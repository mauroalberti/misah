//! Intersection of a triangular surface mesh with a raster grid.
//!
//! Where `raster::intersection` cuts the grid with a single unbounded plane,
//! this cuts it with a surface of arbitrary shape, given as a triangle mesh --
//! the general case a folded or faulted geological surface needs. Each
//! intersection point carries the attitude of the mesh triangle that produced
//! it, so the result is a set of located attitudes, not just a trace.
//!
//! Ported from geoSurfDEM's `IntersectDEM` (C++). Three things are done
//! differently, and deliberately:
//!
//! - **Where the points come from.** The original intersected the two triangle
//!   planes, obtaining a line, then intersected that line with each DEM side in
//!   turn. Here a DEM side is cut against the mesh triangle's plane directly, by
//!   linear interpolation between the signed distances of its two endpoints --
//!   the same construction `raster::intersection` uses at cell edges. It is the
//!   same set of points, reached without the intermediate line, whose
//!   construction in the original needed an arbitrary 100-unit displacement and
//!   an unguarded division, and could return an uninitialised point when all
//!   three of its determinants fell below tolerance. Every point emitted here
//!   lies on the plane and on the side, by construction rather than by luck.
//!
//! - **What gets compared.** The original tested every DEM triangle against
//!   every mesh triangle: on the Malpi example, about 112 million pairs for a
//!   247x200 grid. A DEM *is* a uniform spatial index, so here each mesh
//!   triangle is mapped through the inverse geotransform onto the range of
//!   cells its footprint covers, and only those are visited. No index is built,
//!   because the grid already is one.
//!
//! - **Shared sides.** Adjacent DEM triangles share a side, and the original
//!   emitted the crossing once per triangle: about half the points of its
//!   reference output are duplicates. Crossings are keyed here on the side that
//!   produced them, so each is emitted once, and the count of suppressed
//!   repeats is reported rather than hidden.
//!
//! One consequence of dropping the two-plane construction: the original skipped
//! a pair whose planes were within 0.1 degrees of parallel, since its line was
//! ill-conditioned there. Nothing here is ill-conditioned at small angles, so
//! such a pair is cut normally. Only the exactly coplanar case has no answer to
//! give, and it is counted and skipped.

use std::collections::HashSet;

use crate::geometry::mesh::TriangleMesh;
use crate::geometry::point::Point3D;
use crate::raster::grid::Grid;
use crate::structural::geol_plane::GeologicalPlane;

/// One point where the mesh cuts the topography, with the attitude of the mesh
/// there.
#[derive(Debug, Clone)]
pub struct MeshIntersection {
    pub point: Point3D,
    /// Attitude of the mesh triangle that produced the point. Read off an
    /// upward-oriented normal, so it does not depend on the mesh winding.
    pub attitude: GeologicalPlane,
    /// Index of that triangle in the mesh, for tracing a point back to it.
    pub mesh_triangle: usize,
}

/// What the run did, beyond what it found.
///
/// The C++ original wrote a line to stdout for every degenerate or parallel
/// pair it met, which at a hundred million pairs is not something a kernel can
/// do. These counters carry the same information back to a caller that wants
/// it, and cost nothing to ignore.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct MeshGridIntersectionStats {
    /// Mesh triangles examined, degenerate ones included.
    pub mesh_triangles: usize,
    /// Mesh triangles with no plane of their own: two vertices coincident, or
    /// the three collinear.
    pub degenerate_mesh_triangles: usize,
    /// Mesh triangles skipped outright because their bounding box misses the
    /// grid's. The difference between this and `mesh_triangles` is what the
    /// cell walk actually visited.
    pub mesh_triangles_outside_grid: usize,
    /// DEM triangles compared against a mesh triangle. Against the original's
    /// exhaustive pairing, the ratio to `mesh_triangles` times the DEM triangle
    /// count is what the spatial restriction saved.
    pub dem_triangle_pairs: usize,
    /// DEM sides lying within the mesh triangle's plane, which have no single
    /// crossing point and are skipped.
    pub coplanar_sides: usize,
    /// Crossings suppressed because the DEM side had already been cut by this
    /// same mesh triangle, through its other adjacent triangle.
    pub duplicate_crossings: usize,
}

/// Points where a mesh meets a grid, and a note of how the run went.
#[derive(Debug, Clone, Default)]
pub struct MeshGridIntersection {
    pub intersections: Vec<MeshIntersection>,
    pub stats: MeshGridIntersectionStats,
}

impl MeshGridIntersection {
    pub fn is_empty(&self) -> bool {
        self.intersections.is_empty()
    }

    pub fn len(&self) -> usize {
        self.intersections.len()
    }

    /// Just the points, for feeding an attitude inversion or a writer.
    pub fn points(&self) -> impl Iterator<Item = &Point3D> + '_ {
        self.intersections.iter().map(|i| &i.point)
    }
}

/// Intersect a triangle mesh with a grid.
///
/// The grid is read as the triangulated surface through its node elevations,
/// two triangles per cell, split along the diagonal running from the cell's
/// top-right node to its bottom-left one -- the same decomposition as the C++
/// original. Nodes holding `nodata`, and any non-finite elevation, drop every
/// cell that touches them.
pub fn intersect_mesh_grid(
    mesh: &TriangleMesh,
    grid: &Grid,
    nodata: Option<f64>,
) -> MeshGridIntersection {

    let mut out = MeshGridIntersection::default();

    let (nrows, ncols) = grid.data.dim();

    // A single row or column of nodes spans no cell, so there is no surface to
    // cut.
    if nrows < 2 || ncols < 2 {
        return out;
    }

    let grid_bounds = match grid_bounds(grid, nodata) {
        Some(bounds) => bounds,
        None => return out, // every node is nodata
    };

    let elevation = |row: usize, col: usize| -> f64 {
        let z = grid.data[[row, col]];
        match nodata {
            Some(nd) if z == nd => f64::NAN,
            _ => z,
        }
    };

    let node = |row: usize, col: usize| -> Point3D {
        let (x, y) = grid.transform.node_xy(row, col);
        Point3D::from([x, y, elevation(row, col)])
    };

    // Sides shared between the two triangles of a cell, or between neighbouring
    // cells, are cut once. The key is the side, as its two node indices in a
    // fixed order; the set is cleared per mesh triangle, since a side may
    // legitimately be cut again by the next one.
    let mut cut_sides: HashSet<[usize; 4]> = HashSet::new();

    for (mesh_index, mesh_triangle) in mesh.triangles().enumerate() {

        out.stats.mesh_triangles += 1;

        let plane = match mesh_triangle.to_plane() {
            Some(plane) => plane,
            None => {
                out.stats.degenerate_mesh_triangles += 1;
                continue;
            }
        };

        let triangle_bounds = mesh_triangle.bounds();

        if !bounds_overlap(&triangle_bounds, &grid_bounds) {
            out.stats.mesh_triangles_outside_grid += 1;
            continue;
        }

        let attitude = GeologicalPlane::from_plane(&plane);

        let (row_range, col_range) = match cell_range(grid, &triangle_bounds, nrows, ncols) {
            CellRange::Cells(rows, cols) => (rows, cols),
            CellRange::Outside => {
                out.stats.mesh_triangles_outside_grid += 1;
                continue;
            }
            // Nothing can be ruled out, so the whole grid is walked rather than
            // the triangle being dropped on the quiet.
            CellRange::Unbounded => ((0, nrows - 2), (0, ncols - 2)),
        };

        cut_sides.clear();

        for row in row_range.0..=row_range.1 {
            for col in col_range.0..=col_range.1 {

                // Cell corners, in the order the two DEM triangles refer to
                // them: 0 = (row, col), 1 = (row, col+1), 2 = (row+1, col),
                // 3 = (row+1, col+1).
                let corners = [(row, col), (row, col + 1), (row + 1, col), (row + 1, col + 1)];
                let points = corners.map(|(r, c)| node(r, c));

                if points.iter().any(|p| !p.z().is_finite()) {
                    continue;
                }

                let (cell_z_min, cell_z_max) = points.iter().fold(
                    (f64::INFINITY, f64::NEG_INFINITY),
                    |(lo, hi), p| (lo.min(p.z()), hi.max(p.z())),
                );

                // The accepted points lie inside the mesh triangle, so within
                // its elevation span; a cell wholly above or below it cannot
                // contribute one.
                if cell_z_max < triangle_bounds.0.z() || cell_z_min > triangle_bounds.1.z() {
                    continue;
                }

                let raw_distances = points.map(|p| plane.signed_distance_to_point(&p));

                // A node landing exactly on the plane leaves the sign test with
                // no answer, and would have every side meeting that node report
                // the node itself as a crossing. Nudging it off by the smallest
                // positive double settles the sign the same way
                // `raster::intersection` settles it at cell corners; the point
                // is still emitted, once, from the sides that reach the other
                // side of the plane. A node touched by the plane with no
                // neighbour beyond it is thereby dropped, tangency not being a
                // trace.
                let distances = raw_distances.map(|d| if d == 0.0 { f64::MIN_POSITIVE } else { d });

                for dem_triangle in [[0usize, 1, 2], [2, 1, 3]] {

                    out.stats.dem_triangle_pairs += 1;

                    for side in 0..3 {

                        let (a, b) = (dem_triangle[side], dem_triangle[(side + 1) % 3]);
                        let (da, db) = (distances[a], distances[b]);

                        if raw_distances[a] == 0.0 && raw_distances[b] == 0.0 {
                            out.stats.coplanar_sides += 1;
                            continue;
                        }

                        // Both endpoints on the same side of the plane.
                        if (da > 0.0 && db > 0.0) || (da < 0.0 && db < 0.0) {
                            continue;
                        }

                        if !cut_sides.insert(side_key(corners[a], corners[b])) {
                            out.stats.duplicate_crossings += 1;
                            continue;
                        }

                        // da - db cannot vanish here: that would need da == db,
                        // and equal values are either both zero or of one sign,
                        // both already returned above.
                        let t = da / (da - db);
                        let point = points[a] + points[a].vector_to(&points[b]) * t;

                        if !mesh_triangle.contains_coplanar_point(&point) {
                            continue;
                        }

                        out.intersections.push(MeshIntersection {
                            point,
                            attitude: attitude.clone(),
                            mesh_triangle: mesh_index,
                        });
                    }
                }
            }
        }
    }

    out
}

/// A DEM side, identified by its two nodes with the lower index pair first, so
/// that the two triangles sharing it produce the same key.
#[inline]
fn side_key(a: (usize, usize), b: (usize, usize)) -> [usize; 4] {
    if a <= b {
        [a.0, a.1, b.0, b.1]
    } else {
        [b.0, b.1, a.0, a.1]
    }
}

/// The grid's bounding box: horizontally from its four corner nodes, which
/// covers a rotated transform as well as a north-up one; vertically from the
/// valid elevations. `None` when no node is valid.
fn grid_bounds(grid: &Grid, nodata: Option<f64>) -> Option<(Point3D, Point3D)> {

    let (nrows, ncols) = grid.data.dim();

    if nrows == 0 || ncols == 0 {
        return None;
    }

    let corners = [
        (0, 0),
        (0, ncols - 1),
        (nrows - 1, 0),
        (nrows - 1, ncols - 1),
    ];

    let (mut x_min, mut x_max) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut y_min, mut y_max) = (f64::INFINITY, f64::NEG_INFINITY);

    for (r, c) in corners {
        let (x, y) = grid.transform.node_xy(r, c);
        x_min = x_min.min(x);
        x_max = x_max.max(x);
        y_min = y_min.min(y);
        y_max = y_max.max(y);
    }

    let (mut z_min, mut z_max) = (f64::INFINITY, f64::NEG_INFINITY);

    for &z in grid.data.iter() {
        if !z.is_finite() || nodata == Some(z) {
            continue;
        }
        z_min = z_min.min(z);
        z_max = z_max.max(z);
    }

    if z_min > z_max {
        return None;
    }

    Some((
        Point3D::from([x_min, y_min, z_min]),
        Point3D::from([x_max, y_max, z_max]),
    ))
}

/// Whether two axis-aligned boxes, each given as its lower and upper corner,
/// share any volume. Touching counts.
fn bounds_overlap(a: &(Point3D, Point3D), b: &(Point3D, Point3D)) -> bool {
    (0..3).all(|k| a.0.coords[k] <= b.1.coords[k] && b.0.coords[k] <= a.1.coords[k])
}

/// Which cells a mesh triangle's footprint can reach.
enum CellRange {
    /// Inclusive ranges of cell rows and of cell columns.
    Cells((usize, usize), (usize, usize)),
    /// The footprint falls off the grid entirely.
    Outside,
    /// The transform cannot be inverted, so the footprint cannot be located and
    /// no cell may be ruled out.
    Unbounded,
}

/// Locate a box's horizontal footprint among the grid cells.
///
/// Cell `(r, c)` is the quadrilateral joining node centres `(r, c)` to
/// `(r+1, c+1)`, so a position whose fractional row lies in `[r, r+1)` belongs
/// to cell `r`, and flooring the footprint's fractional bounds gives the range.
/// All four horizontal corners are mapped, not two, because under a rotated
/// transform the extreme row and column need not come from the same corners as
/// the extreme x and y.
///
/// The range is then widened by one cell all round. A footprint ending exactly
/// on a node line has its crossings on the far side of that line, in the cell
/// beyond the one the arithmetic names -- a vertical surface passing through a
/// column of nodes is precisely that case, and without the margin its trace
/// comes out empty.
fn cell_range(grid: &Grid, bounds: &(Point3D, Point3D), nrows: usize, ncols: usize) -> CellRange {

    let (min, max) = bounds;

    let corners = [
        (min.x(), min.y()),
        (min.x(), max.y()),
        (max.x(), min.y()),
        (max.x(), max.y()),
    ];

    let (mut row_min, mut row_max) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut col_min, mut col_max) = (f64::INFINITY, f64::NEG_INFINITY);

    for (x, y) in corners {
        let (row, col) = match grid.transform.node_rc(x, y) {
            Some(rc) => rc,
            None => return CellRange::Unbounded,
        };
        row_min = row_min.min(row);
        row_max = row_max.max(row);
        col_min = col_min.min(col);
        col_max = col_max.max(col);
    }

    let (last_row, last_col) = ((nrows - 1) as f64, (ncols - 1) as f64);

    if row_max < 0.0 || row_min > last_row || col_max < 0.0 || col_min > last_col {
        return CellRange::Outside;
    }

    // The last node row and column start no cell of their own, hence the clamp
    // one short of them.
    let span = |lo: f64, hi: f64, last: f64| {
        let last_cell = last as usize - 1;
        let lo = (lo.floor().clamp(0.0, last - 1.0) as usize).saturating_sub(1);
        let hi = (hi.floor().clamp(0.0, last - 1.0) as usize + 1).min(last_cell);
        (lo, hi)
    };

    CellRange::Cells(
        span(row_min, row_max, last_row),
        span(col_min, col_max, last_col),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    use ndarray::Array2;

    use crate::algebra::vector::Vector3D;
    use crate::geometry::plane::Plane;
    use crate::geometry::triangle::Triangle3D;
    use crate::raster::geotransform::GeoTransform;
    use crate::raster::intersection::intersect_plane_grid;

    /// A north-up grid of unit cells with its lower-left corner at the origin.
    fn unit_grid(nrows: usize, ncols: usize, z: f64) -> Grid {
        Grid {
            transform: GeoTransform {
                data: [0.0, 1.0, 0.0, nrows as f64, 0.0, -1.0],
            },
            epsg_code: 32633,
            data: Array2::from_elem((nrows, ncols), z),
        }
    }

    /// A grid sloping up towards the east, one metre of rise per column.
    fn sloping_grid(nrows: usize, ncols: usize) -> Grid {
        let mut grid = unit_grid(nrows, ncols, 0.0);
        for r in 0..nrows {
            for c in 0..ncols {
                grid.data[[r, c]] = c as f64;
            }
        }
        grid
    }

    /// A triangle of circumradius `radius` lying in `plane` and centred on
    /// `centre`, which must be a point of the plane. Big enough, it covers the
    /// grid outright, so what the kernel finds is decided by the plane and not
    /// by where the patch happens to end.
    fn patch_in_plane(plane: &Plane, centre: Point3D, radius: f64) -> Triangle3D {

        // Along strike, that is horizontal; any in-plane direction would do,
        // and a horizontal plane, having no strike, falls back to east.
        let strike = Vector3D::from([0.0, 0.0, 1.0])
            .cross(&plane.normal)
            .normalize()
            .unwrap_or(Vector3D::from([1.0, 0.0, 0.0]));
        let down_dip = plane.normal.cross(&strike);

        let vertex = |degrees: f64| {
            let a: f64 = degrees.to_radians();
            centre + strike * (radius * a.cos()) + down_dip * (radius * a.sin())
        };

        Triangle3D::new(vertex(90.0), vertex(210.0), vertex(330.0))
    }

    /// A big vertical triangle across the grid at x = `x`, facing east.
    fn vertical_wall(x: f64) -> TriangleMesh {
        TriangleMesh::from_triangle(&Triangle3D::new(
            Point3D::from([x, -100.0, -100.0]),
            Point3D::from([x, 100.0, -100.0]),
            Point3D::from([x, 0.0, 100.0]),
        ))
    }

    #[test]
    fn an_empty_mesh_finds_nothing() {
        let out = intersect_mesh_grid(&TriangleMesh::default(), &unit_grid(10, 10, 0.0), None);

        assert!(out.is_empty());
        assert_eq!(out.stats, MeshGridIntersectionStats::default());
    }

    #[test]
    fn a_mesh_above_the_topography_finds_nothing() {
        let mesh = TriangleMesh::from_triangle(&Triangle3D::new(
            Point3D::from([0.0, 0.0, 500.0]),
            Point3D::from([10.0, 0.0, 500.0]),
            Point3D::from([0.0, 10.0, 500.0]),
        ));

        let out = intersect_mesh_grid(&mesh, &unit_grid(10, 10, 0.0), None);

        assert!(out.is_empty());
    }

    #[test]
    fn a_mesh_beside_the_grid_is_never_walked() {
        // Far to the east of the grid: the bounding-box test alone must settle
        // it, without a single DEM triangle being touched.
        let out = intersect_mesh_grid(&vertical_wall(5000.0), &unit_grid(10, 10, 0.0), None);

        assert!(out.is_empty());
        assert_eq!(out.stats.mesh_triangles_outside_grid, 1);
        assert_eq!(out.stats.dem_triangle_pairs, 0);
    }

    #[test]
    fn a_degenerate_mesh_triangle_is_counted_and_skipped() {
        let mesh = TriangleMesh::new(
            vec![
                Point3D::from([0.0, 0.0, 0.0]),
                Point3D::from([5.0, 5.0, 0.0]),
                Point3D::from([10.0, 10.0, 0.0]),
            ],
            vec![[0, 1, 2]],
        )
        .unwrap();

        let out = intersect_mesh_grid(&mesh, &unit_grid(10, 10, 0.0), None);

        assert!(out.is_empty());
        assert_eq!(out.stats.degenerate_mesh_triangles, 1);
        assert_eq!(out.stats.dem_triangle_pairs, 0);
    }

    #[test]
    fn a_vertical_wall_cuts_a_flat_grid_along_its_own_plane() {
        let out = intersect_mesh_grid(&vertical_wall(4.5), &unit_grid(10, 10, 0.0), None);

        assert!(!out.is_empty());
        for i in &out.intersections {
            assert!((i.point.x() - 4.5).abs() < 1e-9, "point off the wall: {:?}", i.point);
        }
    }

    #[test]
    fn the_reported_attitude_is_that_of_the_source_triangle() {
        // A wall facing east dips 90 degrees towards 90 or 270; from_plane
        // resolves the pair by orienting the normal upward, and for an exactly
        // vertical plane that leaves the winding to decide. Only the dip angle
        // is therefore asserted, being the part that means something here.
        let out = intersect_mesh_grid(&vertical_wall(4.5), &unit_grid(10, 10, 0.0), None);

        assert!(!out.is_empty());
        for i in &out.intersections {
            assert!((i.attitude.dip_angle - 90.0).abs() < 1e-9);
            assert_eq!(i.mesh_triangle, 0);
        }
    }

    #[test]
    fn every_point_lies_on_the_plane_of_its_own_triangle() {
        let mesh = TriangleMesh::from_triangle(&Triangle3D::new(
            Point3D::from([-10.0, -10.0, 2.0]),
            Point3D::from([20.0, -10.0, 8.0]),
            Point3D::from([-10.0, 20.0, 2.0]),
        ));
        let plane = mesh.triangle(0).unwrap().to_plane().unwrap();

        let out = intersect_mesh_grid(&mesh, &sloping_grid(12, 12), None);

        assert!(!out.is_empty());
        for i in &out.intersections {
            assert!(
                plane.distance_to_point(&i.point) < 1e-9,
                "point off plane: {:?}",
                i.point
            );
        }
    }

    #[test]
    fn every_point_falls_within_its_own_triangle() {
        // A patch smaller than the grid, so the containment test does real work
        // rather than accepting everything the plane cuts.
        let triangle = Triangle3D::new(
            Point3D::from([2.0, 2.0, -5.0]),
            Point3D::from([6.0, 2.0, 15.0]),
            Point3D::from([2.0, 6.0, -5.0]),
        );
        let out = intersect_mesh_grid(&TriangleMesh::from_triangle(&triangle), &sloping_grid(12, 12), None);

        assert!(!out.is_empty());
        for i in &out.intersections {
            assert!(triangle.contains_coplanar_point(&i.point));
            let (min, max) = triangle.bounds();
            assert!(i.point.x() >= min.x() - 1e-9 && i.point.x() <= max.x() + 1e-9);
            assert!(i.point.y() >= min.y() - 1e-9 && i.point.y() <= max.y() + 1e-9);
        }
    }

    #[test]
    fn a_shared_side_is_cut_once() {
        // At x = 4.7 the wall passes between node columns, so no cut lands on a
        // node and every crossing belongs to one side alone. The wall spans
        // whole rows of cells, so most of those sides are shared between two
        // DEM triangles: without the keying each point would come out twice, as
        // it does in the C++ original.
        let out = intersect_mesh_grid(&vertical_wall(4.7), &unit_grid(10, 10, 0.0), None);

        assert!(out.stats.duplicate_crossings > 0);

        for (i, a) in out.intersections.iter().enumerate() {
            for b in &out.intersections[i + 1..] {
                assert!(
                    !a.point.approx_eq(&b.point, 1e-9),
                    "point emitted twice: {:?}",
                    a.point
                );
            }
        }
    }

    #[test]
    fn a_wall_running_through_the_nodes_still_traces_the_whole_grid() {
        // x = 4.5 is a node column, so every cut lands on a node: the case the
        // sign test has no answer for. The trace must still cross all ten rows,
        // once each, rather than doubling up or dropping out.
        let out = intersect_mesh_grid(&vertical_wall(4.5), &unit_grid(10, 10, 0.0), None);

        for i in &out.intersections {
            assert!((i.point.x() - 4.5).abs() < 1e-9);
        }

        let ys: HashSet<i64> = out.points().map(|p| (p.y() * 1e6).round() as i64).collect();
        assert_eq!(ys.len(), 10, "one crossing expected per node row");
    }

    #[test]
    fn nodata_cells_are_skipped() {
        let grid = unit_grid(10, 10, -9999.0);

        let out = intersect_mesh_grid(&vertical_wall(4.5), &grid, Some(-9999.0));

        assert!(out.is_empty());
    }

    #[test]
    fn a_nodata_node_removes_only_the_cells_touching_it() {
        let mut grid = sloping_grid(12, 12);
        grid.data[[5, 5]] = -9999.0;

        let with_hole = intersect_mesh_grid(&vertical_wall(5.5), &grid, Some(-9999.0));
        let intact = intersect_mesh_grid(&vertical_wall(5.5), &sloping_grid(12, 12), Some(-9999.0));

        assert!(!with_hole.is_empty());
        assert!(with_hole.len() < intact.len(), "the hole changed nothing");
    }

    #[test]
    fn the_search_is_restricted_to_the_footprint_of_the_mesh() {
        // A patch spanning about two cells of a 40x40 grid. Exhaustive pairing,
        // as the C++ original did it, would test 2 * 39 * 39 = 3042 DEM
        // triangles against it; the cell walk must stay near the handful its
        // footprint covers, plus the one-cell margin around them.
        let mesh = TriangleMesh::from_triangle(&Triangle3D::new(
            Point3D::from([10.0, 10.0, -5.0]),
            Point3D::from([12.0, 10.0, 5.0]),
            Point3D::from([10.0, 12.0, -5.0]),
        ));

        let out = intersect_mesh_grid(&mesh, &unit_grid(40, 40, 0.0), None);

        assert!(!out.is_empty());
        assert!(
            out.stats.dem_triangle_pairs < 100,
            "the cell restriction is not biting: {} pairs",
            out.stats.dem_triangle_pairs
        );
    }

    #[test]
    fn a_planar_mesh_agrees_with_the_plane_kernel() {
        // The two kernels cut different sides -- marching squares walks cell
        // edges, this one also walks the diagonals -- so their point sets are
        // not the same. What must agree is the trace they describe: both on the
        // plane, and over the same ground.
        let grid = sloping_grid(20, 20);

        let centre = Point3D::from([10.0, 10.0, 9.5]);
        let plane = GeologicalPlane::new(90.0, 45.0).to_plane(centre).unwrap();
        let triangle = patch_in_plane(&plane, centre, 1000.0);

        let from_mesh = intersect_mesh_grid(&TriangleMesh::from_triangle(&triangle), &grid, None);
        let from_plane = intersect_plane_grid(&plane, &grid, None);

        assert!(!from_mesh.is_empty() && !from_plane.is_empty());

        let span = |xs: Vec<f64>| {
            xs.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &v| {
                (lo.min(v), hi.max(v))
            })
        };

        let mesh_y = span(from_mesh.points().map(|p| p.y()).collect());
        let plane_y = span(from_plane.points.iter().map(|p| p.y()).collect());

        // Within a cell: the two kernels place their vertices on different
        // sides, so the extremes need not coincide exactly.
        assert!((mesh_y.0 - plane_y.0).abs() <= 1.0, "{:?} vs {:?}", mesh_y, plane_y);
        assert!((mesh_y.1 - plane_y.1).abs() <= 1.0, "{:?} vs {:?}", mesh_y, plane_y);

        for p in from_mesh.points() {
            assert!(plane.distance_to_point(p) < 1e-9);
        }
    }

    #[test]
    fn a_grid_too_small_to_hold_a_cell_finds_nothing() {
        let grid = Grid {
            transform: GeoTransform { data: [0.0, 1.0, 0.0, 1.0, 0.0, -1.0] },
            epsg_code: 32633,
            data: Array2::from_elem((1, 5), 0.0),
        };

        assert!(intersect_mesh_grid(&vertical_wall(2.5), &grid, None).is_empty());
    }

    #[test]
    fn a_side_lying_in_the_plane_is_counted_and_skipped() {
        // A horizontal patch at exactly the elevation of a flat grid: every DEM
        // side lies in its plane, and none has a single crossing point.
        let mesh = TriangleMesh::from_triangle(&Triangle3D::new(
            Point3D::from([-10.0, -10.0, 0.0]),
            Point3D::from([100.0, -10.0, 0.0]),
            Point3D::from([-10.0, 100.0, 0.0]),
        ));

        let out = intersect_mesh_grid(&mesh, &unit_grid(6, 6, 0.0), None);

        assert!(out.is_empty());
        assert!(out.stats.coplanar_sides > 0);
    }

    #[test]
    fn results_carry_the_index_of_the_triangle_that_made_them() {
        // Two disjoint walls, so the index is the only way to tell which of
        // them a point came from.
        let mesh = TriangleMesh::new(
            vec![
                Point3D::from([3.0, -100.0, -100.0]),
                Point3D::from([3.0, 100.0, -100.0]),
                Point3D::from([3.0, 0.0, 100.0]),
                Point3D::from([7.0, -100.0, -100.0]),
                Point3D::from([7.0, 100.0, -100.0]),
                Point3D::from([7.0, 0.0, 100.0]),
            ],
            vec![[0, 1, 2], [3, 4, 5]],
        )
        .unwrap();

        let out = intersect_mesh_grid(&mesh, &unit_grid(10, 10, 0.0), None);

        for i in &out.intersections {
            let expected_x = if i.mesh_triangle == 0 { 3.0 } else { 7.0 };
            assert!((i.point.x() - expected_x).abs() < 1e-9);
        }
        assert!(out.intersections.iter().any(|i| i.mesh_triangle == 0));
        assert!(out.intersections.iter().any(|i| i.mesh_triangle == 1));
    }
}
