//! Intersection of a plane with a raster grid.
//!
//! The intersection is the zero-level set of the plane's signed distance,
//! sampled at the grid nodes and extracted by marching squares. Working on the
//! signed distance rather than on a `z = z(x, y)` plane expression keeps vertical
//! planes -- which have no single-valued z -- as an ordinary case, and costs one
//! pass over the cells.

use std::collections::HashMap;

use crate::geometry::plane::Plane;
use crate::geometry::point::Point3D;
use crate::raster::grid::Grid;

/// Corner numbering: 0 = (r, c), 1 = (r, c+1), 2 = (r+1, c+1), 3 = (r+1, c).
/// Edges are named after the corner pair they join.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Edge {
    Top,
    Right,
    Bottom,
    Left,
}

impl Edge {
    #[inline]
    fn corners(self) -> (usize, usize) {
        match self {
            Edge::Top => (0, 1),
            Edge::Right => (1, 2),
            Edge::Bottom => (3, 2),
            Edge::Left => (0, 3),
        }
    }
}

/// Marching-squares chords per case index; bit i is set when corner i is positive.
/// Complementary cases (i and 15-i) share their chords.
#[inline]
fn case_edges(case: u8) -> &'static [(Edge, Edge)] {
    match case {
        0 | 15 => &[],
        1 | 14 => &[(Edge::Left, Edge::Top)],
        2 | 13 => &[(Edge::Top, Edge::Right)],
        3 | 12 => &[(Edge::Left, Edge::Right)],
        4 | 11 => &[(Edge::Right, Edge::Bottom)],
        6 | 9 => &[(Edge::Top, Edge::Bottom)],
        7 | 8 => &[(Edge::Left, Edge::Bottom)],
        5 => &[(Edge::Left, Edge::Top), (Edge::Right, Edge::Bottom)],
        10 => &[(Edge::Left, Edge::Bottom), (Edge::Top, Edge::Right)],
        _ => unreachable!("case index is four bits"),
    }
}

/// The alternative pairing for a saddle cell.
#[inline]
fn saddle_alternative(case: u8) -> &'static [(Edge, Edge)] {
    match case {
        5 => &[(Edge::Left, Edge::Bottom), (Edge::Top, Edge::Right)],
        10 => &[(Edge::Left, Edge::Top), (Edge::Right, Edge::Bottom)],
        _ => unreachable!("only cases 5 and 10 are saddles"),
    }
}

/// Trace of a plane over a grid.
#[derive(Debug, Clone, Default)]
pub struct PlaneGridIntersection {
    /// Intersection vertices.
    pub points: Vec<Point3D>,
    /// Index pairs into `points`, one per marching-squares chord. They let a
    /// caller rebuild polylines; the bare vertices feed attitude inversion.
    pub segments: Vec<[usize; 2]>,
}

impl PlaneGridIntersection {
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }
}

/// Intersect a plane with a grid.
///
/// Nodes holding `nodata`, and any non-finite elevation, are excluded together
/// with every cell that touches them.
pub fn intersect_plane_grid(
    plane: &Plane,
    grid: &Grid,
    nodata: Option<f64>,
) -> PlaneGridIntersection {
    let (nrows, ncols) = grid.data.dim();

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

    let signed_distance =
        |row: usize, col: usize| -> f64 { plane.signed_distance_to_point(&node(row, col)) };

    let mut out = PlaneGridIntersection::default();
    // Vertices are shared between adjacent cells; key them by edge identity so a
    // crossing is emitted once and segments stay topologically connected.
    let mut vertex_of_edge: HashMap<[usize; 4], usize> = HashMap::new();

    for r in 0..nrows.saturating_sub(1) {
        for c in 0..ncols.saturating_sub(1) {
            let corners = [(r, c), (r, c + 1), (r + 1, c + 1), (r + 1, c)];

            let mut vals = [0.0f64; 4];
            let mut usable = true;
            for (i, &(cr, cc)) in corners.iter().enumerate() {
                let f = signed_distance(cr, cc);
                if !f.is_finite() {
                    usable = false;
                    break;
                }
                // A node exactly on the plane would make the sign test ill-defined;
                // nudging it off by an ulp-scale epsilon keeps the topology consistent.
                vals[i] = if f == 0.0 { f64::MIN_POSITIVE } else { f };
            }
            if !usable {
                continue;
            }

            let mut case = 0u8;
            for (i, &v) in vals.iter().enumerate() {
                if v > 0.0 {
                    case |= 1 << i;
                }
            }

            let mut chords = case_edges(case);
            if chords.is_empty() {
                continue;
            }
            if case == 5 || case == 10 {
                // Saddle: the cell centre tells which of the two pairings is right.
                let centre = vals.iter().sum::<f64>() / 4.0;
                if (centre > 0.0) != (case == 5) {
                    chords = saddle_alternative(case);
                }
            }

            for &(e0, e1) in chords {
                let mut ends = [0usize; 2];
                for (slot, edge) in [e0, e1].iter().enumerate() {
                    let (a, b) = edge.corners();
                    let (r0, c0) = corners[a];
                    let (r1, c1) = corners[b];

                    let key = [r0, c0, r1, c1];
                    ends[slot] = match vertex_of_edge.get(&key) {
                        Some(&idx) => idx,
                        None => {
                            let pa = node(r0, c0);
                            let pb = node(r1, c1);
                            let fa = signed_distance(r0, c0);
                            let fb = signed_distance(r1, c1);
                            let t = fa / (fa - fb);

                            let mut coords = [0.0f64; 3];
                            for k in 0..3 {
                                coords[k] = pa.coords[k] + t * (pb.coords[k] - pa.coords[k]);
                            }

                            let idx = out.points.len();
                            out.points.push(Point3D::from(coords));
                            vertex_of_edge.insert(key, idx);
                            idx
                        }
                    };
                }
                out.segments.push(ends);
            }
        }
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    use ndarray::Array2;

    use crate::algebra::Vector;
    use crate::raster::geotransform::GeoTransform;
    use crate::structural::geol_plane::GeologicalPlane;

    fn unit_grid(nrows: usize, ncols: usize, z: f64) -> Grid {
        Grid {
            transform: GeoTransform {
                data: [0.0, 1.0, 0.0, nrows as f64, 0.0, -1.0],
            },
            epsg_code: 32633,
            data: Array2::from_elem((nrows, ncols), z),
        }
    }

    fn plane_at(x: f64, y: f64, z: f64, dip_dir: f64, dip_angle: f64) -> Plane {
        GeologicalPlane::new(dip_dir, dip_angle)
            .to_plane(Point3D::from([x, y, z]))
            .expect("a geological plane always has a non-zero normal")
    }

    #[test]
    fn horizontal_plane_normal_points_up() {
        let plane = plane_at(0.0, 0.0, 0.0, 0.0, 0.0);
        let n = plane.normal.coords;
        assert!(n[0].abs() < 1e-12 && n[1].abs() < 1e-12 && (n[2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn vertical_plane_normal_is_horizontal() {
        let plane = plane_at(0.0, 0.0, 0.0, 90.0, 90.0);
        let n = plane.normal.coords;
        assert!((n[0] - 1.0).abs() < 1e-12 && n[2].abs() < 1e-12);
    }

    #[test]
    fn grid_entirely_below_the_plane_yields_nothing() {
        let grid = unit_grid(10, 10, 0.0);
        let out = intersect_plane_grid(&plane_at(0.0, 0.0, 100.0, 0.0, 0.0), &grid, None);
        assert!(out.is_empty());
        assert!(out.segments.is_empty());
    }

    #[test]
    fn tilted_plane_cuts_a_flat_grid_along_a_line() {
        let grid = unit_grid(10, 10, 0.0);
        let plane = plane_at(5.0, 5.0, 0.0, 90.0, 45.0);
        let out = intersect_plane_grid(&plane, &grid, None);

        assert!(!out.is_empty());
        for p in &out.points {
            assert!(
                plane.signed_distance_to_point(p).abs() < 1e-9,
                "point off plane: {:?}",
                p
            );
        }
    }

    #[test]
    fn horizontal_plane_cuts_a_slope_at_its_own_elevation() {
        let (nrows, ncols) = (10, 10);
        let mut grid = unit_grid(nrows, ncols, 0.0);
        for r in 0..nrows {
            for c in 0..ncols {
                grid.data[[r, c]] = c as f64;
            }
        }
        let out = intersect_plane_grid(&plane_at(0.0, 0.0, 4.5, 0.0, 0.0), &grid, None);

        assert!(!out.is_empty());
        for p in &out.points {
            assert!((p.z() - 4.5).abs() < 1e-9, "contour off level: {:?}", p);
        }
    }

    #[test]
    fn nodata_cells_are_skipped() {
        let grid = unit_grid(10, 10, -9999.0);
        let out = intersect_plane_grid(
            &plane_at(5.0, 5.0, 0.0, 90.0, 45.0),
            &grid,
            Some(-9999.0),
        );
        assert!(out.is_empty());
    }

    #[test]
    fn segments_index_existing_vertices() {
        let grid = unit_grid(10, 10, 0.0);
        let out = intersect_plane_grid(&plane_at(5.0, 5.0, 0.0, 90.0, 45.0), &grid, None);
        for s in &out.segments {
            assert!(s[0] < out.points.len() && s[1] < out.points.len());
        }
    }

    #[test]
    fn chords_stay_within_a_cell() {
        let grid = unit_grid(10, 10, 0.0);
        let out = intersect_plane_grid(&plane_at(5.0, 5.0, 0.0, 90.0, 45.0), &grid, None);
        for s in &out.segments {
            let d = out.points[s[0]].distance(&out.points[s[1]]);
            assert!(d < 3.0, "chord longer than a cell: {}", d);
        }
    }

    #[test]
    fn plane_from_three_points_agrees_with_the_geological_one() {
        // An independent construction of the same plane, as a cross-check on the
        // dip-direction convention.
        let geological = plane_at(0.0, 0.0, 0.0, 90.0, 45.0);
        let strike = Point3D::from([0.0, 1.0, 0.0]);
        let down_dip = Point3D::from([1.0, 0.0, -1.0]);
        let built = Plane::from_three_points(Point3D::from([0.0, 0.0, 0.0]), strike, down_dip)
            .expect("three distinct non-collinear points");

        let dot = Vector::new(geological.normal.coords).dot(&Vector::new(built.normal.coords));
        assert!(dot.abs() > 1.0 - 1e-12, "normals not parallel: {}", dot);
    }
}
