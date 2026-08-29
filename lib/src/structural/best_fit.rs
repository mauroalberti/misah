//! Best-fit geological planes from located points.
//!
//! The consumer `raster::mesh_intersection` was written for. That kernel cuts
//! a geological surface against a DEM and returns where the two meet; this one
//! reads an attitude back out of those points, by fitting a plane to each
//! neighbourhood of them. Between the two, a mapped surface goes to a field of
//! measured orientations without either step knowing about the other -- the
//! input here is bare coordinates, so a trace digitised from a map or a set of
//! GPS readings along an outcrop feeds it just as well.
//!
//! A port of geoSurfDEM's `BestFitGeoplanes`, whose numerical core is the
//! Fortran `invert_main_attitude` in `GeoInversions`: centre the points of a
//! cell, take the SVD of the resulting m-by-3 matrix through LAPACK, and read
//! the plane normal off the right singular vector of the smallest singular
//! value. Unlike `ForwardStress.f95` -- the other Fortran this crate has
//! ported, where a missing `implicit none` silently truncated a double to a
//! single and moved the answer -- `GeoInversions` compiles clean under
//! `-fimplicit-none`, so the reference numbers can be trusted as they stand.
//! `algebra::eigen` explains why the port reaches the same vectors through an
//! eigen-decomposition of the 3x3 scatter matrix rather than through LAPACK.
//!
//! Three things depart from the C++ deliberately, and one of them changes
//! results rather than only cost:
//!
//! - **Nearly collinear cells are counted, not returned as attitudes.** This
//!   is the substantive difference. A plane is determined by a cloud only if
//!   the cloud spreads in two directions: the test is `s2/s1`, the second
//!   spread against the first. Where the points fall along a line -- the
//!   ordinary condition for a contact crossing a smooth slope, not a rare one
//!   -- every plane through that line fits equally well, and the one returned
//!   is whichever direction rounding happened to favour.
//!
//!   The reference computes this ratio, as `logratio_s(1)`, alongside two
//!   others, writes all three to its results file, and publishes an attitude
//!   regardless of any of them. On the Timpa San Lorenzo dataset that costs it
//!   three of its 311 cells: two carrying attitudes 13.2 and 12.5 degrees away
//!   from the plane whose intersection generated the points, both at `s2/s1`
//!   near 3e-6 where a sound cell sits near 3e-2; and a third where `s2` and
//!   `s3` are both zero, so its own diagnostic column is `-log10(0/0)` -- a
//!   NaN printed beside a confident dip direction.
//!
//!   Which of the three ratios to test on is not obvious and is worth stating,
//!   since the plausible choice is the wrong one. `-log10(s3/s2)` reads like
//!   the measure of a good plane, and on this dataset it has no power to
//!   separate at all: the two 13-degree errors score 3.8 and 3.5 on it,
//!   better than most cells that are correct to a thousandth of a degree.
//!   When `s2` is itself noise, a ratio taken against it is a ratio of two
//!   noises.
//! - **Duplicate points are found through the grid, not by scanning.** The
//!   reference compares each candidate against every point already kept, which
//!   is quadratic; on the 1900-point Timpa San Lorenzo dataset that is under
//!   two million comparisons and costs nothing, but the mesh kernel can hand
//!   over far more than 1900 points. Binning at the coincidence distance makes
//!   two points that could be coincident always fall in the same or an
//!   adjacent bin, so a 3x3 neighbourhood of bins holds every candidate and
//!   the answer is unchanged.
//! - **Coincidence is judged in three dimensions, not two.** The reference
//!   compares only x and y. For points that came off a DEM the two agree,
//!   since a single-valued surface cannot put two elevations at one map
//!   location; for points from anywhere else -- an overturned contact, a cave
//!   survey, a set of readings up a cliff -- the 2D test throws away exactly
//!   the measurements that constrain a steep plane.

use std::collections::HashMap;

use crate::algebra::eigen::symmetric_eigen;
use crate::algebra::vector::Vector;
use crate::geometry::plane::Plane;
use crate::geometry::point::Point3D;
use crate::structural::geol_plane::GeologicalPlane;

/// Below this ratio of the largest singular value, a singular value is taken
/// as zero.
///
/// Not a tolerance on the data but on the arithmetic, and it is looser than it
/// looks for a reason worth stating. Forming the scatter matrix before
/// decomposing it puts eigenvalues at best within about `EPSILON` of each
/// other relatively -- and a singular value is the square root of an
/// eigenvalue, which halves the exponent. A ratio resolvable to 1e-15 among
/// the eigenvalues is resolvable only to about 3e-8 among the singular values.
/// Measured on exactly collinear points, the second singular value comes out
/// around 1e-7 of the first rather than at zero, so a floor tighter than this
/// would call a straight line two-dimensional.
const SINGULAR_VALUE_FLOOR: f64 = 1.0e-6;

/// How well the fitted plane is determined by the points it came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FitQuality {
    /// The points spread in two directions and lie close to one plane: the
    /// normal is determined, and the attitude means something.
    Determined,
    /// The points fall along a line. Every plane containing that line fits
    /// them equally well, so the attitude returned is arbitrary -- not
    /// imprecise, arbitrary -- and says nothing about the surface.
    Collinear,
}

/// A plane fitted to a set of points, with what the fit is worth.
#[derive(Debug, Clone)]
pub struct BestFitPlane {
    pub plane: GeologicalPlane,
    /// The centroid the fit was taken about, and the point the plane passes
    /// through.
    pub centre: Point3D,
    /// The three singular values of the centred coordinates, descending: the
    /// spread of the cloud along its own principal directions, in the units
    /// the points were given in.
    pub singular_values: [f64; 3],
    pub points: usize,
}

impl BestFitPlane {

    /// `-log10(s2/s1)`: how far the second spread falls below the first, in
    /// orders of magnitude.
    ///
    /// The quantity that says whether the normal is pinned down, and note the
    /// direction -- **large is bad**. Zero means the cloud is equally wide in
    /// its two largest directions, which determines a plane as firmly as
    /// points can; large means the cloud is a line, and the fitted plane is
    /// free to rotate about it. `f64::INFINITY` where the second spread is
    /// nothing but rounding, which is the exactly collinear case the reference
    /// turns into a NaN.
    ///
    /// This is the reference's `logratio_s(1)`. Its `logratio_s(3)`,
    /// `-log10(s3/s2)`, looks like it should serve here and does not: see the
    /// module documentation.
    pub fn collinearity(&self) -> f64 {

        let [s1, s2, _] = self.singular_values;

        if s2 <= s1 * SINGULAR_VALUE_FLOOR {
            return f64::INFINITY;
        }

        -(s2 / s1).log10()
    }

    /// Whether the points determine a plane at all.
    ///
    /// `max_collinearity` is in the orders of magnitude `collinearity`
    /// returns, and a fit is kept when it is at or below that;
    /// `DEFAULT_MAX_COLLINEARITY` is the value `best_fit_geoplanes` uses.
    pub fn quality(&self, max_collinearity: f64) -> FitQuality {
        if self.collinearity() <= max_collinearity {
            FitQuality::Determined
        } else {
            FitQuality::Collinear
        }
    }

    /// Root-mean-square distance of the points from the fitted plane.
    ///
    /// `s3 / sqrt(n)`: the smallest singular value is the total scatter
    /// perpendicular to the plane, in the same units as the coordinates. The
    /// residual of the fit, where `collinearity` is its conditioning -- and
    /// the two are independent, which is the trap. A nearly collinear cell has
    /// an excellent residual and a meaningless attitude: the points do lie on
    /// the plane returned, and so do they on every other plane through the
    /// same line. A residual alone never detects this.
    pub fn rms_distance(&self) -> f64 {
        self.singular_values[2] / (self.points as f64).sqrt()
    }
}

/// Three orders of magnitude: a fit is kept while its second spread is above a
/// thousandth of its first.
///
/// Chosen by measurement rather than by taste. On the Timpa San Lorenzo
/// dataset, where the true attitude of every point is known because the points
/// were generated by intersecting known planes with a DEM, cells correct to
/// better than a tenth of a degree run up to a collinearity of 3.3 and the
/// wrong ones begin at 3.1 -- the two populations overlap, so no threshold
/// separates them cleanly. At 3.0 the field keeps 288 of 308 cells and its
/// worst attitude is 0.089 degrees off; at 4.0 it keeps 303 and the worst is
/// 0.873; with no threshold at all, 38.4.
///
/// Tenths of a degree is far below what a compass reading is worth, so this
/// buys a margin nobody will measure. It is set here rather than lower because
/// the points along a contact are genuinely elongated -- a cell wide enough to
/// make them equidimensional is wide enough to average across a fold hinge --
/// and discarding sound cells to chase a thousandth of a degree would trade
/// coverage for a precision the data does not have.
pub const DEFAULT_MAX_COLLINEARITY: f64 = 3.0;

/// Fit a plane to a set of points by least squares perpendicular distance.
///
/// `None` for fewer than three points, or for points so nearly coincident that
/// they have no spread at all -- neither is a badly conditioned fit, both are
/// the absence of one. A collinear set does return a fit, whose `quality` says
/// what it is worth; the decision to use it or drop it belongs to the caller,
/// who knows what the points are.
pub fn fit_plane(points: &[Point3D]) -> Option<BestFitPlane> {

    if points.len() < 3 {
        return None;
    }

    let n = points.len() as f64;
    let mut centre = [0.0f64; 3];
    for point in points {
        for (sum, coord) in centre.iter_mut().zip(point.coords) {
            *sum += coord;
        }
    }
    centre = centre.map(|sum| sum / n);

    // The scatter matrix of the centred coordinates: X^T X, whose eigenvectors
    // are the principal directions of the cloud and whose eigenvalues are the
    // squares of the singular values of X.
    let mut scatter = [[0.0f64; 3]; 3];
    for point in points {
        let d = [
            point.coords[0] - centre[0],
            point.coords[1] - centre[1],
            point.coords[2] - centre[2],
        ];
        for i in 0..3 {
            for j in i..3 {
                scatter[i][j] += d[i] * d[j];
            }
        }
    }

    let (eigenvalues, eigenvectors) = symmetric_eigen(scatter);

    // A scatter matrix is positive semi-definite, so a negative eigenvalue is
    // rounding about zero and clamping it is the honest reading, not a repair.
    let singular_values = eigenvalues.map(|v| v.max(0.0).sqrt());

    if singular_values[0] <= 0.0 {
        // No spread in any direction: every point is the same point.
        return None;
    }

    // The direction of least scatter, which is the normal of the plane the
    // points sit closest to.
    let normal = Vector::new([eigenvectors[0][2], eigenvectors[1][2], eigenvectors[2][2]]);

    let centre = Point3D::from(centre);
    let plane = Plane::new(centre, normal)?;

    Some(BestFitPlane {
        plane: GeologicalPlane::from_plane(&plane),
        centre,
        singular_values,
        points: points.len(),
    })
}

/// The regular grid a field of fits is laid out on.
///
/// Anchored at the top-left of the input points' own extent, as the reference
/// does: `row` runs south from `y_max`, `col` east from `x_min`.
#[derive(Debug, Clone, Copy)]
pub struct CellGrid {
    pub x_min: f64,
    pub y_max: f64,
    pub cell_size: f64,
    pub rows: usize,
    pub columns: usize,
}

impl CellGrid {

    /// The grid covering these points at this cell size.
    ///
    /// `None` for no points, there being no extent to cover.
    pub fn covering(points: &[Point3D], cell_size: f64) -> Option<Self> {

        let first = points.first()?;
        let (mut x_min, mut x_max) = (first.coords[0], first.coords[0]);
        let (mut y_min, mut y_max) = (first.coords[1], first.coords[1]);

        for point in points {
            x_min = x_min.min(point.coords[0]);
            x_max = x_max.max(point.coords[0]);
            y_min = y_min.min(point.coords[1]);
            y_max = y_max.max(point.coords[1]);
        }

        Some(Self {
            x_min,
            y_max,
            cell_size,
            rows: ((y_max - y_min) / cell_size) as usize + 1,
            columns: ((x_max - x_min) / cell_size) as usize + 1,
        })
    }

    /// The cell a point falls in.
    pub fn cell_of(&self, point: &Point3D) -> (usize, usize) {
        (
            ((self.y_max - point.coords[1]) / self.cell_size) as usize,
            ((point.coords[0] - self.x_min) / self.cell_size) as usize,
        )
    }

    /// The map coordinates of a cell's centre.
    pub fn centre_of(&self, row: usize, col: usize) -> (f64, f64) {
        (
            self.x_min + (col as f64 + 0.5) * self.cell_size,
            self.y_max - (row as f64 + 0.5) * self.cell_size,
        )
    }
}

/// A fit, and where on the grid it came from.
#[derive(Debug, Clone)]
pub struct CellFit {
    pub row: usize,
    pub column: usize,
    /// The cell's centre, which is where the reference posts the result --
    /// not the centroid of the points, which is in `fit.centre`.
    pub cell_centre: (f64, f64),
    pub fit: BestFitPlane,
}

/// What a run did, including what it declined to do.
#[derive(Debug, Clone, Copy, Default)]
pub struct FieldStats {
    pub input_points: usize,
    /// After coincident points were merged.
    pub distinct_points: usize,
    pub non_empty_cells: usize,
    pub cells_fitted: usize,
    /// Cells holding points but fewer than three.
    pub cells_too_few_points: usize,
    /// Cells whose points were collinear, so that no attitude follows from
    /// them. Not a failure: on a smooth slope this is most of them.
    pub cells_collinear: usize,
}

/// A field of attitudes, and the grid it is posted on.
#[derive(Debug, Clone)]
pub struct GeoplaneField {
    pub fits: Vec<CellFit>,
    pub grid: CellGrid,
    /// Points per cell, row-major over the whole grid -- the count raster the
    /// reference writes out beside its results, and the thing to look at when
    /// a field comes back thinner than expected.
    pub counts: Vec<usize>,
    pub stats: FieldStats,
}

/// Fit a plane in every cell of a regular grid over the points.
///
/// `cell_size` sets both the grid and the scale the attitudes are averaged
/// over: too small and no cell holds three points, too large and a fold is
/// flattened into a single meaningless plane. `coincidence_distance` merges
/// points closer than it, which matters because an intersection kernel emits a
/// point per crossing and a shared edge is crossed twice -- duplicates would
/// otherwise weight one location more heavily than its neighbours. Points are
/// merged, not summed: the first of a coincident group is kept.
///
/// `max_collinearity` is the threshold on `BestFitPlane::collinearity` above
/// which a cell is counted as collinear and left out of `fits`.
///
/// `None` when there are no points, or when `cell_size` or
/// `coincidence_distance` is not a positive finite number.
pub fn best_fit_geoplanes(
    points: &[Point3D],
    cell_size: f64,
    coincidence_distance: f64,
    max_collinearity: f64,
) -> Option<GeoplaneField> {

    if !(cell_size.is_finite() && cell_size > 0.0) {
        return None;
    }
    if !(coincidence_distance.is_finite() && coincidence_distance > 0.0) {
        return None;
    }

    let distinct = distinct_points(points, coincidence_distance);
    let grid = CellGrid::covering(&distinct, cell_size)?;

    let mut by_cell: HashMap<(usize, usize), Vec<Point3D>> = HashMap::new();
    for point in &distinct {
        by_cell.entry(grid.cell_of(point)).or_default().push(*point);
    }

    let mut counts = vec![0usize; grid.rows * grid.columns];
    for (&(row, column), cell_points) in &by_cell {
        counts[row * grid.columns + column] = cell_points.len();
    }

    let mut stats = FieldStats {
        input_points: points.len(),
        distinct_points: distinct.len(),
        non_empty_cells: by_cell.len(),
        ..Default::default()
    };

    let mut fits = Vec::new();

    for (&(row, column), cell_points) in &by_cell {

        let Some(fit) = fit_plane(cell_points) else {
            stats.cells_too_few_points += 1;
            continue;
        };

        if fit.quality(max_collinearity) == FitQuality::Collinear {
            stats.cells_collinear += 1;
            continue;
        }

        stats.cells_fitted += 1;
        fits.push(CellFit {
            row,
            column,
            cell_centre: grid.centre_of(row, column),
            fit,
        });
    }

    // A HashMap hands its entries back in whatever order it likes, and a field
    // of attitudes that reshuffles between runs is a nuisance to diff, to plot
    // and to test. Row-major, as the grid is read.
    fits.sort_by_key(|f| (f.row, f.column));

    Some(GeoplaneField { fits, grid, counts, stats })
}

/// The points with coincident ones merged, keeping the first of each group.
///
/// Binned at the coincidence distance, so that two points within it of each
/// other always share a bin or lie in adjacent ones and the nine bins around a
/// candidate hold every point it could be coincident with.
fn distinct_points(points: &[Point3D], coincidence_distance: f64) -> Vec<Point3D> {

    let mut kept: Vec<Point3D> = Vec::new();
    let mut bins: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();

    let bin_of = |p: &Point3D| {
        (
            (p.coords[0] / coincidence_distance).floor() as i64,
            (p.coords[1] / coincidence_distance).floor() as i64,
            (p.coords[2] / coincidence_distance).floor() as i64,
        )
    };

    let threshold_squared = coincidence_distance * coincidence_distance;

    'candidates: for point in points {

        let (bx, by, bz) = bin_of(point);

        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let Some(neighbours) = bins.get(&(bx + dx, by + dy, bz + dz)) else {
                        continue;
                    };
                    for &index in neighbours {
                        let other: &Point3D = &kept[index];
                        let separation: f64 = (0..3)
                            .map(|i| (point.coords[i] - other.coords[i]).powi(2))
                            .sum();
                        if separation < threshold_squared {
                            continue 'candidates;
                        }
                    }
                }
            }
        }

        bins.entry((bx, by, bz)).or_default().push(kept.len());
        kept.push(*point);
    }

    kept
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Points on a plane of the given attitude through the origin, sampled at
    /// the given map positions.
    fn on_plane(azimuth: f64, dip: f64, xy: &[(f64, f64)]) -> Vec<Point3D> {

        let plane = GeologicalPlane::new(azimuth, dip);
        let n = plane.normal_vector().coords;

        // n . (p - origin) = 0, solved for z.
        xy.iter()
            .map(|&(x, y)| Point3D::from([x, y, -(n[0] * x + n[1] * y) / n[2]]))
            .collect()
    }

    const SPREAD: [(f64, f64); 6] = [
        (0.0, 0.0), (10.0, 3.0), (4.0, 11.0), (-7.0, 5.0), (-2.0, -9.0), (8.0, -6.0),
    ];

    #[test]
    fn a_plane_is_recovered_from_points_lying_on_it() {
        for &(azimuth, dip) in &[(0.0, 30.0), (135.0, 35.0), (270.0, 70.0), (45.0, 5.0)] {
            let fit = fit_plane(&on_plane(azimuth, dip, &SPREAD)).expect("a fit");

            assert!(
                (fit.plane.azimuth - azimuth).abs() < 1e-6,
                "azimuth {} against {}",
                fit.plane.azimuth,
                azimuth
            );
            assert!((fit.plane.dip_angle - dip).abs() < 1e-6);
            assert_eq!(fit.quality(DEFAULT_MAX_COLLINEARITY), FitQuality::Determined);
        }
    }

    #[test]
    fn points_exactly_on_a_plane_have_no_residual() {
        let fit = fit_plane(&on_plane(135.0, 35.0, &SPREAD)).expect("a fit");

        // Not zero, and this is the floor the scatter-matrix route imposes
        // rather than anything about the data: on coordinates of order ten,
        // a residual below a millionth cannot be told from none at all.
        assert!(fit.rms_distance() < 1e-6, "residual {}", fit.rms_distance());

        // Spread in two directions, so the plane is well determined: the
        // second singular value is within an order of magnitude of the first.
        assert!(fit.collinearity() < 1.0, "collinearity {}", fit.collinearity());
    }

    #[test]
    fn a_horizontal_plane_comes_back_horizontal() {
        let points: Vec<Point3D> = SPREAD.iter().map(|&(x, y)| Point3D::from([x, y, 17.0])).collect();

        let fit = fit_plane(&points).expect("a fit");

        assert!(fit.plane.dip_angle < 1e-9);
        assert!((fit.centre.coords[2] - 17.0).abs() < 1e-9);
    }

    #[test]
    fn collinear_points_are_reported_as_such() {
        // The failure the reference does not catch: every plane through this
        // line fits perfectly, so the attitude returned means nothing, and
        // the residual -- excellent -- does not say so.
        let points: Vec<Point3D> = (0..8)
            .map(|i| Point3D::from([i as f64, 2.0 * i as f64, 3.0 * i as f64]))
            .collect();

        let fit = fit_plane(&points).expect("a fit is returned, with its quality");

        assert!(fit.rms_distance() < 1e-6, "the residual looks perfect");
        assert!(fit.collinearity().is_infinite(), "and yet nothing is determined");
        assert_eq!(fit.quality(DEFAULT_MAX_COLLINEARITY), FitQuality::Collinear);
    }

    #[test]
    fn the_plausible_diagnostic_does_not_detect_collinearity() {
        // Why `collinearity` is s2/s1 and not s3/s2. On a cloud that is a
        // line with a whisker of noise, s3/s2 -- the ratio that reads like a
        // measure of flatness -- is a ratio between two quantities that are
        // both noise, and comes out looking excellent. Only s2/s1 sees that
        // there was never a second direction to fit a plane in.
        let mut points: Vec<Point3D> = (0..12)
            .map(|i| Point3D::from([i as f64, 2.0 * i as f64, 3.0 * i as f64]))
            .collect();
        for (i, point) in points.iter_mut().enumerate() {
            point.coords[2] += if i % 2 == 0 { 1e-5 } else { -1e-5 };
        }

        let fit = fit_plane(&points).expect("a fit");
        let [s1, s2, s3] = fit.singular_values;

        let plausible = -(s3 / s2).log10();
        assert!(
            plausible > DEFAULT_MAX_COLLINEARITY,
            "s3/s2 was supposed to look good here, and reads {plausible}"
        );

        assert!(fit.collinearity() > DEFAULT_MAX_COLLINEARITY, "s2/s1 missed it");
        assert_eq!(fit.quality(DEFAULT_MAX_COLLINEARITY), FitQuality::Collinear);
        assert!(s2 / s1 < 1e-5, "the cloud really is a line: s2/s1 = {}", s2 / s1);
    }

    #[test]
    fn scatter_about_a_plane_shows_up_in_the_residual() {
        let mut points = on_plane(90.0, 45.0, &SPREAD);
        for (i, point) in points.iter_mut().enumerate() {
            point.coords[2] += if i % 2 == 0 { 0.1 } else { -0.1 };
        }

        let fit = fit_plane(&points).expect("a fit");

        assert!(fit.rms_distance() > 0.01, "residual {}", fit.rms_distance());
        // Displacing along z alone leaves the fitted attitude close to the
        // truth, since the displacements cancel; the residual is what moves.
        assert!((fit.plane.dip_angle - 45.0).abs() < 5.0);
    }

    #[test]
    fn fewer_than_three_points_do_not_define_a_plane() {
        let points = on_plane(0.0, 30.0, &SPREAD[..2]);

        assert!(fit_plane(&points).is_none());
        assert!(fit_plane(&[]).is_none());
    }

    #[test]
    fn coincident_points_do_not_define_one_either() {
        let points = vec![Point3D::from([5.0, 5.0, 5.0]); 10];

        assert!(fit_plane(&points).is_none());
    }

    #[test]
    fn coincident_points_are_merged_once() {
        let mut points = on_plane(135.0, 35.0, &SPREAD);
        let duplicates = points.clone();
        points.extend(duplicates);

        let distinct = distinct_points(&points, 0.1);

        assert_eq!(distinct.len(), SPREAD.len());
    }

    #[test]
    fn points_further_apart_than_the_threshold_are_kept() {
        let points = vec![
            Point3D::from([0.0, 0.0, 0.0]),
            Point3D::from([0.05, 0.0, 0.0]),
            Point3D::from([0.5, 0.0, 0.0]),
        ];

        assert_eq!(distinct_points(&points, 0.1).len(), 2);
        assert_eq!(distinct_points(&points, 0.01).len(), 3);
    }

    #[test]
    fn merging_is_judged_in_three_dimensions() {
        // Same map position, ten metres apart vertically: the reference's 2D
        // test would merge these, and with them the only evidence of how
        // steeply the surface dips.
        let points = vec![
            Point3D::from([0.0, 0.0, 0.0]),
            Point3D::from([0.0, 0.0, 10.0]),
        ];

        assert_eq!(distinct_points(&points, 0.1).len(), 2);
    }

    #[test]
    fn a_field_recovers_one_attitude_in_every_cell() {
        // A plane sampled densely enough that each 50-unit cell holds a
        // properly two-dimensional patch of it.
        let mut xy = Vec::new();
        for i in 0..40 {
            for j in 0..40 {
                xy.push((i as f64 * 5.0, j as f64 * 5.0));
            }
        }
        let points = on_plane(135.0, 35.0, &xy);

        let field = best_fit_geoplanes(&points, 50.0, 0.1, DEFAULT_MAX_COLLINEARITY)
            .expect("a field");

        assert!(field.stats.cells_fitted > 10);
        assert_eq!(field.stats.cells_collinear, 0);

        for cell in &field.fits {
            assert!(
                (cell.fit.plane.azimuth - 135.0).abs() < 1e-6
                    && (cell.fit.plane.dip_angle - 35.0).abs() < 1e-6,
                "cell ({}, {}): {:?}",
                cell.row,
                cell.column,
                cell.fit.plane
            );
        }
    }

    #[test]
    fn a_trace_along_a_line_yields_no_attitudes() {
        // What a contact crossing a uniform slope looks like: the points in
        // any one cell are strung along a line, and the reference would post
        // an attitude for each of them.
        let xy: Vec<(f64, f64)> = (0..200).map(|i| (i as f64, i as f64)).collect();
        let points = on_plane(90.0, 20.0, &xy);

        let field = best_fit_geoplanes(&points, 20.0, 0.1, DEFAULT_MAX_COLLINEARITY)
            .expect("a field");

        assert!(field.fits.is_empty(), "collinear cells were fitted anyway");
        assert!(field.stats.cells_collinear > 0);
        assert_eq!(
            field.stats.cells_collinear + field.stats.cells_too_few_points,
            field.stats.non_empty_cells
        );
    }

    #[test]
    fn the_counts_raster_covers_the_whole_grid() {
        let points = on_plane(135.0, 35.0, &SPREAD);
        let field = best_fit_geoplanes(&points, 5.0, 0.1, DEFAULT_MAX_COLLINEARITY).expect("a field");

        assert_eq!(field.counts.len(), field.grid.rows * field.grid.columns);
        assert_eq!(field.counts.iter().sum::<usize>(), field.stats.distinct_points);
        assert_eq!(
            field.counts.iter().filter(|&&c| c > 0).count(),
            field.stats.non_empty_cells
        );
    }

    #[test]
    fn every_cell_is_accounted_for() {
        let mut xy = Vec::new();
        for i in 0..30 {
            for j in 0..3 {
                xy.push((i as f64 * 7.0, j as f64 * 7.0));
            }
        }
        let field = best_fit_geoplanes(&on_plane(45.0, 60.0, &xy), 25.0, 0.1, DEFAULT_MAX_COLLINEARITY)
            .expect("a field");

        assert_eq!(
            field.stats.cells_fitted
                + field.stats.cells_collinear
                + field.stats.cells_too_few_points,
            field.stats.non_empty_cells
        );
        assert_eq!(field.stats.cells_fitted, field.fits.len());
    }

    #[test]
    fn a_cell_centre_is_where_the_reference_puts_it() {
        // Anchored top-left at (x_min, y_max), so row 0 column 0 is centred
        // half a cell east and half a cell south of that corner.
        let points = vec![
            Point3D::from([100.0, 500.0, 0.0]),
            Point3D::from([200.0, 400.0, 0.0]),
            Point3D::from([150.0, 450.0, 1.0]),
        ];
        let grid = CellGrid::covering(&points, 50.0).expect("a grid");

        assert_eq!(grid.centre_of(0, 0), (125.0, 475.0));
        assert_eq!(grid.cell_of(&points[0]), (0, 0));
        assert_eq!(grid.rows, 3);
        assert_eq!(grid.columns, 3);
    }

    #[test]
    fn a_field_needs_a_usable_cell_size() {
        let points = on_plane(135.0, 35.0, &SPREAD);

        for cell_size in [0.0, -10.0, f64::NAN, f64::INFINITY] {
            assert!(best_fit_geoplanes(&points, cell_size, 0.1, 1.0).is_none(), "{cell_size}");
        }
        for coincidence in [0.0, -1.0, f64::NAN] {
            assert!(best_fit_geoplanes(&points, 50.0, coincidence, 1.0).is_none(), "{coincidence}");
        }
        assert!(best_fit_geoplanes(&[], 50.0, 0.1, 1.0).is_none());
    }

    #[test]
    fn the_fits_come_back_in_grid_order() {
        let mut xy = Vec::new();
        for i in 0..20 {
            for j in 0..20 {
                xy.push((i as f64 * 9.0, j as f64 * 9.0));
            }
        }
        let field = best_fit_geoplanes(&on_plane(200.0, 25.0, &xy), 40.0, 0.1, DEFAULT_MAX_COLLINEARITY)
            .expect("a field");

        let order: Vec<(usize, usize)> = field.fits.iter().map(|f| (f.row, f.column)).collect();
        let mut sorted = order.clone();
        sorted.sort();

        assert_eq!(order, sorted, "a HashMap's order leaked into the output");
    }
}
