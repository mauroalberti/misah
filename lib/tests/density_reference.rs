//! The density estimator against the program it descends from.
//!
//! `spatial::density` is a reimplementation of `InterpDensity3D`, a C++ program
//! that reads a list of points and writes a VTK volume of density. This test
//! runs the two on the same points with the same bandwidth and asks them to
//! agree, which covers three things the unit tests cannot reach on their own:
//! the normalizing constants in three dimensions, the registration of the grid,
//! and the order the nodes come out in.
//!
//! It is worth having because this lineage has been wrong about all three. The
//! C++ carried planar kernel constants in a volume, so its densities were not
//! per unit volume and two runs at different bandwidths were not comparable;
//! and it counted rows downwards from `y_max` while writing a header from
//! `y_min`, displacing the volume by 495 m on a 500 m cell. Both are fixed
//! there now, and the numbers here come from the fixed program -- which is the
//! only version this could be checked against, the older one being wrong in
//! ways an agreement would have propagated.
//!
//! ## Where the data come from
//!
//! `density_points.csv` is synthetic: two Gaussian clusters and a uniform
//! scatter, from a seeded generator, rounded to whole metres so that both
//! languages read exactly the same numbers. It is not a real catalogue.
//! `InterpDensity3D`'s own sample dataset would have been the obvious choice
//! and is deliberately not used: it is somebody's earthquake catalogue,
//! arriving here without a statement of its terms, and this repository's
//! `example_data/NOTICE.md` sets the rule that vendored data travels under
//! licence terms that can be named. These points are the repository's own.
//!
//! The extents are arranged so the grid comes out 8 by 6 by 10 -- three
//! different counts, so that a transposed traversal cannot pass by
//! coincidence.
//!
//! `density_interpdensity3d.vtk` is the C++ program's output on those points,
//! unedited, so nothing stands between the reference and the program that
//! produced it. To regenerate it, build `interp_dens.cpp` from
//! <https://gitlab.com/mauroalberti/InterpDensity3D> at 0.9 or later and feed
//! it a parameter file naming `density_points.csv`, a cell size of 1000, the
//! kernels `qn`, and a bandwidth of 2000.

use std::fs;
use std::path::PathBuf;

use misah::geometry::point::Point;
use misah::spatial::density::point_density_field;
use misah::spatial::grid::SamplingGrid;
use misah::spatial::kernel::{Bandwidth, Kernel};

/// The parameters the vendored reference was produced with.
const CELL_SIZE: f64 = 1000.0;
const BANDWIDTH: f64 = 2000.0;

/// `%.6e` carries seven significant digits, so two implementations that agree
/// exactly still differ by up to half a unit in the last one: five parts in
/// ten million. Anything under this means the reference file has no digit left
/// to disagree in.
const PRINTED_PRECISION: f64 = 1.0e-6;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data").join(name)
}

fn points() -> Vec<Point<3>> {

    fs::read_to_string(data("density_points.csv"))
        .expect("the vendored point set")
        .lines()
        .skip(1)
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            let mut fields = line.split(',').map(|f| {
                f.trim().parse::<f64>().expect("three numbers per line")
            });
            Point::from([
                fields.next().expect("an easting"),
                fields.next().expect("a northing"),
                fields.next().expect("an elevation"),
            ])
        })
        .collect()
}

/// The reference volume: its declared registration, and its two scalar blocks
/// in file order.
struct Reference {
    origin: Point<3>,
    spacing: [f64; 3],
    counts: [usize; 3],
    normal: Vec<f64>,
    quartic: Vec<f64>,
}

fn reference() -> Reference {

    let text = fs::read_to_string(data("density_interpdensity3d.vtk"))
        .expect("the vendored reference volume");

    let mut origin = None;
    let mut spacing = None;
    let mut counts = None;
    let mut blocks: Vec<Vec<f64>> = Vec::new();
    let mut in_values = false;

    let numbers = |line: &str| -> Vec<f64> {
        line.split_whitespace().skip(1).map(|f| f.parse().expect("a number")).collect()
    };

    for line in text.lines() {

        let line = line.trim();

        if let Some(rest) = line.strip_prefix("ORIGIN") {
            let v = numbers(&format!("_ {rest}"));
            origin = Some(Point::from([v[0], v[1], v[2]]));
        } else if line.starts_with("SPACING") {
            let v = numbers(line);
            spacing = Some([v[0], v[1], v[2]]);
        } else if line.starts_with("DIMENSIONS") {
            let v = numbers(line);
            counts = Some([v[0] as usize, v[1] as usize, v[2] as usize]);
        } else if line.starts_with("SCALARS") {
            blocks.push(Vec::new());
            in_values = false;
        } else if line.starts_with("LOOKUP_TABLE") {
            in_values = true;
        } else if in_values && !line.is_empty() {
            blocks.last_mut().expect("a scalar block").push(line.parse().expect("a value"));
        }
    }

    assert_eq!(blocks.len(), 2, "the reference should carry both kernels");
    let quartic = blocks.pop().expect("the quartic block");
    let normal = blocks.pop().expect("the normal block");

    Reference {
        origin: origin.expect("an ORIGIN"),
        spacing: spacing.expect("a SPACING"),
        counts: counts.expect("a DIMENSIONS"),
        normal,
        quartic,
    }
}

/// The grid the C++ sampled, rebuilt from the points alone.
///
/// Asserted against the reference's own header rather than read from it: this
/// is the registration check, and reading the answer would defeat it. The
/// program puts its first node half a cell in from the minimum of the data and
/// steps by whole cells from there, which is what `SamplingGrid` means by an
/// origin.
fn grid_from(points: &[Point<3>]) -> SamplingGrid<3> {

    let mut lower = [f64::INFINITY; 3];
    let mut upper = [f64::NEG_INFINITY; 3];

    for point in points {
        for axis in 0..3 {
            lower[axis] = lower[axis].min(point.coords[axis]);
            upper[axis] = upper[axis].max(point.coords[axis]);
        }
    }

    let origin = Point::from(std::array::from_fn(|axis| lower[axis] + CELL_SIZE / 2.0));
    // Truncating division, one more node than whole cells. The C++ takes the
    // range through a single-precision float first, which can shift a count by
    // one where the extent falls within a float's resolution of a multiple of
    // the cell size; that is a wart of the C++ rather than a convention to
    // reproduce, and this dataset sits nowhere near such a boundary.
    let counts = std::array::from_fn(|axis| {
        ((upper[axis] - lower[axis]) / CELL_SIZE) as usize + 1
    });

    SamplingGrid::new(origin, [CELL_SIZE; 3], counts).expect("a well-formed grid")
}

/// Every node compared, with the zeros required to be exactly zero.
fn agrees(name: &str, mine: &[f64], theirs: &[f64]) {

    assert_eq!(mine.len(), theirs.len(), "{name}: different node counts");

    let mut worst = 0.0f64;
    let mut worst_at = 0usize;

    for (node, (&ours, &reference)) in mine.iter().zip(theirs).enumerate() {

        if reference == 0.0 {
            // A compact kernel is zero beyond its reach, and the boundary is
            // part of what is being checked: a value that ought to be nothing
            // must be nothing, not merely small.
            assert_eq!(
                ours, 0.0,
                "{name}: node {node} is zero in the reference and {ours:e} here"
            );
            continue;
        }

        let relative = (ours - reference).abs() / reference;
        if relative > worst {
            worst = relative;
            worst_at = node;
        }
    }

    assert!(
        worst < PRINTED_PRECISION,
        "{name}: worst relative difference {worst:e} at node {worst_at} \
         ({} against {}), past the {PRINTED_PRECISION:e} the reference is printed to",
        mine[worst_at],
        theirs[worst_at]
    );
}

#[test]
fn the_density_matches_interpdensity3d_on_the_same_points() {

    let points = points();
    let reference = reference();
    let grid = grid_from(&points);

    // The registration, before any density is compared: two volumes that
    // disagree about where their nodes are would still both be smooth.
    assert_eq!(grid.origin(), &reference.origin, "the grid origins differ");
    assert_eq!(grid.spacing(), &reference.spacing, "the grid spacings differ");
    assert_eq!(grid.counts(), &reference.counts, "the grid shapes differ");
    assert_eq!(grid.node_count(), reference.normal.len());

    let bandwidth = Bandwidth::isotropic(BANDWIDTH).expect("a positive bandwidth");

    // Untruncated, because the C++ sums every point at every node and a
    // truncation here would be a real difference rather than a rounding one.
    agrees(
        "normal",
        &point_density_field(&points, &grid, Kernel::gaussian(bandwidth)),
        &reference.normal,
    );
    agrees(
        "quartic",
        &point_density_field(&points, &grid, Kernel::quartic(bandwidth)),
        &reference.quartic,
    );
}

#[test]
fn truncating_the_gaussian_removes_a_bounded_amount_and_nothing_else() {

    // The kernel a field will actually be run with, held against the one just
    // checked against the C++.
    //
    // What truncation costs is an *absolute* quantity, and it is worth being
    // clear that it is not the tail fraction of the kernel's mass. That
    // fraction -- three per cent at three bandwidths in three dimensions --
    // says what is lost from the integral over the whole volume. It says
    // nothing about the loss at any one node, because a node far from every
    // point holds a density made entirely of distant tails, and cutting the
    // tails takes essentially all of it. On this dataset, cutting at three
    // bandwidths moves the worst node by twenty-nine per cent of its own
    // value.
    //
    // The bound that does hold is per point and unconditional. The profile
    // falls off monotonically, so a point beyond the cut was contributing less
    // than the profile at the cut; with `n` points, nothing anywhere can move
    // by more than `n` times that.
    let points = points();
    let grid = grid_from(&points);
    let bandwidth = Bandwidth::isotropic(BANDWIDTH).expect("a positive bandwidth");
    let kernel = Kernel::gaussian(bandwidth);

    let origin = Point::from([0.0, 0.0, 0.0]);
    let peak_of_one_kernel = kernel.weight(&origin, &origin);

    let exact = point_density_field(&points, &grid, kernel);
    let peak_density = exact.iter().cloned().fold(0.0f64, f64::max);

    for radii in [3.0, 4.0, 6.0] {

        let cut = Kernel::truncated_gaussian(bandwidth, radii).expect("a positive cut");
        let truncated = point_density_field(&points, &grid, cut);

        let bound = points.len() as f64 * peak_of_one_kernel * (-0.5 * radii * radii).exp();

        for (node, (&whole, &part)) in exact.iter().zip(&truncated).enumerate() {
            // Truncation only ever removes -- to within the rounding of the
            // sum, which is not a formality here. The two fields add the same
            // positive terms in different orders: the untruncated kernel has
            // no reach, so it walks the points as given, while the truncated
            // one walks them bin by bin. Summing `n` positive terms carries an
            // error of about `n` epsilons of the total either way, and at a
            // node where the truncation removes nothing at all that rounding
            // is the entire difference and can fall on either side of zero.
            let rounding = 4.0 * points.len() as f64 * f64::EPSILON * whole;
            assert!(
                part <= whole + rounding,
                "cutting at {radii} bandwidths *added* {:e} at node {node}, past rounding",
                part - whole
            );
            assert!(
                whole - part <= bound + rounding,
                "cutting at {radii} bandwidths cost {:e} at node {node}, over the bound {bound:e}",
                whole - part
            );
        }

        // And the consequence worth stating: where a density map is actually
        // read -- the nodes carrying a tenth of the peak or more -- that same
        // bound is a small fraction of the value.
        let worst_where_it_counts = exact
            .iter()
            .zip(&truncated)
            .filter(|(whole, _)| **whole >= 0.1 * peak_density)
            .map(|(whole, part)| (whole - part).max(0.0) / whole)
            .fold(0.0f64, f64::max);

        assert!(
            worst_where_it_counts <= bound / (0.1 * peak_density),
            "the derived bound does not hold: {worst_where_it_counts:e}"
        );

        if radii >= 4.0 {
            assert!(
                worst_where_it_counts < 1.0e-3,
                "cutting at {radii} bandwidths moved a substantial node by \
                 {worst_where_it_counts:e}"
            );
        }
    }
}
