//! Python bindings for the fields: density on a grid, and a stress tensor at
//! every node of one.
//!
//! Kept apart from `kernels.rs` for length rather than for kind. From Python
//! these register into the same `misah.kernels` module as everything else, and
//! there is no second import path to learn.
//!
//! ## The dimension is a run-time value here and a type over there
//!
//! The Rust side is generic over `N`, the number of coordinates that locate an
//! observation, and `N` is a const parameter -- one compiled copy per
//! dimension, chosen when the code is written. Python has no such thing: the
//! dimension arrives as the width of an array, at run time, and something has
//! to bridge the two. That is all `at_dimension!` does, and it is why a field
//! here is one, two or three dimensional and not any number: each is a
//! monomorphisation somebody has to name.
//!
//! One dimension is not padding for the sake of a tidy range. A histogram of
//! hypocentral depths is a one-dimensional density, and estimating it with a
//! kernel rather than with bins is the same improvement here as anywhere --
//! no edges chosen by hand, and no answer that changes when they move.
//!
//! ## Flat order, and the reshape
//!
//! Every field comes back as a flat array in the grid's own order, **first
//! axis fastest**: x runs, then y, then z. Not reshaped, because the reshape
//! is where the convention would go unsaid. To lay one out as an array, ask
//! numpy for Fortran order explicitly:
//!
//! ```python
//! volume = values.reshape(stats["grid_counts"], order="F")   # volume[i, j, k]
//! ```
//!
//! The same order a VTK `STRUCTURED_POINTS` file wants, which is not a
//! coincidence: it is the order the grid was defined in so that writing one
//! needs no transposition.

use numpy::{IntoPyArray, PyArray1, PyArrayMethods, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use misah::geometry::located::Located;
use misah::geometry::point::Point;
use misah::spatial::density::density_field as density_kernel;
use misah::spatial::grid::SamplingGrid;
use misah::spatial::kernel::{Bandwidth, Kernel};
use misah::structural::fault::FaultPlane;
use misah::structural::geol_axis::GeologicalAxis;
use misah::structural::inversion::SearchGrid;
use misah::structural::stress_field::{
    FieldNode, field_cost as field_cost_kernel, stress_field as stress_field_kernel,
};

use crate::kernels::{checked_grid, faults_from_rows};

/// Call a function that is generic over `N` at the dimension a call turned out
/// to carry.
///
/// `N` is a const generic, so there is one compiled copy per dimension and the
/// choice cannot be made by any ordinary run-time means -- something has to
/// name each copy. This names them once rather than at each entry point, and
/// gives the same refusal for a fourth dimension everywhere.
macro_rules! at_dimension {
    ($dimension:expr, $call:ident ( $($arg:expr),* $(,)? )) => {
        match $dimension {
            1 => $call::<1>($($arg),*),
            2 => $call::<2>($($arg),*),
            3 => $call::<3>($($arg),*),
            other => Err(PyValueError::new_err(format!(
                "a field must have 1, 2 or 3 coordinates per observation, not {other}"
            ))),
        }
    };
}

/// A sequence of the length the dimension calls for, as a fixed array.
///
/// Named so the refusal can say which argument was the wrong length, which is
/// the whole difficulty of an interface where five arguments must agree on a
/// number none of them states.
fn fixed<const N: usize>(values: &[f64], name: &str) -> PyResult<[f64; N]> {

    if values.len() != N {
        return Err(PyValueError::new_err(format!(
            "{name} has {} entries for a {N}-dimensional field",
            values.len()
        )));
    }

    Ok(std::array::from_fn(|i| values[i]))
}

/// The same for the node counts, which are whole numbers and must be positive.
fn fixed_counts<const N: usize>(values: &[usize]) -> PyResult<[usize; N]> {

    if values.len() != N {
        return Err(PyValueError::new_err(format!(
            "counts has {} entries for a {N}-dimensional field",
            values.len()
        )));
    }

    Ok(std::array::from_fn(|i| values[i]))
}

/// The kernel a name and a bandwidth ask for.
///
/// The default elsewhere is `quartic`, and deliberately: it is compact by
/// construction, which is what keeps a node's cost proportional to its own
/// neighbourhood. A plain `gaussian` reaches everywhere and makes every node
/// cost the whole dataset -- correct, sometimes wanted, and never the thing to
/// reach for first.
fn kernel_of<const N: usize>(
    bandwidth: [f64; N],
    shape: &str,
    truncation: f64,
) -> PyResult<Kernel<N>> {

    let bandwidth = Bandwidth::new(bandwidth).ok_or_else(|| {
        PyValueError::new_err(
            "bandwidth must be positive and finite on every axis: it is a width in map units, \
             and a zero one is a kernel of no width and infinite height",
        )
    })?;

    match shape {
        "quartic" => Ok(Kernel::quartic(bandwidth)),
        "gaussian" => Ok(Kernel::gaussian(bandwidth)),
        "truncated_gaussian" => Kernel::truncated_gaussian(bandwidth, truncation).ok_or_else(|| {
            PyValueError::new_err(format!(
                "truncation must be a positive finite number of bandwidths, not {truncation}"
            ))
        }),
        other => Err(PyValueError::new_err(format!(
            "unknown kernel {other:?}: expected 'quartic', 'gaussian' or 'truncated_gaussian'"
        ))),
    }
}

/// The grid, once its parts are known to agree on a dimension.
fn grid_of<const N: usize>(
    origin: &[f64],
    spacing: &[f64],
    counts: &[usize],
) -> PyResult<SamplingGrid<N>> {

    let origin = fixed::<N>(origin, "origin")?;
    let spacing = fixed::<N>(spacing, "spacing")?;
    let counts = fixed_counts::<N>(counts)?;

    SamplingGrid::new(Point::from(origin), spacing, counts).ok_or_else(|| {
        PyValueError::new_err(
            "a grid needs a finite origin, a positive finite spacing on every axis, at least one \
             node on each, and a node count that fits in memory",
        )
    })
}

/// The rows of an (M, N) array as points.
fn positions_of<const N: usize>(view: &ndarray::ArrayView2<'_, f64>) -> Vec<Point<N>> {

    view.rows()
        .into_iter()
        .map(|row| Point::from(std::array::from_fn(|i| row[i])))
        .collect()
}

/// Trend and plunge of an axis, the pair every orientation crosses this
/// boundary as.
fn trend_plunge(versor: &misah::algebra::versor::Versor3D) -> (f64, f64) {
    let axis = GeologicalAxis::from_versor(versor);
    (axis.trend, axis.plunge)
}

/// A grid covering a set of points, with room around them.
///
/// The first step of nearly every field, and worth having rather than done by
/// hand in numpy, because it settles the one question a raster header cannot:
/// `origin` is the **first sample point**, not the corner of a cell around it.
/// Half a spacing either way is the classic silent error in this lineage, and
/// there is nothing here to be half of.
///
/// `margin` extends the grid past the outermost observation on every side, in
/// map units. About one bandwidth is the useful choice: it is where a compact
/// kernel has finished falling away, so the field shows the whole of its own
/// decay rather than being cut off mid-slope.
///
/// Returns a dict of `origin`, `spacing`, `counts` and `node_count`, which is
/// exactly what `density_field` and `stress_field` take. `counts` is the
/// number of nodes along each axis, not the number of intervals.
#[pyfunction]
#[pyo3(signature = (points, spacing, margin = 0.0))]
fn covering_grid<'py>(
    py: Python<'py>,
    points: PyReadonlyArray2<'py, f64>,
    spacing: Vec<f64>,
    margin: f64,
) -> PyResult<Bound<'py, PyDict>> {

    let view = points.as_array();
    let dimension = view.ncols();

    if view.nrows() == 0 {
        return Err(PyValueError::new_err("no points to cover"));
    }
    if spacing.len() != dimension {
        return Err(PyValueError::new_err(format!(
            "spacing has {} entries for points with {dimension} coordinates",
            spacing.len()
        )));
    }

    fn build<const N: usize>(
        view: &ndarray::ArrayView2<'_, f64>,
        spacing: &[f64],
        margin: f64,
    ) -> PyResult<([f64; N], [f64; N], [usize; N])> {

        let spacing = fixed::<N>(spacing, "spacing")?;
        let points = positions_of::<N>(view);

        let grid = SamplingGrid::covering(&points, spacing, margin).ok_or_else(|| {
            PyValueError::new_err(
                "could not cover those points: spacing and margin must be positive and finite, \
                 and every coordinate must be a number",
            )
        })?;

        Ok((grid.origin().coords, *grid.spacing(), *grid.counts()))
    }

    let dict = PyDict::new(py);

    // Written out per dimension rather than through `at_dimension!`, the
    // return type being a different tuple in each arm.
    match dimension {
        1 => {
            let (origin, spacing, counts) = build::<1>(&view, &spacing, margin)?;
            dict.set_item("origin", origin.to_vec())?;
            dict.set_item("spacing", spacing.to_vec())?;
            dict.set_item("counts", counts.to_vec())?;
            dict.set_item("node_count", counts.iter().product::<usize>())?;
        }
        2 => {
            let (origin, spacing, counts) = build::<2>(&view, &spacing, margin)?;
            dict.set_item("origin", origin.to_vec())?;
            dict.set_item("spacing", spacing.to_vec())?;
            dict.set_item("counts", counts.to_vec())?;
            dict.set_item("node_count", counts.iter().product::<usize>())?;
        }
        3 => {
            let (origin, spacing, counts) = build::<3>(&view, &spacing, margin)?;
            dict.set_item("origin", origin.to_vec())?;
            dict.set_item("spacing", spacing.to_vec())?;
            dict.set_item("counts", counts.to_vec())?;
            dict.set_item("node_count", counts.iter().product::<usize>())?;
        }
        other => {
            return Err(PyValueError::new_err(format!(
                "a field must have 1, 2 or 3 coordinates per observation, not {other}"
            )));
        }
    }

    Ok(dict)
}

/// What a kernel does, before any data are put through it.
///
/// Three numbers worth knowing before a field is run rather than after.
/// `reach` is how far one observation is felt along each axis, or `None` for
/// the untruncated Gaussian, which never quite stops -- and a `None` here is
/// the warning that every node will cost the whole dataset. `peak` is the
/// density a single observation produces at its own position, which is the
/// scale a colour ramp should be read against. `retained_mass` is the fraction
/// of the kernel that survives its own truncation: one for the quartic and the
/// plain Gaussian, and about 0.999 for a Gaussian cut at four bandwidths.
///
/// The missing mass is quoted and not corrected for. A density short by a
/// known thousandth is easier to reason about than one silently rescaled, and
/// a caller who wants it back can divide.
#[pyfunction]
#[pyo3(signature = (bandwidth, kernel = "quartic", truncation = 4.0))]
fn kernel_profile<'py>(
    py: Python<'py>,
    bandwidth: Vec<f64>,
    kernel: &str,
    truncation: f64,
) -> PyResult<Bound<'py, PyDict>> {

    fn describe<const N: usize>(
        bandwidth: &[f64],
        shape: &str,
        truncation: f64,
    ) -> PyResult<(Option<Vec<f64>>, f64, f64)> {

        let extents = fixed::<N>(bandwidth, "bandwidth")?;
        let kernel = kernel_of::<N>(extents, shape, truncation)?;

        let origin = Point::<N>::from([0.0; N]);

        Ok((
            kernel.reach().map(|r| r.to_vec()),
            kernel.weight(&origin, &origin),
            kernel.retained_mass(),
        ))
    }

    let (reach, peak, retained) =
        at_dimension!(bandwidth.len(), describe(&bandwidth, kernel, truncation))?;

    let dict = PyDict::new(py);
    dict.set_item("reach", reach)?;
    dict.set_item("peak", peak)?;
    dict.set_item("retained_mass", retained)?;

    Ok(dict)
}

/// Kernel density over located points, at every node of a grid.
///
/// `points` is an (M, D) array of coordinates with D of 1, 2 or 3, and
/// `origin`, `spacing`, `counts` and `bandwidth` are each D long -- the same
/// grid `covering_grid` returns, which is where those three usually come from.
///
/// The value at a node is **observations per unit volume** there, per unit area
/// in two dimensions, per unit length in one. So the field summed over the grid
/// and multiplied by the cell volume approaches the number of observations,
/// less whatever fell outside the grid and whatever a truncated kernel dropped.
/// That is the property the normalization exists for, and it is what makes two
/// runs at different bandwidths comparable with each other.
///
/// `bandwidth` is in **map units, one per axis**. Not in cells: a bandwidth in
/// cells changes meaning when the grid is refined, which is exactly when a
/// reader is looking for the answer to stop changing. One per axis, because
/// depth is not interchangeable with easting even when both are in metres --
/// a seismogenic layer is far wider than it is thick, and an isotropic kernel
/// over one mixes its top with its bottom before it mixes two neighbours.
///
/// Returns a flat (M,) array in the grid's order, first axis fastest. To lay it
/// out: `values.reshape(counts, order="F")`.
///
/// Runs across every core rayon is given, with the interpreter released. The
/// result does not depend on how many cores that was, bit for bit.
#[pyfunction]
#[pyo3(signature = (
    points,
    origin,
    spacing,
    counts,
    bandwidth,
    kernel = "quartic",
    truncation = 4.0,
))]
#[allow(clippy::too_many_arguments)]
fn density_field<'py>(
    py: Python<'py>,
    points: PyReadonlyArray2<'py, f64>,
    origin: Vec<f64>,
    spacing: Vec<f64>,
    counts: Vec<usize>,
    bandwidth: Vec<f64>,
    kernel: &str,
    truncation: f64,
) -> PyResult<Bound<'py, PyArray1<f64>>> {

    let view = points.as_array();
    let dimension = view.ncols();

    fn run<const N: usize>(
        py: Python<'_>,
        view: &ndarray::ArrayView2<'_, f64>,
        origin: &[f64],
        spacing: &[f64],
        counts: &[usize],
        bandwidth: &[f64],
        shape: &str,
        truncation: f64,
    ) -> PyResult<Vec<f64>> {

        let grid = grid_of::<N>(origin, spacing, counts)?;
        let kernel = kernel_of::<N>(fixed::<N>(bandwidth, "bandwidth")?, shape, truncation)?;

        let observations: Vec<Located<(), N>> = positions_of::<N>(view)
            .into_iter()
            .map(|p| Located::new(p, ()))
            .collect();

        // The whole cost of the call is here, and it touches nothing Python
        // owns, so the interpreter is let go for the duration.
        Ok(py.detach(|| density_kernel(&observations, &grid, kernel)))
    }

    let values = at_dimension!(
        dimension,
        run(py, &view, &origin, &spacing, &counts, &bandwidth, kernel, truncation)
    )?;

    Ok(values.into_pyarray(py))
}

/// A stress tensor at every node of a grid, with the density of the data that
/// produced it.
///
/// Hardebeck and Michael's (2006) spatially varying inversion in kernel form.
/// Rather than cutting the dataset into bins and damping neighbouring solutions
/// towards each other, every fault contributes to every node it can reach, by
/// an amount that falls off with distance. Where the data are dense the
/// estimate is local; where they thin out the same kernel reaches further into
/// what there is, and `support` says so.
///
/// `positions` is an (M, D) array saying where each fault was found, and
/// `faults` the (M, 4) array `invert_stress` takes -- right-hand-rule strike,
/// dip, then the trend and plunge of the slickenline. `senses` is optional and
/// works exactly as it does there. `positions` says where; `faults` says what.
/// D is the dimension of the **field**, not of the geology: a fault plane is
/// three-dimensional whether it is posted on a map or in a volume.
///
/// `min_support` is the effective sample size below which a node is left
/// without a tensor. Zero inverts wherever a single fault reaches, which is
/// almost never wanted: four parameters are being estimated, so a node held up
/// by three faults' worth of weight has a solution in the sense that the search
/// returns one and in no other sense. Ten is a defensible floor for real data,
/// and is the default here for that reason.
///
/// **Call `field_cost` first.** It walks the same nodes with the same kernel in
/// a fraction of a second and reports the work this would do. A node costs a
/// whole exhaustive search -- of the order of a tenth of a second, times
/// however many faults reach it -- against microseconds for the density at that
/// same node. A density grid worth looking at is finer than a stress grid can
/// afford by several orders of magnitude, which is why `density_field` stays
/// available on its own.
///
/// Returns `(field, stats)`. `field` is a dict of arrays, all M long and
/// aligned with each other and with the grid's flat order:
///
/// - `positions` (M, D), the node coordinates
/// - `density` (M,), faults per unit area or volume, from the same kernel
/// - `support` (M,), Kish's `(sum w)^2 / sum w^2` -- the number to read beside
///   a tensor. Twenty faults at equal weight give twenty; twenty of which
///   nineteen sit on the far edge of the kernel give barely one
/// - `faults_within_reach` (M,), how many reached at all, at any weight
/// - `has_solution` (M,) boolean
/// - `s1`, `s2`, `s3` (M, 2) as trend and plunge, `phi` (M,),
///   `misfit_degrees` (M,), `faults_scored` (M,) and
///   `effective_sample_size` (M,)
/// - `runner_up_misfit_degrees` (M,), the best runner-up's misfit, and
///   `runner_up_s1_degrees` (M,), the largest angle between the winning S1 and
///   any runner-up's
///
/// Everything from the solution is `NaN`, or zero for the counts, at a node
/// without one. Read `has_solution` rather than testing for `NaN`.
///
/// That last pair is what says whether a minimum stood by itself. On a grid the
/// runners-up are usually the winner's own neighbours, a step away in S1, and
/// `runner_up_s1_degrees` stays small; a runner-up fitting nearly as well from
/// right across the sphere is two tectonic phases meeting at that node, and it
/// shows up as a large angle at a low misfit. A stress map read without it will
/// show a smooth rotation through a region where the data are in fact saying
/// two different things.
///
/// `stats` carries `nodes`, `nodes_inverted`, `nodes_below_threshold`,
/// `nodes_without_data`, and the grid back as `grid_origin`, `grid_spacing` and
/// `grid_counts` for the reshape.
///
/// Nodes with no fault in reach are not an error. A field is entitled to holes,
/// and they are counted rather than filled -- as are nodes whose support fell
/// short, which keep their density, so a thin patch stays legible instead of
/// carrying an answer the search will always produce.
#[pyfunction]
#[pyo3(signature = (
    positions,
    faults,
    origin,
    spacing,
    counts,
    bandwidth,
    senses = None,
    kernel = "quartic",
    truncation = 4.0,
    min_support = 10.0,
    angle_step_degrees = 10.0,
    phi_step = 0.1,
))]
#[allow(clippy::too_many_arguments)]
fn stress_field<'py>(
    py: Python<'py>,
    positions: PyReadonlyArray2<'py, f64>,
    faults: PyReadonlyArray2<'py, f64>,
    origin: Vec<f64>,
    spacing: Vec<f64>,
    counts: Vec<usize>,
    bandwidth: Vec<f64>,
    senses: Option<PyReadonlyArray1<'py, bool>>,
    kernel: &str,
    truncation: f64,
    min_support: f64,
    angle_step_degrees: f64,
    phi_step: f64,
) -> PyResult<(Bound<'py, PyDict>, Bound<'py, PyDict>)> {

    let search = checked_grid(angle_step_degrees, phi_step)?;
    let planes = faults_from_rows(&faults, senses.as_ref())?;

    let view = positions.as_array();
    let dimension = view.ncols();

    if view.nrows() != planes.len() {
        return Err(PyValueError::new_err(format!(
            "positions has {} rows but {} faults were given",
            view.nrows(),
            planes.len()
        )));
    }

    // Not folded into the kernel's own refusal: `min_support` is compared with
    // `<`, so a NaN would leave every node above the threshold and invert the
    // whole grid at whatever a single fault happens to say. Silently, and
    // expensively.
    if min_support.is_nan() {
        return Err(PyValueError::new_err("min_support must be a number"));
    }

    #[allow(clippy::type_complexity)]
    fn run<const N: usize>(
        py: Python<'_>,
        view: &ndarray::ArrayView2<'_, f64>,
        planes: Vec<FaultPlane>,
        origin: &[f64],
        spacing: &[f64],
        counts: &[usize],
        bandwidth: &[f64],
        shape: &str,
        truncation: f64,
        min_support: f64,
        search: SearchGrid,
    ) -> PyResult<(Vec<FieldNode<N>>, [usize; 3])> {

        let grid = grid_of::<N>(origin, spacing, counts)?;
        let kernel = kernel_of::<N>(fixed::<N>(bandwidth, "bandwidth")?, shape, truncation)?;

        let located: Vec<Located<FaultPlane, N>> = positions_of::<N>(view)
            .into_iter()
            .zip(planes)
            .map(|(position, fault)| Located::new(position, fault))
            .collect();

        let field =
            py.detach(|| stress_field_kernel(&located, &grid, kernel, search, min_support));

        let counters = [
            field.nodes_inverted,
            field.nodes_below_threshold,
            field.nodes_without_data,
        ];

        Ok((field.nodes, counters))
    }

    // Each arm builds its own dicts, the node type differing by dimension.
    macro_rules! collect_at {
        ($n:literal) => {{
            let (nodes, counters) = run::<$n>(
                py,
                &view,
                planes,
                &origin,
                &spacing,
                &counts,
                &bandwidth,
                kernel,
                truncation,
                min_support,
                search,
            )?;
            (field_dict::<$n>(py, &nodes)?, counters, nodes.len())
        }};
    }

    let (field, counters, node_count) = match dimension {
        1 => collect_at!(1),
        2 => collect_at!(2),
        3 => collect_at!(3),
        other => {
            return Err(PyValueError::new_err(format!(
                "a field must have 1, 2 or 3 coordinates per observation, not {other}"
            )));
        }
    };

    let stats = PyDict::new(py);
    stats.set_item("nodes", node_count)?;
    stats.set_item("nodes_inverted", counters[0])?;
    stats.set_item("nodes_below_threshold", counters[1])?;
    stats.set_item("nodes_without_data", counters[2])?;
    stats.set_item("grid_origin", origin)?;
    stats.set_item("grid_spacing", spacing)?;
    stats.set_item("grid_counts", counts)?;

    Ok((field, stats))
}

/// The per-node arrays a field comes back as.
///
/// One pass over the nodes filling every column, rather than one pass per
/// column: the columns have to stay aligned with each other and with the
/// grid's flat order, and building them together is what makes that structural
/// instead of a thing to remember.
fn field_dict<'py, const N: usize>(
    py: Python<'py>,
    nodes: &[FieldNode<N>],
) -> PyResult<Bound<'py, PyDict>> {

    let count = nodes.len();

    let mut positions = Vec::with_capacity(count * N);
    let mut density = Vec::with_capacity(count);
    let mut support = Vec::with_capacity(count);
    let mut within_reach = Vec::with_capacity(count);
    let mut has_solution = Vec::with_capacity(count);
    let mut s1 = Vec::with_capacity(count * 2);
    let mut s2 = Vec::with_capacity(count * 2);
    let mut s3 = Vec::with_capacity(count * 2);
    let mut phi = Vec::with_capacity(count);
    let mut misfit = Vec::with_capacity(count);
    let mut scored = Vec::with_capacity(count);
    let mut effective = Vec::with_capacity(count);
    let mut runner_up_misfit = Vec::with_capacity(count);
    let mut runner_up_angle = Vec::with_capacity(count);

    for node in nodes {

        positions.extend_from_slice(&node.position.coords);
        density.push(node.density);
        support.push(node.support);
        within_reach.push(node.faults_within_reach as i64);

        let Some(solution) = &node.solution else {

            has_solution.push(false);
            // NaN and not zero: a node without a tensor has no S1 at north and
            // horizontal, no misfit of nothing, and no shape ratio of zero. A
            // zero would plot.
            s1.extend_from_slice(&[f64::NAN; 2]);
            s2.extend_from_slice(&[f64::NAN; 2]);
            s3.extend_from_slice(&[f64::NAN; 2]);
            phi.push(f64::NAN);
            misfit.push(f64::NAN);
            scored.push(0i64);
            effective.push(f64::NAN);
            runner_up_misfit.push(f64::NAN);
            runner_up_angle.push(f64::NAN);

            continue;
        };

        let best = &solution.best;
        let best_s1 = best.tensor.s1_versor();

        has_solution.push(true);
        let (trend, plunge) = trend_plunge(&best_s1);
        s1.extend_from_slice(&[trend, plunge]);
        let (trend, plunge) = trend_plunge(&best.tensor.s2_versor());
        s2.extend_from_slice(&[trend, plunge]);
        let (trend, plunge) = trend_plunge(&best.tensor.s3_versor());
        s3.extend_from_slice(&[trend, plunge]);
        phi.push(best.tensor.phi);
        misfit.push(best.mean_misfit_degrees);
        scored.push(best.faults_scored as i64);
        effective.push(best.effective_sample_size);

        runner_up_misfit.push(
            solution
                .runners_up
                .first()
                .map(|r| r.mean_misfit_degrees)
                .unwrap_or(f64::NAN),
        );

        // The furthest runner-up rather than the nearest: the question is
        // whether anything fitting nearly as well sits somewhere else
        // entirely, and one that does is not made harmless by four others
        // clustered on the winner. Axes are undirected, so the angle runs from
        // zero to ninety and `abs` on the dot product is what makes it so.
        runner_up_angle.push(
            solution
                .runners_up
                .iter()
                .map(|r| {
                    best_s1
                        .dot(&r.tensor.s1_versor())
                        .abs()
                        .clamp(-1.0, 1.0)
                        .acos()
                        .to_degrees()
                })
                .fold(f64::NAN, f64::max),
        );
    }

    let dict = PyDict::new(py);

    dict.set_item("positions", positions.into_pyarray(py).reshape([count, N])?)?;
    dict.set_item("density", density.into_pyarray(py))?;
    dict.set_item("support", support.into_pyarray(py))?;
    dict.set_item("faults_within_reach", within_reach.into_pyarray(py))?;
    dict.set_item("has_solution", has_solution.into_pyarray(py))?;
    dict.set_item("s1", s1.into_pyarray(py).reshape([count, 2])?)?;
    dict.set_item("s2", s2.into_pyarray(py).reshape([count, 2])?)?;
    dict.set_item("s3", s3.into_pyarray(py).reshape([count, 2])?)?;
    dict.set_item("phi", phi.into_pyarray(py))?;
    dict.set_item("misfit_degrees", misfit.into_pyarray(py))?;
    dict.set_item("faults_scored", scored.into_pyarray(py))?;
    dict.set_item("effective_sample_size", effective.into_pyarray(py))?;
    dict.set_item("runner_up_misfit_degrees", runner_up_misfit.into_pyarray(py))?;
    dict.set_item("runner_up_s1_degrees", runner_up_angle.into_pyarray(py))?;

    Ok(dict)
}

/// What `stress_field` would cost, measured rather than guessed.
///
/// Takes exactly what `stress_field` takes and runs the cheap half of it: the
/// neighbour search and the support at every node. Not an estimate from a
/// formula -- it visits the same nodes with the same kernel and counts the
/// faults that actually reach each, so a dataset clustered in one corner of its
/// own grid is costed as such rather than as if it were spread evenly.
///
/// Returns `nodes`, `nodes_to_invert` and `forward_solutions`. That last is the
/// unit the problem is made of: one candidate tensor resolved onto one fault,
/// summed over the nodes that will be inverted. Divide it by the rate the
/// machine sustains -- of the order of a few million a second per core -- for
/// the answer to "how long".
///
/// The number is usually far smaller than a first guess, because a compact
/// kernel means a node sees its own neighbourhood and not the dataset. In the
/// measurement this library's defaults were tuned against, each node saw about
/// 22 of 300 faults, so a hundred-node field cost four times a single inversion
/// over the whole set rather than a hundred times it.
#[pyfunction]
#[pyo3(signature = (
    positions,
    faults,
    origin,
    spacing,
    counts,
    bandwidth,
    senses = None,
    kernel = "quartic",
    truncation = 4.0,
    min_support = 10.0,
    angle_step_degrees = 10.0,
    phi_step = 0.1,
))]
#[allow(clippy::too_many_arguments)]
fn field_cost<'py>(
    py: Python<'py>,
    positions: PyReadonlyArray2<'py, f64>,
    faults: PyReadonlyArray2<'py, f64>,
    origin: Vec<f64>,
    spacing: Vec<f64>,
    counts: Vec<usize>,
    bandwidth: Vec<f64>,
    senses: Option<PyReadonlyArray1<'py, bool>>,
    kernel: &str,
    truncation: f64,
    min_support: f64,
    angle_step_degrees: f64,
    phi_step: f64,
) -> PyResult<Bound<'py, PyDict>> {

    let search = checked_grid(angle_step_degrees, phi_step)?;
    let planes = faults_from_rows(&faults, senses.as_ref())?;

    let view = positions.as_array();
    let dimension = view.ncols();

    if view.nrows() != planes.len() {
        return Err(PyValueError::new_err(format!(
            "positions has {} rows but {} faults were given",
            view.nrows(),
            planes.len()
        )));
    }
    if min_support.is_nan() {
        return Err(PyValueError::new_err("min_support must be a number"));
    }

    #[allow(clippy::too_many_arguments)]
    fn run<const N: usize>(
        py: Python<'_>,
        view: &ndarray::ArrayView2<'_, f64>,
        planes: &[FaultPlane],
        origin: &[f64],
        spacing: &[f64],
        counts: &[usize],
        bandwidth: &[f64],
        shape: &str,
        truncation: f64,
        min_support: f64,
        search: SearchGrid,
    ) -> PyResult<(usize, usize, u128)> {

        let grid = grid_of::<N>(origin, spacing, counts)?;
        let kernel = kernel_of::<N>(fixed::<N>(bandwidth, "bandwidth")?, shape, truncation)?;

        let located: Vec<Located<FaultPlane, N>> = positions_of::<N>(view)
            .into_iter()
            .zip(planes.iter().cloned())
            .map(|(position, fault)| Located::new(position, fault))
            .collect();

        let cost = py.detach(|| field_cost_kernel(&located, &grid, kernel, search, min_support));

        Ok((cost.nodes, cost.nodes_to_invert, cost.forward_solutions))
    }

    let (nodes, to_invert, solutions) = at_dimension!(
        dimension,
        run(
            py,
            &view,
            &planes,
            &origin,
            &spacing,
            &counts,
            &bandwidth,
            kernel,
            truncation,
            min_support,
            search,
        )
    )?;

    let dict = PyDict::new(py);
    dict.set_item("nodes", nodes)?;
    dict.set_item("nodes_to_invert", to_invert)?;
    dict.set_item("forward_solutions", solutions)?;

    Ok(dict)
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(covering_grid, m)?)?;
    m.add_function(wrap_pyfunction!(kernel_profile, m)?)?;
    m.add_function(wrap_pyfunction!(density_field, m)?)?;
    m.add_function(wrap_pyfunction!(stress_field, m)?)?;
    m.add_function(wrap_pyfunction!(field_cost, m)?)?;
    Ok(())
}
