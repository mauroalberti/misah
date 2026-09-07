//! Python bindings for focal mechanisms and the rotations between them.
//!
//! Kept apart from `kernels.rs` and `fields.rs` for length. From Python these
//! register into the same `misah.kernels` module as everything else.
//!
//! ## How a mechanism crosses the boundary
//!
//! As four numbers: **P trend, P plunge, T trend, T plunge**, in degrees. An
//! `(N, 4)` array of them is a catalogue. B is not carried, being `T x P` and
//! therefore not an independent measurement -- a fourth column pair would be a
//! thing that can disagree with the other two.
//!
//! Nodal planes are the other way a mechanism is usually quoted, and are
//! deliberately not what this takes. Which of the two planes slipped is not
//! something the seismology says, so a pair of planes is an ambiguous input
//! where P and T are not. `ptb_axes` converts from a fault plane and its slip,
//! which is the unambiguous form of the same thing.
//!
//! ## Why the matrix function exists
//!
//! `kagan_angle_matrix` is the reason any of this is in Rust. Comparing a
//! catalogue with itself is `N^2 / 2` rotations, each a handful of quaternion
//! products -- trivial once and ruinous in a Python loop at N in the
//! thousands. It runs across every core and releases the interpreter.

use numpy::{IntoPyArray, PyArray1, PyArray2, PyArrayMethods, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use rayon::iter::{IntoParallelIterator, ParallelIterator};

use misah::structural::focal_mechanism::PTBAxes;
use misah::structural::geol_axis::GeologicalAxis;
use misah::structural::rotation::{
    RotationAxis, focal_mechanism_rotations as rotations_kernel, kagan_angle,
    rotate_focal_mechanism as rotate_kernel,
};

use crate::kernels::faults_from_rows;

/// The mechanisms an (N, 4) array of P and T axes describes.
///
/// Row by row rather than through an iterator, so that a pair of axes that is
/// not a right angle apart names the row it came from. In a catalogue of a
/// thousand that is the difference between a fixable complaint and a search.
fn mechanisms_from_rows(array: &PyReadonlyArray2<'_, f64>) -> PyResult<Vec<PTBAxes>> {

    let view = array.as_array();

    if view.ncols() != 4 {
        return Err(PyValueError::new_err(
            "mechanisms must be an (N, 4) array of P trend, P plunge, T trend, T plunge",
        ));
    }

    let mut out = Vec::with_capacity(view.nrows());

    for (index, row) in view.rows().into_iter().enumerate() {

        let axes = PTBAxes::new(
            GeologicalAxis::new(row[0], row[1]),
            GeologicalAxis::new(row[2], row[3]),
        )
        .map_err(|e| PyValueError::new_err(format!("mechanism {index}: {e}")))?;

        out.push(axes);
    }

    Ok(out)
}

/// The P, T and B axes as three flat columns, in the order the mechanisms came.
fn axes_dict<'py>(py: Python<'py>, mechanisms: &[PTBAxes]) -> PyResult<Bound<'py, PyDict>> {

    let count = mechanisms.len();

    let mut p = Vec::with_capacity(count * 2);
    let mut t = Vec::with_capacity(count * 2);
    let mut b = Vec::with_capacity(count * 2);

    for mechanism in mechanisms {
        let axis = mechanism.p_axis();
        p.extend_from_slice(&[axis.trend, axis.plunge]);
        let axis = mechanism.t_axis();
        t.extend_from_slice(&[axis.trend, axis.plunge]);
        let axis = mechanism.b_axis();
        b.extend_from_slice(&[axis.trend, axis.plunge]);
    }

    let dict = PyDict::new(py);
    dict.set_item("p", p.into_pyarray(py).reshape([count, 2])?)?;
    dict.set_item("t", t.into_pyarray(py).reshape([count, 2])?)?;
    dict.set_item("b", b.into_pyarray(py).reshape([count, 2])?)?;

    Ok(dict)
}

/// The P, T and B kinematic axes of a set of faults.
///
/// `faults` is the same (N, 4) array `invert_stress` takes -- right-hand-rule
/// strike, dip angle, then the trend and plunge of the slickenline -- and
/// `senses` the same optional (N,) boolean array.
///
/// Here, though, **the sense must be known for every fault**. Without it the
/// slip is a line rather than a vector, and reversing it exchanges P with T:
/// the same fault comes back as shortening where it was extension. A fault
/// whose sense nobody could read is refused by name rather than given a
/// direction. That is a stricter rule than `invert_stress`'s, which handles the
/// same missing information by taking its misfit modulo 180 degrees -- an
/// option that does not exist here.
///
/// P and T are **kinematic** axes, at 45 degrees to the fault plane by
/// construction. They are not principal stress axes, and coincide with them
/// only under Anderson's assumption. Where you want stress, `invert_stress` is
/// the function that solves for it.
///
/// Returns a dict of `p`, `t` and `b`, each an (N, 2) array of trend and
/// plunge in degrees.
#[pyfunction]
#[pyo3(signature = (faults, senses = None))]
fn ptb_axes<'py>(
    py: Python<'py>,
    faults: PyReadonlyArray2<'py, f64>,
    senses: Option<PyReadonlyArray1<'py, bool>>,
) -> PyResult<Bound<'py, PyDict>> {

    let planes = faults_from_rows(&faults, senses.as_ref())?;

    let mut mechanisms = Vec::with_capacity(planes.len());

    for (index, fault) in planes.iter().enumerate() {
        mechanisms.push(
            PTBAxes::from_fault_slickenline(fault, 0)
                .map_err(|e| PyValueError::new_err(format!("fault {index}: {e}")))?,
        );
    }

    axes_dict(py, &mechanisms)
}

/// The Kagan angle between mechanisms, pair by pair.
///
/// Both arguments are (N, 4) arrays of P trend, P plunge, T trend, T plunge,
/// of the same length, and the answer is the (N,) array of angles between
/// corresponding rows.
///
/// The Kagan angle is the **smallest** of the four rotations that carry one
/// mechanism onto the other, and is the standard measure of how far apart two
/// focal mechanisms are. It runs from 0 to 120 degrees -- 120 and not 180,
/// because a double couple is unchanged by a half turn about any of its own
/// axes, so no two of them can be further apart than that.
#[pyfunction]
fn kagan_angles<'py>(
    py: Python<'py>,
    first: PyReadonlyArray2<'py, f64>,
    second: PyReadonlyArray2<'py, f64>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {

    let first = mechanisms_from_rows(&first)?;
    let second = mechanisms_from_rows(&second)?;

    if first.len() != second.len() {
        return Err(PyValueError::new_err(format!(
            "{} mechanisms against {}: the two arrays must be the same length",
            first.len(),
            second.len()
        )));
    }

    let angles: Result<Vec<f64>, _> = py.detach(|| {
        first
            .iter()
            .zip(second.iter())
            .map(|(a, b)| kagan_angle(a, b))
            .collect()
    });

    let angles = angles.map_err(|e| PyValueError::new_err(e.to_string()))?;

    Ok(angles.into_pyarray(py))
}

/// Every Kagan angle in a catalogue, as an (N, N) matrix.
///
/// `mechanisms` is an (N, 4) array of P and T axes. The result is symmetric
/// with a zero diagonal, which is what makes it a distance matrix: it can go
/// straight into a clustering routine, and clustering focal mechanisms by
/// Kagan angle is the usual reason to want one.
///
/// Only the upper triangle is computed and it is mirrored, so the cost is
/// `N^2 / 2` rotations rather than `N^2`. That halving matters at the sizes
/// this is for: a catalogue of two thousand mechanisms is two million pairs.
///
/// Runs across every core with the interpreter released.
#[pyfunction]
fn kagan_angle_matrix<'py>(
    py: Python<'py>,
    mechanisms: PyReadonlyArray2<'py, f64>,
) -> PyResult<Bound<'py, PyArray2<f64>>> {

    let mechanisms = mechanisms_from_rows(&mechanisms)?;
    let count = mechanisms.len();

    // One row at a time in parallel, each row computing only the entries to
    // the right of the diagonal and reading the rest from the symmetry.
    let rows: Result<Vec<Vec<f64>>, _> = py.detach(|| {
        (0..count)
            .into_par_iter()
            .map(|i| {
                let mut row = vec![0.0; count];
                for j in (i + 1)..count {
                    row[j] = kagan_angle(&mechanisms[i], &mechanisms[j])?;
                }
                Ok(row)
            })
            .collect()
    });

    let mut rows = rows.map_err(|e: misah::structural::focal_mechanism::FocalMechanismError| {
        PyValueError::new_err(e.to_string())
    })?;

    // Mirror the upper triangle down. Done here rather than inside the
    // parallel pass, where a row cannot see the rows above it.
    //
    // Both indices are load-bearing -- this is a transpose, and it reads
    // `[j][i]` while writing `[i][j]` -- so clippy's advice to iterate one of
    // them away does not apply.
    #[allow(clippy::needless_range_loop)]
    for i in 0..count {
        for j in 0..i {
            rows[i][j] = rows[j][i];
        }
    }

    let flat: Vec<f64> = rows.into_iter().flatten().collect();

    flat.into_pyarray(py).reshape([count, count])
}

/// All four rotations that carry one focal mechanism onto another.
///
/// One pair, given as `(p_trend, p_plunge, t_trend, t_plunge)` twice. Returns a
/// dict of four-element arrays -- `trend`, `plunge` and `angle_degrees` --
/// sorted by the size of the turn, so element 0 is the Kagan angle.
///
/// There are four and not one because a double couple is unchanged by a half
/// turn about any of its own P, T or B axes, so each of those symmetries
/// composed with the base rotation is an equally correct answer. The three
/// beyond the minimum are returned rather than dropped: a minimum quoted
/// without its alternatives is a number nobody can check, and their spread says
/// how well determined the smallest one is.
///
/// The angle is signed and lies in `[-180, 180]` -- the shorter way round about
/// the axis as reported. The axis is directed and its plunge may be negative;
/// the opposite end with the angle negated is the same rotation.
#[pyfunction]
#[pyo3(signature = (
    first_p_trend, first_p_plunge, first_t_trend, first_t_plunge,
    second_p_trend, second_p_plunge, second_t_trend, second_t_plunge,
))]
#[allow(clippy::too_many_arguments)]
fn focal_mechanism_rotations<'py>(
    py: Python<'py>,
    first_p_trend: f64,
    first_p_plunge: f64,
    first_t_trend: f64,
    first_t_plunge: f64,
    second_p_trend: f64,
    second_p_plunge: f64,
    second_t_trend: f64,
    second_t_plunge: f64,
) -> PyResult<Bound<'py, PyDict>> {

    let build = |p_trend, p_plunge, t_trend, t_plunge, which: &str| {
        PTBAxes::new(
            GeologicalAxis::new(p_trend, p_plunge),
            GeologicalAxis::new(t_trend, t_plunge),
        )
        .map_err(|e| PyValueError::new_err(format!("{which} mechanism: {e}")))
    };

    let first = build(first_p_trend, first_p_plunge, first_t_trend, first_t_plunge, "first")?;
    let second = build(
        second_p_trend,
        second_p_plunge,
        second_t_trend,
        second_t_plunge,
        "second",
    )?;

    let rotations =
        rotations_kernel(&first, &second).map_err(|e| PyValueError::new_err(e.to_string()))?;

    let dict = PyDict::new(py);
    dict.set_item(
        "trend",
        rotations.iter().map(|r| r.axis.trend).collect::<Vec<f64>>().into_pyarray(py),
    )?;
    dict.set_item(
        "plunge",
        rotations.iter().map(|r| r.axis.plunge).collect::<Vec<f64>>().into_pyarray(py),
    )?;
    dict.set_item(
        "angle_degrees",
        rotations.iter().map(|r| r.angle_degrees).collect::<Vec<f64>>().into_pyarray(py),
    )?;

    Ok(dict)
}

/// Turn a focal mechanism about an axis.
///
/// The operation `focal_mechanism_rotations` inverts: apply any of the four
/// rotations it returns to the first mechanism and the second comes back.
///
/// Returns a dict of `p`, `t` and `b`, each a `(trend, plunge)` pair, matching
/// the shape `ptb_axes` returns for a single mechanism.
#[pyfunction]
#[pyo3(signature = (p_trend, p_plunge, t_trend, t_plunge, axis_trend, axis_plunge, angle_degrees))]
#[allow(clippy::too_many_arguments)]
fn rotate_focal_mechanism<'py>(
    py: Python<'py>,
    p_trend: f64,
    p_plunge: f64,
    t_trend: f64,
    t_plunge: f64,
    axis_trend: f64,
    axis_plunge: f64,
    angle_degrees: f64,
) -> PyResult<Bound<'py, PyDict>> {

    let mechanism = PTBAxes::new(
        GeologicalAxis::new(p_trend, p_plunge),
        GeologicalAxis::new(t_trend, t_plunge),
    )
    .map_err(|e| PyValueError::new_err(e.to_string()))?;

    let rotation = RotationAxis::new(axis_trend, axis_plunge, angle_degrees);

    let turned = rotate_kernel(&mechanism, &rotation)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;

    let dict = PyDict::new(py);
    let axis = turned.p_axis();
    dict.set_item("p", (axis.trend, axis.plunge))?;
    let axis = turned.t_axis();
    dict.set_item("t", (axis.trend, axis.plunge))?;
    let axis = turned.b_axis();
    dict.set_item("b", (axis.trend, axis.plunge))?;

    Ok(dict)
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(ptb_axes, m)?)?;
    m.add_function(wrap_pyfunction!(kagan_angles, m)?)?;
    m.add_function(wrap_pyfunction!(kagan_angle_matrix, m)?)?;
    m.add_function(wrap_pyfunction!(focal_mechanism_rotations, m)?)?;
    m.add_function(wrap_pyfunction!(rotate_focal_mechanism, m)?)?;
    Ok(())
}
