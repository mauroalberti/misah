//! Python bindings for the numerical kernels.

use numpy::{IntoPyArray, PyArray2, PyArrayMethods, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use misah::geometry::point::Point3D;
use misah::raster::geotransform::GeoTransform;
use misah::raster::grid::Grid;
use misah::raster::intersection::intersect_plane_grid as kernel;
use misah::structural::geol_axis::GeologicalAxis;
use misah::structural::geol_plane::GeologicalPlane;
use misah::structural::stress::ReducedStressTensor;

/// What a grid intersection hands back to Python: the vertices, as an (N, 3)
/// array of coordinates, and the chords, as an (M, 2) array of indices into
/// them.
///
/// Named rather than spelled out at the signature, where four nested generics
/// take longer to read than they say.
type VerticesAndSegments<'py> = (Bound<'py, PyArray2<f64>>, Bound<'py, PyArray2<i64>>);

/// Intersect an unbounded geological plane with a DEM.
///
/// `geotransform` is the GDAL six-element affine transform. Returns
/// `(points, segments)`: an (N, 3) array of intersection vertices and an (M, 2)
/// array of indices into it, one row per marching-squares chord.
#[pyfunction]
#[pyo3(signature = (
    dem,
    geotransform,
    src_pt,
    dip_dir_degr,
    dip_angle_degr,
    nodata = None,
))]
#[allow(clippy::too_many_arguments)]
fn intersect_plane_grid<'py>(
    py: Python<'py>,
    dem: PyReadonlyArray2<'py, f64>,
    geotransform: [f64; 6],
    src_pt: [f64; 3],
    dip_dir_degr: f64,
    dip_angle_degr: f64,
    nodata: Option<f64>,
) -> PyResult<VerticesAndSegments<'py>> {
    let view = dem.as_array();
    if !view.is_standard_layout() {
        return Err(PyValueError::new_err(
            "DEM must be C-contiguous; pass numpy.ascontiguousarray(dem)",
        ));
    }

    let plane = GeologicalPlane::new(dip_dir_degr, dip_angle_degr)
        .to_plane(Point3D::from(src_pt))
        .ok_or_else(|| PyValueError::new_err("degenerate plane normal"))?;

    // Grid owns its elevations, so the array is copied here. For a DEM that is a
    // few megabytes at most, and it keeps the kernel free of lifetimes.
    let grid = Grid {
        transform: GeoTransform { data: geotransform },
        data: view.to_owned(),
    };

    let out = py.detach(|| kernel(&plane, &grid, nodata));

    let n_pts = out.points.len();
    let flat_pts: Vec<f64> = out.points.iter().flat_map(|p| p.coords).collect();
    let points = flat_pts.into_pyarray(py).reshape([n_pts, 3])?;

    let n_segs = out.segments.len();
    let flat_segs: Vec<i64> = out
        .segments
        .iter()
        .flat_map(|s| [s[0] as i64, s[1] as i64])
        .collect();
    let segments = flat_segs.into_pyarray(py).reshape([n_segs, 2])?;

    Ok((points, segments))
}

/// Upward-pointing unit normal of a geological plane, in (East, North, Up).
#[pyfunction]
fn plane_normal(dip_dir_degr: f64, dip_angle_degr: f64) -> [f64; 3] {
    GeologicalPlane::new(dip_dir_degr, dip_angle_degr)
        .normal_vector()
        .coords
}

/// The trend/plunge of a slickenline of the given rake, on a plane of the
/// given strike (right-hand rule) and dip.
///
/// Aki & Richards (1980)'s convention: rake 0 is left-lateral, 90 reverse,
/// +/-180 right-lateral, -90 normal.
#[pyfunction]
fn rake_to_slickenline(strike_rhr_degr: f64, dip_angle_degr: f64, rake_degr: f64) -> (f64, f64) {
    let plane = GeologicalPlane::from_rhr_strike(strike_rhr_degr, dip_angle_degr);
    let axis = GeologicalAxis::from_versor(&plane.rake_to_versor(rake_degr));
    (axis.trend, axis.plunge)
}

/// The direct (forward) Wallace-Bott problem: resolve a reduced stress
/// tensor onto one fault plane, predicting the slip it drives.
///
/// `s1`/`s3` are the principal stress axes as (trend, plunge) in degrees,
/// sub-orthogonal to within a degree; `phi` is the shape ratio
/// `(sigma2 - sigma3) / (sigma1 - sigma3)`. `sigma1`/`sigma3` default to 1/0,
/// the usual normalization when only the tensor's shape is known, as from a
/// fault-slip inversion -- which leaves the predicted rake correct while
/// making `slip_tendency`/`deformation_index` meaningless; pass the true
/// magnitudes when they are known and those are wanted.
///
/// Returns a dict always carrying `is_valid`, `traction`, `traction_magnitude`,
/// `normal_stress`, `normal_stress_magnitude`, `shear_stress` and
/// `shear_stress_magnitude` (the last three vectors in (East, North, Up)),
/// plus `theoretical_rake`, `theoretical_slickenline` (a `(trend, plunge)`
/// pair), `slip_tendency` and `deformation_index`, `None` on all four where
/// the shear stress does not clear `shear_threshold` -- the plane sits on, or
/// acutely close to, a principal stress axis, and no slip direction is
/// defined.
#[pyfunction]
#[pyo3(signature = (
    s1_trend_degr,
    s1_plunge_degr,
    s3_trend_degr,
    s3_plunge_degr,
    phi,
    strike_rhr_degr,
    dip_angle_degr,
    sigma1 = 1.0,
    sigma3 = 0.0,
    shear_threshold = misah::structural::stress::SHEAR_MAGNITUDE_THRESHOLD,
))]
#[allow(clippy::too_many_arguments)]
fn solve_stress<'py>(
    py: Python<'py>,
    s1_trend_degr: f64,
    s1_plunge_degr: f64,
    s3_trend_degr: f64,
    s3_plunge_degr: f64,
    phi: f64,
    strike_rhr_degr: f64,
    dip_angle_degr: f64,
    sigma1: f64,
    sigma3: f64,
    shear_threshold: f64,
) -> PyResult<Bound<'py, PyDict>> {

    let s1 = GeologicalAxis::new(s1_trend_degr, s1_plunge_degr);
    let s3 = GeologicalAxis::new(s3_trend_degr, s3_plunge_degr);

    let tensor = ReducedStressTensor::new(s1, s3, phi, sigma1, sigma3)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;

    let plane = GeologicalPlane::from_rhr_strike(strike_rhr_degr, dip_angle_degr);
    let solution = tensor.solve_with_threshold(&plane, shear_threshold);

    let slickenline = solution.theoretical_slickenline.map(|v| {
        let axis = GeologicalAxis::from_versor(&v);
        (axis.trend, axis.plunge)
    });

    let dict = PyDict::new(py);
    dict.set_item("is_valid", solution.is_valid)?;
    dict.set_item("traction", (solution.traction[0], solution.traction[1], solution.traction[2]))?;
    dict.set_item("traction_magnitude", solution.traction_magnitude)?;
    dict.set_item("normal_stress", (solution.normal_stress[0], solution.normal_stress[1], solution.normal_stress[2]))?;
    dict.set_item("normal_stress_magnitude", solution.normal_stress_magnitude)?;
    dict.set_item("shear_stress", (solution.shear_stress[0], solution.shear_stress[1], solution.shear_stress[2]))?;
    dict.set_item("shear_stress_magnitude", solution.shear_stress_magnitude)?;
    dict.set_item("theoretical_rake", solution.theoretical_rake)?;
    dict.set_item("theoretical_slickenline", slickenline)?;
    dict.set_item("slip_tendency", solution.slip_tendency)?;
    dict.set_item("deformation_index", solution.deformation_index)?;

    Ok(dict)
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(intersect_plane_grid, m)?)?;
    m.add_function(wrap_pyfunction!(plane_normal, m)?)?;
    m.add_function(wrap_pyfunction!(rake_to_slickenline, m)?)?;
    m.add_function(wrap_pyfunction!(solve_stress, m)?)?;
    Ok(())
}
