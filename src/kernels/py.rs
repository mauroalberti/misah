//! Python bindings for the kernels.
//!
//! Kept apart from the kernels themselves so those stay free of pyo3 and remain
//! testable with `cargo test` alone.

use numpy::{IntoPyArray, PyArray2, PyArrayMethods, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use super::plane_dem;

/// Intersect an unbounded geological plane with a DEM.
///
/// `geotransform` is the GDAL six-element affine transform. Returns
/// `(points, segments)`: an (N, 3) array of intersection vertices and an (M, 2)
/// array of indices into it, one row per marching-squares chord.
#[pyfunction]
#[pyo3(signature = (dem, geotransform, src_pt, dip_dir_degr, dip_angle_degr, nodata=None))]
#[allow(clippy::too_many_arguments)]
fn intersect_plane_dem<'py>(
    py: Python<'py>,
    dem: PyReadonlyArray2<'py, f64>,
    geotransform: [f64; 6],
    src_pt: [f64; 3],
    dip_dir_degr: f64,
    dip_angle_degr: f64,
    nodata: Option<f64>,
) -> PyResult<(Bound<'py, PyArray2<f64>>, Bound<'py, PyArray2<i64>>)> {
    let view = dem.as_array();
    let (nrows, ncols) = view.dim();
    let slice = view.as_slice().ok_or_else(|| {
        PyValueError::new_err("DEM must be C-contiguous; pass numpy.ascontiguousarray(dem)")
    })?;

    let gt = plane_dem::GeoTransform {
        x_origin: geotransform[0],
        pixel_width: geotransform[1],
        row_rotation: geotransform[2],
        y_origin: geotransform[3],
        col_rotation: geotransform[4],
        pixel_height: geotransform[5],
    };

    let out = py.detach(|| {
        plane_dem::intersect_plane_dem(
            slice, nrows, ncols, &gt, src_pt, dip_dir_degr, dip_angle_degr, nodata,
        )
    });

    let n_pts = out.points.len();
    let flat_pts: Vec<f64> = out.points.into_iter().flatten().collect();
    let points = flat_pts
        .into_pyarray(py)
        .reshape([n_pts, 3])?;

    let n_segs = out.segments.len();
    let flat_segs: Vec<i64> = out
        .segments
        .into_iter()
        .flat_map(|s| [s[0] as i64, s[1] as i64])
        .collect();
    let segments = flat_segs
        .into_pyarray(py)
        .reshape([n_segs, 2])?;

    Ok((points, segments))
}

/// Upward-pointing unit normal of a geological plane, in (East, North, Up).
#[pyfunction]
fn plane_normal(dip_dir_degr: f64, dip_angle_degr: f64) -> [f64; 3] {
    plane_dem::plane_normal(dip_dir_degr, dip_angle_degr)
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(intersect_plane_dem, m)?)?;
    m.add_function(wrap_pyfunction!(plane_normal, m)?)?;
    Ok(())
}
