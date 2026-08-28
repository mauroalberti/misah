//! Python bindings for the numerical kernels.

use numpy::{IntoPyArray, PyArray2, PyArrayMethods, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use misah::geometry::point::Point3D;
use misah::raster::geotransform::GeoTransform;
use misah::raster::grid::Grid;
use misah::raster::intersection::intersect_plane_grid as kernel;
use misah::structural::geol_plane::GeologicalPlane;

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
    epsg_code = 0,
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
    epsg_code: i32,
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
        epsg_code,
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

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(intersect_plane_grid, m)?)?;
    m.add_function(wrap_pyfunction!(plane_normal, m)?)?;
    Ok(())
}
