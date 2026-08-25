use pyo3::prelude::*;
use pyo3::wrap_pymodule;

use crate::geometry::geom2d::{Point2D, Segment2D};
use crate::geometry::geom3d::{Point3D, Segment3D};

use crate::georeferenced::georef2d::ras2d::GeoTransform;

use crate::orientations::orien3d::{Axis, GeolPlane};

pub mod geometry;
pub mod georeferenced;
pub mod kernels;
pub mod orientations;

// geometry submodules

#[pymodule]
fn geom2d(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Point2D>()?;
    m.add_class::<Segment2D>()?;
    Ok(())
}

#[pymodule]
fn geom3d(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Point3D>()?;
    m.add_class::<Segment3D>()?;
    Ok(())
}

#[pymodule]
#[pyo3(name = "geometry")]
fn geometry_module(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_wrapped(wrap_pymodule!(geom2d))?;
    m.add_wrapped(wrap_pymodule!(geom3d))?;
    Ok(())
}

// georeferenced submodules

#[pymodule]
fn ras2d(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<GeoTransform>()?;
    Ok(())
}

#[pymodule]
fn georef2d(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_wrapped(wrap_pymodule!(ras2d))?;
    Ok(())
}

#[pymodule]
#[pyo3(name = "georeferenced")]
fn georeferenced_module(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_wrapped(wrap_pymodule!(georef2d))?;
    Ok(())
}

// orientations submodules

#[pymodule]
fn orien3d(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Axis>()?;
    m.add_class::<GeolPlane>()?;
    Ok(())
}

#[pymodule]
#[pyo3(name = "orientations")]
fn orientations_module(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_wrapped(wrap_pymodule!(orien3d))?;
    Ok(())
}

// kernels submodule

// Bound to Python only with `extension-module`, which is what pulls in numpy.
#[cfg(feature = "extension-module")]
#[pymodule]
#[pyo3(name = "kernels")]
fn kernels_module(m: &Bound<'_, PyModule>) -> PyResult<()> {
    crate::kernels::py::register(m)
}

// misah module

#[pymodule]
fn misah(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_wrapped(wrap_pymodule!(geometry_module))?;
    m.add_wrapped(wrap_pymodule!(orientations_module))?;
    m.add_wrapped(wrap_pymodule!(georeferenced_module))?;
    #[cfg(feature = "extension-module")]
    m.add_wrapped(wrap_pymodule!(kernels_module))?;

    Ok(())
}
