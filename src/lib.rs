
use pyo3::prelude::*;
use pyo3::{PyResult, Python, py_run, wrap_pyfunction, wrap_pymodule};

use crate::geometry::geom2d::{Point2D, Segment2D};
use crate::geometry::geom3d::{Point3D, Segment3D};

use crate::georeferenced::georef2d::ras2d::{GeoArray, GeoTransform};

use crate::orientations::orien3d::{Axis, GeolPlane};

pub mod geometry;
pub mod georeferenced;
pub mod kernels;
pub mod orientations;


// geometry submodules

#[pymodule]
fn geom2d(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<Point2D>()?;
    m.add_class::<Segment2D>()?;
    Ok(())
}

#[pymodule]
fn geom3d(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<Point3D>()?;
    m.add_class::<Segment3D>()?;
    Ok(())
}

#[pymodule]
fn geometry(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pymodule!(geom2d))?;
    m.add_wrapped(wrap_pymodule!(geom3d))?;
    Ok(())
}

// georeferenced submodules

#[pymodule]
fn ras2d(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<GeoTransform>()?;
    Ok(())
}

#[pymodule]
fn georef2d(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pymodule!(ras2d))?;
    Ok(())
}

#[pymodule]
fn georeferenced(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pymodule!(georef2d))?;
    Ok(())
}

// orientations submodules

#[pymodule]
fn orien3d(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<Axis>()?;
    m.add_class::<GeolPlane>()?;
    Ok(())
}

#[pymodule]
fn orientations(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pymodule!(orien3d))?;
    Ok(())
}

// misah module

#[pymodule]
fn misah(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pymodule!(geometry))?;
    m.add_wrapped(wrap_pymodule!(orientations))?;
    m.add_wrapped(wrap_pymodule!(georeferenced))?;

    Ok(())
}


