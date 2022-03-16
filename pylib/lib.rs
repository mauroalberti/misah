
use crate::features::{Point3D, Segment3D};

use crate::rasters::{GeoArray, GeoTransform};

use crate::orientations::{Axis, GeolPlane};

pub mod features;
pub mod rasters;
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


