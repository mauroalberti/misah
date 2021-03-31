
use pyo3::prelude::*;
use pyo3::{Py, PyResult, Python};

use crate::geom2d::Point2D;
use crate::geom2d::Segment2D;

use crate::geom3d::Point3D;

mod geom2d;
mod geom3d;


#[pymodule]
fn misah(py: Python, m: &PyModule) -> PyResult<()> {

    let geom2d = PyModule::new(py, "geom2d")?;

    m.add_submodule(geom2d)?;

    let geom3d = PyModule::new(py, "geom3d")?;

    m.add_submodule(geom3d)?;

    m.add_class::<Point3D>()?;

    Ok(())
}


