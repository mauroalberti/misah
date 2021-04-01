
use pyo3::prelude::*;
use pyo3::{PyResult, Python};

use crate::geom2d::Point2D;
use crate::geom2d::Segment2D;

use crate::geom3d::Point3D;
use crate::geom3d::Segment3D;

mod geom2d;
mod geom3d;


#[pymodule]
fn misah(py: Python, m: &PyModule) -> PyResult<()> {

    let geom2d = PyModule::new(py, "geom2d")?;

    m.add_submodule(geom2d)?;

    geom2d.add_class::<Point2D>()?;
    geom2d.add_class::<Segment2D>()?;

    let geom3d = PyModule::new(py, "geom3d")?;

    m.add_submodule(geom3d)?;

    geom3d.add_class::<Point3D>()?;
    geom3d.add_class::<Segment3D>()?;

    Ok(())
}


