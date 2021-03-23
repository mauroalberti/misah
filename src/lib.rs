
use pyo3::prelude::*;
use pyo3::{wrap_pyfunction, wrap_pymodule};
use pyo3::types::IntoPyDict;

use crate::geom2d::Point2D;
use crate::geom2d::Segment2D;

use crate::geom3d::Point3D;

mod geom2d;
mod geom3d;


#[pyfunction]
fn dist2d(x1: f64, y1: f64, x2: f64, y2: f64) -> PyResult<f64> {

    let p1 = Point2D{x: x1, y: y1};
    let p2 = Point2D{x: x2, y: y2};

    let out = p1.dist2d(&p2);
    Ok(out)

}

fn init_submod_geom2d(module: &PyModule) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(dist2d, module)?)?;
    Ok(())
}

#[pyfunction]
fn dist3d(x1: f64, y1: f64, z1: f64, x2: f64, y2: f64, z2: f64) -> PyResult<f64> {

    let p1 = Point3D{x: x1, y: y1, z: z1};
    let p2 = Point3D{x: x2, y: y2, z: z2};

    let out = p1.dist3d(&p2);
    Ok(out)

}

fn init_submod_geom3d(module: &PyModule) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(dist3d, module)?)?;
    Ok(())
}

#[pymodule]
fn misah(py: Python, m: &PyModule) -> PyResult<()> {

    let geom2d = PyModule::new(py, "geom2d")?;
    init_submod_geom2d(geom2d)?;
    m.add_submodule(geom2d)?;

    let geom3d = PyModule::new(py, "geom3d")?;
    init_submod_geom3d(geom3d)?;
    m.add_submodule(geom3d)?;

    m.add_class::<Point3D>()?;

    Ok(())
}


