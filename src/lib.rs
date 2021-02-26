
use pyo3::prelude::*;
use pyo3::{wrap_pyfunction, wrap_pymodule};
use pyo3::types::IntoPyDict;

use crate::geometry2d::Point2D;
use crate::geometry2d::Segment2D;

mod geometry2d;


#[pyfunction]
fn distance(x1: f64, y1: f64, x2: f64, y2: f64) -> PyResult<f64> {

    let p1 = Point2D{x: x1, y: y1};
    let p2 = Point2D{x: x2, y: y2};

    let out = p1.distance(&p2);
    Ok(out)

}

fn init_submod_geom2d(module: &PyModule) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(distance, module)?)?;
    Ok(())
}

#[pymodule]
fn misah(py: Python, misah_module: &PyModule) -> PyResult<()> {
    let geom2d = PyModule::new(py, "geometry2d")?;
    init_submod_geom2d(geom2d)?;
    misah_module.add_submodule(geom2d)?;
    Ok(())
}


