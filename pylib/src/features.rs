use misah::features;

use std::thread;
use pyo3::prelude::*;
use pyo3::{PyResult, Python, py_run, wrap_pyfunction, wrap_pymodule};

//* https://github.com/PyO3/pyo3/issues/1517#issuecomment-808664021

#[pymodule]
fn features(py: Python, module: &PyModule) -> PyResult<()> {
    let points = PyModule::new(py, "features.points")?;
    py_run!(py, module, "import sys; sys.modules['features.points'] = points");
    module.add_submodule(points)?;

    Ok(())
}





