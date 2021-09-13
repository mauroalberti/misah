use pyo3::prelude::*;

#[pyclass]
#[derive(Clone, Copy)]
pub struct Axis {
    #[pyo3(get)]
    pub trend: f64,
    #[pyo3(get)]
    pub plunge: f64
}

#[pymethods]
impl Axis {

    #[new]
    fn new(tr: f64, pl: f64) -> Self { Axis{ trend: tr, plunge: pl }}

}