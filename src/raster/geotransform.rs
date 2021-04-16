use pyo3::prelude::*;

#[pyclass]
#[derive(Clone, Copy)]
pub struct GeoTransform {
    #[pyo3(get, set)]
    pub data: [f64; 6]
}