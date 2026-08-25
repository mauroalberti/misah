use pyo3::prelude::*;

#[pyclass(from_py_object)]
#[derive(Clone, Copy)]
pub struct GeoTransform {
    #[pyo3(get, set)]
    pub data: [f64; 6]
}