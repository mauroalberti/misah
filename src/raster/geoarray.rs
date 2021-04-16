use pyo3::prelude::*;
use ndarray::Array2;

// modified from: https://stackoverflow.com/questions/13212212/creating-two-dimensional-arrays-in-rust

#[pyclass]
#[derive(Clone, Copy)]
pub struct GeoArray {
    #[pyo3(get, set)]
    pub epsg_code: i32,
    #[pyo3(get, set)]
    pub data: Array2<f64>
}