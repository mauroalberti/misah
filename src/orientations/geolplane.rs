use pyo3::prelude::*;

use super::{Axis};
use pyo3::PyResult;

#[pyclass]
#[derive(Clone, Copy)]
pub struct GeolPlane {
    #[pyo3(get)]
    pub azimuth: f64,
    #[pyo3(get)]
    pub dip_angle: f64
}

#[pymethods]
impl GeolPlane {

    #[new]
    fn new(az: f64, dip: f64) -> Self { GeolPlane{ azimuth: az, dip_angle: dip }}

    pub fn normal_axis(&self) -> PyResult<Axis> {
        Ok(Axis{trend: (self.azimuth + 180.0) % 360.0, plunge: 90.0 - self.dip_angle})
    }

}