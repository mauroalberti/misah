use pyo3::prelude::*;

use super::{Point3D};
use pyo3::PyResult;

#[pyclass]
pub struct Segment3D {
    pub start_pt: Point3D,
    pub end_pt: Point3D,
}

#[pymethods]
impl Segment3D {

    #[new]
    fn new(start_pt: Point3D, end_pt: Point3D) -> Self {
        Segment3D { start_pt, end_pt }
    }

    pub fn length(&self) -> PyResult<f64> {
        Point3D::distance(&self.start_pt, &self.end_pt)
    }
}