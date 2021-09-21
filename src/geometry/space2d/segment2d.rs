use pyo3::prelude::*;

use super::{Point2D};
use pyo3::PyResult;

#[pyclass]
pub struct Segment2D {
    pub start_pt: Point2D,
    pub end_pt: Point2D,
}

#[pymethods]
impl Segment2D {

    #[new]
    fn new(start_pt: Point2D, end_pt: Point2D) -> Self {
        Segment2D { start_pt, end_pt }
    }

    pub fn length(&self) -> PyResult<f64> {
        Point2D::distance(&self.start_pt, &self.end_pt)
    }
}