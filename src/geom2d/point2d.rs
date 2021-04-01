use pyo3::prelude::*;

#[pyclass]
#[derive(Clone, Copy)]
pub struct Point2D {
    pub x: f64,
    pub y: f64,
}

#[pymethods]
impl Point2D {

    #[new]
    fn new(x: f64, y: f64) -> Self {
        Point2D { x, y}
    }

    pub fn delta_x(&self, other: &Self) -> f64 {
        other.x - self.x
    }

    pub fn delta_y(&self, other: &Self) -> f64 {
        other.y - self.y
    }

    pub fn distance(&self, other: &Self) -> PyResult<f64> {
        Ok((self.delta_x(other) * self.delta_x(other) + self.delta_y(other) * self.delta_y(other)).sqrt())
    }

}
