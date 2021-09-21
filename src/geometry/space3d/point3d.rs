use pyo3::prelude::*;

#[pyclass]
#[derive(Clone, Copy)]
pub struct Point3D {
    #[pyo3(get, set)]
    pub x: f64,
    #[pyo3(get, set)]
    pub y: f64,
    #[pyo3(get, set)]
    pub z: f64
}

#[pymethods]
impl Point3D {

    #[new]
    fn new(x: f64, y: f64, z: f64) -> Self {
        Point3D { x, y, z}
    }

    pub fn delta_x(&self, other: &Self) -> PyResult<f64> {
        Ok(other.x - self.x)
    }

    pub fn delta_y(&self, other: &Self) -> PyResult<f64> {
        Ok(other.y - self.y)
    }

    pub fn delta_z(&self, other: &Self) -> PyResult<f64> {
        Ok(other.z - self.z)
    }

    pub fn distance(&self, other: &Self) -> PyResult<f64> {
        Ok((self.delta_x(other)? * self.delta_x(other)? + self.delta_y(other)? * self.delta_y(other)? + self.delta_z(other)? * self.delta_z(other)?).sqrt())
    }


}