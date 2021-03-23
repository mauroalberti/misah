use pyo3::prelude::*;

#[pyclass]
#[derive(Default)]
pub struct Point3D {
    pub x: f64,
    pub y: f64,
    pub z: f64
}

#[pymethods]
impl Point3D {

    #[new]
    fn new() -> Self {
        Point3D::default()
    }

    pub fn delta_x(&self, other: &Self) -> f64 {
        other.x - self.x
    }

    pub fn delta_y(&self, other: &Self) -> f64 {
        other.y - self.y
    }

    pub fn delta_z(&self, other: &Self) -> f64 {
        other.z - self.z
    }

    pub fn dist3d(&self, other: &Self) -> f64 {
        (self.delta_x(other) * self.delta_x(other) + self.delta_y(other) * self.delta_y(other) + self.delta_z(other) * self.delta_z(other)).sqrt()
    }

}