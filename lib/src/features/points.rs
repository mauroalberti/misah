
pub struct Point3D {
    pub x: f64,
    pub y: f64,
    pub z: f64
}

impl Point3D {

    fn new(x: f64, y: f64, z: f64) -> Self {
        Point3D { x, y, z}
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

    pub fn distance(&self, other: &Self) -> f64 {
        (self.delta_x(other) * self.delta_x(other) + self.delta_y(other) * self.delta_y(other) + self.delta_z(other) * self.delta_z(other)).sqrt()
    }


}