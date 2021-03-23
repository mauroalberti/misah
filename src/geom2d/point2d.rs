
pub struct Point2D {
    pub x: f64,
    pub y: f64,
}

impl Point2D {

    pub fn dist2d(&self, other: &Self) -> f64 {
        ((self.x - other.x) * (self.x - other.x) + (self.y - other.y) * (self.y - other.y)).sqrt()
    }

}