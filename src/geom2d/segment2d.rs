
use super::{Point2D};

pub struct Segment2D {
    pub start_pt: Point2D,
    pub end_pt: Point2D,
}

impl Segment2D {

    pub fn length(&self) -> f64 {
        Point2D::dist2d(&self.start_pt, &self.end_pt)
    }
}