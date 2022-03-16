
use super::point3d::{Point3D};

pub struct Segment3D {
    pub start_pt: Point3D,
    pub end_pt: Point3D,
}

impl Segment3D {

    fn new(start_pt: Point3D, end_pt: Point3D) -> Self {
        Segment3D { start_pt, end_pt }
    }

    pub fn length(&self) -> f64 {
        Point3D::distance(&self.start_pt, &self.end_pt)
    }
}