
use super::points::Point;

pub struct Segment<const N: usize> {
    pub start_pt: Point<N>,
    pub end_pt: Point<N>,
}

pub type Segment2D = Segment<2>;
pub type Segment3D = Segment<3>;

impl<const N: usize> Segment<N> {

    pub fn new(start_pt: Point<N>, end_pt: Point<N>) -> Self {
        Self { start_pt, end_pt }
    }

    pub fn length(&self) -> f64 {
        self.start_pt.distance(&self.end_pt)
    }
}
