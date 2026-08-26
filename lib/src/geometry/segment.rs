
use super::point::Point;

#[derive(Debug, Clone, Copy, PartialEq)]
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

    pub fn midpoint(&self) -> Point<N> {
        self.start_pt.midpoint(&self.end_pt)
    }

    /// The point at parameter `t`: `t = 0` is `start_pt`, `t = 1` is `end_pt`.
    ///
    /// Not clamped to `[0, 1]`, matching `Line::point_at` -- a caller wanting
    /// the extrapolation blocked has to say so explicitly.
    pub fn point_at(&self, t: f64) -> Point<N> {
        self.start_pt + self.start_pt.vector_to(&self.end_pt) * t
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn length_of_a_3_4_5_segment_is_5() {
        let s = Segment::new(Point::<2>::from([0.0, 0.0]), Point::from([3.0, 4.0]));
        assert_eq!(s.length(), 5.0);
    }

    #[test]
    fn midpoint_is_halfway_between_the_endpoints() {
        let s = Segment::new(Point::<2>::from([0.0, 0.0]), Point::from([4.0, 2.0]));
        assert_eq!(s.midpoint(), Point::from([2.0, 1.0]));
    }

    #[test]
    fn point_at_zero_and_one_are_the_endpoints() {
        let s = Segment::new(Point::<2>::from([1.0, 1.0]), Point::from([4.0, 5.0]));
        assert_eq!(s.point_at(0.0), s.start_pt);
        assert_eq!(s.point_at(1.0), s.end_pt);
    }

    #[test]
    fn point_at_half_is_the_midpoint() {
        let s = Segment::new(Point::<2>::from([0.0, 0.0]), Point::from([4.0, 2.0]));
        assert_eq!(s.point_at(0.5), s.midpoint());
    }

    #[test]
    fn point_at_extrapolates_past_the_endpoints() {
        let s = Segment::new(Point::<2>::from([0.0, 0.0]), Point::from([1.0, 0.0]));
        assert_eq!(s.point_at(2.0), Point::from([2.0, 0.0]));
    }
}
