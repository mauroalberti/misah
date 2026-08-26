
use crate::algebra::EPSILON;

use super::point::Point;
use super::segment::Segment;


#[derive(Debug, Clone, PartialEq)]
pub struct Linestring<const N: usize> {
    pub points: Vec<Point<N>>,
}

pub type Linestring2D = Linestring<2>;
pub type Linestring3D = Linestring<3>;

impl<const N: usize> Linestring<N> {

    pub fn new(points: Vec<Point<N>>) -> Self {
        Self { points }
    }

    pub fn from_points(points: Vec<Point<N>>) -> Self {
        Self { points }
    }

    pub fn empty() -> Self {
        Self { points: Vec::new() }
    }

    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    pub fn num_points(&self) -> usize {
        self.points.len()
    }

    pub fn len(&self) -> usize {
        self.points.len()
    }

    pub fn length(&self) -> f64 {
        self.segments().map(|s| s.length()).sum()
    }

    pub fn num_segments(&self) -> usize {
        self.points.len().saturating_sub(1)
    }

    pub fn segments(&self) -> impl Iterator<Item = Segment<N>> + '_ {
        self.points
            .windows(2)
            .map(|w| Segment::new(w[0], w[1]))
    }

    pub fn segment_lengths(&self) -> Vec<f64> {
        self.segments().map(|s| s.length()).collect()
    }

    /// The point at arc-length distance `s` from the first vertex, walking the
    /// segments in order.
    ///
    /// `None` if `s` is negative, if it exceeds the total length by more than
    /// `EPSILON`, or if the linestring has no points. `s` beyond the length by
    /// up to `EPSILON` still resolves to the last point, since summing segment
    /// lengths one at a time does not agree bit-for-bit with reaching the last
    /// vertex by walking them.
    pub fn point_at_distance(&self, s: f64) -> Option<Point<N>> {
        if s < 0.0 {
            return None;
        }

        let first = self.first()?;

        if s == 0.0 {
            return Some(*first);
        }

        let mut remaining = s;
        for seg in self.segments() {
            let seg_len = seg.length();
            if remaining <= seg_len {
                let t = if seg_len > EPSILON { remaining / seg_len } else { 0.0 };
                return Some(seg.point_at(t));
            }
            remaining -= seg_len;
        }

        if remaining <= EPSILON {
            self.last().copied()
        } else {
            None
        }
    }

    pub fn point(&self, i: usize) -> Option<&Point<N>> {
        self.points.get(i)
    }

    pub fn first(&self) -> Option<&Point<N>> {
        self.points.first()
    }

    pub fn last(&self) -> Option<&Point<N>> {
        self.points.last()
    }

    pub fn push(&mut self, point: Point<N>) {
        self.points.push(point);
    }

    /// Whether the first and last vertex coincide, within `EPSILON`.
    ///
    /// A ring read back from storage or through a reprojection rarely has its
    /// closing vertex bit-identical to its first one; comparing with `==` (exact
    /// `Point` equality) called such a ring open.
    pub fn is_closed(&self) -> bool {
        match (self.first(), self.last()) {
            (Some(first), Some(last)) => first.approx_eq(last, EPSILON),
            _ => false,
        }
    }
}

impl<const N: usize> From<Vec<Point<N>>> for Linestring<N> {
    fn from(points: Vec<Point<N>>) -> Self {
        Self::new(points)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_linestring_is_not_closed() {
        assert!(!Linestring::<2>::empty().is_closed());
    }

    #[test]
    fn open_ring_is_not_closed() {
        let ls = Linestring::from(vec![
            Point::<2>::from([0.0, 0.0]),
            Point::from([1.0, 0.0]),
            Point::from([1.0, 1.0]),
        ]);

        assert!(!ls.is_closed());
    }

    #[test]
    fn exactly_matching_endpoints_are_closed() {
        let ls = Linestring::from(vec![
            Point::<2>::from([0.0, 0.0]),
            Point::from([1.0, 0.0]),
            Point::from([1.0, 1.0]),
            Point::from([0.0, 0.0]),
        ]);

        assert!(ls.is_closed());
    }

    #[test]
    fn endpoints_within_epsilon_are_closed() {
        // The closing vertex is a fraction of a nanometre off the first one --
        // exactly the kind of mismatch a reprojection round trip leaves behind.
        // Exact Point equality used to call this ring open.
        let ls = Linestring::from(vec![
            Point::<2>::from([0.0, 0.0]),
            Point::from([1.0, 0.0]),
            Point::from([1.0, 1.0]),
            Point::from([1e-10, 0.0]),
        ]);

        assert!(ls.is_closed());
    }

    #[test]
    fn endpoints_past_epsilon_are_not_closed() {
        let ls = Linestring::from(vec![
            Point::<2>::from([0.0, 0.0]),
            Point::from([1.0, 0.0]),
            Point::from([1e-3, 0.0]),
        ]);

        assert!(!ls.is_closed());
    }

    fn l_shaped() -> Linestring<2> {
        // (0,0) -> (3,0) -> (3,4): legs of length 3 and 4, total length 7.
        Linestring::from(vec![
            Point::from([0.0, 0.0]),
            Point::from([3.0, 0.0]),
            Point::from([3.0, 4.0]),
        ])
    }

    #[test]
    fn segments_yields_one_segment_per_consecutive_pair() {
        let ls = l_shaped();
        let segs: Vec<Segment<2>> = ls.segments().collect();

        assert_eq!(segs.len(), 2);
        assert_eq!(segs[0], Segment::new(Point::from([0.0, 0.0]), Point::from([3.0, 0.0])));
        assert_eq!(segs[1], Segment::new(Point::from([3.0, 0.0]), Point::from([3.0, 4.0])));
    }

    #[test]
    fn length_sums_the_segment_lengths() {
        let ls = l_shaped();
        assert_eq!(ls.length(), 7.0);
        assert_eq!(ls.segment_lengths(), vec![3.0, 4.0]);
    }

    #[test]
    fn point_at_distance_zero_is_the_first_point() {
        let ls = l_shaped();
        assert_eq!(ls.point_at_distance(0.0), Some(Point::from([0.0, 0.0])));
    }

    #[test]
    fn point_at_distance_within_the_first_segment() {
        let ls = l_shaped();
        assert_eq!(ls.point_at_distance(1.5), Some(Point::from([1.5, 0.0])));
    }

    #[test]
    fn point_at_distance_at_a_vertex() {
        let ls = l_shaped();
        assert_eq!(ls.point_at_distance(3.0), Some(Point::from([3.0, 0.0])));
    }

    #[test]
    fn point_at_distance_within_the_second_segment() {
        let ls = l_shaped();
        assert_eq!(ls.point_at_distance(5.0), Some(Point::from([3.0, 2.0])));
    }

    #[test]
    fn point_at_the_total_length_is_the_last_point() {
        let ls = l_shaped();
        assert_eq!(ls.point_at_distance(ls.length()), Some(Point::from([3.0, 4.0])));
    }

    #[test]
    fn point_at_negative_distance_is_none() {
        assert!(l_shaped().point_at_distance(-1.0).is_none());
    }

    #[test]
    fn point_at_distance_past_the_end_is_none() {
        assert!(l_shaped().point_at_distance(100.0).is_none());
    }

    #[test]
    fn point_at_distance_on_an_empty_linestring_is_none() {
        assert!(Linestring::<2>::empty().point_at_distance(0.0).is_none());
    }

    #[test]
    fn point_at_distance_zero_on_a_single_point_linestring() {
        let ls = Linestring::from(vec![Point::<2>::from([1.0, 1.0])]);
        assert_eq!(ls.point_at_distance(0.0), Some(Point::from([1.0, 1.0])));
        assert!(ls.point_at_distance(1.0).is_none());
    }
}
