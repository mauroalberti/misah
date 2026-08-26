
use crate::algebra::EPSILON;

use super::point::Point;


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
        self.segments()
            .map(|(p0, p1)| p0.distance(p1))
            .sum()
    }

    pub fn num_segments(&self) -> usize {
        self.points.len().saturating_sub(1)
    }

    pub fn segments(&self) -> impl Iterator<Item = (&Point<N>, &Point<N>)> {
        self.points
        .windows(2)
        .map(|w| (&w[0], &w[1]))
    }

    pub fn segment_lengths(&self) -> Vec<f64> {
        self.points
            .windows(2)
            .map(|pair| pair[0].distance(&pair[1]))
            .collect()
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
}
