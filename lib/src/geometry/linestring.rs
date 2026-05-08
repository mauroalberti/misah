
use super::point::Point;
use super::segment::Segment;


#[derive(Debug, Clone, PartialEq)]
pub struct Linestring<const N: usize> {
    pub points: Vec<Point<N>>,
}

pub type Linestring2D = Linestring<2>;
pub type Linestring3D = Linestring<3>;

impl<const N: usize> Linestring<N> {

    pub fn from_points(points: Vec<Point<N>>) -> Self {
        Self { points }
    }

    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    pub fn num_points(&self) -> usize {
        self.points.len()
    }

    pub fn num_segments(&self) -> usize {
        self.points.len().saturating_sub(1)
    }

    pub fn segments(&self) -> impl Iterator<Item = (&Point<N>, &Point<N>)> {
        self.points
        .windows(2)
        .map(|w| (&w[0], &w[1]))
    }

    pub fn length(&self) -> f64 {
        self.segments()
            .map(|(p0, p1)| p0.distance(p1))
            .sum()
    }

    pub fn segment_lengths(&self) -> Vec<f64> {
        self.points
            .windows(2)
            .map(|pair| pair[0].distance(&pair[1]))
            .collect()
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

    pub fn is_closed(&self) -> bool {
        match (self.first(), self.last()) {
            (Some(first), Some(last)) => first == last,
            _ => false,
        }
    }
}

impl<const N: usize> From<Vec<Point<N>>> for Line<N> {
    fn from(points: Vec<Point<N>>) -> Self {
        Self::new(points)
    }
}
