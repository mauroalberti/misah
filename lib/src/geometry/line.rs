
use std::convert::TryInto;

use crate::algebra::{AlgebraError, Versor};
use crate::geometry::point::Point;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Line<const N: usize> {
    pub origin: Point<N>,
    pub direction: Versor<N>,
}

pub type Line2D = Line<2>;
pub type Line3D = Line<3>;

impl<const N: usize> Line<N> {

    pub fn new<D>(origin: Point<N>, direction: D) -> Result<Self, D::Error>
    where
        D: TryInto<Versor<N>>,
    {

        let direction = direction.try_into()?;

        Ok(Self { origin, direction })
    }

    pub fn from_two_points(p0: Point<N>, p1: Point<N>) -> Result<Self, AlgebraError> {
        let direction = p0.vector_to(&p1);
        Self::new(p0, direction)
    }

    pub fn point_at(&self, t: f64) -> Point<N> {
        self.origin + self.direction.as_vector() * t
    }

}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra::vector::Vector;
    use crate::geometry::point::Point;

    #[test]
    fn line_point_at_zero_is_origin() {
        let origin = Point::<3>::from([1.0, 2.0, 3.0]);
        let direction = Versor::new( [1.0, 0.0, 0.0] ).unwrap();

        let line = Line::new(origin, direction).unwrap();
        let p = line.point_at(0.0);

        assert_eq!(p, origin);
    }

    #[test]
    fn line_point_at_parameter_is_correct() {
        let origin = Point::<3>::from([0.0, 0.0, 0.0]);
        let direction = Versor::<3>::new([10.0, 0.0, 0.0]).unwrap();

        let line = Line::new(origin, direction).unwrap();
        let p = line.point_at(2.0);

        assert_eq!(p, Point::<3>::from([2.0, 0.0, 0.0]));
    }

}
