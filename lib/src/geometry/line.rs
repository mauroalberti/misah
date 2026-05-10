
use std::convert::TryInto;

use crate::algebra::{Versor, AlgebraError};
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

}
