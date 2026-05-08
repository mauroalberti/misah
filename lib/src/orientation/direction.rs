
use crate::math::{MathError, Vector, Versor};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Direction<const N: usize> {
    versor: Versor<N>,
}

impl<const N: usize> Direction<N> {

    pub fn new(versor: Versor<N>) -> Self {
        Self { versor }
    }

    pub fn from_vector(vector: Vector<N>) -> Result<Self, VectorError> {
        Ok(Self {
            versor: vector.try_into()?,
        })
    }

    pub fn versor(&self) -> Versor<N> {
        self.versor
    }

    pub fn as_vector(&self) -> Vector<N> {
        self.versor.as_vector()
    }

    pub fn opposite(&self) -> Self {
        Self {
            versor: self.versor.opposite()
        }
    }

    pub fn dot(&self, other: &Self) -> f64 {
        self.versor.dot(&other.versor)
    }

    pub fn angle_rad(&self, other: &Self) -> f64 {
        self.dor(other).clamp(-1.0, 1.0).acos()
    }

    pub fn angle_deg(&self, other: &Self) -> f64 {
        self.angle_rad(other).to_degrees()
    }
}

pub type Direction2D = Direction<2>;
pub type Direction3D = Direction<3>;


