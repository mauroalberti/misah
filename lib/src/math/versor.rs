
use std::ops::Neg;

use super::{Vector, VectorError};

// Unit vector
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Versor<const N: usize> {
    coords: [f64, N],
}

impl<const N: usize> Versor<N> {

    pub fn new(coords: [f64; N]) -> Result<Self, VectorError> {
        Vector::new(coords).try_into()
    }

    pub fn coords(&self) -> &[f64; N] {
        &self.coords
    }

    pub fn as_vector(&self) -> Vector<N> {
        Vector::new(self.coords)
    }

    pub fn opposite(&self) -> Self {
        Self {
            coords: self.coords.map(|c| -c),
        }
    }

    pub fn dot(&self, other: &Self) -> f64 {
        self.as_vector().dot(&other.as_vector())
    }
}

impl<const N: usize> TryFrom<Vector<N>> for Versor<N> {

    type Error = VectorError;

    fn try_from(v: Vector<N>) -> Result<Self, Self::Error> {

        let norm = v.norm();
        if norm == 0.0 {
            return Err(VectorError::ZeroVector);
        }

        Ok(Self {
            coords: v.coords.map(|c| c / norm),
        })

    }
}

impl<const N: usize> Neg for Versor<N> {

    type Output = Self;

    fn neg(self) -> Self::Output {
        self.opposite()
    }
}
