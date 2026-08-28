
use std::convert::TryInto;
use std::convert::TryFrom;
use std::ops::Neg;

use super::constants::EPSILON;
use super::{Vector, AlgebraError};

// Unit vector
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Versor<const N: usize> {
    coords: [f64; N],
}

pub type Versor2D = Versor<2>;
pub type Versor3D = Versor<3>;

impl<const N: usize> Versor<N> {

    pub fn new(coords: [f64; N]) -> Result<Self, AlgebraError> {
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

    type Error = AlgebraError;

    /// `AlgebraError::ZeroVector` when the vector is too small to normalize.
    ///
    /// The threshold is `EPSILON`, the same one `Vector::normalize` and
    /// `Vector::is_zero` use, and for the same reason: before this it was
    /// `norm == 0.0` here, so a vector that `is_zero` called zero and
    /// `normalize` refused was still accepted as a direction through this
    /// path. Dividing by a norm that small yields components dominated by
    /// whatever rounding noise produced the vector, which is a direction only
    /// in type.
    fn try_from(v: Vector<N>) -> Result<Self, Self::Error> {

        let norm = v.norm();

        if norm < EPSILON {
            return Err(AlgebraError::ZeroVector);
        }

        let coords = v.coords.map(|c| c / norm);

        Ok(Self { coords })

    }
}

impl<const N: usize> Neg for Versor<N> {

    type Output = Self;

    fn neg(self) -> Self::Output {
        self.opposite()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_scales_to_unit_length() {
        let v = Versor::<3>::new([0.0, 3.0, 4.0]).expect("a non-zero vector");

        assert!((v.as_vector().norm() - 1.0).abs() < 1e-12);
        assert!((v.coords()[1] - 0.6).abs() < 1e-12);
        assert!((v.coords()[2] - 0.8).abs() < 1e-12);
    }

    #[test]
    fn the_zero_vector_is_not_a_direction() {
        assert_eq!(Versor::<3>::new([0.0, 0.0, 0.0]), Err(AlgebraError::ZeroVector));
    }

    #[test]
    fn a_sub_epsilon_vector_is_not_a_direction_either() {
        // The case that used to tell the three paths apart: is_zero() called
        // this zero and normalize() refused it, while Versor took it and
        // returned a unit vector built from rounding noise.
        let v = Vector::<3>::from([EPSILON / 10.0, 0.0, 0.0]);

        assert!(v.is_zero());
        assert!(v.normalize().is_none());
        assert_eq!(Versor::try_from(v), Err(AlgebraError::ZeroVector));
    }

    #[test]
    fn versor_agrees_with_normalize_wherever_both_answer() {
        for coords in [
            [1.0, 0.0, 0.0],
            [1.0, 2.0, 3.0],
            [-4.0, 0.5, 7.25],
            [EPSILON * 10.0, 0.0, 0.0],
        ] {
            let v = Vector::<3>::from(coords);

            let normalized = v.normalize().expect("above the threshold");
            let versor = Versor::try_from(v).expect("and so accepted here too");

            assert_eq!(&normalized.coords, versor.coords(), "{:?}", coords);
        }
    }

    #[test]
    fn opposite_and_neg_are_the_same_operation() {
        let v = Versor::<3>::new([1.0, 2.0, 2.0]).unwrap();

        assert_eq!(v.opposite(), -v);
        assert_eq!((-v).opposite(), v);
    }

    #[test]
    fn dot_of_a_versor_with_itself_is_one() {
        let v = Versor::<3>::new([1.0, -2.0, 0.5]).unwrap();

        assert!((v.dot(&v) - 1.0).abs() < 1e-12);
    }
}
