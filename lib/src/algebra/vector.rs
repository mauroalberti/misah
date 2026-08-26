
use std::ops::{Add, Neg, Sub, Mul, Div};

use super::constants::EPSILON;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vector<const N: usize> {
    pub coords: [f64; N],
}

pub type Vector2D = Vector<2>;
pub type Vector3D = Vector<3>;

impl<const N: usize> Vector<N> {

    pub fn new(coords: [f64; N]) -> Self {
        Self { coords }
    }

    pub fn coord(&self, i: usize) -> Option<f64> {
        self.coords.get(i).copied()
    }

    pub fn norm(&self) -> f64 {
        self.coords
            .iter()
            .map(|c| c*c)
            .sum::<f64>()
            .sqrt()
    }

    pub fn normalize(&self) -> Option<Self> {

        let n = self.norm();

        if n < EPSILON {
            return None;
        }

        let coords = self.coords.map(|v| v / n);
        Some(Self {coords })

    }

    pub fn dot(&self, other: &Self) -> f64 {
        self.coords
            .iter()
            .zip(other.coords.iter())
            .map(|(a, b)| a * b)
            .sum()
    }

    /// Whether this vector is too small to normalize.
    ///
    /// Kept consistent with `normalize`'s own threshold on purpose: before this,
    /// a vector with norm below `EPSILON` but no single component exactly `0.0`
    /// answered `false` here while `normalize` still returned `None` on it, so
    /// `if !v.is_zero() { v.normalize().unwrap() }` could panic on the exact case
    /// this method exists to guard against.
    pub fn is_zero(&self) -> bool {
        self.norm() < EPSILON
    }

}

impl<const N: usize> From<[f64; N]> for Vector<N> {
    fn from(coords: [f64; N]) -> Self {
        Self::new(coords)
    }
}

impl<const N: usize> Neg for Vector<N> {

    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            coords: self.coords.map(|c| -c ),
        }
    }
}

impl<const N: usize> Mul<f64> for Vector<N> {

    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            coords: self.coords.map(|c| c * rhs ),
        }
    }
}

impl<const N: usize> Div<f64> for Vector<N> {

    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Self {
            coords: self.coords.map(|c| c / rhs ),
        }
    }
}

impl<const N: usize> Add for Vector<N> {

    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        Self {
            coords: std::array::from_fn(|i| self.coords[i] + other.coords[i]),
        }
    }
}

impl<const N: usize> Sub for Vector<N> {

    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {
        Self {
            coords: std::array::from_fn(|i| self.coords[i] - other.coords[i]),
        }
    }
}

impl Vector<3> {

    pub fn cross(&self, other: &Self) -> Self {

        let [a1, a2, a3] = self.coords;
        let [b1, b2, b3] = other.coords;

        Self {
            coords: [
              a2 * b3 - a3 * b2,
              a3 * b1 - a1 * b3,
              a1 * b2 - a2 * b1,
            ],
        }


    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coord_returns_the_requested_component() {
        let v = Vector::<3>::from([1.0, 2.0, 3.0]);
        assert_eq!(v.coord(1), Some(2.0));
    }

    #[test]
    fn coord_out_of_range_is_none() {
        let v = Vector::<3>::from([1.0, 2.0, 3.0]);
        assert_eq!(v.coord(3), None);
    }

    #[test]
    fn norm_of_a_3_4_5_vector_is_5() {
        let v = Vector::<2>::from([3.0, 4.0]);
        assert_eq!(v.norm(), 5.0);
    }

    #[test]
    fn norm_of_the_zero_vector_is_zero() {
        let v = Vector::<3>::from([0.0, 0.0, 0.0]);
        assert_eq!(v.norm(), 0.0);
    }

    #[test]
    fn dot_of_orthogonal_vectors_is_zero() {
        let x = Vector::<3>::from([1.0, 0.0, 0.0]);
        let y = Vector::<3>::from([0.0, 1.0, 0.0]);
        assert_eq!(x.dot(&y), 0.0);
    }

    #[test]
    fn dot_of_parallel_vectors_is_the_product_of_norms() {
        let a = Vector::<2>::from([3.0, 0.0]);
        let b = Vector::<2>::from([5.0, 0.0]);
        assert_eq!(a.dot(&b), 15.0);
    }

    #[test]
    fn cross_of_x_and_y_axes_is_z_axis() {
        let x = Vector::<3>::from([1.0, 0.0, 0.0]);
        let y = Vector::<3>::from([0.0, 1.0, 0.0]);
        assert_eq!(x.cross(&y), Vector::from([0.0, 0.0, 1.0]));
    }

    #[test]
    fn cross_of_a_vector_with_itself_is_zero() {
        let v = Vector::<3>::from([1.0, 2.0, 3.0]);
        assert!(v.cross(&v).is_zero());
    }

    #[test]
    fn normalize_scales_to_unit_length() {
        let v = Vector::<2>::from([3.0, 4.0]);
        let n = v.normalize().expect("non-zero vector");
        assert!((n.norm() - 1.0).abs() < 1e-12);
        assert!((n.coords[0] - 0.6).abs() < 1e-12);
        assert!((n.coords[1] - 0.8).abs() < 1e-12);
    }

    #[test]
    fn normalize_of_the_zero_vector_is_none() {
        let v = Vector::<3>::from([0.0, 0.0, 0.0]);
        assert!(v.normalize().is_none());
    }

    #[test]
    fn normalize_of_a_vector_below_epsilon_is_none() {
        // Below EPSILON but no component is exactly 0.0 -- the case that used to
        // tell is_zero() and normalize() apart.
        let v = Vector::<3>::from([EPSILON / 10.0, 0.0, 0.0]);
        assert!(v.normalize().is_none());
    }

    #[test]
    fn is_zero_agrees_with_normalize_on_a_sub_epsilon_vector() {
        let v = Vector::<3>::from([EPSILON / 10.0, 0.0, 0.0]);
        assert!(v.is_zero());
        assert!(v.normalize().is_none());
    }

    #[test]
    fn is_zero_is_false_for_a_clearly_nonzero_vector() {
        let v = Vector::<3>::from([1.0, 0.0, 0.0]);
        assert!(!v.is_zero());
    }

    #[test]
    fn add_sub_and_neg_are_consistent() {
        let a = Vector::<3>::from([1.0, 2.0, 3.0]);
        let b = Vector::<3>::from([4.0, -1.0, 0.5]);

        assert_eq!(a + b, Vector::from([5.0, 1.0, 3.5]));
        assert_eq!(a - b, Vector::from([-3.0, 3.0, 2.5]));
        assert_eq!(-a, Vector::from([-1.0, -2.0, -3.0]));
        assert_eq!(a - b, -(b - a));
    }

    #[test]
    fn mul_and_div_by_scalar_are_inverse() {
        let v = Vector::<2>::from([2.0, -3.0]);
        assert_eq!(v * 2.0, Vector::from([4.0, -6.0]));
        assert_eq!((v * 2.0) / 2.0, v);
    }
}

