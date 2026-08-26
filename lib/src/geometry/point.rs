
use std::ops::{Add, Sub};

use crate::algebra::Vector;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point<const N: usize> {
    pub coords: [f64; N],
}

pub type Point2D = Point<2>;
pub type Point3D = Point<3>;

// conversione
impl<const N: usize> From<[f64; N]> for Point<N> {

    fn from(coords: [f64; N]) -> Self {
        Self { coords }
    }

}

impl<const N: usize> Point<N> {

    pub fn coord(&self, i: usize) -> Option<f64> {
        self.coords.get(i).copied()
    }

    pub fn distance(&self, other: &Self) -> f64 {
        self.vector_to(other).norm()
    }

    /// Whether `other` lies within `tolerance` of `self`, Euclidean distance.
    ///
    /// `PartialEq` on `Point` is exact bit equality, which two points that
    /// should count as the same rarely satisfy once they have been through a
    /// reprojection or a round trip through storage.
    pub fn approx_eq(&self, other: &Self, tolerance: f64) -> bool {
        self.distance(other) <= tolerance
    }

    pub fn vector_to(&self, other: &Self) -> Vector<N> {
        Vector::from(other.coords) - Vector::from(self.coords)
    }

    /// The point halfway between `self` and `other`.
    pub fn midpoint(&self, other: &Self) -> Self {
        *self + self.vector_to(other) * 0.5
    }
}

impl<const N: usize> Add<Vector<N>> for Point<N> {
    type Output = Self;

    fn add(self, rhs: Vector<N>) -> Self::Output {
        Self {
            coords: std::array::from_fn(|i| self.coords[i] + rhs.coords[i]),
        }
    }
}

impl<const N: usize> Sub<Vector<N>> for Point<N> {
    type Output = Self;

    fn sub(self, rhs: Vector<N>) -> Self::Output {
        Self {
            coords: std::array::from_fn(|i| self.coords[i] - rhs.coords[i]),
        }
    }
}

impl Point<1> {
    pub fn x(&self) -> f64 { self.coords[0] }
}

impl Point<2> {
    pub fn x(&self) -> f64 {
        self.coords[0]
    }
    pub fn y(&self) -> f64 {
        self.coords[1]
    }

    /// Embeds this point in 3D by giving it the elevation `z`.
    pub fn with_z(&self, z: f64) -> Point3D {
        Point3D::from([self.coords[0], self.coords[1], z])
    }
}

impl Point<3> {
    pub fn x(&self) -> f64 {
        self.coords[0]
    }
    pub fn y(&self) -> f64 {
        self.coords[1]
    }
    pub fn z(&self) -> f64 {
        self.coords[2]
    }

    /// The plan-view projection of this point, elevation dropped.
    pub fn drop_z(&self) -> Point2D {
        Point2D::from([self.coords[0], self.coords[1]])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn point_distance_2d_is_correct() {
        let p1 = Point::<2>::from([0.0, 0.0]);
        let p2 = Point::<2>::from([3.0, 4.0]);

        assert_eq!(p1.distance(&p2), 5.0);
    }

    #[test]
    fn approx_eq_accepts_points_within_tolerance() {
        let p1 = Point::<2>::from([0.0, 0.0]);
        let p2 = Point::<2>::from([1e-10, 0.0]);

        assert!(p1.approx_eq(&p2, 1e-9));
    }

    #[test]
    fn approx_eq_rejects_points_past_tolerance() {
        let p1 = Point::<2>::from([0.0, 0.0]);
        let p2 = Point::<2>::from([1.0, 0.0]);

        assert!(!p1.approx_eq(&p2, 1e-9));
    }

    #[test]
    fn approx_eq_is_exact_at_zero_tolerance() {
        let p1 = Point::<2>::from([0.0, 0.0]);

        assert!(p1.approx_eq(&p1, 0.0));
    }

    #[test]
    fn vector_to_matches_the_coordinatewise_difference() {
        let p1 = Point::<3>::from([1.0, 2.0, 3.0]);
        let p2 = Point::<3>::from([4.0, 0.0, 3.0]);

        assert_eq!(p1.vector_to(&p2), Vector::from([3.0, -2.0, 0.0]));
    }

    #[test]
    fn adding_a_vector_translates_the_point() {
        let p = Point::<2>::from([1.0, 1.0]);
        let v = Vector::from([2.0, -3.0]);

        assert_eq!(p + v, Point::from([3.0, -2.0]));
    }

    #[test]
    fn subtracting_a_vector_is_the_inverse_of_adding_it() {
        let p = Point::<2>::from([1.0, 1.0]);
        let v = Vector::from([2.0, -3.0]);

        assert_eq!(p + v - v, p);
    }

    #[test]
    fn point_plus_vector_to_another_point_reaches_it() {
        let p1 = Point::<3>::from([1.0, 2.0, 3.0]);
        let p2 = Point::<3>::from([4.0, 0.0, 3.0]);

        assert_eq!(p1 + p1.vector_to(&p2), p2);
    }

    #[test]
    fn midpoint_is_equidistant_from_both_points() {
        let p1 = Point::<2>::from([0.0, 0.0]);
        let p2 = Point::<2>::from([4.0, 2.0]);

        let m = p1.midpoint(&p2);

        assert_eq!(m, Point::from([2.0, 1.0]));
        assert_eq!(p1.distance(&m), p2.distance(&m));
    }

    #[test]
    fn midpoint_with_self_is_self() {
        let p = Point::<2>::from([5.0, -1.0]);

        assert_eq!(p.midpoint(&p), p);
    }

    #[test]
    fn drop_z_keeps_the_planar_coordinates() {
        let p = Point3D::from([1.0, 2.0, 3.0]);

        assert_eq!(p.drop_z(), Point2D::from([1.0, 2.0]));
    }

    #[test]
    fn with_z_embeds_a_2d_point_at_the_given_elevation() {
        let p = Point2D::from([1.0, 2.0]);

        assert_eq!(p.with_z(5.0), Point3D::from([1.0, 2.0, 5.0]));
    }

    #[test]
    fn drop_z_then_with_z_round_trips() {
        let p = Point3D::from([1.0, 2.0, 3.0]);

        assert_eq!(p.drop_z().with_z(p.z()), p);
    }
}





