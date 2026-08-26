
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

    /// Whether `self` and `other` describe the same infinite line, regardless
    /// of which point each uses as its origin or which way its direction
    /// points -- `Line`'s derived `PartialEq` compares those choices, not the
    /// line itself, so two representations of the same line answer `false`
    /// to `==`.
    ///
    /// `angle_tolerance` bounds how far from parallel (or anti-parallel) the
    /// two directions may be, as the sine of the angle between them: both
    /// directions are already unit (`Versor`), so `sin(angle) = sqrt(1 -
    /// dot^2)` needs no cross product and so works at any `N`, unlike the
    /// collinearity check in `Plane::from_three_points`, which is 3D-only.
    /// `distance_tolerance` bounds how far `other.origin` may sit from `self`.
    pub fn is_same_line_as(&self, other: &Self, angle_tolerance: f64, distance_tolerance: f64) -> bool {

        let dot = self.direction.dot(&other.direction);
        let sin_angle = (1.0 - dot * dot).max(0.0).sqrt();

        if sin_angle > angle_tolerance {
            return false;
        }

        let t = self.origin.vector_to(&other.origin).dot(&self.direction.as_vector());
        let closest = self.point_at(t);

        closest.distance(&other.origin) <= distance_tolerance
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::point::{Point, Point3D};

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

    fn x_axis() -> Line3D {
        Line3D::new(Point3D::from([0.0, 0.0, 0.0]), Versor::new([1.0, 0.0, 0.0]).unwrap()).unwrap()
    }

    #[test]
    fn same_line_with_a_different_origin_is_recognized() {
        let other = Line3D::new(Point3D::from([5.0, 0.0, 0.0]), Versor::new([1.0, 0.0, 0.0]).unwrap())
            .unwrap();

        assert!(x_axis().is_same_line_as(&other, 1e-9, 1e-9));
    }

    #[test]
    fn same_line_with_opposite_direction_is_recognized() {
        let other = Line3D::new(Point3D::from([3.0, 0.0, 0.0]), Versor::new([-1.0, 0.0, 0.0]).unwrap())
            .unwrap();

        assert!(x_axis().is_same_line_as(&other, 1e-9, 1e-9));
    }

    #[test]
    fn parallel_offset_lines_are_not_the_same_line() {
        let other = Line3D::new(Point3D::from([0.0, 1.0, 0.0]), Versor::new([1.0, 0.0, 0.0]).unwrap())
            .unwrap();

        assert!(!x_axis().is_same_line_as(&other, 1e-9, 1e-9));
    }

    #[test]
    fn crossing_lines_with_different_directions_are_not_the_same_line() {
        let other = Line3D::new(Point3D::from([0.0, 0.0, 0.0]), Versor::new([0.0, 1.0, 0.0]).unwrap())
            .unwrap();

        assert!(!x_axis().is_same_line_as(&other, 1e-9, 1e-9));
    }

    #[test]
    fn skew_lines_are_not_the_same_line() {
        let other = Line3D::new(Point3D::from([0.0, 0.0, 1.0]), Versor::new([0.0, 1.0, 0.0]).unwrap())
            .unwrap();

        assert!(!x_axis().is_same_line_as(&other, 1e-9, 1e-9));
    }

    #[test]
    fn offset_within_distance_tolerance_is_accepted() {
        let other = Line3D::new(Point3D::from([0.0, 1e-10, 0.0]), Versor::new([1.0, 0.0, 0.0]).unwrap())
            .unwrap();

        assert!(x_axis().is_same_line_as(&other, 1e-9, 1e-9));
    }

    #[test]
    fn offset_past_distance_tolerance_is_rejected() {
        let other = Line3D::new(Point3D::from([0.0, 1e-3, 0.0]), Versor::new([1.0, 0.0, 0.0]).unwrap())
            .unwrap();

        assert!(!x_axis().is_same_line_as(&other, 1e-9, 1e-9));
    }

    #[test]
    fn tilt_within_angle_tolerance_is_accepted() {
        let other = Line3D::new(Point3D::from([0.0, 0.0, 0.0]), Versor::new([1.0, 1e-10, 0.0]).unwrap())
            .unwrap();

        assert!(x_axis().is_same_line_as(&other, 1e-9, 1e-9));
    }

    #[test]
    fn tilt_past_angle_tolerance_is_rejected() {
        let other = Line3D::new(Point3D::from([0.0, 0.0, 0.0]), Versor::new([1.0, 1e-3, 0.0]).unwrap())
            .unwrap();

        assert!(!x_axis().is_same_line_as(&other, 1e-9, 1e-9));
    }
}
