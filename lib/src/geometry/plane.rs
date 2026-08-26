
use super::point::Point3D;
use crate::algebra::constants::EPSILON;
use crate::algebra::vector::Vector3D;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Plane {
    pub point: Point3D,
    pub normal: Vector3D,
}

impl Plane {

    pub fn new(point: Point3D, normal: Vector3D) -> Option<Self> {

        let normal = normal.normalize()?;
        Some(Self { point, normal })
    }

    pub fn from_three_points(
        p0: Point3D,
        p1: Point3D,
        p2: Point3D,
    ) -> Option<Self> {

        let v1 = p0.vector_to(&p1);
        let v2 = p0.vector_to(&p2);

        let len1 = v1.norm();
        let len2 = v2.norm();

        // Two of the three points coincide.
        if len1 < EPSILON || len2 < EPSILON {
            return None;
        }

        let cross = v1.cross(&v2);

        // |cross| / (len1 * len2) is the sine of the angle between v1 and v2: a
        // dimensionless, scale-invariant measure of collinearity. Thresholding
        // the cross product itself would not do, since its magnitude scales
        // with len1 * len2 -- a fixed cutoff would accept anything at the UTM
        // scale (eastings and northings in the 1e5-1e7 range) while rejecting
        // valid triangles a metre across. Points that are nearly but not
        // exactly collinear pass this check as-is: without it, a normal whose
        // direction is mostly rounding noise was returned as if it meant
        // something.
        if cross.norm() <= EPSILON * len1 * len2 {
            return None;
        }

        let normal = cross.normalize()?;

        Some(Self {
            point: p0,
            normal,
        })
    }

    pub fn signed_distance_to_point(&self, point: &Point3D) -> f64 {

        let v = self.point.vector_to(point);
        self.normal.dot(&v)
    }

    pub fn distance_to_point(&self, point: &Point3D) -> f64 {
        self.signed_distance_to_point(point).abs()
    }

    pub fn contains_point(&self, point: &Point3D, tolerance: f64) -> bool {
        self.distance_to_point(point) <= tolerance
    }

    pub fn coefficients(&self) -> (f64, f64, f64, f64) {

        let [a, b, c] = self.normal.coords;
        let [x0, y0, z0] = self.point.coords;

        let d = -(a * x0 + b * y0 + c * z0);

        (a, b, c, d)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn three_well_spread_points_give_a_valid_plane() {
        let plane = Plane::from_three_points(
            Point3D::from([0.0, 0.0, 0.0]),
            Point3D::from([1.0, 0.0, 0.0]),
            Point3D::from([0.0, 1.0, 0.0]),
        )
        .expect("a well-conditioned triangle");

        assert!((plane.normal.norm() - 1.0).abs() < 1e-12);
        assert!(plane.contains_point(&Point3D::from([1.0, 0.0, 0.0]), 1e-12));
    }

    #[test]
    fn exactly_collinear_points_are_rejected() {
        let plane = Plane::from_three_points(
            Point3D::from([0.0, 0.0, 0.0]),
            Point3D::from([1.0, 0.0, 0.0]),
            Point3D::from([2.0, 0.0, 0.0]),
        );

        assert!(plane.is_none());
    }

    #[test]
    fn nearly_collinear_points_are_rejected() {
        // p2 sits a fraction of a nanometre off the line through p0-p1: not bit-
        // identical to collinear, but well past what any real survey, or the
        // float arithmetic computing it, could tell apart from a straight line.
        // Before the tolerance check this returned Some(plane) with a normal
        // that was mostly rounding noise.
        let plane = Plane::from_three_points(
            Point3D::from([0.0, 0.0, 0.0]),
            Point3D::from([1.0, 0.0, 0.0]),
            Point3D::from([2.0, 1e-10, 0.0]),
        );

        assert!(plane.is_none());
    }

    #[test]
    fn a_small_but_real_angle_still_gives_a_plane() {
        // About 0.1 degrees between the two edges: a shallow, real-world sliver
        // of a triangle, not a degenerate one. The fix must not reject this.
        let angle = 0.1_f64.to_radians();
        let plane = Plane::from_three_points(
            Point3D::from([0.0, 0.0, 0.0]),
            Point3D::from([1.0, 0.0, 0.0]),
            Point3D::from([angle.cos(), angle.sin(), 0.0]),
        );

        assert!(plane.is_some());
    }

    #[test]
    fn coincident_points_are_rejected() {
        let plane = Plane::from_three_points(
            Point3D::from([0.0, 0.0, 0.0]),
            Point3D::from([0.0, 0.0, 0.0]),
            Point3D::from([1.0, 1.0, 0.0]),
        );

        assert!(plane.is_none());
    }
}
