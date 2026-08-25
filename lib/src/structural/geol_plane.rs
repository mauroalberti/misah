
use super::geol_axis::{GeologicalAxis};

use crate::algebra::Vector;
use crate::geometry::plane::Plane;
use crate::geometry::point::Point3D;

#[derive(Debug, Clone)]
pub struct GeologicalPlane {
    pub azimuth: f64,
    pub dip_angle: f64
}

impl GeologicalPlane {

    pub fn new(az: f64, dip: f64) -> Self { GeologicalPlane{ azimuth: az, dip_angle: dip }}

    pub fn normal_axis(&self) -> GeologicalAxis {

        GeologicalAxis{
            trend: (self.azimuth + 180.0) % 360.0,
            plunge: 90.0 - self.dip_angle

        }
    }

    /// Upward-pointing unit normal, in (East, North, Up).
    ///
    /// Sanity anchors: a horizontal plane gives (0, 0, 1); a vertical one gives a
    /// horizontal normal pointing along the dip direction.
    pub fn normal_vector(&self) -> Vector<3> {

        let az = self.azimuth.to_radians();
        let dip = self.dip_angle.to_radians();

        Vector::new([
            az.sin() * dip.sin(),
            az.cos() * dip.sin(),
            dip.cos(),
        ])
    }

    /// The geometric plane with this attitude passing through `point`.
    ///
    /// `None` only if the normal degenerates, which an attitude cannot produce;
    /// the signature follows `Plane::new`.
    pub fn to_plane(&self, point: Point3D) -> Option<Plane> {
        Plane::new(point, self.normal_vector())
    }

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normal_of_a_horizontal_plane_points_up() {
        let n = GeologicalPlane::new(0.0, 0.0).normal_vector().coords;

        assert!(n[0].abs() < 1e-12);
        assert!(n[1].abs() < 1e-12);
        assert!((n[2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn normal_of_a_vertical_plane_lies_horizontal() {
        let n = GeologicalPlane::new(90.0, 90.0).normal_vector().coords;

        assert!((n[0] - 1.0).abs() < 1e-12);
        assert!(n[2].abs() < 1e-12);
    }

    #[test]
    fn normal_is_a_unit_vector() {
        let n = GeologicalPlane::new(135.0, 35.0).normal_vector();

        assert!((n.norm() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn plane_passes_through_its_point() {
        let point = Point3D::from([10.0, 20.0, 30.0]);
        let plane = GeologicalPlane::new(135.0, 35.0).to_plane(point).unwrap();

        assert!(plane.signed_distance_to_point(&point).abs() < 1e-12);
    }
}
