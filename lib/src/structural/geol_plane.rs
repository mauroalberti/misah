
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

    /// The attitude of a geometric plane: the inverse of `to_plane`.
    ///
    /// A `Plane`'s normal has no inherent up/down sense -- `from_three_points`
    /// gives one or the other depending only on the order the points were
    /// listed in, and both describe the same physical surface. This picks
    /// whichever of `plane.normal`/`-plane.normal` points upward before
    /// reading azimuth and dip off it, matching the convention
    /// `normal_vector` always produces; for a horizontal plane the two
    /// coincide and the dip direction is undefined, reported here as 0
    /// (`atan2`'s own answer at the origin). For an exactly vertical plane
    /// the upward pick is itself undefined (the z-component is ~0 either
    /// way), so the reported azimuth there depends on which of the two
    /// normals `plane` happened to carry.
    pub fn from_plane(plane: &Plane) -> Self {
        let normal = plane.normal;
        let n = if normal.coords[2] < 0.0 { -normal } else { normal };

        let dip_angle = n.coords[2].clamp(-1.0, 1.0).acos().to_degrees();

        let az = n.coords[0].atan2(n.coords[1]).to_degrees();
        let azimuth = if az < 0.0 { az + 360.0 } else { az };

        GeologicalPlane { azimuth, dip_angle }
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

    #[test]
    fn from_plane_recovers_the_attitude_that_built_it() {
        for &(az, dip) in &[(135.0, 35.0), (0.0, 45.0), (270.0, 60.0), (359.0, 10.0)] {
            let plane = GeologicalPlane::new(az, dip)
                .to_plane(Point3D::from([1.0, 2.0, 3.0]))
                .unwrap();

            let recovered = GeologicalPlane::from_plane(&plane);

            assert!(
                (recovered.azimuth - az).abs() < 1e-9,
                "azimuth: expected {az}, got {}",
                recovered.azimuth
            );
            assert!(
                (recovered.dip_angle - dip).abs() < 1e-9,
                "dip: expected {dip}, got {}",
                recovered.dip_angle
            );
        }
    }

    #[test]
    fn from_plane_of_a_horizontal_plane_has_zero_dip() {
        let plane = GeologicalPlane::new(200.0, 0.0)
            .to_plane(Point3D::from([0.0, 0.0, 0.0]))
            .unwrap();

        // The dip direction of a horizontal plane is undefined; only the dip
        // angle -- the part that actually means something -- is checked.
        assert!(GeologicalPlane::from_plane(&plane).dip_angle.abs() < 1e-9);
    }

    #[test]
    fn from_plane_ignores_which_way_the_normal_happened_to_point() {
        let upward = GeologicalPlane::new(135.0, 35.0).normal_vector();
        let point = Point3D::from([5.0, -5.0, 2.0]);

        let from_upward = GeologicalPlane::from_plane(&Plane::new(point, upward).unwrap());
        let from_downward = GeologicalPlane::from_plane(&Plane::new(point, -upward).unwrap());

        assert!((from_upward.azimuth - from_downward.azimuth).abs() < 1e-9);
        assert!((from_upward.dip_angle - from_downward.dip_angle).abs() < 1e-9);
    }
}
