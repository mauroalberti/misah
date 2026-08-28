
use super::geol_axis::{GeologicalAxis};

use crate::algebra::Vector;
use crate::algebra::versor::Versor3D;
use crate::geometry::plane::Plane;
use crate::geometry::point::Point3D;

#[derive(Debug, Clone)]
pub struct GeologicalPlane {
    pub azimuth: f64,
    pub dip_angle: f64
}

impl GeologicalPlane {

    pub fn new(az: f64, dip: f64) -> Self { GeologicalPlane{ azimuth: az, dip_angle: dip }}

    /// Build a plane from its strike by the right-hand rule rather than its
    /// dip azimuth: walking `strike` with the dip on the right hand reaches
    /// `strike + 90`, which is what this stores as `azimuth`.
    pub fn from_rhr_strike(strike: f64, dip: f64) -> Self {
        Self::new((strike + 90.0).rem_euclid(360.0), dip)
    }

    /// The strike, by the right-hand rule: the inverse of `from_rhr_strike`.
    pub fn rhr_strike(&self) -> f64 {
        (self.azimuth - 90.0).rem_euclid(360.0)
    }

    /// The direction a slickenline of the given rake points on this plane,
    /// Aki & Richards (1980)'s convention: rake 0 is left-lateral, 90
    /// reverse, +/-180 right-lateral, -90 normal.
    ///
    /// Unit length is an algebraic identity of the formula (the strike terms
    /// cancel by `sin^2 + cos^2 = 1`, then so do the dip ones), holding for
    /// every strike, dip and rake -- there is no input this can fail on.
    pub fn rake_to_versor(&self, rake_degrees: f64) -> Versor3D {

        let strike = self.rhr_strike().to_radians();
        let dip = self.dip_angle.to_radians();
        let rake = rake_degrees.to_radians();

        let coords = [
            rake.cos() * strike.sin() - rake.sin() * dip.cos() * strike.cos(),
            rake.cos() * strike.cos() + rake.sin() * dip.cos() * strike.sin(),
            rake.sin() * dip.sin(),
        ];

        Versor3D::new(coords)
            .expect("a unit vector for any strike, dip and rake, by construction")
    }

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
    fn rhr_strike_and_from_rhr_strike_round_trip() {
        for &(strike, dip) in &[(0.0, 60.0), (30.0, 70.0), (200.0, 50.0), (350.0, 10.0)] {
            let plane = GeologicalPlane::from_rhr_strike(strike, dip);
            assert!((plane.rhr_strike() - strike).abs() < 1e-9);
            assert_eq!(plane.dip_angle, dip);
        }
    }

    #[test]
    fn rhr_strike_is_ninety_degrees_behind_the_dip_azimuth() {
        // strike + 90 = dip azimuth, wrapped into [0, 360) rather than left
        // negative -- the case a plain Rust `%` would get wrong.
        assert_eq!(GeologicalPlane::new(45.0, 30.0).rhr_strike(), 315.0);
    }

    #[test]
    fn a_rake_of_minus_90_is_pure_normal_dip_slip() {
        // Rake -90 (Aki & Richards) points straight down the dip vector,
        // whatever the strike: the classic Andersonian normal-fault check.
        let plane = GeologicalPlane::from_rhr_strike(0.0, 60.0);

        let slick = plane.rake_to_versor(-90.0);
        let dip_vector = GeologicalAxis::new(plane.azimuth, plane.dip_angle).as_versor();

        for i in 0..3 {
            assert!(
                (slick.coords()[i] - dip_vector.coords()[i]).abs() < 1e-12,
                "rake -90 should coincide with the dip vector"
            );
        }
    }

    #[test]
    fn rake_to_versor_is_always_unit_length() {
        for &(strike, dip, rake) in &[
            (0.0, 60.0, -90.0), (30.0, 70.0, -2.26), (200.0, 50.0, 43.56), (10.0, 0.5, 179.0),
        ] {
            let v = GeologicalPlane::from_rhr_strike(strike, dip).rake_to_versor(rake);
            let norm_sq: f64 = v.coords().iter().map(|c| c * c).sum();
            assert!((norm_sq - 1.0).abs() < 1e-12, "strike {} dip {} rake {}", strike, dip, rake);
        }
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
