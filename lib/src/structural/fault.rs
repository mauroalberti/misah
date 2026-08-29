
use thiserror::Error;

use super::geol_plane::GeologicalPlane;
use super::slickenline::Slickenline;

/// How far, in degrees, a slickenline may sit off its fault plane and still be
/// accepted.
///
/// Field measurements of a plane and of a lineation on it are taken
/// separately, so they never agree exactly; a degree is the order of what a
/// compass-clinometer reading is worth. Rejecting anything closer would refuse
/// real data, and accepting anything further would let a transcription error
/// through as geometry.
pub const SLICKENLINE_COPLANARITY_TOLERANCE: f64 = 1.0;

#[derive(Debug, Error)]
pub enum FaultError {
    #[error(
        "slickenline {index} lies {departure:.3} degrees off the fault plane, \
         more than the {tolerance:.3} allowed"
    )]
    SlickenlineNotInPlane { index: usize, departure: f64, tolerance: f64 },
}

/// A fault plane and the slip lineations observed on it.
///
/// The pair is what a stress inversion consumes: the plane says which
/// orientation the stress tensor is resolved onto, the slickenlines say what
/// the rock actually did, and the misfit between predicted and observed slip
/// is what such an inversion minimises. `ReducedStressTensor::misfit_on` takes
/// one of these.
#[derive(Debug, Clone)]
pub struct FaultPlane {
    plane: GeologicalPlane,
    slickenlines: Vec<Slickenline>,
}

impl FaultPlane {

    /// Every slickenline is checked against the plane it is said to lie on.
    ///
    /// A lineation off the plane is not a fault with a small error in it; it
    /// is two measurements that do not belong to each other, and the forward
    /// model would happily resolve a stress tensor onto the plane and compare
    /// the result with a direction the fault could never have slipped along.
    pub fn new(
        plane: GeologicalPlane,
        slickenlines: Vec<Slickenline>,
    ) -> Result<Self, FaultError> {

        let normal = plane.normal_vector();

        for (index, slickenline) in slickenlines.iter().enumerate() {

            // The slickenline lies in the plane when it is perpendicular to
            // the normal, so the departure is how far the angle between them
            // falls from a right angle.
            let cosine = normal.dot(&slickenline.direction().as_vector()).clamp(-1.0, 1.0);
            let departure = (cosine.acos().to_degrees() - 90.0).abs();

            if departure > SLICKENLINE_COPLANARITY_TOLERANCE {
                return Err(FaultError::SlickenlineNotInPlane {
                    index,
                    departure,
                    tolerance: SLICKENLINE_COPLANARITY_TOLERANCE,
                });
            }
        }

        Ok(Self { plane, slickenlines })
    }

    /// A fault plane with no slip observed on it.
    ///
    /// Of no use to an inversion, which has nothing to score it against, but
    /// the forward model still predicts what it would do under a given stress.
    pub fn without_slickenlines(plane: GeologicalPlane) -> Self {
        Self { plane, slickenlines: Vec::new() }
    }

    /// The plane, with the given rake as its single slickenline.
    ///
    /// The rake is Aki & Richards (1980)'s, so the lineation is built by the
    /// same formula the forward model predicts one with, and cannot fail the
    /// coplanarity check: `rake_to_versor` returns a direction lying in the
    /// plane by construction.
    pub fn from_rake(
        plane: GeologicalPlane,
        rake_degrees: f64,
        sense: Option<super::slickenline::SlipSense>,
    ) -> Self {
        let direction =
            crate::orientation::direction::Direction3D::new(plane.rake_to_versor(rake_degrees));

        Self {
            plane,
            slickenlines: vec![Slickenline::new(direction, sense)],
        }
    }

    pub fn plane(&self) -> &GeologicalPlane {
        &self.plane
    }

    pub fn slickenlines(&self) -> &[Slickenline] {
        &self.slickenlines
    }

    pub fn slickenline(&self, index: usize) -> Option<&Slickenline> {
        self.slickenlines.get(index)
    }

    pub fn has_slickenlines(&self) -> bool {
        !self.slickenlines.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::orientation::direction::Direction3D;
    use crate::structural::geol_axis::GeologicalAxis;
    use crate::structural::slickenline::SlipSense;

    fn along(trend: f64, plunge: f64) -> Direction3D {
        Direction3D::new(GeologicalAxis::new(trend, plunge).as_versor())
    }

    #[test]
    fn a_dip_slip_lineation_lies_in_its_plane() {
        // A plane dipping 60 towards 090, and the slip straight down its dip.
        let plane = GeologicalPlane::new(90.0, 60.0);
        let slick = Slickenline::new(along(90.0, 60.0), Some(SlipSense::Down));

        assert!(FaultPlane::new(plane, vec![slick]).is_ok());
    }

    #[test]
    fn a_strike_slip_lineation_lies_in_its_plane_too() {
        // Horizontal slip along the strike of the same plane.
        let plane = GeologicalPlane::new(90.0, 60.0);
        let slick = Slickenline::new(along(0.0, 0.0), Some(SlipSense::Right));

        assert!(FaultPlane::new(plane, vec![slick]).is_ok());
    }

    #[test]
    fn a_lineation_off_the_plane_is_refused() {
        // Pointing along the plane's own normal: as far off as a lineation
        // can be, and exactly what a mismatched pair of measurements looks
        // like.
        let plane = GeologicalPlane::new(90.0, 60.0);
        let slick = Slickenline::new(along(270.0, 30.0), None);

        let err = FaultPlane::new(plane, vec![slick]).unwrap_err();

        assert!(
            matches!(err, FaultError::SlickenlineNotInPlane { index: 0, .. }),
            "unexpected error: {}",
            err
        );
    }

    #[test]
    fn the_offending_slickenline_is_named() {
        let plane = GeologicalPlane::new(90.0, 60.0);
        let good = Slickenline::new(along(90.0, 60.0), None);
        let bad = Slickenline::new(along(270.0, 30.0), None);

        let err = FaultPlane::new(plane, vec![good, bad]).unwrap_err();

        assert!(matches!(err, FaultError::SlickenlineNotInPlane { index: 1, .. }));
    }

    #[test]
    fn a_lineation_within_measurement_error_is_accepted() {
        // Half a degree off: two instrument readings that do not quite agree,
        // which is every real pair of them.
        let plane = GeologicalPlane::new(90.0, 60.0);
        let slick = Slickenline::new(along(90.0, 60.5), None);

        assert!(FaultPlane::new(plane, vec![slick]).is_ok());
    }

    #[test]
    fn from_rake_builds_a_fault_that_passes_its_own_check() {
        // rake_to_versor lies in the plane by construction, so every rake must
        // survive the coplanarity test new() applies.
        for rake in [-180.0, -90.0, -45.0, 0.0, 30.0, 90.0, 179.0] {
            let plane = GeologicalPlane::from_rhr_strike(30.0, 70.0);
            let fault = FaultPlane::from_rake(plane.clone(), rake, None);

            let rebuilt = FaultPlane::new(plane, fault.slickenlines().to_vec());
            assert!(rebuilt.is_ok(), "rake {} was refused", rake);
        }
    }

    #[test]
    fn a_fault_without_slickenlines_says_so() {
        let fault = FaultPlane::without_slickenlines(GeologicalPlane::new(90.0, 60.0));

        assert!(!fault.has_slickenlines());
        assert!(fault.slickenline(0).is_none());
    }
}
