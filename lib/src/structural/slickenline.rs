
use crate::orientation::direction::Direction3D;
use crate::structural::geol_axis::GeologicalAxis;

/// Which way the hanging wall moved, where that is known.
///
/// A qualitative label beside the geometry, not a substitute for it: the
/// direction a `Slickenline` carries already points the way of movement
/// whenever the sense is known at all. This records how a field observation
/// named it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SlipSense {
    Up,
    Down,
    Left,
    Right,
}

/// A slip lineation on a fault plane.
///
/// When the sense of movement is known, `direction` points the way the hanging
/// wall moved and the slickenline is a vector. When it is not, the slip is
/// known only as a line: `direction` is then one of the two antipodal
/// representatives, chosen arbitrarily, and every comparison against it has to
/// be taken modulo 180 degrees -- which is what `angle_to` does, and the
/// reason it exists rather than callers reaching for `Direction3D::angle_deg`.
#[derive(Debug, Clone)]
pub struct Slickenline {
    direction: Direction3D,
    sense: Option<SlipSense>,
}

impl Slickenline {

    /// `sense` is `None` when the movement sense was not determined. There is
    /// deliberately no `SlipSense::Unknown`: `None` already says it, and two
    /// spellings of the same state is one too many.
    pub fn new(direction: Direction3D, sense: Option<SlipSense>) -> Self {
        Self { direction, sense }
    }

    /// From a measured trend and plunge, which is how a slickenline is read in
    /// the field and how a GeoProfiler export stores one.
    pub fn from_axis(axis: GeologicalAxis, sense: Option<SlipSense>) -> Self {
        Self::new(Direction3D::new(axis.as_versor()), sense)
    }

    pub fn direction(&self) -> Direction3D {
        self.direction
    }

    pub fn sense(&self) -> Option<SlipSense> {
        self.sense
    }

    pub fn sense_is_known(&self) -> bool {
        self.sense.is_some()
    }

    /// The trend and plunge this slickenline points along.
    pub fn to_axis(&self) -> GeologicalAxis {
        GeologicalAxis::from_versor(&self.direction.versor())
    }

    /// The angle, in degrees, between this slickenline and a direction.
    ///
    /// In `[0, 180]` when the sense is known, and in `[0, 90]` when it is not:
    /// an undirected slip is as well matched by a prediction pointing the
    /// other way, so the supplementary angle is taken whenever it is smaller.
    /// Getting this wrong would cost a stress inversion every fault whose
    /// sense was not read, scoring a perfect antiparallel fit as the worst
    /// possible one.
    pub fn angle_to(&self, other: &Direction3D) -> f64 {

        let angle = self.direction.angle_deg(other);

        if self.sense_is_known() {
            angle
        } else {
            angle.min(180.0 - angle)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn along(trend: f64, plunge: f64) -> Direction3D {
        Direction3D::new(GeologicalAxis::new(trend, plunge).as_versor())
    }

    #[test]
    fn from_axis_round_trips_through_to_axis() {
        for &(trend, plunge) in &[(0.0, 0.0), (90.0, 45.0), (217.0, -30.0)] {
            let slick = Slickenline::from_axis(GeologicalAxis::new(trend, plunge), None);
            let back = slick.to_axis();

            assert!((back.trend - trend).abs() < 1e-9);
            assert!((back.plunge - plunge).abs() < 1e-9);
        }
    }

    #[test]
    fn a_known_sense_is_reported_as_known() {
        let slick = Slickenline::from_axis(GeologicalAxis::new(90.0, 30.0), Some(SlipSense::Down));

        assert!(slick.sense_is_known());
        assert_eq!(slick.sense(), Some(SlipSense::Down));
    }

    #[test]
    fn angle_to_the_same_direction_is_zero_either_way() {
        let axis = GeologicalAxis::new(90.0, 30.0);

        for sense in [None, Some(SlipSense::Down)] {
            let slick = Slickenline::from_axis(axis, sense);
            assert!(slick.angle_to(&along(90.0, 30.0)) < 1e-9);
        }
    }

    #[test]
    fn an_unknown_sense_matches_the_opposite_direction_just_as_well() {
        // The distinction the modulo exists for: with the sense unread, a
        // prediction pointing the other way along the same line is a perfect
        // fit, not the worst one.
        let slick = Slickenline::from_axis(GeologicalAxis::new(90.0, 30.0), None);
        let opposite = along(270.0, -30.0);

        assert!(slick.angle_to(&opposite) < 1e-9);
    }

    #[test]
    fn a_known_sense_does_tell_the_two_apart() {
        let slick =
            Slickenline::from_axis(GeologicalAxis::new(90.0, 30.0), Some(SlipSense::Down));
        let opposite = along(270.0, -30.0);

        assert!((slick.angle_to(&opposite) - 180.0).abs() < 1e-9);
    }

    #[test]
    fn an_unknown_sense_never_reports_more_than_a_right_angle() {
        let slick = Slickenline::from_axis(GeologicalAxis::new(0.0, 0.0), None);

        for trend in [0.0, 30.0, 89.0, 91.0, 150.0, 180.0, 270.0] {
            let angle = slick.angle_to(&along(trend, 0.0));
            assert!(
                (-1e-9..=90.0 + 1e-9).contains(&angle),
                "trend {}: {} degrees",
                trend,
                angle
            );
        }
    }
}
