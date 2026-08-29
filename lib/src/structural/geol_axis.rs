
use crate::algebra::versor::Versor3D;

/// Structural axis, in trend/plunge: trend in `[0, 360)` clockwise from North,
/// plunge in `[-90, 90]`, negative upward.
///
/// `new` stores what it is given, folding neither range -- unlike a full
/// stereographic `Axis`, whose antipodal identification (a trend/plunge pair
/// and its `trend + 180, -plunge` twin name the same line) this type does not
/// carry. `as_versor` and `from_versor` need no such folding to round-trip: a
/// versor already picks one representative, and the one `from_versor` returns
/// is always in the two ranges above.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeologicalAxis {
    pub trend: f64,
    pub plunge: f64
}

impl GeologicalAxis {

    pub fn new(tr: f64, pl: f64) -> Self { GeologicalAxis{ trend: tr, plunge: pl }}

    /// The unit vector this trend/plunge points along, in (East, North, Up).
    ///
    /// Same construction as `GeologicalPlane::normal_vector`, at plunge
    /// `dip_angle - 90` from a plane's dip azimuth/angle -- the two are not
    /// independent formulas kept in sync by hand, they are the same trig
    /// applied to different angles.
    pub fn as_versor(&self) -> Versor3D {

        let trend = self.trend.to_radians();
        let plunge = self.plunge.to_radians();

        let coords = [
            trend.sin() * plunge.cos(),
            trend.cos() * plunge.cos(),
            -plunge.sin(),
        ];

        Versor3D::new(coords)
            .expect("sin^2 + cos^2 = 1: never the zero vector, whatever trend and plunge are")
    }

    /// The trend/plunge a unit vector points along: the inverse of `as_versor`.
    ///
    /// Always the canonical pair -- trend in `[0, 360)`, plunge in `[-90, 90]`
    /// -- regardless of what produced the versor.
    pub fn from_versor(v: &Versor3D) -> Self {

        let &[east, north, up] = v.coords();

        let plunge = (-up).clamp(-1.0, 1.0).asin().to_degrees();

        let trend = east.atan2(north).to_degrees();
        let trend = if trend < 0.0 { trend + 360.0 } else { trend };

        Self::new(trend, plunge)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn north_horizontal_is_the_north_versor() {
        let v = GeologicalAxis::new(0.0, 0.0).as_versor();
        let c = v.coords();

        assert!(c[0].abs() < 1e-12, "east: {:?}", c);
        assert!((c[1] - 1.0).abs() < 1e-12, "north: {:?}", c);
        assert!(c[2].abs() < 1e-12, "up: {:?}", c);
    }

    #[test]
    fn straight_down_has_no_trend_dependence() {
        for trend in [0.0, 90.0, 217.0] {
            let v = GeologicalAxis::new(trend, 90.0).as_versor();
            let c = v.coords();
            assert!((c[2] - (-1.0)).abs() < 1e-12, "trend {}: {:?}", trend, c);
        }
    }

    #[test]
    fn as_versor_then_from_versor_round_trips() {
        for &(trend, plunge) in &[(0.0, 0.0), (90.0, 45.0), (217.0, -30.0), (359.0, 89.0)] {
            let axis = GeologicalAxis::new(trend, plunge);
            let back = GeologicalAxis::from_versor(&axis.as_versor());

            assert!((back.trend - trend).abs() < 1e-9, "trend: {} vs {}", back.trend, trend);
            assert!((back.plunge - plunge).abs() < 1e-9, "plunge: {} vs {}", back.plunge, plunge);
        }
    }

    #[test]
    fn from_versor_of_straight_up_reports_zero_trend() {
        // Trend is undefined at the poles; atan2(0, 0) settles on 0, matching
        // GeologicalPlane::from_plane's own convention at its analogous
        // degenerate case.
        let axis = GeologicalAxis::from_versor(&Versor3D::new([0.0, 0.0, 1.0]).unwrap());

        assert_eq!(axis.trend, 0.0);
        assert!((axis.plunge - (-90.0)).abs() < 1e-9);
    }
}
