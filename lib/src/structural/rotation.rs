//! Rotations between focal mechanisms, after Kagan (1991).
//!
//! ## Why there are four answers and not one
//!
//! A double-couple source is unchanged by a half turn about any of its own P,
//! T or B axes: turn the mechanism over and it radiates the same thing. So the
//! question "what rotation carries mechanism 1 onto mechanism 2" has **four**
//! answers, one for each of those symmetries composed with the base rotation,
//! and all four are equally correct descriptions of the same pair.
//!
//! The one usually wanted is the smallest, and it has a name: the **Kagan
//! angle**, the standard measure of how far apart two focal mechanisms are. It
//! is the first element of what `focal_mechanism_rotations` returns, the four
//! being sorted by the size of their turn. The other three are kept because
//! their spread is informative in itself -- and because a minimum reported
//! without saying it was a minimum over four is a number nobody can check.
//!
//! ## Sign, and the shorter way round
//!
//! An axis and its opposite end describe the same rotation, one turning by `a`
//! and the other by `-a`; and a turn of `a` is a turn of `a - 360` the other
//! way. Neither is more true than the other. What is settled here is only that
//! the angle reported is the **shorter** of the two ways round, so it lies in
//! `[-180, 180]` -- see `RotationAxis::to_minimum`. The axis is then whatever
//! that choice implies, and may point up or down.

use crate::algebra::quaternion::{QUATERNION_MAGNITUDE_THRESHOLD, Quaternion};
use crate::algebra::versor::Versor3D;
use crate::structural::focal_mechanism::{FocalMechanismError, PTBAxes};
use crate::structural::geol_axis::GeologicalAxis;

/// An axis to turn about, and how far to turn.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RotationAxis {
    /// Where the axis points, as trend and plunge in degrees. The plunge may
    /// be negative: this axis is directed, and its direction is what fixes the
    /// sign of the angle.
    pub axis: GeologicalAxis,
    /// The turn about it, in degrees, right-handed about the axis direction.
    pub angle_degrees: f64,
}

impl RotationAxis {

    pub fn new(trend: f64, plunge: f64, angle_degrees: f64) -> Self {
        Self { axis: GeologicalAxis::new(trend, plunge), angle_degrees }
    }

    /// The rotation a quaternion stands for.
    ///
    /// `None` for a quaternion too small to describe an orientation at all. A
    /// quaternion that is unit but has no vector part is *not* that case: it is
    /// the null rotation, and comes back as an angle of zero about an axis
    /// reported as north-horizontal, which is **arbitrary** -- turning by
    /// nothing has no axis, and any direction serves as well as another.
    pub fn from_quaternion(quaternion: &Quaternion) -> Option<Self> {

        let unit = quaternion.normalized()?;
        let angle = unit.rotation_angle_degrees()?;

        // Below the threshold the vector part is rounding rather than a
        // direction, which is the null rotation reached from the scalar side.
        if unit.vector.norm() < QUATERNION_MAGNITUDE_THRESHOLD {
            return Some(Self::new(0.0, 0.0, 0.0));
        }

        let direction = unit.vector.normalize()?;
        let versor = Versor3D::try_from(direction).ok()?;

        Some(Self { axis: GeologicalAxis::from_versor(&versor), angle_degrees: angle })
    }

    /// The unit quaternion of this rotation: `cos(a/2)`, with `sin(a/2)` along
    /// the axis.
    pub fn to_quaternion(&self) -> Quaternion {

        let half = self.angle_degrees.to_radians() / 2.0;

        Quaternion::new(half.cos(), self.axis.as_versor().as_vector() * half.sin())
    }

    /// The same rotation described from the other end of the axis.
    ///
    /// Turning by `a` about an axis is turning by `360 - a` about its
    /// opposite: the same movement, written the other way round.
    pub fn specular(&self) -> Self {

        let opposite = self.axis.as_versor().opposite();
        let angle = 360.0 - self.angle_degrees;
        // Back into (-360, 360], so that a negative angle specularises to a
        // positive one of the complementary size rather than to one past a
        // full turn.
        let angle = if angle > 360.0 { angle - 360.0 } else { angle };

        Self { axis: GeologicalAxis::from_versor(&opposite), angle_degrees: angle }
    }

    /// The same rotation, taking the shorter way round: `|angle| <= 180`.
    pub fn to_minimum(&self) -> Self {
        if self.angle_degrees.abs() <= 180.0 {
            *self
        } else {
            self.specular()
        }
    }
}

/// The four rotations that carry `from` onto `to`, smallest turn first.
///
/// The first is the **Kagan angle** between the two mechanisms, which is what
/// most callers want; the rest are the same rotation composed with the double
/// couple's own symmetries, and are returned rather than discarded because the
/// spread among them is part of the answer.
///
/// Identical mechanisms need no special case here. The base rotation comes out
/// as nothing, and the other three come out as the half turns about P, T and B
/// that a double couple is genuinely invariant under -- which is the honest
/// answer, and more use than the empty list geogst returns or the four zeroes
/// the Fortran does.
pub fn focal_mechanism_rotations(
    from: &PTBAxes,
    to: &PTBAxes,
) -> Result<[RotationAxis; 4], FocalMechanismError> {

    let first = from.to_quaternion();
    let second = to.to_quaternion();

    let first_inverse = first.inverse().ok_or(FocalMechanismError::DegenerateQuaternion)?;
    let second_inverse = second.inverse().ok_or(FocalMechanismError::DegenerateQuaternion)?;

    // q' = q2 * q1^-1, the rotation from the first frame to the second.
    let base = second * first_inverse;

    // The three half turns about the second mechanism's own axes, carried into
    // the reference frame as a(i) = q2 * i * q2^-1. Composed with the base
    // rotation they give the other three solutions.
    let mut quaternions = [base; 4];
    for (slot, generator) in [Quaternion::i(), Quaternion::j(), Quaternion::k()]
        .into_iter()
        .enumerate()
    {
        quaternions[slot] = (second * (generator * second_inverse)) * base;
    }

    let mut rotations = [RotationAxis::new(0.0, 0.0, 0.0); 4];
    for (slot, quaternion) in quaternions.iter().enumerate() {
        rotations[slot] = RotationAxis::from_quaternion(quaternion)
            .ok_or(FocalMechanismError::DegenerateQuaternion)?
            .to_minimum();
    }

    // Smallest turn first. `total_cmp` rather than `partial_cmp().unwrap()`:
    // these angles cannot be NaN, and saying so with a total order costs
    // nothing and keeps a panic off the page.
    rotations.sort_by(|a, b| a.angle_degrees.abs().total_cmp(&b.angle_degrees.abs()));

    Ok(rotations)
}

/// The Kagan angle: how far apart two focal mechanisms are, in degrees.
///
/// The smallest of the four rotations, unsigned. Between 0 and 120 degrees for
/// any pair of double couples -- 120 and not 180, which follows from the
/// fourfold symmetry and is the cheapest sanity check there is on an
/// implementation of this.
pub fn kagan_angle(from: &PTBAxes, to: &PTBAxes) -> Result<f64, FocalMechanismError> {
    Ok(focal_mechanism_rotations(from, to)?[0].angle_degrees.abs())
}

/// Turn a focal mechanism about an axis.
///
/// The operation `focal_mechanism_rotations` inverts: applying any of the four
/// rotations it returns to `from` gives back `to`.
pub fn rotate_focal_mechanism(
    mechanism: &PTBAxes,
    rotation: &RotationAxis,
) -> Result<PTBAxes, FocalMechanismError> {

    PTBAxes::from_quaternion(&(rotation.to_quaternion() * mechanism.to_quaternion()))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The pair of focal mechanisms Kagan (1991) works his example on, as
    /// `(P trend, P plunge)` and `(T trend, T plunge)`.
    const KAGAN_FIRST: ((f64, f64), (f64, f64)) = ((232.0, 41.0), (120.0, 24.0));
    const KAGAN_SECOND: ((f64, f64), (f64, f64)) = ((51.0, 17.0), (295.0, 55.0));

    /// Kagan's four published solutions, as `(azimuth, colatitude from the
    /// bottom, rotation angle)` in degrees -- his own parameterisation, left in
    /// it rather than converted here, so that what is written down is what the
    /// paper prints.
    const KAGAN_SOLUTIONS: [(f64, f64, f64); 4] = [
        (24.8, 101.2, 102.8),
        (257.5, 79.7, 104.3),
        (144.8, 105.2, 124.1),
        (96.8, 16.7, 165.9),
    ];

    fn mechanism(spec: ((f64, f64), (f64, f64))) -> PTBAxes {
        PTBAxes::new(
            GeologicalAxis::new(spec.0.0, spec.0.1),
            GeologicalAxis::new(spec.1.0, spec.1.1),
        )
        .expect("Kagan's own axes are orthogonal to within a degree")
    }

    #[test]
    fn the_four_solutions_are_the_ones_kagan_published() {
        // The test this whole lineage exists to pass, against numbers printed
        // in 1991 rather than against any code.
        let solutions = focal_mechanism_rotations(&mechanism(KAGAN_FIRST), &mechanism(KAGAN_SECOND))
            .expect("two proper mechanisms");

        for (found, &(azimuth, colatitude, angle)) in solutions.iter().zip(&KAGAN_SOLUTIONS) {

            // Kagan gives colatitude measured from the downward vertical, so
            // the plunge is its complement.
            let plunge = 90.0 - colatitude;

            // An axis and its opposite end describe the same rotation, so
            // either representation is accepted.
            let deviation = [(azimuth, plunge, angle), (azimuth + 180.0, -plunge, -angle)]
                .iter()
                .map(|&(trend, plunge, angle)| {
                    let trend_apart = ((found.axis.trend - trend + 180.0).rem_euclid(360.0) - 180.0).abs();
                    trend_apart
                        .max((found.axis.plunge - plunge).abs())
                        .max((found.angle_degrees - angle).abs())
                })
                .fold(f64::INFINITY, f64::min);

            // The published values are quoted to a tenth of a degree, and the
            // worst component deviation over the four is 0.20 degrees.
            assert!(
                deviation < 0.25,
                "solution {found:?} departs from Kagan's ({azimuth}, {plunge}, {angle}) by \
                 {deviation:.4} degrees"
            );
        }
    }

    #[test]
    fn every_solution_really_carries_one_mechanism_onto_the_other() {
        // Stronger than matching the published table, and independent of it:
        // all four are rotations that do the job, which is what makes them
        // four solutions rather than one solution and three artefacts.
        let first = mechanism(KAGAN_FIRST);
        let second = mechanism(KAGAN_SECOND);

        for rotation in focal_mechanism_rotations(&first, &second).expect("two mechanisms") {

            let turned = rotate_focal_mechanism(&first, &rotation).expect("a mechanism");

            assert!(
                turned.almost_equal(&second, 1e-6),
                "rotating by {rotation:?} gave P {:?} T {:?}, not P {:?} T {:?}",
                turned.p_axis(),
                turned.t_axis(),
                second.p_axis(),
                second.t_axis()
            );
        }
    }

    #[test]
    fn the_kagan_angle_is_the_smallest_of_the_four() {
        let first = mechanism(KAGAN_FIRST);
        let second = mechanism(KAGAN_SECOND);

        let solutions = focal_mechanism_rotations(&first, &second).unwrap();
        let angle = kagan_angle(&first, &second).unwrap();

        assert!((angle - solutions[0].angle_degrees.abs()).abs() < 1e-12);

        // Sorted, and the published minimum is 102.8 degrees.
        for pair in solutions.windows(2) {
            assert!(pair[0].angle_degrees.abs() <= pair[1].angle_degrees.abs());
        }
        assert!((angle - 102.8).abs() < 0.25, "Kagan angle {angle}");
    }

    #[test]
    fn no_pair_of_double_couples_is_more_than_120_degrees_apart() {
        // The cheapest check there is on an implementation of this, and one
        // that a wrong set of symmetry generators fails immediately: the
        // fourfold symmetry bounds the minimum rotation at 120 degrees, not at
        // the 180 a general pair of frames would allow.
        let mut seed = 424242u64;
        let mut next = || {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (seed >> 33) as f64 / (1u64 << 31) as f64
        };

        let mut worst: f64 = 0.0;

        for _ in 0..2000 {

            // A random plane and a random rake give a mechanism that is
            // certainly a double couple, which random P and T axes would not.
            let build = |trend: f64, plunge: f64, roll: f64| {
                let t = GeologicalAxis::new(trend, plunge);
                // Any axis at a right angle to T serves as P, and rolling
                // about T sweeps them all.
                let t_versor = t.as_versor().as_vector();
                let helper = if t_versor.coords[2].abs() < 0.9 {
                    crate::algebra::vector::Vector::new([0.0, 0.0, 1.0])
                } else {
                    crate::algebra::vector::Vector::new([1.0, 0.0, 0.0])
                };
                let u = t_versor.cross(&helper).normalize().unwrap();
                let v = t_versor.cross(&u);
                let angle = roll.to_radians();
                let p = u * angle.cos() + v * angle.sin();

                PTBAxes::from_vectors(p, t_versor).expect("orthogonal by construction")
            };

            let a = build(next() * 360.0, next() * 180.0 - 90.0, next() * 360.0);
            let b = build(next() * 360.0, next() * 180.0 - 90.0, next() * 360.0);

            let angle = kagan_angle(&a, &b).expect("two mechanisms");

            assert!(angle.is_finite(), "not a number");
            worst = worst.max(angle);
        }

        assert!(worst <= 120.0 + 1e-9, "worst Kagan angle {worst} exceeds 120 degrees");
        // And the bound is approached, or the test would pass on a function
        // that always returned zero.
        assert!(worst > 110.0, "worst Kagan angle {worst} never approached the bound");
    }

    #[test]
    fn a_mechanism_is_no_rotation_from_itself() {
        // The case the Fortran and geogst both guard with a special branch and
        // this one does not need. The base rotation is nothing, and the other
        // three are the half turns the double couple is invariant under, which
        // is a truthful answer rather than a degenerate one.
        let ptb = mechanism(KAGAN_FIRST);

        let solutions = focal_mechanism_rotations(&ptb, &ptb).expect("a mechanism");

        assert!(solutions[0].angle_degrees.abs() < 1e-6, "{:?}", solutions[0]);

        for solution in &solutions[1..] {
            assert!(
                (solution.angle_degrees.abs() - 180.0).abs() < 1e-6,
                "expected a half turn, got {solution:?}"
            );
        }

        // And every one of them leaves the mechanism where it was, half turns
        // about P, T and B being exactly the double couple's own symmetries.
        for solution in solutions {
            let turned = rotate_focal_mechanism(&ptb, &solution).expect("a mechanism");
            assert!(turned.almost_equal(&ptb, 1e-6), "{solution:?} moved it");
        }
    }

    #[test]
    fn the_rotation_is_reported_the_shorter_way_round() {
        let long_way = RotationAxis::new(90.0, 30.0, 300.0);
        let short = long_way.to_minimum();

        assert!(short.angle_degrees.abs() <= 180.0, "{short:?}");
        // 300 degrees one way is 60 the other, about the opposite end.
        assert!((short.angle_degrees - 60.0).abs() < 1e-12, "{short:?}");
        assert!((short.axis.trend - 270.0).abs() < 1e-12, "{short:?}");
        assert!((short.axis.plunge + 30.0).abs() < 1e-12, "{short:?}");

        // Already short: left alone.
        let already = RotationAxis::new(90.0, 30.0, 120.0);
        assert_eq!(already.to_minimum(), already);
    }

    #[test]
    fn a_specular_rotation_is_the_same_rotation() {
        // Not merely a bookkeeping identity: turning a mechanism by either
        // description has to land it in the same place.
        let ptb = mechanism(KAGAN_FIRST);
        let rotation = RotationAxis::new(35.0, 20.0, 250.0);

        let one = rotate_focal_mechanism(&ptb, &rotation).expect("a mechanism");
        let other = rotate_focal_mechanism(&ptb, &rotation.specular()).expect("a mechanism");

        assert!(one.almost_equal(&other, 1e-6));
    }

    #[test]
    fn a_null_quaternion_has_an_angle_and_no_axis() {
        let nothing = RotationAxis::from_quaternion(&Quaternion::identity()).expect("a rotation");

        assert_eq!(nothing.angle_degrees, 0.0);

        // Too small to be an orientation at all is a different case, and is
        // declined rather than answered.
        assert!(RotationAxis::from_quaternion(&Quaternion::from_scalar_first(
            [0.0, 0.0, 0.0, 0.0]
        ))
        .is_none());
    }

    #[test]
    fn a_rotation_survives_the_round_trip_through_a_quaternion() {
        for angle in [-170.0, -90.0, -1.0, 0.5, 45.0, 179.0] {

            let rotation = RotationAxis::new(217.0, 33.0, angle);
            let back = RotationAxis::from_quaternion(&rotation.to_quaternion())
                .expect("a rotation")
                .to_minimum();

            // The axis may come back reversed with the angle negated, which is
            // the same rotation; compare through what it does rather than
            // through how it is written.
            let ptb = mechanism(KAGAN_FIRST);
            let one = rotate_focal_mechanism(&ptb, &rotation).expect("a mechanism");
            let other = rotate_focal_mechanism(&ptb, &back).expect("a mechanism");

            assert!(one.almost_equal(&other, 1e-6), "{angle}: {rotation:?} against {back:?}");
        }
    }
}
