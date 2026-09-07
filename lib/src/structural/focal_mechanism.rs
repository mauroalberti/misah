//! A focal mechanism as its kinematic axes, and the rotation matrix behind
//! them.
//!
//! A double-couple focal mechanism is usually quoted as two nodal planes, and
//! the two are interchangeable: nothing in the seismology says which of them
//! slipped. What *is* unambiguous is the orthogonal triad the pair defines --
//! **P** for the axis of maximum shortening, **T** for maximum extension, and
//! **B**, the null axis, perpendicular to both. That triad is what this module
//! carries, and it is the form Kagan (1991) works in, because a triad is an
//! orthonormal frame and the rotation between two frames is a quaternion.
//!
//! P and T are not stress axes. They are kinematic, fixed at 45 degrees to the
//! fault plane by construction, and coincide with the principal stresses only
//! under an assumption -- Anderson's -- that real rock keeps only roughly.
//! Where this crate means stress it says stress: `structural::stress` and
//! `structural::inversion` are that, and they take fault-slip data rather than
//! focal mechanisms. Conflating the two is the standard way of reporting a
//! stress field that is really a strain field.
//!
//! ## Fixing T and recomputing P
//!
//! Two measured axes are never exactly orthogonal, and every formula
//! downstream assumes they are. So the constructor takes both, refuses the pair
//! if it departs from orthogonality by more than `PTB_ORTHOGONALITY_TOLERANCE`,
//! and then makes the triad exact by keeping T as given and recomputing P from
//! B and T.
//!
//! Keeping T rather than splitting the difference is a decision copied from
//! geogst and from the Fortran before it, and it is worth naming as a decision:
//! the result is not symmetric in P and T, so a triad built one way and one
//! built with the roles exchanged do not agree to the last bit. They agree to
//! well within the tolerance that let them in.

use thiserror::Error;

use crate::algebra::quaternion::Quaternion;
use crate::algebra::vector::Vector;
use crate::algebra::versor::Versor3D;
use crate::structural::fault::FaultPlane;
use crate::structural::geol_axis::GeologicalAxis;

/// How far a measured P/T pair may depart from orthogonality, in degrees.
///
/// One degree, which is what a pair of compass readings of the same structure
/// is worth, and the same figure `fault::SLICKENLINE_COPLANARITY_TOLERANCE`
/// uses for the same reason.
pub const PTB_ORTHOGONALITY_TOLERANCE: f64 = 1.0;

#[derive(Debug, Error, PartialEq)]
pub enum FocalMechanismError {

    #[error(
        "P and T axes depart from orthogonality by {departure:.3} degrees, more than the \
         {tolerance:.3} allowed"
    )]
    NotOrthogonal { departure: f64, tolerance: f64 },

    #[error("a P or T axis was given as a zero-length vector, which has no direction")]
    ZeroLengthAxis,

    #[error("the fault carries no slickenline at index {index}")]
    NoSuchSlickenline { index: usize },

    #[error(
        "the slickenline has no determined sense of movement, so P and T cannot be told apart"
    )]
    UnknownSlipSense,

    #[error("the quaternion is too small to describe an orientation")]
    DegenerateQuaternion,
}

/// The same axis, pointing into the ground.
///
/// P and T are **axes**: undirected, so a vector and its negation describe the
/// same one. Without a rule for choosing between them the triad would depend on
/// how the axis happened to be handed over rather than on the mechanism, and
/// the two constructors here would disagree about the same fault --
/// `from_fault_slickenline` computes `normal - slip`, whose sign follows from
/// the arithmetic, while `new` takes a trend and plunge a geologist wrote down,
/// which is downward by convention. Two triads differing by that sign are a
/// half turn apart, and a rotation between them would come back as 180 degrees
/// where the answer is zero.
///
/// Downward is the geological convention and geogst's. A horizontal axis is
/// left alone: both its ends are equally horizontal, and either serves.
///
/// "Horizontal" has to be a tolerance and not a test against zero, which is
/// the second thing this function got wrong. A fault dipping 45 degrees with
/// down-dip slip has `normal + slip` exactly horizontal in exact arithmetic
/// and about 6e-15 above horizontal in floating point, so a strict `> 0` puts
/// the sign of that rounding in charge of whether the axis turns over --
/// reintroducing, at the last bit, the very dependence on accident this
/// function exists to remove. Anything within `HORIZONTAL_TOLERANCE` of
/// horizontal is horizontal, and is left as it came.
fn downward(v: Vector<3>) -> Vector<3> {

    // Index 2 is Up, so a positive third component points at the sky.
    if v.coords[2] > HORIZONTAL_TOLERANCE { -v } else { v }
}

/// How far the vertical component of a unit vector may stray from zero and
/// still count as horizontal, in `downward`.
///
/// A thousandth of a millionth of a degree of plunge. Far below any
/// measurement, and far above the rounding of a normalized cross product,
/// which is where the ambiguity actually comes from.
const HORIZONTAL_TOLERANCE: f64 = 1.0e-12;

/// The P, T and B kinematic axes of a focal mechanism.
///
/// Stored as two versors, T and P, with B derived as their cross product --
/// there being no third independent quantity, and a stored B being a thing that
/// can fall out of step with the other two.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PTBAxes {
    t: Versor3D,
    p: Versor3D,
}

impl PTBAxes {

    /// From the two axes as trend and plunge.
    pub fn new(
        p_axis: GeologicalAxis,
        t_axis: GeologicalAxis,
    ) -> Result<Self, FocalMechanismError> {
        Self::from_vectors(p_axis.as_versor().as_vector(), t_axis.as_versor().as_vector())
    }

    /// From the two axes as vectors, which need not be unit vectors but must be
    /// sub-orthogonal.
    pub fn from_vectors(
        p_vector: Vector<3>,
        t_vector: Vector<3>,
    ) -> Result<Self, FocalMechanismError> {

        let p = downward(p_vector.normalize().ok_or(FocalMechanismError::ZeroLengthAxis)?);
        let t = downward(t_vector.normalize().ok_or(FocalMechanismError::ZeroLengthAxis)?);

        // The departure from a right angle, taken unsigned, so a pair 180
        // degrees apart is refused as loudly as one 0 degrees apart. Neither
        // is a triad. Unsigned is also what makes the check indifferent to the
        // canonicalization just applied, which can flip the sign of the dot
        // product without changing the geometry.
        let departure = (p.dot(&t).clamp(-1.0, 1.0).acos().to_degrees() - 90.0).abs();

        if departure > PTB_ORTHOGONALITY_TOLERANCE {
            return Err(FocalMechanismError::NotOrthogonal {
                departure,
                tolerance: PTB_ORTHOGONALITY_TOLERANCE,
            });
        }

        // T is kept and P is rebuilt, which is what makes the triad exactly
        // orthonormal. See the module note on why it is T that is kept.
        let b = t.cross(&p);
        let p = b.cross(&t);

        Ok(Self {
            t: Versor3D::try_from(t).map_err(|_| FocalMechanismError::ZeroLengthAxis)?,
            p: Versor3D::try_from(p).map_err(|_| FocalMechanismError::ZeroLengthAxis)?,
        })
    }

    /// The P and T axes a fault plane and its observed slip imply.
    ///
    /// `P = normal - slip` and `T = normal + slip`, which places both at 45
    /// degrees to the fault plane. That 45 degrees is an assumption of the
    /// construction and not a measurement -- it is what makes these kinematic
    /// axes rather than stress axes.
    ///
    /// The **sense of movement must be known**. Without it the slip is a line
    /// rather than a vector, and reversing it exchanges P with T: the same
    /// fault would come back as shortening where it was extension. A
    /// slickenline of undetermined sense is refused here rather than given a
    /// direction -- a different answer from `inversion`'s, where the same
    /// missing information is handled by taking the misfit modulo 180 degrees.
    /// There is no modulo available here.
    pub fn from_fault_slickenline(
        fault: &FaultPlane,
        index: usize,
    ) -> Result<Self, FocalMechanismError> {

        let slickenline = fault
            .slickenline(index)
            .ok_or(FocalMechanismError::NoSuchSlickenline { index })?;

        if !slickenline.sense_is_known() {
            return Err(FocalMechanismError::UnknownSlipSense);
        }

        let normal = fault.plane().normal_vector();
        let slip = slickenline.direction().as_vector();

        Self::from_vectors(normal - slip, normal + slip)
    }

    /// The triad a quaternion describes.
    ///
    /// Eq. 10 in Kagan (1991), which is the first two columns of the rotation
    /// matrix -- so it is read off `Quaternion::to_rotation_matrix` rather than
    /// written out a second time. Two copies of the same trigonometry are two
    /// things to keep in step, and this way `to_quaternion` and this one are
    /// inverse by construction rather than by inspection.
    pub fn from_quaternion(quaternion: &Quaternion) -> Result<Self, FocalMechanismError> {

        let matrix = quaternion
            .to_rotation_matrix()
            .ok_or(FocalMechanismError::DegenerateQuaternion)?;

        let column = |j: usize| Vector::new([matrix[0][j], matrix[1][j], matrix[2][j]]);

        Self::from_vectors(column(1), column(0))
    }

    pub fn t_versor(&self) -> Versor3D {
        self.t
    }

    pub fn p_versor(&self) -> Versor3D {
        self.p
    }

    /// The null axis, `T x P`. Derived rather than stored.
    pub fn b_versor(&self) -> Versor3D {
        Versor3D::try_from(self.t.as_vector().cross(&self.p.as_vector()))
            .expect("the cross product of two orthonormal versors is a unit vector")
    }

    pub fn t_axis(&self) -> GeologicalAxis {
        GeologicalAxis::from_versor(&self.t)
    }

    pub fn p_axis(&self) -> GeologicalAxis {
        GeologicalAxis::from_versor(&self.p)
    }

    pub fn b_axis(&self) -> GeologicalAxis {
        GeologicalAxis::from_versor(&self.b_versor())
    }

    /// The rotation matrix of this triad, indexed `[row][column]`, with T, P
    /// and B as its **columns**.
    ///
    /// Eq. 3 in Kagan (2007) and eq. 31 in Kagan (2008). Columns and not rows:
    /// the matrix carries the reference frame into this mechanism's frame, and
    /// transposing it would silently give the rotation the other way round.
    pub fn to_matrix(&self) -> [[f64; 3]; 3] {

        let t = self.t.coords();
        let p = self.p.coords();
        let b = self.b_versor();
        let b = b.coords();

        [
            [t[0], p[0], b[0]],
            [t[1], p[1], b[1]],
            [t[2], p[2], b[2]],
        ]
    }

    /// The quaternion of this triad.
    pub fn to_quaternion(&self) -> Quaternion {
        Quaternion::from_rotation_matrix(&self.to_matrix())
    }

    /// Whether two mechanisms have their P and T axes within
    /// `tolerance_degrees` of each other.
    ///
    /// Axes and not directions: an axis and its opposite end are the same axis,
    /// so the comparison is on the absolute dot product.
    pub fn almost_equal(&self, other: &Self, tolerance_degrees: f64) -> bool {

        let apart =
            |a: &Versor3D, b: &Versor3D| a.dot(b).abs().clamp(-1.0, 1.0).acos().to_degrees();

        apart(&self.t, &other.t) <= tolerance_degrees
            && apart(&self.p, &other.p) <= tolerance_degrees
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::orientation::direction::Direction3D;
    use crate::structural::geol_plane::GeologicalPlane;
    use crate::structural::slickenline::{SlipSense, Slickenline};

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9
    }

    #[test]
    fn a_normal_fault_puts_p_up_and_t_along_the_slip() {
        // A plane dipping 45 to the east, slipping down-dip: extension
        // east-west, so T is horizontal and east, P is vertical.
        //
        // The same case geogst's own doctest uses, so the two can be compared
        // directly rather than by reasoning about sign conventions.
        let plane = GeologicalPlane::new(90.0, 45.0);
        let slip = Slickenline::from_axis(GeologicalAxis::new(90.0, 45.0), Some(SlipSense::Down));
        let fault = FaultPlane::new(plane, vec![slip]).expect("the slip lies in the plane");

        let ptb = PTBAxes::from_fault_slickenline(&fault, 0).expect("a triad");

        let t = ptb.t_axis();
        let p = ptb.p_axis();

        assert!(close(t.trend, 90.0) && close(t.plunge, 0.0), "T {t:?}");
        // Plunge 90 and not -90: `normal - slip` comes out pointing at the
        // sky here, and the canonicalization turns it over. Reading the axis
        // as vertical-up would be the same geometry and the wrong convention.
        assert!(close(p.plunge, 90.0), "P {p:?}");

        // And the versors, against geogst's printed values for this fault.
        assert!(close(ptb.t_versor().coords()[0], 1.0));
        assert!(close(ptb.p_versor().coords()[2], -1.0));
        assert!(close(ptb.b_versor().coords()[1], 1.0));
    }

    #[test]
    fn the_triad_is_orthonormal_whatever_went_in() {
        let ptb = PTBAxes::new(
            GeologicalAxis::new(232.0, 41.0),
            GeologicalAxis::new(120.0, 24.0),
        )
        .expect("within tolerance");

        let t = ptb.t_versor();
        let p = ptb.p_versor();
        let b = ptb.b_versor();

        assert!(t.dot(&p).abs() < 1e-15, "T.P = {}", t.dot(&p));
        assert!(t.dot(&b).abs() < 1e-15, "T.B = {}", t.dot(&b));
        assert!(p.dot(&b).abs() < 1e-15, "P.B = {}", p.dot(&b));

        for versor in [t, p, b] {
            let norm = versor.as_vector().norm();
            assert!((norm - 1.0).abs() < 1e-15, "norm {norm}");
        }

        // T came through untouched, which is the documented asymmetry.
        let given = GeologicalAxis::new(120.0, 24.0).as_versor();
        for i in 0..3 {
            assert!(close(t.coords()[i], given.coords()[i]));
        }
    }

    #[test]
    fn axes_that_are_not_a_right_angle_apart_are_refused() {
        let err = PTBAxes::new(
            GeologicalAxis::new(0.0, 0.0),
            GeologicalAxis::new(80.0, 0.0),
        )
        .unwrap_err();

        assert!(matches!(err, FocalMechanismError::NotOrthogonal { .. }), "{err:?}");

        // A pair 180 degrees apart is refused just as loudly as one 0 degrees
        // apart: taking the deviation from 90 unsigned is what catches both.
        assert!(PTBAxes::new(
            GeologicalAxis::new(0.0, 0.0),
            GeologicalAxis::new(180.0, 0.0)
        )
        .is_err());

        // Half a degree is inside the tolerance, being what a compass reading
        // is worth.
        assert!(PTBAxes::new(
            GeologicalAxis::new(0.0, 0.0),
            GeologicalAxis::new(90.5, 0.0)
        )
        .is_ok());
    }

    #[test]
    fn a_slip_of_unknown_sense_is_refused_rather_than_guessed() {
        // Reversing an undetermined slip exchanges P with T, turning
        // shortening into extension. There is no defensible default.
        let plane = GeologicalPlane::new(90.0, 45.0);
        let slip = Slickenline::new(
            Direction3D::new(GeologicalAxis::new(90.0, 45.0).as_versor()),
            None,
        );
        let fault = FaultPlane::new(plane, vec![slip]).expect("the slip lies in the plane");

        assert_eq!(
            PTBAxes::from_fault_slickenline(&fault, 0).unwrap_err(),
            FocalMechanismError::UnknownSlipSense
        );

        assert_eq!(
            PTBAxes::from_fault_slickenline(&fault, 3).unwrap_err(),
            FocalMechanismError::NoSuchSlickenline { index: 3 }
        );
    }

    #[test]
    fn a_triad_survives_the_round_trip_through_a_quaternion() {
        // The two mechanisms Kagan (1991) works his example on.
        for (p, t) in [((232.0, 41.0), (120.0, 24.0)), ((51.0, 17.0), (295.0, 55.0))] {

            let ptb = PTBAxes::new(
                GeologicalAxis::new(p.0, p.1),
                GeologicalAxis::new(t.0, t.1),
            )
            .expect("within tolerance");

            let back = PTBAxes::from_quaternion(&ptb.to_quaternion()).expect("a triad");

            // Compared component by component and not through `almost_equal`,
            // which is the right tool at degree tolerances and the wrong one
            // here. It reads an angle out of `acos`, and `acos` near 1 loses
            // half the significant digits: one unit in the last place of the
            // dot product comes out as about 1e-6 of a degree. A round trip
            // that is exact to 1e-14 in the components therefore cannot be
            // asserted to better than 1e-6 through an angle. The components
            // are the well-conditioned comparison.
            for (a, b) in [
                (back.t_versor(), ptb.t_versor()),
                (back.p_versor(), ptb.p_versor()),
                (back.b_versor(), ptb.b_versor()),
            ] {
                for i in 0..3 {
                    assert!(
                        (a.coords()[i] - b.coords()[i]).abs() < 1e-13,
                        "P/T {p:?} {t:?}: component {i} came back as {} against {}",
                        a.coords()[i],
                        b.coords()[i]
                    );
                }
            }
        }
    }

    #[test]
    fn the_same_mechanism_gives_the_same_triad_however_it_was_specified() {
        // The property the downward canonicalization exists for. A fault gives
        // P as `normal - slip`, whose sign falls out of the arithmetic; read
        // that triad's own axes back and hand them to `new`, which is the
        // route a geologist's notebook takes, and the two must agree.
        //
        // Without the canonicalization these differ by the sign of P, which is
        // a half turn -- so a Kagan rotation between a mechanism and itself
        // would come back as 180 degrees rather than as nothing.
        let plane = GeologicalPlane::new(90.0, 45.0);
        let slip = Slickenline::from_axis(GeologicalAxis::new(90.0, 45.0), Some(SlipSense::Down));
        let fault = FaultPlane::new(plane, vec![slip]).expect("the slip lies in the plane");

        let from_fault = PTBAxes::from_fault_slickenline(&fault, 0).expect("a triad");
        let from_axes = PTBAxes::new(from_fault.p_axis(), from_fault.t_axis()).expect("a triad");

        for (a, b) in [
            (from_fault.t_versor(), from_axes.t_versor()),
            (from_fault.p_versor(), from_axes.p_versor()),
            (from_fault.b_versor(), from_axes.b_versor()),
        ] {
            for i in 0..3 {
                assert!(
                    (a.coords()[i] - b.coords()[i]).abs() < 1e-13,
                    "component {i}: {} against {}",
                    a.coords()[i],
                    b.coords()[i]
                );
            }
        }
    }

    #[test]
    fn the_quaternion_is_the_one_geogst_computes() {
        // Not a round trip: the actual components, against the verified Python
        // port. A sign convention that had drifted would survive a round trip
        // and fail here.
        let first = PTBAxes::new(
            GeologicalAxis::new(232.0, 41.0),
            GeologicalAxis::new(120.0, 24.0),
        )
        .unwrap()
        .to_quaternion();

        let second = PTBAxes::new(
            GeologicalAxis::new(51.0, 17.0),
            GeologicalAxis::new(295.0, 55.0),
        )
        .unwrap()
        .to_quaternion();

        for (got, expected) in first
            .components()
            .iter()
            .zip([-0.41567, 0.85017, -0.31120, -0.08706].iter())
        {
            assert!((got - expected).abs() < 1e-5, "{:?}", first.components());
        }

        for (got, expected) in second
            .components()
            .iter()
            .zip([0.38380, 0.30459, 0.80853, -0.32588].iter())
        {
            assert!((got - expected).abs() < 1e-5, "{:?}", second.components());
        }
    }

    #[test]
    fn the_matrix_carries_the_axes_in_its_columns() {
        let ptb = PTBAxes::new(
            GeologicalAxis::new(232.0, 41.0),
            GeologicalAxis::new(120.0, 24.0),
        )
        .unwrap();

        let matrix = ptb.to_matrix();
        let (t, p, b) = (ptb.t_versor(), ptb.p_versor(), ptb.b_versor());

        for (row, line) in matrix.iter().enumerate() {
            assert!(close(line[0], t.coords()[row]), "T in column 0");
            assert!(close(line[1], p.coords()[row]), "P in column 1");
            assert!(close(line[2], b.coords()[row]), "B in column 2");
        }

        // A rotation matrix, so the determinant is +1 and not -1: the triad is
        // right-handed, and a left-handed one would still be orthonormal.
        let d = matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
            - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
            + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);

        assert!(close(d, 1.0), "determinant {d}");
    }

    #[test]
    fn the_reference_triad_is_the_identity_matrix() {
        // P north-horizontal and T east-horizontal give the frame everything
        // else is measured against.
        let ptb = PTBAxes::new(
            GeologicalAxis::new(0.0, 0.0),
            GeologicalAxis::new(90.0, 0.0),
        )
        .unwrap();

        let matrix = ptb.to_matrix();
        for (row, line) in matrix.iter().enumerate() {
            for (column, value) in line.iter().enumerate() {
                let expected = if row == column { 1.0 } else { 0.0 };
                assert!(close(*value, expected), "[{row}][{column}] = {value}");
            }
        }
    }
}
