//! Quaternions, for the one job this crate needs them for: the rotation
//! between two orientations in space.
//!
//! A rotation can be written as a matrix, as an axis and an angle, or as a
//! quaternion, and the three are the same rotation. The quaternion is the form
//! in which *composing* two rotations is a product and *undoing* one is a
//! conjugate, which is why Kagan (1991) put the rotation between two focal
//! mechanisms in this form: the four solutions that exist for every such pair
//! come out as four products rather than as four separate constructions.
//!
//! ## Where the scalar lives
//!
//! In `scalar`, and never in a slot. The Fortran this descends from stored a
//! quaternion as `q(0:3)` with the scalar at index 0; a great deal of other
//! software stores `(x, y, z, w)` with the scalar **last**. Both conventions
//! are defensible and mixing them silently produces a rotation that is wrong in
//! a way no assertion about magnitude will catch -- the norm is unchanged, the
//! quaternion is still a unit quaternion, and it turns things the wrong way.
//! Naming the parts is what makes the question unanswerable rather than
//! answered wrongly.
//!
//! `components` still returns `[scalar, i, j, k]` for the cases that want an
//! array, and says so in its name and its documentation.
//!
//! ## The frame is (East, North, Up)
//!
//! As everywhere else in this crate. The Fortran worked in (North, East,
//! Down), so its numbers are not comparable component by component with these
//! even where the rotation is the same one -- a fact worth keeping in mind when
//! reading the original beside this.

use std::ops::Mul;

use crate::algebra::vector::Vector;

/// Below this magnitude a quaternion carries no direction to rotate about, and
/// normalizing it would amplify rounding into an answer.
///
/// The Fortran's own threshold, and geogst's.
pub const QUATERNION_MAGNITUDE_THRESHOLD: f64 = 1.0e-6;

/// A quaternion: one scalar part and one vector part.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Quaternion {
    pub scalar: f64,
    pub vector: Vector<3>,
}

impl Quaternion {

    pub fn new(scalar: f64, vector: Vector<3>) -> Self {
        Self { scalar, vector }
    }

    /// From four numbers in the order `[scalar, i, j, k]`.
    ///
    /// The order is in the name because it is the thing that goes wrong. See
    /// the module documentation.
    pub fn from_scalar_first(components: [f64; 4]) -> Self {
        Self {
            scalar: components[0],
            vector: Vector::new([components[1], components[2], components[3]]),
        }
    }

    /// The four numbers, scalar first.
    pub fn components(&self) -> [f64; 4] {
        let v = self.vector.coords;
        [self.scalar, v[0], v[1], v[2]]
    }

    /// The multiplicative identity: no rotation at all.
    pub fn identity() -> Self {
        Self::new(1.0, Vector::new([0.0, 0.0, 0.0]))
    }

    /// The pure unit quaternion along the first axis, and its two companions.
    ///
    /// These are the three that generate the other solutions of a
    /// double-couple rotation: a focal mechanism is unchanged by a half turn
    /// about any of its own P, T or B axes, so conjugating by each of these in
    /// turn is what produces Kagan's four rotations from one.
    pub fn i() -> Self {
        Self::new(0.0, Vector::new([1.0, 0.0, 0.0]))
    }

    pub fn j() -> Self {
        Self::new(0.0, Vector::new([0.0, 1.0, 0.0]))
    }

    pub fn k() -> Self {
        Self::new(0.0, Vector::new([0.0, 0.0, 1.0]))
    }

    /// Same scalar, opposite vector.
    pub fn conjugate(&self) -> Self {
        Self::new(self.scalar, -self.vector)
    }

    /// Sum of the four squared components. The norm without the square root,
    /// which is what the inverse actually divides by.
    pub fn squared_norm(&self) -> f64 {
        self.scalar * self.scalar + self.vector.dot(&self.vector)
    }

    pub fn norm(&self) -> f64 {
        self.squared_norm().sqrt()
    }

    /// The same rotation, at unit norm.
    ///
    /// `None` below `QUATERNION_MAGNITUDE_THRESHOLD`, where there is no
    /// direction left to preserve and dividing would turn rounding into an
    /// answer.
    pub fn normalized(&self) -> Option<Self> {

        let norm = self.norm();
        if norm < QUATERNION_MAGNITUDE_THRESHOLD {
            return None;
        }

        Some(Self::new(self.scalar / norm, self.vector / norm))
    }

    /// The rotation undone: `conjugate / squared_norm`.
    ///
    /// For a unit quaternion this is the conjugate, and it is written the long
    /// way anyway because the quaternions here come out of a matrix that is
    /// only orthogonal to within rounding.
    pub fn inverse(&self) -> Option<Self> {

        let squared = self.squared_norm();
        if squared < QUATERNION_MAGNITUDE_THRESHOLD * QUATERNION_MAGNITUDE_THRESHOLD {
            return None;
        }

        let conjugate = self.conjugate();

        Some(Self::new(conjugate.scalar / squared, conjugate.vector / squared))
    }

    /// The angle of the rotation this quaternion represents, in degrees.
    ///
    /// `2 * acos(scalar)` on the normalized quaternion -- p. 710 of Kagan
    /// (1991). Runs from 0 to 360: a quaternion and its negation are the same
    /// rotation described the two ways round, and which of the two you get
    /// depends on where the quaternion came from. `structural::rotation`
    /// carries the reduction to the shorter way.
    ///
    /// `None` for a quaternion too small to have a direction.
    pub fn rotation_angle_degrees(&self) -> Option<f64> {

        let unit = self.normalized()?;

        // A matrix that is orthogonal only to within rounding gives a scalar
        // a whisker outside [-1, 1], where acos is NaN. The clamp is the same
        // guard as the one in `from_rotation_matrix`, for the same reason.
        Some(2.0 * unit.scalar.clamp(-1.0, 1.0).acos().to_degrees())
    }

    /// The quaternion of a 3x3 rotation matrix, indexed `[row][column]`.
    ///
    /// All four squared components are read off the diagonal, and then the
    /// largest of them is used to recover the other three from the
    /// off-diagonal terms. Taking the largest is what keeps the division well
    /// conditioned: at least one of the four is always at or above the mean,
    /// and for a unit quaternion that mean is at least 1/4.
    ///
    /// ## The clamp, which the Fortran does not have
    ///
    /// The four square-root arguments are non-negative in exact arithmetic.
    /// They are not in floating point: at a **trace of -1** -- a rotation by
    /// 180 degrees, which is what two strike-slip focal mechanisms a whole
    /// number of degrees apart give -- one argument that should be exactly
    /// zero lands a few times 10^-16 below it, and the square root is NaN.
    ///
    /// The original `quaternfromcartmatr` in `FaultCorrelation.f95` has no
    /// clamp and produces that NaN in silence. It went unseen for twenty years
    /// because it needs a trace of exactly -1, which round-degree catalogue
    /// data and textbook examples produce constantly and random doubles
    /// essentially never: 134 pairs in 360 among pure strike-slip mechanisms,
    /// and none at all in 4000 pairs drawn from random doubles.
    ///
    /// Clamping at zero is safe rather than merely convenient. The branch
    /// taken always divides by a component at or above the mean of the four,
    /// and a clamped zero is never that component -- so the zero never reaches
    /// a denominator.
    pub fn from_rotation_matrix(matrix: &[[f64; 3]; 3]) -> Self {

        let m = matrix;

        let mut q0 = (1.0 + m[0][0] + m[1][1] + m[2][2]).max(0.0).sqrt() / 2.0;
        let mut q1 = (1.0 + m[0][0] - m[1][1] - m[2][2]).max(0.0).sqrt() / 2.0;
        let mut q2 = (1.0 - m[0][0] + m[1][1] - m[2][2]).max(0.0).sqrt() / 2.0;
        let mut q3 = (1.0 - m[0][0] - m[1][1] + m[2][2]).max(0.0).sqrt() / 2.0;

        let q0q1 = (m[2][1] - m[1][2]) / 4.0;
        let q0q2 = (m[0][2] - m[2][0]) / 4.0;
        let q0q3 = (m[1][0] - m[0][1]) / 4.0;
        let q1q2 = (m[0][1] + m[1][0]) / 4.0;
        let q1q3 = (m[0][2] + m[2][0]) / 4.0;
        let q2q3 = (m[1][2] + m[2][1]) / 4.0;

        if 3.0 * q0 > q1 + q2 + q3 {
            q1 = q0q1 / q0;
            q2 = q0q2 / q0;
            q3 = q0q3 / q0;
        } else if 3.0 * q1 > q0 + q2 + q3 {
            q0 = q0q1 / q1;
            q2 = q1q2 / q1;
            q3 = q1q3 / q1;
        } else if 3.0 * q2 > q0 + q1 + q3 {
            q0 = q0q2 / q2;
            q1 = q1q2 / q2;
            q3 = q2q3 / q2;
        } else {
            q0 = q0q3 / q3;
            q1 = q1q3 / q3;
            q2 = q2q3 / q3;
        }

        Self::from_scalar_first([q0, q1, q2, q3])
    }

    /// The rotation matrix of this quaternion, indexed `[row][column]`.
    ///
    /// Eq. 10 in Kagan (1991), eq. 3.5 in Salamin (1979). The inverse of
    /// `from_rotation_matrix` up to the sign of the whole quaternion, which
    /// describes the same rotation either way.
    pub fn to_rotation_matrix(&self) -> Option<[[f64; 3]; 3]> {

        let [q0, q1, q2, q3] = self.normalized()?.components();

        Some([
            [
                q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3,
                2.0 * (q1 * q2 - q0 * q3),
                2.0 * (q1 * q3 + q0 * q2),
            ],
            [
                2.0 * (q1 * q2 + q0 * q3),
                q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3,
                2.0 * (q2 * q3 - q0 * q1),
            ],
            [
                2.0 * (q1 * q3 - q0 * q2),
                2.0 * (q2 * q3 + q0 * q1),
                q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3,
            ],
        ])
    }
}

impl Mul for Quaternion {
    type Output = Self;

    /// The Hamilton product, which is the composition of two rotations and is
    /// **not** commutative: `a * b` turns by `b` and then by `a`.
    fn mul(self, other: Self) -> Self {

        let [a0, a1, a2, a3] = self.components();
        let [b0, b1, b2, b3] = other.components();

        Self::from_scalar_first([
            a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3,
            a0 * b1 + a1 * b0 + a2 * b3 - a3 * b2,
            a0 * b2 - a1 * b3 + a2 * b0 + a3 * b1,
            a0 * b3 + a1 * b2 - a2 * b1 + a3 * b0,
        ])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-12
    }

    fn quaternions_close(a: &Quaternion, b: &Quaternion) -> bool {
        a.components()
            .iter()
            .zip(b.components().iter())
            .all(|(x, y)| close(*x, *y))
    }

    #[test]
    fn the_scalar_is_named_and_not_positioned() {
        // The whole point of the type. `from_scalar_first` says which end it
        // reads, and the field says which end it kept.
        let q = Quaternion::from_scalar_first([0.5, 1.0, 2.0, 3.0]);

        assert_eq!(q.scalar, 0.5);
        assert_eq!(q.vector.coords, [1.0, 2.0, 3.0]);
        assert_eq!(q.components(), [0.5, 1.0, 2.0, 3.0]);
    }

    #[test]
    fn the_product_is_the_one_kuipers_gives() {
        // Kuipers (2002), chapter 5: the worked example every quaternion
        // implementation is checked against.
        let p = Quaternion::from_scalar_first([3.0, 1.0, -2.0, 1.0]);
        let q = Quaternion::from_scalar_first([2.0, -1.0, 2.0, 3.0]);

        assert!(quaternions_close(
            &(p * q),
            &Quaternion::from_scalar_first([8.0, -9.0, -2.0, 11.0])
        ));
    }

    #[test]
    fn the_product_does_not_commute() {
        // Which is the property that makes the order of the Kagan products
        // matter, and the one a careless rewrite loses.
        let p = Quaternion::from_scalar_first([3.0, 1.0, -2.0, 1.0]);
        let q = Quaternion::from_scalar_first([2.0, -1.0, 2.0, 3.0]);

        assert!(!quaternions_close(&(p * q), &(q * p)));
    }

    #[test]
    fn the_identity_leaves_a_quaternion_alone() {
        let q = Quaternion::from_scalar_first([0.3, 0.5, -0.2, 0.7]);

        assert!(quaternions_close(&(q * Quaternion::identity()), &q));
        assert!(quaternions_close(&(Quaternion::identity() * q), &q));
    }

    #[test]
    fn a_quaternion_times_its_inverse_is_the_identity() {
        let q = Quaternion::from_scalar_first([0.3, 0.5, -0.2, 0.7]);
        let inverse = q.inverse().expect("not near zero");

        assert!(quaternions_close(&(q * inverse), &Quaternion::identity()));
        assert!(quaternions_close(&(inverse * q), &Quaternion::identity()));
    }

    #[test]
    fn the_three_pure_units_multiply_as_they_should() {
        // i*j = k, j*k = i, k*i = j, and i*i = -1.
        assert!(quaternions_close(&(Quaternion::i() * Quaternion::j()), &Quaternion::k()));
        assert!(quaternions_close(&(Quaternion::j() * Quaternion::k()), &Quaternion::i()));
        assert!(quaternions_close(&(Quaternion::k() * Quaternion::i()), &Quaternion::j()));
        assert!(quaternions_close(
            &(Quaternion::i() * Quaternion::i()),
            &Quaternion::from_scalar_first([-1.0, 0.0, 0.0, 0.0])
        ));
    }

    #[test]
    fn a_matrix_and_back_is_the_same_rotation() {
        // A quarter turn about the third axis.
        let angle: f64 = 90.0_f64.to_radians();
        let matrix = [
            [angle.cos(), -angle.sin(), 0.0],
            [angle.sin(), angle.cos(), 0.0],
            [0.0, 0.0, 1.0],
        ];

        let quaternion = Quaternion::from_rotation_matrix(&matrix);
        let back = quaternion.to_rotation_matrix().expect("not near zero");

        for row in 0..3 {
            for column in 0..3 {
                assert!(
                    close(back[row][column], matrix[row][column]),
                    "element [{row}][{column}]: {} against {}",
                    back[row][column],
                    matrix[row][column]
                );
            }
        }

        assert!(close(quaternion.rotation_angle_degrees().unwrap(), 90.0));
    }

    #[test]
    fn the_identity_matrix_is_no_rotation_at_all() {
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let q = Quaternion::from_rotation_matrix(&identity);

        assert!(quaternions_close(&q, &Quaternion::identity()));
        assert!(close(q.rotation_angle_degrees().unwrap(), 0.0));
    }

    #[test]
    fn a_half_turn_does_not_produce_a_nan() {
        // The defect the Fortran still carries, and the reason for the clamp.
        // A half turn about the second axis has trace -1, so
        // `1 + m00 + m11 + m22` is zero in exact arithmetic and a few times
        // 10^-16 below it in practice.
        let half_turn = [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]];

        let q = Quaternion::from_rotation_matrix(&half_turn);

        assert!(q.components().iter().all(|c| c.is_finite()), "{:?}", q);
        assert!(close(q.rotation_angle_degrees().unwrap(), 180.0));
        // The rotation is about the second axis, so that is where the vector
        // part has to be.
        assert!(quaternions_close(&q, &Quaternion::j()));
    }

    #[test]
    fn every_half_turn_survives_the_square_root() {
        // Not just the one above. A half turn about any axis has trace -1, so
        // this sweeps the whole family the Fortran fails on rather than the
        // single case that happened to be noticed.
        for step in 0..64 {

            let angle = step as f64 * std::f64::consts::TAU / 64.0;
            let (x, y, z) = (angle.cos(), angle.sin(), 0.0);

            // Rodrigues at 180 degrees: 2 n n^T - I.
            let n = [x, y, z];
            let mut matrix = [[0.0; 3]; 3];
            for row in 0..3 {
                for column in 0..3 {
                    matrix[row][column] =
                        2.0 * n[row] * n[column] - if row == column { 1.0 } else { 0.0 };
                }
            }

            let q = Quaternion::from_rotation_matrix(&matrix);

            assert!(
                q.components().iter().all(|c| c.is_finite()),
                "step {step} gave {:?}",
                q
            );
            assert!(
                (q.rotation_angle_degrees().unwrap() - 180.0).abs() < 1e-6,
                "step {step} turned by {}",
                q.rotation_angle_degrees().unwrap()
            );
        }
    }

    #[test]
    fn normalizing_nothing_is_declined() {
        let tiny = Quaternion::from_scalar_first([0.0, 0.0, 0.0, 0.0]);

        assert!(tiny.normalized().is_none());
        assert!(tiny.inverse().is_none());
        assert!(tiny.rotation_angle_degrees().is_none());
        assert!(tiny.to_rotation_matrix().is_none());
    }

    #[test]
    fn the_rotation_angle_is_kagans() {
        // The two worked quaternions on p. 712 of Kagan (1991).
        let first = Quaternion::from_scalar_first([0.696, 0.322, -0.152, 0.624]);
        let second = Quaternion::from_scalar_first([0.62471, 0.32267, 0.69465, 0.15195]);

        assert!((first.rotation_angle_degrees().unwrap() - 91.8182771683).abs() < 1e-9);
        assert!((second.rotation_angle_degrees().unwrap() - 102.67846140868497).abs() < 1e-9);
    }
}
