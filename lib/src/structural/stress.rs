//! The direct (forward) Wallace-Bott problem: given a reduced stress tensor,
//! predict the slip a fault plane would show under it.
//!
//! The inverse of what a fault-slip inversion (right-dihedra, Angelier, ...)
//! solves for. Ported from `ForwardStress.f95`, a Fortran tool of Alberti
//! (2010) that once lived in this repository's own history under
//! `GeoFaults`, before being dropped for a Rust port that was never carried
//! through -- this is that port, arriving by way of a parallel one in geogst
//! first, where `Plane.rake_to_direct` turned out to already implement the
//! Fortran's rake-to-slickenline formula exactly (both are Aki & Richards,
//! 1980), and where the whole algorithm was cross-checked against the
//! Fortran source itself, compiled and run unmodified alongside a corrected
//! build: `vector_processing`, in the original, carries no `implicit none`,
//! so `vector1_magn` inside `vector_normalization` -- reached through
//! `vector_projection`, used for the normal-stress vector every solution goes
//! through -- is an undeclared name, silently typed as a single-precision
//! `REAL` by Fortran's implicit-typing rule rather than the double precision
//! the rest of the tool uses. The truncation is usually too small to matter,
//! but the predicted rake reaches it through an `acos` close to the edge of
//! its domain in near-strike-slip cases, where the error is four orders of
//! magnitude larger in the angle than in the cosine. This port -- and the
//! corrected Fortran build -- agree with a 50-digit independent
//! recomputation of the same cases to 12+ significant digits; the original
//! binary does not.
//!
//! `misah` has no per-fault performance problem to justify a compiled kernel
//! here on its own -- this direct calculation is O(1) per fault, not the
//! O(n) grid search of the inverse problem. It is here as the inner loop that
//! search would call, evaluated once per candidate tensor per fault.

use thiserror::Error;

use crate::algebra::versor::Versor3D;
use crate::structural::geol_axis::GeologicalAxis;
use crate::structural::geol_plane::GeologicalPlane;

/// Below this, a shear stress is numerical noise rather than a driving force
/// -- the fault plane sits on, or acutely close to, a principal stress axis,
/// and no slip direction is defined. Xu (2004) and the Fortran original both
/// treat such a plane as having no valid solution rather than an arbitrary
/// one.
pub const SHEAR_MAGNITUDE_THRESHOLD: f64 = 1.0e-5;

/// Tolerance, in degrees, on how far from exactly orthogonal `S1` and `S3`
/// may be and still be accepted.
pub const AXIS_ORTHOGONALITY_TOLERANCE: f64 = 1.0;

#[derive(Debug, Error)]
pub enum StressError {
    #[error(
        "S1 and S3 axes must be sub-orthogonal: {angle:.3} degrees apart, \
         expected within {tolerance:.3} of 90"
    )]
    AxesNotOrthogonal { angle: f64, tolerance: f64 },

    #[error("Phi must be between 0 and 1, got {0}")]
    PhiOutOfRange(f64),

    #[error("Sigma1 must be greater than Sigma3, got sigma1={sigma1}, sigma3={sigma3}")]
    InvalidMagnitudes { sigma1: f64, sigma3: f64 },
}

/// A reduced stress tensor (Angelier, 1984): the orientations of the S1 and
/// S3 principal stress axes, plus the shape ratio
/// `Phi = (sigma2 - sigma3) / (sigma1 - sigma3)`, which between them are what
/// a fault-slip inversion recovers. S2 completes a right-handed orthonormal
/// triad with S1 and S3, cyclically -- `S1 x S2 = S3`, `S2 x S3 = S1` -- and
/// its magnitude follows from Phi.
///
/// Absolute magnitudes are optional: `sigma1`/`sigma3` default to 1 and 0,
/// the usual normalization when only the shape of the tensor is known, which
/// leaves ratios -- the predicted rake among them -- correct, while making
/// the absolute ones (slip tendency, deformation index) meaningless. Supply
/// the true magnitudes when they are known and those are wanted.
#[derive(Debug, Clone)]
pub struct ReducedStressTensor {
    s1: Versor3D,
    s2: Versor3D,
    s3: Versor3D,
    pub sigma1: f64,
    pub sigma3: f64,
    pub phi: f64,
}

impl ReducedStressTensor {

    pub fn new(
        s1_axis: GeologicalAxis,
        s3_axis: GeologicalAxis,
        phi: f64,
        sigma1: f64,
        sigma3: f64,
    ) -> Result<Self, StressError> {

        let s1 = s1_axis.as_versor();
        let s3 = s3_axis.as_versor();

        let angle = s1.dot(&s3).clamp(-1.0, 1.0).acos().to_degrees();
        if (angle - 90.0).abs() > AXIS_ORTHOGONALITY_TOLERANCE {
            return Err(StressError::AxesNotOrthogonal {
                angle,
                tolerance: AXIS_ORTHOGONALITY_TOLERANCE,
            });
        }

        if !(0.0..=1.0).contains(&phi) {
            return Err(StressError::PhiOutOfRange(phi));
        }

        if sigma1 <= sigma3 {
            return Err(StressError::InvalidMagnitudes { sigma1, sigma3 });
        }

        // S3 cross S1, not S1 cross S3: this order is what makes S1, S2, S3 a
        // right-handed cyclic triad, matching the Fortran original. A versor
        // rather than the raw cross product: exact orthogonality (the
        // ordinary case) makes no difference, but input merely
        // sub-orthogonal within the tolerance above would otherwise leave S2
        // imperceptibly short, and the triad this builds is meant to be
        // orthonormal.
        let s2 = s3.as_vector().cross(&s1.as_vector());
        let s2 = Versor3D::new(s2.coords).expect(
            "S1 and S3 were just confirmed sub-orthogonal, so their cross product is non-zero",
        );

        Ok(Self { s1, s2, s3, sigma1, sigma3, phi })
    }

    /// Build directly from `sigma1 = 1`, `sigma3 = 0`: the usual
    /// normalization when only the shape of the tensor -- the orientation of
    /// S1/S3 and Phi -- is known, as from a fault-slip inversion.
    pub fn normalized(
        s1_axis: GeologicalAxis,
        s3_axis: GeologicalAxis,
        phi: f64,
    ) -> Result<Self, StressError> {
        Self::new(s1_axis, s3_axis, phi, 1.0, 0.0)
    }

    pub fn s1_versor(&self) -> Versor3D { self.s1 }

    /// S3 cross S1, completing a right-handed triad.
    pub fn s2_versor(&self) -> Versor3D { self.s2 }

    pub fn s3_versor(&self) -> Versor3D { self.s3 }

    /// The intermediate principal stress magnitude, from sigma1, sigma3 and phi.
    pub fn sigma2(&self) -> f64 {
        self.phi * self.sigma1 + (1.0 - self.phi) * self.sigma3
    }

    /// The stress tensor, as a 3x3 array in (East, North, Up) components.
    ///
    /// `R . diag(sigma1, sigma2, sigma3) . R^T`, with `R`'s columns the S1,
    /// S2, S3 versors.
    pub fn tensor(&self) -> [[f64; 3]; 3] {

        let s = [self.s1.coords(), self.s2.coords(), self.s3.coords()];
        let sigma = [self.sigma1, self.sigma2(), self.sigma3];

        let mut t = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                t[i][j] = (0..3).map(|k| s[k][i] * sigma[k] * s[k][j]).sum();
            }
        }
        t
    }

    /// Apply this tensor to a fault plane, predicting the slip it drives.
    ///
    /// Follows Xu (2004): the traction on the plane is resolved into a
    /// normal and a shear component; where the shear is not negligible, its
    /// direction gives the predicted rake and slickenline. Where it is (the
    /// plane sits on, or acutely close to, a principal stress axis), the
    /// solution carries the stresses but no direction, since none is
    /// defined.
    pub fn solve(&self, plane: &GeologicalPlane) -> StressSolution {
        self.solve_with_threshold(plane, SHEAR_MAGNITUDE_THRESHOLD)
    }

    pub fn solve_with_threshold(
        &self,
        plane: &GeologicalPlane,
        shear_threshold: f64,
    ) -> StressSolution {

        let t = self.tensor();

        // The forward-pointing normal -- upward for a shallow-dipping plane,
        // horizontal for a vertical one -- rather than the downward
        // `normal_axis`: the sign of every signed quantity below (traction,
        // normal stress) is relative to this choice, so it is the one the
        // theoretical rake further down is built consistently against.
        let n = plane.normal_vector();
        let n = [n.coords[0], n.coords[1], n.coords[2]];

        let traction = matvec(&t, n).map(|c| -c);
        let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        let norm = |a: [f64; 3]| dot(a, a).sqrt();

        let mut traction_magn = norm(traction);
        if dot(n, traction) < 0.0 {
            traction_magn = -traction_magn;
        }

        let normal_scal = dot(n, traction);
        let normal_stress = n.map(|c| c * normal_scal);
        let normal_stress_magn = normal_scal;

        let shear_stress = std::array::from_fn(|i| traction[i] - normal_stress[i]);
        let shear_stress_magn = norm(shear_stress);

        let is_valid = shear_stress_magn > shear_threshold;

        let mut theoretical_rake = None;
        let mut theoretical_slickenline = None;
        let mut slip_tendency = None;
        let mut deformation_index = None;

        if is_valid {
            let shear_versor = shear_stress.map(|c| c / shear_stress_magn);

            let strike_versor = GeologicalAxis::new(plane.rhr_strike(), 0.0).as_versor();
            let dip_versor = GeologicalAxis::new(plane.azimuth, plane.dip_angle).as_versor();

            let scal_strike = dot(shear_versor, *strike_versor.coords()).clamp(-1.0, 1.0);
            let scal_dip = dot(shear_versor, *dip_versor.coords());

            let mut rake = scal_strike.acos().to_degrees();
            if scal_dip > 0.0 {
                rake = -rake;
            }

            theoretical_slickenline = Some(plane.rake_to_versor(rake));
            theoretical_rake = Some(rake);

            slip_tendency = Some(shear_stress_magn / traction_magn.abs());
            // abs() in the numerator alone, not the denominator: Xu (2004)'s
            // index as the Fortran original carried it. The asymmetry is
            // kept rather than "fixed" -- deformation_index changes sign
            // along with traction_magn, which flips only when the traction
            // opposes the forward normal, a distinction the index is
            // apparently meant to preserve.
            deformation_index = Some((traction_magn.abs() - shear_stress_magn) / traction_magn);
        }

        StressSolution {
            traction,
            traction_magnitude: traction_magn,
            normal_stress,
            normal_stress_magnitude: normal_stress_magn,
            shear_stress,
            shear_stress_magnitude: shear_stress_magn,
            is_valid,
            theoretical_rake,
            theoretical_slickenline,
            slip_tendency,
            deformation_index,
        }
    }
}

fn matvec(m: &[[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| (0..3).map(|j| m[i][j] * v[j]).sum())
}

/// The result of resolving a `ReducedStressTensor` onto one fault plane: the
/// traction and its normal/shear components always, in (East, North, Up),
/// and, where the shear is large enough to drive a defined slip, the
/// predicted rake, slickenline, slip tendency and deformation index.
#[derive(Debug, Clone)]
pub struct StressSolution {
    pub traction: [f64; 3],
    pub traction_magnitude: f64,
    pub normal_stress: [f64; 3],
    pub normal_stress_magnitude: f64,
    pub shear_stress: [f64; 3],
    pub shear_stress_magnitude: f64,
    pub is_valid: bool,
    pub theoretical_rake: Option<f64>,
    pub theoretical_slickenline: Option<Versor3D>,
    pub slip_tendency: Option<f64>,
    pub deformation_index: Option<f64>,
}

impl StressSolution {

    /// The angle between the predicted and an observed slickenline -- what a
    /// fault-slip inversion minimizes, evaluated here in the forward
    /// direction. `None` when this solution has no predicted direction to
    /// compare against.
    pub fn angular_misfit(&self, observed: &Versor3D) -> Option<f64> {
        self.theoretical_slickenline
            .as_ref()
            .map(|predicted| predicted.dot(observed).clamp(-1.0, 1.0).acos().to_degrees())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Five cases checked against ForwardStress.f95 itself: the modules up
    /// to `stress_processing`, compiled unmodified with gfortran and driven
    /// by a small program built for the purpose (calling `S2_calc`,
    /// `sigma2_calc`, `rotmatr`, `stresstensorcalc` and `stresssolution_calc`
    /// directly, bypassing the interactive/file-based main program).
    ///
    /// Two of the five (cases 2 and 3) match the *corrected* binary, not the
    /// one as it stands -- see the module doc comment for the
    /// `implicit none` bug that separates them, confirmed independently by a
    /// 50-digit `mpmath` recomputation of the same five inputs. These are
    /// the same five cases and the same values geogst's port of this tool
    /// was checked against.
    mod against_the_fortran_original {
        use super::*;

        fn are_close(a: f64, b: f64) -> bool {
            (a - b).abs() <= 1e-9_f64.max(1e-9 * a.abs().max(b.abs()))
        }

        #[test]
        fn case_1_anderson_normal() {
            let tensor = ReducedStressTensor::new(
                GeologicalAxis::new(0.0, 90.0), GeologicalAxis::new(90.0, 0.0), 0.5, 30.0, 10.0,
            ).unwrap();
            let plane = GeologicalPlane::from_rhr_strike(0.0, 60.0);

            let sol = tensor.solve(&plane);

            assert!(sol.is_valid);
            assert!(are_close(sol.traction_magnitude, -17.320508075689));
            assert!(are_close(sol.normal_stress_magnitude, -15.0));
            assert!(are_close(sol.shear_stress_magnitude, 8.660254037844));
            assert!(are_close(sol.theoretical_rake.unwrap(), -90.0));
            assert!(are_close(sol.slip_tendency.unwrap(), 0.5));
            assert!(are_close(sol.deformation_index.unwrap(), -0.5));

            let slick = GeologicalAxis::from_versor(&sol.theoretical_slickenline.unwrap());
            assert!(are_close(slick.trend, 90.0));
            assert!(are_close(slick.plunge, 60.0));
        }

        #[test]
        fn case_2_oblique() {
            let tensor = ReducedStressTensor::new(
                GeologicalAxis::new(0.0, 0.0), GeologicalAxis::new(90.0, 0.0), 0.3, 40.0, -10.0,
            ).unwrap();
            let plane = GeologicalPlane::from_rhr_strike(30.0, 70.0);

            let sol = tensor.solve(&plane);

            assert!(sol.is_valid);
            assert!(are_close(sol.traction_magnitude, -20.551398971889));
            assert!(are_close(sol.normal_stress_magnitude, -2.792444446101));
            assert!(are_close(sol.shear_stress_magnitude, 20.360801892784));
            assert!(are_close(sol.theoretical_rake.unwrap(), -2.261611727892));
            assert!(are_close(sol.slip_tendency.unwrap(), 0.990725834316));
            assert!(are_close(sol.deformation_index.unwrap(), -0.009274165684));

            let slick = GeologicalAxis::from_versor(&sol.theoretical_slickenline.unwrap());
            assert!(are_close(slick.trend, 30.773871690337));
            assert!(are_close(slick.plunge, 2.125155259309));
        }

        #[test]
        fn case_3_plunging_axes() {
            let tensor = ReducedStressTensor::new(
                GeologicalAxis::new(125.0, 35.0), GeologicalAxis::new(125.0, -55.0), 0.6, 35.0, 5.0,
            ).unwrap();
            let plane = GeologicalPlane::from_rhr_strike(200.0, 50.0);

            let sol = tensor.solve(&plane);

            assert!(sol.is_valid);
            assert!(are_close(sol.traction_magnitude, -34.425635706498));
            assert!(are_close(sol.normal_stress_magnitude, -34.215382549860));
            assert!(are_close(sol.shear_stress_magnitude, 3.798946006885));
            assert!(are_close(sol.theoretical_rake.unwrap(), 43.558963074529));
            assert!(are_close(sol.slip_tendency.unwrap(), 0.110352239804));
            assert!(are_close(sol.deformation_index.unwrap(), -0.889647760196));

            let slick = GeologicalAxis::from_versor(&sol.theoretical_slickenline.unwrap());
            assert!(are_close(slick.trend, 168.565012691676));
            assert!(are_close(slick.plunge, -31.862445102249));
        }

        #[test]
        fn case_4_fault_normal_parallel_to_a_principal_axis() {
            let tensor = ReducedStressTensor::new(
                GeologicalAxis::new(10.0, 5.0), GeologicalAxis::new(10.0, -85.0), 0.4, 50.0, -5.0,
            ).unwrap();
            let plane = GeologicalPlane::from_rhr_strike(100.0, 85.0);

            let sol = tensor.solve(&plane);

            assert!(!sol.is_valid);
            assert!(are_close(sol.traction_magnitude, -50.0));
            assert!(are_close(sol.normal_stress_magnitude, -50.0));
            assert!(sol.shear_stress_magnitude.abs() < 1e-9);
            assert!(sol.theoretical_rake.is_none());
            assert!(sol.theoretical_slickenline.is_none());
            assert!(sol.slip_tendency.is_none());
            assert!(sol.deformation_index.is_none());
        }

        #[test]
        fn case_5_horizontal_fault_under_vertical_s1() {
            let tensor = ReducedStressTensor::new(
                GeologicalAxis::new(0.0, 90.0), GeologicalAxis::new(90.0, 0.0), 0.5, 30.0, 10.0,
            ).unwrap();
            let plane = GeologicalPlane::from_rhr_strike(0.0, 0.0);

            let sol = tensor.solve(&plane);

            assert!(!sol.is_valid);
            assert!(are_close(sol.traction_magnitude, -30.0));
            assert!(are_close(sol.normal_stress_magnitude, -30.0));
            assert!(sol.shear_stress_magnitude.abs() < 1e-9);
        }
    }

    #[test]
    fn sigma2_from_phi_sigma1_sigma3() {
        let tensor = ReducedStressTensor::new(
            GeologicalAxis::new(0.0, 90.0), GeologicalAxis::new(90.0, 0.0), 0.5, 30.0, 10.0,
        ).unwrap();
        assert_eq!(tensor.sigma2(), 20.0);
    }

    #[test]
    fn default_magnitudes_leave_the_rake_unaffected() {
        // Sigma1=1, sigma3=0 is only a normalization: the predicted rake and
        // slickenline, which depend on the shear *direction* and not on the
        // traction's absolute scale, must not move when the magnitudes do.
        let s1 = GeologicalAxis::new(125.0, 35.0);
        let s3 = GeologicalAxis::new(125.0, -55.0);
        let plane = GeologicalPlane::from_rhr_strike(200.0, 50.0);

        let normalized = ReducedStressTensor::normalized(s1, s3, 0.6).unwrap().solve(&plane);
        let scaled = ReducedStressTensor::new(s1, s3, 0.6, 35.0, 5.0).unwrap().solve(&plane);

        let a = normalized.theoretical_rake.unwrap();
        let b = scaled.theoretical_rake.unwrap();
        assert!((a - b).abs() < 1e-9, "{} vs {}", a, b);

        let misfit = normalized
            .angular_misfit(&scaled.theoretical_slickenline.unwrap())
            .unwrap();
        assert!(misfit < 1e-9);

        // But slip tendency does depend on the absolute stress state.
        assert!((normalized.slip_tendency.unwrap() - scaled.slip_tendency.unwrap()).abs() > 1e-6);
    }

    #[test]
    fn axes_must_be_sub_orthogonal() {
        let err = ReducedStressTensor::new(
            GeologicalAxis::new(0.0, 0.0), GeologicalAxis::new(45.0, 0.0), 0.5, 1.0, 0.0,
        ).unwrap_err();
        assert!(matches!(err, StressError::AxesNotOrthogonal { .. }));
    }

    #[test]
    fn phi_must_be_in_unit_range() {
        let err = ReducedStressTensor::new(
            GeologicalAxis::new(0.0, 90.0), GeologicalAxis::new(90.0, 0.0), 1.5, 1.0, 0.0,
        ).unwrap_err();
        assert!(matches!(err, StressError::PhiOutOfRange(_)));
    }

    #[test]
    fn sigma1_must_exceed_sigma3() {
        let err = ReducedStressTensor::new(
            GeologicalAxis::new(0.0, 90.0), GeologicalAxis::new(90.0, 0.0), 0.5, 10.0, 10.0,
        ).unwrap_err();
        assert!(matches!(err, StressError::InvalidMagnitudes { .. }));
    }

    #[test]
    fn s2_completes_a_right_handed_orthonormal_triad() {
        // S1 x S2 = S3, for any orthogonal S1/S3 pair -- the cyclic order the
        // whole construction relies on, not just for one axis-aligned case.
        for (s1_axis, s3_axis) in [
            (GeologicalAxis::new(0.0, 90.0), GeologicalAxis::new(90.0, 0.0)),
            (GeologicalAxis::new(125.0, 35.0), GeologicalAxis::new(125.0, -55.0)),
            (GeologicalAxis::new(40.0, 10.0), GeologicalAxis::new(40.0, -80.0)),
        ] {
            let tensor = ReducedStressTensor::normalized(s1_axis, s3_axis, 0.5).unwrap();
            let (s1, s2, s3) = (tensor.s1_versor(), tensor.s2_versor(), tensor.s3_versor());

            let cross = s1.as_vector().cross(&s2.as_vector());
            let cross = Versor3D::new(cross.coords).unwrap();

            assert!(cross.dot(&s3).clamp(-1.0, 1.0).acos().to_degrees() < 1e-6);
        }
    }

    #[test]
    fn angular_misfit_of_an_invalid_solution_is_none() {
        let tensor = ReducedStressTensor::new(
            GeologicalAxis::new(0.0, 90.0), GeologicalAxis::new(90.0, 0.0), 0.5, 30.0, 10.0,
        ).unwrap();
        let sol = tensor.solve(&GeologicalPlane::from_rhr_strike(0.0, 0.0));

        assert!(sol.angular_misfit(&GeologicalAxis::new(0.0, 0.0).as_versor()).is_none());
    }

    #[test]
    fn angular_misfit_of_the_predicted_direction_against_itself_is_zero() {
        let tensor = ReducedStressTensor::new(
            GeologicalAxis::new(0.0, 90.0), GeologicalAxis::new(90.0, 0.0), 0.5, 30.0, 10.0,
        ).unwrap();
        let sol = tensor.solve(&GeologicalPlane::from_rhr_strike(0.0, 60.0));

        let predicted = sol.theoretical_slickenline.unwrap();
        assert!(sol.angular_misfit(&predicted).unwrap() < 1e-9);
    }
}
