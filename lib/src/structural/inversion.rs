//! Fault-slip inversion: recover the reduced stress tensor that best explains
//! a set of measured faults.
//!
//! The inverse of `structural::stress`, and built directly on it. A candidate
//! tensor is scored by resolving it onto every fault and comparing the slip it
//! predicts against the slip observed, which is one forward solution per
//! candidate per fault -- the loop the forward model's own documentation calls
//! the reason it lives in a compiled crate rather than only in geogst. A grid
//! of a few tens of thousands of candidates over a hundred faults is a few
//! million forward solutions, which is where that stops being a claim.
//!
//! The search is exhaustive rather than iterative, and deliberately so: the
//! misfit surface of this problem is not convex, a fault set with two
//! superposed tectonic phases has two minima by construction, and a descent
//! from an arbitrary start would report whichever one it happened to fall
//! into without ever saying the other existed. An exhaustive pass over a known
//! grid returns a result whose worst case is the grid step, and `best` and
//! `runners_up` together let a caller see whether the minimum stood alone.

use crate::structural::fault::FaultPlane;
use crate::structural::geol_axis::GeologicalAxis;
use crate::structural::stress::ReducedStressTensor;

/// How finely the space of candidate tensors is sampled.
///
/// The whole space is three angles and a ratio: the orientation of S1 (trend
/// and plunge), the roll of S3 about it, and Phi. S2 follows from the other
/// two axes, and the magnitudes do not matter -- the predicted slip direction
/// depends only on the shape of the tensor, which is why the search runs on
/// normalized ones.
#[derive(Debug, Clone, Copy)]
pub struct SearchGrid {
    /// Step, in degrees, for the S1 trend and plunge and for the S3 roll.
    pub angle_step_degrees: f64,
    /// Step for Phi, over `[0, 1]`.
    pub phi_step: f64,
}

impl Default for SearchGrid {
    /// Ten degrees and 0.1: about 21 000 candidates, coarse enough to run in a
    /// moment and fine enough to place a tensor within the scatter of the
    /// measurements it is inverted from, which no fault-slip dataset beats.
    fn default() -> Self {
        Self { angle_step_degrees: 10.0, phi_step: 0.1 }
    }
}

impl SearchGrid {

    /// Candidate tensors this grid will visit, before any are rejected as
    /// degenerate. Useful for sizing a run before starting it.
    pub fn candidate_count(&self) -> usize {

        let trends = (360.0 / self.angle_step_degrees).ceil() as usize;
        // Plunge spans a quarter turn: an axis pointing up is the same axis
        // pointing down, so the downward hemisphere holds every distinct S1.
        let plunges = (90.0 / self.angle_step_degrees).floor() as usize + 1;
        // S3 is likewise an axis, so half a turn about S1 exhausts it.
        let rolls = (180.0 / self.angle_step_degrees).ceil() as usize;
        let phis = (1.0 / self.phi_step).floor() as usize + 1;

        trends * plunges * rolls * phis
    }
}

/// One scored candidate.
#[derive(Debug, Clone)]
pub struct ScoredTensor {
    pub tensor: ReducedStressTensor,
    /// Mean angular misfit, in degrees, over the faults it could be scored on.
    pub mean_misfit_degrees: f64,
    /// How many faults that was. A candidate scored on fewer faults than
    /// another is not straightforwardly better for having a lower mean.
    pub faults_scored: usize,
}

/// What a search found.
#[derive(Debug, Clone)]
pub struct InversionResult {
    pub best: ScoredTensor,
    /// The next few candidates by misfit, worst last. A minimum that stands
    /// alone looks quite different here from one on a plateau, or from two
    /// tectonic phases the mean has averaged into a tensor belonging to
    /// neither.
    pub runners_up: Vec<ScoredTensor>,
    pub candidates_tried: usize,
    pub candidates_scored: usize,
}

/// How many runners-up `invert` keeps.
const RUNNERS_UP: usize = 5;

/// Search the grid for the tensor of least mean misfit.
///
/// `None` when no candidate could be scored on any fault -- an empty set, or
/// one whose faults carry no slickenlines.
pub fn invert(faults: &[FaultPlane], grid: SearchGrid) -> Option<InversionResult> {

    let mut scored: Vec<ScoredTensor> = Vec::new();
    let mut candidates_tried = 0usize;

    for (s1_axis, s3_axis) in candidate_axes(grid) {

        let mut phi: f64 = 0.0;
        while phi <= 1.0 + 1e-9 {

            candidates_tried += 1;

            // A candidate whose axes fall outside the tensor's own
            // orthogonality tolerance is not a failure of the search: it is a
            // pair the grid produced and the tensor declines, and it simply
            // does not compete.
            if let Ok(tensor) = ReducedStressTensor::normalized(s1_axis, s3_axis, phi.min(1.0)) {
                if let Some((mean_misfit_degrees, faults_scored)) = tensor.mean_misfit(faults) {
                    scored.push(ScoredTensor {
                        tensor,
                        mean_misfit_degrees,
                        faults_scored,
                    });
                }
            }

            phi += grid.phi_step;
        }
    }

    if scored.is_empty() {
        return None;
    }

    let candidates_scored = scored.len();

    // Sorted by misfit, and by fault count where two candidates tie: a tensor
    // that explains more of the dataset equally well is the better answer.
    scored.sort_by(|a, b| {
        a.mean_misfit_degrees
            .partial_cmp(&b.mean_misfit_degrees)
            .expect("misfits are finite angles")
            .then(b.faults_scored.cmp(&a.faults_scored))
    });

    let best = scored.remove(0);
    scored.truncate(RUNNERS_UP);

    Some(InversionResult {
        best,
        runners_up: scored,
        candidates_tried,
        candidates_scored,
    })
}

/// The (S1, S3) pairs the grid visits.
///
/// S3 is generated by rolling a reference perpendicular about S1, so the pair
/// is orthogonal by construction rather than by search: sampling two axes
/// independently and discarding the pairs that miss orthogonality would spend
/// most of the grid on candidates that never compete.
fn candidate_axes(grid: SearchGrid) -> Vec<(GeologicalAxis, GeologicalAxis)> {

    let step = grid.angle_step_degrees;
    let mut pairs = Vec::new();

    let mut plunge: f64 = 0.0;
    while plunge <= 90.0 + 1e-9 {

        let mut trend: f64 = 0.0;
        while trend < 360.0 - 1e-9 {

            let s1 = GeologicalAxis::new(trend, plunge);
            let s1_versor = s1.as_versor();

            // Any vector not along S1 gives a starting perpendicular; the
            // vertical serves unless S1 is itself vertical, where east does.
            let seed = if plunge.abs() > 89.0 {
                GeologicalAxis::new(90.0, 0.0).as_versor()
            } else {
                GeologicalAxis::new(0.0, 90.0).as_versor()
            };

            // The seed is chosen never to be parallel to S1, so the cross
            // product cannot vanish; if it somehow did, this S1 contributes no
            // candidates rather than skipping the trend increment below, which
            // would leave the loop spinning on it forever.
            if let Some(u) = s1_versor.as_vector().cross(&seed.as_vector()).normalize() {

                let v = s1_versor.as_vector().cross(&u);

                let mut roll: f64 = 0.0;
                while roll < 180.0 - 1e-9 {

                    let a = roll.to_radians();
                    let s3_vector = u * a.cos() + v * a.sin();

                    if let Ok(s3_versor) = crate::algebra::versor::Versor3D::new(s3_vector.coords) {
                        pairs.push((s1, GeologicalAxis::from_versor(&s3_versor)));
                    }

                    roll += step;
                }
            }

            trend += step;
        }

        // At a vertical S1 every trend names the same axis, so the trend loop
        // above has already covered the pole once; stepping past it would
        // repeat the whole roll sweep for nothing.
        if plunge >= 90.0 - 1e-9 {
            break;
        }
        plunge += step;
    }

    pairs
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::structural::geol_plane::GeologicalPlane;
    use crate::structural::slickenline::{SlipSense, Slickenline};

    /// Faults slipping exactly as a known tensor says they should.
    ///
    /// Built through the forward model, so the inversion is being asked to
    /// recover a tensor from data that tensor itself generated -- the one
    /// case where the right answer is known independently of the search.
    fn faults_from(tensor: &ReducedStressTensor, planes: &[(f64, f64)]) -> Vec<FaultPlane> {

        planes
            .iter()
            .filter_map(|&(strike, dip)| {
                let plane = GeologicalPlane::from_rhr_strike(strike, dip);
                let solution = tensor.solve(&plane);

                solution.theoretical_slickenline.map(|predicted| {
                    let slick = Slickenline::new(
                        crate::orientation::direction::Direction3D::new(predicted),
                        Some(SlipSense::Down),
                    );
                    FaultPlane::new(plane, vec![slick])
                        .expect("a predicted slip lies in the plane it was predicted on")
                })
            })
            .collect()
    }

    /// A spread of orientations, so the inversion is constrained rather than
    /// left free by a dataset that all dips one way.
    const SPREAD: [(f64, f64); 8] = [
        (0.0, 60.0), (45.0, 70.0), (90.0, 50.0), (135.0, 65.0),
        (180.0, 55.0), (225.0, 75.0), (270.0, 45.0), (315.0, 60.0),
    ];

    #[test]
    fn a_tensor_is_recovered_from_the_slip_it_predicts() {
        // Andersonian normal faulting: S1 vertical, S3 horizontal east-west.
        let truth = ReducedStressTensor::normalized(
            GeologicalAxis::new(0.0, 90.0),
            GeologicalAxis::new(90.0, 0.0),
            0.5,
        )
        .unwrap();

        let faults = faults_from(&truth, &SPREAD);
        assert!(faults.len() >= 6, "the synthetic set is too small to constrain anything");

        let result = invert(&faults, SearchGrid::default()).expect("a scorable set");

        // The data were generated without error, so the true tensor is on the
        // grid and the misfit at it is zero; anything above a grid step means
        // the search found a different tensor, not a nearby one.
        assert!(
            result.best.mean_misfit_degrees < 10.0,
            "misfit {} degrees",
            result.best.mean_misfit_degrees
        );
        assert_eq!(result.best.faults_scored, faults.len());

        // S1 recovered to within the grid step, as an axis: an axis and its
        // antipode name the same one.
        let found = result.best.tensor.s1_versor();
        let angle = found.dot(&truth.s1_versor()).abs().clamp(-1.0, 1.0).acos().to_degrees();
        assert!(angle <= 10.0 + 1e-6, "S1 off by {} degrees", angle);
    }

    #[test]
    fn a_strike_slip_tensor_is_recovered_too() {
        // S1 and S3 both horizontal: a setting the vertical-S1 case above
        // could pass without ever exercising.
        let truth = ReducedStressTensor::normalized(
            GeologicalAxis::new(0.0, 0.0),
            GeologicalAxis::new(90.0, 0.0),
            0.5,
        )
        .unwrap();

        let faults = faults_from(&truth, &SPREAD);
        let result = invert(&faults, SearchGrid::default()).expect("a scorable set");

        assert!(
            result.best.mean_misfit_degrees < 10.0,
            "misfit {} degrees",
            result.best.mean_misfit_degrees
        );

        let found = result.best.tensor.s1_versor();
        let angle = found.dot(&truth.s1_versor()).abs().clamp(-1.0, 1.0).acos().to_degrees();
        assert!(angle <= 10.0 + 1e-6, "S1 off by {} degrees", angle);
    }

    #[test]
    fn a_finer_grid_does_not_fit_worse() {
        let truth = ReducedStressTensor::normalized(
            GeologicalAxis::new(20.0, 75.0),
            // Perpendicular to S1 exactly: in the same vertical plane, a
            // quarter turn past it. (110, 5) reads plausibly and is 4.8
            // degrees off, which ReducedStressTensor rightly refuses.
            GeologicalAxis::new(200.0, 15.0),
            0.4,
        )
        .unwrap();
        let faults = faults_from(&truth, &SPREAD);

        let coarse = invert(&faults, SearchGrid { angle_step_degrees: 30.0, phi_step: 0.5 })
            .expect("a scorable set");
        let fine = invert(&faults, SearchGrid::default()).expect("a scorable set");

        assert!(
            fine.best.mean_misfit_degrees <= coarse.best.mean_misfit_degrees + 1e-9,
            "the finer grid fitted worse: {} against {}",
            fine.best.mean_misfit_degrees,
            coarse.best.mean_misfit_degrees
        );
        assert!(fine.candidates_tried > coarse.candidates_tried);
    }

    #[test]
    fn faults_without_slickenlines_cannot_be_inverted() {
        let planes = SPREAD
            .iter()
            .map(|&(strike, dip)| {
                FaultPlane::without_slickenlines(GeologicalPlane::from_rhr_strike(strike, dip))
            })
            .collect::<Vec<_>>();

        assert!(invert(&planes, SearchGrid::default()).is_none());
    }

    #[test]
    fn an_empty_set_inverts_to_nothing() {
        assert!(invert(&[], SearchGrid::default()).is_none());
    }

    #[test]
    fn the_result_reports_what_it_searched() {
        let truth = ReducedStressTensor::normalized(
            GeologicalAxis::new(0.0, 90.0),
            GeologicalAxis::new(90.0, 0.0),
            0.5,
        )
        .unwrap();
        let faults = faults_from(&truth, &SPREAD);

        let grid = SearchGrid { angle_step_degrees: 30.0, phi_step: 0.5 };
        let result = invert(&faults, grid).expect("a scorable set");

        assert!(result.candidates_scored <= result.candidates_tried);
        assert!(result.candidates_tried <= grid.candidate_count());
        assert!(!result.runners_up.is_empty());
        // Sorted: every runner-up fits no better than the winner.
        for other in &result.runners_up {
            assert!(other.mean_misfit_degrees >= result.best.mean_misfit_degrees - 1e-9);
        }
    }

    #[test]
    fn an_unread_slip_sense_does_not_penalise_a_fault() {
        // The same faults, once with the sense recorded and once without, and
        // with half the lineations deliberately stored pointing the wrong way
        // along their own line -- which is exactly what an unread sense means.
        let truth = ReducedStressTensor::normalized(
            GeologicalAxis::new(0.0, 90.0),
            GeologicalAxis::new(90.0, 0.0),
            0.5,
        )
        .unwrap();

        let flipped: Vec<FaultPlane> = SPREAD
            .iter()
            .enumerate()
            .filter_map(|(i, &(strike, dip))| {
                let plane = GeologicalPlane::from_rhr_strike(strike, dip);
                let predicted = truth.solve(&plane).theoretical_slickenline?;

                let direction = crate::orientation::direction::Direction3D::new(predicted);
                let direction = if i % 2 == 0 { direction } else { direction.opposite() };

                FaultPlane::new(plane, vec![Slickenline::new(direction, None)]).ok()
            })
            .collect();

        let misfit = truth.mean_misfit(&flipped).expect("a scorable set").0;

        assert!(
            misfit < 1e-6,
            "reversed lineations with no sense should fit perfectly, got {} degrees",
            misfit
        );
    }
}
