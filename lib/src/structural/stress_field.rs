//! A stress tensor at every node of a grid, with the density of the data that
//! produced it.
//!
//! The pass the rest of this was built for. `inversion::invert_weighted` gives
//! one tensor from a weighted fault set; `spatial::kernel` turns a distance
//! into a weight; `spatial::density` finds what is near a place. Put together
//! they are Hardebeck and Michael's (2006) spatially varying inversion in its
//! kernel form: rather than cutting the dataset into bins and damping the
//! neighbouring solutions towards each other, every fault contributes to every
//! node it can reach, by an amount that falls off with distance. Where the data
//! are dense the estimate is local; where they thin out the same kernel reaches
//! further into what there is, and the effective sample size says so.
//!
//! ## Density and tensor in the same traversal
//!
//! They are computed together because they answer each other. A tensor at a
//! node where four faults reached is not the same object as one where two
//! hundred did, and reading a stress field without the density beside it is
//! reading a smoothed picture of the data distribution and calling it tectonics
//! -- the ordinary failure of this kind of map. Sharing the traversal is also
//! what makes it cheap: the neighbourhood is found once and both passes read
//! it.
//!
//! It does *not* follow that both belong on the same grid. Density costs a
//! kernel evaluation per fault per node; a tensor costs a whole grid search,
//! tens of thousands of candidate tensors each resolved onto every fault in
//! reach. Measured on 300 faults over 40 km of map, a 6 km quartic kernel and
//! the default search grid, that came to about 350 ms per node inverted --
//! against microseconds for the density at the same node. A density field
//! worth looking at may want a spacing a stress field could not afford by
//! several orders of magnitude. So `density_field` remains available on its
//! own, to be run finely, and this pass is for the coarse grid where a tensor
//! is actually wanted -- carrying the density at those same nodes, which is
//! the number that belongs beside those tensors.
//!
//! The compact kernel is what makes it affordable at all. In that same run a
//! node saw about 22 of the 300 faults, so a 100-node field cost four times a
//! single inversion over the whole set rather than a hundred times it. Call
//! `field_cost` first: it walks the same nodes with the same kernel in well
//! under a millisecond and reports the work the real pass would do.
//!
//! ## What `N` is
//!
//! The dimension of the *field*, not of the geology. A fault plane and its
//! slickenlines are three-dimensional whatever grid they are posted on. `N` is
//! how many coordinates say where the fault was found: two for a map of surface
//! measurements, three for hypocentres or for a dataset with depth. The same
//! function serves both.

use crate::geometry::located::Located;
use crate::geometry::point::Point;
use crate::spatial::density::Neighbourhood;
use crate::spatial::grid::SamplingGrid;
use crate::spatial::kernel::Kernel;
use crate::structural::fault::FaultPlane;
use crate::structural::inversion::{InversionResult, SearchGrid, invert_weighted};

/// What was found at one node.
#[derive(Debug, Clone)]
pub struct FieldNode<const N: usize> {
    pub position: Point<N>,
    /// Faults per unit volume -- per unit area in two dimensions -- from the
    /// same kernel that produced the weights. Independent of whether a tensor
    /// was inverted here, and worth mapping on its own: it is the honest
    /// picture of where the data are.
    pub density: f64,
    /// Kish's effective sample size of the kernel weights: `(sum w)^2 / sum
    /// w^2`, the number of faults' worth of information that reached this
    /// node. Twenty faults at equal weight give twenty; twenty of which
    /// nineteen are on the far edge of the kernel give barely one.
    ///
    /// This is the number the threshold is applied to, and the one to read
    /// beside a tensor. Note it is not the `effective_sample_size` inside
    /// `solution`: that one counts only the faults the winning tensor could
    /// actually be scored on, which is fewer wherever a candidate resolved no
    /// shear on a plane.
    pub support: f64,
    /// How many faults reached this node at all, at any weight.
    pub faults_within_reach: usize,
    /// The inversion at this node, or `None` where it was not run or found
    /// nothing: no fault within reach, support below the threshold, or a
    /// neighbourhood that could not be scored.
    ///
    /// The whole result rather than the winning tensor alone, because the
    /// runners-up are what say whether the minimum stood by itself. On a grid
    /// they are usually its own neighbours, a step away in S1; a runner-up
    /// that fits nearly as well from right across the sphere is two tectonic
    /// phases meeting at this node, and dropping it to save a few bytes per
    /// node would hide the one thing a stress field is read for.
    pub solution: Option<InversionResult>,
}

impl<const N: usize> FieldNode<N> {

    /// Mean angular misfit of the winning tensor, where there is one.
    pub fn misfit_degrees(&self) -> Option<f64> {
        self.solution.as_ref().map(|s| s.best.mean_misfit_degrees)
    }

    /// Trend and plunge of S1, where there is a tensor.
    pub fn s1(&self) -> Option<(f64, f64)> {
        self.solution.as_ref().map(|s| {
            let axis = crate::structural::geol_axis::GeologicalAxis::from_versor(
                &s.best.tensor.s1_versor(),
            );
            (axis.trend, axis.plunge)
        })
    }
}

/// A field, and what became of the nodes.
#[derive(Debug, Clone)]
pub struct StressField<const N: usize> {
    /// One per grid node, in the grid's own flat order -- first axis fastest.
    pub nodes: Vec<FieldNode<N>>,
    pub grid: SamplingGrid<N>,
    /// Nodes that came back with a tensor.
    pub nodes_inverted: usize,
    /// Nodes where faults reached but not enough of them: `support` below the
    /// threshold. They keep their density, which is what makes the hole
    /// legible rather than merely blank.
    pub nodes_below_threshold: usize,
    /// Nodes no fault reached at all.
    pub nodes_without_data: usize,
}

/// What a field will cost, measured rather than guessed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FieldCost {
    pub nodes: usize,
    /// Nodes that will actually be inverted, the rest being below the support
    /// threshold or empty.
    pub nodes_to_invert: usize,
    /// Candidate tensors resolved onto a fault, summed over those nodes: the
    /// unit of work this problem is made of, and the number that says whether
    /// a run is minutes or days.
    ///
    /// Measured at about four and a half million a second on one core of an
    /// ordinary desktop, in release build -- and measured on both a single
    /// inversion and a whole field, which came out within seven per cent of
    /// each other. That agreement is what makes this the right unit to count:
    /// the work per forward solution does not depend on how the solutions were
    /// arranged into nodes.
    pub forward_solutions: u128,
}

/// Invert a stress tensor at every node of `grid`, and report the data density
/// there.
///
/// `kernel` turns the distance from a node to a fault into that fault's weight
/// in the inversion at that node, and its bandwidth is the scale the field is
/// smoothed over -- the choice that matters most and the one no data can make
/// for you. A compact kernel, `Kernel::quartic` or a truncated Gaussian, is
/// what keeps the cost of a node proportional to its own neighbourhood; an
/// untruncated Gaussian will work and will resolve every candidate tensor onto
/// every fault in the dataset at every node.
///
/// `min_support` is the effective sample size below which a node is left
/// without a tensor. Zero inverts wherever a single fault reaches, which is
/// almost never what is wanted: four parameters are being estimated, so a node
/// held up by three faults' worth of weight has a solution in the sense that
/// the search returns one and in no other sense. Something in the range of ten
/// is a defensible floor for real data.
///
/// Nodes with no fault in reach are not an error; a field is entitled to holes,
/// and they are counted rather than filled.
pub fn stress_field<const N: usize>(
    faults: &[Located<FaultPlane, N>],
    grid: &SamplingGrid<N>,
    kernel: Kernel<N>,
    search: SearchGrid,
    min_support: f64,
) -> StressField<N> {

    let neighbourhood = Neighbourhood::new(faults, kernel);

    let mut nodes = Vec::with_capacity(grid.node_count());
    let mut nodes_inverted = 0usize;
    let mut nodes_below_threshold = 0usize;
    let mut nodes_without_data = 0usize;

    // Reused across nodes: a field is thousands of these, and the shapes are
    // the same each time.
    let mut found: Vec<(usize, f64)> = Vec::new();
    let mut nearby: Vec<FaultPlane> = Vec::new();
    let mut weights: Vec<f64> = Vec::new();

    for (_, position) in grid.nodes() {

        neighbourhood.near_into(&position, &mut found);

        let density = found.iter().map(|(_, weight)| weight).sum();
        let support = effective_sample_size(&found);

        if found.is_empty() {
            nodes_without_data += 1;
            nodes.push(FieldNode {
                position,
                density,
                support,
                faults_within_reach: 0,
                solution: None,
            });
            continue;
        }

        if support < min_support {
            nodes_below_threshold += 1;
            nodes.push(FieldNode {
                position,
                density,
                support,
                faults_within_reach: found.len(),
                solution: None,
            });
            continue;
        }

        // Only the faults that reach this node are handed to the search.
        // `invert_weighted` would skip the zero-weighted ones before solving
        // anything, so passing the whole dataset would give the same answer --
        // but it would walk it once per candidate tensor, and there are tens of
        // thousands of those per node. Copying a neighbourhood is the cheaper
        // side of that trade by a wide margin.
        nearby.clear();
        weights.clear();
        for &(index, weight) in found.iter() {
            nearby.push(faults[index].value.clone());
            weights.push(weight);
        }

        let solution = invert_weighted(&nearby, &weights, search);
        if solution.is_some() {
            nodes_inverted += 1;
        }

        nodes.push(FieldNode {
            position,
            density,
            support,
            faults_within_reach: found.len(),
            solution,
        });
    }

    StressField {
        nodes,
        grid: *grid,
        nodes_inverted,
        nodes_below_threshold,
        nodes_without_data,
    }
}

/// Size a field before running it.
///
/// Runs the cheap half of `stress_field` -- the neighbour search and the
/// support at every node -- and reports how much of the expensive half would
/// follow. Not an estimate from a formula: it visits the same nodes with the
/// same kernel and counts the faults that actually reach each one, so a
/// dataset clustered in one corner of its own grid is costed as such rather
/// than as if it were spread evenly.
///
/// The arguments are `stress_field`'s, less the search grid's effect on
/// anything but the count.
pub fn field_cost<const N: usize>(
    faults: &[Located<FaultPlane, N>],
    grid: &SamplingGrid<N>,
    kernel: Kernel<N>,
    search: SearchGrid,
    min_support: f64,
) -> FieldCost {

    let neighbourhood = Neighbourhood::new(faults, kernel);
    let candidates = search.candidate_count() as u128;

    let mut nodes_to_invert = 0usize;
    let mut forward_solutions = 0u128;
    let mut found: Vec<(usize, f64)> = Vec::new();

    for (_, position) in grid.nodes() {

        neighbourhood.near_into(&position, &mut found);

        if found.is_empty() || effective_sample_size(&found) < min_support {
            continue;
        }

        nodes_to_invert += 1;
        forward_solutions += candidates * found.len() as u128;
    }

    FieldCost {
        nodes: grid.node_count(),
        nodes_to_invert,
        forward_solutions,
    }
}

/// Kish's effective sample size of a set of weights.
///
/// Zero for an empty neighbourhood, which is the honest answer: no faults'
/// worth of information reached here.
fn effective_sample_size(found: &[(usize, f64)]) -> f64 {

    let mut sum = 0.0;
    let mut sum_of_squares = 0.0;

    for &(_, weight) in found {
        sum += weight;
        sum_of_squares += weight * weight;
    }

    if sum_of_squares <= 0.0 {
        return 0.0;
    }

    sum * sum / sum_of_squares
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::orientation::direction::Direction3D;
    use crate::spatial::kernel::Bandwidth;
    use crate::structural::geol_axis::GeologicalAxis;
    use crate::structural::geol_plane::GeologicalPlane;
    use crate::structural::slickenline::{SlipSense, Slickenline};
    use crate::structural::stress::ReducedStressTensor;

    const SPREAD: [(f64, f64); 8] = [
        (0.0, 60.0), (45.0, 70.0), (90.0, 50.0), (135.0, 65.0),
        (180.0, 55.0), (225.0, 75.0), (270.0, 45.0), (315.0, 60.0),
    ];

    fn normal_phase() -> ReducedStressTensor {
        ReducedStressTensor::normalized(
            GeologicalAxis::new(0.0, 90.0),
            GeologicalAxis::new(90.0, 0.0),
            0.5,
        )
        .unwrap()
    }

    fn strike_slip_phase() -> ReducedStressTensor {
        ReducedStressTensor::normalized(
            GeologicalAxis::new(0.0, 0.0),
            GeologicalAxis::new(90.0, 0.0),
            0.5,
        )
        .unwrap()
    }

    /// Faults slipping exactly as `tensor` says they should, all placed at
    /// `at` -- so the only thing the field has to get right is which of them
    /// reaches which node.
    fn faults_at(tensor: &ReducedStressTensor, at: [f64; 2]) -> Vec<Located<FaultPlane, 2>> {

        SPREAD
            .iter()
            .filter_map(|&(strike, dip)| {
                let plane = GeologicalPlane::from_rhr_strike(strike, dip);
                let predicted = tensor.solve(&plane).theoretical_slickenline?;
                let slick = Slickenline::new(Direction3D::new(predicted), Some(SlipSense::Down));
                let fault = FaultPlane::new(plane, vec![slick]).ok()?;

                Some(Located::new(Point::from(at), fault))
            })
            .collect()
    }

    fn s1_angle(result: &InversionResult, truth: &ReducedStressTensor) -> f64 {
        result
            .best
            .tensor
            .s1_versor()
            .dot(&truth.s1_versor())
            .abs()
            .clamp(-1.0, 1.0)
            .acos()
            .to_degrees()
    }

    #[test]
    fn each_node_recovers_the_phase_that_is_near_it() {
        // The case the whole module exists for. Two tectonic phases, one in
        // the west and one in the east, far enough apart that a compact kernel
        // never sees both. A single inversion over this dataset returns a
        // tensor belonging to neither; a field returns each in its own place.
        let extension = normal_phase();
        let transcurrence = strike_slip_phase();

        let mut faults = faults_at(&extension, [0.0, 0.0]);
        faults.extend(faults_at(&transcurrence, [10_000.0, 0.0]));

        let kernel = Kernel::quartic(Bandwidth::isotropic(2_000.0).unwrap());
        let grid =
            SamplingGrid::<2>::new(Point::from([0.0, 0.0]), [10_000.0, 1.0], [2, 1]).unwrap();

        let field = stress_field(&faults, &grid, kernel, SearchGrid::default(), 4.0);

        assert_eq!(field.nodes.len(), 2);
        assert_eq!(field.nodes_inverted, 2);

        let west = field.nodes[0].solution.as_ref().expect("a tensor in the west");
        let east = field.nodes[1].solution.as_ref().expect("a tensor in the east");

        assert!(
            s1_angle(west, &extension) <= 10.0 + 1e-6,
            "the western node found S1 {} degrees off the extensional phase",
            s1_angle(west, &extension)
        );
        assert!(
            s1_angle(east, &transcurrence) <= 10.0 + 1e-6,
            "the eastern node found S1 {} degrees off the strike-slip phase",
            s1_angle(east, &transcurrence)
        );

        // Each node saw only its own cluster.
        assert_eq!(field.nodes[0].faults_within_reach, 8);
        assert_eq!(field.nodes[1].faults_within_reach, 8);
    }

    #[test]
    fn a_node_out_of_reach_keeps_its_place_in_the_field() {
        let faults = faults_at(&normal_phase(), [0.0, 0.0]);
        let kernel = Kernel::quartic(Bandwidth::isotropic(1_000.0).unwrap());

        // Three nodes: one on the data, two well away from it.
        let grid =
            SamplingGrid::<2>::new(Point::from([0.0, 0.0]), [50_000.0, 1.0], [3, 1]).unwrap();

        let field = stress_field(&faults, &grid, kernel, SearchGrid::default(), 4.0);

        assert_eq!(field.nodes.len(), 3);
        assert_eq!(field.nodes_inverted, 1);
        assert_eq!(field.nodes_without_data, 2);

        // The holes are holes, not zeros passed off as answers.
        assert!(field.nodes[1].solution.is_none());
        assert_eq!(field.nodes[1].density, 0.0);
        assert_eq!(field.nodes[1].support, 0.0);
        assert_eq!(field.nodes[1].faults_within_reach, 0);
    }

    #[test]
    fn thin_support_is_declined_rather_than_answered() {
        let faults = faults_at(&normal_phase(), [0.0, 0.0]);
        let kernel = Kernel::quartic(Bandwidth::isotropic(1_000.0).unwrap());
        let grid = SamplingGrid::<2>::new(Point::from([0.0, 0.0]), [1.0, 1.0], [1, 1]).unwrap();

        // Eight faults reach the node, so a threshold above eight refuses it
        // and one below accepts it -- the same node, the same data.
        let refused = stress_field(&faults, &grid, kernel, SearchGrid::default(), 9.0);
        assert_eq!(refused.nodes_below_threshold, 1);
        assert_eq!(refused.nodes_inverted, 0);
        assert!(refused.nodes[0].solution.is_none());
        // And it still knows what it saw.
        assert_eq!(refused.nodes[0].faults_within_reach, 8);
        assert!(refused.nodes[0].density > 0.0);

        let accepted = stress_field(&faults, &grid, kernel, SearchGrid::default(), 4.0);
        assert_eq!(accepted.nodes_inverted, 1);
    }

    #[test]
    fn support_is_the_weighted_count_and_not_the_plain_one() {
        // Faults at one place, a node offset far enough that the kernel has
        // fallen a long way: every fault still reaches, so the count is
        // unchanged, but they all reach weakly and equally -- which is the
        // case the effective size cannot distinguish from full weight, and
        // should not.
        let faults = faults_at(&normal_phase(), [0.0, 0.0]);
        let kernel = Kernel::quartic(Bandwidth::isotropic(1_000.0).unwrap());
        let grid = SamplingGrid::<2>::new(Point::from([900.0, 0.0]), [1.0, 1.0], [1, 1]).unwrap();

        let field = stress_field(&faults, &grid, kernel, SearchGrid::default(), 0.0);
        let node = &field.nodes[0];

        assert_eq!(node.faults_within_reach, 8);
        assert!((node.support - 8.0).abs() < 1e-9, "support {}", node.support);
        // Weak but present: the density has fallen, which is the number that
        // does notice.
        assert!(node.density > 0.0);

        // Now a node between two clusters at very different distances. Here
        // the count says sixteen and the support says close to eight.
        let mut mixed = faults_at(&normal_phase(), [0.0, 0.0]);
        mixed.extend(faults_at(&strike_slip_phase(), [990.0, 0.0]));
        let grid = SamplingGrid::<2>::new(Point::from([990.0, 0.0]), [1.0, 1.0], [1, 1]).unwrap();

        let field = stress_field(&mixed, &grid, kernel, SearchGrid::default(), 0.0);
        let node = &field.nodes[0];

        assert_eq!(node.faults_within_reach, 16);
        assert!(
            node.support < 9.0,
            "support {} should be near the eight faults that are actually close",
            node.support
        );
    }

    #[test]
    fn the_density_is_the_same_number_the_density_field_would_give() {
        // The two passes share a kernel and a neighbourhood, and this is the
        // assertion that keeps them sharing it: a stress field's density
        // column must be the density field, not a variant of it.
        let mut faults = faults_at(&normal_phase(), [0.0, 0.0]);
        faults.extend(faults_at(&strike_slip_phase(), [1_500.0, 500.0]));

        let kernel = Kernel::quartic(Bandwidth::isotropic(2_000.0).unwrap());
        let grid =
            SamplingGrid::<2>::new(Point::from([-500.0, -500.0]), [500.0, 500.0], [6, 4]).unwrap();

        let field = stress_field(&faults, &grid, kernel, SearchGrid::default(), 1e9);
        let separately = crate::spatial::density::density_field(&faults, &grid, kernel);

        assert_eq!(field.nodes.len(), separately.len());
        for (node, expected) in field.nodes.iter().zip(separately) {
            assert_eq!(node.density, expected);
        }
    }

    #[test]
    fn the_cost_counts_the_nodes_the_field_goes_on_to_invert() {
        let mut faults = faults_at(&normal_phase(), [0.0, 0.0]);
        faults.extend(faults_at(&strike_slip_phase(), [3_000.0, 0.0]));

        let kernel = Kernel::quartic(Bandwidth::isotropic(1_200.0).unwrap());
        let grid =
            SamplingGrid::<2>::new(Point::from([-2_000.0, 0.0]), [1_000.0, 1.0], [8, 1]).unwrap();
        let search = SearchGrid { angle_step_degrees: 30.0, phi_step: 0.5 };

        let cost = field_cost(&faults, &grid, kernel, search, 4.0);
        let field = stress_field(&faults, &grid, kernel, search, 4.0);

        assert_eq!(cost.nodes, grid.node_count());
        assert_eq!(cost.nodes_to_invert, field.nodes_inverted);

        // And the work is the candidates times the faults in reach, node by
        // node -- the sum the field itself performs.
        let expected: u128 = field
            .nodes
            .iter()
            .filter(|n| n.solution.is_some())
            .map(|n| search.candidate_count() as u128 * n.faults_within_reach as u128)
            .sum();
        assert_eq!(cost.forward_solutions, expected);

        // A run that would invert nothing is costed at nothing, which is the
        // answer worth having before starting.
        let none = field_cost(&faults, &grid, kernel, search, 1e9);
        assert_eq!(none.nodes_to_invert, 0);
        assert_eq!(none.forward_solutions, 0);
    }

    #[test]
    fn a_field_over_hypocentral_depths_works_the_same_way() {
        // Three dimensions, with a bandwidth that reaches ten times as far
        // across the map as it does down: a seismogenic layer is wide and
        // thin, and an isotropic kernel would mix faults from top to bottom of
        // it before mixing two neighbouring ones.
        let extension = normal_phase();

        let shallow: Vec<Located<FaultPlane, 3>> = faults_at(&extension, [0.0, 0.0])
            .into_iter()
            .map(|f| Located::new(Point::from([0.0, 0.0, -2_000.0]), f.value))
            .collect();
        let deep: Vec<Located<FaultPlane, 3>> = faults_at(&strike_slip_phase(), [0.0, 0.0])
            .into_iter()
            .map(|f| Located::new(Point::from([0.0, 0.0, -12_000.0]), f.value))
            .collect();

        let mut faults = shallow;
        faults.extend(deep);

        let kernel = Kernel::quartic(Bandwidth::new([20_000.0, 20_000.0, 3_000.0]).unwrap());
        let grid = SamplingGrid::<3>::new(
            Point::from([0.0, 0.0, -12_000.0]),
            [1.0, 1.0, 10_000.0],
            [1, 1, 2],
        )
        .unwrap();

        let field = stress_field(&faults, &grid, kernel, SearchGrid::default(), 4.0);

        assert_eq!(field.nodes_inverted, 2);
        // Wide open across the map, each node still saw only its own depth.
        assert_eq!(field.nodes[0].faults_within_reach, 8);
        assert_eq!(field.nodes[1].faults_within_reach, 8);

        let at_depth = field.nodes[0].solution.as_ref().unwrap();
        let near_surface = field.nodes[1].solution.as_ref().unwrap();

        assert!(s1_angle(at_depth, &strike_slip_phase()) <= 10.0 + 1e-6);
        assert!(s1_angle(near_surface, &extension) <= 10.0 + 1e-6);
    }

    #[test]
    fn a_field_over_no_faults_is_all_holes() {
        let faults: Vec<Located<FaultPlane, 2>> = Vec::new();
        let kernel = Kernel::quartic(Bandwidth::isotropic(1_000.0).unwrap());
        let grid = SamplingGrid::<2>::new(Point::from([0.0, 0.0]), [100.0, 100.0], [3, 3]).unwrap();

        let field = stress_field(&faults, &grid, kernel, SearchGrid::default(), 4.0);

        assert_eq!(field.nodes.len(), 9);
        assert_eq!(field.nodes_without_data, 9);
        assert_eq!(field.nodes_inverted, 0);
        assert!(field.nodes.iter().all(|n| n.solution.is_none()));
    }

    #[test]
    fn faults_without_slickenlines_leave_a_node_without_a_tensor() {
        // They reach, so they count towards the density and the support; they
        // just cannot be scored. A node like this is not a hole in the data
        // and should not be reported as one.
        let faults: Vec<Located<FaultPlane, 2>> = SPREAD
            .iter()
            .map(|&(strike, dip)| {
                Located::new(
                    Point::from([0.0, 0.0]),
                    FaultPlane::without_slickenlines(GeologicalPlane::from_rhr_strike(strike, dip)),
                )
            })
            .collect();

        let kernel = Kernel::quartic(Bandwidth::isotropic(1_000.0).unwrap());
        let grid = SamplingGrid::<2>::new(Point::from([0.0, 0.0]), [1.0, 1.0], [1, 1]).unwrap();

        let field = stress_field(&faults, &grid, kernel, SearchGrid::default(), 4.0);

        assert_eq!(field.nodes_inverted, 0);
        assert_eq!(field.nodes_without_data, 0);
        assert_eq!(field.nodes_below_threshold, 0);
        assert!(field.nodes[0].solution.is_none());
        assert_eq!(field.nodes[0].faults_within_reach, 8);
        assert!(field.nodes[0].density > 0.0);
    }

    #[test]
    fn the_node_helpers_read_off_the_solution() {
        let faults = faults_at(&normal_phase(), [0.0, 0.0]);
        let kernel = Kernel::quartic(Bandwidth::isotropic(1_000.0).unwrap());
        let grid = SamplingGrid::<2>::new(Point::from([0.0, 0.0]), [1.0, 1.0], [1, 1]).unwrap();

        let field = stress_field(&faults, &grid, kernel, SearchGrid::default(), 4.0);
        let node = &field.nodes[0];

        let (_, plunge) = node.s1().expect("a tensor");
        // The extensional phase has S1 vertical, and the grid steps ten
        // degrees.
        assert!(plunge >= 80.0, "S1 plunge {}", plunge);
        assert!(node.misfit_degrees().expect("a tensor") < 10.0);
    }

    #[test]
    fn a_field_keeps_the_grid_it_was_given() {
        let faults = faults_at(&normal_phase(), [0.0, 0.0]);
        let kernel = Kernel::quartic(Bandwidth::isotropic(1_000.0).unwrap());
        let grid = SamplingGrid::<2>::new(Point::from([10.0, 20.0]), [5.0, 5.0], [4, 3]).unwrap();

        let field = stress_field(&faults, &grid, kernel, SearchGrid::default(), 1e9);

        assert_eq!(field.grid, grid);
        assert_eq!(field.nodes.len(), grid.node_count());
        // In the grid's order, so a node can be found by its flat index.
        for (flat, place) in grid.nodes() {
            assert_eq!(field.nodes[flat].position, place);
        }
    }
}
