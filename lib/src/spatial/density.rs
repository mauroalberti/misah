//! Kernel density over located observations, and the neighbour search that
//! makes a field of it affordable.
//!
//! A port, in the loose sense, of `InterpDensity3D` -- a program that read a
//! list of hypocentres and wrote a VTK volume of events per unit volume. What
//! is kept is the estimator; what is not is the shape of the program. That one
//! summed every observation at every node, which for its own datasets was
//! perfectly sensible: a few thousand hypocentres against a few hundred
//! thousand nodes runs in seconds and no indexing would repay its own
//! complexity.
//!
//! The reason for indexing here is not density at all. It is that the same
//! neighbourhood is about to be handed to `structural::inversion`, where each
//! node costs a grid search rather than an exponential, and the difference
//! between scoring the fifty faults that reach a node and scoring the ten
//! thousand that do not is the difference between a field that runs and one
//! that does not. So the search is built once, in terms of the kernel's own
//! reach, and both passes read it.
//!
//! Which is also why a kernel that admits no reach -- the untruncated Gaussian
//! -- is not refused here but falls back to visiting everything. It gives the
//! same answer, more slowly, and it is the kernel to reach for when a result
//! is being checked against something that did the same.

use std::collections::HashMap;

use crate::geometry::located::Located;
use crate::geometry::point::Point;
use crate::spatial::grid::SamplingGrid;
use crate::spatial::kernel::Kernel;

/// Observations indexed by position, ready to be asked what lies near a place.
///
/// Borrows the observations rather than taking them: the intended caller holds
/// one dataset and walks a grid over it, and a field that cloned the data per
/// node would spend more on copying than on arithmetic.
pub struct Neighbourhood<'a, T, const N: usize> {
    observations: &'a [Located<T, N>],
    kernel: Kernel<N>,
    /// Observation indices by bin, or `None` where the kernel reaches
    /// everywhere and there is nothing a bin could exclude.
    bins: Option<HashMap<[i64; N], Vec<usize>>>,
    /// Bin edge along each axis: the kernel's reach, so that anything within
    /// reach of a place lies in that place's own bin or one adjacent to it
    /// along every axis, and the block of `3^N` bins around it holds
    /// everything. Meaningless, and unused, when `bins` is `None`.
    bin_size: [f64; N],
}

impl<'a, T, const N: usize> Neighbourhood<'a, T, N> {

    /// Index `observations` at the reach of `kernel`.
    ///
    /// The kernel is kept rather than taken again per query, so the bins and
    /// the weights cannot come to disagree about how far the kernel reaches --
    /// which would not fail loudly. It would return a field that is slightly
    /// too thin near the edges of each bin.
    pub fn new(observations: &'a [Located<T, N>], kernel: Kernel<N>) -> Self {

        let Some(bin_size) = kernel.reach() else {
            return Self { observations, kernel, bins: None, bin_size: [0.0; N] };
        };

        let mut bins: HashMap<[i64; N], Vec<usize>> = HashMap::new();

        for (index, observation) in observations.iter().enumerate() {
            // A non-finite coordinate has no bin, and would land every query
            // in the same one. It is dropped, and `indexed` reports the loss.
            if observation.position.coords.iter().any(|c| !c.is_finite()) {
                continue;
            }
            let key = bin_of(&observation.position, &bin_size);
            bins.entry(key).or_default().push(index);
        }

        Self { observations, kernel, bins: Some(bins), bin_size }
    }

    pub fn observations(&self) -> &'a [Located<T, N>] {
        self.observations
    }

    pub fn kernel(&self) -> &Kernel<N> {
        &self.kernel
    }

    /// How many observations the index actually holds.
    ///
    /// Below `observations().len()` when some carried a coordinate that was
    /// not finite. Equal to it for a kernel with no reach, which indexes
    /// nothing and skips nothing.
    pub fn indexed(&self) -> usize {

        match &self.bins {
            Some(bins) => bins.values().map(Vec::len).sum(),
            None => self.observations.len(),
        }
    }

    /// Index and kernel weight of every observation that reaches `place`,
    /// appended to `found` after clearing it.
    ///
    /// Takes the vector rather than returning one because a field asks this
    /// once per node and the allocation is worth keeping. Only strictly
    /// positive weights are reported, so the caller need not test for zero
    /// and, for a compact kernel, the length of the answer is the size of the
    /// neighbourhood rather than of the dataset.
    ///
    /// The order is deterministic -- bins are reached by key, never by
    /// iterating the map, so nothing here depends on how a `HashMap` felt like
    /// arranging itself. It is not, however, ascending by index: it is bin
    /// order, and within a bin the order the observations were given in.
    /// Which matters in one place only. Summing these weights adds the same
    /// positive terms a plain scan would, in a different sequence, so a
    /// density computed through the index can differ from one computed without
    /// it in the last bits. Sorting would buy that back and cost a sort per
    /// node, on the one loop where a pure density field spends all its time.
    pub fn near_into(&self, place: &Point<N>, found: &mut Vec<(usize, f64)>) {

        found.clear();

        let Some(bins) = &self.bins else {
            for (index, observation) in self.observations.iter().enumerate() {
                let weight = self.kernel.weight(&observation.position, place);
                if weight > 0.0 {
                    found.push((index, weight));
                }
            }
            return;
        };

        if place.coords.iter().any(|c| !c.is_finite()) {
            return;
        }

        let home = bin_of(place, &self.bin_size);

        // The block of bins touching this one, in every dimension at once:
        // nine on a plane, twenty-seven in a volume. Written as a mixed-radix
        // count rather than as nested loops, which would have to be rewritten
        // for each `N`.
        for step in 0..3usize.pow(N as u32) {

            let mut remaining = step;
            let key: [i64; N] = std::array::from_fn(|i| {
                let offset = (remaining % 3) as i64 - 1;
                remaining /= 3;
                home[i] + offset
            });

            let Some(neighbours) = bins.get(&key) else {
                continue;
            };

            for &index in neighbours {
                let weight = self.kernel.weight(&self.observations[index].position, place);
                if weight > 0.0 {
                    found.push((index, weight));
                }
            }
        }
    }

    /// The same, into a fresh vector.
    pub fn near(&self, place: &Point<N>) -> Vec<(usize, f64)> {

        let mut found = Vec::new();
        self.near_into(place, &mut found);

        found
    }

    /// Summed kernel weight at `place`: the density estimate there, in
    /// observations per unit volume.
    pub fn density_at(&self, place: &Point<N>) -> f64 {

        let mut found = Vec::new();
        self.near_into(place, &mut found);

        found.iter().map(|(_, weight)| weight).sum()
    }
}

/// Which bin a position falls in, at the given edge lengths.
fn bin_of<const N: usize>(place: &Point<N>, bin_size: &[f64; N]) -> [i64; N] {
    std::array::from_fn(|i| (place.coords[i] / bin_size[i]).floor() as i64)
}

/// Density at every node of `grid`, in the grid's flat order.
///
/// The value at a node is observations per unit volume there -- per unit area
/// in two dimensions, and so on -- so that summed over the grid and multiplied
/// by the cell volume it approaches the number of observations, less whatever
/// the kernel's own truncation drops and whatever falls outside the grid.
pub fn density_field<T, const N: usize>(
    observations: &[Located<T, N>],
    grid: &SamplingGrid<N>,
    kernel: Kernel<N>,
) -> Vec<f64> {

    let neighbourhood = Neighbourhood::new(observations, kernel);
    let mut found = Vec::new();

    grid.nodes()
        .map(|(_, place)| {
            neighbourhood.near_into(&place, &mut found);
            found.iter().map(|(_, weight)| weight).sum()
        })
        .collect()
}

/// The same for bare positions, which is the hypocentre case.
pub fn point_density_field<const N: usize>(
    points: &[Point<N>],
    grid: &SamplingGrid<N>,
    kernel: Kernel<N>,
) -> Vec<f64> {

    let located: Vec<Located<(), N>> =
        points.iter().map(|p| Located::new(*p, ())).collect();

    density_field(&located, grid, kernel)
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::spatial::kernel::Bandwidth;

    fn at<const N: usize>(coords: [f64; N]) -> Located<(), N> {
        Located::new(Point::from(coords), ())
    }

    #[test]
    fn a_single_observation_makes_the_kernel_its_own_shape() {
        let observations = [at([0.0, 0.0, 0.0])];
        let kernel = Kernel::quartic(Bandwidth::isotropic(100.0).unwrap());
        let neighbourhood = Neighbourhood::new(&observations, kernel);

        let origin = Point::<3>::from([0.0, 0.0, 0.0]);

        assert_eq!(
            neighbourhood.density_at(&origin),
            kernel.weight(&origin, &origin)
        );
        // And nothing at all beyond the bandwidth.
        assert_eq!(neighbourhood.density_at(&Point::from([150.0, 0.0, 0.0])), 0.0);
    }

    #[test]
    fn the_index_finds_exactly_what_a_scan_would() {
        // The property the binning has to have, and the one it would lose
        // quietly: a bin edge is not a boundary of the kernel.
        let mut observations = Vec::new();
        let mut seed = 12345u64;
        let mut next = || {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((seed >> 33) as f64 / (1u64 << 31) as f64) * 2000.0 - 1000.0
        };
        for _ in 0..400 {
            observations.push(at([next(), next(), next()]));
        }

        let kernel = Kernel::quartic(Bandwidth::new([300.0, 300.0, 120.0]).unwrap());
        let indexed = Neighbourhood::new(&observations, kernel);

        for step in 0..200 {
            let place = Point::<3>::from([next(), next(), next()]);

            let mut by_scan: Vec<(usize, f64)> = observations
                .iter()
                .enumerate()
                .filter_map(|(i, o)| {
                    let w = kernel.weight(&o.position, &place);
                    (w > 0.0).then_some((i, w))
                })
                .collect();

            let mut by_index = indexed.near(&place);

            by_scan.sort_by_key(|(i, _)| *i);
            by_index.sort_by_key(|(i, _)| *i);

            assert_eq!(by_index, by_scan, "the index disagreed at step {}", step);
        }
    }

    #[test]
    fn an_unbounded_kernel_visits_everything() {
        let observations = [at([0.0, 0.0, 0.0]), at([50_000.0, 0.0, 0.0])];
        let kernel = Kernel::gaussian(Bandwidth::isotropic(1000.0).unwrap());
        let neighbourhood = Neighbourhood::new(&observations, kernel);

        // Both are reported at the near one, the far contribution being
        // vanishing rather than absent -- until it underflows, which is the
        // only thing that ever silences a Gaussian.
        assert_eq!(neighbourhood.near(&Point::from([0.0, 0.0, 0.0])).len(), 1);
        assert!(neighbourhood.bins.is_none());
        assert_eq!(neighbourhood.indexed(), 2);
    }

    #[test]
    fn the_density_integrates_to_the_count() {
        // A quartic kernel keeps all its mass, so a grid fine enough to
        // resolve it and wide enough to hold it should recover the number of
        // observations to the accuracy of the quadrature.
        let observations = [
            at([0.0, 0.0, 0.0]),
            at([400.0, -200.0, 100.0]),
            at([-300.0, 250.0, -150.0]),
        ];

        let kernel = Kernel::quartic(Bandwidth::isotropic(500.0).unwrap());
        let spacing = [25.0, 25.0, 25.0];
        let grid = SamplingGrid::covering(
            &observations.iter().map(|o| o.position).collect::<Vec<_>>(),
            spacing,
            600.0,
        )
        .unwrap();

        let values = density_field(&observations, &grid, kernel);
        let cell_volume: f64 = spacing.iter().product();
        let total: f64 = values.iter().sum::<f64>() * cell_volume;

        assert!(
            (total - observations.len() as f64).abs() < 0.01,
            "the field integrated to {} observations, not {}",
            total,
            observations.len()
        );
    }

    #[test]
    fn the_field_is_returned_in_the_grids_own_order() {
        let observations = [at([0.0, 0.0])];
        let kernel = Kernel::quartic(Bandwidth::isotropic(10.0).unwrap());

        let grid =
            SamplingGrid::<2>::new(Point::from([0.0, 0.0]), [5.0, 5.0], [3, 2]).unwrap();
        let values = density_field(&observations, &grid, kernel);

        assert_eq!(values.len(), grid.node_count());
        for (flat, place) in grid.nodes() {
            assert_eq!(
                values[flat],
                kernel.weight(&observations[0].position, &place),
                "node {} did not hold its own value",
                flat
            );
        }
    }

    #[test]
    fn point_density_agrees_with_located_density() {
        let points = [
            Point::<2>::from([0.0, 0.0]),
            Point::<2>::from([30.0, 10.0]),
        ];
        let located: Vec<Located<(), 2>> = points.iter().map(|p| Located::new(*p, ())).collect();

        let kernel = Kernel::quartic(Bandwidth::isotropic(50.0).unwrap());
        let grid = SamplingGrid::covering(&points, [10.0, 10.0], 50.0).unwrap();

        assert_eq!(
            point_density_field(&points, &grid, kernel),
            density_field(&located, &grid, kernel)
        );
    }

    #[test]
    fn an_observation_with_no_finite_position_is_left_out_rather_than_spread() {
        let observations = [at([0.0, 0.0, 0.0]), at([f64::NAN, 0.0, 0.0])];
        let kernel = Kernel::quartic(Bandwidth::isotropic(100.0).unwrap());
        let neighbourhood = Neighbourhood::new(&observations, kernel);

        assert_eq!(neighbourhood.indexed(), 1);
        assert_eq!(neighbourhood.near(&Point::from([0.0, 0.0, 0.0])).len(), 1);
    }

    #[test]
    fn the_answer_does_not_depend_on_how_the_map_arranged_itself() {
        // Two indexes built from the same data are two different HashMaps.
        // Bins are reached by key rather than iterated, so the results should
        // match element for element and in order -- and a density built from
        // them should be equal bit for bit, not merely close.
        let mut observations = Vec::new();
        for i in 0..300 {
            let t = i as f64;
            observations.push(at([t * 7.0 % 900.0, t * 13.0 % 900.0, t * 3.0 % 900.0]));
        }

        let kernel = Kernel::quartic(Bandwidth::isotropic(250.0).unwrap());
        let first = Neighbourhood::new(&observations, kernel);
        let second = Neighbourhood::new(&observations, kernel);

        for step in 0..50 {
            let place = Point::<3>::from([step as f64 * 17.0, step as f64 * 11.0, 100.0]);

            assert_eq!(first.near(&place), second.near(&place));
            assert_eq!(first.density_at(&place), second.density_at(&place));
        }
    }

    #[test]
    fn an_empty_dataset_has_no_density_anywhere() {
        let observations: [Located<(), 3>; 0] = [];
        let kernel = Kernel::quartic(Bandwidth::isotropic(100.0).unwrap());
        let neighbourhood = Neighbourhood::new(&observations, kernel);

        assert_eq!(neighbourhood.density_at(&Point::from([0.0, 0.0, 0.0])), 0.0);
        assert!(neighbourhood.near(&Point::from([0.0, 0.0, 0.0])).is_empty());
    }
}
