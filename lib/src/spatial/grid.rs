//! A regular grid of the places a field is evaluated at.
//!
//! Deliberately a grid of **sample points** and not of cells. The distinction
//! sounds pedantic until it costs a half cell: a raster format that says
//! `ORIGIN` may mean the corner of the first cell or the centre of it, and code
//! that declares one and computes the other produces a field displaced by half
//! a spacing in every direction -- shifted, not wrong-looking, and so
//! survivable for years. It is the defect this module's C++ ancestor carried.
//!
//! There are no cells here to be ambiguous about. `origin` *is* the first node,
//! the rest follow at `spacing`, and a VTK `STRUCTURED_POINTS` header built
//! from `origin`, `spacing` and `counts` describes exactly the points that were
//! computed.
//!
//! Nodes are ordered with the **first axis varying fastest** -- x, then y, then
//! z -- which is what VTK and Fortran readers expect. Note that
//! `structural::best_fit` numbers its cells row-major instead, rows being the
//! natural reading order of a map sheet; the two are different conventions on
//! purpose, and neither is a grid of the other's kind.

use crate::geometry::point::Point;

/// Where a field is sampled: a first node, a step along each axis, and how
/// many nodes each axis carries.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SamplingGrid<const N: usize> {
    origin: Point<N>,
    spacing: [f64; N],
    counts: [usize; N],
}

impl<const N: usize> SamplingGrid<N> {

    /// `None` unless every spacing is finite and positive, every count is at
    /// least one, and the origin is finite.
    pub fn new(origin: Point<N>, spacing: [f64; N], counts: [usize; N]) -> Option<Self> {

        if N == 0 {
            return None;
        }
        if spacing.iter().any(|s| !s.is_finite() || *s <= 0.0) {
            return None;
        }
        if counts.contains(&0) {
            return None;
        }
        if origin.coords.iter().any(|c| !c.is_finite()) {
            return None;
        }

        // A grid large enough to overflow `usize` when flattened cannot be
        // indexed, let alone evaluated; better to refuse it here than to wrap
        // an index somewhere in the middle of a long run.
        counts.iter().try_fold(1usize, |total, c| total.checked_mul(*c))?;

        Some(Self { origin, spacing, counts })
    }

    /// A grid covering `points`, extended by `margin` map units on every side.
    ///
    /// The margin is what keeps a density field from being cut off at the
    /// outermost observation, where the estimate is still falling away; a
    /// margin of about one bandwidth shows the whole of that decay.
    ///
    /// `None` for an empty set, a non-positive or non-finite spacing or
    /// margin, or coordinates that are not finite.
    pub fn covering(points: &[Point<N>], spacing: [f64; N], margin: f64) -> Option<Self> {

        if points.is_empty() || !margin.is_finite() || margin < 0.0 {
            return None;
        }
        if spacing.iter().any(|s| !s.is_finite() || *s <= 0.0) {
            return None;
        }
        if points.iter().any(|p| p.coords.iter().any(|c| !c.is_finite())) {
            return None;
        }

        let mut lower = [f64::INFINITY; N];
        let mut upper = [f64::NEG_INFINITY; N];

        for point in points {
            for i in 0..N {
                lower[i] = lower[i].min(point.coords[i]);
                upper[i] = upper[i].max(point.coords[i]);
            }
        }

        let mut origin = [0.0; N];
        let mut counts = [0usize; N];

        for i in 0..N {
            origin[i] = lower[i] - margin;
            let span = (upper[i] + margin) - origin[i];
            // Nodes, not intervals: one more than the number of steps, and at
            // least one for a set that is flat along this axis.
            counts[i] = (span / spacing[i]).floor() as usize + 1;
        }

        Self::new(Point::from(origin), spacing, counts)
    }

    pub fn origin(&self) -> &Point<N> {
        &self.origin
    }

    pub fn spacing(&self) -> &[f64; N] {
        &self.spacing
    }

    /// Nodes along each axis.
    pub fn counts(&self) -> &[usize; N] {
        &self.counts
    }

    /// Nodes in the whole grid -- how many times a field will evaluate
    /// whatever it is evaluating.
    pub fn node_count(&self) -> usize {
        self.counts.iter().product()
    }

    /// The node at a set of per-axis indices.
    pub fn node(&self, index: [usize; N]) -> Option<Point<N>> {

        if (0..N).any(|i| index[i] >= self.counts[i]) {
            return None;
        }

        Some(Point::from(std::array::from_fn(|i| {
            self.origin.coords[i] + index[i] as f64 * self.spacing[i]
        })))
    }

    /// The per-axis indices of a node given its position in the flat ordering,
    /// first axis fastest.
    pub fn index_of(&self, flat: usize) -> Option<[usize; N]> {

        if flat >= self.node_count() {
            return None;
        }

        let mut remaining = flat;

        Some(std::array::from_fn(|i| {
            let along = remaining % self.counts[i];
            remaining /= self.counts[i];
            along
        }))
    }

    /// The node at a position in the flat ordering.
    pub fn node_at(&self, flat: usize) -> Option<Point<N>> {
        self.node(self.index_of(flat)?)
    }

    /// Every node, in the flat ordering, paired with its position.
    pub fn nodes(&self) -> impl Iterator<Item = (usize, Point<N>)> + '_ {
        (0..self.node_count()).map(move |flat| {
            (flat, self.node_at(flat).expect("flat index is within the grid"))
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_origin_is_the_first_node_and_not_a_corner() {
        // The whole reason this type exists: no half spacing anywhere, so a
        // header written from `origin` and `spacing` names the points that
        // were computed.
        let grid =
            SamplingGrid::<2>::new(Point::from([100.0, 200.0]), [10.0, 10.0], [3, 2]).unwrap();

        assert_eq!(grid.node([0, 0]), Some(Point::from([100.0, 200.0])));
        assert_eq!(grid.node([2, 1]), Some(Point::from([120.0, 210.0])));
    }

    #[test]
    fn the_first_axis_varies_fastest() {
        let grid =
            SamplingGrid::<2>::new(Point::from([0.0, 0.0]), [1.0, 1.0], [3, 2]).unwrap();

        let positions: Vec<Point<2>> = grid.nodes().map(|(_, p)| p).collect();

        assert_eq!(
            positions,
            vec![
                Point::from([0.0, 0.0]),
                Point::from([1.0, 0.0]),
                Point::from([2.0, 0.0]),
                Point::from([0.0, 1.0]),
                Point::from([1.0, 1.0]),
                Point::from([2.0, 1.0]),
            ]
        );
    }

    #[test]
    fn flat_indices_and_per_axis_indices_are_inverses() {
        let grid =
            SamplingGrid::<3>::new(Point::from([0.0, 0.0, 0.0]), [1.0, 2.0, 4.0], [4, 3, 2])
                .unwrap();

        assert_eq!(grid.node_count(), 24);

        for flat in 0..grid.node_count() {
            let index = grid.index_of(flat).expect("within the grid");
            assert_eq!(
                grid.node(index),
                grid.node_at(flat),
                "flat index {} disagreed with its own decomposition",
                flat
            );
        }

        assert!(grid.index_of(24).is_none());
        assert!(grid.node_at(24).is_none());
        assert!(grid.node([4, 0, 0]).is_none());
    }

    #[test]
    fn a_covering_grid_holds_every_point_it_was_built_from() {
        let points = vec![
            Point::<2>::from([12.0, -3.0]),
            Point::<2>::from([57.0, 41.0]),
            Point::<2>::from([30.0, 20.0]),
        ];

        let grid = SamplingGrid::covering(&points, [5.0, 5.0], 10.0).unwrap();

        let last = grid
            .node(std::array::from_fn(|i| grid.counts()[i] - 1))
            .unwrap();

        for point in &points {
            for i in 0..2 {
                assert!(
                    point.coords[i] >= grid.origin().coords[i],
                    "point {:?} fell below the grid on axis {}",
                    point,
                    i
                );
                assert!(
                    point.coords[i] <= last.coords[i],
                    "point {:?} fell past the grid on axis {}",
                    point,
                    i
                );
            }
        }

        // And the margin is honoured on the low side exactly.
        assert_eq!(grid.origin(), &Point::from([2.0, -13.0]));
    }

    #[test]
    fn a_set_flat_along_one_axis_still_gets_a_node_on_it() {
        // Hypocentres all at one depth, or a map-view field: the degenerate
        // axis collapses to a single node rather than to none.
        let points = vec![
            Point::<3>::from([0.0, 0.0, -8000.0]),
            Point::<3>::from([100.0, 100.0, -8000.0]),
        ];

        let grid = SamplingGrid::covering(&points, [50.0, 50.0, 50.0], 0.0).unwrap();

        assert_eq!(grid.counts()[2], 1);
        assert!(grid.node_count() > 0);
    }

    #[test]
    fn a_grid_needs_a_positive_spacing_and_a_node_on_every_axis() {
        let origin = Point::<2>::from([0.0, 0.0]);

        assert!(SamplingGrid::new(origin, [0.0, 1.0], [2, 2]).is_none());
        assert!(SamplingGrid::new(origin, [-1.0, 1.0], [2, 2]).is_none());
        assert!(SamplingGrid::new(origin, [f64::NAN, 1.0], [2, 2]).is_none());
        assert!(SamplingGrid::new(origin, [1.0, 1.0], [0, 2]).is_none());
        assert!(SamplingGrid::new(Point::from([f64::NAN, 0.0]), [1.0, 1.0], [2, 2]).is_none());
        assert!(SamplingGrid::new(origin, [1.0, 1.0], [1, 1]).is_some());
    }

    #[test]
    fn a_grid_that_could_not_be_indexed_is_refused() {
        // Rather than wrapping a flat index somewhere in the middle of a run
        // that was going to take days anyway.
        assert!(SamplingGrid::<3>::new(
            Point::from([0.0, 0.0, 0.0]),
            [1.0, 1.0, 1.0],
            [usize::MAX, 2, 2]
        )
        .is_none());
    }

    #[test]
    fn covering_nothing_is_not_a_grid() {
        assert!(SamplingGrid::<2>::covering(&[], [1.0, 1.0], 0.0).is_none());
        assert!(SamplingGrid::covering(
            &[Point::<2>::from([0.0, 0.0])],
            [1.0, 1.0],
            -1.0
        )
        .is_none());
        assert!(SamplingGrid::covering(
            &[Point::<2>::from([f64::NAN, 0.0])],
            [1.0, 1.0],
            0.0
        )
        .is_none());
    }
}
