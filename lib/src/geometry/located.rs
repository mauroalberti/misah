
use crate::geometry::point::Point;

/// A value together with the place it was observed.
///
/// Position is deliberately not a field on the observations themselves. A
/// fault plane means the same thing wherever it was measured, and the forward
/// model that consumes it never asks where it came from; giving every
/// measurement a coordinate would push a concern that belongs to one kind of
/// caller into a type that most callers use without it. Where an observation
/// sits starts to matter only when several of them are combined into a field,
/// and that is what this pairing is for.
///
/// Generic over the value as well as the dimension, so the same type carries a
/// fault plane now and a focal mechanism later without being touched.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Located<T, const N: usize> {
    pub position: Point<N>,
    pub value: T,
}

pub type Located2D<T> = Located<T, 2>;
pub type Located3D<T> = Located<T, 3>;

impl<T, const N: usize> Located<T, N> {

    pub fn new(position: Point<N>, value: T) -> Self {
        Self { position, value }
    }

    /// Distance from where this was observed to somewhere else -- the argument
    /// a distance-decaying weight is computed from.
    pub fn distance_to(&self, place: &Point<N>) -> f64 {
        self.position.distance(place)
    }

    /// The same position, holding a borrow of the value.
    ///
    /// Lets a caller pass the value on without cloning it, which matters when
    /// the value is a fault with its slickenlines and the caller is walking a
    /// grid of nodes over the whole dataset.
    pub fn as_ref(&self) -> Located<&T, N> {
        Located { position: self.position, value: &self.value }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn distance_to_measures_from_the_observation_to_the_place() {
        let obs = Located::new(Point::<2>::from([0.0, 0.0]), "a fault");

        assert_eq!(obs.distance_to(&Point::from([3.0, 4.0])), 5.0);
    }

    #[test]
    fn distance_to_its_own_position_is_zero() {
        let obs = Located::new(Point::<3>::from([1.0, 2.0, 3.0]), 7u8);

        assert_eq!(obs.distance_to(&obs.position), 0.0);
    }

    #[test]
    fn the_third_coordinate_counts_in_the_distance() {
        // Depth is not decoration: two hypocentres under the same map point
        // are far apart, and a weight built on this must say so.
        let obs = Located::new(Point::<3>::from([0.0, 0.0, 0.0]), ());

        assert_eq!(obs.distance_to(&Point::from([0.0, 0.0, -12.0])), 12.0);
    }

    #[test]
    fn as_ref_keeps_the_position_and_borrows_the_value() {
        let obs = Located::new(Point::<2>::from([5.0, 6.0]), String::from("slickenline"));

        let borrowed = obs.as_ref();

        assert_eq!(borrowed.position, obs.position);
        assert_eq!(borrowed.value, &obs.value);
    }

    #[test]
    fn a_located_value_is_the_value_it_was_given() {
        let obs = Located2D::new(Point::from([0.0, 0.0]), 42);

        assert_eq!(obs.value, 42);
    }
}
