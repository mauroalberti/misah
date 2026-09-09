//! Distance-decaying kernels, in whatever dimension the field has.
//!
//! One kernel serves two jobs that are worth keeping apart in the reader's
//! mind, because only one of them cares whether the arithmetic is normalized.
//!
//! Summed over observations, a kernel estimates a **density**: hypocentres per
//! cubic kilometre, faults per square kilometre. That number means what it says
//! only if the kernel integrates to one over the space, which is what the
//! constants below are for and what makes them depend on the dimension.
//!
//! Used one observation at a time, the same kernel supplies the **weights**
//! `structural::inversion::invert_weighted` takes, and there the normalization
//! cancels: that search reads only the relative sizes. Nothing is lost by
//! normalizing anyway, and a great deal is gained by having a single kernel
//! answer both questions, so that the tensor at a node and the density at the
//! same node were computed from the same neighbourhood with the same reach.
//!
//! ## The constants are derived, not tabulated
//!
//! A quartic kernel is `(1 - u^2)^2` inside the unit scaled radius and zero
//! outside, and the constant in front of it is whatever makes the integral one
//! -- which is a different number in each dimension: `15/16` on a line, `3/pi`
//! on a plane, `105/(32 pi)` in a volume. This is exactly the sort of constant
//! that gets copied from a two-dimensional textbook into three-dimensional code
//! and produces a field that is smooth, plausible, correctly shaped and wrong
//! by a factor. The port this module descends from, `InterpDensity3D`, carried
//! the planar constant in a volume kernel for that reason. Here the constant is
//! computed from `N`, so there is no version of it to pick the wrong one of.
//!
//! ## Bandwidth is in map units, and per axis
//!
//! Not in cells. A bandwidth expressed in cells changes meaning when the grid
//! is refined, which is precisely when a reader is looking for the answer to
//! stop changing. And it is one number per axis rather than one overall,
//! because a depth axis is not interchangeable with a horizontal one even when
//! both are in metres: hypocentral depth is the worst-determined coordinate of
//! an earthquake, and a seismogenic layer is a thing far wider than it is
//! thick. An isotropic bandwidth remains one call away for the cases that want
//! it.

use std::f64::consts::PI;

use crate::geometry::point::Point;

/// How far the kernel reaches along each axis, in the units the coordinates
/// are in.
///
/// Every extent must be finite and strictly positive: a zero bandwidth is a
/// kernel of infinite height and no width, which is not a smoothing but a
/// division by zero.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bandwidth<const N: usize> {
    extents: [f64; N],
}

impl<const N: usize> Bandwidth<N> {

    /// `None` if any extent is zero, negative or not finite.
    pub fn new(extents: [f64; N]) -> Option<Self> {

        // A zero-dimensional field has no volume to spread density over, and
        // the normalizing constants below divide by a gamma function that is
        // not defined there. Nobody writes `Bandwidth<0>` on purpose; refusing
        // it is cheaper than letting it return a NaN that travels.
        if N == 0 {
            return None;
        }

        if extents.iter().any(|e| !e.is_finite() || *e <= 0.0) {
            return None;
        }

        Some(Self { extents })
    }

    /// The same extent along every axis.
    pub fn isotropic(extent: f64) -> Option<Self> {
        Self::new([extent; N])
    }

    pub fn extents(&self) -> &[f64; N] {
        &self.extents
    }

    /// The product of the extents: the volume the normalization divides by.
    fn volume(&self) -> f64 {
        self.extents.iter().product()
    }

    /// Separation between two points measured in bandwidths rather than in map
    /// units -- the only form the kernel shapes below ever see, which is what
    /// makes an anisotropic bandwidth cost nothing to support.
    fn scaled_radius(&self, a: &Point<N>, b: &Point<N>) -> f64 {

        let mut sum = 0.0;
        for i in 0..N {
            let separation = (a.coords[i] - b.coords[i]) / self.extents[i];
            sum += separation * separation;
        }

        sum.sqrt()
    }
}

/// Which kernel, and how far it is allowed to reach.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum KernelShape {
    /// `exp(-u^2 / 2)`, unbounded. Smooth everywhere and never exactly zero,
    /// so every observation contributes to every node and the cost of a field
    /// is the grid times the whole dataset. Use it when that is affordable or
    /// when a reference implementation is being matched.
    Gaussian,
    /// The same, cut off at `radii` scaled radii and **not** renormalized, so
    /// the tail beyond the cut is simply missing; `retained_mass` says how
    /// much. At four radii that is about a thousandth in three dimensions,
    /// which buys a compact support and with it a cost per node that follows
    /// the neighbourhood rather than the dataset.
    TruncatedGaussian { radii: f64 },
    /// `(1 - u^2)^2` within one scaled radius and zero beyond it. Compact by
    /// construction rather than by truncation, so nothing is dropped and the
    /// integral is exactly one.
    Quartic,
}

/// A kernel ready to weigh separations: a shape, a bandwidth, and the constant
/// that makes the two integrate to one.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Kernel<const N: usize> {
    shape: KernelShape,
    bandwidth: Bandwidth<N>,
    /// Shape constant divided by the bandwidth volume, folded together at
    /// construction because it multiplies every one of the millions of weights
    /// a field asks for.
    scale: f64,
}

impl<const N: usize> Kernel<N> {

    /// An untruncated Gaussian.
    pub fn gaussian(bandwidth: Bandwidth<N>) -> Self {
        Self::of_shape(KernelShape::Gaussian, bandwidth)
    }

    /// A Gaussian cut off at `radii` bandwidths.
    ///
    /// `None` for a cut at or below zero, or one that is not finite -- the
    /// latter being a way of asking for the untruncated kernel that would
    /// leave the neighbour search believing it had a bounded reach.
    pub fn truncated_gaussian(bandwidth: Bandwidth<N>, radii: f64) -> Option<Self> {

        if !radii.is_finite() || radii <= 0.0 {
            return None;
        }

        Some(Self::of_shape(KernelShape::TruncatedGaussian { radii }, bandwidth))
    }

    /// A quartic, or biweight, kernel.
    pub fn quartic(bandwidth: Bandwidth<N>) -> Self {
        Self::of_shape(KernelShape::Quartic, bandwidth)
    }

    fn of_shape(shape: KernelShape, bandwidth: Bandwidth<N>) -> Self {

        let constant = match shape {
            KernelShape::Gaussian | KernelShape::TruncatedGaussian { .. } => {
                (2.0 * PI).powf(-(N as f64) / 2.0)
            }
            KernelShape::Quartic => quartic_constant(N),
        };

        Self { shape, bandwidth, scale: constant / bandwidth.volume() }
    }

    pub fn shape(&self) -> KernelShape {
        self.shape
    }

    pub fn bandwidth(&self) -> &Bandwidth<N> {
        &self.bandwidth
    }

    /// What one observation at `observed` contributes at `place`.
    ///
    /// Density per unit volume, so summing this over a set of observations
    /// gives a density and not a count. Zero beyond the reach, exactly, for
    /// the shapes that have one.
    pub fn weight(&self, observed: &Point<N>, place: &Point<N>) -> f64 {

        let u = self.bandwidth.scaled_radius(observed, place);

        let profile = match self.shape {
            KernelShape::Gaussian => (-0.5 * u * u).exp(),
            KernelShape::TruncatedGaussian { radii } => {
                if u > radii {
                    return 0.0;
                }
                (-0.5 * u * u).exp()
            }
            KernelShape::Quartic => {
                if u > 1.0 {
                    return 0.0;
                }
                let falloff = 1.0 - u * u;
                falloff * falloff
            }
        };

        self.scale * profile
    }

    /// How far along each axis an observation can still be felt, or `None`
    /// where the kernel never quite stops.
    ///
    /// This is what lets a neighbour search bin the observations and visit
    /// only the bins around a node. `None` is not a failure -- it is the
    /// untruncated Gaussian saying, correctly, that no bin can be skipped.
    pub fn reach(&self) -> Option<[f64; N]> {

        let radii = match self.shape {
            KernelShape::Gaussian => return None,
            KernelShape::TruncatedGaussian { radii } => radii,
            KernelShape::Quartic => 1.0,
        };

        Some(std::array::from_fn(|i| radii * self.bandwidth.extents[i]))
    }

    /// The fraction of the kernel's mass that lies within its reach.
    ///
    /// One for every shape but the truncated Gaussian, which is the only one
    /// that throws anything away. Quoted rather than corrected for: a density
    /// short by a known thousandth is easier to reason about than one silently
    /// rescaled, and a caller who wants the missing mass back can divide.
    pub fn retained_mass(&self) -> f64 {

        let radii = match self.shape {
            KernelShape::TruncatedGaussian { radii } => radii,
            _ => return 1.0,
        };

        // The mass inside a radius is a one-dimensional integral once the
        // angular part is done: the shell area at radius r times the profile,
        // integrated outwards. Simpson over a few thousand steps settles this
        // far below the precision anyone reads it at, and it is computed on
        // demand rather than at construction because most callers never ask.
        let steps = 4096;
        let width = radii / steps as f64;
        let profile = |r: f64| (-0.5 * r * r).exp() * r.powi(N as i32 - 1);

        let mut total = profile(0.0) + profile(radii);
        for step in 1..steps {
            let r = step as f64 * width;
            total += profile(r) * if step % 2 == 0 { 2.0 } else { 4.0 };
        }
        let integral = total * width / 3.0;

        (2.0 * PI).powf(-(N as f64) / 2.0) * unit_sphere_area(N) * integral
    }
}

/// Surface area of the unit sphere in `N` dimensions: `2 pi^(N/2) / Gamma(N/2)`.
///
/// Two points on a line, a circumference of `2 pi` on a plane, `4 pi` in a
/// volume, and it keeps going.
fn unit_sphere_area(n: usize) -> f64 {
    2.0 * PI.powf(n as f64 / 2.0) / half_integer_gamma(n)
}

/// `Gamma(n/2)` for a positive integer `n`.
///
/// Only half-integer arguments arise here, and they have closed forms reached
/// by stepping the recurrence `Gamma(x+1) = x Gamma(x)` up from `Gamma(1/2) =
/// sqrt(pi)` or `Gamma(1) = 1` -- so no general gamma function is needed, and
/// no approximation enters a constant that ought to be exact.
fn half_integer_gamma(n: usize) -> f64 {

    let even = n.is_multiple_of(2);

    let mut value = if even { 1.0 } else { PI.sqrt() };
    let mut argument = if even { 2 } else { 1 };

    while argument + 2 <= n {
        value *= argument as f64 / 2.0;
        argument += 2;
    }

    value
}

/// The constant that makes `(1 - u^2)^2` integrate to one over the unit ball
/// in `N` dimensions.
///
/// The radial integral `∫ (1 - r^2)^2 r^(N-1) dr` from zero to one comes out
/// as `1/N - 2/(N+2) + 1/(N+4)` term by term, and the angular part is the
/// sphere area. In three dimensions that is `1/(4 pi * 8/105)`, or
/// `105/(32 pi)`.
fn quartic_constant(n: usize) -> f64 {

    let dimension = n as f64;
    let radial = 1.0 / dimension - 2.0 / (dimension + 2.0) + 1.0 / (dimension + 4.0);

    1.0 / (unit_sphere_area(n) * radial)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Sphere areas written out rather than computed, so a test of the
    /// constants is not a test against the same gamma routine that produced
    /// them.
    const AREAS: [f64; 3] = [2.0, 2.0 * PI, 4.0 * PI];

    fn area(n: usize) -> f64 {
        AREAS[n - 1]
    }

    #[test]
    fn the_gamma_recurrence_reaches_the_values_it_should() {
        assert!((half_integer_gamma(1) - PI.sqrt()).abs() < 1e-15);
        assert_eq!(half_integer_gamma(2), 1.0);
        assert!((half_integer_gamma(3) - PI.sqrt() / 2.0).abs() < 1e-15);
        assert_eq!(half_integer_gamma(4), 1.0);
        assert!((half_integer_gamma(5) - 3.0 * PI.sqrt() / 4.0).abs() < 1e-15);
        // Gamma(3) = 2! = 2
        assert!((half_integer_gamma(6) - 2.0).abs() < 1e-15);
    }

    #[test]
    fn the_sphere_areas_come_out_right() {
        for n in 1..=3 {
            assert!(
                (unit_sphere_area(n) - area(n)).abs() < 1e-13,
                "unit sphere area in {} dimensions: {} against {}",
                n,
                unit_sphere_area(n),
                area(n)
            );
        }
    }

    #[test]
    fn the_quartic_constant_is_the_textbook_one_in_each_dimension() {
        // 15/16 on a line, 3/pi on a plane, 105/(32 pi) in a volume. The third
        // is the one the C++ this descends from had wrong.
        assert!((quartic_constant(1) - 15.0 / 16.0).abs() < 1e-14);
        assert!((quartic_constant(2) - 3.0 / PI).abs() < 1e-14);
        assert!((quartic_constant(3) - 105.0 / (32.0 * PI)).abs() < 1e-14);
    }

    /// Total mass of an isotropic kernel of unit bandwidth, by radial
    /// quadrature with the sphere area taken from the table above.
    fn total_mass<const N: usize>(kernel: &Kernel<N>, out_to: f64) -> f64 {

        let steps = 20_000;
        let width = out_to / steps as f64;
        let origin = Point::<N>::from([0.0; N]);

        let shell = |r: f64| {
            let mut coords = [0.0; N];
            coords[0] = r;
            kernel.weight(&origin, &Point::from(coords)) * r.powi(N as i32 - 1)
        };

        let mut total = 0.0;
        for step in 0..steps {
            let a = step as f64 * width;
            let b = a + width;
            let m = 0.5 * (a + b);
            total += (shell(a) + 4.0 * shell(m) + shell(b)) * width / 6.0;
        }

        total * area(N)
    }

    #[test]
    fn a_quartic_kernel_integrates_to_one() {
        // In every dimension it is offered in, and this is the property the
        // constants exist for: a density that integrates to the number of
        // observations is a density, and one that does not is a picture.
        let one_d = Kernel::<1>::quartic(Bandwidth::isotropic(1.0).unwrap());
        let two_d = Kernel::<2>::quartic(Bandwidth::isotropic(1.0).unwrap());
        let three_d = Kernel::<3>::quartic(Bandwidth::isotropic(1.0).unwrap());

        assert!((total_mass(&one_d, 1.0) - 1.0).abs() < 1e-9);
        assert!((total_mass(&two_d, 1.0) - 1.0).abs() < 1e-9);
        assert!((total_mass(&three_d, 1.0) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn a_gaussian_kernel_integrates_to_one() {
        let two_d = Kernel::<2>::gaussian(Bandwidth::isotropic(1.0).unwrap());
        let three_d = Kernel::<3>::gaussian(Bandwidth::isotropic(1.0).unwrap());

        // Out to twelve bandwidths, past which the integrand underflows.
        assert!((total_mass(&two_d, 12.0) - 1.0).abs() < 1e-9);
        assert!((total_mass(&three_d, 12.0) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn the_bandwidth_volume_scales_the_density_and_not_the_count() {
        // Doubling the bandwidth in three dimensions divides the peak by
        // eight, the mass staying at one.
        let tight = Kernel::<3>::quartic(Bandwidth::isotropic(1.0).unwrap());
        let loose = Kernel::<3>::quartic(Bandwidth::isotropic(2.0).unwrap());

        let origin = Point::<3>::from([0.0, 0.0, 0.0]);
        let peak_ratio = tight.weight(&origin, &origin) / loose.weight(&origin, &origin);

        assert!((peak_ratio - 8.0).abs() < 1e-12, "peak ratio {}", peak_ratio);
    }

    #[test]
    fn an_anisotropic_bandwidth_is_felt_axis_by_axis() {
        // A kernel reaching ten units east and one down should weigh a point
        // ten units east exactly as it weighs one a single unit down: both sit
        // at one bandwidth, whatever the map says.
        let kernel = Kernel::<3>::quartic(Bandwidth::new([10.0, 10.0, 1.0]).unwrap());
        let node = Point::<3>::from([0.0, 0.0, 0.0]);

        let along_strike = kernel.weight(&Point::from([5.0, 0.0, 0.0]), &node);
        let in_depth = kernel.weight(&Point::from([0.0, 0.0, 0.5]), &node);

        assert!((along_strike - in_depth).abs() < 1e-15);

        // And the isotropic reading of the same separations disagrees, which
        // is the whole point of carrying one extent per axis.
        let isotropic = Kernel::<3>::quartic(Bandwidth::isotropic(10.0).unwrap());
        assert!(
            isotropic.weight(&Point::from([0.0, 0.0, 0.5]), &node)
                > isotropic.weight(&Point::from([5.0, 0.0, 0.0]), &node)
        );
    }

    #[test]
    fn a_quartic_kernel_stops_at_its_bandwidth() {
        let kernel = Kernel::<3>::quartic(Bandwidth::isotropic(2500.0).unwrap());
        let node = Point::<3>::from([0.0, 0.0, 0.0]);

        assert!(kernel.weight(&Point::from([2499.0, 0.0, 0.0]), &node) > 0.0);
        // Exactly zero at the edge, not merely small: the neighbour search is
        // entitled to skip everything beyond the reach without changing an
        // answer.
        assert_eq!(kernel.weight(&Point::from([2500.0, 0.0, 0.0]), &node), 0.0);
        assert_eq!(kernel.weight(&Point::from([2500.1, 0.0, 0.0]), &node), 0.0);
    }

    #[test]
    fn an_untruncated_gaussian_never_reaches_zero() {
        let kernel = Kernel::<3>::gaussian(Bandwidth::isotropic(1.0).unwrap());
        let node = Point::<3>::from([0.0, 0.0, 0.0]);

        assert!(kernel.weight(&Point::from([20.0, 0.0, 0.0]), &node) > 0.0);
        assert!(kernel.reach().is_none());
    }

    #[test]
    fn truncation_bounds_the_reach_and_says_what_it_cost() {
        let kernel =
            Kernel::<3>::truncated_gaussian(Bandwidth::isotropic(1000.0).unwrap(), 3.0).unwrap();

        assert_eq!(kernel.reach(), Some([3000.0, 3000.0, 3000.0]));

        let node = Point::<3>::from([0.0, 0.0, 0.0]);
        assert!(kernel.weight(&Point::from([2999.0, 0.0, 0.0]), &node) > 0.0);
        assert_eq!(kernel.weight(&Point::from([3001.0, 0.0, 0.0]), &node), 0.0);

        // The Maxwell distribution puts 97.07 per cent of a three-dimensional
        // Gaussian within three bandwidths: erf(3/sqrt 2) - sqrt(2/pi) 3 e^-4.5.
        assert!(
            (kernel.retained_mass() - 0.970706_f64).abs() < 1e-5,
            "retained mass {}",
            kernel.retained_mass()
        );

        // Four bandwidths leave a thousandth outside, which is the reason it
        // is the reach worth reaching for.
        let wider =
            Kernel::<3>::truncated_gaussian(Bandwidth::isotropic(1000.0).unwrap(), 4.0).unwrap();
        assert!((wider.retained_mass() - 0.998866_f64).abs() < 1e-5);
    }

    #[test]
    fn the_shapes_that_keep_everything_say_so() {
        assert_eq!(
            Kernel::<3>::quartic(Bandwidth::isotropic(1.0).unwrap()).retained_mass(),
            1.0
        );
        assert_eq!(
            Kernel::<3>::gaussian(Bandwidth::isotropic(1.0).unwrap()).retained_mass(),
            1.0
        );
    }

    #[test]
    fn a_bandwidth_must_be_positive_and_finite_on_every_axis() {
        assert!(Bandwidth::<3>::new([1.0, 1.0, 0.0]).is_none());
        assert!(Bandwidth::<3>::new([1.0, -1.0, 1.0]).is_none());
        assert!(Bandwidth::<3>::new([1.0, f64::NAN, 1.0]).is_none());
        assert!(Bandwidth::<3>::new([f64::INFINITY, 1.0, 1.0]).is_none());
        assert!(Bandwidth::<3>::isotropic(0.0).is_none());
        assert!(Bandwidth::<3>::new([1.0, 2.0, 3.0]).is_some());
    }

    #[test]
    fn a_zero_dimensional_bandwidth_is_refused() {
        // There is no volume to spread density over, and the constants would
        // divide by a gamma function that is not defined there.
        assert!(Bandwidth::<0>::new([]).is_none());
    }

    #[test]
    fn a_truncation_must_be_a_positive_distance() {
        let bandwidth = Bandwidth::<3>::isotropic(1.0).unwrap();

        assert!(Kernel::truncated_gaussian(bandwidth, 0.0).is_none());
        assert!(Kernel::truncated_gaussian(bandwidth, -1.0).is_none());
        // An infinite cut is a way of asking for the plain Gaussian that would
        // leave `reach` reporting a bound no bin could be tested against.
        assert!(Kernel::truncated_gaussian(bandwidth, f64::INFINITY).is_none());
        assert!(Kernel::truncated_gaussian(bandwidth, 4.0).is_some());
    }

    #[test]
    fn a_truncated_gaussian_agrees_with_the_plain_one_inside_the_cut() {
        let bandwidth = Bandwidth::<3>::isotropic(500.0).unwrap();
        let plain = Kernel::gaussian(bandwidth);
        let cut = Kernel::truncated_gaussian(bandwidth, 4.0).unwrap();

        let node = Point::<3>::from([0.0, 0.0, 0.0]);
        let near = Point::<3>::from([300.0, 200.0, -100.0]);

        // Truncation removes a tail; it does not rescale what is left.
        assert_eq!(plain.weight(&near, &node), cut.weight(&near, &node));
    }
}
