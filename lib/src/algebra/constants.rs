
/// Geometric tolerance shared across the algebra and geometry modules.
///
/// Matches the working precision already assumed elsewhere in the crate (the
/// intersection tests compare against 1e-9): at the metre scale of a projected
/// CRS, a difference below this is rounding noise, not a real geometric
/// distinction. It is an absolute threshold, so it is not appropriate for data
/// at a very different physical scale.
pub const EPSILON: f64 = 1e-9;
