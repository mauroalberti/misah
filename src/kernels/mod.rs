//! Numerical kernels.
//!
//! Kernels are deliberately free of pyo3 and of any I/O: they take plain slices
//! and return plain vectors, so they can be unit-tested with `cargo test` alone
//! and bound to Python separately, in `py`.

pub mod plane_dem;

#[cfg(feature = "extension-module")]
pub mod py;

pub use self::plane_dem::{intersect_plane_dem, plane_normal, GeoTransform, Intersection};
