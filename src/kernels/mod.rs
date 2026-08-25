//! Numerical kernels.
//!
//! Kernels are deliberately free of pyo3 and of any I/O: they take plain slices
//! and return plain vectors, so they can be unit-tested with `cargo test` alone
//! and bound to Python separately.

pub mod plane_dem;

pub use self::plane_dem::{intersect_plane_dem, plane_normal, GeoTransform, Intersection};
