
pub mod constants;
pub mod eigen;
pub mod quaternion;
pub mod vector;
pub mod versor;
pub mod error;

pub use constants::EPSILON;
pub use error::AlgebraError;
pub use quaternion::Quaternion;
pub use vector::Vector;
pub use versor::Versor;


