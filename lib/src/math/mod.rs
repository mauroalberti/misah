
pub mod vector;
pub mod versor;
pub mod error;

pub use error::VectorError;
pub use vector::Vector;
pub use versor::Versor;

pub type Vector2D = Vector<2>;
pub type Vector3D = Vector<3>;

pub type Versor2D = Versor<2>;
pub type Versor3D = Versor<3>;

