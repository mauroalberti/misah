// private sub-module defined in other files
mod axis;
mod geolplane;

// exports identifiers from private sub-modules in the current module namespace
pub use self::orientations::space3d::axis::Axis;
pub use self::orientations::space3d::geolplane::GeolPlane;


