// private sub-module defined in other files
mod point3d;
mod segment3d;

// exports identifiers from private sub-modules in the current module namespace
pub use self::point3d::Point3D;
pub use self::segment3d::Segment3D;


