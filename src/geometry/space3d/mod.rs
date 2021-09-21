// private sub-module defined in other files
mod point3d;
mod segment3d;

// exports identifiers from private sub-modules in the current module namespace
pub use self::geometry::space3d::point3d::Point3D;
pub use self::geometry::space3d::segment3d::Segment3D;


