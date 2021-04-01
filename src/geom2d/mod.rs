// private sub-module defined in other files
mod point2d;
mod segment2d;

// exports identifiers from private sub-modules in the current module namespace
pub use self::point2d::Point2D;
pub use self::segment2d::Segment2D;



