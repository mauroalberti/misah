// private sub-module defined in other files
mod point2d;
mod segment2d;

// exports identifiers from private sub-modules in the current module namespace
pub use self::geometry::space2d::point2d::Point2D;
pub use self::geometry::space2d::segment2d::Segment2D;



