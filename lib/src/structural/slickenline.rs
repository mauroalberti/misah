
use crate::orientation::direction::Direction3D;

#[derive(Debug, Clone)]
pub enum SlipSense {
    Up,
    Down,
    Left,
    Right,
    Unknown,
}

#[derive(Debug, Clone)]
pub struct Slickenline {
    pub lineation: Direction3D,
    pub sense: Option<SlipSense>,
}
