
use crate::orientation::direction::Direction3D;

pub enum SlipSense {
    Up,
    Down,
    Left,
    Right,
    Unknown,
}

pub struct Slickenline {
    pub lineation: Direction3D,
    pub sense: Option<SlipSense>,
}
