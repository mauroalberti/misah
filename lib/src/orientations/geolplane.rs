
use super::axis::{Axis};

pub struct GeolPlane {
    pub azimuth: f64,
    pub dip_angle: f64
}

impl GeolPlane {

    fn new(az: f64, dip: f64) -> Self { GeolPlane{ azimuth: az, dip_angle: dip }}

    pub fn normal_axis(&self) -> Axis {
        Axis{trend: (self.azimuth + 180.0) % 360.0, plunge: 90.0 - self.dip_angle}
    }

}