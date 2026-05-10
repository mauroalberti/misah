
use super::geol_axis::{GeologicalAxis};

pub struct GeologicalPlane {
    pub azimuth: f64,
    pub dip_angle: f64
}

impl GeologicalPlane {

    fn new(az: f64, dip: f64) -> Self { GeologicalPlane{ azimuth: az, dip_angle: dip }}

    pub fn normal_axis(&self) -> GeologicalAxis {

        GeologicalAxis{
            trend: (self.azimuth + 180.0) % 360.0,
            plunge: 90.0 - self.dip_angle

        }
    }

}
