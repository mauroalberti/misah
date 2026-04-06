
#[derive(Debug, Clone)]
pub struct Attitude {
    pub fid: i64,
    pub profile_id: i64,
    pub s: f64,
    pub z: Option<f64>,
    pub slope_degr: f64,
    pub down_sense: Option<String>,
    pub dist: Option<f64>,
    pub src_dip_dir: Option<f64>,
    pub src_dip_ang: Option<f64>,
}
