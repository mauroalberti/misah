

#[derive(Debug, Clone)]
pub struct ProjectedAttitudeRecord {
    pub rec_id: i64,
    pub profile_id: i64,
    pub label: String,
    pub s: f64,
    pub z: f64,
    pub slope_degr: f64,
    pub down_sense: String,
    pub src_dip_dir: f64>,
    pub src_dip_ang: f64,
    pub dist_to_profile: Option<f64>,
    pub src_x: Option<f64>,
    pub src_y: Option<f64>,
    pub src_z: Option<f64>,
    pub src_fid: Option<i64>,
    pub extra_json: Option<String>,
}
