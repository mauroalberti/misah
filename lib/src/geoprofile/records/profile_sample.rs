
#[derive(Debug, Clone)]
pub struct ProfileSampleRecord {
    pub sample_id: i64,
    pub profile_id: i64,
    pub idx: i64,
    pub s: f64,
    pub z: f64,
    pub x: Option<f64>,
    pub y: Option<f64>,
    pub extra_json: Option<String>,
}
