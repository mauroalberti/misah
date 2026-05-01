
#[derive(Debug, Clone)]
pub struct ProfileRecord {
    pub profile_id: i64,
    pub result_set_id: i64,
    pub source_profile_fid: Option<i64>,
    pub profile_name: String,
    pub wkt_crs: Option<String>,
    pub s_max: Option<f64>,
    pub z_min: Option<f64>,
    pub z_max: Option<f64>,
    pub extra_json: Option<String>,
}
