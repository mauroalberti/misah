
#[derive(Debug, Clone)]
pub struct IntersectionRecord {
    pub rec_id: i64,
    pub profile_id: i64,
    pub category: Option<String>,
    pub s_from: f64,
    pub s_to: f64,
    pub extra_json: Option<String>,
}
