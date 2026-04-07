
#[derive(Debug, Clone)]
pub struct LineIntersectionRecord {
    pub rec_id: i64,
    pub profile_id: i64,
    pub feat_category: String,
    pub s: f64,
    pub extra_json: Option<String>,
}
