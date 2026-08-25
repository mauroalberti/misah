
#[derive(Debug, Clone)]
pub struct ResultSetRecord {
    pub result_set_id: i64,
    pub name: String,
    pub description: String,
    pub created_utc: String,
    pub schema_version: i64,
    pub distance_units: String,
    pub extra_json: Option<String>,
}
