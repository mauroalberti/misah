
#[derive(Debug, Clone)]
pub struct SourceRecord {
    pub source_id: i64,
    pub result_set_id: i64,
    pub source_type: String,
    pub source_name: String,
    pub source_uri: String,
    pub wkt_crs: Option<String>,
    pub parameters_json: Option<String>,
}
