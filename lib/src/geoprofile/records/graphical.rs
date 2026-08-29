
/// The plot styles in force when a result set was exported.
///
/// Kept as the JSON the exporter wrote: these are qgSurf's own presentation
/// settings, and a kernel has no business giving them a schema of its own. They
/// travel with the data so a section can be redrawn as it was.
#[derive(Debug, Clone)]
pub struct GraphicalParamsRecord {
    pub params_id: i64,
    pub result_set_id: i64,
    /// The profile these apply to, or `None` when they apply to the whole
    /// result set.
    pub profile_id: Option<i64>,
    pub params_json: String,
}

/// One class of a source layer's legend, with the colour it was drawn in.
///
/// Exported whole, rather than only for the categories that happen to occur in
/// a section, so that a redrawn section keeps the classification and the
/// ordering of the map it came from.
#[derive(Debug, Clone)]
pub struct SourceCategoryRecord {
    pub category_id: i64,
    pub result_set_id: i64,
    /// Which family of data the class belongs to: `polygon_intersections` or
    /// `line_intersections`.
    pub data_kind: String,
    /// The position the class holds in the source legend.
    pub position: i64,
    pub category: String,
    /// `#rrggbb`, when the source symbology carried one.
    pub color: Option<String>,
}
