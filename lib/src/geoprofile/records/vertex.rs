
/// Where a profile starts, turns, and ends.
///
/// A profile drawn as a polyline has a vertex at each change of direction;
/// these are those, with their position along the section as well as on the
/// ground, so a section can be annotated with the bends of its own trace.
#[derive(Debug, Clone)]
pub struct ProfileVertexRecord {
    pub vertex_id: i64,
    pub profile_id: i64,
    /// Position along the profile's own vertex sequence, from zero.
    pub vertex_ndx: i64,
    pub kind: VertexKind,
    /// Distance from the profile start.
    pub s: f64,
    pub x: f64,
    pub y: f64,
    /// Geographic coordinates, written beside the projected ones when the
    /// export could compute them.
    pub lon: Option<f64>,
    pub lat: Option<f64>,
    pub z: Option<f64>,
    pub extra_json: Option<String>,
}

/// What a vertex is to its profile. The export constrains the column to these
/// three, so anything else is a file that did not come from GeoProfiler.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VertexKind {
    Start,
    Break,
    End,
}
