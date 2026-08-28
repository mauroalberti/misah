
use crate::geoprofile::records::{
    graphical::{GraphicalParamsRecord, SourceCategoryRecord},
    intersection::IntersectionRecord,
    profile::ProfileRecord,
    profile_sample::ProfileSampleRecord,
    projection::{ProjectedAttitudeRecord, ProjectedFocalMechanismRecord, ProjectedPointRecord},
    result_set::ResultSetRecord,
    source::SourceRecord,
    vertex::ProfileVertexRecord,
};

/// Everything one GeoProfiler export holds, as far as this reader goes.
///
/// The last four are absent from older exports -- qgSurf added them at schema
/// v2 to v4 -- and come back empty from those, rather than making the file
/// unreadable.
#[derive(Debug, Clone, Default)]
pub struct GeoProfileDataset {
    pub profiles: Vec<ProfileRecord>,
    pub profile_samples: Vec<ProfileSampleRecord>,
    pub projected_points: Vec<ProjectedPointRecord>,
    pub projected_attitudes: Vec<ProjectedAttitudeRecord>,
    pub line_intersections: Vec<IntersectionRecord>,
    pub polygon_intersections: Vec<IntersectionRecord>,
    pub result_sets: Vec<ResultSetRecord>,
    pub sources: Vec<SourceRecord>,
    pub projected_focal_mechanisms: Vec<ProjectedFocalMechanismRecord>,
    pub profile_vertices: Vec<ProfileVertexRecord>,
    pub graphical_params: Vec<GraphicalParamsRecord>,
    pub source_categories: Vec<SourceCategoryRecord>,
}
