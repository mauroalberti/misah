
use crate::geoprofile::records::{
    line_intersection::LineIntersectionRecord,
    polygon_intersection::PolygonIntersectionRecord,
    profile::ProfileRecord,
    profile_sample::ProfileSampleRecord,
    projected_attitude::ProjectedAttitudeRecord,
    projected_point::ProjectedPointRecord,
    result_set::ResultSetRecord,
    source::SourceRecord,
}

#[derive(Debug, Clone, Default)]
pub struct GeoProfileDataset {
    pub result_sets: Vec<ResultSetRecord>,
    pub profiles: Vec<ProfileRecord>,
    pub profile_samples: Vec<ProfileSampleRecord>,
    pub projected_points: Vec<ProjectedPointRecord>,
    pub projected_attitudes: Vec<ProjectedAttitudeRecord>,
    pub line_intersections: Vec<LineIntersectionRecord>,
    pub polygon_intersections: Vec<PolygonIntersectionRecord>,
    pub sources: Vec<SourceRecord>,
}
