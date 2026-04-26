
use crate::geoprofile::records::{
    profile::ProfileRecord,
    profile_sample::ProfileSampleRecord,
    projected_point::ProjectedPointRecord,
    projected_attitude::ProjectedAttitudeRecord,
    line_intersection::LineIntersectionRecord,
    polygon_intersection::PolygonIntersectionRecord,
    result_set::ResultSetRecord,
    source::SourceRecord,
};

#[derive(Debug, Clone, Default)]
pub struct GeoProfileDataset {
    pub profiles: Vec<ProfileRecord>,
    pub profile_samples: Vec<ProfileSampleRecord>,
    pub projected_points: Vec<ProjectedPointRecord>,
    pub projected_attitudes: Vec<ProjectedAttitudeRecord>,
    pub line_intersections: Vec<LineIntersectionRecord>,
    pub polygon_intersections: Vec<PolygonIntersectionRecord>,
    pub result_sets: Vec<ResultSetRecord>,
    pub sources: Vec<SourceRecord>,
}
