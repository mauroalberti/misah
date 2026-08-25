
use crate::geoprofile::records::{
    intersection::IntersectionRecord,
    profile::ProfileRecord,
    profile_sample::ProfileSampleRecord,
    projection::{ProjectedAttitudeRecord, ProjectedPointRecord},
    result_set::ResultSetRecord,
    source::SourceRecord,
};

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
}
