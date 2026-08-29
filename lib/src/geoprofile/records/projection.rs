
use crate::geometry::point::Point3D;

#[derive(Debug, Clone)]
pub struct ProjectedRecord<T> {
    pub base: ProjectedBase,
    pub data: T,
}

#[derive(Debug, Clone)]
pub struct ProjectedBase {
    pub rec_id: i64,
    pub profile_id: i64,
    pub category: Option<String>,
    pub s: f64,
    pub z: f64,
    pub dist_to_profile: Option<f64>,
    pub src_point: Option<Point3D>,
    pub src_fid: Option<i64>,
    pub extra_json: Option<String>,
}

#[derive(Debug, Clone)]
pub struct ProjectedPointData;

#[derive(Debug, Clone)]
pub struct ProjectedAttitudeData {
    pub slope_degr: f64,
    pub down_sense: DownSense,
    pub src_dip_dir: Option<f64>,
    pub src_dip_ang: Option<f64>,
}

#[derive(Debug, Clone)]
pub enum DownSense {
    Left,
    Right,
}

#[derive(Debug, Clone)]
pub struct ProjectedFocalMechanismData {
    pub strike: f64,
    pub dip: f64,
    pub rake: f64,
    /// Azimuth of the profile where the mechanism was projected, which is what
    /// decides how the mechanism appears on that section rather than another.
    pub profile_azimuth: f64,
}

pub type ProjectedPointRecord = ProjectedRecord<ProjectedPointData>;

pub type ProjectedAttitudeRecord = ProjectedRecord<ProjectedAttitudeData>;

pub type ProjectedFocalMechanismRecord = ProjectedRecord<ProjectedFocalMechanismData>;
