
use crate::features::points::Point3D;


pub struct ProjectedRecord<T> {
    pub base: ProjectedBase,
    pub data: T,
}

pub struct ProjectedBase {
    pub rec_id: i64,
    pub profile_id: i64,
    pub category: Optional<String>,
    pub s: f64,
    pub z: f64,
    pub dist_to_profile: Option<f64>,
    pub src_point: Option<Point3D>,
    pub src_fid: Option<i64>,
    pub extra_json: Option<String>,
}

pub struct ProjectedPointData;

pub struct ProjectedAttitudeData {
    pub slope_degr: f64,
    pub down_sense: DownSense,
    pub src_dip_dir: Optional<f64>,
    pub src_dip_ang: Optional<f64>,
}

pub enum DownSense {
    Left,
    Right,
}

pub type ProjectedPointRecord = ProjectedRecord<ProjectedPointData>;

pub type ProjectedAttitudeRecord = ProjectedRecord<ProjectedAttitudeData>;
