
use std::path::Path;
use rusqlite::Connection;

use crate::error::GeoProfileError;
use crate::profiles::attitudes::Attitude;

pub struct SqliteGeoProfileReader {
    conn: Connection,
}

impl SqliteGeoProfileReader {

    pub fn open(path: impl AsRef<Path>) -> Result<Self, GeoProfileError> {
        let conn = Connection::open(path)?;
        Ok(Self { conn })
    }

    pub fn validate_schema(&self) -> Result<(), GeoProfileError> {
        // chiamate a schema::require_table / require_columns
        Ok(())
    }

    pub fn read_result_sets(&self) -> Result<Vec<ResultSetRecord>, GeoProfileError>;

    pub fn read_profiles(&self) -> Result<Vec<ProfileRecord>, GeoProfileError> {

    }

    pub fn read_profile_samples(&self) -> Result<Vec<ProfileSampleRecord>, GeoProfileError>;

    pub fn read_projected_points(&self) -> Result<Vec<ProjectedPointRecord>, GeoProfileError>

    pub fn read_projected_attitudes(&self) -> Result<Vec<ProjectedAttitudeRecord>, GeoProfileError> {

        let mut stmt = self.conn.prepare(
            "SELECT fid, profile_id, s, z, slope_degr, down_sense, dist, src_dip_dir, src_dip_ang FROM attitudes"
        )?;

        let rows = stmt.query_map([], |row| {
            Ok(Attitude {
                fid: row.get(0)?,
                profile_id: row.get(1)?,
                s: row.get(2)?,
                z: row.get(3)?,
                slope_degr: row.get(4)?,
                down_sense: row.get(5)?,
                dist: row.get(6)?,
                src_dip_dir: row.get(7)?,
                src_dip_ang: row.get(8)?,
            })
        })?;

        let mut out = Vec::new();
        for item in rows {
            out.push(item?);
        }
        Ok(out)
    }

    pub fn read_plane_traces(&self) -> Result<Vec<PlaneTrace>, GeoProfileError>;

    pub fn read_line_intersections(&self) -> Result<Vec<LineIntersectionRecord>, GeoProfileError>;

    pub fn read_polygon_intersections(&self) -> Result<Vec<PolygonIntersectionRecord>, GeoProfileError>;

    pub fn read_all(path: impl AsRef<Path>) -> Result<GeoProfileDataset, GeoProfileError>;

}
