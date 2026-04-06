
use std::path::Path;
use rusqlite::Connection;

use crate::error::GxsurfError;
use crate::profiles::attitudes::Attitude;

pub struct SqliteGeoprofileReader {
    conn: Connection,
}

impl SqliteGeoprofileReader {

    pub fn open(path: impl AsRef<Path>) -> Result<Self, GxsurfError> {
        let conn = Connection::open(path)?;
        Ok(Self { conn })
    }

    pub fn validate_schema(&self) -> Result<(), GxsurfError> {
        // chiamate a schema::require_table / require_columns
        Ok(())
    }

    pub fn read_attitudes(&self) -> Result<Vec<Attitude>, GxsurfError> {

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

    pub fn read_profiles(&self) -> Result<Vec<Profile>, GxsurfError> {

    }


    pub fn read_plane_traces(&self) -> Result<Vec<PlaneTrace>, GxsurfError>;
    pub fn read_polygon_intersections(&self) -> Result<Vec<PolygonIntersection>, GxsurfError>;

    pub fn read_all(path: impl AsRef<Path>) -> Result<GeoprofileDataset, GxsurfError>;

}
