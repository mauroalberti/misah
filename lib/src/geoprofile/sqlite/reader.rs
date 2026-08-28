
use crate::geometry::point::Point3D;

use rusqlite::{Connection, OpenFlags};
use std::path::Path;

use crate::geoprofile::{
    dataset::GeoProfileDataset,
    error::GeoProfileError,
    records::{
        intersection::IntersectionRecord,
        profile::ProfileRecord,
        profile_sample::ProfileSampleRecord,
        result_set::ResultSetRecord,
        source::SourceRecord,
    },
    sqlite::schema::validate_geoprofile_schema,
};

use crate::geoprofile::records::projection::{
    DownSense,
    ProjectedAttitudeData,
    ProjectedAttitudeRecord,
    ProjectedBase,
    ProjectedPointData,
    ProjectedPointRecord,
};

pub struct SqliteGeoProfileReader {
    conn: Connection,
}

impl SqliteGeoProfileReader {

    pub fn open(path: impl AsRef<Path>) -> Result<Self, GeoProfileError> {
        let conn = Connection::open_with_flags(
            path,
            OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_URI,
        )?;
        Ok(Self { conn })
    }

    pub fn validate_schema(&self) -> Result<(), GeoProfileError> {
        validate_geoprofile_schema(&self.conn)
    }

    fn make_src_point(
        x: Option<f64>,
        y: Option<f64>,
        z: Option<f64>,
    ) -> Option<Point3D> {

        match (x, y, z) {
            (Some(x), Some(y), Some(z)) => Some(Point3D::from([x, y, z])),
            _ => None,
        }
    }

    pub fn read_result_sets(&self) -> Result<Vec<ResultSetRecord>, GeoProfileError> {
        let mut stmt = self.conn.prepare(
            r#"
            SELECT
                result_set_id,
                name,
                description,
                created_utc,
                schema_version,
                distance_units,
                extra_json
            FROM gp_result_sets
            ORDER BY result_set_id
            "#,
        )?;

        let rows = stmt.query_map([], |row| {
            Ok(ResultSetRecord {
                result_set_id: row.get(0)?,
                name: row.get(1)?,
                // DEFAULT '' in the export schema, so the format's own way of
                // saying "nothing here" is the empty string. NULL is still
                // storable, though, and reading it straight into a String would
                // fail the whole read on one stray value; it maps to the empty
                // string the column would have held anyway.
                description: row.get::<_, Option<String>>(2)?.unwrap_or_default(),
                created_utc: row.get(3)?,
                schema_version: row.get(4)?,
                distance_units: row.get::<_, Option<String>>(5)?.unwrap_or_default(),
                extra_json: row.get(6)?,
            })
        })?;

        rows.collect::<Result<Vec<_>, _>>().map_err(Into::into)

    }

    pub fn read_profiles(&self) -> Result<Vec<ProfileRecord>, GeoProfileError> {
        let mut stmt = self.conn.prepare(
            r#"
            SELECT
                profile_id,
                result_set_id,
                source_profile_fid,
                profile_name,
                wkt_crs,
                s_max,
                z_min,
                z_max,
                extra_json
            FROM gp_profiles
            ORDER BY profile_id
            "#,
        )?;

        let rows = stmt.query_map(
            [],
            |row| {
                Ok(ProfileRecord {
                    profile_id: row.get(0)?,
                    result_set_id: row.get(1)?,
                    source_profile_fid: row.get(2)?,
                    // DEFAULT '' in the export schema; see read_result_sets.
                    profile_name: row.get::<_, Option<String>>(3)?.unwrap_or_default(),
                    wkt_crs: row.get(4)?,
                    s_max: row.get(5)?,
                    z_min: row.get(6)?,
                    z_max: row.get(7)?,
                    extra_json: row.get(8)?,
                    }
                )
            }
        )?;

        rows.collect::<Result<Vec<_>, _>>().map_err(Into::into)

    }

    pub fn read_profile_samples(&self) -> Result<Vec<ProfileSampleRecord>, GeoProfileError> {
        let mut stmt = self.conn.prepare(
            r#"
            SELECT
                sample_id,
                profile_id,
                idx,
                s,
                z,
                x,
                y,
                extra_json
            FROM gp_profile_samples
            ORDER BY profile_id, idx
            "#,
        )?;

        let rows = stmt.query_map(
            [],
            |row| {
                Ok(
                    ProfileSampleRecord {
                        sample_id: row.get(0)?,
                        profile_id: row.get(1)?,
                        idx: row.get(2)?,
                        s: row.get(3)?,
                        z: row.get(4)?,
                        x: row.get(5)?,
                        y: row.get(6)?,
                        extra_json: row.get(7)?,
                    }
                )
            }
        )?;

        rows.collect::<Result<Vec<_>, _>>().map_err(Into::into)
    }

    pub fn read_projected_points(&self) -> Result<Vec<ProjectedPointRecord>, GeoProfileError> {

        let mut stmt = self.conn.prepare(
            r#"
            SELECT
                rec_id,
                profile_id,
                label,
                s,
                z,
                dist_to_profile,
                src_x,
                src_y,
                src_z,
                src_fid,
                extra_json
            FROM gp_projected_points
            ORDER BY profile_id, s, rec_id
            "#,
        )?;

        let rows = stmt.query_map(
            [],
            |row| {
                let src_x: Option<f64> = row.get(6)?;
                let src_y: Option<f64> = row.get(7)?;
                let src_z: Option<f64> = row.get(8)?;

                Ok( ProjectedPointRecord {
                    base: ProjectedBase {
                        rec_id: row.get(0)?,
                        profile_id: row.get(1)?,
                        category: row.get(2)?,
                        s: row.get(3)?,
                        z: row.get(4)?,
                        dist_to_profile: row.get(5)?,
                        src_point: Self::make_src_point(src_x, src_y, src_z),
                        src_fid: row.get(9)?,
                        extra_json: row.get(10)?,
                    },
                    data: ProjectedPointData,
                })
            }
        )?;

        rows.collect::<Result<Vec<_>, _>>().map_err(Into::into)

    }

    pub fn read_projected_attitudes(
        &self,
    ) -> Result<Vec<ProjectedAttitudeRecord>, GeoProfileError> {

        let mut stmt = self.conn.prepare(
            r#"SELECT
                rec_id,
                profile_id,
                label,
                s,
                z,
                slope_degr,
                down_sense,
                src_dip_dir,
                src_dip_ang,
                dist_to_profile,
                src_x,
                src_y,
                src_z,
                src_fid,
                extra_json
            FROM gp_projected_attitudes
            ORDER BY profile_id, s, rec_id
            "#,
        )?;

        let rows = stmt.query_map([], |row| {

            let down_sense_raw: String = row.get(6)?;

            let down_sense = match down_sense_raw.to_lowercase().as_str() {
                "left" => DownSense::Left,
                "right" => DownSense::Right,
                other => {
                    return Err(rusqlite::Error::InvalidColumnType(
                        6,
                        format!("invalid down_sense value: {other}"),
                        rusqlite::types::Type::Text,
                    ));
                }
            };

            let src_x: Option<f64> = row.get(10)?;
            let src_y: Option<f64> = row.get(11)?;
            let src_z: Option<f64> = row.get(12)?;

            Ok(ProjectedAttitudeRecord {
                base: ProjectedBase {
                    rec_id: row.get(0)?,
                    profile_id: row.get(1)?,
                    category: row.get(2)?,
                    s: row.get(3)?,
                    z: row.get(4)?,
                    dist_to_profile: row.get(9)?,
                    src_point: Self::make_src_point(src_x, src_y, src_z),
                    src_fid: row.get(13)?,
                    extra_json: row.get(14)?,
                },
                data: ProjectedAttitudeData {
                    slope_degr: row.get(5)?,
                    down_sense,
                    src_dip_dir: row.get(7)?,
                    src_dip_ang: row.get(8)?,
                },
            })
        })?;

        rows.collect::<Result<Vec<_>, _>>().map_err(Into::into)

    }

    pub fn read_line_intersections(&self) -> Result<Vec<IntersectionRecord>, GeoProfileError> {

        let mut stmt = self.conn.prepare(
            r#"
            SELECT
                rec_id,
                profile_id,
                feat_category,
                s,
                extra_json
            FROM gp_intersected_lines
            ORDER BY profile_id, s, rec_id
            "#,
        )?;

        let rows = stmt.query_map(
            [],
            |row| {
                // A line crosses the section at one point, so the span is
                // degenerate and both ends carry that distance. The exception is
                // a segment lying along the section trace, which is crossed over
                // a stretch -- but the export gives no way to recover it: the
                // exporter writes the two ends of such a stretch as two ordinary
                // rows, indistinguishable from two separate crossings of the
                // same category. Reading a row as a span from itself is
                // therefore all this table supports.
                let s: f64 = row.get(3)?;

                Ok(
                    IntersectionRecord {
                        rec_id: row.get(0)?,
                        profile_id: row.get(1)?,
                        category: row.get(2)?,
                        s_from: s,
                        s_to: s,
                        extra_json: row.get(4)?,
                    }
                )
            }
        )?;

        rows.collect::<Result<Vec<_>, _>>().map_err(Into::into)
    }

    pub fn read_polygon_intersections(&self) -> Result<Vec<IntersectionRecord>, GeoProfileError> {

        let mut stmt = self.conn.prepare(
            r#"
            SELECT
                rec_id,
                profile_id,
                unit_name,
                s_from,
                s_to,
                extra_json
            FROM gp_intersected_polygons
            ORDER BY profile_id, s_from, s_to, rec_id
            "#,
        )?;

        let rows = stmt.query_map(
            [],
            |row| {
                Ok(
                    IntersectionRecord {
                        rec_id: row.get(0)?,
                        profile_id: row.get(1)?,
                        category: row.get(2)?,
                        s_from: row.get(3)?,
                        s_to: row.get(4)?,
                        extra_json: row.get(5)?,
                    }
                )
            }
        )?;

        rows.collect::<Result<Vec<_>, _>>().map_err(Into::into)

    }

    pub fn read_sources(&self) -> Result<Vec<SourceRecord>, GeoProfileError> {
        let mut stmt = self.conn.prepare(
            r#"
            SELECT
                source_id,
                result_set_id,
                source_type,
                source_name,
                source_uri,
                wkt_crs,
                parameters_json
            FROM gp_sources
            ORDER BY source_id
           "#,
        )?;

        let rows = stmt.query_map(
           [],
           |row| {
               Ok(
                   SourceRecord {
                    source_id: row.get(0)?,
                    result_set_id: row.get(1)?,
                    source_type: row.get(2)?,
                    // Both DEFAULT '' in the export schema; see read_result_sets.
                    source_name: row.get::<_, Option<String>>(3)?.unwrap_or_default(),
                    source_uri: row.get::<_, Option<String>>(4)?.unwrap_or_default(),
                    wkt_crs: row.get(5)?,
                    parameters_json: row.get(6)?,
                })
            })?;

        rows.collect::<Result<Vec<_>, _>>().map_err(Into::into)

    }

    pub fn read_all(&self) -> Result<GeoProfileDataset, GeoProfileError> {

        Ok(GeoProfileDataset {
            result_sets: self.read_result_sets()?,
            profiles: self.read_profiles()?,
            profile_samples: self.read_profile_samples()?,
            polygon_intersections: self.read_polygon_intersections()?,
            line_intersections: self.read_line_intersections()?,
            projected_attitudes: self.read_projected_attitudes()?,
            projected_points: self.read_projected_points()?,
            sources: self.read_sources()?,
        })
    }

}

