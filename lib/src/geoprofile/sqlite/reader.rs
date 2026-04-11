
use std::path::Path;
use rusqlite::{Connection, OpenFlags};

use crate::geoprofile::{
    dataset::GeoProfileDataset,
    error::GeoProfileError,
    records::{
        line_intersection::LineIntersectionRecord,
        polygon_intersection::PolygonIntersectionRecord,
        profile::ProfileRecord,
        profile_sample::ProfileSampleRecord,
        projected_attitude::ProjectedAttitudeRecord,
        projected_point::ProjectedPointRecord,
        result_set::ResultSetRecord,
        source::SourceRecord,
    },
    sqlite::schema::validate_geoprofile_schema,
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
                description: row.get(2)?,
                created_utc: row.get(3)?,
                schema_version: row.get(4)?,
                distance_units: row.get(5)?,
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
                    profile_name: row.get(3)?,
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
            "#
        )?;

        let rows = stmt.query_map(
            [],
            |row| {
                Ok(
                    ProjectedPointRecord {
                        rec_id: row.get(0)?,
                        profile_id: row.get(1)?,
                        label: row.get(2)?,
                        s: row.get(3)?,
                        z: row.get(4)?,
                        dist_to_profile: row.get(5)?,
                        src_x: row.get(6)?,
                        src_y: row.get(7)?,
                        src_z: row.get(8)?,
                        src_fid: row.get(9)?,
                        extra_json: row.get(10)?,
                    }
                )
            }

        )?;

        rows.collect::<Result<Vec<_>, _>>().map_err(Into::into)

    }

    pub fn read_projected_attitudes(&self) -> Result<Vec<ProjectedAttitudeRecord>, GeoProfileError> {

        let mut stmt = self.conn.prepare(
            r#"SELECT
                fid,
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
            Ok(ProjectedAttitudeRecord {
                rec_id: row.get(0)?,
                profile_id: row.get(1)?,
                label: row.get(2)?,
                s: row.get(3)?,
                z: row.get(4)?,
                slope_degr: row.get(5)?,
                down_sense: row.get(6)?,
                src_dip_dir: row.get(7)?,
                src_dip_ang: row.get(8)?,
                dist_to_profile: row.get(9)?,
                src_x: row.get(10)?,
                src_y: row.get(11)?,
                src_z: row.get(12)?,
                src_fid: row.get(13)?,
                extra_json: row.get(14)?,
            })
        })?;

        rows.collect::<Result<Vec<_>, _>>().map_err(Into::into)

    }

    pub fn read_line_intersections(&self) -> Result<Vec<LineIntersectionRecord>, GeoProfileError> {

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
                Ok(
                    LineIntersectionRecord {
                        rec_id: row.get(0)?,
                        profile_id: row.get(1)?,
                        feat_category: row.get(2)?,
                        s: row.get(3)?,
                        extra_json: row.get(4)?,
                    }
                )
            }
        )?;

        rows.collect::<Result<Vec<_>, _>>().map_err(Into::into)
    }

    pub fn read_polygon_intersections(&self) -> Result<Vec<PolygonIntersectionRecord>, GeoProfileError> {

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
            "#
        )?;

        let rows = stmt.query_map(
            [],
            |row| {
                Ok(
                    PolygonIntersectionRecord {
                        rec_id: row.get(0)?,
                        profile_id: row.get(1)?,
                        unit_name: row.get(2)?,
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
           "#
        )?;

        let rows = stmt.query_map(
           [],
           |row| {
               Ok(
                   SourceRecord {
                    source_id: row.get(0)?,
                    result_set_id: row.get(1)?,
                    source_type: row.get(2)?,
                    source_name: row.get(3)?,
                    source_uri: row.get(4)?,
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

