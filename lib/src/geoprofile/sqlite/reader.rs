
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
    sqlite::schema::{table_columns, table_exists, validate_geoprofile_schema},
};

use crate::geoprofile::records::graphical::{GraphicalParamsRecord, SourceCategoryRecord};
use crate::geoprofile::records::vertex::{ProfileVertexRecord, VertexKind};

use crate::geoprofile::records::projection::{
    DownSense,
    ProjectedAttitudeData,
    ProjectedAttitudeRecord,
    ProjectedBase,
    ProjectedFocalMechanismData,
    ProjectedFocalMechanismRecord,
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

    /// Line intersections, each as the span it covers along the profile.
    ///
    /// A line generally crosses a section at a point, and the span is then
    /// degenerate; a segment lying along the section trace crosses it over a
    /// stretch, and the span is real. Schema v5 records the pair. An export
    /// older than that holds a single distance, which is read as the degenerate
    /// span it stands for -- and in such a file a stretch is not recoverable at
    /// all, its two ends having been written as two ordinary rows.
    pub fn read_line_intersections(&self) -> Result<Vec<IntersectionRecord>, GeoProfileError> {

        let columns = table_columns(&self.conn, "gp_intersected_lines")?;
        let has_span = columns.iter().any(|c| c == "s_from") && columns.iter().any(|c| c == "s_to");

        // Aliased rather than branched on further down, so that the row mapping
        // and the ordering below are written once.
        let distances = if has_span { "s_from, s_to" } else { "s AS s_from, s AS s_to" };

        let mut stmt = self.conn.prepare(
            &format!(r#"
            SELECT
                rec_id,
                profile_id,
                feat_category,
                {distances},
                extra_json
            FROM gp_intersected_lines
            ORDER BY profile_id, s_from, s_to, rec_id
            "#, distances = distances),
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

    /// Focal mechanisms projected onto the profiles.
    ///
    /// Added to the export at schema v2; absent from anything older, which
    /// answers with an empty list rather than an error, since a file written
    /// before the table existed is not a malformed one.
    pub fn read_projected_focal_mechanisms(
        &self,
    ) -> Result<Vec<ProjectedFocalMechanismRecord>, GeoProfileError> {

        if !table_exists(&self.conn, "gp_projected_focal_mechanisms")? {
            return Ok(Vec::new());
        }

        let mut stmt = self.conn.prepare(
            r#"SELECT
                rec_id,
                profile_id,
                label,
                s,
                z,
                strike,
                dip,
                rake,
                profile_azimuth,
                dist_to_profile,
                src_x,
                src_y,
                src_z,
                src_fid,
                extra_json
            FROM gp_projected_focal_mechanisms
            ORDER BY profile_id, s, rec_id
            "#,
        )?;

        let rows = stmt.query_map([], |row| {

            let src_x: Option<f64> = row.get(10)?;
            let src_y: Option<f64> = row.get(11)?;
            let src_z: Option<f64> = row.get(12)?;

            Ok(ProjectedFocalMechanismRecord {
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
                data: ProjectedFocalMechanismData {
                    strike: row.get(5)?,
                    dip: row.get(6)?,
                    rake: row.get(7)?,
                    profile_azimuth: row.get(8)?,
                },
            })
        })?;

        rows.collect::<Result<Vec<_>, _>>().map_err(Into::into)
    }

    /// The start, break and end vertices of each profile.
    ///
    /// Added at schema v2; older exports answer with an empty list.
    pub fn read_profile_vertices(&self) -> Result<Vec<ProfileVertexRecord>, GeoProfileError> {

        if !table_exists(&self.conn, "gp_profile_vertices")? {
            return Ok(Vec::new());
        }

        let mut stmt = self.conn.prepare(
            r#"SELECT
                vertex_id,
                profile_id,
                vertex_ndx,
                kind,
                s,
                x,
                y,
                lon,
                lat,
                z,
                extra_json
            FROM gp_profile_vertices
            ORDER BY profile_id, vertex_ndx
            "#,
        )?;

        // The kind is carried out as text and turned into the enum below, so
        // that an unrecognised one is reported as the malformed export it is
        // rather than as a column of the wrong type.
        let rows = stmt.query_map([], |row| {
            Ok((
                ProfileVertexRecord {
                    vertex_id: row.get(0)?,
                    profile_id: row.get(1)?,
                    vertex_ndx: row.get(2)?,
                    kind: VertexKind::Start,
                    s: row.get(4)?,
                    x: row.get(5)?,
                    y: row.get(6)?,
                    lon: row.get(7)?,
                    lat: row.get(8)?,
                    z: row.get(9)?,
                    extra_json: row.get(10)?,
                },
                row.get::<_, String>(3)?,
            ))
        })?;

        let mut vertices = Vec::new();
        for row in rows {
            let (mut vertex, kind) = row?;
            vertex.kind = match kind.to_lowercase().as_str() {
                "start" => VertexKind::Start,
                "break" => VertexKind::Break,
                "end" => VertexKind::End,
                other => {
                    return Err(GeoProfileError::InvalidData(format!(
                        "gp_profile_vertices.kind is '{}', not one of start, break, end",
                        other
                    )))
                }
            };
            vertices.push(vertex);
        }

        Ok(vertices)
    }

    /// The plot styles recorded with the export.
    ///
    /// Added at schema v3; older exports answer with an empty list. The JSON is
    /// handed over as written: it is qgSurf's own presentation state, and
    /// giving it a schema here would only fix it in place.
    pub fn read_graphical_params(&self) -> Result<Vec<GraphicalParamsRecord>, GeoProfileError> {

        if !table_exists(&self.conn, "gp_graphical_params")? {
            return Ok(Vec::new());
        }

        let mut stmt = self.conn.prepare(
            r#"SELECT
                params_id,
                result_set_id,
                profile_id,
                params_json
            FROM gp_graphical_params
            ORDER BY result_set_id, profile_id, params_id
            "#,
        )?;

        let rows = stmt.query_map([], |row| {
            Ok(GraphicalParamsRecord {
                params_id: row.get(0)?,
                result_set_id: row.get(1)?,
                // NULL where the styles belong to the whole result set rather
                // than to one profile.
                profile_id: row.get(2)?,
                params_json: row.get(3)?,
            })
        })?;

        rows.collect::<Result<Vec<_>, _>>().map_err(Into::into)
    }

    /// The classification of the source layers, whole.
    ///
    /// Added at schema v4; older exports answer with an empty list. It holds
    /// every class of the source legend, not only those a section happens to
    /// cross, which is what lets a redrawn section keep the ordering and the
    /// colours of the map it came from.
    pub fn read_source_categories(&self) -> Result<Vec<SourceCategoryRecord>, GeoProfileError> {

        if !table_exists(&self.conn, "gp_source_categories")? {
            return Ok(Vec::new());
        }

        let mut stmt = self.conn.prepare(
            r#"SELECT
                category_id,
                result_set_id,
                data_kind,
                position,
                category,
                color
            FROM gp_source_categories
            ORDER BY result_set_id, data_kind, position
            "#,
        )?;

        let rows = stmt.query_map([], |row| {
            Ok(SourceCategoryRecord {
                category_id: row.get(0)?,
                result_set_id: row.get(1)?,
                data_kind: row.get(2)?,
                position: row.get(3)?,
                category: row.get(4)?,
                color: row.get(5)?,
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
            projected_focal_mechanisms: self.read_projected_focal_mechanisms()?,
            profile_vertices: self.read_profile_vertices()?,
            graphical_params: self.read_graphical_params()?,
            source_categories: self.read_source_categories()?,
            sources: self.read_sources()?,
        })
    }

}

