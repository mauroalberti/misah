
use rusqlite::Connection;

use crate::geoprofile::error::GeoProfileError;

pub fn table_exists(conn: &Connection, table_name: &str) -> Result<bool, GeoProfileError> {
    let mut stmt = conn.prepare(
        "SELECT EXISTS(
            SELECT 1
            FROM sqlite_master
            WHERE type='table' AND name=?1
        )"
    )?;

    let exists: i64 = stmt.query_row([table_name], |row| row.get(0))?;
    Ok(exists != 0)
}

pub fn table_columns(conn: &Connection, table_name: &str) -> Result<Vec<String>, GeoProfileError> {

    let pragma = format!("PRAGMA table_info({})", table_name);
    let mut stmt = conn.prepare(&pragma)?;
    let rows = stmt.query_map([], |row| row.get::<_, String>(1))?;

    let mut cols = Vec::new();
    for c in rows {
        cols.push(c?);
    }
    Ok(cols)
}

pub fn require_table(conn: &Connection, table_name: &str) -> Result<(), GeoProfileError> {
    if !table_exists(conn, table_name)? {
        return Err(GeoProfileError::MissingTable(table_name.to_string()));
    }
    Ok(())
}

pub fn require_columns(
    conn: &Connection,
    table_name: &str,
    required_columns: &[&str],
) -> Result<(), GeoProfileError> {
    require_table(conn, table_name)?;
    let cols = table_columns(conn, table_name)?;

    for required in required_columns {
        if !cols.iter().any(|c| c == required) {
            return Err(GeoProfileError::MissingColumn {
                table: table_name.to_string(),
                column: (*required).to_string(),
            });
        }
    }

    Ok(())
}

/// Check that an export carries the tables and columns the reader needs.
///
/// Only the nine tables present since schema v1 are required. qgSurf has since
/// added `gp_profile_vertices`, `gp_projected_focal_mechanisms`,
/// `gp_graphical_params`, `gp_source_categories` (v2 to v4) and the span form
/// of `gp_intersected_lines` (v5); all of those are read when present and
/// passed over when not, so that one reader serves every export written so far.
/// Requiring them would refuse the older files this validation exists to admit.
pub fn validate_geoprofile_schema(conn: &Connection) -> Result<(), GeoProfileError> {
    require_columns(
        conn,
        "gp_result_sets",
        &[
            "result_set_id",
            "name",
            "description",
            "created_utc",
            "schema_version",
            "distance_units",
            "extra_json",
        ],
    )?;

    require_columns(
        conn,
        "gp_profiles",
        &[
            "profile_id",
            "result_set_id",
            "source_profile_fid",
            "profile_name",
            "wkt_crs",
            "s_max",
            "z_min",
            "z_max",
            "extra_json",
        ],
    )?;

    require_columns(
        conn,
        "gp_profile_samples",
        &[
            "sample_id",
            "profile_id",
            "idx",
            "s",
            "z",
            "x",
            "y",
            "extra_json",
        ],
    )?;

    require_columns(
        conn,
        "gp_intersected_polygons",
        &[
            "rec_id",
            "profile_id",
            "unit_name",
            "s_from",
            "s_to",
            "extra_json",
        ],
    )?;

    // The distance columns are checked apart: schema v5 replaced the single s
    // with the span s_from/s_to, and a migrated database carries both.
    require_columns(
        conn,
        "gp_intersected_lines",
        &[
            "rec_id",
            "profile_id",
            "feat_category",
            "extra_json",
        ],
    )?;

    let line_columns = table_columns(conn, "gp_intersected_lines")?;
    let has_span = line_columns.iter().any(|c| c == "s_from")
        && line_columns.iter().any(|c| c == "s_to");
    if !has_span && !line_columns.iter().any(|c| c == "s") {
        return Err(GeoProfileError::MissingColumn {
            table: "gp_intersected_lines".to_string(),
            column: "s_from/s_to, or s in an export older than schema v5".to_string(),
        });
    }


    require_columns(
        conn,
        "gp_projected_attitudes",
        &[
            "rec_id",
            "profile_id",
            "label",
            "s",
            "z",
            "slope_degr",
            "down_sense",
            "src_dip_dir",
            "src_dip_ang",
            "dist_to_profile",
            "src_x",
            "src_y",
            "src_z",
            "src_fid",
            "extra_json",
        ],
    )?;


    require_columns(
        conn,
        "gp_projected_points",
        &[
            "rec_id",
            "profile_id",
            "label",
            "s",
            "z",
            "dist_to_profile",
            "src_x",
            "src_y",
            "src_z",
            "src_fid",
            "extra_json",
        ],
    )?;

    require_columns(
        conn,
        "gp_sources",
        &[
            "source_id",
            "result_set_id",
            "source_type",
            "source_name",
            "source_uri",
            "wkt_crs",
            "parameters_json",
        ],
    )?;

    require_columns(
        conn,
        "gp_meta",
        &[
            "key",
            "value",
        ],
    )?;

    Ok(())

}
