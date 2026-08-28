//! The GeoProfiler SQLite reader, against a real qgSurf export.
//!
//! `tests/data/export_timpasanlorenzo.gpkg` is a GeoProfiler export of the
//! Timpa San Lorenzo section, the same dataset the geogst tutorials use. It was
//! committed with the reader and then left unused, which is how
//! `read_line_intersections` came to query two columns that table has never
//! had: schema validation asks for the right ones, the SELECT asked for the
//! polygon table's, and nothing ever ran the two against each other.
//!
//! The fixture declares `schema_version = 1`; qgSurf writes 5 now. Versions 2
//! to 4 only added tables, and v5 changed how a line intersection records its
//! position, so several tests here are about reading both shapes.

use std::path::PathBuf;

use misah::geoprofile::records::vertex::VertexKind;
use misah::geoprofile::sqlite::reader::SqliteGeoProfileReader;

fn fixture() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/export_timpasanlorenzo.gpkg")
}

fn reader() -> SqliteGeoProfileReader {
    SqliteGeoProfileReader::open(fixture()).expect("the fixture opens read-only")
}

#[test]
fn a_real_export_passes_schema_validation() {
    reader().validate_schema().expect("a genuine GeoProfiler export");
}

#[test]
fn every_reader_method_runs_against_a_real_export() {
    // Each is called on its own, so a failure names the table it came from.
    // Four are empty here, which still exercises the statement and its column
    // names -- SQLite resolves those when the statement is prepared, so a query
    // naming a column that does not exist fails on an empty table just as it
    // would on a full one. That is exactly how the s_from/s_to fault behaved.
    let r = reader();

    assert_eq!(r.read_result_sets().expect("result sets").len(), 1);
    assert_eq!(r.read_profiles().expect("profiles").len(), 3);
    assert_eq!(r.read_profile_samples().expect("samples").len(), 1737);
    assert_eq!(r.read_polygon_intersections().expect("polygons").len(), 26);
    assert_eq!(r.read_line_intersections().expect("lines").len(), 0);
    assert_eq!(r.read_projected_attitudes().expect("attitudes").len(), 0);
    assert_eq!(r.read_projected_points().expect("points").len(), 0);
    assert_eq!(r.read_sources().expect("sources").len(), 0);
}

#[test]
fn read_all_returns_the_whole_dataset() {
    // The method a caller actually uses, and the one that stayed broken while
    // seven of its eight parts worked.
    let dataset = reader().read_all().expect("a genuine GeoProfiler export");

    assert_eq!(dataset.result_sets.len(), 1);
    assert_eq!(dataset.profiles.len(), 3);
    assert_eq!(dataset.profile_samples.len(), 1737);
    assert_eq!(dataset.polygon_intersections.len(), 26);
}

#[test]
fn profiles_carry_their_extent_and_their_crs() {
    for p in reader().read_profiles().expect("readable") {
        let s_max = p.s_max.expect("an exported profile knows its length");
        assert!(s_max > 0.0, "profile {} has no length", p.profile_id);

        let (z_min, z_max) = (
            p.z_min.expect("and its lowest point"),
            p.z_max.expect("and its highest"),
        );
        assert!(z_min <= z_max);

        assert!(
            p.wkt_crs.as_deref().map(|c| !c.is_empty()).unwrap_or(false),
            "profile {} carries no CRS",
            p.profile_id
        );
    }
}

#[test]
fn samples_run_in_order_and_agree_with_the_profile_they_belong_to() {
    let dataset = reader().read_all().expect("readable");

    for profile in &dataset.profiles {
        let mut samples: Vec<_> = dataset
            .profile_samples
            .iter()
            .filter(|s| s.profile_id == profile.profile_id)
            .collect();
        samples.sort_by_key(|s| s.idx);

        assert!(!samples.is_empty(), "profile {} has no samples", profile.profile_id);

        let mut previous = f64::NEG_INFINITY;
        for s in &samples {
            assert!(s.s >= previous, "s went backwards at idx {}", s.idx);
            previous = s.s;
        }

        // The profile's declared extent and its samples are written by
        // different parts of the exporter, so their agreeing is worth checking
        // rather than assuming. Tolerances are a metre and a centimetre against
        // a section kilometres long.
        let s_max = profile.s_max.expect("length");
        let last = samples.last().unwrap();
        assert!(
            (last.s - s_max).abs() <= 1.0,
            "profile {} ends at s={} but declares s_max={}",
            profile.profile_id,
            last.s,
            s_max
        );

        let (z_min, z_max) = (profile.z_min.unwrap(), profile.z_max.unwrap());
        for s in &samples {
            assert!(
                s.z >= z_min - 0.01 && s.z <= z_max + 0.01,
                "sample z={} outside the profile's own {}..{}",
                s.z,
                z_min,
                z_max
            );
        }
    }
}

#[test]
fn polygon_intersections_are_ordered_spans_within_their_profile() {
    let dataset = reader().read_all().expect("readable");

    for rec in &dataset.polygon_intersections {
        assert!(
            rec.s_from <= rec.s_to,
            "span runs backwards: {} to {}",
            rec.s_from,
            rec.s_to
        );

        let profile = dataset
            .profiles
            .iter()
            .find(|p| p.profile_id == rec.profile_id)
            .expect("every intersection names a profile that exists");

        assert!(rec.s_to <= profile.s_max.unwrap() + 1.0, "span past the profile end");
        assert!(
            rec.category.as_deref().map(|c| !c.is_empty()).unwrap_or(false),
            "an intersected polygon with no unit name"
        );
    }
}

#[test]
fn a_pre_v5_line_intersection_reads_as_a_degenerate_span() {
    // This fixture predates schema v5, so its gp_intersected_lines holds one
    // distance per row and the reader must present it as the span it stands
    // for. The table is empty here, hence the inserted row -- and its emptiness
    // is exactly what let a broken query survive in the first place.
    let path = temp_copy("misah_geoprofile_lines_v1.gpkg");

    {
        let conn = rusqlite::Connection::open(&path).expect("the copy opens");
        conn.execute(
            "INSERT INTO gp_intersected_lines(profile_id, feat_category, s, extra_json)
             SELECT profile_id, 'faglia', 1234.5, NULL FROM gp_profiles LIMIT 1",
            [],
        )
        .expect("the copy is writable");
    }

    let lines = SqliteGeoProfileReader::open(&path)
        .expect("opens")
        .read_line_intersections()
        .expect("reads");

    assert_eq!(lines.len(), 1);
    assert_eq!(lines[0].s_from, 1234.5);
    assert_eq!(lines[0].s_to, 1234.5);
    assert_eq!(lines[0].category.as_deref(), Some("faglia"));

    let _ = std::fs::remove_file(&path);
}

#[test]
fn a_v5_line_intersection_keeps_the_span_it_covers() {
    // From v5 the table records s_from and s_to. A crossing is a degenerate
    // span; a segment lying along the section trace covers a real one, which is
    // what the older single column could not say.
    let path = temp_copy("misah_geoprofile_lines_v5.gpkg");

    {
        let conn = rusqlite::Connection::open(&path).expect("the copy opens");
        conn.execute("DROP TABLE gp_intersected_lines", []).expect("writable");
        conn.execute(
            "CREATE TABLE gp_intersected_lines (
                 rec_id INTEGER PRIMARY KEY AUTOINCREMENT,
                 profile_id INTEGER NOT NULL,
                 feat_category TEXT DEFAULT '',
                 s_from REAL NOT NULL,
                 s_to REAL NOT NULL,
                 extra_json TEXT DEFAULT NULL)",
            [],
        )
        .expect("writable");
        conn.execute(
            "INSERT INTO gp_intersected_lines(profile_id, feat_category, s_from, s_to)
             SELECT profile_id, 'faglia', 200.0, 250.0 FROM gp_profiles LIMIT 1",
            [],
        )
        .expect("writable");
        conn.execute(
            "INSERT INTO gp_intersected_lines(profile_id, feat_category, s_from, s_to)
             SELECT profile_id, 'contatto', 300.0, 300.0 FROM gp_profiles LIMIT 1",
            [],
        )
        .expect("writable");
    }

    let r = SqliteGeoProfileReader::open(&path).expect("opens");
    r.validate_schema().expect("the span form is a valid export");

    let lines = r.read_line_intersections().expect("reads");

    assert_eq!(lines.len(), 2);
    assert_eq!((lines[0].s_from, lines[0].s_to), (200.0, 250.0));
    assert_eq!((lines[1].s_from, lines[1].s_to), (300.0, 300.0));

    let _ = std::fs::remove_file(&path);
}

#[test]
fn the_tables_added_after_v1_come_back_empty_rather_than_failing() {
    // gp_projected_focal_mechanisms, gp_profile_vertices, gp_graphical_params
    // and gp_source_categories arrived with schema v2 to v4. A file written
    // before them is not a malformed one.
    let r = reader();

    assert!(r.read_projected_focal_mechanisms().expect("absent, not broken").is_empty());
    assert!(r.read_profile_vertices().expect("absent, not broken").is_empty());
    assert!(r.read_graphical_params().expect("absent, not broken").is_empty());
    assert!(r.read_source_categories().expect("absent, not broken").is_empty());

    let dataset = r.read_all().expect("and read_all does not fall over them");
    assert!(dataset.projected_focal_mechanisms.is_empty());
    assert!(dataset.profile_vertices.is_empty());
}

#[test]
fn the_tables_added_after_v1_are_read_when_present() {
    let path = temp_copy("misah_geoprofile_v4_tables.gpkg");

    {
        let conn = rusqlite::Connection::open(&path).expect("the copy opens");
        conn.execute_batch(
            "CREATE TABLE gp_projected_focal_mechanisms (
                 rec_id INTEGER PRIMARY KEY AUTOINCREMENT, profile_id INTEGER NOT NULL,
                 label TEXT DEFAULT '', s REAL NOT NULL, z REAL NOT NULL,
                 strike REAL NOT NULL, dip REAL NOT NULL, rake REAL NOT NULL,
                 profile_azimuth REAL NOT NULL, dist_to_profile REAL, src_x REAL,
                 src_y REAL, src_z REAL, src_fid INTEGER, extra_json TEXT);
             CREATE TABLE gp_profile_vertices (
                 vertex_id INTEGER PRIMARY KEY AUTOINCREMENT, profile_id INTEGER NOT NULL,
                 vertex_ndx INTEGER NOT NULL, kind TEXT NOT NULL, s REAL NOT NULL,
                 x REAL NOT NULL, y REAL NOT NULL, lon REAL, lat REAL, z REAL,
                 extra_json TEXT);
             CREATE TABLE gp_graphical_params (
                 params_id INTEGER PRIMARY KEY AUTOINCREMENT, result_set_id INTEGER NOT NULL,
                 profile_id INTEGER, params_json TEXT NOT NULL);
             CREATE TABLE gp_source_categories (
                 category_id INTEGER PRIMARY KEY AUTOINCREMENT, result_set_id INTEGER NOT NULL,
                 data_kind TEXT NOT NULL, position INTEGER NOT NULL, category TEXT NOT NULL,
                 color TEXT);

             INSERT INTO gp_projected_focal_mechanisms
                 (profile_id, label, s, z, strike, dip, rake, profile_azimuth,
                  dist_to_profile, src_x, src_y, src_z, src_fid)
             SELECT profile_id, 'ML 4.3', 1500.0, -8000.0, 135.0, 60.0, -90.0, 45.0,
                    120.0, 600000.0, 4440000.0, -8000.0, 7 FROM gp_profiles LIMIT 1;

             INSERT INTO gp_profile_vertices (profile_id, vertex_ndx, kind, s, x, y, lon, lat, z)
             SELECT profile_id, 0, 'start', 0.0, 600000.0, 4440000.0, 16.0, 40.1, 850.0
             FROM gp_profiles LIMIT 1;
             INSERT INTO gp_profile_vertices (profile_id, vertex_ndx, kind, s, x, y, lon, lat, z)
             SELECT profile_id, 1, 'break', 500.0, 600400.0, 4440300.0, NULL, NULL, NULL
             FROM gp_profiles LIMIT 1;

             INSERT INTO gp_graphical_params (result_set_id, profile_id, params_json)
             VALUES (1, NULL, '{\"vertical_exaggeration\": 2}');

             INSERT INTO gp_source_categories (result_set_id, data_kind, position, category, color)
             VALUES (1, 'line_intersections', 0, 'faglia', '#ff0000'),
                    (1, 'line_intersections', 1, 'contatto', NULL);",
        )
        .expect("the copy is writable");
    }

    let r = SqliteGeoProfileReader::open(&path).expect("opens");

    let mechanisms = r.read_projected_focal_mechanisms().expect("reads");
    assert_eq!(mechanisms.len(), 1);
    assert_eq!(mechanisms[0].base.category.as_deref(), Some("ML 4.3"));
    assert_eq!(mechanisms[0].data.strike, 135.0);
    assert_eq!(mechanisms[0].data.rake, -90.0);
    assert_eq!(mechanisms[0].data.profile_azimuth, 45.0);
    // src_x/y/z are present, so they become a point rather than three options.
    assert!(mechanisms[0].base.src_point.is_some());

    let vertices = r.read_profile_vertices().expect("reads");
    assert_eq!(vertices.len(), 2);
    assert_eq!(vertices[0].kind, VertexKind::Start);
    assert_eq!(vertices[1].kind, VertexKind::Break);
    assert_eq!(vertices[0].lon, Some(16.0));
    // The break vertex was written without geographic coordinates.
    assert_eq!(vertices[1].lon, None);

    let params = r.read_graphical_params().expect("reads");
    assert_eq!(params.len(), 1);
    // NULL profile_id: the styles belong to the whole result set.
    assert_eq!(params[0].profile_id, None);
    assert!(params[0].params_json.contains("vertical_exaggeration"));

    let categories = r.read_source_categories().expect("reads");
    assert_eq!(categories.len(), 2);
    assert_eq!(categories[0].category, "faglia");
    assert_eq!(categories[0].color.as_deref(), Some("#ff0000"));
    assert_eq!(categories[1].color, None);

    let _ = std::fs::remove_file(&path);
}

#[test]
fn an_unknown_vertex_kind_is_reported_as_bad_data() {
    // The export constrains kind to start/break/end, so anything else did not
    // come from GeoProfiler and should say so rather than be quietly mapped.
    let path = temp_copy("misah_geoprofile_bad_vertex.gpkg");

    {
        let conn = rusqlite::Connection::open(&path).expect("the copy opens");
        conn.execute_batch(
            "CREATE TABLE gp_profile_vertices (
                 vertex_id INTEGER PRIMARY KEY AUTOINCREMENT, profile_id INTEGER NOT NULL,
                 vertex_ndx INTEGER NOT NULL, kind TEXT NOT NULL, s REAL NOT NULL,
                 x REAL NOT NULL, y REAL NOT NULL, lon REAL, lat REAL, z REAL,
                 extra_json TEXT);
             INSERT INTO gp_profile_vertices (profile_id, vertex_ndx, kind, s, x, y)
             VALUES (1, 0, 'middle', 0.0, 1.0, 2.0);",
        )
        .expect("writable");
    }

    let err = SqliteGeoProfileReader::open(&path)
        .expect("opens")
        .read_profile_vertices()
        .unwrap_err();

    assert!(format!("{}", err).contains("middle"), "unexpected error: {}", err);

    let _ = std::fs::remove_file(&path);
}

#[test]
fn a_null_in_a_defaulted_text_column_does_not_fail_the_read() {
    // The export schema gives these columns DEFAULT '', so the empty string is
    // how the format says "nothing here" -- but NULL remains storable, and read
    // straight into a String it would fail the whole read over one stray value.
    let path = temp_copy("misah_geoprofile_nulls.gpkg");

    {
        let conn = rusqlite::Connection::open(&path).expect("the copy opens");
        conn.execute("UPDATE gp_result_sets SET description = NULL, distance_units = NULL", [])
            .expect("writable");
        conn.execute("UPDATE gp_profiles SET profile_name = NULL", [])
            .expect("writable");
    }

    let r = SqliteGeoProfileReader::open(&path).expect("opens");

    let sets = r.read_result_sets().expect("a NULL description is not a broken export");
    assert_eq!(sets[0].description, "");
    assert_eq!(sets[0].distance_units, "");

    let profiles = r.read_profiles().expect("nor a NULL profile name");
    assert!(profiles.iter().all(|p| p.profile_name.is_empty()));

    let _ = std::fs::remove_file(&path);
}

#[test]
fn a_newer_schema_is_still_readable() {
    // qgSurf writes schema_version 4; this fixture is 1. Versions 2 to 4 only
    // added tables -- gp_profile_vertices, gp_graphical_params,
    // gp_source_categories, gp_projected_focal_mechanisms -- without touching a
    // column of the nine read here. Validation must therefore not turn on the
    // version number, or every current export would be refused.
    let path = temp_copy("misah_geoprofile_v4.gpkg");

    {
        let conn = rusqlite::Connection::open(&path).expect("the copy opens");
        conn.execute("UPDATE gp_meta SET value = '4' WHERE key = 'schema_version'", [])
            .expect("writable");
        conn.execute("UPDATE gp_result_sets SET schema_version = 4", [])
            .expect("writable");
    }

    let r = SqliteGeoProfileReader::open(&path).expect("the copy opens read-only");
    r.validate_schema()
        .expect("a version bump alone must not make an export unreadable");
    let dataset = r.read_all().expect("and its rows must still come back");
    assert_eq!(dataset.profiles.len(), 3);

    let _ = std::fs::remove_file(&path);
}

/// A writable copy of the fixture, under a name of its own so that tests
/// running in parallel do not share one.
fn temp_copy(name: &str) -> PathBuf {
    let bytes = std::fs::read(fixture()).expect("the fixture is readable");
    let mut path = std::env::temp_dir();
    path.push(name);
    std::fs::write(&path, bytes).expect("a writable temp directory");
    path
}
