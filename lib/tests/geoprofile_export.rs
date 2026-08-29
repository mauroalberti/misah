//! The GeoProfiler SQLite reader, against a real qgSurf export.
//!
//! `tests/data/export_timpasanlorenzo.gpkg` is a GeoProfiler export of the
//! Timpa San Lorenzo section, the same dataset the geogst tutorials use. It was
//! committed with the reader and then left unused, which is how
//! `read_line_intersections` came to query two columns that table has never
//! had: schema validation asks for the right ones, the SELECT asked for the
//! polygon table's, and nothing ever ran the two against each other.
//!
//! That fixture declares `schema_version = 1`. `export_synthetic_v5.gpkg` is
//! the other end: one profile, every one of the thirteen tables filled, written
//! by qgSurf's own exporter through its public insert functions -- see
//! `data/make_export_synthetic_v5.py`, which regenerates it. The two projects
//! agree on a file format that nothing else spans, so the fixture is the
//! agreement, and a hand-built copy of the DDL would drift from it unnoticed.
//!
//! Versions 2 to 4 only added tables; v5 changed how a line intersection
//! records its position. Both shapes are read, and both are tested.

use std::path::PathBuf;

use misah::geoprofile::records::vertex::VertexKind;
use misah::geoprofile::sqlite::reader::SqliteGeoProfileReader;

fn fixture() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/export_timpasanlorenzo.gpkg")
}

fn reader() -> SqliteGeoProfileReader {
    SqliteGeoProfileReader::open(fixture()).expect("the fixture opens read-only")
}

fn v5_fixture() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/export_synthetic_v5.gpkg")
}

fn v5_reader() -> SqliteGeoProfileReader {
    SqliteGeoProfileReader::open(v5_fixture()).expect("the v5 fixture opens read-only")
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
fn a_v5_export_is_valid_and_reads_whole() {
    let r = v5_reader();
    r.validate_schema().expect("the span form is a valid export");

    let dataset = r.read_all().expect("every table reads");

    assert_eq!(dataset.result_sets.len(), 1);
    assert_eq!(dataset.profiles.len(), 1);
    assert_eq!(dataset.profile_samples.len(), 11);
    assert_eq!(dataset.sources.len(), 2);
    assert_eq!(dataset.projected_points.len(), 1);
    assert_eq!(dataset.projected_attitudes.len(), 2);
    assert_eq!(dataset.polygon_intersections.len(), 2);
}

#[test]
fn a_v5_line_intersection_keeps_the_span_it_covers() {
    // The case schema v5 exists for. The exporter wrote a crossing at 300, a
    // stretch run along from 600 to 750, and another crossing at 880; before
    // v5 the stretch would have arrived as two rows no different from the two
    // crossings.
    let lines = v5_reader().read_line_intersections().expect("reads");

    assert_eq!(lines.len(), 3);

    assert_eq!((lines[0].s_from, lines[0].s_to), (300.0, 300.0));
    assert_eq!((lines[1].s_from, lines[1].s_to), (600.0, 750.0));
    assert_eq!((lines[2].s_from, lines[2].s_to), (880.0, 880.0));

    assert_eq!(lines[0].category.as_deref(), Some("faglia"));
    assert_eq!(lines[2].category.as_deref(), Some("contatto"));

    // One of the three covers ground rather than a point: that is the whole
    // distinction, and it must survive the file.
    assert_eq!(lines.iter().filter(|l| l.s_to > l.s_from).count(), 1);
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
    let r = v5_reader();

    let mechanisms = r.read_projected_focal_mechanisms().expect("reads");
    assert_eq!(mechanisms.len(), 1);
    assert_eq!(mechanisms[0].base.category.as_deref(), Some("ML 4.3 2026-03-11"));
    assert_eq!(mechanisms[0].data.strike, 135.0);
    assert_eq!(mechanisms[0].data.dip, 60.0);
    assert_eq!(mechanisms[0].data.rake, -90.0);
    assert_eq!(mechanisms[0].data.profile_azimuth, 45.0);
    // A hypocentre eight kilometres down, so the projected z is well below the
    // topography the same profile carries.
    assert!(mechanisms[0].base.z < 0.0);
    assert!(mechanisms[0].base.src_point.is_some());

    let vertices = r.read_profile_vertices().expect("reads");
    assert_eq!(vertices.len(), 3);
    assert_eq!(vertices[0].kind, VertexKind::Start);
    assert_eq!(vertices[1].kind, VertexKind::Break);
    assert_eq!(vertices[2].kind, VertexKind::End);
    assert_eq!(vertices[0].s, 0.0);
    assert_eq!(vertices[0].lon, Some(16.0));
    // The break vertex was written without geographic coordinates.
    assert_eq!(vertices[1].lon, None);
    assert_eq!(vertices[1].z, None);

    let params = r.read_graphical_params().expect("reads");
    assert_eq!(params.len(), 1);
    // NULL profile_id: the styles belong to the whole result set.
    assert_eq!(params[0].profile_id, None);
    assert!(params[0].params_json.contains("vertical_exaggeration"));

    let categories = r.read_source_categories().expect("reads");
    assert_eq!(categories.len(), 4);
    // Ordered by kind, then by the position each class holds in its own legend.
    assert_eq!(categories[0].data_kind, "line_intersections");
    assert_eq!(categories[0].category, "faglia");
    assert_eq!(categories[0].color.as_deref(), Some("#d62728"));
    // A class the source symbology gave no colour.
    assert!(categories.iter().any(|c| c.color.is_none()));
}

#[test]
fn the_two_fixtures_are_the_two_schema_versions_they_claim() {
    // The tests above lean on which version each fixture is, so that is checked
    // rather than assumed -- regenerating one against a newer exporter would
    // otherwise quietly change what they cover.
    let version = |r: &SqliteGeoProfileReader| r.read_result_sets().expect("readable")[0].schema_version;

    assert_eq!(version(&reader()), 1);
    assert_eq!(version(&v5_reader()), 5);
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
