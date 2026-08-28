//! The GeoProfiler SQLite reader, against a real qgSurf export.
//!
//! `tests/data/export_timpasanlorenzo.gpkg` is a GeoProfiler export of the
//! Timpa San Lorenzo section, the same dataset the geogst tutorials use. It was
//! committed with the reader and then left unused, which is how
//! `read_line_intersections` came to query two columns that table has never
//! had: schema validation asks for the right ones, the SELECT asked for the
//! polygon table's, and nothing ever ran the two against each other.
//!
//! The fixture declares `schema_version = 1`; qgSurf writes 4 now, which is
//! what the last test here is about.

use std::path::PathBuf;

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
fn a_line_intersection_reads_as_a_span_from_its_own_distance() {
    // gp_intersected_lines stores one distance per row, so both ends of the
    // record hold it. The fixture has no rows, so the mapping is checked on a
    // copy with one inserted -- also the only coverage this table gets, since
    // its emptiness is what let a broken query survive.
    let path = temp_copy("misah_geoprofile_lines.gpkg");

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
