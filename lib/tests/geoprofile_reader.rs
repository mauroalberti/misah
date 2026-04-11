
use misah::geoprofile::sqlite::reader::SqliteGeoProfileReader;

#[test]
fn reads_real_geoprofile_dataset() {
    let path = "tests/data/export_timpa_sanlorenzo.gpkg";

    let reader = SqliteGeoProfileReader::open(path).unwrap();
    reader.validate_schema().unwrap();

    let dataset = reader.read_all().unwrap();

    assert_eq!(dataset.result_sets.len(), 1);
    assert_eq!(dataset.profiles.len(), 3);
    assert_eq!(dataset.line_intersections.len(), 0);

    assert!(!dataset.profile_samples.is_empty());
    assert!(!dataset.polygon_intersections.is_empty());

    let profile_ids: std::collections::HashSet<i64> =
        dataset.profiles.iter().map(|p| p.profile_id).collect();

    for rec in &dataset.profile_samples {
        assert!(profile_ids.contains(&rec.profile_id));
    }

    for rec in &dataset.polygon_intersections {
        assert!(profile_ids.contains(&rec.profile_id));
        assert!(rec.s_from <= rec.s_to);
    }
}
