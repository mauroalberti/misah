//! Reading rasters from file.
//!
//! Only the ESRI ASCII grid is handled: it is what the geoSurfDEM test data and
//! the qgSurf exports come in, and it needs no library to parse. Anything
//! wider -- GeoTIFF and the rest -- means GDAL, which a crate meant to build as
//! a self-contained wheel should not take on lightly.

use std::collections::HashMap;
use std::fs;
use std::path::Path;

use ndarray::Array2;
use thiserror::Error;

use crate::raster::geotransform::GeoTransform;
use crate::raster::grid::Grid;

#[derive(Debug, Error)]
pub enum RasterIoError {
    #[error("cannot read {path}: {source}")]
    Unreadable { path: String, source: std::io::Error },

    #[error("{source_name}: invalid header: {message}")]
    InvalidHeader { source_name: String, message: String },

    #[error("{source_name}: invalid grid body: {message}")]
    InvalidBody { source_name: String, message: String },
}

/// The header keys an ESRI ASCII grid may carry, all of which take one number.
///
/// Reading the header by key rather than by line position is what lets the
/// optional ones be absent: the two coordinate keys come in a corner and a
/// centre spelling, and `NODATA_VALUE` need not be there at all.
const HEADER_KEYS: [&str; 8] = [
    "NCOLS",
    "NROWS",
    "XLLCORNER",
    "YLLCORNER",
    "XLLCENTER",
    "YLLCENTER",
    "CELLSIZE",
    "NODATA_VALUE",
];

/// Read an ESRI ASCII grid, returning it with its nodata value if it declared
/// one.
///
/// The nodata value is `Option` rather than the -9999 the format nominally
/// defaults to. Applying that default to a file that never asked for it would
/// silently punch holes in bathymetry, where -9999 m is an ordinary depth
/// rather than a marker; a caller wanting the default can supply it knowingly.
pub fn read_esri_ascii_grid(
    path: impl AsRef<Path>,
) -> Result<(Grid, Option<f64>), RasterIoError> {

    let path = path.as_ref();
    let name = path.display().to_string();

    let text = fs::read_to_string(path).map_err(|source| RasterIoError::Unreadable {
        path: name.clone(),
        source,
    })?;

    parse_esri_ascii_grid(&text, &name)
}

/// Parse an ESRI ASCII grid already in memory.
///
/// `source_name` names the origin in any error message, and is not otherwise
/// used; it exists so a grid pulled out of an archive or a download can report
/// itself as something more useful than "input".
///
/// Parsing runs over the whitespace-separated tokens rather than line by line,
/// so nothing depends on how the writer chose to wrap the values -- one row per
/// line and one value per line are both common, and both are read.
pub fn parse_esri_ascii_grid(
    text: &str,
    source_name: &str,
) -> Result<(Grid, Option<f64>), RasterIoError> {

    let header_error = |message: String| RasterIoError::InvalidHeader {
        source_name: source_name.to_string(),
        message,
    };
    let body_error = |message: String| RasterIoError::InvalidBody {
        source_name: source_name.to_string(),
        message,
    };

    let mut tokens = text.split_whitespace().peekable();
    let mut header: HashMap<String, f64> = HashMap::new();

    // The header ends at the first token that is not one of its keys, which is
    // the first elevation.
    while let Some(token) = tokens.peek() {

        let key = token.to_uppercase();
        if !HEADER_KEYS.contains(&key.as_str()) {
            break;
        }
        tokens.next();

        let raw = tokens
            .next()
            .ok_or_else(|| header_error(format!("{} has no value", key)))?;
        let value: f64 = raw
            .parse()
            .map_err(|_| header_error(format!("{} is not a number: {:?}", key, raw)))?;

        header.insert(key, value);
    }

    let required = |key: &str| -> Result<f64, RasterIoError> {
        header
            .get(key)
            .copied()
            .ok_or_else(|| header_error(format!("{} is missing", key)))
    };

    let count = |key: &str| -> Result<usize, RasterIoError> {
        let value = required(key)?;
        if value < 1.0 || value.fract() != 0.0 {
            return Err(header_error(format!(
                "{} must be a positive whole number, found {}",
                key, value
            )));
        }
        Ok(value as usize)
    };

    let ncols = count("NCOLS")?;
    let nrows = count("NROWS")?;

    let cellsize = required("CELLSIZE")?;
    // Finiteness is tested apart from the sign: a NaN compares false against
    // everything, so `cellsize <= 0.0` on its own would let one through and
    // leave every coordinate the transform produces a NaN.
    if !cellsize.is_finite() || cellsize <= 0.0 {
        return Err(header_error(format!(
            "CELLSIZE must be a positive finite number, found {}",
            cellsize
        )));
    }

    // The corner spelling gives the outer edge of the lower-left cell, the
    // centre spelling its sampled point, which is half a cell in from that.
    let edge = |corner: &str, centre: &str| -> Result<f64, RasterIoError> {
        match (header.get(corner), header.get(centre)) {
            (Some(&v), _) => Ok(v),
            (None, Some(&v)) => Ok(v - cellsize / 2.0),
            (None, None) => Err(header_error(format!(
                "neither {} nor {} is present",
                corner, centre
            ))),
        }
    };

    let x_left = edge("XLLCORNER", "XLLCENTER")?;
    let y_bottom = edge("YLLCORNER", "YLLCENTER")?;

    let expected = nrows * ncols;
    let mut values = Vec::with_capacity(expected);

    for (index, raw) in tokens.enumerate() {
        let value: f64 = raw
            .parse()
            .map_err(|_| body_error(format!("value {} is not a number: {:?}", index, raw)))?;
        values.push(value);
    }

    if values.len() != expected {
        return Err(body_error(format!(
            "header declares {} values ({} rows x {} columns) but the body holds {}",
            expected,
            nrows,
            ncols,
            values.len()
        )));
    }

    let grid = Grid {
        // An ESRI header is written from the lower-left corner, a GDAL
        // transform from the upper-left one, so the origin moves up by the full
        // height of the grid and the row step is negative.
        transform: GeoTransform {
            data: [
                x_left,
                cellsize,
                0.0,
                y_bottom + nrows as f64 * cellsize,
                0.0,
                -cellsize,
            ],
        },
        data: Array2::from_shape_vec((nrows, ncols), values)
            .expect("the value count was just checked against the declared shape"),
    };

    Ok((grid, header.get("NODATA_VALUE").copied()))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two rows of three cells, one metre across, lower-left corner at the
    /// origin.
    const SMALL: &str = "\
NCOLS 3
NROWS 2
XLLCORNER 0.0
YLLCORNER 0.0
CELLSIZE 1.0
NODATA_VALUE -9999
10 11 12
20 21 22
";

    fn parse(text: &str) -> Result<(Grid, Option<f64>), RasterIoError> {
        parse_esri_ascii_grid(text, "test")
    }

    #[test]
    fn a_grid_is_read_in_row_major_order_from_the_top() {
        let (grid, _) = parse(SMALL).expect("a well-formed grid");

        assert_eq!(grid.data.dim(), (2, 3));
        // The first row of the body is the northern one.
        assert_eq!(grid.data[[0, 0]], 10.0);
        assert_eq!(grid.data[[1, 2]], 22.0);
    }

    #[test]
    fn the_transform_puts_the_first_node_at_the_top_left_cell_centre() {
        let (grid, _) = parse(SMALL).expect("a well-formed grid");

        // Two rows of unit cells above y = 0, so the top edge is at y = 2 and
        // the first cell centre half a cell in from the top-left corner.
        assert_eq!(grid.transform.node_xy(0, 0), (0.5, 1.5));
        assert_eq!(grid.transform.node_xy(1, 2), (2.5, 0.5));
    }

    #[test]
    fn the_centre_spelling_shifts_the_origin_by_half_a_cell() {
        // Same grid, its lower-left cell named by its centre instead of its
        // corner: the cell centres must land in the same places.
        let by_centre = SMALL
            .replace("XLLCORNER 0.0", "XLLCENTER 0.5")
            .replace("YLLCORNER 0.0", "YLLCENTER 0.5");

        let (from_corner, _) = parse(SMALL).unwrap();
        let (from_centre, _) = parse(&by_centre).unwrap();

        assert_eq!(
            from_corner.transform.node_xy(0, 0),
            from_centre.transform.node_xy(0, 0)
        );
        assert_eq!(from_corner.transform.data, from_centre.transform.data);
    }

    #[test]
    fn a_declared_nodata_value_is_returned() {
        let (_, nodata) = parse(SMALL).unwrap();

        assert_eq!(nodata, Some(-9999.0));
    }

    #[test]
    fn an_absent_nodata_value_stays_absent() {
        // No -9999 is assumed on the file's behalf: at sea that is a depth, not
        // a hole.
        let without = SMALL.replace("NODATA_VALUE -9999\n", "");
        let (_, nodata) = parse(&without).expect("nodata is optional");

        assert_eq!(nodata, None);
    }

    #[test]
    fn the_body_may_be_wrapped_any_way_at_all() {
        let one_per_line = "\
NCOLS 3
NROWS 2
XLLCORNER 0.0
YLLCORNER 0.0
CELLSIZE 1.0
10
11 12 20
21
22
";
        let (grid, _) = parse(one_per_line).expect("wrapping is not meaningful");

        assert_eq!(grid.data[[0, 1]], 11.0);
        assert_eq!(grid.data[[1, 0]], 20.0);
    }

    #[test]
    fn header_keys_are_read_case_insensitively() {
        let lower = SMALL.replace("NCOLS", "ncols").replace("NROWS", "nRows");

        assert!(parse(&lower).is_ok());
    }

    #[test]
    fn a_missing_dimension_is_reported_as_a_header_fault() {
        let without = SMALL.replace("NROWS 2\n", "");
        let err = parse(&without).unwrap_err();

        assert!(
            matches!(&err, RasterIoError::InvalidHeader { message, .. } if message.contains("NROWS")),
            "unexpected error: {}", err
        );
    }

    #[test]
    fn a_missing_pair_of_coordinate_keys_is_reported() {
        let without = SMALL.replace("XLLCORNER 0.0\n", "");
        let err = parse(&without).unwrap_err();

        assert!(
            matches!(&err, RasterIoError::InvalidHeader { message, .. } if message.contains("XLLCENTER")),
            "unexpected error: {}", err
        );
    }

    #[test]
    fn a_zero_cell_size_is_rejected() {
        // It would make the transform singular, and every inverse lookup with
        // it.
        let flat = SMALL.replace("CELLSIZE 1.0", "CELLSIZE 0.0");
        let err = parse(&flat).unwrap_err();

        assert!(
            matches!(&err, RasterIoError::InvalidHeader { message, .. } if message.contains("CELLSIZE")),
            "unexpected error: {}", err
        );
    }

    #[test]
    fn a_cell_size_that_is_not_a_number_is_rejected() {
        // NaN compares false against every bound, so a sign test alone would
        // pass it straight into the transform.
        let broken = SMALL.replace("CELLSIZE 1.0", "CELLSIZE NaN");
        let err = parse(&broken).unwrap_err();

        assert!(
            matches!(&err, RasterIoError::InvalidHeader { message, .. } if message.contains("CELLSIZE")),
            "unexpected error: {}", err
        );
    }

    #[test]
    fn a_fractional_row_count_is_rejected() {
        let odd = SMALL.replace("NROWS 2", "NROWS 2.5");
        let err = parse(&odd).unwrap_err();

        assert!(matches!(err, RasterIoError::InvalidHeader { .. }));
    }

    #[test]
    fn a_header_key_without_a_value_is_reported() {
        let truncated = "NCOLS 3\nNROWS\n";
        let err = parse(truncated).unwrap_err();

        assert!(
            matches!(&err, RasterIoError::InvalidHeader { message, .. } if message.contains("no value")),
            "unexpected error: {}", err
        );
    }

    #[test]
    fn a_body_shorter_than_the_header_declares_is_reported() {
        let short = SMALL.replace("20 21 22\n", "20 21\n");
        let err = parse(&short).unwrap_err();

        assert!(
            matches!(&err, RasterIoError::InvalidBody { message, .. } if message.contains("holds 5")),
            "unexpected error: {}", err
        );
    }

    #[test]
    fn a_body_longer_than_the_header_declares_is_reported() {
        let long = format!("{SMALL}30 31 32\n");
        let err = parse(&long).unwrap_err();

        assert!(matches!(err, RasterIoError::InvalidBody { .. }));
    }

    #[test]
    fn a_non_numeric_elevation_is_reported_with_its_position() {
        let broken = SMALL.replace("20 21 22", "20 NaN? 22");
        let err = parse(&broken).unwrap_err();

        assert!(
            matches!(&err, RasterIoError::InvalidBody { message, .. } if message.contains("value 4")),
            "unexpected error: {}", err
        );
    }

    #[test]
    fn a_missing_file_is_reported_with_its_path() {
        let err = read_esri_ascii_grid("/nonexistent/dem.asc").unwrap_err();

        assert!(
            matches!(&err, RasterIoError::Unreadable { path, .. } if path.contains("dem.asc")),
            "unexpected error: {}", err
        );
    }
}
