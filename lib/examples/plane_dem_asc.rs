//! Run the plane-grid intersection over an ESRI ASCII DEM and write the trace as
//! CSV on stdout. Exists so the kernel can be checked against a reference output
//! before any Python binding is in place.
//!
//! ```text
//! cargo run --release --example plane_dem_asc -- \
//!     <dem.asc> <x0> <y0> <z0> <dip_dir> <dip_angle>
//! ```
//!
//! It carries its own reader because `raster::io` does not compile yet and is not
//! wired into the module tree.

use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

use ndarray::Array2;

use misah::geometry::point::Point3D;
use misah::raster::geotransform::GeoTransform;
use misah::raster::grid::Grid;
use misah::raster::intersection::intersect_plane_grid;
use misah::structural::geol_plane::GeologicalPlane;

fn read_esri_ascii(path: &str, epsg_code: i32) -> (Grid, f64) {
    let file = File::open(path).unwrap_or_else(|e| panic!("cannot open {}: {}", path, e));
    let mut lines = BufReader::new(file).lines();

    let mut header: HashMap<String, f64> = HashMap::new();
    for _ in 0..6 {
        let line = lines.next().expect("truncated header").expect("bad header line");
        let mut it = line.split_whitespace();
        let key = it.next().expect("missing header key").to_uppercase();
        let value: f64 = it
            .next()
            .expect("missing header value")
            .parse()
            .expect("header value is not a number");
        header.insert(key, value);
    }

    let nrows = header["NROWS"] as usize;
    let ncols = header["NCOLS"] as usize;
    let cellsize = header["CELLSIZE"];

    let mut values = Vec::with_capacity(nrows * ncols);
    for line in lines {
        for token in line.expect("bad data line").split_whitespace() {
            values.push(token.parse::<f64>().expect("elevation is not a number"));
        }
    }
    assert_eq!(values.len(), nrows * ncols, "grid body does not match header");

    let grid = Grid {
        // ESRI ASCII gives the lower-left corner; a GDAL transform starts from the
        // upper-left, so the origin moves up by the full height of the grid.
        transform: GeoTransform {
            data: [
                header["XLLCORNER"],
                cellsize,
                0.0,
                header["YLLCORNER"] + nrows as f64 * cellsize,
                0.0,
                -cellsize,
            ],
        },
        epsg_code,
        data: Array2::from_shape_vec((nrows, ncols), values).expect("shape matches header"),
    };

    (grid, header["NODATA_VALUE"])
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 7 {
        eprintln!(
            "usage: {} <dem.asc> <x0> <y0> <z0> <dip_dir> <dip_angle>",
            args[0]
        );
        std::process::exit(2);
    }

    let parse = |i: usize| args[i].parse::<f64>().expect("bad numeric argument");

    let (grid, nodata) = read_esri_ascii(&args[1], 32633);
    let src_pt = Point3D::from([parse(2), parse(3), parse(4)]);
    let plane = GeologicalPlane::new(parse(5), parse(6))
        .to_plane(src_pt)
        .expect("an attitude always yields a valid normal");

    let out = intersect_plane_grid(&plane, &grid, Some(nodata));

    let stdout = std::io::stdout();
    let mut w = BufWriter::new(stdout.lock());
    writeln!(w, "x,y,z").unwrap();
    for p in &out.points {
        writeln!(w, "{:.6},{:.6},{:.6}", p.x(), p.y(), p.z()).unwrap();
    }
    eprintln!(
        "points: {}  segments: {}",
        out.points.len(),
        out.segments.len()
    );
}
