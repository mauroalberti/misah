//! Run the plane-DEM kernel over an ESRI ASCII grid and write the intersection
//! as CSV on stdout. Exists so the Rust kernel can be cross-checked against the
//! numpy reference before any Python binding is in place.
//!
//! ```text
//! cargo run --release --no-default-features --example plane_dem_asc -- \
//!     <dem.asc> <x0> <y0> <z0> <dip_dir> <dip_angle>
//! ```

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

use misah::kernels::{intersect_plane_dem, GeoTransform};

struct Grid {
    data: Vec<f64>,
    nrows: usize,
    ncols: usize,
    geotransform: GeoTransform,
    nodata: f64,
}

fn read_esri_ascii(path: &str) -> Grid {
    let file = File::open(path).unwrap_or_else(|e| panic!("cannot open {}: {}", path, e));
    let mut lines = BufReader::new(file).lines();

    let mut header = std::collections::HashMap::new();
    for _ in 0..6 {
        let line = lines.next().expect("truncated header").expect("bad header line");
        let mut it = line.split_whitespace();
        let key = it.next().expect("missing header key").to_uppercase();
        let value: f64 = it.next().expect("missing header value").parse().expect("bad value");
        header.insert(key, value);
    }

    let nrows = header["NROWS"] as usize;
    let ncols = header["NCOLS"] as usize;
    let cellsize = header["CELLSIZE"];

    let mut data = Vec::with_capacity(nrows * ncols);
    for line in lines {
        for token in line.expect("bad data line").split_whitespace() {
            data.push(token.parse::<f64>().expect("bad elevation"));
        }
    }
    assert_eq!(data.len(), nrows * ncols, "grid body does not match header");

    Grid {
        data,
        nrows,
        ncols,
        geotransform: GeoTransform {
            x_origin: header["XLLCORNER"],
            pixel_width: cellsize,
            row_rotation: 0.0,
            y_origin: header["YLLCORNER"] + nrows as f64 * cellsize,
            col_rotation: 0.0,
            pixel_height: -cellsize,
        },
        nodata: header["NODATA_VALUE"],
    }
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

    let grid = read_esri_ascii(&args[1]);
    let parse = |i: usize| args[i].parse::<f64>().expect("bad numeric argument");
    let src_pt = [parse(2), parse(3), parse(4)];

    let out = intersect_plane_dem(
        &grid.data,
        grid.nrows,
        grid.ncols,
        &grid.geotransform,
        src_pt,
        parse(5),
        parse(6),
        Some(grid.nodata),
    );

    let stdout = std::io::stdout();
    let mut w = BufWriter::new(stdout.lock());
    writeln!(w, "x,y,z").unwrap();
    for p in &out.points {
        writeln!(w, "{:.6},{:.6},{:.6}", p[0], p[1], p[2]).unwrap();
    }
    eprintln!(
        "points: {}  segments: {}",
        out.points.len(),
        out.segments.len()
    );
}
