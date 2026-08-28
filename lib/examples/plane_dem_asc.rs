//! Run the plane-grid intersection over an ESRI ASCII DEM and write the trace as
//! CSV on stdout. Exists so the kernel can be checked against a reference output
//! before any Python binding is in place.
//!
//! ```text
//! cargo run --release --example plane_dem_asc -- \
//!     <dem.asc> <x0> <y0> <z0> <dip_dir> <dip_angle>
//! ```
//!
//! The DEM comes through `raster::io`.

use std::env;
use std::io::{BufWriter, Write};

use misah::geometry::point::Point3D;
use misah::raster::intersection::intersect_plane_grid;
use misah::raster::io::read_esri_ascii_grid;
use misah::structural::geol_plane::GeologicalPlane;

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

    let (grid, nodata) = read_esri_ascii_grid(&args[1], 32633)
        .unwrap_or_else(|e| { eprintln!("{}", e); std::process::exit(1) });
    let src_pt = Point3D::from([parse(2), parse(3), parse(4)]);
    let plane = GeologicalPlane::new(parse(5), parse(6))
        .to_plane(src_pt)
        .expect("an attitude always yields a valid normal");

    let out = intersect_plane_grid(&plane, &grid, nodata);

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
