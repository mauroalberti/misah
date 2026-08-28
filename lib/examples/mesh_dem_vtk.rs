//! Run the mesh-grid intersection over an ESRI ASCII DEM and a VTK surface,
//! writing the located attitudes as CSV on stdout. Exists so the kernel can be
//! checked against geoSurfDEM's reference output, whose columns and precision
//! this reproduces.
//!
//! ```text
//! cargo run --release --example mesh_dem_vtk -- <dem.asc> <surface.vtk>
//! ```
//!
//! It carries its own readers: `raster::io` does not compile yet and is not
//! wired into the module tree, and reading VTK is no business of a raster
//! kernel.

use std::env;
use std::fs;
use std::io::{BufWriter, Write};

use ndarray::Array2;

use misah::geometry::mesh::TriangleMesh;
use misah::geometry::point::Point3D;
use misah::raster::geotransform::GeoTransform;
use misah::raster::grid::Grid;
use misah::raster::mesh_intersection::intersect_mesh_grid;

fn read_esri_ascii(path: &str, epsg_code: i32) -> (Grid, f64) {

    let text = fs::read_to_string(path).unwrap_or_else(|e| panic!("cannot read {}: {}", path, e));
    let mut tokens = text.split_whitespace();

    let mut header = std::collections::HashMap::new();
    for _ in 0..6 {
        let key = tokens.next().expect("truncated header").to_uppercase();
        let value: f64 = tokens
            .next()
            .expect("header key without a value")
            .parse()
            .expect("header value is not a number");
        header.insert(key, value);
    }

    let nrows = header["NROWS"] as usize;
    let ncols = header["NCOLS"] as usize;
    let cellsize = header["CELLSIZE"];

    let values: Vec<f64> = tokens
        .map(|t| t.parse::<f64>().expect("elevation is not a number"))
        .collect();
    assert_eq!(values.len(), nrows * ncols, "grid body does not match header");

    let grid = Grid {
        // ESRI ASCII gives the lower-left corner; a GDAL transform starts from
        // the upper-left, so the origin moves up by the full height of the grid.
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

/// Read an ASCII VTK `POLYDATA` surface made of triangle strips -- what qgSurf
/// writes for a geological surface.
///
/// Parsed over the whitespace-separated token stream rather than line by line,
/// so nothing depends on how the writer chose to wrap its numbers.
fn read_vtk_triangle_strips(path: &str) -> TriangleMesh {

    let text = fs::read_to_string(path).unwrap_or_else(|e| panic!("cannot read {}: {}", path, e));
    let tokens: Vec<&str> = text.split_whitespace().collect();

    let number = |i: usize| -> f64 {
        tokens[i]
            .parse()
            .unwrap_or_else(|_| panic!("expected a number at token {}, found {:?}", i, tokens[i]))
    };
    let count = |i: usize| -> usize {
        tokens[i]
            .parse()
            .unwrap_or_else(|_| panic!("expected a count at token {}, found {:?}", i, tokens[i]))
    };

    let points_at = tokens
        .iter()
        .position(|t| *t == "POINTS")
        .expect("no POINTS section: not an ASCII VTK POLYDATA file");
    let num_points = count(points_at + 1);

    // POINTS <n> <datatype>, then 3n coordinates.
    let first_coord = points_at + 3;
    let vertices: Vec<Point3D> = (0..num_points)
        .map(|n| {
            let i = first_coord + 3 * n;
            Point3D::from([number(i), number(i + 1), number(i + 2)])
        })
        .collect();

    let strips_at = tokens.iter().position(|t| *t == "TRIANGLE_STRIPS").expect(
        "no TRIANGLE_STRIPS section: only the triangle-strip form of POLYDATA is read here",
    );
    let num_strips = count(strips_at + 1);

    // TRIANGLE_STRIPS <n> <total>, then each strip as its length followed by
    // that many vertex indices.
    let mut cursor = strips_at + 3;
    let mut strips = Vec::with_capacity(num_strips);
    for _ in 0..num_strips {
        let len = count(cursor);
        cursor += 1;
        strips.push((0..len).map(|k| count(cursor + k)).collect::<Vec<usize>>());
        cursor += len;
    }

    TriangleMesh::from_triangle_strips(vertices, &strips).expect("VTK indices within the point pool")
}

fn main() {

    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        eprintln!("usage: {} <dem.asc> <surface.vtk>", args[0]);
        std::process::exit(2);
    }

    let (grid, nodata) = read_esri_ascii(&args[1], 32633);
    let mesh = read_vtk_triangle_strips(&args[2]);

    let out = intersect_mesh_grid(&mesh, &grid, Some(nodata));

    let stdout = std::io::stdout();
    let mut w = BufWriter::new(stdout.lock());
    writeln!(w, "x,y,z,dipdir,dipang").unwrap();
    for i in &out.intersections {
        writeln!(
            w,
            "{:.2},{:.2},{:.2},{:.2},{:.2}",
            i.point.x(),
            i.point.y(),
            i.point.z(),
            i.attitude.azimuth,
            i.attitude.dip_angle
        )
        .unwrap();
    }

    let s = out.stats;
    eprintln!(
        "mesh triangles: {} ({} degenerate, {} off the grid)",
        s.mesh_triangles, s.degenerate_mesh_triangles, s.mesh_triangles_outside_grid
    );
    eprintln!(
        "DEM triangles compared: {}  coplanar sides: {}  repeats suppressed: {}",
        s.dem_triangle_pairs, s.coplanar_sides, s.duplicate_crossings
    );
    eprintln!("intersection points: {}", out.len());
}
