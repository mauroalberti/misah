//! Run the mesh-grid intersection over an ESRI ASCII DEM and a VTK surface,
//! writing the located attitudes as CSV on stdout. Exists so the kernel can be
//! checked against geoSurfDEM's reference output, whose columns and precision
//! this reproduces.
//!
//! ```text
//! cargo run --release --example mesh_dem_vtk -- <dem.asc> <surface.vtk>
//! ```
//!
//! The DEM comes through `raster::io`; the VTK reader lives here, reading a
//! mesh being no business of a raster module.

use std::env;
use std::fs;
use std::io::{BufWriter, Write};

use misah::geometry::mesh::TriangleMesh;
use misah::geometry::point::Point3D;
use misah::raster::io::read_esri_ascii_grid;
use misah::raster::mesh_intersection::intersect_mesh_grid;

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

    let (grid, nodata) = read_esri_ascii_grid(&args[1], 32633)
        .unwrap_or_else(|e| { eprintln!("{}", e); std::process::exit(1) });
    let mesh = read_vtk_triangle_strips(&args[2]);

    let out = intersect_mesh_grid(&mesh, &grid, nodata);

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
