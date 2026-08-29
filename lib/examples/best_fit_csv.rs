//! Fit geological planes to a field of located points read from a CSV.
//!
//! The consumer of `mesh_dem_vtk`: that example writes where a surface meets
//! the topography, this one reads an attitude back out of those points. It
//! takes any CSV whose first three columns are x, y and z, which is also the
//! form geoSurfDEM's `BestFitGeoplanes` reads, so the two can be run on the
//! same file and their outputs compared column for column.
//!
//! ```sh
//! cargo run --release -p misah --example best_fit_csv -- \
//!     <points.csv> <cell_size> [coincidence_distance] [max_collinearity]
//! ```
//!
//! Columns are written in the reference's order and units. Parsing CSV is left
//! here rather than put in the library for the same reason reading VTK is left
//! to `mesh_dem_vtk`: it is no business of a structural-geology module.

use std::env;
use std::fs;
use std::process;

use misah::geometry::point::Point3D;
use misah::structural::best_fit::{best_fit_geoplanes, DEFAULT_MAX_COLLINEARITY};

fn main() {

    let args: Vec<String> = env::args().collect();

    if args.len() < 3 {
        eprintln!(
            "usage: {} <points.csv> <cell_size> [coincidence_distance] [max_collinearity]",
            args[0]
        );
        process::exit(2);
    }

    let path = &args[1];
    let cell_size: f64 = args[2].parse().unwrap_or_else(|_| {
        eprintln!("cell size '{}' is not a number", args[2]);
        process::exit(2);
    });
    let coincidence: f64 = args.get(3).map_or(0.1, |a| a.parse().unwrap_or(0.1));
    let max_collinearity: f64 =
        args.get(4).map_or(DEFAULT_MAX_COLLINEARITY, |a| a.parse().unwrap_or(DEFAULT_MAX_COLLINEARITY));

    let text = fs::read_to_string(path).unwrap_or_else(|e| {
        eprintln!("cannot read '{}': {}", path, e);
        process::exit(1);
    });

    let points = read_xyz(&text);
    if points.is_empty() {
        eprintln!("no points read from '{}'", path);
        process::exit(1);
    }

    let field = best_fit_geoplanes(&points, cell_size, coincidence, max_collinearity)
        .unwrap_or_else(|| {
            eprintln!("no field: check the cell size and coincidence distance");
            process::exit(1);
        });

    let s = field.stats;
    eprintln!(
        "read {} points, {} distinct, {} non-empty cells of {}x{}",
        s.input_points, s.distinct_points, s.non_empty_cells, field.grid.rows, field.grid.columns
    );
    eprintln!(
        "fitted {}, set aside {} collinear and {} with fewer than three points",
        s.cells_fitted, s.cells_collinear, s.cells_too_few_points
    );

    println!("x, y, pt_num, dip_dir, dip_ang, s1, s2, s3, definition, rms");
    for cell in &field.fits {
        let [s1, s2, s3] = cell.fit.singular_values;
        println!(
            "{:.6}, {:.6}, {}, {:.6}, {:.6}, {:.6}, {:.6}, {:.6}, {:.3}, {:.6}",
            cell.cell_centre.0,
            cell.cell_centre.1,
            cell.fit.points,
            cell.fit.plane.azimuth,
            cell.fit.plane.dip_angle,
            s1,
            s2,
            s3,
            cell.fit.collinearity(),
            cell.fit.rms_distance(),
        );
    }
}

/// The first three comma-separated fields of every line that parses as three
/// numbers. A header is skipped by failing to parse rather than by being
/// counted, which the reference makes the caller do.
fn read_xyz(text: &str) -> Vec<Point3D> {

    text.lines()
        .filter_map(|line| {
            let mut fields = line.split(',').map(|f| f.trim().parse::<f64>());
            match (fields.next(), fields.next(), fields.next()) {
                (Some(Ok(x)), Some(Ok(y)), Some(Ok(z))) => Some(Point3D::from([x, y, z])),
                _ => None,
            }
        })
        .collect()
}
