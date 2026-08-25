//! Plane-DEM intersection.
//!
//! The intersection is the zero-level set of the signed distance
//!
//! ```text
//! f(x, y) = n . (X - P0),    X = (x, y, z_dem(x, y))
//! ```
//!
//! sampled at the DEM nodes and extracted by marching squares. Working on the
//! signed distance rather than on a `z = z(x, y)` plane expression keeps vertical
//! planes -- which have no single-valued z -- as an ordinary case.
//!
//! Mirrors `python/misah_ref/plane_dem.py`; `tests/test_plane_dem.py` validates
//! both against the geoSurfDEM golden dataset.

use std::collections::HashMap;

/// GDAL-style affine transform, mapping (col, row) offsets to (x, y).
#[derive(Clone, Copy, Debug)]
pub struct GeoTransform {
    pub x_origin: f64,
    pub pixel_width: f64,
    pub row_rotation: f64,
    pub y_origin: f64,
    pub col_rotation: f64,
    pub pixel_height: f64,
}

impl GeoTransform {
    /// DEM values are point samples at cell centres, hence the +0.5 offsets.
    #[inline]
    pub fn node_xy(&self, row: usize, col: usize) -> (f64, f64) {
        let (r, c) = (row as f64 + 0.5, col as f64 + 0.5);
        (
            self.x_origin + c * self.pixel_width + r * self.row_rotation,
            self.y_origin + c * self.col_rotation + r * self.pixel_height,
        )
    }
}

/// Upward-pointing unit normal of a geological plane, in (East, North, Up).
///
/// Sanity anchors: a horizontal plane (dip 0) gives (0, 0, 1); a vertical plane
/// (dip 90) gives a horizontal normal pointing along the dip direction.
pub fn plane_normal(dip_dir_degr: f64, dip_angle_degr: f64) -> [f64; 3] {
    let az = dip_dir_degr.to_radians();
    let dip = dip_angle_degr.to_radians();
    [az.sin() * dip.sin(), az.cos() * dip.sin(), dip.cos()]
}

/// Corner numbering: 0 = (r, c), 1 = (r, c+1), 2 = (r+1, c+1), 3 = (r+1, c).
/// Edges are named after the corner pair they join.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Edge {
    Top,
    Right,
    Bottom,
    Left,
}

impl Edge {
    #[inline]
    fn corners(self) -> (usize, usize) {
        match self {
            Edge::Top => (0, 1),
            Edge::Right => (1, 2),
            Edge::Bottom => (3, 2),
            Edge::Left => (0, 3),
        }
    }
}

/// Marching-squares chords per case index; bit i is set when corner i is positive.
/// Complementary cases (i and 15-i) share their chords.
#[inline]
fn case_edges(case: u8) -> &'static [(Edge, Edge)] {
    match case {
        0 | 15 => &[],
        1 | 14 => &[(Edge::Left, Edge::Top)],
        2 | 13 => &[(Edge::Top, Edge::Right)],
        3 | 12 => &[(Edge::Left, Edge::Right)],
        4 | 11 => &[(Edge::Right, Edge::Bottom)],
        6 | 9 => &[(Edge::Top, Edge::Bottom)],
        7 | 8 => &[(Edge::Left, Edge::Bottom)],
        5 => &[(Edge::Left, Edge::Top), (Edge::Right, Edge::Bottom)],
        10 => &[(Edge::Left, Edge::Bottom), (Edge::Top, Edge::Right)],
        _ => unreachable!(),
    }
}

/// The two alternative pairings for a saddle cell.
#[inline]
fn saddle_alternative(case: u8) -> &'static [(Edge, Edge)] {
    match case {
        5 => &[(Edge::Left, Edge::Bottom), (Edge::Top, Edge::Right)],
        10 => &[(Edge::Left, Edge::Top), (Edge::Right, Edge::Bottom)],
        _ => unreachable!(),
    }
}

pub struct Intersection {
    /// Intersection vertices, as (x, y, z) triples.
    pub points: Vec<[f64; 3]>,
    /// Index pairs into `points`, one per marching-squares chord.
    pub segments: Vec<[usize; 2]>,
}

/// Intersect an unbounded geological plane with a DEM.
///
/// `dem` is row-major with `nrows * ncols` elements, rows running north to south.
/// Non-finite elevations, and those matching `nodata`, are excluded.
pub fn intersect_plane_dem(
    dem: &[f64],
    nrows: usize,
    ncols: usize,
    geotransform: &GeoTransform,
    src_pt: [f64; 3],
    dip_dir_degr: f64,
    dip_angle_degr: f64,
    nodata: Option<f64>,
) -> Intersection {
    assert_eq!(dem.len(), nrows * ncols, "DEM length does not match its shape");

    let normal = plane_normal(dip_dir_degr, dip_angle_degr);

    let elevation = |row: usize, col: usize| -> f64 {
        let z = dem[row * ncols + col];
        match nodata {
            Some(nd) if z == nd => f64::NAN,
            _ => z,
        }
    };

    let signed_distance = |row: usize, col: usize| -> f64 {
        let (x, y) = geotransform.node_xy(row, col);
        let z = elevation(row, col);
        normal[0] * (x - src_pt[0]) + normal[1] * (y - src_pt[1]) + normal[2] * (z - src_pt[2])
    };

    let mut points: Vec<[f64; 3]> = Vec::new();
    let mut segments: Vec<[usize; 2]> = Vec::new();
    // Vertices are shared between adjacent cells; key them by edge identity so a
    // crossing is emitted once and segments stay topologically connected.
    let mut vertex_of_edge: HashMap<[usize; 4], usize> = HashMap::new();

    for r in 0..nrows.saturating_sub(1) {
        for c in 0..ncols.saturating_sub(1) {
            let corners = [(r, c), (r, c + 1), (r + 1, c + 1), (r + 1, c)];
            let mut vals = [0.0f64; 4];
            let mut usable = true;
            for (i, &(cr, cc)) in corners.iter().enumerate() {
                let f = signed_distance(cr, cc);
                if !f.is_finite() {
                    usable = false;
                    break;
                }
                // A node exactly on the plane would make the sign test ill-defined;
                // nudging it off by an ulp-scale epsilon keeps the topology consistent.
                vals[i] = if f == 0.0 { f64::MIN_POSITIVE } else { f };
            }
            if !usable {
                continue;
            }

            let mut case = 0u8;
            for (i, &v) in vals.iter().enumerate() {
                if v > 0.0 {
                    case |= 1 << i;
                }
            }

            let mut chords = case_edges(case);
            if chords.is_empty() {
                continue;
            }
            if case == 5 || case == 10 {
                // Saddle: the cell centre tells which of the two pairings is right.
                let centre = vals.iter().sum::<f64>() / 4.0;
                if (centre > 0.0) != (case == 5) {
                    chords = saddle_alternative(case);
                }
            }

            for &(e0, e1) in chords {
                let mut ends = [0usize; 2];
                for (slot, edge) in [e0, e1].iter().enumerate() {
                    let (a, b) = edge.corners();
                    let (r0, c0) = corners[a];
                    let (r1, c1) = corners[b];
                    let key = [r0, c0, r1, c1];
                    ends[slot] = match vertex_of_edge.get(&key) {
                        Some(&idx) => idx,
                        None => {
                            let (xa, ya) = geotransform.node_xy(r0, c0);
                            let (xb, yb) = geotransform.node_xy(r1, c1);
                            let (za, zb) = (elevation(r0, c0), elevation(r1, c1));
                            let (fa, fb) = (signed_distance(r0, c0), signed_distance(r1, c1));
                            let t = fa / (fa - fb);
                            let idx = points.len();
                            points.push([
                                xa + t * (xb - xa),
                                ya + t * (yb - ya),
                                za + t * (zb - za),
                            ]);
                            vertex_of_edge.insert(key, idx);
                            idx
                        }
                    };
                }
                segments.push(ends);
            }
        }
    }

    Intersection { points, segments }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn flat_dem(nrows: usize, ncols: usize, z: f64) -> Vec<f64> {
        vec![z; nrows * ncols]
    }

    fn unit_transform() -> GeoTransform {
        GeoTransform {
            x_origin: 0.0,
            pixel_width: 1.0,
            row_rotation: 0.0,
            y_origin: 10.0,
            col_rotation: 0.0,
            pixel_height: -1.0,
        }
    }

    #[test]
    fn horizontal_plane_normal_points_up() {
        let n = plane_normal(0.0, 0.0);
        assert!((n[0]).abs() < 1e-12 && (n[1]).abs() < 1e-12 && (n[2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn vertical_plane_normal_is_horizontal() {
        let n = plane_normal(90.0, 90.0);
        assert!((n[0] - 1.0).abs() < 1e-12 && n[2].abs() < 1e-12);
    }

    #[test]
    fn flat_dem_below_plane_yields_nothing() {
        let dem = flat_dem(10, 10, 0.0);
        let out = intersect_plane_dem(
            &dem, 10, 10, &unit_transform(), [0.0, 0.0, 100.0], 0.0, 0.0, None,
        );
        assert!(out.points.is_empty());
        assert!(out.segments.is_empty());
    }

    #[test]
    fn tilted_plane_cuts_flat_dem_along_a_line() {
        let dem = flat_dem(10, 10, 0.0);
        let out = intersect_plane_dem(
            &dem, 10, 10, &unit_transform(), [5.0, 5.0, 0.0], 90.0, 45.0, None,
        );
        assert!(!out.points.is_empty());
        let n = plane_normal(90.0, 45.0);
        for p in &out.points {
            let d = n[0] * (p[0] - 5.0) + n[1] * (p[1] - 5.0) + n[2] * (p[2] - 0.0);
            assert!(d.abs() < 1e-9, "point off plane: {}", d);
        }
    }

    #[test]
    fn nodata_cells_are_skipped() {
        let mut dem = flat_dem(10, 10, 0.0);
        for v in dem.iter_mut() {
            *v = -9999.0;
        }
        let out = intersect_plane_dem(
            &dem, 10, 10, &unit_transform(), [5.0, 5.0, 0.0], 90.0, 45.0, Some(-9999.0),
        );
        assert!(out.points.is_empty());
    }

    #[test]
    fn segments_index_existing_vertices() {
        let dem = flat_dem(10, 10, 0.0);
        let out = intersect_plane_dem(
            &dem, 10, 10, &unit_transform(), [5.0, 5.0, 0.0], 90.0, 45.0, None,
        );
        for s in &out.segments {
            assert!(s[0] < out.points.len() && s[1] < out.points.len());
        }
    }
}
