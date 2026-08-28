
/// GDAL-style affine transform: origin, then the two rows of the matrix, as
/// `[x_origin, pixel_width, row_rotation, y_origin, col_rotation, pixel_height]`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeoTransform {
    pub data: [f64; 6]
}

impl GeoTransform {

    pub fn x_origin(&self) -> f64 { self.data[0] }
    pub fn pixel_width(&self) -> f64 { self.data[1] }
    pub fn row_rotation(&self) -> f64 { self.data[2] }
    pub fn y_origin(&self) -> f64 { self.data[3] }
    pub fn col_rotation(&self) -> f64 { self.data[4] }
    pub fn pixel_height(&self) -> f64 { self.data[5] }

    /// Map a cell index to the ground coordinates of its centre.
    ///
    /// Raster values are point samples taken at cell centres, hence the half-cell
    /// offsets; passing the raw indices would shift every result by half a cell.
    pub fn node_xy(&self, row: usize, col: usize) -> (f64, f64) {

        let (r, c) = (row as f64 + 0.5, col as f64 + 0.5);

        (
            self.data[0] + c * self.data[1] + r * self.data[2],
            self.data[3] + c * self.data[4] + r * self.data[5],
        )
    }

    /// Map ground coordinates back to a fractional cell index: the inverse of
    /// `node_xy`, and so likewise measured to cell centres, meaning
    /// `node_rc(node_xy(r, c))` returns `(r, c)` and a point halfway between
    /// two node centres comes back as a half.
    ///
    /// The result is neither clamped nor rounded: a point outside the raster
    /// yields indices outside its extent, negative ones included. That is what
    /// lets a caller tell "outside" from "on the edge", instead of having the
    /// distinction quietly clamped away.
    ///
    /// `None` when the transform is singular -- a raster with a zero-sized or
    /// collapsed cell -- since then there is no inverse to return.
    pub fn node_rc(&self, x: f64, y: f64) -> Option<(f64, f64)> {

        let [x0, a, b, y0, d, e] = self.data;

        let det = a * e - b * d;

        if det == 0.0 {
            return None;
        }

        let (dx, dy) = (x - x0, y - y0);

        // Solved for the half-shifted indices node_xy works in, hence the half
        // taken back off each at the end.
        let col = (dx * e - dy * b) / det;
        let row = (dy * a - dx * d) / det;

        Some((row - 0.5, col - 0.5))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn node_xy_returns_cell_centres() {
        // North-up grid of unit cells, origin at the top-left corner.
        let gt = GeoTransform { data: [0.0, 1.0, 0.0, 10.0, 0.0, -1.0] };

        assert_eq!(gt.node_xy(0, 0), (0.5, 9.5));
        assert_eq!(gt.node_xy(1, 2), (2.5, 8.5));
    }

    #[test]
    fn node_rc_inverts_node_xy() {
        let gt = GeoTransform { data: [0.0, 1.0, 0.0, 10.0, 0.0, -1.0] };

        for &(r, c) in &[(0usize, 0usize), (1, 2), (7, 3)] {
            let (x, y) = gt.node_xy(r, c);
            let (rr, cc) = gt.node_rc(x, y).expect("a north-up transform is invertible");

            assert!((rr - r as f64).abs() < 1e-12);
            assert!((cc - c as f64).abs() < 1e-12);
        }
    }

    #[test]
    fn node_rc_inverts_a_rotated_transform() {
        // 30 degrees off north, and cells that are not square, so nothing can
        // pass by accident of symmetry.
        let (s, k) = (30.0_f64.to_radians().sin(), 30.0_f64.to_radians().cos());
        let gt = GeoTransform {
            data: [1000.0, 2.0 * k, 3.0 * s, 5000.0, 2.0 * s, -3.0 * k],
        };

        let (x, y) = gt.node_xy(11, 4);
        let (rr, cc) = gt.node_rc(x, y).expect("a rotated transform is still invertible");

        assert!((rr - 11.0).abs() < 1e-9);
        assert!((cc - 4.0).abs() < 1e-9);
    }

    #[test]
    fn node_rc_reports_positions_outside_the_raster() {
        let gt = GeoTransform { data: [0.0, 1.0, 0.0, 10.0, 0.0, -1.0] };

        let (row, col) = gt.node_rc(-5.5, 15.5).expect("invertible");

        assert!(row < 0.0 && col < 0.0);
    }

    #[test]
    fn node_rc_of_a_singular_transform_is_none() {
        // Zero-width cells: the transform collapses the plane onto a line.
        let gt = GeoTransform { data: [0.0, 0.0, 0.0, 10.0, 0.0, -1.0] };

        assert!(gt.node_rc(1.0, 1.0).is_none());
    }
}
