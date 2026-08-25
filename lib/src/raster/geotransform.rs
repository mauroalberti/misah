
/// GDAL-style affine transform: origin, then the two rows of the matrix, as
/// `[x_origin, pixel_width, row_rotation, y_origin, col_rotation, pixel_height]`.
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
}
