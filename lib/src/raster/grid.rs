
use ndarray::Array2;
use crate::raster::geotransform::GeoTransform;

/// A raster: its elevations, and the affine transform placing them on the
/// ground.
///
/// No CRS. There was an `epsg_code` here, threaded through `raster::io` and out
/// to the Python binding as a parameter, and nothing ever read it: it was not
/// validated, not propagated to any output, and did not stop two grids in
/// different systems from being intersected with each other. A field that
/// nothing reads is not metadata, and as an argument in a public signature it
/// advertised a CRS awareness this crate does not have. Callers who need to
/// track the system can do so beside the grid exactly as easily; when the
/// checks are worth writing, they will need a real CRS and a place to apply
/// it, not an integer parked here in advance.
#[derive(Debug)]
pub struct Grid {
    pub transform: GeoTransform,
    pub data: Array2<f64>
}
