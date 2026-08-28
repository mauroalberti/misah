
use ndarray::Array2;
use crate::raster::geotransform::GeoTransform;

// modified from: https://stackoverflow.com/questions/13212212/creating-two-dimensional-arrays-in-rust

#[derive(Debug)]
pub struct Grid {
    pub transform: GeoTransform,
    pub epsg_code: i32,
    pub data: Array2<f64>
}
