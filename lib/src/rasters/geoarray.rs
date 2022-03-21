
use super::*;

use ndarray::Array2;
use crate::rasters::geotransform::GeoTransform;

// modified from: https://stackoverflow.com/questions/13212212/creating-two-dimensional-arrays-in-rust

pub struct GeoArray {
    pub transform: GeoTransform,
    pub epsg_code: i32,
    pub data: Array2<f64>
}