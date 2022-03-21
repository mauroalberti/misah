
// open.rs
// from:
//   - http://rustbyexample.com/std_misc/file/open.html
//   - https://users.rust-lang.org/t/read-a-file-line-by-line/1585/7
//    whitebox-tools

use super::*;

use std::io::prelude::*;
use std::io::BufReader;
use std::error::Error;
use std::fs::File;
use std::path::Path;


pub fn read_ascii_dem(
    path_str: &string,
    epsg_code: i32
) -> Result<GeoArray, Error> {

    // Create a path to the desired file
    let path = Path::new(path_str);

    // Open the path in read-only mode

    let file = match File::open(&path) {
        Err(err) => return Err,
        Ok(file) => file,
    };

    let mut lines = BufReader::new(&file).lines();

    /*
    NCOLS 525
    NROWS 544
    XLLCORNER 317567.8563212453
    YLLCORNER 4730903.4029195514
    CELLSIZE 28.030813
    NODATA_VALUE 170141000918782800000000000000000000000.0000
     */

    // ncols
    let line = match lines.next().unwrap() {
        Err(err) => return Err,
        Ok(line) => line,
    };
    let ncols = line.split_whitespace().nth(1).unwrap().parse::<i32>().unwrap();
    println!("ncols: {}", ncols);

    // nrows
    let line = match lines.next().unwrap() {
        Err(err) => return Err,
        Ok(line) => line,
    };
    let nrows = line.split_whitespace().nth(1).unwrap().parse::<i32>().unwrap();
    println!("nrows: {}", nrows);

    // xllcorner
    let line = match lines.next().unwrap() {
        Err(err) => return Err,
        Ok(line) => line,
    };
    let xllcorner = line.split_whitespace().nth(1).unwrap().parse::<f32>().unwrap();
    println!("xllcorner: {}", xllcorner);

    // yllcorner
    let line = match lines.next().unwrap() {
        Err(err) => return Err,
        Ok(line) => line,
    };
    let yllcorner = line.split_whitespace().nth(1).unwrap().parse::<f32>().unwrap();
    println!("yllcorner: {}", yllcorner);

    // cellsize
    let line = match lines.next().unwrap() {
        Err(err) => return Err,
        Ok(line) => line,
    };
    let cellsize = line.split_whitespace().nth(1).unwrap().parse::<f32>().unwrap();
    println!("cellsize: {}", cellsize);

    // nodata_value
    let line = match lines.next().unwrap() {
        Err(err) => return Err,
        Ok(line) => line,
    };
    let nodata_value = line.split_whitespace().nth(1).unwrap().parse::<f32>().unwrap();
    println!("nodata_value: {}", nodata_value);


}
