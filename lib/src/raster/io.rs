
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
    file_path: &string,
    epsg_code: i32
) -> Result<Grid, Error> {

    // Open the path in read-only mode

    let file = File::open(&file_path)?;
    let lines = BufReader::new(&file).lines();

    // ncols
    let mut line = lines.next()?;
    let ncols = line.split_whitespace().nth(1).unwrap().parse::<i32>().unwrap();
    println!("ncols: {}", ncols);

    // nrows
    line = lines.next()?;
    let nrows = line.split_whitespace().nth(1).unwrap().parse::<i32>().unwrap();
    println!("nrows: {}", nrows);

    // xllcorner
    line = lines.next()?;
    let xllcorner = line.split_whitespace().nth(1).unwrap().parse::<f32>().unwrap();
    println!("xllcorner: {}", xllcorner);

    // yllcorner
    line = lines.next()?;
    let yllcorner = line.split_whitespace().nth(1).unwrap().parse::<f32>().unwrap();
    println!("yllcorner: {}", yllcorner);

    // cellsize
    line = lines.next()?;
    let cellsize = line.split_whitespace().nth(1).unwrap().parse::<f32>().unwrap();
    println!("cellsize: {}", cellsize);

    // nodata_value
    line = lines.next()?;
    let nodata_value = line.split_whitespace().nth(1).unwrap().parse::<f32>().unwrap();
    println!("nodata_value: {}", nodata_value);

    while let line = lines.next() {

    };
}
