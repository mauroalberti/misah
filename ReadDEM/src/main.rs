// open.rs
// from: 
//   - http://rustbyexample.com/std_misc/file/open.html
//   - https://users.rust-lang.org/t/read-a-file-line-by-line/1585/7

use std::io::prelude::*;
use std::io::BufReader;
use std::error::Error;
use std::fs::File;
use std::path::Path;


fn main() {

    let path_str = "/home/mauro/Documents/GeoDati/italia/locali/example_data/srtm/E010_N40.asc";

    // Create a path to the desired file
    
    let path = Path::new(path_str);

    // Open the path in read-only mode
    
    let file = match File::open(&path) {
        Err(err) => panic!("couldn't open {}: {}", path_str, err.description()),
        Ok(file) => file,
    };

    let mut lines = BufReader::new(&file).lines();

    let line = match lines.next().unwrap() {
        Err(err) => panic!("couldn't process first line: {}", err.description()),
        Ok(line) => line,
        };
    let ncols = line.split_whitespace().nth(1).unwrap().parse::<i32>().unwrap();

    println!("column number: {}", ncols);

 
}
