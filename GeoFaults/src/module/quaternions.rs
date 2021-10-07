
use libc::{c_bool, c_int, c_double};

extern "C"{
    pub fn quat_product_(quat1: * [c_double; 4], quat2: * [c_double; 4]);
    pub fn quat_conjugate_(d: *mut c_double);
    pub fn quat_squarednorm_(d: *mut c_double);
    pub fn quat_scalardivision_(d: *mut c_double);
    pub fn quat_inverse_(d: *mut c_double);
    pub fn quat_normaliztest_();
    pub fn quat_normalization_();
    pub fn quaternfromcartmatr();
}