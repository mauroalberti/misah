
use libc::{c_bool, c_int, c_double};

extern "C"{
    pub fn quat_product_(quat1: * [c_double; 4], quat2: * [c_double; 4], product: *mut [c_double; 4]);
    pub fn quat_conjugate_(quat: * [c_double; 4], conjugate: *mut [c_double; 4]);
    pub fn quat_squarednorm_(quat: * [c_double; 4], squarednorm: *mut [c_double; 4]);
    pub fn quat_scalardivision_(quat: * [c_double; 4], scaldiv: * c_double, quat_scaldiv: *mut [c_double; 4]);
    pub fn quat_inverse_(quat: * [c_double; 4], conjugate: *mut [c_double; 4]);
    pub fn is_quat_normalized_(quat: * [c_double; 4], is_normalized: *mut c_bool);
    pub fn quat_normalization_(quat: * [c_double; 4], normalized: *mut [c_double; 4]);
    pub fn quaternfromcartmatr(focmec_matrix: * [[c_double; 3]; 3], quat: *mut [c_double; 4]);
}