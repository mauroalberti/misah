
use std::ops::{Add, Neg, Sub, Mul, Div};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vector<const N: usize> {
    pub coords: [f64; N],
}

impl<const N: usize> Vector<N> {

    pub fn new(coords: [f64; N]) -> Self {
        Self { coords }
    }

    pub fn norm(&self) -> f64 {
        self.coords
            .iter()
            .map(|c| c*c)
            .sum::<f64>()
            .sqrt()
    }

    pub fn dot(&self, other: &Self) -> f64 {
        self.coords
            .iter()
            .zip(other.coords.iter())
            .map(|(a, b)| a * b)
            .sum()
    }

    pub fn is_zero(&self) -> bool {
        self.coords
            .iter()
            .all(|c| *c == 0.0)
    }

}

impl<const N: usize> Neg for Vector<N> {

    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            coords: self.coords.map(|c| -c ),
        }
    }
}

impl<const N: usize> Mul<f64> for Vector<N> {

    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Self {
            coords: self.coords.map(|c| c * rhs ),
        }
    }
}

impl<const N: usize> Div<f64> for Vector<N> {

    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Self {
            coords: self.coords.map(|c| c / rhs ),
        }
    }
}

impl<const N: usize> Add for Vector<N> {

    type Output = Self;

    fn add(self, other: Self) -> Self::Output {

        let mut coords = [0.0; N];

        for i in 0..N {
            coords[i] = self.coords[i] + other.coords[i];
        }

        Self { coords }
    }
}

impl<const N: usize> Sub for Vector<N> {

    type Output = Self;

    fn sub(self, other: Self) -> Self::Output {

        let mut coords = [0.0; N];

        for i in 0..N {
            coords[i] = self.coords[i] - other.coords[i];
        }

        Self { coords }
    }
}

