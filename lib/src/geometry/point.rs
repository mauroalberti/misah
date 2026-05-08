
use super::math::Vector

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point<const N: usize> {
    pub coords: [f64; N],
}

pub type Point2D = Point<2>;
pub type Point3D = Point<3>;

// conversione
impl<const N: usize> From<[f64; N]> for Point<N> {

    fn from(coords: [f64; N]) -> Self {
        Self { coords }
    }

}

impl<const N: usize> Point<N> {

    pub fn coord(&self, i: usize) -> Option<f64> {
        self.coords.get(i).copied()
    }

    pub fn distance(&self, other: &Self) -> f64 {
        self.coords
            .iter()
            .zip(other.coords.iter())
            .map(|(a, b)| {
                let d = b - a;
                d*d
            })
            .sum::<f64>()
            .sqrt()
    }

    pub fn vector_to(&self, other: &Self) -> Vector<N> {

        let mut coords = [0.0; N];

        for i in 0..N {
            coords[i] = other.coords[i] - self.coords[i];
        }

        Vector { coords }
    }
}

impl Point<1> {
    pub fn x(&self) -> f64 { self.coords[0] }
}

impl Point<2> {
    pub fn x(&self) -> f64 { self.coords[0] }
    pub fn y(&self) -> f64 { self.coords[1] }
}

impl Point<3> {
    pub fn x(&self) -> f64 { self.coords[0] }
    pub fn y(&self) -> f64 { self.coords[1] }
    pub fn z(&self) -> f64 { self.coords[2] }
}



