
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vector<const N: usize> {
    pub coords: [f64; N],
}

impl<const N: usize> Vector<N> {

    pub fn new(coords: [f64; N]) -> Self {
        Self { coords }
    }

    pub fn norm(&self) -> f64 {
        self.coords.iter().map(|c| c*c).sum::<f64>().sqrt()
    }

    pub fn dot(&self, other: &Self) -> f64 {
        self.coords
            .iter()
            .zip(other.coords.iter())
            .map(|(a, b)| a * b)
            .sum()
    }

    pub fn is_zero(&self) -> bool {
        self.coords.iter().all(|c| *c == 0.0)
    }

}
