
pub struct Axis {
    pub trend: f64,
    pub plunge: f64
}

impl Axis {

    fn new(tr: f64, pl: f64) -> Self { Axis{ trend: tr, plunge: pl }}

}
