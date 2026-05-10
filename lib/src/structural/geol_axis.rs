
pub struct GeologicalAxis {
    pub trend: f64,
    pub plunge: f64
}

impl GeologicalAxis {

    fn new(tr: f64, pl: f64) -> Self { GeologicalAxis{ trend: tr, plunge: pl }}

}
