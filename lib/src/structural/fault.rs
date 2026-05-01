
use super::geolplane::GeologicalPlane;
use super::slickenline::Slickenline;


#[derive(Debug, Clone)]
pub struct FaultPlane {
    pub plane GeologicalPlane,
    pub slickenlines Vec<Slickenline>,
}
