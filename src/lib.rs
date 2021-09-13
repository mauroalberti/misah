
use pyo3::prelude::*;
use pyo3::{PyResult, Python, py_run};

use crate::geom2d::Point2D;
use crate::geom2d::Segment2D;

use crate::geom3d::Point3D;
use crate::geom3d::Segment3D;

use crate::orientations::Axis;
use crate::orientations::GeolPlane;

pub mod geom2d;
pub mod geom3d;
pub mod orientations;


#[pymodule]
fn misah(python: Python, module: &PyModule) -> PyResult<()> {

    let geom2d = PyModule::new(python, "misah.geom2d")?;
    py_run!(python, geom2d, "import sys; sys.modules['misah.geom2d'] = geom2d");
    module.add_submodule(geom2d)?;

    geom2d.add_class::<Point2D>()?;
    geom2d.add_class::<Segment2D>()?;

    let geom3d = PyModule::new(python, "geom3d")?;
    py_run!(python, geom3d, "import sys; sys.modules['misah.geom3d'] = geom3d");
    module.add_submodule(geom3d)?;

    geom3d.add_class::<Point3D>()?;
    geom3d.add_class::<Segment3D>()?;

    let orientations = PyModule::new(python, "orientations")?;
    py_run!(python, orientations, "import sys; sys.modules['misah.orientations'] = orientations");
    module.add_submodule(orientations)?;

    orientations.add_class::<Axis>()?;
    orientations.add_class::<GeolPlane>()?;

    Ok(())
}


