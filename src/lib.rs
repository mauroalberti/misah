
use pyo3::prelude::*;
use pyo3::{PyResult, Python, py_run};

use crate::geometry::space2d::{Point2D, Segment2D};

use crate::geometry::space3d::{Point3D, Segment3D};

use crate::georeferenced::space2d::raster::{GeoArray, GeoTransform};

use crate::orientations::space3d::{Axis, GeolPlane};

pub mod geometry;
pub mod georeferenced;
pub mod orientations;


#[pymodule]
fn misah(python: Python, module: &PyModule) -> PyResult<()> {

    let geometry = PyModule::new(python, "geometry")?;
    py_run!(python, geometry, "import sys; sys.modules['misah.geometry'] = geometry");
    module.add_submodule(geometry)?;

    let geometry_space2d = PyModule::new(python, "geometry.space2d")?;
    py_run!(python, geometry_space2d, "import sys; sys.modules['misah.geometry.space2d'] = geometry.space2d");
    geometry.add_submodule(geometry_space2d)?;

    geometry_space2d.add_class::<Point2D>()?;
    geometry_space2d.add_class::<Segment2D>()?;

    let geometry_space3d = PyModule::new(python, "geometry.space3d")?;
    py_run!(python, geometry_space3d, "import sys; sys.modules['misah.geometry.space3d'] = geometry.space3d");
    geometry.add_submodule(geometry_space3d)?;

    geometry_space3d.add_class::<Point3D>()?;
    geometry_space3d.add_class::<Segment3D>()?;

    let georeferenced = PyModule::new(python, "georeferenced")?;
    py_run!(python, georeferenced, "import sys; sys.modules['misah.georeferenced'] = georeferenced");
    module.add_submodule(georeferenced)?;

    let georeferenced_space2d = PyModule::new(python, "georeferenced.space2d")?;
    py_run!(python, georeferenced_space2d, "import sys; sys.modules['misah.georeferenced.space2d'] = georeferenced.space2d");
    georeferenced.add_submodule(georeferenced_space2d)?;

    let georeferenced_raster = PyModule::new(python, "georeferenced.space2d.raster")?;
    py_run!(python, georeferenced_raster, "import sys; sys.modules['misah.georeferenced.space2d.raster'] = georeferenced.space2d.raster");
    georeferenced.add_submodule(georeferenced_raster)?;

    georeferenced_raster.add_class::<GeoTransform>()?;

    let orientations = PyModule::new(python, "orientations")?;
    py_run!(python, orientations, "import sys; sys.modules['misah.orientations'] = orientations");
    module.add_submodule(orientations)?;

    let orientations_space3d = PyModule::new(python, "orientations.space3d")?;
    py_run!(python, orientations_space3d, "import sys; sys.modules['misah.orientations.space3d'] = orientations.space3d");
    orientations.add_submodule(orientations_space3d)?;

    //orientations_space3d.add_class::<Axis>()?;
    orientations_space3d.add_class::<GeolPlane>()?;

    Ok(())
}


