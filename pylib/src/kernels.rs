//! Python bindings for the numerical kernels.

use numpy::{IntoPyArray, PyArray1, PyArray2, PyArrayMethods, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use misah::geometry::mesh::TriangleMesh;
use misah::geometry::point::Point3D;
use misah::raster::geotransform::GeoTransform;
use misah::raster::grid::Grid;
use misah::raster::intersection::intersect_plane_grid as kernel;
use misah::raster::mesh_intersection::intersect_mesh_grid as mesh_kernel;
use misah::structural::best_fit::best_fit_geoplanes;
use misah::structural::fault::FaultPlane;
use misah::structural::geol_axis::GeologicalAxis;
use misah::structural::geol_plane::GeologicalPlane;
use misah::structural::inversion::{invert, invert_weighted, ScoredTensor, SearchGrid};
use misah::structural::slickenline::{SlipSense, Slickenline};
use misah::structural::stress::ReducedStressTensor;

/// What a grid intersection hands back to Python: the vertices, as an (N, 3)
/// array of coordinates, and the chords, as an (M, 2) array of indices into
/// them.
///
/// Named rather than spelled out at the signature, where four nested generics
/// take longer to read than they say.
type VerticesAndSegments<'py> = (Bound<'py, PyArray2<f64>>, Bound<'py, PyArray2<i64>>);

/// What a mesh-grid intersection hands back: the points as an (N, 3) array of
/// coordinates, the attitude of the mesh triangle that produced each as an
/// (N, 2) array of dip direction and dip angle, the index of that triangle as
/// an (N,) array, and the run's own statistics as a dict.
type MeshIntersections<'py> = (
    Bound<'py, PyArray2<f64>>,
    Bound<'py, PyArray2<f64>>,
    Bound<'py, PyArray1<i64>>,
    Bound<'py, PyDict>,
);

/// Intersect an unbounded geological plane with a DEM.
///
/// `geotransform` is the GDAL six-element affine transform. Returns
/// `(points, segments)`: an (N, 3) array of intersection vertices and an (M, 2)
/// array of indices into it, one row per marching-squares chord.
#[pyfunction]
#[pyo3(signature = (
    dem,
    geotransform,
    src_pt,
    dip_dir_degr,
    dip_angle_degr,
    nodata = None,
))]
#[allow(clippy::too_many_arguments)]
fn intersect_plane_grid<'py>(
    py: Python<'py>,
    dem: PyReadonlyArray2<'py, f64>,
    geotransform: [f64; 6],
    src_pt: [f64; 3],
    dip_dir_degr: f64,
    dip_angle_degr: f64,
    nodata: Option<f64>,
) -> PyResult<VerticesAndSegments<'py>> {
    let view = dem.as_array();
    if !view.is_standard_layout() {
        return Err(PyValueError::new_err(
            "DEM must be C-contiguous; pass numpy.ascontiguousarray(dem)",
        ));
    }

    let plane = GeologicalPlane::new(dip_dir_degr, dip_angle_degr)
        .to_plane(Point3D::from(src_pt))
        .ok_or_else(|| PyValueError::new_err("degenerate plane normal"))?;

    // Grid owns its elevations, so the array is copied here. For a DEM that is a
    // few megabytes at most, and it keeps the kernel free of lifetimes.
    let grid = Grid {
        transform: GeoTransform { data: geotransform },
        data: view.to_owned(),
    };

    let out = py.detach(|| kernel(&plane, &grid, nodata));

    let n_pts = out.points.len();
    let flat_pts: Vec<f64> = out.points.iter().flat_map(|p| p.coords).collect();
    let points = flat_pts.into_pyarray(py).reshape([n_pts, 3])?;

    let n_segs = out.segments.len();
    let flat_segs: Vec<i64> = out
        .segments
        .iter()
        .flat_map(|s| [s[0] as i64, s[1] as i64])
        .collect();
    let segments = flat_segs.into_pyarray(py).reshape([n_segs, 2])?;

    Ok((points, segments))
}

/// Upward-pointing unit normal of a geological plane, in (East, North, Up).
#[pyfunction]
fn plane_normal(dip_dir_degr: f64, dip_angle_degr: f64) -> [f64; 3] {
    GeologicalPlane::new(dip_dir_degr, dip_angle_degr)
        .normal_vector()
        .coords
}

/// The trend/plunge of a slickenline of the given rake, on a plane of the
/// given strike (right-hand rule) and dip.
///
/// Aki & Richards (1980)'s convention: rake 0 is left-lateral, 90 reverse,
/// +/-180 right-lateral, -90 normal.
#[pyfunction]
fn rake_to_slickenline(strike_rhr_degr: f64, dip_angle_degr: f64, rake_degr: f64) -> (f64, f64) {
    let plane = GeologicalPlane::from_rhr_strike(strike_rhr_degr, dip_angle_degr);
    let axis = GeologicalAxis::from_versor(&plane.rake_to_versor(rake_degr));
    (axis.trend, axis.plunge)
}

/// The direct (forward) Wallace-Bott problem: resolve a reduced stress
/// tensor onto one fault plane, predicting the slip it drives.
///
/// `s1`/`s3` are the principal stress axes as (trend, plunge) in degrees,
/// sub-orthogonal to within a degree; `phi` is the shape ratio
/// `(sigma2 - sigma3) / (sigma1 - sigma3)`. `sigma1`/`sigma3` default to 1/0,
/// the usual normalization when only the tensor's shape is known, as from a
/// fault-slip inversion -- which leaves the predicted rake correct while
/// making `slip_tendency`/`deformation_index` meaningless; pass the true
/// magnitudes when they are known and those are wanted.
///
/// Returns a dict always carrying `is_valid`, `traction`, `traction_magnitude`,
/// `normal_stress`, `normal_stress_magnitude`, `shear_stress` and
/// `shear_stress_magnitude` (the last three vectors in (East, North, Up)),
/// plus `theoretical_rake`, `theoretical_slickenline` (a `(trend, plunge)`
/// pair), `slip_tendency` and `deformation_index`, `None` on all four where
/// the shear stress does not clear `shear_threshold` -- the plane sits on, or
/// acutely close to, a principal stress axis, and no slip direction is
/// defined.
#[pyfunction]
#[pyo3(signature = (
    s1_trend_degr,
    s1_plunge_degr,
    s3_trend_degr,
    s3_plunge_degr,
    phi,
    strike_rhr_degr,
    dip_angle_degr,
    sigma1 = 1.0,
    sigma3 = 0.0,
    shear_threshold = misah::structural::stress::SHEAR_MAGNITUDE_THRESHOLD,
))]
#[allow(clippy::too_many_arguments)]
fn solve_stress<'py>(
    py: Python<'py>,
    s1_trend_degr: f64,
    s1_plunge_degr: f64,
    s3_trend_degr: f64,
    s3_plunge_degr: f64,
    phi: f64,
    strike_rhr_degr: f64,
    dip_angle_degr: f64,
    sigma1: f64,
    sigma3: f64,
    shear_threshold: f64,
) -> PyResult<Bound<'py, PyDict>> {

    let s1 = GeologicalAxis::new(s1_trend_degr, s1_plunge_degr);
    let s3 = GeologicalAxis::new(s3_trend_degr, s3_plunge_degr);

    let tensor = ReducedStressTensor::new(s1, s3, phi, sigma1, sigma3)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;

    let plane = GeologicalPlane::from_rhr_strike(strike_rhr_degr, dip_angle_degr);
    let solution = tensor.solve_with_threshold(&plane, shear_threshold);

    let slickenline = solution.theoretical_slickenline.map(|v| {
        let axis = GeologicalAxis::from_versor(&v);
        (axis.trend, axis.plunge)
    });

    let dict = PyDict::new(py);
    dict.set_item("is_valid", solution.is_valid)?;
    dict.set_item("traction", (solution.traction[0], solution.traction[1], solution.traction[2]))?;
    dict.set_item("traction_magnitude", solution.traction_magnitude)?;
    dict.set_item("normal_stress", (solution.normal_stress[0], solution.normal_stress[1], solution.normal_stress[2]))?;
    dict.set_item("normal_stress_magnitude", solution.normal_stress_magnitude)?;
    dict.set_item("shear_stress", (solution.shear_stress[0], solution.shear_stress[1], solution.shear_stress[2]))?;
    dict.set_item("shear_stress_magnitude", solution.shear_stress_magnitude)?;
    dict.set_item("theoretical_rake", solution.theoretical_rake)?;
    dict.set_item("theoretical_slickenline", slickenline)?;
    dict.set_item("slip_tendency", solution.slip_tendency)?;
    dict.set_item("deformation_index", solution.deformation_index)?;

    Ok(dict)
}

/// Intersect a triangulated surface with a DEM.
///
/// The general case of `intersect_plane_grid`: a surface of arbitrary shape
/// rather than one unbounded plane, which is what a folded or faulted
/// geological surface needs. `vertices` is a (V, 3) array of coordinates and
/// `faces` a (F, 3) array of indices into it -- the shape a VTK `POLYDATA`
/// file gives, reading which is left to the caller.
///
/// Returns `(points, attitudes, mesh_triangles, stats)`: the intersection
/// points as (N, 3), the dip direction and dip angle of the mesh triangle that
/// produced each as (N, 2), the index of that triangle as (N,), and a dict of
/// counters describing the run. The attitudes are what make the result a set
/// of located measurements rather than a bare trace.
#[pyfunction]
#[pyo3(signature = (dem, geotransform, vertices, faces, nodata = None))]
fn intersect_mesh_grid<'py>(
    py: Python<'py>,
    dem: PyReadonlyArray2<'py, f64>,
    geotransform: [f64; 6],
    vertices: PyReadonlyArray2<'py, f64>,
    faces: PyReadonlyArray2<'py, i64>,
    nodata: Option<f64>,
) -> PyResult<MeshIntersections<'py>> {

    let dem_view = dem.as_array();
    if !dem_view.is_standard_layout() {
        return Err(PyValueError::new_err(
            "DEM must be C-contiguous; pass numpy.ascontiguousarray(dem)",
        ));
    }

    let vertex_view = vertices.as_array();
    if vertex_view.ncols() != 3 {
        return Err(PyValueError::new_err("vertices must be an (N, 3) array"));
    }
    let face_view = faces.as_array();
    if face_view.ncols() != 3 {
        return Err(PyValueError::new_err("faces must be an (N, 3) array"));
    }

    let points: Vec<Point3D> = vertex_view
        .rows()
        .into_iter()
        .map(|r| Point3D::from([r[0], r[1], r[2]]))
        .collect();

    // Checked here rather than left to TriangleMesh's own error, so that a
    // negative index -- which numpy allows and a usize cast would turn into an
    // enormous positive one -- is refused as the mistake it is.
    let mut triples: Vec<[usize; 3]> = Vec::with_capacity(face_view.nrows());
    for row in face_view.rows() {
        let mut triple = [0usize; 3];
        for (slot, &index) in triple.iter_mut().zip(row.iter()) {
            if index < 0 {
                return Err(PyValueError::new_err(format!(
                    "face index {index} is negative"
                )));
            }
            *slot = index as usize;
        }
        triples.push(triple);
    }

    let mesh = TriangleMesh::new(points, triples)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;

    let grid = Grid {
        transform: GeoTransform { data: geotransform },
        data: dem_view.to_owned(),
    };

    let out = py.detach(|| mesh_kernel(&mesh, &grid, nodata));

    let n = out.intersections.len();

    let flat_points: Vec<f64> = out
        .intersections
        .iter()
        .flat_map(|i| i.point.coords)
        .collect();
    let points = flat_points.into_pyarray(py).reshape([n, 3])?;

    let flat_attitudes: Vec<f64> = out
        .intersections
        .iter()
        .flat_map(|i| [i.attitude.azimuth, i.attitude.dip_angle])
        .collect();
    let attitudes = flat_attitudes.into_pyarray(py).reshape([n, 2])?;

    let triangles: Vec<i64> = out
        .intersections
        .iter()
        .map(|i| i.mesh_triangle as i64)
        .collect();
    let mesh_triangles = triangles.into_pyarray(py);

    let s = out.stats;
    let stats = PyDict::new(py);
    stats.set_item("mesh_triangles", s.mesh_triangles)?;
    stats.set_item("degenerate_mesh_triangles", s.degenerate_mesh_triangles)?;
    stats.set_item("mesh_triangles_outside_grid", s.mesh_triangles_outside_grid)?;
    stats.set_item("dem_triangle_pairs", s.dem_triangle_pairs)?;
    stats.set_item("coplanar_sides", s.coplanar_sides)?;
    stats.set_item("duplicate_crossings", s.duplicate_crossings)?;

    Ok((points, attitudes, mesh_triangles, stats))
}

/// Build the faults an inversion consumes out of the plain rows Python holds
/// them in.
///
/// Row by row rather than by iterator, so that a refusal names the fault it
/// came from: `FaultPlane`'s own error indexes the slickenline within its
/// fault, which is always 0 here and would leave a caller with two hundred
/// rows and no idea which one to look at.
fn faults_from_rows(
    faults: &PyReadonlyArray2<'_, f64>,
    senses: Option<&PyReadonlyArray1<'_, bool>>,
) -> PyResult<Vec<FaultPlane>> {

    let view = faults.as_array();
    if view.ncols() != 4 {
        return Err(PyValueError::new_err(
            "faults must be an (N, 4) array of strike, dip, slickenline trend, slickenline plunge",
        ));
    }

    let sense_view = senses.map(|s| s.as_array());
    if let Some(ref s) = sense_view {
        if s.len() != view.nrows() {
            return Err(PyValueError::new_err(format!(
                "senses has {} entries for {} faults",
                s.len(),
                view.nrows()
            )));
        }
    }

    let mut out = Vec::with_capacity(view.nrows());

    for (index, row) in view.rows().into_iter().enumerate() {

        let plane = GeologicalPlane::from_rhr_strike(row[0], row[1]);

        // Only whether the sense is known ever enters the arithmetic -- the
        // trend and plunge already say which way the hanging wall moved, and
        // `Slickenline::angle_to` reads nothing from the variant but its
        // presence. So the boundary carries a flag, and `Down` stands in for
        // "determined" without claiming to be the sense that was read.
        let sense = match sense_view {
            Some(ref s) if !s[index] => None,
            _ => Some(SlipSense::Down),
        };

        let slickenline = Slickenline::from_axis(GeologicalAxis::new(row[2], row[3]), sense);

        let fault = FaultPlane::new(plane, vec![slickenline])
            .map_err(|e| PyValueError::new_err(format!("fault {index}: {e}")))?;

        out.push(fault);
    }

    Ok(out)
}

/// The search grid, once its steps are known to be usable.
///
/// Both loops in the search walk by adding their step, so a zero would spin
/// forever rather than return nothing, and a NaN would leave every comparison
/// false and end the loop before its first iteration. Neither is a grid, and
/// both are worth refusing at the boundary where the number arrives.
fn checked_grid(angle_step_degrees: f64, phi_step: f64) -> PyResult<SearchGrid> {

    let usable = |step: f64| step.is_finite() && step > 0.0;

    if !usable(angle_step_degrees) || !usable(phi_step) {
        return Err(PyValueError::new_err(
            "angle_step_degrees and phi_step must be positive finite numbers",
        ));
    }

    Ok(SearchGrid { angle_step_degrees, phi_step })
}

/// The per-fault weights, once they are known to pair with the faults.
///
/// The core answers a mismatched length or a negative weight with `None`,
/// which is right inside Rust and wrong to pass on here: at this boundary it
/// would arrive as "no solution", and read as something the data did rather
/// than something the call got wrong. Refused with a message instead, naming
/// the offending entry.
///
/// Copied out of the array as well as checked, so that the search can run with
/// the interpreter released rather than holding a borrow of a Python object
/// across it.
fn checked_weights(
    weights: Option<&PyReadonlyArray1<'_, f64>>,
    fault_count: usize,
) -> PyResult<Option<Vec<f64>>> {

    let Some(weights) = weights else {
        return Ok(None);
    };

    let view = weights.as_array();

    if view.len() != fault_count {
        return Err(PyValueError::new_err(format!(
            "weights has {} entries but {} faults were given",
            view.len(),
            fault_count
        )));
    }

    for (index, weight) in view.iter().enumerate() {
        if weight.is_nan() || *weight < 0.0 {
            return Err(PyValueError::new_err(format!(
                "weight {index} is {weight}: weights must be non-negative numbers"
            )));
        }
    }

    Ok(Some(view.to_vec()))
}

/// One scored candidate, as a dict.
fn scored_tensor_dict<'py>(
    py: Python<'py>,
    scored: &ScoredTensor,
) -> PyResult<Bound<'py, PyDict>> {

    let axis = |v| {
        let a = GeologicalAxis::from_versor(&v);
        (a.trend, a.plunge)
    };

    let dict = PyDict::new(py);
    dict.set_item("s1", axis(scored.tensor.s1_versor()))?;
    dict.set_item("s2", axis(scored.tensor.s2_versor()))?;
    dict.set_item("s3", axis(scored.tensor.s3_versor()))?;
    dict.set_item("phi", scored.tensor.phi)?;
    dict.set_item("mean_misfit_degrees", scored.mean_misfit_degrees)?;
    dict.set_item("faults_scored", scored.faults_scored)?;
    dict.set_item("effective_sample_size", scored.effective_sample_size)?;

    Ok(dict)
}

/// Fit a geological plane in every cell of a regular grid over a set of
/// located points.
///
/// The consumer of `intersect_mesh_grid`: that kernel says where a surface
/// meets the topography, this one reads an attitude back out of those points.
/// The two are not coupled -- the input here is bare `(N, 3)` coordinates, so
/// a trace digitised from a map, or readings taken along an outcrop, feed it
/// just as well.
///
/// `cell_size` sets both the grid and the scale attitudes are averaged over:
/// too small and no cell holds three points, too large and a fold is flattened
/// into one meaningless plane. `coincidence_distance` merges points closer
/// than it, which matters because an intersection kernel emits a point per
/// crossing and a shared edge is crossed twice. `max_collinearity` is the
/// threshold on `-log10(s2/s1)` above which a cell's points are taken as
/// lying along a line, where every plane through that line fits equally well
/// and the attitude that comes back is arbitrary rather than imprecise.
///
/// Returns `(field, stats)`, or `None` where there is nothing to fit. `field`
/// is a dict of arrays, all of length M and aligned with each other, one entry
/// per fitted cell in row-major order: `cell_centres` (M, 2) in map
/// coordinates, `attitudes` (M, 2) as dip direction and dip angle,
/// `centroids` (M, 3), `cells` (M, 2) as row and column, `point_counts` (M,),
/// `singular_values` (M, 3), `collinearity` (M,) and `rms_distance` (M,). A
/// dict rather than a tuple because eight aligned arrays positionally is not
/// something anyone should have to read.
///
/// `stats` carries the counters, including `cells_collinear` -- the cells set
/// aside. On a smooth slope that is most of them, and it is the first number
/// to look at when a field comes back thinner than expected.
#[pyfunction]
#[pyo3(signature = (
    points,
    cell_size,
    coincidence_distance = 0.1,
    max_collinearity = misah::structural::best_fit::DEFAULT_MAX_COLLINEARITY,
))]
fn best_fit_planes<'py>(
    py: Python<'py>,
    points: PyReadonlyArray2<'py, f64>,
    cell_size: f64,
    coincidence_distance: f64,
    max_collinearity: f64,
) -> PyResult<Option<(Bound<'py, PyDict>, Bound<'py, PyDict>)>> {

    let view = points.as_array();
    if view.ncols() != 3 {
        return Err(PyValueError::new_err("points must be an (N, 3) array"));
    }

    let coordinates: Vec<Point3D> = view
        .rows()
        .into_iter()
        .map(|r| Point3D::from([r[0], r[1], r[2]]))
        .collect();

    let Some(result) = py.detach(|| {
        best_fit_geoplanes(&coordinates, cell_size, coincidence_distance, max_collinearity)
    }) else {
        // `best_fit_geoplanes` declines a non-positive cell size or
        // coincidence distance the same way it declines an empty input. Only
        // the first is a mistake worth naming.
        if !(cell_size.is_finite() && cell_size > 0.0)
            || !(coincidence_distance.is_finite() && coincidence_distance > 0.0)
        {
            return Err(PyValueError::new_err(
                "cell_size and coincidence_distance must be positive finite numbers",
            ));
        }
        return Ok(None);
    };

    let n = result.fits.len();

    let flat = |values: Vec<f64>, width: usize| -> PyResult<Bound<'py, PyArray2<f64>>> {
        values.into_pyarray(py).reshape([n, width])
    };

    let field = PyDict::new(py);
    field.set_item(
        "cell_centres",
        flat(
            result.fits.iter().flat_map(|f| [f.cell_centre.0, f.cell_centre.1]).collect(),
            2,
        )?,
    )?;
    field.set_item(
        "attitudes",
        flat(
            result
                .fits
                .iter()
                .flat_map(|f| [f.fit.plane.azimuth, f.fit.plane.dip_angle])
                .collect(),
            2,
        )?,
    )?;
    field.set_item(
        "centroids",
        flat(result.fits.iter().flat_map(|f| f.fit.centre.coords).collect(), 3)?,
    )?;
    field.set_item(
        "singular_values",
        flat(result.fits.iter().flat_map(|f| f.fit.singular_values).collect(), 3)?,
    )?;
    field.set_item(
        "cells",
        result
            .fits
            .iter()
            .flat_map(|f| [f.row as i64, f.column as i64])
            .collect::<Vec<i64>>()
            .into_pyarray(py)
            .reshape([n, 2])?,
    )?;
    field.set_item(
        "point_counts",
        result.fits.iter().map(|f| f.fit.points as i64).collect::<Vec<i64>>().into_pyarray(py),
    )?;
    field.set_item(
        "collinearity",
        result.fits.iter().map(|f| f.fit.collinearity()).collect::<Vec<f64>>().into_pyarray(py),
    )?;
    field.set_item(
        "rms_distance",
        result.fits.iter().map(|f| f.fit.rms_distance()).collect::<Vec<f64>>().into_pyarray(py),
    )?;

    let s = result.stats;
    let stats = PyDict::new(py);
    stats.set_item("input_points", s.input_points)?;
    stats.set_item("distinct_points", s.distinct_points)?;
    stats.set_item("non_empty_cells", s.non_empty_cells)?;
    stats.set_item("cells_fitted", s.cells_fitted)?;
    stats.set_item("cells_too_few_points", s.cells_too_few_points)?;
    stats.set_item("cells_collinear", s.cells_collinear)?;
    stats.set_item("grid_rows", result.grid.rows)?;
    stats.set_item("grid_columns", result.grid.columns)?;
    stats.set_item("grid_x_min", result.grid.x_min)?;
    stats.set_item("grid_y_max", result.grid.y_max)?;
    stats.set_item("cell_size", result.grid.cell_size)?;

    Ok(Some((field, stats)))
}

/// The reduced stress tensor itself, as a 3x3 array in (East, North, Up).
///
/// `R . diag(sigma1, sigma2, sigma3) . R^T`, with `R`'s columns the principal
/// axes and sigma2 following from the other two magnitudes and Phi. The
/// axes-and-Phi form is what a fault-slip inversion works in and what
/// `invert_stress` returns, so this is the bridge to anything that wants the
/// matrix instead -- a plot, an export, a rotation into another frame -- and
/// exists because rebuilding it in Python means repeating the construction of
/// the right-handed triad by hand, S3 cross S1 and not the other way, which is
/// a sign error waiting to happen.
///
/// `sigma1`/`sigma3` default to 1/0, matching `solve_stress` and the
/// normalization `invert_stress` searches in. Unlike the misfit, the matrix
/// does scale with them: pass the true magnitudes when they are known, or the
/// result is the shape of the tensor rather than the tensor.
#[pyfunction]
#[pyo3(signature = (
    s1_trend_degr,
    s1_plunge_degr,
    s3_trend_degr,
    s3_plunge_degr,
    phi,
    sigma1 = 1.0,
    sigma3 = 0.0,
))]
#[allow(clippy::too_many_arguments)]
fn stress_tensor<'py>(
    py: Python<'py>,
    s1_trend_degr: f64,
    s1_plunge_degr: f64,
    s3_trend_degr: f64,
    s3_plunge_degr: f64,
    phi: f64,
    sigma1: f64,
    sigma3: f64,
) -> PyResult<Bound<'py, PyArray2<f64>>> {

    let tensor = ReducedStressTensor::new(
        GeologicalAxis::new(s1_trend_degr, s1_plunge_degr),
        GeologicalAxis::new(s3_trend_degr, s3_plunge_degr),
        phi,
        sigma1,
        sigma3,
    )
    .map_err(|e| PyValueError::new_err(e.to_string()))?;

    let flat: Vec<f64> = tensor.tensor().into_iter().flatten().collect();

    flat.into_pyarray(py).reshape([3, 3])
}

/// Score one given tensor against a set of faults, without searching for it.
///
/// The quantity `invert_stress` minimises, evaluated once: how far the slip
/// this tensor predicts falls from the slip observed. Useful where the tensor
/// is not being looked for but tested -- a published one against your own
/// faults, or a runner-up from a search against the subset you suspect belongs
/// to a second tectonic phase.
///
/// `faults` and `senses` are exactly as `invert_stress` takes them, and the
/// same coplanarity check applies, naming the row it refuses.
///
/// There are deliberately no `sigma1`/`sigma3` here. The predicted slip
/// direction is unchanged by them -- scaling a tensor and adding a multiple of
/// the identity both leave the shear direction where it was, the second
/// contributing a purely normal traction -- so a misfit depends only on the
/// tensor's shape. Offering magnitudes would invite the reading that they
/// change the answer.
///
/// Returns a dict of `mean_misfit_degrees` (`None` when no fault could be
/// scored at all), the `faults_scored` that mean came from, and `misfits`, an
/// (N,) array in the order the faults were given. A fault the tensor cannot
/// speak about -- its plane on a principal stress axis, where no slip
/// direction is predicted -- comes back as `NaN` rather than as zero or 180:
/// it is not a fault the tensor fits badly, and it is left out of the mean for
/// the same reason.
#[pyfunction]
#[pyo3(signature = (
    s1_trend_degr,
    s1_plunge_degr,
    s3_trend_degr,
    s3_plunge_degr,
    phi,
    faults,
    senses = None,
))]
#[allow(clippy::too_many_arguments)]
fn score_stress<'py>(
    py: Python<'py>,
    s1_trend_degr: f64,
    s1_plunge_degr: f64,
    s3_trend_degr: f64,
    s3_plunge_degr: f64,
    phi: f64,
    faults: PyReadonlyArray2<'py, f64>,
    senses: Option<PyReadonlyArray1<'py, bool>>,
) -> PyResult<Bound<'py, PyDict>> {

    let tensor = ReducedStressTensor::normalized(
        GeologicalAxis::new(s1_trend_degr, s1_plunge_degr),
        GeologicalAxis::new(s3_trend_degr, s3_plunge_degr),
        phi,
    )
    .map_err(|e| PyValueError::new_err(e.to_string()))?;

    let planes = faults_from_rows(&faults, senses.as_ref())?;

    // `mean_misfit` unrolled, so that the per-fault values it averages over
    // can come back too rather than being computed twice. The same
    // `misfit_on` it calls, and the same arithmetic on the result.
    let misfits: Vec<f64> = planes
        .iter()
        .map(|f| tensor.misfit_on(f).unwrap_or(f64::NAN))
        .collect();

    let scored: Vec<f64> = misfits.iter().copied().filter(|m| !m.is_nan()).collect();
    let mean = (!scored.is_empty()).then(|| scored.iter().sum::<f64>() / scored.len() as f64);

    let dict = PyDict::new(py);
    dict.set_item("mean_misfit_degrees", mean)?;
    dict.set_item("faults_scored", scored.len())?;
    dict.set_item("misfits", misfits.into_pyarray(py))?;

    Ok(dict)
}

/// The inverse Wallace-Bott problem: recover the reduced stress tensor that
/// best explains a set of faults and the slip observed on them.
///
/// `solve_stress` run backwards, and the reason the forward model was worth
/// compiling: the search resolves one forward solution per candidate tensor
/// per fault, which on the default grid over a hundred faults is a couple of
/// million of them.
///
/// `faults` is an (N, 4) array -- fault strike (right-hand rule) and dip
/// angle, then the trend and plunge of the slickenline on it. Each lineation
/// must lie in the plane it is recorded on to within a degree, which is what a
/// pair of compass readings of the same fault is worth; a row that misses by
/// more is refused by name rather than inverted, since a lineation off its
/// plane is not a fault with a small error but two measurements that do not
/// belong to each other.
///
/// `senses` is an optional (N,) boolean array saying, per fault, whether the
/// sense of movement was determined. Where it was not, the slip is known only
/// as a line and the misfit is taken modulo 180 degrees, so a prediction
/// pointing the other way along the same lineation still fits perfectly.
/// Omitting it claims every sense was read, which is the reading that uses the
/// data as given; pass `False` for the faults whose sense nobody could tell,
/// or those faults will be scored at 180 degrees for being right.
///
/// `weights` is an optional (N,) array of non-negative numbers saying how much
/// each fault counts for. Omitting it weights them all alike. The caller this
/// exists for is a stress field on a grid: pass the whole dataset once and
/// recompute only the weights at each node, from a kernel of the distance
/// between the node and where each fault was measured. Faults weighted at zero
/// are dropped before any forward solution is computed, so with a kernel of
/// finite reach a node costs what its own neighbourhood costs rather than what
/// the dataset costs. What the weights mean is not this function's business:
/// they need not sum to anything, and only their relative size is read.
///
/// Weighting is also how two superposed tectonic phases are separated, since
/// nothing else in this search can be told to prefer one over the other.
///
/// The search is exhaustive over the grid rather than a descent, because the
/// misfit surface is not convex: a fault set carrying two superposed tectonic
/// phases has two minima by construction, and a descent would report whichever
/// one it fell into without ever saying the other was there.
///
/// Returns `None` when no candidate could be scored on any fault -- an empty
/// set, or every weight zero, which for a field on a grid means a node with no
/// data within reach and is a result rather than a failure. Otherwise a dict
/// with `best`, `runners_up` (the next five by misfit, worst last, so that a
/// minimum standing alone can be told from one on a plateau),
/// `candidates_tried` and `candidates_scored`. Each scored tensor is itself a
/// dict of `s1`, `s2`, `s3` as `(trend, plunge)` pairs, `phi`,
/// `mean_misfit_degrees`, the `faults_scored` that mean came from -- a
/// candidate lying on a principal axis of half the dataset predicts no slip
/// there and is scored on the rest, so the count is what says whether two
/// misfits are comparable -- and `effective_sample_size`.
///
/// That last is the count in the form that survives weighting: Kish's
/// `(sum w)^2 / sum w^2` over the faults that contributed. Unweighted it is the
/// count exactly. Weighted it is the number to read instead, because a count
/// cannot tell twenty faults contributing equally from nineteen weighted at a
/// millionth of the twentieth, and at the edge of a field most nodes are the
/// second kind.
///
/// A `weights` of the wrong length, or holding a negative or a NaN, is refused
/// with a message rather than returning `None`: it is a fault in the call, and
/// silently pairing weights with the wrong faults would return a number that
/// looks like an answer.
#[pyfunction]
#[pyo3(signature = (faults, senses = None, weights = None, angle_step_degrees = 10.0, phi_step = 0.1))]
fn invert_stress<'py>(
    py: Python<'py>,
    faults: PyReadonlyArray2<'py, f64>,
    senses: Option<PyReadonlyArray1<'py, bool>>,
    weights: Option<PyReadonlyArray1<'py, f64>>,
    angle_step_degrees: f64,
    phi_step: f64,
) -> PyResult<Option<Bound<'py, PyDict>>> {

    let grid = checked_grid(angle_step_degrees, phi_step)?;
    let planes = faults_from_rows(&faults, senses.as_ref())?;
    let weights = checked_weights(weights.as_ref(), planes.len())?;

    // The whole point of the exercise is here, so the interpreter is let go
    // for the duration rather than held through a few million solutions.
    let Some(result) = py.detach(|| match weights {
        Some(ref weights) => invert_weighted(&planes, weights, grid),
        None => invert(&planes, grid),
    }) else {
        return Ok(None);
    };

    let runners_up = PyList::new(
        py,
        result
            .runners_up
            .iter()
            .map(|s| scored_tensor_dict(py, s))
            .collect::<PyResult<Vec<_>>>()?,
    )?;

    let dict = PyDict::new(py);
    dict.set_item("best", scored_tensor_dict(py, &result.best)?)?;
    dict.set_item("runners_up", runners_up)?;
    dict.set_item("candidates_tried", result.candidates_tried)?;
    dict.set_item("candidates_scored", result.candidates_scored)?;

    Ok(Some(dict))
}

/// How many candidate tensors `invert_stress` would visit on a given grid.
///
/// Before any are rejected as degenerate, so an upper bound rather than the
/// `candidates_tried` a run reports. Exists because the count is not the
/// product a caller would guess: an axis pointing up is the axis pointing
/// down, so plunge spans a quarter turn and the roll of S3 about S1 half a
/// one. The default grid holds 71 280 candidates, and halving both steps
/// multiplies that by about fourteen, which is worth knowing before starting
/// rather than after.
#[pyfunction]
#[pyo3(signature = (angle_step_degrees = 10.0, phi_step = 0.1))]
fn inversion_candidate_count(angle_step_degrees: f64, phi_step: f64) -> PyResult<usize> {
    Ok(checked_grid(angle_step_degrees, phi_step)?.candidate_count())
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(intersect_plane_grid, m)?)?;
    m.add_function(wrap_pyfunction!(intersect_mesh_grid, m)?)?;
    m.add_function(wrap_pyfunction!(best_fit_planes, m)?)?;
    m.add_function(wrap_pyfunction!(plane_normal, m)?)?;
    m.add_function(wrap_pyfunction!(rake_to_slickenline, m)?)?;
    m.add_function(wrap_pyfunction!(solve_stress, m)?)?;
    m.add_function(wrap_pyfunction!(stress_tensor, m)?)?;
    m.add_function(wrap_pyfunction!(score_stress, m)?)?;
    m.add_function(wrap_pyfunction!(invert_stress, m)?)?;
    m.add_function(wrap_pyfunction!(inversion_candidate_count, m)?)?;
    Ok(())
}
