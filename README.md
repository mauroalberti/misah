# misah

*misah* is an experiment in creating a Rust-based geospatial library with Python bindings (thanks to the *pyo3* library).


Currently it is in alpha mode.

## Plane-grid intersection

`lib/src/raster/intersection.rs` intersects a plane with a raster grid, as the
zero-level set of `Plane::signed_distance_to_point` sampled at the grid nodes and
extracted by marching squares. Building on the signed distance instead of a
`z = z(x, y)` plane expression keeps vertical planes an ordinary case, and the
cost is one pass over the cells. The general case, a surface of arbitrary shape,
is `raster::mesh_intersection` below.

`GeologicalPlane::to_plane` bridges an attitude to the geometric plane, so a
caller works in dip direction and dip angle rather than in normals.

```sh
cargo test -p misah
cargo run --release -p misah --example plane_dem_asc -- \
    <dem.asc> <x0> <y0> <z0> <dip_dir> <dip_angle>
```

The kernel was checked against the geoSurfDEM golden dataset (plane 135/35 over
the Malpi ASTER DEM): 404 vertices, all on the plane to 1e-9, coincident with the
C++ trace to well within a cell. Two properties of that reference output are
worth recording, since they bound how closely anything can agree with it: about
half of its 1228 points are duplicates (618 unique at centimetre precision),
because shared triangle edges are emitted once per triangle, and four points —
CSV rows 644-647, in two pairs sharing an elevation — lie up to 3.56 m off the
plane. A plane fitted to the whole set returns 134.9999/35.0005 and passes 8 mm
from the nominal source point, so those four are a defect of that run rather than
a disagreement over conventions.

## Mesh-grid intersection

`lib/src/raster/mesh_intersection.rs` cuts the grid with a triangulated surface
rather than a single plane — the case a folded or faulted geological surface
needs. It is a port of geoSurfDEM's `IntersectDEM`, the one part of that suite
with no equivalent in geogst. Each point comes out with the attitude of the mesh
triangle that produced it, so the result is a set of located attitudes, ready
for `BestFitGeoplanes`-style inversion, and not merely a trace.

```sh
cargo run --release -p misah --example mesh_dem_vtk -- <dem.asc> <surface.vtk>
```

`geometry::triangle::Triangle3D` and `geometry::mesh::TriangleMesh` carry the
surface; the mesh is indexed (a vertex pool plus index triplets) because that is
how a VTK `POLYDATA` file gives it. The DEM is read by `raster::io`; reading VTK
is left to the example, being no business of a raster module.

Three things depart from the C++ deliberately:

- **Where the points come from.** The original intersected the two triangle
  planes into a line, then intersected that line with each DEM side. Here a DEM
  side is cut against the mesh triangle's plane directly, interpolating between
  the signed distances of its two endpoints — the construction the plane kernel
  already uses. Same points, without the intermediate line, whose construction
  needed an arbitrary 100-unit displacement, an unguarded division, and could
  return an uninitialised point when all three of its determinants fell below
  tolerance.
- **What gets compared.** A DEM *is* a uniform spatial index, so each mesh
  triangle is mapped through the inverse geotransform (`GeoTransform::node_rc`,
  new) onto the cells its footprint covers, widened by one cell so that a
  surface passing exactly through a column of nodes still finds its crossings.
  No index is built, because the grid already is one.
- **Shared sides.** Crossings are keyed on the DEM side that produced them, so
  the side two adjacent DEM triangles have in common is cut once. This is the
  source of the duplication in the reference output noted above; the count of
  suppressed repeats is reported in the run's statistics rather than hidden.

Against the same golden dataset (surface 135/35 over the Malpi ASTER DEM), the
port returns 614 points against the reference's 618 distinct ones, every one
reported at 135.00/35.00, in 29 ms. Two independent checks, both at the 1 cm the
CSV is written to:

| | on the plane, max residual | off the DEM surface, median / max |
| --- | --- | --- |
| misah | 0.008 m | 0.003 m / 0.010 m |
| C++ reference | 3.55 m | 0.038 m / 0.213 m |

Both figures for misah are the CSV rounding itself, so the points sit on the
plane and on the topography as exactly as the output can express; the reference
misses the DEM surface by more than a centimetre at 511 of its 617 testable
points, which is where the ~0.18 m median offset between the two point sets
comes from. The spatial restriction takes the comparisons from the 49 million
pairs the C++ makes after its own volume prefilter (505 mesh triangles × 97,908
DEM triangles) down to 23,606.

Note that `test_data/src_planes/malpi_045_90.vtk` in geoSurfDEM is a byte-copy
of `malpi_135_35.vtk`, so that suite has no vertical-surface case despite the
name; the vertical and node-aligned configurations are covered by unit tests
instead.

## Forward stress on a fault plane

`lib/src/structural/stress.rs` is the direct Wallace-Bott problem: given a
reduced stress tensor (Angelier, 1984 — S1/S3 orientation and the shape ratio
Phi), predict the slip a fault plane would show under it. It is the inverse of
what a fault-slip inversion solves for, and it is a port of `ForwardStress.f95`,
a Fortran tool of Alberti (2010) that once lived in this repository's own
history under `GeoFaults`, before being dropped for a Rust port that was never
carried through.

This one arrived by way of geogst first (`geogst.core.geology.stress`), where
`Plane.rake_to_direct` turned out to already implement the Fortran's own
rake-to-slickenline formula exactly — both are Aki & Richards (1980) — and
where the whole algorithm was checked against the Fortran source itself:
compiled unmodified with gfortran, driven by a small program calling
`stresssolution_calc` directly, over five geological cases. Three matched
outright. Two did not, and the reason was a real bug in the 2010 source:
`module vector_processing` carries no `implicit none`, so `vector1_magn` inside
`vector_normalization` — reached through `vector_projection`, used for the
normal-stress vector every solution goes through — is an undeclared name,
silently typed as single-precision `REAL` by Fortran's implicit-typing rule
rather than the double precision the rest of the tool uses. The truncation is
usually too small to matter, but the predicted rake reaches it through an
`acos` close to the edge of its domain in near-strike-slip cases, where an
error below float32 precision in the cosine becomes one four orders of
magnitude larger in the angle — 1.4e-5 degrees on one case, arrived at through
a chain of double-precision arithmetic that looks correct at a glance.
Recompiling with `-fimplicit-none`, declaring `vector1_magn`'s type, moves
those two cases onto what the geogst port already computed, agreeing to 12+
significant digits; an independent 50-digit `mpmath` recomputation of the same
five cases confirms the corrected numbers, not the original binary's, are the
mathematically correct ones.

This Rust port is checked against the same five cases and the same corrected
values — `structural::stress::tests::against_the_fortran_original` — so the
same bug does not need rediscovering here.

`ReducedStressTensor::new` takes S1 and S3 as `GeologicalAxis`, sub-orthogonal,
plus Phi; S2 completes a right-handed cyclic triad (`S1 x S2 = S3`,
`S2 x S3 = S1`), matching the sign the Fortran original builds. `sigma1`/`sigma3`
default to 1/0 through `ReducedStressTensor::normalized` — the usual
normalization when only the tensor's shape is known, as from an inversion —
which leaves the predicted rake correct while making slip tendency and
deformation index meaningless; supply the true magnitudes when they are known
and those are wanted. `solve` resolves the tensor onto a `GeologicalPlane` and
returns a `StressSolution`: traction, normal and shear stress always, and,
where the shear is not numerical noise, the predicted rake, slickenline, slip
tendency and deformation index. `angular_misfit` compares the prediction
against an observed slickenline — what an inversion would minimize, run
forward.

Two small primitives came out of this that the rest of the structural code
will reuse: `structural::geol_axis::GeologicalAxis::as_versor`/`from_versor`
(trend/plunge in `(East, North, Up)` and back), and
`GeologicalPlane::rhr_strike`/`from_rhr_strike`/`rake_to_versor`, the last
being the Aki & Richards formula itself, exposed rather than kept private to
`solve`, since it is the piece any future fault-slip *inversion* would need to
call once per candidate tensor per fault — which is also the reason this went
into `misah` and not only geogst: the direct problem is O(1) per fault, with no
performance case of its own, but it is the exact inner loop a grid search over
candidate tensors would evaluate, at whatever count of faults and grid points
made that search worth compiling.

`solve_stress` and `rake_to_slickenline`, in `pylib`, are the Python surface of
this — one call each, taking and returning plain numbers rather than exposing
`ReducedStressTensor`/`GeologicalAxis`/`GeologicalPlane` themselves, matching
how `intersect_plane_grid` already builds and consumes a `GeologicalPlane`
without handing the type to the caller. Both are checked from Python against
the same five Fortran-verified cases.

## Fault-slip inversion

`lib/src/structural/inversion.rs` runs the forward model backwards: given
faults with their observed slip, it searches for the reduced stress tensor that
best explains them. That search is what made the forward problem worth
compiling, and what the sentence above was promising — one forward solution per
candidate tensor per fault, which for the default grid of some 21 000
candidates over a hundred faults is a couple of million of them.

`FaultPlane` is the input, and is what connects a group of types that until now
led nowhere: `Direction3D` was reached only by `Slickenline`, `Slickenline` only
by `FaultPlane`, and `FaultPlane` by nothing at all. It now validates the
invariant its absence let through — that a slickenline must lie *in* the plane
it is recorded on, within a degree, since a lineation off the plane is not a
fault with a small error but two measurements that do not belong together.

Two details of the domain that the geometry alone would get wrong:

- **An unread slip sense is a line, not a vector.** Where the sense of movement
  was not determined, a prediction pointing the other way along the same
  lineation fits perfectly, so `Slickenline::angle_to` takes the misfit modulo
  180 degrees. Scoring such a fault at 180 instead of 0 would penalise an
  inversion for exactly the faults nobody could read the sense of.
- **A fault the model cannot speak about is not a fault it fits badly.** A plane
  lying on a principal stress axis carries no shear, so no slip direction is
  predicted; `misfit_on` returns `None` there rather than zero, and
  `mean_misfit` reports how many faults its average actually came from, so a
  candidate that silences half the dataset cannot win on the strength of the
  remainder.

The search is exhaustive rather than a descent, because the misfit surface is
not convex — a fault set carrying two superposed tectonic phases has two minima
by construction, and a descent would report whichever it fell into without ever
saying the other was there. `runners_up` is returned alongside `best` for the
same reason: a minimum that stands alone reads differently from one on a
plateau.

The test that matters generates faults from a known tensor through the forward
model and asks the inversion to recover it, which it does to within the grid
step for a vertical-S1 and for a strike-slip setting alike.

Not yet exposed to Python — the one part of the structural side that is not.
The search is where being compiled matters most, so it is the obvious next
binding; what it needs first is a decision about its own shape, since a caller
handing over a hundred faults wants to pass arrays rather than build a hundred
objects across the boundary.

## Python bindings

`pylib` builds the `misah` Python distribution with maturin, against pyo3 0.29
and `abi3-py39`, so one wheel serves every Python from 3.9 on. That matters for
the QGIS side of the suite: a plugin cannot choose the interpreter it is loaded
into, nor count on a compiler being available to build a version-specific
extension.

Three names, deliberately distinct. The distribution installed from PyPI is
`misah`; the crate under `pylib` is `misah-py`, because the core crate in `lib`
already holds `misah` and a workspace admits one package per name; the compiled
extension is `misah._misah`, private because callers reach its contents through
the aliases `misah/__init__.py` sets up.

```sh
cd pylib
maturin develop --release   # build and install into the active environment
maturin build --release     # produce a wheel under target/wheels
python3 tests/test_kernels.py
```

```python
import numpy as np
from misah.kernels import intersect_plane_grid

points, segments = intersect_plane_grid(
    np.ascontiguousarray(dem),   # C-contiguous, else ValueError
    geotransform,                # the six GDAL elements
    (x0, y0, z0),
    135.0,                       # dip direction
    35.0,                        # dip angle
    nodata,
)
```

```python
from misah.kernels import intersect_mesh_grid

# vertices (V, 3) and faces (F, 3) as a VTK POLYDATA file gives them
points, attitudes, mesh_triangles, stats = intersect_mesh_grid(
    np.ascontiguousarray(dem),
    geotransform,
    vertices,
    faces,
    nodata,
)
attitudes[0]        # (dip direction, dip angle) of the triangle that cut here
stats["duplicate_crossings"]
```

```python
from misah.kernels import solve_stress

solution = solve_stress(
    0.0, 90.0,    # S1: trend, plunge
    90.0, 0.0,    # S3: trend, plunge
    0.5,          # Phi
    0.0, 60.0,    # fault: RHR strike, dip angle
    sigma1=30.0, sigma3=10.0,   # optional; default 1/0 leaves the rake correct
)
solution["theoretical_rake"]           # -90.0: pure normal dip-slip
solution["theoretical_slickenline"]    # (90.0, 60.0), down the dip vector
```

The extension is built as `misah._misah`, so its submodules register themselves
under that name; `misah/__init__.py` aliases them, which is what makes
`import misah.kernels` work rather than only `misah._misah.kernels`.

There is no pure-Python fallback. There was one, `misah/_reference.py`,
mirroring every exposed kernel so that `misah.kernels` answered either way and
`misah.is_compiled` said which; it is gone, deliberately.

It cost a second independent implementation of every function exposed, and that
price is what kept the mesh-grid intersection unexposed for as long as it was:
marching triangles over a spatially indexed DEM is not something worth writing
twice, so the rule that was meant to protect the bindings was instead deciding
which half of the library Python could reach. Dropping it is what let
`intersect_mesh_grid` above be exposed at all.

What it bought was real and is worth naming: an oracle. Two implementations
written from the same formulas, compared elementwise, can say they *agree*
rather than that one of them runs — and it earned its keep, catching
`solve_stress` returning its vectors as lists from one side and tuples from the
other. In its place the bindings are checked against the datasets the kernels
were ported from: the geoSurfDEM golden trace for the raster ones, the five
Fortran-verified cases for the stress ones, and expected values rather than
mutual agreement throughout.

What it also bought, and this is the loss: a QGIS plugin that cannot install a
binary wheel could still use misah, slowly. It now cannot use it at all. That
is the same problem as the wheels, and belongs there — building for the
platforms QGIS runs on — rather than in a second implementation of everything
kept in reserve.

The test files split by what they need, not by what they cover.
`test_stress.py` and `test_mesh.py` are numbers in, numbers out, so both run in
CI; `test_kernels.py` reads the Malpi crop from the geoSurfDEM repository and
runs by hand. The DEM committed under `example_data` is a wider crop of the same
ASTER tile — 213x260 against 200x247 — so it cannot stand in without rewriting
the vertex counts those tests assert.

Run the tests from anywhere except `pylib` itself, whose source tree shadows the
installed package.

### Status

`raster::io` reads ESRI ASCII grids, which is what both examples now use; it is
the only raster format handled, anything wider meaning GDAL. It returns the
nodata value as an `Option` rather than defaulting to -9999 on a file's behalf,
since at sea that is a depth and not a hole.

`structural::focal_mechanism` is still an empty struct, and
`gp_projected_focal_mechanisms` is read from a GeoProfiler export into records
that nothing yet computes with — the next thing in this direction, and the one
`structural::stress` would extend to, focal mechanisms being fault-slip data
that arrives already reduced to a plane and a rake.

The Python surface is five functions: `intersect_plane_grid`,
`intersect_mesh_grid` and `plane_normal` for the raster side, `solve_stress`
and `rake_to_slickenline` for the structural side. Each builds and consumes the
types it needs within the one call rather than handing them across the
boundary, so the surface stays plain arrays and numbers in, arrays and dicts
out.

What remains Rust-only is `structural::inversion`, and the GeoProfiler SQLite
reader. The reader is a different case from the others: Python has `sqlite3` in
its standard library and can read a GeoProfiler export directly, so the Rust
reader exists for Rust callers rather than standing between Python and the
data.

`geoprofile::sqlite` reads a qgSurf GeoProfiler export — all thirteen tables of
it — and `lib/tests` runs that against two: a schema-v1 export of the Timpa San
Lorenzo section, and a small v5 one written by qgSurf's own exporter through its
public insert functions, regenerable with `lib/tests/data/make_export_synthetic_v5.py`.
The two projects agree on a file format that nothing else spans, so the fixture
is that agreement rather than a transcription of it. The four added after schema
v1 (`gp_projected_focal_mechanisms`, `gp_profile_vertices`,
`gp_graphical_params`, `gp_source_categories`) are read when present and come
back empty when not, so one reader serves every export written so far.

A line intersection needs both shapes. A line generally crosses a section at a
point, but a segment lying along the section trace crosses it over a stretch;
until qgSurf's schema v5 the table stored a single distance and such a stretch
was written as two ordinary rows, indistinguishable afterwards from two separate
crossings. v5 records `s_from` and `s_to`, as `gp_intersected_polygons` always
has. Both are read: the older column becomes the degenerate span it stands for.

The setuptools-rust leftovers are gone: `features.rs`, the empty
`georeferenced.rs` and `orientations.rs`, and `setup.py` with `MANIFEST.in`,
superseded by maturin and broken besides — `setup.py` used an `install_requires`
it never defined.

`.gitlab-ci.yml` runs the tests on Linux at every push, builds the workspace and
the examples, gates on clippy, and installs the built wheel to run
`pylib/tests/test_stress.py` against it — the Python bindings imported and
called, not merely compiled. It stops short of building wheels for anything but
that one Linux target.

Nothing yet builds wheels for anything but the host. A wheel built here is
tagged `manylinux_2_34`, since the extension picks up the glibc it is compiled
against, and will not install on Ubuntu 20.04, Debian 11 or RHEL 8. Wheels for
macOS, Windows and aarch64 need runners for those platforms, which is a question
about the GitLab plan rather than about the code, and until it is answered
nothing in `lib` can reach a QGIS installation that is not this one.

`docs/notebooks/misah.ipynb` is gone rather than repaired. Every path in it was
dead, not merely the one previously named here: `misah.geometry.geom2d`,
`misah.geometry.geom3d`, `misah.orientations.orien3d`,
`misah.georeferenced.georef2d` — the last of which was a file deleted two
commits before this line was written — and it compared results against `pygsf`,
which geogst itself superseded. Nothing in it touched any of the four functions
the package exposes today. A notebook covering those would be worth having, and
would be a new one.

`lib/examples/test.rs` is gone with it: a nalgebra hello-world computing the
distance between two points, referring to nothing in this crate, and the sole
reason `nalgebra` appeared in `lib/Cargo.toml` — 97 lines of lockfile and a
whole dependency tree for an example that demonstrated nothing about misah.

Both crates are now on edition 2021. `lib` was on 2018, which is what let
`assert!(cond, "…{err}")` pass as a plain string rather than a format string,
and so print the placeholder instead of the error at the moment a test failed;
ten such assertions had to be corrected by hand once clippy pointed at them.
The edition closes the class rather than the instances. `cargo fix --edition`
needed no source changes to make the move.
