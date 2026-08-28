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

The extension is built as `misah._misah`, so its submodules register themselves
under that name; `misah/__init__.py` aliases them, which is what makes
`import misah.kernels` work rather than only `misah._misah.kernels`.

Where the extension cannot be built or installed, `misah.kernels` is instead the
pure-Python `misah/_reference.py`, under the same names and returning the same
arrays; `misah.is_compiled` says which one answered. QGIS is the case this exists
for: a plugin neither picks the interpreter it is loaded into nor can count on a
binary wheel installing, and vendoring — which is how qgSurf carries geogst —
works for Python and not for a compiled extension. On the Malpi DEM the fallback
takes about 300 ms against the 1 ms of the compiled kernel, the gap being the
per-cell Python loop rather than the algorithm.

Keeping a second implementation also buys the tests an oracle:
`tests/test_reference.py` compares the two elementwise over four attitudes, and
so can say they agree rather than merely that one of them runs.

Run the tests from anywhere except `pylib` itself, whose source tree shadows the
installed package.

### Status

`raster::io` reads ESRI ASCII grids, which is what both examples now use; it is
the only raster format handled, anything wider meaning GDAL. It returns the
nodata value as an `Option` rather than defaulting to -9999 on a file's behalf,
since at sea that is a depth and not a hole.

The Python surface is two functions, `intersect_plane_grid` and `plane_normal`.
Everything else in `lib` — the geometries, the orientations, the mesh-grid
intersection, the GeoProfiler SQLite reader — is reachable from Rust only.

The setuptools-rust leftovers are gone: `features.rs`, the empty
`georeferenced.rs` and `orientations.rs`, and `setup.py` with `MANIFEST.in`,
superseded by maturin and broken besides — `setup.py` used an `install_requires`
it never defined.

`.gitlab-ci.yml` runs the tests on Linux at every push, and builds the workspace
and the examples. It stops there: it does not build wheels.

Nothing yet builds wheels for anything but the host. A wheel built here is
tagged `manylinux_2_34`, since the extension picks up the glibc it is compiled
against, and will not install on Ubuntu 20.04, Debian 11 or RHEL 8. Wheels for
macOS, Windows and aarch64 need runners for those platforms, which is a question
about the GitLab plan rather than about the code, and until it is answered
nothing in `lib` can reach a QGIS installation that is not this one.

`docs/notebooks/misah.ipynb` calls an API that no longer exists
(`misah.orientations.orien3d.Axis`) and does not run.
