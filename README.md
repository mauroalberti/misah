# misah

*misah* is an experiment in creating a Rust-based geospatial library with Python bindings (thanks to the *pyo3* library).


Currently it is in alpha mode.

## Plane-grid intersection

`lib/src/raster/intersection.rs` intersects a plane with a raster grid, as the
zero-level set of `Plane::signed_distance_to_point` sampled at the grid nodes and
extracted by marching squares. Building on the signed distance instead of a
`z = z(x, y)` plane expression keeps vertical planes an ordinary case, and the
cost is one pass over the cells — geoSurfDEM's `IntersectDEM` solves the more
general surface-mesh problem by testing every DEM triangle against every mesh
triangle.

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

`raster::io` does not compile and is not wired into the module tree, so the
example carries its own ESRI ASCII reader.

The Python surface is two functions, `intersect_plane_grid` and `plane_normal`.
Everything else in `lib` — the geometries, the orientations, the GeoProfiler
SQLite reader — is reachable from Rust only.

The setuptools-rust leftovers are gone: `features.rs`, the empty
`georeferenced.rs` and `orientations.rs`, and `setup.py` with `MANIFEST.in`,
superseded by maturin and broken besides — `setup.py` used an `install_requires`
it never defined.

Nothing yet builds wheels for anything but the host. A wheel built here is
tagged `manylinux_2_34`, since the extension picks up the glibc it is compiled
against, and will not install on Ubuntu 20.04, Debian 11 or RHEL 8. Wheels for
macOS, Windows and aarch64 need CI, and there is none: the Travis configuration
was dead and has been removed, and no GitLab equivalent has replaced it yet.

`docs/notebooks/misah.ipynb` calls an API that no longer exists
(`misah.orientations.orien3d.Axis`) and does not run.
