# misah

*misah* is an experiment in creating a Rust-based geospatial library with Python bindings (thanks to the *pyo3* library).


Currently it is in alpha mode.

## Kernels

`src/kernels/` holds the numerical routines. They take plain slices and return
plain vectors, with no pyo3 and no I/O, so that `cargo test` exercises them
without a Python interpreter and the bindings can evolve separately.

### plane_dem

Intersection of an unbounded geological plane with a DEM, as the zero-level set
of the signed distance `f(x, y) = n . (X - P0)` sampled at the DEM nodes and
extracted by marching squares. Using the signed distance instead of a
`z = z(x, y)` plane expression keeps vertical planes an ordinary case, and the
cost is linear in the number of cells — geoSurfDEM's `IntersectDEM` solves the
more general surface-mesh problem by testing every DEM triangle against every
mesh triangle.

Two implementations are kept deliberately in step:

| | |
|---|---|
| `python/misah_ref/plane_dem.py` | numpy reference; also the fallback for environments where the compiled extension cannot be installed, QGIS plugins in particular |
| `src/kernels/plane_dem.rs` | Rust port |

`tests/test_plane_dem.py` validates the reference against the geoSurfDEM golden
dataset (plane 135/35 over the Malpi ASTER DEM, in the sibling `geoSurfDEM`
checkout) and pins the conventions both implementations must share. Run it with
`python3 tests/test_plane_dem.py`.

Two properties of that reference output are worth recording, since they set the
tolerances the tests use: about half of its 1228 points are duplicates (618
unique at centimetre precision), because shared triangle edges are emitted once
per triangle, and four points — CSV rows 644-647, in two pairs sharing an
elevation — lie up to 3.56 m off the plane. A plane fitted to the whole set
returns 134.9999/35.0005 and passes 8 mm from the nominal source point, so those
four are a defect of that run rather than a disagreement over conventions.

### Status

The Rust kernel has never been compiled: this machine has no Rust toolchain and
no access to crates.io. It also cannot be built as it stands, because the crate
is pinned to pyo3 0.13, which predates Python 3.11. Bringing the bindings up to a
current pyo3 and a maturin build is the next step, and needs a compiler in the
loop.
