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
checkout) and pins the conventions both implementations must share.
`tests/test_rust_equivalence.py` then runs the Rust kernel over the same DEM
through `examples/plane_dem_asc` and compares the two outputs elementwise. Both
walk cells in row-major order and key vertices by edge, so their point sequences
correspond one to one; on the Malpi DEM they agree to better than 1e-6 m on an
inclined and on a vertical plane alike.

```sh
python3 tests/test_plane_dem.py         # conventions, against the golden dataset
python3 tests/test_rust_equivalence.py  # Rust and the extension, against numpy
cargo test --no-default-features        # kernel unit tests
```

On the Malpi DEM (200 x 247 cells) the compiled kernel runs in about 1 ms
against 300 ms for the reference — the reference walks cells in a Python loop,
so the gap is the loop, not the algorithm.

Two properties of that reference output are worth recording, since they set the
tolerances the tests use: about half of its 1228 points are duplicates (618
unique at centimetre precision), because shared triangle edges are emitted once
per triangle, and four points — CSV rows 644-647, in two pairs sharing an
elevation — lie up to 3.56 m off the plane. A plane fitted to the whole set
returns 134.9999/35.0005 and passes 8 mm from the nominal source point, so those
four are a defect of that run rather than a disagreement over conventions.

## Building

```sh
maturin develop --release   # build and install into the active environment
maturin build --release     # produce a wheel under target/wheels
```

The crate builds against pyo3 0.29 with `abi3-py39`, so one wheel serves every
Python from 3.9 on. That matters for the QGIS side of the suite: a plugin cannot
choose the interpreter it is loaded into, and cannot count on a compiler being
present to build a version-specific extension.

`extension-module` deliberately leaves libpython unlinked, which a test binary
cannot do, so it is a default feature that the Rust tests turn off:
`cargo test --no-default-features`. It also gates `kernels::py`, and with it the
numpy dependency, keeping the kernels themselves buildable on their own.

### Status

Alpha. The nested modules are attributes of the extension rather than importable
submodules, so it is `misah.geometry.geom2d.Point2D`, not
`from misah.geometry.geom2d import Point2D`. `GeoTransform` has no constructor
exposed and cannot be instantiated from Python.
