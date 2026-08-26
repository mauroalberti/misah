# misah

Python bindings for [misah](https://gitlab.com/mauroalberti/misah), a Rust crate
for structural geology processings.

**Alpha.** One kernel is exposed so far — the intersection of a geological plane
with a DEM — and the API may still change without notice.

## Install

```sh
pip install misah
```

Wheels are `abi3-py39`, so one of them serves every CPython from 3.9 on. Where no
wheel matches the platform, pip falls back to the source distribution and builds
it, which needs a Rust toolchain.

## Use

```python
import numpy as np
from misah.kernels import intersect_plane_grid

points, segments = intersect_plane_grid(
    np.ascontiguousarray(dem),   # C-contiguous, else ValueError
    geotransform,                # the six GDAL elements
    (x0, y0, z0),                # a point the plane passes through
    135.0,                       # dip direction, degrees
    35.0,                        # dip angle, degrees
    nodata,
)
```

`points` is an (N, 3) array of intersection vertices, `segments` an (M, 2) array
of indices into it, one row per marching-squares chord.

`misah.kernels.plane_normal(dip_dir, dip_angle)` returns the upward-pointing unit
normal in (East, North, Up).

## Compiled or not

`misah.kernels` is the compiled extension where it can be built, and the pure
Python `misah/_reference.py` where it cannot. Both expose the same names and
return the same arrays; `misah.is_compiled` says which one answered.

QGIS is the case this exists for: a plugin neither picks the interpreter it is
loaded into nor can count on a binary wheel installing, and vendoring the source
tree works for Python and not for a compiled extension. On the Malpi DEM the
fallback takes about 300 ms against the 1 ms of the compiled kernel, the gap
being the per-cell Python loop rather than the algorithm.

Note the fallback does not rescue a `pip install`: with no matching wheel, pip
builds the sdist and that needs Rust. It covers the vendored-source case.

Keeping a second implementation also buys the tests an oracle: `tests/test_reference.py`
compares the two elementwise over four attitudes, and so can say they agree
rather than merely that one of them runs.

## Accuracy

The kernel was checked against the geoSurfDEM golden dataset — plane 135/35 over
the Malpi ASTER DEM — giving 404 vertices, all on the plane to 1e-9 and
coincident with the C++ trace to well within a cell.

## License

GPL-3.0-or-later. See [LICENSE](https://gitlab.com/mauroalberti/misah/-/blob/main/LICENSE).
