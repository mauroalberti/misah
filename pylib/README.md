# misah

Python bindings for [misah](https://gitlab.com/mauroalberti/misah), a Rust crate
for structural geology processings: cutting geological surfaces against a DEM,
reading attitudes back out of the result, and solving the Wallace-Bott problem
forwards and backwards.

**Alpha.** Ten functions are exposed and the API may still change without
notice.

## Install

```sh
pip install misah
```

That installs the alpha, and no `--pre` is needed for it: pip passes over a
pre-release only when there is a stable version to prefer instead, and misah has
none yet. From the first stable release on, the alphas stop arriving that way,
and `pip install --pre misah` is what reaches them.

Wheels are `abi3-py39`, so one of them serves every CPython from 3.9 on. They
are built for **Linux only** so far — `manylinux_2_17`, x86-64 and aarch64,
which reaches back to glibc 2.14 and 2.17 respectively. Everywhere else pip
falls back to the source distribution and builds it, which needs a Rust
toolchain.

That is a real limitation rather than an oversight, and it is worth stating
where it bites: QGIS, which is the reason these bindings exist, mostly runs on
Windows and macOS. Windows needs LLVM and the MSVC target to cross-compile —
`rusqlite` is bundled, so the build wants a C compiler for Windows and not only
a Rust one — and macOS needs the Apple SDK, a licensing question before a
technical one. Until those are solved, a QGIS plugin can depend on misah on
Linux and nowhere else. The repository README records both attempts in detail.

## Use

Every function takes and returns plain arrays, numbers and dicts. The Rust types
— `GeologicalPlane`, `ReducedStressTensor`, `FaultPlane` — are built and consumed
inside the one call rather than handed across the boundary.

### Surfaces against a DEM

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

For a folded or faulted surface rather than a single plane, `intersect_mesh_grid`
takes the triangulation and returns the attitude of the triangle that produced
each point, which makes the result a set of located measurements and not a bare
trace:

```python
from misah.kernels import intersect_mesh_grid

points, attitudes, mesh_triangles, stats = intersect_mesh_grid(
    np.ascontiguousarray(dem),
    geotransform,
    vertices,                    # (V, 3), as a VTK POLYDATA file gives them
    faces,                       # (F, 3) indices into the vertex pool
    nodata,
)
attitudes[0]                     # (dip direction, dip angle)
stats["duplicate_crossings"]
```

`misah.kernels.plane_normal(dip_dir, dip_angle)` returns the upward-pointing unit
normal in (East, North, Up).

### Attitudes back out of points

```python
from misah.kernels import best_fit_planes

field, stats = best_fit_planes(points, cell_size=50.0)

field["attitudes"]       # (M, 2) dip direction and dip angle, one per cell
field["cell_centres"]    # (M, 2) where to post them
field["collinearity"]    # (M,) how well each is determined; large is bad
stats["cells_collinear"] # cells set aside: on a smooth slope, most of them
```

`None` where there is nothing to fit. Any (N, 3) array will do — a contact
digitised from a map, or readings along an outcrop, as readily as the output of
the kernel above.

The cells set aside are the point of it. A plane is determined by a cloud of
points only if the cloud spreads in two directions, and a contact crossing a
smooth hillside spreads in one, where every plane through that line fits equally
well and the one returned is whichever direction rounding favoured. Those cells
are counted in `stats` rather than published with a confident attitude.

### Stress on a fault plane

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

`rake_to_slickenline` is the Aki & Richards rake-to-lineation formula on its own,
for callers that want it without a tensor.

### And the same problem backwards

```python
from misah.kernels import invert_stress, inversion_candidate_count

# faults: (N, 4) -- RHR strike and dip angle, then the trend and plunge of the
# slickenline on that plane. Every lineation must lie in its own plane to within
# a degree, or the offending row is named and refused.
result = invert_stress(
    faults,
    senses=None,             # (N,) bool: False where the sense was not read
    weights=None,            # (N,) float >= 0: how much each fault counts for
    angle_step_degrees=10.0,
    phi_step=0.1,
)

result["best"]["s1"]                   # (trend, plunge)
result["best"]["phi"]
result["best"]["mean_misfit_degrees"]
result["best"]["faults_scored"]        # how many faults that mean came from
result["best"]["effective_sample_size"]  # the count, once weighting is allowed for
result["runners_up"]                   # the next five, worst last
```

`senses` is the one argument worth reading twice. Omitting it claims every sense
of movement was determined, which is right for data taken as given and wrong for
a dataset where the sense went unrecorded: those faults would be scored at 180
degrees for fitting perfectly the other way round. Pass `False` for them and the
misfit is taken modulo 180.

`weights` is what makes a stress *field* possible. Each fault counts for as much
as its weight says, so pass the whole dataset once and recompute only the weights
at each node of a grid, from a kernel of the distance between that node and where
each fault was measured. Faults weighted at zero are dropped before any forward
solution is computed, so with a kernel of finite reach a node costs what its own
neighbourhood costs rather than what the dataset costs. Only the relative sizes
are read; they need not sum to anything.

It is also the only handle for separating two superposed tectonic phases, which
nothing else in this search can be told to prefer between.

Weighting changes which number says whether two misfits are comparable.
`faults_scored` counts twenty faults whether they contributed equally or whether
nineteen were weighted at a millionth of the twentieth; `effective_sample_size`
is Kish's `(sum w)^2 / sum w^2`, which says twenty in the first case and barely
one in the second. Unweighted the two are equal, exactly. At the edge of a field
most nodes are the second kind, so read the effective size there.

The search is exhaustive rather than a descent, because a fault set carrying two
superposed tectonic phases has two minima by construction and a descent would
report whichever it fell into. `inversion_candidate_count` sizes a run before
starting it — the default grid is 71 280 candidates, and halving both steps
multiplies the work by about fourteen. The GIL is released for the search, so a
hundred faults on the default grid holds it for none of the 1.5 s they take.
`None` comes back when no candidate could be scored on any fault, and also when
every weight is zero — for a field, a node with no data within reach, which is a
result rather than a failure. A `weights` of the wrong length, or holding a
negative or a NaN, raises instead: that is a fault in the call, not in the data.

```python
from misah.kernels import score_stress, stress_tensor

# The same misfit, for a tensor being tested rather than looked for.
score = score_stress(*best["s1"], *best["s3"], best["phi"], faults)
score["misfits"]            # (N,) in the order given; NaN where none predicted

# And the tensor as a matrix, in (East, North, Up).
matrix = stress_tensor(*best["s1"], *best["s3"], best["phi"],
                       sigma1=30.0, sigma3=10.0)
```

## There is no pure-Python fallback

There was one, `misah/_reference.py`, mirroring every exposed kernel so that
`misah.kernels` answered either way. It was dropped deliberately: it cost a
second independent implementation of every function exposed, and that price was
what kept the mesh-grid intersection — marching triangles over a spatially
indexed DEM — unexposed for as long as it was. What it bought was an oracle, two
implementations that could be said to *agree* rather than one that ran. In its
place the bindings are checked against the datasets the kernels were ported
from.

The loss is the QGIS case above: a plugin that cannot install a binary wheel now
cannot use misah at all, rather than falling back to something slower. That is
the same problem as the wheels and belongs there, in building for the platforms
QGIS runs on.

## Accuracy

The kernels are ports, each checked against what it was ported from.

The raster ones against geoSurfDEM's golden dataset, plane 135/35 over the Malpi
ASTER DEM: 404 vertices for the plane kernel, all on the plane to 1e-9; 614
points for the mesh kernel against the reference's 618 distinct ones, every one
reported at 135.00/35.00. On geoSurfDEM's Timpa San Lorenzo dataset the best-fit
kernel returns 288 cells whose worst error against the plane that generated the
points is 0.089°, where the reference publishes 311 and its worst is 38.4°.

The stress ones against `ForwardStress.f95`, Alberti (2010), over five
geological cases — and against the *corrected* Fortran: two of the five moved
when a missing `implicit none` in the original was fixed, and an independent
50-digit recomputation confirms the corrected numbers.

The repository README carries the full comparisons, including what they cost and
where the two disagree.

## License

GPL-3.0-or-later. See
[LICENSE](https://gitlab.com/mauroalberti/misah/-/blob/master/LICENSE).
