"""Pure-Python implementation of the kernels, used when the extension is absent.

`mispy/__init__.py` picks this up in place of the compiled `kernels` module, so
callers get the same names and the same results either way. That matters for
QGIS, whose plugins cannot rely on a binary wheel being installable and cannot
choose the interpreter they are loaded into.

It is also the oracle the Rust kernel is checked against: keeping a second,
independent implementation is what makes `tests/test_kernels.py` able to say the
two agree rather than merely that one of them runs.

Mirrors `lib/src/raster/intersection.rs`. The intersection is the zero-level set
of the signed distance

    f(x, y) = n . (X - P0),    X = (x, y, z_dem(x, y))

sampled at the grid nodes and extracted by marching squares. Working on the
signed distance rather than on a z = z(x, y) plane expression keeps vertical
planes an ordinary case.
"""

from __future__ import annotations

from typing import List, Optional, Sequence, Tuple

import numpy as np

# Corners are numbered 0 = top-left, 1 = top-right, 2 = bottom-right,
# 3 = bottom-left; edges are named after the corner pair they join.
_EDGE_CORNERS = {
    "top": (0, 1),
    "right": (1, 2),
    "bottom": (3, 2),
    "left": (0, 3),
}

# Chords per case index; bit i is set when corner i is positive. Complementary
# cases (i and 15-i) share their chords.
_CASE_EDGES = {
    0: (),
    1: (("left", "top"),),
    2: (("top", "right"),),
    3: (("left", "right"),),
    4: (("right", "bottom"),),
    5: (("left", "top"), ("right", "bottom")),      # saddle, resolved on centre
    6: (("top", "bottom"),),
    7: (("left", "bottom"),),
    8: (("left", "bottom"),),
    9: (("top", "bottom"),),
    10: (("left", "bottom"), ("top", "right")),     # saddle, resolved on centre
    11: (("right", "bottom"),),
    12: (("left", "right"),),
    13: (("top", "right"),),
    14: (("left", "top"),),
    15: (),
}

_SADDLE_ALTERNATIVE = {
    5: (("left", "bottom"), ("top", "right")),
    10: (("left", "top"), ("right", "bottom")),
}


def plane_normal(dip_dir_degr: float, dip_angle_degr: float) -> List[float]:
    """Upward-pointing unit normal of a geological plane, in (East, North, Up).

    Sanity anchors: a horizontal plane gives (0, 0, 1); a vertical one gives a
    horizontal normal pointing along the dip direction.
    """
    az = np.radians(dip_dir_degr)
    dip = np.radians(dip_angle_degr)
    return [
        float(np.sin(az) * np.sin(dip)),
        float(np.cos(az) * np.sin(dip)),
        float(np.cos(dip)),
    ]


def _node_coords(
    geotransform: Sequence[float], nrows: int, ncols: int
) -> Tuple[np.ndarray, np.ndarray]:
    """Ground coordinates of the cell centres, hence the half-cell offsets."""
    rows, cols = np.meshgrid(
        np.arange(nrows, dtype=float) + 0.5,
        np.arange(ncols, dtype=float) + 0.5,
        indexing="ij",
    )
    x = geotransform[0] + cols * geotransform[1] + rows * geotransform[2]
    y = geotransform[3] + cols * geotransform[4] + rows * geotransform[5]
    return x, y


def intersect_plane_grid(
    dem: np.ndarray,
    geotransform: Sequence[float],
    src_pt: Sequence[float],
    dip_dir_degr: float,
    dip_angle_degr: float,
    nodata: Optional[float] = None,
    epsg_code: int = 0,
) -> Tuple[np.ndarray, np.ndarray]:
    """Intersect an unbounded geological plane with a DEM.

    `geotransform` is the GDAL six-element affine transform. Returns
    `(points, segments)`: an (N, 3) array of intersection vertices and an (M, 2)
    array of indices into it, one row per marching-squares chord.
    """
    dem = np.asarray(dem, dtype=float)
    if dem.ndim != 2:
        raise ValueError("DEM must be two-dimensional")
    # The compiled kernel reads the grid as packed memory and refuses anything
    # else; refusing it here too keeps the fallback from quietly accepting what
    # the extension would reject.
    if not dem.flags["C_CONTIGUOUS"]:
        raise ValueError(
            "DEM must be C-contiguous; pass numpy.ascontiguousarray(dem)"
        )
    if len(geotransform) != 6:
        raise ValueError("geotransform must hold six elements")

    z = dem if nodata is None else np.where(dem == nodata, np.nan, dem)

    nrows, ncols = z.shape
    x, y = _node_coords(geotransform, nrows, ncols)

    n = plane_normal(dip_dir_degr, dip_angle_degr)
    x0, y0, z0 = src_pt
    f = n[0] * (x - x0) + n[1] * (y - y0) + n[2] * (z - z0)

    points: list = []
    segments: list = []
    # Vertices are shared between adjacent cells; key them by edge identity so a
    # crossing is emitted once and segments stay topologically connected.
    vertex_of_edge: dict = {}

    tiny = np.finfo(float).tiny

    def edge_vertex(r0, c0, r1, c1):
        key = (r0, c0, r1, c1)
        hit = vertex_of_edge.get(key)
        if hit is not None:
            return hit
        fa, fb = f[r0, c0], f[r1, c1]
        t = fa / (fa - fb)
        points.append(
            (
                x[r0, c0] + t * (x[r1, c1] - x[r0, c0]),
                y[r0, c0] + t * (y[r1, c1] - y[r0, c0]),
                z[r0, c0] + t * (z[r1, c1] - z[r0, c0]),
            )
        )
        vertex_of_edge[key] = len(points) - 1
        return len(points) - 1

    for r in range(nrows - 1):
        for c in range(ncols - 1):
            corners = ((r, c), (r, c + 1), (r + 1, c + 1), (r + 1, c))
            vals = [f[i, j] for i, j in corners]
            if not all(np.isfinite(v) for v in vals):
                continue
            # A node exactly on the plane would make the sign test ill-defined;
            # nudging it off by an ulp-scale epsilon keeps the topology consistent.
            vals = [v if v != 0.0 else tiny for v in vals]

            case = sum((1 << i) for i, v in enumerate(vals) if v > 0.0)
            chords = _CASE_EDGES[case]
            if not chords:
                continue

            if case in _SADDLE_ALTERNATIVE:
                # Saddle: the cell centre tells which of the two pairings is right.
                centre = sum(vals) / 4.0
                if (centre > 0.0) != (case == 5):
                    chords = _SADDLE_ALTERNATIVE[case]

            for e0, e1 in chords:
                ends = []
                for edge in (e0, e1):
                    a, b = _EDGE_CORNERS[edge]
                    (r0, c0), (r1, c1) = corners[a], corners[b]
                    ends.append(edge_vertex(r0, c0, r1, c1))
                segments.append(tuple(ends))

    return (
        np.array(points, dtype=float).reshape(-1, 3),
        np.array(segments, dtype=np.int64).reshape(-1, 2),
    )
