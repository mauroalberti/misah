"""Reference implementation of the plane-DEM intersection kernel.

This is the numpy oracle against which the Rust kernel in ``src/kernels/plane_dem.rs``
is validated, and the fallback used when the compiled extension is unavailable
(QGIS plugins cannot rely on binary wheels being installable).

The intersection is computed as the zero-level set of the signed distance

    f(x, y) = n . (X - P0),    X = (x, y, z_dem(x, y))

sampled at the DEM nodes, extracted by marching squares. Working on the signed
distance rather than on a z = z(x, y) plane expression keeps vertical planes
(where the plane has no single-valued z) as an ordinary case.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

# Marching-squares edge crossings per cell configuration. Corners are numbered
# 0=top-left, 1=top-right, 2=bottom-right, 3=bottom-left; edges are named after
# the corner pair they join. Case index has bit i set when corner i is positive.
_EDGE_CORNERS = {
    "top": (0, 1),
    "right": (1, 2),
    "bottom": (3, 2),
    "left": (0, 3),
}

_CASE_EDGES = {
    0: (),
    1: (("left", "top"),),
    2: (("top", "right"),),
    3: (("left", "right"),),
    4: (("right", "bottom"),),
    5: (("left", "top"), ("right", "bottom")),      # ambiguous: resolved on centre
    6: (("top", "bottom"),),
    7: (("left", "bottom"),),
    8: (("left", "bottom"),),
    9: (("top", "bottom"),),
    10: (("left", "bottom"), ("top", "right")),     # ambiguous: resolved on centre
    11: (("right", "bottom"),),
    12: (("left", "right"),),
    13: (("top", "right"),),
    14: (("left", "top"),),
    15: (),
}


@dataclass(frozen=True)
class GeoTransform:
    """GDAL-style affine transform, mapping (col, row) offsets to (x, y)."""

    x_origin: float
    pixel_width: float
    row_rotation: float
    y_origin: float
    col_rotation: float
    pixel_height: float

    @classmethod
    def from_esri_ascii_header(
        cls, xllcorner: float, yllcorner: float, cellsize: float, nrows: int
    ) -> "GeoTransform":
        return cls(
            x_origin=xllcorner,
            pixel_width=cellsize,
            row_rotation=0.0,
            y_origin=yllcorner + nrows * cellsize,
            col_rotation=0.0,
            pixel_height=-cellsize,
        )

    def node_coords(self, nrows: int, ncols: int) -> Tuple[np.ndarray, np.ndarray]:
        """Return (x, y) arrays of shape (nrows, ncols) for the DEM node centres.

        DEM values are point samples at cell centres, hence the +0.5 offsets.
        """
        rows, cols = np.meshgrid(
            np.arange(nrows, dtype=float) + 0.5,
            np.arange(ncols, dtype=float) + 0.5,
            indexing="ij",
        )
        x = self.x_origin + cols * self.pixel_width + rows * self.row_rotation
        y = self.y_origin + cols * self.col_rotation + rows * self.pixel_height
        return x, y


def plane_normal(dip_dir_degr: float, dip_angle_degr: float) -> np.ndarray:
    """Upward-pointing unit normal of a geological plane, in (East, North, Up).

    Sanity anchors: a horizontal plane (dip 0) gives (0, 0, 1); a vertical plane
    (dip 90) gives a horizontal normal pointing along the dip direction.
    """
    az = np.radians(dip_dir_degr)
    dip = np.radians(dip_angle_degr)
    return np.array(
        [np.sin(az) * np.sin(dip), np.cos(az) * np.sin(dip), np.cos(dip)],
        dtype=float,
    )


def signed_distance(
    dem: np.ndarray,
    geotransform: GeoTransform,
    src_pt: Tuple[float, float, float],
    dip_dir_degr: float,
    dip_angle_degr: float,
    nodata: Optional[float] = None,
) -> np.ndarray:
    """Signed distance from the plane to each DEM node, in ground units.

    NODATA nodes and non-finite elevations become NaN and are excluded downstream.
    """
    z = np.asarray(dem, dtype=float)
    if nodata is not None:
        z = np.where(np.isclose(z, nodata), np.nan, z)

    nrows, ncols = z.shape
    x, y = geotransform.node_coords(nrows, ncols)
    n = plane_normal(dip_dir_degr, dip_angle_degr)
    x0, y0, z0 = src_pt

    return n[0] * (x - x0) + n[1] * (y - y0) + n[2] * (z - z0)


def _interp_edge(fa, fb, xa, ya, za, xb, yb, zb):
    """Zero crossing along an edge, by linear interpolation on the signed distance."""
    t = fa / (fa - fb)
    return (xa + t * (xb - xa), ya + t * (yb - ya), za + t * (zb - za))


def intersect_plane_dem(
    dem: np.ndarray,
    geotransform: GeoTransform,
    src_pt: Tuple[float, float, float],
    dip_dir_degr: float,
    dip_angle_degr: float,
    nodata: Optional[float] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Intersect an unbounded geological plane with a DEM.

    Returns ``(points, segments)`` where ``points`` is an (N, 3) array of
    intersection vertices and ``segments`` an (M, 2) array of indices into
    ``points``, one row per marching-squares chord. Segments let callers rebuild
    polylines; the bare vertices feed attitude inversion directly.
    """
    f = signed_distance(
        dem, geotransform, src_pt, dip_dir_degr, dip_angle_degr, nodata
    )
    nrows, ncols = f.shape
    x, y = geotransform.node_coords(nrows, ncols)
    z = np.asarray(dem, dtype=float)
    if nodata is not None:
        z = np.where(np.isclose(z, nodata), np.nan, z)

    points: list = []
    segments: list = []
    # Vertices are shared between adjacent cells; key them by edge identity so a
    # crossing is emitted once and segments stay topologically connected.
    vertex_of_edge: dict = {}

    def edge_vertex(r0, c0, r1, c1):
        key = ((r0, c0), (r1, c1))
        hit = vertex_of_edge.get(key)
        if hit is not None:
            return hit
        pt = _interp_edge(
            f[r0, c0], f[r1, c1],
            x[r0, c0], y[r0, c0], z[r0, c0],
            x[r1, c1], y[r1, c1], z[r1, c1],
        )
        idx = len(points)
        points.append(pt)
        vertex_of_edge[key] = idx
        return idx

    for r in range(nrows - 1):
        for c in range(ncols - 1):
            # corner order: 0=(r,c) 1=(r,c+1) 2=(r+1,c+1) 3=(r+1,c)
            corners = ((r, c), (r, c + 1), (r + 1, c + 1), (r + 1, c))
            vals = [f[i, j] for i, j in corners]
            if not all(np.isfinite(v) for v in vals):
                continue
            # A node exactly on the plane would make the sign test ill-defined;
            # nudging it off by an ulp-scale epsilon keeps the topology consistent.
            vals = [v if v != 0.0 else np.finfo(float).tiny for v in vals]

            case = sum((1 << i) for i, v in enumerate(vals) if v > 0.0)
            edges = _CASE_EDGES[case]
            if not edges:
                continue

            if case in (5, 10):
                # Saddle: the cell centre tells which of the two pairings is right.
                centre = sum(vals) / 4.0
                if (centre > 0.0) != (case == 5):
                    edges = (
                        (("left", "bottom"), ("top", "right"))
                        if case == 5
                        else (("left", "top"), ("right", "bottom"))
                    )

            for e0, e1 in edges:
                idx = []
                for edge in (e0, e1):
                    a, b = _EDGE_CORNERS[edge]
                    (r0, c0), (r1, c1) = corners[a], corners[b]
                    idx.append(edge_vertex(r0, c0, r1, c1))
                segments.append(tuple(idx))

    pts = np.array(points, dtype=float).reshape(-1, 3)
    segs = np.array(segments, dtype=np.int64).reshape(-1, 2)
    return pts, segs


def read_esri_ascii(path: str) -> Tuple[np.ndarray, GeoTransform, float]:
    """Minimal ESRI ASCII grid reader, sufficient for the test fixtures."""
    header = {}
    with open(path) as src:
        for _ in range(6):
            key, value = src.readline().split()
            header[key.upper()] = float(value)
        data = np.loadtxt(src)

    nrows, ncols = int(header["NROWS"]), int(header["NCOLS"])
    data = data.reshape(nrows, ncols)
    gt = GeoTransform.from_esri_ascii_header(
        header["XLLCORNER"], header["YLLCORNER"], header["CELLSIZE"], nrows
    )
    return data, gt, header["NODATA_VALUE"]
