"""Pure-Python implementation of the kernels, used when the extension is absent.

`misah/__init__.py` picks this up in place of the compiled `kernels` module, so
callers get the same names and the same results either way. That matters for
QGIS, whose plugins cannot rely on a binary wheel being installable and cannot
choose the interpreter they are loaded into.

It is also the oracle the Rust kernel is checked against: keeping a second,
independent implementation is what makes `tests/test_kernels.py` able to say the
two agree rather than merely that one of them runs.

`intersect_plane_grid`/`plane_normal` mirror `lib/src/raster/intersection.rs`.
The intersection is the zero-level set of the signed distance

    f(x, y) = n . (X - P0),    X = (x, y, z_dem(x, y))

sampled at the grid nodes and extracted by marching squares. Working on the
signed distance rather than on a z = z(x, y) plane expression keeps vertical
planes an ordinary case.

`rake_to_slickenline`/`solve_stress` mirror `lib/src/structural/stress.rs` and
the two structural primitives it drew out of `geol_axis.rs`/`geol_plane.rs` --
trend/plunge <-> (East, North, Up), and the Aki & Richards (1980) rake formula.
`solve_stress` is the direct (forward) Wallace-Bott problem: given a reduced
stress tensor, predict the slip a fault plane would show under it.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

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


# Below this, a shear stress is numerical noise rather than a driving force --
# the fault plane sits on, or acutely close to, a principal stress axis, and no
# slip direction is defined.
_SHEAR_MAGNITUDE_THRESHOLD = 1.0e-5

# Tolerance, in degrees, on how far from exactly orthogonal S1 and S3 may be
# and still be accepted.
_AXIS_ORTHOGONALITY_TOLERANCE = 1.0


def _axis_versor(trend_degr: float, plunge_degr: float) -> np.ndarray:
    """Trend/plunge to a unit vector in (East, North, Up).

    Same construction as `plane_normal`, at plunge `dip_angle - 90` from a
    plane's dip azimuth/angle -- the two are not independent formulas kept in
    sync by hand, they are the same trig applied to different angles.
    """
    trend, plunge = np.radians(trend_degr), np.radians(plunge_degr)
    return np.array([
        np.sin(trend) * np.cos(plunge),
        np.cos(trend) * np.cos(plunge),
        -np.sin(plunge),
    ])


def _versor_to_axis(v: Sequence[float]) -> Tuple[float, float]:
    """The trend/plunge a unit vector points along: the inverse of `_axis_versor`."""
    east, north, up = v
    plunge = float(np.degrees(np.arcsin(np.clip(-up, -1.0, 1.0))))
    trend = float(np.degrees(np.arctan2(east, north)))
    if trend < 0.0:
        trend += 360.0
    return trend, plunge


def _rhr_strike(dip_azimuth_degr: float) -> float:
    return (dip_azimuth_degr - 90.0) % 360.0


def _from_rhr_strike(strike_rhr_degr: float) -> float:
    return (strike_rhr_degr + 90.0) % 360.0


def rake_to_slickenline(
    strike_rhr_degr: float, dip_angle_degr: float, rake_degr: float
) -> Tuple[float, float]:
    """The trend/plunge of a slickenline of the given rake, on a plane of the
    given strike (right-hand rule) and dip.

    Aki & Richards (1980)'s convention: rake 0 is left-lateral, 90 reverse,
    +/-180 right-lateral, -90 normal. Unit length is an algebraic identity of
    the formula (the strike terms cancel by sin^2 + cos^2 = 1, then so do the
    dip ones), holding for every strike, dip and rake.
    """
    strike = np.radians(strike_rhr_degr)
    dip = np.radians(dip_angle_degr)
    rake = np.radians(rake_degr)

    versor = np.array([
        np.cos(rake) * np.sin(strike) - np.sin(rake) * np.cos(dip) * np.cos(strike),
        np.cos(rake) * np.cos(strike) + np.sin(rake) * np.cos(dip) * np.sin(strike),
        np.sin(rake) * np.sin(dip),
    ])

    return _versor_to_axis(versor)


def solve_stress(
    s1_trend_degr: float,
    s1_plunge_degr: float,
    s3_trend_degr: float,
    s3_plunge_degr: float,
    phi: float,
    strike_rhr_degr: float,
    dip_angle_degr: float,
    sigma1: float = 1.0,
    sigma3: float = 0.0,
    shear_threshold: float = _SHEAR_MAGNITUDE_THRESHOLD,
) -> Dict[str, object]:
    """The direct (forward) Wallace-Bott problem: resolve a reduced stress
    tensor onto one fault plane, predicting the slip it drives.

    `s1`/`s3` are the principal stress axes as (trend, plunge) in degrees,
    sub-orthogonal to within a degree; `phi` is the shape ratio
    (sigma2 - sigma3) / (sigma1 - sigma3). `sigma1`/`sigma3` default to 1/0,
    the usual normalization when only the tensor's shape is known, as from a
    fault-slip inversion -- which leaves the predicted rake correct while
    making `slip_tendency`/`deformation_index` meaningless; pass the true
    magnitudes when they are known and those are wanted.

    Returns a dict always carrying `is_valid`, `traction`, `traction_magnitude`,
    `normal_stress`, `normal_stress_magnitude`, `shear_stress` and
    `shear_stress_magnitude` (the last three vectors in (East, North, Up)),
    plus `theoretical_rake`, `theoretical_slickenline` (a `(trend, plunge)`
    pair), `slip_tendency` and `deformation_index`, `None` on all four where
    the shear stress does not clear `shear_threshold`.

    Follows Xu (2004), as `structural::stress::ReducedStressTensor::solve`
    does; see that module for how this was checked against ForwardStress.f95,
    the Fortran tool of Alberti (2010) it and this port both descend from.
    """
    s1 = _axis_versor(s1_trend_degr, s1_plunge_degr)
    s3 = _axis_versor(s3_trend_degr, s3_plunge_degr)

    angle = float(np.degrees(np.arccos(np.clip(np.dot(s1, s3), -1.0, 1.0))))
    if abs(angle - 90.0) > _AXIS_ORTHOGONALITY_TOLERANCE:
        raise ValueError(
            f"S1 and S3 axes must be sub-orthogonal: {angle:.3f} degrees apart, "
            f"expected within {_AXIS_ORTHOGONALITY_TOLERANCE:.3f} of 90"
        )
    if not (0.0 <= phi <= 1.0):
        raise ValueError(f"Phi must be between 0 and 1, got {phi}")
    if sigma1 <= sigma3:
        raise ValueError(f"Sigma1 must be greater than Sigma3, got sigma1={sigma1}, sigma3={sigma3}")

    # S3 cross S1, not S1 cross S3: this order is what makes S1, S2, S3 a
    # right-handed cyclic triad (S1 x S2 = S3, S2 x S3 = S1), matching the
    # Fortran original and the Rust port. The tensor itself does not depend
    # on which way S2 points -- it enters only as S2 (x) S2 -- but a versor
    # here still keeps the triad genuinely orthonormal for input merely
    # sub-orthogonal within the tolerance above.
    s2 = np.cross(s3, s1)
    s2 = s2 / np.linalg.norm(s2)

    sigma2 = phi * sigma1 + (1.0 - phi) * sigma3

    r = np.column_stack([s1, s2, s3])
    tensor = r @ np.diag([sigma1, sigma2, sigma3]) @ r.T

    dip_azimuth = _from_rhr_strike(strike_rhr_degr)

    # The forward-pointing normal -- upward for a shallow-dipping plane,
    # horizontal for a vertical one -- rather than the downward-pointing
    # normal_axis: the sign of every signed quantity below (traction, normal
    # stress) is relative to this choice, so it is the one the theoretical
    # rake further down is built consistently against.
    n = np.array(plane_normal(dip_azimuth, dip_angle_degr))

    traction = -(tensor @ n)
    traction_magnitude = float(np.linalg.norm(traction))
    if np.dot(n, traction) < 0.0:
        traction_magnitude = -traction_magnitude

    normal_scal = float(np.dot(n, traction))
    normal_stress = n * normal_scal
    normal_stress_magnitude = normal_scal

    shear_stress = traction - normal_stress
    shear_stress_magnitude = float(np.linalg.norm(shear_stress))

    is_valid = shear_stress_magnitude > shear_threshold

    theoretical_rake = None
    theoretical_slickenline = None
    slip_tendency = None
    deformation_index = None

    if is_valid:
        shear_versor = shear_stress / shear_stress_magnitude

        strike_versor = _axis_versor(strike_rhr_degr, 0.0)
        dip_versor = _axis_versor(dip_azimuth, dip_angle_degr)

        scal_strike = float(np.clip(np.dot(shear_versor, strike_versor), -1.0, 1.0))
        scal_dip = float(np.dot(shear_versor, dip_versor))

        rake = float(np.degrees(np.arccos(scal_strike)))
        if scal_dip > 0.0:
            rake = -rake

        theoretical_slickenline = rake_to_slickenline(strike_rhr_degr, dip_angle_degr, rake)
        theoretical_rake = rake

        slip_tendency = shear_stress_magnitude / abs(traction_magnitude)
        # abs() in the numerator alone, not the denominator: Xu (2004)'s index
        # as the Fortran original carried it, and the asymmetry is kept
        # rather than "fixed" -- see structural::stress in the Rust crate.
        deformation_index = (abs(traction_magnitude) - shear_stress_magnitude) / traction_magnitude

    return {
        "is_valid": bool(is_valid),
        "traction": tuple(float(c) for c in traction),
        "traction_magnitude": traction_magnitude,
        "normal_stress": tuple(float(c) for c in normal_stress),
        "normal_stress_magnitude": normal_stress_magnitude,
        "shear_stress": tuple(float(c) for c in shear_stress),
        "shear_stress_magnitude": shear_stress_magnitude,
        "theoretical_rake": theoretical_rake,
        "theoretical_slickenline": theoretical_slickenline,
        "slip_tendency": slip_tendency,
        "deformation_index": deformation_index,
    }
