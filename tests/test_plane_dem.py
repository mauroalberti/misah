"""Validation of the plane-DEM kernel against the geoSurfDEM C++ golden dataset.

The reference output (``inters_malpi_135_35.csv``) was produced by the C++
IntersectDEM, which intersects the DEM against a *triangulated* planar mesh. Its
vertices therefore also fall on cell diagonals, while marching squares emits
crossings on grid edges only. The two point sets are consequently not in
one-to-one correspondence, so the traces are compared geometrically rather than
element-wise.
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

from misah_ref.plane_dem import (  # noqa: E402
    intersect_plane_dem,
    plane_normal,
    read_esri_ascii,
    signed_distance,
)

GEOSURFDEM = os.path.expanduser("~/Documenti/projects/geoSurfDEM/test_data")
DEM_PATH = os.path.join(GEOSURFDEM, "IntersectDEM/dem_malpi_aster_wgs84utm33.asc")
REF_PATH = os.path.join(GEOSURFDEM, "IntersectDEM/inters_malpi_135_35.csv")

# From malpi_plane_135_35.json: the displacement applied to the synthetic surface.
SRC_PT = (583657.237626, 4440640.5549, 1526.78799803)
DIP_DIR, DIP_ANGLE = 135.0, 35.0


def load_reference():
    return np.loadtxt(REF_PATH, delimiter=",", skiprows=1, usecols=(0, 1, 2))


def test_plane_normal_conventions():
    horizontal = plane_normal(0.0, 0.0)
    assert np.allclose(horizontal, [0.0, 0.0, 1.0])

    vertical = plane_normal(90.0, 90.0)
    assert np.allclose(vertical, [1.0, 0.0, 0.0], atol=1e-12)

    assert np.isclose(np.linalg.norm(plane_normal(135.0, 35.0)), 1.0)


def test_reference_points_lie_on_our_plane():
    """The C++ trace must satisfy our plane equation, or our conventions differ.

    1224 of the 1228 reference points sit on the plane within the centimetre the
    CSV is rounded to. The remaining four (rows 644-647, in two pairs sharing an
    elevation, up to 3.56 m off) are a defect of the C++ run, not of the plane
    definition: a plane fitted to the reference points returns 134.9999/35.0005
    and passes 8 mm from SRC_PT.
    """
    ref = load_reference()
    n = plane_normal(DIP_DIR, DIP_ANGLE)
    residuals = np.abs((ref - np.array(SRC_PT)) @ n)
    assert (residuals <= 0.01).mean() > 0.99, (residuals > 0.01).sum()
    assert (residuals > 0.01).sum() == 4, sorted(np.where(residuals > 0.01)[0])


def test_reference_best_fit_recovers_nominal_attitude():
    """Independent check that our normal convention matches geoSurfDEM's."""
    ref = load_reference()
    centred = ref - ref.mean(axis=0)
    _, _, vt = np.linalg.svd(centred, full_matrices=False)
    normal = vt[-1] * (1.0 if vt[-1][2] >= 0 else -1.0)
    angle = np.degrees(np.arccos(np.clip(normal @ plane_normal(DIP_DIR, DIP_ANGLE), -1, 1)))
    assert angle < 0.01, angle


def test_our_points_lie_on_the_plane():
    dem, gt, nodata = read_esri_ascii(DEM_PATH)
    pts, _ = intersect_plane_dem(dem, gt, SRC_PT, DIP_DIR, DIP_ANGLE, nodata)
    n = plane_normal(DIP_DIR, DIP_ANGLE)
    residuals = (pts - np.array(SRC_PT)) @ n
    assert np.abs(residuals).max() < 1e-6, np.abs(residuals).max()


def test_traces_coincide():
    """Every reference vertex must sit on our trace, and conversely."""
    dem, gt, nodata = read_esri_ascii(DEM_PATH)
    pts, _ = intersect_plane_dem(dem, gt, SRC_PT, DIP_DIR, DIP_ANGLE, nodata)
    ref = load_reference()

    cellsize = gt.pixel_width
    d_ref_to_ours = nearest_distances(ref, pts)
    d_ours_to_ref = nearest_distances(pts, ref)

    # A vertex on a cell diagonal can be up to half a diagonal from any grid-edge
    # crossing, so the tolerance is scaled to the cell, not to zero.
    assert np.median(d_ref_to_ours) < 0.5 * cellsize, np.median(d_ref_to_ours)
    assert d_ref_to_ours.max() < 1.5 * cellsize, d_ref_to_ours.max()
    assert np.median(d_ours_to_ref) < 0.5 * cellsize, np.median(d_ours_to_ref)


def nearest_distances(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Distance from each row of ``a`` to the closest row of ``b``."""
    out = np.empty(len(a))
    for i, pt in enumerate(a):
        out[i] = np.sqrt(((b - pt) ** 2).sum(axis=1)).min()
    return out


def test_segments_reference_valid_vertices():
    dem, gt, nodata = read_esri_ascii(DEM_PATH)
    pts, segs = intersect_plane_dem(dem, gt, SRC_PT, DIP_DIR, DIP_ANGLE, nodata)
    assert len(segs) > 0
    assert segs.min() >= 0
    assert segs.max() < len(pts)
    # Chords are sub-cell by construction; a longer one means mispaired edges.
    lengths = np.sqrt(((pts[segs[:, 0]] - pts[segs[:, 1]]) ** 2).sum(axis=1))
    assert lengths.max() < 3.0 * gt.pixel_width, lengths.max()


def test_horizontal_plane_is_a_contour():
    """A horizontal plane must cut the DEM exactly at its own elevation."""
    dem, gt, nodata = read_esri_ascii(DEM_PATH)
    valid = dem[dem != nodata]
    level = float(np.median(valid))
    pts, _ = intersect_plane_dem(dem, gt, (0.0, 0.0, level), 0.0, 0.0, nodata)
    assert len(pts) > 0
    assert np.abs(pts[:, 2] - level).max() < 1e-9


def test_vertical_plane_is_handled():
    """Vertical planes have no single-valued z; the signed-distance form must cope."""
    dem, gt, nodata = read_esri_ascii(DEM_PATH)
    pts, _ = intersect_plane_dem(dem, gt, SRC_PT, 135.0, 90.0, nodata)
    assert len(pts) > 0
    n = plane_normal(135.0, 90.0)
    residuals = (pts - np.array(SRC_PT)) @ n
    assert np.abs(residuals).max() < 1e-6
    # The trace of a vertical plane is a straight line in map view.
    xy = pts[:, :2] - pts[:, :2].mean(axis=0)
    _, sv, _ = np.linalg.svd(xy, full_matrices=False)
    assert sv[1] / sv[0] < 1e-9, sv


if __name__ == "__main__":
    failures = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"PASS  {name}")
            except AssertionError as exc:
                failures += 1
                print(f"FAIL  {name}: {exc}")
    print("\nfailures:", failures)
    sys.exit(1 if failures else 0)
