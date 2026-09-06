"""Check the misah raster bindings against the geoSurfDEM golden dataset.

The binding must not perturb what the kernel computes, so the comparison is
against a trace independently known to be correct: 404 vertices for plane 135/35
over the Malpi ASTER DEM, every one of them on the plane.

The grid that trace was computed on is in this repository. `example_data` ships
a 213 x 260 crop of the ASTER tile; geoSurfDEM's own crop is 200 x 247 and its
corner sits exactly seven cells in from this one's, in both directions, so the
golden grid is a *window* of what is committed here and not a second copy of it
-- the elevations in that window are identical value for value. Reading the
window is what lets this suite run in CI alongside the others, against data
whoever clones the repository already has.

The two crops carry different nodata values, -99999.0 here against -9999.9
there. No cell in the window is nodata, so the trace does not depend on which
is passed; the committed header's value is used, being the one that describes
this file.

Run from anywhere but the `pylib` directory, which shadows the installed package
with its source tree:

    python3 pylib/tests/test_kernels.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

DEM_PATH = (
    Path(__file__).resolve().parents[2]
    / "example_data" / "raster" / "monte_alpi" / "malpi_aster_w4u3.asc"
)

# The geoSurfDEM crop, in the georeferencing of its own header rather than as
# row and column indices into the file above: the window is then derived, and a
# committed DEM that no longer contains it fails here instead of silently
# yielding a trace over the wrong ground.
GOLDEN_XLLCORNER, GOLDEN_YLLCORNER = 581005.5143976583, 4437398.6422322076
GOLDEN_NCOLS, GOLDEN_NROWS = 200, 247

SRC_PT = (583657.237626, 4440640.5549, 1526.78799803)
DIP_DIR, DIP_ANGLE = 135.0, 35.0


def read_esri_ascii(path):
    header = {}
    with open(path) as src:
        for _ in range(6):
            key, value = src.readline().split()
            header[key.upper()] = float(value)
        data = np.loadtxt(src)

    nrows, ncols = int(header["NROWS"]), int(header["NCOLS"])
    cellsize = header["CELLSIZE"]
    geotransform = (
        header["XLLCORNER"], cellsize, 0.0,
        header["YLLCORNER"] + nrows * cellsize, 0.0, -cellsize,
    )
    return data.reshape(nrows, ncols), geotransform, header["NODATA_VALUE"]


def golden_grid():
    """The geoSurfDEM crop, cut out of the wider one committed here."""
    dem, geotransform, nodata = read_esri_ascii(DEM_PATH)
    x0, cellsize, _, ytop, _, _ = geotransform

    col_offset = (GOLDEN_XLLCORNER - x0) / cellsize
    yll = ytop - dem.shape[0] * cellsize
    rows_from_bottom = (GOLDEN_YLLCORNER - yll) / cellsize
    for name, offset in (("column", col_offset), ("row", rows_from_bottom)):
        assert abs(offset - round(offset)) < 1e-6, (
            f"the golden crop is not on this grid: {name} offset {offset}"
        )

    col0 = round(col_offset)
    row0 = dem.shape[0] - round(rows_from_bottom) - GOLDEN_NROWS
    assert row0 >= 0 and col0 >= 0 \
        and row0 + GOLDEN_NROWS <= dem.shape[0] \
        and col0 + GOLDEN_NCOLS <= dem.shape[1], (
            f"the committed DEM {dem.shape} does not contain the golden crop: "
            f"rows {row0}:{row0 + GOLDEN_NROWS}, columns {col0}:{col0 + GOLDEN_NCOLS}"
        )

    window = np.ascontiguousarray(
        dem[row0:row0 + GOLDEN_NROWS, col0:col0 + GOLDEN_NCOLS]
    )
    window_geotransform = (
        x0 + col0 * cellsize, cellsize, 0.0,
        ytop - row0 * cellsize, 0.0, -cellsize,
    )
    return window, window_geotransform, nodata


def plane_normal_reference(dip_dir, dip_angle):
    az, dip = np.radians(dip_dir), np.radians(dip_angle)
    return np.array(
        [np.sin(az) * np.sin(dip), np.cos(az) * np.sin(dip), np.cos(dip)]
    )


def test_golden_window_is_the_grid_the_trace_was_computed_on():
    dem, geotransform, _ = golden_grid()

    assert dem.shape == (GOLDEN_NROWS, GOLDEN_NCOLS), dem.shape
    assert np.isclose(geotransform[0], GOLDEN_XLLCORNER, atol=1e-4), geotransform[0]
    yll = geotransform[3] - GOLDEN_NROWS * geotransform[1]
    assert np.isclose(yll, GOLDEN_YLLCORNER, atol=1e-4), yll


def test_normal_conventions():
    from misah.kernels import plane_normal

    assert np.allclose(plane_normal(0.0, 0.0), [0.0, 0.0, 1.0])
    assert np.allclose(plane_normal(90.0, 90.0), [1.0, 0.0, 0.0], atol=1e-12)
    for dip_dir, dip_angle in ((135.0, 35.0), (17.0, 62.0), (300.0, 5.0)):
        assert np.allclose(
            plane_normal(dip_dir, dip_angle),
            plane_normal_reference(dip_dir, dip_angle),
        )


def test_malpi_trace_matches_the_known_result():
    from misah.kernels import intersect_plane_grid

    dem, geotransform, nodata = golden_grid()
    points, segments = intersect_plane_grid(
        dem, geotransform, SRC_PT, DIP_DIR, DIP_ANGLE, nodata
    )

    assert len(points) == 404, len(points)
    assert len(segments) == 403, len(segments)

    residuals = np.abs((points - np.array(SRC_PT)) @ plane_normal_reference(DIP_DIR, DIP_ANGLE))
    assert residuals.max() < 1e-9, residuals.max()

    assert segments.min() >= 0 and segments.max() < len(points)


def test_vertical_and_horizontal_planes():
    from misah.kernels import intersect_plane_grid

    dem, geotransform, nodata = golden_grid()

    # A vertical plane has no single-valued z; it must still work.
    points, _ = intersect_plane_grid(
        dem, geotransform, SRC_PT, 135.0, 90.0, nodata
    )
    assert len(points) == 399, len(points)
    xy = points[:, :2] - points[:, :2].mean(axis=0)
    _, sv, _ = np.linalg.svd(xy, full_matrices=False)
    assert sv[1] / sv[0] < 1e-9, "the map trace of a vertical plane is a straight line"

    # A horizontal plane cuts the DEM at its own elevation, i.e. as a contour.
    level = SRC_PT[2]
    points, _ = intersect_plane_grid(
        dem, geotransform, SRC_PT, 0.0, 0.0, nodata
    )
    assert len(points) == 432, len(points)
    assert np.abs(points[:, 2] - level).max() < 1e-9


def test_non_contiguous_dem_is_refused():
    """A strided view would be read as if packed; it must fail loudly."""
    from misah.kernels import intersect_plane_grid

    dem, geotransform, nodata = golden_grid()
    try:
        intersect_plane_grid(
            dem[:, ::2], geotransform, SRC_PT, DIP_DIR, DIP_ANGLE, nodata
        )
    except ValueError:
        return
    raise AssertionError("non-contiguous DEM was accepted")


def test_package_exposes_the_raster_kernels():
    import misah

    assert hasattr(misah, "kernels")
    for name in ("intersect_plane_grid", "plane_normal", "intersect_mesh_grid"):
        assert hasattr(misah.kernels, name), name


if __name__ == "__main__":
    if not DEM_PATH.exists():
        print(f"missing the example DEM: {DEM_PATH}")
        sys.exit(1)

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
