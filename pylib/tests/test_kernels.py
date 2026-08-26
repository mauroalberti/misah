"""Check the misah bindings against the geoSurfDEM golden dataset.

The binding must not perturb what the kernel computes, so the comparison is
against a trace independently known to be correct: 404 vertices for plane 135/35
over the Malpi ASTER DEM, every one of them on the plane.

Run from anywhere but the `pylib` directory, which shadows the installed package
with its source tree:

    python3 pylib/tests/test_kernels.py
"""

from __future__ import annotations

import os
import sys

import numpy as np

DEM_PATH = os.path.expanduser(
    "~/Documenti/projects/geoSurfDEM/test_data/IntersectDEM/dem_malpi_aster_wgs84utm33.asc"
)
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


def plane_normal_reference(dip_dir, dip_angle):
    az, dip = np.radians(dip_dir), np.radians(dip_angle)
    return np.array(
        [np.sin(az) * np.sin(dip), np.cos(az) * np.sin(dip), np.cos(dip)]
    )


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

    dem, geotransform, nodata = read_esri_ascii(DEM_PATH)
    points, segments = intersect_plane_grid(
        np.ascontiguousarray(dem), geotransform, SRC_PT, DIP_DIR, DIP_ANGLE, nodata
    )

    assert len(points) == 404, len(points)
    assert len(segments) == 403, len(segments)

    residuals = np.abs((points - np.array(SRC_PT)) @ plane_normal_reference(DIP_DIR, DIP_ANGLE))
    assert residuals.max() < 1e-9, residuals.max()

    assert segments.min() >= 0 and segments.max() < len(points)


def test_vertical_and_horizontal_planes():
    from misah.kernels import intersect_plane_grid

    dem, geotransform, nodata = read_esri_ascii(DEM_PATH)

    # A vertical plane has no single-valued z; it must still work.
    points, _ = intersect_plane_grid(
        np.ascontiguousarray(dem), geotransform, SRC_PT, 135.0, 90.0, nodata
    )
    assert len(points) == 399, len(points)
    xy = points[:, :2] - points[:, :2].mean(axis=0)
    _, sv, _ = np.linalg.svd(xy, full_matrices=False)
    assert sv[1] / sv[0] < 1e-9, "the map trace of a vertical plane is a straight line"

    # A horizontal plane cuts the DEM at its own elevation, i.e. as a contour.
    level = SRC_PT[2]
    points, _ = intersect_plane_grid(
        np.ascontiguousarray(dem), geotransform, SRC_PT, 0.0, 0.0, nodata
    )
    assert len(points) == 432, len(points)
    assert np.abs(points[:, 2] - level).max() < 1e-9


def test_non_contiguous_dem_is_refused():
    """A strided view would be read as if packed; it must fail loudly."""
    from misah.kernels import intersect_plane_grid

    dem, geotransform, nodata = read_esri_ascii(DEM_PATH)
    try:
        intersect_plane_grid(
            dem[:, ::2], geotransform, SRC_PT, DIP_DIR, DIP_ANGLE, nodata
        )
    except ValueError:
        return
    raise AssertionError("non-contiguous DEM was accepted")


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
