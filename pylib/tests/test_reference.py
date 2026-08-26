"""Check the pure-Python fallback against the compiled kernel.

The two are meant to be interchangeable, so the comparison is elementwise and
tight rather than statistical: same vertex count, same order, same segments.

Run from anywhere but the `pylib` directory, which shadows the installed package
with its source tree:

    python3 pylib/tests/test_reference.py
"""

from __future__ import annotations

import os
import sys

import numpy as np

DEM_PATH = os.path.expanduser(
    "~/Documenti/projects/geoSurfDEM/test_data/IntersectDEM/dem_malpi_aster_wgs84utm33.asc"
)
SRC_PT = (583657.237626, 4440640.5549, 1526.78799803)

ORIENTATIONS = ((135.0, 35.0), (135.0, 90.0), (0.0, 0.0), (17.0, 62.0))


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


def _both():
    """The compiled kernel and the fallback, or None when only one is available."""
    from misah import _reference

    try:
        from misah._misah import kernels as compiled
    except ImportError:
        return None, _reference
    return compiled, _reference


def test_normals_agree():
    compiled, reference = _both()
    if compiled is None:
        print("SKIP  extension not installed")
        return

    for dip_dir, dip_angle in ORIENTATIONS:
        assert np.allclose(
            compiled.plane_normal(dip_dir, dip_angle),
            reference.plane_normal(dip_dir, dip_angle),
            atol=1e-15,
        ), (dip_dir, dip_angle)


def test_traces_agree_on_the_malpi_dem():
    compiled, reference = _both()
    if compiled is None:
        print("SKIP  extension not installed")
        return

    dem, geotransform, nodata = read_esri_ascii(DEM_PATH)
    dem = np.ascontiguousarray(dem)

    for dip_dir, dip_angle in ORIENTATIONS:
        pts_c, segs_c = compiled.intersect_plane_grid(
            dem, geotransform, SRC_PT, dip_dir, dip_angle, nodata
        )
        pts_r, segs_r = reference.intersect_plane_grid(
            dem, geotransform, SRC_PT, dip_dir, dip_angle, nodata
        )

        assert pts_c.shape == pts_r.shape, (dip_dir, dip_angle, pts_c.shape, pts_r.shape)
        # Both walk cells in row-major order and key vertices by edge, so the
        # sequences correspond one to one.
        assert np.abs(pts_c - pts_r).max() < 1e-9, (dip_dir, dip_angle)
        assert np.array_equal(segs_c, segs_r), (dip_dir, dip_angle)


def test_fallback_reproduces_the_known_malpi_result():
    """The fallback alone must reach the figures the C++ trace established."""
    from misah import _reference

    dem, geotransform, nodata = read_esri_ascii(DEM_PATH)
    dem = np.ascontiguousarray(dem)

    points, segments = _reference.intersect_plane_grid(
        dem, geotransform, SRC_PT, 135.0, 35.0, nodata
    )
    assert len(points) == 404, len(points)
    assert len(segments) == 403, len(segments)

    normal = np.array(_reference.plane_normal(135.0, 35.0))
    residuals = np.abs((points - np.array(SRC_PT)) @ normal)
    assert residuals.max() < 1e-9, residuals.max()


def test_fallback_refuses_a_non_contiguous_dem():
    """It must reject what the extension rejects, not quietly accept it."""
    from misah import _reference

    dem, geotransform, nodata = read_esri_ascii(DEM_PATH)
    try:
        _reference.intersect_plane_grid(
            dem[:, ::2], geotransform, SRC_PT, 135.0, 35.0, nodata
        )
    except ValueError:
        return
    raise AssertionError("non-contiguous DEM was accepted")


def test_package_exposes_a_kernels_module_either_way():
    import misah

    assert hasattr(misah, "kernels")
    assert hasattr(misah.kernels, "intersect_plane_grid")
    assert hasattr(misah.kernels, "plane_normal")
    assert isinstance(misah.is_compiled, bool)


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
