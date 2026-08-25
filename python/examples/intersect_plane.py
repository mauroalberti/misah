"""Intersect a geological plane with a DEM and write the trace to CSV.

Run it as:

    python3 python/examples/intersect_plane.py

The kernel comes from the compiled extension when it is installed and from the
numpy reference otherwise, so the same script works inside QGIS, where a binary
wheel may not be installable.
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

DEM_PATH = os.path.expanduser(
    "~/Documenti/projects/geoSurfDEM/test_data/IntersectDEM/dem_malpi_aster_wgs84utm33.asc"
)
SRC_PT = (583657.237626, 4440640.5549, 1526.78799803)
DIP_DIR, DIP_ANGLE = 135.0, 35.0
OUT_PATH = "/tmp/trace_135_35.csv"


def load_dem(path):
    """Return (elevations, GDAL geotransform, nodata).

    Any reader will do -- rasterio and GDAL both hand back the six-element
    geotransform the kernel expects. The ESRI ASCII reader in `misah_ref` is used
    here only to keep the example free of a raster dependency.
    """
    try:
        import rasterio
    except ImportError:
        from misah_ref.plane_dem import read_esri_ascii

        dem, gt, nodata = read_esri_ascii(path)
        return dem, (
            gt.x_origin, gt.pixel_width, gt.row_rotation,
            gt.y_origin, gt.col_rotation, gt.pixel_height,
        ), nodata

    with rasterio.open(path) as src:
        dem = src.read(1).astype(float)
        t = src.transform
        # rasterio orders the affine as (a, b, c, d, e, f); GDAL wants
        # (c, a, b, f, d, e) -- origin first, then the two rows of the matrix.
        return dem, (t.c, t.a, t.b, t.f, t.d, t.e), src.nodata


def intersect(dem, geotransform, src_pt, dip_dir, dip_angle, nodata):
    """Call the compiled kernel, falling back to the numpy reference."""
    try:
        from misah import kernels
    except ImportError:
        from misah_ref.plane_dem import GeoTransform, intersect_plane_dem

        gt = GeoTransform(*geotransform)
        print("kernel: numpy reference (misah extension not installed)")
        return intersect_plane_dem(dem, gt, src_pt, dip_dir, dip_angle, nodata)

    print("kernel: compiled misah extension")
    # The kernel reads the grid as packed memory, so hand it a contiguous array.
    return kernels.intersect_plane_dem(
        np.ascontiguousarray(dem, dtype=float),
        geotransform,
        src_pt,
        dip_dir,
        dip_angle,
        nodata,
    )


def main():
    dem, geotransform, nodata = load_dem(DEM_PATH)
    print(f"DEM {dem.shape[0]}x{dem.shape[1]}, nodata={nodata}")

    points, segments = intersect(
        dem, geotransform, SRC_PT, DIP_DIR, DIP_ANGLE, nodata
    )
    print(f"{len(points)} intersection points, {len(segments)} segments")

    if len(points) == 0:
        print("the plane does not meet the DEM")
        return

    print(
        "elevation range: {:.1f} to {:.1f} m".format(
            points[:, 2].min(), points[:, 2].max()
        )
    )

    # Residuals confirm the points really do lie on the plane; they should be at
    # the level of floating-point noise.
    normal = np.array(
        [
            np.sin(np.radians(DIP_DIR)) * np.sin(np.radians(DIP_ANGLE)),
            np.cos(np.radians(DIP_DIR)) * np.sin(np.radians(DIP_ANGLE)),
            np.cos(np.radians(DIP_ANGLE)),
        ]
    )
    residuals = np.abs((points - np.array(SRC_PT)) @ normal)
    print(f"max distance from the plane: {residuals.max():.2e} m")

    np.savetxt(
        OUT_PATH,
        np.column_stack([points, np.full(len(points), DIP_DIR), np.full(len(points), DIP_ANGLE)]),
        delimiter=",",
        header="x,y,z,dipdir,dipang",
        comments="",
        fmt="%.3f",
    )
    print(f"written: {OUT_PATH}")


if __name__ == "__main__":
    main()
