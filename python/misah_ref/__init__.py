"""Numpy reference implementations of the misah kernels.

These mirror `src/kernels/` and act both as the oracle the Rust ports are
validated against and as the fallback where a compiled extension cannot be
installed.
"""

from .plane_dem import (
    GeoTransform,
    intersect_plane_dem,
    plane_normal,
    read_esri_ascii,
    signed_distance,
)

__all__ = [
    "GeoTransform",
    "intersect_plane_dem",
    "plane_normal",
    "read_esri_ascii",
    "signed_distance",
]
