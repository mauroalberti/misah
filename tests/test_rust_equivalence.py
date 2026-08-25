"""Cross-validation of the Rust kernel against the numpy reference.

Builds and runs `examples/plane_dem_asc`, then compares its output with
`misah_ref` on the same input. The two implementations are meant to stay in
step, so the comparison is elementwise and tight rather than statistical.

Skipped when cargo is unavailable.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))

from misah_ref.plane_dem import intersect_plane_dem, read_esri_ascii  # noqa: E402

REPO = os.path.join(os.path.dirname(__file__), "..")
DEM_PATH = os.path.expanduser(
    "~/Documenti/projects/geoSurfDEM/test_data/IntersectDEM/dem_malpi_aster_wgs84utm33.asc"
)
SRC_PT = (583657.237626, 4440640.5549, 1526.78799803)
DIP_DIR, DIP_ANGLE = 135.0, 35.0


def run_rust(dem_path, src_pt, dip_dir, dip_angle) -> np.ndarray:
    cmd = [
        "cargo", "run", "--quiet", "--release", "--no-default-features",
        "--example", "plane_dem_asc", "--",
        dem_path, *(repr(v) for v in src_pt), repr(dip_dir), repr(dip_angle),
    ]
    proc = subprocess.run(
        cmd, cwd=REPO, capture_output=True, text=True, check=True
    )
    rows = [line.split(",") for line in proc.stdout.strip().splitlines()[1:]]
    return np.array(rows, dtype=float).reshape(-1, 3)


def test_rust_matches_reference():
    if shutil.which("cargo") is None:
        print("SKIP  cargo unavailable")
        return

    rust = run_rust(DEM_PATH, SRC_PT, DIP_DIR, DIP_ANGLE)

    dem, gt, nodata = read_esri_ascii(DEM_PATH)
    ref, _ = intersect_plane_dem(dem, gt, SRC_PT, DIP_DIR, DIP_ANGLE, nodata)

    assert len(rust) == len(ref), (len(rust), len(ref))
    # Both walk cells in row-major order and key vertices by edge, so the
    # sequences correspond one to one; only the CSV's 6 decimals separate them.
    delta = np.abs(rust - ref).max()
    assert delta < 1e-6, delta


def test_rust_matches_reference_on_a_vertical_plane():
    if shutil.which("cargo") is None:
        print("SKIP  cargo unavailable")
        return

    rust = run_rust(DEM_PATH, SRC_PT, 135.0, 90.0)

    dem, gt, nodata = read_esri_ascii(DEM_PATH)
    ref, _ = intersect_plane_dem(dem, gt, SRC_PT, 135.0, 90.0, nodata)

    assert len(rust) == len(ref), (len(rust), len(ref))
    assert np.abs(rust - ref).max() < 1e-6


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
