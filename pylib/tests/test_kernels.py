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


# Five cases checked against ForwardStress.f95 itself -- the same ones, and the
# same values, structural::stress's own Rust tests and geogst's port of the
# same tool are checked against; see either for how they were obtained (the
# modules up to stress_processing compiled unmodified with gfortran, a small
# driver calling stresssolution_calc directly) and for the implicit-none bug
# in the 2010 source that separates cases 2 and 3 here from the binary as it
# stands rather than as corrected.
def test_solve_stress_case_1_anderson_normal():
    from misah.kernels import solve_stress

    sol = solve_stress(0.0, 90.0, 90.0, 0.0, 0.5, 0.0, 60.0, sigma1=30.0, sigma3=10.0)

    assert sol["is_valid"] is True
    assert np.isclose(sol["traction_magnitude"], -17.320508075689)
    assert np.isclose(sol["normal_stress_magnitude"], -15.0)
    assert np.isclose(sol["shear_stress_magnitude"], 8.660254037844)
    assert np.isclose(sol["theoretical_rake"], -90.0)
    assert np.allclose(sol["theoretical_slickenline"], (90.0, 60.0))
    assert np.isclose(sol["slip_tendency"], 0.5)
    assert np.isclose(sol["deformation_index"], -0.5)


def test_solve_stress_case_2_oblique():
    from misah.kernels import solve_stress

    sol = solve_stress(0.0, 0.0, 90.0, 0.0, 0.3, 30.0, 70.0, sigma1=40.0, sigma3=-10.0)

    assert sol["is_valid"] is True
    assert np.isclose(sol["traction_magnitude"], -20.551398971889)
    assert np.isclose(sol["normal_stress_magnitude"], -2.792444446101)
    assert np.isclose(sol["shear_stress_magnitude"], 20.360801892784)
    assert np.isclose(sol["theoretical_rake"], -2.261611727892)
    assert np.allclose(sol["theoretical_slickenline"], (30.773871690337, 2.125155259309))
    assert np.isclose(sol["slip_tendency"], 0.990725834316)
    assert np.isclose(sol["deformation_index"], -0.009274165684)


def test_solve_stress_case_3_plunging_axes():
    from misah.kernels import solve_stress

    sol = solve_stress(125.0, 35.0, 125.0, -55.0, 0.6, 200.0, 50.0, sigma1=35.0, sigma3=5.0)

    assert sol["is_valid"] is True
    assert np.isclose(sol["traction_magnitude"], -34.425635706498)
    assert np.isclose(sol["normal_stress_magnitude"], -34.215382549860)
    assert np.isclose(sol["shear_stress_magnitude"], 3.798946006885)
    assert np.isclose(sol["theoretical_rake"], 43.558963074529)
    assert np.allclose(sol["theoretical_slickenline"], (168.565012691676, -31.862445102249))
    assert np.isclose(sol["slip_tendency"], 0.110352239804)
    assert np.isclose(sol["deformation_index"], -0.889647760196)


def test_solve_stress_case_4_and_5_are_degenerate():
    """The fault normal falls on a principal axis: pure normal loading, no
    direction defined."""
    from misah.kernels import solve_stress

    sol4 = solve_stress(10.0, 5.0, 10.0, -85.0, 0.4, 100.0, 85.0, sigma1=50.0, sigma3=-5.0)
    sol5 = solve_stress(0.0, 90.0, 90.0, 0.0, 0.5, 0.0, 0.0, sigma1=30.0, sigma3=10.0)

    for sol, traction_magn in ((sol4, -50.0), (sol5, -30.0)):
        assert sol["is_valid"] is False
        assert np.isclose(sol["traction_magnitude"], traction_magn)
        assert np.isclose(sol["normal_stress_magnitude"], traction_magn)
        assert abs(sol["shear_stress_magnitude"]) < 1e-9
        assert sol["theoretical_rake"] is None
        assert sol["theoretical_slickenline"] is None
        assert sol["slip_tendency"] is None
        assert sol["deformation_index"] is None


def test_solve_stress_rejects_non_orthogonal_axes():
    from misah.kernels import solve_stress

    try:
        solve_stress(0.0, 0.0, 45.0, 0.0, 0.5, 0.0, 60.0)
    except ValueError:
        return
    raise AssertionError("non-orthogonal S1/S3 axes were accepted")


def test_rake_to_slickenline_minus_90_is_pure_normal_dip_slip():
    """Rake -90 (Aki & Richards) points straight down the dip vector."""
    from misah.kernels import rake_to_slickenline

    trend, plunge = rake_to_slickenline(0.0, 60.0, -90.0)
    assert np.isclose(trend, 90.0)
    assert np.isclose(plunge, 60.0)


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
