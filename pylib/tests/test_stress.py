"""Check the forward-stress bindings against the Fortran-verified values.

Apart from the other test file in that it needs no data: everything here is
numbers in, numbers out, so it is the one Python suite the CI can run without a
DEM from another repository. `test_kernels.py` reads the geoSurfDEM Malpi crop,
which lives outside this tree, and so still runs only by hand.

The five cases are the ones checked against ForwardStress.f95 itself, and carry
the same values `structural::stress`'s own Rust tests and geogst's port of the
same tool are checked against -- see either for how they were obtained (the
modules up to `stress_processing` compiled unmodified with gfortran, a small
driver calling `stresssolution_calc` directly) and for the `implicit none` bug
in the 2010 source that separates cases 2 and 3 here from the binary as it
stands rather than as corrected.

Run from anywhere but the `pylib` directory, which shadows the installed
package with its source tree:

    python3 pylib/tests/test_stress.py
"""

from __future__ import annotations

import sys

import numpy as np

# Keyword form, so a case reads as the geology it stands for.
STRESS_CASES = [
    dict(s1_trend_degr=0.0, s1_plunge_degr=90.0, s3_trend_degr=90.0, s3_plunge_degr=0.0,
         phi=0.5, strike_rhr_degr=0.0, dip_angle_degr=60.0, sigma1=30.0, sigma3=10.0),
    dict(s1_trend_degr=0.0, s1_plunge_degr=0.0, s3_trend_degr=90.0, s3_plunge_degr=0.0,
         phi=0.3, strike_rhr_degr=30.0, dip_angle_degr=70.0, sigma1=40.0, sigma3=-10.0),
    dict(s1_trend_degr=125.0, s1_plunge_degr=35.0, s3_trend_degr=125.0, s3_plunge_degr=-55.0,
         phi=0.6, strike_rhr_degr=200.0, dip_angle_degr=50.0, sigma1=35.0, sigma3=5.0),
    dict(s1_trend_degr=10.0, s1_plunge_degr=5.0, s3_trend_degr=10.0, s3_plunge_degr=-85.0,
         phi=0.4, strike_rhr_degr=100.0, dip_angle_degr=85.0, sigma1=50.0, sigma3=-5.0),
    dict(s1_trend_degr=0.0, s1_plunge_degr=90.0, s3_trend_degr=90.0, s3_plunge_degr=0.0,
         phi=0.5, strike_rhr_degr=0.0, dip_angle_degr=0.0, sigma1=30.0, sigma3=10.0),
]


def test_case_1_anderson_normal():
    from misah.kernels import solve_stress

    sol = solve_stress(**STRESS_CASES[0])

    assert sol["is_valid"] is True
    assert np.isclose(sol["traction_magnitude"], -17.320508075689)
    assert np.isclose(sol["normal_stress_magnitude"], -15.0)
    assert np.isclose(sol["shear_stress_magnitude"], 8.660254037844)
    assert np.isclose(sol["theoretical_rake"], -90.0)
    assert np.allclose(sol["theoretical_slickenline"], (90.0, 60.0))
    assert np.isclose(sol["slip_tendency"], 0.5)
    assert np.isclose(sol["deformation_index"], -0.5)


def test_case_2_oblique():
    from misah.kernels import solve_stress

    sol = solve_stress(**STRESS_CASES[1])

    assert sol["is_valid"] is True
    assert np.isclose(sol["traction_magnitude"], -20.551398971889)
    assert np.isclose(sol["normal_stress_magnitude"], -2.792444446101)
    assert np.isclose(sol["shear_stress_magnitude"], 20.360801892784)
    assert np.isclose(sol["theoretical_rake"], -2.261611727892)
    assert np.allclose(sol["theoretical_slickenline"], (30.773871690337, 2.125155259309))
    assert np.isclose(sol["slip_tendency"], 0.990725834316)
    assert np.isclose(sol["deformation_index"], -0.009274165684)


def test_case_3_plunging_axes():
    from misah.kernels import solve_stress

    sol = solve_stress(**STRESS_CASES[2])

    assert sol["is_valid"] is True
    assert np.isclose(sol["traction_magnitude"], -34.425635706498)
    assert np.isclose(sol["normal_stress_magnitude"], -34.215382549860)
    assert np.isclose(sol["shear_stress_magnitude"], 3.798946006885)
    assert np.isclose(sol["theoretical_rake"], 43.558963074529)
    assert np.allclose(sol["theoretical_slickenline"], (168.565012691676, -31.862445102249))
    assert np.isclose(sol["slip_tendency"], 0.110352239804)
    assert np.isclose(sol["deformation_index"], -0.889647760196)


def test_cases_4_and_5_are_degenerate():
    """The fault normal falls on a principal axis: pure normal loading, no
    direction defined."""
    from misah.kernels import solve_stress

    for case, traction_magn in ((STRESS_CASES[3], -50.0), (STRESS_CASES[4], -30.0)):
        sol = solve_stress(**case)

        assert sol["is_valid"] is False
        assert np.isclose(sol["traction_magnitude"], traction_magn)
        assert np.isclose(sol["normal_stress_magnitude"], traction_magn)
        assert abs(sol["shear_stress_magnitude"]) < 1e-9
        assert sol["theoretical_rake"] is None
        assert sol["theoretical_slickenline"] is None
        assert sol["slip_tendency"] is None
        assert sol["deformation_index"] is None


def test_non_orthogonal_axes_are_refused():
    from misah.kernels import solve_stress

    try:
        solve_stress(0.0, 0.0, 45.0, 0.0, 0.5, 0.0, 60.0)
    except ValueError:
        return
    raise AssertionError("non-orthogonal S1/S3 axes were accepted")


def test_rake_minus_90_is_pure_normal_dip_slip():
    """Rake -90 (Aki & Richards) points straight down the dip vector."""
    from misah.kernels import rake_to_slickenline

    trend, plunge = rake_to_slickenline(0.0, 60.0, -90.0)
    assert np.isclose(trend, 90.0)
    assert np.isclose(plunge, 60.0)


def test_package_exposes_the_stress_functions():
    import misah

    assert hasattr(misah.kernels, "solve_stress")
    assert hasattr(misah.kernels, "rake_to_slickenline")


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
