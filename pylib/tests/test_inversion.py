"""Check the fault-slip inversion bindings, and the two that go with them, by
inverting what the forward ones predict.

`invert_stress` and `inversion_candidate_count` are the search; `score_stress`
evaluates one tensor without searching for it, and `stress_tensor` is the
matrix behind the axes-and-Phi form the other three speak in. The last is
checked by taking it apart again rather than against stored numbers -- a stress
tensor is symmetric, and its principal axes are its eigenvectors.

Numbers in, numbers out like `test_stress.py`, so this runs in CI: the faults
are built here, through `solve_stress`, from a tensor chosen in advance, and
the inversion is then asked to find that tensor again. That is the one case
where the right answer is known independently of the search -- and it exercises
both halves of the Python surface against each other, the forward model
producing the data the inverse model consumes, rather than checking either
against numbers copied from the Rust tests.

The data are generated without error, so the true tensor sits exactly on the
default grid and the misfit at it is zero. Anything above one grid step means
the search found a different tensor, not a nearby one.

Run from anywhere but the `pylib` directory, which shadows the installed
package with its source tree:

    python3 pylib/tests/test_inversion.py
"""

from __future__ import annotations

import sys

import numpy as np

# A spread of orientations, so the inversion is constrained rather than left
# free by a dataset that all dips one way. The same eight the Rust tests use.
SPREAD = [
    (0.0, 60.0), (45.0, 70.0), (90.0, 50.0), (135.0, 65.0),
    (180.0, 55.0), (225.0, 75.0), (270.0, 45.0), (315.0, 60.0),
]

# Andersonian normal faulting: S1 vertical, S3 horizontal east-west.
NORMAL_TENSOR = (0.0, 90.0, 90.0, 0.0, 0.5)

# S1 and S3 both horizontal: a setting the vertical-S1 case could pass without
# ever exercising.
STRIKE_SLIP_TENSOR = (0.0, 0.0, 90.0, 0.0, 0.5)


def faults_from(tensor, planes=SPREAD):
    """The (N, 4) array of faults a tensor would drive on the given planes.

    Planes the tensor predicts no slip on -- lying on a principal axis, where
    the shear vanishes -- are dropped rather than passed with a missing
    lineation, there being nothing to invert them against.
    """
    from misah.kernels import solve_stress

    rows = []
    for strike, dip in planes:
        solution = solve_stress(*tensor, strike, dip)
        slickenline = solution["theoretical_slickenline"]
        if slickenline is not None:
            rows.append((strike, dip, slickenline[0], slickenline[1]))

    return np.array(rows, dtype=np.float64)


def versor(trend_plunge):
    """(East, North, Up), as `GeologicalAxis::as_versor` builds it."""
    trend, plunge = np.radians(trend_plunge)
    return np.array([
        np.sin(trend) * np.cos(plunge),
        np.cos(trend) * np.cos(plunge),
        -np.sin(plunge),
    ])


def axes_apart(a, b):
    """Degrees between two axes, which an antipode does not change."""
    cosine = abs(float(np.dot(versor(a), versor(b))))
    return float(np.degrees(np.arccos(min(cosine, 1.0))))


def test_a_normal_tensor_is_recovered_from_the_slip_it_predicts():
    from misah.kernels import invert_stress

    faults = faults_from(NORMAL_TENSOR)
    assert len(faults) >= 6, "the synthetic set is too small to constrain anything"

    result = invert_stress(faults)

    assert result is not None
    assert result["best"]["mean_misfit_degrees"] < 10.0
    assert result["best"]["faults_scored"] == len(faults)
    assert axes_apart(result["best"]["s1"], NORMAL_TENSOR[:2]) <= 10.0 + 1e-6
    assert np.isclose(result["best"]["phi"], 0.5, atol=0.1)


def test_a_strike_slip_tensor_is_recovered_too():
    from misah.kernels import invert_stress

    result = invert_stress(faults_from(STRIKE_SLIP_TENSOR))

    assert result is not None
    assert result["best"]["mean_misfit_degrees"] < 10.0
    assert axes_apart(result["best"]["s1"], STRIKE_SLIP_TENSOR[:2]) <= 10.0 + 1e-6


def test_the_principal_axes_come_back_orthogonal():
    from misah.kernels import invert_stress

    best = invert_stress(faults_from(NORMAL_TENSOR))["best"]

    for one, other in (("s1", "s2"), ("s2", "s3"), ("s3", "s1")):
        apart = axes_apart(best[one], best[other])
        assert abs(apart - 90.0) < 1e-6, f"{one} and {other} are {apart} degrees apart"


def test_the_result_reports_what_it_searched():
    from misah.kernels import inversion_candidate_count, invert_stress

    result = invert_stress(faults_from(NORMAL_TENSOR), angle_step_degrees=30.0, phi_step=0.5)

    assert result["candidates_scored"] <= result["candidates_tried"]
    assert result["candidates_tried"] <= inversion_candidate_count(30.0, 0.5)

    # Sorted, worst last: every runner-up fits no better than the winner.
    misfits = [r["mean_misfit_degrees"] for r in result["runners_up"]]
    assert misfits == sorted(misfits)
    assert misfits[0] >= result["best"]["mean_misfit_degrees"] - 1e-9


def test_a_finer_grid_does_not_fit_worse():
    from misah.kernels import invert_stress

    faults = faults_from(NORMAL_TENSOR)

    coarse = invert_stress(faults, angle_step_degrees=30.0, phi_step=0.5)
    fine = invert_stress(faults)

    assert fine["best"]["mean_misfit_degrees"] <= coarse["best"]["mean_misfit_degrees"] + 1e-9
    assert fine["candidates_tried"] > coarse["candidates_tried"]


def test_an_unread_sense_is_not_penalised():
    """Half the lineations stored pointing the wrong way along their own line,
    which is exactly what an undetermined sense means: with `senses` saying so,
    they still fit perfectly."""
    from misah.kernels import invert_stress

    faults = faults_from(NORMAL_TENSOR)

    flipped = faults.copy()
    flipped[::2, 2] = (flipped[::2, 2] + 180.0) % 360.0
    flipped[::2, 3] = -flipped[::2, 3]

    senses = np.zeros(len(flipped), dtype=bool)
    unread = invert_stress(flipped, senses=senses)
    assert unread["best"]["mean_misfit_degrees"] < 1e-6

    # Claiming those senses were read is what the flag exists to prevent: the
    # reversed half is then scored as wrong rather than as undirected.
    claimed = invert_stress(flipped)
    assert claimed["best"]["mean_misfit_degrees"] > 10.0


def test_a_lineation_off_its_plane_is_refused_by_row():
    from misah.kernels import invert_stress

    faults = faults_from(NORMAL_TENSOR)
    # Row 2's lineation moved well off the plane it is recorded on: not a
    # fault measured imprecisely, but two measurements that do not go together.
    faults[2, 2] += 40.0

    try:
        invert_stress(faults)
    except ValueError as exc:
        assert "fault 2" in str(exc), f"the offending row is not named: {exc}"
        return
    raise AssertionError("a slickenline off its plane was accepted")


def test_a_slickenline_within_measurement_error_is_accepted():
    """A degree is what a pair of compass readings of the same fault is worth,
    so the tolerance has to admit one."""
    from misah.kernels import invert_stress

    faults = faults_from(NORMAL_TENSOR)
    faults[:, 3] += 0.4

    assert invert_stress(faults) is not None


def test_an_empty_set_inverts_to_nothing():
    from misah.kernels import invert_stress

    assert invert_stress(np.zeros((0, 4), dtype=np.float64)) is None


def test_the_input_shape_is_checked():
    from misah.kernels import invert_stress

    for bad in (np.zeros((3, 3)), np.zeros((3, 5))):
        try:
            invert_stress(bad)
        except ValueError:
            continue
        raise AssertionError(f"an array of shape {bad.shape} was accepted")


def test_senses_must_match_the_faults():
    from misah.kernels import invert_stress

    faults = faults_from(NORMAL_TENSOR)

    try:
        invert_stress(faults, senses=np.ones(len(faults) - 1, dtype=bool))
    except ValueError:
        return
    raise AssertionError("a senses array of the wrong length was accepted")


def test_a_zero_grid_step_is_refused_rather_than_walked():
    """Both loops advance by adding their step, so a zero would not return a
    coarse answer -- it would not return."""
    from misah.kernels import inversion_candidate_count, invert_stress

    faults = faults_from(NORMAL_TENSOR)

    for kwargs in (dict(angle_step_degrees=0.0), dict(phi_step=0.0),
                   dict(angle_step_degrees=float("nan"))):
        try:
            invert_stress(faults, **kwargs)
        except ValueError:
            pass
        else:
            raise AssertionError(f"{kwargs} was accepted")

    try:
        inversion_candidate_count(0.0, 0.1)
    except ValueError:
        return
    raise AssertionError("a zero step was accepted for the candidate count")


def test_the_candidate_count_grows_as_the_grid_is_refined():
    from misah.kernels import inversion_candidate_count

    assert inversion_candidate_count() == inversion_candidate_count(10.0, 0.1)
    assert inversion_candidate_count(5.0, 0.1) > inversion_candidate_count(10.0, 0.1)
    assert inversion_candidate_count(10.0, 0.05) > inversion_candidate_count(10.0, 0.1)


def test_the_tensor_matrix_is_symmetric_with_the_right_eigenpairs():
    """`R . diag(sigma) . R^T` is checked by taking it apart again: a stress
    tensor is symmetric, and its principal axes are its eigenvectors with the
    principal magnitudes for eigenvalues."""
    from misah.kernels import stress_tensor

    matrix = stress_tensor(*NORMAL_TENSOR, sigma1=30.0, sigma3=10.0)

    assert matrix.shape == (3, 3)
    assert np.allclose(matrix, matrix.T)

    # sigma2 = phi * sigma1 + (1 - phi) * sigma3, so 20 at phi 0.5.
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    assert np.allclose(sorted(eigenvalues), [10.0, 20.0, 30.0])

    for magnitude, axis in ((30.0, NORMAL_TENSOR[:2]), (10.0, NORMAL_TENSOR[2:4])):
        column = eigenvectors[:, np.argmin(abs(eigenvalues - magnitude))]
        # An eigenvector is defined up to sign, which is what an axis is too.
        assert abs(abs(np.dot(column, versor(axis))) - 1.0) < 1e-9


def test_the_tensor_scales_with_its_magnitudes():
    """Unlike the misfit, the matrix is not shape alone: the default 1/0 gives
    the normalized tensor `invert_stress` searches in."""
    from misah.kernels import stress_tensor

    normalized = stress_tensor(*NORMAL_TENSOR)
    assert np.allclose(sorted(np.linalg.eigvalsh(normalized)), [0.0, 0.5, 1.0])

    # An affine map of the magnitudes is an affine map of the tensor: scaled by
    # the range, shifted along the identity. This is also why the misfit cannot
    # depend on them -- the identity part loads the plane normally, and a
    # purely normal traction leaves the shear direction alone.
    scaled = stress_tensor(*NORMAL_TENSOR, sigma1=30.0, sigma3=10.0)
    assert np.allclose(scaled, 20.0 * normalized + 10.0 * np.eye(3))


def test_non_orthogonal_axes_are_refused_by_both_new_functions():
    from misah.kernels import score_stress, stress_tensor

    for call in (lambda: stress_tensor(0.0, 0.0, 45.0, 0.0, 0.5),
                 lambda: score_stress(0.0, 0.0, 45.0, 0.0, 0.5, faults_from(NORMAL_TENSOR))):
        try:
            call()
        except ValueError:
            continue
        raise AssertionError("non-orthogonal S1/S3 axes were accepted")


def test_scoring_the_generating_tensor_costs_nothing():
    from misah.kernels import score_stress

    faults = faults_from(NORMAL_TENSOR)
    score = score_stress(*NORMAL_TENSOR, faults)

    assert score["faults_scored"] == len(faults)
    assert score["misfits"].shape == (len(faults),)
    assert np.all(score["misfits"] < 1e-9)
    assert score["mean_misfit_degrees"] < 1e-9


def test_scoring_agrees_with_the_search_that_found_the_tensor():
    """The binding averages the per-fault misfits itself, to hand them back;
    this is what says that average is still `mean_misfit`'s."""
    from misah.kernels import invert_stress, score_stress

    faults = faults_from(STRIKE_SLIP_TENSOR)
    best = invert_stress(faults)["best"]

    score = score_stress(*best["s1"], *best["s3"], best["phi"], faults)

    assert score["faults_scored"] == best["faults_scored"]
    assert np.isclose(score["mean_misfit_degrees"], best["mean_misfit_degrees"])


def test_a_wrong_tensor_scores_worse_than_the_right_one():
    from misah.kernels import score_stress

    faults = faults_from(NORMAL_TENSOR)

    right = score_stress(*NORMAL_TENSOR, faults)["mean_misfit_degrees"]
    wrong = score_stress(*STRIKE_SLIP_TENSOR, faults)["mean_misfit_degrees"]

    assert wrong > right + 10.0


def test_a_fault_the_tensor_cannot_speak_about_is_nan_not_zero():
    """A plane normal to a principal axis carries no shear, so no slip is
    predicted -- which is not the same as a slip predicted and matched."""
    from misah.kernels import score_stress

    faults = faults_from(NORMAL_TENSOR)

    # S1 is vertical here, so a horizontal plane has its normal on it. Its
    # slickenline is arbitrary: nothing is predicted to compare it against.
    horizontal = np.array([[0.0, 0.0, 90.0, 0.0]])
    mixed = np.vstack([faults, horizontal])

    score = score_stress(*NORMAL_TENSOR, mixed)

    assert np.isnan(score["misfits"][-1])
    assert score["faults_scored"] == len(faults)
    assert score["mean_misfit_degrees"] < 1e-9


def test_scoring_an_unscorable_set_gives_no_mean():
    from misah.kernels import score_stress

    score = score_stress(*NORMAL_TENSOR, np.array([[0.0, 0.0, 90.0, 0.0]]))

    assert score["mean_misfit_degrees"] is None
    assert score["faults_scored"] == 0
    assert np.isnan(score["misfits"]).all()


def test_scoring_takes_the_sense_flag_too():
    from misah.kernels import score_stress

    faults = faults_from(NORMAL_TENSOR)
    flipped = faults.copy()
    flipped[:, 2] = (flipped[:, 2] + 180.0) % 360.0
    flipped[:, 3] = -flipped[:, 3]

    unread = score_stress(*NORMAL_TENSOR, flipped, senses=np.zeros(len(flipped), dtype=bool))
    claimed = score_stress(*NORMAL_TENSOR, flipped)

    assert unread["mean_misfit_degrees"] < 1e-9
    assert np.isclose(claimed["mean_misfit_degrees"], 180.0)


def test_package_exposes_the_inversion_functions():
    import misah

    for name in ("invert_stress", "inversion_candidate_count",
                 "stress_tensor", "score_stress"):
        assert hasattr(misah.kernels, name), name


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
