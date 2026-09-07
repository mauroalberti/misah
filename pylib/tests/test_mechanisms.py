"""Check the focal-mechanism bindings against Kagan (1991)'s published table.

The rotation between two double-couple mechanisms has four solutions, not one,
because a double couple is unchanged by a half turn about any of its own P, T
or B axes. Kagan (1991) works one pair through and prints all four, and that
table is what this suite checks against: numbers from a 1991 paper, not from
any code this repository has ever run.

Everything else here follows from properties rather than from stored values:

- every one of the four rotations must actually carry the first mechanism onto
  the second, which is what makes them four solutions rather than one solution
  and three artefacts;
- no two double couples can be more than 120 degrees apart, which is a
  consequence of that same fourfold symmetry and the cheapest check there is on
  an implementation -- a wrong set of symmetry generators fails it at once;
- the matrix must be symmetric with a zero diagonal, or it is not a distance
  matrix and cannot be handed to a clustering routine.

Run from anywhere but the `pylib` directory, which shadows the installed
package with its source tree:

    python3 pylib/tests/test_mechanisms.py
"""

from __future__ import annotations

import sys

import numpy as np

# Kagan (1991)'s pair, as (P trend, P plunge, T trend, T plunge) in degrees.
KAGAN_FIRST = [232.0, 41.0, 120.0, 24.0]
KAGAN_SECOND = [51.0, 17.0, 295.0, 55.0]

# His four published solutions: azimuth, colatitude measured from the downward
# vertical, and rotation angle. Left in his parameterisation rather than
# converted here, so that what is written down is what the paper prints.
KAGAN_SOLUTIONS = [
    (24.8, 101.2, 102.8),
    (257.5, 79.7, 104.3),
    (144.8, 105.2, 124.1),
    (96.8, 16.7, 165.9),
]


def random_mechanisms(count, seed=7):
    """A catalogue of proper double couples.

    Built from a T axis and a roll about it rather than from two random axes,
    which would almost never be orthogonal and would be refused.
    """
    rng = np.random.default_rng(seed)
    out = np.zeros((count, 4))

    for i in range(count):

        t_trend, t_plunge = rng.uniform(0.0, 360.0), rng.uniform(-90.0, 90.0)

        t = np.array([
            np.sin(np.radians(t_trend)) * np.cos(np.radians(t_plunge)),
            np.cos(np.radians(t_trend)) * np.cos(np.radians(t_plunge)),
            -np.sin(np.radians(t_plunge)),
        ])

        helper = np.array([0.0, 0.0, 1.0]) if abs(t[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
        u = np.cross(t, helper)
        u /= np.linalg.norm(u)
        v = np.cross(t, u)

        roll = np.radians(rng.uniform(0.0, 360.0))
        p = u * np.cos(roll) + v * np.sin(roll)

        p_plunge = np.degrees(np.arcsin(-p[2]))
        p_trend = np.degrees(np.arctan2(p[0], p[1])) % 360.0

        out[i] = [p_trend, p_plunge, t_trend, t_plunge]

    return out


def test_the_four_rotations_are_the_ones_kagan_published():
    from misah.kernels import focal_mechanism_rotations

    found = focal_mechanism_rotations(*KAGAN_FIRST, *KAGAN_SECOND)

    assert len(found["angle_degrees"]) == 4

    for i, (azimuth, colatitude, angle) in enumerate(KAGAN_SOLUTIONS):

        trend, plunge = found["trend"][i], found["plunge"][i]
        turn = found["angle_degrees"][i]

        # Kagan gives colatitude from the downward vertical; the plunge is its
        # complement. An axis and its opposite end describe the same rotation,
        # so either representation is accepted.
        expected_plunge = 90.0 - colatitude

        deviation = min(
            max(abs((trend - t + 180.0) % 360.0 - 180.0), abs(plunge - p), abs(turn - a))
            for t, p, a in [
                (azimuth, expected_plunge, angle),
                (azimuth + 180.0, -expected_plunge, -angle),
            ]
        )

        # The published values are quoted to a tenth of a degree, and the worst
        # component deviation over the four is 0.20.
        assert deviation < 0.25, (
            f"solution {i} ({trend:.3f}, {plunge:.3f}, {turn:.3f}) departs from "
            f"Kagan's ({azimuth}, {expected_plunge}, {angle}) by {deviation:.4f} degrees"
        )


def test_every_rotation_really_carries_one_mechanism_onto_the_other():
    from misah.kernels import focal_mechanism_rotations, rotate_focal_mechanism

    found = focal_mechanism_rotations(*KAGAN_FIRST, *KAGAN_SECOND)

    for trend, plunge, angle in zip(
        found["trend"], found["plunge"], found["angle_degrees"]
    ):
        turned = rotate_focal_mechanism(*KAGAN_FIRST, trend, plunge, angle)

        # Against the target's own axes, allowing the third of a degree by
        # which Kagan's P and T are not exactly orthogonal and are made so.
        assert abs(turned["t"][0] - KAGAN_SECOND[2]) < 0.1, turned
        assert abs(turned["t"][1] - KAGAN_SECOND[3]) < 0.1, turned
        assert abs(turned["p"][0] - KAGAN_SECOND[0]) < 0.1, turned
        assert abs(turned["p"][1] - KAGAN_SECOND[1]) < 0.1, turned


def test_the_kagan_angle_is_the_smallest_of_the_four():
    from misah.kernels import focal_mechanism_rotations, kagan_angles

    found = focal_mechanism_rotations(*KAGAN_FIRST, *KAGAN_SECOND)
    angles = np.abs(found["angle_degrees"])

    # Sorted smallest first.
    assert (np.diff(angles) >= 0).all(), angles

    elementwise = kagan_angles(np.array([KAGAN_FIRST]), np.array([KAGAN_SECOND]))

    assert abs(elementwise[0] - angles[0]) < 1e-12
    assert abs(elementwise[0] - 102.8) < 0.25, elementwise


def test_no_pair_of_double_couples_is_more_than_120_degrees_apart():
    from misah.kernels import kagan_angle_matrix

    # 120 and not 180: the fourfold symmetry is what bounds it, and a wrong set
    # of generators fails this immediately.
    matrix = kagan_angle_matrix(random_mechanisms(120))

    assert np.isfinite(matrix).all()
    assert matrix.max() <= 120.0 + 1e-9, matrix.max()
    # And the bound is approached, or this would pass on a function returning
    # zeroes.
    assert matrix.max() > 110.0, matrix.max()


def test_the_matrix_is_a_distance_matrix():
    from misah.kernels import kagan_angle_matrix, kagan_angles

    mechanisms = random_mechanisms(40)
    matrix = kagan_angle_matrix(mechanisms)

    assert matrix.shape == (40, 40)
    # Symmetric with a zero diagonal, or it cannot go into a clustering
    # routine, which is the reason to want one.
    assert np.array_equal(matrix, matrix.T)
    assert (np.diag(matrix) == 0.0).all()

    # And its entries are the same numbers the pairwise call gives, which is
    # what says the mirroring did not transpose anything.
    first = np.repeat(mechanisms[:5], 40, axis=0)
    second = np.tile(mechanisms, (5, 1))
    pairwise = kagan_angles(first, second).reshape(5, 40)

    assert np.allclose(pairwise, matrix[:5], atol=1e-12)


def test_a_mechanism_is_no_distance_from_itself():
    from misah.kernels import focal_mechanism_rotations, kagan_angles

    same = kagan_angles(np.array([KAGAN_FIRST]), np.array([KAGAN_FIRST]))
    assert same[0] < 1e-6, same

    # The other three solutions are the half turns a double couple is
    # genuinely invariant under, and are reported rather than suppressed.
    found = focal_mechanism_rotations(*KAGAN_FIRST, *KAGAN_FIRST)
    rest = np.abs(found["angle_degrees"][1:])

    assert np.allclose(rest, 180.0, atol=1e-6), found["angle_degrees"]


def test_p_and_t_come_out_of_a_fault_and_its_slip():
    from misah.kernels import ptb_axes

    # A plane striking north and dipping 45 east, slipping down-dip: extension
    # east-west, so T is horizontal and east and P is vertical.
    faults = np.array([[0.0, 45.0, 90.0, 45.0]])
    axes = ptb_axes(faults)

    assert axes["p"].shape == (1, 2)
    assert axes["t"].shape == (1, 2)
    assert axes["b"].shape == (1, 2)

    assert abs(axes["t"][0][0] - 90.0) < 1e-9, axes["t"]
    assert abs(axes["t"][0][1]) < 1e-9, axes["t"]
    assert abs(axes["p"][0][1] - 90.0) < 1e-9, axes["p"]

    # B is the null axis, horizontal and along strike here. Its trend must be
    # a bearing: 0 and not 360, which is what a half-open range means.
    assert 0.0 <= axes["b"][0][0] < 360.0, axes["b"]


def test_a_slip_of_unknown_sense_is_refused_rather_than_guessed():
    from misah.kernels import ptb_axes

    faults = np.array([[0.0, 45.0, 90.0, 45.0]])

    # Reversing an undetermined slip exchanges P with T, turning shortening
    # into extension. Stricter than invert_stress, which handles the same
    # missing information by taking its misfit modulo 180 -- an option that
    # does not exist here.
    try:
        ptb_axes(faults, senses=np.array([False]))
    except ValueError as exc:
        assert "sense" in str(exc), exc
    else:
        raise AssertionError("a slickenline of unknown sense was accepted")


def test_axes_that_are_not_orthogonal_are_named_and_refused():
    from misah.kernels import kagan_angles

    # Ten degrees off is not a triad with a small error in it, and in a
    # catalogue of a thousand the row number is the difference between a
    # fixable complaint and a search.
    bad = np.array([[0.0, 0.0, 80.0, 0.0]])

    try:
        kagan_angles(bad, np.array([KAGAN_SECOND]))
    except ValueError as exc:
        assert "mechanism 0" in str(exc), exc
    else:
        raise AssertionError("a non-orthogonal pair of axes was accepted")


def test_mismatched_catalogues_are_refused():
    from misah.kernels import kagan_angles

    try:
        kagan_angles(random_mechanisms(3), random_mechanisms(4, seed=8))
    except ValueError as exc:
        assert "same length" in str(exc), exc
    else:
        raise AssertionError("two catalogues of different length were accepted")


def test_package_exposes_the_mechanism_functions():
    import misah

    for name in ("ptb_axes", "kagan_angles", "kagan_angle_matrix",
                 "focal_mechanism_rotations", "rotate_focal_mechanism"):
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
