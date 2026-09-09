"""Check the field bindings: a density on a grid, and a tensor at every node.

The properties tested here are the ones the boundary could break without the
Rust suite noticing, because they are properties *of the boundary*: that the
dimension is read off the data rather than assumed, that the flat order coming
back is the order the grid was defined in, that the arrays in the returned dict
are aligned with each other, and that a node without a tensor says so rather
than carrying a plausible zero.

Two of them are worth stating as claims rather than as assertions, since they
are what the whole exercise is for.

The first is that a density is a density. A kernel that integrates to one means
the field summed over the grid, times the cell volume, comes back as the number
of observations -- so the test is the count, not a stored array. That also
catches the error this lineage actually made: a normalizing constant right for
two dimensions used in three is off by a factor, and the integral is where a
factor shows.

The second is that a field finds locally what a single inversion cannot find at
all. Two tectonic phases are placed far enough apart that a compact kernel never
sees both; one inversion over the pair returns a tensor belonging to neither,
and the field returns each in its own place. Everything is built here from
`solve_stress`, so the right answer is known independently of the search rather
than copied from the Rust tests.

Run from anywhere but the `pylib` directory, which shadows the installed
package with its source tree:

    python3 pylib/tests/test_fields.py
"""

from __future__ import annotations

import sys

import numpy as np

# The same eight orientations the other suites use: a spread, so an inversion
# is constrained rather than left free by a dataset that all dips one way.
SPREAD = [
    (0.0, 60.0), (45.0, 70.0), (90.0, 50.0), (135.0, 65.0),
    (180.0, 55.0), (225.0, 75.0), (270.0, 45.0), (315.0, 60.0),
]

# Andersonian normal faulting: S1 vertical.
NORMAL_TENSOR = (0.0, 90.0, 90.0, 0.0, 0.5)
# S1 and S3 both horizontal, which the vertical-S1 case could pass without
# ever exercising.
STRIKE_SLIP_TENSOR = (0.0, 0.0, 90.0, 0.0, 0.5)


def faults_from(tensor, planes=SPREAD):
    """The (N, 4) array of faults a tensor would drive on the given planes."""
    from misah.kernels import solve_stress

    rows = []
    for strike, dip in planes:
        solution = solve_stress(*tensor, strike, dip)
        slickenline = solution["theoretical_slickenline"]
        if slickenline is None:
            # The plane lies on a principal axis; no slip is predicted and
            # there is nothing to invert it against.
            continue
        rows.append([strike, dip, slickenline[0], slickenline[1]])

    return np.array(rows)


def two_phases(separation=10_000.0):
    """Two fault sets at two places, obeying two different tensors."""

    west = faults_from(NORMAL_TENSOR)
    east = faults_from(STRIKE_SLIP_TENSOR)

    faults = np.vstack([west, east])
    positions = np.vstack([
        np.zeros((len(west), 2)),
        np.tile([separation, 0.0], (len(east), 1)),
    ])

    return positions, faults


def clustered_points(count=400, seed=20260907):
    """A blob of points in three dimensions, away from the origin.

    Away from it deliberately: a grid whose origin is zero would pass a test
    that a grid built from the data would fail, and the offset is what makes
    the registration mean something.
    """
    rng = np.random.default_rng(seed)

    return rng.normal(0.0, 300.0, size=(count, 3)) + np.array([12_000.0, 34_000.0, -5_000.0])


def test_a_covering_grid_holds_the_points_and_starts_on_a_node():
    from misah.kernels import covering_grid

    points = clustered_points()
    grid = covering_grid(points, [100.0, 100.0, 100.0], 1_000.0)

    origin = np.array(grid["origin"])
    spacing = np.array(grid["spacing"])
    counts = np.array(grid["counts"])

    # The origin is the first sample point, not the corner of a cell around
    # it, so it sits exactly one margin below the lowest observation. Half a
    # spacing either way is the classic silent error here, and this is where it
    # would show.
    assert np.allclose(origin, points.min(axis=0) - 1_000.0), origin

    last = origin + (counts - 1) * spacing
    assert (points.min(axis=0) >= origin).all()
    assert (points.max(axis=0) <= last).all()

    assert grid["node_count"] == int(np.prod(counts))


def test_a_density_integrates_to_the_number_of_observations():
    from misah.kernels import covering_grid, density_field

    points = clustered_points()
    grid = covering_grid(points, [100.0, 100.0, 100.0], 1_200.0)

    values = density_field(
        points, grid["origin"], grid["spacing"], grid["counts"], [500.0] * 3
    )

    assert values.shape == (grid["node_count"],)

    cell_volume = float(np.prod(grid["spacing"]))
    total = values.sum() * cell_volume

    # The quadrature is a Riemann sum over a grid twenty times finer than the
    # bandwidth, so a part in a thousand is the accuracy on offer. A wrong
    # normalizing constant would be out by a factor, not by a part in a
    # thousand.
    assert abs(total - len(points)) < 0.5, total


def test_the_dimension_is_read_off_the_data():
    from misah.kernels import covering_grid, density_field

    rng = np.random.default_rng(11)

    # One dimension is not padding for a tidy range: a histogram of
    # hypocentral depths is a one-dimensional density, and this is how it is
    # estimated without choosing bin edges by hand.
    for dimension in (1, 2, 3):
        points = rng.normal(0.0, 1_000.0, size=(300, dimension))
        grid = covering_grid(points, [200.0] * dimension, 3_000.0)

        values = density_field(
            points,
            grid["origin"],
            grid["spacing"],
            grid["counts"],
            [800.0] * dimension,
        )

        cell = float(np.prod(grid["spacing"]))
        total = values.sum() * cell

        assert abs(total - len(points)) < 1.0, (dimension, total)

    # And a fourth is refused by name rather than silently truncated to three.
    try:
        points = rng.normal(size=(10, 4))
        density_field(points, [0.0] * 4, [1.0] * 4, [2] * 4, [1.0] * 4)
    except ValueError as exc:
        assert "not 4" in str(exc), exc
    else:
        raise AssertionError("a four-dimensional field was accepted")


def test_the_field_comes_back_in_the_grids_own_order():
    from misah.kernels import density_field

    # One observation, on a grid with three different counts, so a transposed
    # traversal cannot pass by coincidence. The peak has to land at the node
    # the observation sits on.
    origin, spacing, counts = [0.0, 0.0, 0.0], [100.0, 100.0, 100.0], [4, 5, 6]
    observation = np.array([[100.0, 200.0, 300.0]])

    values = density_field(observation, origin, spacing, counts, [400.0] * 3)

    # First axis fastest, which is what `order="F"` means to numpy and what a
    # VTK STRUCTURED_POINTS reader expects.
    volume = values.reshape(counts, order="F")
    assert volume.shape == (4, 5, 6)

    peak = np.unravel_index(np.argmax(volume), volume.shape)
    assert peak == (1, 2, 3), peak


def test_a_truncated_kernel_says_what_it_threw_away():
    from misah.kernels import kernel_profile

    quartic = kernel_profile([500.0] * 3)
    assert quartic["reach"] == [500.0, 500.0, 500.0]
    assert quartic["retained_mass"] == 1.0

    # An untruncated Gaussian never stops, and reporting no reach is the
    # warning that every node will cost the whole dataset.
    assert kernel_profile([500.0] * 3, "gaussian")["reach"] is None

    cut = kernel_profile([500.0] * 3, "truncated_gaussian", 4.0)
    assert cut["reach"] == [2_000.0] * 3
    # A thousandth outside four bandwidths in three dimensions, which is why
    # four is the cut worth reaching for.
    assert abs(cut["retained_mass"] - 0.998866) < 1e-5, cut

    # Truncation removes a tail; it does not rescale what is left, so the peak
    # is the plain Gaussian's.
    assert cut["peak"] == kernel_profile([500.0] * 3, "gaussian")["peak"]

    for bad in ([0.0, 1.0, 1.0], [1.0, -1.0, 1.0], [1.0, float("nan"), 1.0]):
        try:
            kernel_profile(bad)
        except ValueError:
            continue
        raise AssertionError(f"a bandwidth of {bad} was accepted")


def test_each_node_recovers_the_phase_that_is_near_it():
    from misah.kernels import stress_field

    positions, faults = two_phases()

    # Two nodes, one on each cluster, and a kernel too tight to see both.
    field, stats = stress_field(
        positions,
        faults,
        origin=[0.0, 0.0],
        spacing=[10_000.0, 1.0],
        counts=[2, 1],
        bandwidth=[2_000.0, 2_000.0],
        min_support=4.0,
    )

    assert stats["nodes"] == 2
    assert stats["nodes_inverted"] == 2
    assert field["has_solution"].all()

    # Each node saw only its own cluster.
    assert (field["faults_within_reach"] == 8).all(), field["faults_within_reach"]
    assert np.allclose(field["support"], 8.0)

    # S1 vertical in the west, horizontal in the east. Plunge is the reading
    # that does not depend on a trend that is undetermined when an axis is
    # vertical.
    assert field["s1"][0][1] > 80.0, field["s1"][0]
    assert field["s1"][1][1] < 10.0, field["s1"][1]

    # The data were generated without error and the truth sits on the search
    # grid, so the misfit is zero rather than small.
    assert (field["misfit_degrees"] < 1e-9).all(), field["misfit_degrees"]


def test_a_node_out_of_reach_is_a_hole_and_not_a_zero():
    from misah.kernels import stress_field

    positions, faults = two_phases()

    field, stats = stress_field(
        positions,
        faults,
        origin=[0.0, 0.0],
        spacing=[100_000.0, 1.0],
        counts=[3, 1],
        bandwidth=[2_000.0, 2_000.0],
        min_support=4.0,
    )

    assert stats["nodes_inverted"] == 1
    assert stats["nodes_without_data"] == 2

    empty = ~field["has_solution"]
    assert empty.sum() == 2

    # Every column from the solution is NaN where there is none. A zero would
    # plot: S1 at north and horizontal, a misfit of nothing, a shape ratio of
    # zero -- all of them readings, and none of them true.
    for key in ("phi", "misfit_degrees", "effective_sample_size",
                "runner_up_misfit_degrees", "runner_up_s1_degrees"):
        assert np.isnan(field[key][empty]).all(), key
    for key in ("s1", "s2", "s3"):
        assert np.isnan(field[key][empty]).all(), key

    # The density and the count are not from the solution and are real
    # numbers: zero, because nothing reached.
    assert (field["density"][empty] == 0.0).all()
    assert (field["faults_within_reach"][empty] == 0).all()


def test_thin_support_is_declined_rather_than_answered():
    from misah.kernels import stress_field

    positions, faults = two_phases()

    # One node on the western cluster, which eight faults reach.
    grid = dict(origin=[0.0, 0.0], spacing=[1.0, 1.0], counts=[1, 1],
                bandwidth=[2_000.0, 2_000.0])

    refused, stats = stress_field(positions, faults, min_support=9.0, **grid)
    assert stats["nodes_below_threshold"] == 1
    assert stats["nodes_inverted"] == 0
    assert not refused["has_solution"][0]

    # And it still says what it saw, which is what keeps a thin patch legible
    # rather than merely blank.
    assert refused["faults_within_reach"][0] == 8
    assert refused["density"][0] > 0.0

    accepted, stats = stress_field(positions, faults, min_support=4.0, **grid)
    assert stats["nodes_inverted"] == 1
    assert accepted["has_solution"][0]


def test_support_is_the_weighted_count_and_not_the_plain_one():
    from misah.kernels import stress_field

    # Two clusters at very different distances from one node. Every fault
    # reaches, so the count says sixteen; the support should say close to the
    # eight that are actually near.
    positions, faults = two_phases(separation=1_900.0)

    field, _ = stress_field(
        positions,
        faults,
        origin=[1_900.0, 0.0],
        spacing=[1.0, 1.0],
        counts=[1, 1],
        bandwidth=[2_000.0, 2_000.0],
        min_support=0.0,
    )

    assert field["faults_within_reach"][0] == 16
    assert field["support"][0] < 9.0, field["support"][0]


def test_the_cost_counts_the_nodes_the_field_goes_on_to_invert():
    from misah.kernels import field_cost, inversion_candidate_count, stress_field

    positions, faults = two_phases(separation=3_000.0)

    grid = dict(
        origin=[-2_000.0, 0.0],
        spacing=[1_000.0, 1.0],
        counts=[8, 1],
        bandwidth=[1_200.0, 1_200.0],
        min_support=4.0,
        angle_step_degrees=30.0,
        phi_step=0.5,
    )

    cost = field_cost(positions, faults, **grid)
    field, stats = stress_field(positions, faults, **grid)

    assert cost["nodes"] == 8
    assert cost["nodes_to_invert"] == stats["nodes_inverted"]

    # The work is the candidates times the faults in reach, node by node --
    # the sum the field itself performs.
    candidates = inversion_candidate_count(30.0, 0.5)
    expected = int(
        (candidates * field["faults_within_reach"][field["has_solution"]]).sum()
    )
    assert cost["forward_solutions"] == expected, cost

    # A run that would invert nothing is costed at nothing, which is the
    # answer worth having before starting rather than after.
    none = field_cost(positions, faults, **{**grid, "min_support": 1e9})
    assert none["nodes_to_invert"] == 0
    assert none["forward_solutions"] == 0


def test_the_density_column_is_the_density_field():
    from misah.kernels import density_field, stress_field

    positions, faults = two_phases(separation=3_000.0)

    grid = dict(origin=[-1_000.0, -1_000.0], spacing=[1_000.0, 1_000.0], counts=[6, 3])
    bandwidth = [2_500.0, 2_500.0]

    # A threshold nothing can clear, so this costs a neighbour search and no
    # inversions at all.
    field, _ = stress_field(
        positions, faults, bandwidth=bandwidth, min_support=1e9, **grid
    )
    separately = density_field(positions, bandwidth=bandwidth, **grid)

    # The same number, not a variant of it: both come from one kernel over one
    # neighbourhood, and a stress map whose density column had drifted from the
    # density field would be the exact failure the two-in-one-pass design
    # exists to prevent.
    assert np.array_equal(field["density"], separately)


def test_the_columns_stay_aligned_with_the_grid():
    from misah.kernels import stress_field

    positions, faults = two_phases(separation=4_000.0)

    counts = [5, 3]
    origin = [-1_000.0, -1_000.0]
    spacing = [1_500.0, 1_500.0]

    field, stats = stress_field(
        positions, faults, origin, spacing, counts, [3_000.0, 3_000.0],
        min_support=4.0,
    )

    node_count = counts[0] * counts[1]
    assert stats["nodes"] == node_count
    for key, values in field.items():
        assert len(values) == node_count, key

    # The positions are the grid's own nodes, first axis fastest.
    expected = np.array([
        [origin[0] + i * spacing[0], origin[1] + j * spacing[1]]
        for j in range(counts[1])
        for i in range(counts[0])
    ])
    assert np.array_equal(field["positions"], expected)

    # And the grid comes back in the stats, so a caller can reshape without
    # having kept it.
    assert stats["grid_origin"] == origin
    assert stats["grid_spacing"] == spacing
    assert stats["grid_counts"] == counts


def test_a_mismatch_between_places_and_faults_is_refused():
    from misah.kernels import stress_field

    positions, faults = two_phases()

    grid = dict(origin=[0.0, 0.0], spacing=[1.0, 1.0], counts=[1, 1],
                bandwidth=[2_000.0, 2_000.0])

    # Refused rather than answered: pairing faults with the wrong places would
    # return a field that looks like a result.
    try:
        stress_field(positions[:-1], faults, **grid)
    except ValueError as exc:
        assert "rows" in str(exc), exc
    else:
        raise AssertionError("positions of the wrong length were accepted")

    # A NaN threshold would leave every node above it -- `<` being false
    # against a NaN -- and invert the whole grid on whatever a single fault
    # said. Silently, and expensively.
    try:
        stress_field(positions, faults, min_support=float("nan"), **grid)
    except ValueError:
        pass
    else:
        raise AssertionError("a NaN min_support was accepted")

    for bad in ({"spacing": [0.0, 1.0]}, {"counts": [0, 1]}, {"origin": [float("nan"), 0.0]}):
        try:
            stress_field(positions, faults, **{**grid, **bad})
        except ValueError:
            continue
        raise AssertionError(f"a grid with {bad} was accepted")


def test_an_unknown_kernel_is_named_rather_than_defaulted():
    from misah.kernels import density_field

    points = clustered_points(count=20)

    try:
        density_field(points, [0.0] * 3, [1.0] * 3, [2] * 3, [1.0] * 3, kernel="epanechnikov")
    except ValueError as exc:
        assert "epanechnikov" in str(exc), exc
    else:
        raise AssertionError("an unknown kernel silently became something else")


def test_package_exposes_the_field_functions():
    import misah

    for name in ("covering_grid", "kernel_profile", "density_field",
                 "stress_field", "field_cost"):
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
