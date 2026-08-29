"""Check the best-fit plane bindings by closing the loop on the mesh kernel.

The pair `intersect_mesh_grid` and `best_fit_planes` are the two halves of one
operation: cut a surface against a DEM, then read the surface's attitude back
out of where the cut fell. Nothing here is checked against stored numbers --
the DEM and the surface are both built from a plane of known attitude, and the
test is that the attitude survives the round trip.

Numbers in, numbers out, so this runs in CI. The comparison against
geoSurfDEM's own `BestFitGeoplanes` output for Timpa San Lorenzo lives outside
this tree, like the Malpi golden trace, and is run by hand through
`cargo run --example best_fit_csv`.

Run from anywhere but the `pylib` directory, which shadows the installed
package with its source tree:

    python3 pylib/tests/test_best_fit.py
"""

from __future__ import annotations

import sys

import numpy as np

CELL = 1.0
NODATA = -9999.0


def plane_z(azimuth, dip, x, y, through=(0.0, 0.0, 0.0)):
    """Elevation of a plane of the given attitude at the given map positions."""
    from misah.kernels import plane_normal

    n = np.asarray(plane_normal(azimuth, dip))
    x0, y0, z0 = through
    return z0 - (n[0] * (x - x0) + n[1] * (y - y0)) / n[2]


def dem_of(azimuth, dip, relief=10.0, wavelength=15.0, rows=60, columns=60):
    """A tilted plane with sinusoidal relief on it, and its geotransform.

    `relief` is what makes this method work at all, and setting it to zero is
    the test of the opposite case. A plane cut by a *smooth* slope meets it
    along a straight line, and no straight line determines a plane; a plane cut
    by rough ground meets it along a trace that runs up a spur and back down a
    gully, covering a genuinely two-dimensional patch of each cell. That is why
    an attitude can be recovered from a mapped contact in real terrain, and why
    it cannot on a planar hillside.
    """
    col, row = np.meshgrid(np.arange(columns), np.arange(rows))
    x = col * CELL
    y = -row * CELL

    z = plane_z(azimuth, dip, x, y)
    if relief:
        z = z + relief * np.sin(2 * np.pi * x / wavelength) * np.sin(2 * np.pi * y / wavelength)

    return np.ascontiguousarray(z), (0.0, CELL, 0.0, 0.0, 0.0, -CELL)


def surface_points(azimuth, dip, dem, geotransform):
    """Where a plane of the given attitude cuts the DEM, through the mesh kernel."""
    from misah.kernels import intersect_mesh_grid

    # A large triangulated square lying in the plane: two triangles are enough
    # to cover the DEM, and the mesh kernel takes vertices and faces.
    span = 200.0
    corners = [(-span, -span), (span, -span), (span, span), (-span, span)]
    vertices = np.array(
        [[cx, cy, plane_z(azimuth, dip, cx, cy, through=(30.0, -30.0, 0.0))] for cx, cy in corners],
        dtype=np.float64,
    )
    faces = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int64)

    points, _, _, _ = intersect_mesh_grid(dem, geotransform, vertices, faces, NODATA)
    return points


def test_the_mesh_kernel_output_feeds_straight_in():
    """The whole point: points out of one kernel, attitudes out of the next."""
    from misah.kernels import best_fit_planes

    dem, geotransform = dem_of(90.0, 20.0)
    points = surface_points(135.0, 35.0, dem, geotransform)
    assert len(points) > 100, f"the surface barely cut the DEM: {len(points)} points"

    field, stats = best_fit_planes(points, cell_size=10.0)

    assert stats["cells_fitted"] > 10
    assert stats["cells_fitted"] == len(field["attitudes"])
    assert stats["cells_collinear"] == 0, "rough ground should leave no cell undetermined"

    # Every cell must return the attitude of the surface that was cut.
    assert np.allclose(field["attitudes"][:, 0], 135.0, atol=1e-3)
    assert np.allclose(field["attitudes"][:, 1], 35.0, atol=1e-3)


def test_a_plane_cut_by_a_smooth_slope_yields_nothing():
    """Two planes meet along one straight line, and no line determines a plane.
    Every cell along it is collinear and the field comes back empty -- which is
    the honest answer, and the one the reference does not give."""
    from misah.kernels import best_fit_planes

    dem, geotransform = dem_of(90.0, 20.0, relief=0.0)
    points = surface_points(135.0, 35.0, dem, geotransform)

    field, stats = best_fit_planes(points, cell_size=10.0)

    assert stats["cells_fitted"] == 0, "a straight line was fitted with a plane"
    assert stats["cells_collinear"] > 0
    assert len(field["attitudes"]) == 0
    assert (
        stats["cells_fitted"] + stats["cells_collinear"] + stats["cells_too_few_points"]
        == stats["non_empty_cells"]
    )


def test_a_plane_is_recovered_from_scattered_points_on_it():
    from misah.kernels import best_fit_planes

    rng = np.random.default_rng(7)
    for azimuth, dip in ((0.0, 30.0), (135.0, 35.0), (250.0, 65.0), (40.0, 8.0)):
        x = rng.uniform(0.0, 100.0, 600)
        y = rng.uniform(0.0, 100.0, 600)
        points = np.column_stack([x, y, plane_z(azimuth, dip, x, y)])

        field, _ = best_fit_planes(points, cell_size=25.0)

        assert len(field["attitudes"]) > 4, f"{azimuth}/{dip}: too few cells"
        assert np.allclose(field["attitudes"][:, 0], azimuth, atol=1e-3)
        assert np.allclose(field["attitudes"][:, 1], dip, atol=1e-3)


def test_every_array_in_the_field_has_the_same_length():
    from misah.kernels import best_fit_planes

    rng = np.random.default_rng(3)
    x, y = rng.uniform(0.0, 100.0, 400), rng.uniform(0.0, 100.0, 400)
    points = np.column_stack([x, y, plane_z(45.0, 50.0, x, y)])

    field, stats = best_fit_planes(points, cell_size=20.0)
    n = stats["cells_fitted"]

    assert field["cell_centres"].shape == (n, 2)
    assert field["attitudes"].shape == (n, 2)
    assert field["centroids"].shape == (n, 3)
    assert field["singular_values"].shape == (n, 3)
    assert field["cells"].shape == (n, 2)
    assert field["point_counts"].shape == (n,)
    assert field["collinearity"].shape == (n,)
    assert field["rms_distance"].shape == (n,)


def test_the_singular_values_descend_and_the_residual_follows_the_smallest():
    from misah.kernels import best_fit_planes

    rng = np.random.default_rng(11)
    x, y = rng.uniform(0.0, 100.0, 400), rng.uniform(0.0, 100.0, 400)
    points = np.column_stack([x, y, plane_z(45.0, 50.0, x, y)])

    field, _ = best_fit_planes(points, cell_size=20.0)
    s = field["singular_values"]

    assert np.all(s[:, 0] >= s[:, 1]) and np.all(s[:, 1] >= s[:, 2])
    # rms = s3 / sqrt(n), which is the definition and worth pinning.
    assert np.allclose(field["rms_distance"], s[:, 2] / np.sqrt(field["point_counts"]))
    # Points exactly on a plane: the residual is the arithmetic floor, not data.
    assert np.all(field["rms_distance"] < 1e-6)


def test_the_cells_come_back_in_grid_order():
    from misah.kernels import best_fit_planes

    rng = np.random.default_rng(5)
    x, y = rng.uniform(0.0, 200.0, 900), rng.uniform(0.0, 200.0, 900)
    points = np.column_stack([x, y, plane_z(200.0, 25.0, x, y)])

    field, _ = best_fit_planes(points, cell_size=40.0)
    cells = field["cells"]

    order = [tuple(c) for c in cells]
    assert order == sorted(order), "a hash map's order leaked through the boundary"


def test_the_collinearity_threshold_is_what_sets_cells_aside():
    from misah.kernels import best_fit_planes

    dem, geotransform = dem_of(90.0, 20.0)
    points = surface_points(135.0, 35.0, dem, geotransform)

    strict, strict_stats = best_fit_planes(points, cell_size=10.0, max_collinearity=0.5)
    loose, loose_stats = best_fit_planes(points, cell_size=10.0, max_collinearity=50.0)

    assert loose_stats["cells_fitted"] >= strict_stats["cells_fitted"]
    assert strict_stats["cells_collinear"] >= loose_stats["cells_collinear"]
    # Whatever survives the strict threshold is below it.
    assert np.all(strict["collinearity"] <= 0.5)
    assert len(loose["attitudes"]) == loose_stats["cells_fitted"]


def test_coincident_points_are_merged():
    from misah.kernels import best_fit_planes

    rng = np.random.default_rng(13)
    x, y = rng.uniform(0.0, 50.0, 200), rng.uniform(0.0, 50.0, 200)
    points = np.column_stack([x, y, plane_z(90.0, 45.0, x, y)])
    doubled = np.vstack([points, points])

    _, stats = best_fit_planes(doubled, cell_size=25.0)

    assert stats["input_points"] == 400
    assert stats["distinct_points"] == 200


def test_the_grid_is_reported_with_the_field():
    from misah.kernels import best_fit_planes

    points = np.array(
        [[100.0, 500.0, 0.0], [200.0, 400.0, 10.0], [150.0, 450.0, 5.0], [120.0, 480.0, 2.0]]
    )

    field, stats = best_fit_planes(points, cell_size=50.0)

    assert stats["cell_size"] == 50.0
    assert stats["grid_x_min"] == 100.0
    assert stats["grid_y_max"] == 500.0
    assert stats["grid_rows"] == 3 and stats["grid_columns"] == 3

    # Cell centres follow from the grid, half a cell in from the top-left.
    for (row, col), (cx, cy) in zip(field["cells"], field["cell_centres"]):
        assert np.isclose(cx, 100.0 + (col + 0.5) * 50.0)
        assert np.isclose(cy, 500.0 - (row + 0.5) * 50.0)


def test_too_few_points_give_no_field():
    from misah.kernels import best_fit_planes

    assert best_fit_planes(np.zeros((0, 3)), cell_size=10.0) is None

    # Points but no cell with three of them: a field with nothing in it.
    field, stats = best_fit_planes(np.array([[0.0, 0.0, 0.0], [100.0, 100.0, 5.0]]), cell_size=1.0)
    assert len(field["attitudes"]) == 0
    assert stats["cells_too_few_points"] == 2


def test_the_inputs_are_checked():
    from misah.kernels import best_fit_planes

    points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])

    try:
        best_fit_planes(np.zeros((4, 2)), cell_size=10.0)
    except ValueError:
        pass
    else:
        raise AssertionError("an (N, 2) array was accepted")

    for kwargs in (dict(cell_size=0.0), dict(cell_size=-5.0), dict(cell_size=float("nan")),
                   dict(cell_size=10.0, coincidence_distance=0.0)):
        try:
            best_fit_planes(points, **kwargs)
        except ValueError:
            continue
        raise AssertionError(f"{kwargs} was accepted")


def test_package_exposes_the_best_fit_function():
    import misah

    assert hasattr(misah.kernels, "best_fit_planes")


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
