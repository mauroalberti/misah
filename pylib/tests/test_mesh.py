"""Check the mesh-grid binding: a surface of arbitrary shape cut against a DEM.

Data-free, like `test_stress.py`, so the CI can run it: the grids and meshes
here are built in a few lines rather than read from the geoSurfDEM crop that
`test_kernels.py` needs. The kernel itself is checked against that crop by the
Rust tests and by `lib/examples/mesh_dem_vtk.rs`; what this file is for is the
binding -- array shapes, index handling, and the statistics dict.

Run from anywhere but the `pylib` directory, which shadows the installed
package with its source tree:

    python3 pylib/tests/test_mesh.py
"""

from __future__ import annotations

import sys

import numpy as np

# A north-up grid of unit cells with its lower-left corner at the origin, so a
# node (row, col) sits at x = col + 0.5, y = nrows - row - 0.5.
NROWS, NCOLS = 10, 10
GEOTRANSFORM = (0.0, 1.0, 0.0, float(NROWS), 0.0, -1.0)


def flat_grid(z=0.0):
    return np.ascontiguousarray(np.full((NROWS, NCOLS), float(z)))


def vertical_wall(x):
    """A big triangle standing vertically at the given easting, facing east."""
    vertices = np.array(
        [[x, -100.0, -100.0], [x, 100.0, -100.0], [x, 0.0, 100.0]], dtype=float
    )
    faces = np.array([[0, 1, 2]], dtype=np.int64)
    return vertices, faces


def test_a_wall_cuts_a_flat_grid_along_its_own_plane():
    from misah.kernels import intersect_mesh_grid

    vertices, faces = vertical_wall(4.7)
    points, attitudes, triangles, stats = intersect_mesh_grid(
        flat_grid(), GEOTRANSFORM, vertices, faces
    )

    assert len(points) > 0
    assert points.shape[1] == 3
    assert np.abs(points[:, 0] - 4.7).max() < 1e-9, "every point should lie on the wall"

    # One triangle in the mesh, so every point traces back to index 0.
    assert triangles.shape == (len(points),)
    assert set(triangles.tolist()) == {0}

    # A vertical plane dips 90; its dip direction is the pair the winding
    # settles, so only the angle is asserted here.
    assert attitudes.shape == (len(points), 2)
    assert np.abs(attitudes[:, 1] - 90.0).max() < 1e-9


def test_a_dipping_surface_reports_its_own_attitude():
    from misah.kernels import intersect_mesh_grid, plane_normal

    # A patch dipping 45 towards 090, built large enough to cross the grid.
    big = 100.0
    centre = np.array([5.0, 5.0, 0.0])
    normal = np.array(plane_normal(90.0, 45.0))
    strike = np.array([0.0, 1.0, 0.0])           # horizontal, along strike
    down_dip = np.cross(normal, strike)

    vertices = np.array(
        [centre + big * (np.cos(a) * strike + np.sin(a) * down_dip)
         for a in np.radians([90.0, 210.0, 330.0])],
        dtype=float,
    )
    faces = np.array([[0, 1, 2]], dtype=np.int64)

    points, attitudes, _, _ = intersect_mesh_grid(
        flat_grid(), GEOTRANSFORM, vertices, faces
    )

    assert len(points) > 0
    assert np.abs(attitudes[:, 0] - 90.0).max() < 1e-6, "dip direction"
    assert np.abs(attitudes[:, 1] - 45.0).max() < 1e-6, "dip angle"


def test_a_mesh_beside_the_grid_finds_nothing():
    from misah.kernels import intersect_mesh_grid

    vertices, faces = vertical_wall(5000.0)
    points, attitudes, triangles, stats = intersect_mesh_grid(
        flat_grid(), GEOTRANSFORM, vertices, faces
    )

    assert points.shape == (0, 3)
    assert attitudes.shape == (0, 2)
    assert triangles.shape == (0,)
    # The bounding-box test settles it without walking a single DEM triangle.
    assert stats["mesh_triangles_outside_grid"] == 1
    assert stats["dem_triangle_pairs"] == 0


def test_the_statistics_describe_the_run():
    from misah.kernels import intersect_mesh_grid

    vertices, faces = vertical_wall(4.7)
    _, _, _, stats = intersect_mesh_grid(flat_grid(), GEOTRANSFORM, vertices, faces)

    for key in (
        "mesh_triangles",
        "degenerate_mesh_triangles",
        "mesh_triangles_outside_grid",
        "dem_triangle_pairs",
        "coplanar_sides",
        "duplicate_crossings",
    ):
        assert key in stats, key
        assert isinstance(stats[key], int), key

    assert stats["mesh_triangles"] == 1
    assert stats["degenerate_mesh_triangles"] == 0
    # Sides shared between adjacent DEM triangles are cut once, and the repeats
    # suppressed are counted rather than hidden.
    assert stats["duplicate_crossings"] > 0


def test_nodata_removes_the_cells_that_touch_it():
    from misah.kernels import intersect_mesh_grid

    vertices, faces = vertical_wall(4.7)

    dem = flat_grid()
    intact, _, _, _ = intersect_mesh_grid(dem, GEOTRANSFORM, vertices, faces, -9999.0)

    holed = flat_grid()
    holed[5, 4] = -9999.0
    with_hole, _, _, _ = intersect_mesh_grid(holed, GEOTRANSFORM, vertices, faces, -9999.0)

    assert len(with_hole) < len(intact), "the hole changed nothing"


def test_a_non_contiguous_dem_is_refused():
    from misah.kernels import intersect_mesh_grid

    vertices, faces = vertical_wall(4.7)
    strided = np.zeros((NROWS, NCOLS * 2))[:, ::2]

    try:
        intersect_mesh_grid(strided, GEOTRANSFORM, vertices, faces)
    except ValueError:
        return
    raise AssertionError("non-contiguous DEM was accepted")


def test_a_face_index_past_the_vertex_pool_is_refused():
    from misah.kernels import intersect_mesh_grid

    vertices, _ = vertical_wall(4.7)
    faces = np.array([[0, 1, 9]], dtype=np.int64)

    try:
        intersect_mesh_grid(flat_grid(), GEOTRANSFORM, vertices, faces)
    except ValueError:
        return
    raise AssertionError("an out-of-range face index was accepted")


def test_a_negative_face_index_is_refused():
    """numpy allows a negative index and a usize cast would turn it into an
    enormous positive one, so it has to be caught rather than converted."""
    from misah.kernels import intersect_mesh_grid

    vertices, _ = vertical_wall(4.7)
    faces = np.array([[0, 1, -1]], dtype=np.int64)

    try:
        intersect_mesh_grid(flat_grid(), GEOTRANSFORM, vertices, faces)
    except ValueError:
        return
    raise AssertionError("a negative face index was accepted")


def test_arrays_of_the_wrong_width_are_refused():
    from misah.kernels import intersect_mesh_grid

    good_vertices, good_faces = vertical_wall(4.7)

    for vertices, faces in (
        (np.zeros((3, 2)), good_faces),
        (good_vertices, np.zeros((1, 4), dtype=np.int64)),
    ):
        try:
            intersect_mesh_grid(flat_grid(), GEOTRANSFORM, vertices, faces)
        except ValueError:
            continue
        raise AssertionError("an array of the wrong width was accepted")


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
