# misah

*misah* is an experiment in creating a Rust-based geospatial library with Python bindings (thanks to the *pyo3* library).


Currently it is in alpha mode.

## Plane-grid intersection

`lib/src/raster/intersection.rs` intersects a plane with a raster grid, as the
zero-level set of `Plane::signed_distance_to_point` sampled at the grid nodes and
extracted by marching squares. Building on the signed distance instead of a
`z = z(x, y)` plane expression keeps vertical planes an ordinary case, and the
cost is one pass over the cells — geoSurfDEM's `IntersectDEM` solves the more
general surface-mesh problem by testing every DEM triangle against every mesh
triangle.

`GeologicalPlane::to_plane` bridges an attitude to the geometric plane, so a
caller works in dip direction and dip angle rather than in normals.

```sh
cargo test -p misah
cargo run --release -p misah --example plane_dem_asc -- \
    <dem.asc> <x0> <y0> <z0> <dip_dir> <dip_angle>
```

The kernel was checked against the geoSurfDEM golden dataset (plane 135/35 over
the Malpi ASTER DEM): 404 vertices, all on the plane to 1e-9, coincident with the
C++ trace to well within a cell. Two properties of that reference output are
worth recording, since they bound how closely anything can agree with it: about
half of its 1228 points are duplicates (618 unique at centimetre precision),
because shared triangle edges are emitted once per triangle, and four points —
CSV rows 644-647, in two pairs sharing an elevation — lie up to 3.56 m off the
plane. A plane fitted to the whole set returns 134.9999/35.0005 and passes 8 mm
from the nominal source point, so those four are a defect of that run rather than
a disagreement over conventions.

### Status

`raster::io` does not compile and is not wired into the module tree, so the
example carries its own ESRI ASCII reader. `pylib` is still on pyo3 0.13, which
predates Python 3.11: it builds, but no extension made from it will import into a
current interpreter, and the kernel is not exposed there yet.
