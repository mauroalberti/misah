# misah

*misah* is an experiment in creating a Rust-based geospatial library with Python bindings (thanks to the *pyo3* library).


Currently it is in alpha mode.

## Plane-grid intersection

`lib/src/raster/intersection.rs` intersects a plane with a raster grid, as the
zero-level set of `Plane::signed_distance_to_point` sampled at the grid nodes and
extracted by marching squares. Building on the signed distance instead of a
`z = z(x, y)` plane expression keeps vertical planes an ordinary case, and the
cost is one pass over the cells. The general case, a surface of arbitrary shape,
is `raster::mesh_intersection` below.

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

## Mesh-grid intersection

`lib/src/raster/mesh_intersection.rs` cuts the grid with a triangulated surface
rather than a single plane — the case a folded or faulted geological surface
needs. It is a port of geoSurfDEM's `IntersectDEM`, the one part of that suite
with no equivalent in geogst. Each point comes out with the attitude of the mesh
triangle that produced it, so the result is a set of located attitudes and not
merely a trace. What consumes them is **Best-fit geoplanes** below.

```sh
cargo run --release -p misah --example mesh_dem_vtk -- <dem.asc> <surface.vtk>
```

`geometry::triangle::Triangle3D` and `geometry::mesh::TriangleMesh` carry the
surface; the mesh is indexed (a vertex pool plus index triplets) because that is
how a VTK `POLYDATA` file gives it. The DEM is read by `raster::io`; reading VTK
is left to the example, being no business of a raster module.

Three things depart from the C++ deliberately:

- **Where the points come from.** The original intersected the two triangle
  planes into a line, then intersected that line with each DEM side. Here a DEM
  side is cut against the mesh triangle's plane directly, interpolating between
  the signed distances of its two endpoints — the construction the plane kernel
  already uses. Same points, without the intermediate line, whose construction
  needed an arbitrary 100-unit displacement, an unguarded division, and could
  return an uninitialised point when all three of its determinants fell below
  tolerance.
- **What gets compared.** A DEM *is* a uniform spatial index, so each mesh
  triangle is mapped through the inverse geotransform (`GeoTransform::node_rc`,
  new) onto the cells its footprint covers, widened by one cell so that a
  surface passing exactly through a column of nodes still finds its crossings.
  No index is built, because the grid already is one.
- **Shared sides.** Crossings are keyed on the DEM side that produced them, so
  the side two adjacent DEM triangles have in common is cut once. This is the
  source of the duplication in the reference output noted above; the count of
  suppressed repeats is reported in the run's statistics rather than hidden.

Against the same golden dataset (surface 135/35 over the Malpi ASTER DEM), the
port returns 614 points against the reference's 618 distinct ones, every one
reported at 135.00/35.00, in 29 ms. Two independent checks, both at the 1 cm the
CSV is written to:

| | on the plane, max residual | off the DEM surface, median / max |
| --- | --- | --- |
| misah | 0.008 m | 0.003 m / 0.010 m |
| C++ reference | 3.55 m | 0.038 m / 0.213 m |

Both figures for misah are the CSV rounding itself, so the points sit on the
plane and on the topography as exactly as the output can express; the reference
misses the DEM surface by more than a centimetre at 511 of its 617 testable
points, which is where the ~0.18 m median offset between the two point sets
comes from. The spatial restriction takes the comparisons from the 49 million
pairs the C++ makes after its own volume prefilter (505 mesh triangles × 97,908
DEM triangles) down to 23,606.

Note that `test_data/src_planes/malpi_045_90.vtk` in geoSurfDEM is a byte-copy
of `malpi_135_35.vtk`, so that suite has no vertical-surface case despite the
name; the vertical and node-aligned configurations are covered by unit tests
instead.

## Best-fit geoplanes

`lib/src/structural/best_fit.rs` reads an attitude back out of located points,
by fitting a plane to each cell of a regular grid over them. It is the consumer
the mesh kernel was written for: that one says where a surface meets the
topography, this one says what the surface's orientation was. The two are not
coupled — the input here is bare coordinates, so a contact digitised from a map
or a set of readings along an outcrop feeds it just as well.

```sh
cargo run --release -p misah --example best_fit_csv -- \
    <points.csv> <cell_size> [coincidence_distance] [max_collinearity]
```

A port of geoSurfDEM's `BestFitGeoplanes`, whose numerical core is the Fortran
`invert_main_attitude` in `GeoInversions`: centre a cell's points, take the SVD
of the m-by-3 matrix through LAPACK, read the normal off the right singular
vector of the smallest singular value. Unlike `ForwardStress.f95` above, this
Fortran compiles clean under `-fimplicit-none` — checked before porting, since
that is exactly the bug that moved two of the five stress cases — so its numbers
can be trusted as they stand.

`algebra::eigen` is new and is why LAPACK does not come with it: the right
singular vectors of the centred coordinates are the eigenvectors of their 3×3
scatter matrix, which a page of Jacobi rotations decomposes exactly. The cost is
real and worth naming — forming `XᵀX` before decomposing it squares the
condition number, so about half the available significant digits go — and it is
paid in a place where nothing needs them: an attitude is wanted to a hundredth
of a degree, and the residual floor this imposes is around 1e-7 of the coordinate
scale. Jacobi rather than the analytic solution of the characteristic cubic,
which is faster and loses accuracy precisely on the nearly degenerate matrices
this sees most.

### What a straight trace cannot tell you

The substantive departure from the reference. A plane is determined by a cloud
of points only if the cloud spreads in **two** directions; the test is `s2/s1`,
the second singular value against the first. Where the points fall along a line
every plane through that line fits equally well, and the one returned is
whichever direction rounding happened to favour — arbitrary, not imprecise, and
with a residual that looks perfect either way.

This is not an edge case. It is the ordinary condition for a contact crossing a
smooth slope: two planes meet along one straight line. The method works on real
terrain because rough ground makes the trace wander — up a spur, back down a
gully — so a cell catches a genuinely two-dimensional patch of it. On a planar
hillside it cannot work at all, and the honest output is nothing.

The reference computes `s2/s1` as `logratio_s(1)`, writes it to its results file
alongside two other ratios, and publishes an attitude regardless of any of them.
Here it is `BestFitPlane::collinearity`, and cells above `max_collinearity` are
counted in the statistics instead of entering the field.

**Which ratio to test on is not obvious, and the plausible choice is wrong.**
`-log10(s3/s2)` reads like the measure of a good plane and has, on the reference
dataset, no power to separate at all: the two worst cells score *better* on it
than most cells correct to a thousandth of a degree. When `s2` is itself
rounding noise, a ratio taken against it is a ratio between two noises. That
mistake was made here first and caught by measuring, which is the only reason
this paragraph exists.

### Against the reference

geoSurfDEM's own Timpa San Lorenzo dataset — 1900 points, from intersecting two
known planes (67.298/39 and 77.292/40) with a DEM, at 50 m cells:

| | cells | error vs. the source plane, median | worst |
| --- | --- | --- | --- |
| misah | 288 | 0.0006° | **0.089°** |
| C++ reference | 311 | 0.0007° | **38.4°** |

On the 288 cells both accept, dip direction and dip angle agree to every one of
the six decimals either writes — the Jacobi-on-scatter-matrix route reproduces
LAPACK exactly where the fit is conditioned at all — and the point count per cell
is identical, so the deduplication and cell assignment match too. The 23 misah
sets aside are the difference. Three of them the reference publishes with
attitudes 13.2°, 12.5° and 0.7° away from the plane that generated their own
points; on the last, `s2` and `s3` are both zero, so its diagnostic column holds
`-log10(0/0)` — a NaN printed beside a confident dip direction. The run takes
6 ms.

`DEFAULT_MAX_COLLINEARITY` is 3.0 for measured reasons rather than tasteful
ones: on that dataset the sound cells run up to 3.3 and the wrong ones start at
3.1, so the populations overlap and no threshold separates them cleanly. At 3.0
the worst surviving attitude is 0.089° out; at 4.0, 0.873°; with no threshold,
38.4°. Tenths of a degree are far below what a compass reading is worth, so this
buys a margin nobody will measure, and buying more would discard sound cells to
chase a precision the data does not have.

Two smaller departures. Duplicate points are found by binning at the coincidence
distance rather than by scanning every kept point, which is what the reference
does — quadratic, invisible at 1900 points, not invisible at what the mesh kernel
can produce. And coincidence is judged in three dimensions rather than two: for
points off a DEM the two agree, since a single-valued surface cannot put two
elevations at one map location, but for an overturned contact or a set of
readings up a cliff the 2D test discards exactly the measurements that constrain
a steep plane.

## Forward stress on a fault plane

`lib/src/structural/stress.rs` is the direct Wallace-Bott problem: given a
reduced stress tensor (Angelier, 1984 — S1/S3 orientation and the shape ratio
Phi), predict the slip a fault plane would show under it. It is the inverse of
what a fault-slip inversion solves for, and it is a port of `ForwardStress.f95`,
a Fortran tool of Alberti (2010) that once lived in this repository's own
history under `GeoFaults`, before being dropped for a Rust port that was never
carried through.

This one arrived by way of geogst first (`geogst.core.geology.stress`), where
`Plane.rake_to_direct` turned out to already implement the Fortran's own
rake-to-slickenline formula exactly — both are Aki & Richards (1980) — and
where the whole algorithm was checked against the Fortran source itself:
compiled unmodified with gfortran, driven by a small program calling
`stresssolution_calc` directly, over five geological cases. Three matched
outright. Two did not, and the reason was a real bug in the 2010 source:
`module vector_processing` carries no `implicit none`, so `vector1_magn` inside
`vector_normalization` — reached through `vector_projection`, used for the
normal-stress vector every solution goes through — is an undeclared name,
silently typed as single-precision `REAL` by Fortran's implicit-typing rule
rather than the double precision the rest of the tool uses. The truncation is
usually too small to matter, but the predicted rake reaches it through an
`acos` close to the edge of its domain in near-strike-slip cases, where an
error below float32 precision in the cosine becomes one four orders of
magnitude larger in the angle — 1.4e-5 degrees on one case, arrived at through
a chain of double-precision arithmetic that looks correct at a glance.
Recompiling with `-fimplicit-none`, declaring `vector1_magn`'s type, moves
those two cases onto what the geogst port already computed, agreeing to 12+
significant digits; an independent 50-digit `mpmath` recomputation of the same
five cases confirms the corrected numbers, not the original binary's, are the
mathematically correct ones.

This Rust port is checked against the same five cases and the same corrected
values — `structural::stress::tests::against_the_fortran_original` — so the
same bug does not need rediscovering here.

`ReducedStressTensor::new` takes S1 and S3 as `GeologicalAxis`, sub-orthogonal,
plus Phi; S2 completes a right-handed cyclic triad (`S1 x S2 = S3`,
`S2 x S3 = S1`), matching the sign the Fortran original builds. `sigma1`/`sigma3`
default to 1/0 through `ReducedStressTensor::normalized` — the usual
normalization when only the tensor's shape is known, as from an inversion —
which leaves the predicted rake correct while making slip tendency and
deformation index meaningless; supply the true magnitudes when they are known
and those are wanted. `solve` resolves the tensor onto a `GeologicalPlane` and
returns a `StressSolution`: traction, normal and shear stress always, and,
where the shear is not numerical noise, the predicted rake, slickenline, slip
tendency and deformation index. `angular_misfit` compares the prediction
against an observed slickenline — what an inversion would minimize, run
forward.

Two small primitives came out of this that the rest of the structural code
will reuse: `structural::geol_axis::GeologicalAxis::as_versor`/`from_versor`
(trend/plunge in `(East, North, Up)` and back), and
`GeologicalPlane::rhr_strike`/`from_rhr_strike`/`rake_to_versor`, the last
being the Aki & Richards formula itself, exposed rather than kept private to
`solve`, since it is the piece any future fault-slip *inversion* would need to
call once per candidate tensor per fault — which is also the reason this went
into `misah` and not only geogst: the direct problem is O(1) per fault, with no
performance case of its own, but it is the exact inner loop a grid search over
candidate tensors would evaluate, at whatever count of faults and grid points
made that search worth compiling.

`solve_stress` and `rake_to_slickenline`, in `pylib`, are the Python surface of
this — one call each, taking and returning plain numbers rather than exposing
`ReducedStressTensor`/`GeologicalAxis`/`GeologicalPlane` themselves, matching
how `intersect_plane_grid` already builds and consumes a `GeologicalPlane`
without handing the type to the caller. Both are checked from Python against
the same five Fortran-verified cases.

## Fault-slip inversion

`lib/src/structural/inversion.rs` runs the forward model backwards: given
faults with their observed slip, it searches for the reduced stress tensor that
best explains them. That search is what made the forward problem worth
compiling, and what the sentence above was promising — one forward solution per
candidate tensor per fault, which for the default grid of 71 280 candidates
over a hundred faults is seven million of them — 1.5 s.

`FaultPlane` is the input, and is what connects a group of types that until now
led nowhere: `Direction3D` was reached only by `Slickenline`, `Slickenline` only
by `FaultPlane`, and `FaultPlane` by nothing at all. It now validates the
invariant its absence let through — that a slickenline must lie *in* the plane
it is recorded on, within a degree, since a lineation off the plane is not a
fault with a small error but two measurements that do not belong together.

Two details of the domain that the geometry alone would get wrong:

- **An unread slip sense is a line, not a vector.** Where the sense of movement
  was not determined, a prediction pointing the other way along the same
  lineation fits perfectly, so `Slickenline::angle_to` takes the misfit modulo
  180 degrees. Scoring such a fault at 180 instead of 0 would penalise an
  inversion for exactly the faults nobody could read the sense of.
- **A fault the model cannot speak about is not a fault it fits badly.** A plane
  lying on a principal stress axis carries no shear, so no slip direction is
  predicted; `misfit_on` returns `None` there rather than zero, and
  `mean_misfit` reports how many faults its average actually came from, so a
  candidate that silences half the dataset cannot win on the strength of the
  remainder.

The search is exhaustive rather than a descent, because the misfit surface is
not convex — a fault set carrying two superposed tectonic phases has two minima
by construction, and a descent would report whichever it fell into without ever
saying the other was there. `runners_up` is returned alongside `best` for the
same reason: a minimum that stands alone reads differently from one on a
plateau.

The test that matters generates faults from a known tensor through the forward
model and asks the inversion to recover it, which it does to within the grid
step for a vertical-S1 and for a strike-slip setting alike.

`invert_stress`, in `pylib`, is the Python surface of this, and it is where the
question the other bindings did not have to answer came up: a caller handing
over a hundred faults wants to pass arrays, not build a hundred objects across
the pyo3 boundary. So the faults arrive as one `(N, 4)` array — strike, dip,
and the trend and plunge of the slickenline — and the result comes back as a
dict of dicts rather than as `InversionResult`, `ScoredTensor` and
`ReducedStressTensor` exposed as types. Same rule as everywhere else on this
surface: plain arrays and numbers in, arrays and dicts out.

The slip sense is the one thing that could not be a column of that array
without lying about itself. `SlipSense` has four variants and the search reads
none of them — only whether the sense is known at all, the trend and plunge
already saying which way the hanging wall moved. So the boundary carries a
boolean per fault, `senses`, and not a code: `False` where nobody could read
the sense, which is what puts that fault's misfit modulo 180 degrees. Omitting
the array claims every sense was read. That is the right default for data taken
as given, and the wrong one for a dataset where the sense went unrecorded —
those faults would be scored at 180 degrees for fitting perfectly the other way
round, which is the whole failure `Slickenline::angle_to` exists to prevent.
`test_an_unread_sense_is_not_penalised` asserts both halves of that.

A refusal names the row. `FaultPlane::new` indexes the offending slickenline
within its own fault, which is always 0 here and tells a caller holding two
hundred rows nothing, so the binding prefixes the fault's index — and it does
refuse, rather than projecting the lineation onto the plane, because a
lineation off its plane is two measurements that do not belong together and
repairing it silently would invert data nobody measured.

`inversion_candidate_count` sizes a run before starting it. The count is not
the product a caller would guess — an axis pointing up is the axis pointing
down, so plunge spans a quarter turn and the roll of S3 about S1 half a one —
and halving both steps multiplies the work by about fourteen, which is worth
knowing before starting rather than after.

`score_stress` is the same misfit evaluated once, for a tensor that is not
being looked for but tested: a published one against your own faults, or a
runner-up against the subset you suspect belongs to a second phase. It returns
the per-fault misfits as well as their mean, and a fault whose plane sits on a
principal stress axis comes back `NaN` — no slip is predicted there, which is
not the same as a slip predicted and missed, and is why it is left out of the
mean rather than counted as zero.

It takes no `sigma1`/`sigma3`, and that absence is the statement: scaling a
tensor, or adding a multiple of the identity to it, leaves the shear direction
where it was, the second because a purely normal traction has no shear to
turn. So a misfit depends on the tensor's shape and nothing else, and offering
magnitudes would invite the reading that they change the answer.

`stress_tensor` is the matrix, `R . diag(sigma1, sigma2, sigma3) . R^T` in
(East, North, Up) — the bridge from the axes-and-Phi form everything here works
in to anything that wants the tensor itself. It exists because rebuilding it in
Python means repeating the construction of the right-handed triad by hand, S3
cross S1 and not the other way, which is a sign error waiting to happen. Unlike
the misfit, this one does scale with the magnitudes: the 1/0 default gives the
normalized tensor the search runs on, and the true values give the tensor.

## Density, and a stress field

`lib/src/spatial/` estimates a density over located observations, and
`lib/src/structural/stress_field.rs` uses the same machinery to invert a
tensor at every node of a grid. Both are generic over the dimension of the
field: two coordinates for a map of surface measurements, three for
hypocentres. The geology stays three-dimensional either way — `N` says how many
numbers locate an observation, not how many a fault plane needs.

`spatial::kernel` is the estimator. A `Bandwidth` carries one extent per axis,
in map units, and a `Kernel` is a Gaussian, a Gaussian truncated at so many
bandwidths, or a quartic. Two decisions there are worth stating, because the
alternative to each is what the C++ this descends from actually did.

The normalizing constant is **derived from the dimension** rather than
tabulated. A quartic kernel needs `15/16` on a line, `3/pi` on a plane and
`105/(32 pi)` in a volume, and `InterpDensity3D` carried the planar constant in
a volume kernel — so its output was not a density per unit volume, and two runs
at different bandwidths were not comparable with each other. Here the constant
comes out of a half-integer gamma recurrence at construction, so there is no
version of it to pick the wrong one of.

The bandwidth is in **map units and per axis**. In cells, it changes meaning
when the grid is refined, which is exactly when a reader is looking for the
answer to stop changing. Per axis, because a depth coordinate is not
interchangeable with a horizontal one even when both are in metres: a
seismogenic layer is far wider than it is thick, and an isotropic kernel over
one mixes the top of it with the bottom before it mixes two neighbours.

`spatial::grid::SamplingGrid` is a grid of sample points and deliberately not
of cells. `origin` *is* the first node. A raster header that says `ORIGIN` may
mean the corner of the first cell or its centre, and code that declares one
while computing the other produces a field displaced by half a spacing —
shifted rather than wrong-looking, and so survivable for years. There is
nothing here to be ambiguous about.

`spatial::density::Neighbourhood` bins the observations at the kernel's own
reach, so a node's cost follows its neighbourhood rather than the dataset. A
kernel with no reach — the untruncated Gaussian — is not refused but falls back
to visiting everything, which is the kernel to use when checking a result
against something else that did the same.

### The two together

`stress_field` walks a grid once, and at each node reports the density, the
weight of data that reached it, and the inverted tensor:

```rust
let kernel = Kernel::quartic(Bandwidth::isotropic(6_000.0).unwrap());
let cost = field_cost(&faults, &grid, kernel, SearchGrid::default(), 10.0);
// -> nodes, nodes_to_invert, forward_solutions

let field = stress_field(&faults, &grid, kernel, SearchGrid::default(), 10.0);
field.nodes[0].density        // faults per unit area, from the same kernel
field.nodes[0].support        // Kish's (sum w)^2 / sum w^2
field.nodes[0].solution       // the whole InversionResult, runners-up included
```

This is Hardebeck and Michael's (2006) spatially varying inversion in its
kernel form: rather than cutting the dataset into bins and damping neighbouring
solutions towards each other, every fault contributes to every node it reaches,
by an amount falling off with distance. Weighting is the only handle this search
has for separating two superposed tectonic phases, and in a field the weights
come from geography rather than from someone already knowing the answer.

The density is computed in the same traversal because the two answer each
other. A tensor from four faults is not the same object as one from two
hundred, and a stress map read without the density beside it is a smoothed
picture of where the data happen to be, presented as tectonics. `support` is
what to read: it counts twenty faults contributing equally as twenty, and
twenty of which nineteen sit on the far edge of the kernel as barely one. Nodes
below `min_support` keep their density and are left without a tensor, so a thin
patch stays legible instead of being filled with an answer the search will
always produce.

It does not follow that both belong on the same grid. On 300 faults over 40 km
of map, a 6 km quartic kernel and the default search, a node cost about 350 ms
to invert and microseconds to take a density at. So `density_field` stays
available on its own for the fine grid, and this pass is for the coarse one
where a tensor is actually wanted. The compact kernel is what makes even that
affordable: each node saw about 22 of the 300 faults, so a hundred-node field
cost four times a single inversion over the whole set rather than a hundred
times it. `field_cost` walks the same nodes with the same kernel in well under
a millisecond and reports the work the real pass would do — measured, not
estimated from a formula, so a dataset clustered in one corner of its own grid
is costed as such.

Both passes run across rayon's thread pool, and both give **bit for bit** what
one thread would. That is worth stating because it is the property parallel
numerics usually lose: a reduction that lets threads combine partial sums adds
the same terms in an order that depends on how the work was split, so the last
bits move when the machine changes. Nothing of the kind happens here. The split
is over nodes; a node's value is a sequential sum over its own neighbourhood in
an order fixed by the bin layout; threads never share an accumulator; and the
field's counters are summed from per-node outcomes, integer addition not caring
what order it happens in.

The split also keeps the cores fed, which is not automatic when nodes cost
anything from nothing to a full search: on the run above, 3.93 threads' worth of
CPU was consumed on four threads and 6.92 on eight. Wall-clock went from 21.6 s
to between five and seven, three to four times rather than seven, and the
shortfall is the machine and not the schedule — a 15 W mobile part holding its
all-core clock far below its single-core turbo, two hardware threads to a
physical core. Where the clock holds, the wall time follows the occupancy.
`rayon::ThreadPoolBuilder::install` bounds the pool where taking the machine is
not acceptable.

### Against the reference

`lib/tests/density_reference.rs` runs the estimator against `InterpDensity3D`
on the same points, checking the three-dimensional constants, the grid
registration and the node ordering at once — the three things this lineage has
been wrong about. It agrees to within five parts in ten million, which is half
a unit in the last digit the reference file prints: there is no digit left in
which the two could disagree.

The committed data are synthetic, from a seeded generator, with the extents
arranged so the grid comes out 8 by 6 by 10 — three different counts, so a
transposed traversal cannot pass by coincidence. `InterpDensity3D`'s own sample
dataset would have been the obvious choice and is deliberately not used: it is
somebody's earthquake catalogue arriving without a statement of its terms, and
`example_data/NOTICE.md` sets the rule that vendored data travels under licence
terms that can be named. That full volume was checked too, off to one side —
91 686 nodes, 1 525 hypocentres, agreeing to the same five parts in ten million
— but it is not this repository's to commit.

Note what the reference is: the *corrected* `InterpDensity3D`, at 0.9 or later.
The older one is wrong in ways an agreement would have propagated.

## Focal mechanisms, and the Kagan angle

`algebra::quaternion` and `structural::focal_mechanism` carry a double-couple
focal mechanism as its **P, T and B** kinematic axes, and
`structural::rotation` gives the rotation between two of them after Kagan
(1991).

Nodal planes are the usual way a mechanism is quoted and deliberately not what
this takes: nothing in the seismology says which of the two planes slipped, so
a pair of planes is an ambiguous input where P and T are not.
`PTBAxes::from_fault_slickenline` converts from a fault plane and its slip,
which is the unambiguous form of the same thing — and requires the sense of
movement to be known, because reversing an undetermined slip exchanges P with
T and returns shortening where the rock recorded extension.

**P and T are not stress axes.** They are kinematic, pinned at 45 degrees to
the fault plane by construction, and coincide with the principal stresses only
under Anderson's assumption. Where this crate means stress it says stress, and
`structural::inversion` is what solves for it.

### Four answers, not one

A double couple is unchanged by a half turn about any of its own axes, so the
question "what rotation carries mechanism 1 onto mechanism 2" has four equally
correct answers. `focal_mechanism_rotations` returns all four, sorted by the
size of the turn; the smallest is the **Kagan angle**, the standard measure of
how far apart two mechanisms are. The other three are kept because a minimum
quoted without its alternatives is a number nobody can check.

No two double couples are more than 120 degrees apart — 120 and not the 180 a
general pair of frames would allow — and that bound is the cheapest check there
is on an implementation of this. A wrong set of symmetry generators fails it at
once. The tests sweep 2000 random pairs against it, and against Kagan's own
published table of four solutions.

### What the Fortran had wrong here

The reference is `GeoFaults/` in this repository's own history — Alberti
(2010)'s `FaultCorrelation.f95` and the modules under it, removed in `961d8b9`
and recoverable with `git show 961d8b9^:GeoFaults/…`. There is no separate
repository, and no submodule.

`quaternfromcartmatr` there takes four square roots whose arguments are
non-negative in exact arithmetic and not in floating point. At a **trace of
-1** — a rotation by 180 degrees, which is what two strike-slip mechanisms a
whole number of degrees apart give — one argument that should be exactly zero
lands a few times 10⁻¹⁶ below it, and the Fortran returns NaN in silence. It
went unseen for twenty years because it needs the trace to be *exactly* -1:
round-degree catalogue data and textbook examples produce that constantly, and
random doubles essentially never. `Quaternion::from_rotation_matrix` clamps,
which is safe because the branch taken always divides by a component at or
above the mean of the four, so a clamped zero never reaches a denominator.

Everything else in that Fortran's Kagan path was verified rather than assumed:
`triadvectors_rotsolution` gives bit-identical output compiled with and without
`-fimplicit-none`, which is not true of its sibling `ForwardStress.f95`.

## Python bindings

`pylib` builds the `misah` Python distribution with maturin, against pyo3 0.29
and `abi3-py39`, so one wheel serves every Python from 3.9 on. That matters for
the QGIS side of the suite: a plugin cannot choose the interpreter it is loaded
into, nor count on a compiler being available to build a version-specific
extension.

Three names, deliberately distinct. The distribution installed from PyPI is
`misah`; the crate under `pylib` is `misah-py`, because the core crate in `lib`
already holds `misah` and a workspace admits one package per name; the compiled
extension is `misah._misah`, private because callers reach its contents through
the aliases `misah/__init__.py` sets up.

```sh
cd pylib
maturin develop --release   # build and install into the active environment
maturin build --release     # produce a wheel under target/wheels
cd ..
python3 pylib/tests/test_kernels.py   # never from pylib, which shadows the install
```

```python
import numpy as np
from misah.kernels import intersect_plane_grid

points, segments = intersect_plane_grid(
    np.ascontiguousarray(dem),   # C-contiguous, else ValueError
    geotransform,                # the six GDAL elements
    (x0, y0, z0),
    135.0,                       # dip direction
    35.0,                        # dip angle
    nodata,
)
```

```python
from misah.kernels import intersect_mesh_grid

# vertices (V, 3) and faces (F, 3) as a VTK POLYDATA file gives them
points, attitudes, mesh_triangles, stats = intersect_mesh_grid(
    np.ascontiguousarray(dem),
    geotransform,
    vertices,
    faces,
    nodata,
)
attitudes[0]        # (dip direction, dip angle) of the triangle that cut here
stats["duplicate_crossings"]
```

```python
from misah.kernels import best_fit_planes

# The points above, read back as attitudes. Any (N, 3) array will do.
field, stats = best_fit_planes(points, cell_size=50.0)

field["attitudes"]       # (M, 2) dip direction and dip angle, one per cell
field["cell_centres"]    # (M, 2) where to post them
field["collinearity"]    # (M,) how well each is determined; large is bad
stats["cells_collinear"] # cells set aside: on a smooth slope, most of them
```

`field` is a dict of eight aligned arrays rather than a tuple, which eight
positional returns would not survive being read. `None` where there is nothing
to fit.

```python
from misah.kernels import solve_stress

solution = solve_stress(
    0.0, 90.0,    # S1: trend, plunge
    90.0, 0.0,    # S3: trend, plunge
    0.5,          # Phi
    0.0, 60.0,    # fault: RHR strike, dip angle
    sigma1=30.0, sigma3=10.0,   # optional; default 1/0 leaves the rake correct
)
solution["theoretical_rake"]           # -90.0: pure normal dip-slip
solution["theoretical_slickenline"]    # (90.0, 60.0), down the dip vector
```

```python
from misah.kernels import invert_stress

# (N, 4): fault RHR strike and dip, then the trend and plunge of the
# slickenline on it. Every lineation must lie in its own plane to within a
# degree, or the row is named and refused.
result = invert_stress(
    faults,
    senses=None,             # (N,) bool: False where the sense was not read
    angle_step_degrees=10.0,
    phi_step=0.1,
)

result["best"]["s1"]                   # (trend, plunge)
result["best"]["phi"]
result["best"]["mean_misfit_degrees"]
result["best"]["faults_scored"]        # how many faults that mean came from
result["runners_up"]                   # the next five, worst last
```

`None` rather than a result when no candidate could be scored on any fault —
an empty set. The GIL is released for the search, so a hundred faults on the
default grid holds it for none of the 1.5 s they take.

```python
from misah.kernels import score_stress, stress_tensor

# The same tensor, evaluated rather than searched for.
score = score_stress(*best["s1"], *best["s3"], best["phi"], faults)
score["mean_misfit_degrees"]
score["misfits"]            # (N,) in the order given; NaN where none predicted

# And as a matrix, in (East, North, Up).
matrix = stress_tensor(*best["s1"], *best["s3"], best["phi"],
                       sigma1=30.0, sigma3=10.0)
```

```python
from misah.kernels import covering_grid, density_field, field_cost, stress_field

# A grid around the data. `origin` is the first sample point, not a cell
# corner, and `margin` should be about the kernel's reach.
grid = covering_grid(hypocentres, spacing=[500.0] * 3, margin=1500.0)

# Observations per unit volume at every node, flat, first axis fastest.
values = density_field(hypocentres, grid["origin"], grid["spacing"],
                       grid["counts"], bandwidth=[1500.0] * 3)
volume = values.reshape(grid["counts"], order="F")

# Summed over the grid and multiplied by the cell volume, that comes back as
# the number of observations. It is one line and it is worth running.
values.sum() * float(np.prod(grid["spacing"]))
```

```python
# What the field will cost, in milliseconds, before committing to it.
field_cost(positions, faults, origin, spacing, counts, [4000.0] * 2,
           min_support=10.0)
# -> nodes, nodes_to_invert, forward_solutions

field, stats = stress_field(positions, faults, origin, spacing, counts,
                            bandwidth=[4000.0] * 2, min_support=10.0)

field["density"]                  # faults per unit area, same kernel
field["support"]                  # Kish's (sum w)^2 / sum w^2 -- read this one
field["s1"]                       # (M, 2) trend and plunge, NaN where no tensor
field["runner_up_s1_degrees"]     # how far the nearest rival answer sits
stats["nodes_below_threshold"]    # holes left as holes
```

`positions` is an (M, D) array with D of 1, 2 or 3 — the dimension of the
*field*, not of the geology, read off the width of the array. Everything comes
back flat in the grid's own order, first axis fastest, so the columns align
with each other and with `density_field` on the same grid; `reshape(counts,
order="F")` lays one out.

Both run across every core and release the interpreter while they do. The
answer does not depend on the thread count, bit for bit: the split is over
nodes, each node's arithmetic stays sequential inside it, and no accumulator is
shared. `RAYON_NUM_THREADS` bounds the pool where taking the machine is not
acceptable — a QGIS plugin, say. Measured on a four-core laptop, one field took
7.54 s at one thread, 3.73 s at two and 1.51 s unbounded.

```python
from misah.kernels import (ptb_axes, kagan_angles, kagan_angle_matrix,
                           focal_mechanism_rotations)

# P, T and B axes from faults and their slip. Same (N, 4) array
# invert_stress takes -- but the sense of movement must be known here, or
# P and T swap and shortening is reported as extension.
axes = ptb_axes(faults)          # -> p, t, b, each (N, 2) trend and plunge

# A mechanism crosses as four numbers: P trend, P plunge, T trend, T plunge.
kagan_angles(first, second)      # (N,) elementwise, 0 to 120 degrees
kagan_angle_matrix(catalogue)    # (N, N), symmetric, zero diagonal

# All four rotations for one pair, smallest turn first.
focal_mechanism_rotations(232.0, 41.0, 120.0, 24.0,
                          51.0, 17.0, 295.0, 55.0)
# -> trend, plunge, angle_degrees, each (4,)
```

The matrix is why any of this is in Rust: comparing a catalogue with itself is
`N²/2` rotations, trivial once and ruinous in a Python loop at N in the
thousands. It runs across every core with the interpreter released — 600
mechanisms, 180 000 pairs, 32 ms. Symmetric with a zero diagonal, so it goes
straight into a clustering routine, which is the usual reason to want one.

Worked examples are in [`docs/notebooks`](docs/notebooks): kernel density in
`01_kernel_density.ipynb`, the stress field in `02_stress_field.ipynb`, focal
mechanisms and the Kagan angle in `03_focal_mechanisms.ipynb`. Each generates
its own data — through the forward model, or against Kagan's published table —
so the answer is known without trusting the code, and each carries its output so
it reads without being run.

The extension is built as `misah._misah`, so its submodules register themselves
under that name; `misah/__init__.py` aliases them, which is what makes
`import misah.kernels` work rather than only `misah._misah.kernels`.

There is no pure-Python fallback. There was one, `misah/_reference.py`,
mirroring every exposed kernel so that `misah.kernels` answered either way and
`misah.is_compiled` said which; it is gone, deliberately.

It cost a second independent implementation of every function exposed, and that
price is what kept the mesh-grid intersection unexposed for as long as it was:
marching triangles over a spatially indexed DEM is not something worth writing
twice, so the rule that was meant to protect the bindings was instead deciding
which half of the library Python could reach. Dropping it is what let
`intersect_mesh_grid` above be exposed at all.

What it bought was real and is worth naming: an oracle. Two implementations
written from the same formulas, compared elementwise, can say they *agree*
rather than that one of them runs — and it earned its keep, catching
`solve_stress` returning its vectors as lists from one side and tuples from the
other. In its place the bindings are checked against the datasets the kernels
were ported from: the geoSurfDEM golden trace for the raster ones, the five
Fortran-verified cases for the stress ones, and expected values rather than
mutual agreement throughout.

What it also bought, and this is the loss: a QGIS plugin that cannot install a
binary wheel could still use misah, slowly. It now cannot use it at all. That
is the same problem as the wheels, and belongs there — building for the
platforms QGIS runs on — rather than in a second implementation of everything
kept in reserve. Which is now what has happened: with wheels for Windows, macOS
and both Linux architectures, a plugin that cannot install one is a far narrower
case than when this was written. Narrow enough that if a fallback returns it
should be a selective one — the functions cheap to write twice, not the grid
search — and not the whole mirror again.

All seven Python suites run in CI. `test_stress.py`, `test_mesh.py`,
`test_inversion.py`, `test_best_fit.py`, `test_fields.py` and
`test_mechanisms.py` are numbers in, numbers out;
`test_inversion.py` needs no fixture at all, generating its faults through
`solve_stress` from a tensor chosen in advance and asking the inversion to find
that tensor again — the two halves of the Python surface checked against each
other, rather than either against numbers copied over from the Rust tests.

`test_kernels.py` needs a DEM, and used to read the Malpi crop from a checkout
of geoSurfDEM sitting beside this one, which meant the raster bindings — the
half of the surface with a golden trace to be measured against — were the half
nobody but the author could check. They are now checked by anyone with this
repository. The crop committed under `example_data` is a wider window on the
same ASTER tile, 213 x 260 against 200 x 247, and *wider on the same grid*:
the golden crop's corner is exactly seven cells in from it, east and north, and
the elevations agree value for value over the overlap. So the suite cuts the
200 x 247 window out of the committed file and intersects that, which is not a
substitute for the grid the trace was computed on but that grid itself — 404
vertices, unchanged, from data in this tree. The window is derived from the
golden crop's own corner coordinates rather than written down as four array
indices, so a committed DEM that stopped containing it would fail the suite
instead of quietly moving the trace onto different ground.

Run the tests from anywhere except `pylib` itself, whose source tree shadows the
installed package.

### Status

`raster::io` reads ESRI ASCII grids, which is what both examples now use; it is
the only raster format handled, anything wider meaning GDAL. It returns the
nodata value as an `Option` rather than defaulting to -9999 on a file's behalf,
since at sea that is a depth and not a hole.

`gp_projected_focal_mechanisms` is read from a GeoProfiler export into records
that nothing yet computes with. `structural::focal_mechanism` now carries the
P/T/B triad those records would feed, so the join is a matter of wiring rather
than of arithmetic.

The Python surface is twenty functions: `intersect_plane_grid`,
`intersect_mesh_grid`, `best_fit_planes` and `plane_normal` for the raster
side; `solve_stress`, `rake_to_slickenline`, `stress_tensor`, `score_stress`,
`invert_stress` and `inversion_candidate_count` for the structural side;
`covering_grid`, `kernel_profile`, `density_field`, `stress_field` and
`field_cost` for the fields; and `ptb_axes`, `kagan_angles`,
`kagan_angle_matrix`, `focal_mechanism_rotations` and
`rotate_focal_mechanism` for focal mechanisms. Each builds and consumes the types it needs within
the one call rather than handing them across the boundary, so the surface stays
plain arrays and numbers in, arrays and dicts out.

The field five are where that costs something worth naming. `N`, the dimension
an observation is located in, is a const generic on the Rust side — one compiled
copy per dimension, chosen when the code is written — and on the Python side it
is the width of an array, at run time. The bridge is a match over 1, 2 and 3,
and it is why a field here is one, two or three dimensional rather than any
number: each is a monomorphisation somebody has to name. One dimension is not
padding for a tidy range, either. A histogram of hypocentral depths is a
one-dimensional density, and estimating it with a kernel rather than with bins
removes both of a histogram's troubles at once — edges chosen by hand, and an
answer that changes when they move.

What remains Rust-only is the GeoProfiler SQLite reader, and deliberately:
Python has `sqlite3` in its standard library and can read a GeoProfiler export
directly, so the Rust reader exists for Rust callers rather than standing
between Python and the data.

`geoprofile::sqlite` reads a qgSurf GeoProfiler export — all thirteen tables of
it — and `lib/tests` runs that against two: a schema-v1 export of the Timpa San
Lorenzo section, and a small v5 one written by qgSurf's own exporter through its
public insert functions, regenerable with `lib/tests/data/make_export_synthetic_v5.py`.
The two projects agree on a file format that nothing else spans, so the fixture
is that agreement rather than a transcription of it. The four added after schema
v1 (`gp_projected_focal_mechanisms`, `gp_profile_vertices`,
`gp_graphical_params`, `gp_source_categories`) are read when present and come
back empty when not, so one reader serves every export written so far.

A line intersection needs both shapes. A line generally crosses a section at a
point, but a segment lying along the section trace crosses it over a stretch;
until qgSurf's schema v5 the table stored a single distance and such a stretch
was written as two ordinary rows, indistinguishable afterwards from two separate
crossings. v5 records `s_from` and `s_to`, as `gp_intersected_polygons` always
has. Both are read: the older column becomes the degenerate span it stands for.

The setuptools-rust leftovers are gone: `features.rs`, the empty
`georeferenced.rs` and `orientations.rs`, and `setup.py` with `MANIFEST.in`,
superseded by maturin and broken besides — `setup.py` used an `install_requires`
it never defined.

The pipeline is in two files, and the split is not arbitrary. `.gitlab-ci.yml`
runs everything on Linux at every push: the workspace and the examples built,
clippy as a gate, and the wheel installed so that all seven Python suites run
against it — the bindings imported and called, not merely compiled.
`.github/workflows/CI.yml` runs those same seven suites on macOS and Windows,
which is the one thing a GitLab Free namespace structurally cannot do, and on
tags builds every artifact and uploads it.

GitHub is a push mirror and nothing is developed there; it holds that job
because its macOS and Windows runners are free and unmetered on public
repositories. The two Linux wheels went with the rest rather than staying
behind, and not for tidiness: a version is final once uploaded, so two pipelines
able to upload the same one is a hazard rather than a redundancy. One place
builds every artifact and one place uploads them, and GitLab keeps what only it
has — the whole Rust side, and `cargo publish`.

Both uploads are manual, for that same finality: the tag builds and checks the
artifacts, and someone then decides. `git push --tags` stays reversible;
pressing the button does not. PyPI goes through Trusted Publishing, so no token
is stored in either repository — the runner is minted a short-lived OIDC token
and the publishing action trades it for one PyPI will accept. A publisher is
registered per PyPI project rather than per account, so this one is its own, and
what it names is the workflow *file*, `CI.yml`, rather than the `name:` inside
it. The manual button has no equivalent on GitHub; a *required reviewer* on the
`release` environment is what stands in for it, and without one `git push
--tags` becomes the irreversible act. crates.io takes the kernel crate alone;
`misah-py` is the extension, depends on `misah` by path rather than by version,
and is not something anyone adds to a Cargo.toml.

The two registries are independent, and nothing here makes one wait for the
other: the sdist carries `lib/` inside it, so the Python distribution builds
without `misah` ever reaching crates.io.

### One version number

`[workspace.package] version` in the root `Cargo.toml`, inherited by both
members. They were maintained separately and had drifted — the crate at `0.2.0`
while the bindings were at `0.2.0-alpha.0` — which a tag cannot express, a
release being one commit. Both are now `0.2.0-alpha.2`, which maturin normalizes
to `0.2.0a2` for PyPI.

The alpha is deliberate, and what it costs an installer is less than it looks —
which is worth recording, because this file first claimed the opposite. pip is
said to skip pre-releases unless asked with `--pre`, and the rule is narrower
than that: it skips them only when a stable version exists to prefer instead.
misah has none, so `pip install misah` installs the alpha, verified by doing it
from PyPI into a clean virtualenv. The flag becomes necessary the day a stable
release exists and an alpha after it is wanted.

Cargo is the stricter of the two, and there the usual statement does hold: a
requirement like `misah = "0.2"` does not match `0.2.0-alpha.2`, since a version
range admits a pre-release only when the range itself names one. Until a stable
version is published, a dependent has to write the pre-release out in full.

**The number moves before a release, not at it.** `0.2.0-alpha.1` was published
on 29 August 2026 and then stood still while the density fields, the stress
fields and the focal mechanisms went in behind it. Nothing in the tooling
objects to that, which is what makes it worth naming: a wheel built from the
tree afterwards calls itself `0.2.0a1` as well, installs cleanly over the
published one, and answers to `density_field` — which the version on PyPI does
not have. The version string is the only thing that distinguishes them, so
testing against a local build says nothing about what an installer receives
unless the two numbers are known to differ.

### Wheels

A wheel built plainly on a development machine takes the glibc it happens to
find: here that produced a `manylinux_2_34` tag, which pip will refuse to
install on Ubuntu 20.04, Debian 11 or RHEL 8. That used to be the whole story,
and the conclusion drawn from it — that nothing could reach a QGIS other than
this one without macOS, Windows and aarch64 runners — was right about needing
the runners and wrong in supposing they could not be had.

There are five artifacts: `manylinux_2_17` wheels for x86-64 and aarch64, a
universal2 wheel for macOS, an MSVC wheel for Windows, and the source
distribution. `abi3-py39` is what keeps that number down — the extension links
no version-specific CPython symbol, so it is one wheel per *platform* rather
than one per platform and interpreter version, five artifacts instead of thirty.

**The compatibility floor is a decision and not an accident.** Which glibc a
Linux wheel demands is otherwise whatever the build image happened to carry.
Four ways of building the same extension, and the floor each one produces:

| built by | tag | glibc actually required |
| --- | --- | --- |
| the host, plainly | `manylinux_2_34_x86_64` | 2.34 |
| a `manylinux_2_28` container | `manylinux_2_28_x86_64` | 2.28 |
| zig, x86-64 | `manylinux_2_17_x86_64` | 2.14 |
| zig, aarch64 | `manylinux_2_17_aarch64` | 2.17 |

`maturin --zig` is how that floor was first reached, cross-compiling both Linux
wheels from the single x86-64 runner a GitLab Free namespace gets, and for that
situation it was the right tool rather than a workaround for lacking runners.
It buys nothing where the wheels are built now — a `manylinux_2_17` container
reaches the same floor on its own — but the floor was kept when the jobs moved,
which is the point of the table: chosen once, and inherited afterwards on
purpose.

Kept, and also checked. `ci/check_wheel.py` opens the extension inside the wheel
and reads its own header rather than its file name: the machine type, so that a
build which quietly produced a host binary is caught, and on Linux the versioned
glibc symbols it references, so the tag becomes a claim the binary is held to.
On macOS it reads Mach-O, thin and fat, and reports universal2 only when the
pair of slices is exactly that — which matters more there than anywhere else,
since that wheel is built on the arm64 runner and its x86-64 half is compiled
and never executed. An arm64-only build would otherwise pass every test there
is.

**Windows and macOS were the two holes, and both turned out to be problems of
cross-compilation rather than of compilation.** That is why native runners
closed them at a stroke, and why they had stayed open for so long against so
much effort. The failed attempts are still worth recording, because neither
failed for the reason one would guess.

`rusqlite` is pulled in with the `bundled` feature, so every build compiles
SQLite's C amalgamation — which makes a Windows *cross*-build need a C compiler
for Windows, not only a Rust target. That is invisible on the Linux side
because zig ships a C compiler; it was the whole difficulty on the Windows one.

- **`x86_64-pc-windows-msvc` through cargo-xwin** downloads the Microsoft CRT
  and SDK — around 2.5 GB — and then stops at `failed to find tool "clang-cl"`.
  cargo-xwin supplies the headers and libraries but not the toolchain, which
  wants `clang-cl`, `lld-link` and `llvm-lib` from LLVM. The download alone is
  a real cost against a metered budget, and large enough to sit awkwardly in a
  CI cache.
- **`x86_64-pc-windows-gnu` through zig** gets further — SQLite compiles — and
  fails at the link with `undefined symbol: PyInit_misah._misah`. That is this
  package's `module-name = "misah._misah"` reaching the export-definition file
  with its dot intact, where the initialiser Python actually looks for is
  `PyInit__misah`. It is also the wrong target to want: CPython on Windows is
  built with MSVC, and a GNU-ABI extension is not the supported configuration.

On `windows-latest` neither arises. MSVC is the native toolchain there, so the
target is the ordinary one and `rusqlite` finds the C compiler it wanted all
along — no CRT to download, no LLVM to install, no ABI to argue about. macOS
was a licensing question before a technical one, and there too the question is
narrower than it sounds: cross-compiling needs the Apple SDK, and what the
licence speaks to is *extracting* one in order to build from Linux. It does not
arise when the compiler runs on Apple hardware. Both holes were closed by moving
the build rather than by solving what had been attempted.

**The source distribution is what everything else gets.** With no wheel
matching, pip falls back to it and builds, which needs a Rust toolchain on the
user's machine but no runner here — a harder install rather than no install.
macOS and Windows no longer depend on that, having wheels of their own, but the
sdist job earns its place for a second reason: `pip install` on the tarball is
the only thing anywhere that exercises maturin's rewriting of the manifests when
it packs a project whose crate lives outside the Python directory, and with the
version inherited from the workspace that rewriting is a claim to be checked
rather than assumed. So the job installs *from the tarball* and never from the
tree, and runs all seven suites against that.

`docs/notebooks/misah.ipynb` was removed rather than repaired. Every path in it
was dead, not merely the one previously named here: `misah.geometry.geom2d`,
`misah.geometry.geom3d`, `misah.orientations.orien3d`,
`misah.georeferenced.georef2d` — the last of which was a file deleted two
commits before that line was written — and it compared results against `pygsf`,
which geogst itself superseded. Nothing in it touched any function the package
exposes.

`01_kernel_density.ipynb` and `02_stress_field.ipynb` are the new ones that
line asked for. They are not a repair of the old one and share nothing with it:
they call the functions that exist, they read no files, and they build their
data through the forward model so that what the inversion should return is
known without trusting the inversion. `docs/README.md` says how they are kept
from going the way of the first.

`lib/examples/test.rs` is gone with it: a nalgebra hello-world computing the
distance between two points, referring to nothing in this crate, and the sole
reason `nalgebra` appeared in `lib/Cargo.toml` — 97 lines of lockfile and a
whole dependency tree for an example that demonstrated nothing about misah.

Both crates are now on edition 2021. `lib` was on 2018, which is what let
`assert!(cond, "…{err}")` pass as a plain string rather than a format string,
and so print the placeholder instead of the error at the moment a test failed;
ten such assertions had to be corrected by hand once clippy pointed at them.
The edition closes the class rather than the instances. `cargo fix --edition`
needed no source changes to make the move.

## Licence, and the data

GPL-3.0-or-later, and the text sits in three places: the repository root, `lib/`
and `pylib/`. The copies are not redundancy but the consequence of how the
artifacts are cut. `cargo package` takes one directory and nothing above it, so
the SPDX expression in `lib/Cargo.toml` would have reached crates.io without the
licence it names; the wheel is the same case, which is why `pylib` already
carried one.

The data under `example_data` is not the repository's to license and travels
under its own terms, recorded in `example_data/NOTICE.md`. The Monte Alpi DEM is
derived from the global ASTER GDEM, which permits redistribution and asks for an
acknowledgement in return: **ASTER GDEM is a product of METI and NASA**.
Committing the file is redistribution, so the acknowledgement belongs in the
tree rather than in the citation of a paper that has not been written.
