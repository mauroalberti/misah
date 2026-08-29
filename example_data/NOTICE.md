# Data in this directory

The code in this repository is GPL-3.0-or-later. The data here is not the
repository's to license, and travels under its own terms.

## `raster/monte_alpi/malpi_aster_w4u3.asc`

A 213 x 260 crop over Monte Alpi (Basilicata, southern Italy), 27.28 m cells,
WGS 84 / UTM zone 33N (EPSG:32633), derived from the global ASTER GDEM and
reprojected. It reaches this repository through
[geoSurfDEM](https://github.com/mauroalberti/geoSurfDEM), whose own crop of the
same tile — 200 x 247 — is the one the golden-trace comparisons in the README
were run against; the two are not interchangeable, since the expected vertex
counts depend on the row and column counts.

ASTER GDEM is a product of METI and NASA.

That sentence is the acknowledgement the METI/NASA terms ask of anyone
redistributing the data, which is what committing this file does. The data
itself is distributed at no cost and may be redistributed; see
<https://asterweb.jpl.nasa.gov/gdem.asp> for the current statement of the
policy, which governs this file rather than anything written here.

The elevations are used here as a surface to cut, not as a measurement: what
the tests assert is where a plane of known attitude meets this grid, which the
grid's own vertical accuracy does not enter.
