"""Python bindings for misah.

`misah.kernels` is the compiled extension. There is no pure-Python fallback:
there was one, mirroring every exposed kernel, and it was dropped deliberately.
It cost a second implementation of each function -- written independently, so
that the tests could say the two agreed rather than that one of them ran -- and
that price is what kept the mesh-grid intersection unexposed for as long as it
was, marching triangles over a spatially indexed DEM not being something worth
writing twice. The kernels are verified against the datasets they were ported
from instead: the geoSurfDEM golden trace, and the five cases checked against
ForwardStress.f95.

What that gives up is the QGIS case the fallback existed for -- a plugin that
cannot install a binary wheel now cannot use misah at all, rather than falling
back to something slower. That is the same problem as the wheels themselves,
and was to be solved there, by building for the platforms QGIS runs on, not by
keeping a second implementation of everything in reserve. It has been: there
are wheels now for Windows, macOS and both Linux architectures, which leaves a
far narrower case than the one this paragraph was written about.

The extension is built as `misah._misah`, so its submodules register under that
name. Aliasing here is what lets callers write `import misah.kernels`.
"""

import sys

from ._misah import kernels

sys.modules[f"{__name__}.kernels"] = kernels

__all__ = ["kernels"]
