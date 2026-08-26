"""Python bindings for misah.

`misah.kernels` is the compiled extension where it can be built, and the pure
Python implementation in `_reference` where it cannot -- QGIS being the case that
matters, since a plugin neither picks the interpreter it is loaded into nor can
count on a binary wheel installing. Both expose the same names and return the
same arrays, so callers need not know which one answered; `is_compiled` tells
them if they care.

The extension is built as `misah._misah`, so its submodules register under that
name. Aliasing here is what lets callers write `import misah.kernels`.
"""

import sys

try:
    from ._misah import kernels

    is_compiled = True
except ImportError:  # pragma: no cover - exercised only without the extension
    from . import _reference as kernels

    is_compiled = False

sys.modules[f"{__name__}.kernels"] = kernels

__all__ = ["kernels", "is_compiled"]
