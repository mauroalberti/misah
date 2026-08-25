"""Python bindings for misah.

The compiled extension is built as `mispy.mispy`, so its submodules register
themselves under that name. Aliasing them here is what lets callers write
`import mispy.kernels` instead of `import mispy.mispy.kernels`.
"""

import sys

from .mispy import kernels

sys.modules[f"{__name__}.kernels"] = kernels

__all__ = ["kernels"]
