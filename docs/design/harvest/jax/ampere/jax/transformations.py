# ============================================================================
# ASPIRATIONAL SKETCH -- harvested read-only from origin/jax (W0.7, 2026-09).
# DOES NOT IMPORT. MUST NOT BE MERGED AS-IS.
# ----------------------------------------------------------------------------
# ampere/jax/module.py (this file or one it imports) does
#     from equinox import Dataset as EqxDataset
#     from equinox import Parameter as EqxParameter
# Neither `equinox.Dataset` nor `equinox.Parameter` exists in the equinox
# public API (checked against equinox on PyPI as of this harvest) -- this
# import fails immediately. Every other module in ampere/jax/ imports
# (directly or transitively) from ampere/jax/module.py, so the whole
# ampere.jax sketch package fails to import.
# This is reference-only design intent for a future jax backend
# (see DEVELOPMENT_PLAN.md); treat it as sketch notes, not implementation.
# Full provenance: docs/design/harvest/jax/README.md
# ============================================================================

""" This module defines the base class for transformations from the output of a
model to the output of an instrument model.

The `Transformation` class is a subclass of `Module`, and is intended to be used
as a base class for all transformations in the `ampere.jax` namespace. It
provides a `__call__` method that should be implemented in a subclass.

We also define a few example transformations as subclasses of `Transformation`,
including `SyntheticPhotometry`, `LineSpreadFunction`, `ResampleSpectrum` and
`CalibrateSpectrum`. These are intended to be used in instrument models, and to
be extended to include more complex transformations as needed.
"""

from typing import Union, List, Dict, Any
import jax
import jax.numpy as jnp
from .module import Module

class Transformation(Module):
    """
    Base class for all transformations in `ampere.jax`. This class is intended to
    be used as a base class for all transformations, and provides a `__call__`
    method that should be implemented in a subclass.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __call__(self, *args, **kwargs):
        raise NotImplementedError("The __call__ method must be implemented in a subclass.")
    
class SyntheticPhotometry(Transformation):
    pass

class LineSpreadFunction(Transformation):
    pass