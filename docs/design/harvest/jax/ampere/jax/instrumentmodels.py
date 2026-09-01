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

""" This module defines the base class for instrument models in `ampere.jax`,
along with a few example instrument models. """ 

from typing import Union, List, Dict, Any
import jax
import jax.numpy as jnp
from .module import Module
from .transformations import Transformation

class InstrumentModel(Module):
    """
    Base class for all instrument models in `ampere.jax`. This class is intended
    to be used as a base class for all instrument models, and provides a
    `__call__` method that should be implemented in a subclass.
    """

    transformations: List[Transformation]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __call__(self, *args, **kwargs):
        raise NotImplementedError("The __call__ method must be implemented in a subclass.")
    
    def add_transformation(self, transformation: Transformation):
        self.transformations.append(transformation)