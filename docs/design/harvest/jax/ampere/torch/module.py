# ============================================================================
# ASPIRATIONAL SKETCH -- harvested read-only from origin/jax (W0.7, 2026-09).
# DOES NOT IMPORT AS A PACKAGE. MUST NOT BE MERGED AS-IS.
# ----------------------------------------------------------------------------
# ampere/torch/data.py (source branch filename: "ampere/torch/data,py" --
# comma typo, renamed here for clarity) does `from .likelihoods import
# Likelihood`, but ampere/torch/likelihoods.py in this same sketch is an
# empty file (0 bytes, no `Likelihood` symbol) -- that import fails.
# ampere/torch/module.py depends on `gpytorch` and `pyro`, optional
# dependencies not declared anywhere in this project's extras.
# This is reference-only design intent for a possible future torch backend
# (see DEVELOPMENT_PLAN.md); treat it as sketch notes, not implementation.
# Full provenance: docs/design/harvest/jax/README.md
# ============================================================================

""" This module defines the base class for all `ampere.torch` objects. """

from typing import Union, List, Dict, Any
import torch
from gpytorch import Module as GPyTorchModule
import pyro
from pyro.nn.module import PyroModule, PyroSample, PyroParam

class Module(GPyTorchModule, PyroModule):
    """
    Base class for all `ampere.torch` objects, including models, likelihoods,
    datasets and inference algorithms.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def forward(self, *args, **kwargs):
        raise NotImplementedError("The forward method must be implemented in a subclass.")

    def __call__(self, *args, **kwargs):
        return self.forward(*args, **kwargs)
    

class Parameter(torch.nn.Parameter):
    """
    Base class for all `ampere.torch` parameters. This class provides an interface
    to pyro primitives, and should be used as a base class for all parameters in
    the `ampere.torch` namespace.
    """

    _name: str
    _prior: Any
    _constraint: Any

    def __init__(self, name: str, prior: Any, constraint: Any = None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._name = name
        self._prior = prior
        if constraint is not None:
            self._constraint = constraint

    def sample(self):
        return pyro.sample(self._name, self._prior, obs=None)