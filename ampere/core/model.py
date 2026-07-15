"""
ampere.core.model
=================

Abstract base class for all physical models in ampere.

This module defines the stable interface that every model must satisfy.
Concrete implementations live in ``ampere.models``.

Classes
-------
Model
    Abstract base; declares the interface only.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from .parameter import ParameterSet
    from ..models.results import ModelResults


class Model(ABC):
    """Abstract base class for all ampere physical models.

    Subclasses must:

    1. Declare a :attr:`parameters` attribute of type
       :class:`~ampere.core.parameter.ParameterSet`.
    2. Implement :meth:`__call__` accepting keyword arguments whose names
       match ``self.parameters.names``.

    The default implementations of :meth:`lnprior` and
    :meth:`prior_transform` delegate to ``self.parameters``; subclasses only
    need to override them for non-standard behaviour.

    Notes
    -----
    Legacy models that still use positional ``*args`` should inherit from
    :class:`~ampere.models.models.LegacyModel` (defined in
    ``ampere.models.models``) which provides a backward-compatibility shim.
    """

    # Subclasses must assign a ParameterSet here (or in __init__).
    parameters: "ParameterSet"

    @abstractmethod
    def __call__(self, **params) -> "ModelResults":
        """Evaluate the model and return a :class:`~ampere.models.results.ModelResults`.

        Parameters
        ----------
        **params
            One keyword argument per parameter in ``self.parameters``.
        """

    def lnprior(self, theta: np.ndarray) -> float:
        """Evaluate the joint log-prior at *theta*.

        Delegates to ``self.parameters.lnprior``.
        """
        return self.parameters.lnprior(theta)

    def prior_transform(self, u: np.ndarray) -> np.ndarray:
        """Map unit-hypercube values to parameter values.

        Delegates to ``self.parameters.prior_transform``.
        """
        return self.parameters.prior_transform(u)

    # ------------------------------------------------------------------
    # Extension point for hierarchical / population inference
    # ------------------------------------------------------------------

    def population_model(self):
        """Return a population-level numpyro model, or ``None``.

        This hook is reserved for future hierarchical inference support
        (Issue 8).  The default implementation returns ``None``.
        """
        return None
