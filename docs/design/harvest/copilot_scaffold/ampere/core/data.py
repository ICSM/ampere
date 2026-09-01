"""
ampere.core.data
================

Abstract base class for all data objects in ampere.

This module defines the stable interface that every data container must
satisfy.  Concrete implementations live in ``ampere.data``.

Classes
-------
Data
    Abstract base; declares the interface only.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Dict, Union

import numpy as np

if TYPE_CHECKING:
    from .parameter import ParameterSet
    from ..models.results import ModelResults


class Data(ABC):
    """Abstract base class for all ampere data containers.

    Subclasses must:

    1. Declare a :attr:`parameters` attribute of type
       :class:`~ampere.core.parameter.ParameterSet` for nuisance parameters
       (e.g. noise-model hyperparameters, calibration offsets).
    2. Implement :meth:`lnlike`.

    The default implementations of :meth:`lnprior` and
    :meth:`prior_transform` delegate to ``self.parameters``.
    """

    # Subclasses must assign a ParameterSet here (or in __init__).
    parameters: "ParameterSet"

    @abstractmethod
    def lnlike(
        self,
        theta: Union[np.ndarray, Dict[str, float]],
        result: "ModelResults",
    ) -> float:
        """Evaluate the log-likelihood given nuisance parameters and model output.

        Parameters
        ----------
        theta : np.ndarray or dict
            Nuisance parameter values.  May be supplied either as a flat
            numpy array (for legacy samplers) or as a name→value dict (for
            the new dict-based API).
        result : ModelResults
            Output of the physical model evaluation.
        """

    def lnprior(self, theta: np.ndarray) -> float:
        """Evaluate the joint log-prior over nuisance parameters.

        Delegates to ``self.parameters.lnprior``.
        """
        return self.parameters.lnprior(theta)

    def prior_transform(self, u: np.ndarray) -> np.ndarray:
        """Map unit-hypercube values to nuisance parameter values.

        Delegates to ``self.parameters.prior_transform``.
        """
        return self.parameters.prior_transform(u)
