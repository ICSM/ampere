"""
ampere.core.parameter
=====================

Named, prior-equipped parameter containers that replace manual ``npars``
counting and positional ``theta`` array slicing throughout ampere.

Classes
-------
Parameter
    A single named scalar parameter with a prior distribution and an optional
    constraint slot (reserved for future JAX/bijector support).
ParameterSet
    An ordered collection of :class:`Parameter` objects for one model or data
    object.  Provides pack/unpack helpers for numpy-based samplers and a
    ``as_paramax()`` method for JAX pytree support.

Notes
-----
``scipy.stats`` frozen distributions are the default prior representation.
A ``numpyro`` distribution may be supplied instead; it will be used
automatically by the future numpyro backend.

JAX / `paramax <https://github.com/danielward27/paramax>`_ are **optional**.
All functionality works with plain numpy; only ``ParameterSet.as_paramax()``
requires JAX and raises ``ImportError`` when it is absent.
"""

from __future__ import annotations

import dataclasses
from typing import Any, Dict, List, Optional

import numpy as np


@dataclasses.dataclass
class Parameter:
    """A single named scalar parameter with a prior distribution.

    Parameters
    ----------
    name : str
        Human-readable identifier, e.g. ``"temperature"`` or ``"log_tau"``.
    prior : object
        A ``scipy.stats`` frozen distribution (or any object that exposes
        ``.logpdf(x)`` and ``.ppf(u)``).  A ``numpyro`` distribution may be
        supplied instead and will be used by the future numpyro backend.
    constraint : object, optional
        A bijector / constraint object (e.g. from ``paramax.constraints``).
        Reserved for future use; has no effect in v1.
    value : float, optional
        Default or initial value.  ``None`` means uninitialised.
    hyperpriors : ParameterSet, optional
        Reserved extension point for hierarchical inference.  Leave as
        ``None`` for standard single-level inference.
    """

    name: str
    prior: Any
    constraint: Optional[Any] = None
    value: Optional[float] = None
    hyperpriors: Optional["ParameterSet"] = None

    def lnprior(self, x: float) -> float:
        """Evaluate the log-prior density at *x*."""
        return float(self.prior.logpdf(x))

    def ppf(self, u: float) -> float:
        """Percent-point function (inverse CDF) used by prior transforms."""
        return float(self.prior.ppf(u))


class ParameterSet:
    """Ordered, named collection of :class:`Parameter` objects.

    Used by every :class:`~ampere.core.model.Model` and
    :class:`~ampere.core.data.Data` subclass to declare the parameters they
    own.  Replaces the manual ``self.npars`` integer and the associated
    positional ``theta`` slicing in inference backends.

    Parameters
    ----------
    parameters : list of Parameter
        Parameters in declaration order.

    Examples
    --------
    >>> from scipy.stats import uniform, norm
    >>> pset = ParameterSet([
    ...     Parameter("temperature", uniform(loc=100, scale=9900)),
    ...     Parameter("log_tau",     norm(loc=0, scale=1)),
    ... ])
    >>> theta = pset.prior_transform(np.array([0.5, 0.5]))
    >>> lp = pset.lnprior(theta)
    """

    def __init__(self, parameters: List[Parameter]) -> None:
        self._params: List[Parameter] = list(parameters)
        self._index: Dict[str, int] = {
            p.name: i for i, p in enumerate(self._params)
        }

    # ------------------------------------------------------------------
    # Sequence helpers
    # ------------------------------------------------------------------

    @property
    def size(self) -> int:
        """Number of parameters in this set."""
        return len(self._params)

    # Alias used by legacy code that reads ``npars``.
    @property
    def npars(self) -> int:
        return self.size

    def __len__(self) -> int:
        return self.size

    def __iter__(self):
        return iter(self._params)

    def __getitem__(self, key):
        if isinstance(key, str):
            return self._params[self._index[key]]
        return self._params[key]

    @property
    def names(self) -> List[str]:
        """Parameter names in declaration order."""
        return [p.name for p in self._params]

    # ------------------------------------------------------------------
    # Pack / unpack
    # ------------------------------------------------------------------

    def pack(self, values: Dict[str, float]) -> np.ndarray:
        """Convert a name→value mapping into a flat numpy array.

        Parameters
        ----------
        values : dict
            Mapping from parameter name to scalar value.

        Returns
        -------
        np.ndarray
            1-D array of length ``self.size``.
        """
        return np.array([values[p.name] for p in self._params], dtype=float)

    def unpack(self, theta: np.ndarray) -> Dict[str, float]:
        """Convert a flat numpy array into a name→value mapping.

        Parameters
        ----------
        theta : np.ndarray
            1-D array of length ``self.size``.

        Returns
        -------
        dict
            Mapping from parameter name to scalar float.
        """
        if len(theta) != self.size:
            raise ValueError(
                f"Expected array of length {self.size}, got {len(theta)}"
            )
        return {p.name: float(theta[i]) for i, p in enumerate(self._params)}

    # ------------------------------------------------------------------
    # Prior evaluation
    # ------------------------------------------------------------------

    def lnprior(self, theta: np.ndarray) -> float:
        """Evaluate the joint log-prior at *theta*.

        Returns ``-inf`` as soon as any individual prior is non-finite.
        """
        lp = 0.0
        for i, p in enumerate(self._params):
            lp += p.lnprior(theta[i])
            if not np.isfinite(lp):
                return -np.inf
        return lp

    def prior_transform(self, u: np.ndarray) -> np.ndarray:
        """Map unit-hypercube values to parameter values.

        Used by nested-sampling backends (dynesty, …).

        Parameters
        ----------
        u : np.ndarray
            1-D array of values in ``[0, 1]``, one per parameter.

        Returns
        -------
        np.ndarray
            Corresponding parameter values.
        """
        return np.array(
            [p.ppf(float(u[i])) for i, p in enumerate(self._params)]
        )

    # ------------------------------------------------------------------
    # JAX / Paramax integration (optional)
    # ------------------------------------------------------------------

    def as_paramax(self):
        """Return a Paramax pytree representation of this parameter set.

        Requires ``paramax`` and ``jax`` to be installed.

        Raises
        ------
        ImportError
            If ``paramax`` or ``jax`` are not available.
        NotImplementedError
            Until the Paramax integration is implemented in Issue 2.
        """
        try:
            import paramax  # noqa: F401
        except ImportError as exc:
            raise ImportError(
                "paramax is required for JAX pytree support. "
                "Install it with:  pip install paramax"
            ) from exc
        raise NotImplementedError(
            "as_paramax() will be completed in Issue 2 "
            "(branch feature/core-parameter)."
        )

    # ------------------------------------------------------------------
    # Merge helper (used by CompositeModel)
    # ------------------------------------------------------------------

    @staticmethod
    def merge(
        sets: List["ParameterSet"], prefixes: List[str]
    ) -> "ParameterSet":
        """Merge multiple ParameterSets with name-scoping.

        Each parameter name becomes ``"{prefix}.{original_name}"``, avoiding
        collisions when two component models share a parameter name.

        Parameters
        ----------
        sets : list of ParameterSet
        prefixes : list of str
            One prefix per set.

        Returns
        -------
        ParameterSet
        """
        merged: List[Parameter] = []
        for pset, prefix in zip(sets, prefixes):
            for p in pset:
                merged.append(dataclasses.replace(p, name=f"{prefix}.{p.name}"))
        return ParameterSet(merged)
