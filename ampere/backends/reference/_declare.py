"""Coercing a constructor argument into a :class:`~ampere.core.Parameter`.

Every model, step and noise model in this backend takes each of its physical
quantities the same three ways — a prior to fit it, a number to hold it fixed,
or a ready-made :class:`~ampere.core.Parameter` for anything more specific
(a tie, a declared unit, a non-default bijection). Doing that in one place is
what keeps the classes themselves short enough to read as physics.

``ampere.core.likelihood`` has a private helper of the same shape for kernel
hyperparameters. This is deliberately not that one: it is a backend concern,
its error messages talk about model parameters rather than kernel
hyperparameters, and importing a private name across the core boundary would
make the backend depend on an implementation detail of the contract.
"""

from __future__ import annotations

from typing import Any

import astropy.units as u
import numpy as np

from ampere.core import Parameter, ParameterError

__all__ = ["as_parameter"]


def as_parameter(
    name: str,
    given: Any,
    *,
    unit: u.UnitBase | None = None,
) -> Parameter:
    """Coerce *given* into a :class:`~ampere.core.Parameter` called *name*.

    Parameters
    ----------
    name
        The declared parameter name. A supplied ``Parameter`` must already
        carry it: names are part of a model's interface — they arrive at
        ``evaluate`` as keyword arguments and appear in the posterior — so
        silently renaming one would make the declaration a lie.
    given
        A frozen ``scipy.stats`` distribution (fitted), a number (held fixed),
        or a ``Parameter`` (used as-is).
    unit
        Declared unit for the parameter, if it has one.

    Notes
    -----
    No bijection is passed: :func:`~ampere.core.default_bijection_for` reads it
    off the prior's own support, which is the contract's machinery for exactly
    this and is right more often than a per-class guess. A caller wanting
    something else builds the ``Parameter`` itself and passes that.

    Raises
    ------
    ParameterError
        If *given* is none of the three, or is a ``Parameter`` under another
        name.
    """
    if isinstance(given, Parameter):
        if given.name != name:
            raise ParameterError(
                f"parameter {name!r} was given a Parameter named {given.name!r}. A model's "
                f"parameter names are part of its interface — they arrive at evaluate() as "
                f"keyword arguments — so they are fixed; rename it with .rename({name!r})."
            )
        return given
    if isinstance(given, (int, float, np.floating, np.integer)) and not isinstance(given, bool):
        return Parameter(name, value=float(given), fixed=True, unit=unit)
    if hasattr(given, "ppf"):
        return Parameter(name, given, unit=unit)
    raise ParameterError(
        f"parameter {name!r} must be a frozen scipy.stats distribution (to fit it), a number "
        f"(to hold it fixed), or an ampere Parameter — got {type(given).__name__}."
    )
