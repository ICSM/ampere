"""Coercing a constructor argument into a :class:`~ampere.core.Parameter`.

The jax backend's copy of the reference backend's helper, and deliberately a
copy rather than an import: ``ampere.backends.reference._declare`` is a private
name in a *different* backend, and a backend that reached across for it would
make the two subpackages share an implementation detail neither contract says
they share. The behaviour is identical because the declaration is identical —
what differs between the backends is the arithmetic, not the vocabulary.

Every model, step and noise model in this backend takes each of its physical
quantities the same three ways — a prior to fit it, a number to hold it fixed,
or a ready-made :class:`~ampere.core.Parameter` for anything more specific
(a tie, a declared unit, a non-default bijection).
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
    off the prior's own support, which is what makes the reference path and this
    backend's ``biject_to`` lowering agree by construction rather than by
    coincidence (``lowering.md`` §4).

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
