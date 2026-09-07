"""Coercing a constructor argument into a :class:`~ampere.core.Parameter`.

Every model, step and noise model in this backend takes each of its physical
quantities the same three ways — a prior to fit it, a number to hold it fixed,
or a ready-made :class:`~ampere.core.Parameter` for anything more specific (a
tie, a declared unit, a non-default bijection). Doing that in one place is what
keeps the classes themselves short enough to read as physics.

``ampere.backends.reference`` has a helper of the same shape, and this is
deliberately not an import of it. That module's own docstring makes the
argument for why a helper like this is a *backend* concern rather than a shared
one, and importing another backend's private module would make the torch
package fail to import the moment the reference backend reorganised something
it never promised about. The extra cost is forty lines; the alternative is a
coupling between two backends that ``architecture.md`` §1 says are rungs of a
ladder, not collaborators.

A ``torch.Tensor`` is accepted where a number is, because a caller working in
tensors should not have to remember to convert one back before fixing a
parameter with it. It is converted immediately: a declaration is neutral data
(``lowering.md`` §1.4 and §5), so nothing tensor-shaped survives into the
:class:`~ampere.core.Parameter` itself.
"""

from __future__ import annotations

from typing import Any

import astropy.units as u
import numpy as np
import torch

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
        A frozen ``scipy.stats`` distribution (fitted), a number or scalar
        tensor (held fixed), or a ``Parameter`` (used as-is).
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
    if isinstance(given, torch.Tensor):
        if given.numel() != 1:
            raise ParameterError(
                f"parameter {name!r} was given a tensor of {given.numel()} element(s); a fixed "
                f"value must be a scalar."
            )
        return Parameter(name, value=float(given.detach().reshape(())), fixed=True, unit=unit)
    if isinstance(given, (int, float, np.floating, np.integer)) and not isinstance(given, bool):
        return Parameter(name, value=float(given), fixed=True, unit=unit)
    if hasattr(given, "ppf"):
        return Parameter(name, given, unit=unit)
    raise ParameterError(
        f"parameter {name!r} must be a frozen scipy.stats distribution (to fit it), a number "
        f"(to hold it fixed), or an ampere Parameter — got {type(given).__name__}."
    )
