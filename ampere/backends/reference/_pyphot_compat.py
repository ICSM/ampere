"""The v2 home of the pyphot unit-API compatibility layer.

This is what ``ampere/legacy/utils/pyphot_compat.py`` provides, lifted so that
v2 imports no legacy module (W7.5, W6.0's carried item); the legacy module
stays untouched, as legacy is frozen, and ``ampere/legacy/data/photometry.py``
keeps importing it.

pyphot 1.x exposed a module-level unit registry, ``pyphot.unit``; pyphot >= 2
removed it in favour of a pluggable adapter at ``pyphot.config.units``, whose
``.U(name)`` returns an ``astropy.units.Unit`` when astropy is importable.
:func:`get_unit` hides the difference: write ``value * get_unit("micron")`` and
both major versions work.
"""

from __future__ import annotations

from typing import Any

import pyphot

#: True once the installed pyphot has dropped the legacy pint-based
#: ``pyphot.unit`` registry. Prefer :func:`get_unit` to testing this directly.
PYPHOT_V2: bool = not hasattr(pyphot, "unit")


def get_unit(name: str) -> Any:
    """A unit object for *name*, in whichever idiom the installed pyphot uses.

    Multiply it onto an array to attach units before passing the result to a
    pyphot ``Filter`` method, or give it to a quantity's ``.to()``.

    Parameters
    ----------
    name : str
        A unit name pyphot recognises, such as ``"micron"``, ``"AA"``,
        ``"flam"``, ``"fnu"`` or ``"Jy"``.

    Returns
    -------
    unit : astropy.units.Unit or pint unit
        ``pyphot.config.units.U(name)`` under pyphot >= 2, ``pyphot.unit[name]``
        under pyphot 1.x.
    """
    if PYPHOT_V2:
        return pyphot.config.units.U(name)
    return pyphot.unit[name]  # type: ignore[attr-defined]
