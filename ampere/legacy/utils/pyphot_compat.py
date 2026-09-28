"""Compatibility layer for pyphot's unit-handling API across the pyphot
1.x -> 2.x transition (W0.9).

**Background.** pyphot 1.x exposed a module-level, pint-backed unit
registry, ``pyphot.unit``: quantities usable with ``Filter.get_flux`` /
``Filter.lpivot`` etc. were built as ``value * pyphot.unit['micron']``.
pyphot >= 2 removed ``pyphot.unit`` entirely as part of a rework that makes
the unit system pluggable ("unit adapters" -- astropy, pint, or pyphot's own
legacy pint wrapper). The active adapter lives at ``pyphot.config.units``
and, when astropy is importable (always true here -- it is a core ampere
dependency), **defaults to an astropy-units-backed adapter**: its ``.U(name)``
returns a plain ``astropy.units.Unit``, and quantities are built exactly as
``value * pyphot.config.units.U('micron')``. This is pyphot's own documented
idiom for >= 2 (see its ``QuickStart`` notebook), not an ampere invention.

**Primary surface (pyphot >= 2 idiom) -- use this in all new code,
including Phase 2's synthetic-photometry Transformation:**

    from ampere.utils.pyphot_compat import get_unit

    wave = model_wavelength_um * get_unit("micron")
    flux = model_flux_flam * get_unit("flam")
    synthetic = a_filter.get_flux(wave, flux).value

``get_unit(name)`` returns a unit object that is a drop-in replacement for
the old ``pyphot.unit[name]`` in every place ampere used it (multiplication
to build a quantity, and ``Quantity.to(unit)`` conversions -- both astropy
and pint quantities accept a plain unit object, or a unit name string,
interchangeably there). Under pyphot >= 2 it *is* the astropy idiom
directly; under pyphot 1.x it falls back to the legacy ``pyphot.unit``
registry so old environments keep working during the transition. Callers
never need to branch on the installed pyphot version themselves -- write
against ``get_unit`` and both major versions work.

Phase 2 code should only ever import ``get_unit`` (or, if it ever needs to
branch explicitly, the ``PYPHOT_V2`` flag below) from this module -- never
reach for ``pyphot.unit`` directly -- so it targets the >= 2 API from its
first line rather than inheriting the 1.x idiom this module exists to
retire.
"""

from __future__ import annotations

import pyphot

#: True once the installed pyphot has moved to the >= 2 unit-adapter
#: rework (i.e. no longer exposes the legacy pint-based ``pyphot.unit``
#: registry). Exposed for the rare case calling code needs to branch
#: explicitly; prefer ``get_unit`` over testing this directly.
PYPHOT_V2 = not hasattr(pyphot, "unit")


def get_unit(name: str):
    """Return a unit object for ``name``, in whichever unit-handling idiom
    the installed pyphot understands.

    The result is meant to be used the same way ``pyphot.unit[name]`` used
    to be under pyphot 1.x: multiply it onto a plain array/scalar to attach
    units (``value * get_unit('micron')``) before passing the result to a
    pyphot ``Filter`` method, or pass it to a quantity's ``.to()`` for unit
    conversion.

    Parameters
    ----------
    name : str
        A unit name pyphot recognises, e.g. ``"micron"``, ``"AA"``,
        ``"flam"``, ``"fnu"``, ``"Jy"``.

    Returns
    -------
    unit : astropy.units.Unit or pint unit
        Under pyphot >= 2 (``PYPHOT_V2`` is ``True``), this is
        ``pyphot.config.units.U(name)`` -- an ``astropy.units.Unit`` for the
        default (astropy) unit adapter, which is what pyphot itself uses
        when astropy is installed. Under pyphot 1.x, this is
        ``pyphot.unit[name]`` from the legacy pint-based registry.
    """
    if PYPHOT_V2:
        return pyphot.config.units.U(name)
    return pyphot.unit[name]
