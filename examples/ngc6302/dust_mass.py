"""Dust masses from a fit's posterior -- the v2 continuation of
``examples/NGC6302-calculate-dust-mass.py`` (W6.13 (2), ruling 4).

Equations 4 and 5 of Kemper et al. (2002, A&A 394, 679), transcribed
exactly. Legacy's own script hard-codes a *single* previous fit's ``n0``
(abundance), ``Tin`` and ``Tout`` per species and component and prints one
number per species; here those three become the posterior draws --
:func:`dust_masses_at` takes one point (a name-to-value mapping, in the
model's own unqualified parameter names), and :func:`dust_masses` takes a
run's ``DataTree`` and reports quantiles over every draw. Both call the same
per-point equations (:func:`_component_masses`).

The overall size normalisation, ``r0``
---------------------------------------
Legacy's script anchors ``r0`` (the shell's characteristic radius, equation
4) to Kemper et al. (2002)'s own *published* cold and warm amorphous-olivine
solution -- a fixed mass, ``Tin``, ``Tout`` and ``n0`` for each component,
external to this model and never fitted (``examples/NGC6302-calculate-dust-mass.py``
lines 96-121, reproduced as :data:`COLD_ANCHOR`/:data:`WARM_ANCHOR`) -- then
rescales it to *this* fit's own ``Tin`` via the ``Tin ** (-1/q)`` factor
equation 4 gives. Only that rescaling step uses a fit-dependent quantity;
the anchor itself stays fixed regardless of which draw is being scored, for
both :func:`dust_masses_at` and :func:`dust_masses`.

The species correspondence
---------------------------
This script's own ``names``/``rhod`` arrays (lines 17-19) list the eight
species in a *different* order from :mod:`examples.ngc6302.ngc6302`'s model
index (0 = calcite, ..., 7 = olivine); :data:`COLD_PARAMETER` and
:data:`WARM_PARAMETER` state the correspondence explicitly (also documented,
in the other direction, in :mod:`examples.ngc6302.generators`'s module
docstring) rather than relying on the two scripts' index orders lining up.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numpy as np

__all__ = [
    "COLD_ANCHOR",
    "COLD_PARAMETER",
    "RHOD",
    "SPECIES_NAMES",
    "WARM_ANCHOR",
    "WARM_PARAMETER",
    "dust_masses",
    "dust_masses_at",
    "format_table",
]

#: Kemper et al. (2002)'s own species order and grain densities, g cm^-3
#: (``examples/NGC6302-calculate-dust-mass.py`` lines 16-19).
SPECIES_NAMES = ("am. oliv.", "forst.", "calcite", "ice", "diopside", "c-enst.", "dolomite", "iron")
RHOD = (3.71, 3.33, 2.71, 1.00, 3.4, 2.80, 2.84, 7.874)

#: Which of :mod:`.ngc6302`'s fifteen fitted parameters supplies each
#: species' n0 in each component -- ``None`` where legacy's own model never
#: fits that species/component combination (``examples/NGC6302.py`` lines
#: 158-178's ``acold``/``awarm`` construction).
COLD_PARAMETER: dict[str, str | None] = {
    "am. oliv.": "logacold7",
    "forst.": "logacold2",
    "calcite": "logacold0",
    "ice": "logacold6",
    "diopside": "logacold3",
    "c-enst.": "logacold1",
    "dolomite": "logacold4",
    "iron": None,
}
WARM_PARAMETER: dict[str, str | None] = {
    "am. oliv.": "logawarm7",
    "forst.": "logawarm2",
    "calcite": None,
    "ice": None,
    "diopside": None,
    "c-enst.": "logawarm1",
    "dolomite": None,
    "iron": "logawarm5",
}

# Fixed constants, ``examples/NGC6302-calculate-dust-mass.py`` lines 91-94.
_P_INDEX = 0.5
_Q_INDEX = 0.5
_GRAIN_RADIUS_CM = 0.1e-4
_MSUN_G = 1.9885e33
_RHOD_ANCHOR = 3.71  # amorphous olivine, both components.

#: The 2002 cold/warm amorphous-olivine anchor (mass in Msun, Tin/Tout in K,
#: n0 dimensionless) -- ``examples/NGC6302-calculate-dust-mass.py`` lines
#: 100-121. Fixed; never drawn from a posterior.
COLD_ANCHOR: dict[str, float] = {"mass": 4.7e-2, "tin": 60.0, "tout": 30.0, "n0": 3.9e-2}
WARM_ANCHOR: dict[str, float] = {"mass": 6.1e-6, "tin": 118.0, "tout": 100.0, "n0": 1.2e-4}


def _r0_anchor(anchor: Mapping[str, float]) -> float:
    """``r0`` implied by the 2002 anchor alone, before rescaling to a fit's own Tin."""
    r0_cubed = (anchor["mass"] * _MSUN_G) / (
        math.pi
        / (3.0 - _P_INDEX)
        * ((anchor["tout"] / anchor["tin"]) ** (-4.0) - 1.0)
        * (4.0 / 3.0)
        * math.pi
        * _GRAIN_RADIUS_CM**3
        * _RHOD_ANCHOR
        * anchor["n0"]
    )
    return r0_cubed ** (1.0 / 3.0)


def _r0(anchor: Mapping[str, float], tin: np.ndarray) -> np.ndarray:
    """Equation (4): rescale the anchor's own ``r0`` to a fit's own ``Tin``."""
    return _r0_anchor(anchor) * (tin / anchor["tin"]) ** (-1.0 / _Q_INDEX)


def _masses(tin: np.ndarray, tout: np.ndarray, r0: np.ndarray, n0: Any, rhod: float) -> np.ndarray:
    """Equation (5), integrated over the shell -- solar masses."""
    return (
        (math.pi * r0**3)
        / (3.0 - _P_INDEX)
        * ((tout / tin) ** (-4.0) - 1.0)
        * (4.0 / 3.0)
        * math.pi
        * _GRAIN_RADIUS_CM**3
        * rhod
        * n0
    ) / _MSUN_G


def _n0(theta: Mapping[str, Any], parameter: str | None) -> Any:
    if parameter is None:
        return 0.0
    return 10.0 ** np.asarray(theta[parameter], dtype=float)


def _component_masses(theta: Mapping[str, Any]) -> dict[str, dict[str, np.ndarray]]:
    """Per-species mass in each component, at *theta* (scalars or draw arrays)."""
    tcold_in = np.asarray(theta["Tcold1"], dtype=float)
    tcold_out = np.asarray(theta["Tcold0"], dtype=float)
    twarm_in = np.asarray(theta["Twarm1"], dtype=float)
    twarm_out = np.asarray(theta["Twarm0"], dtype=float)
    r0_cold = _r0(COLD_ANCHOR, tcold_in)
    r0_warm = _r0(WARM_ANCHOR, twarm_in)

    cold: dict[str, np.ndarray] = {}
    warm: dict[str, np.ndarray] = {}
    for name, rhod in zip(SPECIES_NAMES, RHOD, strict=True):
        cold[name] = _masses(tcold_in, tcold_out, r0_cold, _n0(theta, COLD_PARAMETER[name]), rhod)
        warm[name] = _masses(twarm_in, twarm_out, r0_warm, _n0(theta, WARM_PARAMETER[name]), rhod)
    return {"cold": cold, "warm": warm}


def dust_masses_at(theta: Mapping[str, float]) -> dict[str, Any]:
    """Dust mass per species and component at one parameter point, solar masses.

    *theta* is a name-to-value mapping using :class:`~examples.ngc6302.ngc6302.KemperTwoShell`'s
    own (unqualified) parameter names. At the legacy script's own hard-coded
    inputs (:data:`examples.ngc6302.generators.TRUTH`) this reproduces
    ``examples/NGC6302-calculate-dust-mass.py``'s printed ``mdust`` numbers
    -- checked in the smoke test to ``rtol=1e-6``.
    """
    masses = _component_masses(theta)
    return {
        "cold": {name: float(value) for name, value in masses["cold"].items()},
        "warm": {name: float(value) for name, value in masses["warm"].items()},
        "totals": {
            component: float(sum(values.values())) for component, values in masses.items()
        },
    }


def dust_masses(tree: Any, *, quantiles: tuple[float, ...] = (0.16, 0.5, 0.84)) -> dict[str, Any]:
    """Dust mass per species and component, over a run's posterior draws.

    *tree* is a fit's ``DataTree`` (as :func:`examples.ngc6302.ngc6302.fit`
    returns); *quantiles* are taken over the flattened chain x draw
    posterior. Same equations as :func:`dust_masses_at`, vectorised over
    draws -- legacy's hard-coded ``n0``/``Tin``/``Tout`` become the
    posterior draws (ruling 4).
    """
    posterior = tree["posterior"].dataset
    names = {"Tcold0", "Tcold1", "Twarm0", "Twarm1"} | {
        parameter
        for parameter in (*COLD_PARAMETER.values(), *WARM_PARAMETER.values())
        if parameter is not None
    }
    theta = {name: np.asarray(posterior[f"model.{name}"], dtype=float).ravel() for name in names}

    masses = _component_masses(theta)
    table: dict[str, Any] = {"cold": {}, "warm": {}}
    for component in ("cold", "warm"):
        for name, draws in masses[component].items():
            table[component][name] = {
                f"q{quantile:.2f}": float(np.quantile(draws, quantile)) for quantile in quantiles
            }
    table["totals"] = {
        component: {
            f"q{quantile:.2f}": float(np.quantile(sum(masses[component].values()), quantile))
            for quantile in quantiles
        }
        for component in ("cold", "warm")
    }
    return table


def format_table(table: Mapping[str, Any]) -> str:
    """A human-readable rendering of :func:`dust_masses_at`'s or :func:`dust_masses`' table."""
    lines = ["dust mass (Msun), by species and component:"]

    def _render(value: Any) -> str:
        if isinstance(value, Mapping):
            return ", ".join(f"{key}={val:.4g}" for key, val in value.items())
        return f"{float(value):.4g}"

    for component in ("cold", "warm"):
        lines.append(f"  {component}:")
        for name, value in table[component].items():
            lines.append(f"    {name:12s} {_render(value)}")
    totals = table.get("totals")
    if totals:
        for component, value in totals.items():
            lines.append(f"  total {component:6s} {_render(value)}")
    return "\n".join(lines)
