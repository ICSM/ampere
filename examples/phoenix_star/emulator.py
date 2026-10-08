"""The PHOENIX emulator as a v2 :class:`~ampere.core.Model` -- W6.13 (4).

``phoenix_emulator.npz`` (beside this file, written once by
:mod:`.train_emulator`) holds a PCA of ``log10`` unit-bolometric PHOENIX-ACES
spectra and one Gaussian process per PCA weight over (Teff, log g, [Fe/H]).
This module loads it (:func:`load_emulator`) and evaluates it
(:func:`log_shape`): the GP predictive mean ``k(x, X) alpha`` of every weight at
the standardised input, then ``mean + weights @ basis``. That arithmetic is
written once, against an array namespace (:class:`ampere.core.ArrayOps`, the
public protocol since W7.13) -- :class:`~ampere.core.NumpyOps` here,
:class:`~ampere.backends.torch.TorchOps` in :mod:`.emulator_torch`,
:class:`~ampere.backends.jax.JaxOps` in :mod:`.emulator_jax` --
so the three backends run literally the same sequence of ``exp``, ``sum`` and
``matmul``, and cannot drift. It uses only operations all three spell alike.

What lives here, shared by :mod:`.phoenix_star` and :mod:`examples.star_disc`:

* :class:`PhoenixEmulator`, the bare emulator as a reference-backend model:
  parameters ``teff`` (K), ``logg`` and ``feh``; two channels, ``"sed"`` on
  segment A (``np.geomspace(0.3, 5.5, 582)`` micron) and ``"rvs"`` on segment
  B (the legacy Gaia-RVS grid, 387 points at R = 11 000), each the
  unit-bolometric ``F_lambda`` shape, per angstrom.
* :func:`channel_plan` and :func:`star_log_flux`: the star -- the emulator's
  shape turned into ``F_nu`` in Jy for a luminosity and a distance, placed on
  an arbitrary grid by linear interpolation in ``log F_nu``-``log lambda``
  (the legacy ``interp1d`` on log-log axes), extrapolated beyond 5.5 micron
  as ``F_nu ~ lambda**-2`` from the last emulated pixel (the legacy's own
  choice; ``QuickSED``'s ``F_lambda ~ lambda**-4`` is the same law).
* :func:`ccm89_ab`: Cardelli, Clayton & Mathis (1989) ``a(x)`` and ``b(x)``
  (their equations 2a-3b) -- the infrared power law for ``x < 1.1`` per
  micron and the optical/near-infrared polynomial for ``1.1 <= x <= 3.3``;
  ``A_lambda / A_V = a(x) + b(x) / R_V``. The infrared branch is used below
  its nominal ``x = 0.3`` too (the legacy model grid reaches 10 micron,
  ``x = 0.1``), where the extinction is under 2 % of ``A_V`` anyway. These
  depend on the grid only, so each channel precomputes them and the model
  applies ``10 ** (-0.4 * a_v * (a + b / r_v))`` in backend arithmetic.

The inverse-square law
----------------------
``F_nu(Jy) = shape x f_bol(1 L_sun, 1 pc) x 10**log_luminosity / d_pc**2 x
3.34e5 lambda_A**2``, where ``f_bol(1 L_sun, 1 pc) = L_sun / (4 pi pc**2)``
already carries the ``4 pi``. The legacy scripts (``phoenixstar.py`` line 53,
``QuickSED.py`` line 148) divide by ``4 pi d**2`` a *second* time, making their
fluxes ``4 pi`` too faint; the twins do not reproduce that.

Where the emulator errs
-----------------------
The GP smooths rather than interpolates (see :mod:`.train_emulator`, "The
jitter"), so it misses the grid by a few per cent at its worst pixels. The
three maxima the training reports all sit at the cool edge of the grid and the
blue end of segment A, outside every fitted passband:

* at the nodes (check (ii)), 7.9 % against the PCA reconstruction (7.8 %
  against the true spectrum) at 0.3015 micron, node Teff 5400 K, log g 4.5,
  [Fe/H] +0.5;
* leave-one-out on segment A (check (i)), 14 % at 0.300 micron, node 5000 K,
  4.0, 0.0 -- a corner of the grid, so an extrapolation when held out;
* leave-one-out on segment B, 2.0 % at 0.8668 micron, the same node -- inside
  the RVS window and the Gaia G and RP passbands.

The eight filters of :mod:`.phoenix_star` start at 0.328 micron (Gaia BP, at
1 % of peak), so the two segment-A maxima lie blueward of every passband and
outside the RVS window (0.842-0.872 micron). Both nodes lie below the twins'
Teff priors (6000 K for :mod:`.phoenix_star`, 5500 K for
:mod:`examples.star_disc`). Errors above 1 % are not confined there, though:
over all nodes, 297 of the 371 segment-A pixels that exceed 1 % somewhere
against the true spectrum lie inside a passband, and in segment B every such
pixel is inside the RVS window (41 pixels at the nodes, 24 leave-one-out, at
most 2.5 % at a node). The RMS errors are 0.38 % (A) and 0.21 % (B) at the
nodes, and 0.61 % and 0.29 % leave-one-out.
"""

from __future__ import annotations

import functools
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, ClassVar

import astropy.units as u
import numpy as np
import scipy.stats as st
from astropy import constants

from ampere.core import DTYPE, ArrayOps, Model, ModelResult, NumpyOps, Parameter, Spectrum

__all__ = [
    "EMULATOR_FILE",
    "FBOL_1L1P",
    "FNU_FACTOR",
    "NUMPY_OPS",
    "ChannelPlan",
    "EmulatorArrays",
    "PhoenixEmulator",
    "as_parameter",
    "ccm89_ab",
    "channel_plan",
    "load_emulator",
    "log_shape",
    "star_log_flux",
]

EMULATOR_FILE = Path(__file__).resolve().parent / "phoenix_emulator.npz"

#: Bolometric flux of one solar luminosity at one parsec, erg s^-1 cm^-2 --
#: the legacy ``fbol_1l1p`` (``phoenixstar.py`` lines 23-26).
FBOL_1L1P = float((constants.L_sun / (4 * np.pi * (1 * u.pc) ** 2)).to(u.erg / u.s / u.cm**2).value)

#: ``F_lambda`` (erg s^-1 cm^-2 A^-1) to ``F_nu`` (Jy) is ``x 3.34e5 lambda_A**2``,
#: the legacy's own constant.
FNU_FACTOR = 3.34e5

#: The channel each emulator segment serves.
SEGMENTS = ("sed", "rvs")


@dataclass(frozen=True)
class EmulatorArrays:
    """The committed file's arrays (numpy), and its provenance record."""

    segment_a: np.ndarray
    segment_b: np.ndarray
    mean: np.ndarray
    basis: np.ndarray
    nodes: np.ndarray
    input_centre: np.ndarray
    input_spread: np.ndarray
    alpha: np.ndarray
    hyper: np.ndarray
    grid_parameters: np.ndarray
    stored_index: np.ndarray
    stored_spectra: np.ndarray
    provenance: dict[str, Any]

    @property
    def components(self) -> int:
        return int(self.basis.shape[0])

    def segment(self, channel: str) -> tuple[np.ndarray, slice]:
        """The wavelength grid (micron) and the slice of the emulated vector for *channel*."""
        n_a = self.segment_a.size
        if channel == "sed":
            return self.segment_a, slice(0, n_a)
        if channel == "rvs":
            return self.segment_b, slice(n_a, n_a + self.segment_b.size)
        raise KeyError(f"the PHOENIX emulator has channels {SEGMENTS}, not {channel!r}")


@functools.lru_cache(maxsize=4)
def load_emulator(path: str | Path = EMULATOR_FILE) -> EmulatorArrays:
    """Read ``phoenix_emulator.npz`` (cached per path)."""
    with np.load(path) as data:
        arrays = {name: np.asarray(data[name]) for name in data.files if name != "provenance"}
        provenance = json.loads(str(data["provenance"]))
    return EmulatorArrays(provenance=provenance, **arrays)


#: The reference namespace: the public :class:`~ampere.core.NumpyOps` (W7.13
#: folded this example's own ``Ops`` class onto it). ``asarray`` makes a
#: float64 array, ``asindex`` an integer one for gathers; ``exp``, ``log10``,
#: ``sum`` and ``stack`` are the arithmetic.
NUMPY_OPS: ArrayOps = NumpyOps()


def log_shape(ops: ArrayOps, arrays: dict[str, Any], teff: Any, logg: Any, feh: Any) -> Any:
    """``log10`` of the unit-bolometric shape on segments A and B concatenated.

    *arrays* holds ``mean``, ``basis``, ``nodes``, ``alpha``, ``lengths``,
    ``amplitude``, ``centre`` and ``spread`` already converted by *ops*.
    """
    x = (
        ops.stack([ops.asarray(teff), ops.asarray(logg), ops.asarray(feh)], axis=0)
        - arrays["centre"]
    ) / arrays["spread"]
    # (K, m, 3): every component's scaled distance to every node.
    d = (x - arrays["nodes"])[None, :, :] / arrays["lengths"][:, None, :]
    k = arrays["amplitude"][:, None] ** 2 * ops.exp(-0.5 * ops.sum(d**2, -1))
    weights = ops.sum(k * arrays["alpha_t"], -1)
    return arrays["mean"] + weights @ arrays["basis"]


def emulator_arrays(ops: ArrayOps, data: EmulatorArrays) -> dict[str, Any]:
    """*data*'s arrays in *ops*' kind, in the layout :func:`log_shape` reads."""
    return {
        "mean": ops.asarray(data.mean),
        "basis": ops.asarray(data.basis),
        "nodes": ops.asarray(data.nodes),
        "alpha_t": ops.asarray(data.alpha.T),
        "lengths": ops.asarray(data.hyper[:, :3]),
        "amplitude": ops.asarray(data.hyper[:, 3]),
        "centre": ops.asarray(data.input_centre),
        "spread": ops.asarray(data.input_spread),
    }


def ccm89_ab(wavelength: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """CCM89's ``a(x)`` and ``b(x)`` at *wavelength* (micron), ``x = 1/lambda``.

    Infrared (eqs 2a, 2b) below ``x = 1.1``; optical/near-infrared (eqs 3a,
    3b, ``y = x - 1.82``) for ``1.1 <= x <= 3.3``. Beyond ``x = 3.3`` (the
    ultraviolet) is refused: no grid here reaches it.
    """
    x = 1.0 / np.asarray(wavelength, dtype=float)
    if np.any(x > 3.3):
        raise ValueError("ccm89_ab covers x <= 3.3 per micron (lambda >= 0.303 micron) only")
    a = np.empty_like(x)
    b = np.empty_like(x)
    ir = x < 1.1
    a[ir] = 0.574 * x[ir] ** 1.61
    b[ir] = -0.527 * x[ir] ** 1.61
    y = x[~ir] - 1.82
    a[~ir] = (
        1
        + 0.17699 * y
        - 0.50447 * y**2
        - 0.02427 * y**3
        + 0.72085 * y**4
        + 0.01979 * y**5
        - 0.77530 * y**6
        + 0.32999 * y**7
    )
    b[~ir] = (
        1.41338 * y
        + 2.28305 * y**2
        + 1.07233 * y**3
        - 5.38434 * y**4
        - 0.62251 * y**5
        + 5.30260 * y**6
        - 2.09002 * y**7
    )
    return a, b


@dataclass(frozen=True)
class ChannelPlan:
    """How one channel's grid reads the emulator: gathers, weights, offsets.

    ``log F_nu(target) = w0 * S[i0] + w1 * S[i1] + offset`` where ``S`` is
    ``log10 F_nu`` of the unit-luminosity star at one parsec on the source
    segment (``log_shape`` plus :attr:`source_offset`). Interpolation rows have
    ``offset = 0``; the extrapolation rows beyond the segment's red end have
    ``i0 = i1 = last``, ``w0 = 1``, ``w1 = 0`` and ``offset = -2 log10(lambda /
    lambda_last)``. :attr:`ccm_a` and :attr:`ccm_b` are CCM89's ``a(x)`` and
    ``b(x)`` on the target grid (None when the model applies no extinction).
    """

    wavelength: np.ndarray
    segment: slice
    source_offset: np.ndarray
    i0: np.ndarray
    i1: np.ndarray
    w0: np.ndarray
    w1: np.ndarray
    offset: np.ndarray
    ccm_a: np.ndarray | None
    ccm_b: np.ndarray | None


def channel_plan(
    data: EmulatorArrays,
    channel: str,
    wavelength: Any,
    *,
    extinction: bool = True,
    extrapolate: bool | None = None,
) -> ChannelPlan:
    """The :class:`ChannelPlan` for *channel*'s *wavelength* grid (micron).

    The ``"sed"`` channel extrapolates beyond segment A's red end by default;
    ``"rvs"`` does not (its grid must lie inside segment B). Anything bluer
    than the segment's first pixel is refused.
    """
    target = np.asarray(wavelength, dtype=float)
    source, segment = data.segment(channel)
    extrapolate = channel == "sed" if extrapolate is None else extrapolate
    if target.min() < source[0] * (1 - 1e-12):
        raise ValueError(
            f"channel {channel!r}: the grid starts at {target.min():.6g} micron, bluer than "
            f"the emulator's {source[0]:.6g}."
        )
    beyond = target > source[-1]
    if np.any(beyond) and not extrapolate:
        raise ValueError(
            f"channel {channel!r}: the grid reaches {target.max():.6g} micron, beyond the "
            f"emulator's {source[-1]:.6g}, and this channel does not extrapolate."
        )
    log_source = np.log10(source)
    log_target = np.log10(np.minimum(target, source[-1]))
    i0 = np.clip(np.searchsorted(log_source, log_target, side="right") - 1, 0, source.size - 2)
    i1 = i0 + 1
    w1 = (log_target - log_source[i0]) / (log_source[i1] - log_source[i0])
    w0 = 1.0 - w1
    offset = np.zeros_like(target)
    last = source.size - 1
    i0 = np.where(beyond, last, i0)
    i1 = np.where(beyond, last, i1)
    w0 = np.where(beyond, 1.0, w0)
    w1 = np.where(beyond, 0.0, w1)
    offset = np.where(beyond, -2.0 * np.log10(target / source[-1]), offset)
    source_offset = np.log10(FNU_FACTOR * (source * 1e4) ** 2 * FBOL_1L1P)
    a, b = ccm89_ab(target) if extinction else (None, None)
    return ChannelPlan(target, segment, source_offset, i0, i1, w0, w1, offset, a, b)


def plan_arrays(ops: ArrayOps, plan: ChannelPlan) -> dict[str, Any]:
    """*plan*'s arrays in *ops*' kind."""
    out = {
        "source_offset": ops.asarray(plan.source_offset),
        "i0": ops.asindex(plan.i0),
        "i1": ops.asindex(plan.i1),
        "w0": ops.asarray(plan.w0),
        "w1": ops.asarray(plan.w1),
        "offset": ops.asarray(plan.offset),
    }
    if plan.ccm_a is not None and plan.ccm_b is not None:
        out["ccm_a"] = ops.asarray(plan.ccm_a)
        out["ccm_b"] = ops.asarray(plan.ccm_b)
    return out


def star_log_flux(
    ops: ArrayOps,
    plan: ChannelPlan,
    arrays: dict[str, Any],
    shape: Any,
    log_luminosity: Any,
    distance: Any,
    a_v: Any = None,
    r_v: Any = None,
) -> Any:
    """``log10 F_nu`` (Jy) of the star on *plan*'s grid.

    *shape* is :func:`log_shape`'s output; *arrays* :func:`plan_arrays`'.
    ``log_luminosity`` is ``log10(L / L_sun)``, ``distance`` in parsec. The
    extinction term is added when *plan* carries CCM89 coefficients and *a_v*
    is given.
    """
    source = shape[plan.segment] + arrays["source_offset"]
    log_flux = arrays["w0"] * source[arrays["i0"]] + arrays["w1"] * source[arrays["i1"]]
    log_flux = log_flux + arrays["offset"]
    log_flux = log_flux + ops.asarray(log_luminosity) - 2.0 * ops.log10(ops.asarray(distance))
    if a_v is not None and "ccm_a" in arrays:
        extinction = ops.asarray(a_v) * (arrays["ccm_a"] + arrays["ccm_b"] / ops.asarray(r_v))
        log_flux = log_flux - 0.4 * extinction
    return log_flux


def as_parameter(name: str, spec: Any) -> Parameter:
    """A prior (a frozen scipy distribution) or a bare number, as a Parameter."""
    if isinstance(spec, Parameter):
        return spec
    if hasattr(spec, "logpdf") or hasattr(spec, "rvs"):
        return Parameter(name, spec)
    return Parameter(name, None, value=float(spec), fixed=True)


#: The emulator's own box: Teff 5000-8000 K, log g 4-5, [Fe/H] 0-0.5.
DEFAULT_PRIORS: dict[str, Any] = {
    "teff": st.uniform(5000.0, 3000.0),
    "logg": st.uniform(4.0, 1.0),
    "feh": st.uniform(0.0, 0.5),
}


class PhoenixEmulator(Model):
    """The bare emulator: ``teff``, ``logg``, ``feh`` -> two unit-bolometric shapes.

    Channels ``"sed"`` (segment A) and ``"rvs"`` (segment B), each ``F_lambda``
    normalised to unit bolometric flux, per angstrom. :attr:`OPS` is the array
    namespace; the torch and jax twins subclass this and change only that and
    the capability flags.
    """

    OPS: ClassVar[ArrayOps] = NUMPY_OPS
    UNIT: ClassVar[Any] = u.AA**-1

    def __init__(
        self,
        *,
        teff: Any = DEFAULT_PRIORS["teff"],
        logg: Any = DEFAULT_PRIORS["logg"],
        feh: Any = DEFAULT_PRIORS["feh"],
        path: str | Path = EMULATOR_FILE,
    ) -> None:
        self.data = load_emulator(path)
        for name, given in {"teff": teff, "logg": logg, "feh": feh}.items():
            self.register_parameter(as_parameter(name, given))
        self.arrays = emulator_arrays(self.OPS, self.data)
        self.grids = {
            channel: self.OPS.asarray(self.data.segment(channel)[0]) for channel in SEGMENTS
        }
        self.templates = {
            channel: Spectrum(
                self.data.segment(channel)[0] * u.um,
                np.zeros(self.data.segment(channel)[0].size, dtype=DTYPE),
                unit=self.UNIT,
            )
            for channel in SEGMENTS
        }

    def grid(self, channel: str) -> Any:
        """*channel*'s wavelength grid, micron, as a backend array."""
        return self.grids[channel]

    def flux(self, channel: str, values: Any = None) -> Any:
        """*channel*'s unit-bolometric shape as a backend array."""
        ctx = self.context(values)
        full = log_shape(self.OPS, self.arrays, ctx["teff"], ctx["logg"], ctx["feh"])
        return 10.0 ** full[self.data.segment(channel)[1]]

    def evaluate(self, **values: Any) -> ModelResult:
        emitted = {}
        for channel in SEGMENTS:
            flux = np.asarray(self.OPS.to_numpy(self.flux(channel, values)), dtype=DTYPE)
            emitted[channel] = self.templates[channel].with_values(flux)
        return ModelResult(emitted)
