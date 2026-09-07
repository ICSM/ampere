"""The standard instrument steps on the torch path.

``DEVELOPMENT_PLAN.md`` §5's Phase 2 list — "resampling, LSF, synthetic
photometry, calibration factor" — with the shapes ``transformations.md`` §10's
table fixes, computed in ``torch``. Same declarations, same buffer names, same
published requirements as :mod:`ampere.backends.reference.instrument`: the two
are held to each other by the conformance battery, so a differently named
buffer or a differently shaped requirement is a failure rather than a detail.

Where the two backends genuinely differ
----------------------------------------
Every step here is a **linear operator applied to the flux**, and that
factorisation is what makes a torch version worth having:

* building the operator — bin edges, overlap integrals, a Gaussian kernel, a
  response integral — happens **once**, at construction or at
  ``configure_from`` time, from constant data;
* applying it happens **per evaluation**, and is one ``matmul``.

So the operator is built in torch and registered as a **torch buffer**
(``lowering.md`` §7: traced, moved by ``.to()``, never differentiated), and the
per-evaluation application is ``weights @ values`` on float64 tensors. That is
a real second implementation, not a rename: the arithmetic is torch's, and the
conformance battery's cross-backend rows compare it against numpy's.

The mask is not part of that. A mask is boolean bookkeeping, not arithmetic —
``propagate_mask``'s ANY rule is a contract about *which* samples an output
depends on, identical on every backend — so it is computed with
``ampere.core``'s own function on the numpy view of the same weights, rather
than being reimplemented in torch for the sake of it.

Two rules from the contracts govern everything here, exactly as they do on the
reference path.

**Coordinates come from the container, never from arithmetic**
(``transformations.md`` §10, W1.11 gap I-2): a step reproducing the observed
coordinates takes them *from* the observed container, because
``check_alignment`` compares axes with ``np.array_equal`` and a recomputed grid
differs in its last bits.

**A step whose buffer is tabulated on particular coordinates publishes
``points=`` and looks them up** (``spectrum_photometry.md`` Gap 1):
:class:`SyntheticPhotometry` is this backend's instance of that shape, and
:meth:`~ampere.core.Axis.locate` is how it stays correct once a second
instrument binds the same channel.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, ClassVar

import astropy.units as u
import numpy as np
import torch

from ampere.core import (
    DTYPE,
    AxisRequirement,
    PhotometricPoints,
    Spectrum,
    Transformation,
    TransformationError,
    propagate_mask,
)

from ._config import BACKEND, DEFAULT_DEVICE, DEFAULT_DTYPE, as_tensor, to_numpy
from ._declare import as_parameter
from .parameters import LoweredParameters

__all__ = [
    "DETECTORS",
    "CalibrationScale",
    "LSFConvolution",
    "Resample",
    "SyntheticPhotometry",
    "TorchStep",
    "bin_edges",
    "bundled_filter_library",
]

COORDINATE_UNIT = u.micron

#: The two synthetic-photometry conventions, mapped to the power of lambda
#: dividing the weight — the same table, and the same ruling (Peter,
#: 2026-09-05), as the reference backend's. ``"photon"`` is a detector that
#: counts photons (weight ``R dlambda / lambda``); ``"energy"`` is one that
#: measures energy (weight ``R dlambda / lambda**2``). They are different
#: numbers, not different spellings, so there is no default and no third
#: option. See :mod:`ampere.backends.reference.instrument` for the derivation.
DETECTORS: dict[str, float] = {"photon": 1.0, "energy": 2.0}

#: How far beyond the output range an LSF needs coverage, in kernel sigmas.
_LSF_PADDING_SIGMAS = 5.0
#: How much finer than the kernel the input grid must be.
_LSF_OVERSAMPLING = 3.0


def bin_edges(centres: torch.Tensor) -> torch.Tensor:
    """Bin edges bracketing *centres*: midpoints, with the outer two reflected.

    The usual convention for a spectrum given as sample centres. A single
    sample has no midpoints to work with, so it is given a unit-width bin — the
    width then cancels out of every normalised weight that uses it.
    """
    if centres.numel() == 1:
        return torch.stack([centres[0] - 0.5, centres[0] + 0.5])
    middle = 0.5 * (centres[1:] + centres[:-1])
    first = centres[0] - (middle[0] - centres[0])
    last = centres[-1] + (centres[-1] - middle[-1])
    return torch.cat([first.reshape(1), middle, last.reshape(1)])


class TorchStep(Transformation):
    """Shared plumbing: the declared buffers, their tensors, and the four flags.

    A step's constant arrays live twice, and both are load-bearing (see
    :class:`~ampere.backends.torch.models.TorchSpectralModel` for the same
    split): once in ``ampere.core``'s :class:`~ampere.core.BufferSet`, which is
    the *declared* surface the contract, the provenance record and the model
    fingerprint read; and once as a torch buffer on :attr:`tensors`, which is
    what a device move, a ``state_dict`` and the per-evaluation ``matmul`` see.
    """

    #: The four capability flags (``inference.md`` §18, W2.12's fourth), stated
    #: rather than inherited. ``BATCHABLE`` is ``False`` honestly: nothing here
    #: accepts a stack of parameter vectors yet.
    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = False
    DEVICE: ClassVar[str] = "cpu"
    BACKEND: ClassVar[str] = BACKEND

    def __init__(
        self,
        *,
        label: str | None = None,
        dtype: torch.dtype = DEFAULT_DTYPE,
        device: torch.device = DEFAULT_DEVICE,
    ) -> None:
        super().__init__(label=label)
        self.dtype = dtype
        self.device = device
        self.tensors = LoweredParameters()

    def _declare(self, name: str, value: Any, *, unit: u.UnitBase | None = None) -> torch.Tensor:
        """Register one constant array as both a declared buffer and a torch buffer."""
        array = np.asarray(value, dtype=DTYPE)
        self.register_buffer(name, array, unit=unit)
        tensor = as_tensor(array, dtype=self.dtype, device=self.device)
        self.tensors.register_buffer(name, tensor, persistent=True)
        return tensor

    def _tensor(self, name: str) -> torch.Tensor:
        return self.tensors.get_buffer(name)

    def _data(self, name: str) -> np.ndarray:
        """The declared buffer as a plain array.

        A buffer's name may not shadow a class attribute, so a step whose
        buffer is called ``target`` cannot also expose a ``target`` property;
        this is the readable spelling of ``self.buffers[name].value``.
        """
        return np.asarray(self.buffers[name].value, dtype=DTYPE)


class CalibrationScale(TorchStep):
    """A multiplicative calibration factor: ``flux -> scale * flux``.

    The commonest instrumental nuisance there is — an uncertain absolute
    calibration, or the relative normalisation between two instruments
    observing the same source. Kind-preserving, and it publishes no
    requirements: scaling does not care what grid it is on.

    Parameters
    ----------
    scale
        A prior to fit the factor, or a number to hold it fixed. The natural
        prior is log-normal about 1, since a calibration factor is positive and
        multiplicative.
    label
        Component label for this step's parameters within a chain.
    """

    ACCEPTS: ClassVar[tuple[type, ...]] = (Spectrum,)

    def __init__(
        self,
        scale: Any = 1.0,
        *,
        label: str | None = None,
        dtype: torch.dtype = DEFAULT_DTYPE,
        device: torch.device = DEFAULT_DEVICE,
    ) -> None:
        super().__init__(label=label, dtype=dtype, device=device)
        self.register_parameter(as_parameter("scale", scale))

    def apply_tensor(self, values: torch.Tensor, factor: torch.Tensor) -> torch.Tensor:
        """The differentiable step: one multiply, gradient in both arguments."""
        return values * factor

    def apply(self, samples: Any, values: Any) -> Any:
        factor = as_tensor(self.context(values)["scale"], dtype=self.dtype, device=self.device)
        flux = as_tensor(samples.values, dtype=self.dtype, device=self.device)
        # with_values inherits the axes, the unit and the mask, so the
        # calibration costs one multiply and no revalidation.
        return samples.with_values(to_numpy(self.apply_tensor(flux, factor)))


class Resample(TorchStep):
    """Resample a spectrum onto a coarser grid, conserving mean flux density.

    Each output sample is the overlap-weighted mean of the input samples whose
    bins intersect its own — the standard binning integral, written as an
    explicit ``(n_out, n_in)`` matrix so that the same weights drive both the
    values and the mask. A mean rather than a sum because the containers carry
    a flux *density* (Jy): summing would make the result depend on the input
    sampling.

    The target grid is given at construction and is expected to be the observed
    container's own coordinates — ``transformations.md`` §10 requires a step
    reproducing observed coordinates to take them from the observed container
    rather than recompute them, so build one of these with
    ``observed.spectral_axis.values``.

    Parameters
    ----------
    target
        The output coordinates, micron. Strictly increasing.
    label
        Component label for this step.
    """

    ACCEPTS: ClassVar[tuple[type, ...]] = (Spectrum,)

    def __init__(
        self,
        target: Any,
        *,
        label: str | None = None,
        dtype: torch.dtype = DEFAULT_DTYPE,
        device: torch.device = DEFAULT_DEVICE,
    ) -> None:
        super().__init__(label=label, dtype=dtype, device=device)
        grid = _to_micron(target)
        if grid.ndim != 1 or grid.size == 0:
            raise TransformationError(
                f"Resample needs a one-dimensional, non-empty target grid, got shape {grid.shape}."
            )
        if grid.size > 1 and not bool(np.all(np.diff(grid) > 0.0)):
            raise TransformationError(
                "Resample's target grid must be strictly increasing; a Spectrum's spectral axis "
                "is strictly increasing, so an unsorted target could not be used to build one."
            )
        self._declare("target", grid, unit=COORDINATE_UNIT)

    def requirements(self) -> tuple[AxisRequirement, ...]:
        """Cover the target range, sampled finer than its finest output bin.

        A density requirement rather than ``points=``: the weights are built
        from whatever grid arrives, so this step has no positional buffer to
        keep aligned and does not need the exact coordinates back.
        """
        grid = self._data("target")
        step = float(np.diff(grid).min()) / 2.0 if grid.size > 1 else None
        return (
            AxisRequirement(
                "spectral_axis",
                unit=COORDINATE_UNIT,
                intervals=(float(grid[0]), float(grid[-1])),
                max_step=step,
            ),
        )

    def influence_tensor(self, source: torch.Tensor) -> torch.Tensor:
        """The ``(n_out, n_in)`` weight matrix taking *source* onto the target."""
        target = self._tensor("target")
        source_edges = bin_edges(source)
        target_edges = bin_edges(target)
        lower = torch.maximum(target_edges[:-1, None], source_edges[None, :-1])
        upper = torch.minimum(target_edges[1:, None], source_edges[None, 1:])
        overlap = torch.clamp(upper - lower, min=0.0)
        total = overlap.sum(dim=1, keepdim=True)
        # A target bin that overlaps nothing (a gap in the input, or a target
        # point outside its range) falls back to the single nearest input
        # sample, so no output row is empty and no output value is silently
        # zero -- which would read as a real flux of zero downstream.
        empty = (total <= 0.0).reshape(-1)
        if bool(torch.any(empty)):
            nearest = torch.argmin(torch.abs(target[empty, None] - source[None, :]), dim=1)
            overlap = overlap.clone()
            overlap[torch.nonzero(empty, as_tuple=True)[0], nearest] = 1.0
            total = overlap.sum(dim=1, keepdim=True)
        return overlap / total

    def influence(self, source: Any) -> np.ndarray:
        """:meth:`influence_tensor` on the numpy side of the boundary."""
        return to_numpy(
            self.influence_tensor(as_tensor(source, dtype=self.dtype, device=self.device))
        )

    def apply_tensor(self, weights: torch.Tensor, values: torch.Tensor) -> torch.Tensor:
        """The differentiable step: the binning integral as one ``matmul``."""
        return weights @ values

    def apply(self, samples: Any, values: Any) -> Spectrum:
        source = as_tensor(samples.spectral_axis.values, dtype=self.dtype, device=self.device)
        weights = self.influence_tensor(source)
        flux = as_tensor(samples.values, dtype=self.dtype, device=self.device)
        return Spectrum(
            self._data("target") * COORDINATE_UNIT,
            to_numpy(self.apply_tensor(weights, flux)),
            unit=samples.unit,
            mask=propagate_mask(samples, to_numpy(weights)),
        )


class LSFConvolution(TorchStep):
    """Convolve with a Gaussian line-spread function.

    Either a constant ``resolving_power`` (``R = lambda / FWHM``, the usual
    description of a grating spectrograph, so the kernel widens with
    wavelength) or a constant ``fwhm`` in micron. Exactly one of the two.

    The kernel is built as an explicit matrix against the input grid rather
    than as an FFT, because the input grid is whatever negotiation produced and
    need not be evenly spaced.

    The width is a **buffer**, not a parameter: an LSF known from the
    instrument's own calibration is constant data. Promoting it to a fitted
    parameter later is a configuration change (``parameters.md`` §10), not a
    rewrite — call ``promote_buffer``.

    Parameters
    ----------
    resolving_power
        ``lambda / FWHM``, dimensionless. Mutually exclusive with *fwhm*.
    fwhm
        Constant full width at half maximum, micron. Mutually exclusive with
        *resolving_power*.
    label
        Component label for this step.
    """

    ACCEPTS: ClassVar[tuple[type, ...]] = (Spectrum,)

    #: FWHM of a Gaussian in units of its standard deviation.
    FWHM_PER_SIGMA: ClassVar[float] = float(2.0 * np.sqrt(2.0 * np.log(2.0)))

    def __init__(
        self,
        *,
        resolving_power: float | None = None,
        fwhm: float | None = None,
        label: str | None = None,
        dtype: torch.dtype = DEFAULT_DTYPE,
        device: torch.device = DEFAULT_DEVICE,
    ) -> None:
        super().__init__(label=label, dtype=dtype, device=device)
        if (resolving_power is None) == (fwhm is None):
            raise TransformationError(
                "LSFConvolution takes exactly one of resolving_power (constant lambda/FWHM) or "
                "fwhm (constant width in micron); it was given "
                + ("both" if fwhm is not None else "neither")
                + "."
            )
        if resolving_power is not None:
            value = float(resolving_power)
            if not np.isfinite(value) or value <= 0.0:
                raise TransformationError(
                    f"resolving_power must be finite and positive, got {resolving_power!r}."
                )
            self._declare("resolving_power", value)
        else:
            value = float(_to_micron(fwhm))
            if not np.isfinite(value) or value <= 0.0:
                raise TransformationError(f"fwhm must be finite and positive, got {fwhm!r}.")
            self._declare("fwhm", value, unit=COORDINATE_UNIT)
        self._span: tuple[float, float] | None = None

    def power(self) -> float | None:
        """``lambda / FWHM``, or ``None`` for a constant-width kernel."""
        if "resolving_power" not in self.buffers:
            return None
        return float(self.buffers["resolving_power"].value)

    def width(self) -> float | None:
        """Constant FWHM in micron, or ``None`` for a constant-R kernel."""
        if "fwhm" not in self.buffers:
            return None
        return float(self.buffers["fwhm"].value)

    def sigma_tensor(self, wavelength: torch.Tensor) -> torch.Tensor:
        """Kernel standard deviation, micron, at each of *wavelength*."""
        power = self.power()
        if power is not None:
            widths = wavelength / self._tensor("resolving_power")
        else:
            widths = torch.full_like(wavelength, float(self.width()))
        return widths / self.FWHM_PER_SIGMA

    def sigma(self, wavelength: Any) -> np.ndarray:
        """:meth:`sigma_tensor` on the numpy side of the boundary."""
        return to_numpy(
            self.sigma_tensor(as_tensor(wavelength, dtype=self.dtype, device=self.device))
        )

    def configure_from(self, downstream: Sequence[Transformation]) -> None:
        """Learn the output range from the steps after this one.

        A convolution needs its input to extend *beyond* the range anyone
        actually wants, by several kernel widths, or the outermost output
        samples are convolved against an edge instead of against data. An LSF
        has no idea what that range is: it is set by the resampler or the
        photometry step downstream.

        ``Instrument.__init__`` calls this once, **last step first** (ruled
        2026-09-07), with this step's successors, before the hot loop. Nothing
        is inferred backwards through the chain — the successors' own published
        requirements are read, and if none of them constrains the spectral axis
        this step stays silent.
        """
        low: float | None = None
        high: float | None = None
        for step in downstream:
            for requirement in step.requirements():
                if requirement.axis != "spectral_axis":
                    continue
                if not requirement.intervals and requirement.points is None:
                    continue
                coordinates = np.asarray(
                    requirement.convert_to(COORDINATE_UNIT).coordinates().to_value(COORDINATE_UNIT),
                    dtype=DTYPE,
                )
                if coordinates.size == 0:
                    continue
                low = coordinates.min() if low is None else min(low, float(coordinates.min()))
                high = coordinates.max() if high is None else max(high, float(coordinates.max()))
        self._span = None if low is None or high is None else (float(low), float(high))

    def requirements(self) -> tuple[AxisRequirement, ...]:
        """Cover the downstream range padded by several kernel widths, sampled finely.

        Empty until :meth:`configure_from` has found a downstream range — a
        bare LSF on its own genuinely does not know what to ask for, and
        publishing a guess would over-constrain negotiation.
        """
        if self._span is None:
            return ()
        low, high = self._span
        pad_low = _LSF_PADDING_SIGMAS * float(self.sigma(np.array([low]))[0])
        pad_high = _LSF_PADDING_SIGMAS * float(self.sigma(np.array([high]))[0])
        power = self.power()
        if power is not None:
            return (
                AxisRequirement(
                    "spectral_axis",
                    unit=COORDINATE_UNIT,
                    intervals=(max(low - pad_low, 0.0), high + pad_high),
                    min_resolving_power=_LSF_OVERSAMPLING * power * self.FWHM_PER_SIGMA,
                ),
            )
        return (
            AxisRequirement(
                "spectral_axis",
                unit=COORDINATE_UNIT,
                intervals=(max(low - pad_low, 0.0), high + pad_high),
                max_step=float(self.width()) / _LSF_OVERSAMPLING,
            ),
        )

    def influence_tensor(self, source: torch.Tensor) -> torch.Tensor:
        """The normalised ``(n, n)`` convolution matrix on the grid *source*."""
        sigma = self.sigma_tensor(source)
        separation = source[:, None] - source[None, :]
        kernel = torch.exp(-0.5 * (separation / sigma[:, None]) ** 2)
        # Integrate against the input's own bin widths, so an unevenly sampled
        # grid does not weight its dense regions more heavily.
        widths = torch.diff(bin_edges(source))
        weighted = kernel * widths[None, :]
        return weighted / weighted.sum(dim=1, keepdim=True)

    def influence(self, source: Any) -> np.ndarray:
        """:meth:`influence_tensor` on the numpy side of the boundary."""
        return to_numpy(
            self.influence_tensor(as_tensor(source, dtype=self.dtype, device=self.device))
        )

    def apply_tensor(self, weights: torch.Tensor, values: torch.Tensor) -> torch.Tensor:
        """The differentiable step: the convolution as one ``matmul``."""
        return weights @ values

    def apply(self, samples: Any, values: Any) -> Spectrum:
        source = as_tensor(samples.spectral_axis.values, dtype=self.dtype, device=self.device)
        weights = self.influence_tensor(source)
        flux = as_tensor(samples.values, dtype=self.dtype, device=self.device)
        return samples.with_values(
            to_numpy(self.apply_tensor(weights, flux)),
            mask=propagate_mask(samples, to_numpy(weights)),
        )


class SyntheticPhotometry(TorchStep):
    """Integrate a spectrum through filter response curves.

    ``Spectrum -> PhotometricPoints``: the kind-changing step. Each output is a
    normalised response-weighted mean of the spectrum over one filter's
    support.

    Which mean depends on **what the detector counts**, and the two answers are
    different numbers, not different spellings — see :data:`DETECTORS` and the
    ``detector`` argument. There is no default (ruled by Peter 2026-09-05): a
    silently chosen convention is exactly the failure this argument exists to
    prevent.

    The response curves are **buffers** tabulated on a particular wavelength
    grid, so this class is ``spectrum_photometry.md`` Gap 1's shape exactly: it
    publishes ``points=`` naming that tabulation, and reads the compiled
    container through :meth:`~ampere.core.Axis.locate` rather than
    positionally. Without that, binding a second instrument to the same channel
    silently changes which flux each response column multiplies.

    Parameters
    ----------
    filters
        Filter names, one per output point. Must be unique.
    wavelength
        The tabulation grid the responses are given on, micron.
    response
        ``(n_filters, n_wavelength)`` transmission. Need not be normalised.
    detector
        **Required.** ``"photon"`` or ``"energy"``, applied to every filter; or
        one such string per filter, since a real filter set mixes types.
    pivots
        Pivot wavelengths, micron, one per filter. Computed from the responses
        when omitted; the pivot is convention-independent.
    label
        Component label for this step.
    """

    ACCEPTS: ClassVar[tuple[type, ...]] = (Spectrum,)
    PRODUCES: ClassVar[type] = PhotometricPoints

    def __init__(
        self,
        filters: Sequence[str],
        wavelength: Any,
        response: Any,
        *,
        detector: str | Sequence[str],
        pivots: Any = None,
        label: str | None = None,
        dtype: torch.dtype = DEFAULT_DTYPE,
        device: torch.device = DEFAULT_DEVICE,
    ) -> None:
        super().__init__(label=label, dtype=dtype, device=device)
        names = tuple(str(name) for name in filters)
        if len(set(names)) != len(names):
            raise TransformationError(
                f"SyntheticPhotometry was given repeated filter names {sorted(names)}; a filter "
                f"name is the identity of a photometric point, so they must be unique."
            )
        kinds = _detector_kinds(detector, names)
        grid = as_tensor(_to_micron(wavelength), dtype=dtype, device=device)
        curves = as_tensor(np.asarray(response, dtype=DTYPE), dtype=dtype, device=device)
        if curves.ndim != 2 or tuple(curves.shape) != (len(names), int(grid.numel())):
            raise TransformationError(
                f"SyntheticPhotometry needs a ({len(names)}, {int(grid.numel())}) response array "
                f"for {len(names)} filter(s) on a {int(grid.numel())}-point grid, got shape "
                f"{tuple(curves.shape)}."
            )
        if not bool(torch.all(torch.isfinite(curves))) or bool(torch.any(curves < 0.0)):
            raise TransformationError(
                "SyntheticPhotometry's response curves must be finite and non-negative."
            )
        if not bool(torch.all(grid > 0.0)):
            raise TransformationError(
                "SyntheticPhotometry needs strictly positive wavelengths: both detector "
                "conventions divide by the wavelength."
            )
        widths = torch.diff(bin_edges(grid))
        # The whole of the convention, in one exponent: 1 for a photon
        # counter's R dlambda / lambda, 2 for an energy detector's
        # R dlambda / lambda**2. See DETECTORS for the derivation.
        exponent = as_tensor(
            [DETECTORS[kind] for kind in kinds], dtype=dtype, device=device
        ).reshape(-1, 1)
        raw = curves * widths[None, :] / grid[None, :] ** exponent
        norms = raw.sum(dim=1)
        if not bool(torch.all(norms > 0.0)):
            empty = [names[i] for i in to_numpy(torch.nonzero(norms <= 0.0).reshape(-1))]
            raise TransformationError(
                f"filter(s) {empty} have a response that integrates to zero over the grid they "
                f"were tabulated on, so no flux could ever be measured through them."
            )
        self._names = names
        self._detectors = kinds
        self._declare("wavelength", to_numpy(grid), unit=COORDINATE_UNIT)
        self._declare("response", to_numpy(curves))
        # The weights are registered as constant data in their own right, not
        # merely derived: two steps with identical responses and *different*
        # conventions compute different numbers, and `response` alone would
        # make them indistinguishable to provenance (results.md §13.13).
        self._declare("weights", to_numpy(raw / norms[:, None]))
        self._declare(
            "pivot",
            _to_micron(pivots) if pivots is not None else to_numpy(_pivot(grid, curves, widths)),
            unit=COORDINATE_UNIT,
        )

    def describe(self) -> dict[str, list[str]]:
        """The detector convention per filter — configuration, not data.

        ``results.md`` §13.13's hook, for exactly the case it was landed for: a
        plain Python attribute that changes what this step computes and is
        neither a parameter nor an array.
        """
        return {"detector": list(self._detectors)}

    @property
    def detectors(self) -> tuple[str, ...]:
        """The detector convention used for each filter, in output order."""
        return self._detectors

    @property
    def filters(self) -> tuple[str, ...]:
        """The filter names, in output order."""
        return self._names

    def tabulation(self) -> np.ndarray:
        """The grid the responses are tabulated on, micron."""
        return self._data("wavelength")

    def pivots(self) -> np.ndarray:
        """Pivot wavelengths, micron — the output container's coordinates."""
        return self._data("pivot")

    def requirements(self) -> tuple[AxisRequirement, ...]:
        """``points=`` at the exact tabulation the response buffer was built on."""
        return (AxisRequirement("spectral_axis", unit=COORDINATE_UNIT, points=self.tabulation()),)

    def influence_tensor(self, axis: Any) -> torch.Tensor:
        """The ``(n_filters, n_samples)`` matrix against the whole of *axis*.

        *axis* is an :class:`~ampere.core.Axis`, not a bare array, because the
        lookup is its job: the columns are placed at the indices ``locate``
        reports and every other column is zero, so the same matrix is correct
        for the values and for the mask however much larger than this step's
        own tabulation the negotiated grid turned out to be.
        """
        index = torch.as_tensor(
            np.asarray(axis.locate(self.tabulation()), dtype=np.int64), device=self.device
        )
        full = torch.zeros((len(self._names), int(axis.size)), dtype=self.dtype, device=self.device)
        full[:, index] = self._tensor("weights")
        return full

    def influence(self, axis: Any) -> np.ndarray:
        """:meth:`influence_tensor` on the numpy side of the boundary."""
        return to_numpy(self.influence_tensor(axis))

    def apply_tensor(self, weights: torch.Tensor, values: torch.Tensor) -> torch.Tensor:
        """The differentiable step: the filter integrals as one ``matmul``."""
        return weights @ values

    def apply(self, samples: Any, values: Any) -> PhotometricPoints:
        weights = self.influence_tensor(samples.spectral_axis)
        flux = as_tensor(samples.values, dtype=self.dtype, device=self.device)
        return PhotometricPoints(
            self._names,
            self.pivots() * COORDINATE_UNIT,
            to_numpy(self.apply_tensor(weights, flux)),
            unit=samples.unit,
            mask=propagate_mask(samples, to_numpy(weights)),
        )

    @classmethod
    def from_library(
        cls,
        filters: Sequence[str],
        wavelength: Any,
        *,
        detector: str | Sequence[str] | None = None,
        library: Any = None,
        label: str | None = None,
        dtype: torch.dtype = DEFAULT_DTYPE,
        device: torch.device = DEFAULT_DEVICE,
    ) -> SyntheticPhotometry:
        """Build a step from a pyphot filter library, interpolated onto *wavelength*.

        ``detector`` is optional here, unlike in the constructor: each pyphot
        ``Filter`` records its own convention in ``dtype``, and the bundled
        library genuinely mixes them (2MASS is ``photon``; AKARI and IRAS are
        ``energy``). A filter whose metadata is missing or unrecognised is
        refused rather than guessed at.
        """
        # Imported here rather than at module level: pyphot is a base
        # dependency, but it opens HDF5 filter libraries and is much heavier
        # than the rest of this module, which a chain with no photometry in it
        # should not pay for.
        import pyphot

        from ampere.utils.pyphot_compat import get_unit

        grid = _to_micron(wavelength)
        if library is None:
            library = pyphot.get_library(fname=str(bundled_filter_library()))
        elif isinstance(library, (str, bytes)) or hasattr(library, "__fspath__"):
            library = pyphot.get_library(fname=str(library))
        names = list(filters)
        curves = library.load_filters(names, interp=True, lamb=grid * get_unit("micron"))
        response = np.vstack([np.asarray(curve.transmit, dtype=DTYPE) for curve in curves])
        pivots = np.array(
            [float(curve.lpivot.to(get_unit("micron")).value) for curve in curves], dtype=DTYPE
        )
        if detector is None:
            detector = [
                _library_detector(name, curve) for name, curve in zip(names, curves, strict=True)
            ]
        return cls(
            names,
            grid,
            response,
            detector=detector,
            pivots=pivots,
            label=label,
            dtype=dtype,
            device=device,
        )


def _library_detector(name: str, curve: Any) -> str:
    """The convention a pyphot ``Filter`` declares, or a loud refusal."""
    declared = getattr(curve, "dtype", None)
    if declared not in DETECTORS:
        raise TransformationError(
            f"filter {name!r} does not declare a usable detector type (its pyphot dtype is "
            f"{declared!r}, expected one of {sorted(DETECTORS)}). The photon and energy "
            f"conventions give different fluxes, so this is refused rather than guessed; "
            f"pass detector= explicitly to say which one this filter uses."
        )
    return str(declared)


def _detector_kinds(detector: str | Sequence[str], names: Sequence[str]) -> tuple[str, ...]:
    """Normalise *detector* into one convention per filter, refusing anything else."""
    if isinstance(detector, str):
        kinds = (detector,) * len(names)
    else:
        kinds = tuple(str(kind) for kind in detector)
        if len(kinds) != len(names):
            raise TransformationError(
                f"SyntheticPhotometry was given {len(kinds)} detector type(s) for "
                f"{len(names)} filter(s). Pass one string to use the same convention for "
                f"all of them, or exactly one per filter."
            )
    unknown = {kind for kind in kinds if kind not in DETECTORS}
    if unknown:
        raise TransformationError(
            f"unknown detector type(s) {sorted(unknown)}. A detector either counts photons "
            f"('photon': weight R dlambda/lambda) or measures energy ('energy': weight "
            f"R dlambda/lambda**2), and the two give different numbers, so there is no "
            f"default and no third option."
        )
    return kinds


def bundled_filter_library() -> Any:
    """Path to ``ampere_allfilters.hd5``, the filter set ampere ships.

    Located through :mod:`importlib.resources` rather than by string surgery on
    ``ampere.__file__``, which breaks on any installation layout that is not a
    plain directory.
    """
    from importlib import resources

    return resources.files("ampere") / "ampere_allfilters.hd5"


def _pivot(grid: torch.Tensor, response: torch.Tensor, widths: torch.Tensor) -> torch.Tensor:
    """Response-weighted pivot wavelength of each filter, micron.

    ``lambda_pivot**2 = int(R lambda dlambda) / int(R dlambda / lambda)`` — the
    wavelength at which a flux density in ``f_nu`` and one in ``f_lambda``
    agree, and the conventional coordinate to report a broadband point at.
    """
    weighted = response * widths[None, :]
    numerator = (weighted * grid[None, :]).sum(dim=1)
    denominator = (weighted / grid[None, :]).sum(dim=1)
    return torch.sqrt(numerator / denominator)


def _to_micron(coordinates: Any) -> np.ndarray:
    """Coordinates as bare micron, whether or not they arrive as a Quantity."""
    if isinstance(coordinates, u.Quantity):
        return np.asarray(coordinates.to_value(COORDINATE_UNIT), dtype=DTYPE)
    if isinstance(coordinates, torch.Tensor):
        return to_numpy(coordinates).astype(DTYPE, copy=False)
    return np.asarray(coordinates, dtype=DTYPE)
