"""The standard instrument steps on the numpy path.

``DEVELOPMENT_PLAN.md`` §5's Phase 2 list — "resampling, LSF, synthetic
photometry, calibration factor" — with the shapes ``transformations.md`` §10's
table fixes. ``ampere.core`` deliberately ships no concrete transformation:
core is the vocabulary, and this is the library.

Two rules from the contracts govern everything here.

**Coordinates come from the container, never from arithmetic.**
``transformations.md`` §10 (W1.11 gap I-2): a step reproducing the observed
coordinates takes them *from* the observed container, because
``check_alignment`` compares axes with ``np.array_equal`` and a recomputed grid
differs in its last bits — failing a check whose message is about axes rather
than about arithmetic.

**A step whose buffer is tabulated on particular coordinates publishes
``points=`` and looks them up.** ``spectrum_photometry.md`` Gap 1: once a
second instrument binds the same channel, negotiation's union hands the step a
larger, possibly reordered grid, and reading ``samples.values`` positionally
against a response matrix is then silently wrong. :class:`SyntheticPhotometry`
is this backend's instance of that shape, and
:meth:`~ampere.core.Axis.locate` — landed with this item — is how it stays
correct.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, ClassVar

import astropy.units as u
import numpy as np

from ampere.core import (
    DTYPE,
    AxisRequirement,
    PhotometricPoints,
    Spectrum,
    Transformation,
    TransformationError,
    propagate_mask,
)

from ._declare import as_parameter

__all__ = [
    "CalibrationScale",
    "LSFConvolution",
    "Resample",
    "SyntheticPhotometry",
    "bin_edges",
    "bundled_filter_library",
]

COORDINATE_UNIT = u.micron

#: How far beyond the output range an LSF needs coverage, in kernel sigmas.
#: Truncating a Gaussian at 5 sigma loses about 6e-7 of its mass, which is well
#: below anything a spectral fit can see.
_LSF_PADDING_SIGMAS = 5.0
#: How much finer than the kernel the input grid must be, so the convolution is
#: not itself limited by sampling.
_LSF_OVERSAMPLING = 3.0


def bin_edges(centres: np.ndarray) -> np.ndarray:
    """Bin edges bracketing *centres*: midpoints, with the outer two reflected.

    The usual convention for a spectrum given as sample centres. A single
    sample has no midpoints to work with, so it is given a unit-width bin —
    the width then cancels out of every normalised weight that uses it.
    """
    grid = np.asarray(centres, dtype=DTYPE)
    if grid.size == 1:
        return np.array([grid[0] - 0.5, grid[0] + 0.5], dtype=DTYPE)
    middle = 0.5 * (grid[1:] + grid[:-1])
    first = grid[0] - (middle[0] - grid[0])
    last = grid[-1] + (grid[-1] - middle[-1])
    return np.concatenate([[first], middle, [last]])


class _Step(Transformation):
    """Shared plumbing: read a declared buffer as a plain array.

    Buffers are the declared surface for constant data (``parameters.md`` §10),
    and a buffer's name may not shadow a class attribute — so a step whose
    buffer is called ``target`` cannot also expose a ``target`` property. This
    is the readable spelling of ``self.buffers[name].value`` that the steps use
    instead.
    """

    def _data(self, name: str) -> np.ndarray:
        return np.asarray(self.buffers[name].value, dtype=DTYPE)


class CalibrationScale(_Step):
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

    def __init__(self, scale: Any = 1.0, *, label: str | None = None) -> None:
        super().__init__(label=label)
        self.register_parameter(as_parameter("scale", scale))

    def apply(self, samples: Any, values: Any) -> Any:
        factor = float(self.context(values)["scale"])
        # with_values inherits the axes, the unit and the mask, so the
        # calibration costs one multiply and no revalidation.
        return samples.with_values(samples.values * factor)


class Resample(_Step):
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

    def __init__(self, target: Any, *, label: str | None = None) -> None:
        super().__init__(label=label)
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
        self.register_buffer("target", grid, unit=COORDINATE_UNIT)

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

    def influence(self, source: Any) -> np.ndarray:
        """The ``(n_out, n_in)`` weight matrix taking *source* onto the target."""
        incoming = np.asarray(source, dtype=DTYPE)
        source_edges = bin_edges(incoming)
        target_edges = bin_edges(self._data("target"))
        lower = np.maximum(target_edges[:-1, None], source_edges[None, :-1])
        upper = np.minimum(target_edges[1:, None], source_edges[None, 1:])
        overlap = np.clip(upper - lower, 0.0, None)
        total = overlap.sum(axis=1, keepdims=True)
        # A target bin that overlaps nothing (a gap in the input, or a target
        # point outside its range) falls back to the single nearest input
        # sample, so no output row is empty and no output value is silently
        # zero -- which would read as a real flux of zero downstream.
        empty = (total <= 0.0).ravel()
        if bool(np.any(empty)):
            target = self._data("target")
            nearest = np.argmin(np.abs(target[empty, None] - incoming[None, :]), axis=1)
            overlap[np.flatnonzero(empty), nearest] = 1.0
            total = overlap.sum(axis=1, keepdims=True)
        return overlap / total

    def apply(self, samples: Any, values: Any) -> Spectrum:
        weights = self.influence(samples.spectral_axis.values)
        return Spectrum(
            self._data("target") * COORDINATE_UNIT,
            weights @ np.asarray(samples.values, dtype=DTYPE),
            unit=samples.unit,
            mask=propagate_mask(samples, weights),
        )


class LSFConvolution(_Step):
    """Convolve with a Gaussian line-spread function.

    Either a constant ``resolving_power`` (``R = lambda / FWHM``, the usual
    description of a grating spectrograph, so the kernel widens with
    wavelength) or a constant ``fwhm`` in micron. Exactly one of the two.

    The kernel is built as an explicit matrix against the input grid rather
    than as an FFT, because the input grid is whatever negotiation produced and
    need not be evenly spaced. ``Axis.regular``/``log_regular`` advertise when
    a faster path would be valid; the reference backend deliberately does not
    take it (``architecture.md`` §2: correctness is its only goal).

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
    FWHM_PER_SIGMA: ClassVar[float] = 2.0 * np.sqrt(2.0 * np.log(2.0))

    def __init__(
        self,
        *,
        resolving_power: float | None = None,
        fwhm: float | None = None,
        label: str | None = None,
    ) -> None:
        super().__init__(label=label)
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
            self.register_buffer("resolving_power", value)
        else:
            value = float(_to_micron(fwhm))
            if not np.isfinite(value) or value <= 0.0:
                raise TransformationError(f"fwhm must be finite and positive, got {fwhm!r}.")
            self.register_buffer("fwhm", value, unit=COORDINATE_UNIT)
        self._span: tuple[float, float] | None = None

    def power(self) -> float | None:
        """``lambda / FWHM``, or ``None`` for a constant-width kernel.

        A method rather than a property because ``resolving_power`` is a
        *buffer* name, and a parameter or buffer may not shadow a class
        attribute (``parameters.md`` §10). The declared surface is
        ``self.buffers``; these are the readable spellings of it.
        """
        if "resolving_power" not in self.buffers:
            return None
        return float(self.buffers["resolving_power"].value)

    def width(self) -> float | None:
        """Constant FWHM in micron, or ``None`` for a constant-R kernel."""
        if "fwhm" not in self.buffers:
            return None
        return float(self.buffers["fwhm"].value)

    def sigma(self, wavelength: Any) -> np.ndarray:
        """Kernel standard deviation, micron, at each of *wavelength*."""
        grid = np.asarray(wavelength, dtype=DTYPE)
        power = self.power()
        width = grid / power if power is not None else np.full(grid.shape, float(self.width()))
        return width / self.FWHM_PER_SIGMA

    def configure_from(self, downstream: Sequence[Transformation]) -> None:
        """Learn the output range from the steps after this one.

        A convolution needs its input to extend *beyond* the range anyone
        actually wants, by several kernel widths, or the outermost output
        samples are convolved against an edge instead of against data. But an
        LSF has no idea what that range is: it is set by the resampler or the
        photometry step downstream.

        This is what ``configure_from`` is for (``transformations.md`` §5, gap
        I-3's mechanism): ``Instrument.__init__`` calls it once with this
        step's successors, before the hot loop, and the padded requirement is
        published from what they ask for. Nothing is inferred backwards through
        the chain — the successors' own published requirements are read, and if
        none of them constrains the spectral axis this step stays silent.
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

    def influence(self, source: Any) -> np.ndarray:
        """The normalised ``(n, n)`` convolution matrix on the grid *source*."""
        grid = np.asarray(source, dtype=DTYPE)
        sigma = self.sigma(grid)
        separation = grid[:, None] - grid[None, :]
        kernel = np.exp(-0.5 * (separation / sigma[:, None]) ** 2)
        # Integrate against the input's own bin widths, so an unevenly sampled
        # grid does not weight its dense regions more heavily.
        widths = np.diff(bin_edges(grid))
        weighted = kernel * widths[None, :]
        total = weighted.sum(axis=1, keepdims=True)
        return weighted / total

    def apply(self, samples: Any, values: Any) -> Spectrum:
        weights = self.influence(samples.spectral_axis.values)
        return samples.with_values(
            weights @ np.asarray(samples.values, dtype=DTYPE),
            mask=propagate_mask(samples, weights),
        )


class SyntheticPhotometry(_Step):
    """Integrate a spectrum through filter response curves.

    ``Spectrum -> PhotometricPoints``: the kind-changing step. Each output is
    the response-weighted mean of the spectrum over one filter's support, which
    for a flux density in Jy is the photometric convention.

    The response curves are **buffers** tabulated on a particular wavelength
    grid, so this class is ``spectrum_photometry.md`` Gap 1's shape exactly: it
    publishes ``points=`` naming that tabulation, and reads the compiled
    container through :meth:`~ampere.core.Axis.locate` rather than
    positionally. Without that, binding a second instrument to the same channel
    silently changes which flux each response column multiplies.

    Build one from ampere's bundled filter library with :meth:`from_library`,
    which is the pyphot route; the constructor itself takes plain arrays, so a
    response measured in the lab needs no filter library at all.

    Parameters
    ----------
    filters
        Filter names, one per output point. Must be unique.
    wavelength
        The tabulation grid the responses are given on, micron.
    response
        ``(n_filters, n_wavelength)`` transmission. Need not be normalised.
    pivots
        Pivot wavelengths, micron, one per filter. Computed from the responses
        when omitted.
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
        pivots: Any = None,
        label: str | None = None,
    ) -> None:
        super().__init__(label=label)
        names = tuple(str(name) for name in filters)
        if len(set(names)) != len(names):
            raise TransformationError(
                f"SyntheticPhotometry was given repeated filter names {sorted(names)}; a filter "
                f"name is the identity of a photometric point, so they must be unique."
            )
        grid = _to_micron(wavelength)
        curves = np.asarray(response, dtype=DTYPE)
        if curves.ndim != 2 or curves.shape != (len(names), grid.size):
            raise TransformationError(
                f"SyntheticPhotometry needs a ({len(names)}, {grid.size}) response array for "
                f"{len(names)} filter(s) on a {grid.size}-point grid, got shape {curves.shape}."
            )
        if not bool(np.all(np.isfinite(curves))) or bool(np.any(curves < 0.0)):
            raise TransformationError(
                "SyntheticPhotometry's response curves must be finite and non-negative."
            )
        widths = np.diff(bin_edges(grid))
        norms = (curves * widths[None, :]).sum(axis=1)
        if not bool(np.all(norms > 0.0)):
            empty = [names[i] for i in np.flatnonzero(norms <= 0.0)]
            raise TransformationError(
                f"filter(s) {empty} have a response that integrates to zero over the grid they "
                f"were tabulated on, so no flux could ever be measured through them."
            )
        self._names = names
        self.register_buffer("wavelength", grid, unit=COORDINATE_UNIT)
        self.register_buffer("response", curves)
        self.register_buffer(
            "pivot",
            _to_micron(pivots) if pivots is not None else _pivot(grid, curves, widths),
            unit=COORDINATE_UNIT,
        )
        self._weights = curves * widths[None, :] / norms[:, None]

    @property
    def filters(self) -> tuple[str, ...]:
        """The filter names, in output order."""
        return self._names

    def tabulation(self) -> np.ndarray:
        """The grid the responses are tabulated on, micron.

        A method rather than a property: ``wavelength`` is a *buffer* name, and
        a buffer may not shadow a class attribute (``parameters.md`` §10).
        """
        return self._data("wavelength")

    def pivots(self) -> np.ndarray:
        """Pivot wavelengths, micron — the output container's coordinates."""
        return self._data("pivot")

    def requirements(self) -> tuple[AxisRequirement, ...]:
        """``points=`` at the exact tabulation the response buffer was built on.

        A density requirement would only guarantee *coverage*, and the union is
        free to satisfy that with any grid dense enough — which need not have
        this buffer's length or spacing (``spectrum_photometry.md`` §3).
        """
        return (
            AxisRequirement(
                "spectral_axis",
                unit=COORDINATE_UNIT,
                points=self.tabulation(),
            ),
        )

    def influence(self, axis: Any) -> np.ndarray:
        """The ``(n_filters, n_samples)`` matrix against the whole of *axis*.

        *axis* is an :class:`~ampere.core.Axis`, not a bare array, because the
        lookup is its job: the columns are placed at the indices ``locate``
        reports, and every other column is zero. That makes the same matrix
        correct for the values and for the mask, however much larger than this
        step's own tabulation the negotiated grid turned out to be.
        """
        index = axis.locate(self.tabulation())
        full = np.zeros((len(self._names), int(axis.size)), dtype=DTYPE)
        full[:, index] = self._weights
        return full

    def apply(self, samples: Any, values: Any) -> PhotometricPoints:
        weights = self.influence(samples.spectral_axis)
        return PhotometricPoints(
            self._names,
            self.pivots() * COORDINATE_UNIT,
            weights @ np.asarray(samples.values, dtype=DTYPE),
            unit=samples.unit,
            mask=propagate_mask(samples, weights),
        )

    @classmethod
    def from_library(
        cls,
        filters: Sequence[str],
        wavelength: Any,
        *,
        library: Any = None,
        label: str | None = None,
    ) -> SyntheticPhotometry:
        """Build a step from a pyphot filter library, interpolated onto *wavelength*.

        Targets pyphot >= 2's unit-adapter API through
        :mod:`ampere.utils.pyphot_compat`, whose ``get_unit`` is the supported
        way to attach units to an array for pyphot — never the ``pyphot.unit``
        registry, which pyphot 2 removed.

        Parameters
        ----------
        filters
            Names as the library knows them.
        wavelength
            The grid to tabulate the responses on, micron. This becomes the
            step's ``points=`` requirement, so choose it to resolve the
            narrowest filter in the set.
        library
            An open pyphot library, or a path to one. Defaults to ampere's
            bundled ``ampere_allfilters.hd5``.
        label
            Component label for this step.
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
        curves = library.load_filters(list(filters), interp=True, lamb=grid * get_unit("micron"))
        response = np.vstack([np.asarray(curve.transmit, dtype=DTYPE) for curve in curves])
        pivots = np.array(
            [float(curve.lpivot.to(get_unit("micron")).value) for curve in curves], dtype=DTYPE
        )
        return cls(list(filters), grid, response, pivots=pivots, label=label)


def bundled_filter_library() -> Any:
    """Path to ``ampere_allfilters.hd5``, the filter set ampere ships.

    Located through :mod:`importlib.resources` rather than by string surgery on
    ``ampere.__file__``, which is the idiom the legacy code uses and which
    breaks on any installation layout that is not a plain directory.
    """
    from importlib import resources

    return resources.files("ampere") / "ampere_allfilters.hd5"


def _pivot(grid: np.ndarray, response: np.ndarray, widths: np.ndarray) -> np.ndarray:
    """Response-weighted pivot wavelength of each filter, micron.

    ``lambda_pivot**2 = int(R lambda dlambda) / int(R dlambda / lambda)`` — the
    wavelength at which a flux density in ``f_nu`` and one in ``f_lambda``
    agree, and the conventional coordinate to report a broadband point at.
    """
    weighted = response * widths[None, :]
    numerator = (weighted * grid[None, :]).sum(axis=1)
    denominator = (weighted / grid[None, :]).sum(axis=1)
    return np.sqrt(numerator / denominator)


def _to_micron(coordinates: Any) -> np.ndarray:
    """Coordinates as bare micron, whether or not they arrive as a Quantity."""
    if isinstance(coordinates, u.Quantity):
        return np.asarray(coordinates.to_value(COORDINATE_UNIT), dtype=DTYPE)
    return np.asarray(coordinates, dtype=DTYPE)
