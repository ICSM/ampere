"""The standard instrument steps on the jax path.

``transformations.md`` §10's table — calibration factor, resampling, LSF,
synthetic photometry — with the arithmetic in ``jax.numpy`` and the geometry
unchanged.

Why these subclass the reference steps
--------------------------------------
Each class here derives from its ``ampere.backends.reference`` counterpart and
overrides exactly two things: the four capability flags, and the arithmetic.
Everything else — the constructor validation, the published
:class:`~ampere.core.AxisRequirement`\\ s, ``configure_from``'s downstream
range-finding, ``SyntheticPhotometry``'s response tabulation, pivot
wavelengths and ``from_library`` loading — is inherited unchanged, because
none of it is arithmetic a gradient passes through. It is *geometry*: numpy
data preparation done once at composition time, off the hot loop, whose output
is a constant matrix.

That is a deliberate structural choice and worth defending rather than
assuming. Two backends whose published requirements were written twice would
eventually negotiate differently for the same declaration, and the conformance
suite would then be comparing two different problems and calling the
disagreement a numerical one. Sharing the declaration and forking the
arithmetic is what makes the cross-backend rows mean something. (``mirror.py``
in the conformance suite makes the same judgement for the same reason: "the
arithmetic is still the reference backend's, honestly, because there is nothing
different to offer there".) The import is safe in both directions
``architecture.md`` §4 cares about: the reference backend needs only the base
install, so importing it here adds no dependency, and nothing in the reference
backend imports this.

Two evaluation surfaces
-----------------------
As for the models (see :mod:`ampere.backends.jax.models`):

``apply``
    the contract surface. Returns an ``ampere.core`` container, whose values
    are coerced to numpy by the container itself — so a gradient stops here,
    and that is a property of the frozen container contract, not of this step.
``apply_flux``
    the native surface: a jax array of fluxes and its coordinates in, the
    transformed pair out. Pure, traceable, and what
    :mod:`ampere.backends.jax.problem` composes into the differentiable
    log-density.

Every influence matrix in this module is a function of the *coordinates* alone
— an LSF's width is a buffer, not a parameter (``parameters.md`` §10), and so
are a resampler's target grid and a filter's response curve — so each is built
once, converted to jax once, and cached against the grid it was built for.
``CalibrationScale`` is the one step whose action carries a free parameter, and
its action is a scalar multiply.
"""

from __future__ import annotations

from typing import Any, ClassVar

import astropy.units as u
import jax
import jax.numpy as jnp
import numpy as np

from ampere.backends.reference.instrument import (
    CalibrationScale as _ReferenceCalibrationScale,
)
from ampere.backends.reference.instrument import (
    LSFConvolution as _ReferenceLSFConvolution,
)
from ampere.backends.reference.instrument import (
    Resample as _ReferenceResample,
)
from ampere.backends.reference.instrument import (
    SyntheticPhotometry as _ReferenceSyntheticPhotometry,
)
from ampere.core import Axis, PhotometricPoints, Spectrum, propagate_mask

from ._config import BACKEND, require_x64

__all__ = [
    "COORDINATE_UNIT",
    "CalibrationScale",
    "LSFConvolution",
    "Resample",
    "SyntheticPhotometry",
]

COORDINATE_UNIT = u.micron


class _JaxStep:
    """The four capability flags and the influence-matrix cache, once.

    A mixin rather than a base class, so each concrete step still derives from
    its reference counterpart and inherits that class's declaration wholesale.
    """

    #: The gradient the whole backend exists for passes through
    #: :meth:`apply_flux`; see this module's docstring for what the flag means
    #: and where it stops.
    DIFFERENTIABLE: ClassVar[bool] = True
    #: **True since slice 2** (W2.5). The flag means one thing on this
    #: backend: ``jax.vmap`` over the realised density
    #: (:meth:`~ampere.backends.jax.problem.LoweredProblem.log_prob_unconstrained_batched`)
    #: evaluates a stack of parameter vectors in one call, and it is measured
    #: rather than asserted -- ``tests/backends/test_jax.py`` compares a vmapped
    #: density against the same density in a loop. It is true here because every
    #: operation in this class is whole-array ``jax.numpy``: nothing branches on
    #: a value, nothing indexes by one, so vmap maps it as it maps any pure
    #: function. ``QuasisepGP`` is the one part of this backend that still says
    #: False, and says why.
    BATCHABLE: ClassVar[bool] = True
    #: Never auto-detected (``architecture.md`` §5).
    DEVICE: ClassVar[str] = "cpu"
    BACKEND: ClassVar[str] = BACKEND

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        require_x64(f"a jax {type(self).__name__}")
        super().__init__(*args, **kwargs)
        self._cached_influence: tuple[bytes, jax.Array] | None = None

    def _influence_jax(self, source: Any) -> jax.Array:
        """This step's influence matrix on *source*, as a jax array, cached.

        The cache key is the coordinates' exact bytes. That is not an
        optimisation for its own sake: the matrix is a pure function of the
        grid, the grid is fixed for a run, and rebuilding an ``(n_out, n_in)``
        matrix per likelihood evaluation is precisely the hot-loop waste W2.3
        measured in the reference backend and told this track not to reproduce.
        An exact key means a changed grid rebuilds and a repeated one does not;
        there is no tolerance, so a lookup can never return the matrix of a
        nearby grid.
        """
        coordinates = np.asarray(_coordinates_of(source), dtype=float)
        key = coordinates.tobytes()
        cached = self._cached_influence
        if cached is not None and cached[0] == key:
            return cached[1]
        matrix = jnp.asarray(self.influence(source), dtype=jnp.float64)  # type: ignore[attr-defined]
        self._cached_influence = (key, matrix)
        return matrix


def _coordinates_of(source: Any) -> Any:
    """Bare coordinates, whether *source* is an ``Axis`` or already an array."""
    return getattr(source, "values", source)


class CalibrationScale(_JaxStep, _ReferenceCalibrationScale):
    """A multiplicative calibration factor: ``flux -> scale * flux``, in jax.

    The commonest instrumental nuisance there is, and the one step here whose
    action carries a free parameter. Kind-preserving, publishes no
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

    def apply_flux(
        self, flux: jax.Array, grid: jax.Array, values: Any
    ) -> tuple[jax.Array, jax.Array]:
        """The native surface: scale the fluxes, leave the coordinates alone."""
        factor = jnp.asarray(self.context(values)["scale"], dtype=jnp.float64)
        return flux * factor, grid

    def apply(self, samples: Any, values: Any) -> Any:
        factor = jnp.asarray(self.context(values)["scale"], dtype=jnp.float64)
        # with_values inherits the axes, the unit and the mask, so the
        # calibration costs one multiply and no revalidation.
        return samples.with_values(
            np.asarray(jnp.asarray(samples.values, dtype=jnp.float64) * factor)
        )


class Resample(_JaxStep, _ReferenceResample):
    """Resample a spectrum onto a coarser grid, conserving mean flux density.

    Each output sample is the overlap-weighted mean of the input samples whose
    bins intersect its own — the standard binning integral, as an explicit
    ``(n_out, n_in)`` matrix so that the same weights drive both the values and
    the mask. The matrix is built by the inherited (numpy) ``influence`` and
    then converted once; a resampling matrix is a function of two grids and of
    nothing a sampler proposes.

    Parameters
    ----------
    target
        The output coordinates, micron. Strictly increasing. Expected to be the
        observed container's own coordinates: ``transformations.md`` §10
        requires a step reproducing observed coordinates to take them from the
        observed container rather than recompute them, because a recomputed
        grid differs in its last bits and ``check_alignment`` compares axes
        with ``np.array_equal``.
    label
        Component label for this step.
    """

    def target_grid(self) -> jax.Array:
        """The output coordinates, micron, as a jax array."""
        return jnp.asarray(self.buffers["target"].value, dtype=jnp.float64)

    def apply_flux(
        self, flux: jax.Array, grid: jax.Array, values: Any
    ) -> tuple[jax.Array, jax.Array]:
        """The native surface: ``W @ flux`` onto the declared target grid."""
        return self._influence_jax(np.asarray(grid)) @ flux, self.target_grid()

    def apply(self, samples: Any, values: Any) -> Spectrum:
        weights = self._influence_jax(samples.spectral_axis.values)
        resampled = weights @ jnp.asarray(samples.values, dtype=jnp.float64)
        return Spectrum(
            np.asarray(self.buffers["target"].value) * COORDINATE_UNIT,
            np.asarray(resampled),
            unit=samples.unit,
            mask=propagate_mask(samples, np.asarray(weights)),
        )


class LSFConvolution(_JaxStep, _ReferenceLSFConvolution):
    """Convolve with a Gaussian line-spread function, in jax.

    Either a constant ``resolving_power`` (``R = lambda / FWHM``, so the kernel
    widens with wavelength) or a constant ``fwhm`` in micron; exactly one of
    the two, and the width is a **buffer** rather than a parameter because an
    LSF known from the instrument's own calibration is constant data. The
    kernel is an explicit matrix rather than an FFT, because the input grid is
    whatever negotiation produced and need not be evenly spaced.

    Promoting the width to a fitted parameter later is a configuration change
    (``parameters.md`` §10) — ``promote_buffer`` — but it would also make the
    matrix parameter-dependent, so the cache below would have to go with it.

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

    def apply_flux(
        self, flux: jax.Array, grid: jax.Array, values: Any
    ) -> tuple[jax.Array, jax.Array]:
        """The native surface: convolve on the incoming grid, which is unchanged."""
        return self._influence_jax(np.asarray(grid)) @ flux, grid

    def apply(self, samples: Any, values: Any) -> Spectrum:
        weights = self._influence_jax(samples.spectral_axis.values)
        convolved = weights @ jnp.asarray(samples.values, dtype=jnp.float64)
        return samples.with_values(
            np.asarray(convolved),
            mask=propagate_mask(samples, np.asarray(weights)),
        )


class SyntheticPhotometry(_JaxStep, _ReferenceSyntheticPhotometry):
    """Integrate a spectrum through filter response curves, in jax.

    ``Spectrum -> PhotometricPoints``: the kind-changing step. Each output is a
    normalised response-weighted mean of the spectrum over one filter's
    support, and *which* mean depends on what the detector counts — see
    :data:`ampere.backends.reference.DETECTORS`. There is no default detector
    convention, here or there: a silently chosen one is exactly the failure
    that argument exists to prevent.

    The response tabulation, the ``points=`` requirement published at it, the
    pivot wavelengths and :meth:`from_library` are all inherited. So is the
    column placement: ``influence`` builds its matrix against the *whole* of
    the negotiated axis using :meth:`~ampere.core.Axis.locate`, which is what
    keeps a response matrix correct after negotiation's union hands the step a
    larger, possibly reordered grid than the one it tabulated on
    (``spectrum_photometry.md`` Gap 1). Reading ``samples.values``
    positionally instead would be silently wrong, and this backend inherits the
    fix rather than re-deriving it.
    """

    def pivot_grid(self) -> jax.Array:
        """Pivot wavelengths, micron, as a jax array — the output coordinates."""
        return jnp.asarray(self.pivots(), dtype=jnp.float64)

    def apply_flux(
        self, flux: jax.Array, grid: jax.Array, values: Any
    ) -> tuple[jax.Array, jax.Array]:
        """The native surface: the filter integrals, onto the pivot wavelengths.

        The incoming coordinates are a bare array here, and the inherited
        :meth:`influence` needs an :class:`~ampere.core.Axis` — the lookup is
        the axis's job (``spectrum_photometry.md`` Gap 1: the columns go where
        ``Axis.locate`` reports, never where the tabulation happens to sit in a
        larger or reordered negotiated grid). **So the axis is rebuilt here
        rather than the lookup being skipped.**

        Slice 1 passed the bare array straight through and was simply broken —
        ``influence`` called ``axis.locate`` on an ``ndarray``, which has no
        such method, so every native photometry chain raised
        ``AttributeError``. It was invisible because no conformance
        ``ProblemSpec`` had a photometry step; slice 2 adds one (the battery's
        ``PHOTOMETRIC`` shape) and this method is what it exercises. The torch
        backend rebuilt the axis from the start, and this is that fix.
        """
        coordinates = np.asarray(grid, dtype=float)
        axis = Axis.build("spectral_axis", coordinates, COORDINATE_UNIT)
        return self._influence_jax(axis) @ flux, self.pivot_grid()

    def apply(self, samples: Any, values: Any) -> PhotometricPoints:
        weights = np.asarray(self.influence(samples.spectral_axis), dtype=float)
        integrated = jnp.asarray(weights, dtype=jnp.float64) @ jnp.asarray(
            samples.values, dtype=jnp.float64
        )
        return PhotometricPoints(
            # `filters` is a property on the reference class, not a method;
            # calling it raised `'tuple' object is not callable`. Same cause as
            # `apply_flux`'s own defect above -- nothing composed this step
            # into a problem, so neither surface had ever been evaluated.
            self.filters,
            self.pivots() * COORDINATE_UNIT,
            np.asarray(integrated),
            unit=samples.unit,
            mask=propagate_mask(samples, weights),
        )
