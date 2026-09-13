"""Interferometric visibilities and closure phases on the jax path.

The native twins of :mod:`ampere.backends.reference.interferometry` (W4.3), as
:mod:`ampere.backends.torch.interferometry` is on the other modern backend:
the same five steps and the same six source models, with the arithmetic in
``jax.numpy`` and the geometry unchanged.

Why these subclass the reference classes
----------------------------------------
The reason this module's siblings give (:mod:`ampere.backends.jax.instrument`),
and one more that is specific to this modality. Each class here derives from
its ``ampere.backends.reference`` counterpart and overrides exactly two things:
the four capability flags, and the arithmetic. Everything else — the coverage
buffers and ``from_observed``, the published field-of-view and Nyquist
requirements, ``configure_from``'s expansion product, the canonical-ordering
check, ``compile_for``'s adoption of the negotiated grid and its refusal when a
model holds its own — is inherited, because none of it is arithmetic a gradient
passes through. It is geometry: numpy data preparation done once at composition
time, off the hot loop.

The modality-specific reason is that here the declaration is the *dangerous*
half. ``FourierSample`` publishes ``max_step = 1/(2 s u_max)`` on both image
axes, and an image sampled more coarsely than that **aliases**, folding power
from beyond the limit onto the sampled baselines where it is indistinguishable
from real source structure (gap I-4). Two backends that wrote that requirement
twice would eventually negotiate slightly different grids for one declaration,
and the conformance suite would then be comparing two different problems and
calling the disagreement a numerical one. ``ampere.backends.torch.interferometry``
makes the same judgement for the same reason, so the two modern backends share
one declaration here and fork only the sums.

Two evaluation surfaces, as everywhere in this backend
------------------------------------------------------
``apply`` is the contract surface and returns an ``ampere.core`` container, so a
gradient stops at the container boundary; it is overridden here rather than
inherited so that the conformance rows — which drive ``apply`` — measure *this*
backend's transform against the closed forms in ``tests/conformance/oracles.py``
rather than numpy against numpy. ``apply_flux`` is the native surface:
coordinates and values in, the transformed pair out, pure and traceable, and
what :mod:`ampere.backends.jax.problem` composes into the differentiable
log-density.

The models' native surface is spelled ``native_flux`` / ``native_grid`` rather
than ``flux`` / ``grid``, and not by preference: every source model here
declares a *parameter* called ``flux`` — the source's total flux density, which
is the thing one fits — and ``Parameterised._check_free_name`` refuses a
parameter whose name shadows a class attribute. :mod:`ampere.backends.jax.problem`
accepts either spelling.

The transform as a native operator
----------------------------------
The reference backend's separable form, ported: the DFT kernel factorises as
``exp(-2 pi i u x) exp(-2 pi i v y)``, so one ``(n_uv, ny)`` contraction and one
``(n_uv, nx)`` row product compute the same sum as an ``(n_uv, nx, ny)`` array
would. Both phase matrices and the cell solid angles are functions of the
negotiated grid and the observed coverage alone — nothing a sampler proposes —
so each is built once, placed on this piece's device, and cached against the
grid's exact bytes. The contraction carries a leading ellipsis, so the same
expression evaluates one image or a ``vmap``\\ ped stack of them.

Complex128 throughout, which is why every class here calls ``require_x64``: a
visibility is complex, and ``complex64`` would make the phase of a long
baseline meaningless.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any, ClassVar

import jax
import jax.numpy as jnp
import numpy as np
from jax.scipy.special import bessel_jn

from ampere.backends.reference.interferometry import (
    _FWHM_PER_SIGMA as FWHM_PER_SIGMA,  # the same constant, not a second copy
)
from ampere.backends.reference.interferometry import (
    MAS_PER_RAD,
    SPECTRAL_UNIT,
    cell_solid_angle,
)
from ampere.backends.reference.interferometry import (
    Amplitude as _ReferenceAmplitude,
)
from ampere.backends.reference.interferometry import (
    BandwidthSmearing as _ReferenceBandwidthSmearing,
)
from ampere.backends.reference.interferometry import (
    Binary as _ReferenceBinary,
)
from ampere.backends.reference.interferometry import (
    BinaryVisibilities as _ReferenceBinaryVisibilities,
)
from ampere.backends.reference.interferometry import (
    ClosurePhase as _ReferenceClosurePhase,
)
from ampere.backends.reference.interferometry import (
    FourierSample as _ReferenceFourierSample,
)
from ampere.backends.reference.interferometry import (
    GaussianSource as _ReferenceGaussianSource,
)
from ampere.backends.reference.interferometry import (
    GaussianSourceVisibilities as _ReferenceGaussianSourceVisibilities,
)
from ampere.backends.reference.interferometry import (
    TimeSmearing as _ReferenceTimeSmearing,
)
from ampere.backends.reference.interferometry import (
    UniformDisc as _ReferenceUniformDisc,
)
from ampere.backends.reference.interferometry import (
    UniformDiscVisibilities as _ReferenceUniformDiscVisibilities,
)
from ampere.core import (
    ChannelRequirements,
    ClosurePhases,
    Model,
    Transformation,
    TransformationError,
    VisibilitySet,
    propagate_mask,
)

from ._config import BACKEND, require_x64
from ._device import DEVICE, device_flag, place_on, resolve_device

__all__ = [
    "Amplitude",
    "BandwidthSmearing",
    "Binary",
    "BinaryVisibilities",
    "ClosurePhase",
    "FourierSample",
    "GaussianSource",
    "GaussianSourceVisibilities",
    "TimeSmearing",
    "UniformDisc",
    "UniformDiscVisibilities",
]

#: The three coordinate arrays a ``VisibilitySet``-valued step hands on as its
#: "grid": ``(u, v, wavelength)``. A one-axis kind's grid is one array; a
#: three-axis kind's is three, and the realisation threads whatever it is given
#: from step to step without looking inside it.
Coverage = tuple[jax.Array, jax.Array, jax.Array]

#: Below this argument ``2 J1(z)/z`` is taken from its ascending series rather
#: than from ``jax.scipy.special.bessel_jn``, which returns a NaN as ``z``
#: approaches zero (its downward recurrence normalises by a quantity that
#: underflows). Four series terms are good to 7e-16 over ``[0, 0.1]``, measured.
_SMALL_DISC_ARGUMENT = 0.1


# ---------------------------------------------------------------------------
# Shared plumbing
# ---------------------------------------------------------------------------


class _JaxInterferometryStep:
    """The four capability flags, this backend's placement, and the caches.

    A **mixin** rather than a base class, so each concrete step still derives
    from its reference counterpart and inherits that class's declaration
    wholesale — exactly the shape :class:`ampere.backends.jax.instrument._JaxStep`
    has, and for the same reasons.
    """

    #: The gradient the whole backend exists for passes through
    #: :meth:`apply_flux`; the declaration surface is the reference class's and
    #: carries no gradient at all.
    DIFFERENTIABLE: ClassVar[bool] = True
    #: ``True``, and measured rather than asserted: the batched rows compare a
    #: ``jax.vmap``\\ ped forward model against the same model in a loop. Every
    #: expression below is whole-array ``jax.numpy`` with a leading ellipsis, so
    #: nothing branches on a value and nothing indexes by one.
    BATCHABLE: ClassVar[bool] = True
    #: Never auto-detected (``architecture.md`` §5); per instance, so a step
    #: placed on an accelerator beside CPU models is a device disagreement
    #: ``ampere.core.declared_capabilities`` refuses at composition.
    DEVICE: ClassVar[str] = DEVICE
    BACKEND: ClassVar[str] = BACKEND

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        require_x64(f"a jax {type(self).__name__}")
        # Popped before the reference constructor sees it: ``device=`` is this
        # backend's placement keyword and means nothing to the numpy step whose
        # declaration is inherited wholesale.
        device = kwargs.pop("device", DEVICE)
        super().__init__(*args, **kwargs)
        resolved = resolve_device(device, f"a jax {type(self).__name__}", error=TransformationError)
        self._device = resolved
        object.__setattr__(self, "DEVICE", device_flag(device, resolved))

    def _place(self, values: Any) -> jax.Array:
        """*values* as a float64 jax array on this step's device."""
        return place_on(
            jnp.asarray(np.asarray(values, dtype=float)), getattr(self, "_device", None)
        )

    def _observed_array(self, values: Any) -> jax.Array:
        """A container's own values, complex if the container's are.

        The container decides, not this step: a ``VisibilitySet`` is legally
        real (after :class:`Amplitude`) or complex, and coercing one to the
        other would be either a lost imaginary part or an invented zero.
        """
        array = np.asarray(values)
        dtype = jnp.complex128 if array.dtype.kind == "c" else jnp.float64
        return place_on(jnp.asarray(array, dtype=dtype), getattr(self, "_device", None))


def _unit_phase(angle: jax.Array) -> jax.Array:
    """``exp(i * angle)`` as an explicitly complex128 array.

    Written as ``cos + i sin`` rather than ``jnp.exp(1j * angle)`` so that a
    float64 array becomes complex128 by construction rather than by jax's
    promotion of a Python complex scalar — which is correct only while x64 is
    on, and a silently complex64 phase factor would be a precision loss nothing
    in this module would report.
    """
    return jnp.cos(angle) + 1j * jnp.sin(angle)


def _dense_mask(samples: Any, n_out: int) -> Any:
    """``propagate_mask`` with a Fourier transform's block-of-ones influence.

    A DFT is dense: one masked pixel touches every sample, so the honest
    influence matrix is a block of ones and the honest consequence — a single
    masked pixel masks the whole visibility set — is correct rather than
    convenient. ``None`` for an unmasked image, so the usual case pays nothing.
    """
    if samples.mask is None:
        return None
    return propagate_mask(samples, np.ones((n_out, np.asarray(samples.values).size), dtype=float))


# ---------------------------------------------------------------------------
# The steps
# ---------------------------------------------------------------------------


class FourierSample(_JaxInterferometryStep, _ReferenceFourierSample):
    """Sample a model image at the observed ``(u, v)`` points, in jax.

    ``Image -> VisibilitySet``: the kind-changing, coordinate-changing step
    whose derivative is what makes a gradient-based fit of a sky model possible.
    A direct discrete Fourier transform — no FFT, no gridding, no
    interpolation — in the reference backend's separable form::

        V_n = sum_jk I_jk dOmega_jk exp(-2 pi i (u_n x_j + v_n y_k))

    Parameters
    ----------
    u_pts, v_pts, wavelength, field_of_view, oversampling, label
        As :class:`ampere.backends.reference.interferometry.FourierSample`.
    device
        Where this step's constant arrays live. Never detected.
    """

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        self._coverage_cache: Coverage | None = None
        self._operator_cache: tuple[bytes, tuple[jax.Array, jax.Array, jax.Array]] | None = None

    # -- chain-internal negotiation ------------------------------------------

    def configure_from(self, downstream: Sequence[Transformation]) -> None:
        """The reference's expansion product, then drop the derived arrays."""
        super().configure_from(downstream)
        self._coverage_cache = None
        self._operator_cache = None

    # -- the operator, built once --------------------------------------------

    def coverage_arrays(self) -> Coverage:
        """The expanded ``(u, v, lambda)`` this step evaluates, as jax arrays."""
        found = self._coverage_cache
        if found is not None:
            return found
        u_pts, v_pts, waves = self.expanded_coverage
        built = (self._place(u_pts), self._place(v_pts), self._place(waves))
        self._coverage_cache = built
        return built

    def _operator(self, x_mas: Any, y_mas: Any) -> tuple[jax.Array, jax.Array, jax.Array]:
        """``(phase_x, phase_y, solid_angle)`` for this image grid, cached.

        Keyed on the coordinates' exact bytes, for the reason
        ``instrument._JaxStep._influence_jax`` gives: the operator is a pure
        function of the grid, the grid is fixed for a run, and rebuilding an
        ``(n_uv, nx)`` complex matrix per likelihood evaluation is the hot-loop
        waste W2.3 measured. An exact key means a changed grid rebuilds and a
        repeated one does not, so a lookup can never return the operator of a
        nearby grid.
        """
        x_array = np.asarray(x_mas, dtype=float)
        y_array = np.asarray(y_mas, dtype=float)
        key = x_array.tobytes() + b"|" + y_array.tobytes()
        cached = self._operator_cache
        if cached is not None and cached[0] == key:
            return cached[1]
        u_pts, v_pts, _ = self.coverage_arrays()
        # The solid angles are the reference quadrature's, by value: cell widths
        # from the midpoints of whatever grid negotiation produced. Geometry
        # computed once at composition time, which is the same judgement the
        # other steps in this backend make about their influence matrices.
        omega = self._place(cell_solid_angle(x_array, y_array))
        phase_x = _unit_phase(-2.0 * math.pi * jnp.outer(u_pts, self._place(x_array / MAS_PER_RAD)))
        phase_y = _unit_phase(-2.0 * math.pi * jnp.outer(v_pts, self._place(y_array / MAS_PER_RAD)))
        built = (phase_x, phase_y, omega)
        self._operator_cache = (key, built)
        return built

    # -- the native surface a realisation composes ---------------------------

    def apply_flux(self, flux: jax.Array, grid: Any, values: Any) -> tuple[jax.Array, Any]:
        """The DFT as one batched contraction: ``(..., nx, ny) -> (..., n_uv)``.

        *grid* is the image model's ``(x, y)`` pair in mas; what comes back is
        the ``(u, v, lambda)`` triple of the emitted visibilities, because this
        step changes the coordinates as well as the kind.
        """
        x_mas, y_mas = grid
        phase_x, phase_y, omega = self._operator(x_mas, y_mas)
        weighted = (flux * omega).astype(jnp.complex128)
        along_y = jnp.einsum("nk,...jk->...nj", phase_y, weighted)
        return jnp.einsum("nj,...nj->...n", phase_x, along_y), self.coverage_arrays()

    # -- the contract surface ------------------------------------------------

    def apply(self, samples: Any, values: Any) -> VisibilitySet:
        """The same arithmetic, back on the numpy side of the container boundary."""
        grid = (self._place(samples.x.values), self._place(samples.y.values))
        visibility, _ = self.apply_flux(self._place(samples.values), grid, values)
        template = self._resolve_template(samples.unit, *self.expanded_coverage)
        return template.with_values(
            np.asarray(visibility), mask=_dense_mask(samples, self.expanded_coverage[0].size)
        )


class ClosurePhase(_JaxInterferometryStep, _ReferenceClosurePhase):
    """Three visibilities to one angle, in jax: ``VisibilitySet -> ClosurePhases``.

    ``arg(V_ij V_jk V_ki)`` for a triangle whose three baselines occupy
    positions ``3t``, ``3t + 1``, ``3t + 2`` — the canonical ordering
    :class:`~ampere.core.ClosurePhases` fixes and
    :meth:`FourierSample.from_observed` produces. ``jnp.angle`` of a complex
    product is differentiable and maps under ``vmap``, which is what lets the
    observable that survives the atmosphere carry a gradient.
    """

    def apply_flux(self, flux: jax.Array, grid: Any, values: Any) -> tuple[jax.Array, Any]:
        """``(..., 3n) -> (..., n)``: the argument of each triangle's product."""
        total = int(flux.shape[-1])
        if total % 3:
            raise TransformationError(
                f"ClosurePhase reads three visibilities per triangle, so it needs a multiple of "
                f"three, got {total}. Build the Fourier step before it with "
                f"FourierSample.from_observed(closure_phases, field_of_view=...), which lays the "
                f"three baselines of each triangle out in the canonical order."
            )
        n_out = total // 3
        triples = flux.reshape(*flux.shape[:-1], n_out, 3)
        product = triples[..., 0] * triples[..., 1] * triples[..., 2]
        coverage = tuple(axis.reshape(n_out, 3)[:, 0] for axis in grid)
        return jnp.angle(product), coverage

    def apply(self, samples: Any, values: Any) -> ClosurePhases:
        """The same arithmetic, with the reference class's template and mask."""
        coverage = tuple(
            self._place(samples.axis(name).values) for name in ("u", "v", "spectral_axis")
        )
        phase, _ = self.apply_flux(self._observed_array(samples.values), coverage, values)
        n_out = int(samples.n_samples) // 3
        template = self._resolve_template(samples, n_out)
        mask = propagate_mask(samples, self._influence(n_out) if samples.mask is not None else None)
        return template.with_values(np.asarray(phase), mask=mask)


class Amplitude(_JaxInterferometryStep, _ReferenceAmplitude):
    """Take the modulus, in jax: a complex ``VisibilitySet`` to a real one.

    Kind-preserving, coordinate-preserving, parameter-free — a real scientific
    decision made as a named step rather than as a convention buried in a
    family (``interferometry.md`` §6).
    """

    def apply_flux(self, flux: jax.Array, grid: Any, values: Any) -> tuple[jax.Array, Any]:
        """``|V|``, with the coordinates untouched."""
        return jnp.abs(flux), grid

    def apply(self, samples: Any, values: Any) -> VisibilitySet:
        return samples.with_values(np.asarray(jnp.abs(self._observed_array(samples.values))))


class _JaxAveragingStep(_JaxInterferometryStep):
    """The block reduction both smearing steps share, in jax.

    The geometry — which extra ``(u, v)`` samples to average over
    (``expand_uv``), and how many sub-samples the steps after this one added
    (``configure_from``) — is the reference class's, read by
    :class:`FourierSample` at composition time. What is here is the weighted
    mean, as one ``einsum`` over the flat layout that step fixes.
    """

    def _block(self, total: int) -> tuple[int, int, int]:
        """``(n_out, count, trailing)``, with the reference refusal by name."""
        count = int(self._weights.size)
        trailing = int(self._trailing)
        if total % (count * trailing):
            raise TransformationError(
                f"{type(self).__name__} was given {total} visibilities, which is not a multiple "
                f"of the {count} sub-sample(s) it averages over times the {trailing} the steps "
                f"after it added. That happens when the chain's Fourier step was not the one "
                f"configured for this chain: build the step and the instrument together, so "
                f"configure_from runs on the pair ampere will evaluate."
            )
        return total // count, count, trailing

    def apply_flux(self, flux: jax.Array, grid: Any, values: Any) -> tuple[jax.Array, Any]:
        """The weighted mean of each block, and the un-smeared coordinates.

        Index ``0`` along every sub-axis is the input sample itself, so the
        output coordinates come from *slicing* the input rather than from
        recomputing them — which is what keeps the bit-identical coordinate rule
        on the very container a likelihood compares against the data.
        """
        n_out, count, trailing = self._block(int(flux.shape[-1]))
        lead = tuple(flux.shape[:-1])
        shaped = flux.reshape(*lead, -1, count, trailing)
        weights = self._place(self._weights).astype(shaped.dtype)
        averaged = jnp.einsum("a,...nar->...nr", weights, shaped)
        coverage = tuple(axis.reshape(-1, count, trailing)[:, 0, :].reshape(n_out) for axis in grid)
        return averaged.reshape(*lead, n_out), coverage

    def apply(self, samples: Any, values: Any) -> VisibilitySet:
        """The same average, with the reference class's template and mask rule."""
        coverage = tuple(
            self._place(samples.axis(name).values) for name in ("u", "v", "spectral_axis")
        )
        averaged, reduced = self.apply_flux(self._observed_array(samples.values), coverage, values)
        n_out, count, trailing = self._block(int(samples.n_samples))
        mask = None
        if samples.mask is not None:
            mask = propagate_mask(samples, self._influence(n_out, count, trailing))
        return VisibilitySet(
            np.asarray(reduced[0]),
            np.asarray(reduced[1]),
            np.asarray(reduced[2]) * SPECTRAL_UNIT,
            np.asarray(averaged),
            unit=samples.unit,
            mask=mask,
        )


class BandwidthSmearing(_JaxAveragingStep, _ReferenceBandwidthSmearing):
    """Average the visibility across a spectral channel of finite width, in jax.

    The baseline is fixed in metres, so the sampled spatial frequency
    ``B/lambda`` sweeps radially as ``lambda`` crosses the channel. The effect
    reduces the amplitude of a resolved source at long baselines and looks
    exactly like a larger source if it is not modelled.

    Parameters
    ----------
    resolving_power, nodes, label
        As the reference step.
    device
        Where this step's constant arrays live.
    """


class TimeSmearing(_JaxAveragingStep, _ReferenceTimeSmearing):
    """Average the visibility over a finite integration, along the uv track, in jax.

    Unlike bandwidth smearing the direction is not radial — it is wherever the
    track goes. The uv rates are per-sample buffers given at construction or
    read from the observed container's ``extra_coords`` by the inherited
    ``from_observed``.

    Parameters
    ----------
    integration, du_dt, dv_dt, nodes, label
        As the reference step.
    device
        Where this step's constant arrays live.
    """


# ---------------------------------------------------------------------------
# The image-emitting models
# ---------------------------------------------------------------------------


class _JaxImageModel:
    """Shared plumbing for the image-emitting trio, in jax.

    A mixin over the reference model, for the reason
    :class:`_JaxInterferometryStep` is one. Each concrete model writes one
    method, :meth:`_brightness_native`, returning a surface brightness in Jy/sr
    on the ``(x, y)`` arrays it is handed, in mas; :meth:`_brightness` routes the
    contract path through the same expression, so ``evaluate`` and
    :meth:`native_flux` cannot disagree.
    """

    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = True
    DEVICE: ClassVar[str] = DEVICE
    BACKEND: ClassVar[str] = BACKEND

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        require_x64(f"a jax {type(self).__name__}")
        device = kwargs.pop("device", DEVICE)
        super().__init__(*args, **kwargs)
        resolved = resolve_device(device, f"a jax {type(self).__name__}", error=TransformationError)
        object.__setattr__(self, "_device", resolved)
        object.__setattr__(self, "DEVICE", device_flag(device, resolved))
        self._grids: dict[str, tuple[jax.Array, jax.Array]] = {}
        self._own_grid = tuple(self._place(self._data(name)) for name in ("x", "y"))

    def _place(self, values: Any) -> jax.Array:
        return place_on(
            jnp.asarray(np.asarray(values, dtype=float)), getattr(self, "_device", None)
        )

    # -- negotiation ---------------------------------------------------------

    def compile_for(self, requirements: Mapping[str, ChannelRequirements]) -> Model:
        """Adopt the negotiated image grid, then cache it as jax arrays."""
        super().compile_for(requirements)
        for channel, template in self.templates.items():
            self._grids[channel] = (
                self._place(template.x.values),
                self._place(template.y.values),
            )
        return self

    # -- the native surface a realisation composes ---------------------------

    def native_grid(self, channel: str) -> tuple[jax.Array, jax.Array]:
        """The ``(x, y)`` arrays this channel evaluates on, mas.

        Two arrays rather than one, because an image lives on two axes; the
        realisation threads whatever this returns through the chain's
        ``apply_flux`` calls without looking inside it. ``native_*`` rather than
        ``grid``/``flux`` because ``flux`` is one of these models' own
        parameters — see the module docstring.
        """
        found = self._grids.get(channel)
        if found is not None:
            return found
        return (self._own_grid[0], self._own_grid[1])

    def native_flux(self, channel: str, values: Mapping[str, Any] | None = None) -> jax.Array:
        """This model's surface brightness on *channel*, Jy/sr, as a jax array."""
        x_mas, y_mas = self.native_grid(channel)
        return self._brightness_native(x_mas, y_mas, self.context(values))

    # -- the contract surface ------------------------------------------------

    def _brightness(
        self, x_mas: np.ndarray, y_mas: np.ndarray, context: Mapping[str, Any]
    ) -> np.ndarray:
        """The jax brightness, back on the numpy side of the container boundary."""
        self._validate(context)
        return np.asarray(self._brightness_native(self._place(x_mas), self._place(y_mas), context))

    def _validate(self, context: Mapping[str, Any]) -> None:
        """Refuse a geometrically impossible shape, on the **contract** path only.

        The reference models raise for a non-positive diameter or width. Such a
        refusal cannot live in :meth:`_brightness_native`: under a trace the
        value is a ``Tracer`` and a Python ``if`` on one is an error, and
        ``inference.md`` §10a's rule is that nothing on a native path branches
        on a value. So the check stays where the values are plain numbers, and
        the native path relies on the prior's own support instead.
        """

    def _brightness_native(
        self, x_mas: jax.Array, y_mas: jax.Array, context: Mapping[str, Any]
    ) -> jax.Array:
        raise NotImplementedError


def _gaussian_brightness(
    x_mas: jax.Array,
    y_mas: jax.Array,
    centre_x: Any,
    centre_y: Any,
    sigma_mas: Any,
    flux: Any,
) -> jax.Array:
    """A circular Gaussian of total *flux* (Jy) as a surface brightness in Jy/sr.

    The reference expression, in jax: the outer sum is formed by indexing the
    two coordinate axes, which are never batched (they are the negotiated
    grid), while every parameter may be — so one expression serves one theta and
    a ``vmap``\\ ped stack of them.
    """
    sigma_rad = sigma_mas / MAS_PER_RAD
    squared = ((x_mas[:, None] - centre_x) ** 2 + (y_mas[None, :] - centre_y) ** 2) / sigma_mas**2
    return flux / (2.0 * math.pi * sigma_rad**2) * jnp.exp(-0.5 * squared)


class UniformDisc(_JaxImageModel, _ReferenceUniformDisc):
    """A uniformly bright circular disc, in jax — **not** differentiable.

    ``flux`` is the total flux density in Jy and ``diameter`` the angular
    diameter in mas; the brightness inside the disc is ``flux / (pi
    (theta/2)**2)`` in Jy/sr and zero outside.

    ``DIFFERENTIABLE = False``, declared rather than inherited, and the reason
    is mathematical rather than a missing implementation. A hard-edged disc
    rendered onto a grid is a **piecewise-constant** function of its diameter:
    the ``where`` moves the edge only when a pixel centre crosses it, so
    automatic differentiation returns the derivative of the ``1/r**2`` amplitude
    alone and silently omits the boundary term, which is the larger half. That
    is not a small error in a gradient — it points the wrong way over most of
    the parameter range — and a sampler handed it would converge confidently to
    the wrong diameter. Fit a disc on a gradient-free engine, or fit
    :class:`GaussianSource`, whose gradient is exact. ``BATCHABLE`` is ``False``
    for the same reason and not a second one: a model that declares no
    derivative has no realised density to ``vmap``.

    Parameters
    ----------
    x, y, diameter, flux, channels, adopt_grid
        As the reference model.
    device
        Where this model's grids live.
    """

    DIFFERENTIABLE: ClassVar[bool] = False
    BATCHABLE: ClassVar[bool] = False

    def _validate(self, context: Mapping[str, Any]) -> None:
        radius = 0.5 * float(context["diameter"])
        if radius <= 0.0:
            raise ValueError(f"a uniform disc needs a positive diameter, got {radius * 2.0!r} mas.")

    def _brightness_native(
        self, x_mas: jax.Array, y_mas: jax.Array, context: Mapping[str, Any]
    ) -> jax.Array:
        radius = 0.5 * context["diameter"]
        distance = jnp.hypot(x_mas[:, None], y_mas[None, :])
        area = math.pi * (radius / MAS_PER_RAD) ** 2
        return jnp.where(distance <= radius, context["flux"] / area, 0.0)


class GaussianSource(_JaxImageModel, _ReferenceGaussianSource):
    """A circular Gaussian source, in jax. Band-limited, and differentiable.

    ``flux`` is the total flux density in Jy and ``fwhm`` the full width at half
    maximum in mas. Band-limited, so the direct DFT of its image agrees with its
    closed-form visibility to the solver tolerance on any grid that satisfies
    the Nyquist requirement and covers the source.

    Parameters
    ----------
    x, y, fwhm, flux, channels, adopt_grid
        As the reference model.
    device
        Where this model's grids live.
    """

    def _validate(self, context: Mapping[str, Any]) -> None:
        if float(context["fwhm"]) / FWHM_PER_SIGMA <= 0.0:
            raise ValueError(f"a Gaussian source needs a positive fwhm, got {context['fwhm']!r}.")

    def _brightness_native(
        self, x_mas: jax.Array, y_mas: jax.Array, context: Mapping[str, Any]
    ) -> jax.Array:
        return _gaussian_brightness(
            x_mas, y_mas, 0.0, 0.0, context["fwhm"] / FWHM_PER_SIGMA, context["flux"]
        )


class Binary(_JaxImageModel, _ReferenceBinary):
    """Two Gaussian components, in jax — the source the phase is proved on.

    The primary sits at the phase centre and the secondary at ``separation`` and
    ``position_angle`` (radians, east of north); ``flux_ratio`` is the
    secondary's flux over the primary's and ``flux`` is the pair's total. Both
    components share one ``component_fwhm``, a buffer rather than a parameter:
    it is there to make the source band-limited, and one would not put a prior
    on a numerical device.

    This is the model the NUTS rows fit: its closure phases are non-zero,
    analytic and sensitive to exactly the sign conventions that would otherwise
    go unnoticed, and every one of its four parameters enters the brightness
    through ``sin``, ``cos`` and ``exp``, so the gradient is exact rather than a
    discretisation of one.

    Parameters
    ----------
    x, y, separation, position_angle, flux_ratio, flux, component_fwhm,
    channels, adopt_grid
        As the reference model.
    device
        Where this model's grids live.
    """

    def _brightness_native(
        self, x_mas: jax.Array, y_mas: jax.Array, context: Mapping[str, Any]
    ) -> jax.Array:
        separation = context["separation"]
        angle = context["position_angle"]
        offset_x = separation * jnp.sin(angle)
        offset_y = separation * jnp.cos(angle)
        ratio = context["flux_ratio"]
        sigma = context["component_fwhm"] / FWHM_PER_SIGMA
        primary = context["flux"] / (1.0 + ratio)
        return _gaussian_brightness(x_mas, y_mas, 0.0, 0.0, sigma, primary) + _gaussian_brightness(
            x_mas, y_mas, offset_x, offset_y, sigma, primary * ratio
        )


# ---------------------------------------------------------------------------
# The analytic route: the same three sources, emitting visibilities directly
# ---------------------------------------------------------------------------


class _JaxVisibilityModel:
    """Shared plumbing for the trio that skips the Fourier step, in jax.

    ``interferometry.md`` §1: a model that computes visibilities analytically
    emits a :class:`~ampere.core.VisibilitySet` channel directly, and its
    instrument is the no-steps one. Both routes are supported and neither is
    privileged — and here each is the other's oracle, which is why these are
    worth twinning even though nothing fits with them: a closed form in jax is
    what says this backend's direct transform is right for a reason other than
    "numpy agrees with numpy".
    """

    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = True
    DEVICE: ClassVar[str] = DEVICE
    BACKEND: ClassVar[str] = BACKEND

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        require_x64(f"a jax {type(self).__name__}")
        device = kwargs.pop("device", DEVICE)
        super().__init__(*args, **kwargs)
        resolved = resolve_device(device, f"a jax {type(self).__name__}", error=TransformationError)
        object.__setattr__(self, "_device", resolved)
        object.__setattr__(self, "DEVICE", device_flag(device, resolved))
        self._coverage = tuple(
            self._place(self.buffers[name].value) for name in ("u_pts", "v_pts", "wavelength")
        )

    def _place(self, values: Any) -> jax.Array:
        return place_on(
            jnp.asarray(np.asarray(values, dtype=float)), getattr(self, "_device", None)
        )

    # -- the native surface a realisation composes ---------------------------

    def native_grid(self, channel: str) -> Coverage:
        """The ``(u, v, lambda)`` arrays this model already emits on."""
        return (self._coverage[0], self._coverage[1], self._coverage[2])

    def native_flux(self, channel: str, values: Mapping[str, Any] | None = None) -> jax.Array:
        """This model's complex visibilities on *channel*, Jy, as a jax array."""
        u_pts, v_pts, _ = self.native_grid(channel)
        return self._visibility_native(jnp.hypot(u_pts, v_pts), self.context(values))

    # -- the contract surface ------------------------------------------------

    def _visibility(self, rho: np.ndarray, context: Mapping[str, Any]) -> np.ndarray:
        """The jax closed form, back on the numpy side of the boundary."""
        return np.asarray(self._visibility_native(self._place(rho), context))

    def _visibility_native(self, rho: jax.Array, context: Mapping[str, Any]) -> jax.Array:
        raise NotImplementedError


class UniformDiscVisibilities(_JaxVisibilityModel, _ReferenceUniformDiscVisibilities):
    """``V = flux * 2 J1(pi theta rho) / (pi theta rho)`` in jax — **not** differentiable.

    The disc's closed form, real-valued and negative beyond the first null,
    which is the sign flip an amplitude-only fit throws away. It is here as the
    conformance battery's jax-side oracle and as the contract-path closed form,
    and ``DIFFERENTIABLE = False`` for a **library** reason worth recording
    rather than hiding, since it is the only place in this module where jax's
    own special functions are not good enough.

    ``jax.scipy.special.bessel_jn`` evaluates ``J_v`` by a downward recurrence
    with a fixed, static ``n_iter``, and neither the default nor any other fixed
    choice is usable across a realistic argument range. Measured on jax 0.11.1,
    float64, against ``scipy.special.j1``:

    =============  ========================  =======================
    ``n_iter``     worst error, ``z <= 6``   worst error, ``z <= 200``
    =============  ========================  =======================
    default        7e-16                     3.5e-2
    100            7e-16                     1.5e-2
    400            2e-16 (but NaN for most   2e-16
                   ``z`` below about 20)
    =============  ========================  =======================

    and the default returns a **NaN** as ``z`` approaches zero. So the small-
    argument limit is taken from the ascending series here (four terms, good to
    7e-16 below :data:`_SMALL_DISC_ARGUMENT`), and the rest is ``bessel_jn``'s,
    accurate to 7e-16 for ``pi theta rho <= 6`` — which covers a source
    comparable to the beam, and is the regime the conformance rows live in —
    and degrading beyond it. Writing a Bessel function with a documented
    accuracy over the whole range would make its accuracy ampere's problem
    rather than a library's, which is the same judgement ``_families.py`` makes
    about the incomplete beta. The remedies, in order: the image route
    (:class:`UniformDisc` with :class:`FourierSample`, which needs no Bessel
    function at all), or the reference backend, whose ``scipy.special.j1`` is
    exact everywhere.

    Parameters
    ----------
    u_pts, v_pts, wavelength, diameter, flux, channels
        As the reference model.
    device
        Where this model's coverage buffers live.
    """

    DIFFERENTIABLE: ClassVar[bool] = False
    BATCHABLE: ClassVar[bool] = False

    def _visibility_native(self, rho: jax.Array, context: Mapping[str, Any]) -> jax.Array:
        argument = math.pi * context["diameter"] / MAS_PER_RAD * rho
        small = argument < _SMALL_DISC_ARGUMENT
        # Both branches are evaluated, so each needs an argument the *other*
        # branch's function is finite at: jnp.where propagates a NaN from the
        # branch it discards.
        safe = jnp.where(small, jnp.ones_like(argument), argument)
        half = jnp.where(small, 0.5 * argument, jnp.zeros_like(argument)) ** 2
        series = 1.0 - half / 2.0 + half**2 / 12.0 - half**3 / 144.0 + half**4 / 2880.0
        envelope = jnp.where(small, series, 2.0 * bessel_jn(safe, v=1)[1] / safe)
        return (context["flux"] * envelope).astype(jnp.complex128)


class GaussianSourceVisibilities(_JaxVisibilityModel, _ReferenceGaussianSourceVisibilities):
    """``V = flux * exp(-2 pi**2 sigma**2 rho**2)`` — the Gaussian's closed form, in jax."""

    def _visibility_native(self, rho: jax.Array, context: Mapping[str, Any]) -> jax.Array:
        sigma = context["fwhm"] / FWHM_PER_SIGMA / MAS_PER_RAD
        return (context["flux"] * jnp.exp(-2.0 * math.pi**2 * sigma**2 * rho**2)).astype(
            jnp.complex128
        )


class BinaryVisibilities(_JaxVisibilityModel, _ReferenceBinaryVisibilities):
    """``V = (flux/(1+f)) G(rho) (1 + f exp(-2 pi i (u dx + v dy)))`` in jax.

    The primary at the phase centre and the secondary at ``(dx, dy)`` radians,
    the same geometry :class:`Binary` renders; ``G`` is the components' common
    Gaussian envelope. The exponent carries the ``exp(-2 pi i)`` sign
    convention, which is what makes the closure phases non-zero — and what the
    conformance rows check the sign of against a closed form rather than against
    this module.
    """

    def _visibility_native(self, rho: jax.Array, context: Mapping[str, Any]) -> jax.Array:
        separation = context["separation"]
        angle = context["position_angle"]
        offset_x = separation * jnp.sin(angle)
        offset_y = separation * jnp.cos(angle)
        sigma = context["component_fwhm"] / FWHM_PER_SIGMA / MAS_PER_RAD
        envelope = jnp.exp(-2.0 * math.pi**2 * sigma**2 * rho**2)
        angle_uv = (
            -2.0
            * math.pi
            * (context["u_pts"] * offset_x + context["v_pts"] * offset_y)
            / MAS_PER_RAD
        )
        ratio = context["flux_ratio"]
        return (context["flux"] / (1.0 + ratio) * envelope).astype(jnp.complex128) * (
            1.0 + ratio * _unit_phase(angle_uv)
        )
