"""Interferometric visibilities and closure phases on the torch path.

The native twins of :mod:`ampere.backends.reference.interferometry` (W4.3): the
same five steps and the same six source models, with every number a gradient
passes through computed in ``torch``. What that buys is the thing Phase 4 set
out to prove — a sky model fitted to visibilities and closure phases **under
NUTS**, which needs a derivative of a discrete Fourier transform with respect
to the source parameters, and under ``vmap``, which needs the same transform as
one batched contraction.

Which twin pattern, and why
---------------------------
There is no step registry (``phase4_placement_memo.md`` §1.2): a twin
associates with its reference by **the same class name in this backend's own
module**, the duck-typed surface :mod:`ampere.backends.torch.problem` composes
(``apply_flux`` on a step; values and coordinates on a model, under the
``native_flux`` / ``native_grid`` spelling these use because ``flux`` is one of
their own parameters — see :meth:`_TorchImageModel.native_grid`), and its
``BACKEND`` declaration. All three hold here. What the memo left open is whether these
follow :mod:`ampere.backends.torch.instrument`, where the torch steps
**re-declare** everything and are held to the reference by the conformance
battery, or :mod:`ampere.backends.jax.instrument`, where each step
**inherits** its reference counterpart and overrides only the four capability
flags and the arithmetic.

**These follow the inheriting pattern, on both modern backends**, and the
reason is specific to this modality rather than a preference. In an LSF or a
resampler the declaration is small and the arithmetic is the whole of the
step; here it is the other way round. :class:`FourierSample` publishes
``intervals = (-fov/2, +fov/2)`` and ``max_step = 1/(2 s u_max)`` over the
*expanded* coverage, and its ``configure_from`` computes the expansion product
that fixes the flat sub-sample layout two other steps read back. None of that
is arithmetic a gradient passes through: the coverage is a buffer taken from
the observed container and never recomputed (``transformations.md`` §10, gap
I-2), the requirement is computed once at composition time, and the canonical
closure check is a refusal. And the consequence of writing it twice is worse
here than anywhere else in ampere: two backends that negotiated a *slightly*
different pixel scale for one declaration would not look like a disagreement —
an under-sampled image aliases, and aliased visibilities look like real source
structure rather than like an error (gap I-4,
:meth:`FourierSample.requirements`). Sharing the declaration and forking the
arithmetic is what makes the cross-backend rows mean something, which is
exactly the argument :mod:`ampere.backends.jax.instrument` makes for itself.

What is forked is every per-evaluation number, and the conformance battery
holds all of it to **numpy and to closed forms** rather than to this module:
:meth:`~ampere.core.Transformation.apply` — the container surface the battery
drives — is overridden here too, so a row that compares a transform against
:mod:`tests.conformance.oracles` is comparing torch arithmetic against a
Bessel function, not numpy against numpy. :class:`TorchInterferometryStep`
keeps this backend's own plumbing (``dtype``, ``device``, the ``DEVICE``
capability flag, the ``tensors`` buffer home and ``to()``), so a step here
behaves like every other torch step under a device move.

The transform as a native operator
----------------------------------
The separable form the reference backend uses is what is ported: the DFT kernel
factorises as ``exp(-2 pi i u x) exp(-2 pi i v y)``, so one ``(n_uv, ny)``
contraction followed by one ``(n_uv, nx)`` row product computes the same sum as
an ``(n_uv, nx, ny)`` array would, with none of its memory. Both phase factors
and the cell solid angles are **built once** — they are functions of the
negotiated grid and of the observed coverage, neither of which a sampler
proposes — cached against the grid's exact bytes, and registered as torch
buffers so a device move takes them along. The contraction itself is written
with a leading ellipsis, so the *same* expression evaluates one image or a
stack of them: that is what ``BATCHABLE = True`` means here, and it is measured
(``simulate_many(native=True)`` against the loop) rather than asserted.

Complex tensors end to end: the visibilities are ``torch.complex128``, the
complex partner of this backend's float64 policy (``_config.complex_dtype``),
which is the dtype :mod:`ampere.backends.torch.problem` already keeps a
``VisibilitySet``'s residual in for the ``complex_gaussian`` family.

What is **not** here, for the reasons the reference module gives: no FFT, no
chromatic sky, no reader. One thing is new: :class:`UniformDisc` and this
backend's :class:`UniformDiscVisibilities` declare ``DIFFERENTIABLE = False``,
and say why on themselves.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any, ClassVar

import numpy as np
import torch

from ampere.backends.reference.interferometry import (
    _FWHM_PER_SIGMA as FWHM_PER_SIGMA,  # the same constant, not a second copy
)
from ampere.backends.reference.interferometry import (
    _DEFAULT_SMEARING_NODES,
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

from ._config import (
    BACKEND,
    DEFAULT_DEVICE,
    DEFAULT_DTYPE,
    as_tensor,
    complex_dtype,
    move,
    place,
    to_numpy,
)
from .parameters import LoweredParameters

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
    "TorchInterferometryStep",
    "UniformDisc",
    "UniformDiscVisibilities",
]

#: The three coordinate tensors a ``VisibilitySet``-valued step hands on as its
#: "grid": ``(u, v, wavelength)``. A one-axis kind's grid is one tensor
#: (:mod:`ampere.backends.torch.instrument`); a three-axis kind's is three, and
#: the realisation threads whatever it is given from step to step without
#: looking inside it.
Coverage = tuple[torch.Tensor, torch.Tensor, torch.Tensor]

#: What :meth:`FourierSample._operator` caches: the two phase matrices and the
#: cell solid angles, all three functions of the grid and the coverage alone.
Operator = tuple[torch.Tensor, torch.Tensor, torch.Tensor]


# ---------------------------------------------------------------------------
# Shared plumbing
# ---------------------------------------------------------------------------


class TorchInterferometryStep:
    """The four capability flags, this backend's placement, and the cache.

    A **mixin** rather than a base class, so that each concrete step still
    derives from its reference counterpart and inherits that class's whole
    declaration — the constructor validation, the published requirements,
    ``configure_from``, ``from_observed`` and the canonical-ordering check. See
    the module docstring for why that is the right split for this modality.

    What it adds is exactly what :class:`~ampere.backends.torch.instrument.TorchStep`
    adds to a spectral step: ``dtype`` and ``device`` chosen at construction and
    never detected (``architecture.md`` §5), the ``DEVICE`` capability flag as
    an *instance* attribute so ``ampere.core.declared_capabilities`` composes a
    whole problem on one named device, a
    :class:`~ampere.backends.torch.parameters.LoweredParameters` home for the
    constant tensors so ``.to()`` moves them, and one cache slot per derived
    operator.
    """

    #: The four capability flags (``inference.md`` §18, W2.12's fourth), stated
    #: rather than inherited — the reference classes these derive from declare
    #: ``BACKEND = "reference"`` and ``DIFFERENTIABLE = False``, and inheriting
    #: that would be this backend claiming to be another one.
    DIFFERENTIABLE: ClassVar[bool] = True
    #: ``True``, and measured: ``tests/conformance/test_inference.py``'s batched
    #: rows and ``simulate_many(native=True)`` on the two-dataset interferometric
    #: problem compare a ``vmap``\\ ped forward model against the same model in a
    #: loop. Every expression below is whole-tensor and carries a leading
    #: ellipsis, so nothing branches on a value and nothing indexes by one.
    BATCHABLE: ClassVar[bool] = True
    #: Class default only; ``__init__`` shadows it per instance.
    DEVICE: ClassVar[str] = "cpu"
    BACKEND: ClassVar[str] = BACKEND

    def __init__(
        self,
        *args: Any,
        dtype: torch.dtype = DEFAULT_DTYPE,
        device: torch.device = DEFAULT_DEVICE,
        **kwargs: Any,
    ) -> None:
        # The reference constructor first: it is what validates the arguments
        # and registers the declared buffers these tensors are built from.
        super().__init__(*args, **kwargs)
        place(self, dtype, device)
        #: The torch-side home for this step's constant arrays (``lowering.md``
        #: §7). Every declared buffer is paired with a torch buffer, as
        #: :meth:`~ampere.backends.torch.instrument.TorchStep._declare` pairs
        #: them: that is what makes ``.to(...)`` a real move — and what makes a
        #: device this process has no hardware for a refusal at ``.to()`` rather
        #: than at the first evaluation.
        self.tensors = LoweredParameters()
        for name in self.buffers.names:
            self.tensors.register_buffer(
                name,
                as_tensor(
                    np.asarray(self.buffers[name].value, dtype=float),
                    dtype=self.dtype,
                    device=self.device,
                ),
                persistent=True,
            )

    def to(self, *, dtype: Any = None, device: Any = None) -> TorchInterferometryStep:
        """Move this step's constant tensors, and re-declare where they live.

        ``architecture.md`` §5's "buffers move with parameters under
        ``.to(...)``". Every derived operator is discarded rather than moved:
        it is cheap to rebuild, it is keyed on the grid it was built for, and a
        half-moved cache is the kind of bug that shows up as a device error
        three calls later.
        """
        self._forget()
        return move(self, dtype=dtype, device=device)

    def _forget(self) -> None:
        """Drop every cached operator. Overridden where there is one to drop."""

    def _real(self, value: Any) -> torch.Tensor:
        """*value* as a real tensor where this step's arithmetic happens."""
        return as_tensor(value, dtype=self.dtype, device=self.device)

    def _complex(self, value: Any) -> torch.Tensor:
        """*value* as a complex tensor where this step's arithmetic happens."""
        return as_tensor(value, dtype=complex_dtype(self.dtype), device=self.device)

    def _observed_tensor(self, values: Any) -> torch.Tensor:
        """A container's own values as a tensor, complex if they are complex.

        The container decides, not this step: a ``VisibilitySet`` is legally
        real (after :class:`Amplitude`) or complex, and coercing one to the
        other would be either a lost imaginary part or an invented zero.
        """
        array = np.asarray(values)
        return self._complex(array) if array.dtype.kind == "c" else self._real(array)


def _unit_phase(angle: torch.Tensor) -> torch.Tensor:
    """``exp(i * angle)`` as an explicitly complex tensor.

    Written as ``cos + i sin`` rather than ``torch.exp(1j * angle)`` so that a
    real float64 tensor becomes complex128 by construction rather than by
    torch's promotion of a Python complex scalar — the promotion is correct
    today, but a silently float32-complex phase factor would be a precision
    loss nothing in this module would report.
    """
    return torch.complex(torch.cos(angle), torch.sin(angle))


def _dense_mask(samples: Any, n_out: int) -> Any:
    """``propagate_mask`` with a Fourier transform's block-of-ones influence.

    A DFT is dense: one masked pixel touches every sample, so the honest
    influence matrix is a block of ones and the honest consequence — a single
    masked pixel masks the whole visibility set — is correct rather than
    convenient. ``None`` for an unmasked image, so the usual case pays nothing
    (the reference module's ``_image_influence``, by value).
    """
    if samples.mask is None:
        return None
    return propagate_mask(samples, np.ones((n_out, np.asarray(samples.values).size), dtype=float))


# ---------------------------------------------------------------------------
# The steps
# ---------------------------------------------------------------------------


class FourierSample(TorchInterferometryStep, _ReferenceFourierSample):
    """Sample a model image at the observed ``(u, v)`` points, in torch.

    ``Image -> VisibilitySet``: the kind-changing, coordinate-changing step in
    the middle of the chain, and the one whose derivative makes a
    gradient-based fit of a sky model possible at all. A direct discrete
    Fourier transform — no FFT, no gridding, no interpolation — in the
    separable form the reference backend fixes::

        V_n = sum_jk I_jk dOmega_jk exp(-2 pi i (u_n x_j + v_n y_k))

    computed as ``einsum("nj,...nj->...n", phase_x, einsum("nk,...jk->...nj",
    phase_y, I dOmega))``. The two phase matrices and the solid angles are
    built once from the negotiated grid and the observed coverage and cached
    against the grid's exact bytes; the contraction is all that happens per
    evaluation.

    Everything declarative is the reference class's, unchanged and deliberately
    so: the ``(u, v, lambda)`` buffers and :meth:`from_observed`, the published
    field-of-view and Nyquist requirements, and ``configure_from``'s expansion
    product. See the module docstring.

    Parameters
    ----------
    u_pts, v_pts, wavelength, field_of_view, oversampling, label
        As :class:`ampere.backends.reference.interferometry.FourierSample`.
    dtype, device
        Where this step's constant tensors live, and in what precision. Never
        detected (``architecture.md`` §5).
    """

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)
        self._coverage_cache: Coverage | None = None
        self._operator_cache: tuple[bytes, Operator] | None = None

    # -- construction from the data ------------------------------------------

    @classmethod
    def from_observed(
        cls,
        container: Any,
        *,
        field_of_view: Any,
        oversampling: float = 1.0,
        label: str | None = None,
        **placement: Any,
    ) -> Any:
        """The reference route, with this backend's placement keywords forwarded.

        ``from_observed`` is **the supported route** — the coverage is taken
        from the observed container and never recomputed (gap I-2) — and the
        reference classmethod's signature predates a backend having a device to
        be placed on, so it would drop a ``device=`` silently. The coverage is
        still derived by the reference's own rule, through a probe, so the
        canonical ordering a ``ClosurePhases`` container unfolds into lives in
        exactly one place; this adds the placement and nothing else. (The same
        probe-and-rebuild shape the reference ``_VisibilityModel.from_observed``
        uses, for the same reason.)
        """
        probe = _ReferenceFourierSample.from_observed(
            container, field_of_view=field_of_view, oversampling=oversampling
        )
        u_pts, v_pts, waves = probe.expanded_coverage
        return cls(
            u_pts,
            v_pts,
            waves * SPECTRAL_UNIT,
            field_of_view=field_of_view,
            oversampling=oversampling,
            label=label,
            **placement,
        )

    # -- chain-internal negotiation ------------------------------------------

    def configure_from(self, downstream: Sequence[Transformation]) -> None:
        """The reference's expansion product, then drop the derived tensors.

        The expanded coverage is what both cached operators are built from, so
        a step reconfigured for a different chain must not keep the tensors it
        built for the previous one.
        """
        super().configure_from(downstream)
        self._forget()

    def _forget(self) -> None:
        self._coverage_cache = None
        self._operator_cache = None

    # -- the operator, built once --------------------------------------------

    def coverage_tensors(self) -> Coverage:
        """The expanded ``(u, v, lambda)`` this step evaluates, as tensors."""
        found = self._coverage_cache
        if found is not None:
            return found
        u_pts, v_pts, waves = self.expanded_coverage
        built = (self._real(u_pts), self._real(v_pts), self._real(waves))
        self._coverage_cache = built
        return built

    def _operator(self, x_mas: torch.Tensor, y_mas: torch.Tensor) -> Operator:
        """``(phase_x, phase_y, solid_angle)`` for this image grid, cached.

        Keyed on the coordinates' exact bytes, for the reason
        :meth:`ampere.backends.jax.instrument._JaxStep._influence_jax` gives:
        the operator is a pure function of the grid, the grid is fixed for a
        run, and rebuilding an ``(n_uv, nx)`` complex matrix per likelihood
        evaluation is the hot-loop waste W2.3 measured. An exact key means a
        changed grid rebuilds and a repeated one does not, so a lookup can
        never return the operator of a nearby grid.
        """
        x_array = to_numpy(x_mas)
        y_array = to_numpy(y_mas)
        key = x_array.tobytes() + b"|" + y_array.tobytes()
        cached = self._operator_cache
        if cached is not None and cached[0] == key:
            return cached[1]
        u_pts, v_pts, _ = self.coverage_tensors()
        # The solid angles are the reference quadrature's, by value: cell
        # widths from the midpoints of whatever grid negotiation produced. It
        # is geometry computed once at composition time, not arithmetic a
        # gradient passes through, which is the same judgement the jax steps
        # make about their influence matrices.
        omega = self._real(cell_solid_angle(x_array, y_array))
        x_rad = self._real(x_array) / MAS_PER_RAD
        y_rad = self._real(y_array) / MAS_PER_RAD
        phase_x = _unit_phase(-2.0 * math.pi * torch.outer(u_pts, x_rad))
        phase_y = _unit_phase(-2.0 * math.pi * torch.outer(v_pts, y_rad))
        built = (phase_x, phase_y, omega)
        self._operator_cache = (key, built)
        return built

    # -- the native surface a realisation composes ---------------------------

    def apply_flux(self, flux: torch.Tensor, grid: Any, values: Any) -> tuple[torch.Tensor, Any]:
        """The DFT as one batched contraction: ``(..., nx, ny) -> (..., n_uv)``.

        *grid* is the image model's ``(x, y)`` pair in mas; what comes back is
        the ``(u, v, lambda)`` triple of the emitted visibilities, because this
        step changes the coordinates as well as the kind.
        """
        x_mas, y_mas = grid
        phase_x, phase_y, omega = self._operator(x_mas, y_mas)
        weighted = (flux * omega).to(phase_x.dtype)
        along_y = torch.einsum("nk,...jk->...nj", phase_y, weighted)
        return torch.einsum("nj,...nj->...n", phase_x, along_y), self.coverage_tensors()

    # -- the contract surface ------------------------------------------------

    def apply(self, samples: Any, values: Any) -> VisibilitySet:
        """The same arithmetic, back on the numpy side of the container boundary.

        Overridden rather than inherited so that the conformance battery — which
        drives ``apply`` — measures *this* backend's transform against the
        closed forms, not the reference backend's against itself.
        """
        grid = (self._real(samples.x.values), self._real(samples.y.values))
        visibility, _ = self.apply_flux(self._real(samples.values), grid, values)
        template = self._resolve_template(samples.unit, *self.expanded_coverage)
        return template.with_values(
            to_numpy(visibility), mask=_dense_mask(samples, self.expanded_coverage[0].size)
        )


class ClosurePhase(TorchInterferometryStep, _ReferenceClosurePhase):
    """Three visibilities to one angle, in torch: ``VisibilitySet -> ClosurePhases``.

    ``arg(V_ij V_jk V_ki)`` for a triangle whose three baselines occupy
    positions ``3t``, ``3t + 1``, ``3t + 2`` — the canonical ordering
    :class:`~ampere.core.ClosurePhases` fixes and
    :meth:`FourierSample.from_observed` produces. ``torch.angle`` of a complex
    product is differentiable and maps under ``vmap``, which is the whole
    reason this step can be in a NUTS chain: the closure phase is where the
    astrometric information survives the atmosphere, so a fit that could not
    differentiate it would be throwing away the signal.

    The ordering check, the mask rule and the output template are the reference
    class's.
    """

    def apply_flux(self, flux: torch.Tensor, grid: Any, values: Any) -> tuple[torch.Tensor, Any]:
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
        return torch.angle(product), coverage

    def apply(self, samples: Any, values: Any) -> ClosurePhases:
        """The same arithmetic, with the reference class's template and mask."""
        total = int(samples.n_samples)
        coverage = tuple(
            self._real(samples.axis(name).values) for name in ("u", "v", "spectral_axis")
        )
        phase, _ = self.apply_flux(self._observed_tensor(samples.values), coverage, values)
        n_out = total // 3
        template = self._resolve_template(samples, n_out)
        mask = propagate_mask(samples, self._influence(n_out) if samples.mask is not None else None)
        return template.with_values(to_numpy(phase), mask=mask)


class Amplitude(TorchInterferometryStep, _ReferenceAmplitude):
    """Take the modulus, in torch: a complex ``VisibilitySet`` to a real one.

    Kind-preserving, coordinate-preserving, parameter-free — and a real
    scientific decision made as a named step rather than as a convention buried
    in a family (``interferometry.md`` §6). ``torch.abs`` of a complex tensor is
    differentiable away from the origin, which is where an amplitude fit lives.
    """

    def apply_flux(self, flux: torch.Tensor, grid: Any, values: Any) -> tuple[torch.Tensor, Any]:
        """``|V|``, with the coordinates untouched."""
        return torch.abs(flux), grid

    def apply(self, samples: Any, values: Any) -> VisibilitySet:
        return samples.with_values(to_numpy(torch.abs(self._observed_tensor(samples.values))))


class _TorchAveragingStep(TorchInterferometryStep):
    """The block reduction both smearing steps share, in torch.

    The *geometry* — which extra ``(u, v)`` samples to average over
    (``expand_uv``), and how many sub-samples the steps after this one added
    (``configure_from``) — is the reference class's, unchanged: it is read by
    :class:`FourierSample` at composition time and is not arithmetic a gradient
    passes through. What is here is the weighted mean itself, as one
    ``einsum`` over the flat layout ``FourierSample.configure_from`` fixes.
    """

    def _weight_tensor(self) -> torch.Tensor:
        """The quadrature weights, as a tensor. Constant; built per call is cheap."""
        return self._real(self._weights)

    def _block(self, total: int) -> tuple[int, int, int]:
        """``(n_out, count, trailing)``, with the reference refusal by name."""
        weights = self._weights
        count = int(weights.size)
        trailing = int(self._trailing)
        block = count * trailing
        if total % block:
            raise TransformationError(
                f"{type(self).__name__} was given {total} visibilities, which is not a multiple "
                f"of the {count} sub-sample(s) it averages over times the {trailing} the steps "
                f"after it added. That happens when the chain's Fourier step was not the one "
                f"configured for this chain: build the step and the instrument together, so "
                f"configure_from runs on the pair ampere will evaluate."
            )
        return total // count, count, trailing

    def apply_flux(self, flux: torch.Tensor, grid: Any, values: Any) -> tuple[torch.Tensor, Any]:
        """The weighted mean of each block, and the un-smeared coordinates.

        Index ``0`` along every sub-axis is the input sample itself
        (``_quadrature``'s centre-first ordering), so the output coordinates
        come from *slicing* the input rather than from recomputing them — which
        is what keeps the bit-identical coordinate rule on the very container a
        likelihood compares against the data.
        """
        n_out, count, trailing = self._block(int(flux.shape[-1]))
        lead = tuple(flux.shape[:-1])
        shaped = flux.reshape(*lead, -1, count, trailing)
        averaged = torch.einsum("a,...nar->...nr", self._weight_tensor().to(shaped.dtype), shaped)
        coverage = tuple(axis.reshape(-1, count, trailing)[:, 0, :].reshape(n_out) for axis in grid)
        return averaged.reshape(*lead, n_out), coverage

    def apply(self, samples: Any, values: Any) -> VisibilitySet:
        """The same average, with the reference class's template and mask rule."""
        coverage = tuple(
            self._real(samples.axis(name).values) for name in ("u", "v", "spectral_axis")
        )
        averaged, reduced = self.apply_flux(self._observed_tensor(samples.values), coverage, values)
        n_out, count, trailing = self._block(int(samples.n_samples))
        mask = None
        if samples.mask is not None:
            mask = propagate_mask(samples, self._influence(n_out, count, trailing))
        return VisibilitySet(
            to_numpy(reduced[0]),
            to_numpy(reduced[1]),
            to_numpy(reduced[2]) * SPECTRAL_UNIT,
            to_numpy(averaged),
            unit=samples.unit,
            mask=mask,
        )


class BandwidthSmearing(_TorchAveragingStep, _ReferenceBandwidthSmearing):
    """Average the visibility across a spectral channel of finite width, in torch.

    The baseline is fixed in metres, so the sampled spatial frequency
    ``B/lambda`` sweeps radially as ``lambda`` crosses the channel and what is
    recorded is the mean visibility over that sweep — a smearing that reduces
    the amplitude of a resolved source at long baselines and looks exactly like
    a larger source if it is not modelled.

    Parameters
    ----------
    resolving_power, nodes, label
        As the reference step.
    dtype, device
        Where this step's constant tensors live.
    """


class TimeSmearing(_TorchAveragingStep, _ReferenceTimeSmearing):
    """Average the visibility over a finite integration, along the uv track, in torch.

    Unlike bandwidth smearing the direction is not radial — it is wherever the
    track goes — so the two are genuinely different steps. The uv rates are
    per-sample buffers given at construction or read from the observed
    container's ``extra_coords`` by :meth:`from_observed`; a
    ``VisibilitySet`` records where a sample was taken, not how fast the
    baseline was moving.

    Parameters
    ----------
    integration, du_dt, dv_dt, nodes, label
        As the reference step.
    dtype, device
        Where this step's constant tensors live.
    """

    @classmethod
    def from_observed(
        cls,
        container: Any,
        *,
        integration: Any,
        nodes: int = _DEFAULT_SMEARING_NODES,
        label: str | None = None,
        **placement: Any,
    ) -> Any:
        """The reference route, with this backend's placement keywords forwarded.

        The rates still come from the observed container's ``extra_coords``, and
        a container without them is still refused by name with the two key names
        it needs: that refusal is the reference classmethod's and is reached by
        calling it. See :meth:`FourierSample.from_observed` for why the override
        exists at all.
        """
        probe = _ReferenceTimeSmearing.from_observed(
            container, integration=integration, nodes=nodes
        )
        return cls(
            integration=integration,
            du_dt=probe.buffers["du_dt"].value,
            dv_dt=probe.buffers["dv_dt"].value,
            nodes=nodes,
            label=label,
            **placement,
        )


# ---------------------------------------------------------------------------
# The image-emitting models
# ---------------------------------------------------------------------------


class _TorchImageModel:
    """Shared plumbing for the image-emitting trio, in torch.

    A mixin over the reference model, for the reason
    :class:`TorchInterferometryStep` is one: the negotiation
    (``compile_for``'s adoption of the union grid, and gap I-4's refusal when
    the model holds its own) is the reference class's, and what is forked is
    the brightness itself — which is what a gradient with respect to a
    separation, a position angle or a flux ratio passes through.

    Each concrete model writes one method, :meth:`_brightness_native`,
    returning a surface brightness in Jy/sr on the ``(x, y)`` tensors it is
    handed, in mas. :meth:`_brightness` routes the contract path through the
    same expression, so ``evaluate`` and ``flux`` cannot disagree.
    """

    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = True
    DEVICE: ClassVar[str] = "cpu"
    BACKEND: ClassVar[str] = BACKEND

    def __init__(
        self,
        *args: Any,
        dtype: torch.dtype = DEFAULT_DTYPE,
        device: torch.device = DEFAULT_DEVICE,
        **kwargs: Any,
    ) -> None:
        super().__init__(*args, **kwargs)
        place(self, dtype, device)
        self.tensors = LoweredParameters()
        for name in ("x", "y"):
            self.tensors.register_buffer(
                name,
                as_tensor(self._data(name), dtype=self.dtype, device=self.device),
                persistent=True,
            )
        self._grids: dict[str, tuple[torch.Tensor, torch.Tensor]] = {}

    # -- negotiation ---------------------------------------------------------

    def compile_for(self, requirements: Mapping[str, ChannelRequirements]) -> Model:
        """Adopt the negotiated image grid, then cache it as tensors."""
        super().compile_for(requirements)
        for channel, template in self.templates.items():
            self._grids[channel] = (
                as_tensor(template.x.values, dtype=self.dtype, device=self.device),
                as_tensor(template.y.values, dtype=self.dtype, device=self.device),
            )
        return self

    def to(self, *, dtype: Any = None, device: Any = None) -> Any:
        """Move this model's grids and re-declare where they live."""
        move(self, dtype=dtype, device=device)
        self._grids = {
            channel: (
                grid[0].to(dtype=self.dtype, device=self.device),
                grid[1].to(dtype=self.dtype, device=self.device),
            )
            for channel, grid in self._grids.items()
        }
        return self

    # -- the native surface a realisation composes ---------------------------

    def native_grid(self, channel: str) -> tuple[torch.Tensor, torch.Tensor]:
        """The ``(x, y)`` tensors this channel evaluates on, mas.

        Two tensors rather than one, because an image lives on two axes; the
        realisation threads whatever this returns through the chain's
        ``apply_flux`` calls without looking inside it, which is what lets one
        value-and-coordinates protocol serve a spectrum and an image.

        **Why ``native_grid`` and not ``grid``**: its partner cannot be called
        ``flux``. Every source model here declares a *parameter* named ``flux``
        (it is the source's total flux density, which is the thing one fits),
        and ``Parameterised._check_free_name`` refuses a parameter whose name
        shadows a class attribute — correctly, since ``self.flux`` would then
        mean two different things. So these models expose the pair under the
        ``native_*`` spelling, which :mod:`ampere.backends.torch.problem`
        accepts beside the original one, and the two names are moved together
        so that a reader never has to wonder which half is which.
        """
        found = self._grids.get(channel)
        if found is not None:
            return found
        return (self.tensors.get_buffer("x"), self.tensors.get_buffer("y"))

    def native_flux(self, channel: str, values: Mapping[str, Any] | None = None) -> torch.Tensor:
        """This model's surface brightness on *channel*, Jy/sr, as a tensor."""
        x_mas, y_mas = self.native_grid(channel)
        return self._brightness_native(x_mas, y_mas, self._context_tensors(self.context(values)))

    def _context_tensors(self, context: Mapping[str, Any]) -> dict[str, torch.Tensor]:
        """This model's vocabulary as tensors, keeping any graph it arrived with."""
        return {
            name: as_tensor(value, dtype=self.dtype, device=self.device)
            for name, value in context.items()
        }

    # -- the contract surface ------------------------------------------------

    def _brightness(
        self, x_mas: np.ndarray, y_mas: np.ndarray, context: Mapping[str, Any]
    ) -> np.ndarray:
        """The torch brightness, back on the numpy side of the container boundary.

        The reference class's ``evaluate`` calls this, so the contract path —
        and therefore every conformance row that compares a transform against a
        closed form — runs torch's arithmetic rather than numpy's.
        """
        self._validate(context)
        return to_numpy(
            self._brightness_native(
                as_tensor(x_mas, dtype=self.dtype, device=self.device),
                as_tensor(y_mas, dtype=self.dtype, device=self.device),
                self._context_tensors(context),
            )
        )

    def _validate(self, context: Mapping[str, Any]) -> None:
        """Refuse a geometrically impossible shape, on the **contract** path only.

        The reference models raise for a non-positive diameter or width. A
        refusal like that cannot live in :meth:`_brightness_native`: under
        ``vmap`` the value is a traced batch and ``float()`` of one is an
        error, and ``inference.md`` §10a's rule is that nothing on a native
        path branches on a value. So the check stays where the values are
        plain numbers — here — and the native path relies on the prior's own
        support instead, which is where a positive quantity is declared
        positive.
        """

    def _brightness_native(
        self, x_mas: torch.Tensor, y_mas: torch.Tensor, context: Mapping[str, torch.Tensor]
    ) -> torch.Tensor:
        raise NotImplementedError


def _gaussian_brightness(
    x_mas: torch.Tensor,
    y_mas: torch.Tensor,
    centre_x: Any,
    centre_y: Any,
    sigma_mas: torch.Tensor,
    flux: torch.Tensor,
) -> torch.Tensor:
    """A circular Gaussian of total *flux* (Jy) as a surface brightness in Jy/sr.

    The reference expression, in torch: the outer sum is formed by indexing the
    two coordinate axes, which are never batched (they are the negotiated
    grid), while every parameter may be — so the same expression serves one θ
    and a ``vmap``\\ ped stack of them.
    """
    sigma_rad = sigma_mas / MAS_PER_RAD
    squared = ((x_mas[:, None] - centre_x) ** 2 + (y_mas[None, :] - centre_y) ** 2) / sigma_mas**2
    return flux / (2.0 * math.pi * sigma_rad**2) * torch.exp(-0.5 * squared)


class UniformDisc(_TorchImageModel, _ReferenceUniformDisc):
    """A uniformly bright circular disc, in torch — **not** differentiable.

    ``flux`` is the total flux density in Jy and ``diameter`` the angular
    diameter in mas; the brightness inside the disc is ``flux / (pi
    (theta/2)**2)`` in Jy/sr and zero outside. The arithmetic is torch's and
    the conformance rows hold it to the same closed form and the same
    convergence measurement as the reference backend's.

    ``DIFFERENTIABLE = False``, declared rather than inherited, and the reason
    is mathematical rather than a missing implementation. A hard-edged disc
    rendered onto a grid is a **piecewise-constant** function of its diameter:
    ``torch.where(distance <= radius, ...)`` moves the edge only when a pixel
    centre crosses it, so autograd returns the derivative of the ``1/r**2``
    amplitude alone and silently omits the boundary term, which is the larger
    half. That is not a small error in a gradient — it points the wrong way
    over most of the parameter range — and a sampler handed it would converge
    confidently to the wrong diameter. Fit a disc on a gradient-free engine, or
    fit :class:`GaussianSource`, whose gradient is exact.

    ``BATCHABLE`` is ``False`` for the same reason and not a second one: the
    flag means "a stack of θ is evaluated in one call over the *realised*
    density", and a model that declares no derivative has no realised density
    to vmap.

    Parameters
    ----------
    x, y, diameter, flux, channels, adopt_grid
        As the reference model.
    dtype, device
        Where this model's grids live.
    """

    DIFFERENTIABLE: ClassVar[bool] = False
    BATCHABLE: ClassVar[bool] = False

    def _validate(self, context: Mapping[str, Any]) -> None:
        radius = 0.5 * float(context["diameter"])
        if radius <= 0.0:
            raise ValueError(f"a uniform disc needs a positive diameter, got {radius * 2.0!r} mas.")

    def _brightness_native(
        self, x_mas: torch.Tensor, y_mas: torch.Tensor, context: Mapping[str, torch.Tensor]
    ) -> torch.Tensor:
        radius = 0.5 * context["diameter"]
        distance = torch.hypot(x_mas[:, None], y_mas[None, :])
        area = math.pi * (radius / MAS_PER_RAD) ** 2
        return torch.where(
            distance <= radius,
            context["flux"] / area,
            torch.zeros((), dtype=self.dtype, device=self.device),
        )


class GaussianSource(_TorchImageModel, _ReferenceGaussianSource):
    """A circular Gaussian source, in torch. Band-limited, and differentiable.

    ``flux`` is the total flux density in Jy and ``fwhm`` the full width at
    half maximum in mas. Band-limited, so the direct DFT of its image agrees
    with its closed-form visibility to the solver tolerance on any grid that
    satisfies the Nyquist requirement and covers the source — which makes it
    the cheapest honest check that this backend's transform is the reference's.

    Parameters
    ----------
    x, y, fwhm, flux, channels, adopt_grid
        As the reference model.
    dtype, device
        Where this model's grids live.
    """

    def _validate(self, context: Mapping[str, Any]) -> None:
        if float(context["fwhm"]) / FWHM_PER_SIGMA <= 0.0:
            raise ValueError(f"a Gaussian source needs a positive fwhm, got {context['fwhm']!r}.")

    def _brightness_native(
        self, x_mas: torch.Tensor, y_mas: torch.Tensor, context: Mapping[str, torch.Tensor]
    ) -> torch.Tensor:
        return _gaussian_brightness(
            x_mas, y_mas, 0.0, 0.0, context["fwhm"] / FWHM_PER_SIGMA, context["flux"]
        )


class Binary(_TorchImageModel, _ReferenceBinary):
    """Two Gaussian components, in torch — the source the phase is proved on.

    The primary sits at the phase centre and the secondary at ``separation``
    and ``position_angle`` (radians, east of north); ``flux_ratio`` is the
    secondary's flux over the primary's and ``flux`` is the pair's total. Both
    components share one ``component_fwhm``, a buffer rather than a parameter:
    it is there to make the source band-limited, and one would not put a prior
    on a numerical device.

    This is the model the NUTS rows fit: its closure phases are non-zero,
    analytic and sensitive to exactly the sign conventions that would otherwise
    go unnoticed, and every one of its four parameters enters the brightness
    through ``sin``, ``cos`` and ``exp``, so the gradient is exact rather than
    a discretisation of one.

    Parameters
    ----------
    x, y, separation, position_angle, flux_ratio, flux, component_fwhm, channels, adopt_grid
        As the reference model.
    dtype, device
        Where this model's grids live.
    """

    def _brightness_native(
        self, x_mas: torch.Tensor, y_mas: torch.Tensor, context: Mapping[str, torch.Tensor]
    ) -> torch.Tensor:
        separation = context["separation"]
        angle = context["position_angle"]
        offset_x = separation * torch.sin(angle)
        offset_y = separation * torch.cos(angle)
        ratio = context["flux_ratio"]
        sigma = context["component_fwhm"] / FWHM_PER_SIGMA
        primary = context["flux"] / (1.0 + ratio)
        return _gaussian_brightness(x_mas, y_mas, 0.0, 0.0, sigma, primary) + _gaussian_brightness(
            x_mas, y_mas, offset_x, offset_y, sigma, primary * ratio
        )


# ---------------------------------------------------------------------------
# The analytic route: the same three sources, emitting visibilities directly
# ---------------------------------------------------------------------------


class _TorchVisibilityModel:
    """Shared plumbing for the trio that skips the Fourier step, in torch.

    ``interferometry.md`` §1: a model that computes visibilities analytically
    emits a :class:`~ampere.core.VisibilitySet` channel directly, and its
    instrument is the no-steps one. Both routes are supported and neither is
    privileged — and here each is the other's oracle, which is why these are
    worth twinning even though nothing fits with them: a closed form in torch
    is what says this backend's direct transform is right for a reason other
    than "numpy agrees with numpy".

    The coverage buffers and ``from_observed`` are the reference class's; each
    concrete model writes :meth:`_visibility_native`.
    """

    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = True
    DEVICE: ClassVar[str] = "cpu"
    BACKEND: ClassVar[str] = BACKEND

    def __init__(
        self,
        *args: Any,
        dtype: torch.dtype = DEFAULT_DTYPE,
        device: torch.device = DEFAULT_DEVICE,
        **kwargs: Any,
    ) -> None:
        super().__init__(*args, **kwargs)
        place(self, dtype, device)
        self.tensors = LoweredParameters()
        for name in ("u_pts", "v_pts", "wavelength"):
            self.tensors.register_buffer(
                name,
                as_tensor(
                    np.asarray(self.buffers[name].value, dtype=float),
                    dtype=self.dtype,
                    device=self.device,
                ),
                persistent=True,
            )

    def to(self, *, dtype: Any = None, device: Any = None) -> Any:
        """Move this model's coverage buffers and re-declare where they live."""
        return move(self, dtype=dtype, device=device)

    # -- the native surface a realisation composes ---------------------------

    def native_grid(self, channel: str) -> Coverage:
        """The ``(u, v, lambda)`` tensors this model already emits on.

        ``native_*`` rather than ``grid``/``flux`` for the reason
        :meth:`_TorchImageModel.native_grid` gives: ``flux`` is one of this
        model's own parameters.
        """
        return (
            self.tensors.get_buffer("u_pts"),
            self.tensors.get_buffer("v_pts"),
            self.tensors.get_buffer("wavelength"),
        )

    def native_flux(self, channel: str, values: Mapping[str, Any] | None = None) -> torch.Tensor:
        """This model's complex visibilities on *channel*, Jy, as a tensor."""
        context = self._context_tensors(self.context(values))
        u_pts, v_pts, _ = self.native_grid(channel)
        return self._visibility_native(torch.hypot(u_pts, v_pts), context)

    def _context_tensors(self, context: Mapping[str, Any]) -> dict[str, torch.Tensor]:
        return {
            name: as_tensor(value, dtype=self.dtype, device=self.device)
            for name, value in context.items()
        }

    # -- the contract surface ------------------------------------------------

    def _visibility(self, rho: np.ndarray, context: Mapping[str, Any]) -> np.ndarray:
        """The torch closed form, back on the numpy side of the boundary."""
        return to_numpy(
            self._visibility_native(
                as_tensor(rho, dtype=self.dtype, device=self.device),
                self._context_tensors(context),
            )
        )

    def _visibility_native(
        self, rho: torch.Tensor, context: Mapping[str, torch.Tensor]
    ) -> torch.Tensor:
        raise NotImplementedError

    def _complex(self, value: torch.Tensor) -> torch.Tensor:
        return value.to(complex_dtype(self.dtype))


class UniformDiscVisibilities(_TorchVisibilityModel, _ReferenceUniformDiscVisibilities):
    """``V = flux * 2 J1(pi theta rho) / (pi theta rho)`` in torch — **not** differentiable.

    The disc's closed form, real-valued and negative beyond the first null,
    which is the sign flip an amplitude-only fit throws away. The arithmetic is
    ``torch.special.bessel_j1``, which — unlike ``jax.scipy.special.bessel_jn``
    — **has no backward**: it produces a tensor with no ``grad_fn``, so a
    density built on it would report a gradient of exactly nothing rather than
    raise. Declaring ``DIFFERENTIABLE = False`` is what turns that into a
    refusal by name at composition; it is a torch library gap, not a decision,
    and the remedy in the message is the jax twin or a gradient-free engine.

    Parameters
    ----------
    u_pts, v_pts, wavelength, diameter, flux, channels
        As the reference model.
    dtype, device
        Where this model's coverage buffers live.
    """

    DIFFERENTIABLE: ClassVar[bool] = False
    BATCHABLE: ClassVar[bool] = False

    def _visibility_native(
        self, rho: torch.Tensor, context: Mapping[str, torch.Tensor]
    ) -> torch.Tensor:
        argument = math.pi * context["diameter"] / MAS_PER_RAD * rho
        zero = argument == 0.0
        # 2 J1(z)/z -> 1 as z -> 0, and the ratio is 0/0 there; the guard is
        # written in rather than clipped, and the safe argument keeps the
        # Bessel call itself away from the singular point.
        safe = torch.where(zero, torch.ones_like(argument), argument)
        envelope = torch.where(
            zero, torch.ones_like(argument), 2.0 * torch.special.bessel_j1(safe) / safe
        )
        return self._complex(context["flux"] * envelope)


class GaussianSourceVisibilities(_TorchVisibilityModel, _ReferenceGaussianSourceVisibilities):
    """``V = flux * exp(-2 pi**2 sigma**2 rho**2)`` — the Gaussian's closed form, in torch."""

    def _visibility_native(
        self, rho: torch.Tensor, context: Mapping[str, torch.Tensor]
    ) -> torch.Tensor:
        sigma = context["fwhm"] / FWHM_PER_SIGMA / MAS_PER_RAD
        return self._complex(context["flux"] * torch.exp(-2.0 * math.pi**2 * sigma**2 * rho**2))


class BinaryVisibilities(_TorchVisibilityModel, _ReferenceBinaryVisibilities):
    """``V = (flux/(1+f)) G(rho) (1 + f exp(-2 pi i (u dx + v dy)))`` in torch.

    The primary at the phase centre and the secondary at ``(dx, dy)`` radians,
    the same geometry :class:`Binary` renders; ``G`` is the components' common
    Gaussian envelope. The exponent carries the ``exp(-2 pi i)`` sign
    convention, which is what makes the closure phases non-zero — and what the
    conformance rows check the sign of against a closed form rather than
    against this module.
    """

    def _visibility_native(
        self, rho: torch.Tensor, context: Mapping[str, torch.Tensor]
    ) -> torch.Tensor:
        separation = context["separation"]
        angle = context["position_angle"]
        offset_x = separation * torch.sin(angle)
        offset_y = separation * torch.cos(angle)
        sigma = context["component_fwhm"] / FWHM_PER_SIGMA / MAS_PER_RAD
        envelope = torch.exp(-2.0 * math.pi**2 * sigma**2 * rho**2)
        angle_uv = (
            -2.0
            * math.pi
            * (context["u_pts"] * offset_x + context["v_pts"] * offset_y)
            / MAS_PER_RAD
        )
        ratio = context["flux_ratio"]
        return self._complex(context["flux"] / (1.0 + ratio) * envelope) * (
            1.0 + self._complex(ratio) * _unit_phase(angle_uv)
        )
