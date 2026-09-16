"""The native spectral models on the torch path: blackbody, modified blackbody, power law.

``DEVELOPMENT_PLAN.md`` §5's Phase 2 list again, computed in ``torch`` this
time. Same physics, same parameter names, same units, same channels as
:mod:`ampere.backends.reference.models` — deliberately, because the two are
held to each other by ``tests/conformance/test_cross_backend.py``, and a
renamed parameter or a different reference wavelength would be a failure rather
than a detail.

What is different is the arithmetic and the two things that follow from it.

**The grid is a torch buffer** (``lowering.md`` §7). ``register_buffer`` gives
exactly the properties ``architecture.md`` §5 asks of a buffer: traced through
computations, moved by ``.to(device)``/``.to(dtype)`` alongside parameters,
never differentiated. The same array is *also* registered with
``ampere.core``'s :meth:`~ampere.core.Parameterised.register_buffer`, because
that is the declared surface the contract, the provenance record and the
neutral model identity all read. The two are the same numbers under two
ownerships — a declaration and a runtime — and neither can be dropped: core's
``BufferSet`` is what ``describe()`` and the model fingerprint see, and torch's
is what a device move and a checkpoint see.

**Two evaluation entry points, and why.** :meth:`~TorchSpectralModel.evaluate`
is ``ampere.core``'s contract: it returns a
:class:`~ampere.core.ModelResult` of :class:`~ampere.core.Spectrum` containers,
which hold **numpy** arrays — the containers are backend-neutral and validate
their axes with ``numpy``, so a tensor carrying an autograd graph cannot go in
one. :meth:`~TorchSpectralModel.evaluate_tensor` is the same computation
without that boundary: it returns the channel tensors themselves, so a caller
that wants a gradient through the forward model has one.

That split is the honest statement of where this backend's differentiability
currently reaches, and it is worth being exact about it because
``DIFFERENTIABLE = True`` is a load-bearing claim. These models *are*
differentiable — ``evaluate_tensor`` is a differentiable function of tensor
parameters, and ``tests/backends/test_torch_models.py`` proves it by taking the
gradient — but an end-to-end ``d log_prob / d theta`` additionally needs the
instrument chain and the likelihood to run on tensors, and ``ampere.core``'s
containers are the place that path stops today. See this package's
``__init__`` docstring for the full statement; it is W2.4's principal finding.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, ClassVar

import astropy.constants as const
import astropy.units as u
import numpy as np
import torch

from ampere.core import (
    DTYPE,
    ChannelRequirements,
    Model,
    ModelResult,
    Spectrum,
)
from ._config import BACKEND, DEFAULT_DEVICE, DEFAULT_DTYPE, as_tensor, move, place, to_numpy
from ._declare import as_parameter
from .parameters import LoweredParameters

__all__ = [
    "COORDINATE_UNIT",
    "FLUX_UNIT",
    "BlackBody",
    "ModifiedBlackBody",
    "PowerLaw",
    "TorchSpectralModel",
    "planck_jy",
]

#: The spectral coordinate unit every model here works in.
COORDINATE_UNIT = u.micron
#: The flux unit every model here emits.
FLUX_UNIT = u.Jy

# Converted once, at import, into the bare-float constants the hot loop uses
# (DEVELOPMENT_PLAN.md §7: units are converted at composition time, never per
# evaluation). Identical to the reference backend's, and deliberately taken
# from the same astropy constants rather than transcribed, so the two backends
# cannot drift in the tenth digit for a reason that has nothing to do with
# arithmetic.
_TWO_H_OVER_C2 = float((2.0 * const.h / const.c**2).to_value(u.J * u.s**3 / u.m**2))
_H_OVER_K = float((const.h / const.k_B).to_value(u.s * u.K))
_C_MICRON_HZ = float(const.c.to_value(u.micron * u.Hz))
#: 1 Jy in SI (W m^-2 Hz^-1).
_JY = 1.0e-26


def planck_jy(wavelength: torch.Tensor, temperature: torch.Tensor) -> torch.Tensor:
    """``B_nu(T)`` in Jy/sr for *wavelength* in micron, differentiable in both.

    Written in frequency because the flux unit is Jy — a per-frequency density
    — so evaluating ``B_lambda`` and converting back would only add a
    ``lambda**2`` round trip and its rounding. ``expm1`` rather than
    ``exp(x) - 1`` is what keeps the Rayleigh-Jeans tail honest: at long
    wavelengths ``h nu / kT`` is tiny and the subtraction loses every
    significant digit to cancellation.

    The Wien-side guard differs from the reference backend's in mechanism and
    not in result. numpy can compute ``inf`` and replace it afterwards;
    ``torch.where`` evaluates **both** branches and would propagate a ``nan``
    gradient from the discarded one, so the exponent is clamped *before* the
    ``expm1`` instead. ``expm1(700)`` is about ``1e304``, comfortably finite in
    float64, and the radiance it gives is ``~1e-300`` — indistinguishable from
    the zero the reference backend substitutes, and reached without ever
    forming an ``inf``.
    """
    if not bool(torch.all(wavelength > 0.0)):
        raise ValueError("the Planck function needs strictly positive wavelengths (micron).")
    if not bool(torch.all(torch.isfinite(temperature)) and torch.all(temperature > 0.0)):
        raise ValueError(f"blackbody temperature must be finite and positive, got {temperature!r}.")
    frequency = _C_MICRON_HZ / wavelength
    exponent = torch.clamp(_H_OVER_K * frequency / temperature, max=700.0)
    radiance = _TWO_H_OVER_C2 * frequency**3 / torch.expm1(exponent)
    return radiance / _JY


class TorchSpectralModel(Model):
    """Shared plumbing: the grid buffer, the channels, and the template cache.

    The reference backend's ``_SpectralModel`` with a torch interior. The
    template cache is the same idea and for the same reason
    (``transformations.md`` §14): build one container per channel from the
    negotiated grid at :meth:`compile_for` time, and refill it with
    ``with_values`` on every evaluation, so the axes are validated once rather
    than per draw and successive results share them by identity.
    """

    #: Name of the spectral-axis requirement these models answer.
    AXIS: ClassVar[str] = "spectral_axis"

    #: The four capability flags (``inference.md`` §18, W2.12's fourth), stated
    #: rather than inherited. ``DIFFERENTIABLE`` is the claim
    #: :meth:`evaluate_tensor` makes good on.
    #: ``BATCHABLE`` is ``True`` since W2.4 slice 2, and it means one specific
    #: thing: a stack of parameter vectors is evaluated in one call through
    #: ``torch.func.vmap``, over the *realised* density
    #: (:meth:`ampere.backends.torch.LoweredProblem.log_prob_unconstrained_batched`).
    #: It is a claim about this class only, aggregated conjunctively with every
    #: other part's, so one non-batchable piece — ``QuasisepGP``, whose solve
    #: happens inside a compiled extension — makes the whole problem report
    #: ``batchable=False``, and the batched call then refuses by name.
    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = True
    #: The class-level default only. Since W2.4 slice 3 ``__init__`` sets an
    #: **instance** attribute of the same name from its ``device=`` keyword, so
    #: a model built on ``"cuda"`` declares ``"cuda"`` and
    #: :func:`~ampere.core.declared_capabilities` — which reads ``part.DEVICE``
    #: by attribute access — composes the whole problem on that device or
    #: refuses. Never auto-detected (``architecture.md`` §5).
    DEVICE: ClassVar[str] = "cpu"
    BACKEND: ClassVar[str] = BACKEND

    def __init__(
        self,
        wavelength: Any,
        *,
        channels: str | Sequence[str] = "default",
        dtype: torch.dtype = DEFAULT_DTYPE,
        device: torch.device = DEFAULT_DEVICE,
    ) -> None:
        grid = _to_micron(wavelength)
        if grid.ndim != 1 or grid.size == 0:
            raise ValueError(
                f"a spectral model needs a one-dimensional, non-empty wavelength grid, got shape "
                f"{grid.shape}."
            )
        names = (channels,) if isinstance(channels, str) else tuple(str(c) for c in channels)
        if not names or len(set(names)) != len(names):
            raise ValueError(
                f"a spectral model needs distinct, non-empty channel names, got {names!r}."
            )
        self.channels = names
        # dtype, device and the DEVICE capability flag, set together: see
        # ``_config.place``. This is the line that makes the device a property
        # of the instance rather than of the class.
        place(self, dtype, device)
        #: The torch-side home for this model's constant arrays (``lowering.md``
        #: §7). Buffers rather than parameters: no gradient, but they move with
        #: ``.to()`` and they are in the ``state_dict``.
        self.tensors = LoweredParameters()
        self.register_buffer("wavelength", grid, unit=COORDINATE_UNIT)
        self.tensors.register_buffer(
            "wavelength", as_tensor(grid, dtype=self.dtype, device=self.device), persistent=True
        )
        self.templates: dict[str, Spectrum] = {}
        self._grids: dict[str, torch.Tensor] = {}

    # -- grids ----------------------------------------------------------------

    def compile_for(self, requirements: Mapping[str, ChannelRequirements]) -> Model:
        """Adopt the negotiated grid for each of this model's channels, once."""
        for channel in self.channels:
            asked = requirements.get(channel)
            if asked is None or self.AXIS not in asked:
                continue
            grid = _to_micron(asked[self.AXIS].coordinates())
            self.templates[channel] = Spectrum(
                grid * COORDINATE_UNIT, np.zeros(grid.size, dtype=DTYPE), unit=FLUX_UNIT
            )
            self._grids[channel] = as_tensor(grid, dtype=self.dtype, device=self.device)
        return self

    def to(self, *, dtype: Any = None, device: Any = None) -> TorchSpectralModel:
        """Move this model's buffers, and re-declare where they live.

        ``architecture.md`` §5's "buffers move with parameters under
        ``.to(...)``". Every constant this model owns is registered on
        :attr:`tensors`, an ``nn.Module``, so torch's own recursion does the
        move; the negotiated per-channel grids :meth:`compile_for` cached are
        moved beside them, because they are the same kind of thing under a
        different owner. In place, returning ``self``.
        """
        move(self, dtype=dtype, device=device)
        self._grids = {
            channel: grid.to(dtype=self.dtype, device=self.device)
            for channel, grid in self._grids.items()
        }
        return self

    def grid_tensor(self, channel: str) -> torch.Tensor:
        """The tensor grid this channel evaluates on: negotiated, or the model's own."""
        found = self._grids.get(channel)
        if found is not None:
            return found
        buffer = self.tensors.get_buffer("wavelength")
        return buffer

    # -- the native surface a realisation composes (W2.13) --------------------

    def native_grid(self, channel: str) -> torch.Tensor:
        """:meth:`grid_tensor` under the name the realisation looks for.

        :mod:`ampere.backends.torch.problem` walks a chain through
        ``model.native_grid`` / ``model.native_flux`` / ``step.apply_flux``,
        and so does
        :mod:`ampere.backends.jax.problem`. The two spellings are the same
        surface deliberately: the two lowered problems are the same walk in
        two libraries, and a reviewer reading them side by side should see
        that rather than have to establish it.
        """
        return self.grid_tensor(channel)

    def native_flux(self, channel: str, values: Mapping[str, Any] | None = None) -> torch.Tensor:
        """This model's flux on *channel*, in Jy, as a differentiable tensor.

        The surface a realisation composes, and the one the gradient passes
        through. :meth:`evaluate_tensor` is the same computation for every
        declared channel at once, which is what the contract path wants; this
        is one channel, which is what a lowered chain wants.
        """
        return self._flux(self.grid_tensor(channel), self._context_tensors(values))

    def _emit(self, channel: str, grid: torch.Tensor, flux: torch.Tensor) -> Spectrum:
        values = to_numpy(flux).astype(DTYPE, copy=False)
        template = self.templates.get(channel)
        if template is None:
            return Spectrum(to_numpy(grid) * COORDINATE_UNIT, values, unit=FLUX_UNIT)
        return template.with_values(values)

    # -- evaluation -----------------------------------------------------------

    def _flux(self, grid: torch.Tensor, context: Mapping[str, torch.Tensor]) -> torch.Tensor:
        raise NotImplementedError

    def _context_tensors(self, values: Mapping[str, Any] | None) -> dict[str, torch.Tensor]:
        """This model's full vocabulary, as tensors on this model's dtype/device."""
        context = self.context(values)
        tensors = {
            name: as_tensor(value, dtype=self.dtype, device=self.device)
            for name, value in context.items()
        }
        tensors["wavelength"] = self.tensors.get_buffer("wavelength")
        return tensors

    def evaluate_tensor(self, **values: Any) -> dict[str, torch.Tensor]:
        """The forward model as tensors — the differentiable entry point.

        Accepts tensors (keeping their graph) as readily as floats, so
        ``torch.autograd.grad(model.evaluate_tensor(**tensors)["default"].sum(),
        tensor)`` is the model's Jacobian-vector product. Returns one tensor per
        declared channel, on this model's own dtype and device.
        """
        context = self._context_tensors(values)
        return {
            channel: self._flux(self.grid_tensor(channel), context) for channel in self.channels
        }

    def evaluate(self, **values: Any) -> ModelResult:
        """``ampere.core``'s contract: one :class:`~ampere.core.Spectrum` per channel.

        The tensors come back to numpy here, at the container boundary, because
        that is where they must: a :class:`~ampere.core.Spectrum` validates its
        axes in numpy and cannot hold a tensor carrying an autograd graph.
        :meth:`evaluate_tensor` is the same computation without the boundary.
        """
        emitted = {
            channel: self._emit(channel, self.grid_tensor(channel), flux)
            for channel, flux in self.evaluate_tensor(**values).items()
        }
        return ModelResult(emitted)


class BlackBody(TorchSpectralModel):
    """``F_nu = scale * B_nu(T)`` — an isothermal blackbody in Jy.

    ``scale`` absorbs the solid angle (and any distance dilution): the
    dimensionless factor turning the Planck radiance into an observed flux
    density.

    Parameters
    ----------
    wavelength
        The model's own grid, micron (bare, or a :class:`~astropy.units.Quantity`).
    temperature
        Kelvin. A prior to fit it, a number to hold it fixed.
    scale
        Dimensionless multiplier.
    channels
        Name (or names) of the channel the emitted spectrum appears under.
    dtype, device
        Threaded explicitly (``lowering.md`` §10.1); float64 on the CPU by
        default, and never taken from ``torch.get_default_dtype()``.
    """

    def __init__(
        self,
        wavelength: Any,
        *,
        temperature: Any = 1000.0,
        scale: Any = 1.0,
        channels: str | Sequence[str] = "default",
        dtype: torch.dtype = DEFAULT_DTYPE,
        device: torch.device = DEFAULT_DEVICE,
    ) -> None:
        super().__init__(wavelength, channels=channels, dtype=dtype, device=device)
        self.register_parameter(as_parameter("temperature", temperature, unit=u.K))
        self.register_parameter(as_parameter("scale", scale))

    def _flux(self, grid: torch.Tensor, context: Mapping[str, torch.Tensor]) -> torch.Tensor:
        return context["scale"] * planck_jy(grid, context["temperature"])


class ModifiedBlackBody(TorchSpectralModel):
    """``F_nu = scale * (lambda_0 / lambda)**beta * B_nu(T)`` — optically thin dust.

    ``reference_wavelength`` is where the emissivity is unity, so ``scale``
    stays interpretable as the flux the source would have if it radiated as a
    pure blackbody there; without it ``scale`` and ``beta`` are degenerate in a
    way that makes the posterior hard to read. It is a **buffer**, not a
    parameter — one would not put a prior on the definition of one's own
    normalisation.
    """

    def __init__(
        self,
        wavelength: Any,
        *,
        temperature: Any = 100.0,
        beta: Any = 1.5,
        scale: Any = 1.0,
        reference_wavelength: float = 250.0,
        channels: str | Sequence[str] = "default",
        dtype: torch.dtype = DEFAULT_DTYPE,
        device: torch.device = DEFAULT_DEVICE,
    ) -> None:
        super().__init__(wavelength, channels=channels, dtype=dtype, device=device)
        reference = float(_to_micron(reference_wavelength))
        if not np.isfinite(reference) or reference <= 0.0:
            raise ValueError(
                f"reference_wavelength must be finite and positive (micron), got {reference!r}."
            )
        self.register_buffer("reference_wavelength", reference, unit=COORDINATE_UNIT)
        self.tensors.register_buffer(
            "reference_wavelength", as_tensor(reference, dtype=dtype, device=device)
        )
        self.register_parameter(as_parameter("temperature", temperature, unit=u.K))
        self.register_parameter(as_parameter("beta", beta))
        self.register_parameter(as_parameter("scale", scale))

    def _flux(self, grid: torch.Tensor, context: Mapping[str, torch.Tensor]) -> torch.Tensor:
        emissivity = (context["reference_wavelength"] / grid) ** context["beta"]
        return context["scale"] * emissivity * planck_jy(grid, context["temperature"])


class PowerLaw(TorchSpectralModel):
    """``F_nu = norm * (lambda / lambda_ref)**index`` — a spectral power law.

    Written in wavelength, which is the axis the v1 slice's containers use, so
    a spectral index quoted in frequency changes sign. ``reference_wavelength``
    is a buffer for the same reason as in :class:`ModifiedBlackBody`: it pins
    what ``norm`` means.
    """

    def __init__(
        self,
        wavelength: Any,
        *,
        norm: Any = 1.0,
        index: Any = -1.0,
        reference_wavelength: float = 1.0,
        channels: str | Sequence[str] = "default",
        dtype: torch.dtype = DEFAULT_DTYPE,
        device: torch.device = DEFAULT_DEVICE,
    ) -> None:
        super().__init__(wavelength, channels=channels, dtype=dtype, device=device)
        reference = float(_to_micron(reference_wavelength))
        if not np.isfinite(reference) or reference <= 0.0:
            raise ValueError(
                f"reference_wavelength must be finite and positive (micron), got {reference!r}."
            )
        self.register_buffer("reference_wavelength", reference, unit=COORDINATE_UNIT)
        self.tensors.register_buffer(
            "reference_wavelength", as_tensor(reference, dtype=dtype, device=device)
        )
        self.register_parameter(as_parameter("norm", norm))
        self.register_parameter(as_parameter("index", index))

    def _flux(self, grid: torch.Tensor, context: Mapping[str, torch.Tensor]) -> torch.Tensor:
        ratio = grid / context["reference_wavelength"]
        return context["norm"] * ratio ** context["index"]


def _to_micron(coordinates: Any) -> np.ndarray:
    """Coordinates as bare micron, whether or not they arrive as a Quantity."""
    if isinstance(coordinates, u.Quantity):
        return np.asarray(coordinates.to_value(COORDINATE_UNIT), dtype=DTYPE)
    if isinstance(coordinates, torch.Tensor):
        return to_numpy(coordinates).astype(DTYPE, copy=False)
    return np.asarray(coordinates, dtype=DTYPE)
