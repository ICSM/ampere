"""The opt-in astropy→native translation hook: ``contracts/astropy_compat.md`` §5.

The other half of the 2026-09-01 ruling in ``DEVELOPMENT_PLAN.md`` §2. The
default route, :func:`ampere.core.from_astropy`, evaluates the user's **actual**
``astropy.modeling`` model on the numpy path, and therefore has no gradient.
This module is how a caller asks for a *native* equivalent instead — explicitly,
by backend, in one call — and it raises for anything it has no curated row for
rather than falling back to the black box.

That refusal is the whole point. Silently returning the adapter would hand
somebody who asked for a differentiable model one that is not, and they would
discover it when ``NUTSEngine`` refused, several composition steps away from the
call that caused it. Silently substituting a lookalike for a model that *is* in
the table, without being asked, is the mirror-image failure and is what the
ruling forbids outright: a curated ``BlackBody`` is not guaranteed numerically
identical to astropy's, and a model whose implementation changed without its
author knowing is exactly the silent downgrade this architecture exists to
prevent.

**W4.6 landed the hook and an empty table; W4.7 fills it**
(``BlackBody``, ``PowerLaw1D``, ``BrokenPowerLaw1D``, ``Polynomial1D``,
``Gaussian1D``, ``Const1D`` and their compound sums/products/differences/
ratios). The physics is shared with the torch backend's twin of this module —
:mod:`ampere.core.astropy_translations` writes each curated leaf's formula
once, against a two-method namespace plus this backend's own
:func:`~ampere.backends.jax.models.planck_jy`, and both backends' tables are
built from the same dict (see that module's docstring, "one table serves both
backends"). What is genuinely backend-specific, and lives here, is turning the
translated declaration into an ``ampere.core.Model``: jax arrays, this
backend's parameter/buffer registration, and the pure ``flux`` surface
(:mod:`ampere.backends.jax.models`'s convention) alongside the contract
``evaluate``.

A compound model carrying an astropy ``tied=`` parameter, or combined with an
operator outside ``+ - * /`` (astropy's ``|``/``&``), is refused by name here
too — see :mod:`ampere.core.astropy_translations`'s docstring for why neither
is given a native meaning (a tied parameter's derived value cannot even be
computed while jax is tracing).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, ClassVar

import jax
import jax.numpy as jnp
import numpy as np

from ampere.core import (
    DEFAULT_CHANNEL,
    ChannelRequirements,
    FunctionSamples,
    Model,
    ModelResult,
    astropy_components,
    from_astropy as adapt_astropy,
    translation_refusal,
)
from ampere.core.astropy_translations import (
    LEAF_BUILDERS,
    LeafPlan,
    leaf_grid_values,
    plan_translation,
    tied_parameter_refusal,
)
from ampere.core.exceptions import CompositionError, TransformationError

from ._config import BACKEND, require_x64
from ._device import DEVICE, device_flag, place_on, resolve_device
from .models import planck_jy

__all__ = ["TRANSLATIONS", "NativeAstropyModel", "from_astropy"]

#: The curated table: astropy model class → native constructor, consulted leaf
#: by leaf so that a compound model is translatable exactly when every one of
#: its components is.
#:
#: Built from :data:`ampere.core.astropy_translations.LEAF_BUILDERS`, the
#: backend-neutral physics — this backend contributes nothing to the table
#: itself, only the array type and the buffer/parameter plumbing
#: :class:`NativeAstropyModel` wraps it in. A copy rather than an alias, so
#: nothing outside this module can grow the neutral table through this one's
#: name.
TRANSLATIONS: dict[type, Any] = dict(LEAF_BUILDERS)


class _JaxOps:
    """:class:`~ampere.core.astropy_translations.TranslationOps` in jax."""

    def exp(self, array: Any) -> Any:
        return jnp.exp(array)

    def where(self, condition: Any, if_true: Any, if_false: Any) -> Any:
        return jnp.where(condition, if_true, if_false)

    def blackbody(self, wavelength: Any, temperature: Any) -> Any:
        return planck_jy(wavelength, temperature)


_OPS = _JaxOps()


class NativeAstropyModel(Model):
    """A curated astropy model (or compound), computed natively in jax.

    Built by :func:`from_astropy`, the documented entry point. Wraps an
    internal :class:`~ampere.core.astropy_compat.AdaptedAstropyModel` — the
    "probe" — purely for the bookkeeping ``ampere.core`` already gets right
    and that a second implementation must not be allowed to drift from: kind
    and grid resolution, the negotiation protocol, the output-unit factor
    (§4's solid-angle rule included), and the container this model refills on
    every evaluation. The *flux* is the one thing this class computes
    natively, replacing the probe's call into the actual astropy model with
    the curated, differentiable formula tree
    :func:`~ampere.core.astropy_translations.plan_translation` builds.

    Parameters, priors, frozen-ness and channel all come from the probe's own
    :func:`~ampere.core.astropy_compat.translate_astropy_parameters` — the
    same table :func:`ampere.core.from_astropy` uses — so the two routes
    cannot drift on what a bound or a fixed value means. A tied parameter is
    refused (:func:`~ampere.core.astropy_translations.tied_parameter_refusal`):
    see :mod:`ampere.core.astropy_translations`'s docstring.
    """

    DIFFERENTIABLE: ClassVar[bool] = True
    #: Not verified batchable (astropy's own ``where``-branching formula and
    #: the per-leaf grid conversion are not exercised under ``vmap`` here), so
    #: this declares the conservative, honest default rather than an untested
    #: claim — a capability lie is worse than a missing capability.
    BATCHABLE: ClassVar[bool] = False
    DEVICE: ClassVar[str] = DEVICE
    BACKEND: ClassVar[str] = BACKEND

    def __init__(
        self,
        model: Any,
        *,
        kind: type[FunctionSamples] | None = None,
        priors: Mapping[str, Any] | None = None,
        channel: str = DEFAULT_CHANNEL,
        grid: Any = None,
        output_unit: Any = None,
        equivalencies: Sequence[Any] = (),
    ) -> None:
        require_x64("a jax NativeAstropyModel translation")
        self._probe = adapt_astropy(
            model,
            kind=kind,
            priors=priors,
            channel=channel,
            grid=grid,
            output_unit=output_unit,
            equivalencies=equivalencies,
        )
        if self._probe.ties:
            raise tied_parameter_refusal(BACKEND, sorted(self._probe.ties))
        self._leaves, self._combine = plan_translation(self._probe.astropy_model, BACKEND)
        for parameter in self._probe.parameters:
            self.register_parameter(parameter)
        self.channel = self._probe.channel
        self.kind = self._probe.kind
        self._axis_name = self.kind.AXES[0].name
        resolved = resolve_device(DEVICE, "a jax NativeAstropyModel", error=TransformationError)
        object.__setattr__(self, "_resolved_device", resolved)
        object.__setattr__(self, "DEVICE", device_flag(DEVICE, resolved))
        self._axis_grid_array: jax.Array = jnp.zeros(0, dtype=jnp.float64)
        self._leaf_grids: tuple[jax.Array, ...] = ()
        self._factor_array: jax.Array = jnp.asarray(1.0, dtype=jnp.float64)
        if self._probe.template is not None:
            self._adopt_grid()

    def compile_for(self, requirements: Mapping[str, ChannelRequirements]) -> Model:
        """Delegate the negotiation to the probe, then adopt what it resolved."""
        self._probe.compile_for(requirements)
        if self._probe.template is not None:
            self._adopt_grid()
        return self

    def _adopt_grid(self) -> None:
        """Mirror the probe's axis buffer, and build this leaf's own grid arrays, once."""
        for buffer in self._probe.buffers:
            if buffer.name in self.buffers:
                self._buffers = self.buffers.without(buffer.name)
            self.register_buffer(
                buffer.name, buffer.array, unit=buffer.unit, description=buffer.description
            )
        grid = self._probe.grids[0]
        resolved = getattr(self, "_resolved_device", None)
        self._axis_grid_array = place_on(
            jnp.asarray(self.buffers[self._axis_name].array, dtype=jnp.float64), resolved
        )
        self._leaf_grids = tuple(
            place_on(jnp.asarray(leaf_grid_values(plan.formula, grid), dtype=jnp.float64), resolved)
            for plan in self._leaves
        )
        self._factor_array = place_on(jnp.asarray(self._probe.factor, dtype=jnp.float64), resolved)

    # -- the native surface a realisation composes (W2.13) -------------------

    def _check_grid(self) -> None:
        if not self._leaf_grids:
            raise CompositionError(
                f"{type(self).__name__} has no grid to evaluate on. Either give one — "
                f"from_astropy(model, grid=wavelength * u.micron, ...) — or compose it into a "
                f"FittingProblem, whose negotiation calls compile_for() with the coordinates "
                f"the instruments asked for."
            )

    def _check_channel(self, channel: str) -> None:
        if channel != self.channel:
            raise CompositionError(
                f"{type(self).__name__} was asked for channel {channel!r}, but it only emits "
                f"{self.channel!r} — one astropy model is one channel (§6, "
                f"docs/design/contracts/astropy_compat.md)."
            )

    def grid(self, channel: str) -> jax.Array:
        """The jax coordinates *channel* evaluates on: the negotiated axis, as buffered.

        :mod:`ampere.backends.jax.problem` walks a chain through
        ``model.grid`` / ``model.flux`` / ``step.apply_flux``, and this is the
        same surface every other native jax model here offers.
        """
        self._check_channel(channel)
        self._check_grid()
        return self._axis_grid_array

    def flux(self, channel: str, values: Mapping[str, Any] | None = None) -> jax.Array:
        """This model's flux on *channel*, pure and traceable.

        The surface a gradient passes through.
        """
        self._check_channel(channel)
        self._check_grid()
        context = self.context(values)
        leaf_arrays = [
            _evaluate_leaf(plan, grid, context)
            for plan, grid in zip(self._leaves, self._leaf_grids, strict=True)
        ]
        return self._combine(leaf_arrays) * self._factor_array

    def evaluate(self, **values: Any) -> ModelResult:
        """``ampere.core``'s contract: the composed flux, filed under this model's channel.

        The array comes back to numpy here, at the container boundary
        (:mod:`ampere.backends.jax.models`'s own rule): a gradient does not
        survive it.
        """
        flux = np.asarray(self.flux(self.channel, values))
        emitted = self._probe.template.with_values(flux)
        return ModelResult({self.channel: emitted})

    # -- provenance ------------------------------------------------------------

    def describe(self) -> Mapping[str, Any]:
        """The probe's own ``describe()``, plus this backend's identity and any leaf extras."""
        info: dict[str, Any] = {"adapter": "astropy-native", "backend": BACKEND}
        info.update(self._probe.describe())
        extra = {
            str(plan.index): dict(plan.formula.extra) for plan in self._leaves if plan.formula.extra
        }
        if extra:
            info["leaf_configuration"] = extra
        return info

    def __repr__(self) -> str:
        return (
            f"{type(self).__name__}({type(self._probe.astropy_model).__name__}, "
            f"kind={self.kind.__name__}, channel={self.channel!r})"
        )


def _evaluate_leaf(plan: LeafPlan, grid: jax.Array, context: Mapping[str, Any]) -> jax.Array:
    params = dict(
        zip(plan.formula.param_names, (context[name] for name in plan.context_names), strict=True)
    )
    return plan.formula.compute(_OPS, grid, **params)


def from_astropy(
    model: Any,
    *,
    kind: type[FunctionSamples] | None = None,
    priors: Mapping[str, Any] | None = None,
    channel: str = DEFAULT_CHANNEL,
    grid: Any = None,
    output_unit: Any = None,
    equivalencies: Sequence[Any] = (),
) -> Model:
    """A **native** equivalent of an ``astropy.modeling`` model, or a refusal.

    The signature mirrors :func:`ampere.core.from_astropy` exactly, so that
    changing the import is the whole of the change a caller makes; what differs
    is the promise. This one returns a model built from this backend's own
    pieces — differentiable, and usable by ``NUTSEngine`` and ``VIEngine`` —
    or raises. It never returns the black-box adapter.

    Raises
    ------
    CapabilityError
        If *model*, or any component of a compound *model*, is not in
        :data:`TRANSLATIONS`; if a component is tied
        (:func:`~ampere.core.astropy_translations.tied_parameter_refusal`); or
        if a compound model combines two components with an operator outside
        ``+ - * /``
        (:func:`~ampere.core.astropy_translations.unsupported_operator_refusal`).
        The first is caught here, before construction; the other two are
        raised inside :class:`NativeAstropyModel`'s constructor, which this
        propagates.
    """
    components = astropy_components(model)
    if any(type(part) not in TRANSLATIONS for part in components):
        raise translation_refusal(BACKEND, model, TRANSLATIONS)
    return _compose(
        model,
        components,
        kind=kind,
        priors=priors,
        channel=channel,
        grid=grid,
        output_unit=output_unit,
        equivalencies=equivalencies,
    )


def _compose(
    model: Any,
    components: tuple[Any, ...],
    **options: Any,
) -> Model:
    """Build the native model from the curated rows.

    *components* is unused beyond the caller having already checked every one
    of them is in :data:`TRANSLATIONS`; :class:`NativeAstropyModel` re-derives
    the same leaves itself, in the course of also building the combinator
    over them (:func:`~ampere.core.astropy_translations.plan_translation`),
    which needs the tree structure this function's caller does not have.
    """
    del components
    return NativeAstropyModel(model, **options)
