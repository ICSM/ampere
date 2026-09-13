"""The astropy interop adapter: ``DEVELOPMENT_PLAN.md`` §4.7, ``contracts/astropy_compat.md``.

The casual user's route into ampere. Somebody has an ``astropy.modeling`` model
— their own, a compound sum of astropy's library, something a collaborator sent
— and wants to fit it to data with ampere's likelihood. :func:`from_astropy`
wraps it as an :class:`~ampere.core.transform.Model` and asks for nothing else:
the parameters, their bounds, their fixed flags and their ties are already
declared on the astropy side, and this module translates that declaration
rather than making the user restate it.

What the adapter promises, and what it does not
-----------------------------------------------
It promises to evaluate **the user's actual astropy model**. It never
substitutes a lookalike. That is the whole of the 2026-09-01 ruling recorded in
``DEVELOPMENT_PLAN.md`` §2 ("Curated astropy→native translation: opt-in only,
never silent"): a model that changed implementation without its author knowing
is the silent-downgrade class this architecture exists to forbid, and it cuts
in the other direction too — a *faster* substitute is still a different model.
A native torch/jax equivalent is reachable only through the backend-scoped
``from_astropy()`` hook (:func:`translation_refusal` is what that hook refuses
with), which raises for anything not in its curated table.

It does not promise gradients. A Python callable is a black box to torch and
jax, so an adapted model declares ``DIFFERENTIABLE = False`` and
``BACKEND = "reference"``, and a problem composed from it reaches the
**gradient-free engines** — :class:`~ampere.inference.EmceeEngine`,
:class:`~ampere.inference.DynestyEngine`, :class:`~ampere.inference.ZeusEngine`
— and :class:`~ampere.inference.SBIEngine`, which is the engine a wrapped
external model exists for. It never reaches ``NUTSEngine`` or ``VIEngine``:
:func:`~ampere.core.realise` refuses a problem it cannot lower, by name. That
is the honest answer rather than a limitation to work around.

``BATCHABLE`` is ``False``, and the precision matters: astropy models *are*
vectorised over their **input grid** — that is why one evaluation covers the
whole spectrum — but nothing in ``astropy.modeling`` takes a stack of
*parameter* vectors in one call, and ``BATCHABLE`` is a claim about θ
(``transformations.md`` §4). ``simulate_many`` therefore parallelises an
adapted model by process, not by vectorisation, which is exactly what a black
box wants.

Units
-----
Two conversions happen, both **once, at configuration time**, never per
evaluation (``DEVELOPMENT_PLAN.md`` §7's units trap):

* the **input** grid. A model that declares ``input_units`` is handed a
  :class:`~astropy.units.Quantity`, so astropy does its own conversion through
  its ``input_units_equivalencies`` — this is not optional politeness.
  ``BlackBody.input_units`` is ``Hz``, so handing it bare micron numbers
  evaluates a blackbody at 1 to 30 Hz and returns a plausible-looking array. The
  adapter refuses a unitless grid for a model with declared input units rather
  than let that happen. The built Quantity is kept and reused.
* the **output**. The astropy model's own output unit is probed once, and the
  factor taking it to the caller's ``output_unit`` is computed once against the
  grid, checked to be linear (so that hoisting it out of the loop is exact) and
  then multiplied in. Only two equivalencies are ever enabled: the caller's own,
  and :func:`~astropy.units.spectral_density` where the kind has a spectral
  axis — the exact, axis-determined conversion between flux-density
  conventions. **No solid angle is ever invented**: an ``astropy`` surface
  brightness (anything per steradian, which is everything
  ``astropy.modeling.physical_models.BlackBody`` emits) reaches a flux density
  only through an equivalency the caller states, and the refusal says so.

The worked case, which ``tests/core/test_astropy_compat.py`` holds to the
``exact`` tolerance class against :func:`~ampere.backends.reference.planck_jy`:
``BlackBody(temperature=T*u.K, scale=1.0*u.Jy/u.sr)`` emits ``B_nu`` in Jy/sr,
and ``equivalencies=[u.dimensionless_angles()]`` states a solid angle of
**exactly one steradian** — which is the same convention
:class:`ampere.backends.reference.BlackBody` carries in its dimensionless
``scale``, the factor that absorbs the solid angle and the distance dilution
together. The two models are then the same physics written twice, and the
conformance row compares them digit for digit rather than approximately.
"""

from __future__ import annotations

import copy
import math
from collections.abc import Callable, Iterator, Mapping, Sequence
from typing import TYPE_CHECKING, Any, ClassVar

import astropy.units as u
import numpy as np
import scipy.stats

from .exceptions import CapabilityError, CompositionError, ParameterError
from .parameter import Parameter, Value
from .results_schema import (
    DEFAULT_CHANNEL,
    FunctionSamples,
    Image,
    Layout,
    ModelResult,
    Spectrum,
    TimeSeries,
)
from .transform import ChannelRequirements, Model

if TYPE_CHECKING:  # pragma: no cover - typing only
    from numpy.typing import ArrayLike

__all__ = [
    "ADAPTABLE_KINDS",
    "AdaptedAstropyModel",
    "AstropyTie",
    "astropy_components",
    "from_astropy",
    "translation_refusal",
]

#: The :class:`~ampere.core.results_schema.FunctionSamples` kinds this adapter
#: can build, and the number of astropy inputs each one needs. Keyed by kind so
#: that a refusal can list them in a stable order.
#:
#: Deliberately short. These are the kinds an ``astropy.modeling`` model maps
#: onto without the adapter having to invent structure: one input axis for the
#: two one-dimensional point kinds, two for a separable image.
#: :class:`~ampere.core.results_schema.PhotometricPoints` needs a filter list
#: the model does not have, and :class:`~ampere.core.results_schema.Cube` and
#: :class:`~ampere.core.results_schema.VisibilitySet` arrive with their own
#: modality rather than through a generic adapter.
ADAPTABLE_KINDS: Mapping[type[FunctionSamples], int] = {
    Spectrum: 1,
    TimeSeries: 1,
    Image: 2,
}

#: Physical types that identify each one-input kind when the caller lets the
#: adapter infer. Anything else is ambiguous and is refused by name.
_SPECTRAL_TYPES = ("length", "frequency", "energy")
_TIME_TYPES = ("time",)


class AstropyTie:
    """One astropy ``tied`` parameter, recorded as a callable and applied on evaluation.

    ``astropy.modeling`` lets a parameter be *derived*: ``p.tied`` is a callable
    taking the model and returning that parameter's value. This is not
    :class:`~ampere.core.parameter.Tie`, and the difference is worth stating
    because the two words collide. §4.1's ``Tie`` collapses several parameter
    **sites** into one free parameter — an equality, costing one sampler
    dimension and carrying one prior. An astropy tie is an arbitrary Python
    function, so it is neither free (no prior, no dimension), nor fixed (its
    value moves with the fit), nor deferred (no tie group supplies it).
    Registering it as a :class:`~ampere.core.parameter.Parameter` would
    therefore be a lie in every one of the three states §4.1 allows.

    So it is recorded here instead: a **derived** quantity, not a parameter. It
    occupies no sampler dimension, appears in no posterior, and is recomputed
    from the astropy model on every evaluation, in ``param_names`` order —
    which is the order ``astropy.modeling.fitting`` itself applies ties in, and
    means a tie that reads another tied parameter sees that one's value from
    this same evaluation only if it comes earlier in the declaration.

    A Python callable is a black box: an adapted model carrying one is exactly
    as differentiable as one without (which is to say, not at all), so nothing
    about the capability declaration changes.

    Attributes
    ----------
    name
        The astropy parameter this tie determines.
    function
        The callable astropy stores in ``Parameter.tied``; called with the
        astropy model and returning the value.
    """

    __slots__ = ("function", "name")

    def __init__(self, name: str, function: Callable[[Any], Any]) -> None:
        self.name = name
        self.function = function

    def __repr__(self) -> str:
        where = getattr(self.function, "__qualname__", repr(self.function))
        return f"AstropyTie({self.name!r}, {where})"

    def describe(self) -> dict[str, str]:
        """A JSON-normalisable record of this tie, for provenance."""
        return {
            "name": self.name,
            "function": str(getattr(self.function, "__qualname__", repr(self.function))),
        }


# ---------------------------------------------------------------------------
# The adapted model
# ---------------------------------------------------------------------------


class AdaptedAstropyModel(Model):
    """An ``astropy.modeling`` model as an ampere :class:`~ampere.core.Model`.

    Built by :func:`from_astropy`, which is the documented entry point; this
    class is public so that ``isinstance`` and the contract page have something
    to name.

    Capabilities, stated rather than implied: ``DIFFERENTIABLE = False``,
    ``BATCHABLE = False``, ``BACKEND = "reference"``. The gradient-free engines
    and :class:`~ampere.inference.SBIEngine` fit a problem composed from one;
    ``NUTSEngine`` and ``VIEngine`` never will, because there is no gradient to
    take through a Python callable. See the module docstring for why
    ``BATCHABLE`` is ``False`` even though astropy models are vectorised.

    Attributes
    ----------
    astropy_model
        A **deep copy** of the model handed to :func:`from_astropy`. The
        adapter writes parameter values into it on every evaluation, and
        mutating the caller's own object would make a wrapped model and its
        original disagree after the first draw.
    kind
        The :class:`~ampere.core.results_schema.FunctionSamples` subclass the
        single channel emits.
    channel
        The channel name the result is filed under.
    ties
        ``name -> AstropyTie`` for every derived parameter (§ :class:`AstropyTie`).
    """

    #: The four capability flags, declared rather than inherited, as every
    #: piece of the reference path declares them (W2.12).
    DIFFERENTIABLE: ClassVar[bool] = False
    #: Vectorised over the *grid*, never over a stack of θ — see the module
    #: docstring. ``simulate_many`` parallelises this model by process.
    BATCHABLE: ClassVar[bool] = False
    DEVICE: ClassVar[str] = "cpu"
    #: A wrapped astropy model runs on the numpy path, which is what
    #: ``architecture.md`` §2 point 3 calls the reference backend's third job.
    BACKEND: ClassVar[str] = "reference"

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
        astropy_model = _check_astropy_model(model)
        self.astropy_model = copy.deepcopy(astropy_model)
        self.channel = str(channel)
        self.output_unit = _check_output_unit(output_unit)
        self.equivalencies = _flatten_equivalencies(equivalencies)

        # The kind first, and the parameters after it. The order is the order a
        # reader would want the refusals in: what this model *is* comes before
        # what its parameters mean, so a model with no declared kind says so
        # rather than complaining about the first unbounded parameter it meets.
        supplied = _normalise_grid(grid, self.astropy_model.n_inputs)
        self.kind = _resolve_kind(kind, self.astropy_model, supplied)
        _check_kind_arity(self.kind, self.astropy_model)

        self._param_names: tuple[str, ...] = tuple(self.astropy_model.param_names)
        self.ties: dict[str, AstropyTie] = {}
        self._declare(priors or {})

        #: Per-axis coordinates, in declaration order; ``None`` until either the
        #: caller supplies a grid or ``compile_for`` negotiates one.
        self._grids: tuple[Any, ...] | None = None
        self._inputs: tuple[Any, ...] = ()
        self._template: FunctionSamples | None = None
        self._factor: np.ndarray | float = 1.0
        self._raw_unit: u.UnitBase | None = None
        if supplied is not None:
            self._configure(_order_grid(self.kind, supplied))

    # -- declaration -------------------------------------------------------

    def _declare(self, priors: Mapping[str, Any]) -> None:
        """Translate every astropy parameter into an ampere one, or a tie."""
        unknown = sorted(set(priors) - set(self._param_names))
        if unknown:
            raise ParameterError(
                f"from_astropy() was given priors for {unknown}, which "
                f"{type(self.astropy_model).__name__} does not declare. Its parameters are "
                f"{list(self._param_names)}. In a compound model astropy suffixes them by "
                f"submodel ('temperature_0', 'temperature_1'), and the adapter keeps astropy's "
                f"names exactly — renaming them would break the tie callables, which read the "
                f"astropy model by attribute."
            )
        for name in self._param_names:
            declared = getattr(self.astropy_model, name)
            if declared.tied:
                if name in priors:
                    raise ParameterError(
                        f"from_astropy() was given a prior for {name!r}, but astropy declares it "
                        f"tied: its value is computed by {declared.tied!r} on every evaluation, "
                        f"so a prior would never be consulted. Clear the tie on the astropy "
                        f"model, or drop the prior."
                    )
                self.ties[name] = AstropyTie(name, declared.tied)
                continue
            self.register_parameter(self._translate(name, declared, priors))

    def _translate(self, name: str, declared: Any, priors: Mapping[str, Any]) -> Parameter:
        """One astropy ``Parameter`` as an ampere :class:`~ampere.core.Parameter`."""
        unit = declared.unit
        value = np.asarray(declared.value, dtype=float)
        initial: Value = float(value) if value.ndim == 0 else value
        if name in priors:
            return _as_parameter(name, priors[name], unit=unit, value=initial)
        if declared.fixed:
            return Parameter(name, value=initial, fixed=True, unit=unit)
        low, high = declared.bounds
        if low is not None and high is not None:
            lo, hi = float(low), float(high)
            if not (math.isfinite(lo) and math.isfinite(hi)) or hi <= lo:
                raise ParameterError(
                    f"astropy parameter {name!r} declares bounds {declared.bounds}, which are not "
                    f"a usable interval. Fix them on the astropy model, or give this parameter a "
                    f"prior of its own through from_astropy(priors={{{name!r}: ...}})."
                )
            return Parameter(
                name,
                scipy.stats.uniform(lo, hi - lo),
                value=initial,
                unit=unit,
                description=f"uniform over the astropy bounds {declared.bounds}",
            )
        raise ParameterError(
            f"astropy parameter {name!r} is free and unbounded"
            + (f" (astropy declares bounds={declared.bounds})" if any(declared.bounds) else "")
            + ", so there is no prior to translate. ampere will not invent one: a default "
            f"improper prior is how a fit silently becomes a different fit. Give it one — "
            f"from_astropy(model, priors={{{name!r}: scipy.stats.norm(0.0, 1.0)}}) — or bound it "
            f"on the astropy model (parameter.bounds = (lo, hi)), which translates to a uniform "
            f"prior over that interval, or fix it (parameter.fixed = True)."
        )

    # -- configuration -----------------------------------------------------

    def compile_for(self, requirements: Mapping[str, ChannelRequirements]) -> Model:
        """Adopt the negotiated grid for this model's channel.

        The template pattern ``transformations.md`` §14 asks a compiled model
        for, and the one :class:`ampere.backends.reference.BlackBody` and its
        siblings use: one container built here, refilled with ``with_values``
        on every evaluation, so successive results share their axes by
        identity and the unit factor is computed once rather than per draw.

        A channel whose requirements name **some** of this kind's axes but not
        all of them, with no grid already in hand for the rest, is
        :class:`~ampere.core.exceptions.CompositionError` — the loud option
        ruled 2026-09-03 for a model that engages with a requirement it cannot
        honour.
        """
        asked = requirements.get(self.channel)
        if asked is None:
            return self
        axes = tuple(spec.name for spec in self.kind.AXES)
        if not any(axis in asked for axis in axes):
            return self
        existing = self._grids or (None,) * len(axes)
        coordinates: list[Any] = []
        for axis, held in zip(axes, existing, strict=True):
            if axis in asked:
                coordinates.append(asked[axis].coordinates())
            elif held is not None:
                coordinates.append(held)
            else:
                raise CompositionError(
                    f"{type(self).__name__} was asked to compile channel {self.channel!r} for "
                    f"axes {sorted(str(a) for a in axes if a in asked)}, but nothing supplies "
                    f"axis {axis!r}, which a {self.kind.__name__} also needs. Either negotiate "
                    f"every axis of the kind, or give from_astropy(grid=...) coordinates for the "
                    f"ones the instruments do not constrain."
                )
        self._configure(tuple(coordinates))
        return self

    def _configure(self, grids: Sequence[Any]) -> None:
        """Build the input Quantities, the output template and the unit factor. Once."""
        axes = tuple(spec.name for spec in self.kind.AXES)
        prepared = tuple(
            _check_axis(self.kind, spec, coordinates)
            for spec, coordinates in zip(self.kind.AXES, grids, strict=True)
        )
        self._grids = prepared
        self._inputs = _astropy_inputs(self.astropy_model, self.kind, prepared)

        probe = self._call_astropy()
        raw_unit = probe.unit if isinstance(probe, u.Quantity) else None
        bare = np.asarray(probe.value if isinstance(probe, u.Quantity) else probe, dtype=float)
        expected = _expected_shape(self.kind, prepared)
        if bare.shape != expected:
            raise CompositionError(
                f"{type(self.astropy_model).__name__} returned an array of shape {bare.shape} on "
                f"the grid it was given, but a {self.kind.__name__} on these axes needs "
                f"{expected}. A {self.kind.__name__} is "
                f"{'a separable grid' if self.kind.LAYOUT is Layout.GRID else 'a point set'}; "
                f"check the model's n_inputs and n_outputs against the kind you declared."
            )
        self._raw_unit = raw_unit
        self._factor = self._unit_factor(raw_unit, expected, prepared)
        self._template = _build_template(self.kind, axes, prepared, bare, self.output_unit)
        self._register_grids(axes, prepared)

    def _register_grids(self, axes: Sequence[str], grids: Sequence[Any]) -> None:
        """Declare each axis as a buffer, so provenance and ``context`` can see it."""
        for axis, coordinates in zip(axes, grids, strict=True):
            if axis in self.buffers:
                self._buffers = self.buffers.without(axis)
            if axis in self.parameters:
                raise ParameterError(
                    f"{type(self.astropy_model).__name__} declares a parameter called {axis!r}, "
                    f"which is also the name of a {self.kind.__name__} axis, so the grid and the "
                    f"parameter cannot share one context. Rename the astropy parameter, or "
                    f"declare a different kind."
                )
            values = coordinates.value if isinstance(coordinates, u.Quantity) else coordinates
            unit = coordinates.unit if isinstance(coordinates, u.Quantity) else None
            self.register_buffer(axis, np.asarray(values, dtype=float), unit=unit)

    def _unit_factor(
        self,
        raw_unit: u.UnitBase | None,
        shape: tuple[int, ...],
        grids: Sequence[Any],
    ) -> np.ndarray | float:
        """The hoisted conversion from astropy's own output unit to ``output_unit``.

        Computed once, here, and checked to be **linear** in the value — which
        every flux-density equivalency is, and a magnitude one is not — because
        a non-linear conversion cannot be hoisted out of the evaluation loop
        and multiplying by a factor derived from ``1.0`` would silently be a
        different model.
        """
        # Three cases need no factor at all. ``output_unit=None`` keeps whatever
        # astropy returned; a unitless astropy model has its numbers *declared*
        # to be in ``output_unit`` rather than converted into it (said in the
        # page, because it is the one place the adapter takes the caller's word);
        # and a model already emitting the requested unit is the common case.
        if self.output_unit is None or raw_unit is None or raw_unit == self.output_unit:
            return 1.0
        equivalencies = self.equivalencies + _spectral_density_for(self.kind, grids)
        ones = np.ones(shape, dtype=float)
        try:
            factor = np.asarray(
                u.Quantity(ones, raw_unit).to_value(self.output_unit, equivalencies=equivalencies),
                dtype=float,
            )
            doubled = np.asarray(
                u.Quantity(2.0 * ones, raw_unit).to_value(
                    self.output_unit, equivalencies=equivalencies
                ),
                dtype=float,
            )
        except u.UnitConversionError as exc:
            raise CompositionError(
                f"from_astropy() cannot express this model's output, which astropy returns in "
                f"{raw_unit}, in the requested output_unit={self.output_unit}. No equivalency in "
                f"force makes the conversion, and ampere will not invent one. A surface "
                f"brightness (anything per steradian — which is everything astropy's BlackBody "
                f"emits) becomes a flux density only through a solid angle: state it. The usual "
                f"answer is equivalencies=[astropy.units.dimensionless_angles()], which declares "
                f"a solid angle of exactly one steradian and is the convention "
                f"ampere.backends.reference.BlackBody's dimensionless 'scale' carries; put the "
                f"real solid angle in the astropy model's own scale parameter."
            ) from exc
        if not np.allclose(doubled, 2.0 * factor, rtol=0.0, atol=0.0):
            raise CompositionError(
                f"the conversion from {raw_unit} to {self.output_unit} is not linear in the "
                f"value (a magnitude or decibel equivalency, most likely), so it cannot be "
                f"hoisted out of the evaluation loop — and doing the conversion per draw is what "
                f"DEVELOPMENT_PLAN.md §7's units trap forbids. Emit the model in a linear unit "
                f"and convert the posterior afterwards."
            )
        return factor if factor.ndim else float(factor)

    # -- evaluation --------------------------------------------------------

    def _call_astropy(self) -> Any:
        """The astropy model on its configured inputs, with the ties applied."""
        return self.astropy_model(*self._inputs)

    def _apply(self, context: Mapping[str, Value]) -> None:
        """Write this draw's values into the astropy model, then recompute the ties.

        Values are written through ``Parameter.value`` rather than by attribute
        assignment: astropy refuses ``model.temperature = 4000.0`` on a
        parameter that was initialised as a ``Quantity``, and a bare value in
        the parameter's own declared unit is exactly what ampere holds.
        """
        for name in self._param_names:
            if name in self.ties:
                continue
            value = context.get(name)
            if value is not None:
                getattr(self.astropy_model, name).value = value
        for name, tie in self.ties.items():
            getattr(self.astropy_model, name).value = tie.function(self.astropy_model)

    def evaluate(self, **values: Value) -> ModelResult:
        """Evaluate the wrapped astropy model and file the result under this channel."""
        if self._template is None:
            raise CompositionError(
                f"{type(self).__name__} has no grid to evaluate on. Either give one — "
                f"from_astropy(model, grid=wavelength * u.micron, ...) — or compose it into a "
                f"FittingProblem, whose negotiation calls compile_for() with the coordinates the "
                f"instruments asked for."
            )
        self._apply(self.context(values))
        raw = self._call_astropy()
        if isinstance(raw, u.Quantity):
            if raw.unit != self._raw_unit:
                raise CompositionError(
                    f"{type(self.astropy_model).__name__} returned {raw.unit} on this draw but "
                    f"{self._raw_unit} when it was configured. The adapter hoists the unit "
                    f"conversion out of the evaluation loop, so a model whose output unit "
                    f"depends on its parameter values cannot be wrapped."
                )
            bare = np.asarray(raw.value, dtype=float)
        else:
            bare = np.asarray(raw, dtype=float)
        return ModelResult({self.channel: self._template.with_values(bare * self._factor)})

    # -- provenance --------------------------------------------------------

    def describe(self) -> Mapping[str, Any]:
        """What this wrapper computes, beyond its parameters and buffers.

        The astropy class and its submodels, the parameter names in astropy's
        own order, the derived (tied) parameters, the declared kind and channel
        and the unit convention. All of it changes what the model computes, and
        none of it is a parameter or a buffer, so without this hook two fits of
        two different astropy models would share a cache key (``results.md``
        §13.13's limitation 13).
        """
        return {
            "adapter": "astropy",
            "astropy_class": type(self.astropy_model).__name__,
            "astropy_components": [type(part).__name__ for part in astropy_components(self)],
            "astropy_parameters": list(self._param_names),
            "astropy_ties": [self.ties[name].describe() for name in sorted(self.ties)],
            "kind": self.kind.__name__,
            "channel": self.channel,
            "output_unit": None if self.output_unit is None else str(self.output_unit),
            "raw_unit": None if self._raw_unit is None else str(self._raw_unit),
            "equivalencies": len(self.equivalencies),
        }

    def __repr__(self) -> str:
        tied = f", ties={sorted(self.ties)}" if self.ties else ""
        return (
            f"{type(self).__name__}({type(self.astropy_model).__name__}, "
            f"kind={self.kind.__name__}, channel={self.channel!r}{tied})"
        )


# ---------------------------------------------------------------------------
# The entry point
# ---------------------------------------------------------------------------


def from_astropy(
    model: Any,
    *,
    kind: type[FunctionSamples] | None = None,
    priors: Mapping[str, Any] | None = None,
    channel: str = DEFAULT_CHANNEL,
    grid: Any = None,
    output_unit: Any = None,
    equivalencies: Sequence[Any] = (),
) -> AdaptedAstropyModel:
    """Wrap an ``astropy.modeling`` model as an ampere :class:`~ampere.core.Model`.

    ``DEVELOPMENT_PLAN.md`` §4.7 and ``docs/design/contracts/astropy_compat.md``.
    Compound models included: astropy's own ``param_names`` are kept exactly, so
    a sum of two blackbodies declares ``temperature_0`` and ``temperature_1``
    and the posterior says so.

    Parameters
    ----------
    model
        Any ``astropy.modeling.Model`` instance with one output and one or two
        inputs. Deep-copied, never mutated.
    kind
        The :class:`~ampere.core.results_schema.FunctionSamples` subclass the
        channel emits — see :data:`ADAPTABLE_KINDS`. Inferred **only** when the
        grid's own units make it unambiguous (a length, frequency or energy
        axis is a :class:`~ampere.core.results_schema.Spectrum`; a time axis a
        :class:`~ampere.core.results_schema.TimeSeries`; two axes an
        :class:`~ampere.core.results_schema.Image`), and refused by name
        otherwise. It is never guessed from the model's class name.
    priors
        ``astropy parameter name -> prior``, overriding the translation below.
        A frozen ``scipy.stats`` distribution fits it, a number holds it fixed,
        a ready-made :class:`~ampere.core.Parameter` is used as it stands.
    channel
        Name of the emitted channel. The default is
        :data:`~ampere.core.results_schema.DEFAULT_CHANNEL`.
    grid
        Coordinates the model is evaluated on: an array or
        :class:`~astropy.units.Quantity` for a one-axis kind, a mapping of axis
        name to coordinates (or a sequence in axis order) for two. Optional —
        without it the model has no grid until ``compile_for`` negotiates one,
        which is the normal path inside a
        :class:`~ampere.core.FittingProblem`.
    output_unit
        The unit the emitted container carries. ``None`` keeps whatever astropy
        returns. See the module docstring for how the conversion is hoisted and
        what is never invented.
    equivalencies
        Extra ``astropy.units`` equivalencies for that conversion.
        :func:`~astropy.units.spectral_density` is already in force for a kind
        with a spectral axis.

    Returns
    -------
    AdaptedAstropyModel
        Black-box on the reference backend: ``DIFFERENTIABLE = False``,
        ``BATCHABLE = False``, ``BACKEND = "reference"``. Gradient-free engines
        and :class:`~ampere.inference.SBIEngine` fit it; ``NUTSEngine`` and
        ``VIEngine`` refuse it by name, and that refusal is the contract rather
        than a gap.

    Notes
    -----
    **Parameter translation**, in the order the adapter tries it:

    ============================ =========================================
    astropy declaration          ampere parameter
    ============================ =========================================
    ``priors[name]`` given       whatever the caller passed (wins outright)
    ``tied`` (a callable)        **not** a parameter — an :class:`AstropyTie`
    ``fixed = True``             frozen at its value
    finite ``bounds``            ``scipy.stats.uniform(lo, hi - lo)``
    free and unbounded           refused, by name
    ============================ =========================================

    A parameter's unit, and therefore the unit its bounds and its prior are
    declared in, is astropy's own — ``bounds`` are always bare numbers in it.
    ``BlackBody.temperature`` ships with ``bounds = (0, None)``, which is
    half-open and therefore *not* a uniform prior: it is refused, by name, with
    the three ways to give it one.

    Examples
    --------
    >>> import astropy.units as u, numpy as np, scipy.stats as st
    >>> from astropy.modeling.models import PowerLaw1D
    >>> grid = np.array([1.0, 2.0, 4.0]) * u.micron
    >>> power_law = PowerLaw1D(amplitude=2.0 * u.Jy, x_0=1.0 * u.micron, alpha=1.0)
    >>> model = from_astropy(
    ...     power_law,
    ...     grid=grid,
    ...     priors={"amplitude": st.loguniform(0.1, 10.0), "x_0": 1.0, "alpha": st.norm(1.0, 0.5)},
    ...     output_unit=u.Jy,
    ... )
    >>> model.kind.__name__, model.parameters.free_names
    ('Spectrum', ('amplitude', 'alpha'))
    >>> model(amplitude=2.0, alpha=1.0)["default"].values.tolist()
    [2.0, 1.0, 0.5]
    >>> model.DIFFERENTIABLE, model.BACKEND
    (False, 'reference')
    """
    return AdaptedAstropyModel(
        model,
        kind=kind,
        priors=priors,
        channel=channel,
        grid=grid,
        output_unit=output_unit,
        equivalencies=equivalencies,
    )


# ---------------------------------------------------------------------------
# The opt-in translation hook (the backends' half)
# ---------------------------------------------------------------------------


def astropy_components(model: Any) -> tuple[Any, ...]:
    """The leaf submodels of *model*, or ``(model,)`` if it is not compound.

    An :class:`AdaptedAstropyModel` is unwrapped to the astropy model it holds,
    so a caller may pass either. A compound model's leaves are what a
    translation table is keyed on: a sum of a blackbody and a power law is
    translatable exactly when **both** halves are, which is the question
    :func:`translation_refusal` answers.
    """
    astropy_model = model.astropy_model if isinstance(model, AdaptedAstropyModel) else model
    if getattr(astropy_model, "n_submodels", 1) <= 1:
        return (astropy_model,)
    return tuple(_iter_leaves(astropy_model))


def _iter_leaves(compound: Any) -> Iterator[Any]:
    for part in compound:
        if getattr(part, "n_submodels", 1) > 1:
            yield from _iter_leaves(part)
        else:
            yield part


def translation_refusal(
    backend: str,
    model: Any,
    table: Mapping[type, Any],
) -> CapabilityError:
    """The refusal a backend's ``from_astropy()`` raises, naming what it cannot translate.

    The 2026-09-01 ruling in ``DEVELOPMENT_PLAN.md`` §2 in one function:
    translation to a native equivalent is **opt-in and never silent**, so a
    backend asked for a native model it has no curated row for raises rather
    than quietly handing back the black-box adapter. Silently degrading would
    give the caller a model that is not differentiable when they asked for the
    one that is; silently substituting would give them a model they did not
    write.

    Shared by every backend's hook because the refusal must read the same on
    all of them — one name, one remedy, one place to change it — and because
    the leaf decomposition is astropy's business rather than any backend's.

    Parameters
    ----------
    backend
        The backend's one name (W2.12), for the message.
    model
        The astropy model (or an :class:`AdaptedAstropyModel`) asked for.
    table
        That backend's curated translation table, keyed on astropy model class.
        **Empty in W4.6**: this item lands the hook, W4.7 lands the rows.

    Returns
    -------
    CapabilityError
        Raised by the caller, so the traceback points at the backend's hook.

        :class:`~ampere.core.exceptions.CapabilityError` rather than
        :class:`~ampere.core.exceptions.LoweringError`, and the choice is worth
        a sentence. Nothing about the caller's model is malformed — the
        black-box route consumes it happily — and one specific path cannot
        serve it, which is precisely what ``CapabilityError`` is for; it is
        also a :class:`NotImplementedError`, which is exactly what "W4.7 has
        not landed this row yet" means. ``LoweringError`` would have been the
        other candidate, but its constructor is prior-family-shaped
        (``LoweringError(family, backend=...)``, rendering "prior family
        'BlackBody'"), and bending it to fit would mean editing a frozen
        contract module to say something it was not written to say.
    """
    missing = sorted(
        {type(part).__name__ for part in astropy_components(model) if type(part) not in table}
    )
    rows = sorted(cls.__name__ for cls in table)
    holds = ", ".join(rows) if rows else "nothing yet (W4.7 adds the curated rows)"
    return CapabilityError(
        f"ampere.backends.{backend}.from_astropy() has no native translation for "
        f"{', '.join(missing)}. Translation is opt-in and never silent (DEVELOPMENT_PLAN.md §2, "
        f"ruled 2026-09-01), so this refuses rather than falling back to the black-box adapter: "
        f"you asked for a differentiable model and would have been given one that is not. Use "
        f"ampere.core.from_astropy(), which evaluates your actual astropy model on the reference "
        f"path — gradient-free engines and SBI, never NUTS or VI — or extend the table, which "
        f"currently holds: {holds}."
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _check_astropy_model(model: Any) -> Any:
    """*model* is a single, one-output ``astropy.modeling`` model, or a named refusal.

    The import is **lazy and local**, which is the rule ``architecture.md`` §4
    sets for anything that need not be paid for on ``import ampere.core``:
    ``astropy.modeling`` is a base dependency, so this is a cost question
    rather than an availability one, and ``tests/core/test_astropy_compat.py``
    measures that importing the core package still does not pull it in.
    """
    from astropy.modeling import Model as AstropyModel

    if not isinstance(model, AstropyModel):
        raise CompositionError(
            f"from_astropy() takes an astropy.modeling.Model, got {type(model).__name__}. A "
            f"plain Python callable is not one: write it as an ampere Model directly — that is "
            f"a dozen lines (see ampere.core.Model) and it is the honest declaration."
        )
    if len(model) != 1 or model.model_set_axis is not False:
        raise CompositionError(
            f"from_astropy() was given an astropy model *set* of {len(model)} models. A "
            f"ModelResult channel is one model's output; fit them one at a time, or compose the "
            f"set as several channels by wrapping each member separately."
        )
    if model.n_outputs != 1:
        raise CompositionError(
            f"{type(model).__name__} has n_outputs={model.n_outputs}. from_astropy() builds one "
            f"channel from one output; a multi-output model needs one channel per output, and "
            f"the adapter will not choose the kinds or the names for you. Wrap the outputs "
            f"separately, or write the model against ampere.core.Model, which returns a "
            f"ModelResult of as many named channels as it likes."
        )
    if model.n_inputs not in set(ADAPTABLE_KINDS.values()):
        raise CompositionError(
            f"{type(model).__name__} has n_inputs={model.n_inputs}. from_astropy() builds "
            f"{_kind_list()}, which need one or two input axes."
        )
    return model


def _kind_list() -> str:
    return ", ".join(
        f"{kind.__name__} ({count} axis)" if count == 1 else f"{kind.__name__} ({count} axes)"
        for kind, count in ADAPTABLE_KINDS.items()
    )


def _flatten_equivalencies(given: Sequence[Any]) -> list[tuple[Any, ...]]:
    """``equivalencies=`` as the flat list of tuples astropy wants.

    ``astropy.units.Equivalency`` is itself a list of tuples, so both
    ``equivalencies=u.dimensionless_angles()`` and the more natural-looking
    ``equivalencies=[u.dimensionless_angles()]`` arrive here; astropy accepts
    only the first shape, and the error it gives for the second ("Invalid
    equivalence entry 0") names neither the argument nor the fix. One level of
    flattening accepts both.
    """
    flat: list[tuple[Any, ...]] = []
    for item in given:
        if isinstance(item, tuple):
            flat.append(item)
        else:
            flat.extend(item)
    return flat


def _check_output_unit(output_unit: Any) -> u.UnitBase | None:
    if output_unit is None:
        return None
    try:
        return u.Unit(output_unit)
    except (TypeError, ValueError) as exc:
        raise CompositionError(
            f"from_astropy(output_unit={output_unit!r}) is not an astropy unit."
        ) from exc


def _as_parameter(name: str, given: Any, *, unit: Any, value: Value) -> Parameter:
    """A ``priors=`` entry as a :class:`~ampere.core.Parameter`.

    The same three-way coercion ``ampere.backends.reference._declare`` offers,
    deliberately not imported from there: ``ampere.core`` may not import a
    backend (``architecture.md`` §4), and the messages here talk about
    ``priors=`` entries rather than about constructor arguments.
    """
    if isinstance(given, Parameter):
        if given.name != name:
            raise ParameterError(
                f"priors[{name!r}] is a Parameter named {given.name!r}. The adapter keeps "
                f"astropy's own parameter names, because the tie callables read the astropy "
                f"model by attribute; rename it with .rename({name!r})."
            )
        return given
    if isinstance(given, (int, float, np.floating, np.integer)) and not isinstance(given, bool):
        return Parameter(name, value=float(given), fixed=True, unit=unit)
    if hasattr(given, "ppf"):
        return Parameter(name, given, value=value, unit=unit)
    raise ParameterError(
        f"priors[{name!r}] must be a frozen scipy.stats distribution (to fit it), a number (to "
        f"hold it fixed), or an ampere Parameter — got {type(given).__name__}."
    )


def _normalise_grid(grid: Any, n_inputs: int) -> tuple[Any, ...] | Mapping[str, Any] | None:
    """``grid=`` as one coordinate set per axis — a mapping is kept, and ordered later.

    A mapping is not resolved here because the axis *order* is the kind's, and
    the kind may not be known yet: the grid's own unit is what
    :func:`_inference_reason` reads to infer it. :func:`_order_grid` does the
    ordering once the kind is settled.
    """
    if grid is None:
        return None
    if isinstance(grid, Mapping):
        return grid
    if isinstance(grid, (u.Quantity, np.ndarray)) or not isinstance(grid, Sequence):
        return (grid,)
    values = tuple(grid)
    if n_inputs == 1 or all(isinstance(item, (int, float)) for item in values):
        return (np.asarray(values, dtype=float),)
    return values


def _order_grid(
    kind: type[FunctionSamples], supplied: tuple[Any, ...] | Mapping[str, Any]
) -> tuple[Any, ...]:
    """The supplied coordinates in the kind's own axis order, arity checked."""
    axes = tuple(str(spec.name) for spec in kind.AXES)
    if isinstance(supplied, Mapping):
        missing = [axis for axis in axes if axis not in supplied]
        extra = sorted(set(supplied) - set(axes))
        if missing or extra:
            raise CompositionError(
                f"from_astropy(grid=...) for a {kind.__name__} needs one entry per axis, "
                f"{list(axes)}; it was "
                + (f"missing {missing}" if missing else "")
                + (" and " if missing and extra else "")
                + (f"given unknown {extra}" if extra else "")
                + "."
            )
        return tuple(supplied[axis] for axis in axes)
    if len(supplied) != len(axes):
        raise CompositionError(
            f"from_astropy(grid=...) supplied {len(supplied)} coordinate set(s), but a "
            f"{kind.__name__} has {len(axes)} axes {list(axes)}. Pass a mapping keyed on the "
            f"axis names to be unambiguous."
        )
    return tuple(supplied)


def _first_coordinate(supplied: tuple[Any, ...] | Mapping[str, Any] | None) -> Any:
    """One coordinate set to read an axis unit from, for kind inference."""
    if supplied is None:
        return None
    if isinstance(supplied, Mapping):
        return next(iter(supplied.values()), None)
    return supplied[0] if supplied else None


def _resolve_kind(
    kind: type[FunctionSamples] | None,
    model: Any,
    supplied: tuple[Any, ...] | Mapping[str, Any] | None,
) -> type[FunctionSamples]:
    """The declared kind, or the inferred one, or the refusal."""
    if kind is not None:
        if kind not in ADAPTABLE_KINDS:
            raise CompositionError(
                f"from_astropy(kind={kind.__name__ if isinstance(kind, type) else kind!r}) is not "
                f"a kind this adapter builds. It builds {_kind_list()}. Other kinds arrive with "
                f"their own modality rather than through a generic adapter — a "
                f"PhotometricPoints needs a filter list the astropy model does not have."
            )
        return kind
    reason = _inference_reason(model, supplied)
    if isinstance(reason, str):
        raise CompositionError(
            f"from_astropy() cannot infer the ModelResult kind for "
            f"{type(model).__name__}: {reason}. Declare it — "
            f"from_astropy(model, kind=Spectrum, ...) — choosing from {_kind_list()}. The kind "
            f"is never guessed from the model's class name: §4.7 makes it the caller's "
            f"declaration, and the wrong kind is a fit that runs and is silently wrong."
        )
    return reason


def _inference_reason(
    model: Any, supplied: tuple[Any, ...] | Mapping[str, Any] | None
) -> type[FunctionSamples] | str:
    """The inferred kind, or the sentence saying why there isn't one."""
    if supplied is None:
        return (
            "no grid was given, so there are no axis units to read, and the negotiated grid is "
            "not known until compile_for()"
        )
    if model.n_inputs == 2:
        return Image
    first = _first_coordinate(supplied)
    unit = first.unit if isinstance(first, u.Quantity) else None
    if unit is None:
        return "the grid carries no unit, and the axis unit is the only unambiguous evidence"
    physical = u.get_physical_type(unit)
    if any(physical == name for name in _SPECTRAL_TYPES):
        return Spectrum
    if any(physical == name for name in _TIME_TYPES):
        return TimeSeries
    return (
        f"the grid is in {unit} (physical type {physical}), which is neither a spectral axis "
        f"(length, frequency or energy) nor a time axis"
    )


def _check_kind_arity(kind: type[FunctionSamples], model: Any) -> None:
    wanted = ADAPTABLE_KINDS[kind]
    if model.n_inputs != wanted:
        raise CompositionError(
            f"a {kind.__name__} has {wanted} coordinate "
            f"{'axis' if wanted == 1 else 'axes'}, but {type(model).__name__} takes "
            f"n_inputs={model.n_inputs}. Declare the kind that matches the model's arity "
            f"({_kind_list()})."
        )


def _check_axis(kind: type[FunctionSamples], spec: Any, coordinates: Any) -> Any:
    """One axis's coordinates, with the unit rule the kind declares enforced here."""
    if spec.unit_required and not isinstance(coordinates, u.Quantity):
        raise CompositionError(
            f"a {kind.__name__}'s {spec.name!r} axis must carry a unit, and the grid given for "
            f"it does not. Pass a Quantity — grid=wavelength * astropy.units.micron. This is "
            f"not pedantry: an astropy model that declares input_units converts what it is "
            f"given, so handing BlackBody (whose input_units is Hz) bare micron numbers "
            f"evaluates a blackbody between 1 and 30 Hz and returns a plausible array."
        )
    return coordinates


def _astropy_inputs(
    model: Any,
    kind: type[FunctionSamples],
    grids: Sequence[Any],
) -> tuple[Any, ...]:
    """The arguments the astropy model is called with, built once.

    A model that declares ``input_units`` is handed Quantities, so astropy runs
    its own ``input_units_equivalencies`` (which is how ``BlackBody`` accepts a
    wavelength at all). A model that declares none is handed bare numbers **in
    the grid's own unit**, which the contract page states, because astropy's
    unitless models produce nonsense units from a Quantity input rather than
    refusing.
    """
    declared = model.input_units
    if kind.LAYOUT is Layout.GRID:
        bare = [np.asarray(g.value if isinstance(g, u.Quantity) else g, dtype=float) for g in grids]
        meshed = np.meshgrid(*bare, indexing="ij")
        if declared:
            return tuple(
                u.Quantity(values, grid.unit) if isinstance(grid, u.Quantity) else values
                for values, grid in zip(meshed, grids, strict=True)
            )
        return tuple(meshed)
    if declared:
        for name, grid in zip(model.inputs, grids, strict=True):
            if declared.get(name) is not None and not isinstance(grid, u.Quantity):
                raise CompositionError(
                    f"{type(model).__name__} declares input_units[{name!r}] = {declared[name]}, "
                    f"so it must be given a Quantity: a bare array is interpreted as already "
                    f"being in {declared[name]}, which silently evaluates the model somewhere "
                    f"else entirely. Pass grid=coordinates * <unit>."
                )
        return tuple(grids)
    return tuple(
        np.asarray(g.value if isinstance(g, u.Quantity) else g, dtype=float) for g in grids
    )


def _expected_shape(kind: type[FunctionSamples], grids: Sequence[Any]) -> tuple[int, ...]:
    sizes = tuple(int(np.size(g)) for g in grids)
    return sizes if kind.LAYOUT is Layout.GRID else (sizes[0],)


def _build_template(
    kind: type[FunctionSamples],
    axes: Sequence[str],
    grids: Sequence[Any],
    values: ArrayLike,
    output_unit: u.UnitBase | None,
) -> FunctionSamples:
    """The container this model refills on every evaluation."""
    built = dict(zip(axes, grids, strict=True))
    coordinates = [built[axis] for axis in axes]
    empty = np.zeros_like(np.asarray(values))
    return kind(*coordinates, empty, unit=output_unit)


def _spectral_density_for(kind: type[FunctionSamples], grids: Sequence[Any]) -> list[Any]:
    """``spectral_density`` against the spectral axis, where the kind has one.

    The one equivalency the adapter enables without being asked, because it is
    the exact, axis-determined conversion between per-frequency and
    per-wavelength flux densities rather than an assumption about the source.
    Nothing else is ever added.
    """
    for spec, grid in zip(kind.AXES, grids, strict=True):
        if spec.name == "spectral_axis" and isinstance(grid, u.Quantity):
            return list(u.spectral_density(grid))
    return []
