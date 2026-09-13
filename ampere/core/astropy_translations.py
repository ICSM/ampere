"""The curated astropy→native table's physics: ``contracts/astropy_compat.md`` §5, W4.7.

``ampere.core.astropy_compat`` §5 lands the opt-in hook and an empty table;
this module is what W4.7 fills it with. It is backend-neutral in the sense
``ampere/core/kernels.py``'s :class:`~ampere.core.kernels.ArrayOps` already
established for the GP kernels (W4.5): the mathematics is written **once**,
against a two-method namespace (:class:`TranslationOps`) plus each parameter's
own array type, and a backend supplies the namespace. It turns out to need far
less than :class:`~ampere.core.kernels.ArrayOps` does, because every curated
model here is a plain elementwise function of one coordinate axis, and
``torch.Tensor``/``jax.Array`` both already overload ``+ - * / **`` to do the
right thing with a scalar, a numpy array or another array of the same kind —
so a formula written in ordinary Python arithmetic needs no namespace at all
for that half. What genuinely differs between the two libraries is
``exp``/``where``, and the Planck function, which each backend already owns as
its own tested ``planck_jy`` (:mod:`ampere.backends.torch.models`,
:mod:`ampere.backends.jax.models`) and supplies here through
:meth:`TranslationOps.blackbody` rather than have a third copy of the physics
written against this namespace.

**Six curated leaves** (``DEVELOPMENT_PLAN.md`` §4.7, ``WORK_ITEMS.md``
W4.7): :class:`~astropy.modeling.functional_models.Gaussian1D`,
:class:`~astropy.modeling.functional_models.Const1D`,
:class:`~astropy.modeling.powerlaws.PowerLaw1D`,
:class:`~astropy.modeling.powerlaws.BrokenPowerLaw1D`,
:class:`~astropy.modeling.polynomial.Polynomial1D` and
:class:`~astropy.modeling.physical_models.BlackBody`. :data:`LEAF_BUILDERS`
keys astropy's own classes, so both
:mod:`ampere.backends.torch.astropy` and :mod:`ampere.backends.jax.astropy`
consult the identical table — "one table serves both backends", which the
item's design guidance invites when the physics fits the namespace, and it
does here.

Composition follows astropy's own expression tree
--------------------------------------------------
A compound model is a binary tree of ``.left``, ``.right`` and ``.op``
(``astropy.modeling.core.CompoundModel``); :func:`plan_translation` walks it
once, at construction time, assigning each leaf the same numeric suffix
astropy's own ``param_names`` would (leaf order is left-to-right depth-first
in both, so the two numberings agree without either reading the other), and
builds one Python closure that combines the leaves' evaluated arrays with
``+``, ``-``, ``*`` or ``/`` — exactly the operators astropy's own tree uses
for these — following the same tree. Astropy's other composition operators,
``|`` (output piped into the next model's input) and ``&`` (concatenated
inputs), have no elementwise reading as *one* channel's flux and are refused
by name (:func:`unsupported_operator_refusal`); none of the six curated leaves
is more than one input and one output, so ``|``/``&`` can only ever appear
*between* two already-translatable sums, which is exactly where this catches
them.

Units: a leaf's own declared input unit, honoured once
--------------------------------------------------------
Most of the six take a bare-number grid in whatever unit the caller
negotiated — the case ``docs/design/contracts/astropy_compat.md`` §4.1 calls
"a model that declares no input units". A model built with ``Quantity``
parameters is the exception (``PowerLaw1D(x_0=1.0*u.nm, ...)`` declares
``input_units={'x': u.nm}``), and :func:`leaf_grid_values` converts the
negotiated grid into that unit once, using the leaf's own
``input_units_equivalencies`` — never per evaluation
(``DEVELOPMENT_PLAN.md`` §7). :class:`~astropy.modeling.physical_models.BlackBody`
is a third case, handled explicitly rather than generically: astropy declares
its input unit as Hz, but the backends' own ``planck_jy`` is written in
wavelength (micron, the convention ``ampere.backends.reference.planck_jy``
and W4.6's exact-tolerance row already fix), so :data:`LEAF_BUILDERS` sends it
micron with :func:`astropy.units.spectral` in force — the same physical
conversion astropy's own Hz route makes, reached by feeding a *different*
but physically identical implementation of the same law, and cheap to state
because the conversion happens once, on the buffer, never in the evaluation
loop. Its *output* needs the same care: astropy computes the Planck function
once, internally, in a fixed CGS-style convention and only then performs a
genuine unit conversion into ``scale``'s own declared unit — dimensionless
leaves it in that CGS unit, ``scale=1.0*u.Jy/u.sr`` converts it properly into
Jy/sr — so a bare ``scale * planck_jy(...)`` would not reproduce astropy's
own raw numeric value the way every other curated leaf's bare
parameter-times-formula does. :func:`LEAF_BUILDERS`'s BlackBody row probes
this leaf's own raw unit once, at translation time, and folds in the
constant Jy-to-that-unit correction, so this leaf satisfies the same "my raw
number is astropy's raw number, in astropy's own raw unit" invariant as
every other one — which is what lets the backend hook's single, whole-model
output-unit factor (computed once, against the *real* astropy model) compose
correctly leaf by leaf, whether or not the caller declared an
``output_unit=`` at all.

What is refused, and why
-------------------------
A **tied** parameter (:func:`tied_parameter_refusal`): astropy's ``tied=`` is
an arbitrary Python callable of the whole model, evaluated on plain floats
read off astropy's own ``Parameter`` objects
(``ampere.core.astropy_compat.AstropyTie``). There is no gradient through it
on either native backend, and on jax there is no way to evaluate it at all
while a value is a traced abstract array rather than a concrete float — a
Python-level ``float(traced_value)`` is exactly the "printing a tracer"
mistake jax raises ``ConcretizationTypeError`` for. Rather than support it on
one backend and not the other (a capability that would then depend on which
backend a user happened to pick, which is precisely the kind of surprise
``architecture.md`` commits this project against), a compound model carrying
a tie is refused on **both** native backends, by name, pointing at
``ampere.core.from_astropy()`` — the black-box route this project already
holds to applying a tie exactly as astropy defines it.
"""

from __future__ import annotations

import dataclasses
import operator
from collections.abc import Callable, Mapping, Sequence
from typing import Any, Protocol

import astropy.units as u
import numpy as np
from astropy.modeling.models import (
    BlackBody,
    BrokenPowerLaw1D,
    Const1D,
    Gaussian1D,
    Polynomial1D,
    PowerLaw1D,
)

from .exceptions import CapabilityError

__all__ = [
    "LEAF_BUILDERS",
    "LeafFormula",
    "LeafPlan",
    "TranslationOps",
    "leaf_grid_values",
    "plan_translation",
    "tied_parameter_refusal",
    "unsupported_operator_refusal",
]


class TranslationOps(Protocol):
    """The two elementwise operations, plus the Planck function, a curated formula needs.

    Deliberately this small — see the module docstring for why ordinary
    ``+ - * / **`` needs no namespace at all here, unlike
    :class:`~ampere.core.kernels.ArrayOps`, which exists because a kernel's
    algebra is richer.
    """

    def exp(self, array: Any) -> Any:
        """Elementwise exponential."""
        ...

    def where(self, condition: Any, if_true: Any, if_false: Any) -> Any:
        """Elementwise select, both branches already evaluated (traceable on both backends)."""
        ...

    def blackbody(self, wavelength: Any, temperature: Any) -> Any:
        """``B_nu(temperature)`` in Jy/sr, for *wavelength* already in this backend's own micron.

        Delegates to this backend's own ``planck_jy`` — the module docstring's
        third case — rather than restate the physics against this namespace.
        """
        ...


@dataclasses.dataclass(frozen=True)
class LeafFormula:
    """One curated astropy leaf type, translated to a pure array function.

    Attributes
    ----------
    param_names
        This leaf's **own** (unsuffixed) astropy parameter names, in the order
        :attr:`compute` takes them as keywords. :func:`plan_translation`
        applies the compound model's numeric suffix, if any.
    compute
        ``(ops, grid, **params) -> array``: this leaf's contribution, in
        astropy's own declared unit for the parameters it was given.
    grid_unit
        The unit :attr:`compute` needs the grid in, or ``None`` to use the
        negotiated grid's own unit unconverted (``astropy_compat.md`` §4.1's
        "no input units" rule).
    grid_equivalencies
        Flattened equivalency tuples for that conversion (already flat, in the
        shape :meth:`astropy.units.Quantity.to_value` accepts directly — see
        ``ampere.core.astropy_compat._flatten_equivalencies`` for why a nested
        ``Equivalency`` cannot be passed as one list entry).
    extra
        Non-parameter configuration that changes what :attr:`compute`
        computes (:class:`~astropy.modeling.polynomial.Polynomial1D`'s
        ``domain``/``window`` rescaling), folded into the composed model's
        ``describe()`` so two leaves differing only here do not share a cache
        key (``results.md`` §13.13's limitation 13). Empty for a leaf with
        none.
    """

    param_names: tuple[str, ...]
    compute: Callable[..., Any]
    grid_unit: u.UnitBase | None = None
    grid_equivalencies: tuple[Any, ...] = ()
    extra: Mapping[str, Any] = dataclasses.field(default_factory=dict)


def _generic_grid_unit(leaf: Any) -> tuple[u.UnitBase | None, tuple[Any, ...]]:
    """The unit (and equivalencies) *leaf* itself declares for its one input, if any.

    ``None`` for the common case — a leaf built from bare numbers, which
    astropy leaves ``input_units`` unset for. A leaf built with ``Quantity``
    parameters (``PowerLaw1D(x_0=1.0*u.nm, ...)``) declares one, and this is
    the same fact ``ampere.core.astropy_compat._astropy_inputs`` reads for the
    black-box route — read here rather than re-derived, so the two agree.
    """
    declared = leaf.input_units
    if not declared:
        return None, ()
    (input_name,) = leaf.inputs
    unit = declared.get(input_name)
    if unit is None:
        return None, ()
    equivalencies = tuple((leaf.input_units_equivalencies or {}).get(input_name, ()))
    return unit, equivalencies


def _blackbody(leaf: Any) -> LeafFormula:
    # Every other curated leaf's compute() reproduces astropy's raw numeric
    # value in astropy's own raw unit simply by using each Parameter's bare
    # value (which *is* that value, in its own declared unit) — the invariant
    # a compound expression's single, whole-model output-unit factor (the
    # probe's own, computed once against the *real* astropy model) needs to
    # apply correctly leaf by leaf. BlackBody is the exception: astropy
    # computes the Planck function once, internally, in a fixed CGS-style
    # convention, and only *then* performs a genuine unit conversion into
    # `scale`'s own declared unit (dimensionless leaves it in that CGS unit;
    # `scale=1.0*u.Jy/u.sr` converts it properly into Jy/sr) — so a bare
    # multiplication by `scale` alone would not reproduce it.
    #
    # `planck_jy` computes the same physics directly in Jy (the
    # `dimensionless_angles()` one-steradian convention). Probing this leaf's
    # own raw unit once, here, at translation time (the *unit* only — never
    # the value, so what `scale` happens to be worth does not matter), and
    # converting Jy to it once gives a constant correction that makes this
    # leaf's compute() satisfy the same invariant as every other one, so the
    # probe's factor composes correctly whether or not the caller declared an
    # ``output_unit=`` at all.
    raw_unit = leaf(1.0 * u.Hz).unit
    to_raw_unit = u.Jy.to(raw_unit, equivalencies=list(u.dimensionless_angles()))

    def compute(ops: TranslationOps, grid: Any, *, temperature: Any, scale: Any) -> Any:
        return scale * ops.blackbody(grid, temperature) * to_raw_unit

    # Explicit, not _generic_grid_unit: astropy declares BlackBody's own input
    # unit as Hz, but the backends' planck_jy is written in wavelength, so the
    # grid this leaf wants is micron with `spectral` in force — see the module
    # docstring's third case.
    return LeafFormula(
        ("temperature", "scale"),
        compute,
        grid_unit=u.micron,
        grid_equivalencies=tuple(u.spectral()),
    )


def _power_law(leaf: Any) -> LeafFormula:
    unit, equivalencies = _generic_grid_unit(leaf)

    def compute(ops: TranslationOps, grid: Any, *, amplitude: Any, x_0: Any, alpha: Any) -> Any:
        return amplitude * (grid / x_0) ** (-alpha)

    return LeafFormula(
        ("amplitude", "x_0", "alpha"), compute, grid_unit=unit, grid_equivalencies=equivalencies
    )


def _broken_power_law(leaf: Any) -> LeafFormula:
    unit, equivalencies = _generic_grid_unit(leaf)

    def compute(
        ops: TranslationOps,
        grid: Any,
        *,
        amplitude: Any,
        x_break: Any,
        alpha_1: Any,
        alpha_2: Any,
    ) -> Any:
        ratio = grid / x_break
        below = amplitude * ratio ** (-alpha_1)
        above = amplitude * ratio ** (-alpha_2)
        return ops.where(grid < x_break, below, above)

    return LeafFormula(
        ("amplitude", "x_break", "alpha_1", "alpha_2"),
        compute,
        grid_unit=unit,
        grid_equivalencies=equivalencies,
    )


def _polynomial(leaf: Any) -> LeafFormula:
    degree = int(leaf.degree)
    names = tuple(f"c{i}" for i in range(degree + 1))
    unit, equivalencies = _generic_grid_unit(leaf)
    domain, window = leaf.domain, leaf.window
    extra: dict[str, Any] = {}
    scale, offset = 1.0, 0.0
    if domain is not None:
        domain_lo, domain_hi = float(domain[0]), float(domain[1])
        window_lo, window_hi = float(window[0]), float(window[1])
        scale = (window_hi - window_lo) / (domain_hi - domain_lo)
        offset = (window_lo * domain_hi - window_hi * domain_lo) / (domain_hi - domain_lo)
        # astropy's own default is domain == window == (-1, 1) (not None), which
        # maps identically and is not a distinguishing fact about this leaf, so
        # it is not recorded — only a domain that actually rescales is.
        if (scale, offset) != (1.0, 0.0):
            extra = {"domain": (domain_lo, domain_hi), "window": (window_lo, window_hi)}

    def compute(ops: TranslationOps, grid: Any, **coefficients: Any) -> Any:
        x = offset + scale * grid
        # Horner's method in ascending order: the same recursion astropy's own
        # PolynomialModel.horner uses, checked against it directly.
        result = coefficients[names[-1]]
        for name in reversed(names[:-1]):
            result = result * x + coefficients[name]
        return result

    return LeafFormula(
        names, compute, grid_unit=unit, grid_equivalencies=equivalencies, extra=extra
    )


def _gaussian(leaf: Any) -> LeafFormula:
    unit, equivalencies = _generic_grid_unit(leaf)

    def compute(ops: TranslationOps, grid: Any, *, amplitude: Any, mean: Any, stddev: Any) -> Any:
        return amplitude * ops.exp(-0.5 * ((grid - mean) / stddev) ** 2)

    return LeafFormula(
        ("amplitude", "mean", "stddev"), compute, grid_unit=unit, grid_equivalencies=equivalencies
    )


def _const(leaf: Any) -> LeafFormula:
    unit, equivalencies = _generic_grid_unit(leaf)

    def compute(ops: TranslationOps, grid: Any, *, amplitude: Any) -> Any:
        # Broadcasts to the grid's own shape without a dedicated op: 0.0 times
        # an array of that shape, plus the (possibly tensor) amplitude.
        return amplitude + 0.0 * grid

    return LeafFormula(("amplitude",), compute, grid_unit=unit, grid_equivalencies=equivalencies)


#: astropy class -> builder(leaf instance) -> :class:`LeafFormula`. The
#: curated table's physics, keyed the same way
#: ``ampere.backends.{torch,jax}.astropy.TRANSLATIONS`` are, so both backends
#: build theirs from this one dict (see the module docstring: "one table
#: serves both backends").
LEAF_BUILDERS: dict[type, Callable[[Any], LeafFormula]] = {
    BlackBody: _blackbody,
    PowerLaw1D: _power_law,
    BrokenPowerLaw1D: _broken_power_law,
    Polynomial1D: _polynomial,
    Gaussian1D: _gaussian,
    Const1D: _const,
}


def leaf_grid_values(formula: LeafFormula, grid: Any) -> np.ndarray:
    """*grid* (the negotiated axis, an astropy ``Quantity``) as the bare numbers *formula* wants.

    Computed once, at configuration time, and never in the evaluation loop
    (``DEVELOPMENT_PLAN.md`` §7): the backend converts the result to its own
    array type once and keeps it as a buffer.
    """
    if formula.grid_unit is None:
        return np.asarray(grid.value if isinstance(grid, u.Quantity) else grid, dtype=float)
    return np.asarray(
        grid.to_value(formula.grid_unit, equivalencies=list(formula.grid_equivalencies)),
        dtype=float,
    )


# ---------------------------------------------------------------------------
# Composition: astropy's own expression tree, walked once
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class LeafPlan:
    """One leaf of a (possibly compound) astropy model, ready to translate.

    Attributes
    ----------
    index
        Position in left-to-right leaf order — astropy's own numbering for a
        compound model's parameter suffixes.
    model
        The astropy leaf instance itself (for :func:`leaf_grid_values`'s
        input-unit probe and, indirectly, nothing else — parameter values
        never come from here, only from the composed model's own context).
    suffix
        ``""`` for a model that is not compound (astropy gives its parameters
        no suffix either), ``f"_{index}"`` otherwise.
    formula
        This leaf's translated physics.
    """

    index: int
    model: Any
    suffix: str
    formula: LeafFormula

    @property
    def context_names(self) -> tuple[str, ...]:
        """This leaf's parameters, under the names they are registered on the composed model."""
        return tuple(f"{name}{self.suffix}" for name in self.formula.param_names)


_OPERATORS: dict[str, Callable[[Any, Any], Any]] = {
    "+": operator.add,
    "-": operator.sub,
    "*": operator.mul,
    "/": operator.truediv,
}


def unsupported_operator_refusal(backend: str, model: Any, op: str) -> CapabilityError:
    """The refusal for an astropy composition operator with no elementwise flux meaning.

    ``+``, ``-``, ``*`` and ``/`` are the operators every curated leaf's own
    formula is built from; astropy's other two, ``|`` (pipe an output into the
    next model's input) and ``&`` (concatenate separate inputs), have no such
    reading as one channel's flux and are refused by name rather than given
    one.
    """
    return CapabilityError(
        f"ampere.backends.{backend}.from_astropy() cannot compose {type(model).__name__} "
        f"natively: its expression tree combines two submodels with the {op!r} operator, and "
        f"the curated table only gives a native, elementwise meaning to '+', '-', '*' and '/' "
        f"— the operators every curated leaf's own formula is built from. astropy's other "
        f"composition operators ('|' pipelining one model's output into the next model's "
        f"input, '&' concatenating separate inputs) have no such reading as one channel's "
        f"flux. Use ampere.core.from_astropy(), which evaluates your actual astropy model on "
        f"the reference path, or restructure the expression as a sum/product/difference/ratio "
        f"of the translatable pieces."
    )


def tied_parameter_refusal(backend: str, names: Sequence[str]) -> CapabilityError:
    """The refusal for a compound model carrying an astropy ``tied=`` parameter.

    Refused on **both** native backends rather than attempted on one: see the
    module docstring's "What is refused, and why".
    """
    listed = ", ".join(sorted(names))
    return CapabilityError(
        f"ampere.backends.{backend}.from_astropy() cannot give {listed} a native, "
        f"differentiable meaning. astropy's tied= is an arbitrary Python callable of the whole "
        f"model, evaluated on plain floats read off the astropy model's own Parameter objects "
        f"(ampere.core.astropy_compat.AstropyTie) — there is no gradient through it on this "
        f"backend, and no way to evaluate it at all while a value being sampled is a traced "
        f"array rather than a concrete float. Clear the tie on the astropy model, or use "
        f"ampere.core.from_astropy(), whose black-box route applies it exactly as astropy "
        f"defines it."
    )


def plan_translation(
    model: Any, backend: str
) -> tuple[tuple[LeafPlan, ...], Callable[[Sequence[Any]], Any]]:
    """The leaves of *model* (in astropy's own traversal order) and how to combine them.

    A single recursive walk of astropy's ``.left``/``.right``/``.op`` tree
    (``n_submodels`` distinguishes an internal node from a leaf, exactly as
    ``ampere.core.astropy_compat._iter_leaves`` reads it), building
    :class:`LeafPlan`\\ s and a matching combinator in one pass, so the two
    cannot number the leaves differently. A leaf's astropy class must already
    be known to be in :data:`LEAF_BUILDERS` (the backend hook's job, before
    this is called); an operator this table cannot give a native meaning to
    raises :func:`unsupported_operator_refusal` instead.

    Returns
    -------
    tuple
        ``(leaves, combine)``: the leaves in traversal order, and
        ``combine(leaf_arrays) -> array``, which applies the same tree of
        ``+``/``-``/``*``/``/`` astropy's own expression uses.
    """
    leaves: list[LeafPlan] = []

    def walk(node: Any) -> Callable[[Sequence[Any]], Any]:
        if getattr(node, "n_submodels", 1) <= 1:
            suffix = f"_{len(leaves)}" if model.n_submodels > 1 else ""
            formula = LEAF_BUILDERS[type(node)](node)
            plan = LeafPlan(len(leaves), node, suffix, formula)
            leaves.append(plan)
            position = plan.index
            return lambda values: values[position]
        if node.op not in _OPERATORS:
            raise unsupported_operator_refusal(backend, model, node.op)
        left = walk(node.left)
        right = walk(node.right)
        combine_pair = _OPERATORS[node.op]
        return lambda values: combine_pair(left(values), right(values))

    combine = walk(model)
    return tuple(leaves), combine
