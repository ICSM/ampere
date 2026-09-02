"""The parameter and prior contract (``DEVELOPMENT_PLAN.md`` §4.1).

This module is the backend-neutral vocabulary in which every ampere model,
instrument, likelihood and dataset declares *what can vary* and *what is
constant*. Nothing here knows about torch, jax, numpyro or paramax: it is
numpy/scipy/astropy.units/stdlib only, per ``architecture.md`` §3-4. Lowering
those neutral declarations onto a backend is W1.9's spec and Phase 2's code.

The narrative version of this contract, with worked examples, is
``docs/design/contracts/parameters.md``; that document's examples are executed
as doctests by ``tests/core/test_spec_doctests.py``, so it cannot drift from
this module without the suite going red.

What lives here
---------------
:class:`Parameter`
    One named quantity: prior (declared neutrally, scipy-style), optional
    fixed value, shape, unit, unconstraining bijection, tying label, plate
    membership.
:class:`Buffer`
    One named constant array — the *other* half of ``architecture.md`` §6's
    parameters-vs-buffers contract. Declared explicitly; never inferred.
:class:`ParameterSet`, :class:`BufferSet`
    Ordered, named, immutable collections. ``ParameterSet`` owns pack/unpack,
    ``lnprior``, ``prior_transform``, prior sampling, the unconstrained-space
    round trip, neutral (de)serialisation, and :meth:`ParameterSet.merge`.
:class:`Plate`
    Declarative plate-aware grouping: N members sharing hyperparameters. The
    population-model hook the plan's design horizon (d) asks to keep open.
:class:`Tie`, :class:`Binding`, :class:`ParameterMapping`
    Tying/sharing across models and datasets, and the result of merging.
:class:`Parameterised`
    The declaration mixin a model-like object uses: ``register_parameter`` /
    ``register_buffer``, plus :meth:`Parameterised.context`, the single
    namespace through which an evaluation reads both.

Three conventions worth reading before the API
----------------------------------------------
**Values are plain numbers in the declared unit.** A parameter's unit is
metadata; it is applied once, at composition time (``Parameter.to_value``
accepts a :class:`~astropy.units.Quantity` and returns a float). No
:class:`~astropy.units.Quantity` ever enters the hot loop — see the units trap
in ``DEVELOPMENT_PLAN.md`` §7.

**The flat vector is free parameters only; the value mapping is everything.**
:meth:`ParameterSet.pack` writes the free parameters into a flat array for
samplers; :meth:`ParameterSet.unpack` returns *all* parameters, fixed ones
included. This asymmetry is deliberate: it is what lets a user fix a parameter
(or promote a buffer to one) without the model's evaluation code changing at
all, which is ``architecture.md`` §6's stated test of this contract.

**"Fixed" is not a delta prior.** A fixed parameter has a value and no prior.
It contributes nothing to ``lnprior``, occupies no sampler dimension, and is
still a named part of the model's vocabulary. A delta-function prior would be
a (degenerate, badly behaved) *free* parameter; the two are not
interchangeable and this module refuses to conflate them.
"""

from __future__ import annotations

import dataclasses
import math
import types
from collections.abc import Iterable, Iterator, Mapping, Sequence
from typing import Any, Protocol, runtime_checkable

import astropy.units as u
import numpy as np
import scipy.stats as _stats

from .exceptions import OptionalDependencyError, ParameterError, TyingError

__all__ = [
    "SEPARATOR",
    "Bijection",
    "Binding",
    "Buffer",
    "BufferSet",
    "HierarchicalPrior",
    "Identity",
    "Log",
    "Logit",
    "Parameter",
    "ParameterMapping",
    "ParameterSet",
    "Parameterised",
    "Plate",
    "PlateBinding",
    "Prior",
    "PriorSpec",
    "Tie",
    "default_bijection_for",
    "describe_prior",
    "log_density",
    "prior_from_spec",
]

#: Separator between a component label and a local parameter name in a merged
#: :class:`ParameterSet` (``"spectrum.temperature"``), and between a plate name
#: and its members (``"objects.theta"``). Local names may never contain it.
SEPARATOR = "."

_EMPTY_STR_MAP: Mapping[str, str] = types.MappingProxyType({})
_EMPTY_NUM_MAP: Mapping[str, float] = types.MappingProxyType({})

ArrayLike = Any
Value = Any


# ---------------------------------------------------------------------------
# Names
# ---------------------------------------------------------------------------


def _check_name(name: object, kind: str = "parameter") -> str:
    """Validate a (possibly qualified) declaration name.

    Every dot-separated segment must be a Python identifier, because unqualified
    parameter names are handed to models as keyword arguments and qualified ones
    are built by joining identifiers with :data:`SEPARATOR`.
    """
    if not isinstance(name, str) or not name:
        raise ParameterError(f"a {kind} name must be a non-empty string, got {name!r}")
    if not all(segment.isidentifier() for segment in name.split(SEPARATOR)):
        raise ParameterError(
            f"{kind} name {name!r} is not usable: every {SEPARATOR!r}-separated segment must "
            f"be a valid Python identifier, because parameter values are passed to models as "
            f"keyword arguments and qualified names are built by joining identifiers."
        )
    return name


def _check_local_name(name: object, kind: str = "parameter") -> str:
    """Validate a name that must *not* already be qualified."""
    checked = _check_name(name, kind)
    if SEPARATOR in checked:
        raise ParameterError(
            f"{kind} name {checked!r} must not contain {SEPARATOR!r}: qualified names are "
            f"produced by ParameterSet.merge() and Plate expansion, never declared directly."
        )
    return checked


# ---------------------------------------------------------------------------
# Priors: neutral declaration and neutral description
# ---------------------------------------------------------------------------


@runtime_checkable
class Prior(Protocol):
    """Structural type of a frozen, one-dimensional prior distribution.

    A frozen ``scipy.stats`` distribution (``scipy.stats.norm(0, 1)``,
    ``scipy.stats.loguniform(1e-3, 1e3)``, ``scipy.stats.poisson(3.0)``) is the
    **canonical** declaration and the only form ampere can lower to
    ``torch.distributions`` / numpyro equivalents (W1.9's mapping table works
    from :class:`PriorSpec`, which is read off a frozen scipy distribution).

    Anything else that satisfies this protocol will *evaluate* correctly on the
    reference path but is opaque to lowering and to serialisation:
    :func:`describe_prior` raises for it, and so therefore does
    :meth:`ParameterSet.to_spec`. That is deliberate — an un-lowerable prior
    should fail at declaration/serialisation time with a message naming the
    parameter, not deep inside a backend.

    Notes
    -----
    ``logpdf`` is not part of the protocol because discrete distributions
    expose ``logpmf`` instead; use :func:`log_density`, which accepts either.
    """

    def ppf(self, q: Any) -> Any:
        """Percent-point function (inverse CDF)."""
        ...

    def support(self) -> tuple[Any, Any]:
        """``(lower, upper)`` bounds of the distribution's support."""
        ...


def log_density(prior: object, x: ArrayLike) -> np.ndarray:
    """Log density (continuous) or log mass (discrete) of *prior* at *x*.

    Exists so that callers never have to care whether a prior is continuous or
    discrete; ``scipy`` names those ``logpdf`` and ``logpmf`` respectively.
    """
    fn = getattr(prior, "logpdf", None) or getattr(prior, "logpmf", None)
    if fn is None:
        raise ParameterError(
            f"prior {prior!r} exposes neither 'logpdf' nor 'logpmf'; a prior must be able "
            f"to evaluate its own log density (a frozen scipy.stats distribution can)."
        )
    return np.asarray(fn(x), dtype=float)


def _numeric(value: object, what: str) -> float:
    """Coerce a prior argument to a plain float, or explain why it cannot be."""
    try:
        array = np.asarray(value)
        if array.ndim != 0:
            raise TypeError(what)
        return float(array)
    except (TypeError, ValueError) as exc:
        raise ParameterError(
            f"{what} must be a numeric scalar to be described neutrally and serialised, "
            f"got {value!r}. Array-valued distribution parameters are not part of the "
            f"v1.3 contract; use a Plate or a HierarchicalPrior instead."
        ) from exc


@dataclasses.dataclass(frozen=True)
class PriorSpec:
    """A neutral, serialisable description of a frozen scipy distribution.

    This is the handle W1.9's lowering table works from: a family name plus the
    numbers it was frozen with, in scipy's parametrisation. Translating scipy's
    ``loc``/``scale`` convention into torch's or numpyro's is W1.9's job; this
    contract's job is to make sure the information needed to do so survives
    declaration, merging and round-tripping.

    Parameters
    ----------
    family
        The ``scipy.stats`` distribution name, e.g. ``"norm"``, ``"uniform"``,
        ``"loguniform"``, ``"poisson"``.
    args, kwds
        Positional and keyword arguments the distribution was frozen with.
        :func:`describe_prior` always emits the canonical keyword-only form
        (``args`` empty); ``args`` remains accepted for hand-written specs.
    discrete
        Whether the family is discrete (``rv_discrete``).

    Examples
    --------
    >>> import scipy.stats
    >>> spec = describe_prior(scipy.stats.norm(loc=1.0, scale=2.0))
    >>> spec.family, dict(spec.kwds)
    ('norm', {'loc': 1.0, 'scale': 2.0})
    >>> prior_from_spec(spec).mean()
    np.float64(1.0)
    """

    family: str
    args: tuple[float, ...] = ()
    kwds: Mapping[str, float] = _EMPTY_NUM_MAP
    discrete: bool = False

    def __post_init__(self) -> None:
        if not isinstance(self.family, str) or not self.family:
            raise ParameterError(f"a prior family must be a non-empty string, got {self.family!r}")
        object.__setattr__(self, "args", tuple(self.args))
        object.__setattr__(self, "kwds", types.MappingProxyType(dict(self.kwds)))

    @classmethod
    def of(cls, family: str, *args: float, discrete: bool = False, **kwds: float) -> PriorSpec:
        """Construct a spec the way one would call the distribution itself.

        >>> PriorSpec.of("uniform", loc=100.0, scale=9900.0).family
        'uniform'
        """
        return cls(family=family, args=args, kwds=kwds, discrete=discrete)

    def to_dict(self) -> dict[str, Any]:
        """A JSON-compatible dictionary (the serialisation form)."""
        return {
            "family": self.family,
            "args": list(self.args),
            "kwds": dict(self.kwds),
            "discrete": self.discrete,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> PriorSpec:
        """Inverse of :meth:`to_dict`."""
        return cls(
            family=data["family"],
            args=tuple(data.get("args", ())),
            kwds=dict(data.get("kwds", {})),
            discrete=bool(data.get("discrete", False)),
        )


def _argument_names(dist: object) -> list[str]:
    """The ordered freezing-argument names of a scipy distribution family.

    Shape parameters first (``dist.shapes``, e.g. ``"a, b"`` for
    ``loguniform``), then ``loc``, then — for continuous families — ``scale``;
    the order scipy itself assigns positional arguments when freezing.
    """
    shapes = getattr(dist, "shapes", None)
    names = [] if not shapes else [part.strip() for part in shapes.split(",")]
    names.append("loc")
    if isinstance(dist, _stats.rv_continuous):
        names.append("scale")
    return names


def describe_prior(prior: object) -> PriorSpec:
    """Describe a frozen ``scipy.stats`` distribution as a :class:`PriorSpec`.

    The description is **canonical**: scipy accepts the same freezing either
    positionally or by keyword, so positional arguments are mapped onto their
    names here (via the family's shape names plus ``loc``/``scale``). Two
    declarations of one distribution therefore describe — and compare, tie
    and serialise — identically, and W1.9's lowering table has a single form
    to translate.

    Raises
    ------
    ParameterError
        If *prior* is not a frozen ``scipy.stats`` distribution, i.e. if it
        cannot be lowered onto a backend or written to a run's provenance
        record. The message says so explicitly rather than letting the failure
        surface later, somewhere less informative.
    """
    dist = getattr(prior, "dist", None)
    family = getattr(dist, "name", None)
    if not isinstance(family, str):
        raise ParameterError(
            f"prior {prior!r} is not a frozen scipy.stats distribution, so it cannot be "
            f"described neutrally (and therefore cannot be lowered to torch/numpyro or "
            f"recorded in a run's provenance). Declare priors as e.g. "
            f"scipy.stats.norm(0, 1); a custom prior object may still be evaluated on the "
            f"reference path, but it must not be serialised or lowered."
        )
    names = _argument_names(dist)
    args = tuple(getattr(prior, "args", ()))
    if len(args) > len(names):
        raise ParameterError(
            f"prior family {family!r} was frozen with {len(args)} positional argument(s) "
            f"but only takes {names}; cannot describe it neutrally."
        )
    kwds = {
        name: _numeric(value, f"argument {name!r} of prior family {family!r}")
        for name, value in zip(names, args, strict=False)
    }
    kwds |= {
        str(k): _numeric(v, f"keyword {k!r} of prior family {family!r}")
        for k, v in dict(getattr(prior, "kwds", {})).items()
    }
    return PriorSpec(family=family, kwds=kwds, discrete=isinstance(dist, _stats.rv_discrete))


def _distribution_factory(family: str) -> Any:
    factory = getattr(_stats, family, None)
    if not isinstance(factory, (_stats.rv_continuous, _stats.rv_discrete)):
        raise ParameterError(
            f"{family!r} is not a scipy.stats distribution; known prior families are the "
            f"distributions exported by scipy.stats."
        )
    return factory


def prior_from_spec(spec: PriorSpec) -> Prior:
    """Rebuild a frozen scipy distribution from a :class:`PriorSpec`.

    Together with :func:`describe_prior` this is the prior round trip that
    :meth:`ParameterSet.to_spec` / :meth:`ParameterSet.from_spec` are built on.
    """
    return _distribution_factory(spec.family)(*spec.args, **dict(spec.kwds))


def _priors_equal(left: object, right: object) -> bool:
    """Compare two priors by neutral description, falling back to identity."""
    if left is right:
        return True
    if isinstance(left, HierarchicalPrior) or isinstance(right, HierarchicalPrior):
        return left == right
    if left is None or right is None:
        return False
    try:
        return describe_prior(left) == describe_prior(right)
    except ParameterError:
        return False


# ---------------------------------------------------------------------------
# Hierarchical priors
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class HierarchicalPrior:
    """A prior whose distribution parameters are themselves ampere parameters.

    This is what makes "each object has its own ``theta_i ~ Normal(mu, sigma)``"
    a *declaration* rather than a callback: the family is a neutral name (it
    goes through exactly the same W1.9 lowering table as any other prior) and
    the hyperparameters are recorded as **references by name**, which is what a
    numpyro/pyro model does anyway (``dist.Normal(mu, sigma)`` where ``mu`` and
    ``sigma`` are earlier sample sites).

    Parameters
    ----------
    family
        ``scipy.stats`` distribution name, as for :class:`PriorSpec`.
    hyperparameters
        Mapping from the distribution's keyword argument to the *name of the
        parameter* supplying it, e.g. ``{"loc": "mu", "scale": "sigma"}``.
        Names are resolved within the enclosing :class:`ParameterSet`; a
        :class:`Plate` rewrites references to its own hyperparameters into
        qualified names when it expands.
    args, kwds
        Any remaining, constant arguments of the family.

    Examples
    --------
    >>> hp = HierarchicalPrior("norm", {"loc": "mu", "scale": "sigma"})
    >>> hp.references
    ('mu', 'sigma')
    >>> float(hp.bind({"mu": 2.0, "sigma": 0.5}).mean())
    2.0
    """

    family: str
    hyperparameters: Mapping[str, str] = _EMPTY_STR_MAP
    args: tuple[float, ...] = ()
    kwds: Mapping[str, float] = _EMPTY_NUM_MAP

    def __post_init__(self) -> None:
        if not isinstance(self.family, str) or not self.family:
            raise ParameterError(f"a prior family must be a non-empty string, got {self.family!r}")
        hyper = dict(self.hyperparameters)
        if not hyper:
            raise ParameterError(
                f"HierarchicalPrior({self.family!r}) declares no hyperparameters; a prior with "
                f"no references to other parameters is an ordinary prior — declare it as a "
                f"frozen scipy.stats distribution instead."
            )
        for keyword, reference in hyper.items():
            if not isinstance(keyword, str) or not keyword.isidentifier():
                raise ParameterError(
                    f"hyperparameter key {keyword!r} must be a keyword argument name of "
                    f"scipy.stats.{self.family}"
                )
            _check_name(reference, "hyperparameter reference")
        overlap = set(hyper) & set(self.kwds)
        if overlap:
            raise ParameterError(
                f"HierarchicalPrior({self.family!r}) supplies {sorted(overlap)} both as a "
                f"constant and as a hyperparameter reference; pick one."
            )
        object.__setattr__(self, "hyperparameters", types.MappingProxyType(hyper))
        object.__setattr__(self, "args", tuple(self.args))
        object.__setattr__(self, "kwds", types.MappingProxyType(dict(self.kwds)))

    @property
    def references(self) -> tuple[str, ...]:
        """Names of the parameters this prior depends on, in declaration order."""
        return tuple(self.hyperparameters.values())

    def bind(self, values: Mapping[str, Value]) -> Prior:
        """Freeze the distribution using hyperparameter values from *values*."""
        try:
            resolved = {kw: values[ref] for kw, ref in self.hyperparameters.items()}
        except KeyError as exc:
            raise ParameterError(
                f"hierarchical prior {self.family!r} references parameter {exc.args[0]!r}, "
                f"which is not available; known names are {sorted(values)}."
            ) from exc
        return _distribution_factory(self.family)(*self.args, **dict(self.kwds), **resolved)

    def rename_references(self, rename: Mapping[str, str]) -> HierarchicalPrior:
        """Return a copy with references remapped (used by plate expansion and merging)."""
        return dataclasses.replace(
            self,
            hyperparameters={kw: rename.get(ref, ref) for kw, ref in self.hyperparameters.items()},
        )

    def to_dict(self) -> dict[str, Any]:
        """A JSON-compatible dictionary (the serialisation form)."""
        return {
            "family": self.family,
            "hyperparameters": dict(self.hyperparameters),
            "args": list(self.args),
            "kwds": dict(self.kwds),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> HierarchicalPrior:
        """Inverse of :meth:`to_dict`."""
        return cls(
            family=data["family"],
            hyperparameters=dict(data["hyperparameters"]),
            args=tuple(data.get("args", ())),
            kwds=dict(data.get("kwds", {})),
        )


AnyPrior = Any


# ---------------------------------------------------------------------------
# Bijections to and from unconstrained space
# ---------------------------------------------------------------------------


@runtime_checkable
class Bijection(Protocol):
    """A smooth, invertible map between constrained and unconstrained space.

    Gradient-based samplers (NUTS/HMC) and gradient-based optimisers need the
    parameter vector to live in :math:`\\mathbb{R}^n`. A :class:`Bijection`
    says how a constrained parameter gets there and back.

    The method names are deliberately *not* ``forward``/``inverse``. Both torch
    and numpyro use ``forward`` for unconstrained-to-constrained while a reader
    coming from the statistics literature usually expects the opposite, and a
    reversed log-determinant is a silent, hard-to-find bug. ``constrain`` and
    ``unconstrain`` cannot be read backwards.

    Notes
    -----
    :meth:`log_abs_det_jacobian` is evaluated at the **unconstrained** point,
    matching ``constrain``'s argument; this is the term added to the log density
    when changing variables, i.e.

    .. math:: \\log p(y) = \\log p(x = c(y)) + \\log|\\mathrm{d}c/\\mathrm{d}y|.
    """

    def constrain(self, y: ArrayLike) -> np.ndarray:
        """Map an unconstrained value into the parameter's support."""
        ...

    def unconstrain(self, x: ArrayLike) -> np.ndarray:
        """Map a value in the parameter's support to unconstrained space."""
        ...

    def log_abs_det_jacobian(self, y: ArrayLike) -> np.ndarray:
        """Log absolute determinant of ``d constrain / d y``, evaluated at *y*."""
        ...


@dataclasses.dataclass(frozen=True)
class Identity:
    """The trivial bijection, for parameters already supported on the real line."""

    def constrain(self, y: ArrayLike) -> np.ndarray:
        return np.asarray(y, dtype=float)

    def unconstrain(self, x: ArrayLike) -> np.ndarray:
        return np.asarray(x, dtype=float)

    def log_abs_det_jacobian(self, y: ArrayLike) -> np.ndarray:
        return np.zeros_like(np.asarray(y, dtype=float))


@dataclasses.dataclass(frozen=True)
class Log:
    """Bijection for a half-line support :math:`[l, \\infty)`: ``x = l + exp(y)``.

    Examples
    --------
    >>> b = Log()
    >>> float(b.unconstrain(np.e))
    1.0
    >>> float(b.log_abs_det_jacobian(1.0))
    1.0
    """

    lower: float = 0.0

    def constrain(self, y: ArrayLike) -> np.ndarray:
        return self.lower + np.exp(np.asarray(y, dtype=float))

    def unconstrain(self, x: ArrayLike) -> np.ndarray:
        return np.log(np.asarray(x, dtype=float) - self.lower)

    def log_abs_det_jacobian(self, y: ArrayLike) -> np.ndarray:
        return np.asarray(y, dtype=float)


@dataclasses.dataclass(frozen=True)
class Logit:
    """Bijection for a bounded support :math:`[l, h]` via the logistic map.

    Examples
    --------
    >>> b = Logit(100.0, 10000.0)
    >>> float(b.constrain(0.0))
    5050.0
    >>> float(b.unconstrain(b.constrain(1.25)))
    1.25
    """

    lower: float
    upper: float

    def __post_init__(self) -> None:
        if not self.upper > self.lower:
            raise ParameterError(
                f"Logit bounds must satisfy upper > lower, got {self.lower} and {self.upper}"
            )

    def constrain(self, y: ArrayLike) -> np.ndarray:
        yy = np.asarray(y, dtype=float)
        return self.lower + (self.upper - self.lower) / (1.0 + np.exp(-yy))

    def unconstrain(self, x: ArrayLike) -> np.ndarray:
        xx = np.asarray(x, dtype=float)
        return np.log(xx - self.lower) - np.log(self.upper - xx)

    def log_abs_det_jacobian(self, y: ArrayLike) -> np.ndarray:
        yy = np.asarray(y, dtype=float)
        width = math.log(self.upper - self.lower)
        return width - np.logaddexp(0.0, -yy) - np.logaddexp(0.0, yy)


def _bijection_for_support(lower: float, upper: float, what: str) -> Bijection:
    finite_lower, finite_upper = math.isfinite(lower), math.isfinite(upper)
    if not finite_lower and not finite_upper:
        return Identity()
    if finite_lower and not finite_upper:
        return Log(lower=lower)
    if finite_lower and finite_upper:
        return Logit(lower=lower, upper=upper)
    raise ParameterError(
        f"cannot infer an unconstraining bijection for {what}: its support "
        f"({lower}, {upper}) is bounded above but not below, which ampere v1.3 has no "
        f"built-in bijection for. Declare one explicitly via Parameter(bijection=...); "
        f"any object satisfying the Bijection protocol will do."
    )


def default_bijection_for(prior: AnyPrior, what: str = "this prior") -> Bijection:
    """Infer the unconstraining bijection implied by a prior's support.

    Real-line support gives :class:`Identity`, a half-line gives :class:`Log`,
    a bounded interval gives :class:`Logit` — the same rule numpyro's
    ``biject_to`` applies, stated here so the reference path and the eventual
    torch/jax paths agree by construction rather than by coincidence.

    For a :class:`HierarchicalPrior` the support may depend on values that are
    not known until sampling time. Inference is therefore attempted only when
    it is provably safe (a location-scale family with no shape arguments, whose
    bound-determining hyperparameters are constants); otherwise this raises and
    asks for an explicit declaration, rather than guessing.

    Examples
    --------
    >>> import scipy.stats
    >>> default_bijection_for(scipy.stats.norm(0, 1))
    Identity()
    >>> default_bijection_for(scipy.stats.uniform(loc=100.0, scale=9900.0))
    Logit(lower=100.0, upper=10000.0)
    >>> default_bijection_for(scipy.stats.halfnorm(0.0, 1.0))
    Log(lower=0.0)
    """
    if isinstance(prior, HierarchicalPrior):
        return _default_bijection_for_hierarchical(prior, what)
    support = getattr(prior, "support", None)
    if support is None:
        raise ParameterError(
            f"cannot infer an unconstraining bijection for {what}: the prior does not expose "
            f"support(). Declare one explicitly via Parameter(bijection=...)."
        )
    lower, upper = support()
    return _bijection_for_support(float(lower), float(upper), what)


def _default_bijection_for_hierarchical(prior: HierarchicalPrior, what: str) -> Bijection:
    dist = _distribution_factory(prior.family)
    referenced = set(prior.hyperparameters)
    if getattr(dist, "numargs", 0):
        raise ParameterError(
            f"cannot infer an unconstraining bijection for {what}: the hierarchical family "
            f"{prior.family!r} takes shape arguments, so its support is not determined by "
            f"loc/scale alone. Declare one explicitly via Parameter(bijection=...)."
        )
    standard_lower, standard_upper = float(dist.a), float(dist.b)
    if not math.isfinite(standard_lower) and not math.isfinite(standard_upper):
        # loc/scale cannot make an unbounded support bounded; safe regardless of values.
        return Identity()
    constants = dict(prior.kwds)
    positional = ("loc", "scale")
    for index, keyword in enumerate(positional):
        if index < len(prior.args):
            constants[keyword] = prior.args[index]
    loc = constants.get("loc", 0.0)
    scale = constants.get("scale", 1.0)
    if "loc" in referenced or ("scale" in referenced and standard_upper != standard_lower):
        needed = sorted(referenced & {"loc", "scale"})
        if standard_lower == 0.0 and not math.isfinite(standard_upper) and "loc" not in referenced:
            # A lower bound at loc is unaffected by scale > 0.
            return Log(lower=float(loc))
        raise ParameterError(
            f"cannot infer an unconstraining bijection for {what}: family {prior.family!r} has "
            f"a bounded support whose position depends on hyperparameter(s) {needed}, so the "
            f"bijection would change from draw to draw. Declare one explicitly via "
            f"Parameter(bijection=...)."
        )
    return _bijection_for_support(loc + scale * standard_lower, loc + scale * standard_upper, what)


# ---------------------------------------------------------------------------
# Units and values
# ---------------------------------------------------------------------------


def _check_unit(unit: object) -> u.UnitBase | None:
    if unit is None:
        return None
    if not isinstance(unit, u.UnitBase):
        raise ParameterError(
            f"unit must be an astropy unit (e.g. astropy.units.micron) or None, got {unit!r}. "
            f"Strings are not accepted: use astropy.units.Unit('micron') at the call site so "
            f"that a typo fails there rather than inside ampere."
        )
    return unit


def _readonly(array: np.ndarray) -> np.ndarray:
    array.setflags(write=False)
    return array


def _normalise_shape(shape: object) -> tuple[int, ...]:
    if isinstance(shape, (int, np.integer)):
        dims = (int(shape),)
    elif isinstance(shape, Iterable):
        try:
            dims = tuple(int(d) for d in shape)
        except (TypeError, ValueError) as exc:
            raise ParameterError(f"shape must be an int or a tuple of ints, got {shape!r}") from exc
    else:
        raise ParameterError(f"shape must be an int or a tuple of ints, got {shape!r}")
    if any(d < 1 for d in dims):
        raise ParameterError(f"every shape dimension must be >= 1, got {dims!r}")
    return dims


# ---------------------------------------------------------------------------
# Parameter
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True, eq=False)
class Parameter:
    """One named quantity a model, instrument, noise model or dataset owns.

    A ``Parameter`` is in exactly one of three states:

    ``free``
        Has a ``prior``, is not ``fixed``. Occupies :attr:`size` dimensions of
        the sampler's flat vector and contributes to ``lnprior``.
    ``fixed``
        Has a ``value`` and no ``prior``. Occupies no sampler dimension and
        contributes nothing to ``lnprior``, but is still handed to the model on
        every evaluation. **Not** the same thing as a delta-function prior.
    ``deferred``
        Has ``shared_as`` set and no prior of its own: the prior is contributed
        by another site in the same tie group, and the parameter becomes usable
        once :meth:`ParameterSet.merge` has resolved the group. A set holding
        one is :attr:`~ParameterSet.is_resolved` ``False`` and refuses to
        evaluate priors.

    Parameters
    ----------
    name
        A Python identifier. Parameter values reach models as keyword
        arguments, so a name that is not an identifier is rejected at
        declaration.
    prior
        A frozen ``scipy.stats`` distribution (canonical; see :class:`Prior`),
        or a :class:`HierarchicalPrior`, or ``None`` for a fixed or deferred
        parameter.
    value
        Initial value for a free parameter, or *the* value for a fixed one. May
        be given as a :class:`~astropy.units.Quantity`, which is converted to
        the declared unit once, here, and stored as a plain float or array.
    fixed
        Whether the parameter is held at ``value`` for this fit.
    shape
        ``()`` for a scalar (the default); ``(n,)`` or richer for an
        array-valued parameter such as a per-channel calibration offset. The
        prior applies **independently and identically to each element**; a
        genuinely multivariate prior is out of scope for v1.3.
    unit
        An ``astropy.units`` unit, or ``None`` for a dimensionless quantity.
        The prior is declared in this unit, numerically.
    bijection
        Map to unconstrained space. ``None`` means "infer from the prior's
        support" (see :func:`default_bijection_for`).
    shared_as
        Tie label. Two parameters — in the same set or in different ones —
        carrying the same label are one free parameter with two binding sites
        after :meth:`ParameterSet.merge`.
    plate
        Name of the :class:`Plate` this parameter belongs to, if any. Set by
        :meth:`Plate.expand`; not normally declared by hand.
    description
        Free text for documentation and for results provenance.

    Examples
    --------
    >>> import scipy.stats, astropy.units as u
    >>> t = Parameter("temperature", scipy.stats.uniform(100.0, 9900.0), unit=u.K)
    >>> t.is_free, t.size, t.unit
    (True, 1, Unit("K"))
    >>> fixed = t.fix(2500.0 * u.K)
    >>> fixed.is_fixed, fixed.value, fixed.prior is None
    (True, 2500.0, True)
    """

    name: str
    prior: AnyPrior = None
    _: dataclasses.KW_ONLY
    value: Value = None
    fixed: bool = False
    shape: tuple[int, ...] = ()
    unit: u.UnitBase | None = None
    bijection: Bijection | None = None
    shared_as: str | None = None
    plate: str | None = None
    description: str = ""

    def __post_init__(self) -> None:
        object.__setattr__(self, "name", _check_name(self.name))
        object.__setattr__(self, "unit", _check_unit(self.unit))
        object.__setattr__(self, "shape", _normalise_shape(self.shape))
        if self.shared_as is not None:
            object.__setattr__(self, "shared_as", _check_local_name(self.shared_as, "tie label"))
        if self.plate is not None:
            object.__setattr__(self, "plate", _check_local_name(self.plate, "plate"))
        object.__setattr__(self, "fixed", bool(self.fixed))
        self._normalise_value()
        self._validate_state()
        if self.bijection is not None and not isinstance(self.bijection, Bijection):
            raise ParameterError(
                f"parameter {self.name!r}: bijection {self.bijection!r} does not satisfy the "
                f"Bijection protocol (constrain/unconstrain/log_abs_det_jacobian)."
            )

    def _normalise_value(self) -> None:
        if self.value is None:
            return
        raw = self.to_value(self.value)
        array = np.asarray(raw, dtype=float)
        if not self.shape and array.ndim:
            object.__setattr__(self, "shape", _normalise_shape(array.shape))
        if self.shape:
            try:
                array = np.broadcast_to(array, self.shape).astype(float, copy=True)
            except ValueError as exc:
                raise ParameterError(
                    f"parameter {self.name!r}: value of shape {array.shape} is not "
                    f"broadcastable to the declared shape {self.shape}"
                ) from exc
            object.__setattr__(self, "value", _readonly(array))
        else:
            object.__setattr__(self, "value", float(array))

    def _validate_state(self) -> None:
        if self.fixed:
            if self.prior is not None:
                raise ParameterError(
                    f"parameter {self.name!r} is declared fixed *and* given a prior. A fixed "
                    f"parameter has a value and no prior; if you meant a parameter that varies, "
                    f"drop fixed=True. A delta-function prior is not the same thing and is not "
                    f"supported (see docs/design/contracts/parameters.md)."
                )
            if self.value is None:
                raise ParameterError(
                    f"parameter {self.name!r} is declared fixed but has no value; a fixed "
                    f"parameter is defined by its value."
                )
        elif self.prior is None and self.shared_as is None:
            raise ParameterError(
                f"parameter {self.name!r} has neither a prior nor fixed=True nor a shared_as "
                f"tie label, so nothing determines it. Give it a prior, fix it at a value, or "
                f"tie it to a parameter that has one."
            )

    # -- state -------------------------------------------------------------

    @property
    def is_fixed(self) -> bool:
        """Whether this parameter is held at :attr:`value`."""
        return self.fixed

    @property
    def is_free(self) -> bool:
        """Whether this parameter varies *and* knows its own prior."""
        return not self.fixed and self.prior is not None

    @property
    def is_deferred(self) -> bool:
        """Whether this parameter is waiting for a tie group to supply its prior."""
        return not self.fixed and self.prior is None

    @property
    def is_hierarchical(self) -> bool:
        """Whether this parameter's prior references other parameters."""
        return isinstance(self.prior, HierarchicalPrior)

    @property
    def size(self) -> int:
        """Number of scalar entries this parameter occupies in a flat vector."""
        return math.prod(self.shape) if self.shape else 1

    @property
    def references(self) -> tuple[str, ...]:
        """Parameter names this one's prior depends on (empty unless hierarchical)."""
        return self.prior.references if isinstance(self.prior, HierarchicalPrior) else ()

    # -- units -------------------------------------------------------------

    def to_value(self, value: Value) -> Value:
        """Strip units from *value*, converting to the declared unit.

        This is the "convert once, at composition time" step; nothing
        downstream ever sees a :class:`~astropy.units.Quantity`.
        """
        if isinstance(value, u.Quantity):
            target = self.unit if self.unit is not None else u.dimensionless_unscaled
            try:
                return value.to_value(target)
            except u.UnitConversionError as exc:
                raise ParameterError(
                    f"parameter {self.name!r}: cannot convert {value!r} to the declared unit "
                    f"{target}"
                ) from exc
        return value

    def quantity(self, value: Value = None) -> u.Quantity:
        """Re-attach the declared unit to *value* (default: :attr:`value`)."""
        raw = self.value if value is None else value
        if raw is None:
            raise ParameterError(f"parameter {self.name!r} has no value to attach a unit to")
        return u.Quantity(raw, self.unit if self.unit is not None else u.dimensionless_unscaled)

    # -- derived declarations ---------------------------------------------

    def unconstraining_bijection(self) -> Bijection:
        """The declared bijection, or the one implied by the prior's support."""
        if self.bijection is not None:
            return self.bijection
        if self.prior is None:
            return Identity()
        return default_bijection_for(self.prior, f"parameter {self.name!r}")

    def fix(self, value: Value = None) -> Parameter:
        """Return a fixed copy of this parameter (configuration, not code, changes)."""
        target = self.value if value is None else value
        if target is None:
            raise ParameterError(
                f"parameter {self.name!r} cannot be fixed without a value: it has no declared "
                f"value to fall back on."
            )
        return dataclasses.replace(self, prior=None, value=target, fixed=True)

    def release(self, prior: AnyPrior, *, bijection: Bijection | None = None) -> Parameter:
        """Return a free copy of this parameter with *prior* (the inverse of :meth:`fix`)."""
        return dataclasses.replace(
            self,
            prior=prior,
            fixed=False,
            bijection=self.bijection if bijection is None else bijection,
        )

    def rename(self, name: str) -> Parameter:
        """Return a copy under a new (possibly qualified) name."""
        return dataclasses.replace(self, name=name)

    # -- equality ----------------------------------------------------------

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Parameter):
            return NotImplemented
        if not _values_equal(self.value, other.value):
            return False
        if (self.prior is None) != (other.prior is None):
            return False
        if self.prior is not None and not _priors_equal(self.prior, other.prior):
            return False
        return (
            self.name == other.name
            and self.fixed == other.fixed
            and self.shape == other.shape
            and self.unit == other.unit
            and self.bijection == other.bijection
            and self.shared_as == other.shared_as
            and self.plate == other.plate
            and self.description == other.description
        )

    def __hash__(self) -> int:
        # Consistent with __eq__ (equal parameters agree on all of these), while
        # staying hashable in the presence of an ndarray value.
        return hash((self.name, self.fixed, self.shape, str(self.unit), self.shared_as, self.plate))

    def __repr__(self) -> str:
        bits = [repr(self.name)]
        if self.is_fixed:
            bits.append(f"fixed={self.value!r}")
        elif self.is_deferred:
            bits.append(f"deferred, shared_as={self.shared_as!r}")
        else:
            bits.append(f"prior={_prior_repr(self.prior)}")
            if self.shared_as is not None:
                bits.append(f"shared_as={self.shared_as!r}")
        if self.shape:
            bits.append(f"shape={self.shape}")
        if self.unit is not None:
            bits.append(f"unit={self.unit}")
        if self.plate is not None:
            bits.append(f"plate={self.plate!r}")
        return f"Parameter({', '.join(bits)})"


def _prior_repr(prior: AnyPrior) -> str:
    if isinstance(prior, HierarchicalPrior):
        return f"{prior.family}({', '.join(prior.hyperparameters.values())})"
    try:
        spec = describe_prior(prior)
    except ParameterError:
        return repr(prior)
    inner = ", ".join(
        [*(repr(a) for a in spec.args), *(f"{k}={v!r}" for k, v in spec.kwds.items())]
    )
    return f"{spec.family}({inner})"


def _values_equal(left: Value, right: Value) -> bool:
    if left is None or right is None:
        return left is None and right is None
    return bool(np.array_equal(np.asarray(left), np.asarray(right)))


# ---------------------------------------------------------------------------
# Buffers
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True, eq=False)
class Buffer:
    """A named constant array a model needs but would never put a prior on.

    ``architecture.md`` §6 gives the distinguishing question: *would you ever
    put a prior on it, sample it, or take a gradient with respect to it?* If
    no, it is a buffer — wavelength grids, opacity tables, filter curves,
    response matrices. Buffers are declared **explicitly**
    (:meth:`Parameterised.register_buffer`), never inferred from "array
    attributes that are not parameters", because a mistyped or forgotten
    parameter registration silently becoming an untracked buffer is exactly
    the silent-drift bug class these contracts exist to end.

    Buffers move with device and dtype changes (``architecture.md`` §5) and are
    traced through computations, but they are never sampled, optimised or
    differentiated. On the torch side they lower to ``register_buffer``; on the
    jax side to array fields excluded from the trainable partition — **never**
    to equinox static fields, which hash array contents into the JIT cache key
    (``DEVELOPMENT_PLAN.md`` §7). Those lowerings are W1.9's table.

    The stored array is made read-only: a buffer mutated in place after
    registration would invalidate exactly the caches its constancy justifies.

    Examples
    --------
    >>> import numpy as np, astropy.units as u
    >>> wl = Buffer("wavelength", np.geomspace(1.0, 100.0, 5), unit=u.micron)
    >>> wl.shape, wl.size
    ((5,), 5)
    >>> wl.value[0] = 2.0
    Traceback (most recent call last):
        ...
    ValueError: assignment destination is read-only
    """

    name: str
    value: ArrayLike
    _: dataclasses.KW_ONLY
    unit: u.UnitBase | None = None
    description: str = ""

    def __post_init__(self) -> None:
        object.__setattr__(self, "name", _check_name(self.name, "buffer"))
        object.__setattr__(self, "unit", _check_unit(self.unit))
        raw = self.value
        if isinstance(raw, u.Quantity):
            target = self.unit if self.unit is not None else raw.unit
            object.__setattr__(self, "unit", _check_unit(target))
            raw = raw.to_value(target)
        object.__setattr__(self, "value", _readonly(np.array(raw, copy=True)))

    @property
    def array(self) -> np.ndarray:
        """The stored array (read-only). Always an :class:`numpy.ndarray`."""
        return np.asarray(self.value)

    @property
    def shape(self) -> tuple[int, ...]:
        """Shape of the stored array (``()`` for a scalar constant)."""
        return self.array.shape

    @property
    def size(self) -> int:
        """Number of scalar entries in the stored array."""
        return int(self.array.size)

    @property
    def dtype(self) -> np.dtype[Any]:
        """dtype of the stored array; buffers keep the dtype they were given."""
        return self.array.dtype

    def quantity(self) -> u.Quantity:
        """The buffer with its unit re-attached."""
        return u.Quantity(
            self.value, self.unit if self.unit is not None else u.dimensionless_unscaled
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Buffer):
            return NotImplemented
        return (
            self.name == other.name
            and self.unit == other.unit
            and self.description == other.description
            and _values_equal(self.value, other.value)
        )

    def __hash__(self) -> int:
        return hash((self.name, self.shape, str(self.unit)))

    def __repr__(self) -> str:
        unit = "" if self.unit is None else f", unit={self.unit}"
        return f"Buffer({self.name!r}, shape={self.shape}, dtype={self.dtype}{unit})"


class BufferSet:
    """An ordered, named, immutable collection of :class:`Buffer`\\ s.

    Deliberately *not* a :class:`ParameterSet`: there is no ``lnprior``, no
    ``pack``/``unpack`` and no prior transform, because a buffer has no prior
    and occupies no sampler dimension. The two collections meeting only at
    :meth:`Parameterised.context` is the contract's structural statement that
    parameters and buffers are different kinds of thing.

    Examples
    --------
    >>> import numpy as np
    >>> buffers = BufferSet([Buffer("wavelength", np.array([1.0, 2.0, 3.0]))])
    >>> buffers.names
    ('wavelength',)
    >>> "wavelength" in buffers, len(buffers)
    (True, 1)
    """

    __slots__ = ("_buffers", "_index")

    def __init__(self, buffers: Iterable[Buffer] = ()) -> None:
        items = tuple(buffers)
        for item in items:
            if not isinstance(item, Buffer):
                raise ParameterError(f"BufferSet takes Buffer instances, got {item!r}")
        index: dict[str, int] = {}
        for position, item in enumerate(items):
            if item.name in index:
                raise ParameterError(f"duplicate buffer name {item.name!r}")
            index[item.name] = position
        self._buffers = items
        self._index = index

    @property
    def buffers(self) -> tuple[Buffer, ...]:
        """The buffers, in declaration order."""
        return self._buffers

    @property
    def names(self) -> tuple[str, ...]:
        """Buffer names, in declaration order."""
        return tuple(b.name for b in self._buffers)

    def __len__(self) -> int:
        return len(self._buffers)

    def __iter__(self) -> Iterator[Buffer]:
        return iter(self._buffers)

    def __contains__(self, name: object) -> bool:
        return name in self._index

    def __getitem__(self, key: str | int) -> Buffer:
        if isinstance(key, str):
            try:
                return self._buffers[self._index[key]]
            except KeyError:
                raise KeyError(
                    f"no buffer named {key!r}; declared buffers are {list(self.names)}"
                ) from None
        return self._buffers[key]

    def values(self) -> dict[str, np.ndarray]:
        """A name-to-array mapping of every buffer."""
        return {b.name: b.array for b in self._buffers}

    def with_buffer(self, buffer: Buffer) -> BufferSet:
        """A copy with *buffer* appended (or replacing a same-named one)."""
        if buffer.name in self._index:
            replaced = list(self._buffers)
            replaced[self._index[buffer.name]] = buffer
            return BufferSet(replaced)
        return BufferSet([*self._buffers, buffer])

    def without(self, name: str) -> BufferSet:
        """A copy with the named buffer removed."""
        if name not in self._index:
            raise KeyError(f"no buffer named {name!r}; declared buffers are {list(self.names)}")
        return BufferSet(b for b in self._buffers if b.name != name)

    def __repr__(self) -> str:
        return f"BufferSet({list(self.names)})"


# ---------------------------------------------------------------------------
# Plates
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class Plate:
    """N members sharing hyperparameters — the population-model declaration.

    A ``Plate`` is a *constructor*, not a container: :meth:`expand` turns it
    into ordinary :class:`Parameter`\\ s, which is how it stays inspectable and
    how everything else in this module (pack/unpack, priors, merging) works on
    it without special cases. Its hyperparameters become
    ``"{plate}.{name}"``; each member becomes a single array-valued parameter
    of shape ``(size,) + member.shape``, tagged with :attr:`Parameter.plate`.

    That shape is chosen because it *is* the numpyro lowering: one
    ``numpyro.sample`` inside one ``numpyro.plate(name, size)`` produces a
    batch of ``size`` draws, not ``size`` separate sites. W1.9 owns the actual
    lowering; this contract's obligation is to hand it a declaration with the
    plate name, the plate size, and per-member priors whose hyperparameter
    references are already resolved to names in the enclosing set.

    Ampere v1.3 supports one plate dimension per parameter; nested plates are
    out of scope (see ``docs/design/contracts/parameters.md``).

    Examples
    --------
    >>> import scipy.stats
    >>> objects = Plate(
    ...     "objects",
    ...     size=3,
    ...     hyperparameters=[
    ...         Parameter("mu", scipy.stats.norm(0.0, 5.0)),
    ...         Parameter("sigma", scipy.stats.halfnorm(0.0, 2.0)),
    ...     ],
    ...     members=[
    ...         Parameter("theta", HierarchicalPrior("norm", {"loc": "mu", "scale": "sigma"})),
    ...     ],
    ... )
    >>> [p.name for p in objects.expand()]
    ['objects.mu', 'objects.sigma', 'objects.theta']
    >>> theta = objects.expand()[-1]
    >>> theta.shape, theta.plate, theta.references
    ((3,), 'objects', ('objects.mu', 'objects.sigma'))
    """

    name: str
    size: int
    hyperparameters: Sequence[Parameter] = ()
    members: Sequence[Parameter] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "name", _check_local_name(self.name, "plate"))
        if not isinstance(self.size, (int, np.integer)) or int(self.size) < 1:
            raise ParameterError(f"plate {self.name!r}: size must be an integer >= 1")
        object.__setattr__(self, "size", int(self.size))
        object.__setattr__(self, "hyperparameters", tuple(self.hyperparameters))
        object.__setattr__(self, "members", tuple(self.members))
        if not self.members:
            raise ParameterError(
                f"plate {self.name!r} declares no members; a plate with only hyperparameters "
                f"is just a group of ordinary parameters."
            )
        seen: set[str] = set()
        for parameter in (*self.hyperparameters, *self.members):
            if not isinstance(parameter, Parameter):
                raise ParameterError(f"plate {self.name!r} takes Parameters, got {parameter!r}")
            _check_local_name(parameter.name, "plate member")
            if parameter.name in seen:
                raise ParameterError(f"plate {self.name!r}: duplicate name {parameter.name!r}")
            seen.add(parameter.name)

    @property
    def hyperparameter_names(self) -> tuple[str, ...]:
        """Local hyperparameter names (before qualification)."""
        return tuple(p.name for p in self.hyperparameters)

    @property
    def member_names(self) -> tuple[str, ...]:
        """Local member names (before qualification)."""
        return tuple(p.name for p in self.members)

    def qualified(self, local_name: str) -> str:
        """The name *local_name* takes in the enclosing :class:`ParameterSet`."""
        return f"{self.name}{SEPARATOR}{local_name}"

    def expand(self) -> tuple[Parameter, ...]:
        """Materialise this plate as ordinary parameters."""
        rename = {p.name: self.qualified(p.name) for p in self.hyperparameters}
        expanded: list[Parameter] = []
        for hyper in self.hyperparameters:
            prior = hyper.prior
            if isinstance(prior, HierarchicalPrior):
                prior = prior.rename_references(rename)
            expanded.append(
                dataclasses.replace(hyper, name=self.qualified(hyper.name), prior=prior)
            )
        for member in self.members:
            if member.plate is not None and member.plate != self.name:
                raise ParameterError(
                    f"plate {self.name!r}: member {member.name!r} is already tagged with plate "
                    f"{member.plate!r}; nested plates are out of scope for v1.3."
                )
            prior = member.prior
            if isinstance(prior, HierarchicalPrior):
                prior = prior.rename_references(rename)
            expanded.append(
                dataclasses.replace(
                    member,
                    name=self.qualified(member.name),
                    prior=prior,
                    shape=(self.size, *member.shape),
                    plate=self.name,
                    value=None if member.value is None else member.value,
                )
            )
        return tuple(expanded)

    def to_parameter_set(self) -> ParameterSet:
        """This plate on its own, as a :class:`ParameterSet`."""
        return ParameterSet(self.expand())


# ---------------------------------------------------------------------------
# Tying and merging
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class Tie:
    """An instruction to collapse several parameter sites into one free parameter.

    Declaration-time tying (:attr:`Parameter.shared_as`) covers the case where
    the person writing the model knows the quantity is shared. ``Tie`` covers
    the other, equally common case: the models were written independently — or
    came from a library — and the person *composing the fit* is the one who
    knows that two distances are the same distance.

    Parameters
    ----------
    name
        The name the collapsed parameter takes in the merged set. It lives in
        the merged set's global namespace, not under any component prefix.
    sites
        Qualified ``"component.local_name"`` names to collapse.
    prior
        Optional override. If omitted, the sites' own priors must agree.

    Examples
    --------
    >>> Tie("distance", ("sed.distance", "spectrum.distance")).sites
    ('sed.distance', 'spectrum.distance')
    """

    name: str
    sites: Sequence[str]
    prior: AnyPrior = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "name", _check_name(self.name, "tie"))
        sites = tuple(self.sites)
        if len(sites) < 2:
            raise TyingError(
                f"tie {self.name!r} lists {len(sites)} site(s); a tie collapses two or more."
            )
        if len(set(sites)) != len(sites):
            raise TyingError(f"tie {self.name!r} lists a site more than once: {sites}")
        for site in sites:
            _check_name(site, "tie site")
        object.__setattr__(self, "sites", sites)


@dataclasses.dataclass(frozen=True)
class Binding:
    """Where one merged parameter is consumed.

    A tied parameter has several bindings; an untied one has exactly one. This
    is the record that makes tying *structural* rather than a constraint
    applied after the fact: the merged set has one free dimension, and the
    bindings say which component sees it under which local name.

    ``index`` (ruled 2026-09-02, ``hierarchical_population.md`` gap H-2) is
    the optional element address: when set, :meth:`ParameterMapping.distribute`
    routes ``resolved[global_name][index]`` — one element of an array-valued
    parameter — rather than the whole array, so a component consuming one
    member of a :class:`Plate` receives a scalar under its own local name.
    The merged set is unchanged (one array-valued parameter, one sample
    site — the lowering is unaffected); only the routing changes. This is
    *addressing*, not tying: each member remains its own draw, which is why
    limitation §12.3's refusal of ties across plate members stands.
    """

    global_name: str
    component: str
    local_name: str
    index: int | tuple[int, ...] | None = None


@dataclasses.dataclass(frozen=True)
class PlateBinding:
    """Declare, at merge time, that a component consumes one plate element.

    The composition-time counterpart of :class:`Binding.index` (ruled
    2026-09-02): ``ParameterSet.merge(..., plate_bindings=[...])`` turns each
    of these into an element :class:`Binding`. ``parameter`` is the **fully
    qualified merged name** of the array-valued parameter (after
    qualification and tie collapse, e.g. ``"population.objects.theta"``) —
    fully qualified rather than bare, because two components may each hold a
    plate of the same local name. ``local_name`` is the bare name under which
    the element arrives in the receiving component's ``distribute`` output;
    it must not collide with any name that component already receives, and
    the receiving component's own :class:`ParameterSet` does **not** declare
    it — consuming the routed element is the composing caller's contract.
    """

    parameter: str
    component: str
    local_name: str
    index: int | tuple[int, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "parameter", _check_name(self.parameter, "plate binding"))
        _check_local_name(self.component, "plate binding component")
        _check_local_name(self.local_name, "plate binding local name")
        raw = self.index
        index: int | tuple[int, ...]
        if isinstance(raw, (int, np.integer)) and not isinstance(raw, bool):
            index = int(raw)
        elif (
            isinstance(raw, tuple)
            and raw
            and all(isinstance(i, (int, np.integer)) and not isinstance(i, bool) for i in raw)
        ):
            index = tuple(int(i) for i in raw)
        else:
            raise ParameterError(
                f"plate binding for {self.parameter!r} has index {raw!r}; an element address "
                f"is an int or a non-empty tuple of ints."
            )
        object.__setattr__(self, "index", index)


@dataclasses.dataclass(frozen=True)
class ParameterMapping:
    """The result of :meth:`ParameterSet.merge`: a joint set plus its wiring.

    Attributes
    ----------
    merged
        The joint :class:`ParameterSet` an inference engine samples.
    routing
        One ``(global_name, component, local_name)`` triple per immediate
        consumer — the table :meth:`distribute` walks, one level deep.
    components
        Component labels, in merge order.
    inner
        The retained :class:`ParameterMapping` of every component that was
        merged *as a mapping* (lossless nesting, ruled 2026-09-02). Empty for
        components merged as plain sets.

    The public :attr:`bindings` view **composes** :attr:`routing` through
    :attr:`inner`, so it enumerates the ultimate leaves: a consumer walking
    ``bindings`` for provenance or labelling sees a parameter collapsed by an
    inner merge as the several leaf sites it really drives, not as the single
    local name the routing table needs. Routing and introspection are thereby
    two views of one structure — value flow stays level-by-level (each
    component receives its immediate-level names, so an inner composite
    re-distributes with its own retained mapping), while introspection tells
    the leaf-level truth.
    """

    merged: ParameterSet
    routing: tuple[Binding, ...]
    components: tuple[str, ...]
    inner: Mapping[str, ParameterMapping] = dataclasses.field(
        default_factory=lambda: types.MappingProxyType({})
    )

    @property
    def bindings(self) -> tuple[Binding, ...]:
        """Every leaf consumer of each merged parameter, composed through ``inner``.

        For a mapping with no inner mappings this is exactly :attr:`routing`.
        For a nested one, each routing entry whose component was merged as a
        mapping expands to that component's own (recursively composed) leaf
        bindings, with dotted local paths — so an inner ``shared_as`` collapse
        surfaces here as several bindings under the outer merged name. An
        element binding (``index`` set) does not descend: the element is the
        leaf.
        """
        out: list[Binding] = []
        for binding in self.routing:
            mapping = self.inner.get(binding.component)
            if mapping is None or binding.index is not None:
                out.append(binding)
                continue
            leaves = [b for b in mapping.bindings if b.global_name == binding.local_name]
            if not leaves:
                out.append(binding)
                continue
            out.extend(
                Binding(
                    global_name=binding.global_name,
                    component=binding.component,
                    local_name=f"{leaf.component}{SEPARATOR}{leaf.local_name}",
                    index=leaf.index,
                )
                for leaf in leaves
            )
        return tuple(out)

    def sites_of(self, global_name: str) -> tuple[Binding, ...]:
        """Every leaf binding fed by the merged parameter *global_name*."""
        return tuple(b for b in self.bindings if b.global_name == global_name)

    def global_name_for(self, component: str, local_name: str) -> str:
        """The merged name a component's local parameter was collapsed into.

        Accepts either the immediate local name (a routing entry) or a fully
        qualified leaf path (a composed binding).
        """
        for binding in (*self.routing, *self.bindings):
            if binding.component == component and binding.local_name == local_name:
                return binding.global_name
        raise KeyError(f"no binding for {component!r}.{local_name!r}")

    @property
    def tied_names(self) -> tuple[str, ...]:
        """Merged names driving more than one leaf site, in merged order.

        Composed through ``inner`` (lossless nesting), so a parameter
        collapsed by an inner merge is reported here too. Element bindings
        (``index`` set) do not count towards sharing: each element is its own
        draw, and a plate routed to N components is addressing, not tying.
        """
        counts: dict[str, int] = {}
        for binding in self.bindings:
            if binding.index is not None:
                continue
            counts[binding.global_name] = counts.get(binding.global_name, 0) + 1
        return tuple(name for name in self.merged.names if counts.get(name, 0) > 1)

    def distribute(self, values: Mapping[str, Value] | np.ndarray) -> dict[str, dict[str, Value]]:
        """Route merged values back to each component's local names.

        Accepts either a flat free-parameter vector or a name-to-value mapping
        over the merged set, and returns ``{component: {local_name: value}}``
        — exactly the keyword arguments each component model expects, with
        tied parameters appearing (identically) in every component that binds
        them, and fixed parameters included. Routing is one level deep by
        design (the nesting rule): a component merged as a mapping receives
        its own merged names and re-distributes with its retained mapping. An
        element binding hands the component the addressed element alone.
        """
        resolved = values if isinstance(values, Mapping) else self.merged.unpack(values)
        routed: dict[str, dict[str, Value]] = {component: {} for component in self.components}
        for binding in self.routing:
            try:
                value = resolved[binding.global_name]
            except KeyError as exc:
                raise ParameterError(
                    f"no value supplied for merged parameter {binding.global_name!r}"
                ) from exc
            if binding.index is not None:
                value = np.asarray(value)[binding.index]
            routed[binding.component][binding.local_name] = value
        return routed


# ---------------------------------------------------------------------------
# ParameterSet
# ---------------------------------------------------------------------------


class ParameterSet:
    """An ordered, named, immutable collection of :class:`Parameter`\\ s.

    This is the object an inference engine talks to. It owns the flat-vector
    layout (free parameters only, in declaration order, arrays flattened in C
    order), the joint prior, the unit-cube transform used by nested samplers,
    the unconstrained-space round trip used by gradient-based samplers, and
    neutral (de)serialisation.

    Parameters
    ----------
    parameters
        The parameters, in declaration order.
    plates
        :class:`Plate` declarations, expanded and appended after
        *parameters*.

    Raises
    ------
    ParameterError
        On duplicate names, on a hierarchical prior referencing a name that is
        not in the set, or on a cycle among hierarchical references.

    Examples
    --------
    >>> import scipy.stats
    >>> pset = ParameterSet([
    ...     Parameter("temperature", scipy.stats.uniform(100.0, 9900.0)),
    ...     Parameter("log_tau", scipy.stats.norm(0.0, 1.0)),
    ...     Parameter("distance", value=1.5, fixed=True),
    ... ])
    >>> pset.names
    ('temperature', 'log_tau', 'distance')
    >>> pset.free_names, pset.free_size
    (('temperature', 'log_tau'), 2)
    >>> pset.unpack(pset.prior_transform(np.array([0.5, 0.5])))
    {'temperature': 5050.0, 'log_tau': 0.0, 'distance': 1.5}
    """

    __slots__ = ("_free_size", "_index", "_order", "_parameters", "_slices")

    def __init__(
        self, parameters: Iterable[Parameter] = (), *, plates: Iterable[Plate] = ()
    ) -> None:
        items = list(parameters)
        for plate in plates:
            if not isinstance(plate, Plate):
                raise ParameterError(f"plates must be Plate instances, got {plate!r}")
            items.extend(plate.expand())
        for item in items:
            if not isinstance(item, Parameter):
                raise ParameterError(f"ParameterSet takes Parameter instances, got {item!r}")

        index: dict[str, int] = {}
        for position, item in enumerate(items):
            if item.name in index:
                raise ParameterError(
                    f"duplicate parameter name {item.name!r}; names must be unique within a "
                    f"ParameterSet (use ParameterSet.merge to compose sets that share names)."
                )
            index[item.name] = position

        slices: dict[str, slice] = {}
        offset = 0
        for item in items:
            if item.is_fixed:
                continue
            slices[item.name] = slice(offset, offset + item.size)
            offset += item.size

        self._parameters = tuple(items)
        self._index = index
        self._slices = slices
        self._free_size = offset
        self._order = _evaluation_order(self._parameters, index)

    # -- collection surface ------------------------------------------------

    @property
    def parameters(self) -> tuple[Parameter, ...]:
        """The parameters, in declaration order."""
        return self._parameters

    @property
    def names(self) -> tuple[str, ...]:
        """Every parameter name, in declaration order."""
        return tuple(p.name for p in self._parameters)

    @property
    def free_names(self) -> tuple[str, ...]:
        """Names occupying flat-vector dimensions (free *and* deferred)."""
        return tuple(p.name for p in self._parameters if not p.is_fixed)

    @property
    def fixed_names(self) -> tuple[str, ...]:
        """Names of parameters held at a value."""
        return tuple(p.name for p in self._parameters if p.is_fixed)

    @property
    def deferred_names(self) -> tuple[str, ...]:
        """Names still waiting for a tie group to supply their prior."""
        return tuple(p.name for p in self._parameters if p.is_deferred)

    @property
    def is_resolved(self) -> bool:
        """Whether every non-fixed parameter knows its own prior."""
        return not self.deferred_names

    @property
    def free_size(self) -> int:
        """Length of the flat vector: the number of sampler dimensions.

        Not the same as ``len(self)``, which counts *parameters* (fixed ones
        included, array-valued ones once each).
        """
        return self._free_size

    @property
    def plates(self) -> dict[str, int]:
        """Plate name to plate size, derived from the parameters that carry one."""
        found: dict[str, int] = {}
        for parameter in self._parameters:
            if parameter.plate is not None and parameter.shape:
                found.setdefault(parameter.plate, parameter.shape[0])
        return found

    def __len__(self) -> int:
        return len(self._parameters)

    def __iter__(self) -> Iterator[Parameter]:
        return iter(self._parameters)

    def __contains__(self, name: object) -> bool:
        return name in self._index

    def __getitem__(self, key: str | int) -> Parameter:
        if isinstance(key, str):
            try:
                return self._parameters[self._index[key]]
            except KeyError:
                raise KeyError(
                    f"no parameter named {key!r}; declared parameters are {list(self.names)}"
                ) from None
        return self._parameters[key]

    def __repr__(self) -> str:
        return (
            f"ParameterSet({len(self._parameters)} parameters, "
            f"{self._free_size} free dimensions: {list(self.names)})"
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ParameterSet):
            return NotImplemented
        return self._parameters == other._parameters

    # -- layout ------------------------------------------------------------

    def free_slice(self, name: str) -> slice:
        """The flat-vector slice a free parameter occupies."""
        try:
            return self._slices[name]
        except KeyError:
            if name in self._index:
                raise KeyError(f"parameter {name!r} is fixed and occupies no free slice") from None
            raise KeyError(f"no parameter named {name!r}") from None

    def free_labels(self) -> tuple[str, ...]:
        """One label per flat-vector entry: ``"offset[2]"`` for array elements.

        Intended for corner plots, ArviZ coordinates and any other consumer that
        needs a name per sampler dimension (W1.8).
        """
        labels: list[str] = []
        for parameter in self._parameters:
            if parameter.is_fixed:
                continue
            if not parameter.shape:
                labels.append(parameter.name)
            else:
                for flat in range(parameter.size):
                    index = np.unravel_index(flat, parameter.shape)
                    joined = ",".join(str(int(i)) for i in index)
                    labels.append(f"{parameter.name}[{joined}]")
        return tuple(labels)

    # -- pack / unpack -----------------------------------------------------

    def pack(self, values: Mapping[str, Value]) -> np.ndarray:
        """Flatten a name-to-value mapping into the free-parameter vector.

        Entries for fixed parameters are accepted and ignored, so that
        ``pack(unpack(theta))`` round-trips. Unknown names are an error (they
        are almost always typos), and every free parameter must be present.
        """
        unknown = set(values) - set(self._index)
        if unknown:
            raise ParameterError(
                f"pack() got value(s) for unknown parameter(s) {sorted(unknown)}; this set "
                f"declares {list(self.names)}."
            )
        flat = np.empty(self._free_size, dtype=float)
        for parameter in self._parameters:
            if parameter.is_fixed:
                continue
            try:
                raw = values[parameter.name]
            except KeyError:
                raise ParameterError(
                    f"pack() is missing a value for free parameter {parameter.name!r}"
                ) from None
            array = np.asarray(parameter.to_value(raw), dtype=float)
            if array.size != parameter.size:
                raise ParameterError(
                    f"parameter {parameter.name!r} expects {parameter.size} value(s) "
                    f"(shape {parameter.shape}), got {array.size}"
                )
            flat[self._slices[parameter.name]] = array.reshape(-1)
        return flat

    def unpack(self, theta: ArrayLike) -> dict[str, Value]:
        """Expand a free-parameter vector into a mapping over **all** parameters.

        Fixed parameters are injected at their declared values, so a model
        always receives its full vocabulary and never has to know which
        parameters this particular fit chose to vary.
        """
        flat = np.asarray(theta, dtype=float).reshape(-1)
        if flat.size != self._free_size:
            raise ParameterError(
                f"expected a flat vector of length {self._free_size} "
                f"({list(self.free_names)}), got length {flat.size}"
            )
        values: dict[str, Value] = {}
        for parameter in self._parameters:
            if parameter.is_fixed:
                values[parameter.name] = parameter.value
            else:
                chunk = flat[self._slices[parameter.name]]
                values[parameter.name] = (
                    chunk.reshape(parameter.shape).copy() if parameter.shape else float(chunk[0])
                )
        return values

    def complete(self, values: Mapping[str, Value] | ArrayLike) -> dict[str, Value]:
        """Normalise *values* into a complete name-to-value mapping.

        Accepts a flat free-parameter vector (delegating to :meth:`unpack`) or
        a partial mapping, in which case fixed parameters are filled in from
        their declared values. Anything still missing is an error.
        """
        if isinstance(values, Mapping):
            unknown = set(values) - set(self._index)
            if unknown:
                raise ParameterError(
                    f"got value(s) for unknown parameter(s) {sorted(unknown)}; this set "
                    f"declares {list(self.names)}."
                )
            filled = {p.name: p.value for p in self._parameters if p.is_fixed}
            filled.update(values)
            missing = [p.name for p in self._parameters if p.name not in filled]
            if missing:
                raise ParameterError(f"no value supplied for parameter(s) {missing}")
            return filled
        return self.unpack(values)

    # -- priors ------------------------------------------------------------

    def _require_resolved(self, operation: str) -> None:
        if not self.is_resolved:
            raise TyingError(
                f"cannot {operation}: parameter(s) {list(self.deferred_names)} are declared "
                f"shared (shared_as=...) without a prior, so their priors come from another "
                f"site. Merge this set with the one that supplies them "
                f"(ParameterSet.merge) before evaluating priors."
            )

    def lnprior(self, values: Mapping[str, Value] | ArrayLike) -> float:
        """Joint log prior at *values* (a flat vector or a name-to-value mapping).

        Fixed parameters contribute nothing. Array-valued parameters contribute
        the sum over their elements (the prior is i.i.d. across elements).
        Hierarchical priors are bound to the current values of the parameters
        they reference. Returns ``-inf`` as soon as any contribution is
        non-finite, without evaluating the rest.
        """
        self._require_resolved("evaluate lnprior")
        resolved = self.complete(values)
        total = 0.0
        for name in self._order:
            parameter = self[name]
            if parameter.is_fixed:
                continue
            prior = parameter.prior
            frozen = prior.bind(resolved) if isinstance(prior, HierarchicalPrior) else prior
            contribution = float(np.sum(log_density(frozen, resolved[name])))
            if not math.isfinite(contribution):
                return -math.inf
            total += contribution
        return total

    def prior_transform(self, unit_cube: ArrayLike) -> np.ndarray:
        """Map a unit-hypercube point to a free-parameter vector (nested sampling).

        Hierarchical priors are handled by evaluating parameters in dependency
        order, so a hyperparameter is always transformed before anything whose
        prior references it.
        """
        self._require_resolved("evaluate prior_transform")
        cube = np.asarray(unit_cube, dtype=float).reshape(-1)
        if cube.size != self._free_size:
            raise ParameterError(
                f"expected a unit-cube vector of length {self._free_size}, got {cube.size}"
            )
        flat = np.empty(self._free_size, dtype=float)
        resolved: dict[str, Value] = {p.name: p.value for p in self._parameters if p.is_fixed}
        for name in self._order:
            parameter = self[name]
            if parameter.is_fixed:
                continue
            prior = parameter.prior
            frozen = prior.bind(resolved) if isinstance(prior, HierarchicalPrior) else prior
            chunk = cube[self._slices[name]]
            drawn = np.asarray(frozen.ppf(chunk), dtype=float)
            flat[self._slices[name]] = drawn.reshape(-1)
            resolved[name] = (
                drawn.reshape(parameter.shape).copy()
                if parameter.shape
                else float(drawn.reshape(-1)[0])
            )
        return flat

    def sample(self, rng: np.random.Generator | None = None) -> dict[str, Value]:
        """Draw one point from the joint prior, as a name-to-value mapping.

        Used for initialisation, for SBI's simulation budget, and for tests.
        Dependency order is respected, so hierarchical draws are correct.
        """
        self._require_resolved("sample from the prior")
        generator = np.random.default_rng() if rng is None else rng
        resolved: dict[str, Value] = {p.name: p.value for p in self._parameters if p.is_fixed}
        for name in self._order:
            parameter = self[name]
            if parameter.is_fixed:
                continue
            prior = parameter.prior
            frozen = prior.bind(resolved) if isinstance(prior, HierarchicalPrior) else prior
            size = parameter.shape if parameter.shape else None
            drawn = np.asarray(frozen.rvs(size=size, random_state=generator), dtype=float)
            resolved[name] = drawn if parameter.shape else float(drawn)
        return {p.name: resolved[p.name] for p in self._parameters}

    # -- unconstrained space ----------------------------------------------

    def bijections(self) -> tuple[Bijection, ...]:
        """The unconstraining bijection of each free parameter, in flat order."""
        return tuple(p.unconstraining_bijection() for p in self._parameters if not p.is_fixed)

    def unconstrain(self, values: Mapping[str, Value] | ArrayLike) -> np.ndarray:
        """Map a free-parameter vector (or mapping) into unconstrained space."""
        self._require_resolved("map to unconstrained space")
        flat = self.pack(self.complete(values))
        out = np.empty_like(flat)
        for parameter in self._parameters:
            if parameter.is_fixed:
                continue
            where = self._slices[parameter.name]
            bijection = parameter.unconstraining_bijection()
            out[where] = np.asarray(bijection.unconstrain(flat[where]), dtype=float).reshape(-1)
        return out

    def constrain(self, unconstrained: ArrayLike) -> np.ndarray:
        """Map an unconstrained vector back to a free-parameter vector."""
        self._require_resolved("map from unconstrained space")
        y = np.asarray(unconstrained, dtype=float).reshape(-1)
        if y.size != self._free_size:
            raise ParameterError(
                f"expected an unconstrained vector of length {self._free_size}, got {y.size}"
            )
        out = np.empty_like(y)
        for parameter in self._parameters:
            if parameter.is_fixed:
                continue
            where = self._slices[parameter.name]
            bijection = parameter.unconstraining_bijection()
            out[where] = np.asarray(bijection.constrain(y[where]), dtype=float).reshape(-1)
        return out

    def lnprior_unconstrained(self, unconstrained: ArrayLike) -> float:
        """Log prior density in unconstrained space, Jacobian term included.

        This is the density a gradient-based sampler (NUTS/HMC) works with. It
        is stated here, on the reference path, so that the torch and jax
        lowerings have an oracle to agree with rather than each rediscovering
        the change-of-variables term.
        """
        y = np.asarray(unconstrained, dtype=float).reshape(-1)
        constrained = self.constrain(y)
        total = self.lnprior(constrained)
        if not math.isfinite(total):
            return -math.inf
        for parameter in self._parameters:
            if parameter.is_fixed:
                continue
            where = self._slices[parameter.name]
            bijection = parameter.unconstraining_bijection()
            total += float(np.sum(bijection.log_abs_det_jacobian(y[where])))
        return total if math.isfinite(total) else -math.inf

    # -- composition -------------------------------------------------------

    def with_parameter(self, parameter: Parameter) -> ParameterSet:
        """A copy with *parameter* appended, or replacing a same-named one."""
        if parameter.name in self._index:
            replaced = list(self._parameters)
            replaced[self._index[parameter.name]] = parameter
            return ParameterSet(replaced)
        return ParameterSet([*self._parameters, parameter])

    def without(self, name: str) -> ParameterSet:
        """A copy with the named parameter removed."""
        if name not in self._index:
            raise KeyError(f"no parameter named {name!r}")
        return ParameterSet(p for p in self._parameters if p.name != name)

    @classmethod
    def merge(
        cls,
        sets: Mapping[str, ParameterSet | ParameterMapping],
        *,
        ties: Sequence[Tie] = (),
        plate_bindings: Sequence[PlateBinding] = (),
    ) -> ParameterMapping:
        """Compose several parameter sets into one joint set.

        Names are qualified with their component label
        (``"spectrum.temperature"``) so independently written models never
        collide. Tied sites — declared either at declaration time via
        :attr:`Parameter.shared_as` or at composition time via :class:`Tie` —
        collapse into a **single** merged parameter with several
        :class:`Binding`\\ s, so a tied pair costs one sampler dimension, not
        two-plus-a-constraint.

        Parameters
        ----------
        sets
            Component label to :class:`ParameterSet` — or to a
            :class:`ParameterMapping` (**lossless nesting**, ruled
            2026-09-02): the mapping's ``merged`` set joins the merge exactly
            as a plain set would, and the mapping is retained so the result's
            :attr:`~ParameterMapping.bindings` compose down to the leaves.
            Routing is unchanged either way. Labels must be Python
            identifiers; iteration order fixes the merged declaration order.
        ties
            Composition-time ties, in addition to any ``shared_as`` labels.
        plate_bindings
            :class:`PlateBinding` declarations (ruled 2026-09-02): each routes
            one element of a merged array-valued parameter to a component
            under a bare local name, via :attr:`Binding.index`. Validated
            here — the parameter must exist and be array-valued, the index in
            range, the component present, and the local name free.

        Returns
        -------
        ParameterMapping
            The joint set plus the wiring needed to route values back to
            components. Note that this is *not* a bare ``ParameterSet``: the
            bindings are the point, and losing them would reduce tying to a
            naming convention.

        Raises
        ------
        TyingError
            If tied sites disagree about shape, unit, prior or fixed value; if
            a tie names a site that does not exist; if a site is claimed by two
            ties; or if a tie group has no prior and no fixed value at all.
        ParameterError
            If a plate binding names an absent or scalar parameter, an
            out-of-range index, an unknown component, or a colliding local
            name.
        """
        return _merge(sets, ties, plate_bindings)

    # -- serialisation -----------------------------------------------------

    def to_spec(self) -> dict[str, Any]:
        """A JSON-compatible description of this set.

        Every prior must be describable (:func:`describe_prior`); an opaque
        prior raises here rather than producing a provenance record that
        silently omits it. Bijections are recorded only when they were declared
        explicitly, since an inferred one is a function of the prior and would
        be redundant.
        """
        entries: list[dict[str, Any]] = []
        for parameter in self._parameters:
            entry: dict[str, Any] = {"name": parameter.name}
            if isinstance(parameter.prior, HierarchicalPrior):
                entry["hierarchical_prior"] = parameter.prior.to_dict()
            elif parameter.prior is not None:
                entry["prior"] = describe_prior(parameter.prior).to_dict()
            if parameter.value is not None:
                entry["value"] = (
                    parameter.value.tolist()
                    if isinstance(parameter.value, np.ndarray)
                    else float(parameter.value)
                )
            if parameter.fixed:
                entry["fixed"] = True
            if parameter.shape:
                entry["shape"] = list(parameter.shape)
            if parameter.unit is not None:
                entry["unit"] = parameter.unit.to_string()
            if parameter.bijection is not None:
                entry["bijection"] = _bijection_to_dict(parameter.bijection)
            if parameter.shared_as is not None:
                entry["shared_as"] = parameter.shared_as
            if parameter.plate is not None:
                entry["plate"] = parameter.plate
            if parameter.description:
                entry["description"] = parameter.description
            entries.append(entry)
        return {"version": 1, "parameters": entries}

    @classmethod
    def from_spec(cls, spec: Mapping[str, Any]) -> ParameterSet:
        """Rebuild a :class:`ParameterSet` from :meth:`to_spec`'s output."""
        version = spec.get("version", 1)
        if version != 1:
            raise ParameterError(f"unsupported ParameterSet spec version {version!r}")
        parameters: list[Parameter] = []
        for entry in spec["parameters"]:
            prior: AnyPrior = None
            if "hierarchical_prior" in entry:
                prior = HierarchicalPrior.from_dict(entry["hierarchical_prior"])
            elif "prior" in entry:
                prior = prior_from_spec(PriorSpec.from_dict(entry["prior"]))
            parameters.append(
                Parameter(
                    entry["name"],
                    prior,
                    value=entry.get("value"),
                    fixed=bool(entry.get("fixed", False)),
                    shape=tuple(entry.get("shape", ())),
                    unit=u.Unit(entry["unit"]) if "unit" in entry else None,
                    bijection=(
                        _bijection_from_dict(entry["bijection"]) if "bijection" in entry else None
                    ),
                    shared_as=entry.get("shared_as"),
                    plate=entry.get("plate"),
                    description=entry.get("description", ""),
                )
            )
        return cls(parameters)

    # -- lowering extension point -----------------------------------------

    def as_paramax(self) -> Any:
        """Not implemented in core, and deliberately so.

        The harvested copilot scaffold carried an ``as_paramax()`` stub here.
        It does not belong in ``ampere.core``: ``DEVELOPMENT_PLAN.md`` §6 is
        explicit that "Paramax is a lowering mechanism, not a user API", and
        ``architecture.md`` §4 forbids core from knowing about optional
        dependencies at all. The jax backend lowers a ``ParameterSet`` from the
        neutral information this class already exposes — :class:`PriorSpec`,
        :class:`Bijection`, shape, fixed/free state, plate membership and
        buffers — per W1.9's lowering table.

        This method is kept, raising, only so that anyone porting the scaffold
        finds an explanation instead of an ``AttributeError``.
        """
        raise OptionalDependencyError(
            "paramax",
            extra="jax",
            context="lowering a ParameterSet to a paramax pytree (ampere.backends.jax's job, "
            "not ampere.core's — see docs/design/contracts/parameters.md)",
        )


def _bijection_to_dict(bijection: Bijection) -> dict[str, Any]:
    if isinstance(bijection, Identity):
        return {"kind": "identity"}
    if isinstance(bijection, Log):
        return {"kind": "log", "lower": bijection.lower}
    if isinstance(bijection, Logit):
        return {"kind": "logit", "lower": bijection.lower, "upper": bijection.upper}
    raise ParameterError(
        f"bijection {bijection!r} is not one of ampere's built-ins, so it cannot be serialised "
        f"neutrally. Built-ins are Identity, Log and Logit."
    )


def _bijection_from_dict(data: Mapping[str, Any]) -> Bijection:
    kind = data["kind"]
    if kind == "identity":
        return Identity()
    if kind == "log":
        return Log(lower=float(data.get("lower", 0.0)))
    if kind == "logit":
        return Logit(lower=float(data["lower"]), upper=float(data["upper"]))
    raise ParameterError(f"unknown bijection kind {kind!r}")


def _evaluation_order(parameters: Sequence[Parameter], index: Mapping[str, int]) -> tuple[str, ...]:
    """Topologically order parameters so hyperparameters precede their dependants.

    Declaration order is preserved among independent parameters, so the result
    is deterministic. The flat-vector layout is *not* reordered — only the
    order in which priors are evaluated.
    """
    dependencies: dict[str, tuple[str, ...]] = {}
    for parameter in parameters:
        for reference in parameter.references:
            if reference not in index:
                raise ParameterError(
                    f"parameter {parameter.name!r} has a hierarchical prior referencing "
                    f"{reference!r}, which is not in this set (it declares "
                    f"{sorted(index)}). Hierarchical references must resolve within the set "
                    f"that will evaluate them."
                )
            if reference == parameter.name:
                raise ParameterError(
                    f"parameter {parameter.name!r} has a hierarchical prior referencing itself"
                )
        dependencies[parameter.name] = parameter.references

    ordered: list[str] = []
    state: dict[str, int] = {}  # 0 = unvisited, 1 = in progress, 2 = done

    def visit(name: str, trail: tuple[str, ...]) -> None:
        mark = state.get(name, 0)
        if mark == 2:
            return
        if mark == 1:
            cycle = " -> ".join([*trail, name])
            raise ParameterError(
                f"hierarchical priors form a cycle: {cycle}. A parameter's prior cannot "
                f"depend (even indirectly) on itself."
            )
        state[name] = 1
        for reference in dependencies[name]:
            visit(reference, (*trail, name))
        state[name] = 2
        ordered.append(name)

    for parameter in parameters:
        visit(parameter.name, ())
    return tuple(ordered)


# ---------------------------------------------------------------------------
# Merging implementation
# ---------------------------------------------------------------------------


def _merge(
    sets: Mapping[str, ParameterSet | ParameterMapping],
    ties: Sequence[Tie],
    plate_bindings: Sequence[PlateBinding] = (),
) -> ParameterMapping:
    components = tuple(sets)
    for component in components:
        _check_local_name(component, "component label")

    # Lossless nesting (ruled 2026-09-02): a component given as a mapping
    # contributes its merged set to the algorithm below unchanged, and the
    # mapping itself is retained so the public bindings view composes to the
    # leaves. Nothing downstream of this loop knows the difference.
    inner: dict[str, ParameterMapping] = {}
    resolved_sets: dict[str, ParameterSet] = {}
    for component, given in sets.items():
        if isinstance(given, ParameterMapping):
            inner[component] = given
            resolved_sets[component] = given.merged
        else:
            resolved_sets[component] = given
    sets = resolved_sets

    qualified: dict[str, tuple[str, Parameter]] = {}
    for component in components:
        for parameter in sets[component]:
            qualified[f"{component}{SEPARATOR}{parameter.name}"] = (component, parameter)

    # Which merged parameter does each site belong to?
    group_of: dict[str, str] = {}
    claimed_by: dict[str, str] = {}

    for tie in ties:
        for site in tie.sites:
            if site not in qualified:
                raise TyingError(
                    f"tie {tie.name!r} names site {site!r}, which does not exist; available "
                    f"sites are {sorted(qualified)}."
                )
            if site in claimed_by:
                raise TyingError(
                    f"site {site!r} is claimed by both tie {claimed_by[site]!r} and tie "
                    f"{tie.name!r}; a site may belong to at most one tie."
                )
            claimed_by[site] = tie.name
            group_of[site] = tie.name

    for site, (_component, parameter) in qualified.items():
        if parameter.shared_as is None:
            continue
        if site in claimed_by:
            raise TyingError(
                f"site {site!r} is declared shared_as={parameter.shared_as!r} *and* named by "
                f"tie {claimed_by[site]!r}; that is ambiguous — remove one."
            )
        group_of[site] = parameter.shared_as

    for site in qualified:
        group_of.setdefault(site, site)

    tie_by_name = {tie.name: tie for tie in ties}

    # Preserve component-then-declaration order, emitting each group once.
    order: list[str] = []
    members: dict[str, list[tuple[str, Parameter]]] = {}
    for component in components:
        for parameter in sets[component]:
            site = f"{component}{SEPARATOR}{parameter.name}"
            group = group_of[site]
            if group not in members:
                members[group] = []
                order.append(group)
            members[group].append((component, parameter))

    rename: dict[str, dict[str, str]] = {component: {} for component in components}
    for group, sites in members.items():
        for component, parameter in sites:
            rename[component][parameter.name] = group

    merged: list[Parameter] = []
    bindings: list[Binding] = []
    for group in order:
        sites = members[group]
        merged.append(_collapse(group, sites, rename, tie_by_name.get(group)))
        bindings.extend(
            Binding(global_name=group, component=component, local_name=parameter.name)
            for component, parameter in sites
        )

    merged_set = ParameterSet(merged)
    taken = {(b.component, b.local_name) for b in bindings}
    for pb in plate_bindings:
        if pb.component not in resolved_sets:
            raise ParameterError(
                f"plate binding routes {pb.parameter!r} to component {pb.component!r}, which "
                f"this merge does not have; components are {sorted(resolved_sets)}."
            )
        if pb.parameter not in merged_set.names:
            raise ParameterError(
                f"plate binding names merged parameter {pb.parameter!r}, which does not exist. "
                f"The name is the fully qualified merged one (after qualification and tie "
                f"collapse), e.g. 'population.objects.theta'."
            )
        shape = merged_set[pb.parameter].shape
        if not shape:
            raise ParameterError(
                f"plate binding addresses element {pb.index!r} of {pb.parameter!r}, which is "
                f"scalar. Element routing is for array-valued (plate) parameters."
            )
        index = (pb.index,) if isinstance(pb.index, int) else pb.index
        if len(index) != len(shape) or any(
            not 0 <= i < extent for i, extent in zip(index, shape, strict=True)
        ):
            raise ParameterError(
                f"plate binding index {pb.index!r} is out of range for {pb.parameter!r}, whose "
                f"shape is {shape}."
            )
        if (pb.component, pb.local_name) in taken:
            raise ParameterError(
                f"plate binding would deliver {pb.parameter!r} to "
                f"{pb.component}{SEPARATOR}{pb.local_name}, but component {pb.component!r} "
                f"already receives a value under {pb.local_name!r}. Element routing must not "
                f"shadow a component's own parameter; pick a different local name."
            )
        taken.add((pb.component, pb.local_name))
        bindings.append(
            Binding(
                global_name=pb.parameter,
                component=pb.component,
                local_name=pb.local_name,
                index=pb.index,
            )
        )

    return ParameterMapping(
        merged=merged_set,
        routing=tuple(bindings),
        components=components,
        inner=types.MappingProxyType(inner),
    )


def _qualified_prior(parameter: Parameter, rename: Mapping[str, str]) -> AnyPrior:
    """A site's prior with hierarchical references mapped to merged names."""
    prior = parameter.prior
    if isinstance(prior, HierarchicalPrior):
        return prior.rename_references(rename)
    return prior


def _collapse(
    group: str,
    sites: Sequence[tuple[str, Parameter]],
    rename: Mapping[str, Mapping[str, str]],
    tie: Tie | None,
) -> Parameter:
    """Reduce one tie group to the single :class:`Parameter` the sampler sees."""
    first_component, first = sites[0]

    for component, parameter in sites[1:]:
        if parameter.shape != first.shape:
            raise TyingError(
                f"tied parameters disagree about shape: {first_component}.{first.name} has "
                f"{first.shape} but {component}.{parameter.name} has {parameter.shape}."
            )
        if parameter.unit != first.unit:
            raise TyingError(
                f"tied parameters disagree about units: {first_component}.{first.name} is in "
                f"{first.unit} but {component}.{parameter.name} is in {parameter.unit}. "
                f"Priors are declared numerically in the parameter's unit, so ampere will not "
                f"guess a conversion; declare both in the same unit."
            )
        # This loop body only runs for multi-site groups, i.e. genuine ties.
        if first.plate is not None or parameter.plate is not None:
            plated = first if first.plate is not None else parameter
            raise TyingError(
                f"tie {group!r} spans plate members ({plated.name} is in plate "
                f"{plated.plate!r}); tying across plates is out of scope for v1.3. Express "
                f"the shared quantity as a hyperparameter referenced by a HierarchicalPrior "
                f"instead."
            )

    fixed_sites = [(c, p) for c, p in sites if p.is_fixed]
    prior_sites = [(c, p) for c, p in sites if p.prior is not None]

    if tie is not None and tie.prior is not None:
        if fixed_sites:
            raise TyingError(
                f"tie {group!r} supplies a prior, but site "
                f"{fixed_sites[0][0]}.{fixed_sites[0][1].name} is fixed; a tie group is either "
                f"free or fixed, not both."
            )
        prior: AnyPrior = tie.prior
        if isinstance(prior, HierarchicalPrior):
            # A composer-supplied hierarchical prior is written against local
            # names; qualify it with the first component's mapping
            # (merged-global references pass through unchanged).
            prior = prior.rename_references(rename[first_component])
        fixed = False
        value: Value = first.value
    elif fixed_sites:
        if prior_sites:
            raise TyingError(
                f"tie {group!r} mixes fixed and prior-equipped sites "
                f"({fixed_sites[0][0]}.{fixed_sites[0][1].name} is fixed, "
                f"{prior_sites[0][0]}.{prior_sites[0][1].name} has a prior); decide whether the "
                f"shared quantity varies."
            )
        reference = fixed_sites[0][1]
        for component, parameter in fixed_sites[1:]:
            if not _values_equal(parameter.value, reference.value):
                raise TyingError(
                    f"tie {group!r} fixes the same quantity at different values: "
                    f"{fixed_sites[0][0]} says {reference.value!r}, {component} says "
                    f"{parameter.value!r}."
                )
        prior = None
        fixed = True
        value = reference.value
    else:
        if not prior_sites:
            raise TyingError(
                f"tie {group!r} has no prior: every site is declared shared without one. One "
                f"site must carry the prior, or the Tie must supply it."
            )
        # Hierarchical priors are compared *after* qualifying their references
        # with each site's own component mapping: two sites are only the same
        # prior if their hyperparameters resolve to the same merged parameters.
        # Comparing the raw declarations instead would let a tie silently wire
        # the collapsed parameter to the first component's hyperparameters.
        reference_component, reference = prior_sites[0]
        prior = _qualified_prior(reference, rename[reference_component])
        for component, parameter in prior_sites[1:]:
            candidate = _qualified_prior(parameter, rename[component])
            if not _priors_equal(candidate, prior):
                hint = (
                    " If the hyperparameters are themselves the same quantity, tie them too:"
                    " hierarchical references are compared after qualification."
                    if isinstance(prior, HierarchicalPrior)
                    and isinstance(candidate, HierarchicalPrior)
                    else " Pass an explicit Tie(prior=...) if you meant to override."
                )
                raise TyingError(
                    f"tie {group!r} has disagreeing priors: {reference_component}."
                    f"{reference.name} declares {_prior_repr(prior)} but "
                    f"{component}.{parameter.name} declares {_prior_repr(candidate)}. "
                    f"Tied parameters are one parameter and must have one prior.{hint}"
                )
        fixed = False
        value = reference.value

    declared_bijections = [(c, p) for c, p in sites if p.bijection is not None]
    bijection = declared_bijections[0][1].bijection if declared_bijections else None
    for component, parameter in declared_bijections[1:]:
        if parameter.bijection != bijection:
            raise TyingError(
                f"tied parameters declare different bijections: "
                f"{declared_bijections[0][0]}.{declared_bijections[0][1].name} declares "
                f"{bijection!r} but {component}.{parameter.name} declares "
                f"{parameter.bijection!r}. Declare it on one site, or identically on all."
            )
    description = next((p.description for _c, p in sites if p.description), "")

    return Parameter(
        group,
        prior,
        value=value,
        fixed=fixed,
        shape=first.shape,
        unit=first.unit,
        bijection=bijection,
        shared_as=None,
        plate=first.plate,
        description=description,
    )


# ---------------------------------------------------------------------------
# The declaration mixin
# ---------------------------------------------------------------------------


class Parameterised:
    """Mixin giving an object explicit parameters and explicit buffers.

    Any model-like object — a model, an instrument transformation, a noise
    model, a dataset — inherits this to declare what it owns. Registration is
    imperative (``register_parameter`` / ``register_buffer``) rather than
    class-level, because buffers usually come from data that only exists at
    construction time (a wavelength grid, an opacity table read from disk).

    The important method is :meth:`context`. An evaluation reads **both**
    parameters and buffers out of the single mapping it returns::

        def __call__(self, **values):
            ctx = self.context(values)
            return blackbody(ctx["wavelength"], ctx["temperature"])

    Written that way, moving a quantity between buffer, fixed parameter and
    free parameter is a configuration change (:meth:`promote_buffer`,
    :meth:`Parameter.fix`) and never a code change — which is precisely the
    test ``architecture.md`` §6 sets for this contract.

    A future ``Model`` ABC (W1.5/W1.7, once :class:`ModelResult` exists to type
    its return) inherits this rather than redeclaring it.

    Examples
    --------
    >>> import numpy as np, scipy.stats
    >>> class Blackbody(Parameterised):
    ...     def __init__(self, wavelength):
    ...         self.register_buffer("wavelength", wavelength)
    ...         self.register_parameter(
    ...             Parameter("temperature", scipy.stats.uniform(100.0, 9900.0))
    ...         )
    ...     def __call__(self, **values):
    ...         ctx = self.context(values)
    ...         return ctx["temperature"] * ctx["wavelength"]
    >>> model = Blackbody(np.array([1.0, 2.0]))
    >>> model.parameters.free_names, model.buffers.names
    (('temperature',), ('wavelength',))
    >>> model(temperature=2.0)
    array([2., 4.])
    """

    _parameters: ParameterSet
    _buffers: BufferSet

    # -- declared state ----------------------------------------------------

    @property
    def parameters(self) -> ParameterSet:
        """Everything this object declares as varying or fixed."""
        try:
            return self._parameters
        except AttributeError:
            self._parameters = ParameterSet()
            return self._parameters

    @property
    def buffers(self) -> BufferSet:
        """Everything this object declares as constant data."""
        try:
            return self._buffers
        except AttributeError:
            self._buffers = BufferSet()
            return self._buffers

    # -- declaration -------------------------------------------------------

    def _check_free_name(self, name: str) -> None:
        if name in self.parameters:
            raise ParameterError(f"{type(self).__name__} already declares a parameter {name!r}")
        if name in self.buffers:
            raise ParameterError(f"{type(self).__name__} already declares a buffer {name!r}")
        if hasattr(type(self), name):
            raise ParameterError(
                f"{name!r} shadows an attribute of {type(self).__name__}; pick another name."
            )

    def register_parameter(self, parameter: Parameter) -> Parameter:
        """Declare *parameter*. Returns it, for convenience."""
        if not isinstance(parameter, Parameter):
            raise ParameterError(f"register_parameter takes a Parameter, got {parameter!r}")
        self._check_free_name(parameter.name)
        self._parameters = self.parameters.with_parameter(parameter)
        return parameter

    def register_parameters(self, *parameters: Parameter) -> None:
        """Declare several parameters, in order."""
        for parameter in parameters:
            self.register_parameter(parameter)

    def register_buffer(
        self,
        name: str,
        value: ArrayLike,
        *,
        unit: u.UnitBase | None = None,
        description: str = "",
    ) -> np.ndarray:
        """Declare a constant array. Returns the stored (read-only) array.

        Explicit by design (``architecture.md`` §6): an array attribute that is
        never registered is simply not part of this object's declared state,
        and will not move with device/dtype changes or appear in provenance.
        Nothing here silently adopts stray attributes as buffers.
        """
        self._check_free_name(_check_name(name, "buffer"))
        buffer = Buffer(name, np.asarray(value), unit=unit, description=description)
        self._buffers = self.buffers.with_buffer(buffer)
        return buffer.array

    # -- promotion / demotion ---------------------------------------------

    def promote_buffer(
        self,
        name: str,
        *,
        prior: AnyPrior = None,
        bijection: Bijection | None = None,
    ) -> Parameter:
        """Turn a declared buffer into a parameter, in place.

        With no *prior* the result is a **fixed** parameter holding the
        buffer's value: the quantity joins the model's parameter vocabulary and
        its provenance record without yet varying. With a *prior* it becomes a
        free parameter. Either way the model's evaluation code is untouched, so
        long as it reads through :meth:`context`.
        """
        buffer = self.buffers[name]
        self._buffers = self.buffers.without(name)
        parameter = Parameter(
            buffer.name,
            prior,
            value=buffer.value,
            fixed=prior is None,
            unit=buffer.unit,
            bijection=bijection,
            description=buffer.description,
        )
        self._parameters = self.parameters.with_parameter(parameter)
        return parameter

    def demote_parameter(self, name: str, *, value: Value = None) -> Buffer:
        """Turn a declared parameter back into a buffer (the inverse of promotion)."""
        parameter = self.parameters[name]
        raw = parameter.value if value is None else parameter.to_value(value)
        if raw is None:
            raise ParameterError(
                f"parameter {name!r} has no value, so it cannot become a buffer; pass value=..."
            )
        self._parameters = self.parameters.without(name)
        buffer = Buffer(
            parameter.name,
            np.asarray(raw),
            unit=parameter.unit,
            description=parameter.description,
        )
        self._buffers = self.buffers.with_buffer(buffer)
        return buffer

    # -- evaluation --------------------------------------------------------

    def context(self, values: Mapping[str, Value] | ArrayLike | None = None) -> dict[str, Value]:
        """One namespace holding every buffer and every parameter value.

        Parameters
        ----------
        values
            A flat free-parameter vector, a name-to-value mapping (fixed
            parameters may be omitted), or ``None`` to use each parameter's
            declared value.

        Returns
        -------
        dict
            Buffer arrays and parameter values under their declared names.
            Buffer and parameter names are guaranteed disjoint, so the merge is
            unambiguous.
        """
        if values is None:
            missing = [p.name for p in self.parameters if p.value is None]
            if missing:
                raise ParameterError(
                    f"context() was given no values and parameter(s) {missing} have no declared "
                    f"value to fall back on."
                )
            resolved: dict[str, Value] = {p.name: p.value for p in self.parameters}
        else:
            resolved = self.parameters.complete(values)
        merged: dict[str, Value] = self.buffers.values()
        merged.update(resolved)
        return merged
