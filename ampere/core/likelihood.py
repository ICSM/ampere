"""Likelihood families and noise models: ``DEVELOPMENT_PLAN.md`` §4.4.

This is the contract in which ampere says *how well a prediction matches an
observation*. It has three separable pieces, and keeping them separable is the
whole point (``prior_art.md`` Tension 3: 3ML fuses instrument response and
likelihood into one opaque plugin, and buys extensibility at the cost of all
reuse; ampere factors them and pays for it with this spec):

- a :class:`LikelihoodFamily` — the sampling distribution of a single datum
  given a prediction and its noise. Registered through one method,
  ``log_prob(predicted, observed, noise)`` (``prior_art.md`` lesson 3M3:
  copy 3ML's minimalism, not its opacity), so a user adds a family without
  touching ampere;
- a :class:`NoiseModel` — what the noise *is*: independent per-sample
  uncertainties (:class:`IndependentNoise`) or a Gaussian process over the
  residuals (:class:`GaussianProcessNoise`), which is ampere's flexible,
  misspecification-robust likelihood;
- a :class:`GPSolver` — *how* the GP algebra is done. A swappable strategy:
  :class:`DenseGP` (exact, O(N³), the correctness anchor implemented here),
  :class:`QuasisepGP` (exact, O(N), Phase 2) and the approximate strategies
  for the cases neither covers.

:class:`Likelihood` composes the three and is the object a ``Dataset`` (W1.7)
holds.

Two structural facts drive everything below.

**Float64, always.** Every array this module touches is cast to
:data:`~numpy.float64` (or ``complex128``). ``DEVELOPMENT_PLAN.md`` §7 and
``architecture.md`` §5 are explicit that GP Cholesky and quasiseparable solves
in float32 fail in ways that read as science bugs rather than numerical ones,
and that this is not a per-backend default to inherit.

**Analytic marginalisation is a property of the family/noise *combination*,
not of either alone.** The "add the GP covariance to the noise and integrate"
trick only closes for a Gaussian noise process. For Poisson, Student-t and the
rest, robustness needs a latent-GP formulation (``counts ~ Poisson(rate ·
exp(f))``, ``f ~ GP``) marginalised numerically — so a :class:`Likelihood`
declares its :class:`Marginalisation`, and refuses an engine that cannot
deliver it (``DEVELOPMENT_PLAN.md`` §4.4, "the Likelihood interface must
therefore declare whether it marginalises its noise process analytically or
introduces latent variables").

Nothing here imports torch, jax or any optional dependency
(``architecture.md`` §4 rule 1): it is numpy, scipy, ``astropy.units`` and
stdlib. The narrative spec is ``docs/design/contracts/likelihoods.md``, whose
every example runs as a doctest.
"""

from __future__ import annotations

import abc
import dataclasses
import enum
import math
from collections.abc import Mapping, Sequence
from typing import Any, ClassVar

import astropy.units as u
import numpy as np
import scipy.linalg
import scipy.stats as st

from .exceptions import LikelihoodError
from .parameter import (
    Identity,
    Log,
    Parameter,
    Parameterised,
    ParameterSet,
)
from .results_schema import FunctionSamples, Layout

__all__ = [
    "DTYPE",
    "CauchyFamily",
    "Censoring",
    "ComplexGaussianFamily",
    "DenseGP",
    "GPConditional",
    "GPSolver",
    "GaussianFamily",
    "GaussianProcessNoise",
    "IndependentNoise",
    "InducingPointGP",
    "Kernel",
    "KernelSpec",
    "LatentDeclaration",
    "Likelihood",
    "LikelihoodFamily",
    "LimitKind",
    "Marginalisation",
    "Matern32",
    "NoiseModel",
    "NoiseParams",
    "PoissonFamily",
    "QuasisepGP",
    "RiceFamily",
    "SquaredExponential",
    "StructuredGridGP",
    "StudentTFamily",
    "VecchiaGP",
    "VonMisesFamily",
    "WindowedSparseGP",
    "family_named",
    "latent_parameter",
    "list_families",
    "register_family",
]

#: The dtype every array in this module is held in. ``DEVELOPMENT_PLAN.md`` §7:
#: GP linear algebra in float32 fails in ways that look like science problems.
DTYPE = np.float64

_LOG_2PI = math.log(2.0 * math.pi)
_SQRT3 = math.sqrt(3.0)


def _check_finite(array: np.ndarray, what: str) -> np.ndarray:
    """Refuse NaN or inf. Applies to complex arrays as well as real ones."""
    if not np.all(np.isfinite(array)):
        raise LikelihoodError(
            f"{what} contains non-finite entries. A likelihood cannot be evaluated on NaN or "
            f"inf; mask the affected samples (mask=True excludes them entirely) rather than "
            f"threading sentinels through the arithmetic."
        )
    return array


def _as_float64(array: Any, what: str) -> np.ndarray:
    """Cast to a contiguous float64 array, loudly."""
    return _check_finite(np.ascontiguousarray(np.asarray(array), dtype=DTYPE), what)


def _as_points(array: Any, what: str, *, dimensions: int | None = None) -> np.ndarray:
    """Coerce coordinates to an ``(n, d)`` float64 array of *n points*.

    ``np.atleast_2d`` is the wrong tool here and the reason this helper exists:
    it turns a shape ``(m,)`` array into ``(1, m)`` — **one m-dimensional
    point** — when what a caller passing a bare list of wavelengths means is m
    one-dimensional points. That reading silently broadcasts through the kernel
    and yields a one-element answer instead of an error, so a 1-D input is
    interpreted here as a column, explicitly, and anything ambiguous raises.
    """
    values = np.asarray(array, dtype=DTYPE)
    if values.ndim == 1:
        values = values[:, None]
    elif values.ndim != 2:
        raise LikelihoodError(
            f"{what} must be a 1-D array of coordinates or an (n, d) array of points, but it "
            f"has shape {values.shape}."
        )
    if dimensions is not None and values.shape[1] != dimensions:
        raise LikelihoodError(
            f"{what} has {values.shape[1]} coordinate(s) per point, but the data it is being "
            f"compared against have {dimensions}. Pass an (n, {dimensions}) array"
            + (", or a bare 1-D array of coordinates." if dimensions == 1 else ".")
        )
    return _check_finite(np.ascontiguousarray(values), what)


# ---------------------------------------------------------------------------
# Declarations: marginalisation and censoring
# ---------------------------------------------------------------------------


class Marginalisation(enum.Enum):
    """How a family/noise combination deals with its noise process.

    ``DEVELOPMENT_PLAN.md`` §4.4 requires this to be *declared*, because it
    decides which inference engines can run the problem at all.
    """

    #: The noise process is integrated out in closed form. Any engine works.
    ANALYTIC = "analytic"
    #: The noise process enters as latent variables that inference must sample
    #: or approximate. Effectively a modern-backend (HMC/VI) capability.
    LATENT = "latent"


class LimitKind(enum.IntEnum):
    """What a single observed sample asserts about the truth.

    An :class:`enum.IntEnum` so a per-sample classification is an ordinary
    integer array — which is exactly the aligned-auxiliary-array hook
    ``results_schema.md`` §8 left open for this contract (issue #11).
    """

    #: An ordinary measurement: the value is the datum.
    DETECTION = 0
    #: A non-detection: the truth is *below* the recorded value.
    UPPER_LIMIT = 1
    #: The truth is *above* the recorded value (a saturated pixel, say).
    LOWER_LIMIT = 2


def _readonly(array: np.ndarray) -> np.ndarray:
    array.flags.writeable = False
    return array


@dataclasses.dataclass(frozen=True)
class Censoring:
    """A per-sample statement of which observations are limits, not detections.

    Censoring is *not* masking. ``results_schema.md`` §8 draws the line and
    hands the other side to this contract: a masked sample carries exactly zero
    information, whereas a 3sigma non-detection says something quite definite about
    the source. Discarding it throws that away; recording it as a datum with a
    large error bar is simply the wrong likelihood.

    A ``Censoring`` is declared *alongside* an observed container rather than
    inside it, because it is a statement about how the observation constrains
    the model — a likelihood concern — and the same container may legitimately
    be analysed with and without it.

    Parameters
    ----------
    kinds
        Integer array of :class:`LimitKind` codes, aligned index-by-index with
        the observed container's ``values``.

    Examples
    --------
    >>> import numpy as np
    >>> import astropy.units as u
    >>> from ampere.core import PhotometricPoints
    >>> phot = PhotometricPoints(
    ...     ["WISE_W3", "WISE_W4"],
    ...     [12.1, 22.2] * u.um,
    ...     [0.4, 0.9] * u.Jy,
    ...     uncertainty=[0.05, 0.30] * u.Jy,
    ...     extra_coords={"limit_kind": np.array([0, 1])},
    ... )
    >>> censoring = Censoring.from_extra_coord(phot, "limit_kind")
    >>> censoring.n_censored, censoring.any_censored
    (1, True)
    """

    kinds: np.ndarray

    def __post_init__(self) -> None:
        raw = np.asarray(self.kinds)
        if raw.dtype.kind == "b" or raw.dtype.kind not in "iu":
            raise LikelihoodError(
                f"Censoring's kinds must be an integer array of LimitKind codes "
                f"({[k.name for k in LimitKind]}), got dtype {raw.dtype}. A boolean array is "
                f"ambiguous — say which kind of limit it means with "
                f"Censoring.upper_limits(flags) or Censoring.lower_limits(flags)."
            )
        allowed = {int(kind) for kind in LimitKind}
        unknown = sorted(set(np.unique(raw).tolist()) - allowed)
        if unknown:
            raise LikelihoodError(
                f"Censoring was given unknown limit codes {unknown}. The declared codes are "
                f"{ {kind.name: int(kind) for kind in LimitKind} }."
            )
        object.__setattr__(self, "kinds", _readonly(np.array(raw, dtype=np.int8, copy=True)))

    # -- construction --------------------------------------------------------

    @classmethod
    def from_extra_coord(cls, container: FunctionSamples, name: str = "limit_kind") -> Censoring:
        """Read a censoring declaration off a container's ``extra_coords``.

        This is the hook ``results_schema.md`` §8 reserved: per-sample aligned
        arrays are validated against the container's shape exactly as values
        and mask are, "which is precisely the machinery a per-point
        ``limit_kind`` array needs".
        """
        if name not in container.extra_coords:
            available = sorted(container.extra_coords)
            raise LikelihoodError(
                f"{type(container).__name__} has no extra coordinate {name!r} to read a censoring "
                f"declaration from; it carries {available}. Attach one at construction, e.g. "
                f"extra_coords={{{name!r}: np.array([0, 1, 0])}} with LimitKind codes."
            )
        return cls(np.asarray(container.extra_coords[name]).ravel())

    @classmethod
    def upper_limits(cls, flags: Any) -> Censoring:
        """Build a declaration from a boolean "is an upper limit" array."""
        return cls._from_flags(flags, LimitKind.UPPER_LIMIT)

    @classmethod
    def lower_limits(cls, flags: Any) -> Censoring:
        """Build a declaration from a boolean "is a lower limit" array."""
        return cls._from_flags(flags, LimitKind.LOWER_LIMIT)

    @classmethod
    def _from_flags(cls, flags: Any, kind: LimitKind) -> Censoring:
        booleans = np.asarray(flags)
        if booleans.dtype.kind != "b":
            raise LikelihoodError(
                f"Censoring.{kind.name.lower()}s() expects a boolean array of flags, got dtype "
                f"{booleans.dtype}."
            )
        codes = np.where(booleans, int(kind), int(LimitKind.DETECTION))
        return cls(codes.ravel())

    # -- introspection -------------------------------------------------------

    @property
    def n_samples(self) -> int:
        """Number of samples this declaration covers."""
        return int(self.kinds.size)

    @property
    def n_censored(self) -> int:
        """Number of samples that are limits rather than detections."""
        return int(np.count_nonzero(self.kinds != int(LimitKind.DETECTION)))

    @property
    def any_censored(self) -> bool:
        """Whether this declaration contains any limit at all."""
        return self.n_censored > 0

    def check_against(self, container: FunctionSamples) -> None:
        """Verify alignment with the container this declaration describes."""
        if self.n_samples != container.n_samples:
            raise LikelihoodError(
                f"Censoring covers {self.n_samples} samples but the "
                f"{type(container).__name__} it was given has {container.n_samples}. A censoring "
                f"declaration is aligned index-by-index with the observed container."
            )

    def __repr__(self) -> str:
        counts = {kind.name: int(np.count_nonzero(self.kinds == int(kind))) for kind in LimitKind}
        body = ", ".join(f"{name}={n}" for name, n in counts.items() if n)
        return f"Censoring({body})"


# ---------------------------------------------------------------------------
# Kernels: neutral, declarative, hyperparameters are ordinary Parameters
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class KernelSpec:
    """The neutral, serialisable description of a kernel.

    Family name plus ordered hyperparameter names is the minimum W1.9's
    lowering table needs to emit a ``celerite2``/``tinygp``/GPyTorch term, and
    the maximum that translates across all three. It deliberately carries no
    values: the values are :class:`~ampere.core.parameter.Parameter`\\ s, which
    have their own declaration contract.
    """

    family: str
    hyperparameters: tuple[str, ...]
    quasiseparable: bool

    def to_dict(self) -> dict[str, Any]:
        """A plain-data form for provenance attrs and spec hashing."""
        return {
            "family": self.family,
            "hyperparameters": list(self.hyperparameters),
            "quasiseparable": self.quasiseparable,
        }


def _as_hyperparameter(
    name: str,
    given: Any,
    unit: u.UnitBase | None,
    *,
    positive: bool = True,
) -> Parameter:
    """Coerce a prior, a fixed number or a ready-made Parameter into a Parameter.

    GP hyperparameters are *ordinary parameters*, with ``Log`` bijections —
    ``parameters.md`` §13's instruction to this contract, discharged here in
    one place so no kernel can quietly do it differently.
    """
    if isinstance(given, Parameter):
        if given.name != name:
            raise LikelihoodError(
                f"kernel hyperparameter {name!r} was given a Parameter named {given.name!r}. "
                f"A kernel's hyperparameter names are part of its KernelSpec (W1.9 lowers them "
                f"to term keywords), so they are fixed; rename it with .rename({name!r})."
            )
        return given
    if isinstance(given, (int, float, np.floating, np.integer)) and not isinstance(given, bool):
        return Parameter(name, value=float(given), fixed=True, unit=unit)
    if hasattr(given, "ppf"):
        bijection = Log() if positive else Identity()
        return Parameter(name, given, unit=unit, bijection=bijection)
    raise LikelihoodError(
        f"kernel hyperparameter {name!r} must be a frozen scipy.stats distribution (a prior), a "
        f"number (held fixed), or an ampere Parameter — got {type(given).__name__}."
    )


class Kernel(Parameterised, abc.ABC):
    """A stationary covariance function, declared neutrally.

    A kernel is a *declaration*: a family name (:attr:`FAMILY`) and its
    hyperparameters as ordinary :class:`~ampere.core.parameter.Parameter`\\ s.
    It knows how to build its own dense covariance matrix — that is what
    :class:`DenseGP` needs, and what the conformance suite compares every other
    solver against — but a solver is free to ignore ``matrix`` entirely and
    lower :meth:`spec` to a state-space term instead.

    Subclasses declare :attr:`FAMILY`, :attr:`HYPERPARAMETERS`,
    :attr:`QUASISEPARABLE`, and implement :meth:`_covariance`, which receives
    non-negative separations. A non-stationary kernel would override
    :meth:`matrix` instead; the machinery does not assume stationarity anywhere
    outside :meth:`matrix`'s default implementation.
    """

    #: Neutral family name; W1.9's lowering table is keyed on it.
    FAMILY: ClassVar[str] = ""
    #: Hyperparameter names, in declaration order.
    HYPERPARAMETERS: ClassVar[tuple[str, ...]] = ()
    #: Whether this kernel has an exact quasiseparable (celerite-class)
    #: representation, and so admits an exact O(N) solve on ordered 1D data.
    QUASISEPARABLE: ClassVar[bool] = False

    def spec(self) -> KernelSpec:
        """The neutral description a lowering rule consumes."""
        return KernelSpec(
            family=self.FAMILY,
            hyperparameters=self.HYPERPARAMETERS,
            quasiseparable=self.QUASISEPARABLE,
        )

    # -- evaluation ----------------------------------------------------------

    def resolve(self, values: Mapping[str, Any] | None = None) -> dict[str, Any]:
        """This kernel's own hyperparameter values, out of a wider mapping.

        A :class:`NoiseModel` hands round one flat mapping covering itself and
        its kernel; the kernel picks out its own names rather than requiring
        the caller to split them.
        """
        if values is None:
            return self.context(None)
        selected = {name: values[name] for name in self.parameters.names if name in values}
        return self.context(selected)

    @abc.abstractmethod
    def _covariance(self, separation: np.ndarray, values: Mapping[str, Any]) -> np.ndarray:
        """Covariance at non-negative separations, given resolved values."""

    def matrix(self, left: np.ndarray, right: np.ndarray, values: Mapping[str, Any]) -> np.ndarray:
        """Dense covariance between two coordinate sets.

        ``left`` and ``right`` are ``(n, d)`` and ``(m, d)`` float64 arrays;
        a bare 1-D array is read as a column of ``n`` one-dimensional points.
        The result is ``(n, m)``. Separation is Euclidean in the coordinate
        space, which is why :meth:`GPSolver.check_compatible` requires every
        coordinate axis to share one unit.
        """
        resolved = self.resolve(values)
        points = _as_points(left, "kernel coordinates")
        other = _as_points(right, "kernel coordinates", dimensions=points.shape[1])
        difference = points[:, None, :] - other[None, :, :]
        separation = np.sqrt(np.einsum("ijk,ijk->ij", difference, difference))
        return self._covariance(separation, resolved)

    def diagonal(self, coordinates: np.ndarray, values: Mapping[str, Any]) -> np.ndarray:
        """The prior variance at each coordinate; ``k(0)`` for a stationary kernel."""
        resolved = self.resolve(values)
        n = int(_as_points(coordinates, "kernel coordinates").shape[0])
        return self._covariance(np.zeros(n, dtype=DTYPE), resolved)

    def __repr__(self) -> str:
        declared = ", ".join(repr(self.parameters[name]) for name in self.parameters.names)
        return f"{type(self).__name__}({declared})"


class Matern32(Kernel):
    r"""Matérn-3/2: ampere's canonical flexible-likelihood kernel.

    .. math::
        k(r) = a^2 \left(1 + \frac{\sqrt{3}\,r}{\ell}\right)
               \exp\!\left(-\frac{\sqrt{3}\,r}{\ell}\right)

    ``DEVELOPMENT_PLAN.md`` §2 and §4.4 make this the default throughout,
    replacing legacy's hardcoded RBF, for two independent reasons. It
    represents structured residuals better than a squared exponential (a
    once-differentiable sample path, not an analytic one — real model
    deficiencies are not infinitely smooth); and it is **exactly
    quasiseparable**, expressible as a sum of celerite/SHO terms, which is what
    makes :class:`QuasisepGP` an *exact* O(N) solve rather than an
    approximation. ``prior_art.md`` lesson S1 records that Starfish (Czekala et
    al. 2015) independently arrived at exactly this kernel, in velocity
    separation, for exactly this purpose.

    ``amplitude`` is the marginal **standard deviation** — ``k(0) ==
    amplitude²`` — following celerite2's ``Matern32Term(sigma=..., rho=...)``
    convention, so a prior on it is a prior in the data's own units.

    Parameters
    ----------
    amplitude, length_scale
        A frozen ``scipy.stats`` prior (fitted, with a :class:`Log` bijection),
        a number (held fixed), or a ready-made ``Parameter``.
    amplitude_unit
        Unit of ``amplitude``; must match the observed values' unit.
    length_scale_unit
        Unit of ``length_scale``; must match the coordinate axis's unit.

    Examples
    --------
    >>> import scipy.stats as st
    >>> import astropy.units as u
    >>> kernel = Matern32(st.loguniform(1e-3, 1e1), st.loguniform(1e-2, 1e2))
    >>> kernel.spec()
    KernelSpec(family='matern32', hyperparameters=('amplitude', 'length_scale'),
               quasiseparable=True)
    >>> kernel.parameters.free_names
    ('amplitude', 'length_scale')
    >>> kernel.parameters.bijections()
    (Log(lower=0.0), Log(lower=0.0))
    >>> float(kernel.matrix([[0.0]], [[0.0]], {"amplitude": 2.0, "length_scale": 1.0})[0, 0])
    4.0
    """

    FAMILY: ClassVar[str] = "matern32"
    HYPERPARAMETERS: ClassVar[tuple[str, ...]] = ("amplitude", "length_scale")
    QUASISEPARABLE: ClassVar[bool] = True

    def __init__(
        self,
        amplitude: Any,
        length_scale: Any,
        *,
        amplitude_unit: u.UnitBase | None = None,
        length_scale_unit: u.UnitBase | None = None,
    ) -> None:
        self.register_parameters(
            _as_hyperparameter("amplitude", amplitude, amplitude_unit),
            _as_hyperparameter("length_scale", length_scale, length_scale_unit),
        )

    def _covariance(self, separation: np.ndarray, values: Mapping[str, Any]) -> np.ndarray:
        amplitude = _positive(values["amplitude"], "amplitude", self.FAMILY, allow_zero=True)
        length_scale = _positive(values["length_scale"], "length_scale", self.FAMILY)
        scaled = _SQRT3 * separation / length_scale
        return amplitude * amplitude * (1.0 + scaled) * np.exp(-scaled)


class SquaredExponential(Kernel):
    r"""Squared exponential (RBF): legacy's kernel, kept for comparison.

    .. math::
        k(r) = a^2 \exp\!\left(-\frac{r^2}{2\ell^2}\right)

    Provided so the misspecification study in milestone M2 can compare the new
    default against what legacy ampere actually did, and so a user who wants it
    can have it. It is **not quasiseparable**: :class:`QuasisepGP` will refuse
    it, and it therefore does not scale past :class:`DenseGP`. That asymmetry
    is the concrete reason ``DEVELOPMENT_PLAN.md`` §2 made Matérn the default.
    """

    FAMILY: ClassVar[str] = "squared_exponential"
    HYPERPARAMETERS: ClassVar[tuple[str, ...]] = ("amplitude", "length_scale")
    QUASISEPARABLE: ClassVar[bool] = False

    def __init__(
        self,
        amplitude: Any,
        length_scale: Any,
        *,
        amplitude_unit: u.UnitBase | None = None,
        length_scale_unit: u.UnitBase | None = None,
    ) -> None:
        self.register_parameters(
            _as_hyperparameter("amplitude", amplitude, amplitude_unit),
            _as_hyperparameter("length_scale", length_scale, length_scale_unit),
        )

    def _covariance(self, separation: np.ndarray, values: Mapping[str, Any]) -> np.ndarray:
        amplitude = _positive(values["amplitude"], "amplitude", self.FAMILY, allow_zero=True)
        length_scale = _positive(values["length_scale"], "length_scale", self.FAMILY)
        return amplitude * amplitude * np.exp(-0.5 * (separation / length_scale) ** 2)


def _positive(value: Any, name: str, owner: str, *, allow_zero: bool = False) -> float:
    number = float(np.asarray(value, dtype=DTYPE))
    if not math.isfinite(number) or number < 0.0 or (number == 0.0 and not allow_zero):
        bound = ">= 0" if allow_zero else "> 0"
        raise LikelihoodError(
            f"{owner}'s {name!r} must be finite and {bound}, got {number!r}. Declare it with a "
            f"prior supported on the positive half-line and a Log bijection so no sampler can "
            f"propose a value outside it."
        )
    return number


# ---------------------------------------------------------------------------
# GP solver strategies
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class GPConditional:
    """The GP conditioned on the residuals: mean and variance at each location.

    This is the object ``DEVELOPMENT_PLAN.md`` §4.8 and the diagnostics spec's
    family C consume — "the conditioned GP mean ... already localises where the
    model is deficient". The mean is **signed**, so it shows the direction of
    the local deficiency and not merely its size.
    """

    mean: np.ndarray
    variance: np.ndarray

    @property
    def standard_deviation(self) -> np.ndarray:
        """Pointwise 1sigma band on :attr:`mean`."""
        return np.sqrt(np.clip(self.variance, 0.0, None))


class GPSolver(abc.ABC):
    """How the GP algebra is done. A strategy, chosen per problem.

    ``DEVELOPMENT_PLAN.md`` §4.4 makes the solver swappable behind the noise
    model precisely so the scaling story can change without the science code
    changing. Every strategy computes the same quantity — the Gaussian marginal
    log-likelihood of the residuals under ``K(θ) + diag(sigma²)`` — and
    :class:`DenseGP` is the definition of the right answer.

    A concrete strategy declares its own applicability through the class
    attributes below and enforces it in :meth:`check_compatible`, at
    composition time, loudly. It may branch on ``Axis.regular`` or
    ``Axis.log_regular`` for a fast path but must have a path when both are
    false (``results_schema.md`` §16).
    """

    #: Neutral strategy name, for provenance and error messages.
    NAME: ClassVar[str] = ""
    #: Whether the strategy computes the exact marginal likelihood of the
    #: declared kernel, as opposed to an approximation of it.
    EXACT: ClassVar[bool] = True
    #: Whether it needs a single, ordered, one-dimensional coordinate axis.
    REQUIRES_ORDERED_1D: ClassVar[bool] = False
    #: Whether it needs the kernel to have a quasiseparable representation.
    REQUIRES_QUASISEPARABLE: ClassVar[bool] = False
    #: Whether an implementation exists in the reference (numpy) path.
    IMPLEMENTED: ClassVar[bool] = False

    def check_compatible(self, kernel: Kernel, observed: FunctionSamples) -> None:
        """Composition-time check that this strategy can run this problem."""
        kind = type(observed).__name__
        # Declarative incompatibilities first: they are permanent facts about
        # the choice, whereas "not implemented yet" is temporary, and a user
        # who paired QuasisepGP with a non-quasiseparable kernel needs to hear
        # about the kernel rather than about Phase 2's schedule.
        if observed.LAYOUT is not Layout.POINTS:
            raise LikelihoodError(
                f"{self.NAME} was given a {kind}, whose layout is {observed.LAYOUT.value}. The "
                f"v1 GP solvers work on point-set containers (Spectrum, TimeSeries, "
                f"PhotometricPoints, VisibilitySet); gridded 2D+ data are the subject of the "
                f"SVGP / SKI / Vecchia strategy slots (DEVELOPMENT_PLAN.md §4.4, Phase 5)."
            )
        units = {axis.unit for axis in observed.axes}
        if len(units) > 1:
            named = sorted(str(unit) for unit in units)
            raise LikelihoodError(
                f"{self.NAME} measures separation as a Euclidean distance across a "
                f"{kind}'s coordinate axes, but they carry different units {named}. A single "
                f"isotropic length-scale is meaningless across mixed units; use one axis, or "
                f"declare a kernel that takes a length-scale per axis."
            )
        if self.REQUIRES_QUASISEPARABLE and not kernel.QUASISEPARABLE:
            raise LikelihoodError(
                f"{self.NAME} needs a kernel with an exact quasiseparable representation, but "
                f"{type(kernel).__name__} ({kernel.FAMILY}) has none. Use Matern32 — which is why "
                f"DEVELOPMENT_PLAN.md §2 made it the default — or switch to DenseGP."
            )
        if self.REQUIRES_ORDERED_1D and len(observed.axes) != 1:
            raise LikelihoodError(
                f"{self.NAME} needs one ordered coordinate axis, but a {kind} has "
                f"{len(observed.axes)}: {[axis.name for axis in observed.axes]}."
            )
        if not self.IMPLEMENTED:
            raise LikelihoodError(self._unimplemented_message())

    def _unimplemented_message(self) -> str:
        return (
            f"{self.NAME} is a declared strategy slot with no implementation yet "
            f"(DEVELOPMENT_PLAN.md §4.4). Its interface is fixed so backends can fill it in; "
            f"until then use DenseGP, which computes the same quantity exactly."
        )

    @abc.abstractmethod
    def log_marginal_likelihood(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
    ) -> float:
        """log N(residual; 0, K(values) + diag(variance))."""

    @abc.abstractmethod
    def condition(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
        at: np.ndarray | None = None,
    ) -> GPConditional:
        """The GP posterior at ``at`` (default: the data coordinates)."""

    @abc.abstractmethod
    def latent_transform(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        whitened: np.ndarray,
        values: Mapping[str, Any],
        *,
        jitter: float = 1e-10,
    ) -> np.ndarray:
        """Map whitened latent draws ``z`` to a GP draw ``f = L(θ) z``.

        This is the deterministic half of the latent-GP declaration (see
        :func:`latent_parameter`): the *prior* stays i.i.d. standard normal and
        all the correlation lives here, where the solver can impose it in
        whatever representation it uses.
        """

    def __repr__(self) -> str:
        return f"{type(self).__name__}()"


@dataclasses.dataclass(frozen=True)
class DenseGP(GPSolver):
    """Exact O(N³) dense Cholesky. The correctness anchor for every other solver.

    ``DEVELOPMENT_PLAN.md`` §4.4 calls this "current behaviour, any kernel,
    O(N³); correctness reference", and §4.6 makes "DenseGP↔QuasisepGP agreement
    on quasiseparable kernels" a conformance row. It has no performance goal
    whatsoever: it exists so that "what should the answer be?" has an
    implementation, not an argument.

    Parameters
    ----------
    jitter
        A standard deviation, in the data's own units, added in quadrature to
        the diagonal. Defaults to zero: a covariance that will not factorise is
        a fact about the model, and this contract does not hide it behind a
        silent numerical fudge. The error message tells you to set this when
        setting it is the right answer.
    """

    jitter: float = 0.0

    NAME: ClassVar[str] = "DenseGP"
    EXACT: ClassVar[bool] = True
    IMPLEMENTED: ClassVar[bool] = True

    def __post_init__(self) -> None:
        if not math.isfinite(self.jitter) or self.jitter < 0.0:
            raise LikelihoodError(f"DenseGP's jitter must be finite and >= 0, got {self.jitter!r}.")

    def _factor(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
    ) -> tuple[np.ndarray, bool]:
        covariance = kernel.matrix(coordinates, coordinates, values)
        total = covariance + np.diag(variance + self.jitter**2)
        if not np.all(np.isfinite(total)):
            raise LikelihoodError(
                "the covariance matrix K + diag(sigma^2) contains non-finite entries, so it "
                "cannot be factorised. The usual cause is a kernel amplitude large enough that "
                "amplitude**2 overflows float64; constrain the amplitude prior to the data's "
                "own scale."
            )
        try:
            return scipy.linalg.cho_factor(total, lower=True)
        except scipy.linalg.LinAlgError as error:
            raise LikelihoodError(
                f"the covariance matrix K + diag(sigma^2) is not positive definite, so its "
                f"Cholesky factorisation failed ({error}). Usual causes: a zero or near-duplicate "
                f"observational uncertainty, coordinates that are closer together than float64 "
                f"can separate at this length-scale, or a kernel amplitude far above the data "
                f"scale. Pass DenseGP(jitter=...) — a standard deviation in the data's units — "
                f"if the matrix is merely ill-conditioned rather than wrong."
            ) from error

    def log_marginal_likelihood(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
    ) -> float:
        factor = self._factor(kernel, coordinates, variance, values)
        alpha = scipy.linalg.cho_solve(factor, residual)
        log_determinant = 2.0 * float(np.sum(np.log(np.abs(np.diag(factor[0])))))
        quadratic = float(residual @ alpha)
        return -0.5 * (quadratic + log_determinant + residual.size * _LOG_2PI)

    def condition(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
        at: np.ndarray | None = None,
    ) -> GPConditional:
        points = _as_points(coordinates, "data coordinates")
        factor = self._factor(kernel, points, variance, values)
        if at is None:
            target = points
        else:
            target = _as_points(at, "conditioning grid", dimensions=points.shape[1])
        cross = kernel.matrix(target, points, values)
        mean = cross @ scipy.linalg.cho_solve(factor, residual)
        solved = scipy.linalg.cho_solve(factor, cross.T)
        prior_variance = kernel.diagonal(target, values)
        posterior = prior_variance - np.einsum("ij,ji->i", cross, solved)
        return GPConditional(mean=mean, variance=posterior)

    def latent_transform(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        whitened: np.ndarray,
        values: Mapping[str, Any],
        *,
        jitter: float = 1e-10,
    ) -> np.ndarray:
        points = _as_points(coordinates, "data coordinates")
        covariance = kernel.matrix(points, points, values)
        if not np.all(np.isfinite(covariance)):
            raise LikelihoodError(
                "the kernel matrix K contains non-finite entries, so the whitening transform "
                "f = L z is undefined. The usual cause is a kernel amplitude large enough that "
                "amplitude**2 overflows float64."
            )
        scale = float(np.mean(np.diag(covariance))) or 1.0
        stabilised = covariance + np.eye(covariance.shape[0]) * (jitter * scale)
        try:
            lower = scipy.linalg.cholesky(stabilised, lower=True)
        except scipy.linalg.LinAlgError as error:
            raise LikelihoodError(
                f"the kernel matrix K is not positive definite, so the whitening transform "
                f"f = L z is undefined ({error}). Increase the jitter argument, or check the "
                f"kernel hyperparameters."
            ) from error
        return lower @ _as_float64(whitened, "whitened latent draws")


class _SolverSlot(GPSolver):
    """Base for the declared-but-unimplemented strategies."""

    IMPLEMENTED: ClassVar[bool] = False

    def log_marginal_likelihood(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
    ) -> float:
        raise LikelihoodError(self._unimplemented_message())

    def condition(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
        at: np.ndarray | None = None,
    ) -> GPConditional:
        raise LikelihoodError(self._unimplemented_message())

    def latent_transform(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        whitened: np.ndarray,
        values: Mapping[str, Any],
        *,
        jitter: float = 1e-10,
    ) -> np.ndarray:
        raise LikelihoodError(self._unimplemented_message())


class QuasisepGP(_SolverSlot):
    """Exact O(N) for ordered 1D data via a quasiseparable (celerite-class) solve.

    The strategy ``DEVELOPMENT_PLAN.md`` §4.4 names as the scaling answer:
    celerite2 (numpy and jax), ``tinygp.solvers.QuasisepSolver``, GPyTorch or
    celerite2-torch on the torch side. It is **exact**, not approximate: a
    Matérn-3/2 kernel has an exact representation as a sum of celerite/SHO
    terms, so the state-space recursion computes the same marginal likelihood
    :class:`DenseGP` does, in linear time. That equivalence is a conformance
    row (§4.6), which is why :class:`DenseGP` exists at all.

    Implementation lands in Phase 2, in the backends; the interface is fixed
    here.
    """

    NAME: ClassVar[str] = "QuasisepGP"
    EXACT: ClassVar[bool] = True
    REQUIRES_ORDERED_1D: ClassVar[bool] = True
    REQUIRES_QUASISEPARABLE: ClassVar[bool] = True


class WindowedSparseGP(_SolverSlot):
    """Approximate O(N) by tapering the kernel to zero beyond a cutoff radius.

    Starfish's actual strategy (``prior_art.md`` lesson S2): apply a Hann
    window to the kernel so it is *exactly* zero beyond ``r₀ ≈ 4l``, giving a
    genuinely sparse banded covariance that a sparse Cholesky factorises in
    roughly linear time. This is a qualitatively different family from
    :class:`QuasisepGP` — approximate but simple, needing no state-space
    machinery, and applicable to kernels with no quasiseparable form at all
    (including :class:`SquaredExponential`) and to more than one dimension.

    It is a named slot rather than the default because :class:`QuasisepGP` is
    exact for the kernel ampere actually recommends, and buys that exactness
    without the cutoff radius — one more approximation-controlling
    hyperparameter that nothing chooses for the user. See the likelihoods
    contract for the full argument.
    """

    NAME: ClassVar[str] = "WindowedSparseGP"
    EXACT: ClassVar[bool] = False


class InducingPointGP(_SolverSlot):
    """Sparse variational / inducing-point GP (SVGP): a Phase 5 slot for 2D+."""

    NAME: ClassVar[str] = "InducingPointGP"
    EXACT: ClassVar[bool] = False


class StructuredGridGP(_SolverSlot):
    """SKI / KISS-GP: structured-kernel interpolation. A Phase 5 slot for 2D+."""

    NAME: ClassVar[str] = "StructuredGridGP"
    EXACT: ClassVar[bool] = False


class VecchiaGP(_SolverSlot):
    """Vecchia / nearest-neighbour approximation. A Phase 5 slot for 2D+."""

    NAME: ClassVar[str] = "VecchiaGP"
    EXACT: ClassVar[bool] = False


# ---------------------------------------------------------------------------
# Noise models
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class NoiseParams:
    """What a :class:`NoiseModel` hands a family for one evaluation.

    Every array here covers the **retained** samples only: masking has already
    been applied by excision (see :meth:`Likelihood.log_prob`), so a family
    never sees a mask and cannot forget to honour one.
    """

    #: Per-sample standard deviation, float64. ``None`` when the family
    #: supplies its own dispersion (Poisson).
    sigma: np.ndarray | None
    #: Resolved values of every parameter the likelihood declares.
    values: Mapping[str, Any]
    #: Coordinates, ``(n, d)`` float64, when a correlated solve is needed.
    coordinates: np.ndarray | None = None
    kernel: Kernel | None = None
    solver: GPSolver | None = None
    #: Latent function values, when the combination is
    #: :attr:`Marginalisation.LATENT` and inference supplied them.
    latent: np.ndarray | None = None
    #: :class:`LimitKind` codes for the retained samples, or ``None``.
    limits: np.ndarray | None = None
    #: Boolean inclusion indicator over the *full* containers (ruled
    #: 2026-09-02, W1.11 gap X-2), so a family carrying its own aligned
    #: per-sample data — a background spectrum, an instrumental template — can
    #: excise it the same way every array already in this record was excised.
    retain: np.ndarray | None = None

    @property
    def correlated(self) -> bool:
        """Whether a covariance beyond the diagonal is in play."""
        return self.kernel is not None

    @property
    def variance(self) -> np.ndarray:
        """``sigma**2``, or a loud failure if this noise model has no sigma."""
        if self.sigma is None:
            raise LikelihoodError(
                "the noise model supplied no per-sample sigma, because the observed container "
                "carries no uncertainties, so the diagonal of K + diag(sigma^2) is undefined. "
                "Attach uncertainties to the observed container, or give the noise model a "
                "'jitter' parameter. (A count-based family such as Poisson defines its own "
                "dispersion and never reaches this path.) "
                "Likelihood.check_alignment() catches this at composition time."
            )
        return self.sigma**2


class NoiseModel(Parameterised, abc.ABC):
    """What the noise is. A :class:`Parameterised`, so its knobs are parameters.

    ``parameters.md`` §13 instructs this contract that GP hyperparameters
    "are ordinary parameters on a ``Parameterised`` noise model, with ``Log``
    bijections", and that is what this class is: nothing about noise-model
    parameters is special, so tying, fixing, priors, plates, serialisation and
    W1.9's lowering all work on them unchanged.
    """

    #: Whether this model induces correlations between samples.
    CORRELATED: ClassVar[bool] = False

    @abc.abstractmethod
    def sigma(
        self,
        observed: FunctionSamples,
        retain: np.ndarray,
        values: Mapping[str, Any],
    ) -> np.ndarray | None:
        """Per-sample standard deviation on the retained samples."""

    def check_compatible(self, family: LikelihoodFamily, observed: FunctionSamples) -> None:
        """Composition-time check. Subclasses extend; this checks uncertainties."""
        if family.REQUIRES_UNCERTAINTY and observed.uncertainty is None:
            raise LikelihoodError(
                f"the {family.NAME} family needs per-sample uncertainties, but the "
                f"{type(observed).__name__} it was given has none. Attach them at construction "
                f"(uncertainty=...), or use a family that models its own dispersion."
            )

    def noise_params(
        self,
        observed: FunctionSamples,
        retain: np.ndarray,
        values: Mapping[str, Any],
        *,
        coordinates: np.ndarray | None = None,
        latent: np.ndarray | None = None,
        limits: np.ndarray | None = None,
    ) -> NoiseParams:
        """Assemble the per-evaluation noise description.

        The base implementation describes *uncorrelated* noise. A subclass that
        sets :attr:`CORRELATED` must override it, and is told so rather than
        being allowed to inherit a description that contradicts its own
        declaration — the resulting ``NoiseParams.correlated`` would be
        ``False`` and every family would take the uncorrelated branch.
        """
        if self.CORRELATED:
            raise LikelihoodError(
                f"{type(self).__name__} declares CORRELATED = True but inherits NoiseModel's "
                f"uncorrelated noise_params(), which reports no kernel and no solver. Every "
                f"family would then take its uncorrelated branch and the correlations would "
                f"silently do nothing. Override noise_params() to supply the kernel, solver and "
                f"coordinates (see GaussianProcessNoise)."
            )
        return NoiseParams(
            sigma=self.sigma(observed, retain, values),
            values=values,
            coordinates=None,
            kernel=None,
            solver=None,
            latent=latent,
            limits=limits,
            retain=retain,
        )

    def __repr__(self) -> str:
        declared = ", ".join(repr(self.parameters[name]) for name in self.parameters.names)
        return f"{type(self).__name__}({declared})"


def _observed_sigma(observed: FunctionSamples, retain: np.ndarray, owner: str) -> np.ndarray:
    if observed.uncertainty is None:
        raise LikelihoodError(
            f"{owner} needs the observed {type(observed).__name__}'s per-sample uncertainties, "
            f"but it has none. Attach them at construction (uncertainty=...)."
        )
    sigma = _as_float64(np.asarray(observed.uncertainty).ravel()[retain], "observed uncertainties")
    if np.any(sigma <= 0.0):
        raise LikelihoodError(
            f"{owner} was given zero or negative uncertainties on {int(np.sum(sigma <= 0.0))} "
            f"retained sample(s). A zero uncertainty is an infinitely precise measurement, which "
            f"no likelihood can normalise; mask the sample, or give it a real error bar."
        )
    return sigma


class IndependentNoise(NoiseModel):
    """Uncorrelated per-sample noise: the container's own uncertainties.

    The simple path, and it must stay simple —
    ``IndependentNoise()`` declares no parameters at all and does nothing but
    read ``observed.uncertainty``.

    Two optional knobs cover the cases a real fit needs. ``scale`` multiplies
    every uncertainty (the classic "the catalogue's error bars are
    underestimated by a factor" nuisance parameter); ``jitter`` adds a
    standard deviation in quadrature (an unmodelled noise floor). Both are
    ordinary parameters, so either may be fitted, fixed or tied.

    Examples
    --------
    >>> import numpy as np
    >>> import astropy.units as u
    >>> import scipy.stats as st
    >>> from ampere.core import Spectrum
    >>> noise = IndependentNoise()
    >>> noise.parameters.free_size
    0
    >>> inflated = IndependentNoise(scale=st.loguniform(0.5, 5.0))
    >>> inflated.parameters.free_names
    ('scale',)
    """

    CORRELATED: ClassVar[bool] = False

    def __init__(self, *, scale: Any = None, jitter: Any = None) -> None:
        if scale is not None:
            self.register_parameter(_as_hyperparameter("scale", scale, None))
        if jitter is not None:
            self.register_parameter(_as_hyperparameter("jitter", jitter, None))

    def sigma(
        self,
        observed: FunctionSamples,
        retain: np.ndarray,
        values: Mapping[str, Any],
    ) -> np.ndarray | None:
        resolved = self.context({k: v for k, v in values.items() if k in self.parameters})
        if observed.uncertainty is None:
            if "jitter" not in resolved:
                return None
            floor = _positive(resolved["jitter"], "jitter", "IndependentNoise")
            return np.full(int(np.count_nonzero(retain)), floor, dtype=DTYPE)
        sigma = _observed_sigma(observed, retain, "IndependentNoise")
        if "scale" in resolved:
            sigma = sigma * _positive(resolved["scale"], "scale", "IndependentNoise")
        if "jitter" in resolved:
            floor = _positive(resolved["jitter"], "jitter", "IndependentNoise", allow_zero=True)
            sigma = np.sqrt(sigma**2 + floor**2)
        return sigma


class GaussianProcessNoise(NoiseModel):
    """Ampere's flexible likelihood: a GP over the residuals.

    The distinguishing feature of the package. The model's residuals are given
    a covariance ``K(θ) + diag(sigma²)``; structure the physical model cannot
    explain is absorbed by the GP instead of biasing the physical parameters
    (``prior_art.md`` §4: Starfish is the direct scientific ancestor of this
    idea, and ampere generalises and scales it).

    The kernel's hyperparameters become this object's parameters, flatly: the
    ``Parameter`` objects are the kernel's own, shared by name and not by
    identity (``prior_art.md`` Tension 1 — object-identity sharing is the
    gammapy footgun this contract avoids everywhere).

    Parameters
    ----------
    kernel
        The covariance function. :class:`Matern32` unless you have a reason.
    solver
        The strategy that does the algebra. :class:`DenseGP` by default,
        because it is the one that is right rather than the one that is fast.
    scale, jitter
        As :class:`IndependentNoise`: applied to the *diagonal* term before
        the kernel is added.

    Examples
    --------
    >>> import scipy.stats as st
    >>> noise = GaussianProcessNoise(Matern32(st.loguniform(1e-3, 1e1), st.loguniform(0.1, 10.0)))
    >>> noise.parameters.free_names
    ('amplitude', 'length_scale')
    >>> noise.solver
    DenseGP(jitter=0.0)
    >>> noise.CORRELATED
    True
    """

    CORRELATED: ClassVar[bool] = True

    def __init__(
        self,
        kernel: Kernel,
        solver: GPSolver | None = None,
        *,
        scale: Any = None,
        jitter: Any = None,
    ) -> None:
        if not isinstance(kernel, Kernel):
            raise LikelihoodError(
                f"GaussianProcessNoise needs a Kernel, got {type(kernel).__name__}. Kernels are "
                f"declared neutrally (family name plus Parameter hyperparameters) so that W1.9 "
                f"can lower them to celerite2 / tinygp / GPyTorch terms."
            )
        self._kernel = kernel
        self._solver = DenseGP() if solver is None else solver
        for parameter in kernel.parameters:
            self.register_parameter(parameter)
        if scale is not None:
            self.register_parameter(_as_hyperparameter("scale", scale, None))
        if jitter is not None:
            self.register_parameter(_as_hyperparameter("jitter", jitter, None))

    @property
    def kernel(self) -> Kernel:
        """The covariance function this noise model applies."""
        return self._kernel

    @property
    def solver(self) -> GPSolver:
        """The strategy that performs the GP algebra."""
        return self._solver

    def sigma(
        self,
        observed: FunctionSamples,
        retain: np.ndarray,
        values: Mapping[str, Any],
    ) -> np.ndarray | None:
        resolved = self.context({k: v for k, v in values.items() if k in self.parameters})
        if observed.uncertainty is None:
            # A latent-GP combination (Poisson and friends) puts the GP on the
            # latent function, not on a residual, so there is no diagonal noise
            # term at all. Families that *do* need one say so through
            # REQUIRES_UNCERTAINTY and are refused at composition time.
            return None
        sigma = _observed_sigma(observed, retain, "GaussianProcessNoise")
        if "scale" in resolved:
            sigma = sigma * _positive(resolved["scale"], "scale", "GaussianProcessNoise")
        if "jitter" in resolved:
            floor = _positive(resolved["jitter"], "jitter", "GaussianProcessNoise", allow_zero=True)
            sigma = np.sqrt(sigma**2 + floor**2)
        return sigma

    def check_compatible(self, family: LikelihoodFamily, observed: FunctionSamples) -> None:
        super().check_compatible(family, observed)
        self._solver.check_compatible(self._kernel, observed)
        self._check_hyperparameter_units(observed)

    def _check_hyperparameter_units(self, observed: FunctionSamples) -> None:
        axis_unit = observed.axes[0].unit if observed.axes else None
        declared: Mapping[str, u.UnitBase | None] = {
            name: self.parameters[name].unit for name in self.parameters.names
        }
        length_scale = declared.get("length_scale")
        if length_scale is not None and length_scale != axis_unit:
            raise LikelihoodError(
                f"the kernel's 'length_scale' is declared in {length_scale} but the "
                f"{type(observed).__name__}'s coordinate axis is in {axis_unit}. Priors are "
                f"numeric in the declared unit and rescaling a distribution correctly is "
                f"family-specific, so this contract requires an exact match rather than a "
                f"conversion (parameters.md §8 makes the same ruling for tying)."
            )
        amplitude = declared.get("amplitude")
        if amplitude is not None and amplitude != observed.unit:
            raise LikelihoodError(
                f"the kernel's 'amplitude' is declared in {amplitude} but the "
                f"{type(observed).__name__}'s values are in {observed.unit}. The amplitude is a "
                f"marginal standard deviation in the data's own units; declare it in "
                f"{observed.unit} or leave its unit unset."
            )

    def noise_params(
        self,
        observed: FunctionSamples,
        retain: np.ndarray,
        values: Mapping[str, Any],
        *,
        coordinates: np.ndarray | None = None,
        latent: np.ndarray | None = None,
        limits: np.ndarray | None = None,
    ) -> NoiseParams:
        return NoiseParams(
            sigma=self.sigma(observed, retain, values),
            values=values,
            coordinates=coordinates,
            kernel=self._kernel,
            solver=self._solver,
            latent=latent,
            limits=limits,
            retain=retain,
        )


# ---------------------------------------------------------------------------
# The latent-GP declaration
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class LatentDeclaration:
    """What a latent-GP formulation adds to a fitting problem's parameters.

    The latent values are declared **whitened**: ``z ~ Normal(0, 1)``
    independently, with the correlation supplied deterministically by
    :meth:`GPSolver.latent_transform` as ``f = L(θ) z``. See
    :func:`latent_parameter` for why this, and not a
    :class:`~ampere.core.parameter.HierarchicalPrior`, is the declaration.
    """

    parameter: Parameter
    size: int
    solver: GPSolver

    def as_parameter_set(self) -> ParameterSet:
        """The declaration as a set, ready for ``ParameterSet.merge``."""
        return ParameterSet([self.parameter])


def latent_parameter(name: str, size: int) -> Parameter:
    """Declare ``size`` whitened latent GP values as one array-valued parameter.

    ``parameters.md`` §13 asked this contract to answer whether an
    array-valued ``HierarchicalPrior`` scales as the latent-GP declaration for
    ``N ~ 10⁵``. The answer is that a ``HierarchicalPrior`` is not the right
    declaration at *any* N, for a reason that has nothing to do with size:
    ``parameters.md`` §12.1 states that an array-valued parameter's prior is
    i.i.d. across its elements, and a GP prior is precisely *not* i.i.d. — its
    whole content is the correlation between elements. Declaring ``f`` with a
    hierarchical Normal referencing the GP hyperparameters would silently
    describe white noise with a fitted variance.

    The declaration is therefore the standard non-centred (whitened)
    parameterisation, which needs no extension to ``parameters.md`` at all:

    * ``z`` is one array-valued parameter of shape ``(N,)`` with an i.i.d.
      standard-normal prior and an ``Identity`` bijection — exactly what
      §12.1 already supports, and exactly what a PPL wants for HMC geometry;
    * the covariance enters through ``f = L(θ) z``, a deterministic transform
      owned by the :class:`GPSolver`, where it can be a Cholesky factor
      (:class:`DenseGP`) or a state-space recursion (:class:`QuasisepGP`)
      without the declaration changing.

    On scale, ``N ~ 10⁵`` is fine as a *declaration* — one ``Parameter``, one
    shape tuple, one ``PriorSpec`` — and the flat-vector pack/unpack is a
    single reshape. Three things downstream are not fine, and the likelihoods
    contract records them as obligations rather than leaving them to be
    discovered: ``ParameterSet.free_labels()`` materialises one string per
    element (W1.8 must give ArviZ a dimension and a coordinate, not 10⁵ scalar
    names); the sampler dimension becomes ``N + k``, which no gradient-free
    engine can address (hence :meth:`Likelihood.check_engine`); and forming
    ``L`` densely is O(N³), so the latent path at that scale requires
    :class:`QuasisepGP`, not :class:`DenseGP`.
    """
    if not isinstance(size, (int, np.integer)) or int(size) < 1:
        raise LikelihoodError(
            f"a latent GP declaration needs a positive number of latent values, got {size!r}."
        )
    return Parameter(
        name,
        st.norm(0.0, 1.0),
        shape=(int(size),),
        bijection=Identity(),
        description=(
            "whitened latent GP values; the correlated draw is f = L(theta) @ z, applied by "
            "GPSolver.latent_transform"
        ),
    )


# ---------------------------------------------------------------------------
# Likelihood families and their registry
# ---------------------------------------------------------------------------

_FAMILIES: dict[str, type[LikelihoodFamily]] = {}


def register_family(cls: type[LikelihoodFamily]) -> type[LikelihoodFamily]:
    """Register a family under its :attr:`~LikelihoodFamily.NAME`.

    Usable as a decorator. ``DEVELOPMENT_PLAN.md`` §4.4 requires that "users
    can add their own without touching ampere", so this is public and the
    registry is an ordinary mutable mapping — a third-party family is a class
    with one method and a decorator.
    """
    name = cls.NAME
    if not name:
        raise LikelihoodError(
            f"{cls.__name__} must declare a non-empty NAME before it can be registered."
        )
    existing = _FAMILIES.get(name)
    if existing is not None and existing is not cls:
        raise LikelihoodError(
            f"a likelihood family named {name!r} is already registered "
            f"({existing.__module__}.{existing.__qualname__}). Family names are the key W1.8's "
            f"provenance records and W1.9's lowering table use; pick another."
        )
    _FAMILIES[name] = cls
    return cls


def family_named(name: str) -> type[LikelihoodFamily]:
    """The registered family class of that name."""
    try:
        return _FAMILIES[name]
    except KeyError:
        raise LikelihoodError(
            f"no likelihood family named {name!r} is registered; the registered families are "
            f"{list(list_families())}. Add your own with @register_family."
        ) from None


def list_families() -> tuple[str, ...]:
    """Every registered family name, in registration order."""
    return tuple(_FAMILIES)


class LikelihoodFamily(Parameterised, abc.ABC):
    """The sampling distribution of one datum, given a prediction and its noise.

    One method — :meth:`log_prob` — is the whole evaluation interface, with the
    signature ``DEVELOPMENT_PLAN.md`` §4.4 specifies. ``prior_art.md`` lesson
    3M3 is the reason it stays that small: 3ML's ``get_log_like()`` is the
    right *shape* of contract, and its opacity (the plugin also owning the
    instrument response) is the part not to copy.

    Everything else a family declares is a class attribute, so declaring costs
    nothing at evaluation time and composition can be checked before a sampler
    starts. A family may also be :class:`Parameterised` — Student-t's degrees
    of freedom are an ordinary fitted parameter.
    """

    #: Registry key. Also what appears in provenance and error messages.
    NAME: ClassVar[str] = ""
    #: Whether a GP covariance can be folded into this family's own noise
    #: process and marginalised in closed form. True for Gaussian noise only.
    ANALYTIC_WITH_GP: ClassVar[bool] = False
    #: Whether :meth:`log_prob` actually implements the latent-conditional form
    #: — i.e. whether it reads ``noise.latent`` and refuses to proceed without
    #: it. **False by default, deliberately**: a family that declares
    #: :attr:`Marginalisation.LATENT` but ignores the latent values would
    #: silently return the *uncorrelated* likelihood, so every GP
    #: hyperparameter and every latent value an engine sampled would leave the
    #: log-probability untouched — a fit that runs, converges, and is wrong.
    #: :class:`Likelihood` therefore refuses a latent combination whose family
    #: has not opted in, exactly as it refuses an unimplemented family.
    CONSUMES_LATENT_GP: ClassVar[bool] = False
    #: Whether this family can consume a :class:`Censoring` declaration.
    SUPPORTS_CENSORING: ClassVar[bool] = False
    #: Whether it needs per-sample uncertainties on the observed container.
    REQUIRES_UNCERTAINTY: ClassVar[bool] = True
    #: Whether it accepts complex-valued containers (VisibilitySet).
    ALLOWS_COMPLEX: ClassVar[bool] = False
    #: Whether an implementation exists, as opposed to a declared slot.
    IMPLEMENTED: ClassVar[bool] = True

    @abc.abstractmethod
    def log_prob(
        self,
        predicted: np.ndarray,
        observed: np.ndarray,
        noise: NoiseParams,
    ) -> float:
        """Total log-probability of ``observed`` given ``predicted`` and ``noise``.

        Both arrays cover the retained samples only, in float64 (or complex128
        for a complex family): masking has already been applied by excision.
        """

    def sample(
        self,
        predicted: np.ndarray,
        noise: NoiseParams,
        rng: np.random.Generator,
    ) -> np.ndarray:
        """One draw of the retained observed values, given *predicted* and *noise*.

        The generative counterpart of :meth:`log_prob` (ruled 2026-09-02,
        ``inference.md`` §19 R3): a draw from the same distribution
        ``log_prob`` scores, over the retained samples, using the same
        ``NoiseParams``. ``simulate(observe=True)`` (W1.7) is the consumer.

        The default **refuses, specifically**: ``log_prob`` defines how a
        datum is scored, not how one is generated, and ampere will not guess a
        sampling distribution — a wrong guess would silently train an SBI
        posterior on the wrong forward model. A family for which the
        observation process is well defined overrides this with it; a user
        family may do the same, which is the supported route to
        ``observe=True`` draws for an exotic observation process.
        """
        raise LikelihoodError(
            f"the {self.NAME or type(self).__name__} family does not implement sample(): its "
            f"log_prob defines how a datum is scored, not how one is generated, and ampere will "
            f"not guess a sampling distribution. Subclass {type(self).__name__} and override "
            f"sample(predicted, noise, rng) with the observation process, or use "
            f"simulate(observe=False) and draw observations from the predicted containers "
            f"yourself."
        )

    def check_observed(self, observed: FunctionSamples) -> None:
        """Composition-time precondition on the observed data. Default: none.

        The family's half of the obligation ``NoiseModel.check_compatible``
        already has (ruled 2026-09-02, W1.11 gap I-5): some preconditions are
        properties of the *sampling distribution* rather than of the noise — a
        circular family needs angles in radians, ``PoissonFamily`` needs
        integer counts, a Rice family needs non-negative amplitudes. Called by
        :meth:`Likelihood.check_alignment` on the **observed** container only:
        the unit check has already forced predicted and observed to agree on
        everything a container carries, and value-range properties genuinely
        differ between the two (a Poisson *rate* is not an integer). Checks
        that belong here are the ones otherwise forced into the hot loop or
        nowhere at all. Only retained samples should be held to a
        precondition — a masked sample carries zero information.
        """

    def marginalisation_with(
        self,
        noise: NoiseModel,
        censoring: Censoring | None = None,
    ) -> Marginalisation:
        """Whether this family plus that noise model closes analytically.

        The default rule implements ``DEVELOPMENT_PLAN.md`` §4.4 exactly: an
        uncorrelated noise model is always analytic; a correlated one is
        analytic only for a family whose own noise process is Gaussian; and
        censoring under a correlated noise model is latent regardless, because
        a censored multivariate Gaussian likelihood is an orthant probability
        with no closed form beyond a handful of dimensions.
        """
        if not noise.CORRELATED:
            return Marginalisation.ANALYTIC
        if censoring is not None and censoring.any_censored:
            return Marginalisation.LATENT
        return Marginalisation.ANALYTIC if self.ANALYTIC_WITH_GP else Marginalisation.LATENT

    def _unimplemented(self) -> LikelihoodError:
        return LikelihoodError(
            f"the {self.NAME} family is declared but not implemented "
            f"(DEVELOPMENT_PLAN.md §4.4 puts it in the interface design and stages the "
            f"implementation). Its declaration is live — list_families() reports it and "
            f"composition checks against it — but it cannot be evaluated yet."
        )

    def __repr__(self) -> str:
        declared = ", ".join(repr(self.parameters[name]) for name in self.parameters.names)
        return f"{type(self).__name__}({declared})"


def _independent_sigma(noise: NoiseParams, family: str) -> np.ndarray:
    if noise.sigma is None:
        raise LikelihoodError(
            f"the {family} family needs per-sample uncertainties, but the noise model supplied "
            f"none — the observed container has no 'uncertainty'. Attach them at construction "
            f"(uncertainty=...), or give the noise model a 'jitter' parameter to stand in for "
            f"them. Likelihood.check_alignment() catches this at composition time."
        )
    return noise.sigma


def _location_scale_log_prob(
    distribution: Any,
    standardised: np.ndarray,
    sigma: np.ndarray,
    limits: np.ndarray | None,
    *shape_args: float,
) -> float:
    """Sum a location-scale family's log-density, honouring a censoring declaration.

    This is the whole of the censored-likelihood mathematics, written once for
    every family whose standardised residual has a scipy distribution: a
    detection contributes ``log f(z) - log sigma``; an upper limit contributes
    ``log F(z)``, the probability that the truth lies below the recorded value;
    a lower limit contributes ``log (1 - F(z))``. The Tobit construction, and
    the reason issue #11's "upper limits are fairly straightforward in other
    settings" is true for uncorrelated noise.
    """
    if limits is None:
        return float(np.sum(distribution.logpdf(standardised, *shape_args) - np.log(sigma)))
    total = 0.0
    detection = limits == int(LimitKind.DETECTION)
    if np.any(detection):
        kept = standardised[detection]
        total += float(np.sum(distribution.logpdf(kept, *shape_args) - np.log(sigma[detection])))
    upper = limits == int(LimitKind.UPPER_LIMIT)
    if np.any(upper):
        total += float(np.sum(distribution.logcdf(standardised[upper], *shape_args)))
    lower = limits == int(LimitKind.LOWER_LIMIT)
    if np.any(lower):
        total += float(np.sum(distribution.logsf(standardised[lower], *shape_args)))
    return total


@register_family
class GaussianFamily(LikelihoodFamily):
    """The Gaussian family: the one whose GP marginalisation closes.

    With :class:`IndependentNoise` this is the ordinary chi-square-plus-
    normalisation log-likelihood. With :class:`GaussianProcessNoise` it is the
    GP marginal likelihood ``log N(residual; 0, K(θ) + diag(sigma²))`` — the
    flexible likelihood, delegated to the :class:`GPSolver` strategy.

    It is the only family here with ``ANALYTIC_WITH_GP = True``, and that is
    not an accident of implementation effort: the "add the GP covariance to the
    noise and integrate" trick is a property of Gaussian conjugacy, and every
    other family needs the latent path.

    Examples
    --------
    >>> import numpy as np
    >>> import astropy.units as u
    >>> import scipy.stats as st
    >>> from ampere.core import Spectrum
    >>> data = Spectrum([1.0, 2.0, 3.0] * u.um, [1.0, 2.0, 3.0] * u.Jy,
    ...                 uncertainty=[0.1, 0.1, 0.1] * u.Jy)
    >>> model = data.with_values([1.05, 2.0, 2.9])
    >>> like = Likelihood(GaussianFamily(), IndependentNoise())
    >>> round(like.log_prob(model, data), 6)
    3.52594
    >>> like.marginalisation
    <Marginalisation.ANALYTIC: 'analytic'>
    """

    NAME: ClassVar[str] = "gaussian"
    ANALYTIC_WITH_GP: ClassVar[bool] = True
    SUPPORTS_CENSORING: ClassVar[bool] = True

    def log_prob(
        self,
        predicted: np.ndarray,
        observed: np.ndarray,
        noise: NoiseParams,
    ) -> float:
        residual = observed - predicted
        if noise.correlated:
            if noise.limits is not None and np.any(noise.limits != int(LimitKind.DETECTION)):
                raise LikelihoodError(
                    "censored data under a correlated (GP) noise model do not marginalise "
                    "analytically: the likelihood is a multivariate-normal orthant probability, "
                    "which has no closed form beyond a few dimensions. This combination is "
                    "declared LATENT (see Likelihood.marginalisation) and needs a modern-backend "
                    "engine that can sample the truncated latent values."
                )
            assert noise.solver is not None and noise.kernel is not None  # narrowed by .correlated
            assert noise.coordinates is not None
            return noise.solver.log_marginal_likelihood(
                noise.kernel, noise.coordinates, residual, noise.variance, noise.values
            )
        sigma = _independent_sigma(noise, self.NAME)
        if noise.limits is not None:
            return _location_scale_log_prob(st.norm, residual / sigma, sigma, noise.limits)
        # The uncensored case is written out rather than delegated: it is the
        # hot loop of every ordinary fit ampere runs, and scipy's generic
        # machinery costs several times what the closed form does.
        return float(np.sum(-0.5 * ((residual / sigma) ** 2 + _LOG_2PI) - np.log(sigma)))

    def sample(
        self,
        predicted: np.ndarray,
        noise: NoiseParams,
        rng: np.random.Generator,
    ) -> np.ndarray:
        """A draw from the same distribution :meth:`log_prob` scores.

        * uncorrelated noise — ``x = mu + sigma z``, with the noise model's own
          ``sigma``, so a fitted ``scale`` or ``jitter`` is already in it;
        * correlated (GP) noise — ``x = mu + L z1 + sigma z2``, where ``L``
          comes from :meth:`GPSolver.latent_transform`, the same whitening the
          latent declaration uses. That is a draw from
          ``N(mu, K + diag(sigma^2))`` **up to the solver's numerical
          stabiliser**: the solver's own jitter is part of the covariance it
          scores, so it is folded into the draw here — omitting it would draw
          from a narrower distribution than the likelihood evaluates, and the
          error is not small at the jitter values the library's own error
          message tells a user to raise.
        """
        realisation = np.asarray(predicted, dtype=DTYPE).copy()
        sigma = None if noise.sigma is None else np.asarray(noise.sigma, dtype=DTYPE)
        if noise.correlated:
            assert noise.solver is not None and noise.kernel is not None  # narrowed by .correlated
            assert noise.coordinates is not None
            whitened = rng.standard_normal(realisation.shape)
            realisation = realisation + noise.solver.latent_transform(
                noise.kernel, noise.coordinates, whitened, noise.values
            )
            stabiliser = float(getattr(noise.solver, "jitter", 0.0) or 0.0)
            if stabiliser:
                floor = np.full(realisation.shape, stabiliser)
                sigma = floor if sigma is None else np.sqrt(sigma**2 + floor**2)
        if sigma is not None:
            realisation = realisation + sigma * rng.standard_normal(realisation.shape)
        return realisation


@register_family
class StudentTFamily(LikelihoodFamily):
    """Heavy-tailed robustness to outlying samples.

    ``DEVELOPMENT_PLAN.md`` §4.4 lists Student-t for outlier robustness. Note
    that this is robustness to *individual bad samples*, which is a different
    problem from the correlated-residual misspecification
    :class:`GaussianProcessNoise` addresses — the two compose, but only through
    the latent path, because a Student-t is a scale mixture of Gaussians and
    the mixture does not commute with a GP covariance in closed form.

    ``nu`` (degrees of freedom) is an ordinary parameter: fit it, fix it, or
    tie it across datasets.
    """

    NAME: ClassVar[str] = "student_t"
    ANALYTIC_WITH_GP: ClassVar[bool] = False
    SUPPORTS_CENSORING: ClassVar[bool] = True

    def __init__(self, nu: Any = 4.0) -> None:
        self.register_parameter(_as_hyperparameter("nu", nu, None))

    def log_prob(
        self,
        predicted: np.ndarray,
        observed: np.ndarray,
        noise: NoiseParams,
    ) -> float:
        sigma = _independent_sigma(noise, self.NAME)
        resolved = self.context({k: v for k, v in noise.values.items() if k in self.parameters})
        nu = _positive(resolved["nu"], "nu", self.NAME)
        standardised = (observed - predicted) / sigma
        return _location_scale_log_prob(st.t, standardised, sigma, noise.limits, nu)


@register_family
class CauchyFamily(LikelihoodFamily):
    """Student-t with ``nu = 1``: maximal outlier tolerance, no finite variance.

    Registered separately because it takes no parameters and because a user
    reaching for "the Cauchy likelihood" should find it by that name.
    """

    NAME: ClassVar[str] = "cauchy"
    ANALYTIC_WITH_GP: ClassVar[bool] = False
    SUPPORTS_CENSORING: ClassVar[bool] = True

    def log_prob(
        self,
        predicted: np.ndarray,
        observed: np.ndarray,
        noise: NoiseParams,
    ) -> float:
        sigma = _independent_sigma(noise, self.NAME)
        standardised = (observed - predicted) / sigma
        return _location_scale_log_prob(st.cauchy, standardised, sigma, noise.limits)


@register_family
class ComplexGaussianFamily(LikelihoodFamily):
    """Circular complex Gaussian: the standard interferometric noise model.

    ``results_schema.md`` §16 hands this contract the statement that a
    ``VisibilitySet``'s real-valued uncertainty "encodes the circular complex
    Gaussian only; a non-circular or amplitude/phase noise model must supply
    its own structure". This family is that model, and nothing more: the real
    and imaginary parts are independent ``Normal(0, sigma²)``, so

    ``log p = Σ [ -log(2π sigma²) - |y - μ|² / (2sigma²) ]``.

    A non-circular model (unequal real/imaginary variances, or a correlation
    between them) is a *different noise model* supplying a 2x2 structure per
    sample, and belongs to whoever needs it; the amplitude/phase formulations
    are :class:`RiceFamily` and :class:`VonMisesFamily`.

    :attr:`ANALYTIC_WITH_GP` is deliberately ``False`` even though a *particular*
    complex GP — one real kernel applied independently to the real and
    imaginary parts — would marginalise perfectly well in closed form. Which
    complex GP is the right model for correlated visibility noise (are the two
    components' covariances equal? is there a non-zero pseudo-covariance?) is a
    Phase-4 modelling question this contract will not settle unilaterally, and
    declaring ANALYTIC would amount to answering it. See the likelihoods
    contract's open questions.
    """

    NAME: ClassVar[str] = "complex_gaussian"
    ALLOWS_COMPLEX: ClassVar[bool] = True

    def log_prob(
        self,
        predicted: np.ndarray,
        observed: np.ndarray,
        noise: NoiseParams,
    ) -> float:
        sigma = _independent_sigma(noise, self.NAME)
        residual = np.abs(observed - predicted)
        variance = sigma**2
        return float(np.sum(-(residual**2) / (2.0 * variance) - _LOG_2PI - np.log(variance)))


@register_family
class PoissonFamily(LikelihoodFamily):
    """Counting statistics — and the family that exercises the latent path.

    With :class:`IndependentNoise` this is the plain Poisson log-pmf: the
    prediction is the expected count, the observation is an integer, and there
    is no sigma anywhere (:attr:`REQUIRES_UNCERTAINTY` is ``False``).

    With :class:`GaussianProcessNoise` it is the case
    ``DEVELOPMENT_PLAN.md`` §4.4 singles out. There is no covariance matrix to
    add: robustness needs a latent-GP formulation, ``counts ~ Poisson(rate ·
    exp(f))`` with ``f ~ GP``, marginalised numerically. So the combination
    declares :attr:`Marginalisation.LATENT`, :meth:`Likelihood.check_engine`
    refuses a gradient-free engine, and ``log_prob`` requires the latent values
    to be supplied — it is the *conditional* likelihood given ``f``, which is
    what an HMC/VI backend evaluates per gradient step.

    Examples
    --------
    >>> import numpy as np
    >>> import astropy.units as u
    >>> import scipy.stats as st
    >>> from ampere.core import Spectrum
    >>> counts = Spectrum([1.0, 2.0, 3.0] * u.um, [4.0, 7.0, 2.0])
    >>> rate = counts.with_values([3.5, 6.0, 2.5])
    >>> analytic = Likelihood(PoissonFamily(), IndependentNoise())
    >>> analytic.marginalisation
    <Marginalisation.ANALYTIC: 'analytic'>
    >>> round(analytic.log_prob(rate, counts), 6)
    -5.010413
    >>> flexible = Likelihood(PoissonFamily(), GaussianProcessNoise(Matern32(0.3, 1.0)))
    >>> flexible.marginalisation
    <Marginalisation.LATENT: 'latent'>
    >>> flexible.check_engine(differentiable=False, engine="emcee")
    Traceback (most recent call last):
        ...
    ampere.core.exceptions.LikelihoodError: emcee cannot run this likelihood: the poisson family
    with a GaussianProcessNoise noise model introduces N latent values ...
    """

    NAME: ClassVar[str] = "poisson"
    ANALYTIC_WITH_GP: ClassVar[bool] = False
    CONSUMES_LATENT_GP: ClassVar[bool] = True
    REQUIRES_UNCERTAINTY: ClassVar[bool] = False

    def check_observed(self, observed: FunctionSamples) -> None:
        """Integer counts, checked once at composition (ruled 2026-09-02, X-3).

        This ran in ``log_prob`` before the ``check_observed`` hook existed —
        an O(N) scan of the *data* on every evaluation, which §13's
        compile-once/evaluate-many split says is the wrong place. The
        ``rate > 0`` guard stays in the hot loop, because it is a property of
        the prediction and genuinely varies per draw.
        """
        counts = np.asarray(observed.values).ravel()
        kept = counts[np.asarray(observed.valid).ravel()]
        if np.any(kept < 0.0) or not np.all(kept == np.round(kept)):
            raise LikelihoodError(
                "the poisson family needs non-negative integer counts, but the observed values "
                "are not integral. Counts are counts; if the data are rates, multiply by the "
                "exposure in the instrument chain (W1.5) rather than here."
            )

    def log_prob(
        self,
        predicted: np.ndarray,
        observed: np.ndarray,
        noise: NoiseParams,
    ) -> float:
        counts = np.asarray(observed)
        rate = np.asarray(predicted, dtype=DTYPE)
        if noise.correlated:
            if noise.latent is None:
                raise LikelihoodError(
                    "a Poisson likelihood with a GaussianProcessNoise model is a latent-variable "
                    "model: it is only defined given the latent function f, as "
                    "counts ~ Poisson(rate * exp(f)). Supply latent=... (an inference engine "
                    "that samples f does), or use IndependentNoise for the analytic case."
                )

            latent = _as_float64(noise.latent, "latent GP values")
            if latent.shape != rate.shape:
                raise LikelihoodError(
                    f"the latent GP values have shape {latent.shape} but there are {rate.shape} "
                    f"retained samples. One latent value per retained sample."
                )
            rate = rate * np.exp(latent)
        if np.any(rate <= 0.0):
            raise LikelihoodError(
                "the poisson family needs a strictly positive expected count; the model "
                "predicted a value <= 0."
            )
        return float(np.sum(st.poisson.logpmf(counts, rate)))


@register_family
class RiceFamily(LikelihoodFamily):
    """Rician amplitude noise — polarised intensity, debiased visibility amplitudes.

    Declared, not implemented. The interface question this contract fixes is
    that a Rice family consumes the *same* per-sample sigma as the Gaussian one
    (it is the amplitude of a circular complex Gaussian), so a
    ``VisibilitySet`` needs no extra structure to support it. What is *not*
    fixed, and is an open question in the spec, is whether the model's
    prediction is the true amplitude or the underlying complex value.
    """

    NAME: ClassVar[str] = "rice"
    IMPLEMENTED: ClassVar[bool] = False

    def log_prob(
        self,
        predicted: np.ndarray,
        observed: np.ndarray,
        noise: NoiseParams,
    ) -> float:
        raise self._unimplemented()


@register_family
class VonMisesFamily(LikelihoodFamily):
    """Wrapped/von Mises phase noise — closure phases, position angles.

    Declared, not implemented. The interface question this contract fixes is
    that the observed and predicted values are angles in radians and the
    residual must be wrapped, not subtracted; the open question is the mapping
    from a per-sample sigma to a concentration κ.
    """

    NAME: ClassVar[str] = "von_mises"
    IMPLEMENTED: ClassVar[bool] = False

    def log_prob(
        self,
        predicted: np.ndarray,
        observed: np.ndarray,
        noise: NoiseParams,
    ) -> float:
        raise self._unimplemented()


# ---------------------------------------------------------------------------
# Likelihood: the composition
# ---------------------------------------------------------------------------


class Likelihood(Parameterised):
    """A family, a noise model and (optionally) a censoring declaration.

    This is the object a ``Dataset`` holds (W1.7) and the thing that turns
    "predicted vs observed containers" into a number. It owns three
    responsibilities the pieces below it deliberately do not:

    1. **Masking.** The effective mask is the union of the predicted and
       observed containers' masks — a sample is used only if *both* are valid —
       read through ``weights()`` and applied by **excision**. See
       :meth:`log_prob`.
    2. **Composition-time checking.** :meth:`check_alignment` is the O(N) pass
       that belongs once, not per evaluation.
    3. **The marginalisation declaration** and the engine capability check
       built on it.

    Its parameters are the family's and the noise model's, in one flat
    :class:`~ampere.core.parameter.ParameterSet` — not a nested merge, because
    ``parameters.md`` §12.4 records that ``merge`` is not associative and W1.7
    must merge every dataset in one call. A name collision between the family
    and the noise model raises here rather than being silently qualified.

    Examples
    --------
    >>> import numpy as np
    >>> import astropy.units as u
    >>> import scipy.stats as st
    >>> from ampere.core import Spectrum
    >>> flexible = Likelihood(
    ...     GaussianFamily(),
    ...     GaussianProcessNoise(Matern32(st.loguniform(1e-3, 1e1), st.loguniform(0.1, 10.0))),
    ... )
    >>> flexible.parameters.free_names
    ('amplitude', 'length_scale')
    >>> flexible.marginalisation
    <Marginalisation.ANALYTIC: 'analytic'>
    """

    def __init__(
        self,
        family: LikelihoodFamily,
        noise: NoiseModel | None = None,
        *,
        censoring: Censoring | None = None,
    ) -> None:
        if not isinstance(family, LikelihoodFamily):
            raise LikelihoodError(
                f"Likelihood needs a LikelihoodFamily, got {type(family).__name__}. Register a "
                f"custom family with @register_family; it needs one method, log_prob."
            )
        if not family.IMPLEMENTED:
            raise LikelihoodError(
                f"the {family.NAME} family is declared but not implemented, so it cannot be "
                f"composed into a Likelihood yet (DEVELOPMENT_PLAN.md §4.4 stages the "
                f"implementation). It is visible in list_families() so the target set is on "
                f"record."
            )
        noise = IndependentNoise() if noise is None else noise
        if not isinstance(noise, NoiseModel):
            raise LikelihoodError(f"Likelihood needs a NoiseModel, got {type(noise).__name__}.")
        if censoring is not None and not family.SUPPORTS_CENSORING:
            raise LikelihoodError(
                f"the {family.NAME} family does not consume a censoring declaration, but one was "
                f"given. Families that do: "
                f"{sorted(n for n, c in _FAMILIES.items() if c.SUPPORTS_CENSORING)}."
            )
        self._family = family
        self._noise = noise
        self._censoring = censoring
        # Latency that follows from the family and noise model alone is a fact
        # about the pair, so it is refused here. Latency that only censoring
        # induces depends on whether a limit survives the mask, so it is
        # refused in check_alignment, which has the data.
        if family.marginalisation_with(noise) is Marginalisation.LATENT:
            self._refuse_unconsumed_latent(family, noise, censored=False)
        for parameter in noise.parameters:
            self.register_parameter(parameter)
        for parameter in family.parameters:
            if parameter.name in self.parameters:
                raise LikelihoodError(
                    f"the {family.NAME} family and the {type(noise).__name__} noise model both "
                    f"declare a parameter named {parameter.name!r}. A Likelihood holds one flat "
                    f"namespace (parameters.md §12.4: merge is not associative, so this contract "
                    f"does not nest one); rename one of them."
                )
            self.register_parameter(parameter)

    @staticmethod
    def _refuse_unconsumed_latent(
        family: LikelihoodFamily, noise: NoiseModel, *, censored: bool
    ) -> None:
        """Refuse a latent combination whose family does not implement one.

        The failure this prevents is the worst kind available to a likelihood
        contract: a family that declares :attr:`Marginalisation.LATENT` but
        whose ``log_prob`` never reads ``noise.latent`` returns the
        *uncorrelated* likelihood, so every GP hyperparameter and every latent
        value an engine samples leaves the log-probability untouched. Their
        posteriors come back as their priors and the physical parameters are
        biased exactly as they would have been under the rigid likelihood — a
        fit that runs, converges, and is wrong, which is what this whole
        contract exists to make impossible.
        """
        if family.CONSUMES_LATENT_GP:
            return
        because = (
            "with a censoring declaration on correlated noise"
            if censored
            else f"with a {type(noise).__name__} noise model"
        )
        raise LikelihoodError(
            f"the {family.NAME} family {because} is a latent-variable model (see "
            f"Marginalisation.LATENT), but {type(family).__name__} does not implement the "
            f"latent-conditional log_prob — its CONSUMES_LATENT_GP is False. Evaluating it would "
            f"return the *uncorrelated* likelihood, leaving the GP hyperparameters and every "
            f"latent value with no effect on the log-probability at all: a fit that runs, "
            f"converges, and is wrong. Use IndependentNoise with this family, or the gaussian "
            f"family, whose GP marginalises in closed form. "
            f"family.marginalisation_with(noise, censoring) reports the declaration without "
            f"composing anything."
        )

    # -- declarations --------------------------------------------------------

    @property
    def family(self) -> LikelihoodFamily:
        """The sampling distribution."""
        return self._family

    @property
    def noise(self) -> NoiseModel:
        """The noise model."""
        return self._noise

    @property
    def censoring(self) -> Censoring | None:
        """The censoring declaration, if any."""
        return self._censoring

    @property
    def marginalisation(self) -> Marginalisation:
        """Whether this combination integrates its noise process in closed form.

        Answered from the *declaration* alone, before any data are seen, and so
        deliberately conservative about censoring: a limit that a container's
        mask happens to exclude still counts here, because this property does
        not know which container it will meet. Use
        :meth:`marginalisation_for` once the observed container is in hand —
        which is the answer W1.7 should act on.
        """
        return self._family.marginalisation_with(self._noise, self._censoring)

    def marginalisation_for(self, observed: FunctionSamples) -> Marginalisation:
        """The marginalisation this combination actually needs for these data.

        Identical to :attr:`marginalisation` except that censoring is read
        through the mask, so a limit on a masked sample does not force the
        latent path. §9's rule is that masking beats censoring — a masked
        sample contributes nothing whatever its limit kind — and that has to
        hold for the declaration as well as for the arithmetic, or a problem
        that is analytic in fact gets refused a gradient-free engine.
        """
        if self._censoring is None:
            return self._family.marginalisation_with(self._noise, None)
        self._censoring.check_against(observed)
        retained = Censoring(np.asarray(self._censoring.kinds)[np.asarray(observed.valid).ravel()])
        return self._family.marginalisation_with(self._noise, retained)

    def latent_declaration(self, size: int, name: str = "latent") -> LatentDeclaration:
        """The latent parameters this combination adds, for ``size`` samples.

        Raises when the combination is analytic: an analytic likelihood has no
        latent values, and quietly returning an empty declaration would let a
        caller build a sampler dimension that does nothing.
        """
        if self.marginalisation is Marginalisation.ANALYTIC:
            raise LikelihoodError(
                f"the {self._family.NAME} family with a {type(self._noise).__name__} noise model "
                f"marginalises analytically, so it declares no latent values. "
                f"Likelihood.marginalisation tells you which path you are on."
            )
        if not isinstance(self._noise, GaussianProcessNoise):
            raise LikelihoodError(
                f"a latent declaration needs a GaussianProcessNoise model to supply the "
                f"whitening transform, but this Likelihood has a {type(self._noise).__name__}."
            )
        return LatentDeclaration(
            parameter=latent_parameter(name, size),
            size=int(size),
            solver=self._noise.solver,
        )

    def check_engine(
        self,
        *,
        differentiable: bool,
        engine: str = "this engine",
        observed: FunctionSamples | None = None,
    ) -> None:
        """Refuse an engine that cannot deliver this likelihood's marginalisation.

        ``DEVELOPMENT_PLAN.md`` §4.4: "gradient-free samplers cannot
        realistically handle hundreds of latent values, so non-Gaussian +
        flexible-GP robustness is effectively a modern-backend capability". A
        silent downgrade here would be a fit that runs, converges to something,
        and is wrong.

        Pass ``observed`` when the data are to hand — W1.7's ``Dataset`` always
        has them — so that a censored sample the mask excludes is not counted
        against a gradient-free engine.
        """
        required = self.marginalisation if observed is None else self.marginalisation_for(observed)
        if required is Marginalisation.ANALYTIC or differentiable:
            return
        raise LikelihoodError(
            f"{engine} cannot run this likelihood: the {self._family.NAME} family with a "
            f"{type(self._noise).__name__} noise model introduces N latent values that inference "
            f"must marginalise numerically, one per retained sample. That needs a "
            f"gradient-based engine (HMC/NUTS or VI on the torch/jax backends). Options: use "
            f"IndependentNoise with this engine, or a Gaussian family, whose GP marginalises in "
            f"closed form."
        )

    # -- composition-time checking -------------------------------------------

    def check_alignment(self, predicted: FunctionSamples, observed: FunctionSamples) -> None:
        """Verify predicted and observed describe the same samples. O(N), once.

        This is deliberately not called from :meth:`log_prob`: coordinate
        comparison is a composition-time obligation, matching
        ``results_schema.md`` §10's compile-once/evaluate-many split. W1.7's
        ``Dataset`` calls it when the problem is assembled.
        """
        self._check_kinds(predicted, observed)
        if predicted.shape != observed.shape:
            raise LikelihoodError(
                f"the predicted {type(predicted).__name__} has shape {predicted.shape} but the "
                f"observed one has {observed.shape}. A likelihood compares samples index by "
                f"index; a resampling step belongs in the instrument chain (W1.5)."
            )
        if predicted.unit != observed.unit:
            raise LikelihoodError(
                f"the predicted {type(predicted).__name__} is in {predicted.unit} but the "
                f"observed one is in {observed.unit}. Convert once at composition time with "
                f".to_unit(...); units never reach the hot loop (DEVELOPMENT_PLAN.md §7)."
            )
        # Ruled 2026-09-02 (W1.11 gap I-1): kind equality does not imply
        # comparability — VisibilitySet is legal with real or complex values,
        # so a complex prediction could otherwise be fitted against real
        # amplitudes with no warning anywhere.
        if predicted.values.dtype.kind != observed.values.dtype.kind:
            raise LikelihoodError(
                f"the predicted {type(predicted).__name__} holds "
                f"{'complex' if predicted.values.dtype.kind == 'c' else 'real'} values but the "
                f"observed one holds "
                f"{'complex' if observed.values.dtype.kind == 'c' else 'real'} values. A "
                f"likelihood compares like with like: if you mean to fit amplitudes, take the "
                f"modulus in the instrument chain (W1.5) so both sides are real."
            )
        for left, right in zip(predicted.axes, observed.axes, strict=True):
            if left != right:
                raise LikelihoodError(
                    f"the predicted and observed {left.name!r} axes differ. The likelihood "
                    f"assumes they index the same samples; make the instrument chain (W1.5) "
                    f"produce the observed container's own coordinates."
                )
        if self._censoring is not None:
            self._censoring.check_against(observed)
            if not self._family.SUPPORTS_CENSORING:
                raise LikelihoodError(
                    f"the {self._family.NAME} family cannot consume a censoring declaration."
                )
            # Censoring on correlated noise turns the likelihood into a
            # multivariate-normal orthant probability, which is latent. Whether
            # it actually does depends on whether a limit survives the mask, so
            # this is the earliest point the question can honestly be asked.
            if self.marginalisation_for(observed) is Marginalisation.LATENT:
                self._refuse_unconsumed_latent(self._family, self._noise, censored=True)
        self._noise.check_compatible(self._family, observed)
        self._family.check_observed(observed)
        if observed.values.dtype.kind == "c" and not self._family.ALLOWS_COMPLEX:
            raise LikelihoodError(
                f"the {self._family.NAME} family holds real values, but the observed "
                f"{type(observed).__name__} is complex. Use the complex_gaussian family, or an "
                f"explicit real projection (amplitude, phase, real, imag)."
            )

    def _check_kinds(self, predicted: FunctionSamples, observed: FunctionSamples) -> None:
        related = isinstance(predicted, type(observed)) or isinstance(observed, type(predicted))
        # A subclass may legitimately be compared with its base, but only if it
        # kept the same axis signature; otherwise the axis loop below would be
        # comparing different things under the same name.
        if not related or len(predicted.axes) != len(observed.axes):
            raise LikelihoodError(
                f"the predicted container is a {type(predicted).__name__} but the observed one "
                f"is a {type(observed).__name__}. A likelihood compares like with like; channel "
                f"binding (ModelResult.require) is where a kind mismatch should have been caught."
            )

    # -- evaluation ----------------------------------------------------------

    def log_prob(
        self,
        predicted: FunctionSamples,
        observed: FunctionSamples,
        values: Mapping[str, Any] | Sequence[float] | np.ndarray | None = None,
        *,
        latent: np.ndarray | None = None,
    ) -> float:
        """log p(observed | predicted, parameters).

        **The mask convention, stated once.** Masks arrive through
        ``weights()`` — ``results_schema.md`` §16 asked this contract to pick
        one of ``weights()`` and ``masked_uncertainty()`` and say which, and it
        is ``weights()``. The effective weight is the product of the predicted
        and observed containers' weights, so a sample counts only if both call
        it valid, and because those weights are exactly 0 or 1 the product is
        an inclusion indicator. Retained samples are then **excised**: they are
        the only rows and columns that enter the covariance at all.

        Excision, not the infinite-variance limit, is the rule — for *every*
        noise model, not only the correlated one. ``masked_uncertainty()`` and
        a zero weight are equivalent statements about a **chi-square term**,
        which is exactly what ``results_schema.md`` §7 claims for them; they
        are not equivalent statements about a normalised log-density, and a
        likelihood is a normalised log-density. Each Gaussian term carries a
        ``-log sigma`` alongside its chi-square, so
        ``-1/2 log(2 pi) - log sigma - 1/2 (r/sigma)**2`` tends to ``-inf``
        rather than to zero as ``sigma`` grows. The GP case is the same
        statement with a log-determinant in place of the ``-log sigma``: the
        limit approaches the excised value minus ``1/2 log(2 pi sigma**2)``,
        and so diverges. Since §7 defines a masked sample as carrying *exactly*
        zero information, only excision delivers that, and a 0/1 weight
        multiplying a whole term is excision written arithmetically.

        A fully masked pair returns ``0.0``: no data, no information, no
        contribution — which is the correct limit of the same rule, not a
        special case.

        Masked coordinates are still available as evaluation locations for
        :meth:`conditional`, so a GP-localisation diagnostic can show the
        conditioned mean where the data were excluded.
        """
        if latent is not None and self.marginalisation is Marginalisation.ANALYTIC:
            raise LikelihoodError(
                "latent values were supplied, but this likelihood marginalises analytically and "
                "has none. Check Likelihood.marginalisation before allocating them."
            )
        self._check_shapes(predicted, observed)
        resolved = self.context(values)

        weights = np.asarray(observed.weights()).ravel() * np.asarray(predicted.weights()).ravel()
        retain = weights > 0.0
        if not np.any(retain):
            return 0.0

        observed_values = self._retained(observed, retain, "observed")
        predicted_values = self._retained(predicted, retain, "predicted")
        limits = self._retained_limits(retain)
        coordinates = self._coordinates(observed, retain) if self._noise.CORRELATED else None

        noise = self._noise.noise_params(
            observed,
            retain,
            resolved,
            coordinates=coordinates,
            latent=latent,
            limits=limits,
        )
        return float(self._family.log_prob(predicted_values, observed_values, noise))

    def conditional(
        self,
        predicted: FunctionSamples,
        observed: FunctionSamples,
        values: Mapping[str, Any] | Sequence[float] | np.ndarray | None = None,
        *,
        at: np.ndarray | None = None,
    ) -> GPConditional:
        """The conditioned GP mean and variance — the localisation diagnostic.

        ``DEVELOPMENT_PLAN.md`` §4.8's third diagnostic family: "the
        conditioned GP mean ... already localises where the model is
        deficient". Defaults to the *full* coordinate set including masked
        samples, because "what would the GP have said here?" is exactly the
        question a user asks about a region they excluded.
        """
        if not isinstance(self._noise, GaussianProcessNoise):
            raise LikelihoodError(
                f"a conditioned GP mean needs a GaussianProcessNoise model; this Likelihood has "
                f"a {type(self._noise).__name__}, which induces no correlations to condition."
            )
        self._check_shapes(predicted, observed)
        resolved = self.context(values)
        weights = np.asarray(observed.weights()).ravel() * np.asarray(predicted.weights()).ravel()
        retain = weights > 0.0
        if not np.any(retain):
            raise LikelihoodError(
                "every sample is masked, so there is nothing to condition the GP on."
            )
        coordinates = self._coordinates(observed, retain)
        residual = self._retained(observed, retain, "observed") - self._retained(
            predicted, retain, "predicted"
        )
        sigma = self._noise.sigma(observed, retain, resolved)
        if sigma is None:
            raise LikelihoodError(
                f"conditioning the GP on residuals needs the observed "
                f"{type(observed).__name__}'s per-sample uncertainties, but it has none. A "
                f"latent-GP combination localises through its latent posterior instead, which "
                f"is inference's output rather than this method's."
            )
        # Default to the *full* axis, masked samples included; otherwise hand
        # `at` to the solver untouched, so one place (``_as_points``) decides
        # what a bare 1-D array of coordinates means.
        target = (
            self._coordinates(observed, np.ones(observed.n_samples, dtype=bool))
            if at is None
            else at
        )
        return self._noise.solver.condition(
            self._noise.kernel, coordinates, residual, sigma**2, resolved, at=target
        )

    # -- internals -----------------------------------------------------------

    def _check_shapes(self, predicted: FunctionSamples, observed: FunctionSamples) -> None:
        if predicted.shape != observed.shape:
            raise LikelihoodError(
                f"the predicted {type(predicted).__name__} has shape {predicted.shape} but the "
                f"observed one has {observed.shape}. Call check_alignment() at composition time "
                f"to catch this before a sampler starts."
            )

    def _retained(self, container: FunctionSamples, retain: np.ndarray, role: str) -> np.ndarray:
        flat = np.asarray(container.values).ravel()[retain]
        if flat.dtype.kind == "c":
            if not self._family.ALLOWS_COMPLEX:
                raise LikelihoodError(
                    f"the {self._family.NAME} family holds real values, but the {role} container "
                    f"is complex. Use the complex_gaussian family."
                )
            complex_values = np.ascontiguousarray(flat, dtype=np.complex128)
            return _check_finite(complex_values, f"{role} values")
        return _as_float64(flat, f"{role} values")

    def _retained_limits(self, retain: np.ndarray) -> np.ndarray | None:
        if self._censoring is None:
            return None
        if self._censoring.n_samples != retain.size:
            raise LikelihoodError(
                f"the censoring declaration covers {self._censoring.n_samples} samples but the "
                f"containers hold {retain.size}."
            )
        return np.asarray(self._censoring.kinds)[retain]

    def _coordinates(self, observed: FunctionSamples, retain: np.ndarray) -> np.ndarray:
        if observed.LAYOUT is not Layout.POINTS:
            raise LikelihoodError(
                f"a correlated noise model needs point-set coordinates, but a "
                f"{type(observed).__name__} has a {observed.LAYOUT.value} layout. Gridded 2D+ "
                f"data are the SVGP / SKI / Vecchia strategy slots (Phase 5)."
            )
        stacked = np.column_stack([np.asarray(axis.values, dtype=DTYPE) for axis in observed.axes])
        return np.ascontiguousarray(stacked[retain])

    def __repr__(self) -> str:
        censored = "" if self._censoring is None else f", {self._censoring!r}"
        return (
            f"Likelihood({type(self._family).__name__}, {type(self._noise).__name__}"
            f"{censored}, marginalisation={self.marginalisation.value})"
        )
