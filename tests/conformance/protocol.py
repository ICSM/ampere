"""The backend fixture protocol the conformance battery is parametrised over.

``DEVELOPMENT_PLAN.md`` §4.6 makes this suite "the contract that keeps two
lockstep backends aligned across agent tracks". For that to be more than an
aspiration, the battery has to be written once and run once per backend — so
every test body in ``tests/conformance/`` takes the ``backend`` fixture and
never names a concrete implementation.

This module is the whole of what a backend must supply. It is deliberately
small, and it is derived from what the rows actually need rather than from a
guess at what a backend might one day want:

* the **declarations** a row builds a problem from (:class:`ModelSpec`,
  :class:`TransformationSpec`, :class:`CovarianceSpec`) are backend-neutral
  data; the backend turns each into a concrete object;
* the **contracts** (``ParameterSet``, ``Likelihood``, ``Dataset``,
  ``FittingProblem``) are reused unchanged, which is the claim
  ``inference.md`` §18 and ``results.md`` §14 both make for them — "a backend
  supplies models and transformations […] and nothing else", "a backend
  supplies draws and ``Evaluation``s and nothing else";
* the **capabilities and tolerances** (:class:`BackendCapabilities`,
  :class:`Tolerances`) are the declaration that lets a row skip cleanly rather
  than fail spuriously, and the per-comparison tolerance table
  ``architecture.md`` §2 point 1 asks W1.10 to own.

``tests/conformance/README.md`` is the prose version, written for a Phase-2
backend author. Read it before implementing :class:`ConformanceBackend`.
"""

from __future__ import annotations

import dataclasses
import enum
from collections.abc import Mapping, Sequence
from typing import Any, Protocol, runtime_checkable

import numpy as np

from ampere.core import (
    GPSolver,
    Kernel,
    Model,
    NoiseModel,
    ParameterSet,
    Transformation,
)

__all__ = [
    "BackendCapabilities",
    "ConformanceBackend",
    "CountingModel",
    "CovarianceSpec",
    "KernelFamily",
    "ModelKind",
    "ModelSpec",
    "ParameterSpace",
    "SolverKind",
    "Tolerances",
    "TransformationKind",
    "TransformationSpec",
    "complex_axes",
]


# ---------------------------------------------------------------------------
# Tolerances
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class Tolerances:
    """Per-comparison-class absolute tolerances, as ``architecture.md`` §2 asks.

    That section commits the project to "documented and specific per
    comparison, never a single blanket tolerance", and hands the table to
    W1.10. One field per *kind* of comparison, so a row names the class of
    agreement it is asserting rather than inventing a number:

    ``exact``
        Identities that hold in exact arithmetic up to representation error
        only: ``pack``∘``unpack``, ``constrain``∘``unconstrain``, a spec round
        trip, a mask surviving a transformation. Anything looser here is a bug,
        not a tolerance.
    ``analytic``
        Agreement with a closed form evaluated on the *same* arithmetic path —
        the i.i.d. Gaussian written out, ``scipy.stats.norm.logcdf`` for a
        Tobit limit, a prior quantile from ``scipy``'s own ``ppf``.
    ``linear_algebra``
        Agreement with a closed form that takes a *different* factorisation of
        the same matrix: ``DenseGP`` against
        ``scipy.stats.multivariate_normal.logpdf``, or ``L Lᵀ`` against ``K``.
    ``cross_solver``
        Two ``GPSolver`` strategies computing the same marginal likelihood by
        genuinely different recursions — dense Cholesky against the
        quasiseparable state-space solve. Not the same as ``linear_algebra``:
        the operation count differs by orders of magnitude and so does the
        accumulated rounding.
    ``cross_backend``
        Two backends' independent arithmetic, plus (in Phase 2) a different
        autodiff accumulation order. A float32 backend widens this one, which
        is why it lives on :class:`BackendCapabilities` and not in a constant.
    ``monte_carlo_sigmas``
        Not a tolerance but a multiple: how many standard errors an empirical
        moment may sit from its analytic value before the row fails. The
        standard error itself is computed from the estimator, so the assertion
        stays honest as the draw count changes.
    """

    exact: float = 1e-12
    analytic: float = 1e-9
    linear_algebra: float = 1e-8
    cross_solver: float = 1e-6
    cross_backend: float = 1e-9
    monte_carlo_sigmas: float = 5.0


DEFAULT_TOLERANCES = Tolerances()


# ---------------------------------------------------------------------------
# Declarations a backend realises
# ---------------------------------------------------------------------------


class SolverKind(enum.StrEnum):
    """The GP solve strategies a row may ask a backend for."""

    DENSE = "dense"
    QUASISEP = "quasisep"


class KernelFamily(enum.StrEnum):
    """Kernel families, split by whether they have a quasiseparable form."""

    MATERN32 = "matern32"
    SQUARED_EXPONENTIAL = "squared_exponential"


class ModelKind(enum.StrEnum):
    """The analytic model forms the battery composes problems from.

    All are closed forms, so every likelihood row has an oracle that never
    calls ampere. Between them they cover the three bijections
    ``default_bijection_for`` can pick.

    ``LINEAR``
        ``f(x) = offset + slope * x``, with ``offset ~ Normal(0, 1)`` and
        ``slope ~ Normal(1, 0.5)`` — both unbounded, so both lower to
        ``Identity``.
    ``POWER_LAW``
        ``f(x) = norm * (x / x_ref) ** index``, with
        ``norm ~ LogUniform(0.1, 10)`` (bounded below and above, so ``Logit``)
        and ``index ~ Normal(-1, 0.5)`` (``Identity``).
    ``COMPLEX``
        ``V(x) = norm * exp(2 pi i * index * x)`` — a **complex**-valued
        channel, emitted as a :class:`~ampere.core.VisibilitySet` rather than a
        ``Spectrum``. Added at W2.4 slice 3, for the ``complex_gaussian``
        family: it is the smallest form whose values are complex and whose
        parameters still move both the modulus and the phase, so a backend that
        dropped the imaginary part somewhere would disagree with the oracle
        rather than merely lose precision. The parameter names and priors are
        deliberately ``POWER_LAW``'s, so nothing else in the battery has to
        learn a third vocabulary.

        Only backends declaring
        :attr:`BackendCapabilities.complex_models` are asked for one; every
        row that uses it skips elsewhere, so it is additive in exactly the way
        :attr:`BackendCapabilities.solvers` is.
    """

    LINEAR = "linear"
    POWER_LAW = "power_law"
    COMPLEX = "complex"


class TransformationKind(enum.StrEnum):
    """The instrument steps the battery composes chains from.

    ``SCALE``
        A multiplicative calibration, ``ACCEPTS = (Spectrum,)`` and
        ``PRODUCES = None``, carrying one free parameter ``scale`` with a
        ``LogNormal(0.2)`` prior (support ``(0, ∞)``, so a ``Log`` bijection).
        Publishes no requirements.
    ``REBIN``
        Resamples a ``Spectrum`` onto the coarser grid in
        :attr:`TransformationSpec.target` with an explicit influence matrix.
        It publishes an :class:`~ampere.core.AxisRequirement` covering the
        target grid, and it must propagate the input mask with
        :func:`~ampere.core.propagate_mask` — this is the step the ANY-rule and
        dropped-mask rows exercise.
    ``PHOTOMETRY``
        ``Spectrum → PhotometricPoints``: ``ACCEPTS = (Spectrum,)`` and
        ``PRODUCES = PhotometricPoints``, integrating onto the pivots in
        :attr:`TransformationSpec.target` under the names in
        :attr:`TransformationSpec.filters`. It exists so a row can build a
        chain whose kinds do not compose.
    """

    SCALE = "scale"
    REBIN = "rebin"
    PHOTOMETRY = "photometry"


@dataclasses.dataclass(frozen=True)
class ModelSpec:
    """A backend-neutral declaration of a forward model.

    ``coordinates`` are micron and the emitted values are Jy, on every channel
    and in every backend; fixing the units here keeps the containers
    comparable across backends without making unit negotiation a backend
    author's problem.
    """

    kind: ModelKind = ModelKind.LINEAR
    channels: tuple[str, ...] = ("default",)
    coordinates: tuple[float, ...] = (1.0, 1.6, 2.5, 4.0, 6.3, 10.0)
    reference_coordinate: float = 1.0
    plated: bool = False
    """Whether the model additionally declares a plate of per-channel offsets.

    ``False`` — the default — is the scalar-only declaration. ``True`` asks
    for one extra ``Plate`` named ``"objects"``, of size ``len(channels)``,
    with hyperparameters ``mu ~ Normal(0, 1)`` and ``sigma ~ HalfNormal(1)``
    and one array-valued member ``offsets`` drawn from
    ``Normal(mu, sigma)``; channel *i* adds ``offsets[i]`` to its flux. It
    exists because ``results.md`` §14 asks for a row proving an array-valued
    parameter emits **one** variable carrying a **named** dimension, and a
    plate is the declaration that names one.
    """

    compiled: bool = False
    """Whether the model honours :meth:`~ampere.core.Model.compile_for`.

    ``False`` — the default — leaves ``compile_for`` inherited, so it returns
    ``self`` and the model evaluates on its own declared grid. ``True`` asks
    for the template-caching behaviour ``transformations.md`` §14 names: build
    one container per channel from ``requirements[channel].coordinates()``,
    keep it, and refill it with ``with_values`` on every evaluation, so
    successive results share their axes by identity.
    """


@dataclasses.dataclass(frozen=True)
class TransformationSpec:
    """A backend-neutral declaration of one instrument step."""

    kind: TransformationKind = TransformationKind.SCALE
    label: str | None = None
    target: tuple[float, ...] = ()
    filters: tuple[str, ...] = ()


@dataclasses.dataclass(frozen=True)
class CovarianceSpec:
    """A backend-neutral declaration of a GP kernel.

    ``amplitude`` is the marginal standard deviation (``k(0) == amplitude**2``)
    and ``length_scale`` is in the coordinate's own units — micron, matching
    :attr:`ModelSpec.coordinates`. Both are given as plain numbers, held
    fixed: the battery's GP rows compare numbers against ``scipy``, and a
    fitted hyperparameter would only add a sampling dimension they do not use.
    """

    family: KernelFamily = KernelFamily.MATERN32
    amplitude: float = 0.4
    length_scale: float = 2.0


# ---------------------------------------------------------------------------
# Adapters
# ---------------------------------------------------------------------------


@runtime_checkable
class ParameterSpace(Protocol):
    """One backend's realisation of a parameter declaration.

    ``ampere.core.ParameterSet`` already satisfies every method below, so the
    reference adapter is a thin wrapper. A Phase-2 backend that lowers a
    declaration into native sample sites (numpyro, pyro) returns its own
    object here, and the ``test_parameters.py`` rows then hold it to the same
    behaviour — which is the point: ``parameters.md`` §13 hands those rows over
    precisely because the lowered path must reproduce them.

    Every method takes and returns plain numpy on the boundary. A backend that
    works in its own array type converts at the edge; the battery is asserting
    agreement of *values*, not of container types.
    """

    @property
    def declaration(self) -> ParameterSet:
        """The backend-neutral declaration this space realises."""

    @property
    def free_size(self) -> int:
        """Number of flat free dimensions — the engine's dimension."""

    def free_labels(self) -> tuple[str, ...]:
        """One label per flat free dimension, array elements included."""

    def pack(self, values: Mapping[str, Any]) -> np.ndarray:
        """Flatten a name → value mapping into the free vector."""

    def unpack(self, theta: Any) -> dict[str, Any]:
        """Expand a free vector into every parameter, fixed ones included."""

    def prior_transform(self, unit_cube: Any) -> np.ndarray:
        """Map ``[0, 1]^n`` to the constrained free vector (nested sampling)."""

    def lnprior(self, values: Any) -> float:
        """Log prior density in the constrained parameterisation."""

    def constrain(self, unconstrained: Any) -> np.ndarray:
        """Map the unconstrained vector into the priors' support."""

    def unconstrain(self, values: Any) -> np.ndarray:
        """Map values in the priors' support to the unconstrained vector."""

    def lnprior_unconstrained(self, unconstrained: Any) -> float:
        """Log prior density in unconstrained space, Jacobian term included."""


@runtime_checkable
class CountingModel(Protocol):
    """A :class:`~ampere.core.Model` that records how often it was evaluated.

    ``inference.md`` §18 asks for a row proving an out-of-support θ returns
    ``-inf`` *without evaluating the model* — the single most valuable thing
    the contract does for an expensive simulator. There is no way to assert
    that from outside without the model saying so, so every model a backend
    returns from :meth:`ConformanceBackend.model` must maintain this counter:
    increment ``evaluations`` once per call to ``evaluate``, and reset it to
    zero on :meth:`reset_evaluations`.
    """

    evaluations: int

    def reset_evaluations(self) -> None:
        """Set :attr:`evaluations` back to zero."""


@dataclasses.dataclass(frozen=True)
class BackendCapabilities:
    """What a backend can do, and how closely it is expected to agree.

    ``differentiable``, ``batchable`` and ``device`` mirror
    :class:`ampere.core.Capabilities` — they are the three flags
    ``inference.md`` §18 says a backend declares and nothing else. The rest is
    what the battery needs in order to skip a row honestly instead of failing
    it:

    ``solvers``
        Which :class:`SolverKind` values :meth:`ConformanceBackend.gp_solver`
        can actually return an implemented solver for. ``ampere.core``'s
        ``QuasisepGP`` has been real since W2.3 (celerite2's numpy solver over
        an exact rank-2 Matérn-3/2 representation), so the in-repo fixtures
        declare ``QUASISEP`` and the ``DenseGP``↔``QuasisepGP`` row runs. A
        backend that has not yet supplied one omits it, and that row skips
        with a reason naming exactly what is missing.
    ``float64``
        Whether the backend's likelihood linear algebra runs in double
        precision. ``architecture.md`` §5 makes float64 the policy for GP
        solves; a backend that opts out for GPU throughput declares it here
        and widens ``tolerances.cross_backend`` accordingly.
    ``complex_models``
        Whether :meth:`ConformanceBackend.model` can realise
        :attr:`ModelKind.COMPLEX` — a channel of complex values, which the
        ``complex_gaussian`` rows need. ``False`` by default, so a backend that
        has not written one is unaffected and those rows skip with a reason;
        the alternative, a required fixture method, would have made a
        Phase-2 track's own slice a change to every other track's fixture.
    ``tolerances``
        The per-comparison table (:class:`Tolerances`). It lives on the
        backend, not in a module constant, because a float32 or
        different-accumulation-order backend legitimately needs a looser
        ``cross_backend`` without loosening anyone else's.
    """

    differentiable: bool = False
    batchable: bool = False
    device: str = "cpu"
    float64: bool = True
    complex_models: bool = False
    solvers: frozenset[SolverKind] = frozenset({SolverKind.DENSE})
    tolerances: Tolerances = DEFAULT_TOLERANCES


@runtime_checkable
class ConformanceBackend(Protocol):
    """Everything the conformance battery asks of a backend.

    Nine members. Implement them and every row in ``tests/conformance/``
    runs against your backend; register the instance in
    ``tests/conformance/backends/__init__.py`` and nothing else changes —
    which is W1.10's acceptance criterion ("adding a backend requires only a
    fixture") and the reason no test body may name a concrete backend.
    """

    @property
    def name(self) -> str:
        """The pytest fixture id. Short, lower-case, stable across runs.

        **And the backend's name everywhere else** (W2.12): every model and
        instrument step :meth:`model` and :meth:`transformation` return must
        declare ``BACKEND`` equal to this string, so that a composed problem
        reports it, a run's ``ampere_backend`` carries it, and ``lowering.md``
        §12.8's registry is keyed on it. ``TestBackendIdentity`` asserts it.
        """

    @property
    def capabilities(self) -> BackendCapabilities:
        """What this backend can do, and its tolerance table."""

    def model(self, spec: ModelSpec) -> Model:
        """Realise a :class:`ModelSpec` as a concrete model.

        The returned object must be a :class:`~ampere.core.Model` *and* a
        :class:`CountingModel`. It declares exactly the parameters
        :class:`ModelKind` names, under exactly those names — the ``lnprior``,
        ``prior_transform`` and spec-hash rows compare against ``scipy`` and
        against other backends by name, so a renamed parameter is a failure,
        not a detail.
        """

    def transformation(self, spec: TransformationSpec) -> Transformation:
        """Realise a :class:`TransformationSpec` as a concrete instrument step.

        ``ACCEPTS``/``PRODUCES``, the published requirements and the mask
        propagation are all part of the declaration — see
        :class:`TransformationKind` for what each kind owes.
        """

    def kernel(self, spec: CovarianceSpec) -> Kernel:
        """Realise a :class:`CovarianceSpec` as a concrete kernel.

        The kernel's ``QUASISEPARABLE`` class flag must be truthful: the
        ``DenseGP``↔``QuasisepGP`` row selects on it, and
        ``GPSolver.check_compatible`` refuses a quasiseparable solver a kernel
        that has no such representation.
        """

    def gp_solver(self, kind: SolverKind) -> GPSolver:
        """Return this backend's solver for *kind*.

        Only called for a kind present in ``capabilities.solvers``; a backend
        may raise for anything else.
        """

    def independent_noise(self) -> NoiseModel:
        """This backend's uncorrelated noise model, with no scale or jitter.

        **Added at W2.13.** A :class:`~ampere.core.NoiseModel` carries the four
        capability flags from that item on (``inference.md`` §10a, fold-in 7)
        and :attr:`ampere.core.Likelihood.capability_parts` puts it on the
        composed problem, so a fixture that kept composing
        ``ampere.core.IndependentNoise`` under its own name would be declaring
        two backends and refused at composition — correctly, because a
        nominally native problem whose noise arithmetic ran in numpy is
        exactly what that widening exists to catch.
        """

    def gp_noise(self, kernel: Kernel, solver: GPSolver) -> NoiseModel:
        """This backend's GP noise composition over *kernel* and *solver*.

        **Added at W2.13**, for the reason :meth:`independent_noise` gives.
        *kernel* and *solver* are this fixture's own, from :meth:`kernel` and
        :meth:`gp_solver`.
        """

    def parameter_space(self, declaration: ParameterSet) -> ParameterSpace:
        """Realise a parameter declaration as this backend's engine-facing view."""

    def to_numpy(self, values: Any) -> np.ndarray:
        """Bring a backend array back to numpy, for comparison against oracles.

        The battery's oracles are ``scipy`` and closed forms, which speak
        numpy. A backend working in torch or jax arrays converts here — and
        only here, so no row has to know.
        """


def complex_axes(coordinates: Any) -> tuple[np.ndarray, np.ndarray]:
    """The ``(u, v)`` point set a :attr:`ModelKind.COMPLEX` channel lives on.

    A :class:`~ampere.core.VisibilitySet` has two coordinate axes and a
    ``ModelSpec`` declares one sequence of coordinates, so the second has to
    come from somewhere. It comes from here rather than from either side,
    because ``Likelihood.check_alignment`` compares the predicted and observed
    axes for equality: the fixture's model and the battery's observed container
    must derive them by the *same* rule, and a rule written twice is a rule
    that eventually differs.

    ``v = u / 2`` is arbitrary and deliberately so — the battery's complex rows
    are about the family's arithmetic, not about (u,v) coverage — but it is
    distinct from ``u``, which keeps the two axes from being accidentally
    interchangeable in a comparison.
    """
    u_axis = np.asarray(coordinates, dtype=float)
    return u_axis, 0.5 * u_axis


def solver_kinds(capabilities: BackendCapabilities) -> tuple[SolverKind, ...]:
    """The solver kinds *capabilities* declares, in :class:`SolverKind` order."""
    return tuple(kind for kind in SolverKind if kind in capabilities.solvers)


def as_sequence(values: Any) -> Sequence[float]:
    """Coerce *values* to a plain float sequence, for spec fields."""
    return tuple(float(value) for value in np.asarray(values, dtype=float).ravel())
