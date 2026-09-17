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

import astropy.units as astropy_units
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
    "COMPLEX_WAVELENGTH",
    "AstrometryPieces",
    "BackendCapabilities",
    "ConformanceBackend",
    "CountingModel",
    "CovarianceSpec",
    "InterferometryPieces",
    "KernelFamily",
    "ModelKind",
    "ModelSpec",
    "ParameterSpace",
    "SolverKind",
    "Tolerances",
    "TransformationKind",
    "TransformationSpec",
    "approximation_envelope",
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

    ``approximation_order``, ``approximation_floor``, ``approximation_final``
        **W5.4's class, for an ``EXACT = False`` solver.** These three do not
        describe an agreement at all; they describe a *convergence*, which is
        the only honest thing to assert about an approximation. See
        :func:`approximation_envelope` and ``README.md`` §3.
    """

    exact: float = 1e-12
    analytic: float = 1e-9
    linear_algebra: float = 1e-8
    cross_solver: float = 1e-6
    cross_backend: float = 1e-9
    monte_carlo_sigmas: float = 5.0
    approximation_order: float = 1.0
    approximation_floor: float = 1e-2
    approximation_final: float = 5e-2


def approximation_envelope(
    coarsest_error: float,
    coarsest_size: int,
    size: int,
    tolerances: Tolerances,
    order: float | None = None,
) -> float:
    """The tolerance an approximate solver's error must sit inside at *size*.

    **The ``EXACT = False`` comparison class** (W5.4; ``horizon_notes.md`` §2
    question (a)). An approximate solver cannot be held to a number: how close
    it comes to :class:`~ampere.core.DenseGP` depends on the kernel, the data
    and the approximation's own parameters, and a fixed tolerance would either
    be so loose it asserted nothing or so tight it encoded one fixture's
    arithmetic. What *is* a property of the method — and what a wrong
    implementation breaks — is that the error **falls as the approximation is
    refined**, at a rate the method's own analysis predicts.

    So the envelope is anchored on the error measured at the **coarsest**
    setting of the same row, and tightens from there as
    ``(coarsest_size / size) ** order``: refining by a factor of two must buy
    at least ``2 ** order``. Below that sits ``approximation_floor``, because
    a Hilbert-space basis converges to the kernel *on a finite box* and the
    box's own truncation error does not go away with more basis members — an
    envelope with no floor would assert something false about the method
    rather than something demanding about the implementation. The floor also
    protects against the other direction: a coarsest setting that happens to
    agree well by accident would otherwise set an envelope nothing could meet.

    **The rate is a property of the kernel, not of the solver**, which is why
    ``approximation_order`` defaults to the *slowest* rate any supported
    family has and a row may tighten it. For an isotropic Matérn-``nu`` in
    ``d`` axes the spectral density decays as ``omega ** -(2 nu + d)``, so the
    truncated tail — and with it the reduced-rank error — falls as
    ``m ** -2nu``: order 1 for Matérn-1/2, 3 for Matérn-3/2, 5 for Matérn-5/2,
    and exponentially for a squared exponential. The default is therefore
    ``1.0``.

    ``approximation_final`` is the separate, absolute claim: whatever the
    envelope allowed on the way, the finest setting a row uses must actually
    be close to the exact answer. A row that asserted only the envelope would
    pass on a solver that converged beautifully to the wrong number. It too is
    kernel-dependent — a Matérn-1/2 is the roughest process this method
    supports and a few hundred basis members is not many for it — so a row
    that sweeps a rough kernel states its own and says why.
    """
    if size <= coarsest_size:
        return float(coarsest_error)
    rate = tolerances.approximation_order if order is None else float(order)
    ratio = float(coarsest_size) / float(size)
    return float(coarsest_error) * ratio**rate + tolerances.approximation_floor


DEFAULT_TOLERANCES = Tolerances()


# ---------------------------------------------------------------------------
# Declarations a backend realises
# ---------------------------------------------------------------------------


class SolverKind(enum.StrEnum):
    """The GP solve strategies a row may ask a backend for."""

    DENSE = "dense"
    QUASISEP = "quasisep"
    #: W5.4's reduced-rank spectral solver. The first ``EXACT = False`` kind
    #: the battery carries, and the reason :class:`Tolerances` gained an
    #: approximation class: a row that asks for it compares against ``DENSE``
    #: inside an envelope that tightens with the basis size rather than at a
    #: fixed number. ``ConformanceBackend.gp_solver`` takes ``basis_size`` and
    #: ``boundary_factor`` for this kind and ignores them for the other two.
    HILBERT = "hilbert"


class KernelFamily(enum.StrEnum):
    """Kernel families a row may ask a backend for.

    Every one but ``SQUARED_EXPONENTIAL`` and ``PRODUCT`` has an exact
    semiseparable representation and so reaches the O(N) path; the two that do
    not are here precisely so the battery can hold the refusals (W4.5).
    """

    MATERN12 = "matern12"
    MATERN32 = "matern32"
    MATERN52 = "matern52"
    SHO = "sho"
    ROTATION = "rotation"
    SQUARED_EXPONENTIAL = "squared_exponential"
    SUM = "sum"
    PRODUCT = "product"
    SPECTRAL_MIXTURE = "spectral_mixture"
    #: **W5.7.** A non-stationary wrapper: ``a(x) k(w(x), w(x')) a(x')`` over
    #: the single term in ``CovarianceSpec.terms``, with the knot locations and
    #: knot values fixed in ``tests/conformance/backends/_kernels.py``. Like
    #: ``SUM``, it is built from ``ampere.core`` on every fixture, because a
    #: wrapping kernel has no arithmetic of its own and adopts its child's
    #: namespace, device and capability flags.
    WARPED = "warped"
    #: A kernel declared **outside** ampere, whose celerite representation a
    #: user registers with ``register_quasiseparable_term``. Every backend
    #: builds it from its own ``Matern12`` under a family name of its own, so
    #: one registration must carry it onto all three O(N) paths.
    USER = "conformance_user_term"


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
        ``Spectrum``, on the three axes :func:`complex_axes` fixes. Added at
        W2.4 slice 3, for the ``complex_gaussian``
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

    **W4.5 made the declaration recursive.** ``terms`` is non-empty for a
    :attr:`KernelFamily.SUM`, :attr:`KernelFamily.PRODUCT` or
    :attr:`KernelFamily.SPECTRAL_MIXTURE`, and a backend's ``kernel()``
    builds its children the same way it builds a leaf. ``period`` and
    ``quality`` belong to the oscillators; ``axes`` is the selector, ``None``
    meaning every axis, exactly as in ``ampere.core``.
    """

    family: KernelFamily = KernelFamily.MATERN32
    amplitude: float = 0.4
    length_scale: float = 2.0
    period: float = 1.5
    quality: float = 3.0
    delta_quality: float = 0.5
    fraction: float = 0.3
    axes: tuple[str, ...] | None = None
    length_scale_unit: Any = None
    terms: tuple[CovarianceSpec, ...] = ()


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
    ``interferometry``
        Whether :meth:`ConformanceBackend.interferometry` can return this
        backend's Fourier, closure-phase and smearing steps and its three
        source models. ``False`` by default, so a backend that has not written
        them is unaffected and the interferometric rows skip with a reason
        naming what is owed — the same shape as ``complex_models``, and for
        the same reason: one track's slice must not be a change to every other
        track's fixture. Added at W4.1; the native twins are W4.3's.
    ``astrometry``
        Whether :meth:`ConformanceBackend.astrometry` can return this
        backend's epoch-sampling step and reflex-orbit model. ``False`` by
        default, for the same reason ``interferometry`` is: one modality's
        slice must not be a change to every other track's fixture. Added at
        W4.9.
    ``image``
        Whether :meth:`ConformanceBackend.image` can return this backend's
        PSF-convolution step and the image-emitting source models it convolves.
        ``False`` by default, for the same reason ``interferometry`` and
        ``astrometry`` are: one modality's slice must not be a change to every
        other track's fixture. Added at **W5.5**, the item that made an
        ``Image`` an observation rather than only an intermediate.
    ``picklable``
        Whether a problem composed from this backend's pieces can be sent to a
        worker process (W3.1). ``True`` by default, because a backend whose
        arrays are plain numpy or plain torch tensors pickles and the
        process-pool executor is therefore available to it. jax declares
        ``False``: a jax array carries a ``jaxlib`` ``Device`` handle, which is
        process-local and has no pickle reduction, and jax itself warns that
        forking a jax process is likely to deadlock. That is a fact about the
        backend rather than a gap in it — throughput on a device backend is
        W3.1 slice 2's per-chunk ``vmap``, not a pool of processes — so the
        battery asserts the *refusal* for such a backend rather than skipping
        the row.
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
    interferometry: bool = False
    astrometry: bool = False
    image: bool = False
    picklable: bool = True
    solvers: frozenset[SolverKind] = frozenset({SolverKind.DENSE})
    tolerances: Tolerances = DEFAULT_TOLERANCES
    # W5.9 -- appended, not interleaved.
    #: Whether :meth:`ConformanceBackend.joint_gp_noise` can return this
    #: backend's ``JointGaussianProcessNoise`` -- one correlated process over
    #: several channels of one model on a shared grid. ``False`` by default,
    #: for the same reason ``interferometry``, ``astrometry`` and ``image``
    #: are: one item's slice must not be a change to every other track's
    #: fixture.
    joint_noise: bool = False


@dataclasses.dataclass(frozen=True)
class InterferometryPieces:
    """The interferometric classes one backend supplies (*W4.1*).

    A record rather than eleven protocol methods, because they arrive as one
    slice and are asked for as one: a backend either has an interferometric
    vocabulary or it has none. The names are the reference backend's, which is
    the convention every other twin follows (there is no step registry — a
    twin associates with its reference by class name in the backend's own
    module, the native surface it exposes, and its ``BACKEND`` declaration).

    The steps
    ---------
    ``fourier_sample``
        ``Image -> VisibilitySet``, built with
        ``from_observed(container, field_of_view=..., oversampling=...)``.
    ``closure_phase``
        ``VisibilitySet -> ClosurePhases``, three baselines to one angle.
    ``amplitude``
        ``VisibilitySet -> VisibilitySet``, the modulus.
    ``bandwidth_smearing``, ``time_smearing``
        The two averaging steps, which ask ``fourier_sample`` for the extra
        ``(u, v)`` samples they average over.

    The models
    ----------
    ``uniform_disc``, ``gaussian_source``, ``binary``
        Emit an ``Image`` on a negotiated ``(x, y)`` grid in mas, in Jy/sr.
    ``uniform_disc_visibilities``, ``gaussian_source_visibilities``,
    ``binary_visibilities``
        The same three sources emitting a ``VisibilitySet`` directly from
        their closed forms, in Jy. Both routes are supported and neither is
        privileged; the rows hold both to the oracles in ``oracles.py``.
    """

    fourier_sample: type
    closure_phase: type
    amplitude: type
    bandwidth_smearing: type
    time_smearing: type
    uniform_disc: type
    gaussian_source: type
    binary: type
    uniform_disc_visibilities: type
    gaussian_source_visibilities: type
    binary_visibilities: type


@dataclasses.dataclass(frozen=True)
class AstrometryPieces:
    """The astrometric classes one backend supplies (*W4.9*).

    A two-field record, the same shape as :class:`InterferometryPieces` and
    for the same reason: there is no step registry, so association with the
    reference implementation is by class name in the backend's own
    ``astrometry`` module, the native surface it exposes, and its ``BACKEND``
    declaration.

    ``epoch_sample``
        ``TimeSeries -> TimeSeries``, kind-preserving, built with
        ``from_observed(container)``.
    ``reflex_orbit``
        Emits the two ``TimeSeries`` channels ``"ra"`` and ``"dec"``.
    """

    epoch_sample: type
    reflex_orbit: type


@dataclasses.dataclass(frozen=True)
class ImagePieces:
    """The gridded-image classes one backend supplies (*W5.5*).

    The same shape as :class:`InterferometryPieces` and :class:`AstrometryPieces`,
    and for the same reason: there is no step registry, so association with the
    reference implementation is by class name in the backend's own ``image``
    module, the native surface it exposes, and its ``BACKEND`` declaration.

    ``psf_convolution``
        ``Image -> Image``, kind-preserving and grid-cropping, built with
        ``from_observed(container, fwhm=...)`` or
        ``from_observed(container, kernel=...)``.
    ``gaussian_source``, ``binary``
        Two of W4.1's three image-emitting source models, reused unchanged as
        the things a PSF is convolved *with*. Band-limited, both of them, which
        is what makes a convolution of them comparable against a direct sum to
        the solver tolerance rather than to a convergence rate — the
        uniform disc's sharp edge is the exception W4.1's own rows record, and
        it is deliberately not in this record.
    """

    psf_convolution: type
    gaussian_source: type
    binary: type


@runtime_checkable
class ConformanceBackend(Protocol):
    """Everything the conformance battery asks of a backend.

    Thirteen members (W5.5 added :meth:`image`; W4.9 added :meth:`astrometry`,
    and W4.1 :meth:`interferometry` before it). Implement them and every row in
    ``tests/conformance/``
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

    def gp_solver(
        self,
        kind: SolverKind,
        *,
        basis_size: int | Sequence[int] = 32,
        boundary_factor: float = 2.0,
    ) -> GPSolver:
        """Return this backend's solver for *kind*.

        Only called for a kind present in ``capabilities.solvers``; a backend
        may raise for anything else.

        **W5.4** adds the two keywords, which describe the *approximation* and
        so mean nothing to the two exact kinds: a backend ignores them unless
        *kind* is :attr:`SolverKind.HILBERT`. They are on this one method
        rather than on a second one because the convergence rows ask for the
        same solver at four basis sizes, and a row that had to know which
        method to call for which kind would be naming solvers rather than
        kinds.
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

    def gp_noise(self, kernel: Kernel, solver: GPSolver, *, jitter: Any = None) -> NoiseModel:
        """This backend's GP noise composition over *kernel* and *solver*.

        **Added at W2.13**, for the reason :meth:`independent_noise` gives.
        *kernel* and *solver* are this fixture's own, from :meth:`kernel` and
        :meth:`gp_solver`.

        *jitter* — **added at W3.1 slice 2** — is the optional diagonal floor
        every ``GaussianProcessNoise`` takes at construction, and it is a
        keyword here because a backend that had quietly dropped it would be a
        backend on which ``sigma_eff² = (scale·sigma_data)² + jitter²`` means
        something else. ``None`` (the default) registers no such parameter,
        which is what every row written before this keyword existed asks for.
        """

    def joint_gp_noise(
        self,
        kernel: Kernel,
        solver: GPSolver,
        *,
        datasets: Sequence[str],
        coupling: Any,
    ) -> NoiseModel:
        """This backend's joint noise model over *datasets* (**W5.9**).

        Only called when :attr:`BackendCapabilities.joint_noise` is declared.
        *kernel* and *solver* are this fixture's own, from :meth:`kernel` and
        :meth:`gp_solver`; *coupling* is ``ampere.core``'s, because a
        :class:`~ampere.core.ChannelCoupling` declares parameters and an
        eigendecomposition written once for every array namespace and so has no
        per-backend twin to supply.
        """

    def parameter_space(self, declaration: ParameterSet) -> ParameterSpace:
        """Realise a parameter declaration as this backend's engine-facing view."""

    def interferometry(self) -> InterferometryPieces:
        """This backend's interferometric steps and source models (*W4.1*).

        Only called when :attr:`BackendCapabilities.interferometry` is
        declared; a backend that has not written them may raise, and every
        interferometric row then skips with a reason naming what is owed. The
        two container kinds these speak in are ``ampere.core``'s, not a
        backend's — the placement D1 ruled on 2026-09-11 — so a backend
        supplies only the arithmetic, which is the same claim
        ``inference.md`` §18 makes for every other piece here.
        """

    def astrometry(self) -> AstrometryPieces:
        """This backend's epoch-sampling step and reflex-orbit model (*W4.9*).

        Only called when :attr:`BackendCapabilities.astrometry` is declared;
        a backend that has not written them may raise, and every astrometric
        row then skips with a reason naming what is owed. The kind these speak
        in, :class:`~ampere.core.TimeSeries`, is ``ampere.core``'s own — the
        same placement D1 ruled for interferometry, applied to the modality
        that was chosen to test the ruling.
        """

    def image(self) -> ImagePieces:
        """This backend's PSF-convolution step and its image sources (*W5.5*).

        Only called when :attr:`BackendCapabilities.image` is declared; a
        backend that has not written them may raise, and every gridded row then
        skips with a reason naming what is owed. The kind these speak in,
        :class:`~ampere.core.Image`, is ``ampere.core``'s own — D1's placement
        again, applied to the first modality whose observed container has a
        ``Layout.GRID``.
        """

    def to_numpy(self, values: Any) -> np.ndarray:
        """Bring a backend array back to numpy, for comparison against oracles.

        The battery's oracles are ``scipy`` and closed forms, which speak
        numpy. A backend working in torch or jax arrays converts here — and
        only here, so no row has to know.
        """


#: Wavelength every :attr:`ModelKind.COMPLEX` sample is measured at, micron.
#: A single number because the battery's complex rows are monochromatic: the
#: spectral axis a ``VisibilitySet`` has carried since W4.1 is a constant
#: column here, which is the amendment's own statement of the monochromatic
#: case and keeps the rows about the family's arithmetic.
COMPLEX_WAVELENGTH = 2.2


def complex_axes(coordinates: Any) -> tuple[np.ndarray, np.ndarray, Any]:
    """The ``(u, v, spectral_axis)`` point set a :attr:`ModelKind.COMPLEX` channel lives on.

    A :class:`~ampere.core.VisibilitySet` has three coordinate axes and a
    ``ModelSpec`` declares one sequence of coordinates, so the other two have
    to come from somewhere. They come from here rather than from either side,
    because ``Likelihood.check_alignment`` compares the predicted and observed
    axes for equality: the fixture's model and the battery's observed container
    must derive them by the *same* rule, and a rule written twice is a rule
    that eventually differs.

    ``v = u / 2`` is arbitrary and deliberately so — the battery's complex rows
    are about the family's arithmetic, not about (u,v) coverage — but it is
    distinct from ``u``, which keeps the two axes from being accidentally
    interchangeable in a comparison. The wavelength is
    :data:`COMPLEX_WAVELENGTH` on every sample, a :class:`~astropy.units.Quantity`
    so that the axis arrives with its unit rather than as a bare array the
    kind would refuse.
    """
    u_axis = np.asarray(coordinates, dtype=float)
    wavelength = np.full(u_axis.shape, COMPLEX_WAVELENGTH) * astropy_units.micron
    return u_axis, 0.5 * u_axis, wavelength


def solver_kinds(capabilities: BackendCapabilities) -> tuple[SolverKind, ...]:
    """The solver kinds *capabilities* declares, in :class:`SolverKind` order."""
    return tuple(kind for kind in SolverKind if kind in capabilities.solvers)


def as_sequence(values: Any) -> Sequence[float]:
    """Coerce *values* to a plain float sequence, for spec fields."""
    return tuple(float(value) for value in np.asarray(values, dtype=float).ravel())
