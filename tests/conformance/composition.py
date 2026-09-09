"""Turning a backend-neutral :class:`ProblemSpec` into a ``FittingProblem``.

A conformance row never builds a problem by hand: it declares one and asks the
suite to realise it on whichever backend the fixture supplied. That is what
makes the cross-backend rows possible at all — two backends handed the *same*
declaration must produce the same numbers — and it is what keeps a Phase-2
author's obligation down to :class:`~tests.conformance.protocol.ConformanceBackend`'s
seven members.

Everything the backend does not supply is assembled here, from
``ampere.core``: ``Instrument``, ``Likelihood``, ``Dataset`` and
``FittingProblem`` are reused unchanged in every backend, which is the claim
``inference.md`` §18 makes for them.

The module also carries the **analytic oracles** — :func:`analytic_flux` is the
closed form :class:`~tests.conformance.protocol.ModelKind` defines, written out
from the definition rather than lifted from any implementation, so a row
comparing against it is comparing against mathematics and not against a second
copy of ampere.
"""

from __future__ import annotations

import dataclasses
import enum
from collections.abc import Callable, Mapping
from typing import Any

import astropy.units as u
import numpy as np

from ampere.core import (
    Censoring,
    Dataset,
    FittingProblem,
    FunctionSamples,
    Instrument,
    Kernel,
    Likelihood,
    LimitKind,
    NoiseModel,
    PhotometricPoints,
    Spectrum,
    Tie,
    VisibilitySet,
    family_named,
)

from .protocol import (
    ConformanceBackend,
    CovarianceSpec,
    ModelKind,
    ModelSpec,
    SolverKind,
    TransformationKind,
    TransformationSpec,
    complex_axes,
)

__all__ = [
    "COARSE_GRID",
    "COUNT_EXPOSURE",
    "GP_GRID",
    "CensoringKind",
    "DatasetSpec",
    "KernelFactory",
    "NoiseKind",
    "ProblemSpec",
    "analytic_flux",
    "build_problem",
    "censoring_codes",
    "model_context",
    "observed_container",
    "observed_grid",
    "observed_mask",
    "observed_values",
    "terminal_photometry",
    "truth",
]

FLUX_UNIT = u.Jy
COORDINATE_UNIT = u.micron

#: A finer grid than :class:`ModelSpec`'s default, for the GP rows: a
#: correlated likelihood on six points is a weak test of a covariance.
GP_GRID: tuple[float, ...] = (
    1.0,
    1.4,
    1.9,
    2.5,
    3.2,
    4.0,
    5.0,
    6.1,
    7.3,
    8.6,
    10.0,
    11.5,
)

#: A coarser grid for the resampling rows, strictly inside :data:`GP_GRID`.
COARSE_GRID: tuple[float, ...] = (1.5, 3.5, 6.5, 9.5)


class NoiseKind(enum.StrEnum):
    """Which noise model a dataset composes."""

    IID = "iid"
    GP = "gp"


class CensoringKind(enum.StrEnum):
    """Which Tobit limit a censored sample carries."""

    NONE = "none"
    UPPER = "upper"
    LOWER = "lower"


@dataclasses.dataclass(frozen=True)
class DatasetSpec:
    """One dataset's declaration.

    The observed container is a ``Spectrum`` unless the instrument chain ends
    in a ``PHOTOMETRY`` step, in which case it is a ``PhotometricPoints``,
    because that step changes kind and ``Likelihood.check_alignment`` compares
    like with like. Those are the only two: the battery's per-kind coverage
    lives in ``test_schema.py``, which needs no problem.

    ``censoring`` may only be combined with :attr:`NoiseKind.IID`. Censored
    data under a correlated noise model is a multivariate-normal orthant
    probability with no closed form, and ``GaussianFamily`` refuses it by
    design (``likelihoods.md`` §9); a spec that asks for both is a bug in the
    row, not a backend failure.

    ``family`` also decides the *shape of the data*, and one family changes it:
    ``"poisson"`` needs non-negative integer counts and no uncertainties at
    all (``PoissonFamily`` defines its own dispersion, and
    ``REQUIRES_UNCERTAINTY`` is ``False``). That is derived from the family
    rather than declared separately so the two cannot disagree — a Poisson
    spec carrying Gaussian uncertainties would be refused at composition, and
    a flag able to express it buys nothing. Combined with
    :attr:`NoiseKind.GP` it is the **latent-GP** shape (W2.14): the dataset
    declares one whitened latent value per retained sample, and the family
    scores at ``f = L(θ) z``.
    """

    label: str = "sed"
    channel: str = "default"
    instrument: tuple[TransformationSpec, ...] = ()
    grid: tuple[float, ...] | None = None
    uncertainty: float = 0.1
    masked: tuple[int, ...] = ()
    family: str = "gaussian"
    noise: NoiseKind = NoiseKind.IID
    covariance: CovarianceSpec = dataclasses.field(default_factory=CovarianceSpec)
    solver: SolverKind = SolverKind.DENSE
    gp_jitter: float | None = None
    """A fixed diagonal floor on a :attr:`NoiseKind.GP` noise model (W3.1 slice 2).

    ``None`` registers no such parameter, which is what every row written
    before this field existed asks for. A number registers a held-fixed
    ``jitter``, so ``sigma_eff² = sigma_data² + jitter²`` — the keyword the
    torch class has taken since W2.4 and the jax class had not, which is why
    it is asserted here rather than in one backend's own suite.
    """

    censoring: CensoringKind = CensoringKind.NONE
    censored: tuple[int, ...] = ()
    data_seed: int = 20260902


@dataclasses.dataclass(frozen=True)
class ProblemSpec:
    """A whole fitting problem, declared without naming a backend."""

    model: ModelSpec = dataclasses.field(default_factory=ModelSpec)
    datasets: tuple[DatasetSpec, ...] = dataclasses.field(default_factory=lambda: (DatasetSpec(),))
    ties: tuple[Tie, ...] = ()
    seed: int = 20260902


# ---------------------------------------------------------------------------
# Analytic oracles and deterministic data
# ---------------------------------------------------------------------------


def analytic_flux(spec: ModelSpec, values: Mapping[str, Any], grid: np.ndarray) -> np.ndarray:
    """The closed form :class:`ModelKind` declares, written out from scratch.

    This is an oracle, not a reimplementation: it is the definition in
    :class:`ModelKind`'s docstring transcribed into numpy. Rows that need the
    predicted mean — the GP covariance row, the Tobit row — use it so that
    nothing in the comparison came from the object under test.
    """
    if spec.kind is ModelKind.LINEAR:
        return float(values["offset"]) + float(values["slope"]) * grid
    if spec.kind is ModelKind.COMPLEX:
        # V(x) = norm * exp(i * index * x): constant modulus, winding phase.
        # Written this way on purpose -- ``index`` moves the *phase only*, so a
        # backend that quietly dropped the imaginary part would compute a
        # density independent of it and disagree with this oracle at every
        # point rather than in the last digit.
        return float(values["norm"]) * np.exp(1j * float(values["index"]) * grid)
    return float(values["norm"]) * (grid / spec.reference_coordinate) ** float(values["index"])


def truth(spec: ModelSpec) -> dict[str, float]:
    """A well-inside-the-prior parameter vector for *spec*'s model.

    ``COMPLEX`` shares ``POWER_LAW``'s, because it shares its parameter names
    and priors — see :class:`~tests.conformance.protocol.ModelKind`.
    """
    if spec.kind is ModelKind.LINEAR:
        return {"offset": 0.2, "slope": 1.0}
    return {"norm": 1.5, "index": -1.0}


def model_context(problem: FittingProblem, values: Any = None) -> dict[str, Any]:
    """The model's own parameter values, de-qualified, from a merged vector.

    A row that wants to call :func:`analytic_flux` needs ``{"slope": …}``, not
    ``{"model.slope": …}``. Going through ``distribute`` rather than string
    surgery means the helper keeps working when the composition changes.
    """
    resolved = (
        dict(problem.reference_values)
        if values is None
        else dict(problem.parameters.complete(values))
    )
    routed = problem.mapping.distribute(resolved)
    return dict(routed["model"])


def observed_grid(spec: ProblemSpec, dataset: DatasetSpec) -> np.ndarray:
    """The coordinates the observed container lives on.

    An explicit ``grid`` wins; otherwise a resampling step's target grid does
    (the instrument chain moves the prediction onto it, and the data must be
    there to meet it); otherwise the model's own coordinates.
    """
    if dataset.grid is not None:
        return np.asarray(dataset.grid, dtype=float)
    for step in reversed(dataset.instrument):
        if step.kind in (TransformationKind.REBIN, TransformationKind.PHOTOMETRY):
            return np.asarray(step.target, dtype=float)
    return np.asarray(spec.model.coordinates, dtype=float)


#: How many counts a unit of the oracle's flux is worth, for the Poisson
#: shapes. Chosen so the smallest count on ``GP_GRID`` is still comfortably
#: above one — a shape whose counts were mostly zero would make every row that
#: uses it a test of the ``k = 0`` term rather than of the composition.
COUNT_EXPOSURE = 20.0


def observed_values(spec: ProblemSpec, dataset: DatasetSpec, grid: np.ndarray) -> np.ndarray:
    """Fixed, deterministic pseudo-data on *grid*.

    Built from :func:`analytic_flux` at :func:`truth` plus a seeded draw, so
    the residuals are the size of the uncertainties and every row works on a
    well-posed problem. Crucially it goes through the *oracle*, never through
    a backend's model: the data are then identical in every backend, which is
    what makes ``ampere_data_hash`` comparable across them and keeps the data
    from being a function of the arithmetic under test.

    A Poisson dataset gets **counts** instead: the same oracle mean, scaled by
    :data:`COUNT_EXPOSURE` and drawn from a Poisson. Still seeded, still
    computed here rather than in a backend, so the two claims above hold
    unchanged.
    """
    rng = np.random.default_rng(seed=dataset.data_seed)
    mean = analytic_flux(spec.model, truth(spec.model), grid)
    if dataset.family == "poisson":
        return rng.poisson(np.clip(mean, 1e-6, None) * COUNT_EXPOSURE).astype(float)
    if dataset.family == "complex_gaussian":
        # The circular complex Gaussian is exactly this draw: independent
        # Normal(0, sigma) on each component, one real ``uncertainty`` for
        # both (``results_schema.md`` §16). Generating the data the way the
        # family defines them keeps the residuals the size of the error bars
        # here as everywhere else.
        real = rng.normal(0.0, dataset.uncertainty, grid.size)
        imaginary = rng.normal(0.0, dataset.uncertainty, grid.size)
        return mean + real + 1j * imaginary
    return mean + rng.normal(0.0, dataset.uncertainty, grid.size)


def observed_mask(dataset: DatasetSpec, size: int) -> np.ndarray | None:
    """``True`` where a sample is excluded, or ``None`` for an unmasked set."""
    if not dataset.masked:
        return None
    mask = np.zeros(size, dtype=bool)
    mask[list(dataset.masked)] = True
    return mask


def censoring_codes(dataset: DatasetSpec, size: int) -> np.ndarray | None:
    """The Tobit limit codes for *dataset*, or ``None`` if nothing is censored."""
    if dataset.censoring is CensoringKind.NONE or not dataset.censored:
        return None
    kind = (
        LimitKind.UPPER_LIMIT if dataset.censoring is CensoringKind.UPPER else LimitKind.LOWER_LIMIT
    )
    codes = np.zeros(size, dtype=np.int8)
    codes[list(dataset.censored)] = int(kind)
    return codes


def terminal_photometry(dataset: DatasetSpec) -> TransformationSpec | None:
    """The photometry step a chain *ends* with, or ``None``.

    Which container kind the observed data must be is decided by the last
    kind-changing step, and ``PHOTOMETRY`` is the only one the battery
    declares: a chain ending in it predicts ``PhotometricPoints``, and
    ``Likelihood.check_alignment`` compares like with like, so the observed
    side has to be one too.
    """
    for step in reversed(dataset.instrument):
        if step.kind is TransformationKind.PHOTOMETRY:
            return step
        if step.kind is TransformationKind.REBIN:
            return None
    return None


def observed_container(spec: ProblemSpec, dataset: DatasetSpec) -> FunctionSamples:
    """The observed container for *dataset*, identical in every backend.

    A ``Spectrum`` in every case but one: a chain ending in ``PHOTOMETRY``
    produces :class:`~ampere.core.PhotometricPoints`, so the data must be
    those too (W2.5 slice 2, adding the battery's first end-to-end photometry
    shape). The values, the seed and the uncertainties are the same function
    of the grid either way, so the two shapes stay comparable and
    ``ampere_data_hash`` still depends on the declaration rather than on the
    arithmetic under test.
    """
    grid = observed_grid(spec, dataset)
    values = observed_values(spec, dataset, grid)
    # Counts carry no uncertainties: the family supplies its own dispersion,
    # and attaching a Gaussian sigma to them would be a statement about the
    # data that is not true.
    uncertainty = (
        None if dataset.family == "poisson" else np.full(grid.size, dataset.uncertainty) * FLUX_UNIT
    )
    mask = observed_mask(dataset, grid.size)
    if np.iscomplexobj(values):
        # W2.4 slice 3. ``VisibilitySet`` is the only container kind that takes
        # complex values (``results_schema.md`` §16), and its uncertainty is
        # real: the per-component standard deviation of the circular complex
        # Gaussian, which is the noise model the family is.
        u_axis, v_axis = complex_axes(grid)
        return VisibilitySet(
            u_axis,
            v_axis,
            values * FLUX_UNIT,
            uncertainty=uncertainty,
            mask=mask,
        )
    photometry = terminal_photometry(dataset)
    if photometry is not None:
        return PhotometricPoints(
            photometry.filters,
            grid * COORDINATE_UNIT,
            values * FLUX_UNIT,
            uncertainty=uncertainty,
            mask=mask,
        )
    return Spectrum(
        grid * COORDINATE_UNIT,
        values * FLUX_UNIT,
        uncertainty=uncertainty,
        mask=mask,
    )


# ---------------------------------------------------------------------------
# Composition
# ---------------------------------------------------------------------------


#: How a dataset's kernel is built. ``ConformanceBackend.kernel`` for every row
#: but W3.8's, which passes ``ampere.core``'s own constructors to compose the
#: **foreign** kernel a native problem may opt into.
KernelFactory = Callable[[CovarianceSpec], Kernel]


def build_noise(
    backend: ConformanceBackend,
    dataset: DatasetSpec,
    *,
    kernel_factory: KernelFactory | None = None,
) -> NoiseModel:
    """The noise model *dataset* declares, from the backend's own noise classes.

    The noise model and the solver were ``ampere.core``'s here until W2.13,
    when the capability flags widened to cover both (``inference.md`` §10a,
    fold-in 7). They carry a ``BACKEND`` now, and
    ``Likelihood.capability_parts`` puts them on the composed problem, so a
    fixture composing the core classes under its own name would declare two
    backends and be refused — which is the point of the widening, not a
    casualty of it.
    """
    if dataset.noise is NoiseKind.IID:
        return backend.independent_noise()
    build_kernel = backend.kernel if kernel_factory is None else kernel_factory
    return backend.gp_noise(
        build_kernel(dataset.covariance),
        backend.gp_solver(dataset.solver),
        jitter=dataset.gp_jitter,
    )


def build_likelihood(
    backend: ConformanceBackend,
    dataset: DatasetSpec,
    observed: FunctionSamples,
    *,
    kernel_factory: KernelFactory | None = None,
) -> Likelihood:
    """The likelihood *dataset* declares."""
    codes = censoring_codes(dataset, observed.n_samples)
    censoring = None if codes is None else Censoring(codes)
    return Likelihood(
        family_named(dataset.family)(),
        build_noise(backend, dataset, kernel_factory=kernel_factory),
        censoring=censoring,
    )


def build_instrument(backend: ConformanceBackend, dataset: DatasetSpec) -> Instrument:
    """The instrument chain *dataset* declares, from the backend's steps.

    ``input_kind`` is the kind the model emits, so a complex dataset declares
    :class:`~ampere.core.VisibilitySet` rather than ``Spectrum``. The battery
    declares no step that accepts one — the instrument vocabulary is
    spectrum-shaped, and Phase 4's is where a visibility chain belongs — so a
    complex dataset's chain is always empty and the declaration is only there
    to keep ``Instrument``'s own kind check truthful.
    """
    return Instrument(
        tuple(backend.transformation(step) for step in dataset.instrument),
        channel=dataset.channel,
        input_kind=VisibilitySet if dataset.family == "complex_gaussian" else Spectrum,
        label=dataset.label,
    )


def build_problem(
    backend: ConformanceBackend,
    spec: ProblemSpec,
    *,
    kernel_factory: KernelFactory | None = None,
    **problem_kwargs: Any,
) -> FittingProblem:
    """Realise *spec* on *backend*.

    The only backend-supplied pieces are the model, the instrument steps, the
    kernel and the solver. Everything else is ``ampere.core``, unchanged.

    ``kernel_factory`` replaces :meth:`ConformanceBackend.kernel`, and exists
    for W3.8's rows alone: they compose the *same* declaration with a kernel
    from another backend, which is the composition the ruling refuses by
    default and accepts under ``allow_foreign_parts=True``. ``problem_kwargs``
    reach :class:`~ampere.core.FittingProblem` unchanged, which is how those
    rows pass that flag.
    """
    datasets = []
    for dataset in spec.datasets:
        observed = observed_container(spec, dataset)
        datasets.append(
            Dataset(
                observed,
                build_instrument(backend, dataset),
                build_likelihood(backend, dataset, observed, kernel_factory=kernel_factory),
                label=dataset.label,
            )
        )
    return FittingProblem(
        backend.model(spec.model),
        datasets,
        ties=spec.ties,
        seed=spec.seed,
        **problem_kwargs,
    )
