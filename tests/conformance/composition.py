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
from collections.abc import Mapping
from typing import Any

import astropy.units as u
import numpy as np

from ampere.core import (
    Censoring,
    Dataset,
    FittingProblem,
    GaussianProcessNoise,
    IndependentNoise,
    Instrument,
    Likelihood,
    LimitKind,
    NoiseModel,
    Spectrum,
    Tie,
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
)

__all__ = [
    "COARSE_GRID",
    "GP_GRID",
    "CensoringKind",
    "DatasetSpec",
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

    The observed container is a ``Spectrum`` in every case: the battery's
    per-kind coverage lives in ``test_schema.py``, which needs no problem, and
    keeping one kind here keeps ``build_problem`` honest about alignment.

    ``censoring`` may only be combined with :attr:`NoiseKind.IID`. Censored
    data under a correlated noise model is a multivariate-normal orthant
    probability with no closed form, and ``GaussianFamily`` refuses it by
    design (``likelihoods.md`` §9); a spec that asks for both is a bug in the
    row, not a backend failure.
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
    return float(values["norm"]) * (grid / spec.reference_coordinate) ** float(values["index"])


def truth(spec: ModelSpec) -> dict[str, float]:
    """A well-inside-the-prior parameter vector for *spec*'s model."""
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


def observed_values(spec: ProblemSpec, dataset: DatasetSpec, grid: np.ndarray) -> np.ndarray:
    """Fixed, deterministic pseudo-data on *grid*.

    Built from :func:`analytic_flux` at :func:`truth` plus a seeded draw, so
    the residuals are the size of the uncertainties and every row works on a
    well-posed problem. Crucially it goes through the *oracle*, never through
    a backend's model: the data are then identical in every backend, which is
    what makes ``ampere_data_hash`` comparable across them and keeps the data
    from being a function of the arithmetic under test.
    """
    rng = np.random.default_rng(seed=dataset.data_seed)
    mean = analytic_flux(spec.model, truth(spec.model), grid)
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


def observed_container(spec: ProblemSpec, dataset: DatasetSpec) -> Spectrum:
    """The observed ``Spectrum`` for *dataset*, identical in every backend."""
    grid = observed_grid(spec, dataset)
    values = observed_values(spec, dataset, grid)
    return Spectrum(
        grid * COORDINATE_UNIT,
        values * FLUX_UNIT,
        uncertainty=np.full(grid.size, dataset.uncertainty) * FLUX_UNIT,
        mask=observed_mask(dataset, grid.size),
    )


# ---------------------------------------------------------------------------
# Composition
# ---------------------------------------------------------------------------


def build_noise(backend: ConformanceBackend, dataset: DatasetSpec) -> NoiseModel:
    """The noise model *dataset* declares, with the backend's kernel and solver."""
    if dataset.noise is NoiseKind.IID:
        return IndependentNoise()
    return GaussianProcessNoise(
        backend.kernel(dataset.covariance), backend.gp_solver(dataset.solver)
    )


def build_likelihood(
    backend: ConformanceBackend, dataset: DatasetSpec, observed: Spectrum
) -> Likelihood:
    """The likelihood *dataset* declares."""
    codes = censoring_codes(dataset, observed.n_samples)
    censoring = None if codes is None else Censoring(codes)
    return Likelihood(
        family_named(dataset.family)(), build_noise(backend, dataset), censoring=censoring
    )


def build_instrument(backend: ConformanceBackend, dataset: DatasetSpec) -> Instrument:
    """The instrument chain *dataset* declares, from the backend's steps."""
    return Instrument(
        tuple(backend.transformation(step) for step in dataset.instrument),
        channel=dataset.channel,
        input_kind=Spectrum,
        label=dataset.label,
    )


def build_problem(backend: ConformanceBackend, spec: ProblemSpec) -> FittingProblem:
    """Realise *spec* on *backend*.

    The only backend-supplied pieces are the model, the instrument steps, the
    kernel and the solver. Everything else is ``ampere.core``, unchanged.
    """
    datasets = []
    for dataset in spec.datasets:
        observed = observed_container(spec, dataset)
        datasets.append(
            Dataset(
                observed,
                build_instrument(backend, dataset),
                build_likelihood(backend, dataset, observed),
                label=dataset.label,
            )
        )
    return FittingProblem(backend.model(spec.model), datasets, ties=spec.ties, seed=spec.seed)
