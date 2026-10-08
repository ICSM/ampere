"""The portable model (W7.13): one ``_flux``, three backends, the same numbers.

``ModelKind.PORTABLE_LINEAR`` is the guide's model —
``examples/portable_model/model.py``, written once against
:class:`ampere.core.PortableModel` — built by each fixture as its own
backend's one-line twin. Every other kind is the fixture's own transcription;
this one is not, and the claim these rows test is the item's: the same source
runs on numpy, jax and torch and agrees

* with a numpy oracle written here, never read from the object, at
  ``tolerances.exact`` (rows (a));
* with itself across every pair of fixtures at ``tolerances.cross_backend``
  (rows (b));
* on its native surface with its contract surface (rows (c));
* under NUTS on each fixture with a differentiable path (rows (d));
* in the capability flags it declares (rows (e)).

A fixture that builds no ``PORTABLE_LINEAR`` model (the in-repo mirror, whose
models are deliberately its own transcriptions) skips with that reason. No row
names a concrete backend.
"""

from __future__ import annotations

import warnings
from typing import Any

import astropy.units as u
import numpy as np
import pytest

from ampere.core import (
    Dataset,
    FittingProblem,
    Model,
    Spectrum,
    negotiate,
    realise,
    registered_realisations,
)

from .composition import GP_GRID, DatasetSpec, build_instrument, build_likelihood
from .conftest import cross_backend_pairs, pair_ids
from .protocol import (
    ConformanceBackend,
    ModelKind,
    ModelSpec,
    Tolerances,
    TransformationKind,
    TransformationSpec,
)

SPEC = ModelSpec(kind=ModelKind.PORTABLE_LINEAR, coordinates=GP_GRID)
COMPILED = ModelSpec(kind=ModelKind.PORTABLE_LINEAR, coordinates=GP_GRID, compiled=True)

#: Three parameter points: inside the flat priors on [-10, 10], one negative
#: slope, one near the boundary.
POINTS: tuple[dict[str, float], ...] = (
    {"slope": 1.0, "intercept": 0.2},
    {"slope": -2.5, "intercept": 3.75},
    {"slope": 9.5, "intercept": -9.0},
)

#: The NUTS row's injected truth and its data's noise.
TRUTH = {"slope": 1.0, "intercept": 0.2}
SIGMA = 0.1

#: A resampling target strictly inside ``GP_GRID``, for the rows that need a
#: negotiated grid and an instrument of the fixture's own.
TARGET: tuple[float, ...] = (1.2, 2.2, 3.6, 5.5, 7.9, 10.8)
REBIN = DatasetSpec(
    instrument=(TransformationSpec(kind=TransformationKind.REBIN, target=TARGET),),
)


def oracle(values: dict[str, float], grid: Any) -> np.ndarray:
    """``slope * x + intercept``, transcribed from the kind's docstring."""
    return values["slope"] * np.asarray(grid, dtype=float) + values["intercept"]


def qualified(values: dict[str, float]) -> dict[str, float]:
    """*values* under the merged names a problem reads (the model component is ``model``)."""
    return {f"model.{name}": value for name, value in values.items()}


def portable(backend: ConformanceBackend, spec: ModelSpec = SPEC) -> Model:
    """*backend*'s ``PORTABLE_LINEAR`` model, or a skip naming why it has none."""
    try:
        return backend.model(spec)
    except KeyError:
        pytest.skip(
            f"the {backend.name!r} fixture builds no {ModelKind.PORTABLE_LINEAR} model "
            "(its models are its own transcriptions; this kind's subject is one shared source)"
        )


def flux(backend: ConformanceBackend, model: Model, values: dict[str, float]) -> np.ndarray:
    """The contract surface's single channel, as numpy."""
    return backend.to_numpy(model(**values).single().values)


def rebinned_problem(backend: ConformanceBackend, *, seed: int = 20261008) -> FittingProblem:
    """A problem on the fixture's own resampling step and independent noise.

    The data are this fixture's own instrument applied to the model at
    :data:`TRUTH`, plus a seeded draw: the row's subject is recovery through
    the native path, not the resampling operator's definition.
    """
    model = portable(backend, COMPILED)
    instrument = build_instrument(backend, REBIN)
    model.compile_for(negotiate([instrument]))
    clean = instrument(model(**TRUTH))
    rng = np.random.default_rng(seed)
    values = np.asarray(clean.values, dtype=float)
    grid = np.asarray(clean.spectral_axis.values, dtype=float)
    observed = Spectrum(
        grid * u.micron,
        (values + rng.normal(0.0, SIGMA, values.size)) * u.Jy,
        uncertainty=np.full(values.size, SIGMA) * u.Jy,
    )
    likelihood = build_likelihood(backend, REBIN, observed)
    return FittingProblem(
        portable(backend, COMPILED),
        [Dataset(observed, build_instrument(backend, REBIN), likelihood, label="sed")],
        seed=seed,
    )


# ---------------------------------------------------------------------------
# (a) the numpy oracle, per fixture
# ---------------------------------------------------------------------------


class TestAgainstTheOracle:
    @pytest.mark.parametrize("values", POINTS, ids=["truth", "negative", "edge"])
    def test_the_flux_is_the_closed_form(
        self, backend: ConformanceBackend, tolerances: Tolerances, values: dict[str, float]
    ) -> None:
        model = portable(backend)
        np.testing.assert_allclose(
            flux(backend, model, values),
            oracle(values, GP_GRID),
            rtol=tolerances.exact,
            atol=tolerances.exact,
        )

    def test_a_negotiated_grid_is_adopted(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        model = portable(backend, COMPILED)
        model.compile_for(negotiate([build_instrument(backend, REBIN)]))
        result = model(**TRUTH).single()
        grid = np.asarray(result.spectral_axis.values, dtype=float)
        assert grid.size >= len(TARGET)
        np.testing.assert_allclose(
            backend.to_numpy(result.values),
            oracle(TRUTH, grid),
            rtol=tolerances.exact,
            atol=tolerances.exact,
        )


# ---------------------------------------------------------------------------
# (b) every pair of fixtures
# ---------------------------------------------------------------------------


def shared_tolerance(first: ConformanceBackend, second: ConformanceBackend) -> float:
    return max(
        first.capabilities.tolerances.cross_backend, second.capabilities.tolerances.cross_backend
    )


@pytest.mark.parametrize(("first", "second"), cross_backend_pairs(), ids=pair_ids())
class TestAcrossBackends:
    def test_the_flux_agrees(self, first: ConformanceBackend, second: ConformanceBackend) -> None:
        left, right = portable(first), portable(second)
        for values in POINTS:
            np.testing.assert_allclose(
                flux(first, left, values),
                flux(second, right, values),
                rtol=0.0,
                atol=shared_tolerance(first, second),
            )

    def test_the_composed_density_agrees(
        self, first: ConformanceBackend, second: ConformanceBackend
    ) -> None:
        left, right = rebinned_problem(first), rebinned_problem(second)
        # Each fixture's data came through its own instrument; compare the
        # densities on one set of data so the row tests the model alone.
        right_on_left = FittingProblem(
            portable(second, COMPILED),
            [
                Dataset(
                    left.datasets["sed"].observed,
                    build_instrument(second, REBIN),
                    build_likelihood(second, REBIN, left.datasets["sed"].observed),
                    label="sed",
                )
            ],
            seed=20261008,
        )
        assert right.free_labels() == left.free_labels()
        for values in POINTS:
            merged = qualified(values)
            assert right_on_left.log_prob(merged) == pytest.approx(
                left.log_prob(merged), abs=shared_tolerance(first, second)
            )


# ---------------------------------------------------------------------------
# (c) the native surface against the contract surface
# ---------------------------------------------------------------------------


class TestTheNativeSurface:
    def test_native_flux_and_grid_agree_with_evaluate(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        model = portable(backend, COMPILED)
        model.compile_for(negotiate([build_instrument(backend, REBIN)]))
        channel = SPEC.channels[0]
        for values in POINTS:
            result = model(**values).single()
            native_grid = backend.to_numpy(model.native_grid(channel))  # type: ignore[attr-defined]
            native_flux = backend.to_numpy(model.native_flux(channel, values))  # type: ignore[attr-defined]
            np.testing.assert_allclose(
                native_grid,
                np.asarray(result.spectral_axis.values, dtype=float),
                rtol=0.0,
                atol=0.0,
            )
            np.testing.assert_allclose(
                native_flux, backend.to_numpy(result.values), rtol=0.0, atol=tolerances.exact
            )

    def test_an_unknown_channel_is_refused_by_name(self, backend: ConformanceBackend) -> None:
        model = portable(backend)
        with pytest.raises(Exception, match="elsewhere"):
            model.native_flux("elsewhere", TRUTH)  # type: ignore[attr-defined]


# ---------------------------------------------------------------------------
# (d) NUTS on each differentiable fixture
# ---------------------------------------------------------------------------


class TestNUTS:
    def test_it_recovers_the_injected_truth(self, backend: ConformanceBackend) -> None:
        problem = rebinned_problem(backend)
        if problem.backend not in registered_realisations():
            pytest.skip(
                f"the {problem.backend!r} backend registers no realisation "
                f"(inference.md §10a: the reference backend has no differentiable path)"
            )
        from ampere.inference import NUTSEngine

        assert problem.differentiable is True
        assert problem.backend == backend.name
        # The realised density is the numpy path's, at the fixture's tolerance.
        realisation = realise(problem)
        y = problem.unconstrain(qualified(TRUTH))
        assert float(backend.to_numpy(realisation.log_prob_unconstrained(y))) == pytest.approx(  # type: ignore[attr-defined]
            problem.log_prob_unconstrained(y), abs=backend.capabilities.tolerances.cross_backend
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            run = NUTSEngine(problem).run(draws=200, warmup=200, chains=2)
        for name, truth in TRUTH.items():
            draws = np.asarray(run["posterior"][f"model.{name}"])
            assert abs(float(draws.mean()) - truth) < 5.0 * float(draws.std()), name
        assert run.attrs["ampere_backend"] == backend.name
        assert run.attrs["ampere_engine"] == "nuts"


# ---------------------------------------------------------------------------
# (e) the capability flags
# ---------------------------------------------------------------------------


class TestTheFlags:
    def test_the_twin_declares_its_fixture(self, backend: ConformanceBackend) -> None:
        model = portable(backend)
        assert model.BACKEND == backend.name
        assert model.DIFFERENTIABLE is backend.capabilities.differentiable
        assert model.BATCHABLE is backend.capabilities.batchable
        assert model.DEVICE == "cpu"
