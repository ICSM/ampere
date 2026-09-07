"""The VI driver: a fitted guide, and the honesty about what it is.

``test_nuts.py`` holds the gradient-based *sampler* to a posterior written down
in closed form. This file holds the gradient-based *optimiser* to the same
posterior, and to three things a sampler is not asked about:

1. **It finds the right location.** The agreement problem is deliberately
   conjugate — a power law with its index held fixed is linear in ``norm``, so
   a Gaussian prior and Gaussian noise give a Gaussian posterior whose mean and
   variance are arithmetic — and a Gaussian guide can represent that posterior
   *exactly*. So on this problem VI is not an approximation at all, and the row
   can be tight: an optimiser that converged must land on the closed form.
2. **It is honest about the approximation everywhere else.** The mean-field
   guide understates the variance of a correlated posterior, and the run
   records the guide family, the step count and the ELBO trace so that a reader
   of the archived file can see which family was fitted and whether the
   optimisation converged. Those attrs are asserted, because an approximation
   that does not say so in provenance is the failure mode this driver is most
   likely to cause.
3. **It refuses, by name, what it cannot fit** — a backend with no realisation,
   a problem that declares itself non-differentiable, an unknown guide family,
   and a density that disagrees with the problem it was handed.

Budgets are small and seeds fixed, as in ``test_nuts.py``: this belongs in the
per-PR gate. VI is cheap enough that the budgets here cost less than one NUTS
row.
"""

from __future__ import annotations

import importlib
import json
import math
import warnings
from typing import Any

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.core import (
    Dataset,
    FittingProblem,
    GaussianFamily,
    Likelihood,
    Spectrum,
)
from ampere.inference import EngineError, VIEngine
from ampere.inference._vi import GUIDES, VI_LIBRARIES, supported_backends

REFERENCE_WAVELENGTH = 1.0
SEED = 20260907

#: jax only, today. ``VI_LIBRARIES`` is the table and this module follows it
#: rather than hard-coding a name, so the day a pyro route lands this file runs
#: on both without an edit.
BACKENDS = sorted(VI_LIBRARIES)


def _installed() -> list[str]:
    found: list[str] = []
    for name in BACKENDS:
        try:
            module = importlib.import_module(f"ampere.backends.{name}")
        except ImportError:  # the extra is not installed in this environment
            continue
        if name == "jax":
            # ``lowering.md`` §10.2(a): the *application* turns the flag on.
            module.configure_x64()
        found.append(name)
    return found


INSTALLED = _installed()

pytestmark = pytest.mark.skipif(
    not INSTALLED,
    reason="no backend with a variational route installed; VI needs a registered realisation",
)


@pytest.fixture(scope="module", params=INSTALLED)
def backend(request: Any) -> Any:
    return importlib.import_module(f"ampere.backends.{request.param}")


# ---------------------------------------------------------------------------
# The conjugate problem, whose posterior is arithmetic
# ---------------------------------------------------------------------------

GRID = np.geomspace(1.0, 10.0, 20)
SIGMA = 0.2
INDEX = -1.0
PRIOR = (2.0, 0.5)  # (mean, sd) of the Gaussian prior on `norm`


def _noisy(grid: np.ndarray, values: np.ndarray, sigma: float, seed: int) -> Spectrum:
    rng = np.random.default_rng(seed)
    return Spectrum(
        grid * u.micron,
        (values + rng.normal(0.0, sigma, values.size)) * u.Jy,
        uncertainty=np.full(values.size, sigma) * u.Jy,
    )


DATA = _noisy(GRID, 2.0 * (GRID / REFERENCE_WAVELENGTH) ** INDEX, SIGMA, seed=7)


def analytic() -> tuple[float, float]:
    """The posterior on ``norm``, in closed form: ``(mean, sd)``.

    The same conjugate construction ``test_nuts.py`` uses, and written out for
    the same reason: an oracle produced by a sampler is not an oracle. It
    matters more here, because a Gaussian guide can represent this posterior
    exactly — so the closed form is not merely the right answer, it is an
    answer this driver has no excuse to miss.
    """
    x = (GRID / REFERENCE_WAVELENGTH) ** INDEX
    y = np.asarray(DATA.values)
    mu, tau = PRIOR
    precision = 1.0 / tau**2 + float(np.sum(x**2)) / SIGMA**2
    mean = (mu / tau**2 + float(np.sum(x * y)) / SIGMA**2) / precision
    return mean, 1.0 / math.sqrt(precision)


def agreement_problem(backend: Any, seed: int | None = SEED) -> FittingProblem:
    """One dataset, one free parameter, an exactly known Gaussian posterior."""
    return FittingProblem(
        backend.PowerLaw(
            GRID,
            norm=st.norm(*PRIOR),
            index=INDEX,
            reference_wavelength=REFERENCE_WAVELENGTH,
        ),
        [Dataset(DATA, likelihood=Likelihood(GaussianFamily(), backend.IndependentNoise()))],
        seed=seed,
    )


def a_different_problem(backend: Any, seed: int | None = SEED) -> FittingProblem:
    """The *same shape* as :func:`agreement_problem`, and a different posterior.

    Same free dimension deliberately: a density of the wrong *length* would be
    caught by the first thing that evaluated it, and the check this exercises
    is the one for a density that runs perfectly well and describes something
    else — which is the failure that would otherwise be silent.
    """
    return FittingProblem(
        backend.PowerLaw(
            GRID,
            norm=st.norm(*PRIOR),
            index=-0.4,
            reference_wavelength=REFERENCE_WAVELENGTH,
        ),
        [Dataset(DATA, likelihood=Likelihood(GaussianFamily(), backend.IndependentNoise()))],
        seed=seed,
    )


def fit(problem: FittingProblem, **settings: Any) -> Any:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return VIEngine(problem).run(**settings)


@pytest.fixture(scope="module")
def agreement_run(backend: Any) -> Any:
    return fit(agreement_problem(backend), draws=2000, steps=4000)


# ---------------------------------------------------------------------------
# 1. It finds the posterior it can represent exactly
# ---------------------------------------------------------------------------


class TestTheConjugateFit:
    def test_it_recovers_the_analytic_mean(self, agreement_run: Any) -> None:
        mean, sd = analytic()
        drawn = float(agreement_run["posterior"]["model.norm"].mean())
        assert abs(drawn - mean) < 0.3 * sd

    def test_it_recovers_the_analytic_width(self, agreement_run: Any) -> None:
        """A Gaussian guide can represent a Gaussian posterior *exactly*, so the
        usual mean-field variance deficit has nowhere to come from here — and
        that is what makes this a real test of the fit rather than of the
        approximation."""
        _, sd = analytic()
        drawn = float(agreement_run["posterior"]["model.norm"].std())
        assert drawn == pytest.approx(sd, rel=0.15)

    def test_the_elbo_rose(self, agreement_run: Any) -> None:
        """The optimisation is meant to have optimised something."""
        attrs = agreement_run.attrs
        assert attrs["ampere_vi_elbo_final"] > attrs["ampere_vi_elbo_initial"]

    def test_a_full_covariance_guide_fits_the_same_posterior(self, backend: Any) -> None:
        mean, sd = analytic()
        run = fit(
            agreement_problem(backend),
            draws=2000,
            steps=4000,
            guide="multivariate_normal",
        )
        assert abs(float(run["posterior"]["model.norm"].mean()) - mean) < 0.3 * sd
        assert run.attrs["ampere_vi_guide"] == "multivariate_normal"


# ---------------------------------------------------------------------------
# 2. The run says what it is
# ---------------------------------------------------------------------------


class TestTheRunRecordsWhatItIs:
    def test_the_engine_and_backend_are_recorded(self, agreement_run: Any, backend: Any) -> None:
        assert agreement_run.attrs["ampere_engine"] == "vi"
        assert agreement_run.attrs["ampere_backend"] == backend.BACKEND

    def test_the_guide_family_and_step_count_are_recorded(self, agreement_run: Any) -> None:
        """An approximation that does not say which family it approximated
        within is not a recorded approximation."""
        assert agreement_run.attrs["ampere_vi_guide"] in GUIDES
        assert agreement_run.attrs["ampere_vi_steps"] == 4000
        assert agreement_run.attrs["ampere_vi_optimiser"] == "adam"

    def test_the_elbo_trace_is_recorded_with_its_stride(self, agreement_run: Any) -> None:
        """The trace is what says whether the optimisation converged, and a
        thinned trace with no stride would be a mis-labelled x axis."""
        trace = json.loads(agreement_run.attrs["ampere_vi_elbo_trace"])
        stride = agreement_run.attrs["ampere_vi_elbo_trace_stride"]
        assert isinstance(trace, list) and trace
        assert len(trace) * stride >= agreement_run.attrs["ampere_vi_steps"] - stride
        assert trace[-1] == pytest.approx(agreement_run.attrs["ampere_vi_elbo_final"])

    def test_the_draws_are_one_chain(self, agreement_run: Any) -> None:
        """A guide is a distribution, not a chain: the draws are independent and
        there is nothing for a second chain to mean."""
        assert agreement_run["posterior"]["model.norm"].shape == (1, 2000)

    def test_it_records_that_the_draws_came_through_a_realisation(self, agreement_run: Any) -> None:
        assert agreement_run.attrs["ampere_realised"] == 1

    def test_no_stored_draw_was_re_evaluated(self, agreement_run: Any) -> None:
        """``inference.md`` §10a's optional decomposition, consumed: the driver
        already holds every draw's per-dataset log-likelihood in the arithmetic
        that produced it, so ``Engine.finish`` is handed it rather than
        recomputing 2000 model evaluations on the numpy path."""
        assert agreement_run.attrs["ampere_engine_draws_recomputed"] == 0

    def test_the_per_dataset_decomposition_is_there(self, agreement_run: Any) -> None:
        assert list(agreement_run["log_likelihood"].dataset.data_vars) == ["default"]


# ---------------------------------------------------------------------------
# 3. It refuses what it cannot fit
# ---------------------------------------------------------------------------


class TestRefusals:
    def test_an_unknown_guide_family_is_refused_by_name(self, backend: Any) -> None:
        with pytest.raises(EngineError, match="guide family"):
            VIEngine(agreement_problem(backend)).run(draws=10, guide="planar_flow")

    def test_a_backend_with_no_variational_route_is_refused(self) -> None:
        """The reference backend registers no realisation, deliberately."""
        from ampere.backends.reference import PowerLaw as ReferencePowerLaw

        # ``Dataset``'s default noise model is ``ampere.core``'s, which declares
        # the reference backend -- so this problem is entirely on rung 1 and
        # registers no realisation, by W2.13's fifth sub-decision.
        problem = FittingProblem(
            ReferencePowerLaw(GRID, norm=st.norm(*PRIOR), index=INDEX),
            [Dataset(DATA)],
            seed=SEED,
        )
        assert problem.backend not in supported_backends()
        with pytest.raises(EngineError, match="cannot fit a problem on the 'reference' backend"):
            VIEngine(problem)

    def test_a_nonsensical_step_count_is_refused(self, backend: Any) -> None:
        with pytest.raises(EngineError, match="at least one optimisation step"):
            VIEngine(agreement_problem(backend)).run(draws=10, steps=0)

    def test_a_nonsensical_learning_rate_is_refused(self, backend: Any) -> None:
        with pytest.raises(EngineError, match="learning rate"):
            VIEngine(agreement_problem(backend)).run(draws=10, learning_rate=-1.0)

    def test_a_density_that_disagrees_with_the_problem_is_refused(self, backend: Any) -> None:
        """The failure this check exists for is silent otherwise: a guide fitted
        to the wrong problem looks perfectly healthy."""
        problem = agreement_problem(backend)
        other = a_different_problem(backend)
        with pytest.raises(EngineError, match="disagrees with the problem"):
            VIEngine(problem, backend.lower_problem(other).log_prob_unconstrained)


# ---------------------------------------------------------------------------
# 4. Reproducibility
# ---------------------------------------------------------------------------


class TestReproducibility:
    def test_the_same_seed_gives_the_same_fit(self, backend: Any) -> None:
        """Every stream is derived from the problem's own seed
        (``inference.md`` §12): the optimiser's, and the draws from the guide."""
        first = fit(agreement_problem(backend), draws=200, steps=500)
        second = fit(agreement_problem(backend), draws=200, steps=500)
        assert np.asarray(first["posterior"]["model.norm"]) == pytest.approx(
            np.asarray(second["posterior"]["model.norm"])
        )

    def test_a_different_seed_gives_a_different_fit(self, backend: Any) -> None:
        first = fit(agreement_problem(backend, seed=1), draws=200, steps=500)
        second = fit(agreement_problem(backend, seed=2), draws=200, steps=500)
        assert not np.allclose(
            np.asarray(first["posterior"]["model.norm"]),
            np.asarray(second["posterior"]["model.norm"]),
        )
