"""W4.6's engine rows: a wrapped astropy model fitted, and the capability ladder honoured.

The adapter's whole claim about inference is a claim about *capabilities*, and
the only way to check a capability claim is to run the engines it names and to
watch the ones it disowns refuse. So:

* :class:`TestAnEmceeFitOfAWrappedModel` — an emcee fit of a wrapped
  ``astropy.modeling`` model recovers the injected parameters. Short and
  seeded, like every run in ``tests/inference``.
* :class:`TestTheSBIEngineRunsOnAWrappedModel` — the black-box route:
  :class:`~ampere.inference.SBIEngine` simulating from the same problem through
  ``simulate_many`` under a **process pool**, which is the executor a black box
  wants and is also the row that proves the wrapper (astropy model, its ties,
  its template and its hoisted unit factor) survives pickling to a worker.
  Skipped without the ``sbi`` extra (``pixi run -e sbi``).
* :class:`TestTheGradientEnginesRefuse` — ``NUTSEngine`` and ``VIEngine`` never
  reach a wrapped model, and the refusal is by name. That is the honest answer
  §4.7 promises, asserted rather than described.

They live in ``tests/core/`` beside ``test_astropy_compat.py`` because W4.6
owns ``tests/core/test_astropy*`` and because the subject is the adapter rather
than the drivers; ``ampere.inference`` is a base-install import, so the
directory pays nothing for it.
"""

from __future__ import annotations

import importlib.util

import astropy.units as u
import numpy as np
import pytest
from astropy.modeling.models import PowerLaw1D

from ampere.core import (
    Dataset,
    FittingProblem,
    ProcessExecutor,
    Spectrum,
    from_astropy,
    realise,
)
from ampere.core.exceptions import CapabilityError, LoweringError
from ampere.inference import EmceeEngine, SBIEngine

SEED = 20260911
WAVELENGTH = np.geomspace(1.0, 20.0, 16)
TRUTH = {"model.amplitude": 2.0, "model.alpha": 1.2}
SIGMA = 0.05

HAS_SBI = importlib.util.find_spec("sbi") is not None
needs_sbi = pytest.mark.skipif(not HAS_SBI, reason="needs the 'sbi' extra (pixi run -e sbi ...)")


def wrapped_power_law() -> object:
    """``PowerLaw1D`` in Jy, with ``x_0`` fixed and both bounds finite.

    Bounded priors throughout, so the problem is one ``SBIEngine``'s prior
    bridge is happy with and one emcee can be started in.
    """
    model = PowerLaw1D(amplitude=2.0 * u.Jy, x_0=1.0 * u.micron, alpha=1.2)
    model.x_0.fixed = True
    model.amplitude.bounds = (0.5, 5.0)
    model.alpha.bounds = (0.2, 2.5)
    return from_astropy(model, grid=WAVELENGTH * u.micron, output_unit=u.Jy)


def problem(seed: int | None = SEED) -> FittingProblem:
    """One dataset generated from the wrapped model itself, plus Gaussian noise."""
    model = wrapped_power_law()
    truth = model(amplitude=TRUTH["model.amplitude"], alpha=TRUTH["model.alpha"]).single()
    rng = np.random.default_rng(SEED)
    observed = Spectrum(
        WAVELENGTH * u.micron,
        (truth.values + rng.normal(0.0, SIGMA, truth.values.size)) * u.Jy,
        uncertainty=np.full(truth.values.size, SIGMA) * u.Jy,
    )
    return FittingProblem(wrapped_power_law(), [Dataset(observed)], seed=seed)


@pytest.fixture(scope="module")
def emcee_run() -> object:
    """One short, seeded emcee fit. Module-scoped: it happens once."""
    return EmceeEngine(problem(), walkers=12).run(steps=250, burn_in=80)


class TestAnEmceeFitOfAWrappedModel:
    """The gradient-free route the adapter promises, exercised end to end."""

    def test_it_runs_and_emits_the_run(self, emcee_run: object) -> None:
        posterior = emcee_run["posterior"].dataset  # type: ignore[index]
        assert set(posterior.data_vars) == {"model.amplitude", "model.alpha"}
        assert (posterior.sizes["chain"], posterior.sizes["draw"]) == (12, 170)
        assert emcee_run.attrs["ampere_backend"] == "reference"  # type: ignore[attr-defined]

    def test_it_recovers_the_injected_parameters(self, emcee_run: object) -> None:
        """Inside the central 95 %, which is what a correct fit of clean data does."""
        posterior = emcee_run["posterior"].dataset  # type: ignore[index]
        for name, truth in TRUTH.items():
            draws = np.asarray(posterior[name]).ravel()
            low, high = np.quantile(draws, (0.025, 0.975))
            assert low <= truth <= high, f"{name}: {truth} outside [{low}, {high}]"

    def test_the_astropy_model_is_still_the_one_being_evaluated(self) -> None:
        """No substitution anywhere on this path: the wrapper holds the real model."""
        composed = problem()
        assert type(composed.model.astropy_model).__name__ == "PowerLaw1D"  # type: ignore[attr-defined]


@needs_sbi
class TestTheSBIEngineRunsOnAWrappedModel:
    """The black-box route, under the process pool a black box wants."""

    @pytest.fixture(scope="class")
    def sbi_run(self) -> object:
        engine = SBIEngine(
            problem(),
            method="npe",
            budget=120,
            executor=ProcessExecutor(2),
            chunk_size=60,
        )
        return engine.run(draws=40, training={"max_num_epochs": 15})

    def test_it_fits(self, sbi_run: object) -> None:
        posterior = sbi_run["posterior"].dataset  # type: ignore[index]
        assert set(posterior.data_vars) == {"model.amplitude", "model.alpha"}
        assert (posterior.sizes["chain"], posterior.sizes["draw"]) == (1, 40)
        assert np.all(np.isfinite(np.asarray(posterior["model.amplitude"])))

    def test_the_budget_went_through_the_process_pool(self, sbi_run: object) -> None:
        """Which is also the proof the wrapped astropy model pickles to a worker."""
        attrs = sbi_run.attrs  # type: ignore[attr-defined]
        assert attrs["ampere_sbi_executor"] == "ProcessExecutor"
        assert attrs["ampere_sbi_simulations"] == 120
        assert attrs["ampere_backend"] == "reference"


class TestTheGradientEnginesRefuse:
    """``NUTSEngine`` and ``VIEngine`` never reach a wrapped model, and say why."""

    def test_realising_a_wrapped_problem_refuses_by_name(self) -> None:
        with pytest.raises((LoweringError, CapabilityError)) as raised:
            realise(problem())
        assert "reference" in str(raised.value)

    def test_the_problem_declares_itself_undifferentiable(self) -> None:
        composed = problem()
        assert composed.capabilities.differentiable is False
        assert composed.capabilities.backend == "reference"
