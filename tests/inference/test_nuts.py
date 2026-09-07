"""The NUTS driver, end to end: gradients, refusals, and the same emitted run.

``tests/inference/test_engines.py`` holds the three gradient-free drivers to
five things; this file holds the fourth driver to the same ones, minus what is
structurally shared and already covered there (the import-graph check, the
netCDF round trip, the failure path). What is *new* is what NUTS is for:

1. **It recovers the toy joint problem's posterior.** ``inference.md`` §15's
   shape — one power law on two channels, two genuinely different instrument
   chains, the calibration factors tied — built entirely from
   ``ampere.backends.jax``, sampled through a gradient, and held to the truth
   the data were generated at.
2. **It agrees with a posterior written down in closed form.** The agreement
   problem is deliberately conjugate: a power law with its index held fixed is
   *linear* in ``norm``, so a Gaussian prior and Gaussian noise give a Gaussian
   posterior whose mean and variance are arithmetic rather than a long
   reference run. Comparing a sampler only against other samplers would pass
   several samplers that are wrong in the same way.
3. **It refuses, by name, what it cannot sample.** A problem on another
   backend, a problem that declares itself non-differentiable, and — the one
   that would otherwise be silent — a density that disagrees with the problem
   it was handed.
4. **Every run emits the run**, through the same
   :meth:`~ampere.inference.engine.Engine.finish`, with ``ampere_backend`` read
   off the problem rather than declared by the driver.

Budgets are small and seeds are fixed, as in ``test_engines.py``: this module
belongs in the per-PR gate beside the other new-namespace suites rather than on
a nightly. The tolerance is stated in units of the analytic posterior's own
standard deviation and is several times the Monte Carlo standard error at these
budgets — loose enough that a correct sampler passes essentially always, and
far tighter than the errors a driver bug produces, since a mis-transposed draw
array or an unconstrained-space posterior recorded as a constrained one moves a
summary by whole standard deviations rather than fractions of one.
"""

from __future__ import annotations

import math
import warnings
from typing import Any

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

pytest.importorskip("jax")
pytest.importorskip("numpyro")

from ampere.backends.jax import (
    CalibrationScale,
    PowerLaw,
    Resample,
    configure_x64,
    lower_problem,
)
from ampere.backends.reference import PowerLaw as ReferencePowerLaw
from ampere.core import (
    Dataset,
    DatasetCollection,
    FittingProblem,
    Instrument,
    Spectrum,
    Tie,
)
from ampere.inference import EngineError, NUTSEngine

REFERENCE_WAVELENGTH = 1.0
TRUTH = {"norm": 2.0, "index": -1.2, "calibration": 1.0}
SEED = 20260907

FINE = np.geomspace(1.0, 20.0, 40)
COARSE = np.geomspace(1.5, 15.0, 12)


@pytest.fixture(autouse=True, scope="session")
def _x64() -> None:
    """``lowering.md`` §10.2(a): the application turns the flag on, not ampere."""
    configure_x64()


def power_law(grid: np.ndarray, norm: float, index: float) -> np.ndarray:
    return norm * (grid / REFERENCE_WAVELENGTH) ** index


def noisy(grid: np.ndarray, values: np.ndarray, sigma: float, seed: int) -> Spectrum:
    """One observed spectrum: the truth plus Gaussian noise from a fixed stream."""
    rng = np.random.default_rng(seed)
    return Spectrum(
        grid * u.micron,
        (values + rng.normal(0.0, sigma, values.size)) * u.Jy,
        uncertainty=np.full(values.size, sigma) * u.Jy,
    )


BLUE_DATA = noisy(FINE, power_law(FINE, TRUTH["norm"], TRUTH["index"]), 0.08, seed=11)
RED_DATA = noisy(COARSE, power_law(COARSE, TRUTH["norm"], TRUTH["index"]), 0.05, seed=13)


def joint_problem(seed: int | None = SEED) -> FittingProblem:
    """``inference.md`` §15's shape, on the shipped **jax** backend.

    Deliberately the same declaration ``test_engines.py``'s ``joint_problem``
    builds on the reference backend, so that "NUTS recovers the toy joint
    problem's posterior" means the same problem the other three drivers are
    held to. The tie costs a dimension: three free parameters, not four.
    """
    model = PowerLaw(
        FINE,
        norm=st.lognorm(0.4, scale=2.0),
        index=st.norm(-1.2, 0.3),
        reference_wavelength=REFERENCE_WAVELENGTH,
        channels=("blue", "red"),
    )
    return FittingProblem(
        model,
        DatasetCollection(
            {
                "blue": Dataset(
                    BLUE_DATA,
                    Instrument(
                        [CalibrationScale(st.lognorm(0.05), label="calibration")],
                        channel="blue",
                    ),
                ),
                "red": Dataset(
                    RED_DATA,
                    Instrument(
                        [
                            Resample(COARSE),
                            CalibrationScale(st.lognorm(0.05), label="calibration"),
                        ],
                        channel="red",
                    ),
                ),
            }
        ),
        ties=[
            Tie(
                "calibration",
                ("blue.instrument.calibration.scale", "red.instrument.calibration.scale"),
            )
        ],
        seed=seed,
    )


# -- the conjugate problem the driver is held to ----------------------------

AGREEMENT_GRID = np.geomspace(1.0, 10.0, 20)
AGREEMENT_SIGMA = 0.2
AGREEMENT_INDEX = -1.0
AGREEMENT_PRIOR = (2.0, 0.5)  # (mean, sd) of the Gaussian prior on `norm`
AGREEMENT_DATA = noisy(
    AGREEMENT_GRID, power_law(AGREEMENT_GRID, 2.0, AGREEMENT_INDEX), AGREEMENT_SIGMA, seed=7
)


def analytic() -> tuple[float, float]:
    """The posterior on ``norm``, in closed form: ``(mean, sd)``.

    With ``index`` fixed the power law is *linear* in ``norm`` — the prediction
    is ``norm * x`` for a known design vector ``x`` — so a Gaussian prior and
    known Gaussian noise conjugate exactly. Written out here rather than
    estimated from one long reference run, because an oracle a sampler produced
    is not an oracle.
    """
    x = (AGREEMENT_GRID / REFERENCE_WAVELENGTH) ** AGREEMENT_INDEX
    y = np.asarray(AGREEMENT_DATA.values)
    mu, tau = AGREEMENT_PRIOR
    precision = 1.0 / tau**2 + float(np.sum(x**2)) / AGREEMENT_SIGMA**2
    mean = (mu / tau**2 + float(np.sum(x * y)) / AGREEMENT_SIGMA**2) / precision
    return mean, 1.0 / math.sqrt(precision)


def agreement_problem(seed: int | None = SEED) -> FittingProblem:
    """One dataset, one free parameter, an exactly known posterior."""
    return FittingProblem(
        PowerLaw(
            AGREEMENT_GRID,
            norm=st.norm(*AGREEMENT_PRIOR),
            index=AGREEMENT_INDEX,
            reference_wavelength=REFERENCE_WAVELENGTH,
        ),
        [Dataset(AGREEMENT_DATA)],
        seed=seed,
    )


def sample(problem: FittingProblem, **settings: Any) -> Any:
    """Build the driver from the backend's lowering and run it."""
    engine = NUTSEngine(problem, lower_problem(problem).log_prob_unconstrained)
    with warnings.catch_warnings():
        # A short chain on a three-dimensional posterior occasionally records a
        # failure at a start point; the failure path itself is test_engines.py's.
        warnings.simplefilter("ignore")
        return engine.run(**settings)


@pytest.fixture(scope="module")
def joint_run() -> Any:
    return sample(joint_problem(), draws=400, warmup=400, chains=2)


@pytest.fixture(scope="module")
def agreement_run() -> Any:
    return sample(agreement_problem(), draws=600, warmup=400, chains=2)


# ---------------------------------------------------------------------------
# 1. The toy joint problem samples end to end
# ---------------------------------------------------------------------------


class TestTheJointProblemSamples:
    def test_the_problem_is_the_one_the_contract_describes(self) -> None:
        """Three free dimensions, not four: the tie costs one (``inference.md`` §15)."""
        problem = joint_problem()
        assert problem.free_size == 3
        assert problem.shared_names == ("calibration",)
        assert problem.backend == "jax"
        assert problem.differentiable is True

    def test_the_run_has_the_shape_it_was_asked_for(self, joint_run: Any) -> None:
        assert joint_run["posterior"]["model.norm"].shape == (2, 400)

    def test_every_free_parameter_is_in_the_posterior_under_its_merged_name(
        self, joint_run: Any
    ) -> None:
        assert set(joint_run["posterior"].dataset.data_vars) == {
            "model.norm",
            "model.index",
            "calibration",
        }

    def test_it_recovers_the_truth(self, joint_run: Any) -> None:
        """The gradient-based half of ``test_engines.py``'s criterion 1."""
        for name, key in (
            ("model.norm", "norm"),
            ("model.index", "index"),
            ("calibration", "calibration"),
        ):
            draws = np.asarray(joint_run["posterior"][name])
            spread = float(draws.std())
            assert abs(float(draws.mean()) - TRUTH[key]) < 4.0 * spread

    def test_the_run_carries_the_decomposition_every_run_carries(self, joint_run: Any) -> None:
        stats = joint_run["sample_stats"].dataset
        assert {"lp", "log_prior", "log_likelihood"} <= set(stats.data_vars)
        total = np.asarray(stats["log_prior"]) + np.asarray(stats["log_likelihood"])
        assert np.asarray(stats["lp"]) == pytest.approx(total, abs=1e-9)

    def test_the_per_dataset_terms_are_there(self, joint_run: Any) -> None:
        assert set(joint_run["log_likelihood"].dataset.data_vars) == {"blue", "red"}

    def test_the_backend_is_derived_not_declared(self, joint_run: Any) -> None:
        """W2.12: no driver takes a ``backend=``; the flag is read off the problem."""
        assert joint_run.attrs["ampere_engine"] == "nuts"
        assert joint_run.attrs["ampere_backend"] == "jax"

    def test_the_nuts_diagnostics_are_recorded(self, joint_run: Any) -> None:
        """Divergences are the diagnostic no gradient-free engine can offer."""
        assert "ampere_nuts_divergences" in joint_run.attrs
        assert joint_run.attrs["ampere_nuts_chains"] == 2
        assert joint_run.attrs["ampere_nuts_warmup"] == 400
        assert 0.0 < float(joint_run.attrs["ampere_nuts_step_size"])
        assert 0.0 < float(joint_run.attrs["ampere_nuts_mean_accept_prob"]) <= 1.0

    def test_the_posterior_is_in_the_constrained_space(self, joint_run: Any) -> None:
        """NUTS works in unconstrained space; a run that stored *that* would look
        entirely healthy and be a posterior over the wrong quantity."""
        norm = np.asarray(joint_run["posterior"]["model.norm"])
        assert bool(np.all(norm > 0.0))  # lognormal support


# ---------------------------------------------------------------------------
# 2. Agreement with a posterior that is arithmetic
# ---------------------------------------------------------------------------


class TestAgreementWithTheClosedForm:
    def test_the_posterior_mean_and_width_match_the_conjugate_answer(
        self, agreement_run: Any
    ) -> None:
        mean, sd = analytic()
        draws = np.asarray(agreement_run["posterior"]["model.norm"]).ravel()
        assert abs(float(draws.mean()) - mean) < 0.25 * sd
        assert abs(float(draws.std()) - sd) < 0.25 * sd

    def test_the_fixed_index_takes_no_sampler_dimension(self) -> None:
        problem = agreement_problem()
        assert problem.free_size == 1
        assert problem.parameters.fixed_names == ("model.index",)


# ---------------------------------------------------------------------------
# 3. Refusals
# ---------------------------------------------------------------------------


class TestRefusals:
    def test_a_problem_on_another_backend_is_refused_by_name(self) -> None:
        problem = FittingProblem(
            ReferencePowerLaw(
                AGREEMENT_GRID,
                norm=st.norm(*AGREEMENT_PRIOR),
                index=AGREEMENT_INDEX,
                reference_wavelength=REFERENCE_WAVELENGTH,
            ),
            [Dataset(AGREEMENT_DATA)],
        )
        with pytest.raises(EngineError, match="'reference' backend"):
            NUTSEngine(problem, lambda y: 0.0)

    def test_a_density_that_is_not_callable_is_refused(self) -> None:
        with pytest.raises(EngineError, match="callable"):
            NUTSEngine(agreement_problem(), 3.0)  # type: ignore[arg-type]

    def test_a_density_for_a_different_problem_is_refused(self) -> None:
        """The refusal that would otherwise be silent.

        Sampling one problem's density while recording another's produces a run
        that passes every structural check and describes a posterior nobody
        asked for. One point is enough to catch it and costs nothing.
        """
        problem = agreement_problem()
        other = lower_problem(joint_problem())
        with pytest.raises(EngineError, match="disagrees with the problem"):
            NUTSEngine(problem, lambda y: other.log_prob_unconstrained(np.zeros(3)))

    def test_a_problem_with_nothing_to_sample_is_refused(self) -> None:
        problem = FittingProblem(
            PowerLaw(
                AGREEMENT_GRID, norm=2.0, index=-1.0, reference_wavelength=REFERENCE_WAVELENGTH
            ),
            [Dataset(AGREEMENT_DATA)],
        )
        with pytest.raises(EngineError, match="nothing to sample"):
            NUTSEngine(problem, lambda y: 0.0)


# ---------------------------------------------------------------------------
# 4. Reproducibility
# ---------------------------------------------------------------------------


class TestReproducibility:
    def test_the_same_seed_gives_the_same_draws(self) -> None:
        """``inference.md`` §12: every stream comes from the problem's own seed.

        Initialisation and the sampler's own randomness are separate labels
        under this engine's name, so neither can make the other
        irreproducible.
        """
        first = sample(joint_problem(SEED), draws=30, warmup=30, chains=1)
        second = sample(joint_problem(SEED), draws=30, warmup=30, chains=1)
        assert np.asarray(first["posterior"]["model.norm"]) == pytest.approx(
            np.asarray(second["posterior"]["model.norm"]), abs=0.0
        )

    def test_a_different_seed_gives_different_draws(self) -> None:
        first = sample(joint_problem(SEED), draws=30, warmup=30, chains=1)
        other = sample(joint_problem(SEED + 1), draws=30, warmup=30, chains=1)
        assert not np.allclose(
            np.asarray(first["posterior"]["model.norm"]),
            np.asarray(other["posterior"]["model.norm"]),
        )
