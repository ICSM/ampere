"""W7.12: the ensemble engines' default start, its fallbacks, the verdict, the size warning.

The default start of :class:`~ampere.inference.EmceeEngine` and
:class:`~ampere.inference.ZeusEngine` is a ball at the one-start scipy
optimum; ``initial="prior"`` is the former default, bit for bit; a problem
the optimiser cannot start falls back to the prior with a
:class:`~ampere.inference.DefaultStartWarning`; one whose optimum sits on a
prior bound starts mixed (W7.17: prior draws on the pinned coordinates, the ball
on the rest), silently.
:func:`~ampere.results.check_convergence` is the verdict ``emit`` and
``summary`` warn with. The regression row is the user-journeys memo's
Appendix A problem (``docs/design/user_journeys_memo.md``): the first fit a
user writes on the beta, which reached R-hat 1.5 from the prior.
"""

from __future__ import annotations

import json
import math
import warnings
from typing import Any

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.backends.reference import ModifiedBlackBody, PowerLaw, SyntheticPhotometry
from ampere.core import (
    Dataset,
    FittingProblem,
    GaussianFamily,
    IndependentNoise,
    Instrument,
    Likelihood,
    Model,
    Parameter,
    PhotometricPoints,
    Spectrum,
    negotiate,
)
from ampere.inference import (
    ENSEMBLE_SIZE_WARNING,
    DefaultStartWarning,
    EmceeEngine,
    EngineError,
    EnsembleSizeWarning,
    ZeusEngine,
)
from ampere.results import ConvergenceVerdict, ResultsWarning, check_convergence, summary

SEED = 20261007
#: Twice the measured spread (0.030) of the offset column, which starts in the ball.
OFFSET_BALL_BOUND = 0.06

# ---------------------------------------------------------------------------
# Problems
# ---------------------------------------------------------------------------

GRID = np.geomspace(1.0, 10.0, 20)
SIGMA = 0.2


def conjugate_problem(seed: int | None = SEED) -> FittingProblem:
    """A one-parameter linear model with a normal prior: the posterior is a Gaussian."""
    rng = np.random.default_rng(11)
    values = 2.0 * (GRID / 3.0) ** -1.0
    observed = Spectrum(
        GRID * u.micron,
        (values + rng.normal(0.0, SIGMA, GRID.size)) * u.Jy,
        uncertainty=np.full(GRID.size, SIGMA) * u.Jy,
    )
    model = PowerLaw(GRID, norm=st.norm(2.0, 0.5), index=-1.0, reference_wavelength=3.0)
    return FittingProblem(model, [Dataset(observed)], seed=seed)


EDGE_GRID = np.array([1.0, 2.0, 3.0, 4.0])


class AtTheEdge(Model):
    """``slope * x + offset`` with ``slope`` in ``[0, 1]`` and data at slope 2; ``width`` unused.

    ``tests/inference/test_optimise.py``'s bound case (W7.6): the mode sits at
    ``slope = 1``, a box's upper bound, and at ``width = 0``, a half-line's floor.
    """

    def __init__(self) -> None:
        self.register_buffer("wavelength", EDGE_GRID, unit=u.um)
        self.register_parameter(Parameter("slope", st.uniform(0.0, 1.0)))
        self.register_parameter(Parameter("width", st.halfnorm(scale=1.0)))
        self.register_parameter(Parameter("offset", st.norm(0.0, 1.0)))

    def evaluate(self, **values: Any) -> Spectrum:
        ctx = self.context(values)
        flux = ctx["slope"] * ctx["wavelength"] + 0.0 * ctx["width"] + ctx["offset"]
        return Spectrum(ctx["wavelength"] * u.um, flux * u.Jy)


def edge_problem() -> FittingProblem:
    observed = Spectrum(
        EDGE_GRID * u.um, 2.0 * EDGE_GRID * u.Jy, uncertainty=np.full(4, 0.1) * u.Jy
    )
    return FittingProblem(AtTheEdge(), [Dataset(observed)], seed=SEED)


class Polynomial(Model):
    """A polynomial with :data:`ENSEMBLE_SIZE_WARNING` coefficients, for the size row."""

    def __init__(self) -> None:
        self.register_buffer("wavelength", GRID, unit=u.um)
        for index in range(ENSEMBLE_SIZE_WARNING):
            self.register_parameter(Parameter(f"c{index}", st.norm(0.0, 1.0)))

    def evaluate(self, **values: Any) -> Spectrum:
        ctx = self.context(values)
        x = ctx["wavelength"] / 10.0
        flux = sum(ctx[f"c{index}"] * x**index for index in range(ENSEMBLE_SIZE_WARNING))
        return Spectrum(ctx["wavelength"] * u.um, flux * u.Jy)


def size_problem() -> FittingProblem:
    observed = Spectrum(
        GRID * u.um, np.ones(GRID.size) * u.Jy, uncertainty=np.ones(GRID.size) * u.Jy
    )
    return FittingProblem(Polynomial(), [Dataset(observed)], seed=SEED)


# -- the user-journeys memo's Appendix A, copied from its walkthrough ---------
# (docs/design/walkthroughs/persona_a.py: nine bands, a modified blackbody at
# T = 180, beta = 1.6, scale = 5, 8 % Gaussian noise, three free parameters)

APPENDIX_A_FILTERS = [
    "2MASS_Ks",
    "WISE_RSR_W1",
    "WISE_RSR_W2",
    "WISE_RSR_W3",
    "WISE_RSR_W4",
    "SPITZER_MIPS_24",
    "SPITZER_MIPS_70",
    "HERSCHEL_PACS_100",
    "HERSCHEL_PACS_160",
]
APPENDIX_A_GRID = np.geomspace(1.5, 250.0, 600)
APPENDIX_A_SEED = 20261006
#: The memo's first budget: 24 walkers, 1500 steps, 500 of them burn-in.
APPENDIX_A_WALKERS, APPENDIX_A_STEPS, APPENDIX_A_BURN_IN = 24, 1500, 500


def _appendix_a_model() -> ModifiedBlackBody:
    return ModifiedBlackBody(
        APPENDIX_A_GRID,
        temperature=st.loguniform(30, 1500),
        beta=st.uniform(0.5, 2.0),
        scale=st.loguniform(1e-3, 1e3),
        reference_wavelength=100.0,
    )


def appendix_a_problem() -> FittingProblem:
    step = SyntheticPhotometry.from_library(APPENDIX_A_FILTERS, APPENDIX_A_GRID)
    instrument = Instrument([step], channel="default", label="phot")
    truth = _appendix_a_model().compile_for(negotiate([instrument]))
    clean = instrument(truth(temperature=180.0, beta=1.6, scale=5.0))
    rng = np.random.default_rng(1)
    sigma = 0.08 * np.abs(clean.values)
    observed = PhotometricPoints(
        clean.filters,
        clean.spectral_axis.values * u.um,
        (clean.values + rng.normal(0, sigma)) * u.Jy,
        uncertainty=sigma * u.Jy,
    )
    likelihood = Likelihood(GaussianFamily(), IndependentNoise())
    return FittingProblem(
        _appendix_a_model(),
        {"phot": Dataset(observed, instrument, likelihood=likelihood)},
        seed=APPENDIX_A_SEED,
    )


def _caught(kind: type[Warning], record: list[warnings.WarningMessage]) -> list[str]:
    return [str(item.message) for item in record if issubclass(item.category, kind)]


# ---------------------------------------------------------------------------
# The regression row: Appendix A
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def appendix_a_default() -> Any:
    """The memo's first budget from the default start (the optimiser's mode)."""
    engine = EmceeEngine(appendix_a_problem(), walkers=APPENDIX_A_WALKERS)
    return engine.run(APPENDIX_A_STEPS, burn_in=APPENDIX_A_BURN_IN)


@pytest.fixture(scope="module")
def appendix_a_prior() -> Any:
    """The same budget from the prior: the former default."""
    engine = EmceeEngine(appendix_a_problem(), walkers=APPENDIX_A_WALKERS)
    with pytest.warns(ResultsWarning, match="initial='prior'"):
        return engine.run(APPENDIX_A_STEPS, burn_in=APPENDIX_A_BURN_IN, initial="prior")


def _bulk_ess(run: Any) -> dict[str, float]:
    verdict = check_convergence(run, rhat=np.inf, ess=np.inf)  # every variable listed
    return {name: ess for name, (_, ess) in verdict.failing.items()}


def _max_rhat(run: Any) -> dict[str, float]:
    verdict = check_convergence(run, rhat=np.inf, ess=np.inf)
    return {name: rhat for name, (rhat, _) in verdict.failing.items()}


@pytest.mark.heavy
class TestAppendixA:
    """The memo's first fit: R-hat 1.50 from the prior, below 1.1 from the default.

    What is pinned is what does not move with the platform: the attrs, the
    verdicts, the prior start's identity with ``initial_positions`` **in the
    same process**, and the ESS **ratio** between the two starts — not the
    draws' bytes (a stored digest of them failed on CI's runners, whose BLAS
    does not round as the machine that made it does; W7.0's row learnt the
    same lesson for a divergence count).
    """

    def test_the_default_start_converges_at_the_first_budget(self, appendix_a_default: Any) -> None:
        run = appendix_a_default
        assert run.attrs["ampere_start_kind"] == "optimum"
        assert run.attrs["ampere_start_route"] == "scipy"
        assert json.loads(run.attrs["ampere_start"])["route"] == "scipy"
        assert "ampere_start_fallback" not in run.attrs
        verdict = check_convergence(run)
        assert verdict.passed, str(verdict)

    def test_the_prior_start_is_the_former_default_bit_for_bit(self, appendix_a_prior: Any) -> None:
        """``initial="prior"`` is ``initial_positions(walkers)`` on the same stream, nothing else.

        The former default drew the walkers with ``initial_positions`` on the
        engine's initialisation stream; a run handed those draws as an array
        therefore walks the identical chain, in the same process.
        """
        expected = EmceeEngine(appendix_a_problem(), walkers=APPENDIX_A_WALKERS).initial_positions(
            APPENDIX_A_WALKERS
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ResultsWarning)
            supplied = EmceeEngine(appendix_a_problem(), walkers=APPENDIX_A_WALKERS).run(
                APPENDIX_A_STEPS, burn_in=APPENDIX_A_BURN_IN, initial=expected
            )
        posterior = appendix_a_prior["posterior"].dataset
        for name in posterior.data_vars:
            np.testing.assert_array_equal(
                np.asarray(posterior[name]), np.asarray(supplied["posterior"][name])
            )
        assert appendix_a_prior.attrs["ampere_start_kind"] == "prior"
        assert appendix_a_prior.attrs["ampere_start_route"] == "prior"
        assert "ampere_start" not in appendix_a_prior.attrs
        verdict = check_convergence(appendix_a_prior)
        assert not verdict.passed
        assert set(verdict.failing) == {"model.beta", "model.scale", "model.temperature"}

    def test_the_default_start_beats_the_prior_by_a_margin(
        self, appendix_a_default: Any, appendix_a_prior: Any
    ) -> None:
        """The memo's 1.50 / 44 against 1.06 / 300+, pinned as a ratio (measured 380 / 44 here)."""
        default_ess, prior_ess = _bulk_ess(appendix_a_default), _bulk_ess(appendix_a_prior)
        default_rhat, prior_rhat = _max_rhat(appendix_a_default), _max_rhat(appendix_a_prior)
        for name in prior_ess:
            assert default_ess[name] >= 3.0 * prior_ess[name], (name, default_ess, prior_ess)
            assert default_rhat[name] < prior_rhat[name], (name, default_rhat, prior_rhat)


# ---------------------------------------------------------------------------
# The default start and its attrs
# ---------------------------------------------------------------------------


class TestTheDefaultStart:
    @pytest.mark.parametrize("engine", [EmceeEngine, ZeusEngine])
    def test_the_default_is_the_optimisers_mode(self, engine: Any) -> None:
        run = engine(conjugate_problem(), walkers=8).run(20, burn_in=5)
        assert run.attrs["ampere_start_kind"] == "optimum"
        assert run.attrs["ampere_start_route"] == "scipy"
        record = json.loads(run.attrs["ampere_start"])
        assert record["route"] == "scipy"
        assert record["converged"] is True

    def test_the_default_start_is_reproducible_from_the_seed(self) -> None:
        one = EmceeEngine(conjugate_problem(), walkers=8).run(20, burn_in=5)
        two = EmceeEngine(conjugate_problem(), walkers=8).run(20, burn_in=5)
        np.testing.assert_array_equal(
            np.asarray(one["posterior"]["model.norm"]), np.asarray(two["posterior"]["model.norm"])
        )
        assert one.attrs["ampere_start"] == two.attrs["ampere_start"]

    def test_the_prior_start_draws_what_initial_positions_draws(self) -> None:
        """``initial="prior"`` is ``initial_positions(walkers)`` on the same stream."""
        expected = EmceeEngine(conjugate_problem(), walkers=8).initial_positions(8)
        prior = EmceeEngine(conjugate_problem(), walkers=8).run(6, initial="prior")
        supplied = EmceeEngine(conjugate_problem(), walkers=8).run(6, initial=expected)
        np.testing.assert_array_equal(
            np.asarray(prior["posterior"]["model.norm"]),
            np.asarray(supplied["posterior"]["model.norm"]),
        )
        assert prior.attrs["ampere_start_kind"] == "prior"
        assert "ampere_start" not in prior.attrs

    def test_a_supplied_array_and_optimum_say_supplied(self) -> None:
        from ampere.inference import optimise

        problem = conjugate_problem()
        engine = EmceeEngine(problem, walkers=8)
        run = engine.run(5, initial=np.full((8, 1), 2.0) + 0.01 * np.arange(8)[:, None])
        assert run.attrs["ampere_start_kind"] == "supplied"
        assert run.attrs["ampere_start_route"] == "user"
        run = engine.run(5, initial=optimise(problem, method="scipy", starts=2))
        assert run.attrs["ampere_start_kind"] == "supplied"
        assert run.attrs["ampere_start_route"] == "scipy"

    def test_an_unknown_start_name_is_refused(self) -> None:
        with pytest.raises(EngineError, match="'prior'"):
            EmceeEngine(conjugate_problem(), walkers=8).run(5, initial="optimum")

    def test_the_dynesty_run_says_prior(self) -> None:
        from ampere.inference import DynestyEngine

        run = DynestyEngine(conjugate_problem(), live_points=30).run(maxcall=400)
        assert run.attrs["ampere_start_kind"] == "prior"


class TestTheFallbacks:
    @pytest.mark.parametrize("engine", [EmceeEngine, ZeusEngine])
    def test_a_bound_pinned_optimum_starts_mixed_and_silent(self, engine: Any) -> None:
        with warnings.catch_warnings():
            warnings.simplefilter("error", DefaultStartWarning)
            from ampere.inference import BoundSaturationWarning

            warnings.simplefilter("error", BoundSaturationWarning)
            run = engine(edge_problem(), walkers=8).run(10)
        assert run.attrs["ampere_start_kind"] == "mixed"
        assert run.attrs["ampere_start_route"] == "scipy"
        assert json.loads(run.attrs["ampere_start"])["route"] == "scipy"
        assert run.attrs["ampere_start_fallback"] == "model.slope, model.width"

    def test_the_mixed_positions_are_prior_draws_on_the_saturated_columns_only(self) -> None:
        engine = EmceeEngine(edge_problem(), walkers=64)
        positions = engine._ensemble_start(None, 64, engine._checked_initial)
        labels = engine.problem.parameters.free_labels()
        column = {label: positions[:, labels.index(label)] for label in labels}
        assert np.ptp(column["model.slope"]) > 0.5
        assert np.std(column["model.slope"]) > 0.15
        assert np.ptp(column["model.width"]) > 0.5
        assert np.std(column["model.offset"]) < OFFSET_BALL_BOUND
        assert all(math.isfinite(engine.log_prob(row)) for row in positions)

    def test_a_mixed_run_is_reproducible_from_the_seed(self) -> None:
        first = EmceeEngine(edge_problem(), walkers=8).run(10)
        second = EmceeEngine(edge_problem(), walkers=8).run(10)
        np.testing.assert_array_equal(
            first["posterior"].dataset["model.slope"].values,
            second["posterior"].dataset["model.slope"].values,
        )

    def test_no_finite_start_falls_back_to_the_prior(self, monkeypatch: pytest.MonkeyPatch) -> None:
        from ampere.inference import _optimise

        def refuse(*args: Any, **kwargs: Any) -> Any:
            raise EngineError(
                "optimise('scipy'): none of the 1 starts ended at a point the problem can score."
            )

        monkeypatch.setattr(_optimise, "optimise", refuse)
        with pytest.warns(DefaultStartWarning, match="no start it could score") as record:
            run = EmceeEngine(conjugate_problem(), walkers=8).run(10)
        assert len(_caught(DefaultStartWarning, list(record))) == 1
        assert run.attrs["ampere_start_kind"] == "prior"
        assert "no start it could score" in run.attrs["ampere_start_fallback"]

    def test_the_prior_start_asked_for_is_silent(self) -> None:
        with warnings.catch_warnings():
            warnings.simplefilter("error", DefaultStartWarning)
            run = EmceeEngine(edge_problem(), walkers=8).run(10, initial="prior")
        assert "ampere_start_fallback" not in run.attrs


class TestTheSizeWarning:
    def test_a_large_problem_warns_once_per_run(self) -> None:
        problem = size_problem()
        assert problem.free_size == ENSEMBLE_SIZE_WARNING
        engine = EmceeEngine(problem)
        for _ in range(2):
            with warnings.catch_warnings(record=True) as record:
                warnings.simplefilter("always")
                engine.run(3, initial="prior")
            messages = _caught(EnsembleSizeWarning, record)
            assert len(messages) == 1
            assert "NUTSEngine" in messages[0] and "16 free coordinates" in messages[0]

    def test_a_small_problem_does_not(self) -> None:
        with warnings.catch_warnings():
            warnings.simplefilter("error", EnsembleSizeWarning)
            EmceeEngine(conjugate_problem(), walkers=8).run(3, initial="prior")


# ---------------------------------------------------------------------------
# The verdict
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def long_run() -> Any:
    return EmceeEngine(conjugate_problem(), walkers=16).run(2500, burn_in=500)


@pytest.fixture(scope="module")
def short_prior_run() -> Any:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ResultsWarning)
        return EmceeEngine(conjugate_problem(), walkers=8).run(12, burn_in=2, initial="prior")


class TestTheVerdict:
    def test_a_converged_run_passes_silently(self, long_run: Any) -> None:
        verdict = check_convergence(long_run)
        assert isinstance(verdict, ConvergenceVerdict)
        assert verdict.applicable and verdict.passed
        assert verdict.failing == {}
        assert str(verdict).startswith("Converged")
        with warnings.catch_warnings():
            warnings.simplefilter("error", ResultsWarning)
            summary(long_run)

    def test_one_failing_variable_is_named_with_its_numbers(self, short_prior_run: Any) -> None:
        verdict = check_convergence(short_prior_run)
        assert verdict.applicable and not verdict.passed
        assert list(verdict.failing) == ["model.norm"]
        rhat, ess = verdict.failing["model.norm"]
        assert not (rhat < 1.1 and ess >= 100)
        assert "initial='prior'" in verdict.remedy
        assert "model.norm (R-hat" in str(verdict)
        # one variable in all is not "one fails while the rest pass"
        assert "reparameterisation" not in verdict.remedy

    def test_summary_and_emission_warn_with_the_verdict(self) -> None:
        with pytest.warns(ResultsWarning, match="Not converged"):
            run = EmceeEngine(conjugate_problem(), walkers=8).run(12, burn_in=2, initial="prior")
        with pytest.warns(ResultsWarning, match="Not converged"):
            summary(run)

    def test_the_thresholds_are_the_callers(self, short_prior_run: Any) -> None:
        assert check_convergence(short_prior_run, rhat=np.inf, ess=0).passed

    def test_one_variable_failing_among_several_suggests_a_reparameterisation(
        self, long_run: Any
    ) -> None:
        tree = long_run.copy()
        posterior = tree["posterior"].dataset
        stuck = posterior["model.norm"].copy(
            data=np.cumsum(np.ones(posterior["model.norm"].shape), axis=1)
        )
        tree["posterior"] = posterior.assign(stuck=stuck)
        verdict = check_convergence(tree)
        assert list(verdict.failing) == ["stuck"]
        assert "reparameterisation of stuck" in verdict.remedy

    def test_an_approximate_run_is_never_warned(self, short_prior_run: Any) -> None:
        tree = short_prior_run.copy()
        tree.attrs["ampere_approximation"] = "mean_field"
        verdict = check_convergence(tree)
        assert not verdict.applicable
        assert verdict.passed
        with warnings.catch_warnings(record=True) as record:
            warnings.simplefilter("always")
            check_convergence(tree)
        assert not _caught(ResultsWarning, record)

    def test_an_approximate_emission_is_never_warned(self) -> None:
        """Every engine emits through ``Engine.finish``; a VI or SBI driver says it is approximate.

        Ten copies of eight prior draws are a "chain" that never moved: emitted
        as exact it fails the verdict and warns, emitted with the
        approximation VI writes it is not a chain and is silent.
        """
        engine = EmceeEngine(conjugate_problem(), walkers=8)
        draws = engine.initial_positions(8)[:, None, :].repeat(10, axis=1)
        with warnings.catch_warnings(record=True) as record:
            warnings.simplefilter("always")
            engine.finish(draws, extra_attrs={"approximation": "mean_field"})
        assert not [m for m in _caught(ResultsWarning, record) if "Not converged" in m]
        with warnings.catch_warnings(record=True) as record:
            warnings.simplefilter("always")
            engine.finish(draws)
        assert any("Not converged" in m for m in _caught(ResultsWarning, record))
