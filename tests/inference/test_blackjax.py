"""The blackjax route on jax: MCLMC as a sampler, Pathfinder twice over.

**W5.14**, the third of the item's three parts. The battery is
``test_nested.py``'s, asked of an engine whose two methods answer different
questions, so the rows say which method each claim is about:

1. **Each method emits the run.** The groups, the per-draw split, the
   provenance attrs, and the one root attribute that differs between them —
   ``ampere_approximation`` is ``"none"`` for MCLMC (a Markov chain, whose
   R-hat and ESS mean what they always mean) and ``"pathfinder"`` for
   Pathfinder, which also carries ``sample_stats.proposal_log_density``.
2. **Neither estimates an evidence, and neither pretends to.** The
   engine-neutral triple is *absent*, which is asserted rather than assumed:
   a placeholder under ``ampere_log_evidence`` would be read by a comparison
   that had no business making it. Pathfinder's ELBO is a single path's lower
   bound and is kept under the engine's own name.
3. **Both are SBC-ranked** through ``ampere.results.calibration.sbc``. This
   is where MCLMC's discretisation bias would show if the tuner's default
   energy-variance target were too loose for a problem this small — the
   unadjusted chain is the one engine here with no Metropolis correction, so
   ranking it is measurement rather than ceremony.
4. **One cost record per run**, design horizon (g): ``ampere_engine_
   evaluations`` from ``Engine.finish``, and the library's own integrator or
   L-BFGS count beside it.
5. **Bitwise reproducibility from the problem's seed.** jax has no global
   RNG, so there is nothing to fork and nothing to restore; two runs of the
   same seeded problem must agree exactly, and "exactly" here means
   ``array_equal``, not ``allclose``.

The oracle is the same one ``test_vi.py`` and ``test_nuts.py`` use: a
conjugate problem whose posterior is arithmetic, so a disagreement is a
disagreement with a closed form rather than with another approximation. The
two-parameter rows use an emcee run as the independent reference, which
reaches the same posterior by the contract path with no gradient at all.
"""

from __future__ import annotations

import importlib.util
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
from ampere.inference import BlackjaxEngine, EmceeEngine
from ampere.inference.exceptions import EngineError
from ampere.results.calibration import sbc


def _installed(module: str) -> bool:
    return importlib.util.find_spec(module) is not None


HAVE = _installed("blackjax") and _installed("jax")
if HAVE:
    import ampere.backends.jax as backend

    backend.configure_x64()

pytestmark = pytest.mark.skipif(
    not HAVE,
    reason="blackjax is not installed here; it is behind the 'blackjax' extra (pixi env 'jax')",
)

SEED = 20260910
REFERENCE_WAVELENGTH = 1.0
GRID = np.geomspace(1.0, 10.0, 20)
SIGMA = 0.2
INDEX = -1.0
PRIOR = (2.0, 0.5)


def power_law(grid: np.ndarray, norm: float, index: float) -> np.ndarray:
    return norm * (grid / REFERENCE_WAVELENGTH) ** index


def noisy(grid: np.ndarray, truth: np.ndarray, sigma: float, seed: int) -> Spectrum:
    rng = np.random.default_rng(seed)
    return Spectrum(
        grid * u.um,
        (truth + rng.normal(0.0, sigma, grid.size)) * u.Jy,
        uncertainty=np.full(grid.size, sigma) * u.Jy,
    )


DATA = noisy(GRID, power_law(GRID, 2.0, INDEX), SIGMA, seed=7)


def analytic() -> tuple[float, float]:
    """The posterior on ``norm`` with ``index`` held fixed, in closed form."""
    x = (GRID / REFERENCE_WAVELENGTH) ** INDEX
    y = np.asarray(DATA.values)
    mu, tau = PRIOR
    precision = 1.0 / tau**2 + float(np.sum(x**2)) / SIGMA**2
    mean = (mu / tau**2 + float(np.sum(x * y)) / SIGMA**2) / precision
    return mean, 1.0 / math.sqrt(precision)


def likelihood() -> Likelihood:
    return Likelihood(GaussianFamily(), backend.IndependentNoise())


def conjugate_problem(seed: int | None = SEED) -> FittingProblem:
    """One free parameter, and a posterior that is arithmetic."""
    return FittingProblem(
        backend.PowerLaw(
            GRID,
            norm=st.norm(*PRIOR),
            index=INDEX,
            reference_wavelength=REFERENCE_WAVELENGTH,
        ),
        [Dataset(DATA, likelihood=likelihood())],
        seed=seed,
    )


def correlated_problem(seed: int | None = SEED) -> FittingProblem:
    """Both parameters free, so the posterior is strongly correlated."""
    return FittingProblem(
        backend.PowerLaw(
            GRID,
            norm=st.norm(2.0, 1.0),
            index=st.norm(-1.0, 1.0),
            reference_wavelength=REFERENCE_WAVELENGTH,
        ),
        [Dataset(DATA, likelihood=likelihood())],
        seed=seed,
    )


#: Budgets. Small and fixed, so this belongs in the per-PR gate: MCLMC's cost
#: is its integrator steps and Pathfinder's is one L-BFGS path.
DRAWS = 1500
WARMUP = 1500
PATHFINDER_DRAWS = 4000

METHODS = ("mclmc", "pathfinder")


def problem_for(method: str) -> FittingProblem:
    """The problem each method's rows are made on, and why they differ.

    MCLMC needs **two** free parameters — blackjax refuses one and the driver
    refuses it first, by name — so the conjugate one-parameter problem, whose
    posterior is arithmetic, is available to Pathfinder and not to the
    sampler. That is not a gap in the battery: MCLMC's agreement is checked
    against an emcee run of the same two-parameter problem, which is an
    independent engine reaching the same posterior by the contract path with
    no gradient at all, and its *calibration* is checked by SBC, which needs
    no oracle.
    """
    return conjugate_problem() if method == "pathfinder" else correlated_problem()


def fit(problem: FittingProblem, method: str, **options: Any) -> Any:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        engine = BlackjaxEngine(problem, method=method)
        if method == "pathfinder":
            return engine.run(draws=PATHFINDER_DRAWS, **options)
        return engine.run(draws=DRAWS, warmup=WARMUP, **options)


@pytest.fixture(scope="module")
def runs() -> dict[str, Any]:
    """One run per method, on the problem that method can be run on."""
    return {method: fit(problem_for(method), method) for method in METHODS}


# ---------------------------------------------------------------------------
# 1. The run is a run
# ---------------------------------------------------------------------------


class TestEachMethodEmitsTheRun:
    @pytest.mark.parametrize("method", METHODS)
    def test_the_groups_and_the_per_draw_split_are_there(
        self, method: str, runs: dict[str, Any]
    ) -> None:
        run = runs[method]
        assert set(run.children) >= {"posterior", "sample_stats", "observed_data"}
        stats = run["sample_stats"].dataset
        for name in ("lp", "log_prior", "log_likelihood"):
            assert name in stats
        drawn = np.asarray(run["posterior"]["model.norm"])
        assert drawn.shape == (1, DRAWS if method == "mclmc" else PATHFINDER_DRAWS)

    @pytest.mark.parametrize("method", METHODS)
    def test_the_run_knows_what_it_was_a_run_of(self, method: str, runs: dict[str, Any]) -> None:
        attrs = runs[method].attrs
        assert attrs["ampere_engine"] == "blackjax"
        assert attrs["ampere_backend"] == "jax"
        assert attrs["ampere_blackjax_method"] == method
        assert attrs["ampere_realised"] == 1
        assert isinstance(attrs["ampere_blackjax_version"], str)
        assert attrs["ampere_blackjax_version"] != "unknown"

    def test_pathfinders_posterior_matches_the_closed_form(self, runs: dict[str, Any]) -> None:
        """The oracle row, on the one problem a closed form exists for."""
        mean, sd = analytic()
        drawn = np.asarray(runs["pathfinder"]["posterior"]["model.norm"]).ravel()
        assert float(drawn.mean()) == pytest.approx(mean, abs=0.3 * sd)
        assert float(drawn.std(ddof=1)) == pytest.approx(sd, rel=0.25)

    def test_mclmc_is_a_chain_and_says_so(self, runs: dict[str, Any]) -> None:
        """``ampere_approximation`` answers "do chain diagnostics apply?"."""
        attrs = runs["mclmc"].attrs
        assert attrs["ampere_approximation"] == "none"
        # ...and the bias the key is *not* about is recorded separately.
        # Stored as an integer, because that is what a netCDF attribute
        # is: the emission layer writes a bool as 0/1 like every other
        # engine's flag (``ampere_nuts_dense_mass``), so a reader tests
        # its truth rather than its identity.
        assert not attrs["ampere_blackjax_adjusted"]
        assert attrs["ampere_blackjax_desired_energy_var"] > 0.0
        assert attrs["ampere_blackjax_step_size"] > 0.0
        assert attrs["ampere_blackjax_L"] > 0.0

    def test_pathfinder_is_an_approximation_and_carries_its_density(
        self, runs: dict[str, Any]
    ) -> None:
        """W5.0's contract for an approximate engine, applied to Pathfinder."""
        run = runs["pathfinder"]
        assert run.attrs["ampere_approximation"] == "pathfinder"
        proposal = np.asarray(run["sample_stats"]["proposal_log_density"]).ravel()
        assert proposal.shape == (PATHFINDER_DRAWS,)
        assert np.all(np.isfinite(proposal))
        # A genuine per-draw density is never constant across thousands of
        # independent draws -- the guard W5.0's review found the need for.
        assert np.ptp(proposal) > 0.0

    def test_the_importance_weights_correct_the_approximation(self, runs: dict[str, Any]) -> None:
        """``results.md`` §9's formula, from the stored groups alone.

        Pathfinder's Gaussian is fitted along one L-BFGS path and is not the
        posterior; the weight built from ``log_prior + log_likelihood -
        proposal_log_density`` is what makes the difference correctable, and a
        density that was wrong by anything but an additive constant would move
        the reweighted mean away from the closed form rather than towards it.
        """
        run = runs["pathfinder"]
        stats = run["sample_stats"].dataset
        log_weight = (
            np.asarray(stats["log_prior"]).ravel()
            + np.asarray(stats["log_likelihood"]).ravel()
            - np.asarray(stats["proposal_log_density"]).ravel()
        )
        log_weight -= log_weight.max()
        weight = np.exp(log_weight)
        weight /= weight.sum()
        assert 1.0 / np.sum(weight**2) > PATHFINDER_DRAWS / 10.0

        mean, sd = analytic()
        drawn = np.asarray(run["posterior"]["model.norm"]).ravel()
        assert float(np.sum(weight * drawn)) == pytest.approx(mean, abs=0.2 * sd)

    def test_mclmc_agrees_with_emcee_on_a_correlated_posterior(self, runs: dict[str, Any]) -> None:
        """Two parameters, no closed form, and an independent gradient-free oracle."""
        run = runs["mclmc"]
        reference = EmceeEngine(correlated_problem(), walkers=16).run(steps=1500, burn_in=500)
        for name in ("model.norm", "model.index"):
            values = np.asarray(reference["posterior"][name]).ravel()
            drawn = np.asarray(run["posterior"][name]).ravel()
            assert float(drawn.mean()) == pytest.approx(
                float(values.mean()), abs=0.5 * float(values.std())
            )
            assert float(drawn.std(ddof=1)) == pytest.approx(float(values.std()), rel=0.3)


# ---------------------------------------------------------------------------
# 2. Evidence: not applicable, and said so
# ---------------------------------------------------------------------------


class TestNeitherMethodClaimsAnEvidence:
    """``test_nested.py``'s evidence section has no counterpart here, deliberately.

    MCLMC samples an unnormalised density like every other MCMC in this
    package, and Pathfinder's ELBO is a *lower bound* from a single L-BFGS
    path whose tightness is unknown. Writing either under
    ``ampere_log_evidence`` would put a number a reader is entitled to compare
    across engines beside a nested sampler's estimate, so the triple is absent
    and this class is the assertion that it stays absent.
    """

    @pytest.mark.parametrize("method", METHODS)
    def test_the_engine_neutral_triple_is_absent(self, method: str, runs: dict[str, Any]) -> None:
        attrs = runs[method].attrs
        for name in (
            "ampere_log_evidence",
            "ampere_log_evidence_err",
            "ampere_evidence_method",
        ):
            assert name not in attrs

    def test_pathfinders_elbo_is_kept_under_its_own_name(self, runs: dict[str, Any]) -> None:
        elbo = runs["pathfinder"].attrs["ampere_blackjax_pathfinder_elbo"]
        assert np.isfinite(float(elbo))


# ---------------------------------------------------------------------------
# 3. Calibration: both methods SBC-ranked
# ---------------------------------------------------------------------------

#: A deliberately tiny simulating problem: SBC costs one full fit per
#: simulation, and on jax it costs one compilation too.
SBC_GRID = np.linspace(1.0, 3.0, 6)
SBC_DATA = noisy(SBC_GRID, power_law(SBC_GRID, 2.0, INDEX), SIGMA, seed=11)

SBC_COUNT = 6
SBC_FULL_COUNT = 100
SBC_DRAWS = 30
SBC_RUN_OPTIONS: dict[str, dict[str, Any]] = {
    "mclmc": {"draws": 400, "warmup": 400},
    "pathfinder": {"draws": 400},
}


def calibration_problem(seed: int | None = SEED) -> FittingProblem:
    """Two free parameters, because MCLMC cannot be run over one."""
    return FittingProblem(
        backend.PowerLaw(
            SBC_GRID,
            norm=st.norm(*PRIOR),
            index=st.norm(INDEX, 0.3),
            reference_wavelength=REFERENCE_WAVELENGTH,
        ),
        [Dataset(SBC_DATA, likelihood=likelihood())],
        seed=seed,
    )


def calibrate(method: str, *, count: int) -> Any:
    def factory(problem: FittingProblem) -> Any:
        return BlackjaxEngine(problem, method=method)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return sbc(
            calibration_problem(),
            factory,
            count=count,
            draws=SBC_DRAWS,
            run_options=dict(SBC_RUN_OPTIONS[method]),
            seed=515151,
            label=f"W5.14 engine battery, blackjax {method}",
        )


@pytest.fixture(scope="module")
def calibrations() -> dict[str, Any]:
    """One SBC study per method, computed once for the whole module."""
    return {method: calibrate(method, count=SBC_COUNT) for method in METHODS}


class TestEachMethodIsRanked:
    @pytest.mark.parametrize("method", METHODS)
    def test_the_ranks_are_produced_and_shaped(
        self, method: str, calibrations: dict[str, Any]
    ) -> None:
        calibration = calibrations[method]
        ranks = np.asarray(calibration["ranks"].values)
        assert ranks.shape == (SBC_COUNT, 2)
        assert ranks.min() >= 0
        assert ranks.max() <= SBC_DRAWS
        assert "BlackjaxEngine" in calibration.attrs["ampere_calibration_engine"]

    @pytest.mark.parametrize("method", METHODS)
    def test_nothing_failed_and_the_uniformity_test_ran(
        self, method: str, calibrations: dict[str, Any]
    ) -> None:
        calibration = calibrations[method]
        assert calibration.attrs["ampere_calibration_failures"] == 0
        pvalues = np.asarray(calibration["ks_pvalue"].values)
        assert np.all((pvalues >= 0.0) & (pvalues <= 1.0))

    @pytest.mark.engines_full
    @pytest.mark.parametrize("method", METHODS)
    def test_the_ranks_are_uniform_at_a_budget_with_power(self, method: str) -> None:
        """The row the per-PR budget cannot make: Talts et al.'s actual test."""
        calibration = calibrate(method, count=SBC_FULL_COUNT)
        pvalues = np.asarray(calibration["ks_pvalue"].values)
        assert np.all(pvalues > 0.01), f"{method}'s ranks are not uniform: KS p-values {pvalues}"


# ---------------------------------------------------------------------------
# 4. The cost record (design horizon (g))
# ---------------------------------------------------------------------------


class TestTheCostRecord:
    @pytest.mark.parametrize("method", METHODS)
    def test_every_run_records_its_evaluations(self, method: str, runs: dict[str, Any]) -> None:
        attrs = runs[method].attrs
        assert attrs["ampere_engine_evaluations"] > 0
        assert attrs["ampere_engine_draws_recomputed"] >= 0

    def test_mclmc_records_the_integrator_steps_it_really_spent(self, runs: dict[str, Any]) -> None:
        """Tuning included: those were gradient evaluations that were spent."""
        attrs = runs["mclmc"].attrs
        assert attrs["ampere_blackjax_tuning_steps"] > 0
        assert attrs["ampere_blackjax_integrator_steps"] >= attrs["ampere_blackjax_tuning_steps"]
        assert attrs["ampere_blackjax_integrator_steps"] >= DRAWS

    def test_pathfinder_records_its_path_length(self, runs: dict[str, Any]) -> None:
        assert runs["pathfinder"].attrs["ampere_blackjax_lbfgs_iterations"] > 0


# ---------------------------------------------------------------------------
# 5. Reproducibility, the initialiser, and the refusals
# ---------------------------------------------------------------------------


class TestReproducibility:
    @pytest.mark.parametrize("method", METHODS)
    def test_the_same_seed_gives_bitwise_identical_draws(self, method: str) -> None:
        """Bitwise, not merely close: jax has no global RNG to leak through."""
        first = fit(problem_for(method), method)
        second = fit(problem_for(method), method)
        assert np.array_equal(
            np.asarray(first["posterior"]["model.norm"]),
            np.asarray(second["posterior"]["model.norm"]),
        )

    def test_a_different_seed_gives_different_draws(self) -> None:
        """The other half of the claim: the seed is doing the work."""
        first = fit(correlated_problem(seed=SEED), "mclmc")
        second = fit(correlated_problem(seed=SEED + 1), "mclmc")
        assert not np.array_equal(
            np.asarray(first["posterior"]["model.norm"]),
            np.asarray(second["posterior"]["model.norm"]),
        )


class TestPathfinderAsAnInitialiser:
    """The memo's second use for Pathfinder, and the reason it is not only a method."""

    def test_the_chains_start_from_the_approximation(self) -> None:
        run = fit(correlated_problem(), "mclmc", chains=2, initial="pathfinder")
        assert np.asarray(run["posterior"]["model.norm"]).shape == (2, DRAWS)
        reference = EmceeEngine(correlated_problem(), walkers=16).run(steps=1500, burn_in=500)
        values = np.asarray(reference["posterior"]["model.norm"]).ravel()
        drawn = np.asarray(run["posterior"]["model.norm"]).ravel()
        assert float(drawn.mean()) == pytest.approx(
            float(values.mean()), abs=0.5 * float(values.std())
        )

    def test_it_is_still_reproducible_from_the_problems_seed(self) -> None:
        first = fit(correlated_problem(), "mclmc", initial="pathfinder")
        second = fit(correlated_problem(), "mclmc", initial="pathfinder")
        assert np.array_equal(
            np.asarray(first["posterior"]["model.norm"]),
            np.asarray(second["posterior"]["model.norm"]),
        )


class TestRefusals:
    def test_an_unknown_method_is_refused_by_name(self) -> None:
        with pytest.raises(EngineError, match="does not know the method"):
            BlackjaxEngine(conjugate_problem(), method="nuts")

    def test_a_problem_on_another_backend_is_refused_by_name(self) -> None:
        from ampere.backends.reference import PowerLaw

        problem = FittingProblem(
            PowerLaw(GRID, norm=st.norm(*PRIOR), index=INDEX),
            [Dataset(DATA, likelihood=Likelihood(GaussianFamily()))],
            seed=SEED,
        )
        with pytest.raises(EngineError, match="cannot fit a problem"):
            BlackjaxEngine(problem)

    def test_mclmc_refuses_a_one_dimensional_problem_by_name(self) -> None:
        """blackjax refuses it one layer down; ampere says so first."""
        engine = BlackjaxEngine(conjugate_problem(), method="mclmc")
        with pytest.raises(EngineError, match="direction to decorrelate"):
            engine.run(draws=10)

    def test_pathfinder_refuses_a_second_chain(self) -> None:
        engine = BlackjaxEngine(conjugate_problem(), method="pathfinder")
        with pytest.raises(EngineError, match="no such thing as a second chain"):
            engine.run(draws=10, chains=2)

    def test_pathfinder_refuses_to_initialise_itself(self) -> None:
        engine = BlackjaxEngine(conjugate_problem(), method="pathfinder")
        with pytest.raises(EngineError, match="fit the same approximation twice"):
            engine.run(draws=10, initial="pathfinder")

    def test_an_unknown_named_initialiser_is_refused(self) -> None:
        # Two free parameters: MCLMC's own floor is checked before the
        # start point is, and a one-parameter problem would be refused
        # for the other reason.
        engine = BlackjaxEngine(correlated_problem(), method="mclmc")
        with pytest.raises(EngineError, match="only named initialiser"):
            engine.run(draws=10, initial="laplace")

    def test_a_zero_tuning_budget_is_refused(self) -> None:
        engine = BlackjaxEngine(correlated_problem(), method="mclmc")
        with pytest.raises(EngineError, match="tuning budget"):
            engine.run(draws=10, warmup=0)


class TestTheDriverStaysBackendNeutral:
    def test_the_module_imports_no_backend(self) -> None:
        """``architecture.md`` §4 rule 2: nothing at module level is a backend.

        ``test_nested.py``'s version of this row scans the whole file for
        ``ampere.backends``, which it can because that module carries no
        worked example. This one cannot: ``BlackjaxEngine``'s docstring
        *demonstrates a fit*, and a fit needs a problem, and a problem needs
        a backend — the same shape ``_nuts.py``'s and ``_vi.py``'s examples
        have. So the rule is checked where it lives, on the import lines,
        and the run-time half of it is ``tests/inference/test_engines.py``'s
        subprocess probe that importing the namespace imports no backend.
        """
        import ampere.inference._blackjax as module

        source = module.__file__
        assert source is not None
        text = open(source, encoding="utf-8").read()
        for line in text.splitlines():
            if line.startswith(("import ", "from ")):
                assert "blackjax" not in line
                assert "jax" not in line
