"""The engine battery, and its first customers: nautilus and ultranest.

**W5.14.** The inference-extensions memo (``docs/design/
inference_extensions_memo.md`` §6) names the battery that should come with any
new engine, and this module is it:

1. **the run is a run** — every engine emits the ArviZ ``DataTree`` with the
   per-draw split, the provenance attrs and ``ampere_approximation = "none"``
   (these are exact samplers), and a nested sampler's weighted output is
   resampled to equal weight with the original count recorded
   (``results.md`` §9, W5.0);
2. **the evidence is checked against a closed form** — not against another
   sampler. The problem is the two-parameter conjugate linear-Gaussian case
   W3.6's calibration rows use in one dimension, widened by one parameter
   because nautilus refuses a one-dimensional problem: a straight line in the
   flux with Gaussian priors on both coefficients and known Gaussian noise has
   a marginal likelihood in closed form, and :func:`analytic` writes it out
   rather than estimating it from a long reference run;
3. **every engine is SBC-ranked** through :func:`ampere.results.calibration.
   sbc` — the honest, expensive route (Talts et al. 2018), at a per-PR budget
   here and at a budget with power behind ``-m engines_full``;
4. **one cost record per run** — the memo's design horizon (g). It is
   ``ampere_engine_evaluations``, which :meth:`ampere.inference.engine.Engine.
   finish` writes for every run from the evaluation cache, beside the
   sampler's own count of likelihood calls under the engine's own name.

``DynestyEngine`` is in the table as the **control**, not because it is new.
Three nested samplers agreeing with each other proves nothing if all three are
wrong in the same way; three agreeing with a closed form, one of which has
been in the gate since W2.2, is a check on the two that have not.

Budgets and skips
-----------------
Every run here is small and seeded, and the module is a couple of minutes, for
the reason ``tests/inference/test_engines.py``'s own header gives. The
tolerances are stated in units of each engine's *own* reported uncertainty
with an absolute floor, and are several times the run-to-run spread measured
while the drivers were written (recorded in W5.14's report): loose enough that
a correct sampler passes essentially always, and far tighter than a wrong
evidence, which moves ``log Z`` by whole nats rather than by tenths.

nautilus and ultranest are each behind their own extra, so each engine's rows
skip independently when that package is absent — which is what the ``dev``
environment's ``test-all`` does, while the ``nested`` environment runs them.
"""

from __future__ import annotations

import dataclasses
import importlib.util
import math
import warnings
from collections.abc import Callable, Mapping
from typing import Any

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.core import Dataset, FittingProblem, Model, Parameter, Spectrum
from ampere.inference import (
    DynestyEngine,
    EngineError,
    NautilusEngine,
    UltranestEngine,
)
from ampere.results import sbc

# ---------------------------------------------------------------------------
# The conjugate problem, and its closed form
# ---------------------------------------------------------------------------

SEED = 20260922
GRID = np.linspace(1.0, 4.0, 16)
SIGMA = 0.3
#: (mean, sd) of the independent Gaussian priors on the two coefficients.
PRIOR_MEAN = np.array([2.0, 0.0])
PRIOR_SD = np.array([0.5, 0.5])
TRUTH = np.array([2.0, 0.2])

DATA = TRUTH[0] * GRID + TRUTH[1] + np.random.default_rng(SEED).normal(0.0, SIGMA, GRID.size)
OBSERVED = Spectrum(GRID * u.um, DATA * u.Jy, uncertainty=np.full(GRID.size, SIGMA) * u.Jy)


class Line(Model):
    """``flux = slope * wavelength + offset`` — linear in both free parameters.

    Linearity is the whole point: with Gaussian priors and known Gaussian
    noise the posterior *and* the marginal likelihood are both conjugate, so
    :func:`analytic` is an oracle rather than another sampler's opinion.
    """

    def __init__(self, wavelength: np.ndarray) -> None:
        self.register_buffer("wavelength", wavelength, unit=u.um)
        self.register_parameter(Parameter("slope", st.norm(PRIOR_MEAN[0], PRIOR_SD[0])))
        self.register_parameter(Parameter("offset", st.norm(PRIOR_MEAN[1], PRIOR_SD[1])))

    def evaluate(self, **values: Any) -> Spectrum:
        ctx = self.context(values)
        return Spectrum(
            ctx["wavelength"] * u.um,
            (ctx["slope"] * ctx["wavelength"] + ctx["offset"]) * u.Jy,
        )


@dataclasses.dataclass(frozen=True)
class Closed:
    """The closed-form answers: the marginal likelihood and the posterior."""

    log_evidence: float
    mean: np.ndarray
    sd: np.ndarray


def analytic(grid: np.ndarray = GRID, data: np.ndarray = DATA) -> Closed:
    r"""The conjugate answers, written out.

    With design matrix ``X = [x, 1]``, prior ``beta ~ N(mu, Lambda)`` and
    likelihood ``y ~ N(X beta, sigma**2 I)``:

    * the posterior precision is ``Lambda**-1 + X' X / sigma**2`` and the
      posterior mean is that inverse times ``Lambda**-1 mu + X' y / sigma**2``;
    * the marginal likelihood is ``y ~ N(X mu, sigma**2 I + X Lambda X')``,
      which is the evidence every nested sampler here is held to.
    """
    design = np.stack([grid, np.ones_like(grid)], axis=1)
    prior_cov = np.diag(PRIOR_SD**2)
    marginal = SIGMA**2 * np.eye(grid.size) + design @ prior_cov @ design.T
    residual = data - design @ PRIOR_MEAN
    _, logdet = np.linalg.slogdet(marginal)
    log_evidence = -0.5 * (
        grid.size * math.log(2.0 * math.pi)
        + logdet
        + float(residual @ np.linalg.solve(marginal, residual))
    )
    precision = np.linalg.inv(prior_cov) + design.T @ design / SIGMA**2
    covariance = np.linalg.inv(precision)
    mean = covariance @ (np.linalg.solve(prior_cov, PRIOR_MEAN) + design.T @ data / SIGMA**2)
    return Closed(log_evidence=log_evidence, mean=mean, sd=np.sqrt(np.diag(covariance)))


def evidence_problem(seed: int | None = SEED) -> FittingProblem:
    """Two free parameters, an exactly known posterior and an exact evidence."""
    return FittingProblem(Line(GRID), [Dataset(OBSERVED)], seed=seed)


# A second, deliberately tiny problem for the calibration rows: SBC costs one
# full fit per simulation, so its problem must be cheap rather than precise.
SBC_GRID = np.linspace(1.0, 3.0, 6)
SBC_OBSERVED = Spectrum(
    SBC_GRID * u.um,
    (TRUTH[0] * SBC_GRID + TRUTH[1]) * u.Jy,
    uncertainty=np.full(SBC_GRID.size, SIGMA) * u.Jy,
)


def calibration_problem(seed: int | None = SEED) -> FittingProblem:
    """The simulating problem the SBC rows draw truths and data from."""
    return FittingProblem(Line(SBC_GRID), [Dataset(SBC_OBSERVED)], seed=seed)


# ---------------------------------------------------------------------------
# The engines under test, and the budgets they run at
# ---------------------------------------------------------------------------


def _installed(module: str) -> bool:
    return importlib.util.find_spec(module) is not None


@dataclasses.dataclass(frozen=True)
class Case:
    """One engine, at the two budgets the battery runs it at."""

    name: str
    build: Callable[..., Any]
    #: Constructor keywords at the (larger) evidence budget.
    settings: Mapping[str, Any]
    #: ``run()`` keywords at the evidence budget.
    run_options: Mapping[str, Any]
    #: Constructor keywords at the (tiny) calibration budget.
    sbc_settings: Mapping[str, Any]
    #: ``run()`` keywords at the calibration budget.
    sbc_run_options: Mapping[str, Any]
    #: The engine's own record of how many times it called the likelihood.
    calls_attr: str
    module: str

    def engine(self, problem: FittingProblem) -> Any:
        return self.build(problem, **self.settings)

    def sbc_engine(self, problem: FittingProblem) -> Any:
        return self.build(problem, **self.sbc_settings)


CASES: tuple[Case, ...] = (
    Case(
        name="dynesty",
        build=DynestyEngine,
        settings={"live_points": 200},
        run_options={},
        sbc_settings={"live_points": 40},
        sbc_run_options={"maxcall": 4000},
        calls_attr="ampere_dynesty_ncall",
        module="dynesty",
    ),
    Case(
        name="nautilus",
        build=NautilusEngine,
        # `n_networks=0` turns the neural boundary off and falls back to plain
        # ellipsoidal sampling. It is the right setting for a two-dimensional
        # conjugate toy — the network has nothing to learn that an ellipsoid
        # does not already describe — and it halves the module's wall time.
        # The network path is exercised once, in its own row below.
        settings={"live_points": 200, "options": {"n_networks": 0}},
        run_options={"n_eff": 2000},
        sbc_settings={"live_points": 40, "options": {"n_networks": 0}},
        sbc_run_options={"n_eff": 150},
        calls_attr="ampere_nautilus_likelihood_calls",
        module="nautilus",
    ),
    Case(
        name="ultranest",
        build=UltranestEngine,
        settings={"live_points": 200},
        run_options={"min_ess": 1000},
        sbc_settings={"live_points": 40},
        sbc_run_options={"min_ess": 100, "dlogz": 1.0},
        calls_attr="ampere_ultranest_likelihood_calls",
        module="ultranest",
    ),
)

#: The two W5.14 added, for the rows that are about *them* rather than about
#: nested sampling in general.
NEW = tuple(case for case in CASES if case.name != "dynesty")


def parameters(cases: tuple[Case, ...] = CASES) -> list[Any]:
    """*cases* as pytest params, each skipping on its own missing extra."""
    return [
        pytest.param(
            case,
            id=case.name,
            marks=pytest.mark.skipif(
                not _installed(case.module),
                reason=(
                    f"{case.module} is not installed here; it is behind the "
                    f"'{case.module}' extra (pixi environment 'nested')"
                ),
            ),
        )
        for case in cases
    ]


#: Tolerances. In units of the engine's own reported evidence uncertainty,
#: with an absolute floor in nats — the floor is what keeps nautilus's small
#: reported uncertainty (`1/sqrt(n_eff)`, of order 0.02 at this budget) from
#: making the row a test of that estimator rather than of the evidence. The
#: deviations measured while writing the drivers were -0.15, -0.26 and +1.46
#: sigma, i.e. 0.03, 0.01 and 0.22 nats.
EVIDENCE_SIGMA = 4.0
EVIDENCE_FLOOR = 0.4
#: Posterior agreement, in units of the analytic posterior's own standard
#: deviation, and as a ratio for the standard deviation itself. Measured
#: agreement was within 0.05 sd and 2 per cent respectively.
MEAN_TOLERANCE = 0.5
SD_TOLERANCE = 0.25


@pytest.fixture(scope="module")
def battery_runs() -> dict[str, Any]:
    """One run per available engine on the conjugate problem, computed once."""
    runs: dict[str, Any] = {}
    for case in CASES:
        if not _installed(case.module):
            continue
        runs[case.name] = case.engine(evidence_problem()).run(**case.run_options)
    return runs


def summary(run: Any, name: str) -> tuple[float, float]:
    column = np.asarray(run["posterior"][name]).reshape(-1)
    return float(column.mean()), float(column.std(ddof=1))


# ---------------------------------------------------------------------------
# 1. The run is a run
# ---------------------------------------------------------------------------


class TestEveryEngineEmitsTheRun:
    """The emission obligations, asked of each engine identically."""

    @pytest.mark.parametrize("case", parameters())
    def test_the_groups_and_the_per_draw_split_are_there(
        self, case: Case, battery_runs: dict[str, Any]
    ) -> None:
        run = battery_runs[case.name]
        assert sorted(run.children) == [
            "constant_data",
            "log_likelihood",
            "observed_data",
            "posterior",
            "sample_stats",
        ]
        stats = run["sample_stats"].dataset
        for variable in ("lp", "log_prior", "log_likelihood"):
            assert variable in stats.data_vars
        assert run["posterior"]["model.slope"].sizes["chain"] == 1
        assert run["posterior"]["model.slope"].sizes["draw"] > 0

    @pytest.mark.parametrize("case", parameters())
    def test_the_run_knows_what_it_was_a_run_of(
        self, case: Case, battery_runs: dict[str, Any]
    ) -> None:
        attrs = battery_runs[case.name].attrs
        assert attrs["ampere_engine"] == case.name
        assert attrs["ampere_backend"] == "reference"
        assert len(attrs["ampere_spec_hash"]) == 32

    @pytest.mark.parametrize("case", parameters())
    def test_a_nested_sampler_is_an_exact_sampler(
        self, case: Case, battery_runs: dict[str, Any]
    ) -> None:
        """W5.0's ``ampere_approximation``: ``"none"``, and by construction."""
        assert battery_runs[case.name].attrs["ampere_approximation"] == "none"

    @pytest.mark.parametrize("case", parameters(NEW))
    def test_the_weighted_output_is_resampled_and_the_original_count_kept(
        self, case: Case, battery_runs: dict[str, Any]
    ) -> None:
        """``results.md`` §9's weighted-draw rule (a), obeyed identically."""
        run = battery_runs[case.name]
        attrs = run.attrs
        dead = attrs[f"ampere_{case.name}_dead_points"]
        equal = attrs[f"ampere_{case.name}_equal_weight_draws"]
        assert dead > 0
        assert equal == run["posterior"]["model.slope"].sizes["draw"]
        assert attrs[f"ampere_{case.name}_live_points"] == case.settings["live_points"]

    @pytest.mark.parametrize("case", parameters(NEW))
    def test_the_library_object_is_kept_for_the_raw_weighted_output(self, case: Case) -> None:
        """``engine.sampler`` is the library's own object, as for every driver."""
        engine = case.sbc_engine(evidence_problem())
        assert engine.sampler is None
        engine.run(**case.sbc_run_options)
        assert engine.sampler is not None


# ---------------------------------------------------------------------------
# 2. The evidence, against the closed form
# ---------------------------------------------------------------------------


class TestTheEvidenceMatchesTheClosedForm:
    """Accept criterion: "the evidence rows within the stated error"."""

    @pytest.mark.parametrize("case", parameters())
    def test_the_engine_neutral_triple_is_written(
        self, case: Case, battery_runs: dict[str, Any]
    ) -> None:
        attrs = battery_runs[case.name].attrs
        assert math.isfinite(attrs["ampere_log_evidence"])
        assert attrs["ampere_log_evidence_err"] > 0.0
        assert attrs["ampere_evidence_method"] == "nested_sampling"

    @pytest.mark.parametrize("case", parameters())
    def test_the_evidence_agrees_with_the_closed_form(
        self, case: Case, battery_runs: dict[str, Any]
    ) -> None:
        truth = analytic()
        attrs = battery_runs[case.name].attrs
        got = float(attrs["ampere_log_evidence"])
        err = float(attrs["ampere_log_evidence_err"])
        tolerance = max(EVIDENCE_SIGMA * err, EVIDENCE_FLOOR)
        assert abs(got - truth.log_evidence) < tolerance, (
            f"{case.name} reported log Z = {got:.4f} +- {err:.4f} against the closed-form "
            f"{truth.log_evidence:.4f}"
        )

    @pytest.mark.parametrize("case", parameters())
    @pytest.mark.parametrize("name", ["model.slope", "model.offset"])
    def test_the_posterior_agrees_with_the_closed_form(
        self, case: Case, name: str, battery_runs: dict[str, Any]
    ) -> None:
        """An evidence can be right for the wrong reason; the posterior cannot."""
        truth = analytic()
        index = ["model.slope", "model.offset"].index(name)
        mean, sd = summary(battery_runs[case.name], name)
        assert abs(mean - truth.mean[index]) < MEAN_TOLERANCE * truth.sd[index]
        assert abs(sd / truth.sd[index] - 1.0) < SD_TOLERANCE

    @pytest.mark.parametrize("case", parameters(NEW))
    def test_a_new_engine_agrees_with_dynesty_on_the_evidence(
        self, case: Case, battery_runs: dict[str, Any]
    ) -> None:
        """The cross-check the oracle cannot make: two independent estimators."""
        mine = float(battery_runs[case.name].attrs["ampere_log_evidence"])
        theirs = float(battery_runs["dynesty"].attrs["ampere_log_evidence"])
        assert abs(mine - theirs) < 2.0 * EVIDENCE_FLOOR


# ---------------------------------------------------------------------------
# 3. Calibration: every engine SBC-ranked
# ---------------------------------------------------------------------------

#: Per-PR budget. Below ``sbc``'s own goodness-of-fit power floor, and it says
#: so — the warning is allowed through, as every other reduced-budget
#: calibration row in this repository does. ``-m engines_full`` runs the
#: budget that has power.
SBC_COUNT = 8
SBC_FULL_COUNT = 100
SBC_DRAWS = 30


@pytest.fixture(scope="module")
def calibrations() -> dict[str, Any]:
    """One SBC study per available engine, computed once for the whole module.

    A fixture rather than a call per row: each study is *count* full nested-
    sampling fits, so computing it twice would double the module's wall time
    to assert two things about one object.
    """
    return {
        case.name: calibrate(case, count=SBC_COUNT) for case in CASES if _installed(case.module)
    }


def calibrate(case: Case, *, count: int) -> Any:
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=r".*goodness-of-fit test.*")
        return sbc(
            calibration_problem(),
            case.sbc_engine,
            count=count,
            draws=SBC_DRAWS,
            run_options=dict(case.sbc_run_options),
            seed=515151,
            label=f"W5.14 engine battery, {case.name}",
        )


class TestEveryEngineIsRanked:
    """The memo's battery: "every new engine SBC-ranked through calibration.sbc"."""

    @pytest.mark.parametrize("case", parameters())
    def test_the_ranks_are_produced_and_shaped(
        self, case: Case, calibrations: dict[str, Any]
    ) -> None:
        calibration = calibrations[case.name]
        ranks = np.asarray(calibration["ranks"].values)
        assert ranks.shape == (SBC_COUNT, 2)
        assert ranks.min() >= 0
        assert ranks.max() <= SBC_DRAWS
        assert sorted(str(name) for name in calibration["parameter"].values) == [
            "model.offset",
            "model.slope",
        ]
        assert case.build.__name__ in calibration.attrs["ampere_calibration_engine"]

    @pytest.mark.parametrize("case", parameters())
    def test_nothing_failed_and_the_uniformity_test_ran(
        self, case: Case, calibrations: dict[str, Any]
    ) -> None:
        calibration = calibrations[case.name]
        assert calibration.attrs["ampere_calibration_failures"] == 0
        pvalues = np.asarray(calibration["ks_pvalue"].values)
        assert pvalues.shape == (2,)
        assert np.all((pvalues >= 0.0) & (pvalues <= 1.0))

    @pytest.mark.engines_full
    @pytest.mark.parametrize("case", parameters())
    def test_the_ranks_are_uniform_at_a_budget_with_power(self, case: Case) -> None:
        """The row the per-PR budget cannot make: Talts et al.'s actual test."""
        calibration = calibrate(case, count=SBC_FULL_COUNT)
        pvalues = np.asarray(calibration["ks_pvalue"].values)
        assert np.all(pvalues > 0.01), f"{case.name}'s ranks are not uniform: KS p-values {pvalues}"


# ---------------------------------------------------------------------------
# 4. The cost record (design horizon (g))
# ---------------------------------------------------------------------------


class TestTheCostRecord:
    """One cost record per run, which is what horizon (g)'s hook asks for."""

    @pytest.mark.parametrize("case", parameters())
    def test_every_run_records_its_evaluations(
        self, case: Case, battery_runs: dict[str, Any]
    ) -> None:
        attrs = battery_runs[case.name].attrs
        assert attrs["ampere_engine_evaluations"] > 0
        assert attrs["ampere_engine_draws_recomputed"] >= 0

    @pytest.mark.parametrize("case", parameters())
    def test_the_sampler_own_call_count_is_recorded_beside_it(
        self, case: Case, battery_runs: dict[str, Any]
    ) -> None:
        """The two are not the same number, and the difference is the cache."""
        attrs = battery_runs[case.name].attrs
        assert attrs[case.calls_attr] > 0

    @pytest.mark.parametrize("case", parameters(NEW))
    def test_the_engine_records_its_library_version(
        self, case: Case, battery_runs: dict[str, Any]
    ) -> None:
        version = battery_runs[case.name].attrs[f"ampere_{case.name}_version"]
        assert isinstance(version, str) and version != "unknown"


# ---------------------------------------------------------------------------
# 5. Reproducibility, and the refusals
# ---------------------------------------------------------------------------


class TestReproducibility:
    @pytest.mark.parametrize("case", parameters(NEW))
    def test_the_same_seed_gives_the_same_draws(self, case: Case) -> None:
        """Including ultranest, whose library draws from numpy's global state."""
        first = case.sbc_engine(evidence_problem()).run(**case.sbc_run_options)
        second = case.sbc_engine(evidence_problem()).run(**case.sbc_run_options)
        assert np.allclose(
            np.asarray(first["posterior"]["model.slope"]),
            np.asarray(second["posterior"]["model.slope"]),
        )

    @pytest.mark.skipif(not _installed("ultranest"), reason="ultranest is not installed here")
    def test_ultranest_restores_the_global_generator_it_seeded(self) -> None:
        """``global_seed``'s contract: seed, run, put back."""
        np.random.seed(4242)
        before = np.random.get_state()[1][:8].copy()  # type: ignore[index]
        UltranestEngine(evidence_problem(), live_points=40).run(min_ess=100, dlogz=1.0)
        after = np.random.get_state()[1][:8].copy()  # type: ignore[index]
        assert np.array_equal(before, after)


class TestRefusals:
    @pytest.mark.parametrize("case", parameters(NEW))
    def test_too_few_live_points_is_refused(self, case: Case) -> None:
        with pytest.raises(EngineError, match="at least 2 live points"):
            case.build(evidence_problem(), live_points=1)

    @pytest.mark.skipif(not _installed("nautilus"), reason="nautilus is not installed here")
    def test_nautilus_refuses_a_one_dimensional_problem_by_name(self) -> None:
        """The library raises a bare ValueError; this driver says what to do."""

        class OneParameter(Model):
            def __init__(self, wavelength: np.ndarray) -> None:
                self.register_buffer("wavelength", wavelength, unit=u.um)
                self.register_parameter(Parameter("slope", st.norm(2.0, 0.5)))

            def evaluate(self, **values: Any) -> Spectrum:
                ctx = self.context(values)
                return Spectrum(ctx["wavelength"] * u.um, (ctx["slope"] * ctx["wavelength"]) * u.Jy)

        problem = FittingProblem(OneParameter(GRID), [Dataset(OBSERVED)], seed=SEED)
        with pytest.raises(EngineError, match="fewer than 2 free parameters"):
            NautilusEngine(problem)

    @pytest.mark.skipif(not _installed("nautilus"), reason="nautilus is not installed here")
    def test_the_neural_boundary_path_runs_too(self) -> None:
        """The battery's budget turns it off; something has to turn it on."""
        engine = NautilusEngine(evidence_problem(), live_points=60)
        run = engine.run(n_eff=200)
        assert run.attrs["ampere_nautilus_n_networks"] == 4
        assert math.isfinite(run.attrs["ampere_log_evidence"])


class TestTheDriversStayBackendNeutral:
    """The rule ``tests/inference/test_engines.py`` enforces for the namespace."""

    @pytest.mark.parametrize("case", parameters(NEW))
    def test_the_module_imports_no_backend(self, case: Case) -> None:
        import ampere.inference._nested as module

        source = module.__file__
        assert source is not None
        text = open(source, encoding="utf-8").read()
        assert "ampere.backends" not in text
