"""``SBIEngine`` end to end: the prior bridge, the summary, and three families.

Five things are under test, and the first is the one that would be easiest to
get wrong invisibly.

1. **The prior bridge is a real change of variables.** ``sbi`` is handed a
   torch ``Distribution`` over ℝⁿ whose ``sample`` draws θ in the constrained
   space and unconstrains it, and whose ``log_prob`` is
   ``ParameterSet.lnprior_unconstrained``. Those two must be the same
   distribution, and a driver that got the Jacobian wrong would still train,
   still sample, and still produce a posterior — a plausible one, and a wrong
   one. So the Jacobian is checked against a **numerically differentiated**
   ``constrain`` rather than against the analytic term the implementation
   itself uses, and the pushforward is checked against a closed form on a
   problem where one exists.

2. **The summary layout is fixed and consistent.** The observation the
   posterior is conditioned on and the rows the network trained on must have
   the same columns in the same order, masked samples dropped from both.

3. **NPE recovers a posterior an established sampler agrees with.** The
   ``test_engines`` joint problem, against an emcee reference run of the same
   problem — thresholds in units of the *reference's* own width, since the
   reference is where "right" is defined for a problem with no closed form.
   NLE and NRE are held to shape and finiteness at CI budgets, not to
   accuracy: both sample by MCMC, and an MCMC run short enough for a per-PR
   gate has no business being asked about its posterior's width.

4. **A black-box model on the reference backend is first class**, including
   once under a process pool — the case the whole engine exists for, run
   through ``examples/sbi/external_simulator.py`` rather than through a
   fixture that flatters it.

5. **The extra is optional, and the namespace stays clean.** Everything that
   does not need ``sbi`` runs in ``dev``: the refusals, the layout function,
   and the proof that importing ``ampere.inference`` imports neither torch nor
   ``sbi``. The rest skips.

Budgets and the trade-off
-------------------------
Every run here is small and seeded. The NPE agreement run is the expensive one
(a 2 000-draw budget and a full training run, about 40 s) and it is a
module-scoped fixture, so it happens once; everything else is tens of draws
and capped epochs. That costs statistical power, which is why the agreement
thresholds are stated in units of the emcee reference's own standard deviation
and set several times the Monte Carlo error at these budgets — loose enough
that a correct fit passes essentially always, and far tighter than the errors a
bridge bug produces (a dropped Jacobian moves a lognormal posterior's mean by
whole widths).
"""

from __future__ import annotations

import importlib
import importlib.util
import math
import subprocess
import sys
from pathlib import Path
from typing import Any

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.backends.reference import Resample
from ampere.core import (
    Dataset,
    FittingProblem,
    Instrument,
    Model,
    Parameter,
    ProcessExecutor,
    Spectrum,
)
from ampere.core.exceptions import OptionalDependencyError
from ampere.inference import EmceeEngine, EngineError, SBIEngine
from ampere.inference._sbi import SUMMARY_LAYOUT, _summary_of, _thinned

HAS_SBI = importlib.util.find_spec("sbi") is not None

#: One decorator rather than a skipif on every row, so that the reason a
#: reader sees is the same everywhere and the module still *collects* in
#: ``dev`` — which is what lets the refusal rows below run there.
needs_sbi = pytest.mark.skipif(
    not HAS_SBI,
    reason="needs the 'sbi' extra (pixi run -e sbi ...)",
)

SEED = 20260909
HERE = Path(__file__).resolve().parent
EXAMPLES = HERE.parents[1] / "examples" / "sbi"


def _sibling(name: str) -> Any:
    """Import a module beside this one (or an example) by *name*, not by path.

    The directory is left on ``sys.path`` afterwards, and that is the whole
    point rather than an oversight: a worker process has to be able to
    **re-import** ``external_simulator`` to unpickle the model it is handed, and
    a module loaded under a synthetic name from a file path cannot be
    re-imported anywhere. ``tests/examples/test_external_simulator.py`` records
    the same reasoning for the same module.
    """
    directory = str(HERE if name == "test_engines" else EXAMPLES)
    if directory not in sys.path:
        sys.path.insert(0, directory)
    return importlib.import_module(name)


def joint_problem() -> FittingProblem:
    """``test_engines``'s two-dataset joint fit -- imported, never copied."""
    return _sibling("test_engines").joint_problem()


# ---------------------------------------------------------------------------
# Problems
# ---------------------------------------------------------------------------


class Line(Model):
    """A one-parameter model with a *bounded* parameter, so the bridge matters."""

    def __init__(self, wavelength: np.ndarray, prior: Any = None) -> None:
        self.register_buffer("wavelength", wavelength, unit=u.um)
        self.register_parameter(
            Parameter("slope", st.lognorm(0.4, scale=2.0) if prior is None else prior)
        )

    def evaluate(self, **values: Any) -> Spectrum:
        context = self.context(values)
        return Spectrum(
            context["wavelength"] * u.um,
            context["slope"] * context["wavelength"] * u.Jy,
        )


GRID = np.array([1.0, 2.0, 3.0, 4.0])


def bounded_problem(seed: int | None = SEED, *, mask: Any = None) -> FittingProblem:
    """One positive-parameter dataset: ``slope`` is lognormal, so the reals are not its support."""
    observed = Spectrum(
        GRID * u.um,
        (2.0 * GRID) * u.Jy,
        uncertainty=np.full(GRID.size, 0.3) * u.Jy,
        mask=mask,
    )
    return FittingProblem(Line(GRID), [Dataset(observed)], seed=seed)


def two_dataset_problem(seed: int | None = SEED) -> FittingProblem:
    """Two datasets of genuinely different lengths, the longer one partly masked.

    The short one is reached through a ``Resample`` step, which is where a
    length change belongs (``likelihoods.md``: a likelihood compares samples
    index by index), so this is a real two-instrument problem rather than two
    views of one array.
    """
    short = Spectrum(GRID[:2] * u.um, (2.0 * GRID[:2]) * u.Jy, uncertainty=np.full(2, 0.3) * u.Jy)
    long = Spectrum(
        GRID * u.um,
        (2.0 * GRID) * u.Jy,
        uncertainty=np.full(GRID.size, 0.3) * u.Jy,
        mask=np.array([False, True, False, False]),
    )
    return FittingProblem(
        Line(GRID),
        {
            "short": Dataset(short, Instrument([Resample(GRID[:2])]), label="short"),
            "long": Dataset(long, label="long"),
        },
        seed=seed,
    )


# ---------------------------------------------------------------------------
# 1. The prior bridge
# ---------------------------------------------------------------------------


@needs_sbi
class TestThePriorBridge:
    """``sample`` and ``log_prob`` must be one distribution, not two."""

    @staticmethod
    def built(problem: FittingProblem) -> Any:
        import torch

        from ampere.inference._sbi import _prior_class

        return _prior_class(torch)(
            problem, problem.rng("test.prior"), dtype=torch.float32, device="cpu"
        )

    def test_it_has_the_shapes_a_torch_distribution_must_have(self) -> None:
        import torch

        problem = joint_problem()
        prior = self.built(problem)
        assert tuple(prior.event_shape) == (problem.free_size,)
        assert tuple(prior.batch_shape) == ()
        assert tuple(prior.sample().shape) == (problem.free_size,)
        assert tuple(prior.sample((5,)).shape) == (5, problem.free_size)
        drawn = prior.sample((7,))
        assert tuple(prior.log_prob(drawn).shape) == (7,)
        assert bool(torch.isfinite(prior.log_prob(drawn)).all())

    def test_sbi_accepts_it_unchanged(self) -> None:
        """``process_prior`` must not have to wrap it — the bridge is the contract."""
        from sbi.utils import process_prior

        problem = joint_problem()
        prior = self.built(problem)
        processed, dimension, returns_numpy = process_prior(prior)
        assert processed is prior
        assert dimension == problem.free_size
        assert returns_numpy is False

    def test_the_density_carries_a_jacobian_measured_independently(self) -> None:
        """The change-of-variables term, against a numerical derivative of ``constrain``.

        ``lnprior_unconstrained`` adds the analytic ``log_abs_det_jacobian`` of
        each parameter's bijection. Checking it against that same analytic term
        would check nothing, so it is checked against a central difference of
        ``constrain`` itself: the bijections are per-parameter and therefore
        diagonal, so the log determinant is the sum of the logs of the diagonal.
        """
        problem = joint_problem()
        prior = self.built(problem)
        rng = np.random.default_rng(4)
        for _ in range(6):
            y = np.asarray(prior.sample().numpy(), dtype=float) + rng.normal(0.0, 0.2, 3)
            step = 1e-6
            diagonal = []
            for index in range(problem.free_size):
                plus, minus = y.copy(), y.copy()
                plus[index] += step
                minus[index] -= step
                derivative = (problem.constrain(plus)[index] - problem.constrain(minus)[index]) / (
                    2.0 * step
                )
                diagonal.append(math.log(abs(derivative)))
            expected = problem.parameters.lnprior(problem.constrain(y)) + sum(diagonal)
            assert float(prior.log_prob(_as_tensor(y))) == pytest.approx(expected, rel=1e-4)

    def test_the_draws_are_distributed_as_the_density_says(self) -> None:
        """The pushforward, against a closed form rather than against itself.

        ``slope ~ lognorm(s=0.4, scale=2.0)`` unconstrains to
        ``log(slope) ~ Normal(log 2.0, 0.4)`` exactly, so both halves of the
        bridge have an independent oracle: the density is checked pointwise
        against ``scipy``'s normal, and the draws against that normal's mean
        and standard deviation. A dropped Jacobian passes neither.
        """
        problem = bounded_problem()
        prior = self.built(problem)
        oracle = st.norm(math.log(2.0), 0.4)
        for y in (-0.5, 0.0, 0.3, 1.1):
            assert float(prior.log_prob(_as_tensor([y]))) == pytest.approx(
                float(oracle.logpdf(y)), rel=1e-5
            )
        drawn = np.asarray(prior.sample((4000,)).numpy(), dtype=float).reshape(-1)
        # 4000 draws: the standard error on the mean is 0.4/sqrt(4000) = 0.006,
        # so 0.05 is eight of them and a dropped bijection is many more.
        assert float(drawn.mean()) == pytest.approx(math.log(2.0), abs=0.05)
        assert float(drawn.std()) == pytest.approx(0.4, abs=0.05)

    def test_a_draw_maps_back_into_the_constrained_support(self) -> None:
        """Every unconstrained draw constrains to a point the prior scores finitely."""
        problem = bounded_problem()
        prior = self.built(problem)
        for row in np.asarray(prior.sample((50,)).numpy(), dtype=float):
            constrained = problem.constrain(row)
            assert float(constrained[0]) > 0.0
            assert math.isfinite(problem.parameters.lnprior(constrained))


def _as_tensor(values: Any) -> Any:
    import torch

    return torch.as_tensor(np.asarray(values, dtype=float).reshape(-1), dtype=torch.float32)


# ---------------------------------------------------------------------------
# 2. The summary layout (no extra needed)
# ---------------------------------------------------------------------------


class TestTheSummaryLayout:
    """``_summary_of``: the one function W3.3 replaces, checked without ``sbi``."""

    def test_the_observation_is_the_datasets_concatenated_in_order(self) -> None:
        problem = two_dataset_problem()
        summary = _summary_of(
            {label: problem.datasets[label].observed for label in problem.datasets},
            problem.datasets,
            batched=False,
        )
        # 'short' has 2 samples; 'long' has 4 with one masked out.
        assert summary.shape == (1, 5)
        assert list(problem.datasets) == ["short", "long"]
        assert summary[0, :2].tolist() == [2.0, 4.0]

    def test_the_mask_drops_the_same_column_from_a_batch(self) -> None:
        """The consistency the layout's name promises: same columns, both sides."""
        problem = two_dataset_problem()
        batch = problem.simulate_many(3, observe=True)
        rows = _summary_of(batch.observations, problem.datasets, batched=True)
        observation = _summary_of(
            {label: problem.datasets[label].observed for label in problem.datasets},
            problem.datasets,
            batched=False,
        )
        assert rows.shape == (3, observation.shape[1])
        # The masked sample of 'long' is index 1 of that dataset, so the
        # feature vector is short[0:2] + long[0], long[2], long[3].
        assert rows.shape[1] == 5

    def test_the_layout_has_a_name_a_run_can_record(self) -> None:
        assert SUMMARY_LAYOUT == "flat"

    def test_a_complex_container_is_refused_by_name(self) -> None:
        """W3.3's contract owns the real/imaginary encoding, so this one refuses."""
        problem = bounded_problem()

        class Complex:
            values = np.array([1.0 + 2.0j, 3.0 + 0.0j])

        with pytest.raises(EngineError, match="complex"):
            _summary_of({"default": Complex()}, problem.datasets, batched=False)

    def test_a_missing_observation_is_refused_by_name(self) -> None:
        problem = bounded_problem()
        with pytest.raises(EngineError, match="no observation"):
            _summary_of({}, problem.datasets, batched=False)

    def test_a_fully_masked_problem_is_refused_rather_than_summarised_to_nothing(
        self,
    ) -> None:
        problem = bounded_problem(mask=np.ones(GRID.size, dtype=bool))
        with pytest.raises(EngineError, match="empty"):
            _summary_of(
                {label: problem.datasets[label].observed for label in problem.datasets},
                problem.datasets,
                batched=False,
            )


class TestTheThinnedTrace:
    """The loss trace is a provenance attribute, so it has to stay small."""

    def test_a_short_trace_is_kept_whole_at_stride_one(self) -> None:
        values, stride = _thinned([1.0, 2.0, 3.0])
        assert (values, stride) == ([1.0, 2.0, 3.0], 1)

    def test_a_long_trace_is_strided_and_keeps_its_last_value(self) -> None:
        raw = [float(index) for index in range(5000)]
        values, stride = _thinned(raw)
        assert stride == 25
        assert len(values) <= 202
        assert values[0] == 0.0
        assert values[-1] == 4999.0


# ---------------------------------------------------------------------------
# 3. Refusals and the optional extra (no extra needed)
# ---------------------------------------------------------------------------


class TestWhatItRefuses:
    """Every refusal names what was wrong and what to do — checked, not asserted."""

    def test_an_unknown_method_is_refused_with_the_three_that_exist(self) -> None:
        with pytest.raises(EngineError, match="npe"):
            SBIEngine(bounded_problem(), method="snpe")

    def test_a_context_is_refused_by_name_as_a_later_item(self) -> None:
        """The reserved slot: the signature exists, the machinery does not."""
        with pytest.raises(EngineError, match="reserved"):
            SBIEngine(bounded_problem(), context={"sigma": 0.1})

    def test_a_budget_below_one_is_refused(self) -> None:
        with pytest.raises(EngineError, match="budget"):
            SBIEngine(bounded_problem(), budget=0)

    def test_a_round_count_below_one_is_refused(self) -> None:
        with pytest.raises(EngineError, match="round"):
            SBIEngine(bounded_problem(), rounds=0)

    def test_a_builder_and_an_embedding_together_are_refused(self) -> None:
        """Two statements about one network is a mistake, not a configuration."""
        with pytest.raises(EngineError, match="builder"):
            SBIEngine(bounded_problem(), density_estimator=lambda *a, **k: None, embedding="FC")

    def test_a_nonsense_density_estimator_is_refused(self) -> None:
        with pytest.raises(EngineError, match="density_estimator"):
            SBIEngine(bounded_problem(), density_estimator=17)

    def test_zero_draws_are_refused(self) -> None:
        with pytest.raises(EngineError, match="draw"):
            SBIEngine(bounded_problem(), budget=10).run(draws=0)


class TestTheOptionalExtra:
    """``architecture.md`` §4 rules 2 and 3: on use, never on import."""

    def test_constructing_the_engine_needs_no_extra(self) -> None:
        """The construction is pure ampere; only ``run`` reaches for ``sbi``."""
        engine = SBIEngine(bounded_problem(), method="nle", budget=10)
        assert engine.method == "nle"
        assert engine.backend == "reference"

    @pytest.mark.skipif(HAS_SBI, reason="the extra is installed here")
    def test_running_without_the_extra_refuses_by_name(self) -> None:
        with pytest.raises(OptionalDependencyError) as raised:
            SBIEngine(bounded_problem(), budget=10).run(draws=2)
        assert raised.value.extra == "sbi"

    def test_importing_the_namespace_imports_neither_torch_nor_sbi(self) -> None:
        """The rule this module's engine could most easily have broken.

        ``test_engines`` proves no backend is imported; this proves the two
        heavy libraries behind the ``sbi`` extra are not either, which is what
        keeps ``import ampere.inference`` cheap in the base install. Run as a
        subprocess because this interpreter may well have imported both
        already.
        """
        probe = subprocess.run(
            [
                sys.executable,
                "-c",
                (
                    "import sys, ampere.inference; "
                    "print(sorted(n for n in ('torch', 'sbi') if n in sys.modules))"
                ),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        assert probe.stdout.strip() == "[]"


# ---------------------------------------------------------------------------
# 4. NPE against an emcee reference
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def emcee_reference() -> dict[str, tuple[float, float]]:
    """Where "right" is, for a problem with no closed form.

    Long enough that its own Monte Carlo error is small beside the tolerances
    below, and seeded, so the reference is the same number every run.
    """
    run = EmceeEngine(joint_problem(), walkers=16).run(steps=600, burn_in=200)
    posterior = run["posterior"].dataset
    return {
        name: (float(posterior[name].mean()), float(posterior[name].std()))
        for name in ("model.norm", "model.index", "calibration")
    }


@pytest.fixture(scope="module")
def npe_run() -> Any:
    """One NPE fit of the joint problem. The expensive fixture; module-scoped."""
    engine = SBIEngine(joint_problem(), method="npe", budget=2000)
    return engine.run(draws=1000)


@needs_sbi
class TestNPERecoversTheJointPosterior:
    """Accept criterion 1: mean and width against the emcee reference."""

    NAMES = ("model.norm", "model.index", "calibration")

    def test_the_posterior_is_in_the_constrained_space_under_merged_names(
        self, npe_run: Any
    ) -> None:
        posterior = npe_run["posterior"].dataset
        assert set(posterior.data_vars) == set(self.NAMES)
        assert posterior.sizes == {"chain": 1, "draw": 1000}
        # `norm` and `calibration` are lognormal, so the constrained space is
        # positive; a run that forgot to `constrain` would show negatives.
        assert float(np.asarray(posterior["model.norm"]).min()) > 0.0
        assert float(np.asarray(posterior["calibration"]).min()) > 0.0

    @pytest.mark.parametrize("name", NAMES)
    def test_the_mean_matches_the_reference(
        self, name: str, npe_run: Any, emcee_reference: dict[str, tuple[float, float]]
    ) -> None:
        drawn = np.asarray(npe_run["posterior"][name])
        mean, width = emcee_reference[name]
        assert abs(float(drawn.mean()) - mean) < 0.5 * width

    @pytest.mark.parametrize("name", NAMES)
    def test_the_width_matches_the_reference(
        self, name: str, npe_run: Any, emcee_reference: dict[str, tuple[float, float]]
    ) -> None:
        """Width, not only location: an over-confident estimator passes the mean row."""
        drawn = np.asarray(npe_run["posterior"][name])
        _, width = emcee_reference[name]
        assert 0.6 < float(drawn.std()) / width < 1.5


@needs_sbi
class TestTheRunItEmits:
    """``results.md``'s obligations, discharged by this driver like any other."""

    def test_it_emits_the_arviz_groups(self, npe_run: Any) -> None:
        assert sorted(npe_run.children) == [
            "constant_data",
            "log_likelihood",
            "observed_data",
            "posterior",
            "sample_stats",
        ]

    def test_every_draw_carries_the_true_split_scored_on_the_numpy_path(self, npe_run: Any) -> None:
        """The point of scoring through ``Engine.finish``: these are ``evaluate``'s."""
        stats = npe_run["sample_stats"].dataset
        prior = np.asarray(stats["log_prior"])
        likelihood = np.asarray(stats["log_likelihood"])
        assert np.all(np.isfinite(prior))
        assert np.all(np.isfinite(likelihood))
        assert np.asarray(stats["lp"]) == pytest.approx(prior + likelihood)

    def test_the_estimators_own_log_density_is_beside_it_per_draw(self, npe_run: Any) -> None:
        """What makes SBC and importance reweighting possible later."""
        stats = npe_run["sample_stats"].dataset
        estimator = np.asarray(stats["ampere_sbi_log_prob"])
        assert estimator.shape == np.asarray(stats["lp"]).shape
        assert np.all(np.isfinite(estimator))
        # It is the network's density, not ampere's: the two must not be equal.
        assert not np.allclose(estimator, np.asarray(stats["lp"]))
        assert npe_run.attrs["ampere_sbi_log_prob_kind"] == "normalised"

    def test_the_log_likelihood_group_decomposes_per_dataset(self, npe_run: Any) -> None:
        group = npe_run["log_likelihood"].dataset
        assert sorted(group.data_vars) == ["blue", "red"]
        total = sum(np.asarray(group[label]) for label in group.data_vars)
        assert total == pytest.approx(np.asarray(npe_run["sample_stats"]["log_likelihood"]))

    def test_the_attrs_say_what_was_fitted_and_how(self, npe_run: Any) -> None:
        attrs = npe_run.attrs
        assert attrs["ampere_engine"] == "sbi"
        assert attrs["ampere_backend"] == "reference"
        assert attrs["ampere_schema_version"] == 5
        assert attrs["ampere_sbi_method"] == "npe"
        assert attrs["ampere_sbi_trainer"].startswith("NPE")
        assert attrs["ampere_sbi_density_estimator"] == "maf"
        assert attrs["ampere_sbi_budget"] == 2000
        assert attrs["ampere_sbi_rounds"] == 1
        assert attrs["ampere_sbi_simulations"] == 2000
        assert attrs["ampere_sbi_failures"] == 0
        assert attrs["ampere_sbi_embedding"] == "none"
        assert attrs["ampere_sbi_summary_layout"] == "flat"
        assert attrs["ampere_sbi_summary_features"] == 52
        assert attrs["ampere_sbi_parameterisation"] == "unconstrained"
        assert attrs["ampere_sbi_context"] == "none"
        assert attrs["ampere_sbi_device"] == "cpu"
        assert int(attrs["ampere_sbi_epochs_trained"]) > 0
        assert attrs["ampere_sbi_version"]
        assert attrs["ampere_torch_version"]

    def test_the_training_loss_trace_is_recorded_with_its_stride(self, npe_run: Any) -> None:
        import json

        trace = json.loads(npe_run.attrs["ampere_sbi_training_loss"])
        assert len(trace) > 1
        assert int(npe_run.attrs["ampere_sbi_training_loss_stride"]) >= 1
        # A training loss that never fell is a fit that did not happen.
        assert trace[-1] < trace[0]

    def test_it_survives_the_netcdf_round_trip(self, npe_run: Any, tmp_path: Path) -> None:
        """Including the extra sample statistic, which is why it is a variable."""
        from ampere.results import from_netcdf, to_netcdf

        path = to_netcdf(npe_run, tmp_path / "sbi.nc")
        back = from_netcdf(path)
        assert back.attrs["ampere_spec_hash"] == npe_run.attrs["ampere_spec_hash"]
        assert back.attrs["ampere_sbi_method"] == "npe"
        assert np.asarray(back["sample_stats"]["ampere_sbi_log_prob"]) == pytest.approx(
            np.asarray(npe_run["sample_stats"]["ampere_sbi_log_prob"])
        )


# ---------------------------------------------------------------------------
# 5. NLE and NRE: shape and finiteness at CI budgets
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module", params=["nle", "nre"])
def mcmc_run(request: Any) -> tuple[str, Any]:
    """One NLE and one NRE fit, each done once for the rows below.

    Both build an MCMC posterior (``sbi`` has no direct sampler for either),
    so the draw count and the chain settings are what keeps this affordable:
    forty draws off four short chains says everything a shape-and-finiteness
    row can say, and nothing this budget could not support.
    """
    if not HAS_SBI:  # pragma: no cover - the module-level skip covers this
        pytest.skip("needs the 'sbi' extra")
    method = str(request.param)
    engine = SBIEngine(bounded_problem(), method=method, budget=400)
    run = engine.run(
        draws=40,
        training={"max_num_epochs": 30},
        posterior_options={"num_chains": 4, "warmup_steps": 15, "thin": 1},
    )
    return method, run


@needs_sbi
class TestTheOtherTwoFamiliesRunEndToEnd:
    """Accept criterion 2. Both sample by MCMC, so the budgets are tiny."""

    def test_it_produces_a_run_of_the_right_shape(self, mcmc_run: tuple[str, Any]) -> None:
        _, run = mcmc_run
        posterior = run["posterior"].dataset
        assert (posterior.sizes["chain"], posterior.sizes["draw"]) == (1, 40)
        assert set(posterior.data_vars) == {"model.slope"}
        assert np.all(np.isfinite(np.asarray(posterior["model.slope"])))
        # Constrained space: `slope` is lognormal and cannot be negative.
        assert float(np.asarray(posterior["model.slope"]).min()) > 0.0

    def test_the_true_split_is_finite_and_the_attrs_name_the_method(
        self, mcmc_run: tuple[str, Any]
    ) -> None:
        method, run = mcmc_run
        stats = run["sample_stats"].dataset
        assert np.all(np.isfinite(np.asarray(stats["lp"])))
        assert run.attrs["ampere_sbi_method"] == method
        # NLE and NRE know their posterior only up to the evidence, and say so.
        assert run.attrs["ampere_sbi_log_prob_kind"] == "unnormalised"
        assert np.all(np.isfinite(np.asarray(stats["ampere_sbi_log_prob"])))

    def test_the_trainer_recorded_is_the_variant_that_ran(self, mcmc_run: tuple[str, Any]) -> None:
        """``sbi`` aliases ``NLE`` to ``NLE_A`` and ``NRE`` to ``NRE_B``."""
        method, run = mcmc_run
        assert run.attrs["ampere_sbi_trainer"].startswith(method.upper())


# ---------------------------------------------------------------------------
# 6. The embedding vocabulary (legacy parity)
# ---------------------------------------------------------------------------


@needs_sbi
class TestTheEmbeddingVocabulary:
    """``ampere/infer/sbi.py``'s four spellings, carried over and refused precisely."""

    #: A realistic summary width. Not 8: ``sbi``'s ``CNNEmbedding`` defaults to
    #: two convolutions with kernel 5 and a pool of 2, which drives an
    #: eight-feature input to zero width and asserts — the legacy class had the
    #: same trap, and the remedy (fewer layers, or a smaller kernel through the
    #: hyperparameter dict) is in ``sbi``'s own message.
    FEATURES = 52

    @staticmethod
    def resolved(embedding: Any, features: int = FEATURES, free_size: int = 2) -> Any:
        import torch

        from ampere.inference._sbi import _embedding_of

        return _embedding_of(embedding, torch=torch, features=features, free_size=free_size)

    def test_none_is_no_embedding_at_all(self) -> None:
        resolved = self.resolved(None)
        assert (resolved.module, resolved.name, resolved.output_dim) == (None, "none", 0)

    @pytest.mark.parametrize("named", ["FC", "FullyConnected"])
    def test_the_fully_connected_spellings(self, named: str) -> None:
        resolved = self.resolved(named)
        assert resolved.name == "FC"
        assert resolved.output_dim == 4  # 2 * free_size, the legacy default
        assert resolved.module is not None

    @pytest.mark.parametrize("named", [True, "default", "Conv", "CNN"])
    def test_the_convolutional_spellings(self, named: Any) -> None:
        resolved = self.resolved(named)
        assert resolved.name == "CNN"
        assert resolved.output_dim == 4

    def test_a_dict_of_legacy_hyperparameters(self) -> None:
        """The legacy key names, including the one ``sbi`` 0.27 renamed."""
        resolved = self.resolved(
            {
                "type": "CNN",
                "output_dim": 6,
                "n_conv_layers": 1,
                "out_channels_per_layer": [4],
                "kernel_size_per_layer": 3,
                "n_linear_layers": 1,
                "num_linear_units": 8,
            }
        )
        assert (resolved.name, resolved.output_dim) == ("CNN", 6)

    def test_a_dict_for_the_fully_connected_net(self) -> None:
        resolved = self.resolved({"type": "FC", "output_dim": 5, "n_layers": 1, "num_hiddens": 7})
        assert (resolved.name, resolved.output_dim) == ("FC", 5)

    def test_a_user_module_is_taken_and_measured(self) -> None:
        """A user's own net: accepted as-is, and its width measured rather than asked for."""
        import torch

        module = torch.nn.Sequential(torch.nn.Linear(self.FEATURES, 3), torch.nn.ReLU())
        resolved = self.resolved(module)
        assert resolved.module is module
        assert (resolved.name, resolved.output_dim) == ("custom", 3)

    def test_an_unknown_name_is_refused_with_the_vocabulary(self) -> None:
        with pytest.raises(EngineError, match="does not know the embedding"):
            self.resolved("transformer")

    def test_an_unknown_hyperparameter_is_refused_by_name(self) -> None:
        """The legacy code fell through silently here; this does not."""
        with pytest.raises(EngineError, match="hyperparameter"):
            self.resolved({"type": "FC", "num_hidden": 7})

    def test_an_embedding_reaches_the_run_and_is_recorded(self) -> None:
        engine = SBIEngine(bounded_problem(), method="npe", budget=300, embedding="FC")
        run = engine.run(draws=30, training={"max_num_epochs": 10})
        assert run.attrs["ampere_sbi_embedding"] == "FC"
        assert run.attrs["ampere_sbi_embedding_output_dim"] == 2
        assert run["posterior"].dataset.sizes["draw"] == 30


# ---------------------------------------------------------------------------
# 7. The black-box external simulator, including under a process pool
# ---------------------------------------------------------------------------


def external_problem() -> FittingProblem:
    """``examples/sbi/external_simulator.py``'s problem, imported not reimplemented."""
    return _sibling("external_simulator").build_problem(seed=SEED)


@pytest.fixture(scope="module")
def pooled_external_run() -> Any:
    """One fit of the subprocess simulator, its budget run under a process pool.

    Both halves of accept criterion 3 in one run, because each draw of this
    problem is a real ``subprocess.run`` and a second budget would double the
    slowest thing in the module for nothing: the fit exercises the black-box
    path, and ``executor=``/``chunk_size=`` exercise the pool.
    """
    if not HAS_SBI:  # pragma: no cover - the module-level skip covers this
        pytest.skip("needs the 'sbi' extra")
    engine = SBIEngine(
        external_problem(),
        method="npe",
        budget=120,
        executor=ProcessExecutor(2),
        chunk_size=60,
    )
    return engine.run(draws=50, training={"max_num_epochs": 20})


@needs_sbi
class TestABlackBoxSimulatorOnTheReferenceBackend:
    """Accept criterion 3, and the case the whole engine exists for.

    The model is a subprocess call to a program that knows nothing about
    ampere, declares its own failure class, and crashes on part of the prior.
    Nothing about this driver is told any of that.
    """

    def test_it_fits_a_subprocess_simulator(self, pooled_external_run: Any) -> None:
        posterior = pooled_external_run["posterior"].dataset
        assert set(posterior.data_vars) == {"model.index", "model.norm"}
        assert (posterior.sizes["chain"], posterior.sizes["draw"]) == (1, 50)
        assert np.all(np.isfinite(np.asarray(posterior["model.norm"])))
        assert pooled_external_run.attrs["ampere_backend"] == "reference"

    def test_the_budget_went_through_the_process_pool(self, pooled_external_run: Any) -> None:
        """``executor=`` and ``chunk_size=`` straight through to ``simulate_many``."""
        attrs = pooled_external_run.attrs
        assert attrs["ampere_sbi_executor"] == "ProcessExecutor"
        assert attrs["ampere_sbi_chunk_size"] == 60
        assert attrs["ampere_sbi_simulations"] == 120

    def test_the_crashing_draws_are_dropped_and_counted(self, pooled_external_run: Any) -> None:
        """Reject-and-record: the network never sees them, and the run says how many.

        ``toy_powerlaw.py`` exits non-zero for ``norm`` above 9.5 and the prior
        is ``uniform(0.5, 9.5)`` — which ``scipy`` reads as ``[0.5, 10.0]``, so
        about a twentieth of the budget crashes. At 120 draws the chance of
        seeing none is under a per cent, and the run is seeded, so this is a
        fixed number rather than a hopeful one.

        The counts are the *parent's* even though the draws ran in workers,
        which is ``inference.md`` §13's *Execution* paragraph: ``simulate_many``
        owns both ends of the pool and replays the records in draw order.
        """
        attrs = pooled_external_run.attrs
        assert attrs["ampere_sbi_failures"] > 0
        assert (
            attrs["ampere_sbi_usable_simulations"] + attrs["ampere_sbi_failures"]
            == attrs["ampere_sbi_simulations"]
        )
        assert "model_failed" in pooled_external_run.attrs["ampere_failure_counts"]


# ---------------------------------------------------------------------------
# 8. The training-set path and multiple rounds
# ---------------------------------------------------------------------------


@needs_sbi
class TestTheTrainingSetPath:
    """``training_set=``: the pairs on disk, written a chunk at a time."""

    def test_the_pairs_reach_the_file_and_can_be_read_back(self, tmp_path: Path) -> None:
        from ampere.results.training import read_training_set

        path = tmp_path / "budget.nc"
        engine = SBIEngine(
            bounded_problem(), method="npe", budget=60, chunk_size=20, training_set=path
        )
        run = engine.run(draws=20, training={"max_num_epochs": 10})
        assert run.attrs["ampere_sbi_training_set"] == str(path)
        assert path.exists()
        training = read_training_set(path)
        # Every draw is written, failures included; three chunks of twenty.
        assert len(training) == 60


@needs_sbi
class TestTheShippedExample:
    """``examples/sbi/fit_external_simulator.py`` runs, and reports what it fitted.

    An example only the documentation build exercises is an example that rots
    (W2.11's reason for ``tests/examples``). This one is run at the smallest
    budget that still goes through every step it demonstrates.
    """

    def test_it_fits_and_reports(self) -> None:
        example = _sibling("fit_external_simulator")
        run = example.fit(budget=40, workers=2, draws=20, chunk_size=20)
        report = example.report(run)
        assert "npe on reference" in report
        assert "posterior:" in report
        assert "model.norm" in report
        assert "unconstrained coordinates" in report
        assert run["posterior"].dataset.sizes["draw"] == 20


@needs_sbi
class TestMultipleRounds:
    """``rounds=``: each round after the first proposes from the last posterior."""

    def test_two_rounds_run_and_are_recorded(self) -> None:
        engine = SBIEngine(bounded_problem(), method="npe", budget=200, rounds=2)
        run = engine.run(draws=40, training={"max_num_epochs": 15})
        assert run.attrs["ampere_sbi_rounds"] == 2
        assert run.attrs["ampere_sbi_simulations"] == 400
        assert run["posterior"].dataset.sizes["draw"] == 40
        assert np.all(np.isfinite(np.asarray(run["sample_stats"]["lp"])))
