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
import itertools
import json
import math
import subprocess
import sys
import warnings
from pathlib import Path
from typing import Any, ClassVar

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.backends.reference import Resample
from ampere.core import (
    Dataset,
    EncodingError,
    EncodingLayout,
    FittingProblem,
    Instrument,
    Model,
    Parameter,
    ProcessExecutor,
    ScaledSigma,
    Spectrum,
    VisibilitySet,
    encode_observations,
)
from ampere.core.exceptions import OptionalDependencyError
from ampere.inference import DEFAULT_TRUNCATION_EPSILON, EmceeEngine, EngineError, SBIEngine
from ampere.inference._sbi import METHODS, SUMMARY_LAYOUT, _summary_of, _thinned
from ampere.inference._tmnre import (
    GRID_POINTS_1D,
    GRID_POINTS_2D,
    TruncationBox,
    grid_between,
    interval_above,
    marginal_indices,
    marginal_log_density,
    pair_labels,
    pair_mesh,
)
from ampere.results import PROVENANCE_SCHEMA_VERSION
from ampere.results.artefacts import ArtefactCacheWarning


class _Observed:
    """A stand-in dataset: the encoding reads ``observed`` and nothing else.

    Used only where composing a whole problem would mean choosing a likelihood
    family the test has no opinion about -- a complex ``VisibilitySet``, say.
    """

    def __init__(self, observed: Any) -> None:
        self.observed = observed
        self.effective_mask = None if observed.mask is None else np.asarray(observed.mask).ravel()


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
        """The ``"set"`` layout encodes real and imaginary parts; ``"flat"`` will not.

        The refusal moved with W3.3 from this driver into the layout, which is
        where it belongs: a complex container is a fact about the *problem*, so
        it can be refused before any observation is packed rather than at the
        first draw.
        """
        visibilities = VisibilitySet(
            [10.0, -30.0] * u.dimensionless_unscaled,
            [-20.0, 40.0] * u.dimensionless_unscaled,
            [2.2, 2.2] * u.um,
            np.array([1.0 + 2.0j, 3.0 + 0.0j]) * u.Jy,
        )
        datasets = {"vis": _Observed(visibilities)}
        with pytest.raises(EncodingError, match="complex"):
            EncodingLayout.from_datasets(datasets, kind="flat")
        # The set layout takes it: two value columns, real and imaginary.
        assert EncodingLayout.from_datasets(datasets).group("value").width == 2

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

    def test_a_context_that_is_not_a_prior_is_refused_by_name(self) -> None:
        """**W5.10** filled the reserved slot; it did not widen it to anything.

        The row this replaces asserted that every value but ``None`` was
        refused, which was right while the slot was reserved. What it takes
        now is a ``ContextPrior``, and a bare settings mapping -- the shape a
        user would most plausibly reach for -- is still refused, by a message
        that says what one is.
        """
        with pytest.raises(EngineError, match="ContextPrior"):
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

    @pytest.mark.skipif(HAS_SBI, reason="the extra is installed here")
    def test_calibrating_without_the_extra_refuses_by_name(self) -> None:
        """W3.6's fast path is behind the same door as the fit itself.

        ``sbi.diagnostics`` is what does the arithmetic, so ``calibrate`` must
        reach for the extra exactly as ``run`` does — and must say so with the
        extra named, rather than failing on ``self.posterior is None`` and
        blaming the user for not having run a fit that could never have run.
        """
        with pytest.raises(OptionalDependencyError) as raised:
            SBIEngine(bounded_problem(), budget=10).calibrate(count=4, posterior_draws=4)
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

    def test_the_contract_name_is_beside_it_too(self, npe_run: Any) -> None:
        """W5.0: the engine-neutral name, not a bare alias of the SBI-specific one.

        Both are the estimator's own density at the same draw, but in
        different coordinates (module docstring): ``ampere_sbi_log_prob`` is
        unconstrained, ``proposal_log_density`` is moved into the same
        constrained coordinates the stored ``log_prior``/``log_likelihood``
        are, so the two are finite together and need not be numerically
        equal.
        """
        stats = npe_run["sample_stats"].dataset
        contract = np.asarray(stats["proposal_log_density"])
        estimator = np.asarray(stats["ampere_sbi_log_prob"])
        assert contract.shape == estimator.shape
        assert np.all(np.isfinite(contract))

    def test_it_writes_the_engine_neutral_approximation_family(self, npe_run: Any) -> None:
        assert npe_run.attrs["ampere_approximation"] == "density_estimator"

    def test_the_log_likelihood_group_decomposes_per_dataset(self, npe_run: Any) -> None:
        group = npe_run["log_likelihood"].dataset
        assert sorted(group.data_vars) == ["blue", "red"]
        total = sum(np.asarray(group[label]) for label in group.data_vars)
        assert total == pytest.approx(np.asarray(npe_run["sample_stats"]["log_likelihood"]))

    def test_the_attrs_say_what_was_fitted_and_how(self, npe_run: Any) -> None:
        attrs = npe_run.attrs
        assert attrs["ampere_engine"] == "sbi"
        assert attrs["ampere_backend"] == "reference"
        assert attrs["ampere_schema_version"] == PROVENANCE_SCHEMA_VERSION
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
    def resolved(
        embedding: Any,
        features: int = FEATURES,
        free_size: int = 2,
        layout: Any = None,
    ) -> Any:
        import torch

        from ampere.inference._sbi import _embedding_of

        return _embedding_of(
            embedding, torch=torch, features=features, free_size=free_size, layout=layout
        )

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
            self.resolved("convnext")

    def test_the_set_and_transformer_default_width_is_raised_but_flat_keeps_its_own(
        self,
    ) -> None:
        """W3.11's two defaults: ``max(2 * free_size, 32)`` for a pooled set,
        ``2 * free_size`` (legacy parity, unchanged) for ``"flat"``.
        """
        set_layout = EncodingLayout.from_datasets(two_dataset_problem().datasets, kind="set")

        # Accept criterion: a one-parameter problem's set embedding gets 32,
        # not 2 -- the legacy default would leave an untrained ReLU net a
        # coin flip away from emitting all zeros.
        for named in ("set", "transformer"):
            resolved = self.resolved(named, free_size=1, layout=set_layout)
            assert resolved.output_dim == 32

        # Above the floor, 2 * free_size wins for a set embedding too.
        resolved = self.resolved("set", free_size=20, layout=set_layout)
        assert resolved.output_dim == 40

        # "flat" (CNN/FC) is untouched: 2 * free_size, even below 32.
        resolved = self.resolved("FC", free_size=1, layout=None)
        assert resolved.output_dim == 2

        # An explicit output_dim= always wins over either default.
        resolved = self.resolved({"type": "set", "output_dim": 5}, free_size=1, layout=set_layout)
        assert resolved.output_dim == 5

    def test_a_set_embedding_under_a_flat_layout_is_refused_by_name(self) -> None:
        """W3.3's two need column groups, and the flat layout has none."""
        flat = EncodingLayout.from_datasets(two_dataset_problem().datasets, kind="flat")
        for named in ("set", "transformer"):
            with pytest.raises(EngineError, match="column groups"):
                self.resolved(named, layout=flat)

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


# ---------------------------------------------------------------------------
# 9. W3.3: the coordinate-value-mask layout and the two set embeddings
# ---------------------------------------------------------------------------


class TestTheLayoutArgument:
    """``layout=``: what the network is trained on, named and hashed. No extra."""

    def test_the_default_is_the_flat_summary_w3_2_shipped(self) -> None:
        engine = SBIEngine(two_dataset_problem(), budget=10)
        assert engine.layout == SUMMARY_LAYOUT == "flat"
        assert engine.encoding is None  # resolved at run(), not at construction

    def test_a_set_layout_packs_one_row_per_sample_of_every_dataset(self) -> None:
        problem = two_dataset_problem()
        layout = EncodingLayout.from_datasets(problem.datasets, kind="set")
        # 'short' has 2 samples, 'long' 4 -- masked ones are rows, not absences.
        assert layout.row_cap == 6
        assert layout.bounds == ((0, 2), (2, 6))
        encoded = encode_observations(problem.datasets, layout=layout)
        mask = np.asarray(encoded.values)[0, :, layout.group("mask").offset]
        assert mask.tolist() == [1.0, 1.0, 1.0, 0.0, 1.0, 1.0]

    def test_an_unknown_layout_name_is_refused_with_the_two_that_exist(self) -> None:
        from ampere.inference._sbi import _layout_of

        with pytest.raises(EngineError, match="does not know the layout"):
            _layout_of(two_dataset_problem(), "sequence")

    def test_a_layout_from_another_problem_is_refused_saying_which_field_differs(
        self,
    ) -> None:
        """Accept criterion: a mismatch is refused **by name**, field by field."""
        from ampere.inference._sbi import _layout_of

        other = EncodingLayout.from_datasets(bounded_problem().datasets)
        with pytest.raises(EngineError) as raised:
            _layout_of(two_dataset_problem(), other)
        assert "dataset labels" in str(raised.value)
        assert other.hash in str(raised.value)

    def test_a_row_cap_the_observation_exceeds_is_refused_by_name(self) -> None:
        with pytest.raises(EncodingError, match="row cap"):
            EncodingLayout.from_datasets(two_dataset_problem().datasets, row_cap=2)

    def test_something_that_is_neither_a_name_nor_a_layout_is_refused(self) -> None:
        from ampere.inference._sbi import _layout_of

        with pytest.raises(EngineError, match="layout="):
            _layout_of(two_dataset_problem(), 17)

    def test_a_set_layout_with_no_embedding_is_refused_by_name(self) -> None:
        """A flow over a padded matrix is not a posterior over an observation.

        Checked through ``_embedding_of`` rather than through ``run``, so that it
        runs in ``dev`` like every other refusal in this section: the refusal is
        reached before the ``torch`` argument is ever touched, which is why
        ``None`` is a legitimate thing to pass here and an
        ``OptionalDependencyError`` is not what a caller should meet first.
        """
        from ampere.inference._sbi import _embedding_of

        layout = EncodingLayout.from_datasets(two_dataset_problem().datasets, kind="set")
        with pytest.raises(EngineError, match="no embedding"):
            _embedding_of(None, torch=None, features=0, free_size=1, layout=layout)


@needs_sbi
class TestTheSetEmbeddingWrapper:
    """The masked-pooling wrapper, checked on its own before any training.

    Three properties, and each is a thing a network would otherwise have to
    *learn* -- badly, and only for the layout it saw: padded rows contribute
    nothing, masked rows contribute nothing, and row order means nothing.
    """

    #: A wide output and a fixed seed, so the untrained net is *alive*: the two
    #: shipped nets both end in a ReLU, and a narrow randomly initialised one
    #: emits all zeros, under which every invariance below would hold of
    #: nothing at all. ``test_a_real_row_does_reach_it`` is the row that would
    #: catch that, and it is why it is here.
    WIDTH = 16

    @staticmethod
    def built(problem: Any, *, row_cap: int | None = None) -> Any:
        import torch

        from ampere.inference._sbi import _embedding_of

        torch.manual_seed(20260909)
        layout = EncodingLayout.from_datasets(problem.datasets, row_cap=row_cap)
        resolved = _embedding_of(
            {"type": "set", "output_dim": TestTheSetEmbeddingWrapper.WIDTH},
            torch=torch,
            features=0,
            free_size=problem.free_size,
            layout=layout,
        )
        module = resolved.module
        module.eval()
        return layout, module

    @staticmethod
    def tensor(problem: Any, layout: Any) -> Any:
        import torch

        return torch.as_tensor(
            np.asarray(encode_observations(problem.datasets, layout=layout).values),
            dtype=torch.float32,
        )

    def test_padded_rows_contribute_nothing_whatever_they_hold(self) -> None:
        """Accept criterion: adding padded rows changes nothing."""
        import torch

        problem = two_dataset_problem()
        layout, module = self.built(problem, row_cap=11)
        x = self.tensor(problem, layout)
        with torch.no_grad():
            plain = module(x)
        junk = x.clone()
        junk[:, 6:, :] = 3.7  # every padded row, every column, including the mask
        junk[:, 6:, layout.group("mask").offset] = 0.0  # ... except the mask
        with torch.no_grad():
            noisy = module(junk)
        assert torch.allclose(plain, noisy, atol=1e-6)

    def test_a_masked_rows_value_does_not_reach_the_embedding(self) -> None:
        import torch

        problem = two_dataset_problem()
        layout, module = self.built(problem)
        x = self.tensor(problem, layout)
        with torch.no_grad():
            plain = module(x)
        # Row 3 is 'long' sample 1, the masked one.
        moved = x.clone()
        moved[:, 3, layout.group("value").offset] = 99.0
        moved[:, 3, layout.group("value_asinh").offset] = -99.0
        with torch.no_grad():
            changed = module(moved)
        assert torch.allclose(plain, changed, atol=1e-6)

    def test_reversing_the_rows_changes_nothing(self) -> None:
        """Permutation invariance: the property the whole packing is for."""
        import torch

        problem = two_dataset_problem()
        layout, module = self.built(problem)
        x = self.tensor(problem, layout)
        with torch.no_grad():
            plain = module(x)
            reversed_rows = module(torch.flip(x, dims=(1,)))
        assert torch.allclose(plain, reversed_rows, atol=1e-6)

    def test_a_real_row_does_reach_it(self) -> None:
        """The three invariances above would all hold of a constant function."""
        import torch

        problem = two_dataset_problem()
        layout, module = self.built(problem)
        x = self.tensor(problem, layout)
        moved = x.clone()
        moved[:, 0, layout.group("value").offset] = 42.0
        with torch.no_grad():
            assert not torch.allclose(module(x), module(moved), atol=1e-6)

    def test_it_pools_with_a_mean_rather_than_a_sum(self) -> None:
        """Trap 5: a summed embedding scales with the row count."""
        _, module = self.built(two_dataset_problem())
        assert module.net.aggregation_fn == "mean"


@needs_sbi
class TestTheTransformerWrapper:
    """Non-causal, no positional embedding, explicit dropout, and masked."""

    @staticmethod
    def built(problem: Any) -> Any:
        return TestTheTransformerWrapper.built_with_row_cap(problem, row_cap=None)

    @staticmethod
    def built_with_row_cap(problem: Any, *, row_cap: int | None) -> Any:
        import torch

        from ampere.inference._sbi import _embedding_of

        torch.manual_seed(20260909)
        layout = EncodingLayout.from_datasets(problem.datasets, row_cap=row_cap)
        resolved = _embedding_of(
            {"type": "transformer", "output_dim": 16},
            torch=torch,
            features=0,
            free_size=problem.free_size,
            layout=layout,
        )
        module = resolved.module
        module.eval()
        return layout, module

    def test_the_defaults_sbi_ships_are_all_overridden(self) -> None:
        """``is_causal``, ``pos_emb`` and both dropouts (trap 4)."""
        _, module = self.built(two_dataset_problem())
        config = module.net.config
        assert module.net.is_causal is False
        assert config["is_causal"] is False
        assert config["pos_emb"] == "none"
        assert config["attention_dropout"] == 0.1
        assert config["vit_dropout"] == 0.1

    def test_a_masked_rows_value_does_not_change_the_output(self) -> None:
        """It cannot: the wrapper zeroes a masked row's token after projection."""
        import torch

        problem = two_dataset_problem()
        layout, module = self.built(problem)
        x = torch.as_tensor(
            np.asarray(encode_observations(problem.datasets, layout=layout).values),
            dtype=torch.float32,
        )
        moved = x.clone()
        moved[:, 3, layout.group("value").offset] = 99.0
        moved[:, 3, layout.group("value_asinh").offset] = -99.0
        with torch.no_grad():
            plain = module(x)
            assert torch.allclose(plain, module(moved), atol=1e-6)
        # The positive control: a *retained* row does change it, so the row
        # above is excluded rather than the whole net being deaf.
        retained = x.clone()
        retained[:, 4, layout.group("value").offset] = 99.0
        with torch.no_grad():
            assert not torch.allclose(plain, module(retained), atol=1e-6)

    def test_padded_rows_contribute_nothing_whatever_they_hold(self) -> None:
        """Accept criterion: padded rows do not change the pooled output."""
        import torch

        problem = two_dataset_problem()
        layout, module = self.built_with_row_cap(problem, row_cap=11)
        x = torch.as_tensor(
            np.asarray(encode_observations(problem.datasets, layout=layout).values),
            dtype=torch.float32,
        )
        with torch.no_grad():
            plain = module(x)
        junk = x.clone()
        junk[:, 6:, :] = 3.7  # every padded row, every column, including the mask
        junk[:, 6:, layout.group("mask").offset] = 0.0  # ... except the mask
        with torch.no_grad():
            noisy = module(junk)
        assert torch.allclose(plain, noisy, atol=1e-6)

    def test_reversing_the_rows_changes_nothing(self) -> None:
        """W3.11's accept criterion: this is the assertion the last-token read fails.

        Built and checked exactly as the set embedding's own permutation test,
        against the *same* module twice (plain and row-reversed) so that a
        regression to the last-token read -- which reads a different row after
        a reversal, and so a different summary -- would be caught here rather
        than only in a slower end-to-end fit.
        """
        import torch

        problem = two_dataset_problem()
        layout, module = self.built(problem)
        x = torch.as_tensor(
            np.asarray(encode_observations(problem.datasets, layout=layout).values),
            dtype=torch.float32,
        )
        with torch.no_grad():
            plain = module(x)
            reversed_rows = module(torch.flip(x, dims=(1,)))
        assert torch.allclose(plain, reversed_rows, atol=1e-6)

    def test_the_readout_depends_on_named_sbi_attributes(self) -> None:
        """The masked-mean readout has no hook onto sbi's forward; it reads
        the net's body modules by name, and must fail loudly, not silently
        wrongly, if a future sbi renames one.
        """
        import torch

        from ampere.inference._sbi import _transformer_masked_mean

        problem = two_dataset_problem()
        layout, module = self.built(problem)
        tokens = torch.zeros(1, layout.row_cap, module.net.aggregator.in_features)
        keep = torch.ones(1, layout.row_cap, dtype=torch.bool)

        # The real net has every attribute this depends on, and is causal=False.
        _transformer_masked_mean(module.net, tokens, keep)  # does not raise

        class Renamed:
            """A stand-in missing one of the attributes the readout reads."""

            preprocess = staticmethod(lambda x: x)
            norm = staticmethod(lambda x: x)
            aggregator = staticmethod(lambda x: x)
            is_causal = False
            # 'layers' is deliberately absent.

        with pytest.raises(AttributeError, match=r"\['layers'\]"):
            _transformer_masked_mean(Renamed(), tokens, keep)

        class Causal:
            """A stand-in with every attribute, but the wrong value for one."""

            preprocess = staticmethod(lambda x: x)
            layers: tuple[Any, ...] = ()
            norm = staticmethod(lambda x: x)
            aggregator = staticmethod(lambda x: x)
            is_causal = True

        with pytest.raises(AttributeError, match="is_causal=True"):
            _transformer_masked_mean(Causal(), tokens, keep)


@needs_sbi
class TestTheSetLayoutEndToEnd:
    """Accept criterion: it trains and samples on a two-dataset toy problem."""

    def test_the_set_embedding_trains_and_samples(self) -> None:
        engine = SBIEngine(
            two_dataset_problem(), method="npe", budget=200, embedding="set", layout="set"
        )
        run = engine.run(draws=25, training={"max_num_epochs": 8})
        posterior = run["posterior"].dataset
        assert (posterior.sizes["chain"], posterior.sizes["draw"]) == (1, 25)
        assert set(posterior.data_vars) == {"model.slope"}
        assert float(np.asarray(posterior["model.slope"]).min()) > 0.0
        assert np.all(np.isfinite(np.asarray(run["sample_stats"]["lp"])))
        assert run.attrs["ampere_sbi_embedding"] == "set"

    def test_the_transformer_embedding_trains_at_a_tiny_budget(self) -> None:
        engine = SBIEngine(
            two_dataset_problem(),
            method="npe",
            budget=80,
            embedding="transformer",
            layout="set",
        )
        run = engine.run(draws=10, training={"max_num_epochs": 3})
        assert run["posterior"].dataset.sizes["draw"] == 10
        assert run.attrs["ampere_sbi_embedding"] == "transformer"

    def test_the_run_records_the_layout_by_name_and_by_hash(self) -> None:
        """Accept criterion: the name and the hash are in the run's attrs."""
        import json

        problem = two_dataset_problem()
        layout = EncodingLayout.from_datasets(problem.datasets)
        engine = SBIEngine(problem, method="npe", budget=80, embedding="set", layout=layout)
        run = engine.run(draws=10, training={"max_num_epochs": 3})
        assert run.attrs["ampere_sbi_summary_layout"] == "set"
        assert run.attrs["ampere_encoding_hash"] == layout.hash
        assert run.attrs["ampere_sbi_encoding_rows"] == layout.row_cap
        assert run.attrs["ampere_sbi_encoding_columns"] == layout.columns_total
        recorded = json.loads(run.attrs["ampere_encoding_layout"])
        assert EncodingLayout.from_dict(recorded) == layout
        assert engine.encoding is layout

    def test_a_training_set_records_the_layout_too(self, tmp_path: Path) -> None:
        """Accept criterion: and in a training set's attrs, so W3.5 can key on it."""
        from ampere.results.training import read_training_set

        path = tmp_path / "set_budget.nc"
        problem = two_dataset_problem()
        engine = SBIEngine(
            problem,
            method="npe",
            budget=40,
            chunk_size=20,
            embedding="set",
            layout="set",
            training_set=path,
        )
        run = engine.run(draws=8, training={"max_num_epochs": 3})
        training = read_training_set(path)
        assert len(training) == 40
        import xarray

        stored = xarray.open_datatree(path)
        try:
            assert stored.attrs["ampere_encoding_hash"] == run.attrs["ampere_encoding_hash"]
            assert stored.attrs["ampere_encoding_layout"]
        finally:
            stored.close()

    def test_sbi_is_told_not_to_z_score_a_set_tensor(self) -> None:
        """Trap 1: sbi's column-wise z-scoring would standardise the mask column."""
        import sbi.neural_nets

        seen: dict[str, Any] = {}
        original = sbi.neural_nets.posterior_nn

        def spy(*args: Any, **kwargs: Any) -> Any:
            seen.update(kwargs)
            return original(*args, **kwargs)

        problem = two_dataset_problem()
        engine = SBIEngine(problem, budget=20, embedding="set", layout="set")
        sbi.neural_nets.posterior_nn = spy
        try:
            engine.run(draws=2, training={"max_num_epochs": 1})
        finally:
            sbi.neural_nets.posterior_nn = original
        assert seen["z_score_x"] == "none"

    def test_the_flat_layout_keeps_sbis_own_z_scoring(self) -> None:
        """Where it is right: one row of real values, and nothing structural."""
        import sbi.neural_nets

        seen: dict[str, Any] = {}
        original = sbi.neural_nets.posterior_nn

        def spy(*args: Any, **kwargs: Any) -> Any:
            seen.update(kwargs)
            return original(*args, **kwargs)

        engine = SBIEngine(bounded_problem(), budget=20, embedding="FC")
        sbi.neural_nets.posterior_nn = spy
        try:
            engine.run(draws=2, training={"max_num_epochs": 1})
        finally:
            sbi.neural_nets.posterior_nn = original
        assert "z_score_x" not in seen


@needs_sbi
class TestTheArtefactCache:
    """W3.5's wiring: a second identical run trains nothing and says so."""

    def test_a_second_identical_run_is_a_hit_that_trains_nothing(self, tmp_path: Any) -> None:
        from ampere.results.artefacts import ArtefactStore

        store = ArtefactStore(tmp_path / "artefacts")
        first = SBIEngine(joint_problem(), method="npe", budget=150, cache=store).run(
            draws=20, training={"max_num_epochs": 3}
        )
        assert first.attrs["ampere_sbi_cache_hit"] == 0
        assert first.attrs["ampere_sbi_simulations"] > 0

        engine = SBIEngine(joint_problem(), method="npe", budget=150, cache=store)
        second = engine.run(draws=20, training={"max_num_epochs": 3})
        assert second.attrs["ampere_sbi_cache_hit"] == 1
        assert second.attrs["ampere_sbi_cache_key"] == first.attrs["ampere_sbi_cache_key"]
        # Nothing was simulated or trained on the hit, and the run says so
        # rather than repeating the first run's numbers.
        assert second.attrs["ampere_sbi_simulations"] == 0
        assert engine.posterior is not None

    def test_a_different_budget_is_a_miss(self, tmp_path: Any) -> None:
        from ampere.results.artefacts import ArtefactStore

        store = ArtefactStore(tmp_path / "artefacts")
        SBIEngine(joint_problem(), method="npe", budget=150, cache=store).run(
            draws=10, training={"max_num_epochs": 2}
        )
        run = SBIEngine(joint_problem(), method="npe", budget=160, cache=store).run(
            draws=10, training={"max_num_epochs": 2}
        )
        assert run.attrs["ampere_sbi_cache_hit"] == 0

    def test_without_a_store_the_attrs_say_nothing_about_a_cache(self) -> None:
        run = SBIEngine(joint_problem(), method="npe", budget=120).run(
            draws=10, training={"max_num_epochs": 2}
        )
        assert "ampere_sbi_cache_hit" not in run.attrs

    def test_a_set_embedding_posterior_survives_the_store_too(self, tmp_path: Any) -> None:
        from ampere.results.artefacts import ArtefactStore

        store = ArtefactStore(tmp_path / "artefacts")
        settings = dict(method="npe", budget=120, embedding="set", layout="set", cache=store)
        SBIEngine(two_dataset_problem(), **settings).run(draws=10, training={"max_num_epochs": 2})
        run = SBIEngine(two_dataset_problem(), **settings).run(
            draws=10, training={"max_num_epochs": 2}
        )
        assert run.attrs["ampere_sbi_cache_hit"] == 1


@needs_sbi
class TestServingANamedArtefact:
    """W5.23: ``serve_artefact=`` restores one named digest regardless of match.

    ``joint_problem()``'s ``budget=`` is the single key ingredient moved here
    (as ``TestTheArtefactCache.test_a_different_budget_is_a_miss`` above
    already does) rather than a ``HilbertSpaceGP`` ``basis_size``: this
    module's fixtures have no HSGP problem, and budget is already the
    established "one differing ingredient" case in this file.
    """

    def test_a_served_artefact_from_a_different_budget_is_restored_and_sampled(
        self, tmp_path: Path
    ) -> None:
        from ampere.results.artefacts import ArtefactStore

        store = ArtefactStore(tmp_path / "artefacts")
        trained = SBIEngine(joint_problem(), method="npe", budget=150, cache=store).run(
            draws=20, training={"max_num_epochs": 3}
        )
        served_digest = trained.attrs["ampere_sbi_cache_key"]

        engine = SBIEngine(
            joint_problem(),
            method="npe",
            budget=160,
            cache=store,
            serve_artefact=served_digest,
        )
        with pytest.warns(ArtefactCacheWarning, match="budget"):
            served_run = engine.run(draws=20, training={"max_num_epochs": 3})

        # Restored, not trained: nothing simulated, and the posterior is set.
        assert engine.posterior is not None
        assert served_run.attrs["ampere_sbi_simulations"] == 0
        assert served_run.attrs["ampere_sbi_artefact_served"] == served_digest
        assert served_run.attrs["ampere_sbi_cache_key"] != served_digest
        mismatch = json.loads(served_run.attrs["ampere_sbi_artefact_mismatch"])
        assert mismatch == {"budget": [150, 160]}

        # Calibration works on the served posterior exactly as on a trained one.
        # ``calibrate()`` used to sample the posterior for its SBC ranks
        # through ``sbi.diagnostics.run_sbc`` without going through
        # ``self._seeded`` (unlike ``run()``, which reseeds torch's *global*
        # generator before every draw and restores it afterwards), leaving
        # torch's global RNG advanced by however many draws SBC took --
        # caught here first, by ``TestTheCalibrationFastPath``'s TARP
        # thresholds intermittently moving when this test ran earlier in the
        # same session. W5.28(j) closed the gap: ``calibrate()`` now saves
        # and restores torch's (and numpy's legacy global) state itself, so
        # this test no longer has to do it by hand to protect its neighbours.
        report = engine.calibrate(count=8, posterior_draws=8, tarp=False)
        assert report is not None

    def test_a_matching_digest_serves_with_an_empty_mismatch_and_no_warning(
        self, tmp_path: Path
    ) -> None:
        from ampere.results.artefacts import ArtefactStore

        store = ArtefactStore(tmp_path / "artefacts")
        trained = SBIEngine(joint_problem(), method="npe", budget=150, cache=store).run(
            draws=10, training={"max_num_epochs": 2}
        )
        served_digest = trained.attrs["ampere_sbi_cache_key"]

        engine = SBIEngine(
            joint_problem(), method="npe", budget=150, cache=store, serve_artefact=served_digest
        )
        with warnings.catch_warnings():
            warnings.simplefilter("error", ArtefactCacheWarning)
            served_run = engine.run(draws=10, training={"max_num_epochs": 2})
        assert served_run.attrs["ampere_sbi_artefact_mismatch"] == "{}"

    def test_a_wrong_digest_is_refused_by_name(self, tmp_path: Path) -> None:
        from ampere.results.artefacts import ArtefactStore

        store = ArtefactStore(tmp_path / "artefacts")
        SBIEngine(joint_problem(), method="npe", budget=150, cache=store).run(
            draws=10, training={"max_num_epochs": 2}
        )
        engine = SBIEngine(
            joint_problem(),
            method="npe",
            budget=150,
            cache=store,
            serve_artefact="0" * 32,
        )
        with pytest.raises(EngineError, match="serve_artefact"):
            engine.run(draws=10, training={"max_num_epochs": 2})

    def test_serve_artefact_without_a_cache_is_refused_by_name(self) -> None:
        with pytest.raises(EngineError, match="cache"):
            SBIEngine(joint_problem(), method="npe", budget=150, serve_artefact="0" * 32)


# ---------------------------------------------------------------------------
# 12. Family D's fast path: SBC and TARP on the trained posterior (W3.6)
# ---------------------------------------------------------------------------


class _Narrowed:
    """A posterior wrapper whose draws are squeezed towards their own mean.

    The deliberately miscalibrated arm, and the reason it has to be a *wrapper*
    rather than a badly trained network: a network trained on too little data is
    miscalibrated by an amount nobody controls, so a test built on one would be
    a test of the training budget. Temperature-scaling a **correct** posterior
    is miscalibration of a known size and a known kind — over-confident, unbiased
    — which is exactly what SBC's U-shaped rank histogram and TARP's negative
    area-to-curve are defined to detect.

    Only the two methods ``sbi.diagnostics`` reaches for are implemented, which
    is the whole surface ``run_sbc``/``run_tarp`` use: ``sample_batched`` on the
    batched path and ``sample`` on the fallback.
    """

    def __init__(self, inner: Any, temperature: float = 0.3) -> None:
        self.inner = inner
        self.temperature = temperature

    def _squeeze(self, drawn: Any) -> Any:
        centre = drawn.mean(dim=0, keepdim=True)
        return centre + (drawn - centre) * self.temperature

    def sample_batched(self, sample_shape: Any, x: Any = None, **kwargs: Any) -> Any:
        return self._squeeze(self.inner.sample_batched(sample_shape, x=x, **kwargs))

    def sample(self, sample_shape: Any, x: Any = None, **kwargs: Any) -> Any:
        return self._squeeze(self.inner.sample(sample_shape, x=x, **kwargs))


#: The calibration smoke budget. A hundred simulations is ``sbi``'s own floor
#: for ``check_sbc`` and ``run_tarp`` (both warn below it), and a hundred
#: posterior draws is enough resolution for a rank on a one-parameter problem.
#: Re-conditioning an amortised posterior is free, so the whole check costs one
#: simulation batch — about a second beside the fit's own thirteen.
CALIBRATION_COUNT = 100
CALIBRATION_DRAWS = 100


@pytest.fixture(scope="module")
def calibration_engine() -> Any:
    """One trained NPE engine on the one-parameter bounded problem.

    ``bounded_problem`` rather than the joint one: a single lognormal parameter
    over four points is a problem NPE learns *well* at a budget a per-PR gate
    can afford, and family D's null needs a posterior that really is calibrated
    — an under-trained network would fail the check for a reason that has
    nothing to do with the code under test.
    """
    engine = SBIEngine(bounded_problem(), method="npe", budget=1000)
    engine.run(200, training={"max_num_epochs": 200})
    return engine


@pytest.fixture(scope="module")
def calibration(calibration_engine: Any) -> Any:
    return calibration_engine.calibrate(count=CALIBRATION_COUNT, posterior_draws=CALIBRATION_DRAWS)


@needs_sbi
class TestTheCalibrationFastPath:
    """``diagnostics.md`` §11 on the route where re-conditioning is free."""

    def test_the_group_has_the_shape_results_md_promises(self, calibration: Any) -> None:
        assert calibration.sizes["simulation"] == CALIBRATION_COUNT
        assert [str(name) for name in calibration.coords["parameter"].values] == ["model.slope"]
        assert list(calibration["ranks"].dims) == ["simulation", "parameter"]
        assert list(calibration["coverage"].dims) == ["level", "parameter"]
        assert "c2st_ranks" in calibration
        assert "tarp_coverage" in calibration
        ranks = np.asarray(calibration["ranks"].values)
        assert ranks.min() >= 0 and ranks.max() <= CALIBRATION_DRAWS

    def test_it_records_the_layout_the_batch_was_encoded_under(
        self, calibration: Any, calibration_engine: Any
    ) -> None:
        """The one thing that could go silently wrong, asserted.

        A calibration batch packed under a layout rebuilt from the *simulated*
        containers would standardise its columns differently from the ones the
        network trained on, and would then report a different network as
        calibrated. The recorded hash is what makes that checkable after the
        fact, and it must be the run's own.
        """
        assert calibration.attrs["ampere_calibration_route"] == "sbi"
        assert calibration.attrs["ampere_calibration_method"] == "npe"
        assert calibration.attrs["ampere_calibration_parameterisation"] == "unconstrained"
        assert calibration.attrs["ampere_calibration_encoding_hash"] == (
            calibration_engine.encoding.hash
        )
        assert calibration.attrs["ampere_calibration_uniformity_check"] == (
            "sbi.diagnostics.check_sbc"
        )

    def test_a_trained_npe_posterior_passes_check_sbc(self, calibration: Any) -> None:
        """Accept criterion: uniform within ``check_sbc``'s own thresholds.

        W5.26 (7): the threshold was ``> 0.05`` until a GitHub-hosted runner
        (CI run 35684651060, reproduced byte-for-byte on 35688582595 at
        0.01983926) failed it on a genuinely calibrated posterior — a p-value
        threshold set at 0.05 fails 5 % of calibrated fits *by construction*,
        so any new machine's floats re-roll that die. ``> 0.001`` keeps two
        orders of magnitude clear of that false-positive rate while staying
        far below the miscalibrated arm's own ``< 0.01`` a few lines down
        (and the ``> 0.01`` precedent at
        ``TestAmortisationOverTheObservationContext.test_it_stays_calibrated_
        at_a_rescale_the_prior_covers`` above), so the two arms remain
        separated by more than an order of magnitude either way.
        """
        assert float(np.min(calibration["ks_pvalue"].values)) > 0.001
        # C2ST between the ranks and a uniform baseline: 0.5 is "indistinguishable".
        assert float(np.max(calibration["c2st_ranks"].values)) < 0.65

    def test_tarps_expected_coverage_passes_check_tarp(self, calibration: Any) -> None:
        """The joint diagnostic, which the marginal ranks cannot stand in for."""
        assert abs(float(calibration.attrs["ampere_calibration_tarp_atc"])) < 0.03
        assert float(calibration.attrs["ampere_calibration_tarp_ks_pvalue"]) > 0.05

    def test_a_temperature_scaled_posterior_fails_the_same_check(
        self, calibration_engine: Any, calibration: Any
    ) -> None:
        """The arm that proves the check can fail — same engine, same batch size."""
        narrowed = calibration_engine.calibrate(
            count=CALIBRATION_COUNT,
            posterior_draws=CALIBRATION_DRAWS,
            posterior=_Narrowed(calibration_engine.posterior),
        )
        assert float(np.min(narrowed["ks_pvalue"].values)) < 0.01
        # Under-dispersion is a *negative* area-to-curve in TARP's convention,
        # and the honest posterior's is not.
        atc = float(narrowed.attrs["ampere_calibration_tarp_atc"])
        assert atc < -0.02
        assert atc < float(calibration.attrs["ampere_calibration_tarp_atc"])
        # And the coverage curve says it in the units a reader quotes.
        levels = np.asarray(narrowed.coords["level"].values, dtype=float)
        squeezed = float(np.interp(0.68, levels, np.asarray(narrowed["coverage"].values)[:, 0]))
        honest = float(np.interp(0.68, levels, np.asarray(calibration["coverage"].values)[:, 0]))
        assert squeezed < 0.45 < honest

    def test_it_attaches_to_a_run_and_survives_netcdf(
        self, calibration_engine: Any, calibration: Any, tmp_path: Any
    ) -> None:
        from ampere.results import CALIBRATION_GROUP, attach_calibration, from_netcdf, to_netcdf

        run = calibration_engine.run(20, training={"max_num_epochs": 2})
        attach_calibration(run, calibration)
        assert CALIBRATION_GROUP in run.children
        path = tmp_path / "calibrated.nc"
        to_netcdf(run, path)
        stored = from_netcdf(path)[CALIBRATION_GROUP].dataset
        np.testing.assert_array_equal(
            stored["ranks"].values, np.asarray(calibration["ranks"].values)
        )
        np.testing.assert_allclose(
            stored["tarp_coverage"].values, np.asarray(calibration["tarp_coverage"].values)
        )

    def test_both_plots_render_from_the_group(self, calibration: Any) -> None:
        import matplotlib

        matplotlib.use("Agg")
        from ampere.results import figure_metadata, plot_coverage, plot_sbc_ranks

        figure = plot_sbc_ranks(calibration)
        assert any(key.endswith(".ks_pvalue") for key in figure_metadata(figure))
        axes = plot_coverage(calibration)
        metadata = figure_metadata(axes.get_figure())
        assert "tarp.area_to_curve" in metadata

    def test_skipping_tarp_leaves_the_group_without_its_curve(
        self, calibration_engine: Any
    ) -> None:
        result = calibration_engine.calibrate(count=20, posterior_draws=20, tarp=False)
        assert "tarp_coverage" not in result
        assert "ampere_calibration_tarp_atc" not in result.attrs

    def test_an_untrained_engine_refuses_and_says_to_run_first(self) -> None:
        engine = SBIEngine(bounded_problem(), method="npe", budget=20)
        with pytest.raises(EngineError, match="has not trained"):
            engine.calibrate(count=4, posterior_draws=4)


# ---------------------------------------------------------------------------
# 13. W3.4: truncated marginal ratio estimation
# ---------------------------------------------------------------------------


#: The three-round fit's settings, in one place because the report quotes them.
#: 800 simulations a round on the three-parameter joint problem, and an epoch
#: cap high enough that the *joint* estimator — the one the run's draws come
#: from — gets most of the way down rather than stopping early in its descent.
#: The marginal estimators reach ``sbi``'s own early stop well before it.
TMNRE_BUDGET = 800
TMNRE_EPOCHS = 120
TMNRE_DRAWS = 200

#: How far a **marginal** estimator's mean may sit from the emcee reference's,
#: in units of the reference's own standard deviation. Five three-round fits of
#: this problem, at budgets 800/1 000/2 000 and epoch caps 40/120/150/200, put
#: the worst marginal at 0.54; 1.0 is comfortably outside that and still far
#: inside the errors a bridge or truncation bug produces (a lost Jacobian moves
#: a lognormal mean by whole widths).
TMNRE_MEAN_TOLERANCE = 1.0

#: The central posterior mass the emcee reference's mean must fall inside, for
#: the **joint** draws.
#:
#: A coverage statement rather than a distance in reference widths, and the
#: reason is measured. The run's draws come from a *ratio* estimator — this
#: module's own section 5 already declines to hold ``method="nre"`` to accuracy
#: at CI budgets — and, since neither ``sbi`` nor this driver seeds torch's
#: global generator, two runs of the same seeded problem differ by however much
#: two random network initialisations differ. Across the five fits above the
#: joint estimator's worst parameter landed between 0.05 and 1.15 reference
#: standard deviations from the emcee mean, always on ``model.index``, whose
#: reference width (0.03) is the tightest of the three. A fixed multiple of
#: that width is therefore either flaky or vacuous.
#:
#: What is *not* noisy, in every one of those fits, is the pair of statements
#: below: the run's width sits inside W3.2's own 0.6-1.5 band, and the
#: reference's mean lies inside the run's own central interval. Together they
#: are strictly stronger than a location tolerance alone would be — an
#: over-confident estimator centred perfectly fails the width row, and a
#: correctly-wide estimator centred somewhere else fails this one — and neither
#: depends on how tight the parameter happens to be.
TMNRE_COVERAGE = 0.95

#: The truth ``test_engines`` generated the joint problem's data at, in the
#: order ``free_labels()`` gives. Every truncation box must contain it: a box
#: that does not has cut away the answer, and no later round can put it back.
TMNRE_TRUTH = np.array([2.0, -1.2, 1.0])


class TestTheTruncationPieces:
    """``_tmnre``'s arithmetic, checked without ``sbi`` so ``dev`` runs it.

    Every one of these is a claim the round loop rests on and none of them
    needs a trained network to state: that boxes nest, that a threshold
    crossing keeps both modes of a bimodal marginal, that a degenerate column
    does not raise.
    """

    def test_an_unbounded_box_accepts_everything_and_has_no_volume(self) -> None:
        box = TruncationBox.unbounded(3)
        assert box.size == 3
        assert box.accepts_everything
        assert box.log_volume == math.inf
        assert box.contains(np.zeros(3))

    def test_intersection_never_grows_a_box(self) -> None:
        wide = TruncationBox((-1.0, -1.0), (1.0, 1.0))
        narrow = TruncationBox((-0.5, -2.0), (0.25, 2.0))
        both = wide.intersect(narrow)
        assert both.lower == (-0.5, -1.0)
        assert both.upper == (0.25, 1.0)
        assert both.log_volume < wide.log_volume

    def test_a_box_knows_what_is_inside_it(self) -> None:
        box = TruncationBox((-1.0, 0.0), (1.0, 2.0))
        assert box.contains([0.0, 1.0])
        assert box.contains([-1.0, 2.0])  # edges are inside
        assert not box.contains([0.0, 2.5])

    def test_the_indicator_is_a_picklable_object_rather_than_a_closure(self) -> None:
        """A trained posterior carries its prior, and the store pickles it."""
        import pickle

        indicator = TruncationBox((-1.0,), (1.0,)).indicator()
        assert pickle.loads(pickle.dumps(indicator)).lower.tolist() == [-1.0]

    def test_the_box_records_both_parameterisations(self) -> None:
        problem = bounded_problem()
        recorded = TruncationBox((-0.5,), (0.5,)).to_dict(problem.constrain)
        # `slope` is positive, so its bijection is the log and a constrained
        # edge is exp(edge): a reader sees the box in their own coordinates.
        assert recorded["lower_constrained"][0] == pytest.approx(math.exp(-0.5), rel=1e-6)
        assert recorded["upper_constrained"][0] == pytest.approx(math.exp(0.5), rel=1e-6)

    def test_the_interval_spans_both_modes_of_a_bimodal_marginal(self) -> None:
        """A truncation that kept only the peak's own mode would be a silent bug."""
        grid = np.linspace(-5.0, 5.0, 401)
        density = np.log(
            np.exp(-0.5 * ((grid + 3.0) / 0.3) ** 2) + np.exp(-0.5 * ((grid - 3.0) / 0.3) ** 2)
        )
        lower, upper = interval_above(grid, density, 1e-4)
        assert lower < -3.0 < 3.0 < upper
        assert lower > -5.0 and upper < 5.0

    def test_a_flat_marginal_keeps_the_whole_grid(self) -> None:
        grid = np.linspace(0.0, 1.0, 11)
        assert interval_above(grid, np.zeros(11), 1e-4) == (0.0, 1.0)

    def test_a_degenerate_column_gives_a_flat_density_rather_than_raising(self) -> None:
        """Every draw identical — a tied or effectively fixed parameter."""
        samples = np.full((50, 1), 1.5)
        values = marginal_log_density(samples, np.linspace(1.0, 2.0, 7).reshape(-1, 1))
        assert values.shape == (7,)
        assert np.allclose(values, 0.0)

    def test_a_kernel_density_peaks_where_the_draws_are(self) -> None:
        rng = np.random.default_rng(7)
        samples = rng.normal(2.0, 0.5, 500).reshape(-1, 1)
        grid = np.linspace(0.0, 4.0, 81).reshape(-1, 1)
        values = marginal_log_density(samples, grid)
        assert float(grid[int(np.argmax(values)), 0]) == pytest.approx(2.0, abs=0.2)

    def test_the_marginals_are_enumerated_in_a_fixed_order(self) -> None:
        assert marginal_indices(3, 1) == ((0,), (1,), (2,))
        assert marginal_indices(3, 2) == ((0, 1), (0, 2), (1, 2))
        assert pair_labels(("a", "b", "c"), marginal_indices(3, 2)) == ("a|b", "a|c", "b|c")

    def test_a_degenerate_grid_is_widened_rather_than_left_with_no_width(self) -> None:
        nodes = grid_between(1.0, 1.0, 5)
        assert nodes.size == 5
        assert nodes[-1] > nodes[0]

    def test_the_pair_mesh_reshapes_back_to_its_two_axes(self) -> None:
        mesh = pair_mesh(np.array([0.0, 1.0]), np.array([10.0, 20.0, 30.0]))
        assert mesh.shape == (6, 2)
        assert mesh[:, 0].reshape(2, 3)[1].tolist() == [1.0, 1.0, 1.0]
        assert mesh[:, 1].reshape(2, 3)[0].tolist() == [10.0, 20.0, 30.0]


class TestWhatTMNRERefuses:
    """The new arguments belong to one method, and say so. No ``sbi`` needed."""

    def test_an_unknown_marginal_order_names_what_the_two_are_for(self) -> None:
        with pytest.raises(EngineError, match="corner plot"):
            SBIEngine(bounded_problem(), method="tmnre", marginals=3)

    def test_a_threshold_outside_zero_to_one_is_refused_by_name(self) -> None:
        with pytest.raises(EngineError, match="strictly between 0 and 1"):
            SBIEngine(bounded_problem(), method="tmnre", truncation_epsilon=2.0)

    def test_an_unknown_sampler_names_both(self) -> None:
        with pytest.raises(EngineError, match="rejection, mcmc"):
            SBIEngine(bounded_problem(), method="tmnre", sample_with="vi")

    def test_the_default_sampler_is_mcmc_and_rejection_stays_selectable(self) -> None:
        """Ruled 2026-09-10 on W3.4's measurement (353.6 s against 43.8 s)."""
        assert SBIEngine(bounded_problem(), method="tmnre").sample_with == "mcmc"
        assert (
            SBIEngine(bounded_problem(), method="tmnre", sample_with="rejection").sample_with
            == "rejection"
        )
        assert SBIEngine(bounded_problem(), method="npe").sample_with is None

    @pytest.mark.parametrize(
        ("options", "match"),
        [
            ({"marginals": 2}, "belongs to method='tmnre'"),
            ({"truncation_epsilon": 0.1}, "does not truncate anything"),
            ({"sample_with": "mcmc"}, "posterior_options= is where"),
        ],
    )
    def test_the_truncation_arguments_are_refused_for_the_other_methods(
        self, options: dict[str, Any], match: str
    ) -> None:
        """Accepted-and-ignored is the failure mode this refusal exists to stop."""
        with pytest.raises(EngineError, match=match):
            SBIEngine(bounded_problem(), method="npe", **options)

    def test_the_method_list_names_the_fourth(self) -> None:
        assert set(METHODS) == {"npe", "nle", "nre", "tmnre"}
        with pytest.raises(EngineError, match="no longer amortised"):
            SBIEngine(bounded_problem(), method="tmnr")


@pytest.fixture(scope="module")
def tmnre_run() -> Any:
    """One three-round TMNRE fit of the joint problem. The expensive fixture."""
    if not HAS_SBI:  # pragma: no cover - the class-level skip covers this
        pytest.skip("needs the 'sbi' extra")
    engine = SBIEngine(joint_problem(), method="tmnre", rounds=3, budget=TMNRE_BUDGET)
    return engine.run(draws=TMNRE_DRAWS, training={"max_num_epochs": TMNRE_EPOCHS})


@needs_sbi
class TestTMNRERecoversTheJointPosterior:
    """Accept criterion 1: the run recovers the toy joint posterior.

    The draws come from the joint estimator trained across the rounds and
    multiplied by the *truncated* prior, so this is a statement about the whole
    loop: a box that cut posterior mass, a ratio multiplied by the wrong
    proposal density, or a rejection sampler proposing outside the box would
    all show up here.

    The width band is W3.2's NPE band exactly. Location is a *coverage*
    statement rather than a distance in reference widths — see
    :data:`TMNRE_COVERAGE`, which records the five fits the choice was made
    from — and the marginal estimators, which are what TMNRE exists to train,
    are held to a distance (:data:`TMNRE_MEAN_TOLERANCE`) as well.
    """

    NAMES = ("model.norm", "model.index", "calibration")

    def test_the_posterior_is_in_the_constrained_space_under_merged_names(
        self, tmnre_run: Any
    ) -> None:
        posterior = tmnre_run["posterior"].dataset
        assert set(posterior.data_vars) == set(self.NAMES)
        assert posterior.sizes == {"chain": 1, "draw": TMNRE_DRAWS}
        assert float(np.asarray(posterior["model.norm"]).min()) > 0.0
        assert float(np.asarray(posterior["calibration"]).min()) > 0.0

    def test_it_writes_the_engine_neutral_approximation_family_and_proposal_density(
        self, tmnre_run: Any
    ) -> None:
        """W5.0: TMNRE runs the same code path as NPE/NLE/NRE, so it gets both too."""
        assert tmnre_run.attrs["ampere_approximation"] == "density_estimator"
        stats = tmnre_run["sample_stats"].dataset
        assert np.all(np.isfinite(np.asarray(stats["proposal_log_density"])))

    @pytest.mark.parametrize("name", NAMES)
    def test_the_posterior_covers_the_reference(
        self, name: str, tmnre_run: Any, emcee_reference: dict[str, tuple[float, float]]
    ) -> None:
        """Where emcee says the answer is, this run's posterior puts mass.

        See :data:`TMNRE_COVERAGE` for why this rather than a distance in
        reference widths. Read together with the width row below: right width
        *and* covering the reference is a stronger pair of claims than either
        alone, and neither is sensitive to how tight the parameter is.
        """
        drawn = np.asarray(tmnre_run["posterior"][name]).reshape(-1)
        mean, _ = emcee_reference[name]
        tail = 0.5 * (1.0 - TMNRE_COVERAGE)
        low, high = np.quantile(drawn, [tail, 1.0 - tail])
        assert low < mean < high

    @pytest.mark.parametrize("name", NAMES)
    def test_the_width_matches_the_reference(
        self, name: str, tmnre_run: Any, emcee_reference: dict[str, tuple[float, float]]
    ) -> None:
        """W3.2's own band: an over-confident estimator passes the mean row."""
        drawn = np.asarray(tmnre_run["posterior"][name])
        _, width = emcee_reference[name]
        assert 0.6 < float(drawn.std()) / width < 1.5

    @pytest.mark.parametrize("name", NAMES)
    def test_the_marginal_estimators_agree_with_the_reference_too(
        self, name: str, tmnre_run: Any, emcee_reference: dict[str, tuple[float, float]]
    ) -> None:
        """The method's *own* product, checked against the same reference.

        The joint estimator supplies the draws; the marginal ones are what
        TMNRE exists to train, and they are a separate estimate of the same
        quantity — so an implementation that trained them on the wrong columns,
        or reconstructed ``p(θ_i|x) ∝ r(θ_i, x) · q(θ_i)`` against the wrong
        ``q``, would disagree here while the ``posterior`` group looked fine.
        Read on the group's own grid, in the constrained coordinates the
        reference is in.
        """
        group = tmnre_run["marginals"].dataset
        index = list(group.coords["marginal_parameter"].values).index(name)
        grid = np.asarray(group["grid_constrained"])[index]
        density = np.exp(np.asarray(group["log_density"])[index])
        density /= np.trapezoid(density, grid)
        estimate = float(np.trapezoid(density * grid, grid))
        mean, width = emcee_reference[name]
        assert abs(estimate - mean) < TMNRE_MEAN_TOLERANCE * width

    def test_every_draw_still_carries_the_true_split_scored_on_the_numpy_path(
        self, tmnre_run: Any
    ) -> None:
        """Nothing about the emitted run is special: ``finish`` scored these."""
        stats = tmnre_run["sample_stats"].dataset
        prior = np.asarray(stats["log_prior"])
        likelihood = np.asarray(stats["log_likelihood"])
        assert np.all(np.isfinite(prior))
        assert np.all(np.isfinite(likelihood))
        assert np.asarray(stats["lp"]) == pytest.approx(prior + likelihood)
        assert np.all(np.isfinite(np.asarray(stats["ampere_sbi_log_prob"])))


@needs_sbi
class TestTheTruncationHistory:
    """Accept criterion 2: the boxes nest, they contain the truth, and say so."""

    @staticmethod
    def history(run: Any) -> list[dict[str, Any]]:
        import json

        return json.loads(run.attrs["ampere_sbi_truncation"])

    def test_there_is_one_record_per_round(self, tmnre_run: Any) -> None:
        records = self.history(tmnre_run)
        assert [record["round"] for record in records] == [1, 2, 3]
        assert [record["simulations"] for record in records] == [TMNRE_BUDGET] * 3

    def test_the_boxes_never_grow(self, tmnre_run: Any) -> None:
        volumes = [record["log_volume"] for record in self.history(tmnre_run)]
        assert all(later <= earlier + 1e-9 for earlier, later in itertools.pairwise(volumes))
        # And they genuinely shrank rather than merely not growing: three
        # rounds that changed nothing would pass a monotonicity check alone.
        assert volumes[-1] < volumes[0]

    def test_every_box_contains_the_truth_the_data_were_generated_at(self, tmnre_run: Any) -> None:
        truth = joint_problem().unconstrain(TMNRE_TRUTH)
        for record in self.history(tmnre_run):
            box = TruncationBox(tuple(record["lower"]), tuple(record["upper"]))
            assert box.contains(truth), f"round {record['round']} lost the truth"

    def test_the_boxes_are_nested_round_by_round(self, tmnre_run: Any) -> None:
        records = self.history(tmnre_run)
        for earlier, later in itertools.pairwise(records):
            assert np.all(np.asarray(later["lower"]) >= np.asarray(earlier["lower"]) - 1e-9)
            assert np.all(np.asarray(later["upper"]) <= np.asarray(earlier["upper"]) + 1e-9)

    def test_the_attrs_name_the_method_the_rounds_and_the_threshold(self, tmnre_run: Any) -> None:
        attrs = tmnre_run.attrs
        assert attrs["ampere_sbi_method"] == "tmnre"
        assert attrs["ampere_sbi_trainer"].startswith("NRE")
        assert attrs["ampere_sbi_rounds"] == 3
        assert attrs["ampere_sbi_marginals"] == 1
        assert attrs["ampere_sbi_truncation_epsilon"] == DEFAULT_TRUNCATION_EPSILON
        # The default sampler since 2026-09-10 (Peter's ruling on W3.4's timing).
        assert attrs["ampere_sbi_truncation_sampler"] == "mcmc"
        assert attrs["ampere_sbi_marginal_estimators"] == 3
        assert attrs["ampere_sbi_log_prob_kind"] == "unnormalised"

    def test_a_truncated_run_says_it_is_not_amortised(self, tmnre_run: Any) -> None:
        """The cost of truncation, recorded rather than left to be inferred."""
        assert tmnre_run.attrs["ampere_sbi_amortised"] == 0

    def test_a_single_round_fit_of_another_method_still_is(self) -> None:
        run = SBIEngine(bounded_problem(), method="npe", budget=120).run(
            draws=10, training={"max_num_epochs": 2}
        )
        assert run.attrs["ampere_sbi_amortised"] == 1


@pytest.fixture(scope="module")
def pair_run() -> Any:
    """``marginals=2`` at a smoke budget: the group's *shape* is what is tested.

    MCMC rather than rejection, and two rounds rather than three, because
    nothing below asks about accuracy — the accuracy claim is
    :func:`tmnre_run`'s — and a slice sampler over a box is the cheap way to
    get a posterior of the right shape.
    """
    if not HAS_SBI:  # pragma: no cover - the class-level skip covers this
        pytest.skip("needs the 'sbi' extra")
    engine = SBIEngine(
        joint_problem(), method="tmnre", rounds=2, budget=300, marginals=2, sample_with="mcmc"
    )
    run = engine.run(
        draws=40,
        training={"max_num_epochs": 8},
        posterior_options={"num_chains": 4, "warmup_steps": 10, "thin": 1},
    )
    return engine, run


@needs_sbi
class TestTheMarginalsGroup:
    """Accept criterion 3: the group carries 1-D and 2-D and survives netCDF."""

    def test_the_run_carries_a_marginals_group_beside_the_usual_five(
        self, pair_run: tuple[Any, Any]
    ) -> None:
        _, run = pair_run
        assert sorted(run.children) == [
            "constant_data",
            "log_likelihood",
            "marginals",
            "observed_data",
            "posterior",
            "sample_stats",
        ]

    def test_the_one_dimensional_marginals_are_one_row_per_parameter(
        self, pair_run: tuple[Any, Any]
    ) -> None:
        _, run = pair_run
        group = run["marginals"].dataset
        assert [str(name) for name in group.coords["marginal_parameter"].values] == [
            "model.norm",
            "model.index",
            "calibration",
        ]
        assert group["log_ratio"].dims == ("marginal_parameter", "marginal_node")
        assert group["log_density"].shape == (3, GRID_POINTS_1D)
        assert np.all(np.isfinite(np.asarray(group["log_ratio"])))
        # ``log_density`` is normalised to a maximum of zero, so it is the
        # marginal posterior up to its own constant and a plot can draw it.
        assert np.asarray(group["log_density"]).max(axis=1) == pytest.approx(np.zeros(3))

    def test_each_grid_spans_the_final_box_in_both_parameterisations(
        self, pair_run: tuple[Any, Any]
    ) -> None:
        engine, run = pair_run
        group = run["marginals"].dataset
        grid = np.asarray(group["grid"])
        assert engine.truncation is not None
        assert grid[:, 0] == pytest.approx(np.asarray(engine.truncation.lower))
        assert grid[:, -1] == pytest.approx(np.asarray(engine.truncation.upper))
        # `model.norm` and `calibration` are lognormal, so their constrained
        # grids are positive and monotone; `model.index` is unbounded and its
        # two grids are the same numbers.
        constrained = np.asarray(group["grid_constrained"])
        assert np.all(constrained[0] > 0.0)
        assert np.all(np.diff(constrained[0]) > 0.0)
        assert constrained[1] == pytest.approx(grid[1])

    def test_the_pairs_are_there_with_their_own_two_axes(self, pair_run: tuple[Any, Any]) -> None:
        _, run = pair_run
        group = run["marginals"].dataset
        assert [str(name) for name in group.coords["marginal_pair"].values] == [
            "model.norm|model.index",
            "model.norm|calibration",
            "model.index|calibration",
        ]
        assert group["pair_log_ratio"].dims == ("marginal_pair", "marginal_row", "marginal_column")
        assert group["pair_log_ratio"].shape == (3, GRID_POINTS_2D, GRID_POINTS_2D)
        assert np.all(np.isfinite(np.asarray(group["pair_log_density"])))
        assert run.attrs["ampere_sbi_marginals"] == 2
        assert run.attrs["ampere_sbi_marginal_estimators"] == 6

    def test_the_group_says_what_parameterisation_it_is_in(self, pair_run: tuple[Any, Any]) -> None:
        _, run = pair_run
        attrs = run["marginals"].attrs
        assert attrs["ampere_marginals_parameterisation"] == "unconstrained"
        assert attrs["ampere_marginals_order"] == 2
        assert attrs["ampere_marginals_schema_version"] == 1

    def test_it_survives_the_netcdf_round_trip(
        self, pair_run: tuple[Any, Any], tmp_path: Path
    ) -> None:
        from ampere.results import from_netcdf, to_netcdf

        _, run = pair_run
        back = from_netcdf(to_netcdf(run, tmp_path / "tmnre.nc"))
        assert "marginals" in back.children
        group, original = back["marginals"].dataset, run["marginals"].dataset
        for name in ("grid", "log_ratio", "log_density", "pair_log_ratio", "pair_log_density"):
            assert np.asarray(group[name]) == pytest.approx(np.asarray(original[name]))
        assert back["marginals"].attrs["ampere_marginals_order"] == 2

    def test_a_one_dimensional_run_carries_no_pairs(self, tmnre_run: Any) -> None:
        group = tmnre_run["marginals"].dataset
        assert "pair_log_ratio" not in group.data_vars
        assert "marginal_pair" not in group.coords


@pytest.fixture(scope="module")
def bounded_tmnre() -> Any:
    """A one-parameter TMNRE engine on the *rejection* sampler, for calibration.

    ``bounded_problem`` rather than the joint one, and ``sample_with=
    "rejection"`` explicitly -- it was the default until 2026-09-10 and is
    now the selectable alternative -- because what these rows are about is
    the rebuild path: a calibration check re-conditions the posterior at
    ``count`` fresh observations, and a rejection posterior pays its
    find-the-maximum stage at every one of them, so ``calibrate()`` rebuilds
    it under MCMC. One parameter keeps the whole thing inside a per-PR budget
    while still exercising the rebuild.
    """
    if not HAS_SBI:  # pragma: no cover - the class-level skip covers this
        pytest.skip("needs the 'sbi' extra")
    engine = SBIEngine(
        bounded_problem(), method="tmnre", rounds=2, budget=300, sample_with="rejection"
    )
    engine.run(draws=30, training={"max_num_epochs": 20})
    return engine


@needs_sbi
class TestTMNRECalibrationAndCaching:
    """Accept criterion 4, and W3.5's store against a run with more in it."""

    def test_calibrate_runs_on_a_tmnre_run_against_the_truncated_prior(
        self, bounded_tmnre: Any
    ) -> None:
        """The reference distribution is the box, and the group says so.

        Calibrating a truncated posterior against the *full* prior would report
        miscalibration for every truth the box excludes — an artefact of the
        method rather than a property of the estimator — so the fresh batch and
        the reference draws both come from the truncated prior.
        """
        result = bounded_tmnre.calibrate(count=10, posterior_draws=20, tarp=False)
        assert result.sizes["simulation"] == 10
        assert result.attrs["ampere_calibration_reference"] == "truncated_prior"
        assert result.attrs["ampere_calibration_method"] == "tmnre"
        assert np.all(np.isfinite(np.asarray(result["ks_pvalue"])))

    def test_it_checks_a_rejection_run_through_an_mcmc_posterior(self, bounded_tmnre: Any) -> None:
        """And says which, because it is not the object the run sampled with.

        Running SBC against the rejection posterior itself does not merely take
        longer: its fixed 10 000-proposal-draw maximisation is paid once per
        conditioning observation, out of a prior whose draws are themselves
        rejected against the box, so the check does not finish.
        """
        assert bounded_tmnre.sample_with == "rejection"
        result = bounded_tmnre.calibrate(count=6, posterior_draws=10, tarp=False)
        assert result.attrs["ampere_calibration_sampler"] == "mcmc"

    def test_an_untruncated_method_still_calibrates_against_the_prior(self) -> None:
        """And through its own posterior: nothing is rebuilt for NPE, NLE or NRE."""
        engine = SBIEngine(bounded_problem(), method="npe", budget=200)
        engine.run(draws=20, training={"max_num_epochs": 5})
        result = engine.calibrate(count=10, posterior_draws=10, tarp=False)
        assert result.attrs["ampere_calibration_reference"] == "prior"
        assert result.attrs["ampere_calibration_sampler"] == "as_run"

    def test_a_cache_hit_restores_the_marginals_group_and_the_history(self, tmp_path: Any) -> None:
        """A TMNRE run's answer is more than its posterior, so the store holds more."""
        from ampere.results.artefacts import ArtefactStore

        store = ArtefactStore(tmp_path / "artefacts")
        settings = dict(method="tmnre", rounds=2, budget=120, sample_with="mcmc", cache=store)
        options = dict(
            draws=10,
            training={"max_num_epochs": 2},
            posterior_options={"num_chains": 2, "warmup_steps": 5, "thin": 1},
        )
        first = SBIEngine(joint_problem(), **settings).run(**options)
        assert first.attrs["ampere_sbi_cache_hit"] == 0

        engine = SBIEngine(joint_problem(), **settings)
        second = engine.run(**options)
        assert second.attrs["ampere_sbi_cache_hit"] == 1
        assert second.attrs["ampere_sbi_simulations"] == 0
        # The group and the history come back, not only the posterior.
        assert "marginals" in second.children
        assert engine.truncation is not None
        assert len(engine.truncation_history) == 2
        assert second.attrs["ampere_sbi_truncation"] == first.attrs["ampere_sbi_truncation"]

    def test_a_different_threshold_is_a_miss(self, tmp_path: Any) -> None:
        """ε changes the box, the box changes the estimator, so it must not hit.

        ``artefact_key`` has no ``truncation_epsilon`` ingredient yet, so the
        engine folds TMNRE's settings into the architecture ingredient
        (``SBIEngine._key_architecture``) rather than letting two different
        runs collide on one digest. This row is what makes that a fact.
        """
        from ampere.results.artefacts import ArtefactStore

        store = ArtefactStore(tmp_path / "artefacts")
        settings = dict(method="tmnre", rounds=1, budget=100, sample_with="mcmc", cache=store)
        options = dict(
            draws=6,
            training={"max_num_epochs": 2},
            posterior_options={"num_chains": 2, "warmup_steps": 4, "thin": 1},
        )
        SBIEngine(joint_problem(), **settings).run(**options)
        again = SBIEngine(joint_problem(), truncation_epsilon=1e-3, **settings).run(**options)
        assert again.attrs["ampere_sbi_cache_hit"] == 0


@needs_sbi
class TestTheShippedTMNREExample:
    """``examples/sbi/tmnre_fit.py`` runs, at the budget it is asked for."""

    def test_it_runs_and_reports_the_box(self) -> None:
        script = EXAMPLES / "tmnre_fit.py"
        assert script.exists()
        finished = subprocess.run(
            [sys.executable, str(script), "--budget", "120", "--rounds", "2", "--draws", "20"],
            check=True,
            capture_output=True,
            text=True,
        )
        report = finished.stdout
        assert "tmnre on reference" in report
        assert "truncation box" in report
        assert "not amortised" in report


# ---------------------------------------------------------------------------
# 14. W3.15: torch is seeded from the problem too
# ---------------------------------------------------------------------------


@needs_sbi
class TestTorchIsSeededFromTheProblem:
    """Before this, ``problem.seed`` fixed the budget and nothing else.

    Neither ``ampere`` nor ``sbi`` seeded torch's own global generator, so a
    network's initial weights, its trainer's batch order and (for an
    MCMC-sampled posterior) its sampler's momentum draws varied between two
    runs of the *same* seeded problem — W3.4's carried finding. ``run()`` now
    seeds torch from the problem's own sub-stream before every step that
    spends it, so a run repeats bitwise, the network included, and records
    the seed it used as ``ampere_sbi_torch_seed``.
    """

    NPE_SETTINGS: ClassVar[dict[str, Any]] = {"method": "npe", "budget": 120}
    NPE_OPTIONS: ClassVar[dict[str, Any]] = {"draws": 10, "training": {"max_num_epochs": 2}}
    TMNRE_SETTINGS: ClassVar[dict[str, Any]] = {
        "method": "tmnre",
        "rounds": 1,
        "budget": 100,
        "sample_with": "mcmc",
    }
    TMNRE_OPTIONS: ClassVar[dict[str, Any]] = {
        "draws": 6,
        "training": {"max_num_epochs": 2},
        "posterior_options": {"num_chains": 2, "warmup_steps": 4, "thin": 1},
    }

    def test_npe_repeats_bitwise_from_the_same_problem_seed(self) -> None:
        first = SBIEngine(bounded_problem(SEED), **self.NPE_SETTINGS).run(**self.NPE_OPTIONS)
        second = SBIEngine(bounded_problem(SEED), **self.NPE_SETTINGS).run(**self.NPE_OPTIONS)
        assert np.asarray(first["posterior"]["model.slope"]) == pytest.approx(
            np.asarray(second["posterior"]["model.slope"]), abs=0.0
        )
        assert first.attrs["ampere_sbi_torch_seed"] == second.attrs["ampere_sbi_torch_seed"]

    def test_tmnre_repeats_bitwise_at_a_tiny_budget(self) -> None:
        first = SBIEngine(bounded_problem(SEED), **self.TMNRE_SETTINGS).run(**self.TMNRE_OPTIONS)
        second = SBIEngine(bounded_problem(SEED), **self.TMNRE_SETTINGS).run(**self.TMNRE_OPTIONS)
        assert np.asarray(first["posterior"]["model.slope"]) == pytest.approx(
            np.asarray(second["posterior"]["model.slope"]), abs=0.0
        )
        assert first.attrs["ampere_sbi_torch_seed"] == second.attrs["ampere_sbi_torch_seed"]

    def test_a_different_problem_seed_gives_a_different_torch_seed_and_draws(self) -> None:
        first = SBIEngine(bounded_problem(SEED), **self.NPE_SETTINGS).run(**self.NPE_OPTIONS)
        other = SBIEngine(bounded_problem(SEED + 1), **self.NPE_SETTINGS).run(**self.NPE_OPTIONS)
        assert first.attrs["ampere_sbi_torch_seed"] != other.attrs["ampere_sbi_torch_seed"]
        assert not np.allclose(
            np.asarray(first["posterior"]["model.slope"]),
            np.asarray(other["posterior"]["model.slope"]),
        )

    def test_a_cache_hit_samples_reproducibly_too(self, tmp_path: Path) -> None:
        """Training is skipped on a hit, so sampling must be pinned on its own."""
        from ampere.results.artefacts import ArtefactStore

        store = ArtefactStore(tmp_path / "artefacts")
        settings = dict(self.NPE_SETTINGS, cache=store)
        SBIEngine(bounded_problem(SEED), **settings).run(**self.NPE_OPTIONS)
        first = SBIEngine(bounded_problem(SEED), **settings).run(**self.NPE_OPTIONS)
        second = SBIEngine(bounded_problem(SEED), **settings).run(**self.NPE_OPTIONS)
        assert first.attrs["ampere_sbi_cache_hit"] == 1
        assert second.attrs["ampere_sbi_cache_hit"] == 1
        assert np.asarray(first["posterior"]["model.slope"]) == pytest.approx(
            np.asarray(second["posterior"]["model.slope"]), abs=0.0
        )
        assert first.attrs["ampere_sbi_torch_seed"] == second.attrs["ampere_sbi_torch_seed"]

    def test_the_attr_is_absent_for_an_unseeded_problem(self) -> None:
        """``seed=None`` asks for fresh randomness, torch included -- not a fixed one."""
        run = SBIEngine(bounded_problem(seed=None), **self.NPE_SETTINGS).run(**self.NPE_OPTIONS)
        assert "ampere_sbi_torch_seed" not in run.attrs

    def test_the_attr_is_an_int_for_a_seeded_problem(self) -> None:
        run = SBIEngine(bounded_problem(SEED), **self.NPE_SETTINGS).run(**self.NPE_OPTIONS)
        assert isinstance(run.attrs["ampere_sbi_torch_seed"], int)


@needs_sbi
class TestCalibrateReseedsTorch:
    """W5.28(j): ``calibrate()`` reseeds torch as ``run()`` does, for the same reason.

    Before this, the SBC batch's truncated-prior draws and ``sbi``'s own
    ``run_sbc``/``run_tarp`` sampled through torch's (and, for a slice-MCMC
    posterior, numpy's legacy global) generator with no reseed of their own --
    found live in ``TestServingANamedArtefact``'s workaround, which had to
    save and restore torch's RNG state by hand around a ``calibrate()`` call
    to stop it perturbing this class's own TARP thresholds when run earlier in
    the same session. Two calibrations of the *same trained posterior*, with
    *different* torch global state walking in, must now agree bitwise --
    calibrate() must not care what state it inherits, the same guarantee
    ``run()``'s own final draw already has.
    """

    SETTINGS: ClassVar[dict[str, Any]] = {"method": "npe", "budget": 120}
    OPTIONS: ClassVar[dict[str, Any]] = {"draws": 10, "training": {"max_num_epochs": 2}}

    def test_two_calibrations_of_the_same_posterior_repeat_bitwise(self) -> None:
        import torch

        first_engine = SBIEngine(bounded_problem(SEED), **self.SETTINGS)
        first_engine.run(**self.OPTIONS)
        torch.manual_seed(1)  # an arbitrary global state walking into calibrate()
        first = first_engine.calibrate(count=12, posterior_draws=8, tarp=True)

        second_engine = SBIEngine(bounded_problem(SEED), **self.SETTINGS)
        second_engine.run(**self.OPTIONS)
        torch.manual_seed(999)  # a *different* arbitrary state
        second = second_engine.calibrate(count=12, posterior_draws=8, tarp=True)

        np.testing.assert_array_equal(np.asarray(first["ranks"]), np.asarray(second["ranks"]))
        np.testing.assert_allclose(
            np.asarray(first["tarp_coverage"]), np.asarray(second["tarp_coverage"]), atol=0.0
        )

    def test_calibrate_restores_torchs_global_state_on_exit(self) -> None:
        """The other half of the contract ``_seeded`` states: nothing leaks out."""
        import torch

        engine = SBIEngine(bounded_problem(SEED), **self.SETTINGS)
        engine.run(**self.OPTIONS)
        torch.manual_seed(20260922)
        before = torch.get_rng_state().clone()
        engine.calibrate(count=12, posterior_draws=8, tarp=True)
        after = torch.get_rng_state()
        torch.testing.assert_close(before, after)


# ---------------------------------------------------------------------------
# W5.10: amortisation over the observation context
# ---------------------------------------------------------------------------

#: The context prior the amortised fixture trains under: half to twice the
#: observed error bars, log-uniformly. Stated once because three rows quote it.
CONTEXT_LOW, CONTEXT_HIGH = 0.5, 2.0

#: Where the acceptance row asks the question. ``1.5`` is inside
#: ``[0.5, 2.0]``; ``12`` is a long way outside it, and deliberately so — the
#: claim under test is not "coverage degrades the moment you leave the range"
#: (it should not, and a posterior that brittle would be useless) but "a noise
#: level the prior never showed the network is not one it can be believed at".
COVERED_FACTOR = 1.5
UNCOVERED_FACTOR = 12.0

#: The same floor ``sbi``'s own checks warn below (see ``CALIBRATION_COUNT``),
#: at the smaller budget this pair of checks can afford twice over.
CONTEXT_CALIBRATION_COUNT = 100
CONTEXT_CALIBRATION_DRAWS = 100


class TestTheContextPriorRefusals:
    """The vocabulary, in ``dev``: what ``context=`` accepts and what it does not."""

    def test_something_that_is_not_a_context_prior_is_refused_by_name(self) -> None:
        with pytest.raises(EngineError, match="ContextPrior"):
            SBIEngine(bounded_problem(), context="noisier")

    def test_a_shipped_prior_is_accepted_and_kept(self) -> None:
        prior = ScaledSigma(CONTEXT_LOW, CONTEXT_HIGH)
        engine = SBIEngine(bounded_problem(), context=prior)
        assert engine.context is prior

    def test_no_context_is_still_the_default(self) -> None:
        assert SBIEngine(bounded_problem()).context is None


@pytest.fixture(scope="module")
def amortised() -> Any:
    """One NPE engine trained under a sigma-pattern context prior (**W5.10**).

    ``layout="set"`` with the set embedding, and that is the point rather than
    a detail: the ``"flat"`` summary is the observed *values* and nothing else,
    so a network trained under it cannot see which noise level it is looking
    at however hard the simulator varies one. The set packing carries each
    sample's ``log sigma``, so this network can — with no layout change, which
    is the claim the item makes.

    The fitted run is returned beside the engine rather than left to a row to
    produce, because ``run()`` **retrains**: a test calling it again would
    quietly replace the network every other row here is about.
    """
    engine = SBIEngine(
        bounded_problem(),
        method="npe",
        budget=1500,
        embedding="set",
        layout="set",
        context=ScaledSigma(CONTEXT_LOW, CONTEXT_HIGH),
    )
    run = engine.run(200, training={"max_num_epochs": 150})
    return engine, run


@pytest.fixture(scope="module")
def covered_calibration(amortised: Any) -> Any:
    """SBC and TARP at a rescale the training context prior covers."""
    return amortised[0].calibrate(
        count=CONTEXT_CALIBRATION_COUNT,
        posterior_draws=CONTEXT_CALIBRATION_DRAWS,
        context=ScaledSigma(COVERED_FACTOR, COVERED_FACTOR),
    )


@pytest.fixture(scope="module")
def uncovered_calibration(amortised: Any) -> Any:
    """The same check at a rescale it does not."""
    return amortised[0].calibrate(
        count=CONTEXT_CALIBRATION_COUNT,
        posterior_draws=CONTEXT_CALIBRATION_DRAWS,
        context=ScaledSigma(UNCOVERED_FACTOR, UNCOVERED_FACTOR),
    )


@needs_sbi
class TestAmortisationOverTheObservationContext:
    """W5.10's acceptance row: calibrated where the prior reaches, and not beyond.

    Both arms are pinned as **inequalities** rather than as numbers, and the
    pair is the claim: a posterior trained under a sigma-pattern prior stays
    calibrated at a noise level inside that prior, and its coverage degrades
    measurably at one outside it. Either half alone proves nothing — a
    posterior that failed everywhere would pass the second, and one that
    ignored the data would pass the first.
    """

    def test_the_run_records_the_context_prior_and_its_digest(self, amortised: Any) -> None:
        run = amortised[1]
        recorded = json.loads(run.attrs["ampere_sbi_context"])
        assert recorded == ScaledSigma(CONTEXT_LOW, CONTEXT_HIGH).describe()
        assert len(run.attrs["ampere_sbi_context_hash"]) == 32

    def test_a_contextless_run_still_records_none(self) -> None:
        engine = SBIEngine(bounded_problem(), method="npe", budget=40)
        run = engine.run(5, training={"max_num_epochs": 2})
        assert run.attrs["ampere_sbi_context"] == "none"
        assert run.attrs["ampere_sbi_context_hash"] == ""

    def test_the_budget_was_actually_drawn_under_the_prior(self, amortised: Any) -> None:
        """The simulated observations carry the drawn sigma, not the observed one."""
        batch = amortised[0].batch
        assert batch is not None
        sigmas = {
            float(np.asarray(draw.observations["default"].uncertainty)[0]) for draw in batch.usable
        }
        assert len(sigmas) > 1
        assert min(sigmas) < 0.3 < max(sigmas)

    def test_it_stays_calibrated_at_a_rescale_the_prior_covers(
        self, covered_calibration: Any
    ) -> None:
        """Arm one: SBC ranks uniform and TARP flat at ``COVERED_FACTOR``.

        Measured at this budget: ``atc = -0.007``, TARP's own KS ``p = 1.0``.
        The thresholds are several times that, for the reason this file's
        header gives about every quantitative row here — loose enough that a
        correct fit passes essentially always, and far tighter than the
        failure the second arm produces.
        """
        assert covered_calibration.attrs["ampere_calibration_context"] != "none"
        assert float(np.min(covered_calibration["ks_pvalue"].values)) > 0.01
        assert abs(float(covered_calibration.attrs["ampere_calibration_tarp_atc"])) < 0.02
        assert float(covered_calibration.attrs["ampere_calibration_tarp_ks_pvalue"]) > 0.05

    def test_its_coverage_degrades_at_one_the_prior_does_not(
        self, covered_calibration: Any, uncovered_calibration: Any
    ) -> None:
        """Arm two: the same network, the same check, a noise level it never saw.

        Under-dispersion is a **negative** area-to-curve in TARP's convention:
        the posterior is too narrow, because it is reading error bars an order
        of magnitude smaller than the ones the observation actually has.

        The two thresholds are set from the measurement rather than chosen a
        priori, the way this file's other quantitative rows are. At the
        budget above the covered arm scores ``atc = -0.007`` and the uncovered
        one ``-0.046`` — a separation of ``0.039`` — so ``0.02`` is a little
        under half the effect and comfortably above the Monte Carlo error of
        a 100-simulation TARP curve. What is *not* asserted is a particular
        size of failure: how badly a network extrapolates outside its training
        context is a property of that network, and pinning it would be pinning
        noise.
        """
        covered = float(covered_calibration.attrs["ampere_calibration_tarp_atc"])
        uncovered = float(uncovered_calibration.attrs["ampere_calibration_tarp_atc"])
        # It degrades, by a margin, and in the direction under-dispersion
        # takes -- both halves matter, because a posterior that went *wider*
        # outside its training range would also move the number.
        assert uncovered < covered - 0.02
        assert uncovered < -0.02

    def test_calibrate_inherits_the_runs_prior_by_default(self, amortised: Any) -> None:
        report = amortised[0].calibrate(count=20, posterior_draws=20, tarp=False)
        recorded = json.loads(report.attrs["ampere_calibration_context"])
        assert recorded == ScaledSigma(CONTEXT_LOW, CONTEXT_HIGH).describe()

    def test_calibrate_at_no_context_is_spelled_none(self, amortised: Any) -> None:
        report = amortised[0].calibrate(count=20, posterior_draws=20, tarp=False, context=None)
        assert report.attrs["ampere_calibration_context"] == "none"


@needs_sbi
class TestTheFilmConditioningRoute:
    """W5.10's opt-in second route: off by default, and it trains when on."""

    @staticmethod
    def built(spec: Any) -> Any:
        import torch

        from ampere.inference._sbi import _embedding_of

        torch.manual_seed(SEED)
        problem = two_dataset_problem()
        layout = EncodingLayout.from_datasets(problem.datasets)
        resolved = _embedding_of(
            spec, torch=torch, features=0, free_size=problem.free_size, layout=layout
        )
        return layout, resolved.module

    def test_it_is_off_unless_asked_for(self) -> None:
        _, module = self.built("set")
        assert module.film is None
        _, transformer = self.built("transformer")
        assert transformer.film is None

    def test_asking_for_it_builds_a_conditioner(self) -> None:
        _, module = self.built({"type": "set", "film": True})
        assert module.film is not None
        _, transformer = self.built({"type": "transformer", "film": True})
        assert transformer.film is not None

    def test_it_starts_as_the_identity(self) -> None:
        """Switching the option on must not perturb a run before it has learnt.

        The conditioner's output layer is zero-initialised, so ``gamma`` and
        ``beta`` are zero and ``x * (1 + 0) + 0`` is ``x``: an untrained FiLM
        embedding produces exactly what the same embedding without one does.
        """
        import torch

        problem = two_dataset_problem()
        layout, plain = self.built("set")
        _, conditioned = self.built({"type": "set", "film": True})
        conditioned.net = plain.net
        x = torch.as_tensor(
            np.asarray(encode_observations(problem.datasets, layout=layout).values),
            dtype=torch.float32,
        )
        plain.eval()
        conditioned.eval()
        with torch.no_grad():
            assert torch.allclose(plain(x), conditioned(x), atol=1e-6)

    def test_the_set_embedding_trains_with_film(self) -> None:
        engine = SBIEngine(
            two_dataset_problem(),
            method="npe",
            budget=80,
            layout="set",
            embedding={"type": "set", "film": True},
        )
        run = engine.run(draws=10, training={"max_num_epochs": 5})
        assert run["posterior"].dataset.sizes["draw"] == 10
        assert run.attrs["ampere_sbi_embedding"] == "set"

    def test_the_transformer_embedding_trains_with_film(self) -> None:
        engine = SBIEngine(
            two_dataset_problem(),
            method="npe",
            budget=60,
            layout="set",
            embedding={"type": "transformer", "film": True},
        )
        run = engine.run(draws=8, training={"max_num_epochs": 3})
        assert run["posterior"].dataset.sizes["draw"] == 8

    def test_a_film_embedding_pickles(self) -> None:
        """The wrappers carry the conditioner through ``__reduce__``, as W3.5's cache needs."""
        import pickle

        import torch

        _, module = self.built({"type": "set", "film": True})
        restored = pickle.loads(pickle.dumps(module))
        assert restored.film is not None
        assert restored.layout.hash == module.layout.hash
        assert isinstance(restored.film.net, torch.nn.Sequential)


@needs_sbi
class TestTheContextIsInTheCacheKey:
    """W5.10: the store must not serve a network trained under another context.

    The three rows are the three ways the key could have got this wrong, and
    each of them was true before ``context=`` became an ingredient: a
    context-amortised run served a context-free network, a context-free run
    served a context-amortised one, and a run under one prior served a
    network trained under another. Cheap to check — the digest is computed
    before any simulation, so these run at a two-epoch budget.
    """

    @staticmethod
    def settings(**overrides: Any) -> dict[str, Any]:
        base: dict[str, Any] = {"method": "npe", "budget": 40, "layout": "set", "embedding": "set"}
        base.update(overrides)
        return base

    def test_a_context_run_misses_a_context_free_entry(self, tmp_path: Any) -> None:
        from ampere.results.artefacts import ArtefactStore

        store = ArtefactStore(tmp_path / "artefacts")
        plain = SBIEngine(bounded_problem(), cache=store, **self.settings()).run(
            draws=5, training={"max_num_epochs": 2}
        )
        amortised = SBIEngine(
            bounded_problem(),
            cache=store,
            context=ScaledSigma(CONTEXT_LOW, CONTEXT_HIGH),
            **self.settings(),
        )
        run = amortised.run(draws=5, training={"max_num_epochs": 2})
        assert run.attrs["ampere_sbi_cache_hit"] == 0
        # Filed apart, which is the property that stops the second run from
        # ever being handed the first's network. *Which* field differs is
        # named by ArtefactStore.diff -- asserted, with the pair of values, in
        # tests/results/test_artefacts.py.
        assert run.attrs["ampere_sbi_cache_key"] != plain.attrs["ampere_sbi_cache_key"]

    def test_the_same_context_prior_is_a_hit(self, tmp_path: Any) -> None:
        from ampere.results.artefacts import ArtefactStore

        store = ArtefactStore(tmp_path / "artefacts")
        settings = self.settings(cache=store, context=ScaledSigma(CONTEXT_LOW, CONTEXT_HIGH))
        SBIEngine(bounded_problem(), **settings).run(draws=5, training={"max_num_epochs": 2})
        run = SBIEngine(bounded_problem(), **settings).run(draws=5, training={"max_num_epochs": 2})
        assert run.attrs["ampere_sbi_cache_hit"] == 1

    def test_two_different_priors_do_not_share_an_entry(self, tmp_path: Any) -> None:
        from ampere.results.artefacts import ArtefactStore

        store = ArtefactStore(tmp_path / "artefacts")
        SBIEngine(
            bounded_problem(), cache=store, context=ScaledSigma(0.5, 2.0), **self.settings()
        ).run(draws=5, training={"max_num_epochs": 2})
        wider = SBIEngine(
            bounded_problem(), cache=store, context=ScaledSigma(0.1, 10.0), **self.settings()
        )
        run = wider.run(draws=5, training={"max_num_epochs": 2})
        assert run.attrs["ampere_sbi_cache_hit"] == 0
