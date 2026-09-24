"""W5.13: population inference by reweighting archived fits, end to end.

The toy: 200 objects, each with a truth :math:`\\theta_i \\sim \\mathcal{N}(\\mu=1.0,
\\tau=0.3)` and one noisy Gaussian datum :math:`d_i = \\theta_i +
\\mathcal{N}(0, \\sigma=0.2)`, each fit independently and archived under a
*wide* interim prior with :class:`~ampere.inference.EmceeEngine`. Every claim
this module's docstring makes is exercised against that toy:

1. The reweighted :math:`(\\mu, \\tau)` posterior recovers the population truth
   and is narrower on :math:`\\mu` than any single object's own posterior.
2. The in-memory (:class:`~ampere.results.DataTreeRunColumns`) and
   netCDF-directory (:func:`~ampere.results.runs_from_netcdf_directory`)
   readers agree bitwise.
3. An object whose draws sit far from the population collapses its effective
   sample size, and :func:`~ampere.results.fit_population` refuses by name.
4. Runs that do not share ``ampere_spec_hash`` are refused by name.
5. :class:`~ampere.results.RunColumns` is a structural protocol: a plain
   dataclass that is not an :class:`xarray.DataTree` satisfies it.
6. An object fitted by :class:`~ampere.inference.VIEngine` (mean-field guide,
   sbi environment) is reweighted through its ``proposal_log_density`` and
   agrees with the same object fitted by ``EmceeEngine``.
7. A two-parameter object model (a nuisance parameter with a narrow,
   non-flat prior, alongside the reweighted one) recovers the same
   population posterior as the one-parameter case on the same objects, and
   the truth inside its own central 95 % -- the regression test for
   dividing by the named parameter's marginal interim prior rather than
   the run's stored *joint* one (§ "why it needs the marginal interim
   prior" above), a fix a single-parameter model's toy cannot exercise
   because the two priors coincide there.
8. **W5.22**: :func:`~ampere.results.fit_population` with no ``interim_prior``
   reads the named parameter's marginal prior back off the archive's own
   provenance (schema 8, ``ampere_free_priors``) and reproduces the same
   recovery as claim 1; a supplied ``interim_prior`` that disagrees with the
   stored one is refused by name; an archive with no stored priors and no
   supplied ``interim_prior`` is refused by name too.

Fits are deliberately tiny (a few hundred steps, few walkers): the 200
per-object fits are timed in :class:`TestTheTwoHundredObjectsFit` and
asserted to complete inside a generous ceiling -- not a budget from
``WORK_ITEMS.md`` (W5.13 sets none) and not tight even under the CPU
contention this machine's concurrent gates create, but tight enough that a
regression which made per-object fitting expensive by orders of magnitude
would fail a test rather than only a stopwatch.

**W5.22**: fitting and timing the full 200-object archive is minutes of
wall clock a per-PR gate should not pay for by default, so
:class:`TestTheTwoHundredObjectsFit` carries a ``population_full`` marker
(``pyproject.toml``, skipped by default -- ``tests/results/conftest.py``;
run with ``pytest -m population_full``), on the ``image_full``/``m2_full``
pattern. Every other test in this module -- the reweighted-posterior
recovery, the netCDF round trip, the ESS refusal, the marginal-interim-prior
regression -- shares a *reduced* archive of :data:`N_OBJECTS_REDUCED` (50)
objects instead, a bitwise prefix of the same 200-draw realisation, since
none of their assertions actually need the full count; only the timing
regression itself does, and :class:`TestTheReducedObjectsFit` pins the same
guard at the reduced count, always in the default gate, at a
proportionally looser ceiling.

The population truth (:math:`\\mu=1.0,\\ \\tau=0.3`) is the *generating*
hyperparameter, not the empirical mean/std of any one realised sample of
200 draws from it; SEED is fixed at a value whose realised sample is close
to both (mean 0.998, std 0.298 for the 200 draws :func:`_truths_and_data`
makes) so that the single-realisation "truth inside the central 95 %"
check is not a coin flip against sampling variance in the population
draws themselves -- confirmed analytically (bypassing ``emcee`` with exact
conjugate per-object posteriors) before being fixed here, and cross-checked
against several other seeds where the same central-95 % check narrowly and
correctly fails precisely because that seed's realised sample sits further
from the generating hyperparameters.
"""

from __future__ import annotations

import importlib
import time
import warnings
from typing import Any, ClassVar

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.core import (
    Dataset,
    FittingProblem,
    GaussianFamily,
    Likelihood,
    Model,
    ModelResult,
    Parameter,
    Spectrum,
)
from ampere.core.exceptions import ResultsError
from ampere.inference import EmceeEngine
from ampere.results import (
    ATTR_PREFIX,
    DataTreeRunColumns,
    GaussianPopulationModel,
    RunColumns,
    fit_population,
    from_netcdf,
    runs_from_netcdf_directory,
    to_netcdf,
)
from ampere.results.population import _self_normalised_log_weights

pytest.importorskip("arviz", reason="ampere.results needs arviz")

SEED = 28  # realised sample close to (MU_TRUTH, TAU_TRUTH); see module docstring
N_OBJECTS = 200  # the full archive; population_full only (TestTheTwoHundredObjectsFit)
#: The default-gate archive size: a bitwise prefix of the N_OBJECTS realisation
#: (same seed, same per-object seeds), so every test below that does not
#: itself need the full count shares this smaller one instead (W5.22).
#: TestTheMarginalInterimPriorFix's own N reuses this constant, since its
#: 50-object subset predates this item and is exactly this size already.
N_OBJECTS_REDUCED = 50
MU_TRUTH = 1.0
TAU_TRUTH = 0.3
SIGMA = 0.2
INTERIM_PRIOR = st.norm(0.0, 10.0)  # ConstantModel's theta prior -- fit_population's marginal pi_0


# ---------------------------------------------------------------------------
# The toy: one free scalar, one Gaussian datum
# ---------------------------------------------------------------------------


class ConstantModel(Model):
    """The smallest model this contract can fit: ``theta``, read straight back."""

    def __init__(self, prior: Any = None) -> None:
        self.register_buffer("grid", np.array([1.0]), unit=u.micron)
        self.register_parameter(
            Parameter("theta", prior if prior is not None else st.norm(0.0, 10.0))
        )

    def evaluate(self, **values: Any) -> ModelResult:
        ctx = self.context(values)
        return ModelResult(
            Spectrum(ctx["grid"] * u.micron, np.full_like(ctx["grid"], ctx["theta"]) * u.Jy)
        )


NUISANCE_PRIOR = st.norm(0.0, 0.05)


class TwoParameterModel(Model):
    """``theta`` (reweighted) plus ``phi``, a pure nuisance the likelihood never touches.

    ``phi``'s prior is narrow and non-flat (``NUISANCE_PRIOR``) precisely so
    its marginal log-density varies materially from draw to draw: that
    per-draw variation is the uncancelled factor the pre-fix code left in
    the reweighting ratio whenever an object had more than one free
    parameter (module docstring, "why it needs the marginal interim
    prior"). ``phi`` never entering ``evaluate`` means its exact posterior
    *is* its declared prior -- the likelihood carries no information about
    it at all -- which is what makes the toy a clean regression check
    rather than merely "a second parameter exists".
    """

    def __init__(self, theta_prior: Any = None, phi_prior: Any = None) -> None:
        self.register_buffer("grid", np.array([1.0]), unit=u.micron)
        self.register_parameter(
            Parameter("theta", theta_prior if theta_prior is not None else st.norm(0.0, 10.0))
        )
        self.register_parameter(
            Parameter("phi", phi_prior if phi_prior is not None else NUISANCE_PRIOR)
        )

    def evaluate(self, **values: Any) -> ModelResult:
        ctx = self.context(values)
        return ModelResult(
            Spectrum(ctx["grid"] * u.micron, np.full_like(ctx["grid"], ctx["theta"]) * u.Jy)
        )


def _two_parameter_problem(datum: float, seed: int) -> FittingProblem:
    observed = Spectrum(
        np.array([1.0]) * u.micron,
        np.array([datum]) * u.Jy,
        uncertainty=np.array([SIGMA]) * u.Jy,
    )
    return FittingProblem(TwoParameterModel(), [Dataset(observed)], seed=seed)


def _fit_two_parameter_object(
    datum: float, seed: int, *, walkers: int = 8, steps: int = 300, burn_in: int = 100
) -> Any:
    problem = _two_parameter_problem(datum, seed)
    return EmceeEngine(problem, walkers=walkers).run(steps=steps, burn_in=burn_in)


def _truths_and_data(seed: int, n: int) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    truths = rng.normal(MU_TRUTH, TAU_TRUTH, n)
    data = truths + rng.normal(0.0, SIGMA, n)
    return truths, data


def _object_problem(datum: float, seed: int, *, prior: Any = None) -> FittingProblem:
    observed = Spectrum(
        np.array([1.0]) * u.micron,
        np.array([datum]) * u.Jy,
        uncertainty=np.array([SIGMA]) * u.Jy,
    )
    return FittingProblem(ConstantModel(prior), [Dataset(observed)], seed=seed)


def _fit_object(
    datum: float,
    seed: int,
    *,
    prior: Any = None,
    walkers: int = 8,
    steps: int = 300,
    burn_in: int = 100,
) -> Any:
    problem = _object_problem(datum, seed, prior=prior)
    return EmceeEngine(problem, walkers=walkers).run(steps=steps, burn_in=burn_in)


def _population_model() -> GaussianPopulationModel:
    return GaussianPopulationModel(
        Parameter("mu", st.norm(0.0, 5.0)),
        Parameter("tau", st.halfnorm(scale=2.0)),
    )


def _fit_archive(n: int) -> tuple[list[Any], float]:
    """Fit and time *n* per-object runs, timed as one archive.

    :func:`_truths_and_data` draws from a freshly-seeded generator, so its
    first *n* truths/data are bitwise the same as the first *n* of a larger
    draw from the same seed (e.g. :data:`N_OBJECTS`'s 200) -- fitting *n*
    directly, rather than fitting 200 and slicing, is the identical data for
    a fraction of the cost (W5.22).
    """
    _, data = _truths_and_data(SEED, n)
    start = time.perf_counter()
    runs = [_fit_object(float(datum), seed=SEED + index) for index, datum in enumerate(data)]
    elapsed = time.perf_counter() - start
    return runs, elapsed


@pytest.fixture(scope="module")
def timed_object_runs() -> tuple[list[Any], float]:
    """The default-gate archive: :data:`N_OBJECTS_REDUCED` per-object runs, timed.

    A bitwise prefix of the full :data:`N_OBJECTS` realisation used under
    ``population_full`` -- see :class:`TestTheTwoHundredObjectsFit` and
    :class:`TestTheReducedObjectsFit`.
    """
    return _fit_archive(N_OBJECTS_REDUCED)


@pytest.fixture(scope="module")
def object_runs(timed_object_runs: tuple[list[Any], float]) -> list[Any]:
    runs, _ = timed_object_runs
    return runs


@pytest.fixture(scope="module")
def full_timed_object_runs() -> tuple[list[Any], float]:
    """The full 200-object archive, and the wall clock it took to build it.

    ``population_full`` only (skipped by default -- W5.22): only
    :class:`TestTheTwoHundredObjectsFit` needs the full count, so nothing
    else in this module requests this fixture.
    """
    return _fit_archive(N_OBJECTS)


@pytest.mark.population_full
class TestTheTwoHundredObjectsFit:
    """The population_full row: the full archive, fit and timed (W5.22).

    Skipped by default (``pyproject.toml``'s ``population_full`` marker,
    ``tests/results/conftest.py``'s hook); run with
    ``pytest -m population_full``. :class:`TestTheReducedObjectsFit` pins
    the same guard, unmarked, at the default-gate archive size.
    """

    def test_two_hundred_objects_fit_without_a_gross_regression(
        self, full_timed_object_runs: tuple[list[Any], float]
    ) -> None:
        """A coarse regression guard, not a wall-clock SLA.

        No budget for these fits is set by ``WORK_ITEMS.md``. The ceiling
        here is loose enough to absorb this machine's own CPU contention
        (concurrent five-suite gates in the main checkout routinely halve a
        single-threaded process's throughput) and is meant to catch an
        order-of-magnitude regression in per-object fitting cost, not to
        enforce a specific wall clock.
        """
        runs, elapsed = full_timed_object_runs
        assert len(runs) == N_OBJECTS
        print(f"\n200 emcee single-object fits: {elapsed:.2f} s wall clock.")
        assert elapsed < 600.0, (
            f"200 tiny emcee fits took {elapsed:.2f} s -- more than ten minutes, which is no "
            f"longer explainable by ordinary machine contention alone."
        )


@pytest.mark.study
class TestTheReducedObjectsFit:
    """The default-gate sibling of :class:`TestTheTwoHundredObjectsFit` (W5.22).

    :data:`N_OBJECTS_REDUCED` objects instead of the full :data:`N_OBJECTS`:
    pins the same regression guard, always in the default gate, at a
    looser ceiling stated as a fraction of the full row's -- 1/4 the
    objects, a little more than 1/4 the time budget, to absorb this
    machine's own contention the same way the full row's ceiling does.
    """

    def test_fifty_objects_fit_without_a_gross_regression(
        self, timed_object_runs: tuple[list[Any], float]
    ) -> None:
        runs, elapsed = timed_object_runs
        assert len(runs) == N_OBJECTS_REDUCED
        print(f"\n{N_OBJECTS_REDUCED} emcee single-object fits: {elapsed:.2f} s wall clock.")
        assert elapsed < 200.0, (
            f"{N_OBJECTS_REDUCED} tiny emcee fits took {elapsed:.2f} s -- more than the "
            f"proportionally-scaled ceiling, which is no longer explainable by ordinary "
            f"machine contention alone."
        )


# ---------------------------------------------------------------------------
# 1. The reweighted population posterior
# ---------------------------------------------------------------------------


@pytest.mark.study
class TestReweightedPopulationPosterior:
    @pytest.fixture(scope="class")
    def result(self, object_runs: list[Any]) -> Any:
        columns = [DataTreeRunColumns(run) for run in object_runs]
        return fit_population(
            columns,
            "model.theta",
            _population_model(),
            INTERIM_PRIOR,
            walkers=8,
            steps=2000,
            burn_in=500,
            seed=SEED,
        )

    def test_the_central_95_percent_interval_contains_the_truth(self, result: Any) -> None:
        mu = np.asarray(result["posterior"]["mu"]).ravel()
        tau = np.asarray(result["posterior"]["tau"]).ravel()
        mu_lo, mu_hi = np.percentile(mu, [2.5, 97.5])
        tau_lo, tau_hi = np.percentile(tau, [2.5, 97.5])
        print(
            f"\nreweighted population posterior: mu 95% = ({mu_lo:.4f}, {mu_hi:.4f}) "
            f"[truth {MU_TRUTH}], tau 95% = ({tau_lo:.4f}, {tau_hi:.4f}) [truth {TAU_TRUTH}]"
        )
        assert mu_lo <= MU_TRUTH <= mu_hi, (mu_lo, MU_TRUTH, mu_hi)
        assert tau_lo <= TAU_TRUTH <= tau_hi, (tau_lo, TAU_TRUTH, tau_hi)

    def test_mu_is_narrower_than_any_single_objects_own_posterior(
        self, result: Any, object_runs: list[Any]
    ) -> None:
        mu = np.asarray(result["posterior"]["mu"]).ravel()
        population_std = float(mu.std())
        object_stds = [
            float(np.asarray(run["posterior"]["model.theta"]).std()) for run in object_runs
        ]
        assert population_std < min(object_stds), (population_std, min(object_stds))

    def test_root_attrs_name_the_derived_product(self, result: Any) -> None:
        import json

        assert result.attrs["ampere_population_parameter"] == "model.theta"
        runs_attr = json.loads(result.attrs["ampere_population_runs"])
        assert len(runs_attr) == N_OBJECTS_REDUCED
        assert all("problem_hash" in row and "spec_hash" in row for row in runs_attr)
        ess_attr = json.loads(result.attrs["ampere_population_ess"])
        assert len(ess_attr["per_object"]) == N_OBJECTS_REDUCED
        assert ess_attr["min"] <= ess_attr["median"]
        model_attr = json.loads(result.attrs["ampere_population_model"])
        assert model_attr["kind"] == "GaussianPopulationModel"
        assert model_attr["hyperparameters"] == ["mu", "tau"]


# ---------------------------------------------------------------------------
# 8. The stored interim prior (W5.22): optional, verified, refused when absent
# ---------------------------------------------------------------------------


@pytest.mark.study
class TestTheStoredInterimPrior:
    """``fit_population`` reads ``model.theta``'s prior back off the archive's
    own provenance (schema 8, ``ampere_free_priors``) when ``interim_prior``
    is omitted, on the same :data:`N_OBJECTS_REDUCED`-object archive
    :class:`TestReweightedPopulationPosterior` uses -- module docstring
    claim 8.
    """

    SETTINGS: ClassVar[dict[str, int]] = dict(walkers=8, steps=2000, burn_in=500, seed=SEED)

    @pytest.fixture(scope="class")
    def columns(self, object_runs: list[Any]) -> list[Any]:
        return [DataTreeRunColumns(run) for run in object_runs]

    @pytest.fixture(scope="class")
    def without_prior(self, columns: list[Any]) -> Any:
        """One ``fit_population`` call with no ``interim_prior``, shared by
        both tests below -- the same settings and archive
        :class:`TestReweightedPopulationPosterior`'s own ``result`` fixture
        uses, so a second, independent population fit is not paid for
        twice."""
        return fit_population(columns, "model.theta", _population_model(), **self.SETTINGS)

    def test_omitting_interim_prior_reproduces_the_supplied_result(
        self, columns: list[Any], without_prior: Any
    ) -> None:
        """No ``interim_prior`` reproduces the explicit-prior row bitwise.

        ``ConstantModel``'s declared prior (every object in ``object_runs``)
        and :data:`INTERIM_PRIOR` are the same ``norm(0.0, 10.0)``, so the
        stored prior :func:`~ampere.results.fit_population` reads back and
        the one :class:`TestReweightedPopulationPosterior` supplies
        explicitly are one declaration -- the two runs are the same
        ``emcee`` ensemble at the same seed and must agree exactly, not just
        within a tolerance.
        """
        with_prior = fit_population(
            columns, "model.theta", _population_model(), INTERIM_PRIOR, **self.SETTINGS
        )
        for name in ("mu", "tau"):
            np.testing.assert_array_equal(
                np.asarray(with_prior["posterior"][name]),
                np.asarray(without_prior["posterior"][name]),
            )

    def test_omitting_interim_prior_recovers_the_truth(self, without_prior: Any) -> None:
        """The same central-95 % claim :class:`TestReweightedPopulationPosterior` pins,
        with the prior read back from provenance instead of supplied."""
        mu = np.asarray(without_prior["posterior"]["mu"]).ravel()
        tau = np.asarray(without_prior["posterior"]["tau"]).ravel()
        mu_lo, mu_hi = np.percentile(mu, [2.5, 97.5])
        tau_lo, tau_hi = np.percentile(tau, [2.5, 97.5])
        assert mu_lo <= MU_TRUTH <= mu_hi, (mu_lo, MU_TRUTH, mu_hi)
        assert tau_lo <= TAU_TRUTH <= tau_hi, (tau_lo, TAU_TRUTH, tau_hi)

    def test_a_disagreeing_supplied_prior_is_refused(self, columns: list[Any]) -> None:
        disagreeing = st.norm(0.0, 20.0)  # every object's declared prior is norm(0.0, 10.0)
        with pytest.raises(ResultsError, match="disagrees"):
            fit_population(
                columns, "model.theta", _population_model(), disagreeing, **self.SETTINGS
            )

    def test_no_supplied_and_no_stored_prior_is_refused(self, object_runs: list[Any]) -> None:
        """A schema-7-style archive (no ``ampere_free_priors``) with no
        ``interim_prior`` supplied either: nothing here to read back and
        nothing to fall back on, refused by name (W3.12's append-refusal
        precedent for a file that predates the attribute it needs)."""
        pre_schema_8 = []
        for run in object_runs[:5]:
            columns_run = DataTreeRunColumns(run)
            attrs = dict(columns_run.attrs)
            del attrs[f"{ATTR_PREFIX}free_priors"]
            pre_schema_8.append(
                ManualRunColumns(
                    theta=columns_run.parameter_draws("model.theta"),
                    log_prior=columns_run.log_prior,
                    log_likelihood=columns_run.log_likelihood,
                    proposal_log_density=columns_run.proposal_log_density,
                    attrs=attrs,
                )
            )
        with pytest.raises(ResultsError, match="no stored prior"):
            fit_population(
                pre_schema_8,
                "model.theta",
                _population_model(),
                walkers=8,
                steps=100,
                burn_in=20,
                seed=SEED,
            )


# ---------------------------------------------------------------------------
# 2. The netCDF-directory reader agrees bitwise with the in-memory one
# ---------------------------------------------------------------------------


@pytest.mark.study
class TestTheFileBackedReaderAgreesWithTheInMemoryOne:
    def test_bitwise_equal_result(self, object_runs: list[Any], tmp_path: Any) -> None:
        model = _population_model()
        settings = dict(walkers=8, steps=200, burn_in=50, seed=SEED)

        in_memory = [DataTreeRunColumns(run) for run in object_runs]
        result_memory = fit_population(in_memory, "model.theta", model, INTERIM_PRIOR, **settings)

        directory = tmp_path / "archived_runs"
        directory.mkdir()
        for index, run in enumerate(object_runs):
            to_netcdf(run, directory / f"object_{index:03d}.nc")
        file_backed = runs_from_netcdf_directory(directory)
        result_file = fit_population(file_backed, "model.theta", model, INTERIM_PRIOR, **settings)

        for name in ("mu", "tau"):
            np.testing.assert_array_equal(
                np.asarray(result_memory["posterior"][name]),
                np.asarray(result_file["posterior"][name]),
            )
        np.testing.assert_array_equal(
            np.asarray(result_memory["sample_stats"]["lp"]),
            np.asarray(result_file["sample_stats"]["lp"]),
        )
        assert (
            result_memory.attrs["ampere_population_runs"]
            == result_file.attrs["ampere_population_runs"]
        )

    def test_round_trip_reads_back_the_same_columns(
        self, object_runs: list[Any], tmp_path: Any
    ) -> None:
        run = object_runs[0]
        path = tmp_path / "one_run.nc"
        to_netcdf(run, path)
        reread = from_netcdf(path)
        in_memory = DataTreeRunColumns(run)
        from_disk = DataTreeRunColumns(reread)
        np.testing.assert_array_equal(
            in_memory.parameter_draws("model.theta"), from_disk.parameter_draws("model.theta")
        )
        np.testing.assert_array_equal(in_memory.log_prior, from_disk.log_prior)
        np.testing.assert_array_equal(in_memory.log_likelihood, from_disk.log_likelihood)


# ---------------------------------------------------------------------------
# 3. A minimal RunColumns implementation that is not a DataTree
# ---------------------------------------------------------------------------


class ManualRunColumns:
    """A bare structural :class:`RunColumns` -- no ``xarray`` anywhere in it.

    Exists to prove horizon (b)'s "buildable entirely on stored files" claim
    is not secretly "buildable entirely on ``xarray.DataTree``": this class
    never constructs one, and :func:`~ampere.results.fit_population` cannot
    tell the difference.
    """

    def __init__(
        self,
        theta: np.ndarray,
        log_prior: np.ndarray,
        log_likelihood: np.ndarray,
        proposal_log_density: np.ndarray | None,
        attrs: dict[str, Any],
    ) -> None:
        self._theta = theta
        self._log_prior = log_prior
        self._log_likelihood = log_likelihood
        self._proposal_log_density = proposal_log_density
        self._attrs = attrs

    @property
    def attrs(self) -> dict[str, Any]:
        return self._attrs

    def parameter_draws(self, name: str) -> np.ndarray:
        assert name == "model.theta"
        return self._theta

    @property
    def log_prior(self) -> np.ndarray:
        return self._log_prior

    @property
    def log_likelihood(self) -> np.ndarray:
        return self._log_likelihood

    @property
    def proposal_log_density(self) -> np.ndarray | None:
        return self._proposal_log_density


def test_the_protocol_accepts_a_non_datatree_implementation() -> None:
    manual = ManualRunColumns(
        theta=np.array([1.0, 1.1]),
        log_prior=np.array([-1.0, -1.1]),
        log_likelihood=np.array([-2.0, -2.1]),
        proposal_log_density=None,
        attrs={"ampere_spec_hash": "abc"},
    )
    assert isinstance(manual, RunColumns)
    assert not isinstance(object(), RunColumns)


# ---------------------------------------------------------------------------
# 4. The ESS refusal
# ---------------------------------------------------------------------------


@pytest.mark.study
class TestTheEffectiveSampleSizeRefusal:
    def test_a_collapsed_object_is_refused_by_name(self, object_runs: list[Any]) -> None:
        good = [DataTreeRunColumns(run) for run in object_runs[:5]]
        spec_hash = good[0].attrs["ampere_spec_hash"]
        n_draws = good[0].log_prior.size

        # 99% of the draws sit miles from the population; one lands exactly on
        # the truth. Whatever alpha the other five objects pull the fit
        # towards, this object's weight collapses onto that one draw.
        theta = np.full(n_draws, 100.0)
        theta[0] = MU_TRUTH
        interim = st.norm(0.0, 10.0)
        adversarial = ManualRunColumns(
            theta=theta,
            log_prior=interim.logpdf(theta),
            log_likelihood=np.zeros(n_draws),
            proposal_log_density=None,
            attrs={"ampere_spec_hash": spec_hash, "ampere_approximation": "none"},
        )

        with pytest.raises(ResultsError, match="effective sample size"):
            fit_population(
                [*good, adversarial],
                "model.theta",
                _population_model(),
                INTERIM_PRIOR,
                walkers=8,
                steps=300,
                burn_in=100,
                ess_floor=20.0,
                seed=SEED,
            )


# ---------------------------------------------------------------------------
# 5. The mixed-spec-hash refusal
# ---------------------------------------------------------------------------


@pytest.mark.study
def test_runs_that_do_not_share_a_spec_hash_are_refused() -> None:
    run_a = _fit_object(1.0, seed=1, prior=st.norm(0.0, 10.0))
    run_b = _fit_object(1.0, seed=2, prior=st.norm(0.0, 20.0))  # a different declared prior

    assert (
        DataTreeRunColumns(run_a).attrs["ampere_spec_hash"]
        != DataTreeRunColumns(run_b).attrs["ampere_spec_hash"]
    )

    with pytest.raises(ResultsError, match="spec_hash"):
        fit_population(
            [DataTreeRunColumns(run_a), DataTreeRunColumns(run_b)],
            "model.theta",
            _population_model(),
            INTERIM_PRIOR,
            walkers=8,
            steps=100,
            burn_in=20,
            seed=1,
        )


# ---------------------------------------------------------------------------
# 6. An object fitted by VIEngine, reweighted through proposal_log_density
# ---------------------------------------------------------------------------


def _vi_kit() -> tuple[str, Any] | tuple[None, None]:
    try:
        from ampere.inference._vi import VARIATIONAL_LIBRARIES
    except ImportError:
        return None, None
    for name in sorted(VARIATIONAL_LIBRARIES):
        try:
            module = importlib.import_module(f"ampere.backends.{name}")
        except ImportError:
            continue
        if name == "jax":
            module.configure_x64()
        return name, module
    return None, None


_VI_BACKEND_NAME, _VI_BACKEND_MODULE = _vi_kit()

requires_vi = pytest.mark.skipif(
    _VI_BACKEND_MODULE is None,
    reason="no backend with a variational library installed (sbi environment)",
)

_VI_GRID = np.geomspace(1.0, 10.0, 20)
_VI_SIGMA = 0.2
_VI_INDEX = -1.0
_VI_REFERENCE_WAVELENGTH = 1.0
_VI_TRUE_NORM = 2.0
_VI_INTERIM_PRIOR = st.norm(2.0, 0.5)  # PowerLaw's declared norm prior; index is fixed, not free


def _vi_observed(seed: int = 7) -> Spectrum:
    rng = np.random.default_rng(seed)
    truth = _VI_TRUE_NORM * (_VI_GRID / _VI_REFERENCE_WAVELENGTH) ** _VI_INDEX
    return Spectrum(
        _VI_GRID * u.um,
        (truth + rng.normal(0.0, _VI_SIGMA, _VI_GRID.size)) * u.Jy,
        uncertainty=np.full(_VI_GRID.size, _VI_SIGMA) * u.Jy,
    )


def _vi_problem(seed: int) -> FittingProblem:
    module = _VI_BACKEND_MODULE
    likelihood = Likelihood(GaussianFamily(), module.IndependentNoise())
    return FittingProblem(
        module.PowerLaw(
            _VI_GRID,
            norm=st.norm(2.0, 0.5),
            index=_VI_INDEX,
            reference_wavelength=_VI_REFERENCE_WAVELENGTH,
        ),
        [Dataset(_vi_observed(), likelihood=likelihood)],
        seed=seed,
    )


@requires_vi
class TestTheApproximateEngineRow:
    def test_vi_fitted_object_agrees_with_emcee(self) -> None:
        from ampere.inference import VIEngine

        emcee_run = EmceeEngine(_vi_problem(seed=SEED), walkers=16).run(steps=800, burn_in=200)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            vi_run = VIEngine(_vi_problem(seed=SEED)).run(draws=4000, steps=3000, guide="normal")

        assert vi_run.attrs["ampere_approximation"] != "none"
        vi_columns = DataTreeRunColumns(vi_run)
        assert vi_columns.proposal_log_density is not None

        # A direct check of the W5.0 importance weight, before it ever reaches
        # a population fit: the self-normalised weighted mean of the VI
        # object's draws should land near the emcee posterior's own mean.
        weights = np.exp(_self_normalised_log_weights(vi_columns, vi_columns.log_prior))
        weighted_mean = float(np.sum(weights * vi_columns.parameter_draws("model.norm")))
        emcee_mean = float(np.asarray(emcee_run["posterior"]["model.norm"]).mean())
        assert weighted_mean == pytest.approx(emcee_mean, abs=0.1)

        model = GaussianPopulationModel(
            Parameter("mu", st.norm(2.0, 0.5)),
            Parameter("tau", st.uniform(0.005, 0.1)),
        )
        settings = dict(walkers=8, steps=800, burn_in=200, seed=SEED)
        result_emcee = fit_population(
            [DataTreeRunColumns(emcee_run)], "model.norm", model, _VI_INTERIM_PRIOR, **settings
        )
        result_vi = fit_population([vi_columns], "model.norm", model, _VI_INTERIM_PRIOR, **settings)

        mu_emcee = float(np.asarray(result_emcee["posterior"]["mu"]).mean())
        mu_vi = float(np.asarray(result_vi["posterior"]["mu"]).mean())
        assert mu_vi == pytest.approx(mu_emcee, abs=0.1)


# ---------------------------------------------------------------------------
# 7. The marginal-vs-joint interim-prior fix: a two-parameter regression
# ---------------------------------------------------------------------------


@pytest.mark.study
class TestTheMarginalInterimPriorFix:
    """A run's stored ``log_prior`` is the *joint* interim prior over every
    free parameter; the reweighting ratio needs the named parameter's own
    *marginal* interim prior alone (module docstring, "why it needs the
    marginal interim prior"). ``ConstantModel`` has one free parameter, so
    the joint and the marginal coincide and the earlier version of this
    module -- which divided by the joint column -- passed every row above
    without ever exercising the bug. ``TwoParameterModel`` adds ``phi``, a
    pure nuisance with a narrow, non-flat prior that never enters the
    likelihood, so its marginal log-density varies materially draw to draw
    -- exactly the factor a joint-column denominator leaves uncancelled.

    Reuses :data:`N_OBJECTS_REDUCED` of :func:`_truths_and_data`'s 200
    (truth, datum) pairs -- the same realisation :data:`object_runs`
    (``ConstantModel``) already fitted, and since W5.22 that fixture *is*
    exactly this size rather than a larger archive sliced down -- so the
    one-parameter control costs no extra fitting, and only the
    :data:`N_OBJECTS_REDUCED` two-parameter fits are new. Kept well under
    ten minutes even alongside a concurrent gate.
    """

    N = N_OBJECTS_REDUCED
    POPULATION_SETTINGS: ClassVar[dict[str, int]] = dict(
        walkers=8, steps=1500, burn_in=400, seed=SEED
    )

    @pytest.fixture(scope="class")
    def two_parameter_runs(self) -> list[Any]:
        _, data = _truths_and_data(SEED, N_OBJECTS)
        return [
            DataTreeRunColumns(_fit_two_parameter_object(float(datum), seed=SEED + index))
            for index, datum in enumerate(data[: self.N])
        ]

    @pytest.fixture(scope="class")
    def one_parameter_result(self, object_runs: list[Any]) -> Any:
        """The one-parameter control, over the *same* N_OBJECTS_REDUCED objects."""
        columns = [DataTreeRunColumns(run) for run in object_runs[: self.N]]
        return fit_population(
            columns, "model.theta", _population_model(), INTERIM_PRIOR, **self.POPULATION_SETTINGS
        )

    @pytest.fixture(scope="class")
    def two_parameter_result(self, two_parameter_runs: list[Any]) -> Any:
        return fit_population(
            two_parameter_runs,
            "model.theta",
            _population_model(),
            INTERIM_PRIOR,
            **self.POPULATION_SETTINGS,
        )

    def test_truth_inside_the_central_95_percent(self, two_parameter_result: Any) -> None:
        mu = np.asarray(two_parameter_result["posterior"]["mu"]).ravel()
        tau = np.asarray(two_parameter_result["posterior"]["tau"]).ravel()
        mu_lo, mu_hi = np.percentile(mu, [2.5, 97.5])
        tau_lo, tau_hi = np.percentile(tau, [2.5, 97.5])
        print(
            f"\ntwo-parameter reweighted posterior: mu 95% = ({mu_lo:.4f}, {mu_hi:.4f}) "
            f"[truth {MU_TRUTH}], tau 95% = ({tau_lo:.4f}, {tau_hi:.4f}) [truth {TAU_TRUTH}]"
        )
        assert mu_lo <= MU_TRUTH <= mu_hi, (mu_lo, MU_TRUTH, mu_hi)
        assert tau_lo <= TAU_TRUTH <= tau_hi, (tau_lo, TAU_TRUTH, tau_hi)

    def test_agrees_with_the_one_parameter_case(
        self, one_parameter_result: Any, two_parameter_result: Any
    ) -> None:
        """The nuisance parameter is irrelevant to the likelihood, so a
        correct reweighting must recover essentially the same (mu, tau) as
        the one-parameter control on the same 50 objects. The tolerance is
        generous (a fraction of the 50-object posterior's own width) because
        this compares two independent short ``emcee`` runs on two
        independently-sampled per-object fits, not a bitwise check.
        """
        mu_one = float(np.asarray(one_parameter_result["posterior"]["mu"]).mean())
        mu_two = float(np.asarray(two_parameter_result["posterior"]["mu"]).mean())
        tau_one = float(np.asarray(one_parameter_result["posterior"]["tau"]).mean())
        tau_two = float(np.asarray(two_parameter_result["posterior"]["tau"]).mean())
        print(
            f"\none-parameter control: mu={mu_one:.4f} tau={tau_one:.4f}; "
            f"two-parameter: mu={mu_two:.4f} tau={tau_two:.4f}"
        )
        assert mu_two == pytest.approx(mu_one, abs=0.15)
        assert tau_two == pytest.approx(tau_one, abs=0.1)
