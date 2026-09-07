"""The NUTS driver, end to end, on **every** backend that registers a realisation.

``tests/inference/test_engines.py`` holds the three gradient-free drivers to
five things; this file holds the fourth driver to the same ones, minus what is
structurally shared and already covered there (the import-graph check, the
netCDF round trip, the failure path). What is *new* is what NUTS is for:

1. **It recovers the toy joint problem's posterior.** ``inference.md`` §15's
   shape — one power law on two channels, two genuinely different instrument
   chains, the calibration factors tied — built entirely from one backend's
   own pieces, sampled through a gradient, and held to the truth the data were
   generated at.
2. **It agrees with a posterior written down in closed form.** The agreement
   problem is deliberately conjugate: a power law with its index held fixed is
   *linear* in ``norm``, so a Gaussian prior and Gaussian noise give a Gaussian
   posterior whose mean and variance are arithmetic rather than a long
   reference run. Comparing a sampler only against other samplers would pass
   several samplers that are wrong in the same way.
3. **It refuses, by name, what it cannot sample.** A problem on a backend with
   no realisation, a problem that declares itself non-differentiable, and — the
   one that would otherwise be silent — a density that disagrees with the
   problem it was handed.
4. **Every run emits the run**, through the same
   :meth:`~ampere.inference.engine.Engine.finish`, with ``ampere_backend`` read
   off the problem rather than declared by the driver, and (W2.13)
   ``ampere_realised`` saying the draws came through a realisation.

**Parametrised over the backends installed here** (W2.13). ``inference.md``
§10a makes the differentiable form of a problem a registered object, and this
driver dispatches between numpyro and pyro on the problem's backend — so the
claim "NUTS recovers the toy joint posterior" is one claim per backend, and
writing it once per backend by hand would let the two drift. Each row runs on
whichever of jax and torch this environment has; both, in the ``torch`` and
``jax`` pixi environments, is not possible today (they are separate extras)
and is not assumed anywhere here.

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

import dataclasses
import importlib
import math
import warnings
from typing import Any

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.backends.reference import PowerLaw as ReferencePowerLaw
from ampere.core import (
    Dataset,
    DatasetCollection,
    FittingProblem,
    GaussianFamily,
    Instrument,
    Likelihood,
    Spectrum,
    Tie,
)
from ampere.inference import EngineError, NUTSEngine

REFERENCE_WAVELENGTH = 1.0
TRUTH = {"norm": 2.0, "index": -1.2, "calibration": 1.0}
SEED = 20260907

FINE = np.geomspace(1.0, 20.0, 40)
COARSE = np.geomspace(1.5, 15.0, 12)


# ---------------------------------------------------------------------------
# The backends this environment can run
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class Kit:
    """One backend's pieces, by name, so a test body never names a library.

    Everything a problem here is built from — the model, the two instrument
    steps, the noise model — comes off this object, which is exactly the set
    ``inference.md`` §18 says a backend supplies. The rows below are then the
    same rows for every backend, which is what makes "NUTS works on this
    backend" a claim the suite can check rather than one each track asserts
    about itself.
    """

    name: str
    module: Any

    def likelihood(self) -> Likelihood:
        """Gaussian, with **this backend's** uncorrelated noise.

        Not ``Dataset``'s default: since W2.13 a noise model carries the four
        capability flags and ``ampere.core.IndependentNoise`` declares
        ``"reference"``, so the default would make every problem here a
        two-backend problem and be refused at composition.
        """
        return Likelihood(GaussianFamily(), self.module.IndependentNoise())


def _installed_kits() -> list[Kit]:
    found: list[Kit] = []
    for name in ("jax", "torch"):
        try:
            module = importlib.import_module(f"ampere.backends.{name}")
        except ImportError:  # the extra is not installed in this environment
            continue
        if name == "jax":
            # ``lowering.md`` §10.2(a): the application turns the flag on, not
            # ampere. A test suite is an application.
            module.configure_x64()
        found.append(Kit(name, module))
    return found


KITS = _installed_kits()

pytestmark = pytest.mark.skipif(
    not KITS,
    reason="no differentiable backend installed; NUTS needs a registered realisation",
)


@pytest.fixture(scope="module", params=[kit.name for kit in KITS])
def kit(request: Any) -> Kit:
    return next(found for found in KITS if found.name == request.param)


# ---------------------------------------------------------------------------
# The data (backend-neutral: containers are ampere.core's on every backend)
# ---------------------------------------------------------------------------


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


def joint_problem(kit: Kit, seed: int | None = SEED) -> FittingProblem:
    """``inference.md`` §15's shape, on *kit*'s backend.

    Deliberately the same declaration ``test_engines.py``'s ``joint_problem``
    builds on the reference backend, so that "NUTS recovers the toy joint
    problem's posterior" means the same problem the other three drivers are
    held to. The tie costs a dimension: three free parameters, not four.
    """
    model = kit.module.PowerLaw(
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
                        [kit.module.CalibrationScale(st.lognorm(0.05), label="calibration")],
                        channel="blue",
                    ),
                    kit.likelihood(),
                ),
                "red": Dataset(
                    RED_DATA,
                    Instrument(
                        [
                            kit.module.Resample(COARSE),
                            kit.module.CalibrationScale(st.lognorm(0.05), label="calibration"),
                        ],
                        channel="red",
                    ),
                    kit.likelihood(),
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


def agreement_problem(kit: Kit, seed: int | None = SEED) -> FittingProblem:
    """One dataset, one free parameter, an exactly known posterior."""
    return FittingProblem(
        kit.module.PowerLaw(
            AGREEMENT_GRID,
            norm=st.norm(*AGREEMENT_PRIOR),
            index=AGREEMENT_INDEX,
            reference_wavelength=REFERENCE_WAVELENGTH,
        ),
        [Dataset(AGREEMENT_DATA, likelihood=kit.likelihood())],
        seed=seed,
    )


def sample(problem: FittingProblem, kit: Kit, **settings: Any) -> Any:
    """Build the driver from the backend's lowering, explicitly, and run it.

    The explicit route, which is the one that exercises the ``density=``
    argument and the driver's own agreement check. ``NUTSEngine(problem)`` —
    the registered realisation — is what the rest of the suite uses and what
    ``TestTheRegisteredRealisation`` compares this against.
    """
    engine = NUTSEngine(problem, kit.module.lower_problem(problem).log_prob_unconstrained)
    with warnings.catch_warnings():
        # A short chain on a three-dimensional posterior occasionally records a
        # failure at a start point; the failure path itself is test_engines.py's.
        warnings.simplefilter("ignore")
        return engine.run(**settings)


def realised_sample(problem: FittingProblem, **settings: Any) -> Any:
    """``NUTSEngine(problem)`` with no density argument: §10a's whole point."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return NUTSEngine(problem).run(**settings)


@pytest.fixture(scope="module")
def joint_run(kit: Kit) -> Any:
    return realised_sample(joint_problem(kit), draws=400, warmup=400, chains=2)


@pytest.fixture(scope="module")
def agreement_run(kit: Kit) -> Any:
    return realised_sample(agreement_problem(kit), draws=600, warmup=400, chains=2)


# ---------------------------------------------------------------------------
# 1. The toy joint problem samples end to end
# ---------------------------------------------------------------------------


class TestTheJointProblemSamples:
    def test_the_problem_is_the_one_the_contract_describes(self, kit: Kit) -> None:
        """Three free dimensions, not four: the tie costs one (``inference.md`` §15)."""
        problem = joint_problem(kit)
        assert problem.free_size == 3
        assert problem.shared_names == ("calibration",)
        assert problem.backend == kit.name
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

    def test_the_backend_is_derived_not_declared(self, joint_run: Any, kit: Kit) -> None:
        """W2.12: no driver takes a ``backend=``; the flag is read off the problem."""
        assert joint_run.attrs["ampere_engine"] == "nuts"
        assert joint_run.attrs["ampere_backend"] == kit.name

    def test_the_run_records_that_it_was_realised(self, joint_run: Any) -> None:
        """W2.13, ``inference.md`` §10a's "Provenance" paragraph."""
        assert joint_run.attrs["ampere_realised"] == 1
        assert "ampere_registered_lowerings" in joint_run.attrs

    def test_the_sampler_that_drove_it_is_recorded(self, joint_run: Any, kit: Kit) -> None:
        from ampere.inference._nuts import SAMPLER_LIBRARIES

        assert joint_run.attrs["ampere_nuts_sampler"] == SAMPLER_LIBRARIES[kit.name]

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

    def test_the_fixed_index_takes_no_sampler_dimension(self, kit: Kit) -> None:
        problem = agreement_problem(kit)
        assert problem.free_size == 1
        assert problem.parameters.fixed_names == ("model.index",)


# ---------------------------------------------------------------------------
# 3. The registered realisation, and the refusals
# ---------------------------------------------------------------------------


class TestTheRegisteredRealisation:
    """``inference.md`` §10a: the density comes through ``ampere.core.realise``.

    Importing the backend registered its realisation; the driver asks core for
    it and never imports the backend itself.
    """

    def test_importing_the_backend_registered_its_realisation(self, kit: Kit) -> None:
        from ampere.core import registered_realisations

        assert registered_realisations().get(kit.name) is True

    def test_the_driver_needs_no_density_argument(self, kit: Kit) -> None:
        problem = agreement_problem(kit)
        engine = NUTSEngine(problem)
        reference = problem.unconstrain(problem.reference_values)
        assert float(np.asarray(engine.density(reference))) == pytest.approx(
            float(problem.log_prob_unconstrained(reference))
        )

    def test_the_realisation_is_kept_so_provenance_can_read_it(self, kit: Kit) -> None:
        engine = NUTSEngine(agreement_problem(kit))
        assert engine.realisation is not None
        assert engine.realisation.backend == kit.name

    def test_an_explicit_density_is_not_recorded_as_realised(self, kit: Kit) -> None:
        """``ampere_realised`` is a fact about the run, not about the backend."""
        problem = agreement_problem(kit)
        engine = NUTSEngine(problem, kit.module.lower_problem(problem).log_prob_unconstrained)
        assert engine.realisation is None
        run = sample(agreement_problem(kit), kit, draws=20, warmup=20, chains=1)
        assert run.attrs["ampere_realised"] == 0

    def test_the_registered_route_samples_the_same_posterior(self, kit: Kit) -> None:
        # Two problems with one seed, as TestReproducibility does: every
        # stream comes from the problem's seed, so the only difference between
        # the two runs is where the density came from -- and there is none.
        explicit = sample(agreement_problem(kit, SEED), kit, draws=30, warmup=30, chains=1)
        registered = realised_sample(agreement_problem(kit, SEED), draws=30, warmup=30, chains=1)
        assert np.asarray(registered["posterior"]["model.norm"]) == pytest.approx(
            np.asarray(explicit["posterior"]["model.norm"]), abs=0.0
        )

    def test_a_backend_without_a_realisation_is_refused_by_name(self, kit: Kit) -> None:
        problem = agreement_problem(kit)
        from ampere.core import realisation as registry

        saved = dict(registry._REALISATIONS)
        try:
            registry._REALISATIONS.clear()
            with pytest.raises(EngineError, match="cannot sample a problem"):
                NUTSEngine(problem)
        finally:
            registry._REALISATIONS.clear()
            registry._REALISATIONS.update(saved)

    def test_supported_backends_is_answered_from_the_registry(self, kit: Kit) -> None:
        from ampere.inference._nuts import supported_backends

        assert kit.name in supported_backends()
        assert "reference" not in supported_backends()


class TestRefusals:
    def test_a_problem_on_another_backend_is_refused_by_name(self, kit: Kit) -> None:
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

    def test_a_density_that_is_not_callable_is_refused(self, kit: Kit) -> None:
        with pytest.raises(EngineError, match="callable"):
            NUTSEngine(agreement_problem(kit), 3.0)  # type: ignore[arg-type]

    def test_a_density_for_a_different_problem_is_refused(self, kit: Kit) -> None:
        """The refusal that would otherwise be silent.

        Sampling one problem's density while recording another's produces a run
        that passes every structural check and describes a posterior nobody
        asked for. One point is enough to catch it and costs nothing.
        """
        problem = agreement_problem(kit)
        other = kit.module.lower_problem(joint_problem(kit))
        with pytest.raises(EngineError, match="disagrees with the problem"):
            NUTSEngine(problem, lambda y: other.log_prob_unconstrained(np.zeros(3)))

    def test_a_problem_with_nothing_to_sample_is_refused(self, kit: Kit) -> None:
        problem = FittingProblem(
            kit.module.PowerLaw(
                AGREEMENT_GRID, norm=2.0, index=-1.0, reference_wavelength=REFERENCE_WAVELENGTH
            ),
            [Dataset(AGREEMENT_DATA, likelihood=kit.likelihood())],
        )
        with pytest.raises(EngineError, match="nothing to sample"):
            NUTSEngine(problem, lambda y: 0.0)


# ---------------------------------------------------------------------------
# 4. Reproducibility
# ---------------------------------------------------------------------------


class TestReproducibility:
    def test_the_same_seed_gives_the_same_draws(self, kit: Kit) -> None:
        """``inference.md`` §12: every stream comes from the problem's own seed.

        Initialisation and the sampler's own randomness are separate labels
        under this engine's name, so neither can make the other
        irreproducible.
        """
        first = realised_sample(joint_problem(kit, SEED), draws=30, warmup=30, chains=1)
        second = realised_sample(joint_problem(kit, SEED), draws=30, warmup=30, chains=1)
        assert np.asarray(first["posterior"]["model.norm"]) == pytest.approx(
            np.asarray(second["posterior"]["model.norm"]), abs=0.0
        )

    def test_a_different_seed_gives_different_draws(self, kit: Kit) -> None:
        first = realised_sample(joint_problem(kit, SEED), draws=30, warmup=30, chains=1)
        other = realised_sample(joint_problem(kit, SEED + 1), draws=30, warmup=30, chains=1)
        assert not np.allclose(
            np.asarray(first["posterior"]["model.norm"]),
            np.asarray(other["posterior"]["model.norm"]),
        )
