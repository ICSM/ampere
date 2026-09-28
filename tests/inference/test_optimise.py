"""``ampere.inference.optimise`` and ``warm_start_gp`` — point estimates and starts (W6.7).

The pinned rows of the item's ruling 4, on the numpy path here and on torch
and jax where those are installed (each native row skips by ``find_spec`` and
runs in ``pixi run -e torch`` / ``-e jax``):

(a) each route's estimate inside the sampled posterior's central 50 % on every
    free parameter — exact to ``1e-3`` on the conjugate ``agreement_problem``
    (MAP = posterior mean there), and against an emcee posterior on
    ``examples/sed_composition``;
(f) the refusals by name, and a covariance refusal on a saddle.

The bridge, provenance, the native routes, the warm start and the NUTS/emcee
rows join as their units land.
"""

from __future__ import annotations

import math
from importlib.util import find_spec
from typing import Any

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.backends.reference import PowerLaw
from ampere.core import Dataset, FittingProblem, Model, Parameter, Spectrum
from ampere.inference import EmceeEngine, EngineError, optimise
from ampere.inference._optimise import (
    constrained_objective,
    covariance_from_hessian,
    finite_difference_hessian,
)
from ampere.results import Optimum

SEED = 20260929

# -- the conjugate problem, copied from tests/inference/test_engines.py -------

REFERENCE_WAVELENGTH = 1.0
AGREEMENT_GRID = np.geomspace(1.0, 10.0, 20)
AGREEMENT_SIGMA = 0.2
AGREEMENT_INDEX = -1.0
AGREEMENT_PRIOR = (2.0, 0.5)


def _agreement_data() -> Spectrum:
    rng = np.random.default_rng(7)
    values = 2.0 * (AGREEMENT_GRID / REFERENCE_WAVELENGTH) ** AGREEMENT_INDEX
    return Spectrum(
        AGREEMENT_GRID * u.micron,
        (values + rng.normal(0.0, AGREEMENT_SIGMA, values.size)) * u.Jy,
        uncertainty=np.full(values.size, AGREEMENT_SIGMA) * u.Jy,
    )


AGREEMENT_DATA = _agreement_data()


def analytic() -> tuple[float, float]:
    """The conjugate posterior's (mean, sd) on ``norm`` — also its mode."""
    x = (AGREEMENT_GRID / REFERENCE_WAVELENGTH) ** AGREEMENT_INDEX
    y = np.asarray(AGREEMENT_DATA.values)
    mu, tau = AGREEMENT_PRIOR
    precision = 1.0 / tau**2 + float(np.sum(x**2)) / AGREEMENT_SIGMA**2
    mean = (mu / tau**2 + float(np.sum(x * y)) / AGREEMENT_SIGMA**2) / precision
    return mean, 1.0 / math.sqrt(precision)


def agreement_problem(seed: int | None = SEED) -> FittingProblem:
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


# -- a saddle: the likelihood is flat in one direction ------------------------

GRID = np.array([1.0, 2.0, 3.0, 4.0])


class OnlySlope(Model):
    """``slope * x``; ``ignored`` is declared, has a flat prior, and is never used."""

    def __init__(self) -> None:
        self.register_buffer("wavelength", GRID, unit=u.um)
        self.register_parameter(Parameter("slope", st.norm(1.0, 1.0)))
        self.register_parameter(Parameter("ignored", st.uniform(0.0, 1.0)))

    def evaluate(self, **values: Any) -> Spectrum:
        ctx = self.context(values)
        return Spectrum(ctx["wavelength"] * u.um, ctx["slope"] * ctx["wavelength"] * u.Jy)


def flat_problem() -> FittingProblem:
    observed = Spectrum(GRID * u.um, 2.0 * GRID * u.Jy, uncertainty=np.full(4, 0.1) * u.Jy)
    return FittingProblem(OnlySlope(), [Dataset(observed)], seed=SEED)


# -- the four-parameter composition, sampled once -----------------------------

HAS_TORCH = find_spec("torch") is not None
HAS_JAX = find_spec("jax") is not None

#: emcee budget for the sampled posterior the routes are held to: the
#: example's own default (16 walkers, 450 steps, 180 burn-in), about 40 s on
#: dev — the sampled central 50 % is the yardstick, so it must be a posterior.
SED_WALKERS, SED_STEPS, SED_BURN_IN = 16, 450, 180


@pytest.fixture(scope="module")
def sed_posterior() -> dict[str, np.ndarray]:
    """``examples/sed_composition``'s reference problem, sampled once by emcee."""
    from examples.sed_composition.sed_composition import build_problem

    run = EmceeEngine(build_problem("reference"), walkers=SED_WALKERS).run(
        SED_STEPS, burn_in=SED_BURN_IN
    )
    posterior = run["posterior"].dataset
    return {name: np.asarray(posterior[name]).ravel() for name in posterior.data_vars}


def central_50(optimum: Optimum, posterior: dict[str, np.ndarray]) -> dict[str, bool]:
    """Whether each free parameter's estimate lies in the sampled central 50 %."""
    inside: dict[str, bool] = {}
    for name in optimum.free_names:
        lower, upper = np.percentile(posterior[name], [25.0, 75.0])
        inside[name] = bool(lower <= float(optimum.constrained[name]) <= upper)
    return inside


# ---------------------------------------------------------------------------
# The objective and its curvature
# ---------------------------------------------------------------------------


class TestTheObjective:
    def test_it_is_the_constrained_density_without_the_jacobian(self) -> None:
        problem = flat_problem()
        objective = constrained_objective(problem)
        u0 = np.array([0.3, -0.7])
        theta = problem.constrain(u0)
        assert objective(u0) == pytest.approx(problem.log_prob(theta))
        # the uniform prior's logit bijection makes the two conventions differ
        assert objective(u0) != pytest.approx(problem.log_prob_unconstrained(u0))

    def test_finite_difference_hessian_of_a_quadratic(self) -> None:
        matrix = np.array([[4.0, 1.0], [1.0, 9.0]])

        def quadratic(x: np.ndarray) -> float:
            return 0.5 * float(x @ matrix @ x)

        hessian = finite_difference_hessian(quadratic, np.array([0.2, -0.1]))
        np.testing.assert_allclose(hessian, matrix, rtol=1e-6)

    def test_an_indefinite_hessian_is_refused_by_name(self) -> None:
        covariance, refusal = covariance_from_hessian(np.diag([1.0, -1.0]))
        assert covariance is None
        assert refusal is not None and "not positive definite" in refusal


# ---------------------------------------------------------------------------
# (a) The scipy route
# ---------------------------------------------------------------------------


class TestTheScipyRoute:
    @pytest.mark.parametrize("minimiser", ["Powell", "L-BFGS-B"])
    def test_the_conjugate_mode_is_exact(self, minimiser: str) -> None:
        mean, sd = analytic()
        optimum = optimise(agreement_problem(), method="scipy", starts=3, minimiser=minimiser)
        assert isinstance(optimum, Optimum)
        assert optimum.route == "scipy" and optimum.backend == "reference"
        assert optimum.constrained["model.norm"] == pytest.approx(mean, abs=1e-3)
        assert optimum.covariance is not None
        # identity bijection on a real-support prior: the unconstrained
        # covariance is the posterior's own variance
        assert math.sqrt(optimum.covariance[0, 0]) == pytest.approx(sd, rel=1e-3)
        assert len(optimum.starts) == 3 and optimum.converged
        assert optimum.log_prob_constrained == pytest.approx(
            agreement_problem().log_prob(np.array([optimum.constrained["model.norm"]]))
        )

    def test_auto_on_the_reference_backend_is_scipy(self) -> None:
        assert optimise(agreement_problem(), starts=1).route == "scipy"

    def test_seeded_problems_optimise_reproducibly(self) -> None:
        first = optimise(agreement_problem(), method="scipy", starts=2)
        second = optimise(agreement_problem(), method="scipy", starts=2)
        assert first.identity == second.identity
        assert [s.start_hash for s in first.starts] == [s.start_hash for s in second.starts]
        other = optimise(agreement_problem(), method="scipy", starts=2, seed=1)
        assert [s.start_hash for s in other.starts] != [s.start_hash for s in first.starts]

    def test_inside_the_sampled_central_50_percent(
        self, sed_posterior: dict[str, np.ndarray]
    ) -> None:
        from examples.sed_composition.sed_composition import build_problem

        optimum = optimise(build_problem("reference"), method="scipy", starts=2)
        assert central_50(optimum, sed_posterior) == dict.fromkeys(optimum.free_names, True)

    def test_provenance_rides_on_the_optimum(self) -> None:
        optimum = optimise(agreement_problem(), method="scipy", starts=1)
        assert optimum.provenance["ampere_engine"] == "optimise.scipy"
        assert optimum.provenance["ampere_seed"] == SEED


# ---------------------------------------------------------------------------
# (f) Refusals by name
# ---------------------------------------------------------------------------


class TestRefusals:
    def test_map_on_the_reference_backend_names_scipy(self) -> None:
        with pytest.raises(EngineError, match="method='scipy'"):
            optimise(agreement_problem(), method="map")

    def test_vi_on_the_reference_backend_names_scipy(self) -> None:
        with pytest.raises(EngineError, match="method='scipy'"):
            optimise(agreement_problem(), method="vi")

    def test_an_unknown_method(self) -> None:
        with pytest.raises(EngineError, match="does not know the method 'bayesopt'"):
            optimise(agreement_problem(), method="bayesopt")

    def test_an_unknown_option(self) -> None:
        with pytest.raises(EngineError, match=r"does not take the option\(s\) \['steps'\]"):
            optimise(agreement_problem(), method="scipy", steps=10)

    def test_a_saddle_refuses_its_covariance(self) -> None:
        optimum = optimise(flat_problem(), method="scipy", starts=2)
        assert optimum.covariance is None
        assert optimum.covariance_refusal is not None
        assert "not positive definite" in optimum.covariance_refusal
        assert optimum.constrained["model.slope"] == pytest.approx(2.0, abs=1e-2)


# ---------------------------------------------------------------------------
# The bridge: around= and run(initial=optimum)
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def conjugate_optimum() -> Optimum:
    return optimise(agreement_problem(), method="scipy", starts=2)


class TestTheBridge:
    def test_the_ball_is_tighter_than_the_posterior(self, conjugate_optimum: Optimum) -> None:
        mean, sd = analytic()
        engine = EmceeEngine(agreement_problem(), walkers=8)
        positions = engine.initial_positions(400, around=conjugate_optimum)
        assert positions.shape == (400, 1)
        assert np.mean(positions) == pytest.approx(mean, abs=0.1 * sd)
        assert np.std(positions) == pytest.approx(0.5 * sd, rel=0.15)

    def test_without_around_nothing_changed(self) -> None:
        first = EmceeEngine(agreement_problem(), walkers=8).initial_positions(8)
        second = EmceeEngine(agreement_problem(), walkers=8).initial_positions(8)
        np.testing.assert_array_equal(first, second)
        assert np.std(first) > analytic()[1]  # prior draws: far wider than the posterior

    def test_a_refused_covariance_gives_the_diagonal_ball(self) -> None:
        problem = flat_problem()
        optimum = optimise(problem, method="scipy", starts=1)
        assert optimum.covariance is None
        positions = EmceeEngine(problem, walkers=8).initial_positions(200, around=optimum)
        spread = np.std(np.stack([problem.unconstrain(t) for t in positions]), axis=0)
        expected = 0.01 * np.abs(optimum.unconstrained) + 1e-3
        np.testing.assert_allclose(spread, expected, rtol=0.2)

    def test_mismatched_names_are_refused_by_name(self, conjugate_optimum: Optimum) -> None:
        engine = EmceeEngine(flat_problem(), walkers=8)
        with pytest.raises(EngineError, match=r"missing \['model.ignored', 'model.slope'\]"):
            engine.initial_positions(8, around=conjugate_optimum)

    @pytest.mark.parametrize("engine", ["emcee", "zeus"])
    def test_ensembles_run_from_an_optimum(self, conjugate_optimum: Optimum, engine: str) -> None:
        from ampere.inference import ZeusEngine

        mean, sd = analytic()
        factory = {"emcee": EmceeEngine, "zeus": ZeusEngine}[engine]
        run = factory(agreement_problem(), walkers=8).run(20, initial=conjugate_optimum)
        first = np.asarray(run["posterior"].dataset["model.norm"])[:, 0]
        assert np.all(np.abs(first - mean) < 3 * sd)

    def test_dynesty_refuses_an_optimum_by_name(self, conjugate_optimum: Optimum) -> None:
        from ampere.inference import DynestyEngine

        with pytest.raises(EngineError, match=r"nested sampler .* no start point"):
            DynestyEngine(agreement_problem(), live_points=30).run(initial=conjugate_optimum)

    @pytest.mark.parametrize("library", ["nautilus", "ultranest"])
    def test_the_other_nested_samplers_refuse_it(
        self, conjugate_optimum: Optimum, library: str
    ) -> None:
        if find_spec(library) is None:
            pytest.skip(f"{library} is not installed here")
        from ampere.inference import NautilusEngine, UltranestEngine

        factory = {"nautilus": NautilusEngine, "ultranest": UltranestEngine}[library]
        with pytest.raises(EngineError, match=r"nested sampler .* no start point"):
            factory(agreement_problem()).run(initial=conjugate_optimum)


# ---------------------------------------------------------------------------
# (e) The start in provenance (schema 9)
# ---------------------------------------------------------------------------


class TestTheStartInProvenance:
    def test_a_run_from_an_optimum_records_it(self, conjugate_optimum: Optimum) -> None:
        import json

        run = EmceeEngine(agreement_problem(), walkers=8).run(10, initial=conjugate_optimum)
        assert run.attrs["ampere_start_route"] == "scipy"
        assert json.loads(run.attrs["ampere_start"]) == conjugate_optimum.start_record()
        assert run.attrs["ampere_schema_version"] == 9

    def test_a_prior_started_run_says_prior(self) -> None:
        run = EmceeEngine(agreement_problem(), walkers=8).run(10)
        assert run.attrs["ampere_start_route"] == "prior"
        assert "ampere_start" not in run.attrs

    def test_a_callers_array_says_user(self) -> None:
        engine = EmceeEngine(agreement_problem(), walkers=8)
        run = engine.run(10, initial=np.full((8, 1), 2.0) + 0.01 * np.arange(8)[:, None])
        assert run.attrs["ampere_start_route"] == "user"

    def test_the_optimum_carries_the_schema_too(self, conjugate_optimum: Optimum) -> None:
        from ampere.results import PROVENANCE_SCHEMA_VERSION

        assert PROVENANCE_SCHEMA_VERSION == 9
        assert conjugate_optimum.provenance["ampere_schema_version"] == 9
        tree = conjugate_optimum.to_datatree()
        assert tree.attrs["ampere_schema_version"] == 9
        assert "posterior" not in tree.children
