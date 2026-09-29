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
from ampere.inference import EmceeEngine, EngineError, optimise, warm_start_gp
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


# ---------------------------------------------------------------------------
# (a) The native routes, on sed_composition's torch and jax twins
# ---------------------------------------------------------------------------

NATIVE = [
    pytest.param(
        "torch", marks=pytest.mark.skipif(not HAS_TORCH, reason="torch is not installed here")
    ),
    pytest.param("jax", marks=pytest.mark.skipif(not HAS_JAX, reason="jax is not installed here")),
]


@pytest.fixture(scope="module")
def sed_scipy() -> Optimum:
    """The scipy route's mode on the reference twin, shared by the native rows."""
    from examples.sed_composition.sed_composition import build_problem

    return optimise(build_problem("reference"), method="scipy", starts=2)


@pytest.mark.parametrize("backend", NATIVE)
class TestTheMapRoute:
    def test_inside_the_sampled_central_50_percent(
        self, backend: str, sed_posterior: dict[str, np.ndarray], sed_scipy: Optimum
    ) -> None:
        """The native twin composes the same declaration on the same data (same
        seed), so the reference emcee posterior is its posterior too."""
        from examples.sed_composition.sed_composition import build_problem

        problem = build_problem(backend)
        optimum = optimise(problem, method="map", starts=2)
        assert optimum.route == "map" and optimum.backend == backend and optimum.converged
        assert central_50(optimum, sed_posterior) == dict.fromkeys(optimum.free_names, True)
        # the objective convention is one convention: same mode as scipy's
        np.testing.assert_allclose(optimum.unconstrained, sed_scipy.unconstrained, atol=1e-3)
        # and the autodiff curvature is the finite-difference one
        assert optimum.covariance is not None and sed_scipy.covariance is not None
        np.testing.assert_allclose(
            np.sqrt(np.diag(optimum.covariance)), np.sqrt(np.diag(sed_scipy.covariance)), rtol=0.02
        )

    def test_auto_takes_the_map_route(self, backend: str) -> None:
        from examples.sed_composition.sed_composition import build_problem

        assert optimise(build_problem(backend), starts=1).route == "map"


@pytest.mark.parametrize("backend", NATIVE)
class TestTheVIRoute:
    def test_inside_the_sampled_central_50_percent(
        self, backend: str, sed_posterior: dict[str, np.ndarray]
    ) -> None:
        from examples.sed_composition.sed_composition import build_problem

        optimum = optimise(build_problem(backend), method="vi", starts=8)
        assert optimum.route == "vi" and optimum.backend == backend
        assert optimum.covariance is not None
        assert central_50(optimum, sed_posterior) == dict.fromkeys(optimum.free_names, True)

    def test_an_unknown_option_is_refused(self, backend: str) -> None:
        from examples.sed_composition.sed_composition import build_problem

        with pytest.raises(EngineError, match=r"does not take the option\(s\) \['minimiser'\]"):
            optimise(build_problem(backend), method="vi", minimiser="Powell")


# ---------------------------------------------------------------------------
# (b) The warm start, on test_hsgp.py's solver fixture
# ---------------------------------------------------------------------------

#: ``tests/core/test_hsgp.py``'s ``TestTheSolver`` fixture — a Matérn-3/2 of
#: amplitude 0.4 and length scale 2.0 over [1, 11.5], white noise of 0.1, a
#: 24-member basis at boundary factor 2 — made into a fitting problem: the same
#: kernel and noise drawn onto 80 points (the fixture's 12 are too few to pin
#: three hyperparameters to a factor of two) on top of a fixed power law.
HSGP_AMPLITUDE, HSGP_LENGTH_SCALE, HSGP_SIGMA = 0.4, 2.0, 0.1
HSGP_GRID = np.linspace(1.0, 11.5, 80)
HSGP_BASIS = 24


def _gp_data() -> np.ndarray:
    from ampere.core import Matern32

    kernel = Matern32(HSGP_AMPLITUDE, HSGP_LENGTH_SCALE)
    values = kernel.resolve({})
    separation = np.abs(HSGP_GRID[:, None] - HSGP_GRID[None, :])
    covariance = np.asarray(kernel.value(separation, values), dtype=float)
    rng = np.random.default_rng(11)
    factor = np.linalg.cholesky(covariance + 1e-10 * np.eye(HSGP_GRID.size))
    truth = 2.0 * HSGP_GRID**-1.0
    return (
        truth
        + factor @ rng.standard_normal(HSGP_GRID.size)
        + rng.normal(0, HSGP_SIGMA, HSGP_GRID.size)
    )


GP_DATA = _gp_data()


def gp_problem(solver: Any = None) -> FittingProblem:
    from ampere.core import (
        GaussianFamily,
        GaussianProcessNoise,
        HilbertSpaceGP,
        Likelihood,
        Matern32,
    )

    chosen = (
        HilbertSpaceGP(basis_size=HSGP_BASIS, boundary_factor=2.0) if solver is None else solver
    )
    likelihood = Likelihood(
        GaussianFamily(),
        GaussianProcessNoise(
            Matern32(st.loguniform(0.05, 5.0), st.loguniform(0.3, 20.0)),
            chosen,
            scale=st.loguniform(0.3, 3.0),
        ),
    )
    observed = Spectrum(
        HSGP_GRID * u.micron,
        GP_DATA * u.Jy,
        uncertainty=np.full(HSGP_GRID.size, HSGP_SIGMA) * u.Jy,
    )
    return FittingProblem(
        PowerLaw(HSGP_GRID, norm=2.0, index=-1.0, reference_wavelength=1.0),
        [Dataset(observed, likelihood=likelihood)],
        seed=SEED,
    )


HYPERPARAMETERS = (
    "default.likelihood.amplitude",
    "default.likelihood.length_scale",
    "default.likelihood.scale",
)

#: The emcee run the warm start is held to: 12 walkers, 400 steps, 150
#: burned — about 6 s on dev over three hyperparameters.
GP_WALKERS, GP_STEPS, GP_BURN_IN = 12, 400, 150


@pytest.fixture(scope="module")
def gp_posterior_medians() -> dict[str, float]:
    run = EmceeEngine(gp_problem(), walkers=GP_WALKERS).run(GP_STEPS, burn_in=GP_BURN_IN)
    posterior = run["posterior"].dataset
    return {name: float(np.median(np.asarray(posterior[name]))) for name in HYPERPARAMETERS}


class TestTheWarmStart:
    def test_within_a_factor_of_two_of_the_posterior_median(
        self, gp_posterior_medians: dict[str, float]
    ) -> None:
        optimum = warm_start_gp(gp_problem())["default"]
        assert optimum.route == "empirical_bayes"
        assert optimum.free_names == HYPERPARAMETERS
        ratios = {n: optimum.constrained[n] / gp_posterior_medians[n] for n in HYPERPARAMETERS}
        print("warm-start / posterior-median ratios:", ratios)
        assert all(0.5 <= ratio <= 2.0 for ratio in ratios.values()), ratios

    def test_the_dense_marginal_likelihood_agrees_at_the_answer(self) -> None:
        """At the returned hyperparameters, ``DenseGP`` and the reduced-rank
        solver score the same residual alike, to the approximation's tolerance."""
        from ampere.core import DenseGP, HilbertSpaceGP, Matern32

        optimum = warm_start_gp(gp_problem())["default"]
        a, ell, s = (optimum.constrained[n] for n in HYPERPARAMETERS)
        kernel = Matern32(a, ell)
        values = kernel.resolve({})
        coordinates = HSGP_GRID.reshape(-1, 1)
        residual = GP_DATA - 2.0 * HSGP_GRID**-1.0
        variance = np.full(HSGP_GRID.size, (s * HSGP_SIGMA) ** 2)
        dense = DenseGP().log_marginal_likelihood(kernel, coordinates, residual, variance, values)

        def reduced(m: int) -> float:
            return HilbertSpaceGP(basis_size=m, boundary_factor=2.0).log_marginal_likelihood(
                kernel, coordinates, residual, variance, values
            )

        coarse, fine = reduced(HSGP_BASIS), reduced(4 * HSGP_BASIS)
        print(f"dense {dense:.6f}, reduced-rank m=24 {coarse:.6f}, m=96 {fine:.6f}")
        # The approximation's tolerance at m=24: a Matérn-3/2 spectrum decays
        # only as omega**-4, so the truncated basis is a few per cent short of
        # the dense marginal likelihood; the claim is agreement to that
        # tolerance and convergence towards it as the basis grows.
        assert coarse == pytest.approx(dense, rel=0.05)
        assert abs(fine - dense) < abs(coarse - dense)

    def test_a_brute_force_grid_does_not_beat_the_root_find(self) -> None:
        from ampere.core import HilbertSpaceGP, Matern32

        optimum = warm_start_gp(gp_problem())["default"]
        a, ell, s = (optimum.constrained[n] for n in HYPERPARAMETERS)
        coordinates = HSGP_GRID.reshape(-1, 1)
        residual = GP_DATA - 2.0 * HSGP_GRID**-1.0
        solver = HilbertSpaceGP(basis_size=HSGP_BASIS, boundary_factor=2.0)

        def score(amplitude: float, scale: float) -> float:
            kernel = Matern32(amplitude, ell)
            variance = np.full(HSGP_GRID.size, (scale * HSGP_SIGMA) ** 2)
            return solver.log_marginal_likelihood(
                kernel, coordinates, residual, variance, kernel.resolve({})
            )

        found = score(a, s)
        best = max(
            score(x, y)
            for x in np.geomspace(a / 3.0, a * 3.0, 61)
            for y in np.geomspace(s / 3.0, s * 3.0, 61)
        )
        print(f"root find {found:.9f}, best of the 61x61 grid {best:.9f}")
        assert best <= found + 1e-6

    def test_the_dense_solver_takes_the_same_start(self) -> None:
        from ampere.core import DenseGP

        dense = warm_start_gp(gp_problem(DenseGP()))["default"]
        reduced = warm_start_gp(gp_problem())["default"]
        assert "temporary HilbertSpaceGP(basis_size=32)" in dense.message
        for name in HYPERPARAMETERS:
            assert dense.constrained[name] == pytest.approx(reduced.constrained[name], rel=0.1)

    def test_it_combines_with_a_model_optimum_and_seeds_a_run(self) -> None:
        import json

        problem = gp_problem()
        warm = warm_start_gp(problem)["default"]
        run = EmceeEngine(problem, walkers=8).run(10, initial=warm)
        assert run.attrs["ampere_start_route"] == "empirical_bayes"
        assert json.loads(run.attrs["ampere_start"])["identity"] == warm.identity

    def test_no_gp_dataset_is_refused_by_name(self) -> None:
        with pytest.raises(
            EngineError, match="no dataset in this problem has a GaussianProcessNoise"
        ):
            warm_start_gp(agreement_problem())


# ---------------------------------------------------------------------------
# (c) NUTS from the MAP against NUTS from the prior; (d) emcee's burn-in
# ---------------------------------------------------------------------------

#: The pinned NUTS budgets, per backend: short enough that the prior-started
#: run has not converged, the same for both runs, at the same seed, as
#: ``(warmup, draws, chains, max_tree_depth)``. The tree depth is capped for
#: both runs, so a prior-started chain far out in a tail costs at most
#: ``2**depth`` leapfrog steps an iteration rather than 1024 — the comparison
#: stays fair and the row stays affordable.
#:
#: jax keeps 4 chains of 100 draws at depth 6 (the ruling's 200 draws halved,
#: orchestrator's note): about three minutes, and a wide margin — the
#: MAP-started R-hat near 1.1 against the prior-started run's 2.3 to 3.2.
#: torch is cut to 2 chains of 50 draws at depth 4 (W6.7's successor): at the
#: jax budget the torch row ran past fifty minutes on a shared machine, since
#: pyro costs about 0.75 s an iteration from the prior and 1.7 s from the MAP
#: at depth 4, 5.8 s at depth 6 (the MAP-started chains build full trees at
#: their small adapted step size). Two chains are the fewest R-hat can
#: compare; the torch row takes about nine minutes and its margin is still
#: clear. The same cut on jax left too thin a margin (MAP-started R-hat up to
#: 2.09 against the prior's 2.09 to 3.09), hence the split. At warm-up 50 the
#: prior-started run is far from converged on both backends, so the warm-up
#: stays at 50. torch is held to four threads for the row.
NUTS_BUDGET: dict[str, tuple[int, int, int, int]] = {
    "jax": (50, 100, 4, 6),
    "torch": (50, 50, 2, 4),
}


def _nuts_summary(run: Any) -> dict[str, Any]:
    import arviz

    posterior = run["posterior"]
    rhat = arviz.rhat(posterior)
    ess = arviz.ess(posterior)
    names = [str(n) for n in posterior.dataset.data_vars]
    return {
        "rhat": {n: float(np.asarray(rhat[n])) for n in names},
        "ess": {n: float(np.asarray(ess[n])) for n in names},
        "divergences": int(run.attrs["ampere_nuts_divergences"]),
    }


@pytest.mark.parametrize("backend", NATIVE)
def test_nuts_from_the_map_adapts_faster_than_from_the_prior(backend: str) -> None:

    if backend == "torch":
        import torch  # pyrefly: ignore[missing-import]

        threads = torch.get_num_threads()
        torch.set_num_threads(4)
    try:
        runs = _nuts_pair(backend)
    finally:
        if backend == "torch":
            torch.set_num_threads(threads)
    print(f"NUTS on {backend}, (warmup, draws, chains, depth)={NUTS_BUDGET[backend]}:", runs)
    prior, started = runs["prior"], runs["map"]
    for name in prior["rhat"]:
        assert started["rhat"][name] < prior["rhat"][name], (name, runs)
        assert started["ess"][name] >= prior["ess"][name], (name, runs)
    assert started["divergences"] <= prior["divergences"], runs


def _nuts_pair(backend: str) -> dict[str, Any]:
    from examples.sed_composition.sed_composition import build_problem

    from ampere.inference import NUTSEngine

    warmup, draws, chains, depth = NUTS_BUDGET[backend]
    optimum = optimise(build_problem(backend), method="map", starts=2)
    runs = {}
    for label, initial in (("prior", None), ("map", optimum)):
        runs[label] = _nuts_summary(
            NUTSEngine(build_problem(backend)).run(
                draws,
                warmup=warmup,
                chains=chains,
                max_tree_depth=depth,
                initial=initial,
            )
        )
    return runs


#: emcee's budget for the burn-in comparison: the example's 16 walkers, 600
#: steps, the plateau the mean over the last quarter.
BURN_WALKERS, BURN_STEPS = 16, 600


def burn_in_step(log_prob: np.ndarray) -> int:
    """The first step at which the ensemble's mean log-probability is within
    one unit of its plateau (the mean over the last quarter of the run)."""
    mean = np.mean(log_prob, axis=1)
    plateau = float(np.mean(mean[-(mean.size // 4) :]))
    return int(np.argmax(mean >= plateau - 1.0))


def test_emcee_from_the_optimum_burns_in_faster(sed_scipy: Optimum) -> None:
    from examples.sed_composition.sed_composition import build_problem

    steps = {}
    for label, initial in (("prior", None), ("optimum", sed_scipy)):
        engine = EmceeEngine(build_problem("reference"), walkers=BURN_WALKERS)
        engine.run(BURN_STEPS, initial=initial)
        steps[label] = burn_in_step(np.asarray(engine.sampler.get_log_prob()))
    print("emcee burn-in equivalent (steps):", steps)
    assert steps["optimum"] < steps["prior"], steps
