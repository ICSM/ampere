"""The variational driver, end to end, on every backend it can fit.

``test_nuts.py`` holds the gradient-based *sampler*; this file holds the
gradient-based *optimiser*, and the two are checked against the same oracle for
the same reason: comparing an approximation only against another approximation
would pass two methods that are wrong in the same way.

What this file claims, and how each claim is checked:

1. **The fit recovers a posterior written down in closed form.** The agreement
   problem is deliberately conjugate — a power law with its index held fixed is
   *linear* in ``norm``, so a Gaussian prior and Gaussian noise give a Gaussian
   posterior whose mean and standard deviation are arithmetic. A Gaussian
   posterior is also the one case where a Gaussian guide is not an
   approximation at all, which is what makes a **tight** tolerance legitimate
   here: mean-field VI on a one-dimensional Gaussian target is exact in the
   limit, so a wide tolerance would hide a real error rather than absorb an
   honest one.
2. **The approximation shows up where it should.** On a *correlated* posterior
   the mean-field guide underestimates the marginal variances and the
   full-covariance guide does not — the textbook failure, asserted rather than
   described, so that ``vi_guide`` in a run's attrs means something a reader
   can act on.
3. **A run is a run.** The same ``DataTree`` every other driver emits, with the
   ELBO trace and the guide family in the attrs, ``ampere_realised = 1``, and
   ``engine_draws_recomputed = 0`` — the last because this driver consumes
   ``inference.md`` §10a's ``log_likelihood_terms`` rather than recomputing the
   per-dataset split on the numpy path.
4. **It refuses, by name, what it cannot fit**: a backend with no variational
   library, a non-differentiable problem, an unknown guide family, and (W5.14)
   a flow over a single parameter.
5. **Each guide family is what it says it is** (W5.14). ``laplace`` is a
   curvature approximation, so it is checked where a curvature approximation
   is *exact* — a Gaussian posterior — against the same closed form claim 1
   uses; ``flow`` keeps no fitted parameters at all, so it is checked against
   the one thing every guide has independently of its shape, the library's own
   ELBO.
6. **Every guide is SBC-ranked** through ``ampere.results.calibration.sbc``,
   which is W5.14's engine battery applied to a driver whose "engine" is the
   guide family. At the per-PR budget that is a check of the machinery (ranks
   produced, nothing failed); the budget at which Talts et al.'s uniformity
   test has power is behind ``-m engines_full``, as it is for the nested
   samplers in ``test_nested.py``.

Parametrised over the backends installed here that this driver supports, for
the reason ``test_nuts.py`` gives: writing the claim once per backend by hand
is how two backends drift apart. Since W2.5 slice 2 that is both of them —
pyro on torch, numpyro on jax — so every claim above is one claim per route,
including the guide-family one, which is the row most likely to depend on a
library's autoguide implementation rather than on the mathematics. Only one
row is route-specific, and it says why.

Budgets are small and seeds fixed, so this belongs in the per-PR gate. VI is
cheap — a fit is a few thousand cheap gradient steps, not a chain — which is
part of why it exists.
"""

from __future__ import annotations

import dataclasses
import importlib
import json
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
    FittingProblem,
    GaussianFamily,
    Likelihood,
    Spectrum,
)
from ampere.inference import EmceeEngine, VIEngine
from ampere.results.calibration import sbc
from ampere.inference.engine import unconstrained_jacobian_correction
from ampere.inference.exceptions import EngineError

SEED = 20260907
REFERENCE_WAVELENGTH = 1.0


def power_law(grid: np.ndarray, norm: float, index: float) -> np.ndarray:
    return norm * (grid / REFERENCE_WAVELENGTH) ** index


def noisy(grid: np.ndarray, truth: np.ndarray, sigma: float, seed: int) -> Spectrum:
    rng = np.random.default_rng(seed)
    return Spectrum(
        grid * u.um,
        (truth + rng.normal(0.0, sigma, grid.size)) * u.Jy,
        uncertainty=np.full(grid.size, sigma) * u.Jy,
    )


@dataclasses.dataclass(frozen=True)
class Kit:
    """One backend's pieces, by name, so a test body never names a library."""

    name: str
    module: Any

    def likelihood(self) -> Likelihood:
        return Likelihood(GaussianFamily(), self.module.IndependentNoise())


def _installed_kits() -> list[Kit]:
    from ampere.inference._vi import VARIATIONAL_LIBRARIES

    found: list[Kit] = []
    for name in sorted(VARIATIONAL_LIBRARIES):
        try:
            module = importlib.import_module(f"ampere.backends.{name}")
        except ImportError:  # the extra is not installed in this environment
            continue
        if name == "jax":
            module.configure_x64()
        found.append(Kit(name, module))
    return found


KITS = _installed_kits()

pytestmark = pytest.mark.skipif(
    not KITS,
    reason="no backend with a variational library installed",
)


@pytest.fixture(scope="module", params=[kit.name for kit in KITS])
def kit(request: Any) -> Kit:
    return next(found for found in KITS if found.name == request.param)


# ---------------------------------------------------------------------------
# The conjugate problem, and its exact posterior
# ---------------------------------------------------------------------------

GRID = np.geomspace(1.0, 10.0, 20)
SIGMA = 0.2
INDEX = -1.0
PRIOR = (2.0, 0.5)  # (mean, sd) of the Gaussian prior on `norm`
DATA = noisy(GRID, power_law(GRID, 2.0, INDEX), SIGMA, seed=7)


def analytic() -> tuple[float, float]:
    """The posterior on ``norm``, in closed form: ``(mean, sd)``.

    The same conjugate construction ``test_nuts.py`` uses, and deliberately the
    same numbers: the two drivers are then checked against one oracle, so a
    disagreement between them is a disagreement with arithmetic rather than
    with each other.
    """
    x = (GRID / REFERENCE_WAVELENGTH) ** INDEX
    y = np.asarray(DATA.values)
    mu, tau = PRIOR
    precision = 1.0 / tau**2 + float(np.sum(x**2)) / SIGMA**2
    mean = (mu / tau**2 + float(np.sum(x * y)) / SIGMA**2) / precision
    return mean, 1.0 / math.sqrt(precision)


def conjugate_problem(kit: Kit, seed: int | None = SEED) -> FittingProblem:
    return FittingProblem(
        kit.module.PowerLaw(
            GRID,
            norm=st.norm(*PRIOR),
            index=INDEX,
            reference_wavelength=REFERENCE_WAVELENGTH,
        ),
        [Dataset(DATA, likelihood=kit.likelihood())],
        seed=seed,
    )


def correlated_problem(kit: Kit, seed: int | None = SEED) -> FittingProblem:
    """Both power-law parameters free, so the posterior is strongly correlated.

    ``norm`` and ``index`` trade off against each other on a short lever arm:
    raising the index and lowering the normalisation describes almost the same
    spectrum. That is the geometry a mean-field guide cannot represent, and it
    is not contrived — it is the geometry of every SED fit ampere runs.
    """
    return FittingProblem(
        kit.module.PowerLaw(
            GRID,
            norm=st.norm(2.0, 1.0),
            index=st.norm(-1.0, 1.0),
            reference_wavelength=REFERENCE_WAVELENGTH,
        ),
        [Dataset(DATA, likelihood=kit.likelihood())],
        seed=seed,
    )


def fit(problem: FittingProblem, **settings: Any) -> Any:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return VIEngine(problem).run(**settings)


@pytest.fixture(scope="module")
def conjugate_run(kit: Kit) -> Any:
    return fit(conjugate_problem(kit), draws=2000, steps=3000)


# ---------------------------------------------------------------------------
# 1. Agreement with the closed form
# ---------------------------------------------------------------------------


class TestAgreementWithTheClosedForm:
    def test_the_posterior_mean_matches_the_conjugate_answer(self, conjugate_run: Any) -> None:
        mean, sd = analytic()
        drawn = np.asarray(conjugate_run["posterior"]["model.norm"]).ravel()
        assert float(drawn.mean()) == pytest.approx(mean, abs=0.2 * sd)

    def test_the_posterior_width_matches_the_conjugate_answer(self, conjugate_run: Any) -> None:
        """The half a Gaussian guide *can* get right, so it must.

        A guide family that had collapsed — the classic VI failure, where the
        ELBO is optimised into a spike — would pass the mean row and fail this
        one, which is why both are here.
        """
        _, sd = analytic()
        drawn = np.asarray(conjugate_run["posterior"]["model.norm"]).ravel()
        assert float(drawn.std(ddof=1)) == pytest.approx(sd, rel=0.15)

    def test_the_elbo_climbed(self, conjugate_run: Any) -> None:
        import json

        trace = json.loads(conjugate_run.attrs["ampere_vi_elbo_trace"])
        assert len(trace) > 1
        # Not monotone -- it is a stochastic estimate -- but the end must be
        # well above the start, or the optimiser did nothing.
        assert trace[-1] > trace[0]
        assert conjugate_run.attrs["ampere_vi_final_elbo"] == pytest.approx(trace[-1])


# ---------------------------------------------------------------------------
# 2. The approximation is where the guide family says it is
# ---------------------------------------------------------------------------


class TestWhatTheGuideFamilyAssumes:
    def test_mean_field_underestimates_a_correlated_posteriors_width(self, kit: Kit) -> None:
        """The textbook failure, asserted rather than described.

        A mean-field guide fits the *conditional* widths rather than the
        marginal ones, so on a correlated posterior it is too narrow — and a
        user who reads ``vi_guide == "normal"`` in an archived run needs that
        to be a fact about the code rather than a caution in a docstring.
        """
        mean_field = fit(correlated_problem(kit), draws=3000, steps=4000, guide="normal")
        full = fit(correlated_problem(kit), draws=3000, steps=4000, guide="multivariate")
        for name in ("model.norm", "model.index"):
            narrow = float(np.asarray(mean_field["posterior"][name]).std(ddof=1))
            wide = float(np.asarray(full["posterior"][name]).std(ddof=1))
            assert narrow < wide

    def test_the_full_covariance_guide_recovers_the_correlation(self, kit: Kit) -> None:
        """And the mean-field one reports none, because it cannot."""
        mean_field = fit(correlated_problem(kit), draws=3000, steps=4000, guide="normal")
        full = fit(correlated_problem(kit), draws=3000, steps=4000, guide="multivariate")

        def correlation(run: Any) -> float:
            a = np.asarray(run["posterior"]["model.norm"]).ravel()
            b = np.asarray(run["posterior"]["model.index"]).ravel()
            return float(np.corrcoef(a, b)[0, 1])

        assert abs(correlation(mean_field)) < 0.15
        assert abs(correlation(full)) > 0.5

    def test_the_full_covariance_guide_names_its_own_approximation_family(self, kit: Kit) -> None:
        """W5.0: ``"multivariate"``, not the mean-field default — cheap budget, attrs only."""
        full = fit(correlated_problem(kit), draws=10, steps=10, guide="multivariate")
        assert full.attrs["ampere_approximation"] == "multivariate"


# ---------------------------------------------------------------------------
# 3. Every run emits the run
# ---------------------------------------------------------------------------


class TestTheRunItEmits:
    def test_the_shape_is_one_chain_of_independent_draws(self, conjugate_run: Any) -> None:
        assert np.asarray(conjugate_run["posterior"]["model.norm"]).shape == (1, 2000)

    def test_the_attrs_name_the_engine_the_backend_and_the_guide(
        self, conjugate_run: Any, kit: Kit
    ) -> None:
        assert conjugate_run.attrs["ampere_engine"] == "vi"
        assert conjugate_run.attrs["ampere_backend"] == kit.name
        assert conjugate_run.attrs["ampere_vi_guide"] == "normal"
        assert conjugate_run.attrs["ampere_vi_guide_class"] == "AutoNormal"
        assert conjugate_run.attrs["ampere_vi_steps"] == 3000
        assert conjugate_run.attrs["ampere_vi_library"] in {"pyro", "numpyro"}

    def test_the_attrs_name_the_approximation_family(self, conjugate_run: Any) -> None:
        """W5.0: ``ampere_approximation`` is the guide family, not ``"none"``."""
        assert conjugate_run.attrs["ampere_approximation"] == "mean_field"

    def test_the_proposal_log_density_is_beside_the_true_split(self, conjugate_run: Any) -> None:
        """W5.0: one per draw, finite, in ``sample_stats`` beside ``lp``."""
        stats = conjugate_run["sample_stats"].dataset
        proposal = np.asarray(stats["proposal_log_density"])
        assert proposal.shape == np.asarray(stats["lp"]).shape
        assert np.all(np.isfinite(proposal))

    def test_the_draws_came_through_the_realisation(self, conjugate_run: Any) -> None:
        assert conjugate_run.attrs["ampere_realised"] == 1

    def test_nothing_was_recomputed_on_the_numpy_path(self, conjugate_run: Any) -> None:
        """``inference.md`` §10a's optional member, consumed (W2.4 slice 2).

        The decomposition comes from the realisation, so no stored draw is
        re-scored through ``FittingProblem.evaluate``. A non-zero count here
        would mean the driver had silently fallen back — which is allowed by
        the contract but is not what this driver does, and the difference is a
        full model evaluation per draw.
        """
        assert conjugate_run.attrs["ampere_engine_draws_recomputed"] == 0

    def test_the_per_dataset_decomposition_is_present_and_sums_correctly(
        self, conjugate_run: Any
    ) -> None:
        groups = list(conjugate_run["log_likelihood"].data_vars)
        assert groups == ["default"]
        total = np.asarray(conjugate_run["log_likelihood"]["default"])
        joint = np.asarray(conjugate_run["sample_stats"]["lp"]) - np.asarray(
            conjugate_run["sample_stats"]["log_prior"]
        )
        assert total == pytest.approx(joint, abs=1e-9)

    def test_it_repeats_exactly_from_the_problems_seed(self, kit: Kit) -> None:
        first = fit(conjugate_problem(kit, SEED), draws=50, steps=100)
        second = fit(conjugate_problem(kit, SEED), draws=50, steps=100)
        assert np.asarray(first["posterior"]["model.norm"]) == pytest.approx(
            np.asarray(second["posterior"]["model.norm"]), abs=0.0
        )

    def test_it_leaves_no_global_state_behind(self, kit: Kit) -> None:
        """pyro's parameter store is process-global; a fit must not write into it.

        A guide registers its parameters under site names, so a second fit in
        the same interpreter would find the first one's and start from them —
        a silently different answer rather than a crash. The fit runs inside
        ``pyro.get_param_store().scope()`` for exactly this, and the row proves
        the scope closes.

        **A claim about pyro, so it is asked of the pyro route only** (W2.5
        slice 2). numpyro has no process-global parameter store at all: the
        fitted parameters come back in the ``SVIRunResult``, so there is no
        state for a jax fit to leave behind and nothing here to assert. The
        row is skipped rather than deleted or generalised, because the hazard
        it guards is real on one route and structurally absent on the other,
        and saying which is more useful than a row that passes vacuously.
        """
        if kit.name != "torch":
            pytest.skip(
                f"the {kit.name!r} route uses no process-global parameter store; numpyro returns "
                f"its fitted parameters in the SVIRunResult, so there is nothing to leak"
            )
        pyro = importlib.import_module("pyro")
        before = set(pyro.get_param_store().keys())
        fit(conjugate_problem(kit), draws=10, steps=20)
        assert set(pyro.get_param_store().keys()) == before


# ---------------------------------------------------------------------------
# 4. Refusals
# ---------------------------------------------------------------------------


class TestRefusals:
    def test_a_problem_on_another_backend_is_refused_by_name(self, kit: Kit) -> None:
        problem = FittingProblem(
            ReferencePowerLaw(
                GRID, norm=st.norm(*PRIOR), index=INDEX, reference_wavelength=REFERENCE_WAVELENGTH
            ),
            [Dataset(DATA, likelihood=Likelihood(GaussianFamily()))],
            seed=SEED,
        )
        with pytest.raises(EngineError, match="cannot fit a problem"):
            VIEngine(problem)

    def test_an_unknown_guide_family_is_refused_by_name(self, kit: Kit) -> None:
        """``"laplace"`` was this row's unknown family until W5.14 added it.

        Which is the point of keeping the row parametrised on a name no
        library here implements: the refusal must name what *is* available,
        and the set of available families is now something that grows.
        """
        engine = VIEngine(conjugate_problem(kit))
        with pytest.raises(EngineError, match="guide family"):
            engine.run(draws=10, guide="student_t")

    def test_a_non_positive_step_count_is_refused(self, kit: Kit) -> None:
        engine = VIEngine(conjugate_problem(kit))
        with pytest.raises(EngineError, match="optimiser step"):
            engine.run(draws=10, steps=0)

    def test_a_non_positive_learning_rate_is_refused(self, kit: Kit) -> None:
        engine = VIEngine(conjugate_problem(kit))
        with pytest.raises(EngineError, match="learning rate"):
            engine.run(draws=10, steps=10, learning_rate=0.0)

    def test_a_wrongly_shaped_start_point_is_refused(self, kit: Kit) -> None:
        engine = VIEngine(conjugate_problem(kit))
        with pytest.raises(EngineError, match="element"):
            engine.run(draws=10, steps=10, initial=np.zeros(7))

    def test_supported_backends_is_answered_from_the_registry(self, kit: Kit) -> None:
        from ampere.inference._vi import supported_backends

        assert kit.name in supported_backends()
        assert "reference" not in supported_backends()


# ---------------------------------------------------------------------------
# 5. The proposal density is really the guide's (review-caught regression)
# ---------------------------------------------------------------------------


def _unconstrained_draws(problem: FittingProblem, run: Any) -> np.ndarray:
    """A run's stored (constrained) draws, mapped back to unconstrained space.

    Nothing about a stored run keeps the unconstrained vector VI actually
    fitted in — only ``unconstrain``'s own inverse of ``constrain`` recovers
    it, exact up to floating point since the two are one bijection's forward
    and inverse maps.
    """
    posterior = run["posterior"].dataset
    names = list(problem.parameters.free_names)
    constrained = np.stack([np.asarray(posterior[name]).ravel() for name in names], axis=1)
    return np.stack(
        [
            problem.unconstrain(problem.parameters.pack(dict(zip(names, row, strict=True))))
            for row in constrained
        ]
    )


class TestTheProposalDensityIsReallyTheGuides:
    """A review of this item caught a real bug here, twice, before this landed.

    Both autoguide classes route a fitted draw through an **auxiliary**
    sample site under the real ``Normal``/``MultivariateNormal`` the guide
    optimised, and only then report ``_SITE`` itself — as a ``Delta`` at the
    identity-transformed value, whose ``log_prob`` is the change-of-variables
    term between the two sites (zero here, since the transform is the
    identity), **not** the guide's density. Reading ``trace.nodes[_SITE]
    ["fn"].log_prob(...)`` alone therefore stored a constant zero for every
    draw on both the pyro and the numpyro route — the importance-reweighting
    test below happened to still pass, for the wrong reason, because
    ``proposal_log_density`` differed from the truth by exactly the same
    (then-missing) constant on every draw and a self-normalised weight is
    invariant to an additive constant in the log. ``np.ptp(proposal) > 0`` is
    the one-line guard that would have caught it outright — a genuine
    per-draw density is never constant across 500 independent draws from a
    continuous guide — and the direct comparison against an independently
    built ``scipy`` density is the check that the *value*, not merely its
    variation, is right.
    """

    DRAWS = 500
    STEPS = 200

    def _check(self, kit: Kit, guide: str) -> None:
        problem = correlated_problem(kit)
        engine = VIEngine(problem)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            run = engine.run(draws=self.DRAWS, steps=self.STEPS, guide=guide)
        stats = run["sample_stats"].dataset
        proposal = np.asarray(stats["proposal_log_density"]).ravel()

        # The cheap, always-on guard: a genuine per-draw density is not a
        # constant, which is exactly what the Delta bug produced.
        assert np.ptp(proposal) > 0.0

        unconstrained = _unconstrained_draws(problem, run)
        jacobian = unconstrained_jacobian_correction(problem, unconstrained)
        if guide in {"multivariate", "laplace"}:
            covariance = engine.guide_scale_tril @ engine.guide_scale_tril.T
            independent = st.multivariate_normal(engine.guide_loc, covariance).logpdf(unconstrained)
        else:
            independent = (
                st.norm(engine.guide_loc, engine.guide_scale).logpdf(unconstrained).sum(axis=1)
            )
        assert (proposal + jacobian) == pytest.approx(independent, abs=1e-6)

    def test_the_mean_field_guides_stored_density_agrees_with_scipy(self, kit: Kit) -> None:
        self._check(kit, "normal")

    def test_the_full_covariance_guides_stored_density_agrees_with_scipy(self, kit: Kit) -> None:
        self._check(kit, "multivariate")

    def test_the_laplace_guides_stored_density_agrees_with_scipy(self, kit: Kit) -> None:
        """W5.14: the same check for the family whose Gaussian is a *curvature*.

        Worth making separately rather than trusting the row above, because
        the two families arrive at their ``scale_tril`` by different roads and
        only one of them is a fitted parameter: pyro hands back a whole
        ``AutoMultivariateNormal`` from ``laplace_approximation()`` (so its
        covariance is the row-scaled correlation split, as for the fitted
        guide), while numpyro has no such object and the factor must be read
        off ``get_transform``'s affine. A driver that read either wrongly
        would store a density for a Gaussian it did not draw from.
        """
        self._check(kit, "laplace")


# ---------------------------------------------------------------------------
# 6. W5.0's accept criterion: reweighting from stored groups alone
# ---------------------------------------------------------------------------


class TestImportanceCorrectedPosteriorAgreesWithEmcee:
    """The contract's own check: not merely that the column exists, but that

    ``exp(log_prior + log_likelihood - proposal_log_density)``, built from
    nothing but ``run["sample_stats"]`` and ``run["posterior"]``, is a valid
    importance weight that a genuinely biased VI fit can be corrected with.

    The guide is deliberately **under-trained** (a handful of SVI steps) on
    the *conjugate* problem, whose posterior is otherwise the one case a
    mean-field guide fits essentially exactly (``TestAgreementWithTheClosedForm``
    above) — so few steps in, its mean is measurably off, which is exactly
    the room a real importance correction needs to demonstrate it does
    something rather than passing because the raw fit was already right.
    Widths are deliberately **not** checked here: self-normalised importance
    sampling from an under-dispersed proposal is well known to underestimate
    a target's second moment (finite draws in the tails), so a width
    assertion would be testing a sampling-theory limitation, not this
    formula's correctness — the mean is where the correction's arithmetic is
    checkable against an independent reference.
    """

    def test_the_reweighted_mean_matches_an_emcee_reference(self, kit: Kit) -> None:
        vi_run = fit(conjugate_problem(kit), draws=8000, steps=15, guide="normal")
        reference = EmceeEngine(conjugate_problem(kit), walkers=16).run(steps=800, burn_in=200)

        stats = vi_run["sample_stats"].dataset
        log_prior = np.asarray(stats["log_prior"]).ravel()
        log_likelihood = np.asarray(stats["log_likelihood"]).ravel()
        proposal_log_density = np.asarray(stats["proposal_log_density"]).ravel()
        log_weight = log_prior + log_likelihood - proposal_log_density
        log_weight -= log_weight.max()
        weight = np.exp(log_weight)
        weight /= weight.sum()
        # A degenerate weighting (all mass on a handful of draws) would pass
        # a mean-agreement check for the wrong reason; guard against it.
        effective_sample_size = 1.0 / np.sum(weight**2)
        assert effective_sample_size > 500

        drawn = np.asarray(vi_run["posterior"]["model.norm"]).ravel()
        raw_mean = float(drawn.mean())
        reweighted_mean = float(np.sum(weight * drawn))
        reference_values = np.asarray(reference["posterior"]["model.norm"]).ravel()
        reference_mean = float(reference_values.mean())
        reference_sd = float(reference_values.std())

        # The under-trained guide's raw mean is measurably off...
        assert abs(raw_mean - reference_mean) > 0.5 * reference_sd
        # ...the correction moves it substantially closer to the reference...
        assert abs(reweighted_mean - reference_mean) < 0.5 * abs(raw_mean - reference_mean)
        # ...and lands within one reference standard deviation of it.
        assert reweighted_mean == pytest.approx(reference_mean, abs=reference_sd)


# ---------------------------------------------------------------------------
# 7. W5.14's two guides: what each one is, and what each one keeps
# ---------------------------------------------------------------------------

#: The flow's budget. It is longer than the Gaussian families' because a flow
#: has more to learn than a location and a scale, not because it is started
#: badly -- ``_vi.py`` moves its base to the same start point the others are
#: initialised at, and the measurement that made that necessary is in the
#: module docstring there.
FLOW_STEPS = 2000
FLOW_LEARNING_RATE = 0.05

#: Draws taken from a fitted flow, and deliberately far fewer than the two
#: thousand the Gaussian families are asked for elsewhere in this file. A
#: guide draw is cheap only when the guide is cheap: both routes take a flow
#: draw by *tracing* the guide, once per draw in Python, and a trace of an
#: inverse-autoregressive flow evaluates three autoregressive networks and
#: their inverses where a trace of ``AutoNormal`` evaluates one Gaussian.
#: Measured on the numpyro route, two thousand flow draws dominated this
#: file's wall time; five hundred is well inside the Monte-Carlo error every
#: claim below is asserted at.
FLOW_DRAWS = 500


@pytest.fixture(scope="module")
def laplace_run(kit: Kit) -> Any:
    """The laplace guide on the **conjugate** problem, where it is exact."""
    return fit(conjugate_problem(kit), draws=2000, steps=1500, guide="laplace")


@pytest.fixture(scope="module")
def flow_run(kit: Kit) -> Any:
    """The flow guide on the correlated problem, at a converged budget."""
    return fit(
        correlated_problem(kit),
        draws=FLOW_DRAWS,
        steps=FLOW_STEPS,
        guide="flow",
        learning_rate=FLOW_LEARNING_RATE,
    )


class TestTheLaplaceGuide:
    """A Laplace approximation is the curvature at the mode, and nothing else.

    Which is why the claim is made where it can be made exactly. On a
    conjugate problem the posterior *is* Gaussian, so the inverse Hessian at
    the MAP point is the posterior covariance to within arithmetic: a tight
    agreement with ``analytic()`` here is therefore a check of the second
    derivative the driver asked the backend for, and not merely of an
    optimiser having got somewhere reasonable. A driver that took the Hessian
    at the wrong point, or that forgot that pyro's ``AutoLaplaceApproximation``
    is a ``Delta`` guide until ``laplace_approximation()`` is called on it,
    fails this row rather than quietly returning the MAP point 2000 times.
    """

    def test_the_posterior_mean_matches_the_conjugate_answer(self, laplace_run: Any) -> None:
        mean, _ = analytic()
        drawn = np.asarray(laplace_run["posterior"]["model.norm"]).ravel()
        assert float(drawn.mean()) == pytest.approx(mean, abs=0.02)

    def test_the_posterior_width_matches_the_conjugate_answer(self, laplace_run: Any) -> None:
        """The row the Hessian has to be right for: a MAP point has no width."""
        _, sd = analytic()
        drawn = np.asarray(laplace_run["posterior"]["model.norm"]).ravel()
        assert float(drawn.std(ddof=1)) == pytest.approx(sd, rel=0.1)

    def test_the_draws_are_not_the_mode_repeated(self, laplace_run: Any) -> None:
        """The failure mode this family invites, asserted against directly."""
        drawn = np.asarray(laplace_run["posterior"]["model.norm"]).ravel()
        assert float(np.ptp(drawn)) > 0.0
        assert len(np.unique(drawn)) > drawn.size // 2

    def test_the_attrs_name_the_guide_and_its_approximation(self, laplace_run: Any) -> None:
        attrs = laplace_run.attrs
        assert attrs["ampere_vi_guide"] == "laplace"
        assert attrs["ampere_vi_guide_class"] == "AutoLaplaceApproximation"
        assert attrs["ampere_approximation"] == "laplace"
        assert attrs["ampere_vi_guide_parameters"] == "loc, scale_tril"

    def test_the_fit_is_kept_as_a_full_covariance(self, kit: Kit) -> None:
        engine = VIEngine(correlated_problem(kit))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            engine.run(draws=200, steps=800, guide="laplace")
        assert engine.guide_loc is not None
        assert engine.guide_scale_tril is not None
        assert engine.guide_scale is None
        assert engine.guide_scale_tril.shape == (2, 2)

    def test_the_curvature_sees_the_correlation(self, kit: Kit) -> None:
        """On the correlated problem the off-diagonal is not an accident."""
        engine = VIEngine(correlated_problem(kit))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            engine.run(draws=200, steps=800, guide="laplace")
        covariance = engine.guide_scale_tril @ engine.guide_scale_tril.T
        correlation = covariance[0, 1] / math.sqrt(covariance[0, 0] * covariance[1, 1])
        assert correlation < -0.4


class TestTheFlowGuide:
    """The only family here that is not a Gaussian, and the only one that

    leaves nothing behind that a ``scipy`` distribution could be rebuilt from.
    That is the whole reason it is worth having and the whole reason it needs
    its own rows: every check the Gaussian families get for free through
    :attr:`VIEngine.guide_loc` has to be made some other way.

    The way is the library's own ELBO. ``vi_elbo_trace`` is
    ``E_q[log p - log q]`` as the *library* computed it, from its own guide
    object and its own draws; the stored groups give the same expectation
    from nothing but ``log_prior + log_likelihood`` and
    ``proposal_log_density``, through the driver's own code path. The two
    agreeing is a check of ``proposal_log_density``'s absolute value -- not
    merely of its variation, which ``ptp > 0`` covers, and not merely up to a
    constant, which an importance weight would be blind to.
    """

    def test_it_refuses_a_one_dimensional_problem_by_name(self, kit: Kit) -> None:
        engine = VIEngine(conjugate_problem(kit))
        with pytest.raises(EngineError, match="autoregressive"):
            engine.run(draws=10, steps=10, guide="flow")

    def test_the_attrs_name_the_guide_and_its_approximation(self, flow_run: Any) -> None:
        attrs = flow_run.attrs
        assert attrs["ampere_vi_guide"] == "flow"
        assert attrs["ampere_vi_guide_class"] == "AutoIAFNormal"
        assert attrs["ampere_approximation"] == "normalising_flow"

    def test_it_keeps_no_fitted_parameters_and_says_so(self, kit: Kit) -> None:
        """The contract W5.14 adds: all three attributes stay ``None``."""
        engine = VIEngine(correlated_problem(kit))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            run = engine.run(
                draws=50,
                steps=FLOW_STEPS,
                guide="flow",
                learning_rate=FLOW_LEARNING_RATE,
            )
        assert engine.guide_loc is None
        assert engine.guide_scale is None
        assert engine.guide_scale_tril is None
        assert run.attrs["ampere_vi_guide_parameters"] == "none"

    def test_it_finds_the_correlated_posterior(self, kit: Kit, flow_run: Any) -> None:
        """A flow that had not reached the mass would fail here, loudly."""
        reference = EmceeEngine(correlated_problem(kit), walkers=16).run(steps=800, burn_in=200)
        for name in ("model.norm", "model.index"):
            values = np.asarray(reference["posterior"][name]).ravel()
            drawn = np.asarray(flow_run["posterior"][name]).ravel()
            assert float(drawn.mean()) == pytest.approx(
                float(values.mean()), abs=float(values.std())
            )

    def test_the_stored_density_is_the_guides_own(self, flow_run: Any) -> None:
        """``ptp > 0``, and the absolute value against the library's ELBO."""
        stats = flow_run["sample_stats"].dataset
        proposal = np.asarray(stats["proposal_log_density"]).ravel()
        assert np.ptp(proposal) > 0.0

        log_prior = np.asarray(stats["log_prior"]).ravel()
        log_likelihood = np.asarray(stats["log_likelihood"]).ravel()
        # E_q[log p - log q]. The Jacobian correction is the same on both
        # terms -- `proposal_log_density` was moved into the constrained
        # coordinates the stored log_prior/log_likelihood live in -- so it
        # cancels in the difference and the estimate is the unconstrained
        # ELBO the library reports.
        estimate = float(np.mean(log_prior + log_likelihood - proposal))
        trace = json.loads(flow_run.attrs["ampere_vi_elbo_trace"])
        tail = float(np.mean(np.asarray(trace[-40:], dtype=float)))
        assert estimate == pytest.approx(tail, abs=2.0)


# ---------------------------------------------------------------------------
# 8. W5.14's battery: every guide SBC-ranked
# ---------------------------------------------------------------------------

#: A deliberately small simulating problem: SBC costs one full fit per
#: simulation, so it must be cheap rather than precise.
SBC_GRID = np.geomspace(1.0, 10.0, 8)
#: The reference wavelength at the grid's geometric centre, which is what
#: makes ``norm`` and ``index`` nearly orthogonal rather than the strongly
#: traded-off pair ``correlated_problem`` deliberately builds. The choice is
#: the difference between asking "are these guides calibrated?" and asking
#: "does a mean-field guide underestimate a correlation?" -- the second
#: question has a known answer (section 2 above measures it) and would make
#: the mean-field row fail by construction.
SBC_REFERENCE = float(np.sqrt(SBC_GRID[0] * SBC_GRID[-1]))
SBC_DATA = noisy(SBC_GRID, power_law(SBC_GRID, 2.0, INDEX), SIGMA, seed=11)

#: Per-PR budget: below ``sbc``'s own goodness-of-fit power floor, and it says
#: so. ``-m engines_full`` runs the budget that has power, exactly as
#: ``test_nested.py`` does for the nested samplers.
#:
#: Four rather than the nested battery's eight, because the cost of a
#: *refit* is not the same on the two routes and this was measured rather
#: than assumed: a replica fit is a few seconds on the pyro route and some
#: twenty on the numpyro one, where every replica problem is a fresh
#: realisation and therefore a fresh jax compilation, which no budget inside
#: the fit can shorten. Four simulations per guide over four guides is what
#: keeps this section's share of the jax leg in minutes rather than a
#: quarter of an hour, and at either count the row's claim is the same one:
#: that ranks come out, well shaped, with nothing failed.
SBC_COUNT = 4
SBC_FULL_COUNT = 100
SBC_DRAWS = 30

#: What each guide is run at during calibration. The Gaussian families need
#: only a short optimisation on a two-parameter problem; the flow needs the
#: longer one section 7 uses, for the reason given there.
SBC_RUN_OPTIONS: dict[str, dict[str, Any]] = {
    "normal": {"steps": 600},
    "multivariate": {"steps": 600},
    "laplace": {"steps": 600},
    "flow": {"steps": FLOW_STEPS, "learning_rate": FLOW_LEARNING_RATE},
}

#: The per-PR row's override, and the other half of this section's wall
#: time. The reduced row asserts only that ranks come out well shaped with
#: nothing failed -- a claim about the machinery, not about the fit -- which
#: a short optimisation shows exactly as well. The ``engines_full`` row,
#: which *is* a claim about the fit, uses the converged budgets above.
SBC_REDUCED_OVERRIDES: dict[str, dict[str, Any]] = {
    "normal": {"steps": 250},
    "multivariate": {"steps": 250},
    "laplace": {"steps": 250},
    "flow": {"steps": 600},
}

#: Draws per fit. ``sbc``'s own ``draws`` is how many the rank is taken
#: against, thinned out of what the fit produced, so the fit has to produce
#: at least that many -- so this is that many with a little room, and no
#: more: a draw is cheap for a Gaussian guide and not for a flow (see
#: ``FLOW_DRAWS``), and nothing here ranks against more than ``SBC_DRAWS``.
SBC_FIT_DRAWS = 40
GUIDES = tuple(SBC_RUN_OPTIONS)


def calibration_problem(kit: Kit, seed: int | None = SEED) -> FittingProblem:
    """The simulating problem the SBC rows draw truths and data from."""
    return FittingProblem(
        kit.module.PowerLaw(
            SBC_GRID,
            norm=st.norm(2.0, 0.5),
            index=st.norm(-1.0, 0.3),
            reference_wavelength=SBC_REFERENCE,
        ),
        [Dataset(SBC_DATA, likelihood=kit.likelihood())],
        seed=seed,
    )


def calibrate(kit: Kit, guide: str, *, count: int) -> Any:
    options = dict(SBC_RUN_OPTIONS[guide])
    if count <= SBC_COUNT:
        options.update(SBC_REDUCED_OVERRIDES.get(guide, {}))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return sbc(
            calibration_problem(kit),
            VIEngine,
            count=count,
            draws=SBC_DRAWS,
            run_options={"draws": SBC_FIT_DRAWS, "guide": guide, **options},
            seed=515151,
            label=f"W5.14 guide battery, {guide}",
        )


@pytest.fixture(scope="module")
def calibrations(kit: Kit) -> dict[str, Any]:
    """One SBC study per guide family, computed once for the whole module."""
    return {guide: calibrate(kit, guide, count=SBC_COUNT) for guide in GUIDES}


class TestEveryGuideIsRanked:
    """W5.14's battery row: "every new engine SBC-ranked through ``calibration.sbc``".

    A guide family is not an engine, but it is what an approximate run is an
    approximation *by*, so the battery is applied per family rather than once
    to the driver: a study that ranked ``normal`` only would say nothing about
    the Hessian ``laplace`` takes or the network ``flow`` trains, which are
    the two pieces of new machinery this item adds.
    """

    @pytest.mark.parametrize("guide", GUIDES)
    def test_the_ranks_are_produced_and_shaped(
        self, guide: str, calibrations: dict[str, Any]
    ) -> None:
        calibration = calibrations[guide]
        ranks = np.asarray(calibration["ranks"].values)
        assert ranks.shape == (SBC_COUNT, 2)
        assert ranks.min() >= 0
        assert ranks.max() <= SBC_DRAWS
        assert sorted(str(name) for name in calibration["parameter"].values) == [
            "model.index",
            "model.norm",
        ]
        assert "VIEngine" in calibration.attrs["ampere_calibration_engine"]

    @pytest.mark.parametrize("guide", GUIDES)
    def test_nothing_failed_and_the_uniformity_test_ran(
        self, guide: str, calibrations: dict[str, Any]
    ) -> None:
        calibration = calibrations[guide]
        assert calibration.attrs["ampere_calibration_failures"] == 0
        pvalues = np.asarray(calibration["ks_pvalue"].values)
        assert pvalues.shape == (2,)
        assert np.all((pvalues >= 0.0) & (pvalues <= 1.0))

    @pytest.mark.engines_full
    @pytest.mark.parametrize("guide", GUIDES)
    def test_the_ranks_are_uniform_at_a_budget_with_power(self, kit: Kit, guide: str) -> None:
        """The row the per-PR budget cannot make: Talts et al.'s actual test.

        It is a fair test of *these* four families only because the simulating
        problem's posterior is nearly Gaussian and nearly uncorrelated (see
        ``SBC_REFERENCE``): each family can represent it, so a failure here is
        a failure of the driver rather than the known and separately measured
        limitation of a mean-field guide.
        """
        calibration = calibrate(kit, guide, count=SBC_FULL_COUNT)
        pvalues = np.asarray(calibration["ks_pvalue"].values)
        assert np.all(pvalues > 0.01), f"{guide}'s ranks are not uniform: KS p-values {pvalues}"
