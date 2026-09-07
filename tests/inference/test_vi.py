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
   library, a non-differentiable problem, an unknown guide family.

Parametrised over the backends installed here that this driver supports, for
the reason ``test_nuts.py`` gives: writing the claim once per backend by hand
is how two backends drift apart.

Budgets are small and seeds fixed, so this belongs in the per-PR gate. VI is
cheap — a fit is a few thousand cheap gradient steps, not a chain — which is
part of why it exists.
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
    FittingProblem,
    GaussianFamily,
    Likelihood,
    Spectrum,
)
from ampere.inference import VIEngine
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
        """
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
        engine = VIEngine(conjugate_problem(kit))
        with pytest.raises(EngineError, match="guide family"):
            engine.run(draws=10, guide="laplace")

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
