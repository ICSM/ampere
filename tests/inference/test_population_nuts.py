"""NUTS fits a population jointly (**W5.12**).

``hierarchical_population.md``'s design horizon (a): the *joint* fit, as
opposed to horizon (b)'s importance reweighting of archived single-object fits
(``ampere.results.population``, W5.13). §11 Q5 left both routes open and
deferred the choice to Phase 5; this module is the joint route's evidence that
it works — one ``Population`` declaration, fifty member components, lowered to
a plate by the backend's realisation and sampled by NUTS through
``ampere.core.realise``.

The claim is a recovery claim, and it is the one the work item names: the
population hyperparameters of a fifty-member synthetic sample lie inside the
central 95 % of the posterior. Fifty members is the point — it is past where
the tie-based pattern of ``parameters.md`` §9 is comfortable (the flat layout
refuses above ``MAX_FLAT_MEMBERS``), and it is the scale the plate exists for.

Every member's likelihood is informative (three points at a tenth of the
population's own scatter), which is deliberate: a centred hierarchical
parameterisation funnels when the data say little about each θ_i, and a row
that failed for *that* reason would be testing Neal's funnel rather than
ampere's lowering. The non-centred reparameterisation is a user-level
declaration and a separate question.
"""

from __future__ import annotations

import dataclasses
import importlib
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
    HierarchicalPrior,
    Instrument,
    Likelihood,
    Parameter,
    Population,
    Spectrum,
)
from ampere.inference import NUTSEngine

#: The synthetic population. Fifty members, as the work item asks.
MEMBERS = 50

#: Truth: the hyperparameters NUTS has to find.
TRUE_MU = -1.30
TRUE_SIGMA = 0.35

#: Each member is observed at three wavelengths with a tenth of the
#: population's scatter, so its own index is well determined and the posterior
#: geometry is benign (see the module docstring).
GRID = np.array([1.0, 3.0, 9.0])
NOISE = 0.02
NORM = 2.0

SEED = 20260919


@dataclasses.dataclass(frozen=True)
class Kit:
    """One differentiable backend's pieces, so no test body names a library."""

    name: str
    module: Any


def _installed_kits() -> list[Kit]:
    found: list[Kit] = []
    for name in ("jax", "torch"):
        try:
            module = importlib.import_module(f"ampere.backends.{name}")
        except ImportError:  # the extra is not installed in this environment
            continue
        if name == "jax":
            # ``lowering.md`` §10.2(a): the application turns x64 on, not ampere.
            module.configure_x64()
        found.append(Kit(name, module))
    return found


KITS = _installed_kits()

pytestmark = pytest.mark.skipif(
    not KITS,
    reason="no differentiable backend installed; a joint population fit needs a realisation",
)


@pytest.fixture(scope="module", params=[kit.name for kit in KITS])
def kit(request: Any) -> Kit:
    return next(found for found in KITS if found.name == request.param)


def synthetic_indices() -> np.ndarray:
    """The fifty spectral indices the members were generated with.

    Drawn from ``Normal(TRUE_MU, TRUE_SIGMA)`` and then **standardised back
    onto those two numbers**, so the sample's own mean is ``TRUE_MU`` and its
    own scatter is ``TRUE_SIGMA`` exactly. That is deliberate, and it is what
    makes this row a test of the fit rather than of sampling luck: with fifty
    well-measured members the posterior on ``mu`` is about
    ``TRUE_SIGMA / sqrt(50) = 0.05`` wide and concentrates on the **sample**
    mean, so a draw whose sample mean happens to sit two standard errors from
    the parent would fail a 95 % interval check while the estimator was
    working perfectly. Standardising removes that coin-flip without weakening
    the claim: the fifty values still scatter exactly as a population does,
    and recovering their mean and spread is precisely what a joint population
    fit is for.
    """
    raw = np.random.default_rng(SEED).normal(0.0, 1.0, MEMBERS)
    return TRUE_MU + TRUE_SIGMA * (raw - raw.mean()) / raw.std()


def population_problem(kit: Kit) -> FittingProblem:
    """Fifty power-law members, declared as one population.

    ``norm`` is held fixed on every member: the population is about the
    *index*, and fifty nuisance amplitudes would double the sampler's
    dimension to say nothing extra about it.
    """
    indices = synthetic_indices()
    rng = np.random.default_rng(SEED + 1)
    models: dict[str, Any] = {}
    datasets: list[Dataset] = []
    for index, value in enumerate(indices):
        label = f"obj{index}"
        truth = NORM * GRID**value
        observed = Spectrum(
            GRID * u.micron,
            (truth + rng.normal(0.0, NOISE, GRID.size)) * u.Jy,
            uncertainty=np.full(GRID.size, NOISE) * u.Jy,
        )
        models[label] = kit.module.PowerLaw(
            GRID,
            norm=Parameter("norm", value=NORM, fixed=True),
            index=st.norm(TRUE_MU, 1.0),
        )
        datasets.append(
            Dataset(
                observed,
                Instrument([], channel="default", input_kind=Spectrum, label=f"scope{index}"),
                Likelihood(GaussianFamily(), kit.module.IndependentNoise()),
                model=label,
                label=f"d{index}",
            )
        )
    declaration = Population(
        "objects",
        members=[Parameter("index", HierarchicalPrior("norm", {"loc": "mu", "scale": "sigma"}))],
        hyperpriors=[
            Parameter("mu", st.norm(-1.0, 1.0)),
            Parameter("sigma", st.halfnorm(0.0, 1.0)),
        ],
        over=[f"obj{index}" for index in range(MEMBERS)],
    )
    return FittingProblem(models, datasets, populations=[declaration], seed=SEED)


@pytest.fixture(scope="module")
def problem(kit: Kit) -> FittingProblem:
    return population_problem(kit)


@pytest.fixture(scope="module")
def run(problem: FittingProblem) -> Any:
    with warnings.catch_warnings():
        # A hierarchical posterior sampled in its centred parameterisation
        # records the occasional divergence, and the driver warns about it; the
        # recovery claim is what this module asserts.
        warnings.simplefilter("ignore")
        return NUTSEngine(problem).run(draws=300, warmup=300, chains=1)


class TestTheDeclaration:
    """The problem NUTS is handed: one plate, fifty element bindings."""

    def test_it_is_one_array_valued_site_over_fifty_members(self, problem: FittingProblem) -> None:
        plated = problem.parameters["objects.index"]
        assert plated.shape == (MEMBERS,)
        assert plated.plate == "objects"
        # Two hyperparameters and one plate of fifty: fifty-two dimensions, not
        # fifty-two parameters. That is what the plate buys.
        assert problem.free_size == MEMBERS + 2
        assert len(problem.parameters.free_names) == 3

    def test_every_member_is_addressed_by_its_own_index(self, problem: FittingProblem) -> None:
        elements = {
            binding.component: binding.index
            for binding in problem.mapping.bindings
            if binding.index is not None
        }
        assert elements == {f"obj{index}": index for index in range(MEMBERS)}


class TestRecovery:
    """The claim: the truth is inside the posterior's central 95 %."""

    @staticmethod
    def draws(run: Any, name: str) -> np.ndarray:
        return np.asarray(run.posterior[name].values).reshape(-1)

    @pytest.mark.parametrize(
        ("name", "truth"), [("objects.mu", TRUE_MU), ("objects.sigma", TRUE_SIGMA)]
    )
    def test_the_hyperparameter_is_inside_the_central_95_percent(
        self, run: Any, name: str, truth: float
    ) -> None:
        samples = self.draws(run, name)
        lower, upper = np.quantile(samples, [0.025, 0.975])
        assert lower <= truth <= upper, (
            f"{name}: truth {truth} outside the central 95 % [{lower}, {upper}] "
            f"of {samples.size} draws (mean {samples.mean()})"
        )

    def test_the_posterior_concentrates_where_the_members_are(self, run: Any) -> None:
        """A sharper statement than the interval.

        The failure it catches is the plausible one: a population fitted
        against the *hyperprior's* own location (``Normal(-1.0, 1.0)``) rather
        than against its members would sit between the two, and 0.3 is six
        standard errors away.
        """
        standard_error = TRUE_SIGMA / np.sqrt(MEMBERS)
        assert self.draws(run, "objects.mu").mean() == pytest.approx(
            TRUE_MU, abs=3.0 * standard_error
        )

    def test_the_member_draws_recover_their_own_indices(self, run: Any) -> None:
        """The plate is fifty *different* draws, not one shared value.

        A population that had silently collapsed to a single θ — the failure
        ``hierarchical_population.md`` §4 warns ``Tie(..., prior=...)`` would
        produce, "silently and plausibly" — would pass every hyperparameter row
        above and fail here.
        """
        drawn = np.asarray(run.posterior["objects.index"].values)
        posterior_means = drawn.reshape(-1, MEMBERS).mean(axis=0)
        assert posterior_means.std() > 0.5 * TRUE_SIGMA
        assert np.corrcoef(posterior_means, synthetic_indices())[0, 1] > 0.9
