"""M2 over many spectra: the per-dataset GP amplitude as a population (**W7.2**).

The claim, stated so that it can fail. Eight spectra of one object class, the
``strong_sharp`` line injected in three of them and one clean member at five
times the noise (:mod:`examples.m2_populations.populations`), fitted three
ways — independent amplitudes, one tied amplitude, and the design memo's §4
non-centred population over ``d*.likelihood``:

1. **The declaration is the memo's.** ``z`` is internal and routed nowhere,
   ``gp.amplitude`` is derived and reaches spectrum *i*'s noise model as
   element *i*; the three arms differ in the amplitude declaration and nothing
   else (48, 41 and 50 free dimensions over the same data). Reference backend,
   no sampling: the per-PR rows every environment runs.
2. **The population stays honest.** Its model parameters stay within
   :data:`~examples.m2_misspecification.study.FLEXIBLE_MAX_BIAS_WIDTHS` on the
   deviated spectra (M2's own claim, re-made per member) and within
   :data:`POPULATION_MAX_BIAS_WIDTHS` on every spectrum.
3. **It localises.** GP localisation under the population peaks at the
   injected line on every deviated spectrum, and its score there exceeds every
   clean spectrum's by :data:`~examples.m2_misspecification.study.LOCALISATION_CONTRAST`.
4. **It recovers the spread.** ``gp.sigma``'s central 95 % contains the sample
   spread of the members' log-amplitudes (the mapping:
   :func:`~examples.m2_populations.populations.sample_spread`).
5. **It shrinks the member that needs it, and only that one.** The poorly
   constrained member's amplitude is narrower under the population than under
   independent amplitudes by :data:`POPULATION_SHRINKAGE_FACTOR`, while the
   deviated members' amplitudes stay above :data:`POPULATION_SIGNAL_FLOOR`.
6. **It is calibrated.** SBC over refits of the population arm: the
   hyperparameters' coverage at the study's ``NOMINAL`` level does not fall
   more than ``COVERAGE_ALLOWANCE`` below it, and their ranks pass the
   uniformity check — beside the independent arm's per-member coverage on the
   same replicas, the contrast.

Engines and budgets
-------------------
Claims 2-6 sample a fifty-dimensional problem. On the reference path one
evaluation costs about 6.5 ms (eight GP solves and fifty scipy prior
densities), so emcee — which needs more than a hundred walkers here — costs
about twenty minutes per arm at a budget worth asserting on; NUTS on jax costs
about three and a half minutes for all three arms at the study's own
:data:`~examples.m2_misspecification.study.TEST_NUTS`. So the science rows are
the **gradient path's**, parametrised over torch and jax as
``test_backend_agreement.py``'s are, and marked ``heavy`` (out of
``test-fast``; ``test-all`` runs them in the torch and jax environments, which
is where CI's per-push gate runs them). The reference path's emcee is the
``m2_full`` row, at :data:`~examples.m2_populations.populations.POPULATION_EMCEE`.
The module is one xdist group, so a ``-n`` run builds each fixture once.
"""

from __future__ import annotations

import warnings
from typing import Any

import numpy as np
import pytest

from examples.m2_misspecification import study
from examples.m2_misspecification.generators import scenario_named
from examples.m2_populations import populations
from examples.m2_populations.populations import (
    ARMS,
    DEVIATED,
    N_SPECTRA,
    POORLY_CONSTRAINED,
    POOR_NOISE_FRACTION,
    POPULATION_LABEL,
    dataset_label,
)

pytestmark = pytest.mark.xdist_group("m2-populations")

NATIVE = ["torch", "jax"]

# ---------------------------------------------------------------------------
# The pinned margins. Every one is measured at POPULATION_NUTS on jax over three
# data seeds (spectrum 0 at 42, 142 and 242; the run seed the study's), and set
# a third below the worst of the three, the study's convention.
# ---------------------------------------------------------------------------

#: Every model parameter of every spectrum under the population, in posterior
#: widths. Not the study's 1.5: that was measured on one spectrum, and with
#: thirty-two parameters the largest |median - truth| / width of a *calibrated*
#: posterior is itself about 2 (the expected maximum of 32 half-normal draws is
#: 2.3). The independent arm reaches the same numbers on the same spectra —
#: they are the noise realisations', not the declaration's. Measured worst:
#: 1.87, 2.57 and 1.92 (the independent arm 1.76, 2.63 and 2.07, on the same
#: spectra); the bound sits a third above the worst.
POPULATION_MAX_BIAS_WIDTHS = 3.5

#: The poorly constrained member's amplitude 68 % half-width, independent over
#: population. Measured: 1.45, 1.60 and 1.72 (0.0080 Jy against 0.0055 at the
#: per-PR seed). The margin is a third of the worst one's **excess over one** —
#: ``1 + 0.45 x 2/3`` — rather than a third of the ratio itself, which would be
#: 0.97 and assert that the population may *widen* the member: the claim is
#: the narrowing, and the narrowing is the excess.
POPULATION_SHRINKAGE_FACTOR = 1.3

#: ... and the deviated members' amplitudes must survive it, Jy. A population
#: that shrank every member towards its centre would satisfy the factor above
#: and be useless. Measured: the smallest deviated member's median under the
#: population is 0.0180, 0.0177 and 0.0176 Jy (0.0185-0.0213 under independent
#: amplitudes), and the largest *clean* member's is 0.0036, 0.0037 and 0.0034
#: (the poorly constrained one). The floor sits a third below the deviated
#: worst and three times above every clean member, so it separates the two.
POPULATION_SIGNAL_FLOOR = 0.011

#: The tied arm's cost, pinned where it is unambiguous: its one amplitude over
#: the largest clean (1 %-noise) member's independent amplitude. Measured: 9.6,
#: 5.7 and 5.1 (0.016 Jy tied, against at most 0.0031 for a clean spectrum on
#: its own); the factor sits a third below the worst. The deviated side —
#: whether the tie's too-small amplitude biases the deviated spectra — had the
#: same sign at all three seeds (the tied arm's mean worst bias over the three
#: deviated spectra 0.98, 0.96 and 1.15 widths against the population's 0.84,
#: 0.80 and 0.97) but a difference of 0.14-0.18 widths, below the study's
#: measured median scatter at this budget (0.31), so it is reported, not pinned.
TIED_INFLATION_FACTOR = 3.0


def _quiet_run(backend: str) -> populations.Comparison:
    with warnings.catch_warnings():
        # NUTS' own chatter (divergence counts, the realisation's notices) is
        # not what these rows assert on.
        warnings.simplefilter("ignore")
        return populations.compare(backend=backend, budget=populations.POPULATION_NUTS)


@pytest.fixture(scope="module")
def comparisons() -> dict[str, populations.Comparison]:
    """One three-arm comparison per native backend, built on first use."""
    return {}


def _comparison(cache: dict[str, populations.Comparison], backend: str) -> populations.Comparison:
    pytest.importorskip(backend)
    if backend not in cache:
        cache[backend] = _quiet_run(backend)
    return cache[backend]


# ---------------------------------------------------------------------------
# 1. The declaration (reference backend, no sampling: every environment)
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def data() -> tuple[Any, ...]:
    return populations.generate_population()


class TestTheDeclaration:
    def test_the_collection_is_the_scenario(self, data: tuple[Any, ...]) -> None:
        assert len(data) == N_SPECTRA
        deviated = [i for i, spectrum in enumerate(data) if spectrum.scenario.amplitude]
        assert deviated == list(DEVIATED)
        noisy = data[POORLY_CONSTRAINED]
        np.testing.assert_allclose(noisy.uncertainty, POOR_NOISE_FRACTION * np.abs(noisy.truth))
        # Distinct noise realisations: one object class, eight observations.
        assert len({float(spectrum.observed[0]) for spectrum in data}) == N_SPECTRA

    def test_the_arms_differ_only_in_the_amplitude(self, data: tuple[Any, ...]) -> None:
        problems = {arm: populations.build_arm(data, arm) for arm in ARMS}
        assert {arm: p.free_size for arm, p in problems.items()} == {
            "independent": 8 * 6,
            "tied": 8 * 5 + 1,
            "population": 8 * 5 + N_SPECTRA + 2,
        }
        shared = {
            name
            for name in problems["independent"].parameters.free_names
            if not name.endswith(".amplitude")
        }
        for problem in problems.values():
            assert shared <= set(problem.parameters.free_names)

    def test_the_population_is_the_memos_declaration(self, data: tuple[Any, ...]) -> None:
        problem = populations.build_arm(data, "population")
        free = problem.parameters.free_names
        assert {f"{POPULATION_LABEL}.{n}" for n in ("mu", "sigma", "z")} <= set(free)
        assert not any(name.endswith("likelihood.amplitude") for name in free)
        assert problem.parameters.derived_names == (f"{POPULATION_LABEL}.amplitude",)
        routed = {
            (b.global_name, b.index): (b.component, b.local_name)
            for b in problem.mapping.bindings
            if b.index is not None
        }
        # z is internal and routed nowhere; the derived amplitude reaches
        # spectrum i's noise model as element i.
        assert routed == {
            (f"{POPULATION_LABEL}.amplitude", i): (dataset_label(i), "likelihood.amplitude")
            for i in range(N_SPECTRA)
        }

    def test_the_population_scores_what_the_independent_arm_scores(
        self, data: tuple[Any, ...]
    ) -> None:
        """At matched amplitudes the two arms' likelihoods are the same number.

        The routing claim made numerically: element *i* of the derived
        amplitude is what spectrum *i*'s GP sees, so setting the independent
        arm's eight amplitudes to the population's eight derived values gives
        the same log-likelihood to rounding.
        """
        population = populations.build_arm(data, "population")
        independent = populations.build_arm(data, "independent")
        values = dict(
            population.parameters.complete(population.sample_prior(np.random.default_rng(7)))
        )
        amplitudes = np.asarray(values[f"{POPULATION_LABEL}.amplitude"], dtype=float)
        matched = {
            name: value
            for name, value in values.items()
            if not name.startswith(f"{POPULATION_LABEL}.")
        }
        matched.update(
            {f"{dataset_label(i)}.likelihood.amplitude": amplitudes[i] for i in range(N_SPECTRA)}
        )
        left = population.evaluate(values).log_likelihood
        right = independent.evaluate(matched).log_likelihood
        assert left == pytest.approx(right, rel=1e-10)

    def test_the_tied_arm_has_one_amplitude(self, data: tuple[Any, ...]) -> None:
        problem = populations.build_arm(data, "tied")
        amplitudes = [n for n in problem.parameters.free_names if n.endswith("amplitude")]
        assert amplitudes == [populations.TIED_NAME]


# ---------------------------------------------------------------------------
# 2-5. The science (NUTS on torch and jax; heavy)
# ---------------------------------------------------------------------------


@pytest.mark.heavy
@pytest.mark.parametrize("backend", NATIVE)
class TestThePopulationOnTheM2Pattern:
    def test_the_deviated_spectra_stay_inside_the_flexible_bound(
        self, comparisons: Any, backend: str
    ) -> None:
        """M2's claim per deviated member. Measured worst: 0.93, 0.99 and 1.17."""
        comparison = _comparison(comparisons, backend)
        worst = {i: comparison.worst_bias("population", i) for i in comparison.deviated}
        print(f"\n{backend} population, deviated worst bias: {worst}")
        assert max(worst.values()) <= study.FLEXIBLE_MAX_BIAS_WIDTHS, worst

    def test_every_spectrum_stays_inside_the_population_bound(
        self, comparisons: Any, backend: str
    ) -> None:
        comparison = _comparison(comparisons, backend)
        worst = {i: comparison.worst_bias("population", i) for i in range(N_SPECTRA)}
        print(f"\n{backend} population, worst bias per spectrum: {worst}")
        assert max(worst.values()) <= POPULATION_MAX_BIAS_WIDTHS, worst

    def test_the_gp_localises_the_line_on_the_deviated_spectra(
        self, comparisons: Any, backend: str
    ) -> None:
        """Peak within the study's tolerance, and far above every clean spectrum.

        Measured: every deviated peak within 0.7 grid spacings of the line, and
        the weakest deviated score over the strongest clean one 29.9, 14.9 and
        13.0, against the study's contrast of 5.
        """
        comparison = _comparison(comparisons, backend)
        found = populations.localise(comparison.fits["population"])
        centre = scenario_named(populations.DEVIATION).localised_at
        assert centre is not None
        spacing = float(np.diff(comparison.data[0].wavelength)[0])
        clean = max(found[dataset_label(i)].largest for i in comparison.clean)
        for index in comparison.deviated:
            hit = found[dataset_label(index)]
            print(f"\n{backend} {hit.dataset}: peak {hit.peak:.5f} um, score {hit.largest:.2f}")
            assert abs(hit.peak - centre) <= study.LOCALISATION_TOLERANCE_POINTS * spacing
            assert hit.largest >= study.LOCALISATION_CONTRAST * clean

    def test_sigma_recovers_the_amplitude_spread(self, comparisons: Any, backend: str) -> None:
        """Measured: 1.33 in [0.95, 2.80], 1.09 in [0.73, 2.60], 1.11 in [0.80, 2.58]."""
        comparison = _comparison(comparisons, backend)
        spread = populations.sample_spread(comparison)
        low, high = np.quantile(comparison.hyperparameter("sigma"), [0.025, 0.975])
        print(f"\n{backend} spread {spread:.3f} in gp.sigma's central 95 % [{low:.3f}, {high:.3f}]")
        assert low <= spread <= high

    def test_the_poorly_constrained_member_is_shrunk(self, comparisons: Any, backend: str) -> None:
        comparison = _comparison(comparisons, backend)
        factor = comparison.shrinkage(POORLY_CONSTRAINED)
        print(f"\n{backend} shrinkage of {dataset_label(POORLY_CONSTRAINED)}: {factor:.2f}")
        assert factor >= POPULATION_SHRINKAGE_FACTOR

    def test_the_deviated_members_are_not_shrunk_away(self, comparisons: Any, backend: str) -> None:
        comparison = _comparison(comparisons, backend)
        medians = {i: comparison.amplitude("population", i).median for i in comparison.deviated}
        print(f"\n{backend} deviated amplitudes under the population: {medians}")
        assert min(medians.values()) >= POPULATION_SIGNAL_FLOOR

    def test_the_tie_inflates_the_clean_members(self, comparisons: Any, backend: str) -> None:
        """One amplitude for all: about ten times what a clean spectrum needs."""
        comparison = _comparison(comparisons, backend)
        tied = comparison.amplitude("tied", 0).median
        clean = max(
            comparison.amplitude("independent", i).median
            for i in comparison.clean
            if i != POORLY_CONSTRAINED
        )
        print(f"\n{backend} tied amplitude {tied:.4f} Jy, largest clean {clean:.4f} Jy")
        assert tied >= TIED_INFLATION_FACTOR * clean


# ---------------------------------------------------------------------------
# 6. Calibration (SBC over refits; NUTS on torch and jax; heavy)
# ---------------------------------------------------------------------------

#: The study's SBC thresholds (``test_many_lines_calibration.py``): the nominal
#: level coverage is read at, and how far below it a calibrated fit may fall at
#: ten refits — two binomial standard errors, rounded up.
NOMINAL = 0.9
COVERAGE_ALLOWANCE = 0.2

#: Refits at the per-PR level, and behind ``-m m2_full``.
SIMULATIONS = 10
FULL_SIMULATIONS = 24

#: The uniformity check's floor: Talts et al.'s ranks of a calibrated fit fail
#: a KS test at this level one time in a hundred per column.
UNIFORMITY_FLOOR = 0.01

HYPERPARAMETERS = (f"{POPULATION_LABEL}.mu", f"{POPULATION_LABEL}.sigma")


def _calibrate(arm: str, count: int, backend: str) -> Any:
    pytest.importorskip(backend)
    with warnings.catch_warnings():
        # sbc's own warning that ten simulations is a smoke-test count, and
        # NUTS' chatter: the coverage is the assertion.
        warnings.simplefilter("ignore")
        return populations.calibrate(arm, count=count, backend=backend)


@pytest.fixture(scope="module")
def calibrations() -> dict[tuple[str, str], Any]:
    return {}


def _calibration(cache: dict[tuple[str, str], Any], arm: str, backend: str) -> Any:
    if (arm, backend) not in cache:
        cache[(arm, backend)] = _calibrate(arm, SIMULATIONS, backend)
    return cache[(arm, backend)]


@pytest.mark.heavy
@pytest.mark.parametrize("backend", NATIVE)
class TestThePopulationIsCalibrated:
    """What is simulated and what is ranked: :func:`~examples.m2_populations.populations.calibrate`.

    Each simulation draws the population's hyperparameters, the members'
    ``z``, the length scales and the model parameters from the population
    arm's own prior and observes every spectrum as the model plus a Matérn-3/2
    realisation at ``exp(mu + sigma z_i)`` plus the white noise — the
    deviation, drawn per member at a population-drawn size. The truth ranked
    against is that draw.
    """

    def test_the_hyperparameters_cover_at_nominal(self, calibrations: Any, backend: str) -> None:
        """Measured on jax: 0.8 on both (ten refits, 236 s), against a floor of 0.7."""
        coverage = populations.coverage_at(
            _calibration(calibrations, "population", backend), NOMINAL
        )
        print(f"\n{backend} population coverage at {NOMINAL}: {coverage}")
        for name in HYPERPARAMETERS:
            assert coverage[name] >= NOMINAL - COVERAGE_ALLOWANCE, (name, coverage)

    def test_the_hyperparameters_ranks_are_uniform(self, calibrations: Any, backend: str) -> None:
        """Measured on jax: KS p-values 0.53 (``mu``) and 0.83 (``sigma``)."""
        calibration = _calibration(calibrations, "population", backend)
        names = [str(name) for name in calibration["parameter"].values]
        pvalues = dict(zip(names, np.asarray(calibration["ks_pvalue"].values), strict=True))
        print(f"\n{backend} uniformity p-values: {pvalues}")
        for name in HYPERPARAMETERS:
            assert pvalues[name] >= UNIFORMITY_FLOOR, (name, pvalues)

    def test_the_contrast_is_reported(self, calibrations: Any, backend: str) -> None:
        """The independent arm's per-member coverage on the same replicas.

        Reported beside the population's, and held only to having produced a
        usable fit on every replica: whether N independent half-normal
        amplitudes over- or under-cover a population-drawn truth depends on
        where the draw lands against the half-normal's scale, so the direction
        is not a claim this budget can pin. Measured on jax at ten refits: the
        population's per-member coverage 0.7, 0.7, 0.8, 0.8; the independent
        arm's 0.8, 0.7, 0.9, 0.9 — no contrast this count resolves, which is
        itself the finding: N independent half-normal amplitudes are not
        miscalibrated here, they are only less informed (the shrinkage rows).
        """
        population = _calibration(calibrations, "population", backend)
        independent = _calibration(calibrations, "independent", backend)
        print(
            f"\n{backend} member amplitudes' coverage at {NOMINAL}:"
            f"\n  population  {populations.coverage_at(population, NOMINAL)}"
            f"\n  independent {populations.coverage_at(independent, NOMINAL)}"
        )
        for calibration in (population, independent):
            assert int(calibration.attrs["ampere_calibration_failures"]) == 0
            assert calibration["ranks"].shape[0] == SIMULATIONS


# ---------------------------------------------------------------------------
# The marked budget: the reference path's emcee, and SBC at the full count
# ---------------------------------------------------------------------------


@pytest.mark.m2_full
@pytest.mark.study
def test_the_reference_path_reaches_the_same_conclusions() -> None:
    """emcee on the numpy path, the independent and population arms, at POPULATION_EMCEE.

    The claims the per-PR rows make on the gradient path, made again by the
    ensemble sampler on the reference backend: the spread recovered, the
    shrinkage factor, the floor. About forty minutes. Measured: M2FULL_PLACEHOLDER.
    """
    with warnings.catch_warnings():
        # EnsembleSizeWarning says exactly what the module docstring says.
        warnings.simplefilter("ignore")
        comparison = populations.compare(
            backend="reference",
            budget=populations.POPULATION_EMCEE,
            arms=("independent", "population"),
        )
    print("\n" + comparison.table())
    low, high = np.quantile(comparison.hyperparameter("sigma"), [0.025, 0.975])
    assert low <= populations.sample_spread(comparison) <= high
    assert comparison.shrinkage(POORLY_CONSTRAINED) >= POPULATION_SHRINKAGE_FACTOR
    for index in comparison.deviated:
        assert comparison.amplitude("population", index).median >= POPULATION_SIGNAL_FLOOR
        assert comparison.worst_bias("population", index) <= study.FLEXIBLE_MAX_BIAS_WIDTHS


@pytest.mark.m2_full
@pytest.mark.parametrize("backend", NATIVE)
def test_the_calibration_at_the_full_count(backend: str) -> None:
    """SBC at FULL_SIMULATIONS, where a coverage estimate has some resolution."""
    calibration = _calibrate("population", FULL_SIMULATIONS, backend)
    coverage = populations.coverage_at(calibration, NOMINAL)
    print(f"\n{backend} population coverage at {NOMINAL}, {FULL_SIMULATIONS} refits: {coverage}")
    for name in HYPERPARAMETERS:
        assert coverage[name] >= NOMINAL - COVERAGE_ALLOWANCE, (name, coverage)
