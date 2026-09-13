"""W4.2: the M2 question asked of visibilities, by simulation-based calibration.

M2's claim is that the flexible likelihood keeps the *physical* parameters
honest when the model is wrong. W4.2 implements the circular complex GP, which
is that likelihood on the Phase-4 proof modality, so the item's last acceptance
criterion is the M2 experiment in miniature: an interferometric binary observed
with a **correlated calibration error** nobody models, fitted twice — under
independent noise and under the GP — with simulation-based calibration (W3.6's
:func:`ampere.results.calibration.sbc`) deciding which set of credible
intervals can be believed.

What is injected, and why it is not a GP draw
---------------------------------------------
The systematic is a fixed, smooth, complex gain error across the ``(u, v)``
plane — a product of sinusoids on a scale of four mega-wavelengths, with a
per-component RMS about twice the thermal sigma. Deliberately **not** a draw
from the kernel the GP is then given: a systematic drawn from the fitted
covariance function would make the GP fit correctly specified, and "a correctly
specified model is calibrated" is a statement about SBC rather than about
ampere. A deterministic function the GP's family does not contain is the
harder case, and the one a calibration error actually is.

What is asserted, and what is merely reported
---------------------------------------------
The assertion is about **coverage**, and it is the directional claim plus a
margin, in the pattern ``tests/m2/test_fringing_kernel.py`` established: the
independent fit's central-90% intervals contain the truth *less* often than
they claim, the GP's contain it at least as often, and the gap between them is
wide enough not to be a coincidence of the budget.

Three things are reported rather than pinned, and the third is why the
uniformity p-value is not what this row asserts.

1. The budget is small — :data:`SIMULATIONS` simulations, where Talts et al.
   (2018) want a hundred or more — so the numbers are evidence of a gross
   effect and nothing finer. ``sbc`` warns about exactly that and the warning
   is allowed through.
2. The measured numbers at the pinned seed: independent 0.75 / 0.75, GP
   1.00 / 1.00, against a nominal 0.90.
3. **The GP over-covers.** Its intervals are conservative rather than exact at
   this amplitude prior, so a Kolmogorov—Smirnov test of its ranks against the
   uniform flags it too — in the opposite direction from the independent fit's
   under-coverage. Over-coverage is not a calibration failure in the sense
   that matters for a published error bar: an interval that contains the truth
   more often than it claims does not invite a false conclusion, and an
   interval that contains it less often does. So coverage is what is pinned
   and the p-value is not, and that choice is stated here rather than left to
   be inferred from which assertion happens to be present.

Cost: about eighty seconds on the reference backend, which is the same order as
the rest of this suite's session fixture. It earns it by being the only row
anywhere that runs the circular complex GP through a whole fit, many times, and
looks at what came out.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any, ClassVar

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.core import (
    ComplexGaussianFamily,
    Dataset,
    DatasetCollection,
    DenseGP,
    FittingProblem,
    GaussianProcessNoise,
    IndependentNoise,
    Instrument,
    Likelihood,
    Marginalisation,
    Matern32,
    Transformation,
    VisibilitySet,
)
from ampere.backends.reference.interferometry import BinaryVisibilities
from ampere.inference import EmceeEngine
from ampere.results.calibration import sbc

#: The observing wavelength, micron. Monochromatic: the chromatic case is
#: W4.4's, and the claim here is about (u, v) correlation.
WAVELENGTH = 2.2

#: A four-telescope array, metres east and north of the array centre. The same
#: stations ``tests/conformance/test_interferometry.py`` uses, so the two files
#: describe one instrument.
STATIONS = np.array([[0.0, 0.0], [-32.4, 39.8], [6.5, -30.1], [-14.7, -48.9]])

#: Six hour angles, which fills the uv plane rather than sampling it at six
#: points: two parameters need more than one snapshot to be well determined,
#: and a GP needs enough neighbouring baselines for "correlated" to mean
#: anything.
HOUR_ANGLES = np.linspace(-0.7, 0.7, 6)

#: The binary every simulation is drawn around.
BINARY = {
    "separation": 12.0,
    "position_angle": 0.7,
    "flux_ratio": 0.42,
    "flux": 1.7,
    "component_fwhm": 2.0,
}

#: Thermal noise, per component, Jy — the sigma a ``VisibilitySet`` carries.
SIGMA = 0.03

#: Amplitude of the injected gain error, Jy. About twice :data:`SIGMA` in RMS
#: once the sinusoids are averaged over the coverage, which is the regime the
#: question is interesting in: large enough to bias a rigid fit, small enough
#: that the binary is still there to be found.
GAIN_AMPLITUDE = 0.12

#: Scale of the injected error, in wavelengths. Comparable to the array's own
#: extent, so the error is smooth *across* the coverage rather than
#: sample-to-sample — which is what makes it invisible to independent noise and
#: visible to a kernel.
GAIN_SCALE = 4.0e6

#: How many simulations. See the module docstring on what this budget does and
#: does not buy.
SIMULATIONS = 16

#: Posterior draws each rank is taken against, and the emcee budget per fit.
RANK_DRAWS = 150
WALKERS = 16
STEPS = 300
BURN_IN = 150

#: The nominal level the assertion is made at.
NOMINAL = 0.9

#: How far apart the two fits' coverages must be. Informational in the sense
#: that the number is a judgement rather than a prediction: the measured gap is
#: 0.25 and one simulation is worth 1/16 = 0.0625, so 0.15 leaves about two
#: simulations of drift before the row fails.
COVERAGE_MARGIN = 0.15


def coverage() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(u, v, lambda)`` for every baseline at every hour angle, in wavelengths."""
    metres_per_wavelength = WAVELENGTH * 1e-6
    pairs = np.array(
        [
            (STATIONS[j] - STATIONS[i]) / metres_per_wavelength
            for i in range(len(STATIONS))
            for j in range(i + 1, len(STATIONS))
        ]
    )
    u_points = np.concatenate(
        [pairs[:, 0] * math.cos(h) - pairs[:, 1] * math.sin(h) for h in HOUR_ANGLES]
    )
    v_points = np.concatenate(
        [pairs[:, 0] * math.sin(h) + pairs[:, 1] * math.cos(h) for h in HOUR_ANGLES]
    )
    return u_points, v_points, np.full(u_points.size, WAVELENGTH)


U_POINTS, V_POINTS, WAVELENGTHS = coverage()
N_SAMPLES = U_POINTS.size


def injected_gain() -> np.ndarray:
    """The systematic: a deterministic, smooth complex gain across the plane.

    A product of sinusoids in ``u`` and ``v``, so it is smooth on
    :data:`GAIN_SCALE` and belongs to no kernel family — see the module
    docstring on why it is not a draw from the Matérn the GP is given.
    """
    x = U_POINTS / GAIN_SCALE
    y = V_POINTS / GAIN_SCALE
    real = np.cos(x + 0.3) * np.cos(y - 0.2)
    imaginary = np.sin(x - 0.5) * np.cos(y + 0.7)
    return GAIN_AMPLITUDE * (real + 1j * imaginary)


GAIN = injected_gain()


class CalibrationError(Transformation):
    """Add a fixed complex gain error to a ``VisibilitySet``.

    The misspecification, expressed where a misspecification belongs: an
    instrument step in the *simulating* problem's chain, absent from both fitted
    problems. Nothing here is a parameter — the error is a fact about the
    experiment, not something either fit is allowed to learn.
    """

    ACCEPTS: ClassVar[tuple[type, ...]] = (VisibilitySet,)
    PRODUCES: ClassVar[type] = VisibilitySet

    def __init__(self, gain: np.ndarray, *, label: str | None = None) -> None:
        super().__init__(label=label)
        # Two real buffers rather than one complex one: a buffer is declared
        # data and the container vocabulary is real-valued throughout.
        self.register_buffer("gain_real", np.real(gain))
        self.register_buffer("gain_imag", np.imag(gain))

    def apply(self, samples: Any, values: Mapping[str, Any]) -> VisibilitySet:
        context = self.context(values)
        gain = np.asarray(context["gain_real"]) + 1j * np.asarray(context["gain_imag"])
        return samples.with_values(np.asarray(samples.values) + gain)


def binary_model() -> BinaryVisibilities:
    """The analytic binary, with the two parameters the study ranks left free.

    The analytic visibility twin rather than an image plus ``FourierSample``,
    and the reason is cost: this row pays for thirty-two full fits, and a direct
    DFT of a grid per likelihood evaluation would put it out of reach of a gate.
    W4.1's conformance rows are what say the two routes agree.
    """
    return BinaryVisibilities(
        U_POINTS,
        V_POINTS,
        WAVELENGTHS * u.micron,
        channels="vis",
        separation=st.uniform(6.0, 14.0),
        position_angle=BINARY["position_angle"],
        flux_ratio=st.uniform(0.15, 0.7),
        flux=BINARY["flux"],
        component_fwhm=BINARY["component_fwhm"],
    )


def observed_template() -> VisibilitySet:
    """An empty container carrying the coverage and the thermal sigma."""
    return VisibilitySet(
        U_POINTS,
        V_POINTS,
        WAVELENGTHS * u.micron,
        np.zeros(N_SAMPLES, dtype=complex) * u.Jy,
        uncertainty=np.full(N_SAMPLES, SIGMA) * u.Jy,
    )


def simulating_problem() -> FittingProblem:
    """What the data are drawn from: the binary **plus** the calibration error."""
    return FittingProblem(
        binary_model(),
        DatasetCollection(
            {
                "vis": Dataset(
                    observed_template(),
                    Instrument(
                        [CalibrationError(GAIN, label="gain")],
                        channel="vis",
                        input_kind=VisibilitySet,
                        label="array",
                    ),
                    likelihood=Likelihood(ComplexGaussianFamily(), IndependentNoise()),
                    label="vis",
                )
            }
        ),
        seed=20260913,
    )


def circular_gp_noise() -> GaussianProcessNoise:
    """The flexible likelihood's noise model for this modality.

    ``axes=("u", "v")`` because it must be: a ``VisibilitySet``'s third axis is
    a wavelength in micron, and W4.5's per-leaf unit rule refuses one isotropic
    length-scale across mixed units. The length-scale prior spans two decades
    around :data:`GAIN_SCALE` and the amplitude prior two decades around the
    injected RMS — wide enough that the fit is finding the correlation rather
    than being told it.
    """
    return GaussianProcessNoise(
        Matern32(st.loguniform(1e-2, 0.3), st.loguniform(5e5, 5e7), axes=("u", "v")),
        DenseGP(),
    )


def fitted_problem(observed: VisibilitySet, *, flexible: bool) -> FittingProblem:
    """The binary alone, against *observed*, with or without the GP.

    Neither fit carries the calibration step, which is the whole point: both
    are misspecified, and the question is which one's credible intervals
    survive it.
    """
    noise: Any = circular_gp_noise() if flexible else IndependentNoise()
    return FittingProblem(
        binary_model(),
        DatasetCollection(
            {
                "vis": Dataset(
                    observed,
                    Instrument([], channel="vis", input_kind=VisibilitySet, label="array"),
                    likelihood=Likelihood(ComplexGaussianFamily(), noise),
                    label="vis",
                )
            }
        ),
        seed=20260913,
    )


def run_sbc(*, flexible: bool) -> Any:
    """``sbc`` over :data:`SIMULATIONS` refits of one of the two formulations.

    The factory reads the simulated container off the replica and builds its
    *own* problem, which is the route ``sbc``'s docstring describes for putting
    a statistic on trial: the simulating problem and the fitted one are
    different formulations of one experiment.
    """

    def factory(replica: FittingProblem) -> Any:
        observed = replica.datasets["vis"].observed
        return EmceeEngine(fitted_problem(observed, flexible=flexible), walkers=WALKERS)

    with pytest.warns(UserWarning, match="goodness-of-fit test"):
        return sbc(
            simulating_problem(),
            factory,
            count=SIMULATIONS,
            draws=RANK_DRAWS,
            run_options={"steps": STEPS, "burn_in": BURN_IN},
            parameters=["model.separation", "model.flux_ratio"],
            seed=4242,
            label=f"{'flexible' if flexible else 'rigid'} visibility fit",
        )


def coverage_at_nominal(calibration: Any) -> np.ndarray:
    """The empirical coverage of the central :data:`NOMINAL` interval, per parameter."""
    levels = np.asarray(calibration["level"].values)
    index = int(np.argmin(np.abs(levels - NOMINAL)))
    assert abs(float(levels[index]) - NOMINAL) < 1e-9, "the nominal level is not on the grid"
    return np.asarray(calibration["coverage"].values)[index]


class TestTheFlexibleLikelihoodComposesOnVisibilities:
    """The structural half: cheap, and it has to hold before the numbers mean anything."""

    def test_the_pair_marginalises_analytically(self) -> None:
        """``complex_gaussian`` + GP is ANALYTIC, so no latent block is declared.

        The alternative — a latent formulation — would put one whitened value
        per retained sample into the parameter vector and put this study out of
        reach of an emcee budget entirely. That the closed form exists is what
        makes the rest of this module affordable.
        """
        problem = fitted_problem(observed_template(), flexible=True)
        dataset = problem.datasets["vis"]
        assert dataset.likelihood.marginalisation is Marginalisation.ANALYTIC
        assert dataset.latent is None
        # Two binary parameters plus the kernel's two, and nothing else.
        assert problem.free_size == 4

    def test_the_kernel_acts_on_the_two_dimensionless_axes(self) -> None:
        """The selector is bound to columns of the coordinate block, not assumed."""
        problem = fitted_problem(observed_template(), flexible=True)
        noise = problem.datasets["vis"].likelihood.noise
        observed = problem.datasets["vis"].observed
        bound = noise.kernel_for(observed)
        assert bound.selected_axes([axis.name for axis in observed.axes]) == ("u", "v")

    def test_the_injected_error_dominates_the_thermal_noise(self) -> None:
        """The experiment is in the regime the question is about.

        If the systematic were small compared with sigma there would be nothing
        for either fit to get wrong, and a row asserting that the GP wins would
        be asserting noise.
        """
        assert np.std(np.real(GAIN)) > SIGMA
        assert np.std(np.imag(GAIN)) > SIGMA


class TestTheGpKeepsTheIntervalsHonest:
    """The M2 question on visibilities: which fit's credible intervals hold up.

    One SBC study per formulation, both at the same seed and the same budget,
    so the only difference between them is the noise model.
    """

    @pytest.fixture(scope="class")
    @classmethod
    def rigid(cls) -> Any:
        return run_sbc(flexible=False)

    @pytest.fixture(scope="class")
    @classmethod
    def flexible(cls) -> Any:
        return run_sbc(flexible=True)

    def test_the_rigid_fit_under_covers(self, rigid: Any) -> None:
        """Independent noise on a correlated error: the intervals are too narrow.

        Measured 0.75 on both parameters against a nominal 0.90. The assertion
        is "below nominal", which at this budget is ``<= 14/16``.
        """
        measured = coverage_at_nominal(rigid)
        print(f"\nrigid coverage at {NOMINAL}: {measured}")
        assert np.all(measured < NOMINAL)

    def test_the_flexible_fit_does_not(self, flexible: Any) -> None:
        """The GP's intervals contain the truth at least as often as they claim.

        Measured 1.00 on both parameters: conservative rather than exact, which
        the module docstring states and explains. What is asserted is that they
        are not *narrow*, because a credible interval that is too narrow is the
        failure a published error bar cannot survive.
        """
        measured = coverage_at_nominal(flexible)
        print(f"\nflexible coverage at {NOMINAL}: {measured}")
        assert np.all(measured >= NOMINAL)

    def test_the_gap_between_them_is_wide(self, rigid: Any, flexible: Any) -> None:
        """The pinned comparison, with :data:`COVERAGE_MARGIN` of slack.

        The margin rather than two absolute thresholds, for the reason W4.5's
        fringing row gives for pinning a margin: it is the difference the
        experiment is about, and it is the quantity that does not move when the
        budget does.
        """
        rigid_coverage = coverage_at_nominal(rigid)
        flexible_coverage = coverage_at_nominal(flexible)
        gap = float(np.min(flexible_coverage) - np.min(rigid_coverage))
        print(f"\ncoverage gap (flexible - rigid) at {NOMINAL}: {gap:.3f}")
        assert gap > COVERAGE_MARGIN

    def test_every_simulation_produced_a_usable_fit(self, rigid: Any, flexible: Any) -> None:
        """A row that silently ranked four of sixteen would still pass the above."""
        for calibration in (rigid, flexible):
            assert int(calibration.attrs["ampere_calibration_failures"]) == 0
            assert calibration["ranks"].shape == (SIMULATIONS, 2)
