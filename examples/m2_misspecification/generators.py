"""The misspecified data: four scenarios, one truth, any number of points.

Milestone M2's controlled experiment needs data whose *deviation from the
model* is known exactly, because the whole claim under test is about what a
likelihood does when the model is wrong. So the observed spectra here are
built in three explicit stages, and each stage is kept:

1. the **truth** — the fitting model itself, evaluated at
   :data:`~examples.m2_misspecification.model.TRUTH`. Nothing about it is
   unknown to the fit;
2. the **deviation** — a multiplicative factor ``1 + delta(lambda)`` the
   fitting model has no parameter for, and cannot produce at any point of its
   parameter space. This is the misspecification, and :class:`Scenario` is
   the only place it is defined;
3. the **noise** — Gaussian, 1 % of the *undeviated* truth, seeded.

Keeping all three (:class:`SyntheticSpectrum` holds each) is what lets a test
assert that the GP localised the deviation *where the deviation actually is*
rather than merely somewhere.

The four scenarios are the paper study's
(``examples/examples_paper/flexible_likelihood_comparison.py``), transcribed
rather than reinvented so the reproduction is one:

======================  ==========================================  =============
key                     deviation                                   character
======================  ==========================================  =============
``none``                nothing                                     control
``mild``                2.5 % sinusoidal fringing, 0.0028 µm period smooth, global
``strong_smooth``       7 % sinusoidal fringing, same period        smooth, global
``strong_sharp``        12 % Gaussian emission line at 0.8630 µm    sharp, local
======================  ==========================================  =============

The ``none`` row is not a formality. It is the control that stops the flexible
likelihood being praised for something it did not do: if the GP inflated the
posterior when there was nothing to absorb, the comparison against the standard
likelihood would be a comparison of a wide posterior with a narrow one rather
than of a wrong answer with a right one.

Parametrised by size
--------------------
:func:`wavelength_grid` spans the same physical range (a Gaia-RVS-like
0.842 to 0.872 µm) at any number of points, so the size ladder 200 → 2 000 →
20 000 samples *the same spectrum* more finely rather than observing more of
it. That is deliberate: it holds the physics, the deviation and the noise level
fixed so that a change in the answer along the ladder is a change in the
inference and not in the experiment. It also means the posterior narrows as
:math:`\\sqrt{N}`, which is exactly the regime in which an unmodelled 7 %
ripple is most damaging — a bias that was two standard deviations at 200 points
is twenty at 20 000.

Everything is seeded (:data:`SEED`, 42, the paper study's) and every array is
float64.
"""

from __future__ import annotations

import dataclasses
from typing import Any

import numpy as np

from ampere.core import Spectrum

from .model import FLUX_UNIT, TRUTH, WAVELENGTH_UNIT, flux_at

__all__ = [
    "NOISE_FRACTION",
    "SCENARIOS",
    "SEED",
    "SIZES",
    "WAVELENGTH_MAX",
    "WAVELENGTH_MIN",
    "Scenario",
    "SyntheticSpectrum",
    "generate",
    "scenario_named",
    "wavelength_grid",
]

#: The Gaia-RVS-like range the paper study works in, micron.
WAVELENGTH_MIN = 0.842
WAVELENGTH_MAX = 0.872

#: Per-point 1-sigma uncertainty, as a fraction of the **undeviated** truth.
#: Taking it from the truth rather than from the deviated flux is what keeps
#: the four scenarios' error bars identical, so a difference between two fits
#: cannot be a difference in how much they were told to trust the data.
NOISE_FRACTION = 0.01

#: The paper study's seed. Fixed here rather than passed around, so "the
#: 200-point strong-fringing spectrum" names exactly one array.
SEED = 42

#: The size ladder (``DEVELOPMENT_PLAN.md`` §5's M2 line: 10-100x the paper
#: study's data). The dense solver is only affordable at the first two.
SIZES: tuple[int, ...] = (200, 2_000, 20_000)


@dataclasses.dataclass(frozen=True)
class Scenario:
    """One controlled misspecification: a multiplicative deviation from truth.

    A scenario is *either* a sinusoidal fringe (``period`` set) or a Gaussian
    line (``centre`` and ``width`` set) or nothing at all (``amplitude`` zero),
    and :meth:`deviation` dispatches on which. Two shapes rather than an
    open-ended callable because these two are the ones the study compares — a
    smooth, global deviation a stationary kernel can absorb, and a sharp, local
    one it can only partly absorb — and a scenario the tests can reason about
    (:attr:`localised_at`) is worth more here than a general one.

    Attributes
    ----------
    key
        The scenario's name, used in run labels and figure panels.
    description
        Human-readable label for figures and tables.
    amplitude
        Peak fractional deviation: 0.025 is a 2.5 % ripple.
    period
        Fringe period in micron, for a sinusoidal deviation.
    phase
        Phase offset of the sinusoid, radians.
    centre, width
        Centre and Gaussian sigma in micron, for a line-shaped deviation.
    """

    key: str
    description: str
    amplitude: float = 0.0
    period: float | None = None
    phase: float = 0.5
    centre: float | None = None
    width: float | None = None

    def __post_init__(self) -> None:
        if self.period is not None and self.centre is not None:
            raise ValueError(
                f"scenario {self.key!r} declares both a fringe period and a line centre; a "
                f"scenario is one deviation, so that the tests can say where it is."
            )
        if self.centre is not None and self.width is None:
            raise ValueError(f"scenario {self.key!r} declares a line centre but no width.")

    @property
    def localised_at(self) -> float | None:
        """Where the deviation lives, micron — or ``None`` if it is everywhere.

        A fringe has no location: it is the same amplitude across the whole
        band, so "did the GP localise it?" is not a question with an answer.
        A line does, and
        :attr:`examples.m2_misspecification.study.Diagnosis.localisation_peak`
        is held to it.
        """
        return self.centre

    def deviation(self, wavelength: Any) -> np.ndarray:
        """The fractional deviation ``delta(lambda)``; the flux is ``truth * (1 + delta)``."""
        grid = np.asarray(wavelength, dtype=float)
        if self.amplitude == 0.0:
            return np.zeros_like(grid)
        if self.period is not None:
            return self.amplitude * np.sin(2.0 * np.pi * grid / self.period + self.phase)
        assert self.centre is not None and self.width is not None  # __post_init__
        return self.amplitude * np.exp(-0.5 * ((grid - self.centre) / self.width) ** 2)


#: The four scenarios, in the order every figure and table uses them.
SCENARIOS: tuple[Scenario, ...] = (
    Scenario("none", "No misspecification"),
    Scenario("mild", "Mild (smooth fringing)", amplitude=0.025, period=0.0028),
    Scenario("strong_smooth", "Strong (smooth fringing)", amplitude=0.07, period=0.0028),
    Scenario(
        "strong_sharp",
        "Strong (sharp line)",
        amplitude=0.12,
        centre=0.8630,
        width=0.00035,
    ),
)


def scenario_named(key: str) -> Scenario:
    """The :class:`Scenario` called *key*, or a ``KeyError`` naming the four."""
    for scenario in SCENARIOS:
        if scenario.key == key:
            return scenario
    raise KeyError(f"unknown scenario {key!r}; the four are {[s.key for s in SCENARIOS]}.")


def wavelength_grid(size: int = 200) -> np.ndarray:
    """*size* points spanning the same 0.842-0.872 µm band, micron.

    The same band at every size, never a longer one: the ladder samples one
    spectrum more finely, so the deviation, the physics and the noise level are
    held fixed along it.
    """
    if size < 2:
        raise ValueError(f"a spectrum needs at least two points, got {size!r}.")
    return np.linspace(WAVELENGTH_MIN, WAVELENGTH_MAX, int(size), dtype=float)


@dataclasses.dataclass(frozen=True)
class SyntheticSpectrum:
    """One generated spectrum, with every stage of its construction kept.

    Attributes
    ----------
    scenario
        Which deviation was applied.
    wavelength
        The grid, micron.
    truth
        The fitting model at :data:`~examples.m2_misspecification.model.TRUTH`
        — what a correctly specified experiment would have observed in the
        absence of noise.
    deviated
        ``truth * (1 + delta)``: the noise-free flux that was actually
        observed through, and the thing the fitting model cannot reproduce.
    observed
        ``deviated`` plus one seeded Gaussian noise realisation.
    uncertainty
        The per-point 1-sigma the fit is told about.
    seed
        The seed the noise came from.
    """

    scenario: Scenario
    wavelength: np.ndarray
    truth: np.ndarray
    deviated: np.ndarray
    observed: np.ndarray
    uncertainty: np.ndarray
    seed: int

    @property
    def size(self) -> int:
        """Number of points."""
        return int(self.wavelength.size)

    @property
    def deviation(self) -> np.ndarray:
        """The fractional deviation applied, on this grid."""
        return self.scenario.deviation(self.wavelength)

    def container(self) -> Spectrum:
        """The observed data as the :class:`~ampere.core.Spectrum` a dataset takes."""
        return Spectrum(
            self.wavelength * WAVELENGTH_UNIT,
            self.observed * FLUX_UNIT,
            uncertainty=self.uncertainty * FLUX_UNIT,
        )

    def truth_container(self) -> Spectrum:
        """The undeviated truth as a container, for overplotting on a figure."""
        return Spectrum(self.wavelength * WAVELENGTH_UNIT, self.truth * FLUX_UNIT)


def generate(
    scenario: Scenario | str,
    *,
    size: int = 200,
    seed: int = SEED,
    noise_fraction: float = NOISE_FRACTION,
) -> SyntheticSpectrum:
    """Generate the observed spectrum for *scenario* at *size* points.

    Deterministic in ``(scenario, size, seed)``: the noise comes from
    ``numpy.random.default_rng(seed)`` alone, so two calls with the same
    arguments return the same numbers, and the four scenarios at one size share
    a noise realisation (the generator is drawn after the deviation is applied,
    from a generator that has done nothing else). Sharing it is deliberate —
    the four rows of the paper's first figure then differ only by their
    deviation, which is what makes them comparable by eye.

    Examples
    --------
    >>> data = generate("strong_smooth", size=200)
    >>> data.size
    200
    >>> bool(np.allclose(data.deviated, data.truth * (1.0 + data.deviation)))
    True
    >>> float(np.max(np.abs(data.deviation))).__round__(3)
    0.07
    """
    chosen = scenario_named(scenario) if isinstance(scenario, str) else scenario
    grid = wavelength_grid(size)
    truth = flux_at(grid, **TRUTH)
    deviated = truth * (1.0 + chosen.deviation(grid))
    uncertainty = float(noise_fraction) * np.abs(truth)
    rng = np.random.default_rng(seed)
    observed = deviated + rng.normal(0.0, uncertainty)
    return SyntheticSpectrum(
        scenario=chosen,
        wavelength=grid,
        truth=truth,
        deviated=deviated,
        observed=observed,
        uncertainty=uncertainty,
        seed=int(seed),
    )
