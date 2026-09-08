"""The data are what the study says they are.

Cheap, exact assertions about the experiment's inputs. They matter more than
their cost suggests: every claim in ``test_science.py`` is a claim about what a
likelihood did to a *particular* deviation, so a generator that silently
changed its amplitude, its location or its noise would turn a science result
into a coincidence.
"""

from __future__ import annotations

import astropy.units as u
import numpy as np
import pytest

from examples.m2_misspecification import generators
from examples.m2_misspecification.model import TRUTH, flux_at


def test_the_four_scenarios_are_the_paper_studys() -> None:
    keys = [scenario.key for scenario in generators.SCENARIOS]
    assert keys == ["none", "mild", "strong_smooth", "strong_sharp"]
    by_key = {scenario.key: scenario for scenario in generators.SCENARIOS}
    assert by_key["none"].amplitude == 0.0
    assert by_key["mild"].amplitude == pytest.approx(0.025)
    assert by_key["mild"].period == pytest.approx(0.0028)
    assert by_key["strong_smooth"].amplitude == pytest.approx(0.07)
    assert by_key["strong_smooth"].period == pytest.approx(0.0028)
    assert by_key["strong_sharp"].amplitude == pytest.approx(0.12)
    assert by_key["strong_sharp"].centre == pytest.approx(0.8630)
    assert by_key["strong_sharp"].width == pytest.approx(0.00035)


def test_only_the_line_scenario_is_localised() -> None:
    """A fringe is everywhere, so it has no location to be found at."""
    by_key = {scenario.key: scenario for scenario in generators.SCENARIOS}
    assert by_key["none"].localised_at is None
    assert by_key["mild"].localised_at is None
    assert by_key["strong_smooth"].localised_at is None
    assert by_key["strong_sharp"].localised_at == pytest.approx(0.8630)


def test_a_scenario_is_one_deviation() -> None:
    with pytest.raises(ValueError, match="both a fringe period and a line centre"):
        generators.Scenario("bad", "two deviations", 0.1, period=0.003, centre=0.86, width=1e-4)
    with pytest.raises(ValueError, match="no width"):
        generators.Scenario("bad", "no width", 0.1, centre=0.86)


@pytest.mark.parametrize("size", [200, 2_000, 20_000])
def test_the_grid_is_the_same_band_at_every_size(size: int) -> None:
    grid = generators.wavelength_grid(size)
    assert grid.size == size
    assert grid[0] == pytest.approx(generators.WAVELENGTH_MIN)
    assert grid[-1] == pytest.approx(generators.WAVELENGTH_MAX)
    assert np.all(np.diff(grid) > 0.0)


def test_generation_is_deterministic_and_decomposed() -> None:
    first = generators.generate("strong_smooth", size=200)
    second = generators.generate("strong_smooth", size=200)
    assert np.array_equal(first.observed, second.observed)
    # The three stages really are the three stages.
    assert np.allclose(first.truth, flux_at(first.wavelength, **TRUTH))
    assert np.allclose(first.deviated, first.truth * (1.0 + first.deviation))
    assert np.allclose(first.uncertainty, generators.NOISE_FRACTION * np.abs(first.truth))


def test_the_control_deviates_by_nothing() -> None:
    data = generators.generate("none", size=200)
    assert np.array_equal(data.deviated, data.truth)
    assert float(np.max(np.abs(data.deviation))) == 0.0


def test_the_four_scenarios_share_a_noise_realisation() -> None:
    """So the rows of figure 1 differ by their deviation and by nothing else."""
    residuals = {
        scenario.key: generators.generate(scenario, size=200).observed
        - generators.generate(scenario, size=200).deviated
        for scenario in generators.SCENARIOS
    }
    control = residuals["none"]
    for key, values in residuals.items():
        assert np.allclose(values, control), key


@pytest.mark.parametrize("key", ["none", "mild", "strong_smooth", "strong_sharp"])
def test_the_deviation_has_the_amplitude_it_claims(key: str) -> None:
    scenario = generators.scenario_named(key)
    # 20 000 points resolves both the 0.0028 um fringe and the 0.00035 um line.
    data = generators.generate(scenario, size=20_000)
    assert float(np.max(np.abs(data.deviation))) == pytest.approx(scenario.amplitude, rel=1e-3)


def test_the_container_carries_the_right_units() -> None:
    data = generators.generate("mild", size=200)
    container = data.container()
    assert container.spectral_axis.unit == u.micron
    assert container.unit == u.Jy
    assert container.uncertainty is not None
    assert np.allclose(container.values, data.observed)


def test_an_unknown_scenario_names_the_four() -> None:
    with pytest.raises(KeyError, match="strong_sharp"):
        generators.scenario_named("nope")
