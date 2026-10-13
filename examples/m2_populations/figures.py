"""The population scenario's figure (**W7.2**).

One figure, two panels when a comparison is given and one when it is not:

* **the collection** — each spectrum's departure from the undeviated truth,
  ``observed - truth``, offset vertically, the deviated members and the
  poorly constrained one marked. This is the panel the documentation page
  shows, drawn from :func:`~examples.m2_populations.populations.generate_population`
  directly: it is what the three declarations are fitted *to*, and it costs
  no sampling, so a docs build can afford it.
* **the amplitudes** — each member's GP amplitude posterior (median and 68 %
  interval) under each arm, on a log axis, from a
  :class:`~examples.m2_populations.populations.Comparison`. Written by
  ``python -m examples.m2_populations --figures DIR``.
"""

from __future__ import annotations

import pathlib
from collections.abc import Sequence
from typing import Any

import numpy as np

from examples.m2_misspecification.generators import SyntheticSpectrum, scenario_named

from .populations import (
    ARMS,
    DEVIATION,
    POORLY_CONSTRAINED,
    Comparison,
    dataset_label,
)

__all__ = ["FIGURE_NAME", "figure_population", "save_figure"]

#: The file stem :func:`save_figure` writes.
FIGURE_NAME = "m2_populations"

_ARM_MARKERS = {"independent": "o", "tied": "s", "population": "D"}


def _pyplot() -> Any:
    """``matplotlib.pyplot``, imported late, for the study's reason."""
    import matplotlib.pyplot as plt

    return plt


def figure_population(
    data: Sequence[SyntheticSpectrum], comparison: Comparison | None = None
) -> Any:
    """The collection and, given *comparison*, the members' amplitudes per arm."""
    plt = _pyplot()
    panels = 1 if comparison is None else 2
    figure, axes = plt.subplots(1, panels, figsize=(5.2 * panels, 5.0), squeeze=False)
    left = axes[0, 0]
    step = 0.06
    centre = scenario_named(DEVIATION).localised_at
    for index, spectrum in enumerate(data):
        offset = -step * index
        deviated = spectrum.scenario.amplitude != 0.0
        colour = "C3" if deviated else ("C1" if index == POORLY_CONSTRAINED else "C0")
        left.plot(
            spectrum.wavelength,
            spectrum.observed - spectrum.truth + offset,
            color=colour,
            lw=0.8,
        )
        note = spectrum.scenario.key if deviated else "none"
        if index == POORLY_CONSTRAINED and not deviated:
            note = "none, 5x noise"
        left.text(
            spectrum.wavelength[-1] + 0.0005,
            offset,
            f"{dataset_label(index)} ({note})",
            va="center",
            fontsize=8,
            color=colour,
        )
    if centre is not None:
        left.axvline(centre, color="0.6", ls=":", lw=0.8)
    left.set_xlabel("wavelength (µm)")
    left.set_ylabel("observed - truth (Jy), offset")
    left.set_title("The collection: one object class, eight spectra")
    left.set_xlim(spectrum.wavelength[0], spectrum.wavelength[-1] + 0.009)
    if comparison is not None:
        right = axes[0, 1]
        arms = [arm for arm in ARMS if arm in comparison.fits]
        positions = np.arange(len(comparison.data))
        for shift, arm in enumerate(arms):
            summaries = [comparison.amplitude(arm, i) for i in positions]
            medians = np.array([s.median for s in summaries])
            lower = medians - np.array([s.low for s in summaries])
            upper = np.array([s.high for s in summaries]) - medians
            right.errorbar(
                positions + 0.2 * (shift - 1),
                medians,
                yerr=[lower, upper],
                fmt=_ARM_MARKERS.get(arm, "o"),
                ms=4,
                capsize=2,
                label=arm,
            )
        right.set_yscale("log")
        right.set_xticks(positions, [dataset_label(i) for i in positions])
        right.set_ylabel("GP amplitude (Jy)")
        right.set_title(f"Amplitudes per arm ({comparison.backend})")
        right.legend(fontsize=8)
    figure.tight_layout()
    return figure


def save_figure(
    data: Sequence[SyntheticSpectrum],
    directory: str | pathlib.Path,
    comparison: Comparison | None = None,
    *,
    suffix: str = "pdf",
) -> list[pathlib.Path]:
    """Write :func:`figure_population` into *directory* and return its path."""
    plt = _pyplot()
    target = pathlib.Path(directory)
    target.mkdir(parents=True, exist_ok=True)
    figure = figure_population(data, comparison)
    path = target / f"{FIGURE_NAME}.{suffix}"
    figure.savefig(path, bbox_inches="tight", dpi=150)
    plt.close(figure)
    return [path]
