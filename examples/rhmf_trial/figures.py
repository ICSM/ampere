"""Figures for the RHMF trial: one ``plot_anomaly_score`` panel per scenario, and the image panels.

Matplotlib's Agg backend, no display. Nothing here is written anywhere but the
run directory the caller names.
"""

from __future__ import annotations

import pathlib
from collections.abc import Mapping, Sequence
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from ampere.results import plot_anomaly_score

from .trial import (
    ImageCollection,
    RHMFFit,
    SpectraCollection,
    anomaly_score,
)


def save_spectra_figures(
    collection: SpectraCollection,
    fits: Mapping[tuple[int, float], RHMFFit],
    best: Mapping[str, Mapping[str, Any]],
    directory: pathlib.Path,
) -> list[pathlib.Path]:
    """One figure per deviated scenario, at that scenario's best grid point.

    Each draws the first deviated row's per-feature score through
    :func:`ampere.results.plot_anomaly_score` (so the provenance and the notes
    travel with it), shades the injected band, and overlays the first control
    row's score in grey. "Best" is the oracle choice of
    :func:`~examples.rhmf_trial.trial.best_by_scenario`.
    """
    written = []
    for key, record in best.items():
        fit = fits[(record["rank"], record["robust_scale"])]
        row = collection.rows(key)[0]
        control = collection.rows("none")[0]
        axes = plot_anomaly_score(anomaly_score(fit, collection.wavelength, row=row), color="C3")
        figure = axes.get_figure()
        axes.plot(
            collection.wavelength,
            fit.score()[control],
            color="0.45",
            linewidth=0.9,
            label="a control row",
        )
        band = collection.bands[key]
        axes.fill_between(
            collection.wavelength,
            0.0,
            1.0,
            where=band,
            color="C0",
            alpha=0.12,
            transform=axes.get_xaxis_transform(),
            label="injected deviation (|delta| >= half its peak)",
        )
        axes.set_xlabel("wavelength (micron)")
        axes.legend(loc="upper right", fontsize="x-small")
        figure.set_size_inches(9, 4.2)
        figure.suptitle(
            f"{key}: rank {record['rank']}, robust_scale {record['robust_scale']:g} "
            f"(excess contrast {record['excess']:.2f}, best of the grid with the answer known)",
            fontsize="small",
        )
        path = directory / f"spectra_{key}.png"
        figure.savefig(path, dpi=110, bbox_inches="tight")
        plt.close(figure)
        written.append(path)
    return written


def save_image_figure(
    collection: ImageCollection,
    collection_fit: RHMFFit,
    single_fits: Mapping[str, RHMFFit],
    directory: pathlib.Path,
    *,
    name: str = "image_scores.png",
) -> pathlib.Path:
    """The injected deviation and the score under each flattening, deviated and control."""
    pixels = collection.pixels
    deviated, control = collection.rows(True)[0], collection.rows(False)[0]
    weights = collection_fit.score().reshape(-1, pixels, pixels)
    panels: Sequence[tuple[str, np.ndarray]] = [
        ("injected deviation (noiseless)", collection.injected),
        ("collection: deviated image", weights[deviated]),
        ("collection: control image", weights[control]),
        ("single image: deviated", single_fits["deviated"].score()),
        ("single image: control", single_fits["control"].score()),
    ]
    figure, axes = plt.subplots(1, len(panels), figsize=(3.0 * len(panels), 3.2))
    for axis, (title, data) in zip(axes, panels, strict=True):
        shown = axis.imshow(data, origin="lower", cmap="viridis")
        axis.contour(collection.band, levels=[0.5], colors="w", linewidths=0.8)
        axis.set_title(title, fontsize="x-small")
        axis.set_xticks([])
        axis.set_yticks([])
        figure.colorbar(shown, ax=axis, fraction=0.046)
    figure.suptitle(
        f"image trial: collection fit rank {collection_fit.rank}, scale "
        f"{collection_fit.robust_scale:g}; white contour = half-peak of the injected background",
        fontsize="small",
    )
    path = directory / name
    figure.savefig(path, dpi=110, bbox_inches="tight")
    plt.close(figure)
    return path
