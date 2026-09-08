"""The figures: the paper's three, plus the ``ampere.results`` renderers.

Two kinds of figure, and the distinction is worth being explicit about.

**The paper's three** — spectra with residuals, parameter recovery, the 1-D
posteriors — are *comparisons*: each panel holds two fits of the same data side
by side, and no single-run renderer can draw that, because a run does not know
what it is being compared with. They are composed here, but every number in
them comes out of the stored runs and the groups ``ampere.results`` derived
(``residuals``, ``gp_localisation``), never from a re-fit or a re-computation.
:func:`~examples.m2_misspecification.study.prepare` is their precondition.

**The shipped renderers** — :func:`ampere.results.plot_corner`,
:func:`~ampere.results.plot_trace`,
:func:`~ampere.results.plot_posterior_predictive`,
:func:`~ampere.results.plot_residuals`,
:func:`~ampere.results.plot_gp_localisation` and
:func:`~ampere.results.plot_anomaly_score` — are what a *user* of this study
would reach for, one run at a time, and :func:`save_result_figures` renders
each of them so that the study exercises W2.8's surface rather than only
describing it. Two of them are scoped: ``plot_residuals`` warns on a GP fit
(its whiteness family is for standard-likelihood fits) and
``plot_gp_localisation`` needs a GP, so each is rendered for the runs it is
for.

Nothing here is committed. Every function takes the directory to write into,
and the caller decides — ``tests/m2/test_figures.py`` gives it a ``tmp_path``,
``python -m examples.m2_misspecification`` gives it whatever ``--figures``
named. AGENTS.md ground rule 7: no binary artefacts in git.
"""

from __future__ import annotations

import pathlib
from collections.abc import Mapping, MutableMapping
from typing import Any

import numpy as np

from ampere.results import (
    GP_LOCALISATION_GROUP,
    plot_anomaly_score,
    plot_corner,
    plot_gp_localisation,
    plot_residuals,
    plot_trace,
)

from .generators import SCENARIOS
from .model import PARAMETER_NAMES, TRUTH, flux_at
from .study import PHYSICAL_NAMES, summarise

__all__ = [
    "COLOUR_FLEXIBLE",
    "COLOUR_STANDARD",
    "PAPER_FIGURES",
    "figure_parameter_recovery",
    "figure_posteriors_1d",
    "figure_spectra_residuals",
    "save_paper_figures",
    "save_result_figures",
]

#: The paper study's two colours, kept so the reproduction is recognisable
#: beside the original: blue for the standard likelihood, red for the flexible.
COLOUR_STANDARD = "#2166ac"
COLOUR_FLEXIBLE = "#ff1d38"

#: The three figures :func:`save_paper_figures` writes, in order.
PAPER_FIGURES: tuple[str, ...] = (
    "fig_spectra_residuals",
    "fig_parameter_recovery",
    "fig_posteriors_1d",
)

_LATEX = {"A": r"$A$", "B": r"$B$", "d1": r"$d_1$", "d2": r"$d_2$"}
_SHORT = {
    "none": "None",
    "mild": "Mild\n(smooth)",
    "strong_smooth": "Strong\n(smooth)",
    "strong_sharp": "Strong\n(sharp)",
}


def _pyplot() -> Any:
    """``matplotlib.pyplot``, imported here rather than at module scope.

    Importing pyplot chooses a backend, and a module that does so on import
    decides for its importer. ``tests/m2`` selects ``Agg`` before it imports
    anything; a user running the study in a notebook keeps whatever they had.
    """
    import matplotlib.pyplot as plt

    return plt


def _draws(run: Any, name: str) -> np.ndarray:
    return np.asarray(run["posterior"][name].values, dtype=float).ravel()


def _localisation(entry: Mapping[str, Any]) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """The conditioned GP mean's median and 68 % band, if the group is there."""
    run = entry["run"]
    if GP_LOCALISATION_GROUP not in run.children:
        return None
    group = run[GP_LOCALISATION_GROUP]
    label = f"{next(iter(entry['problem'].datasets))}_mean"
    if label not in group.dataset.data_vars:
        return None
    values = np.asarray(group[label].values, dtype=float)
    flat = values.reshape(-1, values.shape[-1])
    return (
        np.nanmedian(flat, axis=0),
        np.nanpercentile(flat, 16.0, axis=0),
        np.nanpercentile(flat, 84.0, axis=0),
    )


def figure_spectra_residuals(results: Mapping[tuple[str, str], Mapping[str, Any]]) -> Any:
    """Figure 1: the data and both fits, and their residuals with the GP band.

    Four rows, one per scenario; the left column holds the observed spectrum,
    the truth and the two posterior-median fits, the right column the residuals
    about each fit. On the flexible row the conditioned GP mean is drawn over
    its residuals with the 68 % band across posterior draws — that curve is
    what the flexible likelihood *says* the model is missing, and seeing it
    trace the injected fringe (or spike at the injected line) is the figure's
    whole point.
    """
    plt = _pyplot()
    figure, axes = plt.subplots(
        len(SCENARIOS),
        2,
        figsize=(9.0, 10.5),
        sharex="col",
        gridspec_kw={"width_ratios": [1.25, 1.0], "hspace": 0.10, "wspace": 0.28},
    )
    for row, scenario in enumerate(SCENARIOS):
        standard = results[(scenario.key, "standard")]
        flexible = results[(scenario.key, "flexible")]
        data = standard["data"]
        grid = data.wavelength
        fits = {}
        for kind, entry in (("standard", standard), ("flexible", flexible)):
            medians = {
                name: float(np.median(_draws(entry["run"], f"model.{name}")))
                for name in PARAMETER_NAMES
            }
            fits[kind] = flux_at(grid, **medians)

        left = axes[row, 0]
        left.errorbar(
            grid,
            data.observed,
            yerr=data.uncertainty,
            fmt=".",
            color="0.6",
            ms=2.0,
            lw=0.7,
            alpha=0.8,
            zorder=1,
            rasterized=True,
            label="observed" if row == 0 else None,
        )
        left.plot(grid, data.truth, color="0.35", lw=1.6, ls="--", zorder=2, label="truth")
        left.plot(grid, fits["standard"], color=COLOUR_STANDARD, lw=1.8, zorder=3, label="standard")
        left.plot(grid, fits["flexible"], color=COLOUR_FLEXIBLE, lw=1.0, zorder=4, label="flexible")
        left.set_ylabel("flux (Jy)")
        left.text(
            0.02,
            0.94,
            f"({'abcd'[row]}) {scenario.description}",
            transform=left.transAxes,
            fontsize=8,
            fontweight="bold",
            va="top",
            bbox={"boxstyle": "round,pad=0.3", "fc": "white", "ec": "0.8", "alpha": 0.9},
        )
        if row == 0:
            left.legend(loc="upper right", fontsize=7, framealpha=0.9)
        if row == len(SCENARIOS) - 1:
            left.set_xlabel(r"wavelength ($\mu$m)")

        right = axes[row, 1]
        right.axhline(0.0, color="0.7", lw=0.6)
        right.fill_between(
            grid,
            -data.uncertainty,
            data.uncertainty,
            color="0.7",
            alpha=0.5,
            zorder=0,
            label=r"$\pm 1\sigma$ noise" if row == 0 else None,
        )
        for kind, colour in (("standard", COLOUR_STANDARD), ("flexible", COLOUR_FLEXIBLE)):
            right.plot(
                grid,
                data.observed - fits[kind],
                ".",
                color=colour,
                ms=2.2,
                alpha=0.75,
                zorder=2,
                rasterized=True,
                label=kind if row == 0 else None,
            )
        band = _localisation(flexible)
        if band is not None:
            centre, lower, upper = band
            right.fill_between(grid, lower, upper, color=COLOUR_FLEXIBLE, alpha=0.18, zorder=1)
            right.plot(
                grid,
                centre,
                color=COLOUR_FLEXIBLE,
                lw=1.1,
                zorder=3,
                label="GP mean" if row == 0 else None,
            )
        right.set_ylabel("residual (Jy)")
        if row == 0:
            right.legend(loc="upper right", fontsize=6.5, framealpha=0.9)
        if row == len(SCENARIOS) - 1:
            right.set_xlabel(r"wavelength ($\mu$m)")
    return figure


def figure_parameter_recovery(results: Mapping[tuple[str, str], Mapping[str, Any]]) -> Any:
    """Figure 2: median and 68 % interval per parameter, per scenario, per likelihood.

    The milestone's headline picture. The dashed line is the truth, and the
    claim to read off the figure is not that the red points sit on it — they do
    not always — but that the red *error bars* cross it in every scenario while
    the blue ones stop crossing it as soon as the model is wrong.
    """
    plt = _pyplot()
    figure, axes = plt.subplots(1, len(PARAMETER_NAMES), figsize=(10.0, 3.0))
    positions = np.arange(len(SCENARIOS))
    offset = 0.14
    for column, parameter in enumerate(PARAMETER_NAMES):
        panel = axes[column]
        name = f"model.{parameter}"
        for index, scenario in enumerate(SCENARIOS):
            for kind, colour, shift, marker in (
                ("standard", COLOUR_STANDARD, -offset, "o"),
                ("flexible", COLOUR_FLEXIBLE, offset, "s"),
            ):
                summary = summarise(results[(scenario.key, kind)]["run"], names=[name])[name]
                panel.errorbar(
                    index + shift,
                    summary.median,
                    yerr=[[summary.median - summary.low], [summary.high - summary.median]],
                    fmt=marker,
                    color=colour,
                    ms=4.5,
                    capsize=2.5,
                    lw=1.2,
                    zorder=3,
                    label=kind if (index == 0 and column == 0) else None,
                )
        panel.axhline(TRUTH[parameter], color="0.35", ls="--", lw=0.9, zorder=1)
        panel.set_xticks(positions)
        panel.set_xticklabels([_SHORT[s.key] for s in SCENARIOS], fontsize=6.5)
        panel.set_title(_LATEX[parameter], fontsize=11)
        panel.yaxis.grid(True, alpha=0.3, lw=0.5)
        panel.set_axisbelow(True)
        if column == 0:
            panel.legend(loc="best", fontsize=7, framealpha=0.9)
    figure.tight_layout()
    return figure


def figure_posteriors_1d(results: Mapping[tuple[str, str], Mapping[str, Any]]) -> Any:
    """Figure 3: the marginal posteriors themselves, standard against flexible.

    What figure 2 summarises into a median and an interval, this one shows
    whole: the standard likelihood's narrow spike well away from the truth, and
    the flexible likelihood's broader distribution containing it.
    """
    plt = _pyplot()
    figure, axes = plt.subplots(len(SCENARIOS), len(PARAMETER_NAMES), figsize=(10.0, 8.0))
    for row, scenario in enumerate(SCENARIOS):
        for column, parameter in enumerate(PARAMETER_NAMES):
            panel = axes[row, column]
            name = f"model.{parameter}"
            standard = _draws(results[(scenario.key, "standard")]["run"], name)
            flexible = _draws(results[(scenario.key, "flexible")]["run"], name)
            both = np.concatenate([standard, flexible])
            bins = np.linspace(*np.percentile(both, (0.5, 99.5)), 45)
            panel.hist(
                standard,
                bins=bins,
                density=True,
                alpha=0.45,
                color=COLOUR_STANDARD,
                label="standard" if (row == 0 and column == 0) else None,
            )
            panel.hist(
                flexible,
                bins=bins,
                density=True,
                alpha=0.45,
                color=COLOUR_FLEXIBLE,
                label="flexible" if (row == 0 and column == 0) else None,
            )
            panel.axvline(TRUTH[parameter], color="0.3", ls="--", lw=0.9)
            panel.set_yticks([])
            panel.tick_params(labelsize=6)
            if row == 0:
                panel.set_title(_LATEX[parameter], fontsize=10)
            if column == 0:
                panel.set_ylabel(scenario.description, fontsize=6.5, fontweight="bold")
            if row == 0 and column == 0:
                panel.legend(fontsize=6, loc="upper left")
    figure.tight_layout()
    return figure


def save_paper_figures(
    results: Mapping[tuple[str, str], Mapping[str, Any]],
    directory: str | pathlib.Path,
    *,
    suffix: str = "pdf",
) -> list[pathlib.Path]:
    """Write the three comparison figures into *directory* and return their paths."""
    plt = _pyplot()
    target = pathlib.Path(directory)
    target.mkdir(parents=True, exist_ok=True)
    builders = (
        figure_spectra_residuals,
        figure_parameter_recovery,
        figure_posteriors_1d,
    )
    written: list[pathlib.Path] = []
    for name, builder in zip(PAPER_FIGURES, builders, strict=True):
        figure = builder(results)
        path = target / f"{name}.{suffix}"
        figure.savefig(path, bbox_inches="tight", dpi=150)
        plt.close(figure)
        written.append(path)
    return written


def save_result_figures(
    results: MutableMapping[tuple[str, str], MutableMapping[str, Any]],
    directory: str | pathlib.Path,
    *,
    suffix: str = "png",
    scenarios: tuple[str, ...] | None = None,
) -> list[pathlib.Path]:
    """Render W2.8's shipped plotting surface, one figure per run.

    Each renderer is called on the runs it is *for*: the corner and trace plots
    on every run, ``plot_residuals`` on the standard fits (its whiteness family
    is scoped to them), ``plot_gp_localisation`` and ``plot_anomaly_score`` on
    the flexible ones. That scoping is the point rather than a convenience —
    ``diagnostics.md`` §3.1 and §4 say which question each answers, and running
    the wrong one produces a figure that means something else.
    """
    from ampere.results import gp_localisation_score

    plt = _pyplot()
    target = pathlib.Path(directory)
    target.mkdir(parents=True, exist_ok=True)
    written: list[pathlib.Path] = []

    def _save(figure: Any, stem: str) -> None:
        path = target / f"{stem}.{suffix}"
        figure.savefig(path, bbox_inches="tight", dpi=120)
        plt.close(figure)
        written.append(path)

    for (key, kind), entry in results.items():
        if scenarios is not None and key not in scenarios:
            continue
        run = entry["run"]
        label = next(iter(entry["problem"].datasets))
        _save(
            plot_corner(run, var_names=list(PHYSICAL_NAMES), truths=list(TRUTH.values())),
            f"corner_{key}_{kind}",
        )
        _save(plot_trace(run, var_names=list(PHYSICAL_NAMES)), f"trace_{key}_{kind}")
        if kind == "standard":
            _save(plot_residuals(run, datasets=[label]), f"residuals_{key}_{kind}")
        else:
            _save(plot_gp_localisation(run, datasets=[label]), f"localisation_{key}_{kind}")
            score = gp_localisation_score(run, dataset=label)
            figure = plot_anomaly_score(score).get_figure()
            _save(figure, f"anomaly_{key}_{kind}")
    return written
