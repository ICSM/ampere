"""The six shipped plots, plus the calibration figures, for a stored run.

Following :mod:`examples.m2_misspecification.figures`'s split, but simpler:
this study has no side-by-side "paper" comparison panels, only
:func:`save_arm_figures` (:func:`~ampere.results.plot_corner`,
:func:`~ampere.results.plot_trace`,
:func:`~ampere.results.plot_posterior_predictive`, and — scoped by arm, as
``diagnostics.md`` §3.1 and §4 say — :func:`~ampere.results.plot_residuals`
on the standard-likelihood arms or
:func:`~ampere.results.plot_gp_localisation` and
:func:`~ampere.results.plot_anomaly_score` on the flexible one) and
:func:`save_calibration_figures` (:func:`~ampere.results.plot_sbc_ranks` and
:func:`~ampere.results.plot_coverage`, for :func:`.study.run_calibration`'s
output).

Nothing here is committed — every function takes the directory to write
into, and the caller decides (AGENTS.md ground rule 7).
"""

from __future__ import annotations

import pathlib
from collections.abc import Mapping
from typing import Any

from .generators import TRUTH

__all__ = ["save_arm_figures", "save_calibration_figures"]


def _pyplot() -> Any:
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    return plt


def save_arm_figures(
    results: Mapping[str, Mapping[str, Any]],
    directory: str | pathlib.Path,
    *,
    suffix: str = "png",
    thin: int = 20,
) -> list[pathlib.Path]:
    """Render the six-plot surface for every arm in a :func:`.study.run_study` result."""
    from ampere.results import (
        add_posterior_predictive,
        gp_localisation_score,
        plot_anomaly_score,
        plot_corner,
        plot_gp_localisation,
        plot_posterior_predictive,
        plot_residuals,
        plot_trace,
    )

    plt = _pyplot()
    target = pathlib.Path(directory)
    target.mkdir(parents=True, exist_ok=True)
    written: list[pathlib.Path] = []

    def _save(figure: Any, stem: str) -> None:
        path = target / f"{stem}.{suffix}"
        figure.savefig(path, bbox_inches="tight", dpi=120)
        plt.close(figure)
        written.append(path)

    names = list(TRUTH)
    for arm, entry in results.items():
        run = entry["run"]
        _save(plot_corner(run, var_names=names, truths=list(TRUTH.values())), f"corner_{arm}")
        _save(plot_trace(run, var_names=names), f"trace_{arm}")
        entry["run"] = run = add_posterior_predictive(run, entry["problem"], thin=thin)
        _save(plot_posterior_predictive(run, datasets=["vis", "t3"]), f"posterior_predictive_{arm}")
        if arm == "flexible":
            _save(plot_gp_localisation(run, datasets=["vis"]), f"localisation_{arm}")
            score = gp_localisation_score(run, dataset="vis")
            _save(plot_anomaly_score(score).get_figure(), f"anomaly_{arm}")
        else:
            _save(plot_residuals(run, datasets=["vis"]), f"residuals_{arm}")
    return written


def save_calibration_figures(
    calibrations: Mapping[str, Any],
    directory: str | pathlib.Path,
    *,
    suffix: str = "png",
) -> list[pathlib.Path]:
    """Render :func:`~ampere.results.plot_sbc_ranks` and
    :func:`~ampere.results.plot_coverage` for each of :func:`.study.run_calibration`'s outputs.
    """
    from ampere.results import plot_coverage, plot_sbc_ranks

    plt = _pyplot()
    target = pathlib.Path(directory)
    target.mkdir(parents=True, exist_ok=True)
    written: list[pathlib.Path] = []

    def _save(figure: Any, stem: str) -> None:
        path = target / f"{stem}.{suffix}"
        figure.savefig(path, bbox_inches="tight", dpi=120)
        plt.close(figure)
        written.append(path)

    for arm, calibration in calibrations.items():
        _save(plot_sbc_ranks(calibration).get_figure(), f"sbc_ranks_{arm}")
        _save(plot_coverage(calibration).get_figure(), f"coverage_{arm}")
    return written
