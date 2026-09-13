"""The plots this study can actually render, plus the calibration figures.

Found running this item, not assumed: four of the "six shipped plots"
(:func:`~ampere.results.plot_posterior_predictive`,
:func:`~ampere.results.plot_residuals`,
:func:`~ampere.results.plot_gp_localisation`,
:func:`~ampere.results.plot_anomaly_score`) read a stored group's single
coordinate axis (:func:`ampere.results._plotting.coordinate_of`), and refuse
by name for a "point kind with several axes" — which both
:class:`~ampere.core.VisibilitySet` and
:class:`~ampere.core.ClosurePhases` are, deliberately: ``results.md`` §4 and
``DEVELOPMENT_PLAN.md`` §4.4/§4.8 stage "gridded and multi-axis kinds" for
Phase 5, and the refusal is the same one whether the derived variable is
complex (W4.2's own carried note — "W4.4's figures must pick a component or
modulus" — undersold the gap: even a real, single-component view still has
no *ordered* coordinate to plot against). So this study's two datasets
cannot use any of those four today, on either arm.

What *does* render, and what this module renders: :func:`~ampere.results.plot_corner`
and :func:`~ampere.results.plot_trace` (posterior-only, no data coordinate
needed) for every arm, and — for :func:`.study.run_calibration`'s output,
which lives in *parameter* space rather than on a data axis and so does not
hit the same wall — :func:`~ampere.results.plot_sbc_ranks` and
:func:`~ampere.results.plot_coverage`.

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
) -> list[pathlib.Path]:
    """Render what this study's datasets support: corner and trace, per arm.

    See the module docstring for why ``plot_posterior_predictive``,
    ``plot_residuals``, ``plot_gp_localisation`` and ``plot_anomaly_score``
    are not called here: all four refuse a multi-axis point kind by name, and
    both of this study's datasets are one.
    """
    from ampere.results import plot_corner, plot_trace

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
