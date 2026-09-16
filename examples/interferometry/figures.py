"""The plots this study renders, plus the calibration figures.

**W5.3 lifts the limitation found running this item.** Four of the "six
shipped plots" (:func:`~ampere.results.plot_posterior_predictive`,
:func:`~ampere.results.plot_residuals`, :func:`~ampere.results.plot_gp_localisation`,
:func:`~ampere.results.plot_anomaly_score`) used to refuse a point kind with
several axes by name — both :class:`~ampere.core.VisibilitySet` and
:class:`~ampere.core.ClosurePhases` are that, deliberately (``results.md``
§4 and ``DEVELOPMENT_PLAN.md`` §4.4/§4.8 staged "gridded and multi-axis
kinds" for Phase 5) — so this study's two datasets could not use any of
those four. W5.3 gives ``ampere.results._plotting.coordinate_of`` a
``coordinate=`` argument and lets a kind declare a default
(``VisibilitySet``'s is baseline length; ``ClosurePhases``'s is its longest
baseline), and gives the complex half its own answer, ``component=``. This
module now renders all six shipped plots, per arm: ``vis`` (complex) with
``component="abs"`` and its default coordinate, ``t3`` (real) with its
default coordinate and no ``component=``, and — only on the ``"flexible"``
arm, the one dataset any of this study's models fits with a
:class:`~ampere.core.GaussianProcessNoise` — ``plot_gp_localisation`` and
``plot_anomaly_score`` for ``vis`` as well.
:doc:`/interferometry`'s "The plots: a found limitation, lifted at W5.3"
section is the worked account.

Nothing here is committed — every function takes the directory to write
into, and the caller decides (AGENTS.md ground rule 7).
"""

from __future__ import annotations

import pathlib
import warnings
from collections.abc import Mapping
from typing import Any

from .generators import TRUTH

__all__ = ["save_arm_figures", "save_calibration_figures"]

#: The complex visibility dataset is always plotted as its amplitude — the
#: same view :meth:`~ampere.core.VisibilitySet.amplitude` exposes on the
#: container itself, and a natural default for a study whose flexible arm's
#: GP is over the amplitude-dominated correlated structure the disc leaves.
_VISIBILITY_COMPONENT = "abs"


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
    """Render all six shipped plots this study's data now supports, per arm.

    ``plot_corner`` and ``plot_trace`` need nothing beyond the run itself.
    The four data-coordinate diagnostics (**W5.3**) are computed on demand
    per ``results.md`` §7 — ``add_posterior_predictive``/``add_residuals``
    first, ``gp_localisation`` only where a dataset actually carries a GP —
    and each entry's ``"run"`` is updated in place with the groups this adds,
    so a caller inspecting ``results`` afterwards sees them too.
    """
    from ampere.results import (
        add_posterior_predictive,
        add_residuals,
        gp_localisation,
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
        problem = entry["problem"]
        _save(plot_corner(run, var_names=names, truths=list(TRUTH.values())), f"corner_{arm}")
        _save(plot_trace(run, var_names=names), f"trace_{arm}")

        # Posterior-predictive replicates: vis as its amplitude, t3 as itself
        # (real-valued; closure phases have no component to pick).
        run = add_posterior_predictive(
            run, problem, datasets=["vis"], component=_VISIBILITY_COMPONENT
        )
        _save(
            plot_posterior_predictive(run, datasets=["vis"], component=_VISIBILITY_COMPONENT),
            f"posterior_predictive_vis_{arm}",
        )
        run = add_posterior_predictive(run, problem, datasets=["t3"])
        _save(
            plot_posterior_predictive(run, datasets=["t3"]),
            f"posterior_predictive_t3_{arm}",
        )

        # Signed residuals. The "flexible" arm's vis is GP-fitted, so
        # plot_residuals warns (diagnostics.md §3.1: whiteness is close to
        # circular for a GP-augmented fit) rather than refusing — expected
        # here, since that is exactly the arm plot_gp_localisation covers
        # below, and not a defect this study's figures need to surface twice.
        run = add_residuals(run, problem, datasets=["vis"], component=_VISIBILITY_COMPONENT)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            residuals_vis = plot_residuals(run, datasets=["vis"], component=_VISIBILITY_COMPONENT)
        _save(residuals_vis, f"residuals_vis_{arm}")
        run = add_residuals(run, problem, datasets=["t3"])
        _save(plot_residuals(run, datasets=["t3"]), f"residuals_t3_{arm}")

        if arm == "flexible":
            run = gp_localisation(run, problem, datasets=["vis"], component=_VISIBILITY_COMPONENT)
            _save(
                plot_gp_localisation(run, datasets=["vis"], component=_VISIBILITY_COMPONENT),
                f"gp_localisation_vis_{arm}",
            )
            score = gp_localisation_score(run, dataset="vis", component=_VISIBILITY_COMPONENT)
            axes = plot_anomaly_score(score)
            _save(axes.get_figure(), f"anomaly_score_vis_{arm}")

        entry["run"] = run
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
