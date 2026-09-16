"""The plots this study can render, and the ones it has to draw by hand.

The same wall :mod:`examples.interferometry.figures` hit, reached from the
other side. Four of the six shipped plots —
:func:`~ampere.results.plot_posterior_predictive`,
:func:`~ampere.results.plot_residuals`,
:func:`~ampere.results.plot_gp_localisation` and
:func:`~ampere.results.plot_anomaly_score` — read a stored group through
``ampere.results._plotting.coordinate_of``, which needs **exactly one ordered
coordinate axis**. An :class:`~ampere.core.Image` has two, so all four refuse
it, on every arm.

That refusal has a different flavour here, and it is worth saying which.
Interferometry's kinds are multi-axis *point* sets with no natural order at
all, so there is genuinely nothing to plot a residual against — a
one-dimensional panel of a visibility set would be a line through an
unordered cloud. An image's two axes *are* ordered, and the natural picture is
not a line at all: it is a **panel**. So the gap for a gridded kind is not "no
coordinate to use" but "the renderer's output shape is wrong" — a
one-dimensional figure where a two-dimensional one is wanted. Those are
different gaps behind one refusal message, and a future item that widens
``coordinate_of`` to two axes will still have to decide what a residual
*picture* is before any of the four means anything here.

So this module renders :func:`~ampere.results.plot_corner` and
:func:`~ampere.results.plot_trace` (posterior-only, no data coordinate needed),
the two calibration figures (parameter space, same reason), and — because the
interesting thing about an image fit is the residual map — a hand-rolled
``matplotlib`` panel row of data, model and residual per arm. The panels are
this module's own code, deliberately: they are what the shipped renderer would
have to produce, drawn once here so that the gap is demonstrated rather than
described.

Nothing is committed — every function takes the directory to write into, and
the caller decides (AGENTS.md ground rule 7).
"""

from __future__ import annotations

import pathlib
from collections.abc import Mapping
from typing import Any

import numpy as np

from .generators import TRUTH

__all__ = ["save_arm_figures", "save_calibration_figures", "save_image_panels"]


def _pyplot() -> Any:
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    return plt


def _writer(directory: str | pathlib.Path, suffix: str) -> tuple[Any, list[pathlib.Path], Any]:
    plt = _pyplot()
    target = pathlib.Path(directory)
    target.mkdir(parents=True, exist_ok=True)
    written: list[pathlib.Path] = []

    def save(figure: Any, stem: str) -> None:
        path = target / f"{stem}.{suffix}"
        figure.savefig(path, bbox_inches="tight", dpi=120)
        plt.close(figure)
        written.append(path)

    return plt, written, save


def save_arm_figures(
    results: Mapping[str, Mapping[str, Any]],
    directory: str | pathlib.Path,
    *,
    suffix: str = "png",
) -> list[pathlib.Path]:
    """Corner and trace per arm: the two shipped plots a gridded kind can use."""
    from ampere.results import plot_corner, plot_trace

    _, written, save = _writer(directory, suffix)
    names = list(TRUTH)
    for arm, entry in results.items():
        run = entry["run"]
        save(plot_corner(run, var_names=names, truths=list(TRUTH.values())), f"corner_{arm}")
        save(plot_trace(run, var_names=names), f"trace_{arm}")
    return written


def save_image_panels(
    results: Mapping[str, Mapping[str, Any]],
    directory: str | pathlib.Path,
    *,
    suffix: str = "png",
) -> list[pathlib.Path]:
    """Data, posterior-median model and residual, as three panels per arm.

    Hand-rolled, for the reason the module docstring gives. The model panel is
    the prediction at the posterior **median** rather than a band: a band over
    a two-dimensional field is a third dimension, and choosing how to show it
    is exactly the design question a future gridded ``plot_residuals`` owes.
    """
    _, written, save = _writer(directory, suffix)
    plt = _pyplot()
    for arm, entry in results.items():
        problem = entry["problem"]
        run = entry["run"]
        observed = problem.datasets["image"].observed
        posterior = run["posterior"]
        free = set(problem.parameters.free_names)
        median = {
            str(name): float(np.median(np.asarray(posterior[name].values, dtype=float)))
            for name in posterior.dataset.data_vars
            if str(name) in free
        }
        # ``simulate(observe=False)`` is the public route to a prediction: it
        # runs the model and the instrument and stops before the noise draw,
        # which is exactly the noiseless map this panel wants.
        predicted = problem.simulate(median, observe=False).predicted["image"]
        data = np.asarray(observed.values, dtype=float)
        model = np.asarray(predicted.values, dtype=float)
        residual = data - model
        extent = (
            float(observed.y.values[0]),
            float(observed.y.values[-1]),
            float(observed.x.values[0]),
            float(observed.x.values[-1]),
        )
        figure, axes = plt.subplots(1, 3, figsize=(11.0, 3.6), constrained_layout=True)
        for axis, panel, title in zip(
            axes, (data, model, residual), ("observed", "model (median)", "residual"), strict=True
        ):
            limit = float(np.max(np.abs(panel)))
            image = axis.imshow(
                panel,
                origin="lower",
                extent=extent,
                cmap="RdBu_r" if title == "residual" else "viridis",
                vmin=-limit if title == "residual" else None,
                vmax=limit if title == "residual" else None,
            )
            axis.set_title(f"{title} — {arm}")
            axis.set_xlabel("y (mas)")
            axis.set_ylabel("x (mas)")
            figure.colorbar(image, ax=axis, shrink=0.85)
        save(figure, f"panels_{arm}")
    return written


def save_calibration_figures(
    calibrations: Mapping[str, Any],
    directory: str | pathlib.Path,
    *,
    suffix: str = "png",
) -> list[pathlib.Path]:
    """SBC ranks and coverage per arm: parameter space, so no coordinate wall."""
    from ampere.results import plot_coverage, plot_sbc_ranks

    _, written, save = _writer(directory, suffix)
    for arm, calibration in calibrations.items():
        save(plot_sbc_ranks(calibration).get_figure(), f"sbc_ranks_{arm}")
        save(plot_coverage(calibration).get_figure(), f"coverage_{arm}")
    return written
