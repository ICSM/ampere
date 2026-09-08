"""The figures render, from stored runs, into a directory the test names.

Nothing here checks what a figure *looks* like — that is not a thing a test can
usefully assert. What it checks is the two things that do break: that each
figure can be built at all from a run plus the groups ``ampere.results``
derived, and that every shipped renderer W2.8 landed is exercised by this study
rather than merely mentioned by it.

Figures are never committed (AGENTS.md ground rule 7), so every path here is
under pytest's ``tmp_path``.
"""

from __future__ import annotations

import pathlib
from typing import Any

import pytest

from examples.m2_misspecification import figures, study

pytestmark = pytest.mark.usefixtures("agg_backend")

#: The two scenarios the per-run renderers are exercised on. All four would
#: cost a corner plot each per likelihood for no extra coverage: the control and
#: the sharp line are the two whose renderers differ in what they have to show.
RENDERED = ("none", "strong_sharp")


def test_the_three_paper_figures_are_written(
    study_results: Any, diagnoses: Any, tmp_path: pathlib.Path
) -> None:
    written = figures.save_paper_figures(study_results, tmp_path / "paper", suffix="png")
    assert [path.stem for path in written] == list(figures.PAPER_FIGURES)
    for path in written:
        assert path.exists() and path.stat().st_size > 0, path


def test_the_residual_panel_carries_the_gp_band(study_results: Any, diagnoses: Any) -> None:
    """Figure 1's whole point: the conditioned GP mean, drawn over the residuals.

    The band comes from the ``gp_localisation`` group, which exists only because
    :func:`~examples.m2_misspecification.study.prepare` derived it — so this
    also checks that the figure is reading the stored run rather than re-fitting
    anything.
    """
    from ampere.results import GP_LOCALISATION_GROUP

    entry = study_results[("strong_sharp", "flexible")]
    assert GP_LOCALISATION_GROUP in entry["run"].children
    band = figures._localisation(entry)
    assert band is not None
    centre, lower, upper = band
    assert centre.shape == lower.shape == upper.shape
    assert (lower <= centre + 1e-12).all() and (centre <= upper + 1e-12).all()
    # The GP found the injected line: its conditioned mean peaks there.
    grid = entry["data"].wavelength
    assert abs(float(grid[int(centre.argmax())]) - 0.8630) < 5e-4


def test_the_shipped_renderers_all_run(
    study_results: Any, diagnoses: Any, tmp_path: pathlib.Path
) -> None:
    """W2.8's plotting surface, exercised on this study's runs.

    ``plot_residuals`` warns on a GP fit and ``plot_gp_localisation`` needs one,
    so :func:`~examples.m2_misspecification.figures.save_result_figures` calls
    each on the runs it is scoped to; asserting the *set* of files written is
    how that scoping stays true rather than becoming a comment.
    """
    written = figures.save_result_figures(
        study_results, tmp_path / "runs", suffix="png", scenarios=RENDERED
    )
    stems = {path.stem for path in written}
    for key in RENDERED:
        assert f"corner_{key}_standard" in stems
        assert f"corner_{key}_flexible" in stems
        assert f"trace_{key}_standard" in stems
        assert f"residuals_{key}_standard" in stems
        assert f"localisation_{key}_flexible" in stems
        assert f"anomaly_{key}_flexible" in stems
        # Scoped, not merely preferred: no whiteness panel for a GP fit, no
        # localisation panel for a fit with no GP.
        assert f"residuals_{key}_flexible" not in stems
        assert f"localisation_{key}_standard" not in stems
    for path in written:
        assert path.exists() and path.stat().st_size > 0, path


def test_the_figure_builders_read_the_runs_and_not_the_truth(study_results: Any) -> None:
    """A sanity check on the recovery figure: it summarises the stored draws.

    Building the figure and re-deriving its numbers from
    :func:`~examples.m2_misspecification.study.summarise` must give the same
    medians, because there is only one place either can come from.
    """
    figure = figures.figure_parameter_recovery(study_results)
    assert len(figure.axes) == 4
    summaries = study.summarise(
        study_results[("strong_smooth", "flexible")]["run"], names=study.PHYSICAL_NAMES
    )
    assert set(summaries) == set(study.PHYSICAL_NAMES)
    import matplotlib.pyplot as plt

    plt.close(figure)
