"""``examples/modified_blackbody`` builds, negotiates and fits end to end (W6.13 (3)).

Fast and always on: the ten AKARI/Herschel bands all load from ampere's
bundled filter library; the problem builds and negotiates one dataset; a
tiny-budget fit runs on each of the four engines (the ``sbi`` arm skips
without the ``sbi`` extra installed, exactly as
:mod:`tests.examples.test_npe_native` does); ``--all`` at tiny budgets writes
the posterior-predictive overlay to a ``tmp_path``; and ``main`` prints a
report. The full-budget 95 % recovery this item is accepted on is **not**
part of this suite, for the reason :mod:`tests.examples.test_sed_composition`
gives at length: it is verified by running
``python -m examples.modified_blackbody`` directly, and
:mod:`examples.modified_blackbody.modified_blackbody`'s own docstring has the
numbers from doing exactly that, once per engine.
"""

from __future__ import annotations

import importlib.util

import pytest

import examples.modified_blackbody.modified_blackbody as modified_blackbody_module
from examples.modified_blackbody.modified_blackbody import (
    ENGINES,
    FILTERS,
    QUALIFIED_TRUTH,
    available_engines,
    build_instrument,
    build_model,
    build_problem,
    fit,
    main,
    overlay_figure,
    recovers_truth,
    report,
    run_all,
)

HAS_SBI = importlib.util.find_spec("sbi") is not None
needs_sbi = pytest.mark.skipif(
    not HAS_SBI, reason="needs the 'sbi' extra (pixi run -e sbi ...), which brings in torch too"
)

TINY_ENSEMBLE = {"walkers": 8, "steps": 20, "burn_in": 5}
#: dynesty needs enough live points to bound a four-dimensional posterior
#: efficiently; fewer than this makes the ellipsoidal bootstrap pathologically
#: inefficient (dynesty's own warning) rather than merely coarse.
TINY_NESTED = {"live_points": 25, "dlogz": 5.0}


class TestTheFiltersLoad:
    """All ten AKARI/Herschel bands load from ampere's bundled filter library."""

    def test_every_filter_loads(self) -> None:
        instrument = build_instrument()
        assert instrument.label == "catalogue"

    def test_ten_bands_named(self) -> None:
        assert len(FILTERS) == 10
        assert all(name.startswith(("AKARI", "HERSCHEL")) for name in FILTERS)


class TestTheProblemBuilds:
    """One model, one dataset, the four parameters as free parameters."""

    def test_build_problem_composes_one_dataset(self) -> None:
        problem = build_problem()
        assert problem.backend == "reference"
        assert set(problem.parameters.free_names) == set(QUALIFIED_TRUTH)

    def test_model_alone_evaluates_on_its_own_grid(self) -> None:
        from examples.modified_blackbody.modified_blackbody import GRID

        model = build_model()
        result = model(temperature=30.0, logmass=1.0, beta=-2.0, distance=0.11)
        assert result.single().values.size == GRID.size

    def test_available_engines_excludes_sbi_without_the_extra(self) -> None:
        engines = available_engines()
        assert set(engines).issubset(set(ENGINES))
        if not HAS_SBI:
            assert "sbi" not in engines


class TestTheFitRuns:
    """A tiny-budget fit, end to end, on each of the three always-on engines."""

    def test_emcee_arm_returns_the_expected_parameters(self) -> None:
        problem = build_problem()
        run = fit(problem, engine="emcee", **TINY_ENSEMBLE)
        posterior = run["posterior"].dataset
        assert set(posterior.data_vars) == set(QUALIFIED_TRUTH)
        assert posterior.sizes["chain"] == TINY_ENSEMBLE["walkers"]

    def test_report_and_recovers_truth_name_every_parameter(self) -> None:
        problem = build_problem()
        run = fit(problem, engine="emcee", **TINY_ENSEMBLE)
        text = report(run)
        covered = recovers_truth(run)
        assert set(covered) == set(QUALIFIED_TRUTH)
        for name in QUALIFIED_TRUTH:
            assert name in text

    def test_zeus_arm_runs(self) -> None:
        problem = build_problem()
        run = fit(problem, engine="zeus", **TINY_ENSEMBLE)
        posterior = run["posterior"].dataset
        assert set(posterior.data_vars) == set(QUALIFIED_TRUTH)

    def test_dynesty_arm_runs(self) -> None:
        problem = build_problem()
        run = fit(problem, engine="dynesty", **TINY_NESTED)
        posterior = run["posterior"].dataset
        assert set(posterior.data_vars) == set(QUALIFIED_TRUTH)


@needs_sbi
class TestTheSBIArm:
    """The sbi engine, at a tiny simulation budget."""

    def test_sbi_arm_runs(self) -> None:
        problem = build_problem()
        run = fit(problem, engine="sbi", budget=32, draws=16)
        posterior = run["posterior"].dataset
        assert set(posterior.data_vars) == set(QUALIFIED_TRUTH)


class TestTheOverlay:
    """``--all`` at tiny budgets, and the posterior-predictive overlay figure.

    dynesty's own smoke row (:class:`TestTheFitRuns`) is the one place in
    this file a nested-sampling run actually happens -- ``available_engines``
    is monkeypatched to ``("emcee", "zeus")`` here, both a few seconds each,
    so the overlay mechanism (several engines' runs, one composite figure) is
    checked without paying for a second and third dynesty run under a minute
    budget. ``run_all``/``available_engines`` are the module's own dispatch,
    so patching the module-level name is what actually takes effect.
    """

    def test_run_all_and_the_overlay_figure(self, monkeypatch, tmp_path) -> None:
        monkeypatch.setattr(
            modified_blackbody_module, "available_engines", lambda: ("emcee", "zeus")
        )
        problem = build_problem()
        runs = run_all(problem, **TINY_ENSEMBLE, **TINY_NESTED)
        assert set(runs) == {"emcee", "zeus"}

        figure = overlay_figure(runs, problem)
        assert len(figure.axes) == len(runs)
        destination = tmp_path / "overlay.png"
        figure.savefig(destination)
        assert destination.exists()
        assert destination.stat().st_size > 0


class TestMain:
    """The CLI entry point, at the tiny budget, plain and with ``--all``."""

    def test_main_prints_a_report(self, capsys: pytest.CaptureFixture[str]) -> None:
        assert (
            main(
                [
                    "--engine",
                    "emcee",
                    "--walkers",
                    str(TINY_ENSEMBLE["walkers"]),
                    "--steps",
                    str(TINY_ENSEMBLE["steps"]),
                    "--burn-in",
                    str(TINY_ENSEMBLE["burn_in"]),
                ]
            )
            == 0
        )
        out = capsys.readouterr().out
        assert "free parameters:" in out
        assert "emcee on reference" in out
        assert "wall clock" in out

    def test_main_all_writes_the_overlay(
        self, monkeypatch, tmp_path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        # As TestTheOverlay: emcee only, so --all's CLI plumbing (mkdir,
        # savefig, the "N engine(s)" summary line) is checked without a
        # second dynesty run under this file's minute budget.
        monkeypatch.setattr(modified_blackbody_module, "available_engines", lambda: ("emcee",))
        assert (
            main(
                [
                    "--all",
                    "--walkers",
                    str(TINY_ENSEMBLE["walkers"]),
                    "--steps",
                    str(TINY_ENSEMBLE["steps"]),
                    "--burn-in",
                    str(TINY_ENSEMBLE["burn_in"]),
                    "--figures",
                    str(tmp_path),
                ]
            )
            == 0
        )
        out = capsys.readouterr().out
        assert "engine(s)" in out
        assert (tmp_path / "modified_blackbody_overlay.png").exists()
