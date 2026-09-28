"""``examples/linear_sed`` builds, negotiates and fits end to end (W6.13 (1)).

Fast and always on, like :mod:`tests.examples.test_sed_composition`: the IRS
grids this package's :mod:`~examples.linear_sed.generators` reproduces from
the tracked CASSIS file are checked against the legacy reader's own; the
problem builds and negotiates three datasets on one channel; a tiny-budget
fit runs on each of the four engines (the ``sbi`` arms skip without the
``sbi`` extra installed, exactly as :mod:`tests.examples.test_npe_native`
does); and ``main`` prints a report. The full-budget 95 % recovery this item
is accepted on is **not** part of this suite, for the same reason
:mod:`tests.examples.test_sed_composition` gives: it is one to two orders of
magnitude slower than everything else here, and is instead verified by
running ``python -m examples.linear_sed`` directly -- the branch report and
:mod:`examples.linear_sed.linear_sed`'s own docstring have the numbers from
doing exactly that.
"""

from __future__ import annotations

import importlib.util

import numpy as np
import pytest

from examples.linear_sed import generators
from examples.linear_sed.linear_sed import (
    EMBEDDING,
    FILTERS,
    GRID,
    QUALIFIED_TRUTH,
    build_instruments,
    build_model,
    build_problem,
    fit,
    main,
    recovers_truth,
    report,
)

HAS_SBI = importlib.util.find_spec("sbi") is not None
needs_sbi = pytest.mark.skipif(
    not HAS_SBI, reason="needs the 'sbi' extra (pixi run -e sbi ...), which brings in torch too"
)

TINY = {"walkers": 8, "steps": 20, "burn_in": 5}


class TestTheIRSGrids:
    """:func:`generators.irs_wavelength_grids` reproduces the legacy reader exactly."""

    def test_grids_match_the_legacy_reader(self) -> None:
        from ampere.legacy.data import Spectrum

        sl, ll = generators.irs_wavelength_grids()
        legacy = Spectrum.fromFile(str(generators.IRS_FILE), format="SPITZER-YAAAR")
        assert np.array_equal(sl, legacy[0].wavelength)
        assert np.array_equal(ll, legacy[1].wavelength)

    def test_deduplicated_grids_are_strictly_increasing(self) -> None:
        sl, ll = generators.deduplicated_grids()
        assert sl.size < 200  # the raw SL grid has 17 repeated samples
        assert ll.size < 187  # the raw LL grid has 10 repeated samples
        assert np.all(np.diff(sl) > 0.0)
        assert np.all(np.diff(ll) > 0.0)


class TestTheProblemBuilds:
    """One model, three datasets (catalogue, sl, ll), distinct labels."""

    def test_build_problem_composes_three_datasets(self) -> None:
        problem = build_problem()
        assert problem.backend == "reference"
        for name in QUALIFIED_TRUTH:
            assert name in problem.parameters.free_names

    def test_build_instruments_have_distinct_labels(self) -> None:
        catalogue, sl_instrument, ll_instrument = build_instruments()
        labels = {catalogue.label, sl_instrument.label, ll_instrument.label}
        assert labels == {"catalogue", "sl", "ll"}

    def test_no_gp_problem_has_no_gp_parameters(self) -> None:
        problem = build_problem(gp=False)
        assert not any("likelihood" in name for name in problem.parameters.free_names)

    def test_model_alone_evaluates_on_its_own_grid(self) -> None:
        model = build_model()
        result = model(slope=1.0, intercept=1.0)
        assert result.single().values.size == GRID.size

    def test_catalogue_filters_are_the_legacy_two_bands(self) -> None:
        assert FILTERS == ("WISE_RSR_W1", "SPITZER_MIPS_70")


class TestTheFitRuns:
    """A tiny-budget fit, end to end, returns the four expected parameters.

    ``gp=False`` throughout: the default (GP) problem has eight free
    parameters (the four qualified ones plus a Matern32 amplitude and length
    scale per spectrum), which needs at least sixteen walkers
    (``ampere.inference``'s ``2 * n_dim`` ensemble floor) --
    :data:`TINY`'s eight is sized for the four-dimensional ``gp=False``
    problem, which these tests are about (the four qualified parameters'
    coverage), not the GP arm.
    """

    def test_fit_returns_the_expected_parameters(self) -> None:
        problem = build_problem(gp=False)
        run = fit(problem, engine="emcee", **TINY)
        posterior = run["posterior"].dataset
        assert set(QUALIFIED_TRUTH).issubset(set(posterior.data_vars))
        assert posterior.sizes["chain"] == TINY["walkers"]
        assert posterior.sizes["draw"] == TINY["steps"] - TINY["burn_in"]

    def test_report_names_every_qualified_parameter(self) -> None:
        problem = build_problem(gp=False)
        run = fit(problem, engine="emcee", **TINY)
        text = report(run)
        assert "emcee on reference" in text
        for name in QUALIFIED_TRUTH:
            assert name in text

    def test_recovers_truth_reports_one_bool_per_qualified_parameter(self) -> None:
        problem = build_problem(gp=False)
        run = fit(problem, engine="emcee", **TINY)
        covered = recovers_truth(run)
        assert set(covered) == set(QUALIFIED_TRUTH)
        assert all(isinstance(value, bool) for value in covered.values())

    def test_embedding_without_sbi_is_refused(self) -> None:
        problem = build_problem(gp=False)
        with pytest.raises(SystemExit, match="--embedding only applies to --engine sbi"):
            fit(problem, engine="emcee", embedding=True, **TINY)


class TestTheOtherEngines:
    """The zeus and dynesty arms, both in ``dev``, at tiny budgets."""

    def test_zeus_arm_runs(self) -> None:
        problem = build_problem(gp=False)
        run = fit(problem, engine="zeus", walkers=8, steps=20, burn_in=5)
        posterior = run["posterior"].dataset
        assert set(QUALIFIED_TRUTH).issubset(set(posterior.data_vars))

    def test_dynesty_arm_runs(self) -> None:
        problem = build_problem(gp=False)
        run = fit(problem, engine="dynesty", live_points=25)
        posterior = run["posterior"].dataset
        assert set(QUALIFIED_TRUTH).issubset(set(posterior.data_vars))


@needs_sbi
class TestTheSBIArm:
    """The sbi engine, at a tiny simulation budget, with and without --embedding."""

    def test_sbi_arm_runs(self) -> None:
        problem = build_problem(gp=False)
        run = fit(problem, engine="sbi", budget=32, draws=16)
        posterior = run["posterior"].dataset
        assert set(QUALIFIED_TRUTH).issubset(set(posterior.data_vars))

    def test_sbi_arm_with_embedding_runs(self) -> None:
        problem = build_problem(gp=False)
        run = fit(problem, engine="sbi", budget=32, draws=16, embedding=True)
        assert run.attrs["ampere_engine"] == "sbi"

    def test_the_embedding_dict_is_the_legacy_one(self) -> None:
        assert EMBEDDING == {
            "type": "FC",
            "num_hiddens": 100,
            "n_layers": 3,
            "output_dim": 20,
        }


class TestMain:
    """The CLI entry point, at the tiny budget."""

    def test_main_prints_a_report(self, capsys: pytest.CaptureFixture[str]) -> None:
        assert (
            main(
                [
                    "--no-gp",
                    "--walkers",
                    str(TINY["walkers"]),
                    "--steps",
                    str(TINY["steps"]),
                    "--burn-in",
                    str(TINY["burn_in"]),
                ]
            )
            == 0
        )
        out = capsys.readouterr().out
        assert "free parameters:" in out
        assert "emcee on reference" in out
        assert "wall clock" in out
