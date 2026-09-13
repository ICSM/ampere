"""``examples/sed_composition`` builds, negotiates, and fits end to end (W4.11).

Full-budget recovery — the central-95 % coverage W4.11 is accepted on — is
**not** part of this suite. Verifying it means sampling at the script's
default budget, which is designed to finish under two minutes on ``dev`` but
is still one to two orders of magnitude slower than everything else this
suite runs, and milestone M2's own full ladder sets the precedent for keeping
that kind of run out of the per-PR gate (`tests/m2/conftest.py`). Unlike that
ladder, gating it behind an opt-in pytest marker would mean registering the
marker in ``pyproject.toml``, which is outside this item's file ownership —
so instead of a slow, always-collected-but-skipped test, there is simply no
pytest coverage of the full budget: it is verified by running
``python -m examples.sed_composition`` directly, and the branch report quotes
the recovered intervals and the wall-clock time from doing exactly that.

What *is* here, and fast: the composition builds and negotiates; the
two-instruments-one-channel label collision the tutorial page walks through
is reproduced word for word; a tiny-budget fit runs end to end on the
reference backend and returns the four expected parameters; ``main`` prints
the negotiated requirements and a report.
"""

from __future__ import annotations

import astropy.units as u
import numpy as np
import pytest

from ampere.backends.reference import LSFConvolution, Resample, SyntheticPhotometry
from ampere.core import (
    Dataset,
    DatasetCollection,
    DatasetError,
    Instrument,
    PhotometricPoints,
    Spectrum,
)
from examples.sed_composition import generators
from examples.sed_composition.sed_composition import (
    FILTERS,
    PHOTOMETRY_TABULATION,
    QUALIFIED_TRUTH,
    build_instruments,
    build_model,
    build_problem,
    fit,
    main,
    recovers_truth,
    report,
)

TINY = {"walkers": 8, "steps": 20, "burn_in": 5}


class TestTheCompositionBuilds:
    """One model, two datasets, distinct labels, one negotiated channel."""

    def test_build_problem_composes_two_datasets(self) -> None:
        problem = build_problem("reference")
        assert problem.backend == "reference"
        assert set(problem.datasets) == {"irs", "catalogue"}
        assert set(problem.parameters.free_names) == set(QUALIFIED_TRUTH)

    def test_the_channel_is_negotiated_from_both_instruments(self) -> None:
        problem = build_problem("reference")
        requirements = problem.requirements["model"]["sed"]
        assert set(requirements.sources) == {"irs", "catalogue"}
        assert "spectral_axis" in requirements.axes
        # The photometry step's exact tabulation must have made it into the
        # union grid -- this is finding 1's whole point.
        axis = requirements.axes["spectral_axis"]
        assert axis.points is not None
        assert axis.points.size == PHOTOMETRY_TABULATION.size

    def test_build_model_and_instruments_are_reusable_alone(self) -> None:
        model = build_model("reference")
        spectrograph, camera = build_instruments("reference")
        assert spectrograph.label == "irs"
        assert camera.label == "catalogue"
        assert spectrograph.channel == camera.channel == "sed"
        observed_spectrum, observed_photometry = generators.synthetic_data(
            model, spectrograph, camera, seed=generators.SEED
        )
        assert observed_spectrum.values.size == spectrograph.steps[1].buffers["target"].value.size
        assert tuple(observed_photometry.filters) == FILTERS


class TestTheLabelCollision:
    """Two instruments on one channel, left to their default labels, collide by name.

    This is the exact scenario :doc:`the tutorial page </sed_composition>`
    walks through -- reproduced here so the quoted error text cannot drift
    from what the code actually says.
    """

    def test_two_default_labelled_instruments_collide(self) -> None:
        grid = np.linspace(5.0, 35.0, 12)
        spectrograph = Instrument(
            [LSFConvolution(resolving_power=100.0), Resample(grid)], channel="sed"
        )
        camera = Instrument(
            [SyntheticPhotometry.from_library(["2MASS_J"], np.geomspace(1.0, 10.0, 50))],
            channel="sed",
        )
        assert spectrograph.label == camera.label == "sed"

        observed_spectrum = Spectrum(
            grid * u.um, np.ones_like(grid) * u.Jy, uncertainty=0.1 * np.ones_like(grid) * u.Jy
        )
        observed_photometry = PhotometricPoints(
            ["2MASS_J"], [1.235] * u.um, [1.0] * u.Jy, uncertainty=[0.1] * u.Jy
        )

        with pytest.raises(DatasetError, match="two datasets are labelled 'sed'") as excinfo:
            DatasetCollection(
                [
                    Dataset(observed_spectrum, spectrograph),
                    Dataset(observed_photometry, camera),
                ]
            )
        assert "Instrument's own label" in str(excinfo.value)

    def test_the_fix_is_distinct_labels(self) -> None:
        spectrograph, camera = build_instruments("reference")
        assert spectrograph.label != camera.label
        # Should not raise: build_problem already does exactly this.
        build_problem("reference")


class TestTheFitRuns:
    """A tiny-budget fit, end to end, returns the four expected parameters."""

    def test_fit_returns_the_expected_parameters(self) -> None:
        problem = build_problem("reference")
        run = fit(problem, backend="reference", **TINY)
        posterior = run["posterior"].dataset
        assert set(posterior.data_vars) == set(QUALIFIED_TRUTH)
        assert posterior.sizes["chain"] == TINY["walkers"]
        assert posterior.sizes["draw"] == TINY["steps"] - TINY["burn_in"]

    def test_report_names_every_parameter(self) -> None:
        problem = build_problem("reference")
        run = fit(problem, backend="reference", **TINY)
        text = report(run)
        assert "emcee on reference" in text
        for name in QUALIFIED_TRUTH:
            assert name in text

    def test_recovers_truth_reports_one_bool_per_parameter(self) -> None:
        problem = build_problem("reference")
        run = fit(problem, backend="reference", **TINY)
        covered = recovers_truth(run)
        assert set(covered) == set(QUALIFIED_TRUTH)
        assert all(isinstance(value, bool) for value in covered.values())


class TestMain:
    """The CLI entry point, at the tiny budget."""

    def test_main_prints_requirements_and_a_report(
        self, capsys: pytest.CaptureFixture[str]
    ) -> None:
        assert (
            main(
                [
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
        assert "negotiated channels:" in out
        assert "sources asking of channel 'sed':" in out
        assert "free parameters:" in out
        assert "emcee on reference" in out
        assert "wall clock" in out
