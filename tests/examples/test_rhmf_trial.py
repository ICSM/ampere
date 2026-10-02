"""``examples/rhmf_trial``: the adapter's shapes, and the ``--quick`` preset end to end (W6.9).

The trial needs ``robusta-hmf`` (the non-default ``rhmf`` extra, run in the
``rhmf`` pixi environment), so every row that fits skips cleanly where it is
absent — which is every other environment, including ``dev``. The one row that
runs *only* where it is absent pins the other half of the contract: without the
extra the adapter names it rather than failing with a bare ``ImportError``.

Nothing here claims the screen *works*; ``diagnostics.md`` §7 ("Amended W6.9")
carries the findings. The rows pin the adapter's shapes and the one direction
the quick preset is stable on: a sharp injected line is flagged, against the
controls, at its best grid point.
"""

from __future__ import annotations

import importlib.util
import pathlib

import numpy as np
import pytest

from ampere.core import AnomalyScore
from ampere.core.exceptions import OptionalDependencyError

from examples.rhmf_trial import trial

HAVE_RHMF = importlib.util.find_spec("robusta_hmf") is not None
needs_rhmf = pytest.mark.skipif(not HAVE_RHMF, reason="needs the rhmf extra (robusta-hmf)")


def _collection() -> trial.SpectraCollection:
    return trial.spectra_collection(controls=4, copies=2, seed=7)


@pytest.mark.skipif(HAVE_RHMF, reason="the extra is installed")
def test_without_the_extra_the_adapter_names_it() -> None:
    with pytest.raises(OptionalDependencyError, match="rhmf"):
        trial.require_rhmf()


@needs_rhmf
class TestAdapter:
    def test_to_matrix_shapes_and_weights(self) -> None:
        collection = _collection()
        assert collection.Y.shape == (4 + 2 * 4, 200)
        assert collection.W.shape == collection.Y.shape
        sigma = np.asarray(collection.spectra[0].uncertainty)
        np.testing.assert_allclose(collection.W[0], 1.0 / sigma**2)
        assert collection.wavelength.shape == (200,)

    def test_a_mask_becomes_zero_weight(self) -> None:
        from ampere.core import Spectrum

        spectrum = _collection().spectra[0]
        mask = np.zeros(spectrum.values.shape, dtype=bool)
        mask[10:13] = True
        masked = Spectrum(
            spectrum.spectral_axis.values * spectrum.spectral_axis.unit,
            spectrum.values * spectrum.unit,
            uncertainty=spectrum.uncertainty * spectrum.unit,
            mask=mask,
        )
        _, weights, _ = trial.to_matrix([masked])
        assert np.all(weights[0, 10:13] == 0.0)
        assert np.all(weights[0, :10] > 0.0)

    def test_misaligned_spectra_are_refused_not_resampled(self) -> None:
        from examples.m2_misspecification.generators import generate

        with pytest.raises(ValueError, match="common grid"):
            trial.to_matrix(
                [generate("none", size=200).container(), generate("none", size=100).container()]
            )

    def test_rank_and_scale_have_no_default(self) -> None:
        collection = _collection()
        with pytest.raises(TypeError):
            trial.fit_rhmf(collection.Y, collection.W)  # type: ignore[call-arg]

    def test_the_score_carries_its_provenance_and_caveat(self) -> None:
        collection = _collection()
        fit = trial.fit_rhmf(collection.Y, collection.W, rank=2, robust_scale=3.0, max_iter=60)
        score = trial.anomaly_score(fit, collection.wavelength, row=0)
        assert isinstance(score, AnomalyScore)
        assert score.provenance == "rhmf_prefit"
        assert "rank 2" in score.interpretation_notes
        assert trial.COMMIT[:12] in score.interpretation_notes
        assert "Nothing has been fitted" in score.interpretation_notes
        assert score.values.shape == (200,)
        assert np.all((score.values >= 0.0) & (score.values < 1.0))
        assert score.mask is not None and not score.mask.any()

    def test_a_per_object_score_is_one_value_per_row(self) -> None:
        collection = _collection()
        fit = trial.fit_rhmf(collection.Y, collection.W, rank=1, robust_scale=3.0, max_iter=60)
        score = trial.anomaly_score(fit, collection.wavelength, aggregate=0.1)
        assert score.values.shape == (collection.Y.shape[0],)
        with pytest.raises(ValueError, match="exactly one"):
            trial.anomaly_score(fit, collection.wavelength)

    def test_the_renderer_accepts_it(self) -> None:
        from matplotlib import pyplot

        from ampere.results import plot_anomaly_score

        collection = _collection()
        fit = trial.fit_rhmf(collection.Y, collection.W, rank=1, robust_scale=3.0, max_iter=60)
        axes = plot_anomaly_score(trial.anomaly_score(fit, collection.wavelength, row=0))
        try:
            assert "rhmf_prefit" in axes.get_legend_handles_labels()[1]
        finally:
            pyplot.close(axes.get_figure())


@needs_rhmf
def test_the_quick_preset_runs_and_flags_the_sharp_line(tmp_path: pathlib.Path) -> None:
    from examples.rhmf_trial import __main__ as runner

    assert runner.main(["--quick", "--out", str(tmp_path)]) == 0
    assert (tmp_path / "spectra_grid.csv").exists() and (tmp_path / "image_grid.csv").exists()
    assert (tmp_path / "spectra_strong_sharp.png").exists()
    import csv

    with (tmp_path / "spectra_grid.csv").open() as handle:
        sharp = [row for row in csv.DictReader(handle) if row["scenario"] == "strong_sharp"]
    assert max(float(row["excess"]) for row in sharp) > 2.0
