"""``examples/ngc6302`` builds, negotiates and fits end to end (W6.13 (2)).

This first commit's rows: the opacity buffers load with the right shape and
the two unit conversions, and one evaluation equals the legacy
``SpectrumNGC6302`` at a fixed theta to ``rtol=1e-10``. The data, problem,
engine and dust-mass rows follow in later commits on this branch (see
:mod:`examples.ngc6302.ngc6302`'s own module docstring).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from examples.ngc6302 import generators
from examples.ngc6302.ngc6302 import DEFAULT_GRID, build_model


class TestOpacityBuffers:
    """The eight opacity tables load as sixteen buffers, converted once."""

    def test_every_species_has_a_wavelength_and_opacity_buffer(self) -> None:
        model = build_model()
        for name in generators.SPECIES:
            assert f"{name}_wavelength" in model.buffers
            assert f"{name}_opacity" in model.buffers
            wavelength = model.buffers[f"{name}_wavelength"].array
            opacity = model.buffers[f"{name}_opacity"].array
            assert wavelength.shape == opacity.shape
            assert wavelength.ndim == 1 and wavelength.size > 0

    def test_enstatite_and_diopside_are_converted_by_1e_4(self) -> None:
        model = build_model()
        filenames = np.loadtxt(generators.OPACITY_FILE_LIST, dtype=str)
        for name, filename in zip(generators.SPECIES, filenames, strict=True):
            if name not in ("enstatite", "diopside"):
                continue
            raw_opacity = np.loadtxt(
                generators.OPACITY_DIRECTORY / str(filename), comments="#", unpack=True
            )[1]
            np.testing.assert_allclose(model.buffers[f"{name}_opacity"].array, raw_opacity * 1e-4)

    def test_calcite_and_dolomite_are_converted_by_the_density_factor(self) -> None:
        model = build_model()
        factors = {"calcite": 2.71 * (4.0 / 3.0) * 1e-4, "dolomite": 2.87 * (4.0 / 3.0) * 1e-4}
        filenames = np.loadtxt(generators.OPACITY_FILE_LIST, dtype=str)
        for name, filename in zip(generators.SPECIES, filenames, strict=True):
            if name not in factors:
                continue
            raw_opacity = np.loadtxt(
                generators.OPACITY_DIRECTORY / str(filename), comments="#", unpack=True
            )[1]
            np.testing.assert_allclose(
                model.buffers[f"{name}_opacity"].array, raw_opacity * factors[name]
            )


class TestModelEqualsLegacy:
    """One evaluation, at a fixed theta, equals ``examples.NGC6302.SpectrumNGC6302``."""

    def test_evaluate_matches_legacy_model(self, monkeypatch: pytest.MonkeyPatch) -> None:
        import examples.NGC6302 as legacy_module

        # The legacy __init__ resolves its opacity directory from os.getcwd(),
        # and __call__ references a bare global `wavelengths` (not
        # self.wavelength) that only exists when the script runs as
        # __main__ -- both reproduced here, inside the test only (ruling 1).
        monkeypatch.chdir(Path(__file__).resolve().parents[2] / "examples")
        legacy_module.wavelengths = DEFAULT_GRID
        legacy_model = legacy_module.SpectrumNGC6302(DEFAULT_GRID)
        legacy_model(
            generators.TRUTH["logacold0"],
            generators.TRUTH["logacold1"],
            generators.TRUTH["logacold2"],
            generators.TRUTH["logacold3"],
            generators.TRUTH["logacold4"],
            generators.TRUTH["logacold6"],
            generators.TRUTH["logacold7"],
            generators.TRUTH["logawarm1"],
            generators.TRUTH["logawarm2"],
            generators.TRUTH["logawarm5"],
            generators.TRUTH["logawarm7"],
            generators.TRUTH["Tcold0"],
            generators.TRUTH["Tcold1"],
            generators.TRUTH["Twarm0"],
            generators.TRUTH["Twarm1"],
        )
        legacy_flux = legacy_model.modelFlux

        model = build_model()
        result = model(**generators.TRUTH).single().values
        np.testing.assert_allclose(result, legacy_flux, rtol=1e-10)
