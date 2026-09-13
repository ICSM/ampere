"""``BinaryWithDisc`` composes the way :mod:`examples.interferometry.model` claims to.

Cheap, structural checks that have to hold before any of the science in
``tests/interferometry/test_calibration.py`` means anything: the composite
declares the union of its children's capabilities and exactly the binary's
free parameters, negotiates both children onto one grid, and its brightness
really is the two children's summed.
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.stats as st

import ampere.backends.reference.interferometry as itf
from ampere.core import negotiate
from examples.interferometry import generators as gen
from examples.interferometry import model as m


def _instrument():
    import interferometry_fixtures as fixtures

    return fixtures.chain(itf, fixtures.visibilities(), "vis")


class TestComposition:
    def test_free_parameters_are_the_binarys_alone(self) -> None:
        grid = gen.seed_grid()
        composite = m.binary_with_disc(
            itf,
            grid,
            grid,
            disc_flux=gen.DISC_FLUX,
            disc_fwhm=gen.DISC_FWHM,
            separation=st.uniform(6.0, 14.0),
            flux_ratio=st.uniform(0.1, 0.7),
            **{k: v for k, v in gen.BINARY.items() if k not in ("separation", "flux_ratio")},
            component_fwhm=gen.COMPONENT_FWHM,
        )
        assert composite.parameters.free_names == ("separation", "flux_ratio")
        # The disc's own parameters are declared nowhere on the composite: it
        # would otherwise be a third and fourth thing the sampler could move.
        assert "flux" not in [n for n in composite.parameters.names if "disc" in n]

    def test_capability_flags_are_the_binarys(self) -> None:
        grid = gen.seed_grid()
        binary = itf.Binary(grid, grid, **gen.BINARY, component_fwhm=gen.COMPONENT_FWHM)
        composite = gen.truth_model(itf)
        assert composite.DIFFERENTIABLE == binary.DIFFERENTIABLE
        assert composite.BATCHABLE == binary.BATCHABLE
        assert composite.BACKEND == binary.BACKEND

    def test_evaluate_is_the_sum_of_both_children(self) -> None:
        """The composite's Image really is binary-brightness plus disc-brightness."""
        grid = gen.seed_grid()
        composite = gen.truth_model(itf)
        instrument = _instrument()
        compiled = composite.compile_for(negotiate([instrument]))
        result = compiled.evaluate()
        combined_image = result["sky"]

        binary_alone = itf.Binary(grid, grid, **gen.BINARY, component_fwhm=gen.COMPONENT_FWHM)
        disc_alone = itf.GaussianSource(grid, grid, fwhm=gen.DISC_FWHM, flux=gen.DISC_FLUX)
        requirements = negotiate([instrument])
        binary_alone = binary_alone.compile_for(requirements)
        disc_alone = disc_alone.compile_for(requirements)
        expected = binary_alone.evaluate()["sky"].values + disc_alone.evaluate()["sky"].values

        np.testing.assert_allclose(combined_image.values, expected, rtol=1e-12)

    def test_both_children_adopt_the_same_negotiated_grid(self) -> None:
        composite = gen.truth_model(itf)
        instrument = _instrument()
        compiled = composite.compile_for(negotiate([instrument]))
        np.testing.assert_array_equal(
            compiled._binary.templates["sky"].x.values, compiled._disc.templates["sky"].x.values
        )
        np.testing.assert_array_equal(
            compiled._binary.templates["sky"].y.values, compiled._disc.templates["sky"].y.values
        )

    def test_disc_is_fixed_at_its_declared_truth(self) -> None:
        composite = gen.truth_model(itf)
        assert composite._disc.parameters["flux"].fixed
        assert composite._disc.parameters["fwhm"].fixed
        assert composite._disc.parameters["flux"].value == pytest.approx(gen.DISC_FLUX)
        assert composite._disc.parameters["fwhm"].value == pytest.approx(gen.DISC_FWHM)
