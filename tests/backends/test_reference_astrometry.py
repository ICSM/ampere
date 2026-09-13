"""``ampere.backends.reference.astrometry``: the refusals, the messages, and a fit.

The conformance battery (``tests/conformance/test_astrometry.py``) holds the
arithmetic to the closed-form ephemeris, once per registered backend. This
file is the reference path's own: the construction-time refusals a user
actually meets, and the end-to-end fit — an emcee run on the synthetic orbit
that has to recover the injected proper motion, period, phase and reflex
amplitude it was made from (the item's own acceptance criterion, exercised
directly on this backend's pieces rather than through the example script).
"""

from __future__ import annotations

from typing import Any

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

import astrometry_fixtures as kit_data
from ampere.backends.reference import astrometry as reference_astrometry
from ampere.backends.reference.astrometry import EpochSample, ReflexOrbit
from ampere.core import (
    Dataset,
    DatasetCollection,
    FittingProblem,
    GaussianFamily,
    IndependentNoise,
    Instrument,
    Likelihood,
    TransformationError,
    negotiate,
)
from ampere.inference import EmceeEngine


class TestEpochSampleConstruction:
    """The instrument that takes its coordinates from the data, by construction."""

    def test_from_observed_takes_a_time_series_unchanged(self) -> None:
        observed = kit_data.observed_channel()
        step = EpochSample.from_observed(observed)
        (requirement,) = step.requirements()
        assert np.array_equal(requirement.points, observed.time.values)

    def test_an_empty_epoch_array_is_refused(self) -> None:
        with pytest.raises(TransformationError, match="at least one observed epoch"):
            EpochSample(np.array([]))

    def test_unsorted_epochs_are_refused_by_name(self) -> None:
        with pytest.raises(TransformationError, match="strictly increasing"):
            EpochSample(np.array([0.0, 50.0, 20.0]))

    def test_apply_is_the_identity(self) -> None:
        observed = kit_data.observed_channel()
        step = EpochSample.from_observed(observed)
        assert step(observed, None) is observed


class TestReflexOrbitConstruction:
    """The model's own default grid, checked at construction like an image model's."""

    def test_an_empty_time_grid_is_refused(self) -> None:
        with pytest.raises(ValueError, match="strictly increasing"):
            ReflexOrbit(np.array([]))

    def test_an_unsorted_time_grid_is_refused(self) -> None:
        with pytest.raises(ValueError, match="strictly increasing"):
            ReflexOrbit(np.array([0.0, 30.0, 10.0]))

    def test_evaluate_without_compile_for_uses_the_own_grid(self) -> None:
        model = ReflexOrbit(kit_data.EPOCHS * u.day, **kit_data.ORBIT)
        result = model.evaluate()
        assert result["ra"].time.values.size == kit_data.EPOCHS.size
        assert result["dec"].time.values.size == kit_data.EPOCHS.size

    def test_each_channel_may_be_compiled_onto_its_own_epochs(self) -> None:
        """A strict generalisation of the design sketch's shared-epoch example."""
        ra_epochs = kit_data.EPOCHS
        dec_epochs = kit_data.EPOCHS[:-2]
        ra_instrument = Instrument([EpochSample(ra_epochs)], channel="ra", label="astrom_ra")
        dec_instrument = Instrument([EpochSample(dec_epochs)], channel="dec", label="astrom_dec")
        model = ReflexOrbit(ra_epochs * u.day, **kit_data.ORBIT)
        compiled = model.compile_for(negotiate([ra_instrument, dec_instrument]))
        result = compiled.evaluate()
        assert result["ra"].time.values.size == ra_epochs.size
        assert result["dec"].time.values.size == dec_epochs.size


class TestTheSyntheticOrbitIsRecovered:
    """An emcee run on the synthetic orbit recovers what it was made from."""

    @pytest.fixture(scope="class")
    @classmethod
    def run(cls) -> Any:
        compiled, observed_ra, observed_dec = kit_data.synthetic(reference_astrometry)
        ra_instrument = kit_data.chain(reference_astrometry, "ra", "astrom_ra")
        dec_instrument = kit_data.chain(reference_astrometry, "dec", "astrom_dec")
        model = ReflexOrbit(
            compiled.buffers["time"].value * u.day,
            pmra=st.norm(0.0, 5.0),
            pmdec=st.norm(0.0, 5.0),
            period=st.norm(400.0, 30.0),
            phase=st.uniform(0.0, 2.0 * np.pi),
            amp_ra=st.uniform(0.0, 2.0),
            amp_dec=st.uniform(0.0, 2.0),
        )
        datasets = DatasetCollection(
            {
                "ra": Dataset(
                    observed_ra,
                    ra_instrument,
                    likelihood=Likelihood(GaussianFamily(), IndependentNoise()),
                    label="ra",
                ),
                "dec": Dataset(
                    observed_dec,
                    dec_instrument,
                    likelihood=Likelihood(GaussianFamily(), IndependentNoise()),
                    label="dec",
                ),
            }
        )
        problem = FittingProblem(model, datasets, seed=kit_data.SEED)
        return EmceeEngine(problem, walkers=24).run(steps=800, burn_in=300)

    @pytest.mark.parametrize(
        ("name", "injected"),
        [
            ("model.pmra", kit_data.ORBIT["pmra"]),
            ("model.pmdec", kit_data.ORBIT["pmdec"]),
            ("model.period", kit_data.ORBIT["period"]),
            ("model.amp_ra", kit_data.ORBIT["amp_ra"]),
            ("model.amp_dec", kit_data.ORBIT["amp_dec"]),
        ],
    )
    def test_the_injected_value_is_inside_the_central_95_per_cent(
        self, run: Any, name: str, injected: float
    ) -> None:
        draws = np.asarray(run["posterior"][name].values).reshape(-1)
        low, high = np.percentile(draws, [2.5, 97.5])
        assert low <= injected <= high
