"""The squared-amplitude step: ``|V|**2``, the observable an ``OI_VIS2`` table holds (W6.12).

The conformance row ``interferometry.rst`` §7 says a new step owes: every
registered backend's :class:`SquaredAmplitude` against a closed form that is
not ampere's (:func:`.oracles.binary_visibility`, squared), on the same array,
binary and field of view as ``test_interferometry.py``, and — where the backend
registers a realisation — the native path held to the numpy contract path on a
composed problem, which is the row that exercises the twins' ``apply_flux``.

The step is found **through the backend's own module** — the module its
``Amplitude`` lives in — rather than through a new :class:`InterferometryPieces`
field: the protocol is not this item's to change, and a backend whose
interferometry module has no ``SquaredAmplitude`` is skipped by name rather
than passed.
"""

from __future__ import annotations

import importlib
from typing import Any

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.core import (
    Dataset,
    DatasetCollection,
    FittingProblem,
    GaussianFamily,
    Instrument,
    Likelihood,
    VisibilitySet,
    negotiate,
    realise,
    registered_realisations,
)

from .oracles import binary_visibility, summed_log_abs_det
from .protocol import ConformanceBackend, InterferometryPieces, Tolerances
from .test_interferometry import (
    BINARY,
    FIELD_OF_VIEW,
    OVERSAMPLING,
    WAVELENGTH,
    image_model,
    observed_visibilities,
    pieces_or_skip,
    transformed,
    visibility_coverage,
)

#: The squared-visibility uncertainty the composed rows score with.
SIGMA_SQUARED = 0.01


def squared_amplitude(pieces: InterferometryPieces) -> type:
    """This backend's ``SquaredAmplitude``, from the module its ``Amplitude`` lives in."""
    module = importlib.import_module(pieces.amplitude.__module__)
    step = getattr(module, "SquaredAmplitude", None)
    if step is None:
        pytest.skip(
            f"{module.__name__} has no SquaredAmplitude: the twin W6.12 added beside Amplitude "
            f"is owed by this backend."
        )
    return step


def squared_observed(values: Any, sigma: Any = SIGMA_SQUARED) -> VisibilitySet:
    """A real ``VisibilitySet`` of squared visibilities on the array's coverage.

    Plain arrays make the normalised, unitless container an ``OI_VIS2`` table
    becomes; Quantities in ``Jy**2`` make the un-normalised one.
    """
    u_pts, v_pts, _, _ = visibility_coverage()
    return VisibilitySet(
        u_pts,
        v_pts,
        np.full(u_pts.size, WAVELENGTH) * u.micron,
        values,
        uncertainty=sigma * np.ones(u_pts.size),
        meta={"observable": "squared_visibility"},
    )


class TestSquaredAmplitude:
    """``|V|**2`` on every backend, against the binary's closed form."""

    def test_the_step_matches_the_squared_closed_form(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """Real values, ``|V|**2`` of the oracle, and the square of ``Amplitude``'s output."""
        pieces = pieces_or_skip(backend)
        step = squared_amplitude(pieces)
        u_pts, v_pts, _, _ = visibility_coverage()
        got = transformed(pieces, image_model(pieces, "binary"), steps=(step(),))
        values = np.asarray(backend.to_numpy(got.values))
        assert values.dtype.kind == "f"
        expected = np.abs(binary_visibility(u_pts, v_pts, **BINARY)) ** 2
        assert np.max(np.abs(values - expected)) / np.max(expected) < tolerances.cross_solver
        modulus = transformed(pieces, image_model(pieces, "binary"), steps=(pieces.amplitude(),))
        np.testing.assert_allclose(
            values, np.asarray(backend.to_numpy(modulus.values)) ** 2, rtol=1e-12, atol=0.0
        )

    def test_the_step_declares_the_backend_its_siblings_declare(
        self, backend: ConformanceBackend
    ) -> None:
        """One name per backend: the twin is this backend's piece, not the reference's."""
        pieces = pieces_or_skip(backend)
        step = squared_amplitude(pieces)
        for flag in ("BACKEND", "DEVICE"):
            assert getattr(step, flag, None) == getattr(pieces.amplitude, flag, None), flag

    def test_the_square_is_in_the_square_of_the_unit(self, backend: ConformanceBackend) -> None:
        """``Jy`` in, ``Jy**2`` out, and a Gaussian on ``Jy**2`` data aligns with it."""
        pieces = pieces_or_skip(backend)
        step = squared_amplitude(pieces)
        got = transformed(pieces, image_model(pieces, "binary"), steps=(step(),))
        assert got.unit == u.Jy**2
        assert got.uncertainty is None
        observed = squared_observed(
            np.asarray(backend.to_numpy(got.values)) * u.Jy**2, SIGMA_SQUARED * u.Jy**2
        )
        Likelihood(GaussianFamily(), backend.independent_noise()).check_alignment(got, observed)

    def test_the_normalised_square_is_unitless_and_divides_by_the_flux(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """``normalisation=F`` gives ``|V/F|**2`` with no unit — an ``OI_VIS2``'s representation."""
        pieces = pieces_or_skip(backend)
        step = squared_amplitude(pieces)
        flux = BINARY["flux"]
        got = transformed(
            pieces,
            image_model(pieces, "binary"),
            steps=(step(normalisation=flux * 1000.0 * u.mJy),),
        )
        assert got.unit is None
        u_pts, v_pts, _, _ = visibility_coverage()
        expected = np.abs(binary_visibility(u_pts, v_pts, **BINARY) / flux) ** 2
        values = np.asarray(backend.to_numpy(got.values))
        assert np.max(np.abs(values - expected)) < tolerances.cross_solver
        observed = squared_observed(values)
        Likelihood(GaussianFamily(), backend.independent_noise()).check_alignment(got, observed)

    def test_a_normalisation_that_is_not_a_flux_is_refused(
        self, backend: ConformanceBackend
    ) -> None:
        """By name, at construction, rather than as a wrong number later."""
        step = squared_amplitude(pieces_or_skip(backend))
        with pytest.raises(Exception, match="zero-spacing flux"):
            step(normalisation=3.0 * u.m)
        with pytest.raises(Exception, match="finite and positive"):
            step(normalisation=-1.0)

    def test_the_mask_is_inherited_one_to_one(self, backend: ConformanceBackend) -> None:
        """A flagged visibility stays flagged after squaring, and nothing else is."""
        pieces = pieces_or_skip(backend)
        step = squared_amplitude(pieces)
        container = observed_visibilities()
        mask = np.zeros(container.n_samples, dtype=bool)
        mask[[1, 4]] = True
        flagged = container.with_values(container.values, mask=mask)
        got = step().apply(flagged, {})
        np.testing.assert_array_equal(np.asarray(got.mask), mask)

    def test_the_native_path_matches_the_numpy_contract_path(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """The twin's ``apply_flux`` through :func:`ampere.core.realise`, at several points."""
        pieces = pieces_or_skip(backend)
        step = squared_amplitude(pieces)
        problem = self._problem(backend, pieces, step)
        at_truth = problem.log_prob({"model.separation": 12.0, "model.flux_ratio": 0.42})
        assert np.isfinite(at_truth)
        assert at_truth > problem.log_prob({"model.separation": 16.0, "model.flux_ratio": 0.42})
        if problem.backend not in registered_realisations():
            pytest.skip(f"the {problem.backend!r} backend registers no realisation")
        realised = realise(problem)
        rng = np.random.default_rng(20261005)
        for row in rng.uniform(0.05, 0.95, size=(5, problem.free_size)):
            theta = problem.prior_transform(row)
            y = problem.unconstrain(theta)
            total = float(np.asarray(backend.to_numpy(realised.log_prob_unconstrained(y))))
            assert total == pytest.approx(
                problem.log_prob_unconstrained(y), abs=tolerances.cross_backend
            )
            assert np.isfinite(total - summed_log_abs_det(problem.parameters, y))

    @staticmethod
    def _problem(backend: ConformanceBackend, pieces: InterferometryPieces, step: type) -> Any:
        """Normalised ``V**2`` of the binary at its truth, separation and flux ratio free.

        The route a fit to an ``OI_VIS2`` table takes: the model's total flux
        fixed, and the same number as the step's ``normalisation``.
        """
        instrument = Instrument(
            [
                pieces.fourier_sample.from_observed(
                    observed_visibilities(),
                    field_of_view=FIELD_OF_VIEW * u.mas,
                    oversampling=OVERSAMPLING,
                ),
                step(normalisation=BINARY["flux"]),
            ],
            channel="sky",
            label="array_vis2",
        )
        compiled = image_model(pieces, "binary").compile_for(negotiate([instrument]))
        observed = squared_observed(
            np.asarray(backend.to_numpy(instrument(compiled.evaluate()).values))
        )
        model = pieces.binary(
            compiled.buffers["x"].value * u.mas,
            compiled.buffers["y"].value * u.mas,
            channels="sky",
            **{
                **BINARY,
                "separation": st.uniform(4.0, 20.0),
                "flux_ratio": st.uniform(0.05, 0.9),
            },
        )
        datasets = DatasetCollection(
            {
                "vis2": Dataset(
                    observed,
                    instrument,
                    likelihood=Likelihood(GaussianFamily(), backend.independent_noise()),
                    label="vis2",
                )
            }
        )
        return FittingProblem(model, datasets, seed=20261005)
