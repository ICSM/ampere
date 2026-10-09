"""The squared-amplitude step: ``|V|**2``, the observable an ``OI_VIS2`` table holds (W6.12).

The conformance row ``interferometry.rst`` §7 says a new step owes: every
registered backend's :class:`SquaredAmplitude` against a closed form that is
not ampere's (:func:`.oracles.binary_visibility`, squared), on the same array,
binary and field of view as ``test_interferometry.py``, and — where the backend
registers a realisation — the native path held to the numpy contract path on a
composed problem, which is the row that exercises the twins' ``apply_flux``.

The step is :attr:`InterferometryPieces.squared_amplitude` (the field W7.4
added; W6.12 found the step through the module its ``Amplitude`` lives in,
because the protocol was not its item's to change). W7.4's rows follow the
original ones: the model form ``normalisation="model"`` against the analytic
uniform disc, against the buffer form it must equal when the buffer is the
model's total flux, on the native path, and refused on the analytic route.
"""

from __future__ import annotations

from typing import Any

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.core import (
    ClosurePhases,
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

from .oracles import binary_visibility, summed_log_abs_det, uniform_disc_visibility
from .protocol import ConformanceBackend, InterferometryPieces, Tolerances
from .test_interferometry import (
    BINARY,
    DISC,
    FIELD_OF_VIEW,
    OVERSAMPLING,
    WAVELENGTH,
    image_model,
    observed_visibilities,
    pieces_or_skip,
    transformed,
    triangle_coverage,
    visibility_coverage,
)

#: The squared-visibility uncertainty the composed rows score with.
SIGMA_SQUARED = 0.01


def squared_amplitude(pieces: InterferometryPieces) -> type:
    """This backend's ``SquaredAmplitude`` (the protocol's field since W7.4)."""
    return pieces.squared_amplitude


def squared_observed(
    values: Any, sigma: Any = SIGMA_SQUARED, wavelength: Any = None
) -> VisibilitySet:
    """A real ``VisibilitySet`` of squared visibilities on the array's coverage.

    Plain arrays make the normalised, unitless container an ``OI_VIS2`` table
    becomes; Quantities in ``Jy**2`` make the un-normalised one.
    """
    u_pts, v_pts, _, _ = visibility_coverage()
    axis = WAVELENGTH * u.micron if wavelength is None else wavelength
    return VisibilitySet(
        u_pts,
        v_pts,
        np.full(u_pts.size, axis.value) * axis.unit,
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
    def _problem(
        backend: ConformanceBackend,
        pieces: InterferometryPieces,
        step: type | None = None,
        normalisation: Any = None,
    ) -> Any:
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
                (step or pieces.squared_amplitude)(
                    normalisation=BINARY["flux"] if normalisation is None else normalisation
                ),
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
                **({"flux": st.uniform(0.5, 3.0)} if normalisation == "model" else {}),
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


# ---------------------------------------------------------------------------
# W7.4: the model form, and the observed container's own spectral unit
# ---------------------------------------------------------------------------


class TestModelNormalisation:
    """``SquaredAmplitude(normalisation="model")``: ``|V / V(0)|**2`` with the model's own flux."""

    def test_the_uniform_disc_matches_the_normalised_closed_form(
        self, backend: ConformanceBackend
    ) -> None:
        """``(2 J1(z) / z)**2``, converging as the disc's sharp edge is resolved.

        A uniform disc is not band-limited (``test_a_sharp_edged_image_converges_rather_than_
        agreeing`` in ``test_interferometry.py`` says why), so the row asserts the shape of the
        error instead of ``cross_solver``: it falls by more than a factor of twenty between
        oversampling 2 and 32 and lands inside ``5e-3``. A wrong divisor does not converge at
        all.
        """
        pieces = pieces_or_skip(backend)
        u_pts, v_pts, _, _ = visibility_coverage()
        expected = np.abs(uniform_disc_visibility(u_pts, v_pts, **DISC)) ** 2 / DISC["flux"] ** 2
        errors = []
        for oversampling in (2.0, 32.0):
            got = transformed(
                pieces,
                image_model(pieces, "disc"),
                steps=(pieces.squared_amplitude(normalisation="model"),),
                oversampling=oversampling,
            )
            assert got.unit is None
            assert got.n_samples == u_pts.size
            values = np.asarray(backend.to_numpy(got.values))
            errors.append(float(np.max(np.abs(values - expected))))
        assert errors[1] < 5.0e-3
        assert errors[1] < errors[0] / 20.0

    def test_a_fixed_buffer_and_the_model_form_agree_when_the_buffer_is_the_flux(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """The binary is band-limited here, so ``V(0)`` is its flux to ``cross_solver``."""
        pieces = pieces_or_skip(backend)
        model = image_model(pieces, "binary")
        buffer = transformed(
            pieces, model, steps=(pieces.squared_amplitude(normalisation=BINARY["flux"] * u.Jy),)
        )
        free = transformed(pieces, model, steps=(pieces.squared_amplitude(normalisation="model"),))
        np.testing.assert_allclose(
            np.asarray(backend.to_numpy(free.values)),
            np.asarray(backend.to_numpy(buffer.values)),
            rtol=0.0,
            atol=tolerances.cross_solver,
        )
        np.testing.assert_array_equal(free.u.values, buffer.u.values)

    def test_it_composes_with_a_smearing_step_in_the_chain(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """The expansions multiply: the twin of a smeared sample is smeared at zero spacing."""
        pieces = pieces_or_skip(backend)
        model = image_model(pieces, "binary")
        buffer = transformed(
            pieces,
            model,
            steps=(
                pieces.bandwidth_smearing(resolving_power=30.0, nodes=3),
                pieces.squared_amplitude(normalisation=BINARY["flux"] * u.Jy),
            ),
        )
        free = transformed(
            pieces,
            model,
            steps=(
                pieces.bandwidth_smearing(resolving_power=30.0, nodes=3),
                pieces.squared_amplitude(normalisation="model"),
            ),
        )
        np.testing.assert_allclose(
            np.asarray(backend.to_numpy(free.values)),
            np.asarray(backend.to_numpy(buffer.values)),
            rtol=0.0,
            atol=tolerances.cross_solver,
        )

    def test_the_analytic_route_is_refused_by_name(self, backend: ConformanceBackend) -> None:
        """No ``FourierSample``, nothing to expand, so the step refuses."""
        pieces = pieces_or_skip(backend)
        u_pts, v_pts, waves, _ = visibility_coverage()
        model = pieces.binary_visibilities(u_pts, v_pts, waves * u.micron, channels="vis", **BINARY)
        instrument = Instrument(
            [pieces.squared_amplitude(normalisation="model")],
            channel="vis",
            input_kind=VisibilitySet,
            label="analytic",
        )
        with pytest.raises(Exception, match="needs a FourierSample earlier in the chain"):
            instrument(model.evaluate())

    def test_a_string_that_is_not_model_is_refused(self, backend: ConformanceBackend) -> None:
        with pytest.raises(Exception, match="'model'"):
            pieces_or_skip(backend).squared_amplitude(normalisation="mode1")

    def test_the_buffer_form_does_not_expand(self, backend: ConformanceBackend) -> None:
        """Only the model form carries ``expand_uv``; the default and the buffer form do not."""
        step = pieces_or_skip(backend).squared_amplitude
        assert getattr(step(), "expand_uv", None) is None
        assert getattr(step(normalisation=1.0), "expand_uv", None) is None
        assert getattr(step(normalisation="model"), "expand_uv", None) is not None

    def test_the_native_path_matches_the_numpy_path_and_is_flux_free(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """The twin's ``apply_flux`` with the total flux a free parameter, which drops out."""
        pieces = pieces_or_skip(backend)
        problem = TestSquaredAmplitude._problem(backend, pieces, normalisation="model")
        theta = {"model.separation": 12.0, "model.flux_ratio": 0.42}
        low = problem.log_prob({**theta, "model.flux": 0.6})
        high = problem.log_prob({**theta, "model.flux": 2.9})
        elsewhere = problem.log_prob(
            {"model.separation": 16.0, "model.flux_ratio": 0.42, "model.flux": 1.0}
        )
        assert np.isfinite(low)
        assert low == pytest.approx(high, abs=1e-6)
        assert low > elsewhere
        if problem.backend not in registered_realisations():
            pytest.skip(f"the {problem.backend!r} backend registers no realisation")
        realised = realise(problem)
        rng = np.random.default_rng(20261009)
        for row in rng.uniform(0.05, 0.95, size=(5, problem.free_size)):
            y = problem.unconstrain(problem.prior_transform(row))
            total = float(np.asarray(backend.to_numpy(realised.log_prob_unconstrained(y))))
            assert total == pytest.approx(
                problem.log_prob_unconstrained(y), abs=tolerances.cross_backend
            )


class TestSpectralUnit:
    """``FourierSample.from_observed`` emits the observed container's own spectral unit (W7.4)."""

    def test_a_container_in_nanometres_composes_and_is_predicted_in_nanometres(
        self, backend: ConformanceBackend
    ) -> None:
        pieces = pieces_or_skip(backend)
        u_pts, _, _, _ = visibility_coverage()
        observed = squared_observed(np.ones(u_pts.size), wavelength=WAVELENGTH * 1e3 * u.nm)
        instrument = Instrument(
            [
                pieces.fourier_sample.from_observed(
                    observed, field_of_view=FIELD_OF_VIEW * u.mas, oversampling=OVERSAMPLING
                ),
                pieces.squared_amplitude(normalisation=BINARY["flux"]),
            ],
            channel="sky",
            label="nm_vis2",
        )
        compiled = image_model(pieces, "binary").compile_for(negotiate([instrument]))
        got = instrument(compiled.evaluate())
        assert got.spectral_axis.unit == u.nm
        np.testing.assert_allclose(
            got.spectral_axis.values, observed.spectral_axis.values, rtol=1e-12, atol=0.0
        )
        Likelihood(GaussianFamily(), backend.independent_noise()).check_alignment(got, observed)
        # And through composition: the Dataset and the problem build without refusal.
        sampled = observed.with_values(np.asarray(backend.to_numpy(got.values)))
        problem = FittingProblem(
            pieces.binary.on_field(
                FIELD_OF_VIEW * u.mas,
                4,
                channels="sky",
                **{**BINARY, "separation": st.uniform(4.0, 20.0)},
            ),
            DatasetCollection(
                {
                    "vis2": Dataset(
                        sampled,
                        instrument,
                        likelihood=Likelihood(GaussianFamily(), backend.independent_noise()),
                        label="vis2",
                    )
                }
            ),
            seed=20261009,
        )
        assert np.isfinite(problem.log_prob({"model.separation": 12.0}))

    def test_a_smeared_chain_keeps_the_unit(self, backend: ConformanceBackend) -> None:
        pieces = pieces_or_skip(backend)
        u_pts, _, _, _ = visibility_coverage()
        observed = squared_observed(np.ones(u_pts.size), wavelength=WAVELENGTH * 1e3 * u.nm)
        instrument = Instrument(
            [
                pieces.fourier_sample.from_observed(
                    observed, field_of_view=FIELD_OF_VIEW * u.mas, oversampling=OVERSAMPLING
                ),
                pieces.bandwidth_smearing(resolving_power=30.0, nodes=3),
                pieces.squared_amplitude(normalisation="model"),
            ],
            channel="sky",
        )
        compiled = image_model(pieces, "binary").compile_for(negotiate([instrument]))
        assert instrument(compiled.evaluate()).spectral_axis.unit == u.nm

    def test_closure_phases_in_nanometres_are_predicted_in_nanometres(
        self, backend: ConformanceBackend
    ) -> None:
        pieces = pieces_or_skip(backend)
        u1, v1, u2, v2, waves, _ = triangle_coverage()
        observed = ClosurePhases(u1, v1, u2, v2, waves * 1e3 * u.nm, np.zeros(u1.size) * u.rad)
        instrument = Instrument(
            [
                pieces.fourier_sample.from_observed(
                    observed, field_of_view=FIELD_OF_VIEW * u.mas, oversampling=OVERSAMPLING
                ),
                pieces.closure_phase(),
            ],
            channel="sky",
        )
        compiled = image_model(pieces, "binary").compile_for(negotiate([instrument]))
        got = instrument(compiled.evaluate())
        assert got.spectral_axis.unit == u.nm
        np.testing.assert_allclose(
            got.spectral_axis.values, observed.spectral_axis.values, rtol=1e-12, atol=0.0
        )

    def test_explicit_construction_still_emits_micron(self, backend: ConformanceBackend) -> None:
        pieces = pieces_or_skip(backend)
        u_pts, v_pts, waves, _ = visibility_coverage()
        step = pieces.fourier_sample(
            u_pts, v_pts, waves * u.micron, field_of_view=FIELD_OF_VIEW * u.mas
        )
        instrument = Instrument([step], channel="sky")
        compiled = image_model(pieces, "binary").compile_for(negotiate([instrument]))
        assert instrument(compiled.evaluate()).spectral_axis.unit == u.micron
