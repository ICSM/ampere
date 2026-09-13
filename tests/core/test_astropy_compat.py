"""W4.6's acceptance rows: the astropy interop adapter (``DEVELOPMENT_PLAN.md`` §4.7).

The battery is organised around what the item asks for, and the first class is
the one that matters most.

1. **The oracle rows** (:class:`TestItAgreesWithTheReferenceModels`). A wrapped
   ``astropy.modeling.BlackBody`` and a wrapped ``PowerLaw1D`` are composed into
   fitting problems beside :class:`ampere.backends.reference.BlackBody` and
   :class:`~ampere.backends.reference.PowerLaw` over the *same* data, and the
   two problems' ``log_prob`` are compared at the ``exact`` tolerance class
   (1e-12 relative; ``tests/conformance/protocol.py``) — digits, not
   approximation. That is only possible because the solid-angle convention is
   stated rather than fudged: astropy's ``BlackBody`` emits a surface
   brightness, ``scale=1*u.Jy/u.sr`` puts it in Jy/sr, and
   ``u.dimensionless_angles()`` states a solid angle of exactly one steradian,
   which is the convention the reference model's dimensionless ``scale``
   carries. If those two conventions ever drift apart, this class fails.

   The rows live here rather than in ``tests/conformance/`` for a structural
   reason worth recording: every test body in that directory takes the
   ``backend`` fixture and runs once per registered backend, and none may name
   a concrete backend (``tests/conformance/README.md`` §1). The adapter is
   reference-only by construction — ``BACKEND = "reference"`` is the honest
   answer for a Python callable — so a row comparing it with the reference
   backend's own models would have to name one, which the battery forbids.

2. **Parameter translation** (:class:`TestParameterTranslation`): bounds to a
   uniform, ``fixed`` to frozen, ``tied`` to an :class:`AstropyTie`, the
   ``priors=`` override, and the refusals for what cannot be translated.
3. **The tie round trip** (:class:`TestACompoundModelWithATie`), on a compound
   model, which is where astropy's ``temperature_0``/``temperature_1`` naming
   and its ``tied`` callables actually turn up.
4. **Units and the grid** (:class:`TestUnitsAndTheGrid`), including the trap
   that motivates the unit rule: ``BlackBody.input_units`` is Hz.
5. **Kinds** (:class:`TestTheKind`), the inference and the refusal.
6. **Composition** (:class:`TestItComposesLikeAnyOtherModel`): negotiation
   through ``compile_for``, the capability declaration, provenance, pickling.

``tests/core/test_astropy_engines.py`` carries the two engine rows (an emcee
fit and an SBI run), which need ``ampere.inference``.
"""

from __future__ import annotations

import ast
import copy
import pickle
import subprocess
import sys
from pathlib import Path

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st
from astropy.modeling.models import BlackBody, Const1D, Gaussian1D, Gaussian2D, PowerLaw1D

import ampere.core.astropy_compat
from ampere.backends.reference import BlackBody as ReferenceBlackBody
from ampere.backends.reference import PowerLaw as ReferencePowerLaw
from ampere.core import (
    AdaptedAstropyModel,
    AstropyTie,
    AxisRequirement,
    ChannelRequirements,
    Dataset,
    FittingProblem,
    Image,
    Instrument,
    PhotometricPoints,
    Spectrum,
    TimeSeries,
    from_astropy,
)
from ampere.core.exceptions import CompositionError, ParameterError

#: ``tests/conformance/protocol.py``'s ``Tolerances.exact``: identities that
#: hold in exact arithmetic up to representation error only.
EXACT = 1e-12

WAVELENGTH = np.geomspace(1.0, 30.0, 12)
TEMPERATURE_PRIOR = st.uniform(100.0, 9900.0)
SCALE_PRIOR = st.loguniform(0.1, 10.0)
SEED = 20260911

#: One steradian, stated. See the module docstring and ``astropy_compat.md`` §4.
ONE_STERADIAN = u.dimensionless_angles()


def observed(values: np.ndarray, sigma: float = 0.05) -> Spectrum:
    """The truth plus a fixed noise draw, with uncertainties."""
    rng = np.random.default_rng(SEED)
    noise = rng.normal(0.0, sigma * np.abs(values))
    return Spectrum(
        WAVELENGTH * u.micron,
        (values + noise) * u.Jy,
        uncertainty=(sigma * np.abs(values)) * u.Jy,
    )


def wrapped_blackbody() -> AdaptedAstropyModel:
    """``astropy``'s blackbody, in Jy, on one steradian."""
    return from_astropy(
        BlackBody(temperature=3000.0 * u.K, scale=1.0 * u.Jy / u.sr),
        grid=WAVELENGTH * u.micron,
        priors={"temperature": TEMPERATURE_PRIOR, "scale": SCALE_PRIOR},
        output_unit=u.Jy,
        equivalencies=ONE_STERADIAN,
    )


def wrapped_power_law() -> AdaptedAstropyModel:
    """``PowerLaw1D``, whose ``alpha`` is the *negative* of ampere's ``index``."""
    model = PowerLaw1D(amplitude=2.0 * u.Jy, x_0=1.0 * u.micron, alpha=1.2)
    model.x_0.fixed = True
    return from_astropy(
        model,
        grid=WAVELENGTH * u.micron,
        priors={"amplitude": SCALE_PRIOR, "alpha": st.norm(1.0, 0.5)},
        output_unit=u.Jy,
    )


# ---------------------------------------------------------------------------
# 1. The oracle rows
# ---------------------------------------------------------------------------


class TestItAgreesWithTheReferenceModels:
    """The headline acceptance criterion: same physics, same ``log_prob``, same digits."""

    def test_a_wrapped_blackbody_emits_what_the_reference_one_does(self) -> None:
        native = ReferenceBlackBody(WAVELENGTH, temperature=TEMPERATURE_PRIOR, scale=SCALE_PRIOR)
        adapted = wrapped_blackbody()
        for temperature, scale in ((3000.0, 2.0), (450.0, 0.3), (9000.0, 7.5)):
            theirs = native(temperature=temperature, scale=scale).single()
            ours = adapted(temperature=temperature, scale=scale).single()
            assert ours.unit == theirs.unit == u.Jy
            np.testing.assert_allclose(ours.values, theirs.values, rtol=EXACT, atol=0.0)

    def test_a_wrapped_blackbody_scores_the_same_log_prob(self) -> None:
        """The criterion as written: ``log_prob`` on the *same* data."""
        native = ReferenceBlackBody(WAVELENGTH, temperature=TEMPERATURE_PRIOR, scale=SCALE_PRIOR)
        data = observed(native(temperature=3000.0, scale=2.0).single().values)
        their_problem = FittingProblem(native, [Dataset(data)], seed=SEED)
        our_problem = FittingProblem(wrapped_blackbody(), [Dataset(data)], seed=SEED)
        assert our_problem.parameters.free_names == their_problem.parameters.free_names
        for temperature, scale in ((3000.0, 2.0), (2500.0, 1.7), (5000.0, 0.4)):
            values = {"model.temperature": temperature, "model.scale": scale}
            assert np.isfinite(their_problem.log_prob(values))
            np.testing.assert_allclose(
                our_problem.log_prob(values), their_problem.log_prob(values), rtol=EXACT, atol=0.0
            )

    def test_a_wrapped_power_law_scores_the_same_log_prob(self) -> None:
        """``PowerLaw1D(alpha)`` is ``PowerLaw(index=-alpha)``; the likelihood cannot tell."""
        native = ReferencePowerLaw(
            WAVELENGTH, norm=SCALE_PRIOR, index=st.norm(-1.0, 0.5), reference_wavelength=1.0
        )
        data = observed(native(norm=2.0, index=-1.2).single().values)
        their_problem = FittingProblem(native, [Dataset(data)], seed=SEED)
        our_problem = FittingProblem(wrapped_power_law(), [Dataset(data)], seed=SEED)
        for norm, index in ((2.0, -1.2), (0.5, -0.4), (7.0, -2.0)):
            theirs = their_problem.log_prob({"model.norm": norm, "model.index": index})
            ours = our_problem.log_prob({"model.amplitude": norm, "model.alpha": -index})
            # The priors differ by the sign flip, so compare the likelihoods,
            # which is the half the two models are supposed to share.
            their_likelihood = their_problem.evaluate(
                {"model.norm": norm, "model.index": index}
            ).log_likelihood
            our_likelihood = our_problem.evaluate(
                {"model.amplitude": norm, "model.alpha": -index}
            ).log_likelihood
            assert np.isfinite(theirs) and np.isfinite(ours)
            np.testing.assert_allclose(our_likelihood, their_likelihood, rtol=EXACT, atol=0.0)

    def test_the_solid_angle_convention_is_one_steradian_and_nothing_else(self) -> None:
        """Drop the equivalency and the adapter refuses rather than guessing."""
        with pytest.raises(CompositionError, match="solid angle"):
            from_astropy(
                BlackBody(temperature=3000.0 * u.K, scale=1.0 * u.Jy / u.sr),
                grid=WAVELENGTH * u.micron,
                priors={"temperature": TEMPERATURE_PRIOR, "scale": SCALE_PRIOR},
                output_unit=u.Jy,
            )


# ---------------------------------------------------------------------------
# 2. Parameter translation
# ---------------------------------------------------------------------------


class TestParameterTranslation:
    """§4.7's table: bounds, ``fixed``, ``tied``, and the ``priors=`` override."""

    def test_finite_bounds_become_a_uniform_prior_in_the_parameters_own_unit(self) -> None:
        model = Gaussian1D(1.0, 2.0, 0.5)
        model.amplitude.bounds = (0.5, 4.5)
        model.mean.bounds = (0.0, 10.0)
        model.stddev.fixed = True
        adapted = from_astropy(model, grid=WAVELENGTH * u.micron)
        amplitude = adapted.parameters["amplitude"]
        assert amplitude.is_free
        assert amplitude.prior.ppf(0.0) == pytest.approx(0.5)
        assert amplitude.prior.ppf(1.0) == pytest.approx(4.5)
        assert "0.5" in amplitude.description

    def test_a_unit_ful_parameters_bounds_are_read_in_that_unit(self) -> None:
        model = PowerLaw1D(amplitude=2.0 * u.Jy, x_0=1.0 * u.micron, alpha=1.0)
        model.amplitude.bounds = (0.5, 5.0)
        model.x_0.fixed = True
        model.alpha.bounds = (0.0, 3.0)
        adapted = from_astropy(model, grid=WAVELENGTH * u.micron, output_unit=u.Jy)
        assert adapted.parameters["amplitude"].unit == u.Jy
        assert adapted.parameters["amplitude"].prior.ppf(1.0) == pytest.approx(5.0)

    def test_fixed_becomes_frozen_and_costs_no_sampler_dimension(self) -> None:
        adapted = wrapped_power_law()
        assert adapted.parameters.fixed_names == ("x_0",)
        assert adapted.parameters["x_0"].value == pytest.approx(1.0)
        assert adapted.parameters.free_size == 2

    def test_a_priors_entry_overrides_the_translated_prior(self) -> None:
        model = Gaussian1D(1.0, 2.0, 0.5)
        for name in model.param_names:
            getattr(model, name).bounds = (0.01, 10.0)
        adapted = from_astropy(
            model, grid=WAVELENGTH * u.micron, priors={"mean": st.norm(2.0, 0.1), "stddev": 0.5}
        )
        assert adapted.parameters["mean"].prior.mean() == pytest.approx(2.0)
        assert adapted.parameters["stddev"].is_fixed

    def test_a_free_unbounded_parameter_is_refused_by_name(self) -> None:
        """``BlackBody.temperature`` ships ``bounds=(0, None)``: half-open is not a prior."""
        with pytest.raises(ParameterError) as raised:
            from_astropy(
                BlackBody(temperature=3000.0 * u.K), grid=WAVELENGTH * u.micron, output_unit=None
            )
        message = str(raised.value)
        assert "'temperature' is free and unbounded" in message
        assert "priors=" in message
        assert "parameter.bounds = (lo, hi)" in message

    def test_a_prior_for_a_parameter_the_model_does_not_have_is_refused(self) -> None:
        with pytest.raises(ParameterError, match="does not declare"):
            from_astropy(
                BlackBody(temperature=3000.0 * u.K, scale=1.0 * u.Jy / u.sr),
                grid=WAVELENGTH * u.micron,
                priors={"temperature": TEMPERATURE_PRIOR, "scale": SCALE_PRIOR, "radius": 1.0},
                output_unit=u.Jy,
                equivalencies=ONE_STERADIAN,
            )

    def test_a_prior_for_a_tied_parameter_is_refused(self) -> None:
        model = Gaussian1D(1.0, 2.0, 0.5) + Const1D(0.3)
        model.mean_0.tied = _twice_the_amplitude
        for name in ("amplitude_0", "stddev_0", "amplitude_1"):
            getattr(model, name).bounds = (0.01, 10.0)
        with pytest.raises(ParameterError, match="astropy declares it tied"):
            from_astropy(
                model,
                grid=WAVELENGTH * u.micron,
                kind=Spectrum,
                priors={"mean_0": st.norm(0.0, 1.0)},
            )

    def test_a_nonsense_prior_is_refused(self) -> None:
        model = Gaussian1D(1.0, 2.0, 0.5)
        for name in model.param_names:
            getattr(model, name).bounds = (0.01, 10.0)
        with pytest.raises(ParameterError, match=r"frozen scipy\.stats distribution"):
            from_astropy(model, grid=WAVELENGTH * u.micron, priors={"mean": "flat"})


# ---------------------------------------------------------------------------
# 3. The tie round trip
# ---------------------------------------------------------------------------


def _twice_the_amplitude(model: object) -> float:
    """A module-level tie, so the wrapped model stays picklable."""
    return float(model.amplitude_0.value) * 2.0  # type: ignore[attr-defined]


def tied_compound() -> tuple[object, AdaptedAstropyModel]:
    """``Gaussian1D + Const1D``, with ``mean_0`` tied to twice ``amplitude_0``."""
    model = Gaussian1D(3.0, 6.0, 1.5) + Const1D(0.3)
    model.mean_0.tied = _twice_the_amplitude
    for name in ("amplitude_0", "stddev_0", "amplitude_1"):
        getattr(model, name).bounds = (0.01, 12.0)
    return model, from_astropy(model, grid=WAVELENGTH * u.micron, kind=Spectrum, output_unit=u.Jy)


class TestACompoundModelWithATie:
    """A compound model round-trips its tie: recorded, applied, and reported."""

    def test_the_tie_is_recorded_rather_than_declared_as_a_parameter(self) -> None:
        _, adapted = tied_compound()
        assert set(adapted.ties) == {"mean_0"}
        assert isinstance(adapted.ties["mean_0"], AstropyTie)
        assert adapted.ties["mean_0"].function is _twice_the_amplitude
        assert "mean_0" not in adapted.parameters
        assert adapted.parameters.free_names == ("amplitude_0", "stddev_0", "amplitude_1")
        assert adapted.parameters.free_size == 3

    def test_astropys_compound_names_are_kept_exactly(self) -> None:
        """``temperature_0``/``temperature_1`` is astropy's convention, and it survives."""
        model = BlackBody(temperature=3000.0 * u.K, scale=1.0 * u.Jy / u.sr) + BlackBody(
            temperature=300.0 * u.K, scale=1.0 * u.Jy / u.sr
        )
        adapted = from_astropy(
            model,
            grid=WAVELENGTH * u.micron,
            kind=Spectrum,
            priors={
                "temperature_0": TEMPERATURE_PRIOR,
                "scale_0": SCALE_PRIOR,
                "temperature_1": TEMPERATURE_PRIOR,
                "scale_1": SCALE_PRIOR,
            },
            output_unit=u.Jy,
            equivalencies=ONE_STERADIAN,
        )
        assert adapted.parameters.names == (
            "temperature_0",
            "scale_0",
            "temperature_1",
            "scale_1",
        )

    def test_the_tie_is_applied_on_every_evaluation(self) -> None:
        _, adapted = tied_compound()
        for amplitude in (1.0, 4.0, 2.5):
            result = adapted(amplitude_0=amplitude, stddev_0=1.5, amplitude_1=0.3)
            expected = (Gaussian1D(amplitude, 2.0 * amplitude, 1.5) + Const1D(0.3))(WAVELENGTH)
            assert adapted.astropy_model.mean_0.value == pytest.approx(2.0 * amplitude)
            np.testing.assert_allclose(result.single().values, expected, rtol=EXACT, atol=0.0)

    def test_the_callers_own_astropy_model_is_never_mutated(self) -> None:
        original, adapted = tied_compound()
        before = np.array(original.parameters, copy=True)
        adapted(amplitude_0=5.0, stddev_0=1.5, amplitude_1=0.3)
        np.testing.assert_array_equal(np.asarray(original.parameters), before)
        assert adapted.astropy_model is not original

    def test_the_tie_survives_provenance(self) -> None:
        _, adapted = tied_compound()
        description = adapted.describe()
        assert description["astropy_ties"] == [
            {"name": "mean_0", "function": "_twice_the_amplitude"}
        ]
        assert description["astropy_components"] == ["Gaussian1D", "Const1D"]


# ---------------------------------------------------------------------------
# 4. Units and the grid
# ---------------------------------------------------------------------------


class TestUnitsAndTheGrid:
    """The trap the unit rule exists for, and the conversions that are hoisted."""

    def test_a_bare_grid_is_refused_for_a_model_with_declared_input_units(self) -> None:
        """``BlackBody.input_units`` is Hz; bare micron numbers would be read as Hz."""
        with pytest.raises(CompositionError, match="must carry a unit"):
            from_astropy(
                BlackBody(temperature=3000.0 * u.K, scale=1.0 * u.Jy / u.sr),
                grid=WAVELENGTH,
                kind=Spectrum,
                priors={"temperature": TEMPERATURE_PRIOR, "scale": SCALE_PRIOR},
                output_unit=u.Jy,
                equivalencies=ONE_STERADIAN,
            )

    def test_a_wavelength_grid_reaches_a_frequency_model_correctly(self) -> None:
        """The proof the Quantity route works: agreement with the frequency call."""
        adapted = wrapped_blackbody()
        direct = BlackBody(temperature=3000.0 * u.K, scale=1.0 * u.Jy / u.sr)
        frequency = (WAVELENGTH * u.micron).to(u.Hz, equivalencies=u.spectral())
        expected = 2.0 * direct(frequency).to_value(u.Jy, equivalencies=ONE_STERADIAN)
        np.testing.assert_allclose(
            adapted(temperature=3000.0, scale=2.0).single().values,
            expected,
            rtol=EXACT,
            atol=0.0,
        )

    def test_a_unitless_model_has_its_output_unit_declared_not_converted(self) -> None:
        model = Gaussian1D(1.0, 2.0, 0.5)
        for name in model.param_names:
            getattr(model, name).bounds = (0.01, 10.0)
        adapted = from_astropy(model, grid=WAVELENGTH * u.micron, output_unit=u.Jy)
        result = adapted(amplitude=1.0, mean=2.0, stddev=0.5).single()
        assert result.unit == u.Jy
        np.testing.assert_allclose(result.values, model(WAVELENGTH), rtol=EXACT, atol=0.0)

    def test_spectral_density_is_enabled_for_a_spectral_axis(self) -> None:
        """The one equivalency the adapter turns on unasked, because it is exact."""
        model = PowerLaw1D(
            amplitude=2.0 * u.erg / (u.cm**2 * u.s * u.AA), x_0=1.0 * u.micron, alpha=1.0
        )
        model.x_0.fixed = True
        model.alpha.bounds = (0.0, 3.0)
        adapted = from_astropy(
            model,
            grid=WAVELENGTH * u.micron,
            priors={"amplitude": SCALE_PRIOR},
            output_unit=u.Jy,
        )
        raw = model(WAVELENGTH * u.micron)
        expected = raw.to_value(u.Jy, equivalencies=u.spectral_density(WAVELENGTH * u.micron))
        np.testing.assert_allclose(
            adapted(amplitude=2.0, alpha=1.0).single().values, expected, rtol=EXACT, atol=0.0
        )

    def test_an_impossible_conversion_names_the_two_units_and_the_remedy(self) -> None:
        with pytest.raises(CompositionError) as raised:
            from_astropy(
                BlackBody(temperature=3000.0 * u.K, scale=1.0 * u.Jy / u.sr),
                grid=WAVELENGTH * u.micron,
                priors={"temperature": TEMPERATURE_PRIOR, "scale": SCALE_PRIOR},
                output_unit=u.Jy,
            )
        message = str(raised.value)
        assert "Jy / sr" in message
        assert "dimensionless_angles" in message
        assert "will not invent one" in message

    def test_the_unit_factor_is_computed_once_not_per_evaluation(self) -> None:
        """``DEVELOPMENT_PLAN.md`` §7's units trap, asserted rather than described."""
        adapted = wrapped_blackbody()
        factor = adapted._factor
        adapted(temperature=1000.0, scale=1.0)
        adapted(temperature=8000.0, scale=3.0)
        assert adapted._factor is factor

    def test_successive_results_share_their_axes_by_identity(self) -> None:
        adapted = wrapped_blackbody()
        first = adapted(temperature=1000.0, scale=1.0).single()
        second = adapted(temperature=8000.0, scale=3.0).single()
        assert first.spectral_axis is second.spectral_axis


# ---------------------------------------------------------------------------
# 5. The kind
# ---------------------------------------------------------------------------


class TestTheKind:
    """Declared by the caller; inferred only where the axis unit settles it."""

    def test_an_undeclared_kind_with_no_grid_is_refused_word_for_word(self) -> None:
        """The acceptance criterion's refusal, checked sentence by sentence."""
        with pytest.raises(CompositionError) as raised:
            from_astropy(Gaussian1D(1.0, 2.0, 0.5))
        assert str(raised.value) == (
            "from_astropy() cannot infer the ModelResult kind for Gaussian1D: no grid was "
            "given, so there are no axis units to read, and the negotiated grid is not known "
            "until compile_for(). Declare it — from_astropy(model, kind=Spectrum, ...) — "
            "choosing from Spectrum (1 axis), TimeSeries (1 axis), Image (2 axes). The kind is "
            "never guessed from the model's class name: §4.7 makes it the caller's declaration, "
            "and the wrong kind is a fit that runs and is silently wrong."
        )

    def test_a_spectral_axis_infers_a_spectrum(self) -> None:
        assert wrapped_power_law().kind is Spectrum

    def test_a_time_axis_infers_a_time_series(self) -> None:
        model = Gaussian1D(1.0, 2.0, 0.5)
        for name in model.param_names:
            getattr(model, name).bounds = (0.01, 10.0)
        adapted = from_astropy(model, grid=np.linspace(0.1, 5.0, 6) * u.s)
        assert adapted.kind is TimeSeries
        assert adapted(amplitude=1.0, mean=2.0, stddev=0.5).single().time.size == 6

    def test_two_inputs_infer_an_image(self) -> None:
        model = Gaussian2D(1.0, 0.0, 0.0, 1.0, 1.0)
        for name in model.param_names:
            getattr(model, name).bounds = (-5.0, 5.0)
        adapted = from_astropy(
            model,
            grid={"x": np.linspace(-2.0, 2.0, 5) * u.deg, "y": np.linspace(-1.0, 1.0, 3) * u.deg},
        )
        assert adapted.kind is Image
        assert adapted().single().values.shape == (5, 3)

    def test_an_unrecognisable_axis_unit_is_refused_by_name(self) -> None:
        model = Gaussian1D(1.0, 2.0, 0.5)
        for name in model.param_names:
            getattr(model, name).bounds = (0.01, 10.0)
        with pytest.raises(CompositionError, match="physical type"):
            from_astropy(model, grid=np.linspace(1.0, 5.0, 4) * u.kg)

    def test_a_kind_the_adapter_does_not_build_is_refused_with_the_list(self) -> None:
        with pytest.raises(CompositionError, match="PhotometricPoints needs a filter list"):
            from_astropy(
                Gaussian1D(1.0, 2.0, 0.5),
                grid=WAVELENGTH * u.micron,
                kind=PhotometricPoints,
            )

    def test_a_kind_that_disagrees_with_the_models_arity_is_refused(self) -> None:
        with pytest.raises(CompositionError, match="n_inputs=2"):
            from_astropy(
                Gaussian2D(1.0, 0.0, 0.0, 1.0, 1.0), grid=WAVELENGTH * u.micron, kind=Spectrum
            )

    def test_a_multi_output_model_is_refused(self) -> None:
        model = Gaussian1D(1.0, 2.0, 0.5) & Const1D(0.3)
        with pytest.raises(CompositionError, match="n_outputs=2"):
            from_astropy(model, grid=WAVELENGTH * u.micron, kind=Spectrum)

    def test_something_that_is_not_an_astropy_model_is_refused(self) -> None:
        with pytest.raises(CompositionError, match=r"takes an astropy\.modeling\.Model"):
            from_astropy(lambda x: x, grid=WAVELENGTH * u.micron, kind=Spectrum)


# ---------------------------------------------------------------------------
# 6. Composition
# ---------------------------------------------------------------------------


class TestItComposesLikeAnyOtherModel:
    """An adapted model is a ``Model``: negotiation, capabilities, provenance, pickling."""

    def test_the_capability_declaration_is_the_honest_one(self) -> None:
        adapted = wrapped_blackbody()
        assert (adapted.DIFFERENTIABLE, adapted.BATCHABLE, adapted.BACKEND) == (
            False,
            False,
            "reference",
        )
        problem = FittingProblem(adapted, [Dataset(observed(np.ones(WAVELENGTH.size)))], seed=SEED)
        assert problem.capabilities.differentiable is False
        assert problem.backend == "reference"
        assert problem.foreign_parts == ()

    def test_compile_for_adopts_the_negotiated_grid(self) -> None:
        """The template pattern, on a model built without a grid of its own."""
        adapted = from_astropy(
            BlackBody(temperature=3000.0 * u.K, scale=1.0 * u.Jy / u.sr),
            kind=Spectrum,
            priors={"temperature": TEMPERATURE_PRIOR, "scale": SCALE_PRIOR},
            output_unit=u.Jy,
            equivalencies=ONE_STERADIAN,
        )
        with pytest.raises(CompositionError, match="no grid to evaluate on"):
            adapted(temperature=3000.0, scale=1.0)
        asked = AxisRequirement(
            "spectral_axis", unit=u.micron, intervals=[(2.0, 10.0)], max_step=2.0
        )
        adapted.compile_for(
            {"default": ChannelRequirements("default", Spectrum, {"spectral_axis": asked})}
        )
        result = adapted(temperature=3000.0, scale=2.0).single()
        np.testing.assert_allclose(
            result.spectral_axis.values, asked.coordinates().to_value(u.micron), rtol=EXACT
        )
        assert result.unit == u.Jy

    def test_a_partially_negotiated_image_channel_is_refused_loudly(self) -> None:
        """``transformations.md`` §15 Q3's loud option: engaged but unable to honour."""
        model = Gaussian2D(1.0, 0.0, 0.0, 1.0, 1.0)
        for name in model.param_names:
            getattr(model, name).bounds = (-5.0, 5.0)
        adapted = from_astropy(model, kind=Image)
        asked = AxisRequirement("x", unit=u.deg, intervals=[(-2.0, 2.0)], max_step=1.0)
        with pytest.raises(CompositionError, match="nothing supplies axis 'y'"):
            adapted.compile_for({"default": ChannelRequirements("default", Image, {"x": asked})})

    def test_it_emits_under_the_channel_it_was_given(self) -> None:
        model = PowerLaw1D(amplitude=2.0 * u.Jy, x_0=1.0 * u.micron, alpha=1.2)
        model.x_0.fixed = True
        model.alpha.bounds = (0.0, 3.0)
        adapted = from_astropy(
            model,
            grid=WAVELENGTH * u.micron,
            priors={"amplitude": SCALE_PRIOR},
            output_unit=u.Jy,
            channel="continuum",
        )
        result = adapted(amplitude=2.0, alpha=1.2)
        assert list(result) == ["continuum"]
        Instrument(channel="continuum").bind(result)

    def test_the_grid_is_declared_as_a_buffer(self) -> None:
        adapted = wrapped_blackbody()
        assert adapted.buffers.names == ("spectral_axis",)
        assert adapted.buffers["spectral_axis"].unit == u.micron

    def test_describe_names_what_changes_the_computation(self) -> None:
        description = wrapped_blackbody().describe()
        assert description["adapter"] == "astropy"
        assert description["astropy_class"] == "BlackBody"
        assert description["kind"] == "Spectrum"
        assert description["output_unit"] == "Jy"
        assert description["raw_unit"] == "Jy / sr"

    def test_an_adapted_model_pickles_for_the_process_pool(self) -> None:
        """``simulate_many``'s ``ProcessExecutor`` sends the whole problem to a worker."""
        adapted = wrapped_blackbody()
        revived = pickle.loads(pickle.dumps(adapted))
        np.testing.assert_allclose(
            revived(temperature=3000.0, scale=2.0).single().values,
            adapted(temperature=3000.0, scale=2.0).single().values,
            rtol=EXACT,
            atol=0.0,
        )

    def test_a_compound_model_with_a_module_level_tie_pickles_too(self) -> None:
        _, adapted = tied_compound()
        revived = pickle.loads(pickle.dumps(adapted))
        assert set(revived.ties) == {"mean_0"}
        np.testing.assert_allclose(
            revived(amplitude_0=4.0, stddev_0=1.5, amplitude_1=0.3).single().values,
            adapted(amplitude_0=4.0, stddev_0=1.5, amplitude_1=0.3).single().values,
            rtol=EXACT,
            atol=0.0,
        )

    def test_deep_copying_the_wrapper_keeps_the_two_independent(self) -> None:
        adapted = wrapped_blackbody()
        clone = copy.deepcopy(adapted)
        clone(temperature=9000.0, scale=1.0)
        adapted(temperature=1000.0, scale=1.0)
        assert clone.astropy_model.temperature.value == pytest.approx(9000.0)
        assert adapted.astropy_model.temperature.value == pytest.approx(1000.0)


class TestTheImportCost:
    """``astropy.modeling`` is imported on use, not at module scope.

    Checked by **parsing the module**, the way ``tests/inference/test_engines.py``
    checks ``ampere.inference``'s import graph, rather than by grepping the
    prose. A runtime probe — ``'astropy.modeling' in sys.modules`` after
    ``import ampere.core`` — cannot be written today and the reason is worth
    recording: ``ampere/__init__.py`` imports the frozen legacy
    ``ampere.models``, and ``ampere/models/starScreen.py`` imports
    ``astropy.modeling`` at module scope, so the submodule is already present
    before ``ampere.core`` is reached. That is a legacy cost this item does not
    touch; what it must not do is *add* one, and this row is what holds it to
    that.
    """

    def test_astropy_modeling_is_imported_inside_a_function_not_at_module_scope(self) -> None:
        source = Path(ampere.core.astropy_compat.__file__).read_text(encoding="utf-8")
        tree = ast.parse(source)
        top_level = {
            alias.name.split(".")[0]
            for node in tree.body
            if isinstance(node, (ast.Import, ast.ImportFrom))
            for alias in getattr(node, "names", [])
        }
        top_level |= {
            node.module.split(".")[0]
            for node in tree.body
            if isinstance(node, ast.ImportFrom) and node.module
        }
        assert "astropy" in top_level  # astropy.units, which is cheap and always needed
        modeling_imports = [
            node
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom)
            and (node.module or "").startswith("astropy.modeling")
        ]
        assert modeling_imports, "the adapter must import astropy.modeling somewhere"
        assert all(node not in tree.body for node in modeling_imports), (
            "astropy.modeling must be imported inside a function, so that importing "
            "ampere.core does not pay for it"
        )

    def test_importing_the_module_alone_does_not_import_astropy_modeling(self) -> None:
        """The same claim at run time, with the legacy import stubbed out of the way.

        The frozen legacy subpackages are replaced by empty stand-ins before
        ``ampere.core`` is imported, so what is measured is this package's own
        import cost rather than ``ampere/models/starScreen.py``'s.
        """
        probe = subprocess.run(
            [
                sys.executable,
                "-c",
                (
                    "import sys, types\n"
                    "for name in ('ampere.models', 'ampere.data', 'ampere.infer', "
                    "'ampere.utils'):\n"
                    "    sys.modules[name] = types.ModuleType(name)\n"
                    "import ampere.core\n"
                    "print('astropy.modeling' in sys.modules)\n"
                    "print(hasattr(ampere.core, 'from_astropy'))\n"
                ),
            ],
            capture_output=True,
            text=True,
            check=True,
        )
        assert probe.stdout.split() == ["False", "True"]
