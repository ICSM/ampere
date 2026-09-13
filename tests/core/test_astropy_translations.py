"""W4.7's backend-neutral physics (``ampere.core.astropy_translations``).

The curated table's mathematics — six leaf formulas plus the composition rule
that walks astropy's own expression tree — is written once, against
:class:`~ampere.core.astropy_translations.TranslationOps`, so both
``ampere.backends.torch.astropy`` and ``ampere.backends.jax.astropy`` build
their tables from it (the module's own docstring: "one table serves both
backends"). This file holds that shared physics to the astropy formulas
directly, in plain numpy, without either extra installed — the backend rows
(conformance against the black-box adapter, NUTS, the two refusals) are
``tests/core/test_astropy_translations_backends.py``, parametrised over
whichever of ``torch``/``jax`` this environment has.

Also covers the two small extensions ``ampere.core.astropy_compat`` grew for
this item: :func:`~ampere.core.astropy_compat.translate_astropy_parameters`
(exposed so a native translation and the black-box adapter cannot drift on
what a bound, a fixed value or a tie means) and the ``grids``/``factor``/
``template`` read-only views a native translation reuses rather than
re-deriving.
"""

from __future__ import annotations

import astropy.units as u
import numpy as np
import pytest
from astropy.modeling.models import (
    BlackBody,
    BrokenPowerLaw1D,
    Const1D,
    Gaussian1D,
    Polynomial1D,
    PowerLaw1D,
)

from ampere.core import Spectrum, from_astropy, translate_astropy_parameters
from ampere.core.astropy_translations import (
    LEAF_BUILDERS,
    leaf_grid_values,
    plan_translation,
    tied_parameter_refusal,
    unsupported_operator_refusal,
)
from ampere.core.exceptions import CapabilityError

GRID = np.array([0.5, 1.0, 2.0, 5.0])


class _NumpyOps:
    """A minimal, independent :class:`TranslationOps`, for testing the formulas alone."""

    def exp(self, array):
        return np.exp(array)

    def where(self, condition, if_true, if_false):
        return np.where(condition, if_true, if_false)

    def blackbody(self, wavelength, temperature):
        import astropy.constants as const

        frequency = (const.c / (wavelength * u.micron)).to(u.Hz)
        radiance = (
            2.0
            * const.h
            * frequency**3
            / const.c**2
            / np.expm1((const.h * frequency / (const.k_B * temperature * u.K)).to_value(""))
        )
        return radiance.to_value(u.Jy)


OPS = _NumpyOps()


def _evaluate(formula, leaf, grid):
    grid_values = leaf_grid_values(formula, grid * u.micron)
    params = {name: getattr(leaf, name).value for name in formula.param_names}
    return formula.compute(OPS, grid_values, **params)


class TestTheCuratedFormulas:
    """Each of the six leaves' formula against astropy's own ``evaluate``."""

    def test_power_law(self) -> None:
        leaf = PowerLaw1D(amplitude=2.0, x_0=1.0, alpha=1.5)
        formula = LEAF_BUILDERS[PowerLaw1D](leaf)
        assert _evaluate(formula, leaf, GRID) == pytest.approx(leaf(GRID))

    def test_broken_power_law(self) -> None:
        leaf = BrokenPowerLaw1D(amplitude=2.0, x_break=1.5, alpha_1=1.0, alpha_2=2.0)
        formula = LEAF_BUILDERS[BrokenPowerLaw1D](leaf)
        assert _evaluate(formula, leaf, GRID) == pytest.approx(leaf(GRID))

    def test_gaussian(self) -> None:
        leaf = Gaussian1D(amplitude=3.0, mean=2.0, stddev=0.7)
        formula = LEAF_BUILDERS[Gaussian1D](leaf)
        assert _evaluate(formula, leaf, GRID) == pytest.approx(leaf(GRID))

    def test_const(self) -> None:
        leaf = Const1D(amplitude=0.4)
        formula = LEAF_BUILDERS[Const1D](leaf)
        assert _evaluate(formula, leaf, GRID) == pytest.approx(leaf(GRID))

    def test_polynomial(self) -> None:
        leaf = Polynomial1D(degree=2, c0=1.0, c1=2.0, c2=3.0)
        formula = LEAF_BUILDERS[Polynomial1D](leaf)
        assert _evaluate(formula, leaf, GRID) == pytest.approx(leaf(GRID))

    def test_polynomial_with_a_rescaled_domain(self) -> None:
        """``domain``/``window`` rescaling — Horner's own affine map, matched."""
        leaf = Polynomial1D(degree=1, c0=0.0, c1=1.0, domain=(0.0, 10.0), window=(-1.0, 1.0))
        formula = LEAF_BUILDERS[Polynomial1D](leaf)
        assert _evaluate(formula, leaf, GRID) == pytest.approx(leaf(GRID))
        assert formula.extra == {"domain": (0.0, 10.0), "window": (-1.0, 1.0)}

    def test_a_leaf_with_no_domain_has_no_extra_configuration(self) -> None:
        leaf = Polynomial1D(degree=1, c0=0.0, c1=1.0)
        formula = LEAF_BUILDERS[Polynomial1D](leaf)
        assert formula.extra == {}

    def test_blackbody_against_the_reference_convention(self) -> None:
        """``scale * B_nu(T)``, matched to astropy's own raw answer, digit for digit.

        The leaf's ``compute()`` reproduces astropy's own raw numeric value in
        astropy's own raw unit (here, dimensionless ``scale`` puts that in a
        CGS-style unit) — the invariant
        ``test_blackbody_reproduces_astropys_own_raw_unit_regardless_of_scales_unit``
        states explicitly; this row is the same claim at the convention
        ``docs/design/contracts/astropy_compat.md`` §4.3 fixes.
        """
        leaf = BlackBody(temperature=3000.0 * u.K, scale=1.0)
        formula = LEAF_BUILDERS[BlackBody](leaf)
        got = _evaluate(formula, leaf, GRID)
        expected = leaf(GRID * u.micron).value
        assert got == pytest.approx(expected, rel=1e-9)

    def test_blackbody_reproduces_astropys_own_raw_unit_regardless_of_scales_unit(self) -> None:
        """The invariant every leaf's ``compute()`` must hold for a whole-model factor to compose.

        Dimensionless ``scale`` puts astropy's raw output in a CGS-style unit
        (``erg / (Hz s sr cm2)``); ``scale`` declared as ``Jy / sr`` puts it in
        ``Jy / sr`` directly. Both must come out of this leaf's ``compute()``
        exactly as astropy's own raw evaluate would, in astropy's own raw
        unit, or the backend hook's single output-unit factor (computed once,
        against the *real* astropy model) would not apply correctly.
        """
        dimensionless = BlackBody(temperature=3000.0 * u.K, scale=1.0)
        with_unit = BlackBody(temperature=3000.0 * u.K, scale=1.0 * u.Jy / u.sr)
        for leaf in (dimensionless, with_unit):
            formula = LEAF_BUILDERS[BlackBody](leaf)
            got = _evaluate(formula, leaf, GRID)
            expected = leaf(GRID * u.micron).value
            assert got == pytest.approx(expected, rel=1e-9)


class TestGridUnits:
    """A leaf built with ``Quantity`` parameters wants its own declared grid unit."""

    def test_a_bare_number_leaf_takes_the_grids_own_unit(self) -> None:
        leaf = PowerLaw1D(amplitude=2.0, x_0=1.0, alpha=1.5)
        formula = LEAF_BUILDERS[PowerLaw1D](leaf)
        assert formula.grid_unit is None
        values = leaf_grid_values(formula, GRID * u.micron)
        assert values == pytest.approx(GRID)

    def test_a_quantity_leaf_converts_the_grid_once(self) -> None:
        leaf = PowerLaw1D(amplitude=2.0 * u.Jy, x_0=100.0 * u.nm, alpha=1.5)
        formula = LEAF_BUILDERS[PowerLaw1D](leaf)
        assert formula.grid_unit == u.nm
        values = leaf_grid_values(formula, GRID * u.micron)
        assert values == pytest.approx(GRID * 1000.0)

    def test_blackbody_always_wants_micron_with_spectral_in_force(self) -> None:
        leaf = BlackBody(temperature=3000.0 * u.K, scale=1.0)
        formula = LEAF_BUILDERS[BlackBody](leaf)
        assert formula.grid_unit == u.micron
        values = leaf_grid_values(formula, GRID * u.nm)
        assert values == pytest.approx(GRID / 1000.0)


class TestComposition:
    """:func:`plan_translation` follows astropy's own expression tree."""

    def test_a_single_leaf_has_no_suffix(self) -> None:
        model = PowerLaw1D(amplitude=2.0, x_0=1.0, alpha=1.0)
        leaves, _ = plan_translation(model, "torch")
        assert [plan.suffix for plan in leaves] == [""]
        assert leaves[0].context_names == ("amplitude", "x_0", "alpha")

    def test_a_compound_sum_numbers_leaves_left_to_right(self) -> None:
        model = Gaussian1D(1.0, 2.0, 0.5) + Const1D(0.3)
        leaves, combine = plan_translation(model, "torch")
        assert [plan.suffix for plan in leaves] == ["_0", "_1"]
        assert leaves[0].context_names == ("amplitude_0", "mean_0", "stddev_0")
        assert leaves[1].context_names == ("amplitude_1",)
        combined = combine([np.array([1.0]), np.array([2.0])])
        assert combined == pytest.approx([3.0])

    def test_a_nested_compound_numbers_all_the_way_down(self) -> None:
        model = (Gaussian1D(1.0, 2.0, 0.5) + Const1D(0.3)) * PowerLaw1D(1.0, 1.0, 1.0)
        leaves, combine = plan_translation(model, "torch")
        assert [type(plan.model).__name__ for plan in leaves] == [
            "Gaussian1D",
            "Const1D",
            "PowerLaw1D",
        ]
        assert [plan.suffix for plan in leaves] == ["_0", "_1", "_2"]
        combined = combine([np.array([1.0]), np.array([2.0]), np.array([3.0])])
        assert combined == pytest.approx([9.0])  # (1 + 2) * 3

    @pytest.mark.parametrize(
        ("op", "expected"),
        [("+", 5.0), ("-", -1.0), ("*", 6.0), ("/", 2.0 / 3.0)],
    )
    def test_every_supported_operator_composes(self, op: str, expected: float) -> None:
        if op == "+":
            model = PowerLaw1D(1.0, 1.0, 1.0) + PowerLaw1D(1.0, 1.0, 1.0)
        elif op == "-":
            model = PowerLaw1D(1.0, 1.0, 1.0) - PowerLaw1D(1.0, 1.0, 1.0)
        elif op == "*":
            model = PowerLaw1D(1.0, 1.0, 1.0) * PowerLaw1D(1.0, 1.0, 1.0)
        else:
            model = PowerLaw1D(1.0, 1.0, 1.0) / PowerLaw1D(1.0, 1.0, 1.0)
        _, combine = plan_translation(model, "torch")
        combined = combine([np.array([2.0]), np.array([3.0])])
        assert combined == pytest.approx([expected])

    def test_an_operator_with_no_elementwise_meaning_is_refused_by_name(self) -> None:
        piped = Gaussian1D(1.0, 2.0, 0.5) | Const1D(0.3)
        with pytest.raises(CapabilityError, match="'\\|'"):
            plan_translation(piped, "torch")


class TestTheRefusals:
    """The two refusals :mod:`ampere.core.astropy_translations` adds beyond W4.6's."""

    def test_unsupported_operator_refusal_names_the_operator_and_the_backend(self) -> None:
        model = Gaussian1D(1.0, 2.0, 0.5) & Const1D(0.3)
        message = str(unsupported_operator_refusal("torch", model, "&"))
        assert "ampere.backends.torch.from_astropy()" in message
        assert "'&'" in message
        assert "ampere.core.from_astropy()" in message

    def test_unsupported_operator_refusal_is_a_capability_error(self) -> None:
        refusal = unsupported_operator_refusal("jax", Const1D(0.3), "|")
        assert isinstance(refusal, CapabilityError)
        assert isinstance(refusal, NotImplementedError)

    def test_tied_parameter_refusal_names_the_parameter_and_both_reasons(self) -> None:
        message = str(tied_parameter_refusal("jax", ["mean_0"]))
        assert "ampere.backends.jax.from_astropy()" in message
        assert "mean_0" in message
        assert "traced" in message
        assert "ampere.core.from_astropy()" in message

    def test_tied_parameter_refusal_is_a_capability_error(self) -> None:
        refusal = tied_parameter_refusal("torch", ["mean_0"])
        assert isinstance(refusal, CapabilityError)
        assert isinstance(refusal, NotImplementedError)


WAVELENGTH = np.geomspace(1.0, 20.0, 8)


class TestTranslateAstropyParametersIsExposed:
    """The helper W4.6's adapter already had, exposed so W4.7 reuses it verbatim."""

    def test_it_agrees_with_the_adapter_it_backs(self) -> None:
        """Building it directly gives the same parameters ``from_astropy`` registers."""
        model = Gaussian1D(1.0, 2.0, 0.5)
        for name in model.param_names:
            getattr(model, name).bounds = (0.01, 10.0)
        parameters, ties = translate_astropy_parameters(model)
        adapted = from_astropy(model, grid=WAVELENGTH * u.micron, kind=Spectrum)
        assert [p.name for p in parameters] == list(adapted.parameters.free_names)
        assert ties == {}

    def test_a_tied_parameter_comes_back_as_a_tie_not_a_parameter(self) -> None:
        compound = Gaussian1D(3.0, 6.0, 1.5) + Const1D(0.3)
        compound.mean_0.tied = lambda m: float(m.amplitude_0.value) * 2.0
        for name in ("amplitude_0", "stddev_0", "amplitude_1"):
            getattr(compound, name).bounds = (0.01, 12.0)
        parameters, ties = translate_astropy_parameters(compound)
        assert "mean_0" not in [p.name for p in parameters]
        assert sorted(ties) == ["mean_0"]

    def test_an_unknown_prior_name_is_refused_by_name(self) -> None:
        from ampere.core.exceptions import ParameterError

        model = PowerLaw1D(amplitude=2.0, x_0=1.0, alpha=1.0)
        with pytest.raises(ParameterError, match="does not declare"):
            translate_astropy_parameters(model, {"nonexistent": 1.0})


class TestTheAdapterExposesItsBookkeeping:
    """The three read-only views a native translation reuses (§5, W4.7)."""

    def test_grids_factor_and_template_are_none_before_a_grid_is_known(self) -> None:
        model = PowerLaw1D(amplitude=2.0, x_0=1.0, alpha=1.0)
        model.amplitude.bounds = (0.1, 10.0)
        model.x_0.fixed = True
        model.alpha.bounds = (0.1, 5.0)
        adapted = from_astropy(model, kind=Spectrum)
        assert adapted.grids is None
        assert adapted.template is None
        assert adapted.factor == 1.0

    def test_grids_factor_and_template_are_populated_once_a_grid_is_given(self) -> None:
        model = PowerLaw1D(amplitude=2.0 * u.Jy, x_0=1.0 * u.micron, alpha=1.0)
        model.amplitude.bounds = (0.1, 10.0)
        model.x_0.fixed = True
        model.alpha.bounds = (0.1, 5.0)
        adapted = from_astropy(model, grid=WAVELENGTH * u.micron, kind=Spectrum, output_unit=u.Jy)
        assert adapted.grids is not None
        assert adapted.template is not None
        assert adapted.factor == pytest.approx(1.0)
