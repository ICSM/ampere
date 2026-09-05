"""Tests for W2.6 — the lowering registry (``ampere.core.lowering``).

``docs/design/lowering.md`` §12.8 rules the mechanism; ``WORK_ITEMS.md``'s
W2.6 entry fixes the acceptance criteria this file exercises directly:

1. a user-registered prior family lowers natively on a backend, and its
   provenance record says so (:class:`TestProvenanceStamping`);
2. the overwrite refusal and the ``override=True`` escape are tested
   (:class:`TestOverwriteRefusal`);
3. the registrant battery runs "from the docs example" —
   ``tests/core/test_spec_doctests.py``'s
   ``test_module_docstring_examples_run[lowering]`` executes
   :mod:`ampere.core.lowering`'s module docstring, which registers a prior
   family on a stub backend, runs the battery against a correct and a wrong
   registration, and stamps the result in provenance. This file adds
   focused, non-doctest coverage of the same behaviour so failures are
   easy to localise.

No real backend exists yet (W2.4/W2.5 land the torch/jax rows); every test
here uses backend names it invents itself (``"stub-*"``), exactly as the item
brief anticipates, and scipy as the reference-vs-native agreement baseline.
"""

from __future__ import annotations

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.core import Dataset, FittingProblem, Model, Parameter, Spectrum
from ampere.core.exceptions import LoweringError, ParameterError
from ampere.core.lowering import (
    BatteryReport,
    LoweringResolution,
    lookup_bijection_lowering,
    lookup_lowering,
    provenance_entries,
    register_bijection_lowering,
    register_lowering,
    registered_lowerings,
    run_registrant_battery,
)
from ampere.core.parameter import Bijection, Logit
from ampere.results.provenance import provenance_attrs


class _CustomPowerLaw(st.rv_continuous):
    """A prior family with no native torch/numpyro equivalent -- the motivating
    case ``lowering.md`` §3.4/§12.8 write the registry for: ``pdf(x) = 2x`` on
    ``[0, 1]``.
    """

    def _pdf(self, x: np.ndarray) -> np.ndarray:
        return 2.0 * x


def _custom_power_law() -> object:
    """A frozen instance of :class:`_CustomPowerLaw`, describable and lowerable."""
    return _CustomPowerLaw(a=0.0, b=1.0, name="w26_custom_power")()


class _NativePowerLaw:
    """A correct "native" stand-in: recomputes the exact log-density from scratch
    rather than delegating back to scipy, so the battery is checking something.
    """

    def log_prob(self, x: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        return np.log(2.0) + np.log(x)


class _WrongNativePowerLaw:
    """A deliberately wrong "native" stand-in (uniform density, not 2x)."""

    def log_prob(self, x: np.ndarray) -> np.ndarray:
        return np.zeros_like(np.asarray(x, dtype=float))


def _build_native_power_law(spec: object) -> _NativePowerLaw:
    del spec  # the family has no free parameters
    return _NativePowerLaw()


def _build_wrong_power_law(spec: object) -> _WrongNativePowerLaw:
    del spec
    return _WrongNativePowerLaw()


# ---------------------------------------------------------------------------
# Registration and lookup
# ---------------------------------------------------------------------------


class TestRegistrationAndLookup:
    def test_register_then_lookup_round_trips(self) -> None:
        resolution = register_lowering("w26_test_family_a", "stub-a", _build_native_power_law)
        assert isinstance(resolution, LoweringResolution)
        assert resolution.kind == "prior"
        assert resolution.name == "w26_test_family_a"
        assert resolution.backend == "stub-a"
        assert resolution.builtin is False
        assert resolution.constructor is _build_native_power_law

        found = lookup_lowering("w26_test_family_a", "stub-a")
        assert found == resolution

    def test_lookup_missing_names_family_and_backend_with_remedies(self) -> None:
        with pytest.raises(LoweringError) as excinfo:
            lookup_lowering("w26_never_registered", "nowhere")
        message = str(excinfo.value)
        assert "w26_never_registered" in message
        assert "nowhere" in message
        # LoweringError's standard remedies, so a caller catching by type
        # always sees the same three options (lowering.md §3.4).
        assert "register your own lowering" in message
        assert excinfo.value.family == "w26_never_registered"
        assert excinfo.value.backend == "nowhere"

    def test_registered_lowerings_filters_by_kind_and_backend(self) -> None:
        register_lowering("w26_filter_family", "stub-filter", _build_native_power_law)
        rows = registered_lowerings(kind="prior", backend="stub-filter")
        assert any(r.name == "w26_filter_family" for r in rows)
        assert all(r.kind == "prior" and r.backend == "stub-filter" for r in rows)

    def test_register_lowering_rejects_bad_arguments(self) -> None:
        with pytest.raises(ParameterError):
            register_lowering("", "stub", _build_native_power_law)
        with pytest.raises(ParameterError):
            register_lowering("w26_x", "", _build_native_power_law)
        with pytest.raises(ParameterError):
            register_lowering("w26_x", "stub", constructor="not callable")  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# The overwrite refusal, and override=True
# ---------------------------------------------------------------------------


class TestOverwriteRefusal:
    def test_second_registration_without_override_is_refused(self) -> None:
        register_lowering("w26_overwrite_family", "stub-ow", _build_native_power_law)
        with pytest.raises(LoweringError) as excinfo:
            register_lowering("w26_overwrite_family", "stub-ow", _build_wrong_power_law)
        message = str(excinfo.value)
        assert "already registered" in message
        assert "override=True" in message
        # the still-registered row is the original, not the rejected one
        assert (
            lookup_lowering("w26_overwrite_family", "stub-ow").constructor
            is _build_native_power_law
        )

    def test_override_true_replaces_the_row(self) -> None:
        register_lowering("w26_override_family", "stub-ov", _build_native_power_law)
        replaced = register_lowering(
            "w26_override_family", "stub-ov", _build_wrong_power_law, override=True
        )
        assert replaced.constructor is _build_wrong_power_law
        assert lookup_lowering("w26_override_family", "stub-ov").constructor is (
            _build_wrong_power_law
        )

    def test_a_builtin_row_also_refuses_overwrite_without_override(self) -> None:
        """``builtin=True`` is reserved for the backends ampere ships (W2.4/W2.5),
        but the refusal must cover it too -- "no overwrite of a built-in *or*
        existing row" (WORK_ITEMS.md W2.6) -- and the message should say which
        kind of row it refused to clobber.
        """
        register_lowering(
            "w26_builtin_family", "stub-builtin", _build_native_power_law, builtin=True
        )
        with pytest.raises(LoweringError) as excinfo:
            register_lowering("w26_builtin_family", "stub-builtin", _build_wrong_power_law)
        assert "built-in" in str(excinfo.value)

        # override=True still works on a built-in row.
        register_lowering(
            "w26_builtin_family",
            "stub-builtin",
            _build_wrong_power_law,
            override=True,
            builtin=True,
        )
        assert (
            lookup_lowering("w26_builtin_family", "stub-builtin").constructor
            is _build_wrong_power_law
        )


# ---------------------------------------------------------------------------
# The bijection slot
# ---------------------------------------------------------------------------


class _DoubledLogit(Bijection):
    """A custom Bijection: an ordinary Logit with the bounds doubled -- deliberately
    not one of ``Identity``/``Log``/``Logit``, so it exercises the "unsupported --
    raises at lowering" table row (``lowering.md`` §4) this registration closes.
    """

    def __init__(self, lower: float, upper: float) -> None:
        self._inner = Logit(lower=lower, upper=upper)

    def constrain(self, y: object) -> np.ndarray:
        return self._inner.constrain(y)

    def unconstrain(self, x: object) -> np.ndarray:
        return self._inner.unconstrain(x)

    def log_abs_det_jacobian(self, y: object) -> np.ndarray:
        return self._inner.log_abs_det_jacobian(y)


def _build_doubled_logit_native(bijection: _DoubledLogit) -> dict[str, float]:
    return {"lower": bijection._inner.lower, "upper": bijection._inner.upper}


class TestBijectionSlot:
    def test_register_and_lookup_by_class_or_instance(self) -> None:
        resolution = register_bijection_lowering(
            _DoubledLogit, "stub-bij", _build_doubled_logit_native
        )
        assert resolution.kind == "bijection"
        assert resolution.name == "_DoubledLogit"

        instance = _DoubledLogit(0.0, 1.0)
        by_class = lookup_bijection_lowering(_DoubledLogit, "stub-bij")
        by_instance = lookup_bijection_lowering(instance, "stub-bij")
        assert by_class == by_instance == resolution

        native = by_instance.constructor(instance)
        assert native == {"lower": 0.0, "upper": 1.0}

    def test_overwrite_refused_for_bijection_slot_too(self) -> None:
        register_bijection_lowering(_DoubledLogit, "stub-bij-ow", _build_doubled_logit_native)
        with pytest.raises(LoweringError):
            register_bijection_lowering(_DoubledLogit, "stub-bij-ow", _build_doubled_logit_native)
        register_bijection_lowering(
            _DoubledLogit, "stub-bij-ow", _build_doubled_logit_native, override=True
        )

    def test_lookup_missing_bijection_names_the_class(self) -> None:
        with pytest.raises(LoweringError) as excinfo:
            lookup_bijection_lowering(_DoubledLogit, "nowhere-bij")
        assert "_DoubledLogit" in str(excinfo.value)

    def test_register_bijection_lowering_rejects_non_bijection(self) -> None:
        with pytest.raises(ParameterError):
            register_bijection_lowering(object, "stub-bij", lambda b: b)  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# The registrant battery
# ---------------------------------------------------------------------------


class TestRegistrantBattery:
    def test_battery_passes_for_a_correct_registration(self) -> None:
        register_lowering(
            "w26_custom_power", "stub-battery-ok", _build_native_power_law, override=True
        )
        report = run_registrant_battery(_custom_power_law(), "stub-battery-ok")
        assert isinstance(report, BatteryReport)
        assert report.passed
        assert report.family == "w26_custom_power"
        assert report.max_abs_error < 1e-9

    def test_battery_fails_for_a_wrong_registration_with_a_legible_message(self) -> None:
        register_lowering(
            "w26_custom_power", "stub-battery-bad", _build_wrong_power_law, override=True
        )
        with pytest.raises(LoweringError) as excinfo:
            run_registrant_battery(_custom_power_law(), "stub-battery-bad")
        message = str(excinfo.value)
        assert "registrant battery failed" in message
        # The message names the point and both the native and reference values,
        # not just "disagreement" -- a registrant should be able to see why.
        assert "reference (scipy)" in message
        assert "w26_custom_power" in message

    def test_battery_reports_a_legible_error_without_log_prob(self) -> None:
        register_lowering(
            "w26_custom_power", "stub-battery-nolp", lambda spec: object(), override=True
        )
        with pytest.raises(LoweringError) as excinfo:
            run_registrant_battery(_custom_power_law(), "stub-battery-nolp")
        assert "log_prob" in str(excinfo.value)

    def test_battery_honours_explicit_log_prob_override(self) -> None:
        register_lowering(
            "w26_custom_power", "stub-battery-explicit", lambda spec: "opaque", override=True
        )
        report = run_registrant_battery(
            _custom_power_law(),
            "stub-battery-explicit",
            log_prob=lambda native, xs: np.log(2.0) + np.log(xs),
        )
        assert report.passed


# ---------------------------------------------------------------------------
# Provenance stamping (Accept: "its provenance record says so")
# ---------------------------------------------------------------------------


class _PowerLawModel(Model):
    def __init__(self) -> None:
        self.register_buffer("grid", np.array([0.2, 0.5, 0.8]), unit=u.um)
        self.register_parameter(Parameter("norm", _custom_power_law()))

    def evaluate(self, **values: float) -> Spectrum:
        ctx = self.context(values)
        return Spectrum(ctx["grid"] * u.um, ctx["norm"] * ctx["grid"] * u.Jy)


def _toy_problem() -> FittingProblem:
    model = _PowerLawModel()
    grid = np.array([0.2, 0.5, 0.8])
    observed = Spectrum(grid * u.um, [0.1, 0.25, 0.4] * u.Jy, uncertainty=[0.05, 0.05, 0.05] * u.Jy)
    return FittingProblem(model, [Dataset(observed)], seed=20260906)


class TestProvenanceStamping:
    def test_a_user_registered_family_lowers_natively_and_is_stamped(self) -> None:
        # "lowers natively on a backend": resolve the registry and invoke the
        # registered constructor, exactly as a backend (W2.4/W2.5) will.
        register_lowering("w26_custom_power", "stub-prov", _build_native_power_law, override=True)
        prior = _custom_power_law()
        from ampere.core.parameter import describe_prior

        spec = describe_prior(prior)
        resolution = lookup_lowering(spec.family, "stub-prov")
        native = resolution.constructor(spec)
        assert isinstance(native, _NativePowerLaw)
        assert np.isclose(native.log_prob(0.5), np.log(2.0) + np.log(0.5))

        problem = _toy_problem()
        entries = provenance_entries([resolution])
        attrs = provenance_attrs(
            problem, backend="stub-prov", extra={"registered_lowerings": entries}
        )

        assert "ampere_registered_lowerings" in attrs
        import json

        recorded = json.loads(attrs["ampere_registered_lowerings"])
        assert recorded == [
            {
                "kind": "prior",
                "name": "w26_custom_power",
                "backend": "stub-prov",
                "builtin": False,
                "constructor": (
                    f"{_build_native_power_law.__module__}.{_build_native_power_law.__qualname__}"
                ),
            }
        ]
        assert attrs["ampere_backend"] == "stub-prov"

    def test_builtin_rows_are_not_stamped_as_user_registered(self) -> None:
        register_lowering(
            "w26_custom_power",
            "stub-prov-builtin",
            _build_native_power_law,
            builtin=True,
            override=True,
        )
        resolution = lookup_lowering("w26_custom_power", "stub-prov-builtin")
        assert provenance_entries([resolution]) == []
