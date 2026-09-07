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

Two findings from the first W2.6 review are addressed and covered here too:
``_REGISTRY`` is module-global state exactly like
``ampere.core.likelihood._FAMILIES`` (W0.10 finding (c)'s failure class), so
``tests/core/conftest.py``'s autouse fixture now snapshots and restores it
the same way (:class:`TestFixtureRestoresRegistryBetweenTests` below proves
it); and the bijection slot is keyed on the module-qualified class name, not
the bare one, so two unrelated same-named classes cannot cross-hit each
other's row on lookup
(:meth:`TestBijectionSlot.test_same_named_classes_from_different_modules_do_not_collide`).
"""

from __future__ import annotations

import warnings
from typing import ClassVar

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.core import (
    Dataset,
    FittingProblem,
    GaussianFamily,
    IndependentNoise,
    Likelihood,
    Model,
    Parameter,
    Spectrum,
)
from ampere.core.exceptions import (
    AmpereError,
    LoweringError,
    LoweringFallbackWarning,
    ParameterError,
)
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


def _make_bijection_class(module_name: str) -> type[Bijection]:
    """A fresh class named ``_SameName``, "from" a synthetic module.

    Stands in for two unrelated third-party packages that both happen to
    define a class called ``_SameName`` -- W2.6 review Finding 2's collision
    scenario -- without needing two real module files. Only ``__module__``
    is synthetic; the class is otherwise an ordinary, independently created
    type, so ``is``-identity and qualified-name collision behave exactly as
    they would for two real packages.
    """
    cls = type(
        "_SameName",
        (Bijection,),
        {
            "__init__": lambda self, tag: setattr(self, "tag", tag),
            "constrain": lambda self, y: np.asarray(y),
            "unconstrain": lambda self, x: np.asarray(x),
            "log_abs_det_jacobian": lambda self, y: np.zeros_like(np.asarray(y, dtype=float)),
        },
    )
    cls.__module__ = module_name
    return cls


class TestBijectionSlot:
    def test_register_and_lookup_by_class_or_instance(self) -> None:
        resolution = register_bijection_lowering(
            _DoubledLogit, "stub-bij", _build_doubled_logit_native
        )
        assert resolution.kind == "bijection"
        # Keyed (and recorded) on the module-qualified name, not the bare
        # class name -- W2.6 review Finding 2.
        assert resolution.name == f"{_DoubledLogit.__module__}.{_DoubledLogit.__qualname__}"

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

    def test_same_named_classes_from_different_modules_do_not_collide(self) -> None:
        """W2.6 review Finding 2: a bare-name key would let this lookup for
        ``package_b``'s ``_SameName`` silently return ``package_a``'s row and
        hand its constructor the wrong instance. Keying on the module-qualified
        name means the second registration is not even refused as a duplicate
        -- the two rows are genuinely independent -- and each lookup finds
        only its own.
        """
        package_a = _make_bijection_class("w26_package_a")
        package_b = _make_bijection_class("w26_package_b")
        assert package_a.__name__ == package_b.__name__ == "_SameName"
        assert package_a is not package_b

        def _build_a(bijection: object) -> str:
            return "native-a"

        def _build_b(bijection: object) -> str:
            return "native-b"

        register_bijection_lowering(package_a, "stub-collision", _build_a)
        # Not a duplicate: different qualified key, so no override needed.
        register_bijection_lowering(package_b, "stub-collision", _build_b)

        resolution_a = lookup_bijection_lowering(package_a, "stub-collision")
        resolution_b = lookup_bijection_lowering(package_b, "stub-collision")
        assert resolution_a.constructor is _build_a
        assert resolution_b.constructor is _build_b
        assert resolution_a.name == "w26_package_a._SameName"
        assert resolution_b.name == "w26_package_b._SameName"

        # Human-facing messages still show the bare, ambiguous name --
        # that convention is about legibility, not storage.
        with pytest.raises(LoweringError) as excinfo:
            lookup_bijection_lowering(_make_bijection_class("w26_package_c"), "nowhere")
        assert "_SameName" in str(excinfo.value)

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
    # W2.12: the backend is a capability flag the *pieces* declare, so a
    # problem that ran on "stub-prov" is one whose model says so. Asserting it
    # through ``provenance_attrs(backend=...)`` is no longer possible, and that
    # is the point: ``ampere_backend`` is a fact now.
    BACKEND: ClassVar[str] = "stub-prov"

    def __init__(self) -> None:
        self.register_buffer("grid", np.array([0.2, 0.5, 0.8]), unit=u.um)
        self.register_parameter(Parameter("norm", _custom_power_law()))

    def evaluate(self, **values: float) -> Spectrum:
        ctx = self.context(values)
        return Spectrum(ctx["grid"] * u.um, ctx["norm"] * ctx["grid"] * u.Jy)


class _StubNoise(IndependentNoise):
    # W2.13: a noise model is a capability part now, so the dataset's default
    # (which declares "reference") would make this a two-backend problem.
    BACKEND: ClassVar[str] = "stub-prov"


def _toy_problem() -> FittingProblem:
    model = _PowerLawModel()
    grid = np.array([0.2, 0.5, 0.8])
    observed = Spectrum(grid * u.um, [0.1, 0.25, 0.4] * u.Jy, uncertainty=[0.05, 0.05, 0.05] * u.Jy)
    dataset = Dataset(observed, likelihood=Likelihood(GaussianFamily(), _StubNoise()))
    return FittingProblem(model, [dataset], seed=20260906)


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
        attrs = provenance_attrs(problem, extra={"registered_lowerings": entries})

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


# ---------------------------------------------------------------------------
# The autouse registry-restoration fixture (W2.6 review Finding 1)
# ---------------------------------------------------------------------------


class TestFixtureRestoresRegistryBetweenTests:
    """Proves ``tests/core/conftest.py``'s ``_restore_lowering_registry`` fixture
    actually restores ``ampere.core.lowering._REGISTRY`` between tests --
    the same class of leak W0.10 finding (c) found in
    ``ampere.core.likelihood._FAMILIES`` before that module got its own
    autouse snapshot/restore fixture.

    Two steps, relying on pytest's default (source-order, non-randomised)
    collection within one module -- there is no ``pytest-randomly`` or
    similar in this project's dev dependencies, and the two methods are
    adjacent and named so nothing reorders them independently of this class.
    If the fixture in ``conftest.py`` were removed (or only cleared
    ``_FAMILIES`` and not ``_REGISTRY``), step 2 would fail: it would find
    step 1's row still registered.
    """

    def test_step1_registers_a_marker_row(self) -> None:
        register_lowering("w26_fixture_proof_row", "stub-fixture-proof", _build_native_power_law)
        assert lookup_lowering("w26_fixture_proof_row", "stub-fixture-proof").constructor is (
            _build_native_power_law
        )

    def test_step2_the_marker_row_is_gone_afterwards(self) -> None:
        assert registered_lowerings(backend="stub-fixture-proof") == ()
        with pytest.raises(LoweringError):
            lookup_lowering("w26_fixture_proof_row", "stub-fixture-proof")


class TestTheSharedFallbackWarning:
    """W2.13, fold-in 8: one class, in the core, for both backends.

    W2.4 and W2.5 each invented their own (``IcdfFallbackWarning`` and a
    backend-local ``LoweringFallbackWarning``) and each said in its own
    docstring that a shared core class would be better, because a test that
    wants to assert "no native run silently computed in numpy" must be able to
    ``pytest.warns`` on one type across every backend. It lives here.
    """

    def test_it_is_a_warning_and_not_an_ampere_error(self) -> None:
        # Deliberately separate hierarchies: the sanctioned §3.6 fallback is
        # not an error, and ``except AmpereError`` must not swallow it.
        assert issubclass(LoweringFallbackWarning, UserWarning)
        assert not issubclass(LoweringFallbackWarning, AmpereError)

    def test_it_is_exported_from_ampere_core(self) -> None:
        import ampere.core as core

        assert core.LoweringFallbackWarning is LoweringFallbackWarning
        assert "LoweringFallbackWarning" in core.__all__

    def test_one_pytest_warns_catches_whatever_a_backend_raises(self) -> None:
        with pytest.warns(LoweringFallbackWarning, match="ppf"):
            warnings.warn(
                "prior_transform for 'gamma' falls back to scipy's ppf",
                LoweringFallbackWarning,
                stacklevel=1,
            )
