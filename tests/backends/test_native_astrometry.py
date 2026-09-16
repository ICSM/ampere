"""The astrometric twins on torch and jax: the rows the neutral battery cannot reach.

``tests/conformance/test_astrometry.py`` already holds both twins to the §4
contracts and to the closed-form ephemeris, once per registered backend. What
this file covers is what that backend-neutral battery has no vocabulary for:

* the **gradient**, which is the whole reason a native twin exists — a
  reflex-orbit offset differentiated with respect to its own orbital
  parameters, through each backend's own autodiff;
* the **inherited declaration**, checked by direct comparison rather than by
  trusting the inheriting pattern: both twins' published requirement and
  epoch buffer are the reference class's, unchanged;
* a **NUTS fit of the two-channel composition**, including the flexible
  likelihood (``GaussianProcessNoise`` + ``QuasisepGP``) — the O(N) path this
  modality is chosen to exercise;
* the ``EpochSample`` identity surviving a device move (``.to()``).

**Parametrised over the backends installed here**, in
``tests/backends/test_native_interferometry.py``'s shape: each row runs on
whichever of torch and jax this environment has.
"""

from __future__ import annotations

import dataclasses
import importlib
import importlib.util
from typing import Any

import astropy.units as u
import numpy as np
import pytest

import astrometry_fixtures as kit_data
from ampere.backends.reference import astrometry as reference
from ampere.core import AxisRequirement

#: How closely a native twin must reproduce the reference offset. Both
#: backends compute the same four-line expression in float64.
TRANSLATION_RTOL = 1e-12


@dataclasses.dataclass(frozen=True)
class Kit:
    """One backend's astrometric vocabulary, by name, so no row names a library."""

    name: str
    backend: Any
    astro: Any

    def numpy(self, values: Any) -> np.ndarray:
        """*values* back on the numpy side, for comparison against the reference."""
        if self.name == "torch":
            return self.backend.to_numpy(values)
        return np.asarray(values)

    def gradient(self, function: Any, at: float) -> float:
        """``d function / d x`` at *at*, through this backend's own autodiff."""
        if self.name == "torch":
            torch = importlib.import_module("torch")
            argument = torch.tensor(float(at), dtype=torch.float64, requires_grad=True)
            function(argument).backward()
            assert argument.grad is not None
            return float(argument.grad)
        jax = importlib.import_module("jax")
        return float(jax.grad(function)(float(at)))


def _installed_kits() -> list[Kit]:
    found: list[Kit] = []
    for name in ("jax", "torch"):
        try:
            backend = importlib.import_module(f"ampere.backends.{name}")
        except ImportError:  # the extra is not installed in this environment
            continue
        if name == "jax":
            backend.configure_x64()
        twins = importlib.import_module(f"ampere.backends.{name}.astrometry")
        found.append(Kit(name, backend, twins))
    return found


KITS = _installed_kits()

pytestmark = pytest.mark.skipif(
    not KITS, reason="no modern backend installed; the native astrometry twins need one"
)


@pytest.fixture(scope="module", params=[kit.name for kit in KITS])
def kit(request: Any) -> Kit:
    return next(found for found in KITS if found.name == request.param)


# ---------------------------------------------------------------------------
# Agreement with the reference path
# ---------------------------------------------------------------------------


class TestAgreementWithTheReferencePath:
    """One orbit, evaluated by two backends, must agree to translation precision."""

    def test_both_channels_agree(self, kit: Kit) -> None:
        native = kit_data.truth_model(kit.astro).evaluate()
        expected = kit_data.truth_model(reference).evaluate()
        for channel in ("ra", "dec"):
            got = kit.numpy(native[channel].values)
            want = np.asarray(expected[channel].values)
            assert np.allclose(got, want, rtol=TRANSLATION_RTOL, atol=1e-12)

    def test_the_epoch_sample_requirement_is_inherited_unchanged(self, kit: Kit) -> None:
        native = kit.astro.EpochSample(kit_data.EPOCHS)
        (requirement,) = native.requirements()
        assert isinstance(requirement, AxisRequirement)
        assert requirement.axis == "time"
        assert np.array_equal(requirement.points, kit_data.EPOCHS)

    def test_epoch_sample_apply_is_still_the_identity(self, kit: Kit) -> None:
        observed = kit_data.observed_channel()
        step = kit.astro.EpochSample.from_observed(observed)
        assert step(observed, None) is observed


# ---------------------------------------------------------------------------
# The gradient
# ---------------------------------------------------------------------------


class TestTheGradient:
    """The whole reason a native twin exists."""

    def test_the_ra_offset_gradient_with_respect_to_pmra_is_exact(self, kit: Kit) -> None:
        model = kit.astro.ReflexOrbit(kit_data.EPOCHS * u.day, **{**kit_data.ORBIT, "pmra": 0.0})
        compiled = model.compile_for({})

        def offset_sum(pmra: Any) -> Any:
            return compiled.native_flux("ra", {**kit_data.ORBIT, "pmra": pmra}).sum()

        got = kit.gradient(offset_sum, 1.0)
        expected = float(np.sum(kit_data.EPOCHS / 365.25))
        assert abs(got - expected) < 1e-9

    def test_the_dec_offset_gradient_with_respect_to_amp_dec_is_exact(self, kit: Kit) -> None:
        model = kit.astro.ReflexOrbit(kit_data.EPOCHS * u.day, **kit_data.ORBIT)
        compiled = model.compile_for({})

        def offset_sum(amp_dec: Any) -> Any:
            return compiled.native_flux("dec", {**kit_data.ORBIT, "amp_dec": amp_dec}).sum()

        got = kit.gradient(offset_sum, kit_data.ORBIT["amp_dec"])
        cycle = 2.0 * np.pi * kit_data.EPOCHS / kit_data.ORBIT["period"] + kit_data.ORBIT["phase"]
        expected = float(np.sum(np.cos(cycle)))
        assert abs(got - expected) < 1e-9


# ---------------------------------------------------------------------------
# NUTS on the two-channel composition
# ---------------------------------------------------------------------------


class TestNutsOnTheComposition:
    """The two-channel problem fits under NUTS, independent noise and the flexible likelihood."""

    def test_independent_noise_recovers_the_proper_motion(self, kit: Kit) -> None:
        from ampere.inference import NUTSEngine

        problem = kit_data.two_channel_problem(kit.backend, kit.astro, gp=False)
        run = NUTSEngine(problem).run(80, warmup=120, chains=1, progress=False)
        posterior = run["posterior"]
        for name in kit_data.TRUTH:
            values = np.asarray(posterior[name].values, dtype=float)
            assert values.size == 80
            assert np.all(np.isfinite(values))

    def test_the_flexible_likelihood_also_lowers_natively(self, kit: Kit) -> None:
        from ampere.inference import NUTSEngine

        problem = kit_data.two_channel_problem(kit.backend, kit.astro, gp=True)
        run = NUTSEngine(problem).run(40, warmup=60, chains=1, progress=False)
        posterior = run["posterior"]
        for name in kit_data.TRUTH:
            values = np.asarray(posterior[name].values, dtype=float)
            assert np.all(np.isfinite(values))


needs_torch = pytest.mark.skipif(
    importlib.util.find_spec("torch") is None, reason="needs ampere[torch]"
)


@needs_torch
class TestTorchDeviceMove:
    """``.to()`` re-declares placement; the identity step has nothing to move.

    Torch-only: jax's device placement is immutable, fixed at construction
    (``device=``), with no ``.to()`` verb — a fact about jax's own convention
    (:mod:`ampere.backends.jax.instrument`), not about this modality.
    """

    def test_epoch_sample_to_returns_itself(self) -> None:
        from ampere.backends.torch import astrometry as torch_astrometry

        step = torch_astrometry.EpochSample(kit_data.EPOCHS)
        assert step.to() is step
