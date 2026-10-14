"""The line-integration twins on torch and jax: the rows the neutral battery cannot reach.

``tests/conformance/test_line_fluxes.py`` already holds both twins to the §4
contracts and to the closed form, once per registered backend. What this file
covers is what that backend-neutral battery has no vocabulary for:

* the **inherited declaration**, checked by direct comparison rather than by
  trusting the inheriting pattern: both twins' tabulation, weights and rest
  wavelengths are the reference class's, bit for bit;
* the **native surface** ``apply_flux`` on bare arrays — on a negotiated grid
  wider than the step's own and in another order — against the reference
  step's ``apply`` on the same numbers;
* the **gradient**, which is the reason a native twin exists: the line
  integrals differentiated with respect to a scale on the flux, through each
  backend's own autodiff, against the reference integral;
* a **device declaration**: torch's ``.to()`` moves the weights tensor, jax's
  ``device=`` keyword is accepted and declared.

**Parametrised over the backends installed here**, in
``tests/backends/test_native_astrometry.py``'s shape.
"""

from __future__ import annotations

import dataclasses
import importlib
from typing import Any

import astropy.units as u
import numpy as np
import pytest

from ampere.backends.reference import LineIntegration as ReferenceLineIntegration
from ampere.core import Spectrum

#: How closely a native twin must reproduce the reference integrals: both are
#: one float64 ``matmul`` of the same weights.
TRANSLATION_RTOL = 1e-12

TRANSITIONS = ["[Ne II] 12.81", "[S IV] 10.51", "[Ne III] 15.56", "H I 7-6"]
REST = np.array([12.8135, 10.5105, 15.5551, 12.3719])


@dataclasses.dataclass(frozen=True)
class Kit:
    """One backend's twin, by name, so no row names a library."""

    name: str
    backend: Any

    def array(self, values: Any) -> Any:
        if self.name == "torch":
            torch = importlib.import_module("torch")
            return torch.as_tensor(np.array(values, dtype=float), dtype=torch.float64)
        jnp = importlib.import_module("jax.numpy")
        return jnp.asarray(np.asarray(values, dtype=float))

    def numpy(self, values: Any) -> np.ndarray:
        if self.name == "torch":
            return self.backend.to_numpy(values)
        return np.asarray(values)

    def gradient(self, function: Any, at: float) -> np.ndarray:
        """``d function / d x`` at *at* for a vector-valued *function*, summed."""
        if self.name == "torch":
            torch = importlib.import_module("torch")
            argument = torch.tensor(float(at), dtype=torch.float64, requires_grad=True)
            function(argument).sum().backward()
            assert argument.grad is not None
            return np.asarray(float(argument.grad))
        jax = importlib.import_module("jax")
        return np.asarray(jax.grad(lambda x: function(x).sum())(float(at)))


def _installed_kits() -> list[Kit]:
    found: list[Kit] = []
    for name in ("jax", "torch"):
        try:
            backend = importlib.import_module(f"ampere.backends.{name}")
        except ImportError:  # the extra is not installed in this environment
            continue
        if name == "jax":
            backend.configure_x64()
        found.append(Kit(name, backend))
    return found


KITS = _installed_kits()

pytestmark = pytest.mark.skipif(not KITS, reason="neither torch nor jax is installed here")


@pytest.fixture(params=KITS, ids=[kit.name for kit in KITS])
def kit(request: pytest.FixtureRequest) -> Kit:
    return request.param


def build(cls: Any, **overrides: Any) -> Any:
    arguments: dict[str, Any] = {
        "window": 300.0,
        "continuum": "linear",
        "samples_per_window": 33,
        **overrides,
    }
    return cls(TRANSITIONS, REST, **arguments)


def wide_grid(step: Any) -> np.ndarray:
    """The step's tabulation inside a wider grid, shuffled: what a union hands it."""
    extra = np.linspace(9.0, 17.0, 41)
    grid = np.concatenate([step.tabulation(), extra])
    return np.random.default_rng(3).permutation(grid)


def flux_on(grid: np.ndarray) -> np.ndarray:
    """A continuum plus a bump at each line: something the continuum does not remove."""
    flux = 1.0 + 0.05 * grid
    for rest in REST:
        flux = flux + 3.0 * np.exp(-0.5 * ((grid - rest) / (rest * 2.0e-4)) ** 2)
    return flux


class TestTheInheritedDeclaration:
    def test_the_tabulation_weights_and_rest_wavelengths_are_the_references(self, kit: Kit) -> None:
        twin = build(kit.backend.LineIntegration)
        reference = build(ReferenceLineIntegration)
        assert twin.BACKEND == kit.name
        assert twin.DIFFERENTIABLE and twin.BATCHABLE
        assert np.array_equal(twin.tabulation(), reference.tabulation())
        assert np.array_equal(twin._data("weights"), reference._data("weights"))
        assert np.array_equal(twin.rest_wavelengths(), reference.rest_wavelengths())
        assert twin.describe() == reference.describe()
        (requirement,) = twin.requirements()
        assert np.array_equal(requirement.points, reference.tabulation())


class TestTheNativeSurface:
    def test_apply_flux_on_a_wider_shuffled_grid_is_the_reference_apply(self, kit: Kit) -> None:
        twin = build(kit.backend.LineIntegration)
        reference = build(ReferenceLineIntegration)
        grid = wide_grid(reference)
        flux = flux_on(grid)
        order = np.argsort(grid)
        expected = reference(Spectrum(grid[order] * u.um, flux[order] * u.Jy), None)
        got, coordinates = twin.apply_flux(kit.array(flux), kit.array(grid), {})
        np.testing.assert_allclose(kit.numpy(got), expected.values, rtol=TRANSLATION_RTOL)
        np.testing.assert_array_equal(kit.numpy(coordinates), REST)

    def test_the_contract_apply_is_the_reference_apply(self, kit: Kit) -> None:
        twin = build(kit.backend.LineIntegration)
        reference = build(ReferenceLineIntegration)
        grid = np.sort(wide_grid(reference))
        spectrum = Spectrum(grid * u.um, flux_on(grid) * u.Jy)
        got = twin(spectrum, None)
        expected = reference(spectrum, None)
        assert list(got.transitions) == TRANSITIONS
        assert got.unit == expected.unit
        np.testing.assert_allclose(got.values, expected.values, rtol=TRANSLATION_RTOL)

    def test_the_line_integrals_take_a_gradient(self, kit: Kit) -> None:
        """``d/da Σ_k F_k(a·f) = Σ_k F_k(f)``: the step is linear, so the gradient is exact."""
        twin = build(kit.backend.LineIntegration)
        grid = twin.tabulation()
        flux = flux_on(grid)
        reference = build(ReferenceLineIntegration)(Spectrum(grid * u.um, flux * u.Jy), None)
        base = kit.array(flux)
        coordinates = kit.array(grid)

        def integrals(scale: Any) -> Any:
            return twin.apply_flux(scale * base, coordinates, {})[0]

        gradient = kit.gradient(integrals, 2.0)
        np.testing.assert_allclose(gradient, reference.values.sum(), rtol=1e-10)


class TestTheDevice:
    def test_the_device_is_declared(self, kit: Kit) -> None:
        if kit.name == "torch":
            twin = build(kit.backend.LineIntegration)
            assert twin.to(device="cpu") is twin
            assert twin.DEVICE == "cpu"
            assert str(twin.tensors.get_buffer("weights").device) == "cpu"
        else:
            twin = build(kit.backend.LineIntegration, device="cpu")
            assert twin.DEVICE == "cpu"
