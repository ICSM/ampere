"""The toy model: one physics, three libraries, and a gradient where claimed.

The three variants are held to each other here rather than trusted to agree.
They share their declaration by import — the priors, the line positions, the
reference wavelength all come from :mod:`examples.m2_misspecification.model` —
so this module's job is the half that cannot be shared: the arithmetic, and the
capabilities each variant declares because of it.

It also records W2.10's structural finding as an executable statement: the
differentiable backends lower a model by **duck typing** on ``grid`` and
``flux``, so a hand-written model gets NUTS without shipping in
``ampere.backends``.
"""

from __future__ import annotations

import numpy as np
import pytest

from ampere.core import Model
from examples.m2_misspecification import generators, study
from examples.m2_misspecification.model import (
    PARAMETER_NAMES,
    PRIOR_LIMITS,
    TRUTH,
    AbsorptionLines,
    flux_at,
)

WAVELENGTH = generators.wavelength_grid(200)
OFF_TRUTH = {"A": 0.93, "B": 1.7, "d1": 0.21, "d2": 0.06}


def test_the_declaration_is_the_paper_studys() -> None:
    model = AbsorptionLines(WAVELENGTH)
    assert model.parameters.free_names == PARAMETER_NAMES
    assert PRIOR_LIMITS == {"A": (0.5, 1.5), "B": (-5.0, 12.0), "d1": (0.0, 0.5), "d2": (0.0, 0.5)}
    assert "wavelength" in model.buffers


def test_evaluate_matches_the_free_function() -> None:
    model = AbsorptionLines(WAVELENGTH)
    result = model.evaluate(**TRUTH)
    assert np.allclose(result["default"].values, flux_at(WAVELENGTH, **TRUTH))
    assert np.allclose(model.flux("default", TRUTH), flux_at(WAVELENGTH, **TRUTH))


def test_the_lines_are_absorption_and_the_continuum_is_linear() -> None:
    """The physics, stated as two properties rather than as a golden array."""
    flux = flux_at(WAVELENGTH, **TRUTH)
    continuum = flux_at(WAVELENGTH, A=TRUTH["A"], B=TRUTH["B"], d1=0.0, d2=0.0)
    assert np.all(flux <= continuum + 1e-12)
    assert np.allclose(np.diff(continuum, n=2), 0.0, atol=1e-12)
    # Deeper lines remove more flux, and each line acts where it is.
    deeper = flux_at(WAVELENGTH, A=TRUTH["A"], B=TRUTH["B"], d1=0.3, d2=TRUTH["d2"])
    assert float(np.sum(continuum - deeper)) > float(np.sum(continuum - flux))


def test_the_reference_variant_declares_the_reference_backend() -> None:
    model = AbsorptionLines(WAVELENGTH)
    assert model.BACKEND == "reference"
    assert model.DIFFERENTIABLE is False
    assert isinstance(model, Model)


def test_a_bad_grid_is_refused() -> None:
    with pytest.raises(ValueError, match="one-dimensional, non-empty"):
        AbsorptionLines(np.zeros((3, 3)))


# ---------------------------------------------------------------------------
# The accelerated variants. Skipped where the extra is absent, which is how
# ``dev`` runs this module at all.
# ---------------------------------------------------------------------------


def test_torch_variant_agrees_with_numpy_and_differentiates() -> None:
    torch = pytest.importorskip("torch")
    from examples.m2_misspecification.model_torch import (
        AbsorptionLines as TorchAbsorptionLines,
        flux_tensor,
    )

    model = TorchAbsorptionLines(WAVELENGTH)
    assert model.BACKEND == "torch"
    assert model.DIFFERENTIABLE is True
    assert model.parameters.free_names == PARAMETER_NAMES
    tensor = model.flux("default", OFF_TRUTH)
    expected = flux_at(WAVELENGTH, **OFF_TRUTH)
    assert np.allclose(tensor.detach().numpy(), expected, rtol=0, atol=1e-12)
    assert np.allclose(model.evaluate(**OFF_TRUTH)["default"].values, expected)

    grid = torch.as_tensor(WAVELENGTH, dtype=torch.float64)
    values = {
        name: torch.tensor(OFF_TRUTH[name], dtype=torch.float64, requires_grad=True)
        for name in PARAMETER_NAMES
    }
    flux_tensor(grid, **values).sum().backward()
    for name, value in values.items():
        assert value.grad is not None and bool(torch.isfinite(value.grad)), name


def test_jax_variant_agrees_with_numpy_and_differentiates() -> None:
    jax = pytest.importorskip("jax")
    from ampere.backends.jax import configure_x64

    configure_x64()
    from examples.m2_misspecification.model_jax import (
        AbsorptionLines as JaxAbsorptionLines,
        flux_jax,
    )

    model = JaxAbsorptionLines(WAVELENGTH)
    assert model.BACKEND == "jax"
    assert model.DIFFERENTIABLE is True
    assert model.parameters.free_names == PARAMETER_NAMES
    assert np.allclose(
        np.asarray(model.flux("default", OFF_TRUTH)),
        flux_at(WAVELENGTH, **OFF_TRUTH),
        rtol=0,
        atol=1e-12,
    )

    def total(values: dict[str, float]) -> object:
        return flux_jax(WAVELENGTH, **values).sum()

    gradient = jax.grad(total)({name: OFF_TRUTH[name] for name in PARAMETER_NAMES})
    for name, value in gradient.items():
        assert np.isfinite(float(value)), name


# ---------------------------------------------------------------------------
# The finding: a hand-written model lowers, because the lowering duck-types.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("backend", ["torch", "jax"])
def test_a_hand_written_model_lowers_into_a_differentiable_density(backend: str) -> None:
    """W2.10's structural finding, as an assertion rather than as prose.

    Neither ``ampere.backends.torch.problem`` nor its jax twin knows this class
    exists. Both find its native surface by ``hasattr(model, "flux")``, so the
    problem lowers, the density is finite, and the gradient is finite — which is
    the whole reason the toy model did not have to ship inside a backend.
    """
    pytest.importorskip(backend)
    module = pytest.importorskip(f"ampere.backends.{backend}")
    if backend == "jax":
        module.configure_x64()
    data = generators.generate("strong_smooth", size=200)
    problem = study.build_problem(data, backend=backend, likelihood="flexible")
    assert problem.backend == backend
    assert problem.differentiable is True

    lowered = module.lower_problem(problem)
    unconstrained = np.zeros(lowered.free_size)
    if backend == "jax":
        import jax

        value = float(lowered.log_prob_unconstrained(unconstrained))
        gradient = np.asarray(jax.grad(lowered.log_prob_unconstrained)(unconstrained))
    else:
        import torch

        point = torch.zeros(lowered.free_size, dtype=torch.float64, requires_grad=True)
        density = lowered.log_prob_unconstrained(point)
        density.backward()
        value = float(density.detach())
        gradient = point.grad.numpy()
    assert np.isfinite(value)
    assert gradient.shape == (lowered.free_size,)
    assert np.all(np.isfinite(gradient))
    assert np.any(gradient != 0.0)


@pytest.mark.parametrize("backend", ["torch", "jax"])
def test_the_three_backends_score_the_same_density(backend: str) -> None:
    """The contract path's ``log_prob`` is one number, whatever computes it.

    Bit-for-bit is not asserted — different libraries reduce in different orders
    — but the tolerance is tight enough that any real disagreement in the
    physics, the kernel or the solve would fail it.
    """
    pytest.importorskip(backend)
    data = generators.generate("strong_smooth", size=200)
    for kind, solver in (("standard", "quasisep"), ("flexible", "quasisep"), ("flexible", "dense")):
        reference = study.build_problem(data, backend="reference", likelihood=kind, solver=solver)
        native = study.build_problem(data, backend=backend, likelihood=kind, solver=solver)
        expected = float(reference.log_prob(reference.reference_values))
        got = float(native.log_prob(native.reference_values))
        assert got == pytest.approx(expected, rel=1e-10), f"{kind}/{solver}"
