"""``examples/portable_model``: one model source, fitted on three backends (W7.13).

The reference row fits under emcee at a tiny budget. The jax row runs the
guide's NUTS cell exactly -- :func:`examples.portable_model.fit.fit_nuts` at
its own :data:`~examples.portable_model.fit.NUTS_BUDGET` -- so CI's jax leg
checks the result ``docs/source/notebooks/portable_model.ipynb`` quotes. The
torch row runs the same function at a smaller budget. Every row also checks
that the twins it can import agree on a flux with the one numpy source.
"""

from __future__ import annotations

import importlib.util
import warnings

import numpy as np
import pytest

from examples.portable_model import LinearModel
from examples.portable_model.fit import (
    GRID,
    TRUTH,
    backend_module,
    build_problem,
    fit_emcee,
    fit_nuts,
    model_class,
)

needs_torch = pytest.mark.skipif(
    importlib.util.find_spec("torch") is None, reason="needs ampere[torch]"
)
needs_jax = pytest.mark.skipif(importlib.util.find_spec("jax") is None, reason="needs ampere[jax]")

POINT = {"slope": -0.7, "intercept": 4.25}


def _agreement(backend: str) -> None:
    """*backend*'s twin and the numpy source give the same flux, on both surfaces."""
    backend_module(backend)  # jax: float64 on before any jax work
    reference = LinearModel(GRID)(**POINT)["sed"].values
    twin = model_class(backend)(GRID)
    np.testing.assert_allclose(twin(**POINT)["sed"].values, reference, rtol=1e-12, atol=0.0)
    native = np.asarray(twin.native_flux("sed", POINT))
    np.testing.assert_allclose(native, reference, rtol=1e-12, atol=0.0)
    np.testing.assert_allclose(reference, POINT["slope"] * GRID + POINT["intercept"], rtol=1e-12)


def _recovers(run: object) -> None:
    for name, truth in TRUTH.items():
        draws = np.asarray(run["posterior"][f"model.{name}"])  # type: ignore[index]
        assert abs(float(draws.mean()) - truth) < 5.0 * float(draws.std()), name


class TestTheReferenceTwin:
    def test_the_package_imports_without_a_native_backend(self) -> None:
        assert LinearModel.BACKEND == "reference"
        assert LinearModel.DIFFERENTIABLE is False

    def test_a_tiny_emcee_fit(self) -> None:
        problem = build_problem("reference")
        assert problem.backend == "reference"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            run = fit_emcee(problem, walkers=8, steps=40, burn_in=10)
        assert run.attrs["ampere_engine"] == "emcee"
        assert {"model.slope", "model.intercept"} <= set(run["posterior"].dataset.data_vars)
        _agreement("reference")


@needs_jax
class TestJax:
    def test_the_guide_nuts_cell(self) -> None:
        _agreement("jax")
        problem = build_problem("jax")
        assert problem.backend == "jax" and problem.differentiable is True
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            run = fit_nuts(problem)
        assert run.attrs["ampere_backend"] == "jax"
        assert run.attrs["ampere_engine"] == "nuts"
        _recovers(run)


@needs_torch
class TestTorch:
    def test_a_small_nuts_fit(self) -> None:
        import torch

        torch.set_num_threads(4)
        _agreement("torch")
        problem = build_problem("torch")
        assert problem.backend == "torch" and problem.differentiable is True
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            run = fit_nuts(problem, draws=150, warmup=150, chains=1)
        assert run.attrs["ampere_backend"] == "torch"
        _recovers(run)
