"""The item's pinned claim: an emcee and a NUTS fit recover an injected orbit.

Measured at the pinned seed (``examples.astrometry.generators.SEED``) and the
budgets below (a per-PR budget, well short of ``python -m
examples.astrometry``'s own default): every one of the six orbital
parameters (``pmra``, ``pmdec``, ``period``, ``phase``, ``amp_ra``,
``amp_dec``) lands inside its central 95 % credible interval, on both the
reference backend under emcee and the torch/jax backends under NUTS. See
``docs/source/astrometry.rst``'s closing section for why the period prior is
informative rather than a wide search: a reflex orbit is periodic, unlike
every model the interferometry template fits, and a period search wide
enough to reach past the observed epochs' own baseline is a genuinely
multi-modal problem an MCMC chain can alias onto and never leave.
"""

from __future__ import annotations

import importlib.util
from typing import Any

import numpy as np
import pytest

from examples.astrometry.astrometry import QUALIFIED_TRUTH, build_problem, fit

needs_torch = pytest.mark.skipif(
    importlib.util.find_spec("torch") is None, reason="needs ampere[torch]"
)
needs_jax = pytest.mark.skipif(importlib.util.find_spec("jax") is None, reason="needs ampere[jax]")


def _covered(run: Any, *, level: float = 0.95) -> dict[str, bool]:
    posterior = run["posterior"].dataset
    tail = (1.0 - level) / 2.0 * 100.0
    covered: dict[str, bool] = {}
    for name, truth in QUALIFIED_TRUTH.items():
        draws = np.asarray(posterior[name], dtype=float).ravel()
        lower, upper = np.percentile(draws, [tail, 100.0 - tail])
        covered[name] = bool(lower <= truth <= upper)
    return covered


class TestEmceeRecoversTheOrbit:
    """The reference backend, gradient-free, at a per-PR budget."""

    @pytest.fixture(scope="class")
    @classmethod
    def run(cls) -> Any:
        problem = build_problem("reference")
        return fit(problem, backend="reference", walkers=24, steps=600, burn_in=200)

    @pytest.mark.parametrize("name", list(QUALIFIED_TRUTH))
    def test_every_parameter_is_inside_the_central_95_percent(self, run: Any, name: str) -> None:
        assert _covered(run)[name], f"{name} missed its central 95% interval"


class TestEmceeRecoversTheOrbitUnderHeteroscedasticJointNoise:
    """W5.24: the joint arm on channels with their own per-epoch sigmas, same floor.

    ``build_problem(joint=True, heteroscedastic=True)`` injects the correlated
    centroiding systematic *and* gives each channel its own error bar per
    epoch, so the joint group is scored on the dense route; the coupling's
    three parameters are fitted beside the orbit (nine dimensions, hence the
    walker count). Every orbital parameter must land inside its central 95 %
    interval, the same claim the rows above make.
    """

    @pytest.fixture(scope="class")
    @classmethod
    def run(cls) -> Any:
        problem = build_problem("reference", joint=True, heteroscedastic=True)
        return fit(problem, backend="reference", walkers=24, steps=600, burn_in=200)

    @pytest.mark.parametrize("name", list(QUALIFIED_TRUTH))
    def test_every_parameter_is_inside_the_central_95_percent(self, run: Any, name: str) -> None:
        assert _covered(run)[name], f"{name} missed its central 95% interval"


@needs_torch
class TestNutsOnTorchRecoversTheOrbit:
    """The torch backend, at a per-PR budget far short of a global period search."""

    @pytest.fixture(scope="class")
    @classmethod
    def run(cls) -> Any:
        problem = build_problem("torch")
        return fit(problem, backend="torch", draws=150, warmup=150, chains=2)

    @pytest.mark.parametrize("name", list(QUALIFIED_TRUTH))
    def test_every_parameter_is_inside_the_central_95_percent(self, run: Any, name: str) -> None:
        assert _covered(run)[name], f"{name} missed its central 95% interval"


@needs_jax
class TestNutsOnJaxRecoversTheOrbit:
    """The jax backend, same shape as the torch row above."""

    @pytest.fixture(scope="class")
    @classmethod
    def run(cls) -> Any:
        from ampere.backends import jax as _jax

        _jax.configure_x64()
        problem = build_problem("jax")
        return fit(problem, backend="jax", draws=150, warmup=150, chains=2)

    @pytest.mark.parametrize("name", list(QUALIFIED_TRUTH))
    def test_every_parameter_is_inside_the_central_95_percent(self, run: Any, name: str) -> None:
        assert _covered(run)[name], f"{name} missed its central 95% interval"
