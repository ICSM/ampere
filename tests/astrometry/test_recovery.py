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
import math
from typing import Any

import numpy as np
import pytest

from examples.astrometry.astrometry import QUALIFIED_TRUTH, build_problem, fit, period_modes

needs_torch = pytest.mark.skipif(
    importlib.util.find_spec("torch") is None, reason="needs ampere[torch]"
)
needs_jax = pytest.mark.skipif(importlib.util.find_spec("jax") is None, reason="needs ampere[jax]")
needs_nautilus = pytest.mark.skipif(
    importlib.util.find_spec("nautilus") is None, reason="needs ampere[nautilus]"
)
needs_ultranest = pytest.mark.skipif(
    importlib.util.find_spec("ultranest") is None, reason="needs ampere[ultranest]"
)
#: The `study` marker's convention (W5.26), applied locally because this
#: directory has no conftest hook for it: a dev-only budget row runs where
#: neither modern backend is installed, so the dev gate leg pays its two and
#: a half minutes once and the torch, jax and sbi legs do not repeat it.
dev_only = pytest.mark.skipif(
    importlib.util.find_spec("torch") is not None or importlib.util.find_spec("jax") is not None,
    reason="a dev-only budget row: the torch, jax and sbi legs skip it (the `study` convention)",
)


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


# ---------------------------------------------------------------------------
# The wide-prior arm: nested sampling resolves the aliasing (W5.15)
# ---------------------------------------------------------------------------


def _assert_resolves_the_true_mode(run: Any) -> None:
    """The pinned claim: the true period is inside the largest-mass mode."""
    modes = period_modes(run)
    assert modes, "period_modes returned no mode at all"
    largest = max(modes, key=lambda mode: mode["mass_fraction"])
    assert largest["lower"] <= 400.0 <= largest["upper"], (
        f"the true period (400) is not inside the largest mode's interval "
        f"[{largest['lower']}, {largest['upper']}]"
    )
    assert all(largest["mass_fraction"] >= mode["mass_fraction"] for mode in modes)
    assert math.isfinite(run.attrs["ampere_log_evidence"])
    assert run.attrs["ampere_log_evidence_err"] > 0.0


class TestNestedSamplingResolvesThePeriodModes:
    """W5.15: the wide-prior posterior's mode structure, pinned rather than narrated.

    Measured on this branch: dynesty at its default live points (175) with
    the study's ``DYNESTY_SAMPLE`` ("rslice"; ``examples.astrometry.astrometry``
    says why the engine's own ``sample="auto"`` did not converge here in 38
    minutes) resolves the true mode in about four minutes -- over the
    per-PR bar, so that row is ``astrometry_full``; nautilus and ultranest
    (W5.14) resolve the same claim in one to two and a half minutes each at
    their library defaults and carry no marker, skipping by ``find_spec``
    exactly as ``tests/inference/test_nested.py`` does where the package is
    absent. Every row is the same call ``python -m examples.astrometry
    --wide-prior --engine <name>`` makes.
    """

    @pytest.mark.astrometry_full
    def test_dynesty_resolves_the_true_mode(self) -> None:
        """The CLI's own call, at dynesty's default live points: 241 s measured."""
        problem = build_problem("reference", wide_prior=True)
        run = fit(problem, backend="reference", engine="dynesty")
        _assert_resolves_the_true_mode(run)

    @dev_only
    def test_dynesty_resolves_the_true_mode_at_a_per_pr_budget(self) -> None:
        """The per-PR pin: 100 live points, 142 s measured, the same mode structure.

        Below ``default_live_points``'s ``25 (n_dim + 1)`` floor the evidence
        is less trustworthy (``+91.63 +- 0.67`` here against ``+92.25 +-
        0.55`` at 175), which is why the CLI does not run at this budget; the
        mode structure -- one dominant mode, the truth inside it -- is what
        this row pins, and it is unchanged.
        """
        problem = build_problem("reference", wide_prior=True)
        run = fit(problem, backend="reference", engine="dynesty", live_points=100)
        _assert_resolves_the_true_mode(run)

    @needs_nautilus
    def test_nautilus_resolves_the_true_mode(self) -> None:
        problem = build_problem("reference", wide_prior=True)
        run = fit(problem, backend="reference", engine="nautilus")
        _assert_resolves_the_true_mode(run)

    @needs_ultranest
    def test_ultranest_resolves_the_true_mode(self) -> None:
        problem = build_problem("reference", wide_prior=True)
        run = fit(problem, backend="reference", engine="ultranest")
        _assert_resolves_the_true_mode(run)
