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

    Measured on this branch (see the branch report): dynesty's ``bound="multi"``
    ellipsoidal decomposition needs on the order of half an hour on this
    problem's likelihood surface at dynesty's default live points, and it is
    the bound's own machinery that is slow, not the wide prior's
    multi-modality specifically -- the *informed* prior (``norm(400, 30)``,
    unimodal) is measured no faster. No live-point/``dlogz`` combination
    tried resolved the true mode reliably in under four minutes, so this
    row is ``astrometry_full`` rather than a per-PR dev row (the item's own
    escape valve for exactly this finding). Nautilus and ultranest (W5.14)
    resolve the same claim in one to two minutes each at their own library
    defaults and carry no such marker; they skip by ``find_spec`` exactly as
    ``tests/inference/test_nested.py`` skips them where the package is absent.
    """

    @pytest.mark.astrometry_full
    def test_dynesty_resolves_the_true_mode(self) -> None:
        problem = build_problem("reference", wide_prior=True)
        run = fit(problem, backend="reference", engine="dynesty", dlogz=200.0)
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
