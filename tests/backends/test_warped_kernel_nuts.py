"""NUTS over a warped kernel's knots, on torch and jax (W5.7).

``tests/conformance/test_likelihoods.py`` already holds :class:`WarpedKernel`
to the §6 contract on every registered fixture — the warp is what the explicit
warp gives, the O(N) path reproduces the dense one, the identity warp is the
base kernel bit for bit. What the neutral battery has no vocabulary for is the
reason the knots were declared as ordinary
:class:`~ampere.core.parameter.Parameter`\\ s in the first place: **a gradient
flows to them**, so a warp is *fitted* rather than chosen, on both
differentiable backends, with no bespoke sampler support at all.

Three rows, all at small budgets and fixed seeds, so this file belongs in the
per-PR gate beside the other backend suites rather than on a nightly:

1. the **non-centred default** — ``u_k = s · z_k`` with ``z_k ~ Normal(0, 1)``
   — sampled through :class:`~ampere.inference.NUTSEngine`, with every warp
   variable and every base hyperparameter present in the posterior and finite;
2. the **centred form**, declared with
   :class:`~ampere.core.parameter.HierarchicalPrior`, which is what the plan
   names and what a lowering rule has to be able to carry: the same fit, the
   same finiteness, so the shrinkage declaration is not a form that only
   type-checks;
3. the **gradient itself**, taken directly through each backend's autodiff at
   a fixed point, and compared *between* the backends — the two agree to
   float64 because both differentiate the same generators through the same
   celerite recursion, so a disagreement would mean one of the two threadings
   of ``Kernel.warped_coordinate`` is wrong.

Parametrised over whichever of torch and jax this environment has, in
``tests/backends/test_native_astrometry.py``'s shape.
"""

from __future__ import annotations

import dataclasses
import importlib
from typing import Any

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.core import (
    Dataset,
    DatasetCollection,
    FittingProblem,
    GaussianFamily,
    HierarchicalPrior,
    Instrument,
    Likelihood,
    Spectrum,
    WarpedKernel,
)

SEED = 20260916
#: Knot locations spanning the grid. Four for the input warp (three free
#: increments) and three for the amplitude warp: "few knots" is the
#: degrees-of-freedom guard's first clause, not a budget compromise.
INPUT_KNOTS = (1.0, 4.0, 7.0, 10.0)
AMPLITUDE_KNOTS = (1.0, 5.5, 10.0)
GRID = np.linspace(1.0, 10.0, 48)


@dataclasses.dataclass(frozen=True)
class Kit:
    """One backend's pieces, by name, so no row here names a library."""

    name: str
    module: Any

    def gradient(self, function: Any, at: np.ndarray) -> np.ndarray:
        """``d function / d x`` at *at*, through this backend's own autodiff."""
        if self.name == "torch":
            torch = importlib.import_module("torch")
            argument = torch.tensor(at, dtype=torch.float64, requires_grad=True)
            function(argument).backward()
            assert argument.grad is not None
            return np.asarray(argument.grad.detach().cpu().numpy(), dtype=float)
        jax = importlib.import_module("jax")
        jnp = importlib.import_module("jax.numpy")
        return np.asarray(jax.grad(function)(jnp.asarray(at, dtype=jnp.float64)), dtype=float)


def _installed_kits() -> list[Kit]:
    found: list[Kit] = []
    for name in ("jax", "torch"):
        try:
            module = importlib.import_module(f"ampere.backends.{name}")
        except ImportError:  # the extra is not installed in this environment
            continue
        if name == "jax":
            # ``lowering.md`` §10.2(a): the application turns x64 on, not
            # ampere. A test suite is an application.
            module.configure_x64()
        found.append(Kit(name, module))
    return found


KITS = _installed_kits()

pytestmark = pytest.mark.skipif(
    not KITS,
    reason="no differentiable backend installed; a fitted warp needs one",
)


@pytest.fixture(scope="module", params=[kit.name for kit in KITS])
def kit(request: Any) -> Kit:
    return next(found for found in KITS if found.name == request.param)


# ---------------------------------------------------------------------------
# The data: a smooth model with a genuinely non-stationary residual
# ---------------------------------------------------------------------------


def observations() -> Spectrum:
    """A line whose residual wiggles fast at the blue end and slowly at the red.

    Deliberately the situation a stationary Matérn cannot express: one
    correlation length cannot be both. The warp is not asked to *recover* this
    to a tolerance here — that is W5.8's study, with margins — only to be
    sampled over without the geometry falling apart.
    """
    rng = np.random.default_rng(SEED)
    frequency = 6.0 / GRID
    wiggle = 0.25 * np.sin(2.0 * np.pi * frequency * GRID)
    sigma = 0.05
    values = 1.5 * GRID**-0.3 + wiggle + rng.normal(0.0, sigma, GRID.size)
    return Spectrum(
        GRID * u.micron,
        values * u.Jy,
        uncertainty=np.full(GRID.size, sigma) * u.Jy,
    )


DATA = observations()


def warped_kernel(kit: Kit, *, non_centred: bool = True) -> WarpedKernel:
    """A fitted Matérn-3/2 under a fitted warp, all knots free.

    The centred form is spelled out rather than taken from a keyword so that
    the row which uses it is visibly declaring a
    :class:`~ampere.core.parameter.HierarchicalPrior`, references and all.
    """
    base = kit.module.Matern32(st.loguniform(1e-2, 1.0), st.loguniform(0.3, 6.0))
    if non_centred:
        return WarpedKernel(
            base,
            input_warp=INPUT_KNOTS,
            amplitude_warp=AMPLITUDE_KNOTS,
            non_centred=True,
        )
    return WarpedKernel(
        base,
        input_warp=INPUT_KNOTS,
        amplitude_warp=AMPLITUDE_KNOTS,
        increments=HierarchicalPrior("norm", {"scale": "input_warp.scale"}, kwds={"loc": 0.0}),
        levels=HierarchicalPrior("norm", {"scale": "amplitude_warp.scale"}, kwds={"loc": 0.0}),
        non_centred=False,
    )


def warped_problem(kit: Kit, *, non_centred: bool = True) -> FittingProblem:
    """One line, one channel, and the flexible likelihood under a fitted warp."""
    model = kit.module.PowerLaw(GRID, norm=st.lognorm(0.3, scale=1.5), index=st.norm(-0.3, 0.2))
    # **This backend's** noise model, not ``ampere.core``'s: since W2.13 a
    # noise model carries the four capability flags, and the core's declares
    # "reference", which would make this a two-backend problem and be refused
    # at composition. That refusal names ``WarpedKernel`` as *torch*, which is
    # the adoption rule working.
    noise = kit.module.GaussianProcessNoise(
        warped_kernel(kit, non_centred=non_centred), kit.module.QuasisepGP()
    )
    return FittingProblem(
        model,
        DatasetCollection(
            {"line": Dataset(DATA, Instrument([]), Likelihood(GaussianFamily(), noise))}
        ),
        seed=SEED,
    )


# ---------------------------------------------------------------------------
# Rows
# ---------------------------------------------------------------------------


class TestNUTSOverTheKnots:
    """The knots are ordinary parameters, so a gradient sampler gets them free."""

    @pytest.mark.parametrize("non_centred", [True, False], ids=["non_centred", "centred"])
    def test_nuts_samples_the_warp_and_the_hyperparameters(
        self, kit: Kit, non_centred: bool
    ) -> None:
        from ampere.inference import NUTSEngine

        problem = warped_problem(kit, non_centred=non_centred)
        free = set(problem.mapping.merged.free_names)
        # The declaration is what this row is really about: every knot
        # variable and every shrinkage scale is a *sampled* coordinate, not a
        # setting, and the base kernel's own hyperparameters are still there.
        wanted = {
            "line.likelihood.input_warp.scale",
            "line.likelihood.input_warp.increment0",
            "line.likelihood.input_warp.increment2",
            "line.likelihood.amplitude_warp.scale",
            "line.likelihood.amplitude_warp.level0",
            "line.likelihood.amplitude_warp.level2",
            "line.likelihood.base.amplitude",
            "line.likelihood.base.length_scale",
        }
        assert wanted <= free, sorted(free)

        run = NUTSEngine(problem).run(30, warmup=50, chains=1, progress=False)
        posterior = run["posterior"]
        for name in sorted(wanted):
            values = np.asarray(posterior[name].values, dtype=float)
            assert values.size == 30
            assert np.all(np.isfinite(values)), name
        # A sampled warp must actually move: a chain pinned at its initial
        # point would pass every finiteness check above and mean nothing.
        moved = np.asarray(posterior["line.likelihood.input_warp.increment0"].values, dtype=float)
        assert float(np.std(moved)) > 0.0


class TestTheGradientReachesTheKnots:
    """Counted rather than trusted: ``d log p / d u_k`` is non-zero and backend-agnostic."""

    #: A fixed point at which both backends are asked for the same derivative.
    KNOT_VALUES = np.array([0.6, -0.4, 0.3, 0.2, -0.3, 0.4])

    #: What torch and jax agree to. Both differentiate the same generators
    #: through the same compiled celerite recursion in float64, so the only
    #: difference is accumulation order.
    CROSS_BACKEND = 1e-8

    def _gradient(self, kit: Kit) -> np.ndarray:
        module = kit.module
        base = module.Matern32(0.3, 2.0)
        kernel = WarpedKernel(
            base,
            input_warp=INPUT_KNOTS,
            amplitude_warp=AMPLITUDE_KNOTS,
            increments=1.0,
            levels=1.0,
            input_scale=1.0,
            amplitude_scale=1.0,
        )
        solver = module.QuasisepGP()
        points = GRID[:, None]
        residual = np.asarray(DATA.values, dtype=float) - 1.5 * GRID**-0.3
        variance = np.full(GRID.size, 0.05**2)

        def density(free: Any) -> Any:
            values = dict(kernel.resolve(None))
            for index in range(3):
                values[f"input_warp.increment{index}"] = free[index]
                values[f"amplitude_warp.level{index}"] = free[3 + index]
            if kit.name == "torch":
                return solver.log_marginal_likelihood_tensor(
                    kernel, points, residual, variance, values
                )
            return solver.log_marginal_likelihood_jax(kernel, points, residual, variance, values)

        return kit.gradient(density, self.KNOT_VALUES)

    def test_every_knot_variable_receives_a_gradient(self, kit: Kit) -> None:
        gradient = self._gradient(kit)
        assert gradient.shape == self.KNOT_VALUES.shape
        assert np.all(np.isfinite(gradient))
        # Not merely finite: a detached graph gives exact zeros, which is the
        # failure W2.13 found in the torch kernels and the one this row exists
        # to make impossible for the warp.
        assert np.all(np.abs(gradient) > 1e-6)

    def test_the_two_backends_agree_on_it(self, kit: Kit) -> None:
        """Both columns against one recorded reference, produced on **jax**.

        Only one differentiable backend is installed per pixi environment
        (``torch`` and ``jax`` are separate extras), so the two cannot be
        compared live in one process. Recording jax's value and holding torch
        to it makes the torch column a genuine cross-backend row; the jax
        column is then the regression pin that keeps the recorded value
        honest. Measured, the two agree to 6e-14 — accumulation order in the
        same float64 recursion — and the tolerance is set well inside what a
        wrongly threaded warp would cost (the same derivative with the warp
        omitted from the recursion differs in the first digit).
        """
        expected = np.array(
            [
                -3.299189679476724,
                -2.1256652990843388,
                -3.2399221450942477,
                -4.628824456726996,
                -5.737920860919208,
                -4.2409871676306645,
            ]
        )
        assert self._gradient(kit) == pytest.approx(expected, abs=self.CROSS_BACKEND)
