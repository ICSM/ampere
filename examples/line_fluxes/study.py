"""Does the flexible likelihood survive a misspecified excitation ladder? (W7.9).

Three arms, as :mod:`examples.m2_misspecification` has them, fitted to one
synthetic line table (:mod:`.generators`) through one instrument,
``LineIntegration.from_observed(table, window=300 km/s)``:

``"correct"``
    The two-component ladder the truth was drawn from, independent noise — the
    control.
``"misspecified"``
    The one-component ladder, independent noise: confidently wrong, because
    one Boltzmann term cannot follow the truth's curvature along the ladder.
``"flexible"``
    The same one-component ladder under the flexible likelihood — a Matérn-3/2
    over **rest wavelength** with ``DenseGP`` (ten points) — which absorbs the
    residual the missing component leaves, correlated across neighbouring
    transitions.

The parameters compared are the cool component's ``temperature`` and
``column``, which every arm fits. ``run`` dispatches on the problem's backend:
emcee on the reference path, NUTS through the realisation on torch and jax.
"""

from __future__ import annotations

import dataclasses
from typing import Any

import numpy as np
import scipy.stats as st

from ampere.core import Dataset, DatasetCollection, FittingProblem, GaussianFamily, Likelihood

from . import generators as gen
from .model import FLUX_SCALE

__all__ = [
    "ARMS",
    "BACKENDS",
    "CI_BUDGET",
    "DOC_BUDGET",
    "FLEXIBLE_MAX_BIAS_WIDTHS",
    "MISSPECIFIED_MIN_BIAS_WIDTHS",
    "NUTS_BUDGET",
    "PARAMETERS",
    "EmceeBudget",
    "NutsBudget",
    "bias_widths",
    "build_problem",
    "instrument_for",
    "localisation_peak",
    "main",
    "model_for",
    "run",
    "truth_model",
]

ARMS: tuple[str, ...] = ("correct", "misspecified", "flexible")
BACKENDS: tuple[str, ...] = ("reference", "torch", "jax")
#: The parameters every arm fits, and against which the arms are compared.
PARAMETERS: tuple[str, ...] = ("model.temperature", "model.column")


@dataclasses.dataclass(frozen=True)
class EmceeBudget:
    """Walkers, steps and burn-in for an ensemble run."""

    walkers: int
    steps: int
    burn_in: int


@dataclasses.dataclass(frozen=True)
class NutsBudget:
    """Draws, warm-up and chains for a NUTS run."""

    draws: int
    warmup: int
    chains: int


#: The per-PR budget, ``tests/examples/test_line_fluxes.py``'s.
CI_BUDGET = EmceeBudget(walkers=16, steps=600, burn_in=300)
#: The documentation / ``__main__`` budget.
DOC_BUDGET = EmceeBudget(walkers=32, steps=2000, burn_in=800)
#: The gradient path's, for the torch and jax rows (marked ``heavy``).
NUTS_BUDGET = NutsBudget(draws=300, warmup=300, chains=2)

# ---------------------------------------------------------------------------
# The thresholds the tests assert, each measured over three seeds at
# CI_BUDGET and set a third beyond the worst (the study's margin convention).
# ---------------------------------------------------------------------------

#: The flexible arm's ``temperature`` and ``column`` must stay inside this many
#: posterior widths of the cool component's truth (the control is held to it
#: too). Measured at seeds 20261014, 1 and 2: the flexible arm's worst is 1.93,
#: 1.84 and 1.87 widths, the control's 1.85, 0.18 and 0.40; 1.93 x 4/3 = 2.6.
FLEXIBLE_MAX_BIAS_WIDTHS = 2.6
#: The misspecified, independent-noise arm must miss by at least this many on
#: at least one of them. Measured: 16.27, 16.83 and 15.99 widths (always on
#: ``temperature``); 15.99 x 3/4 = 12.
MISSPECIFIED_MIN_BIAS_WIDTHS = 12.0

#: The prior on the cool component, every arm.
TEMPERATURE_PRIOR = st.uniform(200.0, 2800.0)
COLUMN_PRIOR = st.loguniform(300.0, 30000.0)
#: The flexible arm's kernel: an amplitude prior in units of the median
#: tabulated flux (the residual a missing component leaves is a fraction of
#: the flux, not of its error bar), and a length scale in micron spanning a
#: few neighbouring lines to the whole ladder.
GP_AMPLITUDE_FRACTION = 0.5
GP_LENGTH_PRIOR = (1.0, 30.0)


def _backend_module(backend: str) -> Any:
    """Where this backend's ``LineIntegration`` and noise models live."""
    if backend == "reference":
        from ampere.backends import reference as module

        return module
    if backend == "torch":
        from ampere.backends import torch as module  # type: ignore[assignment]

        return module
    if backend == "jax":
        from ampere.backends import jax as module  # type: ignore[assignment]

        module.configure_x64()
        return module
    raise ValueError(f"unknown backend {backend!r}; the three are {list(BACKENDS)}.")


def _noise_module(backend: str) -> Any:
    """Where this backend's ``IndependentNoise``/``GaussianProcessNoise`` live."""
    if backend == "reference":
        import ampere.core as module

        return module
    return _backend_module(backend)


def _model_class(backend: str) -> Any:
    if backend == "torch":
        from .model_torch import TorchEmissionLadder

        return TorchEmissionLadder
    if backend == "jax":
        _backend_module("jax")
        from .model_jax import JaxEmissionLadder

        return JaxEmissionLadder
    from .model import EmissionLadder

    return EmissionLadder


def model_for(backend: str, arm: str, *, fitted: bool = True) -> Any:
    """The ladder *arm* fits with on *backend*: two components for ``"correct"``, else one.

    ``fitted=False`` fixes every parameter at the truth (the generator's model).
    """
    if arm not in ARMS:
        raise ValueError(f"unknown arm {arm!r}; the three are {list(ARMS)}.")
    truth = gen.TRUTH
    kwargs: dict[str, Any] = {
        "temperature": TEMPERATURE_PRIOR if fitted else truth["temperature"],
        "column": COLUMN_PRIOR if fitted else truth["column"],
    }
    if arm == "correct":
        # Carried at its true, fixed value, as the image study's control
        # carries its background: the control asks whether a fit *told* the
        # hot component is there recovers the cool one honestly. Fitting it
        # too would spend the budget on the two components' near-degeneracy
        # over a ladder this short — a question the study is not asking.
        kwargs["hot_temperature"] = truth["hot_temperature"]
        kwargs["hot_column"] = truth["hot_column"]
    return _model_class(backend)(
        np.linspace(6.0, 27.0, 50),
        rest=gen.REST,
        energy=gen.ENERGY,
        strength=gen.STRENGTH,
        width=gen.WIDTH,
        **kwargs,
    )


def truth_model(backend: str = "reference") -> Any:
    """The two-component ladder at the truth, every parameter fixed."""
    return model_for(backend, "correct", fitted=False)


def instrument_for(backend: str, observed: Any = None) -> Any:
    """``LineIntegration`` over the line list, as the one-step instrument.

    From the observed table when there is one (the supported route,
    ``from_observed``); from the line list itself to generate one.
    """
    from ampere.core import Instrument

    step_class = _backend_module(backend).LineIntegration
    if observed is None:
        step = step_class(list(gen.LABELS), gen.REST, window=gen.WINDOW)
    else:
        step = step_class.from_observed(observed, window=gen.WINDOW)
    return Instrument([step], channel="lines", label="lines")


def _noise_for(backend: str, arm: str, *, scale: float, solver: Any = None) -> Any:
    module = _noise_module(backend)
    if arm != "flexible":
        return module.IndependentNoise()
    kernel = module.Matern32(
        st.halfnorm(scale=GP_AMPLITUDE_FRACTION * scale),
        st.uniform(GP_LENGTH_PRIOR[0], GP_LENGTH_PRIOR[1] - GP_LENGTH_PRIOR[0]),
        axes=("spectral_axis",),
    )
    return module.GaussianProcessNoise(kernel, solver or module.DenseGP())


def build_problem(
    backend: str, arm: str, *, seed: int = gen.SEED, solver: Any = None
) -> FittingProblem:
    """One arm's problem on the synthetic table drawn at *seed*."""
    table = gen.synthetic_table(truth_model(), instrument_for("reference"), seed=seed)
    instrument = instrument_for(backend, table)
    scale = float(np.median(np.asarray(table.flux, dtype=float)))
    noise = _noise_for(backend, arm, scale=scale, solver=solver)
    datasets = DatasetCollection(
        {
            "lines": Dataset(
                table,
                instrument,
                likelihood=Likelihood(GaussianFamily(), noise),
                label="lines",
            )
        }
    )
    return FittingProblem(model_for(backend, arm), datasets, seed=seed)


def run(
    problem: FittingProblem,
    budget: EmceeBudget | NutsBudget = CI_BUDGET,
    *,
    progress: bool = False,
) -> Any:
    """Sample *problem* with the engine its backend can feed (emcee, or NUTS natively)."""
    from ampere.inference import EmceeEngine, NUTSEngine

    if problem.backend == "reference":
        chosen = budget if isinstance(budget, EmceeBudget) else CI_BUDGET
        engine = EmceeEngine(problem, walkers=chosen.walkers)
        return engine.run(chosen.steps, burn_in=chosen.burn_in, progress=progress)
    nuts = budget if isinstance(budget, NutsBudget) else NUTS_BUDGET
    return NUTSEngine(problem).run(
        nuts.draws, warmup=nuts.warmup, chains=nuts.chains, progress=progress
    )


def bias_widths(tree: Any) -> dict[str, float]:
    """``|posterior mean - truth| / posterior sd`` for the cool component's two parameters."""
    truth = {"model.temperature": gen.TRUTH["temperature"], "model.column": gen.TRUTH["column"]}
    out: dict[str, float] = {}
    for name in PARAMETERS:
        draws = np.asarray(tree.posterior[name], dtype=float).ravel()
        out[name] = float(abs(draws.mean() - truth[name]) / draws.std(ddof=1))
    return out


def localisation_peak(tree: Any, problem: FittingProblem) -> str:
    """The transition at which the flexible arm's conditioned GP mean is largest in magnitude."""
    from ampere.results import gp_localisation

    group = gp_localisation(tree, problem)["gp_localisation"]
    names = [name for name in group.data_vars if str(name).startswith("lines")]
    mean = np.asarray(group[names[0]], dtype=float)
    mean = mean.reshape(-1, mean.shape[-1]).mean(axis=0)
    observed = problem.datasets["lines"].observed
    return str(observed.transitions[int(np.argmax(np.abs(mean)))])


def main(argv: list[str] | None = None) -> int:
    """Fit the three arms and print each one's bias in posterior widths."""
    import argparse

    parser = argparse.ArgumentParser(prog="python -m examples.line_fluxes")
    parser.add_argument("--backend", choices=BACKENDS, default="reference")
    parser.add_argument("--seed", type=int, default=gen.SEED)
    parser.add_argument("--quick", action="store_true", help="the per-PR budget")
    arguments = parser.parse_args(argv)
    budget: EmceeBudget | NutsBudget = CI_BUDGET if arguments.quick else DOC_BUDGET
    if arguments.backend != "reference":
        budget = NUTS_BUDGET
    print(f"line fluxes: {len(gen.LABELS)} transitions, seed {arguments.seed}")
    print(
        f"truth (cool component): T = {gen.TRUTH['temperature']:.0f} K, "
        f"N = {gen.TRUTH['column']:.0f} x {FLUX_SCALE:.0e} W m^-2"
    )
    for arm in ARMS:
        problem = build_problem(arguments.backend, arm, seed=arguments.seed)
        tree = run(problem, budget)
        widths = bias_widths(tree)
        summary = ", ".join(
            f"{name.split('.')[-1]} {float(np.mean(tree.posterior[name])):.4g} "
            f"({widths[name]:.2f} widths)"
            for name in PARAMETERS
        )
        print(f"  {arm:13s} {summary}")
        if arm == "flexible":
            print(f"  {'':13s} GP localisation peaks at {localisation_peak(tree, problem)}")
    return 0
