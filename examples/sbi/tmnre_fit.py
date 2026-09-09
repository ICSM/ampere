#!/usr/bin/env python
"""Truncated marginal ratio estimation: rounds, a shrinking box, marginals.

W3.4's ``SBIEngine(method="tmnre", ...)`` end to end on a small hand-written
problem, printed so that every claim the method makes is visible in the output
rather than asserted in a docstring:

* the **truncation box** after each round, in the user's own coordinates, with
  the fraction of untruncated prior draws it accepts — the box's own prior
  mass, and the honest measure of how much the round's budget was concentrated;
* the **marginals group**, one estimated 1-D posterior per parameter over the
  final box, printed as a credible interval and as a sparkline so that a
  reader can see the shape without a plotting backend;
* the **posterior**, which is an ordinary posterior: a joint ratio estimator
  trained across the rounds and multiplied by the final truncated prior,
  sampled by rejection, and scored on the numpy contract path like every other
  ampere run's draws.

Run it::

    python examples/sbi/tmnre_fit.py                         # ~6 min: see below
    python examples/sbi/tmnre_fit.py --sample-with mcmc      # the same fit in ~45 s
    python examples/sbi/tmnre_fit.py --rounds 4 --budget 1500 --marginals 2

Requires the ``sbi`` extra (``pixi run -e sbi python examples/sbi/tmnre_fit.py``).

**Where the time goes, and what to do about it.** Almost all of a default run
is the final *rejection* sampling, not the training or the simulating: the
sampler proposes from the truncated prior and accepts through the ratio, and
*both* of those rejections get harder as the box tightens — which is to say,
as the method works. On the fit below the third-round box holds about 1 % of
the prior's mass and the ratio accepts about 1 % of what the box proposes, so
every stored draw costs of order 10⁴ draws from the untruncated prior; ``sbi``
says so itself, in a warning naming the remedy. Measured here: 353.6 s with
``--sample-with rejection`` (the default, whose draws are i.i.d.) against
43.8 s with ``--sample-with mcmc`` for the same three rounds. The run records
which sampler produced its draws either way.

What to look at
---------------
**The boxes nest, and they contain the truth.** The problem below is generated
at a known θ, printed beside every box. A box that lost the truth would be the
method failing in the one way that matters — no later round can put back mass
the truncation threw away — so the script says explicitly whether each one
still holds it.

**The estimate is not amortised, and the run says so.** ``ampere_sbi_amortised``
is ``0`` for every TMNRE run: the box is chosen at *this* observation, so the
trained estimator must not be re-conditioned on another one. That is the price
of the concentration the boxes show, and it is why ``method="nre"`` at one
round — which stays amortised — is still the right choice when the same network
is meant to serve many observations.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from typing import Any

import astropy.units as u
import numpy as np
import scipy.stats as st

from ampere.core import Dataset, FittingProblem, Model, Parameter, Spectrum
from ampere.inference import SBIEngine

#: The θ the observed spectrum below is generated at. Printed beside every box.
TRUTH = {"model.norm": 3.0, "model.index": -0.8}

SEED = 20260909

#: Where the model is evaluated. Short on purpose: the point of the example is
#: the round loop, not the forward model.
GRID = np.geomspace(1.0, 12.0, 24)

#: The observation's per-sample uncertainty, in Jy.
SIGMA = 0.05


class PowerLaw(Model):
    """``norm * (λ / λ₀) ** index``, with priors whose support is real.

    Both priors are chosen the way ``ampere.inference``'s docstring asks: a
    lognormal on the normalisation, whose support is the positive half-line the
    parameter actually lives on, and a normal on the index, which is genuinely
    unbounded. The unconstrained parameterisation ``SBIEngine`` trains in then
    covers ℝ² exactly, which is what lets a truncation box be a plain
    hyperrectangle.
    """

    def __init__(self, wavelength: np.ndarray) -> None:
        self.register_buffer("wavelength", wavelength, unit=u.um)
        self.register_parameter(Parameter("norm", st.lognorm(0.5, scale=3.0)))
        self.register_parameter(Parameter("index", st.norm(-0.8, 0.5)))

    def evaluate(self, **values: Any) -> Spectrum:
        context = self.context(values)
        wavelength = context["wavelength"]
        return Spectrum(
            wavelength * u.um,
            (context["norm"] * wavelength ** context["index"]) * u.Jy,
        )


def problem(seed: int = SEED) -> FittingProblem:
    """One power law observed once, at :data:`TRUTH` plus Gaussian noise."""
    rng = np.random.default_rng(seed)
    values = TRUTH["model.norm"] * GRID ** TRUTH["model.index"]
    observed = Spectrum(
        GRID * u.um,
        (values + rng.normal(0.0, SIGMA, values.size)) * u.Jy,
        uncertainty=np.full(values.size, SIGMA) * u.Jy,
    )
    return FittingProblem(PowerLaw(GRID), [Dataset(observed)], seed=seed)


def sparkline(values: np.ndarray) -> str:
    """A one-line picture of a marginal, so the shape needs no plotting backend."""
    blocks = " ▁▂▃▄▅▆▇█"
    finite = np.asarray(values, dtype=float)
    span = float(finite.max() - finite.min())
    if span <= 0.0:
        return blocks[0] * finite.size
    scaled = (finite - finite.min()) / span
    return "".join(blocks[round(value * (len(blocks) - 1))] for value in scaled)


def interval(grid: np.ndarray, log_density: np.ndarray, mass: float = 0.9) -> tuple[float, float]:
    """The central *mass* interval of a marginal given on a grid."""
    density = np.exp(np.asarray(log_density, dtype=float))
    weight = np.cumsum(density)
    weight /= weight[-1]
    low = float(np.interp(0.5 * (1.0 - mass), weight, grid))
    high = float(np.interp(1.0 - 0.5 * (1.0 - mass), weight, grid))
    return low, high


def report(run: Any, engine: SBIEngine, elapsed: float) -> None:
    """Everything the run says about its own truncation, printed."""
    attrs = run.attrs
    print(
        f"{attrs['ampere_sbi_method']} on {attrs['ampere_backend']}: "
        f"{attrs['ampere_sbi_rounds']} round(s) x {attrs['ampere_sbi_budget']} simulations, "
        f"{attrs['ampere_sbi_marginal_estimators']} marginal estimator(s), "
        f"{elapsed:.1f}s"
    )
    print(
        f"  amortised: {'yes' if attrs['ampere_sbi_amortised'] else 'no'} "
        f"-- a truncated estimator is chosen at this observation and is "
        f"{'' if attrs['ampere_sbi_amortised'] else 'not amortised across others'}"
    )

    labels = list(engine.problem.free_labels())
    truth = engine.problem.unconstrain(np.array([TRUTH[name] for name in labels]))
    print(f"  truncation box per round (constrained coordinates; truth {TRUTH}):")
    for record in json.loads(attrs["ampere_sbi_truncation"]):
        edges = ", ".join(
            f"{name} [{low:.3f}, {high:.3f}]"
            for name, low, high in zip(
                labels, record["lower_constrained"], record["upper_constrained"], strict=True
            )
        )
        holds = np.all(
            (truth >= np.asarray(record["lower"])) & (truth <= np.asarray(record["upper"]))
        )
        print(
            f"    round {record['round']}: {edges}  "
            f"(log-volume {record['log_volume']:+.3f}, "
            f"prior mass accepted {record['proposal_acceptance']:.1%}, "
            f"contains the truth: {'yes' if holds else 'NO'})"
        )

    group = run["marginals"].dataset
    print("  marginals (estimated 1-D posteriors over the final box):")
    for position, name in enumerate(labels):
        grid = np.asarray(group["grid_constrained"])[position]
        density = np.asarray(group["log_density"])[position]
        low, high = interval(grid, density)
        print(f"    {name:<12} 90% [{low:8.3f}, {high:8.3f}]  {sparkline(density)}")

    posterior = run["posterior"].dataset
    print(f"  posterior (joint estimator, sampled by {attrs['ampere_sbi_truncation_sampler']}):")
    for name in labels:
        drawn = np.asarray(posterior[name]).reshape(-1)
        offset = (drawn.mean() - TRUTH[name]) / drawn.std() if drawn.std() > 0 else math.nan
        print(
            f"    {name:<12} {drawn.mean():8.3f} +/- {drawn.std():.3f}  "
            f"(truth {TRUTH[name]:.3f}, {offset:+.2f} sigma)"
        )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rounds", type=int, default=3, help="rounds of simulate-and-train")
    parser.add_argument("--budget", type=int, default=500, help="simulations per round")
    parser.add_argument("--draws", type=int, default=100, help="posterior draws to store")
    parser.add_argument("--epochs", type=int, default=80, help="training epoch cap")
    parser.add_argument(
        "--marginals", type=int, default=1, choices=(1, 2), help="1-D only, or 1-D and 2-D"
    )
    parser.add_argument(
        "--truncation-epsilon",
        type=float,
        default=None,
        help="threshold on each 1-D marginal, as a fraction of its own maximum",
    )
    parser.add_argument(
        "--sample-with",
        default="rejection",
        choices=("rejection", "mcmc"),
        help="how the final posterior draws; mcmc is the narrow-box fallback",
    )
    args = parser.parse_args(sys.argv[1:] if argv is None else argv)

    settings: dict[str, Any] = dict(
        method="tmnre",
        rounds=args.rounds,
        budget=args.budget,
        marginals=args.marginals,
        sample_with=args.sample_with,
    )
    if args.truncation_epsilon is not None:
        settings["truncation_epsilon"] = args.truncation_epsilon

    engine = SBIEngine(problem(), **settings)
    options: dict[str, Any] = {}
    if args.sample_with == "mcmc":
        options = {"num_chains": 4, "warmup_steps": 20, "thin": 1}
    started = time.perf_counter()
    run = engine.run(
        draws=args.draws,
        training={"max_num_epochs": args.epochs},
        posterior_options=options,
    )
    report(run, engine, time.perf_counter() - started)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
