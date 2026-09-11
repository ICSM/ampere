#!/usr/bin/env python
"""The bare case: an in-process NPE fit of a native problem, nothing else going on.

``external_simulator.py``/``fit_external_simulator.py`` show the case
:class:`~ampere.inference.SBIEngine` exists for — a compiled simulator with no
likelihood to write down, fitted through a pool of subprocesses.
``cached_fit.py`` shows the training-cache seam. Both are true stories, and
both obscure the smallest one: an ordinary ampere
:class:`~ampere.core.FittingProblem`, composed from one of the *native*
backends, fitted with :class:`~ampere.inference.SBIEngine` directly — no
subprocess, no ``cache=``, no truncation, one call.

What makes the problem "native" is nothing this script does. It is that
:func:`build_problem` builds its model from :mod:`ampere.backends.torch`
pieces, so ``problem.backend`` is ``"torch"`` and ``problem.batchable`` is
``True``. :meth:`~ampere.core.dataset.FittingProblem.simulate_many` — which
:meth:`~ampere.inference.SBIEngine.run` calls, on its own training budget —
then takes the native batched forward path (W3.1 slice 2) **by default**:
``native=None``, "use it if it works". The whole budget goes through one
``vmap``-ped forward call per chunk instead of a Python loop, and
:class:`~ampere.inference.SBIEngine` never asks for that by name — it calls
``problem.simulate_many(...)`` exactly as it would of a reference-backend
problem or a subprocess-wrapped one. Swap the ``ampere.backends.torch``
import in :func:`build_problem` for ``ampere.backends.jax`` and nothing else
in this script changes.

Run it::

    python examples/sbi/npe_native.py                 # budget 2000, 500 draws
    python examples/sbi/npe_native.py --budget 5000

Requires the ``sbi`` extra, which brings in torch too
(``pixi run -e sbi python examples/sbi/npe_native.py``).
"""

from __future__ import annotations

import argparse
import sys
import time
from typing import Any

import astropy.units as u
import numpy as np
import scipy.stats as st

from ampere.core import Dataset, FittingProblem, GaussianFamily, Likelihood, Spectrum
from ampere.inference import SBIEngine

GRID = np.geomspace(1.0, 10.0, 12)
TRUTH = {"model.norm": 2.0, "model.index": -1.0}
UNCERTAINTY = 0.05


def build_problem(*, seed: int = 20260911) -> FittingProblem:
    """One :class:`~ampere.backends.torch.PowerLaw`, fitted to synthetic data.

    Small on purpose: this script is about the seam between ``SBIEngine`` and
    a native realisation, not about the model. ``norm`` and ``index`` keep
    their reference-backend priors — a native problem is fitted exactly like
    any other, priors included. The likelihood's noise model is torch's own
    :class:`~ampere.backends.torch.IndependentNoise`, not
    :class:`~ampere.core.Dataset`'s reference-backend default: every piece of
    a native problem must declare the same backend (``DatasetError`` refuses
    a mixed one by name), and the family (:class:`~ampere.core.GaussianFamily`)
    is the one backend-neutral part of a likelihood.
    """
    from ampere.backends.torch import IndependentNoise, PowerLaw

    rng = np.random.default_rng(seed)
    truth = TRUTH["model.norm"] * GRID**TRUTH["model.index"]
    observed = Spectrum(
        GRID * u.um,
        (truth + rng.normal(0.0, UNCERTAINTY, GRID.size)) * u.Jy,
        uncertainty=np.full(GRID.size, UNCERTAINTY) * u.Jy,
    )
    model = PowerLaw(GRID, norm=st.lognorm(0.3, scale=2.0), index=st.norm(-1.0, 0.3))
    likelihood = Likelihood(GaussianFamily(), IndependentNoise())
    problem = FittingProblem(model, [Dataset(observed, likelihood=likelihood)], seed=seed)
    assert problem.backend == "torch" and problem.batchable  # the point of the example
    return problem


def fit(problem: FittingProblem, *, budget: int, draws: int) -> Any:
    """One line: no ``executor=``, no ``cache=``, ``method="npe"``'s one round."""
    return SBIEngine(problem, budget=int(budget)).run(draws=int(draws))


def report(run: Any) -> str:
    attrs = run.attrs
    posterior = run["posterior"].dataset
    lines = [
        f"{attrs['ampere_sbi_method']} on {attrs['ampere_backend']} (native, batchable): "
        f"{attrs['ampere_sbi_simulations']} simulation(s), "
        f"{attrs['ampere_sbi_usable_simulations']} usable "
        f"({attrs['ampere_sbi_failures']} failed)",
        "  posterior (truth in brackets):",
    ]
    for name in sorted(posterior.data_vars):
        values = np.asarray(posterior[name], dtype=float).ravel()
        lower, upper = np.percentile(values, [16.0, 84.0])
        truth = TRUTH.get(name)
        bracket = f"  (truth {truth:+.3f})" if truth is not None else ""
        lines.append(
            f"    {name:14s} {values.mean():+.4f} +- {values.std():.4f}   "
            f"[{lower:+.4f}, {upper:+.4f}]{bracket}"
        )
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--budget", type=int, default=2000, help="simulations to train on")
    parser.add_argument("--draws", type=int, default=500, help="posterior draws to return")
    args = parser.parse_args(sys.argv[1:] if argv is None else argv)

    problem = build_problem()
    started = time.perf_counter()
    run = fit(problem, budget=args.budget, draws=args.draws)
    elapsed = time.perf_counter() - started

    print(report(run))
    print(f"  {elapsed:.1f} s wall clock")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
