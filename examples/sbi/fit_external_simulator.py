#!/usr/bin/env python
"""Fitting a compiled simulator with SBI: the case the whole layer exists for.

``external_simulator.py`` stops at the ``(θ, x)`` pairs, because that is where
W3.1 slice 1 stops. This script picks them up: it composes the *same*
black-box model — a subprocess call to a program that has never heard of
ampere — and hands the problem to :class:`~ampere.inference.SBIEngine`, which
trains a neural posterior estimator on a simulation budget and evaluates it at
the observed data.

Nothing about the fit knows that the forward model is a subprocess. That is
the claim ``DEVELOPMENT_PLAN.md`` §4.5 makes about the engine-facing surface,
and this file is what it looks like when the claim is true: the same
``SBIEngine(problem, ...)`` that fits a hand-written numpy model fits this,
and ``executor=`` is the only line that mentions parallelism.

Run it::

    python examples/sbi/fit_external_simulator.py                 # 400 draws, 4 workers
    python examples/sbi/fit_external_simulator.py 2000 8          # 2000 draws, 8 workers
    python examples/sbi/fit_external_simulator.py 400 4 --training-set pairs.nc

Requires the ``sbi`` extra (``pixi run -e sbi python ...``); the refusal names
it if it is missing.

What it demonstrates that ``external_simulator.py`` does not
------------------------------------------------------------
**A budget spent through a pool is a fit.** ``executor=ProcessExecutor(n)``
and ``chunk_size=`` go straight through to
:meth:`~ampere.core.dataset.FittingProblem.simulate_many`; the engine never
sees a worker. Failed draws — this simulator crashes for ``norm`` above 9.5 —
are dropped from what the network trains on and counted in the run's attrs,
which is ``inference.md`` §13's reject-and-record arriving where a user can
read it.

**The posterior is scored honestly.** Every stored draw is evaluated on the
numpy contract path, so the run's ``lp``, ``log_prior`` and ``log_likelihood``
are the true ones and not the network's opinion of them. The network's own
log-density is stored beside them as ``ampere_sbi_log_prob``. For this toy the
two can be compared directly, because the toy's likelihood *is* writable — a
Gaussian on a power law — which is exactly why it makes a good demonstration
and a bad advertisement: for the real codes this exists for, the second column
is all there would be.

**The pairs can be kept.** ``--training-set`` writes the whole budget to a
netCDF training set (``results.md`` §11), failures included, a chunk at a
time, so a budget larger than memory reaches the file without being held.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

from ampere.core import ProcessExecutor
from ampere.inference import SBIEngine

# This directory goes on sys.path and stays there. A worker process has to be
# able to *re-import* `external_simulator` to unpickle the model it is handed,
# and a module loaded from a file path under a synthetic name cannot be
# re-imported anywhere -- the same reasoning
# tests/examples/test_external_simulator.py records.
_HERE = str(Path(__file__).resolve().parent)
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from external_simulator import build_problem  # noqa: E402


def fit(
    budget: int = 400,
    workers: int = 4,
    *,
    draws: int = 1000,
    method: str = "npe",
    chunk_size: int | None = None,
    training_set: str | None = None,
) -> Any:
    """Simulate *budget* draws through *workers* subprocesses, then fit.

    ``chunk_size`` bounds how many simulations exist at once; with a
    ``training_set`` it is also the size of one write, and few large chunks
    beat many small ones (``results.md`` limitation 13.9's append is
    ``O(existing + new)``).
    """
    problem = build_problem()
    engine = SBIEngine(
        problem,
        method=method,
        budget=int(budget),
        executor=ProcessExecutor(int(workers)),
        chunk_size=chunk_size,
        training_set=training_set,
    )
    return engine.run(draws=int(draws))


def report(run: Any) -> str:
    """What a driver should print after an SBI fit: the posterior, and its provenance."""
    attrs = run.attrs
    posterior = run["posterior"].dataset
    lines = [
        f"{attrs['ampere_sbi_method']} on {attrs['ampere_backend']}: "
        f"{attrs['ampere_sbi_simulations']} simulation(s), "
        f"{attrs['ampere_sbi_usable_simulations']} usable "
        f"({attrs['ampere_sbi_failures']} failed)",
        f"  network: {attrs['ampere_sbi_density_estimator']}, "
        f"embedding {attrs['ampere_sbi_embedding']}, "
        f"{attrs['ampere_sbi_epochs_trained']} epoch(s), "
        f"summary layout {attrs['ampere_sbi_summary_layout']} "
        f"({attrs['ampere_sbi_summary_features']} feature(s))",
        f"  prior in {attrs['ampere_sbi_parameterisation']} coordinates; "
        f"context {attrs['ampere_sbi_context']}",
        "  posterior:",
    ]
    for name in sorted(posterior.data_vars):
        values = np.asarray(posterior[name], dtype=float).ravel()
        lower, upper = np.percentile(values, [16.0, 84.0])
        lines.append(
            f"    {name:14s} {values.mean():+.4f} +- {values.std():.4f}   "
            f"[{lower:+.4f}, {upper:+.4f}]"
        )
    stats = run["sample_stats"].dataset
    true_lp = np.asarray(stats["lp"], dtype=float).ravel()
    estimator_lp = np.asarray(stats["ampere_sbi_log_prob"], dtype=float).ravel()
    lines.append(
        f"  true log p(theta, x) over the draws: {true_lp.mean():+.3f} +- {true_lp.std():.3f}"
    )
    lines.append(
        f"  estimator log-density ({attrs['ampere_sbi_log_prob_kind']}): "
        f"{estimator_lp.mean():+.3f} +- {estimator_lp.std():.3f}"
    )
    if "ampere_sbi_training_set" in attrs:
        lines.append(f"  pairs written to {attrs['ampere_sbi_training_set']}")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    arguments = list(sys.argv[1:] if argv is None else argv)
    training_set = None
    if "--training-set" in arguments:
        index = arguments.index("--training-set")
        training_set = arguments[index + 1]
        del arguments[index : index + 2]
    positional = [item for item in arguments if not item.startswith("-")]
    budget = int(positional[0]) if positional else 400
    workers = int(positional[1]) if len(positional) > 1 else 4

    started = time.perf_counter()
    run = fit(
        budget,
        workers,
        chunk_size=max(1, budget // 2),
        training_set=training_set,
    )
    print(report(run))
    print(f"  {time.perf_counter() - started:.1f} s wall clock")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
