#!/usr/bin/env python
"""Caching a trained SBI posterior: train once, reuse on every identical rerun.

W3.5's :class:`~ampere.results.ArtefactStore` / :func:`~ampere.results.artefact_key`
seam, demonstrated end to end against a real, small, seeded ``sbi`` 0.27 NPE
posterior. ``SBIEngine`` does not take a ``cache=`` argument yet: W3.3 (the
coordinate-value-mask encoding) had not merged when this item landed, and the
two touch ``ampere/inference/_sbi.py`` at the same time
(``WORK_ITEMS.md``'s dispatch order: "W3.3 first and W3.5 appends"). So this
script demonstrates the *equivalent* — the exact seam ``SBIEngine(cache=store)``
will call, wired here around a small standalone training routine that plays
the role ``SBIEngine.run`` otherwise would.

Run it::

    python examples/sbi/cached_fit.py                    # first run: trains
    python examples/sbi/cached_fit.py                    # second run: a hit
    python examples/sbi/cached_fit.py --budget 400        # a miss; names why

Requires the ``sbi`` extra (``pixi run -e sbi python examples/sbi/cached_fit.py``).

What it demonstrates
---------------------
**A second, identical run trains nothing.** The cache lives at ``--cache-dir``
(default: a directory under the platform cache directory, *not* inside the
package or the repository — ``AGENTS.md`` ground rule 7: no binary artefacts
in git). Delete it, or pass ``--budget`` a different number, to force a miss.

**A miss says which ingredient changed.** :meth:`~ampere.results.ArtefactStore.diff`
compares the requested key against the store's last write and reports the
ingredient(s) that moved — printed here rather than left for a caller to
notice on their own.

**The stored posterior is not the training routine's return value handed back
by reference.** The second run's posterior comes back through :mod:`pickle`,
and its ``sample()`` is exercised to prove the round trip actually happened,
not merely trusted.

**The cache never looks inside what it stores.** :func:`train` returns a real
``sbi`` ``NeuralPosterior``, but :class:`~ampere.results.ArtefactStore` treats
it as opaque bytes — the unconstrained-prior bridge and the coordinate-value-
mask summary layout stay ``ampere/inference/_sbi.py``'s business, not this
module's or the cache's.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import Any

import astropy.units as u
import numpy as np
import scipy.stats as st

from ampere.core import Dataset, FittingProblem, Model, Parameter, Spectrum
from ampere.results import ArtefactStore, artefact_key

GRID = np.array([1.0, 2.0, 3.0, 4.0, 5.0])

#: Outside the repository, as every cache must be (``AGENTS.md`` ground rule 7).
DEFAULT_CACHE_DIR = Path.home() / ".cache" / "ampere-examples" / "sbi-cached-fit"


class Line(Model):
    """A one-parameter linear model — the same toy ``_sbi.py``'s own docstring uses."""

    def __init__(self, wavelength: np.ndarray) -> None:
        self.register_buffer("wavelength", wavelength, unit=u.um)
        self.register_parameter(Parameter("slope", st.norm(2.0, 1.0)))

    def evaluate(self, **values: Any) -> Spectrum:
        ctx = self.context(values)
        return Spectrum(ctx["wavelength"] * u.um, ctx["slope"] * ctx["wavelength"] * u.Jy)


def build_problem(*, seed: int = 20260909) -> FittingProblem:
    """One dataset, one free parameter — small enough to train in seconds."""
    observed = Spectrum(
        GRID * u.um, (2.0 * GRID) * u.Jy, uncertainty=np.full(GRID.size, 0.3) * u.Jy
    )
    return FittingProblem(Line(GRID), [Dataset(observed)], seed=seed)


def train(problem: FittingProblem, *, budget: int, epochs: int) -> Any:
    """Simulate *budget* draws from the prior and train a tiny sbi NPE posterior.

    Stands in for the training half of ``SBIEngine.run``: the unconstrained
    prior bridge and the ``"flat"`` summary layout are ``_sbi.py``'s to own
    until W3.3 lands, so this uses ``sbi``'s own prior directly and a bare
    observed-values vector as the summary — exactly the kind of detail
    :class:`~ampere.results.ArtefactStore` is designed not to need to know:
    it caches whatever this function returns, unopened.
    """
    import torch
    from sbi.inference import NPE
    from sbi.utils import BoxUniform

    rng = problem.rng("cached_fit.train")
    torch.manual_seed(int(rng.integers(0, 2**31 - 1)))

    label = next(iter(problem.datasets))
    observed_values = np.array(problem.datasets[label].observed.values, dtype=float, copy=True)

    prior = BoxUniform(low=torch.tensor([-1.0]), high=torch.tensor([5.0]))
    theta = prior.sample((budget,))
    rows = [
        np.asarray(Line(GRID).evaluate(slope=float(value)).values, dtype=float)
        for value in theta.reshape(-1)
    ]
    x = torch.as_tensor(np.stack(rows), dtype=torch.float32)

    trainer = NPE(prior=prior)
    trainer.append_simulations(theta, x)
    estimator = trainer.train(max_num_epochs=epochs, show_train_summary=False)
    posterior = trainer.build_posterior(estimator)
    posterior.set_default_x(torch.as_tensor(observed_values, dtype=torch.float32))
    return posterior


def fit_cached(
    *,
    budget: int = 200,
    epochs: int = 20,
    cache_dir: Path = DEFAULT_CACHE_DIR,
) -> tuple[Any, bool, dict[str, tuple[Any, Any]]]:
    """Train, or load, a posterior for :func:`build_problem` — the W3.5 seam.

    Returns the posterior, whether it was a cache hit, and
    :meth:`~ampere.results.ArtefactStore.diff` against the store's last write
    (empty on a hit, or when nothing was ever cached here before).
    """
    problem = build_problem()
    store = ArtefactStore(cache_dir)
    key = artefact_key(
        problem,
        layout="flat",
        method="npe",
        architecture="maf",
        budget=budget,
        rounds=1,
    )
    diff = store.diff(key)
    posterior, hit = store.train_or_load(key, lambda: train(problem, budget=budget, epochs=epochs))
    return posterior, hit, diff


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--budget", type=int, default=200, help="simulations to train on")
    parser.add_argument("--epochs", type=int, default=20, help="training epoch cap")
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR)
    args = parser.parse_args(sys.argv[1:] if argv is None else argv)

    started = time.perf_counter()
    posterior, hit, diff = fit_cached(
        budget=args.budget, epochs=args.epochs, cache_dir=args.cache_dir
    )
    elapsed = time.perf_counter() - started

    if hit:
        print(f"cache hit: trained nothing ({elapsed:.2f}s to load and verify)")
    else:
        print(f"cache miss: trained a new posterior ({elapsed:.2f}s)")
        if diff:
            named = ", ".join(
                f"{name} ({old!r} -> {new!r})" for name, (old, new) in sorted(diff.items())
            )
            print(f"  ingredient(s) that changed since the last write here: {named}")
        else:
            print("  nothing was cached here before")
    print(f"  cache directory: {args.cache_dir}")

    draws = posterior.sample((5,), show_progress_bars=False)
    print(f"  five posterior draws of slope: {[round(v, 3) for v in draws.reshape(-1).tolist()]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
