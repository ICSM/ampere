"""One :class:`~ampere.backends.reference.ReflexOrbit`, two ``TimeSeries`` datasets:
right ascension and declination, sharing one orbital-parameter evaluation —
the second modality W4.4's template was written to prove (W4.9). See
:doc:`the tutorial page </astrometry>` for the walk-through; this module is
the code it walks through.

The physics
-----------
A source drifting under proper motion and wobbling under a companion's reflex
motion (a circular, node-aligned orbit projected onto the sky — see
:mod:`ampere.backends.reference.astrometry`'s module docstring for why this
is the design sketch's own choice rather than an omission) emits two
channels, ``"ra"`` and ``"dec"``, from one model evaluation: ``period`` and
``phase`` are shared at the language level, needing no ``Tie``.

The two instruments
--------------------
**astrom_ra**, **astrom_dec** — each one ``EpochSample`` step, pinning the
model onto the observed epochs; no parameters, no arithmetic. Both bind
different channels (``"ra"``, ``"dec"``), so — unlike
:mod:`examples.sed_composition` and :mod:`examples.interferometry`, whose two
instruments share one channel — there is no label collision to walk through
here; distinct labels are given anyway, for symmetry with the other pages.

The likelihood
--------------
Two arms, chosen with ``--gp``: independent Gaussian noise (the rigid
comparison), or :class:`~ampere.core.GaussianProcessNoise` with
:class:`~ampere.core.Matern32` on the ``QuasisepGP`` solver — the O(N) path a
:class:`~ampere.core.TimeSeries`'s one ordered ``time`` axis is exactly the
right shape for. Independent noise is the default because the recovery this
item is accepted on is about the orbital parameters, not about the flexible
likelihood's own calibration claim (that is M2's question, asked of a
misspecified sky model, not of this modality's simple truth).

Backends
--------
:func:`build_problem` takes a backend name and builds every piece from that
backend's own namespace (the composition's one-backend rule). On
``"reference"`` the model has no gradient, so :func:`fit` samples with
:class:`~ampere.inference.EmceeEngine`; on ``"torch"`` or ``"jax"`` it samples
with :class:`~ampere.inference.NUTSEngine`.

Run it::

    python -m examples.astrometry                      # reference, emcee
    python -m examples.astrometry --backend torch       # NUTS
    python -m examples.astrometry --gp                  # the flexible likelihood

:mod:`tests.examples.test_astrometry_example` is this module's own coverage:
a fast, always-on suite at a tiny budget. The full-budget recovery this item
is accepted on is verified by running this module directly and by
``tests/astrometry``'s own pinned rows at a per-PR budget larger than the
smoke test's but far short of the full one.
"""

from __future__ import annotations

import argparse
import sys
import time
from typing import Any

import numpy as np
import scipy.stats as st

from ampere.core import (
    Dataset,
    DatasetCollection,
    FittingProblem,
    GaussianFamily,
    Instrument,
    Likelihood,
)

from . import generators

__all__ = [
    "BACKENDS",
    "DEFAULT_BURN_IN",
    "DEFAULT_CHAINS",
    "DEFAULT_DRAWS",
    "DEFAULT_STEPS",
    "DEFAULT_WALKERS",
    "DEFAULT_WARMUP",
    "QUALIFIED_TRUTH",
    "backend_module",
    "build_instruments",
    "build_model",
    "build_problem",
    "fit",
    "main",
    "noise_module",
    "recovers_truth",
    "report",
]

BACKENDS = ("reference", "torch", "jax")

#: The truth, qualified by the names a built ``FittingProblem`` actually
#: samples.
QUALIFIED_TRUTH: dict[str, float] = {
    "model.pmra": generators.TRUTH["pmra"],
    "model.pmdec": generators.TRUTH["pmdec"],
    "model.period": generators.TRUTH["period"],
    "model.phase": generators.TRUTH["phase"],
    "model.amp_ra": generators.TRUTH["amp_ra"],
    "model.amp_dec": generators.TRUTH["amp_dec"],
}

# Budgets. Twelve epochs and six free parameters is a small problem, so a
# reference/emcee fit is fast even at a generous walker/step count; torch/jax
# NUTS needs far fewer gradient-based draws for the same coverage.
DEFAULT_WALKERS = 24
DEFAULT_STEPS = 1500
DEFAULT_BURN_IN = 500
DEFAULT_DRAWS = 500
DEFAULT_WARMUP = 500
DEFAULT_CHAINS = 2

#: The flexible likelihood's own hyperparameters, held fixed here (the point
#: of this study is recovering the orbit, not fitting the noise process).
GP_AMPLITUDE = 0.03
GP_LENGTH_SCALE = 120.0


def backend_module(backend: str) -> Any:
    """The namespace *backend*'s model and instrument pieces live in."""
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
    raise ValueError(f"unknown backend {backend!r}; the three are {BACKENDS}.")


def noise_module(backend: str) -> Any:
    """Where ``IndependentNoise`` comes from for *backend*.

    ``ampere.core.IndependentNoise`` already declares ``BACKEND =
    "reference"``, so on the reference backend it is imported from
    ``ampere.core`` directly rather than through the backend package, exactly
    as :mod:`examples.sed_composition.sed_composition` does.
    """
    if backend == "reference":
        import ampere.core as module

        return module
    return backend_module(backend)


def build_model(backend: str) -> Any:
    """The one :class:`~ampere.backends.reference.ReflexOrbit`, on *backend*."""
    module = backend_module(backend)
    return module.ReflexOrbit(
        generators.EPOCHS,
        pmra=st.norm(0.0, 5.0),
        pmdec=st.norm(0.0, 5.0),
        # A period search over many decades is a genuinely multi-modal
        # problem at this modality's sparse, irregular sampling -- unlike
        # every model the interferometry template fits, a reflex orbit is
        # periodic, and a period prior wide enough to reach past the epochs'
        # own baseline lets a sampler alias onto a spurious cycle indefinitely
        # (see docs/source/astrometry.rst's closing section). A period search
        # informed to within a few tens of days -- the realistic case once a
        # periodogram or a previous epoch has suggested roughly where to look
        # -- is what this study fits, rather than a global period search,
        # which is a different (and harder) problem this item does not claim
        # to solve.
        period=st.norm(400.0, 30.0),
        phase=st.uniform(0.0, 2.0 * np.pi),
        amp_ra=st.uniform(0.0, 2.0),
        amp_dec=st.uniform(0.0, 2.0),
    )


def build_instruments(backend: str) -> tuple[Instrument, Instrument]:
    """The two epoch-sampling instruments, ``"astrom_ra"`` and ``"astrom_dec"``."""
    module = backend_module(backend)
    ra_instrument = Instrument(
        [module.EpochSample(generators.EPOCHS)], channel="ra", label="astrom_ra"
    )
    dec_instrument = Instrument(
        [module.EpochSample(generators.EPOCHS)], channel="dec", label="astrom_dec"
    )
    return ra_instrument, dec_instrument


def build_problem(
    backend: str = "reference", *, gp: bool = False, seed: int = generators.SEED
) -> FittingProblem:
    """The composed problem: one model, two channels, distinct labels.

    ``gp=True`` swaps independent noise for the flexible likelihood
    (``GaussianProcessNoise(Matern32(...), QuasisepGP())``) on both channels.
    """
    model = build_model(backend)
    ra_instrument, dec_instrument = build_instruments(backend)
    observed_ra, observed_dec = generators.synthetic_data(
        model, ra_instrument, dec_instrument, seed=seed
    )
    noise = noise_module(backend)

    def likelihood() -> Likelihood:
        if not gp:
            return Likelihood(GaussianFamily(), noise.IndependentNoise())
        kernel = noise.Matern32(GP_AMPLITUDE, GP_LENGTH_SCALE, axes=("time",))
        return Likelihood(GaussianFamily(), noise.GaussianProcessNoise(kernel, noise.QuasisepGP()))

    datasets = DatasetCollection(
        {
            "ra": Dataset(observed_ra, ra_instrument, likelihood=likelihood(), label="ra"),
            "dec": Dataset(observed_dec, dec_instrument, likelihood=likelihood(), label="dec"),
        }
    )
    return FittingProblem(model, datasets, seed=seed)


def fit(
    problem: FittingProblem,
    *,
    backend: str,
    walkers: int | None = None,
    steps: int | None = None,
    burn_in: int | None = None,
    draws: int | None = None,
    warmup: int | None = None,
    chains: int | None = None,
    progress: bool = False,
) -> Any:
    """Sample *problem* with the engine its backend can feed."""
    if backend == "reference":
        from ampere.inference import EmceeEngine

        engine = EmceeEngine(problem, walkers=DEFAULT_WALKERS if walkers is None else walkers)
        return engine.run(
            DEFAULT_STEPS if steps is None else steps,
            burn_in=DEFAULT_BURN_IN if burn_in is None else burn_in,
            progress=progress,
        )
    from ampere.inference import NUTSEngine

    engine = NUTSEngine(problem)
    return engine.run(
        DEFAULT_DRAWS if draws is None else draws,
        warmup=DEFAULT_WARMUP if warmup is None else warmup,
        chains=DEFAULT_CHAINS if chains is None else chains,
        progress=progress,
    )


def recovers_truth(run: Any, *, level: float = 0.95) -> dict[str, bool]:
    """Whether each parameter's central *level* interval contains its truth."""
    posterior = run["posterior"].dataset
    tail = (1.0 - level) / 2.0 * 100.0
    covered: dict[str, bool] = {}
    for name, truth in QUALIFIED_TRUTH.items():
        draws = np.asarray(posterior[name], dtype=float).ravel()
        lower, upper = np.percentile(draws, [tail, 100.0 - tail])
        covered[name] = bool(lower <= truth <= upper)
    return covered


def report(run: Any) -> str:
    """A human-readable posterior summary, truth in brackets, 95 % coverage flagged."""
    attrs = run.attrs
    posterior = run["posterior"].dataset
    covered = recovers_truth(run)
    lines = [
        (
            f"{attrs['ampere_engine']} on {attrs['ampere_backend']}: "
            f"{posterior.sizes['chain']} chain(s) x {posterior.sizes['draw']} draw(s)"
        ),
        "  posterior (truth in brackets, 95 % coverage flagged):",
    ]
    for name in sorted(posterior.data_vars):
        values = np.asarray(posterior[name], dtype=float).ravel()
        lower, upper = np.percentile(values, [2.5, 97.5])
        truth = QUALIFIED_TRUTH.get(name)
        flag = "" if truth is None else ("  ok" if covered[name] else "  MISS")
        bracket = "" if truth is None else f"  (truth {truth:+.6g})"
        lines.append(
            f"    {name:20s} {values.mean():+.6g} +- {values.std():.3g}   "
            f"95%[{lower:+.6g}, {upper:+.6g}]{bracket}{flag}"
        )
    return "\n".join(lines)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend", default="reference", choices=list(BACKENDS))
    parser.add_argument("--gp", action="store_true", help="the flexible likelihood")
    parser.add_argument("--seed", type=int, default=generators.SEED)
    parser.add_argument("--walkers", type=int, default=None, help="emcee only")
    parser.add_argument("--steps", type=int, default=None, help="emcee only")
    parser.add_argument("--burn-in", type=int, default=None, help="emcee only")
    parser.add_argument("--draws", type=int, default=None, help="NUTS only")
    parser.add_argument("--warmup", type=int, default=None, help="NUTS only")
    parser.add_argument("--chains", type=int, default=None, help="NUTS only")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(sys.argv[1:] if argv is None else argv)

    problem = build_problem(args.backend, gp=args.gp, seed=args.seed)
    print(f"negotiated channels: {list(problem.requirements['model'])}")
    for channel_name in problem.requirements["model"]:
        req = problem.requirements["model"][channel_name]
        print(f"sources asking of channel {channel_name!r}: {req.sources}")
    print(f"free parameters: {problem.parameters.free_names}")

    started = time.perf_counter()
    run = fit(
        problem,
        backend=args.backend,
        walkers=args.walkers,
        steps=args.steps,
        burn_in=args.burn_in,
        draws=args.draws,
        warmup=args.warmup,
        chains=args.chains,
    )
    elapsed = time.perf_counter() - started

    print(report(run))
    print(f"  {elapsed:.1f} s wall clock")
    return 0
