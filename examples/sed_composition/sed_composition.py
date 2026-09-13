"""One :class:`~ampere.backends.reference.ModifiedBlackBody`, fitted to a spectrum through
:class:`~ampere.core.Instrument`'s LSF-then-resample chain and to photometry
through synthetic photometry, both bound to the same model channel — memo
§7.1, made runnable (W4.11). See :doc:`the tutorial page </sed_composition>`
for the walk-through; this module is the code it walks through.

The physics
-----------
:class:`~ampere.backends.reference.ModifiedBlackBody` (or its ``torch``/``jax``
twin) is a 180 K, optically thin dust greybody, emissivity index 1.6 — cool
enough that it has essentially nothing to say in the near infrared and most of
its flux beyond 20 micron, which is realistic for the kind of source this
instrument combination targets (a debris disk or an embedded protostar, not a
photosphere) and is why :data:`FILTERS` is mid/far-infrared only.
``scale`` folds in the source's solid angle, so it is a
dimensionless number of order :math:`10^{-15}` — not a flux — chosen so that
the flux the two instruments actually *see* comes out at a few Jy, a sensible
size for a catalogued mid/far-infrared source (see :data:`generators.TRUTH`).

The two instruments
--------------------
**irs** — a slit spectrograph: :class:`~ampere.backends.reference.LSFConvolution` at
constant resolving power 100, :class:`~ampere.backends.reference.Resample` onto a 120-point
observed grid from 5 to 35 micron, and a :class:`~ampere.backends.reference.CalibrationScale`
nuisance parameter for the instrument's own flux calibration (deliberately
mis-set 2 % high in the synthetic data — :data:`generators.CALIBRATION_TRUTH`).

**catalogue** — :meth:`~ampere.backends.reference.SyntheticPhotometry.from_library`, five
bundled mid/far-infrared filters (WISE, IRAS, Spitzer MIPS, Herschel PACS —
see :data:`FILTERS`) tabulated on a 500-point grid.

Both bind the model's one channel, ``"sed"``, so — per
``transformations.md`` §15 Q4 — they need distinct labels; this module gives
them ``"irs"`` and ``"catalogue"`` explicitly, and :doc:`the tutorial page
</sed_composition>` shows what happens (and why) if you leave that out.

Backends
--------
:func:`build_problem` takes a backend name and builds every piece from that
backend's own namespace, noise model included (the composition's one-backend
rule; see ``docs/source/overview.rst``). On ``"reference"`` the model has no
gradient, so :func:`fit` samples with :class:`~ampere.inference.EmceeEngine`;
on ``"torch"`` or ``"jax"`` it samples with
:class:`~ampere.inference.NUTSEngine`. Nothing about the composition itself
changes — same model, same two instruments, same two datasets — only which
module supplied the pieces and which engine can consume the result.

Run it::

    python -m examples.sed_composition                      # reference, emcee, ~1 minute
    python -m examples.sed_composition --backend torch       # NUTS
    python -m examples.sed_composition --walkers 32 --steps 200

:mod:`tests.examples.test_sed_composition` is this module's own coverage: a
fast, always-on suite at a tiny budget. The full-budget 95 % recovery this
item is accepted on is not part of that suite — see its docstring for why —
and is instead verified by running this module directly; the branch report
for W4.11 has the numbers from doing exactly that.
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
from ampere.inference import EmceeEngine, NUTSEngine

from . import generators

__all__ = [
    "BACKENDS",
    "DEFAULT_BURN_IN",
    "DEFAULT_CHAINS",
    "DEFAULT_DRAWS",
    "DEFAULT_STEPS",
    "DEFAULT_WALKERS",
    "DEFAULT_WARMUP",
    "FILTERS",
    "GRID",
    "OBSERVED_WAVELENGTH",
    "PHOTOMETRY_TABULATION",
    "QUALIFIED_TRUTH",
    "RESOLVING_POWER",
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

#: The model's own grid — a fallback only: negotiation replaces it the moment
#: an instrument is bound (see :mod:`.generators`).
GRID = np.geomspace(1.0, 200.0, 2000)

#: The spectrograph's own observed grid, micron. A ``Resample`` step takes
#: this from the observed container it targets, never recomputes it
#: (``transformations.md`` §10) — this constant *is* that container's grid.
OBSERVED_WAVELENGTH = np.linspace(5.0, 35.0, 120)

#: The photometry step's response tabulation, micron — becomes its exact
#: ``points=`` requirement.
PHOTOMETRY_TABULATION = np.geomspace(1.0, 100.0, 500)

RESOLVING_POWER = 100.0

#: Chosen so every filter actually sees the 180 K source (a few tenths to a
#: few Jy — see the tutorial page): 2MASS/near-IR bands were tried first and
#: dropped, because a 180 K greybody puts essentially nothing there
#: (down to 1e-19 Jy) and NUTS's gradients through such an extreme dynamic
#: range made warmup pathologically slow on the torch/jax variant, for no
#: statistical benefit (a point with no signal constrains nothing).
FILTERS = (
    "WISE_RSR_W3",
    "WISE_RSR_W4",
    "IRAS_60",
    "SPITZER_MIPS_70",
    "HERSCHEL_PACS_100",
)

#: The truth, qualified by the names a built :class:`~ampere.core.dataset.FittingProblem`
#: actually samples — ``problem.parameters.free_names`` once the two datasets
#: and the model are merged.
QUALIFIED_TRUTH: dict[str, float] = {
    "model.temperature": generators.TRUTH["temperature"],
    "model.beta": generators.TRUTH["beta"],
    "model.scale": generators.TRUTH["scale"],
    "irs.instrument.calibration_scale.scale": generators.CALIBRATION_TRUTH,
}

# Budgets. Reference/emcee: measured in the high teens of milliseconds per
# evaluation on this composition after W4.0's LSF operator cache (twenty
# calls to ``problem.log_prob`` at the truth; see the tutorial page for the
# number measured there, and __main__.py for why single-threaded BLAS
# matters here), so 16 walkers x 450 steps (7 200 evaluations) is
# comfortably under two minutes and the 270 kept steps recover all four
# parameters inside their central 95 %. torch/jax NUTS needs far fewer
# gradient-based draws for the same coverage.
DEFAULT_WALKERS = 16
DEFAULT_STEPS = 450
DEFAULT_BURN_IN = 180
DEFAULT_DRAWS = 500
DEFAULT_WARMUP = 500
DEFAULT_CHAINS = 2


def backend_module(backend: str) -> Any:
    """The namespace *backend*'s model and instrument pieces live in.

    ``"jax"`` is configured for float64 here, before anything of its is
    touched, matching ``docs/source/overview.rst``'s "before any jax work"
    rule.
    """
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

    ``ampere.core.IndependentNoise`` is ``IndependentNoise``'s actual home and
    already declares ``BACKEND = "reference"`` (it is "the reference path's
    honest answer", per ``ampere/core/likelihood.py``'s ``NoiseModel``
    docstring) — :mod:`ampere.backends.reference` does not re-export it, so it
    is imported from ``ampere.core`` directly. Composing it into a torch or
    jax problem instead would be exactly the "numpy island" the one-backend
    rule refuses by name (``docs/source/overview.rst``): each native backend
    supplies its own.
    """
    if backend == "reference":
        import ampere.core as module

        return module
    return backend_module(backend)


def build_model(backend: str) -> Any:
    """The one :class:`~ampere.backends.reference.ModifiedBlackBody`, on *backend*."""
    module = backend_module(backend)
    return module.ModifiedBlackBody(
        GRID,
        temperature=st.uniform(50.0, 400.0),
        beta=st.uniform(0.5, 2.5),
        scale=st.loguniform(1e-16, 1e-14),
        channels="sed",
    )


def build_instruments(backend: str) -> tuple[Instrument, Instrument]:
    """The spectrograph (``"irs"``) and the catalogue camera (``"catalogue"``).

    Both bind channel ``"sed"``; distinct labels are passed explicitly rather
    than left to default, which is exactly what a second instrument on one
    channel needs — see :doc:`the tutorial page </sed_composition>` for what
    happens without them.
    """
    module = backend_module(backend)
    spectrograph = Instrument(
        [
            module.LSFConvolution(resolving_power=RESOLVING_POWER),
            module.Resample(OBSERVED_WAVELENGTH),
            module.CalibrationScale(st.lognorm(0.05, scale=1.0)),
        ],
        channel="sed",
        label="irs",
    )
    camera = Instrument(
        [module.SyntheticPhotometry.from_library(list(FILTERS), PHOTOMETRY_TABULATION)],
        channel="sed",
        label="catalogue",
    )
    return spectrograph, camera


def build_problem(backend: str = "reference", *, seed: int = generators.SEED) -> FittingProblem:
    """The composed problem: one model, two datasets, distinct labels.

    Builds the model and the two instruments from *backend*'s own pieces,
    generates the synthetic data through :func:`generators.synthetic_data`
    (the negotiate-then-compile_for dance — see that module), and binds
    everything into one :class:`~ampere.core.dataset.FittingProblem`.
    """
    model = build_model(backend)
    spectrograph, camera = build_instruments(backend)
    observed_spectrum, observed_photometry = generators.synthetic_data(
        model, spectrograph, camera, seed=seed
    )
    noise = noise_module(backend)

    def gaussian() -> Likelihood:
        return Likelihood(GaussianFamily(), noise.IndependentNoise())

    datasets = DatasetCollection(
        {
            "irs": Dataset(observed_spectrum, spectrograph, likelihood=gaussian()),
            "catalogue": Dataset(observed_photometry, camera, likelihood=gaussian()),
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
    """Sample *problem* with the engine its backend can feed.

    Dispatches on *backend* rather than discovering it from the problem, so a
    caller who asked for ``"torch"`` and got a reference-backend problem by
    accident is told so by :class:`~ampere.inference.NUTSEngine` itself
    (``differentiable=False``) rather than silently sampling the wrong thing.
    """
    if backend == "reference":
        engine = EmceeEngine(problem, walkers=DEFAULT_WALKERS if walkers is None else walkers)
        return engine.run(
            DEFAULT_STEPS if steps is None else steps,
            burn_in=DEFAULT_BURN_IN if burn_in is None else burn_in,
            progress=progress,
        )
    engine = NUTSEngine(problem)
    return engine.run(
        DEFAULT_DRAWS if draws is None else draws,
        warmup=DEFAULT_WARMUP if warmup is None else warmup,
        chains=DEFAULT_CHAINS if chains is None else chains,
        progress=progress,
    )


def recovers_truth(run: Any, *, level: float = 0.95) -> dict[str, bool]:
    """Whether each parameter's central *level* interval contains its truth.

    The acceptance statistic W4.11 is measured against: not the point
    estimate, the interval.
    """
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
            f"    {name:45s} {values.mean():+.6g} +- {values.std():.3g}   "
            f"95%[{lower:+.6g}, {upper:+.6g}]{bracket}{flag}"
        )
    return "\n".join(lines)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend", default="reference", choices=list(BACKENDS))
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

    problem = build_problem(args.backend, seed=args.seed)
    req = problem.requirements["model"]["sed"]
    print(f"negotiated channels: {list(problem.requirements['model'])}")
    print(f"sources asking of channel 'sed': {req.sources}")
    for axis, axis_req in req.axes.items():
        print(f"  axis {axis}: {axis_req}")
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
