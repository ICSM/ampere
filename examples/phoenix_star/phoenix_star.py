"""The v2 twin of ``examples/examples_paper/phoenixstar.py`` -- W6.13 (4).

The legacy script stays exactly as it is. It fits a PHOENIX stellar model,
emulated by Starfish, to eight-band photometry and a synthetic Gaia-RVS
spectrum with SBI. This package is the same model, data and question
against :mod:`ampere.core`, with one substitution, ruled 2026-09-28
(``docs/design/example_dependencies_memo.md`` section 2): **no Starfish**. The
emulator is ampere's own, a PCA plus one Gaussian process per weight, trained
once by :mod:`.train_emulator` and committed as ``phoenix_emulator.npz``
(:mod:`.emulator`). For the Starfish route, see the legacy script's lines
84-110. Starfish on PyPI still needs Python < 3.10, and master carries Peter's
July 2026 fix, unreleased.

The model
---------
:class:`PhoenixStar` is ``StarfishStellarModel``'s twin. Its parameters, with
the legacy labels and ``lims`` as flat priors:

========================  ============  ==================
twin                      legacy        prior
========================  ============  ==================
``log_luminosity``        ``log(L)``    U(-1, 1), log10 L_sun
``teff``                  ``T_eff``     U(6000, 8000) K
``logg``                  ``log(g)``    U(4, 5)
``a_v``                   ``A(V)``      U(0, 3) mag
``r_v``                   ``R(V)``      U(2.9, 4.0)
========================  ============  ==================

``feh`` is a **buffer** at 0.0, because the legacy fixes [Fe/H] at 0.0.
``model.promote_buffer("feh", prior=st.uniform(0.0, 0.5))`` frees it without
touching the evaluation code. ``distance`` is a buffer at 100 pc, the legacy's
default, which its ``__main__`` never changes.

The ``"sed"`` channel is on the legacy grid ``10 ** np.linspace(-0.5, 1.0,
2000)`` micron. It is the emulator's shape on segment A, scaled to ``F_nu``
in Jy, placed on the grid by log-log interpolation and extrapolated beyond
5.5 micron as ``lambda**-2``, as the legacy does. CCM89 is then applied with
``a_v`` and ``r_v``. The ``"rvs"`` channel is the same on segment B, the
legacy ``specwaves``. See :mod:`.emulator` for the arithmetic.

What changed in the translation
-------------------------------
* **The inverse-square law.** The twin computes ``F_nu = shape x f_bol(1
  L_sun, 1 pc) x L / d**2 x 3.34e5 lambda_A**2``. The legacy divides by ``4 pi
  d**2`` *again* after that (line 53), although ``f_bol(1 L_sun, 1 pc)``
  already carries the ``4 pi``, so its fluxes are ``4 pi`` too faint. The twin
  does not reproduce that. A luminosity fitted with the legacy script is
  therefore ``4 pi`` times the twin's for the same photometry.
* **Vacuum wavelengths.** PHOENIX HiRes is tabulated in vacuum. The legacy's
  Starfish grid converted to air and the twin does not. The synthetic data
  come from the same emulator, so this is self-consistent.
* **CCM89 in the model.** The legacy calls ``extinction.ccm89``. The twin
  implements CCM89's ``a(x)``, ``b(x)`` itself (:func:`.emulator.ccm89_ab`), so
  the extinction is differentiable on every backend. A smoke row checks it
  against ``dust_extinction``'s ``CCM89`` to 1e-6.
* **The RVS spectrum** is resampled by :class:`~ampere.backends.reference.Resample`
  (legacy: ``spectres``). ``calUnc=0.0025`` becomes ``CalibrationScale(
  st.lognorm(0.0025, scale=1.0))``. ``scaleLengthPrior=0.01`` becomes the
  ``Matern32`` length-scale prior ``st.halfnorm(scale=0.01)`` micron, as in
  :mod:`examples.linear_sed`. The amplitude has no legacy number and gets a
  weakly informative ``st.halfnorm(scale=1e-3)`` Jy, a few per cent of the
  RVS continuum. ``--no-gp`` uses ``IndependentNoise`` instead.

Engines
-------
``--engine sbi|nuts|emcee`` and ``--backend reference|torch|jax``. ``sbi`` is
the legacy's ``SBI_SNPE`` (NPE, 10 000 simulations, one round) on the reference
backend, cached under ``~/.cache/ampere-phoenix/``. ``nuts`` runs on torch or
jax through the native twins (:mod:`.emulator_torch`, :mod:`.emulator_jax`)
and refuses ``reference`` by name. ``emcee`` runs on any backend.

::

    python -m examples.phoenix_star --engine sbi                 # pixi run -e sbi
    python -m examples.phoenix_star --engine nuts --backend torch  # pixi run -e torch

Coverage runs (Accept criterion)
--------------------------------
Pending: recorded here once run.
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
import time
from collections.abc import Mapping
from pathlib import Path
from typing import Any, ClassVar

import astropy.units as u
import numpy as np
import scipy.stats as st

from ampere.core import (
    DTYPE,
    ChannelRequirements,
    Dataset,
    DatasetCollection,
    FittingProblem,
    GaussianFamily,
    Instrument,
    Likelihood,
    Model,
    ModelResult,
    Spectrum,
)

from . import generators
from .emulator import (
    EMULATOR_FILE,
    NUMPY_OPS,
    ChannelPlan,
    Ops,
    _to_numpy,
    as_parameter,
    channel_plan,
    emulator_arrays,
    load_emulator,
    log_shape,
    plan_arrays,
    star_log_flux,
)

__all__ = [
    "BACKENDS",
    "CACHE_DIR",
    "FILTERS",
    "GRID",
    "PRIORS",
    "QUALIFIED_TRUTH",
    "PhoenixStar",
    "backend_module",
    "build_instruments",
    "build_model",
    "build_problem",
    "fit",
    "main",
    "recovers_truth",
    "report",
    "star_class",
]

BACKENDS = ("reference", "torch", "jax")

#: The legacy model grid (``phoenixstar.py`` line 166), micron.
GRID = 10 ** np.linspace(-0.5, 1.0, 2000)

#: The legacy's eight filters (line 145), all in the bundled library.
FILTERS = (
    "Gaia_BP",
    "Gaia_G",
    "Gaia_RP",
    "2MASS_J",
    "2MASS_H",
    "2MASS_Ks",
    "WISE_RSR_W1",
    "WISE_RSR_W2",
)

#: The legacy ``lims`` as flat priors (``scipy.stats.uniform(lo, width)``).
PRIORS: dict[str, Any] = {
    "log_luminosity": st.uniform(-1.0, 2.0),
    "teff": st.uniform(6000.0, 2000.0),
    "logg": st.uniform(4.0, 1.0),
    "a_v": st.uniform(0.0, 3.0),
    "r_v": st.uniform(2.9, 1.1),
}

CALIBRATION_PRIOR = st.lognorm(0.0025, scale=1.0)
GP_LENGTH_SCALE_PRIOR = st.halfnorm(scale=0.01)
GP_AMPLITUDE_PRIOR = st.halfnorm(scale=1e-3)

QUALIFIED_TRUTH: dict[str, float] = {f"model.{k}": v for k, v in generators.TRUTH.items()}

CACHE_DIR = Path.home() / ".cache" / "ampere-phoenix"

# Budgets: the legacy SBI budget; NUTS as the dispatch; emcee as linear_sed.
SBI_DEFAULT_BUDGET = 10_000
SBI_DEFAULT_DRAWS = 10_000
DEFAULT_WARMUP = 300
DEFAULT_DRAWS = 500
DEFAULT_CHAINS = 4
DEFAULT_WALKERS = 32
DEFAULT_STEPS = 3000
DEFAULT_BURN_IN = 1000


class PhoenixStar(Model):
    """A PHOENIX star: the emulator, a luminosity, a distance and CCM89.

    Parameters
    ----------
    wavelength
        The ``"sed"`` channel's grid, micron (default :data:`GRID`).
    rvs_wavelength
        The ``"rvs"`` channel's grid (default: the emulator's segment B).
    log_luminosity, teff, logg, a_v, r_v
        A prior to fit it, or a number to hold it fixed.
    feh, distance
        Buffers (see the module docstring); ``promote_buffer`` frees either.
    """

    OPS: ClassVar[Ops] = NUMPY_OPS
    CHANNELS: ClassVar[tuple[str, ...]] = ("sed", "rvs")
    AXIS: ClassVar[str] = "spectral_axis"

    def __init__(
        self,
        wavelength: Any = GRID,
        *,
        rvs_wavelength: Any = None,
        log_luminosity: Any = PRIORS["log_luminosity"],
        teff: Any = PRIORS["teff"],
        logg: Any = PRIORS["logg"],
        a_v: Any = PRIORS["a_v"],
        r_v: Any = PRIORS["r_v"],
        feh: float = 0.0,
        distance: float = 100.0,
        path: str | Path = EMULATOR_FILE,
    ) -> None:
        self.data = load_emulator(path)
        declared = {
            "log_luminosity": log_luminosity,
            "teff": teff,
            "logg": logg,
            "a_v": a_v,
            "r_v": r_v,
        }
        for name, given in declared.items():
            self.register_parameter(as_parameter(name, given))
        self.register_buffer("feh", np.asarray(feh, dtype=DTYPE))
        self.register_buffer("distance", np.asarray(distance, dtype=DTYPE), unit=u.pc)
        self.arrays = emulator_arrays(self.OPS, self.data)
        grids = {
            "sed": np.asarray(wavelength, dtype=DTYPE),
            "rvs": self.data.segment_b if rvs_wavelength is None else np.asarray(rvs_wavelength),
        }
        self.plans: dict[str, ChannelPlan] = {}
        self.plan_arrays: dict[str, dict[str, Any]] = {}
        self.templates: dict[str, Spectrum] = {}
        self.grids: dict[str, Any] = {}
        for channel, grid in grids.items():
            self._adopt(channel, grid)

    def _adopt(self, channel: str, grid: np.ndarray) -> None:
        plan = channel_plan(self.data, channel, grid)
        self.plans[channel] = plan
        self.plan_arrays[channel] = plan_arrays(self.OPS, plan)
        # Built here, outside any trace: a jax array made inside a jitted
        # function is a tracer, and the photometry step reads the grid in numpy.
        self.grids[channel] = self.OPS.asarray(plan.wavelength)
        self.templates[channel] = Spectrum(
            plan.wavelength * u.um, np.zeros(plan.wavelength.size, dtype=DTYPE), unit=u.Jy
        )

    def compile_for(self, requirements: Mapping[str, ChannelRequirements]) -> Model:
        """Adopt the negotiated grid for each channel, once."""
        for channel in self.CHANNELS:
            asked = requirements.get(channel)
            if asked is None or self.AXIS not in asked:
                continue
            self._adopt(channel, np.asarray(asked[self.AXIS].coordinates(), dtype=DTYPE))
        return self

    # -- the native surface (found by hasattr on torch and jax) ------------

    def grid(self, channel: str) -> Any:
        """*channel*'s wavelength grid, micron, as a backend array."""
        return self.grids[channel]

    def flux(self, channel: str, values: Mapping[str, Any] | None = None) -> Any:
        """*channel*'s ``F_nu`` in Jy, as a backend array."""
        return 10.0 ** self._log_flux(channel, self.context(values))

    def _log_flux(self, channel: str, ctx: Mapping[str, Any]) -> Any:
        shape = log_shape(self.OPS, self.arrays, ctx["teff"], ctx["logg"], ctx["feh"])
        return star_log_flux(
            self.OPS,
            self.plans[channel],
            self.plan_arrays[channel],
            shape,
            ctx["log_luminosity"],
            ctx["distance"],
            ctx["a_v"],
            ctx["r_v"],
        )

    # -- the contract surface ------------------------------------------------

    def evaluate(self, **values: Any) -> ModelResult:
        ctx = self.context(values)
        shape = log_shape(self.OPS, self.arrays, ctx["teff"], ctx["logg"], ctx["feh"])
        emitted = {}
        for channel in self.CHANNELS:
            log_flux = star_log_flux(
                self.OPS,
                self.plans[channel],
                self.plan_arrays[channel],
                shape,
                ctx["log_luminosity"],
                ctx["distance"],
                ctx["a_v"],
                ctx["r_v"],
            )
            flux = np.asarray(_to_numpy(10.0**log_flux), dtype=DTYPE)
            emitted[channel] = self.templates[channel].with_values(flux)
        return ModelResult(emitted)


def backend_module(backend: str) -> Any:
    """Where *backend*'s instruments and noise models live (jax set to float64)."""
    if backend == "reference":
        import ampere.backends.reference as module

        return module
    if backend == "torch":
        import ampere.backends.torch as module  # type: ignore[no-redef]

        return module
    if backend == "jax":
        import ampere.backends.jax as module  # type: ignore[no-redef]

        module.configure_x64()
        return module
    raise ValueError(f"unknown backend {backend!r}; the three are {BACKENDS}.")


def noise_module(backend: str) -> Any:
    """``IndependentNoise`` and ``GaussianProcessNoise``'s home for *backend*."""
    if backend == "reference":
        import ampere.core as module

        return module
    return backend_module(backend)


def star_class(backend: str) -> type[PhoenixStar]:
    """:class:`PhoenixStar` for *backend* (the native twins import torch/jax)."""
    if backend == "reference":
        return PhoenixStar
    backend_module(backend)
    if backend == "torch":
        from .emulator_torch import PhoenixStar as TorchStar

        return TorchStar
    from .emulator_jax import PhoenixStar as JaxStar

    return JaxStar


def build_model(backend: str = "reference") -> PhoenixStar:
    """The star with the legacy priors, on *backend*."""
    return star_class(backend)()


def build_instruments(backend: str = "reference") -> tuple[Instrument, Instrument]:
    """The eight-band catalogue on ``"sed"``, and the RVS spectrograph on ``"rvs"``."""
    module = backend_module(backend)
    data = load_emulator()
    catalogue = Instrument(
        [module.SyntheticPhotometry.from_library(list(FILTERS), GRID)],
        channel="sed",
        label="catalogue",
    )
    rvs = Instrument(
        [module.Resample(data.segment_b), module.CalibrationScale(CALIBRATION_PRIOR)],
        channel="rvs",
        label="rvs",
    )
    return catalogue, rvs


def _rvs_likelihood(backend: str, *, gp: bool) -> Likelihood:
    noise = noise_module(backend)
    if not gp:
        return Likelihood(GaussianFamily(), noise.IndependentNoise())
    module = noise_module(backend)
    kernel = module.Matern32(
        GP_AMPLITUDE_PRIOR,
        GP_LENGTH_SCALE_PRIOR,
        amplitude_unit=u.Jy,
        length_scale_unit=u.um,
        axes=("spectral_axis",),
    )
    return Likelihood(GaussianFamily(), noise.GaussianProcessNoise(kernel, module.DenseGP()))


def build_problem(
    backend: str = "reference", *, gp: bool = True, seed: int = generators.SEED
) -> FittingProblem:
    """The composed problem, entirely from *backend*'s pieces."""
    model = build_model(backend)
    catalogue, rvs = build_instruments(backend)
    # The observation is drawn on the reference backend, so every backend
    # fits the same numbers.
    observed_photometry, observed_rvs = generators.synthetic_data(
        build_model("reference"), *build_instruments("reference"), seed=seed
    )
    noise = noise_module(backend)
    datasets = DatasetCollection(
        {
            "catalogue": Dataset(
                observed_photometry,
                catalogue,
                likelihood=Likelihood(GaussianFamily(), noise.IndependentNoise()),
            ),
            "rvs": Dataset(observed_rvs, rvs, likelihood=_rvs_likelihood(backend, gp=gp)),
        }
    )
    return FittingProblem(model, datasets, seed=seed)


def _needs(extra: str, command: str) -> None:
    if importlib.util.find_spec(extra) is None:
        raise SystemExit(
            f"this engine needs {extra!r}, which is not installed in this environment. "
            f"Run `pixi install -e {command} --frozen` and "
            f"`pixi run -e {command} python -m examples.phoenix_star ...`."
        )


def fit(
    problem: FittingProblem,
    *,
    engine: str = "sbi",
    backend: str = "reference",
    budget: int | None = None,
    draws: int | None = None,
    warmup: int | None = None,
    chains: int | None = None,
    walkers: int | None = None,
    steps: int | None = None,
    burn_in: int | None = None,
    cache: bool = True,
    progress: bool = False,
) -> Any:
    """Sample *problem* with the named engine."""
    if engine == "emcee":
        from ampere.inference import EmceeEngine

        emcee_engine = EmceeEngine(problem, walkers=walkers or DEFAULT_WALKERS)
        return emcee_engine.run(
            steps or DEFAULT_STEPS,
            burn_in=DEFAULT_BURN_IN if burn_in is None else burn_in,
            progress=progress,
        )
    if engine == "nuts":
        if backend == "reference":
            raise SystemExit(
                "engine 'nuts' needs gradients: run it with --backend torch or --backend jax "
                "(the reference backend is not differentiable)."
            )
        from ampere.inference import NUTSEngine

        return NUTSEngine(problem).run(
            draws or DEFAULT_DRAWS,
            warmup=DEFAULT_WARMUP if warmup is None else warmup,
            chains=chains or DEFAULT_CHAINS,
            progress=progress,
        )
    if engine == "sbi":
        if backend != "reference":
            raise SystemExit("engine 'sbi' runs on the reference backend; drop --backend.")
        _needs("sbi", "sbi")
        from ampere.inference import SBIEngine
        from ampere.results import ArtefactStore

        store = ArtefactStore(CACHE_DIR) if cache else None
        sbi_engine = SBIEngine(problem, budget=budget or SBI_DEFAULT_BUDGET, rounds=1, cache=store)
        return sbi_engine.run(draws or SBI_DEFAULT_DRAWS, progress=progress)
    raise SystemExit(f"unknown engine {engine!r}; choose sbi, nuts or emcee.")


def recovers_truth(run: Any, *, level: float = 0.95) -> dict[str, bool]:
    """Whether each of the five parameters' central *level* interval covers its truth."""
    posterior = run["posterior"].dataset
    tail = (1.0 - level) / 2.0 * 100.0
    covered: dict[str, bool] = {}
    for name, truth in QUALIFIED_TRUTH.items():
        values = np.asarray(posterior[name], dtype=float).ravel()
        lower, upper = np.percentile(values, [tail, 100.0 - tail])
        covered[name] = bool(lower <= truth <= upper)
    return covered


def report(run: Any, truth: Mapping[str, float] = QUALIFIED_TRUTH) -> str:
    """A posterior summary, truth in brackets, 95 % coverage flagged, R-hat where defined."""
    attrs = run.attrs
    posterior = run["posterior"].dataset
    rhat: dict[str, float] = {}
    if posterior.sizes["chain"] > 1:
        try:
            import arviz as az

            summary = az.rhat(posterior)
            rhat = {name: float(summary[name]) for name in summary.data_vars}
        except Exception:
            rhat = {}
    tail = 2.5
    lines = [
        (
            f"{attrs['ampere_engine']} on {attrs['ampere_backend']}: "
            f"{posterior.sizes['chain']} chain(s) x {posterior.sizes['draw']} draw(s)"
        ),
        "  posterior (truth in brackets, 95 % coverage flagged):",
    ]
    for name in sorted(posterior.data_vars):
        values = np.asarray(posterior[name], dtype=float).ravel()
        lower, upper = np.percentile(values, [tail, 100.0 - tail])
        value = truth.get(name)
        flag = "" if value is None else ("  ok" if lower <= value <= upper else "  MISS")
        bracket = "" if value is None else f"  (truth {value:+.6g})"
        extra = f"  R-hat {rhat[name]:.3f}" if name in rhat else ""
        lines.append(
            f"    {name:42s} {values.mean():+.6g} +- {values.std():.3g}   "
            f"95%[{lower:+.6g}, {upper:+.6g}]{bracket}{flag}{extra}"
        )
    return "\n".join(lines)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--engine", default="sbi", choices=["sbi", "nuts", "emcee"])
    parser.add_argument("--backend", default=None, choices=list(BACKENDS))
    parser.add_argument("--no-gp", dest="gp", action="store_false", help="IndependentNoise")
    parser.add_argument("--seed", type=int, default=generators.SEED)
    parser.add_argument("--budget", type=int, default=None, help="sbi only")
    parser.add_argument("--draws", type=int, default=None, help="sbi/nuts")
    parser.add_argument("--no-cache", dest="cache", action="store_false", help="sbi only")
    parser.add_argument("--warmup", type=int, default=None, help="nuts only")
    parser.add_argument("--chains", type=int, default=None, help="nuts only")
    parser.add_argument("--walkers", type=int, default=None, help="emcee only")
    parser.add_argument("--steps", type=int, default=None, help="emcee only")
    parser.add_argument("--burn-in", type=int, default=None, help="emcee only")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(sys.argv[1:] if argv is None else argv)
    backend = args.backend or ("torch" if args.engine == "nuts" else "reference")
    if backend == "torch" or args.engine == "sbi":
        import torch

        # Four threads, as the BLAS pools in __main__: a shared machine.
        torch.set_num_threads(4)
    problem = build_problem(backend, gp=args.gp, seed=args.seed)
    print(f"free parameters: {problem.parameters.free_names}")
    started = time.perf_counter()
    run = fit(
        problem,
        engine=args.engine,
        backend=backend,
        budget=args.budget,
        draws=args.draws,
        warmup=args.warmup,
        chains=args.chains,
        walkers=args.walkers,
        steps=args.steps,
        burn_in=args.burn_in,
        cache=args.cache,
    )
    elapsed = time.perf_counter() - started
    print(report(run))
    print(f"  {elapsed:.1f} s wall clock")
    return 0
