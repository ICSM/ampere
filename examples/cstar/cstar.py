"""The v2 twin of ``examples/cstar_model_test_sbi_v2.py`` and its ``_embedding``
variant -- W6.13 (5), the carbon star on Hyperion.

The legacy scripts (unchanged, as is ``examples/cstar_data/``) fit an OGLE LMC
carbon-rich AGB star -- eleven photometric points and a two-chunk Spitzer IRS
spectrum -- with ``HyperionCStarRTModel`` (``ampere/legacy/models/Hyperion.py``):
a spherical power-law dust shell of amorphous carbon and silicon carbide round
a tabulated photosphere, computed by the Hyperion Monte Carlo radiative
transfer code, and inferred by sequential neural posterior estimation (two
rounds, 10 000 simulations, ``SBI_SNPE``). The ``_embedding`` variant adds a
fully connected embedding network, ``{"type": "FC", "num_hiddens": 100,
"n_layers": 3, "output_dim": 28}``, and an emcee block the twin does not
reproduce: one Hyperion model costs minutes of CPU, so a likelihood-based
sampler needing ~10^5 of them is not an example anyone can run; SBI, which
spends a fixed simulation budget up front, is the reason this model is an
SBI example at all.

The model
---------
:class:`CarbonStarShell` is a black-box :class:`~ampere.core.Model` on the
reference backend -- an external Monte Carlo simulator, not differentiable,
declaring its failures through ``simulator_failures=`` -- in the shape of
:mod:`examples.sbi.external_simulator`. Its seven free parameters are the
legacy's six with the legacy boxes and names, plus one abundance:

======================  =======================================  ==================
twin                    meaning                                  prior
======================  =======================================  ==================
``envelope_mass``       log10 shell dust mass, M_sun             U(-10, -6)
``envelope_rin``        log10 inner radius, stellar radii        U(-2, 2)
``envelope_rout``       log10 outer radius, stellar radii        U(2, 4)
``envelope_r0``         log10 density reference radius, R_star   U(-2, 2)
``stellar_mass``        M_sun, linear                            U(1, 3)
``stellar_luminosity``  L_sun, **linear** (legacy)               U(1000, 10000)
``sic_fraction``        SiC mass fraction (carbon: 1 - it)       U(0, 1)
======================  =======================================  ==================

All three radii are in **stellar radii**, as the legacy ``__call__`` applies
them (``10**x * R_star``); the legacy ``__init__`` comment calls
``envelope_rout`` "inner radii", which ``__call__`` does not do. The legacy
exposed two abundances, ``abundance_1`` (carbon) and ``abundance_2`` (SiC),
with a Dirichlet(1, 1) prior and a prior transform that normalised two gamma
quantiles. A Dirichlet(1, 1) over two fractions that sum to one *is* a
uniform on either of them, so ``sic_fraction ~ U(0, 1)`` is the legacy prior
with the redundant coordinate removed: the legacy's eight parameters were
seven degrees of freedom.

Fixed, as buffers, from the legacy ``__init__`` defaults: the photosphere
table (:func:`.generators.read_photosphere`), ``stellar_radius = 2.056479e13``
cm, ``stellar_temperature = 3000`` K (recorded only -- the spectrum is
tabulated, so Hyperion never uses it, exactly as in the legacy), ``distance
= 50`` kpc, the envelope power ``-2``, 251 radial cells, the SED range
0.1-200 micron at 101 points, the one 45 degree viewing angle, and the
photon numbers (``photons=``: ``"legacy"`` is the legacy 1e4 initial and
1e5 imaging / raytracing-sources / raytracing-dust photons with 3
temperature iterations; ``"quick"`` is 1e3 / 1e4 / 1e4 / 1e4 with 2).

:meth:`CarbonStarShell.evaluate` rebuilds the ``AnalyticalYSOModel`` exactly
as the legacy ``__call__`` does -- the star's spectrum, luminosity, mass and
radius; the power-law envelope's mass ``10**envelope_mass`` M_sun and radii;
the spherical-polar grid; raytracing; the modified random walk with gamma =
2; the peeled SED with ``track_origin="detailed"``; the photon numbers; the
convergence criterion (99th percentile, absolute 2, relative 1.1) -- with the
dust from :mod:`.dust`. It writes ``model.rtin`` in a **per-worker scratch
directory keyed on the process id**, runs Hyperion, reads the SED at the
45 degree view, ``aperture=-1``, the star's distance, ``component="total"``,
in Jy, flips it to ascending wavelength and interpolates onto the legacy grid
:data:`GRID` (``np.geomspace(0.2, 200, 1000)`` micron), deletes the output
(tens of MB each), and returns a :class:`~ampere.core.Spectrum` on the
``"sed"`` channel. A non-zero Hyperion exit, a missing output, or an error
from Hyperion's own set-up raises :class:`SimulatorFailed` with the tail of
Hyperion's log.

What changed in the translation
-------------------------------
* **No ``bhmie``.** The legacy ran the unpackaged ``bhmie`` Fortran program
  per call to make the dust; the twin computes Mie opacities with
  ``miepython`` once per process and mixes them per call (:mod:`.dust`, ruled
  2026-09-28 on ``docs/design/example_dependencies_memo.md`` §3). The phase
  function is Henyey-Greenstein at miepython's asymmetry parameter rather than
  ``bhmie``'s tabulated Mie matrix (see :mod:`.dust`).
* **Serial per model, parallel across models.** The legacy ran each model on
  70 MPI processes (``nproc=70``). The twin runs Hyperion serially
  (``mpi=False``, one process) and parallelises *across* simulations with
  ``SBIEngine(executor=ProcessExecutor(workers))``: that is what a simulation
  bank wants (no MPI start-up per model, no idle ranks while a model's
  photons run out), and an mpich inside a pixi environment under WSL2 is the
  brittle route. There is no ``--mpi`` switch.
* **Seven parameters, not eight** (the Dirichlet above).
* **1 000 posterior draws by default, not 10 000.** The legacy's
  ``nsamples_post=10000`` cost nothing beyond the network. ``SBIEngine``
  scores every *stored* draw with the true ``log_prior`` and
  ``log_likelihood`` through ``problem.evaluate`` (what makes importance
  reweighting and calibration possible later), and here each such evaluation
  is a Hyperion run, serial, in the driving process: 10 000 draws would be
  some ten hours at ``photons="quick"`` and two days at ``"legacy"``.
  ``--draws`` sets it; the training budget is unaffected.
* **The data** as v2 containers, the IRS chunks sorted (see
  :mod:`.generators`): the legacy stored them in file order, which is not
  monotonic.
* **Noise.** The photometry is ``IndependentNoise`` at the votable's own
  uncertainties; each IRS chunk carries a
  ``CalibrationScale(st.lognorm(0.05, scale=1.0))`` and, by default, a
  ``GaussianProcessNoise(Matern32)`` with ``st.halfnorm(scale=0.01)`` micron on
  the length scale (the legacy ``scaleLengthPrior`` default family) and
  ``st.halfnorm(scale=0.01)`` Jy on the amplitude (a few tens of per cent of
  the IRS flux); ``--no-gp`` drops the GP.
* **Hyperion 0.9.11 on NumPy 2** needs the one alias
  :func:`.dust.hyperion_numpy_alias` restores (``np.string_``); see
  :mod:`.dust`.
* **The legacy's own defects** are not carried over and not fixed there: the
  undefined ``temperature`` name in the legacy's no-spectrum branch (lines 445
  and 613) is unreachable here because the spectrum is always given.

Requirements
------------
Hyperion and miepython are **example-only requirements** (D12 (b)): the pixi
``hyperion`` feature and environment (``pixi install -e hyperion``;
conda-forge ``hyperion`` and ``hyperion-fortran`` 0.9.11 plus ``miepython``,
on top of the ``sbi`` environment), not a ``pyproject`` extra, and absent
from CI. The tests skip where ``find_spec`` finds neither.

::

    pixi run -e hyperion python -m examples.cstar                       # the data, legacy budget
    pixi run -e hyperion python -m examples.cstar --embedding           # the _embedding variant
    pixi run -e hyperion python -m examples.cstar --synthetic --photons quick \\
        --rounds 1 --simulations 12800 --workers 8                      # the coverage run

The default ``--simulations 10000 --rounds 2 --photons legacy`` is the legacy
budget, documented here and not run: at the measured cost below it is about
90 CPU-hours of Hyperion (20 000 simulations at ~16 s), some eleven hours on
eight workers. ``--cache DIR`` (default ``~/.cache/ampere-cstar``) is an
:class:`~ampere.results.ArtefactStore` handed to ``SBIEngine(cache=)``: a
rerun at the same settings restores the trained posterior instead of
simulating or training again (``ampere_sbi_cache_hit`` in the run's attrs),
and :meth:`~ampere.inference.SBIEngine.calibrate` on the same engine reuses the
trained posterior without retraining. The simulated pairs themselves are
**not** written as a netCDF training set (``training_set=``), which the item
asked for: ``ampere.results.training``'s writer stores every
``extra_coords`` entry as float64, and a ``PhotometricPoints`` observation's
filter names are strings, so the first write fails (``could not convert
string to float: 'MCPS_B'``). That is a bug under ``ampere/``, outside this
example's scope, recorded in the W6.13 (5) report; ``training_set=`` joins the
call once it is fixed.

``--synthetic`` replaces the observed fluxes by **one simulation** at
:data:`.generators.SYNTHETIC_TRUTH` (the legacy defaults), pushed through the
instruments with noise at the data's own uncertainties. It uses the same
``--photons`` preset as the bank, so the "truth" carries the same Monte Carlo
noise every simulation does. That is the point of doing this by SBI: the
network learns the simulator's own scatter as part of the noise, whereas a
Gaussian likelihood on a noisy simulator would treat its Monte Carlo noise as
signal.

Measured cost and the coverage run (Accept criterion)
-----------------------------------------------------
One simulation at :data:`.generators.SYNTHETIC_TRUTH`, serially, in the
``hyperion`` environment on the 16-core development machine (2026-09-30; the
once-per-process Mie tables, 2.0 s, excluded): **16.1 s at
``photons="legacy"``, 3.7 s at ``photons="quick"``** -- of which about 2.2 s
is Hyperion's own mean-opacity and LTE-emissivity tabulation for the mixed
dust, the same at either preset. Pooled, 32 quick simulations on eight
workers took 17.9 s of wall clock, **4.5 s per simulation per worker**
(process start-up and the per-process Mie tables included), with no
failures.

The coverage budget is the largest that fits two hours of wall clock on
eight workers at the *pooled* cost, ``8 x 7200 / 4.5 = 12 800`` simulations
(the serial 3.7 s would give 15 500, which the pool does not achieve)::

    pixi run -e hyperion python -m examples.cstar --synthetic --photons quick \
        --rounds 1 --simulations 12800 --workers 8

RESULT-PENDING
"""

from __future__ import annotations

import argparse
import atexit
import importlib.util
import os
import shutil
import sys
import tempfile
import time
from pathlib import Path
from typing import Any

import astropy.units as u
import numpy as np
import scipy.stats as st

from ampere.core import (
    Dataset,
    DatasetCollection,
    DenseGP,
    FittingProblem,
    GaussianFamily,
    GaussianProcessNoise,
    IndependentNoise,
    Instrument,
    Likelihood,
    Matern32,
    Model,
    ModelResult,
    Parameter,
    PhotometricPoints,
    ProcessExecutor,
    Spectrum,
)
from ampere.backends.reference import CalibrationScale, Resample, SyntheticPhotometry

from . import dust, generators

__all__ = [
    "CALIBRATION_PRIOR",
    "DEFAULT_CACHE",
    "EMBEDDING",
    "GRID",
    "PHOTONS",
    "PRIORS",
    "CarbonStarShell",
    "SimulatorFailed",
    "build_instruments",
    "build_model",
    "build_problem",
    "fit",
    "main",
    "recovers_truth",
    "report",
    "working_directory",
]

#: The legacy model grid, ``np.logspace(log10(0.2), log10(200), 1000)`` micron.
GRID = np.geomspace(0.2, 200.0, 1000)

#: The legacy ``lims`` as flat priors (``scipy.stats.uniform(loc, width)``).
PRIORS: dict[str, Any] = {
    "envelope_mass": st.uniform(-10.0, 4.0),
    "envelope_rin": st.uniform(-2.0, 4.0),
    "envelope_rout": st.uniform(2.0, 2.0),
    "envelope_r0": st.uniform(-2.0, 4.0),
    "stellar_mass": st.uniform(1.0, 2.0),
    "stellar_luminosity": st.uniform(1.0e3, 9.0e3),
    "sic_fraction": st.uniform(0.0, 1.0),
}

#: Photon numbers ``(initial, imaging, raytracing sources, raytracing dust,
#: temperature iterations)``.
PHOTONS: dict[str, tuple[float, float, float, float, int]] = {
    "legacy": (1e4, 1e5, 1e5, 1e5, 3),
    "quick": (1e3, 1e4, 1e4, 1e4, 2),
}

#: The IRS chunks' calibration prior, 5 % lognormal.
CALIBRATION_PRIOR = st.lognorm(0.05, scale=1.0)
#: The GP length-scale prior, micron (the legacy ``scaleLengthPrior`` family).
GP_LENGTH_SCALE_PRIOR = st.halfnorm(scale=0.01)
#: The GP amplitude prior, Jy -- no legacy number; weakly informative.
GP_AMPLITUDE_PRIOR = st.halfnorm(scale=0.01)

#: The legacy ``_embedding`` variant's dict, verbatim (line 60 there).
EMBEDDING: dict[str, Any] = {"type": "FC", "num_hiddens": 100, "n_layers": 3, "output_dim": 28}

#: The legacy budget: two rounds of 10 000 simulations, 10 000 posterior draws.
LEGACY_SIMULATIONS = 10_000
LEGACY_ROUNDS = 2
LEGACY_DRAWS = 10_000
#: The twin's default posterior draws -- a tenth of the legacy's, because every
#: stored draw costs one Hyperion run here (see "What changed").
DEFAULT_DRAWS = 1_000

DEFAULT_CACHE = Path.home() / ".cache" / "ampere-cstar"

#: How much of Hyperion's log a failure keeps.
LOG_TAIL = 600

#: One scratch directory per worker process (module state, not pickled).
_WORKING_DIRECTORIES: dict[int, Path] = {}


class SimulatorFailed(RuntimeError):
    """Hyperion declined to produce an SED; declared through ``simulator_failures=``."""


@atexit.register
def _clean_working_directories() -> None:
    for path in _WORKING_DIRECTORIES.values():
        shutil.rmtree(path, ignore_errors=True)


def working_directory() -> Path:
    """This process's scratch directory, made once, on first use."""
    pid = os.getpid()
    existing = _WORKING_DIRECTORIES.get(pid)
    if existing is None:
        existing = Path(tempfile.mkdtemp(prefix=f"ampere-cstar-{pid}-"))
        _WORKING_DIRECTORIES[pid] = existing
    return existing


def _tail(path: Path) -> str:
    try:
        return path.read_text(errors="replace").strip()[-LOG_TAIL:] or "(empty)"
    except OSError:
        return "(no log)"


class CarbonStarShell(Model):
    """A dusty power-law shell round a carbon star, computed by Hyperion.

    Parameters
    ----------
    wavelength
        The output grid, micron (default the legacy :data:`GRID`).
    photons
        ``"legacy"`` or ``"quick"`` (:data:`PHOTONS`).
    carbon
        The amorphous-carbon constants, ``"rouleau91"`` (legacy) or ``"zubko96"``.
    channel
        The channel the emitted spectrum appears under.
    **priors
        Override any of :data:`PRIORS` (a prior, or a number to fix it).
    """

    def __init__(
        self,
        wavelength: Any = GRID,
        *,
        photons: str = "legacy",
        carbon: str = "rouleau91",
        channel: str = "sed",
        **priors: Any,
    ) -> None:
        unknown = set(priors) - set(PRIORS)
        if unknown:
            raise TypeError(f"CarbonStarShell has no parameter(s) {sorted(unknown)}")
        if photons not in PHOTONS:
            raise ValueError(f"unknown photons preset {photons!r}; choose {sorted(PHOTONS)}.")
        if carbon not in dust.CARBON:
            raise ValueError(f"unknown carbon {carbon!r}; choose {sorted(dust.CARBON)}.")
        self.channel = str(channel)
        self.photons = photons
        self.carbon = carbon
        grid = np.asarray(wavelength, dtype=float)
        for name, default in PRIORS.items():
            spec = priors.get(name, default)
            if hasattr(spec, "rvs"):
                self.register_parameter(Parameter(name, spec))
            else:
                self.register_parameter(Parameter(name, None, value=float(spec), fixed=True))
        nu, fnu = generators.read_photosphere()
        self.register_buffer("wavelength", grid, unit=u.um)
        self.register_buffer("photosphere_nu", nu, unit=u.Hz)
        self.register_buffer("photosphere_fnu", fnu)
        self.register_buffer("stellar_radius", np.asarray(2.056479e13), unit=u.cm)
        self.register_buffer("stellar_temperature", np.asarray(3000.0), unit=u.K)
        self.register_buffer("distance", np.asarray(50.0), unit=u.kpc)
        self.register_buffer("envelope_power", np.asarray(-2.0))
        self.register_buffer("grid_cells", np.asarray(251))
        self.register_buffer("photon_numbers", np.asarray(PHOTONS[photons], dtype=float))
        self.register_buffer("sed_range", np.asarray([0.1, 200.0, 101.0]), unit=u.um)
        self.register_buffer("viewing_angle", np.asarray(45.0), unit=u.deg)
        self._template = Spectrum(grid * u.um, np.zeros(grid.size), unit=u.Jy)

    def hyperion_model(self, ctx: Any) -> Any:
        """The ``AnalyticalYSOModel``, built exactly as the legacy ``__call__`` builds it."""
        dust.hyperion_numpy_alias()
        from hyperion.model import AnalyticalYSOModel

        radius = float(ctx["stellar_radius"])
        model = AnalyticalYSOModel()
        model.star.spectrum = (
            np.asarray(ctx["photosphere_nu"]),
            np.asarray(ctx["photosphere_fnu"]),
        )
        model.star.luminosity = (float(ctx["stellar_luminosity"]) * u.solLum).to(u.erg / u.s).value
        model.star.mass = (float(ctx["stellar_mass"]) * u.solMass).to(u.g).value
        model.star.radius = radius

        envelope = model.add_power_law_envelope()
        envelope.mass = ((10 ** float(ctx["envelope_mass"])) * u.solMass).to(u.g).value
        envelope.rmin = 10 ** float(ctx["envelope_rin"]) * radius
        envelope.rmax = 10 ** float(ctx["envelope_rout"]) * radius
        envelope.r_0 = 10 ** float(ctx["envelope_r0"]) * radius
        envelope.power = float(ctx["envelope_power"])
        envelope.dust = dust.build_dust(float(ctx["sic_fraction"]), self.carbon)

        model.set_spherical_polar_grid_auto(int(ctx["grid_cells"]), 1, 1)
        model.set_raytracing(True)
        model.set_mrw(True, gamma=2)
        angle = float(ctx["viewing_angle"])
        sed = model.add_peeled_images(sed=True, image=False)
        sed.set_viewing_angles([angle], [angle])
        sed.set_uncertainties(uncertainties=True)
        lmin, lmax, nl = (float(v) for v in ctx["sed_range"])
        sed.set_wavelength_range(int(nl), lmin, lmax)
        sed.set_track_origin("detailed")
        initial, imaging, sources, dust_photons, iterations = (
            float(v) for v in ctx["photon_numbers"]
        )
        model.set_n_photons(
            initial=initial,
            imaging=imaging,
            raytracing_sources=sources,
            raytracing_dust=dust_photons,
        )
        model.set_n_initial_iterations(int(iterations))
        model.set_convergence(True, percentile=99.0, absolute=2.0, relative=1.1)
        return model

    def evaluate(self, **values: Any) -> ModelResult:
        dust.hyperion_numpy_alias()
        from hyperion.model import ModelOutput

        ctx = self.context(values)
        scratch = working_directory()
        rtin = scratch / "model.rtin"
        rtout = scratch / "model.rtout"
        log = scratch / "hyperion.log"
        rtout.unlink(missing_ok=True)
        where = ", ".join(f"{name}={float(ctx[name]):.6g}" for name in PRIORS)
        try:
            model = self.hyperion_model(ctx)
            model.write(str(rtin), overwrite=True)
            model.run(str(rtout), logfile=str(log), mpi=False, n_processes=1, overwrite=True)
        except (Exception, SystemExit) as error:
            raise SimulatorFailed(
                f"Hyperion failed at {where}: {type(error).__name__}: {error}; log: {_tail(log)}"
            ) from None
        if not rtout.exists():
            raise SimulatorFailed(f"Hyperion wrote no output at {where}; log: {_tail(log)}")
        try:
            distance = (float(ctx["distance"]) * u.kpc).to(u.cm).value
            sed = ModelOutput(str(rtout)).get_sed(
                inclination=0, aperture=-1, distance=distance, component="total", units="Jy"
            )
            wav = np.asarray(sed.wav, dtype=float)
            val = np.asarray(sed.val, dtype=float).ravel()
        finally:
            rtout.unlink(missing_ok=True)
        grid = self._template.spectral_axis.values
        flux = np.interp(grid, np.flip(wav), np.flip(val))
        return ModelResult({self.channel: self._template.with_values(flux)})


def build_model(*, photons: str = "legacy", carbon: str = "rouleau91") -> CarbonStarShell:
    """The shell with the legacy ``lims`` as flat priors."""
    return CarbonStarShell(photons=photons, carbon=carbon)


def photometry_step(names: list[str]) -> SyntheticPhotometry:
    """The votable's eleven filters from the bundled library, on :data:`GRID`."""
    return SyntheticPhotometry.from_library(list(names), GRID)


def build_instruments(
    names: list[str], chunks: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]]
) -> dict[str, Instrument]:
    """The photometry and one resampled, calibrated instrument per IRS chunk."""
    instruments = {
        "photometry": Instrument([photometry_step(names)], channel="sed", label="photometry")
    }
    for chunk, (wavelength, _, _) in chunks.items():
        label = f"irs_{chunk}"
        instruments[label] = Instrument(
            [Resample(wavelength), CalibrationScale(CALIBRATION_PRIOR)],
            channel="sed",
            label=label,
        )
    return instruments


def _gp_likelihood(*, gp: bool) -> Likelihood:
    if not gp:
        return Likelihood(GaussianFamily(), IndependentNoise())
    kernel = Matern32(
        GP_AMPLITUDE_PRIOR,
        GP_LENGTH_SCALE_PRIOR,
        amplitude_unit=u.Jy,
        length_scale_unit=u.um,
        axes=("spectral_axis",),
    )
    return Likelihood(GaussianFamily(), GaussianProcessNoise(kernel, DenseGP()))


def observed_data(
    names: list[str],
    flux: np.ndarray,
    error: np.ndarray,
    chunks: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]],
    step: SyntheticPhotometry,
) -> dict[str, Any]:
    """The tracked data as v2 containers, labelled as :func:`build_instruments` labels."""
    observed: dict[str, Any] = {
        "photometry": PhotometricPoints(
            names, step.pivots() * u.um, flux * u.Jy, uncertainty=error * u.Jy
        )
    }
    for chunk, (wavelength, values, sigma) in chunks.items():
        observed[f"irs_{chunk}"] = Spectrum(
            wavelength * u.um, values * u.Jy, uncertainty=sigma * u.Jy
        )
    return observed


def build_problem(
    *,
    synthetic: bool = False,
    gp: bool = True,
    photons: str = "legacy",
    carbon: str = "rouleau91",
    seed: int = generators.SEED,
    validate: bool = False,
) -> FittingProblem:
    """The composed problem: the votable and the IRS chunks (or their synthetic twin).

    ``synthetic=True`` runs Hyperion once, at :data:`.generators.SYNTHETIC_TRUTH`.
    ``validate=False`` by default: ``FittingProblem``'s composition-time check
    evaluates the model once at the prior median, which here is a whole
    Hyperion run (and impossible where Hyperion is absent, which is where the
    data-loading smoke rows run). The check is still available as
    ``problem.validate()``, and the ``hyperion`` rows call it.
    """
    names, flux, error = generators.read_votable()
    chunks = generators.read_irs()
    instruments = build_instruments(names, chunks)
    step = instruments["photometry"].steps[0]
    observed = observed_data(names, flux, error, chunks, step)
    if synthetic:
        sigmas = {label: np.asarray(container.uncertainty) for label, container in observed.items()}
        observed = generators.synthetic_observations(
            build_model(photons=photons, carbon=carbon), instruments, sigmas, seed=seed
        )
    datasets = {
        label: Dataset(
            observed[label],
            instrument,
            likelihood=Likelihood(GaussianFamily(), IndependentNoise())
            if label == "photometry"
            else _gp_likelihood(gp=gp),
        )
        for label, instrument in instruments.items()
    }
    return FittingProblem(
        build_model(photons=photons, carbon=carbon),
        DatasetCollection(datasets),
        seed=seed,
        validate=validate,
        simulator_failures=(SimulatorFailed,),
    )


def _needs_sbi() -> None:
    if importlib.util.find_spec("sbi") is None:
        raise SystemExit(
            "engine='sbi' needs the 'sbi' extra, which is not installed in this environment. "
            "This example also needs Hyperion: install both with `pixi install -e hyperion` "
            "and run `pixi run -e hyperion python -m examples.cstar`."
        )


def _needs_hyperion() -> None:
    missing = [name for name in ("hyperion", "miepython") if importlib.util.find_spec(name) is None]
    if missing:
        raise SystemExit(
            f"examples.cstar needs {' and '.join(missing)}, which this environment lacks: they "
            f"are example-only requirements in the pixi 'hyperion' feature, not an ampere extra. "
            f"Install with `pixi install -e hyperion` and run `pixi run -e hyperion python -m "
            f"examples.cstar`."
        )


def fit(
    problem: FittingProblem,
    *,
    engine: str = "sbi",
    embedding: bool = False,
    simulations: int = LEGACY_SIMULATIONS,
    rounds: int = LEGACY_ROUNDS,
    draws: int = DEFAULT_DRAWS,
    workers: int = 4,
    cache: Path | str | None = DEFAULT_CACHE,
    training: dict[str, Any] | None = None,
    progress: bool = False,
) -> Any:
    """Fit *problem* by NPE, ``rounds`` rounds of ``simulations`` each, pooled over *workers*."""
    if engine != "sbi":
        raise SystemExit(f"examples.cstar fits with --engine sbi only, got {engine!r}.")
    _needs_sbi()
    from ampere.inference import SBIEngine
    from ampere.results import ArtefactStore

    store = None
    if cache is not None:
        cache = Path(cache).expanduser()
        cache.mkdir(parents=True, exist_ok=True)
        store = ArtefactStore(cache)
    run_engine = SBIEngine(
        problem,
        budget=int(simulations),
        rounds=int(rounds),
        embedding=EMBEDDING if embedding else None,
        executor=ProcessExecutor(int(workers)) if int(workers) > 1 else None,
        cache=store,
    )
    return run_engine.run(int(draws), training=training, progress=progress)


def recovers_truth(run: Any, *, level: float = 0.95) -> dict[str, bool]:
    """Whether each of the seven model parameters' central *level* interval covers its truth."""
    posterior = run["posterior"].dataset
    tail = (1.0 - level) / 2.0 * 100.0
    covered: dict[str, bool] = {}
    for name, truth in generators.SYNTHETIC_TRUTH.items():
        draws = np.asarray(posterior[f"model.{name}"], dtype=float).ravel()
        lower, upper = np.percentile(draws, [tail, 100.0 - tail])
        covered[name] = bool(lower <= truth <= upper)
    return covered


def report(run: Any, *, synthetic: bool = True) -> str:
    """A posterior summary: 95 % intervals, their width against the prior's, the truth."""
    attrs = run.attrs
    posterior = run["posterior"].dataset
    lines = [
        (
            f"{attrs['ampere_engine']} on {attrs['ampere_backend']}: "
            f"{posterior.sizes['chain']} chain(s) x {posterior.sizes['draw']} draw(s)"
        ),
        f"    {'parameter':36s} {'mean':>7s}  95 % interval{'':18s}width/prior  truth",
    ]
    for name in sorted(posterior.data_vars):
        values = np.asarray(posterior[name], dtype=float).ravel()
        lower, upper = np.percentile(values, [2.5, 97.5])
        short = name.removeprefix("model.")
        prior = PRIORS.get(short) if name.startswith("model.") else None
        ratio = ""
        if prior is not None:
            lo, hi = prior.ppf([0.025, 0.975])
            ratio = f"{(upper - lower) / (hi - lo):11.2f}"
        truth = generators.SYNTHETIC_TRUTH.get(short) if name.startswith("model.") else None
        flag = ""
        if synthetic and truth is not None:
            flag = f"  {truth:+.5g} {'ok' if lower <= truth <= upper else 'MISS'}"
        lines.append(
            f"    {name:36s} {values.mean():+.5g}  [{lower:+.5g}, {upper:+.5g}]{ratio:>13s}{flag}"
        )
    return "\n".join(lines)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m examples.cstar",
        description="The carbon star on Hyperion (W6.13 (5)); see the module docstring.",
    )
    parser.add_argument("--engine", default="sbi", choices=["sbi"])
    parser.add_argument("--embedding", action="store_true", help="the _embedding variant's net")
    parser.add_argument("--synthetic", action="store_true", help="fit one simulation at the truth")
    parser.add_argument("--no-gp", dest="gp", action="store_false", help="IndependentNoise on IRS")
    parser.add_argument("--photons", default="legacy", choices=sorted(PHOTONS))
    parser.add_argument("--carbon", default="rouleau91", choices=sorted(dust.CARBON))
    parser.add_argument("--simulations", type=int, default=LEGACY_SIMULATIONS)
    parser.add_argument("--rounds", type=int, default=LEGACY_ROUNDS)
    parser.add_argument("--draws", type=int, default=DEFAULT_DRAWS, help="each costs a run")
    parser.add_argument("--workers", type=int, default=4, help="Hyperion runs at once")
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--seed", type=int, default=generators.SEED)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(sys.argv[1:] if argv is None else argv)
    _needs_sbi()
    _needs_hyperion()
    started = time.perf_counter()
    problem = build_problem(
        synthetic=args.synthetic,
        gp=args.gp,
        photons=args.photons,
        carbon=args.carbon,
        seed=args.seed,
    )
    print(f"free parameters: {problem.parameters.free_names}")
    run = fit(
        problem,
        embedding=args.embedding,
        simulations=args.simulations,
        rounds=args.rounds,
        draws=args.draws,
        workers=args.workers,
        cache=args.cache,
    )
    elapsed = time.perf_counter() - started
    print(report(run, synthetic=args.synthetic))
    if args.synthetic:
        covered = recovers_truth(run)
        print(f"  95 % coverage of the seven truths: {sum(covered.values())}/7")
    print(f"  failed simulations: {dict(problem.failure_counts)}")
    print(f"  cache hit: {run.attrs.get('ampere_sbi_cache_hit')}")
    print(f"  {elapsed:.1f} s wall clock")
    return 0
