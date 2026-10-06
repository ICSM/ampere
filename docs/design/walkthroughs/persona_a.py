"""Persona A walkthrough: photometry of one object, a modified blackbody, emcee, the table and plots.

Synthetic photometry (a truth model plus 8 % noise) so the numbers are honest; the
measurement is friction and time, not science. Every step timed; every refusal caught.
"""
import os as _os
OUT = _os.environ.get("AMPERE_WALK_OUT", "/tmp/ampere-walkthroughs"); _os.makedirs(OUT, exist_ok=True)
import time, traceback, sys, os
T0 = time.perf_counter()
def stamp(label):
    print(f"[{time.perf_counter()-T0:7.2f} s] {label}", flush=True)

import numpy as np, astropy.units as u, scipy.stats as st
stamp("numpy/astropy/scipy imported")
import ampere
stamp("import ampere")
from ampere.core import (Dataset, FittingProblem, GaussianFamily, IndependentNoise, Instrument,
                         Likelihood, Parameter, PhotometricPoints, GaussianProcessNoise, Matern32, DenseGP)
stamp("import ampere.core names")
from ampere.backends.reference import ModifiedBlackBody, SyntheticPhotometry
stamp("import ampere.backends.reference")
from ampere.inference import EmceeEngine, DynestyEngine
stamp("import ampere.inference")
from ampere.results import summary, plot_corner, plot_posterior_predictive, plot_residuals, add_posterior_predictive, to_netcdf
stamp("import ampere.results")

# --- the user's data: filter names and fluxes. Where do the wavelengths come from? ---
filters = ["2MASS_Ks", "WISE_RSR_W1", "WISE_RSR_W2", "WISE_RSR_W3", "WISE_RSR_W4",
           "SPITZER_MIPS_24", "SPITZER_MIPS_70", "HERSCHEL_PACS_100", "HERSCHEL_PACS_160"]
grid = np.geomspace(1.5, 250.0, 600)
step = SyntheticPhotometry.from_library(filters, grid)
stamp("SyntheticPhotometry.from_library built")
# FRICTION PROBE 1: does the library tell the user the effective wavelengths?
try:
    import pyphot
    lib = pyphot.get_library(os.path.join(os.path.dirname(ampere.__file__), "ampere_allfilters.hd5"))
    eff = np.array([float(lib[f].lpivot.to("micron").magnitude) if hasattr(lib[f].lpivot, "to") else float(lib[f].lpivot) for f in filters])
    print("  pivot wavelengths from pyphot (um):", np.round(eff, 2))
except Exception as e:
    print("  could not get pivot wavelengths from pyphot:", repr(e)); eff = None

instrument = Instrument([step], channel="default", label="phot")
truth_model = ModifiedBlackBody(grid, temperature=st.loguniform(30, 1500), beta=st.uniform(0.5, 2.0),
                                scale=st.loguniform(1e-3, 1e3), reference_wavelength=100.0)
stamp("model declared")
from ampere.core import negotiate
clean = instrument(truth_model.compile_for(negotiate([instrument]))(temperature=180.0, beta=1.6, scale=5.0))
print("  the step's own output wavelengths (um):", np.round(clean.spectral_axis.values, 2))
rng = np.random.default_rng(1)
sigma = 0.08 * np.abs(clean.values)
observed = PhotometricPoints(clean.filters, clean.spectral_axis.values * u.um,
                             (clean.values + rng.normal(0, sigma)) * u.Jy, uncertainty=sigma * u.Jy)
stamp("observed PhotometricPoints built (from the step's own wavelengths)")

# FRICTION PROBE 2: observed built with the user's own wavelengths (pyphot pivot) -- accepted?
if eff is not None:
    try:
        observed_user = PhotometricPoints(filters, eff * u.um, observed.values * u.Jy, uncertainty=sigma * u.Jy)
        like = Likelihood(GaussianFamily(), IndependentNoise())
        p = FittingProblem(truth_model, {"phot": Dataset(observed_user, instrument, likelihood=like)}, seed=1)
        p.evaluate(temperature=180.0, beta=1.6, scale=5.0)
        print("  PROBE 2: user-supplied pivot wavelengths ACCEPTED")
    except Exception as e:
        print("  PROBE 2: user-supplied pivot wavelengths REFUSED:", type(e).__name__, str(e)[:400])

# FRICTION PROBE 3: fluxes in mJy -- what happens?
try:
    observed_mjy = PhotometricPoints(clean.filters, clean.spectral_axis.values * u.um,
                                     (observed.values * 1000) * u.mJy, uncertainty=sigma * 1000 * u.mJy)
    like = Likelihood(GaussianFamily(), IndependentNoise())
    p = FittingProblem(truth_model, {"phot": Dataset(observed_mjy, instrument, likelihood=like)}, seed=1)
    lp = p.evaluate(temperature=180.0, beta=1.6, scale=5.0)
    print("  PROBE 3: mJy observation accepted; log_prob =", lp.log_prob if hasattr(lp, 'log_prob') else lp)
except Exception as e:
    print("  PROBE 3: mJy observation REFUSED:", type(e).__name__, str(e)[:400])

# FRICTION PROBE 4: wavelengths in Angstrom
try:
    observed_aa = PhotometricPoints(clean.filters, clean.spectral_axis.values * 1e4 * u.AA,
                                    observed.values * u.Jy, uncertainty=sigma * u.Jy)
    p = FittingProblem(truth_model, {"phot": Dataset(observed_aa, instrument, likelihood=like)}, seed=1)
    p.evaluate(temperature=180.0, beta=1.6, scale=5.0)
    print("  PROBE 4: Angstrom axis accepted")
except Exception as e:
    print("  PROBE 4: Angstrom axis REFUSED:", type(e).__name__, str(e)[:400])

# --- the fit ---
model = ModifiedBlackBody(grid, temperature=st.loguniform(30, 1500), beta=st.uniform(0.5, 2.0),
                          scale=st.loguniform(1e-3, 1e3), reference_wavelength=100.0)
likelihood = Likelihood(GaussianFamily(), IndependentNoise())
problem = FittingProblem(model, {"phot": Dataset(observed, instrument, likelihood=likelihood)}, seed=20261006)
stamp("FittingProblem built")
t = time.perf_counter()
run = EmceeEngine(problem, walkers=24).run(1500, burn_in=500)
stamp(f"emcee 24 walkers x 1500 steps done ({time.perf_counter()-t:.1f} s)")
table = summary(run)
print(table.to_string())
stamp("summary table")
# convergence: what does the user see?
print("  r_hat max:", float(table["r_hat"].max()), " ess_bulk min:", float(table["ess_bulk"].min()))
out = OUT
for name, fn in [("corner", lambda: plot_corner(run)),
                 ("ppc", lambda: plot_posterior_predictive(add_posterior_predictive(run, problem))),
                 ("residuals", lambda: plot_residuals(run))]:
    try:
        t = time.perf_counter(); fig = fn()
        (fig[0] if isinstance(fig, list) else fig).savefig(f"{out}/persona_a_{name}.png", dpi=80)
        stamp(f"plot {name} ({time.perf_counter()-t:.1f} s)")
    except Exception as e:
        print(f"  plot {name} FAILED: {type(e).__name__}: {str(e)[:300]}")
try:
    t = time.perf_counter(); path = to_netcdf(run, f"{out}/persona_a_run.nc"); stamp(f"to_netcdf ({time.perf_counter()-t:.1f} s): {path}")
except Exception as e:
    print("  to_netcdf FAILED:", type(e).__name__, str(e)[:300])
# evidence
try:
    t = time.perf_counter(); nested = DynestyEngine(problem).run()
    stamp(f"dynesty done ({time.perf_counter()-t:.1f} s): lnZ = {nested.attrs.get('ampere_log_evidence')} +- {nested.attrs.get('ampere_log_evidence_err')}")
except Exception as e:
    print("  dynesty FAILED:", type(e).__name__, str(e)[:300])
# the flexible likelihood on photometry: does it refuse, warn, or run?
try:
    kernel = Matern32(st.halfnorm(scale=0.5), st.halfnorm(scale=20.0), amplitude_unit=u.Jy, length_scale_unit=u.um, axes=("spectral_axis",))
    flexible = Likelihood(GaussianFamily(), GaussianProcessNoise(kernel, DenseGP()))
    problem_gp = FittingProblem(model, {"phot": Dataset(observed, instrument, likelihood=flexible)}, seed=20261006)
    t = time.perf_counter(); run_gp = EmceeEngine(problem_gp, walkers=24).run(600, burn_in=200)
    stamp(f"emcee with GP noise on 9 photometric points ({time.perf_counter()-t:.1f} s)")
    print(summary(run_gp).to_string())
except Exception as e:
    print("  GP on photometry FAILED:", type(e).__name__, str(e)[:400])
stamp("END")
