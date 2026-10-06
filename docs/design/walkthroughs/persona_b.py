"""Persona B: a real Spitzer IRS spectrum (CASSIS, PG 1011-040) read by hand, a power-law continuum,
the flexible likelihood with QuasisepGP, the diagnostics a user would reach for. Timed; refusals caught."""
import os as _os
OUT = _os.environ.get("AMPERE_WALK_OUT", "/tmp/ampere-walkthroughs"); _os.makedirs(OUT, exist_ok=True)
import time, numpy as np, astropy.units as u, scipy.stats as st
from astropy.io import fits
from astropy.table import Table
T0 = time.perf_counter()
def stamp(l): print(f"[{time.perf_counter()-T0:7.2f} s] {l}", flush=True)
from ampere.core import (Dataset, FittingProblem, GaussianFamily, IndependentNoise, Instrument, Likelihood,
                         Spectrum, GaussianProcessNoise, Matern32, QuasisepGP, DenseGP)
from ampere.backends.reference import PowerLaw, Resample
from ampere.inference import EmceeEngine, optimise
from ampere.results import summary, add_residuals, plot_residuals, residual_whiteness, add_posterior_predictive, plot_posterior_predictive, plot_gp_localisation, gp_localisation_score
stamp("imports")
# --- reading the file by hand (what a reader would replace) ---
with fits.open("examples/test_data/cassis_yaaar_spcfw_14191360t.fits") as h:
    hd = h[0].header; names = [hd[f"COL{i:02d}DEF"] for i in range(1, 16)] + ["DUMMY"]
    t = Table(h[0].data, names=names)
t = t[np.isfinite(t["flux"]) & np.isfinite(t["error (RMS+SYS)"]) & (t["error (RMS+SYS)"] > 0)]
t.sort("wavelength")
# IRS orders overlap: keep unique wavelengths (the quickstart's trick)
_, idx = np.unique(t["wavelength"], return_index=True); t = t[idx]
wl = np.asarray(t["wavelength"], float); fl = np.asarray(t["flux"], float); er = np.asarray(t["error (RMS+SYS)"], float)
print(f"  {wl.size} samples, {wl.min():.1f}-{wl.max():.1f} um, flux {np.median(fl):.3f} Jy median, median S/N {np.median(fl/er):.0f}")
observed = Spectrum(wl * u.um, fl * u.Jy, uncertainty=er * u.Jy)
stamp("file read by hand: 9 lines of astropy")
instrument = Instrument([Resample(wl)], channel="default", label="irs")
def model():
    return PowerLaw(np.geomspace(5.0, 40.0, 800), norm=st.loguniform(1e-3, 10.0), index=st.uniform(-3.0, 6.0))
independent = Likelihood(GaussianFamily(), IndependentNoise())
p0 = FittingProblem(model(), {"irs": Dataset(observed, instrument, likelihood=independent)}, seed=1)
t0 = time.perf_counter(); opt = optimise(p0); stamp(f"optimise (independent) {time.perf_counter()-t0:.1f} s: {opt}")
t0 = time.perf_counter(); run0 = EmceeEngine(p0, walkers=24).run(1500, burn_in=500, initial=opt)
tab = summary(run0); print(tab.to_string()); stamp(f"emcee independent {time.perf_counter()-t0:.1f} s; r_hat max {tab['r_hat'].max():.2f}")
# the flexible likelihood
kernel = Matern32(st.halfnorm(scale=0.05), st.halfnorm(scale=5.0), amplitude_unit=u.Jy, length_scale_unit=u.um, axes=("spectral_axis",))
for solver_name, solver in [("QuasisepGP", QuasisepGP()), ("DenseGP", DenseGP())]:
    try:
        flexible = Likelihood(GaussianFamily(), GaussianProcessNoise(kernel, solver))
        p1 = FittingProblem(model(), {"irs": Dataset(observed, instrument, likelihood=flexible)}, seed=1)
        t0 = time.perf_counter(); opt1 = optimise(p1); stamp(f"optimise ({solver_name}) {time.perf_counter()-t0:.1f} s: {opt1}")
        t0 = time.perf_counter(); run1 = EmceeEngine(p1, walkers=24).run(1500, burn_in=500, initial=opt1)
        tab = summary(run1); print(tab.to_string()); stamp(f"emcee {solver_name} {time.perf_counter()-t0:.1f} s; r_hat max {tab['r_hat'].max():.2f}")
        break
    except Exception as e:
        print(f"  {solver_name} FAILED: {type(e).__name__}: {str(e)[:400]}")
# diagnostics as a user would try them
for name, fn in [("add_residuals+plot_residuals", lambda: plot_residuals(add_residuals(run1, p1))),
                 ("residual_whiteness", lambda: print("   whiteness:", residual_whiteness(run1))),
                 ("gp_localisation_score", lambda: print("   localisation:", gp_localisation_score(run1, p1))),
                 ("plot_gp_localisation", lambda: plot_gp_localisation(run1, p1)),
                 ("plot_posterior_predictive", lambda: plot_posterior_predictive(add_posterior_predictive(run1, p1)))]:
    try:
        t0 = time.perf_counter(); fig = fn(); stamp(f"{name} ({time.perf_counter()-t0:.1f} s)")
        if fig is not None and hasattr(fig, "savefig"): fig.savefig(f"{OUT}/persona_b_{name.split('+')[-1]}.png", dpi=80)
    except Exception as e:
        print(f"  {name} FAILED: {type(e).__name__}: {str(e)[:400]}")
stamp("END")
