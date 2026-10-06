"""Persona B, second pass: separate the kernel-prior choice from the package's behaviour.
Five fits of the same IRS spectrum, the start (optimiser vs prior) and the kernel priors (wide vs M2-like) crossed."""
import time, os, numpy as np, astropy.units as u, scipy.stats as st
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from astropy.io import fits; from astropy.table import Table
OUT = os.environ.get("AMPERE_WALK_OUT", "/tmp/ampere-walkthroughs"); os.makedirs(OUT, exist_ok=True)
T0 = time.perf_counter()
def stamp(l): print(f"[{time.perf_counter()-T0:7.2f} s] {l}", flush=True)
from ampere.core import Dataset, FittingProblem, GaussianFamily, IndependentNoise, Instrument, Likelihood, Spectrum, GaussianProcessNoise, Matern32, QuasisepGP
from ampere.backends.reference import PowerLaw, Resample
from ampere.inference import EmceeEngine, optimise
from ampere.results import summary, add_residuals, plot_residuals, gp_localisation, plot_gp_localisation, gp_localisation_score
with fits.open("examples/test_data/cassis_yaaar_spcfw_14191360t.fits") as h:
    hd = h[0].header; names = [hd[f"COL{i:02d}DEF"] for i in range(1, 16)] + ["DUMMY"]; t = Table(h[0].data, names=names)
t = t[np.isfinite(t["flux"]) & np.isfinite(t["error (RMS+SYS)"]) & (t["error (RMS+SYS)"] > 0)]; t.sort("wavelength")
_, idx = np.unique(t["wavelength"], return_index=True); t = t[idx]
wl = np.asarray(t["wavelength"], float); fl = np.asarray(t["flux"], float); er = np.asarray(t["error (RMS+SYS)"], float)
observed = Spectrum(wl * u.um, fl * u.Jy, uncertainty=er * u.Jy)
instrument = Instrument([Resample(wl)], channel="default", label="irs")
def model(): return PowerLaw(np.geomspace(5.0, 40.0, 800), norm=st.loguniform(1e-3, 10.0), index=st.uniform(-3.0, 6.0))
def kernel(amp, ls): return Matern32(st.halfnorm(scale=amp), st.halfnorm(scale=ls), amplitude_unit=u.Jy, length_scale_unit=u.um, axes=("spectral_axis",))
PRIORS = {"wide": (0.05, 5.0), "m2like": (0.01, 2.0)}
runs = {}
def fit(name, likelihood, start):
    p = FittingProblem(model(), {"irs": Dataset(observed, instrument, likelihood=likelihood)}, seed=1)
    eng = EmceeEngine(p, walkers=24)
    if start == "optimiser":
        t0 = time.perf_counter(); opt = optimise(p); stamp(f"{name}: optimise {time.perf_counter()-t0:.1f} s -> {opt}")
        init = eng.initial_positions(24, around=opt); print(f"   initial ball, norm column: min {init[:,0].min():.5f} max {init[:,0].max():.5f} unique {np.unique(init[:,0]).size}")
    else:
        init = None
    t0 = time.perf_counter(); r = eng.run(1500, burn_in=500, initial=init)
    tab = summary(r); print(tab.to_string()); stamp(f"{name}: emcee {time.perf_counter()-t0:.1f} s; r_hat max {np.nanmax(tab['r_hat']):.2f}; ess min {tab['ess_bulk'].min():.0f}")
    runs[name] = (r, p); return r, p
fit("independent/optimiser", Likelihood(GaussianFamily(), IndependentNoise()), "optimiser")
for pname, (amp, ls) in PRIORS.items():
    for start in ("optimiser", "prior"):
        fit(f"gp-{pname}/{start}", Likelihood(GaussianFamily(), GaussianProcessNoise(kernel(amp, ls), QuasisepGP())), start)
# --- figures ---
from ampere.results import to_netcdf
for name, (r, p) in runs.items(): to_netcdf(r, f"{OUT}/persona_b2_{name.replace('/', '_')}.nc")
def predict(r, p):
    post = r.posterior; med = {v: float(np.median(post[v].values)) for v in post.data_vars}
    ev = p.evaluate(med)
    for attr in ("predictions", "predicted", "prediction", "model_result"):
        obj = getattr(ev, attr, None)
        if obj is not None:
            try:
                c = obj["irs"] if "irs" in obj else list(obj.values())[0]
                return med, np.asarray(c.values, float)
            except Exception: pass
    print("   Evaluation attributes:", [a for a in dir(ev) if not a.startswith("_")])
    return med, med["model.norm"] * wl ** med["model.index"]
fig, ax = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
ax[0].errorbar(wl, fl, er, fmt=".", ms=2, color="0.5", alpha=0.6, label="IRS, PG 1011-040")
for name, (r, p) in runs.items():
    med, pred = predict(r, p)
    ax[0].plot(wl, pred, lw=1.2, label=f"{name}: norm {med['model.norm']:.4f}, index {med['model.index']:.2f}")
    ax[1].plot(wl, (fl - pred) / er, lw=0.8, label=name)
ax[0].set_yscale("log"); ax[0].set_ylabel("flux (Jy)"); ax[0].legend(fontsize=7); ax[1].set_ylabel("(data - model)/sigma"); ax[1].set_xlabel("wavelength (um)"); ax[1].axhline(0, color="k", lw=0.5); ax[1].legend(fontsize=7)
fig.savefig(f"{OUT}/persona_b2_overlay.png", dpi=90); stamp("overlay saved")
r, p = runs["independent/optimiser"]; plot_residuals(add_residuals(r, p)).savefig(f"{OUT}/persona_b2_residuals_independent.png", dpi=90); stamp("independent residual plot saved")
for name in ("gp-wide/prior", "gp-m2like/prior", "gp-m2like/optimiser"):
    try:
        r, p = runs[name]; t0 = time.perf_counter(); tree = gp_localisation(r, p); f = plot_gp_localisation(tree)
        (f[0] if isinstance(f, list) else f).savefig(f"{OUT}/persona_b2_localisation_{name.replace('/', '_')}.png", dpi=90)
        stamp(f"localisation {name} ({time.perf_counter()-t0:.1f} s): score {gp_localisation_score(tree)}")
    except Exception as e:
        print(f"  localisation {name} FAILED: {type(e).__name__}: {str(e)[:300]}")
stamp("END")
