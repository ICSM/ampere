"""Persona A, second pass: (1) a clean fit with nothing shared; (2) the same after a refused
composition shared the instrument; (3) the optimiser warm start; (4) a longer budget; (5) pivot wavelengths."""
import time, os, numpy as np, astropy.units as u, scipy.stats as st
T0 = time.perf_counter()
def stamp(label): print(f"[{time.perf_counter()-T0:7.2f} s] {label}", flush=True)
import ampere
from ampere.core import Dataset, FittingProblem, GaussianFamily, IndependentNoise, Instrument, Likelihood, PhotometricPoints, negotiate
from ampere.backends.reference import ModifiedBlackBody, SyntheticPhotometry
from ampere.inference import EmceeEngine, DynestyEngine, optimise
from ampere.results import summary
filters = ["2MASS_Ks", "WISE_RSR_W1", "WISE_RSR_W2", "WISE_RSR_W3", "WISE_RSR_W4",
           "SPITZER_MIPS_24", "SPITZER_MIPS_70", "HERSCHEL_PACS_100", "HERSCHEL_PACS_160"]
grid = np.geomspace(1.5, 250.0, 600)
TRUTH = dict(temperature=180.0, beta=1.6, scale=5.0)
def fresh():
    inst = Instrument([SyntheticPhotometry.from_library(filters, grid)], channel="default", label="phot")
    mdl = ModifiedBlackBody(grid, temperature=st.loguniform(30, 1500), beta=st.uniform(0.5, 2.0),
                            scale=st.loguniform(1e-3, 1e3), reference_wavelength=100.0)
    return inst, mdl
inst0, truth_model = fresh()
clean = inst0(truth_model.compile_for(negotiate([inst0]))(**TRUTH))
rng = np.random.default_rng(1); sigma = 0.08 * np.abs(clean.values)
observed = PhotometricPoints(clean.filters, clean.spectral_axis.values * u.um,
                             (clean.values + rng.normal(0, sigma)) * u.Jy, uncertainty=sigma * u.Jy)
like = lambda: Likelihood(GaussianFamily(), IndependentNoise())
def lp_at_truth(problem):
    ev = problem.evaluate({f"model.{k}": v for k, v in TRUTH.items()})
    return float(problem.log_prob({f"model.{k}": v for k, v in TRUTH.items()}))

# (1) clean
inst, mdl = fresh()
problem = FittingProblem(mdl, {"phot": Dataset(observed, inst, likelihood=like())}, seed=20261006)
print("(1) clean: log_prob at truth =", lp_at_truth(problem))
t = time.perf_counter(); run = EmceeEngine(problem, walkers=24).run(1500, burn_in=500)
tab = summary(run); print(tab.to_string()); stamp(f"(1) emcee 24x1500 {time.perf_counter()-t:.1f} s; r_hat max {tab['r_hat'].max():.2f}")
t = time.perf_counter(); nested = DynestyEngine(problem).run()
stamp(f"(1) dynesty {time.perf_counter()-t:.1f} s: lnZ {nested.attrs['ampere_log_evidence']:.2f} +- {nested.attrs['ampere_log_evidence_err']:.2f}")

# (2) the instrument first used in a refused composition
inst, mdl = fresh()
try:
    bad = PhotometricPoints(clean.filters, clean.spectral_axis.values * 1e4 * u.AA, observed.values * u.Jy, uncertainty=sigma * u.Jy)
    pbad = FittingProblem(mdl, {"phot": Dataset(bad, inst, likelihood=like())}, seed=1); pbad.log_prob({f"model.{k}": v for k, v in TRUTH.items()})
except Exception as e:
    print("(2) refused as expected:", type(e).__name__)
problem2 = FittingProblem(mdl, {"phot": Dataset(observed, inst, likelihood=like())}, seed=20261006)
print("(2) after the refused composition, same instrument and model: log_prob at truth =", lp_at_truth(problem2))
inst, _ = fresh(); _, mdl = fresh()
# model fresh, instrument reused from a refused composition
try:
    pbad = FittingProblem(mdl, {"phot": Dataset(bad, inst, likelihood=like())}, seed=1); pbad.log_prob({f"model.{k}": v for k, v in TRUTH.items()})
except Exception as e: pass
_, mdl2 = fresh()
problem3 = FittingProblem(mdl2, {"phot": Dataset(observed, inst, likelihood=like())}, seed=20261006)
print("(2b) fresh model, reused instrument: log_prob at truth =", lp_at_truth(problem3))
inst2, _ = fresh()
problem4 = FittingProblem(mdl, {"phot": Dataset(observed, inst2, likelihood=like())}, seed=20261006)
print("(2c) reused model, fresh instrument: log_prob at truth =", lp_at_truth(problem4))

# (3) the optimiser warm start
inst, mdl = fresh()
problem = FittingProblem(mdl, {"phot": Dataset(observed, inst, likelihood=like())}, seed=20261006)
t = time.perf_counter(); opt = optimise(problem); stamp(f"(3) optimise {time.perf_counter()-t:.1f} s: {opt}")
try:
    t = time.perf_counter(); run = EmceeEngine(problem, walkers=24).run(1500, burn_in=500, initial=opt)
    tab = summary(run); print(tab.to_string()); stamp(f"(3) emcee from the optimum {time.perf_counter()-t:.1f} s; r_hat max {tab['r_hat'].max():.2f}")
except Exception as e:
    print("(3) initial=opt failed:", type(e).__name__, str(e)[:300])
    eng = EmceeEngine(problem, walkers=24)
    t = time.perf_counter(); run = eng.run(1500, burn_in=500, initial=eng.initial_positions(24, around=opt))
    tab = summary(run); print(tab.to_string()); stamp(f"(3b) emcee from initial_positions(around=opt) {time.perf_counter()-t:.1f} s; r_hat max {tab['r_hat'].max():.2f}")

# (4) a longer budget from the prior
inst, mdl = fresh()
problem = FittingProblem(mdl, {"phot": Dataset(observed, inst, likelihood=like())}, seed=20261006)
t = time.perf_counter(); run = EmceeEngine(problem, walkers=48).run(6000, burn_in=2000)
tab = summary(run); print(tab.to_string()); stamp(f"(4) emcee 48x6000 {time.perf_counter()-t:.1f} s; r_hat max {tab['r_hat'].max():.2f}; ess min {tab['ess_bulk'].min():.0f}")

# (5) pivot wavelengths from the library, accepted?
import pyphot
lib = pyphot.get_library(os.path.join(os.path.dirname(ampere.__file__), "ampere_allfilters.hd5"))
eff = []
for f in filters:
    lp = lib[f].lpivot
    eff.append(float(lp.to("micron").value) if hasattr(lp, "to") else float(lp))
eff = np.array(eff); print("(5) pivot (um):", np.round(eff, 2), " step's own:", np.round(clean.spectral_axis.values, 2))
inst, mdl = fresh()
try:
    obs_user = PhotometricPoints(filters, eff * u.um, observed.values * u.Jy, uncertainty=sigma * u.Jy)
    p = FittingProblem(mdl, {"phot": Dataset(obs_user, inst, likelihood=like())}, seed=1)
    print("(5) user pivot wavelengths accepted; log_prob at truth =", lp_at_truth(p))
except Exception as e:
    print("(5) user pivot wavelengths REFUSED:", type(e).__name__, str(e)[:500])
stamp("END")
