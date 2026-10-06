"""Persona E: twenty objects as a Population on the numpy path (emcee), then the same twenty one by one
with to_netcdf per object, as the catalogue loop a user would write today."""
import os as _os
OUT = _os.environ.get("AMPERE_WALK_OUT", "/tmp/ampere-walkthroughs"); _os.makedirs(OUT, exist_ok=True)
import time, os, numpy as np, astropy.units as u, scipy.stats as st
T0 = time.perf_counter()
def stamp(l): print(f"[{time.perf_counter()-T0:7.2f} s] {l}", flush=True)
from ampere.core import Dataset, FittingProblem, GaussianFamily, HierarchicalPrior, IndependentNoise, Instrument, Likelihood, Parameter, Population, Spectrum
from ampere.backends.reference import PowerLaw
from ampere.inference import EmceeEngine, optimise
from ampere.results import summary, to_netcdf
stamp("imports")
N = 20; GRID = np.array([1.0, 3.0, 9.0]); NOISE = 0.02
rng = np.random.default_rng(20260919); indices = rng.normal(-1.30, 0.35, N)
models, datasets, observed = {}, [], {}
for i, index in enumerate(indices):
    truth = 2.0 * GRID**index
    observed[i] = Spectrum(GRID * u.micron, (truth + rng.normal(0.0, NOISE, GRID.size)) * u.Jy, uncertainty=np.full(GRID.size, NOISE) * u.Jy)
    models[f"obj{i}"] = PowerLaw(GRID, norm=Parameter("norm", value=2.0, fixed=True), index=st.norm(-1.3, 1.0))
    datasets.append(Dataset(observed[i], Instrument([], channel="default", input_kind=Spectrum, label=f"scope{i}"),
                            Likelihood(GaussianFamily(), IndependentNoise()), model=f"obj{i}", label=f"d{i}"))
population = Population("objects", members=[Parameter("index", HierarchicalPrior("norm", {"loc": "mu", "scale": "sigma"}))],
                        hyperpriors=[Parameter("mu", st.norm(-1.0, 1.0)), Parameter("sigma", st.halfnorm(0.0, 1.0))],
                        over=[f"obj{i}" for i in range(N)])
problem = FittingProblem(models, datasets, populations=[population], seed=20260919)
stamp(f"population problem built: free_size {problem.free_size}")
t0 = time.perf_counter(); opt = optimise(problem); stamp(f"optimise {time.perf_counter()-t0:.1f} s")
walkers = 2 * problem.free_size + 4
t0 = time.perf_counter(); run = EmceeEngine(problem, walkers=walkers).run(3000, burn_in=1000, initial=opt)
tab = summary(run, var_names=["objects.mu", "objects.sigma"]); print(tab.to_string())
full = summary(run); stamp(f"emcee {walkers}x3000 {time.perf_counter()-t0:.1f} s; r_hat max over all {full['r_hat'].max():.2f}; truth mu -1.30 sigma 0.35")
# the catalogue loop a user writes today
out = f"{OUT}/loop"; os.makedirs(out, exist_ok=True)
t0 = time.perf_counter(); times = []
for i in range(N):
    ti = time.perf_counter()
    m = PowerLaw(GRID, norm=Parameter("norm", value=2.0, fixed=True), index=st.norm(-1.3, 1.0))
    p = FittingProblem(m, {"d": Dataset(observed[i], Instrument([], channel="default", input_kind=Spectrum, label="scope"), Likelihood(GaussianFamily(), IndependentNoise()))}, seed=i)
    r = EmceeEngine(p, walkers=12).run(600, burn_in=200, initial=optimise(p))
    to_netcdf(r, f"{out}/obj{i}.nc"); times.append(time.perf_counter() - ti)
stamp(f"loop over {N} objects: {time.perf_counter()-t0:.1f} s total, {np.median(times):.2f} s median per object; files {sum(os.path.getsize(f'{out}/obj{i}.nc') for i in range(N))/1e6:.1f} MB")
from ampere.results import from_netcdf
r0 = from_netcdf(f"{out}/obj0.nc"); print("  provenance attrs on a reloaded run:", [k for k in r0.attrs if k.startswith("ampere_")][:12])
stamp("END")
