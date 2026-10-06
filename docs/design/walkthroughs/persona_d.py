"""Persona D: the contest OIFITS file through the reader, the binary fit on V2 and closure phases."""
import time, numpy as np, astropy.units as u, scipy.stats as st
T0 = time.perf_counter()
def stamp(l): print(f"[{time.perf_counter()-T0:7.2f} s] {l}", flush=True)
from ampere.core import Dataset, DatasetCollection, FittingProblem, Instrument, Likelihood
from ampere.interferometry import (read_oifits, FourierSample, SquaredAmplitude, ClosurePhase, Binary,
                                   GaussianFamily, VonMisesFamily)
from ampere.inference import EmceeEngine, optimise
from ampere.results import summary
stamp("imports")
t0 = time.perf_counter(); data = read_oifits("tests/data/contest-2008-binary.oifits"); stamp(f"read_oifits {time.perf_counter()-t0:.2f} s: {data.target}, V2 {data.squared_visibilities.values.shape}, T3 {data.closure_phases.values.shape}, channels {len(data.wavelengths)}")
FOV = 16.0 * u.mas
v2, t3 = data.squared_visibilities, data.closure_phases
vis2_instrument = Instrument([FourierSample.from_observed(v2, field_of_view=FOV), SquaredAmplitude(normalisation=1.0 * u.Jy)], channel="sky", label="vis2")
t3_instrument = Instrument([FourierSample.from_observed(t3, field_of_view=FOV), ClosurePhase()], channel="sky", label="t3")
model = Binary.on_field(FOV, 4, channels="sky", component_fwhm=0.9, separation=st.uniform(2.0, 8.0),
                        position_angle=st.uniform(0.0, 180.0), flux_ratio=st.uniform(0.02, 0.5), flux=1.0)
problem = FittingProblem(model, DatasetCollection({
    "vis2": Dataset(v2, vis2_instrument, likelihood=Likelihood(GaussianFamily()), label="vis2"),
    "t3": Dataset(t3, t3_instrument, likelihood=Likelihood(VonMisesFamily()), label="t3")}), seed=20261006)
stamp(f"problem built; free parameters {problem.free_size}; log_prob at the published truth = {problem.log_prob({'model.separation': 5.0, 'model.position_angle': 30.0, 'model.flux_ratio': 1/8.9}):.1f}")
t0 = time.perf_counter(); opt = optimise(problem); stamp(f"optimise {time.perf_counter()-t0:.1f} s: {opt}")
t0 = time.perf_counter(); run = EmceeEngine(problem, walkers=24).run(800, burn_in=300, initial=opt)
tab = summary(run); print(tab.to_string()); stamp(f"emcee {time.perf_counter()-t0:.1f} s; r_hat max {tab['r_hat'].max():.2f}; evaluations/s ~ {24*800/(time.perf_counter()-t0):.0f}")
stamp("END")
