# Ampere v2 — Results, Provenance & Plotting Contract (W1.8)

Status: **frozen at `spec-v1.0`** (the tag created at the W1.13 merge,
2026-09; any later change to a §4 contract requires a decision-log entry in
`DEVELOPMENT_PLAN.md` in the same PR — ground rule 9; the §15 rulings of
2026-09-03 are all landed). Implements `DEVELOPMENT_PLAN.md` §4.6,
and discharges the obligations `inference.md` §18, `likelihoods.md` §16,
`parameters.md` §13, `results_schema.md` §16/§17.6, `lowering.md` §9.2 and
`diagnostics.md` §7 place on this item. Code: `ampere/results/`. Tests:
`tests/results/test_results.py`.

Every worked example below is executed as a doctest by
`tests/results/test_results_doctests.py`, so this document cannot drift from the
implementation without the suite going red. The examples share one namespace and
build on each other; read them in order.

Where this document and `DEVELOPMENT_PLAN.md` disagree, the plan wins and this
document is wrong — file it as a decision-log correction.

---

## 1. What this contract is for

Everything upstream of here produces one number and a handful of arrays per
draw. This contract decides what is *kept*, in what shape, and what may be
concluded from it later.

`DEVELOPMENT_PLAN.md` §4.6 makes the remit small and absolute: **ArviZ is the
single results format**, every engine emits it, all plotting is written once
against it, and serialisation is netCDF. That is a reaction to a specific
failure. Legacy ampere's `mixins.py` gave each sampler driver its own results
handling and its own plotting, so a fix to one never reached the others and
"sampler plotting parity" was a thing that had to be maintained by hand and
therefore was not. The cure is not a tidier monolith; it is that **every
function here takes the emitted run and nothing else**. A plotting function that
knew which sampler produced its input would be reintroducing the bug.

The second half of the remit is provenance, and it is not book-keeping.
`lowering.md` §9.2 establishes that a lowered model's random draws depend on the
*order* ampere emits its parameters in, so a seed alone does not identify a run;
`DEVELOPMENT_PLAN.md` §7 requires trained-artefact cache keys to hash the
model/prior/data spec so stale artefacts invalidate themselves. Both land here,
in the same hash.

It is **backend-neutral**, and its one dependency is optional: `ampere.results`
imports arviz lazily and never at module import (§10).

### Setup for the examples

```pycon
>>> import json
>>> import pathlib
>>> import tempfile
>>> from typing import Any, ClassVar
>>> import numpy as np
>>> import scipy.stats as st
>>> import astropy.units as u
>>> from ampere.core import (
...     Dataset, DatasetCollection, FittingProblem, Instrument, Model, ModelResult,
...     Parameter, Spectrum, Tie, Transformation,
... )
>>> from ampere.results import (
...     DrawRecorder, canonical_json, container_from_dict, container_to_dict,
...     describe_likelihood, emit, from_netcdf, gp_localisation_caveat, hash_of,
...     model_result_to_dict, plot_gp_localisation, problem_fingerprint,
...     provenance_attrs, to_netcdf, training_pair_to_dict,
... )

```

The problem every example below is emitted from is `inference.md` §15's — two
datasets, one two-channel model, an instrument chain each, and a calibration
nuisance tied across them. It is chosen because it is the smallest composition
that exercises everything this contract has to get right: a tie (one merged name
with two leaves), a per-dataset decomposition with more than one term, and two
observed containers with different coordinates and different units.

```pycon
>>> class Powerlaw(Model):
...     def __init__(self, **grids):
...         self._channels = tuple(grids)
...         for name, grid in grids.items():
...             self.register_buffer(name, np.asarray(grid, dtype=float), unit=u.micron)
...         self.register_parameter(Parameter("index", st.norm(-1.0, 0.5)))
...         self.register_parameter(Parameter("norm", st.loguniform(0.1, 10.0)))
...     def evaluate(self, **values):
...         ctx = self.context(values)
...         return ModelResult({
...             name: Spectrum(ctx[name] * u.micron, ctx["norm"] * ctx[name] ** ctx["index"] * u.Jy)
...             for name in self._channels
...         })
>>> class Calibrate(Transformation):
...     ACCEPTS: ClassVar[tuple[type, ...]] = (Spectrum,)
...     def __init__(self, **kwargs):
...         super().__init__(**kwargs)
...         self.register_parameter(Parameter("scale", st.lognorm(0.2)))
...     def apply(self, samples, values):
...         return samples.with_values(samples.values * self.context(values)["scale"])
>>> blue_grid, red_grid = np.array([1.0, 2.0, 4.0]), np.array([10.0, 20.0, 40.0])
>>> blue_data = Spectrum(
...     blue_grid * u.micron, [1.0, 0.5, 0.25] * u.Jy, uncertainty=[0.05] * 3 * u.Jy
... )
>>> red_data = Spectrum(
...     red_grid * u.micron, [0.1, 0.05, 0.025] * u.Jy, uncertainty=[0.005] * 3 * u.Jy
... )
>>> joint = FittingProblem(
...     Powerlaw(blue=blue_grid, red=red_grid),
...     DatasetCollection({
...         "blue": Dataset(blue_data, Instrument([Calibrate()], channel="blue")),
...         "red": Dataset(red_data, Instrument([Calibrate()], channel="red")),
...     }),
...     ties=[Tie(
...         "calibration",
...         ("blue.instrument.calibrate.scale", "red.instrument.calibrate.scale"),
...     )],
...     seed=20260902,
... )

```

## 2. The objects

| Object | Role |
|---|---|
| `emit` | Draws + one `Evaluation` per draw → the run, as a `DataTree` |
| `DrawRecorder` | The same, accumulated a draw at a time, for a driver that does not hold its chain as an array |
| `to_netcdf` / `from_netcdf` | The serialisation `DEVELOPMENT_PLAN.md` §4.6 names |
| `provenance_attrs` | The `ampere_*` attributes every run carries |
| `canonical_json`, `hash_of`, `hash_container`, `problem_fingerprint` | The hashing recipe (§9) |
| `describe_likelihood` | Family + noise + solver + kernel + censoring, which no other spec serialises (§12, Q8) |
| `container_to_dict` / `model_result_to_dict` and their inverses | `results_schema.md` §17.6's unclaimed serialisation (§11) |
| `add_posterior_predictive`, `add_residuals`, `gp_localisation`, `add_pointwise_log_likelihood` | Groups computed on demand, never stored by default (§7; the last of the four is §6's, added W2.8) |
| `write_training_set` / `append_training_set` / `read_training_set` | §11 layer 2's on-disk training set (added W2.8) |
| `plot_*` | The plotting surface (§8) |
| `ResultsError` | This contract's error; a `ContractError`, hence a `ValueError` |

### A note on `InferenceData`

`DEVELOPMENT_PLAN.md` §4.6 says "ArviZ `InferenceData` is the single results
format". ArviZ 1.0 retired the `InferenceData` *class* in favour of
`xarray.DataTree`. **The format is unchanged** — the same named groups, the same
netCDF layout, and `arviz.from_netcdf` still reads files written by ArviZ 0.x —
so the decision stands untouched and only the Python type moved. Everything here
returns and accepts a `DataTree`. §15's ruling request R3 asks Peter to ratify
the wording change rather than let two documents disagree quietly.

## 3. The simple path

An engine driver has a θ and wants the run recorded. That is two lines, and the
recorder evaluates the problem itself when it is not handed an `Evaluation`:

```pycon
>>> recorder = DrawRecorder(joint, chains=2)
>>> rng = np.random.default_rng(20260902)
>>> for chain in (0, 1):
...     for _ in range(4):
...         _ = recorder.record(joint.prior_transform(rng.random(joint.free_size)), chain=chain)
>>> run = recorder.emit(engine="emcee")
>>> sorted(run.children)
['constant_data', 'log_likelihood', 'observed_data', 'posterior', 'sample_stats']

```

A driver that already holds its chains as an array calls `emit` directly with
the array and the evaluations; `DrawRecorder` is a convenience over it, not a
second path.

**`evaluate`, not `log_prob`.** `inference.md` §10 shapes `Evaluation` to be
exactly a per-draw record — the split, the per-dataset terms, the failure — in
one pass. A driver that calls `log_prob` has thrown that away and cannot get it
back without re-evaluating, so this contract refuses the float rather than
silently storing less:

```pycon
>>> emit(joint, np.zeros((1, 3)), [joint.log_prob({"model.index": -1.0, "model.norm": 1.0, "calibration": 1.0})])
Traceback (most recent call last):
    ...
ampere.core.exceptions.ResultsError: every recorded draw must be a FittingProblem.evaluate()
result, got float. ...

```

## 4. The groups, and what is in each

| Group | Contents | Dims |
|---|---|---|
| `posterior` | one variable per **merged parameter name** | `(chain, draw)`, plus a named dimension for an array-valued parameter |
| `sample_stats` | `lp` (= `log_prob`), `log_prior`, `log_likelihood` (the scalar joint), `failed`, `failure_reason`, `failure_where` | `(chain, draw)` |
| `log_likelihood` | one variable per **dataset label** — the per-dataset decomposition (§5) | `(chain, draw)` |
| `observed_data` | one variable per dataset: the observed values | the dataset's own coordinate axis |
| `constant_data` | per dataset: uncertainties, mask, extra coordinates, and any axes that are not dimensions | as above |

The posterior is keyed by the merged name, which is the third of the three
reasons `inference.md` §4.5 gives for the nested merge topology: a merged name
is "a `Tie` site, a `free_labels()` entry, an ArviZ coordinate, a corner-plot
axis label, and a string a user types". A tie is one dimension with one name, in
the results as everywhere else:

```pycon
>>> sorted(run["posterior"].data_vars)
['calibration', 'model.index', 'model.norm']
>>> run["posterior"]["calibration"].dims, run["posterior"]["calibration"].shape
(('chain', 'draw'), (2, 4))
>>> sorted(run["log_likelihood"].data_vars)
['blue', 'red']

```

Observed data keep their own coordinates, scoped by dataset label so that two
spectra of different length can live in one file, with units recorded as CF-style
attributes:

```pycon
>>> run["observed_data"]["blue"].dims
('blue_spectral_axis',)
>>> run["observed_data"].coords["blue_spectral_axis"].values.tolist()
[1.0, 2.0, 4.0]
>>> run["observed_data"]["blue"].attrs["units"]
'Jy'
>>> sorted(run["constant_data"].data_vars)
['blue_uncertainty', 'red_uncertainty']

```

The rule for dimensions, stated once. A **gridded** kind (`Image`, `Cube`) takes
one dimension per axis. A **point** kind with a single axis takes that axis, so a
spectrum plots against wavelength with nothing further asked of the caller. A
point kind with several axes — a `VisibilitySet`'s *u* and *v*, which index
samples jointly rather than separately — takes one sample dimension `<label>
_index`, and its axes become ordinary variables on it. Complex values are split
into `<label>_real` and `<label>_imag` because netCDF has no complex type, and
named so that nothing mistakes one part for the whole.

### Array-valued parameters are never 10⁵ names

`likelihoods.md` §16(a) is explicit: "the latent block cannot be labelled with
`free_labels()` at scale: 10⁵ scalar names is the wrong ArviZ representation,
and a named dimension with a coordinate is the right one." An array-valued
parameter is therefore **one variable with a named dimension**, whatever its
size, and `free_labels()` is not what names it:

```pycon
>>> from ampere.core import GaussianProcessNoise, Likelihood, Matern32, PoissonFamily
>>> counts = Spectrum([1.0, 2.0, 3.0] * u.um, np.array([4.0, 7.0, 2.0]))
>>> class Rate(Model):
...     def __init__(self, grid):
...         self.register_buffer("grid", grid, unit=u.um)
...         self.register_parameter(Parameter("rate", st.loguniform(0.5, 50.0)))
...     def evaluate(self, **values):
...         ctx = self.context(values)
...         return Spectrum(ctx["grid"] * u.um, np.full(ctx["grid"].shape, ctx["rate"]))
>>> latent = FittingProblem(
...     Rate(np.array([1.0, 2.0, 3.0])),
...     [Dataset(
...         counts,
...         likelihood=Likelihood(PoissonFamily(), GaussianProcessNoise(Matern32(0.3, 1.0))),
...         label="counts",
...     )],
... )
>>> draws = np.zeros((1, 2, latent.free_size))
>>> draws[..., 0] = 3.0
>>> latent_run = emit(latent, draws, [[latent.evaluate(draws[0, d]) for d in (0, 1)]])
>>> sorted(latent_run["posterior"].data_vars)
['counts.latent.z', 'model.rate']
>>> latent_run["posterior"]["counts.latent.z"].dims
('chain', 'draw', 'counts.latent.z_dim_0')

```

A **plate** member takes the plate's own name as its dimension, so every member
of one plate shares one dimension and one coordinate — which is what makes
`hierarchical_population.md` §10.5's "`objects.theta` should be one dimension
with the dataset labels as its coordinate" expressible at all. And when the
plate is routed element by element, the coordinate is **read off the mapping's
own `Binding.index` entries** rather than guessed from the dataset ordering:

```pycon
>>> from ampere.core import HierarchicalPrior, ParameterSet, Plate, PlateBinding
>>> from ampere.results.emission import _dimension_names, _index_coordinate
>>> plate = Plate(
...     "objects", size=2,
...     hyperparameters=[Parameter("mu", st.norm(0.0, 5.0))],
...     members=[Parameter("theta", HierarchicalPrior("norm", {"loc": "mu"}))],
... )
>>> mapping = ParameterSet.merge(
...     {
...         "population": ParameterSet([], plates=[plate]),
...         "ngc1": ParameterSet([Parameter("cal", st.lognorm(0.1))]),
...         "ngc2": ParameterSet([Parameter("cal", st.lognorm(0.1))]),
...     },
...     plate_bindings=[
...         PlateBinding("population.objects.theta", "ngc1", "theta", 0),
...         PlateBinding("population.objects.theta", "ngc2", "theta", 1),
...     ],
... )
>>> member = mapping.merged["population.objects.theta"]
>>> _dimension_names(member)
('objects',)
>>> _index_coordinate(mapping.bindings, member)
['ngc1', 'ngc2']

```

`hierarchical_population.md` §10.2 asks for exactly this — "that ordering then
becomes part of the run's provenance, because it defines which posterior element
belongs to which object — W1.8 must record the dataset labels as the plate's
coordinate, not an integer range". Reading the bindings rather than the
collection order is what keeps it true when the wiring is *not* in dataset order.
Where nothing routes the block element by element, the honest label is an integer
range, and `emit(..., coords=...)` is how a caller supplies a better one.

## 5. "Per-sample `log_likelihood`" — the ambiguity, resolved

`DEVELOPMENT_PLAN.md` §4.6 asks that "every run stores per-sample
`log_likelihood` and `log_prior`". `likelihoods.md` §16(b) points out that the
sentence has two legitimate readings and that they differ *for this project*:

- the **scalar joint** log-likelihood per posterior draw, which is what design
  horizon (b)'s population importance reweighting needs, and what
  `FittingProblem.log_likelihood` returns;
- ArviZ's own `log_likelihood` group, which is conventionally
  **per-observation**, because that is what LOO and WAIC consume — and which
  `diagnostics.md` §3.3 reads §4.6 as promising.

**Ampere stores the scalar joint in `sample_stats`, and the per-*dataset*
decomposition in the `log_likelihood` group.** Neither reading is discarded and
neither is silently substituted for the other.

The per-dataset decomposition is the one that is always well defined. A GP
likelihood does not factorise over observations — the samples are not
independent, so the joint value is not a sum of per-point terms — whereas
datasets *are* conditionally independent given θ by `inference.md` §7's single
modelling assumption, which is the whole reason a joint fit is a sum. So the
split adds up exactly, and it does so on a stored run as it does in memory:

```pycon
>>> joint_ll = run["sample_stats"]["log_likelihood"].values
>>> summed = run["log_likelihood"]["blue"].values + run["log_likelihood"]["red"].values
>>> bool(np.allclose(summed, joint_ll))
True

```

The group says which decomposition it is, in its own attributes, so no consumer
has to *assume* ArviZ's convention applies:

```pycon
>>> run["log_likelihood"].attrs["ampere_decomposition"]
'per_dataset'

```

Practical consequence, stated plainly: `arviz.loo` over this group computes
**leave-one-dataset-out**, which is a meaningful quantity for a joint fit and a
useless one for a single-dataset fit of 10⁵ points. That is a reason to label the
group, which this contract does, and not a reason to withhold it.

### A prior-rejected draw stores NaN, and it matters

`inference.md` §10 returns `-inf` for `log_prob` **without evaluating the
model**, `NaN` for `log_likelihood`, and records no failure. §18(c) asks this
contract not to flatten that on the way in, and it does not:

```pycon
>>> outside = {"model.index": -1.0, "model.norm": 1e9, "calibration": 1.0}
>>> rejected = emit(joint, np.array([joint.parameters.pack(outside)]), [joint.evaluate(outside)])
>>> stats = rejected["sample_stats"]
>>> float(stats["log_prior"].values[0, 0]), bool(np.isnan(stats["log_likelihood"].values[0, 0]))
(-inf, True)
>>> bool(stats["failed"].values[0, 0])
False

```

"Not evaluated" and "impossible" are different statements. An importance
reweighting consumer must be able to tell them apart, and coercing NaN to `-inf`
destroys the distinction irreversibly — the point would then be indistinguishable
from one the likelihood genuinely scored as impossible. The per-dataset group
records the same thing the same way: NaN, not zero, because a dataset that was
never scored contributed nothing rather than contributing nothing *to the sum*.

## 6. The per-observation decomposition: named, and deliberately not computed

`likelihoods.md` §16(b) says that if this contract wants ArviZ's per-observation
convention "it must choose a decomposition and name it — the leave-one-out
conditional terms are the standard choice and are computable from the same
Cholesky this contract already forms". The decision is to **name it and defer
computing it**, for three reasons that are worth separating.

1. **It is a `likelihoods.md` API addition, not a results one.** For an
   independent-noise likelihood the per-observation terms are the family's own
   `log_prob` evaluated pointwise; for a GP likelihood they are the leave-one-out
   conditional terms `log N(y_i | μ_i^{-i}, σ_i^{2,-i})`, obtainable in closed
   form from the same Cholesky (`μ_i^{-i} = y_i - [K̃⁻¹y]_i / [K̃⁻¹]_{ii}`,
   `σ_i^{2,-i} = 1/[K̃⁻¹]_{ii}`). Neither is reachable from outside
   `Likelihood`, so emitting them means adding a method to a merged contract.
   This item may not do that unilaterally (`AGENTS.md` ground rule 9); §15's R2
   asks for it.
2. **Two different decompositions must not share a name.** Writing
   independent-noise pointwise terms and GP conditional terms into one group
   called `log_likelihood`, alongside the per-dataset terms, would give
   `arviz.loo` three different meanings depending on the fit. Each needs its own
   group and its own declared decomposition.
3. **The cost is real and the demand is not yet.** Per-observation terms are
   `N_draws × N_obs` per dataset; the per-dataset decomposition is
   `N_draws × N_datasets`. Nothing in Phase 1 consumes the former.

**Reserved now**, as the constant
`ampere.results.POINTWISE_LOG_LIKELIHOOD_GROUP` rather than in prose: the group
name `pointwise_log_likelihood`, with a required `ampere_decomposition`
attribute taking one of `"factorised"` (independent noise; exact) or
`"conditional_loo"` (GP; the leave-one-out conditionals named above). Written
only by an explicit call, never by default. A run emitted today is
forward-compatible with one emitted after it lands, because the group it would
occupy is empty rather than misused.

*(Landed at the freeze — ruled 2026-09-03, §15 R2, with the decision-log
entry: `Likelihood.pointwise_log_prob` now computes both decompositions
under exactly these names — the family's pointwise terms for independent
noise, `GPSolver.conditional_loo` for a GP — see `likelihoods.md` §8. The
not-stored-by-default rule is unchanged: emitting the group remains an
explicit call, and the emission helper itself is Phase 2's, beside the
engine drivers that produce runs worth decomposing.)*

*(Landed W2.8.* `ampere.results.add_pointwise_log_likelihood(tree, problem)`
is that helper, and it is the only thing that writes the group. It stores one
term per retained observation on the container's own coordinate axis, masked
samples as NaN; it declares the decomposition per variable as well as on the
group, since a joint fit that mixes an independent-noise dataset with a GP one
has two of them and no single honest answer; and where the solver does not
implement `conditional_loo` it **refuses by name rather than falling back to a
dense solve**, because at the sizes `QuasisepGP` exists for that substitution
is a different program and not a slower answer.
`ampere.results.pointwise_as_log_likelihood(tree)` is the one-line bridge to
`arviz.loo`, which reads the group *named* `log_likelihood` — where ampere
deliberately keeps the per-dataset split (§15 R6), so the swap has to be
explicit rather than discovered by experiment.*)

## 7. Groups a run does not store, and the rule for getting them

`diagnostics.md` §7 puts two questions to this contract.

**Are posterior-predictive replicates stored by default?** *No.* The cost is
`N_draws × N_obs` per dataset — for a 10⁵-point spectrum and 10⁴ draws, 8 GB of
float64 in a file whose whole point is being cheap enough to archive — and family
B needs them only when a check is actually asked for.
`ampere.results.add_posterior_predictive` computes them into the reserved group
`posterior_predictive`, drawing through `FittingProblem.simulate(observe=True)`
so the replicates come from the same `LikelihoodFamily.sample` the likelihood
scores with, never from a Gaussian assumption bolted on at the diagnostic layer
(`inference.md` §13). It takes its randomness from the named sub-stream
`"posterior_predictive"`, distinct from `"simulate"` on purpose: adding a
predictive check must not change an SBI budget's draws (`lowering.md` §9.2).

**Where do signed per-point residuals live?** In the reserved group `residuals`,
also on demand, computed by `add_residuals`. `diagnostics.md` §3.3 is precise
about why this is not free: §4.6's log-likelihood requirement makes the
posterior-predictive half cheap, but a Ljung-Box whiteness statistic needs the
**sign**, and for a Gaussian noise model `log_likelihood_i` recovers
`|residual_i|` and not its sign. The derivation is fixed here so no
implementation has to invent it:

> The standardised residual of dataset *d* at draw *k* is
> `r_dk = (observed_d − predicted_d(θ_k)) / σ_d`, with `predicted_d(θ_k)` the
> instrument chain's output — `Dataset.predict` on the model result at `θ_k` —
> and `σ_d` the observed container's own `masked_uncertainty()`. Masked samples
> are carried as NaN, never dropped, so the coordinate axis stays the
> container's own and the statistic decides for itself what to do with a gap.

Family C's input is not a group at all: `likelihoods.md` §16 fixes it as
`Likelihood.conditional`, which returns a signed mean and a variance **on the
full coordinate axis, masked samples included** — because "what would the GP have
said here?" is exactly the question asked about an excluded region.
`gp_localisation` evaluates it across posterior draws into the reserved group
`gp_localisation`.

All three are declared surfaces with no implementation. Phase 2 lands them behind
the backends that make them computable; W1.8's job is that their group names,
their inputs and their cost policy are settled first, so two backend tracks and
three diagnostic families write against one already-agreed answer.

## 8. The plotting surface

Six functions, each taking the emitted run and nothing else.

| Function | Family | Notes |
|---|---|---|
| `plot_corner` | — | selects by merged parameter name; an array-valued block is one variable, and a corner plot of a 10⁵-element latent block must be refused loudly rather than attempted |
| `plot_trace` | — | reads `sample_stats` too, so a prior-rejected draw shows as a gap (NaN), not as zero |
| `plot_posterior_predictive` | B | consumes `posterior_predictive`; its precondition is `add_posterior_predictive`, and its refusal must say so rather than reporting a missing group |
| `plot_residuals` | B | consumes `residuals`; warns when the run's likelihood provenance says a GP was fitted, because `diagnostics.md` §3.1 scopes family B to standard-likelihood fits |
| `plot_gp_localisation` | C | carries the degeneracy caveat by construction |
| `plot_anomaly_score` | A and C | one renderer, both provenances |

Two of these are contract rather than style.

**The GP-localisation caveat is mandatory and machine-readable.**
`diagnostics.md` §4.3 makes it "a plotting-function requirement for W1.8, not a
'please remember to mention this' note", and requires it to reach a user
extracting the score programmatically as well as one looking at the figure —
"a caption a user can silently crop out of a screenshot is not durable
protection against over-interpretation". So it is a module constant, it is in the
plotting function's own docstring, and it is reachable as a value:

```pycon
>>> gp_localisation_caveat().startswith("A large fitted GP amplitude localises")
True
>>> "does not say why" in gp_localisation_caveat()
True
>>> "it does not say why" in " ".join(plot_gp_localisation.__doc__.split())
True

```

**An anomaly score is never rendered without its provenance.**
`diagnostics.md` §5 resolves Tension 5 in favour of a shared visual grammar for
pre-fit RHMF flags and post-fit GP localisation, and polices the comparability
risk that creates with mandatory metadata rather than with visual distinctness —
"metadata travels with the data wherever it is consumed; visual distinctness only
protects the one plot function that respects it". `provenance` and
`interpretation_notes` are therefore required fields of the input, not optional
labels, and `show_provenance=False` may not suppress them when two scores of
different provenance are drawn together.

`diagnostics.md` §5 proposes the `AnomalyScore` container itself for
`ampere.core`, so that `ampere.diagnostics` (family A, which carries a JAX
dependency) and `ampere.results` can each produce one without either namespace
depending on the other. *(R4 granted and landed at the freeze, 2026-09-03:
`ampere.core.AnomalyScore` exists — `results_schema.md` §13 — and satisfies
the protocol.)* `plot_anomaly_score` stays typed against the *shape*
(`AnomalyScoreLike`, a runtime-checkable `Protocol`), so a caller may hand
it the real class or anything matching.

Every plotting function is a declared signature that raises `NotImplementedError`
naming Phase 2. That is deliberate: the surface is what two backend tracks and
three diagnostic families need agreed *before* they are written, and a
half-implemented plot would be a worse commitment than an honest refusal.
*(Amended W2.8: none of the six does any longer — W2.7 landed `plot_residuals`,
`plot_gp_localisation` and `plot_anomaly_score`, W2.8 `plot_corner`, `plot_trace`
and `plot_posterior_predictive`, and the surface above was not moved to make
either possible. What they raise now is `ResultsError`, on an input that is not
an emitted run.)*

## 9. Provenance: the recipe, stated once

Every run carries a flat set of `ampere_`-prefixed root attributes. Each value is
a netCDF-safe scalar — an `int`, a finite `float`, or a `str` holding canonical
JSON for anything structured.

**A boolean is not a netCDF type**, and it is worth saying so rather than
leaving it to be discovered: both engines refuse one (`netCDF4`: *illegal data
type for attribute*; `h5netcdf`: *boolean dtypes are not a supported NetCDF
feature*), and a bare `isinstance(value, int)` does not catch it because `bool`
*is* an `int` in Python. So `extra=` coerces booleans to `0`/`1` and numpy
scalars to their Python equivalents at the point of writing, not at the point of
serialising — an engine setting like `adapt=True` is exactly what `extra=` is
for, and it should not fail two steps later inside a backend.

```pycon
>>> attrs = provenance_attrs(joint, engine="emcee")
>>> attrs["ampere_seed"], attrs["ampere_seed_source"]
(20260902, 'explicit')
>>> json.loads(attrs["ampere_free_names"])
['model.index', 'model.norm', 'calibration']
>>> json.loads(attrs["ampere_sites"])["calibration"]
['blue.instrument.calibrate.scale', 'red.instrument.calibrate.scale']
>>> json.loads(attrs["ampere_dataset_labels"])
['blue', 'red']
>>> len(attrs["ampere_spec_hash"]), len(attrs["ampere_data_hash"])
(32, 32)

```

Three things about that list are decisions rather than defaults.

**The composed bindings, not `free_labels()`.** `inference.md` §4.6 says either
the composed bindings or `sites()` "is right for W1.8's provenance — they carry
the same information" since lossless nesting landed, and `sites()` is the
rendered form, so that is what is written. What is *not* written is
`free_labels()`: it is one string per flat-vector element, and a 10⁵-element
latent block's labels are exactly what `likelihoods.md` §16 tells this contract
not to materialise. `ampere_free_names` (one per parameter) and
`ampere_free_size` carry the same information in bounded space.

**The dataset ordering is provenance.** `ampere_dataset_labels` is written in the
collection's own order because that order defines which posterior element belongs
to which object (`hierarchical_population.md` §10.2).

**The joint spec hash and the per-component ones do not share a namespace.**
`spec_hashes` returns the joint entry under `"spec"` and the components under
`"components"`, because a component label is a user's choice and `spec` is an
ordinary word for a dataset. A flat mapping would let a dataset labelled `spec`
overwrite the joint entry silently, leaving `ampere_spec_hash` reporting one
component's declaration instead of the whole run's — which would quietly break
every property below and Phase 5's "may these two archived fits be reweighted
together?".

**The backend is derived, not declared** *(W2.12, decided by Fable 2026-09-07;
the decision-log row is in `DEVELOPMENT_PLAN.md` §2)*. `ampere_backend` used to
be whatever the engine driver was told to write down — `Engine(problem,
backend=...)` — which made it a claim rather than a record, and would have let
two lockstep backend tracks each invent their own spelling. The backend is now
`inference.md`'s **fourth capability flag**: models and transformations declare
`BACKEND` beside `DIFFERENTIABLE`/`BATCHABLE`/`DEVICE`, `declared_capabilities`
aggregates it by the device rule (all parts agree, or it raises), and
`provenance_attrs` reads it off `problem.capabilities`. `backend=` survives as
an optional **cross-check**: an explicit value that disagrees with the problem
raises rather than being recorded, because writing it down would be writing
down a falsehood and silently preferring either side would hide a real
configuration mistake. It is one name per backend everywhere — the same string
`lowering.md` §12.8's registry is keyed on, and the same string that ids that
backend's conformance fixture.

The same change put a `backend` key into `Capabilities.to_dict()`, so the
`ampere_capabilities` payload changed shape and **`PROVENANCE_SCHEMA_VERSION`
went to 4**, per the rule below. `capabilities` is not an input to
`problem_fingerprint` — but the schema constant is, so every
`ampere_problem_hash` moved at that bump exactly as at the previous ones.

**Three attributes joined at W2.13, and the constant is now 5** *(`inference.md`
§10a; the decision-log row "Realisation surface (W2.13)" in
`DEVELOPMENT_PLAN.md` §2, sub-decisions 10 and 11)*.

- **`ampere_realised`** — whether the draws were scored through the backend's
  **realisation**, its differentiable native form, rather than through the
  numpy contract path. Written on every run as `1`/`0` (netCDF has no boolean
  attribute type), including the gradient-free ones, because "this run's
  gradients were real" and "this run had no gradients" are the two answers a
  reader must be able to tell apart and silence distinguishes neither.
- **`ampere_registered_lowerings`** — the *user-registered* rows the run's
  lowering consulted, from `ampere.core.lowering.provenance_entries`. This is
  `lowering.md` §12.8's stamping, promoted from an `extra=` convention to a
  first-class key: W2.6 deferred that "until a real backend drives lowering end
  to end", and a realisation is exactly that. Empty for a run that used only
  ampere's own table, which is the signal — the question is "did this depend on
  something outside the conformance suite's guarantees?".
- **`ampere_solver_config`** — each dataset's `GPSolver.provenance_config()`,
  by dataset label. **Recorded and never hashed.** A solver's `jitter` is a
  *declaration* (it changes the number a given θ scores) and belongs in the
  spec and the spec hash; its dtype, its device and a future deliberate float32
  opt-out are *configuration* — two backends legitimately differ on them, and
  folding them into the spec hash would break §14's cross-backend agreement,
  which is the cheapest detector of a lowering bug this contract has. This
  attribute is what makes `architecture.md` §5's "a reduced-precision run must
  be visible in provenance" satisfiable without that cost.

None of the three is an input to `problem_fingerprint`, but the schema constant
is, so `ampere_problem_hash` moved again at this bump as at every previous one.

**Failures travel.** `ampere_failure_counts` is the unbounded count per
`FailureReason`; `ampere_failures` is the bounded history, each entry
`Failure.to_dict()`. `inference.md` limitation 17.7 notes that the history is
per-process and that a multiprocessing driver must aggregate before emitting —
that aggregation is the driver's, and this contract's part is that there is
somewhere for the answer to go.

### The hashing recipe

Four steps, and each is a decision.

1. **Normalise to a JSON tree.** numpy scalars become Python scalars; arrays
   become a `{dtype, shape, digest}` fingerprint; units become `to_string()`;
   enums become their values; and a non-finite float becomes one of three
   sentinel strings (`"__inf__"`, `"__-inf__"`, `"__nan__"`). Anything with no
   defined normal form is **refused**, because a record that silently dropped a
   field it did not recognise would compare equal to a run that differed in
   exactly that field.
2. **Serialise canonically**: `sort_keys=True`, no whitespace, UTF-8,
   `allow_nan=False`. Mapping order is therefore irrelevant and **list order is
   significant**:

   ```pycon
   >>> canonical_json({"a": 1, "b": 2}) == canonical_json({"b": 2, "a": 1})
   True
   >>> canonical_json([1, 2]) == canonical_json([2, 1])
   False

   ```

3. **Digest with BLAKE2b**, 16-byte output (32 hex characters), personalisation
   `b"ampere-prov"`. `hashlib`, never `hash()` — Python's built-in hash is salted
   per process, so a label-derived digest would differ between two runs of the
   same script, which is precisely the trap `ampere.core.rng` documents for
   seeds.
4. **Array bytes are normalised to little-endian, C-contiguous** before hashing,
   and the dtype is recorded by `dtype.name`, so a digest does not depend on the
   machine's byte order. Object arrays cannot be hashed by their bytes at all —
   those bytes are pointers — so they go through the normal form of `tolist()`.

### Why list order is load-bearing, and not fastidiousness

`lowering.md` §9.2: numpyro's `seed` handler splits its key once per stochastic
site, **in trace order**, rather than folding the site name in. The draws a
lowered model produces therefore depend on the order ampere emits its parameters
in, which means the seed alone does not identify a run — "adding a parameter to a
model changes the draws of every parameter emitted after it, and ... a run's
provenance must record the parameter-set spec (not just the seed) for a result to
be reproducible."

`ParameterSet.to_spec()` is an *ordered* list, so hashing the **merged** set
captures exactly that order. This is why `ampere_spec_hash` is the hash of
`problem.parameters.to_spec()` and not of a sorted union of the components'
specs; `ampere_component_spec_hashes` carries the per-component digests as well,
purely so that when two runs disagree, the record says where.

A changed declaration is a changed run, at the same seed:

```pycon
>>> class Steeper(Powerlaw):
...     def __init__(self, **grids):
...         self._channels = tuple(grids)
...         for name, grid in grids.items():
...             self.register_buffer(name, np.asarray(grid, dtype=float), unit=u.micron)
...         self.register_parameter(Parameter("index", st.norm(-2.0, 0.5)))
...         self.register_parameter(Parameter("norm", st.loguniform(0.1, 10.0)))
>>> steeper = FittingProblem(
...     Steeper(blue=blue_grid, red=red_grid),
...     [Dataset(blue_data, Instrument([Calibrate()], channel="blue"), label="blue")],
...     seed=20260902,
... )
>>> provenance_attrs(steeper)["ampere_spec_hash"] == attrs["ampere_spec_hash"]
False

```

and so is a changed dataset, with the spec hash untouched — which is what makes
the two hashes worth keeping separately:

```pycon
>>> moved = Spectrum(
...     blue_grid * u.micron, [1.0, 0.5, 0.26] * u.Jy, uncertainty=[0.05] * 3 * u.Jy
... )
>>> shifted = FittingProblem(
...     Powerlaw(blue=blue_grid, red=red_grid),
...     DatasetCollection({
...         "blue": Dataset(moved, Instrument([Calibrate()], channel="blue")),
...         "red": Dataset(red_data, Instrument([Calibrate()], channel="red")),
...     }),
...     ties=[Tie(
...         "calibration",
...         ("blue.instrument.calibrate.scale", "red.instrument.calibrate.scale"),
...     )],
...     seed=20260902,
... )
>>> shifted_attrs = provenance_attrs(shifted)
>>> shifted_attrs["ampere_data_hash"] == attrs["ampere_data_hash"]
False
>>> shifted_attrs["ampere_spec_hash"] == attrs["ampere_spec_hash"]
True

```

Free-form `meta` on a container is deliberately excluded from the data hash: a
changed comment is not a different dataset.

### `ampere_problem_hash` — the cache key

`DEVELOPMENT_PLAN.md` §7's trap list: "cache keys must hash the model/prior/data
spec so stale artefacts are invalidated automatically — the recent SBI caching
bugs on master are the evidence this bites." `problem_fingerprint` is what that
key is taken over, and the whole difficulty is that **the parameter spec is not
enough**. Three things change a run's numbers without changing a single
parameter declaration, and all three are covered here deliberately.

**The likelihood's structure.** Family, noise model, kernel, censoring, and the
solver **together with its configuration** — not merely its name. Two
likelihoods differing only in Matérn-3/2 versus squared-exponential, or
`DenseGP` versus `QuasisepGP`, have identical `ParameterSet` specs; and two
`DenseGP`s differing only in `jitter` produce different log-likelihoods at the
same θ, because the jitter is added to the diagonal before the factorisation.
Recording the strategy's name alone would have missed the second. This is what
`likelihoods.md` §17 question 8 asks about, and §15's R7 recommends moving it
onto `Likelihood` itself.

**Buffers.** A buffer is by definition the thing nobody puts a prior on
(`architecture.md` §6), so it appears nowhere in `to_spec()` — and it is the
wavelength grid, the opacity table, the filter curve, the response matrix. A
model rebuilt on a different grid is a different forward model with the same
declaration, and an emulator trained against one response matrix must not be
served for a fit against another. Buffers are hashed by content, for the model
and for every step of every instrument chain.

**Model identity.** Two models of different classes can declare the same
parameters and compute entirely different things, so the class and its module
are part of the fingerprint too.

What the fingerprint reaches is therefore **parameters, buffers, declared class
identity, ampere's own configuration objects, and whatever the model declares
through `describe()`**. Hashing an arbitrary `__dict__` is not a safe general
answer (it would sweep in caches, file handles and unhashable state), so the
last of those is opt-in rather than automatic.

**Amended and landed W2.1** (ruled 2026-09-03 at the freeze's escalations;
decision-log entry in `DEVELOPMENT_PLAN.md` §2). `model_fingerprint` gained a
`describe` key carrying `Parameterised.describe()`'s return value — `None`
unless the model overrides it. A `Redden(law="ccm89")` whose behaviour is set by
a string that is neither a parameter nor a buffer therefore no longer shares a
cache key with `Redden(law="f99")`, which was limitation 13 of §13. The payload
must be normalisable by the recipe above; anything else is refused when the hash
is taken, rather than silently omitted. Because adding a key changes every
`ampere_problem_hash`, this rode a `PROVENANCE_SCHEMA_VERSION` bump to 3.
*(And then to 4 at W2.12, when `Capabilities.to_dict()` gained its `backend`
key — see §9. Adding a key to the recorded `capabilities` payload does not
change the fingerprint's inputs, but the schema constant is one of them, so the
hash values moved anyway.)*

Two smaller consequences, recorded rather than left implicit: `describe` is now
a class attribute of `Parameterised`, so it joins `parameters`, `buffers` and
`context` as a name a parameter or buffer may not shadow; and a `describe()`
payload must be **backend-neutral and deterministic**, because §14's neutral
identity and W1.10's cross-backend equivalence row both compare it across
implementations. A device string or a dtype does not belong in one.

**The derived neutral identity.** `neutral_model_identity(model)` is this
fingerprint minus its `class` and `module` keys, and `model_identity_hash` is
its digest, recorded as `ampere_model_identity_hashes`. §14 states what it
licenses — and, more importantly, what it does not.

```pycon
>>> description = describe_likelihood(joint.datasets["blue"].likelihood)
>>> description["family"], description["noise"], description["marginalisation"]
('gaussian', 'IndependentNoise', 'analytic')
>>> len(hash_of(problem_fingerprint(joint)))
32

```

Moving a buffer leaves the spec hash alone and moves the problem hash, which is
exactly the division of labour the two hashes are for — the spec hash answers
"was this the same declaration?", the problem hash answers "may I reuse what I
computed last time?":

```pycon
>>> def on_grid(grid):
...     observed = Spectrum(
...         grid * u.micron, [1.0, 0.5, 0.25] * u.Jy, uncertainty=[0.05] * 3 * u.Jy
...     )
...     return FittingProblem(
...         Powerlaw(blue=grid),
...         [Dataset(observed, Instrument([Calibrate()], channel="blue"), label="blue")],
...         seed=20260902,
...     )
>>> coarse, fine = on_grid(blue_grid), on_grid(np.array([1.0, 2.0, 4.000001]))
>>> provenance_attrs(coarse)["ampere_spec_hash"] == provenance_attrs(fine)["ampere_spec_hash"]
True
>>> hash_of(problem_fingerprint(coarse)) == hash_of(problem_fingerprint(fine))
False

```

## 10. netCDF, and the dependency decision

Serialisation is netCDF, per `DEVELOPMENT_PLAN.md` §4.6, and the round trip is
the acceptance criterion of this item:

```pycon
>>> with tempfile.TemporaryDirectory() as folder:
...     written = to_netcdf(run, pathlib.Path(folder) / "run.nc")
...     back = from_netcdf(written)
...     restored = (
...         sorted(back.children),
...         int(back.attrs["ampere_seed"]),
...         back.attrs["ampere_spec_hash"] == run.attrs["ampere_spec_hash"],
...         bool(np.allclose(
...             back["log_likelihood"]["blue"].values,
...             run["log_likelihood"]["blue"].values,
...             equal_nan=True,
...         )),
...     )
>>> restored
(['constant_data', 'log_likelihood', 'observed_data', 'posterior', 'sample_stats'], 20260902, True, True)

```

**arviz stays an optional dependency, imported lazily.** `architecture.md` §4
rule 2 names `ampere.results` explicitly among the namespaces that "must import
optional dependencies lazily, inside the function/method/`__init__` that needs
them", and rule 3 requires the failure to be `OptionalDependencyError`, raised on
use and never on package import. That is implemented literally, and the existing
minimal-install CI job (`pip install -e .`, then the import sweep) enforces it
without anybody having to remember: every module under `ampere/` is walked and
must import cleanly with no extras.

`architecture.md` §3's extras table anticipates arviz being "folded into the base
install when `ampere.results` lands". This item lands the namespace but not the
capability — nothing yet *emits* a run, because no engine driver exists until
Phase 2 — so promoting arviz now would enlarge the base install (arviz pulls
xarray and pandas) in exchange for nothing a user can currently do. §15's R1 asks
Peter to time the promotion with Phase 2's engine drivers instead, and to decide
it together with the second half of the problem, which the extras table does not
currently mention:

> **arviz does not require a netCDF engine.** `pip install "ampere[arviz]"`
> installs something that can build a run and cannot write one; the failure
> surfaces as a bare xarray `ValueError` listing engine names. Whichever way R1
> goes, the `arviz` extra should name an engine — `h5netcdf` plus `h5py`, or
> `netCDF4`. In the meantime `ampere.results` catches that `ValueError` and
> re-raises it as `OptionalDependencyError` naming the remedy, and the pixi
> `dev` environment installs both engines so this document's round trip runs.

## 11. Serialising containers, results and training sets

`results_schema.md` §15.7 and §17.6 leave container serialisation unclaimed and
route the decision here, with design horizon (c)'s "serialisable training sets"
for emulators as the requirement to design against, "nailed down before Phase 2
rather than after". The answer is in two layers.

**Layer 1 — a plain-data form per container**, versioned, built only from lists,
strings, numbers, booleans and `None`, and round-tripping by value through the
**base** contract's validation — axis names and physical types, coordinate
ordering, shapes, complex support, non-negative uncertainties, a strictly
boolean mask — so a record tampered with in any of those ways fails at the
boundary rather than inside a likelihood. An invariant a *subclass* declares in
its own `__init__` is not re-checked, which today means exactly one thing and is
limitation 12 of §13:

```pycon
>>> encoded = container_to_dict(blue_data)
>>> encoded["kind"], encoded["unit"], encoded["coordinates"]["spectral_axis"]["unit"]
('Spectrum', 'Jy', 'micron')
>>> container_from_dict(encoded) == blue_data
True
>>> simulation = joint.simulate({"model.index": -1.0, "model.norm": 1.0, "calibration": 1.0})
>>> pair = training_pair_to_dict(simulation.parameters, simulation.results["model"])
>>> sorted(pair), sorted(pair["result"]["channels"])
(['failed', 'result', 'theta', 'version'], ['blue', 'red'])

```

The pair is exactly what `inference.md` §13 already produces —
`Simulation.results` "carries the raw `ModelResult` per model — the
`(θ, ModelResult)` pair design horizon (c) wants for emulator training sets" —
and `failed` travels with it so a rejected draw is recorded rather than lost,
which is the only way `simulate`'s "reject-and-record rather than train on
garbage" survives to the file.

A kind ampere does not know is refused by name, and an out-of-tree kind registers
itself once with `register_kind`, the same way it already registers with
`ampere.core`'s other extension points.

**Layer 2 — the on-disk training-set format: netCDF, the same as everything
else.** Specified here, and **implemented at W2.8** in `ampere.results.training`
(`write_training_set` / `append_training_set` / `read_training_set`); the two
rows marked below are that item's additions, which
`serialisation_review.md` §4 asked the writer for and which this table, written
before the review, does not otherwise have a home for:

| Group | Contents | Dims |
|---|---|---|
| root attrs | `ampere_spec_hash`, `ampere_problem_hash`, `ampere_data_hash`, versions, seed — the same recipe as §9 | — |
| `theta` | one variable per merged parameter name | `(sample,)`, plus the parameter's own dimensions |
| `<model>.<channel>` | one group per model channel: `values`, and `uncertainty`/`mask` where present | `(sample,)` + the channel's coordinate dimensions |
| `coordinates` | each channel's coordinate arrays, with units as attributes | the channel's dimensions |
| `sample_stats` | `failed`, and the whole `Failure` record — `failure_reason`, `failure_message`, `failure_where`, `failure_exception_type`, `failure_values` *(the last four amended W2.8)* | `(sample,)` |
| `observations` *(added W2.8)* | one subgroup per dataset label, holding the noisy draws `simulate(observe=True)` produced: `values`, and `uncertainty`/`mask` where present. Absent from a set whose budget drew none | `(sample,)` + the dataset's coordinate dimensions |

Three properties, each of which is why the format is netCDF and not JSON or a
pickle. NaN is native, so a masked or crashed sample needs no sentinel. The
coordinate arrays are stored once for the whole set rather than per pair, which
is what makes a coordinate-conditioned (neural-operator style) emulator's
training set the same size as a fixed-grid one. And the spec hash sits in the
attributes, so `DEVELOPMENT_PLAN.md` §7's "spec-hash invalidation of trained
artefacts" is a comparison of two strings rather than a convention somebody has
to remember.

The one thing this format deliberately does **not** do is store the model. An
emulator is trained on `(θ, ModelResult)` pairs and validated against the spec
hash; reconstructing the simulator that produced them is out of scope for a
training set and is what the composed problem is for.

## 12. Decisions and their reasoning

| Decision | Reasoning |
|---|---|
| The scalar joint log-likelihood in `sample_stats`, the per-dataset decomposition in `log_likelihood` | Both readings of §4.6 are legitimate and they differ; the per-dataset one is the only decomposition that is well defined for a GP likelihood (§5) |
| The `log_likelihood` group declares its own decomposition in an attribute | `arviz.loo` over it computes leave-one-dataset-out; a consumer must not have to infer that from the shape |
| Per-observation terms are named but not computed | They need a method on `Likelihood`, which is a merged contract; and two decompositions must not share a group name (§6) |
| `log_likelihood` is NaN, never `-inf`, for a prior-rejected draw | `inference.md` §18(c): "not evaluated" ≠ "impossible", and the coercion is irreversible |
| Posterior-predictive replicates and residuals are on demand | `N_draws × N_obs` per dataset per run, for something family B needs only when asked (§7) |
| An array-valued parameter is one variable with a named dimension | `likelihoods.md` §16(a): 10⁵ scalar names is the wrong representation |
| A plate's coordinate is read off `Binding.index`, not the dataset ordering | It is a fact the composition already recorded; reading it stays correct when the wiring is not in collection order (`hierarchical_population.md` §10.2) |
| The spec hash is over the **merged** set, in its own order | numpyro's seeding is trace-order dependent (`lowering.md` §9.2), so order is part of a run's identity |
| The problem hash covers the likelihood's structure, every buffer's contents, and the model class | Each changes a run's numbers without changing one parameter declaration; a parameters-only key never invalidates a stale artefact (§9) |
| `hashlib.blake2b`, never `hash()` | Salted per process; the same argument `ampere.core.rng` makes for seeds |
| Non-finite floats become sentinel strings, so `allow_nan=False` stays on | JSON has no `NaN`; the alternative is a non-portable token in a netCDF attribute |
| An unrecognised object is refused rather than dropped from a record | A record missing a field it did not understand compares equal to a run that lacked it |
| arviz is lazily imported and stays an extra for now | `architecture.md` §4 rule 2 names this namespace; and this item lands the namespace, not the capability (§10) |
| Every plot is a declared signature raising `NotImplementedError` | The surface is what two backend tracks and three diagnostic families need agreed first; a half-drawn plot is a worse commitment than a refusal. *(Amended W2.8: all six are drawn now, against the surface fixed here and unchanged by the drawing — see §8)* |
| The GP-localisation caveat is a constant, a docstring and a function | `diagnostics.md` §4.3 requires it to reach a programmatic consumer, not only a viewer |
| Container serialisation is functions in `ampere.results`, not methods on the containers | A hot-loop object should not carry the one method no evaluation calls; and `results_schema.py` is a merged contract (§15 R5) |
| Training sets are netCDF | NaN is native, coordinates are stored once, and the spec hash sits in the attributes where invalidation can see it |
| `ResultsError` lives in `ampere/core/exceptions.py`, re-exported here | §15 R5 asked for the move and it was made the same day the ruling landed (2026-09-03): one class, two import paths, no call-site changes — pinned by a test at the freeze |

## 13. Deliberate limitations of v1.8

Each is a decision, not an oversight. Each has an extension point.

1. **No per-observation log-likelihood.** §6 names the decomposition and the
   group; what is missing is a method on `Likelihood`. Extension point: R2.
   *(Closed: R2 was granted at the freeze and `Likelihood.pointwise_log_prob`
   landed with it; **W2.8** added `ampere.results.add_pointwise_log_likelihood`,
   which writes the group on an explicit call and never by default. What stays
   deliberately unreached is a solver that does not implement
   `conditional_loo` — `QuasisepGP` and its jax twin — where the emission
   refuses by name rather than falling back to a dense solve.)*
2. **No engine drivers.** Nothing in ampere currently produces the draws `emit`
   consumes; `DrawRecorder` is the shape a driver will use, exercised here by the
   test suite rather than by a sampler. Phase 2's drivers are the consumers.
3. **The plotting functions are signatures.** They raise; Phase 2 draws.
   *(Closed: W2.7 drew three and W2.8 the other three.)*
4. **`warmup_*` groups are not emitted.** ArviZ has a convention for them and
   ampere has no sampler to produce them yet. The group names are ArviZ's and
   cost nothing to adopt when a driver has warmup to store.
5. **`prior` and `prior_predictive` groups are not emitted.**
   `FittingProblem.sample_prior` and `simulate` make them cheap, but they are a
   driver's choice of what to spend a budget on, not this contract's.
6. **One `emit` call is one rectangular run.** Chains of different length are
   refused rather than padded, because ArviZ's `(chain, draw)` layout is
   rectangular and padding would invent draws. A driver trims first.
7. **Multiprocessing failure counts are not aggregated here.**
   `inference.md` limitation 17.7 puts the aggregation on the driver; this
   contract provides the attribute it goes into.
8. **`meta` is refused rather than dropped when it will not serialise.** A
   container whose metadata holds a live object cannot be written; the remedy is
   to keep such things outside the container. Silent stripping is the failure
   mode this refuses.
9. **The training-set writer is specified, not implemented.** §11's layer 2 is a
   table and a rationale; the writer lands with Phase 2's SBI and emulator work,
   against the format fixed here. *(Closed at **W2.8**:
   `ampere.results.training` writes, appends to and reads the format, and
   `training_pair_from_dict` completes the in-memory round trip.
   `serialisation_review.md` §4's two named losses are closed with it — θ keeps
   its dtype (`CONTAINER_SCHEMA_VERSION` 2) and `Simulation.observations` and
   the `Failure` detail are carried. One limitation is created and recorded in
   its place: **append is read-concatenate-rewrite**, `O(existing + new)` per
   call, because an in-place unlimited-dimension resize is a second
   serialisation path through h5netcdf rather than xarray. That is the
   extension point when a budget outgrows memory.)*
10. **`AnomalyScore` is a `Protocol`, not a class.** *(Closed at the freeze:
    R4 was granted and `ampere.core.AnomalyScore` landed 2026-09-03. The
    renderer stays typed against the shape, which the class satisfies.)*
11. **Data-group variable names can collide, and a collision is refused rather
    than resolved.** A name is the dataset label joined to an axis or role name,
    so a dataset `a` with uncertainties and a dataset `a_uncertainty` both want
    `a_uncertainty` — the same flattening collision `inference.md` §4.4 names for
    design A, inherited here because a flat namespace is what netCDF gives. The
    remedy is to rename a dataset, and the reason it is a refusal is that the
    alternative is silently dropping one dataset's data from the stored run. The
    extension point, if it ever bites in practice, is a nested group per dataset
    rather than a flat namespace, which ArviZ's own conventions do not use.
12. **Reconstruction re-checks the base contract, not a subclass's own
    `__init__`.** `container_from_dict` cannot call the kind-specific
    constructors — their signatures differ per kind by design — so it builds the
    base and inherits the base's checks. Today that misses exactly one
    invariant: `PhotometricPoints`' rule that filter names are unique. It
    matters only for a hand-edited or corrupted record, since anything ampere
    wrote had the check applied when it was first built. The extension point is
    a `validate()` classmethod on `FunctionSamples`, which is W1.4's to add and
    which this function would then call (R5).
13. **The problem hash cannot see a plain Python attribute.** Parameters,
    buffers, class identity and ampere's own configuration objects are all
    hashed; a user's `Redden(law="ccm89")`, configured by a bare attribute that
    is neither a parameter nor a buffer, is not — and two such fits share a
    cache key while scoring differently. Hashing an arbitrary `__dict__` is not
    a safe general answer, so the extension point is an opt-in `describe()`
    hook on `Parameterised` that a model or transformation implements when its
    behaviour depends on something the contracts do not model. **Ruled by
    Peter, 2026-09-03** (at the freeze's escalations): the hook is adopted for
    **early Phase 2** — it lands with W2.1, folded into `model_fingerprint`,
    with the decision-log entry ground rule 9 requires, so the cache-key hole
    is closed before any emulator cache exists to poison.

    **Landed W2.1**, so this is no longer a limitation of the shipped
    contract: `Parameterised.describe()` returns `None` by default and a
    normalisable mapping when a model overrides it, and §9 folds it into the
    fingerprint. What stays deliberately unreached is a plain attribute on a
    model that does *not* opt in — the author's choice now, rather than the
    contract's blind spot.

## 14. What this contract hands to the specs downstream

- **W1.10 (Conformance suite)** — candidate rows already implemented and testable
  here: the per-dataset contributions sum to the scalar joint `log_likelihood`;
  `lp == log_prior + log_likelihood` wherever finite; a prior-rejected draw
  stores `-inf`/NaN and no failure; a netCDF round trip preserves every group,
  every NaN and every provenance attribute; the spec hash changes when a prior
  changes and not when the data change, and the data hash vice versa; the digest
  is identical across processes with different `PYTHONHASHSEED`; a container of
  each kind round-trips by value; an array-valued parameter emits one variable
  with a named dimension. The **backend-spanning** rows this contract adds
  (claim corrected 2026-09-02, when W1.10 tested it): two backends emitting
  the same problem must produce the same `ampere_spec_hash` and
  `ampere_data_hash` — those two are functions of the declaration and the
  data, not of the arithmetic, so any disagreement is a lowering bug and the
  hash is the cheapest possible detector of it. `ampere_problem_hash` is
  **deliberately not backend-invariant**: it fingerprints the model's class
  and module (§9) precisely so that two differently implemented forward
  models never share a cache key, and a backend's lowered model is a
  different implementation. W1.10's row asserts that actual equivalence.
  Whether Phase 2 additionally wants a backend-neutral model identity — so
  an emulator trained on the reference backend can be *offered* (never
  silently served) to a torch fit of the same declaration — is a freeze
  question, and it belongs with §13.13's `describe()` hook. **Ruled by
  Peter, 2026-09-03**: yes, in early Phase 2, as one mechanism with that
  hook (W2.1): a *derived* neutral identity — the model fingerprint minus
  its class/module component — used only to **offer** a cross-backend
  emulator with its provenance shown, never to serve one silently; the
  neutral identity cannot pin the mathematics, which is why "offer, never
  serve" is the rule. `ampere_problem_hash` itself stays deliberately
  backend-variant.

  **Landed W2.1.** `provenance.neutral_model_identity(model)` is
  `model_fingerprint` minus `class` and `module` — so the parameter
  declaration, the buffers and the `describe()` configuration remain — and
  `model_identity_hash` is its digest, recorded per model as the
  `ampere_model_identity_hashes` attribute. Two conformance rows now hold
  the pair together: the pre-existing equivalence row (unchanged) asserts
  that `ampere_problem_hash` agrees exactly when the recorded model
  identity does, and a new row asserts that the *neutral* identity agrees
  across backends unconditionally. Both are real assertions on the in-repo
  pair, which compute the same closed forms by different routes in
  different classes and modules.

  Scope, stated because it is narrower than "a backend-neutral problem
  identity": this is the **model's** identity only. `_describe_instrument`
  records each chain step's class too, so a neutral identity for a whole
  `FittingProblem` would need that second site as well. This section words
  the ruling as "the model fingerprint minus its class/module component",
  and that is exactly what landed; the instrument case is unclaimed.
- **W1.12 (Diagnostics)** — the three answers it asked for: the reserved group
  names and the derivation for signed residuals (§7), the ruling that `y_rep` is
  computed on demand rather than stored (§7), and `plot_anomaly_score` as a
  first-class stub with the provenance/caveat requirement built in (§8).
- **W1.13 (Spec assembly & freeze)** — seven requests, all in §15, plus the
  wording correction in §2 (`InferenceData` is now `xarray.DataTree`) and the
  `architecture.md` §3 extras-table amendment in §10.
- **Phase 2 (backends)** — a backend supplies draws and `Evaluation`s and nothing
  else: the whole of this contract is reused unchanged, which is the same claim
  `inference.md` §18 makes for its own surface. The one obligation on a backend
  is that its flat vector is in `FittingProblem.parameters` order, which `emit`
  checks by width and cannot check by meaning.
- **Phase 3 (SBI)** — §11's training-set format is what a simulation budget
  writes, and `ampere_problem_hash` is the cache key that invalidates a stale
  posterior. The batched `simulate_many` that `inference.md` limitation 17.5
  defers is what will fill it efficiently.
- **Phase 5 (population inference)** — design horizon (b) is buildable entirely
  on stored runs: `sample_stats` carries the scalar `log_prior` and
  `log_likelihood` per draw, the NaN convention distinguishes an unevaluated
  draw from an impossible one, and `ampere_spec_hash` says whether two archived
  fits used the same declaration and may be reweighted together.

## 15. Open questions for review

**Ruled by Peter, 2026-09-03 — all seven, in the direction each request
recommends.** R1: arviz joins the base install with **Phase 2's engine
drivers**, promoted together with a netCDF engine, when a user can first
emit a run. R2: `Likelihood.pointwise_log_prob` is **granted at the
freeze but not stored by default** — W1.13 lands the §4.4 addition with
its decision-log entry. R3: the `InferenceData` → `xarray.DataTree`
wording is ratified — `DEVELOPMENT_PLAN.md` §4.6 and `architecture.md`
§3 corrected the same day, with the decision-log entry. R4:
`AnomalyScore` lands in `ampere.core` at W1.13. R5: `ResultsError`
**moved to `ampere/core/exceptions.py`** the same day (`ampere.results`
still re-exports it, so no import site changes); container serialisation
**stays functions** in `ampere.results`. R6: the per-dataset
`log_likelihood` group — the conventional name, keyed by dataset, with
the decomposition declared in the group's attributes — is confirmed. R7
rides with the consolidated cross-contract serialisation review at W1.13
(Peter's same-day `likelihoods.md` §17 Q8 ruling), where
`describe_likelihood`'s promotion to `Likelihood.to_spec()` is the
leading candidate rather than a separately ruled point. *(The review is
done — `docs/design/serialisation_review.md` — and confirmed the
promotion: `Likelihood.to_spec()` landed, `describe_likelihood` composes
it and keeps the content fingerprints, and the review found and fixed
one provenance hole: family- and noise-model-owned buffers were
invisible to `buffer_fingerprint(likelihood)`;
`PROVENANCE_SCHEMA_VERSION` is now 2.)* The original
requests are kept below for the record.

**R1 — when does arviz join the base install, and with which netCDF engine?**
`architecture.md` §3's extras table says arviz is "folded into the base install
when `ampere.results` lands". This item lands the namespace but not the
capability, so it has kept arviz an extra and imported it lazily, per §4 rule 2
(§10). The recommendation is to promote it with **Phase 2's engine drivers**,
when a user can actually emit a run, and to promote a netCDF engine at the same
time — `pip install "ampere[arviz]"` today installs something that can build a
run and cannot write one, which is a real gap whichever way the timing goes.
Confirm, or ask for the promotion now.

**R2 — should `Likelihood` gain `pointwise_log_prob`?** §6 names the
decomposition ArviZ's `log_likelihood` convention wants — pointwise terms for
independent noise, leave-one-out conditionals for a GP, both computable from the
Cholesky `DenseGP` already forms — but computing them needs a method on
`Likelihood`, which is a merged §4.4 contract. It is small and additive, it is
what makes `arviz.loo`/`waic` mean the usual thing, and `diagnostics.md` §3.3
reads `DEVELOPMENT_PLAN.md` §4.6 as already promising it. Against: it is
`N_draws × N_obs` of storage, nothing in Phase 1 consumes it, and the per-dataset
decomposition already answers "which dataset is driving this fit?". This document
recommends **granting it at the freeze but not storing it by default**, so the
group is available to a user who asks for LOO and absent from an archived run
that does not need it.

**R3 — ratify the `InferenceData` → `DataTree` wording.** ArviZ 1.0 retired the
`InferenceData` class in favour of `xarray.DataTree`. The *format* and the
decision are unchanged (§2), but `DEVELOPMENT_PLAN.md` §4.6 and
`architecture.md` §3 both name the class. A one-line correction in each, with a
decision-log entry, keeps the documents from disagreeing with the code. This item
has not made it, because §4.6 is the plan's own text.

**R4 — land `AnomalyScore` in `ampere.core`?** `diagnostics.md` §5 proposes it
and §7 flags it for W1.4, which had already landed. Until it exists,
`plot_anomaly_score` is typed against a `Protocol` here. The class is small
(coordinates, values, mask, `provenance`, `interpretation_notes`), plain numpy,
and belongs in core precisely so neither `ampere.diagnostics` nor
`ampere.results` depends on the other. Recommend landing it at W1.13 with the
freeze.

**R5 — two homes to ratify or move.** `ResultsError` is declared in
`ampere/results/exceptions.py` while every sibling contract error lives in
`ampere/core/exceptions.py`; it is there only because that file is merged and
widening it is a decision-log matter. Likewise, the container serialisation of
§11 is functions in `ampere.results` rather than `to_spec`/`from_spec` methods on
`FunctionSamples`/`ModelResult`, which is where `results_schema.md` §15.7 half
expects them. This document recommends **moving `ResultsError` to core** (it is
two lines and no call site changes) and **leaving the serialisation as
functions** — a hot-loop container should not carry the one method no evaluation
calls, and keeping it outside means a user-defined kind gets serialisation for
free rather than having to implement it.

**R6 — is the per-dataset `log_likelihood` group the right ArviZ citizenship?**
ArviZ's convention is that `log_likelihood` variable names match `observed_data`
variable names and share the observation dimension. Ampere's match by name but
have no observation dimension, because the terms are per dataset. The
alternatives were a differently named group (honest, but invisible to
`arviz.loo`) or no group at all until R2 lands (which would lose the
decomposition that *is* well defined). The chosen middle — the conventional group
name, keyed by dataset, with the decomposition declared in the group's attributes
— is the one this document recommends, but it is worth confirming, because
changing it after runs are archived is not cheap.

**R7 — should `Likelihood` gain `to_spec()`?** `likelihoods.md` §17 question 8
routes this here: "`ParameterSet.to_spec()` covers its parameters, and
`KernelSpec.to_dict()` its kernel, but there is no `Likelihood.to_spec()` for the
family name plus solver plus censoring. W1.8 owns provenance hashing and should
decide whether it wants one."

**It does, and the argument is a correctness one rather than a tidiness one.**
Two likelihoods differing only in Matérn-3/2 versus squared-exponential, or in
`DenseGP` versus `QuasisepGP`, have *identical* `ParameterSet` specs, so a
provenance record built from the parameters alone cannot tell them apart and a
cache key built on it never invalidates (§9). Something has to serialise the
triple, and the only question is where it lives.

`ampere.results.describe_likelihood` is that something today, assembled from the
public surface — `family.NAME`, `type(noise).__name__`, `solver.NAME`,
`kernel.spec().to_dict()`, `marginalisation`, the censoring counts, the
parameters and the buffers. That works, and it is deliberately written so that
promoting it is a move rather than a rewrite. The recommendation is to promote
it: `Likelihood.to_spec()`, returning the same mapping, at W1.13. Three reasons.
A likelihood knows things about itself that an outside reader has to
reverse-engineer (whether a solver's approximation parameters matter, what a
future family's extra state is), and each new family would otherwise need this
module amended in step. The conformance suite (W1.10) wants one definition to
compare two backends' likelihoods against, not two. And every sibling — the
kernel, the parameter set — already serialises itself, so this is the odd one
out rather than a new idea.

Against it: it is an addition to a merged §4.4 contract, needing a decision-log
entry; and if the answer is no, nothing breaks — this module keeps doing it, and
this document records that it is the definition.
