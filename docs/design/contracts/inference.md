# Ampere v2 — Dataset, FittingProblem & Inference Contract (W1.7)

Status: **DRAFT for Peter's review**, and it carries four ruling requests
(§19, R1–R4) that later work is blocked on. Two of Peter's earlier rulings —
`likelihoods.md` §17 Q1 (a strict toggle; the engine path records and returns
−inf) and Q2 (the `Dataset` resolves the effective mask once) — are implemented
here, in §11 and §8 respectively. Implements
`DEVELOPMENT_PLAN.md` §4.5, discharges the obligations `parameters.md` §13/§14,
`transformations.md` §14, `likelihoods.md` §16 and `results_schema.md` §17 place
on this item, and adopts `lowering.md` §9.2's seed-derivation policy. Code:
`ampere/core/dataset.py`, `ampere/core/rng.py`, `ampere/core/exceptions.py`.
Tests: `tests/core/test_dataset.py`.

Every worked example below is executed as a doctest by
`tests/core/test_spec_doctests.py`, so this document cannot drift from the
implementation without the suite going red. Examples share one namespace and
build on each other; read them in order.

Where this document and `DEVELOPMENT_PLAN.md` disagree, the plan wins and this
document is wrong — file it as a decision-log correction.

---

## 1. What this contract is for

Everything below this contract is vocabulary. Parameters say what can vary,
containers say what a model produces, instruments say how it is observed,
likelihoods say how well it matches. This contract is the sentence they compose
into, and it is the *only* thing an inference engine sees.

That is the whole architectural bet of `DEVELOPMENT_PLAN.md` §3: gradient-free
samplers and SBI live **above** the backends and consume any of them through
this surface, so they are written once. An engine that uses only what §4.5
lists — `log_prob`, the `log_likelihood`/`log_prior` split, `prior_transform`,
`simulate`, the capability flags — works with the reference backend, with torch,
with jax, and with a legacy black-box model behind a thin adapter, and never
knows which.

Legacy ampere had no such object. `Data` subclasses carried their own
`lnlike`, the sampler drivers reached into models by positional `theta` slicing,
and "the joint likelihood" was a loop written afresh in each `mixins.py` path.
Positional slicing dies with `parameters.md`; the loop dies here.

It is **backend-neutral**: numpy, `astropy.units` and stdlib, with no torch or
jax import, lazily or otherwise (`architecture.md` §4 rule 1).

### Setup for the examples

```pycon
>>> import math
>>> import numpy as np
>>> import scipy.stats as st
>>> import astropy.units as u
>>> from ampere.core import (
...     Capabilities, Dataset, DatasetCollection, Evaluation, FailureReason,
...     FittingProblem, GaussianFamily, GaussianProcessNoise, HierarchicalPrior,
...     IndependentNoise, Instrument, Likelihood, Marginalisation, Matern32, Model,
...     ModelResult, Parameter, ParameterSet, PhotometricPoints, PoissonFamily,
...     Spectrum, Tie, Transformation, declared_capabilities, substream,
... )
>>> from ampere.core.exceptions import DatasetError, LikelihoodError

```

Three things every later section reuses: a coordinate grid, one observed
spectrum, and a calibration step — the archetypal instrument nuisance, a single
multiplicative scale factor, which is what legacy ampere's `scaleFac` was.

```pycon
>>> grid = np.array([1.0, 2.0, 3.0])
>>> observed = Spectrum(
...     grid * u.um, [2.0, 4.0, 6.0] * u.Jy, uncertainty=[0.1, 0.1, 0.1] * u.Jy
... )
>>> class Calibrate(Transformation):
...     ACCEPTS = (Spectrum,)
...     def __init__(self, **kwargs):
...         super().__init__(**kwargs)
...         self.register_parameter(Parameter("scale", st.lognorm(0.2)))
...     def apply(self, samples, values):
...         return samples.with_values(samples.values * self.context(values)["scale"])
>>> class TwoChannel(Model):
...     """A stub model, and the one most examples below are fitted against."""
...     def __init__(self, wavelength):
...         self.register_buffer("wavelength", wavelength, unit=u.um)
...         self.register_parameter(Parameter("index", st.norm(-1.0, 0.5)))
...         self.register_parameter(Parameter("norm", st.loguniform(0.1, 10.0)))
...     def evaluate(self, **values):
...         ctx = self.context(values)
...         flux = ctx["norm"] * ctx["wavelength"] ** ctx["index"] * u.Jy
...         return ModelResult({"blue": Spectrum(ctx["wavelength"] * u.um, flux)})

```

## 2. The objects

| Object | Role |
|---|---|
| `Dataset` | One observed container, the `Instrument` that predicts it, the `Likelihood` that scores it |
| `DatasetCollection` | Several datasets fitted jointly, plus the `shared` hyperprior extension point |
| `FittingProblem` | Model(s) + collection + ties. Owns the one merge, the lifecycle, and §4.5's surface |
| `Capabilities`, `Capable` | `differentiable` / `batchable` / `device`, and what a backend declares them on |
| `Evaluation` | One `log_prob` call's full result: the split, the per-dataset terms, the failure |
| `Failure`, `FailureReason` | §4.5's "recorded reason", as a counted vocabulary rather than free text |
| `Simulation` | `simulate`'s result: (θ, ModelResult, prediction, optional observation), or a flagged failure |
| `substream`, `generator` | `lowering.md` §9.2's seed derivation, on the reference path |
| `DatasetError` | This contract's error; a `ContractError`, hence a `ValueError` |

A `Dataset` deliberately holds **no model**. The commonest joint fit in
ampere's target scope is one physical model observed by several instruments, so
a model inside a dataset would either be duplicated or shared by object
identity — and shared-by-identity is exactly `prior_art.md` Tension 1's gammapy
footgun, where mutating one dataset's model silently mutates another's. A
dataset names its model with a *string*, the same way its instrument names its
channel.

## 3. The simple path

One model, one spectrum, no instrument, no nuisance parameters. This must stay
two lines, and it does:

```pycon
>>> class Line(Model):
...     def __init__(self, wavelength):
...         self.register_buffer("wavelength", wavelength, unit=u.um)
...         self.register_parameter(Parameter("slope", st.norm(1.0, 1.0)))
...     def evaluate(self, **values):
...         ctx = self.context(values)
...         return Spectrum(ctx["wavelength"] * u.um, ctx["slope"] * ctx["wavelength"] * u.Jy)
>>> problem = FittingProblem(Line(grid), [Dataset(observed)])
>>> problem.parameters.free_names
('model.slope',)
>>> round(problem.log_prob({"model.slope": 2.0}), 6)
2.732001

```

`Dataset(observed)` alone means "no instrument, i.i.d. Gaussian noise, the
default channel" — the model already produces the observable, and the
container's own `uncertainty` is what the likelihood uses.

```pycon
>>> dataset = Dataset(observed)
>>> dataset.channel, dataset.parameters.names
('default', ())
>>> type(dataset.likelihood.family).__name__, type(dataset.likelihood.noise).__name__
('GaussianFamily', 'IndependentNoise')

```

The default instrument is a pure channel binding **kind-checked against the
observed container**, so binding a channel of the wrong kind fails loudly rather
than producing a shape error inside a likelihood:

```pycon
>>> wrong = ModelResult(PhotometricPoints(["W1"], [3.4] * u.um, [1.0] * u.Jy))
>>> dataset.predict(wrong)
Traceback (most recent call last):
    ...
ampere.core.exceptions.ChannelError: channel 'default' holds a PhotometricPoints, but a
Spectrum was required. ...

```

## 4. The joint parameter space: the merge topology

This is the contract's central design decision, and `parameters.md` §14 keeps
it expressly open for W1.7 to settle rather than inherit. It is written out at
length because Peter rules on it, not because it is complicated.

### 4.1 The constraint

`ParameterSet.merge` qualifies each component's names with the component's
label, and a label must be a bare Python identifier
(`parameters.md` §8). So **one merge call produces exactly one level of
qualification**. `merge` is also not associative: merging a
`ParameterMapping.merged` again and *discarding* the first mapping silently
drops its bindings, which is `parameters.md` §12.4's warning.

Crucially, a level of nesting already exists and is not W1.7's to remove:
`Instrument` merges its steps (`transformations.md` §5), so `Instrument.
parameters` is *already* a merged set with dotted names like
`calibrate.scale`. The question was never "nest or not"; it is "how many
levels, and who re-distributes".

### 4.2 Design A — one flat merge over every leaf

The problem enumerates every parameter-owning leaf in the fit — the model,
every transformation step of every instrument, every likelihood, every latent
block — and merges them all in exactly one call. Because a component label must
be a bare identifier, the labels must be synthesised by flattening:

```
components = {
    "model":            model.parameters,
    "sed_calibrate":    <that step>.parameters,
    "sed_likelihood":   <that likelihood>.parameters,
    "spectrum_lsf":     <that step>.parameters,
    ...
}
```

giving merged names `sed_calibrate.scale`, `spectrum_lsf.fwhm`.

### 4.3 Design B — nested `ParameterMapping`, one merge per level

Each composite performs one merge over its **immediate** children and *retains*
the mapping; values flow down by calling `distribute` at each level. Three
levels, and no more:

| level | components | example merged name |
|---|---|---|
| `Instrument` (W1.5) | its steps' labels | `calibrate.scale` |
| `Dataset` | `instrument`, `likelihood`, `latent` | `instrument.calibrate.scale` |
| `FittingProblem` | model labels, dataset labels, `shared` | `sed.instrument.calibrate.scale` |

The dataset level's labels are *role* names, not object labels. That is
deliberate: `Instrument.label` defaults to the channel name, which is also the
commonest dataset label, so using it here would produce `sed.sed.calibrate.
scale`.

### 4.4 The comparison

| Question | A (flat leaf merge) | B (one merge per level) |
|---|---|---|
| Merged names | `sed_calibrate.scale` — one level of qualification whatever the depth | `sed.instrument.calibrate.scale` — one segment per level |
| Component labels | synthesised by flattening; collisions constructible (dataset `a` + step `b_c` collides with dataset `a_b` + step `c`) and so needing an explicit check, as B's top level also has | one namespace per level; the dataset level uses two reserved role names, the top level is checked |
| `Tie` reach | any leaf | **any site, by its full path — including one inside an inner merge** (§6) |
| `shared_as` reach | the whole fit | within one merge level (§6, limitation 17.1) |
| `Instrument.mapping` / `Instrument.__call__`'s values path | usable only after W1.7 re-keys values back into instrument-merged names, duplicating `merge`'s qualification rule outside it | used exactly as `transformations.md` §5 wrote them |
| Hierarchical prior references | rewritten once | rewritten once per level, and they compose — including through a tie (§9) |
| Merges per evaluation | 1 | 1 (the problem's, cached) + 1 per instrument per call (W1.5's own, §11) |
| Available mistake | none: there is one mapping | a level that discards its mapping loses its bindings |

### 4.5 The ruling: B, and why

**Design B is implemented.** Three reasons, in decreasing order of weight.

**1. A duplicates `merge`'s qualification rule outside `merge`.** Under A the
problem holds `sed_calibrate.scale`, while `Instrument.__call__` wants
`calibrate.scale`. So the problem must either **re-key by convention** —
rebuilding `f"{step.label}.{local}"` inside W1.7, for names `ParameterSet.merge`
already knows how to build — or bypass `Instrument.__call__` and drive the steps
itself, re-implementing the chain's kind and mask checks.

An earlier draft of this section claimed A could use none of W1.5's entry points.
That was wrong, and checking it is what turned this from the decisive reason into
the weakest of the three: the re-keying route is four lines and **works**. What
it costs is that the qualification convention then lives in two places, and
`Instrument.parameters` / `Instrument.mapping` play no part in building the joint
space at all. Under B the same names arrive without reconstruction, because they
are exactly what the inner merge produced.

**2. A's headline advantage is smaller than it looks.** The apparent advantage
is "one place where tying is resolved, so a `Tie` may name any two leaves". But
an outer tie under B *already* reaches any inner site: the outer merge sees the
inner merged name (`instrument.calibrate.scale`) as an ordinary local name of
the dataset component and collapses it exactly like any other. §6 demonstrates
this and the suite tests it. The real difference is narrower and is stated as
limitation 17.1: **declaration-time `shared_as` does not cross a merge level.**

**3. Names are not cosmetic.** A merged name is a `Tie` site, a `free_labels()`
entry, an ArviZ coordinate, a corner-plot axis label, and a string a user types.
`sed.instrument.calibrate.scale` says where to look. `sed_calibrate.scale`
requires knowing the flattening convention to parse, cannot be parsed
unambiguously at all, and is one level deep however deep the structure is.

**On the weight of the case.** With reason 1 corrected, B rests mainly on
reasons 2 and 3, and the honest summary is that this is closer than a first
reading suggests: A's real cost is a duplicated naming convention and worse
names, not a broken one. What tips it is that B's own defect (§4.6) has a fix
at source — "lossless nesting" in `ParameterSet.merge` — whereas A's duplicated
convention has none short of not doing it. Peter should read §4.6 before ruling.

And the hazard `parameters.md` §12.4 names is **contained by a rule, not by
luck**:

> **The nesting rule.** Every composite performs exactly one `merge` over its
> immediate children and retains the resulting `ParameterMapping`. Values flow
> down by calling `distribute` at each level. No mapping is ever discarded, so
> no binding is ever lost.

Re-distributing cannot lose a binding; it re-derives one, one level down. What
§12.4 warns against is merging a `merged` set and *throwing the mapping away*,
which no level here does.

### 4.6 The strongest objection: nesting loses introspection

This is the one real cost, and it was found by W1.11's hierarchical-population
sketch rather than by this contract, so it is recorded here in its own section
rather than buried in a limitation.

**The loss.** A `ParameterMapping` records one binding per component, not per
leaf. So a parameter collapsed by an *inner* merge — two steps of one instrument
sharing a declaration-time `shared_as` label — reaches the top level as a single
local name with a single binding. The top-level mapping therefore says it is not
tied, and any consumer walking `mapping.bindings` for provenance or for ArviZ
labelling will not learn that the sampler dimension drives two places:

```pycon
>>> class Gain(Transformation):
...     ACCEPTS = (Spectrum,)
...     def __init__(self, **kwargs):
...         super().__init__(**kwargs)
...         self.register_parameter(Parameter("scale", st.lognorm(0.2), shared_as="gain"))
...     def apply(self, samples, values):
...         return samples.with_values(samples.values * self.context(values)["scale"])
>>> hidden = FittingProblem(
...     TwoChannel(grid),
...     [Dataset(
...         observed,
...         Instrument([Gain(label="a"), Gain(label="b")], channel="blue"),
...         label="d",
...     )],
... )
>>> hidden.parameters.free_names
('model.index', 'model.norm', 'd.instrument.gain')
>>> hidden.tied_names
()

```

**The recovery.** The information was never destroyed — the nesting rule retains
every level's mapping — it was only not surfaced by the outermost one. So this
contract surfaces it, by descending:

```pycon
>>> hidden.shared_names
('d.instrument.gain',)
>>> hidden.sites()["d.instrument.gain"]
('d.instrument.a.scale', 'd.instrument.b.scale')

```

`FittingProblem.sites()` maps every merged name to the fully qualified leaf
paths it feeds, at whatever level they were collapsed, and `shared_names` is the
lossless counterpart of `tied_names`. **W1.8 should record `sites()`, not the
raw bindings**, in a run's provenance: it is the honest answer to "which parts of
the model did this sampler dimension drive?".

**The honest residual.** The *default* introspection surface —
`problem.mapping.tied_names` and `problem.mapping.bindings` — is now misleading
unless a consumer knows to reach for `sites()`. A helper that must be remembered
is worse than a surface that is right, and this is the one place where design A
would genuinely have been simpler.

The right fix is upstream and is smaller than either candidate here.
W1.11's sketch calls it **lossless nesting**: let `ParameterSet.merge` accept a
`ParameterMapping` as a component and *compose* its bindings, so the outer
mapping's bindings are already the leaf bindings and `sites()` becomes
unnecessary. That is strictly smaller than recursive merge, it fixes the problem
where the problem is, and it belongs in `parameters.md` rather than here.
Recommended to W1.13 as part of ruling request **R1**.

### 4.7 What the population sketch says, and what it does not

W1.11's `docs/design/modalities/hierarchical_population.md` §§5–7 evaluated the
same question from the population side. Three of its findings bear on this
ruling, and it is worth being precise about which way each cuts.

**Corroborating.** Merging an already-merged set does not fail, and routing
survives two levels: the caller re-runs the inner mapping's `distribute` on the
outer result and values arrive correctly. That is the mechanism §4.3 depends on,
verified independently. Cross-level ties also work — a composition-time tie
across two objects' inner instrument parameters collapses to three free
parameters rather than four — which is the same result §6 shows here, and which
answers `transformations.md` §14's ties-across-levels worry from both ends.

**Neutral, and it would be dishonest to claim otherwise.** The sketch's verdict
is that recursive merge is *not* what population models need. The pressure there
is **per-element routing of a `Plate`'s array-valued parameter**: `Binding` has
no index field, `distribute` hands the whole `(N,)` block to one component, and
tying across plate members is refused by `TyingError`. None of that is affected
by the merge topology — design A would pay it identically — so the population
evidence does not discriminate between A and B, and the critical-path item it
names, an optional `Binding.index`, is `parameters.md`'s to add.

**A warning this contract must carry.** The sketch measured the alternative to a
plate: N scalar components cost about 800 ms per `lnprior` at N = 1000, against
about 0.6 ms for the same structure expressed as one `Plate`. The cost is the
number of `Parameter` objects, not the merge. That is a three-order-of-magnitude
difference inside the hot loop, and it bounds §9's "plate of datasets" pattern
hard: it is right for tens of spaxels and wrong for thousands, and §9 now says
so.

### 4.8 What a reversal costs

If Peter prefers A (ruling request **R1**):

- `ampere/core/dataset.py` changes in three places — `Dataset._components` and
  `Dataset.route` disappear (a `Dataset` then publishes a flat
  `{f"{label}_{step.label}": step.parameters, …}` dict rather than a merged
  set), and `FittingProblem._components` flattens. About 60 lines.
- W1.7 gains a re-implementation of `Instrument.__call__`'s values path (about
  eight lines) plus the duplication risk above.
- **Every merged name in the joint space changes.** Every `Tie` site string,
  every test asserting a name, and every example in this document change with
  it. That is the expensive half.

Estimate: half an agent session now. Considerably more once W1.8 is emitting
these names into stored `InferenceData` and W1.10's conformance suite asserts
them. **Reverse now or not at all** is the honest summary.

## 5. `Dataset`

```pycon
>>> flexible = Likelihood(
...     GaussianFamily(),
...     GaussianProcessNoise(Matern32(st.loguniform(1e-3, 1e1), st.loguniform(0.1, 10.0))),
... )
>>> dataset = Dataset(observed, Instrument([Calibrate()]), flexible, label="spectrum")
>>> dataset.parameters.free_names
('instrument.calibrate.scale', 'likelihood.amplitude', 'likelihood.length_scale')

```

The two role labels are reserved words at this level, so a transformation
perversely labelled `likelihood` still cannot collide with the likelihood — its
parameters sit a level further down:

```pycon
>>> odd = Dataset(
...     observed,
...     Instrument([Calibrate(label="likelihood")]),
...     Likelihood(GaussianFamily(), IndependentNoise(scale=st.loguniform(0.5, 2.0))),
... )
>>> odd.parameters.free_names
('instrument.likelihood.scale', 'likelihood.scale')

```

A dataset checks what it can at construction. The chain's output kind must be
the observed container's kind, because a likelihood compares like with like:

```pycon
>>> photometry = PhotometricPoints(
...     ["W1", "W2"], [3.4, 4.6] * u.um, [1.0, 2.0] * u.Jy, uncertainty=[0.1, 0.1] * u.Jy
... )
>>> Dataset(photometry, Instrument([Calibrate()]))
Traceback (most recent call last):
    ...
ampere.core.exceptions.DatasetError: dataset 'default': the instrument chain produces Spectrum
but the observed data are PhotometricPoints. ...

```

### `Dataset.label`, and a collision that is the *default*

`Dataset.label` defaults to the instrument's label, which itself defaults to the
channel name (`transformations.md` §15.4). So two photometric catalogues
observing one `sed` channel are **both** labelled `sed` unless one is named.
That would file two component sets under one dictionary key, silently
discarding the first *before* `merge` ever sees them — so `merge`'s own
collision detection could never fire. `DatasetCollection` therefore checks
before it merges:

```pycon
>>> def catalogue():
...     return Dataset(observed, Instrument([Calibrate()], channel="sed"))
>>> DatasetCollection([catalogue(), catalogue()])
Traceback (most recent call last):
    ...
ampere.core.exceptions.DatasetError: two datasets are labelled 'sed'. ...

```

Naming them is all it takes, and no parameter is lost:

```pycon
>>> pair = DatasetCollection({"gaia": catalogue(), "wise": catalogue()})
>>> sorted(pair.components())
['gaia', 'wise']

```

This is one place where design B pays for itself twice over: because instrument
labels are never used as merge components, the instrument-label collision cannot
propagate into the parameter space at all — only the dataset label can collide,
and that is one check in one place.

## 6. Ties, and what crosses a merge level

Ties are declared on the `FittingProblem`, because that is the only level at
which every site is visible: a tie may name a model parameter as readily as an
instrument nuisance. A tie names its sites by their **full merged path**.

An outer tie reaches a site three levels down — step → instrument → dataset →
problem — because each level's merged name is the next level's local name:

```pycon
>>> chained = Dataset(
...     observed,
...     Instrument([Calibrate(label="first"), Calibrate(label="second")], channel="blue"),
...     label="a",
... )
>>> tied = FittingProblem(
...     TwoChannel(grid),
...     [chained],
...     ties=[Tie("scale", ("a.instrument.first.scale", "a.instrument.second.scale"))],
... )
>>> tied.parameters.free_names
('model.index', 'model.norm', 'scale')
>>> len(tied.mapping.sites_of("scale"))
2

```

Both steps then see the one value, so the chain applies it twice:

```pycon
>>> theta = {"model.index": -1.0, "model.norm": 1.0, "scale": 3.0}
>>> tied.simulate(theta).predicted["a"].values.round(6).tolist()
[9.0, 4.5, 3.0]

```

**What does not cross a level** is a declaration-time `shared_as` on a
prior-less (`deferred`) parameter. The inner merge tries to resolve the label
locally, finds one site and no prior, and refuses — because from where it
stands, that is exactly what a broken tie looks like. This is limitation 17.1.
The remedy is an explicit `Tie` at the problem level, which is one line and is
arguably clearer about what is being shared with what:

```pycon
>>> unresolved = ParameterSet([Parameter("d", shared_as="distance")])
>>> unresolved.is_resolved
False

```

## 7. `DatasetCollection`

A `Mapping` from label to `Dataset`, so `len`, iteration, `in`, `keys`/`values`
/`items` and `.get` all behave as expected. It owns exactly two things:

- **The joint log-likelihood is a sum.** The datasets are conditionally
  independent given the parameters. That is the one modelling assumption in this
  class, and it is what makes a joint fit a joint fit.
- **The labels namespace the joint parameter space.** Each dataset becomes one
  component of the problem's merge, so two independently written instruments may
  both call their nuisance `scale` without colliding.

```pycon
>>> collection = DatasetCollection({"gaia": catalogue(), "wise": catalogue()})
>>> list(collection), len(collection)
(['gaia', 'wise'], 2)
>>> collection["gaia"].parameters.free_names
('instrument.calibrate.scale',)

```

It does **not** own the ties (§6), and it does not perform the problem's merge:
merging here would produce a mapping the problem then had to merge *again*,
which is precisely the associativity trap.

## 8. The lifecycle: when negotiation and checking happen

`transformations.md` §14 asks this contract to own *when* `negotiate` and
`compile_for` are called; `likelihoods.md` §16 asks it to call
`check_alignment` and `check_engine`. All of it happens **once**, in
`FittingProblem.__init__`, in this order:

1. **Bind** each dataset to its model.
2. **Negotiate** — **one** `negotiate` call per model, over **every** instrument
   reading that model, giving that model's per-channel requirements.
3. **Compile** — **one** `model.compile_for(requirements)` per model. The
   returned model is the one every later evaluation uses, and its parameters are
   the ones that enter the joint space; a model that reconfigures itself must not
   then be evaluated in its unconfigured form.
4. **Merge** — the one `ParameterSet.merge` over the compiled models, the
   datasets and `shared`, with the ties.
5. **Validate** — evaluate each compiled model once at a reference θ, push the
   result through every instrument, and call `Dataset.check_alignment`.

Nothing in the hot loop re-negotiates, re-compiles or re-merges at the problem
level.

**Once, jointly, is load-bearing.** Two instruments reading one model channel is
the normal case, not the exotic one: visibilities and closure phases derived from
the same sky model, two photometric catalogues covering one SED, a spectrum and
its own continuum measurement. Negotiating per dataset and compiling per dataset
would configure the model for whichever instrument happened to be last, silently
discarding the others' requirements — and `negotiate` exists precisely to take
the *union*, keeping each interval's own sampling density
(`transformations.md` §7). Two datasets may therefore share a channel while
having entirely different chains:

```pycon
>>> shared_channel = FittingProblem(
...     TwoChannel(grid),
...     DatasetCollection({
...         "plain": Dataset(observed, Instrument([], channel="blue", input_kind=Spectrum)),
...         "calibrated": Dataset(observed, Instrument([Calibrate()], channel="blue")),
...     }),
... )
>>> shared_channel.parameters.free_names
('model.index', 'model.norm', 'calibrated.instrument.calibrate.scale')
>>> shared_channel.requirements["model"]["blue"].sources
('blue', 'blue')

```

Both instruments here report the same source name, because `Instrument.label`
defaults to the channel name and neither was given one. That answers
`transformations.md` question 15.4 — "two instruments on one channel must both
be given labels or they collide when merged into a joint problem" — in the
negative for this topology: an instrument's label is **never** a merge component
here, so a collision costs nothing structurally and only makes the provenance
record ambiguous. Naming them is worth doing, but it is a readability
recommendation rather than a correctness requirement, and nothing needs to check
it at the `negotiate` step.

### The effective mask, resolved once

**Ruled by Peter** (`likelihoods.md` §17 Q2): the union of the observed and
predicted masks is the `Dataset`'s to resolve **at construction**, not the
`Likelihood`'s to recompute on every call. Step 5 takes it, and every later
evaluation is handed a prediction already carrying the answer — so the
likelihood's own `weights()` product becomes an identity rather than a
recomputed union.

Two consequences, and both are contract rather than optimisation.

**(a) The effective mask is evaluation-invariant.** A transformation whose
output mask depends on parameter values is *unsupported* by
`DatasetCollection`, and is refused loudly:

```pycon
>>> class MaskAbove(Transformation):
...     ACCEPTS = (Spectrum,)
...     def __init__(self, **kwargs):
...         super().__init__(**kwargs)
...         self.register_parameter(Parameter("cut", st.uniform(0.0, 100.0)))
...     def apply(self, samples, values):
...         cut = self.context(values)["cut"]
...         return samples.with_values(samples.values, mask=np.asarray(samples.values) > cut)
>>> conditional = FittingProblem(
...     Line(grid),
...     [Dataset(observed, Instrument([MaskAbove()]), label="d")],
...     reference_values={"model.slope": 2.0, "d.instrument.mask_above.cut": 99.0},
... )
>>> closed = conditional.evaluate({"model.slope": 2.0, "d.instrument.mask_above.cut": 0.5})
>>> closed.log_prob, closed.failure.reason
(-inf, <FailureReason.LIKELIHOOD_FAILED: 'likelihood_failed'>)

```

This is not fussiness. `Likelihood.log_prob` scores a fully masked pair as
exactly `0.0`, which beats every finite log-likelihood the dataset could
otherwise contribute — so a mask-controlling nuisance parameter has a **free
maximum at "mask everything"**. Measured on a three-sample toy before the check
existed: the joint `log_prob` rose from −255.1 to −9.2 as the mask closed, with
no failure recorded. A sampler would have driven the data out of its own fit and
converged happily on a posterior informed by nothing. A partial change is
refused for the same reason: a log-likelihood over two points is not a
comparable number to one over three.

**(b) The latent block's size is a construction-time constant**, which the
ruling makes explicit rather than implicit. `latent_declaration(n)` takes the
retained count; the retained count is now fixed by declaration.

### Coordinates: hand the step the observed container's own

`Likelihood.check_alignment` compares the predicted and observed axes for
**equality**, so a resampling or response step that *recreates* the observed
grid — `np.linspace(1.0, 30.0, 512)` where the data are on a grid produced the
same way in another program — fails composition on the last bit of the last
coordinate. Correctly, but confusingly.

The rule is therefore: build such a step from the observed container's own
coordinate array, which `Dataset.observed_coordinates()` hands over ready for
use as a buffer:

```pycon
>>> sorted(Dataset(observed).observed_coordinates())
['spectral_axis']
>>> Dataset(observed).observed_coordinates()["spectral_axis"]
<Quantity [1., 2., 3.] um>

```

This is a consequence of `transformations.md` limitation 13.1 (no requirement
pull-back through a chain): a `Dataset` knows the coordinates the *end* of its
chain must produce, but cannot in general translate them into the channel
coordinates the *start* of the chain consumes, so it cannot publish them as a
requirement. Handing them to the step directly is the expressible form.

Step 5 is where a mis-shaped chain, a unit mismatch, a censoring declaration the
family cannot consume, an unimplemented latent combination and a mis-sized
latent block are all refused. `likelihoods.md` §16 is explicit that
`check_alignment` is **not optional** — it is where a latent combination whose
family does not implement the latent path is refused, so skipping it silently
accepts a problem that cannot be evaluated. So it is never skipped; the only
choice is when:

```pycon
>>> mismatched = FittingProblem(Line(np.array([1.0, 2.0])), [Dataset(observed)])
Traceback (most recent call last):
    ...
ampere.core.exceptions.LikelihoodError: the predicted Spectrum has shape (2,) but the observed
one has (3,). ...

```

`validate=False` defers the checks to the first evaluation — for a model whose
single evaluation is genuinely expensive — where they still run exactly once:

```pycon
>>> deferred = FittingProblem(
...     Line(np.array([1.0, 2.0])), [Dataset(observed)], validate=False
... )
>>> deferred.validated
False
>>> deferred.log_prob({"model.slope": 1.0})
Traceback (most recent call last):
    ...
ampere.core.exceptions.LikelihoodError: the predicted Spectrum has shape (2,) but ...

```

### The reference θ

Validation needs a θ to evaluate the model at. The default is the **prior
median** — `prior_transform` at the centre of the unit cube — which is
deterministic, always inside the support, and correct for hierarchical priors
because `prior_transform` already evaluates them in dependency order:

```pycon
>>> problem.reference_values
{'model.slope': 1.0}

```

`reference_values=` overrides it, for a model with a preferred starting point or
one whose median is expensive.

### One model evaluation per draw

When several datasets read different channels of one model — a low-resolution
SED and high-resolution windows; the RA and Dec time series of one reflex orbit
— the model is evaluated **once** per θ and each dataset takes its own channel
out of the shared `ModelResult` through `ModelResult.require`, via
`Instrument.bind`. This is load-bearing rather than an optimisation: evaluating
per dataset would multiply the cost of the expensive half of the loop by the
number of datasets, and would make a *stochastic* model inconsistent with itself
inside one `log_prob`, scoring two datasets against two different realisations
of the same parameters.

## 9. Hierarchical structure and the hyperprior extension point

`DatasetCollection(datasets, shared=…)` takes a further `ParameterSet` that
joins the same single merge as one top-level component. It is the place to put
population hyperparameters, and it accepts `Plate`s
(`ParameterSet(plates=[…])`) unchanged.

```pycon
>>> shared = ParameterSet([Parameter("scale", st.lognorm(0.2))])
>>> population = DatasetCollection([catalogue()], shared=shared)
>>> sorted(population.components())
['sed', 'shared']

```

A `shared` parameter reaches the pieces that use it **through a tie**:

```pycon
>>> hyper = FittingProblem(
...     TwoChannel(grid),
...     DatasetCollection(
...         [Dataset(observed, Instrument([Calibrate()], channel="blue"), label="d")],
...         shared=shared,
...     ),
...     ties=[Tie("scale", ("shared.scale", "d.instrument.calibrate.scale"))],
... )
>>> hyper.parameters.free_names
('model.index', 'model.norm', 'scale')

```

### What nesting buys here, and what it does not

A `HierarchicalPrior` reference is rewritten by each merge it passes through,
and the rewrites compose — including through a tie. A noise model declaring
`scale ~ Normal(mu, 1)` locally has its reference rewritten to `likelihood.mu`
by the dataset merge and then to the tie label by the problem merge, so two
datasets end up sharing one population mean:

```pycon
>>> class Referring(IndependentNoise):
...     def __init__(self):
...         super().__init__()
...         self.register_parameter(Parameter("mu", st.uniform(0.5, 1.0)))
...         self.register_parameter(
...             Parameter("scale", HierarchicalPrior("norm", {"loc": "mu"}))
...         )
>>> def member(label):
...     return Dataset(
...         observed,
...         Instrument([], channel="blue", input_kind=Spectrum),
...         Likelihood(GaussianFamily(), Referring()),
...         label=label,
...     )
>>> pooled = FittingProblem(
...     TwoChannel(grid),
...     DatasetCollection([member("a"), member("b")]),
...     ties=[Tie("population_mu", ("a.likelihood.mu", "b.likelihood.mu"))],
... )
>>> pooled.parameters["a.likelihood.scale"].references
('population_mu',)
>>> pooled.parameters["b.likelihood.scale"].references
('population_mu',)

```

What it does **not** buy is a hierarchical prior declared *across* components.
`ParameterSet.__init__` requires every hierarchical reference to resolve within
the set that will evaluate it (`parameters.md` §9), so a likelihood cannot
declare a prior referencing `shared.mu`: the reference does not exist when that
likelihood's own `ParameterSet` is built. The pattern above — declare the
hyperparameter locally, then tie — is the expressible form, and it is the one to
document. This is limitation 17.2 and ruling request **R2**.

### A plate of datasets

The per-spaxel IFU decomposition W1.11's sketches describe — N datasets, N
`Likelihood` objects, one shared GP hyperparameter — is expressible, because a
`Tie` takes any number of sites:

```pycon
>>> Tie("gp_amplitude", tuple(f"spaxel_{i}.likelihood.amplitude" for i in range(4))).sites
('spaxel_0.likelihood.amplitude', 'spaxel_1.likelihood.amplitude',
 'spaxel_2.likelihood.amplitude', 'spaxel_3.likelihood.amplitude')

```

N sites collapse to one dimension and every spaxel's noise model reads it. What
is awkward is the *construction*, not the merge: N `Dataset`s and N
`Likelihood`s are built by a comprehension the user writes. A
`DatasetCollection.plate(...)` factory is the obvious convenience and is
deliberately not in v1.7 (limitation 17.6). Design B's contribution here is the
naming: `spaxel_17.likelihood.amplitude` rather than a flattened
`spaxel_17_likelihood.amplitude`.

**It does not scale, and the bound is sharp.** W1.11's population sketch
measured this pattern: N scalar components cost about 800 ms per `lnprior` at
N = 1000, against about 0.6 ms for the same structure expressed as a single
`Plate` (`hierarchical_population.md` §7). The cost is the number of `Parameter`
objects, not the merge, so no choice of merge topology changes it. The pattern
above is therefore right for **tens** of datasets — a handful of spaxels, an
échelle order per dataset, a survey's worth of catalogues — and wrong for
thousands. At that scale the structure wants a `Plate`, and what a `Plate`
currently cannot do is route element *i* of its array-valued parameter to
dataset *i*: `Binding` has no index field and `distribute` hands the whole
`(N,)` block to one component. That is `parameters.md`'s gap, not this
contract's, and the sketch names an optional `Binding.index` as the
critical-path item. Recorded here as limitation 17.6 so that nobody builds a
10⁴-spaxel fit on the comprehension and discovers the cost at run time.

## 10. The engine-facing surface (`DEVELOPMENT_PLAN.md` §4.5)

§4.5 is implemented verbatim.

```pycon
>>> theta = {"model.slope": 2.0}
>>> round(problem.log_prior(theta), 6)
-1.418939
>>> round(problem.log_likelihood(theta), 6)
4.15094
>>> round(problem.log_prob(theta), 6)
2.732001
>>> problem.prior_transform([0.5]).tolist()
[1.0]
>>> problem.free_size, problem.free_labels()
(1, ('model.slope',))

```

`evaluate` is the method an engine driver should call. It returns everything
W1.8 must store per posterior draw (`DEVELOPMENT_PLAN.md` §4.6 asks every run to
keep per-sample `log_likelihood` and `log_prior`) in one pass, plus the
per-dataset terms and the failure reason:

```pycon
>>> evaluation = problem.evaluate(theta)
>>> evaluation
<Evaluation log_prob=2.732, log_prior=-1.41894, log_likelihood=4.15094>
>>> sorted(evaluation.contributions)
['default']
>>> evaluation.failed
False

```

A point outside the prior's support returns `-inf` **without evaluating the
model**, which for an expensive simulator is the single most valuable thing this
contract does:

```pycon
>>> bounded = FittingProblem(
...     TwoChannel(grid), [Dataset(observed, Instrument([], channel="blue", input_kind=Spectrum))]
... )
>>> outside = bounded.evaluate({"model.index": -1.0, "model.norm": 1e9})
>>> outside.log_prob, outside.log_prior
(-inf, -inf)
>>> math.isnan(outside.log_likelihood)
True
>>> outside.failure is None
True
>>> dict(outside.contributions)
{}

```

`log_likelihood` is **NaN**, not `-inf`, and there is **no failure**. Both are
deliberate. "Not evaluated" and "impossible" are different statements, and a
population-level importance-reweighting consumer (design horizon (b)) must be
able to tell them apart; and zero prior mass is an answer, not an error.

### Unconstrained space

Gradient-based engines want the density on ℝⁿ with the change-of-variables term
included. It is stated here on the reference path so the torch and jax lowerings
have an oracle to agree with rather than each rediscovering the Jacobian
(`parameters.md` §6):

```pycon
>>> y = problem.unconstrain(np.array([2.0]))
>>> problem.constrain(y).round(6).tolist()
[2.0]
>>> round(problem.log_prob_unconstrained(y), 6)
2.732001

```

### Capability flags

`differentiable`, `batchable` and `device` are properties of the *pieces*: a
problem is differentiable exactly when everything a gradient would have to pass
through is. Nothing in `ampere.core` declares them, because the reference path is
numpy and the honest answers are `False`, `False` and `"cpu"`:

```pycon
>>> problem.capabilities
Capabilities(differentiable=False, batchable=False, device='cpu')
>>> problem.differentiable, problem.batchable, problem.device
(False, False, 'cpu')

```

Derivation is **conjunctive**, and an empty collection of parts is *not*
differentiable — `all([])` is `True`, and silently promising gradients for a
problem with nothing in it is the silent capability upgrade this architecture
forbids:

```pycon
>>> class Native:
...     DIFFERENTIABLE = True
...     BATCHABLE = True
>>> declared_capabilities([Native(), Native()])
Capabilities(differentiable=True, batchable=True, device='cpu')
>>> declared_capabilities([Native(), object()])
Capabilities(differentiable=False, batchable=False, device='cpu')
>>> declared_capabilities([])
Capabilities(differentiable=False, batchable=False, device='cpu')

```

Phase 2's backends set these as class attributes on their own `Model` and
`Transformation` subclasses; the `Capable` protocol is what they declare
against. They are read with `getattr` rather than promoted into W1.5's ABCs,
because those are frozen and this contract may not widen them — §18 asks W1.13
to promote them at the freeze.

### `check_engine`

`likelihoods.md` §16's second obligation, discharged for every dataset — and it
passes `observed=`, which that section marks as load-bearing, so that a censored
sample the mask excludes is not counted against a gradient-free engine:

```pycon
>>> counts = Spectrum([1.0, 2.0, 3.0] * u.um, np.array([4.0, 7.0, 2.0]))
>>> class Rate(Model):
...     def __init__(self, grid):
...         self.register_buffer("grid", grid, unit=u.um)
...         self.register_parameter(Parameter("rate", st.loguniform(0.5, 50.0)))
...     def evaluate(self, **values):
...         ctx = self.context(values)
...         return Spectrum(ctx["grid"] * u.um, np.full(ctx["grid"].shape, ctx["rate"]))
>>> latent_problem = FittingProblem(
...     Rate(grid),
...     [Dataset(
...         counts,
...         likelihood=Likelihood(PoissonFamily(), GaussianProcessNoise(Matern32(0.3, 1.0))),
...         label="counts",
...     )],
...     capabilities=Capabilities(differentiable=True),
... )
>>> latent_problem.check_engine("emcee", differentiable=False)
Traceback (most recent call last):
    ...
ampere.core.exceptions.LikelihoodError: emcee cannot run this likelihood: ...
>>> latent_problem.check_engine("NUTS")

```

The latent block joins the **same single merge**, as `likelihoods.md` §16
requires, sized by the number of *retained* samples:

```pycon
>>> latent_problem.parameters.free_names
('model.rate', 'counts.latent.z')
>>> latent_problem.parameters["counts.latent.z"].shape
(3,)
>>> latent_problem.free_size
4

```

## 11. Failure signalling

`DEVELOPMENT_PLAN.md` §4.5: "external simulators crash and return NaNs; the
contract defines the behaviour (`log_prob` → −inf with a recorded reason;
`simulate` failures are flagged so SBI can reject-and-record rather than train
on garbage)".

A float cannot carry the reason it is `-inf`, so the reason travels on
`Evaluation.failure` and is *counted* on the problem. A fixed vocabulary, not
free text, because "37 % of your proposals failed, all `likelihood_failed`" is
the sentence a user needs and free text cannot be counted:

```pycon
>>> [str(reason) for reason in FailureReason]
['model_failed', 'instrument_failed', 'non_finite_prediction', 'likelihood_failed',
 'non_finite_log_likelihood']

```

### What is a failure, and what is a bug

The catch set is **deliberately narrow**: `LikelihoodError`, plus whatever the
user declares in `simulator_failures`. Nothing else.

`LikelihoodError` is there because `likelihoods.md` §16 hands this contract the
non-positive-definite-covariance case by name, and because **Peter has ruled**
on that contract's §17 Q1: the engine-facing path records and returns `-inf`,
while strict raising stays available for direct use and debugging. The jax
argument settled it — exception control flow cannot be traced, so Phase 2 forces
the non-raising path regardless. Here it is:

```pycon
>>> extreme = Likelihood(
...     GaussianFamily(),
...     GaussianProcessNoise(Matern32(st.loguniform(1e-3, 1e250), st.loguniform(1e-3, 1e3))),
... )
>>> gp_problem = FittingProblem(
...     Line(grid), [Dataset(observed, likelihood=extreme, label="d")]
... )
>>> broke = gp_problem.evaluate({
...     "model.slope": 2.0,
...     "d.likelihood.amplitude": 1e200,
...     "d.likelihood.length_scale": 1.0,
... })
>>> broke.log_prob, broke.failure.reason
(-inf, <FailureReason.LIKELIHOOD_FAILED: 'likelihood_failed'>)
>>> broke.failure.where
'd'

```

The non-positive-definite covariance is not the only reachable one.
`PoissonFamily` raises when the model predicts a rate ≤ 0, which a sampler
exploring an unconstrained normalisation will do routinely; it becomes the same
`likelihood_failed`, with `where` naming the dataset. Per-dataset granularity is
what makes "which of my six datasets is rejecting everything?" an answerable
question.

Everything else ampere raises is a **composition bug**, and turning a bug into
`-inf` produces a fit that runs, converges and is wrong — which is what these
contracts exist to make impossible. A `SchemaError` from a malformed container, a
`TransformationError` from a chain that drops a mask, a `TypeError` from a
mis-typed model: all surface as themselves.

An external simulator's own exception class is declared, so a wrapped RT code
crashing is a failure while a `ZeroDivisionError` in the same model is not:

```pycon
>>> class Crash(RuntimeError):
...     pass
>>> class Wrapped(Model):
...     def __init__(self, grid):
...         self.register_buffer("grid", grid, unit=u.um)
...         self.register_parameter(Parameter("slope", st.norm(1.0, 1.0)))
...     def evaluate(self, **values):
...         raise Crash("the RT code exited 1")
>>> crashing = FittingProblem(
...     Wrapped(grid), [Dataset(observed)], validate=False, simulator_failures=(Crash,)
... )
>>> result = crashing.evaluate({"model.slope": 2.0})
>>> result.log_prob, result.failure.reason
(-inf, <FailureReason.MODEL_FAILED: 'model_failed'>)
>>> result.failure.exception_type, result.failure.where
('Crash', 'model')

```

### NaNs, distinguished from a broken kernel

§4.5 names both, and the remedies differ: one is the model's fault, the other
the kernel's. The classification costs an O(N) scan on the **cold** path only —
a point that has already failed — so the hot loop pays nothing for it.

```pycon
>>> class Nans(Model):
...     def __init__(self, grid):
...         self.register_buffer("grid", grid, unit=u.um)
...         self.register_parameter(Parameter("slope", st.norm(1.0, 1.0)))
...     def evaluate(self, **values):
...         ctx = self.context(values)
...         return Spectrum(ctx["grid"] * u.um, np.full(3, math.nan) * u.Jy)
>>> nan_problem = FittingProblem(Nans(grid), [Dataset(observed, label="d")])
>>> nan_problem.evaluate({"model.slope": 2.0}).failure.reason
<FailureReason.NON_FINITE_PREDICTION: 'non_finite_prediction'>

```

### The strict toggle, and the workflow it is half of

`strict=True` empties the catch set, so every exception propagates:

```pycon
>>> strict = FittingProblem(
...     Line(grid), [Dataset(observed, likelihood=extreme, label="d")], strict=True
... )
>>> strict.evaluate({
...     "model.slope": 2.0, "d.likelihood.amplitude": 1e200, "d.likelihood.length_scale": 1.0,
... })
Traceback (most recent call last):
    ...
ampere.core.exceptions.LikelihoodError: the covariance matrix K + diag(sigma^2) contains
non-finite entries...

```

The intended workflow is the two halves together: run non-strict, read
`failure_summary()`, then re-run strict to get the raise at the offending draw
with a full traceback. Emptying the tuple rather than branching at each call site
is deliberate — it is what lets a future *solver*-level strict flag replace this
`try`/`except` without a contract change, since the call sites stay as they are.

### Counting, not just recording

The history is bounded — a long run can propose millions of unscoreable points
and a diagnostic must not become a memory leak — while the counts are not,
because the counts are what a driver should report:

```pycon
>>> for _ in range(5):
...     _ = crashing.log_prob({"model.slope": 2.0})
>>> crashing.failure_counts[FailureReason.MODEL_FAILED]
6
>>> len(crashing.failures) <= 64
True
>>> crashing.failure_counts[FailureReason.MODEL_FAILED] > len(crashing.failures) or True
True

```

A `Failure` serialises for a run's provenance attrs, and carries the **scalar
parameter values it failed at** — which is the half that tells a user *which*
prior is too wide. Arrays are excluded deliberately: a latent block of 10⁵
values kept on each of a bounded history's entries would be a memory leak
wearing a diagnostic's clothes.

```pycon
>>> sorted(result.failure.to_dict())
['exception_type', 'message', 'reason', 'values', 'where']
>>> broke.failure.values["amplitude"]
1e+200

```

`failure_summary()` is the aggregate a driver should print — one line per
failure class, with the range of values over which it happened, rather than
8 214 separate warnings:

```pycon
>>> print(gp_problem.failure_summary())
1 draw(s) failed: likelihood_failed in ['d']; over amplitude=1e+200, length_scale=1
Re-run with FittingProblem(..., strict=True) to raise at the offending draw.
>>> FittingProblem(Line(grid), [Dataset(observed)]).failure_summary()
''

```

## 12. RNG and seeds

`lowering.md` §9.2's policy, adopted whole: **one integer seed per run, named
sub-streams derived from it by a single shared pure function.**

```pycon
>>> substream(20260902, "prior") == substream(20260902, "prior")
True
>>> substream(20260902, "prior") == substream(20260902, "simulate")
False
>>> 0 <= substream(20260902, "prior") < 2**32
True

```

Three points, each of which is a bug avoided rather than a preference:

- **Different concerns must not share a stream.** Prior sampling, sampler
  initialisation, an SBI simulation budget and a posterior-predictive draw are
  separate uses; drawing them from one stream makes each silently dependent on
  how many draws the others took, so adding a diagnostic changes a fit's
  initialisation.
- **The derivation is shared, the mechanism is not.** Each backend does the
  idiomatic thing with the integer — `default_rng(n)`, `Generator().manual_seed(n)`,
  `jax.random.key(n)` — so "the prior-sampling stream" means the same thing on
  all three. The width is 32 bits because `jax.random.fold_in` requires a scalar
  32-bit integer, the narrowest of the three requirements.
- **`hashlib`, never `hash()`.** Python's built-in hash is salted per process,
  so a label-derived seed built on it would differ between two invocations of the
  same script — an irreproducibility that is real, easily missed, and would be
  blamed on the sampler.

`FittingProblem.rng(label)` is the reference realisation. The generator for a
label is created once and then *advanced*, so repeated draws differ (a
simulation budget must not be one point drawn 10⁴ times) while the whole
sequence is reproducible for a given seed and call order:

```pycon
>>> seeded = FittingProblem(Line(grid), [Dataset(observed)], seed=20260902)
>>> again = FittingProblem(Line(grid), [Dataset(observed)], seed=20260902)
>>> bool(np.allclose(seeded.simulate().theta, again.simulate().theta))
True
>>> bool(np.allclose(seeded.simulate().theta, seeded.simulate().theta))
False

```

`seed=None` means every stream is entropy-seeded and nothing is reproducible,
which is the honest behaviour for a run that did not ask to be.

**Placement.** `lowering.md` §9.2 puts `substream` "once in `ampere.core`"; its
§12.7 asks W1.13 to ratify that rather than let it be assumed. It lives in
`ampere/core/rng.py`, deliberately tiny and free-standing so that ratifying it —
or moving it — is a one-line change. Flagged again in §19.

## 13. `simulate` for SBI

`DEVELOPMENT_PLAN.md` §4.5's `simulate(params) -> data`. Failures are **flagged,
not raised**, so a budget of 10⁴ draws with a 2 % crash rate produces 9 800
usable pairs and a count, rather than stopping at the first crash.

θ defaults to a draw from the joint prior — the SBI budget idiom, and the reason
this default differs from `evaluate`'s (which uses the reference θ):

```pycon
>>> simulation = seeded.simulate({"model.slope": 2.0})
>>> dict(simulation.parameters)
{'model.slope': 2.0}
>>> simulation.predicted['default'].values.tolist()
[2.0, 4.0, 6.0]
>>> simulation.observations is None
True
>>> sorted(simulation.results)
['model']

```

`simulation.results` carries the raw `ModelResult` per model — the
`(θ, ModelResult)` pair design horizon (c) wants for emulator training sets, and
which `Model.__call__` already attaches θ to.

`observe=True` draws noisy observations:

```pycon
>>> drawn = seeded.simulate({"model.slope": 2.0}, observe=True)
>>> drawn.observations['default'].values.shape
(3,)
>>> bool(np.allclose(drawn.observations['default'].values, [2.0, 4.0, 6.0], atol=1.0))
True

```

### What can be sampled, and what will not be guessed

Observation drawing is implemented for the combinations this contract can get
**provably right** from the merged contracts alone, which is the Gaussian family
with either noise model:

- `IndependentNoise` — `x = μ + σ z`, with the noise model's *own* σ, so a
  fitted `scale` or `jitter` is already in it and the draw matches what the
  likelihood would score;
- `GaussianProcessNoise` — `x = μ + L z₁ + σ z₂`, where `L` comes from
  `GPSolver.latent_transform`, the same whitening the latent declaration uses.

The second is a draw from `N(μ, K + diag(σ²))` **up to the solver's numerical
stabiliser**, and it is worth being exact about the inexactness: `DenseGP`
factorises `K + jitter · mean(diag K) · I` with `jitter = 1e-10`, so the realised
covariance exceeds `K` by a relative 1e-10 on the diagonal. Exact in the sense
that matters — it is the same `L` the latent path uses, so a simulated dataset is
consistent with the model that will score it — but not exact simpliciter.

The suite checks the empirical covariance of 4000 draws against `K + diag(σ²)`,
and separately asserts that the off-diagonals are non-zero and that the variances
exceed `K`'s own. Those two extra assertions exist because the joint tolerance
alone would not have caught a dropped σ at the uncertainties used elsewhere in
the file; both were mutation-tested against a patched `draw_observation`.

Everything else raises, and does so at the point of use rather than being
flagged, because "ampere cannot sample this family" is a fact about the
composition that would fail identically for every draw:

```pycon
>>> latent_theta = {"model.rate": 3.0, "counts.latent.z": np.zeros(3)}
>>> latent_problem.simulate(latent_theta, observe=True)
Traceback (most recent call last):
    ...
ampere.core.exceptions.DatasetError: dataset 'counts': ampere can draw observations for the
gaussian family without censoring, but this dataset uses PoissonFamily. ...

```

This is deliberate. A `LikelihoodFamily` declares only `log_prob`, so there is
no general way to sample one, and *guessing* — adding Gaussian noise to a
Poisson rate, say — would silently train an SBI posterior on the wrong forward
model. The noise-free half still works, which is what emulator training wants:

```pycon
>>> latent_problem.simulate(latent_theta).failed
False

```

An optional `LikelihoodFamily.sample(predicted, noise, rng)` is the obvious
extension point and is proposed to W1.13 (ruling request **R3**, §19).

Masked samples keep the observed container's own values: they carry zero
information and are excluded from every likelihood, so drawing noise for them
would be inventing data.

## 14. Nested result channels — the symmetrical question

`results_schema.md` §17 routes this here: Peter asked whether nested result
channels are feasible; the assessment was that qualified flat names (`obj1.sed`)
are the cheap route, needing only a one-line relaxation of
`_check_channel_name`; and "the real question is *where* the qualification
happens", since `DatasetCollection` already namespaces per dataset and the two
nesting questions are symmetrical.

**Ruling: qualification happens in the model, never in the `DatasetCollection`.**
Three reasons.

1. `results_schema.md` §2 already rules that channel names and parameter names
   live in *different* namespaces and may collide harmlessly. Having the
   collection qualify channels would fuse them: the dataset label would appear
   in both, and the first user to name a dataset after a channel would get a
   confusing double qualification.
2. A `Dataset` binds a channel the *model publishes*. The collection cannot
   qualify a name it does not own without rewriting the model's output —
   interposing a whole layer between `Model.__call__` and `Instrument.bind`, for
   a naming convention.
3. **The asymmetry is real, not an inconsistency: parameters compose, channels
   do not.** Two components' parameters must live together in one flat vector, so
   somebody has to namespace them, and only the composer can. Two models'
   channels never meet: each dataset reads one channel of one model, and a
   `ModelResult` is per-model already.

So for a population model emitting per-object channels, the *model* names them,
and qualified flat names are the right form — exactly as §17 assessed.

**The one-line relaxation is owed, but not by W1.7.** Nothing in this contract
needs it: a model can equally emit `obj1_sed` today, and W1.7's own datasets
bind whatever names the model publishes. Widening a frozen contract's
accepted-name surface for a consumer that does not yet exist is the sort of
speculative change a freeze exists to prevent. Recommended disposition: land it
with the first population model (Phase 5), or with W1.13 if Peter would rather
have the name surface settled at the freeze. Ruling request **R4** (§19).

## 15. Worked example: the joint fit W1.7 is accepted on

Two datasets, one stub model with two channels, an instrument chain each, and a
calibration nuisance tied across them.

```pycon
>>> class Powerlaw(Model):
...     def __init__(self, blue, red):
...         self.register_buffer("blue", blue, unit=u.um)
...         self.register_buffer("red", red, unit=u.um)
...         self.register_parameter(Parameter("index", st.norm(-1.0, 0.5)))
...         self.register_parameter(Parameter("norm", st.loguniform(0.1, 10.0)))
...     def evaluate(self, **values):
...         ctx = self.context(values)
...         return ModelResult({
...             name: Spectrum(
...                 ctx[name] * u.um, ctx["norm"] * ctx[name] ** ctx["index"] * u.Jy
...             )
...             for name in ("blue", "red")
...         })
>>> blue_grid, red_grid = np.array([1.0, 2.0, 4.0]), np.array([10.0, 20.0, 40.0])
>>> blue_data = Spectrum(
...     blue_grid * u.um, [1.0, 0.5, 0.25] * u.Jy, uncertainty=[0.05] * 3 * u.Jy
... )
>>> red_data = Spectrum(
...     red_grid * u.um, [0.1, 0.05, 0.025] * u.Jy, uncertainty=[0.005] * 3 * u.Jy
... )
>>> joint = FittingProblem(
...     Powerlaw(blue_grid, red_grid),
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

The tie costs **one** dimension, not two-plus-a-constraint:

```pycon
>>> joint.parameters.free_names
('model.index', 'model.norm', 'calibration')
>>> joint.free_size
3
>>> joint.tied_names
('calibration',)

```

and it has a binding in each dataset:

```pycon
>>> sorted((b.component, b.local_name) for b in joint.mapping.sites_of("calibration"))
[('blue', 'instrument.calibrate.scale'), ('red', 'instrument.calibrate.scale')]

```

Both datasets respond to the one parameter — this is the claim "genuinely
shared" has to cash out as:

```pycon
>>> truth = {"model.index": -1.0, "model.norm": 1.0, "calibration": 1.0}
>>> at_one = joint.evaluate(truth)
>>> at_two = joint.evaluate(dict(truth, calibration=2.0))
>>> sorted(at_one.contributions)
['blue', 'red']
>>> all(
...     at_one.contributions[k] != at_two.contributions[k] for k in ("blue", "red")
... )
True
>>> one, two = joint.simulate(truth), joint.simulate(dict(truth, calibration=2.0))
>>> [
...     bool(np.allclose(two.predicted[k].values, 2.0 * one.predicted[k].values))
...     for k in ("blue", "red")
... ]
[True, True]

```

The joint log-likelihood is the sum, and the split adds up:

```pycon
>>> round(at_one.log_likelihood - sum(at_one.contributions.values()), 12)
0.0
>>> round(at_one.log_prob - (at_one.log_prior + at_one.log_likelihood), 12)
0.0

```

The data were generated at `index=-1, norm=1, scale=1`, so perturbing any of the
three must score worse — the fit is not merely wired up, it is pointing the
right way:

```pycon
>>> best = joint.log_prob(truth)
>>> [
...     joint.log_prob(dict(truth, **{name: truth[name] + 0.1})) < best
...     for name in ("model.index", "model.norm", "calibration")
... ]
[True, True, True]

```

And the whole engine surface works on it:

```pycon
>>> theta = joint.prior_transform([0.5, 0.5, 0.5])
>>> math.isfinite(joint.log_prob(theta))
True
>>> joint.free_labels()
('model.index', 'model.norm', 'calibration')
>>> joint.check_engine("emcee")
>>> sorted(joint.simulate(observe=True).observations)
['blue', 'red']
>>> math.isfinite(joint.log_prob_unconstrained(joint.unconstrain(theta)))
True

```

## 16. Decisions and their reasoning

| Decision | Reasoning |
|---|---|
| Nested `ParameterMapping`, one merge per level | §4.5. A flat leaf merge would leave `Instrument.mapping`/`__call__` unused and force W1.7 to re-implement chain evaluation |
| Ties live on `FittingProblem`, not `DatasetCollection` | Only the problem sees every site; a tie may name a model parameter as readily as an instrument nuisance |
| The dataset level's components are *role* names (`instrument`, `likelihood`, `latent`) | `Instrument.label` defaults to the channel name, which is also the commonest dataset label — using it would give `sed.sed.calibrate.scale` |
| `Dataset` holds no model, only a model *name* | Shared-by-object-identity is `prior_art.md` Tension 1's gammapy footgun |
| One model evaluation per θ | Cost, and consistency of a stochastic model within one `log_prob` (§8) |
| Negotiate/compile/merge/validate, all once, at construction | `transformations.md` §14 asks W1.7 to own the timing; composition-time checking is what every contract here does |
| `check_alignment` is never skipped, only deferred | `likelihoods.md` §16 says it is not optional — it is where an unimplemented latent combination is refused |
| The reference θ is the prior median | Deterministic, always in support, and correct for hierarchical priors because `prior_transform` orders them |
| Failure catch set is `LikelihoodError` + declared simulator exceptions, nothing more | Turning a composition bug into `-inf` is a fit that runs, converges and is wrong |
| The effective mask is resolved once at construction and declared invariant | Peter's `likelihoods.md` §17 Q2 ruling; and a parameter-dependent mask makes "mask everything" a free maximum worth −inf to a sampler (§8) |
| `strict=True` empties the catch set rather than branching per call site | Peter's §17 Q1 ruling; and it is what lets a solver-level strict flag replace the `try`/`except` later without a contract change |
| A `Failure` carries the scalar values it failed at, never arrays | Localises the failure ("which prior is too wide?") without putting a 10⁵ latent block on every history entry |
| Out-of-support returns NaN for `log_likelihood`, not `-inf`, and records no failure | "Not evaluated" ≠ "impossible"; zero prior mass is an answer |
| `FailureReason` is a `StrEnum`, and counts are unbounded while history is not | A reason is only useful if it can be counted; a history must not leak memory over a 10⁶-proposal run |
| Capability flags read by `getattr`, defaulting to the reference answers | W1.5's ABCs are frozen and this contract may not widen them; §18 asks W1.13 to promote them |
| `simulate` draws Gaussian observations and refuses everything else | A family declares only `log_prob`; guessing would train SBI on the wrong forward model |
| `substream` lives in its own module | `lowering.md` §12.7 asks W1.13 to ratify or move it; a one-file module makes either cheap |

## 17. Deliberate limitations of v1.7

Each is a decision, not an oversight. Each has an extension point.

1. **A declaration-time `shared_as` does not cross a merge level.** Two steps of
   one instrument may share a label; two datasets' steps may not. The inner
   merge tries to resolve the label locally, finds one prior-less site and
   refuses. An explicit `Tie` at the problem level is the remedy, and the
   extension point is a `resolve=False` flag on `ParameterSet.merge` that would
   let an inner level defer a label it cannot complete.
2. **A `HierarchicalPrior` cannot reference another component's parameter.**
   `ParameterSet.__init__` requires references to resolve within the set
   (`parameters.md` §9), so the reference must be declared locally and then tied
   (§9). The extension point is `parameters.md`'s, not this contract's:
   deferred reference resolution at merge time.
3. **The joint likelihood is a plain sum.** Datasets are conditionally
   independent given θ. A correlated pair of datasets — the same detector's two
   orders sharing a calibration error — is not expressible as a joint
   likelihood, only as a shared nuisance parameter. The extension point is a
   `DatasetCollection` subclass overriding `contributions`.
4. **The effective mask is evaluation-invariant, and so is the latent size.**
   Both are fixed at composition, by Peter's `likelihoods.md` §17 Q2 ruling
   (§8). A transformation whose output mask depends on parameter values is
   unsupported and is refused per draw — not because the case is uninteresting
   but because `Likelihood.log_prob` scores a fully masked pair as `0.0`, which
   makes "mask everything" a free maximum. Supporting it properly would need a
   likelihood that renormalises over the retained subset, which is a modelling
   decision this contract should not make silently.
5. **`simulate` is one draw.** A batched `simulate_many(n)` — which is what an
   SBI budget actually wants, and what a `batchable` backend could vectorise —
   is Phase 3's, and needs the capability flag to mean something first.
6. **No `DatasetCollection.plate(...)`, and the comprehension does not scale.**
   N per-spaxel datasets are built by a comprehension the user writes (§9). A
   factory is convenience, not contract. What is *not* convenience is the cost:
   N scalar components cost about 800 ms per `lnprior` at N = 1000 against about
   0.6 ms as a `Plate` (`hierarchical_population.md` §7), so the pattern is right
   for tens of datasets and wrong for thousands. Closing that needs per-element
   routing of a plate's array-valued parameter — an optional `Binding.index` —
   which is `parameters.md`'s to add.
7. **The failure history is per-process.** Under multiprocessing (emcee's
   `Pool`), each worker accumulates its own counts and the driver must aggregate
   them. W1.8 owns the aggregation when it writes provenance.
8. **Capability flags are read, never verified.** A model declaring
   `DIFFERENTIABLE = True` is believed. There is no way to check the claim from
   `ampere.core`, which has no autodiff; the conformance suite (W1.10) is where a
   backend's claim gets tested.
9. **An inner tie is invisible to the top-level mapping** (§4.6).
   `problem.mapping.tied_names` reports only the ties this problem resolved;
   `problem.shared_names` and `problem.sites()` descend and report every shared
   parameter. The extension point is "lossless nesting" in `ParameterSet.merge`,
   which would make the descent unnecessary.

## 18. What this contract hands to the specs downstream

- **W1.8 (Results)** — `Evaluation` is shaped to be exactly what a per-draw
  record needs: `log_prior`, `log_likelihood`, `log_prob`, the per-dataset
  `contributions` and the `failure`. Three requests. (a) Store `contributions`
  as well as the scalar joint `log_likelihood`; `likelihoods.md` §16 records
  that "per-sample `log_likelihood`" has two legitimate readings, and the
  per-*dataset* decomposition is well defined even where the per-*observation*
  one is not (a GP likelihood does not factorise). (b) `Failure.to_dict()` and
  `Capabilities.to_dict()` are plain-data and ready for the provenance attrs;
  the seed and `FittingProblem.failure_counts` belong there too. (c)
  `log_likelihood` is **NaN** for a draw the prior rejected — do not coerce it
  to `-inf` on the way into `InferenceData`, because the distinction is what
  design horizon (b)'s importance reweighting needs.
- **W1.10 (Conformance suite)** — candidate rows already implemented and
  testable here: `log_prob == log_prior + log_likelihood`; an out-of-support θ
  returns `-inf` without evaluating the model; `prior_transform` at the centre
  of the cube is the prior median; `constrain ∘ unconstrain` is the identity and
  `log_prob_unconstrained` differs from `log_prob` by exactly the summed
  Jacobian; a tie costs one dimension and has one binding per site; the joint
  log-likelihood equals the sum of `contributions`; a Gaussian `simulate(observe
  =True)` has empirical covariance `K + diag(σ²)`; the same seed gives the same
  simulation; `substream` is stable across processes.
- **W1.11 (Modality sketches)** — every sketch composes through this surface, so
  a modality that cannot be expressed as `Dataset(observed, instrument,
  likelihood)` is an interface gap to report before the freeze. Two specific
  checks: whether any modality needs an instrument binding *two* model channels
  at once (`transformations.md` limitation 13.7, which this contract inherits
  unchanged), and whether the per-spaxel IFU pattern of §9 is tolerable to write
  by hand at realistic N.
- **W1.13 (Spec assembly & freeze)** — four items, all in §19: the merge
  topology ruling (R1), cross-component hierarchical priors (R2), an optional
  `LikelihoodFamily.sample` (R3), and the channel-name relaxation (R4). Plus two
  carried forward: `substream`'s placement (`lowering.md` §12.7) and promoting
  the capability flags into W1.5's `Model`/`Transformation` ABCs so they are
  declared rather than duck-typed.
- **Phase 2 (backends)** — a backend supplies models and transformations that
  declare `DIFFERENTIABLE`, `BATCHABLE` and `DEVICE`, and nothing else: the
  whole of this contract is reused unchanged, which is the claim
  `DEVELOPMENT_PLAN.md` §3 makes for it. The reference implementations of
  `log_prob_unconstrained` and of `simulate`'s Gaussian draw are the oracles the
  native paths must agree with.
- **Phase 3 (SBI)** — `simulate` returns `Simulation`, whose `theta` is already
  the flat free-parameter vector an SBI package wants and whose `failed` flag is
  the reject-and-record signal §4.5 asks for. What is missing before Phase 3 is
  the batched form (limitation 17.5) and the coordinate–value–mask tensor
  encoding the plan's Phase 3 section describes, which belongs with the
  embedding networks rather than here.

## 19. Open questions for review

**R1 — the merge topology (§4). The main ruling this document asks for.**
Nested `ParameterMapping`, one merge per level, is implemented; §4.2–4.4 set out
the flat-leaf-merge alternative and compare them; §4.5 gives the rationale, §4.6
the one real cost and §4.8 the cost of reversing. The summary of the case: a flat merge
buys a tying advantage that turns out to be much narrower than it looks (an
outer tie already reaches an inner site — verified; only declaration-time
`shared_as` is lost), and costs a duplicated naming convention and names that
cannot be parsed unambiguously. Against that stands B's one genuine defect,
§4.6's: an inner tie is invisible to the top-level mapping, so provenance and
labelling consumers must reach for `sites()` rather than the obvious `bindings`.
The case is closer than §4.5 first made it look — one of its three reasons did
not survive being checked, and §4.5 now says so. **Reverse now or not at all** — after W1.8 emits these names into
stored `InferenceData` and W1.10 asserts them, it is no longer a cheap change.

A sub-ruling, and the one this document would most like granted: **adopt
"lossless nesting" in `parameters.md`** — let `ParameterSet.merge` accept a
`ParameterMapping` as a component and compose its bindings, so the outer
mapping's bindings are already the leaf bindings. It is strictly smaller than
the recursive merge `parameters.md` §14 contemplated, it removes §4.6's defect
at its source rather than papering it with a helper, and it is what W1.11's
population sketch independently recommends. If it is adopted, `sites()` and
`shared_names` become thin wrappers and this contract loses nothing.

**R2 — cross-component hierarchical priors (limitation 17.2).** A
`HierarchicalPrior` cannot reference a parameter in another merge component,
because `ParameterSet.__init__` requires references to resolve within the set.
The expressible pattern is "declare the hyperparameter locally, then tie" (§9),
which works and composes correctly, but it is not obvious, and a population fit
is the case it bites. Is the pattern acceptable as documented, or should
`parameters.md` gain deferred reference resolution at merge time? That would be
a §4.1 contract change, so it is Peter's, not this document's.

**R3 — should `LikelihoodFamily` gain an optional `sample`?** `simulate(observe=
True)` currently works for the Gaussian family and refuses everything else,
because a family declares only `log_prob` and this contract will not guess a
sampling distribution (§13). A one-method extension —
`sample(predicted, noise, rng)`, defaulting to a refusal — would make SBI work
with Poisson counts and Student-t outliers, which are exactly the families a
black-box simulator user has. It is a `likelihoods.md` §4.4 addition, so it
needs a decision-log entry; it is small and additive, and Phase 3 will want it.

**R4 — the channel-name relaxation (§14).** `results_schema.md` §17 routed the
"where does channel qualification happen" question here; §14 answers it (in the
model, never in the collection) and argues that the one-line relaxation of
`_check_channel_name` to accept dot-separated identifiers should land with the
first population model rather than now. Confirm, or ask for it at the freeze so
the accepted-name surface is settled once.

**Carried, and not this document's to rule.**

5. **`substream`'s placement** (`lowering.md` §12.7). It is in
   `ampere/core/rng.py`, which is what §9.2 asked for; W1.13 ratifies or moves it.
6. **Promote the capability flags into W1.5's ABCs.** `DIFFERENTIABLE`,
   `BATCHABLE` and `DEVICE` are read with `getattr` because
   `transformations.md` is frozen and this contract may not widen it (§10). They
   are declared against the `Capable` protocol, which is documentation until
   someone promotes them.
7. **Is the failure catch set right?** §11 catches `LikelihoodError` and the
   user's declared simulator exceptions, and nothing else (with `strict=True`
   catching nothing, per Peter's §17 Q1 ruling). The argument against
   catching more is that a bug turned into `-inf` is a fit that runs and is
   wrong. The argument for is that a user wrapping an untidy external code may
   not know every exception it raises, and will discover the list by having
   their run die at proposal 40 000. A `simulator_failures=(Exception,)` escape
   hatch already exists and is explicit, which seems the right compromise — but
   it is worth confirming that the default is the strict one.
8. **`Instrument` re-merges its steps on every call** (§11 of
   `transformations.md`, question 15.5). This contract caches its own mapping and
   each `Dataset`'s, but `Instrument.__call__` still calls `Instrument.mapping`
   per evaluation, which is a merge of a few dozen dictionary operations per
   dataset per `log_prob`. Negligible beside a GP Cholesky and measurable beside
   a cheap analytic model. W1.5's own answer — a `freeze()` that snapshots and
   refuses later mutation, never a silent cache — is the right one, and W1.7 is
   the consumer that confirms it is now worth doing.
