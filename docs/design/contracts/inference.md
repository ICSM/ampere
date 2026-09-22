# Ampere v2 — Dataset, FittingProblem & Inference Contract (W1.7)

Status: **frozen at `spec-v1.0`** (the tag created at the W1.13 merge,
2026-09; any later change to a §4 contract requires a decision-log entry in
`DEVELOPMENT_PLAN.md` in the same PR — ground rule 9). Previously: reviewed
and merged; **all four ruling requests R1–R4 were ruled by
Peter on 2026-09-02** (see §19's preamble for the dispositions — the merge
topology is ratified, `LikelihoodFamily.sample` and the dotted channel-name
surface landed the same day, and the tie-based hierarchical pattern stands as
documented). Two of Peter's earlier rulings — `likelihoods.md` §17 Q1 (a
strict toggle; the engine path records and returns −inf) and Q2 (the
`Dataset` resolves the effective mask once) — are implemented here, in §11
and §8 respectively. Implements
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
...     Capabilities, CauchyFamily, Dataset, DatasetCollection, Evaluation, FailureReason,
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
| `Capabilities`, `Capable` | `differentiable` / `batchable` / `device` / `backend`, and what a backend declares them on |
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

This is the contract's central design decision, and `parameters.md` §14 kept
it expressly open for W1.7 to settle rather than inherit. It is written out at
length because Peter ruled on it, not because it is complicated. **Ruled
2026-09-02: design B is ratified and stays** (§19 R1); the lossless-nesting
sub-proposal was granted later the same day and is landed —
`parameters.md` §8 carries the contract with executed examples, and §4.6's
defect is closed at source.

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

### 4.6 The strongest objection: nesting loses introspection — closed at source

This was the one real cost, found by W1.11's hierarchical-population sketch
rather than by this contract. **It no longer exists**: Peter granted the
lossless-nesting sub-ruling on 2026-09-02 ("lossless, not necessarily
recursive"), `parameters.md` §8 now lets `merge` take a `ParameterMapping` as
a component and compose its bindings, and this contract passes its mappings
accordingly — so the section below records what the defect *was*, and what
the surface says now.

**The loss, as it was.** A `ParameterMapping` used to record one binding per
component, not per leaf. So a parameter collapsed by an *inner* merge — two
steps of one instrument sharing a declaration-time `shared_as` label —
reached the top level as a single local name with a single binding, the
top-level mapping said it was not tied, and a consumer walking
`mapping.bindings` for provenance or ArviZ labelling would not learn that
the sampler dimension drives two places. Today the composed surface tells
the truth by itself:

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
('d.instrument.gain',)
>>> hidden.shared_names
('d.instrument.gain',)
>>> hidden.sites()["d.instrument.gain"]
('d.instrument.a.scale', 'd.instrument.b.scale')

```

`tied_names` now reports the collapse at whatever level it happened, because
the problem's bindings compose through every retained inner mapping.
`FittingProblem.sites()` remains as the convenient rendered form (merged name
to fully qualified leaf paths) and `shared_names` as the established name for
what is now the same statement `tied_names` makes; either the composed
bindings or `sites()` is right for W1.8's provenance — they carry the same
information. Value **routing** is untouched by all of this: `distribute`
stays one level deep, each composite re-distributes with its retained
mapping, and `Instrument.__call__`'s values path is exactly as
`transformations.md` §5 wrote it. The two views — leaf-level introspection,
level-by-level routing — are two readings of one retained structure, which
is what "lossless, not necessarily recursive" means in practice.

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

### 4.9 A third kind of top-level component: a joint noise group (*added W5.9*)

The merge topology above has two kinds of top-level component — one per model,
one per dataset — plus the optional `shared` set §9 uses. **W5.9 adds a third**,
and it is the first thing in this contract that is not one term per dataset.

A `JointGaussianProcessNoise` (`likelihoods.md` §7) is one correlated process
over `T` channels of one model on a shared grid, so it belongs to no single
dataset. It is declared on the `DatasetCollection`:

```python
DatasetCollection({"ra": …, "dec": …}, joint={"astrom": noise})
```

and joins the merge as one further component under its own label, exactly as
`shared` does and for the same reason: it owns parameters no dataset owns (the
kernel's, the coupling's, the group's `scale` and `jitter`). Its merged names
are therefore `astrom.angle`, `astrom.log_variance_0`, `astrom.length_scale`
and so on, and ties, plates, priors and W1.9's lowering reach them exactly as
they reach any other component's.

**The `"joint"` decomposition.** `DatasetCollection.contributions` returns one
log-likelihood per *key*, and a joint group replaces its members with **one**
key, its own:

```python
collection = DatasetCollection({"ra": ..., "dec": ...}, joint={"astrom": noise})
collection.contribution_labels()   # ('astrom',)  -- not ('ra', 'dec')
collection.group_of("ra")          # 'astrom'
```

The sum is unchanged — `log_likelihood` is still the sum of `contributions`'
values — and what changes is the *decomposition*. That widening is the same one
`"mixed"` was for `results.md` §6's pointwise group: a name for the case where
the obvious per-dataset reading does not apply, so that a consumer is told
rather than left to infer it. A run's per-dataset `log_likelihood` group
therefore carries a variable named for the **group** where a joint fit is
concerned, and `results.md` §6 says what the pointwise group holds (the rotated
outputs, one per member label, declared `"joint"`).

**`contribution_labels()` is the API**, not `tuple(collection)`: a caller that
assumed one term per dataset gets the right answer from it for a collection
with no groups and stays right for one with them.

**Simulation.** `DatasetCollection.draw_group` draws a group's `T` channels in
**one correlated call**, and `FittingProblem.simulate` routes grouped datasets
through it instead of through `Dataset.draw_observation`. This is not an
optimisation: drawing each channel from its own marginal produces observations
whose cross-covariance is zero, which is data from a different model than the
one being fitted — the precise failure an SBC study of a joint fit exists to
catch, and one that every marginal check passes.

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

**A third thing, since W5.9**: it owns the **joint noise groups**, passed as
`joint={label: noise}`. A group is checked at construction — the members exist,
no dataset belongs to two groups, the channels share a grid, a mask and a
per-sample uncertainty, and each member's own likelihood is a bare family — and
it is then one further component of the merge and one further entry of
`contributions()`, in place of its members'. See §4.9. The first bullet above
is unchanged in substance and sharper in wording: the joint log-likelihood is
still a **sum**, but a sum over `contribution_labels()` rather than over
datasets, because the datasets a group spans are not conditionally independent
given the parameters — which is exactly what the group says about them.

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

An instrument's label is **never** a merge component here, so two instruments
sharing a label on one channel costs nothing structurally — but it makes the
provenance record ambiguous, and **ruled 2026-09-03** (the
`transformations.md` §15 Q4 residual, landed at the freeze): the label and
the channel are distinct concepts — the channel says which part of the
simulation an instrument consumes, the label is how the user identifies the
instrument and which parameters are constrained by which data — so when more
than one instrument reads a channel, **distinct instrument labels are
required**, checked at problem composition (not at `negotiate`, which still
merges nothing). `Instrument.label` defaults to the channel name, so two
unnamed instruments on one channel are refused:

```pycon
>>> FittingProblem(
...     TwoChannel(grid),
...     DatasetCollection({
...         "plain": Dataset(observed, Instrument([], channel="blue", input_kind=Spectrum)),
...         "calibrated": Dataset(observed, Instrument([Calibrate()], channel="blue")),
...     }),
... )
Traceback (most recent call last):
    ...
ampere.core.exceptions.DatasetError: model 'model': 2 instruments read channel 'blue' but share the instrument label(s) ['blue']. ...

```

Named, the same composition works, and the `sources` tuple — the provenance
record the ruling protects — says which instrument asked for what. The
one-instrument default-to-channel-name case is unchanged:

```pycon
>>> shared_channel = FittingProblem(
...     TwoChannel(grid),
...     DatasetCollection({
...         "plain": Dataset(
...             observed,
...             Instrument([], channel="blue", input_kind=Spectrum, label="direct"),
...         ),
...         "calibrated": Dataset(
...             observed, Instrument([Calibrate()], channel="blue", label="scaled")
...         ),
...     }),
... )
>>> shared_channel.parameters.free_names
('model.index', 'model.norm', 'calibrated.instrument.calibrate.scale')
>>> shared_channel.requirements["model"]["blue"].sources
('direct', 'scaled')

```

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
...         Instrument([], channel="blue", input_kind=Spectrum, label=f"{label}_scope"),
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
`DatasetCollection.plate(...)` factory is the obvious convenience and was
deliberately not in v1.7 (limitation 17.6) — **it landed at W5.12**, below.
Design B's contribution here is the naming:
`spaxel_17.likelihood.amplitude` rather than a flattened
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

*`Binding.index` landed at the freeze (2026-09-02) and the declaration that
produces it landed at W5.12; the next subsection is the route this paragraph
was waiting for. The bound itself has not moved — it has become a refusal.*

### A population of datasets (*W5.12*)

`parameters.md` §9's `Population` is the composition-time declaration of
hierarchy, and it is what closes the paragraph above. A `FittingProblem`
takes them beside the ties they are the counterpart of —
`FittingProblem(models, datasets, populations=[...])` — and
`DatasetCollection.plate(name, datasets, members=..., hyperpriors=...)`
builds one from the datasets' own order, which is the factory limitation 17.6
named.

```pycon
>>> from ampere.core import HierarchicalPrior, Parameter, Population
>>> spaxels = [
...     Dataset(observed, Instrument([], channel="blue", input_kind=Spectrum,
...                                  label=f"scope{i}"),
...             model=f"spaxel_{i}", label=f"spaxel_{i}")
...     for i in range(3)
... ]
>>> plated = DatasetCollection.plate(
...     "spaxels", spaxels,
...     members=[Parameter("index",
...                        HierarchicalPrior("norm", {"loc": "mu", "scale": "sigma"}))],
...     hyperpriors=[Parameter("mu", st.norm(0.0, 1.0)),
...                  Parameter("sigma", st.halfnorm(0.0, 1.0))],
... )
>>> plated.populations[0].over
('spaxel_0', 'spaxel_1', 'spaxel_2')
>>> plated.populations[0].size
3

```

The draws are routed to the **models** the datasets name, not to the datasets
themselves, and that is a rule rather than a default: a population addresses
its members by *bare* local name, while a dataset joins the merge as a
`ParameterMapping` whose names are qualified (`likelihood.scale`), so an
element routed there would never reach a leaf. The merge refuses it by name.
A quantity genuinely shared *across* the spaxels' noise models is still a
`Tie`, as above — sharing and population are alternatives, not layers.

What this buys over the comprehension is the scaling the paragraph above
measured: the population's draws are **one** array-valued parameter, one
sample site, and one `numpyro.plate` when the problem is realised on torch or
jax (§10a). The N-component form remains available as
`Population(layout="flat")`, refused above `MAX_FLAT_MEMBERS` for the reason
this section has just given.

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

### Three nested samplers on one surface (*Added W5.14*)

The claim this section makes is that an engine consuming only §4.5's surface
runs on every backend, and `prior_transform` is the half of it that exists for
nested sampling. W5.14 is the first time that half has had **three**
consumers — `DynestyEngine` since W2.2, and `NautilusEngine` and
`UltranestEngine` from the inference-extensions memo's §6 tier 1 — so it is
the first time the claim has been checked by more than one library of the same
kind.

Nothing here changed to admit them, which is the finding. Both new drivers
consume `prior_transform` and `log_likelihood` and nothing else; both are
written against `ampere.inference.engine.Engine`'s shared machinery; and the
one place they differ from dynesty is in what their libraries *report*, not in
what they ask of a problem. The differences that did have to be absorbed are
all on the driver side and are recorded in `ampere/inference/_nested.py`:

* nautilus refuses a problem with fewer than two free parameters, so the
  driver refuses one by name rather than letting a library `ValueError` out of
  a constructor;
* nautilus reports no evidence uncertainty, so the driver estimates the
  importance-sampling one (`1/sqrt` of the library's own Kish effective sample
  size) and records that it did;
* ultranest draws from numpy's *process-global* generator, like zeus, so its
  run is wrapped in `ampere.inference.engine.global_seed` — the helper W5.14
  moved out of the zeus driver so that two engines share one implementation of
  seed-and-restore rather than two.

`results.md` §9's weighted-draw rule (W5.0) is obeyed by all three through
**one** function, `dynesty.utils.resample_equal` on the engine's own
`resample` stream: three nested samplers must not be able to produce three
slightly different posteriors from the same dead points. Each library's own
equal-weight output is deliberately unused, because each draws from its
library's randomness rather than from the problem's seed.

The engine-neutral evidence triple is written by all three, and
`ampere_evidence_method` is `"nested_sampling"` for all three — W5.0 carried
the question of whether that attribute should name the method family or the
engine, and the second and third evidence engines are the occasion to answer
it: the **family**, because `ampere_engine` already names the engine and what
a reader needs from the second attribute is whether two archived evidences
were estimated the same way.

### Four guide families and a second gradient library (*Added W5.14*)

The same claim, asked of the **realisation** surface (§10a) rather than of
§4.5's: an engine that consumes `log_prob_unconstrained` and nothing else runs
on any backend that registers one. W5.14's other two thirds are the first time
that half has been asked for something other than a sampler.

`VIEngine`'s `GUIDE_FAMILIES` gains `laplace` and `flow`, and nothing in the
contract moved to admit them — the driver already had the one-site model, the
unconstrained density and the emission shape. Two things had to be absorbed on
the driver side, and both are properties of the libraries rather than of this
surface:

* a Laplace autoguide is a `Delta` guide during optimisation in **both**
  libraries, so what SVI fits is the MAP location and the Gaussian exists only
  afterwards — pyro's `laplace_approximation()` returns a whole
  `AutoMultivariateNormal`, numpyro's `get_posterior(params)` builds the
  distribution on demand from the Hessian at that location. The driver draws
  from the Gaussian in each case, never from the `Delta`;
* a flow autoguide has no location parameter, so `init_loc_fn` — the hook the
  other three families are started with — does nothing for it. It is started
  by its **base distribution** instead (`get_base_dist`, which both libraries
  define and both `get_posterior` implementations call), placed at the same
  seeded prior draw every other driver starts from. In ampere's unconstrained
  coordinates that is not optional: a parameter declared with an unbounded
  prior *is* its own coordinate, so a flow left at the origin is starting many
  base standard deviations from the posterior, and the measurement that
  prompted this is in `ampere/inference/_vi.py`'s module docstring.

The contract consequence is one widened vocabulary, not a new key:
`results.md` §9's `ampere_approximation` gains `"laplace"` and
`"normalising_flow"` beside `"mean_field"` and `"multivariate"`. Every
consumer in the repository tests it against `"none"`, so the widening needs no
reader taught anything; what a reader gains is that the *family* is still
recoverable from the one engine-neutral key. Beside it, and new at W5.14,
`ampere_vi_guide_parameters` records which fitted parameters a run left on the
engine (`"loc, scale"`, `"loc, scale_tril"`, or `"none"` for a flow, whose fit
is a neural network's weights) — because "this run's guide has no location and
no scale" is a fact about an archived run, and an absent attribute is not a
statement.

`BlackjaxEngine` is the second consumer of §10a's density on the jax side, and
it is the first engine here that is **backend-specific by nature** rather than
by what happens to be installed: blackjax is a jax library and
`inference_extensions_memo.md` §2.2 records that torch has no counterpart to
borrow, so the driver refuses any other backend by name in its constructor
rather than reporting a missing library. `supported_backends()` is therefore a
one-element intersection rather than a table, which is the honest shape for a
driver with one route.

Its two methods make the same point from opposite ends of the approximation
question, and settle how `ampere_approximation` is to be read:

* **MCLMC writes `"none"`.** An unadjusted microcanonical chain carries a
  discretisation bias — there is no accept/reject step, and the tuner controls
  the bias by holding the energy variance near a target. `ampere_approximation`
  is nevertheless `"none"`, because the question that key answers is *"do
  chain diagnostics mean anything for this run?"* — it is what `plot_trace`
  and `ampere.results.summary` read before reporting an R-hat, an ESS or a
  trace shape — and an MCLMC run is a Markov chain, so they do. The bias is
  recorded as its own fact (`ampere_blackjax_adjusted`,
  `ampere_blackjax_desired_energy_var`) and is *measured* by the engine
  battery's SBC rather than asserted here.
* **Pathfinder writes `"pathfinder"`** and the per-draw
  `sample_stats.proposal_log_density` beside it, in the constrained
  coordinates W5.0 fixed, so the stored groups alone reweight it. blackjax
  returns the approximation's density beside its draws, so this costs the
  driver nothing but the coordinate change.

Neither writes the engine-neutral evidence triple, and that absence is
deliberate rather than pending: MCLMC samples an unnormalised density like
every other MCMC here, and Pathfinder's ELBO is a single L-BFGS path's *lower
bound* on the log evidence, whose tightness is unknown. A bound recorded under
the name a nested sampler's estimate uses would invite a comparison that has
no meaning, so the ELBO is kept under the engine's own name. The triple stays
what `results.md` §9 says it is: written by whichever engine *estimates* a
marginal likelihood.

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

`differentiable`, `batchable`, `device` and `backend` are properties of the
*pieces*: a problem is differentiable exactly when everything a gradient would
have to pass through is, and it is a torch problem exactly when every piece of
it is. On the reference path the honest answers are `False`, `False`, `"cpu"`
and `"reference"`:

```pycon
>>> problem.capabilities
Capabilities(differentiable=False, batchable=False, device='cpu', backend='reference')
>>> problem.differentiable, problem.batchable, problem.device, problem.backend
(False, False, 'cpu', 'reference')

```

**Promoted into W1.5's ABCs at the freeze** (ruled 2026-09-03, §19.6):
`Model` and `Transformation` carry `DIFFERENTIABLE`, `BATCHABLE` and
`DEVICE` as class attributes whose conservative defaults reproduce the
earlier `getattr` semantics exactly, so every piece a problem composes
declares them — silence inherits the reference answers — and
`declared_capabilities` reads the attributes directly. **`BACKEND` joined them
at W2.12** (decided by Fable 2026-09-07, `DEVELOPMENT_PLAN.md` §4.5 and its
decision log), for the same reason and with the same conservative default: the
base install's whole toolkit *is* the reference backend, and a hand-written
numpy model runs on the reference path.

```pycon
>>> Model.DIFFERENTIABLE, Model.BATCHABLE, Model.DEVICE, Model.BACKEND
(False, False, 'cpu', 'reference')
>>> from ampere.core import Transformation
>>> (Transformation.DIFFERENTIABLE, Transformation.BATCHABLE, Transformation.DEVICE)
(False, False, 'cpu')
>>> Transformation.BACKEND
'reference'

```

**One name per backend, everywhere.** The string a piece declares is *the*
name of that backend across the whole project: it is what
`declared_capabilities` aggregates onto `FittingProblem.backend`, what
`ampere.results` records as `ampere_backend`, the `backend` key
`lowering.md` §12.8's registry is consulted with, and the `name` of that
backend's fixture in `tests/conformance/`. There is no translation table
anywhere, and a backend that wanted two spellings would be introducing one.

Derivation of the two booleans is **conjunctive**, and an empty collection of
parts is *not* differentiable — `all([])` is `True`, and silently promising
gradients for a problem with nothing in it is the silent capability upgrade
this architecture forbids:

```pycon
>>> class Native:
...     DIFFERENTIABLE = True
...     BATCHABLE = True
...     DEVICE = "cpu"
...     BACKEND = "reference"
>>> class Conservative(Transformation):
...     def apply(self, samples, values):
...         return samples
>>> declared_capabilities([Native(), Native()])
Capabilities(differentiable=True, batchable=True, device='cpu', backend='reference')
>>> declared_capabilities([Native(), Conservative()])
Capabilities(differentiable=False, batchable=False, device='cpu', backend='reference')
>>> declared_capabilities([])
Capabilities(differentiable=False, batchable=False, device='cpu', backend='reference')

```

`device` and `backend` derive by a different rule, because neither is a
promise that can be weakened: all the parts must **agree**, and a disagreement
raises rather than being resolved. Ampere moves arrays between neither devices
nor array libraries on the user's behalf — the first turns a configuration
mistake into a silent performance collapse, the second into a run that fails
two steps later inside a backend, or that quietly has no gradients. The error
names the values it found and the `capabilities=Capabilities(backend=...)`
override for a caller who genuinely means one of them:

```pycon
>>> class NativeTorch:
...     DIFFERENTIABLE = True
...     BATCHABLE = True
...     DEVICE = "cpu"
...     BACKEND = "torch"
>>> declared_capabilities([NativeTorch(), Native()])
Traceback (most recent call last):
    ...
ampere.core.exceptions.DatasetError: the pieces of this problem declare different backends...

```

Phase 2's backends override these on their own `Model` and `Transformation`
subclasses; the `Capable` protocol remains the statement of the surface for
anything duck-typed into `capability_parts`.

### Non-native parts in a native problem

*(Added by W3.8, ruled by Peter 2026-09-08 on W2.4 slice 3's carried finding —
a post-freeze §4.5 addition; the decision-log row of the same date in
`DEVELOPMENT_PLAN.md` §2 is the record.)*

W2.13's fold-in 7 widened the capability parts to the noise model and the GP
solver, because a nominally differentiable problem whose GP solve ran in scipy
was a problem whose hyperparameters got no gradient. **The kernel joined them
here, for the same reason one level down**: a solver builds its covariance by
*calling* `Kernel.matrix`, so a numpy kernel inside a native solver hands back a
numpy array, the solver converts it, and the conversion detaches the graph in
exactly the amplitude and length scale a native GP fit exists to fit. The
argument that had kept the kernel out — "the solver is the piece that computes"
— does not survive contact with the code, and is withdrawn.

The ruling has three parts.

**Refused by default.** A part whose `BACKEND` is not the problem's is the
backend disagreement `declared_capabilities` already raises, and it names the
offending pieces by their **fully qualified** class — `ampere.core.likelihood.
Matern32` rather than `Matern32`, because a backend's own kernel deliberately
shares the bare name (`Likelihood.to_spec` records the declaration and
`results.md` §14 compares it across backends, so the two *must* agree there).
The remedy has two halves: build the piece from the backend's own classes, or
take the opt-in below and give up the gradient.

**Accepted under an explicit opt-in when no gradient is needed.**
`FittingProblem(..., allow_foreign_parts=True)` — spelt to sit beside `strict`,
and like `strict` an explicit per-problem declaration with no global form. It
exists for one case and it is not the case of a mistake: a piece expressible
**only in Python** — a tabulated kernel, a legacy callback, an external code
with no native twin — that a torch or jax problem wants to call through and
sample gradient-free. Under it the problem's `differentiable` is **`False`**
whatever the parts declare (a foreign part that claimed `True` would otherwise
buy back a gradient that does not exist), the run's provenance records
`ampere_foreign_parts` with the located names of the pieces, and the
gradient-free engines run on the numpy contract path.

Which backend is "the problem's" is resolved by the one asymmetry that is real:
`"reference"` is the name a piece inherits **by silence**, so it can never
identify a problem whose other pieces have deliberately declared a native one.
If the parts declare exactly one non-`"reference"` backend, that one is the
problem's and the `"reference"` pieces are foreign; two native backends have no
foreign/native split to find and are refused with or without the flag. The
check is **structural** — a part whose `BACKEND` is not the problem's — with no
list of substitutable pieces and no assumption about families, so a part stops
being foreign on the day it declares the backend and nothing here is edited
when core grows `sample` implementations or a backend grows a native twin.

**Always refused where a gradient is required.** `realise`, a backend's
differentiable `LoweredProblem`, `NUTSEngine` and `VIEngine` refuse by name
whatever the flag says, and whatever `strict` says: the flag declares that a
gradient-free run may call through a Python-only piece, and it cannot make one
differentiable. On the gradient-free side `realise`'s refusal *is* the routing —
`_EvaluationCache` treats a `LoweringError` as "no usable realisation here" and
scores on the contract path, exactly as it already does for an unlowerable
family — so §10a's fast path falls back and `ampere_realised` is `0`.

Whether a realisation could instead keep the native path and call the foreign
piece through the backend's escape hatch (torch: detach → numpy → tensor; jax:
`jax.pure_callback`) was **measured** rather than assumed, on the 400-point
flexible-likelihood problem of `tests/benchmarks/test_engine_fast_path.py`, one
`Engine.log_prob`: torch 10.7–12.6 ms contract, 7.8–9.1 ms callback, 7.3–7.5 ms
all-native; jax 10.8–13.6 ms contract, 9.3–9.8 ms callback, 2.0–2.2 ms
all-native. On torch the whole fast path is worth only 1.4–1.7×, and on jax —
where it is worth 5–7× — the callback recovers just 1.2–1.4× of it, because the
host round trip is precisely what jit exists to avoid. So the callback route
buys little where the fast path is cheap and almost nothing where it is
valuable, at the price of a second silently non-differentiable path. It is
**not implemented**: the fall-back stands.

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

## 10a. Realisation: the differentiable native form of a problem

*(Added by W2.13, decided 2026-09-07 — a post-freeze §4.5 addition; the
decision-log row in `DEVELOPMENT_PLAN.md` §2 is the record. Prototyped on
`w2.13-realisation-prototype` before the text was written, so every sentence
here has run.)*

**Why this section exists.** §10's surface is written against numpy and it is
**not traceable** on any backend, for three reasons that are each deliberate:
every `results_schema.md` container coerces its values with `numpy.asarray`
(a torch tensor that requires grad cannot enter a `Spectrum`; a jax tracer
raises `TracerArrayConversionError` at the same line); `ParameterSet.lnprior`
short-circuits on `math.isfinite`, which is Python control flow on a value;
and `Kernel.matrix` builds covariances in numpy. So a gradient-based engine
cannot be written against §10 alone, and `ampere.inference` — which may
import nothing but `ampere.core` and `ampere.results` — cannot reach a
backend's native evaluations by import. Both Phase 2 tracks found this
independently (W2.4, W2.5). `lowering.md` §0 already names the step that
closes the gap — lowering "happens once, when a `FittingProblem` is
*realised* on a backend" — and this section specifies it.

**Definition.** A **realisation** is a backend's one-way translation of a
composed `FittingProblem` into a native object that exposes the log density
as a pure function of the unconstrained free vector, built from the
backend's own model, step, noise, kernel and solver evaluations — never from
core containers in the hot loop. It is the *differentiable* form of the
problem; the numpy path of §10 remains the definition of the quantity.

**The registry** (`ampere.core.realisation`), with the shape of
`lowering.md` §12.8's:

- `register_realisation(backend, factory, *, override=False, builtin=False)`
  — one slot per backend, keyed on the backend's one name (the `BACKEND`
  capability flag); no silent overwrites; built-in rows distinguished.
  A backend registers at import: importing `ampere.backends.jax` is the
  user's opt-in to jax, and it is also the moment the jax realisation
  becomes reachable.
- `realise(problem) -> Realisation` — dispatches on `problem.backend` and
  **checks the result on use**: it names the same backend, its `free_size`
  is the problem's, and it agrees with `problem.log_prob_unconstrained` at
  the problem's reference values. A backend with nothing registered (the
  reference backend, which has no differentiable path) is refused by name
  with the remedy: import the backend that supplies one, or use a
  gradient-free engine.
- `registered_realisations()` for provenance and tests.

**The surface** (`Realisation`, a runtime-checkable Protocol). Mandatory:
`backend: str`, `free_size: int`, `log_prob_unconstrained(theta) -> native
scalar` (log prior + log likelihood + the change of variables, exactly §10's
`log_prob_unconstrained`, in the backend's array type). Optional:
`log_likelihood_terms(theta) -> mapping of dataset label to native scalar`,
consumed by a driver when present to emit the per-dataset decomposition
natively; absent, the driver recomputes the split for **stored draws only**
on the numpy path and records that it did (`engine_draws_recomputed`).

**Failure signalling on the native path** narrows §11 deliberately: inside a
trace nothing can raise and no `Failure` record can be built, so a realised
density returns a bare `-inf` (computed with the backend's `where`), and the
reasons §11 promises are recovered post hoc by re-evaluating stored draws on
the numpy path. `strict=True` on a realised problem means the **factory
refuses at construction**, by name, anything it cannot lower exactly (a
family, a noise model, a censoring declaration); runtime evaluation failures
are always `-inf`. This is the trace-purity ruling (`likelihoods.md` §17 Q1)
applied to the whole problem.

**Coverage.** W2.13 fixes the floor both backends must meet: Gaussian
families with independent or dense-GP noise, masks, plates and hierarchical
priors, with native kernels so GP hyperparameters are trainable; anything
else refused by name at construction. Widening (censoring, latent GPs,
non-Gaussian families, the quasiseparable solver) is each track's slice 2.

**What the conformance suite owes** (`tests/conformance`): for every
registered realisation, agreement with the numpy path at many points
including near a support boundary, at `tolerances.cross_backend`; the
one-point check in `realise` is a guard, not the proof.

**Provenance.** A run sampled through a realisation records
`ampere_realised = True` and the registered lowering rows it consulted
(`lowering.md` §12.8's stamping, populated automatically for the first time
— W2.6's deferral ends here, with `PROVENANCE_SCHEMA_VERSION` 5).

**The gradient-free engines are unchanged.** They consume §10 and run on
every backend; the reference backend registers no realisation in v1.

**W3.8**: `realise` refuses, before dispatch, a problem composing a part from
another backend — whether or not it was built with `allow_foreign_parts=True`,
and whether or not it is `strict`. A realisation *is* the differentiable form,
and that flag never buys a gradient. See §4.5, "Non-native parts in a native
problem".

**W4.3: the per-model native surface gained a second spelling.** Every native
model up to Phase 4 exposed its differentiable body and its compiled
coordinate grid as `flux`/`grid`, and the realisation machinery looked those
up by exactly those names. An interferometric source model cannot: every one
of them already declares a *parameter* called `flux` (the object's own
brightness), and `Parameterised._check_free_name`'s rule of the day — that a
parameter may not shadow a class attribute — made `flux` unavailable as a
method name on exactly these classes. `flux`/`grid` therefore stopped being
the only spelling: both backends' realisation machinery resolves either pair
on a model, and asks for both before refusing.

**W5.20: `native_flux`/`native_grid` is the canonical spelling.** *Ruled by
Peter 2026-09-15 (D1); decision-log entry in `DEVELOPMENT_PLAN.md` §2.* The
constraint that forced W4.3's second spelling is gone — `parameters.md` §10's
reserved set is core's, finite and stated, and no longer grows with a class's
own attributes, so an interferometric model could now have a method called
`flux` beside its parameter of that name. It still should not: a surface a
*realisation* composes and a quantity a *user* fits competing for one word is
a collision whether or not the namespace rule catches it. So the ruling keeps
the spelling and promotes it, rather than reverting to one name:

* **`native_flux`/`native_grid` is canonical.** Every shipped model on every
  backend offers that pair — the spectral models, the astrometry twins, the
  native astropy models, the interferometry twins — and both `problem.py`
  lookup tables list it first.
* **`flux`/`grid` remains a supported legacy alias.** A model written before
  the ruling composes, realises and differentiates unchanged; the
  documentation calls the pair *legacy*, and nothing shipped uses it.
* **A model offering both pairs is refused as ambiguous**, naming both
  methods. Before, the lookup took the first match and said nothing, so a
  model carrying a legacy `flux` beside a canonical `native_flux` — the shape
  a half-finished rename leaves behind — composed at whichever the table
  happened to list first. Two methods, one surface, no way to tell which the
  author meant: that is a question for the author.
* **A model offering half of either pair is refused with the missing half
  named**, as before. The two go together: the first supplies the values and
  the second the coordinates the instrument chain transforms them on.

**W5.12: a realisation routes its own values through the merge's table.** Both
backends' lowered problems call `ParameterMapping.distribute` on native
values — the wiring is the merge's, and reproducing it in each backend would
be two more places for the routing to be wrong. That makes `distribute` part
of the traced path in one specific respect: an **element binding**
(`Binding.index`, `parameters.md` §8) has to take its element *in the value's
own array type*. Coercing to numpy first, which is what it did until a
`Population` produced element bindings on a native path, raises on a torch
tensor that requires grad and on a jax tracer, and would detach the graph if
it did not. Nothing else in `distribute` touches values, so nothing else is
affected; this is recorded because it is the one line of `ampere.core` that a
backend's trace runs through by design rather than by accident.

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
 'non_finite_log_likelihood', 'execution_failed']

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

***Amended W5.10***, one sub-stream. A **context prior draws from the run's
own seed derivation and never from a generator of its own**: the contexts of
`simulate_many(count, context=prior)` come off `"<stream>.context"`, spawned
by index exactly as the draws are. Two properties follow, and both are the
reason it is a separate stream rather than the draw stream:

- a budget's θ and its noise are **bit-for-bit what they were without the
  argument**, so adding a context prior to a study does not silently re-draw
  the parameters it was comparing against;
- draw *i*'s context does not depend on the chunking, which is the same
  partition independence §13's batched form promises of everything else.

`rng=` names the *draw* stream, and therefore leaves the context stream alone:
a caller overriding the draws still gets contexts from the problem's own
derivation, because a prior that reached for `default_rng()` would make a
seeded budget irreproducible in the one place nobody would look.

**Placement.** `lowering.md` §9.2 puts `substream` "once in `ampere.core`"; its
§12.7 asked W1.13 to ratify that rather than let it be assumed. It lives in
`ampere/core/rng.py`, and **the home is ratified** (ruled 2026-09-03, §19.5):
pure stdlib+numpy, deliberately free-standing, the alternative having been
three backends agreeing by convention.

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

### What can be sampled, and what will not be guessed (*Amended W3.14*)

**Ruled by Peter, 2026-09-02 (R3)**: `LikelihoodFamily` has the generative
half — an overridable `sample(predicted, noise, rng)`, receiving the same
`NoiseParams` that `log_prob` scores with, with the **default a specific
refusal** that names the family and the override to provide. Observation
drawing delegates to it, so what a family cannot sample it refuses precisely,
and a user with an exotic observation process supplies it by subclassing:

```pycon
>>> class SamplingCauchy(CauchyFamily):
...     def sample(self, predicted, noise, rng):
...         return predicted + noise.sigma * rng.standard_cauchy(size=predicted.shape)

```

**Amended W3.14.** The refusal was the default for every shipped family but the
Gaussian one until Peter's ruling of 2026-09-09, whose use case arrived in W3.6:
a counting experiment could not `simulate(observe=True)`, so neither SBC nor SBI
on count data worked without the user subclassing `PoissonFamily`. The refusal
was written for families whose observation process is *genuinely ambiguous*, and
three of the shipped families are not — `poisson`, `student_t` and
`complex_gaussian` each have exactly one generative form, and each family's own
`log_prob` already fixes which one. Those three implement `sample`; `cauchy` and
every declared-but-unimplemented family keep the refusal **word for word**. What
changed is the list of families the contract can honestly say it knows the
answer for, not the rule: where the observation process is not determined by the
density, ampere still declines to guess it.

- **`poisson`** — `counts ~ Poisson(rate)`, returned as float because that is
  what the containers hold. `predicted` means what it means in `log_prob`: the
  expected count *before* the latent is applied. Under `IndependentNoise` that
  is the whole rate, and mean and variance are both `predicted`. Under
  `GaussianProcessNoise` the family's model is `counts ~ Poisson(predicted ·
  exp(f))`, so the rate is `predicted · exp(f)` and `f` is the latent the caller
  supplied — the *conditional* draw, at the same `f` `log_prob` scores at, which
  is what makes the whitened `z` in θ and the drawn data one model rather than
  two. `Dataset.draw_observation` reads that `z` out of θ and hands it to
  `noise_params`, exactly as `Dataset.log_likelihood_of` does on the scoring
  side. With no latent supplied a fresh GP realisation is drawn through
  `GPSolver.latent_transform`, as the Gaussian GP branch draws its own; the
  marginal is then a log-normal mixture of Poissons and is over-dispersed
  relative to a Poisson, which is what the model says. A negative or non-finite
  rate is refused by name (`log_prob` guards `rate <= 0`; the draw guards
  `rate < 0`, a zero rate being a well-defined point mass the density
  nonetheless declines to score).
- **`student_t`** — `x = μ + σ · t_ν`, with the family's own `ν` and the noise
  model's `σ` used **as the scale**, exactly as `log_prob` standardises by it.
  `σ` is not the standard deviation: the variance is `σ² ν / (ν − 2)` for
  `ν > 2` and undefined below, and matching the variance instead would be a draw
  from a different density from the one that scores it. A correlated (GP) noise
  model is refused by name — a Student-t is a scale mixture of Gaussians, the
  mixture does not commute with a GP covariance, and the marginal is not a
  Student-t at all (the composition is already refused at construction, the
  family declaring no `CONSUMES_LATENT_GP`).
- **`complex_gaussian`** — independent `Normal(0, σ²)` on the real and the
  imaginary part, which is the circular symmetry the family's density assumes.
  **`σ` is the per-component standard deviation, not the total**, and that is
  read off the density rather than assumed: `−|r|²/(2σ²) − log 2π − log σ²` is
  the joint density of two independent `Normal(0, σ²)` components, so
  `E|x − μ|² = 2σ²`. A correlated noise model is refused with the same message
  the density refuses it with — the circular complex GP is declared analytic and
  its implementation is Phase 4's, so there is no marginal to draw from either.

`GaussianFamily.sample` is implemented for both noise models — the
combinations the merged contracts get **provably right**:

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

A family that does not implement `sample` refuses at the point of use rather
than being flagged, because "this family cannot sample" is a fact about the
composition that would fail identically for every draw — and the refusal is
specific, naming the family and the override that provides the observation
process:

```pycon
>>> cauchy_problem = FittingProblem(
...     Rate(grid),
...     [Dataset(
...         counts.with_values([4.0, 7.0, 2.0], uncertainty=[1.0, 1.0, 1.0]),
...         likelihood=Likelihood(CauchyFamily(), IndependentNoise()),
...         label="counts",
...     )],
... )
>>> cauchy_problem.simulate({"model.rate": 3.0}, observe=True)
Traceback (most recent call last):
    ...
ampere.core.exceptions.DatasetError: dataset 'counts': the cauchy family does not implement
sample(): its log_prob defines how a datum is scored, not how one is generated, and ampere
will not guess a sampling distribution. Subclass CauchyFamily and override
sample(predicted, noise, rng) with the observation process, or use simulate(observe=False)
and draw observations from the predicted containers yourself.

```

The default refuses rather than guessing — inventing a symmetric error bar for a
density that fixes none would silently train an SBI posterior on the wrong
forward model. The noise-free half still works, which is what emulator training
wants:

```pycon
>>> cauchy_problem.simulate({"model.rate": 3.0}).failed
False

```

The three families W3.14 added draw instead, at the moments their own densities
imply:

```pycon
>>> latent_theta = {"model.rate": 3.0, "counts.latent.z": np.zeros(3)}
>>> drawn = latent_problem.simulate(latent_theta, observe=True)
>>> values = drawn.observations['counts'].values
>>> bool(np.all(values >= 0.0) and np.all(values == np.round(values)))
True

```

Masked samples keep the observed container's own values: they carry zero
information and are excluded from every likelihood, so drawing noise for them
would be inventing data.

#### Sampling on a backend (*Amended W3.1*, *W3.14*)

**Ruled by Peter, 2026-09-08**: every backend supports observation sampling
natively. W3.1 slice 2 landed it, and the ruling has a **ceiling** which is
stated here because it is the load-bearing half: *a backend samples exactly what
`ampere.core` samples.* A native twin for a family the core declines to guess a
sampling distribution for would be a backend inventing an observation process
the contract has just finished declining to guess — the "silently train an SBI
posterior on the wrong forward model" this whole subsection exists to prevent —
so `cauchy` keeps the refusal above **on every backend**, and a family whose
`sample` a *user* has overridden keeps it too: their override is the observation
process they wrote, and running something else instead would be worse than
running it slowly.

**Amended W3.14**, on both sides of that ceiling. The core now samples four
families, so the twins follow: `gaussian`, `poisson`, `student_t` and
`complex_gaussian` all draw natively on torch and on jax, through the
presence-based dispatch slice 2 left open. The ceiling itself is unchanged — a
backend samples what the core samples, no more — and the *floor* is stated too:
a family the core samples but a backend has no twin for **falls back to the
numpy path**, silently and correctly, rather than refusing. A backend may be
slower than the core; it may not be more restrictive than it, because a refusal
where the contract samples would be a capability withdrawn by an implementation
detail.

The two backends reach the same distributions by different routes, and the
difference is `vmap`'s, not the contract's. `jax.random` takes a traced rate and
a traced `ν` inside `jax.vmap`, so jax's twins are one expression each. Under
`torch.func.vmap` a random operation is either an error or driven by the global
generator, and the second would make a draw depend on how the chunk was
scheduled — the property the per-draw seed exists to protect — so torch's
non-Gaussian twins are **two-stage**: the transform computes each draw's
distribution *parameters* (the Poisson rate, including `exp(f)` for a latent GP;
the Student-t location, scale and `ν`), and the variate is drawn outside it, per
draw, from that draw's own seed. That keeps a draw a pure function of its seed
on both backends, which is what partition independence needs.

A realisation may therefore offer an optional
`sample_observations(theta, predicted, seeds)`, checked by presence exactly as
`log_likelihood_terms` is. It draws the retained values for a whole chunk in the
backend's own arithmetic and from the backend's own random stream; a dataset it
may not sample is refused **before** any draw is made, and `simulate_many` then
runs the numpy path, so what a user meets is `ampere.core`'s own refusal text
rather than a backend's paraphrase of it. Masked samples keep the observed
container's values by the same rule, applied in one place
(`Dataset.place_observation`) rather than reimplemented per backend, and a
censoring declaration that survives the mask blocks a draw on every backend.

**The numpy path stays the oracle, and the comparison is distributional.**
`jax.random` and `torch.Generator` do not reproduce numpy's stream and could not
be made to without reimplementing one library inside another, so a natively
drawn budget is a draw from the same distribution, not the same draw. The
conformance battery therefore holds a native draw to its *moments* — the mean
and the covariance of 4 000 draws against `K + diag(σ²)`, at the Monte-Carlo
standard error of the estimator itself — with the same two extra assertions the
numpy row carries and for the same reason: the joint tolerance alone would not
catch a dropped `K` (the off-diagonals collapse) or a dropped `σ` (the variances
become `K`'s alone), and both look entirely plausible in a plot.

Which stream produced a budget is **recorded, not inferred**:
`SimulationBatch.provenance['sample_backend']` names it, and the training-set
writer puts it in the file's root attributes as `ampere_sample_backend`. Its
value is `"reference"` whenever `LikelihoodFamily.sample` drew — including on a
torch or jax problem whose observations came from the loop.

### Batched form: `simulate_many` (*Amended W3.1*)

Limitation 17.5 said `simulate` is one draw and that `simulate_many(n)` was
Phase 3's. It landed at **W3.1 slice 1**, and the shape it landed in is fixed by
Peter's ruling of 2026-09-08 (`DEVELOPMENT_PLAN.md` §2, *Batched simulation*):
`vmap` is a single-device route and must not be the design, because a budget
routinely wants many devices or many machines and a single simulation may not
fit on one device at all. So **the loop is the reference semantics**, and it is
stated as an equality rather than as a description:

> `problem.simulate_many(n, stream=s)` returns, in order, exactly what
> `[problem.simulate(rng=child) for child in problem.rng(s).spawn(n)]` returns —
> whatever executor ran it, and however it was chunked.

That is the strongest equality available, and the reason it is spelled with
`spawn` rather than "n calls of `simulate`" is §12's advance-the-stream rule:
`simulate` *advances* its sub-stream, so draw 7 of a plain loop depends on how
many draws preceded it **in this process**, which a partitioned budget cannot
promise. Deriving one child generator per draw **by index** removes the
dependency by construction. numpy's spawning is cumulative — `spawn(a)` then
`spawn(b)` gives the children `spawn(a + b)` would — so children are spawned a
chunk at a time without the batch ever holding *n* generators, and `chunk_size`
1 and 7 give bit-identical batches.

θ is drawn **in the parent**, before any work is handed out, and the
already-advanced child travels with it. That is not an optimisation: it is what
lets a draw whose worker was killed still record the θ it was killed on, which
is the difference between reject-and-record and losing the evidence.

```pycon
>>> batch = seeded.simulate_many(5, observe=True)
>>> len(batch), batch.theta.shape
(5, (5, 1))
>>> batch.observations['default'].values.shape
(5, 3)
>>> bool(batch.failed.any())
False

```

A batch is a **sequence of `Simulation`s**, so every consumer of the single-draw
contract above keeps working unchanged, with stacked views over the columns an
SBI trainer wants laid beside it:

```pycon
>>> batch[2].observations['default'].values.tolist() == (
...     batch.observations['default'].values[2].tolist()
... )
True
>>> sorted(batch.results), len(batch.results['model'])
(['model'], 5)

```

The stacks are `ContainerBatch`, **not** `FunctionSamples` with a leading sample
axis. `results_schema.md` fixes a container's value shape to the one its
coordinates imply, so a container carrying a sample axis would be a container
whose values contradict its axes; `ContainerBatch` instead holds what is shared
once (kind, axes, unit, extra coordinates) and what varies as an array, and
indexing rebuilds the draw's own container. It is the same economy `results.md`
§11 gives the on-disk training set, for the same reason.

`values=None` draws each θ from the joint prior on the batch sub-stream — the
budget idiom. An array of shape `(count, free_size)`, or a sequence of that many
mappings or vectors, simulates at **given** θ: the simulation-based-calibration
idiom, where the θ are the ones whose rank statistics are being checked. One
mapping is refused rather than broadcast, because simulating a single point *n*
times is a deliberate thing to ask for and not something to arrive at by
accident.

```pycon
>>> given = seeded.simulate_many(3, values=[[1.0], [2.0], [3.0]])
>>> given.theta.ravel().tolist()
[1.0, 2.0, 3.0]
>>> seeded.simulate_many(3, values={"model.slope": 2.0})
Traceback (most recent call last):
    ...
ampere.core.exceptions.DatasetError: simulate_many's values= is one θ per draw ...

```

**The observation context (*Amended W5.10*).** `context=` was reserved by
`simulate_many`'s signature and refused any value but `None`; it now takes a
**`ContextPrior`** — `draw(rng, observed) -> ObservationContext` and
`describe() -> mapping`, structural like `Executor` — and draws **one context
per simulation** from it. `None` remains the default and means what it always
did: every observation drawn at the observed containers' own uncertainties.

An `ObservationContext` is a σ array per dataset label plus a small JSON-safe
`record` of what the prior drew. The σ reaches exactly one place — the
observed container of that dataset is replaced, **for that draw**, by one
carrying it (`Dataset.contextual_observed`) — so every family's `sample`
receives `NoiseParams` built from the context's σ and the drawn observation
carries it. That is what makes this a draw at a different noise *level* rather
than a draw at the observed level rescaled afterwards, and it is also why
`encoding.md` needs no new column: the packing already carries each sample's
`log σ` and whitened value, so a network reading `unpack` sees the context of
the observation it is conditioned on, and the reserved `context` column group
stays width 0.

Three shipped priors, and the list is the item's rather than a taxonomy:
`ScaledSigma` (scaled copies of the observed σ pattern), `SigmaArchive` (real
error arrays, drawn from uniformly) and `SignalToNoise` (a parametric S/N
model, which needs no observed σ at all). A user's own is any object with the
two methods.

The context is recorded **on every `Simulation`** (`Simulation.context`,
failures included — a draw that crashed in a context is evidence about that
context) and the *prior* in `SimulationBatch.provenance['simulation_context']`,
as canonical JSON of `describe()` or the string `"none"`. `results.md` §11's
training set stores the per-draw records in its optional `context` group and
the prior in `ampere_simulation_context`; the σ arrays are not stored twice,
because a drawn observation carries its own uncertainties.

Three refusals, each a claim the code could not honestly make:

- `observe=False` with a context — the context *is* the level the observations
  are drawn at, so a budget drawing none would record a context that did
  nothing;
- `native=True` with a context — the native path draws its noise from the
  realised problem's own σ, in the backend's arithmetic, which a per-draw
  container σ does not reach. With the default `native=None` the loop runs
  (which draws the context correctly) and the provenance records both facts;
- a problem declaring a **joint noise group** (W5.9) — a group's σ is read
  from the first of its datasets and its channels are drawn in one correlated
  call, so what a per-dataset context means for the cross-covariance is a
  design question rather than a detail, and guessing would put a training set
  on disk whose covariance and whose values came from two different contexts.

**Chunking.** `chunk_size` bounds how many simulations exist at once;
`as_chunks=True` yields `SimulationBatch` chunks lazily, and
`write_training_set`/`append_training_set` take the iterator, so a budget larger
than memory reaches the file without ever being held. The cost moves rather than
disappearing: each chunk after the first pays `results.md` limitation 13.9's
`O(existing + new)` append, which is quadratic in the *number* of chunks, so few
large chunks beat many small ones. `on_chunk` is called with each chunk's index
before that chunk runs — the device-placement hook. A single simulation larger
than one device is **not** partitioned by ampere; a model-parallel simulator is a
`Model` whose `__call__` does its own placement, and `chunk_size=1` plus the hook
is what guarantees such a model is never asked to hold two simulations at once.
W3.1 slice 2 wires the backends into the hook and adds per-chunk `vmap`
*underneath* this contract, never instead of it.

**`BATCHABLE` on the reference backend.** The flag means one thing said in two
dialects — *this part can take a stack of θ*. On torch and jax that is
`log_prob_unconstrained_batched` through `vmap`, which a backend model declares
without implementing anything here; on the reference backend it is
`Model.evaluate_batch(batch)`, the form an external code that takes a *table* of
parameter sets in one call already has. `simulate_many` uses it under the serial
executor only, and only when the model both declares the flag **and** overrides
the hook — a table evaluated in one call is by definition not partitioned across
workers, so the two are alternative ways of spending the same batch, and a torch
model that declares the flag for `vmap` falls to the loop honestly.

**The native path, and exactly what it guarantees (*Amended W3.1*, slice 2).**
A realisation may offer an optional
`simulate_batched(theta, *, chunk_size=None, sharder=None)` — checked by
presence, as `log_likelihood_terms` is — which runs a chunk's whole noise-free
forward model through `torch.func.vmap` or `jax.vmap` and hands back a
`BatchedPrediction`: every model's every channel, and every dataset's
instrument-transformed prediction on the whole observed grid. *Every* channel,
because `results.md` §11 writes one training-set group per `<model>.<channel>`
and a fast path returning fewer would silently write a smaller file than the
loop. θ is the **constrained** free vector — the coordinates
`SimulationBatch.theta` holds, not the unconstrained ones a sampler works in,
because a simulation budget is a set of parameter values rather than a set of
sampler positions.

`simulate_many` uses it when the problem is realised and every part declares
`BATCHABLE`, and the guarantee it then makes is deliberately **two-graded**,
because one sentence could not be true of both halves:

- **the noise-free prediction** still satisfies the equality above, to
  `tolerances.cross_backend` rather than bitwise. Vectorised arithmetic is not
  scalar arithmetic in the last digits, and a `vmap` is free to accumulate in a
  different order; the claim is that it is the same function, not that XLA and
  ATen reassociate identically. The *partition* independence is unweakened, and
  is a stronger claim here than on the loop, because a chunk is what gets
  `vmap`ped and therefore decides the shape of every array in the trace:
  `chunk_size` 1, 7 and the whole budget give the same predictions;
- **the observations**, where they were drawn natively, are a draw from the same
  distribution and *not* the same draw — see "sampling on a backend" above.

`native=False` switches the fast path off and restores the bitwise equality
exactly, which is the way that claim stays checkable; `native=True` requires the
fast path and **refuses by name** when a part is not `BATCHABLE`, when an
executor was given (a pool partitions the draws and a `vmap` evaluates them
together — they are alternative ways of spending one chunk), when the backend
registers no realisation, or — ***W5.10*** — when a `context=` was given as
well. The default, `native=None`, uses it where it is
available and falls back to the loop where it is not, recording which happened
in `SimulationBatch.provenance['simulate_batched']` — written to a training set
as `ampere_simulate_batched`, and merged conservatively across chunks so a mixed
budget never reads as a native one. `Model.evaluate_batch` is recorded
separately (`ampere_evaluate_batch`): a model taking a table of θ and a lowered
problem running through a backend's `vmap` are different claims, and one
attribute could not answer for both.

Three properties are contractual and follow from the design rather than from
care. The `vmap` is **per chunk with the chunks looped**, never over the budget
— that is exactly the single-device memory trap the ruling names. The
`Simulation` objects a native chunk produces are the loop's: containers rebuilt
from templates taken once at the reference values, which is sound because a
predicted container's axes and its effective mask are evaluation-invariant *by
contract* (`Likelihood.check_alignment` and `Dataset._masked_pair`). And the
native path is checked against the contract path at that same reference point,
channel by channel and dataset by dataset, to the tolerances `realise` uses for
the density — a guard rather than a proof, for the same reason §10a gives, with
the conformance battery making the full comparison.

**Sharding a chunk (*Amended W3.1*, slice 2; Peter, 2026-09-09).** There are
three axes a budget can be spread along and the contract keeps them apart:
across processes or machines is the `Executor`; across the draws in one chunk is
the `vmap`; across the accelerators of one host is a
`ChunkSharder` — `devices()` and `shard(fn, stacked)`, where the contract is a
*value* contract: sharding must not change the answer, only where it was
computed. jax ships `SingleDeviceSharder` and a `pmap`-based `MeshSharder`;
torch ships `SingleDeviceSharder` and a `DistributedSharder` that evaluates one
rank's stride of the chunk and gathers. ampere **does not launch a process
group**: a multi-rank run is launched under `torchrun`, and a sharder built
without one refuses by name rather than hanging in a collective. The
single-device case is the reference implementation and is what CPU CI exercises;
the multi-device rows live in `tests/gpu` and skip without hardware, with the
full exercise scheduled with the GPU item.

### Execution: the executor protocol (*Amended W3.1*)

`ampere.core.simulate.Executor` is one method — `map(fn, items)`, results in
item order — and it is deliberately the shape `concurrent.futures.Executor`
already has, so dask's `Client.get_executor()`, ray's executor wrappers and
mpi4py's `MPIPoolExecutor` satisfy it as they stand. A queue-driven cluster
array is the one case that wants a thin adapter, and the adapter is that method.
Two properties are contractual, and both are properties of the *result* rather
than of the schedule: results come back positionally, and how the executor
groups, distributes or retries items must not change the values it returns.
`simulate_many` guarantees the randomness half of the second by handing each
draw its own generator, so an executor only has to avoid reordering.

Three are shipped. `SerialExecutor` is the default and the oracle.
`ProcessExecutor` is the route for an external simulator: one worker process per
simulation in flight, the problem broadcast to each worker once at start-up
rather than pickled per draw, and `max_tasks_per_child=1` available for a
routine that cannot be run twice in one process. `ThreadExecutor` suits an
I/O-bound wrapper.

**The pool outlives a `map` call (*Amended W3.1*, slice 2).** Slice 1 built one
per call, which made `chunk_size` — whose job is bounding *memory* — bound
throughput as well, since every chunk paid for process creation plus one pickle
of the problem per worker. `ProcessExecutor` now keeps its pool, and the
replacement rule is exactly the two states in which the held pool is unusable:
**a pool is replaced only when a worker died or a draw expired.** Broadcasting a
different payload closes it too, since workers are handed the payload at
start-up. `shutdown()` therefore matters, and the context-manager form calls it.

**The default start method is `forkserver` on POSIX (*Amended W3.1*, slice 2;
ruled for that item).** Not the platform's own, for two reasons that point the
same way: forking a process that already holds a threaded runtime may deadlock —
jax warns about it in as many words, and torch's intra-op pools have the same
shape — and Python 3.14 moves the platform default off `fork` on Linux for that
reason, so inheriting the platform default would mean ampere's behaviour
changing at an interpreter upgrade rather than at a decision. `forkserver` over
`spawn` because it keeps most of the start-up saving; `fork` remains available
through `mp_context=`. The price is `spawn`'s: everything a worker touches must
pickle **and unpickle**, which is why `simulate_many`'s picklability check is a
round trip rather than a dump, and why a `timeout` is now taken only after the
workers have been started with an undeadlined no-op — a deadline measured from
submission is a deadline on the *simulation* only if the worker is already
running when the item is submitted.

**What is flagged and what is raised.** A per-simulation `timeout` on either
pool expires as a flagged failure, never an exception, and so does a worker that
dies outright — `FailureReason.EXECUTION_FAILED`, its own code beside
`MODEL_FAILED` because "the simulator said no" and "the executor lost the
simulator" call for different remedies and a single count could not tell a user
which they had. An exception the *simulator* raises propagates, exactly as it
would from `simulate`: a declared `simulator_failures` class is already flagged
inside `simulate`, so an exception that gets that far is an undeclared one, and
turning it into a silently dropped draw would hide a bug inside a plausible
failure rate.

Two limits are stated rather than engineered around. A **thread** cannot be
interrupted, so `ThreadExecutor`'s deadline bounds the answer and not the work:
the draw is flagged at the deadline, and `map` still returns only once the
abandoned thread has finished, because leaving a non-daemon worker thread
running would hang the interpreter at exit instead. A **crashed worker** breaks
its whole pool and the operating system does not say which draw did it, so the
draws that were merely in flight beside it are re-run one per fresh pool, where
whichever draw is guilty convicts itself rather than an innocent neighbour.

**Failure accounting is the parent's.** Recording is suspended while a chunk
runs and replayed on the problem in draw order, whichever process produced the
draws — so `failure_counts` is right after a pooled budget, which limitation
17.7 says it would not otherwise be. That limitation stands for the engines that
drive their own pools; this path is the exception, and it is one because
`simulate_many` owns both ends of the pool.

**Picklability is a precondition, and it is checked.** A process pool sends the
problem to its workers, so every model, instrument, likelihood and prior in it
must pickle; `simulate_many` checks once, before a worker starts, and refuses
with a message naming the usual culprits rather than letting an opaque
`PicklingError` surface from inside a worker's bootstrap. Making that true at
all needed a `copyreg` reduction for `mappingproxy`
(`ampere/core/_pickling.py`), since every frozen mapping in `ampere.core` is one
and none of them could be pickled before. Since W3.1 slice 2 the check is a
**round trip**: a forked worker inherits the problem by memory and never
reconstructs it, whereas a `forkserver` worker does, and a class that writes but
does not read is a failure inside a worker bootstrap rather than a sentence here.

**The process pool is therefore not a universal executor, and that is the right
answer rather than a gap.** A problem composed on the *reference* backend
pickles, which is what matters: it is where a wrapped external simulator is
composed, and an external simulator is the case the pool exists for. A problem
composed on **jax** does not — a jax array carries a `jaxlib` `Device` handle,
which is process-local and has no pickle reduction, and jax warns in its own
right that forking a jax process is likely to deadlock — so `simulate_many`
refuses the pool for it, by name, before any worker starts. Throughput on a
device backend is not more processes; it is W3.1 slice 2's per-chunk `vmap`,
with the serial and thread executors and `chunk_size` still available meanwhile.
The conformance battery carries this as a declared capability
(`BackendCapabilities.picklable`) and asserts *both* halves: a backend that
claims to pickle round-trips, and one that does not must genuinely fail to,
so the refusal can never rest on a stale declaration.

### Truncated proposals over this surface (*Amended W3.4*)

W3.4 lands `SBIEngine(method="tmnre", ...)` — truncated marginal ratio
estimation, Miller et al. (2021) — and it asks this section for exactly one
thing beyond what W3.1 already fixed: **`values=` is how a round proposes.**
The other multi-round methods propose from the *posterior* the previous round
trained, which is an `sbi` object that can be handed to `append_simulations` as
a `proposal=`. TMNRE's proposal is a truncated **prior** — the joint prior
renormalised on the hyperrectangle where the estimated 1-D marginals put their
mass — which is not an `sbi` posterior at all, so the engine draws from it
itself and passes the resulting θ table through `values=`. Nothing here
changes: `values=` already means "simulate at *these* θ", one per draw, and
the per-round sub-stream rule (`"sbi.simulate.<round>"`) is what keeps a round
reproducible whatever the previous round's rejection sampler consumed.

Two consequences worth recording where a reader of this section will meet them.

**The truncated prior needs no importance correction, and that is a property of
`values=` rather than a claim about the estimator.** Because the proposal is
the prior restricted to a subset — not a learned distribution — the target
inside the box is the same target, so the pairs a round produces are drawn from
the same joint the untruncated budget draws from, conditioned on the box. A
proposal that were a *fitted* posterior would not have that property, which is
why `rounds > 1` costs amortisation for every method and an importance weight
for some of them.

**A truncated round's cost is a prior-sampling cost, and it is ampere's.** The
box is enforced by rejection against the joint prior, and `sample_prior` is a
Python loop because §6's ties and §9's hierarchies are resolved per draw. A box
holding a thousandth of the prior's mass therefore needs of order a thousand
`sample_prior` calls per simulation, whatever the simulator costs. The engine
records each round's acceptance rate for exactly this reason and warns below
1e-3; a future vectorised `sample_prior` for the tie-free case would remove
the constant, and nothing in this section forecloses it.

### Training and sampling are reproducible too (*Amended W3.15*)

The per-round sub-stream rule above fixes every *simulation*, but until W3.15
it fixed nothing else: neither `ampere` nor `sbi` seeded torch's own global
generator, so a density estimator's initial weights and a trainer's batch
order came from whatever state torch's process-wide RNG happened to be in.
Two runs of the same seeded problem therefore agreed on every simulated pair
and disagreed on the posterior trained from them, so W3.5's "identical
posterior" claim held only through the artefact cache, not from the
problem's seed itself.

`SBIEngine.run` now seeds torch's global generator from the problem's own
sub-stream — `self.integer_seed("sbi.torch")`, the same idiom `_nuts.py` and
`_vi.py` use for pyro's kernels — immediately before every step that spends
it: building the first network, each round's `train()` (TMNRE's marginal and
joint estimators included), a multi-round proposal's own `sample()` call, and
the run's final posterior draw (a distinct sub-stream, `"sbi.torch.sample"`).
The generator behind a label is created once and then advanced (§ above), so
calling `integer_seed` with the same concern before every round's training
already gives each round its own seed, with no per-round label arithmetic
needed. Reseeding immediately before the final draw — rather than trusting
whatever state training left behind — is what makes a **cache hit** sample
reproducibly too: a hit trains nothing, so without its own seed point the
draw would depend on how much randomness the (still unconditionally
constructed, then discarded) network setup happened to consume. The seed
actually used to start training is recorded as `ampere_sbi_torch_seed`;
absent, not `None`, when `problem.seed is None`, because an unseeded run
asked for fresh randomness and torch is not pinned behind its back.

**A second global generator turned up beside torch's.** `sbi`'s default MCMC
method (`"slice_np_vectorized"`, what an NLE/NRE posterior and a TMNRE
`sample_with="mcmc"` one both sample by) draws its slice proposals through
`numpy`'s *legacy* global generator (`np.random`) directly, not through
anything `torch.manual_seed` reaches and not through `problem.rng`'s own
sub-streams either — found by this item's own TMNRE acceptance check
refusing to repeat once torch alone was seeded. So the same call that seeds
torch also seeds `np.random`, wherever a posterior might be MCMC-sampled; an
NPE posterior's flow never reads it, so doing this unconditionally costs
that path nothing. And both are **saved and restored** around each seeded
block — this driver's own version of the rule `_nuts.py`'s
`torch.random.fork_rng` and `_zeus.py`'s `_global_seed` both state already:
ampere does not leave a library's global random state changed behind it, so
unrelated code drawing from `np.random` or torch after `run()` returns is
exactly as reproducible as it would have been had this engine never run.

Two things this does not reach, deliberately. TMNRE's own round-to-round
proposal (the truncated prior, `_propose`) draws through the same
numpy-backed `_UnconstrainedPrior` every other prior draw in this section
does, so it needs no seed of its own at all — it was already reproducible.
And `calibrate()` never retrains — a rejection-sampled TMNRE posterior is
merely *rebuilt* under MCMC (§ above), not refit — so it has no training step
for this to seed; its own SBC/TARP posterior sampling (`sbi.diagnostics`,
internal to that call) is unaffected and remains as reproducible, or not, as
it was before W3.15.

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

**The relaxation has landed — ruled by Peter, 2026-09-02 (R4).** The name
surface is settled now rather than with the first population model, on his
observation that grouped data that are not hierarchical at all want it too —
one object's several sub-mm CO lines (`co.j3_2`, `co.j2_1`) are the worked
case. `_check_channel_name` accepts `.`-separated identifiers; an
`Instrument`'s channel binding follows the same rule (and so does its label,
which is provenance, never a merge component). What stays bare-identifier is
every **merge component** — step labels, dataset labels — so a dotted
instrument label reaching a `Dataset`'s default is refused loudly and the
user names the dataset explicitly.

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
| Capability flags are class attributes on W1.5's ABCs, defaulting to the reference answers | Promoted at the freeze (ruled 2026-09-03, §19.6), replacing the interim `getattr` reads with identical semantics: every composed piece declares the three flags, silence inherits `False`/`False`/`"cpu"`, and `declared_capabilities` reads them directly |
| `simulate` draws Gaussian observations and refuses everything else | A family declares only `log_prob`; guessing would train SBI on the wrong forward model |
| `substream` lives in its own module | `lowering.md` §12.7 asked W1.13 to ratify or move it; **ratified in place** (ruled 2026-09-03) — pure stdlib+numpy, deliberately free-standing |

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
   (§9). **Ruled by Peter, 2026-09-02 (R2): the tie-based pattern stands as
   the documented route.** For the revisit, if the pattern ever proves too
   awkward in practice, the required change is recorded here so it does not
   have to be re-derived — it is `parameters.md`'s, in three parts:
   (a) `ParameterSet.__init__` accepts a hierarchical reference it cannot
   resolve when (and only when) an explicit flag marks it deferred
   (`HierarchicalPrior(..., defer=True)` or equivalent), so today's loud
   refusal stays the default and a typo'd reference is still caught;
   (b) `ParameterSet.merge` resolves deferred references against the *merged*
   namespace — the reference is then a full merged path (`shared.mu`),
   resolved after qualification and tie collapse in the same pass that
   already rewrites local references, and refused loudly if it still
   dangles; (c) the lowering consequence is confined to naming — a deferred
   reference lowers exactly like a local one once resolved, since the
   topological ordering is over merged names already — so no backend work
   follows. The cost that kept this out of v1.7: error timing moves from
   declaration to merge for deferred references, a real loss of locality,
   and nothing yet needs it — every population case composes with the tie
   pattern.
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

   *(Amended W5.4, 2026-09-16.)* The limitation is unchanged; what changes is
   a claim this item used to make in passing, that the latent size is one
   whitened value **per retained sample**. That was a property of the two
   exact solvers, not of the contract. The whitened block is whatever the
   solver's whitening takes: a reduced-rank strategy factorises `K` as
   `(N, m)` rather than `(N, N)`, so `m` whitened variables produce `N`
   correlated ones and the block is `m`. `GPSolver.latent_size(kernel,
   n_samples)` is the single place that answers it — `n_samples` by default,
   which is what `DenseGP` and `QuasisepGP` return and why nothing about the
   two exact paths changes, and the basis size for `HilbertSpaceGP`
   (`likelihoods.md` §7). It settles `docs/design/horizon_notes.md` §2's
   question (b).

   The invariance the item is *about* is untouched and is what makes the
   ruling safe: the size must be computable from the **declaration alone**,
   before any container is in hand, because `Likelihood.latent_declaration`
   is called at composition. `HilbertSpaceGP`'s `basis_size` is therefore a
   per-axis declaration whose product is `m`, checked against the axes the
   kernel selects at composition rather than discovered as a shape error
   inside a solve. The obligation the horizon note attached to the ruling —
   that `simulate(observe=True)` and the latent-GP likelihood path must agree
   on which whitening they use — is met structurally: every family's `sample`
   draws `solver.latent_size(...)` whitened values and applies
   `solver.latent_transform` to them, which is the same transform
   `GaussianProcessNoise.noise_params` applies to the declared latent block,
   and a conformance row asserts that the covariance the marginal likelihood
   scores is exactly the covariance that whitening draws from.
5. **`simulate` is one draw.** A batched `simulate_many(n)` — which is what an
   SBI budget actually wants, and what a `batchable` backend could vectorise —
   is Phase 3's, and needs the capability flag to mean something first.
   *(Amended W3.13, 2026-09-10)*: **closed at W3.1.** `simulate_many` lands
   above (the "Batched form" and "Execution" subsections of §13), on an
   order-preserving, partition-independent executor protocol; slice 2 adds
   the per-chunk native form (`simulate_batched`) that a `BATCHABLE` backend
   vectorises. Kept in this list, like item 9, because the limitation it
   replaces was load-bearing in earlier discussion.
6. **No `DatasetCollection.plate(...)`, and the comprehension does not scale.**
   ~~N per-spaxel datasets are built by a comprehension the user writes (§9). A
   factory is convenience, not contract.~~ What is *not* convenience is the
   cost: N scalar components cost about 800 ms per `lnprior` at N = 1000
   against about 0.6 ms as a `Plate` (`hierarchical_population.md` §7), so the
   pattern is right for tens of datasets and wrong for thousands. Closing that
   needs per-element routing of a plate's array-valued parameter — an optional
   `Binding.index` — which is `parameters.md`'s to add.

   **Lifted at W5.12.** `Binding.index` landed at the freeze (2026-09-02) and
   `parameters.md` §9's `Population` is the declaration that produces it; the
   factory is `DatasetCollection.plate(...)` (§9). The caller still writes the
   comprehension that builds the N `Dataset`s, because each carries its own
   observations and nothing can guess them — what the factory removes is the
   *wiring*, which was N hand-written `PlateBinding`s or N rewritten parameter
   sets. The cost sentence survives as a refusal: `Population(layout="flat")`
   is the comprehension's own shape and is refused above `MAX_FLAT_MEMBERS`.
7. **The failure history is per-process.** Under multiprocessing (emcee's
   `Pool`), each worker accumulates its own counts and the driver must aggregate
   them. W1.8 owns the aggregation when it writes provenance.
8. **Capability flags are read, never verified.** A model declaring
   `DIFFERENTIABLE = True` is believed. There is no way to check the claim from
   `ampere.core`, which has no autodiff; the conformance suite (W1.10) is where a
   backend's claim gets tested.
9. **Resolved — an inner tie is now visible to the top-level mapping**
   (§4.6). Lossless nesting landed on 2026-09-02: the problem passes its
   retained mappings into the merge, so `tied_names` and the composed
   `bindings` report every collapse at whatever level it happened;
   `shared_names` and `sites()` remain as the established rendered forms of
   the same information. Kept in this list because the limitation it
   replaces was load-bearing in §4's ruling discussion.

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
  declared rather than duck-typed. *(All closed: R1–R4 were ruled 2026-09-02
  and landed then; the carried pair were ruled 2026-09-03 and landed at the
  freeze — §19's preamble is the record.)*
- **Phase 2 (backends)** — a backend supplies models and transformations that
  declare the **four** flags `DIFFERENTIABLE`, `BATCHABLE`, `DEVICE` and
  `BACKEND` *(the fourth added W2.12)*, and nothing else: the whole of this
  contract is reused unchanged, which is the claim `DEVELOPMENT_PLAN.md` §3
  makes for it. `BACKEND` is the backend's one name everywhere — the key
  `lowering.md` §12.8's registry uses and the `name` of its conformance
  fixture — and it is what makes a run's `ampere_backend` a fact rather than
  a declaration by whoever built the engine. A backend that offers gradients
  also registers its **realisation** (§10a) at import; that is the whole of
  what a gradient-based engine needs from it. The reference implementations of
  `log_prob_unconstrained` and of `simulate`'s Gaussian draw are the oracles the
  native paths must agree with.
- **Phase 3 (SBI)** — `simulate` returns `Simulation`, whose `theta` is already
  the flat free-parameter vector an SBI package wants and whose `failed` flag is
  the reject-and-record signal §4.5 asks for. What is missing before Phase 3 is
  the batched form (limitation 17.5) and the coordinate–value–mask tensor
  encoding the plan's Phase 3 section describes, which belongs with the
  embedding networks rather than here. *(Amended W3.13, 2026-09-10)*: **both
  have landed.** The batched form is §13's own "Batched form" and "Execution"
  subsections (*Amended W3.1*) and limitation 17.5 above is closed; the
  encoding is `docs/design/contracts/encoding.md`, frozen at W3.3 as its own
  contract exactly where this sentence said it belonged, and consumed by
  `SBIEngine(layout=...)`.

## 19. Open questions for review

**Ruled by Peter, 2026-09-02 — all four requests.** R1: design B (nested
`ParameterMapping`, one merge per level) is **ratified**; later the same day
the **lossless-nesting sub-ruling was granted too** ("lossless, not
necessarily recursive"), together with `Binding.index` — both landed in
`parameters.md` §8 the same day; `sites()` and `shared_names` are now thin
wrappers and §4.6's defect is closed at source. R2: the tie-based pattern for
cross-component hierarchical structure **stands as documented**; the design
for a future revisit (deferred reference resolution at merge time) is
recorded in limitation 17.2 so it need not be re-derived. R3: **granted, and
landed the same day** — `LikelihoodFamily.sample(predicted, noise, rng)`
exists with a default that refuses specifically, `GaussianFamily` implements
it for both noise models, and a user family overrides it to supply an exotic
observation process (§13; decision-log entry in `DEVELOPMENT_PLAN.md` §2).
R4: **the name surface is settled now, not at Phase 5** — grouped
non-hierarchical data (one object's several sub-mm CO lines) want dotted
channels too, so `_check_channel_name` and the `Instrument` channel binding
accept `.`-separated identifiers (§14; merge-component labels stay bare).
The original requests are kept below for the record. Item 8 was closed on
2026-09-03 by Peter's acceptance of `transformations.md` §15 Q5 —
`freeze()` is confirmed as worth doing, and W1.13 lands it. **Later the
same day items 5–7 were ruled too**: item 5 — `substream` is ratified in
`ampere.core` (closing `lowering.md` §12.7 with it); item 6 — the
capability flags are promoted into W1.5's ABCs at the freeze, as class
attributes whose conservative defaults (`False`/`False`/`"cpu"`)
reproduce the current `getattr` semantics exactly *(landed — §10's
"Capability flags" carries the promoted form, and
`declared_capabilities` reads the attributes directly)*; item 7 — confirmed:
the catch set stays narrow by default (`LikelihoodError` plus declared
`simulator_failures`, with `strict=True` catching nothing), and
`simulator_failures=(Exception,)` remains the explicit escape hatch for
an untidy external simulator.

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
