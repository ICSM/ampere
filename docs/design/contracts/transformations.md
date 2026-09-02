# Ampere v2 — Transformation & Instrument Contract (W1.5)

Status: **DRAFT for Peter's review.** Implements `DEVELOPMENT_PLAN.md` §4.3 and
the composition half of §4.2. Code: `ampere/core/transform.py`,
`ampere/core/exceptions.py`. Tests: `tests/core/test_transform.py`,
`tests/core/thirdparty_polarimeter.py` (the out-of-tree extension),
`tests/core/test_spec_doctests.py`.

Every worked example below is executed as a doctest by
`tests/core/test_spec_doctests.py`, so this document cannot drift from the
implementation without the suite going red. Examples share one namespace and
build on each other; read them in order.

Where this document and `DEVELOPMENT_PLAN.md` disagree, the plan wins and this
document is wrong — file it as a decision-log correction.

---

## 1. What this contract is for

A model computes physics. An instrument records a number. Everything between
the two — the line-spread function, the resampling onto detector pixels, the
integration through a filter, the calibration scale factor, the response
matrix — is what this contract is about.

Legacy ampere has no concept for it. The one instrumental effect it does model,
a photometric scale factor, is a hard-coded `scaleFac` inside the likelihood;
resampling happens ad hoc, on every call, wherever somebody needed it (issue
#12). There is nowhere to put an LSF and nothing to reuse between two
spectrographs.

The replacement is `DEVELOPMENT_PLAN.md` §4.3: an **`Instrument` is an ordered
chain of small, reusable `Transformation`s** bound to one named channel of a
`ModelResult`. This is a deliberate choice against the alternative, and
`prior_art.md` §3 (lesson 3M1) and §6 (Tension 3) record it: 3ML puts the
response *and* the likelihood inside one opaque per-instrument plugin, which
makes extensibility trivial and reuse impossible — no two of its plugins share
a resampling routine. Ampere factors the chain (here) apart from the
likelihood and noise model (W1.6) so that the LSF written for one spectrograph
is the LSF for the next one. The bet is that real instruments decompose;
W1.11's modality sketches carry the deliberately-awkward stress test that
checks the bet before the spec freeze.

Three things make the chain more than a list of callables, and they are the
three halves of this contract: **kinds are checked when the chain is
assembled**, **transformations may own nuisance parameters**, and
**instruments may negotiate the grids the model evaluates on** — the last
being entirely optional, at every level.

It is **backend-neutral**: numpy, `astropy.units` and stdlib, nothing else
(`architecture.md` §3–4).

### Setup for the examples

```pycon
>>> import numpy as np
>>> import astropy.units as u
>>> import scipy.stats as st
>>> from ampere.core import (
...     DEFAULT_CHANNEL, AxisRequirement, Instrument, Model, ModelResult, Parameter,
...     PhotometricPoints, Spectrum, Transformation, negotiate, propagate_mask,
... )
>>> from ampere.core.exceptions import ChannelError, CompositionError, TransformationError

```

## 2. The objects

| Object | Role |
|---|---|
| `Transformation` | One step: a container in, a container out. The ABC a user subclasses |
| `Instrument` | An ordered chain of steps, bound to one channel by name; the chain is kind-checked when it is built, the channel when it is bound |
| `Model` | A `Parameterised` that produces a `ModelResult`; the other side of the negotiation protocol |
| `AxisRequirement` | What one step needs of one coordinate axis: coverage, required points, sampling density |
| `ChannelRequirements` | The union of those, per channel — the output of `negotiate` and the input to `Model.compile_for` |
| `negotiate` | Collects and unions requirements across instruments |
| `propagate_mask` | The ratified mask-propagation rule (§6) |
| `TransformationError`, `CompositionError` | This contract's errors; `CompositionError` is the assemble-time subset |

`Transformation` is this contract's word, and it is available precisely because
W1.3 named the unconstrained-space map `Bijection` rather than `Transform`
(`parameters.md` §13).

## 3. The simple path: a fixed grid and no negotiation

Nothing below this section is needed to fit data. `DEVELOPMENT_PLAN.md` §4.3
requires that "the simple path must remain: fixed grid, no negotiation", so
that path is shown first, in full, before any machinery appears.

A model on a fixed wavelength grid, written exactly as `parameters.md` §10
writes one — buffers and parameters read out of one `context`:

```pycon
>>> class GreyBody(Model):
...     """Flux proportional to T lambda**-beta on a fixed grid."""
...     def __init__(self, wavelength):
...         self.register_buffer("wavelength", wavelength, unit=u.micron)
...         self.register_buffer("beta", 1.8)
...         self.register_parameter(
...             Parameter("temperature", st.uniform(100.0, 900.0), unit=u.K, value=300.0)
...         )
...     def evaluate(self, **values):
...         ctx = self.context(values)
...         flux = ctx["temperature"] * ctx["wavelength"] ** -ctx["beta"]
...         return Spectrum(ctx["wavelength"] * u.micron, flux * u.Jy)
>>> model = GreyBody(np.geomspace(1.0, 100.0, 6))
>>> model(temperature=300.0)
<ModelResult default: Spectrum[6]>

```

An instrument that applies a flux calibration factor — the nuisance parameter
legacy hard-codes as `scaleFac`:

```pycon
>>> class CalibrationScale(Transformation):
...     """Multiply by a free scale factor."""
...     ACCEPTS = (Spectrum,)
...     def __init__(self, **kwargs):
...         super().__init__(**kwargs)
...         self.register_parameter(Parameter("scale", st.lognorm(0.1), value=1.0))
...     def apply(self, samples, values):
...         return samples.with_values(samples.values * self.context(values)["scale"])
>>> photometer = Instrument([CalibrationScale()])
>>> photometer
<Instrument 'default' on channel 'default': Spectrum -> calibration_scale -> Spectrum>

```

That is the whole composition. Evaluating it is one call, and the chain's
parameters are one `ParameterSet` an inference engine can already consume:

```pycon
>>> predicted = photometer(model(temperature=300.0), {"calibration_scale.scale": 2.0})
>>> np.round(predicted.values[:3], 3).tolist()
[600.0, 114.328, 21.785]
>>> photometer.parameters.free_names, photometer.parameters.free_size
(('calibration_scale.scale',), 1)

```

Nothing was negotiated, and nothing published a requirement:

```pycon
>>> photometer.requirements()
()

```

The channel defaults to `DEFAULT_CHANNEL`, so a single-output model needs no
channel names at all — the same "keep the trivial case trivial" rule W1.4
applies to `ModelResult`:

```pycon
>>> photometer.channel == DEFAULT_CHANNEL
True

```

## 4. `Transformation`: one method, and what you get for free

A subclass declares which container kinds it accepts, optionally which kind it
produces, and implements `apply`. That is the entire extension surface —
`DEVELOPMENT_PLAN.md` §4.3 makes "one ABC, one or two methods, no changes
inside ampere" a first-class requirement, and §11 below discharges it with a
test written out of tree.

| Declaration | Meaning |
|---|---|
| `ACCEPTS` | Tuple of container kinds this step may be handed. Default: any |
| `PRODUCES` | The kind returned, or `None` (the default) for "the same kind it was given" |
| `apply(samples, values)` | The transformation itself. `values` are this step's *local* parameter values |
| `requirements()` | Optional; §7. Empty by default |
| `label` | Component label for this step's parameters. Defaults to the class name in snake_case |

`PRODUCES = None` is the common case, because most instrumental effects are
kind-preserving: convolution, resampling and calibration all take a `Spectrum`
to a `Spectrum`. Synthetic photometry is the other sort — it genuinely changes
the kind:

```pycon
>>> class SyntheticPhotometry(Transformation):
...     """Integrate a spectrum through tabulated filter responses."""
...     ACCEPTS = (Spectrum,)
...     PRODUCES = PhotometricPoints
...     def __init__(self, filters, pivots, response, **kwargs):
...         super().__init__(**kwargs)
...         self.register_buffer("pivots", pivots, unit=u.micron)
...         self.register_buffer("response", response)
...         self._filters = tuple(filters)
...     def apply(self, samples, values):
...         ctx = self.context(values)
...         weights = ctx["response"] / ctx["response"].sum(axis=1, keepdims=True)
...         return PhotometricPoints(
...             self._filters,
...             ctx["pivots"] * u.micron,
...             weights @ samples.values * samples.unit,
...             mask=propagate_mask(samples, weights),
...         )

```

The response matrix is a **buffer**, not a parameter — `architecture.md` §6's
distinguishing question ("would you ever put a prior on it?") answers itself,
and `parameters.md` §10 requires buffers to be declared explicitly rather than
inferred from stray array attributes.

```pycon
>>> response = np.array([[1.0, 1.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0, 1.0, 1.0]])
>>> bands = SyntheticPhotometry(["blue", "red"], np.array([1.6, 45.0]), response)
>>> bands.buffers.names, bands.parameters.names
(('pivots', 'response'), ())
>>> observed = bands(model(temperature=300.0).single())
>>> observed.filters.tolist(), np.round(observed.values, 3).tolist()
(['blue', 'red'], [178.582, 0.849])

```

Call the *instance*, never `apply` directly: `__call__` is where this contract's
checks live. Handing a step a kind it does not accept is a `CompositionError`
that names both kinds and both ways out:

```pycon
>>> bands(observed)
Traceback (most recent call last):
    ...
ampere.core.exceptions.CompositionError: SyntheticPhotometry accepts Spectrum, but was given a PhotometricPoints. Either bind this step to a channel of an accepted kind, or put a transformation that produces one earlier in the chain.

```

And a step that does not return what it advertises is caught immediately,
rather than becoming a confusing shape error inside a likelihood:

```pycon
>>> class Miscoded(Transformation):
...     ACCEPTS = (Spectrum,)
...     PRODUCES = PhotometricPoints
...     def apply(self, samples, values):
...         return samples
>>> Miscoded()(model(temperature=300.0).single())
Traceback (most recent call last):
    ...
ampere.core.exceptions.TransformationError: Miscoded.apply returned a Spectrum, but the class declares PRODUCES=PhotometricPoints, so given a Spectrum it must return a PhotometricPoints. Fix apply, or declare PRODUCES = Spectrum.

```

### Nuisance parameters

A transformation carrying parameters inherits `Parameterised` and declares
through it — the obligation `parameters.md` §13 places on this document. There
is no second parameter mechanism for instruments, which means tying, fixing,
buffer promotion, priors, bijections and serialisation all work on a
calibration factor exactly as they work on a temperature:

```pycon
>>> scale = CalibrationScale()
>>> scale.parameters["scale"].is_free, scale.label
(True, 'calibration_scale')
>>> scale.parameters["scale"].unconstraining_bijection()
Log(lower=0.0)

```

A transformation with no parameters at all is perfectly ordinary — an LSF whose
width is known from the instrument's own calibration is a buffer, not a
parameter, and promoting it later is a configuration change (`parameters.md`
§10), not a rewrite.

## 5. `Instrument`: a chain bound to a channel

An `Instrument` is the composition. It owns three facts: the channel it binds,
the kind it expects there, and the ordered steps.

```pycon
>>> spectrograph = Instrument(
...     [CalibrationScale(), bands], channel="sed_lowres", label="wise"
... )
>>> spectrograph
<Instrument 'wise' on channel 'sed_lowres': Spectrum -> calibration_scale -> synthetic_photometry -> PhotometricPoints>
>>> spectrograph.input_kind.__name__, spectrograph.output_kind.__name__
('Spectrum', 'PhotometricPoints')

```

**The chain is checked when it is built.** Reordering those two steps is a
type error, and it is reported at composition time, naming the step, its
position and what it wanted:

```pycon
>>> Instrument([bands, CalibrationScale()], channel="sed_lowres", label="wise")
Traceback (most recent call last):
    ...
ampere.core.exceptions.CompositionError: instrument 'wise' cannot be composed: step 0 ('synthetic_photometry') produces a PhotometricPoints, but step 1 ('calibration_scale', CalibrationScale) accepts Spectrum. Insert a step that produces a Spectrum, or reorder the chain.

```

There are two moments called "composition time" here, and it is worth being
precise about which is which. The chain's *internal* consistency is checked
when the `Instrument` is constructed, above, because everything it needs is
known then. The *channel's* kind can only be checked against a real
`ModelResult`, so it is checked on binding — `Instrument.bind` is public
precisely so that a caller assembling a fit (W1.7) can do it once, at setup,
rather than discovering it a thousand samples in.

The binding goes through `ModelResult.require(name, kind)`, never
`result[name]` — `results_schema.md` §16 asks for exactly that, so a missing
channel or a wrong kind raises `ChannelError` naming what *is* available:

```pycon
>>> spectrograph(model(temperature=300.0))
Traceback (most recent call last):
    ...
ampere.core.exceptions.ChannelError: no channel named 'sed_lowres'. Available channels: 'default' (Spectrum). This result has only the automatic 'default' channel, which means the model returned a bare container; name its channels explicitly if an instrument needs to bind to one by name.

```

Name the channel and it binds:

```pycon
>>> named = ModelResult({"sed_lowres": model(temperature=300.0).single()})
>>> spectrograph(named, {"calibration_scale.scale": 1.0}).filters.tolist()
['blue', 'red']

```

### Parameters compose by merging

A chain's parameters are its steps', merged with each step's label as the
component (`parameters.md` §13's instruction), so two independently written
steps never collide and every parameter is traceable to the step that owns it:

```pycon
>>> spectrograph.parameters.names
('calibration_scale.scale',)
>>> spectrograph.mapping.components
('calibration_scale', 'synthetic_photometry')
>>> spectrograph.mapping.distribute({"calibration_scale.scale": 1.5})
{'calibration_scale': {'scale': 1.5}, 'synthetic_photometry': {}}

```

That is also how values reach the steps at evaluation time, and why the
instrument accepts either a mapping over qualified names or a flat
free-parameter vector in `parameters` order:

```pycon
>>> np.round(spectrograph(named, [2.0]).values, 3).tolist()
[357.164, 1.698]

```

Two steps of the same class in one chain therefore need distinct labels, and
this is enforced rather than patched over with automatic numbering:

```pycon
>>> Instrument([CalibrationScale(), CalibrationScale()], label="twice")
Traceback (most recent call last):
    ...
ampere.core.exceptions.CompositionError: instrument 'twice' has two steps labelled 'calibration_scale' (positions 0 and 1). Step labels become component names when the chain's parameters are merged, so they must be unique: pass label='...' to one of them. They are not numbered automatically, because inserting a step would then silently rename every parameter after it.
>>> both = Instrument(
...     [CalibrationScale(label="detector"), CalibrationScale(label="aperture")], label="twice"
... )
>>> both.parameters.free_names
('detector.scale', 'aperture.scale')

```

An instrument's own label is the component under which its parameters join a
joint fit (W1.7); it defaults to the channel name. `Instrument` is deliberately
**not** `Parameterised` itself — its parameters are its steps' — so an
instrument-level nuisance parameter is a one-step transformation, which is
exactly what `CalibrationScale` already is.

### Inferring the kind, and when you must say

The kind on the channel is inferred from the first step when that step accepts
exactly one kind, which covers nearly everything. When it does not, guessing
would be the wrong answer:

```pycon
>>> class PassThrough(Transformation):
...     ACCEPTS = (Spectrum, PhotometricPoints)
...     def apply(self, samples, values):
...         return samples
>>> Instrument([PassThrough()], channel="sed_lowres")
Traceback (most recent call last):
    ...
ampere.core.exceptions.CompositionError: instrument 'sed_lowres' cannot infer the kind it binds: its first step ('pass_through', PassThrough) accepts Spectrum or PhotometricPoints, so there is no single answer. Pass input_kind=... to say which kind channel 'sed_lowres' holds.
>>> Instrument([PassThrough()], channel="sed_lowres", input_kind=Spectrum).output_kind.__name__
'Spectrum'

```

An instrument with no steps at all is legal: a model that already produces the
observable needs a channel binding and nothing else.

```pycon
>>> Instrument(channel="sed_lowres", input_kind=Spectrum)
<Instrument 'sed_lowres' on channel 'sed_lowres': Spectrum -> (no steps) -> Spectrum>

```

## 6. Masks: the rule this contract ratifies

`results_schema.md` §15.9 leaves mask propagation to this contract, and §16
proposes a rule for the many-to-one case. **This contract ratifies it**: an
output sample is masked if **any** input sample that influences it is masked.

The reasoning is `results_schema.md` §7's definition of what a mask means. A
masked sample carries *exactly zero* information — it is to be treated as if it
had never been observed. An output bin that integrates over such a sample
therefore has no defensible value, and the alternative rule (mask only if
*every* contributing input is masked) silently feeds a partly-invalid number
into a likelihood, where it is indistinguishable from a good one. Being
conservative loses a little data at the edges of a masked region; being
permissive corrupts the fit. `propagate_mask` is that rule:

```pycon
>>> propagate_mask(None) is None
True
>>> propagate_mask(np.array([False, True, False])).tolist()
[False, True, False]

```

For a many-to-one step, pass the same weights the values go through — a
resampling matrix, a response matrix, a boolean adjacency. Non-zero weight
means influence:

```pycon
>>> weights = np.array([[0.5, 0.5, 0.0, 0.0], [0.0, 0.0, 0.5, 0.5]])
>>> propagate_mask(np.array([False, False, True, False]), weights).tolist()
[False, True]

```

An obligation that is only documented is not an obligation, so
`Transformation.__call__` enforces it: a step whose input carried a mask and
whose output carries none is a bug, caught where it happens.

```pycon
>>> class Forgetful(Transformation):
...     ACCEPTS = (Spectrum,)
...     def apply(self, samples, values):
...         return Spectrum(samples.spectral_axis.quantity(), samples.values * u.Jy)
>>> flagged = Spectrum(
...     [1.0, 2.0, 4.0] * u.um, [3.0, 2.5, 1.0] * u.Jy, mask=np.array([False, True, False])
... )
>>> Forgetful()(flagged)
Traceback (most recent call last):
    ...
ampere.core.exceptions.TransformationError: Forgetful.apply dropped the mask: its input excluded 1 of 3 sample(s) and its output carries no mask. A masked sample carries zero information and that must survive the chain (results_schema.md §16). Use propagate_mask(samples, weights) and pass it as mask=..., or pass an explicit all-False mask if the output genuinely does not depend on the excluded samples.

```

The check is cheap to satisfy and free in the common case: `with_values`
inherits the template's mask, so a kind-preserving step gets propagation
without writing anything.

```pycon
>>> CalibrationScale()(flagged, {"scale": 2.0}).mask.tolist()
[False, True, False]

```

The escape hatch is explicit, not implicit. A step whose output genuinely does
not depend on any excluded sample says so by passing an all-`False` mask —
which is what `propagate_mask` returns for a filter that misses the masked
region:

```pycon
>>> narrow = np.array([[1.0, 0.0, 0.0]])
>>> propagate_mask(flagged, narrow).tolist()
[False]

```

## 7. Requirements negotiation

`DEVELOPMENT_PLAN.md` §4.3 decided this feature, and decided its shape: an
instrument may publish what it needs; a one-off compilation step before the hot
loop lets the model configure per-channel evaluation grids from the union;
**models are free to ignore requests**. It exists because evaluating the fine
grid everywhere is prohibitive and evaluating the coarse grid at a line is
useless (the §4.2 pattern, worked in §8 below), because it removes the
per-call resampling of issue #12, and because it is the natural place to cache.

It is optional at every level. A transformation publishes nothing by default; a
model ignores requirements by default; §3 never mentioned any of it.

### What a requirement is

A requirement is a statement about **coordinates**, because a container is
built coordinates-first — the shape `results_schema.md` §16 asks for. It has
two independent halves:

| Half | Fields | Meaning |
|---|---|---|
| *where* | `intervals`, `points` | coverage to span; coordinates that must be present exactly |
| *how finely* | `max_step`, `min_resolving_power` | absolute spacing; λ/Δλ, the spectrograph's own language |

A density belongs to the coverage declared **alongside** it — the `@` in the
repr is that association, and it is what makes the union below usable.

```pycon
>>> band = AxisRequirement(
...     "spectral_axis", intervals=(1.0, 30.0) * u.um, min_resolving_power=40.0, source="wise"
... )
>>> band
<AxisRequirement 'spectral_axis' um [1,30]@R>=40>
>>> band.intervals, band.constrains_density
(((1.0, 30.0),), True)

```

A density with no coverage to apply to is therefore refused, rather than
quietly becoming a constraint on whatever anybody else asks for:

```pycon
>>> AxisRequirement("spectral_axis", min_resolving_power=1000.0)
Traceback (most recent call last):
    ...
ampere.core.exceptions.TransformationError: the requirement on axis 'spectral_axis' declares a sampling density but no coverage for it to apply to. A density belongs to the intervals declared alongside it — otherwise, once requirements are unioned, there is no way to say where it holds. Add intervals=(low, high).

```

Units are converted once, here, exactly as `results_schema.md` §6 converts
them once at construction — a requirement in nanometres and one in microns are
the same requirement:

```pycon
>>> AxisRequirement("spectral_axis", intervals=(2000.0, 30000.0) * u.nm).convert_to(u.um).intervals
((2.0, 30.0),)

```

A requirement can build a reference grid satisfying itself, so that the common
case needs no arithmetic from the model author:

```pycon
>>> grid = band.coordinates()
>>> grid.unit, grid.size
(Unit("um"), 139)
>>> grid.value[0], grid.value[-1]
(np.float64(1.0), np.float64(30.0))
>>> float(np.max(grid[1:] / grid[:-1])) <= 1.0 + 1.0 / 40.0
True

```

### The union

Two instruments wanting different things from one channel is the whole point.
Coverage and required points accumulate, and **each density stays attached to
the coverage that asked for it**:

```pycon
>>> deep = AxisRequirement(
...     "spectral_axis", intervals=(20.0, 200.0) * u.um, max_step=2.0, source="pacs"
... )
>>> both_bands = band.union(deep)
>>> both_bands
<AxisRequirement 'spectral_axis' um [1,30]@R>=40 [20,200]@step<=2>
>>> both_bands.segments()
((1.0, 30.0, None, 40.0), (20.0, 200.0, 2.0, None))
>>> both_bands.source
'wise, pacs'

```

That the two intervals overlap and were *not* merged is the point, not an
oversight. Merging them would mean applying the stricter of R=40 and
Δλ ≤ 2 µm across the whole 1–200 µm span, which is a different — and much more
expensive — request than either instrument made. The pathological case is easy
to reach: a broad SED at R=40 and one 0.1 µm line window at 0.0025 µm sampling
would become 0.0025 µm sampling from 1 to 200 µm, some eighty thousand
coordinates in place of a couple of hundred, defeating the purpose §4.3 gives
negotiation. Per-interval densities are what prevent it:

```pycon
>>> broad = AxisRequirement(
...     "spectral_axis", intervals=(1.0, 200.0) * u.um, min_resolving_power=40.0
... )
>>> window = AxisRequirement("spectral_axis", intervals=(866.9, 867.0) * u.um, max_step=0.0025)
>>> broad.union(window).coordinates().size
258

```

Intervals asking for the *same* density do merge, because then there is nothing
to lose by it:

```pycon
>>> coarse = AxisRequirement("spectral_axis", intervals=(1.0, 10.0) * u.um, max_step=0.5)
>>> more = AxisRequirement("spectral_axis", intervals=(5.0, 20.0) * u.um, max_step=0.5)
>>> coarse.union(more).intervals
((1.0, 20.0),)

```

`max_step` and `min_resolving_power` are **constructor-only** (ruled by Peter,
2026-09-02 — §15.8): they are folded into the per-interval densities at
construction and do not survive as attributes. On a merged requirement a
single scalar could only be a strictest-anywhere summary, and a reader who
trusted it would overestimate what was asked of every other interval.
`segments()` is the canonical — and only — public form, and it is what
`coordinates()` reads:

```pycon
>>> both_bands.max_step
Traceback (most recent call last):
    ...
AttributeError: ...

```

Disjoint intervals stay disjoint — the case that makes negotiation worth having
at all, because a union that closed the gap between two narrow windows would
force the model to evaluate finely across hundreds of microns of nothing:

```pycon
>>> line_32 = AxisRequirement("spectral_axis", intervals=(866.9, 867.0) * u.um, max_step=0.005)
>>> line_21 = AxisRequirement("spectral_axis", intervals=(1300.3, 1300.5) * u.um, max_step=0.005)
>>> lines = line_32.union(line_21)
>>> lines.intervals
((866.9, 867.0), (1300.3, 1300.5))
>>> lines.coordinates().size
64

```

Requirements ampere cannot reconcile are refused rather than guessed at. A
unitless requirement and a unit-bearing one on the same axis are a composition
error, not a pair to fix silently:

```pycon
>>> band.union(AxisRequirement("spectral_axis", intervals=(1.0, 2.0)))
Traceback (most recent call last):
    ...
ampere.core.exceptions.CompositionError: requirements on axis 'spectral_axis' disagree about units: one is in no unit and another in um. Declare a unit on both (or neither) — ampere will not assume a bare number is in the other's unit.

```

### Publishing, and collecting

A step publishes by overriding `requirements()`. A resampler knows the
detector's own wavelengths at construction, so it can say what it needs of the
model in the model's own coordinates:

```pycon
>>> class NearestResampler(Transformation):
...     """A deliberately naive binning resampler; Phase 2 writes the real one."""
...     ACCEPTS = (Spectrum,)
...     def __init__(self, target, **kwargs):
...         super().__init__(**kwargs)
...         self._target = np.asarray(target, dtype=float)
...     def requirements(self):
...         return (
...             AxisRequirement(
...                 "spectral_axis",
...                 unit=u.um,
...                 intervals=(self._target[0], self._target[-1]),
...                 max_step=float(np.diff(self._target).min()) / 2.0,
...             ),
...         )
...     def _weights(self, source):
...         nearest = np.abs(self._target[None, :] - source[:, None]).argmin(axis=1)
...         hit = nearest[None, :] == np.arange(self._target.size)[:, None]
...         return hit / np.maximum(hit.sum(axis=1, keepdims=True), 1)
...     def apply(self, samples, values):
...         weights = self._weights(samples.spectral_axis.values)
...         return Spectrum(
...             self._target * u.micron,
...             weights @ samples.values * samples.unit,
...             mask=propagate_mask(samples, weights),
...         )

```

`negotiate` collects across instruments, groups by channel, checks that
instruments binding one channel agree about its kind, and unions their
requirements per axis:

```pycon
>>> detector = Instrument(
...     [NearestResampler(np.linspace(5.0, 25.0, 21)), CalibrationScale()],
...     channel="sed_lowres",
...     label="irs",
... )
>>> requirements = negotiate([detector, spectrograph])
>>> sorted(requirements)
['sed_lowres']
>>> asked = requirements["sed_lowres"]
>>> asked
<ChannelRequirements 'sed_lowres' Spectrum axes=['spectral_axis'] from ['irs', 'wise']>
>>> asked["spectral_axis"]
<AxisRequirement 'spectral_axis' um [5,25]@step<=0.5>

```

Every bound channel appears, even one with no requirements at all: "this
channel is consumed" is itself worth telling a model. And a requirement naming
an axis the channel's kind does not have is a composition error, because it is
always a mistake:

```pycon
>>> class WrongAxis(Transformation):
...     ACCEPTS = (Spectrum,)
...     def requirements(self):
...         return (AxisRequirement("time", intervals=(0.0, 1.0) * u.day),)
...     def apply(self, samples, values):
...         return samples
>>> negotiate([Instrument([WrongAxis()], channel="sed_lowres", label="clock")])
Traceback (most recent call last):
    ...
ampere.core.exceptions.CompositionError: instrument 'clock' published a requirement on axis 'time', but the Spectrum on channel 'sed_lowres' is indexed by ['spectral_axis']. Requirements name the container's own axes.

```

### Compiling, once

`Model.compile_for` is the one-off step. Its default implementation is one
line — `return self` — and that is the contract, not a placeholder: a model
that ignores every request must remain a working model.

```pycon
>>> model.compile_for(requirements) is model
True
>>> model(temperature=300.0).single().n_samples          # unchanged: 6, as declared
6

```

A model that honours requests builds its evaluation grids from
`coordinates()`, keeps the resulting containers as **templates**, and refills
them with `with_values` on every evaluation. That compile-once/evaluate-many
split is the one `results_schema.md` §16 specifies, and it is where the caching
`DEVELOPMENT_PLAN.md` §4.3 wants naturally lives.

## 8. Worked example: a low-resolution SED and high-resolution CO windows

This is the pattern `DEVELOPMENT_PLAN.md` §4.2 uses to justify named channels
and §4.3 uses to justify negotiation, and it is this item's headline criterion.
A dusty envelope emits a broad continuum *and* two narrow CO lines. A
photometer wants coverage from 1 to 200 µm at R ≈ 40. A heterodyne receiver
wants 0.005 µm sampling in two 0.1 µm windows, 430 µm apart. Neither grid is
usable for the other.

The model produces two channels and lets its grids be configured:

```pycon
>>> class DustyEnvelope(Model):
...     """Continuum on one channel, CO lines on another; grids are negotiable."""
...     LINES = (866.96, 1300.40)
...     def __init__(self):
...         self.register_parameter(
...             Parameter("temperature", st.uniform(50.0, 450.0), unit=u.K, value=300.0)
...         )
...         self.register_parameter(Parameter("line_flux", st.uniform(0.0, 5.0), value=1.0))
...         self._templates = {
...             "sed_lowres": self._template(np.geomspace(1.0, 200.0, 12)),
...             "co_windows": self._template(np.linspace(866.9, 867.0, 8)),
...         }
...     @staticmethod
...     def _template(wavelength):
...         return Spectrum(np.asarray(wavelength) * u.micron, np.zeros(len(wavelength)) * u.Jy)
...     def compile_for(self, requirements):
...         for channel, asked in requirements.items():
...             if "spectral_axis" in asked:
...                 self._templates[channel] = self._template(
...                     asked["spectral_axis"].coordinates().to_value(u.micron)
...                 )
...         return self
...     def evaluate(self, **values):
...         out = {}
...         for channel, template in self._templates.items():
...             lam = template.spectral_axis.values
...             flux = values["temperature"] * lam**-1.8
...             for centre in self.LINES:
...                 flux = flux + values["line_flux"] * np.exp(-0.5 * ((lam - centre) / 0.02) ** 2)
...             out[channel] = template.with_values(flux)
...         return ModelResult(out)

```

Two instruments, each bound to the channel it needs, each publishing what it
needs of that channel:

```pycon
>>> class Coverage(Transformation):
...     """Publishes a requirement and does nothing else; the LSF goes here in Phase 2."""
...     ACCEPTS = (Spectrum,)
...     def __init__(self, requirement, **kwargs):
...         super().__init__(**kwargs)
...         self._requirement = requirement
...     def requirements(self):
...         return (self._requirement,)
...     def apply(self, samples, values):
...         return samples
>>> sed_instrument = Instrument(
...     [
...         Coverage(
...             AxisRequirement(
...                 "spectral_axis", intervals=(1.0, 200.0) * u.um, min_resolving_power=40.0
...             )
...         ),
...         CalibrationScale(),
...     ],
...     channel="sed_lowres",
...     label="photometer",
... )
>>> heterodyne = Instrument(
...     [
...         Coverage(
...             AxisRequirement(
...                 "spectral_axis",
...                 intervals=[(866.9, 867.0), (1300.3, 1300.5)] * u.um,
...                 max_step=0.005,
...             )
...         ),
...         NearestResampler(np.linspace(866.92, 866.99, 15)),
...     ],
...     channel="co_windows",
...     label="heterodyne",
... )

```

Negotiation, once:

```pycon
>>> asked = negotiate([sed_instrument, heterodyne])
>>> sorted(asked)
['co_windows', 'sed_lowres']
>>> asked["sed_lowres"]["spectral_axis"]
<AxisRequirement 'spectral_axis' um [1,200]@R>=40>
>>> asked["co_windows"]["spectral_axis"]
<AxisRequirement 'spectral_axis' um [866.9,867]@step<=0.005 [866.92,866.99]@step<=0.0025 [1300.3,1300.5]@step<=0.005>

```

The heterodyne requirement is itself a union of two steps' needs: the windows
and their nominal sampling come from the receiver's tuning (the `Coverage`
step), and the finer sampling over the part of the first window the detector
actually reads comes from the resampler. Neither step had to know about the
other, and the resampler's 0.0025 µm did not leak onto the second window.

The channel also keeps its **gap**. That is the whole argument for
negotiating coordinates rather than a range and a resolution: the union of two
narrow windows is two narrow windows, and the model never evaluates across the
430 µm of continuum between them.

```pycon
>>> envelope = DustyEnvelope().compile_for(asked)
>>> result = envelope(temperature=300.0, line_flux=2.0)
>>> result
<ModelResult sed_lowres: Spectrum[216], co_windows: Spectrum[94]>
>>> float(np.diff(result["co_windows"].spectral_axis.values).max()) > 400.0
True

```

Both instruments now evaluate against a grid built for them:

```pycon
>>> continuum = sed_instrument(result, {"calibration_scale.scale": 1.0})
>>> line = heterodyne(result)
>>> continuum.n_samples, line.n_samples
(216, 15)

```

And the hot loop refills the templates rather than rebuilding them, so the
coordinate validation, the unit checks and the advertised structure are all
computed once (`results_schema.md` §10):

```pycon
>>> again = envelope(temperature=250.0, line_flux=2.0)
>>> again["sed_lowres"].axes is result["sed_lowres"].axes
True
>>> round(float(again["sed_lowres"].values[0] / result["sed_lowres"].values[0]), 6)
0.833333

```

The same model, uncompiled, still works — on its own declared grids. That is
what "models are free to ignore requests" has to mean in practice:

```pycon
>>> DustyEnvelope()(temperature=300.0, line_flux=2.0)
<ModelResult sed_lowres: Spectrum[12], co_windows: Spectrum[8]>

```

## 9. `Model`: why the ABC is here

`parameters.md` left the question open ("a future `Model` ABC (W1.5/W1.7)").
It is settled here, and the answer is **here, minimally**.

The reason is that negotiation is a *protocol*, and a protocol with one side
missing is a suggestion. `compile_for` is the model's half of §4.3's decided
feature; it has to live with the requirements types it consumes, and it has no
sensible home in W1.7's `Dataset`/`FittingProblem` vocabulary. So `Model` is
defined here as exactly two things — a `Parameterised` that returns a
`ModelResult`, plus the `compile_for` hook — and nothing else.

Everything engine-facing stays W1.7's, per `DEVELOPMENT_PLAN.md` §4.5, which
puts `log_prob`, `prior_transform`, `simulate` and the capability flags on the
*fitting problem* (model + instruments + datasets), not on the model. W1.7
therefore composes this ABC rather than redefining it, and nothing here needs
revisiting when it does.

`Model.__call__` normalises what an engine has to hand — a flat vector, a
mapping, keywords, or nothing — and files a bare container under the default
channel, so the single-output case never has to name anything:

```pycon
>>> model([400.0]).single().n_samples           # a flat free-parameter vector
6
>>> model({"temperature": 400.0}) == model(temperature=400.0)
True
>>> model()                                     # declared values
<ModelResult default: Spectrum[6]>

```

`__call__` also attaches the resolved values to the result it returns
(`results_schema.md` §3 — the `(θ, result)` pairing ruled 2026-09-01), so
every evaluation is a self-contained training-set/provenance pair with no
effort from the model author; a record `evaluate()` attached itself is
respected:

```pycon
>>> model(temperature=400.0).parameters
mappingproxy({'temperature': 400.0})

```

Mixing a vector with keywords is refused rather than resolved by precedence,
because either interpretation would be a silent surprise:

```pycon
>>> model([400.0], temperature=500.0)
Traceback (most recent call last):
    ...
ampere.core.exceptions.TransformationError: GreyBody was called with both a flat parameter vector and the keyword value(s) ['temperature']. A vector already fixes every free parameter; pass one or the other.

```

## 10. The standard library: interfaces now, implementations in Phase 2

`DEVELOPMENT_PLAN.md` §4.3 names the transformations ampere ships. This
contract specifies the *shape* each takes; the physics is Phase 2's, and
`ampere/core/transform.py` deliberately ships **no** concrete transformation —
core is the vocabulary, not the library.

| Transformation | Accepts → produces | Parameters (typical) | Publishes |
|---|---|---|---|
| LSF convolution | `Spectrum` → `Spectrum` | resolving power or FWHM, usually a buffer | coverage padded beyond the output range by several kernel widths; `min_resolving_power` several times the kernel's |
| Spectral resampling | `Spectrum` → `Spectrum` | none | `intervals` spanning the detector, `max_step` finer than its finest bin |
| Synthetic photometry | `Spectrum` → `PhotometricPoints` | none (filter curves are buffers) | `intervals` per filter's non-zero support, `max_step` from the response's own sampling |
| Calibration / scale factor | any → same | the scale factor, or per-channel offsets as an array-valued parameter | nothing |
| Epoch sampling | `TimeSeries` → `TimeSeries` | none | `points` at the observed epochs |
| Fourier sampling | `Image` → `VisibilitySet` | none | `intervals` on `x`/`y` from the field of view, `max_step` from the longest baseline |
| Response matrix (RMF/ARF) | `Spectrum` → `Spectrum` | none (the matrix is a buffer) | `points` at the matrix's own tabulated energies |

Two of them are worth a note. A **response matrix is a matrix multiply**, so
X-ray forward folding is an ordinary `Transformation` and needs nothing special
— `results_schema.md` §16 asks W1.11 to confirm that `Spectrum` with an energy
axis suffices, and this table is the claim it should check. And **Fourier
sampling changes both the kind and the meaning of the axes**, which is why it
is also the one whose requirements cannot be pulled back through a chain
(§13.1).

## 11. User extensibility, proved out of tree

`DEVELOPMENT_PLAN.md` §4.3 states the acceptance criterion: the test suite must
include an extension written *out-of-tree*, imported as if third-party, to
prove the interfaces suffice. `tests/core/thirdparty_polarimeter.py` is that
module. It imports only `ampere.core`'s public API — no private helper, no
internal module — and it defines a container kind, a transformation with a
nuisance parameter and published requirements, and an instrument factory.
`tests/core/test_transform.py::TestOutOfTreeExtension` runs it through a chain.

The shape of it, in miniature: a user's transformation is a class body and one
method.

```pycon
>>> class AtmosphericExtinction(Transformation):
...     """Somebody else's package, subclassing ampere and nothing more."""
...     ACCEPTS = (Spectrum,)
...     def __init__(self, optical_depth, **kwargs):
...         super().__init__(**kwargs)
...         self.register_buffer("tau", optical_depth)
...         self.register_parameter(Parameter("airmass", st.uniform(1.0, 2.0), value=1.2))
...     def apply(self, samples, values):
...         ctx = self.context(values)
...         return samples.with_values(samples.values * np.exp(-ctx["tau"] * ctx["airmass"]))
>>> ground = Instrument([AtmosphericExtinction(np.full(6, 0.1))], label="ground")
>>> ground.parameters.free_names
('atmospheric_extinction.airmass',)
>>> np.round(ground(model(temperature=300.0), {"atmospheric_extinction.airmass": 2.0})
...          .values[0], 4)
np.float64(245.6192)

```

Nothing in `ampere.core` knows that class exists, and nothing had to change to
accommodate it — which is the property `prior_art.md` lesson 3M1 says 3ML's
monolithic plugin buys trivially and a factored design has to earn.

## 12. Decisions and their reasoning

| Decision | Reasoning |
|---|---|
| `Instrument` is a chain of `Transformation`s, not one opaque plugin | `prior_art.md` 3M1/Tension 3: 3ML's plugin buys extensibility at the cost of zero reuse across instruments, and a library of shared LSF/resampling/photometry steps is exactly what ampere wants |
| Kinds are declared (`ACCEPTS`/`PRODUCES`) and checked when the chain is built | §4.2 asks for loud mismatches at composition time; a shape error inside a likelihood ten frames down is the failure mode being avoided |
| `PRODUCES = None` means "same kind in, same kind out" | Kind-preserving is the common case; making it the default keeps the declaration to one line for convolution, resampling and calibration alike |
| Binding goes through `ModelResult.require`, never `result[name]` | `results_schema.md` §16's instruction; `require` is what makes a wrong kind a `ChannelError` instead of a broadcasting surprise |
| `Transformation.__call__` wraps `apply` with the checks | An obligation the base class can check is one the user cannot forget; `apply` stays the one method a subclass writes |
| Transformations declare parameters through `Parameterised` | `parameters.md` §13; one parameter mechanism means tying, fixing, priors, bijections and buffer promotion all work on a calibration factor for free |
| A chain's parameters compose via `ParameterSet.merge`, step label as component | Also `parameters.md` §13. It makes every nuisance parameter traceable to the step that owns it, and independently written steps cannot collide |
| `Instrument` is **not** `Parameterised` | Its parameters are its steps', merged — a `ParameterMapping`, not a `ParameterSet`. An instrument-level nuisance parameter is a one-step transformation, which is what a calibration factor already is |
| The merge is recomputed on access, not cached at construction | A step may be reconfigured (`promote_buffer`) after the chain is built, and a stale snapshot is the silent-drift bug class these contracts exist to end. The hot loop is inside `apply`, not here |
| Duplicate step labels raise; they are not auto-numbered | Numbering makes a parameter's name depend on its position, so inserting a step silently renames everything after it — and those names go into priors, provenance and ArviZ coordinates |
| Mask rule: an output touching **any** masked input is masked | `results_schema.md` §7 defines a masked sample as carrying zero information. Conservative loses a little data at a masked edge; permissive feeds a partly-invalid number to a likelihood, indistinguishable from a good one |
| Dropping a mask raises rather than warning | `results_schema.md` §16 makes propagation an obligation; an unenforced obligation is documentation. The escape hatch — an explicit all-`False` mask — is one argument |
| Requirements are expressed as coordinates | `results_schema.md` §16: containers are built coordinates-first, and `with_values` refills them. Anything else would need translating before it could be used |
| Requirements are declared in the *channel's* coordinates; no pull-back through the chain | A general pull-back needs every step to invert its own coordinate map, and a kind-changing step (photometry, Fourier sampling) cannot in general. The standard library's steps all know their targets at construction, so they can say what they need directly (§13.1) |
| `negotiate` unions coverage; density stays attached to the interval that asked for it | The weakest grid satisfying everybody is the one that lets one evaluation feed several instruments. A single strictest-everywhere density is *not* that grid: a broad SED at R=40 unioned with one 0.1 µm line window at 0.0025 µm sampling would demand eighty thousand coordinates instead of a couple of hundred, defeating the purpose §4.3 gives negotiation. Disjoint intervals stay disjoint for the same reason |
| `max_step`/`min_resolving_power` survive a union as a strictest-anywhere *summary*; `segments()` is canonical | The scalar reads naturally on a freshly declared requirement and is useful for reporting; keeping it as the thing grids are built from is what would reintroduce the over-refinement above |
| A density declared with no coverage raises | It has nowhere to apply. Silently dropping it in a union loses a real constraint; silently applying it to everyone else's coverage is the over-refinement bug by another route |
| `compile_for` defaults to `return self` | §4.3: models are free to ignore requests, and a model that does must still be a working model. Making the default a no-op is what keeps the simple path simple |
| Units on a requirement are converted once, and unitless never mixes with unit-bearing | `results_schema.md` §6's rule, and the units trap in `DEVELOPMENT_PLAN.md` §7. Assuming a bare number shares the other requirement's unit is exactly the silent factor-of-1000 to avoid |
| No spectral equivalencies in requirement unit conversion | Converting a wavelength interval to frequency reverses it and turns even spacing into uneven; a loud refusal beats a subtly wrong grid |
| `Model` ABC lives here, minimally | Negotiation needs a model side, and `compile_for` belongs with the types it consumes. Everything engine-facing (`log_prob`, `simulate`, capability flags) is W1.7's per §4.5 (§9) |
| `Model.__call__` accepts a vector, a mapping, keywords or nothing; refuses vector + keywords | Engines have vectors, users have keywords, and both are wanted. Silently letting one override the other is the kind of precedence rule nobody remembers |
| Core ships **no** concrete transformation | §4.3's standard library is Phase 2 physics; core is the vocabulary the backends share. The interfaces are specified in §10 so Phase 2 has a target |

## 13. Deliberate limitations of v1.5

Each is a decision, not an oversight. Each has an extension point.

1. **No requirement pull-back through a chain.** A step's `requirements()` are
   read as statements about the *channel's* coordinates, not about that step's
   own input. For the standard library this is exact — every one of those steps
   knows its target grid, its filter curves or its response matrix at
   construction — but a chain whose third step needs something only expressible
   in the second step's output coordinates cannot say so. The extension point
   is a `Transformation.pull_back(requirement)` method with an identity default
   for kind-preserving steps and a loud refusal for kind-changing ones.
2. **Density is piecewise-constant, per interval.** A requirement cannot
   express "either of these two grids will do", nor a density that varies
   continuously across an interval. Nor can it express "no coarser than R=40
   *and* no finer than R=200" — there is no upper bound on sampling, only a
   lower one, so a model is always free to over-sample. Splitting an interval
   is the way to vary density, and the union keeps the pieces separate.
   Overlapping intervals of different density are also kept separate rather
   than being split at their boundaries: the resulting grid is the union of
   both, which satisfies both and costs only the coordinates in the overlap.
3. **`negotiate` does not build containers.** It produces coordinates; the
   model builds the container, because only the model knows the value unit,
   the fidelity tag and (for `PhotometricPoints` or a user kind) what else the
   constructor needs. `ChannelRequirements.template()` is the obvious extension
   if that turns out to be boilerplate everybody writes.
4. **No partial-coverage tolerance in the mask rule.** A resampled bin is
   masked if it touches one masked input out of five hundred. A transformation
   is free to implement a coverage threshold and document it, but it must not
   be silently more permissive than `propagate_mask`; a `min_valid_fraction`
   argument is the natural extension if real pipelines want it.
5. **Mask propagation is 1-D.** `propagate_mask`'s influence matrix is
   `(n_out, n_in)`, which fits `Layout.POINTS`. A `Layout.GRID` container's
   mask must be flattened by the transformation itself.
6. **The influence matrix is dense.** `(n_out, n_in)` booleans are built in
   full. For the negotiated grids this contract is designed around (hundreds to
   thousands of samples) that is nothing; a 10⁵-sample resample would want a
   sparse or banded path, which is a Phase 2 backend concern.
7. **Chains are linear.** No branching, no merging of two channels into one
   prediction. A model producing two channels feeds two instruments; an
   instrument that genuinely needs two model channels at once is not
   expressible, and would want an `Instrument` that binds a tuple of channels.
8. **No instrument-level parameters.** Deliberate (§5); the workaround is a
   one-step transformation, which is what a calibration factor is anyway.
9. **`compile_for` has no failure mode.** A model that *cannot* satisfy a
   requirement has no way to say so other than ignoring it or raising. Whether
   a negotiation should be able to fail loudly — "you asked for R=10⁵ and my
   opacity tables stop at 10³" — is an open question (§15.3).

## 14. What this contract hands to the specs downstream

- **W1.6 (Likelihood / NoiseModel)** — the predicted container arriving at a
  likelihood is an ordinary `FunctionSamples` whose mask has already been
  propagated by this contract's rule (§6). The likelihood's own mask handling
  therefore combines *two* masks — predicted and observed — and should say
  which combination it uses (the union is the only defensible one). Note also
  that a calibration nuisance parameter and a noise-model hyperparameter are
  both ordinary parameters on a `Parameterised`, merged by the same mechanism:
  if a real instrument couples them physically (`prior_art.md` 3M1's warning),
  that coupling has to be expressible as a tie or a `HierarchicalPrior`, not by
  merging the two objects.
- **W1.7 (Dataset / FittingProblem)** — a `Dataset` pairs an observed container
  with an `Instrument`. *(Superseded in part, 2026-09-02: this bullet's
  original "the instrument's `label` is the component label to merge its
  `parameters` under" did not survive W1.7's ratified topology —
  `inference.md` §4's nested design merges the instrument's own mapping under
  the reserved role name `instrument` inside the **dataset's** component, so
  an instrument label is provenance and never a merge component, and the
  label-collision worry in §15.4 collapses to one check on dataset labels.)*
  Nesting works because each level re-distributes, and ties that cross levels
  collapse correctly — `inference.md` §6 demonstrates both, answering the
  question this bullet left open. W1.7 also owns *when* `negotiate` and
  `compile_for` are called: once each, jointly across every instrument of a
  model, at `FittingProblem` construction (`inference.md` §8).
- **W1.9 (Lowering)** — a `Transformation` lowers as its parameters and buffers
  do; `apply` is ordinary array code in the backend's array type. The
  declaration forms needing a lowering row are: nothing new for parameters, but
  the chain fold itself must be traceable (no Python-level control flow that
  depends on parameter *values*), and `ACCEPTS`/`PRODUCES`/`label` lower to
  nothing, being composition-time metadata. Negotiation happens entirely before
  tracing.
- **W1.10 (Conformance suite)** — candidate rows already implemented and
  testable here: chain composition rejects a kind mismatch; `require` is used
  for binding; a dropped mask raises; `propagate_mask` matches the ANY rule on
  a random influence matrix; a union of requirements satisfies each of its
  inputs; `compile_for`'s default is the identity; a compiled model's templates
  share axes across evaluations.
- **W1.11 (Modality sketches)** — the stress test `prior_art.md` 3M1 asks for
  belongs here: an instrument whose calibration uncertainty and noise
  covariance are physically coupled, to check that the §4.3/§4.4 split does not
  force such users back into a monolithic plugin. §10's table is the claim to
  check per modality — in particular whether X-ray RMF/ARF really is a matrix
  multiply on an energy-axis `Spectrum`, and whether interferometry's
  `Image` → `VisibilitySet` step can publish anything useful given
  limitation 13.1.
- **W1.12 (Diagnostics)** — a per-sample residual diagnostic sees the predicted
  container's mask, not the model's; the propagation rule (§6) is what relates
  them, and a diagnostic that localises misspecification should not report a
  conservatively masked bin as a data gap.

## 15. Open questions for review

**Ruled by Peter, 2026-09-02**: question 1 — the ANY rule stands as the
default; `min_valid_fraction` remains the named, per-transformation opt-in
extension, to be added only if the conservative rule bites a real pipeline
(none of the W1.11 sketches found it doing so). Question 7 — confirmed:
`Model` belongs to this contract; W1.7 composes it. Question 8 — ruled in
favour of demoting the scalars: `max_step`/`min_resolving_power` are
constructor arguments that do not survive as attributes, and `segments()` is
the only public statement of density (implemented the same day; §7 and the
class docstring updated). Questions 2 and 3 carry recommendations from the
W1.11 interferometry sketch (`docs/design/modalities/interferometry.md`
gaps I-3 and I-4: a chain-internal `configure_from` rather than `pull_back`;
`compile_for` may raise) and await Peter's ruling; questions 4–6 remain open
as written, with question 4's concrete failure mode now recorded in
`docs/design/modalities/spectrum_photometry.md` gap 2 and routed to W1.7.

1. **Is the ANY mask rule too conservative for real resampling?** §6 ratifies
   it, and limitation 13.4 names the extension. The case against: a
   ground-based spectrum with a handful of bad pixels, resampled onto a much
   coarser grid, could lose several output bins entirely. The case for: the
   alternative silently averages an invalid number into a good one. A
   `min_valid_fraction` on the *transformation* (never on `propagate_mask`'s
   default) would settle it if the conservative rule bites in practice.
2. **Should requirements be pulled back through the chain?** Limitation 13.1.
   The current answer says the standard library does not need it. W1.11's
   modality sketches are the place that will find out.
3. **Should a model be able to *refuse* a requirement?** Limitation 13.9. A
   model that cannot reach the requested resolution currently ignores the
   request silently, and the instrument gets a grid it cannot use — a failure
   that surfaces as a bad fit rather than an error. A `compile_for` that may
   raise, or a returned report of what it could and could not honour, would be
   loud; both add API surface to a feature that is meant to be optional.
4. **`Instrument.label` defaults to the channel name.** Convenient, and right
   when one instrument reads one channel. Two instruments on one channel (the
   §7 example) must both be given labels or they collide when merged into a
   joint problem — and nothing currently checks that at the `negotiate` step,
   because `negotiate` does not merge parameters. Should it?
5. **Should `Instrument` cache its merged parameters?** §12 argues no. If W1.7
   finds the repeated merge measurable inside a sampler loop, the answer is a
   `freeze()` that snapshots and refuses later mutation, not a silent cache.
6. **Chains are linear** (limitation 13.7). An instrument that needs two model
   channels — a spectrum plus a continuum-only variant, say, for a line-to-
   continuum ratio — is not expressible. Is that a real case in the target
   scope?
7. **`Model` is defined here, not in W1.7** (§9). Confirm the split: this
   contract owns "produces a `ModelResult`, may be compiled"; W1.7 owns
   everything an engine calls.
8. **`max_step` and `min_resolving_power` mean two things.** On a freshly
   declared requirement they *are* the density; after a union they are a
   strictest-anywhere summary, and `segments()` is what the grid is built from.
   That dual reading is convenient — the scalar is how an instrument author
   naturally writes a requirement — but a reader who trusts `max_step` on a
   merged requirement will overestimate what was asked for. The alternative is
   to make `segments()` the only public form and demote the scalars to
   constructor arguments that do not survive as attributes. Cheap to change
   now, awkward after the freeze.
