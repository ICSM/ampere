# Ampere v2 — ModelResult Schema Contract (W1.4)

Status: **DRAFT for Peter's review.** Implements `DEVELOPMENT_PLAN.md` §4.2 and
the functional-data stance of `architecture.md` §7. Code:
`ampere/core/results_schema.py`, `ampere/core/exceptions.py`. Tests:
`tests/core/test_results_schema.py`, `tests/core/test_spec_doctests.py`.

Every worked example below is executed as a doctest by
`tests/core/test_spec_doctests.py`, so this document cannot drift from the
implementation without the suite going red. Examples share one namespace and
build on each other; read them in order.

Where this document and `DEVELOPMENT_PLAN.md` disagree, the plan wins and this
document is wrong — file it as a decision-log correction.

---

## 1. What this contract is for

Legacy ampere has two conventions for "what the model produced": an attribute
(`self.model.modelFlux`) and a dictionary entry (`result['spectrum']`). Code
that reaches for the wrong one fails late, confusingly, and repeatedly —
`DEVELOPMENT_PLAN.md` §7 lists this first among the traps, and §4.2 exists to
end it. There is now exactly one answer to the question, and it is a
`ModelResult`.

The contract has two halves, and the second is the interesting one:

- **A `ModelResult` is a mapping of named channels to typed containers.** Not a
  single array, not one container per kind: a *name* per output, because one
  model legitimately produces several results of the same kind. Instruments
  bind by name and check the kind (§4.3 of the plan); mismatches fail loudly.
- **A container is a set of samples of a function.** Every observable ampere
  handles *is* a function — flux(λ), V(u,v), position(t) — and a container
  holds explicit coordinates, values, optional uncertainties and a mask. There
  is no regular-grid assumption anywhere, because spectra, light curves and
  (u,v) coverage are all irregular in general.

It is **backend-neutral**. Nothing here imports torch, jax, numpyro or paramax,
lazily or otherwise (`architecture.md` §3–4); it is numpy, `astropy.units` and
stdlib. Turning a container into a backend array is W1.9's lowering spec and
Phase 2's code.

### Setup for the examples

```pycon
>>> import numpy as np
>>> import astropy.units as u
>>> from ampere.core import (
...     DEFAULT_CHANNEL, Axis, AxisSpec, Cube, FunctionSamples, Image, Layout,
...     ModelResult, Order, PhotometricPoints, Spectrum, TimeSeries,
...     VisibilitySet,
... )
>>> from ampere.core.exceptions import ChannelError, SchemaError

```

## 2. The objects

| Object | Role |
|---|---|
| `ModelResult` | Immutable mapping of channel name to container; the single answer to "what did the model produce" |
| `FunctionSamples` | Base container: coordinates + values + optional uncertainty + mask + unit + fidelity + metadata |
| `Spectrum` | Flux at strictly increasing spectral coordinates (wavelength, frequency or energy) |
| `PhotometricPoints` | One value per uniquely named filter; coordinate order explicitly arbitrary |
| `Image` | A 2D map on separable, strictly monotonic spatial axes |
| `Cube` | Two spatial axes plus a spectral axis, separable |
| `TimeSeries` | A quantity at strictly increasing times |
| `VisibilitySet` | **Complex** visibilities at scattered (u,v) points |
| `Axis` | One coordinate axis: values, unit, and *advertised* structure (`regular`, `log_regular`) |
| `AxisSpec`, `Layout`, `Order` | How a container kind declares its axis signature, its layout and its ordering rule |
| `SchemaError`, `ChannelError` | This contract's errors; `ChannelError` is the binding failure a consumer may want to catch on its own |

**Channel names and parameter names are different namespaces.** This is the
obligation `parameters.md` §13 places on this document, and the answer is: they
may collide, harmlessly and deliberately. A model with a parameter `temperature`
may perfectly well emit a channel `temperature` (a temperature map), and nothing
disambiguates them because nothing needs to — parameter names are keyword
arguments *into* a model evaluation, channel names are keys *out of* one, and no
structure in ampere holds both in one flat namespace. Provenance records (W1.8)
must therefore keep them in separate groups rather than merging them into one
attribute dictionary.

Both namespaces do require Python identifiers, for overlapping but not identical
reasons: parameters because their names become keyword arguments, channels
because their names become group names on serialisation and coordinate labels in
`ampere.results`.

```pycon
>>> ModelResult({"2 channels": Spectrum([1.0, 2.0] * u.um, [1.0, 1.0])})
Traceback (most recent call last):
    ...
ampere.core.exceptions.SchemaError: channel name '2 channels' is not usable: a channel name must be a valid Python identifier, because channel names become group names on serialisation and coordinate labels in ampere.results. Rename the channel, e.g. 'channel_2_channels'.

```

## 3. `ModelResult`: named channels, and the trivial case

A `ModelResult` maps names to containers, and implements
`collections.abc.Mapping` — `len`, iteration, `in`, `keys`/`values`/`items` and
`.get` all behave as expected.

```pycon
>>> spectrum = Spectrum([1.0, 2.0, 4.0] * u.um, [3.0, 2.5, 1.0] * u.Jy)
>>> result = ModelResult({"sed": spectrum})
>>> len(result), "sed" in result, list(result)
(1, True, ['sed'])

```

**The single-channel case stays trivial**, as §4.2 requires: a model returning
one bare container gets a default channel name automatically.

```pycon
>>> bare = ModelResult(spectrum)
>>> list(bare), bare.is_single
(['default'], True)
>>> DEFAULT_CHANNEL
'default'
>>> bare.single() is spectrum
True

```

`single()` is sugar for the one-channel case only. On a multi-channel result it
raises rather than picking the first, because silently guessing which output was
meant is precisely the bug class this contract exists to remove.

```pycon
>>> two = ModelResult({"sed": spectrum, "line": spectrum})
>>> two.single()
Traceback (most recent call last):
    ...
ampere.core.exceptions.ChannelError: single() is only meaningful for a one-channel result, but this one has 2 channels: ['sed', 'line'].
Ask for the one you mean by name.

```

### Binding: `require`, and loud mismatches

An instrument binds to a channel **by name** and checks the **kind** at
composition time (`DEVELOPMENT_PLAN.md` §4.2, §4.3). `require` is that entry
point. Both failure modes raise `ChannelError`, so a consumer that wants to fall
back to an alternative channel can catch one type.

```pycon
>>> result.require("sed", Spectrum) is spectrum
True

```

A missing channel names what *is* available, rather than raising a bare
`KeyError` whose message is only the name that was absent:

```pycon
>>> result.require("photometry")
Traceback (most recent call last):
    ...
ampere.core.exceptions.ChannelError: no channel named 'photometry'. Available channels: 'sed' (Spectrum).

```

A kind mismatch says what was found, what was wanted, and both ways to fix it:

```pycon
>>> result.require("sed", Image)
Traceback (most recent call last):
    ...
ampere.core.exceptions.ChannelError: channel 'sed' holds a Spectrum, but an Image was required. Either bind to a channel of the right kind (there are none in this result), or change the model to produce an Image on 'sed'.

```

`ChannelError` is a `KeyError` as well as a `SchemaError`, following
`exceptions.py`'s rule that a contract error stays catchable by the builtin a
caller would naturally reach for. That is load-bearing rather than decorative:
`Mapping`'s `get` and `__contains__` mixins are defined in terms of catching
`KeyError`, so a result would otherwise raise where it should answer.

```pycon
>>> result.get("photometry") is None
True
>>> "photometry" in result, "sed" in result
(False, True)
>>> isinstance(ChannelError("x"), KeyError)
True

```

And asking a default-only result for a named channel says *why* the name is
missing, which is almost always "the model returned a bare container":

```pycon
>>> bare.require("sed")
Traceback (most recent call last):
    ...
ampere.core.exceptions.ChannelError: no channel named 'sed'. Available channels: 'default' (Spectrum). This result has only the automatic 'default' channel, which means the model returned a bare container; name its channels explicitly if an instrument needs to bind to one by name.

```

### Immutable updates

A transformation chain rewrites one channel and passes the rest through, so
updates return a new result rather than mutating:

```pycon
>>> scaled = result.with_channels(sed=spectrum.with_values([6.0, 5.0, 2.0]))
>>> list(scaled), float(scaled["sed"].values[0])
(['sed'], 6.0)
>>> float(result["sed"].values[0])       # the original is untouched
3.0

```

## 4. Containers are coordinate-indexed function samples

Every container holds explicit coordinates. A `Spectrum` is **not** two
positionally aligned arrays with implicit ordering assumptions baked into
whatever reads them next — that is the legacy idiom `architecture.md` §7 names
as the thing to kill.

```pycon
>>> spectrum.spectral_axis
<Axis 'spectral_axis' n=3 um log-regular>
>>> spectrum.spectral_axis.values.tolist()
[1.0, 2.0, 4.0]
>>> spectrum.flux.tolist()
[3.0, 2.5, 1.0]
>>> spectrum.unit
Unit("Jy")

```

Coordinates and values are aligned index by index, and a mismatch is caught
here rather than becoming an unrelated broadcasting error inside a likelihood:

```pycon
>>> Spectrum([1.0, 2.0, 4.0] * u.um, [3.0, 2.5])
Traceback (most recent call last):
    ...
ampere.core.exceptions.SchemaError: Spectrum's values have shape (2,), but its coordinates imply (3,) (spectral_axis=3). Coordinates and values are aligned index-by-index; a mismatch here is a bug that would otherwise surface as an unrelated broadcasting error inside a likelihood.

```

Containers are immutable, and so are the arrays inside them: a consumer cannot
accidentally corrupt a model's output in place.

```pycon
>>> spectrum.values[0] = 99.0
Traceback (most recent call last):
    ...
ValueError: assignment destination is read-only

```

### Two layouts, neither of which is a grid assumption

`Layout.POINTS` — every axis has length N and `values.shape == (N,)`: an
arbitrary point set in as many dimensions as there are axes. `Layout.GRID` —
axes are separable, so `values.shape` is the tuple of axis lengths.

```pycon
>>> Spectrum.LAYOUT, Image.LAYOUT
(<Layout.POINTS: 'points'>, <Layout.GRID: 'grid'>)

```

Separable is **not** the same claim as regular: a grid layout says the axes
factorise, not that either is evenly spaced. Nothing in this contract requires
even spacing anywhere.

## 5. Ordering: validated or explicitly tolerated, never assumed

`architecture.md` §7 requires each container either to *validate* its ordering
assumptions or to *document* that it tolerates arbitrary order. Every kind here
does one or the other, and which one is part of the kind's definition.

| Kind | Axis | Rule | Why |
|---|---|---|---|
| `Spectrum` | `spectral_axis` | strictly increasing | quasiseparable GP solvers need ordered 1D coordinates; duplicates make a covariance singular |
| `TimeSeries` | `time` | strictly increasing | same |
| `Image`, `Cube` | `x`, `y` | strictly monotonic (either direction) | sky axes legitimately run either way |
| `Cube` | `spectral_axis` | strictly increasing | as `Spectrum` |
| `PhotometricPoints` | `spectral_axis` | **any** | a point is identified by its filter, not its position; catalogues arrive in archive order |
| `VisibilitySet` | `u`, `v` | **any** | a (u,v) point set has no natural order; imposing one would be a fiction |

Violating a validated rule fails at construction with a message that names both
the problem and the fix:

```pycon
>>> Spectrum([2.0, 1.0, 3.0] * u.um, [1.0, 1.0, 1.0])
Traceback (most recent call last):
    ...
ampere.core.exceptions.SchemaError: Spectrum requires its 'spectral_axis' coordinates to be strictly increasing, but its coordinates are not sorted. This is checked here rather than assumed, because quasiseparable GP solvers and resampling both silently misbehave on unordered or duplicated coordinates. Use Spectrum.from_unsorted(...) to sort them, or split genuinely overlapping data into separate channels.

```

Duplicates are reported as such, because the fix differs — sorting will not help:

```pycon
>>> Spectrum([1.0, 1.0, 3.0] * u.um, [1.0, 1.0, 1.0])
Traceback (most recent call last):
    ...
ampere.core.exceptions.SchemaError: Spectrum requires its 'spectral_axis' coordinates to be strictly increasing, but it contains 1 repeated coordinate(s). ...

```

`from_unsorted` is the explicit opt-in for data that genuinely arrive in
arbitrary order. It sorts, carrying uncertainties and mask along with the
permutation; duplicates still raise.

```pycon
>>> tidy = Spectrum.from_unsorted(
...     [3.0, 1.0, 2.0] * u.um, [30.0, 10.0, 20.0] * u.Jy, uncertainty=[3.0, 1.0, 2.0] * u.Jy
... )
>>> tidy.spectral_axis.values.tolist(), tidy.flux.tolist(), tidy.uncertainty.tolist()
([1.0, 2.0, 3.0], [10.0, 20.0, 30.0], [1.0, 2.0, 3.0])

```

### Regularity is advertised, never required

An `Axis` reports whether it happens to be evenly spaced, so an implementation
may choose an FFT convolution or a Toeplitz solve when it legitimately can.
Nothing may *require* it.

```pycon
>>> linear = Spectrum(np.linspace(1.0, 10.0, 10) * u.um, np.ones(10))
>>> linear.spectral_axis.regular, round(linear.spectral_axis.step, 3)
(True, 1.0)

```

Constant resolving power λ/Δλ — the case in which a constant-velocity LSF
becomes a convolution — is advertised separately, because it is a different fast
path:

```pycon
>>> logarithmic = Spectrum(np.geomspace(1.0, 100.0, 20) * u.um, np.ones(20))
>>> logarithmic.spectral_axis.regular, logarithmic.spectral_axis.log_regular
(False, True)

```

And an irregular axis simply advertises nothing, which is a perfectly ordinary
state of affairs:

```pycon
>>> irregular = Spectrum([1.0, 1.5, 9.0] * u.um, np.ones(3))
>>> irregular.spectral_axis.regular, irregular.spectral_axis.log_regular
(False, False)
>>> irregular.spectral_axis.step is None
True

```

## 6. Units: converted once, at construction

A container stores plain numpy arrays plus a unit. No `Quantity` ever reaches
the hot loop — the units trap in `DEVELOPMENT_PLAN.md` §7.

```pycon
>>> type(spectrum.values).__name__, spectrum.unit
('ndarray', Unit("Jy"))

```

Passing a `Quantity` adopts its unit; passing a `Quantity` *and* an explicit
`unit=` converts once, there and then:

```pycon
>>> converted = Spectrum([1.0, 2.0] * u.um, [1000.0, 2000.0] * u.mJy, unit=u.Jy)
>>> converted.values.tolist(), converted.unit
([1.0, 2.0], Unit("Jy"))

```

`to_unit` is the composition-time conversion, and it carries uncertainties:

```pycon
>>> in_mjy = spectrum.to_unit(u.mJy)
>>> in_mjy.values.tolist(), in_mjy.unit
([3000.0, 2500.0, 1000.0], Unit("mJy"))

```

An axis whose physical type matters says so, and rejects the wrong one — this
catches transposed arguments, which are otherwise silent:

```pycon
>>> Spectrum(np.arange(3.0) * u.s, np.ones(3))
Traceback (most recent call last):
    ...
ampere.core.exceptions.SchemaError: Spectrum's 'spectral_axis' coordinates are in s (physical type 'time'), which this axis does not accept. It accepts ['length', 'frequency', 'energy']. Check you have not passed the axes in the wrong order.

```

A bare array where a unit is needed is a loud failure, not a guess:

```pycon
>>> Spectrum(np.arange(3.0), np.ones(3))
Traceback (most recent call last):
    ...
ampere.core.exceptions.SchemaError: Spectrum's 'spectral_axis' coordinates need a unit, but a bare array was given. Pass an astropy Quantity of length or frequency or energy physical type, e.g. spectral_axis=values * u.um. Units are converted once here and never in the hot loop.

```

*Value* units are deliberately unconstrained: a model may emit F_λ, F_ν, a
surface brightness or a dimensionless ratio, and this contract has no business
ruling on which. Only *coordinate* axes carry a physical-type requirement,
because there the meaning is unambiguous.

## 7. Masks: first-class, and exactly zero information

`mask` is a boolean array in which **`True` marks an excluded sample**,
following the numpy and `astropy.nddata.NDData` convention. A masked sample
carries **zero information**: it is to be treated exactly as if it had not been
observed.

```pycon
>>> masked = Spectrum(
...     [1.0, 2.0, 4.0] * u.um,
...     [3.0, 2.5, 1.0] * u.Jy,
...     uncertainty=[0.1, 0.2, 0.4] * u.Jy,
...     mask=np.array([False, True, False]),
... )
>>> masked.is_masked, masked.n_samples, masked.n_valid
(True, 3, 2)
>>> masked.valid.tolist()
[True, False, True]

```

`prior_art.md` lesson R1 records that RHMF expresses the same idea as either
zero weight or infinite uncertainty, and notes that ampere's mask composes with
that convention without impedance mismatch. Both representations are available,
so a consumer uses whichever suits its algebra:

```pycon
>>> masked.weights().tolist()
[1.0, 0.0, 1.0]
>>> masked.masked_uncertainty().tolist()
[0.1, inf, 0.4]

```

Masked values stay in place in `values`. Ampere never threads sentinel NaNs
through arithmetic — R1's explicit point — so a masked sample is still a
perfectly ordinary number that simply carries no weight:

```pycon
>>> masked.values.tolist()
[3.0, 2.5, 1.0]

```

A mask must be boolean. An array of *good*-data flags is the complement, and
saying so is the caller's job, not a guess this contract makes:

```pycon
>>> Spectrum([1.0, 2.0] * u.um, [1.0, 1.0], mask=np.array([1, 0]))
Traceback (most recent call last):
    ...
ampere.core.exceptions.SchemaError: Spectrum's mask must be a boolean array, got dtype int64. True marks an *excluded* sample, following the numpy and astropy NDData convention; if you have an array of good-data flags, pass ~good.

```

## 8. Masking is not censoring — and the hook that stays open

An upper limit is **not** a masked point. A non-detection at 3σ says something
quite definite about the source; discarding it as "no information" throws that
away, and encoding it as a data point with a large error bar is simply the wrong
likelihood. `DEVELOPMENT_PLAN.md` §4.2 keeps the two distinct and puts censoring
in the likelihood contract (W1.6, issue #11).

This contract therefore does **not** add a `limit` flag, and deliberately does
not overload `mask` to mean one. What it does is leave the hook open, in two
concrete ways:

1. **`mask` is strictly binary include/exclude, and documented as such.** It has
   no third state and no reserved values, so W1.6 may introduce a censoring
   classification alongside it without a semantic collision or a compatibility
   break.
2. **Per-sample auxiliary arrays are a supported, aligned concept**
   (`extra_coords`). They are validated against the container's shape exactly as
   values and mask are, which is precisely the machinery a per-point
   `limit_kind` array needs. W1.6 can therefore express censoring either as a
   `Dataset`-level array or as an extra coordinate, and this contract's shape
   rules already cover it.

```pycon
>>> photometry = PhotometricPoints(
...     ["WISE_W3", "WISE_W4"],
...     [12.1, 22.2] * u.um,
...     [0.4, 0.9] * u.Jy,
...     uncertainty=[0.05, 0.30] * u.Jy,
...     extra_coords={"detected": np.array([True, False])},
... )
>>> photometry.extra_coords["detected"].tolist()
[True, False]
>>> photometry.is_masked          # a non-detection is NOT masked
False

```

The `detected` array above is an illustration of the shape of the hook, not a
ratified part of the censoring interface — naming and semantics are W1.6's to
fix.

## 9. Fidelity tags

A channel may carry a fidelity tag naming which variant of a model produced it.
This is the hook the plan's design horizon (a) reserves for multi-fidelity
modelling: cheap and expensive variants sharing one parameter space, with
inference engines that mix fidelities slotting in above the backends later.

Ampere attaches **no meaning** to the string. It is a label for an engine that
understands the model's own fidelity ladder; nothing in the core branches on it.

```pycon
>>> cheap = Spectrum([1.0, 2.0] * u.um, [1.0, 0.9] * u.Jy, fidelity="lte")
>>> expensive = Spectrum([1.0, 2.0] * u.um, [1.1, 0.8] * u.Jy, fidelity="full_nlte")
>>> ladder = ModelResult({"sed_lte": cheap, "sed_nlte": expensive})
>>> dict(ladder.fidelities())
{'sed_lte': 'lte', 'sed_nlte': 'full_nlte'}

```

## 10. The hot loop: `with_values`

Full construction validates coordinates, units and ordering in O(N). That is the
right price to pay **once**, at composition time, and the wrong price to pay on
every likelihood evaluation.

`with_values` is the hot-loop constructor: it reuses the already-validated axes
— the same `Axis` objects, with their unit checks and advertised structure
already computed — and only swaps the values.

```pycon
>>> template = Spectrum(np.geomspace(1.0, 100.0, 50) * u.um, np.zeros(50), unit=u.Jy)
>>> evaluated = template.with_values(np.ones(50))
>>> evaluated.axes is template.axes          # axes shared, not revalidated
True
>>> evaluated.unit, evaluated.shape
(Unit("Jy"), (50,))

```

This is the point at which §4.2 meets §4.3: requirements negotiation exists so
that a model can build its per-channel evaluation grids **once**, from the union
of its instruments' requirements, and then evaluate onto them repeatedly.
`with_values` is what "then evaluate onto them" means in this contract.

The template's unit is authoritative. A `Quantity` in a convertible unit is
*converted*, not reinterpreted — getting this wrong would be the units trap in
its most damaging form, a factor of 1000 that never announces itself:

```pycon
>>> jansky = Spectrum([1.0, 2.0] * u.um, [1.0, 1.0] * u.Jy)
>>> jansky.with_values([2000.0, 3000.0] * u.mJy).values.tolist()
[2.0, 3.0]

```

Changing the unit is refused outright, because the unit belongs to the composed
container and was fixed once, at construction:

```pycon
>>> Spectrum([1.0, 2.0] * u.um, [1.0, 1.0]).with_values([1.0, 2.0] * u.Jy)
Traceback (most recent call last):
    ...
ampere.core.exceptions.SchemaError: with_values() cannot change Spectrum's value unit: the template is in no unit but the new values are in Jy. ...

```

Shape is still checked, because a model that changes its output length between
evaluations is a bug worth catching immediately:

```pycon
>>> template.with_values(np.ones(49))
Traceback (most recent call last):
    ...
ampere.core.exceptions.SchemaError: Spectrum's values have shape (49,), but its coordinates imply (50,) (spectral_axis=50). ...

```

## 11. Worked example: a low-resolution SED plus high-resolution CO windows

This is the case `DEVELOPMENT_PLAN.md` §4.2 uses to justify named channels, and
it is the acceptance criterion for this item. A dusty circumstellar envelope
model produces a broad, coarsely sampled SED *and* a handful of narrow, finely
sampled windows around CO rotational lines. Both are `Spectrum`. Neither is a
resampling of the other. Evaluating the fine grid everywhere would be
prohibitive, and evaluating the coarse grid at the lines would be useless.

The SED channel: broad coverage, logarithmically spaced, low resolution.

```pycon
>>> sed_wavelengths = np.geomspace(1.0, 2000.0, 60) * u.um
>>> sed = Spectrum(
...     sed_wavelengths,
...     np.geomspace(50.0, 0.2, 60) * u.Jy,
...     fidelity="continuum",
...     meta={"note": "dust continuum only"},
... )
>>> sed.n_samples, sed.spectral_axis.log_regular
(60, True)

```

The CO channel: two narrow windows around the J=3–2 (866.96 µm) and J=2–1
(1300.40 µm) lines, finely sampled *within* each window and separated by a gap
of several hundred microns. This is one channel, and it is emphatically not on
any grid.

```pycon
>>> window_32 = np.linspace(866.86, 867.06, 21)
>>> window_21 = np.linspace(1300.30, 1300.50, 21)
>>> co_wavelengths = np.concatenate([window_32, window_21]) * u.um
>>> co_windows = Spectrum(
...     co_wavelengths,
...     np.concatenate([np.linspace(1.0, 1.4, 21), np.linspace(0.8, 1.1, 21)]) * u.Jy,
...     fidelity="co_lte",
... )
>>> co_windows.n_samples
42

```

The channel is strictly increasing — so a quasiseparable solver can consume it —
while being nothing like evenly spaced, because of the gap between the windows:

```pycon
>>> co_windows.spectral_axis.regular, co_windows.spectral_axis.log_regular
(False, False)
>>> float(np.diff(co_windows.spectral_axis.values).max()) > 400.0
True

```

Both channels go into one result, under names an instrument can bind to:

```pycon
>>> result = ModelResult(
...     {"sed_lowres": sed, "co_windows": co_windows},
...     meta={"model": "dusty_envelope"},
... )
>>> result
<ModelResult sed_lowres: Spectrum@continuum[60], co_windows: Spectrum@co_lte[42]>
>>> sorted(result)
['co_windows', 'sed_lowres']

```

Two channels of the **same kind** is the whole point: a kind-only schema could
not tell them apart, so the schema is name-first.

```pycon
>>> {name: kind.__name__ for name, kind in result.kinds().items()}
{'sed_lowres': 'Spectrum', 'co_windows': 'Spectrum'}
>>> sorted(result.of_kind(Spectrum))
['co_windows', 'sed_lowres']

```

A photometric instrument binds to the SED; a heterodyne receiver binds to the CO
windows. Each says which by name, and each gets a kind check for free:

```pycon
>>> photometer_input = result.require("sed_lowres", Spectrum)
>>> heterodyne_input = result.require("co_windows", Spectrum)
>>> photometer_input.n_samples, heterodyne_input.n_samples
(60, 42)

```

Binding to a channel that does not exist fails loudly at composition time,
naming what is there — the failure mode §4.2 asks for:

```pycon
>>> result.require("co_lines")
Traceback (most recent call last):
    ...
ampere.core.exceptions.ChannelError: no channel named 'co_lines'. Available channels: 'sed_lowres' (Spectrum), 'co_windows' (Spectrum).

```

The channels carry independent fidelity tags, so a multi-fidelity engine can
later mix them without the model changing:

```pycon
>>> dict(result.fidelities())
{'sed_lowres': 'continuum', 'co_windows': 'co_lte'}
>>> result.meta["model"]
'dusty_envelope'

```

And in the hot loop, each channel is refilled from its own template, with no
coordinate revalidation and no unit arithmetic:

```pycon
>>> next_sample = result.with_channels(
...     sed_lowres=sed.with_values(np.geomspace(55.0, 0.18, 60)),
...     co_windows=co_windows.with_values(np.ones(42)),
... )
>>> next_sample["sed_lowres"].axes is sed.axes
True
>>> round(float(next_sample["sed_lowres"].values[0]), 1)
55.0

```

## 12. The other kinds

`PhotometricPoints` — one value per uniquely named filter. The filter name is
the *identity* of a point, so duplicates are refused; coordinate order is
explicitly arbitrary.

```pycon
>>> phot = PhotometricPoints(
...     ["2MASS_Ks", "WISE_W1", "IRAS_12"],
...     [2.16, 3.35, 12.0] * u.um,
...     [12.0, 8.0, 3.0] * u.Jy,
...     uncertainty=[0.5, 0.4, 0.3] * u.Jy,
... )
>>> phot.filters.tolist()
['2MASS_Ks', 'WISE_W1', 'IRAS_12']
>>> PhotometricPoints(["W1", "W1"], [3.4, 3.4] * u.um, [1.0, 1.0] * u.Jy)
Traceback (most recent call last):
    ...
ampere.core.exceptions.SchemaError: PhotometricPoints was given repeated filter names ['W1']. A filter name is the identity of a photometric point, so duplicates make the channel ambiguous. Put repeat observations of one filter in separate channels, or disambiguate the names (e.g. 'WISE_W1_epoch1').

```

`TimeSeries` — strictly increasing times, arbitrary cadence.

```pycon
>>> lightcurve = TimeSeries([0.0, 1.0, 1.5, 9.0] * u.day, [1.0, 1.2, 1.1, 0.9] * u.Jy)
>>> lightcurve.time.regular, lightcurve.n_samples
(False, 4)

```

`Image` — separable, strictly monotonic spatial axes; note the y axis running
*downwards*, which is legitimate and accepted.

```pycon
>>> sky = Image(
...     np.linspace(-2.0, 2.0, 5) * u.arcsec,
...     np.linspace(2.0, -2.0, 5) * u.arcsec,
...     np.ones((5, 5)) * u.Jy,
... )
>>> sky.shape, sky.x.regular
((5, 5), True)

```

`Cube` — two spatial axes and a spectral axis, in that order.

```pycon
>>> cube = Cube(
...     np.linspace(-1.0, 1.0, 3) * u.arcsec,
...     np.linspace(-1.0, 1.0, 4) * u.arcsec,
...     np.linspace(866.9, 867.0, 6) * u.um,
...     np.ones((3, 4, 6)) * u.Jy,
... )
>>> cube.shape
(3, 4, 6)

```

### `VisibilitySet`: the kind that stresses the schema hardest

Complex values, an irregular point set no grid could describe, and no natural
ordering. It is the Phase 4 proof modality, and it is an ordinary container.

```pycon
>>> vis = VisibilitySet(
...     [120.0, -35.0, 88.0, -210.0],
...     [45.0, 190.0, -66.0, 12.0],
...     [1.0 + 0.2j, 0.6 - 0.3j, 0.4 + 0.0j, 0.1 - 0.05j],
...     uncertainty=[0.02, 0.02, 0.03, 0.05],
...     extra_coords={"frequency_ghz": np.array([230.5, 230.5, 345.8, 345.8])},
... )
>>> vis.n_samples, vis.values.dtype.kind
(4, 'c')
>>> np.round(vis.amplitude(), 3).tolist()
[1.02, 0.671, 0.4, 0.112]

```

The (u,v) coordinates are dimensionless by convention — baselines measured in
wavelengths — so a bare array is accepted without ceremony; spatial frequency in
rad⁻¹ is accepted too:

```pycon
>>> vis.u.unit
Unit(dimensionless)
>>> VisibilitySet([1.0, 2.0] / u.rad, [3.0, 4.0] / u.rad, [1.0 + 0j, 0.5 + 0j]).v.unit
Unit("1 / rad")

```

Uncertainty stays **real** on complex data: it is the per-component standard
deviation of a circular complex Gaussian, the standard interferometric noise
model. Non-circular noise, and the amplitude/phase (Rice, von Mises)
formulations, are noise-model concerns and belong to W1.6 — the container
deliberately does not encode them.

Complex values in a kind that does not admit them is an error, not a silent cast
to the real part:

```pycon
>>> Spectrum([1.0, 2.0] * u.um, [1.0 + 1j, 2.0 + 0j])
Traceback (most recent call last):
    ...
ampere.core.exceptions.SchemaError: Spectrum values include complex numbers, but this container kind holds real values. Complex data belong in a VisibilitySet (or a container kind that declares ALLOW_COMPLEX = True); to keep a real projection of complex data, take its amplitude, phase, real or imaginary part explicitly.

```

## 13. User-defined kinds

`DEVELOPMENT_PLAN.md` §4.3 makes user extensibility a first-class requirement: a
new data container must be definable by subclassing one class, out of tree, with
no changes inside ampere. A kind is three class attributes and (optionally) a
friendlier `__init__`; everything else — validation, masks, units, `with_values`,
equality, `repr` — is inherited.

```pycon
>>> class PolarisationCurve(FunctionSamples):
...     """Degree of polarisation against strictly increasing wavelength."""
...     AXES = (AxisSpec("spectral_axis", physical_types=("length",),
...                      order=Order.STRICTLY_INCREASING),)
...     LAYOUT = Layout.POINTS
...
>>> curve = PolarisationCurve({"spectral_axis": [0.4, 0.6, 0.8] * u.um}, [0.02, 0.03, 0.01])
>>> curve.n_samples, curve.axis("spectral_axis").unit
(3, Unit("um"))
>>> ModelResult({"polarisation": curve}).require("polarisation", PolarisationCurve) is curve
True

```

The inherited validation applies to the new kind unchanged, which is the point:

```pycon
>>> PolarisationCurve({"spectral_axis": [0.6, 0.4] * u.um}, [0.02, 0.03])
Traceback (most recent call last):
    ...
ampere.core.exceptions.SchemaError: PolarisationCurve requires its 'spectral_axis' coordinates to be strictly increasing, but its coordinates are not sorted. ...

```

## 14. Decisions and their reasoning

| Decision | Reasoning |
|---|---|
| A result is a mapping of *named* channels, not one container per kind | One model legitimately produces several outputs of the same kind (the §11 example); a kind-keyed schema cannot express it at all |
| A bare container is filed under `DEFAULT_CHANNEL` | §4.2 requires the one-channel case to stay trivial; the alternative taxes every simple model to serve the complex one |
| `single()` raises on a multi-channel result | Silently returning the first output would reintroduce exactly the guess-which-output ambiguity the contract exists to remove |
| Channel names must be Python identifiers | They become group names on serialisation and coordinate labels in `ampere.results`; identifiers are the portable intersection |
| Channel and parameter names are separate namespaces and may collide | They are keys in opposite directions — arguments in, results out — and nothing holds both in one flat namespace (`parameters.md` §13's question, answered) |
| Binding failures raise `ChannelError`, a `SchemaError` subclass | A consumer may reasonably want to catch "no such channel / wrong kind" and try an alternative, without catching every malformed-container error |
| `ChannelError` is also a `KeyError`, with `__str__` overridden | `exceptions.py`'s rule that contract errors stay catchable by the obvious builtin; and `Mapping.get`/`__contains__` are *defined* in terms of catching `KeyError`, so without it a result would raise where it should answer. `KeyError` alone formats as `repr(args[0])`, which would put quotes round every message |
| Kind checking happens in `require`, at composition time | §4.2 asks for loud mismatches; a shape error inside a likelihood ten frames down is the failure mode being avoided |
| Every container is coordinate + value (+ uncertainty, + mask) | §4.2's functional view; it is also what coordinate-aware SBI embeddings and coordinate-conditioned emulators (design horizon (c)) need to exist |
| Ordering is validated per kind, or documented as tolerated | `architecture.md` §7 forbids the third option — an unstated assumption — which is what legacy had |
| `Spectrum`/`TimeSeries` require strictly increasing coordinates | Quasiseparable GP solvers need ordered 1D coordinates, and duplicates make a covariance singular; both fail in ways that look like science problems |
| `from_unsorted` exists, but sorting is never automatic | Silently reordering a user's arrays is how the alignment between coordinates and values gets quietly broken; opting in is one call |
| `PhotometricPoints`/`VisibilitySet` tolerate arbitrary order | A filter-keyed point and a (u,v) sample have no natural order; requiring one would be a fiction, so it is documented instead |
| Two layouts (`POINTS`, `GRID`), neither implying regularity | Storing an image as 10⁴ scattered coordinates wastes memory for no gain; separability is a real distinction, regularity is not the same claim |
| `regular`/`log_regular` advertised on the axis | `architecture.md` §7: fast paths may be *taken*, never *required*. `log_regular` earns its place because constant-velocity LSF convolution needs exactly it |
| Units converted once at construction; values are plain arrays | The units trap, `DEVELOPMENT_PLAN.md` §7 |
| Physical types constrained on coordinate axes only, never on values | An axis's meaning is unambiguous; a model may legitimately emit F_λ, F_ν, brightness or a ratio, and the schema has no business ruling |
| `mask` is `True` = excluded | numpy and `astropy.nddata.NDData` both do this; inventing the opposite convention would be a permanent papercut |
| Both `weights()` and `masked_uncertainty()` provided | `prior_art.md` R1: zero weight and infinite uncertainty are the same statement, and different consumers' algebra wants different ones |
| Masked values stay in the array; no sentinel NaNs | R1 explicitly; a NaN threaded through arithmetic contaminates everything it touches |
| Censoring is *not* in this contract, but the hook is left open (§8) | §4.2 puts it in the likelihood contract (W1.6, issue #11); overloading `mask` would conflate zero information with real information |
| `uncertainty` stays real for complex values | The circular complex Gaussian is the standard interferometric model; non-circular noise is a `NoiseModel`'s business (W1.6), not a container's |
| `with_values` is the hot-loop constructor | O(N) coordinate validation is right once and wrong per evaluation; it is also where §4.2 meets §4.3's negotiated grids |
| Containers are immutable, with read-only arrays | A consumer must not be able to corrupt a model's output in place, and shared axes (`with_values`) make aliasing routine |
| `extra_coords` for per-sample labels that are not geometry | Filter names, per-visibility frequency, epoch labels; keeping them out of `axes` keeps "what indexes this container" a crisp question |

## 15. Deliberate limitations of v1.4

Each of these is a decision, not an oversight. Each has an extension point.

1. **No scattered `Image`/`Cube`.** Both are `Layout.GRID`, so a mosaic with
   irregular pixel positions is not expressible as one container. The extension
   point is a `ScatteredImage` kind with `Layout.POINTS` and two spatial axes —
   the machinery already supports it, only the kind is absent.
2. **One value array per container.** Stokes I/Q/U/V, or a model's flux and its
   optical depth, are separate *channels*, not extra columns. This keeps
   "coordinates index values" a single unambiguous statement, and channels are
   free.
3. **Only linear and logarithmic regularity are advertised.** A richer
   grid-structure descriptor (the architecture spec's phrasing) would name other
   fast paths; `Axis` gains a field when a consumer needs one.
4. **No unit on `extra_coords`.** They are labels, and a labelled quantity that
   needs a unit is arguably an axis. If W1.6 wants a per-point frequency with
   units, `extra_coords` should become a mapping of `Axis`.
5. **No cross-channel consistency checks.** Nothing verifies that two channels
   claiming to be the same model's output share a unit or overlap sensibly. That
   is composition's job (W1.5) and it needs the instruments to say what they
   want first.
6. **`with_values` cannot change shape or unit.** By design — both are
   properties of the negotiated grid, and changing either mid-loop means the
   template was wrong. A `Quantity` in a *convertible* unit is converted (§10);
   an inconvertible one, or any unit at all on a unitless template, raises.
   Build a new container.
7. **No serialisation.** `ModelResult` has no `to_spec`/`from_spec` yet, unlike
   `ParameterSet`. W1.8 owns results emission and should decide the format
   rather than having one imposed here; design horizon (c)'s "serialisable
   training sets" for emulators is the requirement to design against.
8. **Complex support is per kind, via `ALLOW_COMPLEX`.** Only `VisibilitySet`
   sets it. A complex-valued spectrum (a model amplitude before squaring) would
   need its own kind or a flag change.
9. **Masks do not propagate automatically.** A container knows its own mask;
   nothing here makes a transformation carry it through. That obligation is
   stated for W1.5 below rather than enforced here, because only a
   transformation knows whether its output samples correspond one-to-one with
   its input ones.

## 16. What this contract hands to the specs downstream

- **W1.5 (Transformation / Instrument)** — the binding interface is
  `ModelResult.require(name, kind)`; use it rather than `[]` so kind mismatches
  are loud. Two obligations follow from §7 and §15.9: a transformation **must**
  propagate masks (and say what it does when resampling makes the
  correspondence many-to-one — a resampled bin touching any masked input sample
  is the conservative rule, and W1.5 should ratify or replace it), and
  requirements negotiation should express its output as the *coordinates* a
  channel needs, since a container is built from coordinates first. The
  compile-once/evaluate-many split this contract assumes is
  `with_values`: negotiation produces the template, the hot loop refills it.
- **W1.6 (Likelihood / NoiseModel)** — masks arrive as either `weights()` or
  `masked_uncertainty()`; pick one and state it. Censoring is yours (§8): the
  hook is that `mask` is strictly binary and that per-sample aligned arrays are
  a supported concept. `VisibilitySet`'s real-valued uncertainty encodes the
  circular complex Gaussian only; a non-circular or amplitude/phase noise model
  must supply its own structure. Note also that a container advertises
  `Axis.regular` and `Axis.log_regular` but never requires them — a solver
  strategy may branch on them, and must have a path when both are false.
- **W1.7 (Dataset / FittingProblem)** — observed data are containers of the same
  kinds, so a `Dataset` pairs an observed container with a channel name and an
  instrument. The shape rules here are what makes "predicted and observed are
  aligned index-by-index" checkable at composition time rather than at
  evaluation time.
- **W1.8 (Results)** — channel names are ArviZ group/coordinate names, which is
  why they are identifiers. Keep channel names and parameter names in *separate*
  groups in the provenance attrs (§2): they may legitimately collide.
  Serialisation of containers is unclaimed (§15.7) and should be settled here.
- **W1.9 (Lowering)** — a container lowers to (coordinate arrays, value array,
  uncertainty array, mask) in the backend's array type; `unit` and `fidelity`
  lower to nothing, being composition-time metadata. The `Layout` distinction is
  the one structural thing a backend must respect.
- **W1.10 (Conformance suite)** — candidate rows already implemented and
  testable here: `with_values` preserves axes identity and unit; mask/weights
  round trip; unit conversion is exact; ordering validation fires per kind;
  complex values are confined to `ALLOW_COMPLEX` kinds; a default channel is
  created for a bare container; kind mismatch raises `ChannelError`.
- **W1.11 (Modality sketches)** — `VisibilitySet` is the Phase 4 proof modality
  and is deliberately the hardest case here (complex, scattered, unordered). The
  X-ray modality should check whether `Spectrum` with an energy axis plus a
  response matrix is expressible, since the axis already accepts `energy`.

## 17. Open questions for review

1. **`DEFAULT_CHANNEL` is the string `"default"`.** It is short and obvious, but
   it is also a name a user might plausibly want for a real channel. Reserving
   it (refusing an explicit `"default"`) would be safer and more annoying;
   currently it is not reserved.
2. **Strictly increasing, not merely non-decreasing.** Duplicate spectral
   coordinates are refused outright. Two overlapping échelle orders are the real
   case this affects, and the contract's answer is "two channels". Confirm that
   is acceptable rather than needing a merge helper.
3. **`Image`/`Cube` axis names are `x`/`y`.** Neutral, but a WCS-aware
   alternative (`lon`/`lat`, with a frame) may be wanted once §4.7's astropy
   interop is written. Deferred deliberately; W1.11's imaging sketch should say.
4. **Should `Cube` be `(x, y, spectral)` or `(spectral, x, y)`?** Chosen as
   `(x, y, spectral)` so the spatial axes stay adjacent and the spectral axis is
   the fastest-varying, which suits per-pixel spectral fitting. FITS convention
   would put spectral first. Cheap to change now, expensive later.
5. **`extra_coords` has no unit support** (§15.4). If W1.6 wants per-visibility
   frequency as a first-class quantity, it should become a mapping of `Axis`.
6. **No serialisation yet** (§15.7). W1.8 owns the format, but design horizon
   (c) wants serialisable (θ, `ModelResult`) training sets for emulators; that
   requirement should be nailed down before Phase 2 rather than after.
7. **Should `ModelResult` carry the parameter values that produced it?** It
   would make (θ, result) pairs self-contained for emulator training and for
   provenance, at the cost of coupling this contract to W1.3's. Currently it
   does not; `meta` is free-form and could hold them by convention.
