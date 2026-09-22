# Ampere v2 — Parameter & Prior Contract (W1.3)

Status: **frozen at `spec-v1.0`** (the tag created at the W1.13 merge, 2026-09; any later change to a §4 contract requires a decision-log entry in `DEVELOPMENT_PLAN.md` in the same PR — ground rule 9). Implements `DEVELOPMENT_PLAN.md` §4.1 and
the parameters-vs-buffers half of `architecture.md` §6. Code:
`ampere/core/parameter.py`, `ampere/core/exceptions.py`. Tests:
`tests/core/`.

Every worked example below is executed as a doctest by
`tests/core/test_spec_doctests.py`, so this document cannot drift from the
implementation without the suite going red. Examples share one namespace and
build on each other; read them in order.

Where this document and `DEVELOPMENT_PLAN.md` disagree, the plan wins and this
document is wrong — file it as a decision-log correction.

---

## 1. What this contract is for

Legacy ampere addresses parameters positionally: a model declares `npars`, an
inference driver slices `theta[i:j]`, and the correspondence between slice and
meaning lives in whoever wrote both. That is the source of a recurring bug
class, and it makes tying, fixing, and hierarchical structure impossible to
express at all. This contract replaces it.

It is the vocabulary in which **every** ampere object — model, instrument
transformation, noise model, dataset — declares:

- what varies, under what prior, in what units, with what shape;
- what is held fixed;
- what is constant *data* rather than a parameter at all (buffers);
- which quantities in different models are secretly the same quantity (tying);
- which quantities are drawn from a shared population (plates).

It is **backend-neutral**. Nothing here imports torch, jax, numpyro or paramax,
lazily or otherwise (`architecture.md` §3–4). It is numpy, scipy,
`astropy.units` and stdlib. Turning these declarations into
`torch.distributions` / numpyro sites is W1.9's lowering spec and Phase 2's
code; this contract's obligation is to make sure the information needed to do
so survives declaration, composition and serialisation.

### Setup for the examples

```pycon
>>> import numpy as np
>>> import scipy.stats as st
>>> import astropy.units as u
>>> from ampere.core import (
...     Buffer, HierarchicalPrior, Identity, Log, Logit, Parameter,
...     Parameterised, ParameterSet, Plate, PlateBinding, Population, PriorSpec, Tie,
...     describe_prior, prior_from_spec,
... )
>>> from ampere.core.exceptions import ParameterError, TyingError

```

## 2. The objects

| Object | Role |
|---|---|
| `Parameter` | One named quantity: prior, value, shape, unit, bijection, tie label, plate membership |
| `Buffer` | One named constant array — declared explicitly, never inferred |
| `ParameterSet` | Ordered, named, immutable collection; owns the flat-vector layout, `lnprior`, `prior_transform`, sampling, (de)serialisation, `merge` |
| `BufferSet` | The same for buffers, deliberately without any prior machinery |
| `Plate` | Declarative "N members sharing hyperparameters"; expands into ordinary parameters |
| `HierarchicalPrior` | A prior whose distribution parameters are named references to other parameters |
| `Tie` | Composition-time instruction to collapse several sites into one free parameter |
| `Binding`, `ParameterMapping` | The result of `merge`: the joint set plus its wiring |
| `Bijection`, `Identity`, `Log`, `Logit` | Maps to and from unconstrained space |
| `Parameterised` | The declaration mixin: `register_parameter` / `register_buffer` / `context`, plus the opt-in `describe()` identity hook (§10) |

## 3. `Parameter`: three states, and why "fixed" is not a delta prior

A `Parameter` is in exactly one of three states.

**Free** — has a prior, is not fixed. Occupies `size` dimensions of the flat
vector and contributes to `lnprior`.

**Fixed** — has a value and *no prior*. Occupies no sampler dimension,
contributes nothing to `lnprior`, and is still handed to the model on every
evaluation.

**Deferred** — carries a `shared_as` tie label but no prior of its own: another
site in the tie group supplies it (§7). A set containing one refuses to
evaluate priors until merged.

```pycon
>>> temperature = Parameter("temperature", st.uniform(100.0, 9900.0), unit=u.K)
>>> temperature.is_free, temperature.is_fixed, temperature.size
(True, False, 1)
>>> distance = Parameter("distance", value=1.5 * u.kpc, fixed=True, unit=u.kpc)
>>> distance.is_fixed, distance.value, distance.prior is None
(True, 1.5, True)

```

Declaring both is an error, and the message says why rather than quietly
picking one:

```pycon
>>> Parameter("temperature", st.uniform(100.0, 9900.0), value=2500.0, fixed=True)
Traceback (most recent call last):
    ...
ampere.core.exceptions.ParameterError: parameter 'temperature' is declared fixed *and* given a prior. A fixed parameter has a value and no prior; if you meant a parameter that varies, drop fixed=True. A delta-function prior is not the same thing and is not supported (see docs/design/contracts/parameters.md).

```

**Decision — fixed is a state, not a degenerate prior.** A delta-function prior
would make the quantity a *free* parameter of zero-width support: it would take
a sampler dimension, contribute a divergent or arbitrary density to `lnprior`,
and behave pathologically in every sampler ampere targets. It would also lower
badly (`torch.distributions` and numpyro both have no usable point mass on the
real line). The two concepts are not interchangeable and this contract refuses
to conflate them. Fixing and unfixing are one-liners that never touch a model's
evaluation code:

```pycon
>>> held = temperature.fix(2500.0 * u.K)
>>> held.is_fixed, held.value
(True, 2500.0)
>>> held.release(st.uniform(100.0, 9900.0)).is_free
True

```

Names must be Python identifiers, because parameter values reach models as
keyword arguments. A name that would not work as one fails at declaration, not
at the first evaluation:

```pycon
>>> Parameter("2 temperatures", st.norm(0.0, 1.0))
Traceback (most recent call last):
    ...
ampere.core.exceptions.ParameterError: parameter name '2 temperatures' is not usable: every '.'-separated segment must be a valid Python identifier, because parameter values are passed to models as keyword arguments and qualified names are built by joining identifiers.

```

## 4. Priors: neutral declaration, and what "neutral" has to mean

A prior is declared as a **frozen `scipy.stats` distribution**. That is the
canonical, backend-neutral form: `scipy.stats.norm(0, 1)` names a family and
the numbers it was frozen with, and both survive introspection.

```pycon
>>> spec = describe_prior(st.norm(loc=1.0, scale=2.0))
>>> spec.family, dict(spec.kwds), spec.discrete
('norm', {'loc': 1.0, 'scale': 2.0}, False)
>>> float(prior_from_spec(spec).mean())
1.0

```

The description is **canonical**: scipy accepts the same freezing positionally
or by keyword, and `describe_prior` maps positional arguments onto their names
(the family's shape names, then `loc`, then `scale`). Two declarations of one
distribution therefore describe — and compare, tie and serialise —
identically:

```pycon
>>> describe_prior(st.norm(1.0, 2.0)) == spec
True

```

`PriorSpec` is the handle W1.9's mapping table works from. Translating scipy's
`loc`/`scale` convention into torch's or numpyro's is W1.9's job; this
contract's job is to guarantee the family name and its numeric arguments are
recoverable — from a freshly declared parameter, from a merged set, and from a
set that has been round-tripped through JSON.

**Decision — duck-typed priors evaluate but do not lower.** Any object exposing
`ppf`, `support` and `logpdf`/`logpmf` satisfies the `Prior` protocol and will
evaluate correctly on the reference path. It cannot be described neutrally, so
`describe_prior` — and therefore `ParameterSet.to_spec`, and therefore any
backend lowering and any provenance record — raises for it, naming the problem:

```pycon
>>> describe_prior(object())
Traceback (most recent call last):
    ...
ampere.core.exceptions.ParameterError: prior <object object at ...> is not a frozen scipy.stats distribution, so it cannot be described neutrally (and therefore cannot be lowered to torch/numpyro or recorded in a run's provenance). Declare priors as e.g. scipy.stats.norm(0, 1); a custom prior object may still be evaluated on the reference path, but it must not be serialised or lowered.

```

This is the deliberate middle path between "scipy frozen distributions only"
(which would block a user with a genuinely custom prior from ever running) and
"anything with a `logpdf`" (which would let an un-lowerable prior travel
silently as far as a backend, and fail there instead). Discrete families work
throughout: `log_density` accepts `logpdf` or `logpmf`, so `scipy.stats.poisson`
needs no special case.

## 5. Shape and units

**Array-valued parameters are supported.** `shape=()` is a scalar (the
default); `shape=(n,)` or richer is an array-valued parameter such as a
per-channel calibration offset. The prior applies **independently and
identically to every element**, which is the semantics per-channel offsets,
per-epoch jitters and plate members all actually want.

```pycon
>>> offsets = Parameter("offset", st.norm(0.0, 0.05), shape=(4,), unit=u.mag)
>>> offsets.size, offsets.shape
(4, (4,))

```

An array-valued parameter occupies a contiguous block of the flat vector, in C
order, and gets one label per entry — which is what corner plots and ArviZ
coordinates need (W1.8):

```pycon
>>> ParameterSet([offsets]).free_labels()
('offset[0]', 'offset[1]', 'offset[2]', 'offset[3]')

```

A genuinely *multivariate* prior — one that correlates the elements — is out of
scope for v1.3 (§12).

**Units are metadata, applied once.** A parameter's value is a plain float (or
plain array) in its declared unit; no `astropy.units.Quantity` ever enters the
hot loop, per the units trap in `DEVELOPMENT_PLAN.md` §7. Conversion happens at
declaration and at any explicit call to `to_value`/`quantity`:

```pycon
>>> d = Parameter("distance", value=1500.0 * u.pc, fixed=True, unit=u.kpc)
>>> d.value
1.5
>>> d.quantity()
<Quantity 1.5 kpc>
>>> float(d.to_value(2000.0 * u.pc))
2.0

```

Priors are declared **numerically, in the parameter's declared unit**. Ampere
does not attach units to distributions and does not convert priors; that is why
tying two parameters declared in different units is an error rather than a
silent rescaling (§7).

## 6. Unconstrained space

Gradient-based samplers need the parameter vector on the real line. A
`Bijection` says how a parameter gets there and back.

The methods are named `constrain` / `unconstrain`, not `forward` / `inverse`.
Both torch and numpyro use `forward` for unconstrained-to-constrained while a
reader coming from the statistics literature usually expects the opposite; a
reversed log-determinant is a silent, expensive bug. `constrain` and
`unconstrain` cannot be read backwards.

The default bijection is inferred from the prior's support, by the same rule
numpyro's `biject_to` applies — stated here so that the reference path and the
eventual torch/jax paths agree by construction rather than coincidence:

| Support | Bijection |
|---|---|
| `(-inf, inf)` | `Identity()` |
| `[l, inf)` | `Log(lower=l)`, i.e. `x = l + exp(y)` |
| `[l, h]` | `Logit(lower=l, upper=h)` |
| `(-inf, h]` | **no built-in** — raises, asking for an explicit declaration |

```pycon
>>> from ampere.core import default_bijection_for
>>> default_bijection_for(st.norm(0.0, 1.0))
Identity()
>>> default_bijection_for(st.uniform(loc=100.0, scale=9900.0))
Logit(lower=100.0, upper=10000.0)
>>> default_bijection_for(st.halfnorm(0.0, 1.0))
Log(lower=0.0)

```

**A discrete family is refused, with a typed capability error** (ruled
2026-09-03, `lowering.md` §12.1 — a §4.1 change, recorded in the plan's
decision log). No continuous bijection to unconstrained space can be right
for an integer-supported family, and before this rule the support-based
inference happily returned one (`Log(lower=0.0)` for a Poisson). The
refusal is `CapabilityError` — a `NotImplementedError`, deliberately *not*
a `ValueError` — because nothing is malformed: the door stays ajar for the
non-gradient routes that might eventually support discrete parameters
((variational) EM, numpyro-style enumeration, SBI, nested sampling,
Bayesian optimisation — none in the current plan), so the refusal lives
**only** here. Declaration, prior sampling, constrained-space `lnprior`
and `prior_transform` (scipy's discrete families implement `ppf`) all
work, and discreteness is queryable from the canonical description so an
engine path branches rather than catches:

```pycon
>>> default_bijection_for(st.poisson(3.0))
Traceback (most recent call last):
    ...
ampere.core.exceptions.CapabilityError: no unconstraining bijection exists for this prior: prior family 'poisson' is discrete, ...
>>> describe_prior(st.poisson(3.0)).discrete
True
>>> ParameterSet([Parameter("counts", st.poisson(3.0))]).prior_transform([0.7])
array([4.])

```

A `ParameterSet` exposes the whole round trip, with the change-of-variables
term. Stating it here, on the reference path, gives the torch and jax lowerings
an oracle to agree with rather than each rediscovering the Jacobian:

```pycon
>>> pset = ParameterSet([
...     Parameter("temperature", st.uniform(100.0, 9900.0), unit=u.K),
...     Parameter("log_tau", st.norm(0.0, 1.0)),
... ])
>>> theta = pset.prior_transform(np.array([0.5, 0.5]))
>>> theta
array([5050.,    0.])
>>> y = pset.unconstrain(theta)
>>> y
array([0., 0.])
>>> bool(np.allclose(pset.constrain(y), theta))
True
>>> # lnprior in unconstrained space = lnprior + log|d constrain / d y|
>>> bool(np.isclose(
...     pset.lnprior_unconstrained(y),
...     pset.lnprior(theta) + np.log(9900.0) - 2 * np.log(2.0),
... ))
True

```

## 7. `ParameterSet`: layout, and the asymmetry that makes fixing free

`ParameterSet` is what an inference engine talks to.

**The flat vector holds free parameters only; the value mapping holds
everything.** `pack` writes free parameters into a flat array for samplers;
`unpack` returns *all* parameters, fixed ones injected at their declared
values. This asymmetry is deliberate and load-bearing: it is what lets a user
fix a parameter, or promote a buffer to one, without a model's evaluation code
changing at all — `architecture.md` §6's stated test of this contract.

```pycon
>>> pset = ParameterSet([
...     Parameter("temperature", st.uniform(100.0, 9900.0), unit=u.K),
...     Parameter("log_tau", st.norm(0.0, 1.0)),
...     Parameter("offset", st.norm(0.0, 0.05), shape=(3,), unit=u.mag),
...     Parameter("distance", value=1.5 * u.kpc, fixed=True, unit=u.kpc),
... ])
>>> len(pset), pset.free_size
(4, 5)
>>> pset.names
('temperature', 'log_tau', 'offset', 'distance')
>>> pset.free_names, pset.fixed_names
(('temperature', 'log_tau', 'offset'), ('distance',))

```

`len()` counts *parameters*; `free_size` is the sampler's dimension. The two
differ whenever anything is fixed or array-valued, so the contract never
overloads one name for both.

The round trip is exact, fixed parameters included:

```pycon
>>> theta = pset.prior_transform(np.full(5, 0.5))
>>> values = pset.unpack(theta)
>>> values['distance'], values['offset']
(1.5, array([0., 0., 0.]))
>>> bool(np.array_equal(pset.pack(values), theta))
True

```

`pack` accepts (and ignores) entries for fixed parameters — that is exactly
what makes `pack(unpack(theta))` an identity — but rejects unknown names,
because those are almost always typos:

```pycon
>>> pset.pack({**values, 'temprature': 3000.0})
Traceback (most recent call last):
    ...
ampere.core.exceptions.ParameterError: pack() got value(s) for unknown parameter(s) ['temprature']; this set declares ['temperature', 'log_tau', 'offset', 'distance'].

```

`lnprior` accepts either form, sums over the elements of array-valued
parameters, and short-circuits to `-inf` at the first non-finite contribution:

```pycon
>>> bool(np.isclose(pset.lnprior(theta), pset.lnprior(values)))
True
>>> pset.lnprior({**values, 'temperature': 50.0})
-inf

```

### Serialisation

`to_spec` / `from_spec` produce and consume plain JSON-compatible dictionaries.
This is the round trip a run's provenance record needs (`DEVELOPMENT_PLAN.md`
§4.6 wants a spec hash in the `InferenceData` attrs) and the reason priors must
be describable at all.

```pycon
>>> import json
>>> rebuilt = ParameterSet.from_spec(json.loads(json.dumps(pset.to_spec())))
>>> rebuilt == pset
True
>>> bool(np.isclose(rebuilt.lnprior(theta), pset.lnprior(theta)))
True

```

## 8. Tying and sharing across models and datasets

Two datasets of the same object share its distance. Two models fitted jointly
share a calibration scale. This must cost **one** sampler dimension with
several binding sites — not N dimensions plus a constraint that they be equal,
which wastes the sampler's time and corrupts its geometry.

Ampere supports both directions of declaration, because both situations are
real:

**Declaration-time** (`shared_as="label"`) — the person writing the model knows
the quantity is shared.

**Composition-time** (`Tie`) — the models were written independently, or came
from a library, and the person *composing the fit* is the one who knows two
distances are the same distance.

`ParameterSet.merge` qualifies each component's names with its label
(`"spectrum.temperature"`) so independently written models never collide, and
collapses tie groups into single parameters:

```pycon
>>> sed = ParameterSet([
...     Parameter("temperature", st.uniform(100.0, 9900.0), unit=u.K),
...     Parameter("distance", st.norm(1.5, 0.1), unit=u.kpc, shared_as="distance"),
... ])
>>> spectrum = ParameterSet([
...     Parameter("temperature", st.uniform(100.0, 9900.0), unit=u.K),
...     Parameter("distance", st.norm(1.5, 0.1), unit=u.kpc, shared_as="distance"),
...     Parameter("scale", st.norm(1.0, 0.05)),
... ])
>>> sed.free_size + spectrum.free_size            # five, if merged naively
5
>>> mapping = ParameterSet.merge({"sed": sed, "spectrum": spectrum})
>>> mapping.merged.names
('sed.temperature', 'distance', 'spectrum.temperature', 'spectrum.scale')
>>> mapping.merged.free_size                       # four: the tied pair is one
4
>>> mapping.tied_names
('distance',)
>>> [(b.component, b.local_name) for b in mapping.sites_of("distance")]
[('sed', 'distance'), ('spectrum', 'distance')]

```

`merge` returns a `ParameterMapping`, not a bare `ParameterSet`. The bindings
*are* the point: without them, tying would be a naming convention rather than a
structure, and nothing could route a sampled value back to the components that
consume it. `distribute` does that routing, producing exactly the keyword
arguments each component expects — tied values appearing (identically) in every
component that binds them, fixed values included:

```pycon
>>> theta = mapping.merged.prior_transform(np.full(4, 0.5))
>>> routed = mapping.distribute(theta)
>>> sorted(routed)
['sed', 'spectrum']
>>> routed["sed"]["distance"] == routed["spectrum"]["distance"]
True

```

The composition-time form is equivalent, for models that never anticipated
being shared:

```pycon
>>> plain_sed = ParameterSet([
...     Parameter("temperature", st.uniform(100.0, 9900.0), unit=u.K),
...     Parameter("distance", st.norm(1.5, 0.1), unit=u.kpc),
... ])
>>> plain_spec = ParameterSet([
...     Parameter("temperature", st.uniform(100.0, 9900.0), unit=u.K),
...     Parameter("distance", st.norm(1.5, 0.1), unit=u.kpc),
... ])
>>> composed = ParameterSet.merge(
...     {"sed": plain_sed, "spectrum": plain_spec},
...     ties=[Tie("distance", ("sed.distance", "spectrum.distance"))],
... )
>>> composed.merged.names, composed.merged.free_size
(('sed.temperature', 'distance', 'spectrum.temperature'), 3)

```

### What tied sites must agree about

Shape, unit, prior, and any explicitly declared bijection. Disagreement is an
error with a message naming both sites, never a silent choice:

```pycon
>>> ParameterSet.merge({
...     "a": ParameterSet([Parameter("d", st.norm(1.0, 1.0), unit=u.kpc, shared_as="d")]),
...     "b": ParameterSet([Parameter("d", st.norm(1.0, 1.0), unit=u.pc, shared_as="d")]),
... })
Traceback (most recent call last):
    ...
ampere.core.exceptions.TyingError: tied parameters disagree about units: a.d is in kpc but b.d is in pc. Priors are declared numerically in the parameter's unit, so ampere will not guess a conversion; declare both in the same unit.

```

Units are *equivalent* here, not merely different — kpc and pc convert freely.
Ampere still refuses, because the priors were declared numerically in those
units and rescaling a `scipy.stats` distribution's parameters correctly is
family-specific (`loc` scales, `scale` scales, a shape parameter may not). A
loud error with a one-line fix beats a clever conversion that is right for
`norm` and wrong for `lognorm`.

Hierarchical priors are compared **after their references are qualified**: two
sites are only the same prior if their hyperparameters resolve to the same
merged parameters. Anything else would force the collapsed parameter to adopt
one component's hyperparameters and silently orphan the other's:

```pycon
>>> def hier_component():
...     return ParameterSet([
...         Parameter("mu", st.norm(0.0, 5.0)),
...         Parameter("theta", HierarchicalPrior("norm", {"loc": "mu"}), shared_as="theta"),
...     ])
>>> ParameterSet.merge({"a": hier_component(), "b": hier_component()})
Traceback (most recent call last):
    ...
ampere.core.exceptions.TyingError: tie 'theta' has disagreeing priors: a.theta declares norm(a.mu) but b.theta declares norm(b.mu). Tied parameters are one parameter and must have one prior. If the hyperparameters are themselves the same quantity, tie them too: hierarchical references are compared after qualification.

```

Adding `shared_as="mu"` to both `mu` declarations makes the two references
resolve to one merged `"mu"`, and the tie then collapses cleanly.

A site may also *defer* its prior entirely — useful when one model is the
authority on a quantity and another merely consumes it. The consuming set is
then not self-sufficient, and says so rather than silently treating a missing
prior as improper:

```pycon
>>> authority = ParameterSet([
...     Parameter("distance", st.norm(1.5, 0.1), unit=u.kpc, shared_as="distance")])
>>> consumer = ParameterSet([Parameter("distance", unit=u.kpc, shared_as="distance")])
>>> consumer.is_resolved, consumer.deferred_names
(False, ('distance',))
>>> consumer.lnprior(np.array([1.5]))
Traceback (most recent call last):
    ...
ampere.core.exceptions.TyingError: cannot evaluate lnprior: parameter(s) ['distance'] are declared shared (shared_as=...) without a prior, so their priors come from another site. Merge this set with the one that supplies them (ParameterSet.merge) before evaluating priors.
>>> joint = ParameterSet.merge({"one": authority, "two": consumer})
>>> joint.merged.names, joint.merged.free_size, joint.merged.is_resolved
(('distance',), 1, True)

```

### Lossless nesting: a `ParameterMapping` as a component

**Ruled by Peter, 2026-09-02** ("lossless, not necessarily recursive"):
`merge` accepts an already-merged `ParameterMapping` as a component. Its
`merged` set joins the merge exactly as a plain set would — same names, same
free dimensions, same routing — and the mapping is *retained*, so the
result's public `bindings` view composes down to the ultimate leaves. What
this buys is truthful introspection under the nested topology `inference.md`
§4 ratified: a parameter collapsed by an *inner* merge (two instrument steps
sharing a `shared_as` label, say) used to reach the outer mapping as a
single binding, so `tied_names` under-reported and a provenance consumer
walking `bindings` saw no sharing.

```pycon
>>> instrument = ParameterSet.merge({
...     "a": ParameterSet([Parameter("scale", st.lognorm(0.2), shared_as="gain")]),
...     "b": ParameterSet([Parameter("scale", st.lognorm(0.2), shared_as="gain")]),
... })
>>> nested = ParameterSet.merge({
...     "instrument": instrument,
...     "model": ParameterSet([Parameter("t", st.norm(0.0, 1.0))]),
... })
>>> nested.merged.names
('instrument.gain', 'model.t')
>>> nested.tied_names
('instrument.gain',)
>>> [(b.component, b.local_name) for b in nested.sites_of("instrument.gain")]
[('instrument', 'a.scale'), ('instrument', 'b.scale')]

```

Routing and introspection are two views of one structure. `distribute` stays
**one level deep** — the component receives its own merged names and
re-distributes with its retained mapping, which is what keeps
`Instrument.__call__`'s values path working unchanged (the nesting rule,
`inference.md` §4.5) — while `bindings`/`sites_of`/`tied_names` compose
through the retained mappings and tell the leaf-level truth. The raw
one-level table remains available as `routing`.

```pycon
>>> routed = nested.distribute({"instrument.gain": 2.0, "model.t": 0.1})
>>> routed["instrument"]
{'gain': 2.0}
>>> instrument.distribute(routed["instrument"])
{'a': {'scale': 2.0}, 'b': {'scale': 2.0}}

```

### Plate bindings: routing one element to one component

**Ruled by Peter, 2026-09-02** (the population sketch's gap H-2): `Binding`
carries an optional element `index`, and `merge` accepts `PlateBinding`
declarations that create such bindings — so one member of a `Plate`'s
array-valued parameter reaches one component as a scalar under a bare local
name. The merged set is untouched: one array-valued parameter, one sample
site, so the numpyro lowering in `lowering.md` is unaffected; only the
routing changes. This is *addressing*, not tying — each member remains its
own draw, which is why limitation §12.3's refusal of ties across plate
members stands — and element bindings accordingly do not count towards
`tied_names`.

```pycon
>>> survey = ParameterSet([], plates=[Plate(
...     "objects", size=3,
...     hyperparameters=[Parameter("mu", st.norm(0.0, 5.0)),
...                      Parameter("sigma", st.halfnorm(0.0, 2.0))],
...     members=[Parameter("theta", HierarchicalPrior("norm", {"loc": "mu", "scale": "sigma"}))],
... )])
>>> population = ParameterSet.merge(
...     {"population": survey,
...      "obj0": ParameterSet([Parameter("cal", st.lognorm(0.1))]),
...      "obj1": ParameterSet([Parameter("cal", st.lognorm(0.1))])},
...     plate_bindings=[
...         PlateBinding("population.objects.theta", "obj0", "theta", 0),
...         PlateBinding("population.objects.theta", "obj1", "theta", 1),
...     ],
... )
>>> routed = population.distribute({
...     "population.objects.mu": 0.0, "population.objects.sigma": 1.0,
...     "population.objects.theta": np.array([10.0, 20.0, 30.0]),
...     "obj0.cal": 1.1, "obj1.cal": 0.9,
... })
>>> routed["obj0"]["cal"], float(routed["obj0"]["theta"])
(1.1, 10.0)
>>> routed["obj1"]["cal"], float(routed["obj1"]["theta"])
(0.9, 20.0)

```

(The element arrives as a numpy scalar — dtype-preserving, and arithmetically
a float.)

`PlateBinding.parameter` is the **fully qualified merged name** (after
qualification and tie collapse) — fully qualified rather than the sketch's
bare form, because two components may each hold a plate of the same local
name. Everything is validated at merge: the parameter must exist and be
array-valued, the index in range for its shape, the component present, and
the local name must not shadow anything the component already receives. Note
the receiving component's own `ParameterSet` does **not** declare the local
name: the element arrives as an extra key in `distribute`'s output, and
consuming it is the composing caller's contract — the intended caller being
the plate-of-datasets construction, which builds these bindings from its own
dataset ordering. **Since W5.12 that caller exists**: `Population` (§9) is the
declaration that produces these bindings, and `DatasetCollection.plate`
(`inference.md` §9) is the convenience that derives them from the datasets'
order, so writing them out by hand is the low-level route rather than the
expected one.

## 9. Hierarchical structure: `HierarchicalPrior` and `Plate`

The plan's design horizon (d) asks that population models stay expressible and
lower cleanly to pyro/numpyro. Two constructs do that, and between them they
cover both shapes the problem takes.

### `HierarchicalPrior` — a prior that references other parameters

The family is a neutral name, exactly as for any other prior, so it goes
through the same W1.9 lowering table. The hyperparameters are recorded as
**references by name** — which is what a numpyro model does anyway
(`dist.Normal(mu, sigma)`, where `mu` and `sigma` are earlier sample sites).

```pycon
>>> hp = HierarchicalPrior("norm", {"loc": "mu", "scale": "sigma"})
>>> hp.references
('mu', 'sigma')
>>> float(hp.bind({"mu": 2.0, "sigma": 0.5}).mean())
2.0

```

A `ParameterSet` containing hierarchical priors orders its prior evaluations
topologically, so a hyperparameter is always resolved before anything that
references it. References that do not resolve, and cycles, fail at
construction:

```pycon
>>> ParameterSet([Parameter("theta", HierarchicalPrior("norm", {"loc": "mu"}))])
Traceback (most recent call last):
    ...
ampere.core.exceptions.ParameterError: parameter 'theta' has a hierarchical prior referencing 'mu', which is not in this set (it declares ['theta']). Hierarchical references must resolve within the set that will evaluate them.
>>> ParameterSet([
...     Parameter("a", HierarchicalPrior("norm", {"loc": "b"})),
...     Parameter("b", HierarchicalPrior("norm", {"loc": "a"})),
... ])
Traceback (most recent call last):
    ...
ampere.core.exceptions.ParameterError: hierarchical priors form a cycle: a -> b -> a. A parameter's prior cannot depend (even indirectly) on itself.

```

### `Plate` — N members sharing hyperparameters

This is the canonical population declaration: *N objects share `mu` and
`sigma`, each has its own `theta_i ~ Normal(mu, sigma)`.*

```pycon
>>> objects = Plate(
...     "objects",
...     size=4,
...     hyperparameters=[
...         Parameter("mu", st.norm(0.0, 5.0)),
...         Parameter("sigma", st.halfnorm(0.0, 2.0)),
...     ],
...     members=[
...         Parameter("theta", HierarchicalPrior("norm", {"loc": "mu", "scale": "sigma"})),
...     ],
... )
>>> population = ParameterSet([Parameter("background", st.norm(0.0, 1.0))], plates=[objects])
>>> population.names
('background', 'objects.mu', 'objects.sigma', 'objects.theta')
>>> population.plates
{'objects': 4}
>>> population.free_size                     # 1 + 1 + 1 + 4
7

```

A `Plate` is a *constructor*, not a container: it expands into ordinary
`Parameter`s, which is how everything else in this module — pack/unpack,
priors, merging, serialisation — works on it with no special cases. Its
hyperparameters become `"objects.mu"` and `"objects.sigma"`; each member
becomes a **single array-valued parameter** of shape `(size,) + member.shape`,
tagged with its plate, with its hyperparameter references already rewritten to
the qualified names:

```pycon
>>> theta_i = population["objects.theta"]
>>> theta_i.shape, theta_i.plate, theta_i.references
((4,), 'objects', ('objects.mu', 'objects.sigma'))

```

That shape is chosen because it *is* the numpyro lowering: one
`numpyro.sample` inside one `numpyro.plate("objects", 4)` produces a batch of
four draws, not four separate sites. The declared structure is inspectable —
plate name, plate size, per-member priors with resolved references — which is
everything W1.9 needs to emit that.

The structure is sound numerically, not merely declarative. Each member is
drawn from the population distribution implied by the current hyperparameter
values:

```pycon
>>> drawn = population.unpack(population.prior_transform(np.full(7, 0.6)))
>>> mu, sigma = drawn["objects.mu"], drawn["objects.sigma"]
>>> bool(np.allclose(drawn["objects.theta"], st.norm(mu, sigma).ppf(0.6)))
True
>>> sample = population.sample(np.random.default_rng(0))
>>> sample["objects.theta"].shape
(4,)
>>> population.free_labels()[3:]
('objects.theta[0]', 'objects.theta[1]', 'objects.theta[2]', 'objects.theta[3]')

```

### `regularised_horseshoe` — sparsity over a set of amplitudes (*Added W5.8*)

The one named declaration this contract ships, because it is the one the plan
asks for by name: the sparsity guard on the amplitudes of a `Sum` of noise
components (`likelihoods.md` §6). It is `HierarchicalPrior` used three levels
deep and nothing else — `τ ~ C⁺(0, τ₀)`, `s_j | τ ~ C⁺(0, τ)`, `a_j | s_j ~
N⁺(0, s_j)` — returning the `Parameter`s in that order, because a
`ParameterSet` refuses a reference that is not already in it and the levels
therefore have to be registered outermost first. The horseshoe's multiplicative
`s_j = τ λ_j` is written as one declaration rather than as a product, which it
may be because the half-Cauchy is a **scale family**: `τ·C⁺(0,1)` and
`C⁺(0, τ)` are the same distribution. That is what lets the whole prior be
declared with a construct that references a parameter by name and has no
product node. Piironen & Vehtari's slab is the one part that cannot: `λ̃_j² =
c²λ_j²/(c² + τ²λ_j²)` is a deterministic function of two sampled parameters,
and this contract declares parameters and priors, not deterministic nodes —
`tail="regularised"` gives the declarable form of the same guard instead (a
`gamma(a=½, scale=τ)` local level: the same spike at zero, an exponential tail
rather than a Cauchy one), and a `Derived` node is what would close the gap.

### `Population` — hierarchy declared at composition time (*Added W5.12*)

**Ruled by Peter, 2026-09-03** (`hierarchical_population.md` §11 Q2): the
sketch's gap **H-1** lands with Phase 5, and W5.12 is where it landed. The gap
was an asymmetry. §8 supports tying in *both* directions, at declaration time
and at composition time, "because the models were written independently — or
came from a library — and the person composing the fit is the one who knows".
Hierarchy had only one direction: a `HierarchicalPrior` must be on the site
when its `ParameterSet` is built, and its references must already resolve
there, so a user with N library models had to rewrite every one of them —
each declaring population hyperparameters it never uses — before they could
say "these are draws from a population". `Tie` is not the escape: it collapses
N sites into **one**, which fits one θ for the whole survey, silently and
plausibly.

`Population` is the counterpart `Tie` never had. It names the members once,
and `merge` does the rewriting:

```pycon
>>> def object_set():
...     return ParameterSet([Parameter("theta", st.norm(0.0, 1.0)),
...                          Parameter("cal", st.lognorm(0.1))])
>>> objects = Population(
...     "objects",
...     members=[Parameter("theta", HierarchicalPrior("norm", {"loc": "mu", "scale": "sigma"}))],
...     hyperpriors=[Parameter("mu", st.norm(0.0, 5.0)),
...                  Parameter("sigma", st.halfnorm(0.0, 2.0))],
...     over=["obj0", "obj1", "obj2"],
... )
>>> survey = ParameterSet.merge(
...     {label: object_set() for label in ("obj0", "obj1", "obj2")},
...     populations=[objects],
... )
>>> survey.merged.names
('obj0.cal', 'obj1.cal', 'obj2.cal', 'objects.mu', 'objects.sigma', 'objects.theta')
>>> survey.merged["objects.theta"].shape, survey.merged["objects.theta"].plate
((3,), 'objects')
>>> [(b.component, b.local_name, b.index) for b in survey.bindings if b.index is not None]
[('obj0', 'theta', 0), ('obj1', 'theta', 1), ('obj2', 'theta', 2)]

```

Three things happened there, and each is the answer to one of H-1's
complaints. The hyperpriors joined the merge as one further component, so no
per-object model had to declare them. Each member component's *own* `theta` —
an ordinary `Normal(0, 1)` a library model declared without knowing about this
fit — was replaced by its draw from the population. And the draws are one
array-valued parameter with three element bindings (§8), so `obj0` receives a
scalar `theta` under its own local name and never learns that it is object
zero of three.

`members` is `Plate`'s `members`, and it means the same thing: one entry per
quantity each member has its own value of. An entry whose prior is a
`HierarchicalPrior` is the population draw; an entry with an ordinary prior is
a **per-member nuisance parameter**, an i.i.d. array — "each object has its own
calibration scale, from a common prior" — which is the reading §11 Q4
confirms. `hyperpriors` are ordinary `Parameter`s: nothing about them is
special, which is what lets the same machinery fit them.

#### Two layouts, one density

```pycon
>>> flat = ParameterSet.merge(
...     {label: object_set() for label in ("obj0", "obj1", "obj2")},
...     populations=[Population(
...         "objects", members=objects.members, hyperpriors=objects.hyperpriors,
...         over=objects.over, layout="flat",
...     )],
... )
>>> flat.merged.names  # doctest: +NORMALIZE_WHITESPACE
('obj0.theta', 'obj0.cal', 'obj1.theta', 'obj1.cal', 'obj2.theta', 'obj2.cal',
 'objects.mu', 'objects.sigma')
>>> flat.merged["obj0.theta"].references
('objects.mu', 'objects.sigma')
>>> flat.merged.free_size == survey.merged.free_size
True

```

`layout="plate"` — the default above — is one array-valued site routed by
element, which is what a `numpyro.plate` or a `pyro.plate` *is*, and is
therefore what the torch and jax realisations lower (`inference.md` §10a).
`layout="flat"` is the N-component pattern below: each member keeps its own
scalar parameter, re-priored onto the shared hyperpriors. The two declare the
same joint density over the same number of dimensions, and the conformance
suite holds them to that on every backend, realised and unrealised.

The flat layout is **refused above `MAX_FLAT_MEMBERS` (128)**, naming the
plate layout as the remedy. This is the cost recorded below made into a
guardrail: `lnprior` is O(number of `Parameter` objects), and past a hundred or
so members the flat declaration is a mistake rather than a trade-off. The
plate layout has no such limit — it is one object however many members it
holds.

#### What it refuses

By name, at merge, because each of these is a different model from the one the
caller meant: a component label the merge already has; a member component that
does not exist, or the population's own; a member that disagrees with the
component's own declaration about shape or unit; a **fixed** site (a fixed site
has no draw); a site already tied or `shared_as` (sharing collapses N sites
into one value, a population keeps them N — pick one); and a component merged
as a `ParameterMapping`, because a composite's merged names are qualified
(`likelihood.scale`) while a population addresses its members by bare local
name, so the draw would never reach a leaf. For a plate of per-object
*datasets*, declare the population over the models those datasets name —
which is what `DatasetCollection.plate` does.

### Which construct to use

Three patterns are expressible, and they serve different data layouts:

- **One component, N members vectorised** (`Plate`): the N objects are fitted
  from one dataset, or from data already stacked. Lowers to a numpyro plate.
- **N components, declared at composition time** (`Population`, *W5.12*): the
  N objects were written independently and each has its own dataset, and the
  person composing the fit is the one who knows they are a population. This is
  the one to reach for when the models are somebody else's; it produces the
  first pattern's plate (or, at `layout="flat"`, the third pattern's
  parameters) without any of them being rewritten.
- **N components, declared by hand** (`HierarchicalPrior` + tying): each
  object gets its own `ParameterSet` with a scalar `theta`, whose
  `HierarchicalPrior` references `mu` and `sigma` — themselves tied across all
  components, so they collapse to one shared pair on merge. This is the shape
  a `DatasetCollection` (W1.7) of per-object datasets naturally takes, and it
  is what `Population(layout="flat")` builds for you.

The second pattern in full, for two objects:

```pycon
>>> def object_set():
...     return ParameterSet([
...         Parameter("mu", st.norm(0.0, 5.0), shared_as="mu"),
...         Parameter("sigma", st.halfnorm(0.0, 2.0), shared_as="sigma"),
...         Parameter("theta", HierarchicalPrior("norm", {"loc": "mu", "scale": "sigma"})),
...     ])
>>> pop = ParameterSet.merge({"obj1": object_set(), "obj2": object_set()})
>>> pop.merged.names
('mu', 'sigma', 'obj1.theta', 'obj2.theta')
>>> pop.merged.free_size                     # not 6: mu and sigma are shared
4
>>> pop.merged["obj1.theta"].references      # references follow the collapse
('mu', 'sigma')

```

**Cost, at population scale** (recorded at the freeze — the W1.11
population sketch's gap H-4): `lnprior` and its relatives are O(number of
`Parameter` objects), so the N-component pattern above evaluates N scalar
priors through N objects per call, where a `Plate` evaluates one
vectorised prior — measured at N = 1000 as roughly two orders of magnitude
apart. That is not a defect (`ParameterSet` is a declaration container,
not a hot-loop object, and an engine pays this once per proposal beside a
model evaluation), but at genuinely population scale prefer the `Plate`
layout where the data allow it, and expect the N-component layout's prior
overhead to be visible beside a cheap model. *W5.12* turned that guidance
into a refusal for the one declaration that can always be expressed the other
way: `Population(layout="flat")` is refused above `MAX_FLAT_MEMBERS`.

## 10. Buffers: explicit, and why

`architecture.md` §6 gives the distinguishing question: *would you ever put a
prior on it, sample it, or take a gradient with respect to it?* If no, it is a
buffer — wavelength grids, opacity tables, filter curves, response matrices.

**Buffers are declared explicitly** (decided in `architecture.md` §6, and
implemented here as the binding API). The alternative — treating any
unregistered array attribute as a buffer — saves boilerplate and fails
silently: a mistyped or forgotten parameter registration becomes an untracked
buffer, which is exactly the silent-drift bug class these contracts exist to
end.

`BufferSet` is deliberately *not* a `ParameterSet`. It has no `lnprior`, no
`pack`/`unpack` and no prior transform, because a buffer has no prior and
occupies no sampler dimension. The two collections meeting only at
`Parameterised.context` is the contract's structural statement that parameters
and buffers are different kinds of thing.

Stored buffer arrays are read-only — a buffer mutated in place after
registration would invalidate exactly the caches its constancy justifies:

```pycon
>>> wl = Buffer("wavelength", np.geomspace(1.0, 100.0, 5), unit=u.micron)
>>> wl.shape, wl.dtype
((5,), dtype('float64'))
>>> wl.value[0] = 2.0
Traceback (most recent call last):
    ...
ValueError: assignment destination is read-only

```

### `Parameterised`: one namespace for evaluation

A model-like object declares both through the `Parameterised` mixin, and reads
both out of the single mapping `context()` returns:

```pycon
>>> class ModifiedBlackbody(Parameterised):
...     def __init__(self, wavelength):
...         self.register_buffer("wavelength", wavelength, unit=u.micron)
...         self.register_buffer("beta", 1.8)
...         self.register_parameter(
...             Parameter("temperature", st.uniform(100.0, 9900.0), unit=u.K)
...         )
...     def __call__(self, **values):
...         ctx = self.context(values)
...         return ctx["temperature"] * ctx["wavelength"] ** -ctx["beta"]
>>> model = ModifiedBlackbody(np.array([1.0, 2.0, 4.0]))
>>> model.parameters.names, model.buffers.names
(('temperature',), ('wavelength', 'beta'))
>>> model.parameters.free_size
1
>>> np.round(model(temperature=100.0), 4)
array([100.    ,  28.7175,   8.2469])

```

Written that way, **promoting a buffer to a parameter is a configuration
change, never a code change** — which is the test `architecture.md` §6 sets for
this contract. The dust emissivity index `beta` is a buffer by default; a user
who wants to fit it says so at composition time, and `__call__` is untouched:

```pycon
>>> model.promote_buffer("beta", prior=st.norm(1.8, 0.2))
Parameter('beta', prior=norm(loc=1.8, scale=0.2))
>>> model.parameters.names, model.buffers.names
(('temperature', 'beta'), ('wavelength',))
>>> model.parameters.free_size
2
>>> np.round(model(temperature=100.0, beta=1.8), 4)     # same call, same answer
array([100.    ,  28.7175,   8.2469])

```

Promotion without a prior gives a **fixed** parameter: the quantity joins the
model's declared vocabulary and its provenance record without yet varying. The
inverse, `demote_parameter`, exists too, so the three states form a cycle a
user can move around freely.

Buffers and parameters occupy one namespace and may not collide; the collision
is caught at registration, not at the first ambiguous lookup:

```pycon
>>> model.register_parameter(Parameter("wavelength", st.norm(0.0, 1.0)))
Traceback (most recent call last):
    ...
ampere.core.exceptions.ParameterError: ModifiedBlackbody already declares a buffer 'wavelength'

```

#### The reserved names — *Amended W5.20*

*Ruled by Peter 2026-09-15; decision-log entry in `DEVELOPMENT_PLAN.md` §2.*

A parameter or buffer may not take a name from a **stated, finite reserved
set**, and **that set is core's and identical on every backend**. It has three
sources, spelled once in `ampere.core.parameter.reserved_names()`:

* every public name of `Parameterised` — `parameters`, `buffers`, `context`,
  `describe`, the `register_*`/`promote_buffer`/`demote_parameter` API;
* every public name of `Model` — `evaluate`, `evaluate_batch`, `call_batch`,
  `compile_for`, and the four capability flags;
* `TORCH_MODULE_NAMES`, the public namespace of `torch.nn.Module` (`to`,
  `type`, `float`, `apply`, `train`, `state_dict`, …). Torch's lowering nests
  every parameter as an attribute of an `nn.Module` (`lowering.md` §6.1), so
  these are the names no backend could carry whatever core thought of them.
  Core must not import torch, so the list is pinned as a literal and a test in
  the torch environment holds it to the real class.

The first two are *computed from the classes*, so the set cannot drift from the
code; the third is checked against torch. **A backend may add public attributes
to its model classes freely without changing which parameter names are legal.**

That last sentence is the amendment. The rule used to be
`hasattr(type(self), name)` over the whole MRO, which made the legal parameter
names a property of *which backend's base class a model inherits*: a torch
spectral model reserved `grid`, `grid_tensor`, `to` and `AXIS`, a jax one
reserved a different subset, and an interferometric source model could declare
a parameter called `flux` on the reference class but not on its twin — the same
declaration, legal on one backend and refused on another, which is not a thing
a *core* contract may say. Nothing in ampere reads a parameter as an attribute
(there is no `__getattr__` on `Parameterised`; values arrive as
`context[name]`), so the old rule guarded a convention rather than a live
defect, and the convention is now stated once:

```pycon
>>> class Native(Parameterised):
...     AXIS = "spectral_axis"
...     def flux(self, channel): ...          # a backend's native surface
>>> _ = Native().register_parameter(Parameter("flux", st.norm(0.0, 1.0)))
>>> Native().register_parameter(Parameter("to", st.norm(0.0, 1.0)))
Traceback (most recent call last):
    ...
ampere.core.exceptions.ParameterError: 'to' is a reserved name: it belongs to the core parameter namespace (`parameters.md` §10), which is the same set on every backend. Pick another name.

```

The backends' own native value-and-coordinates surface is spelled
`native_flux`/`native_grid` — canonical since the same ruling, `flux`/`grid`
kept as a legacy alias — which is now a matter of taste rather than of
necessity: `inference.md` §10a has the reasoning.

Buffers are excluded from everything prior-related, and parameters from the
buffer set — the separation this contract exists to make:

```pycon
>>> "beta" in model.buffers, "wavelength" in model.parameters
(False, False)
>>> sorted(model.context({"temperature": 100.0, "beta": 1.8}))
['beta', 'temperature', 'wavelength']

```

Per-backend lowering of buffers (torch `register_buffer`; jax array fields
excluded from the trainable partition, **never** equinox static fields, which
hash array contents into the JIT cache key — `DEVELOPMENT_PLAN.md` §7) is
W1.9's table, not this document's.

### `describe()`: configuration that is neither parameter nor buffer

*Added post-freeze, ruled by Peter 2026-09-03 at the freeze's escalations and
landed with W2.1; decision-log entry in `DEVELOPMENT_PLAN.md` §2. The
motivating problem and the provenance half are `results.md` §9/§13.13/§14.*

Parameters and buffers do not exhaust what changes a model's output. A
`Redden(law="ccm89")` selects its extinction curve with a plain string: not a
quantity anyone would put a prior on, so not a parameter, and not an array, so
not a buffer. Provenance therefore cannot see it, and two such fits share a
cache key while scoring differently.

Hashing an arbitrary `__dict__` is not a safe general answer — it would sweep in
caches, file handles, open file descriptors and unhashable state — so the
answer is an **opt-in** declaration instead:

```
Parameterised.describe() -> Mapping[str, Any] | None
```

The default returns `None`: nothing extra declared. A class whose behaviour
depends on something the contracts do not model overrides it:

```python
def describe(self):
    return {"law": self.law}
```

`results.md` §9's `model_fingerprint` folds the result in. Three obligations on
what a payload may contain:

* it must be **normalisable** by `results.md` §9's recipe — strings, numbers,
  booleans, `None`, arrays, and lists or mappings of those. Anything else is
  refused loudly when the hash is taken, never silently omitted;
* it must be **deterministic**, since it becomes part of a cache key;
* it must be **backend-neutral**. A device string, a dtype, or an array's
  backing type does not belong in one: `results.md` §14's derived neutral
  identity and W1.10's cross-backend equivalence row both compare `describe()`
  payloads across implementations of the same declaration.

One consequence for §10's namespace rule: `describe` is a class attribute of
`Parameterised`, so it joins `parameters`, `buffers` and `context` among the
names a parameter or buffer may not shadow.

## 11. Decisions and their reasoning

| Decision | Reasoning |
|---|---|
| Priors declared as frozen `scipy.stats` distributions; `PriorSpec` is the neutral description | Family name plus numeric arguments is the minimum W1.9 needs and the maximum that translates across scipy/torch/numpyro. Anything richer would bake in one library's semantics |
| `describe_prior` canonicalises to keyword form | scipy accepts the same freezing positionally or by keyword; one distribution must have one description, so equality, tying and W1.9's table see a single form |
| Tied hierarchical priors must reference the same merged hyperparameters | Compared after qualification; the alternative silently wires the collapsed parameter to one component's hyperparameters and orphans the other's |
| Tied sites' explicit bijections must agree | First-declared-wins was a silent choice; disagreement errors like every other tie disagreement |
| Duck-typed priors evaluate but do not lower or serialise | Neither locks out custom priors nor lets an un-lowerable one travel silently to a backend |
| Fixed is a state, not a delta prior | A delta prior is a degenerate free parameter: a wasted sampler dimension with a pathological density, and no clean lowering |
| Flat vector = free parameters; value mapping = all parameters | Makes fixing, and buffer promotion, invisible to model code — `architecture.md` §6's test |
| Array-valued parameters supported, i.i.d. priors across elements | Per-channel offsets and plate members both need it, and it costs one `shape` field plus C-order flattening |
| Bijection methods named `constrain`/`unconstrain` | `forward`/`inverse` are ambiguous across the libraries we lower to; a reversed Jacobian is a silent bug |
| Bijection inferred from support, by numpyro's rule | Reference and native paths agree by construction, not coincidence |
| No built-in bijection for `(-inf, h]`; raises | Better a loud gap with an obvious extension point than a silently wrong default |
| Tying by label (declaration) **and** by `Tie` (composition) | Both situations are real; a library model cannot anticipate the fit it will be composed into |
| Tied sites must agree on unit exactly, not merely be convertible | Priors are numeric in the declared unit; rescaling a distribution correctly is family-specific |
| `merge` returns `ParameterMapping`, not `ParameterSet` | The bindings are the structure; dropping them would reduce tying to a naming convention |
| Plates expand into ordinary parameters | Everything else works on them unchanged; the array-valued member *is* the numpyro plate lowering |
| Hierarchical hyperparameters referenced by name, not by callable | A callable is opaque to lowering and serialisation; a name reference is exactly what a PPL sample site is |
| Buffers explicit, and in a separate collection type | `architecture.md` §6; a typo must never silently become an untracked buffer |
| `as_paramax()` removed from the core API | `DEVELOPMENT_PLAN.md` §6: "Paramax is a lowering mechanism, not a user API", and core may not know about optional dependencies (§4). Kept as a raising stub so a porter finds an explanation, not an `AttributeError` |

## 12. Deliberate limitations of v1.3

Each of these is a decision, not an oversight. Each has an extension point.

1. **No multivariate priors.** An array-valued parameter's prior is i.i.d.
   across elements. A correlated prior (a covariance across calibration
   channels, say) needs a `MultivariatePrior` slot alongside `Prior`; the
   `Parameter.prior` field already accommodates a new type without a
   layout change.
2. **One plate dimension per parameter.** Nested plates (objects within
   surveys) are not expressible. `Parameter.plate` would become a tuple.
3. **Tying across plate members is refused.** A tie group whose sites carry a
   plate raises: tying makes N members one draw, which contradicts what a
   plate is. The index-alignment concept this limitation originally deferred
   **now exists** (ruled 2026-09-02): `Binding.index` and `merge`'s
   `plate_bindings` route one *element* to one component (§8), which is
   addressing rather than tying, so the refusal stands unchanged for genuine
   ties.
4. **`merge` is not associative in its bindings — unless the mapping is
   passed.** Merging a `ParameterMapping.merged` (the bare set) again routes
   correctly but silently drops the first merge's bindings from
   introspection. Since 2026-09-02 the lossless route exists and is the
   recommendation for composing composites: pass the `ParameterMapping`
   itself as the component (§8) and the outer bindings compose to the
   leaves. Merging all leaf components in one flat call remains equivalent
   where the structure allows it.
5. **Bijections are declared per parameter, not per set.** A joint bijection
   (a Cholesky factor over several parameters) is not expressible.
6. **No built-in bijection for supports bounded above only.** §6.
7. **Units must match exactly for tying**, not merely be convertible. §8.
8. **Discrete parameters** evaluate and transform correctly (`log_density`
   handles `logpmf`, `ppf` is well defined) but no sampler in ampere currently
   consumes them; treat the support as declared-but-unexercised. *(Amended at
   the freeze, ruled 2026-09-03: the one thing that could never be right —
   a continuous default bijection — is now refused with a typed
   `CapabilityError` in `default_bijection_for`, and only there; see §6.)*

## 13. What this contract hands to the specs downstream

Obligations and hooks the later contract specs should reconcile against.

- **W1.4 (ModelResult schema)** — nothing structural, but note that channel
  names and parameter names live in *different* namespaces and may collide
  harmlessly. If W1.4 wants a single flat namespace for provenance, say so.
- **W1.5 (Transformation / Instrument)** — transformations carrying nuisance
  parameters should inherit `Parameterised` and declare through it; a chain's
  parameters compose via `ParameterSet.merge` with the transformation's label
  as the component. Note the name clash the contract avoided: `Bijection`, not
  `Transform`, precisely so `Transformation` stays W1.5's word.
- **W1.6 (Likelihood / NoiseModel)** — GP hyperparameters (amplitude,
  length-scale) are ordinary parameters on a `Parameterised` noise model, with
  `Log` bijections. A latent-GP formulation's latent values are naturally an
  array-valued parameter with a `HierarchicalPrior`; whether that scales as a
  declaration for N ~ 10⁵ is an open question for W1.6.
- **W1.7 (Dataset / FittingProblem)** — `ParameterMapping` is the intended
  substrate for `DatasetCollection`'s joint parameter space:
  `merge({dataset_label: dataset.parameters, ...})`, then `distribute` per
  evaluation. `free_size` is the engine's dimension; `free_labels()` is the
  per-dimension naming; `prior_transform` and `lnprior` are already
  engine-shaped.
- **W1.8 (Results)** — `free_labels()` gives ArviZ coordinates;
  `ParameterSet.to_spec()` is the serialisable declaration to hash into the
  provenance attrs.
- **W1.9 (Lowering)** — the declaration forms needing a lowering row are:
  free/fixed/deferred state, `PriorSpec` per family, `HierarchicalPrior`,
  `Plate` (name + size + array-valued member), `Identity`/`Log`/`Logit`,
  array-valued parameters, buffers, and `Parameter.unit` (which lowers to
  nothing — it is composition-time metadata). `lnprior_unconstrained` is the
  reference semantics the native paths must agree with.
- **W1.10 (Conformance suite)** — candidate rows already implemented and
  testable here: pack/unpack identity, spec round trip, `prior_transform`
  against analytic quantiles, `lnprior` against summed scipy `logpdf`,
  `constrain`∘`unconstrain` identity, the Jacobian identity in §6, and merge
  dimension counting.

## 14. Open questions for review

**Ruled by Peter, 2026-09-01**: the recommendations below stand — question 1
was already resolved by the W1.2 review amendment (core's dependency floor is
numpy/scipy/astropy/stdlib); question 3's flat tie-label namespace stands;
question 4's lone `shared_as` stays allowed; question 5 goes to W1.13 as
written *(closed there, 2026-09-03: `OptionalDependencyError` ratified in
place, with `LoweringError` landing beside it — `lowering.md` §12.4)*;
question 6's `npars` removal is confirmed. **Question 2 (recursive
merge) is expressly kept open**, not closed: Peter can see cases where
hierarchical/nested merging is the natural approach, so W1.7 must treat a
nested `ParameterMapping` as a live design option for `DatasetCollection` —
evaluated on its merits, not dismissed because single-call merge is what
exists today. (W1.5's §14 note that ties crossing merge levels are uncovered
is part of the same question.)

**Closed 2026-09-02**: W1.7 evaluated both topologies (`inference.md` §4) and
Peter **ratified the nested design** — one merge per immediate level, every
mapping retained, values re-distributed level by level; cross-level ties are
verified to collapse correctly. Recursive merge is not adopted. **Later the
same day Peter approved the two follow-ons**: *lossless nesting* — `merge`
accepts a `ParameterMapping` as a component and composes its bindings, so an
outer mapping's bindings are the leaf bindings (`inference.md` §4.6, the
population sketch's recommendation; "lossless, not necessarily recursive") —
and the optional **`Binding.index`** for per-element plate routing through
`distribute` (`hierarchical_population.md` gap H-2). **Both are implemented**
(same day): §8's "Lossless nesting" and "Plate bindings" subsections carry
the contract and executed examples, and §12's items 3 and 4 are updated to
match.

1. **`astropy.units` in `ampere.core`.** `architecture.md` §3–4 says core is
   "numpy/scipy/typing/stdlib"; this module imports `astropy.units` at module
   level. That seems right — astropy is a *required* dependency of the base
   install (not an extra), `DEVELOPMENT_PLAN.md` §4.1 asks parameters to carry
   units, and §3's own namespace listing puts `astropy_compat.py` in core — but
   the architecture spec's wording should be amended to say so explicitly
   rather than leaving core's dependency floor ambiguous.
2. **Should `merge` be recursive?** Limitation 4 above. W1.7 will know whether
   `DatasetCollection`s nest in practice.
3. **Tie label namespace.** A tie label is global and unqualified
   (`"distance"`, not `"shared.distance"`). That reads well and matches how
   users think, but it is a flat global namespace across a whole fitting
   problem. Prefixing (`"shared.distance"`) would be safer and uglier.
4. **`shared_as` on a single component.** Currently allowed: a parameter
   declared shared but merged with nothing else simply takes the bare label as
   its merged name. The alternative — raising on a lone `shared_as`, on the
   theory that it is a typo'd label — would catch a real class of mistake but
   break the legitimate one-component case.
5. **`OptionalDependencyError`'s shape** is pinned in
   `ampere/core/exceptions.py` as `architecture.md` §9 asked; W1.13 should
   ratify it or move it.
6. **`npars` is gone.** The harvested scaffold kept it as a legacy alias.
   Nothing in frozen legacy consumes `ampere.core`, so it has been dropped in
   favour of the unambiguous pair `len(pset)` / `pset.free_size`. Confirm no
   migration path needs it.
