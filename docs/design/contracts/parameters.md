# Ampere v2 — Parameter & Prior Contract (W1.3)

Status: **DRAFT for Peter's review.** Implements `DEVELOPMENT_PLAN.md` §4.1 and
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
...     Parameterised, ParameterSet, Plate, PriorSpec, Tie,
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
| `Parameterised` | The declaration mixin: `register_parameter` / `register_buffer` / `context` |

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

### Which construct to use

Both patterns are expressible, and they serve different data layouts:

- **One component, N members vectorised** (`Plate`): the N objects are fitted
  from one dataset, or from data already stacked. Lowers to a numpyro plate.
- **N components, each with its own dataset** (`HierarchicalPrior` +
  tying): each object gets its own `ParameterSet` with a scalar `theta`, whose
  `HierarchicalPrior` references `mu` and `sigma` — themselves tied across all
  components, so they collapse to one shared pair on merge. This is the shape
  a `DatasetCollection` (W1.7) of per-object datasets naturally takes.

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
   plate raises, because cross-component plate identity needs an
   index-alignment concept that belongs with W1.7's `DatasetCollection`. The
   reference-based pattern in §9 covers the case in practice.
4. **`merge` is not associative.** Merging a `ParameterMapping.merged` again
   would silently drop the first merge's bindings. Merge all components in one
   call. A nested `ParameterMapping` is the obvious extension if W1.7 needs it.
5. **Bijections are declared per parameter, not per set.** A joint bijection
   (a Cholesky factor over several parameters) is not expressible.
6. **No built-in bijection for supports bounded above only.** §6.
7. **Units must match exactly for tying**, not merely be convertible. §8.
8. **Discrete parameters** evaluate and transform correctly (`log_density`
   handles `logpmf`, `ppf` is well defined) but no sampler in ampere currently
   consumes them; treat the support as declared-but-unexercised.

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
written; question 6's `npars` removal is confirmed. **Question 2 (recursive
merge) is expressly kept open**, not closed: Peter can see cases where
hierarchical/nested merging is the natural approach, so W1.7 must treat a
nested `ParameterMapping` as a live design option for `DatasetCollection` —
evaluated on its merits, not dismissed because single-call merge is what
exists today. (W1.5's §14 note that ties crossing merge levels are uncovered
is part of the same question.)

**Closed 2026-09-02**: W1.7 evaluated both topologies (`inference.md` §4) and
Peter **ratified the nested design** — one merge per immediate level, every
mapping retained, values re-distributed level by level; cross-level ties are
verified to collapse correctly. Recursive merge is not adopted and `merge`
itself is unchanged. Two related items remain open at W1.13: *lossless
nesting* (`merge` accepting a `ParameterMapping` as a component and composing
its bindings — `inference.md` §4.6, also the population sketch's
recommendation), and the optional `Binding.index` for per-element plate
routing (`hierarchical_population.md` gap H-2).

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
