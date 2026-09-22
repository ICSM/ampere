# Ampere v2 — Lowering Spec (W1.9)

Status: **frozen at `spec-v1.0`** (the tag created at the W1.13 merge, 2026-09; later changes to lowering rules follow ground rule 9). Implements the lowering half of
`DEVELOPMENT_PLAN.md` §4.1 and operationalises `architecture.md` §5. This is a
**document, not code**: Phase 2 writes the backends, and this spec is the thing
they are written against. Nothing here is executable, so — unlike
`contracts/parameters.md`, whose examples are doctests — it cannot be held
honest by the test suite. W1.10's conformance suite is the forcing function
instead; §11 lists the rows it owes.

Where this document and `DEVELOPMENT_PLAN.md` disagree, the plan wins and this
document is wrong — file it as a decision-log correction.

---

## 0. What lowering is, and what it is not

`ampere.core` declares parameters, priors, plates and buffers in a neutral
vocabulary that knows nothing about torch, jax, numpyro or paramax
(`architecture.md` §3–4). **Lowering** is the one-way translation of those
declarations into a specific backend's objects. It happens once, when a
`FittingProblem` is realised on a backend — never per evaluation.
*(What "realised" means was left unspecified at the freeze; `inference.md`
§10a specifies it — W2.13, 2026-09-07 — as the backend's registered,
differentiable native form of the whole problem.)*

Three targets, and the jax column is really two:

| Target | Distributions & sampling | Module structure & arrays |
|---|---|---|
| **reference** | `scipy.stats` (the declaration *is* the implementation) | plain numpy arrays, no module system |
| **torch** | `torch.distributions` | `torch.nn.Module` (`register_parameter` / `register_buffer`) |
| **jax** | **numpyro** (`numpyro.distributions`, `numpyro.sample`, `numpyro.plate`) | **equinox** `Module` pytrees, with **paramax** wrappers for the trainable/non-trainable split. *(**Amended W2.15**: no paramax. §6.2 ranked `eqx.partition` first and the plan's "jax non-trainable mechanism" row ruled for it on 2026-09-01; W2.5 implemented an `eqx.partition` **filter spec** — `ampere/backends/jax/parameters.py` — and the `jax` extra is `jax`, `numpyro`, `equinox` and nothing else.)* |

Splitting the jax column matters: a jax-backed fit that runs NUTS goes through
numpyro's sample sites, while a jax-backed fit that runs a gradient optimiser
or an SBI simulator goes through equinox pytrees and never constructs a
numpyro model at all. Both consume the same `ParameterSet`, and the rows below
say what each does with it.

**Lowering is not a translation of user code.** It translates *declarations*.
A model's `__call__` is written once against the backend it targets; this spec
governs only the parameter/prior/buffer plumbing around it.

## 1. Standing preconditions — what may reach a backend

These are preconditions, not per-row caveats. A lowering implementation may
assume them and should assert them once at entry rather than defensively
re-checking per parameter.

### 1.1 Lower the *merged* set, never a component set

`ParameterSet.merge` is the resolution step. It qualifies component names
(`"spectrum.temperature"`), collapses tie groups, and returns a
`ParameterMapping` whose `.merged` set is the joint parameter space the
sampler sees. **Lowering consumes `mapping.merged` and nothing else.**

This is `prior_art.md` §6 Tension 1's resolution made operational. Gammapy
tied parameters by Python object identity and was bitten because each
downstream tree-walking consumer had to deduplicate independently; ampere has
*more* such consumers than gammapy did — jax pytree flattening, torch
parameter registration, numpyro site naming, InferenceData provenance — and
this rule is what stops each of them needing its own tie logic. By the time a
backend sees a parameter, tying has already happened.

### 1.2 Tie labels do not lower — there is nothing left to lower

After `merge`, every parameter in `mapping.merged` has `shared_as is None`
(`_collapse` constructs it that way). A tie is not a runtime object, a
constraint, or a shared reference: it is a **compile-time collapse**. Two
sites that were tied are, downstream of `merge`, literally one `Parameter`
with one name, one prior, one flat-vector slice, one torch `nn.Parameter`, one
numpyro sample site.

The bindings that record *where the collapsed value is consumed* live in
`ParameterMapping.bindings`, and they are routing metadata for
`mapping.distribute(theta)` — the step that builds each component's keyword
arguments. Routing is backend-neutral and happens above the backend boundary.
**No backend should ever see, emit, or name a tie label.** If a lowering
implementation finds itself wanting to, it has been handed an unmerged set;
that is a bug in the caller, and §1.5's entry assertion should have caught it.

### 1.3 Deferred parameters must not reach lowering

A deferred parameter (`shared_as` set, no prior of its own — `parameters.md`
§3) is *by construction* not evaluable: `ParameterSet` refuses `lnprior`,
`prior_transform` and `sample` on an unresolved set. `merge` resolves it, and
a tie group in which no site supplies a prior raises `TyingError` rather than
producing a prior-less merged parameter.

Therefore: **a merged set contains only free and fixed parameters.** Lowering
has exactly two states to handle, not three. A lowering implementation that
encounters `parameter.is_deferred` has been given an unmerged set and must
raise, not invent an improper prior.

### 1.4 Units lower to nothing

`Parameter.unit` is composition-time metadata. Values are plain numbers in the
declared unit before they reach any backend (`parameters.md` §5, and the
units-in-the-hot-loop trap in `DEVELOPMENT_PLAN.md` §7). No backend
constructs a `Quantity`, converts a unit, or carries a unit on a tensor.

The unit is still *recorded*: it belongs in the provenance attrs and in the
ArviZ coordinate metadata (W1.8), because a posterior in kelvin and a
posterior in log-kelvin are different results. But that is a results-side
concern, and it travels through `to_spec()`, not through a backend.

### 1.5 Opaque priors never reach lowering

`describe_prior` raises for anything that is not a frozen `scipy.stats`
distribution (`parameters.md` §4). That is the gate: a duck-typed prior
evaluates on the reference path and is refused at serialisation *and* at
lowering, with a message naming the parameter. A backend therefore never has
to ask "can I describe this?" — it works from `PriorSpec` throughout.

### 1.6 The entry assertion

One function, run once per lowering, checking all of the above:

```
lower(pset) requires:
    pset.is_resolved                       # §1.3 — no deferred parameters
    all(p.shared_as is None for p in pset)  # §1.2 — merge has happened
    pset.to_spec() succeeds                 # §1.5 — every prior is describable
```

Failing loudly here, once, is worth more than three subtly different errors
raised from three backends halfway through construction.

## 2. The direction convention

**Read this section before writing a single line of transform code in any
backend.** A reversed Jacobian is the named silent-bug class this contract
exists to prevent (`parameters.md` §6, §11): it does not crash, it does not
produce NaNs, it produces a subtly wrong posterior that looks plausible.

Ampere deliberately does not use `forward`/`inverse`, because both target
libraries use `forward` for the direction a statistician often reads the other
way. The mapping is:

| Concept | ampere (`Bijection`) | torch (`torch.distributions.transforms`) | numpyro (`numpyro.distributions.transforms`) |
|---|---|---|---|
| unconstrained → constrained (real line → support) | `constrain(y)` | `t(y)` / `t.forward` | `t(y)` / `t.__call__` |
| constrained → unconstrained | `unconstrain(x)` | `t.inv(x)` | `t.inv(x)` |
| the transform obtained from a support/constraint | `default_bijection_for(prior)` | `biject_to(constraint)` | `biject_to(constraint)` |
| log-abs-det Jacobian of `d(constrain)/dy` | `log_abs_det_jacobian(y)` — **one argument, the unconstrained point** | `log_abs_det_jacobian(x, y)` — **two arguments**, `x` = input = *unconstrained* | `log_abs_det_jacobian(x, y, intermediates=None)` — **two arguments**, same order |
| the registry to use | `default_bijection_for` | `biject_to` — **not** `transform_to` (§2 b′) | `biject_to` — **not** `transform_to` |

Three consequences, each of which is a place an implementation goes wrong:

**(a) `constrain` is `forward`, not `inverse`.** `biject_to(constraint)` in
both libraries returns a transform whose *forward* direction maps the
unconstrained real line **onto** the constrained support. That is the same
direction as ampere's `constrain`. So:

```
ampere:  x = bijection.constrain(y)
torch:   x = biject_to(constraint)(y)
numpyro: x = biject_to(constraint)(y)
```

and the inverse direction is `bijection.unconstrain(x)` ↔ `t.inv(x)`.

**(b) The two-argument Jacobian: the first argument is the transform's
*input*, i.e. the unconstrained point.** Both libraries' signature is
`log_abs_det_jacobian(x, y)` where `x` is the forward transform's input and
`y` its output, and the return value is `log|dy/dx|` — the determinant of the
*forward* direction. Both arguments are taken because different transforms
need different ones: torch's `ExpTransform` returns `x`, its `AffineTransform`
returns `log|scale|` and uses neither. For a transform obtained from
`biject_to`, the input is the unconstrained value and the output is the
constrained one. So the correct call, and the only one that matches ampere's
reference semantics, is:

```
ampere:  ladj = bijection.log_abs_det_jacobian(y_unconstrained)
torch:   ladj = t.log_abs_det_jacobian(y_unconstrained, x_constrained)
numpyro: ladj = t.log_abs_det_jacobian(y_unconstrained, x_constrained)
```

**Beware the variable-name collision, which is the reason to write this out
in full.** Ampere's convention (following the statistics literature, and
`parameters.md` §6) is `x` = constrained, `y` = unconstrained. torch's and
numpyro's transform docstrings use the opposite: `x` = input = unconstrained,
`y` = output = constrained. The two libraries' `x` is ampere's `y`. Copying a
call out of the torch documentation into ampere-named variables — or the
reverse — therefore transposes the arguments while *looking* correct in both
places. This is not a hypothetical: the identifiers match, the shapes match,
the dtypes match, and the result is merely wrong.

Passing the arguments the other way round is the silent bug. **A conformance
row must pin this numerically at a point where the two arguments differ**
(§11), because no type checker will, and a test at a symmetric point will not
either.

**(b′) Use `biject_to`, never `transform_to`.** Both registries map from the
real line to a constraint, in the same direction, and for simple constraints
they often return the same transform — so a substitution passes casual
inspection. They differ in guarantees: `biject_to` returns a transform
guaranteed `bijective = True` and guaranteed to implement
`log_abs_det_jacobian`; `transform_to` guarantees neither, exists for
unconstrained *optimisation*, and may overparameterise. The canonical
divergence is the simplex: `transform_to` gives a softmax (cheap,
non-bijective), `biject_to` gives stick-breaking (bijective, has a Jacobian).
Ampere needs the Jacobian on every path that reaches
`lnprior_unconstrained`, so the answer is always `biject_to`.

**(c) The reference implementation is the oracle.**
`ParameterSet.lnprior_unconstrained(y)` (see `ampere/core/parameter.py`) is
the definition every backend must agree with:

```
lnprior_unconstrained(y) = lnprior(constrain(y)) + Σ_i log|d constrain_i / dy_i|
```

with the Jacobian evaluated **at the unconstrained point `y`**, summed over
every element of every free parameter, and short-circuiting to `-inf` if any
term is non-finite. The native paths must reproduce this to tolerance, not
merely produce "a valid unconstrained density" — a density differing by a
non-constant function of `y` is a different posterior, and one differing by a
constant is fine for MCMC but wrong for evidence.

**(d) Do not double-count.** `TransformedDistribution` (both libraries)
already applies the change of variables internally: torch's implementation
computes `log_prob(y) = base.log_prob(x) − Σ log|dy/dx|`, subtracting the
*forward* determinant because the density is being evaluated in the output
space. Composing that with an additional ampere-computed Jacobian double-counts
the correction, and — because the two terms have the same sign structure —
produces a density that is smooth, finite and wrong by a factor that varies
with position.

The rule is that **exactly one mechanism owns the correction**, and which one
depends on why the transform is there:

- a transform used to *express a family the target lacks* (§3.4's
  `loguniform`) belongs inside a `TransformedDistribution`, which owns its own
  Jacobian; ampere adds nothing;
- a transform used to *unconstrain a parameter for a sampler* (§4) is applied
  by the sampling layer, and its Jacobian is the term
  `lnprior_unconstrained` adds.

These compose without conflict — the first is part of the prior's definition,
the second part of the change of variables — but an implementation that routes
both through the same code path will apply one of them twice.

## 3. Distribution mapping

### 3.1 The source side is canonical and keyword-only

`describe_prior` emits a canonical, **keyword-only** `PriorSpec`: scipy accepts
the same freezing positionally or by keyword, and positional arguments are
already mapped onto their names (shape names, then `loc`, then `scale`) before
lowering sees them. The table therefore has one form per family to translate,
and `PriorSpec.args` is always empty for a described prior.

The keyword names below are the verified output of `describe_prior` on scipy
1.17, not a reconstruction from documentation.

### 3.2 The table

Required families first — these are the ones the Phase 1 specs actually use.

| Family (`PriorSpec.family`) | canonical `kwds` | reference (scipy) | torch.distributions | numpyro |
|---|---|---|---|---|
| `norm` | `loc`, `scale` | `norm(loc, scale)` — identity | `Normal(loc=loc, scale=scale)` | `Normal(loc=loc, scale=scale)` |
| `uniform` | `loc`, `scale` | `uniform(loc, scale)` — identity | `Uniform(low=loc, high=loc+scale)` | `Uniform(low=loc, high=loc+scale)` |
| `halfnorm` | `loc`, `scale` | `halfnorm(loc, scale)` — identity | `HalfNormal(scale=scale)` when `loc == 0`; otherwise shift (§3.3) | `HalfNormal(scale=scale)` when `loc == 0`; otherwise shift (§3.3) |
| `loguniform` | `a`, `b` | `loguniform(a, b)` — identity | **no native equivalent**; `TransformedDistribution(Uniform(log a, log b), ExpTransform())` | `LogUniform(low=a, high=b)` |
| `poisson` | `mu` | `poisson(mu)` — identity | `Poisson(rate=mu)` | `Poisson(rate=mu)` — **but see §3.5, discrete parameters do not lower to a gradient path** |
| `truncnorm` | `a`, `b`, `loc`, `scale` | `truncnorm(a, b, loc, scale)` — identity | **unsupported — raises at lowering** (§3.4); absent from torch *and* from pyro | `TruncatedNormal(loc, scale, low=loc+a*scale, high=loc+b*scale)` — **`low`/`high` are keyword-only**; note the standardisation, §3.3 |

Second tier — not used by the Phase 1 specs, but common enough that leaving
them to be rediscovered per backend invites divergence. Same discipline, same
obligations. These rows were verified less thoroughly than the first tier
(§13): `Gamma` and `Beta` signatures were read in both libraries' sources, but
`Exponential`'s was not, so confirm `rate = 1/scale` against the library before
relying on it.

| Family | canonical `kwds` | reference | torch | numpyro |
|---|---|---|---|---|
| `lognorm` | `s`, `loc`, `scale` | identity | `LogNormal(loc=log(scale), scale=s)`, **`loc` must be 0** | `LogNormal(loc=log(scale), scale=s)`, **`loc` must be 0** |
| `expon` | `loc`, `scale` | identity | `Exponential(rate=1/scale)`, **`loc` must be 0** — otherwise shift | `Exponential(rate=1/scale)`, **`loc` must be 0** — otherwise shift |
| `gamma` | `a`, `loc`, `scale` | identity | `Gamma(concentration=a, rate=1/scale)`, **`loc` must be 0** | `Gamma(concentration=a, rate=1/scale)`, **`loc` must be 0** |
| `beta` | `a`, `b`, `loc`, `scale` | identity | `Beta(concentration1=a, concentration0=b)`, **`loc`=0, `scale`=1** | `Beta(concentration1=a, concentration0=b)`, **`loc`=0, `scale`=1** |
| `halfcauchy` | `loc`, `scale` | `halfcauchy(loc, scale)` — identity | `HalfCauchy(scale=scale)` when `loc == 0`; otherwise shift (§3.3) | `HalfCauchy(scale=scale)` when `loc == 0`; otherwise shift (§3.3) |

*(**Amended W2.15**, 2026-09-08: "must be 0" describes what the *target
library's* constructor takes, and both backends now do what §3.3 rule 1 says
to do about it rather than refusing. W2.4 and W2.5 compose the affine map away
for `lognorm`, `expon`, `gamma` and `beta` alike — a shifted `gamma` lowers,
and a `beta` on `[loc, loc+scale]` lowers to a `Beta` under an affine
transform whose support is the interval scipy's has. What is still refused is
a family with no exact construction from the target's primitives, which is
rule 2. The rows are unchanged as statements about the libraries;
`ampere/backends/torch/lowering.py` and
`ampere/backends/jax/distributions.py` are the implementations.)*

*(**Added W5.25**: `halfcauchy` is the horseshoe's global scale under both of
`shrinkage_horseshoe`'s tails (§3.2.1; *Amended W5.27* — renamed, the old
name kept as a deprecated alias, because it promised Piironen & Vehtari's
slab and the default tail is a gamma-tailed sibling of it) and had no row on
either modern backend, so the recommended prior for any `Sum` of noise terms was
reference-only on NUTS. `torch.distributions.HalfCauchy` and
`numpyro.distributions.HalfCauchy` are both exact, so this is §3.4's fallback
used as intended, on `halfnorm`'s own pattern: native at `loc == 0`, the
§3.3 shift otherwise — for `torch`, `_shifted`'s affine route; for `numpyro`,
`TruncatedCauchy(loc, scale, low=loc)`, preferred over an affine
`TransformedDistribution` for the same reason `_halfnorm` is, in
`ampere/backends/torch/lowering.py` and `ampere/backends/jax/distributions.py`
respectively.)*

### 3.2.1 Hierarchical priors, and the parameterisation NUTS wants (*Added W5.8*)

Everything above is a *frozen* distribution, whose arguments are numbers. A
`HierarchicalPrior` (§8) supplies one or more of them from another parameter's
value instead, and lowers through exactly the same table: the family name is
the same name, the referenced value is what the lowered constructor receives,
and the ordering §8 imposes is what guarantees it has a value by then. Nothing
in §3.2 changes; this subsection exists because *which* of two equivalent
declarations is written changes how well a gradient sampler explores the
result, and that is a lowering concern rather than a modelling one.

**The funnel.** `θ_k | s ~ Normal(0, s)` with `s` itself sampled has a
posterior whose width in `θ_k` is proportional to `s` — a funnel, narrowing to
a point at `s = 0`. NUTS chooses one step size for the whole geometry, so a
step that works in the mouth of the funnel overshoots its neck and a step that
works in the neck cannot cross the mouth: the reported symptom is divergences,
sometimes many, and the silent symptom is a chain that never visits small `s`.
The **non-centred** parameterisation is the same prior written so that the
sampled coordinates are independent of each other:

| centred (what is modelled) | non-centred (what is sampled) |
|---|---|
| `s ~ HalfNormal(σ)`, `θ_k | s ~ Normal(0, s)` | `s ~ HalfNormal(σ)`, `z_k ~ Normal(0, 1)`, `θ_k = s · z_k` |

The right-hand column has no funnel, because `z_k` and `s` are a priori
independent, and it lowers to two ordinary rows of §3.2 (`halfnorm` and `norm`)
with a multiplication in between.

**Where the multiplication lives.** This contract has no deterministic node —
`parameters.md` §9 — so `θ_k = s · z_k` is formed by whatever consumes the
knot variables, not declared. `WarpedKernel` does exactly that and it is the
reason `non_centred=True` is its default: the kernel declares `z_k ~ Normal(0,
1)` and forms `s · z_k` itself, so the sampler sees the right-hand column while
the model is the left-hand one. `non_centred=False` declares the left-hand
column directly with `HierarchicalPrior`, which is what a strongly identified
warp can afford and what the plan names.

**The horseshoe's chain now lowers, in full.** (*Amended W5.25*)
`shrinkage_horseshoe` (`parameters.md` §9) is three hierarchical levels, and
its default `tail="regularised"` is chosen so that every one of them is in the
table above: `halfcauchy` for the global scale, `gamma` for each local scale,
`halfnorm` for each amplitude. `halfcauchy` was in neither backend's
registration until W5.25 added one `register_lowering("halfcauchy", backend,
…)` row per backend, on `halfnorm`'s own pattern — native at `loc == 0`, the
§3.3 shift otherwise, which is §3.4's fallback rule used as intended
(`torch.distributions.HalfCauchy` and `numpyro.distributions.HalfCauchy` are
both exact). `tail="cauchy"`'s local level is the same family and needed no
row of its own as a consequence — it lowers on both backends too. On jax this
closed the chain completely on its own, because a `HierarchicalPrior`'s
dispatch *is* the flat table above (§5): any family §3.2 has, a hierarchical
reference to it gets for free. Torch's dispatch is a second, per-family
registry (`_HIERARCHICAL_BUILDERS` in `ampere/backends/torch/lowering.py`),
which W5.25 extended with both `halfcauchy` and `gamma` — the latter found
only once the former was in place and the default tail's local level still
would not reach NUTS on torch. `gamma`'s shape argument (`a`) is fixed rather
than referenced in every declaration this contract writes, so it lowers as an
ordinary tensor alongside `loc`/`scale`, exactly as §8's "arguments are other
parameters' values" already allows; what it cannot do is acquire a bijection
automatically (`_default_bijection_for_hierarchical` refuses any family with a
shape argument, referenced or not), which is why `shrinkage_horseshoe`
supplies `bijection=Log()` explicitly rather than relying on the default. All
three levels of both tails now reach NUTS on both backends. The horseshoe's
own funnel is the one described above, one level deeper, and it is why
`tests/m2`'s shrinkage rows run at a longer emcee budget than the rest of
that suite and pin a margin a third below what they measure.

Spike-and-slab — the other classical sparsity prior, and the one a reader may
expect here — stays out, for the reason `horizon_notes.md` §1 gives and §12 Q1
rules: its inclusion indicator is **discrete**, so it has no default bijection
(`parameters.md` §6) and no gradient path. A continuous shrinkage prior is what
this contract can lower, which is why it is the one `parameters.md` §9 ships.

### 3.3 The parametrisation traps, spelled out

Every row above that says "must be 0" is a trap, and they share one cause:
**scipy gives every continuous family a `loc`/`scale` shift-and-stretch that
the target libraries' equivalents do not have.** `scipy.stats.halfnorm(3, 2)`
is supported on `[3, ∞)`; `torch.distributions.HalfNormal(2)` is supported on
`[0, ∞)` and has no way to say otherwise. Lowering the former to the latter
silently moves the prior's support by 3 units.

**Rule: a non-default `loc` (or, where noted, `scale`) that the target family
cannot express is not silently dropped.** The lowering has two admissible
responses, in order of preference:

1. **Compose it away.** Where an affine shift is exactly equivalent — `expon`,
   `halfnorm`, `gamma`, and any other family whose `loc` is a pure
   translation of the support — lower to
   `TransformedDistribution(base, AffineTransform(loc=loc, scale=1))`. This is
   exact, not an approximation, and it is the preferred route because it keeps
   the prior a first-class distribution object on the backend.
2. **Raise at lowering** with a message naming the parameter, the family, and
   the offending argument — never warn, never drop.

For the shifted half-normal specifically, both backends have an exact route
and they are different, which is worth writing down rather than rediscovering:

| | shifted `halfnorm(loc=l, scale=s)` |
|---|---|
| torch | `TransformedDistribution(HalfNormal(s), AffineTransform(loc=l, scale=1))` |
| numpyro | `TruncatedNormal(l, s, low=l)` — a normal centred at `l` and truncated below at `l` **is** a half-normal shifted to `l`, exactly (the truncation renormalises by 2, which is the half-normal's factor) |

The numpyro route is preferable there because it produces a distribution whose
`support` is `greater_than(l)`, so `biject_to` yields the right bijection
automatically (§4); the `TransformedDistribution` route requires the bijection
to be reasoned about separately.

`truncnorm` deserves its own paragraph, because it is the one row where the
numbers change rather than merely move. **scipy's `a` and `b` are in
standardised units**: the support is `[loc + a·scale, loc + b·scale]`, not
`[a, b]`. Verified: `scipy.stats.truncnorm(a=-1, b=2, loc=5, scale=3)` has
support `(2, 11)`. numpyro's `TruncatedNormal` takes `low`/`high` in the
**unstandardised** space. Passing scipy's `a`/`b` straight through as
`low`/`high` produces a prior truncated at `(-1, 2)` around a mean of 5 — a
prior that excludes its own mode, and which will read as a badly behaved
sampler rather than as a translation bug. The conversion in the table is
mandatory and belongs in a conformance row.

Two further facts about numpyro's truncation machinery constrain how far this
row generalises. `TruncatedNormal` is a **factory function, not a class**, and
its `low`/`high` are **keyword-only** — a positional call fails. And the
underlying `TruncatedDistribution` accepts only a fixed set of real-supported
base families (`Normal`, `Cauchy`, `Laplace`, `Logistic`, `SoftLaplace`,
`StudentT`). So "truncate the base distribution" is not a general escape hatch
for arbitrary truncated families; outside that set, §3.4's rule applies and
lowering raises.

`lognorm` is the second numeric conversion: scipy parametrises by the shape
`s` (the log-scale standard deviation) and `scale = exp(mu)`; both libraries'
`LogNormal` take `loc = mu` directly. So `loc_target = log(scale_scipy)`.
Passing `scale` through unchanged is wrong by an exponential.

`uniform` is the third, and the most likely to be got right by accident and
wrong under maintenance: scipy's second argument is a **width**, not an upper
bound. `high = loc + scale`.

### 3.4 Families with no native equivalent — the fallback rule

There is no torch `LogUniform` and no torch `TruncatedNormal`, and there will
be families beyond this table that neither target implements. Silence is not
an option — the whole point of the `describe_prior` gate (§1.5) is that an
un-lowerable prior fails early with a name attached, rather than travelling to
a backend and failing there or, worse, being quietly approximated.

The rule, in order:

1. **Exact construction from primitives.** If the family can be built exactly
   out of the target's own pieces — `TransformedDistribution`, an affine or
   exponential transform, a truncation wrapper — do that, and register it in
   this table so it is a decision on the record rather than a per-backend
   improvisation. torch's `loguniform` row is exactly this case.
2. **Explicitly unsupported → raise at lowering.** If no exact construction
   exists, lowering raises a typed error naming the family, the parameter, and
   the backend, and says what the user's options are: change the prior, or run
   on a backend that has it.
3. **Never silently approximate.** A "close enough" substitution (a
   `TruncatedNormal` approximated by a `Normal`, a `LogUniform` approximated
   by a wide `LogNormal`) is forbidden. It changes the posterior in a way no
   diagnostic will attribute to the translation layer.
4. **Fall back to the reference path only when the user asked for it.** The
   reference backend can always evaluate a described prior, because scipy *is*
   the declaration. But falling back automatically would silently drop the
   gradients or the speed the user chose a native backend for. A fallback is
   therefore opt-in and per-run, never a default — the same "opt-in, never
   silent" ruling `DEVELOPMENT_PLAN.md` §2 recorded for curated translation.

The error type should be one shared class across backends, so tests and
tooling assert on it uniformly, in the same spirit as `architecture.md` §4's
rule 3. **Landed at the freeze** (ruled 2026-09-03, §12.4):
`ampere.core.exceptions.LoweringError(AmpereError)` — deliberately **not**
under `ContractError`, because a family torch does not implement is a
capability gap in the backend, not a malformed declaration by the user;
conflating the two would make "your prior is invalid" and "this backend
cannot express your valid prior" indistinguishable to a caller catching the
exception. It carries the family, the parameter name and the backend as
fields (plus a free-text `detail`), and its message names the three options:
change the prior, register your own lowering (§12.8's hook, Phase 2), or
run on a backend that has the family. `OptionalDependencyError` was
ratified in place in the same pass (`parameters.md` §14 Q5), and
`ResultsError` had already moved to core (`results.md` §15 R5).

### 3.5 Discrete families

`poisson` lowers cleanly as a *distribution* on all three targets. It does
**not** lower to a gradient-based sampling path: a discrete parameter has no
meaningful unconstraining bijection, and NUTS/HMC cannot sample it.

*(Amended at the freeze: the upstream fix this section asked for landed —
ruled 2026-09-03, §12.1. `default_bijection_for(scipy.stats.poisson(3.0))`
now raises `CapabilityError`, a typed capability refusal that is
deliberately not a malformed-declaration error; the refusal lives only
there, and discreteness is queryable as `describe_prior(prior).discrete`.
The paragraph below records the pre-freeze state for the history.)*

This was worth stating because ampere's reference implementation used to
hand you one — `default_bijection_for(scipy.stats.poisson(3.0))` returned
`Log(lower=0.0)`, inferred from the support `[0, ∞)` with no regard for
discreteness. `parameters.md` §12.8 declares discrete parameters
"declared-but-unexercised", so that was consistent with the contract rather
than in conflict with it, but the lowering layer is where it would have
bitten.

**Rule: lowering a discrete parameter onto a gradient-requiring path (torch
HMC, numpyro NUTS, any `lnprior_unconstrained` consumer) raises** — now
enforced at the source, in `default_bijection_for` itself. Lowering it
onto a non-gradient path (a reference-backend nested sampler via
`prior_transform`, an SBI simulator, §3.5's lowering-as-distribution) is
fine and untouched, which is exactly the door the ruling keeps ajar.

### 3.6 `icdf` availability, and what it costs the nested-sampling path

`ParameterSet.prior_transform` — the unit-hypercube map nested samplers need —
is an inverse-CDF composition. On the reference path every scipy family has
`ppf`, so it always works. On the native backends it is **partial**, and the
gap is a real constraint rather than a detail:

| family | torch `icdf` | numpyro `icdf` |
|---|---|---|
| `norm`, `uniform`, `halfnorm` | implemented | implemented |
| `lognorm`, `loguniform` | inherited from the base via `TransformedDistribution` | inherited the same way (both *are* `TransformedDistribution`s of `Normal`/`Uniform`) |
| `gamma` | **absent** (`cdf` exists, `icdf` does not) | present, but **needs TensorFlow Probability** — see below |
| `beta` | **absent** (neither `cdf` nor `icdf`) | present, but **needs TensorFlow Probability** — see below |
| `poisson` | **absent** | **absent** (`cdf` only) |

Both libraries declare `icdf` on the base `Distribution` and raise
`NotImplementedError` unless a family overrides it, so the failure is loud
rather than silent.

***Amended W2.13*** *(2026-09-07; the "Realisation surface (W2.13)" row in
`DEVELOPMENT_PLAN.md` §2, spec correction 1).* The original table read
"implemented" for numpyro's `Gamma` and `Beta`, and that is true of the
*method* and false of the *environment*. Both are implemented by delegating to
`gammaincinv`/`betaincinv`, which numpyro does not have and imports from
TensorFlow Probability on use; with the `jax` extra as ampere ships it —
jax, numpyro, equinox, no TFP — calling either raises `ImportError: Please
install 'tensorflow_probability>=0.18' for gammaincinv` (measured on numpyro
0.21.0). So the asymmetry the original paragraph drew between the two backends
is not there in a default install: **neither backend has a native `icdf` for
`gamma` or `beta`**, and both take §3.6's sanctioned reference fallback for
them.

Two consequences, both deliberate. A backend must decide `icdf` availability
by **asking the environment**, not by consulting a table — `has_native_icdf`
in `ampere.backends.jax.distributions` probes the call — because whether these
two rows are available depends on a dependency ampere neither requires nor
forbids. And a user who installs TFP alongside the `jax` extra silently gains
native coverage for them, which is exactly why the fallback is warned about
rather than silent: the run says which path it took.

**Rule: a nested-sampling run whose priors lack a native `icdf` evaluates its
prior transform on the reference path.** This is legitimate and cheap — the
prior transform is called once per proposal on a vector of size `free_size`,
not inside the likelihood's hot loop, and it is exactly the deterministic,
backend-independent map §9.3 wants it to be. It is *not* the silent fallback
§3.4 forbids, because nothing about the posterior changes: `prior_transform`
has one mathematical definition and the reference path computes it exactly.
The distinction worth holding onto is that §3.4 forbids substituting a
*different prior*, not computing the *same* quantity somewhere else.

**Ruled (2026-09-03): sanctioned on these terms, but loud.** Taking the
fallback emits a warning naming the affected families and the backend, so
a nominally torch-backed run that computes its prior transform in numpy is
never a surprise discovered later; a `strict` option turns the warning
into a raise.

**The flag's home is pinned (W1.13, at the freeze): the fitting problem's
existing `strict` toggle** — `FittingProblem(strict=...)`, the run-level
control `inference.md` §11 defines — not a new per-lowering knob. One flag,
one meaning: `strict=True` already says "I would rather fail than have
anything smoothed over" (exceptions propagate instead of becoming −inf
with a recorded reason), and refusing to mix a numpy prior transform into
a nominally native run is the same preference at lowering time. Phase 2
implements exactly this: when a backend lowers `prior_transform` and a
family lacks a native `icdf`, it consults the problem's `strict` — `False`
(the default) takes the reference fallback and warns **once per run**,
naming the families and the backend; `True` raises a `LoweringError`
naming them instead. The warning is per run, not per call, because the
fallback decision is made once at lowering time, before any sampling.

## 4. Bijection mapping

`parameters.md` §6 fixes the inference rule — real line → `Identity`,
half-line → `Log`, bounded interval → `Logit`, bounded-above-only → raise —
and states that it is numpyro's `biject_to` rule, chosen so the reference and
native paths agree by construction rather than by coincidence. This table is
the other half of that promise.

| ampere | maths | reference | torch | numpyro |
|---|---|---|---|---|
| `Identity()` | `x = y` | `Identity` | `biject_to(constraints.real)` | `IdentityTransform()` |
| `Log(lower=0)` | `x = exp(y)` | `Log(lower=0.0)` | `ExpTransform()` | `ExpTransform()` — what `biject_to(constraints.positive)` returns, bare |
| `Log(lower=l)`, `l ≠ 0` | `x = l + exp(y)` | `Log` | `ComposeTransform([ExpTransform(), AffineTransform(loc=l, scale=1)])` | same, with `AffineTransform(l, 1, domain=constraints.positive)` |
| `Logit(lower=0, upper=1)` | `x = σ(y)` | `Logit(0, 1)` | `SigmoidTransform()` | `SigmoidTransform()` — `biject_to(constraints.unit_interval)` |
| `Logit(lower=l, upper=h)` | `x = l + (h−l)·σ(y)` | `Logit` | `ComposeTransform([SigmoidTransform(), AffineTransform(loc=l, scale=h−l)])` | same, with `AffineTransform(l, h−l, domain=constraints.unit_interval)` |
| bounded above only | — | raises at declaration | (never reaches lowering) | (never reaches lowering) |
| user-supplied custom `Bijection` | arbitrary | works (duck-typed) | **unsupported — raises at lowering** | **unsupported — raises at lowering** |

Notes that decide implementations:

- **Compose order.** In both libraries `ComposeTransform([f, g])` applies `f`
  first, then `g`. So the exponential/sigmoid comes first and the affine
  shift second, matching `x = l + exp(y)`. Getting this backwards gives
  `exp(l + y)` — a different, and still perfectly finite, prior.
- **`positive` and `greater_than(0)` are different constraint types**, and
  numpyro registers them separately: a distribution whose support is
  `constraints.positive` (`HalfNormal`, `Gamma`, `LogNormal`) gets a **bare**
  `ExpTransform`, while one whose support is `constraints.greater_than(l)` (a
  left-truncated normal) gets the composed form. At `l = 0` the two are
  mathematically identical, so this is not a correctness issue — but it means
  the *object* you get back from `biject_to` differs between two priors with
  the same support, and any code that pattern-matches on transform type rather
  than on behaviour will be surprised. Do not pattern-match on transform type.
- **numpyro's composed forms carry an explicit `domain=`** on the inner
  `AffineTransform` (`constraints.positive` or `constraints.unit_interval`).
  This matters only if a lowering builds the chain by hand — which the next
  note says it should not.

  ***Amended W2.13*** *(spec correction 3).* Read this as "numpyro's own
  composed forms state a `domain=` and a lowering that builds one by hand must
  too", because `AffineTransform`'s constructor **defaults `domain` to
  `constraints.real`**. An `AffineTransform(loc, scale)` written without it is
  therefore declared to accept the whole line, and the composed
  distribution's `support` — which is derived from the outermost transform's
  codomain, itself derived from that domain — comes back wrong. W2.5 hit this
  building the shifted-family rows.
- **Prefer `biject_to(dist.support)` over reconstructing the transform.** When
  a parameter's bijection is the *inferred default* (`Parameter.bijection is
  None`), the correct lowering is to ask the target library for
  `biject_to(lowered_distribution.support)` rather than to build the compose
  chain by hand from the table. The table then documents what you will get and
  serves as the cross-check; the library stays the authority on its own
  constraint registry, and the two agree because ampere adopted its rule.
  Build by hand **only** when the user declared a bijection explicitly, in
  which case `biject_to` would give the wrong answer.

  ***Amended W2.13*** *(spec correction 2, from W2.4).* This advice holds only
  where `lowered.support` is the distribution's *real* support, and on torch
  it often is not: `torch.distributions.TransformedDistribution.support` is
  **the last transform's codomain**, not the support of the composition. For
  the two rows §3.2/§3.3 need it — `loguniform`, and any shifted family
  lowered as an affine composition — that codomain is strictly larger than the
  distribution's actual support, so `biject_to` of it yields a bijection onto
  the wrong set and an unconstrained sample can land where the density is
  zero. A torch lowering must therefore **declare the support it actually
  has** rather than read it off the composed object; W2.4 does this with a
  `TransformedDistribution` subclass overriding `support`, and the advice
  above then applies to the declared value. numpyro's composed forms do not
  have the problem, because they carry the explicit `domain=` the note above
  now insists on.
- **The rule does agree, and `loguniform` is the check worth citing.** Ampere
  infers `Logit(lower=a, upper=b)` for `loguniform(a, b)` — from its support,
  which is the bounded interval `[a, b]`, not from its log-shaped density.
  numpyro's `LogUniform` declares `support = interval(low, high)`, so
  `biject_to` independently produces sigmoid-plus-affine over the same
  interval. The two paths arrive at the same bijection without either being
  told about the other, which is the property `parameters.md` §6 was after.
  (One might argue a log-space bijection suits a log-uniform prior better; that
  is a question about ampere's inference rule, not about lowering, and the rule
  is fixed by the contract.)
- **A custom `Bijection` is a reference-path-only feature**, exactly parallel
  to a duck-typed prior (§1.5). It satisfies ampere's protocol with numpy
  operations that a torch tensor or a jax tracer will not accept. The
  extension point is the **hardened registration hook** ruled 2026-09-03
  (§12.8): `register_lowering(family, backend, constructor)` for prior
  families and the analogous per-backend slot for a custom `Bijection` — no
  silent overwrites (`override=True` required), user-registered rows stamped
  in provenance, an opt-in conformance battery for registrants, and
  constructors that must return trace-pure objects since the registry
  resolves before any tracing. The plumbing is Phase 2's, beside the
  backends that consume it.

## 5. The declaration-form lowering table

This is the table W1.9's acceptance criterion asks for: every declaration form
`parameters.md` §13 enumerates, with a row per backend or an explicit
"unsupported".

| §13 declaration form | reference | torch | jax |
|---|---|---|---|
| **free** parameter | an entry in the flat float64 array; `scipy` prior evaluated directly | `torch.nn.Parameter(tensor, requires_grad=True)` registered via `register_parameter`; prior as a `torch.distributions` object per §3 | equinox: an ordinary array field on the `Module`, **in** the trainable partition. numpyro: one `numpyro.sample(name, fn)` site |
| **fixed** parameter | not in the flat vector; injected by `unpack` at its declared value | a **buffer**, not a `Parameter` — `register_buffer(name, tensor)`, `requires_grad=False`, no prior object constructed | equinox: an array field **excluded** from the trainable partition (§6.2). numpyro: **not a sample site** — a plain constant closed over, or `numpyro.deterministic` if it must appear in the trace |
| **deferred** parameter | — | — | — (all three: **must not reach lowering**; `merge` resolves it first, §1.3) |
| **`PriorSpec`** per family | identity — scipy is the declaration | §3.2, with the `loc`/`scale` rules of §3.3 and the fallback rule of §3.4 | §3.2 via numpyro; same rules |
| **`HierarchicalPrior`** | evaluated in topological order; `bind()` freezes the family against already-resolved values | distribution object constructed **per evaluation** from the referenced parameter tensors (its arguments are tensors, so it cannot be built once at lowering) | numpyro: the natural form — `dist.Normal(mu, sigma)` where `mu`, `sigma` are earlier sample sites in the same model function. Site order follows the same topological order |
| **`Plate`** (name, size, array-valued member) | already expanded into ordinary parameters by `Plate.expand()`; the member is one array-valued parameter of shape `(size,) + member.shape` | a single `nn.Parameter` of shape `(size,) + member.shape`; the prior is the member distribution **batched** over the leading axis | numpyro: `with numpyro.plate(plate_name, size, dim=-1): numpyro.sample(member_name, fn)` — one site producing `size` draws, which is exactly why the contract chose the array-valued shape |
| **`Identity`** bijection | `Identity` | §4 | §4 (numpyro `biject_to`) |
| **`Log`** bijection | `Log` | §4 | §4 |
| **`Logit`** bijection | `Logit` | §4 | §4 |
| **array-valued** parameter | contiguous C-order block of the flat vector; prior i.i.d. per element, summed | one tensor of the declared `shape`; prior `.expand(shape)` and, for a single scalar log-density, `Independent(base, len(shape))` | numpyro: `fn.expand(shape).to_event(len(shape))`, or a plate — see §5.1 |
| **buffers** | read-only numpy array in `BufferSet`; passed through `context()` | `register_buffer(name, tensor, persistent=True)` | equinox: ordinary array field, **excluded from the trainable partition**; **never** `eqx.field(static=True)` — §7 |
| **`Parameter.unit`** | **lowers to nothing** — composition-time metadata only (§1.4) | **lowers to nothing** | **lowers to nothing** |

### 5.1 `to_event` versus `plate` — they are not interchangeable

Both give `size` draws; they differ in what the model claims about
independence, and numpyro uses that claim.

- **`to_event(n)`** declares the last `n` dimensions to be a single event: one
  log-density scalar for the whole array. This is the right lowering for an
  **array-valued parameter** (`shape=(4,)` per-channel offsets), which the
  contract defines as i.i.d. across elements and sums over. **Always pass `n`
  explicitly**: `to_event()` with no argument consumes *every* batch
  dimension, which for a parameter that has acquired a batch dimension from
  somewhere else silently reinterprets more than intended.
- **`numpyro.plate`** declares a batch dimension of conditionally independent
  draws, which enables subsampling and lets numpyro's machinery reason about
  the structure. This is the right lowering for a **`Plate`**, and it is what
  `parameters.md` §9 says the array-valued member shape was chosen for.

The two produce the same joint density here, because ampere's array-valued
priors are i.i.d. That makes the distinction easy to get wrong and hard to
notice. The rule is **structural, not numerical**: lower a `Plate` to a
`numpyro.plate` because it *is* a plate and downstream tooling (subsampling,
hierarchical diagnostics, ArviZ dimension naming) depends on the declaration;
lower a bare array-valued parameter with `to_event`. `Parameter.plate` is the
flag that tells them apart, and it is set precisely by `Plate.expand()`.

The same distinction on the torch side is less consequential — torch has no
plate primitive — so a plated parameter and an array-valued one lower
identically, and the plate structure survives only in the names and in the
provenance record. That is an accepted asymmetry: torch is not the PPL
backend.

### 5.2 Fixed parameters lower to buffers, and why that is not a category error

`architecture.md` §6 and `parameters.md` §10 are emphatic that a buffer is a
different *kind* of thing from a parameter — the distinguishing question is
"would you ever put a prior on it?" — and a fixed parameter answers "yes, in
principle; not in this fit".

They nonetheless lower to the same mechanism, because the mechanism answers a
narrower question: *does the optimiser take a gradient with respect to it, and
does it move with device and dtype?* A fixed parameter and a buffer give the
same answer to both. Keeping them distinct at the contract level and merging
them at the lowering level is the correct factoring: the user-visible
distinction (a fixed parameter is still in the model's declared vocabulary and
its provenance record; a buffer is not) is preserved where it matters, and the
runtime does the one thing that is actually needed.

The visible consequence is the one `parameters.md` §7 promised: fixing a
parameter changes which lowering row it takes, and changes nothing else — not
the model's `__call__`, not the buffer/parameter namespace, not the results
schema.

## 6. Module and pytree conventions

### 6.1 torch

A `Parameterised` object lowers to a `torch.nn.Module`:

- free parameters → `register_parameter(name, torch.nn.Parameter(t))`
- fixed parameters and buffers → `register_buffer(name, t, persistent=True)`
- names: the merged, qualified name (`"spectrum.temperature"`) is **not** a
  legal torch submodule path component, since torch splits `state_dict` keys
  on `.`. Lowering must either substitute the separator or nest real
  submodules per component. **Recommendation: nest.** One `nn.Module` per
  merge component, so `state_dict` keys come out as
  `spectrum.temperature` naturally and torch's own device/dtype recursion does
  the right thing. Tied parameters, being merged, belong to no single
  component — they go on the root module. This keeps the flat name and the
  module path identical, which matters because the flat name is also the
  ArviZ coordinate (W1.8).

`persistent=True` is the default and the right choice: a fixed parameter's
value is part of the fit's definition and belongs in a checkpoint.

### 6.2 jax: equinox partitions, and where paramax fits

`DEVELOPMENT_PLAN.md` §6 is explicit that **Paramax is a lowering mechanism,
not a user API** — which is exactly this document's remit, and why
`as_paramax()` was removed from the core API (`parameters.md` §11). Whatever
this section concludes, nothing a user writes mentions paramax.

A `Parameterised` object lowers to an `eqx.Module`:

- free parameters → ordinary array fields, **in** the trainable partition
- fixed parameters and buffers → ordinary array fields, **excluded from** the
  trainable partition
- **never** `eqx.field(static=True)` for any array — see §7

Equinox has **no buffer concept** — confirmed: there is no `eqx.Buffer`, and
`eqx.nn.State`/`StateIndex` are stateful-layer machinery (BatchNorm running
statistics), not PyTorch-style buffers. Three documented mechanisms exist for
"this array is not trainable", and they are not equivalent:

1. **An explicit `eqx.partition` filter spec.** Build a boolean pytree
   mirroring the module — the documented idiom is
   `jtu.tree_map(lambda _: False, model)` followed by `eqx.tree_at(...)` to
   flip the trainable leaves to `True` — then
   `trainable, frozen = eqx.partition(model, spec)` and
   `eqx.combine(trainable, frozen)` inside the differentiated function. The
   excluded leaves are literally absent (`None`) from the differentiated
   argument, so `eqx.filter_grad` returns `None` for them. (Equinox's own
   examples call the second half `static`, which is an unfortunate collision:
   it has nothing to do with `eqx.field(static=True)` and the arrays in it are
   still ordinary pytree leaves. Ampere's code should not reuse that name —
   §7 is the reason.)
2. **`jax.lax.stop_gradient`** applied in `__call__`. This is the Equinox
   FAQ's headline answer to "how do I mark arrays as non-trainable, like
   PyTorch's buffers?". The leaf *is* differentiated; its cotangent is zeroed.
3. **`paramax` wrappers.** `paramax.non_trainable(tree)` wraps each inexact
   array leaf in a `NonTrainable`, whose `unwrap()` applies `stop_gradient`.
   `paramax.NonTrainable(tree)` wraps a whole subtree instead; paramax's own
   docs prefer the function.

**Recommendation: (1), the `eqx.partition` filter spec**, with paramax
available as an interop layer if Phase 2's GP library choice demands it.

*(**Amended W2.15**, 2026-09-08: the recommendation was implemented — W2.5's
`ampere/backends/jax/parameters.py` builds an `eqx.partition` filter spec —
and **the interop caveat lapsed unspent**. The jax GP library was chosen by
measurement at W2.5 slice 2 and is `celerite2.jax`, not GPJax, so no
paramax-founded library sits at ampere's boundary and paramax is not a
dependency of the `jax` extra. The ruling therefore stands unconditionally
rather than pending a revisit.)*

This reverses the reading one might take from `DEVELOPMENT_PLAN.md` §4.1,
which names `paramax.NonTrainable` and an `eqx.partition` filter spec as
alternatives without ranking them. The plan's binding requirement is the
*negative* one — never equinox static fields — and both candidates satisfy it.
Ranking them needs one fact about how paramax actually works:

> `NonTrainable` is an `eqx.Module`, so the array it holds is **still an
> ordinary array leaf of the model**. It is *not* filtered out of
> `eqx.filter_grad`. The freezing happens only when `paramax.unwrap(model)` is
> called, which applies `lax.stop_gradient`.

That makes mechanisms (2) and (3) the same mechanism with different
ergonomics, and it means **forgetting to call `unwrap` inside the loss
silently trains the buffers**. Nothing raises; nothing warns; the model simply
learns its wavelength grid. That is precisely the failure signature these
contracts exist to eliminate, and it is the deciding argument: `eqx.partition`
*structurally cannot* leak, because the leaf is not in the differentiated
pytree at all, whereas paramax and `stop_gradient` both rest on a discipline
someone has to remember at every call site.

Two secondary considerations point the same way:

- **Regularisation and optimiser state still see a `NonTrainable` leaf**, because
  it is still there. paramax's own documentation warns about this and
  recommends filtering such leaves out anyway — i.e. falling back to (1) — to
  keep weight decay off them. A mechanism whose documentation tells you to
  add the other mechanism is not the simpler choice.
- **paramax is version 0.0.5, single-maintainer, pre-1.0.** It is alive and it
  is a reasonable dependency, but taking a hard dependency on a `0.0.x` API
  for something ampere can express with equinox core is an avoidable risk in
  the layer that every jax model passes through.

`DEVELOPMENT_PLAN.md` §6's "Paramax is a lowering mechanism, not a user API"
stands unchanged either way: it constrains where paramax may appear, not
whether it must. If Phase 2 selects a GP library founded on paramax wrappers
(GPJax is the candidate named there), ampere's modules must *interoperate*
with `unwrap` — accepting wrapped leaves from the library — and that is an
argument for supporting paramax at the boundary, not for adopting it as
ampere's own freezing mechanism. `paramax.Parameterize` is separately worth
noting as the ecosystem's idiom for applying a bijection at unwrap time; it is
not how ampere lowers §4's bijections, which are applied by the sampling layer.

## 7. Parameters versus buffers — the trap, restated as a rule

`DEVELOPMENT_PLAN.md` §7 lists it and `parameters.md` §10 points here, so this
document owns the operational statement:

> **Never mark an array `eqx.field(static=True)`.**

Equinox's own documentation says as much, in terms: *"`static=True` means that
this field is not a node of the PyTree, so it does not interact with any JAX
transforms, like JIT or grad. **This means that it is usually a bug to make JAX
arrays be static fields.** `static=True` should very rarely be used. It is
preferred to just filter out each field with `eqx.partition`."* (Note in
passing that `eqx.static_field` is deprecated in favour of
`eqx.field(static=True)`.)

The reason is mechanical, not stylistic. In equinox, a static field's *value*
goes into the flattening's auxiliary data — it is part of the pytree
**structure**, not its leaves. JAX's own documentation completes the argument:
auxiliary data "becomes part of the treedef, which JAX compares and hashes …
so it must support meaningful equality and hashing", and metadata fields "are
used to generate JIT cache keys". Equinox's FAQ states the consequence: a
function is recompiled "if any of its static (non-array) inputs change (as
measured by `__eq__`)". Putting an array there means:

- the array's contents are hashed/compared on every traced call;
- any change to the contents — a new wavelength grid, a re-read opacity
  table — is a *different* structure, so it triggers a full recompilation;
- large constant arrays are held alive in the compilation cache, so memory
  grows with the number of distinct arrays ever seen;
- and numpy arrays are not hashable, so depending on the path this is either a
  silent performance collapse or an error a long way from its cause.

Static is for genuinely static, small, hashable metadata: an integer size, a
string name, an enum, a boolean flag. `Plate.size` is static. A `Plate`'s
member values are not.

The correct home for a constant array is a **leaf that is excluded from the
trainable partition** (§6.2). It then behaves as buffers must
(`architecture.md` §5): traced through computations, moved by `device_put` and
dtype changes alongside parameters, never differentiated. torch's
`register_buffer` gives exactly these properties natively; equinox has no
buffer concept, and the partition filter is the replacement.

## 8. Plates and hierarchical priors, end to end

Worked shape, so Phase 2 has a target rather than a description. Given the
`parameters.md` §9 declaration:

```
Plate("objects", size=4,
      hyperparameters=[Parameter("mu", norm(0, 5)),
                       Parameter("sigma", halfnorm(0, 2))],
      members=[Parameter("theta", HierarchicalPrior("norm", {"loc": "mu", "scale": "sigma"}))])
```

`Plate.expand()` has already produced three ordinary parameters —
`objects.mu` (scalar), `objects.sigma` (scalar), `objects.theta` (shape
`(4,)`, `plate="objects"`, references `("objects.mu", "objects.sigma")`) — so
lowering never sees a `Plate` object. It sees parameters, one of which carries
a plate tag. The numpyro lowering is:

```
mu    = numpyro.sample("objects.mu",    dist.Normal(0.0, 5.0))
sigma = numpyro.sample("objects.sigma", dist.HalfNormal(2.0))
with numpyro.plate("objects", 4):
    theta = numpyro.sample("objects.theta", dist.Normal(mu, sigma))
```

Three things this pins:

- **Site names are the merged parameter names, verbatim.** They are already
  unique, already qualified, and already what `free_labels()` will use as
  ArviZ coordinates (W1.8). Inventing a second naming scheme at the numpyro
  boundary would put a translation layer between the sampler's output and the
  results schema, which is precisely the bug class §4.1 exists to end.
- **Emission order is `ParameterSet`'s topological order**, not declaration
  order. A hyperparameter is always sampled before anything referencing it;
  the core already computes that order (`_evaluation_order`) and lowering
  reuses it rather than recomputing.
- **The plate name is the `Plate`'s name**, and its size is static. Both are
  small hashable metadata, so on the equinox side they are legitimately
  `static=True` fields (§7) — unlike the member values.

A numpyro detail that ampere v1.3 is insulated from but a future version is
not: **`numpyro.plate` allocates its `dim` at construction, not on entry, and
the outermost plate takes `dim=-1` with inner plates going leftwards** — the
reverse of the ordering most readers assume. Worse, two plates *constructed*
before either is entered both take `dim=-1` and collide silently. Ampere v1.3
supports one plate dimension per parameter (`parameters.md` §12.2), so a
lowering emits one plate at a time and none of this bites. It becomes live the
moment nested plates are added, and at that point the lowering must pass `dim`
explicitly rather than relying on allocation order.

For the torch path there is no plate primitive: `objects.theta` is one tensor
of shape `(4,)`, and its hierarchical prior is constructed per evaluation as
`Normal(mu, sigma).expand((4,))` from the current `mu`/`sigma` tensors. The
hierarchical structure is real but implicit, surviving in names and provenance
rather than in a framework object.

**A `HierarchicalPrior` cannot be lowered once and cached.** Its arguments are
other parameters' *values*, which change every evaluation. On torch and in
numpyro's model function this is natural. It is worth stating because the
non-hierarchical rows *can* be built once at lowering time, and an
implementation that hoists distribution construction out of the evaluation
loop as an optimisation will break exactly the hierarchical rows, and break
them into a plausible-looking wrong answer (the prior frozen at its
initial-value hyperparameters).

## 9. RNG lowering and seed derivation

Three incompatible RNG models, and a policy that has to survive all three.

| Concept | reference (numpy) | torch | jax |
|---|---|---|---|
| the object | `np.random.Generator` | `torch.Generator` | a PRNG key (array-valued, threaded explicitly) |
| construction from an integer | `np.random.default_rng(seed)` | `g = torch.Generator(); g.manual_seed(seed)` | `jax.random.key(seed)` — the **typed** key, which is what to use; `PRNGKey` is legacy (still supported, not deprecated, but untyped, `uint32`, carries an extra trailing axis and no RNG-implementation information) |
| independent sub-stream | `rng.spawn(n)` / `SeedSequence(seed).spawn(n)` | a further `Generator` seeded from a derived integer | `jax.random.split(key, n)` |
| per-item stream | derive an integer per item | derive an integer per item | `jax.random.fold_in(key, i)` |
| statefulness | stateful, mutated in place | stateful, mutated in place, **and per-device** | **stateless**; keys are values and must be threaded |
| does the distribution API accept it? | **yes** — `frozen.rvs(random_state=rng)` | **no** — see §9.1 | n/a — every `jax.random` call takes a key explicitly |

### 9.1 torch's generator gap

`torch.distributions.Distribution.sample()` and `.rsample()` take
`sample_shape` and **nothing else**: there is no `generator=` parameter, and
the underlying calls (`torch.normal`, `torch.rand`, `_standard_normal`) are
made against the global RNG without one being threaded through. A torch-backed
draw is therefore *not* reproducible by handing a `Generator` to the
distribution, which is the obvious thing to try and which fails silently — the
draws are perfectly valid, just not the seeded ones.

Two working routes, and the choice matters:

1. **`torch.random.fork_rng()`** around a `torch.manual_seed(s)` plus the
   sampling call. Restores the global RNG state on exit, so it does not leak
   into a host process. Works for any distribution.
2. **Draw uniforms explicitly and push them through `icdf`**:
   `torch.rand(shape, generator=g)` does accept a generator. This gives a
   genuinely `Generator`-scoped stream with no global state touched at all.

**Recommendation: (2) where `icdf` exists (§3.6), (1) otherwise.** Route (2)
is the only one that is safe under threading and the only one that never
mutates process-global state — the same objection §10.1 raises against
`set_default_dtype`. Its limitation is exactly the `icdf` gap, so the two
routes are complementary rather than competing, and a backend needs both.

### 9.2 Policy

**One integer seed per run**, recorded in the provenance attrs alongside the
backend name and library versions.

**Named sub-streams, derived by a single shared function.** Different concerns
must not share a stream — prior sampling, sampler initialisation, an SBI
simulation budget and a posterior-predictive draw are separate uses, and
reusing one stream makes any of them silently dependent on how many draws the
others took. The derivation lives once in `ampere.core` as a pure function of
`(seed, label)` returning an integer:

```
substream(seed: int, label: str) -> int
```

implemented over a *stable* hash (`hashlib.blake2b` of the label, not Python's
`hash()`, which is salted per process and would make runs irreproducible
across invocations — a real and easily missed trap). Each backend then does
the idiomatic thing with that integer: `default_rng(substream(...))`,
`Generator().manual_seed(substream(...))`, `jax.random.key(substream(...))`.
Sharing the derivation, not the mechanism, is what makes "the prior-sampling
stream" mean the same thing on all three.

**jax keys are threaded, never stored.** A key on an `eqx.Module` field is a
trap: it makes the module's identity depend on RNG state, and a stale key
silently reuses draws. Keys are function arguments. `jax.random.fold_in`
requires a **scalar 32-bit** integer, which is the shape `substream` must
therefore produce; `split` requires a scalar key and rejects a batched one.

**numpyro's own seeding is order-dependent, and ampere must not fight it.**
The `seed` handler splits its key once per stochastic site, sequentially, in
trace order — it does not derive per-site keys by folding in the site name. So
the draws a numpyro model produces depend on the *order* the sites are
emitted in, which for ampere is the topological order of §8. That order is
deterministic and reproducible, so this is not a defect; it does mean that
adding a parameter to a model changes the draws of every parameter emitted
after it, and that a run's provenance must record the parameter-set spec (not
just the seed) for a result to be reproducible. `ParameterSet.to_spec()` is
already what W1.8 hashes into the InferenceData attrs, so the requirement is
met — but the reason it is *needed* is this, and it should not be discovered
later.

### 9.3 What reproducibility means, and what it does not

**Same seed + same backend + same library versions → identical results.** This
is the contract, and it is testable.

**Same seed across backends → NOT identical, and ampere must not claim
otherwise.** numpy's PCG64, torch's Mersenne/Philox and jax's Threefry are
different bit generators; even where the algorithm matches, the order in which
draws are consumed differs. Cross-backend agreement is asserted
**statistically** — a conformance row compares distributions, moments or
quantiles to tolerance, never sample-by-sample equality.

The temptation to make prior draws bit-identical across backends by generating
uniforms once in numpy and pushing them through each backend's `icdf` is worth
naming and rejecting: it works only for families implementing `icdf`, it does
not survive a sampler (which draws its own randomness internally), and it
buys a property nothing actually needs. The reference implementation's
`prior_transform` already gives a deterministic, backend-independent map from
the unit hypercube to parameter values for the case where determinism genuinely
matters — nested sampling — and that is the right tool for it.

## 10. dtype, device, and the jax x64 call site

`architecture.md` §5 fixes the policy — **float64 for all likelihood and GP
linear algebra, on every backend, always; CPU by default, GPU opt-in and
explicit; buffers move with parameters** — and delegates the mechanics here.

### 10.1 torch

Thread dtype explicitly through construction. Every tensor a lowering creates
takes an explicit `dtype=` (default `torch.float64`) and an explicit
`device=`. **Never call `torch.set_default_dtype`**: it is global mutable
process state, and a library that sets it changes the behaviour of every other
consumer of torch in the same interpreter.

Two consequences worth pinning.

**`.to(dtype)` does not touch integer buffers.** `torch.nn.Module.to` casts
only floating-point and complex parameters and buffers; integral ones are
moved between devices but keep their dtype. That is the behaviour ampere wants
— an index array or a channel-count buffer must not silently become a float —
but it means a lowering cannot rely on a blanket `.to()` to make everything
float64. Buffers holding physical quantities are created float64 at
registration.

**Distribution arguments infer dtype from the first tensor they are given.**
`torch.distributions` normalises its arguments through `broadcast_all`, which
promotes bare Python floats using `torch.get_default_dtype()` — float32 — but
switches to the dtype *and device* of the first `Tensor` argument it finds.
So `Normal(0.0, 1.0)` is float32, while
`Normal(torch.tensor(0.0, dtype=torch.float64), 1.0)` is float64 throughout,
with the bare `1.0` promoted to match. **This is the supported mechanism for
pinning a lowered prior's precision without touching global state**: pass at
least one argument as an explicitly-typed tensor on the intended device.
Passing tensors of *mixed* dtypes is not a supported input to that path and
should be avoided; construct all arguments at the target dtype.

### 10.2 jax, and the x64 activation problem

This is `architecture.md` §9's open item, assigned here.

**The facts the recommendation rests on.** With `jax_enable_x64` off — the
default — `jnp.zeros(5)` is float32. An *explicit* float64 request is not
quite silent: `jnp.arange(5, dtype='float64')` emits a `UserWarning`
("Explicitly requested dtype float64 … is not available, and will be truncated
to dtype float32") and returns float32. Current jax makes even that
configurable through `jax_explicit_x64_dtypes`, whose modes are `ALLOW` /
`WARN` (the default) / `ERROR`. jax's documentation is unambiguous that the
main flag is a whole-program setting:

> "The X64 flag is intended as a **global setting** that should have one value
> for your whole program, **set at the top of your main file**. A common
> feature request is for the flag to be contextually configurable … this turns
> out to be difficult to implement within JAX's programming model, where code
> execution may happen in a different context than code compilation."

Mechanically, `jax.config.update("jax_enable_x64", True)` sets a
**process-global** value, while the `jax.enable_x64(...)` context manager sets
a **thread-local** one; the flag participates in the JIT cache key. The
scoped form is newly non-experimental — it replaced
`jax.experimental.enable_x64`, which has been *removed*, not merely
deprecated — and the experimental version carried an explicit warning that it
was "fundamentally broken and can result in unexpected behavior … particularly
when used in conjunction with JAX transformations". That warning did not
survive the move out of `jax.experimental`; **its absence is a documentation
gap, not evidence the hazard is gone**, and ampere should not build a
numerical guarantee on the scoped form until that is established.

The recommendation:

**(a) Do not set the flag as an import side effect of `ampere.backends.jax`.**
This is the tempting option — `architecture.md` §5 currently reads as though it
were the plan ("must be set at the *earliest* possible point when
`ampere.backends.jax` is first imported") — and it is wrong, for the same
reason §5 itself forbids `torch.set_default_dtype`: a library that mutates
global interpreter state on import changes the numerical behaviour of every
other jax user in the process, silently, from an `import` statement they may
not have written themselves. jax's documentation addresses the flag to the
*application* author ("the top of your main file"), never to a library.

Note that §5's stated *policy* — "x64 on, always, for this backend" — is
unchanged by this. What changes is who turns it on: the application, loudly
prompted, rather than ampere, silently. **This is a recommendation that
narrows §5's implementation sketch, and if accepted it should be reflected
there.**

**(b) Ship an explicit, idempotent activation function.**
`ampere.backends.jax.configure_x64()` calls
`jax.config.update("jax_enable_x64", True)` — the process-global setter, not
the thread-local context manager — logs once at INFO, and is a no-op if the
flag is already set. It is documented as "call before any jax work", because
that is the only condition under which it is safe.

Deliberately **not** named `enable_x64`: jax itself now exports
`jax.enable_x64` as a *context manager* with different (thread-local, scoped)
semantics, and two same-named functions doing different things across a
library boundary is a support burden nobody needs.

**(c) Guard at construction, and raise rather than warn.** The jax backend's
construction path — where a likelihood, GP or model is built — checks the flag.
If it is off, it **raises**, with a message naming the three remedies:

1. set `JAX_ENABLE_X64=1` in the environment (jax reads it at its own import;
   this is the only route that is guaranteed to precede array creation, and
   therefore the one the documentation should lead with);
2. call `ampere.backends.jax.configure_x64()` before creating any arrays;
3. explicitly opt into reduced precision — `precision="float32"` on the run —
   which is `DEVELOPMENT_PLAN.md` §7's sanctioned "per-run opt-out for GPU
   throughput", and which must be recorded in the provenance attrs so a
   surprising posterior can be traced to it.

Raising rather than warning is the load-bearing choice. A warning in a
notebook scrolls away; the failure it precedes is a plausible-looking wrong
answer.

**(d) The multi-library host process, honestly.** If a host application has
already used jax at float32 and ampere is imported into it, there is no
correct automatic behaviour, and the spec should say so rather than pretend.
jax publishes **no guidance at all** for a library in this position, and no
supported per-scope precision isolation: the thread-local context manager is
documented by jax's own maintainers as the thing that "turns out to be
difficult to implement within JAX's programming model", and it carried a
"fundamentally broken" warning in its previous incarnation. Ampere should not
stake a numerical guarantee on it.

**What happens to arrays created before the flag flips is undocumented.** It
is not stated to be an error, and it is not stated to be a silent
re-truncation. The expected behaviour — existing arrays keep float32, new ones
are float64, and mixed-precision promotion follows — is an *inference* from
jax's promotion machinery, not a documented claim, and the design here is
built so that it does not have to be relied upon. If Phase 2 wants to depend
on it, it should be established empirically first.

The supported configurations are therefore exactly two:

- **whole process at x64** — the host sets `JAX_ENABLE_X64=1`, or calls
  ampere's activation before its own jax work; or
- **ampere explicitly at float32**, opted into per run, with the science
  caveat recorded in provenance.

The guard in (c) is what turns an unsupported configuration into an error
message naming these two options, instead of a silently degraded solve.

One partial third route is worth flagging for Phase 2 to evaluate rather than
adopt now: setting `jax_explicit_x64_dtypes` to `ALLOW` makes explicitly
requested float64 dtypes be honoured *without* flipping the global default. If
that reliably holds through a GP solve's intermediate operations — which is
exactly the question, since promotion rules still follow `x64=False` — it
would give a genuine per-library precision island. That is an empirical
question about a recently added flag, and this spec does not assume the answer.

**(e) Device.** CPU by default; `device=` is explicit and never auto-detected
(`architecture.md` §5). On jax this means an explicit `jax.device_put` to a
chosen device rather than relying on the default device, so that a machine
with a GPU present does not silently take a different code path from CI.

## 11. What the conformance suite owes (for W1.10)

Rows this spec creates. Each is a place a backend can be wrong without
crashing, which is the criterion for needing a mechanical check.

1. **Jacobian direction.** For each of `Identity`, `Log`, `Logit`, and at a
   point where `x ≠ y`, assert the native `log_abs_det_jacobian(y, x)` equals
   `ampere.log_abs_det_jacobian(y)`. Argument-order-sensitive by construction:
   the test must fail if the arguments are swapped. §2(b).
2. **`lnprior_unconstrained` agreement.** Native path versus
   `ParameterSet.lnprior_unconstrained`, to tolerance, at several points
   including near a boundary. §2(c).
3. **Distribution parametrisation, family by family.** For every row of §3.2,
   compare log-densities at a handful of points against the scipy original.
   `truncnorm`, `lognorm` and `uniform` are the rows that catch the
   conversions of §3.3; include a `truncnorm` case with `loc ≠ 0` and
   `scale ≠ 1`, since a case with defaults passes under the wrong mapping.
4. **Support preservation.** For every family, assert the lowered
   distribution's support matches the scipy one — the check that catches a
   dropped `loc` on `halfnorm`/`expon`/`gamma`.
5. **Unsupported families raise.** Asking torch for a `truncnorm` raises the
   shared lowering error, names the family and the backend, and does **not**
   silently fall back. §3.4.
6. **Plate structure.** A `Plate` of size N lowers to N draws with the right
   joint density, and the numpyro path produces a plate rather than a bare
   batched site. §5.1.
7. **Hierarchical priors track their hyperparameters.** Change `mu`, and the
   member's log-density changes accordingly — the check that catches a
   distribution object hoisted out of the evaluation loop. §8.
8. **Buffers are not trainable and not static.** No gradient flows to a
   buffer; and on the jax side, changing a buffer's *contents* does not change
   the pytree structure (which would trigger recompilation). §7.
9. **Fixed parameters take no sampler dimension** on any backend, and
   `free_size` agrees across all three.
10. **Seed reproducibility.** Same seed, same backend → identical draws.
    Across backends → agreement in distribution only, never element-wise. §9.3.
11. **x64 guard.** With the flag off, constructing a jax likelihood raises;
    the message names the three remedies. §10.2(c).

## 12. Open questions for review

1. ***Ruled 2026-09-03: raise, with the door left ajar.*** Option (a) is
   accepted — `default_bijection_for` raises for a discrete family — with
   Peter's modification: the route to *eventual* discrete support
   (numpyro-style enumeration, (variational) EM, SBI, nested sampling,
   bare optimisation such as Bayesian optimisation) must not be
   foreclosed, though none of it is in the current development plan.
   Three design constraints make that so, recorded for W1.13's landing
   (*landed at the freeze* — `CapabilityError` in `default_bijection_for`,
   `parameters.md` §6, §3.5 above amended):
   the refusal lives **only** in `default_bijection_for` — declaration,
   prior sampling, constrained-space `log_prob`, `prior_transform` (an
   inverse-CDF composition, and scipy's discrete families implement
   `ppf`) and §3.5's lowering-as-distribution are untouched, so the
   non-gradient routes keep working; the error is a typed *capability*
   refusal naming the reason, not a malformed-declaration error; and
   discreteness becomes queryable from the canonical prior description,
   so a future engine path branches on it rather than catching.
   *(Original question follows for the record.)*
   **Discrete parameters get a continuous default bijection.**
   `default_bijection_for(scipy.stats.poisson(3.0))` returns `Log(lower=0.0)`,
   inferred from the support with no discreteness check. `parameters.md` §12.8
   declares discrete parameters unexercised, so nothing is broken today, but
   the honest fix is upstream — `default_bijection_for` should raise for a
   discrete family rather than return a bijection that cannot be right. That
   is a §4.1 contract change and therefore not this document's to make; §3.5
   states the lowering-side rule in the meantime. **Recommend routing to W1.13.**
2. ***Ruled 2026-09-01: `eqx.partition` wins.*** Recorded in
   `DEVELOPMENT_PLAN.md` §2's decision table with §4.1 amended to match.
   The Phase 2 revisit noted below still applies if the GP library choice
   (GPJax) puts paramax-wrapped leaves at ampere's boundary.
   ***Revisit discharged W2.15 (2026-09-08):*** the choice was made at W2.5
   slice 2 by measuring celerite2.jax against tinygp, and it was
   `celerite2.jax`. GPJax was not adopted, paramax is in neither extra, and
   no wrapped leaf reaches ampere's boundary — so the ruling is now
   unconditional. *(Original
   question follows for the record.)* `DEVELOPMENT_PLAN.md` §4.1 named
   `paramax.NonTrainable` and an `eqx.partition` filter spec as
   alternatives without ranking them. §6.2 ranks them, and comes down on
   `eqx.partition`, because `paramax.NonTrainable` turns out not to remove
   the leaf from the differentiated pytree at all — it applies
   `stop_gradient` at `unwrap()` time, so forgetting to call `unwrap`
   silently trains the buffers. The plan's binding requirement (never
   equinox static fields) is unaffected. Worth
   revisiting at Phase 2 start when the GP library is chosen, since a
   paramax-founded library (GPJax) would put wrapped leaves at ampere's
   boundary regardless of what ampere itself uses internally.
3. ***Closed at the freeze (W1.13):*** the flat tie-label ruling stands —
   Peter ruled `parameters.md` §14.3 on 2026-09-01 (tie labels stay global
   and unqualified) and nothing since has reopened it — so §6.1's
   recommendation is settled as written: one `nn.Module` per merge
   component, tied parameters on the root module. *(Original question
   follows for the record.)*
   **Qualified names versus torch `state_dict` keys.** §6.1 recommends nesting
   one `nn.Module` per merge component so that dotted names come out
   naturally. This is clean for names produced by `merge`, but a tie label is
   deliberately unqualified (`"distance"`, not `"shared.distance"` —
   `parameters.md` §14.3), so tied parameters land on the root module. If that
   open question is ever resolved the other way, this recommendation should be
   revisited with it.
4. ***Ruled 2026-09-03: accepted.*** `LoweringError(AmpereError)` lands in
   `ampere.core.exceptions` at W1.13. Peter's phrasing ("the relevant
   exceptions", plural) supports the same disposition for the other
   homeless exception types on W1.13's list — `OptionalDependencyError`
   ratified in place, `ResultsError` moved to core per `results.md` §15
   R5's recommendation — to be confirmed in that pass. *(Original
   question follows for the record.)* **Where the shared lowering error type lives.** §3.4 proposes
   `ampere.core.exceptions.LoweringError(AmpereError)` — a sibling of
   `ContractError`, not a subclass, for the reason given there. Adding it puts
   a *backend-facing* exception in `ampere.core`, which is defensible (core
   already owns `OptionalDependencyError`, which is equally backend-facing,
   and the whole point is that all three backends raise the same type) but is
   a small widening of core's remit and should be ratified rather than
   assumed. `architecture.md` §9 already hands `OptionalDependencyError` to
   W1.13 to ratify or move; same disposition suggested.
5. ***Ruled 2026-09-01: guard-and-raise accepted.*** Recorded in
   `DEVELOPMENT_PLAN.md` §2's decision table; `architecture.md` §5 amended
   to match and its §9 item closed. One addition from Peter's review,
   verified against numpyro source: `numpyro.enable_x64()` is a thin
   wrapper over `jax.config.update("jax_enable_x64", ...)` with no
   separate numpyro-side dtype state, so the one jax flag governs numpyro
   draws too — the §10.2(c) guard covers the numpyro path, and Phase 2's
   numpyro lowering may call `numpyro.enable_x64()` as its idiom for the
   same switch. *(Original question:)* §10.2(a) narrows `architecture.md`
   §5's "set at the earliest possible point on first import" sketch: same
   policy, opposite mechanism — ampere never sets the flag on import, and
   instead guards and raises.
6. ***Ruled 2026-09-03: sanctioned, loudly.*** The reference fallback for
   a missing native `icdf` is accepted on two conditions: taking it emits
   a loud warning naming the families and the backend, and a `strict`
   switch turns that warning into a raise for the user who would rather
   fail than mix paths. §3.6 amended to match; the flag's home was fixed
   at the freeze — **the fitting problem's existing `strict` toggle**,
   threaded to the lowering path, one flag with one meaning (§3.6 has the
   Phase 2 implementation contract). *(Original
   question follows for the record.)* **Should the reference backend be a sanctioned fallback for a missing
   `icdf`?** §3.6 says yes and argues it is not the silent substitution §3.4
   forbids, because `prior_transform` has one mathematical definition. That
   reasoning is sound but it does put a numpy computation inside a nominally
   torch-backed run, which someone will eventually be surprised by. Worth an
   explicit ruling.
7. ***Ruled 2026-09-03: ratified in place.*** `substream` stays in
   `ampere.core` (`ampere/core/rng.py`) — pure stdlib+numpy, deliberately
   free-standing, and the alternative was three backends agreeing by
   convention. Closes `inference.md` §19 item 5 with it. *(Original
   question follows for the record.)*
   **`substream(seed, label)` belongs in `ampere.core`** (§9.2), which means
   the core gains a small RNG-policy responsibility it does not have today.
   That seems right — it is pure stdlib, and the alternative is three backends
   agreeing by convention — but it is a (small) addition to a frozen contract's
   surface and should be ratified rather than assumed.
8. ***Ruled 2026-09-03: accepted in principle, hardened.*** The
   registration hook sketched below is the mechanism — a module-level
   registry keyed on the *neutral* family name and the backend, plus the
   analogous per-backend slot for a custom `Bijection` — with three
   hardenings: registration refuses to overwrite a built-in row or an
   existing registration without an explicit `override=True`; every
   registered row is stamped user-registered in provenance (as this item
   already demands); and the conformance suite grows an *opt-in* battery
   a registrant can run against their own lowering (reference-vs-native
   agreement), so "bypasses the suite's guarantees" becomes "can
   self-certify". On the jax question that prompted the ruling: the
   registry is consulted at lowering time, before any tracing, so the
   mechanism itself cannot interact with `jax.jit`/`vmap`/gradients —
   what must hold is that a registered constructor returns trace-pure
   objects, the same requirement the built-in table rows meet, and the
   hook's documentation states that as a rule. The plumbing is Phase 2's,
   beside the backends that consume it.

   ***Amended W2.12:*** the `backend` this registry is keyed on is the
   backend's **one name everywhere** — the `BACKEND` capability flag its
   models and transformations declare (`inference.md`'s capability-flags
   section), the `backend` a composed `FittingProblem` reports, the
   `ampere_backend` in a run's provenance, and the `name` of that backend's
   conformance fixture. There is no translation table between them, and a
   backend wanting two spellings would be introducing one.

   ***Landed W2.6:*** `ampere.core.lowering` — `register_lowering`/
   `lookup_lowering` for prior families and `register_bijection_lowering`/
   `lookup_bijection_lowering` for the bijection slot, sharing one
   underlying keyed-by-`(kind, name, backend)` store; `run_registrant_battery`
   is the opt-in self-certification battery; `provenance_entries` renders
   the non-built-in rows a run consulted for
   `ampere.results.provenance.provenance_attrs`'s `extra=`. The module's own
   docstring carries the worked example (registration, the overwrite
   refusal and its `override=True` escape, a passing and a failing battery
   run, and the provenance entry), executed as a doctest by
   `tests/core/test_spec_doctests.py` — this document stays prose-only, per
   §0's "a document, not code". *(Original question follows for
   the record.)*
   **Backend-specific lowerings for user-defined cases** (raised by Peter,
   2026-09-01). A duck-typed prior and a custom `Bijection` currently
   evaluate on the reference path only; §3.4 and §4 both name "supply a
   backend-native equivalent alongside" as the extension point without
   specifying its plumbing. The natural shape is a registration hook keyed
   on the neutral name — `register_lowering(family, backend, constructor)`
   for prior families, and the analogous per-backend transform slot for a
   custom `Bijection` — mirroring how §4.4's likelihood families register.
   That would let a user with a genuinely custom prior run it natively
   (and serialise it, if their registration round-trips) instead of being
   confined to the reference path. Needs design care at exactly one point:
   a user-registered lowering bypasses the conformance suite's guarantees,
   so registered entries should be marked as such in provenance. Route to
   W1.13 for the decision in principle; the plumbing is Phase 2's, beside
   the backends that consume it.

## 13. Confidence and sourcing

In the style of `prior_art.md` §7. torch and jax were **not** installed (the
work item forbids it), so every native-API claim is documentation-sourced or
flagged.

| Area | Confidence | Basis |
|---|---|---|
| ampere-side facts: `PriorSpec` canonical keywords, inferred bijections, merge semantics, deferred/tie resolution | **High.** | Executed directly against the repository at `1d251ad` — `describe_prior` and `default_bijection_for` were run on all ten families in the tables, and the merge/collapse behaviour was read from `ampere/core/parameter.py` source, not inferred from the prose spec. |
| scipy parametrisations, including the `truncnorm` standardisation and `loguniform`/`lognorm` conventions | **High.** | Executed against scipy 1.17.1: supports and shape names printed directly. The `truncnorm(a=-1, b=2, loc=5, scale=3) → support (2, 11)` result in §3.3 is a real observation, not a recollection. |
| §2 direction convention, the two-argument Jacobian, `biject_to` vs `transform_to`, the `TransformedDistribution` sign | **High** for torch. | Verified against the PyTorch v2.13.0 source (`torch/distributions/{transforms,constraint_registry,transformed_distribution}.py`), which is what generates the rendered docstrings. `biject_to`'s docstring states the direction verbatim ("from `constraints.real` to the given `constraint`"); `log_abs_det_jacobian(x, y)`'s states "given input and output". The rendered `distributions.html` page exceeds the fetch limit and truncates before these entries, so the source tag was used rather than the HTML. |
| §3.2 torch rows; absence of `LogUniform` and `TruncatedNormal` | **High.** | Checked against the complete `__all__` of `torch/distributions/__init__.py` at v2.13.0 — an enumeration, not a search. `TruncatedNormal`'s absence from pyro was checked the same way against `pyro/distributions/__init__.py`. The `loguniform` construction is idiomatic-by-precedent (torch's own `LogNormal` *is* `TransformedDistribution(Normal, ExpTransform)`) rather than a documented recipe — **flagged**. |
| §3.6 `icdf` availability per family | **High.** | Read the per-family source: `Normal`, `Uniform`, `HalfNormal` and `TransformedDistribution` implement `icdf`; `Gamma` has `cdf` but not `icdf`; `Beta` and `Poisson` have neither; the base class raises `NotImplementedError`. |
| §4 bijection composition rows | **High** for torch. | `ComposeTransform` applies parts left-to-right in the forward direction (verified in source), `AffineTransform(loc, scale)` computes `loc + scale * x`, `ExpTransform` and `SigmoidTransform` have the stated domains/codomains. |
| §9.1 torch's missing `generator` argument | **High.** | Verified by reading `sample`/`rsample` in `distribution.py`, `normal.py`, `uniform.py` and `transformed_distribution.py`: no `generator` parameter and none threaded into the underlying `torch.normal` / `torch.rand` / `_standard_normal` calls. `torch.rand`'s own `generator=` parameter and `torch.random.fork_rng` were confirmed on their rendered doc pages. |
| §10.1 torch dtype mechanics | **High** for the two stated behaviours. | `Module.to`'s float-only casting is quoted from the rendered `torch.nn.Module` page; the `broadcast_all` dtype/device inference was read from `torch/distributions/utils.py`. **Unverified**: the behaviour when arguments of *different* tensor dtypes are passed — the text avoids relying on it. |
| §3.2 numpyro rows; `LogUniform`, `TruncatedNormal`, `HalfNormal` | **High.** | Verified against the NumPyro **0.21.0 release tag** (not `master`, which was checked and does differ in `continuous.py`/`transforms.py`/`primitives.py`). `LogUniform(low, high)` exists with both arguments required; `TruncatedNormal` is a factory function whose `low`/`high` are keyword-only; `TruncatedDistribution` accepts only a fixed set of real-supported base families. The shifted-half-normal-as-truncated-normal equivalence in §3.3 is my own derivation, checked against the densities — **it is mathematics, not a citation**. |
| §4 numpyro `biject_to` dispatch table | **High.** | Read from `transforms.py`'s registry at 0.21.0, entry by entry, including the `positive` vs `greater_than(0)` type distinction and the `domain=` arguments on the composed forms. |
| §6.2 equinox/paramax mechanisms, and the recommendation against paramax | **High** on the mechanism; the recommendation is judgement. | `paramax.NonTrainable.unwrap` was read in source: it is an `eqx.Module` performing `eqx.partition`/`lax.stop_gradient`/`eqx.combine`, so the leaf remains in the differentiated pytree and freezing depends on `unwrap` being called. paramax's own documentation supplies the regularisation caveat. Version 0.0.5 and its maintenance state come from PyPI/GitHub metadata. **The ranking of the three mechanisms is a design judgement**, and it deliberately narrows `DEVELOPMENT_PLAN.md` §4.1, which lists two of them without ranking. |
| §7 equinox static-field trap | **High**, with one attribution caveat. | The "usually a bug to make JAX arrays be static fields" sentence is quoted verbatim from Equinox's advanced-fields documentation. The *hashing / JIT-cache-key* half of the argument is **JAX's** wording, not Equinox's — jax's pytree documentation states that auxiliary data "becomes part of the treedef, which JAX compares and hashes", and `register_dataclass` states metadata "are used to generate JIT cache keys". Equinox's FAQ supplies the recompilation consequence and the concrete `unhashable type: ArrayImpl` error. The claim holds; it is assembled from two sources rather than quoted from one. |
| §9 RNG APIs | **High.** | `jax.random.key` vs `PRNGKey` (legacy, discouraged, **not** deprecated — no `DeprecationWarning`), `split`'s and `fold_in`'s scalar-only requirements, and numpyro's `seed` handler splitting sequentially per site rather than folding in site names, all read from source or from the rendered `jax.random` page. |
| §10.2 jax x64 semantics | **High** for what is documented; **explicitly flagged** for what is not. | Documented and verified: the flag is a global setting "set at the top of your main file"; contextual configuration is acknowledged by jax as difficult within its programming model; `jax.config.update` is process-global while `jax.enable_x64` is thread-local; the flag is in the JIT cache key; `jax.experimental.enable_x64` has been **removed** (not merely deprecated) in favour of `jax.enable_x64`, and the old one's "fundamentally broken" warning did not survive the move; the float64-truncation `UserWarning` text; the `jax_explicit_x64_dtypes` ALLOW/WARN/ERROR modes. **Not documented anywhere, and flagged as such in the text: what happens to arrays created before the flag flips.** Also not documented: any guidance for a library setting the flag on a host's behalf — jax says nothing, so §10.2(a) is reasoning from the docs' audience, not citing a rule. |
| §10.2(d) multi-library host recommendation | **Recommendation, not a verified fact.** | The *constraint* (process-global flag, must precede array creation) is documented; the *disposition* (two supported configurations, guard-and-raise) is a design judgement this document is making, and is the kind of thing review exists to overturn. |
| §3.2 second-tier rows (`lognorm`, `expon`, `gamma`, `beta`) | **Medium.** | `Gamma` and `Beta` constructor signatures were read in both libraries' sources; `Exponential`'s was **not** in either, so `rate = 1/scale` is asserted from the standard parametrisation rather than verified. The `loc`-must-be-zero constraints follow from the scipy supports, which were executed. Flagged in the table itself. |
| §4 torch `Identity` row | **Medium.** | torch's `biject_to(constraints.real)` is stated to return an identity transform; the registry entry itself was not read, so the table names the call rather than the object it returns. Behaviourally certain, nominally unverified. |
| §5.2, §6.1, §8 structural recommendations (fixed→buffer, module nesting, emission order) | **Design judgement.** | Derived from the contracts rather than from any external source; no external verification applies. Flagged as recommendations in the text. |

**Overall**: the rows most likely to cause a silent, science-affecting error —
the Jacobian argument order (§2), the `truncnorm` standardisation (§3.3), the
equinox static-field trap (§7) and paramax's unwrap dependency (§6.2) — are
the four that were checked hardest, three of them against library source
rather than prose. The weakest areas are the second-tier distribution rows,
which are flagged inline, and the undocumented question of what jax does with
arrays that predate an x64 flip, which §10.2(d) is deliberately designed not
to depend on. torch and jax were not installed, per the work item, so **no
claim here has been executed against the target libraries** — only against
scipy and ampere's own code. W1.10's conformance suite (§11) is where these
rows stop being documentation and start being tested; until it exists, treat
the tables as carefully sourced intent.
