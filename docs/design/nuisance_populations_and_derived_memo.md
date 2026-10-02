# Per-dataset nuisance populations, and the `Derived` parameter node — design memo (W6.11)

*Drafted 2026-10-02 by Fable (session ampere-51), on the plan's two Phase 6
design items (Peter, 2026-09-22, ruling on W5.12's point (3) and W5.8's
point (2)). Opus review: pending at the time of writing; the review's
findings are folded in below the line that names them. The product is a
memo, not an implementation (D8, ruled 2026-09-28: the memo in Phase 6, the
implementation items open Phase 7). §8 drafts the decision-log rows the
Phase 7 items carry, §9 names the conformance rows they owe, §10 drafts the
items themselves. Every claim about the code below was checked against
master `f7c05ba`; §1 reproduces the two probes it rests on.*

## 0. The two questions, and the short answers

**(1) Can a `Population` draw a dataset's own nuisance parameter — a GP
amplitude, a calibration scale — from a shared prior?** Yes, by one rule
change and no new routing. Today a population addresses its members by bare
local name, so `over` may name a model but not a dataset, whose merged names
are qualified (`likelihood.amplitude`). The change is that **an `over` entry
may be a qualified component path** (`"d0.likelihood"`), and the element
binding it produces carries the matching **qualified local path**
(`local_name="likelihood.amplitude"`). Nothing else moves: `distribute` stays
one level deep, the dataset's own retained mapping re-distributes the
element to its noise model exactly as it re-distributes every other value,
and the plate layout's "strip the member's own declaration" step learns to
strip a qualified leaf from a composite's *outer* declaration. §1's second
probe runs this through the live merge by hand and the amplitude arrives at
the noise model under its bare name. Routing stays in `ParameterMapping`,
which is what `hierarchical_population.md` §11 Q1 required; nothing is
sliced in `DatasetCollection`.

**(2) Can a parameter be a pure function of other parameters?** It needs a
fourth state. A `Derived` object sits in a `Parameter`'s prior slot the way
`HierarchicalPrior` does — *"what determines this parameter"* is exactly what
that slot means — holding a small, serialisable **expression** over other
parameters' names. A derived parameter occupies no sampler dimension,
contributes nothing to the prior, is computed after its inputs in the
topological order the set already keeps, may be referenced by a
`HierarchicalPrior`, may be a plate or population member, and lowers to a
`numpyro.deterministic` site on jax and a tensor computation on torch. The
two customers the plan names become declarations: Piironen & Vehtari's slab
as `shrinkage_horseshoe(tail="slab")`, and the non-centred population
`θ_i = μ + σ z_i` as a `Population` whose member `theta` is derived from an
internal standard-normal member `z`. The expression is a closed grammar
rather than a callable, for the same reason `to_spec()` refuses an opaque
prior: provenance, hashing and lowering all need to *read* it.

The two compose: a per-dataset GP amplitude drawn non-centred from a shared
log-normal is `over=["d0.likelihood", …]` with members `z ~ N(0, 1)` and
`amplitude = exp(mu + sigma * z)`. §4 writes it out. Land (2) first: it is
self-contained and the slab is waiting on it; (1) second; the combined
customer and its validation third.

## 1. What is true today, in the code's own words

### 1.1 The per-dataset case

A dataset joins a `FittingProblem`'s merge as a `ParameterMapping`
(`Dataset._components`: the instrument as *its* mapping, the likelihood as a
set), so its names reach the problem qualified one level —
`('likelihood.amplitude', 'likelihood.length_scale')` for a flexible
likelihood with no instrument parameters. `Population` addresses members by
bare local name and `_apply_populations` refuses a composite by name
(`ampere/core/parameter.py`, the `inner` check). The probe
(`probe_w611.py`, three datasets each with a `GaussianProcessNoise(Matern32)`
likelihood, a population over them with `amplitude ~ lognorm(s=spread)`):

```
ParameterError: population 'gp' routes a draw to component 'd0', which was
merged as a ParameterMapping (a Dataset, or another composite). A population
addresses its members by *bare* local name, and a composite's merged names
are qualified ('likelihood.scale'), so the draw would never reach a leaf.
Declare the population over the models those datasets name, or share one
quantity across them with a Tie.
```

The refusal is correct as written: the draw *would* never reach a leaf,
because three things stand in its way, and they are the whole of what item
(1) changes.

1. **`PlateBinding.local_name` must be a bare identifier.**
   `PlateBinding.__post_init__` runs `_check_local_name` on it:
   `plate binding local name name 'likelihood.amplitude' must not contain
   '.': qualified names are produced by ParameterSet.merge() and Plate
   expansion, never declared directly.` `Binding.local_name`, by contrast,
   already carries dotted paths in the composed `bindings` view (`('d1',
   'likelihood.amplitude', 1)` in the probe), so the type is ready and only
   the declaration-time check refuses.
2. **The plate layout strips a bare name.** `_apply_populations` removes
   the component's own declaration of the member so the element can arrive
   under that name without shadowing it; for a composite the declaration to
   strip is the qualified leaf (`likelihood.amplitude`) in the component's
   *outer* declaration list — the copy `_merge` already makes
   (`declarations = {component: list(sets[component])}`), never the
   frozen inner mapping. Left unstripped, `_merge`'s own guard fires:
   `plate binding would deliver 'gp.amplitude' to d0.likelihood.amplitude,
   but component 'd0' already receives a value under 'likelihood.amplitude'.
   Element routing must not shadow a component's own parameter`.
3. **`_apply_populations` refuses `inner` components outright** — the
   message above. With 1 and 2 in place the refusal becomes a requirement
   that a composite in `over` carry a path.

With 1 bypassed and 2 done by hand, the probe's second half shows the rest
of the machinery already composes:

```
merged names (stripped): ('d0.likelihood.length_scale', 'd1.likelihood.length_scale',
  'd2.likelihood.length_scale', 'gp.spread', 'gp.amplitude')   free 7
d1 receives: {'likelihood.length_scale': 1.0, 'likelihood.amplitude': 0.2}
d1's inner re-distribute -> likelihood: {'amplitude': 0.2, 'length_scale': 1.0}
bindings view for gp.amplitude: [('gp', 'amplitude', None),
  ('d0', 'likelihood.amplitude', 0), ('d1', 'likelihood.amplitude', 1), ('d2', 'likelihood.amplitude', 2)]
tied_names: ()
```

Seven free dimensions — three length scales, one spread, three amplitudes as
one array — one sample site for the amplitudes, each dataset's noise model
receiving its own element under `amplitude`, and the bindings view already
naming the leaf path. `distribute` was not touched; the dataset's retained
mapping did the second hop, which is the nesting rule (`inference.md` §4.5)
doing what it was ratified for. The same hop runs natively: both backends'
lowered problems call `self.dataset.instrument.mapping.distribute(chain)`
(`ampere/backends/{torch,jax}/problem.py`) on native values, and element
bindings take their element in the value's own array type since W5.12.

### 1.2 The `Derived` case

Three places record the gap, in these words:

- `parameters.md` §9 (`shrinkage_horseshoe`): the slab `λ̃_j² = c²λ_j²/(c²
  + τ²λ_j²)` "is a deterministic function of two sampled parameters, and
  this contract declares parameters and priors, not deterministic nodes …
  a `Derived` node is what would close the gap."
- `lowering.md` §3.2.1: "This contract has no deterministic node … so `θ_k
  = s · z_k` is formed by whatever consumes the knot variables, not
  declared. `WarpedKernel` does exactly that and it is the reason
  `non_centred=True` is its default."
- The W5.12 status row's carried note: "the centred parameterisation is the
  only one expressible (a non-centred `θ_i = μ + σ z_i` needs the `Derived`
  node already recorded as Phase 6)."

What the code offers a derived node to stand on: `Parameter` has three
states and `_validate_state`'s refusal names the slot's meaning — *"has
neither a prior nor fixed=True nor a shared_as tie label, so nothing
determines it"*; `ParameterSet.evaluation_order` is a public topological
order over `Parameter.references`, which both backends already walk
(`lowering.md` §8, W2.13's fold-in 9); `complete()` is the one place a
name-to-value mapping is filled in before anything evaluates; jax's
`numpyro_model` and `_resolved` and torch's `_resolved` each fill a
`resolved` dict in that order, which is where a derived value is formed and
where a `HierarchicalPrior.bind` would find it; `ArrayOps` (`kernels.py`,
W4.5) is the backend-neutral namespace a kernel's closed form is written
against once — it has `exp`, `log1p`, `absolute`, `cos`, `sin` and not
`sqrt` or `log`. The blast radius of a new prior-slot type is countable:
`isinstance(…, HierarchicalPrior)` or `is_hierarchical` appears at 17 sites
in `core/parameter.py`, 3 in `backends/torch/parameters.py`, 2 in
`core/kernels.py`, 1 in `backends/jax/parameters.py`, 1 in
`results/provenance.py` — 24 branch points, each of which gains one case.

### 1.3 Two findings on the way (recorded, not fixed — ground rule 4)

- **`Population` declarations are not in provenance.** `grep Population
  ampere/results/provenance.py` is empty. The *merged* spec differs between a
  problem with and without a population (`gp.amplitude` versus
  `d0.likelihood.amplitude`), so the spec hash does move, but the declaration
  itself — name, members, hyperpriors, `over`, layout — is recorded nowhere
  a reader can see it. The per-component hashes (`provenance.py` ~649,
  `hash_of(merged.to_spec())` per dataset) hash each dataset's *own*
  mapping, which after item (1) no longer declares the leaf the population
  replaced. Item W7.1 carries the fix: `populations` in the problem record.
- **A routed member a component does not declare fails late.** In the plate
  layout `_apply_populations` adds the element binding whether or not the
  component declares the member; a model receiving an undeclared key fails at
  its first evaluation, in `Parameterised.context` → `ParameterSet.complete`:
  `got value(s) for unknown parameter(s) [...]`. The contract's "consuming it
  is the composing caller's contract" is true but the failure should be at
  merge, by name. Item (2)'s internal-member rule (§3.5) is the occasion to
  make it so.

## 2. Item (1): a population over a dataset's qualified path

### 2.1 The surface

`Population.over` accepts entries of the form `component[.path]`:

```python
Population(
    "gp",
    members=[Parameter("amplitude", HierarchicalPrior("lognorm", {"s": "spread"}))],
    hyperpriors=[Parameter("spread", st.halfnorm(0.0, 1.0))],
    over=["d0.likelihood", "d1.likelihood", "d2.likelihood"],
)
```

and `DatasetCollection.plate` gains `within=`, the convenience that writes
those entries from the datasets' own labels:

```python
DatasetCollection.plate("gp", datasets, members=[...], hyperpriors=[...], within="likelihood")
```

A bare entry keeps today's meaning and today's rule: a plain-set component
(a model) addressed by bare local name. An entry with a path names a
composite and the component *inside* it that declares the member. The two
motivating cases:

| Case | `over` entry | member | leaf stripped / re-priored |
|---|---|---|---|
| each dataset's GP amplitude from one shared prior | `"d0.likelihood"` | `amplitude ~ HierarchicalPrior("lognorm", {"s": "spread"})` | `d0.likelihood.amplitude` |
| each dataset's calibration scale under a fitted spread | `"d0.instrument.calibrate"` | `scale ~ HierarchicalPrior("lognorm", {"s": "spread"})` | `d0.instrument.calibrate.scale` (two levels: the dataset's mapping, then the instrument's) |

### 2.2 The rules, in `_apply_populations`'s terms

For an `over` entry `component.rest`:

1. `component` is a top-level component of this merge (as today).
2. `rest` resolves **through the retained inner mappings**, one segment per
   level: `inner[component].components` contains the first segment; if that
   component was itself merged as a mapping, the next segment resolves in
   *its* `inner`, and so on. A segment that does not resolve is refused by
   name with the components available at that level. A path given for a
   component merged as a plain set is refused ("`d0` is a `ParameterSet`,
   not a composite; address it without a path").
3. The member's leaf is `f"{rest}.{member.name}"`, looked up in the
   component's **outer declaration** (`declarations[component]`, the
   dataset's merged names). The existing checks run against it unchanged —
   shape, unit, fixed, tied or `shared_as` — plus one new one: **the leaf
   must exist** (a composite cannot be *given* a leaf it never declared,
   because its inner routing table has no row for it; a model may, as
   today, since the element arrives as an extra key its own set declares).
4. **Plate layout**: strip the leaf from the outer declaration; add
   `PlateBinding(parameter=population.qualified(member), component=component,
   local_name=f"{rest}.{member.name}", index=i)`. **Flat layout**: replace
   the leaf in the outer declaration with the member, re-priored, its
   hierarchical references renamed onto the qualified hyperpriors — exactly
   the bare case with a longer name.
5. `PlateBinding.local_name` and `Binding.local_name` admit a qualified
   path; `_merge`'s shadow check compares the full path, so a path cannot
   shadow a leaf the composite still declares.

Nothing in `distribute` changes. The outer `distribute` hands `d0` a key
`likelihood.amplitude`; the dataset's own `mapping.distribute` (the second
hop both the reference path and both realisations already take) routes it
to `likelihood` under `amplitude`. The inner mapping's merged set still
*declares* `likelihood.amplitude` with the dataset's original prior; that
prior is never evaluated by a problem (only the outer merged set's
`lnprior` is), which is why the strip can be on the outer copy alone. It
is also why §1.3's provenance finding matters: the dataset's own spec is no
longer the truth about that leaf.

### 2.3 Why routing stays in one place

`hierarchical_population.md` §11 Q1 weighed `Binding.index` against "let
`DatasetCollection` slice before it calls the model" and chose the binding
because slicing "puts parameter routing in two places … and every other
consumer that walks the bindings — provenance, ArviZ labelling, W1.9's
lowering — would see an incomplete picture." The plan's bullet restates it
for this item: the lift is "`PlateBinding`/`Binding` carrying a qualified
path … not slicing in the dataset collection." §2.2 honours it exactly: one
new binding with a longer `local_name`, the routing table the single source,
the bindings view already truthful (§1.1's probe: `sites_of("gp.amplitude")`
names the dataset-level leaf path for every element), `_index_coordinate`
in `results/emission.py` reading the element bindings' components to label
the plate dimension with the dataset labels — which is the labelling
`hierarchical_population.md` §10.2 demanded and which this item gets for
free, since the components *are* the datasets.

One refinement to the bindings view: `ParameterMapping.bindings` does not
descend an element binding ("the element is the leaf"). With a qualified
path the leaf *is* the path, so the view is already right without
descending; `global_name_for("d0", "likelihood.amplitude")` resolves because
the composed view is searched too. No change.

### 2.4 Sharing and population remain alternatives, per quantity

A population over `d*.likelihood` for `amplitude` and a `Tie` over
`d*.likelihood.length_scale` coexist in one problem — the existing refusal
of a tied or `shared_as` *site* as a population member is per leaf, not per
component, and stays. `inference.md` §9's sentence "sharing and population
are alternatives, not layers" is still true at the level it was said.

### 2.5 What this does not do

- Nested plates (`parameters.md` §12.2): still one plate dimension per
  parameter. A population whose `over` entries sit at different depths
  (`"d0.likelihood"` and `"m1"`) is *allowed* — each entry resolves on its
  own — provided the leaves agree in shape and unit, as today.
- Members with dotted names: rejected. `Parameter.name` stays an
  identifier; the path lives on the `over` entry, where it describes *where*
  the member goes rather than *what* it is.
- A per-element path (different members to different leaves within one
  component): not needed by either case; one path per entry.

### 2.6 Contract text this amends

- `parameters.md` §8 "Plate bindings": `PlateBinding.local_name` may be a
  qualified path relative to the component, consumed by the component's
  retained mapping; `Binding.local_name` likewise. §9 "What it refuses": the
  composite refusal rewritten as "a composite in `over` without a path", and
  the `over` grammar stated. §12 gains no new limitation.
- `inference.md` §9 "A population of datasets": the paragraph "the draws are
  routed to the models the datasets name, not to the datasets themselves,
  and that is a rule rather than a default" rewritten — it was a rule for
  want of a path; `DatasetCollection.plate(within=)` documented beside
  `over=`.
- `hierarchical_population.md` §11 Q1: a dated note that the second
  addressing form landed by the same mechanism, routing still in one place.
- `results.md` §9: the problem record gains `populations` (§1.3's finding).

## 3. Item (2): the `Derived` node

### 3.1 The declaration

```python
Parameter("theta", Derived("mu + sigma * z"))
Parameter("s", Derived("sqrt(c**2 * tau**2 * lam**2 / (c**2 + tau**2 * lam**2))"), shape=(3,))
```

`Derived(expression: str)` is a frozen dataclass in `ampere.core.parameter`
beside `HierarchicalPrior`, and `AnyPrior` admits it. A `Parameter` whose
prior slot holds one is in the fourth state, **derived**: `is_derived` true,
`is_free`, `is_fixed` and `is_deferred` false. `_validate_state` accepts it
as "something determines it"; `fixed=True`, `shared_as`, `value` and
`bijection` alongside a `Derived` are refused by name (a derived parameter
has no draw to fix, tie, initialise or transform — fix, tie or transform its
inputs).

### 3.2 The expression grammar, and why not a callable

The expression is parsed once with `ast.parse(mode="eval")` and the tree is
checked against a whitelist: `Name`, numeric `Constant`, `BinOp` over `+ - *
/ **`, `UnaryOp` with `-`, parentheses (free), and `Call` to exactly `sqrt`,
`exp`, `log`, `log1p`, `abs`. Everything else — attribute access,
subscripts, comparisons, conditionals, any other call — is refused at
declaration with the offending node named. The names are the `inputs`
(`Derived.references`, in first-appearance order); `rename_references`
rewrites them exactly as `HierarchicalPrior.rename_references` does, which
is what lets qualification, tie collapse and a population's hyperprior
rename all pass through unchanged. `to_dict` records the normalised source
(`ast.unparse` of the checked tree) so the spec round-trips and hashes
stably. Evaluation is `Derived.evaluate(ops, resolved)`: one walk of the
tree over an `ArrayOps`-shaped namespace, so the same expression runs in
numpy (`complete()`), torch and jax (each backend's `_resolved`) with no
per-backend transcription — W4.5's "one builder or three" answer, reused.
`ArrayOps` gains `sqrt` and `log` (both backends' namespaces are one line
each) and a `derived(expression, resolved)` entry point is **not** added:
the evaluator is core's, the namespace the backend's.

A callable (`Derived(lambda mu, sigma, z: mu + sigma * z)`) was the
obvious alternative and is rejected for this version, for the reason the
contracts already give about opaque priors: `to_spec()` "raises rather than
producing a provenance record that silently omits it" (`parameters.md`
§7). A callable cannot be hashed into `ampere_problem_hash`, cannot be
serialised into a stored run, cannot be lowered without trusting it to be
traceable on each backend, and cannot be shown to a reader of the results.
The grammar above covers both customers and every non-centred form in the
contracts; a conditional (`where`) is the first thing it will want, and it
is one whitelist entry when a customer appears. Recorded as a limitation
(§3.7); the escape hatch is what exists today — the consumer forms the
quantity itself, as `WarpedKernel` does.

Shape: declared as for any parameter (`shape=(3,)` above). The expression's
result must broadcast to it; checked at `ParameterSet` construction when
every input has a declared `value` and otherwise at the first `complete()`,
refused by name either way. Unit: declared, carried into provenance, **not
checked through the arithmetic** — priors are declared numerically in the
parameter's unit and so is a derived expression; recorded as a limitation.

### 3.3 What each consumer does with a derived parameter

| Consumer | Behaviour |
|---|---|
| `free_size`, `free_names`, `free_labels`, `free_slice`, `pack`, `unpack`, `bijections`, `prior_transform`, `constrain`, `unconstrain` | **excluded** — not a sampler dimension, like a fixed parameter |
| `names`, `__iter__`, `__getitem__`, `evaluation_order` | **included** — it is a parameter of the set, ordered after its inputs (its `references` are the expression's names, so `_evaluation_order` needs no new logic; cycles are refused by the same walk) |
| `complete(values)`, `unpack(theta)` | **computed** in evaluation order from the resolved inputs; a derived name *supplied* in the mapping is refused by name ("is computed, not given") |
| `lnprior`, `lnprior_unconstrained` | **skipped** — no density; its inputs carry theirs |
| `sample(rng)` | computed from the sampled inputs |
| `HierarchicalPrior.bind(resolved)` | finds it in `resolved` like any other name — a hierarchical prior may reference a derived parameter, which is the slab's whole point |
| `Tie`, `shared_as`, `PlateBinding` on it | tie/share **refused** ("tie the inputs"); element routing **allowed** (a plate-tagged derived array is the non-centred population's `theta`) |
| `Plate.expand`, `Population` member | allowed; the plate tag and shape `(size,) + member.shape` applied as to any member; references renamed onto the qualified hyperpriors |
| `to_spec` / `from_spec` | `{"derived": {"expression": "..."}}`; opaque never arises |
| `fix`, `release` | refused |
| default bijection, `default_bijection_for` | never asked (not free) |

### 3.4 Lowering

`lowering.md` §5's declaration-form table gains a row: **derived → no site**.
On jax, `numpyro_model` computes the value in `_order` and emits
`numpyro.deterministic(name, value)` so the sample stream carries it under
the merged name (numpyro's own convention for exactly this), and
`_resolved`/`log_prior` compute it into `resolved` before any site that
references it. On torch, `_resolved` does the same; there is no plate
primitive and no deterministic primitive, so the structure survives in
names and provenance as the plate does (`lowering.md` §8). The
`lowering_provenance` rows gain a `kind: "derived"` entry per derived
parameter so a run records which names were computed rather than sampled.
No bijection is lowered. `_lower_site` returns `None` for a derived
parameter as it does for a fixed one, and both backends' `_resolved` gain
the computation loop — the one new piece of backend code, written once per
backend against `ArrayOps`.

### 3.5 Population members, and the internal-member rule

The non-centred population:

```python
Population(
    "objects",
    members=[
        Parameter("z", st.norm(0.0, 1.0)),                    # internal
        Parameter("theta", Derived("mu + sigma * z")),        # routed
    ],
    hyperpriors=[Parameter("mu", st.norm(0.0, 5.0)), Parameter("sigma", st.halfnorm(0.0, 2.0))],
    over=["obj0", "obj1", "obj2"],
)
```

`z` must not be routed: no member component declares it, and today's
behaviour (§1.3) would be a late `unknown parameter` failure at the first
model evaluation. The rule: **a member that no `over` component declares
is routed to none of them; it is permitted only if it is an input of a
derived member of the same population, and refused by name at merge
otherwise.** The first half makes `z` population-internal (it lives under
`objects.z`, one array site, plate-tagged, sampled alongside `mu` and
`sigma`, and reaches the members only through `theta`); the second half
turns the late failure into a merge-time one for the typo'd member, which
closes §1.3's second finding in passing. For the flat layout the same rule
reads: an internal member is declared once on the population's own
component, not per member.

The declared density is the centred one — `θ_i ~ N(μ, σ)` — written so that
the sampled coordinates are `(μ, σ, z_i)`: `lowering.md` §3.2.1's right-hand
column, now declarable rather than formed by a consumer. `WarpedKernel`'s
`non_centred=True` *could* be re-expressed this way and is left alone: it
works, and rewriting a landed kernel for symmetry is not a customer.

### 3.6 Results

The `posterior` group gains **one variable per derived parameter**, under
its merged name, computed from the stored draws at emission
(`_posterior_variables` walks `problem.parameters`; a derived one is
`complete()`d from the draw rather than sliced from it) — on **every**
engine, not only NUTS, so an emcee run of a slab-shrunk `Sum` shows the
effective scales beside the amplitudes. A root attr `ampere_derived` lists
the names, so a diagnostic that must not treat a deterministic function of
draws as a sampled dimension can tell (R-hat and ESS on it are meaningful;
`plot_trace` may show it; the SBC and TARP routes operate on the free vector
and never see it). This is a results-schema change: `PROVENANCE_SCHEMA_VERSION`
→ 10. `Optimum` stays on the free vector (`inference.md` §10b: the
unconstrained and constrained free vectors along `free_parameter`); a
derived value at the optimum is `complete()` of the constrained vector and
is left to the caller.

Everything else is unaffected by construction: the SBI encoding packs the
free vector (`encoding.md`), `simulate` draws free parameters and
`complete()` fills the rest, `fit_population`'s reweighting reads
`log_prior`/`log_likelihood`, and `prior_transform` maps the unit cube onto
free dimensions only.

### 3.7 Deliberate limitations

1. **No callable.** §3.2. The consumer forms it, as today.
2. **No conditional.** `where` is one whitelist entry away; no customer yet.
3. **No unit arithmetic.** The declared unit is recorded, not derived.
4. **Not a sampler dimension, ever.** A user who wants to *sample* `θ` and
   *derive* `z` has written the centred form; declare it that way.
5. **No derived input to a plate size or a kernel's structural argument.**
   A derived value is a value, not configuration.

### 3.8 Contract text this amends

- `parameters.md` §3: three states become four, with the slot's meaning
  stated; a new §9 subsection "`Derived` — a parameter that is a function of
  others" with the grammar, the table of §3.3 and both customers; the
  `shrinkage_horseshoe` subsection's "a `Derived` node is what would close
  the gap" discharged by `tail="slab"`; §12 gains limitations 1–3 of §3.7.
- `lowering.md` §3.2.1: the "where the multiplication lives" paragraph
  rewritten — it lives in a `Derived`; §5's table row; §8's worked shape
  extended with the non-centred population and its `deterministic` site.
- `results.md` §4 (the posterior's derived variables and `ampere_derived`),
  §9 (schema 10).
- `inference.md` §9 (the non-centred population as the recommended form for
  NUTS), §10a coverage ("plates, hierarchical priors and derived
  parameters").

## 4. The two together: a non-centred per-dataset GP amplitude

```python
collection = DatasetCollection.plate(
    "gp", datasets, within="likelihood",
    members=[
        Parameter("z", st.norm(0.0, 1.0)),
        Parameter("amplitude", Derived("exp(mu + sigma * z)")),
    ],
    hyperpriors=[Parameter("mu", st.norm(-3.0, 2.0)), Parameter("sigma", st.halfnorm(0.0, 1.0))],
)
```

Each dataset's `likelihood.amplitude` is stripped from its outer declaration
and replaced by element *i* of a plate-tagged derived array `gp.amplitude`;
the sampler sees `(gp.mu, gp.sigma, gp.z[0..N-1])` plus every length scale;
NUTS on torch or jax walks a funnel-free geometry; the `posterior` group
carries `gp.amplitude` as a derived variable with the dataset labels as its
plate coordinate. This is the flexible likelihood's own hierarchical prior
— the case the plan calls "the motivating case" — and it is the validation
customer for the third Phase 7 item: the M2 pattern over many spectra, the
per-dataset amplitudes' shrinkage towards the population versus independent
amplitudes, pinned as a margin (W4.5's `_period_margin` precedent).

## 5. Alternatives considered and rejected

- **Slicing in `DatasetCollection`** (§11 Q1's alternative): routing in two
  places; rejected at the freeze and again here.
- **Dotted member names**: `Parameter.name` an identifier is load-bearing
  (values reach models as keyword arguments); rejected.
- **A `DerivedParameter` class beside `Parameter`**: every `for parameter in
  set` consumer branches on type instead of on state; the prior-slot object
  keeps the 24 branch points of §1.2 as the complete list of edits and
  leaves `ParameterSet`'s iteration protocol alone. Rejected.
- **A callable `Derived`**: §3.2.
- **Derived values in their own results group** (`derived` beside
  `posterior`): numpyro and ArviZ both put deterministic sites in
  `posterior`; a separate group would make every plot and summary a special
  case for a quantity that is, per draw, a posterior quantity. The attr
  marks them instead. Rejected.
- **Sampling the derived parameter and deriving the input** ("let the user
  pick which is free"): that is a reparameterisation, which is what
  `bijection=` is for when it is one-to-one and what the centred declaration
  is otherwise. Rejected.

## 6. Cost and risk

- (1) touches `PlateBinding.__post_init__`, `_apply_populations`,
  `Population.__post_init__` (the `over` grammar), `DatasetCollection.plate`
  (`within=`), the provenance problem record, and the contract pages of
  §2.6. Nothing in `distribute`, nothing in either backend. Risk: the path
  resolution through nested `inner` mappings — the one new algorithm — and
  the stripped leaf's absence from the dataset's own spec (§1.3).
- (2) touches the 24 branch points of §1.2, adds `Derived`, two `ArrayOps`
  methods on three namespaces, a computation loop in each backend's
  `_resolved`, `numpyro.deterministic` in `numpyro_model`,
  `_posterior_variables`, the `ampere_derived` attr and the schema bump,
  and the contract pages of §3.8. Risk: the evaluation-order dependence
  (`lowering.md` §8's warning that hoisting distribution construction out of
  the loop breaks hierarchical rows applies verbatim to derived values), and
  broadcasting of a plate-tagged derived array against scalar hyperpriors —
  both pinned by conformance rows below.
- Both are §4 changes (ground rule 9): a decision-log row each, in the same
  PR as the conformance rows, the suite green or amended with justification.

## 7. Open questions for Peter

1. **The grammar versus a callable** (§3.2). Recommendation: the grammar,
   with `where` added when a customer appears.
2. **Derived variables in `posterior` with an attr, or a separate group**
   (§3.6). Recommendation: `posterior`.
3. **`over` paths as the low-level form plus `within=` on the factory**
   (§2.1), or `within=` on `Population` itself. Recommendation: paths on
   `over` (one grammar, heterogeneous depths allowed) and the convenience
   on the factory only.
4. **Order of landing** (§10): Derived, then paths, then the customer.
   Recommendation: as listed; (2) is self-contained and (1) is a smaller
   change once the internal-member rule exists.

## 8. The decision-log rows (drafted for the Phase 7 items to carry)

| Topic | Decision |
|---|---|
| Per-dataset nuisance populations: a `Population` over a qualified component path (W7.1) | **`Population.over` entries may be qualified component paths (`"d0.likelihood"`, `"d0.instrument.calibrate"`); `PlateBinding.local_name` and `Binding.local_name` may be qualified paths relative to the component; the plate layout strips the addressed leaf from the composite's outer declaration; `DatasetCollection.plate(within=)` writes the entries from the dataset labels; routing is unchanged.** The population sketch's §11 Q1 ruling (2026-09-02) that routing lives in `ParameterMapping` alone is kept: the one new binding has a longer local name and the dataset's retained mapping takes the second hop it already takes for every other value, so neither `distribute` nor either backend changes. The W5.12 refusal of a composite in `over` becomes a refusal of a composite *without a path*. A leaf the composite does not declare is refused (a composite's inner routing has no row for it); the bare-name case is unchanged. The problem's provenance record gains `populations` (name, members, hyperpriors, `over`, layout): the merged spec already moved with a population but the declaration was recorded nowhere, and after this change a dataset's own spec no longer declares the leaf the population replaced. Motivated by the per-dataset GP amplitude under a shared prior — the flexible likelihood's natural hierarchical prior across many spectra — and the per-dataset calibration scale under a fitted spread, neither declarable before. `parameters.md` §8/§9, `inference.md` §9, `hierarchical_population.md` §11, `results.md` §9 amended; the conformance rows of the memo's §9.1. |
| The `Derived` parameter node (W7.0) | **A fourth `Parameter` state, `derived`: `Parameter(name, Derived("<expression>"))`, the expression a closed grammar (names, numeric literals, `+ - * / **`, unary minus, `sqrt exp log log1p abs`) parsed once, serialised as its normalised source, evaluated over `ArrayOps` on every backend; no sampler dimension, no prior term, computed after its inputs in `evaluation_order`; referenceable by a `HierarchicalPrior`; a legal plate and population member; lowered as `numpyro.deterministic` on jax and a tensor computation on torch; one `posterior` variable per derived parameter on every engine, named in `ampere_derived`; `PROVENANCE_SCHEMA_VERSION` → 10.** A callable is refused for the reason `to_spec()` refuses an opaque prior: provenance, hashing, serialisation and lowering must read it. `ArrayOps` gains `sqrt` and `log`. A population member that no `over` component declares is routed nowhere and is permitted only as an input of a derived member of the same population — the non-centred `θ_i = μ + σ z_i` — and refused by name otherwise, which also moves the undeclared-member failure from the first model evaluation to the merge. Discharges the two recorded gaps: `shrinkage_horseshoe(tail="slab", slab_scale=c)` declares Piironen & Vehtari's slab exactly, and `lowering.md` §3.2.1's "where the multiplication lives" has an answer that is a declaration. Deliberate limitations: no callable, no conditional, no unit arithmetic. `parameters.md` §3/§9/§12, `lowering.md` §3.2.1/§5/§8, `results.md` §4/§9, `inference.md` §9/§10a amended; the conformance rows of the memo's §9.2. |

## 9. The conformance rows (named; `tests/conformance`, once per registered fixture unless stated)

### 9.1 `test_population.py`, a new class `TestAPopulationOverADatasetPath`

1. `test_the_leaf_is_stripped_and_the_element_reaches_the_noise_model` —
   three flexible-likelihood datasets, `over=["d*.likelihood"]`,
   `amplitude ~ lognorm(s=spread)`: merged names contain `gp.amplitude` of
   shape `(3,)` and no `d*.likelihood.amplitude`; `evaluate` at a θ with
   `gp.amplitude = (a0, a1, a2)` equals, per dataset, the log-likelihood of
   the same dataset built alone with its amplitude fixed at `a_i`
   (`tolerances.cross_solver`).
2. `test_the_joint_log_prior_decomposes_over_the_path` — `lnprior` equals
   the hyperprior's density plus the sum of the members' densities at the
   sampled spread.
3. `test_a_two_level_path_reaches_an_instrument_step` — `over=["d*.
   instrument.calibrate"]`, member `scale`: the calibration step receives
   its element (the transformed spectrum scales by it).
4. `test_the_two_layouts_agree_on_the_joint_density` — plate versus flat
   over a path, on the W5.12 row's pattern.
5. `test_the_realised_population_agrees_with_the_numpy_path` — torch and
   jax fixtures, at points including near a support boundary
   (`tolerances.cross_backend`); and the element arrives in the native array
   type (no coercion).
6. `test_refusals_by_name` — a path that does not resolve (names the
   available components at that level); a path on a plain-set component; a
   leaf the composite does not declare; a tied leaf; a composite without a
   path (the W5.12 message, re-worded to name the remedy).
7. `test_the_factory_writes_the_paths` — `DatasetCollection.plate(within=)`
   produces the same `Population` as the explicit `over`.
8. `test_the_plate_coordinate_is_the_dataset_labels` — the emitted run's
   `gp.amplitude` dimension is labelled by the dataset labels (dev only,
   through `ampere.results`).
9. `test_the_population_is_in_provenance` — the problem record carries the
   declaration and the hash moves when `over` changes.

### 9.2 `test_derived.py`, new

1. `test_a_derived_parameter_is_not_a_free_dimension` — `free_size`,
   `free_labels`, `pack`/`unpack` exclude it; `names` and
   `evaluation_order` include it after its inputs.
2. `test_complete_computes_it_and_refuses_it_given` — both halves of the
   §3.3 row.
3. `test_the_grammar_is_closed` — each refused node kind refused by name
   (attribute, subscript, comparison, an unlisted call, a lambda).
4. `test_the_spec_round_trips_and_hashes_stably` — `to_spec` → `from_spec`
   → equal set; the normalised source identical for `a+b` and `a + b`.
5. `test_a_hierarchical_prior_may_reference_a_derived_parameter` —
   `lnprior` equals the hand-computed density with the derived scale.
6. `test_the_slab_is_piironen_and_vehtari` — `shrinkage_horseshoe(tail=
   "slab", slab_scale=c)` on three amplitudes: the effective local scale
   equals `sqrt(c²τ²λ²/(c²+τ²λ²))` at sampled `(τ, λ)`; the plain and
   regularised tails unchanged byte for byte (the existing rows still pass).
7. `test_the_non_centred_population_declares_the_centred_density` —
   `Population` with internal `z` and derived `theta`: `theta_i = mu +
   sigma z_i`; the joint log prior over `(mu, sigma, z)` plus the
   change-of-variables-free identity equals the centred declaration's
   marginal check by Monte Carlo at a loose tolerance; `z` reaches no
   member; an undeclared non-input member is refused at merge.
8. `test_the_realised_derived_agrees_with_the_numpy_path` — torch and jax:
   `log_prob_unconstrained` agreement at `tolerances.cross_backend` on a
   problem whose hierarchical prior references a derived parameter and
   whose model consumes a derived element; the jax sample stream carries the
   deterministic site under the merged name.
9. `test_the_posterior_carries_the_derived_variable` — an emcee run's
   `posterior` has the derived variable, equal to `complete()` of each draw;
   `ampere_derived` lists it; schema version 10.
10. `test_the_evaluation_order_is_honoured_under_tracing` — a derived
    parameter referencing a hierarchical one referencing a derived one:
    identical on numpy, torch and jax (the §6 hoisting hazard).

## 10. The implementation items, drafted for D8's ruling (Phase 7)

### W7.0 — The `Derived` parameter node [M; Opus]
§3 of this memo, whole: `Derived` in the prior slot with the closed grammar;
the state, the §3.3 consumer table, `ArrayOps.sqrt`/`log` on the three
namespaces; the lowering (jax `deterministic`, both `_resolved` loops, the
provenance row kind); `_posterior_variables` and `ampere_derived`,
`PROVENANCE_SCHEMA_VERSION` → 10; the internal-member rule on `Population`
(§3.5) with its merge-time refusal; `shrinkage_horseshoe(tail="slab",
slab_scale=)`; the contract amendments of §3.8 and the decision-log row of
§8. **Depends:** nothing. **Accept:** §9.2's ten rows green on every fixture;
the existing `test_population.py`, shrinkage and horseshoe rows unchanged
and green; a NUTS run on torch and jax of the non-centred fifty-member
population (`tests/inference/test_population_nuts.py`'s problem re-declared
with `z` and a derived `theta`) recovering `mu` and `sigma` inside the
central 95 % with fewer divergences than the centred declaration at the
same budget (pinned as an inequality); lint/format/pyrefly clean; gates dev
+ torch + jax.

### W7.1 — Populations over a qualified component path [M; Opus]
§2 of this memo, whole: the `over` grammar and its resolution through the
retained inner mappings; `PlateBinding`/`Binding` qualified local paths; the
outer-declaration strip and the flat re-prior on a leaf; the new leaf-must-
exist refusal and the re-worded composite refusal;
`DatasetCollection.plate(within=)`; `populations` in the problem's
provenance record; the contract amendments of §2.6 and the decision-log row
of §8. **Depends:** W7.0 (the internal-member rule and the refusal-at-merge
it brings; not the `Derived` type itself). **Accept:** §9.1's nine rows
green on every fixture; `test_population.py`'s existing rows unchanged; the
two-level (instrument-step) path exercised; `inference.md` §9's doctest of
the plate-of-datasets rewritten to the path form and executing; gates dev +
torch + jax.

### W7.2 — The per-dataset GP amplitude as a population: the M2 validation [M; Opus]
§4's customer, measured. An `examples/m2_misspecification` scenario — or a
sibling under `examples/` if the study's `SCENARIOS` must stay as the
milestone fixed them (W5.8's precedent) — with N spectra of one object class
sharing a population of GP amplitudes, non-centred (`z`, derived
`amplitude`), against (a) independent amplitudes and (b) one tied amplitude:
bias, calibration and localisation on the M2 pattern, the shrinkage of a
poorly constrained member's amplitude towards the population pinned as a
margin, SBC over refits with the deviation injected, NUTS on torch and jax
through the realisation, emcee on the reference path. A `docs/source`
section on `population.rst` teaching the three declarations (tie, flat,
non-centred population) and when each is right. **Depends:** W7.0, W7.1.
**Accept:** the scenario in the driver with its `tests/m2` rows and pinned
margins; the SBC contrast pinned; the docs section; gates dev + torch.

---

*Opus review: pending. Findings are appended here when the review lands.*
