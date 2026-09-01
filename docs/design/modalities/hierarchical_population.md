# Modality sketch (f) — hierarchical population inference

Status: **DRAFT for Peter's review.** W1.11 sketch (f). Validates
`DEVELOPMENT_PLAN.md` §4.1 (plate-aware groups), §4.5 and design horizon (d)
against the merged contracts, and is the stress test the W1.11 dispatch asks
for: **does flat, single-call parameter merging scale to populations, or does
the nested-merge question Peter expressly kept open (`parameters.md` §14
preamble) bite here?**

The short answer, argued from measurements in §5 and §6: **the nested-merge
question does not bite — nesting works for routing and is lossy only in the
bindings — but a different limitation bites hard, and it is the one that decides
whether populations are usable at all.** `Plate` vectorises the population and
lowers to a numpyro plate but cannot be routed to per-object datasets;
`HierarchicalPrior` + tying routes to per-object datasets but costs one
`Parameter` per object, which is **400× slower per prior evaluation** at
N = 1000. Neither construct serves the modality, and the missing piece is a
single index on `Binding`.

Every measurement below was taken against the merged `ampere.core` at commit
`8c4e99d`; §8 lists what was run.

---

## 1. The modality

A survey yields N objects of one class — say N post-AGB stars, each with a
photometric SED and one or two spectra. Each object has its own dust
temperature, optical depth and distance; the *population* has a mean and a
spread in each, and the population parameters are what the science is about.
Each object also has its own instrumental nuisance parameters (a calibration
scale per spectrum), which are emphatically not shared.

In contract terms:

```
DatasetCollection                              (W1.7)
├── population hyperparameters  mu, sigma      shared by every object
├── obj0 : Dataset(observed, Instrument, Likelihood)   model params theta_0 ~ N(mu, sigma)
├── obj1 : Dataset(...)                                theta_1 ~ N(mu, sigma)
│   ...
└── objN : Dataset(...)                                theta_N ~ N(mu, sigma)
```

`parameters.md` §9 already anticipates exactly this and names it *"N components,
each with its own dataset … the shape a `DatasetCollection` (W1.7) of per-object
datasets naturally takes"*. This sketch takes that recommendation at its word
and pushes on it.

`DEVELOPMENT_PLAN.md`'s design horizon (b) — populations by post-processing,
via importance reweighting of archived single-object fits — is the *other*
answer to the same science question and is unaffected by anything here; it needs
only the stored per-sample `log_likelihood`/`log_prior` §4.6 already mandates.
This sketch is about the joint fit.

## 2. The two constructs, and which layout each serves

`parameters.md` §9 offers two, and is explicit that they serve different data
layouts:

| Construct | Population lives as | Lowers to | Per-object dataset? |
|---|---|---|---|
| `Plate` | one array-valued `Parameter` of shape `(N,)` | one `numpyro.sample` inside one `numpyro.plate` | **no** (§6) |
| `HierarchicalPrior` + `shared_as` tying | N scalar `Parameter`s, one per component | N separate sample sites | yes |

Both were exercised. The rest of this document is about the fact that the
modality needs the first construct's *cost* and the second construct's
*routing*, and no combination of the two delivers both.

## 3. The composition that works today

Per-object parameter set, following `parameters.md` §9 pattern 2:

```python
def object_set():
    return ParameterSet([
        Parameter("mu", st.norm(0.0, 5.0), shared_as="mu"),
        Parameter("sigma", st.halfnorm(0.0, 2.0), shared_as="sigma"),
        Parameter("theta", HierarchicalPrior("norm", {"loc": "mu", "scale": "sigma"})),
    ])

pop = ParameterSet.merge({f"obj{i}": object_set() for i in range(4)})
# pop.merged.names -> ('mu', 'sigma', 'obj0.theta', 'obj1.theta', 'obj2.theta', 'obj3.theta')
# pop.merged.free_size -> 6          (2 hyperparameters + 4 objects, not 12)
# pop.merged["obj0.theta"].references -> ('mu', 'sigma')
```

The tie collapses the N duplicate `mu` declarations into one merged parameter,
the hierarchical references follow the collapse, and the sampler sees N + 2
dimensions. That much is exactly as advertised.

### The hyperprior can be declared once, and should be

Declaring `st.norm(0.0, 5.0)` N times and relying on the tie to notice they
agree is fragile: a typo in object 137 raises a `TyingError` at composition, which
is at least loud, but there is no reason to write it N times at all.
`parameters.md` §8's *deferred* declaration is the better idiom, and it composes
(verified):

```python
sets = {"population": ParameterSet([
            Parameter("mu", st.norm(0.0, 5.0), shared_as="mu"),
            Parameter("sigma", st.halfnorm(0.0, 2.0), shared_as="sigma"),
        ])}
sets.update({f"obj{i}": ParameterSet([
            Parameter("mu", shared_as="mu"),            # deferred: no prior of its own
            Parameter("sigma", shared_as="sigma"),
            Parameter("theta", HierarchicalPrior("norm", {"loc": "mu", "scale": "sigma"})),
        ]) for i in range(4)})

pop = ParameterSet.merge(sets)
# names -> ('mu', 'sigma', 'obj0.theta', ..., 'obj3.theta');  free_size -> 6;  is_resolved -> True
```

The hyperprior has exactly one authority, each object carries a stub, and the
whole thing is `is_resolved` after the merge and refuses to evaluate priors
before it (`parameters.md` §8's own behaviour, applied at population scale).
**Recommendation: W1.7's `DatasetCollection` should build this shape, not the
N-duplicate one**, and the population component should be a first-class thing it
owns rather than a dictionary key the user remembers to add.

### The consequence nobody will expect

`ParameterSet` requires a hierarchical prior's references to resolve *within the
set that will evaluate them*:

```
ParameterError: parameter 'theta' has a hierarchical prior referencing 'mu',
which is not in this set (it declares ['theta']). Hierarchical references must
resolve within the set that will evaluate them.
```

So every per-object component must declare `mu` and `sigma`, if only as deferred
stubs — and `distribute` therefore hands every object's model values for them:

```python
routed = pop.distribute(theta)
sorted(routed["obj0"])          # -> ['mu', 'sigma', 'theta']
```

and `Model.__call__` **refuses** values it did not declare:

```
ParameterError: got value(s) for unknown parameter(s) ['mu', 'sigma'];
this set declares ['theta'].
```

Therefore: **a per-object `Model` written for a single object cannot be dropped
into a population fit unchanged.** It must declare the population
hyperparameters it does not use. That is a real composability cost and it is
the first half of gap H-1.

## 4. Gap H-1 in full — hierarchy is declaration-time only

`parameters.md` §8 makes a point of supporting tying in *both* directions,
because "the models were written independently — or came from a library — and
the person composing the fit is the one who knows". For **hierarchy, only one
direction exists.** There is no composition-time analogue of `Tie`:

- `Tie` collapses N sites into **one** parameter. A population needs them to
  stay N, drawn from a shared prior. `Tie(..., prior=...)` overrides the
  *collapsed* parameter's prior, which is a different thing entirely — using it
  here would fit one θ for the whole survey, silently and plausibly.
- `HierarchicalPrior` must be present on the site at declaration, and its
  references must already resolve locally (§3).

So a user with N library models has no contract-level way to say "these are
draws from a population". They must rewrite each component's parameter set.

That rewrite is at least *possible* in user code, which caps the severity. A
~15-line helper — inject deferred hyperparameter stubs, `release()` the target
parameter with a `HierarchicalPrior`, add a `population` component — was written
and run against the merged code:

```python
m = populate(plain_sets, target="theta", family="norm",
             mapping={"loc": "mu", "scale": "sigma"},
             hyperparameters=[Parameter("mu", st.norm(0.0, 5.0)),
                              Parameter("sigma", st.halfnorm(0.0, 2.0))])
# merged -> ('mu', 'sigma', 'obj0.theta', 'obj0.scale', 'obj1.theta', ...)
# free_size 8, is_resolved True, lnprior finite, obj0.theta.references ('mu', 'sigma')
```

Note that `obj0.scale` — the per-object calibration nuisance — passes straight
through untouched, which is right: population structure applies to one named
quantity, not to the whole component.

So the amendment is a *contract-level declaration for something already
expressible*, not new machinery. See §9.

## 5. The scaling stress test

This is the measurement the dispatch asked for. `parameters.md` §9 recommends
the N-component pattern for per-object datasets; here is what it costs.
Timings are single-threaded numpy/scipy on the pixi `dev` environment, and the
column that matters is `lnprior`, because it runs **once per log-probability
evaluation**:

| N objects | `free_size` | `merge` (once) | `lnprior` (per evaluation) | `free_labels` | `distribute` |
|---|---|---|---|---|---|
| 100 | 102 | 83 ms | **46 ms** | 0.02 ms | 0.10 ms |
| 1000 | 1002 | 703 ms | **373 ms** | 0.09 ms | 0.73 ms |

373 ms per prior evaluation means a 10⁵-sample run spends **ten hours** in
`lnprior` alone, before the model has been evaluated once. At N = 100 it is
46 ms, still an order of magnitude more than a typical SED model evaluation.
The N-component pattern is unusable beyond a few tens of objects.

The same population declared as a `Plate`:

| N objects | `free_size` | build | `lnprior` | `free_labels` | `len(pset)` |
|---|---|---|---|---|---|
| 100 | 102 | 0.05 ms | **0.68 ms** | 0.13 ms | 3 |
| 1000 | 1002 | 0.04 ms | **0.61 ms** | 1.28 ms | 3 |

Roughly **600× faster** at N = 1000, and flat in N.

### The cost is not the merge, and not the hierarchy

This matters for choosing the fix, so it was isolated. Three sets, all with
1002-ish free dimensions:

| Declaration | `lnprior` |
|---|---|
| 1000 plain scalar `Parameter`s, `st.norm(0, 1)` each | 190 ms |
| 1000 scalar `Parameter`s with `HierarchicalPrior` | 509 ms |
| **1 array-valued `Parameter` of shape `(1000,)` with a `HierarchicalPrior`** | **1.28 ms** |

So: the hierarchy roughly triples the per-parameter cost (topological ordering
plus a `bind()` per site), but the dominant term is simply **the number of
`Parameter` objects**, each of which costs a Python-level `scipy` call. One
array-valued parameter is 150× faster than a thousand scalars carrying the same
information, because scipy vectorises across the elements.

`Plate` already produces exactly that array-valued parameter — `parameters.md`
§9 says so, and says why: *"That shape is chosen because it is the numpyro
lowering."* The vectorisation is not an optimisation to be added later; it is
already there and the modality cannot reach it.

**Nothing in this table is about `merge`.** Merging is a one-off cost (703 ms at
N = 1000 is tolerable at setup). The per-evaluation cost is a property of the
*declaration*, and would be identical without any merging at all.

## 6. Gap H-2 — a `Plate` member cannot be routed to per-object datasets

This is the headline gap, and it is exactly the seam `parameters.md` §12.3
predicts: *"cross-component plate identity needs an index-alignment concept that
belongs with W1.7's `DatasetCollection`."*

A plate expands to one array-valued parameter, and `merge`/`distribute` route
whole parameters to components. Verified:

```python
plate = Plate("objects", size=4,
              hyperparameters=[Parameter("mu", st.norm(0.0, 5.0)),
                               Parameter("sigma", st.halfnorm(0.0, 2.0))],
              members=[Parameter("theta", HierarchicalPrior("norm",
                                                            {"loc": "mu", "scale": "sigma"}))])
survey = ParameterSet([], plates=[plate])
m = ParameterSet.merge({"survey": survey})

[(b.global_name, b.component, b.local_name) for b in m.bindings]
# [('survey.objects.mu', 'survey', 'objects.mu'),
#  ('survey.objects.sigma', 'survey', 'objects.sigma'),
#  ('survey.objects.theta', 'survey', 'objects.theta')]

m.distribute(theta)["survey"]["objects.theta"]
# array([...])          <- the WHOLE (4,) block, to ONE component
```

`Binding` is `(global_name, component, local_name)` — there is no index — so
"object 1's dataset consumes element 1 of `objects.theta`" is not expressible.
`free_labels()` already names the elements (`'objects.theta[0]'`, …), so the
information exists; only the binding does not.

The two workarounds, and why neither is acceptable:

1. **Tie the plate members across components.** Refused, explicitly:

   ```
   TyingError: tie 'theta' spans plate members (objects.theta is in plate
   'objects'); tying across plates is out of scope for v1.3. Express the shared
   quantity as a hyperparameter referenced by a HierarchicalPrior instead.
   ```

   The advice in that message is sound for a *shared* quantity and does not
   apply here: θ_i is not shared, it is the i-th draw.

2. **Give every per-object model the whole `(N,)` array and let it index by
   position.** This works today (verified: the routed value is the full block)
   and is what a determined user would do. It makes every per-object model
   depend on N, which defeats the modularity `DatasetCollection` exists for; it
   hands each of 10³ models 10³ values of which it uses one; and it puts the
   object's identity in an integer offset the model carries by hand, which is
   precisely the positional-`theta`-slicing idiom `parameters.md` §1 exists to
   kill.

## 7. The nested-merge question, answered

Peter kept `parameters.md` §14 question 2 expressly open and asked W1.7 to
evaluate a nested `ParameterMapping` on its merits. This modality is the case
that would need it — three levels of structure (population → object →
{model, instrument, likelihood}) — so here is the evidence, in three parts.

### (a) A flat single merge cannot express three levels

Component labels must be plain identifiers:

```
ParameterError: component label name 'obj1.model' must not contain '.':
qualified names are produced by ParameterSet.merge() and Plate expansion,
never declared directly.
```

So the "flatten it all into one call" route requires labels like `obj1_model`,
`obj1_instr`, `obj1_like` — 3N components whose structure is encoded in a naming
convention rather than in the data, and which then appear in every ArviZ
coordinate and provenance record as `obj1_instr.calibration_scale.scale`. That
is expressible and ugly.

### (b) Nesting works for routing, and loses only the bindings

Merging an already-merged set does **not** fail, and `distribute` round-trips
correctly through two levels — `transformations.md` §14's claim ("each level
re-distributes") is correct. What is lost is the inner *wiring*:

```python
inner = ParameterSet.merge({"model": ..., "instr": ...})     # 'model.t' tied to 'instr.t'
# inner.merged.names -> ('t', 'instr.scale');  inner.bindings -> 3;  inner.tied_names -> ('t',)

outer = ParameterSet.merge({"obj0": inner.merged, "obj1": inner.merged})
# outer.merged.names -> ('obj0.t', 'obj0.instr.scale', 'obj1.t', 'obj1.instr.scale')
# outer.bindings     -> [('obj0.t', 'obj0', 't'), ('obj0.instr.scale', 'obj0', 'instr.scale'), ...]
```

The outer mapping records four bindings, one per merged name. The inner fact
that `obj0.t` feeds *two* sites — the model's `t` and the instrument's `t` — is
gone: `outer.tied_names` is empty, and any consumer that walks `outer.bindings`
to find shared quantities (provenance, an ArviZ coordinate builder, a
plate-lowering pass) sees none. Routing survives because the caller re-runs
`inner.distribute` on the outer result; *introspection* does not.

That is a weaker statement than `parameters.md` §12.4's "silently drops the
first merge's bindings" implies for routing, and exactly as strong as it implies
for structure.

### (c) Cross-level ties do work

`transformations.md` §14 flags this as uncovered and asks W1.7 to say what
happens. Verified — it works. Taking a plain inner merge this time (a model
parameter `t` and an instrument parameter `scale`, no inner tie), and tying the
two objects' calibration scales at the outer level:

```python
inner = ParameterSet.merge({"model": ParameterSet([Parameter("t", st.norm(0, 1))]),
                            "instr": ParameterSet([Parameter("scale", st.norm(1, 0.1))])})
outer = ParameterSet.merge({"obj0": inner.merged, "obj1": inner.merged},
                           ties=[Tie("cal", ("obj0.instr.scale", "obj1.instr.scale"))])
# outer.merged.names -> ('obj0.model.t', 'cal', 'obj1.model.t');  free_size -> 3
# outer.distribute(...) -> {'obj0': {'model.t':…, 'instr.scale':…}, 'obj1': {...}}
# inner.distribute(outer_result['obj0']) -> {'model': {'t': …}, 'instr': {'scale': …}}
```

A composition-time tie across two objects' instruments collapses correctly, and
the two-level re-distribution routes it to both. So the answer to W1.5 §14's
open worry is: **cross-level ties are expressible and correct.**

### Verdict on recursive merge

**Recursive merge is not what this modality needs, and adding it would not fix
anything measured here.** The pressure the population case applies is entirely
about §6 — per-element routing of a vectorised parameter — and about §5's cost,
which is a property of the declaration and unaffected by how sets are combined.
Nesting is already usable for routing; its only real deficiency is the lost
bindings in (b), which is a *reporting* problem with a much cheaper fix than
recursion (gap H-3).

If Peter wants the option kept open anyway, the honest framing is: nesting
should be made **lossless** rather than **recursive**. Those are different
changes, and the first is a few lines.

## 8. What was verified, and how

Executed against `ampere.core` at `8c4e99d` in the pixi `dev` environment.

| Claim | Evidence |
|---|---|
| §9 pattern 2 merges N components to N + 2 dimensions; references follow the collapse | ran at N = 4; `free_size` 6, `references` `('mu','sigma')` |
| The deferred-hyperprior idiom composes and resolves | ran at N = 4; `is_resolved` `True`, same 6 dimensions |
| Hierarchical references must resolve within the component's own set | `ParameterError` quoted in §3 |
| `distribute` hands every object `mu`, `sigma`, `theta`; `Model.__call__` refuses undeclared values | both errors quoted in §3 |
| Timings, N-component pattern and `Plate` | §5's two tables; `time.perf_counter`, single run each |
| The cost is the number of `Parameter` objects, not the merge or the hierarchy | §5's third table |
| A composition-time population rewrite is expressible in user code | ~15-line `populate()` helper; merged, resolved, finite `lnprior` |
| `Plate` routes the whole `(N,)` block to one component; `Binding` has no index | §6, bindings printed |
| Tying across plate members is refused | `TyingError` quoted in §6 |
| Dotted component labels are refused | `ParameterError` quoted in §7(a) |
| Nested merge loses the inner bindings but not the routing | §7(b) |
| Cross-level `Tie` works through two levels | §7(c) |
| Dotted *channel* names are refused | `SchemaError`; see gap H-5 |
| Nested plates are not expressible as plates | §10 question 3 |

## 9. Interface gaps

The composition **does not close cleanly**. It closes structurally — every piece
composes, and a population fit of a few dozen objects would run today — but the
modality the plan scopes (a survey, so hundreds to thousands of objects) is
blocked by H-2, and H-2 is not a performance problem that Phase 5 can optimise
away: it is a missing binding.

The gaps are listed in severity order rather than in numerical order; H-2 is the
one on the critical path.

### H-2 — `Binding` cannot address one element of an array-valued parameter

**Severity: blocking for the modality. Should land in the freeze**, because it
changes a frozen dataclass and `distribute`'s contract.

`Plate` produces the vectorised declaration the population needs and the numpyro
lowering W1.9 wants; `DatasetCollection` needs per-object scalars. One index
closes the gap.

**Proposed amendment** — `parameters.md` §7 (the `merge` section), §12.3, and
`Binding`/`ParameterMapping.distribute`:

> `Binding` gains an optional element index:
>
> ```python
> @dataclasses.dataclass(frozen=True)
> class Binding:
>     global_name: str
>     component: str
>     local_name: str
>     index: int | tuple[int, ...] | None = None
> ```
>
> When `index` is not `None`, `distribute` routes `resolved[global_name][index]`
> rather than the whole array, so a component consuming one member of a plate
> receives a scalar under its own local name. The merged set is unchanged — one
> array-valued parameter, one plate, one numpyro sample site — so the lowering
> in `lowering.md` is unaffected; only the routing changes.
>
> `ParameterSet.merge` gains a way to declare such bindings. The minimal form is
> an explicit argument:
>
> ```python
> ParameterSet.merge(
>     {"population": hyper_set, "obj0": obj0_set, ...},
>     plate_bindings=[PlateBinding("objects.theta", component="obj0", local_name="theta", index=0), ...],
> )
> ```
>
> and the convenient form is for W1.7's `DatasetCollection` to construct them
> from the order of its datasets, which is the natural index.
>
> Limitation 3 ("tying across plate members is refused") is unaffected and
> should stay: this is not tying — each member remains its own draw — it is
> *addressing*, which is what the limitation says is missing.

This also removes limitation 12.3's stated blocker, so §12.3's text should be
rewritten from "needs an index-alignment concept that belongs with W1.7" to
"is provided by `Binding.index`; W1.7 constructs them".

### H-1 — hierarchy has no composition-time declaration

**Severity: real composability cost; a convenience amendment, could follow the
freeze.** §4 has the argument. A user with independently written per-object
models must rewrite their parameter sets to introduce a population; `Tie` covers
the analogous case for sharing and has no counterpart here.

**Proposed amendment** — `parameters.md` §8 and §9:

> Alongside `Tie`, `merge` accepts a composition-time population declaration:
>
> ```python
> @dataclasses.dataclass(frozen=True)
> class Population:
>     """N components' same-named parameter are draws from one shared prior."""
>     parameter: str                      # local name in each component
>     family: str                         # neutral prior family, as HierarchicalPrior
>     hyperparameters: Mapping[str, Parameter]   # e.g. {"loc": Parameter("mu", ...)}
>     components: Sequence[str]
> ```
>
> `merge` then adds the hyperparameters to the merged set under their own names,
> rewrites each named component's parameter to carry the corresponding
> `HierarchicalPrior`, and records the bindings. It is the hierarchical
> counterpart of `Tie` and exists for the same reason §8 gives for `Tie`: a
> library model cannot anticipate the fit it will be composed into, and a
> population is a statement the person composing the fit makes.
>
> Note the consequence this removes: without it, every per-object model must
> declare the population hyperparameters it never uses, because a
> `HierarchicalPrior`'s references must resolve within the component's own set.

### H-3 — a nested merge silently loses the inner bindings

**Severity: reporting; cheap.** §7(b). Routing is fine; introspection is not, and
provenance, ArviZ coordinate naming and any future plate-lowering pass all walk
the bindings.

**Proposed amendment** — `parameters.md` §12.4, replacing "merge is not
associative":

> `merge` is not associative in its **bindings**. Merging a
> `ParameterMapping.merged` again routes correctly — each level re-distributes,
> and cross-level `Tie`s work — but the outer `ParameterMapping` records one
> binding per merged name and therefore cannot see the inner merge's tie
> structure: `tied_names` under-reports, and a consumer walking `bindings` to
> find shared quantities finds none.
>
> Merging all components in one call remains the recommendation. Where nesting
> is unavoidable (W1.7's `DatasetCollection` of `Dataset`s of instruments), a
> caller that needs the full structure should keep the inner `ParameterMapping`s
> alongside the outer one. The extension point, if the reporting gap bites, is
> for `merge` to accept `ParameterMapping`s as well as `ParameterSet`s and
> compose their bindings — which is *lossless nesting*, and a strictly smaller
> change than recursive merge.

### H-4 — `lnprior` is O(number of `Parameter` objects), which the recommended population idiom maximises

**Severity: documentation now, implementation later; no interface change.** §5's
numbers. This is not a defect — `ParameterSet` is a declaration container, not a
hot-loop object — but `parameters.md` §9 currently recommends the N-component
pattern for exactly the layout in which it is 600× slower than the alternative,
with no warning.

**Proposed amendment** — `parameters.md` §9, "Which construct to use", add:

> **Cost.** The second pattern declares one `Parameter` per member, and
> `lnprior` evaluates one scipy call per `Parameter`: measured at 46 ms for 100
> members and 373 ms for 1000, *per log-probability evaluation*. The first
> pattern declares one array-valued parameter and costs ~1 ms at either size,
> because scipy vectorises across the elements. For more than a few tens of
> members, use a `Plate` — and, once `Binding.index` exists (§7), a `Plate` can
> feed per-object datasets too, so the two patterns' data layouts stop being a
> reason to choose the slow one.

### H-5 — channel names cannot be qualified, which blocks the one-model population layout

**Severity: one line; routed here by `results_schema.md` §17.** The other
population layout — a single model that evaluates all N objects and emits one
channel per object — needs channel names like `obj1.sed`. Refused:

```
SchemaError: channel name 'obj1.sed' is not usable: a channel name must be a
valid Python identifier … Rename the channel, e.g. 'obj1_sed'.
```

`results_schema.md` §17's added note already proposes the fix (accept
dot-separated identifiers, as parameter names do) and routes the decision to
W1.7, pairing it with the nested-merge question. This sketch's answer to the
pairing: **the two are not symmetrical and should be decided separately.**
Qualified channel names are a one-line relaxation of `_check_channel_name` with
no structural consequence; nested merge is a structural question, and §7 says it
should be answered "make nesting lossless", not "make it recursive".

**Proposed amendment** — `results_schema.md` §2 and `_check_channel_name`:

> A channel name must be a `.`-separated sequence of Python identifiers, exactly
> as a parameter name is (`parameters.md` §3). Qualification is how a model that
> emits per-object output namespaces its channels (`obj1.sed`), mirroring how
> `ParameterSet.merge` qualifies parameter names; nothing in `ModelResult`
> interprets the dots, and `require` still matches the full name.

## 10. Requirements on W1.7

1. **`DatasetCollection` owns the population component.** It should construct
   the deferred-stub shape of §3 (one authority for each hyperprior, a stub per
   object) rather than leaving the user to remember it, and should expose the
   hyperparameters as a named group so W1.8 can label them separately from the
   per-object block.
2. **`DatasetCollection` constructs the `Binding.index` wiring** (gap H-2) from
   its own dataset ordering. That ordering then becomes part of the run's
   provenance, because it defines which posterior element belongs to which
   object — W1.8 must record the dataset labels as the plate's coordinate, not
   an integer range.
3. **One merge call, over datasets, with the instrument mappings kept.**
   `likelihoods.md` §5 and `parameters.md` §12.4 both require a single merge;
   §7(b) shows what is lost if `Dataset.parameters` is itself a merge. W1.7
   should either merge model + instrument + likelihood + datasets in one call
   (with flat labels), or keep every inner `ParameterMapping` and document that
   `distribute` is applied level by level.
4. **The hyperprior extension point W1.7's item text promises is `Population`
   (H-1) plus `Binding.index` (H-2)**, not a new mechanism. Both live in
   `parameters.md`; `DatasetCollection`'s job is to build them, not to
   re-implement them.
5. **`free_labels()` at population scale.** `likelihoods.md` §16 already tells
   W1.8 that 10⁵ scalar names is the wrong ArviZ representation for a latent
   block. The same applies to a plate member: `objects.theta` should be one
   dimension with the dataset labels as its coordinate. Measured cost is small
   (1.3 ms at N = 1000) so this is about representation, not speed.
6. **Failure signalling is per object.** One object's simulator failing should
   not lose the whole population draw silently. §4.5's "flagged failed
   simulations" needs a per-dataset granularity here, and the recorded reason
   should name the dataset.

## 11. Open questions for review

1. **Is `Binding.index` (H-2) the right shape, or should `DatasetCollection`
   slice before it calls the model?** The alternative to amending `Binding` is
   for W1.7 to hold the plate itself, distribute the whole array to itself, and
   hand each dataset its slice — no core change at all. That works, but it puts
   parameter routing in two places (some through `ParameterMapping`, some
   through `DatasetCollection`), and every other consumer that walks the
   bindings — provenance, ArviZ labelling, W1.9's lowering — would see an
   incomplete picture. This sketch recommends the amendment; the alternative is
   cheaper today and more expensive later.
2. **Should `Population` (H-1) land before the freeze or with Phase 5?** It
   changes `merge`'s signature, which argues for the freeze; nothing depends on
   it before Phase 5's hierarchical implementation, which argues against. The
   H-2 amendment is the one that must not wait.
3. **Nested plates.** `parameters.md` §12.2 defers "objects within surveys".
   Two surveys each holding a plate merge fine — `('survey0.objects.mu',
   'survey0.objects.theta', 'survey1.objects.mu', …)` — but `Parameter.plate` is
   a single label, so the outer level is not tagged as a plate and W1.9 cannot
   lower it to nested numpyro plates. The modality does not need it yet; the
   extension is §12.2's own (`Parameter.plate` becomes a tuple), and it should
   be decided together with H-2 since both are about plate addressing.
4. **Per-object nuisance parameters as plate members.** A `Plate` member with an
   ordinary prior and no hyperparameters is an i.i.d. array (verified:
   `objects.cal` of shape `(5,)`), which is exactly right for "each object has
   its own calibration scale, from a common prior". With H-2 this vectorises the
   nuisance parameters too. Worth confirming that is the intended reading of
   §9's `members` argument, since the spec's example only shows hierarchical
   members.
5. **Does design horizon (b) make any of this less urgent?** Importance
   reweighting of archived single-object fits gets population hyperparameters
   without a joint fit, needs nothing beyond the stored per-sample
   `log_likelihood`/`log_prior`, and scales to any N. If Peter regards it as the
   primary route for large N, H-2 becomes a Phase 5 convenience rather than a
   freeze item. If the joint fit is the primary route, H-2 is on the critical
   path. The contracts should record which.
