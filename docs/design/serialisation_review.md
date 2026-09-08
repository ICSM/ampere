# Ampere v2 — the consolidated serialisation review (W1.13)

Status: **the freeze's one serialisation pass.** Ruled 2026-09-03: instead of
settling serialisation piecemeal per contract, one review inventories every
mechanism, closes the questions routed here — `likelihoods.md` §17 Q8,
`results_schema.md` §17 Q6 (with design-horizon (c)'s training-set
requirement), `results.md` §15 R7 — and states one coherent approach the
Phase 2 backends implement against. Everything below was verified against the
code at the freeze; the inventory tables name files so drift is checkable.

---

## 1. The rule, stated once

Ampere has exactly **three serialisation layers**, and each thing belongs to
one of them:

1. **Declarative specs live on the objects, in `ampere.core`.**
   `to_spec()` / `to_dict()` / `spec().to_dict()`: plain data, JSON-able,
   round-trippable where a round trip is declared. A spec describes *the
   declaration* — names, families, priors, shapes, units, configuration —
   and never bulk content. `ParameterSet.to_spec` ↔ `from_spec`,
   `PriorSpec.to_dict` ↔ `from_dict`, `HierarchicalPrior.to_dict` ↔
   `from_dict`, `KernelSpec.to_dict` (one-way, deliberately: the values are
   `Parameter`s with their own contract), and — landed with this review —
   `Likelihood.to_spec()`.
2. **Content fingerprints live in `ampere.results.provenance`.** One-way by
   design: `hash_array`, `buffer_fingerprint`, `container_fingerprint`,
   `describe_likelihood`, `dataset_fingerprint`, `model_fingerprint`,
   `problem_fingerprint`, `spec_hashes`, `provenance_attrs`. They **compose**
   the layer-1 specs and add content hashes for what specs deliberately omit
   (buffer bytes, data arrays, censoring code positions). This is the
   cache-key and reproducibility layer `DEVELOPMENT_PLAN.md` §7 demands;
   nothing is ever reconstructed from it.
3. **Bulk storage is netCDF, through `ampere.results`.** Containers and
   `ModelResult`s serialise by value through the *functions* in
   `ampere.results.serialisation` (ruled R5: functions, not methods — a
   user-defined kind gets serialisation for free); whole runs go through
   `emission.emit` → `to_netcdf` ↔ `from_netcdf` (returning a `DataTree`;
   there is deliberately no `problem_from_netcdf` — a stored run *records* a
   problem, it does not reconstitute one).

The boundary rule that keeps the three honest: **a spec describes the
declaration; provenance fingerprints content; storage carries values.** Two
likelihoods differing only in kernel family or solver jitter are
distinguishable at layer 1 (that is R7's correctness argument); two runs
differing only in an opacity table are distinguishable at layer 2; and the
data themselves live only at layer 3.

## 2. What this review landed

- **`Likelihood.to_spec()`** (`ampere/core/likelihood.py`) — R7 confirmed.
  The declarative mapping `describe_likelihood` used to assemble from
  outside: family name and class, noise-model class, marginalisation,
  `parameters.to_spec()`, and for a GP the kernel spec plus the solver's
  name/class/exactness/dataclass configuration (`DenseGP.jitter` changes
  the number the same θ scores, so it is identity). Censoring appears as
  counts only — code positions are content. `describe_likelihood` now
  *composes* it and adds the fingerprints, so backends and the conformance
  suite share one definition.
- **A provenance hole, fixed**: `Likelihood.__init__` forwards *parameters*
  from its family and noise model, never buffers, so
  `buffer_fingerprint(likelihood)` was blind to a family's background
  template or a noise model's tabulated data — precisely the stale-cache
  trap §7 warns about. `describe_likelihood` now fingerprints the family's
  and the noise model's buffers too; `PROVENANCE_SCHEMA_VERSION` bumped
  to 2.

## 3. Inventory and dispositions

| Mechanism (file) | Round trip? | Disposition |
|---|---|---|
| `PriorSpec.to_dict/from_dict`; `describe_prior` ↔ `prior_from_spec` (`core/parameter.py`) | yes, tested | Sound. Canonical keyword-only form is the single thing W1.9 lowers. `discrete` field is the queryability the 2026-09-03 discrete ruling requires |
| `ParameterSet.to_spec/from_spec` (`core/parameter.py`) | yes, tested (unit + conformance) | Sound; order-significant, which the spec hash depends on. Inferred bijections deliberately unwritten (a function of the prior); custom bijections/opaque priors refuse loudly rather than serialise wrongly |
| `HierarchicalPrior.to_dict/from_dict` | yes | Sound; references by name, resolved in the enclosing set |
| `KernelSpec.to_dict` (`core/likelihood.py`) | one-way | Correct as-is: values are `Parameter`s; *reconstruction* of a kernel from a spec is Phase 2's lowering business (the §12.8 registry names the route for user kernels/priors/bijections, and "serialise it, if their registration round-trips" is that hook's own bar) |
| `Likelihood.to_spec()` (new) | one-way | Landed; see §2. Reconstruction from spec is deliberately out of scope: a spec identifies, the registry (Phase 2) reconstructs |
| `Censoring` (`core/likelihood.py`) | read-side only (`from_extra_coord`) | Acceptable: the declaration lives beside the container; a writer back into `extra_coords` lands if and when a modality needs to *store* censored observations (none does yet — recorded as a named absence, not a gap) |
| `Capabilities.to_dict`, `Failure.to_dict` (`core/dataset.py`) | one-way | Correct: run *records*, not declarations; nothing should ever be rebuilt from them |
| Container/`ModelResult`/training-pair functions (`results/serialisation.py`) | containers and results round-trip, tested; `training_pair_to_dict` one-way | See §4 (the design-horizon (c) requirement). The subclass-invariant re-validation gap stays `results.md` §13.12's named limitation with the `validate()` classmethod as extension point |
| Provenance fingerprints and hashes (`results/provenance.py`) | one-way by design | Sound; sensitivity is tested (kernel, solver config, buffers, ties, fixed values). Spec/data hashes are the backend-spanning promise; the problem hash stays backend-variant (W1.10's narrowed claim, ratified) |
| Run emission (`results/emission.py`) | file-level round trip, tested | Sound. `from_netcdf` returns data, never a problem — by design |
| `Instrument`/`Transformation`/`Model` | **no `to_spec` at all** | **Deliberate, now recorded here rather than implicit**: chains and models are *code*, and ampere will not pretend a class name plus buffer hashes reconstructs one. Identity for caching is `model_fingerprint`/`_describe_instrument` (class, module, parameter spec, buffer hashes) — sufficient for spec-hash invalidation of trained artefacts, which is the only consumer design horizon (c) names. A declarative instrument-chain spec would only be worth its maintenance the day a *declarative chain builder* exists; nothing in the plan wants one |
| `Dataset`/`DatasetCollection`/`Plate`/`Tie`/`Binding` | no `to_*` | Same disposition: identity via `problem_fingerprint` (which records ties as `{name, sites}` and plates via the per-parameter tag and shape); structural round-trips have no consumer and are not promised |

## 4. The training-set requirement (design-horizon (c)), nailed down

`results_schema.md` §17 Q6 asked that serialisable (θ, `ModelResult`)
training sets be pinned before Phase 2. They are, as follows:

- **The pair format is fixed**: `training_pair_to_dict(theta, result, failed)`
  over `model_result_to_dict` — exactly what `simulate()` produces, θ
  attached to the result per `results_schema.md` §17 Q7's ruling. Containers
  serialise by value with kind, axes, units, masks and `extra_coords`
  complete, so a coordinate-conditioned emulator has everything the
  functional-data stance requires.
- **The invalidation key is fixed**: the spec hash (plus buffer
  fingerprints inside the problem fingerprint) — §7's cache-key trap is
  covered for emulators exactly as for SBI posteriors.
- **What Phase 2 owes** (carried into the Phase 2 work items): the
  training-set *writer* itself (`results.md` §11 layer 2's table is the
  format; a `training_pair_from_dict` completes the round trip when the
  first consumer — the SBI layer or an emulator trainer — lands), and the
  batching/append story for large budgets. Two known losses to fix in that
  writer, cheap while nothing is archived: θ array values go through
  `.tolist()` (dtype lost), and `Simulation.observations`/`Failure` detail
  are not carried (only the `failed` flag).

  **Discharged at W2.8** (decision-log row in `DEVELOPMENT_PLAN.md` §2).
  `ampere.results.training` writes, appends to and reads §11's format;
  `training_pair_from_dict` completes the in-memory round trip; θ keeps its
  dtype, which bumped `CONTAINER_SCHEMA_VERSION` to **2** under §5's rule
  below; and `Simulation.observations` and the whole `Failure` record travel,
  for which `results.md` §11's table gained two rows. Append is a real
  operation on an existing file — the `sample` dimension grows and the spec
  hash is checked, so §7's stale-artefact trap refuses rather than mixing two
  simulators in one file — implemented as read-concatenate-rewrite, which is
  the one limitation the item created and recorded in place of the two it
  closed.

## 5. Versioning policy, restated as the contract

Three version stamps exist and each guards its own layer:
`ParameterSet.to_spec()`'s `"version"` (layer 1),
`CONTAINER_SCHEMA_VERSION` (layer 3's plain-data form), and
`PROVENANCE_SCHEMA_VERSION` (layer 2 — bumped to 2 by this review). The rule
going forward: any change to what a mapping *means* bumps its layer's
version in the same commit, and readers refuse versions they do not know —
already the behaviour of `from_spec`/`container_from_dict`, and the
convention Phase 2's writers inherit.
