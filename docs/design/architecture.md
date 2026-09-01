# Ampere v2 — Architecture Spec (W1.2)

Status: **DRAFT — pending reconciliation with W1.1 (prior-art memo,
in progress) and Peter's review.** This is not yet frozen; §10 lists what
might move once the prior-art memo lands.

Relationship to other documents: `DEVELOPMENT_PLAN.md` §2–3 record the
architectural *decisions*; this document expands them into the connective
tissue — the cross-cutting policies every per-contract spec (`contracts/`,
W1.3–W1.9) must follow, so those specs can each stay focused on their one
contract rather than re-deriving shared ground rules. Where this document
and `DEVELOPMENT_PLAN.md` disagree, `DEVELOPMENT_PLAN.md` wins and this
document is wrong — file it as a decision-log correction.

---

## 1. The capability ladder

Ampere v2 is not four peer backends; it is one backend-neutral core plus an
escalating capability ladder, because the axis that actually matters is
**differentiability**, not array library:

| Rung | What runs here | Inference available |
|---|---|---|
| 0 — black-box | External RT codes (Dusty, Hyperion), legacy numpy models, any user callable ampere cannot see inside | Gradient-free samplers (emcee, zeus, dynesty) + SBI |
| 1 — reference | `backends/reference` (pure numpy/scipy); astropy-adapted models execute here by default | Same as rung 0, plus the O(N) GP contract as a correctness anchor |
| 2 — native | `backends/torch`, `backends/jax` | Everything above, plus NUTS/HMC, VI (pyro/numpyro), gradient-based optimisation, GPU batching |

This is a **ladder, not a matrix**: differentiability is monotonic
capability. Anything rung 0 can do, a rung-2 model can also do (by running
its gradient-free/SBI path); the reverse never holds. A model does not
choose a rung — its capability flags (`differentiable`, `batchable`,
`device`; §4.5 of the plan) are a fact about how it was written, and the
`FittingProblem` surface exposes them so inference engines can self-select.
**An engine requesting a capability the problem doesn't have must fail
loudly and specifically** ("this model is not differentiable; NUTS is
unavailable — see capability flags") — never silently downgrade, and never
fail deep inside a tensor op with an opaque autodiff error.

Astropy-adapted models (§4.7 of the plan) are rung 1 by default: the
adapter's parameter translation and compound-model support make them easy
to use, but the underlying `astropy.modeling` evaluation is not
differentiable. The curated-translation escape hatch (mapping common
analytic astropy models — blackbody, power laws, polynomials — to native
torch/jax equivalents) is precisely a rung-1-to-rung-2 promotion for
specific, recognised models; it is opt-in and the adapter must never
silently substitute a native model for a user's actual astropy definition
without that model being byte-for-byte in the curated table.

## 2. Reference backend: the trade-off, written out

`DEVELOPMENT_PLAN.md` §3 records the decision; here is the analysis behind
it, since W1.2's brief asks for it explicitly.

**Option A — promote legacy to a full backend.** Legacy `ampere.data`/
`ampere.models`/`ampere.infer` already runs, has demonstrated science
output, and covers the v1 slice (spectra + photometry). Promoting it would
mean zero new numpy-side implementation work.

Cost: legacy has zero test coverage before W0.4 (mitigated, not solved —
the characterisation suite protects known flows, not the full surface),
carries in-place mutation as its idiom (`self.covMat`, `self.bestPars` set
as side effects rather than returned), and until this year had two
coexisting model-output conventions (the exact bug class §4.2's schema
exists to end). Making it conform to the frozen §4 contracts means
touching its internals — which is precisely what "frozen in place" (a
deliberate risk-reduction decision, not an oversight) forbids. Option A is
self-contradictory: contract conformance and frozen-in-place cannot both
hold for the same code.

**Option B — a small, new, pure-numpy reference backend.** Cost: a third
implementation of the §4 contracts to write and maintain, alongside
torch and jax. Scoped down deliberately — no gradients, no GPU, no
performance target — so this is materially cheaper than a third full
backend; it is closer in size to a large contract-conformance test fixture
than to a production backend.

Benefit, three-fold, made concrete:

1. **Conformance oracle.** The suite (§4.6 of the plan; W1.10) computes
   reference values — prior round-trips, `log_prob` on analytic cases,
   `DenseGP` outputs at small N — and asserts torch/jax agree with the
   reference to a **documented numerical tolerance**: tight (`rtol` on the
   order of 1e-6 to 1e-8) for closed-form/exact-arithmetic cases, looser
   and backend-pair-specific for GP solves that take genuinely different
   linear-algebra paths (dense Cholesky vs quasiseparable recursion vs
   different autodiff accumulation order). The exact tolerance table is
   W1.10's job; this document only commits to "documented and specific
   per comparison", never a single blanket tolerance.
2. **A useful base install.** `pip install ampere` (no extras) is a
   complete, scalable, numpy-only fitting environment: reference backend +
   astropy adapter + emcee/dynesty + the O(N) flexible likelihood. For
   that last clause to be true rather than aspirational, **celerite2's
   numpy interface is a core dependency of `backends/reference`, not gated
   behind an extra** — the reference backend's entire point is being a
   complete, scalable environment with no heavy deps, and a base install
   that still has the O(N³) problem would defeat it. (zeus stays an
   extra — it is a genuine alternative sampler, not part of the flagship
   scaling story.)
3. **Execution venue for adapted models.** Astropy-adapted and
   legacy-adapted (black-box) models are not evaluated *by* the reference
   backend in the sense of differentiating through them — they're opaque —
   but the reference backend supplies the concrete numpy-side
   Instrument/Likelihood/Dataset/inference-contract machinery that wraps
   their raw output. This is what "rung 1" in §1 means operationally: the
   scaffolding around a black-box or astropy model at rung 1 is the
   reference backend's implementation of §4.

Option B wins. It is the standing decision.

## 3. Namespace layout and extras mapping

```
ampere/
├── core/            # Backend-neutral contracts. STRICT: no torch/jax
│                     # imports, ever, at any point, even lazily. Only
│                     # numpy/scipy/typing/stdlib. Everything else depends
│                     # on core; core depends on nothing backend-specific.
│   ├── parameter.py      # §4.1 — Parameter, ParameterSet, priors, tying
│   ├── results_schema.py # §4.2 — ModelResult, named channels, containers
│   ├── transform.py      # §4.3 — Transformation, Instrument, negotiation
│   ├── likelihood.py     # §4.4 — Likelihood, NoiseModel strategy
│   ├── dataset.py        # §4.5/4.7 seed — Dataset, DatasetCollection
│   ├── exceptions.py     # OptionalDependencyError and friends (§4 below)
│   └── astropy_compat.py # §4.7 — astropy.modeling adapter
├── backends/
│   ├── reference/   # numpy/scipy. Core dependency: celerite2 (numpy).
│   ├── torch/       # requires extra `torch`. May import torch at module
│   │                # top level WITHIN this subpackage only (see §4).
│   └── jax/         # requires extra `jax`. Same top-level-import carve-out.
├── inference/       # Backend-agnostic engines consuming §4.5's contract:
│                     # emcee/dynesty/zeus drivers, the SBI layer, optimisers.
├── results/         # ArviZ InferenceData emission + all plotting (§4.6).
├── data/, models/, infer/, utils/   # LEGACY — frozen; see plan §2.
└── (docs/design/, tests/, examples/ — not package code)
```

| Extra | Adds | Gates |
|---|---|---|
| *(none)* | reference backend, astropy adapter, emcee/dynesty, celerite2-numpy | the base install; must always work |
| `zeus` | zeus-mcmc | `inference`'s zeus driver |
| `sbi` | torch, sbi | `inference`'s SBI layer (also unlocks `backends/torch` incidentally, but does not itself require native-model authoring) |
| `extinction` | dust_extinction | legacy `extinctionModels.F99Extinction` (W0.5) |
| `torch` | torch, (GP solver library — deferred choice, plan §6) | `backends/torch` |
| `jax` | jax, (GP solver library — deferred choice, plan §6) | `backends/jax` |
| `dev` | pytest, ruff, pyrefly, sphinx, … | contributor tooling |
| `all` | everything above | — |

## 4. Extras and lazy-import policy

Rules, binding on all new code:

1. **`ampere.core` never imports a heavy or optional dependency**, lazily
   or otherwise. It is pure numpy/scipy/stdlib. This is what keeps the
   base install light and keeps `backends/torch` and `backends/jax`
   interchangeable consumers of the same contracts.
2. **Outside `ampere.core`, an optional dependency is imported at module
   top level only within the subpackage that exists *because of* that
   dependency** — e.g. `ampere/backends/torch/__init__.py` may `import
   torch` at the top, because nothing imports that subpackage without
   wanting torch. Any module reachable without opting into that
   subpackage (`ampere/backends/__init__.py` itself, `ampere.inference`,
   `ampere.results`) must import optional dependencies lazily, inside the
   function/method/`__init__` that needs them.
3. **One shared exception type carries the failure**:
   `ampere.core.exceptions.OptionalDependencyError(package, extra,
   context)`, raised on *use*, never on package import. This standardises
   what W0.5 did ad hoc per-module in legacy (each module rolled its own
   message) — new code uses the one helper, so tests and tooling can
   assert on it uniformly.
4. **Enforcement is the minimal-install CI job** (already live, W0.6): a
   bare `pip install -e .` followed by `import ampere` and the import
   sweep (`tests/test_imports.py`). New namespaces inherit this guard as
   they land; a PR adding an eager heavy import outside its own backend
   subpackage fails CI, not review-by-eyeball.

## 5. Dtype, device, and precision policy

(Expands the float-precision trap in `DEVELOPMENT_PLAN.md` §7 into an
actual policy; exact per-backend mechanics belong to `lowering.md`, W1.9 —
this section commits to the *policy*, not the incantations.)

- **Default is float64 for all likelihood and GP linear algebra, on every
  backend, always.** GP Cholesky/quasiseparable solves in float32 fail in
  ways that read as science bugs, not numerical ones — this is not a
  performance knob to leave at each backend's default.
- **torch**: dtype is threaded explicitly through construction (models,
  buffers, parameters all declare or inherit an explicit dtype), never
  relied upon via `torch.set_default_dtype` — that call is global mutable
  process state and unsafe to set implicitly from library code that
  another user's process also imports.
- **jax**: `jax.config.update("jax_enable_x64", True)` is process-global
  and must be set at the *earliest* possible point when `ampere.backends.
  jax` is first imported, with a loud one-time log message — jax warns
  (and silently truncates precision) if arrays are created before the
  flag is set. The exact call site and interaction with a host
  application that also uses jax is a W1.9 gate check (flagged there
  already); this document only fixes the policy: **x64 on, always, for
  this backend.**
- **Device is CPU by default, GPU strictly opt-in and explicit**
  (`device="cuda"` or backend-equivalent) — never auto-detected. This
  keeps behaviour reproducible and guarantees CI (CPU-only runners)
  exercises the real default path, not a path nobody ships.
- **Buffers move with parameters under device/dtype changes** (`.to(...)`
  / `device_put` equivalents apply to both) — this is the operational
  reason the buffer/parameter distinction (§6) is a contract requirement
  and not just a JIT-cache-key nicety.

## 6. Parameters vs buffers — the model contract

Expands `DEVELOPMENT_PLAN.md` §4.1's buffer note into a usable test.

**The distinguishing question**: *would you ever put a prior on it,
sample it, or take a gradient with respect to it?* If no, it is a buffer.
Wavelength grids, opacity tables, filter curves: buffers. Temperature,
scale factor, calibration offset: parameters. A fixed cosmological
constant the user has not chosen to fit: a buffer *by default*, but the
contract must make it trivial to re-declare as a parameter — promotion
between buffer and parameter should be a **configuration-level choice**
(what the user declares when building the model), never a code-level
fork requiring a different model implementation for the fit-it vs
fix-it cases. This is a real test of whether the parameter contract
(W1.3) is well designed: if promoting a buffer to a parameter requires
touching the model's `__call__`, the contract has failed.

**Decided**: buffers are declared **explicitly**, alongside parameters,
not inferred from "any array attribute that isn't a parameter". The
alternative (implicit: unregistered arrays default to buffer status) saves
boilerplate but fails silently — a typo'd or forgotten parameter
registration quietly becomes an untracked buffer, which is exactly the
kind of silent-drift bug class §4.2's ModelResult schema exists to end
elsewhere in the contracts; the parameter/buffer boundary deserves the
same discipline. Explicit declaration is also the friendlier shape for
pyrefly, since a buffer's presence is then a checkable structural fact
about the class rather than something only known at runtime. W1.3
implements this as the binding API.

Per-backend lowering mechanics (torch `register_buffer`; jax via
partition filters, never `equinox` static fields — see the trap in plan
§7) are W1.9's table, not repeated here.

## 7. The functional-data stance

Expands `DEVELOPMENT_PLAN.md` §4.2's "containers are coordinate-indexed
function samples" into a concrete illustration (the container class
hierarchy itself is W1.4's contract, not this document's job).

A `Spectrum` channel is **not** two positionally-aligned arrays
(`wavelength`, `flux`) with implicit ordering assumptions baked into
whatever consumes them — that is precisely the legacy idiom, and legacy
code (`spectres`-based resampling, the dense covariance-matrix
construction) implicitly assumes sorted, deduplicated wavelength arrays
without ever stating so. A container is a **coordinate + value (+
uncertainty, + mask)** structure that either validates its ordering
assumptions explicitly at construction, or explicitly documents that it
tolerates arbitrary coordinate order — never an unstated assumption
inherited by whoever reads the array next.

A container may **advertise** structure (`regular: bool`, or a richer
grid-structure descriptor) so implementations can take fast paths — FFT
convolution on a regular grid, Toeplitz structure in a covariance
solve — but nothing in the contracts may *require* regularity. This
connects directly to requirements negotiation (§4.3 of the plan): an
instrument's published requirements and a model's advertised output
structure are exactly the two things that negotiation reconciles.

## 8. What this document does not own

| Detail | Owned by |
|---|---|
| Exact `Parameter`/`ParameterSet` API, tying/plate syntax | W1.3, `contracts/parameters.md` |
| Exact `ModelResult`/container class hierarchy, channel API | W1.4, `contracts/results_schema.md` |
| Exact `Transformation`/`Instrument`/negotiation API | W1.5, `contracts/transformations.md` |
| Exact `Likelihood`/`NoiseModel` API, kernel registry | W1.6, `contracts/likelihoods.md` |
| Exact `Dataset`/`DatasetCollection`/`FittingProblem` API | W1.7, `contracts/inference.md` |
| Exact `InferenceData` emission helpers, plotting API | W1.8, `contracts/results.md` |
| Distribution mapping table, per-backend lowering mechanics, RNG lowering | W1.9, `lowering.md` |
| Per-modality worked compositions | W1.11, `modalities/` |
| Diagnostics API | W1.12, `contracts/diagnostics.md` |

## 9. Open items carried from this document

- **jax x64 activation call site** in a multi-library host process — W1.9
  gate check.
- **GP solver library per backend** — already an open deferred choice in
  `DEVELOPMENT_PLAN.md` §6; unaffected by this document.
- **`OptionalDependencyError` exact shape** (fields, message format) —
  small, but should be pinned once, in `core/exceptions.py`'s own
  docstring, so it isn't reinvented per contract spec.

## 10. Status and what W1.1 might change

This document was drafted in parallel with W1.1 (prior-art memo), per the
project's orchestration policy — not after it, despite the plan's stated
dependency — because the capability-ladder framing above does not depend
on prior-art specifics to state, and starting the specs in parallel saves
a full round trip. It is **not final**. Once W1.1 lands, re-check
specifically:

- Whether bilby's `Likelihood`/`PriorDict`/`Result` separation suggests a
  different shape for the capability-flag / `FittingProblem` framing in
  §1.
- Whether gammapy's `Datasets` joint-fitting pattern changes how §3's
  "execution venue" language should describe multi-dataset composition
  (this bears more on W1.7 than this document, but the ladder framing
  should stay consistent with it).
- Whether 3ML's per-instrument plugin-likelihood pattern argues for
  likelihoods being instrument-owned rather than composed separately —
  this would touch §1's rung description and the boundary between §4.3
  and §4.4 of the plan.
- Whether Starfish's implementation experience surfaces a kernel or
  cost pitfall the dtype/precision policy (§5) or the reference-backend
  scope (§2) should account for.

Peter's review pass is the other half of this item's acceptance
criterion and is independent of W1.1's content.
