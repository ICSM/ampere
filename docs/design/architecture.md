# Ampere v2 — Architecture Spec (W1.2)

Status: **frozen at `spec-v1.0`** (the tag created at the W1.13 merge,
2026-09) — merged 2026-09-01 (Fable-reviewed, reconciled against W1.1);
Peter's review pass remains the formal accept gate. §10 records
the reconciliation outcome; the curated-translation default flagged during
review has been ruled opt-in/never-silent and recorded in
`DEVELOPMENT_PLAN.md` §2's decision table (§1, §9).

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
specific, recognised model classes and parameterisations.

**Decided (2026-09-01, plan §2 decision log): the substitution is opt-in,
never silent.** A curated native implementation is not guaranteed to be
numerically identical to the astropy original, and a fit whose model changed
implementation without the user's knowledge is exactly the silent-downgrade
class this section forbids in the other direction. The default adapter path
always wraps the user's actual astropy model as a black box; promotion to a
native equivalent happens only through an explicit, backend-scoped request —
sketched as a `from_astropy()` constructor on the backend subpackage, so the
call site itself names the intent, and raising if the model (or any
component of a compound model) is not fully in the curated table rather than
silently falling back to black-box. The exact API is fixed with the adapter
contract (Phase 4).

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
│                     # numpy/scipy/astropy/stdlib — the required base
│                     # dependencies (astropy is one: §4.1 of the plan puts
│                     # units on parameters, and astropy_compat.py below is
│                     # core). Everything else depends on core; core depends
│                     # on nothing backend-specific or optional.
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
├── diagnostics/     # Misspecification diagnostics (§4.8) — the peer
│                     # namespace contracts/diagnostics.md §6 fixes; lands
│                     # Phase 2 (with its decision-log entry then), extra
│                     # `diagnostics` for the JAX-carrying RHMF family.
├── data/, models/, infer/, utils/   # LEGACY — frozen; see plan §2.
└── (docs/design/, tests/, examples/ — not package code)
```

| Extra | Adds | Gates |
|---|---|---|
| *(none)* | reference backend, astropy adapter, emcee/dynesty, celerite2-numpy | the base install; must always work |
| `zeus` | zeus-mcmc | `inference`'s zeus driver |
| `arviz` | arviz | legacy postprocessing today; folded into the base install **with Phase 2's engine drivers**, together with a netCDF engine (ruled 2026-09-03, `results.md` §15 R1) — §4.6 of the plan makes ArviZ's format the single results format (`xarray.DataTree` since ArviZ 1.0 retired the `InferenceData` class; same format, corrected wording per R3) |
| `sbi` | torch, sbi | `inference`'s SBI layer (also unlocks `backends/torch` incidentally, but does not itself require native-model authoring) |
| `extinction` | dust_extinction | legacy `extinctionModels.F99Extinction` (W0.5) |
| `torch` | torch, (GP solver library — deferred choice, plan §6) | `backends/torch` |
| `jax` | jax, (GP solver library — deferred choice, plan §6) | `backends/jax` |
| `dev` | pytest, ruff, pyrefly, sphinx, … | contributor tooling |
| `all` | everything above | — |

## 4. Extras and lazy-import policy

Rules, binding on all new code:

1. **`ampere.core` never imports an optional dependency**, lazily or
   otherwise. It is pure numpy/scipy/astropy/stdlib — the required base
   dependencies and nothing else. This is what keeps the base install
   light and keeps `backends/torch` and `backends/jax` interchangeable
   consumers of the same contracts. (Amended 2026-09-01 during the W1.3
   review: the original wording omitted astropy, which the plan requires
   for units on parameters and which is a required dependency of the base
   install — `ampere/core/parameter.py` importing `astropy.units` at
   module level is in-policy.)
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
- **jax**: the policy is **x64 on, always, for this backend** — and the
  mechanism was ruled 2026-09-01 (decision log, amending this bullet's
  original set-on-first-import sketch): ampere **never** flips
  `jax_enable_x64` as an import side effect, for the same reason this
  section forbids `torch.set_default_dtype` — a library mutating global
  interpreter state on import changes every other jax user in the
  process. Instead `ampere.backends.jax` ships an explicit, idempotent
  `configure_x64()` (deliberately not named `enable_x64`, which jax now
  uses for a context manager with different semantics), and construction
  of any jax-backed likelihood, GP or model **raises** when the flag is
  off, naming the three remedies: `JAX_ENABLE_X64=1` in the environment
  (the only route guaranteed to precede array creation), the explicit
  call, or a per-run, provenance-recorded `float32` opt-out. numpyro
  needs no separate switch — its `enable_x64` is a verified thin wrapper
  over the same jax flag. Full analysis in `docs/design/lowering.md`
  §10.2.
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

- **Curated-translation default (§1)**: *resolved 2026-09-01* — ruled
  opt-in/never silent; recorded in `DEVELOPMENT_PLAN.md` §2's decision
  table, §4.7 amended to match. No longer open.
- **jax x64 activation call site**: *resolved 2026-09-01* — ruled
  guard-and-raise, never set-on-import; recorded in the plan's §2 decision
  table, this document's §5 amended to match, analysis in
  `docs/design/lowering.md` §10.2. No longer open.
- **GP solver library per backend** — already an open deferred choice in
  `DEVELOPMENT_PLAN.md` §6; unaffected by this document.
- **`OptionalDependencyError` exact shape** (fields, message format) —
  pinned by W1.3 in `ampere/core/exceptions.py` as this section asked;
  **ratified in place at the freeze** (ruled 2026-09-03 with
  `lowering.md` §12.4's "the relevant exceptions" disposition; W1.13
  landed `LoweringError` beside it and `ResultsError` had already moved
  to core). No longer open.

## 10. Reconciliation with W1.1 (done 2026-09-01)

This document was drafted in parallel with W1.1 (prior-art memo), per the
project's orchestration policy — not after it, despite the plan's stated
dependency. The memo has since landed (`docs/design/prior_art.md`), and the
four questions this section originally posed were checked against it in the
2026-09-01 Fable review:

- **bilby (memo §1)**: lessons B1/B3 corroborate the framing here rather
  than changing it — bilby's own core moved to explicit-argument
  `log_likelihood(parameters)` (re-verified against current `bilby-dev`
  source during the review), and its sampler-parametrised `Result` grab-bag
  is precisely the failure mode the single-InferenceData decision avoids.
  No change to §1's capability-flag / `FittingProblem` framing.
- **gammapy (memo §2)**: G1 (joint log-likelihood as a sum over member
  datasets; stacking is an opt-in change of objective, never a silent
  optimisation) and G2 (identity-based tying broke a downstream consumer;
  tie by explicit declaration instead) bear on W1.7 and W1.3 respectively
  — W1.3 implements name/declaration-based tying. The ladder and
  "execution venue" language here needed no change.
- **3ML (memo §3)**: 3M1's verdict is to keep the factored
  Transformation/NoiseModel split, with an awkward-instrument stress test
  (physically coupled calibration and noise) delegated to W1.11's modality
  sketches before the freeze. §1's rung description and the §4.3/§4.4
  boundary stand.
- **Starfish (memo §4)**: S1 corroborates the Matérn-3/2 default; S2's
  windowed-sparse truncation is a third solver-strategy family for W1.6 to
  name (or explicitly decline); S3's trans-dimensional local kernels are a
  scope boundary for W1.6 to state. Nothing touches the dtype/precision
  policy (§5) or the reference-backend scope (§2).

Peter's review pass remains this item's outstanding acceptance criterion.
