# Ampere v2 Development Plan — DRAFT for refinement

Status: **plan settled; work-item breakdown live** (2026-09-01). All
architectural proposals are confirmed; remaining open items are
implementation-level choices deferred to their natural phase (§6). This
document is the source of truth for the redevelopment of ampere: decisions
taken, target architecture, and phased roadmap. The agent-sized work items
for Phases 0–1 live in `WORK_ITEMS.md`.

---

## 1. Motivation

Ampere's core scientific idea — a flexible, GP-based likelihood that provides
robustness to model misspecification — works well but does not scale. The
current implementation builds a dense covariance matrix with a hardcoded
squared-exponential kernel and evaluates it with dense `slogdet`/solves
(`ampere/data/spectrum.py:417-431`), i.e. O(N²) memory and O(N³) time. At the
same time, the codebase has accumulated structural debt (positional parameter
slicing, two coexisting model-output conventions, a 1,000-line
plotting/postprocessing monolith, no tests, no CI) that blocks extension to
new data types, gradient-based inference, and GPU execution.

Rather than incremental patching, this plan pursues a redesign around a
backend-neutral core plus modern computational backends, targeting:

1. **Scalable GPs** — exact O(N) methods for 1D data (celerite-class
   quasiseparable kernels / state-space GPs), with approximate methods
   (SVGP, SKI, Vecchia) as swappable strategies for 2D+ data later.
2. **Extensibility** — a composition system for data, instruments, and
   likelihoods that supports spectra, photometry, IFU cubes and images,
   interferometric visibilities and closure phases, astrometric time series,
   and hierarchical/population inference.
3. **Modern inference** — gradient-based inference (NUTS, VI) where models
   are differentiable, simulation-based inference everywhere, on GPU where
   available.

## 2. Decisions taken (2026-09-01)

| Topic | Decision |
|---|---|
| Target data scale | All of: long 1D spectra (10³–10⁶ pts), many joint datasets, IFU cubes/images, population/hierarchical, interferometric visibilities & phases, astrometric time series. Interfaces must be designed for all; implementations are staged (see §5). |
| Where the O(N) GP lands | **New architecture only.** The legacy code keeps its dense implementation; the paper revision proceeds on legacy at current scale. (A cheap escape hatch exists — celerite2 wired directly into the paper script — if referee/timeline pressure changes.) |
| Modern backends | **torch and jax in lockstep**, implemented in parallel against a frozen interface spec by separate agent tracks. The shared conformance test suite (§4.6) is the synchronisation mechanism. |
| Legacy compatibility | **Frozen in place.** `ampere.data`, `ampere.models`, `ampere.infer` keep working as-is where they are; only critical bugfixes land there. New code grows in new namespaces alongside. |
| Astropy's role | **Interop adapter in the core, not a peer backend.** The objective is: a user defines their model with `astropy.modeling` and ampere handles it. The adapter translates astropy parameter metadata automatically (§4.7). Astropy-defined models are black-box to torch/jax (gradient-free inference + SBI, not NUTS/VI), with a curated-translation escape hatch for common analytic models. |
| Flexible-likelihood kernel | **Matérn-class (Matérn-3/2 / sums of SHO terms) becomes the default throughout**, replacing the legacy RBF — it represents structured residuals better and is exactly O(N) via quasiseparable solvers. Validated against the misspecification study in milestone M2. |
| Packaging | **Single distribution with extras** (`ampere[torch]`, `ampere[jax]`, `ampere[sbi]`, …), lazy imports; split distributions are not pursued. |
| Phase-4 modality proof | **Interferometric visibilities** — complex-valued data and Fourier sampling stress the result schema and transformation chains hardest; other modalities follow the template it establishes. |
| CI/CD | **Full CI/CD is in scope**, as a cross-cutting workstream (§5) growing from lint+tests (Phase 0) through conformance/backend matrices and benchmark tracking (Phases 1–2) to automated PyPI releases and versioned docs (Phase 6). Includes pyrefly type checking scoped to the new namespaces. |
| Reference backend | **Confirmed**: pure-numpy `backends/reference` as conformance oracle, numpy-only base install, and execution venue for astropy-adapted models (§3). |
| Result schema | **Confirmed**: named channels + requirements negotiation; containers are coordinate-indexed function samples with first-class masks; no regular-grid assumptions anywhere in the contracts (§4.2–4.3). |
| Non-Gaussian robustness | **Accepted**: flexible-GP robustness for non-Gaussian likelihood families (latent-GP path) is promised on the torch/jax rungs only (§4.4). |
| Buffers | The model contract distinguishes parameters from **buffers** (constant arrays); lowering rules per backend in §4.1 (torch `register_buffer`; jax via partition filters, never equinox static fields). |
| Diagnostics | **Misspecification diagnostics are in scope** (§4.8): RHMF-style pre-fit screening (Hilder et al. 2026), post-fit residual tests, GP-based localisation. |

## 3. Architecture: a core and a capability ladder, not four peer backends

The four-backend framing (legacy / astropy / torch / jax) is restated as **one
backend-neutral core plus an escalating capability ladder**, because the real
dividing line is not the array library — it is whether a model is
differentiable:

- **Black-box models** (external RT codes: Dusty, Hyperion, …; legacy numpy
  models) can never be differentiated. They get gradient-free inference
  (emcee, zeus, dynesty) and **SBI** — and SBI only needs `(θ, x)` samples, so
  it is backend-agnostic with respect to the model already.
- **Native models** (written in torch or jax) additionally unlock NUTS/HMC and
  VI (pyro/numpyro), gradient-based optimisation, and batched GPU evaluation.

Consequences:

- Gradient-free samplers and SBI live **above** the backends and consume any
  of them through the inference contracts (§4.5). They are written once.
- "A backend" is precisely: an array library for writing models and
  transformations, a set of GP solver implementations, and the extra inference
  engines its differentiability unlocks. Everything else is shared.
- **Astropy is an interop adapter in the core, not a peer backend**
  (decided): a wrapper turning any `astropy.modeling` model into an ampere
  model — inheriting astropy's compound-model algebra, which retires the
  unfinished `eval`-based `CompositeModel` — plus units-awareness in the
  data layer. Contract in §4.7.

### Decided (2026-09-01): a numpy reference backend

Rather than promoting the legacy code to a true backend (it is untested,
carries side-effect state and two output conventions, and refactoring it
contradicts "frozen in place"), add a **small, clean, pure-numpy/scipy
reference backend** (`backends/reference`) implementing the §4 contracts
with correctness as its only goal. It earns its keep three ways:

1. **Conformance oracle** — the suite tests torch and jax against it (and
   against analytic cases), giving lockstep development a ground truth.
2. **A useful base install** — `pip install ampere` (no extras) becomes a
   complete numpy-only fitting environment: reference backend + astropy
   adapter + emcee/dynesty + O(N) GP via celerite2's numpy interface. This
   is, in effect, the "streamlined, extensible, astropy-compatible" backend
   2 of the original four-backend vision, recovered at low cost.
3. **Execution venue** for astropy-adapted and legacy-adapted (black-box)
   models within the new architecture.

Cost: a third contract implementation to maintain — kept small by having no
performance goals, no gradients, and no GPU support. Legacy remains frozen
with only a thin black-box adapter.

### Target namespace layout

```
ampere/
├── core/            # NEW: backend-neutral contracts (§4). No torch/jax imports.
│   ├── parameter.py     # named params, priors, transforms, tying, hierarchy
│   ├── results_schema.py# typed ModelResult containers per observable kind
│   ├── transform.py     # Transformation & Instrument chain ABCs
│   ├── likelihood.py    # Likelihood + NoiseModel strategy interfaces
│   ├── dataset.py       # Dataset, DatasetCollection containers
│   └── astropy_compat.py# astropy.modeling adapter, units interop (D1)
├── backends/
│   ├── reference/   # NEW: minimal numpy oracle + numpy-only base install
│   ├── torch/       # NEW: torch implementations of core contracts
│   └── jax/         # NEW: jax implementations of core contracts
├── inference/       # NEW: backend-agnostic engines (emcee/dynesty/zeus drivers,
│                    #      SBI layer, optimisers) consuming §4.5 contracts
├── results/         # NEW: ArviZ InferenceData emission + all plotting, written once
├── data/, models/, infer/   # LEGACY: frozen in place, critical fixes only
└── utils/
```

Packaging: single distribution with extras (`ampere[torch]`, `ampere[jax]`,
`ampere[sbi]`, …) and lazy imports; torch and jax are never both required.
CI tests them in separate matrix jobs. (Open decision D3 on eventual split
into separate distributions — deferred until beta.)

## 4. Core contracts (Phase 1 deliverable — "the constitution")

These interfaces are frozen, documented, and covered by a conformance suite
**before** backend implementation starts. Prior art to study while writing
them (≤ half a day each): bilby (likelihood/prior/sampler decoupling),
gammapy (Datasets container), 3ML (per-instrument plugin likelihoods),
Starfish (misspecification GPs for spectra).

### 4.1 Parameters and priors
Named, first-class parameters with: prior distribution (declared neutrally,
scipy-style; **lowered** to `torch.distributions` / numpyro equivalents),
optional transform to unconstrained space, units, shape, and — critically for
the target scope — **tying/sharing across models and datasets** and
hierarchical grouping (plate-like, so population models lower cleanly to
pyro/numpyro). Seed: `ampere/core/parameter.py` on
`copilot/explore-implementation-plan`. Positional `theta`-slicing dies here.

**Parameters vs buffers.** Not every numeric attribute is a parameter:
models carry large constant data (wavelength grids, opacity tables, filter
curves). The model contract distinguishes **buffers** — arrays that move
with device/dtype changes and are traced through computations, but are
never sampled, optimised, or differentiated — from parameters. Lowering:
torch has this natively (`register_buffer`); equinox does not, so on the
jax side buffers lower to ordinary array fields excluded from the trainable
partition (`paramax.NonTrainable` or an `eqx.partition` filter spec) —
**never** equinox static fields, which hash array contents into the JIT
cache key (see §7).

### 4.2 ModelResult schema
Typed containers for what a model produces — `Spectrum`, `PhotometricPoints`,
`Image`, `Cube`, `TimeSeries`, `VisibilitySet` — with units and metadata.
This permanently kills the `self.model.modelFlux`-attribute vs
`result['spectrum']` split that currently causes recurring bugs.

**Decided (2026-09-01) — named channels, not just kinds.** A
`ModelResult` is a mapping of *named channels* to typed containers (e.g.
`{"sed_lowres": Spectrum, "co_windows": Spectrum, "image_850um": Image}`),
because one model may legitimately produce several results of the same kind
at different resolutions or in different regions. Instruments bind to a
channel **by name**, with kind checking at composition time; mismatches
fail loudly. Channels may carry a fidelity tag (see Design horizon). The
simple case (one channel) must stay trivial: a model returning a single
container gets a default channel name automatically.

**Containers are coordinate-indexed function samples.** Every observable
is an underlying function — flux(λ), V(u,v), position(t) — and every
container holds explicit coordinates + values (+ uncertainties) sampled
from it, with **no regular-grid assumption anywhere in the contracts**:
spectra, time series, and (u,v) coverage are all irregular in general.
Regularity is a property a container may *advertise* so implementations
can take fast paths (FFT convolution, Toeplitz structure); it is never
required. Quasiseparable GP solvers handle irregular 1D coordinates
natively. This functional view feeds SBI (coordinate-aware embeddings,
Phase 3) and emulation (coordinate-conditioned emulators, Design horizon).

**Missing data are first-class.** Containers carry masks that propagate
through transformations and likelihoods. Masking (data excluded, zero
information) is distinct from censoring (upper/lower limits, which carry
information — issue #11); both are part of the likelihood contract.

### 4.3 Instruments as transformation chains
Promoted from the `jax` branch sketch: an `Instrument` is a composition of
`Transformation`s mapping model output to predicted-data space. Standard
library: LSF convolution, spectral resampling, synthetic photometry,
calibration/scale factors (absorbing the current `scaleFac` nuisance),
epoch sampling (astrometry), Fourier sampling at (u,v) points
(interferometry), and instrument response matrices (X-ray RMF/ARF are
simply matrix-multiply Transformations). Transformations may carry their
own named nuisance parameters (§4.1).

**User extensibility is a first-class requirement**: a user must be able to
define a new Transformation, Instrument, data container, or likelihood
family by subclassing one ABC and implementing one or two methods, with no
changes inside ampere. Acceptance criterion: a "bring your own instrument
and likelihood" tutorial exists, and the test suite includes one such
extension written *out-of-tree* (imported as if third-party) to prove the
interfaces suffice.

**Decided (2026-09-01) — requirements negotiation.** Instruments
may publish their requirements (wavelength coverage, resolution, epochs,
(u,v) coverage); an optional one-off *compilation* step before the hot loop
lets the model configure per-channel evaluation grids from the union of
requirements. This directly serves the low-res-SED + high-res-windows
pattern (§4.2), avoids evaluating expensively everywhere, removes redundant
per-call resampling (issue #12), and creates the natural caching point.
Models are free to ignore requests (compute what they compute); the simple
path must remain: fixed grid, no negotiation.

### 4.4 Likelihoods and noise models
`Likelihood` compares predicted vs observed data given a `NoiseModel`.
The GP solver is a **swappable strategy** behind `NoiseModel`:

- `DenseGP` — current behaviour, any kernel, O(N³); correctness reference.
- `QuasisepGP` — exact O(N) for 1D ordered data (celerite2 numpy/jax,
  `tinygp.solvers.QuasisepSolver`, GPyTorch/celerite2-torch on the torch
  side). **Decided:** the canonical flexible-likelihood kernel becomes
  Matérn-class (Matérn-3/2 / sums of SHO terms), which is exactly
  quasiseparable and represents structured residuals better than the legacy
  hardcoded RBF; new-architecture posteriors will therefore not be
  bit-comparable to legacy ones. Validated against the misspecification
  study in M2.
- Future strategies for 2D+: SVGP/inducing points, SKI/KISS-GP, Vecchia —
  interface slots designed now, implemented in Phase 5.

Beyond Gaussian: likelihood families needed by the target scope — Poisson
(X-ray counts), Student-t and Cauchy (outlier robustness), Rice (polarised
intensity), complex Gaussian (visibilities), wrapped/von Mises (phases),
censored/upper-limit support (issue #11) — are part of the interface design
even where implementation is staged. Families register through a minimal
interface (`log_prob(predicted, observed, noise_params)`), so users can add
their own without touching ampere.

**Structural consequence of non-Gaussian families**: the "GP as covariance
matrix added to a Gaussian" trick only marginalises analytically for
Gaussian likelihoods. For Poisson and friends, misspecification robustness
needs a **latent-GP formulation** (e.g. counts ~ Poisson(rate·exp(f)),
f ~ GP), marginalised numerically — which is exactly what the PPL backends
are for (HMC/VI over the latent function; still O(N) per gradient
evaluation with a quasiseparable GP prior). The Likelihood interface must
therefore declare whether it marginalises its noise process analytically or
introduces latent variables; gradient-free samplers cannot realistically
handle hundreds of latent values, so non-Gaussian + flexible-GP robustness
is effectively a modern-backend capability (Laplace-type approximations are
a possible later fallback for the reference path).

### 4.5 Inference contracts
A fitting problem (model + instruments + datasets) exposes:

- `log_prob(params) -> float` and `log_likelihood` / `log_prior` split;
- `prior_transform(u)` for nested sampling;
- `simulate(params) -> data` for SBI;
- capability flags: `differentiable`, `batchable`, `device`.

Any engine consuming only this contract works with every backend, including
legacy black-box models (via a thin adapter over the frozen legacy classes).

Also part of this contract:
- **Failure signalling**: external simulators crash and return NaNs; the
  contract defines the behaviour (`log_prob` → −inf with a recorded reason;
  `simulate` failures are flagged so SBI can reject-and-record rather than
  train on garbage).
- **RNG policy**: named seed handling that lowers to each backend's model
  (numpy Generators, torch Generators, jax PRNG keys), so runs are
  reproducible across backends.

### 4.6 Results, plotting, and the conformance suite
- **ArviZ `InferenceData` is the single results format.** Every engine emits
  it; corner/trace/posterior-predictive plotting is written once against it
  in `ampere.results`; serialisation is netCDF. This is how the
  `mixins.py` monolith is retired and how sampler plotting parity stops
  regressing.
- **Every run stores per-sample `log_likelihood` and `log_prior`** (and run
  provenance: package versions, seeds, data hashes, spec hash, in the
  InferenceData attrs). Cheap now; it is the enabling requirement for
  population-level importance reweighting later (Design horizon, item b).
- **Conformance suite**: a pytest suite parametrised over backends that any
  implementation of the contracts must pass (round-trip priors, schema
  validation, log_prob agreement vs analytic cases, DenseGP↔QuasisepGP
  agreement on quasiseparable kernels, InferenceData emission). This is the
  contract that keeps two lockstep backends aligned across agent tracks.

### 4.7 Astropy interop contract
`core/astropy_compat.py` wraps any `astropy.modeling` model — including
compound models — as an ampere model:

- **Parameter translation is automatic**: each astropy `Parameter` (name,
  value, `bounds`, `fixed`, `tied`) maps to an ampere parameter (§4.1);
  `fixed` → frozen, `tied` → tying, `bounds` → a default uniform prior.
  Users override priors via a simple dict; nothing else is required, so a
  casual user's astropy model "just works".
- **Units**: astropy models with units and `Quantity` inputs are supported;
  the data layer accepts Quantities throughout.
- **Output kind**: the user declares (or the adapter infers where
  unambiguous) which ModelResult kind (§4.2) the model produces.
- **Capability consequence (advertise clearly in docs)**: astropy-defined
  models are black-box to torch/jax, so they get gradient-free inference and
  SBI — not NUTS/VI. As a later nicety, a curated translation table can map
  common analytic astropy models (blackbody, power laws, polynomials, …) to
  native torch/jax equivalents, silently restoring differentiability for the
  most frequent cases.

### 4.8 Misspecification diagnostics
Tools that tell users **where** a flexible likelihood is needed, rather
than leaving it as an act of faith:

- **Pre-fit, data-driven screening**: robust low-rank factorisation in the
  style of RHMF (Hilder, Hogg, Casey & Rix 2026, arXiv:2607.08081 —
  robust heteroskedastic matrix factorisation with an implicit Student-t
  likelihood, per-feature uncertainties, missing-data support, and
  automatic per-feature/per-object anomaly flags; JAX implementation
  Robusta-HMF), and extensions thereof, applied across collections of
  spectra to flag the features/regions any smooth model will struggle
  with. Its missing-data and heteroskedastic machinery aligns directly
  with the masks-first-class container design (§4.2).
- **Post-fit residual tests**: whiteness/autocorrelation statistics
  (Ljung-Box-style) and posterior-predictive checks on standard-likelihood
  fits, flagging structure that motivates enabling the GP.
- **Localisation from the flexible fit itself**: the posteriors on GP
  amplitude and length-scale, and the conditioned GP mean, already
  localise where the model is deficient — surfaced as standard plots in
  `ampere.results`.

Design spec in Phase 1 (including an adoptability check on Robusta-HMF);
1D implementation lands with Phase 2.

## 5. Phased work breakdown

Scope discipline: interfaces are designed for the full data scope (§2), but
the **v1 vertical slice is spectra + photometry** — today's working use
cases — delivered end-to-end on both modern backends before any new modality
is built. Phases 0 and 1 are decomposed into agent-sized work items in
`WORK_ITEMS.md`; later phases are decomposed when their prerequisites
freeze.

### Phase 0 — Safety net & hygiene (small, do first)
- Commit the pending `mixins.py` MAP-plot fix.
- Characterisation tests: golden-output tests with fixed seeds wrapping the
  `minimal_working_example*.py` scripts (emcee, dynesty, zeus, SBI), so the
  frozen legacy code is protected without refactoring it.
- CI (GitHub Actions): lint (ruff) + pytest; Python floor **≥3.11**;
  `requires-python` updated.
- Fix-or-quarantine the broken legacy modules (#74 import error, #75
  `extinction` dependency, #76 syntax errors, #77 astropy BlackBody API) —
  lazy imports so broken corners cannot block `import ampere`.
- Branch harvest & archive: extract the useful diffs (swyft TMNRE ~550
  lines; `optim_only` optimisers; `jax` branch design sketches — note its
  `equinox.Dataset`/`equinox.Parameter` imports don't exist, and it carries
  ~50 MB of committed lightning checkpoints that must not be merged), then
  tag and archive stale branches.

### Phase 1 — Core contracts (§4)
- Write the interface specs as **layered design documents** under
  `docs/design/` — this plan stays a roadmap; the specs carry the detail:
  - `architecture.md` — the core + capability ladder, with the
    legacy-vs-reference-backend trade-off analysis written out;
  - `contracts/` — one spec per contract (parameters, result schema,
    transformations, likelihoods, inference, results), each with worked
    examples across multiple modalities and edge cases;
  - `lowering.md` — precisely how neutral declarations lower into each
    backend: the distribution mapping table (scipy ↔ torch.distributions ↔
    numpyro), transform conventions, module/pytree conventions (equinox on
    the jax side, with Paramax a candidate mechanism for constrained
    parameters), the parameter-vs-buffer lowering rules (§4.1), RNG
    lowering, and the device/dtype policy.
- ABCs in `ampere.core` + conformance suite skeleton. Start from the
  copilot-branch scaffold.
- **Paper-validate the contracts against every target modality** (spectra,
  photometry, visibilities+closure phases, astrometric time series, IFU
  cube, hierarchical population): one short design sketch per modality
  showing model→instrument→likelihood composition typechecks conceptually.
  No implementation — this is how the broad scope informs interfaces
  without exploding v1.
- Freeze the spec (version it; changes thereafter require explicit
  decision log entries).

### Phase 2 — Twin modern backends, lockstep (the big one)
Parallel agent tracks for `backends/torch` and `backends/jax`, both against
the frozen spec, both gated by the conformance suite. v1 slice scope:
- Models: blackbody, modified blackbody, power laws (ported natively);
  legacy-model adapter (black-box capability flags).
- Instruments: resampling, LSF, synthetic photometry, calibration factor.
- Likelihoods: iid Gaussian + flexible GP with `DenseGP` and `QuasisepGP`
  strategies (torch: GPyTorch or celerite2-torch; jax: tinygp quasisep).
- Inference: NUTS + VI natively (pyro / numpyro); emcee + dynesty via the
  §4.5 contracts from `ampere.inference`.
- Results: everything emits InferenceData; plotting from `ampere.results`.
- Diagnostics (§4.8): post-fit residual tests and GP-localisation plots;
  RHMF-style pre-fit screening if Robusta-HMF proves adoptable.
- Shared numeric utilities via the Python array-API standard
  (`array_api_compat`) where practical, to reduce duplication.

**Milestone M2 (flagship validation):** reproduce the
`flexible_likelihood_comparison` misspecification study on both backends at
10–100× the current data size, with wall-clock benchmarks vs legacy. This is
paper-grade evidence the redesign delivers its central promise.

### Phase 3 — SBI layer (backend-spanning)
- One SBI module in `ampere.inference` consuming `simulate()` from any
  backend (including legacy black-box models): NPE/NLE/NRE via `sbi`;
  revive the swyft TMNRE implementation from the branch; embedding-network
  support carried over from the existing `infer/sbi.py` work. jax-native
  SBI (sbijax/flowjax) optional, later.
- A canonical **coordinate–value–mask tensor encoding** of containers for
  embedding networks, so set/attention-based embeddings can consume any
  modality — irregular sampling and missing data included. Fixed-size
  summaries remain the simple default for a single fitting problem, where
  the data layout is fixed anyway; the encoding is what makes amortisation
  across differently-sampled datasets possible.

### Phase 4 — Extensibility proof: one new modality end-to-end
- Implement **interferometric visibilities** (decided) through the whole
  stack to prove the composition design: Fourier sampling at (u,v) points as
  a Transformation, complex Gaussian likelihood, closure phases. It stresses
  the schema (complex data) and the transformation chain hardest; other
  modalities (astrometric time series, IFU cubes) then follow the template
  it establishes.
- Astropy interop adapter (`core/astropy_compat.py`, §4.7).

### Phase 5 — Scale-out & advanced inference
- Approximate GP strategies for images/IFU (SVGP, SKI, Vecchia) behind the
  `NoiseModel` interface (GPJax is the natural provider on the jax side,
  GPyTorch on the torch side).
- Hierarchical/population inference: implement the container + hyperprior
  design from Phase 1 (plates in pyro/numpyro).
- **A dedicated benchmark-driven optimisation pass**: profile against the
  CI benchmark baselines established in Phase 2, then attack the levers in
  evidence order — requirements-negotiation compilation and caching (§4.3),
  vmap/batched evaluation, GPU placement, precision policy, solver
  selection, resampling (issues #12, #29, #67). No speculative optimisation
  before profiles exist.

### Phase 6 — Docs, migration, release
- Sphinx docs rebuilt around the new core; example gallery migrated;
  migration guide from legacy; deprecation policy for `ampere.{data,models,
  infer}`; beta release (addresses issues #57–60, #62).

### Cross-cutting workstream — CI/CD (grows with each phase)

Full CI/CD is in scope (decided 2026-09-01). Staged so each phase's gate
exists before the work it protects:

- **Phase 0**: GitHub Actions PR gate — ruff lint + format check, **pyrefly
  type checking** (scoped strictly to the new namespaces `ampere.core`,
  `ampere.backends`, `ampere.inference`, `ampere.results`; legacy is
  excluded and never gets typed), pytest (characterisation + unit tests),
  Python 3.11–3.13 matrix, coverage reporting; pre-commit hooks mirroring
  the linters. New code is typed from the first line — the contracts are
  exactly where static types pay for themselves.
- **Phase 1**: conformance suite wired into the gate; Sphinx docs build
  check with warnings-as-errors.
- **Phase 2**: backend matrix — separate torch, jax, and no-extras jobs
  (the last catches lazy-import breakage, i.e. `import ampere` must never
  require torch or jax); dependency caching; benchmark suite
  (pytest-benchmark or asv) with results tracked as CI artefacts. GPU tests
  are nightly/manual — hosted runners are CPU-only.
- **Phase 3**: SBI smoke tests with tiny simulation budgets.
- **Phase 6**: release automation — tag-driven builds via setuptools_scm
  (already configured in pyproject.toml), PyPI trusted publishing,
  changelog generation, versioned docs deployment (Read the Docs or
  gh-pages).

### Design horizon — capabilities to keep unblocked (hooks reserved now)

Not scheduled, but the contracts must not paint them out:

- **(a) Multi-fidelity modelling**: cheap and expensive variants of a model
  sharing one parameter space. Hook: fidelity tags on result channels
  (§4.2) and on model registrations; inference engines that mix fidelities
  (multi-fidelity BO/SBI) then slot in above the backends.
- **(b) Population models by post-processing**: hierarchical inference from
  archived single-object fits via importance reweighting under population
  hyperpriors (Hogg-style). Hook — **decided**: every run stores per-sample
  `log_likelihood` and `log_prior` (§4.6); a population module can then be
  built entirely on stored InferenceData files.
- **(c) Model emulation**: an emulator is just a Model trained on
  (θ, ModelResult) pairs that `simulate()` already produces; it satisfies
  the same contract (channels included) and is differentiable even when the
  original simulator is not — quietly upgrading black-box models to
  NUTS/VI. Hooks: serialisable training sets; model interchangeability;
  spec-hash invalidation of trained artefacts. Emulators may be
  **coordinate-conditioned** (functional emulation, neural-operator style),
  so an emulator trained once can be evaluated on any sampling of its
  output function — required for emulated models to participate in
  requirements negotiation (§4.3).
- **(d) Hierarchical SBI**: hook is plate-aware parameter groups (§4.1) and
  a `simulate()` that supports batched hierarchical draws.

Dependencies: 0 → 1 → 2 → {3, 4} → 5 → 6, with 3 and 4 parallelisable and
the CI/CD workstream running alongside every phase.

## 6. Deferred implementation choices

Settled at the start of the phase that needs them, not now:

- **Torch GP solver library** (Phase 2): GPyTorch structured solvers vs
  celerite2's experimental torch interface for `QuasisepGP`. Evaluate both
  against the conformance suite; GPyTorch is the safer bet for the broader
  Phase 5 strategies (SVGP/SKI).
- **Jax GP libraries** (Phase 2 gate check): resolved policy is
  *both-where-strongest behind the strategy interface*, not a single
  library. `QuasisepGP` comes from tinygp's `QuasisepSolver` or
  celerite2.jax (verify at Phase 2 start which is better maintained);
  GPJax — attractive for its Equinox/Paramax foundations and variational
  machinery — is the candidate provider for Phase 5 sparse/variational
  strategies, but **verify it actually offers a quasiseparable/state-space
  O(N) exact solver before letting it displace tinygp** (believed absent as
  of early 2026). Ampere's own §4.1 layer remains the only user-facing
  parameter interface regardless of provider, for cross-backend parity;
  Paramax is a lowering mechanism, not a user API.
- **SBI package set beyond `sbi` + swyft** (Phase 3): jax-native SBI
  (sbijax/flowjax) if and when maturity warrants.
- **Benchmark harness** (Phase 2): pytest-benchmark vs asv.

## 7. Known traps (do not rediscover these)

- Two model-output conventions caused a recurring bug class; the schema
  (§4.2) exists to end it. No backend code before the spec freeze.
- Zero test coverage on legacy: characterisation tests **before** anything
  else touches shared files; legacy stays frozen precisely because
  refactoring untested code is how regressions happen.
- Lockstep backends drift without a mechanical forcing function: the
  conformance suite is that function; a feature is "done" only when both
  backends pass it.
- torch + jax in one environment is dependency pain: lazy imports, extras,
  separate CI jobs. Never make `import ampere` require either.
- No binary artefacts in git (the `jax` branch's 50 MB of checkpoints).
- `eval`-based composition (`CompositeModel`) is replaced by astropy
  compound models (§4.7) and native torch/jax module composition — never
  reimplemented.
- **Equinox static fields are not buffers**: marking large constant arrays
  `static=True` hashes their contents into the JIT cache key — silent
  recompiles and memory blow-ups. Buffers lower per §4.1.
- **Float precision**: torch and jax default to float32, and GP Cholesky /
  quasiseparable solves in float32 fail in ways that look like science
  problems. Policy: float64 by default for likelihood linear algebra (jax
  needs the explicit x64 flag), with per-run opt-out for GPU throughput.
- **Trained-artefact caching** (SBI posteriors, embedding nets, emulators):
  cache keys must hash the model/prior/data spec so stale artefacts are
  invalidated automatically — the recent SBI caching bugs on master are the
  evidence this bites.
- **Units in the hot loop**: convert Quantities to canonical internal units
  once at composition time, never per-evaluation.
- Agent handoff discipline: work items sized to one agent session; every
  interface change goes through this document's decision log; PRs against
  frozen specs, not moving targets.
