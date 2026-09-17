# Ampere v2 Development Plan — DRAFT for refinement

Status: **plan settled; Phase 1 spec frozen** (freeze recorded 2026-09-03,
W1.13; the `spec-v1.0` tag is created at that item's merge); **Phases 0–3
complete** (Phase 2 closed 2026-09-08, Phase 3 closed 2026-09-10 — the §5
Phase 3 section carries the landed summary); **Phase 4 complete
(2026-09-11 to 2026-09-13; the §5 Phase 4 paragraph carries the landed
summary; decisions D1–D4 are recorded in `docs/design/phase4_placement_memo.md`);
Phase 5 is next to draft.** All
architectural proposals are confirmed; remaining open items are
implementation-level choices deferred to their natural phase (§6). This
document is the source of truth for the redevelopment of ampere: decisions
taken, target architecture, and phased roadmap. The agent-sized work items
for Phases 0–3 live in `WORK_ITEMS.md`.

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
| Flexible-likelihood kernel | **Matérn-class (Matérn-3/2 / sums of SHO terms) becomes the default throughout**, replacing the legacy RBF — it represents structured residuals better and is exactly O(N) via quasiseparable solvers. Validated against the misspecification study in milestone M2 — **done, W2.10 merged 2026-09-08**: the flexible likelihood stays within one posterior width of the truth with the truth covered on every scenario and every rung of the 200→20 000 ladder while the standard likelihood's bias grows as √N to 113 widths; see the W2.10 status row and `docs/source/m2_misspecification.rst`. (Validation note ratified by Peter 2026-09-08.) |
| Packaging | **Single distribution with extras** (`ampere[torch]`, `ampere[jax]`, `ampere[sbi]`, …), lazy imports; split distributions are not pursued. |
| Phase-4 modality proof | **Interferometric visibilities** — complex-valued data and Fourier sampling stress the result schema and transformation chains hardest; other modalities follow the template it establishes. |
| CI/CD | **Full CI/CD is in scope**, as a cross-cutting workstream (§5) growing from lint+tests (Phase 0) through conformance/backend matrices and benchmark tracking (Phases 1–2) to automated PyPI releases and versioned docs (Phase 6). Includes pyrefly type checking scoped to the new namespaces. |
| Reference backend | **Confirmed**: pure-numpy `backends/reference` as conformance oracle, numpy-only base install, and execution venue for astropy-adapted models (§3). |
| Result schema | **Confirmed**: named channels + requirements negotiation; containers are coordinate-indexed function samples with first-class masks; no regular-grid assumptions anywhere in the contracts (§4.2–4.3). |
| Non-Gaussian robustness | **Accepted**: flexible-GP robustness for non-Gaussian likelihood families (latent-GP path) is promised on the torch/jax rungs only (§4.4). |
| Buffers | The model contract distinguishes parameters from **buffers** (constant arrays); lowering rules per backend in §4.1 (torch `register_buffer`; jax via partition filters, never equinox static fields). |
| Diagnostics | **Misspecification diagnostics are in scope** (§4.8): RHMF-style pre-fit screening (Hilder et al. 2026), post-fit residual tests, GP-based localisation. |
| Dependency management | **pixi replaces manual environment control** (work item W0.8): pyproject.toml stays the single source of truth for dependencies; pixi provides locked environments (features mirroring the extras), tasks, and the CI environment setup. Adopt before/with W0.6 so CI is built on it once. |
| Dependency pins | `pyphot<2` and `sbi<0.28` pinned 2026-09-01 as **temporary** measures (pyphot ≥2 removed `pyphot.unit`; sbi 0.27 changed `posterior.map()` shapes) — **both lifted by W0.9 (2026-09-03)**. pyphot: `ampere/utils/pyphot_compat.py` targets the ≥2 unit-adapter API as its primary surface (astropy-backed `pyphot.config.units`), with 1.x still supported through the same shim; `ampere/data/photometry.py` migrated onto it. sbi: `ampere/infer/sbi.py`/`mixins.py`'s W0.5 fixes (batched `posterior.map()` shape, `_prior_is_normalised` default) hold against the currently-resolved sbi; PyPI's latest sbi release is still 0.27.0 as of the migration (no 0.28+ exists yet), so the upper pin was simply removed rather than validated against a newer release — re-check when sbi actually publishes ≥0.28. **Addendum 2026-09-07**: pyphot 2.1.1 imports `requests` without declaring it, which broke `import pyphot` under a bare `pip install ampere` (the minimal-install CI job's finding, carried from W2.1). Ruled by Peter: `requests` is declared in `[project.dependencies]` on pyphot's behalf, for now, with the pyproject comment recording that ampere itself makes no HTTP requests and that httpx, aiohttp or niquests would be the better choices should it ever need a client of its own. Drop the line when pyphot fixes its metadata. |
| Curated astropy→native translation | **Opt-in only, never silent** (ruled 2026-09-01 during the W1.2 review, amending this document's original §4.7 "silently restoring differentiability" wording). The default adapter path always wraps the user's actual astropy model as a black box; a curated native equivalent is substituted only on an explicit, backend-scoped request — sketched as a `from_astropy()` constructor on the backend subpackage, which raises if the model (or any component of a compound model) is not fully in the curated table, rather than silently falling back to black-box. Exact API fixed with the adapter contract (Phase 4). |
| jax non-trainable mechanism | **`eqx.partition` filter specs** (ruled 2026-09-01 at the W1.9 review, ranking what §4.1 originally left unranked): buffers and fixed parameters lower to ordinary array leaves excluded from the trainable partition via an explicit filter spec — not `paramax.NonTrainable`, whose freezing happens only when `unwrap()` is called, so a forgotten call silently trains the buffers; and never equinox static fields (§7). paramax remains an interop layer at the boundary if Phase 2's GP library choice (GPJax) puts wrapped leaves there. Analysis in `docs/design/lowering.md` §6.2. |
| Joint-space merge topology | **Nested `ParameterMapping`, one merge per level — ratified** (ruled 2026-09-02 at the W1.7 review, closing the question `parameters.md` §14 expressly kept open): every composite performs one merge over its immediate children and retains the mapping; values re-distribute level by level; cross-level ties collapse correctly. The flat-leaf-merge alternative and the full case are `inference.md` §4; the lossless-nesting sub-proposal and `Binding.index` were still open when this was recorded — both approved and landed later the same day (next row). |
| `LikelihoodFamily.sample` | **Granted and landed** (ruled 2026-09-02, W1.7's R3 — a §4.4 interface addition, hence this entry): families gain an optional generative half, `sample(predicted, noise, rng)`, drawing from the same distribution `log_prob` scores with the same `NoiseParams`. The default refuses specifically, naming the family and the override to provide; `GaussianFamily` implements both noise models; a user family overrides it for an exotic observation process. `simulate(observe=True)` delegates to it. |
| Channel-name surface | **Dot-qualified channel names, settled now** (ruled 2026-09-02, W1.7's R4): `_check_channel_name` and the `Instrument` channel binding accept `.`-separated identifiers — qualification happens in the model, never the `DatasetCollection`, and grouped non-hierarchical data (one object's several sub-mm CO lines) motivate it as much as population models. Merge-component labels (steps, datasets) stay bare identifiers. |
| Lossless nesting & plate routing | **Approved 2026-09-02** (the R1 follow-ons): `ParameterSet.merge` accepts a `ParameterMapping` as a component and composes its bindings ("lossless, not necessarily recursive" — closes `inference.md` §4.6's introspection defect at source), and `Binding` gains an optional `index` so `distribute` routes one element of a plate's array-valued parameter to one component (`hierarchical_population.md` H-2 — what makes a `Plate` feed per-object datasets). Both amend §4.1's contract; **implemented and merged the same day** — `parameters.md` §8 carries the contract with executed examples, and W1.7's `Dataset`/`FittingProblem` pass their mappings so bindings compose end to end. |
| Likelihood data checks (I-1, I-5, X-2, X-3) | **Approved as proposed and landed 2026-09-02** (W1.11 amendments): `check_alignment` compares value dtype kinds (a complex prediction cannot silently be fitted against real amplitudes); `LikelihoodFamily.check_observed(observed)` is the family's composition-time hook (observed container only — value-range properties differ between prediction and data), where `PoissonFamily`'s integrality test now lives; `NoiseParams.retain` hands families the caller's inclusion indicator so their own aligned arrays excise correctly under masking. X-1 (prediction-aware `NoiseModel`) **accepted 2026-09-03** — `awkward_instrument.md` §6's detailed design as written: keyword-only `predicted=None` on `sigma`/`noise_params`, passed at every call site (`log_prob`, `conditional`, `draw_observation`); an outright pre-freeze signature change with no shims; marginalisation declarations untouched; `FractionalModelNoise` named in `likelihoods.md` §5 and implemented with the reference backend. Lands at W1.13 with its three conformance rows. |
| Likelihood failure signalling | **Non-strict engine path with a `strict` toggle** (ruled 2026-09-02 at the W1.6 §17 Q1 review): sampling-time evaluation failures reach the engine as −inf with a recorded reason; strict raising remains for direct use and debugging (see the amended §4.5 bullet). Decided in the end by jax — exception control flow does not trace, so the non-raising path is forced by Phase 2 regardless. W1.7 owns the conversion and recording; the `DenseGP`-level flag follows once that recording mechanism exists. |
| Effective-mask resolution | **The `Dataset` resolves the effective (predicted ∪ observed) mask once at construction** (ruled 2026-09-02 at the W1.6 §17 Q2 review), not `Likelihood` per evaluation. The effective mask is thereby a declared evaluation-time invariant — parameter-dependent output masks are unsupported in a fitting problem and W1.7 checks this loudly — making explicit what `latent_declaration(n)`'s fixed shape already assumed. `likelihoods.md` §8's mechanics are unchanged; a pre-resolved pair makes its internal union the identity. |
| Results-contract rulings (W1.8's R1–R7) | **All accepted as recommended** (ruled 2026-09-03, `results.md` §15): arviz + a netCDF engine join the base install with Phase 2's engine drivers (R1); `Likelihood.pointwise_log_prob` granted at the freeze, not stored by default — W1.13 lands the §4.4 addition (R2); ArviZ 1.0's `InferenceData` → `xarray.DataTree` wording corrected in §4.6 and `architecture.md` §3, format unchanged (R3); `AnomalyScore` lands in `ampere.core` at W1.13 (R4); `ResultsError` moved to `ampere/core/exceptions.py` same day with `ampere.results` re-exporting, and container serialisation stays functions (R5); the per-dataset `log_likelihood` group shape confirmed (R6); `Likelihood.to_spec()` folds into W1.13's consolidated cross-contract serialisation review (R7, with `likelihoods.md` §17 Q8). |
| Discrete parameters | **Refuse the bijection, keep the door ajar** (ruled 2026-09-03, `lowering.md` §12 Q1): `default_bijection_for` raises a typed *capability* refusal for a discrete family instead of returning a continuous bijection that cannot be right. Declaration, prior sampling, `prior_transform` (scipy's discrete families implement `ppf`) and lowering-as-distribution are untouched, and discreteness becomes queryable from the canonical prior description — so the non-gradient routes that might eventually support discrete parameters ((variational) EM, numpyro-style enumeration, SBI, nested sampling, Bayesian optimisation) stay reachable, though none is in the current plan. W1.13 lands the §4.1 change. |
| User-registered lowerings | **Registry accepted in principle, hardened** (ruled 2026-09-03, `lowering.md` §12 Q8): `register_lowering(family, backend, constructor)` keyed on the neutral name, plus the analogous per-backend `Bijection` slot — no silent overwrites (`override=True` required), user-registered rows stamped in provenance, and an opt-in conformance battery for registrants. The registry resolves at lowering time before any tracing, so it is inert to `jax.jit`/`vmap`/gradients; constructors must return trace-pure objects. Plumbing lands in Phase 2 beside the backends. **Landed W2.6**: `ampere.core.lowering` (`register_lowering`/`lookup_lowering`, `register_bijection_lowering`/`lookup_bijection_lowering`, `run_registrant_battery`, `provenance_entries`); consumed by W2.4/W2.5 when the torch/jax backends register their own built-in rows through the same mechanism. |
| 2026-09-03 ratification batch | `substream` ratified in `ampere.core` (`lowering.md` §12.7 = `inference.md` §19.5); capability flags promoted into W1.5's ABCs at the freeze, defaults reproducing the `getattr` semantics (§19.6); the narrow failure catch-set default confirmed (§19.7); distinct instrument labels required when more than one instrument reads a channel, checked at problem composition, the singleton default unchanged (`transformations.md` §15 Q4 residual); two-channel instruments deferred, sketch-first (Q6); W1.10's narrowed hash claim stands — spec/data hashes are the backend-spanning promise, the problem hash stays backend-variant, backend-neutral model identity routed to the freeze beside `describe()`; CI expansion authorised as W0.10, with py3.11 droppable if it proves problematic. |
| Transformation-chain rulings (W1.5's §15 Q2/Q3/Q5 + the overlooked slots) | **All landed at W1.13** (ruled 2026-09-03; §4.3 contract changes). No requirement pull-back — gap I-3's chain-internal `Transformation.configure_from(downstream)` is the adopted mechanism, called once per step at `Instrument` construction, default no-op, push-forward-and-raise the posture. `compile_for` engaging with an unachievable requirement raises `CompositionError` by default (gap I-4's loud option); the explicit opt-out is `FittingProblem(lenient_compile=True)`, downgrading the refusal to a warning and proceeding with the unconfigured model. `Instrument.freeze()` snapshots the merged parameters and *refuses* later step reconfiguration (never a silent cache); `Dataset` freezes its instrument at construction, removing the per-`log_prob` re-merge `inference.md` §19.8 evidenced. And the image/1-D-spatial standard-library slots the roadmap overlooked (PSF convolution, spatial resampling, affine + the future WCS carrier on `Image`/`Cube`, Hankel, NUFFT) are recorded in `transformations.md` §10 as Phase 2+ interfaces the freeze precludes none of. |
| Interferometric likelihood interfaces | **Rice, von Mises and the circular complex GP fixed** (ruled 2026-09-03, `likelihoods.md` §17 Q3/Q4/Q6 — the W1.11 interferometry sketch's recommendations accepted as written; landed at W1.13). Rice: the model predicts the underlying complex value and an `Amplitude` instrument step takes the modulus, so the family receives real non-negative amplitudes. von Mises: `κ = 1/σ²` per sample from the container's own uncertainties. `complex_gaussian` + `GaussianProcessNoise` declares **`ANALYTIC`** with the circular (equal-component, zero-pseudo-covariance) complex GP as the fixed meaning — a §4.4 declaration change — while the implementation stays Phase 4's: composition refuses with the schedule named (`GP_ANALYTIC_IMPLEMENTED`, the declared-but-staged discipline), never a silently different model. Family implementations are Phase 4's, with the visibility modality. |
| Consolidated serialisation review | **One review, not piecemeal — done at W1.13** (ruled 2026-09-03; supersedes `likelihoods.md` §17 Q8, `results_schema.md` §17 Q6 and `results.md` §15 R7 as separate questions). The review is `docs/design/serialisation_review.md`. Its rule: **specs describe declarations (on the objects, in core), provenance fingerprints content (one-way, in `ampere.results`), storage carries values (netCDF)**. It confirmed R7's promotion — `Likelihood.to_spec()` landed as a §4.4 addition, `describe_likelihood` now composes it — pinned the design-horizon (c) training-set format (the pair form and the spec-hash invalidation key are fixed; the writer is Phase 2's), recorded the deliberate absence of instrument/model reconstruction, and found and fixed one provenance hole (family/noise-model buffers were invisible to the likelihood fingerprint; `PROVENANCE_SCHEMA_VERSION` bumped to 2). |
| WStat / profile likelihoods | **Not shipped; documented as a compared workaround** (ruled 2026-09-03, `awkward_instrument.md` §9 Q3): ampere's standard library never carries the profiled Cash-with-background statistic — the docs are deliberately opinionated that the two-dataset Bayesian formulation is the correct approach — but a worked example shows the user-family route (safe under masking since X-2's `retain`; its `sample` refuses; per-sample log-likelihood/LOO semantics degrade) beside the joint fit, comparing the two's pros, cons and results. Lands with Phase 2's engine drivers; W1.13 carries it into the Phase 2 breakdown. |
| Freeze escalations (model identity, `Population`, `Axis.locate`) | **All three ruled 2026-09-03**, closing W1.13's escalations: (1) the `describe()` hook and the derived backend-neutral model identity are adopted for **early Phase 2** as one mechanism, landing with W2.1 — the neutral identity *offers* a cross-backend emulator, never serves one silently, and `ampere_problem_hash` stays backend-variant (`results.md` §13.13/§14). (2) `Population` (H-1) lands with **Phase 5**, timing deliberately adaptable — nothing frozen precludes it, and the tie-based pattern is the documented route until then; no primary population route is committed (`hierarchical_population.md` §11 Q2/Q5). (3) `Axis.locate` — the lookup is the approved approach, landing with W2.1 (`spectrum_photometry.md` Gap 1). The W2.1 landings carry their own decision-log entries per ground rule 9. |
| **Phase 1 spec freeze** | **The §4 contracts are frozen** (W1.13, 2026-09-03; the `spec-v1.0` tag is created at the W1.13 merge, after Peter's review). Frozen surface: `docs/design/architecture.md`, the seven contract specs under `docs/design/contracts/`, and `docs/design/lowering.md`, each at v1.x **as amended through the freeze** — the amendments being the recorded rulings W1.13 implemented (X-1's prediction-aware `NoiseModel`; the circular complex GP declared `ANALYTIC` with Rice/von Mises interfaces fixed; `configure_from`, `Instrument.freeze()` and the loud `compile_for` with `lenient_compile` as the opt-out; capability flags promoted into the ABCs; the instrument-label requirement; the discrete-family `CapabilityError`; `LoweringError` + the exception ratifications; the icdf strict-flag home; `pointwise_log_prob`; `AnomalyScore`; `Likelihood.to_spec()` with the consolidated serialisation review) plus the cross-review harmonisation and the W1.11 gap dispositions. From this point ground rule 9 is in force: **any change to a §4 contract requires a decision-log entry in this table in the same PR.** The Phase 2 work-item breakdown (W2.1–W2.11) is written against this frozen spec in `WORK_ITEMS.md`. |
| jax x64 activation | **Guard-and-raise, never set-on-import** (ruled 2026-09-01 at the W1.9 review, amending `architecture.md` §5's original set-on-first-import sketch — the *policy*, float64 always for likelihood/GP linear algebra, is unchanged). Ampere never flips `jax_enable_x64` as an import side effect; `ampere.backends.jax` ships an explicit, idempotent `configure_x64()`, and construction of any jax-backed likelihood/GP/model raises when the flag is off, naming the three remedies (the `JAX_ENABLE_X64=1` environment variable; the explicit call; or a per-run, provenance-recorded `float32` opt-out). numpyro needs no separate switch — its `enable_x64` is a verified thin wrapper over the same jax flag. Analysis in `docs/design/lowering.md` §10.2. |
| `Parameterised.describe()` + backend-neutral model identity | **Landed with W2.1** (the escalation ruled 2026-09-03; a post-freeze §4 addition, hence this entry — `results.md` §9/§13.13/§14, `parameters.md` §2/§10). One mechanism, in two halves. (1) `Parameterised.describe()` is an **opt-in** hook returning a normalisable mapping, or `None` (the default); `provenance.model_fingerprint` folds it in under a `describe` key, closing §13.13's hole — a `Redden(law="ccm89")` configured by a plain attribute that is neither parameter nor buffer no longer shares a cache key with `Redden(law="f99")`. Hashing an arbitrary `__dict__` stays rejected: it would sweep in caches, file handles and unhashable state. (2) `provenance.neutral_model_identity` is the fingerprint **minus `class` and `module`**, with `model_identity_hash` its digest, recorded as `ampere_model_identity_hashes`. It **offers, never serves**: a match licenses proposing a cross-backend emulator with its provenance shown, and never silently substituting one, because the neutral identity cannot pin the mathematics. `ampere_problem_hash` stays deliberately backend-variant and W1.10's equivalence row is unchanged; a new conformance row asserts the neutral identity *agrees* across the two fixtures, which is the other half of the same claim. Two consequences recorded deliberately: `describe` is now a class attribute, so it joins `parameters`/`buffers`/`context` as a name a parameter or buffer may not shadow; and adding a fingerprint key changes every `ampere_problem_hash`, so `PROVENANCE_SCHEMA_VERSION` rides to **3**. Scoped to `model_fingerprint` as §14 words it — `_describe_instrument` records each step's class too, so a whole-*problem* neutral identity would need that second site and is not claimed here. |
| `Axis.locate(values)` | **Landed with W2.1** (approved 2026-09-03; a post-freeze §4 addition, hence this entry — `spectrum_photometry.md` Gap 1, `results_schema.md` §2/§5, `transformations.md` §10). Index lookup on an `Axis`, matching within `COORDINATE_RTOL` and raising `SchemaError` naming any value with no match. It exists because a step whose buffer is tabulated on the coordinates it published as `points=` has no supported way to find them again once `negotiate`'s union hands it a larger, reordered grid — the failure being a bare `matmul` shape error naming neither channel nor negotiation. The tolerance is the point rather than a convenience: the union collapses coordinates coinciding to within `COORDINATE_RTOL` and keeps one representative, so a step's own published value may legitimately differ from the survivor. `COORDINATE_RTOL` therefore **moved from `transform.py` down to `results_schema.py`** (re-exported from its old home, so both import paths keep working): `locate` is the exact inverse of `_dedupe`'s collapsing and the two must not be able to drift apart. The axis is not assumed sorted — `PhotometricPoints` and `VisibilitySet` declare `Order.ANY`, and Gap 1's own scenario is a photometry step. W2.1's `SyntheticPhotometry` is the first consumer; `transformations.md` §10 now records the pattern as **required** for any step whose buffer is tied to specific `points=`. |
| `QuasisepGP` on the reference path | **Landed with W2.3, with its own celerite term and one recorded deferral** (2026-09-05; a post-freeze §4.4 landing, hence this entry — `likelihoods.md` §7). Three decisions. (1) **The implementation lives in `ampere.core`, beside `DenseGP`, not in `backends/reference`.** The class was already declared there; `GPSolver.IMPLEMENTED` is documented as "whether an implementation exists in the reference (numpy) path"; and `backends/reference/__init__.py` states in as many words that "kernels, GP solvers, families, containers and the fitting problem itself are **not** here: they are backend-neutral and live in `ampere.core`", which is `inference.md` §18's "a backend supplies models and transformations […] and nothing else". `architecture.md` §2/§3's "celerite2's numpy interface is a core dependency of `backends/reference`" is read as the *packaging* statement it makes — celerite2 joins `[project.dependencies]`, never an extra — not as a module address. celerite2 is imported **lazily**, inside `QuasisepGP`, so `ampere.core` keeps the numpy/scipy/astropy/stdlib import surface it promises; same idiom as the reference backend's pyphot import, and the ~20 ms is not charged to a problem with no GP in it. (2) **ampere supplies its own celerite term rather than using `celerite2.terms.Matern32Term`,** which is an approximation controlled by `eps` (the celerite basis has no `τe^{−cτ}` member) and misses `tolerances.cross_solver` by three orders of magnitude at its default. Matérn-3/2 has an *exact* rank-2 representation in the semiseparable form celerite2's solver actually factorises — `U_n = (a²(1+f t_n), −a² f)`, `V_m = (1, t_m)`, `c = (f, f)`, `f = √3/ℓ` — which is an algebraic identity, not a limit, and is what makes §4.4's "exact O(N)" true. Its generators grow linearly in the coordinate, so ampere re-references coordinates to the midpoint of their own range and the accuracy degrades gracefully with the number of length scales spanned (~5e-12 over 10, ~2e-8 over 10⁴ against the dense Cholesky). A kernel declaring `QUASISEPARABLE` without a registered exact representation is refused by name, never silently approximated. `QuasisepGP` gains an optional `jitter`, the same knob, meaning and default as `DenseGP`'s, since the same ill-conditioning failure mode exists. (3) **`GPSolver.conditional_loo` is DEFERRED on the quasiseparable path**, and its refusal stays live and tested. Every leave-one-out term needs `A_ii` for `A = (K + diag(σ²))⁻¹`, and celerite2's public numpy interface has no O(N) route to that diagonal — its own conditional variance forms the cross-covariance densely, so it is O(N·M). An O(N) route does exist (a backward accumulation of `Σ_{k>i} U'_k U'^T_k / d_k` over the inverse of celerite's unit-lower semiseparable factor, or equivalently a Kalman smoother on the Matérn-3/2 SDE; and given an O(N) posterior variance `v_i` at the data points, `A_ii = (σ_i² − v_i)/σ_i⁴` closes it), but each reimplements celerite2's internal factorisation convention in numpy — a coupling to undocumented internals, for a decomposition `results.md` §15 R2 does not store by default and which `DenseGP` already computes exactly. Shipping an O(N²) fallback under an O(N) name was rejected outright. `condition` is likewise O(N·M) rather than O(N), because a dense cross-covariance block is unavoidable when M outputs each need all N inputs; the docstring and `likelihoods.md` §7 say so rather than implying otherwise. |
| Detector-aware synthetic photometry | **Both conventions, per filter, one convention across backends** (ruled by Peter 2026-09-05 on the W2.1 branch). `SyntheticPhotometry` weights an `f_nu` container by `R dlambda / lambda` for a **photon**-counting detector and `R dlambda / lambda**2` for an **energy** detector; both are normalised weighted means, so they differ on a sloped spectrum and agree on a flat one. The photon form is pyphot-consistent by construction, not by coincidence: pyphot's photon convention is the lambda-weighted mean of `f_lambda`, and converting that to `f_nu` at the pivot (`lambda_p**2 = INT R lambda dlambda / INT R dlambda/lambda`) collapses exactly to it. `detector=` is a **required** keyword on the raw constructor — a silently chosen convention is the failure being guarded against, and the class is new so there is no compatibility cost — accepting one string for all filters or one per filter, since real filter sets mix types (the bundled library is 819 photon to 64 energy: 2MASS counts photons, AKARI and IRAS measure energy). `from_library` fills it per filter from pyphot's own `Filter.dtype`, with an explicit `detector=` overriding, and refuses loudly rather than guessing when a filter's metadata is missing or unrecognised. The choice is visible to provenance through a registered `weights` buffer and legible through `describe()`. **The reference backend is the conformance oracle, so this is the standard torch and jax must reproduce** — the battery is the enforcement. The plain `R dlambda` weighting an earlier draft of this class used is neither convention and was replaced, not kept as a third option; it was never released.
| arviz + a netCDF engine in the base install | **Landed with W2.2** (the promotion ruled 2026-09-03, `results.md` §15 R1: arviz joins the base install "with Phase 2's engine drivers, promoted together with a netCDF engine, when a user can first emit a run"). That moment is the engine drivers: `ampere.inference` gives every run an ArviZ `DataTree` emitted through `ampere.results`, so a base install that could not build one would be a base install that cannot sample. `arviz` and `h5netcdf` move into `[project.dependencies]` and the `arviz` extra is **removed** rather than retained as an empty alias -- `pip install "ampere[arviz]"` failing loudly with "does not provide the extra" teaches the right thing where a silent no-op teaches the wrong one. h5netcdf rather than netCDF4 because arviz requires neither and one is needed to *write* a run (precisely the gap R1 names), and it is the lighter of the two; the pixi `arviz` feature becomes `netcdf`, carrying only the second engine that `tests/results`' both-engines row `importorskip`s. `architecture.md` §3's extras table is updated to match, which is why this entry exists: the table is part of the frozen surface. Two consequences recorded deliberately. (1) `ampere.results` keeps its lazy import behind `OptionalDependencyError`, but the justification changes from "optional dependency" (`architecture.md` §4 rule 2) to cost plus diagnosis -- arviz pulls xarray and pandas onto every path that touches the namespace, including a caller who only wanted a provenance hash, and an environment assembled without it deserves a named remedy; the error now carries `extra=None`, so it says `pip install arviz`. (2) The CI matrix installs arviz on every leg now, so `tests/results` runs there instead of skipping, and `tests/inference` joins it. |
| `configure_from` call order | **Last step first** (ruled by Peter 2026-09-07 on W2.1's out-of-scope finding; a §4.3 clarification of the 2026-09-03 chain ruling, hence this entry — `transformations.md` §5). `Instrument.__init__` calls each step's `configure_from(successors)` walking the chain **from the end**, so every step reads successors already in their final state. The ruling text left the order unstated and the implementation walked forwards; that under-padded chained same-axis convolutions, because the reference backend's `LSFConvolution` publishes nothing until it has learned the range downstream of it and a first convolution therefore read the second one *unconfigured* — it saw the resampler's bare range and padded for one kernel where two were needed, so the second kernel's outermost outputs were convolved against an edge. The hook's signature, its once-per-step call, the default no-op and the push-forward-and-raise posture are unchanged; only the order is now specified. Regression tests at both levels: a contract-level one in `tests/core` (a step whose requirement derives from its successors', the reach compounding) and the two-LSF case in `tests/backends`. torch and jax inherit the fix through `Instrument`; a backend that overrides chain construction must preserve the order. |
| Backend identity on `FittingProblem` | **The backend is §4.5's fourth capability flag** (decided by Fable 2026-09-07 on W2.2's carried finding; ratified by Peter 2026-09-08 — a post-freeze §4.5 addition, hence this entry). No §4 surface named a problem's backend, so `Engine(problem, backend=...)` could only *declare* one and a run's `ampere_backend` recorded that declaration rather than a fact — and the two lockstep backend tracks would each have invented their own answer. `Capable`, `Model` and `Transformation` now carry `BACKEND: ClassVar[str] = "reference"` beside the three existing flags (the conservative default: the base install's whole toolkit *is* the reference backend, and a hand-written numpy model runs on the reference path), every piece `ampere.backends.reference` ships declares it explicitly, `Capabilities.backend` is validated non-empty like `device`, and `declared_capabilities` aggregates it by the **device rule** — all parts agree or `DatasetError` names the backends and the `capabilities=Capabilities(backend=...)` override, because ampere converts arrays between libraries no more readily than between devices. `FittingProblem.backend` joins the other three properties. **One name per backend, everywhere**: the same string keys `lowering.md` §12.8's registry and ids a conformance fixture. `Engine.__init__` drops `backend=` (pre-release; no shim) and reads it off the problem; `emit`/`provenance_attrs` take it as an optional cross-check that raises on disagreement rather than recording it. `Capabilities.to_dict()` gained a key, so `PROVENANCE_SCHEMA_VERSION` is **4** (`results.md` §9); `capabilities` is not an input to `problem_fingerprint`, but the schema constant is, so problem hashes moved as at every previous bump. W2.6's deferred first-class lowering-provenance key is **not** this bump and stays deferred to W2.4/W2.5. Landed W2.12. |
| Realisation surface (W2.13) | **A registered, checked-on-use, differentiable native form of the problem** (all twelve sub-decisions ruled by Peter 2026-09-07 on Fable's recommendations, after a prototype ran; a post-freeze §4.5 addition — `inference.md` §10a). (1) Mandatory surface is minimal — `backend`, `free_size`, `log_prob_unconstrained` — with an optional `log_likelihood_terms`; absent, drivers recompute the decomposition for stored draws only and say so. (2) A realised density returns a bare −inf; §11's recorded reasons are recovered post hoc for stored draws. (3) `strict=True` on the native path means the factory refuses at construction; runtime failures are always −inf. (4) `realise` checks agreement with the numpy path at one point on use; the conformance suite compares many points per registered realisation. The numpy path stays the oracle. (5) The reference backend registers no realisation in v1 — 'realisation' means the differentiable native form; a numpy fast path for the reference hot loop is a later performance item. (6) W2.13's coverage floor: Gaussian families, independent/dense-GP noise, masks, plates, hierarchical priors, native kernels; the rest refused by name and widened in each track's slice 2. Fold-ins ruled with it: (7) the capability flags widen to `NoiseModel` and `GPSolver` (four ClassVars, reference defaults) and the likelihood contributes its noise model and solver to `capability_parts` — `LikelihoodFamily` stays flag-free — so a numpy solver in a native problem makes it non-differentiable and a backend disagreement, loudly; (8) one shared `ampere.core.exceptions.LoweringFallbackWarning` replaces both backends' local classes; (9) public `ParameterSet.evaluation_order()` and `Dataset.effective_mask`; (10) `GPSolver.provenance_config()` — backend-specific solver configuration (dtype, device, a future precision opt-out) recorded in the run's attrs, never in the spec hash, so cross-backend spec-hash agreement survives; (11) `PROVENANCE_SCHEMA_VERSION` 5, taken once here: `ampere_realised` and the consulted registered lowering rows — W2.6's deferred first-class key lands now that a backend drives lowering end to end; (12) three spec corrections from the tracks ride along: `lowering.md` §3.6's numpyro `Gamma`/`Beta` icdf need TensorFlow Probability; §4's `biject_to(lowered.support)` advice gains the declared-support caveat (torch's `TransformedDistribution.support` is the last transform's codomain); numpyro's `AffineTransform` defaults `domain` to `real`. Also ruled the same day: CI stays CPU-only for torch/jax, with API-level GPU smoke tests later. |
| Torch GP solver library (§6's deferred choice) | **celerite2's compiled kernels under `torch.autograd`, written by ampere; GPyTorch rejected** (decided 2026-09-07 on W2.4 slice 2's measurements — §6 asked for the two candidates to be *measured* against the conformance suite rather than chosen from documentation, and this row is that measurement). Neither candidate existed in the form §6 supposed. (1) **celerite2 ships no torch interface** (0.3.3 has `jax`, `pymc`, `pymc3`, `theano`), but it does ship `celerite2.backprop` — the compiled forward *and reverse* passes of the semiseparable factorisation, which is exactly what all four of those wrappers are built on. So `ampere.backends.torch._celerite` registers them as `torch.autograd.Function`s, the direct counterpart of `celerite2/jax/ops.py`. (2) **GPyTorch has no quasiseparable operator at all**: `linear_operator` 0.6.1's structured operators are Toeplitz (a *regular* grid), Kronecker (a product grid) and SKI / inducing points (approximate). ampere's coordinates are irregular by contract (`results_schema.md` §16 forbids a regular-grid assumption), so GPyTorch's only *exact* route on this problem is a dense Cholesky — i.e. `DenseGP`, already shipped, under another name. The table, at 10³ irregular points against `ampere.core.DenseGP` as oracle: agreement — celerite2+autograd **agrees with `ampere.core.QuasisepGP` to the last bits** (0.0 relative difference at this size; ~1e-15 in general, the two being different summation orders over one factorisation, and 1.7e-13 from the dense oracle, which is the dense solver's own accumulation error), GPyTorch-exact 1.4e-15, GPyTorch-SKI **9e-3 to 2e-2 and non-deterministic between runs** — three to four orders outside `tolerances.cross_solver` (1e-6). Gradients in (amplitude, length scale) match central differences to 8 significant figures on both viable candidates. Wall clock for value **and** gradient: 0.0056 s / 0.021 s / 0.144 s at 10³ / 10⁴ / 10⁵ for celerite2+autograd (exponent ~0.84 over the last decade, i.e. linear), against 0.69 s / 117 s / infeasible for GPyTorch — **123× at 10³ and 5 657× at 10⁴**, and at 10⁵ the dense covariance alone needs 80 GB. Dependency cost: **zero** for celerite2 (already a base dependency since W2.3 — so the `torch` extra is unchanged and no lockfile moves), against 5.2 MB of new pure-Python for gpytorch + linear_operator + jaxtyping. A third, control candidate — the celerite recursion written in pure torch — was measured too and rejected: the factorisation is a sequential scan over N with a J×J state and torch has no scan primitive, so it is a Python loop building an N-node autograd graph, **1 800× slower than the wrapper at 10³ points** and worse with N. The price of the compiled kernels is stated in the capability flags rather than hidden: they are CPU float64 C++, so `ampere.backends.torch.QuasisepGP` declares `DEVICE = "cpu"` and `BATCHABLE = False`, refuses a `device=` or float32 request by name, and reports `{dtype, device, library}` from `provenance_config()`. `DEVELOPMENT_PLAN.md` §6's bullet is struck accordingly. |
| `QuasisepGP.conditional_loo` on the torch path | **The O(N) recursion W2.3 deferred is supplied, on the backend that already pays for the coupling** (decided 2026-09-07 on W2.4 slice 2; a change of circumstances, not a reversal — the 2026-09-05 deferral stands exactly where it was taken). W2.3's reason was precise: the leave-one-out terms need `A_ii` for `A = (K + diag(σ²))⁻¹`; celerite2's **public numpy** interface exposes no O(N) route to that diagonal; and an O(N) route exists but "reimplements celerite2's internal factorisation convention in numpy — a coupling to undocumented internals". The torch backend calls `celerite2.backprop` directly, so it holds `t, c, U, d, W` in hand and the coupling is already paid for and, more importantly, **tested**: `tests/backends` asserts the factorisation against an independent torch transcription and the conformance row asserts the terms against `DenseGP`'s Cholesky. The recursion, recorded so it need not be rederived: with `M = L⁻¹` (unit lower triangular), `A_ii = Σ_{k≥i} M_ki²/d_k`; column *i* of `M` solves `L z = e_i`, and substituting `e_i` into celerite's forward substitution gives `f^(i)_k = G_k f^(i)_{k-1}` with `G_k = diag(p_k)(I − W_{k-1} U_{k-1}ᵀ)` **independent of i** — so every column is the same linear recursion from a different start, and the sum collapses to one backward accumulation of a J×J matrix, `A_ii = 1/d_i + w_iᵀ R_i w_i`, `R_i = U_{i+1}U_{i+1}ᵀ/d_{i+1} + G_{i+2}ᵀ R_{i+1} G_{i+2}`, `w_i = p_{i+1} ⊙ W_i`, `R_{N-1} = 0`. No inverses, so it is as stable as the factorisation. `ampere.core.QuasisepGP` keeps its refusal (the numpy path's circumstances are unchanged), and the conformance row is widened from "it refuses" to **"it refuses by name *or* it agrees with `DenseGP` at `tolerances.cross_solver`"** — the stronger claim, and the one that catches a wrong `A_ii`, which would otherwise be finite, per-sample and plausible. |
| Variational inference, batching, and the float64 opt-out (W2.4 slice 2) | **Three surfaces landed 2026-09-07 on the torch track; none is a §4 change, and this row exists because each answers a question §4 or §5 asked and left open.** (1) **`VIEngine`** in `ampere.inference`, beside `NUTSEngine` and reached the same way — through `ampere.core.realise`, with the library imported lazily inside `run`, so the import-graph rule is untouched. §5's Phase 2 asks each track for "NUTS + VI"; this is torch's half, over pyro's SVI with `AutoNormal`/`AutoMultivariateNormal` guides. The one piece of real work is that SVI takes a *model* while a realisation is a bare density with no site structure (§10a fixes it that way): the model declares one standard-normal site over the whole unconstrained vector and adds a `pyro.factor` of `density(theta) - base.log_prob(theta)`, so the carrier cancels identically and the target is exactly the realisation's. A VI run emits **one chain** of i.i.d. guide draws — R-hat has nothing to say about draws with no Markov structure — with the guide family, step count, optimiser, learning rate and a thinned ELBO trace in the attrs, because the guide family *is* the approximation and an archived run must say which was made. `VARIATIONAL_LIBRARIES` is torch-only today; numpyro's `SVI` is the jax track's row to add. (2) **`log_likelihood_terms` is consumed**: `Engine.finish` gained an optional `log_likelihood_terms=` and builds each stored draw's `Evaluation` from the decomposition the realisation computed, taking only the prior from the declaration (cheap — it touches no model and no data). W2.13 shipped the member on both realisations and left the consumption to whichever slice 2 arrived first; the observable consequence is `engine_draws_recomputed = 0` on a realised run, where it was previously one full numpy model evaluation per stored draw of a quantity the backend had just computed. The fallback stays live and tested, because §10a makes the member optional. (3) **`GPSolver.configured(dtype=, device=)`** is `architecture.md` §5's per-run float64 opt-out and its device opt-in, in the only form that does not disturb the spec hash: a *copy* of the solver whose dtype and device are instance attributes shadowing the ClassVar policy, so `dataclasses.fields` still reports `jitter` alone (`Likelihood.to_spec` records fields; `results.md` §14 requires cross-backend spec-hash agreement) and `provenance_config()` — W2.13 fold-in 10 — reports what actually happened. `QuasisepGP.configured` **refuses** a float32 or non-CPU request by name rather than ignoring it, because celerite2's compiled kernels are float64 CPU and a `provenance_config` recording a policy the run did not follow is worse than a refusal. **`BATCHABLE` is now `True`** on torch's models, steps, kernels, noise models and `DenseGP`, and it means one specific thing: `LoweredProblem.log_prob_unconstrained_batched` evaluates a stack of free vectors in one call through `torch.func.vmap`. Making that true needed one change — the `math.isfinite` short circuit in `TorchParameterSpace.log_prior_tensor` became a `torch.where`, the last data-dependent Python branch on the density's path; the value is unchanged (a prior contribution is never `+inf`) and the cost is a handful of skipped scalar `log_prob` calls. `QuasisepGP` stays `False` and the batched call refuses by name, because a compiled extension reached through a `torch.autograd.Function` is not something `vmap` can rewrite and a loop under a batched name would be a false cost model. **CI stays CPU-only** (ruled with W2.13), so `device=` is tested at `"cpu"` only and the GPU smoke tests remain a later item. |
| jax quasiseparable GP library | **celerite2.jax, chosen by measurement; tinygp rejected on scaling** (decided by Fable on the W2.5 slice-2 branch 2026-09-07, discharging §6's deferred choice and its "verify at Phase 2 start which is better maintained" instruction; ratified by Peter 2026-09-08). Both were implemented as a full `QuasisepGP` against `ampere.core`'s **exact rank-2** Matérn-3/2 representation — never celerite2's approximate `Matern32Term`, which costs ~5e-3 in the log-likelihood at its default `eps` — and measured on one machine, CPU, float64. **Agreement**: both reach `scipy.stats.multivariate_normal` to 2e-11 on a log-likelihood of 4e3, both far inside `tolerances.cross_solver` (1e-6); against the dense Cholesky, tinygp holds 1e-13 to 1e-12 over 10 to 10⁴ length scales while celerite2 degrades from 1.6e-10 to 2.7e-8 over the same span (the cancellation `ampere.core`'s term docstring describes, bounded by the midpoint centring both paths inherit). Gradients agree to 1.7e-11 (tinygp) and 1.4e-10 (celerite2). **Scaling, which decided it**: celerite2 is linear at ~0.25 µs/point — 10⁴ points in 2.3 ms, 10⁶ in 0.26 s, value *and* gradient at 10⁶ in 0.75 s — while tinygp's jax-native `lax.scan` recursions are **quadratic in practice on XLA's CPU backend**: 3.0 µs/point at N=500, 6.4 at 10³, 10.3 at 2×10³, 18.1 at 4×10³, 40.0 at 8×10³, 98.0 at 1.6×10⁴, i.e. a doubling of the per-point cost with every doubling of N, and 0.47 s against 2.3 ms for one evaluation at 10⁴. A solver §2 calls "the scaling answer" cannot stop scaling at 10³ points. **Dependency cost**: none — celerite2 has been a *base* dependency since W2.3; tinygp would have added one package (needing only jax and equinox, both already in the extra). **Maintenance**: not decisive — celerite2 0.3.3 (2026-07-12) and tinygp 0.3.1 (2026-03-15) are both live. **Three costs, taken deliberately.** (1) `celerite2.jax` flips `jax_enable_x64` as an **import side effect**, which is precisely what `lowering.md` §10.2(a) forbids ampere from doing; the import is therefore **deferred** to first use, behind `require_x64`, which every caller has already passed at construction — so the flag is always already on, celerite2's `config.update` is a no-op, its warning never fires, and the side effect can never be what made ampere's guard pass. `_celerite2_jax()` asserts that invariant and a subprocess row checks it. (2) celerite2's primitives register **no `vmap` batching rule**, so `QuasisepGP` declares `BATCHABLE = False` alone in the backend and `declared_capabilities` aggregates that down to a problem that honestly says it cannot be batched — a refusal by name instead of a `NotImplementedError` about a primitive the user never wrote. (3) It exposes **no O(N) route to `diag((K + diag(σ²))⁻¹)`**, so `conditional_loo` stays refused on jax exactly as on the reference path (2026-09-05's deferral); tinygp *would* have given it — `L.inv().transpose() @ L.inv()` is a quasiseparable matrix whose `.diag.d` matched a dense inverse to 7e-15 — and that is the measured price of the factor of 200. **If the balance ever changes** (a tinygp release that lowers to a linear-time primitive, or an XLA:CPU fix for scan output aliasing), the strategy interface is what makes the swap a one-file change; re-run `tests/scaling` before believing it. |
| W2.5 slice 2's surface additions | **The realised path widened, VI added, batching and precision plumbed** (implemented by Fable on the W2.5 slice-2 branch 2026-09-07; ratified by Peter 2026-09-08 — post-freeze §4.5/§4.4 additions, hence this entry). Five things. (1) **The jax realisation lowers everything `ampere.core` implements**: every implemented family (`gaussian`, `student_t`, `cauchy`, `complex_gaussian`, `poisson`) through a new `ampere.backends.jax.families`, which transcribes the core's closed forms and refuses by name anything it does not hold; **censoring**, as the Tobit decomposition, with the CDFs `jax.scipy.stats` lacks written out (Student's t through `betainc`, agreeing with `scipy.stats.t.logcdf` to 1.1e-15); **latent GPs**, for the family that consumes one (`CONSUMES_LATENT_GP`, today Poisson — §4.4's singled-out case, which no gradient-free engine can run); and **both GP solves**, differentiable in the hyperparameters. Slice 1's four refusals (non-Gaussian family, censoring, latent, non-dense solver) are gone; what remains refused is a family `ampere.core` itself does not implement, an unimplemented analytic-GP combination, and another backend's model, step or solver. (2) **`VIEngine`** joins `ampere.inference` beside `NUTSEngine`, reached the same way — through `ampere.core.realise`, importing no backend — over numpyro's SVI with an `AutoNormal`/`AutoMultivariateNormal` guide. The bridge from a bare density to a numpyro *model* is one `ImproperUniform(real_vector)` site plus a `factor`, which adds exactly nothing to the log-density and leaves the guide fitted in the realisation's own coordinates. A run records `ampere_engine = "vi"`, the guide family, the step count, the optimiser and a (thinned, stride-recorded) **ELBO trace**, because draws from a fitted guide are not draws from the posterior and an archived file has to say so. jax only: pyro's route is one row in `VI_LIBRARIES` and one branch away, and is not written because nothing on this branch could exercise it. (3) **`Engine.finish` consumes `log_likelihood_terms`** (§10a sub-decision 1's "a driver uses it when it is there", which W2.13 recorded as owed): a driver that sampled through a realisation hands the per-dataset decomposition it already holds, `Engine.evaluations_from_terms` builds the per-draw records from it plus the model-free `log_prior`, and `engine_draws_recomputed` on a realised run drops from one-per-draw to **zero**. (4) **`BATCHABLE` becomes `True`** on every jax model, instrument step, kernel, noise model and `DenseGP` — it is measured, not asserted (`LoweredProblem.log_prob_unconstrained_batched` is `jax.vmap` of the same function and is compared against the loop), and `QuasisepGP` alone still says `False`. (5) **`architecture.md` §5's float32 opt-out and a `device=` opt-in** land on the jax solvers as `dataclasses.InitVar`s — so they stay out of the `Likelihood.to_spec` fields and therefore out of the cross-backend spec hash — and are reported from `provenance_config()` into `ampere_solver_config`, which is the home fold-in 10 created for them. Neither is auto-detected; the device is validated against `jax.devices()` and an absent platform is a refusal, never a silent fallback. GPU is untested and stays a later item; a non-CPU device today is caught by `declared_capabilities`' device rule against CPU models, which is the honest behaviour and says what a GPU item must widen. |
| `ampere.diagnostics` and the RHMF pre-fit family (W2.7) | **Deferred: the namespace and the `diagnostics` extra do not land, because the recorded adoptability assessment no longer holds on its maturity axis** (re-verified 2026-09-08 by the W2.7 agent, confirmed at the Fable review, discharging `diagnostics.md` §7's obligation that "the code-landing PR still carries the decision-log entry"; ratified by Peter 2026-09-08). W1.12 §2.2 made the namespace conditional on that assessment holding "at implementation time (licence, maturity, API — re-verify)", and two of the three axes do hold. **Licence holds**: MIT on the repository (GitHub reports SPDX `MIT`) with `LICENSE` shipped in the sdist — the one nit is that the wheel/sdist metadata carries only `License-File` and no `License`/`License-Expression`, so PyPI reports no licence and an automated audit sees nothing. **API holds, and is slightly better than recorded**: `Robusta(rank=…, robust=True, robust_scale=…)`, `fit(Y, W, max_iter=…)` and `synthesize()` are all present in 0.0.2, and so is `robust_weights(Y, W, state)` — precisely the IRLS weights §2.3's anomaly-score conversion wants — so the adapter is writable. **Maturity does not hold as recorded.** The 2026-09-01 assessment credited the project with "`tests/` and CI"; the tests exist (7 modules) but the CI does not — the sole GitHub workflow is `release.yml`, a publish-on-tag job that runs no tests on any push or pull request, and the project's own README still has "CI, automated tests, automated releases, and PyPI" unticked in its TODO list. The released artefact is **0.0.2, published 2025-11-24 — 9½ months stale**, two releases ever and no GitHub releases; the repository is alive (last push 2026-08-03) but its recent work is the paper's peer review rather than the library. The README's "Usage" and "Citation" sections are literally "TODO"; `Robusta.mse` raises `NotImplementedError` saying the author does not yet know which metric is appropriate; `fit()` prints to stdout unconditionally; `infer()` carries "NOTE: This implementation is a bit of a mess" and a docstring whose stated return does not match its code. A 0.0.x dependency with no release cadence and nothing running its tests is one that breaks silently, and the failure it would break in ampere is the screening step that exists to catch misspecification going undetected. **Two sequencing reasons point the same way**, and are the reason this is "not yet" rather than "never". §2.7 asks Phase 2 to *check* rather than assume away the JAX-pin compatibility between `ampere[jax]` and `ampere[diagnostics]`; W2.5 has since merged (its `jax` extra pins nothing yet), so there is no ampere JAX pin to check against, and shipping the extra now would fix a `jax>=0.6.0` floor (plus `equinox>=0.13.0`, `optax>=0.2.0`) in ampere's metadata and in `pixi.lock` ahead of the backend track whose job that choice is. And §2.5's Phase-2 obligation — validate `rank`/`robust_scale` against representative ampere SED/spectrum collections before promoting any default — still cannot be discharged, because no such collection exists in this repository; the namespace would ship with no default and no validated heuristic, which is the same reason W1.12 gave for shipping no defaults, applied one level up. **Consequences**: no `ampere/diagnostics/` package, no `diagnostics` extra in `pyproject.toml`, no pixi feature, `pixi.lock` untouched. `architecture.md` §3's namespace diagram still lists `diagnostics/` as the Phase-2 peer namespace and that stays correct — this defers the landing, not the placement, and §2.4's three reasons for a peer namespace are unaffected. **Nothing that landed forecloses it**: `ampere.core.AnomalyScore` and `ampere.results.plot_anomaly_score` already accept and render `provenance="rhmf_prefit"` (tested), so family A's adapter is an addition rather than a change. **Revisit trigger**: a `robusta-hmf` release at ≥ 0.1 with a test job running in CI, taken up after W2.5 has fixed ampere's own JAX pin. The post-fit families of the same item — B (separation-binned, permutation-calibrated residual whiteness, plus the chi-square posterior-predictive p-value) and C (`Likelihood.conditional` → `AnomalyScore` → the shared renderer) — landed in full and are not affected by this row; neither is a §4 change, since `likelihoods.md` §12 had already chosen the separation-binned statistic and named permutation calibration as its price. **Ratification note (Peter, 2026-09-08)**: early testing of RHMF is worth doing before its maturity gate is met, deferred to a later phase — recorded as a Phase 5 exploratory bullet in §5, not a Phase 3 item. |
| The latent-GP likelihood must see its kernel (W2.14) | **`GaussianProcessNoise.noise_params` applies `f = solver.latent_transform(kernel, coordinates, z, values)`, so `noise.latent` is `f` and not the whitened `z`** (ruled by Peter 2026-09-08; a §4.4 *clarification*, not a change — §4.4 always said the GP hyperparameters are ordinary fitted parameters, and `latent_parameter` always said the covariance enters through `f = L(θ) z`; nothing on the scoring path applied it). Found independently by W2.4 and W2.5 slice 2 and reproduced on master: on a Poisson + `GaussianProcessNoise` problem `problem.log_likelihood` was bit-identical for amplitude 0.5/5/50 and length scale 0.1/1/100, so a latent fit sampled both hyperparameters against a flat likelihood and reported nothing wrong. The transform belongs in `noise_params` because that is the only place holding the solver, the kernel, the **retained** coordinates and the resolved hyperparameters at once — no family changes, and every family already reads `noise.latent` as `f`. Consequences taken here: `likelihoods.md` §17 limitation 6 is amended (the remaining limitation is the *inference*, not the mathematics); the torch backend's blanket refusal of the latent path and the test that asserted the flatness are gone, and both backends' realised `_latent` now apply their own native `latent_transform` (torch `latent_transform_native`, jax `latent_transform_jax` — added for this, differentiable in the hyperparameters, because the contract-surface `latent_transform` returns numpy and would cut the graph); the conformance battery gains a latent-GP shape (`LATENT_GP`, a Poisson family over a dense GP) whose realised density is held to the corrected numpy oracle per backend. |
| Results completion: the plotting surface, the pointwise group, the training-set writer (W2.8) | **`results.md`'s three declared-but-unimplemented surfaces are code, and the contract text that said otherwise is amended in place** (implemented 2026-09-08; ratified by Peter 2026-09-08 — three §4 amendments, hence this row, and none of them changes a shape, a name or an obligation the contract fixed). **(1) §8's closing paragraph and §12's matching row** said "every plotting function is a declared signature that raises `NotImplementedError`". After W2.7 drew three and W2.8 drew the other three — `plot_corner`, `plot_trace`, `plot_posterior_predictive` — none does; what they raise on a non-run is `ResultsError`. The surface itself was not moved to make the drawing possible, which is the whole point of having fixed it before anybody drew. Judgement calls inside the surface, recorded because they are the implementation's and not the contract's: `MAX_CORNER_VARIABLES = 20` and `MAX_TRACE_VARIABLES = 40` are where §8's "refused loudly rather than attempted" bites, both overridable by a named `max_variables=`; `plot_corner` **excludes prior-rejected draws** (θ is stored for them, §5, and a point of zero prior mass is not a posterior sample); `plot_trace` renders them as NaN gaps and draws `lp` beside the parameters; and `plot_posterior_predictive`'s default discrepancy is a σ-standardised sum of squares about the **replicate mean**, a function of `y` alone, because the group stores replicates and not the per-draw prediction — a caller who wants `T(y, θ)` passes one. **(2) §6's per-observation group is emitted**, by `add_pointwise_log_likelihood(tree, problem)` and by nothing else: the never-by-default rule is unchanged, and having a computation available is exactly when that rule stops enforcing itself. Two additions beyond the reserved constant, both recorded in §6: the `ampere_decomposition` attribute is written **per variable** as well as on the group, and reads `"mixed"` on a joint fit whose datasets disagree, because reporting one of two decompositions would be reporting a falsehood about the other; and `pointwise_as_log_likelihood(tree)` is a one-line bridge to `arviz.loo`, which reads the group *named* `log_likelihood` where ampere deliberately keeps the per-dataset split (§15 R6). Under a solver that does not implement `conditional_loo` — `QuasisepGP` and its jax twin, deferred at W2.3 — the emission **refuses by name and does not fall back to a dense solve**: at the sizes the quasiseparable path exists for that substitution is a different program, not a slower answer. **(3) §11 layer 2's training-set writer exists** (`ampere.results.training`: `write_training_set` / `append_training_set` / `read_training_set`), and its table gained two rows for what `serialisation_review.md` §4 asked the writer to carry and the table, written first, had nowhere to put — an `observations` group (one subgroup per dataset label, when the budget drew any) and the whole `Failure` record in `sample_stats` rather than the reason alone. The review's two named losses are closed: θ keeps its dtype, which changes what `training_pair_to_dict`/`model_result_to_dict` **write** and therefore bumps **`CONTAINER_SCHEMA_VERSION` to 2** under the review's own §5 rule (the reader still accepts the old bare-list form, since a record that stated less than it could is not a corrupt one); and `training_pair_from_dict` completes the in-memory round trip. Append is real — the `sample` dimension grows and the file's `ampere_spec_hash` is compared before anything is written, so §7's stale-artefact trap refuses rather than leaving one file whose halves came from two simulators — and it is read-concatenate-rewrite, `O(existing + new)` per call, which is the one limitation created in place of the two closed and is recorded as §13.9's extension point. **No version constant moved except that one.** `PROVENANCE_SCHEMA_VERSION` stays **5**: no existing mapping changed meaning, and adding groups (`pointwise_log_likelihood`, `posterior_predictive`) or a new file format is explicitly not a bump under §9's rule — every `ampere_problem_hash` in the repository is unchanged, which the conformance suite checks. `TRAINING_SET_SCHEMA_VERSION` starts at 1 because the format starts here. *(Amended W3.10, 2026-09-09, on Peter's ruling of 2026-09-08)*: `MAX_CORNER_VARIABLES`/`MAX_TRACE_VARIABLES` no longer refuse outright above the cap — `plot_corner`/`plot_trace` **page** into a `list[Figure]`, one page per cap's worth, in merged-name order, with an array-valued block split only when it cannot fit a page whole; a loud `ResultsWarning` names the page count, the cap and the `var_names=` override; `paginate=False` restores this row's original refusal unchanged; a call that fits within the cap still returns a single `Figure`. |
| Benchmark harness (§6's deferred choice), and the docs gate's missing `-W` (W2.11) | **pytest-benchmark, not asv** (decided 2026-09-08 on W2.11; ratified by Peter 2026-09-08). §6 deferred the choice to "the phase that needs them", and this is it — §5's Phase 2 line asks for a "benchmark suite … with results tracked as CI artefacts", which is the requirement that decides it. **Environments decide it first.** asv owns its environments: it builds them itself from a matrix in `asv.conf.json`, out of a wheel it builds from the checkout. Ampere's environments are pixi's, and the three the benchmarks care about (`dev`, `torch`, `jax`) differ in exactly the ways that move the numbers — one of them now pins torch to a CPU-only index — while AGENTS.md's working agreement is that every CI command is a pixi task. An asv matrix would be a second, divergent definition of those three environments, and the first time it drifted the benchmark numbers would silently stop describing what CI runs. pytest-benchmark is a pytest plugin: it inherits the environment it is invoked in, and its rows skip or run on the same `importorskip` as every other suite. **What "tracked" means here decides the rest.** asv's real value is its history database and web front end, which want either a committed `.asv/` results tree or a separate results repository; ground rule 7 forbids run outputs in git, so that history would need infrastructure outside the repository, for a project whose stated requirement today is that the numbers are attached to the PR that changed them. `--benchmark-json` uploaded as a workflow artefact is exactly that, and carries more than a reviewer needs (per-round statistics, machine, interpreter, commit). **Reuse** closes it: these rows need ampere's own fixtures, skip logic and jax `configure_x64`, which in pytest is code already written and in asv would be rewritten in asv's class-based API — so the O(N) claims would be measured by code the conformance suite has never run. Cost of the choice, stated: no automatic regression detection against history and no bisection. Revisit if this project ever wants a benchmark *server* rather than benchmark artefacts; nothing forecloses it, since the measured quantities are ordinary functions in `tests/benchmarks/test_gp_solvers.py`. §6's bullet is struck accordingly. **Ridden along, because it is a documented deviation rather than a decision**: §5 puts "Sphinx docs build check with warnings-as-errors" in Phase 1, and W2.11 lands the docs job **without** `-W`. `pixi run docs` was red on master for four reasons; three are fixed at source (the `dev` pixi feature now declares `pandoc` and `ipykernel`; `docs/source/conf.py` sets `nbsphinx_execute = 'never'` and records, per notebook, why neither legacy notebook can run — one reads a Spitzer/CASSIS FITS file that is not and cannot be in the repository, the other needs torch+sbi and trains on 10 000 simulations; the `automodule` entry for the non-existent `ampere.infer.ptemceesearch` is gone). The fourth is 21 residual warnings, all of them legacy docstrings and legacy `.rst` under frozen paths, where editing to silence a warning would be touching frozen code for a cosmetic reason. So warnings-as-errors waits for the Phase 6 docs rebuild that retires those pages, and what the gate asserts today is the stronger half of the intent anyway: the documentation builds, from a clean clone, in a declared environment, on every PR. |
| `QuasisepGP.conditional_loo` on the jax path (W2.5 slice 3) | **Supplied, in O(N), by the same recursion the torch path runs; the refusal is lifted here too** (implemented 2026-09-08 on the W2.5 slice-3 branch; ratified by Peter 2026-09-08). This is the open path W2.5 slice 2 recorded rather than a new ruling: that row closed with "a jax `conditional_loo` becomes cheap the day this backend calls celerite2's kernels directly, as torch does", and it does now. The 2026-09-05 deferral's reason was a **coupling**, never the mathematics — the leave-one-out terms need `A_ii` for `A = (K + diag(σ²))⁻¹`, celerite2's high-level surface has no O(N) route to it, and the route that exists "reimplements celerite2's internal factorisation convention". Slice 2 reached celerite2 through its `GaussianProcess` wrapper, so that coupling would have been new *and* would have meant reading private `_d`/`_W`. It turns out not to be needed: **`celerite2.jax.ops` is public** and is the same set of compiled kernels `celerite2.backprop` exposes to the torch backend, so `ops.factor(t, c, a, U, V)` hands over `d` and `W` beside the term's own `c` and `U`, and `ops.solve_lower`/`solve_upper` are what `GaussianProcess._do_solve` itself calls. The recursion is W2.4 slice 2's, unchanged and not rederived: `A_ii = 1/d_i + w_iᵀ R_i w_i`, `R_i = U_{i+1}U_{i+1}ᵀ/d_{i+1} + G_{i+2}ᵀ R_{i+1} G_{i+2}`, `w_i = p_{i+1} ⊙ W_i`, `R_{N-1} = 0`. One thing changed in the transcription: torch runs a Python loop, which on jax would stage N copies of the body into the jaxpr, so this is a single `jax.lax.scan` with `reverse=True` — the terms land at their own indices and no `.at[].set()` bookkeeping is needed. Measured against `DenseGP`'s Cholesky on a 250-point problem: the precision diagonal to 1.2e-13 relative, the terms to 1.0e-12 absolute. **`ampere.core.QuasisepGP` keeps its refusal** — the numpy path's circumstances are genuinely unchanged, since it composes celerite2's `GaussianProcess` and does not call the kernels itself — so the conformance row still states both branches ("refuses by name *or* agrees with `DenseGP` at `tolerances.cross_solver`") and what changes is only which branch jax takes. tinygp is still not needed: it was rejected on a factor of 200 in the marginal likelihood, and nothing here touches that. |
| The gradient-free fast path (W2.5 slice 3) | **`Engine`'s evaluation cache scores through the backend's realisation when there is one, `use_realisation=True` by default, and the numpy path stays the oracle and the fallback** (implemented 2026-09-08 on the W2.5 slice-3 branch; ratified by Peter 2026-09-08). An engine-side change under `inference.md` §10a's existing text — no contract shape, name or obligation moves — recorded here because it changes what a default gradient-free run *does*. The problem was measured at W2.10: the jax **contract** path costs ~28 ms flat per `log_prob` whatever the problem's size, because `FittingProblem.evaluate` is numpy and a jax model evaluated an operation at a time pays jax's dispatch on each one — so emcee, dynesty and zeus were *slower* on the backend built for speed. §10a's optional `log_likelihood_terms` already existed and W2.4 slice 2 already consumed it for the per-draw record; slice 3 lets `_EvaluationCache` **score** through it as well: the prior from `problem.parameters.lnprior` (it touches no model and no data, and §10a's mandatory surface has no prior/likelihood split to ask for), the per-dataset terms from the realisation, the same `Evaluation` out. **Nothing in `ampere.inference` imports a backend** — `ampere.core.realise` dispatches on `problem.backend` — so §10's import rule is untouched and a backend with no registered realisation keeps the contract path silently. Four guards, all deliberate: `realise`'s own reference-point agreement check gates the attachment; the dataset labels are checked once at attachment (which also pays the backend's compilation before sampling rather than on the first proposal); any exception at evaluation time falls back to `problem.evaluate` for that θ; and a problem declaring **`strict=True` keeps the contract path**, because a realised density can neither raise nor build a `Failure` (§10a sub-decision 2) and `strict` is precisely the declaration that reasons matter more than speed. **The cost, stated rather than hidden**: on this path a proposal the model cannot score is a bare `-inf` with no reason, so `failure_counts` and `ampere_failure_summary` stay empty; `use_realisation=False` restores them. Provenance gains `engine_realised_evaluations`, and `ampere_realised` is now answered by the cache rather than fixed at 0 for the gradient-free drivers — which is what that attribute has always meant ("whether the draws were scored through the backend's realisation"). **Measured** (400-point flexible likelihood, one `Engine.log_prob`, `tests/benchmarks/test_engine_fast_path.py`): jax `QuasisepGP` 27.3 ms → 0.62 ms (**44x**), jax `DenseGP` 7.9 → 1.8 ms (4.4x), torch `QuasisepGP` 1.07 → 1.08 ms (1.0x), torch `DenseGP` 6.0 → 5.2 ms (1.2x). The jax number needed one change in the backend as well: `LoweredProblem.log_likelihood_terms` is now `jax.jit`-compiled on first use and cached per lowered problem, because it was the one member of the realised surface a caller might invoke many times *without* a surrounding transformation (numpyro jits `potential()` itself) and eager dispatch was the whole of the 28 ms. **Open for Peter**: on torch the fast path is a wash, so it buys nothing there and still costs the failure reasons; leaving it on by default is what the item asked for and keeps one rule for both backends, but a torch-specific default of `False` would also be defensible. `NUTSEngine` and `VIEngine` pass `use_realisation=False`: they already hold a realisation and the base class would otherwise lower the same problem a second time for the start-point search alone. |
| Phase 2 documentation sweep: the frozen specs' stale claims, annotated in place (W2.15) | **Twenty-one annotations across seven frozen documents, additive and marked *Amended W2.15*; no contract semantics change, and the conformance suite is untouched** (swept 2026-09-08 on the W2.15 branch; ratified by Peter 2026-09-08). Every "landed"/"amended"/"Phase 2 will" claim in `architecture.md`, `lowering.md` and the seven `contracts/*.md` was checked against the merged code rather than against the work-item reports. `contracts/parameters.md` and `contracts/inference.md` came back **clean** — §10a's realisation surface, the four capability flags and the engine-side claims all check out line by line. The rest: **`architecture.md`** gains one subsection after §3 recording the tree and the extras table as built — `core/` has `lowering.py`, `realisation.py` and `rng.py` and does *not* have `astropy_compat.py` (Phase 4); `backends/reference/` is not "models and transformations only" since W2.13 fold-in 7 made every backend declare its own noise models and GP solvers; `inference/` ships five drivers, not three, and neither the SBI layer nor the optimisers the comment names; `diagnostics/` is the ruled deferral of 2026-09-08, placement intact; the `torch`/`jax` rows' "(GP solver library — deferred choice)" is spent (celerite2 on all three backends, so neither extra gained a package); `sbi` installs torch but not pyro-ppl, so it unlocks importing the torch backend and not its engines; `all` is the feature extras and deliberately not `dev`; and the table predates pixi, so the environments, the `netcdf` feature and W2.11's CPU-wheel index pin are recorded beside it. §9's "GP solver library per backend" open item is marked resolved. **`lowering.md`**: §0's jax row said paramax and the answer is an `eqx.partition` filter spec; §6.2's and §12's paramax interop caveat was conditional on a GPJax-founded GP library and lapsed unspent when the choice went to `celerite2.jax`; §3.2's "`loc` must be 0" rows are what the target libraries take, and both backends compose the affine map away per §3.3 rule 1 for `lognorm`, `gamma` and `beta` as well as `expon`. **`contracts/likelihoods.md`**: the namespace's dependency list now includes celerite2 (base, lazily imported inside `QuasisepGP`), and the `KernelSpec`→celerite2-term translation landed in `ampere.core` at W2.3 rather than in a backend. **`contracts/results.md`** carries five: arviz is not optional (W2.2 executed §15 R1 — base dependency, extra deleted, lazy import kept for cost); the three derived surfaces landed at W2.7/W2.8 **on the numpy path in `ampere.results`**, not behind a backend; limitation 13.2 (no engine drivers) is closed; limitation 13.4 (`warmup_*`) stands but its stated reason has expired, since samplers with warm-up now exist and choose not to store it; and §15's archival "`PROVENANCE_SCHEMA_VERSION` is now 2" is 5. **`contracts/results_schema.md`**: only the `ampere.results` producer of `AnomalyScore` exists, and W2.8 discharged the training-set writer obligation. **`contracts/transformations.md`**: the capability-flag row is three flags and there are four (`BACKEND`, W2.12). **`contracts/diagnostics.md`**: §2.2's adoptability verdict was re-verified at W2.7 and family A is deferred on the maturity axis; §2.6's and §7's Phase-2 checklists are half-done (B and C landed in `ampere.results`, A did not); limitation 9.3's irregular-coordinate strategy is closed by `likelihoods.md` §12's separation-binned choice, implemented hand-rolled with permutation calibration and no statsmodels. **Two findings recorded rather than fixed, both out of this item's scope**: the name `ampere` on PyPI belongs to an unrelated battery-modelling package, so every `pip install ampere[...]` instruction in the package's docstrings and in `OptionalDependencyError`'s runtime message names the wrong distribution until ampere is published — the README and the install page now warn, but the runtime string is a code change and was left; and `SAMPLE_STATS_GROUP` is defined twice, identically, in `ampere/results/emission.py` and `ampere/results/training.py`. |
| Non-native parts in a native problem (W2.4 slice 3's carried finding) | **Refused by default; accepted under an explicit opt-in when no gradient is needed; always refused where a gradient is required** (ruled by Peter 2026-09-08 on Fable's recommendation; **landed W3.8** 2026-09-09, which supplied the §4.5 text — `inference.md` §4.5, "Non-native parts in a native problem"). The cause was one level below fold-in 7's: a solver builds its covariance by *calling* `Kernel.matrix`, so an `ampere.core` kernel inside a native solver hands back numpy, the solver converts it, and the conversion detaches the graph — the problem declared `differentiable=True` while its amplitude and length scale got no gradient at all. The kernel therefore **joins `Likelihood.capability_parts`** beside the noise model and the solver, and `Kernel` declares the four capability flags; fold-in 7's "the kernel is consumed *by* the solver" is withdrawn. (1) **Default**: the backend-disagreement refusal, naming the offending pieces by their **fully qualified** class — the bare name does not discriminate, since a backend's kernel deliberately shares it so that `Likelihood.to_spec` and `results.md` §14's cross-backend spec hash agree — with both halves of the remedy. (2) **Opt-in**: `FittingProblem(..., allow_foreign_parts=True)` (the working name kept; it sits beside `strict` and like `strict` has no global form), for a piece expressible only in Python. Under it `differentiable` is `False` whatever the parts declare, and the run records `ampere_foreign_parts` with the located names (`"sed: ampere.core.likelihood.Matern32"`) — conditional, and **no `PROVENANCE_SCHEMA_VERSION` bump**, because the attribute can only appear on a composition no earlier ampere could build (accepted at review; W3.12's bump to 6 will fold it in). "The problem's backend" is resolved **structurally**, with no part list and no family assumption: `"reference"` is the name a piece inherits by silence, so exactly one non-`"reference"` backend among the parts makes that one the problem's and the silent ones foreign, while two native backends have no split and stay refused either way — a part stops being foreign the day it declares the backend, which keeps the pathway open for later core `sample` implementations and native twins. `FittingProblem.foreign_parts` reads the *parts*, not `problem.backend`, so `capabilities=Capabilities(backend=...)` keeps its documented meaning while still failing to hide a genuine mixture. (3) **Gradient routes**: `ampere.core.foreign_parts_refusal`, called from `realise` (before dispatch) and from both backends' `LoweredProblem`, and `_refuse_foreign_parts` in `ampere.inference.engine`, called by `NUTSEngine` and `VIEngine` ahead of their other checks — refusing by name whatever `allow_foreign_parts` and `strict` say. `realise`'s refusal *is* the gradient-free routing: `_EvaluationCache` already falls back on a `LoweringError`, so the fast path falls back with no engine change and `ampere_realised` is 0. **Measured, not assumed** (400-point flexible likelihood, one `Engine.log_prob`, medians of two runs): torch 10.7–12.6 ms contract, 7.8–9.1 ms callback (detach→numpy→tensor), 7.3–7.5 ms all-native; jax 10.8–13.6 ms contract, 9.3–9.8 ms callback (`jax.pure_callback`), 2.0–2.2 ms all-native. The callback buys 1.2–1.6× where the whole fast path is worth 1.4–1.7× (torch) and 1.2–1.4× of an available 5–7× (jax), at the price of a second silently non-differentiable path — so it is **not implemented**. Conformance: `tests/conformance/test_foreign_parts.py`, a row per backend for each of the three parts; the mirror fixture gained kernels declared as its own. |
| Batched simulation: execution, external simulators, native sampling (Phase 3, W3.1) | **Three rulings by Peter 2026-09-08 on the first draft of W3.1.** (1) `vmap` is a single-device route and must not be the design: `simulate_many` is built on an order-preserving, partition-independent **executor protocol** (serial default; a process pool for external simulators; user-supplied mappers satisfying `concurrent.futures.Executor`'s `map`/`submit` shape — dask, ray, `MPIPoolExecutor` — for many machines), **chunked** so a budget never has to fit in memory at once and a single simulation is never asked to share a device with another; a single simulation larger than one device is the model's own concern (a model-parallel `__call__`), not something the framework partitions. (2) Slow **external simulators** (compiled Fortran/C/C++/Rust routines behind a Python call, composed as black-box models on the reference backend) are the canonical SBI case and are first class: process-pool execution with per-simulation timeouts and crash capture as flagged failures, and a reference-backend meaning for `BATCHABLE` (a model that takes a table of θ in one call). (3) **Every backend supports observation sampling natively**: `LikelihoodFamily.sample` gains torch and jax twins over the libraries' own distributions, the numpy path remaining the oracle and compared distributionally since RNG streams differ by backend. W3.1 slice 1 lands (1)–(2) on the numpy path; slice 2 lands per-chunk `vmap`, (3), and — Peter's addendum of 2026-09-09 on approving the revised text — API-level hooks for sharding a chunk across devices (`jax.pmap`/`shard_map`, torch's distributed equivalents), smoke-tested without hardware. The §4.5 text lands with the slices, which amend this row. **Slice 1 landed 2026-09-09**, as `inference.md` §13's two new subsections (*batched form*, *execution*) with limitation 17.5 closed. `FittingProblem.simulate_many(count, *, values, observe, rng, stream, executor, chunk_size, as_chunks, on_chunk)` -> `SimulationBatch`, whose contract is an **equality**: `simulate_many(n, stream=s)` returns exactly what `[simulate(rng=child) for child in rng(s).spawn(n)]` returns, in order, whatever executor ran it and however it was chunked. The per-draw generators are `numpy.random.Generator.spawn` children derived by index (spawning is cumulative, so a chunk at a time equals all at once), and theta is drawn in the parent so a draw killed by a timeout still records the theta it died on. `ampere/core/simulate.py` holds the `Executor` protocol (`map(fn, items)`, `concurrent.futures.Executor`'s shape) with `SerialExecutor`/`ThreadExecutor`/`ProcessExecutor`; `FailureReason` gains `EXECUTION_FAILED` for a draw the executor lost (timeout, dead worker), distinct from a simulator that reported its own failure, and failure recording is replayed on the parent in draw order so a pooled budget's counts are right. Stacks are `ContainerBatch`, **not** `FunctionSamples` with a sample axis — `results_schema.md` fixes a container's value shape to the one its coordinates imply, so the item's wording would have contradicted a frozen invariant. `Model.evaluate_batch`/`call_batch` give `BATCHABLE` its reference-backend meaning, used under the serial executor only. `results.md` §11/13.9 amended: the training-set writers take the chunk iterator, which moves the cost of a very large budget onto the `O(existing + new)` append. Enabling all of it needed a `copyreg` reduction for `mappingproxy` (`ampere/core/_pickling.py`) — nothing in `ampere.core` could be pickled before, so no problem could reach a worker. **Slice 2 landed 2026-09-09**, as §13's *sampling on a backend* subsection plus the native paragraphs of *batched form* and *execution*, and `likelihoods.md` §3's native-twin note. A realisation may now offer two optional members, looked up by presence exactly as `log_likelihood_terms` is: `simulate_batched(theta, *, chunk_size, sharder)` -> `BatchedPrediction` (every model's every channel and every dataset's full instrument-transformed prediction, through `torch.func.vmap`/`jax.vmap`, **per chunk with the chunks looped**, on the *constrained* free vector) and `sample_observations(theta, predicted, seeds)`. `simulate_many` gains `native=`/`sharder=`/`context=`; `_NativeBatch` (`ampere/core/dataset.py`) is the single place the fallback lives, taking its container templates and its one-point agreement check from one contract-path evaluation at the reference values — sound because a predicted container's axes and effective mask are evaluation-invariant by contract. **The guarantee is two-graded and stated as such**: `native=False` reproduces the loop bitwise, the native prediction agrees to `tolerances.cross_backend` (with partition independence unweakened — `chunk_size` 1, 7 and whole agree), and native observations are a draw from the same distribution rather than the same draw, since `jax.random` and `torch.Generator` are not numpy's stream. Which path ran is recorded rather than inferred: `SimulationBatch.provenance` carries `simulate_batched`, `sample_backend`, `evaluate_batch` and the reserved `simulation_context`, merged conservatively across chunks so a mixed budget never reads as a native one, and written to the training set's root attributes. **Peter's ruling (3) has a ceiling, and it is the load-bearing half**: a backend samples exactly what `ampere.core` samples, which today is `GaussianFamily.sample` alone — a native twin for `student_t`, `cauchy`, `complex_gaussian` or `poisson` would be a backend guessing an observation process the contract declines to guess, so those refuse on every backend, with the refusal handed back to the numpy path so the text matches exactly. Native draws are checked *distributionally* (mean and covariance of 4 000 draws against `K + diag(sigma^2)` at the estimator's own standard error, plus the numpy row's two mutation-tested assertions). **Sharding** (Peter, 2026-09-09) is `ChunkSharder` — `devices()` and `shard(fn, stacked)`, a *value* contract — with `SingleDeviceSharder` on both backends, `MeshSharder` (`jax.pmap`, padded chunks) and `DistributedSharder` (torch, strided ranks over an **existing** process group, refused by name without one; ampere launches none); CPU CI exercises the degenerate case and `tests/gpu` holds the rest. **Executor carry-overs**: `ProcessExecutor` reuses its pool across `map` calls (replaced only when a worker died or a draw expired) and defaults to **`forkserver`** on POSIX (`fork` still available via `mp_context=`) — forking a threaded runtime may deadlock and Python 3.14 moves the platform default anyway. That default exposed two real defects, both fixed here: a pickled `FittingProblem` could not be *used*, because `Instrument.freeze` fingerprints its steps with `id()` and pickle does not preserve object identity (`Dataset.__setstate__` re-freezes from the steps; the `id()` fingerprint itself is carried as a finding), and `timeout` measured from submission included pool start-up (workers are now started with an undeadlined no-op first). Also landed: the two W3.0 carry-overs — torch's GP-marginal zero-uncertainty refusal, and jax's `GaussianProcessNoise(scale=, jitter=)` parity, with a conformance row. |
| The SBI engine's shape (W3.2) | **`SBIEngine` is an ordinary `ampere.inference` engine over `simulate_many`, and six implementation decisions are recorded here because each answers a question the item left open** (implemented 2026-09-09 on the W3.2 branch; accepted at the Fable review; Peter to ratify at leisure — none is a §4 contract change). (1) The prior sbi trains against is the **unconstrained** one: `sample` draws through `FittingProblem.sample_prior` and maps through `unconstrain`, `log_prob` is `ParameterSet.lnprior_unconstrained` (prior plus Jacobian), so bounded priors need no restriction and the density estimator sees ℝⁿ; draws are mapped back through `constrain` before the numpy path scores them. (2) The estimator's own per-draw log-density lives in **`sample_stats` as `ampere_sbi_log_prob`**, added after `emit()` rather than through a new hook, because an importance weight is one value per stored draw in draw order and belongs in the `(chain, draw)`-shaped group; `ampere_sbi_log_prob_kind` says normalised (NPE) or unnormalised (NLE/NRE). `emit()` gaining a first-class extra-sample-stats hook is deferred to W3.6, which may need it — that would be a `results.md` §4 amendment and its own row. (3) Every sbi posterior is queried through `potential()`, the one surface all of them share (`log_prob` is deprecated on the two that answer only up to the evidence). (4) The `"flat"` summary is one function, `_summary_of`, refusing complex containers by name — the real/imaginary encoding is W3.3's to fix; W3.3 also adds the `layout=` argument rather than W3.2 pre-empting its vocabulary. (5) `run()` takes `training=` and `posterior_options=` pass-throughs to sbi's `train` and the posterior's `sample`, with no ampere defaults interposed, because NLE/NRE's MCMC sampling is not affordable in a per-PR gate otherwise. (6) Each round simulates on its own sub-stream (`sbi.simulate.<round>`), so changing round 1's budget cannot silently change round 2's draws; and the legacy embedding vocabulary's fall-through on an unknown dict `type` (which handed the dict itself to sbi as a module) is a refusal by name — a bug fix, not parity. **Ridden along at the merge**: the `sbi` pixi environment gains the `torch` feature, because the `sbi` extra installs torch without pyro-ppl and the torch NUTS/VI/M2 rows do not skip without it, so the environment's gate had been red since W2.11 created it; and the eight `*_version` provenance attrs in `_nuts.py`/`_vi.py` are `str()`-wrapped, since `TorchVersion` is a `str` subclass h5netcdf cannot write and a torch NUTS/VI run could not be saved. The `context=` slot exists and accepts only `None` (`ampere_sbi_context = "none"`), per the reserved hook. Largest budget run: 20 000 simulations of the joint problem in 152 s end to end; a further 5 000-draw append took 0.9 s, so `results.md` 13.9's append stays deferred. **W3.5 (2026-09-09)**: `SBIEngine(cache=ArtefactStore)` wraps the round loop — a hit skips simulation and training and the run says so (`ampere_sbi_cache_hit`, the key digest; `ampere_sbi_simulations = 0`); the key is `ampere.results.artefacts.artefact_key` with the **encoding hash** as its layout ingredient, plus a `model_hash` beyond the item text so a kernel/solver/family swap with identical parameters is a miss. The lazily-built prior and embedding-wrapper classes gained `__reduce__` hooks so a trained posterior pickles. |
| Encoding contract (post-freeze §4 addition, drafted for W3.3) | **`docs/design/contracts/encoding.md` is drafted 2026-09-09 (Fable, at Peter's request) as the packing contract for embedding networks: rows are samples, column groups in fixed order (dataset index; standardised coordinates; Fourier coordinate features; whitened value `y/σ`; `asinh`-scaled value; `log σ`; one mask column; per-set features `log N`, `has_sigma`, `is_complex`; a reserved zero-width context group), a frozen hashed `Layout` computed from the observed data alone (never from simulations), padding to a row cap with refusal beyond it, sbi's own z-scoring switched off for set layouts, and `unpack` as the one way a network reads the tensor so that every embedding is a wrapper over the unpacked view.** The reasons: sbi's interface is a module over one tensor, so the packing is the only contract there can be; σ and coordinates as columns are what amortisation over noise and over sampling need (horizon notes §4–5); statistics from the observation make the layout a function of the problem, fixed before any simulation and identical at inference. W3.3 lands the code, binds the text, and amends this row; the document is frozen at that merge. **W3.3 landed the code and bound the text 2026-09-09** (`ampere/core/encoding.py`; the wrappers, `layout=` and the attrs in `ampere/inference/_sbi.py`; `tests/core/test_encoding.py`), and the document is frozen at this merge (`1ca558a`). Seven sentences were **amended in place** where the code proved them wrong, each marked *Amended W3.3*. Two change what a reader would otherwise get wrong: (a) a dataset's `rows` is its **total** sample count with `valid` beside it, since a masked sample must be a row with `mask = 0`; (b) **the mask is frozen into the layout** as each dataset's excluded indices, because `Dataset.effective_mask` resolves lazily at the first prediction, so encoding the observation before simulating and the draws after could give two packings under one hash. The rest: the class is `EncodingLayout` (`ampere.core.Layout` already names the container schema's enum) and `EncodingError` subclasses `ContractError`; `Encoded` is always `(count, rows, columns)`; an absent axis contributes zero Fourier features; a masked sample's value is left as held unless non-finite (then zeroed), while a retained non-finite value or non-positive σ is refused; padding sits at the tail and per-dataset row counts are in the layout, so per-dataset capacity stays reserved. One drafted refusal did not land: ignored non-numeric extra coordinates are not recorded, because putting free-form per-sample labels in the hash would make renaming a filter a retrain. **§7 gained two sbi 0.27 facts a later encoder must not re-derive**: `PermutationInvariantEmbedding` pools with a sum by default (built with `"mean"` here) and computes its valid-row count from the first batch element only, correct under this contract solely because the mask is a layout property identical across a batch; and `TransformerEmbedding.forward` discards `attention_mask` unless `is_causal` is true and aggregates by the last token, so the wrapper zeroes masked tokens after its own projection and passes the mask anyway. `_summary_of` survives as a thin adapter over `encode`, so `"flat"` has one implementation. |
| Family D lands: the `calibration` group (W3.6) | **`results.md` §4 gains the `calibration` group; §7 and §8 amended** (2026-09-09, W3.6; accepted at the Fable review; Peter to ratify at leisure). `diagnostics.md` §11's family D lands in `ampere.results.calibration`, with the `sbi` fast path as `SBIEngine.calibrate`. A run may now carry a fourth non-default group, `calibration`, holding `ranks (simulation, parameter)`, `coverage (level, parameter)`, `ks_pvalue (parameter)` and the route's own extras (`c2st_ranks`, TARP's curve on its own `tarp_level` dimension), written by `ampere.results.attach_calibration` and by nothing else; the coverage curve is *derived from the ranks* so the two halves of one figure cannot disagree. Like `posterior_predictive`, `residuals` and `gp_localisation` it is **not stored by default** — it costs a fresh simulation batch at best and a full fit per simulation at worst — so §7's rule extends to it unchanged; its randomness comes from the named sub-stream `"calibration"`, distinct from `"simulate"` (`lowering.md` §9.2), so adding a calibration check cannot change what an SBI budget simulated. The SBI route encodes the calibration batch with **the run's own `EncodingLayout`** (a layout rebuilt from simulated containers would standardise differently and report a different network as calibrated). §8's plotting surface grows from six functions to eight (`plot_sbc_ranks`, `plot_coverage`), both taking the run or the bare group. `PROVENANCE_SCHEMA_VERSION` unchanged (no existing mapping changed meaning); the group carries `ampere_calibration_schema_version = 1`. **W3.2 decision 2's deferred `emit()` extra-sample-stats hook is declined**: the group is attached after emission exactly as `ampere_sbi_log_prob` is. The any-engine route's `engine_factory` is duck-typed on shape (an engine, or an emitted run) because `ampere.inference` imports `ampere.results` and the dependency stays one-way. **The WStat coverage study** (W2.9's request) runs by default under the *faint* prior `log-U(0.5, 3)` with `--broad-prior` as the comparison, because SBC averages over the prior it draws from and the example's broad prior is mostly bright sources where profiling is harmless — a judgement for Peter to confirm (recorded on the item). Conformance unaffected; `dev` 1900/200 and `sbi` 2615/56 green on the branch. |
| TMNRE lands: the `marginals` group and the truncated prior (W3.4) | **`results.md` §4 gains the `marginals` group, §7 gains its exception; `inference.md` §13 gains "Truncated proposals over this surface"** (2026-09-09, W3.4; Peter's ruling of the same date to express TMNRE through sbi rather than revive swyft; accepted at the Fable review; **ratified by Peter 2026-09-10** — the coverage row, ε and the final-round-only pair estimators confirmed, the sampler default ruled below). The results question resolves **against** the marginals-as-posterior reading: a TMNRE run's `posterior` group is ordinary joint draws, from a joint ratio estimator trained across the rounds and multiplied by the *final truncated prior*, scored on the numpy contract path exactly as every other run's are — so §4.6 is untouched and nothing downstream learns a second shape of run. The marginals sit **beside** them as a fifth reserved group, and the first stored **by default**: the other four can be recomputed from a stored run plus the problem, while one ratio estimator per marginal exists only inside the fit — and marginal ratio estimation is half of what the method *is*. What is stored is a bounded **summary** (`(free_size, 129)` for 1-D, `(pairs, 33, 33)` for pairs): each estimator's `log_ratio` and estimated marginal posterior `log_density` on a grid over the final box, in both parameterisations; never the networks. The truncation history rides in the attrs (`ampere_sbi_truncation`, one record per round). Six decisions: (1) the marginal prior the ratio is reconstructed against is a **Gaussian KDE** of the trained-on column, not a histogram — the number it feeds is a threshold crossing, and a KDE's over-smoothing biases the box *outwards*, the safe direction; (2) the box is the product of the 1-D intervals above `ε·max` over **all** crossing nodes (a bimodal marginal keeps both modes), intersected with the incoming box so the sequence is nested by construction; default `ε = 1e-4`; (3) the joint estimator trains on **every** round (`discard_prior_samples=True` measured to cost a third of the budget for no gain); (4) `sample_with="rejection"` is the default because i.i.d. draws are what the single "chain" promises — **not** because it is cheap: measured **353.6 s against 43.8 s** for `"mcmc"` on the example, since rejection's cost rises as truncation succeeds — **ruled 2026-09-10: `"mcmc"` is the default (W3.16), `"rejection"` stays selectable**; (5) `calibrate()` on a TMNRE run draws its batch **and** its reference draws from the *truncated* prior and checks through an MCMC posterior over the same estimator, because a rejection posterior pays a fixed maximisation per conditioning observation; (6) `ampere_sbi_amortised = 0` for every TMNRE run, `= 1` only for a single-round NPE/NLE/NRE fit. **Accuracy row**: the 0.5σ location tolerance is unreachable for a ratio estimator's joint on the toy problem (worst parameter 0.05–1.15 reference σ over five fits); the suite asserts **coverage** (the emcee mean inside the run's central 95 %) plus W3.2's width band — a strictly stronger pair — and holds the marginal estimators to 1.0σ (worst 0.54). **Shortfall carried to W3.12**: `artefact_key` has no `marginals`/`truncation_epsilon` ingredients, so `_key_architecture` folds them into the `architecture` ingredient as a stopgap; the proper diff is in W3.4's report. `sbi` 2706/70, `dev` 1943/248 on the branch. |
| `LikelihoodFamily.sample`: the refusal becomes the exception for three families (W3.14) | **`poisson`, `student_t` and `complex_gaussian` implement `sample` in `ampere.core` and have native twins on torch and jax; `cauchy` and every declared-but-unimplemented family keep the refusal word for word** (§4.4 contract change, 2026-09-09, W3.14; approved by Peter the same day on W3.6's use case — a counting experiment could not `simulate(observe=True)`; accepted at the Fable review; **ratified by Peter 2026-09-10**). §13's refusal was written for families whose observation process is genuinely ambiguous; these three are not — each has exactly one generative form and its own `log_prob` fixes which. Poisson draws `rng.poisson(rate)` with the rate `predicted·exp(f)` under a latent GP, `f` the whitened latent read out of θ (so a latent-consuming family draws at the rate its `log_prob` scores, or SBC ranks are wrong); Student-t is location-scale with σ the scale and correlated noise refused by name (a scale mixture does not commute with a GP covariance); complex Gaussian draws independent components with **σ per component** (total variance 2σ², verified against the density). The backend **ceiling** is unchanged (a backend samples exactly what the core samples) and W3.14 states the **floor**: a family the core samples but a backend has no twin for falls back to the numpy path, never a refusal. torch's twins for Poisson and Student-t are two-stage (parameters under `vmap`, variates after it per draw from that draw's seed, `torch.distributions` inside `fork_rng`) because a random op under `torch.func.vmap` is an error or globally seeded; jax needed nothing. The draw's working dtype comes from the containers so a complex family draws complex observations. `examples/wstat_comparison.py`'s `CountingPoisson` deleted with the study's ranks pinned unchanged. **Ruled 2026-09-10 (Peter)**: by *principle*, not by list — a sampling form is added when a data type needs its likelihood, and then for every inference approach on every backend at once, never for one or a few (`likelihoods.md` §3 states it); `cauchy` waits for the data type that needs it, exactly as polarimetry will bring Rice's sampling form and likelihood together; and `sample`'s guard is aligned to `<= 0`, the more conservative one (W3.16). Gates on the branch: dev 1937/227, jax 2432/160, torch 2608/156, sbi 2684/80. |
| Model identity hash promoted to provenance (W3.12) | **`ampere_model_hash` joins the root attributes of every run and training set; `PROVENANCE_SCHEMA_VERSION` → 6; `append_training_set` refuses a model-hash mismatch and a pre-schema-6 file, both by name; `ArtefactKey` gains `marginals`/`truncation_epsilon`/`sample_with`** (`results.md` §4, §9, §11 amended; 2026-09-09, W3.12; ruled by Peter 2026-09-09, accepted at the Fable review; **ratified by Peter 2026-09-10** — schema 6 and the `model_hash` name confirmed). `ampere_spec_hash` is only the merged parameter declaration, so a likelihood family, noise model, solver or kernel swap that leaves every parameter's name and prior unchanged was invisible to W2.8's append check — plan §7's poisoned cache arriving by the other door. The composition W3.5 built privately (every model fingerprint minus `"parameters"`, every dataset fingerprint minus `"observed"`, plus the bindings) is now `provenance.model_hash(problem)`, called from the cache key and the append check alike. **Named `model_hash`, not the item's `model_identity_hash(problem)`**: `model_identity_hash(model)` already exists (W2.1) with different, cross-backend "offer" semantics, and reusing the name would have been a silent collision. A pre-schema-6 file has no hash to compare and is **refused**, not passed — the conservative reading, one `write_training_set` away from a fresh file. W3.4's carried key diff lands here **with the third folded setting the item text omitted, `sample_with`**: sbi bakes the sampling mode into the built posterior, which is exactly the stored artefact, so two TMNRE runs differing only there must not share a digest. `ARTEFACT_CACHE_SCHEMA_VERSION` unchanged (an optional key absent when unset, so every pre-W3.4 non-TMNRE digest is byte-identical — pinned in the tests). Branch `tests/results` 363/3, lint/format/pyrefly clean; merged-master gates in the status row. |
| SBI runs reproducible from the problem's seed (W3.15) | **`SBIEngine.run` seeds torch's global generator — and numpy's legacy global one — from the problem's own sub-stream before every network build, training round and posterior draw, restoring both on exit; `ampere_sbi_torch_seed` recorded, absent when `problem.seed is None`** (`inference.md` §13 amended; 2026-09-10, W3.15; drafted from W3.4's carried finding and dispatched on Fable's judgement; accepted at the Fable review; **confirmed by Peter 2026-09-10**). Neither ampere nor sbi seeded torch, so two runs of one seeded problem agreed on every simulated pair and disagreed on the posterior trained from them — W3.5's "identical posterior" held only through the cache. The item's own TMNRE acceptance check then refused to repeat with torch alone seeded: sbi's default MCMC method draws its slice proposals through `np.random` directly, so the same context manager seeds that too. Both generators are restored afterwards — the rule `_nuts.py`'s `fork_rng` and `_zeus.py`'s `_global_seed` already state: ampere leaves no library's global random state changed behind it. The final draw is reseeded on its own concern so a cache hit, which trains nothing, samples reproducibly as well. Not reached, deliberately: `calibrate()` never retrains, and its internal `sbi.diagnostics` sampling is unchanged. Branch gates dev 1982/264, sbi 2766/80. |
| Peter's Phase 3 rulings applied (W3.16) | **TMNRE samples by MCMC by default (`TMNRE_DEFAULT_SAMPLER = "mcmc"`, `"rejection"` selectable); `PoissonFamily.sample` guards `rate <= 0` exactly as `log_prob` does; and the principle for sampling forms: a data type's likelihood arrives with its sampling form, for every inference approach on every backend at once** (2026-09-10, ruled by Peter on the Phase 3 questions block; Fable-authored, merged at `fea7ef9`). The sampler default follows W3.4's own measurement (353.6 s against 43.8 s, and rejection's cost rising as the method works); the canonical TMNRE fixture now runs on the default and the calibration fixture asks for rejection explicitly, since the rebuild path it exercises is the rejection posterior's. The guard follows the rule that `simulate(observe=True)` must never hand an SBC or SBI consumer a draw the fitting likelihood cannot score; the native twins draw inside a trace and cannot raise, so a non-positive rate there is scored to −inf at evaluation as before. The principle closes the `cauchy` question by making it wait for its data type. Also ruled the same day, recorded where they act: `actionlint` and a path-gated CI split (W4.10); the overview refresh and a bare NPE example (W4.0); the stored proposal density, `DataTree` optima, BO as acquisition only, samplers behind extras, three new design horizons (§5) and the results-contract item they need (W5.0) — from the inference-extensions memo. Gates on merged master in W3.16's status row. |
| Phase 3 documentation pass: stale annotations corrected in place (W3.13) | **Five stale sentences across `docs/design/architecture.md` §3 and `docs/design/contracts/{inference,diagnostics,results}.md`, plus this plan's own §5 Phase 3 bullets and §6's SBI deferred bullet, corrected in place and marked *Amended W3.13*; no contract semantics change** (2026-09-10; docs-only, `pixi run docs` unchanged at 13 warnings, byte-identical set; branch gates dev 1994/268, sbi 2782/80; accepted at the Fable review with one correction). `architecture.md` §3's "the SBI layer is not here" corrected now that `SBIEngine` landed (W3.2, W3.4); `inference.md` limitation 17.5 and its §18 Phase 3 bullet closed (batched `simulate_many`, W3.1; the encoding, W3.3); `diagnostics.md` row D landed at W3.6; `results.md` §18 corrected rather than dated — it named `ampere_problem_hash` as the SBI cache key, which is wrong: a training set is simulated from the prior and is checked on `ampere_spec_hash` + `ampere_model_hash` (W3.12), while the trained-artefact key adds `data_hash` of the observed containers because a stored posterior is conditioned on its observation (the review's correction of the agent's "same two hashes"). §5 gains the landed-summary paragraph (sbi 0.27, no swyft, the encoding, `forkserver`, schema 6); §6's SBI bullet loses swyft with Peter's ruling of 2026-09-09, keeps jax-native SBI, and gains the embedding-network study with the ε study folded in. |
| Where a shipped observable lives (Phase 4 D1) | **The kind in `core/results_schema.py`; the steps in `backends/{reference,torch,jax}/<observable>.py`; no grouping namespace** (ruled by Peter 2026-09-11 on `docs/design/phase4_placement_memo.md` §2). The draft's `ampere/modalities/` was rejected on the name and on a rule it broke — it moved the reference backend's steps out of `backends/reference/`, so the three backends' step names were no longer parallel. A kind is three class attributes and not a contract, so adding one to core is not a §4 change; the out-of-tree extensibility claim remains proven by `tests/core/thirdparty_polarimeter.py`. A grouping namespace is revisited after realistic usage (end of Phase 4 or later); a per-observable front door (`ampere.interferometry`) arrives only with the first reader (OIFITS, Phase 6). |
| The closure-phase signature, and wavelength as an axis (Phase 4 D2) | **`VisibilitySet` gains a third axis `spectral_axis`; `ClosurePhases` is `(u1, v1, u2, v2, spectral_axis)` with a canonical baseline ordering; kernels gain an `axes` selector so a `Product` of a (u, v) block and a spectral block is expressible** (ruled by Peter 2026-09-11, memo §3.6–3.7, after his question on chromatic misspecification). A missing band in a patch of sky gives δV = S(λ)·F(B/λ): sharp in wavelength, smooth in (u, v); two baselines of different length at the same (u, v) sit at different wavelengths, so the B/λ convention absorbs wavelength only for a grey error, and a kernel sees a container's axes only. The `VisibilitySet` amendment is a frozen-kind change carried by W4.1 with the conformance update in the same PR. A GP on *closure phases* is a latent composition (the family is wrapped) and is Phase 5's (W5.1); Phase 4's flagship GP is on the visibilities. The kernel space is dense-only (3 and 5 axes); in (B, λ) coordinates a dispersed observation is a product structure whose covariance is a Kronecker product with a quasiseparable spectral factor — Phase 5's structured-solver route, kept recoverable by the baseline labels and the wavelength axis. |
| Kernel algebra, the axes selector and one term registry (W4.5) | **The quasiseparable-term registry is keyed on the kernel family alone, not on (family, backend)** (merged 2026-09-11, `6a78a2b`). A semiseparable representation is generators built from the coordinate and the hyperparameters with `exp`, `cos`, `sin` and `stack` — the same mathematics in numpy, torch and jax — so the array namespace became an argument (`ArrayOps`, carried by the kernel instance, rebound by each solver through `Kernel.with_ops`) rather than a second registry key: a user registers **once** and reaches the O(N) path on all three backends, and the three private tables each backend carried collapse into one public registry of `lowering.py`'s shape. The kernels' closed forms moved into `ampere.core.kernels`, written against `ArrayOps`, so the backends' per-family transcriptions are gone. `likelihoods.md` §6–§8 amended; §15.1 and §15.3 lifted, §15.2 half lifted. **§15.3's claim that Matérn-5/2 is not exactly quasiseparable is withdrawn**: true of the celerite basis, false of the semiseparable form the solver factorises, in which Matérn-5/2 is exactly rank 3 (the correction W2.3 made for Matérn-3/2, one order further). Every kernel gains an `axes=` selector (D2, memo §3.6 item 3) defaulting to "every axis" so no existing declaration or spec hash moves; the single-unit rule applies per leaf to the selected subset; `Product` composes kernels on disjoint subsets on the dense path and is refused by name on the O(N) path. Composite hyperparameters are qualified by term label in declaration order (positional labels are not Python identifiers, and sorting would let incomparable chains share an identity). **Confirmed by Peter 2026-09-15.** |
| `VisibilitySet` gains a spectral axis (W4.1) | **A frozen kind amended: `VisibilitySet` is `(u, v, spectral_axis)`, one wavelength per sample, `Order.ANY`** (merged 2026-09-13, `eb7b702`; the D2 ruling of 2026-09-11 executed). A chromatic sky error is sharp in wavelength and smooth in spatial frequency, a kernel sees a container's axes only, and (u, v) in wavelengths absorb λ only for a grey error — so wavelength is a coordinate, not a label. The monochromatic case is a constant column; the alternative of a second dispersed kind would double what every step and family accepts. Every existing `VisibilitySet` conformance row updated in the same PR (ground rule 9). `ClosurePhases` lands beside it with five axes and a canonical baseline ordering in its docstring, asserted by the conformance battery against a closed form. **Confirmed by Peter 2026-09-15.** The interferometric conventions (baseline sign `b_ij = r_j − r_i`, the canonical ordering, the DFT sign with x east and y north, position angle east of north) **accepted by Peter 2026-09-15**. |
| The wrapped family implemented, with its draw (W4.1) | **`VonMisesFamily` is implemented — the normalised density with `scipy.special.i0e`, the residual wrapped into (−π, π], `κ = 1/σ²` per sample, radians enforced at composition — and gains `sample()` (`rng.vonmises(mu, κ)`) under the sampling-form principle of 2026-09-10** (merged 2026-09-13). A GP on a wrapped observable is a latent composition; the family refuses `GaussianProcessNoise` by name and points at Phase 5's W5.1. **Confirmed by Peter 2026-09-15.** |
| The astropy adapter, a post-freeze §4 addition (W4.6) | **`core/astropy_compat.py` implements §4.7 as `contracts/astropy_compat.md` states it** (merged 2026-09-13, `b6d9298`): the adapter is black-box on the reference backend (`DIFFERENTIABLE = False`, `BATCHABLE = False` — astropy vectorises over inputs, not over parameter vectors), bounds become uniform priors, fixed parameters are frozen, ties are recorded and applied rather than sampled, an unbounded free parameter without a prior is refused; **no solid angle is ever invented** — a per-steradian output converts to a flux density only through an equivalency the caller states; the kind is declared or inferred only where the axis unit and arity make it unambiguous. The opt-in native translation of 2026-09-01 is a per-backend `from_astropy` that refuses anything outside its curated table; the table is empty until W4.7. |
| The circular complex GP implemented, and a solver's right-hand side widened (W4.2) | **`complex_gaussian` + `GaussianProcessNoise` is `ANALYTIC` *and implemented*: the circular (equal-component, zero-pseudo-covariance) complex GP, computed as one Cholesky of `S = K(θ) + diag(σ²)` with a two-column right-hand side rather than as a `2N` by `2N` factorisation** (merged 2026-09-13, `4d93323`; the 2026-09-03 ruling of `likelihoods.md` §17 Q6 executed, `GP_ANALYTIC_IMPLEMENTED` now `True`). Circularity says the real `2N` covariance of `(Re r, Im r)` is `diag(S, S)` with a zero off-diagonal block, so the density is `−½[rᵉᵀS⁻¹rᵉ + rⁱᵀS⁻¹rⁱ + 2 log|S| + 2N log 2π]`: one factorisation, two solves, `log|S|` once and counted twice. Materialising the `2N` matrix costs eight times the arithmetic to carry a zero block the model has declared; calling the real solver twice computes the log-determinant twice. **The `GPSolver` contract is amended accordingly** (`likelihoods.md` §7, new subsection): a `residual` may be an `(n, k)` block of `k` realisations independent of one another and sharing one covariance — `log_marginal_likelihood` sums the `k` marginals, `conditional_loo` stays one term per sample (the components share `A_ii`), `condition` gives an `(m, k)` mean and one `(m,)` variance, `latent_transform` maps `(n, k)` to `(n, k)`. `GPSolver.STACKED_RESIDUALS` is the opt-in, `False` by default because both failure modes are silent: flattening scores `2n` residuals against an `n` by `n` covariance, and taking the first column drops the imaginary part of every visibility. **`QuasisepGP` is refused by name, and it is structural rather than scheduling**: `REQUIRES_ORDERED_1D` cannot hold for a point of the `(u, v)` plane at a wavelength, whatever the kernel selects, so the refusal precedes the solver's own check — `interferometry.md` §7's prediction that the O(N) path "will inherit it" is withdrawn. σ and the kernel `amplitude` are both per-component (`results_schema.md` §16), so `E|r|² = 2(K_ii + σ²)`; the draw is two real realisations sharing `L` and nothing else. Two native-path bugs closed in the same PR, exposed by the first multi-axis GP customer: both backends' lowerings took `observed.axes[0]` as the GP coordinates rather than the full `(n, d)` stack, and both used the unbound `noise.kernel` where W4.5 put the axis binding on `kernel_for(observed)`. **Confirmed by Peter 2026-09-15.** The calibration row's coverage pin (rather than the rank test, where the flexible GP over-covers) confirmed the same day. |
| The native model surface has a second spelling (W4.3) | **`ampere.backends.{torch,jax}.problem` resolve a native model's value-and-coordinates pair over two names: `flux`/`grid` (the original) and `native_flux`/`native_grid`** (merged 2026-09-13, `64284d5`). A genuine collision rather than taste: `Parameterised._check_free_name` refuses a parameter whose name shadows a class attribute, and `flux` is exactly what every interferometric source model calls its total flux density, so such a model cannot have a method called `flux` and the convention the placement memo records as a duck-typed `flux`/`grid` surface was unsatisfiable for the first modality that needed it. Either pair composes; a model offering half of either is refused with the missing half named. Not a §4 change (the surface is a backend-internal convention), recorded because it is a documented convention. **Ruled by Peter 2026-09-15: (a) and (b) together, as W5.20 — `native_flux`/`native_grid` canonical with `flux`/`grid` an alias, and the shadow check narrowed to a backend-invariant core reserved set** (provisionally (a) earlier the same day, made final on the consequences analysis he asked for) (recorded in `docs/development.md`'s review block of 2026-09-15). The analysis's finding: `_check_free_name` tests `hasattr(type(self), name)` over the whole MRO, so each backend's base class adds its own reserved parameter names — core `Model` reserves 17, a torch spectral model adds `AXIS`, `evaluate_tensor`, `flux`, `grid`, `grid_tensor`, `to`, a jax one `AXIS`, `flux`, `grid` — and a parameter legal on the reference backend can be illegal on its twin; separately, torch's lowering nests parameters as `nn.Module` attributes, so `nn.Module`'s own namespace (`to`, `type`, `float`, `apply`, …) is a third reserved list enforced only at lowering. Nothing reads a parameter as an attribute (no `__getattr__`; values arrive as `context[name]`), so the rule guards a convention, not a live defect. The lookup order `("flux", "native_flux")` means a class offering both is silently taken at the legacy name. W5.20 (proposed) makes the namespace backend-invariant; the canonical spelling flips the lookup order and refuses a model offering both pairs. Also recorded: both modern backends' interferometric twins use the *inheriting* pattern — the declaration (requirements, `configure_from`, the canonical-closure check) is written once, because a pixel-scale requirement written twice would alias silently. |
| Phase 4 documentation pass: stale claims annotated in place (W4.8) | **Nineteen additive *Amended W4.8* annotations across ten frozen design documents; no contract semantics change; the conformance suite untouched** (merged 2026-09-13, `c8bf2c0`). Every "Phase 4", "landed", "declared, not implemented" and "will" claim in `likelihoods.md`, `results_schema.md`, `results.md`, `transformations.md`, `inference.md`, `encoding.md`, `architecture.md` and the interferometry, spectrum-photometry and astrometry sketches was checked against the merged code. One is a correction: the interferometry sketch's Q2 ruled that `extra_coords` would gain `Axis` support for a per-visibility frequency, and Phase 4 met the need with a `spectral_axis` on the container instead (D2), so `results_schema.md` §15.4's extension point stands unspent. Two new limitations recorded where they had no home: `results.md` §13 item 14 (four of the six plots refuse a point kind with several axes — Phase 5's multi-axis diagnostics) and `encoding.md` §9 item 6 (the packing aligns axes by position, not name). The template page `interferometry.rst` is amended for the four gaps W4.9 found by following it. |
| The results contract for approximate and evidence-producing engines (W5.0) | **`ampere_approximation` in the root attrs of every run (`"none"` for an exact sampler; `"mean_field"`/`"multivariate"` for `VIEngine`, `"density_estimator"` for every `SBIEngine` method), the engine-neutral evidence triple `ampere_log_evidence`/`ampere_log_evidence_err`/`ampere_evidence_method` (dynesty today, `"nested_sampling"`; `ampere_dynesty_logz` kept as the engine's own), and `sample_stats.proposal_log_density` — a proposal's own log-density per draw, in the same constrained coordinates the stored `log_prior`/`log_likelihood` are in, via one shared `unconstrained_jacobian_correction` (`log q_θ = log q_u − log|∂θ/∂u|`); `PROVENANCE_SCHEMA_VERSION` → 7** (ruled by Peter 2026-09-10 on the inference-extensions memo §5, §7.1–7.2; merged 2026-09-16). Written once — the default at `Engine.finish()`, the hook on `emit(sample_stats=)` — so every later engine is one item. SBI's `ampere_sbi_log_prob` (unconstrained) stays beside the new field; the two are not aliases and the module says so. `plot_trace` and the new `ampere.results.summary` warn when the approximation is not `"none"`. The accept criterion is the contract's own check: an importance-corrected VI posterior computed from the stored groups alone agrees with an emcee reference (`tests/inference/test_vi.py`). **Found at the Fable review by running the branch, fixed before merge**: both VI routes stored the wrong density — pyro's autoguide reports the model site as a `Delta` whose `log_prob` is the change-of-variables term (zero here), and numpyro's `AutoNormal` has no `get_posterior` — so the density is now the sum of `log_prob` over every sample site of the traced guide on both routes, with a test against an independently built guide density (to 1e-6) and a non-constancy guard; pyro's `AutoMultivariateNormal` covariance is `scale[..., None] * scale_tril`, not `scale_tril` alone. |
| A backend-invariant parameter namespace (W5.20) | **The reserved parameter names are core's, stated and finite — the public names of `Parameterised` and `Model` (computed from the classes) plus the pinned `torch.nn.Module` namespace (`TORCH_MODULE_NAMES`, 53 names plus `training`, checked against the real class by a torch-environment test) — and identical on every backend; a backend may add public attributes to its model classes without changing which parameter names are legal. The native value-and-coordinates surface is canonically `native_flux`/`native_grid`; `flux`/`grid` is a supported legacy alias; a model offering both is refused as ambiguous** (D1 ruled by Peter 2026-09-15 as options (a) and (b) together; merged 2026-09-16, `acaad79`). `hasattr(type(self), name)` had made a core contract's namespace rule depend on which backend's base class a model inherited — `flux` legal on the reference `UniformDisc` and illegal on its twin, `grid` reserved on torch and free on jax — while nothing in ampere reads a parameter as an attribute, so the rule guarded a convention, not a defect; torch's lowering, which nests parameters as `nn.Module` attributes, now refuses a colliding leaf in ampere's words instead of torch's `KeyError`. Keeping `native_*` after the constraint is gone is deliberate: a surface a realisation composes and a quantity a user fits should not compete for one word. No spec hash moves (a pinned hash asserts it). **Behavioural breaks, for the release note**: a subclass of a shipped model that overrides `flux` beside the inherited `native_flux` is now refused as ambiguous (it fired once, in the jax conformance fixture); and names in `torch.nn.Module`'s namespace (`type`, `eval`, `train`, `float`, `apply`, …) are refused at declaration on every backend where before they failed only at torch lowering. |
| The first approximate solver: `HilbertSpaceGP` (W5.4) | **`GPSolver` admits `EXACT = False`, held to a *convergence* rather than a tolerance, and the latent block is the solver's to size** (merged 2026-09-16, `359386c`; Opus-authored, Fable-reviewed as the second pass; terra owed). Two contract questions `horizon_notes.md` §2 left open are settled. (a) The conformance class for an approximate solver is an envelope anchored on the error measured at the coarsest setting of the same sweep and tightening as the approximation is refined (`coarsest_error · (m₀/m)^order + floor`), plus a separate absolute claim at the finest — because a fixed number would either assert nothing or encode one fixture's arithmetic, whereas a wrong basis, a spectral density with the wrong dimension or a dropped Woodbury factor each give an error that is *stable* under refinement and sails through any fixed tolerance. The rate is the kernel's property (an isotropic Matérn-ν converges as m^−2ν: order 1, 3, 5 for ν = 1/2, 3/2, 5/2), so it is per row with the slowest supported case as the default; Matérn-1/2 genuinely reaches only 1.9e-1 nats at m = 128 on the 12-point grid and says so in its row — the regime the horizon note reserves for a Vecchia solver, now with the measurement to prioritise it. (b) The latent size is `GPSolver.latent_size(kernel, n_samples)` — `n_samples` for both exact solvers, the basis size `m` for a reduced-rank one — so `inference.md` §17.4's "one value per retained sample" was a property of the two exact solvers and not of the contract; `LatentDeclaration.size` keeps meaning the retained count the mask is checked against and `whitened_size` is the block the sampler carries, and one whitening serves `simulate(observe=True)` and the latent path by construction (a latent-GP problem has 15 sampler dimensions at m = 12 where the exact solver has 51). **Against the item's and the plan's wording, `basis_size` and `boundary_factor` are dataclass fields and so enter the spec hash** — fold-in 10's test is whether two backends may legitimately differ on a value: they may on a dtype, they may not on the approximation, and two runs at different `m` are not the same model; they are repeated in `provenance_config()` for the reader. The §5 bullet and `horizon_notes.md` §2 are corrected accordingly. **Merged on Fable's judgement; for Peter to confirm or strike**: the hashing decision, and `SpectralMixture` supported rather than refused (it is a `Sum` of `SHO`s, whose spectral density is exact). `conditional_loo` is not refused: Woodbury gives `diag(A⁻¹)` in closed form at O(N m²). |
| Multi-axis diagnostics (W5.3) | **The family B/C plots and derived groups reach a point kind with several axes through one rule: `coordinate=` (an axis name, or a callable from the kind's axes to `(coordinates, label)`) resolved against a `FunctionSamples.PLOT_COORDINATE` default declared on the kind, else refused by name listing the axes; a complex dataset is derived and plotted as the `component=` the caller names (`real`/`imag`/`abs`/`phase`), stored as `<label>_<component>`, and refused with the four choices when none is given** (merged 2026-09-16, `47988b4`; Sonnet-authored, Fable-reviewed; lifts `results.md` §13 item 14). `VisibilitySet` defaults to baseline length and `ClosurePhases` to the longest of its three baselines; a single-axis kind declares nothing and resolves exactly as before, so every existing plot is byte-identical (pinned). The default is precomputed from the live container when the derived group is built and stored on the group beside the kind's raw axes — plain xarray variables on the `<label>_index` dimension, assigned after `arviz.from_dict` because arviz assumes every named variable carries `(chain, draw)` — so a run loaded from netCDF plots without the container (verified at review by a round trip). `plot_anomaly_score` is unchanged: it takes a built `AnomalyScore` whose coordinates are already one-dimensional, the resolution having happened upstream in `gp_localisation_score`, which gained the arguments instead. Not a §4 shape change (a kind attribute, per Phase 4's D1); the conformance suite is untouched. Two pre-existing bugs fixed on the way: `_observed_axes` mis-zipped a multi-axis kind's axes against the index dimension, and `_draw_whiteness` passed `width=None` explicitly, which a recent matplotlib refuses. |
| The first gridded observation: `PSFConvolution` and `propagate_mask_grid` (W5.5) | **A §4 transformation addition and a §13.5 amendment, both within the frozen surface** (merged 2026-09-16, `614d7c9`; Opus-authored, Fable-reviewed; terra owed). `transformations.md` §10's PSF row is filled by a `Transformation` that publishes *one* constraint in *two* forms — `points=` at the observed pixel centres padded by the kernel's half-support, and `intervals`/`max_step` restating the pixel scale as a density — because a PSF, unlike an LSF, has no resampler downstream to land it on the observed grid. The two must agree to the last bit or their negotiated union is not evenly spaced and the FFT route refuses on whichever backend ran, so the padded grid is built the way `AxisRequirement.coordinates` builds one: **if a step publishes the same constraint twice, one form must compute the other.** §13.5 limitation 5 ("mask propagation is 1-D") is amended rather than lifted: `propagate_mask_grid` changes the *expression* of the ANY rule, not the rule — a local separable kernel's influence is a neighbourhood, carried as a dilation by the kernel's bounding half-support at O(kN) instead of an O(N²) influence matrix — and a grid step whose influence is *not* local (an arbitrary warp, a non-separable resampling) still flattens its own mask. `ifu_cube.md` gap 2 closes by the route its freeze disposition named. Conformance green in all four environments. |
| `Layout.GRID` and correlated noise: a gate, not a limit (W5.5) | **`ampere.core` refuses a correlated noise model on every `Layout.GRID` container, and the refusal is a closed Phase 5 slot rather than mathematics** (`GPSolver.check_compatible`, `likelihood.py:588`, and `Likelihood._coordinates`, `:4047`, both "`LAYOUT is not Layout.POINTS`"). A stationary kernel over `(x, y)` is a function of coordinates, and a grid keeps its coordinates separably, so the `(N, 2)` matrix must be broadcast out of the axes rather than column-stacked — which `ampere.core.encoding` has done since W3.3 and W5.5 lifted into `ampere.core.sample_coordinates`. W5.5's ownership put `likelihood.py` out of scope, so `examples/image/grid_gp.py` lifts the two gates *out of tree* from the public API (three subclasses overriding only the layout check), as `transformations.md` §11 exists to prove is possible. **Proposed, and ruled in by Peter 2026-09-16 as W5.21**: the library change is two edits in `ampere/core/likelihood.py` — allow `Layout.GRID` in `GPSolver.check_compatible` when the kernel's selected axes are the container's, and build `Likelihood._coordinates` with `sample_coordinates` — after which `grid_gp.py` is deleted and the study imports `DenseGP` and `HilbertSpaceGP` directly; it touches `likelihoods.md` §7 and needs its own row when it lands. Nothing here exploits a grid's structure; that is W5.6's bake-off. |
| `sample_coordinates`: one coordinate matrix for both paths (W5.5) | **What a kernel sees on the native path and on the contract path is now the same function.** W4.2 fixed `observed_coordinates` on both modern backends' `LoweredProblem` by column-stacking the axes — right for a three-axis *point* kind, wrong for a grid, where the axes are the grids the samples are the product of. `ampere.core.sample_coordinates` is that stacking made layout-aware by `encoding.py`'s rule, and is now the single implementation used by `Dataset.draw_observation` and `LoweredProblem.observed_coordinates` on torch and jax (and, under the change proposed above, `Likelihood._coordinates`). Two further `Layout.GRID` assumptions were fixed alongside it in both backends' `problem.py`: values and uncertainties masked before being flattened, and `predict()` applying a flat mask to a container-shaped prediction; the point-set branch of `predict()` is deliberately unchanged. |
| Population inference by reweighting archived fits (W5.13) | **`ampere.results.population` closes design horizon (b)'s "buildable entirely on stored runs" promise: `RunColumns` (a `typing.Protocol`) plus two implementations — `DataTreeRunColumns` (in memory) and `NetCDFRunColumns`/`runs_from_netcdf_directory` (a directory of archived `.nc` files) — read the per-draw `log_prior`/`log_likelihood` every run stores and, for approximate engines, W5.0's `proposal_log_density`; `fit_population` self-normalised-importance-reweights them (Hogg, Myers & Bovy 2010) under a declared `PopulationModel` (`GaussianPopulationModel` shipped) into a population-hyperparameter posterior sampled with `emcee`, refusing by name when any object's effective sample size at the posterior mean collapses below a stated floor or when the runs do not share `ampere_spec_hash`** (merged 2026-09-16, `81442e5`; Sonnet-authored, Fable-reviewed). The identity divides each draw by the named parameter's *marginal* interim prior — the nuisance priors cancel — and that prior is a required `interim_prior` argument, because a run's provenance records parameter names and a hash of the declaration, never the `PriorSpec` (found and fixed at review; the joint `log_prior` column is used only where it is right, turning an approximate engine's proposal draws into posterior draws). Not a `PROVENANCE_SCHEMA_VERSION` change: the population posterior is a new derived `DataTree`, not a change to what a single-object run stores; `results.md` §13 item 16 carries the account, including what is deliberately left out (per-object nuisance re-sampling, multi-level populations, a selection function, and the joint-fit sibling route W5.12, which this item cross-checks against rather than merges with). Conformance suite untouched. Proposed for a later row: store each parameter's `PriorSpec` in provenance so the reader can verify the supplied interim prior. |
| Phase 5's second approximate solver, chosen by measurement (W5.6) | **Neither EFGP nor Vecchia is promoted at W5.6; `HilbertSpaceGP` stays the phase's reduced-rank solver, `EquispacedFourierGP` is recommended for a follow-on three-backend item, and `VecchiaGP` stays a slot** (merged 2026-09-16, `8d785e1`; Opus-authored, Fable-reviewed). §5's rule is that the first 2–3D solver landed wins a benchmark at realistic N, and the benchmark exists: `examples/image/bakeoff.py`, attached to `pixi run bench` by `tests/benchmarks/test_solver_bakeoff.py`, with both candidates implemented in the reference path (`ampere/core/efgp.py`, `ampere/core/vecchia.py`) and the full table in `likelihoods.md` §7. **EFGP is real but bounded.** Its equispaced grid makes the Woodbury normal matrix block-Toeplitz, so the assembly drops from O(N m²) to O(N 2^d m) and the `(N, m)` block is never formed — measured at N = 16 384 it is faster than HSGP past m ≈ 576 (0.59 s against 0.92 s at m ≈ 1024; 2.64 against 3.37 at m ≈ 2400) at 4.7× less memory (58.7 MiB against 273.3), and at matched frequency reach it is *more* accurate on every smooth kernel measured (|Δlog p| 1.1e-2 against 2.1e-1 on the image's own Matérn-3/2). But the published O(N + m log m) is a claim about **posterior-mean regression**, and a `GPSolver` is not that: the same normal equations at m = 6561 cost 0.16 s and 2.6 MiB by conjugate gradients on FFT matvecs against 10.95 s and 1.31 GiB by Cholesky, and the entire difference is log|M|, for which a Toeplitz matrix has no FFT; the alternatives are the O(m³) factorisation or stochastic Lanczos quadrature, and a stochastic log-determinant is a *noisy* likelihood, which every engine `inference.md` describes treats as a correctness failure. EFGP therefore buys a constant factor on the same complexity class — worth a follow-on item, not worth pre-empting the phase's other work. **Vecchia stays a slot on evidence, not on schedule.** It is the right method for a short-range rough process (on a 1-D Matérn-1/2 with ℓ one eighth of the span it reaches 1e-6 nats at k = 32 where a 513-point Fourier grid is still 8.5 nats away), but the flexible likelihood absorbs *unmodelled structure*, which is broad by construction, and there the screening effect (Stein 2002) fails: on the image's own kernel it is 24–32 nats from exact at k = 30 with localisation 0.68 against the spectral solvers' 0.90, at 25–35× the cost per evaluation at 65 536 pixels (14.3 s against 0.41–0.64 s); its latent block is N rather than m, the opposite of what the NUTS path needs. Both prototypes stay in the tree, tested and documented, so a modality that *is* short-range and rough reopens the question with the measurement already written. |
| Non-stationary flexible likelihood: `WarpedKernel`, and the coordinate a recursion runs on (W5.7) | **`ampere.core.WarpedKernel(base, input_warp=, amplitude_warp=)` makes the flexible likelihood non-stationary without leaving the exact O(N) path on any backend** (merged 2026-09-17, `97f7662`; Opus-authored, Fable-reviewed). A monotone input warp `w` and a diagonal amplitude warp `a` give `a(x) k(w(x), w(x′)) a(x′)`: the base family's registered generators evaluated at `w(t)` and scaled by `a`, at unchanged rank and unchanged decays, so W4.5's registry is the extension point and nothing new is approximated. **Three additive contract changes, each omitted when unused so every spec hash minted before W5.7 is unchanged** (checked against a literal recorded on `a6a09c3`: `Matern32(0.4, 2.0)` → `2ff5c78e4dfe672df7e0eb81a156eeaf`; re-derived independently at review). (i) `KernelSpec.metadata` carries the knot *locations*, which are declaration rather than hyperparameters — nothing samples them — yet make two warps different models. (ii) **`Kernel.warped_coordinate(axis, values)`, the identity everywhere else, is the one hook the registry lacked and it could not be avoided**: celerite2 forms its propagators `e^{−cⱼ(tₙ−tₘ)}` from the coordinate handed to `compute`/`factor`, never from the term, and folding the difference into `U`,`V` means multiplying by `e^{±ct}` — precisely the overflow the factored form exists to prevent. It is a `Kernel` method rather than a builder-signature change so that `(kernel, values, axis) → CeleriteRepresentation` and the meaning of `axis` (the raw sorted coordinate) are untouched for any representation registered before W5.7; `CeleriteRepresentation.marginal` may now also be `(n,)`, which every consumer already broadcasts. (iii) `GPConditional` gains `coordinates` and `warp`, both `None` for an unwarped kernel, so `Likelihood.conditional` reports where the residuals are *claimed* to be stationary and the whiteness (family B) and localisation (family C) diagnostics run there — the degrees-of-freedom guard's last clause. **The guard is in the declaration, not in advice**: few knots, one shrinkage scale per warp with every knot variable normal about zero under it, non-centred by default (`u_k = s·z_k`, the same prior with the geometry NUTS wants) with `HierarchicalPrior` for the centred form, and the identity warp **bit-for-bit** the base kernel — the offset parameterisation `x ↦ x + δ(x)` with slopes `softplus(u)/softplus(0)` makes `δ` exactly zero at `u = 0`, so "the data must pay to leave the identity warp" is a statement about the base model itself and not about something near it. Monotonicity is structural (every slope is a softplus), so `QuasisepGP`'s ordering precondition cannot be violated and there is no rejection region in the posterior; the one float64 caveat, recorded in `likelihoods.md` §6 beside §7's midpoint re-referencing, is that below `u ≈ −37` a slope falls under machine epsilon and its segment goes flat to within a ULP. **Where a warp may sit**: a warp of a warp composes into one monotone map and is fine; a warp *under* a composite is refused by name on every backend (`refuse_warped_composite`), because two terms warped differently have two coordinates while the propagators come from one axis — `WarpedKernel(Sum(...))` is the representable order and `DenseGP` takes either. Knot placement is **fixed and explicit in the declaration**, with `quantile_knots(coordinate, count)` offered so quantile placement is a user act that ends up recorded rather than a hidden dependence on the data. No backend twin: a wrapping kernel has no arithmetic of its own and adopts its child's namespace, device and flags, as `Sum`, `Product` and `_Composite._adopt` already do — measured, torch's dense and quasisep solves agree to `0.0` and jax's to `5.7e-14` with both warps active, and the two backends' gradients with respect to all six knot variables agree to `6e-14`. Conformance extended (eight rows on every fixture) and green in all three environments. **Merged on Fable's second pass alone (the Codex quota ruling); Peter to confirm the no-twin choice.** |


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
partition via an **`eqx.partition` filter spec** (ruled 2026-09-01 —
preferred over `paramax.NonTrainable`, whose freezing depends on `unwrap()`
being called; see §2 and `docs/design/lowering.md` §6.2) — **never**
equinox static fields, which hash array contents into the JIT cache key
(see §7).

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
  side). **Implemented on the numpy path at W2.3** over celerite2, with
  ampere's own exact rank-2 Matérn-3/2 term and `conditional_loo` deferred
  (see the decision log above). **Decided:** the canonical flexible-likelihood kernel becomes
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
a possible later fallback for the reference path). *Clarified 2026-09-08
(W2.14, §2)*: "f ~ GP" is a statement the **scoring path** has to make good
on. The declaration is whitened — `z ~ N(0, 1)` i.i.d., with the covariance
supplied by `GPSolver.latent_transform` as `f = L(θ) z` — so
`GaussianProcessNoise.noise_params` applies that transform before the family
sees `noise.latent`. Without it the kernel hyperparameters enter the
likelihood nowhere at all and the latent block is white noise wearing a GP's
name.

### 4.5 Inference contracts
A fitting problem (model + instruments + datasets) exposes:

- `log_prob(params) -> float` and `log_likelihood` / `log_prior` split;
- `prior_transform(u)` for nested sampling;
- `simulate(params) -> data` for SBI;
- capability flags: `differentiable`, `batchable`, `device`, `backend`.
- **realisation** (added W2.13): a backend that offers gradients registers, at
  import, a factory turning a composed problem into its differentiable native
  form (`inference.md` §10a); gradient-based engines obtain it through
  `ampere.core.realise(problem)` and never import a backend.

The flags are properties of the *pieces*, not of the problem: models and
transformations declare `DIFFERENTIABLE`, `BATCHABLE`, `DEVICE` and `BACKEND`
as class attributes with conservative defaults (`False`, `False`, `"cpu"`,
`"reference"`), and the problem's `Capabilities` are aggregated from them.
`differentiable` and `batchable` are **conjunctive** — one black-box piece
withdraws the claim for the whole. `device` and `backend` are **identities**,
with no conservative answer available, so all the parts must agree and a
disagreement raises; the `capabilities=` argument is the override for a caller
who genuinely means one of them. `backend` (added W2.12) is a backend's one
name everywhere: the key `lowering.md` §12.8's registry uses, the id of its
conformance fixture, and what a run records as `ampere_backend` — which is
therefore a fact about the problem rather than a declaration by whoever
constructed the engine.

Any engine consuming only this contract works with every backend, including
legacy black-box models (via a thin adapter over the frozen legacy classes).

Also part of this contract:
- **Failure signalling**: external simulators crash and return NaNs; the
  contract defines the behaviour (`log_prob` → −inf with a recorded reason;
  `simulate` failures are flagged so SBI can reject-and-record rather than
  train on garbage). Ruled 2026-09-02 (`likelihoods.md` §17 Q1): the
  engine-facing path is **non-strict** — reachable evaluation failures (a
  non-positive-definite GP covariance, a non-positive Poisson rate) become
  −inf with a recorded reason, never an exception, because jax's traced hot
  loops cannot use exception control flow; a `strict` toggle lets the
  exception propagate instead, so a user who sees the recorded-reason
  warnings can re-run strict and get the raise at the offending draw.
- **RNG policy**: named seed handling that lowers to each backend's model
  (numpy Generators, torch Generators, jax PRNG keys), so runs are
  reproducible across backends.

### 4.6 Results, plotting, and the conformance suite
- **ArviZ's results format is the single results format.** (ArviZ 1.0
  retired the `InferenceData` *class* in favour of `xarray.DataTree`; the
  format and this decision are unchanged, and "InferenceData" elsewhere in
  this document names the format, not the class — corrected 2026-09-03,
  `results.md` §15 R3.) Every engine emits
  it; corner/trace/posterior-predictive plotting is written once against it
  in `ampere.results`; serialisation is netCDF. This is how the
  `mixins.py` monolith is retired and how sampler plotting parity stops
  regressing.
- **Every run stores per-draw `log_likelihood` and `log_prior`** (and run
  provenance: package versions, seeds, data hashes, spec hash, in the
  InferenceData attrs). Cheap now; it is the enabling requirement for
  population-level importance reweighting later (Design horizon, item b).
  (Wording reconciled at the freeze, resolving the ambiguity
  `likelihoods.md` §16 flagged: what every run stores is the **per-draw**
  value — the scalar joint log-likelihood per posterior draw, decomposed
  per dataset, which is what reweighting needs. The **per-observation**
  terms ArviZ's LOO/WAIC convention wants have no GP factorisation, so they
  are a *named* decomposition available on request —
  `Likelihood.pointwise_log_prob`, ruled 2026-09-03, `results.md` §6 —
  and are never stored by default.)
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
  native torch/jax equivalents, restoring differentiability for the most
  frequent cases. **Decided 2026-09-01: translation is opt-in, never
  silent** (see §2) — the default adapter always evaluates the user's actual
  astropy model; a native equivalent is substituted only through an explicit,
  backend-scoped request (sketch: `from_astropy()` on the backend, raising
  when the model is not fully translatable rather than silently degrading to
  black-box). A substituted implementation is not guaranteed numerically
  identical, and a model changing implementation without the user's
  knowledge is exactly the silent-downgrade class the architecture forbids.

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
is built. Phases 0–2 are decomposed into agent-sized work items in
`WORK_ITEMS.md` (Phase 2's at the W1.13 freeze); later phases are decomposed
when their prerequisites freeze.

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
- Branch harvest & archive: DONE (W0.7, merged 2026-09-01) — harvested
  content and per-branch provenance/caveats live in `docs/design/harvest/`
  (swyft TMNRE; `optim_only` optimisers; annotated `jax`-branch sketches;
  the copilot core scaffold that seeds W1.3), with archival
  recommendations for all 15 remote branches in
  `docs/design/harvest/branch_triage.md` awaiting Peter's approval.
  Phase 1+ design work should read the relevant harvest README before
  reinventing or reviving anything.

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
  decision log entries). — **Done** (W1.13, 2026-09-03; the freeze row in
  §2 is the record, and the `spec-v1.0` tag is created at merge).

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

**M2 REACHED 2026-09-08 (W2.10, `956e1ab`).** Reproduced on all three backends over 200 → 2 000 → 20 000 points with the science claim asserted as thresholds, posterior agreement across backends at stated tolerances, and the benchmark table as a CI artefact. The adversarial pass was Fable's (Codex quota-blocked; a retroactive `gpt-5.6-terra` pass is owed on the merged range, along with the earlier ones).

### Phase 3 — SBI layer (backend-spanning) — **complete 2026-09-10** (W3.0–W3.15; the landed summary is the paragraph after the two original bullets)
- One SBI module in `ampere.inference` consuming `simulate()` from any
  backend (including legacy black-box models): NPE/NLE/NRE via `sbi`;
  revive the swyft TMNRE implementation from
  `docs/design/harvest/swyft/` (note its README's latent bug: the lazy
  swyft-import guard never protected the `SwyftNetwork*` class
  definitions — fix during revival); embedding-network support carried
  over from the existing `infer/sbi.py` work. jax-native SBI
  (sbijax/flowjax) optional, later.
- A canonical **coordinate–value–mask tensor encoding** of containers for
  embedding networks, so set/attention-based embeddings can consume any
  modality — irregular sampling and missing data included. Fixed-size
  summaries remain the simple default for a single fitting problem, where
  the data layout is fixed anyway; the encoding is what makes amortisation
  across differently-sampled datasets possible.

**Phase 3 landed 2026-09-10 (W3.0–W3.15; the per-item detail is each status
row in `WORK_ITEMS.md` and the decision-log rows in §2 above — this
paragraph only reconciles the two bullets above against what shipped).**
`ampere.inference.SBIEngine` is the one module, over **sbi 0.27.0** (pinned
CPU torch, no pyro-ppl — the `sbi` pixi environment and extra, W3.2); NPE,
NLE and NRE are `sbi.inference`'s own trainers, unamortised truncated
marginal ratio estimation (TMNRE, Miller et al. 2021) is expressed through
those same `NRE` trainers and a `RestrictedPrior` subclass (W3.4) rather
than a revived swyft — **swyft is not used anywhere in ampere v2**, ruled by
Peter 2026-09-09 on W3.4's finding that swyft's `pytorch-lightning<=1.9.5`
pin cannot install beside the torch the `sbi` extra resolves to; the
harvested swyft material stays exactly what §5's Phase 0 bullet always
called it, an archived reference, not a dependency. Embedding-network
support is `EncodingLayout`'s masked set and transformer wrappers (W3.3,
W3.11), not a port of `infer/sbi.py`'s dict-based vocabulary, though that
vocabulary's shape is what a caller still writes. jax-native SBI did not
land and needed no revisiting (§6 below). The coordinate–value–mask
encoding is `ampere/core/encoding.py` and `docs/design/contracts/encoding.md`
(frozen at W3.3), exactly the second bullet's shape: one frozen `Layout`
per problem, hashed from the observed data alone, that both the flat and
the set/transformer embeddings read through `unpack`. Landed beyond either
bullet's wording, because the batching work Phase 3 needed turned out to be
a prerequisite rather than a detail: `FittingProblem.simulate_many` and its
executor protocol (`SerialExecutor`/`ThreadExecutor`/`ProcessExecutor`,
defaulting to **`forkserver`** on POSIX, Peter's ruling of 2026-09-09,
W3.1), per-chunk native batched prediction and sampling (W3.1 slice 2),
trained-artefact caching (`ArtefactStore`, W3.5) wired into the engine,
non-native ("foreign") parts opt-in (W3.8), and posterior calibration
(SBC/TARP, `ampere.results.calibration`, W3.6). `PROVENANCE_SCHEMA_VERSION`
is **6** (W3.12): every SBI run and training set carries `ampere_model_hash`
beside `ampere_spec_hash`, which is what the artefact cache and the
training-set append check key on, and — since W3.15 — `ampere_sbi_torch_seed`
records the seed that makes a run's network, training and posterior draws
reproducible bitwise from the problem's own seed.

### Phase 4 — Extensibility proof: one new modality end-to-end
- Implement **interferometric visibilities** (decided) through the whole
  stack to prove the composition design: Fourier sampling at (u,v) points as
  a Transformation, complex Gaussian likelihood, closure phases. It stresses
  the schema (complex data) and the transformation chain hardest; other
  modalities (astrometric time series, IFU cubes) then follow the template
  it establishes.
- The photometry + spectrum composition as a documented example (W4.11,
  ruled 2026-09-11): the simplest combined fit, written up for the docs
  before the interferometry template page cites it.
- Astropy interop adapter (`core/astropy_compat.py`, §4.7).
- **Kernel algebra and a public quasiseparable-term registry** (ruled by
  Peter 2026-09-09 as extensibility work: users must be able to compose a
  more expressive noise model). `Sum`/`Product` kernels with `KernelSpec`
  composition and celerite translation (a sum of quasiseparable terms is
  quasiseparable; a product is refused on the O(N) path); a damped
  periodic/SHO term for fringing; Matérn-1/2 and -5/2; a public
  `register_quasiseparable_term` so a user kernel can reach the O(N) path;
  a spectral-mixture kernel as a sum of SHO terms. Conformance rows per
  term against its dense closed form. Background and the sparsity-prior
  companion note: `docs/design/horizon_notes.md` §1 and its follow-up.

**Phase 4 as landed (2026-09-11 to 2026-09-13; W4.0–W4.11, closed at W4.8's
merge `c8bf2c0`).** The
composition design held for a modality nobody wrote the contracts for, and
then for a second one built by following the first's page. **Placement
(D1)**: a shipped observable's kind lives in `core/results_schema.py`, its
steps and models in `backends/{reference,torch,jax}/<observable>.py`, no
grouping namespace. **Interferometry** (W4.1, W4.3): `VisibilitySet` amended
to `(u, v, spectral_axis)` and `ClosurePhases` `(u1, v1, u2, v2, spectral_axis)`
with a canonical baseline ordering (D2, after Peter's chromatic-misspecification
question: a kernel sees axes only, and a missing band in a patch of sky is
sharp in wavelength and smooth in spatial frequency); `FourierSample` (a
direct DFT of a gridded `Image`, coverage from the observed container,
Nyquist requirements over the expanded coverage), `ClosurePhase`, the
bandwidth and time smearing steps (the first cross-kind uses of
`configure_from`), `Amplitude`, three image models and their analytic
visibility twins, inheriting native twins on both modern backends (the
declaration written once because a pixel scale written twice would alias
silently), `VonMisesFamily` implemented with `sample()` and native on both
backends, NUTS/VI/NPE all recovering the binary. **The circular complex GP**
(W4.2): one Cholesky with a two-column right-hand side, the `GPSolver`
contract widened to `(n, k)` residuals behind `STACKED_RESIDUALS`, the O(N)
path refused structurally, two native-path bugs (the coordinate stack and
the unbound kernel) fixed; SBC shows the flexible fit calibrated where the
rigid one is not. **Kernel algebra** (W4.5): seven families with exact
semiseparable representations (Matérn-5/2 rank 3, `likelihoods.md` §15.3
withdrawn), `Sum`/`Product`/`SpectralMixture`, the `axes=` selector with the
per-leaf unit rule, one public `register_quasiseparable_term` registry keyed
on the family with the array namespace as an argument. **The astropy
adapter** (W4.6, W4.7): `from_astropy` black-box on the reference backend
with bounds, fixed values and ties translated and no solid angle invented;
the opt-in native route with six curated models and compound `+ − * /`,
tied parameters and `|`/`&` refused by name. **The examples**: the
photometry + spectrum composition as the simplest combined fit (W4.11, far
infrared filters; the reference LSF operator cached by W4.0 to make it
tractable), the interferometry study with SBC-pinned coverage and a
chromatic arm reported informationally (W4.4), and **astrometry by the
template** (W4.9): a reflex orbit on two `TimeSeries` channels with no
change to core, four gaps in the template page found and handed to W4.8.
**CI** (W4.10): path-gated jobs whose skipped required checks report
success, typecheck split out of the backend legs, `actionlint` with
SHA-pinned actions. **Process**: every wave gated once on the merged master
in all four environments (all green); the token-economy rules of
`docs/orchestration.md` adopted mid-phase after the session limit killed
the first wave twice. **Carried to Phase 5 or the owed list**: a latent GP
on closure phases (W5.1); matrix-free exact GPs for the 3- and 5-axis
kernels; multi-axis point kinds refused by four of the six plots;
`NUTSEngine` tuning knobs; the encoding's positional axis alignment; the
Kronecker structure of dispersed data as the structured-solver route; the
grouping-namespace question, to be revisited after realistic usage; the terra reviews owed on W4.1, W4.2 and W4.5. **The documentation pass** (W4.8)
fixed the template for the gaps astrometry found, added the astropy and
kernel pages, and left nineteen *Amended W4.8* annotations across the
frozen documents — one of them a correction rather than a confirmation:
the interferometry sketch's Q2 (per-visibility frequency as an
`extra_coords` quantity) landed as a container axis instead.

### Phase 5 — Scale-out & advanced inference
- Approximate GP strategies for images/IFU behind the `GPSolver` /
  `NoiseModel` interface — **SVGP, SKI, Vecchia, HSGP, EFGP, chosen by
  measurement as celerite2 was** (*widened 2026-09-10 from
  `docs/design/horizon_notes.md` §2*): each is a `GPSolver` with
  `EXACT = False`, its approximation parameters in `provenance_config()`
  (~~recorded, never hashed~~ — **corrected at W5.4**: they are dataclass
  fields and enter the spec hash, because two runs at different `m` are not
  the same model; `provenance_config()` repeats them), a reduced-dimension
  `latent_transform`, and
  `conditional_loo` exact-in-the-approximation or refused by name. The rule:
  the first 2–3D solver landed is the one that wins a benchmark on the IFU
  sketch's cube at realistic N; HSGP is the cheapest to land and the most
  useful for NUTS, EFGP the one to reach for at image scale, Vecchia the
  one that handles rough processes without a rank limit. Two contract
  questions to settle when the phase opens: an approximation-aware
  conformance tolerance for `EXACT = False` solvers (compare to `DenseGP`
  with a tolerance that tightens with the basis size, and assert the
  convergence rather than a fixed number), and whether the latent size fixed
  at composition (`inference.md` limitation 17.4) may be the reduced rank
  rather than N — the contract permits it, but `simulate(observe=True)` and
  the latent-GP path must agree on the whitening. (GPJax is the natural
  provider on the jax side, GPyTorch on the torch side.)
- **Matrix-free exact GPs in 2+ dimensions** (noted by Peter 2026-09-11,
  not immediate): the interferometric kernels of Phase 4 live in 3 and 5
  axes where the quasiseparable tools do not apply, so exact GPs are
  dense there. `gpytorch` and `gpjax` avoid materialising the covariance
  by treating it as a linear operator and solving with conjugate
  gradients (MVM-based inference, with stochastic trace estimators for
  the log-determinant), which bounds memory at O(N) for exact GPs. This
  is a third solver strategy beside `DenseGP` and `QuasisepGP`, to be
  taken up when a Phase 4 or Phase 5 case actually exceeds the dense
  path's memory; nothing in the solver interface may preclude it.
- **Non-stationary flexible likelihood: input and amplitude warping as
  kernel wrappers preserving quasiseparability** (*added 2026-09-10 from
  `horizon_notes.md` §1 and its follow-up*): a `WarpedKernel(base,
  input_warp=, amplitude_warp=)` whose monotone input warp keeps
  `QuasisepGP`'s ordering precondition and whose amplitude warp is `D K D`
  with `D` diagonal, so both keep the O(N) exact solve; the warp knots are
  ordinary `Parameter`s (NUTS over them on torch/jax for free). Validated
  by an M2 extension with a "many lines / one band" scenario comparing
  stationary Matérn, warped Matérn and a sum of two kernels on bias,
  calibration and localisation. **The degrees-of-freedom guard is part of
  the item, not an afterthought**: few knots with hierarchical shrinkage to
  the identity warp; the whiteness and localisation diagnostics run on the
  *warped* residuals; SBC on injected misspecification. The same guard
  generalises to every sum of noise components W4.5's kernel algebra makes
  possible — a sparsity-inducing prior on component amplitudes (the
  regularised horseshoe, Piironen & Vehtari 2017, expressible today with
  `HierarchicalPrior`; spike-and-slab stays out, `lowering.md` §12 Q1),
  marginalised by the GP as usual, with the lowering rules offering the
  non-centred parameterisation NUTS wants (a `lowering.md` §3 note). Deep
  kernel learning waits for the multi-dimensional GP work above.
- **Joint noise over a tuple of channels — the linear model of
  coregionalisation** (*added 2026-09-10 from `horizon_notes.md` §3 and
  its follow-up; the limitation `likelihoods.md` §15 records*): a
  `NoiseModel` bound to several channels with `K = Σ_q B_q ⊗ k_q`, scoped
  first to the **shared-grid intrinsic model** `B ⊗ K_x`, which is exact
  and stays O(N) — diagonalise `B`, rotate the outputs, solve T scalar GPs
  with `QuasisepGP` — with `B` parameterised physically (a rotation for
  Q/U leakage) rather than a free T(T+1)/2; the general LMC and mismatched
  grids use the dense or a reduced-rank solver. Polarimetry the worked
  modality, the astrometric sketch the second test, interferometry the
  third (design horizon (h)). It is the first real use of
  `DatasetCollection.contributions` as something other than a sum, so it
  touches `inference.md` §4 and the results decomposition needs a
  `"joint"` entry — the same vocabulary widening `"mixed"` was; the
  diagnostics generalise per rotated output.
- **Amortisation over observation context** (*added 2026-09-10 from
  `horizon_notes.md` §4–5; the reserved hook is design horizon (i)*): a
  per-draw context — the σ-pattern drawn from a noise-realisation prior
  (scaled copies of the observed one, an archive of real error arrays, a
  parametric S/N model), the grid, the instrument settings (resolution,
  filter set, exposure) — as an input to `simulate_many`, recorded on the
  `Simulation`, passed to `sample` in place of the container's σ, with the
  chain re-negotiated per context (grouped by `chunk_size` to keep it
  cheap) and the context surfaced to the embedding as row features (σ,
  `log σ`, whitened values) and a per-set conditioning vector (FiLM, or
  extra tokens). SBC per observation (W3.6) is the check that the context
  prior covered the observation at hand. This is the amortisation the
  population use case of `docs/design/inference_extensions_memo.md` §8.1
  depends on.
- Hierarchical/population inference: implement the container + hyperprior
  design from Phase 1 (plates in pyro/numpyro).
- **RHMF exploratory trial** (Peter's ratification note, 2026-09-08, on the
  W2.7 deferral row): early testing of the pre-fit robust-factorisation
  family (`diagnostics.md` §2) against a pinned RHMF commit or a ≥0.1
  release, behind a non-default extra, before the deferred
  `ampere.diagnostics` namespace lands — the §2.2 adoptability re-check is
  the gate, and the trial informs it rather than waiting for it.
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
- **Composition tutorials that the capability already supports but no
  example shows** (Peter, 2026-09-09): a worked, runnable example of a
  problem with several observations — photometry plus one or more spectra
  on one model, `SyntheticPhotometry` beside `Resample`/`LSFConvolution`
  — and how calibration uncertainty is represented (a `CalibrationScale`
  step with a prior as the nuisance parameter, per spectrum or shared via a
  tie; the flexible likelihood as the complement for what calibration does
  not explain). The capability is in `ampere.core` and the reference
  backend and is exercised by `docs/design/modalities/spectrum_photometry.md`
  and the tests; what is missing is the user-facing example, best written
  when the legacy tutorials and `minimal_working_example*.py` are converted
  here. Not needed before then, but a `tests/examples` smoke test of the
  composition may be worth adding earlier so the capability cannot regress
  unnoticed.

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
  check with warnings-as-errors. **Landed at W2.11, minus the `-W`** — the
  build check is a job, warnings-as-errors waits for Phase 6 to retire the
  legacy docstrings that produce them; the decision-log row says why.
- **Phase 2**: backend matrix — separate torch, jax, and no-extras jobs
  (the last catches lazy-import breakage, i.e. `import ampere` must never
  require torch or jax); dependency caching; benchmark suite
  (pytest-benchmark or asv) with results tracked as CI artefacts. GPU tests
  are nightly/manual — hosted runners are CPU-only. **Landed at W2.11**: one
  job per environment (`dev`, `torch`, `jax`) with `fail-fast: false` so a
  broken backend fails only its own leg; `minimal-install` is the no-extras
  leg and stays a bare `pip install -e .`; `setup-pixi`'s cache covers all
  of them; `pixi run bench` (pytest-benchmark) uploads `benchmark.json` per
  environment. The torch environments resolve torch from PyTorch's CPU
  index, which is what makes "hosted runners are CPU-only" cheap as well as
  true.
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
  built entirely on stored InferenceData files. **Constraint added
  2026-09-10 (Peter, on the inference-extensions memo §8.1)**: population
  inference over amortised SBI at 10⁴–10⁹ sources will need a columnar,
  partitioned store rather than one file per source; its design waits for
  the population work, but nothing landed before then may exclude it — no
  results or training-set format may assume one file per source, and the
  per-draw columns reweighting needs (`log_prior`, the proposal's own
  log-density, a log-ratio where one exists) stay separable from per-run
  provenance.
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
- **(e) Model comparison and averaging** (added 2026-09-10, approved by
  Peter on `docs/design/inference_extensions_memo.md` §9): evidences from
  nested sampling, SMC and Laplace, the learnt harmonic-mean estimator on
  archived draws, bridge sampling, Savage–Dickey ratios, and Bayesian model
  averaging over runs sharing a data hash. Hook — **decided**: an
  engine-neutral `ampere_log_evidence`/`ampere_log_evidence_err`/
  `ampere_evidence_method` triple in the root attrs (W5.0), beside the
  per-draw `log_prior`/`log_likelihood` every run already stores; LOO/WAIC
  already ride on the pointwise group.
- **(f) Approximate-inference correction** (added 2026-09-10, same
  approval): importance correction of VI, Laplace, Pathfinder and SBI runs
  toward the true posterior, population reweighting at scale (b), and
  retroactive density emulation on archives. Hook — **decided**: every
  engine whose draws are not from the target stores the proposal's own
  log-density per draw, and `ampere_approximation` names the family in the
  root attrs (`"none"` for an exact sampler) — W5.0; `SBIEngine` stores it
  today, `VIEngine` will.
- **(g) Engines as a registry, with uniform cost accounting** (added
  2026-09-10, same approval): third-party and extra-gated engines
  discoverable by name, stamped into provenance by shared code, and
  enumerated by an engine battery that rank-calibrates every engine
  through `ampere.results.calibration.sbc`. Hook: a registry shaped like
  the realisation registry, and one cost record every engine writes
  (evaluations, simulations, surrogate calls — per fidelity where (a)
  applies). Each new sampler lives behind its own extra (ruled
  2026-09-10): the base install stays quick to start with, and a user
  upgrades for a specific problem.
- **(h) Vector observables with correlated misspecification** (added
  2026-09-10 from `horizon_notes.md` §3's follow-up): one wrong sky model
  produces coherent errors in every quantity derived from it — Stokes
  components, RA/Dec, visibility amplitude and phase, multi-band fluxes —
  so a coregionalised misspecification model is the *default* shape for
  vector data, not a special case for instrumental leakage. Hook:
  `results_schema.md` §15.2's one-value-array-per-container rule stands,
  and **Phase 4's interferometry item leaves the channel pairing between
  visibilities and closure phases visible** (W4.1) so Phase 5's joint
  noise model can bind to it; `contributions` and the `"joint"`
  decomposition entry are the likelihood-side hooks.
- **(i) Amortisation over observation context** (added 2026-09-10 from
  `horizon_notes.md` §4–5 and its "Consequence for Phase 3"): σ-pattern,
  grid and instrument settings as per-draw simulator inputs from a context
  prior, so an amortised posterior is valid across instruments and error
  bars, not only across θ. Hook — **reserved at W3.1/W3.2/W3.3 and
  confirmed by Peter 2026-09-09**: `simulate_many`/`SBIEngine` carry a
  `context=` slot (accepting `None` today, recorded in provenance), and the
  W3.3 encoding carries the uncertainty column by default. The machinery
  itself is the Phase 5 bullet above; nothing before it may fix a call
  signature that leaves no room for a context.

Dependencies: 0 → 1 → 2 → {3, 4} → 5 → 6, with 3 and 4 parallelisable and
the CI/CD workstream running alongside every phase.

## 6. Deferred implementation choices

Settled at the start of the phase that needs them, not now:

- ~~**Torch GP solver library** (Phase 2): GPyTorch structured solvers vs
  celerite2's experimental torch interface for `QuasisepGP`. Evaluate both
  against the conformance suite; GPyTorch is the safer bet for the broader
  Phase 5 strategies (SVGP/SKI).~~ **Settled 2026-09-07 by measurement**
  (W2.4 slice 2; the decision-log row carries the table): celerite2's
  compiled kernels under `torch.autograd`, at zero dependency cost and
  bit-identical to the reference solver, against a GPyTorch that has no
  quasiseparable operator at all. The parenthetical hope — "GPyTorch is the
  safer bet for the broader Phase 5 strategies (SVGP/SKI)" — is **not**
  discharged and is not contradicted: GPyTorch remains the obvious provider
  of Phase 5's *approximate* strategies (`InducingPointGP`,
  `StructuredGridGP`), which are different slots on the same `GPSolver`
  interface. The two choices were never one choice; treating them as one is
  what the measurement corrected.
- **Jax GP libraries** (Phase 2 gate check): resolved policy is
  *both-where-strongest behind the strategy interface*, not a single
  library. ~~`QuasisepGP` comes from tinygp's `QuasisepSolver` or
  celerite2.jax (verify at Phase 2 start which is better maintained)~~ —
  **settled at W2.5 slice 2 (2026-09-07): `QuasisepGP` is celerite2.jax.**
  Both were implemented against ampere's exact rank-2 Matérn-3/2
  representation and measured; the decision-log row above carries the
  table. Maintenance did not decide it and accuracy did not decide it
  (tinygp is the more accurate, by three to five orders, and both are far
  inside `tolerances.cross_solver`). **Scaling decided it**: celerite2's
  XLA custom calls are linear at ~0.25 µs/point to 10⁶ points, while
  tinygp's jax-native `lax.scan` recursions are quadratic in practice on
  XLA's CPU backend — a factor of 200 at 10⁴ points and widening. The
  costs taken with it, all recorded: a deferred import (celerite2.jax
  flips `jax_enable_x64` at import, which §10.2(a) forbids ampere from
  doing), `BATCHABLE = False` for that solver alone (no `vmap` batching
  rule), and — at slice 2 — `conditional_loo` still refused (no O(N) route
  to the inverse diagonal through `celerite2.jax`'s `GaussianProcess`
  wrapper, which tinygp *does* offer). **Slice 3 removed the third cost**
  without disturbing the choice: `celerite2.jax.ops` is public and is the
  same set of compiled kernels torch reaches through `celerite2.backprop`,
  so the O(N) recursion runs on jax too (§2's own row). The strategy
  interface is what would make a later swap a one-file change.

  **Both tracks converged on celerite2, by different routes and on different
  evidence, and that is worth recording once here rather than twice above.**
  These were always two questions — the bullet above is torch's — and neither
  answer was assumed from the other: torch measured celerite2's compiled
  kernels under `torch.autograd` against GPyTorch (which turned out to have no
  quasiseparable operator at all) and jax measured `celerite2.jax` against
  tinygp (which turned out to be quadratic in practice on XLA's CPU backend).
  The convergence has one concrete consequence and one caution. The
  consequence: all three implemented `QuasisepGP`s — numpy, torch and jax —
  now factorise through **one library's** arithmetic over one representation
  ampere wrote three times, so the conformance suite's cross-solver and
  cross-backend rows test ampere's lowering rather than three third-party
  solvers. The caution is the same fact read the other way: celerite2 is now a
  single point of failure for the project's distinguishing feature, and the
  strategy interface is what keeps that a risk rather than a lock-in. **Since W2.5 slice 3 the two tracks
  no longer disagree about anything**: both supply `conditional_loo` in O(N),
  by the same recursion, because both call celerite2's compiled kernels
  directly and so already hold the factorisation it walks — torch through
  `celerite2.backprop`, jax through `celerite2.jax.ops`. `ampere.core`'s numpy
  solver is the one that still refuses, and it refuses for the reason it
  always gave: it composes celerite2's `GaussianProcess` rather than calling
  the kernels, so the coupling would be new there.
  `tests/conformance/README.md` §6 states that as one debt entry, and the
  conformance row states both branches so a backend remains free to choose.

  GPJax — attractive for its Equinox/Paramax foundations and variational
  machinery — remains the candidate provider for Phase 5 sparse/variational
  strategies, but **verify it actually offers a quasiseparable/state-space
  O(N) exact solver before letting it displace this one** (believed absent
  as of early 2026). Ampere's own §4.1 layer remains the only user-facing
  parameter interface regardless of provider, for cross-backend parity;
  Paramax is a lowering mechanism, not a user API.
- **SBI package set beyond `sbi` + ~~swyft~~** (Phase 3): **half settled
  2026-09-09 (Peter's ruling on W3.4's finding, recorded in the TMNRE
  decision-log row above): no swyft.** TMNRE is expressed through `sbi`
  0.27's own `NRE` trainers and a `RestrictedPrior` subclass instead of a
  revived swyft implementation — swyft's `pytorch-lightning<=1.9.5` pin
  cannot install beside the torch the `sbi` extra resolves to, and nothing
  in Phase 3 needed it once that was known. **jax-native SBI (sbijax/flowjax)
  remains open, unchanged**: if and when maturity warrants — nothing in
  Phase 3 needed it, since `simulate_many`'s native `simulate_batched` is the
  jax half that would matter if it ever did.
- **Embedding-network study** (Phase 3, deferred at W3.3/W3.11's review;
  Peter's rider, 2026-09-09: "defaults, not findings"). The set and
  transformer embeddings' default output width
  (`max(2·free_size, 32)`, W3.11) and readout (mask-weighted mean pooling
  over every retained token, replacing `sbi`'s own last-token read, W3.11)
  are principled defaults chosen to make the wrappers correct, not the
  product of a study of which width and readout suit which combination of
  data, model and structure — that guidance for a user choosing between
  `layout="flat"`, `"set"` and `"transformer"` is future work. Folds in the
  **truncation-epsilon (ε) study** carried from W3.4 (`docs/development.md`'s
  "Questions collected for Peter" §3): whether TMNRE's default
  `truncation_epsilon = 1e-4` (a box of roughly ±4σ that stops shrinking
  after round 2 on the worked example) generalises, or needs per-problem
  tuning guidance the way the width does. *(Added 2026-09-10 from
  `horizon_notes.md`'s embedding follow-up)*: the first *architecture*
  experiment worth running now that W3.6's calibration exists is a
  ConvCNP-style encoder (SetConv onto a fine internal grid, then a CNN —
  translation-equivariant along the coordinate and discretisation-invariant)
  as an `examples/sbi/` script through W3.2's `embedding=` slot, not a
  contract change; neural-operator branch networks are the same idea shared
  with design horizon (c)'s emulators.
- ~~**Benchmark harness** (Phase 2): pytest-benchmark vs asv.~~ **Settled
  2026-09-08 (W2.11): pytest-benchmark.** The decision-log row above carries
  the reasoning; in one line, asv owns its own environments and its own
  history store, and this project's environments are pixi's and its history
  is the PR. `tests/benchmarks/` is the suite, `pixi run bench` the one
  command, and `benchmark.json` the artefact CI attaches per environment.

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
