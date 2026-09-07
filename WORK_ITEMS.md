# Ampere v2 — Work Items, Phases 0–2

Companion to `DEVELOPMENT_PLAN.md` (the source of truth for architecture and
decisions — read it first). Each item below is sized for delegation to a
development agent. The Phase 2 items were written at the spec freeze
(W1.13), against `spec-v1.0`; Phase 3+ items are written once their
prerequisites freeze.

Sizes: **S** ≈ half an agent session, **M** ≈ one session, **L** ≈ 1–2
sessions.

## Working agreement (all items)

- One item = one feature branch = one PR against `master`; Peter reviews and
  merges. No direct pushes to `master`; no tag pushes or branch deletions
  without explicit approval.
- CI must be green before review; from W1.10 onwards, contract changes must
  keep the conformance suite green or update it in the same PR with
  justification.
- Legacy modules (`ampere/data`, `ampere/models`, `ampere/infer`) are
  frozen: touch them only where an item explicitly says so.
- Any change to a §4 contract after the spec freeze requires a decision-log
  entry in `DEVELOPMENT_PLAN.md` in the same PR.
- British English in all documentation and prose.
- Update this file's status column in the PR that completes an item.

## Phase 0 — Safety net & hygiene

### W0.1 — Commit the pending MAP-plot fix [S]
Commit the uncommitted `ampere/infer/mixins.py` change (MAP plotting via the
dict-based model return) on `master`.
**Accept:** clean `git status` for tracked files; the fix's example still runs.

### W0.2 — Repository hygiene [S]
Extend `.gitignore` for locally generated artefacts (run logs, PNGs, pickles
under `examples/`, `miniforge.sh`, egg-info). Do **not** delete the user's
untracked paper-revision files.
**Accept:** `git status` shows no noise from a fresh example run.

### W0.3 — Packaging floor [S]
`requires-python >= 3.11`; fix the invalid GPL trove classifier; add extras
skeleton (`torch`, `jax`, `sbi`, `dev` updated; `all`); switch to package
auto-discovery so new namespaces don't need manual listing; add a dev
environment file.
**Accept:** `pip install -e ".[dev]"` clean in a fresh Python 3.11 env;
`python -c "import ampere"` works (jointly with W0.5).

### W0.4 — Characterisation test harness [L]
Golden-output pytest suite wrapping the `examples/minimal_working_example*`
flows (emcee, dynesty, zeus, SBI) with fixed seeds and small iteration
budgets: assert posterior summary statistics within tolerances, output
shapes, and successful plot/pickle generation. These tests define "legacy
still works" for every later PR. Tolerances: loose enough to survive
dependency-version noise; tight enough to catch a broken likelihood —
document the choice per test.
**Depends:** W0.3. **Accept:** `pytest tests/characterisation` green twice in
a row locally (determinism) in < ~10 min.

### W0.5 — Broken-module quarantine and dependency pins (issues #74–77 + W0.4 findings) [M]
Fix the astropy `BlackBody` API usage (#77); fix `CCMExtinctionLaw` syntax
errors (#76); replace the abandoned `extinction` dependency with
`dust_extinction` or inline curves (#75); fix the `extinctionModels` import
error (#74). Make heavy/optional imports lazy so no broken or missing
optional corner can block `import ampere`.

Additional scope from W0.4's findings (2026-09-01):
- (The `pyphot<2` and `sbi<0.28` pins were applied directly on master
  2026-09-01; forward-migration to current versions is W0.9.)
- Fix or quarantine-with-informative-error the two SBI defects W0.4
  isolated: `SBI_SNPE.postProcess()` crashes on sbi 0.27's batched
  `posterior.map()` shapes (`mixins.py:542`), and
  `check_prior_normalisation=False` raises `AttributeError` because
  `_prior_is_normalised` is never set. Consider an upper pin on `sbi` if
  fixing is disproportionate.
**Accept:** `import ampere` succeeds in a minimal-deps env; a test imports
every non-quarantined module; quarantined corners raise informative errors
on use, not on import; **`pytest tests/characterisation` in a fresh
`.[dev,zeus,sbi]` env yields real passes, not pyphot xfails**; the W0.4
xfail markers that the fixes obsolete are removed in the same PR.

### W0.6 — CI + pre-commit [M]
GitHub Actions PR gate: ruff (lint + format), pyrefly (scoped to
`ampere/core`, `ampere/backends`, `ampere/inference`, `ampere/results` —
initially empty scope is fine), pytest + coverage, Python 3.11–3.13 matrix.
Matching pre-commit config.
**Depends:** W0.4, W0.5; coordinate with W0.8 — prefer building the CI
environment setup on pixi directly rather than migrating it afterwards.
**Accept:** green runs on a PR; a deliberately broken test/lint fails the
gate.

### W0.8 — pixi migration [M]
Move environment management to pixi while keeping pyproject.toml the single
source of truth for dependencies: `[tool.pixi]` tables with features/
environments mirroring the extras (`dev`, `torch`, `jax`, `sbi`), a
committed `pixi.lock`, tasks for the common commands (test, lint,
typecheck, docs), and CI switched to `setup-pixi` with caching. Retire
`environment.yml`; update AGENTS.md and `docs/development.md` environment
instructions.
**Depends:** W0.3 merged. Schedule just before or together with W0.6 so CI
is built on pixi once, not twice.
**Accept:** from a clean clone with only pixi installed, `pixi run test`
(and lint/typecheck tasks) work; lockfile committed; CI green via pixi;
`environment.yml` removed; docs updated.

### W0.9 — Dependency forward-migration: pyphot ≥2 and current sbi [L]
The `pyphot<2` and `sbi<0.28` pins are temporary; both ecosystems should
be migrated forwards, not frozen out — we want current pyphot and the
latest sbi inference algorithms available.
- **pyphot**: write a small compat layer (e.g. `ampere/utils/
  pyphot_compat.py`) supporting both pyphot 1.x and ≥2 unit handling;
  migrate `ampere/data/photometry.py` onto it (a permitted critical fix to
  frozen legacy — a breaking dependency counts); lift the pin. MUST land
  before Phase 2's synthetic-photometry Transformation is written, so new
  code targets the pyphot ≥2 API from the start and never inherits the
  legacy idiom.
- **sbi**: update `ampere/infer/sbi.py` for the current sbi API (batched
  `posterior.map()` shapes in `postProcess`, prior-validation behaviour);
  lift the pin. Natural timing: with, or just before, Phase 3, which
  targets latest sbi for the new SBI layer anyway.
- Characterisation suite must pass for real against the migrated versions
  in the same PR(s); remove any obsoleted xfail markers.
**Depends:** W0.4 merged (the suite is the safety net for exactly this
migration); the sbi half also benefits from W0.5's quarantine decisions.
**Accept:** fresh-env suite green with unpinned pyphot ≥2 (and current sbi
for the sbi half); pins removed from pyproject.toml; compat layer
documented.

### W0.10 — CI coverage of the Phase 1 suites [M]
Authorised by Peter 2026-09-03. `ci.yml` never runs `tests/core`,
`tests/results` or `tests/conformance` — the py3.11 core-import break that
survived until W1.10's review (fixed `a8e6f54`) is the proof it bites.
- Run the three Phase-1 suites in CI per PR; different schedules for
  different parts are acceptable if runtime demands it (e.g. the full
  cross-backend battery on a schedule, the fast suites per PR).
- Fix finding (c): `tests/core/test_likelihood.py` registers throwaway
  likelihood families module-globally and never deregisters them, so the
  suites fail in one pytest process — a registry-snapshot fixture around
  the registering tests.
- Fix finding (b): `pixi run <task>` in the *default* environment cannot
  import ampere from a worktree — default the tasks to `dev` or pin the
  default environment's Python.
- Policy (Peter, 2026-09-03): py3.11 stays in the matrix but may be
  dropped if it proves problematic — with 3.15 imminent it is not worth
  fighting for.
**Accept:** a deliberately broken conformance row fails a PR; the three
suites green in a single pytest invocation; pixi tasks work from a
worktree.

### W0.7 — Branch harvest & archive proposal [M]
Extract into `docs/design/harvest/`: the swyft TMNRE diff, the
`optim_only` optimiser module, the `jax`-branch design sketches (annotated:
aspirational, does not import), and the copilot-branch core scaffold +
IMPLEMENTATION_PLAN.md. Produce a table of all remote branches with a
recommendation (archive-tag / delete / keep) for Peter's approval — do not
execute archival.
**Accept:** harvest files in place with one-paragraph provenance notes each;
branch table in the PR description.

## Phase 1 — Core contracts

### W1.1 — Prior-art memo [M]
Study bilby (likelihood/prior/sampler decoupling), gammapy (Datasets
container), 3ML (per-instrument plugin likelihoods), Starfish
(misspecification GPs for spectra), and RHMF/Robusta-HMF (arXiv:2607.08081,
for §4.8). Write `docs/design/prior_art.md`: for each, the two or three
interface lessons ampere should copy or avoid, with pointers.
**Accept:** memo reviewed by Peter; lessons cross-referenced from later specs.

### W1.2 — Architecture spec [M]
`docs/design/architecture.md`: expand DEVELOPMENT_PLAN §3 — capability
ladder, reference-backend rationale, namespace layout, extras/lazy-import
policy, dtype/device/precision policy, the parameters-vs-buffers model
contract, and the functional-data stance (coordinate-indexed containers).
**Depends:** W1.1. **Accept:** review pass by Peter.

### W1.3 — Parameter & prior contract [L]
`docs/design/contracts/parameters.md` + `ampere/core/parameter.py` + unit
tests. Named parameters; neutral (scipy-style) prior declarations;
transforms to unconstrained space; frozen/fixed; tying/sharing across
models and datasets; plate-aware groups; units; buffer declaration.
Start from the copilot-branch scaffold (via W0.7 harvest).
**Depends:** W1.2. **Accept:** typed (pyrefly clean); worked examples in the
spec run as doctests/tests; round-trip prior serialisation.

### W1.4 — ModelResult schema contract [L]
`docs/design/contracts/results_schema.md` + `ampere/core/results_schema.py`
+ tests. Named channels; typed kinds (`Spectrum`, `PhotometricPoints`,
`Image`, `Cube`, `TimeSeries`, `VisibilitySet`); coordinate-indexed
function-sample containers with uncertainties; first-class masks; units;
fidelity tags; default-channel sugar for single-output models.
**Depends:** W1.2. **Accept:** typed; the low-res-SED + CO-windows example
from the plan expressed and validated in tests.

### W1.5 — Transformation & Instrument contract [L]
`docs/design/contracts/transformations.md` + `ampere/core/transform.py` +
tests. Transformation ABC (with nuisance parameters), Instrument as chain,
channel binding with kind checks, and the optional requirements-negotiation
protocol (published requirements → one-off compilation → hot loop). Include
the out-of-tree extension example (a user-defined Transformation imported
as if third-party).
**Depends:** W1.3, W1.4. **Accept:** typed; simple path (fixed grid, no
negotiation) demonstrably trivial; out-of-tree extension test passes.

### W1.6 — Likelihood & NoiseModel contract [L]
`docs/design/contracts/likelihoods.md` + `ampere/core/likelihood.py` +
tests. Family registry (`log_prob(predicted, observed, noise_params)`);
analytic-marginalisation vs latent-variable declaration; NoiseModel
strategy interface (DenseGP / QuasisepGP / future approximate slots);
Matérn-class kernel specification; mask propagation; censoring interface
(#11). Include a numpy DenseGP implementation with Matérn-3/2 as the
correctness anchor for later solvers.
**Depends:** W1.4. **Accept:** typed; DenseGP validated against analytic
marginal-likelihood cases; latent-path declaration exercised by a stub
Poisson family.

### W1.7 — Dataset, FittingProblem & inference contracts [L]
`docs/design/contracts/inference.md` + `ampere/core/dataset.py` +
fitting-problem object + tests. Dataset (data + instrument + likelihood),
DatasetCollection (joint log-likelihood, shared parameters, hyperprior
extension point); the engine-facing surface: `log_prob`, `prior_transform`,
`simulate`, capability flags, failure signalling (−inf + recorded reason;
flagged failed simulations), RNG/seed policy.
**Depends:** W1.3–W1.6. **Accept:** typed; a toy two-dataset joint problem
with a tied parameter works end-to-end against a stub model.

### W1.8 — Results contract [M]
`docs/design/contracts/results.md` + `ampere/results` skeleton.
InferenceData emission helpers (always storing per-sample `log_likelihood`
and `log_prior`, plus provenance attrs: versions, seeds, data/spec hashes);
the plotting API surface (corner, trace, posterior-predictive,
GP-localisation) as stubs against InferenceData.
**Depends:** W1.7. **Accept:** typed; emission round-trips through netCDF.

### W1.9 — Lowering spec [M]
`docs/design/lowering.md`: distribution mapping table (scipy ↔
torch.distributions ↔ numpyro), transform conventions, module/pytree
conventions (equinox + Paramax candidate), parameter-vs-buffer lowering
(torch `register_buffer`; jax partition filters — never static fields),
RNG lowering, dtype/device policy (float64 for likelihood linear algebra,
jax x64 flag). Document, not code — the backends implement it in Phase 2.
**Depends:** W1.3. **Accept:** every §4.1 declaration form has a lowering
row for each backend, or an explicit "unsupported" entry.

### W1.10 — Conformance suite skeleton [L]
`tests/conformance/`: the parametrised-over-backends battery every backend
must pass — prior round-trips, schema validation, log_prob agreement
against analytic oracles and the reference backend, DenseGP↔QuasisepGP
agreement on quasiseparable kernels, mask/censoring behaviour,
InferenceData emission. Runs against the core stubs now; backends plug in
during Phase 2.
**Depends:** W1.3–W1.8. **Accept:** suite runs (and passes) against a
minimal in-repo stub backend; adding a backend requires only a fixture.

### W1.11 — Modality design sketches [L, splittable across agents]
`docs/design/modalities/`: one short sketch each for (a) joint
spectrum+photometry (the v1 slice), (b) the multi-channel low-res +
CO-windows case, (c) interferometric visibilities + closure phases,
(d) astrometric time series, (e) IFU cube, (f) hierarchical population —
each showing the model → channels → instrument → likelihood composition in
contract terms, naming any interface gap found.
**Depends:** drafts of W1.4–W1.6. **Accept:** no unresolved interface gaps,
or gaps fed back as spec amendments before the freeze.

### W1.12 — Diagnostics design spec [M]
`docs/design/contracts/diagnostics.md` per plan §4.8: RHMF-style pre-fit
screening (assess Robusta-HMF's adoptability: licence, maturity, API),
post-fit residual whiteness tests, GP-localisation outputs; where each
lives (`ampere.results` vs a `diagnostics` module) and what lands in
Phase 2.
**Depends:** W1.1. **Accept:** review pass by Peter.

### W1.13 — Spec assembly & freeze [M]
Cross-review all specs for contradictions; resolve W1.11 amendments;
Peter's review pass; tag `spec-v1.0`; record the freeze in
DEVELOPMENT_PLAN's decision log; write the Phase 2 work-item breakdown
against the frozen spec.
**Depends:** everything above. **Accept:** tag exists; Phase 2 items added
to this file.

## Phase 2 — Twin modern backends, lockstep (written at the freeze, W1.13)

Everything below is against the **frozen spec** (`spec-v1.0`): any change a
Phase 2 item needs to a §4 contract requires a decision-log entry in
`DEVELOPMENT_PLAN.md` in the same PR (working agreement), and the
conformance suite (`tests/conformance/`) is the lockstep mechanism — a
feature is done only when every registered backend passes it. Dispatch per
`docs/orchestration.md`: **Opus per backend track** (W2.4, W2.5), Sonnet
for the well-specified rest, Fable review at merge, adversarial
cross-model review at milestone M2. The v1 slice stays spectra +
photometry end-to-end; no new modality before Phase 4.

### W2.1 — `ampere.backends.reference`: the numpy backend and base-install slice [L]
The conformance oracle becomes a real package (`architecture.md` §1/§3):
native models (blackbody, modified blackbody, power law — plan §5's list),
the standard instrument steps (spectral resampling, LSF convolution,
calibration scale; **synthetic photometry only after W0.9's pyphot
migration**, so new code targets the ≥2 API from the first line), and
**`FractionalModelNoise`** — the prediction-aware noise model
`likelihoods.md` §5 names (X-1), with `f` an ordinary fitted parameter —
registered so the X-1 conformance rows run against the shipped class
instead of the in-repo double. Register the backend as a conformance
fixture (one registry line, per W1.10's acceptance).
Also lands here (ruled 2026-09-03 at the freeze's escalations; both are
post-freeze §4 additions, so this item's PR carries their decision-log
entries): the opt-in `describe()` hook on `Parameterised` folded into
`model_fingerprint`, plus the derived backend-neutral model identity
(offer, never serve — `results.md` §14); and `Axis.locate(values)`
matching within `COORDINATE_RTOL` (`spectrum_photometry.md` Gap 1),
which this item's resampling/photometry steps are the first to consume.
**Depends:** spec-v1.0 merged; W0.9 for the photometry step only (the rest
must not wait on it).
**Accept:** conformance suite green with the new backend registered; the
X-1 σ_eff row exercises `FractionalModelNoise` itself; `import ampere`
still requires no extras; a no-extras environment composes and scores a
blackbody + resampling + flexible-GP problem end-to-end.

### W2.2 — Engine drivers: emcee, dynesty, zeus in `ampere.inference` [L]
The gradient-free engines, written once against §4.5's surface and nothing
else (`inference.md` §10) — they must work unchanged with every backend.
Every run emits the ArviZ `DataTree` through `ampere.results` (per-draw
`log_likelihood`/`log_prior`, provenance attrs, netCDF round trip), which
is also the moment **arviz + a netCDF engine join the base install**
(`results.md` §15 R1's ruling: promoted with the engine drivers, when a
user can first emit a run) — lift them from the extras in the same PR.
Non-strict failure signalling consumed as declared (−inf + recorded
reason; `failure_summary` surfaced to the user); `check_engine` called
with `observed=`.
**Depends:** W2.1 (a real backend to drive).
**Accept:** the toy two-dataset joint problem samples end-to-end on all
three engines; posterior summaries within tolerance of each other on a
known problem; the stored run round-trips netCDF with hashes intact;
arviz imports from the base install.

### W2.3 — `QuasisepGP` on the reference path (celerite2) [M]
Fill the declared O(N) solver slot with celerite2's numpy interface and
un-skip the `DenseGP`↔`QuasisepGP` conformance rows (the debt W1.10
recorded by name). Includes `GPSolver.conditional_loo`'s O(N) recursion
for the leave-one-out terms — or, if that recursion is deferred, a
decision-log-recorded deferral with the refusal row kept live. Matérn-3/2
exactness is the whole point: agreement at `tolerances.cross_solver`.
**Depends:** W2.1.
**Accept:** the skipped solver rows run and pass for the reference
backend; the empty-slot refusal row flips to the implemented branch;
10³–10⁵-point scaling demonstrated (wall-clock, not asymptotics claimed).

### W2.4 — The torch backend track [L, multiple sessions; Opus]
`ampere.backends.torch` against the frozen spec, gated by the conformance
suite throughout. Scope per plan §5 and `lowering.md`: distribution and
bijection lowering (the §3 tables; `LoweringError` for rule-2 refusals),
buffers via `register_buffer` (§7), RNG via `substream` →
`torch.Generator.manual_seed` (§9), float64 policy for likelihood linear
algebra (§10), native models (the W2.1 trio), GP solvers (`DenseGP`
natively; `QuasisepGP` via GPyTorch or celerite2-torch — evaluate both
against the conformance suite, plan §6), NUTS + VI via pyro, capability
flags declared on the subclasses. **The icdf-fallback contract**
(`lowering.md` §3.6): a family without native `icdf` takes the reference
fallback with a once-per-run warning naming families and backend, and
`FittingProblem.strict=True` raises `LoweringError` instead.
**Depends:** W2.1, W2.2, W2.6; lockstep with W2.5 — neither track merges a
feature the other cannot pass the suite on without a recorded reason.
**Accept:** full conformance battery green for the torch fixture
(cross-backend rows against reference included); NUTS recovers the toy
joint problem's posterior; the icdf warning and strict raise are tested.

### W2.5 — The jax backend track [L, multiple sessions; Opus]
`ampere.backends.jax`, mirror of W2.4: numpyro distributions and
`biject_to` (never `transform_to` — `lowering.md` §2), buffers and fixed
parameters via `eqx.partition` filter specs (never static fields — plan
§7), `configure_x64()` guard-and-raise (never set-on-import; the plan §2
ruling), RNG via `substream` → `jax.random.key`, `QuasisepGP` via tinygp's
`QuasisepSolver` or celerite2.jax (verify maintenance at track start, plan
§6), NUTS + VI via numpyro, the same icdf-fallback contract. Trace purity
throughout: no exception control flow on the engine path (the non-strict
ruling exists because of exactly this).
**Depends/Accept:** as W2.4, with the jax fixture.

### W2.6 — The lowering registry: `register_lowering` and the bijection slot [M]
The hardened hook ruled 2026-09-03 (`lowering.md` §12.8):
`register_lowering(family, backend, constructor)` keyed on the neutral
name, plus the per-backend custom-`Bijection` slot. Hardenings are part of
the contract: no overwrite of a built-in or existing row without
`override=True`; user-registered rows stamped in provenance; an **opt-in
conformance battery** a registrant runs against their own lowering
(reference-vs-native agreement); constructors must return trace-pure
objects (the registry resolves before tracing, so the mechanism itself is
inert to jit/vmap/grad — document that as a rule).
**Depends:** W2.1 (the reference rows to agree with); consumed by W2.4/W2.5.
**Accept:** a user-registered prior family lowers natively on a backend and
its provenance record says so; the overwrite refusal and `override=True`
path are tested; the registrant battery runs from the docs example.

### W2.7 — Diagnostics, the 1D slice [M]
Plan §4.8 / `diagnostics.md`: post-fit residual whiteness tests
(separation-binned, permutation-calibrated) and GP-localisation output
(`Likelihood.conditional` → `AnomalyScore` → the shared renderer) in
`ampere.results`; the `ampere.diagnostics` namespace with RHMF pre-fit
screening **if** the recorded adoptability assessment of Robusta-HMF
holds at implementation time (licence, maturity, API — re-verify), behind
the `diagnostics` extra. The namespace addition carries its decision-log
entry (`diagnostics.md` §7's note). Carries `diagnostics.md` §2.5's
obligation: validate `rank`/`robust_scale` heuristics on real ampere
collections before promoting any default, with its own decision-log entry
when that happens.
**Depends:** W2.2 (posteriors to diagnose).
**Accept:** both post-fit families produce output on the toy problem;
every `AnomalyScore` renders with provenance shown; RHMF path either
lands behind the extra or its deferral is recorded with reasons.

### W2.8 — Results completion: plots, pointwise emission, training sets [M]
Implement the declared plotting surface (corner, trace,
posterior-predictive, GP-localisation) against the stored `DataTree`;
the explicit-call emission of the `pointwise_log_likelihood` group
(`results.md` §6's reserved names over `Likelihood.pointwise_log_prob`,
never by default); and the training-set writer
(`serialisation_review.md` §4's obligations: `training_pair_from_dict`
completing the round trip, the θ-dtype loss fixed, batching/append for
large budgets).
**Depends:** W2.2.
**Accept:** each plot renders from a stored run with merged names as
labels; a run with the pointwise group survives netCDF and `arviz.loo`
consumes it; a training set written, appended to, and read back
round-trips by value.

### W2.9 — The WStat comparison example [S]
The 2026-09-03 ruling's obligation, carried here by W1.13: ampere ships no
WStat, and the docs take the opinionated line. A worked example builds the
profiled Cash-with-background statistic as a **user family** (safe under
masking via `NoiseParams.retain`; its `sample` refuses; per-sample
log-likelihood semantics degrade — say so), then the **two-dataset
Bayesian formulation** (source + background as a `DatasetCollection` with
a shared background model), runs both, and compares pros, cons and
*results* — posteriors, so this needs the engine drivers.
**Depends:** W2.1, W2.2.
**Accept:** the example runs end-to-end in the docs build; the comparison
states the recommendation and the trade-offs rather than false balance.

### W2.10 — Milestone M2: the flagship misspecification validation [L]
Reproduce the `flexible_likelihood_comparison` study on **both** modern
backends at 10–100× the current data size, with wall-clock benchmarks
against legacy (plan §5's M2). This is the paper-grade evidence the
redesign delivers its central promise, and the adversarial cross-model
review milestone (`docs/orchestration.md`).
**Depends:** W2.3, W2.4, W2.5.
**Accept:** posterior agreement between backends within stated tolerances;
the benchmark table produced by CI-runnable code, not by hand; adversarial
review dispositions recorded.
**Scoping (Fable, 2026-09-07, from reading `examples/examples_paper/flexible_likelihood_comparison.py` — 713 lines of legacy API, untracked in git; leave the script itself alone, it is paper-revision work).** The study is: a 4-parameter toy (linear continuum × two Gaussian absorption lines, uniform priors) on a Gaia-RVS-like grid (0.842–0.872 µm, **200 points**, 1% Gaussian noise, seed 42); four misspecification scenarios applied multiplicatively to the truth (none; 2.5% and 7% sinusoidal fringing at 0.0028 µm period; a 12% Gaussian emission line at 0.863 µm, σ = 0.00035 µm); two likelihoods per scenario (“standard” = GP priors pinned to 1e-10, “flexible” = RBF GP with scale-length prior 0.003 and amplitude prior 0.3); emcee 40 walkers × 5000 steps, burn-in 2000, DE + snooker moves; three figures (spectra + residuals with the GP predictive band, parameter recovery, 1-D posteriors). What v2 already has: the model (a native model per backend — the W2.1 trio does not include it, so W2.10 adds one small model to each of reference/torch/jax, or a single hand-written `Model` subclass per backend in the example), `Spectrum`, `GaussianFamily` + `GaussianProcessNoise` with Matérn-3/2 (the legacy study used RBF; **decide**: Matérn-3/2 is the v2 default and has the O(N) path, so the reproduction should use it and say so, with the RBF run kept only at 200 points as the legacy cross-check), `DenseGP`/`QuasisepGP`, the emcee driver, the DataTree emission. What is missing and must be built by W2.10 or before it: (a) the misspecified-data generators as a small reusable module (the four scenarios, parametrised by size), (b) the “standard likelihood” expressed honestly in v2 — a plain `IndependentNoise` likelihood rather than a GP with priors pinned to 1e-10, (c) the GP predictive band → `Likelihood.conditional` (W2.7's GP-localisation output; **W2.10 depends on W2.7** for the band and for the whiteness test that makes the comparison quantitative, or ships its own dense predictive as a stopgap), (d) the figures → W2.8's plotting surface (posterior-predictive and GP-localisation renderers; corner for recovery), (e) the size ladder 200 → 2 000 → 20 000 points, at which the dense solver is 1500× slower than the quasiseparable one (W2.3's scaling row), so the benchmark table wants both solvers at the small sizes and only `QuasisepGP` at the largest, (f) NUTS on torch/jax vs emcee on reference as the backend comparison, with agreement asserted on the 4 physical parameters' posterior medians and 68% intervals at stated tolerances, (g) the benchmark harness choice (pytest-benchmark vs asv, plan §6) — W2.11's, to be taken before W2.10 starts. **Dispatch order implied: W2.7 and W2.8 before W2.10; W2.11's benchmark half alongside.** Sizing: L holds; two agents (data + reference/legacy benchmark; torch/jax reproduction) is the natural split once W2.4/W2.5 slice 2 has `QuasisepGP` natively — without it the 20 000-point torch/jax runs are dense and the table is not the one the paper wants.

### W2.11 — CI/CD Phase 2 expansion [M]
Plan §5's cross-cutting workstream at this phase: separate torch, jax and
**no-extras** matrix jobs (the last catches lazy-import breakage —
`import ampere` must never require either), dependency caching, and the
benchmark suite (choose pytest-benchmark vs asv here, plan §6) with
results tracked as CI artefacts. GPU tests stay nightly/manual.
**Depends:** W2.4/W2.5 far enough along to have something to gate; the
benchmark half can trail with W2.10.
**Accept:** a deliberately broken backend row fails only its own job; the
no-extras job is green; benchmark results appear as artefacts on a PR.

### W2.12 — Backend identity on `FittingProblem` [M; Opus]
W2.2's carried finding: no §4 surface names a problem's backend, so
`Engine(problem, backend=...)` *declares* it and `ampere_backend` in a run's
provenance records a declaration rather than a fact — and two lockstep
tracks would each invent their own answer. **Decided by Fable 2026-09-07,
Peter to ratify at merge review** (a §4.5 addition — ground rule 9: this PR
carries the decision-log row): the backend is a **fourth capability flag**,
following the three existing ones exactly.
1. `Capable`, `Model` and `Transformation` gain `BACKEND: ClassVar[str] =
   "reference"` — the conservative default, since the base install's whole
   toolkit *is* the reference backend and a hand-written numpy model runs
   on the reference path. Every piece `ampere.backends.reference` ships
   declares it explicitly rather than inheriting. Every other capability
   part that carries the three flags (follow `Dataset.capability_parts`)
   carries the fourth.
2. `Capabilities.backend: str = "reference"`, validated non-empty like
   `device`; `declared_capabilities` aggregates by the **device rule**:
   all parts must agree, and disagreement raises `DatasetError` naming
   the backends and the override, `capabilities=Capabilities(backend=...)`.
   Ampere does not convert arrays between libraries on the user's behalf —
   a mixed problem is a configuration mistake that would otherwise fail
   two steps later inside a backend, or silently drop gradients.
3. `FittingProblem.backend` property beside `differentiable`/`batchable`/
   `device`.
4. **One name per backend, everywhere**: the string is the key
   `lowering.md` §12.8's registry uses for `backend`, and the conformance
   fixtures' `name`. Say so in the spec; a conformance row asserts each
   fixture's composed problem reports that fixture's `name`.
5. `Engine.__init__` **drops** `backend=` (pre-release; no shim) and reads
   `problem.capabilities.backend`. `ampere.results.emit` /
   `provenance_attrs`' `backend=` becomes optional, defaulting to the
   problem's own; an explicit value that disagrees with the problem
   raises rather than being recorded. `ampere_backend` is then a fact.
6. `Capabilities.to_dict` gains the key, so the `capabilities` payload in
   the provenance attrs changes shape: bump `PROVENANCE_SCHEMA_VERSION`
   to 4 per `results.md`'s rule. W2.6's deferred first-class lowering
   provenance key is **not** this bump — it stays deferred to W2.4/W2.5.
   Check whether `capabilities` feeds `ampere_problem_hash`; do not change
   the hash inputs, and state in the report what you found.
Spec amendments in the same PR: `inference.md` (the `Capabilities` row,
§18's Phase 2 bullet — a backend declares the *four* flags), `results.md`
(`ampere_backend` derived, schema version 4), `DEVELOPMENT_PLAN.md` §4.5's
capability-flags bullet, plus the decision-log row.
**Depends:** W2.1, W2.2, W2.6. **Blocks** W2.4/W2.5.
**Accept:** `pixi run test-all` green with new rows in `tests/core`
(agreement, disagreement raise naming the backends, the override, the
property), `tests/conformance` (item 4, per backend fixture — including
the in-repo third-backend demonstration), `tests/inference` (`Engine`
refuses `backend=`; the emitted attrs come from the problem) and
`tests/results` (the schema version; the disagreeing explicit value
raises); `grep -rn "backend" ampere/inference` shows no declared default;
the executed example in `ampere/inference/__init__.py` still prints
`ampere_engine`/`ampere_backend`; lint/format/pyrefly clean.

### W2.13 — The realisation surface, and the fold-ins [L; Opus]
`inference.md` §10a and its decision-log row (ruled 2026-09-07) turned into
code, from the prototype on `w2.13-realisation-prototype` (`ampere/core/
realisation.py`, the jax registration, `NUTSEngine(problem)` with no density
argument — keep all of it, extend it). In scope, in this order:
1. **Core**: `Realisation` gains the optional `log_likelihood_terms`; the
   Protocol and registry docstrings cite §10a; `ampere.core.exceptions.
   LoweringFallbackWarning` replaces `ampere.backends.torch.lowering.
   IcdfFallbackWarning` and `ampere.backends.jax.parameters.
   LoweringFallbackWarning` (both backends import the shared one);
   `ParameterSet.evaluation_order()` and `Dataset.effective_mask` public,
   both backends switched to them; the four capability ClassVars on
   `NoiseModel` and `GPSolver` (defaults `False, False, "cpu", "reference"`),
   `Likelihood.capability_parts` → its noise model and, when a GP is declared,
   its solver, `Dataset.capability_parts` appending them; `GPSolver.
   provenance_config()` (default `{}`) recorded by `provenance_attrs` under
   `ampere_solver_config`, never hashed.
2. **torch**: `ampere.backends.torch.problem.LoweredProblem`/`lower_problem`
   built from the tensor twins W2.4 shipped, at jax's coverage floor
   (Gaussian family, independent/dense-GP noise, masks, plates, hierarchical
   priors), refusing the rest by name at construction; native kernels
   (`Matern32`, `SquaredExponential` subclassing core's) so the torch
   `DenseGP` differentiates in the hyperparameters; `BACKEND`/flags on the
   torch `DenseGP`; registration at import. A pyro-backed `NUTSEngine`
   route: extend `ampere.inference._nuts` to dispatch on `problem.backend`
   between numpyro (jax) and pyro (torch) — both lazily imported inside
   `run`, both via a potential function; `SUPPORTED_BACKENDS` becomes
   "whatever is registered".
3. **jax**: adopt the shared warning, the public accessors, the flags on
   `DenseGP`/kernels/noise; `log_likelihood_terms` on `LoweredProblem`.
4. **Provenance**: `PROVENANCE_SCHEMA_VERSION` → 5; `ampere_realised` and
   the consulted registered lowering rows (`provenance_entries`) stamped by
   the drivers that sample through a realisation.
5. **Conformance**: a row per registered realisation comparing
   `log_prob_unconstrained` with the numpy path at many points including
   near a boundary (`tolerances.cross_backend`), and the identity rows
   extended to the widened parts set; `tests/core` rows for every new core
   surface.
6. **Specs**: the three corrections in the decision-log row (`lowering.md`
   §3.6, §4, and the `AffineTransform` domain note), each marked *Amended
   W2.13*.
**Depends:** W2.4 slice 1, W2.5 slice 1, the prototype. **Blocks** both
slice 2s.
**Accept:** `pixi run -e torch test-all`, `-e jax test-all` and `-e dev
test-all` green (the dev run proving neither backend is imported); NUTS
recovers the toy joint posterior on **both** backends with
`NUTSEngine(problem)` and no density argument; the per-realisation
conformance row runs for both; a numpy `DenseGP` in a native problem is
refused as a backend disagreement; the two fallback-warning classes are
gone; schema 5 with `ampere_realised` in an emitted run; lint/format clean,
pyrefly 0 errors in all three environments.

## Status

| Item | Status |
|---|---|
| W0.1 | done 2026-09-01 (committed directly to master — the pending fix existed only in the local working tree, so the branch-per-item rule was waived for it) |
| W0.2 | merged 2026-09-01 |
| W0.3 | merged 2026-09-01 |
| W0.4 | merged 2026-09-01 |
| W0.7 | merged 2026-09-01 — branch-archival actions in `docs/design/harvest/branch_triage.md` await Peter's approval |
| W0.5 | merged 2026-09-01 |
| W0.6 + W0.8 | merged 2026-09-01 |
| W0.9 | merged 2026-09-03 (Sonnet-authored, Fable-reviewed; no fixes needed — compat layer live-probed under pyphot 2.1.1, no residual `pyphot.unit` in the package, gates re-run independently: dev characterisation 3 passed no-xfail, sbi run 6 real passes on sbi 0.27.0, phase1 1109, lint/format/pyrefly clean). `ampere/utils/pyphot_compat.py` presents the ≥2 idiom as its primary surface — **Phase 2's synthetic-photometry step is now ungated**. Caveat recorded in the plan's pins row: no sbi ≥0.28 exists on PyPI yet, so re-check when one ships. Out-of-scope finding: `examples/minimal_working_example*.py` and the paper scripts still call `pyphot.unit[...]` directly and break standalone under pyphot ≥2 — needs a small follow-up (examples touch-up) in Phase 2 |
| W0.10 | merged 2026-09-03 (Sonnet-authored, Fable-reviewed; no fixes needed — gates re-run independently: single-process phase1 1060 green, py3.11 leg spot-checked, lint/format/pyrefly clean, ci.yml YAML-validated, lockfile diff confined to default/dev). One correction to the item's own text: finding (c)'s leak source is the `likelihoods.md` "laplace" doctest via `test_spec_doctests.py`, not the `test_likelihood.py` sites the item named — the autouse fixture covers both. Live-CI verification deferred until Peter pushes |
| W1.1 | merged 2026-09-01 (Sonnet-authored, Fable-reviewed; two load-bearing claims re-verified against live sources). Peter's review pass remains the formal accept gate |
| W1.2 | merged 2026-09-01 (Sonnet-authored, Fable-reviewed and amended; reconciled against W1.1, see its §10). The curated-translation conflict found in review was ruled **opt-in, never silent** — decision recorded in `DEVELOPMENT_PLAN.md` §2, §4.7 amended. Peter's review pass remains the formal accept gate |
| W1.3 | merged 2026-09-01 (Opus-authored, Fable-reviewed; three tying/serialisation defects fixed on the branch with tests — 112 pass, pyrefly/ruff clean). Remaining open questions for Peter in the spec's §14; obligations on W1.4–W1.10 in its §13 |
| W1.4 | merged 2026-09-01 (Opus-authored, Fable-reviewed; two defects fixed on the branch — `from_unsorted` extra-coordinate misalignment, `to_unit` error type — 247 core tests, pyrefly/ruff clean). Open questions for Peter in the spec's §17 |
| W1.9 | merged 2026-09-01 — **both flagged rulings approved by Peter** (`eqx.partition` over paramax; x64 guard-and-raise) and recorded in the plan's §2 decision table with §4.1/architecture.md §5 amended. Remaining §12 ratification items route to W1.13, now including item 8 (backend-specific lowerings for user-defined priors/bijections, from Peter's review) |
| W1.12 | merged 2026-09-01 — **approved by Peter** (§10 defaults stand). One addition from his review, recorded as the spec's §11: posterior calibration (SBC/coverage) as a fourth, future family landing with Phase 3's SBI layer, plus the extension template for families to come |
| W1.5 | merged 2026-09-01 (Opus-authored, Fable-reviewed; no defects needed fixing — the agent's self-review caught the per-interval-density union bug itself). 376 core tests. Open questions for Peter in the spec's §15; `Model` ABC placed here (§9), mask rule ratified (§6) |
| W1.6 | merged 2026-09-01 (Opus-authored, Fable-reviewed; no fixes needed — the critical silent-GP defect was caught pre-review by the adversarial pass the agent commissioned, and fixed with the `CONSUMES_LATENT_GP` opt-in). DenseGP anchor validated to ~1e-15; 491 core tests post-merge. Open questions for Peter in the spec's §17 |
| W1.11 (a, b, d, e) | merged 2026-09-02 (Sonnet-authored, Fable-reviewed; every gap claim reproduced by execution during review; one fix on the branch — the proposed `Axis.locate` must match within `COORDINATE_RTOL`, not exactly). Five gaps recorded as proposed amendments, all routed to W1.7/W1.13 |
| W1.11 (c, f + stress test) | merged 2026-09-02 (Opus-authored, Fable-reviewed; no fixes needed — structural claims verified against the code, decisive probes re-run, von Mises table checked analytically). All four routed claim-checks answered; fourteen gaps recorded as proposed amendments, five flagged for the freeze (I-1, I-5, X-1, X-2, H-2). Item complete pending W1.13's resolution of the amendments |
| W1.7 | merged 2026-09-02 (Opus-authored, Fable-reviewed; gates verified independently — 671 tests, 174 spec doctests, pyrefly/ruff clean — and adversarial probes beyond the agent's own two passes all held; no fixes needed). Implements plan §4.5 verbatim plus Peter's two same-day rulings; the agent's commissioned review caught and fixed the mask-escape free maximum, the deferred-validation escape, and the dropped solver jitter. **Ruling requests R1–R4 in `inference.md` §19 — R1 (merge topology) must be ruled before W1.8/W1.10 dispatch** |
| W1.8 | merged 2026-09-02 (Opus-authored, Fable-reviewed; the agent's commissioned adversarial review found five provenance/validation defects — all reproduced and fixed pre-review — and the Fable pass found and fixed one more: a string value equal to a float sentinel hashed identically to the real NaN). 834 tests, 87 spec doctests; cross-process hash determinism verified. **Ruling requests R1–R7 in `results.md` §15**; arviz stays a lazy extra with the pixi dev env gaining it plus a netCDF engine |
| W1.10 | merged 2026-09-02 (Opus-authored, Fable-reviewed; kernel oracles verified bit-exactly against hand computations; no fixes needed on the branch). 275 backend-parametrised rows over two in-repo fixtures, a third backend demonstrated as one registry line; QuasisepGP agreement is a skipped skeleton naming Phase 2's debt, with a live row asserting the empty slot refuses. Its review found and led to fixing (on master, `a8e6f54`): the py3.11 core-import break (seven mappingproxy dataclass defaults) and `results.md` §14's overreaching problem-hash claim |
| W2.1 | merged 2026-09-05 (Opus-authored; **Fable-reviewed twice** — the second pass ruled by Peter in lieu of the cross-model review while the Codex quota is exhausted; the retroactive `gpt-5.6-terra` pass on the merged range is owed when quota returns, ~2026-09-30). The reference backend is a real package: three native models, the four standard steps, `FractionalModelNoise`(+GP), and the conformance battery rewired so its `reference` fixture is an adapter onto the shipped package (Peter-confirmed — no test-local duplicate kept). Both post-freeze §4 landings carry their decision-log rows: `describe()` + the neutral model identity (`PROVENANCE_SCHEMA_VERSION` → 3) and `Axis.locate` (`COORDINATE_RTOL` moved to `results_schema.py`, re-exported). Peter's 2026-09-05 photometry ruling implemented on the branch: detector-aware weighting (photon `R dλ/λ`, energy `R dλ/λ²`, per filter), `detector=` required, `from_library` reading pyphot's `Filter.dtype`; the draft `R dλ` weighting was neither convention and was replaced before release. Gates at merge: phase1 1133 in one process, backends 85, lint/format/pyrefly clean. Out-of-scope findings, **both ruled by Peter and landed on master 2026-09-07**: pyphot 2.1.1's undeclared `requests` dependency broke the minimal-install CI job (W0.9 aftermath) — `requests` declared on pyphot's behalf (`a02836d`); `Instrument.__init__` ran `configure_from` in forward chain order, under-padding chained same-axis convolutions (core, pre-existing) — now last step first, with regression tests at both levels and a decision-log row (`8a9534a`). `celerite2` deliberately not added to base deps — W2.3 first uses it |
| W2.6 | merged 2026-09-05 (Sonnet-authored, Fable-reviewed; two findings fixed on the branch and re-verified — the lowering-registry snapshot fixture in `tests/core/conftest.py` (W0.10 finding (c)'s leak class on a new registry) and module-qualified bijection keys (a silent wrong-lookup path for same-named classes closed; bare names kept in messages). `ampere/core/lowering.py` implements `lowering.md` §12.8 as ruled: neutral-name keying, the bijection slot, overwrite hardening, `provenance_entries` via `provenance_attrs`'s `extra=` (first-class schema key deferred until a real backend drives lowering — recorded), the opt-in registrant battery run from the module docstring's executed example. Consumed by W2.4/W2.5 |
| W2.2 | merged 2026-09-05 (Opus-authored, Fable-reviewed; no fixes needed). `ampere.inference`: emcee/dynesty/zeus against §4.5 only — backend neutrality import-graph-tested; every run emits the DataTree with per-dataset log-likelihood decomposition, provenance attrs and a hash-preserving netCDF round trip; three engines agree on an analytically known posterior (oracle guarded by its own row). Carries R1's promotion with its decision-log row: arviz + h5netcdf into base deps, the `arviz` extra REMOVED (loud failure over an empty alias), pixi feature renamed `netcdf`, `tests/results` + `tests/inference` join the CI matrix. Findings carried forward: zeus draws from TWO process-global RNGs (numpy legacy + stdlib `random`; both seeded/restored by the driver — any future zeus/Pool work inherits this); zeus wants a burn-in/re-ball helper (small future item); the reference hot loop pays ~28% of ~4 ms/eval rebuilding containers when no requirement was published (`compile_for` template gap — cheapest visible perf win); no `FittingProblem` surface names its backend, so `Engine(..., backend=)` is declared — wants a small addition before W2.4/W2.5 |
| W2.3 | merged 2026-09-05 (Opus-authored, Fable-reviewed with the GP algebra hand-checked line by line in lieu of the blocked cross-model pass — first in the retroactive terra queue). `QuasisepGP` is real in `ampere.core` (placement recorded in the decision log; `architecture.md`'s module-address wording clarified in place): ampere's own **exact** rank-2 semiseparable Matérn-3/2 term (algebraic identity, verified against the dense kernel matrix to 1.7e-16) rather than celerite2's ε-limit `Matern32Term` (misses `tolerances.cross_solver` at its default); coordinates midpoint-centred against cancellation; internal sort under permutation invariance. W1.10's skipped solver rows now run and pass (agreement 2.7e-15); the refusal row became the stronger not-the-dense-solver-in-disguise row (declared-slot discipline kept in `tests/core`). `conditional_loo` DEFERRED with refusal live and the O(N) route recorded. celerite2 into base deps, imported lazily. Scaling demonstrated: ~57 ms at 10⁵ points, exponent ~0.94, >1500× over dense. Finding for W2.4/W2.5: celerite2 returns quiet NaN where a Cholesky raises — guard preconditions on every backend's celerite2 path |
| W2.12 | merged 2026-09-07 (Opus-authored, Fable-reviewed; no fixes needed — gates re-run independently: test-all 1368, zero skips, lint/format/pyrefly clean). The backend is §4.5's fourth capability flag: `BACKEND` on `Capable`/`Model`/`Transformation` (default `"reference"`, declared explicitly on every shipped reference piece), `Capabilities.backend` aggregated by the device rule with the `capabilities=` override, `FittingProblem.backend`; `Engine` drops `backend=` and reads the problem; `emit`/`provenance_attrs` keep it as an optional cross-check that raises on disagreement; `PROVENANCE_SCHEMA_VERSION` 4 (`capabilities` is not a fingerprint input, the schema constant is). One name per backend everywhere (flag = registry key = fixture `name` = `ampere_backend`), asserted per fixture by the new conformance identity rows — the mirror fixture gained three one-line step subclasses so the row is not a tautology. Decision-log row present; **Peter ratifies the contract decision at his review of the merge**. Carried findings: `Dataset.capability_parts` excludes the likelihood/family/noise model, so a native GP solver contributes nothing to any aggregated flag (pre-existing; widening is a §4 change wanting a ruling); `ZeusEngine`'s `**sampler_settings` remains a silent-typo surface for other misspelt keywords |
| W2.4 (slice 1) | merged 2026-09-07 (Opus-authored, Fable-reviewed; no fixes needed — gates re-run independently on the branch: torch test-all 1745/7 skipped, dev 1368 unchanged, lint/format clean, pyrefly 0 errors in both environments). Items 1–9 of slice 1 landed: `ampere.backends.torch` — every `lowering.md` §3.2/§4 row registered through the core registry (`builtin=True`), `truncnorm` refused by name, exact affine composition for shifted families (`beta` included: §3.3 rule 1 outranks the table's loc/scale note), **declared true supports** because torch's `TransformedDistribution.support` is the last transform's codomain (a spec clarification owed for §4's `biject_to(lowered.support)` advice), nested `nn.Module` per merge component, icdf fallback with once-per-run `IcdfFallbackWarning` and strict `LoweringError`, `substream` → `Generator.manual_seed`, float64 threaded with `set_default_dtype` never called (tested), the model trio and four steps with tensor twins, a torch `DenseGP` agreeing with the oracle to 1e-9, the fixture green on every cross-backend row. **Item 10 (NUTS) deferred** on the verified structural gap (containers coerce to numpy; `ampere.inference` may not import a backend) — the realisation ruling, proposed as W2.13. Findings carried: the torch `DenseGP` covariance comes from core's numpy `Kernel.matrix`, so GP hyperparameters get no gradient (W2.5 solved this with native kernels — torch slice 2 should mirror it); `Likelihood.to_spec` records dataclass solver fields, so backend-specific solver config cannot be dataclass fields without breaking cross-backend spec-hash agreement (ruling before any float32/GPU work); `pip install ampere[torch]` resolves to a ~2 GB CUDA wheel on Linux; no torch noise models yet; `ParameterSet` wants a public `evaluation_order()`; `run_registrant_battery`'s zero-argument form is unusable by a tensor backend. Slice 2: `QuasisepGP` (GPyTorch vs celerite2-torch, measured), VI, batching/GPU, the realisation consumer |
| W2.5 (slice 1) | merged 2026-09-07 (Opus-authored, Fable-reviewed; no fixes needed — gates re-run independently on the branch: jax conformance 475/5 skipped, backend+NUTS+import suites, lint/format/pyrefly clean; the agent's own jax test-all 1662/5 and dev 1369/3). **All ten slice-1 items landed**: `ampere.backends.jax` with `configure_x64()`/`require_x64` (guard-and-raise naming the three remedies; never at import), every §3.2/§4 row through the core registry with `biject_to` and explicit `domain=` on composed affine maps (numpyro's `AffineTransform` defaults its domain to `real` — a §4 note owed), real `numpyro.plate`s and explicit `to_event(n)`, `eqx.partition` filter specs (never static fields), the icdf contract probed behaviourally (**`lowering.md` §3.6 is wrong: numpyro's `Gamma`/`Beta` icdf delegate to TensorFlow Probability and raise without it** — spec correction owed), keys threaded never stored, **native kernels** (`Matern32`, `SquaredExponential` subclassing core's) so GP hyperparameters are trainable, a jax `DenseGP`, the steps subclassing the reference ones so declarations are written once, `lower_problem` → `LoweredProblem` (trace-pure: refusals at lowering time, −inf via `jnp.where`), the fixture green on every row, and **`NUTSEngine` in `ampere.inference`** recovering the toy joint posterior and the conjugate closed form to 0.25σ, taking the lowered density as an argument so the import-graph rule stays literally true. Findings carried: `celerite2.jax` flips `jax_enable_x64` on import (the very thing §10.2 forbids — slice 2's library choice must weigh it); `Likelihood.describe()` records solver `config` only for dataclasses (cross-backend-visible); `Dataset._effective_mask` and `ParameterSet._order` reached privately; the lowering provenance key stays deferred (every consulted row is built-in, so a first-class key would be empty — `lowering_provenance()` exists and is tested with a user-registered family). Slice 2: the GP-library choice (survey in `gp.py`'s docstring: celerite2.jax at zero dependency cost vs tinygp 0.3.1), `QuasisepGP` + the O(N) `conditional_loo`, widening the native path to censoring/latent GPs/non-Gaussian families, VI, `vmap`/`device_put`, benchmarks |
| W2.13 | merged 2026-09-07 (Fable-prototyped and contract-drafted, Opus-implemented, Fable-reviewed; no fixes needed — gates re-run independently in all three environments). `inference.md` §10a is code: `ampere.core.realisation` (`register_realisation`/`realise`/`registered_realisations`, `log_likelihood_terms_of`), the torch `LoweredProblem` mirroring jax's line for line (refusals at construction, −inf via `torch.where` because pyro refuses a potential without a `grad_fn`), native torch kernels (the amplitude now has a gradient — W2.4's principal finding closed), both backends' `IndependentNoise`/`GaussianProcessNoise` declarations under the core names, `NUTSEngine.run` dispatching on `problem.backend` between numpyro and pyro (both lazy, both over a potential function, chains sequential, pyro seeded inside `fork_rng`), `supported_backends()` as the registry ∩ `SAMPLER_LIBRARIES`. Fold-ins landed: four flags on `NoiseModel`/`GPSolver` with `Likelihood.capability_parts`; `ampere.core.exceptions.LoweringFallbackWarning`; public `ParameterSet.evaluation_order()`/`Dataset.effective_mask`; `GPSolver.provenance_config()` → `ampere_solver_config` (recorded, never hashed — tested). Schema **5**: `ampere_realised`, `ampere_registered_lowerings`, `ampere_solver_config` on every run. Conformance: `TestTheRealisation` (24 points × 3 shapes per registered realisation, decomposition, refusals; reference/mirror skip by the ruling), `ConformanceBackend` gained `independent_noise()`/`gp_noise()`. The three spec corrections marked *Amended W2.13*. **Consequence for users (Peter to note)**: `Dataset(observed)` on a native backend is now a `DatasetError` naming the classes — the default noise model is core's and declares `reference`; every native problem passes `Likelihood(GaussianFamily(), backend.IndependentNoise())`. Ruled as written; softening options recorded in the handoff. Carried: jax `SyntheticPhotometry.apply_flux` is broken (`influence` needs an `Axis`; torch's rebuilds it) and no conformance spec has a photometry step — slice 2 (jax) fixes it and adds the spec; torch ships no prediction-aware noise models (`sigma_tensor` hook dormant) — slice 2 (torch); `log_likelihood_terms` is supplied by both realisations but not consumed (`Engine.finish` must learn to accept a decomposition it did not compute) — either slice 2 |
| W1.13 | **merged 2026-09-03 — Phase 1 complete, spec frozen** (Fable-authored across a session-limit resume, Fable-reviewed independently: every contract change read against its recorded ruling, LOO algebra hand-checked, gates re-run — 1109 tests in one process, py3.11 leg, lint/format/pyrefly — and no fixes needed). 18 commits: the twelve ruled contract changes, cross-review harmonisation, the W1.11 gap dispositions, the consolidated serialisation review (`docs/design/serialisation_review.md`), the freeze record, the Phase 2 breakdown (W2.1–W2.11), the `gpt-5.6-terra` adversarial pass (two real pre-existing defects fixed — the censoring/masking union in `draw_observation`, loud seed validation — one draw-epsilon nuance documented for Phase 2), and Peter's rulings on the three escalations (model identity + `describe()` → W2.1; `Population` → Phase 5, adaptable; `Axis.locate` approved → W2.1). Tag `spec-v1.0` at the freeze-content commit `58daa86` |
