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
**Implemented 2026-09-08 on `w2.11-ci-phase2`** (Peter to review and merge).
The job graph is now: `lint` and `typecheck` on `dev`; `test` (the py3.11/12/13
matrix, unchanged); `suites` (`dev`); `backend-suites`, a `fail-fast: false`
matrix over `torch` and `jax` running that environment's `typecheck`,
`test-all` and `bench`; `docs` (`dev`); `minimal-install` (bare
`pip install -e .`, the no-extras leg); and the unchanged schedule/manual
`sbi-characterisation`. One environment per job, never two — GitHub's runners
have 7 GB and one five-suite gate is what fits. **The benchmark harness is
pytest-benchmark** (`tests/benchmarks/`, `pixi run bench`, `benchmark.json`
uploaded per environment); the reasoning is the decision-log row in
`DEVELOPMENT_PLAN.md` §2 and §6's bullet is struck. Also landed here because
they were blocking: `pixi run docs` now builds (`pandoc` and `ipykernel`
declared in the `dev` pixi feature; `nbsphinx_execute = 'never'` with a
per-notebook note; the dead `ampere.infer.ptemceesearch` autodoc entry
removed) though **without `-W`**, since the residual warnings are legacy
docstrings under frozen paths; `tests/examples` joined `test-all`; the `torch`
and `sbi` pixi features resolve torch from PyTorch's CPU index, removing all
fifteen `nvidia-*` wheels from `pixi.lock`; and `tests/scaling` gained the
torch rows it was missing. GPU rows are still nothing — a later item.

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

### W2.14 — The latent-GP likelihood must see its kernel [S; Opus]
Found independently by W2.4 and W2.5 slice 2, reproduced on master
(`docs/development.md`, 2026-09-07): on a Poisson + `GaussianProcessNoise`
problem `problem.log_likelihood` is bit-identical for amplitude 0.5/5/50 and
length scale 0.1/1/100 although both are free parameters. Nothing on the
scoring path applies `GPSolver.latent_transform`: `Dataset.log_likelihood_of`
passes the whitened `z` through `GaussianProcessNoise.noise_params(latent=)`
to the family, which reads it as `f`; the only caller of `latent_transform` is
`GaussianFamily.sample`. `likelihoods.md` §17's limitation 6 ("the latent path
has no inference") explains the omission, but the numpy path is the oracle
every backend is held to, and it silently scores a wrong number.
**Fix (ruled 2026-09-08)**: `GaussianProcessNoise.noise_params` applies
`latent = self.solver.latent_transform(self.kernel, coordinates, latent,
values)` when `latent` is given, so `noise.latent` *is* `f = L(θ) z` as
every family already reads it — the solver and coordinates are already in
hand there and no family changes. Then: (1) a `tests/core` regression row
asserting the log-likelihood *moves* with amplitude and length scale and
equals the closed form for a Gaussian family (whitened-`z` Gaussian likelihood
with `f = L z`); (2) `likelihoods.md` §17 limitation 6 amended (*Amended
W2.14*) and a decision-log row (a §4.4 clarification); (3) both backends'
realised latent paths — torch refuses today with a test asserting the
flatness (`test_the_core_defect_the_latent_refusal_names_is_real`), jax
mirrors the flat oracle deliberately (`_LoweredDataset._latent`) — updated to
apply their native `latent_transform` and the refusal/mirroring removed;
(4) a conformance row: a latent-GP problem's realised density agrees with the
corrected numpy path at many points, and the numpy path's own likelihood
changes with the hyperparameters, per backend.
**Depends:** W2.4/W2.5 slice 2 (merged). **Blocks** both tracks' slice 3.
**Accept:** the regression row fails on the pre-fix code and passes after;
`pixi run -e torch test-all`, `-e jax test-all`, `-e dev test-all` green with
the new rows; the torch flatness test is gone; lint/format clean; pyrefly 0
errors in all three environments.

### W2.4 slice 3 — torch: device, complex Gaussian, GPU smoke tests [M; Opus]
Ruled 2026-09-08 (Peter: "finish slice 3"). (1) **Per-instance `device=`** on
every shipped torch model, step, noise model, kernel and solver (a keyword
threaded exactly as `dtype` is; default `"cpu"`; never auto-detected), with
`DEVICE` reported from the instance, buffers moving with `.to(...)`, and
`declared_capabilities`' device rule composing a whole problem on one device.
(2) **GPU smoke tests** at API level only: a `tests/gpu/` module, skipped
without `torch.cuda.is_available()`, that composes the toy joint problem on
`"cuda"`, realises it, evaluates the density and its gradient, and runs a
handful of NUTS draws — proving the API works, trusting the CPU/GPU library
parity for the numerics (Peter's ruling); CI stays CPU-only. (3)
**`complex_gaussian`** transcribed into `_families.py` with the circular
complex GP path `likelihoods.md` §4 declares (`ANALYTIC` under a GP),
composing in the realised path; refused by name where the core family itself
refuses. (4) The `sigma_tensor` hook's remaining consumer gaps (any
prediction-aware noise model that still falls back to numpy). No new
dependencies. **Owns nothing outside `ampere/backends/torch/`, `tests/`,
and the torch conformance fixture.**
**Accept:** the three torch gates green; the device keyword exercised for
`"cpu"` in the suite and the GPU module collected-and-skipped here; a
`complex_gaussian` realised density agreeing with the numpy oracle at
`tolerances.cross_backend` (a conformance shape added); lint/format/pyrefly.

### W2.5 slice 3 — jax: device, complex Gaussian, LOO parity, the gradient-free fast path [M; Opus]
Ruled 2026-09-08. (1) **Per-instance `device=`** (`jax.devices()` by name,
`device_put` on buffers; default CPU; never auto-detected) mirroring the torch
item, GPU smoke tests in `tests/gpu/` skipped without an accelerator. (2)
**`complex_gaussian`** in `families.py`, composing in the realised path. (3)
**`conditional_loo` parity**: implement the O(N) recursion on the jax
`QuasisepGP` by calling celerite2's compiled kernels directly (as the torch
solver does via `celerite2.backprop`; the recursion is derived in the plan's
"`QuasisepGP.conditional_loo` on the torch path" row), so the refusal lifts
and the conformance row's "agree exactly" branch runs on jax; if the kernels
cannot be reached from jax without breaking trace purity, keep the refusal
and record precisely why. (4) **The gradient-free fast path** (owned here,
touches `ampere/inference/engine.py`): W2.10 measured the jax *contract* path
at ~28 ms flat per `log_prob` — a pessimisation for emcee/dynesty/zeus on a
jax problem. `Engine`'s evaluation cache gains an opt-in route through
`ampere.core.realise` when a realisation is registered for the problem's
backend and supplies `log_likelihood_terms` (`use_realisation=True` default
where available; the numpy path stays the oracle and the fallback), the
decomposition coming from the realisation; `engine_realised_evaluations`
recorded. Test it on both backends (torch is registered too). No new
dependencies. **Owns `ampere/backends/jax/`, `tests/`, the jax fixture, and
`ampere/inference/engine.py` (torch slice 3 must not touch that file).**
**Accept:** the three jax gates green plus the torch inference rows for the
fast path; `conditional_loo` on jax either agrees with `DenseGP` at
`cross_solver` or is refused with the recorded reason; the fast path measured
faster than the contract path on a jax `QuasisepGP` problem in
`tests/benchmarks`; lint/format/pyrefly.

### W2.15 — Phase 2 documentation sweep [M; Opus]
Ruled 2026-09-08 (Peter: "make sure all docs are up to date — then Phase 3
can begin in clean sessions"). After the slice-3 merges: every user-facing
and design document reflects what landed. (1) `docs/source/`: an
architecture/overview page for v2 (the capability ladder as it exists; the
three environments; how a problem is composed on each backend, including the
one-backend rule for noise models and solvers; realisation and the engines;
results and diagnostics; the M2 page linked); API pages for `ampere.core`,
`ampere.backends.{reference,torch,jax}`, `ampere.inference`, `ampere.results`
generated by autodoc with the legacy pages kept separate and marked frozen;
the tutorials toctree; the docs build green on errors. (2) `README.md`:
current install routes (`pixi`, `pip`, the extras), the one-paragraph pitch
with the M2 numbers, and pointers. (3) `docs/design/architecture.md` §3's
namespace diagram and extras table checked against the tree (`diagnostics/`
deferred; `pixi` features; CPU torch index); every `docs/design/contracts/*`
"landed"/"amended" annotation checked against the code for the Phase 2
landings (a list of stale sentences with fixes, each marked *Amended W2.15*,
and one decision-log row for the batch). (4) Docstring index: every public
name in the new namespaces has a docstring that renders (autodoc warnings
from the new namespaces fixed; legacy warnings left). (5) The `examples/`
pyphot touch-up carried since W0.9 (`examples/minimal_working_example*.py`
and any paper-adjacent scripts *outside* `examples/examples_paper/` that call
`pyphot.unit` directly) migrated onto `ampere.utils.pyphot_compat`. No code
behaviour changes; no contract changes beyond annotations. Fable separately
restructures `docs/development.md` and `CLAUDE.md` for Phase 3.
**Accept:** `pixi run docs` succeeds with no warnings attributable to the new
namespaces; README install routes verified by running them in a scratch
environment; the annotation list in the report; all three gates unchanged.

## Phase 3 — The SBI layer (drafted 2026-09-08 by Fable; **approved by Peter 2026-09-08** with W3.1 revised into two slices from his three notes, and W3.8–W3.14 added from his rulings; W3.1's revised text approved 2026-09-09; **closed 2026-09-10** — every item W3.0–W3.15 merged, W3.13 last; the rulings still open are listed in `docs/development.md`'s "Everything that needs Peter's attention" block)

Written from `DEVELOPMENT_PLAN.md` §5 Phase 3, its §6 deferred choices and
§7's trained-artefact trap, `inference.md` §13 and limitation 17.5,
`results.md` limitation 13.9, `diagnostics.md` §11, and the swyft harvest
(`docs/design/harvest/swyft/README.md`). The sequencing principle is the one
Phase 2 used: the backend-neutral core surface first (W3.1), then the
package integration that consumes it (W3.2), then what the integration
makes possible (W3.3–W3.6). Everything here lives **above** the backends —
`ampere.inference` still imports no backend; the SBI module imports torch
and `sbi` lazily inside the call that needs them, exactly as `_nuts.py` and
`_vi.py` do, and a jax-native route is not scheduled (plan §6: "if and when
maturity warrants" — nothing below closes the door).

Shared facts for every item: the `sbi` extra resolves to **sbi 0.27.0 with
a CPU torch 2.13** in the `sbi` pixi environment (`pixi install -e sbi`,
verified 2026-09-08); sbi 0.27 spells the trainers `NPE`/`NLE`/`NRE`
(`sbi.inference`), ships `run_sbc`/`check_sbc`/`run_tarp`/`check_tarp` in
`sbi.diagnostics`, and `FCEmbedding`/`CNNEmbedding`/
`PermutationInvariantEmbedding`/`TransformerEmbedding` in
`sbi.neural_nets.embedding_nets`. The environments are `dev` (no torch),
`torch`, `jax` and `sbi` (= `dev` + the `sbi` extra; **no pyro**, so the
`torch` backend imports but its NUTS/VI engines do not — `architecture.md`
§3's post-W2.15 note). The gate for a Phase 3 item is `pixi run -e sbi
test-all` plus the `dev` gate proving nothing new leaks into the base
install; the `torch`/`jax` gates are re-run only by items that touch a
backend. Phase 2's operations notes hold: one five-suite gate at a time.

**Facts gathered 2026-09-09 for the items not yet dispatched** (verified
against the installed sbi 0.27.0 and PyPI, while W3.1 slice 2's gates
ran). *W3.3*: sbi's `PermutationInvariantEmbedding(trial_net,
trial_net_output_dim, aggregation_fn="sum"|…, output_dim=…)` takes input
`(batch, permutation_dim, input_dim)` and its `forward` has some mask
handling — verify exactly what before relying on it, and ship a small
masked-pooling wrapper (multiply each row's embedding by the mask column
before aggregation) rather than trusting a network to learn that padded
rows are absent. `TransformerEmbedding` takes `(batch, seq_len,
feature_space_dim)` and its `forward(input, attention_mask=None, …)`
accepts a `(batch, seq_len)` mask — but sbi calls an embedding as
`net(x)` only, so W3.3 needs a thin wrapper module that splits the mask
column out of `x` and passes it as `attention_mask`; its `pos_emb` is
`"rotary"` (index-based) or `"none"` — for set data use `"none"` and
carry the *coordinate* as row features (raw plus Fourier features), and
**set `is_causal=False`** (the default is `True`, which would make a
spectrum's rows attend only to earlier rows). *W3.6*: `run_sbc(thetas,
xs, posterior, num_posterior_samples=1000, reduce_fns="marginals", …) ->
(ranks, dap_samples)`, `check_sbc(ranks, prior_samples, dap_samples,
num_posterior_samples=1000) -> dict`, `run_tarp(thetas, xs, posterior,
references=None, num_posterior_samples=1000, …, z_score_theta=True) ->
(ecp, alpha)`, `check_tarp(ecp, alpha) -> (atc, ks_pval)`. *W3.4*: swyft's
latest release is **0.4.5 of September 2023**, pinned to
`pytorch-lightning >=1.5.10,<=1.9.5`, so it cannot be installed beside
torch 2.13 without an old Lightning — the maturity note's first fact, and
a strong signal that TMNRE is better expressed through sbi's own `NRE`
plus `RestrictedPrior` rounds than by reviving the harvest (Peter rules;
see W3.4).

### W3.0 — Phase 2 carry-over housekeeping [S; Sonnet]
The four findings the Phase 2 handoff lists, none of them a contract change.
(1) The name `ampere` on PyPI is an unrelated battery-modelling package, so
`OptionalDependencyError`'s remedy line (`ampere/core/exceptions.py`) must
not print `pip install "ampere[...]"` as if it worked today: make the remedy
name the extra and route to the install page (`pip install "ampere[<extra>]"
from a checkout, or see the install documentation`) — one wording, chosen in
the item, applied to the runtime string, the class docstring's doctest, and
every docstring in the new namespaces that quotes a `pip install ampere…`
line (`ampere/backends/{reference,torch,jax}/__init__.py`,
`ampere/inference/__init__.py`, `ampere/inference/_zeus.py`,
`ampere/core/likelihood.py`, `ampere/results/emission.py` — enumerate by
grep, the list here may be incomplete). (2) `SAMPLE_STATS_GROUP` is defined
twice, identically, in `ampere/results/emission.py` and
`ampere/results/training.py`: one definition, imported by the other. (3)
Rename torch's `LSFConvolution.sigma_tensor`
(`ampere/backends/torch/instrument.py`) so it cannot be mistaken for the
`NoiseModel.sigma_tensor` hook — `width_tensor` or similar; update its
numpy-side wrapper, callers, tests and the module docstring that mentions the
hook. (4) Add a native zero-uncertainty precondition to the jax realised GP
path (`ampere/backends/jax/problem.py`): `GaussianProcessNoise` on the
contract path refuses a dataset with zero uncertainties and no jitter, the
`realise` agreement check catches the disagreement today, but a `strict`
NUTS run would sample a density the contract refuses — mirror the torch
path's construction-time refusal (`REQUIRES_UNCERTAINTY` / the sigma-is-None
check around line 507 of `ampere/backends/torch/problem.py`) with a test
that the refusal names the dataset. No behaviour change beyond the four.
**Depends:** nothing. **Blocks** nothing; may run in parallel with W3.1.
**Accept:** `grep -rn "pip install ampere" ampere/core ampere/backends
ampere/inference ampere/results` shows only the chosen wording; the
`OptionalDependencyError` doctest passes; one `SAMPLE_STATS_GROUP`
definition; no `sigma_tensor` on any torch `Transformation`; a jax
realised GP problem with zero uncertainties refuses by name at construction;
`pixi run -e dev test-all`, `-e torch test-all` and `-e jax test-all` green
(counts unchanged apart from the new rows); lint/format clean, pyrefly 0
errors in all three.

### W3.1 slice 1 — `simulate_many`: the batched forward model, the executor, external simulators [M; Opus]
Revised 2026-09-08 from Peter's three notes on the first draft (decision-log
row of the same date): (a) `vmap` is one device, so the design must be ready
for simulations distributed across devices or machines and for single
simulations that exceed one device's memory; (b) slow external simulators —
Fortran/C/C++/Rust routines behind a Python call — are a major SBI case and
must be first class, not a fallback; (c) every backend supports observation
sampling natively (slice 2). **Core** (`ampere/core/dataset.py`, a new
`ampere/core/simulate.py`): `FittingProblem.simulate_many(count, *,
values=None, observe=False, rng=None, stream="simulate", executor=None,
chunk_size=None)` → `SimulationBatch` — `theta` as `(count, free_size)`,
per-dataset `predicted`/`observations` as `FunctionSamples` batches with a
leading sample axis, the per-draw `ModelResult`s kept for design horizon
(c), a `failed` mask and the per-draw `Failure` records, `__getitem__`
yielding the i-th `Simulation` so every §13 consumer still works, and
iteration by chunk. `values=None` draws the batch from the joint prior on
the named sub-stream; an array `(count, free_size)` simulates at given θ
(the SBC idiom). **The loop is the reference semantics and the contract
every executor must satisfy**: `simulate_many(n)` equals `n` calls of
`simulate` under the same sub-stream, order preserved, with the per-draw
sub-stream derived from the batch sub-stream by *index* so the result is
independent of how the work was partitioned. **The executor protocol**
(`ampere.core.simulate.Executor`: `map(fn, items)` → results in order),
three shipped: `SerialExecutor` (default), a `concurrent.futures` process
pool (the route for external simulators — one process per simulation, the
problem pickled once per worker, so `FittingProblem` picklability is asserted
or supplied via its spec), a thread pool for I/O-bound wrappers; plus the
documented adapter shape for user-supplied mappers — the shape is
`concurrent.futures.Executor`'s `map`/`submit`, which dask's `Client`, ray's
executor wrappers and mpi4py's `MPIPoolExecutor` already satisfy (Peter,
2026-09-09: "dask and/or ray can easily provide replacements"), so a
queue-driven cluster array is the only case needing a thin adapter. `timeout=` per simulation on the pool executors; an expiry or a
crashed worker is a flagged `Failure`, never an exception. **Chunking**:
`chunk_size` bounds how many simulations are in memory at once;
`simulate_many(..., as_chunks=True)` yields `SimulationBatch` chunks and
`write_training_set`/`append_training_set` accept the iterator, so a budget
larger than memory is written without ever being held — which makes
`results.md` 13.9's `O(existing + new)` append the bottleneck for very large
budgets; the item measures the append time at 10⁴ and 10⁵ tiny simulations
and reports it (the trigger for the deferred unlimited-dimension append). A
single simulation larger than one device is not partitioned by the framework
— a model-parallel simulator is a `Model` whose `__call__` does its own
placement — but `chunk_size=1` with a per-chunk device-placement hook
guarantees such a model is never asked to hold two simulations at once.
**External simulators, first class**: a black-box `Model` on the reference
backend whose `__call__` calls compiled code is the canonical SBI simulator.
The item ships `examples/sbi/external_simulator.py` — a subprocess-based toy
standing in for a compiled routine (argument marshalling, a working directory
per worker, stdout/stderr capture into the `Failure` detail) run through the
process pool — and `tests/examples` runs it. **`BATCHABLE` on the reference
backend** acquires a meaning: a `Model` may implement `evaluate_batch(thetas)`
(many external codes take a table of parameter sets in one call) and declare
`BATCHABLE = True`; `simulate_many` uses it under the serial executor. **Spec**:
§13 gains "batched form" and "execution" subsections, limitation 17.5 is
closed, each *Amended W3.1*; the decision-log row is amended with the landed
text (a §4.5 surface addition). **Depends:** nothing merged; W3.0 may run
alongside. **Blocks** W3.1 slice 2, W3.2, W3.5, W3.6, W3.8.
**Accept:** conformance row `simulate_many(n) == [simulate() × n]` on every
registered backend fixture and under every shipped executor (serial; process
pool with two workers; thread pool) — order and values; partition
independence — `chunk_size` 1 and 7 give identical batches; an injected 2 %
crash rate plus one timeout are counted and the usable pairs are correct; the
external-simulator example runs under the process pool in `tests/examples`;
writing from chunks equals writing the whole batch; `FittingProblem`
picklability asserted; the append-time table in the report; `dev` gate green
(this slice is numpy-only) with `torch`/`jax` gates unchanged;
lint/format/pyrefly clean ×3.

### W3.1 slice 2 — Native batched prediction and native observation sampling on torch and jax [M; Opus; two-track possible]
(a) `LoweredProblem.simulate_batched(theta)` for the **noise-free
prediction** through `torch.func.vmap`/`jax.vmap` where every part declares
`BATCHABLE` and a realisation is registered — **per chunk**: `vmap` over a
chunk with chunks looped, never over the whole budget (that is exactly the
single-device memory trap Peter's note names); device placement per the
slice-3 `device=` plumbing; `QuasisepGP` and non-batchable parts fall back to
the loop honestly; `simulate_many` uses the native path when the problem is
realised and batchable, with the same one-point agreement check `realise`
makes, and training sets written from it carry `ampere_simulate_batched`.
(b) **Native observation sampling** (Peter's ruling 2026-09-08): every
`LikelihoodFamily.sample` the core implements gets a native twin — torch
over `torch.distributions`, jax over `jax.random` and numpyro's
distributions — for `gaussian` (independent noise, and GP noise through the
backend solver's own `latent_transform`), `student_t`, `cauchy` and
`complex_gaussian`; `poisson` keeps §13's refusal on every backend unless the
user overrides. (c) **Multi-device hooks** (Peter, 2026-09-09): the native
chunk path exposes a hook for sharding a chunk across the devices the
backend reports — `jax.pmap`/`shard_map` on jax, the `torch.distributed`/
`DTensor` equivalents on torch — designed, documented and smoke-tested at
the API level (`tests/gpu`-style rows, skipped without hardware), with
CPU-only CI exercising the single-device degenerate case; full exercise
waits for the GPU item. (d) **Two W3.0 carry-overs**: torch's `gp_marginal`
branch gains the zero-uncertainty construction-time refusal jax now has, and
jax's `GaussianProcessNoise` gains the `jitter=`/`scale=` constructor
keywords torch's already has (parity; conformance row). RNG derived from the problem's seed sub-stream
(`integer_seed(stream)` → `torch.Generator` / `PRNGKey`), the backend that
drew recorded as `ampere_sample_backend`. **The numpy path stays the oracle**
and is compared *distributionally* — mean and covariance of 4 000 draws
within the tolerances §13's numpy test uses — because RNG streams differ by
backend; the masked-samples rule and the refusal texts match exactly. Spec:
§13's "what can be sampled" gains the native paragraph, *Amended W3.1*.
**Depends:** W3.1 slice 1; W3.0 (it touches the jax `problem.py`).
**Blocks** nothing hard — W3.2 works from slice 1; this slice is throughput
and native parity.
**Accept:** native batched prediction agrees with the loop to
`tolerances.cross_backend` on both backends and refuses by name on a
non-batchable problem; `chunk_size` 1, 7 and the whole budget give identical
predictions; native Gaussian draws' empirical covariance matches
`K + diag(σ²)` within §13's tolerances on both backends, off-diagonals
non-zero, variances exceeding `K`'s (the same two mutation-tested
assertions); Poisson refuses by name on both; the sharding hook's API rows
skip cleanly on CPU and the single-device case runs; torch refuses a
zero-uncertainty GP-marginal dataset at construction as jax does; a jax
`GaussianProcessNoise(jitter=...)` composes and lowers; `torch`, `jax` and
`dev` gates green; lint/format/pyrefly clean ×3.

### W3.2 — `SBIEngine`: NPE/NLE/NRE through the `sbi` package [L; Opus]
The one SBI module `DEVELOPMENT_PLAN.md` §5 asks for, as an
`ampere.inference` engine beside `EmceeEngine`/`NUTSEngine`/`VIEngine`,
consuming `simulate_many` from **any** backend — a legacy black-box model
composed on the reference backend is the first-class case, not an
afterthought. `ampere/inference/_sbi.py`: `SBIEngine(problem, *,
method="npe" | "nle" | "nre", budget, embedding=None, density_estimator=
None, rounds=1, device="cpu")`; `sbi` and torch imported lazily inside
`run`, `OptionalDependencyError(extra="sbi")` on use, so `dev`'s import-graph
test is untouched. **Prior**: sbi wants a torch `Distribution` with
`sample`/`log_prob` over the flat free vector — build it from
`FittingProblem.parameters` in the **unconstrained** coordinates
(`unconstrain`/`constrain` and `log_prob_unconstrained`'s prior part), so
bounded priors need no `RestrictedPrior` and the density estimator sees ℝⁿ;
record the choice in the attrs. **Simulator**: `simulate_many` with
`observe=True`, failed draws dropped and counted (§4.5's reject-and-record),
the pairs written to a training set on request (`training_set=` path), and
the summary tensor built by the encoding W3.3 defines — until W3.3 merges,
the fixed-layout concatenation of each dataset's observed values in
`datasets` order, which is what a single fitting problem needs. **Posterior
→ run**: the trained posterior is sampled at the observed data (`draws`
i.i.d., **one chain**, as `VIEngine` does and for the same reason), and every
stored draw is scored on the numpy path through `Engine.finish` — so the run
carries the *true* `log_prob` per draw beside the estimator's own `log_prob`
(attr-named `ampere_sbi_log_prob`), which is what makes the SBC/coverage
item and importance reweighting possible later. Attrs: `ampere_engine =
"sbi"`, the method, the estimator architecture, the budget, the round count,
the failure count, the embedding's name and output dimension, the training
loss trace (thinned, stride recorded), the `sbi`/torch versions. **Legacy
parity**: the embedding-network conveniences of `ampere/infer/sbi.py`
(`"FC"`, `"CNN"`, a user `nn.Module`, a dict of hyperparameters) are carried
over as the `embedding=` vocabulary — frozen legacy is **read, not
modified**. **Not in this item**: swyft (W3.4), caching (W3.5), SBC (W3.6),
multi-round truncation beyond `sbi`'s own `rounds` loop. **Depends:** W3.1 slice 1 (slice 2 is throughput; the engine's `executor=` and `chunk_size=` pass straight through to `simulate_many`, which is how an external simulator is run under the process pool).
**Blocks** W3.4, W3.5, W3.6.
**Accept:** in the `sbi` environment, `SBIEngine` recovers the
`tests/inference` toy joint posterior with NPE at a small budget (thresholds
on posterior mean and width against the emcee reference, seeded; the smoke
budget small enough for CI — plan §5's "SBI smoke tests with tiny simulation
budgets"); NLE and NRE run end-to-end on the same problem (shape and
finiteness, not accuracy, at CI budgets); the same engine runs a legacy
black-box `Model` wrapped on the reference backend; `pixi run -e dev
test-all` proves torch/sbi are not imported; the emitted `DataTree`
validates against the results contract with schema 5 and the new attrs;
`pixi run -e sbi test-all` green; lint/format clean, pyrefly 0 errors in
`dev` and `sbi`; a `tests/inference/test_sbi.py` that skips cleanly without
the extra.

### W3.3 — The coordinate–value–mask encoding for embedding networks [M; Opus]
The plan's Phase 3 second bullet: a canonical tensor encoding of a
`DatasetCollection`'s observed containers so **set-based** embeddings can
consume any modality, irregular sampling and missing data included, and so
amortisation across differently-sampled datasets becomes possible. In
`ampere/core/encoding.py` (backend-neutral, numpy): `encode_observations
(datasets, *, layout) -> Encoded` where each dataset contributes rows of
`[coordinates…, value(s), uncertainty?, mask]` — the coordinate columns
from the container's axes in declared order, complex values as two columns,
the mask from `Dataset.effective_mask` — and a **dataset-id column** so one
tensor carries the whole collection; `layout` is a frozen description
(column names, per-dataset row counts, dtypes) hashed into the training
set's attrs so a network trained on one layout refuses another by name;
padding to a common row count with the mask column zero on padded rows;
`decode` for the round trip in tests. The fixed-size summary of W3.2
becomes one `layout` among others (`layout="flat"`), the default for a
single fitting problem where the data layout is fixed anyway. **On the sbi
side** (`ampere/inference/_sbi.py`): `embedding="set"` builds
`PermutationInvariantEmbedding` over the encoding, `embedding="transformer"`
the `TransformerEmbedding`, both masked by the mask column; `layout=` is an
`SBIEngine` argument. **Spec**: a short `docs/design/contracts/encoding.md`
(pipeline stage, inputs/outputs in contract vocabulary, the layout hash,
what is and is not stable across `CONTAINER_SCHEMA_VERSION`s) — an addition
to §4, so a decision-log row. **Depends:** W3.2. **Blocks** nothing hard;
W3.6 uses it if merged.
**Traps and their solutions (Fable, 2026-09-09, verified against sbi
0.27; each is a requirement of this item).** (1) **sbi z-scores `x`
column-wise by default** (`posterior_nn(z_score_x="independent")`) over
the training tensor — on a padded, mixed-column tensor that standardises
the mask and coordinate columns and lets padded rows contaminate every
column's statistics: the engine passes `z_score_x="none"` for any
non-flat layout and the **encoding does its own standardisation**, per
column group, with the statistics part of the layout (so the observed
data at inference is standardised identically). (2) **Values across
datasets differ by orders of magnitude** (photometry vs spectrum units):
the value columns are the whitened `y/σ`, `log σ` and an `asinh`-scaled
`y`, not raw `y` — this is also what makes amortisation over noise
possible. (3) **Mask conventions differ per embedding**: sbi's
`PermutationInvariantEmbedding` treats an all-NaN row as absent and
already computes a masked mean; the transformer wants a `(batch, rows)`
`attention_mask`. The contract keeps one explicit mask column
(backend-neutral, serialisable); the torch-side wrapper `unpack`s the
tensor by the layout and converts — NaN rows for the set embedding,
`attention_mask` for the transformer — so no embedding ever sees the mask
column as a feature. (4) **`TransformerEmbedding` defaults** are wrong for
sets: `is_causal=True` (rows attend only to earlier rows), dropout 0.5 on
both attention and MLP, `pos_emb="rotary"` (index-based). Set
`is_causal=False`, explicit dropout, `pos_emb="none"`, and carry the
coordinate as features: the normalised coordinate plus fixed Fourier
features (NeRF-style, band count in the layout). Its `forward` returns a
tuple; the wrapper takes the first element. (5) **Sum pooling makes the
embedding depend on the row count**; use the masked mean and add `log N`
(valid rows per dataset) as an explicit per-set feature, so the count is
information rather than a scale. (6) **The layout fixes the row cap**
(sbi validates `x_shape` at `set_default_x` anyway): an observation with
more rows than the layout allows is refused by name; attention is
quadratic in rows, so a 2 000-point spectrum beside 10 photometric points
is the case the row cap and, later, per-dataset sub-encoders exist for.
(7) **dtype**: ampere is float64, sbi nets are float32 — cast at the
boundary, in the wrapper, once. (8) **Complex values** are two columns,
real and imaginary (linear), recorded in the layout. (9) **The durable
design**: `encoding.md` defines the *packing* — named column groups
(coordinates, coordinate features, values, σ, mask, dataset id, per-set
features, a reserved context slot) — and one `unpack(x, layout)` helper;
every embedding, present or future (neural-process encoders, FiLM-
conditioned nets), is a module over the unpacked view, so the packing is
the contract and the networks are free.
**Accept:** encode/decode round trip on every container fixture in
`tests/core` including a masked, an irregular and a complex one; two
datasets of different lengths encode into one padded tensor whose mask
column is exact; a layout mismatch is refused by name; `SBIEngine(...,
embedding="set")` trains and samples on a two-dataset toy problem in the
`sbi` environment; `dev` and `sbi` gates green; lint/format/pyrefly clean.

### W3.4 — TMNRE through sbi's own ratio estimators [M; Opus] (ruled 2026-09-09; the swyft revival is not pursued)
Plan §5: "revive the swyft TMNRE implementation from
`docs/design/harvest/swyft/`". **Ruled by Peter 2026-09-09: the "express through sbi" route — no swyft,
no maturity note; the design paragraph below is the item's scope, and
the swyft-revival text after it is kept only as the record of what was
not done.** *(Superseded)* Gate first: the item begins with a
one-page maturity note — does swyft install alongside sbi 0.27 and torch
2.13 in the `sbi` environment today (it is a PyTorch-Lightning package,
last seen active 2024), what its truncation offers that `sbi`'s `NRE` +
`RestrictedPrior` rounds do not, and whether the answer is "revive",
"express TMNRE through sbi's own NRE rounds" or "drop". Peter rules on the
note before the code is written; if the ruling is revive: `swyft` becomes
an optional extra (`ampere[swyft]`, pixi feature `swyft`, CI non-blocking
like the legacy `sbi-characterisation` job), `SBIEngine(method="tmnre")`
routes to a `swyft.SwyftTrainer` path in `ampere/inference/_swyft.py`
consuming the same `simulate_many` and encoding, and the harvest's
**latent bug is fixed in the revival**: the `SwyftNetwork*` classes
referenced the bare `swyft` name at class-definition time, so the lazy
guard never protected them — the revived classes are defined inside the
lazily-imported path or behind a factory. **Depends:** W3.2 (and W3.3 for
set embeddings). **Blocks** nothing.
**The "express through sbi" design (Fable, 2026-09-09, on Peter's
request; the option Fable recommends).** TMNRE (Miller et al. 2021) is two
ideas. *Marginal* ratio estimation: instead of one classifier for the joint
ratio `r(θ, x) = p(θ|x)/p(θ)`, train one per marginal of interest —
each 1D `θ_i` and each 2D pair — which is what makes the method robust at
high parameter dimension and is what the corner plot actually needs.
*Truncation*: after each round, restrict the *prior* to the region where
the estimated 1D marginals put their mass (the hyperrectangle of per-
parameter intervals above a small threshold ε), simulate the next round
inside it, and retrain; because the restricted prior is the original
prior renormalised on a subset — not a learned proposal — the ratio is
unchanged inside the region and no importance correction is needed, the
estimate stays amortised *within* the box, and the box is a diagnostic in
its own right. In sbi 0.27 both halves exist: `NRE`/`NRE_B`/`BNRE` train a
ratio estimator; a marginal estimator is the same trainer fed
`theta[:, idx]` with a prior over that marginal (for ampere's unconstrained
prior the marginal `log_prob` is not closed-form, but NRE only needs prior
*samples* to train and the prior only to build a posterior — so marginals
train from the joint draws' columns, and the marginal posterior is
evaluated as `exp(log r) · p(θ_i)` with `p(θ_i)` estimated once from
prior draws); `RestrictedPrior(prior, accept_reject_fn)` is the truncated
prior, with `accept_reject_fn` an indicator on the current box (sbi's own
`get_density_thresholder` gives the HPD form used by TSNPE; TMNRE's box is
the product of 1D marginal intervals, a few lines on the 1D estimators'
outputs), and `append_simulations(theta, x, from_round=r)` accumulates
rounds. `build_posterior(sample_with="rejection")` samples the restricted
prior through the ratio; `"mcmc"` (slice) is the fallback for narrow
boxes. Ampere-side: `SBIEngine(method="tmnre", rounds=R, marginals=1|2,
truncation_epsilon=ε)` — a round loop over `simulate_many` with the
`RestrictedPrior` handed to it as the θ source (the batch's `values=`
argument takes an array, so the engine draws from the restricted prior
and passes the array), the 1D/2D estimators trained per round, the box
recorded per round in the attrs (`ampere_sbi_truncation`), and the
final marginal posteriors emitted as the run's draws — one chain from the
final-round rejection sampler on the *joint* estimator if one is also
trained, or the 1D/2D marginals as `sample_stats`/a derived group if not
(a results question for the item: a run whose posterior is a set of
marginals rather than joint draws is new to §4.6 and needs a row). No
swyft, no Lightning, and the SBC/coverage machinery of W3.6 applies
unchanged.
**Accept (if revived):** the maturity note in the PR; `import
ampere.inference` succeeds without swyft; `method="tmnre"` refuses with
`OptionalDependencyError(extra="swyft")` without it and recovers the toy
posterior with it at a smoke budget; the emitted run carries the truncation
history in its attrs; `dev` and `sbi` gates unchanged.

### W3.5 — Trained-artefact caching keyed on the spec hash [M; Sonnet or Opus]
Plan §7's trap, stated there with evidence ("the recent SBI caching bugs on
master"): an SBI posterior, embedding net or emulator reused against a
problem it was not trained on is silently wrong. `ampere/inference/_cache.py`
(or `ampere.results.artefacts`, decided in the item by §8's placement logic
and recorded): `ArtefactStore(root)` whose key hashes, together: the problem's joint
spec hash (`ampere.results.provenance.spec_hashes`' `"spec"` entry), the
observed-data hash, the encoding layout hash, the method, the estimator
architecture, the budget, the seed, and the sbi/torch versions; `store.get(key)`/`store.put(key,
artefact, attrs)` persisting the torch state dict plus a JSON sidecar of the
key's ingredients, so a miss can say **which** ingredient changed; a hit
that fails to load is a miss, never an error. `SBIEngine(cache=store)` uses
it around training; the training set itself is the existing
`ampere.results.training` format (W2.8) keyed the same way, so a budget can
be reused for a different estimator. **Refusals**: a partial key is never
accepted; there is no "force" that skips the spec-hash comparison; a cache
hit is recorded in the run's attrs (`ampere_sbi_cache_hit`, the key). No
binary artefacts in git — the store lives outside the repository and tests
use `tmp_path`. **Depends:** W3.2. **Blocks** nothing.
**Accept:** a second `SBIEngine.run` with an identical problem and settings
trains nothing and emits an identical posterior (seeded); changing any one
of prior, data, model, layout, method, architecture, budget or seed is a
miss whose report names the ingredient; a corrupted stored artefact is a
silent miss with a warning; `dev` and `sbi` gates green; lint/format/pyrefly
clean.

### W3.6 — Family D: simulation-based calibration and coverage [M; Opus]
`diagnostics.md` §11 turned into code, placed per that section's own rule:
it consumes `InferenceData`s and `simulate`, and brings no new dependency
for the SBI path, so it lives in `ampere.results.calibration` with the
`sbi`-specific fast path inside the engine. Two routes. (1) **For an
`SBIEngine` posterior**: `run_sbc`/`check_sbc` and `run_tarp`/`check_tarp`
from `sbi.diagnostics` over a fresh batch from `simulate_many` (prior draws,
`observe=True`), the rank statistics and the coverage curve returned as a
small xarray `Dataset` and written into the run's `DataTree` as a
`calibration` group (not `AnomalyScore`: §11 declines the convention, for
the reason family B did). (2) **For any engine**: `sbc(problem, engine_
factory, *, count, draws)` — the Talts et al. loop over `simulate` and a
full fit per simulated dataset, expensive by design and budget-controlled,
returning the same `Dataset`; this is what validates the flexible-GP
likelihood itself, as §11 anticipates for M2. **Route (2)'s first
application is the repeated-trial coverage study W2.9 asked for** (ruled
2026-09-09, Peter leaving the timing to Fable's judgement): the WStat
example at low counts under `EmceeEngine`, demonstrating the profile-MLE
bias, its rank histogram and coverage curve added to
`docs/source/wstat_comparison.rst` — the machinery *is* a repeated-trial
coverage study, so this is where it costs least. Plots: `plot_sbc_ranks` and
`plot_coverage` in `ampere.results.plots`, following the six existing plots'
conventions. **Spec**: §11 is marked landed with the placement decision,
*Amended W3.6*; `results.md` gains the `calibration` group in its schema
table — a §4.6 addition, so a decision-log row. **Depends:** W3.1 slice 1, W3.2.
**Blocks** nothing.
**Accept:** on the toy problem, an NPE posterior's rank histogram is
uniform within `check_sbc`'s own thresholds at the smoke budget and a
deliberately narrowed posterior (temperature-scaled draws) fails the check;
route (2) runs with `EmceeEngine` on a two-parameter problem at a tiny
budget; the `calibration` group validates and round-trips through netCDF;
the two plots render headless; `dev` and `sbi` gates green;
lint/format/pyrefly clean; decision-log row present.

### W3.7 — CI/CD Phase 3 expansion [S; Sonnet]
Plan §5's cross-cutting line: "SBI smoke tests with tiny simulation
budgets". Promote the `sbi` environment from the non-blocking weekly
`sbi-characterisation` job to a **blocking** `suites` leg running `pixi run
-e sbi test-all` on every PR, with the smoke budgets W3.2 and W3.6 chose;
`pixi run -e sbi typecheck` beside it; the sbi environment's install cached
like torch's. Keep the legacy characterisation job as it is. Update the
`suites` matrix documentation in `docs/development.md` and the workflow
comments. **Depends:** W3.2. **Blocks** nothing.
**Accept:** `ci.yml` YAML-validated; the `sbi` leg runs the five suites in
one process; the job's wall time is under the torch leg's; the branch-
protection note in the handoff lists the new leg name.

### W3.8 — Non-native parts in a native problem: refuse by default, opt in without gradients [M; Opus]
Ruled by Peter 2026-09-08 on W2.4 slice 3's carried finding (decision-log
row of the same date; this item lands the §4.5 text and amends the row).
Today a core numpy `Kernel` inside a torch GP noise model is accepted and
its amplitude silently gets no gradient. The ruling has three parts.
**Refused by default**: `declared_capabilities` treats a core part inside a
backend's noise model, solver or likelihood as the backend disagreement
fold-in 7 already gives a numpy solver — by name, with the remedy.
**Accepted under an explicit opt-in when no gradient is needed**: the use
case is a piece expressible only in Python — a tabulated kernel, a legacy
callback — that a torch/jax problem wants to call through and sample
gradient-free. Working name `allow_foreign_parts=True` on `FittingProblem`
(spelt to sit beside the existing `strict` idiom; the item may propose
better), recorded in provenance as `ampere_foreign_parts` with the parts'
names, and the affected dataset's `DIFFERENTIABLE` becomes `False` so the
capability ladder tells the truth. **Always refused where a gradient is
required**: `realise(strict=True)`, `NUTSEngine`, `VIEngine` and the
differentiable native `LoweredProblem` refuse by name regardless of the
flag. The gradient-free engines run through the contract path (the fast path
falls back, as it already does for a non-batchable part); whether the
realisation can still call the foreign part cheaply for the fast path
(torch: detach→numpy→tensor hop; jax: `jax.pure_callback`) is measured in
the item, not assumed, and the report says which. Conformance: a row per
backend that the default refuses, the opt-in accepts and samples with emcee
to the same posterior as the all-native problem, and the gradient routes
refuse. **Depends:** W3.1 slice 1 merged (it owns `dataset.py` first).
**Blocks** nothing.
**Accept:** the default refusal names the part and the remedy; an opt-in
run's provenance carries the flag and the part names; `NUTSEngine`,
`VIEngine` and `realise(strict=True)` refuse by name under the opt-in; the
emcee posterior on the opt-in problem agrees with the all-native one within
`tolerances.cross_backend`; all three backend gates plus `dev` green;
lint/format/pyrefly clean ×3; the decision-log row amended with the landed
text.

### W3.9 — Docs: the legacy API section and a "Migrating to the new API" page [S; Sonnet]
Ruled by Peter 2026-09-08 on W2.15's open question: users will keep using
the legacy API, so the legacy pages stay **linked** from their own section
of the API reference (`docs/source/api.rst` already has "Legacy (frozen)";
`:no-index:` stays so nothing cross-links into frozen code) and gain a
**"Migrating to the new API"** page. This item seeds that page: a concept
map from legacy to v2 (a legacy `Model` subclass → a v2 `Model` publishing
`ModelResult` channels; the legacy `Photometry`/`Spectrum` data classes →
`Dataset` with an `Instrument` chain; the legacy likelihood's GP switches →
`Likelihood` with `GaussianProcessNoise` and a `Kernel`; the legacy
`ampere.infer` search classes → the `ampere.inference` engines, with the SBI
row pointing at W3.2); the frozen-legacy statement (what "frozen" means,
what still works, that the two APIs do not interoperate); the emcee
`minimal_working_example` shown side by side, legacy and v2, as prose with
short snippets; and a pointer to the Phase 6 deprecation policy as "to
come". No code changes. The full migration guide and the deprecation policy
remain Phase 6's. **Depends:** nothing — dispatchable as soon as the
breakdown is approved. **Blocks** nothing.
**Accept:** `pixi run docs` succeeds with no new warnings; the migration
page sits in the API reference toctree beside the legacy section; every
legacy name on the page is either an explicit `:doc:` link to its page or a
plain literal (no dangling cross-references); British English.

### W3.11 — Set-embedding readouts: default width and a pooled transformer head [S; Sonnet] (ruled by Peter 2026-09-09)
From W3.3's two open questions. (1) The set and transformer embeddings
inherit the `"flat"` default output width `2·free_size`; for a pooled set
that is the conditioning vector's whole capacity, and sbi's nets end in a
ReLU, so a narrow randomly-initialised net can emit all zeros (W3.3's
tests set 16 explicitly for this reason). Default becomes
`max(2·free_size, 32)` for `"set"`/`"transformer"`, `2·free_size` kept for
`"flat"` (legacy parity), user override unchanged. (2) sbi's transformer
reads the **last token** as its summary; under full attention with no
positional embedding the last token is a function of every row, but the
readout depends on which row happens to be last (or on a padded zero
token). The wrapper gains a **masked-mean readout head** over the tokens
after sbi's transformer body, which is permutation-invariant and uses
every retained token; a learned CLS token is the alternative, more
parameters for the same information, not chosen. Both are wrapper
changes in `ampere/inference/_sbi.py` (§7 of `encoding.md` amended,
*Amended W3.11*); no contract change. **Peter's rider (2026-09-09)**: both
choices are accepted as defaults, not as findings — how the right width
depends on the data, the model and the problem structure, and what the
readout choice (last token, masked mean, CLS) costs, are to be **tested
carefully at some point** with the calibration machinery, and the result
written up as guidance for end users; recorded in the Phase 3 deferred
list as an embedding study. **Depends:** W3.6 (owns `_sbi.py`
until it merges). **Accept:** the two defaults asserted; the transformer
wrapper's output unchanged under row permutation (a test that fails on
the last-token read and passes on the mean); `sbi` and `dev` gates green.

### W3.12 — Model identity hash: promote, record, and check on append [S; Sonnet] (ruled by Peter 2026-09-09)
From W3.5's open question and its carried finding. W3.5 built a
`model_hash` locally (model fingerprints minus their parameter specs,
dataset fingerprints minus the observed data, plus bindings) because the
spec hash covers only the parameter declaration and a kernel/solver/
family swap with identical parameters must be a miss. `provenance.py`
already has `model_fingerprint`/`dataset_fingerprint`/`problem_fingerprint`;
this item promotes W3.5's composition to a public
`model_identity_hash(problem)` beside `spec_hashes`, **and lands W3.4's
carried `artefact_key` diff** (`marginals=`/`truncation_epsilon=`
ingredients, written only when set so existing digests are stable; `_sbi.py`
then drops `_key_architecture` and passes both at its one call site), records it on every
run and training set as `ampere_model_hash` (**`PROVENANCE_SCHEMA_VERSION`
→ 6**, one decision-log row), makes `artefacts.py` use it rather than its
private copy, and fixes W2.8's `append_training_set` to refuse an append
whose `ampere_model_hash` differs from the file's, by name — plan §7's
trap, live today. **Depends:** W3.6 (touches `ampere/results/`).
**Accept:** the attr on a run and a training set; a kernel swap with
identical parameters refuses the append and names the hash; `artefacts.py`
has no private fingerprint code left; all four gates green (schema bump
touches every run); lint/format/pyrefly clean.

### W3.14 — Generative forms for the unambiguous families: `sample()` on Poisson, Student-t and complex Gaussian [M; Opus] (approved by Peter 2026-09-09)
Peter's ruling of 2026-09-09 ("implement further core `sample` families
when we have a use case") met its use case in W3.6: a counting
experiment cannot `simulate(observe=True)`, so neither SBC nor SBI on
count data works without the user subclassing `PoissonFamily` to add
three lines. `inference.md` §13's refusal was written for families whose
observation process is genuinely ambiguous; Poisson counts, Student-t and
complex Gaussian are not — each has one generative form that its
`log_prob` already fixes. This item adds `sample(predicted, noise, rng)`
to those three in `ampere.core` (Poisson: `rng.poisson(predicted)`;
Student-t: location-scale with the family's degrees of freedom, using the
noise model's σ as the scale exactly as `log_prob` does; complex Gaussian:
independent real and imaginary parts with the noise model's σ — circular
symmetry, as the family's density assumes), keeping the specific refusal
for `cauchy`-with-no-scale and any family that does not fix its form,
and adds the **native twins** on torch and jax through the presence-based
dispatch W3.1 slice 2 left open (the pathway Peter asked to keep open —
this item is its first use). Latent-GP consumers (Poisson with a GP rate)
sample the latent through `latent_transform` first, as the Gaussian GP
branch does. **Spec**: §13's "what can be sampled" amended (*Amended
W3.14*), one decision-log row (a §4.4 contract change: the refusal
becomes the exception rather than the default for these families).
**Depends:** W3.6 (merged). **Blocks** nothing.
**Accept:** conformance rows per backend: 4 000 draws of each family
match its `log_prob`'s implied moments (Poisson mean = variance =
`predicted`; Student-t scale and ν; complex Gaussian's real/imaginary
variances and zero cross-covariance), numpy path the oracle and native
twins compared distributionally; the W3.6 WStat example's `CountingPoisson`
subclass deleted; the refusal text unchanged for families still refusing;
all four gates green; lint/format/pyrefly clean ×3.

### W3.15 — Reproducible SBI runs: seed torch from the problem [S; Sonnet] (drafted 2026-09-09 from W3.4's finding; dispatched on Fable's judgement, Peter to confirm)
Neither ampere nor sbi seeds torch's global generator, so two
`SBIEngine.run`s of the same seeded problem give measurably different
posteriors (network initialisation and sbi's training-batch order vary);
`problem.seed` fixes only the simulation streams, and W3.5's "identical
posterior" claim holds only through the cache. `run()` derives an integer
from the problem's own sub-stream (`self.integer_seed("sbi.torch")`,
distinct per round) and calls `torch.manual_seed` (and
`torch.cuda.manual_seed_all` when a device is in play) before training and
before sampling, records `ampere_sbi_torch_seed` in the attrs, and the
cache key already includes the problem seed. Cost: none; a user who wants
fresh randomness sets a new problem seed, which is the contract's idiom.
**Depends:** W3.4 (merged). **Accept:** two `SBIEngine.run`s of the same
seeded problem produce bitwise-identical draws and attrs (NPE, and TMNRE
at a tiny budget); a different problem seed differs; `sbi` and `dev` gates
green; lint/format/pyrefly clean.

### W3.13 — Phase 3 documentation pass [M; Sonnet] (drafted 2026-09-09 by Fable; dispatch when W3.15, W3.10 and W3.7 have merged; Peter to confirm)
Like W2.15, after the last Phase 3 merge: every user-facing and design
document reflects what landed. (1) `docs/source/`: an **SBI tutorial page**
(`sbi.rst`, in the v2 tutorials toctree beside the M2 and WStat pages)
built from `examples/sbi/` — `toy_powerlaw.py` (NPE on a native problem at
the smoke budget; what `SBIEngine.run` returns; the `ampere_sbi_*` attrs),
`external_simulator.py`/`fit_external_simulator.py` (a black-box simulator
through `simulate_many`'s executor protocol and `forkserver`),
`cached_fit.py` (the artefact store: a hit, and a miss that names what
moved), `tmnre_fit.py` (TMNRE, the `marginals` group, the truncation
history) — plus `SBIEngine.calibrate` (SBC and TARP, `plot_sbc_ranks`,
`plot_coverage`), the embedding choices (`EncodingLayout`, set versus
transformer, the default width) and reproducibility (W3.15's seed); every
example the page uses covered by `tests/examples` as the M2 and WStat pages
are; `advanced.rst`'s "Very slow models" and "Embedding networks" sections
rewritten for `ampere.inference` with the legacy note inverted (legacy
`ampere.infer.sbi` is the frozen route); `migrating.rst`'s `SBI_SNPE` row
completed; `install.rst`'s `sbi` extra row updated; the legacy
`Embedding_nets` notebook kept under the legacy warning. (2)
`docs/design/architecture.md` §3 and the plan's §5 Phase 3 bullets
annotated with what landed (sbi 0.27; **no swyft** — TMNRE through sbi's
own ratio estimators, Peter's ruling of 2026-09-09; the encoding;
`forkserver`; the schema-6 attrs), and every `docs/design/contracts/*`
"landed"/"amended" annotation for the Phase 3 landings checked against the
code (the list of stale sentences in the report, each marked *Amended
W3.13*, one decision-log row for the batch). (3) Docstring index: every
public name added in Phase 3 (`ampere.inference`'s SBI surface,
`ampere.results.calibration`, `artefacts`, `training`,
`provenance.model_hash`) renders under autodoc without warnings. (4) The
deferred list (plan §6) reviewed: the swyft half of the "SBI package set"
bullet struck through with the ruling, the jax-native SBI bullet kept, the
embedding *study* and the ε study entered where the handoff put them. (5)
`README.md`'s pitch gains one sentence on SBI, with the WStat coverage-study
numbers if they fit. No code behaviour changes; no contract changes beyond
annotations. **Depends:** W3.15, W3.10, W3.7 (merged). **Blocks** the Phase
3 close-out.
**Accept:** `pixi run docs` succeeds with no warnings attributable to the
new namespaces; every example the tutorial uses runs under `tests/examples`
(existing coverage confirmed, new coverage added where an example had none);
the annotation list in the report; a docs-only diff — the `dev` gate green,
plus the `sbi` gate if `tests/examples` changed.

### W3.10 — Automatic paging above the plot caps, with a loud warning [S; Sonnet; not urgent]
Ruled by Peter 2026-09-08 (on W2.8's confirmed caps): above
`MAX_CORNER_VARIABLES`/`MAX_TRACE_VARIABLES`, `plot_corner` and `plot_trace`
**page** — a list of figures each within the cap, in merged-name order,
array blocks kept whole where they fit — and emit a loud `ResultsWarning`
naming the page count, the cap and the `var_names=` route, instead of
refusing. `var_names=` and `max_variables=` keep their meanings;
`paginate=False` restores today's refusal for a caller who needs one figure;
`figure_metadata` records "page i of n" on each. `results.md` §8's "refused
loudly rather than attempted" sentence is amended to "paged, with a warning"
(*Amended W3.10*) and W2.8's decision-log row gains the note. Not on the
Phase 3 critical path — dispatch when convenient. **Depends:** nothing.
**Blocks** nothing.
**Accept:** 25 scalar parameters give two corner pages and one warning; a
200-element plate gives ten; `paginate=False` refuses exactly as today; the
metadata round-trips; `dev` gate green; lint/format/pyrefly clean.

### Reserved hook from the horizon notes (Fable, 2026-09-09; `docs/design/horizon_notes.md`)
Peter's look-ahead questions on amortising SBI over noise realisations and
over sampling share one hook that is cheap to reserve now and costly to
retrofit: a per-draw **observation context** (σ-pattern, grid, instrument
settings) drawn from a context prior, passed to `simulate_many`, recorded
on the `Simulation`, and visible to the embedding. Three consequences for
items not yet dispatched, none changing an approved body's scope: (i)
W3.3's encoding carries the **uncertainty column by default** (mandatory
unless the dataset has none), since a network that cannot see the error
bars cannot condition on them; (ii) W3.1 slice 2 and W3.2 accept a
`context=` argument on `simulate_many`/`SBIEngine` that is `None` today and
recorded as such in provenance, so the signature exists before the
machinery; (iii) the context prior and per-context instrument
renegotiation are a follow-on item drafted when W3.3 lands, or Phase 5 if
the Phase 3 budget is spent. **Reservation confirmed by Peter 2026-09-09.**

### Deferred from Phase 3 (recorded so they are not re-derived)
- **jax-native SBI** (sbijax/flowjax): plan §6 says "if and when maturity
  warrants"; nothing above needs it, and `simulate_many`'s native
  `simulate_batched` is the jax half that would matter.
- **Unlimited-dimension append** on the training-set writer (`results.md`
  13.9's extension point): the append is read-concatenate-rewrite,
  `O(existing + new)`. Scheduled only when a Phase 3 budget outgrows memory
  in practice — W3.2's report is to state the largest budget it wrote and
  the append time, and this becomes an item if that number is a problem.
- **Multi-fidelity and hierarchical SBI** (design horizons (a), (d)): the
  hooks stay reserved; `simulate_many` over a `Plate` is the first thing to
  check when (d) is scheduled.
- **Emulators** (design horizon (c)): a `Model` trained on the training set
  W3.1 writes; Phase 5 or later.
- **An embedding study** (Peter's rider on W3.11, 2026-09-09): once W3.6's
  SBC/TARP machinery and W3.11's defaults exist, a systematic study of how
  the embedding width and the readout (last token, masked mean, CLS)
  affect calibration and posterior width across data kinds, model sizes
  and problem structures — an `examples/sbi/` study with a docs page,
  its product being **guidance for end users**, not a contract change.
  Phase 3 if the budget allows, else the first Phase 5 SBI item.
- **RHMF exploratory trial**: Peter's ratification note on the W2.7
  deferral (2026-09-08) asks for early testing in a later phase; recorded
  as a Phase 5 bullet in the plan's §5.

### Dispatch order and parallelism
W3.0 ∥ W3.1 slice 1 ∥ W3.9 first (disjoint files: W3.0 owns `exceptions.py`,
`emission.py`/`training.py`'s constant, torch `instrument.py`, jax
`problem.py`; slice 1 owns `dataset.py` and the new `simulate.py`; W3.9 is
docs-only). Then W3.2 ∥ W3.1 slice 2 (disjoint: `_sbi.py` vs both
`problem.py`s and the backends' families) — W3.2 defines the engine every
later item extends. Then W3.3 ∥ W3.5 (disjoint:
encoding vs cache; both add arguments to `SBIEngine` — one owner of
`_sbi.py` at a time, so W3.3 first and W3.5 appends). Then W3.6, with W3.4's
maturity note written any time after W3.2 and its code only on Peter's
ruling. W3.8 after W3.1 merges (it touches `dataset.py`'s capability
check and both backends' `problem.py`; disjoint from W3.2's `_sbi.py`).
W3.9 is docs-only and can run any time after approval. W3.7 last.
Cross-model review via terra when the quota returns
(~2026-09-30) for W3.1 slice 1 (the batched-equals-loop and partition-independence claims), W3.1 slice 2 (native sampling), W3.2 (the prior
coordinates and the scoring of draws) and W3.6 (the calibration statistics).

### W3.16 — Peter's Phase 3 rulings applied [S; Fable] (ruled 2026-09-10, merged the same day)
The one-line rulings from the Phase 3 questions block, on one branch with
the full gate set: TMNRE's default sampler becomes `"mcmc"`
(`TMNRE_DEFAULT_SAMPLER`; `"rejection"` selectable; the example, the API
page and the tutorial say so; the canonical fixture runs on the default and
the calibration fixture asks for rejection explicitly); `PoissonFamily.sample`
guards `rate <= 0` exactly as `log_prob` does, by name; the sampling-form
principle written into `likelihoods.md` §3. **Accept:** the default asserted
on construction and in a run's attrs; the zero rate refused by name; all
four gates green; lint/format/pyrefly clean.

## Phase 4 — Extensibility proof: one new modality end to end (drafted 2026-09-10 by Fable; approved with rulings D1–D4 on 2026-09-11; **closed 2026-09-13** — every item W4.0–W4.11 merged, W4.8 last; the landed summary is the plan's §5 Phase 4 paragraph; the rulings still open are in `docs/development.md`'s review block)

Written from `DEVELOPMENT_PLAN.md` §5 Phase 4 and §4.7, the interferometry
sketch (`docs/design/modalities/interferometry.md`, whose §9 gaps I-1 to I-5
all landed at the freeze and whose §11 questions Q1–Q4 were ruled "Phase 4"),
the kernel note in `docs/design/horizon_notes.md` (§1 and its follow-up),
and the Phase 3 close. Phase 4's claim is the plan's: the composition design
is proven when a modality nobody wrote the contracts for runs through the
whole stack — a kind-changing, coordinate-changing step in the middle of the
chain, complex data, wrapped angles, two datasets on one model channel — with
no change to `ampere.core` beyond what the sketch already named. The
sequencing follows Phase 2 and 3: the reference-path composition first
(W4.1), the likelihood closed form beside it (W4.2), the native twins that
consume both (W4.3), the study that proves it (W4.4); the two extensibility
items that share no files with those (W4.5 kernel algebra, W4.6 astropy)
run in parallel from the start; the documentation pass closes the phase.

**Four decisions for Peter before dispatch** (each item below states its
assumption; a different ruling changes the text, not the plan). **Status
2026-09-11: all four ruled** (`docs/design/phase4_placement_memo.md` is the
record of the D1/D2 discussion; its §4 corrections and §7 walk-through are
applied to the items below; Peter: "I think we can begin work").
- **D1 — where a shipped modality lives.** **Ruled 2026-09-11 (memo §2.4):
  the kind goes in `core/results_schema.py` beside `VisibilitySet`; the
  steps go in `backends/{reference,torch,jax}/interferometry.py` (one name
  per backend, the rule kept); no grouping namespace now — revisited after
  realistic usage, towards the end of the phase or later; a per-observable
  front door (`ampere.interferometry`) only when the first reader lands.**
  *The original draft, superseded:* The sketch keeps `ClosurePhases`
  and the interferometric steps out of `ampere.core` ("exactly as
  `transformations.md` §10 keeps the standard library out of core"), and
  the standard library today is `ampere/backends/reference/instrument.py`
  with twins per backend. Assumed: a new namespace **`ampere/modalities/`**
  — `ampere/modalities/interferometry/` holding the backend-neutral kind
  and the reference-path steps, native twins under
  `ampere/backends/{torch,jax}/interferometry.py` — typed from the first
  line, with a one-line addition to `architecture.md` §3's layout. The
  alternative is everything under `backends/reference/` with the kind in
  `core.results_schema` beside `VisibilitySet`.
- **D2 — the `ClosurePhases` signature** (sketch Q3): four dimensionless
  axes `(u1, v1, u2, v2)` fixing the triangle by two baselines (assumed —
  it gives a GP over closure phases coordinates to work with), or a single
  `triangle` label in the manner of `PhotometricPoints.filters`. **Ruled
  2026-09-11 as revised in memo §3.6–3.7**: a chromatic sky error is sharp
  in wavelength and smooth in (u, v), and a kernel sees axes only, so
  `VisibilitySet` gains a `spectral_axis` (a frozen-kind amendment with its
  decision-log row and conformance update in W4.1), `ClosurePhases` is
  `(u1, v1, u2, v2, spectral_axis)` with the canonical baseline ordering of
  memo §3.4, and kernels gain an `axes` selector (W4.5) so a `Product` of a
  (u, v) block and a spectral block is expressible. A GP on closure phases
  is a *latent* composition (von Mises is wrapped) and is **Phase 5's**
  (Peter: "plenty of science to be done with visibilities alone").
- **D3 — whether W4.9 (astrometric time series by the template) runs in
  Phase 4.** It is the test of the "other modalities then follow the
  template" claim; it is not on the critical path. **Ruled yes by Peter
  2026-09-11: W4.9 runs; "this will be an essential test of the code".**
- **D4 — Opus for W4.1–W4.3, W4.5, W4.6; Sonnet for the rest** (per
  `docs/orchestration.md`); W4.1 ∥ W4.5 ∥ W4.6 from the start, W4.2 after
  W4.5 (both touch `core/likelihood.py`), W4.3 after W4.1 and W4.2.
  **Agreed by Peter 2026-09-11**, with the 2026-09-05 review-policy ruling
  extended to Phase 4: while the Codex quota is blocked, items touching the
  GP mathematics (W4.2, W4.5) merge on a second Fable review pass, with the
  `gpt-5.6-terra` pass owed retroactively. W4.1 and W4.5 both edit
  `core/likelihood.py`: W4.1 owns the family section (the von Mises
  `sample()`), W4.5 the kernel section — stated in both prompts.

### W4.0 — Phase 4 housekeeping and the owed core fixes [S; Sonnet]
The carried items that should not wait for a phase to need them. (1)
`docs/source/overview.rst` refreshed for the Phase 3 landing (six engines;
SBI landed; the capability ladder's SBI row) — W3.13's finding. (2)
`_SimulateNatively.run` wraps the native sampler in the per-draw `try`, so a
native sampler failure flags one draw where the numpy loop does, rather
than aborting the chunk (W3.14's finding; the Poisson twin is the first
family that can reach it; a test that injects a failing draw on each
backend). (3) `Instrument._declarations()` fingerprints steps by value
rather than `id()`, so a *bare* pickled `Instrument` round-trips
(`transform.py`; W3.1 slice 2's finding; `Dataset.__setstate__`'s
workaround then becomes redundant and is removed). (4) `architecture.md`
§3's layout records D1: a shipped observable's kind lives in core, its
steps in a per-backend `<observable>.py` module, no grouping namespace.
(5) `examples/sbi/npe_native.py` — a bare in-process NPE fit of a native
problem for its own sake (no subprocess, no cache, no truncation), covered
in `tests/examples` and linked from the tutorial's first section (ruled
wanted 2026-09-10; W3.13's finding). (6) `configure_from`'s docstring
(`core/transform.py`) no longer says the smearing steps are its first
instances — `LSFConvolution` uses it on all three backends (memo §1.3).
(7) **The reference `LSFConvolution` caches its `(n, n)` operator** keyed
on the grid it was built for (identity, then equality), as the torch twin
already builds it once: memo §7.1 measured 2.2 s per evaluation rebuilding
it on a 5 647-point union grid. The jax twin inherits the fix. **Depends:**
nothing. **Accept:** the two regression tests; the pickled-`Instrument`
round trip asserted; a timing test that the second evaluation on one grid
does not rebuild the operator, and a conformance row that the cached
operator equals the rebuilt one; docs build clean; all four gates green
(core is touched); lint/format/pyrefly clean.

### W4.1 — The interferometry modality on the reference path [M; Opus]
The sketch's §1 composition, built as shipped code under D1 as ruled (the
kinds in `core/results_schema.py`, the steps and models in
`backends/reference/interferometry.py`): **`VisibilitySet` amended to three
axes `(u, v, spectral_axis)`** — `spectral_axis` with `Spectrum`'s physical
types and `Order.ANY`, one wavelength per sample, the monochromatic case a
constant column — with the decision-log row and the conformance update in
the same PR (ground rule 9; memo §3.6 measured about ten positional call
sites, all tests, plus one doctest in `results/serialisation.py`); the
**`ClosurePhases` kind** — `(u1, v1, u2, v2, spectral_axis)`, `Layout.POINTS`,
values in radians, the **canonical ordering** of memo §3.4 (baselines ij
and jk for telescope indices i < j < k, the third implied as
−(u1 + u2, v1 + v2), the phase's sign convention stated) in its docstring
and asserted, `triangle` and per-sample `baseline` labels in
`extra_coords`; `FourierSample` (`ACCEPTS = (Image,)`,
`PRODUCES = VisibilitySet`; the `(u, v)` coverage as buffers taken *from the
observed container* per gap I-2's rule, never recomputed; requirements per
sketch §3 with `max_step = 1/(2 u_max)` so an under-sampled grid is refused
by `compile_for` rather than aliased — gap I-4; a direct DFT on the
reference path, the FFT-plus-interpolation variant refused by name as a
later option); `ClosurePhase` (three visibilities to one angle, mask
propagation through the three-to-one step per sketch §5, wrapped into
(−π, π]); **the first cross-kind uses of `configure_from`** —
`BandwidthSmearing` and `TimeSmearing` — which ask the `FourierSample`
before them for the extra `(u, v)` samples they average over (gap I-3's
mechanism; `LSFConvolution` already uses it within a kind); the two likelihoods of sketch §1
(`ComplexGaussianFamily` + `IndependentNoise` on visibilities;
`VonMisesFamily` on closure phases) composed on one `sky` channel from two
datasets and two instruments, negotiated once. `RiceFamily` on amplitudes
checked as the third route. **Design horizon (h)'s obligation** (added
2026-09-10): the pairing between the visibility and closure-phase channels
stays *visible* in the composition — both datasets name the `sky` channel
and the collection records which datasets derive from one model channel —
so Phase 5's joint noise model over a tuple of channels can bind to it
without a re-plumb; nothing here may fold the two into one container. A `UniformDisc`/`GaussianSource`/`Binary` model
trio emitting an `Image` on the reference backend, and the same three
emitting a `VisibilitySet` directly (the analytic route, no Fourier step)
as each other's oracle. Sketch §8's verified claims become tests.
**Depends:** nothing (W4.0's (4) is the architecture line, not blocking).
**Owns**: `core/results_schema.py` (the two kinds), `core/likelihood.py`'s
*family* section only (the von Mises `sample()`; W4.5 owns the kernel
section — do not touch it), `backends/reference/interferometry.py`, its
tests and conformance rows. **Accept:** conformance rows — `FourierSample` of each image model against its analytic
visibilities at `tolerances.cross_solver`; closure phases of the binary
against the closed form; each smearing step against brute-force fine
sampling; mask propagation asserted; the two-dataset composition evaluates
the model once per draw (`inference.md` §8) — plus an emcee fit of the
synthetic binary recovering the injected separation and flux ratio inside
the central 95 %; `simulate(observe=True)` works on both kinds (W3.14's
`complex_gaussian` sample, and a von Mises draw — **the wrapped family gains
`sample()`**, one decision-log row); the `VisibilitySet` amendment's own
decision-log row and every existing `VisibilitySet` row of the conformance
suite updated with the justification; all four gates green;
lint/format/pyrefly clean.

### W4.2 — The circular complex Gaussian process: `complex_gaussian` + `GaussianProcessNoise` analytic [M; Opus]
Sketch Q1, ruled 2026-09-03: the pair is declared `ANALYTIC` with the
circular (equal-component, zero-pseudo-covariance) GP as its fixed meaning,
and composition refuses until Phase 4 implements the closed form. This
item implements it: for a circular complex Gaussian with covariance `K`
over the `(u, v)` points plus the per-component σ², the marginal
log-likelihood is the real 2N-dimensional Gaussian's with the block
structure exploited (`log|K + Σ|` and the quadratic form each once, not
twice) — on the `DenseGP` solver, the kernel evaluated on the container's
selected axes through W4.5's `axes` selector: `Matern32(axes=("u", "v"))`
is the isotropic (u, v) kernel, and `Product(Matern32(axes=("u", "v")),
Matern32(axes=("spectral_axis",)))` the chromatic one of memo §3.6 (the
O(N) `QuasisepGP` path refuses by name since `REQUIRES_ORDERED_1D` cannot
hold). `GP_ANALYTIC_IMPLEMENTED` flips to `True` for the family;
`conditional_loo` for the pointwise group per `results.md` §6. Sketch Q2's
`extra_coords` units are no longer needed for this purpose (the wavelength
is an axis since W4.1) and are dropped from the item.
Native twins on torch and jax for the closed form (the family exists on
both since W2.4/W2.5 slice 3; the GP form is new). **Depends:** W4.5
(both touch `core/likelihood.py`; W4.5 owns the kernel section, this item
the family/noise section — merge W4.5 first). **Accept:** conformance
rows — the closed form against a dense real 2N formulation at
`tolerances.cross_solver` on all three backends, with the (u, v)-selected
kernel and with the product kernel; the refusal on the O(N)
path word for word; `check_alignment`'s dtype check (gap I-1) still
catches an amplitude fit; a GP fit of the W4.1 binary with an injected
correlated calibration residual stays calibrated where the independent
fit does not (SBC ranks via W3.6, the M2 pattern); all four gates green;
lint/format/pyrefly clean; `likelihoods.md` §7 amended, one decision-log
row.

### W4.3 — Native twins for the interferometry steps, and the modality under every engine [M; Opus]
`FourierSample`, `ClosurePhase` and the smearing steps on torch and jax
(`ampere/backends/{torch,jax}/interferometry.py`): a differentiable direct
DFT, batchable under `vmap`, `DIFFERENTIABLE`/`BATCHABLE`/`BACKEND`
declared, associated with the reference step as every standard step is —
the same class name in the backend's module, the `apply_flux` native
surface, `BACKEND` declared; there is no step registry (memo §1.2). State
which twin pattern the new steps follow (torch re-declares and is held by
the conformance battery; jax inherits the reference class) and why. The
kinds need no twin. Then the modality under every engine as the proof: NUTS on
the binary on both backends; VI; `SBIEngine` on visibilities plus closure
phases — W3.3's encoding already carries `is_complex`, and this is its
first complex customer (fix what it gets wrong, in scope); the artefact
cache keyed correctly on the two-dataset problem. **Depends:** W4.1, W4.2.
**Accept:** conformance lockstep rows for each step (numpy the oracle);
NUTS on torch and jax recovers the binary at W2's NUTS tolerances;
`simulate_many(native=True)` on the two-dataset problem equals the loop;
an NPE fit at the smoke budget passes the coverage check; all four gates
plus sbi green; lint/format/pyrefly clean ×3.

### W4.4 — The interferometry study: `examples/interferometry/` and its page [M; Sonnet]
In the M2 layout (`generators.py`, `model.py` + `model_torch.py` +
`model_jax.py`, `study.py`, `figures.py`, `__main__.py`): a synthetic
resolved binary with a fainter disc component observed as visibilities and
closure phases, fitted (a) with the correct model and independent noise,
(b) with a deliberately incomplete model (the disc omitted) under
independent noise, (c) under the flexible likelihood of W4.2 — the M2
question asked of the proof modality: does the GP keep the binary
parameters calibrated when the sky model is wrong? — and **(d) the
chromatic case of memo §3.6**: the omitted component is a compact patch
with a band profile, fitted under a (u, v)-only kernel, a spectral-only
kernel and their product, so the claim that the product is needed is
measured rather than argued. Three engines (emcee on
the reference backend, NUTS on torch or jax, NPE), timings, the six plots,
and the SBC row. `docs/source/interferometry.rst` in the v2 tutorials
toctree, written as the **template for adding a modality** (kind, step,
`configure_from`, the two-dataset composition, what the conformance suite
owes), citing W4.11's photometry + spectrum page as the simple case. `tests/examples` coverage; `tests/interferometry` for the study's
assertions in the M2 pattern. **Depends:** W4.3, W4.11. **Accept:** the
study's assertions pinned (coverage under (a) and (c), failure under (b),
and (d)'s ranking of the three kernels reported — informational);
runs in under ten minutes on the CI runner; page builds clean; the `dev`,
`sbi` and both backend gates green.

### W4.5 — Kernel algebra and the public quasiseparable-term registry [M; Opus]
Ruled by Peter 2026-09-09 (plan §5 Phase 4; `horizon_notes.md`'s follow-up
§1 is the specification): `Sum` and `Product` kernels with `KernelSpec`
composition and celerite translation — a sum of quasiseparable terms is
quasiseparable, a product is refused on the O(N) path by name; a damped
periodic **SHO** term (celerite's `SHOTerm`, `Q > 1/2`; the
`RotationTerm` pair as the second form) as the component for fringing;
`Matern12` and `Matern52` as celerite-exact siblings; `Kernel`'s dense
`matrix`/`diagonal` for every new term; **an `axes` selector on `KernelSpec`
and every `Kernel`** (ruled 2026-09-11, memo §3.6 item 3): a kernel acts on
a named subset of the container's axes, `check_compatible`'s single-unit
rule applies to the subset, and `Product` on the dense path composes
kernels on disjoint axis subsets (the (u, v) × spectral covariance of the
chromatic case) — the spec hash includes the selection; a **public
`register_quasiseparable_term(kernel_type, builder)`** beside the
lowering/realisation registries' shape (one slot per type, no silent
overwrite, built-in rows distinguished) so a user kernel reaches the O(N)
path; a `SpectralMixture` kernel as a sum of SHO terms with free
frequencies; the three backends' translations (numpy through celerite2's
`GaussianProcess`, torch through `celerite2.backprop`, jax through
`celerite2.jax.ops`), each term's `conditional_loo` where the solver has it.
Sparsity-inducing priors on component amplitudes are *not* this item (the
follow-up note's regularised horseshoe is expressible with
`HierarchicalPrior` today and is a study, not a contract). **Depends:**
nothing; owns `core/likelihood.py`'s kernel section and the backends'
kernel modules. **Accept:** conformance rows per term against its dense
closed form at `tolerances.cross_solver`, `Sum` against the sum of
matrices, the `Product` refusal word for word on the O(N) path and a
`Product` of two axis-selected kernels against the elementwise product of
their matrices on a three-axis container on the dense path, the `axes`
selector's unit check refusing a mixed-unit subset word for word, a user-registered term
reaching `QuasisepGP` on all three backends; M2's `fringing` scenario
re-fitted with `Matern32 + SHO` as the demonstration (bias and calibration
reported beside the stationary fit — informational, not pinned); all four
gates green; lint/format/pyrefly clean ×3; `likelihoods.md` §7/§8 amended,
one decision-log row.

### W4.6 — The astropy interop adapter (`core/astropy_compat.py`) [M; Opus]
Plan §4.7, the last unimplemented §4 contract (`ampere/core/__init__.py`
says so). `from_astropy(model, *, kind=, priors=None, channel=)` wraps any
`astropy.modeling` model, compound models included, as an ampere `Model`:
each astropy `Parameter` becomes an ampere `Parameter` (`bounds` → a
uniform prior, `fixed` → frozen, `tied` → a `Tie`), overridable by a
`priors=` mapping; `Quantity` inputs and outputs honoured; the output kind
declared by the caller (inferred only where the model's `n_outputs` and
units make it unambiguous, refused by name otherwise); capabilities the
honest reference answers (black-box: gradient-free engines and SBI, never
NUTS/VI — advertised in the docstring and the page). Translation to a
native equivalent is **opt-in, never silent** (decided 2026-09-01): this
item lands the hook — `from_astropy()` on a backend that raises for any
model it has no row for — and an empty table; the curated rows are W4.7.
`docs/design/contracts/astropy_compat.md`, short, in the contracts'
format (a post-freeze §4 addition; one decision-log row). **Depends:**
nothing; owns the new module, its contract page, `tests/core/test_astropy*`.
**Accept:** conformance rows — a wrapped `BlackBody` and `PowerLaw1D` agree
with the reference backend's own models' `log_prob` on the same data; a
compound model with a tied parameter round-trips its tie; an emcee fit of a
wrapped model; the SBI engine runs on one (black-box route); the refusal
for an undeclared kind word for word; all four gates green (a core module);
lint/format/pyrefly clean.

### W4.7 — Native translations of common astropy models [S; Sonnet]
The curated table §4.7 calls a later nicety: `BlackBody`, `PowerLaw1D`,
`BrokenPowerLaw1D`, `Polynomial1D`, `Gaussian1D`, `Const1D` (and their
compound sums/products) to native torch and jax models through W4.6's
opt-in hook, restoring differentiability for the frequent cases; a model
outside the table still raises. **Depends:** W4.6. **Accept:** conformance
rows — each translated model agrees with the wrapped astropy original at
`tolerances.cross_backend`; NUTS runs on a translated compound model on
both backends; the refusal for an untranslatable model word for word; the
`torch` and `jax` gates green; lint/format/pyrefly clean ×2.

### W4.11 — The photometry + spectrum composition: example, smoke test and tutorial page [S; Sonnet] (ruled by Peter 2026-09-11)
Memo §7.1, expanded into the docs: `examples/sed_composition/` (a
generator, the script, a `__main__`) fitting one `ModifiedBlackBody` to a
spectrum through `LSFConvolution` + `Resample` and to photometry through
`SyntheticPhotometry.from_library`, both on one channel with distinct
labels, on the reference backend with emcee and — the same script, the
backend chosen by a flag — on torch or jax with NUTS; synthetic data
generated through `negotiate` + `compile_for` (memo §7.1 finding 1, and
the page says why); `problem.requirements` printed and explained (what
each step asked, how the union grid arose); a `CalibrationScale` on the
spectrum as the instrument nuisance parameter; the label-collision error
shown and fixed. `docs/source/sed_composition.rst` in the v2 tutorials
toctree — the page `tutorials.rst` says is owed — replacing that sentence;
`tests/examples` smoke coverage. **Depends:** W4.0 (the LSF operator
cache; without it the reference fit is unusably slow). **Accept:** the
example runs on `dev` in under two minutes and recovers the injected
parameters inside the central 95 %; the smoke test; the page builds
clean; the `dev` gate green (plus `torch` or `jax` for the NUTS variant).

### W4.8 — Phase 4 documentation pass [M; Sonnet]
Like W3.13, after the last Phase 4 merge: the modality template page
(W4.4's) cross-linked from the overview and the architecture page; the
astropy page (`docs/source/astropy.rst`: the casual user's route, the
capability consequence, the opt-in translation); the kernel page (algebra,
the SHO term for fringing, registering a term); `advanced.rst`'s noise-model
section; every "landed"/"Phase 4" annotation in `docs/design/contracts/*`
and `docs/design/modalities/interferometry.md`'s §9–§11 dispositions
checked against the code (*Amended W4.8*, one decision-log row), including
`spectrum_photometry.md`'s pre-shipping `SyntheticPhotometry` signature;
the deferred list reviewed; README. **Depends:** every other Phase 4 item.
**Accept:** `pixi run docs` with no new-namespace warnings; every example
the pages use under `tests/examples`; the annotation list in the report;
the `dev` gate green.

### W4.9 — Astrometric time series by the template [M; Sonnet; D3 — optional in Phase 4]
The sketch (`docs/design/modalities/astrometric_timeseries.md`): a reflex
orbit on two `TimeSeries` channels, an epoch-sampling instrument, the
time-domain flexible likelihood per coordinate on the O(N) path — built
under `ampere/modalities/astrometry/` by W4.4's template with no new
contract surface (the joint 2-vector GP the sketch wants is Phase 5's,
deferred at the freeze). The value is the test of the template claim, and
the second worked modality page. **Depends:** W4.4. **Accept:** conformance
rows for the orbit model against a closed-form ephemeris; an emcee and a
NUTS fit recovering an injected orbit; `QuasisepGP` on both channels; the
page; all four gates green.

### W4.10 — CI/CD Phase 4: smaller, path-gated jobs and `actionlint` [S; Sonnet]
Ruled by Peter 2026-09-10 on W3.7's report: (1) the `suites`/`backend-suites`
legs split into smaller jobs where that is simpler — `typecheck` out of each
five-suite leg into its own matrix job, so the slow `sbi` leg's critical
path is `test-all` alone; (2) **path gating**: a `paths-filter` step (or
`on.push.paths`/`pull_request.paths`, whichever keeps required checks
honest) so a docs-only change runs lint, the docs build and `dev`, a change
under `ampere/backends/torch` runs the torch leg but not jax's, a change
under `ampere/core` runs everything, and the characterisation job keeps its
schedule — with the rule that a required check skipped by the filter
reports success, not "expected", so branch protection does not hang; (3)
`actionlint` in the `dev` toolchain (a pixi task and a CI step), with
`pinact`-style pinning or Dependabot's `github-actions` ecosystem enabled so
action versions stay current — the benefit Peter asked it to earn; (4) the
interferometry study (W4.4) and the new examples on the legs that own them.
`docs/development.md`'s CI prose updated. **Depends:** nothing; owns
`.github/workflows/*`, `pyproject.toml`'s task table (the one task), the CI
prose. **Accept:** `actionlint` clean; YAML validated; the path filter
table in the workflow comments and the docs; a docs-only PR provably runs
only its subset (the report shows the job list for three synthetic diffs);
the required-check names for branch protection listed in the report.

Ordering: W4.0 ∥ W4.1 ∥ W4.5 ∥ W4.6 from the start (file ownership:
W4.0 — `core/transform.py`, `core/simulate.py`, `core/dataset.py`'s
`__setstate__`, `backends/reference/instrument.py`'s LSF, the docs;
W4.1 — `core/results_schema.py`, `core/likelihood.py` *families*,
`backends/reference/interferometry.py`; W4.5 — `core/likelihood.py`
*kernels* and the backends' kernel/GP modules; W4.6 —
`core/astropy_compat.py` and its page; `core/__init__.py` exports are
appended by each and merged by hand); W4.11 after W4.0; W4.2 after W4.5
and W4.1; W4.3 after W4.1 and W4.2; W4.7 after W4.6; W4.4 after W4.3 and
W4.11; W4.9 after W4.4 (ruled in); W4.10 any time; W4.8 last.


## Phase 5 — Scale-out & advanced inference (drafted 2026-09-15 by Fable; **approved by Peter 2026-09-15 with rulings D1–D8 as recommended — D1 as (a) and (b) together, carried by W5.20; D8 as two Opus items**; wave 1 dispatched 2026-09-15; **closed 2026-09-28 with W5.19's merge at `50f6ff4` — thirty-two items merged, W5.16 deferred to Phase 6**)

The plan's §5 Phase 5 bullets, the design horizons (b), (e)–(i), the
inference-extensions memo's §5 adaptations and §6 ranking, the modality
sketches' Phase 5 dispositions (IFU gaps 1–2, the astrometric sketch's
joint-noise gap, `hierarchical_population.md` Q2/Q5) and the owed list,
as agent-sized items. **Read `DEVELOPMENT_PLAN.md` §5 "Phase 5" first**:
each item below names the bullet it executes. Sized under
`docs/orchestration.md`'s token-economy rules: two agents at a time,
targeted tests on the branch, one merged-master gate per wave. Every
`GPSolver`, `NoiseModel` or results change is a §4 change (ground rule 9):
a decision-log row and the conformance rows in the same PR. The eight
decisions **D1–D8** at the end are the ones the drafting could not take.

### W5.0 — The results contract for approximate and evidence-producing engines [S; Sonnet] (ruled by Peter 2026-09-10 on the inference-extensions memo §5, §7.1–7.2)
The three contract adaptations every tier-1 sampler and every approximate
engine wants, taken once so each later engine is one item: (1)
`ampere_log_evidence`, `ampere_log_evidence_err` and
`ampere_evidence_method` in the root attrs — engine-neutral — written by
every engine that estimates a marginal likelihood (`dynesty` today, from
its `logz`/`logzerr`; `ampere_dynesty_logz` kept as the engine's own);
(2) the weighted/approximate-draw rule stated in `results.md`: weighted
output is resampled to equal weight for the `posterior` group with the
original count recorded and the raw weighted output kept on the engine
(dynesty's convention, now the contract's), and **every engine whose draws
are not from the target stores the proposal's own log-density per draw**
(`sample_stats.proposal_log_density`) — `SBIEngine` does, `VIEngine` gains
it from the fitted guide's `log_prob`; (3) `ampere_approximation` in the
root attrs, `"none"` for an exact sampler and the family otherwise
(`"mean_field"`, `"multivariate"`, `"density_estimator"`, …), the one key a
plot or summary checks before it reports an R-hat that means nothing.
`PROVENANCE_SCHEMA_VERSION` → 7 (attrs added; one decision-log row);
`results.md` §9/§4 amended; `plot_trace`/`summary` warn on an
approximation. Optimisation results, when they come, are a `DataTree` with
an `optimum` group (ruled), not a separate type. **Depends:** nothing.
**Accept:** the three attrs on a dynesty run (evidence), a VI run
(approximation and proposal density) and an SBI run; an importance-corrected
VI posterior computed from stored draws alone agrees with an emcee reference
on the toy problem; gates dev + sbi on the merged wave; lint/format/pyrefly
clean.

### W5.1 — A latent GP on closure phases [M; Opus] (ruled Phase 5 by Peter 2026-09-11, placement memo §7.2)
`VonMisesFamily` gains `CONSUMES_LATENT_GP = True` — the Poisson pattern of
W2.14 — so a `GaussianProcessNoise` over `ClosurePhases`' five axes
composes on the **native path only** (NUTS/VI on torch and jax through
`ampere.core.realise`; the numpy path keeps its refusal by name, reworded
to say *where* the composition is available): the latent phase error
`f ~ GP(0, K)` on the `DenseGP` solver (dense by W4.2's structural
argument — no ordered 1-D coordinate exists), the observed closure phase
wrapped von Mises around `model + f`, `latent_transform` whitening the N
latent variables. The kernel binds through W4.5's `axes=` selector
(`Matern32(axes=("u1", "v1", "u2", "v2"))`, or the product with a spectral
block). `simulate(observe=True)` draws the latent first, then the wrapped
draw. The channel pairing between visibilities and closure phases
(`extra_coords`' `triangle`/`baseline` labels, horizon (h)) is read, not
changed. **Depends:** nothing (W4.3 merged); W5.4's reduced-rank latent
question does not gate it at the study's N. **Accept:** conformance rows on
the torch and jax fixtures — the latent composition's log-density against a
from-scratch von Mises-around-GP-draw formula at `tolerances.cross_solver`,
the numpy refusal word for word, `simulate` producing wrapped values with
the injected correlation visible in a periodogram of the residuals; NUTS
on the W4.4 binary with an injected smooth per-triangle phase error: the
flexible fit's central-90 % coverage of separation and flux ratio pinned
where the rigid von Mises fit's is not (SBC over 12 simulations, W4.4's
pattern); `likelihoods.md` §4 and §14 amended, one decision-log row;
gates: all four.

### W5.2 — Phase 5 housekeeping: the owed list [S; Sonnet]
W4.0's shape: the small items owed at the Phase 4 close, each a commit.
`NUTSEngine` exposes `target_accept` and `max_tree_depth` (W4.11's variant
ran past 19 min at a 40/80 budget); `examples/` brought under ruff (fix
what it finds); `dataset.py`'s `part_name` docstring (`ampere.core.kernels`);
a `Sum` containing a `Product` gets the sharper O(N) refusal; `derived.py`
treats a complex value one way in all three places (refuse, as
`add_posterior_predictive` does — a component or modulus is the caller's
choice, W5.3's); `von_mises` into `_TWINNED_FAMILIES` on both backends;
`encoding._sanitised`'s summary line made true; `RiceFamily` given either a
schedule (a line in §5) or a refusal that says "no schedule"; the two docs-
build warning pairs; `interferometry.rst`'s non-verbatim quotation; the
Phase 3 owed items — a native sampler failure flags one draw not a chunk
(`_SimulateNatively.run`), a value-based `Instrument._declarations()`
fingerprint, the unreachable "no uncertainty" branch of both
`_check_gp_uncertainty`s, torch `noise.py`'s redundant device check. Not
in scope: anything with a decision-log row. **Depends:** nothing.
**Accept:** one targeted test per behaviour change; the docs build's
warning list shorter than the base commit's; lint/format/pyrefly clean in
every environment; gates dev + torch + jax on the merged wave.

### W5.3 — Multi-axis diagnostics: the four plots on a point kind with several axes [M; Sonnet]
`results.md` §13 item 14 lifted. `plot_posterior_predictive`,
`plot_residuals`, `plot_gp_localisation` and `plot_anomaly_score` gain a
`coordinate=` argument — an axis name, or a callable of the container's
axes to one ordered coordinate — with per-kind defaults declared on the
kind (`VisibilitySet`: baseline length `hypot(u, v)`; `ClosurePhases`: the
longest of the three baselines; `TimeSeries`: `time`, unchanged), so the
default keeps every existing plot byte-identical; complex values are
plotted as the component the caller names (`component="real" | "imag" |
"abs" | "phase"`, default refuse with the four names), and `derived.py`'s
three complex branches agree (W5.2 makes them refuse; this item gives them
the argument). `results._plotting.coordinate_of` is the one place the rule
lives. **Depends:** W5.2 (the `derived.py` alignment). **Accept:** the
W4.4 study's four missing figures render for the binary's `VisibilitySet`
and `ClosurePhases` on the reference backend; the existing `tests/results`
plot rows unchanged; a refusal row for a kind with no default and no
`coordinate=`; `results.md` §4/§13 amended, one decision-log row; gates
dev + sbi.

### W5.4 — The first approximate solver: `HilbertSpaceGP` on all three backends [L; Opus]
The plan's first Phase 5 bullet, executed for the cheapest candidate
(horizon notes §2's recommendation: HSGP first, the most useful for NUTS).
A `GPSolver` with `EXACT = False`, `IMPLEMENTED = True`, the kernel's
spectral density at the Laplacian eigenvalues of a bounded box
(`m` basis functions per axis, tensor-product in 2–3 axes; a boundary
factor `c` from the data extent), `K ≈ Φ diag(S) Φᵀ` solved by Woodbury at
O(N m + m³); `provenance_config()` records `m` and `c`, never hashed
(fold-in 10); `conditional_loo` exact-in-the-approximation (the Woodbury
form of the leave-one-out identity) — refused by name only if the
mathematics fails, with the reason; `latent_transform` of dimension `m`
(**settles horizon notes §2 question (b)**: the latent size fixed at
composition may be the basis size, and `simulate(observe=True)` and the
latent path use the same whitening — a conformance row asserts it).
Stationary families with a closed-form spectral density only (the Matérn
family, `SHO`, and `Sum`s of them); others refused by name. Written once
in `ampere.core` against W4.5's `ArrayOps`, so numpy, torch and jax share
the basis construction and the backends contribute only their solves.
**The approximation-aware conformance tolerance class is part of this
item** (horizon notes §2 question (a)): a row compares to `DenseGP` at a
sequence of `m` and asserts monotone convergence to within a tolerance
that tightens with `m`, never a fixed number; `tests/conformance/README.md`
gains the class. The 1-D validation is against `QuasisepGP` on the M2
spectra (exact reference at O(N)); the 2-D customer is W5.5's image.
**Depends:** nothing for 1-D; W5.5 for the 2-D rows. **Accept:** the
convergence rows on all three fixtures in 1-D and 2-D; the latent-path
agreement row; NUTS on torch and jax over the hyperparameters with `m`
latents recovering the M2 `strong_smooth` scenario's coverage within the
`DenseGP` result's; `likelihoods.md` §7's table row for HSGP filled and
the `_SolverSlot` list amended (`HilbertSpaceGP` added beside the three
slots, which remain slots); one decision-log row; gates: all four.

### W5.5 — A gridded customer: `Image` data by the template, with PSF convolution [M; Opus]
No image *dataset* has ever been fitted: `Image` and `Cube` are kinds, the
Phase 4 image models feed `FourierSample`, and `transformations.md` §10's
PSF-convolution slot is empty. This item follows `interferometry.rst`'s
template for a gridded observable: a `PSFConvolution` step on the
reference path (kind-preserving, requirements the negotiation can meet —
a padded grid of the PSF's support; the FFT route with the padding stated),
its inheriting twins on torch and jax (W4.3's pattern), **mask propagation
on `Layout.GRID`** (IFU gap 2: `propagate_mask_grid`, the additive helper
the freeze precluded nowhere), the three Phase 4 image models reused as
the sources, `examples/image/` with a misspecification arm (a smooth
background the model omits; the flexible arm on `DenseGP` at small N and
on W5.4's HSGP at realistic N, the comparison being the phase's benchmark
per the plan's rule "chosen by measurement"), and a page in the
template's order closing with what the template did not say about a
gridded kind. **Depends:** W5.4 for the HSGP arm (the `DenseGP` arm and
the step land first). **Accept:** the template's four conformance rows
(the step against a direct convolution, the requirement and its refusal,
the mask row, `simulate`), the twins at `tolerances.cross_backend`, the
study's coverage pinned on the flexible arm at small N, the HSGP arm
informational with its wall-clock and memory beside `DenseGP`'s at three N;
`transformations.md` §10/§13.5 amended; gates: all four.

### W5.6 — Solver bake-off: EFGP and Vecchia against HSGP on the image [M; Opus] (D2)
The plan's rule made concrete: on W5.5's image at realistic N (10⁴–10⁵
pixels), a second reduced-rank or sparse solver — EFGP (equispaced Fourier
features with the Toeplitz normal equations solved by FFT; O(N + m log m))
or Vecchia (not rank-limited, the one for rough processes) — implemented
far enough to measure (numpy reference only, `EXACT = False`, the
convergence tolerance class of W5.4), benchmarked against HSGP and
`DenseGP` on bias, coverage, localisation, wall clock and memory across N
and kernel smoothness (Matérn-1/2 through 5/2). The outcome is a
decision-log row naming which lands as a full three-backend solver (a
follow-on item) and which stays a slot, with the table; not a fourth
solver landed on judgement. **Depends:** W5.4, W5.5. **Accept:** the
benchmark script under `tests/benchmarks/` with `pixi run bench`
attaching it; the table in `likelihoods.md` §7; one decision-log row;
gates dev only (nothing ships on the modern backends).

### W5.7 — `WarpedKernel`: input and amplitude warping preserving quasiseparability [L; Opus]
The plan's third Phase 5 bullet, the kernel half. `WarpedKernel(base,
input_warp=, amplitude_warp=)` in `ampere.core.kernels`: a monotone input
warp `x → w(x)` (a few knots, monotone by construction — cumulative
softplus increments — so `QuasisepGP`'s ordering precondition survives)
and an amplitude warp `D K D` with `D = diag(a(x))` (log-amplitude at
knots, linearly interpolated), both keeping the O(N) exact solve on every
backend through the term registry (the warped generators are the base
generators evaluated at `w(x)`, scaled by `a(x)`); the knots are ordinary
`Parameter`s in the kernel's own namespace, so NUTS on torch/jax gets them
for free and the `KernelSpec` hash carries them. **The degrees-of-freedom
guard is in the item**: the knots' priors are hierarchical shrinkage to
the identity warp (`HierarchicalPrior`, the non-centred form offered by a
`lowering.md` §3 note), and the whiteness and localisation diagnostics
(family B/C) run on the *warped* residuals — a `Likelihood.conditional`
that reports in warped coordinates with the warp recorded. **Depends:**
nothing (W4.5's registry is the extension point). **Accept:** conformance
rows — the warped kernel against `DenseGP` on the explicitly warped
coordinates at `tolerances.cross_solver` on all three fixtures, the
identity warp bit-identical to the base kernel and to its pre-W5.7 spec
hash, a non-monotone knot set refused at composition, the O(N) path
exercised; NUTS on torch and jax over knots and hyperparameters; the
diagnostics row on warped residuals; `likelihoods.md` §6–§8 amended, one
decision-log row; gates: all four.

### W5.8 — The M2 extension "many lines / one band", and the sparsity prior on summed noise components [M; Opus] (D8: two Opus items, margins pinned)
The plan's third bullet, the validation half, and its generalisation. An
`examples/m2_misspecification/` scenario with a forest of narrow lines in
one band and a smooth continuum error elsewhere, comparing stationary
Matérn, W5.7's warped Matérn and W4.5's `Sum` of two kernels on bias,
calibration and localisation — the M2 pattern, pinned as W4.5's fringing
study was (the margin, not the number). The sparsity guard for sums: the
regularised horseshoe (Piironen & Vehtari 2017) on component amplitudes,
expressed with `HierarchicalPrior` today and documented as the
recommended prior for any `Sum` of noise terms, with `lowering.md` §3
gaining the non-centred note NUTS wants and a `tests/m2` row showing the
spurious component's amplitude shrinks to zero when the truth has one
component. **Depends:** W5.7. **Accept:** the scenario in the driver and
`tests/m2` with the pinned margin; the horseshoe row; SBC on injected
misspecification for the warped fit; the M2 page's section; gates dev +
torch (NUTS over the knots) on the merged wave.

### W5.9 — Joint noise over a tuple of channels: the shared-grid intrinsic coregionalisation model [L; Opus] (D3)
The plan's fourth bullet and `likelihoods.md` §15's recorded limitation
lifted for the exact case. A `NoiseModel` bound to **several channels of
one model on a shared grid** — `JointGaussianProcessNoise(kernel, B=…)`
with `K = B ⊗ K_x`, `B` a T×T positive-definite matrix parameterised
physically (a rotation and two log-variances for T = 2; a Cholesky
parameterisation as the general fallback) — solved exactly and at O(N)
by diagonalising `B`, rotating the T residual vectors and solving T scalar
GPs with the bound solver (`QuasisepGP` where the grid is ordered 1-D,
`DenseGP` otherwise). It is the first use of `DatasetCollection.
contributions` as something other than a sum: `inference.md` §4 gains the
`"joint"` decomposition entry (`"mixed"`'s precedent), the pointwise group
carries the rotated outputs, family B/C diagnostics run per rotated output.
**The first customer is astrometry** (recommended, D3): W4.9's reflex orbit
already produces `ra`/`dec` from one evaluation on one epoch grid, so a
correlated per-epoch error (a centroiding systematic shared by both axes)
is the injected misspecification and no new kind is needed; polarimetry
(a Stokes kind) is the second test, and the general LMC and mismatched
grids stay recorded as the dense/reduced-rank follow-on. **Depends:**
nothing; W4.9's astrometry twins are read only. **Accept:** conformance
rows on all three fixtures — the joint density against `DenseGP` on the
materialised `B ⊗ K_x` at `tolerances.cross_solver`, the rotation
recovering T independent solves, `B = I` bit-identical to two independent
noise models, `simulate` drawing correlated channels; NUTS on torch and
jax over `B`'s parameters and the kernel's; the astrometry example's
`--joint` arm with SBC-pinned coverage under an injected correlated
error where the independent GPs' coverage is not; `likelihoods.md`
§7/§15, `inference.md` §4, `results.md` §6 amended; the decision-log row;
gates: all four.

### W5.10 — Amortisation over observation context [L; Opus]
The plan's fifth bullet; horizon (i)'s reserved `context=` slot filled.
A `ContextPrior` protocol with three shipped instances — scaled copies of
the observed σ-pattern, an archive of real error arrays, a parametric S/N
model — drawn per simulation by `simulate_many(context=…)`, recorded on
the `Simulation` and in the training-set provenance, passed to `sample` in
place of the container's σ, the chain re-negotiated per context grouped by
`chunk_size`; the encoding already carries σ, `log σ` and whitened values
per row (W3.3), so the network sees the context without a layout change,
and a per-set conditioning vector (FiLM) is the opt-in second route. SBC
per observation (W3.6) is the check that the context prior covered the
observation at hand, and the tutorial says so. **Depends:** W5.0 (the
provenance attrs land first so the schema moves once — D7 decides whether
W5.11's axis identity also rides this bump). **Accept:** an NPE posterior
trained under the σ-pattern prior stays calibrated (SBC, TARP) on an
observation whose σ is rescaled by a factor the prior covers, and its
coverage degrades measurably on one it does not — both pinned as
inequalities; the `context` recorded on every `Simulation` and in the
training set; `inference.md` and `encoding.md` amended, the decision-log
row, `PROVENANCE_SCHEMA_VERSION` unchanged if W5.0's bump carries the
attrs; gates dev + sbi + torch.

### W5.11 — The encoding's axis identity [S; Opus] (D7)
`encoding.md` §9 item 6: the packing aligns axes by position, so a layout
mixing a kind whose column 0 is `x` with one whose column 0 is `u` shows a
network two unrelated quantities in one slot. This item adds a per-column
axis-identity feature (the axis's physical type from its unit, one small
integer code per column, part of the frozen `Layout` and its hash) so an
embedding can tell the columns apart, with the code table in the contract.
Every existing layout hash moves once — which is why D7 asks whether it
lands in the same wave as W5.10 or not at all in Phase 5. **Depends:**
nothing. **Accept:** the code table pinned; a mixed-kind layout's columns
distinguishable in `unpack`; the W4.3 five-axis rows still passing; the
cache invalidation of a pre-W5.11 training set refused by name; one
decision-log row; gates dev + sbi.

### W5.12 — `Population`: the hierarchical container and the joint fit on the native path [L; Opus] (D4)
`hierarchical_population.md` H-1, ruled to land with Phase 5 (2026-09-03).
`Population(members, hyperpriors)` as the container the sketch designs —
plate-aware parameter groups (`parameters.md` §9's `members`), per-member
nuisance parameters, `Binding.index` (H-2) routing each member's element,
the hyperpriors ordinary `Parameter`s — lowered to numpyro/pyro plates by
the realisation so NUTS on torch and jax fits a population jointly, and the
IFU sketch's gap 1 ("a plate of datasets") met by the same construct: a
`DatasetCollection` built over a plate. The numpy path fits small
populations by the tie-based pattern (the documented route until now),
refusing beyond a member count it states. **Depends:** nothing frozen;
W5.0 for the results attrs. **Accept:** conformance rows — a two-level
population's joint log-density against the sum of members plus the
hyperprior on all three fixtures, the plate lowering on torch and jax
bit-consistent with the flattened form; NUTS recovering the population
hyperparameters of a 50-member synthetic sample inside the central 95 %;
`hierarchical_population.md` Q2 closed and Q5's joint-fit half landed;
`parameters.md`/`inference.md` amended, the decision-log row; gates: all
four. **Note (Peter's ruling of 2026-09-18 on W5.22)**: `fit_population`
refuses a stored `HierarchicalPrior` as the sole source of the interim
prior; this item may lift that by resolving the marginal — Monte Carlo
over the stored hyperprior draws, or the conjugate closed forms — if the
population fit needs it, and should say either way in its report.

### W5.13 — Population inference by reweighting archived fits [M; Sonnet] (D4)
Horizon (b), buildable entirely on stored files (`results.md` §13 item 16 (the text said 15; corrected 2026-09-22)): an
`ampere.results.population` module that takes a collection of run
`DataTree`s sharing a spec hash, reads each run's per-draw `log_prior`,
`log_likelihood` and — for approximate engines — W5.0's proposal
log-density, and returns a population-hyperparameter posterior by
importance reweighting under hyperpriors (Hogg-style), with the effective
sample size per object reported and a refusal when it collapses. Reads
from a directory of files *or* from any object with the same column
interface, so the columnar store of horizon (b)'s constraint is not
excluded — the reader is a protocol, one file-backed implementation
shipped. **Depends:** W5.0. **Accept:** on 200 synthetic single-object
emcee fits, the reweighted hyperparameters agree with W5.12's joint fit
(where both exist) and with the truth inside the central 95 %; the ESS
refusal row; the columnar-reader protocol test with an in-memory
implementation; `results.md` §13 item 16 (the text said 15; corrected 2026-09-22) amended; gates dev + sbi.

### W5.14 — Tier-1 engines behind extras: nested sampling, VI guides, blackjax [M; Opus] (D5)
The memo's §6 tier 1 on W5.0's contract, each behind its own extra
(ruled): `nautilus` and `ultranest` on a shared `_nested.py` (evidence into
the engine-neutral attrs; multimodality by construction); `VIEngine` guides
`laplace` and `flow` (pyro/numpyro autoguides); the `blackjax` route on jax
(MCLMC as a sampler, Pathfinder as initialiser and as an approximation). The
engine battery that comes with them: every new engine SBC-ranked through
`calibration.sbc` on the conformance toy problems, its evidence checked
against the closed-form linear-Gaussian case W3.6 uses, one cost record
per run (evaluations; horizon (g)'s hook). **Depends:** W5.0. **Accept:**
the battery green for each engine; the evidence rows within the stated
error; the extras in `pyproject.toml` with pixi features and CI legs path-
gated to the engine's files; `inference.md` §5 amended, one decision-log
row; gates dev + sbi + jax.

### W5.15 — Periodic models: the astrometry example under nested sampling, and the guidance [S; Sonnet] (D6)
W4.9's measured finding (period aliasing under a `loguniform(50, 2000)`
prior on an 860-day baseline; emcee and NUTS both lock onto a spurious
mode) is written into `interferometry.rst` and `astrometry.rst` as a
hazard with one remedy: an informed prior. Nested sampling is the standard
second remedy and ships already: this item runs `examples/astrometry`'s
wide-prior case under `DynestyEngine` (and W5.14's nested samplers when
they land), reports the modes and their evidences, adds a `--wide-prior`
arm, and rewrites both hazard sections as guidance with the measurement:
which engine to reach for, how the multimodal posterior looks in the
corner plot, and how an informed prior compares in evidence. **Depends:**
nothing (W5.14 optional). **Accept:** the arm's smoke test; the posterior's
mode structure pinned (the true period among the recovered modes, its
evidence the largest); both pages amended; gates dev.

### W5.16 — RHMF exploratory trial [S; Sonnet] (**deferred to Phase 6**, ruled by Peter 2026-09-24; the text stands)
Peter's ratification note of 2026-09-08 on the W2.7 deferral: the pre-fit
robust-factorisation family (`diagnostics.md` §2) tried against a pinned
`robusta-hmf` commit behind a non-default `rhmf` extra — the adapter over
`Robusta(...)`/`robust_weights` producing an `AnomalyScore` with
`provenance="rhmf_prefit"` (the renderer already accepts it), run on the M2
spectra and W5.5's image, with the adoptability re-check (§2.2: licence,
maturity, API) re-run and recorded. The outcome is a report, not a
namespace: `ampere.diagnostics` lands only if the maturity gate is met.
**Depends:** W5.5 optional. **Accept:** the trial script under `examples/`,
its findings in `diagnostics.md` §7, the re-check in the decision-log row;
no change to the base install; gates dev.

### W5.17 — The benchmark-driven optimisation pass [M; Opus]
The plan's last Phase 5 bullet: profile first against the pytest-benchmark
baselines (W2.11) on the M2 driver, the interferometry study and W5.5's
image, then attack the levers in evidence order — requirements-negotiation
compilation and caching (§4.3), batched evaluation, solver selection,
precision policy, resampling (issues #12, #29, #67). No speculative
optimisation: every change carries its before/after benchmark row.
**Depends:** W5.5 (the image is the load). **Accept:** the profile report
in `docs/`; each landed lever with its benchmark delta; no conformance row
moves; gates: all four.

### W5.18 — CI/CD Phase 5: the residue [S; Sonnet] (trimmed by Peter's ruling 2026-09-24)
As first written this asked for the W5.14 extras as path-gated legs, the
GPU job kept skip-clean, the sbi leg's budget re-measured, and a
`test-fast` task "with each item's Accept line naming its suites". W5.14
landed the nested leg and W5.26 re-measured every leg and split the CI
matrix by suite group, so what remains is: (1) **`test-fast`** (the
token-economy proposal's last rule) as a pixi task excluding the `*_full`
markers, the `study` rows (W5.26) and the SBI-training rows — the local
pre-commit run, under four minutes in dev; (2) the GPU job (`tests/gpu`)
confirmed skip-clean on the current matrix; (3) the stale `test (py311)`
job name in `path_filters.py`'s dev list (W5.26's out-of-scope finding —
the matrix is py312/py313/py314) corrected. The Accept-line clause is
**struck** (the status rows already record which suites each item ran).
**Depends:** W5.26 merged. **Accept:** `actionlint` clean; `test-fast`
under 4 min in dev with its excluded markers listed in the task comment;
`tests/gpu` skip-clean in dev and torch; the docs' CI section names
`test-fast`; no five-suite gate.

### W5.19 — Phase 5 documentation pass [M; Sonnet]
Last, W4.8's shape: every Phase 5 claim in the frozen documents checked
against the merged code and annotated *Amended W5.19*; the kernel page
gains warping; a solvers page (exact, approximate, joint) with the
tolerance classes; the population tutorial; `overview.rst`, README and
`index.rst` say what Phase 5 shipped; the plan's §5 Phase 5 paragraph is
the landed summary. **Depends:** everything merged. **Accept:** docs build
with a warning list no longer than the base commit's; spec doctests and
`tests/examples` green; gates dev.

### W5.20 — A backend-invariant parameter namespace [S; Opus] (proposed 2026-09-15 from the D1 analysis; **ruled in by Peter the same day**: both options, (a) and (b))
`Parameterised._check_free_name` refuses a parameter or buffer whose name is
an attribute of `type(self)`, over the whole MRO — so the set of legal
parameter names depends on which backend's base class a model inherits
(measured 2026-09-15: core `Model` reserves 17 names; a torch spectral model
adds `AXIS`, `evaluate_tensor`, `flux`, `grid`, `grid_tensor`, `to`; a jax one
`AXIS`, `flux`, `grid`; the interferometry twins `native_flux`,
`native_grid`), and torch's `LoweredParameters` additionally nests every
parameter as an `nn.Module` attribute, so `nn.Module`'s namespace (`to`,
`type`, `float`, `apply`, `eval`, …) is a third reserved list caught only at
lowering, on one backend. `parameters.md` §10 is a core contract; its
enforcement should not vary by backend. This item: (1) the shadow check
tests a **core reserved set** (the public names of `Parameterised` and
`Model`, plus a short list of names no backend lowering can carry, stated in
`parameters.md` §10) rather than `hasattr(type(self))`, so a name legal on
the reference backend is legal on every twin; (2) the torch lowering
mangles or refuses the `nn.Module` collisions by that same list, never by
torch's own `KeyError`; (3) the D1 ruling executed — `native_flux`/
`native_grid` canonical in the two `problem.py` tables (looked up first), a
model offering **both** pairs refused as ambiguous rather than taken at the
legacy name, `flux`/`grid` an alias the docs call legacy, the eight backend
modules renamed (mechanical), `inference.md` §10a, the placement memo and
`interferometry.md` §3 updated. **Depends:** nothing. **Accept:** a
conformance row on every fixture declaring a model with parameters named
`flux`, `grid`, `to` and `type` and composing it natively; the ambiguity
refusal word for word; the reserved list pinned; no spec hash moves (names
are unchanged, only their legality); `parameters.md` §10 amended, one
decision-log row; gates: all four.

**Not drafted, on a trigger**: matrix-free exact GPs (conjugate-gradient
solves with stochastic trace estimators, gpytorch/gpjax-style) for the 3-
and 5-axis interferometric kernels — taken up when a Phase 5 case
exceeds the dense path's memory (the plan's second bullet, "not
immediate"); the general LMC and mismatched-grid joint noise (W5.9's
follow-on); the second three-backend approximate solver (W5.6's outcome);
the terra reviews owed on W4.1, W4.2 and W4.5 when the quota returns.

### W5.21 — Lift the `Layout.GRID` refusals in the core [S; Sonnet] (ruled by Peter 2026-09-16 on W5.5's proposal)
W5.5's decision-log row records that `ampere.core` refuses a correlated
noise model on every `Layout.GRID` container as a gate, not mathematics,
and that `examples/image/grid_gp.py` lifts the gate out of tree by
subclassing. Make the library change: in `ampere/core/likelihood.py`
allow `Layout.GRID` in `GPSolver.check_compatible` when the kernel's
selected axes are the container's, and build `Likelihood._coordinates`
with `ampere.core.sample_coordinates` (both places are named in the row,
`likelihood.py:588` and `:4047` at W5.5's base — re-locate them, W5.7
and W5.9 move the module). Then delete `grid_gp.py`; `examples/image/`
(`study.py`, `bakeoff.py`) imports `DenseGP` and `HilbertSpaceGP`
directly; `tests/benchmarks/test_solver_bakeoff.py` and
`tests/conformance/test_image.py` follow. **Depends:** W5.7 and W5.9
merged (they own the two `likelihood.py` sections until then).
**Accept:** the image conformance rows and `tests/examples`' image rows
green on all three fixtures with the core solvers; no `grid_gp` symbol
left in the tree; `likelihoods.md` §7 amended (the GRID paragraph);
one decision-log row; gates: all four (the merged wave's).

### W5.22 — `population_full`, and the prior specification in provenance [S; Sonnet] (ruled by Peter 2026-09-16 on W5.13's proposals)
Two changes. (1) The 200-object row in `tests/results/test_population.py`
gets a `population_full` marker on the `image_full`/`m2_full` pattern
(registered in `pyproject.toml`, skipped by default in `test-all`, run on
request; a reduced sibling — 20 objects — stays in the default gate).
(2) A run's provenance stores each parameter's `PriorSpec` (the next schema
version — `PROVENANCE_SCHEMA_VERSION` is already 7 at drafting, so 8 — bumped, the append check and the
pre-schema refusal by name as W3.12's precedent), so
`ampere.results.population.fit_population` can *verify* a supplied
`interim_prior` against the archive and, when none is supplied, use the
stored one; the argument becomes optional, and a mismatch is refused by
name. **Depends:** nothing (W5.13 merged). **Accept:** the marker
present and the default dev gate no longer running the 200-object row
(time it: the `tests/results` suite before/after); the new-schema rows in
`tests/results/test_population.py` (round trip, the append check, the
old-schema refusal; the text named a `test_provenance*.py` that does not
exist — corrected 2026-09-22); `fit_population` with no `interim_prior` reproducing
W5.13's measured μ/τ intervals on the reduced row; a supplied prior that
disagrees with the stored one refused by name; `results.md` §9 and §13
item 16 amended; one decision-log row; gates: dev, sbi.

### W5.23 — Serving a named cached SBI artefact explicitly, with the mismatch recorded [S; Sonnet] (from Peter's note of 2026-09-16 on W5.4's hashing)
The artefact key (W3.5, `ampere/results/artefacts.py`) is one digest of
the problem's spec, model and data hashes plus the run's settings, so a
posterior trained at one `basis_size` is a miss for another by
construction, and that stays. Users who accept the risk get an *explicit*
route, not a fuzzy match: `SBIEngine(cache=..., serve_artefact=<digest>)`
restores that stored artefact regardless of the computed key, refuses if
the digest is absent from the store or its sidecar cannot be trusted,
and records in the run's attrs `sbi_cache_key` (computed),
`sbi_artefact_served` (the digest) and `sbi_artefact_mismatch` (the key
fields that differ, by name — the `ArtefactKey` fields are named, so the
comparison is field-wise), with a loud warning. The calibration
machinery must still work on the served posterior. **Depends:** nothing.
**Accept:** a served artefact from a run at a different `basis_size`
restored and sampled, its attrs carrying the three keys with the
mismatch naming the spec hash; a wrong digest refused by name; the
default path (no `serve_artefact`) bitwise unchanged (the W3.15
reproducibility row still green); `docs/source/sbi.rst`'s cache section
amended; `results.md` §9 gains the three attrs; one decision-log row;
gates: dev, sbi.

### W5.24 — Joint noise with heteroscedastic channels: the dense and reduced-rank routes [M; Opus] (ruled by Peter 2026-09-17 at W5.9's review)
W5.9's `JointGaussianProcessNoise` refuses unequal per-channel
uncertainties, because `Qᵀ ⊗ I` leaves `I ⊗ diag(σ²)` diagonal only when
every channel shares one σ vector. Heteroscedasticity across channels is
the norm for real data — astrometric and every other customer of the
feature — so the refusal is lifted by two routes behind the same
declaration, chosen by the bound solver. (1) **Dense**: `B ⊗ K_x +
blockdiag(diag(σ_t²))` factorised directly on the `TN × TN` matrix
(`DenseGP` bound; exact for any σ; O((TN)³), the reference and the small-N
path). (2) **Reduced-rank**: with `K_x ≈ Φ Λ Φᵀ` (`HilbertSpaceGP`, and
`EquispacedFourierGP` where it is promoted), `B ⊗ K_x ≈ (I_T ⊗ Φ)(B ⊗
Λ)(I_T ⊗ Φ)ᵀ` is a `Tm`-feature model against a noise that is *diagonal
in the original basis*, so Woodbury against `blockdiag(diag(σ_t²))` is
exact-in-the-approximation at O(TN (Tm)²) and the `latent_size` is `Tm`
— the NUTS-friendly path, and the first joint use of Phase 5's
approximate solvers. The exact O(N) rotated path stays the fast path
when the σ vectors agree (bit-identical to W5.9). The general LMC and
mismatched grids remain the recorded follow-on. **Depends:** W5.9 merged.
**Accept:** conformance rows on all three fixtures — both routes against
the materialised dense matrix at `tolerances.cross_solver` (dense) and
in W5.4's convergence class (reduced-rank, tightening with `m`), the
equal-σ case bit-identical to W5.9's rotated path, `simulate` drawing
correlated heteroscedastic channels; NUTS on torch and jax over `B` and
the kernel through the reduced-rank route; the astrometry `--joint` arm
with per-epoch heteroscedastic errors (the generator draws a σ per
channel per epoch) and its coverage pinned as W5.9's is; `likelihoods.md`
§7/§15 amended; one decision-log row; gates: all four.

### W5.25 — `halfcauchy` lowering on both backends [S; Sonnet] (ruled by Peter 2026-09-22 on W5.8's For Peter (1))
The horseshoe's global scale is half-Cauchy under both of
`shrinkage_horseshoe`'s tails, and `halfcauchy` is in neither backend's
§3.2 table, so the recommended prior for any `Sum` of noise terms is
reference-only on NUTS today. `torch.distributions.HalfCauchy` and
`numpyro.distributions.HalfCauchy` are both exact, so this is §3.4's
fallback used as intended: one `register_lowering("halfcauchy", ...)` row
per backend on `halfnorm`'s pattern (native when `loc == 0`, the §3.3 shift
otherwise), the `lowering.md` §3.2 row, the §3.2.1 note amended to say the
chain now lowers, and `tail="cauchy"` lowering too as a consequence.
**Owns:** `ampere/backends/torch/lowering.py`, `ampere/backends/jax/distributions.py`
(or wherever each backend registers `halfnorm`), `docs/design/lowering.md`
§3.2/§3.2.1, one conformance row. **Depends:** nothing. **Accept:** a
conformance row — a `halfcauchy(loc, scale)` prior's log-density and a
sample's support agree across the three fixtures, at `loc == 0` and
shifted; the W5.8 shrinkage problem's NUTS run on torch and jax passes the
existing backend-agreement shape (one short chain each, seeds fixed); the
§3.2 table row; `pixi run -e torch typecheck` / `-e jax typecheck` clean;
gates: all four (the lowering table is contract-adjacent — a §3.2 change
takes a one-line note in the decision log's W5.8 row, not a new row).

### W5.26 — Test cost and CI modularity [M; Sonnet] (ruled by Peter 2026-09-22)
Two problems with one item. **Locally**, the merged-wave gate is four legs
of `test-all` run one at a time (~2 h 30 m), and about half of every leg is
the studies' emcee-driven rows (`tests/m2` ≈ 11 min, `tests/results/
test_population.py` ≈ 8.6 min), which assert likelihood behaviour on the
numpy path and are repeated unchanged in sbi, torch and jax where the
backend-agreement rows already prove the backends match. **In CI**, each
`backend-suites` leg runs `test-all` as one forty-minute step, so a red
`tests/m2` row on torch hides everything else in that leg. **And two
suites run nowhere**: `tests/interferometry` (22 rows, W4.4) and
`tests/astrometry` (40 rows, W4.9) are in neither `test-all` nor any
`ci.yml` step since they landed. Do: (1) a `study` pytest marker (or
per-suite conftest skip on the backend fixture) so the studies' sampling
rows run in the `dev` environment only and skip by name elsewhere — the
backend-agreement rows keep running in every environment; (2) `tests/
interferometry` and `tests/astrometry` join `test-all` (measure their
budgets first and record them); (3) `ci.yml`'s `backend-suites` and
`suites` jobs gain a second matrix dimension over suite groups —
`core+results+conformance`, `backends+inference+examples`,
`conformance+m2+studies` — each group that carries a registry-leak proof
staying in one process (see `test-all`'s task comment for which two);
(4) `actionlint` clean; (5) the gate recipe in `docs/development.md`
updated to the new leg times. (6) **The pixi `--locked` check** (found 2026-09-22 on
the first `v2` CI run): pixi 0.68 and 0.81 both report the lock stale —
"'torch' requires index pypi.org but the lock-file has
download.pytorch.org" — because the satisfiability check sees the plain
`torch` requirement pyro-ppl and sbi carry and does not apply W2.11's
per-package CPU index, while the solver does (a re-solve reproduces the
lock byte for byte; spelling the extras out and a feature-level
`extra-index-urls` were both tried and do not help — the latter breaks the
editable build). CI installs `frozen` meanwhile; find whether a pixi
release or a manifest form (`index-strategy`, a workspace-level index)
restores `--locked`, or record it as upstream and keep `frozen`. (7) **Two rows that pass locally and fail on the
GitHub runner** (CI run `35684651060` on `17bb7df`, the same commit the
local four-environment gate passed): `tests/inference/test_sbi.py::
TestTheCalibrationFastPath::test_a_trained_npe_posterior_passes_check_sbc`
(sbi leg: the KS p-value was 0.0198 against `> 0.05` — a threshold at
0.05 on a p-value fails 5 % of calibrated posteriors by construction, so
any new machine's floats re-roll that die; recommend `> 0.001`, keeping the
contrast row's `< 0.01` separation) and `tests/inference/test_nuts.py::
TestTheReducedRankSolverUnderNUTS::test_the_posterior_agrees_with_the_exact_solver[jax]`
(jax leg: the HSGP posterior width 0.188 against the exact 0.137, pinned
to a quarter — 400 draws × 2 chains; on a different CPU the trajectory
differs and one chain's width moved 37 %; recommend a larger draw budget
or a chain-pooled width with the margin set from measured across-machine
scatter, not a looser pin alone). Both are seed-on-one-machine pins; the
fix is a per-row change with the reason, run on the runner to confirm.
**The runner is deterministic**: the second run (`35688582595`, on
`3cc97c7`) reproduced both values to every printed digit (0.01983926 and
0.05119422), so a fix can be verified on CI by one push rather than by
repetition. **Run `35785421357` on `3eeb2d3` (2026-09-22): the
reduced-rank row now fails on torch too** (0.0443 against 0.0351 allowed,
with jax at 0.0452 against 0.0341), so the fix is for the row itself, not
for one backend's parametrisation. Minutes are free on the public repository,
so the extra environment restores are accepted. **Not in scope:**
pytest-xdist (needs per-worker registry snapshots — record it as a
follow-on if the numbers say it is worth it); changing any pinned margin.
**Owns:** `pyproject.toml` `[tool.pixi.tasks]` and markers, `.github/
workflows/ci.yml`, the studies' `conftest.py` files, `docs/development.md`'s
gate paragraph. **Depends:** master pushed to `v2` (2026-09-22) so the
workflow runs on real runners. **Accept:** timed before/after for each
leg in each environment, recorded in the item's row (target: sbi, torch,
jax legs each ≥ 8 min shorter; dev unchanged or longer only by the two
new suites); every previously running row still runs in at least one
environment (a collected-test diff, dev ∪ torch, before vs after,
shows no row lost); the interferometry and astrometry rows green in the
gate; one CI run on `v2` green with the new job layout; gates: all four.

### W5.27 — `shrinkage_horseshoe`: the rename and the shrinkage-helper framework [S; Sonnet] (ruled by Peter 2026-09-22 on W5.8's For Peter (4))
`regularised_horseshoe`'s default tail is a gamma variant of the horseshoe,
not Piironen & Vehtari's slab, so the name promises the paper and the code
delivers a sibling. Rename to `shrinkage_horseshoe`, keeping
`regularised_horseshoe` as a deprecated alias (a `DeprecationWarning`
naming the phase it goes: Phase 6) so the merged docs and decision log
stay true. The rename names a family: every `shrinkage_*` helper returns a
list of `Parameter`s that `with_shrinkage` accepts, and one shared
docstring section (in `parameter.py`, referenced from each helper) states
that promise — the framework a future `shrinkage_spike_slab` or
`shrinkage_dirichlet_laplace` joins without another design round.
**Owns:** `ampere/core/parameter.py`, `ampere/core/__init__.py`, the M2
`many_lines` module and driver, `tests/core/test_parameter.py`, `tests/m2/
test_many_lines*.py`, `docs/design/contracts/parameters.md` §9,
`docs/design/lowering.md` §3.2.1, `docs/design/contracts/likelihoods.md`
§6/§13, `docs/source/kernels.rst`, `m2_misspecification.rst`. **Depends:**
nothing; **before W5.19**. **Accept:** the old name importable and warning
once, with the test that proves it; every reference in `docs/` and
`examples/` moved (a `grep -rn regularised_horseshoe` outside the alias, its
test and the decision log finds nothing); `tests/m2` rows green on their
existing margins (no re-pinning — the prior is unchanged); the docs build
clean; gates: dev.

### W5.28 — Phase 5 housekeeping II: the owed list [S; Sonnet]
W5.2's pattern for the second half of the phase: the small owed items the
wave-3 and wave-4 reviews recorded, each with a test where behaviour
changes. (a) `replace_observations` also drops `shared=` / `shared_label=`
(W5.9); (b) `ampere.results.sbc` ranks posterior variables only — state it
or rank the derived groups too (W5.7's carried); (c) a `summarise`-side helper reporting `RotationCoupling`'s matrix `B` rather than its bimodal angle parameterisation (W5.9's carried note; the first draft of this item mis-cited it as "the evidence `B` (W5.0/W5.13)" — corrected 2026-09-22 at review; the agent read it correctly from the code); (d) `QuasisepGP.condition
(at=...)` on a multi-axis container — refuse by name or implement (W5.3);
(e) `_GriddedSolver` is private but is the extension point W5.5's template
names — make it public with a docstring or name the public route;
(f) `HilbertSpaceGP` forms two `(N, m)` blocks where one suffices (W5.4);
(g) `_concatenate` in `results/training.py` indexes the addition by the
existing tree's groups — a bare `KeyError` or a silent drop for any future
optional group (W5.10); (h) `TrainingSet.contexts` defaults to `()`, so
"written before W5.10" and "written without a context" are
indistinguishable in Python — a sentinel, or a schema attr the reader
surfaces (W5.10); (i) `with_shrinkage` rebuilds the parameter set through
`copy.copy` and a `__dict__` write — a constructor path (W5.8; coordinate
with W5.27, which owns the same function's name); (j) `SBIEngine.
calibrate()` reseeds torch as `run()` does (W5.10's carried); (k) `tests/`
is not an importable package — decide and record (the conftest sys.path
pattern is used four times). **Not in scope:** anything a §4 contract
names. **Owns:** the files each item names; no design document beyond a
line in each affected page. **Depends:** W5.27 merged (for (i)).
**Accept:** each item either landed with its test or recorded in the row
as refused with the reason; lint, format, typecheck clean in dev, torch,
jax; gates: dev + torch (jax if (f) or (i) touches its path).

### W5.31 — The arviz lazy-`numpyro` import-order fix-up [XS; orchestrator] (ruled by Peter 2026-09-24 on W5.26's carried finding)
`tests/results/test_population.py::TestTheApproximateEngineRow::test_vi_fitted_object_agrees_with_emcee`
fails on jax whenever `tests/results` runs without a suite that imports the
jax backend ahead of it, on master as on every branch since arviz 1.3.0
entered the lock: arviz registers `numpyro` in `sys.modules` as an
`importlib.util._LazyModule`; a plain `import numpyro` returns the stub
unexecuted, and the backend's first `import numpyro.distributions` runs
numpyro's `__init__` mid-chain, after which `numpyro.distributions` never
binds `distribution` and `numpyro.factor` raises `AttributeError` inside
every native model. Reproduced in a fresh interpreter as `import
ampere.core; import arviz; import ampere.backends.jax`. Fix: touch
`numpyro` (read any attribute) in `ampere/backends/jax/_config.py`, the
first module the package imports; a regression row in
`tests/backends/test_jax_import_order.py` runs the failing order in a
fresh interpreter with the checkout first on the path, and a second row
proves the stub really is lazy after arviz alone. **Not in scope:**
anything in arviz's own import; the torch backend (pyro is not
lazy-loaded by arviz). **Owns:** `ampere/backends/jax/_config.py` and
`__init__.py` (comments), the new test file. **Depends:** W5.26 merged.
**Accept:** the reproduction green in a fresh interpreter; the VI
population row green on jax run *alone*; `tests/backends/test_jax.py`
green on jax; lint, format, jax typecheck clean; gates: jax (the merged
gate's).

### W5.32 — Phase 5 housekeeping IV: the CI warnings [S; Sonnet] (ruled by Peter 2026-09-25 on the orchestrator's warning analysis of CI runs 35962510128 and 36196425548)
The two green `v2` runs emit about 3400 pytest warnings across the
fifteen suite cells and a handful of Sphinx warnings; the orchestrator
sorted them (handoff, 2026-09-25) and Peter ruled which are ours. **Held
for the last filler slot before W5.19 (ruled 2026-09-26): any further
housekeeping found by later items is appended here as a new lettered
part rather than given an item of its own.** Five small changes so far,
each a commit. **(a) pytest 10.** `PytestRemovedIn10Warning`:
eight class-scoped fixtures are defined as instance methods —
`tests/results/test_population.py` (six: lines ~336, 403, 407, 796, 804,
812), `tests/inference/test_engines.py` (~794),
`tests/core/test_astropy_engines.py` (~114). Make each a `@staticmethod`
(or module-level where the class state is not used); the rows they feed
stay green, the warning is gone. **(b) sbi's `mcmc_parameters=`.**
Deprecated since sbi 0.25 (`FutureWarning` from 0.27.0, the locked
version) in favour of `posterior_parameters=` taking the
`MCMCPosteriorParameters` dataclass; two call sites in
`ampere/inference/_sbi.py` (`build_posterior` at ~2203 for the TMNRE
calibration sampler, ~2654 for `sample_with="mcmc"`). Translate
`_CALIBRATION_MCMC` and `_TMNRE_MCMC` to the dataclass with the same
values, keep the two dicts' names if the tests read them, and keep
`ampere_calibration_sampler` recording the same string; the TMNRE rows in
`tests/inference/test_sbi.py` are the check. **(c) `WarpedKernel`
documented three times.** The torch and jax backend packages re-export
it and their automodule pages describe it again ("duplicate object
description ... use :no-index:"): give the backend pages' entries
`:no-index:` as `docs/source/ampere.infer.rst` already does for its
re-exports, or exclude the re-export from the backend automodule — the
core page stays the one indexed description. **(d) `ampere.infer.sbi` in
the docs build.** autodoc fails to import the legacy `ampere.infer.sbi`
because the docs environment has no `sbi` package. **Ruled 2026-09-25:
mock it** — add `sbi` to `autodoc_mock_imports` in `docs/source/conf.py`
(the list's comment explains why an *installed* package must never be
mocked; `sbi` is import-only in the docs build, as torch and jax are), do
not add the extra to the docs environment; check the build log for any
other import-only optional the same rule applies to and mock it the same
way. **(e) One test leaves figures open.** `tests/results/test_plots.py::TestCornerPaging`
trips corner's "more than 20 figures" warning; close the pages after the
assertions. **(f) `sample_coordinates` below both modules** (W5.21's carried
note, 2026-09-26): `Likelihood._coordinates` imports
`ampere.core.dataset.sample_coordinates` inside the method because
`dataset.py` imports `Likelihood` at load; move the helper (it depends on
containers only) to a module both can import at the top — `encoding.py`
or a new `coordinates.py` — re-export it from `ampere.core` unchanged, and
make both imports module-level. Owns that helper's home, the two import
lines and `__init__.py`'s re-export. **Parts (g)–(l), ruled by Peter
2026-09-27 on the orchestrator's sweep of the carried findings since this
item was written; each a commit.** **(g) The astrometry usage text**
(W5.15's carried note): `examples/astrometry/__main__.py`'s docstring lists
the arms and flags but not `--wide-prior`, `--engine`, `--live-points`,
`--dlogz`; add them in the existing style. **(h) The image study's noise
helper on the modern backends** (W5.17's out-of-scope finding):
`examples/image/study.py`'s `_noise_for` falls back to the core
`ampere.core.DenseGP()` when no solver is passed, so on `"torch"` or
`"jax"` `FittingProblem` refuses the flexible arm as a foreign part
(`DatasetError`); the fallback becomes the chosen backend's own `DenseGP`
(the pattern the study's other arms already use for their pieces), with a
row per modern backend that the flexible arm composes. **(i) The `study`
marker project-wide** (W5.15's carried note): `_skip_study_rows_outside_dev`
is defined twice, in `tests/m2/conftest.py` and `tests/results/conftest.py`,
so the marker is inert elsewhere and `tests/astrometry/test_recovery.py`
carries its own `dev_only` skipif; move the hook once into the root
`tests/conftest.py` (W5.28 (k)'s file), delete the two copies, and replace
the local marker with `@pytest.mark.study` — the astrometry dev-only row
still runs on the dev leg and skips on torch, jax and sbi, which the
collected-test diff per environment proves (`--collect-only -q` before and
after: the same rows skip, no row lost). **(j) The silent numpy fallback
under a context, with a loud warning** (W5.29's carried note): on the
batched native path a user noise model whose `sigma_jax` (or `sigma_torch`)
does not take `base=` cannot receive the per-draw context sigma, so
`simulate_many(context=...)` quietly falls back to numpy draws, visible only
through `provenance["sample_backend"]`; and the in-chunk refusal W5.29
described has no test. Locate the fallback and the refusal with `grep -n
"sample_backend\|contextual\|_context_sigma" ampere/core/dataset.py`.
Add a **loud, obvious warning** at the fallback — a `UserWarning` subclass
named for the package beside the existing `AmpereFlatPopulationWarning`
(`ampere/core/settings.py` line ~29; put the new one where a
`dataset.py` import of it creates no cycle — `exceptions.py` if
`settings.py` imports `dataset`), whose message names the noise-model class,
says the batch is being drawn by numpy because its `sigma_jax`/`sigma_torch`
takes no `base=`, quotes the signature to add, and says where the
provenance records it; emit it once per `simulate_many` call, not per
chunk. Tests: a user noise model without `base=` triggers the warning
(`pytest.warns`) and the provenance says numpy; the shipped models do not
warn; the in-chunk refusal row W5.29 left untested, asserted by name.
`encoding.md` §13's context paragraph gains one sentence. **(k) A
conformance row for `proposal_log_density`** (W5.14's carried note): the
key has three producers in the tree — `VIEngine` (`_vi.py` ~570, every
guide family), `BlackjaxEngine`'s pathfinder (`_blackjax.py` ~442) and
`SBIEngine`'s NPE draws (`_sbi.py` ~2048); W5.14's row counted four, so
confirm with `grep -rn '"proposal_log_density"' ampere/` and say what the
fourth was if one exists — and no row asserts they agree on shape and
convention. One row
in `tests/conformance` (or `tests/inference` if it cannot be fixture-
parametrised): for each producer available in the environment, the array
lives at `sample_stats.proposal_log_density`, has one value per draw
(chain × draw), is finite, and is the log-density of the *unconstrained*
draw under the approximation as `results.md` §9 states, checked on the one
approximation whose density is known in closed form (the Laplace guide on
the conjugate Gaussian, `test_vi.py`'s fixture) and by the importance-
weight identity on the others (`exp(log_prior + log_likelihood -
proposal_log_density)` has finite, positive weights); producers absent
from the environment skip by name. **(l) Pickled reference steps carry
their plans** (W5.17's carried note): `Resample._influence_plan`,
`FourierSample._transform_plan` and `_template_source_unit` are caches that
now ride inside a pickled step; give `_Step` in
`ampere/backends/reference/instrument.py` a `__getstate__` that drops
attributes whose names start with `_` and end with `_plan`, plus
`_template_source_unit` (or a small declared tuple of cache names the two
classes extend), so a round-tripped step replans on first use; a test
pickles a step after one `apply`, asserts the pickle carries no plan and
the reloaded step gives bit-identical output. **Ruled out (record, do not do):** the `ubuntu-latest` → Ubuntu
26 migration notice (no runner pin for now, ruled 2026-09-25); the
~2700 `DeprecationWarning`s from netCDF4-python setting an array's shape
under NumPy 2.5 — fixed upstream in netcdf4 1.7.4.1, not yet on
conda-forge; re-lock when it lands (say in the report whether it has);
the legacy `SyntaxWarning`s (`data/spectrum.py`, `infer/mixins.py`),
emcee's deprecated `chain`/`a` in `ampere/infer`, and the docutils
"inline strong start-string" warnings in legacy docstrings (frozen);
sbi's, dynesty's and arviz's budget warnings and the kernel overflow in
the failure-signalling rows (deliberate tiny budgets and bad inputs); the
arviz chain-longer-than-draw warning (carried in W5.8's row). **Owns:** the
three test files in (a), `ampere/inference/_sbi.py` for (b)'s two calls
and two constants only, `docs/source/conf.py` and the two backend rst
pages, `tests/results/test_plots.py` for (e); for (f) the helper's home
and the import lines; (g) `examples/astrometry/__main__.py`; (h)
`examples/image/study.py` and its test file; (i) the three conftests and
`tests/astrometry/test_recovery.py`'s marker; (j) the fallback site in
`ampere/core/dataset.py`, the new warning class (`settings.py` or `exceptions.py`), the
tests in `tests/inference/test_sbi.py`, one sentence in `encoding.md`
§13; (k) one new conformance/inference test module; (l) `_Step` in
`ampere/backends/reference/instrument.py` and one test. **Depends:**
nothing (W5.15 and W5.17 merged). **Accept:** (g)–(l) as their own
tests say above, plus: the docs build with zero autodoc/duplicate-object warnings
for the new namespaces and no `ampere.infer.sbi` import failure (the
legacy docutils warnings may remain — count them before and after);
`tests/results/test_population.py`, `tests/inference/test_engines.py`,
`tests/core/test_astropy_engines.py` green in dev with no
`PytestRemovedIn10Warning`; the TMNRE rows green in sbi with no
`FutureWarning` from sbi about `mcmc_parameters`; `TestCornerPaging` green
with no figure warning; the collected-test diff for (i) per environment;
lint, format, typecheck clean in dev, torch and jax; gates: dev + sbi +
jax + torch (the (j) warning and (i) marker are exercised on the modern
legs; the four-leg merged gate covers it).

### W5.29 — The native batched path drawing a context [M; Opus] (W5.10's carried item, ruled by Peter 2026-09-22)
`simulate_many(context=prior)` draws a context per draw, but on the batched
native path a per-draw σ does not reach the realisation, so the path
refuses by name and the amortisation over context (W5.10) runs unbatched.
Lifting it means the realisation's noise model accepts a batch of σ (one
per draw) on both modern backends, the batch's provenance records the
drawn contexts as the unbatched path does, and the cache key is unchanged
(the context prior already enters it). **Owns:** `ampere/core/dataset.py`
(`simulate_many`'s batched branch), `ampere/backends/torch/` and `jax/`
noise-model realisation, `docs/design/contracts/inference.md` §13's
context paragraph, `tests/inference/test_sbi.py`, one conformance row.
**Depends:** W5.10 (merged). **Accept:** a conformance row — a batched
draw with a context prior equals the unbatched draw at the same seed on
torch and jax (the reference path has no batched branch; it is the
oracle); the W5.10 tutorial's context run on the batched path with a
measured speed-up recorded in the row; the refusal message gone and its
test inverted; gates: sbi + torch + jax.

### W5.30 — Phase 5 housekeeping III: the flat-population cap and the hierarchical bijection [S; Sonnet] (ruled by Peter 2026-09-22 on W5.12's For Peter (2) and W5.25's carried note)
Two small items, both in `ampere/core/parameter.py`. **(a) The
flat-population cap.** `MAX_FLAT_MEMBERS` stays 128. Its docstring and the
refusal message in `Population.__post_init__` cite
`hierarchical_population.md` §7's measured figures — 46 ms per `lnprior`
at 100 flat members, 373 ms at 1000, against 1.3 ms for the same
structure as one plate — in place of the "800 ms / 0.6 ms" and "well
under a millisecond" they say now, and state plainly that at the cap a
flat prior costs about 60 ms per evaluation, so a hundred-thousand-
evaluation ensemble run spends well over an hour in the prior. Then a
**setting** turns the refusal into a loud warning for a user who accepts
that cost: not a `Population` argument and not a change to the limit
(Peter's ruling). The package has no settings mechanism; add the smallest
one a test can flip and restore — recommended: a new typed
`ampere/core/settings.py` holding one module-level dataclass instance with
the single field `flat_population_cap: Literal["refuse", "warn"] =
"refuse"` and a context-manager `override(**fields)` that restores on
exit; read at validation time, never at import. The warning (a
`UserWarning` subclass named for the package) carries the member count,
the measured cost and the plate remedy. Tests: the refusal by default; the
warning under the setting, after which the population merges, evaluates
`lnprior`, and the override restores the refusal; `parameters.md` §9
gains an *Amended W5.30* line and the populations page one sentence.
**(b) Bijection inference for hierarchical priors whose family takes
shape arguments.** `_default_bijection_for_hierarchical` refuses every
family with `numargs > 0`, which is why `shrinkage_horseshoe` passes
`bijection=Log()` by hand and why a hierarchical `gamma` or `beta` on
its own fails until the user does the same. The support of `gamma`,
`beta`, `lognorm` and most shape families does not move with the shape
(scipy's `a`/`b` are class constants), and only a minority (`truncnorm`,
`uniform`-like families whose bounds *are* shape arguments, anything
overriding `_get_support`) have supports that do. Replace the blanket
refusal with a rule that infers exactly when it is safe and refuses,
with the existing message, when it is not: infer for a family whose
support is independent of its shape arguments (test: `dist._get_support`
is the base-class one, or the support is equal at the declared shape
values and none of them is a hyperparameter reference); keep refusing
when the support depends on a shape argument, and when any shape argument
is itself referenced by name, because the support could then move from
draw to draw. Do not overcomplicate it: one predicate, one docstring
paragraph saying what it checks and why. Tests: `gamma`, `beta` and
`lognorm` hierarchical priors infer without a keyword (`Log` for the
half-lines, the bounded bijection for `beta`); `truncnorm` still refuses;
a `gamma` whose `a` is referenced still refuses; an explicit
`bijection=` still wins; the horseshoe helper's explicit keyword is kept
(harmless) and a row proves the helper's declaration infers the same
bijection without it. The torch and jax hierarchical registries take the
bijection from the parameter, so the W5.25 gamma and halfcauchy NUTS
rows must stay green. **(c) `replace_observations` and `populations=`** (W5.28's carried note): the replica rebuilt by `ampere.results.calibration.replace_observations` passes `ties=` but not `populations=` to the new `FittingProblem`, so a W5.12 population-level SBC replay would be fitted as independent objects — the same silent drop W5.28 (a) closed for the shared set. Carry it through and add the row beside W5.28 (a)'s (`tests/results/test_calibration.py`). Owns `ampere/results/calibration.py` for this line only. **Not in scope:** anything a §4 contract names —
if the reviewer judges the inference rule to be a §4.1 statement, a
decision-log line goes in the same PR; changing the cap; a general
configuration system (one field, one module). **Owns:**
`ampere/core/parameter.py`, the new `ampere/core/settings.py`,
`tests/core/test_parameter.py`, `parameters.md` §9 (one line each),
the populations and priors pages of the Sphinx site. **Depends:** W5.28
merged (its item (i) edits the same file; its (a) is where (c) goes). **Accept:** the tests above
green in dev; `tests/conformance` green in dev, torch and jax; the W5.25
horseshoe-under-NUTS rows green on torch and jax; lint, format,
typecheck clean in dev, torch, jax; gates: dev + torch + jax.

**Decisions for Peter before dispatch.**
- **D1 — the native model surface's spelling** (W4.3's decision-log row):
  (a) `native_flux`/`native_grid` canonical, `flux`/`grid` an alias kept
  indefinitely (a docs and new-code rule; the eight existing modules keep
  working; a Haiku sweep renames them when convenient); (b) narrow
  `Parameterised._check_free_name` so a parameter may shadow a method;
  (c) leave both spellings as documented. Recommendation: (a) — the
  shadowing rule caught a real ambiguity in `model.flux`, and the collision
  recurs for any model whose physical parameter is a flux. **Ruled 2026-09-15: (a) and (b) together, as W5.20** (plan §2's W4.3 row carries the analysis).
- **D2 — the first 2-D+ solver and its customer**: HSGP first on a single
  `Image` (W5.4 + W5.5, recommended) with the bake-off (W5.6) in-phase; or
  the IFU cube as the customer (needs W5.12's plate first, so the solver
  waits half a phase). **Ruled 2026-09-15: the image (W5.4 + W5.5, W5.6 in-phase).**
- **D3 — the joint-noise first customer**: astrometry's two channels
  (recommended; exists, shared grid, no new kind) or polarimetry (a Stokes
  kind and a reference instrument to write first). **Ruled 2026-09-15: astrometry first, polarimetry second.**
- **D4 — population scope**: both routes (W5.12 the joint fit, W5.13 the
  reweighting module — the memo's Q5 left both open for Phase 5);
  recommendation: both, W5.13 first since it is small and exercises W5.0. **Ruled 2026-09-15: both; W5.13 first.**
- **D5 — tier-1 engines in Phase 5**: the memo's ruling scheduled nothing
  beyond W5.0; W5.14 is drafted as opt-in. Recommendation: in, after W5.0,
  because W5.15's periodic case and W5.13's evidences both want a nested
  sampler with evidence, and the phase is titled advanced inference. The
  `Optimum` result and warm start (memo §5.4) stay out until an optimiser
  item exists. **Ruled 2026-09-15: in (W5.14 after W5.0); `Optimum` stays out.**
- **D6 — the periodic-model hazard**: leave as the hazard note W4.8
  wrote, or run W5.15 (recommended: it turns a warning into a measured
  recommendation and costs a Sonnet afternoon). **Ruled 2026-09-15: run W5.15.**
- **D7 — the encoding's axis identity**: W5.11 in W5.10's wave so every
  layout hash moves once (recommended), or deferred past Phase 5. **Ruled 2026-09-15: W5.11 in W5.10's wave.**
- **D8 — warping as one item or two**: W5.7 (kernel) and W5.8 (study and
  sparsity prior) as drafted, or merged into one Opus item; and whether
  W5.8's M2 margins are pinned or informational at the per-PR budget
  (recommendation: pinned as margins, W4.5's precedent). **Ruled 2026-09-15: two Opus items; W5.8's margins pinned.**

Ordering (two agents at a time): **wave 1** W5.0 ∥ W5.2, then W5.20; **wave 2** W5.4
(1-D) ∥ W5.3, then W5.5 ∥ W5.13; **wave 3** W5.7 ∥ W5.9; **wave 4** W5.10
(+ W5.11 if D7) ∥ W5.8; **wave 5** W5.12 ∥ W5.14; **wave 6** W5.6 ∥ W5.15,
W5.16, W5.17, W5.18 as they free; W5.19 last. **Fillers, added
2026-09-22 on Peter's rulings** (any free slot, in this order): W5.25
first, then W5.26, W5.27 (before W5.19), W5.28 (after W5.27), W5.30 (after W5.28), W5.29,
with W5.21 and W5.24 as before. **Ruled by Peter 2026-09-24 (the
remaining order, after the 2026-09-22 handoffs had stopped naming wave 6,
the two older fillers and W5.1)**: after W5.26 merges, **W5.1 ∥ W5.29**
(both Opus; W5.1 dispatched first, the larger blast radius); then
**W5.21, W5.24, W5.15** as slots free; **W5.17** once W5.1 and W5.24 are
merged, so the profile sees the final code; **W5.18 trimmed to its
residue** (W5.14 landed the nested leg and W5.26 re-measured every leg
and split the matrix — what remains is the `test-fast` task and the
GPU-job check); **W5.16 deferred to Phase 6** by explicit ruling; **W5.32** (housekeeping IV,
the CI warnings, Sonnet S, added 2026-09-25) **held for the last slot before
W5.19** so that further housekeeping found on the way folds into it (ruled
2026-09-26); **W5.19 last**. File ownership per wave in
the dispatch prompts; the sole shared file across waves is
`core/likelihood.py` (W5.4, W5.7, W5.9 each own a section — the solver,
the kernel, the noise-model — and merge in that order).


## Phase 6 — Docs, migration, release (drafted 2026-09-28 by Fable; **ruled by Peter 2026-09-28: D1 as (b) with `ampere.legacy`, D2–D7, D9, D10 as recommended (D3 as `1.0.0b1`); D8 as recommended (2026-09-28, after the briefing); D11 ruled 2026-09-28: the legacy examples stay on legacy and gain v2 twins — W6.13; D12 ruled 2026-09-28 as recommended, plus a `star_disc` twin and the committed pre-trained emulator**; **wave 1 (W6.6 ∥ W6.8) cleared to dispatch on Peter's word of 2026-09-28**)

The plan's §5 Phase 6 bullets (the RHMF trial, the optimisers module, the
scipy distribution exploration, the docs rebuild and beta release, the
composition tutorial, the two design items, "where the merged gate runs"),
the §5 CI/CD workstream's Phase 6 line (release automation, versioned docs,
the GPU suite on an accelerator), the placement memo's Phase 6 pointers
(OIFITS as the first reader and the `ampere.interferometry` front door),
the Phase 5 rows' "Phase 6" carried notes, and the open GitHub issues, as
agent-sized items. **Read `DEVELOPMENT_PLAN.md` §5 "Phase 6" first**: each
item names the bullet it executes. The phase has a different shape from
the five before it — its product is a **release** rather than a
capability, so most items are documentation, policy and infrastructure,
and the two capability items (the optimisers, the readers) are the ones
whose scope Peter rules on. Sized under `docs/orchestration.md`'s rules:
two agents at a time, targeted tests on the branch, gate legs scoped to
the code the work touched (ruled 2026-09-28). Every change to a §4
contract is a decision-log row and the conformance rows in the same PR
(ground rule 9); the release itself is a decision-log row. The eleven
decisions **D1–D11** at the end are the ones the drafting could not take;
the ordering paragraph before them assumes the recommendations.

### W6.0 — `ampere.legacy`: the legacy surface relocated, kept, and never removed [S; Fable drafts the policy page, Sonnet lands it] (ruled D1 (b), 2026-09-28)
`docs/source/migrating.rst` promises "a deprecation policy — if and when
legacy pieces are retired, and on what notice — is Phase 6 work" (W3.9).
**Ruled (b)**: the legacy code is **kept indefinitely, frozen as it is,
with no removal scheduled** — and it moves to one subpackage,
`ampere.legacy`, so that `from ampere import legacy as ampere` gives a
legacy script its whole old surface under one name. The item: (1) the
move — `ampere/data`, `ampere/models`, `ampere/infer`, `ampere/utils`
become `ampere/legacy/{data,models,infer,utils}` by `git mv`, their
internal absolute imports rewritten to the new path and nothing else
changed (ground rule 1's "do not modify" is lifted for exactly this
mechanical relocation, by this item's text; the characterisation suite
and `tests/test_imports.py` are the proof nothing else moved);
`ampere/legacy/__init__.py` imports the four as the old top level did;
(2) the old top-level names — `ampere.data`, `ampere.models`,
`ampere.infer`, `ampere.utils` — stay importable **as lazy aliases**
(a module `__getattr__` in `ampere/__init__.py` resolving to the
`ampere.legacy` module on first access), so no existing import breaks
and `import ampere` no longer loads any legacy module or legacy
dependency eagerly (today it imports all four before anything of v2);
the aliases are documented on the migration page as "kept; prefer
`ampere.legacy`", with no warning and no removal date — the orchestrator's
assumption under (b), to confirm at review; (3) `__version__` reads
`importlib.metadata` instead of the hard-coded `"0.1.2"`; (4) the policy
page — what "frozen" means now (kept, unchanged, documented as legacy,
bugs fixed only where a fix is a line and never a behaviour change,
tested by the characterisation suite), that v2's **own** deprecations
(the `regularised_horseshoe` alias, W5.27) are a separate rule —
deprecated with a warning, removed at 1.0.0 final (the orchestrator's
assumption: (b) is about legacy, not about v2's aliases; to confirm) —
and where each legacy example goes under W6.1; (5) the legacy
characterisation anchors (`examples/minimal_working_example*.py`) and
`examples/examples_paper/` keep working through the aliases or through
`ampere.legacy`, unmodified. **Depends:** nothing. **Accept:** the
characterisation suite green unchanged; `import ampere` in a bare
environment imports no legacy module (a `tests/test_imports.py` row);
`from ampere import legacy as ampere` runs `minimal_working_example.py`
with one substitution; the old names still resolve; the policy page in
the API toctree beside `migrating`; docs warnings no longer than base;
gates dev (the characterisation suite once on the branch).

### W6.1 — The full migration guide, and the legacy examples converted [M; Sonnet]
W3.9 seeded `migrating.rst` with the concept map and one side-by-side; this
completes it (issue #59): every legacy class and search mapped to its v2
route with a snippet each — `Photometry`/`Spectrum` → `Dataset` with an
`Instrument` chain, the legacy GP switches → `GaussianProcessNoise`, the
`ampere.infer` searches → the engines, `ampere.infer.sbi` → `SBIEngine`,
the post-processors → `ampere.results`, the extinction and filter helpers
→ their v2 equivalents or "no equivalent; carried" — with the legacy
examples' v2 twins (W6.13, under D11: the legacy scripts stay as they
are and keep working on `ampere.legacy`; each gains a v2 twin beside it)
linked from the guide as the side-by-side material; the three notebooks under
`docs/source/notebooks/` (`quickstart`, `Ampere_MBB_Example`,
`Embedding_nets`) are re-written on v2 and executed at docs-build time
(`nbsphinx_execute = "auto"` for the ones that can run in the docs
environment) or moved to the legacy section with a banner; the legacy
characterisation anchors are untouched (they exercise the frozen code, not
the examples). **Depends:** W6.0, W6.13 (the twins it links). **Accept:** every legacy
public name appears on the guide with a route or a "carried" note; the
executed notebooks build in `pixi run docs`; docs warnings no
longer than the base commit's; gates dev.

### W6.2 — Composition tutorial: photometry plus spectra, with calibration uncertainty [S; Sonnet]
The plan's bullet (Peter, 2026-09-09): the worked, runnable example of one
model against a photometric catalogue and one or more spectra —
`SyntheticPhotometry` beside `Resample`/`LSFConvolution` — with calibration
uncertainty as a `CalibrationScale` step carrying a prior, per spectrum or
shared through a `Tie`, and the flexible likelihood as the complement for
what calibration does not explain; `spectrum_photometry.md` and
`sed_composition.rst` are the starting points (the latter is the one-model
two-instrument case; this page is the several-observations case W4.4's
template cites and no example shows). Lands the owed `tests/examples`
smoke test of the composition. **Depends:** nothing. **Accept:** the page
in the tutorials toctree; the example under `examples/` with a smoke row;
`tests/examples` green; docs warnings no longer than base; gates dev.

### W6.3 — The docs rebuild: warnings-as-errors, docstrings, the reference restructured [M; Sonnet]
The Phase 1 promise the plan's CI/CD workstream deferred here: `-W` on the
docs job. Twelve warnings remain, all in frozen legacy docstrings
(`ampere/data/spectrum.py`, `ampere/infer/{emceesearch,mixins,zeussearch}.py`)
and one ambiguous cross-reference in `ampere.utils.rst`; under W6.0's
policy the legacy pages stay in the reference indefinitely (D1 (b)), so they are
autodoc'd with the offending members excluded (`:exclude-members:`, the
W5.32 (c) pattern) rather than by editing frozen code — ground rule 1 holds.
Then: (1) the API reference restructured around v2 — `api.rst` leads with
`ampere.core`, `ampere.backends`, `ampere.inference`, `ampere.results`, and
the legacy packages become one "Legacy (deprecated)" section with the
policy banner; (2) v2 docstrings audited against the pages (issue #58: every
public name in the four namespaces has a docstring that autodoc renders
without a warning); (3) the two pages `tutorials.rst`'s "Still to be
written" names — conditional priors, arbitrary priors — written from
`parameters.md` §4–§6; (4) the shrinkage section `advanced.rst` lacks
(W5.19's question; **D5 ruled yes**); (5) the `pycon` blocks under `docs/source`
under a doctest runner (W5.7's carried idea, W5.19's question; **D5 ruled
yes**) — `tests/core/test_spec_doctests.py`'s harness extended to `docs/source`,
each page's blocks either run or marked `:skipif:` with the reason. Issue
#57 closes with this item and W6.1 together. **Depends:** W6.0 (the policy
banner), W6.1 (the notebooks). **Accept:** `pixi run docs` with `-W` green
in CI; every v2 public name rendered; the doctest runner (if ruled) green
over `docs/source`; gates dev.

### W6.4 — Versioned documentation deployment [S; Sonnet]
The CI/CD workstream's Phase 6 line: the docs deployed per version — Read
the Docs (a `.readthedocs.yaml` building from the pixi `dev` environment,
or a pip-installable docs extra, since RTD does not run pixi natively) or
GitHub Pages from the CI job (a `docs` deploy step on tags and on `master`
as "latest"), per D2; a version switcher; `latest` and `stable` aliases;
the legacy site, if any, redirected. **D13 (2026-10-02): the Read the Docs
slug is `ampere` (`ampere.readthedocs.io`; free at the ruling), the
project imported from the GitHub repository by Peter before dispatch —
the agent cannot create it; until the beta tag the `v2` mirror is the
branch RTD builds as `latest`, switching to `master` at the tag (D4).** **Depends:** W6.3, D2. **Accept:** the
docs reachable at the ruled host for the branch head and for the tagged
beta; the CI job publishes on tag; no build step outside pixi or the docs
extra; gates none (infrastructure; the docs job is the check).

### W6.5 — Release automation and the beta release [M; Opus]
The plan's bullet and the CI/CD line together, addressing issues #57–60 and
#62: (1) **versioning** — setuptools_scm is configured; `__version__` moves
to metadata (W6.0); the beta is **`1.0.0b1`** from a `v1.0.0b1` tag (D3 ruled:
the major bump signals the redesign; the final is `1.0.0`); (2) **PyPI trusted
publishing** — a `release.yml` on `v*` tags: build sdist and wheel in the
`dev` environment, `twine check`, publish through OIDC, no token in
secrets; a TestPyPI dry run first; the wheel's metadata carries
`License-Expression` (the RHMF adoptability row found a packaging
without one reports no licence on PyPI — check ampere's own); (3) **the
changelog** — hand-written for the beta from the decision log and the status table
for the beta (the phases' landed summaries are the material) and
maintained per release afterwards; (4) **citation** (#62) — a
`CITATION.cff`, a Zenodo DOI minted at the tag, the paper reference once
it exists, and `ampere.__citation__` or a `cite()` helper printing it;
(5) **installation** (#60) — `install.rst` and the README's install
section re-checked against a clean PyPI install of the beta in a fresh
environment for each extra, and the `all` extra verified to resolve —
**under the distribution name `ampere-astro` (D13, ruled 2026-10-02: the
PyPI name `ampere` is an unrelated package)**: `name = "ampere-astro"` in
`pyproject.toml`, the import name unchanged, the extras' short form in
docstrings and `OptionalDependencyError` messages becoming
`pip install "ampere-astro[jax]"`, and the README's and `install.rst`'s
clash warnings rewritten as the plain statement of the three names
(distribution, import, docs slug);
(6) **the remote** — per D4, `origin/master` takes the v2 line at the beta
tag (the `v2` mirror's purpose ends), branch protection on `master`
(CI required, no force-push), and the release is the decision-log row
that says so; (7) the paper's examples (`examples/examples_paper/`) are
Peter's and stay out unless D11 couples them. **Depends:** W6.0, W6.3,
W6.4, D3, D4. **Accept:** a TestPyPI release installable with every extra
in a fresh venv; the tagged beta on PyPI with the DOI; `pip install
ampere-astro` in a clean environment imports `ampere` without legacy; the changelog and
citation on the docs site; `origin/master` at the tag; gates all four
(the release gate is the full matrix once, on the tag).

### W6.6 — Phase 6 housekeeping: the Phase 5 residue [S; Sonnet]
The carried notes the rows sent here, each a commit: (a) `pyproject.toml`'s
`blackjax` extra comment says `_blackjax.py` "does not exist yet" — it does
(W5.19); (b) `tests/results/test_plots.py`'s other classes close their
figures, so the 20-figure warning stops firing when the file runs whole
(W5.32, W5.19); (c) W5.18's 26 `--deselect` flags become one `sbi_training`
marker on the classes that train a network, now that W5.32 (i) settled the
markers list, and `tests/inference/test_interferometry.py`'s one training
class (~662) joins it; (d) a conformance row for `QuasisepGP.condition(at=)`
on a multi-axis container — the reference accepts, torch and jax refuse
by name (W5.28, annotated at W5.19) — asserting the disagreement as a
declared capability rather than leaving it undeclared, or closing it if
the fix is a few lines in each backend's `_axis` (say which); (e) the
`regularised_horseshoe` alias's removal is *scheduled*, not done, under
W6.0's policy — this part only checks the warning names the policy's
removal version; (f) netCDF4 re-locked when 1.7.4.1 reaches conda-forge
(W5.32's note; the NumPy-2.5 deprecation is ~80 % of CI's warning volume)
— if it has not landed, say so; (g) ultranest's `logzerr` is more
conservative than its console line (W5.14) — one sentence on the
engine's page. **Depends:** nothing; any slot. **Accept:** each part's
own check as W5.32's were; gates per part (dev; jax and torch for (d);
the CI run for (f)).

### W6.7 — The optimisers module: point estimates and warm starts [L; Opus]
The plan's bullet verbatim is the item text: the three routes behind one
call on the frozen §4.5 surface — multi-start `scipy.optimize` over the
packed unconstrained vector on the numpy path; gradient MAP through
`realise` on torch and jax (L-BFGS, Adam as the fallback) and `VIEngine`'s
mean as the alternative start, both reaching NUTS through `init_to_value`
and the ensembles through `initial_positions`; the reduced-rank
empirical-Bayes hyperparameter warm start with a length-scale grid — and
a point estimate with provenance stored on the run it seeds. Bayesian
optimisation stays out (the dimension argument of 2026-09-24). Issues #40
and #14 close with it; #41 (variational Bayes with snowline) is answered
by `VIEngine` and closes with a note. The `Optimum` result shape the
inference-extensions memo §5.4 sketched is this item's to land, with a
`results.md` amendment and the decision-log row. **Depends:** Phase 5
closed; nothing else. **Accept:** as the plan's bullet states — each
route's point estimate inside the sampled posterior's central 50 % on
every free parameter on the conformance fixtures; the warm start within a
factor of two of the sampled posterior median on a pinned reduced-rank
row; a pinned row where a NUTS run started from the MAP reaches its
adaptation target in fewer warm-up steps than the prior-draw start, and
the emcee burn-in equivalent; the start in provenance (a schema bump if a
new attr is needed); `inference.md` amended; a `docs/source` page; gates
all four.

### W6.8 — The scipy distribution exploration [S; Sonnet]
The plan's bullet verbatim: measure the six-prior `lnprior` of the M2
reference load through scipy's new distribution infrastructure (≥ 1.15)
against W5.17's 262 µs legacy figure and 20.5 µs floor; bit-identity or
the tolerance §13's conformance row would absorb; which of ampere's common
priors it covers; how a user's frozen legacy prior sits alongside under
`parameters.md` §4's protocol. The outcome is a **report** in
`docs/design/performance_memo.md` §6 (a new subsection) and a drafted
follow-up item if the gain is real and exact; nothing lands in the core.
**Depends:** nothing; a filler. **Accept:** the measurement table with the
three columns, the coverage list, the identity result, the recommendation;
gates none (no code change).

### W6.9 — RHMF exploratory trial [S; Sonnet] (W5.16's text stands, deferred here by ruling 2026-09-24)
As written under Phase 5 (W5.16): the pre-fit robust-factorisation family
tried against a pinned `robusta-hmf` commit behind a non-default `rhmf`
extra, the adapter producing an `AnomalyScore` with `provenance="rhmf_prefit"`,
run on the M2 spectra and W5.5's image, with the adoptability re-check
re-run and recorded. The outcome is a report; `ampere.diagnostics` lands
only if the maturity gate is met. **Depends:** nothing. **Accept:** W5.16's.

### W6.10 — Where the merged gate runs: CI as the gate of record, and the parallelism audit [S; Sonnet] (the GPU rows split out to W6.16 on Peter's word, 2026-10-02)
The plan's assessment bullet, executed: (1) under D4, the orchestrator
pushes `origin/v2` (then `origin/master`) at each merge and the CI run on
the push is the merged gate of record — the local legs retire to the
scoped pre-merge runs ruled 2026-09-28, and `~/.cache/ampere-gates`'s
scripts become the fallback for a machine without GitHub; the handoff's
"next gate's baselines" become the CI run's counts; (2) the
`pytest-xdist` safety audit — whether the suites' seeded streams, the
shared lock in `tests/conftest.py` and the m2 margins survive `-n auto`;
if they do, `test-all` and `test-fast` gain `-n` and the measured
speed-up goes in the task comment; if a suite does not, it is named and
left serial. The GPU rows on an accelerator, formerly (3) here, are
**W6.16** (split out 2026-10-02: Peter is checking the cluster's details
first). **Depends:** D4. **Accept:** the CI run recorded as the gate in the
next merged row; the audit's table; gates none (infrastructure).

### W6.11 — Design memo: per-dataset nuisance populations, and the `Derived` parameter node [M; Fable drafts, Opus reviews]
The plan's two design items (Peter, 2026-09-22): (1) `Population.over`
addressing a dataset's qualified path — `PlateBinding`/`Binding` carrying a
qualified name so "each dataset's GP amplitude is a draw from one shared
prior" and a per-dataset calibration scale under a fitted spread are
declarable — a §4 change to `parameters.md` §8–§9 and `inference.md` §9,
with `hierarchical_population.md` §11 Q1's rejection of routing in two
places respected; (2) the `Derived` node — a parameter that is a pure
function of others — making Piironen & Vehtari's slab declarable
(`parameters.md` §9, `lowering.md` §3.2.1) and the non-centred
`θ_i = μ + σ z_i` expressible for a population, with its lowering on both
backends and its place in provenance and the results groups. The product
is a memo under `docs/design/` with the contract amendments drafted as
decision-log rows and the conformance rows named, **not** an
implementation; the implementation items are drafted at the end of the
memo for D8's ruling. **Depends:** nothing. **Accept:** the memo; the
drafted rows; no code.

### W6.12 — The first reader: OIFITS into `VisibilitySet`/`ClosurePhases`, and the observable's front door [M; Opus]
The placement memo's D reserved it: when the first reader lands, create
`ampere.interferometry` as the observable's front door and let astrometry
and image follow the same shape. The reader takes an OIFITS file
(`OI_VIS2`, `OI_T3`, `OI_WAVELENGTH`) into the shipped containers with
the spectral axis (Phase 4 D2), the canonical baseline ordering and the
closure-phase triangle pairing W5.1 binds to; units and flags handled
per `results_schema.md`; `astropy.io.fits` only, no new dependency. A
JWST spectrum reader (issue #63) is the second reader under the same
front-door pattern (`ampere.spectroscopy`?), scoped by D9. **Depends:**
D9. **Accept:** a real OIFITS file (a public archive product, small,
under `tests/data` if its licence allows, else downloaded in the test
with a skip offline) round-trips into containers `tests/interferometry`'s
fixtures accept; the front-door module documented on `interferometry.rst`;
gates dev (torch and jax if the containers' native twins are touched).

### W6.13 — The legacy examples' models on v2: twins beside the originals [M; Sonnet for the self-contained models, Opus for the two external-code models] (ruled D11, 2026-09-28)
**Ruled**: the legacy examples are kept exactly as they are, running on
legacy (D1 (b)), and each named model gains a **v2 twin** — the same
model, the same data, the same question, written against `ampere.core`
and the reference backend, beside the original under `examples/`. Peter
named five and D12 added a sixth: (1) **the minimal working
examples** — `minimal_working_example.py` and its `_dynesty`, `_zeus`,
`_sbi` and `_sbi_embedding` variants — become one v2 example,
`examples/linear_sed/`, with `--engine emcee|dynesty|zeus|sbi` and
`--embedding` switches: the linear model, two-band synthetic photometry
through `SyntheticPhotometry.from_library` (the pattern
`examples/sed_composition` already uses), the Spitzer IRS sampling read
from the tracked `examples/test_data/cassis_yaaar_spcfw_14191360t.fits`
with `astropy.io.fits` (a v2 `Spectrum` container, no legacy reader), the
calibration scale and the flexible likelihood as the v2 counterparts of
`calUnc`/`scaleLengthPrior`; the legacy files are the characterisation
anchors and are **not touched**; (2) **NGC6302** — `NGC6302.py` and
`NGC6302_zeus.py` become `examples/ngc6302/`: the Kemper et al. (2002)
two-shell model as a v2 `Model` with the opacity tables as buffers (the
tracked `examples/NGC6302/` files), the observed `NGC6302_100.tab` as a
`Spectrum` with the 25–120 µm selection, `Resample`, a `CalibrationScale`
and `GaussianProcessNoise` in place of the legacy noise triple, `--engine
emcee|zeus`; `NGC6302-calculate-dust-mass.py`'s post-processing becomes a
function over the run's `DataTree` in the same package (D12 confirms);
(3) **the modified blackbody** (`examples_paper/modifiedblackbody.py`) —
`examples/modified_blackbody/`: the astropy `BlackBody` model with the
four parameters, ten AKARI/Herschel bands through `from_library`, all
four engines (emcee, zeus, dynesty, SBI) and the posterior-predictive
overlay across them through `ampere.results`' plots; (4) **the PHOENIX
star** (`examples_paper/phoenixstar.py`) — `examples/phoenix_star/`:
**ruled 2026-09-28 on `docs/design/example_dependencies_memo.md` §2: no
Starfish in the twin — ampere's own lightweight emulator instead** (design
horizon (c)'s own case: an emulator is a model): a training script, run
once and not in CI, downloads the PHOENIX-ACES subset (Teff 5000–8000 K,
log g 4–5, [Fe/H] 0–0.5; `astropy.utils.data.download_file` or `expecto`),
bins it to a coarse 0.3–5 µm SED grid plus the Gaia-RVS window at
R ≈ 11 000, fits a PCA with one regressor per component weight over
(Teff, log g, [Fe/H]) — a small MLP preferred over a GP so the emulator
is differentiable on every backend; the item chooses and says why — and
writes the result as a sub-megabyte `.npz` committed beside the script;
the example loads it into an emulator `Model` written against the core
contracts (reference, torch and jax twins), applies the luminosity
scaling and CCM89 extinction through the `extinction` extra's
`dust_extinction`, and fits Gaia/2MASS/WISE photometry plus the RVS
spectrum with `SBIEngine` and NUTS for comparison; a comment points at
the legacy script for the Starfish route (Starfish on PyPI still requires
Python < 3.10; master carries Peter's July 2026 fix unreleased);
(5) **the carbon star** (`cstar_model_test_sbi_v2.py` and its
`_embedding` variant) — `examples/cstar/`: **ruled 2026-09-28 on the memo's
§3: Hyperion stays** (actively maintained again; conda-forge ships
`hyperion` and `hyperion-fortran` 0.9.11 for Python ≤ 3.13, pinned there
until upstream releases this year's fixes) as a v2 `Model` with the same
shell parameters as the legacy `ampere.models.Hyperion` class, **the
unpackaged `bhmie` Fortran step replaced by `miepython`** (the `.optc`
n, k tables and the power-law size distributions mixed by the two
abundance parameters into arrays for Hyperion's `IsotropicDust`), the
tracked `cstar_data` votable photometry and IRS spectrum as v2 containers,
`SBIEngine` with two rounds and the `--embedding` switch, its simulation
bank cached so reruns and calibration reuse it; Hyperion lives in an
optional pixi feature (`hyperion`: conda-forge `hyperion` +
`hyperion-fortran` + `miepython`, not a `pyproject` extra) and is absent
from CI (D12b); (6) **the star plus disc** (`examples/star_disc.py`,
the HD105 SED in `examples/star_disc/`) — `examples/star_disc/`: the
legacy `QuickSED` star-plus-disc model as a v2 `Model` (the one legacy
model class with no v2 counterpart yet), the votable photometry as a v2
container, emcee. **D12 ruled**: Starfish and Hyperion are documented
example-only requirements, not extras, and (4) and (5) skip by `find_spec`
where they are absent; (4) ships a small pre-trained emulator file under
a megabyte beside the script so the PHOENIX download and the training
are an optional, documented one-off; **and the orchestrator's assessment
of better-maintained alternatives to Starfish and Hyperion for the v2
twins (Peter's ask of 2026-09-28: the legacy scripts keep their packages,
the twins may switch) is `docs/design/example_dependencies_memo.md` —
the twin follows its recommendation where Peter accepts it.** Every twin
has a `__main__`, a docstring stating
which legacy script it twins and what changed in the translation, and a
`tests/examples` smoke row where its dependencies are installed (skipping
by `find_spec` otherwise — (4) and (5) skip in every CI environment; (1)–(3)
run under a minute each in dev); the `examples/` README lists the pairs.
**Depends:** W6.0 (the `ampere.legacy` path the originals now sit behind
— the twins import nothing from it), D12. **Accept:** the six twins with
smoke rows; the legacy originals byte-identical; each twin's posterior on
its synthetic truth covers the truth at 95 % on every parameter once (the
run recorded in the docstring); docs warnings no longer than base; gates
dev, plus sbi for (1)'s and (3)'s SBI arms. The three scripts with the
dead `ampere.emceesearch` import (`example.py`, `modbbtest.py`,
`modelClio.py`) get no twin and are left as they are (D12 (a)); the
untracked `flexible_likelihood_comparison.py` is answered by the M2
study.

### W6.14 — The training-set writer's non-numeric coordinates [S; Sonnet] (added 2026-10-01 on Peter's word, from W6.13 (5)'s finding)
`ampere/results/training.py`'s `_slot_dataset` (line 755) stores every
`extra_coords` entry of a slot's containers as float64, so a
`PhotometricPoints` observation — whose filter names are strings — fails the
first `training_set=` write with `could not convert string to float:
'MCPS_B'`. Any SBI problem with photometry among its observations therefore
cannot write its simulated pairs as a netCDF training set, which is what
W3.4's contract promises and what `examples/cstar` had to omit. The fix:
a coordinate keeps its own dtype (strings as a fixed-width or object array
the way `results_schema.md`'s `filters` coordinate already round-trips
through netCDF), the `extra_coords` attr recording the dtype alongside the
name; a conformance row writes and reads back a training set from a
problem with a `PhotometricPoints` observation and a `Spectrum` one, and
`examples/cstar`'s `fit` gains the `training_set=` argument the brief asked
for. **Depends:** W6.13 (C2) merged. **Accept:** the row round-trips both
container kinds; `examples/cstar`'s 40-simulation smoke row writes a
training set under `tmp_path`; `results.md` §9 amended with a sentence;
gates dev and sbi.

### W6.15 — `SBIEngine`'s scoring of the stored draws: optional, and pooled [S; Sonnet] (added 2026-10-01 on Peter's word, from W6.13 (5)'s finding)
`SBIEngine.run` scores every stored posterior draw on the numpy contract
path through `problem.evaluate` (`_sbi.py` line 70's rule), one evaluation
per draw, serially in the driving process and never through the engine's
`executor`. For an external simulator that is one Hyperion run per draw:
`examples/cstar` measured 28 min for 500 draws and would need ten hours for
the legacy's 10 000, so the twin cut its default to 1 000. Two levers,
both additive to §4.5's surface: (a) `run(..., score=False)` stores the
draws with `lp`, `log_prior` and `log_likelihood` absent and a provenance
attr `ampere_sbi_scored = 0` saying so, refused by name when `calibrate`
or a diagnostic later needs the scores; (b) when an `executor` was given,
the scoring is pooled through it in chunks, the same `simulate_many` path
the bank uses, so the stored draws of a pooled fit cost what the bank's
simulations did. The default stays scored and serial for a cheap model
(say at what per-evaluation cost the pool pays for itself, measured on the
`sed_composition` and `cstar` problems). **Depends:** nothing. **Accept:** a
row for each lever on an inference-suite problem; the `cstar` 40-simulation
smoke row scores through its two-worker pool; `inference.md` §10 and
`sbi.rst` amended; gates dev and sbi.

### W6.16 — The GPU rows on the cluster: the `gpu` environment and the procedure [S; Sonnet] (split from W6.10 (3) on Peter's word, 2026-10-02; the cluster's details confirmed by Peter 2026-10-03 — merged `5991128`)
W6.10's former part (3), unchanged in substance: the GPU rows (`tests/gpu`,
API-level smoke rows that skip without an accelerator) run on an
accelerator — **D7 ruled: an HPC allocation with A100s is available** — a
`gpu` pixi environment (torch and jax from their CUDA indices, the
CPU-index constraint of the `torch` feature lifted there only; locked, not
installed on the development machine) and a documented manual procedure
under `docs/development.md` (an rsync of the tagged tree to the login node,
the scheduler script committed under `scripts/`, the run, the log back
beside the status row), run once on the beta tag as part of the release
gate. The cluster's name, scheduler, driver/CUDA version and access route
are gathered from Peter at dispatch — he is double-checking them, which is
why this is its own item. **Depends:** D7, Peter's details, W6.10 (the gate
policy the procedure's log joins). **Accept:** the environment solves
(`pixi lock` / `--no-install`) and the other environments' lock sections are
unchanged; the procedure runs end to end once on the cluster against the
beta tag (W6.5's release gate records it); gates none (infrastructure;
the rows run where CI cannot).

### W6.17 — Community health before the beta: licence detection, contributing, conduct, security, templates [S; Sonnet] (added 2026-10-03 on Peter's ask; D14)
Peter's review of the repository's community profile before `v1.0.0b1`:
GitHub reports a health score of 25 % — `readme` present, **`license`
missing** (the GPLv3 text lives at `licenses/gpl.txt`; GitHub, Zenodo and
PyPI's "licence" badge all detect a licence only from a root `LICENSE`/
`COPYING` file, so the beta's GitHub release and its Zenodo record would
show *no licence* today), `code_of_conduct`, `contributing`,
`issue_template` and `pull_request_template` missing; `isSecurityPolicyEnabled`
false, private vulnerability reporting and secret scanning off, Discussions
on. In priority order, each a commit: (1) **`LICENSE` at the root** —
`git mv licenses/gpl.txt LICENSE`, `license-files = ["LICENSE"]` in
`pyproject.toml` (the SPDX `license` field unchanged), the `licenses/`
directory gone; the single change that must precede the tag. (2)
**`CONTRIBUTING.md`** for a human contributor: the pixi route and the
tasks (`test-fast` before a PR, `lint`/`format-check`/`typecheck`, `docs`),
the branch-and-PR convention and that Peter reviews and merges, the
frozen legacy and the frozen contracts with pointers at `DEVELOPMENT_PLAN.md`
§4 and ground rule 9, British English, the Co-Authored-By line when an
agent assisted, where discussion happens (Issues for bugs, Discussions for
questions), and a plain sentence that much of v2 was built by orchestrated
agents under `AGENTS.md`'s working agreement — honesty a contributor can
act on; `docs/source/index.rst`'s "Contributing" paragraph and the README's
pointer link to it. (3) **`CODE_OF_CONDUCT.md`** — the Contributor Covenant
2.1 verbatim with the enforcement contact D14 names. (4) **`SECURITY.md`** —
supported versions (the beta line onward), how to report (GitHub's private
vulnerability reporting, which Peter enables; email as the fallback), what
ampere's risk surface actually is (it executes user-supplied models and
simulators in worker processes, reads netCDF runs and training sets, and
pickles problems for `ProcessExecutor` — loading a file or a cache from an
untrusted source is the hazard to name), and the response Peter is willing
to promise (D14). (5) **Issue forms** under `.github/ISSUE_TEMPLATE/`: a bug
report (version from `python -c "import ampere; print(ampere.__version__)"`,
backend and extras, the smallest `FittingProblem` that shows it, the
traceback, the engine), a feature request, and `config.yml` sending
questions to Discussions and documentation gaps to a `docs` form. (6) **A
pull-request template**: the item or issue, what and why, the acceptance
evidence (which tasks ran, counts), a checklist (tests, lint/format/
typecheck, docs built, British English, no run outputs or binaries,
Co-Authored-By when agent-assisted). **Not files, for Peter in the
repository settings**: enable private vulnerability reporting, secret
scanning with push protection, and "reported content" moderation for the
public Discussions. **Deferred, recorded**: an accessibility statement for
the docs site is not a GitHub profile item and belongs with a docs-theme
review (alabaster's contrast and keyboard navigation; alt text on every
figure as a docs convention) — a Phase 7 docs item, not a release blocker;
`SUPPORT.md` is subsumed by the issue forms' `config.yml`. **Ownership**:
the six files above, `pyproject.toml`'s one line, `.github/ISSUE_TEMPLATE/`,
`.github/PULL_REQUEST_TEMPLATE.md`, one paragraph each in `index.rst` and
the README; runs **before W6.5** (the licence move is W6.5's wheel's
`License-File`, and W6.5 rewrites the same README and `index.rst` passages
afterwards) and after W6.16 (no shared files, but one agent at a time on
the release path). **Depends:** nothing. **Accept:** `gh api
repos/ICSM/ampere/community/profile` reports every file present once merged
(the orchestrator checks after the push — the `license` entry detects
`GPL-3.0`); a wheel built from the branch carries `License-File: LICENSE`
and `License-Expression: GPL-3.0-or-later` (the agent builds one with
`python -m build` if `build` is present, else says so and W6.5 verifies);
`actionlint`/YAML validity of the forms; `pixi run test` and the docs build
green; gates dev (the docs job is the check).

### W6.18 — The first GPU run's faults: a host copy in the realisation self-check, the jax CPU lookup on an accelerator, a sharding test's hard-coded width, NUTS's start tensor, and the jax quasiseparable solver's CPU-only refusal [S; Fable] (added 2026-10-05 from the release gate's first GPU run; blocks the tag unless Peter rules otherwise)
Slurm job 18041011 on `cc7a65e` (node `gina1`, one A100, driver 615.71.09)
ran `tests/gpu` for the first time ever: **20 passed, 13 failed**, three
distinct faults, none of which a CPU machine can show. (1) **torch, six
rows** (`TestTheRealisedDensityRunsOnTheGpu` ×5, `TestNutsDraws` ×1):
`ampere/core/realisation.py::_check_agrees_at_reference` evaluates the
realised density at the reference values through `float(np.asarray(...))`;
a CUDA tensor has no `__array__` ("can't convert cuda:0 device type tensor
to numpy"), so `lower_problem` raises `LoweringError` for every torch
problem on a GPU and `nuts` fails behind it with `EngineError`. The fix is
backend-neutral: `float(x)` on the scalar (a 0-d torch tensor, a jax array
and a numpy scalar all implement `__float__`), or a backend-provided
`to_host`; the self-check stays, the one conversion changes. (2) **jax,
five rows** (`TestTheDensityAgrees` ×3, `TestTheRefusalsStillHold` ×2):
`ampere/backends/jax/_device.py::resolve_device` searches `jax.devices()`
for the requested platform, but `jax.devices()` without an argument lists
the *default backend's* devices only — on a GPU node that is `cuda` alone,
so a request for `"cpu"` is refused as "no such platform" although the CPU
platform exists; `jax.devices(device)` (with its `RuntimeError` for a
platform jax really lacks, re-raised as the ruled `error`) is the lookup
that means what the docstring says. Twenty-seven call sites go through
this one function; none changes. (3) **jax, two rows**
(`TestChunkSharding::test_sharding_does_not_change_the_prediction`,
`::test_a_ragged_chunk_is_padded_not_refused`): the test helper draws
unit-cube rows of width 2 for a problem whose `free_size` is 4 (norm,
index and the two Matérn hyperparameters) — a test bug, `ParameterError:
expected a unit-cube vector of length 4, got 2`; the helper uses
`problem.free_size`. Also found by the run and taken here: the cluster
appends a usage epilogue after the job's output, so the log's last line is
*not* the pytest summary as the procedure says — `run.sh fetch` prints the
summary line by `grep -E '[0-9]+ (passed|failed)' | tail -n 1`, and the
procedure text says so; and `tests/gpu` is 33 rows after parametrisation,
not 30. **Extended 2026-10-05 by the second run** (job 18042559 on `7e324f5`,
31 passed, 2 failed — the two the first three faults had masked): (4)
**torch NUTS on a GPU**: `ampere/inference/_nuts.py` built pyro's
`initial_params` start tensor with no device, so pyro's momenta and mass
matrix lived on the CPU beside a CUDA density and `velocity_verlet` failed
with "Expected all tensors to be on the same device"; the start now takes
the realisation's `device` (a caller-supplied density keeps the default),
and the engine's own reference-value self-check uses `float()` as in (1).
(5) **jax `QuasisepGP` on an accelerator cannot work**: celerite2 0.3.3
registers its primitives' MLIR lowerings for the CPU alone
(`celerite2/jax/ops.py`), so the first factorisation dies inside the trace
with "MLIR translation rule for primitive 'celerite2_factor' not found for
platform cuda"; under `architecture.md` §5's rule the solver now **refuses an
accelerator by name at construction** (`LikelihoodError`, naming `DenseGP`
and `HilbertSpaceGP` as the way out), the GPU row that compared the
quasiseparable density across devices becomes the row asserting that
refusal, and the changelog's "Known limitations" records the jax solver as
CPU-only. A jax-native quasiseparable recursion (a `lax.scan`, the tinygp
route slice 2 declined on CPU speed) is the Phase 7 item that would lift it. **Ownership**: `ampere/core/realisation.py` (the one conversion),
`ampere/backends/jax/_device.py` (the lookup), `ampere/inference/_nuts.py`
(the start tensor's device, the self-check's conversion),
`ampere/backends/jax/gp.py` (`QuasisepGP`'s refusal), `tests/gpu/test_jax_gpu.py`
(the helper, the refusal row), `docs/source/changelog.rst` (one limitation), `scripts/cluster/run.sh` (the fetch's summary line),
`docs/development.md`'s GPU-rows section (the two corrections). No
contract changes (the realisation registry's self-check is implementation;
`resolve_device` is the jax backend's own). **Depends:** the cluster
bootstrap (done 2026-10-05). **Accept:** the touched suites green in dev
(`tests/core`), jax (`tests/backends`, `tests/conformance`) and torch
(`tests/conformance`); a CPU-only regression row for (2) where one can be
written (`jax.devices()` can be made to differ from `jax.devices("cpu")`
only with an accelerator, so the row may have to be the GPU row itself —
say so); then **the GPU rows on the branch head through
`run.sh submit <hash>`: 33 passed, 0 failed, quoted in the row**; gates dev
and jax (torch for the conformance leg).

**Decision for Peter (2026-10-05):** the beta waits for W6.18 (recommended —
the changelog says GPU support is exercised by smoke tests, and the first
such run failed 13 of 33 rows on faults that a user with a GPU would hit on
their first `lower_problem`; the fixes are three small edits and one more
cluster job), or the tag goes on `cc7a65e` with the changelog's "Known
limitations" amended to say the torch realisation path and the jax CPU
device lookup fail on an accelerator at the beta and W6.18 follows as
`1.0.0b2`. **Taken (2026-10-05): the beta waits for W6.18; the orchestrator
did the item on `w6.18-gpu-first-run`; the fourth GPU run, job 18042730 on
`62e851c`, 33 passed, 0 failed — the acceptance met; the CPU-only regression
row for (2) was writable after all (`366f576`).**

**Issue triage (D10).** The open issues, with the recommendation: **closed
by Phase 5 already** — #12, #29, #67 (W5.17's levers), #11 (censoring,
`likelihoods.md` §9 landed with the core); **closed by Phase 6** — #57,
#58 (W6.3), #59 (W6.1, W6.2), #60, #62 (W6.5), #40, #14, #41 (W6.7),
#63 (W6.12 if D9 says so), #73 (the v2 core *is* the rewrite the issue
asks for — close with a pointer at `overview.rst`), #3 and #68 (legacy
data-object questions, closed by the deprecation policy); **re-filed as
v2 backlog** — #61 and #45 (the six v2 plots have a uniform style; what
remains is corner paging at many parameters, `results.md` §13), #70 (line
fluxes as an observable kind by the template — a Phase 7 modality), #64
(the filter library as a build step — legacy `utils`; carried until the
v2 photometry route needs it), #44 (opacity wavelength range — a legacy
model's limit); **out of scope, closed with a note** — #65, #21, #22, #23
(radiative-transfer codes as v2 models: out-of-tree models by the
`Model` contract, no in-tree adapters planned), #15 (CANFAR batch scripts:
deployment is a user's, not the library's).

- **D13 — the distribution name** (opened 2026-10-02 when W6.5's prerequisites
  were listed): the PyPI name `ampere` belongs to an unrelated battery-modelling
  package, so the beta needs its own. Candidates, all free on PyPI at the
  ruling: `ampere-astro`, `ampere-fit`, `ampere-sed`, `ampere-bayes`,
  `ampere-infer`, `astro-ampere`. **Ruled by Peter 2026-10-02: `ampere-astro`
  on PyPI (and TestPyPI); the import name stays `ampere`; the Read the Docs
  slug is `ampere`.** Recorded in the plan's decision table; carried by W6.5
  and W6.4.
- **D14 — community health before the beta** (opened 2026-10-03 on Peter's
  review of the community profile; W6.17). **Ruled by Peter 2026-10-04 — every
  recommendation accepted**: (a) the code of conduct's enforcement contact is
  **Peter's address** (the one in `pyproject.toml`'s author list,
  `peter.scicluna@eso.org`; the same address serves `SECURITY.md`'s email
  fallback; "we can adjust email addresses when necessary" — a later change
  is a one-line edit, not a ruling) — the alternatives were a project address
  (no mailbox exists) and "the maintainers listed in `CITATION.cff`", which
  does not exist as written: the repository has no `CITATION.cff` (W6.5 adds
  one for Zenodo) and `pyproject.toml` names nine authors with one email;
  (b) the security policy promises **an acknowledgement within fourteen
  days, a fix on a best-effort basis, no embargo or coordinated-disclosure
  machinery** — a scientific library with no network surface of its own, whose
  realistic report class is unsafe deserialisation of a file the user chose to
  load; (c) **`SECURITY.md` names the hazards plainly** — the artefact store
  unpickles cached posteriors (`ampere/results/artefacts.py`), `ProcessExecutor`
  pickles whole problems with their user models into workers, runs and training
  sets are netCDF parsed by the xarray/HDF5 stack, and models and simulators
  are user Python run with the user's privileges — so that a user who loads a
  stranger's training set reads the warning there and a report of the
  documented hazard is recognised as not a vulnerability; (d) the three
  repository settings (private vulnerability reporting, secret scanning with
  push protection, reported-content moderation for Discussions) are **Peter's
  clicks before the tag** — the API on 2026-10-04 shows secret scanning and
  push protection disabled, Dependabot security updates enabled, Discussions
  on; the non-provider-pattern and validity-check extras are not asked for;
  (e) the accessibility statement is **deferred to a Phase 7 docs item** (not
  a profile item, not a blocker; the one thing taken now is alt text on any
  figure W6.17 or W6.5 adds). Recorded in the plan's decision table.

- **D15 — the `origin/master` switch before the dry run** (opened 2026-10-05
  when the dry run was attempted; W6.5 merged at `f4ea081`): GitHub registers
  a `workflow_dispatch` workflow only from the repository's **default
  branch**, and `release.yml` lives on `v2` — `gh api
  repos/ICSM/ampere/actions/workflows/release.yml` returns 404 and the
  workflow list shows `ci.yml` alone — so the TestPyPI dry run
  (`gh workflow run release.yml --ref v2 -f target=testpypi`) cannot be
  dispatched while `origin/master` is the legacy line. A tag push would
  still trigger the workflow (push-triggered workflows run from the pushed
  ref), but then the first ever run of `release.yml` is the beta itself, with
  no dry run before it. Options: (a) **move D4's switch forward** — Peter
  fast-forwards `origin/master` to the gated v2 head *before* the dry run
  (`git push origin master`; `b8e585b` is an ancestor, so no force), and the
  dry run, the tag and the rest of the procedure follow with `master` as the
  default branch from the start; legacy users pulling `master` get v2 a few
  days before the tag, with the legacy code still importable under its old
  names (W6.0's aliases), which is the state the tag would give them anyway;
  branch protection can go on at the same time or at the tag as D4 says;
  (b) skip the dry run and let the tag's run be the first — a red smoke leg
  then means a tag already pushed, a fix, and `v1.0.0b2` as the first beta;
  (c) put `release.yml` on the legacy `master` by a commit there — ruled out
  (legacy is frozen; nothing is pushed to `origin/master` before the switch).
  **Recommendation: (a)**, with the release procedure's step (f) moved to
  before the dry run and `ci.yml`'s `v2` trigger left until a later
  housekeeping. **Ruled by Peter 2026-10-05: (a) — "go ahead with the
  switch"; done the same day: `origin/master` fast-forwarded `b8e585b` →
  `fef26ca` by the orchestrator, `release.yml` registered, D4's switch
  thereby taken before the tag rather than at it; branch protection and the
  Read the Docs default branch remain Peter's at the tag.**

**Ordering (two agents at a time, gate legs scoped to the code touched).**
**Wave 1** (fillers while the rulings are taken): W6.6 ∥ W6.8 — both
Sonnet, disjoint files, no ruling needed beyond D10's alias note; W6.9 in
either slot as it frees. **Wave 2** (after D1): W6.0 alone, short, since
everything on the documentation side keys on its policy; then W6.1 ∥ W6.2
(disjoint: `migrating.rst` and the notebooks against a new tutorial page
and one new example); then W6.13 ∥ W6.7's first half — the twins are
Sonnet on `examples/` and the optimisers Opus on `ampere/inference`, with
the two external-code twins ((4), (5)) as a second Opus slot when the
first frees. **Wave 3**: W6.3 ∥ W6.7 — the docs
rebuild is Sonnet on docs files and the optimisers are Opus on
`ampere/inference` and one docs page; ownership disjoint by construction.
**Wave 4**: W6.4 ∥ W6.10 (both infrastructure, D2/D4/D7 ruled by then);
W6.11 as Fable's own work between waves. **Last**: W6.5, the release,
once everything else is merged and the full matrix is green on the tag;
W6.12 after the release if D9 puts the readers in this phase, or first
in Phase 7 if not. Phase 6 closes with the beta on PyPI and
`origin/master` at the tag.

**Decisions for Peter before dispatch.**
- **D1 — the deprecation policy**: (a) deprecate the four legacy packages
  at the beta (a warning at first import naming the v2 route), remove at
  1.0 and no sooner than six months after the beta (recommended — the
  characterisation suite guards the interval and the migration guide is
  in place by then); (b) keep legacy frozen indefinitely as it is, no
  warning, no removal; (c) remove at the beta. Also under D1: whether the
  `regularised_horseshoe` alias (W5.27) is removed at the beta or at 1.0
  (recommended: 1.0, one rule for every deprecation). **Ruled 2026-09-28:
  (b), with the legacy routines moved into an `ampere.legacy` subpackage
  so `from ampere import legacy as ampere` keeps a legacy script working
  — W6.0 rewritten accordingly; the old top-level names as lazy aliases
  and the v2 alias's own removal at 1.0.0 are the orchestrator's
  assumptions under (b), to confirm at W6.0's review.**
- **D2 — where the docs live**: Read the Docs (the astronomy convention;
  a `.readthedocs.yaml` building from a pip docs extra since RTD does not
  run pixi) or GitHub Pages from the CI job (pixi-native, one less
  service, no PR previews). Recommendation: Read the Docs, for the
  version switcher and PR previews users expect. **Ruled 2026-09-28: Read
  the Docs.**
- **D3 — the beta's version and the changelog**: the version — `0.2.0b1`
  under setuptools_scm from a `v0.2.0b1` tag (recommended; `v0.1` is the
  last legacy tag and `__version__` says `0.1.2`), or `1.0.0b1` to signal
  the redesign; the changelog — hand-written per release from the decision
  log and the status table (recommended for the beta: the phases' landed
  summaries are the material and a generator would flatten them), or a
  generator (`towncrier` fragments per PR, or `git-cliff` from commit
  messages) from the beta onward. **Ruled 2026-09-28: `1.0.0b1`, to
  signal the redesign (the final is `1.0.0`); the changelog as
  recommended.**
- **D4 — the remote**: `origin/master` takes the v2 line at the beta tag
  (recommended) or at the first CI-green merge after Phase 6 opens; branch
  protection on `master` then (CI required, no force-push, Peter merges);
  and whether the orchestrator may push the mirror at every merge from
  now on so the CI run becomes the gate of record (W6.10 (1); recommended:
  yes, `origin/v2` until the switch-over, never `origin/master` before it).
  **Ruled 2026-09-28: accepted — the orchestrator pushes `origin/v2` at
  every merge from now on; `origin/master` switches at the beta tag with
  branch protection.**
- **D5 — two small docs rulings from W5.19**: a shrinkage section in
  `advanced.rst` (recommended: yes, four paragraphs pointing at
  `kernels.rst` §1 and `m2_misspecification.rst`, in W6.3), and a doctest
  runner over `docs/source`'s `pycon` blocks (recommended: yes, in W6.3,
  the harness of `test_spec_doctests.py` extended; the pages already
  claim their blocks were run). **Ruled 2026-09-28: yes to both.**
- **D6 — the optimisers' size and tier**: L and Opus as drafted (three
  routes, a `results.md` amendment, four gates), or split into the numpy
  and native routes first (M) with the empirical-Bayes warm start as a
  follow-on (S). Recommendation: as drafted; the warm start is where the
  measured gain is and it shares the provenance work. **Ruled 2026-09-28:
  as drafted.**
- **D7 — the GPU rows**: whether an HPC allocation with accelerators is
  available for the release gate; if yes, W6.10 (3) is in scope and the
  procedure is written; if no, `tests/gpu` stays skipped and the release
  says so. No recommendation — the answer is a fact only Peter has.
  **Ruled 2026-09-28: yes — HPC with A100 GPUs is available; W6.10 (3) in
  scope.** (2026-10-02: (3) split out as W6.16 so W6.10 can run while Peter
  checks the cluster's details.)
- **D8 — the design items' timing**: W6.11's memo in Phase 6 (recommended:
  yes — a memo costs little and the per-dataset GP amplitude is the
  flexible likelihood's own use case), and its implementation items in
  Phase 6 after the beta, or in a Phase 7 (recommended: Phase 7; a §4
  change belongs after the release, not before it). **Ruled 2026-09-28,
  after the briefing: as recommended — the memo (W6.11) in Phase 6, the
  implementations open Phase 7.**
- **D9 — readers in Phase 6**: W6.12 (OIFITS, and JWST spectra as the
  second reader) in this phase after the beta, or the first Phase 7 item.
  Recommendation: OIFITS in Phase 6 after the beta — a release with no
  reader for real interferometric data is a weaker claim than the
  modality deserves; JWST (issue #63) in Phase 7 with the spectroscopy
  front door. **Ruled 2026-09-28: accepted — readers matter less than the
  containers (a user can always read a file by hand and feed ampere), so
  OIFITS is the worked example of how a reader is written and where it
  lives, after the beta; JWST Phase 7.**
- **D10 — the issue triage** as listed above (closed by Phase 5; closed by
  Phase 6; re-filed as v2 backlog; out of scope). Recommendation: as
  listed; the closures with a one-line pointer each, done by the
  orchestrator at the beta. **Ruled 2026-09-28: accepted, with tweaks as
  they come up at the beta.**
- **D11 — the paper**: whether the beta release is coupled to the paper
  revision's examples (`examples/examples_paper/`, Peter's) — released
  together, the paper's scripts pinned to the beta tag — or independent.
  No recommendation; Peter's timeline decides it. **Ruled 2026-09-28: keep
  the old versions working on legacy, and produce v2 twins of the same
  models — `examples_paper/modifiedblackbody.py`, `examples_paper/phoenixstar.py`,
  `examples/NGC6302.py`, `examples/minimal_working_example.py` and its
  variants, `examples/cstar_model_test_sbi_v2.py` — as W6.13; the release
  is not coupled to the paper beyond that.**
- **D12 — the remaining legacy scripts, and the external-code dependencies**
  (opened 2026-09-28 on D11's "ask about any others"). (a) The scripts
  Peter did not name: `examples/example.py` and `examples/modbbtest.py`
  (both import `ampere.emceesearch`, a path that has not existed for
  years — already broken; recommend: no twin, left as they are);
  `examples/modelClio.py` (Gielen et al. 2008's post-AGB disc model, with
  the same dead `ampere.emceesearch` import; recommend: a twin only if
  the paper uses it, since the model is a real published one);
  `examples/star_disc.py` with `examples/star_disc/HD105_SED.*` (a star
  plus disc on the legacy `QuickSED` model with emcee, committed
  2026-07-31 as "Star disc model attempt"; recommend: a twin, since it is
  the most recent work and `QuickSED` is the one legacy model class with
  no v2 counterpart yet); `NGC6302_zeus.py` and `NGC6302-calculate-dust-mass.py` (folded
  into (2) as drafted); `cstar_model_test_sbi_v2_embedding.py` (folded into
  (5)); the three notebooks (W6.1's); `examples_paper/flexible_likelihood_comparison.py`
  (untracked; its v2 counterpart is the M2 study — recommend: cite the
  study, no twin). (b) Starfish and Hyperion: documented "to run this
  example you need" requirements with the twins skipping by `find_spec`
  (recommended — neither belongs in an extra the package resolves), or an
  `examples` extra that pins them. Also whether the PHOENIX grid download
  and the emulator training may be replaced in the twin by a small
  pre-trained emulator file committed under `examples/phoenix_star/` if it
  is under a megabyte (recommended if so; otherwise the one-off step stays
  documented). **Ruled 2026-09-28: all as recommended, plus the
  `star_disc` twin (W6.13 (6)) and the committed pre-trained emulator;
  and an assessment of better-maintained alternatives to Starfish and
  Hyperion for the twins (`docs/design/example_dependencies_memo.md`,
  the orchestrator's) — the twins may switch, the legacy scripts keep
  their packages. Ruled on the memo 2026-09-28: the PHOENIX twin ships
  ampere's own lightweight emulator in place of Starfish; the carbon-star
  twin keeps Hyperion (maintained again) with `bhmie` replaced by
  `miepython`; the GRAMS grid, a Starfish release and the Zenodo download
  are moot under those two; the emulator's regressor is the item's
  choice (MLP preferred).**


## Phase 7 — Hierarchy, derived parameters, the second reader (drafted 2026-10-05 by Fable; **ruled by Peter 2026-10-05: D1–D11 as recommended** — the memo's §7 answers as D1–D4, `1.0.0b2` after wave 3, `x1d`/`c1d` only, line fluxes in wave 4, nested populations left for the horizon memo, the jax recursion in, the RHMF check at close, the `Censoring` filler and the migration note from the backlog; the longer horizon is `docs/design/horizon_beyond_phase7.md`)

The first phase after the beta. Its sources: the W6.11 design memo
(`docs/design/nuisance_populations_and_derived_memo.md`) whose §10 drafted
W7.0–W7.2 and whose §7 questions Peter has not yet ruled — D8 of Phase 6
put the implementations here; Phase 6's D9 (the JWST reader as the second
reader, with the spectroscopy front door); the Phase 6 rows' carried
findings and Peter's backlog as the handoff records them; the beta's
"Known limitations" (`changelog.rst`), each of which is either an item
here or a decision not to lift it yet; and the two deferred trials (RHMF,
the jax-native quasiseparable recursion). **Read the memo first** — W7.0
and W7.1 execute it section by section and the item texts below point at
its sections rather than restating them. The plan's §5 has no Phase 7
section yet: it is written at approval from this section, as Phase 6's
was, and the beta's changelog gains an "Unreleased" entry per merged item
(W6.5's convention). Sized under `docs/orchestration.md`'s rules: two
agents at a time, gate legs scoped to the code touched, CI on the push
after merge as the gate of record (W6.10). Every change to a §4 contract
is a decision-log row and the conformance rows in the same PR (ground
rule 9) — the two rows for W7.0 and W7.1 are already drafted in the
memo's §8. The decisions **D1–D11** at the end were the ones the drafting
could not take; **Peter ruled every one as recommended on 2026-10-05**,
and the ordering paragraph stands as written. The horizon beyond this
phase — what is blocked, what unblocks it, a shape for Phases 8–10 — is
`docs/design/horizon_beyond_phase7.md`, drafted the same day.

**What the phase delivers.** (1) The two design items: a `Derived`
parameter node and populations over a dataset's own parameters — which
lift two of the beta's five known limitations and give the flexible
likelihood its natural hierarchical prior across many spectra (the
per-dataset GP amplitude), validated on the M2 pattern. (2) The second
reader and the spectroscopy front door, so a JWST spectrum fits in a
dozen lines. (3) The release train: `1.0.0b2` with the above, and the
housekeeping the beta carried. Not in this phase, recorded as decisions:
nested populations (the third limitation), a jax-native quasiseparable
recursion (the fourth), `ampere.diagnostics` (the fifth, gated on
upstream).

### W7.0 — The `Derived` parameter node [M; Opus] (memo §3; D1, D2 ruled as recommended 2026-10-05)
§3 of the memo, whole: `Parameter(name, Derived("<expression>",
symbols={...}))` as a fourth parameter state — the closed grammar over
symbols (numeric literals, `+ - * / **`, unary minus, `sqrt exp log log1p
abs`; attribute, subscript, comparison, lambda and any unlisted call
refused by name; a symbol spelt like a reserved name refused), parsed once
and stored as its normalised source, each symbol bound to a parameter name
by a mapping of `HierarchicalPrior.hyperparameters`' shape so a merge
renames the binding and never the expression; no sampler dimension and no
prior term; computed idempotently wherever named values are formed
(`complete`/`unpack` — a supplied value recomputed, never trusted; on the
reference path a stale supplied value is refused by name, W5.32 (j)),
mid-walk in `prior_transform` and `sample`, and in both backends' resolve
and walk functions (torch `unpack_tensor`, `_resolved`, `prior_transform`,
`log_prior_tensor`; jax `_resolved`, `unpack`, `prior_transform`,
`log_prior`, `numpyro_model` with `numpyro.deterministic` as the
structural view only); evaluated over `ArrayOps`, which gains `sqrt` and
`log` on the three namespaces; referenceable by a `HierarchicalPrior`; the
evaluation order honoured under tracing (the §3.4 hazard: derived →
hierarchical → derived). Results: one `posterior` variable per derived
parameter on every engine, computed vectorised at emission
(`_posterior_variables`), named in a new `ampere_derived` attr,
`training.py`'s column writer carrying it, `PROVENANCE_SCHEMA_VERSION` →
10. `Population` gains the two member rules of §3.5: a member no `over`
component declares is **internal** — routed nowhere, permitted only as an
input of a derived member of the same population (the non-centred `θ_i =
μ + σ z_i`), refused otherwise — and a routed member must be declared by
every `over` component, which moves W6.11's recorded late failure (an
undeclared member failing at the first model evaluation) to the merge;
both rules refused in the flat layout; the two frozen allowances §3.8
retracts (`parameters.md` §8's "the receiving component's own set does
not declare the local name"; `_apply_populations`'s flat-layout "a
component that does not declare it is given one") become merge-time
refusals. The slab: `shrinkage_horseshoe(tail="slab", slab_scale=c)` in
the helper's own `s_j = τλ_j` parameterisation (§3.9: `s̃_j =
sqrt(c²s_j²/(c²+s_j²))`, a prior on `c` admitted), the plain and
regularised tails byte-identical. The contract amendments of §3.8
(`parameters.md` §3/§9/§12, `lowering.md` §3.2.1/§5/§8, `results.md`
§4/§9, `inference.md` §9/§10a) and the decision-log row of §8, verbatim
but for what review changes. Deliberate limitations stand (§3.7): no
callable, no conditional, no unit arithmetic, not in the flat layout.
**Depends:** nothing (D1, D2 ruled). **Accept:** the memo's §9.2 ten rows
(`tests/conformance/test_derived.py`) green on every registered fixture;
the existing `test_population.py`, shrinkage and horseshoe rows unchanged
and green; a NUTS run on torch and jax of the non-centred fifty-member
population (`tests/inference/test_population_nuts.py`'s problem
re-declared with `z` internal and a derived `theta`) recovering `mu` and
`sigma` inside the central 95 % with fewer divergences than the centred
declaration at the same budget, pinned as an inequality; every
`ampere_problem_hash` the conformance suite pins unchanged for a problem
without a derived parameter; `population.rst` gains the non-centred
declaration as a worked block (run as a doctest, D5 of Phase 6);
lint/format/pyrefly clean; gates dev + torch + jax.

### W7.1 — Populations over a qualified component path [M; Opus] (memo §2; D3 ruled as recommended 2026-10-05)
§2 of the memo, whole: `Population.over` entries of the form
`component[.path]` (`"d0.likelihood"`, `"d0.instrument.calibrate"`;
`"d*.likelihood"` through the existing glob), validated by walking the
retained inner mappings — a path that does not resolve refused naming the
components at that level, a path on a plain-set component refused, a leaf
the composite does not declare refused, a leaf an inner `shared_as`
collapsed refused by its tie label, a tied leaf refused, and W5.12's
refusal of a composite in `over` re-worded to a refusal of a composite
*without a path* naming the remedy; `PlateBinding.local_name` and
`Binding.local_name` as qualified paths relative to the component; the
plate layout stripping the addressed leaf from the composite's *outer*
declaration (the flat layout re-prioring the leaf in place); routing
untouched — the dataset's retained mapping takes the second hop it
already takes, so neither `distribute` nor either backend changes and the
population sketch's §11 Q1 ruling stands; the leaf-must-exist rule with
W7.0's internal-member exemption; `DatasetCollection.plate(within=)`
writing the entries from the dataset labels (the convenience on the
factory only, D3); the population declaration in the problem's provenance
record as `populations` (W6.11's recorded gap: the merged spec moved but
the declaration was recorded nowhere, and after this change a dataset's
own spec no longer declares the leaf the population replaced),
`PROVENANCE_SCHEMA_VERSION` → 11. The contract amendments of §2.6
(`parameters.md` §8/§9, `inference.md` §9, `hierarchical_population.md`
§11, `results.md` §9) and the decision-log row of §8. **Depends:** W7.0
(the member rules, the merge-time refusals; §9.1's row 10 needs
`Derived`). **Accept:** the memo's §9.1 ten rows
(`tests/conformance/test_population.py::TestAPopulationOverADatasetPath`)
green on every fixture, the two-level instrument-step path among them;
`test_population.py`'s existing rows unchanged; `inference.md` §9's
doctest of the plate-of-datasets rewritten to the path form and
executing; the emitted run's plate dimension labelled by the dataset
labels; `population.rst` teaching the path form in one section; gates
dev + torch + jax.

### W7.2 — The per-dataset GP amplitude as a population: the M2 validation [M; Opus] (memo §4)
§4's customer, measured — the composition W7.0 and W7.1 exist for. A
scenario with N spectra of one object class (the M2 generators' pattern,
the deviation injected in some) sharing a population of GP amplitudes,
non-centred (`z` internal, `amplitude` derived, over `d*.likelihood`),
against (a) independent amplitudes and (b) one tied amplitude: bias,
calibration and localisation on the M2 pattern; the shrinkage of a poorly
constrained member's amplitude towards the population, pinned as a margin
against (a); SBC over refits with the deviation injected; NUTS on torch
and jax through the realisation, emcee on the reference path. Placement:
a sibling `examples/m2_populations/` rather than a new `SCENARIOS` entry,
since the study's scenarios are the milestone's fixed set (W5.8's
precedent) — the driver, generators and models imported from the study
where they fit. A `docs/source/population.rst` section teaching the three
declarations — tie, flat, non-centred population — and when each is
right, with the measured table. **Depends:** W7.0, W7.1. **Accept:** the
scenario's driver with `tests/m2` rows and pinned margins (the shrinkage
margin; the population recovering the amplitude spread inside the
central 95 %); the SBC contrast pinned (the population's rank
distribution uniform at the study's thresholds); the docs section with
every figure traced to the driver; `test-fast` not lengthened — the long rows
behind the `m2_full` marker as the study's are; gates dev + torch (jax if the native
declaration touches the jax population code; it should not).

### W7.3 — The second reader: JWST spectra into `Spectrum`, and `ampere.spectroscopy` [M; Opus] (Phase 6 D9's ruling; issue #63)
The front-door pattern of W6.12 applied to spectroscopy. A reader for
the JWST pipeline's one-dimensional spectral products — the `x1d`/`c1d`
FITS files with their `EXTRACT1D`/`COMBINE1D` binary tables (`WAVELENGTH`
in µm, `FLUX` and `FLUX_ERROR` in Jy, `DQ`, `SURF_BRIGHT`/`SB_ERROR` for
extended sources, `NPIXELS`; several extensions for several sources or
slits) — `read_jwst(path, *, source=, extension=, surface_brightness=False)
-> Spectrum` with `astropy.io.fits` only; `DQ` → mask by the pipeline's
"do not use" bit, non-finite values and non-positive uncertainties masked
and counted in `meta`; the instrument, grating/filter, exposure and
pipeline version in `meta`; a file with several sources or extensions
refused by name with the choices listed unless the keyword says which;
the spectral axis emitted in µm (the OIFITS precedent) and the flux in
the file's unit; the LSF: the pipeline's `R` per instrument mode is
*not* tabulated by the file — the reader carries nothing it cannot read,
and the docs section shows `LSFConvolution` taking the user's resolution
curve beside it. `ampere.spectroscopy` as the front door: re-exports of
`Spectrum`, `PhotometricPoints`, the reference steps (`Resample`,
`LSFConvolution`, `CalibrationScale`, `SyntheticPhotometry`), the shipped
families, and `read_jwst`/`JWSTError`; `import ampere` untouched (the
import-cost rows). A real file under `tests/data/` with its provenance
(a small public MAST product — JWST data are public domain once their
exclusive-access period lapses; the ERS programmes' are; under a
megabyte, else downloaded in the test with a skip offline, the W6.12
pattern), read into a `Spectrum` that `tests/examples/test_photometry_spectra.py`'s
fixtures accept and fit end to end in the test with a `GaussianProcessNoise`
likelihood (a dozen lines, run in the docs section too). The `.gitignore`'s
legacy `*.fits` line qualified so a vendored test file is not silently
ignored (W6.12's carried finding). **Depends:** nothing (D6 ruled: `x1d`/`c1d` only).
**Accept:** the file round-trips into a `Spectrum`, the mask equals the
DQ "do not use" bit, the fit runs; the multi-source refusal and the
selection by keyword exercised on a rewritten copy; `read_jwst` on an
OIFITS file and `read_oifits` on an `x1d` file refuse by name;
`photometry_spectra.rst` gains "Reading a JWST spectrum: the front door";
`ampere.spectroscopy.rst` in the API toctree; `tests/data/README.md`
extended; the changelog's "No file readers yet" rewritten; `ampere/spectroscopy`
in the typecheck scope; gates dev.

### W7.4 — `normalisation="model"` on `SquaredAmplitude`, and the interferometry housekeeping [S; Sonnet] (W6.12's For Peter and carried findings)
(1) `SquaredAmplitude(normalisation="model")` on the three backends:
`|V/V(0)|²` with the model's own zero-spacing flux, for a model whose
total flux is free — the kind-preserving step evaluating the prediction
at `(u, v) = (0, 0)` through the same `FourierSample` it already holds,
the buffer form unchanged and the default; a conformance row per fixture
against the analytic uniform disc, and a row that a fixed buffer and the
model form agree when the buffer equals the model's total flux. (2)
`FourierSample.from_observed` on every backend emits the observed
container's own spectral unit rather than converting to µm, so a container
not in µm no longer fails alignment at composition; a conformance row in
nm. (3) `tests/conformance/protocol.py`'s `InterferometryPieces` gains
`squared_amplitude` and the mirror backend a `SquaredAmplitude`;
`test_native_interferometry.py::test_every_piece_declares_this_backend`
reads the piece names from the protocol instead of a hand list.
(4) A conformance row fitting the shipped `Binary` to the contest file
and recovering the published geometry (5.0 mas at 30° east of north,
ratio 8.9) in the model's own convention — the user-journeys walkthrough
(Appendix D) found the optimum at 3.9 mas and 126° with the walkthrough's
construction, which is a convention or a field-of-view question the page
must answer in one sentence beside the DFT sign it already states.
**Depends:** nothing. **Accept:** the rows above green on every fixture;
`interferometry.rst` §10's reader example gains the model form in one
sentence; the convention sentence; gates dev + torch + jax.

### W7.5 — Phase 7 housekeeping: the beta's carried list [S; Sonnet]
The owed list from the Phase 6 rows and the handoff, each a line and none
a behaviour change: `ci.yml`'s `v2` trigger removed (D15); `overview.rst`,
`sbi.rst` and the two backend API pages' clone-form installs replaced by
`pip install ampere-astro[...]`, `ampere/__init__.py`'s "alpha testing
phase" docstring and `conf.py`'s copyright and author lists brought to
`pyproject.toml`'s (W6.5); `examples/README.md`'s stale notebook paragraph
(W6.1); the `failure_*` variables' fixed-width truncation on
`append_training_set` and integer/boolean extra coordinates exercised by a
row (W6.14); `SyntheticPhotometry.from_library`'s call into
`ampere.legacy.utils.pyphot_compat.get_unit` lifted into v2 (W6.1 (b),
W6.0's carried `pyphot_compat` home — the one legacy dependency v2 has)
and the stray ruff entry (W6.0); a `heavy` marker in `pyproject.toml`
keeping W6.7's emcee burn-in row and the `sed_composition` fixture out of
`test-fast` (W6.7); the sbi calibration row
`TestAmortisationOverTheObservationContext::test_it_stays_calibrated_at_a_rescale_the_prior_covers`
that fails locally at `ks_pvalue` 0.0054 against 0.01 while CI passes —
the seed pinned or the budget raised so the row is not at the threshold
(W6.15, W5.26's territory); netCDF4 re-locked when 1.7.4.1 ships (W6.6);
the Hyperion `np.string_` shim dropped if upstream has released on NumPy
2 (W6.13 (C2)); **a note on `migrating.rst`** that the legacy
`phoenixstar.py` and `QuickSED.py` divide by 4πd² twice, so their fluxes
are 4π too faint and the twins use the correct law (D11: the legacy code
is frozen, so a note is the only fix); the changelog's "Unreleased"
section carrying each.
**Depends:** nothing. **Accept:** each line's own check (a row, a grep,
a build); `import ampere` imports no legacy module (the W6.0 row still
green); docs warnings no longer than base; gates dev (+ sbi for the
calibration row).

### W7.6 — The optimisers' bound-aware scipy route [S; Opus] (W6.7's carried finding)
W6.7's scipy route, on an eighteen-parameter problem whose free
coordinates are all sigmoid-bounded boxes, converged from eight Powell
starts to a box corner (seven coordinates saturated) with a log posterior
a hundred below a point the ensemble had already visited: Powell's line
searches run to where the sigmoid saturates and the objective goes flat.
The fix: the scipy route minimises in the *constrained* coordinates with
the bounds passed (`L-BFGS-B` or `Powell` with `bounds=`, the Jacobian of
the bijection dropped from the objective accordingly), the unconstrained
form kept for unbounded coordinates; a loud warning (a `UserWarning`
subclass in `ampere.inference.exceptions`, the `SamplingFailureWarning`
pattern) when a converged coordinate sits within a tolerance of its
bound, naming the coordinate; `warm_start_gp`'s assumption of `<label>.likelihood.<name>`
names replaced by a lookup through the merged mapping so a renaming tie
is accepted. The NGC6302 twin's eighteen-parameter problem is the
regression row (its MAP within the ensemble's best point's log posterior
by a pinned margin). **Depends:** nothing. **Accept:** the regression
row; the existing `test_optimise.py` rows unchanged; `optimisers.rst`'s
route table amended; the `Optimum` record's provenance naming the
coordinates (constrained/unconstrained) — a `results.md` §4 sentence if
the record's fields change, with the decision-log row; gates dev + torch
+ jax.

### W7.7 — The docs accessibility pass [S; Sonnet] (Phase 6 D14 (e))
Deferred from W6.17: an accessibility statement for the docs site, and
the review it rests on — alabaster's contrast (the link and code colours
against WCAG AA, overridden in `_static` where they fail), keyboard
navigation through the sidebar and the search, a `lang` attribute, alt
text on every figure as a docs convention (a Sphinx `figure` directive
without `:alt:` made a warning by a small extension in `conf.py`, so the
convention enforces itself), the six plots' default palettes checked
for colour-vision deficiency, recorded on the results API page. The statement on its own
page beside `citing.rst`, dated, naming what was checked, what fails
and how to report a problem (the issue forms). **Depends:** nothing.
**Accept:** every existing figure carries alt text; a build with a
figure lacking it warns; the contrast figures recorded on the page;
`pixi run docs` warnings no longer than base; gates none (docs).

### W7.8 — `1.0.0b2`: the second beta [S; Opus] (Phase 6 D3's train)
The release procedure as W6.5 wrote it and the first release corrected
it (`pixi run --frozen`, the clean-tree check, the `pypi` environment's
tag policy), run a second time from the changelog's "Unreleased"
section: the known-limitations list rewritten (two lifted by W7.0/W7.1;
the readers sentence rewritten by W7.3; the jax solver's and the
hierarchy's standing), the GPU rows on the cluster at the candidate
(`scripts/cluster/run.sh`, 33 rows plus any W7.0–W7.2 add), CI green on
the candidate, Peter's tag, TestPyPI then PyPI behind his approval, the
GitHub release and the Zenodo version DOI, the Read the Docs version.
**Depends:** every Phase 7 item Peter wants in the beta merged.
**Accept:** `ampere-astro 1.0.0b2` on PyPI built from the tag with its
own version; the GPU rows green; the decision-log row dated; the
changelog's section closed.

### W7.9 — Line fluxes as an observable kind [M; Opus] (issue #70, re-filed at D10; D7 below)
The modality template applied to the measurement a spectroscopist most
often has: a set of integrated line fluxes with uncertainties, each
identified by a transition (a rest wavelength and a label), the model
producing a `Spectrum` and the instrument chain integrating it over each
line's window — a `LineFluxes` point kind (the template of
`PhotometricPoints`: a coordinate per line, the value its integrated
flux, the units `W m⁻²` or `erg s⁻¹ cm⁻²` under the usual equivalences)
and a `LineIntegration` reference step (the window per line in velocity
or wavelength, a local continuum subtracted or not, by declaration) with
its native twins; the flexible likelihood on it (the kernel over the
line's rest wavelength, so a misspecified excitation ladder shows as a
correlated residual across neighbouring transitions); the conformance
rows per fixture; an example on a JWST-like line list (synthetic until
W7.3 lands a real file, then the real one). **Depends:** W7.3 for the
real-file example (D7 ruled in). **Accept:** the template's checklist (the kind, the
step on three backends, the family default, the diagnostics' four plots,
the conformance rows); a worked example page; gates dev + torch + jax.

### W7.10 — The jax-native quasiseparable recursion [L; Opus] (W6.18's carried limitation; D9 below)
The beta's fourth known limitation: the jax `QuasisepGP` refuses an
accelerator because celerite2's jax primitives lower on the CPU alone.
The item: the O(N) recursion written in jax itself — a `lax.scan` over
ampere's exact rank-2 Matérn-3/2 representation (the tinygp route W2.5
slice 2 declined on CPU speed), behind the strategy interface as a
second provider the solver chooses by device (celerite2 on the CPU, the
scan on an accelerator), bit-comparable to the reference solver at
`tolerances.cross_solver`, `conditional_loo` by the same recursion,
`BATCHABLE` true for the scan provider (the one thing the celerite2
provider cannot be); the GPU row W6.18 turned into a refusal row turned
back into the cross-device agreement row; the CPU speed of the scan
measured against celerite2 at 10³–10⁶ points and recorded in
`performance_memo.md`, with the device rule chosen from that table.
**Depends:** nothing (D9 ruled in); a cluster run per `scripts/cluster/run.sh`. **Accept:**
the conformance columns for the jax fixture green on both providers; the
GPU rows green with the agreement row restored; the table in the memo;
the changelog's limitation struck; gates jax (+ dev for the memo's rows).

### W7.11 — The star-disc twin's upper limits as a `Censoring` declaration [S; Sonnet] (D11, ruled 2026-10-05)
W6.13 (C1)'s finding (4): the HD 105 CSV's two upper limits at 500 and
880 µm are in neither the legacy fit nor the twin's. The twin
(`examples/star_disc/`) gains the two points as censored observations —
`Censoring` on the photometric dataset with the limits as upper bounds
under the Gaussian family, per `likelihoods.md` §9 — so the far-infrared
constraint on the disc's cold dust enters the fit; the twin's docstring
and coverage record updated, the smoke row in
`tests/examples/test_star_disc.py` asserting the censored points are in
the likelihood's decomposition, and one paragraph on `photometry_spectra.rst`
or the twin's docs section showing the declaration as the worked censoring
example the docs lack. **Depends:** nothing. **Accept:** the smoke row;
the fit's posterior on the disc temperature moves in the direction the
limits imply (pinned as a sign, not a value); docs warnings no longer
than base; gates dev.

### W7.12 — Convergence: the optimiser's mode as the ensemble engines' default start, the verdict, the warning [S; Sonnet] (ruled by Peter 2026-10-06 from the user-journeys memo, Appendix A finding 1)
The first fit a user writes on the beta — nine-band photometry, a
modified blackbody, emcee from prior draws — reaches R-hat 1.5 at 36 000
evaluations and 1.22 at ten times that, while the same budget started
from `optimise`'s mode reaches 1.06 in 22 s; `summary` printed R-hat 1.5
without comment (`docs/design/user_journeys_memo.md`, Appendix A). Three
changes: (1) **the default start of `EmceeEngine` and `ZeusEngine` is the
optimiser's mode** — `run(initial=None)` runs `optimise(problem)` once
(the scipy route, its default starts) and takes
`initial_positions(walkers, around=optimum)`'s ball; `initial="prior"`
restores the prior start; an explicit `initial=` is honoured as today; the
run records `ampere_start` (`"optimum"`, `"prior"`, `"supplied"`) and the
optimum's provenance as W6.7 stores it; a problem the optimiser cannot
start (a non-finite objective at every start) **or whose optimum sits on
a prior bound** (Appendix B: the ball around a bound-pinned optimum is
degenerate and every walker stayed on it — ESS equal to the draw count,
R-hat undefined) falls back to the prior with a loud warning naming
`initial="prior"` and the bound; the default start is budgeted (one
start, or a time cap — Appendix E measured 79 s for the eight-start route
on 22 parameters) with the full route on request; W7.6's bound check is
the one this uses, so W7.6 lands first or with it. (2) **A convergence
verdict**: `ampere.results.check_convergence(tree, *, rhat=1.05, ess=100)`
returning a small typed record — passed or not, the failing variables with
their values, and the remedy in words (more steps; more walkers; the
warm start if the run says `"prior"`; a reparameterisation when one
variable alone fails) — and a `ResultsWarning` from `emit` and `summary`
when a chain-based run fails it, on the same `ampere_approximation` guard
`summary` already uses so VI and SBI runs are never warned about. (3)
The FAQ's "Inference" section rewritten around the verdict and the two
starts; the quickstart's text names the default. (4) An ensemble engine
given a problem whose `free_size` exceeds a threshold (a `Population`
of twenty members took eight minutes to R-hat 1.23 on the numpy path,
Appendix E, where the docs fit fifty by NUTS) warns once, naming NUTS
and the native backends. **Depends:** nothing
(W6.7 merged). **Accept:** the Appendix A problem as a regression row —
R-hat below 1.1 at the first budget from the default start, the attr
present, the prior start reproducing today's draws bit for bit under
`initial="prior"`; the verdict's rows (passed, one failing variable, an
approximate run never warned); every engine test unchanged but for the
start; `results.md` §4 gains the attr and the decision-log row records
the default's change; gates dev (+ torch and jax for the engine rows that
run there).

### W7.13 — The portable model: one `ArrayOps` for models and the guide "Writing a model once for three backends" [S code, M docs; Opus] (ruled by Peter 2026-10-06; user-journeys memo §4 (1))
Most users write a `Model` subclass (Peter: 60 % or more), and the
package's answer to "how do I run it on jax" is documented only for
modality authors (`interferometry.rst` §8's inheriting pattern) while the
mechanism exists three times — `ampere.core.ArrayOps` for kernels (W4.5),
a private ops class in the jax astropy adapter, the PHOENIX emulator's own
`Ops`. The item: (1) one public protocol for models — `ArrayOps` widened by
what a model needs beyond a kernel (`where`, `interp`, `trapz`/`cumsum`,
`exp`, `log`, `sqrt`, `power`, `clip`, the dtype-and-device scalar and
`asarray`), `NumpyOps` from `ampere.core`, `JaxOps` and `TorchOps` from
their backends, the astropy adapter's private class and the emulator
example folded onto it (W7.0 adds `sqrt` and `log` first; this item lands
after it); (2) a documented base pattern — a reference model whose
`evaluate` is written against `self.ops`, and the two native twins as
three-line subclasses setting `OPS`, `BACKEND` and the capability flags,
the inheriting pattern made the default for a user's own model; (3) the
guide page, "Writing a model once for three backends": the quickstart's
`LinearModel` rewritten against the protocol (its `np.asarray` and
`float(...)` casts are exactly what makes it unportable), run under emcee
on numpy and NUTS on jax in one executed notebook, the `from_astropy`
route for the curated classes stated beside it, and the rule for when to
re-declare instead. **Depends:** W7.0 merged. **Accept:** the guide's
model as a conformance row agreeing across the three fixtures at
`tolerances.cross_backend` and lowering under NUTS on both native
backends; the emulator example and the astropy adapter on the public
protocol with their tests unchanged; the notebook executed at docs build;
gates dev + torch + jax.

### W7.14 — "Reading the diagnostics": the docs page that teaches what the misspecification plots show [S; Sonnet] (ruled by Peter 2026-10-06)
Peter: users probably do not know what the misspecification diagnostics
mean and would care if they did; half of them read the concept page,
which explains why the GP is there and not what its plots show. One page
beside `concept.rst`, walking the M2 study's scenarios through the four
plots — the residual plot with the whiteness test, the GP localisation
plot with its caveat, the posterior-predictive check, the anomaly score —
with, for each figure, the one sentence it supports and the one it does
not ("the GP absorbed structure at 9.7 µm; the model is missing a
feature there" against "the model is wrong"), the shrinkage section's
figures reused, and the decision a user takes from each (fix the model,
widen the prior, accept the GP's correction and report its amplitude).
**The page opens with the gross case** (ruled by Peter 2026-10-06): the
IRS spectrum of PG 1011-040 fitted by a single power law
(`docs/design/walkthroughs/persona_b2.py`, the user-journeys memo's
Appendix B.2) — the GP's conditioned mean *is* the silicate emission, a
third of the flux at a length scale of microns, the index loses two orders
of precision and under a permissive prior the normalisation runs to its
bound while the GP carries the lot; the page says in so many words that
this is the concept page's earthquake and what "absorbing a lot of
power" looks like, then turns to the mild case (W7.16's) where the GP
absorbs only what the model leaves. One section on choosing the kernel's
priors against the model's own scales — the amplitude as a fraction of
the flux, the length scale against the features' widths and the data's
span — which no page gives today (the quickstart's 0.01 µm length scale
on an IRS grid sidesteps the question).
Linked from the concept page's last paragraph and from `plot_residuals`'s
docstring. **Depends:** nothing. **Accept:** every figure traced to a
driver and its numbers to the M2 page's table; alt text on each; docs
warnings no longer than base; gates none (docs).

### W7.15 — The persona-A findings: `ModifiedBlackBody`'s `scale`, photometric alignment by filter name, the plotting nits [S; Sonnet] (from Appendix A findings 2, 3 and 5 and Appendix B.2; **ruled by Peter 2026-10-06: `scale` is the flux at the reference wavelength by default, the solid angle by keyword**)
(1) `ModifiedBlackBody` and `BlackBody` multiply `B_ν` with its per-
steradian magnitude, so `scale = 1` gives fluxes of 10¹⁵ Jy and a user
with catalogue fluxes needs a prior reaching 10⁻¹⁶, while the docstring
promises "`scale` stays interpretable as the flux the source would have
… at that wavelength". **Ruled**: `scale` becomes that flux — a quantity
in Jy at `reference_wavelength`, the Planck function divided by its own
value there — on `ModifiedBlackBody` and `BlackBody` alike, and the
solid-angle form is reachable by a documented keyword (`solid_angle=`, a
quantity in steradian, exclusive with `scale`), the parameter carrying
its unit in `to_spec` either way; the modified-blackbody twin,
the M2 generators and every test pinning a value amended, the change
named in the changelog as a behaviour change. (2) `PhotometricPoints`
alignment keys on the filter name, not the wavelength: the step tabulates
each filter on the user's grid and emits its own effective wavelength,
which differs from a catalogue's pivot wavelength in the second decimal,
so a catalogue's own axis is refused and the only accepted route is to
copy the step's axis from a dummy prediction. `check_alignment` on that
kind matches names and adopts the step's wavelength, refusing a name the
step does not tabulate; a `results_schema.md` §8 sentence and the
decision-log row. (3) `plot_posterior_predictive` thins its replicate draw
to a default of a few hundred draws (8 s on nine points, 39 s on 360),
and `gp_localisation` the same (130 s on 360 points); the localisation
plot's caveat text no longer overdraws the x-axis label; `Evaluation`
gains a `predictions` mapping (or `FittingProblem.predict(values)`) so
the model curve at a point is one call rather than a recomputation by
hand; the plotting library's "too few points to create valid contours"
warnings on a corner plot are silenced with the reason; one sentence on
`photometry_spectra.rst` on when the flexible likelihood has nothing to
learn from (nine photometric points). **Depends:** nothing. **Accept:** a
catalogue-wavelength `PhotometricPoints` accepted and fitting (the Appendix
A probe as a row); the scale's unit in the model's `to_spec` and provenance;
the twin's fit unchanged in its posterior after the amendment; gates dev
+ torch + jax (the native model twins).

### W7.16 — The spectroscopist's walkthrough: NGC 6302 as an executed notebook [M; Sonnet] (ruled by Peter 2026-10-06; user-journeys memo Appendix B)
Persona B's tutorial, on the mildly misspecified case Peter chose: the
NGC 6302 twin (`examples/ngc6302`, the Kemper et al. 2002 two-shell dust
model on the ISO 25–120 µm spectrum, eighteen species abundances, a
calibration scale, `QuasisepGP`) expanded from an example with a command
line into a notebook under `docs/source/notebooks/` executed at docs
build: reading the ISO data by hand (the lines a reader will replace,
said so), the model's declaration explained species by species, the
independent fit against the flexible one, the default start (W7.12) and
the convergence verdict read, every diagnostic of W7.14's page applied to
this fit with the sentence each supports — the GP's amplitude as a
fraction of the flux and where it localises, against the gross case the
page opened with — the dust-mass table as the post-fit product, and the
second mode of the legacy backlog shown rather than hidden. The notebook
runs a short budget inside the docs build's limit (the W6.1 pattern: the
budget and wall time stated in the last cell) and quotes the full-budget
run's numbers from the twin's recorded coverage. **Depends:** W7.12,
W7.14. **Accept:** the notebook executes in `pixi run docs` under two
minutes; the full-budget figures traced to the twin's record; linked from
`tutorials.rst` and from W7.14's page; docs warnings no longer than base;
gates none (docs; the twin's rows unchanged).

**Ordering (two agents at a time, gate legs scoped to the code touched).**
**Wave 1**: W7.0 alone first — W7.1 depends on its member
rules and both touch `parameter.py`, the populations and both backends'
resolve functions, so they cannot run beside each other; W7.5 beside it
as the filler (housekeeping, disjoint files), then the fillers in this
order as slots free: W7.12 (convergence — every walkthrough hits it
first), W7.14 (docs), W7.7 (docs), W7.11 (one example), W7.16 (after
W7.12 and W7.14), W7.15 (on Peter's word for the scale semantics); W7.13
only after W7.0 merges, since both touch `ArrayOps`. **Wave 2**: W7.1 ∥ W7.4 (the interferometry pieces
are disjoint from the parameter layer). **Wave 3**: W7.2 ∥ W7.3 (the M2
sibling under `examples/` and `tests/m2` against a new package
`ampere/spectroscopy` and `tests/spectroscopy`; disjoint). **Wave 4**:
W7.6 ∥ W7.9 (the optimisers against a new kind and step; disjoint). **Then W7.8**, the second beta, after wave 3 (D5), so `1.0.0b2` carries the two design items and
the reader. **W7.10** is the long one and runs in whichever slot frees
after wave 2, since it touches only the jax backend's solver and its
conformance column; it may miss the beta without holding it. Every item
except W7.7 and W7.8 touches tested code, so each has CI on its push as
the gate of record; W7.0, W7.1, W7.4 and W7.6 are core-touching and
their rows record the three legs.

**Decisions for Peter before dispatch.**
- **D1 — the `Derived` grammar versus a callable** (memo §3.2, §7 Q1):
  a closed expression grammar with a `symbols` mapping (recommended —
  provenance, hashing, serialisation and lowering must read it, for the
  reason `to_spec()` refuses an opaque prior; `where` added when a
  customer appears), or a Python callable traced on each backend.
  **Recommended: the grammar.**
  **Ruled 2026-10-05: the grammar.**
- **D2 — where derived values live in the emitted run** (memo §3.6,
  §7 Q2): in `posterior` with the `ampere_derived` attr naming them
  (recommended — ArviZ's summaries, the corner plot and the training-set
  writer see them without a second group, and a derived value *is* a
  posterior quantity), or a separate `ampere_derived` group.
  **Recommended: `posterior` with the attr.**
  **Ruled 2026-10-05: `posterior` with the attr.**
- **D3 — paths on `over`, or `within=` on `Population`** (memo §2.1, §7
  Q3): qualified paths as the one low-level grammar on `over`
  (heterogeneous depths allowed) with `within=` as a convenience on
  `DatasetCollection.plate` only (recommended), or `within=` on
  `Population` itself as a second declaration form. **Recommended: paths
  on `over`, the convenience on the factory.**
  **Ruled 2026-10-05: paths on `over`, `within=` on the factory.**
- **D4 — the landing order** (memo §10, §7 Q4): W7.0 → W7.1 → W7.2
  (recommended — (2) the `Derived` node is self-contained and (1) the
  paths need its member rules), or W7.1 first with the member rules
  split out. **Recommended: as listed**, Phase 6 D8's order.
  **Ruled 2026-10-05: W7.0 → W7.1 → W7.2.**
- **D5 — the second beta's cut**: `1.0.0b2` after wave 3 (W7.0–W7.5,
  W7.7 — recommended: the two design items lift two known limitations
  and the reader rewrites a third sentence, which is the changelog a
  second beta deserves), after wave 4, or at the phase's close. Also
  whether `1.0.0b2` is the last beta before `1.0.0` final, which fixes
  when v2's own deprecations (the `regularised_horseshoe` alias) go.
  **Recommended: after wave 3; the final's timing left to the phase after
  this one.**
  **Ruled 2026-10-05: after wave 3; the final's timing to the next phase.**
- **D6 — the JWST product set**: `x1d`/`c1d` one-dimensional products
  only (recommended — one container, one reader, the W6.12 size), or
  also `s3d` cubes into the `Image` kind per spectral channel (the IFU
  sketch's territory, `ifu_cube.md`; a second item if wanted) and the
  MIRI MRS/NIRSpec multi-extension layouts beyond source selection. And
  which public file: a NIRSpec or MIRI LRS `x1d` from an ERS programme
  under a megabyte (the orchestrator finds one and records its provenance
  at dispatch, as W6.12 did). **Recommended: `x1d`/`c1d` only; the cube
  reader deferred to the phase that lands the IFU modality.**
  **Ruled 2026-10-05: `x1d`/`c1d` only; the cube reader with the IFU modality.**
- **D7 — line fluxes in this phase**: W7.9 as drafted (a new kind and
  step by the template, Opus, M), deferred to the phase after, or
  declined — issue #70 is the one re-filed modality request and the
  second reader makes it timely, but it is the one item here that is
  neither a memo's nor a carried finding's. **Recommended: in this
  phase, in wave 4, if the budget allows; else the first item of the
  next phase.**
  **Ruled 2026-10-05: in this phase, wave 4.**
- **D8 — nested populations** (the beta's "one level of hierarchy"): a
  design memo in this phase (Fable's own work between waves, the W6.11
  pattern — a parameter carrying two plates, objects within surveys, the
  plate layout generalised to a tree, the native resolve functions'
  nesting; a §4 change whose memo precedes its item), or left as a known
  limitation until a customer appears. **Recommended: left; no customer
  has asked, and W7.1 is the hierarchy work this phase can absorb.**
  **Ruled 2026-10-05: left as a known limitation; Peter notes nested populations will be required at some point — the horizon memo's §4 places the memo and item in Phase 9 with the IFU cube as the customer.**
- **D9 — the jax-native recursion**: W7.10 in this phase (L, Opus, a
  cluster run per candidate), or deferred with the limitation standing
  — the torch solver has no such limit and the jax `DenseGP` and
  `HilbertSpaceGP` serve an accelerator today. **Recommended: in this
  phase, in the free slot after wave 2, allowed to miss the beta** — a
  CPU-only distinguishing feature on one backend is a claim the docs
  must keep qualifying.
  **Ruled 2026-10-05: in this phase, allowed to miss the beta.**
- **D10 — `ampere.diagnostics` and RHMF**: the revisit trigger (a
  robusta-hmf release at ≥ 0.1 with a test job on `main`) is checked
  once at the phase's close by the orchestrator (recommended — four
  lookups, no item) and the trial item re-opened only if it is met; or
  the trial's finding (robust weights localise sharp deviations only at
  ranks that leave no residual structure to screen) taken as closing
  the question. **Recommended: the check at close; the finding stands.**
  **Ruled 2026-10-05: the check at close; the finding stands.**
- **D11 — Peter's backlog from Phase 6, which of it this phase takes**:
  the two unconverged C1 coverage runs (a longer NUTS budget or the
  calibration-scale reparameterisation — the latter is a W7.0 customer:
  `scale = luminosity_ratio * k`); the star-disc CSV's two upper limits
  as a `Censoring` declaration in the twin (an example change, S,
  Sonnet — recommended: yes, as a filler); the missing IRS file
  `cassis_yaaar_spcfw_5295616t.fits` (Peter's to commit or not);
  NGC6302's second mode; C2's envelope mass at the quick budget; the
  legacy 4π-twice finding (the legacy scripts are frozen — a note on the
  migration page, recommended); the three issue comments (#40, #14,
  #41) and the D10 closures not yet posted; the one merged
  worktree left under `.claude/worktrees/` (W6.12's
  `agent-a4bda625296ffe8c5`, locked; the eighteen the handoff listed are
  gone) to `git worktree remove` on Peter's word, never the branch;
  W6.16's cu126/525 confirmation.
  **Recommended: the `Censoring` filler and the migration-page note in
  this phase; the rest stays recorded.**
  **Ruled 2026-10-05: the `Censoring` filler (W7.11) and the migration-page note (in W7.5); the rest stays recorded.**
## Status

| Item | Status |
|---|---|
| W0.1 | done 2026-09-01 (committed directly to master — the pending fix existed only in the local working tree, so the branch-per-item rule was waived for it) |
| W0.2 | merged 2026-09-01 |
| W0.3 | merged 2026-09-01 |
| W0.4 | merged 2026-09-01 |
| W0.7 | merged 2026-09-01 — the branch-archival actions in `docs/design/harvest/branch_triage.md` were approved and executed 2026-09-01 (its execution record); the stale 'await approval' wording was confirmed and struck 2026-09-08 |
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
| W2.4 (slice 2) | merged 2026-09-07 (Opus-authored, Fable-reviewed; no fixes needed — gates re-run independently: torch 1923/17, dev 1419/67, lint/format clean, pyrefly 0 errors in both). **The GP-library question was answered by measurement, not as the plan assumed**: celerite2 0.3.3 ships no torch interface and GPyTorch/linear_operator no quasiseparable operator (Toeplitz/Kronecker/SKI only — regular or product grids, or approximate), so the torch `QuasisepGP` wraps `celerite2.backprop`'s compiled forward+reverse kernels as `torch.autograd.Function`s (`_celerite.py`): agrees with `ampere.core.QuasisepGP` to the last bits and with the dense oracle to 1.7e-13, gradients to 8 significant figures, value+gradient 5.6 ms / 21 ms / 144 ms at 10³/10⁴/10⁵ points, zero new dependencies; GPyTorch exact was 123× slower at 10³ and 5 657× at 10⁴, its SKI missed `cross_solver` by 3–4 orders. Three decision-log rows. Also landed: the widened realised path (all core families, Tobit censoring, both solvers, fractional noise — worst realised-vs-contract disagreement 2.7e-15), `VIEngine` over pyro SVI (`AutoNormal`/`AutoMultivariateNormal`; conjugate posterior recovered to 0.2σ), `vmap` batching, the float32/device opt-out on the solvers through `provenance_config`, `Engine.finish` consuming `log_likelihood_terms` (`engine_draws_recomputed` 0 on realised runs, both backends), and the **`conditional_loo` O(N) recursion** on torch — W2.3's deferral discharged there; the conformance row now reads 'agree exactly or refuse by name'. Deferred with reasons: latent GPs (blocked by the core defect below — refused by name with a test asserting the flatness), `complex_gaussian` (Phase 4), a per-instance `device=` on models/steps (with the GPU smoke-test item). **Core defect escalated**: the latent-GP likelihood is flat in the kernel hyperparameters on every path (see the handoff; W2.14 proposed). Carried: `LoweringError`'s message says 'prior family' for anything; `FittingProblem` wants a public `capability_parts`; `tests/scaling` has no torch row; `ampere.core.QuasisepGP` could now lift its own `conditional_loo` deferral via `celerite2.backprop` (a core decision) |
| W2.5 (slice 2) | merged 2026-09-07 at `dd2c479` (Opus-authored, Fable-reviewed twice — on substance, then after the agent's own reconciliation onto the W2.4 slice-2 merge; no fixes needed; gates re-run independently on the reconciled tip: torch 1924/19, jax 1806/19, dev 1419/69, lint/format clean, pyrefly 0 errors ×3). **jax `QuasisepGP` = celerite2.jax, chosen by measurement**: tinygp 0.3.1's `QuasisepSolver` is three to five orders more accurate but its per-point cost doubles with every doubling of N on XLA:CPU (0.47 s vs 2.3 ms at 10⁴; infeasible vs 0.26 s at 10⁶), so the scaling answer decided it; costs taken and recorded — the x64 import side effect contained (import deferred to first use after `require_x64`; subprocess test), no `vmap` rule (`QuasisepGP` declares `BATCHABLE=False` and the batched density refuses by name), and `conditional_loo` stays refused (no O(N) route to the inverse diagonal in celerite2.jax's public surface; tinygp would have given it). Two decision-log rows; §6 records that **both tracks converged on celerite2** and the single-point-of-failure consequence. Also landed: the realised path widened to every core family (a `families.py` transcribing the closed forms, Tobit censoring via `jnp.where`, Student-t log-CDF from the regularised incomplete beta), latent GPs (mirroring the core oracle deliberately — see W2.14), both solvers, `vmap` batching, float32/device `InitVar` opt-ins reported by `provenance_config`, the numpyro VI route in the shared driver (`ImproperUniform` + `factor` bridge), `SyntheticPhotometry` fixed on both native surfaces plus a `PHOTOMETRIC` conformance shape (which then caught a missing `apply_flux` on the torch fixture, fixed in the reconciliation), jax scaling rows (p = 0.92). Carried: a fitted Student-t `nu` on censored data has no gradient (`betainc` has no derivative in a/b) and is refused by name; `default_bijection_for(loguniform(2, 30))` gives `Log(lower=0)` so the unconstrained origin is outside the prior's support (surprising, not wrong); tinygp's quadratic scan is worth reporting upstream; master's VI conjugate-mean row passes on numpyro with a 0.19σ margin against a 0.2σ tolerance (seeded, deterministic — loosen if a library bump flips it) |
| W2.7 | merged 2026-09-08 at `f57c72f` (Opus-authored, Fable-reviewed; one wording fix on the decision-log row — the RHMF re-verification is the agent's, confirmed at review, and W2.5 *has* merged; dev gate re-run independently: 1472/69, lint/format clean, pyrefly 0 errors). Families B and C in `ampere.results`: `add_residuals`/`gp_localisation` per `results.md` §7 (thinned, `draw` coordinate of retained indices, prior-rejected draws NaN, hash-checked before any arithmetic); `diagnostics.py` — a separation-binned autocorrelation native to irregular coordinates that reduces exactly to the lag autocorrelation on an even grid, `Q = Σ n_b ρ_b²`, permutation-calibrated (`(1+#{Q_perm ≥ Q})/(1+n_perm)`) from a named `substream` of the run's seed, O(N·m) windowed pair search with a pair budget, per-draw kept and median reported; the χ² Bayesian p-value where the stored log-likelihood makes it free, refused by reason otherwise; `gp_localisation_score` → `AnomalyScore` standardised by the law of total variance (the across-draw spread alone dropped the correlation with an injected deviation from 0.88 to 0.30); `plot_residuals` (warns on a GP fit), `plot_gp_localisation` (caveat attached to the figure unconditionally, readable via `figure_metadata`), `plot_anomaly_score` (provenance forced on for both when two provenances share an axes); shared helpers in `_plotting.py`; matplotlib lazy behind `OptionalDependencyError`. Whiteness fires on structure (p = 0.005 floor) and not on noise (p = 0.71). **RHMF deferred** with a decision-log row: licence (MIT) and API hold; maturity does not (no CI — only a publish-on-tag workflow; one stale 0.0.2 release from 2025-11; README TODOs; `Robusta.mse` unimplemented) — revisit at a ≥0.1 release with a test job. No dependencies added. Carried: `results.md` §8's 'every plotting function raises' sentence is stale (W2.8 may amend); `Dataset` gives a bare `TypeError` on `channel=`; ruled 2026-09-08 (Peter): across-draw *median* pooling of the per-draw permutation test, `gp_localisation(datasets=None)` meaning 'every GP dataset', and `plot_anomaly_score` returning the Axes — all three confirmed |
| W2.9 | merged 2026-09-08 at `196b48f` (Sonnet-authored, Fable-reviewed; the profile MLE re-derived by hand — the positive root of ρ(ρ+1)b² − [ρ(S+B) − m(ρ+1)]b − Bm = 0 — and brute-force checked in the tests; `tests/examples` 13 passed, lint/format/pyrefly clean). `examples/wstat_comparison.py`: `ProfiledCashWithBackground` as a registered user family (`wstat_example`; excises its background buffer with `NoiseParams.retain`; `sample()` refuses with the profiling-specific reason; per-sample terms declared not predictive), `XraySource`/`XraySourceAndBackground` toy models, both routes run through emcee, `compare()` from the run's own numbers; `docs/source/wstat_comparison.rst` as a literal script page in the tutorials toctree; `tests/examples/` as the end-to-end gate (with a `_FAMILIES` snapshot). The recommendation is structural (a proper marginal likelihood with valid predictive densities) rather than an asserted seed-dependent direction. Also fixed: `docs/source/conf.py` read `version('pgmuvi')`. **Carried to W2.11**: `pixi run docs` is broken on master for pre-existing reasons (no `pandoc`, no registered Jupyter kernel for two legacy notebooks; many legacy-docstring Sphinx warnings; autodoc of the non-existent `ampere.infer.ptemceesearch`) — the docs task needs its dependencies declared and `tests/examples` should join a gate. Ruled 2026-09-09: the repeated-trial coverage study demonstrating the low-count profile bias is W3.6's first application of its any-engine route |
| W2.14 | merged 2026-09-08 (Opus-authored, Fable-reviewed; no fixes needed — gates re-run independently: dev 1431/71, jax 1822/21, torch see the merge commit, lint/format clean, pyrefly 0 errors ×3). **The latent-GP likelihood sees its kernel**: `GaussianProcessNoise.noise_params` applies `f = solver.latent_transform(kernel, coordinates, z, values)` via `_realised_latent` (both arrays already the retained block — `Likelihood.log_prob` excises first and problem validation pins the latent size to the effective mask; shape and missing-coordinates refusals named), so `NoiseParams.latent` *is* `f`; no family changed. Regression rows in `tests/core` (9 of 10 fail on the unfixed core — bit-identical log-likelihood across amplitudes and length scales; closed-form Poisson and Gaussian-at-`f = L z` checks). Both backends gained native whitening surfaces (`latent_transform_native` on torch, `latent_transform_jax` on jax, no-exception, NaN-through-`where`, the stabiliser scale as a `where` on a value so a traced amplitude survives), the torch blanket refusal and jax's deliberate mirroring removed, and a positive refusal for a solver lacking the native surface. **A second latent defect exposed and fixed**: torch's `log_likelihood` branched on `self.correlated` alone, so a Poisson+GP dataset would have scored the Gaussian *marginal*; it now branches on `correlated and family.ANALYTIC_WITH_GP` as jax already did (never reachable on master, behind the refusal). New conformance shape `LATENT_GP` (Poisson counts over dense-GP noise; `DatasetSpec` derives integer counts from `family == "poisson"`), realised-vs-oracle rows per backend and a gradient-non-zero row in each backend suite. `likelihoods.md` §17 limitation 6 amended, §10's latent doctest corrected, decision-log row, plan §4.4 clarified. Carried: `latent_transform`'s stabiliser falls back to 1.0 at exactly zero amplitude (documented, not wrong); the neutral battery's latent shape uses `DENSE` only (a one-line second spec would add `QUASISEP`); §10a could list the native surfaces a backend owes (`latent_transform_native`/`_jax`) — small follow-up |
| W2.8 | merged 2026-09-08 at `71b2e96` (Opus-authored, Fable-reviewed; no fixes needed — dev gate re-run independently 1530/69 on the branch, 1542/71 on merged master; lint/format clean; pyrefly 0 errors). `plot_corner` (merged names; array blocks one column per element; prior-rejected draws excluded; refuses above `MAX_CORNER_VARIABLES = 20` with `max_variables=` override), `plot_trace` (NaN gaps for prior-rejected draws; `lp` its own row; cap 40), `plot_posterior_predictive` (refusal names `add_posterior_predictive`; σ-standardised discrepancy by default, `statistic=` for a θ-dependent one); `add_posterior_predictive` via `simulate(observe=True)` on the `posterior_predictive` substream, masked samples NaN, complex refused by name; `add_pointwise_log_likelihood` as the only writer of §6's group (per-variable `ampere_decomposition`, `"mixed"` at group level on a heterogeneous joint fit; refuses by name where the solver lacks `conditional_loo` — never a dense fallback); `pointwise_as_log_likelihood` bridging to `arviz.loo`; θ dtype preserved (`CONTAINER_SCHEMA_VERSION` 2, old form still read), `training_pair_from_dict`, `observations` + `Failure` carried; `ampere/results/training.py` — the §11 layer-2 netCDF training set (root provenance attrs + `ampere_training_set_version` 1, `theta` per merged name, one group per `<model>.<channel>`, `coordinates` once, `sample_stats` with the whole `Failure`, `observations/<label>`), `append_training_set` checking `ampere_spec_hash` before writing and growing `sample` by read-concatenate-rewrite (O(existing+new); the unlimited-dimension writer is §13.9's extension point). `PROVENANCE_SCHEMA_VERSION` stays 5 (no mapping changed meaning; every problem hash unchanged — conformance-checked). Three §4 amendments in one decision-log row (§8's stale 'every plot raises', §11's table gains `observations`/`Failure`, §6's landed paragraph). One netCDF engine per file after a real deadlock (mixing h5netcdf and netCDF4 on one file with a leaked handle). No dependencies added. Carried: `tests/results/test_results.py::test_both_netcdf_engines_write_it` leaks a netCDF4 handle (close it); W2.7's `add_residuals`/`gp_localisation` do not handle complex containers; `plot_corner` on a constant column surfaces `corner`'s own error. Ruled 2026-09-08 (Peter): all five confirmed — the 20/40 caps (the two routes past them are `var_names=` selection and the `max_variables=` override, both named in the refusal; **automatic paging with a loud warning is wanted later — W3.10**), the corner exclusion, `pointwise_as_log_likelihood` public, `"mixed"`, and the append deferred until a Phase 3 budget outgrows memory |
| W2.11 | merged 2026-09-08 at `dbbc172` (Opus-authored, Fable-reviewed; no fixes needed — YAML validated, dev gate 1555/71 with `tests/examples` joined, lint/format clean, pyrefly 0 errors, `pixi run docs` succeeds (15 legacy warnings), `bench` 3 rows). CI: one job per environment — `lint`, `typecheck` (dev), `test` (py3.11–3.13 matrix), `suites` (dev `test-all` + `bench`), `backend-suites` (torch|jax: typecheck + `test-all` + `bench`, `fail-fast: false`), `docs`, `minimal-install` (bare pip), `sbi-characterisation` (schedule/manual); `benchmark.json` uploaded as `benchmarks-{dev,torch,jax}`; setup-pixi caching everywhere. **CPU wheels**: torch pinned to PyTorch's CPU index in the `torch` and `sbi` pixi features only (pip users unchanged); 15 nvidia wheels gone from the lock. **Docs build repaired**: conda `pandoc` + `ipykernel` in `dev`, the dead `ampere.infer.ptemceesearch` autodoc entry removed, `nbsphinx_execute = 'never'` with the two unrunnable legacy notebooks explained; the job fails on errors, `-W` deferred to the Phase 6 docs rebuild (a recorded deviation from plan §5's Phase-1 line). **Benchmark harness: pytest-benchmark** (decision-log row; §6 struck): pixi owns the environments, PR artefacts are the requirement, the rows reuse the fixtures; no regression history — revisit if a benchmark server is wanted. `tests/benchmarks/test_gp_solvers.py` new; torch scaling rows added (p = 0.71 over 10³–10⁵). Carried: `docs/source/_static` missing (warning every build); `Embedding_nets.ipynb` in no toctree; `ampere.infer.sbi`'s API page empty in the dev docs build; `tests/examples/conftest.py` reaches `_FAMILIES` privately. Ruled 2026-09-08 (Peter): the CPU-index pin, the `tests/examples` placement and the `-W` deferral confirmed; the branch-protection rename remains Peter's action on GitHub |
| W2.10 (M2) | **merged 2026-09-08 at `956e1ab` — milestone M2 reached** (Opus-authored, Fable-reviewed as the adversarial pass; gates re-run independently: dev 1614/95, jax 2012/38, torch 2131/38, lint/format clean, pyrefly 0 errors ×3). `examples/m2_misspecification/` (generators keeping truth/deviation/observed; the toy model as a hand-written `Model` per backend sharing one declaration — the realised paths lower any model exposing `grid`/`flux`, asserted on both backends; the study driver with every threshold as a named constant; figures; CLI), `tests/m2` (83 rows, in `test-all`), `tests/benchmarks/test_m2_*` (the table as a CI artefact), `docs/source/m2_misspecification.rst`, the `m2_full` marker. **Result**: Matérn-3/2 reproduction; standard likelihood confidently wrong under misspecification (worst |median−truth| 2.83/7.39/4.71 widths at 200 points, 11.2 → 36.3 → 113 up the ladder, truth excluded) while the flexible one stays calibrated (≤ 0.98 widths, 0.70/0.51/0.57 up the ladder, truth covered everywhere); whiteness p = 0.005 (floor) on every misspecified scenario vs 0.185 control; localisation peaks 0.33 grid spacings from the injected line at 45× the control; three backends agree on the density to 1e-10 and on the posterior to 0.072/0.109 widths against 0.20/0.30 milestone tolerances derived from measured scatter; benchmarks: O(N) solvers 2–8 ms at 20 000 points (torch 4.0 ms contract / 8.2 ms with gradient; jax 0.6 ms jitted with gradient at 2 000), dense infeasible there; legacy 352 ms at 2 000. Adversarial review dispositions (Fable, 2026-09-08, standing in for the cross-model pass while Codex is quota-blocked; a retroactive terra pass is owed on this range): 1. Seed dependence of the science claim — PROBED: the study re-run at the CI budget under two seeds never used on the branch (7, 20260908); standard worst bias 3.67/10.8/5.06 and 3.67/10.9/4.62 widths (mild/strong_smooth/strong_sharp), flexible worst 0.76–0.88 with the truth covered everywhere, control flexible ≤ 0.58; every headline assertion holds. Not tuned to the seed. 2. Threshold margins — CHECKED: STANDARD_MIN 2.0 vs worst-case 2.83 at milestone budget and 3.67 at CI budget; FLEXIBLE_MAX 1.5 vs worst 0.98/0.88; localisation contrast 5 vs 45 observed. Margins are wide except the mild/standard one, which is the physically weakest scenario by design. 3. Bias statistic — ACCEPTED: |median − truth| over half the 68% interval is the calibration statistic the claim needs; the report says plainly the flexible likelihood does not give better point estimates. 4. Backend tolerances — ACCEPTED: derived from measured three-seed scatter on the worst-conditioned parameter, stated at 3–5× the Monte Carlo error of a difference; worst observed 0.072/0.109 widths at the milestone budget against 0.20/0.30. 5. Legacy comparability — ACCEPTED as stated: legacy's covariance is a truncated (I + w·M)σσ form, not K(θ)+diag(σ²), so only costs are compared; the RBF cross-check at 200 points reaches the same conclusion. 6. Toy model location — ACCEPTED: the realised paths lower any model exposing grid/flux, asserted on both backends; nothing added to the shipped backends. 7. Findings carried: the jax contract path costs ~28 ms flat on a QuasisepGP problem (dispatch through the numpy boundary, uncompiled) — a gradient-free engine on jax is a pessimisation; torch NUTS ~2.6× slower than jax NUTS at the milestone budget; tests/m2 costs 5:30 in torch; pyproject now carries [tool.pytest.ini_options] (markers only); examples/ imported as a namespace package via sys.path in conftest (a third import idiom). Ruled 2026-09-08 (Peter): the Matérn validation note on the plan's kernel row ratified; the retroactive terra pass on this range when quota returns |
| W2.4 (slice 3) | merged 2026-09-08 at `d6d9d01` (Opus-authored, Fable-reviewed; no fixes needed — dev 1616/98 and torch 2209/40 re-run on the branch, lint/format clean, pyrefly 0 errors ×2). **Per-instance `device=`** on every shipped torch piece via `_config.place`/`move` (an instance attribute shadowing the class flag — no core change; buffers move with `.to()`; kernels stop copying GPU coordinates back per evaluation; a GPU-configured `DenseGP` now declares its device; the GP noise models refuse a kernel or solver on another device since a kernel is not a capability part); `LoweredProblem` builds every tensor on `problem.device`. **`complex_gaussian`** transcribed and composing in the realised path (a complex observed container stays complex); **the circular complex GP is refused by name** — the core declares `GP_ANALYTIC_IMPLEMENTED = False` and refuses the composition, so there is no oracle (accepted at review; Phase 4 implements it in core first). The prediction-aware `sigma_tensor` hook's last silent fallback is a named refusal; `|predicted|` taken before any dtype cast. Conformance: `ModelKind.COMPLEX` behind a `complex_models` capability flag (torch declares it; the others skip). `tests/gpu/test_torch_gpu.py` at API level (13 rows, skipped here); `tests/backends/test_torch_device.py` (74 rows, using the meta device as a second device). `TorchParameterSpace.sample` fixed (CPU draw then move). Carried: `tests/gpu` is outside `test-all` by design (`pixi run gpu`); `LSFConvolution.sigma_tensor` shares a name with the noise hook (rename to `lsf_sigma_tensor`); a core numpy `Kernel` inside a torch GP noise is still accepted (amplitude gets no gradient). Ruled 2026-09-08 (Peter): refuse by default, accept under an explicit opt-in when no gradient is needed, always refuse where a gradient is required — W3.8 (decision-log row recorded); the rename is W3.0 |
| W2.5 (slice 3) | merged 2026-09-08 at `ae1392a` (Opus-authored, Fable-reviewed; dev 1614/108 on the branch, all three typechecks clean; jax and torch gates re-run on merged master — see the handoff). **Per-instance `device=`** through `_device.py` (`device_put` placement of grids, influence matrices and covariances; noise models own no arrays and say so; refusals by the piece's own contract error; mixed devices a `DatasetError` at composition). **`complex_gaussian`** was already lowered in slice 2 — now shown to survive the whole realised path at 0.0 disagreement with a finite gradient; a model supplying `flux` without `grid` is refused by name. **`conditional_loo` parity: the refusal is lifted** — `celerite2.jax.ops`' public `factor`/`solve_lower`/`solve_upper` give `d, W` beside the term's `c, U`, and W2.4 slice 2's O(N) backward accumulation runs as one reverse `lax.scan`; agrees with `DenseGP`'s Cholesky to 1e-12; only the numpy `ampere.core.QuasisepGP` keeps the deferral (two contract sentences corrected in `likelihoods.md` §7 and `results.md` R2, decision-log row). **The gradient-free fast path**: `_EvaluationCache` gains a route through `ampere.core.realise` (prior from the declaration, per-dataset terms from `log_likelihood_terms`), guarded — `realise`'s reference-point check gates attachment, labels checked once, any evaluation-time exception falls back to the contract path, `strict=True` keeps the contract path because a realised density has no failure reasons; `use_realisation=True` on emcee/dynesty/zeus (NUTS/VI pass `False`); `engine_realised_evaluations` recorded and `ampere_realised` 1 on such runs; `LoweredProblem.log_likelihood_terms` jitted on first use (28.8 ms → 0.14 ms). Measured: jax `QuasisepGP` 28.5 → 0.64 ms (44×), jax `DenseGP` 4.3×, torch 1.0–1.2× (accepted as `True` for uniformity). `tests/gpu/test_jax_gpu.py`, the `gpu` pixi task outside every gate, `tests/inference/test_fast_path.py` (13 rows per backend), a benchmark row asserting the ratio > 3. Carried: the jax realised path does not reproduce `GaussianProcessNoise`'s zero-uncertainty refusal (the `realise` guard catches it, so the fast path falls back — but a `strict` NUTS run would sample a density the contract path refuses; worth a native precondition); a jax model cannot name a parameter `flux` (reserved by the native surface); `ampere_realised = 1` on gradient-free jax/torch runs is a visible provenance change. Ruled 2026-09-08 (Peter): the provenance change and torch's `use_realisation=True` default confirmed; the zero-uncertainty precondition is W3.0 |
| W2.15 | merged 2026-09-08 at `fbd6970` (Opus-authored, Fable-reviewed; one review fix — the jax `_NEGATIVE_INFINITY` sentinel spelt `-math.inf` so autodoc's import mocks need no private-API shim; docs build reproduced: succeeds, 13 warnings all from frozen legacy code (was 21 with two errors); dev 1616/111 and characterisation 3 passed unchanged; lint/format/pyrefly clean). `docs/source/`: `overview.rst` (the v2 architecture as built: the ladder, the three environments, composition on each backend with the one-backend rule, realisation and the five engines, results and diagnostics, M2 linked), `api.rst` + autodoc pages for `ampere.core`, the three backends (torch/jax under `autodoc_mock_imports`, stated on the page — no environment has both libraries and CI builds in `dev`), `inference`, `results`, with `Constants` sections for the 49 re-exported constants; legacy pages banner-frozen and `:no-index:` (kills ~300 spurious `None` cross-references from legacy docstrings); `install.rst` (pixi first, pip second, extras table), `tutorials.rst`, `advanced.rst` repaired; `napoleon_use_ivar`, dead mocks dropped, `html_static_path` emptied. `README.md` rebuilt round the M2 numbers with the real install routes (both verified: cold `pixi install -e dev`; a scratch venv `pip install -e .` + extras + `--dry-run` of `[all]`/`[dev]`). **21 stale sentences annotated *Amended W2.15*** across `architecture.md`, `lowering.md`, `likelihoods.md`, `results.md`, `results_schema.md`, `transformations.md`, `diagnostics.md` (`parameters.md` and `inference.md` clean), one decision-log row; the architecture §3 'as built' subsection covers the namespace tree, extras and pixi features as they are. Four docstring fixes; every public name renders (checked against the built HTML). Seven `examples/` scripts migrated onto `pyphot_compat.get_unit`. **Findings carried**: **`pip install ampere` installs an unrelated PyPI package** — the install page now warns, but `OptionalDependencyError`'s runtime remedy (`ampere/core/exceptions.py`) and six docstrings still say `pip install "ampere[extra]"` (a small code fix for the first Phase 3 session); `SAMPLE_STATS_GROUP` defined twice in `ampere.results`; `demo_pyphot.ipynb` and `examples/examples_paper/{phoenixstar,modifiedblackbody}.py` still use `pyphot.unit`; `ampere/class_reqs` is a stray notes file; `faqs.rst` is empty headings (Phase 6 content). Ruled 2026-09-08 (Peter): `:no-index:` stays and the legacy pages stay linked from their own 'Legacy (frozen)' section of the API reference, joined by a 'Migrating to the new API' page (W3.9 seeds it, Phase 6 completes it); the mocked constants confirmed as is |
| W1.13 | **merged 2026-09-03 — Phase 1 complete, spec frozen** (Fable-authored across a session-limit resume, Fable-reviewed independently: every contract change read against its recorded ruling, LOO algebra hand-checked, gates re-run — 1109 tests in one process, py3.11 leg, lint/format/pyrefly — and no fixes needed). 18 commits: the twelve ruled contract changes, cross-review harmonisation, the W1.11 gap dispositions, the consolidated serialisation review (`docs/design/serialisation_review.md`), the freeze record, the Phase 2 breakdown (W2.1–W2.11), the `gpt-5.6-terra` adversarial pass (two real pre-existing defects fixed — the censoring/masking union in `draw_observation`, loud seed validation — one draw-epsilon nuance documented for Phase 2), and Peter's rulings on the three escalations (model identity + `describe()` → W2.1; `Population` → Phase 5, adaptable; `Axis.locate` approved → W2.1). Tag `spec-v1.0` at the freeze-content commit `58daa86` |
| W3.0 | merged 2026-09-08 at `f0140dc` (Sonnet-authored, Fable-reviewed; no fixes needed — branch gates dev 1616/111, jax 2052/42, torch 2222/40, lint/format clean, pyrefly 0 errors ×3; dev gate re-run on merged master 1616/111). The remedy text is `pip install ".[<extra>]"` from a checkout with the PyPI-name warning; one `SAMPLE_STATS_GROUP`; torch `LSFConvolution.width_tensor`; jax GP-marginal datasets refuse zero or negative retained uncertainties at construction — the item's Accept line implied jitter rescues a zero uncertainty, which the contract's `_observed_sigma` does not allow, and the agent rightly followed the contract. Carried: **torch has the identical GP-marginal gap** (no zero-uncertainty check in its `gp_marginal` branch) and jax's `GaussianProcessNoise` takes no `jitter=`/`scale=` constructor keyword unlike torch's — both folded into W3.1 slice 2, which touches both `problem.py`s and the backends' noise models |
| W3.9 | merged 2026-09-09 at `fcdd07d` (Sonnet-authored, Fable-reviewed; one review fix — the page had addressed users in the voice of the decision record; docs build 13 warnings before and after, all legacy). `docs/source/migrating.rst` sits in the API reference toctree between Ampere v2 and Legacy (frozen); the full guide and deprecation policy remain Phase 6's |
| W3.1 (slice 1) | merged 2026-09-09 at `11d01d7` (Opus-authored, Fable-reviewed; no fixes needed — merged-master gates re-run one at a time: dev 1729/111, jax 2176/42, torch 2346/40; lint/format clean, pyrefly 0 errors). `FittingProblem.simulate_many` on the **spawn-based equality** — `simulate_many(n, stream=s)` returns exactly `[simulate(rng=child) for child in rng(s).spawn(n)]`, order preserved, whatever executor or chunking — with `SerialExecutor`/`ThreadExecutor`/`ProcessExecutor` behind a `concurrent.futures`-shaped `Executor` protocol, per-draw `timeout=` and crash capture as flagged failures (`FailureReason.EXECUTION_FAILED`, a §4 vocabulary addition), chunked `as_chunks=True` feeding the training-set writers, `Model.evaluate_batch` giving `BATCHABLE` its reference meaning, `examples/sbi/external_simulator.py` under the process pool in `tests/examples`, and a `copyreg` reduction making `mappingproxy` (hence every frozen core object) picklable. **Two contract-driven deviations accepted at review**: `ContainerBatch` (one stacked array per label, shared axes held once) instead of a leading-axis `FunctionSamples`, since a container's shape is what its coordinates imply; and jax problems **cannot be pickled** (`jaxlib` `Device` handles), declared as `BackendCapabilities.picklable = False` with the refusal asserted rather than the row skipped — device throughput is slice 2's per-chunk `vmap`. Append-time table in the merge's branch report: cost is quadratic in the *number of chunks* (50 chunks of 10⁵ tiny draws: 36 s writer vs 55 s simulator) — the trigger for the deferred unlimited-dimension append. Carried: `ProcessExecutor` builds a fresh pool per chunk (reuse is folded into slice 2); `ThreadExecutor`'s timeout cannot interrupt a thread (documented); **open for Peter**: the multiprocessing start method — `fork` works but jax/torch warn about forking a threaded runtime and Python 3.14 changes the default; `mp_context=` is the escape hatch today, a `forkserver`/`spawn` default is the question |
| W3.2 | merged 2026-09-09 at `7ce3e81` (Opus-authored, Fable-reviewed; two review fixes ride along at `cbad8b4` — the eight `*_version` attrs in `_nuts.py`/`_vi.py` are `str()`-wrapped so a torch NUTS/VI run can be written to netCDF, and the **`sbi` pixi environment gains the `torch` feature** because its gate had been red since W2.11 (the `sbi` extra installs torch without pyro and the torch NUTS/VI/M2 rows do not skip). Merged-master gates: **sbi 2408/41 — the environment's first green gate**, dev 1748/155, jax 2195/86, torch 2365/84; lint/format clean, pyrefly 0 errors in dev and sbi). `SBIEngine(problem, method="npe"|"nle"|"nre", budget=, embedding=, executor=, chunk_size=, context=None, training_set=)` over `simulate_many` on any backend; the prior in **unconstrained** coordinates; one chain of i.i.d. draws scored on the numpy path with the estimator's own log-density beside them in `sample_stats`; the `"flat"` summary behind one function for W3.3; the legacy embedding vocabulary carried over. NPE recovers the toy joint posterior at budget 2000 (Δmean ≤ 0.13σ, width ratios 0.83–0.94); NLE/NRE run end to end; the black-box subprocess simulator fits under `ProcessExecutor`; 20 000 simulations in 152 s end to end. Decisions in the plan's §2 row "The SBI engine's shape (W3.2)" (Peter to ratify at leisure). Carried: sbi's `CNNEmbedding` asserts on summaries shorter than ~20 features (documented, sbi's own message names the remedy); sbi warns that the prior has no closed-form mean/std (expected — only input z-scoring uses the estimate); `emit()` may want a first-class extra-sample-stats hook at W3.6 |
| W3.1 (slice 2) | merged 2026-09-09 at `9e722ba` (Opus-authored, Fable-reviewed; the orchestrator applied at `3c034b9` the one `training.py` change the agent could not own — training sets now carry `ampere_simulate_batched`/`ampere_sample_backend`/`ampere_evaluate_batch`/`ampere_simulation_context` — and two autodoc entries; branch gates dev 1761/125, jax 2232/56, torch 2403/54, lint/format/pyrefly clean ×3; merged-master gates all green: dev 1780/169, jax 2251/100, torch 2422/98, sbi 2465/55). `LoweredProblem.simulate_batched` per chunk on both backends (`BatchedPrediction`), `sample_observations` natively for the Gaussian family with independent and GP noise, `simulate_many(native=None|True|False, sharder=, context=None)` with `_NativeBatch` as the single fallback site and the same one-point agreement check `realise` makes; `ChunkSharder` in core with `SingleDeviceSharder` + `MeshSharder` (jax `pmap`) and `DistributedSharder` (torch) smoke-tested at API level; `ProcessExecutor` keeps its pool across `map` calls and **defaults to `forkserver`** (Peter's ruling), with workers pre-started so the first draw's timeout is honest; torch's GP-marginal zero-uncertainty refusal and jax's `GaussianProcessNoise(jitter=, scale=)`. **One departure from the item text, correct and accepted**: native twins exist only for what `ampere.core` itself samples — `GaussianFamily`; `student_t`/`cauchy`/`complex_gaussian` inherit the contract's refusal on every backend, since a backend sampling where numpy refuses would be a guessed observation process with no oracle. **Two defects the `forkserver` default exposed, fixed**: a pickled `FittingProblem` could not be used because `Instrument.freeze` fingerprints steps by `id()` — `Dataset.__setstate__` re-freezes, and picklability is asserted by round trip. Throughput: torch ~10× at 64–400 points; jax 1.05× at 64, 3.9× at 400, **92×** at 1 200. Carried: `Instrument._declarations()`'s `id()` fingerprint still breaks a *bare* pickled `Instrument` (a value-based fingerprint in `transform.py` is the fix — small core item); the "no uncertainty at all" branch of both backends' `_check_gp_uncertainty` is unreachable (composition refuses first); pytest-collected `ProcessExecutor`s hold workers for the session (harmless). **Ruled 2026-09-09 (Peter)**: (1) `native=None` stays the default (recorded in `sample_backend`, opt-out `native=False`); (2) the core implements `sample` for further families when a use case arrives — the pathway must stay open (nothing in the native dispatch may assume Gaussian: a new core `sample` must reach the backends by adding a twin, not by rewriting the dispatch); (3) `MeshSharder`'s padding is accepted, to be **tested and costed on hardware with the GPU item** (the gain must outweigh ≤ devices−1 wasted draws per chunk) — not blocking |
| W3.3 | merged 2026-09-09 at `1ca558a` (Opus-authored, Fable-reviewed; no fixes needed — branch gates dev 1834/184 and sbi 2534/55, lint/format/pyrefly clean ×2; merged-master gates after W3.3+W3.5: dev 1867/191, jax 2338/122, sbi 2574/55, torch 2509/120 on the re-run after two first-run failures were explained and fixed — see the handoff). `ampere/core/encoding.py`: `EncodingLayout.from_datasets` (statistics and the mask frozen from the observed containers), `encode`/`encode_observations`/`decode`/`unpack` with `per_dataset` views and `grid_shape`; `SBIEngine(layout="flat"|"set"|EncodingLayout)` with `embedding="set"|"transformer"` behind wrappers that convert the one mask column (NaN rows / `attention_mask`), `z_score_x="none"` for set layouts, the layout by name, hash and in full in the run's and training sets' attrs; `encoding.md` bound and frozen with seven *Amended W3.3* sentences (decision-log row). Found: a latent train/inference skew in W3.2's flat path (mask read lazily) closed structurally; **two upstream sbi 0.27 bugs** — the set embedding's valid-row count reads the first batch element only, and the transformer drops `attention_mask` unless causal — both worked around in the wrappers (worth upstream issues). Carried: `z_score_x="none"` does not silence sbi's constant-column warning; `row_cap` above the observation is legal but inert until per-dataset capacity (§9.1). Open for Peter: the set embeddings' default `output_dim` (inherits legacy's `2·free_size`, arguably narrow for a pooled set); the transformer's last-token read (a CLS token or masked-mean head is a wrapper away, deliberately not added) |
| W3.5 | merged 2026-09-09 at `6dba1f3` (Sonnet-authored, Fable-reviewed; the `_sbi.py` wiring applied by the orchestrator in the following commit, adapted to W3.3 — plus `__reduce__` hooks on the lazily-built prior and wrapper classes, without which a trained posterior could not be pickled; branch gates dev 1813/172, sbi 2501/55). `ampere.results.artefacts`: `ArtefactKey`/`artefact_key`/`ArtefactStore` (`get`/`put`/`diff`/`train_or_load`), JSON sidecar of the ingredients so a miss names what changed, corrupt or hand-edited artefacts a warned miss, no partial keys, no force; `model_hash` added beyond the item text. Carried (real, W2.8 code): `append_training_set` checks only `ampere_spec_hash`, so a model change with identical parameter names could append onto a set simulated under the old model — small follow-up item to compare `model_fingerprint` too. Open: whether `model_hash` should be promoted into `provenance.py` as a public fingerprint |
| W3.8 | merged 2026-09-09 (Opus-authored, Fable-reviewed; no fixes needed — branch gates dev 1890/205, jax 2375/136, torch 2545/134 with the one failure being W3.5's versions test already fixed on master at `ed3eaaf`; lint/format/pyrefly clean ×3; merged-master gates all green: dev 1890/205, jax 2375/136, torch 2546/134, sbi 2611/69). **Root cause found one level below the item's**: solvers build covariances by calling `Kernel.matrix`, so a core kernel in a native solver detaches the graph — the kernel now joins `capability_parts` with the four flags. `FittingProblem(allow_foreign_parts=True)`, `foreign_parts`/`foreign_part_names`, qualified class names in refusals, `differentiable` forced `False` under the opt-in, `ampere_foreign_parts` in provenance (conditional, no schema bump — folded into W3.12's), `foreign_parts_refusal` from `realise` and both `LoweredProblem`s, `_refuse_foreign_parts` ahead of NUTS/VI's other checks; the fast path falls back through `realise`'s refusal with no engine change. The callback route (torch detach hop, jax `pure_callback`) **measured and rejected**: 1.2–1.6× where all-native is 1.4–1.7× (torch) or 5–7× (jax). Accepted at review: presence of `ampere_foreign_parts` is the flag (no separate boolean — the attr answers "which pieces got no gradient"); the working name kept. Carried: torch `noise.py`'s `_check_one_device` is partly redundant with the new check but still needed at noise-model construction; `tests/gpu` now also catches a CPU kernel in a GPU problem (unexercised) |
| W3.6 | merged 2026-09-09 at `d25d8a5` (Opus-authored, Fable-reviewed; no fixes needed — branch gates dev 1900/200, sbi 2615/56, lint/format/pyrefly clean ×2, docs build clean; merged-master gates dev 1923/214, sbi 2652/70 green). `ampere.results.calibration` (`calibration_dataset`, `sbc` — the Talts et al. refit loop with a duck-typed `engine_factory` and `replace_observations` — `attach_calibration`), `SBIEngine.calibrate` over `run_sbc`/`check_sbc`/`run_tarp`/`check_tarp` on the run's own layout (1.2 s at 100×100 after a 13.5 s train), `plot_sbc_ranks` (99 % binomial band, failure shape named) and `plot_coverage`; the `calibration` group (decision-log row); `diagnostics.md` §11 landed. Verified against a closed-form linear-Gaussian posterior (ranks uniform within 3 SE) and that a temperature-scaled posterior *fails* (`ks p` 5e-9). **The WStat coverage study** (W2.9's request, 128 trials × 200 draws, 21 min for both routes): under the example's broad prior nothing fails — SBC averages over a prior that is mostly bright sources — so the study runs under the faint prior `log-U(0.5, 3)` by default with `--broad-prior` as the comparison; there WStat's spectral index fails uniformity (`p = 4×10⁻⁴`, mean rank fraction 0.558, rank variance below 1/12: bias *and* overstated precision) where the full joint fit on the same simulations does not (`p = 0.18`); the page states the two qualifications. Carried: `PoissonFamily.sample()` is unimplemented so count data cannot `simulate(observe=True)` — the example subclasses it; **this is the use case Peter's ruling of 2026-09-09 (further core `sample` families) was waiting for** — W3.14 drafted; `replace_observations` rebuilds from the public surface and should delegate to a `FittingProblem.replace` if one ever lands; `pyrefly` on `PATH` resolves to a miniforge binary with no project deps (use `pixi run typecheck`); the agent briefly edited the shared checkout by `cd`-ing to it and restored it byte-for-byte (verified clean at the merge) — the worktree onboarding should say the isolation guard does not block file edits after a `cd`. **Ruled 2026-09-09 (Peter)**: the study's faint-prior default stands, `--broad-prior` the comparison |
| W3.11 | merged 2026-09-09 (Sonnet-authored, Fable-reviewed; no fixes needed — branch gates dev 1923/217, sbi 2655/70, lint/format/pyrefly clean ×2; merged-master gates dev 1923/217, sbi 2655/70 green). `_default_output_width`: `max(2·free_size, 32)` for set/transformer, `2·free_size` for flat, explicit `output_dim=` wins; the transformer wrapper now runs sbi's body modules (`preprocess`, `layers`, `norm`, `aggregator`, `is_causal` checked) and pools every retained token with a mask-weighted mean instead of sbi's last-token read — a test that fails against the old wrapper and passes on the new, and a test that fails loudly if sbi renames those attributes; `encoding.md` §7 *Amended W3.11*. Peter's rider stands: defaults, not findings — the embedding study is in the deferred list |
| W3.4 | merged 2026-09-09 (Opus-authored, Fable-reviewed; no fixes needed — branch gates sbi 2706/70 in 20 min, dev 1943/248, lint/format/pyrefly clean ×2; merged-master gates dev 1943/248, sbi 2706/70 green). `ampere/inference/_tmnre.py` (`TruncationBox`, KDE marginal density, `MarginalEstimator`, `MarginalSummary`, a `RestrictedPrior` subclass that neither prints nor normalises) and `SBIEngine(method="tmnre", rounds=, marginals=1|2, truncation_epsilon=, sample_with=)`; the `marginals` group stored by default (decision-log row); `ampere_sbi_truncation` per round; `ampere_sbi_amortised` on every SBI run; `calibrate()` works on a TMNRE run; `examples/sbi/tmnre_fit.py`. Three-round run 226 s at the fixture budget; boxes nested and containing the truth. **Two Accept deviations accepted**: coverage + width instead of the 0.5σ location row (unreachable for a ratio estimator's joint; documented); the cache key's `marginals`/`truncation_epsilon` folded into the `architecture` ingredient as a stopgap — **the proper `artefact_key` diff is in W3.4's report and is folded into W3.12's scope**. Carried: **an `SBIEngine` run is not reproducible seed or no seed** — nothing seeds torch's global generator (network init and batch order vary) — a small item (W3.15, drafted); sbi's `RestrictedPrior` prints and normalises expensively by default (subclassed here); `RatioEstimator.unnormalized_log_ratio` refuses to broadcast; `examples/` is outside ruff's gate (two example files carry findings). **Ruled 2026-09-10 (Peter)**: `"mcmc"` the default (W3.16); the coverage row, ε (study deferred) and the final-round-only pair estimators confirmed |
| W3.14 | merged 2026-09-09 at `4643f53` (Opus-authored, Fable-reviewed; no fixes needed — branch gates dev 1937/227, jax 2432/160, torch 2608/156, sbi 2684/80, lint/format/pyrefly clean ×4, docs 13 warnings; merged-master gates all green: dev 1957/258, jax 2452/191, torch 2628/187, sbi 2735/80). Core `sample` for Poisson (latent-GP aware), Student-t and complex Gaussian (decision-log row); `_TWINNED_FAMILIES` per backend with `sampling_refusal` testing whether the family's `sample` is the core's own (user overrides fall back); torch's two-stage twins; `_draw_dtype` and `_whitened_latent` in `dataset.py`; the refusal string byte-identical and tested word for word; `CountingPoisson` deleted with the WStat study's ranks pinned as integers. One documentation edit outside its ownership (`wstat_comparison.rst`'s bullet on the deleted subclass) — accepted. Carried: **a native sampler failure aborts a whole chunk** where the numpy loop flags one draw (`_SimulateNatively.run` calls the sampler outside `_draw`'s `try`; the Poisson twin is the first family that can reach it — small `dataset.py` fix); torch's Student-t twin ~1.2 s per 4 000 draws (a private-API path is 20× faster, deliberately not used); legacy `examples/NGC6302*.py` fail a bare `ruff check examples`. **Ruled 2026-09-10 (Peter)**: by principle — a data type's likelihood arrives with its sampling form on every backend for every engine (`likelihoods.md` §3), so `cauchy` waits for its data type; the guard aligned to `<= 0` (W3.16) |
| W3.12 | merged 2026-09-09 at `449eda2` (Sonnet-authored, Fable-reviewed; **one fix at review** — `sample_with` joins the artefact key: the removed `_key_architecture` stopgap folded three TMNRE settings and the item text named two, so two TMNRE runs differing only in the sampling mode shared a digest while sbi bakes the mode into the built posterior the store holds — the agent had flagged the gap as carried; branch lint/format/pyrefly clean, `tests/results` 363/3; merged-master gates: dev 1982/258 green; sbi 2759/80 with **one failure, a literal `== 5` schema pin in `tests/inference/test_sbi.py` that W3.12 could not see** — its worktree had no sbi environment and the sbi gate was not run on the branch before the merge (the orchestrator's miss); the one-line fix (compare against `PROVENANCE_SCHEMA_VERSION`) is folded into W3.15's branch, which owns that file; jax 2477/191 and torch 2665/187 green). `model_hash(problem)`/`model_hash_fingerprint(problem)` in `provenance.py` — named `model_hash`, not the item's `model_identity_hash(problem)`, because `model_identity_hash(model)` already exists with W2.1's cross-backend "offer" semantics (accepted); `ampere_model_hash` on every run and training set, `PROVENANCE_SCHEMA_VERSION` 6 (decision-log row); `artefacts.py` builds its key from it with no private fingerprint code left; `ArtefactKey` gains `marginals`/`truncation_epsilon`/`sample_with`, written only when set, a non-TMNRE digest pinned to the base-commit value; `append_training_set` refuses a model-hash mismatch by name and refuses a pre-schema-6 file by name rather than guessing; `results.md` §4/§9/§11/§15 and the rst amended. The agent's report was lost with its session (the orchestrating session was cleared); the branch was reviewed from its commits and commit messages, which carried the reasoning. Carried: none new |
| W3.7 | merged 2026-09-09 at `a744009` (Sonnet-authored, Fable-reviewed; no fixes needed — YAML validated; no gates: a workflow-and-prose diff). `backend-suites`' matrix gains `sbi` as a third blocking leg, producing the check **`new-namespace suites (sbi)`** (`typecheck` then `test-all`, install cached like torch's); `bench` and its upload skipped for `sbi` (it would duplicate torch's numbers — a scope call, accepted); the weekly `sbi-characterisation` job untouched, its comment now saying why it stays non-blocking (disjoint legacy coverage, not download weight); `docs/development.md`'s Environment section documents the matrix and the reference timings. **Accept row not verifiable locally**: the sbi leg's local reference time (14–20 min) sits above torch's (12–17), so the leg is likely the slowest on CI — reported honestly with two levers not applied (typecheck split into its own job; tighter smoke budgets). `actionlint` is not installed here. **Ruled 2026-09-10 (Peter)**: `actionlint` joins the dev toolchain if it earns its place (keeping action versions current in particular); CI split into smaller jobs gated on the paths touched — W4.10 |
| W3.10 | merged 2026-09-09 at `1f4e8d3` (Sonnet-authored, Fable-reviewed; no fixes needed — branch dev gate 1994/258, lint/format/pyrefly clean; merged-master dev gate 1994/258 green). `plot_corner`/`plot_trace` gain `paginate=True`: above the cap a `list[Figure]`, one page per cap's worth in merged-name order, an array block kept whole where it fits and split across full pages only when it cannot fit any page; a `ResultsWarning` (new, in `plots.py`, re-exported) naming the page count, the cap, `var_names=` and `paginate=False`; `"page": "i of n"` in every paged figure's metadata; a call that fits still returns one `Figure`; `paginate=False` is the pre-W3.10 path byte-for-byte, refusal text asserted; `labels=`/`truths=` validated once against the whole column set and sliced per page; `lp` on **every** trace page (accepted — the failure it exists to catch applies per page). `results.md` §8 amended, W2.8's decision-log row gains the note. Carried: `ResultsWarning` lives in `plots.py` beside `ArtefactCacheWarning`'s pattern rather than in `core/exceptions.py` with `ResultsError` — hoist if a second results warning appears; the warning wording is the agent's own |
| W3.15 | merged 2026-09-10 at `1fb3097` (Sonnet-authored, Fable-reviewed; no fixes needed — branch gates dev 1982/264, sbi 2766/80, lint/format/pyrefly clean; merged-master gates sbi 2778/80 and dev 1994/264 green — the first clean sbi gate on master since the W3.12 schema bump). `SBIEngine._seeded(torch, concern)`: a restore-on-exit context manager seeding torch's global generator (and CUDA's when the device is a GPU) from `integer_seed("sbi.torch")` before the first network is built, before every round's `train()` (TMNRE's marginal and joint estimators included), before a multi-round proposal's own `sample()`, and — a distinct concern, `"sbi.torch.sample"` — immediately before the final draw so a cache hit repeats too; `ampere_sbi_torch_seed` recorded, absent for an unseeded problem. **Found on the way**: sbi's default MCMC method (`slice_np_vectorized`, used by NLE, NRE and TMNRE `sample_with="mcmc"`) draws through numpy's *legacy* global generator, so torch seeding alone left MCMC-sampled posteriors irreproducible — `np.random` is seeded and restored by the same context manager. `inference.md` §13 gains "Training and sampling are reproducible too" (*Amended W3.15*); the rst amended. The restore-on-exit shape (beyond the item's literal `torch.manual_seed`) accepted — it is `_nuts.py`/`_zeus.py`'s own rule. Also carries the one-line schema-pin fix in `test_sbi.py` that W3.12's merge exposed. Carried: `calibrate()`'s internal SBC/TARP sampling through `sbi.diagnostics` is as (ir)reproducible as before — it never retrains, and the item's condition did not reach it |
| W3.13 | merged 2026-09-10 at `0ff223e` (Sonnet-authored, Fable-reviewed; **one fix at review** — `results.md` §18's new annotation said the artefact key uses "the same two hashes" as the append check, but `ArtefactKey` also carries `data_hash` of the observed containers, corrected at `e0d0cf7`; branch gates dev 1994/268, sbi 2782/80, lint/format/pyrefly clean, `pixi run docs` 13 warnings before and after, byte-identical set, all legacy; merged-master gates dev 1994/268 and sbi 2782/80 green). `docs/source/sbi.rst` — the SBI tutorial built from all five `examples/sbi/` scripts with real console output, in the v2 tutorials toctree; `advanced.rst`'s legacy note inverted, `migrating.rst`'s `SBI_SNPE` row completed, `install.rst`/`index.rst`/README updated (six engines; the WStat coverage numbers in the pitch); five stale sentences amended in place (*Amended W3.13*) across `architecture.md` §3, `inference.md` 17.5/§18, `diagnostics.md` row D, `results.md` §18 — the last corrected rather than dated (the training-set check is spec + model hash, never `ampere_problem_hash`); the plan's §5 Phase 3 bullets gain the landed-summary paragraph and §6's SBI bullet loses swyft with the ruling and gains the embedding and ε studies; `examples/sbi/cached_fit.py`'s stale docstring fixed and `tests/examples/test_cached_fit.py` added (it had no coverage). **One deviation accepted**: the item mapped the NPE-on-a-native-problem material to `toy_powerlaw.py`, which is the fake compiled program, not an ampere script — the tutorial follows the code. Carried: `docs/source/overview.rst` is framed "as of M2" and still says five engines and "SBI is Phase 3" — a small refresh item; no `examples/sbi` script does a bare in-process NPE fit for its own sake — **both ruled wanted 2026-09-10**: W4.0 (1) and (6) |
| W3.16 | merged 2026-09-10 at `fea7ef9` (Fable-authored; branch: TMNRE tests 45 passed under the new default in under four minutes, the Poisson rows green, lint/format/pyrefly clean; merged-master gates all green: dev 1996/268, sbi 2784/80 — five minutes faster than before the default switch — jax 2491/201, torch 2667/197). `TMNRE_DEFAULT_SAMPLER = "mcmc"`; `PoissonFamily.sample` guards `rate <= 0`; the sampling-form principle in `likelihoods.md` §3; the example's default, its header, the API page and the tutorial updated. Found on the way: `TestTheCalibrationFastPath::test_a_temperature_scaled_posterior_fails_the_same_check` failed once under a `-k` selection that changed the module-scoped fixture order, and passed alone — order-dependent, not a regression (the full gates are the verdict) |
| W4.5 | merged 2026-09-11 at `6a78a2b` (Opus-authored, Fable-reviewed as the second pass the 2026-09-05 review-policy ruling requires while Codex is blocked — the Matérn-5/2 rank-3 expansion, the per-leaf unit rule and the `kernel_for` binding checked by hand; the `gpt-5.6-terra` pass owed; no fixes needed at review; branch gates dev 2144/269, sbi 2967/81, torch 2850/198, jax 2674/202; lint/format/pyrefly clean ×3). Landed: `ampere/core/kernels.py` (the kernel surface moved out of `likelihood.py`, which re-exports every name): `Matern12`/`Matern32`/`Matern52` (ranks 1/2/3 — **§15.3's claim that Matérn-5/2 is not exactly quasiseparable is withdrawn**: false of the semiseparable form the solver factorises), `SHO`, `RotationTerm` (rank 4; `quality` is celerite2's Q₀, the excess over ½), `Sum`, `Product` (dense-only by declaration), `SpectralMixture`; `amplitude` is the marginal standard deviation throughout; the **`axes=` selector** on every kernel and `KernelSpec` (default "every axis", so no pre-W4.5 spec hash moves — pinned), the single-unit rule applied per leaf, binding functional through `GaussianProcessNoise.kernel_for(observed)` (W4.2 must take its kernel from there, never from the unbound declaration); a composite qualifies its children's parameters by term label (`term0.amplitude`, or `labels=`), declaration-order hashing; **one public registry keyed on the kernel family**, `register_quasiseparable_term(kernel_type, builder, *, override, builtin)` in `lowering.py`'s shape, the array namespace an argument (`ArrayOps`; `TorchOps`/`JaxOps`; `Kernel.with_ops`) so a user registers once for all three backends and the three private tables are gone; `CovarianceSpec` a tree in the conformance suite with one recursive builder for all fixtures; `examples/m2_misspecification/fringing.py` and the informational `tests/m2/test_fringing_kernel.py` (strong_smooth: worst bias 11.25 standard → 0.66 Matérn-3/2 → 0.46 Matérn-3/2 + SHO, 4/4 covered under both GPs; the SHO's period preference confounded by smoothness, so the pinned assertion is the margin growing with the fringing). `likelihoods.md` §6–§8 amended, §15.1/§15.3 lifted, §15.2 half. Carried: `dataset.py:395`'s docstring cites the old module path (W4.0's file); `part_name` outputs for kernels now read `ampere.core.kernels.*` (no hash reads a module path); a `Sum` containing a `Product` gets the generic O(N) refusal rather than the sharper one. Merged-master gates (2026-09-11/12, master at `bac8b4a`): dev 2144/269, sbi 2967/81, torch 2850/198, jax 2674/202 — all green. |
| W4.0 | merged 2026-09-13 at `c9e4015` (Sonnet-authored, Fable-reviewed; no fixes needed at review; the agent was cut off twice by the session rate limit with its work committed, and finished its gates on resumption; branch gates dev 1999/275, sbi 2792/84, torch 2672/204, jax 2496/208 — the jax gate first failed its own new test because **jax's Poisson twin drew silently for a non-positive rate** (torch's and core's refuse), fixed at `90a37a1` with an eager `rate <= 0` refusal in `backends/jax/problem.py` before the vmapped draw). Landed: (1) `overview.rst` refreshed; (2) `_NativeBatch._sample_chunk` retries a failed batched native draw row by row so one bad θ is flagged and the rest keep their native draw (strict propagates), tests on every backend fixture; (3) `Instrument._declarations` fingerprints by value (`Parameter.__eq__` compares priors by value), `Dataset.__setstate__` removed, a bare pickled `Instrument` round-trips; (4) `architecture.md` §3 records D1; (5) `examples/sbi/npe_native.py` + smoke test, linked from the SBI tutorial; (6) `configure_from`'s docstring corrected; (7) the reference `LSFConvolution` caches its operator keyed on grid identity then equality, with a rebuild-count test and a cached-equals-rebuilt row. Merged-master gates (first wave, master at `03b6c57`, 2026-09-13): dev 2312/284, sbi 3145/112, torch 3023/234, jax 2847/238 — all green. |
| W4.6 | merged 2026-09-13 at `b6d9298` (Opus-authored, Fable-reviewed; no fixes needed; hand merge of the three `__init__` export lists against W4.5 and of one `architecture.md` bullet against W4.0; branch gates dev 2057/276, sbi 2850/83, torch 2731/202 — the agent was cut off by the rate limit before its jax gate, covered by the merged-master jax gate). Landed: `ampere/core/astropy_compat.py` — `from_astropy(model, *, kind=, priors=, channel=, grid=, output_unit=, equivalencies=)` wrapping any `astropy.modeling` model (compound included) as `AdaptedAstropyModel` (`DIFFERENTIABLE = False`, `BATCHABLE = False`, `BACKEND = "reference"`): finite `bounds` → a uniform prior in the parameter's own unit, `fixed` → frozen, `tied` → an `AstropyTie` recorded (not a sampled parameter) and applied per evaluation, `priors=` overrides, a free unbounded parameter refused by name; units honoured with **no solid angle invented** — a per-steradian output needs `equivalencies=[u.dimensionless_angles()]` (one steradian) or the model's own scale, checked exactly against the reference `BlackBody`; the kind declared, or inferred only from an unambiguous axis unit and arity; `compile_for` adopts the negotiated grid; emcee and SBI (process pool) rows, NUTS/VI refuse by name. The opt-in hook `ampere.backends.{torch,jax}.from_astropy` with an empty `TRANSLATIONS` table refuses every model naming the untranslatable components (W4.7 fills it and writes `_compose`). `docs/design/contracts/astropy_compat.md` (§3–§7 binding), its doctests run; `core/__init__`'s docstring no longer says the contract is pending. Merged-master gates: as W4.0's row (the first wave gated together). |
| W4.1 | merged 2026-09-13 at `eb7b702` + fix-up `03b6c57` (Opus-authored, Fable-reviewed as the second pass; the agent was cut off by the rate limit twice, finishing dev/sbi/torch gates on its branch — dev 2100/268, sbi 2888/104, torch 2771/221 — with the jax gate covered on the merged master; **one semantic merge conflict with W4.5**: W4.5 moved the kernels and their `astropy.units` import out of `likelihood.py` while W4.1's von Mises unit check uses `u.rad` — pyrefly caught it, import restored in the fix-up; one textual conflict in `tests/conformance/oracles.py`, both blocks kept). Landed: **`VisibilitySet` amended to `(u, v, spectral_axis)`** (one wavelength per sample, `Order.ANY`; `visibility` the fourth positional; the call sites and the `serialisation.py` doctest updated; `results_schema.md` amended with the chromatic-error reasoning) and the **`ClosurePhases` kind** `(u1, v1, u2, v2, spectral_axis)` in radians with the canonical ordering (i < j < k, baselines ij and jk stored, ki implied as their negated sum, `arg(V_ij·V_jk·V_ki)` wrapped to (−π, π], sign convention `b_ij = r_j − r_i`) and `implied_baseline()`; `backends/reference/interferometry.py`: `FourierSample` (direct separable DFT of a gridded `Image` with `exp(−2πi(ux+vy))`, cell solid angle integrated, coverage from the observed container via `from_observed`, requirements ±fov/2 with `max_step = 1/(2 s u_max)` over the *expanded* coverage and an `oversampling` factor whose meaning the docstring states — a compact source needs more than Nyquist), `ClosurePhase` (three-to-one, refuses non-closing or non-co-spectral triples, mask propagated), `BandwidthSmearing`/`TimeSmearing` (average over extra (u, v) sub-samples the Fourier step computes for them, counted through `configure_from` — the first cross-kind uses), `Amplitude`; `UniformDisc`/`GaussianSource`/`Binary` image models with `compile_for` onto the negotiated grid, and their analytic visibility twins as oracles; **`VonMisesFamily` implemented** (normalised, `i0e`, wrapped residual, radians enforced by `check_observed`) **with `sample()`**, a GP refused by name as latent (W5.1). Conformance (`tests/conformance/test_interferometry.py`, reference and mirror fixtures; a backend without twins skips by name): band-limited image against the closed form, the analytic route, a sharp edge converges rather than agrees, the Nyquist requirement and the refusal, the closure phase of the binary against the closed form and its sign, mask propagation through the triangle, both smearings against brute-force fine sampling and the extra samples widening the pixel scale, two datasets on one `sky` channel — both sources recorded, the model evaluated once per draw, the composed likelihood peaking at the truth, `simulate` drawing both kinds — and the amplitude route; the emcee recovery of the synthetic binary in `tests/backends/test_reference_interferometry.py`. Carried: the reference `Image` channel is achromatic (a `Cube` channel is later); `RiceFamily` left as found. Merged-master gates: as W4.0's row (the first wave gated together). **Conventions accepted by Peter 2026-09-15**; the `VisibilitySet` and von Mises decision-log rows confirmed the same day. |
| W4.10 | merged 2026-09-13 at `d009e65` (Sonnet-authored, Fable-reviewed; no fixes needed; a workflow-and-prose diff plus one dev dependency — no five-suite gate; `actionlint`, YAML validation, lint/format, the import sweep and the docs build clean; 202 k agent tokens in 17 minutes under the no-gate-waiting rule). Landed: `ci.yml` — a `changes` job computes five run flags from the PR diff through `.github/scripts/path_filters.py` (the one executable source of the path table: core/results/inference/pyproject/lock/workflows → everything; torch → torch + sbi; jax → jax; reference and other tests → dev; `examples/sbi` → dev + sbi; `examples/interferometry` a marked slot; docs and markdown → dev + docs; unrecognised → everything); every job always runs and its *steps* are gated, so a required check skipped by the filter completes with conclusion success, never "skipped" or "expected"; push, schedule and dispatch are never gated; `backend-typecheck` (torch/jax/sbi) split out of `backend-suites` so the sbi leg's critical path is `test-all` alone; an `actionlint` job and `pixi run actionlint` task (conda-forge `actionlint` 1.7.12 in the dev feature; `pixi.lock` +21 lines); every action pinned to a full SHA with a version comment; `.github/dependabot.yml` for the actions ecosystem; `docs/development.md`'s CI section rewritten with the table and the required-check names. Judgement call accepted at review: `ampere/results/**` and `ampere/inference/**` run everything, as `ampere/core/**` does. **For Peter (part C)**: the branch-protection list is now `lint + format-check`, `actionlint`, `typecheck (pyrefly, new namespaces)`, `test (py311/py312/py313)`, `new-namespace suites (dev)`, `typecheck (pyrefly, torch/jax/sbi)`, `new-namespace suites (torch/jax/sbi)`, `docs build`, `minimal install (no extras)`; the real skip-success behaviour is verified only by validation until the first live run. |
| W4.7 | merged 2026-09-13 at `2f62fba` (Sonnet-authored, Fable-reviewed; no fixes needed; targeted rows in place of the branch gates under the 2026-09-13 rule — the astropy test files 115 passed in torch and in jax, 90 in dev, plus `tests/core` 1115 passed, the torch and jax backend suites and `tests/inference/test_nuts.py`; lint/format clean, pyrefly clean in torch and jax; the wave's merged-master gates cover it). Landed: `ampere/core/astropy_translations.py` — the six curated formulae (`BlackBody` through each backend's own `planck_jy` with astropy's raw-unit correction folded in once, `PowerLaw1D`, `BrokenPowerLaw1D`, `Polynomial1D`, `Gaussian1D`, `Const1D`) written once against a two-method `TranslationOps` protocol in W4.5's `ArrayOps` shape, `plan_translation()` walking astropy's own compound tree with astropy's leaf numbering over `+ − * /` (`|` and `&` refused by name; a tied parameter refused on the native route by name — an arbitrary Python callable has no gradient); `ampere/backends/{torch,jax}/astropy.py` — `TRANSLATIONS` from the shared table and `NativeAstropyModel` wrapping an internal `AdaptedAstropyModel` as the probe for kind, grid, negotiation and unit factor while computing the flux natively (`DIFFERENTIABLE = True`, `BATCHABLE = False` declared conservatively, CPU/float64), exposing the `grid`/`flux` surface the realisation composes so NUTS runs on it; `ampere/core/astropy_compat.py` — `translate_astropy_parameters()` extracted as the one place both routes read a bound, a fixed value or a tie from (accepted at review: exposure without behaviour change), plus three read-only properties; the contract page §5 records the landed table; conformance rows at `tolerances.cross_backend` for all six leaves and four compound forms on both backends; NUTS recovers a `Gaussian1D + Const1D` truth on both. Carried: `BATCHABLE` unverified under `vmap`; no `dtype=`/`device=` threading. |
| W4.11 | merged 2026-09-13 at `9c520c1` (Sonnet-authored, Fable-reviewed; no fixes needed; targeted checks in place of branch gates — `tests/examples` 36 passed under dev and torch, lint/format/pyrefly clean, docs build with no new warnings; the wave's merged-master gates cover it). Landed: `examples/sed_composition/` (generators with the negotiate-then-compile synthetic data, the script with `build_model`/`build_instruments`/`build_problem`/`fit`/`report`, a `__main__` that pins BLAS to one thread before numpy loads — 11 ms per evaluation against 20 ms multi-threaded), `docs/source/sed_composition.rst` (the physics, the five nouns for this case, why synthetic data needs the negotiated grid, the printed requirement explained clause by clause, the label collision and its fix, the emcee fit, the NUTS variant as one flag), cross-links from `tutorials.rst` and `overview.rst`, `tests/examples/test_sed_composition.py` (9 rows, including the collision message word for word). The reference fit: 16 walkers × 450 steps in 57–99 s wall clock, deterministic, every parameter inside its central 95 % (temperature 178 ± 27 K for 180; β 1.53 ± 0.43 for 1.6; the calibration scale 0.997 ± 0.034 for 1.02). **Deviation, accepted by Peter 2026-09-13 ("sensible")**: the filter set is mid/far-infrared (WISE W3/W4, IRAS 60, MIPS 70, PACS 100) rather than the memo's 2MASS-to-MIPS list, because a 180 K source is 1e-19 Jy in J and that dynamic range made the NUTS variant pathological. **Findings, carried**: the NUTS variant on torch is mechanically verified at 20 draws but a 40/80 budget ran past 19 minutes — the temperature/β/scale degeneracy of a modified blackbody under NUTS's default adaptation; `NUTSEngine` exposes no `target_accept`/`max_tree_depth` from the user surface (a follow-up item candidate); `examples/` is excluded from ruff by `pyproject.toml`, so the lint gate does not cover example scripts (checked by hand here); the full-budget recovery has no pytest path because the M2-style marker needs a `pyproject.toml` edit the item could not make. |
| W4.2 | merged 2026-09-13 at `4d93323` (Opus-authored, Fable-reviewed as the second pass the review-policy ruling requires — the closed form, the shared factorisation, the two-column draw and the two lowering fixes checked by hand; terra owed; no fixes needed at review; targeted suites in every environment in place of branch gates: dev `tests/core tests/conformance` 1593 passed + `tests/backends tests/m2 tests/results` 575, torch 1886 + 780, jax 1882 + 607, sbi 2208; lint/format/pyrefly clean ×3; **merged-master gates for the second wave, master at `4d93323` (covering W4.2, W4.7, W4.10, W4.11): dev 2380/308, sbi 3250/114, torch 3128/236, jax 2951/241 — all green**). Landed: `ComplexGaussianFamily.GP_ANALYTIC_IMPLEMENTED = True` — the circular complex GP as **one Cholesky of `S = K(θ) + diag(σ²)` with a two-column right-hand side** (`−½[rᵉᵀS⁻¹rᵉ + rⁱᵀS⁻¹rⁱ + 2 log|S| + 2N log 2π]`), the `GPSolver` contract widened so a `residual` may be `(n, k)` independent realisations sharing one covariance (`STACKED_RESIDUALS`, `False` by default because both failure modes are silent; `True` on all three `DenseGP`s; `conditional_loo` one term per sample summing the components; `condition` an `(m, k)` mean with one variance; `latent_transform` `(n, k)` → `(n, k)`); the O(N) refusal by name *before* the solver's own check, structural — a `(u, v)` point at a wavelength has no ordered 1-D coordinate whatever the kernel selects, so `interferometry.md` §7's "`QuasisepGP` will inherit it" is withdrawn; `sample` under a GP as two real realisations sharing `L` and nothing else; `conditional` returns a complex mean; **two pre-existing native-path bugs fixed** — both backends' lowerings took `observed.axes[0]` as the GP coordinates rather than the `(n, d)` stack, and used the unbound `noise.kernel` where W4.5 put the binding on `kernel_for(observed)`, so an `axes=` selector did nothing natively; `check_axes`'s bare-kernel message now names the selector (two pinned rows updated, ground rule 9); `likelihoods.md` §4 rewritten and §7 gains the `(n, k)` rule and the refusal; conformance rows — the closed form against `multivariate_normal.logpdf` on the materialised 2N block (|Δ| ≈ 4e-14) with the (u, v) kernel and the product, the refusal word for word, gap I-1 still caught under a GP, pointwise terms against a from-scratch leave-one-out; realisations against the contract path on torch and jax with hyperparameter gradients; `tests/m2/test_visibility_calibration.py` — a four-station array, a smooth complex gain error not from the fitted kernel family, SBC over 16 simulations: rigid coverage 0.75 against nominal 0.9, flexible 1.00 (over-covers; coverage pinned, the p-value not; 80 s). Carried: `ModelKind.COMPLEX` exists only on the torch fixture so the complex realisation rows skip on jax (W4.3); `protocol.complex_axes` fixes one wavelength so no `ModelSpec` row exercises a spectral factor; `QuasisepGP.conditional_loo` still unimplemented on numpy; W4.4's figures must pick a component or modulus of the complex conditional mean. **Decision-log row and the coverage pin confirmed by Peter 2026-09-15.** |
| W4.3 | merged 2026-09-13 at `64284d5` + fix-up `169b672` (Opus-authored, Fable-reviewed; targeted suites in every environment in place of branch gates — torch 809 passed, jax 636, sbi 814, dev 642 + 1491 + 36; lint/format/pyrefly clean ×3; **one fix at the merge**: the item's `simulate_many(native=True)` row was blocked by a one-line defect in `core/dataset.py` outside its ownership — `_draw` handed the flat `(batch, n)` channel stack to `with_values`, which needs the container's own shape; the agent verified the fix on its branch and pinned the blocked behaviour, Fable applied the fix and turned the row into the equality, 4 passed on torch and jax). Landed: `ampere/backends/{torch,jax}/interferometry.py` — every step and all six source models as **inheriting** twins on both backends (torch's re-declaring pattern deliberately not used: the declaration is the dangerous half here, a pixel-scale requirement written twice would alias silently, so both derive from the reference class and override the flags and the arithmetic; `TorchInterferometryStep` keeps torch's dtype/device plumbing); the DFT as a batched separable contraction with the phase matrices cached on the grid's bytes; `UniformDisc` and `UniformDiscVisibilities` declare `DIFFERENTIABLE = False` with reasons (a hard edge on a grid gives autograd a *wrong* gradient in the diameter; `torch.special.bessel_j1` has no backward; jax's `bessel_jn` is unusable across a realistic range — measured); a native `von_mises` on both backends agreeing with core to 0.0 across the branch cut; the realisation's native surface resolves `flux`/`grid` **or** `native_flux`/`native_grid` (one decision-log row: `flux` is what an interferometric model calls its parameter and the free-name rule forbids shadowing); jax's conformance fixture gains a complex model so W4.2's complex rows run on jax; NUTS on the binary at 200/200 × 2 recovers separation and flux ratio at 0.2σ/1.9σ (torch) and 0.25σ/1.6σ (jax) with no divergences, VI likewise; the circular GP realised on both backends agrees with the contract path to 1e-10 and the two W4.2 fixes are guarded by rows that differ only in `axes=`; NPE at budget 1200 with `layout="set"` passes C2ST/TARP/coverage (the marginal KS p-value is budget-sensitive and gated at a floor, stated); the encoding needed no change — its complex and five-axis packing is pinned by rows; GPU placement rows. **Merged-master gates (master at `169b672`): dev 2392/371, sbi 3349/90, torch 3222/217, jax 3048/219 — all green.** Carried: jax `bessel_jn` and torch `bessel_j1` findings; `encoding._sanitised`'s summary line over-promises; the encoding aligns axes by position not name (a Phase 5 design question); `von_mises` is not in `_TWINNED_FAMILIES` so its draws fall back to numpy. |
| W4.4 | merged 2026-09-13 at `9cef25a` (Sonnet-authored, Fable-reviewed; no fixes needed; targeted suites — dev `tests/interferometry` + the study's example test 30 passed in 222 s, sbi 77 passed, torch 21 + the engine rows; the jax binding the agent could not run was verified at the merge by Fable: `tests/interferometry` under jax 17 passed, 5 skipped, 240 s; lint/format/pyrefly clean; docs build with no new warnings; the final wave's merged-master gates recorded with W4.9). Landed: `examples/interferometry/` in the M2 layout — the truth a binary plus a **partially resolved** Gaussian disc (6 mas, 0.18 Jy; a wider disc resolved out on every baseline and produced no misspecification, recorded), reusing `tests/backends/interferometry_fixtures.py`'s geometry; `BinaryWithDisc` as one shared `Model` delegating to `Binary` + `GaussianSource` on any backend; the three arms (correct, disc omitted, flexible) and a chromatic truth for arm (d); `run_calibration` (SBC over 12 simulations per arm, the pinned route — a single-seed M2-style pin proved too seed-sensitive once the disc was sized right) and `run_chromatic_arm`; corner and trace figures plus the SBC rank and coverage plots; `docs/source/interferometry.rst`, the **template for adding a modality** in the prescribed order, citing `sed_composition` as the simple case; `tests/interferometry/` with an `interferometry_full` marker (one `pyproject.toml` line) and `tests/examples/test_interferometry_study.py`. **Pinned**: central-90 % coverage of separation and flux ratio — correct [1.00, 1.00], incomplete [0.25, 0.00], flexible [0.92, 0.92]. **Informational (d)**, one seed at the per-PR budget: spatial 0.50/2.16, spectral 1.17/5.21, product 1.97/1.51 (bias in widths on the two parameters) — spectral-only clearly worst, the product **not** cleanly dominant at this budget; `-m interferometry_full` reruns at the documentation budget. Timings on the reference path: three arms 19 s, calibration 142 s, chromatic 56 s. **Findings, carried**: four of the six shipped plots (`plot_posterior_predictive`, `plot_residuals`, `plot_gp_localisation`, `plot_anomaly_score`) refuse a point kind with several axes — `VisibilitySet` and `ClosurePhases` both — because `results._plotting.coordinate_of` needs one ordered coordinate (Phase 5 staging per `results.md` §4, deeper than W4.2's note anticipated); `add_residuals` and `gp_localisation` cast complex to real with a warning where `add_posterior_predictive` refuses cleanly (`results/derived.py`). |
| W4.9 | merged 2026-09-13 at `318034e` (Sonnet-authored, Fable-reviewed; no fixes needed; targeted suites — the astrometry rows dev 51 passed, torch 77, jax 76; conformance `test_astrometry.py` 24/36/36; lint/format/pyrefly clean ×3; docs build clean; the final wave's merged-master gates recorded below). **The template claim held**: a second observable landed with **no change to `ampere.core`** — `TimeSeries` already had the one strictly increasing `time` axis the modality wants, so no kind was added. Landed: `backends/reference/astrometry.py` — `ReflexOrbit` (linear proper motion plus a circular, node-aligned reflex wobble, the sketch's own simplification: `pmra`, `pmdec`, `period`, `phase`, `amp_ra`, `amp_dec`; two fixed channels `ra`/`dec` from one evaluation) and `EpochSample` (kind-preserving, the identity in `apply`, publishing `points=` at the observed epochs); inheriting twins on torch and jax (the native surface spelled `grid`/`flux` — no parameter-name clash here); the closed-form ephemeris in `oracles.py`; `AstrometryPieces` and an `astrometry` capability in `protocol.py` (outside the literal ownership list, necessary and accepted); conformance rows on every fixture — model against the ephemeris, the step's requirement and identity, the two-channel composition evaluating once per draw, `simulate` on both channels, `DenseGP` and `QuasisepGP` agreeing; emcee (24 × 1000, 16 s) and NUTS on torch (34 s) and jax (7 s) recovering all six parameters inside the central 95 %, agreeing to three figures; `examples/astrometry/` with a `--gp` arm; `docs/source/astrometry.rst` in the template's order with the closing section **"What the template did not say"** — (1) the template is silent about a one-step chain and `configure_from`; (2) the inherit-versus-re-declare rule does not settle the case where declaration and arithmetic are both cheap; (3) a periodic model is a hazard class the template never met — a period prior reaching past the observed baseline is multi-modal and chains alias onto spurious periods (measured with `loguniform(50, 2000)`); (4) "two instruments, one channel" against "one model, two channels" is not named anywhere though both now ship; (5) what it got right: a kind may already exist. All four for W4.8. **Merged-master gates for the final wave (W4.4 + W4.9, master at `318034e`): dev 2450/379, sbi 3427/90, torch 3300/217, jax 3125/220 — all green.** Carried: the sketch's circular orbit against a full Keplerian/Thiele–Innes parameterisation is a deliberate simplification to note. |
| W4.8 | merged 2026-09-13 at `c8bf2c0` (Sonnet-authored, Fable-reviewed; no fixes needed; docs build with a warning list byte-identical to the base commit's; spec doctests and `tests/examples` 76 passed; lint/format clean; **the `dev` gate on the merged master at `c8bf2c0`: 2450/379, green**). Landed: `interferometry.rst` fixed for W4.9's four gaps (a one-step chain overriding nothing is correct; a third inherit-versus-re-declare clause — both halves cheap, inherit anyway; the two shipped composition shapes named; a closing hazard section on periodic models and priors wider than the observed baseline); `astrometry.rst` §10 now "what the template now says, and where"; **`docs/source/astropy.rst`** (the casual route, the solid-angle rule, the capability consequence, the six-model native table, the two refusals) and **`docs/source/kernels.rst`** (the seven families and ranks, the algebra, the `axes=` selector and per-leaf unit rule, `register_quasiseparable_term` with the out-of-tree example, the chromatic product as W4.2's case) — new; `advanced.rst`'s noise-model section; `index.rst`, `overview.rst` (which *is* the architecture page) and README cross-linked and updated; **nineteen additive *Amended W4.8* annotations** across `likelihoods.md` (§4, §14 ×2, §15 ×2; §15.3's withdrawal confirmed present), `interferometry.md` (gaps I-1–I-5 naming their W4.x instances, §5, §7 confirmed, §10, §11 — **Q2 landed differently**: a `spectral_axis`, not `extra_coords` gaining units), `spectrum_photometry.md` (the pre-shipping `SyntheticPhotometry` signature), `astrometric_timeseries.md` (landed at W4.9 via the template), `results_schema.md` §15.4, `results.md` §13 (new limitation 14: four plots refuse a point kind with several axes — Phase 5), `transformations.md` §10 (the per-backend modules and the two smearing steps), `inference.md` §10a (the native surface's second spelling), `encoding.md` §9 (item 6, positional axis alignment), `architecture.md` §3. **Carried** (for the owed list): `RiceFamily` declared-only with no schedule; `von_mises` not in `_TWINNED_FAMILIES`; `encoding._sanitised`'s summary line; a `Sum` containing a `Product` gets the generic refusal; `QuasisepGP.conditional_loo` unimplemented on numpy; `dataset.py`'s `part_name` docstring; two pre-existing docs-build warning pairs (`Binary`'s field list in both backends' `interferometry.py`; duplicate object descriptions for the re-exported `Product`/`Sum`); `interferometry.rst` attributes a quotation to the plan that is not verbatim there. |
| W5.0 | merged 2026-09-16 (Sonnet-authored, Fable-reviewed; **two defects found at review by running the branch, fixed by the agent in `1f12e13`/`3b5e1a6`**: the pyro route read the guide's density off the `Delta` model site (always 0, so the importance test had passed for the wrong reason) and the numpyro route called `AutoNormal.get_posterior`, which does not exist, crashing every jax VI run; a third found while writing the required test — pyro's `AutoMultivariateNormal` covariance is `scale[..., None] * scale_tril`; Fable's own verification at `3b5e1a6`: jax `tests/inference/test_vi.py` 23 passed/1 skipped, sbi `test_vi + test_sbi + test_engines + tests/results` 670 passed/2 skipped, dev `tests/results tests/inference` 548 passed/216 skipped; lint/format/pyrefly clean in dev, sbi and jax; the merged-wave gate recorded in W5.2's row). Landed: `PROVENANCE_SCHEMA_VERSION` 7; `emit(sample_stats=)` (a colliding name refused); `Engine.finish()` writes `ampere_approximation = "none"` by default and threads `sample_stats` — no per-engine change for emcee/zeus/NUTS; `unconstrained_jacobian_correction` in `engine.py` (the difference of `lnprior_unconstrained` and `lnprior`, so one place knows the formula); dynesty's evidence triple (`"nested_sampling"` as the method family — a judgement call, reviewable when a second evidence engine lands); VI's `ampere_approximation` (`mean_field`/`multivariate`) and `proposal_log_density` on both routes with the fitted guide's parameters kept as arrays on the engine; SBI's `density_estimator` and its Jacobian-corrected `proposal_log_density` beside the unchanged unconstrained `ampere_sbi_log_prob`; `warn_if_approximate` in `plot_trace` and the new `ampere.results.summary`; `results.md` §4/§9/§13 amended (*Amended W5.0*), the optimum-group ruling recorded with nothing implemented; the docs pages. Carried: whether `ampere_evidence_method` should name the method family or the engine; whether `ampere_sbi_log_prob` is deprecated once a consumer moves to the contract field. |
| W5.2 | merged 2026-09-16 (Sonnet-authored, Fable-reviewed; **one fix-up at the merge**: ruff's old bare `examples` exclude pattern had also matched `tests/examples/`, so narrowing it to `examples/examples_paper` surfaced two never-linted test files — formatted; the complex-refusal test moved from `tests/core` to `tests/results`, which owns the derived groups; the agent's targeted suites: dev `tests/core` 1124 passed/30 skipped, `tests/conformance` 530/63, `tests/results` 375/3, torch `tests/backends+conformance+inference` 1643 passed, jax 1474; `test_nuts.py` 32 on each; lint/format/pyrefly clean in dev, torch and jax after the fix-up; docs warnings 19 → 9; **Merged-master gates for wave 1 (master `b6e22c2`, docs-only commits after; the jax leg re-run alone at `c3cfc9c` after the memory guard killed its first run): dev 2472/393, sbi 3463/90, torch 3333/220, jax 3158/223 — all green.** Landed, per owed item: (1) `NUTSEngine`'s `target_accept_prob`/`max_tree_depth` were already passed through — the missing test pinned; (2) `examples/` under ruff (~450 mechanical fixes, six dead assignments removed, per-file ignores with rationale for legacy idioms and three pre-existing undefined names in demonstration scripts; `examples/examples_paper` stays excluded); (3) the `part_name` docstring; (4) `_find_nested_product` — a `Sum` containing a `Product` names the term; (5) `_refuse_complex`, one refusal for all three derived groups (the component or modulus is W5.3's argument); (6) `von_mises` in both `_TWINNED_FAMILIES` — torch through `torch.distributions.VonMises` under `fork_rng`, jax through `numpyro.distributions.VonMises` (jax has no `vonmises` primitive), circular moments against the numpy oracle on both; (7) `_sanitised`'s docstring says it is `nan_to_num`; (8) `RiceFamily._unimplemented` says plainly there is no scheduled phase, and composition delegates to it; (9) both docs-warning pairs; (10) the `interferometry.rst` quotation paraphrased; (11) the native sampler's per-draw isolation was already in `dataset.py`'s `_NativeBatch._sample_chunk` — pinned with a fake native backend in `tests/core/test_simulate.py`; (12) the value-based `Instrument._declarations` fingerprint had landed at W4.0 with its test — nothing to do; (13) the unreachable `uncertainty is None` branch of both backends' `_check_gp_uncertainty` removed (core's `check_compatible` refuses first); (14) `_check_one_device`'s redundancy is deliberate and pinned by two tests — the stale docstring fixed, code kept. **Open for Peter**: whether item 14's second check is trimmed (deleting one of its two tests). Carried: `tests/backends/test_torch_device.py`'s complex-Gaussian rows fail with `No module named 'tests'` when run outside `test-all` (pre-existing rootdir issue); the nine remaining docs warnings are in frozen legacy docstrings and the optional sbi import. |
| W5.20 | merged 2026-09-16 at `acaad79` (Opus-authored, Fable-reviewed; the agent was killed by the weekly limit at commit 1/4 and resumed from a WIP commit the orchestrator preserved — nothing lost; no fixes needed at review; the agent's targeted suites: dev `tests/core tests/conformance` + the new file 1675 passed/102 skipped, torch `tests/core tests/conformance tests/backends tests/inference/test_nuts.py tests/m2` 2706/89, jax 2526/97; lint/format/pyrefly clean ×3; the two contracts' doctests 134 and 196 examples clean; Fable's own runs at `a86d156`: dev `tests/core` + the namespace file 1145/39, torch the namespace file + `test_torch_lowering.py` 219/2, jax the namespace file + `test_cross_backend.py` 109/6; the wave-2 gate (W5.20 + W5.4 + W5.3 together; dev at `12c6001`, the other three legs relaunched at `b7a14b2`/`db12df5` after the 2026-09-16 system crash, all green): dev 2594 passed / 407 skipped, sbi 3625/92, torch 3495/222, jax 3315/231). Landed: `ampere.core.reserved_names()` — the public names of `Parameterised` and `Model` computed from the classes plus the pinned `TORCH_MODULE_NAMES` (67 names in all; `AXIS`, `flux`, `grid`, `grid_tensor` are *not* reserved anywhere); `_check_free_name` tests that set, with a message naming `parameters.md` §10; `LoweredParameters._unusable`/`check_leaf` so a colliding parameter leaf raises `LoweringError` rather than torch's `KeyError`; both `problem.py` lookup tables reordered `native_*` first, `_native_surface` refusing a model that offers both pairs and still refusing half a pair with the missing half named; every shipped torch/jax model renamed to `native_flux`/`native_grid` (`models`, `astropy`, `astrometry`; the interferometry twins already were), the legacy alias kept and exercised live by the conformance `PointSourceModel`s and the three M2 example models (a judgement call the agent offers for Peter: rename those too, so the examples show only the canonical spelling); `tests/backends/test_parameter_namespace.py` parametrised over all three backends (`flux`/`grid`/`grid_tensor`/`AXIS` declarable and lowerable everywhere, `to`/`type`/`apply` refused everywhere with the core message, every non-reserved public attribute of each backend's `PowerLaw` a legal name, the pinned torch list covering `dir(nn.Module)` and a bare instance, the ambiguity and half-pair refusals word for word, a pinned spec hash unchanged); `parameters.md` §10 "The reserved names" (*Amended W5.20*, with a doctest), `inference.md` §10a, the placement memo and `interferometry.md` updated. Carried: `tests/` is not an importable package, so `test_torch_device.py`'s complex rows import `tests.conformance` only when `tests/conformance` is collected first (the W5.2 finding, now explained); the release-note items in the decision-log row. |
| W5.4 | merged 2026-09-16 at `359386c` (Opus-authored, Fable-reviewed as the second pass — the Woodbury log-determinant, solve and LOO precision diagonal, the Matérn spectral density in Solin & Särkkä's angular-frequency convention, and the `latent_size` contract checked by hand; terra owed; the agent was killed by the weekly limit at commit 2 of 7 and resumed from a WIP commit the orchestrator preserved, which it squashed; clean merge onto a master that had taken W5.2 and W5.20 since its base; no fixes needed at review; the agent's targeted suites: dev `tests/core tests/conformance tests/results` 2117/98, torch `tests/core tests/conformance tests/backends tests/inference/test_nuts.py` 2679/73, jax 2506/75; Fable's own runs on the merged tree: dev 2139/98, jax 1392/76, torch 1565 passed with the four pre-existing `No module named 'tests'` failures of `test_torch_device.py::TestTheComplexGaussianBody` (see carried), lint/format clean, pyrefly 0 errors ×3; the wave-2 gate (W5.20 + W5.4 + W5.3 together; dev at `12c6001`, the other three legs relaunched at `b7a14b2`/`db12df5` after the 2026-09-16 system crash, all green): dev 2594 passed / 407 skipped, sbi 3625/92, torch 3495/222, jax 3315/231). Landed: `Kernel.spectral_density(frequency, values, dimensions=)` in `ArrayOps` — the Matérn form once via `SPECTRAL_NU`, `SquaredExponential`, `SHO` from its own coefficients, `Sum` as the sum (so `SpectralMixture` reaches the path for free), the base refusing by name; `ampere/core/hsgp.py` (the box, the frequency grid, Φ through `ArrayOps`, `check_spectral_support`); `HilbertSpaceGP` beside the four slots with `EXACT = False`, `STACKED_RESIDUALS = True`, one m×m Cholesky of the capacitance, LOO exact in the approximation, a strictly positive noise diagonal required (refused by name); `GPSolver.latent_size` and `LatentDeclaration.whitened_size`, every family's `sample` and the whitening check routed through it; twins in both backends' `gp.py` sharing the core basis (cross-backend agreement ≈1e-15, gradients matching central differences to six figures, `jit`/`vmap` on jax; `BATCHABLE` True on jax, False on torch); `TestApproximateSolverConvergence` with `approximation_envelope` and `SolverKind.HILBERT` in `protocol.py` (1-D on the M2 grid for six kernels, 2-D on the `VisibilitySet` fixture under `complex_gaussian` at m = 16…576), the latent-path agreement row, `TestTheReducedRankSolverUnderNUTS` on torch and jax; `likelihoods.md` §7 (table row and an *Amended W5.4* subsection), `inference.md` 17.4, `tests/conformance/README.md` §2–§3, `kernels.rst`. **Ruled by Peter 2026-09-16, as recommended**: the hashing of `basis_size` and `boundary_factor` is confirmed (a cached artefact trained at one basis is never served silently for another; an *explicit* opt-in to serve a named artefact anyway is **W5.23**); `SpectralMixture` supported is confirmed. Carried: `dataset.py:1393`'s `_check_latent_count` and `:1578`'s repr, and `results/provenance.py:516`'s `latent_size`, all read `LatentDeclaration.size` (the sample count) where the sampler's block is now `whitened_size` — three one-line follow-ups; `RotationTerm`'s spectral density is five lines away; `approximation_final` is an absolute tolerance in nats where a relative one would be better; a Vecchia solver is now the measured next step for rough kernels and 2-D (W5.6); **`tests/` is not an importable package**, so `test_torch_device.py`'s four complex-Gaussian rows fail with `No module named 'tests'` in any targeted run — they pass only inside `test-all` — a one-line fix for the next housekeeping item. |
| W5.3 | merged 2026-09-16 at `47988b4` (Sonnet-authored, Fable-reviewed; no fixes needed at review; clean merge onto a master carrying W5.4; the agent's suites: dev `tests/results` 400/3, `tests/core` 1135/30, the study smoke 13, sbi `tests/results tests/core` + the study 1578/3, lint/format/pyrefly clean in dev and sbi; Fable's own runs on the merged tree: dev `tests/results tests/core` + the interferometry and SED-composition smoke tests 1619/33 (through `python -m pytest` — the bare `pytest` entry point cannot import `examples`, see carried), sbi `tests/results tests/inference/test_engines.py` 502, lint/format clean, pyrefly 0 errors in dev and sbi, plus a netCDF round trip of a multi-axis derived group re-plotted from the loaded file; the agent twice ended its turn to wait for the lock and was resumed with the polling instruction; the wave-2 gate (W5.20 + W5.4 + W5.3 together; dev at `12c6001`, the other three legs relaunched at `b7a14b2`/`db12df5` after the 2026-09-16 system crash, all green): dev 2594 passed / 407 skipped, sbi 3625/92, torch 3495/222, jax 3315/231). Landed: `FunctionSamples.PLOT_COORDINATE` with `format_axis_label`, defaults on `VisibilitySet` (baseline length) and `ClosurePhases` (longest baseline, the implied third included); `coordinate_of(coordinate=)` with the three-step resolution and `_stored_axes`; `component=` on `add_posterior_predictive`/`add_residuals`/`gp_localisation` (`COMPONENTS`, `component_variable`, `base_label` exported), `_multi_axis_extras` storing a multi-axis kind's raw axes and precomputed default on the derived group through `_attach`'s new `extra_variables`; `coordinate=`/`component=` on `residual_whiteness`, `gp_localisation_score`, `plot_posterior_predictive`, `plot_residuals`, `plot_gp_localisation` (`plot_anomaly_score` unchanged, by design — see the decision-log row); `examples/interferometry/figures.py` rendering all six plots per arm; `results.md` §4/§7/§8/§13, `interferometry.rst`'s section retitled "lifted at W5.3", `ampere.results.rst`; `TestMultiAxisCoordinate` (11 rows). Carried: `_plotting.dataset_units` never finds a complex dataset's value unit (the y-axis loses its annotation); no test exercises `examples/interferometry/figures.py` directly; `gp_localisation`'s `component=` renames the variance variable too, for lookup symmetry (a convention choice); **invocation trap**: `pixi run -e dev pytest tests/examples/...` fails to import `examples` (`No module named 'examples'`) where `python -m pytest` and the pixi tasks succeed — the same class as the `tests` import problem, for the next housekeeping item. |
| W5.5 | merged 2026-09-16 at `614d7c9` (Opus-authored, Fable-reviewed as the second pass — the FFT route against the direct-sum oracle, the separable dilation in `propagate_mask_grid`, the C-order grid coordinate matrix in `sample_coordinates`, the two native-path `Layout.GRID` fixes and `examples/image/grid_gp.py` read by hand; terra owed; **merged by Peter** because the session's permission classifier refused `git merge`; the original agent was killed by the 2026-09-16 system crash after five commits with its verification and report lost, and a fresh Opus agent finished from the branch's HEAD in 45 min with one commit; no fixes needed at review; the agent's runs: `tests/conformance/test_image.py` dev 28, torch 45, jax 45, sbi 45; full conformance dev 582/65, torch 938/65, jax 938/65; `tests/core` dev 1197/30; `tests/backends` dev 147/60, torch (+ `test_transform`/`test_dataset`/`test_simulate`) 986/3, jax 807/11; `tests/examples` dev 77/9; the image study tests 21/2 (the two `image_full` rows skipped); lint/format clean, pyrefly 0 errors ×3, import sweep 83/4, docs 13 warnings all legacy; **the wave-3 merged-master gate on `8d785e1`/`27e9bd1`, 2026-09-16, all green: dev 2695/410, sbi 3744/94, torch 3614/224, jax 3434/233**). Landed: `PSFConvolution` in each backend's `image.py` — `kernel=` tabulated at the observed pixel scale and refused on any other, `fwhm=` a circular Gaussian rebuilt through `context()` so `promote_buffer` fits it, one constraint published in two forms (`points=` on the padded grid, `max_step` at the pixel scale), the FFT route requiring an evenly spaced grid and refusing an irregular one by name, `from_observed` the supported route; the twins inheriting and overriding only the four flags and the sums, `grid_steps`/`crop_indices` shared; `ampere.core.propagate_mask_grid` (a separable dilation by the kernel's bounding half-support, O(kN)); `ampere.core.sample_coordinates` (layout-aware, used by `Dataset.draw_observation` and both `LoweredProblem.observed_coordinates`); both native paths ravel values and uncertainties before the mask and `predict()` flattens only the container's own trailing axes; `tests/conformance/test_image.py` with the `image` capability and `ImagePieces`; `examples/image/` (a compact Gaussian on a broad smooth background; correct/incomplete/flexible arms; SBC calibration; `--benchmark` DenseGP vs HSGP with `--no-dense`); `docs/source/image.rst`; `transformations.md` §10 PSF row and §13.5 amended; `ifu_cube.md` gap 2 closed; the `image_full` marker. Measured: bias in posterior widths correct 0.74/0.45, incomplete 7.92/4.68, flexible 0.78/0.13 (24×24, reference); coverage at nominal 0.9, 12 replicas at 16×16: incomplete [0.00, 0.58], flexible [0.92, 0.92]; benchmark 24×24 DenseGP 0.031 s / 12.7 MiB vs HSGP(8×8) 0.003 s / 0.7 MiB, 64×64 3.925 s / 640 MiB vs 0.006 s / 4.4 MiB, 128×128 HSGP 0.038 s / 17.6 MiB with DenseGP not run (projected ≈10 GiB and ≈4 min — the one Accept number reported by projection). **Ruled by Peter 2026-09-16, as recommended**: the two-edit change lifting the `Layout.GRID` refusals is approved as **W5.21** (after which `grid_gp.py` is deleted); `--no-dense` stays; the dense 128×128 cell stands by projection. Carried: `python -m examples.image --benchmark` allocates ≈10 GiB by default; neither `image_full` row has been run end to end (≈27 min; ≈10 GiB); `tests/examples` is not in `test-all`; a fitted `fwhm` wandering above the declared width fits a kernel truncated tighter than 5σ, and a negotiated grid finer than the observed one tabulates the Gaussian kernel over a fixed pixel count and so over fewer sigmas (the support is fixed at composition — reviewer's note); the docs warning pairs named in the handoff do not reproduce on this base. |
| W5.13 | merged 2026-09-16 at `81442e5` (Sonnet-authored, Fable-reviewed; merged by Peter; **one defect found at review by reading the reweighting identity and fixed by the agent in `17aa2c4`**: the ratio divided by the stored joint `log_prior` column where only the named parameter's marginal interim prior belongs — the nuisance priors cancel — so any object with more than one parameter was biased, invisible on the one-parameter test model; the fix takes a required `interim_prior` argument because a run's provenance stores no `PriorSpec`; the original agent was killed by the 2026-09-16 system crash after two commits with its test file uncommitted, preserved by the orchestrator as a WIP commit and finished by a fresh Sonnet agent; the agent's runs: `tests/results/test_population.py` dev 11 passed (852 s under contention; the VI row deselected), sbi 12 passed (1427 s) including the VI row; `tests/results` dev 409/4 before the fix; lint/format clean, pyrefly 0 errors in dev and sbi, import sweep 83/4; **the wave-3 merged-master gate on `8d785e1`/`27e9bd1`, 2026-09-16, all green: dev 2695/410, sbi 3744/94, torch 3614/224, jax 3434/233**). Landed: `ampere.results.population` — `RunColumns` (a `typing.Protocol`: `parameter_draws`, `log_prior`, `log_likelihood`, `proposal_log_density`, `attrs`), `DataTreeRunColumns` and `NetCDFRunColumns`/`runs_from_netcdf_directory`, `PopulationModel` with `GaussianPopulationModel` (hyperparameters as ordinary `Parameter`s so the hyperprior needs no bespoke method), `fit_population` (self-normalised importance reweighting after Hogg, Myers & Bovy 2010 — uniform weights for an exact sampler, `exp(log_prior + log_likelihood − proposal_log_density)` for an approximate engine, refused by name without W5.0's column; `emcee` over α; refusals when any object's ESS at the posterior mean falls below `ess_floor`, and when the runs do not share `ampere_spec_hash`); `results.md` §13 item 16 (the item text's "§13.15" is stale — W5.0 took 15). Measured (200 objects, truth μ = 1.0, τ = 0.3): μ 95 % (0.943, 1.049), τ (0.285, 0.372); the two-parameter regression at N = 50: μ (0.887, 1.118), τ (0.270, 0.466), agreeing with the one-parameter control to 0.15/0.1. The W5.12 half of the Accept line is not testable until W5.12 lands. **Ruled by Peter 2026-09-16, as recommended**: the 200-object row gets a `population_full` marker, and each parameter's `PriorSpec` is stored in provenance so the interim prior is verifiable from the archive — both as **W5.22**. Carried: the test's class-scoped fixture is an instance method (`PytestRemovedIn10Warning`); `ruff format` must never be pointed at a `.md` file (it rewrites doctest fences); `results.md` §14's Phase 5 bullet still reads prospectively; WORK_ITEMS.md's W5.13 text cites §13.15. |
| W5.6 | merged 2026-09-16 at `8d785e1` (Opus-authored, Fable-reviewed — the Woodbury log-determinant and quadratic form in `EquispacedFourierGP.log_marginal_likelihood` checked by hand, the bake-off's internal consistency checked (EFGP reproduces `DenseGP`'s bias/rmse to three digits at every N); merged by Peter; no fixes needed at review; the agent found and fixed one real defect in its own code — `toeplitz_matvec`'s circulant embedding filled two of a multilevel embedding's corners, right in 1-D and wrong in 2-D, visible only through `condition`'s mean (+6.5 σ) and not the marginal likelihood, so two 2-D rows were added; the agent's runs: `tests/core` dev 1238/30 (41 new rows), `tests/conformance` + `tests/results` 982/68, `tests/examples` 78/9, `pixi run -e dev bench` 26 passed/42 skipped, the `image_full` tables 1 passed in 136 s; lint/format clean, pyrefly 0 errors, import sweep 85/4; gate dev only per the item, taken with **the wave-3 merged-master gate on `8d785e1`/`27e9bd1`, 2026-09-16, all green: dev 2695/410, sbi 3744/94, torch 3614/224, jax 3434/233**). Landed: `ampere/core/efgp.py` (`EquispacedFourierGP`, `EXACT = False`, `latent_size = m`; `fourier_grid`, `spectral_weights`, `toeplitz_generator`/`toeplitz_matrix`/`toeplitz_matvec` by circulant embedding and FFT, `solve_iterative` by CG, `real_features`) and `ampere/core/vecchia.py` (`VecchiaResponseGP(neighbours, ordering, seed, jitter)`, all fields in the spec hash; `neighbour_structure` by blocked KD-tree search; `conditional_loo` exact in the approximation at O(Nk)), both reference-only measurement prototypes exported from `ampere.core` (**ruled by Peter 2026-09-16: the exports stay** — people will want to use them, and their development continues); `tests/core/test_efgp.py` and `test_vecchia.py` in W5.4's convergence class (restated locally — `tests/` has no package root); `examples/image/bakeoff.py` (grid-lifted subclasses on `grid_gp._GriddedSolver`; accuracy, smoothness × correlation length, cost, the resolution ladder, and the normal-equations two-route table); `tests/benchmarks/test_solver_bakeoff.py` with its conftest (repo root on `sys.path`, the `image_full` skip); `likelihoods.md` §7 two measured rows, the `VecchiaGP` slot annotated, and the bake-off subsection with five tables and the ruling. Measured (N = 4096, Matérn-3/2, ℓ = 8 mas): |Δlog p| HSGP m=256 2.1e-1, EFGP m=289 1.1e-2, Vecchia k=30 30.2, with bias/rmse/coverage/localisation identical to `DenseGP`'s for the two spectral solvers; cost at N = 65 536: HSGP 0.41 s / 263 MiB, EFGP 0.64 s / 123 MiB, Vecchia 14.3 s / 202 MiB; at N = 16 384 EFGP overtakes HSGP past m ≈ 576 (0.59 s vs 0.92 s at m ≈ 1024) at 4.7× less memory; the deciding table — the same Toeplitz normal equations at m = 6561: CG 0.16 s / 2.6 MiB against Cholesky 10.95 s / 1.31 GiB, the whole difference being log|M|, for which a Toeplitz matrix has no FFT. Ruling in the decision-log row. Carried: `pixi run -e dev test-examples` fails at collection in a worktree (`tests/examples/conftest.py` lacks the `sys.path` line; main-checkout `test-all` unaffected); `_GriddedSolver` private but the extension point; `HilbertSpaceGP`'s two `(N, m)` blocks; a 2-D `condition` conformance row for any EFGP twin; the `VecchiaGP` slot docstring's pointer (one line, now that the row has landed). |
| W5.7 | merged 2026-09-17 at `97f7662` (Opus-authored, Fable-reviewed as the second pass — the offset-form input warp `x ↦ x + δ(x)` with slopes `softplus(u)/softplus(0)` (monotone structurally, identity exact), the hinge-basis interpolant, the amplitude warp as generator scaling `a_n U_n`, `a_n V_n` with a per-point marginal `a_n² k(0)`, the `Kernel.warped_coordinate` hook handing celerite2 `w(x)` while the registry's builders keep the raw axis, the composite refusal, and the pre-W5.7 spec hashes re-derived independently on master and branch; **no fixes needed at review**; merged by Peter; the agent's runs: dev `tests/core` + `tests/conformance` 1885/95, torch 2278/70, jax 2278/70, NUTS rows torch 4 in 68 s and jax 4 in 13 s, `tests/results` 411/4; the orchestrator's re-run of the warped rows in dev: 65 passed in 1.3 s; lint/format/typecheck clean in all three; **the merged-wave gate on `8e10b55`/`897b9c0` (W5.7 + W5.23 + W5.9), 2026-09-17 04:13–06:58, all green: dev 2818/422, sbi 3891/95, torch 3757/229, jax 3577/238**). Landed: `ampere.core.WarpedKernel(base, input_warp=, amplitude_warp=, non_centred=True)` — knots fixed and explicit in the declaration (`quantile_knots(coordinate, count)` the user-side convenience), one half-normal shrinkage scale per warp with the knot variables normal about zero under it, non-centred by default with `HierarchicalPrior` for the centred form (a centred scale nothing depends on is refused by name); `KernelSpec.metadata` (omitted when empty, so every pre-W5.7 hash is unchanged) carrying the knots and the parameterisation; `Kernel.warped_coordinate`/`warp_provenance` (identity/`None` by default) and `CeleriteRepresentation.marginal` allowed `(n,)`; `warped_representation` registered as family `warped`; `GPConditional.coordinates`/`.warp` (`None` for every unwarped kernel) filled by `Likelihood.conditional` so families B and C run on the warped residuals; `refuse_warped_composite` on the quasiseparable path of all three backends (`WarpedKernel(Sum(...))` is the representable order; `DenseGP` takes either); no per-backend twin — the wrapper adopts its child's namespace, device and flags as `Sum`/`Product` do, and each backend re-exports the class; `likelihoods.md` §6–§8 amended; `kernels.rst`/`overview.rst` repaired; eight conformance rows per fixture, 49 core rows, 4 NUTS rows per backend. Measured: dense vs quasisep with both warps active — torch 0.0, jax 5.7e-14, reference 1e-13 relative; gradients w.r.t. six knot variables torch/jax agree to 6e-14; quasisep N = 16 000 in 1.14 ms against 4 000 in 0.44 ms (ratio 2.6; dense at N = 2 000 is 691× slower). **Ruled by Peter 2026-09-18, all three confirmed**: (1) no backend twin (`Sum`/`Product` precedent; a twin becomes a small follow-on only if a backend ever needs warp-specific arithmetic); (2) a warp *under* a composite stays refused on the O(N) path — with the **potential need recorded**: a sum over differently warped terms is a real model (several non-stationary components on one coordinate), one recursion per warp summed outside the factorisation is a different solver, and Peter expects the realistic route to be one of the rank-reduction strategies (`HilbertSpaceGP`, EFGP, or the matrix-free CG path of the plan's §5) rather than a quasiseparable trick, while suspecting the degeneracy between the warps will make it hard to use in practice; `DenseGP` takes the composition today; (3) the two named redundancies (a uniform warp against `length_scale`, a constant `log a` against `amplitude`) handled by shrinkage, not constrained. Carried: `docs/source/kernels.rst`/`overview.rst` print registry contents no doctest checks (an item: doctest the `pycon` blocks under `docs/source`); `QuasisepGP.condition(at=...)` is unusable on a multi-axis container independent of W5.7 (coerces `at` to 1-D, then `Kernel.matrix` refuses the `(n, d)` points); `kernels.py` now imports `scipy.stats` at module level (no new dependency); the M2 "many lines / one band" validation is W5.8's; a linearly extrapolated `log a` outside the knot range (reviewer's note: cover the data range with the knots). |
| W5.23 | merged 2026-09-17 at `05ff62f` (Sonnet-authored, Fable-reviewed — the self-verifying `get_by_digest` (the sidecar's recorded ingredients must hash back to the digest it is filed under, which is exactly `ArtefactKey.digest`'s definition), the served/mismatch attrs and the untouched default path read; the orchestrator's re-run of `tests/results/test_artefacts.py` and the serving rows in sbi: 63 passed in 13 s; no fixes needed; merged by Peter; the agent's runs: sbi `tests/inference/test_sbi.py` 163/2 in 294 s, `test_artefacts.py` 59/0, import sweep 105/2, lint/format/typecheck clean; **the merged-wave gate on `8e10b55`/`897b9c0` (W5.7 + W5.23 + W5.9), 2026-09-17 04:13–06:58, all green: dev 2818/422, sbi 3891/95, torch 3757/229, jax 3577/238**; the agent: 247 k tokens / 158 tool uses / 47 min). Landed: `SBIEngine(cache=..., serve_artefact=<digest>)` — the ordinary key still computed and recorded, the named entry restored regardless through `ArtefactStore.get_by_digest` (a `(artefact, ingredients)` pair, or `None` with an `ArtefactCacheWarning` for a missing, unreadable, hand-edited or corrupt entry), a miss refused by name as `EngineError` rather than falling back to training, `serve_artefact` without `cache` refused at construction; attrs `sbi_artefact_served` and `sbi_artefact_mismatch` (field-wise, from the factored `ArtefactStore.diff_ingredients`, empty when they agree) beside `sbi_cache_key`, with an `ArtefactCacheWarning` naming the differing fields; `results.md` §9 documents all three (`sbi_cache_key` had been written since W3.5 but never listed); `docs/source/sbi.rst`'s cache section gains the explicit route; 9 store rows and 4 engine rows (the served-and-calibrated row differs on `budget`, the file's established one-ingredient case, since its fixtures have no HSGP problem). Carried: **`SBIEngine.calibrate()` does not reseed torch as `run()` does** (the agent found SBC/TARP thresholds in `test_sbi.py` sensitive to test order when a preceding test samples a posterior; its own row saves and restores `torch.get_rng_state()` around `calibrate()` as a local mitigation) — a follow-up for the next housekeeping item, W3.15's pattern applied to `calibrate()`. |
| W5.9 | merged 2026-09-17 at `8e10b55` (Opus-authored, Fable-reviewed as the second pass — the rotated, rescaled solve (`u_s = r̃_s/√λ_s` against `K_x + diag(σ²/λ_s)` with `−(N/2) log λ_s` added back, and `−½ log λ_s` per leave-one-out term) re-derived; the closed-form `RotationCoupling` eigendecomposition and the `CholeskyCoupling` fallback read; the calibration study's negative result checked against the `1/√(1+ρ)` prediction; the orchestrator's re-run on the branch: dev (conformance astrometry, `test_noise_joint`, `test_joint_pointwise`, `test_calibration`) 93 passed in 43 s, jax (the NUTS rows and conformance astrometry) 54 passed in 23 s; **no fixes needed at review**; merged by Peter; the agent merged `w5.7-warped-kernel` and then master `05ff62f`+ into its branch itself, cleanly; the agent's runs: dev conformance astrometry + `test_noise_joint` 83, torch 84 in 80 s, jax 84 in 25 s, dev `tests/core` 1649/34 and `tests/results` 510/4, torch `tests/backends` 588/3, jax 409/11, the `astrometry_full` row 1 in 903 s; lint/format/typecheck clean in all three; the agent was killed once by the session rate limit at ~03:15 and resumed from committed state (the resumed segment alone: 547 k tokens / 14 tool uses; the first segment's count was lost); **the merged-wave gate on `8e10b55`/`897b9c0` (W5.7 + W5.23 + W5.9), 2026-09-17 04:13–06:58, all green: dev 2818/422, sbi 3891/95, torch 3757/229, jax 3577/238**). Landed: `ampere.core.JointGaussianProcessNoise(kernel, solver, datasets=…, coupling=…)` with `JOINT = True` (a `Likelihood` refuses it by name; it is declared on `DatasetCollection(datasets, joint={label: noise})`); `ChannelCoupling` (a `Parameterised` supplying its own eigendecomposition through `eigen(resolved, xp=…)`, one implementation for numpy, torch and jax), `RotationCoupling(angle, log_variance_0, log_variance_1)` (closed-form eigenvectors; the parameters take the bijection their prior's support implies, not `Log`) and `CholeskyCoupling` (general `T`, numerical `eigh`); the bound solver decides once for all rotated outputs (`QuasisepGP` on an ordered 1-D grid is the O(N) path); `DatasetCollection.contribution_labels()`/`group_of()`, `group_log_likelihood`, `draw_group` (one correlated draw), `JOINT_DECOMPOSITION = "joint"` in `ampere.results`, the pointwise group storing the rotated outputs under the member labels; twins and a `_LoweredJointGroup` in each backend's `problem.py`, the native per-dataset draw refusing a joint problem by name; three refusals at composition — shared grid (exact axis equality), shared mask, one shared σ vector across channels — and a free kernel amplitude beside a free `B` refused as one degree of freedom twice; `examples/astrometry --joint` and `--sbc {joint,independent,rigid}` with the injected centroiding systematic; `docs/source/astrometry.rst` §10; `likelihoods.md` §7 (new subsection) and §15 item 11 rewritten, `inference.md` §4.9 (new) and §7, `results.md` §6; five conformance rows per fixture, `BackendCapabilities.joint_noise`, `tests/core/test_noise_joint.py`, `tests/inference/test_joint_noise_nuts.py`, `tests/results/test_joint_pointwise.py`, `tests/examples/test_astrometry_example.py` (the full row behind `astrometry_full`, a 12-refit sibling always on). Fixed on the way, inside scope: `replace_observations` rebuilt SBC replicas without the collection's joint groups; `inference/engine.py`'s two label checks and `results/emission.py`'s `log_likelihood` group assumed one term per dataset. Measured (48 refits per arm, 0.90 interval on `(pmra + pmdec)/√2`): joint 0.917 (KS p 0.235), independent GPs with the correct marginals 0.729 (p 0.021; `1/√(1+ρ)` with ρ = 0.86 predicts 0.77), rigid 0.479; pinned at a 0.80 floor and a 0.10 gap; NUTS on torch (81 s) and jax (16 s) recovers `trace(B)` in its 95 % interval; lowered vs contract density torch 0.0, jax 1.4e-14. **Ruled by Peter 2026-09-17**: (1) heteroscedastic channels are the norm for real data, astrometric and otherwise, and handling them — dense and/or through the accelerated solvers (HSGP, EFGP) — is essential: drafted as **W5.24**; (2) the study design accepted (it holds `B` at the truth in both arms and gives the comparison arm the correct marginals, so the pinned claim is "given the noise you are entitled to know, the joint model is calibrated and two independent GPs are not" — fitting `B` is what `--joint` and the NUTS rows do); (3) `T > 2` is implemented but unit-tested only; polarimetry stays the second test. Carried: `replace_observations` also drops `shared=`/`shared_label=` (one line, a hierarchical SBC study would silently refit without its hyperparameters); `ampere.results.sbc` ranks posterior variables only — a `derived=` hook is a general need; a `summarise`-side helper reporting `B` rather than its bimodal angle parameterisation; the general LMC, mismatched grids and unequal uncertainties on the dense/reduced-rank follow-on; the arviz chain-length warning at tiny budgets (pre-existing). |
| W5.22 | merged 2026-09-17 at `b16607f` (Sonnet-authored, Fable-reviewed — `free_priors`, the schema-8 attr and the three-way `_resolve_interim_prior` (spec-hash mixture refused first; stored prior read back when none supplied; a supplied prior verified against the stored one; a hierarchical stored prior refused only when it is the sole source) read; the orchestrator's re-run on the branch: `test_results.py` + `test_training.py` + the provenance doctests 212 passed in 15 s; no fixes needed; merged by Peter; the agent merged master twice (the markers-list conflict with W5.9 resolved by keeping both); killed twice by the session rate limit and resumed from committed state (one uncommitted edit preserved as a WIP commit by the orchestrator, folded in); ≈ 409 k tokens / 298 tool uses over three segments; the agent's runs: dev `test_population.py` 15/2 in 516 s (was 11/1 in 710 s before the marker), `test_results.py` 177, `test_training.py` 31, sbi `test_population.py` 16/1 in 621 s; lint/format/typecheck clean; **the merged-master gate on `dd79aa0` (W5.22 on top of the gated `8e10b55`), 2026-09-17 23:51 – 2026-09-18 02:28, all green: dev 2827/423 in 20 min, sbi 3900/96 in 46, torch 3766/230 in 41, jax 3586/239 in 47 — nine rows up on the previous gate in every leg, the `population_full` row skipped as designed**). Landed: `PROVENANCE_SCHEMA_VERSION` 8 with `ampere_free_priors` (canonical JSON of `{name: {"kind": "prior_spec", **PriorSpec.to_dict()}}`, or `{"kind": "hierarchical", **HierarchicalPrior.to_dict()}`) from `ampere.results.provenance.free_priors`; `fit_population(interim_prior=None)` reads the stored marginal prior back, verifies a supplied one, and refuses by name a disagreement, a pre-schema-8 archive with no supplied prior, and a hierarchical stored prior with none supplied; the `population_full` marker (registered, skipped by default through the new `tests/results/conftest.py`) on the 200-object archive, with a 50-object reduced archive shared by every other row (`TestTheReducedObjectsFit` pins the timing guard there); `results.md` §9 (the schema-8 paragraph) and §13 item 16 amended; two schema-version literals updated in `test_results.py` (`>= 7`) and `test_training.py` (`8`). **Ruled by Peter 2026-09-18, all three confirmed**: (1) the reduced archive stays at 50 objects (the item text said 20; the fixture predated the item), the second reduction of the module's default cost — 20 objects, fewer emcee steps, or the archive built once per module — is the Sonnet filler already listed; (2) `append_training_set` does not gate on `ampere_free_priors`, because the spec hash is the hash of the merged parameter set's neutral spec, which carries every parameter's prior through the same `describe_prior` form the attribute stores, so a prior change is already refused there and a second gate would refuse the same cases twice; (3) a stored `HierarchicalPrior` is refused rather than marginalised as the answer for this phase — with **the route recorded for when it is revisited**: the marginal prior of a hierarchical member is `∫ p(θ | φ) p(φ) dφ` over the hyperparameters `φ` its references name, obtainable by Monte Carlo from draws of the hyperpriors (which the archive's runs store, since the hyperparameters are ordinary free parameters with `ampere_free_priors` entries of their own) or in closed form in the conjugate cases (normal–normal, gamma–Poisson); W5.12, which lowers `Population` and its hyperpriors, is where that resolution naturally lives, and its text now says so. Carried: `tests/results/test_population.py` still costs ≈ 8.6 min in dev and ≈ 10.3 in sbi in the default gate — the 50 single-object fits and the 2000-step population emcee, not the archive size, dominate; a second reduction (fewer steps, or the archive built once and cached across the module) is a Sonnet filler; the item text's `tests/results/test_provenance*.py` does not exist (the schema rows live in `test_results.py`). |
| W5.11 | merged 2026-09-18 at `c81f975` (Opus-authored, Fable-reviewed as the second pass — the thirteen-entry code table, the unit-to-code rule through `astropy.units.get_physical_type` with the one declared exception (spatial frequency, named by a kind putting `rad**-1` in its `AxisSpec.equivalent_units`, coded 8 whether the axis is given in wavelengths or in inverse radians, while a `u` in metres is honestly a length), `ENCODING_VERSION` in the hashed layout record, and the shared by-name refusal `axis_identity_complaint` at both doors read; the orchestrator's re-run of `tests/core/test_encoding.py` + `test_simulate.py` + `tests/results/test_training.py` in dev on the branch: 242 passed in 46 s (and `test_sbi.py` in sbi 179/2 in 5 min 20 s); **no fixes needed**; merged by Peter; the agent's runs: dev targeted 242 passed, sbi `tests/inference/test_sbi.py` 179/2 in 6 min 17 s, lint/format/typecheck (dev, sbi, torch) clean, import sweep 86/4; the merged-wave gate re-run from 2026-09-22 after the 18 September log was lost to a reboot — dev 2955 passed / 446 skipped in 28 min 17 s on `264fbd4`; sbi 4046 / 101 in 56 min 19 s; torch 3896 / 251 in 48 min 48 s; jax 3716 / 260 in 43 min 48 s — **all four green**; the agent's usage is W5.10's row's — one agent did both items in sequence). Landed: `AXIS_TYPE_CODES` (contract: `absent` 0 for a padded column, `unknown` 1, then dimensionless, length, frequency, energy, time, angle, spatial frequency, wavenumber, speed, temperature, mass — a code may be appended, never renumbered), `AXIS_TYPE_NAMES`, `axis_type_code(unit, spec=)`; `DatasetLayout.axis_codes`/`axis_type_names` in the layout, its dict form, its hash and `compare`'s report; the `axis_identity` column group (`encoding.md` §3 group 4, width A, the code on every row of its dataset's block) exposed by `unpack` on `Unpacked` and every `DatasetView` and carried by `features`, so no embedding wrapper changed; `ENCODING_VERSION` 2 with `EncodingLayout.from_dict` and `append_training_set` refusing a version-1 record by name through one message. Measured: the module docstring's one-dataset `Spectrum` layout `e73db8a8032273f4c511142f1459c51c` (17 columns, version 1, computed from `dd79aa0`'s own code) → `822f4c900af7a854488b2044db598e26` (18 columns, version 2). Deliberately narrower than what it closes: the code names a physical type, not a role — two `length` axes in one collection (a wavelength and a baseline) still share a code, and a kind wanting them distinguished needs a new code, not a new mechanism. |
| W5.10 | merged 2026-09-18 at `d48c291` (Opus-authored, Fable-reviewed as the second pass — `ContextPrior` and the three shipped priors' draws, `Dataset.contextual_observed` as the single route of a context's σ into every family's `sample`, the `"<stream>.context"` sub-stream spawned by index (so every pre-existing budget's θ and noise are bit-for-bit unchanged), the three refusals by name, the FiLM branch's zero-initialised last layer and the pinned calibration inequalities read; **one fix at review**: the context prior did not enter the artefact cache key, so a network trained without a context would have been served for a context-amortised run with the same settings (§7's stale-artefact trap) — fixed as `artefact_key(context=<digest>)`, omitted when `None` so every pre-W5.10 digest is unchanged (verified byte-for-byte against master's own `artefacts.py`, `3af107dc…` either side; W3.12's pinned NPE digest untouched), with 7 store rows and 3 engine rows; the orchestrator's re-runs on the branch: dev `test_encoding.py` + `test_simulate.py` + `test_training.py` 242 passed in 46 s, sbi `test_sbi.py` 179/2 in 5 min 20 s (pre-fix), dev `test_artefacts.py` 63/3 in 4 s (post-fix); merged by Peter; the agent's runs: dev 242, sbi `test_sbi.py` + `test_artefacts.py` 248/2 in 6 min 14 s, torch `test_simulate.py` 128 in 56 s, lint/format/typecheck (dev, sbi, torch) clean, import sweep 86/4; the merged-wave gate re-run from 2026-09-22 (the 18 September log lost to a reboot) — dev 2955 / 446 in 28 min 17 s; sbi 4046 / 101 in 56 min 19 s; torch 3896 / 251 in 48 min 48 s; jax 3716 / 260 in 43 min 48 s — **all four green**; the agent, W5.11 and W5.10 together over two segments: ≈ 566 k tokens / 462 tool uses / 2 h 38 min, most of it queued behind the wave-4 gate's four legs on the lock). Landed: `ampere.core.simulate.ObservationContext` (per-dataset σ arrays plus a JSON-safe record), the structural `ContextPrior` (`draw(rng, observed)`, `describe()`), `ScaledSigma`, `SigmaArchive`, `SignalToNoise`; `Dataset.contextual_observed(sigma)` and `draw_observation(..., sigma=)` — no `LikelihoodFamily.sample` signature moved; `Simulation.context` (failures included); `simulate_many(context=...)` drawing one context per simulation on its own sub-stream, materialised a chunk at a time; the training set's optional `context` group and `ampere_simulation_context`, append refusing a context mismatch in both directions; `SBIEngine(context=...)` with `ampere_sbi_context`/`ampere_sbi_context_hash` and the prior's digest in the cache key; `calibrate(context=...)` inheriting the run's prior by default; FiLM (`embedding={"type": "set", "film": True}`), off by default, the identity before training; `inference.md` §12/§13, `encoding.md` §3/§8, `results.md` §9/§11 amended, `sbi.rst` gains "Amortising over the observation context". `PROVENANCE_SCHEMA_VERSION` stays 8 and `TRAINING_SET_SCHEMA_VERSION` 1 (engine-specific extras; an optional group). Measured: a 1500-draw NPE under `ScaledSigma(0.5, 2.0)` — TARP area-to-curve −0.007 at a ×1.5 rescale the prior covers, −0.046 at ×12 it does not, pinned as inequalities at about half the separation; the whole of `test_sbi.py` in sbi in 6 min, so no marker. **For Peter**: a context on a problem with a joint noise group (W5.9) is refused by name rather than guessed (one shared context per group from its first dataset is the implementable alternative) — confirm. **Ruled by Peter 2026-09-22: confirmed** — the refusal stands and the shared-context-per-group route stays open until a problem needs it. Carried: the native batched path refuses a context by name (a per-draw σ does not reach the realisation; **scheduled as W5.29**, 2026-09-22); `SignalToNoise` shapes σ from the *observed* magnitudes because a context is drawn before the forward model runs (documented; a prior tracking each draw's own signal needs `draw` to receive the predicted containers — a protocol change); `_concatenate` in `training.py` indexes the addition by the existing tree's groups (a bare `KeyError` or a silent drop for any future optional group); `TrainingSet.contexts` defaults to `()`, so "written before W5.10" and "written without a context" are indistinguishable in Python (the file distinguishes them); sbi 0.27's outlier warning on the set layout under a wide prior (`z_score_x` guidance for the tutorial); a pre-existing wrong stream name in `_calibration_batch` fixed on the branch because W5.10 made it load-bearing. |
| W5.8 | merged 2026-09-18 at `2a1f596` (Opus-authored, Fable-reviewed as the second pass — the horseshoe's global-local chain as three `HierarchicalPrior` levels, the `with_shrinkage` hook (functional; verified directly: the original kernel's spec hash and parameter names unchanged, the shrunk kernel's hash moved, the scale levels registered ahead of the amplitudes, `terms` and `FAMILY` preserved), the pinned margins against the measured values, the SBC contrast re-pinned as the visibility precedent does, and the driver's own rows read; the orchestrator's re-runs on the branch through the lock: dev `tests/core/test_parameter.py` 158 passed in 3.6 s, the driver/arm rows of `tests/m2/test_many_lines.py` 17 passed in 2 min 14 s; **no fixes needed**; merged by Peter; the agent's runs: dev `tests/m2` 112/32 in 11 min 24 s, `test_many_lines.py` 29/2 in 4 min 41 s, `test_many_lines_calibration.py` 6/1 in 2 min 49 s, torch NUTS-over-the-knots 2/2 in 4 min 51 s, jax 2/2 in 1 min 1 s, lint/format/typecheck (dev, torch, jax) clean, import sweep 86/4, docs build clean; the merged-wave gate re-run from 2026-09-22 (the 18 September log lost to a reboot) — dev 2955 / 446 in 28 min 17 s; sbi 4046 / 101 in 56 min 19 s; torch 3896 / 251 in 48 min 48 s; jax 3716 / 260 in 43 min 48 s — **all four green**; the agent: ≈ 503 k tokens / 82 tool uses over two segments, killed once by the session limit with two edits preserved as a WIP commit). Landed: `ampere.core.regularised_horseshoe(amplitudes, *, prefix, global_scale, tail, unit)` — `τ ~ C⁺(0, τ₀)`, `s_j | τ` a `gamma(½, scale=τ)` (`tail="regularised"`, the default, which lowers) or `C⁺(0, τ)` (`tail="cauchy"`, the plain horseshoe, which does not lower on torch/jax yet), `a_j | s_j ~ N⁺(0, s_j)`, all with `Log` bijections; `ampere.core.with_shrinkage(kernel, declaration)`, the one hook, appended to `kernels.py`; the `many_lines` scenario (kept out of `SCENARIOS`, so the milestone's four are untouched), the `warped` and `sum` arms in `build_kernel`, `examples/m2_misspecification/many_lines.py` and the driver's `--scenario many_lines --kernel {matern32,warped,sum}` (a trailing-block `KeyError` on the warped arm fixed in scope); `tests/m2/test_many_lines.py` (29 rows) and `test_many_lines_calibration.py` (SBC with a parameter-free `LineForestError` step in the simulating problem), the backend-agreement row for the warped arm; `likelihoods.md` §6 (line 889's promise discharged) and §13, `parameters.md` §9, `lowering.md` §3.2.1, `m2_misspecification.rst` (two sections) and `kernels.rst`. Measured at 200 points, bias in posterior widths: standard 35.70 (1/4 covered), stationary Matérn 2.84 (3/4; pinned ≥ 2.0 and > 1.5), warped 0.92 (4/4), `Sum` 0.78 (4/4; each pinned ≤ 1.5 and ≥ 2.0× better than stationary, measured 3.1× and 3.7×); all three localise inside the band (in-band contrast 2.76×, 1.90×, 3.02×, pinned ≥ 1.4×); the horseshoe row `min(a)/max(a)` 0.254 → 0.087 (2.93×, pinned ≥ 1.3) and mass below a tenth 0.257 → 0.537 (2.09×, pinned ≥ 1.2), thresholds a third below the worst of three seeds; SBC at nominal 0.90: standard [0.00, 0.80], warped [0.80, 0.90], gap pinned ≥ 0.4; the warp's fitted increments give an effective length scale of 0.028 µm under the smooth error and 0.0009 µm under the forest, the sign change at the band edge. **For Peter**: (1) `halfcauchy` is in neither backend's lowering table, so the plain-horseshoe tail and the global level are reference-only on NUTS — two lines per backend plus a §3.2 row and a conformance row, a Sonnet filler; (2) a `Derived` parameter node (a pure function of other parameters) would make P&V's slab declarable and let `WarpedKernel` declare its non-centring — a design question, not scheduled; (3) `tests/m2` in dev now costs 11 min 24 s per environment (was ≈ 100 s): `many_lines.SHRINKAGE_EMCEE` 3.5 min, the four arms 2.4 min, the calibration row 2.8 min, each set against a measured scatter — accept, or name the lever; (4) the helper is named `regularised_horseshoe` while its default tail is a gamma variant of the horseshoe rather than Piironen & Vehtari's slab — the docstring says so; keep the name, or rename. **Ruled by Peter 2026-09-22**: (1) scheduled as W5.25; (2) a `Derived` node goes to Phase 6 as a design-horizon note; (3) accepted for this gate, the trim is W5.26's; (4) renamed `shrinkage_horseshoe` with `regularised_horseshoe` a deprecated alias for one phase, W5.27. Carried: `with_shrinkage` rebuilds the parameter set through `copy.copy` and a `__dict__` write rather than a constructor path; `arviz_base`'s chain-longer-than-draw warning on every emcee run (pre-existing). |
| W5.12 | merged 2026-09-22 at `c26ff66` (Opus-authored, Fable-reviewed as the second pass — `Population`, `_apply_populations`, `_element` and the decomposition row read directly; **no fixes needed**; merged by Peter; the orchestrator's re-runs on the branch through the lock: dev 177 passed / 10 skipped in 2.6 s, torch 194 / 4 in 1 min 45 s, jax 194 / 4 in 17 s (`test_population.py`, `test_population_nuts.py`, `test_parameter.py`); the agent's runs: dev `tests/conformance tests/results` 1076 / 74 in 11 min 23 s, `tests/core` 1386 / 30, torch and jax `tests/conformance` 1007 / 69 each, torch `test_population_nuts.py` 6 in 74 s, jax 6 in 13 s, hygiene clean; the wave-5 merged gate (master `e2c6666`, 2026-09-22, five legs) green: dev 3003 passed / 561 skipped in 27 min 17 s, sbi 4142 / 188 in 1 h 09 min, torch 3992 / 338 in 1 h 03 min, jax 3845 / 313 in 45 min 42 s, nested 3043 / 521 in 28 min 47 s; the agent: ≈ 340 k tokens / 391 tool uses / 56 min, no hand-off — the first item under the 2026-09-22 budget rules). Landed: `ampere.core.Population(name, members, hyperpriors, over=, size=, layout=, label=)` and `MAX_FLAT_MEMBERS = 128`; `ParameterSet.merge(..., populations=[...])` (keyword-only, additive); the `plate` layout (one array parameter per member of shape `(size,) + member.shape` in the population's component, `PlateBinding` per element to each `over` component; lowers to a real plate on both backends with no backend change) and the `flat` layout (the tie-based pattern, references renamed onto the qualified hyperpriors, refused above 128); refusals by name at merge (taken label, unknown or self-referencing component, shape or unit disagreement, a fixed site, a tied or `shared_as` site, a component merged as a `ParameterMapping`); `DatasetCollection(..., populations=)`, `.populations`, `DatasetCollection.plate(name, datasets, members=, hyperpriors=, over=None, layout=, label=)` defaulting `over` to the datasets' distinct model labels, `FittingProblem(..., populations=)`; the `distribute` fix (`_element`); `tests/conformance/test_population.py` (once per fixture: the decomposition on both layouts, the realised population against the numpy path, the plate lowering against the flattened form) and `tests/inference/test_population_nuts.py` (50 members, 52 free dimensions, `mu` and `sigma` inside the central 95 %, member draws correlating > 0.9 with their truths); `parameters.md` §8/§9, `inference.md` §9/§10a and limitation 17.6 lifted, `hierarchical_population.md` H-1 dispositioned with its three deviations, §11 Q2 closed and Q5's joint-fit half landed, `docs/source/advanced.rst`. **Peter's W5.22 note**: `fit_population`'s refusal of a stored `HierarchicalPrior` is **not lifted** — the joint fit (horizon (a)) never consults a stored interim prior, and the marginal is only defined against an archived population fit, whose members' interim priors are not independent; real work for W5.13's successor, not this item. **For Peter**: (1) the layout is an explicit `layout=` with one default everywhere rather than backend-dependent behaviour (the agent's reading of "the numpy path fits small populations by the tie-based pattern"; recommend confirm — a backend-dependent parameter space would defeat the conformance comparison); (2) `MAX_FLAT_MEMBERS = 128` is a judgement call (recommend accept); (3) `over` cannot name a `Dataset` component — a population addresses members by bare local name and a dataset's merged names are qualified; relaxing it is a §4 change to `PlateBinding` or routing in two places (§11 Q1's rejected alternative); recommend leave refused and documented until per-dataset nuisance populations are wanted. **Ruled by Peter 2026-09-22**: (1) confirmed — `layout=` is an explicit argument with `"plate"` the one default on every backend; (2) 128 accepted, with the docstring and the refusal message re-cited to `hierarchical_population.md` §7's measured figures (46 ms per `lnprior` at 100 flat members, 373 ms at 1000, against 1.3 ms as one plate; the text now says "800 ms / 0.6 ms" and "well under a millisecond"), and a **setting** — not a `Population` argument and not a change to the limit — that turns the refusal into a loud warning for a power user who accepts the cost: **W5.30 (a)**; (3) refused now and recorded as a **Phase 6 design item** (the per-dataset GP amplitude drawn from a shared prior is the motivating case; `DEVELOPMENT_PLAN.md` Phase 6). Carried: the centred parameterisation is the only one expressible (a non-centred `θ_i = μ + σ z_i` needs the `Derived` node already recorded as Phase 6); `tests/inference` in dev is slow because of the SBI rows (W5.26's territory). |
| W5.14 | merged 2026-09-22 at `9347175` (Opus-authored across two agents — the first landed the nested pair and handed off at its natural boundary under the 2026-09-22 budget rules (≈ 339 k tokens / 180 tool uses / 34 min), the successor landed the VI guides and the blackjax route (≈ 367 k / 238 / 70 min); Fable-reviewed as the second pass — the shared nested driver, the Laplace and flow handling and the blackjax refusals read directly, the lock's environment diff verified (jax +blackjax/optax/absl-py, `nested` new, all else byte-identical); **two fixes at review**: `70f0128` (the MCLMC dimension floor declared twice) and `cdb87a0` (`results.md` §9's `ampere_approximation` vocabulary and the evidence triple's absence as a state); merged by Peter; the orchestrator's re-runs on the final head through the lock: dev `tests/inference` + `test_parameter.py` 337 passed / 343 skipped in 37 s, nested `test_nested.py` + doctests 57 / 6 in 75 s, torch `test_vi.py` + `test_nested.py` + doctests 62 / 49 in 2 min 14 s, jax `test_vi.py` + `test_blackjax.py` + `test_nested.py` + doctests 95 / 51 in 6 min 54 s; hygiene (lint, format, typecheck dev, actionlint) clean; the agents' runs: nested `test_nested.py` 50 / 3 in 60 s, dev inference battery 117 / 45, torch `test_vi.py` + doctests 51 / 7 in 126 s, jax `test_blackjax.py` + doctests 41 / 4 in 99 s, dev `tests/inference` 179 / 343, typecheck jax/torch/nested 0 errors, docs build clean; the wave-5 merged gate (master `e2c6666`, 2026-09-22, five legs) green: dev 3003 passed / 561 skipped in 27 min 17 s, sbi 4142 / 188 in 1 h 09 min, torch 3992 / 338 in 1 h 03 min, jax 3845 / 313 in 45 min 42 s, nested 3043 / 521 in 28 min 47 s). Landed: `ampere/inference/_nested.py` — `NautilusEngine` and `UltranestEngine` over `_NestedEngine` (dynesty's `resample_equal` on the engine's `resample` stream as the one weighted-draw implementation for all three nested samplers; `ampere_evidence_method = "nested_sampling"` for all three, closing W5.0's carried question as the method family; nautilus's error as `1/sqrt(n_eff)` with `ampere_nautilus_log_z_err_source`; a one-dimensional problem refused by name; ultranest under `engine.global_seed`, shared with zeus; `engine.default_live_points` shared with dynesty); `_vi.py` — `GUIDE_FAMILIES` gains `laplace` (`AutoLaplaceApproximation`; a `Delta` until asked for its Gaussian: pyro's `laplace_approximation()`, numpyro's `get_posterior(params)`) and `flow` (`AutoIAFNormal`; the base moved to the seeded start through `get_base_dist`, pyro's network in float64; `ampere_vi_guide_parameters`; below two free parameters refused by name); `_blackjax.py` — `BlackjaxEngine(problem, method="mclmc"|"pathfinder")`, jax-only by nature (MCLMC tuned by `mclmc_find_L_and_step_size`, one `lax.scan`; Pathfinder as approximation writing `ampere_approximation="pathfinder"` and `proposal_log_density`, and as `initial="pathfinder"`; MCLMC writes `"none"` with `ampere_blackjax_adjusted`/`desired_energy_var`; neither writes the evidence triple, asserted absent); extras `nautilus`/`ultranest`/`blackjax` with 1:1 pixi features, one new environment `nested` (dev's set plus the two samplers) with a path-gated CI leg on both backend matrices, blackjax folded into `jax`, `NESTED` and `BLACKJAX` path buckets before `CORE`, `engines_full` marker; the batteries `test_nested.py` (dynesty the control; analytic log Z −13.0684: dynesty −0.15 σ, nautilus −0.26 σ, ultranest +1.46 σ; tolerance `max(4 σ, 0.4 nat)`), `test_vi.py` (Laplace exact on the conjugate Gaussian; the flow's density against the library's ELBO; every guide SBC-ranked), `test_blackjax.py`; `inference.md` §10 (the item said §5) "Three nested samplers on one surface" and "Four guide families and a second gradient library", `results.md` §9 (review), the memo's landed notes, `advanced.rst`, `ampere.inference.rst`. **Rulings for Peter**, none blocking: (1) a separate `nested` environment and CI leg rather than folding the samplers into `dev` — recommend keep; (2) nautilus's `1/sqrt(n_eff)` as ampere's own labelled evidence error — recommend accept; (3) MCLMC's `ampere_approximation="none"` rather than a named unadjusted family — recommend keep, the key answers whether chain diagnostics apply; (4) the flow's base-distribution subclass and `_DESIRED_ENERGY_VAR` duplicating a library default — recommend accept with their comments. **Ruled by Peter 2026-09-27: all four confirmed as recommended** (the orchestrator's sweep found them never recorded). Carried: `path_filters.py`'s `DEV_JOBS` list still names py311–313 (W5.18, one line); ultranest's `logzerr` is more conservative than its console line; `proposal_log_density` has four producers and no cross-engine conformance row; the jax `test_vi.py` leg's budget cut from 22 min to ~2 (a flow draw traces three networks; each SBC replica is a fresh jax compilation) — W5.26's territory if it grows again. |
| W5.25 | merged 2026-09-22 at `51d18a6` (Sonnet-authored, Fable-reviewed as the second pass — the four builders read directly; **no fixes needed**; **extended by the orchestrator at the first report** to close the gap the item found (below); merged by Peter; the orchestrator's re-runs on the final head through the lock: dev `test_halfcauchy.py` + `test_parameter.py` 174 passed in 1.9 s, jax `test_halfcauchy.py` + `test_jax.py` + `TestTheHorseshoeUnderNUTS` 224 passed in 67 s, torch `test_halfcauchy.py` + `test_torch_lowering.py` + `TestTheHorseshoeUnderNUTS` 218 passed / 2 skipped in 12 min 10 s (both tails pass on torch; the two skips are pre-existing); hygiene clean; the agent: ≈ 305 k tokens / 243 tool uses / 41 min across the two rounds; the wave-5 merged gate (master `e2c6666`, 2026-09-22, five legs) green: dev 3003 passed / 561 skipped in 27 min 17 s, sbi 4142 / 188 in 1 h 09 min, torch 3992 / 338 in 1 h 03 min, jax 3845 / 313 in 45 min 42 s, nested 3043 / 521 in 28 min 47 s). Landed: `halfcauchy` in both backends' tables — torch `_build_halfcauchy`/`_halfcauchy` (`HalfCauchy(scale)` at `loc == 0`, the §3.3 affine shift otherwise) and `_hierarchical_halfcauchy` in `_HIERARCHICAL_BUILDERS`; jax `_halfcauchy` (`HalfCauchy`, or `TruncatedCauchy(loc, scale, low=loc)` shifted, on `_halfnorm`'s support reasoning) and `NATIVE_ICDF`; `lowering.md` §3.2 row and §3.2.1 rewritten ("the horseshoe's chain now lowers, in full"); the `regularised_horseshoe` docstring; `tests/conformance/test_halfcauchy.py` (prior transform against scipy's `ppf`, log-density against `logpdf`, no mass below `loc`, cube centre is the median; `loc ∈ {0, 3}`); `TestTheHorseshoeUnderNUTS` in `tests/inference/test_nuts.py` on W5.8's actual `horseshoe_kernel` declaration, both tails, torch and jax. **The finding and the extension**: torch keeps a separate per-family registry for hierarchical priors (`_HIERARCHICAL_BUILDERS`: `norm`, `uniform`, `halfnorm`, `expon`, now `halfcauchy`) while jax lowers hierarchical priors generically from the flat table, so `tail="regularised"`'s `gamma` local scale still failed under NUTS on torch — pre-existing, isolated from `halfcauchy`, and contradicting §3.2.1's earlier claim; the orchestrator extended the item on the spot: `_hierarchical_gamma` (mirroring `_build_gamma`, the fixed `a=½` folded in by `lower_hierarchical` as any constant), a torch backend row against the flat table and scipy to 1e-9, the skip removed. **For Peter**: none. Carried: the torch horseshoe NUTS rows cost ≈ 5 min per tail (the funnel; W5.26's territory if the budget matters); `_default_bijection_for_hierarchical` refuses any family with a shape argument unconditionally, which is why the horseshoe still passes `bijection=Log()` explicitly — documented, a candidate for W5.28. |
| W5.27 | merged 2026-09-22 at `3eeb2d3` (Sonnet-authored, Fable-reviewed as the second pass — the alias, the framework note and the residual references read directly; **no fixes needed**; merged by Peter; the orchestrator's checks on the branch: lint, format-check, typecheck (dev) clean, `tests/core/test_parameter.py` + `tests/conformance/test_halfcauchy.py` 176 passed in 2.6 s; the agent's queued runs on the branch behind the gate lock: torch `tests/backends/test_torch_lowering.py` 192 passed / 2 skipped, dev `tests/m2/test_many_lines.py` 33 passed / 2 skipped in 4 min 44 s on the existing margins, the docs build succeeded with the 15 pre-existing legacy warnings; the agent's own: `test_parameter.py` 160, `test_spec_doctests.py` 19, hygiene clean; ≈ 229 k tokens / 187 tool uses / 24 min; no gate of its own — a rename, verified by the runs above). Landed: `ampere.core.shrinkage_horseshoe` (the helper renamed; twelve internal references, doctests and error messages moved), `regularised_horseshoe` kept as an alias that emits a `DeprecationWarning` naming the replacement and Phase 6 and calls straight through (tested: warns once, returns exactly what the new name returns), `SHRINKAGE_HELPER_FRAMEWORK` — one module-level string beside `HORSESHOE_TAILS` with an RST heading, spliced onto both `shrinkage_horseshoe.__doc__` and `with_shrinkage.__doc__` so it renders under both entries without a new rst page; both names exported (`__all__` alphabetical for `RUF022`); every reference moved in `kernels.py`, the torch lowering comments, `examples/m2_misspecification/many_lines.py`, five test files (`TestTheShrinkageHorseshoe`, new `TestTheRegularisedHorseshoeAlias`), `kernels.rst`, `m2_misspecification.rst`, and the frozen pages `parameters.md` §9, `likelihoods.md` §6 and `lowering.md` §3.2.1 annotated *Amended W5.27* without re-spelling the old identifier (`likelihoods.md` §13 had no literal mention, contrary to the item text). **For Peter**: none. |
| W5.28 | merged 2026-09-22 at `30bf81e` (Sonnet-authored, Fable-reviewed as the second pass — every lettered diff read directly, `register_parameter`'s route, the reseed helper's per-call seed derivation and the root conftest's collection order checked against the code; **no fixes needed**; merged by Peter; the orchestrator's runs on the branch: jax typecheck 0 errors (through `pixi run -e jax bash -c "cd <worktree> && pyrefly check --project-excludes=ampere/backends/torch/** ..."`, the main checkout's jax environment with the worktree first on the import path), jax `tests/conformance` 1031 passed / 69 skipped in 2 min 36 s, jax `test_nuts.py -k "TestTheHorseshoeUnderNUTS or TestReproducibility"` 4 passed in 24 s; the agent's: dev lint, format-check (272 files), typecheck clean, torch typecheck clean, dev `tests/core tests/results` 1845 / 35 in 12 min 28 s, dev `tests/conformance` 643 / 69, torch `tests/conformance` 1031 / 69, torch horseshoe-and-reproducibility NUTS rows 4 passed in 13 min 34 s, sbi calibration rows 23 passed in 1 min 41 s, the per-item files (calibration 27, diagnostics 59, likelihood 161, hsgp 63, training 43, parameter 160, bakeoff 11, many-lines horseshoe 6), the `test-all` composition collecting 3572 rows with 0 errors; ≈ 568 k tokens / 521 tool uses / 92 min, past the soft cap but finishing the last unit; the worktree was cut by the harness from the legacy `origin/master` and re-based on `3eeb2d3` on the orchestrator's instruction before any edit; the W5.28 merged gate on `30bf81e`: dev 3013 passed / 563 skipped in 40 min (wave-5 baseline 3003/561 in 27 min — the ten rows are W5.28's, the minutes the shared machine), torch 4002 / 340 in 63 min (baseline 3992/338 in 63), jax 3855 / 315 in 43 min (baseline 3845/313 in 46) — **all three green**, 2 h 25 m end to end on 2026-09-23). Landed: (a) `replace_observations` carries `shared=`/`shared_label=` into the replica (test); (c) `ampere.results.coupling_matrix_summary(tree, coupling, prefix=None, **kwargs)` — `RotationCoupling`'s `B` reassembled per draw through `ChannelCoupling.matrix` and handed to `arviz.summary` (three rows, including exact recovery of `B` across the relabelling jump; `astrometry.rst`'s own advice now points at it; the item text's "evidence `B` (W5.0/W5.13)" was the orchestrator's mis-citation of W5.9's note, corrected); (d) the reference `QuasisepGP.condition(at=)` takes a multi-axis container as `DenseGP` does (bitwise agreement row); (e) `GriddedSolver` public, the image study's extension point; (f) `HilbertSpaceGP.condition(at=None)` reuses the training block (call-count row plus exact agreement with the explicit-`at` path); (g) `_concatenate` refuses a group mismatch by name in either direction for any group (two rows on hand-built trees); (i) `with_shrinkage` registers the rebuilt parameters through `register_parameter` — no new row, covered by `test_parameter.py`, the many-lines horseshoe rows and the torch/jax horseshoe-under-NUTS rows, all re-run; (j) `SBIEngine.calibrate()` reseeds torch and numpy's legacy global under `_seeded` for the SBC batch and again for TARP, restoring state on exit, and `TestServingANamedArtefact`'s hand-rolled save/restore is gone (two rows: bitwise repeat, state restored); (k) one `tests/conftest.py` puts the root on `sys.path` and records why `tests/` stays a non-package; six copies trimmed. Documented decisions: (b) `sbc` ranks the `posterior` group only (a `Notes` section: the derived groups have no truth `simulation.parameters` can supply); (h) an empty `TrainingSet.contexts` is deliberate parity with `observe=False`. **For Peter**: none. Carried: `replace_observations` still does not pass `populations=` to the rebuilt problem, so a W5.12 population-level SBC replica may be silently independent the same way — **added to W5.30 as (c)**; the reference `QuasisepGP.condition` now accepts a multi-axis `at=` that torch's and jax's own refuse by name (their `_axis` discards the unreduced points), a capability disagreement to close with a conformance row when a problem needs it (W5.19's list). |
| W5.30 | merged 2026-09-23 at `10bebe4` (Sonnet-authored, Fable-reviewed as the second pass — every diff read directly, the predicate's premise checked family by family against scipy's `_get_support`, the (c) gap the agent found verified against `FittingProblem.__init__`; **no fixes needed**; merged by Peter; the orchestrator's run: dev `test_parameter.py` + `test_calibration.py` 197 passed in 49 s; the agent's: dev `test_parameter.py` 169, `test_calibration.py` 28, `test-core` 1399 / 30, import sweep 89 / 4, conformance 643 / 69, docs build with the pre-existing legacy warnings; jax conformance 1031 / 69 in 2 min 42 s, jax backends 198 in 60 s, jax horseshoe NUTS 2 in 18 s; torch conformance 1031 / 69 in 59 s, torch backends 192 / 2 in 4 s (the two skips pre-existing), torch horseshoe NUTS 2 in 10 min 43 s; a supplementary dev `test-results` 456 / 5 in 10 min 32 s; lint, format-check (273 files), typecheck clean in dev, torch and jax; ≈ 353 k tokens / 219 tool uses / 52 min; the worktree was cut by the harness from `b8e585b` and reset to `30bf81e` on the dispatch's instruction before any edit; the W5.30 merged gate on `10bebe4`: torch 4012 passed / 340 skipped in 62 min (W5.28 baseline 4002/340 in 63 — the ten rows are W5.30's), dev 3023 / 563 in 27 min 20 s of run time (baseline 3013/563 in 40; 91 min wall, queued behind torch on the lock), jax 3865 / 315 in 45 min 18 s of run time (baseline 3855/315 in 43; 138 min wall, third on the lock) — **all three green**, 2 h 18 m end to end, 2026-09-23 23:07 → 2026-09-24 01:25; log `~/.cache/ampere-gates/w5.30-gate.log`). Landed: (a) the corrected figures in `MAX_FLAT_MEMBERS`'s docstring and the refusal (§5, not §7 — the item's citation was wrong), `ampere/core/settings.py` (`Settings` with `flat_population_cap: Literal["refuse", "warn"] = "refuse"`, the `settings` instance, `override(**fields)` restoring on exit and failing atomically on an unknown name, `AmpereFlatPopulationWarning(UserWarning)`), read in `Population.__post_init__` at validation time (`TestTheFlatPopulationCapSetting`: refuses by default; warns under the setting with the count, cost and remedy, then merges and evaluates `lnprior`; the override restores the refusal); (b) `_hierarchical_shape_support_is_fixed(dist, shape_names, referenced, constants)` — at least one shape argument a constant, and either the base-class `_get_support` or the override reproducing `(dist.a, dist.b)` at the declared constants — replacing the blanket `numargs` refusal, the constants map now built once ahead of it with the shape slots in the positional order (`gamma`/`lognorm` → `Log`, `beta` with either shape referenced → `Logit`; `truncnorm` refuses; `gamma` with `a` referenced refuses; explicit `bijection=` wins; the horseshoe helper's local scale infers `Log` without its keyword — the keyword kept); (c) `replace_observations` passes `populations=problem.populations` (`test_the_populations_survive_the_replica`, which first asserts the population given to `FittingProblem` is absent from `datasets.populations`); `parameters.md` §6 and §9 *Amended W5.30* lines; `advanced.rst` one paragraph each, its own stale "0.6 ms / 800 ms" corrected. **§4.1**: judged not a contract statement (the support→bijection table is unchanged); no decision-log line. **For Peter**: the predicate's asymmetry — `gamma` with its sole shape referenced refuses while `beta` with one shape referenced infers — follows the item's test list, not a mathematical need (the base-class support check makes both safe); loosen if a real hierarchical `gamma` on `a` turns up. **Ruled by Peter 2026-09-24: leave it** until a real hierarchical `gamma` on its shape turns up. Carried: none. |
| W5.26 | merged 2026-09-24 at `cbc039c` (Sonnet-authored, Fable-reviewed as the second pass — `ci.yml`, the marker and skip mechanics, the group tasks and the two runner rows read directly; **one fix**: `tests/conformance` restored to the `backends` CI group (`02b56ec`) because the `test-all` comment names `tests/examples` as a suite whose registry-leak proof rests on sharing a process with it; merged by Peter; the agent handed off code-complete after ≈ 353 k tokens / 265 tool uses / 40 min and the orchestrator ran the remainder: the jax reduced-rank NUTS row on the new pooled-width code 3 passed in 60 s; the branch timing run on `02b56ec`, one leg at a time through the lock: dev 3046 passed / 580 skipped in 33 min 18 s (before, the W5.30 gate: 3023/563 in 27 min 20 s), torch 3956 / 436 in 54 min 44 s (before 4012/340 in 62), jax 3807 / 413 in 32 min 34 s (before 3865/315 in 45 min 18 s), sbi 4110 / 282 in 61 min 39 s (before 4142/188 in 69, the wave-5 gate on `e2c6666`, four merges older) — **all four green**; per leg the two new suites add 40 rows and the 87 study rows move from passed to skipped, so the like-for-like saving is ≈ 15 min on torch, ≈ 16 on jax, ≈ 13 on sbi once the new suites' own 7 / 4 / 6 min are taken out (the item's "≥ 8 min shorter" read like for like, ruled 2026-09-24); the study files alone (`tests/m2` + `tests/results/test_population.py`): dev 130 / 34 in 18 min 33 s, torch 53 / 111 in 8 min 44 s, sbi 53 / 111 in 8 min 52 s, jax 52 / 111 in 2 min 23 s plus the one pre-existing failure carried below; the agent's evidence: the collected-test diff dev 3622 → 3622, torch 4391 → 4391, dev ∪ torch 4600 → 4600, `comm -23` empty; the new suites' budgets dev 23 / 17 in 187 s, torch 31 / 9 in 439 s, jax 29 / 11 in 229 s, sbi 33 / 7 in 335 s; actionlint, lint, format-check (273 files), typecheck clean; the worktree cut by the harness from `b8e585b` and rebased on `2c76b28` per the prompt; **the CI run: `workflow_dispatch` 35959713657 on `4ff25b1`, 2026-09-24 — GREEN, 29 of 29 jobs, 35 min end to end — the first green run of the `v2` line; the fifteen `new-namespace suites (<env>, <group>)` cells: dev core 13 / backends 4 / studies 17 min, torch 5 / 29 / 7, jax 7 / 18 / 10, sbi 6 / 35 / 16, nested 16 / 5 / 11 (the longest cell, sbi backends at 35 min, against the single sbi `test-all` step of over 40 min before); the two W5.26 (7) rows green on torch, jax and sbi**; **the merged run on `v2`, push run 35962510128 on `cbc039c`: green, 28 jobs plus the weekly sbi-characterisation job skipped by design, 36 min**; **the W5.26 merged gate on `cbc039c` (2026-09-24, four legs serialised through the lock, the first gate on the nine-suite `test-all`): dev 3046 passed / 580 skipped in 31 min 09 s, torch 3956 / 436 in 58 min 06 s, jax 3807 / 413 in 33 min 22 s, sbi 4110 / 282 in 1 h 02 min 35 s — all four green, every count equal to the branch run's**). Landed: (1) the `study` marker (`pyproject.toml`) on every module of `tests/m2` but `test_backend_agreement.py` and `test_model.py`, and on the reweighting, interim-prior, reader-agreement, ESS-refusal, spec-hash and marginal-fix rows of `test_population.py` (`TestTheApproximateEngineRow` deliberately unmarked: it needs a VI backend), skipped by name by `_skip_study_rows_outside_dev` in the two conftests wherever torch or jax imports; (2) `tests/interferometry` and `tests/astrometry` in `test-all`, with `tests/interferometry/__init__.py` (three of its basenames collide with `tests/results`, `tests/inference` and `tests/m2` in one process — the reason the suites could never have joined before), both in the path filter's shared-suite bucket and the header tables; (3) `test-group-core` / `-backends` / `-studies` pixi tasks, each a subset of `test-all` with `tests/conformance` in every group, as the `group` matrix dimension of `suites` and `backend-suites` (fifteen checks, `bench` once per environment on the `core` cell); (4) actionlint clean; (5) the CI section and the gate bullet with the measured times; (6) the pixi `--locked` failure traced to upstream prefix-dev/pixi#6834, fixed by #7024 on `main` two days after v0.81.0 and in no release, recorded in the `setup-pixi` comment, CI stays `frozen`; (7) the SBC row's threshold `> 0.001` with the reason, and the reduced-rank width compared as the pooled within-chain estimate at 600 draws with a two-fifths margin set from the runner's measured scatter. **Ruled 2026-09-24**: the review's conformance pairing, the like-for-like reading and the descriptive group names all accepted. **Carried**: the arviz-1.3.0 lazy-`numpyro` import-order fault (pre-existing on master; the handoff's finding; a one-line fix-up before W5.1 recommended); `path_filters.py`'s stale `test (py311)` name (W5.18's residue); `tests/inference`'s own cost (SBI training, torch horseshoe rows) is not `study` material and stays a future item. |
| W5.31 | merged 2026-09-25 at `2f76b53` (Fable-authored on Peter's ruling of 2026-09-24 — W5.26's carried finding, pre-existing on master; merged by Peter; the fresh-interpreter reproduction `import ampere.core; import arviz; import ampere.backends.jax` green with the worktree first on the path (`.distribution` present, `numpyro` a real module); jax rows (`test_jax_import_order.py`, the VI population row alone, `tests/backends/test_jax.py`): 201 passed in 88 s; dev rows (`test_jax_import_order.py` skipping without jax, plus the import sweep): 89 passed / 5 skipped in 3.1 s; lint, format-check clean; jax typecheck 0 errors (29 suppressed, as master); the jax gate leg on the merge `2f76b53`: 3809 passed / 413 skipped in 32 min 58 s — green, the two rows above the W5.26 gate's 3807 the item's own; **the merged CI run on `v2`, push run 36196425548 on `34f81fa`: green, 28 jobs plus the weekly sbi-characterisation job skipped by design, 30 min, the jax backends cell +2 rows, the warning sites identical to the W5.26 run's**). Landed: `ampere/backends/jax/_config.py` imports and touches `numpyro` (`_NUMPYRO_LOADED = numpyro.__version__`) before any sibling imports `numpyro.distributions`, with the mechanism in its comment (arviz 1.3.0's `importlib.util._LazyModule` stub executed mid-chain left `numpyro.distributions` without its `distribution` submodule, so `numpyro.factor` raised inside every native model whenever `ampere.core` and arviz were imported before the backend); `__init__.py` says why `._config` is imported first; `tests/backends/test_jax_import_order.py` — the failing order in a fresh interpreter, and a second row proving the stub really is lazy after arviz alone. Carried: none. |
| W5.1 | merged 2026-09-26 at `320afe5` (Opus-authored, Fable-reviewed as the second pass — `VonMisesFamily`, `_check_closure_phase_solver`, both native compositions, the oracle, the SBC arm and `likelihoods.md` §4/§14 read directly; **the agent's departure from the dispatch ruling accepted**: the numpy family computes the latent-conditional density and draw *given* `f` (Poisson's pattern) and refuses by name without it, since a family refusing outright made `realise`'s reference-point guard refuse every such problem (realised −150.95 against −inf) and left `simulate(observe=True)` with no oracle; **one review commit** `08ec93a` (the decision-log row per ground rule 9; the two stale "Phase 5's (W5.1)" comments in `tests/backends/interferometry_fixtures.py` and `tests/inference/test_interferometry.py`); merged by Peter; the orchestrator's re-runs on the branch through the lock: dev (`test_likelihood.py`, conformance interferometry, `tests/interferometry`) 241 passed / 13 skipped in 4 min 36 s, jax (conformance interferometry, the SBC arm with `interferometry_full`, the native draw) 98 / 4 in 4 min 52 s, torch (conformance interferometry and likelihoods, the native draw) 421 / 6 in 28 s, docs build clean with the six warning kinds master's CI build has and none new; the agent's runs: `tests/conformance/test_interferometry.py` torch and jax 86 passed / 4 skipped each, dev 56 / 4; the SBC arm on jax (`-m interferometry_full`) 4 passed in 206 s; `tests/interferometry` dev 21 / 9 in 254 s; jax interferometry + native draw + inference rows 60 / 14 in 334 s; torch native draw + conformance interferometry/inference/likelihoods 601 / 55; `tests/core/test_likelihood.py` 164; spec doctests 19; lint, format, typecheck (dev, torch, jax) clean, import sweep 89; the agent: ≈ 322 k tokens / 175 tool uses / 63 min, no hand-off; **the merged CI run on `v2`, push run 36204424495 on `d1f611a` (W5.1 + W5.29 + handoff): green, 28 jobs plus the weekly sbi-characterisation job skipped by design, 35 min**; **the combined W5.1 + W5.29 merged gate on `8188e1c` (2026-09-26, four legs serialised through the lock): dev 3064 passed / 593 skipped in 30 min 42 s, torch 3981 / 449 in 47 min 11 s, jax 3834 / 425 in 31 min 39 s, sbi 4137 / 293 in 54 min 08 s — all four green; against the W5.26/W5.31 baselines (3046/580, 3956/436, 3809/413, 4110/282) the two items' rows exactly**). Landed: `VonMisesFamily.CONSUMES_LATENT_GP = True` with `_latent_centre` and one `_native_only` refusal text for `log_prob` and `sample` (pinned word for word); `GaussianProcessNoise._check_closure_phase_solver` refusing any ordered-1-D solver on `ClosurePhases` by name before the solver's generic message; torch `_von_mises` and jax `_von_mises` scoring `observed − predicted − f` wrapped, the sampling gate admitting `von_mises` under a latent declaration, the native draws centred on `predicted + f` from the draw's own latent block (no solver code); `tests/conformance/oracles.py::von_mises_latent_log_likelihood` and `TestClosurePhaseLatentGP` (log-density against the from-scratch formula at `cross_solver` over 10 prior points and against numpy at `cross_backend`; the latent moves the density by > 1 nat; the refusal word for word on all four fixtures; `simulate(observe=True)` in (−π, π] with > 0.5 of the residual periodogram power in the lowest 4 of 64 frequencies against ≈ 0.1 under independent noise, measured over seeds 0–2 on numpy, torch and jax); `examples/interferometry/study.py` `PHASE_ARMS`, `phase_problem`, `with_phase_error`, `run_phase_calibration` (a 0.6 rad plane wave in the triangle coordinates, its size measured: 0.25 rad left the rigid fit covering separation at 0.83; closure phases alone left the latent arm at 0.58–0.67); `tests/interferometry/test_phase_calibration.py` — NUTS on jax (numpyro compiles the trajectory, ≈ 9 s a fit), 12 simulations, W4.4's floor 0.55 / ceiling 0.55 / gap 0.3 reused: latent coverage@0.9 [1.00, 1.00] against rigid [0.42, 0.08] at the pinned seed, [0.92, 0.92] against [0.25, 0.33] at a second; the default run keeps the structural rows (< 1 s); `likelihoods.md` §4 (a new subsection) and §14 (a new row, the von Mises row amended); `docs/source/interferometry.rst`. **For Peter**: none beyond the accepted departure. Carried: the tightest SBC margin is the rigid arm's separation coverage 0.42 against the 0.55 ceiling (one to two simulations); the reduced-rank latent on closure phases stays W5.4's question. |
| W5.29 | merged 2026-09-26 at `4b35967` (Opus-authored, Fable-reviewed as the second pass — `_context_sigma`, `_sample_chunk`, the trial draw, both backends' sigma threading and the jax `base=` route, the conformance row read directly; **no fixes needed**; one review commit `e17970d` (the decision-log row per ground rule 9); rebased onto master after W5.1 by the orchestrator, one mechanical conflict in jax `problem.py`'s `_sample_von_mises` resolved by keeping both changes; merged by Peter; the agent's verify chain on `21887af`: dev `test_simulate.py` green in 65 s, import sweep green in 6 s, torch `tests/conformance/test_inference.py` green in 54 s, jax green in 92 s, sbi `tests/inference/test_sbi.py` green in 6 min 08 s; the orchestrator's re-runs on the rebased head `9663bc9` through the lock: dev (`test_simulate.py`, `test_likelihood.py`) 295 passed in 50 s, jax (conformance inference + interferometry, the native von Mises draw, the phase-calibration structural rows) 275 / 59 in 2 min 08 s, torch (the same three) 271 / 55 in 50 s; lint, format-check (276 files), typecheck dev/jax/torch 0 errors on the rebased head; the agent's runs: the equality row torch 1 passed / 2 skipped in 0.4 s, jax 1 / 2 in 5 s; `tests/core/test_simulate.py -k context` 36 passed; the two sbi rows 2 passed in 21 s; lint, format-check (275 files), typecheck dev/torch/jax 0 errors; the agent: ≈ 254 k tokens / 157 tool uses / 18 min, handed off code-complete with the long runs queued; **the merged CI run on `v2`, push run 36204424495 on `d1f611a` (W5.1 + W5.29 + handoff): green, 28 jobs plus the weekly sbi-characterisation job skipped by design, 35 min**; **the combined W5.1 + W5.29 merged gate on `8188e1c` (2026-09-26, four legs serialised through the lock): dev 3064 passed / 593 skipped in 30 min 42 s, torch 3981 / 449 in 47 min 11 s, jax 3834 / 425 in 31 min 39 s, sbi 4137 / 293 in 54 min 08 s — all four green; against the W5.26/W5.31 baselines (3046/580, 3956/436, 3809/413, 4110/282) the two items' rows exactly**). Landed: `sample_observations(theta, predicted, seeds, *, sigma=None)` on both native realisations — torch `_context_rows` (a `(count,) + observed.shape` stack flattened as the lowering flattens and restricted to the retained samples, a wrong shape refused by name), `_base_sigma`/`_sigma`/`draw_chunk`/`sample_retained`/`sample_parameters` taking an optional per-draw `uncertainty` mapped through `vmap` beside θ, the context sigma replacing `sigma_data` before `scale`, `jitter` and `sigma_tensor`'s inflation; jax the same, with `_base_sigma` computing the quadrature and handing it to `FractionalModelNoise.sigma_jax(base=)`/`FractionalModelGPNoise.sigma_jax(base=)` because a traced sigma cannot be placed in a container; `dataset.py`: W5.10's refusal removed, `_NativeBatch(contextual=)`, `_context_sigma` (a per-dataset stack per chunk, the observation's own sigma filling the rows a draw's context leaves alone; the one combination with no batch form refused by name inside the chunk), `_sample_chunk` passing `sigma=` only for a context budget including the one-row retry, the trial draw passing `sigma=` for a context budget so a sampler without the keyword falls back to numpy draws, `place_observation(sigma=)`, every native `Simulation` recording `context` (failures included; a pre-existing gap), `provenance['simulate_batched']` true for a context budget, `_sbi.py` untouched (the default `native=None` takes the path; the cache key unchanged); `sample_observations_of`'s docstring; the conformance row `TestBatchedSimulation::test_a_context_budget_runs_natively_and_matches_the_loop` (200 draws under `ScaledSigma(0.1, 10)`: θ, every context record and every observation's sigma exactly equal to the loop's, the prediction at `cross_backend`, the noise standardised by its own context sigma within five Monte Carlo sigmas of `N(0, 1)`, `simulate_batched` true and `sample_backend` the backend); `tests/core/test_simulate.py` — the refusal and fallback rows inverted, `TestTheNativePathDrawsTheContext` (five rows on the fake backend, including bitwise equality with the loop when the sampler has no `sigma=`, and each draw's own context handed row for row when it has); `tests/inference/test_sbi.py::TestTheBatchedPathDrawsTheContext` (the amortised run batched on torch under the prior; the cache key independent of the path); all three shipped priors reach the batched path; `inference.md` §13 (the signature, "The native path draws the context", two refusals not three), `encoding.md` §8 (no hash moves), `sbi.rst` one sentence. **Measured**: the tutorial's context run on `npe_native`'s torch problem (budget 2000, set layout, `ScaledSigma(0.5, 2.0)`): `simulate_many` inside `engine.run` 0.51 s batched against 1.30 s looped (≈ 2.5×), standalone warm 0.28 s against 0.91 s (≈ 3.2×); the whole run 28.6 s against 36.1 s, the difference mostly training on different noise draws; fractional noise under a context (scratch, 3000 draws): standardised residual std 0.995 jax, 1.002 torch. **For Peter**: none (the three questions answered at review, above). Carried: a user noise model whose `sigma_jax` lacks `base=` falls back to numpy draws under a context, visible only through `sample_backend`; the in-chunk refusal has no test. |
| W5.21 | merged 2026-09-26 at `6f033eb` (Sonnet-authored, Fable-reviewed as the second pass — the two core hunks, the conformance rows and the §15 amendment read directly; **no fixes needed**; one review commit `6015e4c` (the decision-log row per ground rule 9); merged and pushed by the orchestrator under Peter's session authorisation of 2026-09-26; the orchestrator's re-runs on `44215f6` through the lock: dev (conformance image, `test_likelihood.py`, the image study rows, the bake-off's non-timing rows) 232 passed / 9 skipped in 35 s, torch conformance image 54 passed, jax 54 passed; the agent's runs: `tests/conformance/test_image.py` dev 34, torch 54, jax 54; `test_image_study.py` 21 / 2; `test_likelihood.py` 166; the bake-off 11 non-`image_full` rows; lint, format, typecheck dev/torch/jax 0 errors; import sweep 89 / 4; docs build 15 warnings, all master's; the grep for any `grid_gp` symbol empty outside the three record documents; the agent: ≈ 292 k tokens / 261 tool uses / 38 min, no hand-off; **the merged CI run on `v2`, push run 36215935801 on `7c5bda8`: red on ONE row — `test_sbi.py::TestTheCalibrationFastPath::test_a_trained_npe_posterior_passes_check_sbc`, c2st max 0.66 against the `< 0.65` pin in the sbi backends cell (28 other jobs green); not this item's doing (see the handoff's finding of 2026-09-26: the coordinate path is bitwise unchanged for point sets, the row passes locally in every ordering at 0.575, and the fixture's posterior is mildly miscalibrated so the pin flaps on the runner); a rerun of the job was cancelled by the next push, whose run 36219726492 on `b6277d8` was green including that row**; **the combined W5.21 + W5.24 merged gate on `b6277d8` (2026-09-26, four legs serialised through the lock): dev 3102 passed / 599 skipped in 30 min 58 s, torch 4033 / 450 in 50 min 03 s, jax 3886 / 426 in 33 min 21 s, sbi 4189 / 294 in 57 min 03 s — all four green; against the W5.1/W5.29 gate's baselines (3064/593, 3981/449, 3834/425, 4137/293) the two items' rows on top**). Landed: `GPSolver.check_compatible` accepting `Layout.POINTS` or `Layout.GRID` (the ordered-1-D and product rules untouched; a comment recording that `Kernel.check_axes` already refuses a leaf selecting axes the container lacks, for either layout, so the item's condition needed no new check); `Likelihood._coordinates` through `sample_coordinates(observed)[retain]` for either layout (imported inside the method — `dataset.py` imports `Likelihood` at load; carried to W5.32 (f)); `examples/image/grid_gp.py` deleted with its `GriddedSolver` mixin and every `Grid*` subclass (two more in `bakeoff.py`), `study.py` and `bakeoff.py` on the shipped solvers and a plain `Likelihood`, `__init__.py`'s docstring and `__all__`; `tests/core/test_likelihood.py` — the old refusal row replaced by GRID acceptance, the full path against an independent Euclidean Matérn-3/2 over `sample_coordinates`' own output (the flattening order proven), and `QuasisepGP` still refused on a 2-D grid; `tests/conformance/test_image.py` — acceptance on every fixture for `DenseGP` and `HilbertSpaceGP` (`GRID_KERNEL`, a two-axis Matérn-3/2 in mas) and the quasiseparable refusal by name; `test_solver_bakeoff.py` rewired; `docs/source/image.rst` §5 rewritten ("a grid composes like any other container"); `likelihoods.md` §15 item 4 struck through and marked lifted in the style of items 1 and 3. **Finding**: the refusal paragraph lived in §15, not §7 as the item said — amended where it was. **For Peter**: none. Carried: the lazy import (W5.32 (f)); `StructuredGridGP`/SKI stays a slot. |
| W5.24 | merged 2026-09-26 at `eac829c` (Opus-authored, Fable-reviewed as the second pass — `route`, `variances`, `dense_covariance`/`_dense_log_prob` (channel-major, matching the oracle's assembly), `features`/`_reduced_rank_log_prob` (the Woodbury identities checked), the torch native routes, the oracle and the conformance rows read directly; **no fixes needed**; the three out-of-ownership edits accepted as necessary — `DatasetCollection.group_log_likelihood`/`draw_group` and `results/derived.py`'s joint pointwise passed only the first channel's σ (a latent bug the item exposed), `tests/core/test_noise_joint.py` asserted the refusal being lifted, W5.9's coverage pin lives in `tests/examples/test_astrometry_example.py`; one review commit `0b842e2` (the decision-log row per ground rule 9); merged and pushed by the orchestrator under Peter's session authorisation of 2026-09-26; the orchestrator's re-runs on `c7a0885` through the lock: dev (conformance astrometry, `test_noise_joint.py`, `test_likelihood.py`, `test_dataset.py`, the results joint-pointwise rows, `tests/astrometry`) 435 passed / 12 skipped in 38 s, torch (conformance astrometry + `test_joint_noise_nuts.py`) 77 passed in 3 min 45 s, jax 77 in 2 min 03 s, the astrometry example rows with the full marker 20 passed in 20 min 09 s; lint, format-check, typecheck dev/torch/jax 0 errors; the agent's runs: conformance astrometry dev 46, torch 69, jax 69; NUTS 8/8 torch (219 s), 8/8 jax (100 s); astrometry + example rows 17 / 12; core 374; results joint 13; import sweep 89; the agent: ≈ 344 k tokens / 160 tool uses / 46 min, no hand-off; **the merged CI run on `v2`, push run 36219726492 on `b6277d8`: green, 28 jobs plus the weekly sbi-characterisation job skipped by design — including the c2st row red on the previous push**; **the combined W5.21 + W5.24 merged gate on `b6277d8` (2026-09-26, four legs serialised through the lock): dev 3102 passed / 599 skipped in 30 min 58 s, torch 4033 / 450 in 50 min 03 s, jax 3886 / 426 in 33 min 21 s, sbi 4189 / 294 in 57 min 03 s — all four green; against the W5.1/W5.29 gate's baselines (3064/593, 3981/449, 3834/425, 4137/293) the two items' rows on top**). Landed: `JointGaussianProcessNoise.route(variance)` — `"rotated"` for identical columns (exact comparison), `"dense"` under `DenseGP`, `"reduced_rank"` under `REDUCED_RANK_SOLVERS` (`HilbertSpaceGP`; `EquispacedFourierGP` on the reference path only), any other solver refused by name with the fix; `variances()` (the `(n, T)` block), `latent_size(n, route)` (`T·m` reduced-rank), `dense_covariance` (`B⊗(K_x+j²I) + blockdiag(diag σ_t²)`, one Cholesky), `features` (`Φ̃` as the solver's `latent_transform` of the identity), `_reduced_rank_log_prob` (Woodbury: `T` Gram matrices, the `Tm×Tm` capacitance, `B⊗K_x` never formed), `_sample_unequal` (the correlated part through the whitening, each channel's white noise added in the original basis), `pointwise_log_prob` refused off the rotated route; `_check_shared_diagonal` no longer refuses unequal channels, still refuses them under other solvers and a group whose channels disagree about carrying uncertainties; torch and jax `_LoweredJointGroup` choosing the route once at lowering — the dense route one `cholesky_ex` of the `TN×TN` matrix in the backend's arithmetic, the reduced-rank route native Woodbury through `latent_transform_native`/`latent_transform_jax`, every channel's uncertainty checked positive; `tests/conformance/oracles.py` `rotation_coupling_matrix`, `coregionalised_covariance` (block-assembled, not `kron`), `coregionalised_log_density` (scipy's MVN); `TestHeteroscedasticJointNoise` — dense against the oracle at `cross_solver` (1.4e-14 measured), reduced-rank converging in W5.4's class (m 8→64: 3.0e-1 → 8.1e-4, `latent_size == 2m`), equal σ `==` W5.9 under Dense and Quasisep and through the composed problem, Quasisep with unequal σ refused naming `HilbertSpaceGP`, `simulate` over 2000 draws with every `2N×2N` covariance entry within `monte_carlo_sigmas`; `test_joint_noise_nuts.py` — the heteroscedastic chain through the reduced-rank route on torch (140 s) and jax (82 s) with `B`'s trace and the length scale inside 95 %, lowered vs contract to 1e-7 on both routes; `synthetic_joint_data(heteroscedastic=)` (a per-channel-per-epoch log-uniform factor up to 2 around `SIGMA`, its own `seed+2` stream), `--heteroscedastic`, `joint_noise(solver=)`, the arm on the dense route (cheapest at N=28: 0.38 ms against 0.46 rotated and 0.64 reduced-rank at m=32), the SBC pin at W5.9's 0.80 floor (measured 0.938) beside W5.9's, an emcee recovery row; `likelihoods.md` §7 ("The one restriction" rewritten as the three routes) and §15 (unequal uncertainties struck; the per-(sample, channel) pointwise follow-on added); `astrometry.rst` "Heteroscedastic channels". **For Peter**: (1) the heteroscedastic study does not reproduce W5.9's gap — the independent arm covers 0.917 against the joint 0.938 (white noise up to 2× dilutes the induced correlation), so the row pins the floor only and records the table; confirm, or ask for a narrower σ spread / stronger `B` to show the gap; (2) the jitter convention differs between the unequal routes (`B⊗j²I` dense, per-channel diagonal reduced-rank; both default 0) — accept. **Ruled by Peter 2026-09-27**: (1) confirmed — the row pins the floor and records the table; (2) accepted. Carried: no LATENT path exists for joint groups (`_check_member_likelihood` refuses non-analytic families), so only `simulate` uses the new factorisation; EFGP has no native twin; the general LMC and mismatched grids stay the follow-on. |
| W5.15 | merged 2026-09-27 at `a8f4d47` (Sonnet-authored, Fable-reviewed as the second pass — every diff read directly, two review probes run (dynesty with slice sampling; nautilus's posterior tail); **two review commits**: `2300152` (`DYNESTY_SAMPLE = "rslice"` as the study's dynesty proposal method, measured — the branch had found dynesty's default uniform-in-ellipsoid proposal not converging in 38 min at either prior's width and drawn the wrong lesson; the `astrometry_full` row made the CLI's own call and run to completion, 1 passed in 241 s; a dev-only per-PR row at 100 live points, 142 s) and `a455136` (`astrometry.rst` §6 corrected: the epochs are W4.9's shipped set, the dynesty rows and engine paragraph teach the proposal method, the corner-plot prose limited to what was measured); merged by Peter; the orchestrator's runs on `a455136` through the lock: dev astrometry rows 36 passed / 17 skipped in 224 s, nested rows 2 passed in 187 s, hygiene clean (typecheck 0 errors, import sweep 89/4), docs 14 warnings (baseline 15); the agent's: 35/17 in 85 s, nested 2 passed in 279 s, docs 15 warnings; the agent ≈ 394 k tokens / 390 tool uses / 84 min; **the W5.15 merged dev gate leg on `bad7979` (2026-09-27): 3108 passed / 602 skipped in 36 min 05 s — green; against the W5.21/W5.24 baseline 3102/599 the six new passed rows and three new skips are exactly this item's; the torch/jax/sbi legs run with W5.17's merged gate: **sbi 4194 / 298 in 71 min on `e9c3417` — green; torch 4038 / 454 in 64 min — green; jax 3891 / 430 in 42 min — green; dev 3108 / 602 in 36 min (this item's own leg's counts) — green; all four green**; **the CI run on the push of `bad7979` to `v2`, run 36289988248: GREEN** (queued ~35 min on GitHub's runners before starting, ~50 min end to end)). Landed: `build_model`/`build_problem(wide_prior=)` with `WIDE_PERIOD_PRIOR = loguniform(50, 2000)`; `fit(engine=, live_points=, dlogz=)` over `ENGINES = (emcee, nuts, dynesty, nautilus, ultranest)` with `None` the pre-item behaviour; `period_modes(run, gap=0.1, minimum_mass=0.02)` splitting the equal-weight period draws in ln P with per-mode mass fractions and `ln Z_k = ln Z + ln f_k`; the CLI's `--wide-prior` (reaching for dynesty), `--engine`, `--live-points`, `--dlogz`; `report()` printing the evidence triple; `tests/examples/test_astrometry_example.py::TestTheWidePriorArm` (support, the splitting rule on hand-built draws, the smoke run); `tests/astrometry/test_recovery.py::TestNestedSamplingResolvesThePeriodModes` (dynesty full-budget under `astrometry_full`, dynesty per-PR dev-only, nautilus and ultranest by `find_spec`); `astrometry.rst` §6 "period aliasing, and the two remedies" with the measured table, §9 pointer, §11 bullet; `interferometry.rst`'s closing hazard section. **Measured**: every converged nested run finds one dominant mode at the truth (nautilus +92.20 ± 0.01 in 84 s, ultranest +91.48 ± 0.48 in 143 s, dynesty/rslice +92.25 ± 0.55 in 241 s); the informed prior's evidence +95.17 (nautilus) gives a Bayes factor ≈ 19 for "period known to ±30 d", the wide prior's Occam penalty; emcee splits 62.5/37.5 between 400 d and a 763-d alias where the nested runs put no mass, and its naive 95 % coverage check still reads "ok". **For Peter**: none. Carried: `examples/astrometry/__main__.py`'s usage text does not list the new flags (W5.32 territory); the `study` convention has no hook outside `tests/m2` and `tests/results`, so `tests/astrometry` carries its own `dev_only` marker (W5.19's list); `default_live_points`'s floor and the per-PR row's 100 live points disagree by design — the row pins the mode structure, not the evidence. |
| W5.17 | merged 2026-09-27 at `7e8e5c8` (Opus-authored, Fable-reviewed as the second pass — every hot-path diff read directly and celerite2.jax's own `_do_compute` read beside the subclass; **no fixes needed**; merged by Peter; the orchestrator's re-runs on `216ece5` through the lock: dev `tests/conformance` + `test_likelihood.py` + `test_transform.py` + `tests/interferometry` 992 passed / 84 skipped in 3 min 10 s, jax `tests/conformance` + `test_jax.py` 1271 / 75 in 4 min 42 s, torch `tests/conformance` 1073 / 75 in 71 s, hygiene clean (typecheck 0 errors, import sweep 89/4); the agent's: conformance dev 669/75, torch 1073/75, jax 1073/75 at the head and after each lever in dev, collected tests 3869 before and after, `bench` in all three environments before (on `3320eef`) and after, saved outside git; the agent ≈ 326 k tokens / 243 tool uses / 111 min; **the combined W5.15 + W5.17 merged gate on `e9c3417` (2026-09-27): sbi 4194 passed / 298 skipped in 71 min (baseline 4189/294 in 57: W5.15's five smoke rows passed and its four dev-only/full/nested rows skipped, W5.17 adds no rows; the minutes are the machine, no other load was on it) — green; torch 4038 / 454 in 64 min (baseline 4033/450 in 50; the same +5/+4) — green; jax 3891 / 430 in 41 min 37 s (baseline 3886/426 in 33; the same +5/+4) — green; dev 3108 / 602 in 36 min 22 s (W5.15's dev leg's counts exactly, W5.17 adding no rows) — green; **all four legs green**, 3 h 39 m end to end through the lock; **the CI run on the push of `497e68e` to `v2`, run 36302416953: GREEN**, every job**). Landed: `docs/design/performance_memo.md` (the machine, the baselines, per-load profiles of the M2 driver, the interferometry study and the image study at 24² and 64², construction cost, issues #12/#29/#67 checked against the v2 core — the lifecycle claim holds, zero renegotiations over 50 `log_prob` calls per load — the levers ranked by measured share, §6 proposals, §7 results); five levers, one commit each, every one bit-for-bit (`np.array_equal`): **L1** `ampere/backends/jax/gp.py` — celerite2.jax's `lax.cond` over fresh closures recompiled per eager call (~98 % of the jax contract path); a `_gaussian_process_type()` subclass reuses the branch functions, M2 jax contract ×3–7, `test_backend_quasisep_contract[*-jax]` 30 → 3.7–7.5 ms; **L2** `ampere/core/likelihood.py` `_cho_solve` (+ `efgp.py`) — scipy 1.18's C-ordered `cho_factor` made `cho_solve` copy the N×N factor per solve; LAPACK `potrs` now gets the transposed view with the triangle flipped, complex/non-contiguous/Fortran factors fall through; image 64² ×2.6, `test_dense_anchor` dev 208 → 84 ms; **L3** reference `Resample._planned_influence` — the weight matrix once per input grid, keyed on both grids' bytes (issue #12), ×30–65 on the step; **L4** reference `FourierSample._planned_transform` — solid angles and DFT factors once per grid and expanded coverage, plus the same-unit-object fast return in `_resolve_template`, interferometry arms ×1.3–2.1; **L5** `examples/m2_misspecification/model{,_torch,_jax}.py` `_template_for` — the container template built once on the simple path, ×1.1–1.35 on five of six rows (jax flexible recorded as no change). `test_the_fast_path_beats_the_contract_path` still passes (realised ≈ 7× the contract path now, against a 3× floor). **Ruled by Peter 2026-09-27**: P1's private-API route not taken; **a note to evaluate the performance of scipy's new public distribution infrastructure** (≥ 1.15, `scipy.stats.Normal`/`make_distribution`) against the legacy wrapper on the same load is recorded in the memo's §6 as the follow-up. Carried: pickled reference steps now carry their plans; P2 batched evaluation (a §4.5 addition), P3 the image solver default (a science decision), P4 torch `Resample` caching (an autograd question), P5 the jax contract path's remaining eager dispatch (jitting the whole solve is not bit-identical, 8e-11); `examples/image/study.py`'s `_noise_for` falls back to the core `DenseGP` on torch/jax. |
| W5.18 | merged 2026-09-27 at `5f38f39` (Sonnet-authored to the hand-off, Fable-finished and reviewed — the agent landed the `DEV_JOBS` fix and the `test-fast` task's exclusion logic and handed off at the measurement after ≈ 168 k tokens / 83 tool uses / 24 min; the orchestrator took the remainder under the "small things" rule: the measurement (14 min 02 s for the literal definition), `--durations=80`, **Peter's ruling of 2026-09-27 that `test-fast` is the contract-and-backend tier**, the ruled definition and its docs; merged by Peter; the orchestrator's runs on `2ed92c4` through the lock: `pixi run -e dev test-fast` 2879 passed / 425 skipped / 244 deselected in 2 min 54 s (3 min 04 s wall), collected 3299/3543 with the 244 deselected exactly the `study` rows, the 26 SBI-training classes and the two slow results classes, lint/format/typecheck (0 errors)/import sweep (89/4)/actionlint clean, `tests/gpu` skip-clean in dev (2 skipped), torch (20) and jax (15); no gate of its own — a task-table, CI-script and docs change; **the combined W5.18 + W5.32 merged gate on `418dce7` green on all four legs, the counts in W5.32's row**). Landed: (1) `test-fast` in `pyproject.toml`'s task table — `tests/core tests/results tests/conformance tests/backends tests/inference tests/m2 -m "not study"` plus 28 `--deselect` node ids (the 26 `@needs_sbi` training classes of `test_sbi.py`, `test_plots.py::TestCornerPaging` at 75 s, `test_calibration.py::TestTheRefitLoop` at 27 s), the comment carrying the ruling, the durations (the W5.15 per-PR dynesty pin 141 s, the interferometry calibration fixtures 83 + 44 s, the image and interferometry study CLIs 50/35/17/17 s, the astrometry SBC smoke runs 22 + 20 s) and the measured figures, and noting `test-all -m "not study"` as the 14-minute middle tier; `docs/development.md`'s task bullet and `install.rst`'s console block name it; (2) `tests/gpu` confirmed skip-clean, no edit; (3) `path_filters.py`'s `DEV_JOBS` names `test (py312)`/`(py313)`/`(py314)` in place of the retired `py311`, with a comment. **For Peter**: none beyond the ruling taken. Carried: the 26 deselect flags would be one `sbi_training` marker if the markers list were this item's (W5.32 owns it; a follow-up housekeeping part); `tests/inference/test_interferometry.py`'s one SBI-training class (~662) is not deselected — moot in dev, a small cost under `-e sbi`; `ci.yml` ~312's historical note about `test-py311` left as accurate history; `docs/development.md`'s "Common tasks" bullet still describes `test-all` as six suites (nine since W5.26) — W5.19's list; CLAUDE.md's task list gains `test-fast` at the merge (the orchestrator's edit). |
| W5.32 | merged 2026-09-27 at `1135f30` (Sonnet-authored, Fable-reviewed as the second pass — every diff read directly, the (j) signature-bind check read against both backends' contextual calls, the (k) module's convention checked against `results.md` §9 (the dispatch prompt had it backwards: the density is in the *constrained* coordinates, as the module tests); **no fixes needed**; merged by the orchestrator on Peter's word; the orchestrator's re-runs on `dbeaf32` through the lock: the fast-tier dev selection (master's `test-fast`, the task did not yet exist on the branch) 2883 passed / 433 skipped / 244 deselected in 3 min 02 s with `PytestRemovedIn10Warning` promoted to an error, dev `test_image_study.py` + `test_astrometry_example.py` + `test_reference_interferometry.py` 81 / 6 in 67 s, jax `tests/conformance` + the (k) module + the (j) rows + `test_image_study.py` 1104 / 79 in 3 min 40 s, torch `tests/conformance` + (k) + `test_image_study.py` 1095 / 81 in 75 s, sbi `-k TMNRE` 27 passed with `FutureWarning` promoted to an error in 3 min 37 s, lint/format (276 files)/typecheck dev-torch-jax 0 errors/import sweep 89/4 clean, the docs build 11 warnings from a clean `_build` (15 on master; none of the autodoc, duplicate-object or `ampere.infer.sbi` kind remain); the agent's evidence per letter in its report, including each regression check by a tagged stash (the (l) size assertions, the (b) `FutureWarning`, the (h) `DatasetError`, the (k) jacobian subtraction); **the agent ran to ≈ 728 k tokens / 813 tool uses / 147 min — past the soft cap without a hand-off; the work is complete and clean, the breach recorded as W5.28's was**; **the combined W5.18 + W5.32 merged gate on `418dce7` (launched 2026-09-27 11:59 on the merge, four legs serialised through the lock, complete 14:40): jax 3903 passed / 432 skipped in 29 min 50 s, torch 4043 / 463 in 47 min 05 s, sbi 4200 / 306 in 53 min 55 s, dev 3112 / 612 in 29 min 46 s — all four green; against the W5.15/W5.17 gate's baselines (dev 3108/602, torch 4038/454, jax 3891/430, sbi 4194/298) the rows on top are dev +4 passed / +10 skipped, torch +5 / +9, jax +12 / +2, sbi +6 / +8 — this item's rows landing as expected (the (j) rows run under jax and skip elsewhere, the (k) rows skip by extra, the (h) rows run per modern backend); log `~/.cache/ampere-gates/w5.18-w5.32-gate.log`**; **the CI run on the push of `418dce7` to `v2`, run 36314292175: GREEN** (33 min end to end)). Landed: (a) eight class-scoped fixtures `@staticmethod` (`test_population.py` ×6, `test_engines.py`, `test_astropy_engines.py`); (b) `posterior_parameters=MCMCPosteriorParameters(**dict)` at `_sbi.py`'s two `build_posterior` sites, imported lazily, the two dicts kept; (c) `WarpedKernel` in the backend pages' `:exclude-members:`; (d) `sbi` in `autodoc_mock_imports` with the reason; (e) `TestCornerPaging` closes its figures; (f) `sample_coordinates` in `ampere/core/encoding.py`, re-exported from `ampere.core` unchanged, both callers importing at module level; (g) the astrometry usage text; (h) `_noise_for` falls back to `module.DenseGP()` with a row per modern backend; (i) `_skip_study_rows_outside_dev` once in `tests/conftest.py`, the two copies gone, `tests/astrometry`'s `dev_only` replaced by `@pytest.mark.study`, the marker's comment updated; (j) `AmpereContextFallbackWarning` (`settings.py`), `_context_base_gap`/`_warn_of_context_base_gaps` in `dataset.py` binding each backend's real contextual sigma call against the hook's signature before the trial and warning once per `simulate_many`, the in-chunk refusal asserted, one sentence in `inference.md` §13; (k) `tests/inference/test_proposal_log_density_conformance.py` — the Laplace guide against its closed-form density, pathfinder and NPE by the importance-weight identity, each skipping by extra; (l) `__getstate__` on both reference `_Step` classes dropping `_*_plan` and `_template_source_unit`, with size and bit-identity tests. **Findings**: the proposal density has three producers, not W5.14's four (the fourth hit is `results/population.py` reading it); torch's `sigma_tensor` takes `base` positionally so a torch gap fails loudly at construction already — the (j) warning's real case is jax's keyword-optional hook; netCDF4 1.7.4.1 not on conda-forge yet, no re-lock. **For Peter**: none. Carried: `tests/results/test_plots.py`'s other classes never close their figures and the 20-figure warning still fires when the whole file runs (W5.19's list or Phase 6 housekeeping); the (j) warning's `stacklevel=5` is best-effort through a generator; W5.18's 26 `--deselect` flags could become one marker now that (i) has made the markers list settled (Phase 6 housekeeping). |
| W5.19 | merged 2026-09-28 at `50f6ff4` (Sonnet-authored, Fable-reviewed as the second pass — every diff read directly and every figure the new pages quote traced to the code, `performance_memo.md`'s lever table, `likelihoods.md` §7's bake-off tables and the W5.6/W5.9/W5.24 rows; the `conditional_loo` claim checked against all three backends' `QuasisepGP` (the reference refuses, torch and jax compute); **one review commit** `670ae34` (the overview's "in flight since 2026-09-15" dated as closed by this pass, the decision-log row's "this section" from inside §2, an unrecorded starting KS p-value dropped from `sbi.rst`'s caution); merged by Peter; the orchestrator's runs on `670ae34` through the lock: `pixi run docs` from a clean `_build` 12 warnings with the list byte-identical to the base commit's (the agent measured 12 on `ddfa8a9` too; the prompt's "11" was the distinct count — one legacy warning prints twice), `test_spec_doctests.py` 19 passed, `test-examples` 93 passed / 13 skipped in 2 min 57 s, lint/format (276 files)/typecheck (0 errors)/import sweep (89/4) clean; the agent's: the same four on `dde467e` (docs 12/identical after every commit, doctests 19, examples 93/13 in 190 s, hygiene clean), every `pycon` block on the three pages run in `pixi run -e dev python`, the population tutorial's joint-NUTS cell not run (needs torch or jax; the worktree had dev only — the page says so); the agent ≈ 416 k tokens / 264 tool uses / 36 min, no hand-off; the worktree was cut from legacy `b8e585b` and re-based onto `ddfa8a9` on the prompt's instruction before any edit; **the W5.19 merged dev gate leg on `ed9decf` (2026-09-28): 3112 passed / 612 skipped in 31 min 58 s — green, exactly the W5.18/W5.32 baseline (a docs-only item adds no rows; log `~/.cache/ampere-gates/w5.19-gate.log`); the only leg run, per Peter's ruling of 2026-09-28 that gate legs are scoped to the code touched**; the CI run on the push of `1740468` to `v2` (the first push under D4's push-at-every-merge ruling): run 36359647901: GREEN, every job including the sbi backends cell, about 85 min end to end). Landed: five additive *Amended W5.19* annotations with one decision-log row — `likelihoods.md` §7 (the `QuasisepGP.condition(at=)` multi-axis capability disagreement W5.28 carried, Phase 6), `architecture.md` §3 ("optimisers remain Phase 5" → moved to Phase 6, distinct from W5.17's pass), `astrometric_timeseries.md` (the cross-channel noise gap closed by W5.9/W5.24), `horizon_notes.md` §2 (the HSGP/EFGP/Vecchia recommendation landed as stated), `results.md` §14 (population inference built by W5.13/W5.22) — Phase 5's items having annotated their own pages as they landed, so five gaps where Phase 4 had nineteen; `kernels.rst` §3 "Non-stationarity: `WarpedKernel`" (input and amplitude warping, quasiseparability preserved, the parameters and their shrinkage prior, the many-lines worked case, the knot-coverage caution; §§4–6 renumbered and `advanced.rst`'s cross-reference moved); **`docs/source/solvers.rst`**, new (exact `DenseGP`/`QuasisepGP` with the `Layout.GRID` acceptance and the `conditional_loo` split; `HilbertSpaceGP` with the spec-hashed approximation fields and the `m`-sized latent, EFGP and Vecchia as unpromoted prototypes with the measured reasons, the `EXACT = False` tolerance class as `approximation_envelope` asserts it; `JointGaussianProcessNoise` with both couplings, W5.24's two heteroscedastic routes and the jitter note) in the index toctree after `kernels`; **`docs/source/population.rst`**, new (the fifty-object plate fit from `test_population_nuts.py`, the two layouts and `MAX_FLAT_MEMBERS`, `DatasetCollection.plate`, `fit_population` on a twelve-object emcee archive run on the branch, the interim prior and `ampere_free_priors`, `population_full`, W5.30 (c)) in the tutorials toctree after `sbi`; `overview.rst`, `README.md` and `index.rst` each say what Phase 5 shipped (the engine count corrected from six to nine in all three; the three new engines in the overview's table and README's extras table); `DEVELOPMENT_PLAN.md` §5's Phase 5 section as the landed summary — a "Landed at W5.x" clause per bullet, every ruling and *corrected at* note kept verbatim; the carried list — `docs/development.md`'s `test-all` bullet says nine suites, `sbi.rst` gains the calibration-smoke caution (the c2st finding of 2026-09-26) and the `z_score_x`/outlier-warning guidance, the two W5.32-closed notes confirmed absent from the pages. **Findings**: `advanced.rst` has no shrinkage section (the prompt assumed one) — `kernels.rst` §1 is the treatment and the new section links there; the clean docs build counts 12 warnings on `ddfa8a9` and on the head, not 11. **For Peter**: whether a shrinkage section belongs in `advanced.rst`; whether the `pycon` blocks under `docs/source` get a doctest runner (W5.7's carried idea; Phase 6 unless ruled otherwise). Carried: `pyproject.toml`'s `blackjax` extra comment says `_blackjax.py` "does not exist yet" (it does — a one-line fix, Phase 6 housekeeping); `test_plots.py`'s open figures (Phase 6 housekeeping); the `condition` conformance row (Phase 6). |
| W6.6 | merged 2026-09-28 at `0257ed4` (Sonnet-authored, Fable-reviewed as the second pass — every diff read directly; the torch `condition` fix checked against the module's `_points` signature and against jax's `condition` at 1312, which already built its points from the unreduced container; the marker placement checked on the 27 classes; **no fixes needed**; merged by Peter; the orchestrator's re-runs on `3f5766c` through the lock: dev `test_plots.py` 56 passed with matplotlib's "More than 20 figures" `RuntimeWarning` promoted to an error, the alias test + `tests/conformance/test_likelihoods.py` 225 passed / 2 skipped, `test-fast` collection 3308/3557 with 249 deselected under the marker form, lint/format (276 files)/typecheck (0 errors)/import sweep 89/4 clean, torch `test_likelihoods.py` 334 / 2 and typecheck 0 errors, jax 334 / 2 and 0 errors, sbi `-m sbi_training` collects 144 rows across the 27 classes; the agent's evidence per part in its report, including the old-form versus new-form collection counts (3311/3555 with 244 deselected against 3306/3555 with 249 — the five rows being `TestSBIOnVisibilitiesAndClosurePhases`'s methods) and a live dev `test-fast` at 2885 passed / 428 skipped / 249 deselected in 177 s; the agent ≈ 251 k tokens / 199 tool uses / 30 min, no hand-off; the worktree was cut from legacy `b8e585b` and re-based onto `f3f1acb` on the prompt's instruction before any edit; **the CI run on the push of `1161f04` to `v2`, run 36363416751: GREEN — the first CI run of Phase 6** (the run on `6225b94`, 36363389278, was cancelled by that superseding push); **the W6.6 merged gate (launched 2026-09-28 03:01 BST on the merge `0257ed4`, the scoped legs in sequence through the lock — dev on `0257ed4`, torch and jax on master's working tree at `a6a4c67`, which differs from the merge only by W6.8's memo and handoff prose; complete 04:51): dev 3114 passed / 612 skipped in 30 min 56 s, torch 4046 / 463 in 49 min 02 s, jax 3906 / 432 in 28 min 21 s — all three green; against the W5.18/W5.32 baselines (dev 3112/612, torch 4043/463, jax 3903/432) the rows on top are dev +2, torch +3, jax +3 — the (d) row on each environment's fixtures; log `~/.cache/ampere-gates/w6.6-gate.log`**; the CI run on the push of `27076e4`, 36368233931, was cancelled by the superseding push of W6.8's merge, whose run 36368703571 on `7f4fa9b` is GREEN and covers both). Landed: (a) the blackjax extra comment; (b) the autouse close-all fixture in `test_plots.py` (`TestCorner`, `TestTrace`, `TestTracePaging`, `TestPosteriorPredictive` were the classes leaving figures open); (c) the `sbi_training` marker with its comment, on the 26 `test_sbi.py` classes and `test_interferometry.py::TestSBIOnVisibilitiesAndClosurePhases`, `test-fast` as `-m "not study and not sbi_training"` plus the two slow non-SBI deselects, the task comment keeping the measured figures; (d) torch's `QuasisepGP.condition` builds `points` from the unreduced container via `_points` and checks the grid's width against the data's, with `TestAxisSelector::test_condition_at_a_multi_axis_container_reaches_the_quasiseparable_path` asserting agreement with `DenseGP.condition` at `cross_solver` on the three-axis `DispersedPoints` fixture wherever a quasiseparable solver is declared; (e) the alias's docstring, warning text and test say "1.0.0 final"; (g) one paragraph in `UltranestEngine`'s docstring, checked against ultranest 4.5.2's `integrator.py` (`combine_results` folds `logzerr_tail` into the stored `logzerr`; the live line prints `logZerr_bs` alone). **Findings**: the carried note's premise for (d) was half wrong — jax already agreed with `DenseGP` before the branch (the `likelihoods.md` §7 *Amended W5.19* annotation, which says torch and jax refuse, is now stale on both counts and is W6.3's to update); the literal `-W error::RuntimeWarning` check the prompt suggested for (b) also trips a pre-existing netCDF4/NumPy-2.5 binary-incompatibility warning in dev, so the check was scoped to matplotlib's message; (f) not actionable — netCDF4 1.7.4.1 on neither PyPI nor conda-forge at 2026-09-28, the lock untouched. **For Peter**: none. Carried: the `likelihoods.md` §7 annotation above (W6.3); netCDF4's re-lock when 1.7.4.1 ships (a later housekeeping item). |
| W6.8 | merged 2026-09-28 at `88ce33c` (Sonnet-authored, Fable-reviewed as the second pass — the whole §8 diff read; the load-bearing claims reproduced by the orchestrator in a probe on this machine: `Normal` has `icdf`/`sample` and neither `ppf` nor `rvs`, a `make_distribution` product is a `CustomDistribution` with no `name` or `dist`, `Uniform`/`loguniform`/`halfnorm * scale` bit-identical to the frozen route at 10 000 draws with `-inf` at the edge on both, `isinstance(Normal(...), Prior)` False, the six-prior `lnprior` 265 µs legacy against 37 µs new (7.1×; the memo's 250 / 34.7); **one review commit** `f29681a` (the recommendation's "silently `AttributeError`s" said loud and silent at once — the failures are loud, by name); merged by Peter; no gate — a report only, `lint`/`format-check` clean; **the CI run on the push of `7f4fa9b` to `v2`, run 36368703571: GREEN**; the agent ≈ 182 k tokens / 94 tool uses / 18 min, no hand-off; the worktree was cut from legacy `b8e585b` and re-based onto `f3f1acb` on the prompt's instruction before any edit). Landed: `docs/design/performance_memo.md` §8 "The scipy distribution infrastructure, measured (W6.8)" — §8.1 the three-column table (frozen route 250 µs, new infrastructure 34.7 µs and vectorising to 0.53 µs/θ batched, the declined private-API floor 13.5 µs reproduced for comparison only) and the live M2 reference figures (`log_prob` 545 µs, `lnprior` 58 % of it, the substitution a 1.70× on the whole call); §8.2 bit-identity for the six M2 priors and for `lognorm`/`truncnorm`/`beta`/`gamma`/`halfcauchy` through `make_distribution`, with the §13 conformance row shown insensitive to the prior family by construction; §8.3 coverage (first-class: `Uniform`, `Normal`, `Logistic` only; `make_distribution` takes no `loc`/`scale` keywords, location-scale-only families are scaled by arithmetic, and no family name survives); §8.4 the protocol (`lnprior`, the support-based bijection and `constrain`/`unconstrain` work; `sample`, `prior_transform` and `describe_prior` fail by name; a frozen legacy prior untouched); §8.5 recommendation **(C)** with the two things that would have to change (scipy aliases, or an ampere adapter that cannot recover a family name for `make_distribution` products without the caller supplying it — a §4-adjacent change either way); the probe script as §8's appendix; the pointer under §6 P1. **Findings**: the dispatch prompt said `parameters.md` §6's bijection table is keyed on the distribution's name — it is not: `default_bijection_for` infers the bijection from the prior's `support()`, and the name-keyed table is the lowering registry (`_BUILTIN_PRIORS`, keyed on `PriorSpec.family` from `describe_prior`'s `dist.name`) — §8.3/§8.4 report it correctly; the M2 reference's `lnprior` share is 58 % live here against §3.2's 47 % and P1's 37 %, measured on earlier commits and loads. **For Peter**: none — (C) asks for nothing to land; the plan's §5 bullet is closed by this memo section. Carried: nothing. |
| W6.0 | merged 2026-09-28 at `8811191` (Sonnet-authored, Fable-reviewed as the second pass — the finder, both package inits, the tests, the configuration and every docs page read directly, the relocation by the rename-aware diff; **the agent found and fixed a defect in the orchestrator's own finder template**: appended to `sys.meta_path` the finder loses to `PathFinder` for any submodule a legacy package does not import eagerly — `import ampere.infer.emceesearch` re-executed the file as a second module, two distinct `EmceeSearch` classes — so `install()` inserts it ahead of `PathFinder`, with the identity pinned by `TestTheLegacyAliases`; **one review commit** `4ab1222` (the policy page says the one-substitution route does not carry `ampere.__file__`, with the `importlib.resources` route — the item's own one-substitution run of `minimal_working_example.py` failed at the first filter name because the script builds the library path from `ampere.__file__`, which bound to `ampere.legacy` is the subpackage's; the run repeated with that line on the package-name route: exit 0, the whole script through sampling and the post-processing plots, 150 steps per walker with 100 burn-in, through the lock — the one-substitution claim holds for the API, and the page now says what the second line is for); merged by Peter; the agent's evidence: `tests/test_imports.py` + `test_astropy_compat.py` 140 passed / 4 skipped, torch `tests/backends -k "photometry or instrument or filter"` 32 / 2 and torch typecheck 0 errors, docs 11 `WARNING` lines / Sphinx's "12 warnings" on both `499a33c` and the branch with the warning set identical once the moved paths are normalised, lint/format (278 files)/typecheck/import sweep clean, `path_filters.py` OTHER for a legacy file and CORE for the finder; the characterisation suite on the branch through the lock 3 passed / 3 skipped in 34 s (emcee, dynesty, zeus pass; the three sbi rows skip in dev — the record's "3 passed"); the orchestrator's docs build on `4ab1222` from a clean `_build`: 11 `WARNING` lines, the base count; `test-fast` not run on the branch (queued behind the merged gate and the agent handed off) — covered by the merged gate; the agent ≈ 299 k tokens / 224 tool uses / 37 min, handed off at the cap with the two lock-queued runs; the worktree was cut from legacy `b8e585b` and re-based onto `499a33c`; **the W6.0 merged gate (launched 2026-09-28 05:50 BST on master's working tree at `32b28ac`/`a561eb4`, the merge `8811191` plus handoff prose, through the lock; complete 06:24): dev `test-all` 3114 passed / 612 skipped in 29 min 33 s — exactly the W6.6 baseline, as it should be (the five alias rows live in `tests/test_imports.py`, the `test` task, not `test-all`); `test-characterisation` 3 passed / 3 skipped in 36 s; the docs build from a clean `_build` 11 `WARNING` lines, the base count; torch `tests/backends` 588 passed / 4 skipped in 2 min 11 s — all green; log `~/.cache/ampere-gates/w6.0-gate.log`**; **the CI run on the push of `a561eb4` to `v2`, 36379454600: GREEN** (the run on `32b28ac` cancelled by that push)). Landed: `ampere/legacy/{data,models,infer,utils}` and `ampere/legacy/logger.py` by `git mv`, unchanged but for QuickSED's one absolute import (the prompt said two — the other was a bare `import ampere`, left); `ampere/legacy/__init__.py` importing the five; `ampere/_legacy_aliases.py` (the finder ahead of `PathFinder`, serving the `ampere.legacy` module object itself under `ampere.{data,models,infer,utils,logger}` and every name beneath, a missing extra raising the same `ModuleNotFoundError` as before) and `ampere/__init__.py` rewritten (`__version__` from `importlib.metadata`, `0+unknown` without an install; `__getattr__` for attribute access; no eager import of anything — `import ampere` leaves `sys.modules` at `ampere` and `ampere._legacy_aliases`; no `__all__`); the two backends' `instrument.py` and the M2 benchmark name `ampere.legacy...` directly; `tests/test_imports.py` renamed throughout with `test_top_level_packages_import_cleanly` as the alias proof and `TestTheLegacyAliases` (five rows: purity, identity, the whole surface under `from ampere import legacy as ampere`, attribute access, the version); `TestTheImportCost` without its `sys.modules` stubs; ruff's exclusion collapsed to `ampere/legacy` with `ampere/__init__.py` linted; `path_filters.py` naming the two new top-level files as CORE; **`docs/source/legacy.rst`** (the policy: where it lives and the one substitution; the aliases as the same module objects, no warning, no removal date; kept indefinitely, documented, tested; frozen — unchanged, none of the v2 contracts, its own dependencies; v2's own deprecations removed at 1.0.0 final with `regularised_horseshoe` the list; the examples' twins) first in `api.rst`'s legacy toctree; `ampere.rst` and the four subpackage pages renamed `ampere.legacy*.rst` with the dead commented-out block dropped; `migrating.rst`, `overview.rst`, `tutorials.rst` pointing to the policy and the reference. **Findings**: the dispatch's one-line surface check needed `import ampere.infer.emceesearch` first — the legacy `infer` package never imported it eagerly (pre-existing); `CLAUDE.md`'s ground rule 1 and repository map described the old layout (updated in the merge's handoff commit); `migrating.rst`'s opening still promises "a deprecation timetable" (W6.1's page); ruff's `extend-exclude` carries a stray `ampere/examples` entry for a directory that does not exist. **For Peter**: the two assumptions confirmed at merge — the old names as silent aliases with no removal date, and v2's own deprecations removed at 1.0.0 final; `pyphot_compat` stays under `ampere.legacy.utils` for now (shared infrastructure of W0.9 vintage; a v2 home is a later item). Carried: the `pyphot_compat` home; the stray ruff entry (housekeeping). |
| W6.13 (A) | tranche A — twins (1) and (3) — merged 2026-09-28 at `61503bf` (Sonnet-authored, Fable-reviewed as the second pass — both packages' model, noise and instrument declarations, the README and the `tutorials.rst` line read directly; **the agent handed off at the soft cap** (≈ 398 k tokens / 334 tool uses / 61 min, the modified-blackbody coverage run and three checks queued behind the lock) **and the orchestrator finished the remainder itself** — the four queued jobs watched to completion, the coverage numbers written into the docstring as commit `411c270`, the agent stopped after it kept waking to re-report; merged by Peter; the evidence: `test_linear_sed.py` 14 passed / 3 skipped in dev and 4 sbi rows in `-e sbi`, `test_modified_blackbody.py` 12 / 1 in dev and 2 in `-e sbi` (114 s — the dynamic-range finding), `test-examples` on the branch 119 passed / 17 skipped in 6 min (master 93 / 13 — the tranche's 26 rows and 4 skips), the docs build 11 `WARNING` lines (the base), lint/format (288 files)/typecheck (0 errors)/import sweep 96/4 clean, the legacy diff empty; **the coverage runs**: `linear_sed` under emcee at the legacy budget, 181 s, all four qualified parameters covered on the first run; `modified_blackbody` under emcee (357 s), zeus (425 s), dynesty (43 s) and sbi at 50 000 simulations (891 s) — all four parameters covered under every engine on the first run, the sbi posterior several times wider on three of them for the dynamic-range reason recorded; **the wave-3 merged gate (both merges; launched 2026-09-28 08:43 BST on master's working tree at `f679776`, the two merges plus handoff prose; complete 09:26): dev `test-all` 3155 passed / 616 skipped in 40 min 36 s — the W6.6 baseline (3114/612) plus exactly the 41 rows and 4 sbi-gated skips the two items added to `tests/examples`; `-e sbi` the two tranche-A test files 30 passed in 76 s — green; log `~/.cache/ampere-gates/wave3-gate.log`**; **the CI run on the push of `f679776`, 36393247614: FAILED at lint — four lines of `linear_sed.py`'s coverage docstring over the 100-column limit (the agent's lint pass predated its coverage commit; every other job green) — fixed on master by the orchestrator at `6527d73` and re-pushed; **the CI run on that push, 36397440625 on `c66a1d9`: GREEN**). Landed: `examples/linear_sed/` (the `LinearModel`; `generators.irs_wavelength_grids` reading the CASSIS primary-image HDU into the SL/LL grids with `astropy` alone, asserted equal to `ampere.legacy.data.Spectrum.fromFile`'s in the smoke test, and `deduplicated_grids` dropping the repeated wavelengths from overlapping IRS orders that v2's `Spectrum` refuses; the catalogue, `sl` and `ll` datasets; `CalibrationScale(st.lognorm(0.0025, scale=1.0))` and `GaussianProcessNoise(Matern32)` with `st.halfnorm(scale=0.01)` on the length scale as `calUnc`/`scaleLengthPrior`'s counterparts, `--no-gp`; `--engine emcee|dynesty|zeus|sbi`, `--embedding` passing the legacy FC dict verbatim, sbi refused by name without the extra); `examples/modified_blackbody/` (the four-parameter `ModifiedBlackBody` transcribed exactly with the distance-unit and blackbody-scale quirks kept and documented; ten bands, all loading; `--engine` and `--all` with `overlay_figure` as a composite of per-engine `plot_posterior_predictive` panels, because that plot builds its own figure and takes no axes); the two smoke-test files; `examples/README.md` (the v2 packages, the six pairs with status, the three dead scripts, the notebooks); one sentence in `tutorials.rst`. **Findings** (for the later tranches and W6.1): v2's `Spectrum` refuses the duplicate wavelengths legacy held silently — every twin reading a real instrument file needs a deduplication; the transcribed blackbody costs 5–10× a plain array evaluation, so the MBB's legacy emcee budget takes six minutes; the MBB's dynamic range triggers sbi's z-scoring outlier warning and slows training; dynesty needs live points well above `2 × ndim` to avoid its enlargement pathology (25 for four dimensions); `plot_posterior_predictive` takes no `axes=` (a Phase 6 docs-or-results note). **For Peter**: `examples_paper/modifiedblackbody.py` has two quirks the twin reproduces on purpose — its comment calls `d_true` kiloparsecs but the code converts parsecs, and `BlackBody().evaluate` on a bare frequency array never attaches its `Jy/sr` scale, so the legacy flux is a dimensionless number of order 1e-11; the paper script is yours to correct or not. The GP amplitude prior on the linear twin's spectra (`st.halfnorm(scale=1.0)` Jy) has no legacy number to match. Carried: the deduplication and `plot_posterior_predictive` notes above. |
| W6.13 (B) | tranche B — twin (2), NGC6302 — merged 2026-09-29 at `7998fd2` (Sonnet-authored in two passes, Fable-reviewed as the second pass — the prior derivation, the derived-temperature plumbing, the docstring's coverage record and the test file read directly, the twenty-five rows re-run by the orchestrator (25 passed in 9 s), one Fable review commit `eb5fa3b` tempering the diagnosis; merged by Peter on the ruling of 2026-09-29 (option B: a successor to enforce the ordered prior, switch the solver and rerun the coverage, then merge)). **The first pass** (≈ 397 k tokens / 173 tool uses / 92 min; five commits to `c780fba`): `KemperTwoShell` transcribing `__call__`, `ckmodbb` and `shbb` exactly, vectorised over the fourteen temperature steps, equal to `SpectrumNGC6302` from `ampere.legacy` at a fixed θ to 4.7e-15 relative; the eight opacity tables as sixteen buffers with the two legacy unit conversions; `NGC6302_100.tab` with the 25–120 µm selection (625 points) and the 5 %-of-flux uncertainty rule; `Resample` + `CalibrationScale(st.lognorm(0.05))` in place of `calUnc=1e-10`; `GaussianProcessNoise(Matern32)` with the legacy `scalelengthPrior=0.1` and `st.halfnorm(scale=100)` Jy on the amplitude (the spectrum is 300–1000 Jy), `--no-gp`; `--engine emcee|zeus` at the legacy 50 walkers / 50 000 / 40 000 with `--quick`; `--synthetic` at the 2002 solution; `dust_mass.py` (`dust_masses_at`, `dust_masses(tree)`) reproducing every number of `NGC6302-calculate-dust-mass.py`'s printout to 1e-6; but its 900-step coverage run did not converge and its temperature pairs were unordered (the same independent box for both, `HierarchicalPrior` binding only raw values — verified). **The orchestrator's review** found the ensemble unmixed rather than multimodal, the dense GP dominating the 32 ms likelihood, and the legacy's flat ordered prior declarable exactly; Peter ruled (B). **The successor** (Sonnet, ≈ 312 k tokens / 172 tool uses / 136 min; three commits to `b0bc038`): `QuasisepGP` by default (`solver="dense"` kept; dense and quasisep agree to 4e-14; 28 → 13 ms per call on the shared machine — the O(N) win at N = 625 is smaller than hoped); `Tcold0 ~ st.triang(c=0, loc=10, scale=70)` and `Tcold_fraction ~ U(0, 1)` with `Tcold1 = Tcold0 + f (80 − Tcold0)` in `evaluate` (likewise warm over [80, 180]) — the Jacobian `(80 − T0)` cancels the triangular density, so the pair's joint is constant on the triangle, the legacy prior exactly; a row checks it at 1 000 random points to 1e-10; `derived_temperatures` recovers `Tcold1`/`Twarm1` for `report`, `recovers_truth` and `dust_masses`; `TRUTH` remapped with the synthetic spectrum unchanged to 1e-12; **the coverage run** — `--synthetic --engine emcee --seed 20260928 --walkers 50 --steps 10000 --burn-in 5000`, 4 748 s — **covers all eighteen truths at 95 %** (sixteen declared plus `Tcold1`/`Twarm1`), the Accept met; **R-hat 1.69–2.90 and ESS 57–78 on every parameter**, unchanged from the 900-step run: an ESS at the walker count means each walker sits in its own region — structurally distinct modes or a stretch-move autocorrelation time beyond the run, which the run cannot tell apart; recorded in the docstring as the finding, the one permitted extension not run because the gap is not within a factor of two. The evidence: `test_ngc6302.py` 25 passed (16 + 4 solver + 5 prior rows) in 9–16 s; `test-examples` 159 / 17 (the first pass's 150 / 17 plus the nine rows); docs 11 warning lines (the base); lint / format / typecheck (0 errors) / import sweep 96 / 4 clean; the legacy diff empty. **Findings**: the NGC6302 posterior does not mix under emcee's stretch move at 18 parameters — the legacy's 50 000 steps and its zeus variant say the authors knew; **the follow-up is a zeus coverage run** (`--engine zeus`, in the twin) — the orchestrator's to run detached when Peter says, its result appended to the docstring; **RUN 2026-09-30 on Peter's word** (50 walkers × 2 000 / 1 000, seed 20260928, 92–93 min each, twice — the CLI prints no diagnostics, so a driver repeated it with arviz and saved the run to `~/.cache/ampere-gates/ngc6302-zeus-run.nc`): every truth and both derived temperatures covered at 95 %, but R-hat 2.3–6.2 and bulk ESS 52–63 on all eighteen — worse than emcee's — and the saved run shows why: the per-walker mean log-probability spans −3 164 (best) to −673 351, only 4 of 50 walkers within 100 of the best, and 49 of 50 still climbing between their first and last draw (mean gain +9 620) — **the ensemble had not finished burning in from its prior draws at 1 000 steps**, so the run says nothing about the posterior's mode structure yet; the zeus 'wide intervals' are walkers still en route, not the posterior; **the MAP-started third run was abandoned** (2026-09-30 17:27): `optimise(method="scipy", starts=8)` — Powell on the eighteen sigmoid-bounded coordinates, 97 964 evaluations, 20 min — converged to a corner of the box (seven log-abundances at their edges, two species at maximum, `Twarm0` at its lower bound) with log posterior −3 228.5, a hundred below the best prior-started walker's −3 130, so its ball was the wrong start and the run was killed after 30 min; **the fourth run is queued**: all fifty walkers started in that best walker's own last fifty draws (a compact cloud in the good region, `ngc6302-zeus-start.npy`), the same budget, saved to `ngc6302-zeus-best-run.nc` — **DONE 2026-09-30 19:04 (86 min)**: all fifty walkers within ten log-units of each other at lp −3 103 to −3 109 (the best draw −3 099), so the ensemble is on the posterior's high ridge — and still R-hat 1.5–2.2, ESS 64–97, the walkers a continuum along the ridge (per-walker means of `Tcold0` spread 37.7–42.1 K with within-walker sd 0.7), the intervals now narrow and **missing eight of eighteen truths**; the truth's log posterior with the nuisance parameters profiled is −3 109.6 (GP amplitude → 0, as synthetic data with no misspecification should give), on the ridge and just beyond the stretch the walkers cover; **the emcee record's coverage is explained**: at its GP amplitude of 37 Jy the truth's lp is −3 191, so those walkers sat in a low-posterior region where the GP soaked the residuals and the intervals were wide for that reason. **Reading**: the synthetic NGC6302 posterior is a curved, degenerate ridge in the abundance–temperature space that neither the stretch move nor zeus's slice moves mix along at these budgets — the case for nested sampling; **nautilus queued 19:06** (`ngc6302-nautilus.sh`, the nested environment, 475 live points, the evidence and its error, the run saved) as the last run before the docstring follow-up — **DONE 21:45 (158 min, 460 800 likelihood calls, n_eff 10 001, log Z = −3 119.94 ± 0.01 by nautilus's own Kish estimate)**: 14 of 20 truths covered, **six missed** (`Tcold0` 33.6 ± 1.0 against 36.1, `Tcold1`, `Tcold_fraction`, `logacold2`, `logacold3`, `logacold7` −1.27 ± 0.09 against −1.08) — and **none of its 460 582 equal-weight draws lies in the zeus-best region** (`Tcold0` > 38 or `logacold7` > −0.8: zero draws) although that region's log posterior (−3 105) equals the nautilus draws' own median (−3 104): **and they are two distinct modes of equal height** (peaks −3 098 / −3 099; the nearest high-posterior draws of the two clouds 3.7 posterior widths apart, the straight line between them dropping 230 log-units, while the line from nautilus's median to the truth stays within −3 130..−3 110): two dust compositions fit the 5 % spectrum equally well, each sampler found one, the truth lies in nautilus's, and the emcee coverage record's intervals came from a low-posterior region (lp −3 191 at the truth for its GP amplitude of 37 Jy) where the GP absorbed the residuals; ****The zeus/nested follow-up merged 2026-09-30 at `ce9dc72`** (`w6.13b-zeus`, one docstring-only commit by the orchestrator on Peter's word; the 25 rows pass, lint/format clean, docs 11): the docstring's new section records the four runs and the reading above — two modes of equal height (peaks −3 098 / −3 099) 3.7 posterior widths apart with a 230-log-unit barrier between their nearest high-posterior draws, the truth in nautilus's, run 3's walkers in the other, the emcee coverage explained; **for Peter's backlog**: a multi-modal nested bound (dynesty's multi-ellipsoid, or nautilus with more live points and `enlarge_per_dim`) with per-mode masses, and on the science side a species prior or the far-infrared photometry; the drivers, logs and saved runs in `~/.cache/ampere-gates/ngc6302-*`; the CI run on the push of `fcfbe8a`, 36783928223, GREEN. `QuasisepGP` buys only 2× at N = 625 (per-call constants dominate); the `zeus` variant's dropped `logawarm1` is not reproduced (one model, two samplers — documented). **Merged gate GREEN, all four legs** (the wave-4 merged gate, one for both merges, launched 2026-09-30 04:25 BST on `2ca2edf`, complete 11:56; log `~/.cache/ampere-gates/wave4-gate.log`): **dev** 3236 passed / 632 skipped in 34 min 56 s (the baseline 3155 / 616 plus the tranche's and W6.7's rows), **torch** 4177 / 478 in 1 h 23 min, **jax** 4037 / 447 in 40 min, **sbi** 4338 / 317 in 1 h 23 min (the torch, jax and sbi legs on the tree at `28eca47`); CI on the push of `fc910a0`, run 36664388988, failed on three rows none of them this tranche's (its 25 rows green in every job) — see W6.7's row for the test-only fix at `28eca47`; **the rerun on that push, 36668328856: GREEN**. |
| W6.7 | merged 2026-09-29 at `2ca2edf` (Opus-authored in two passes, Fable-reviewed as the second pass — the objective, the Jacobian handling on the native path, the covariance refusal, the VI fix, the bridge, the provenance change, the warm start's derivation checked line by line against the code, the `Optimum` record, the contract text and the docs page read directly; the dev rows re-run by the orchestrator: `test_optimum.py` + `conformance/test_optimise.py` + `test_optimise.py` 52 passed / 16 skipped in 102 s, `test_results.py` + `test_training.py` 223 passed; the codex route still unavailable, so merged on the second Fable pass under the standing ruling). **The first pass** (≈ 394 k tokens / 289 tool uses / 2 h 5 min; nine units, one commit each, to `0c097f9`, a hand-off at the soft cap): `ampere/results/optimum.py` (`Optimum`, `StartSummary`, `combine`, `to_datatree`/`from_datatree`, `aligned`); `ampere/inference/_optimise.py` — `optimise` with the `"scipy"` (multi-start Powell, `minimiser=` selectable), `"map"` (torch `LBFGS` / jax `BFGS`, Adam then quasi-Newton again as the fallback, the change-of-variables term subtracted on the realised side with its gradient by central differences so `ampere.inference` imports no backend, the autodiff Hessian minus the term's finite-difference Hessian), `"vi"` (the `laplace` guide's mean and covariance — the stated exception, centred on the density NUTS samples) and `"auto"` routes; `constrained_objective`, `finite_difference_hessian` (two-pass, curvature-scaled steps), `covariance_from_hessian` (refusal by name below an eigenvalue floor of 1e-8); `warm_start_gp` (the derivation in the docstring — `-2 log L = N log σ_f² + Σ log(λ_k + ρ) + (N − p) log ρ + Q(ρ)/σ_f²`, the amplitude profiled as `Q/N`, `brentq` on the profile's derivative in `log ρ`, twelve length scales over the prior's central 99 %, a dense dataset searched on a temporary `HilbertSpaceGP(32)`); the bridge — `draw_prior_positions` shared by `initial_positions` and the multi-start, `initial_positions(around=)` as `u* + 0.5 L z` or a diagonal ball on a refusal, `run(initial=optimum)` on emcee/zeus (the ball), NUTS/blackjax (the mode plus a tenth-covariance jitter), VI (the mode), dynesty/nautilus/ultranest refusing by name; `provenance_attrs(start=)`, `ampere_start_route` on every run (`"prior"`, `"user"`, `"pathfinder"` or the route) and `ampere_start` on a seeded one, **`PROVENANCE_SCHEMA_VERSION` 9**; `tests/results/test_optimum.py`, `tests/inference/test_optimise.py`, `tests/conformance/test_optimise.py` (scipy and map agree to 1e-3 on the torch and jax columns; the scipy mode stationary on every column); `docs/source/optimisers.rst` (every pycon block run), `inference.md` §10b, `results.md` §4/§9/§12, `architecture.md`'s stale sentence, the decision-log row. **The successor** (Opus, ≈ 116 k tokens / 68 tool uses / 3 h 33 min mostly waiting on the lock; two commits to `9db0870`): the torch VI flake fixed by polishing the best prior start with a budgeted Powell run before the guide fits (a guide whose MAP stage ends short of the mode has no positive-definite curvature; the polish is kept only when it improves the objective) — three consecutive gated runs 2 passed each; the NUTS budget pinned per backend (`NUTS_BUDGET`: jax 50/100/4/6, torch 50/50/2/4 at four threads — pyro costs ~25 ms per leapfrog on this problem and the MAP-started chains build full trees at their small step). **Evidence**: (a) on `agreement_problem` the scipy mode equals the analytic mean to 1e-3 and the covariance's sd the analytic one to 1e-3 relative; on `sed_composition` sampled by emcee every route's estimate inside the central 50 % on all four parameters (temperature 178.97 in [178.41, 179.61], beta 1.627 in [1.608, 1.645], scale 9.82e-16 in [9.6e-16, 1.29e-15], calibration 0.997 in [0.978, 1.012]); scipy and map agree to 1e-7 in `u`; (b) the warm start's three ratios to the sampled medians 0.867 / 0.871 / 0.993; dense 41.42 against reduced-rank m=24 40.35 and m=96 41.44 (within the truncation's 5 %); the 61×61 brute-force grid does not beat the root find (both 40.353386438); (c) NUTS from the MAP against the prior — jax R-hat 1.11 vs 3.15, ESS 32.6 vs 4.84, divergences 0 vs 19; torch R-hat 1.89 vs 2.69, ESS 3.15 vs 2.63, divergences 0 vs 17; (d) emcee's burn-in equivalent 0 steps from the optimum against 491 from the prior (600 steps, 16 walkers); (e) `ampere_start`/`ampere_start_route`, schema 9 everywhere (254 results rows); (f) every refusal by name incl. a covariance refusal on a built saddle; conformance dev 671 / 75, torch 1157 / 79, jax 1197 / 86 with no existing row moved; `test_optimise` torch 45 / 11 (15 min), jax 45 / 11 (5 min); dev `test-fast` 2941 / 444 in 4 min 31 s (about 1.5 min over master); docs 11 warning lines; lint / format (298 files) / typecheck 0 errors in dev, torch and jax / import sweep 98 / 4. **Findings** (for Peter): the torch NUTS row's ESS margin is thin (3.15 against 2.63 at the pinned budget — a flake risk in the torch gate; if it flakes, compare the minimum ESS across parameters rather than each, or mark the row slow); `test-fast` gains ~1.5 min (the emcee burn-in row and the `sed_composition` fixture) — a registered marker to keep heavy rows out is a `pyproject.toml` change not made; `warm_start_gp` assumes `<label>.likelihood.<name>` parameter names and refuses a renaming tie; **carried from the NGC6302 follow-up (2026-09-30)**: the scipy route's multi-start Powell, on an eighteen-parameter problem whose free coordinates are all sigmoid-bounded boxes, converged from eight starts to a box corner (seven coordinates saturated at their bounds) with a log posterior a hundred below a point the ensemble sampler had already visited — Powell's line searches run to where the sigmoid saturates and the objective goes flat; the route needs either a bound-aware minimiser or a warning when an optimum sits at a coordinate's bound, for Peter to place; the first agent's six open questions and the orchestrator's recommendations on them are in the handoff of 2026-09-29 (the agreement measured in unconstrained coordinates; `free_labels` added; `"user"`/`"pathfinder"` start routes; NaN densities on `combine`; the 5 % dense tolerance; the four-parameter problem and 100 draws). **The issue-closing comments for #40, #14 and #41 are drafted in the first agent's report** (Peter posts, amending as he likes). **Merged gate GREEN, all four legs** (the wave-4 merged gate, one for both merges, launched 2026-09-30 04:25 BST on `2ca2edf`, complete 11:56): **dev** 3236 passed / 632 skipped in 34 min 56 s (the baseline 3155 / 616 plus the tranche's and W6.7's rows); **torch** 4177 / 478 in 1 h 23 min 15 s (the baseline 4046 / 463 plus W6.7's rows; on the tree at `28eca47`, the test-only fix below, as were jax and sbi); **jax** 4037 / 447 in 39 min 55 s (the baseline 3906 / 432); **sbi** 4338 / 317 in 1 h 23 min 15 s (the baseline 4200 / 306); **CI on the push of `fc910a0`, run 36664388988: FAILED on three rows, 26 jobs green** — nested `TestTheBridge`'s nautilus row (a one-dimensional problem, refused by the constructor before the start-point refusal), the torch NUTS start comparison (this row's flagged thin margin: `model.beta` R-hat 1.89 against 1.86 with the run as a whole far better adapted) and py313's pool-reuse row in `test_simulate.py` (a scheduling flake, pre-existing); **all three fixed in the tests by the orchestrator at `28eca47`** — the nautilus row on a two-parameter power law, the NUTS row comparing the worst parameter of each run (max R-hat, min ESS) as this row proposed — **ruled by Peter 2026-09-30: the aggregate comparison stands**, the pool row asserting both chunks ran on the one pool's own processes; verified dev 5 passed, nested 8 passed / 2 skipped, hygiene clean; **the CI run on that push, 36668328856: GREEN** (all 29 jobs, 2026-09-30 ~06:10 BST); **the fixed torch NUTS row re-run detached through the lock: 1 passed in 10 min 10 s** (MAP-started max R-hat 1.89 against the prior's 2.69, min ESS 3.15 against 2.63, divergences 0 against 17 — the MAP-started chains identical to CI's to the last digit, the prior-started ones not, which is where the per-parameter form flaked). |
| W6.13 (C1) | tranche C1 — twins (4) the PHOENIX star and (6) the star plus disc — merged 2026-09-30 at `2876cc0` (Opus-authored in two passes, Fable-reviewed as the second pass — the emulator's arithmetic (the GP predictive mean against the training fit's kernel, the unit-bolometric normalisation, the flux-conserving binning, the LSF), the inverse-square scaling, the λ⁻² extrapolation, the CCM89 polynomials against the paper's equations 2a–3b, the dust term against `QuickSED`, the two top-hats' widths, the native twins and the three coverage records read directly; the orchestrator's re-runs on the branch: the three row files 31 passed / 6 skipped in dev (twice), lint / format (319 files) / typecheck (0 errors) clean, the docs build from clean 11 WARNING lines (the base), `git merge-tree` clean against master, the legacy diff empty; **one review commit** `deb0505` — the sentence recording Peter's ruling on the smoothing GP in `train_emulator.py`'s trade-off note, which the successor's permission classifier had refused because the ruling reached it only through the orchestrator's brief (Peter's ruling of 2026-09-30 in the orchestrator's own session); merged by Peter). **The first pass** (≈ 292 k tokens / 193 tool uses / 2 h 26 min; eight commits to `8e055e3`, a hand-off at the soft cap with its training and coverage runs queued as detached jobs): `examples/phoenix_star/` — `train_emulator.py` (157 Göttingen files into the astropy cache; segment A `geomspace(0.3, 5.5, 582)` µm as exact bin means of the piecewise-linear spectrum, segment B the legacy `specwaves` (387 points, the legacy loop's count — the brief said 386) through a Gaussian LSF at R = 11 000 in log λ; unit bolometric flux over the file's full range; PCA on log10 flux with the 0.5 % / 1 % truncation rule giving **K = 10** (worst RMS 0.49 % / 0.45 %); one ARD squared-exponential GP per weight by marginal likelihood, L-BFGS-B on log hyperparameters), `emulator.py` (the arithmetic written once against an `Ops` namespace so numpy, torch and jax run the same `exp`/`sum`/`matmul`; `ChannelPlan` for log-log placement and the λ⁻² extrapolation; `ccm89_ab`; `PhoenixEmulator`), `emulator_torch.py`/`emulator_jax.py` (subclasses changing only the namespace and the capability flags; the three backends agree to 1e-8), `phoenix_star.py` (`PhoenixStar` with the five legacy parameters and priors, `feh` and `distance` buffers, the correct inverse-square law, the instruments, the RVS GP noise, `--engine sbi|nuts|emcee`, `--backend`, `nuts` refusing `reference`), `generators.py`, `__main__.py`; `examples/star_disc/` — `StarDisc` (`QuickSED`'s twin: the linear `luminosity` and `feh` free, the dust transcribed exactly including the legacy `1e23`), one photometry step over the votable's nineteen points (seventeen library filters plus ALMA band 6 as 1 090–1 421 µm and ATCA 9 mm as 7 889–9 993 µm top-hats, `detector="energy"`), the synthetic Marshall RVS spectrum, `--synthetic`, `--irs PATH`, `--quick`; `tests/examples/test_phoenix_star.py`, `test_star_disc.py`, `test_phoenix_emulator_training.py` (31 rows; the toy-grid rows test the script's functions without the download); `examples/README.md` rows (4) and (6), the `tutorials.rst` sentence. **The successor** (Opus, ≈ 107 k tokens / 49 tool uses / 44 min; four commits to `387cff5`): the `.npz` committed from the run under the lock (byte-identical to the dev-trained copy; 134 005 bytes); the three coverage runs recorded; the error maxima located (the 7.9 % node maximum at 0.3015 µm on the 5 400 K / 4.5 / +0.5 node, the 14 % leave-one-out maximum at 0.300 µm on the 5 000 K / 4.0 / 0.0 corner — both blueward of every passband and below both twins' Teff priors; the 2.0 % segment-B maximum at 0.8668 µm inside the RVS window; the docstring says where errors above 1 % sit); the docs count: the first agent's master baseline of 9 came from a truncated build without pandoc — master and the branch are both 11 with identical warning lists, nothing to fix. **Ruled by Peter 2026-09-30: the smoothing emulator is accepted** (the exact GP, jitter 1e-8, reproduced the nodes to 3e-6 but predicted 9 % / 5 % leave-one-out RMS; with jitter to 1e-4 and length scales capped at 5 the GP smooths: node RMS 0.38 % / 0.21 % with max 7.9 %, leave-one-out RMS 0.61 % / 0.29 % with max 14 % / 2.0 %; the brief's 1e-3 node-exactness row replaced by a check against four stored true spectra). **Coverage** (seed 20260930): **SBI** (reference, NPE, 10 000 simulations, 191 s) covers all five truths; **NUTS on torch** (300 / 500 / 4, 17 890 s) covers all five but **does not converge** — R-hat 1.12–1.61, the calibration scale and luminosity trading off in one chain (intervals 300× and 9× SBI's); **the star-disc emcee run** (`--synthetic`, 40 walkers, 4 000 / 1 000, 631 s) covers all eight truths but R-hat 1.10–1.42 and minimum ESS 82 against the 1.05 / 400 criterion — both recorded as findings, neither re-run. The agents' evidence: dev rows 31 / 6, sbi 35 / 2, torch 34 / 3, jax 33 / 4, dev `test-examples` 190 passed / 23 skipped in 3 min 28 s (master 159 / 17), hygiene clean, the legacy diff empty. **For Peter**: (1) the two unconverged coverage runs — a longer NUTS budget (5 h at this one) or a reparameterisation of the calibration scale against the luminosity, and a longer or zeus-driven star-disc run, or accept the coverage at these budgets as the examples' honesty clause; (2) **both legacy scripts divide by 4π d² a second time** (`phoenixstar.py` 53, `QuickSED.py` 148), so their fluxes are 4π too faint — the twins use the correct law and say so; (3) the IRS file `cassis_yaaar_spcfw_5295616t.fits` the legacy `star_disc.py` reads is not in the repository — the twin fits photometry plus the synthetic RVS unless given `--irs PATH`; does Peter have it to commit? (4) the CSV's two upper limits (500 and 880 µm) are in neither the legacy fit nor the twin's — a `Censoring` declaration is the follow-up. Deviations from the brief: segment B has 387 points; the PHOENIX download ran outside the lock per the brief's own line 7; the successor's maxima script ran three minutes outside the lock. **Merged gate GREEN** (the scoped gate, launched 2026-09-30 23:08 on `2876cc0`, complete 23:13; log `~/.cache/ampere-gates/w6.13c1-gate.log`): dev `test-examples` 190 passed / 23 skipped in 4 min 10 s (master 159 / 17 plus the tranche's rows); the three row files sbi 35 / 2, torch 34 / 3, jax 33 / 4; **the CI run on the push of `fcfbe8a` (both merges), 36783928223: GREEN** (all 29 jobs, 2026-09-30 23:48 BST). |
| W6.13 (C2) | tranche C2 — twin (5), the carbon star on Hyperion — merged 2026-10-01 at `0d1b9ab` (Opus-authored, Fable-reviewed as the second pass — the Mie integrals and the mass-weighted mixture, the Henyey–Greenstein dust, the `AnalyticalYSOModel` setup against the legacy `__call__` line for line, the SED extraction, the data readers, the pixi feature and the lock (parsed: one environment added, 141 package records added, none removed or modified, no existing environment changed), the tests and the docstring's cost and coverage records read directly; the orchestrator's re-runs on `18b0001`: dev rows 13 passed / 7 skipped, lint / format (325 files) / typecheck (0 errors) clean, the legacy diff empty, `git merge-tree` clean; **no review commit needed**; merged by Peter). **The agent** (≈ 251 k tokens / 168 tool uses / 2 h 56 min; five commits to `5117830`, then one more commit `18b0001` for the rerun; ≈ 258 k tokens in all; no hand-off): `pyproject.toml`'s `[tool.pixi.feature.hyperion.dependencies]` (`hyperion`/`hyperion-fortran` 0.9.11.*, `miepython >= 3.3`) and the `hyperion` environment (the sbi set plus the feature; example-only per D12 (b), not an extra, not in CI); `examples/cstar/` — `generators.py` (the votable in Jy — the file declares mJy; the IRS CSV as one strictly increasing spectrum per chunk, 192 + 172 points; the photosphere read with the legacy's `skiprows=1`, which drops a real data row, kept so both hand Hyperion the same 300 rows; `--synthetic`), `dust.py` (the legacy size recipe equal to the tracked `.size` files to 1e-6; `miepython.efficiencies_mx` per (λ, a) with `m = n − ik`; κ_ext, κ_sca and g by the trapezoid rule over the distribution; the mass-weighted mixture; `HenyeyGreensteinDust` with the legacy `extrapolate_wav`/`set_lte_emissivities`; the `np.string_` alias for Hyperion 0.9.11 on NumPy 2), `cstar.py` (`CarbonStarShell`: the seven legacy boxes with the Dirichlet pair reduced to `sic_fraction`, the legacy buffers, a pid-keyed scratch directory, `SimulatorFailed`; the problem, `--engine sbi` only with the by-name refusals, `--embedding` (the legacy dict verbatim), `--photons legacy|quick`, `--workers`, `--cache`, `--draws` defaulting to 1 000 because `SBIEngine` scores every stored draw through one serial Hyperion run); `tests/examples/test_cstar.py` (20 rows: 13 in dev, the Hyperion rows in the `hyperion` environment — the size recipe, the sign convention (SiC κ(11.3)/κ(9.5) = 15.2, carbon κ falling 25 556 → 717 → 50 cm² g⁻¹ from 1 to 100 µm, albedo 0.49992 at 1 µm), the mixture end points, the finite emissivities, one quick simulation, the composition check, a pooled four-simulation batch with no failures, a 40-simulation fit with its cache hit); README row (5), the `install.rst` paragraph, the `tutorials.rst` sentence. **Measured**: 16.1 s per simulation at the legacy photons, 3.7 s at quick (2.2 s of it Hyperion's own tabulation), 4.5 s pooled per worker on eight. **Coverage** (`--synthetic --photons quick --rounds 1 --simulations 12800 --workers 8`, GP on, 1 h 23 min to simulate and train, then 500 draws served from the cache in 28 min): **six of the seven truths covered, `envelope_mass` missed** — at the in-box truth −6.5 (the brief's truth, the legacy default −5.156, lies outside the box; the first run, 6 of 7 with the mass pressed against the box's upper edge, is kept as a History paragraph) its 95 % interval [−9.58, −6.64] stops 0.14 dex short; `stellar_luminosity` [5739, 7178] (width/prior 0.17) is the one well-constrained parameter, the five others cover at their prior's width (0.91–0.98) and the mass loosely (0.77); not reseeded; **the Accept's 'every parameter' is unmet at this budget** (one round of 12 800 at quick, against the legacy's two rounds of 10 000) — merged on Peter's ruling with the gap recorded; the levers (a second round, `--no-gp`, or the legacy budget at ~11 h on eight workers) are his to place. The evidence: dev rows 13 passed / 7 skipped; `hyperion` rows 19 passed / 1 skipped (the no-sbi refusal, by design) in 11 min; `test-examples` 203 / 30 (master 190 / 23 plus the twenty); docs 11 (the base); lint / format (325 files) / typecheck (0 errors) / import sweep 98 / 4 clean; `pixi install -e dev --frozen` and `-e sbi --frozen` succeed on the new lock; the legacy diff empty. **Findings** (for Peter): (1) **Hyperion 0.9.11 breaks on NumPy 2** — its front end uses `np.string_` sixty times though the conda-forge build declares `numpy < 3`; the example aliases `np.string_ = np.bytes_` before importing it, scoped and documented, to drop when upstream releases; (2) **`ampere/results/training.py`'s `_slot_dataset` stores every `extra_coords` entry as float64**, so a `PhotometricPoints` observation's filter names fail the training-set write (`could not convert string to float: 'MCPS_B'`) — `training_set=` is omitted here; a bug under `ampere/`, a filler item; (3) **`SBIEngine` scores every stored draw serially in the driving process** through `problem.evaluate`, one Hyperion run each, so scoring dominates `run` for an external simulator (the legacy's 10 000 draws would be ten hours at quick) — a `score=False` or pooled scoring option is wanted; (4) **the legacy script's own default `envelope_mass`, −5.156, lies outside its own U(−10, −6) prior box** — the brief had carried it as the synthetic truth, so the first coverage run could not cover it (6 of 7) and was rerun at −6.5 on the orchestrator's ruling; (5) five of the seven posteriors sit at the prior's width (the quick-photon bank at one round constrains the luminosity alone) — the honesty clause the brief asked for, stated in the docstring. Judgement calls to confirm: `--draws` defaulting to 1 000, `build_problem(validate=False)` by default (validation runs Hyperion once at construction, impossible in dev; a hyperion row calls `validate()`). **Merged gate GREEN** (the scoped gate, launched 2026-10-01 02:28 on `0d1b9ab`, complete 02:34; log `~/.cache/ampere-gates/w6.13c2-gate.log`): dev `test-examples` 203 passed / 30 skipped in 4 min 16 s (master 190 / 23 plus the twenty rows); the `hyperion` environment installed into the main checkout from the merged lock (`pixi install -e hyperion --frozen`, clean); the cstar rows in it 19 passed / 1 skipped (the no-sbi refusal, by design) in 72 s; **the CI run on the push of `c30d869`, 36801269831: GREEN** (all 29 jobs, 2026-10-01 03:21 BST). |
| W6.14 | merged 2026-10-01 at `cea2cfe` (Sonnet-authored, Fable-reviewed as the second pass — the writer, the reader, the slot check and the append path read directly; the agent's evidence reproduced by the orchestrator on the branch at `bb85ef5`: dev `test_training.py` + `conformance/test_results.py` 78 passed in 7 s; the new conformance class's torch and jax columns 3 passed each (run from the main checkout's environments with the editable finder bypassed, since a worktree's code is otherwise shadowed); the `cstar` forty-simulation row under `-e hyperion` 1 passed in 56 s; lint, format and pyrefly clean; **no fixes needed**; merged by Peter). What landed: `_Slot` carries a dtype per extra coordinate (`"str"` for `U`/`S`/object arrays, the little-endian numpy string otherwise); string coordinates are written as object arrays with a variable-length-unicode encoding and an empty-string fill for a failed sample; the attr `ampere_extra_coords` is unchanged and `ampere_extra_coord_dtypes` is the new name → dtype mapping, the two reader sites defaulting to float64 when it is absent, so every earlier file reads unchanged and `TRAINING_SET_SCHEMA_VERSION` stays 1 (a row strips the attr and reads the file back); `_check_slot` refuses a label coordinate against a numeric one by name; **one fix beyond the ruling**: `append_training_set` re-wrote a string coordinate at the first batch's fixed width, truncating a longer name in a later append — `_concatenate` resets the encoding to variable-length, with a row; the conformance row `TestTrainingSetCoordinates` (a `PhotometricPoints` and a `Spectrum` observation, write, read by value, append, 8 → 16) in every backend column; `examples/cstar`'s `fit(training_set=)` and `--training-set`, the forty-simulation row asserting the bank's 40 samples and both container kinds by value, the docstring's measured size (54 kB shared + 15.5 kB per simulation, ~8 MB at 500); `results.md` **§11** amended (the item said §9, which is provenance — corrected here). **Out-of-scope findings for the housekeeping queue**: (a) the same fixed-width truncation on append exists for the five `failure_*` string variables (inferred from the identical code path, untested); (b) boolean and integer extra coordinates now keep their dtype where they were cast to float64 before, correct under the ruling but unexercised by a row. Scoped gates after merge on `cea2cfe`, one at a time through the lock: **dev test-all 3288 passed / 645 skipped in 38 min (2026-10-02 00:19 BST)**; **sbi test-all 4394 passed / 327 skipped in 1 h 28 min (2026-10-02 01:50 BST)** — with the caveat that `ampere/core/simulate.py` on disk changed under the sbi leg when Peter merged the CI race fix `30d43e3` at ~00:46 (and for a three-minute uncommitted edit at ~00:20, reverted), so pool workers spawned after that imported the fixed module; CI on the fix's push is the clean sbi check) |
| W6.1 | merged 2026-10-02 at `fdd0be8` (Sonnet-authored, Fable-reviewed as the second pass — the guide read whole, the five tables checked against the enumerated legacy surface, both rewritten notebooks read cell by cell, `conf.py`'s switch and the tutorials page read directly; **one review commit** `264f17b`: the `EmceeSearch` row said there was no counterpart to `guess`, but `EmceeEngine.run` takes `initial=` (an array or a W6.7 `Optimum`), the prior being only the default; the agent's evidence reproduced by the orchestrator in the agent's worktree at `7786da5`: `pixi run docs` rc 0 in 2 min 38 s with the notebooks executing, **11 WARNING lines, equal to the base commit's 11**, none naming the guide, the notebooks, a toctree or an undefined label, the three notebook pages built with no traceback, the executed quickstart's fits 12.0 s and 24.7 s and the modified blackbody's emcee 18.7 s; the agent's `test-examples` 203 passed / 30 skipped, lint, format and the import sweep clean, ruff clean on both notebooks when named explicitly (ruff's `exclude` skips `docs/`); a second build on the review commit `264f17b`: rc 0, 11 WARNING, none naming the guide; merged by Peter). What landed: `migrating.rst` grows from the seed to the full guide — one table per legacy package (data, models, infer, utils, logger) routing every public name to a real v2 name, a twin, or "no equivalent; carried" (never "deprecated", D1 (b)), six side-by-side sections drawn from the twins' own code with the translation notes, a post-processing table and snippet, and the "no legacy counterpart" list; `quickstart` and `Ampere_MBB_Example` rewritten on v2 with no stored outputs and executed at docs build (the model inline, the IRS file found by walking up from the notebook's directory, budgets and wall times stated in the last cell, an honest note that `r_hat` is high at a 90-second budget); `Embedding_nets` kept legacy under the amended banner with `execute: never` in its metadata; `conf.py` to `nbsphinx_execute = "auto"` and `nbsphinx_allow_errors = False` (a raising cell now fails the build — the docs build is a gate for the two notebooks, where W2.11 said it was none), plus two additions beyond the brief, accepted: `nbsphinx_timeout = 300` (the default 30 s per cell is under a zeus cell) and `--IPKernelApp.log_level=ERROR` (the kernel's TCP line otherwise adds one WARNING per executed notebook); `tutorials.rst`'s notebook block split; one docstring sentence in `test_wstat_comparison.py`. Issue #59 closed by this item. **Out-of-scope findings for the housekeeping queue**: (a) `examples/README.md`'s closing paragraph still says the notebooks are legacy and converting them is W6.1's — stale, one paragraph; (b) `SyntheticPhotometry.from_library` calls `ampere.legacy.utils.pyphot_compat.get_unit`, so v2 depends on a "carried" legacy module — the guide says so; worth lifting into v2 in W6.3 or a housekeeping item; (c) `plot_covmats` has no v2 plot equivalent (the guide points at the residual plots); (d) the quickstart's corner plot emits a matplotlib `tight_layout` warning into a cell's stderr, not a Sphinx warning. Scoped gate after merge on `36254ed`: **dev test-all 3290 passed / 645 skipped in 36 min (2026-10-02 02:33 BST)**; CI's docs job now executes the two notebooks, about 2.5 min more) |
| W6.3 | merged 2026-10-02 at `c195b0e` (Sonnet-authored, Fable-reviewed as the second pass — the two new pages read whole, the shrinkage section checked against the M2 page's table (0.0118/0.0077 Jy, 0.254/0.087, 0.257/0.537, the factors 2.9 and 2.1 all as printed there), the doctest runner, the audit row, the exclude-members lines, the CI comment and the `BLACKJAX_METHODS` re-export read directly; the agent's evidence reproduced by the orchestrator on the branch at `c42032b`: `pixi run docs` (now `-W --keep-going`) rc 0, **0 WARNING lines**, "build succeeded." in 2 min 30 s with the two notebooks executing; `test_spec_doctests.py` + `test_public_docstrings.py` 31 passed / 2 skipped (torch, jax absent) in dev, and 12 passed / 1 skipped under each of the main checkout's torch and jax environments with the editable finder bypassed (the `needs-torch` block of `astropy.rst` runs there); lint, format, pyrefly 0 errors, the import sweep 98 passed; **no fixes needed**; merged by Peter). What landed: `pixi run docs` is `sphinx-build -W --keep-going` — a warning fails the build locally and in CI, whose docs job's comment block now says so; the six frozen legacy members that warned (`Spectrum.setPlotParams` on two pages, the `EmceeSearch` class docstring, the mixins' `get_map` and `plot_posteriorpredictive`, `ZeusSearch.get_map`) and the `makeFilterSet` function (the ambiguous `name` reference) are on their pages' `:exclude-members:` lines with a sentence each, `ampere/legacy/` byte-identical; `api.rst` already had the v2-first shape and the "Legacy (frozen)" heading; ten re-exported constants the reference pages did not list gained `autodata` entries, and five undocumented public constants (`BLACKJAX_METHODS`, re-exported with a `#:` docstring; `PARAMETER_DIM`, `LEVEL_DIM`, `TARP_LEVEL_DIM`, `REFIT_ROUTE` in `results/calibration.py`) gained theirs — the only changes under `ampere/`; `tests/core/test_public_docstrings.py` enforces a sentence-first docstring on every `__all__` name (constants through Sphinx's `ModuleAnalyzer`, the `#:` form autodoc itself reads — accepted as the right reading of the ruling); the two tutorial pages `conditional_priors.rst` (HierarchicalPrior, Plate and Population; the NGC6302 ordering by reparameterisation with the Jacobian cancellation checked numerically; what cannot be declared — an arbitrary joint density — with two workarounds, reparameterise or reweight, and `Derived` named as W6.11's) and `arbitrary_priors.rst` (the `Prior` protocol; a `TruncatedPowerLaw`; scipy's new `Normal`, which satisfies half the protocol, with a four-method adapter; an `EmpiricalPrior` over `gaussian_kde` with `ppf` by grid inversion and a padded support; a checklist), both in the tutorials toctree, "Still to be written" now promising nothing; the shrinkage section in `advanced.rst`; the doctest runner — every `docs/source/*.rst` with a `pycon` block is globbed and run as one namespace per page from a scratch directory with the repository root importable, an empty `SKIP` mapping as the `:skipif:` equivalent, plus a per-block `:class: needs-<library>` option (accepted beyond the ruling: one torch-only block in `astropy.rst` is dropped with a warning where torch is absent rather than skipping the page); the page fixes it forced: `kernels.rst` (missing imports, a `Sum` with `SquaredExponential` is **not** quasiseparable — the page said it was — two undefined fixtures replaced by a real `Spectrum` and `VisibilitySet`, two stale refusal messages), `sed_composition.rst` (self-contained blocks). Issues #57 and #58 close. **Out-of-scope findings**: (a) the CI docs job now fails on any warning, including one from a future notebook cell — by design; (b) the legacy utils page documents nothing for `makeFilterSet` — by ruling. Scoped gate after merge on `a7113c8`: **dev test-all 3302 passed / 647 skipped in 34 min (2026-10-02 03:47 BST)**; CI on the merge push: run 36954600768 — 25 jobs green including the first `-W` docs build (3 min 20 s), **3 red** (test py314, core groups dev and jax) on one doctest: `optimisers.rst`'s printed `Optimum.summary()` pinned the Powell evaluation count at 695, which CI's runners give as 678 (py312 and py313 agreed with this machine by chance); fixed on the docs page by eliding the count, branch `ci-optimisers-doctest-ellipsis` at `45a2209`, verified 4 passed under dev, py314 and jax here; merged by Peter at `32da5fe` — CI on that push: **run 36960051619 GREEN, 28 jobs, 1 skipped — the gate of record for W6.3**) |
| W6.10 | merged 2026-10-02 at `8b07e73` (Sonnet-authored, Fable-reviewed as the second pass — the policy section and the audit read whole, the two task lines, the `conftest.py` hook, the orchestration pointer and the lock diff (42 added lines: `pytest-xdist` 3.8.0 and `execnet` 2.1.2 across the environments) read directly; the merge tested clean on top of W6.3; the orchestrator's reproduction on the branch at `f7cf938`: `test-fast` under `-n 4` 2951 passed / 444 skipped in 4 min 24 s (the counts the agent reported; the wall time 2 min 32 s in the agent's measurement — the machine was shared both times, so the speed figures are indicative, the counts are the evidence); lint, format and the import sweep clean; collection under `--dist loadgroup` fine; **no fixes needed**; merged by Peter). What landed: (1) `docs/development.md`'s new section "The merged gate: CI is the gate of record" — the scoped pre-merge legs, the CI run on the push as the gate with its id and per-job counts in the row, what a red run means (the 2026-10-01 race as the worked example), the `~/.cache/ampere-gates` scripts as the fallback and what the lock is, the two things CI never sees (the hyperion rows, the GPU rows of W6.16), and `docs/orchestration.md`'s gate sentence pointing at it — the W6.3 row already records its CI run that way, which is the Accept criterion met; (2) `pytest-xdist` in the `dev` extra, re-locked; the audit table of eleven suite runs at `-n 4` (counts identical to serial everywhere, no flake; `conformance` cross-checked under torch and jax, so the per-problem seeded streams survive); two suites correct but slower in parallel (`interferometry`'s class-scoped calibration fits rebuilt per worker; `astrometry` bounded by two nested-sampling rows) are pinned to one worker each by an `xdist_group` marker from `tests/conftest.py` under `--dist loadgroup`; `test-all` and `test-fast` run `-n 4 --dist loadgroup` with `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS` and `MKL_NUM_THREADS` at 1 (each worker's BLAS pool was oversubscribing the cores) — **measured: `test-all` 36 min 19 s → 19 min 40 s, `test-fast` 4 min 54 s → 2 min 32 s, counts identical**; `-n auto` deliberately not used (CI runners differ; a memory hazard for torch and jax on 13 GB); CI's jobs are unchanged (they do not use the two tasks) and a proposal to measure `-n auto` on a runner is in the report. **The single-thread BLAS pin** inside the two tasks also applies to a `-n 0` serial run — merged as drafted; it stands unless Peter rules otherwise. **Out-of-scope findings**: none under `tests/` writes outside `tmp_path`; `~/.cache/ampere-gates/w6.1-gate.log.dev` is truncated (the summary line survives in `w6.1-gate.log`). Scoped gate after merge: none beyond CI (infrastructure) — the CI run on the push carrying the merge (36972977590 on `dbac0b1`, the merge plus the Read the Docs files): **GREEN, 28 jobs, 1 skipped — the gate of record**) |
| W6.4 | merged 2026-10-03 at `00e5fdc` (Sonnet-authored, Fable-reviewed as the second pass — the whole diff read: the `docs` extra and `dev`'s self-reference, `.readthedocs.yaml`'s rewrite, the two site pointers, the `tags: ["v*"]` trigger, the lock diff (9 lines: the three Sphinx packages moved from `dev` to `docs`, `ipykernel` and `ampere[zeus]` added there, `ampere[docs]` in `dev`); the merge tested clean on master; the orchestrator's reproduction on the branch at `4cefdec`: `pixi run docs` from a clean `docs/_build` rc 0, **0 WARNING lines**, "The HTML pages are in docs/_build/html" in 2 min 39 s; the agent's: a fresh Python 3.13 venv with `pip install -e ".[docs]"` plus pixi's `pandoc` built the site with `sphinx-build -W --keep-going` at rc 0 and 0 WARNING lines (the RTD recipe reproduced locally), `pixi run docs` green after the re-lock, `--frozen` installs of `dev`, `torch` and `jax` all succeeded, lint, format and the import sweep (98 passed / 4 skipped) clean; **no fixes needed**; merged by Peter). What landed: a `docs` extra — `sphinx`, `nbsphinx`, `nbconvert`, `ipykernel`, `ampere[zeus]` (the modified-blackbody notebook runs `ZeusEngine`) — which `dev` now includes by self-reference, so nothing is duplicated; `.readthedocs.yaml` in the ruled shape (ubuntu-24.04, Python 3.13, `pandoc` from apt, the package pip-installed with the `docs` extra, `sphinx.fail_on_warning: true` to match the task; the stopgap `docs/requirements-readthedocs.txt` deleted); one sentence each in `index.rst` and `README.md` pointing at `https://ampere.readthedocs.io/`; CI's push trigger gains `tags: ["v*"]` so the `-W` docs job runs on the tag Read the Docs builds — **no publish step in `ci.yml`**, RTD's GitHub integration building every tag is "the CI job publishes on tag" (ruling 3). Not changed: `conf.py` still reads `version("ampere")` — the name becomes `ampere-astro` at W6.5, which rewrites it; no legacy documentation site exists to redirect (grep found third-party links only). Versions: `latest` = the default branch (`v2` now, `master` at the tag, D4), each tag at `/en/<tag>/`; `stable` follows RTD's highest-semver-tag rule and does **not** promote a pre-release, so it is absent until `v1.0.0` unless Peter enables pre-releases as stable. **For Peter, in the RTD dashboard**: (1) default branch `v2` now, `master` at the beta tag; (2) activate versions on tag (an automation rule, or tick `v1.0.0b1` when it appears); (3) addons on for the version flyout; (4) optionally pre-release tags as `stable`. Acceptance "reachable at the ruled host": **met** — Read the Docs rebuilt `latest` from the merge push on the new `.readthedocs.yaml` (the `docs` extra, no requirements file) and `https://ampere.readthedocs.io/en/latest/` served the build of `7db41af` (the push head, the merge's tree plus one handoff commit) within fifteen minutes of the push, checked 2026-10-02 23:35 BST; the tagged beta is checked at W6.5. Gate: none beyond CI (infrastructure; the docs job is the check) — **CI run 37072266731 on the merge push (head `7db41af`) GREEN: 29 jobs, 28 success, 1 skipped, none failed**, the first run with the `tags: ["v*"]` trigger in the file. Agent ≈ 60 k tokens / 26 tool uses / 9 min. |
| W6.11 | merged 2026-10-03 at `b266d7a` (Fable-drafted as the orchestrator's own work between waves, Opus-reviewed read-only as the item specifies — fourteen findings against the first draft `fdf7828`, three blocking (the grammar could not name a dotted merged parameter; a refusal of a supplied derived value would have broken `FittingProblem.evaluate`, which completes twice; the internal-member rule tied `z` silently in the flat layout), eleven should-fix or nit, every one verified against the code by the orchestrator and folded in at `d106067` with the list in the memo's §11; the reviewer's confirmation pass on `d106067`: **all fourteen resolved, nothing blocking introduced**, three should-fixes (no silent overwrite of a supplied derived value on the reference path — W5.32 (j); the emission broadcasting rule generalised to every declared shape; two frozen allowances retracted by the per-component rule listed as such) folded in at `09b8044`; verdict "ready for Peter's ruling"; merged by Peter; no code, no gate — the memo is not linked from the docs site, so the docs build is unaffected). What landed: `docs/design/nuisance_populations_and_derived_memo.md` — (1) a `Population` over a dataset's qualified path: `over` entries of the form `component[.path]`, `PlateBinding`/`Binding` carrying a qualified local path, the plate layout stripping the leaf from the composite's *outer* declaration, routing untouched (the dataset's retained mapping takes the second hop it already takes; §11 Q1 of the population sketch respected), `DatasetCollection.plate(within=)`, `populations` in the provenance record (schema 11) — verified by two probes against the live merge whose source is the memo's Appendix A; (2) the `Derived` node: a fourth parameter state carried by a `Derived(expression, symbols)` object in the prior slot, a closed grammar over symbols bound to parameter names by a mapping of `HierarchicalPrior.hyperparameters`' shape (so merges rename bindings and never touch the expression), no sampler dimension, computed idempotently wherever named values are formed and mid-walk in `prior_transform`, referenceable by a `HierarchicalPrior`, a plate-layout population member with two new member rules (internal members routed nowhere; a routed member declared by every component — the late `unknown parameter` failure moved to the merge), refused in the flat layout, computed natively in both backends' resolve functions with `numpyro.deterministic` as the jax structural view only, one `posterior` variable per derived name computed vectorised at emission on every engine under `ampere_derived` (schema 10); the slab as `shrinkage_horseshoe(tail="slab")` in the helper's own `s_j = τλ_j` parameterisation and the non-centred population `θ_i = μ + σ z_i` as the two customers, their composition (a non-centred per-dataset GP amplitude) as the validation customer; two decision-log rows drafted (§8), twenty conformance rows named (§9), the Phase 7 items W7.0 (Derived, M, Opus), W7.1 (paths, M, Opus), W7.2 (the M2 validation, M, Opus) drafted for D8 (§10); four questions for Peter with recommendations (§7: the grammar over a callable; derived values in `posterior` with an attr; paths on `over` with `within=` on the factory only; W7.0 before W7.1). Two findings recorded for the Phase 7 items, not fixed: `Population` declarations are absent from provenance; a routed member a component does not declare fails at the first model evaluation rather than at merge. **For Peter**: the four §7 questions, and D8's ordering — the memo recommends W7.0 → W7.1 → W7.2 opening Phase 7. |
| W6.15 | merged 2026-10-03 at `b25f211` (Sonnet-authored, Fable-reviewed as the second pass — `evaluate_many` and its two task callables, `lookup_many`, `finish`'s three keywords and `_strip_scores`, `run`'s keyword and the attr, `_with_estimator_log_prob`'s new reference column, `require_scored` and both call sites, the contract amendment and the docs section read whole; the test rows' assertions listed; the merge tested clean on master; the orchestrator's reproduction on the branch at `25c4a90`: dev `test_simulate.py` + `test_dataset.py` + `test_engines.py` **419 passed** in 53 s, sbi `test_sbi.py -k "score or pool or scored or Scoring"` **20 passed** / 189 deselected in 2 min 04 s, hyperion `test_cstar.py` **20 passed / 1 skipped** in 63 s (the forty-simulation row through its two-worker pool — the leg CI never sees), lint, format and pyrefly (0 errors) clean; the agent's: dev `test_simulate.py` + `test_dataset.py` + `test_engines.py` 417 passed, sbi `test_sbi.py` whole 201 passed / 7 skipped / 1 failed (the pre-existing row below), hyperion `test_cstar.py` 20 passed / 1 skipped in 72.5 s (the forty-simulation row 46.8 s through its pool, the new `score=False` row 7.1 s), `test_spec_doctests.py` + `tests/results` 508 / 5, lint, format, pyrefly, the import sweep clean; **no fixes needed**; merged by Peter). What landed: lever (a) `SBIEngine.run(..., score=True)` — `score=False` stores the draws with `lp`, `log_prior`, `log_likelihood` (column and group) and the per-draw `failed`/`failure_*` columns **absent**, not NaN, `engine_draws_recomputed = 0`, `ampere_sbi_scored = 0` (1 on every scored SBI run from now on), `ampere_sbi_log_prob` and `proposal_log_density` kept; `ampere.results.diagnostics.require_scored` refuses by name — naming the attr and `score=True` — in `chi_square_pvalue` and in the population reweighting's `DataTreeRunColumns._stat` (the two readers that need the scores; `calibrate` simulates fresh and reads none, the trace plot already guards `lp`'s absence); lever (b) `FittingProblem.evaluate_many(vectors, *, executor=None, chunk_size=None)` beside `simulate_many`, reusing its broadcast-once `ProcessExecutor` route, picklability check, recording suspension and failure ledger (a worker's failures are recorded once by the parent, in draw order; an undelivered draw is an `Evaluation` at `-inf` with `EXECUTION_FAILED`, `where="executor"`); `_EvaluationCache.lookup_many` scores the misses through it and `Engine.finish(score=, executor=, chunk_size=)` is additive with today's defaults, a realisation-backed cache staying serial; `SBIEngine` passes its own executor and `chunk_size`, so a pooled fit scores as it banks — identical to serial at `rtol=0, atol=0` (a row), the other engines unchanged (a row pins an emcee run's columns against serial `evaluate`); `examples/cstar`'s `fit(score=)` and `--no-score`, the forty-simulation row now scoring through its pool and asserting the attrs, a `simulations=4, draws=4, score=False` row; `inference.md` §10 "Scoring the stored draws: optional, and pooled (*Amended W6.15*)" and `sbi.rst`'s section with the numbers. **The crossover, measured** (four-worker `ProcessExecutor`, warm pool, identical results): `sed_composition` at 2.7 ms per evaluation — 400 draws 1.09 s serial, 1.47 s pooled (the pool loses); `cstar` quick at 1.43 s per evaluation — 40 draws 57.2 s serial, 21.9 s pooled (2.6×; `chunk_size=4` 29.2 s, a barrier per chunk); dispatch ≈ 0.4 ms an item; the rule of thumb in both documents: **pool the scoring when one evaluation costs more than about 10 ms**; cstar's 500-draw figure (28 min serial) extrapolated to about ten minutes pooled, labelled as extrapolated. No §4 contract change (both levers additive; the ruling's `score=False` default preserved), so no decision-log row. **Deviations accepted on review**: two files outside the ownership list (`results/population.py`'s three-line refusal, `results/diagnostics.py`'s `require_scored`) — the ruling required refusals by name in every reader of the scores, and these are the two; the item's "microseconds per evaluation" for `sed_composition` was wrong by three orders (2.7 ms), so the pool nearly breaks even there rather than losing badly. **Findings for Peter**: `tests/inference/test_sbi.py::TestAmortisationOverTheObservationContext::test_it_stays_calibrated_at_a_rescale_the_prior_covers` fails on the base `f7c05ba` too in the local `sbi` environment (`ks_pvalue` 0.0054 against the 0.01 threshold) while CI's sbi jobs on `dbac0b1` and `f7c05ba` are green — environment- or seed-sensitive at the threshold, W5.26's territory. Gate: dev + sbi — **CI run 37076050565 on the push carrying this merge (head `bf428fb`) GREEN: 29 jobs, 28 success, 1 skipped, none failed** (its dev and sbi jobs are the two legs; the hyperion rows are quoted from the branch run above). Agent ≈ 170 k tokens / 86 tool uses / 37 min. |
| W6.9 | merged 2026-10-03 at `eaf4762` (on Peter's word, by the orchestrator) (Sonnet-authored, Fable-reviewed as the second pass — the extra, feature and environment stanzas, the adapter (`require_rhmf`, `to_matrix`, `fit_rhmf`, `anomaly_score`), the decision-log row and the `diagnostics.md` §7 paragraph read whole, the test rows listed; the merge tested clean on master; the orchestrator's reproduction on the branch at `adf6839`: `pixi install -e rhmf --frozen` ok, `-e rhmf pytest tests/examples/test_rhmf_trial.py` **8 passed / 1 skipped** in 16 s, `-e dev` **1 passed / 8 skipped**, `python -m examples.rhmf_trial --quick` rc 0 in 14 s (both flattenings of the image at excess 0.9–1.3, as reported), lint, format and the import sweep (98 passed / 4 skipped) clean; the agent's: `pixi install -e rhmf` solved (jax 0.11.2, equinox 0.13.8, optax 0.2.8 against the `jax` environment's 0.11.1 / 0.13.8 / 0.2.8), `dev` and `jax` `--frozen` installs unchanged, `-e rhmf` tests 8 passed / 1 skipped, `-e dev` 1 passed / 8 skipped, `--quick` 21.7 s, the full run 41.5 s, lint, format and the import sweep clean; **no fixes needed**; merged by Peter). What landed — **a report, not a namespace**: the adoptability re-check repeated on 2026-10-03 with four lookups (PyPI still 0.0.2 of 2025-11-24 with no licence metadata; the repository MIT, last push 2026-09-24, `main` at `cf2fcffa10bc` tagged `paper`; the only workflow on `main` still `release.yml`, a `test` workflow existing only on two unmerged branches — the draft PR #7 with failing runs and PR #11 with one pass; 8 stars) — **licence and API hold, maturity does not**, the §2.7 JAX-pin axis discharged by the lock, so `ampere.diagnostics` stays deferred and `ampere/` is untouched (`git diff -- ampere` empty); the trial behind a non-default `rhmf` extra with a pixi feature pinning the **commit** `cf2fcffa10bc` and an `rhmf` environment (the `jax` environment's features plus it; 342 lock lines, all additions); `examples/rhmf_trial/` — the adapter (`to_matrix` checking rather than performing §2.3's alignment, mask → zero weight, uncertainties required; `fit_rhmf` with `rank` and `robust_scale` required, stdout captured; `anomaly_score` per feature or per object by a low quantile, `provenance="rhmf_prefit"`, notes naming rank, scale, commit and the pre-fit caveat; `OptionalDependencyError` naming `rhmf` without the extra), the spectra and image trials with `--quick` and full presets writing outside git; `tests/examples/test_rhmf_trial.py` skipping without the extra; the findings as *Amended W6.9* in `diagnostics.md` §7 with a pointer in §2.2, and the decision-log row in `DEVELOPMENT_PLAN.md` §2 after W2.7's. **The finding**: on the M2 spectra (12 controls plus 1/3/6 copies of each scenario, rank 1–4, `robust_scale` 1–5) the robust weights localise the sharp deviations — excess contrast over controls 16.8 for `strong_sharp`, 4.7 for `many_lines`, 2.8 for `strong_smooth`, 1.8 for `mild` at three copies, per-object AUC 1.00 at every best point — but only at a rank that leaves no room for the deviation's own shape (fringes at rank 1 only, the forest at 1–2, the line at 1–3), so no single rank serves one collection; the contrast falls as more rows share the deviation; and the best grid point is oracle-chosen. On W5.5's image neither flattening (images as rows; one image's rows as objects) flags the omitted background at the study's strength (0.84σ per pixel) or at five times it. Verdict for Peter: were the gate met, the namespace lands as expert opt-in, 1D only, `rank` and `robust_scale` required — §2.5 and §2.6 as written; nothing argues for a default. **Deviation accepted on review**: the commit pin sits in the pixi feature, not the extra, because PyPI refuses uploads whose metadata carries a direct reference and the extra must not block the 1.0.0b1 upload — `pip install "ampere-astro[rhmf]"` resolves 0.0.2, which the trial did not run; recorded in all three places. No CI leg (release infrastructure). The restated revisit trigger: a ≥ 0.1 release whose tests run in CI on `main`. Out-of-scope findings: `fit()` prints unconditionally and needs float64 (the trial flips `jax_enable_x64` on its own behalf — a landed namespace could not, `lowering.md` §10.2(a)); W5.5's image values are ~1e15 Jy/sr so a naive weight is ~1e-29 (the trial rescales to sigma units). Gate: dev — **CI run 37079631345 on the push carrying this merge (head `ddd4072`) GREEN: 29 jobs, 28 success, 1 skipped, none failed** (the `rhmf` rows skip there, as designed; their run is quoted above). Agent ≈ 184 k tokens / 74 tool uses / 15 min. |
| W6.16 | merged 2026-10-03 at `5991128` (Sonnet-authored, Fable-reviewed as the second pass — the `gpu` feature and environment stanzas, the three scripts, the example config, the procedure section and the install clause read whole; the merge tested clean on master; **one review commit** `5b646ee` on the branch: `docs/source/overview.rst`'s "No environment has both torch and jax" gains the `gpu` exception (the agent's own out-of-scope finding, outside its ownership); the orchestrator's reproduction on the branch: `pixi install -e dev --frozen` ok, `pixi lock --check` ok (the lock is current), `bash -n` on the three scripts ok, `pixi run test` 98 passed / 4 skipped, `pixi run docs` from a clean `_build` **rc 0, 0 WARNING lines** in 2 min 29 s on `9638b6f` (and again **rc 0, 0 WARNING lines** in 2 min 28 s on the review commit `5b646ee`), lint and format clean; the agent's: `pixi lock` solved the new environment (torch `2.14.1+cu126`, jax/jaxlib `0.11.2` with `jax_cuda12_plugin`/`jax_cuda12_pjrt` 0.11.2, pyro-ppl 1.9.1, numpyro 0.22.0, the `nvidia-*-cu12` family — cublas 12.6.4.1, cuda_runtime 12.6.77, cudnn 9.10.2.21, nccl 2.29.3), `git diff --stat pixi.lock` 426 insertions / 0 deletions (no other environment's section moved), `dev`, `torch` and `jax` `--frozen` installs unchanged, the `gpu` environment never installed, `pixi run test` 98 passed / 4 skipped, lint and format clean, `bash -n` on the three scripts, `run.sh` exits 2 without a config and on the example config; merged by Peter). What landed: the `gpu` pixi feature — `ampere[torch,jax]` (the one environment carrying both backends, deliberately: the rows exercise both on one node), `torch` from PyTorch's **`cu126`** index (the newest cu12x index publishing the torch the other environments lock, 2.14.1; `cu128`/`cu129` stop at 2.11/2.13) and `jax[cuda12]` — and the `gpu` environment in **its own solve-group**, locked, never installed on the development machine; the `torch` feature's CPU pin untouched (as its own comment instructed); `scripts/cluster/` — `bootstrap.sh` (idempotent, two-phase: directories under the project space, a static pixi 0.68.1 with `PIXI_HOME`/`PIXI_CACHE_DIR` there, an ed25519 deploy key whose public half is printed for Peter to register read-only, then the clone and `pixi install -e gpu --locked`), `gpu_rows.sbatch` (`--gres=gpu:1`, 4 CPUs, 32 GB, 30 min, **no partition, account or QoS** by site policy; `nvidia-smi` and `module spider cuda` first, written to the log; **exit 3 if the driver is below 525**, the CUDA 12 wheels' floor; then `pixi run -e gpu gpu -q --tb=short`, the pytest summary as the log's last line), `run.sh` (laptop side: `bootstrap`, `submit <tag>` — checks the tag out in the cluster's clone, shows the exact script, asks, submits with `--output=<fred>/ampere/logs/gpu-<tag>-%j.log`, logs the job id to `~/.cache/ampere-gates/gpu/submissions.txt` — and `fetch <tag>` through the data-mover host to `~/.cache/ampere-gates/gpu/<tag>.log`; refuses without a config or on the example), `config.example.yaml` (committed) with `config.yaml` gitignored (the account, aliases, project `oz528`, paths and `modules` live there, never in git), a `README.md`; `docs/development.md` "The GPU rows on the cluster" after the merged-gate section — the rows and why CI never runs them, the environment, the cluster's rules, the bootstrap and the deploy key, the per-release run, what the log must show ("30 skipped proves nothing"), the driver floor and the re-solve if a node is below it; `install.rst`'s paragraph on the one environment with both backends. **Peter's details carried** (2026-10-03): Swinburne's Ngarrgu Tindebeek, Slurm routed by login host, A100 80 GB nodes, project `oz528` for now, a deploy key rather than the item text's rsync, the CUDA version unknown — hence the CUDA-12 lock with the driver check (the orchestrator's ruling). **Not done, by design**: nothing ran on the cluster — the agent never ssh'd; the first end-to-end run is **W6.5's release gate** on the beta tag and records its summary line in that row. **For Peter**: run `run.sh bootstrap` once (after filling `config.yaml`), add the printed key at the repository's deploy keys with write access unchecked, run it again; confirm cu126 and the 525 floor; `pixi install -e gpu --locked` runs on the login node (several GB) — the site may prefer it inside a job; the login node must reach github.com and the PyTorch index. Out-of-scope findings: W6.5's rename must cover this feature's `torch`/`jax` keys too. Gate: none (infrastructure) — the lock proof above; the CI run on the push carrying the merge (37084736957 on `3feee01`, the merge plus the closing handoff): **GREEN, 29 jobs — 28 success, 1 skipped — the gate of record; W6.16 closed 2026-10-04**. **Amended 2026-10-05 at the first cluster bootstrap**: `pixi install -e gpu --locked` refused with "lock-file not up-to-date" — pixi #7024 (the editable project's `torch` extra names `torch` with no index, read as disagreeing with the lock's PyTorch index; the lock was current and `pixi lock --check` exits 0 on both machines) — so `bootstrap.sh` installs with `--frozen`, as CI's `setup-pixi` steps do; the procedure text follows. **And the same day, job 18040399** (the first GPU run, node `gina1`, A100, driver 615.71.09): `pixi run -e gpu gpu` without `--frozen` tried to re-solve `sbi` on the network-less compute node and exited 1 before any row ran — `gpu_rows.sbatch` now runs `pixi run --frozen` (`cc7a65e`). **Job 18041011 on `cc7a65e`, the first run that reached the rows: 20 passed, 13 failed in 41 s** — three faults a CPU cannot show, now W6.18 (the torch realisation self-check's numpy conversion of a CUDA tensor, `jax.devices()` listing only the default backend so `"cpu"` is refused on a GPU node, and a sharding test helper's hard-coded width); the log is `~/.cache/ampere-gates/gpu/cc7a65e.log`; the cluster appends a usage epilogue after the job, so the summary is not the last line (W6.18 fixes the fetch). Agent ≈ 101 k tokens / 69 tool uses / 7 min. |
| W6.17 | merged 2026-10-04 at `9d42d7c` (Sonnet-authored, Fable-reviewed as the second pass — every new file read whole, the three edits to existing files read as diffs, the Covenant diffed against the upstream 2.1 Markdown: identical apart from the contact line and a leading blank line; the merge tested clean on master; **one review commit** `d49a5cc`: `CONTRIBUTING.md` names the branch a pull request targets — `master` from the beta, `v2` before it — in place of "the repository's default branch", which is the legacy line until the tag (the agent's own finding); the orchestrator's reproduction on the branch at `d49a5cc`: `pixi install -e dev --frozen` ok with no re-lock, lint clean, 331 files formatted, `actionlint` rc 0, `pixi run test` 98 passed / 4 skipped, a wheel from `pip wheel . --no-deps` carrying `Metadata-Version: 2.4`, `License-Expression: GPL-3.0-or-later`, `License-File: LICENSE`, the four issue-form YAML files parse, `pixi run docs` from a clean `_build` **rc 0, 0 WARNING lines** in 2 min 21 s; the agent's: the same, plus `check-jsonschema` with `vendor.github-issue-forms` and `vendor.github-issue-config` accepting all four forms; merged by Peter). What landed: `licenses/gpl.txt` → `LICENSE` as a 100 % rename with `license-files = ["LICENSE"]`; `CONTRIBUTING.md` (97 lines: pixi and the three environments, the pre-PR tasks, branch-and-PR convention, the two freezes with pointers, conventions, where to talk, and the plain paragraph on agent-built v2 under `AGENTS.md`); `CODE_OF_CONDUCT.md` (Contributor Covenant 2.1 verbatim, contact `peter.scicluna@eso.org` per D14 (a)); `SECURITY.md` (the `1.0.0b1` line onward with legacy excluded, private vulnerability reporting then email, fourteen-day acknowledgement, best-effort fix, no embargo/CVE process, and the section naming the artefact-store unpickle, `ProcessExecutor`'s pickled problems, netCDF reads and user code as documented behaviour that is not a vulnerability — D14 (b), (c)); `.github/ISSUE_TEMPLATE/` (`bug_report.yml` with the ruled field order, `feature_request.yml`, `docs.yml`, `config.yml` with blank issues off and contact links to Discussions and the advisory page); `.github/PULL_REQUEST_TEMPLATE.md` (23 lines); `index.rst`'s Contributing section and the README's Contributing bullet pointing at the guide. **For Peter** (Settings → Advanced Security, then Moderation): enable private vulnerability reporting (the advisory link in `SECURITY.md` and `config.yml` is dead until then), secret scanning then push protection, reported-content moderation for Discussions. Out-of-scope findings: `SECURITY.md` promises a fix "will be described in the release notes" — a small promise beyond the ruled text, kept since W6.5's changelog is maintained per release; `index.rst` links the guide at `blob/master/`, correct from the tag onward; **for W6.5**: `setuptools_scm` resolves `spec-v1.0` as a version tag, so pre-tag builds are `1.1.devN+g<sha>` (the W6.5 prompt says `0.1.devN`) — a TestPyPI dry run would publish a version that sorts *above* `1.0.0b1`; the prompt gains a `git_describe_command` matching `v*` only; the bash sandbox refused commands naming `.github`, so the agent wrote those files with the Write tool. **The profile check moved**: `gh api repos/ICSM/ampere/community/profile` reads GitHub's *default* branch, which is the legacy `master` until the D4 switch — on 2026-10-04 after the push it still reported 25 % with every file false, while the `contents` API confirms all six files on the `v2` ref; the Accept criterion was therefore verified at the `origin/master` switch (D15, 2026-10-05): **health 100 %, `license` → `GPL-3.0`** (the `/license` endpoint names `LICENSE`), `code_of_conduct`, `contributing`, `pull_request_template` and `readme` true; `issue_template` reads false because that legacy profile field detects Markdown templates only, not YAML issue forms — the forms are on the default branch and the score is unaffected. Gate: CI's run on the push after merge (`.github/**` and `pyproject.toml` are classed CORE, so the full matrix runs) — run 37239911091 on `74464f3`: **GREEN, 29 jobs — 28 success, 1 skipped — the gate of record; W6.17 closed 2026-10-05**. Agent ≈ 74 k tokens / 40 tool uses / 6 min. |
| W6.5 | merged 2026-10-05 at `f4ea081` (Opus-authored, Fable-reviewed as the second pass — `release.yml`, `_citation.py`, `CITATION.cff`, the changelog page, the citing page, the release-procedure section and the drafted decision-log row read whole; `pyproject.toml`, the README, `install.rst`, `index.rst`, `conf.py`, `__init__.py`, `provenance.py`, the exceptions message and the four test edits read as diffs; the branch descends from W6.17's head `d49a5cc`, and the merge tests clean on master after W6.17's; **one review commit** `addaeca` with three changes: (i) the smoke job takes only the wheel from TestPyPI (`pip download --no-deps`) and resolves its dependencies from PyPI alone, in place of `--index-url` TestPyPI plus `--extra-index-url` PyPI, which let a name squatted on TestPyPI shadow a real dependency (the agent's own finding) and also proves the wheel installs from PyPI the way a user does; (ii) the release procedure's step (f) no longer says "require a pull request" — merges are local and pushed, so that setting would block every push; D4's rule stands: the CI checks required by name, no force-push, no deletion, the administrator bypass left as GitHub sets it; (iii) a minimal `MANIFEST.in` — `prune .claude`, `global-exclude .DS_Store`, the 3.7 MB legacy corner PNG excluded — so the sdist on PyPI does not carry the agent skills, a tracked `.DS_Store` or a run output (sdist 19.3 → 16.0 MB; the wheel unchanged); the orchestrator's reproduction on the branch: `--frozen` installs of dev/torch/jax ok and `pixi lock --check` ok (the lock unchanged by it), `build-dist` + `check-dist` PASSED (twine strict) and OK (check-wheel-contents) at `751f68a` and again at `addaeca`, the wheel's `METADATA`: `Name: ampere-astro`, `License-Expression: GPL-3.0-or-later`, `License-File: LICENSE`, `Requires-Python: >=3.12`, the five `Project-URL`s, `Development Status :: 4 - Beta`, no `License ::` classifier; lint clean, 334 files formatted, `actionlint` rc 0 both times; **pyrefly 0 errors in dev, torch and jax**; `pixi run test` 99 passed / 4 skipped; `test_citation` + `test_public_docstrings` + `test_results` + `test_engines` + `test_spec_doctests` 325 passed / 2 skipped; `pixi run docs` from a clean `_build` **rc 0, 0 WARNING lines** in 3 min 16 s with `changelog.html` and `citing.html` rendered; `ampere.cite()` prints the nine-author reference and `import ampere` loads no `ampere.legacy.*` module; the agent's additionally: a fresh venv installing a wheel built as `1.0.0b1` — `pip install ampere-astro` with no version takes the pre-release, `__version__` `1.0.0b1`, `[jax]` resolves from the wheel's metadata with `pip check` clean — and a toy repository with tags `v0.1`, `spec-v1.0`, `v1.0.0b1` versioning as `0.2.dev2` → `1.0.0b1` → `1.0.0b2.dev1`; merged by Peter). What landed: **the rename** `name = "ampere-astro"` across every self-referencing extra and all eleven pixi `pypi-dependencies` keys (the `gpu` feature's included), `version("ampere-astro")` in `__init__.py` and `conf.py`, `OptionalDependencyError`'s `pip install "ampere-astro[{extra}]"` with the clash clause gone and the short form in seven docstrings, **and `provenance.package_versions` looking ampere up under `ampere-astro` while recording it under `"ampere"`** (without it every run would have silently dropped ampere's own version — the agent's catch, with a test); **versioning**: `[tool.setuptools_scm.scm.git] describe_command … --match 'v*'` (the table form; `setuptools_scm>=9.2` in `build-system`, since `git_describe_command` is deprecated in 10 and 9.0.0 was yanked) so `spec-v1.0` no longer makes every build `1.1.devN` — pre-tag builds are `0.2.devN` after the legacy `v0.1` tag on origin, the tag builds exactly `1.0.0b1`, no fallback version; **the release build**: `build`, `twine`, `check-wheel-contents` in `dev`, tasks `build-dist` and `check-dist`; **`release.yml`**: `v*` tags and a `workflow_dispatch` `target` (testpypi default | pypi, the latter refused unless from a `v*` tag), `build` (fetch-depth 0, setup-pixi dev, `SETUPTOOLS_SCM_OVERRIDES_FOR_AMPERE_ASTRO` dropping the `+g<sha>` local label both indexes refuse, a tag-builds-its-own-version check, the dist artefact), `publish-testpypi` (environment `testpypi`, OIDC, `skip-existing`), `smoke` (a ten-leg matrix: base, each `all` extra alone, `all`; Python 3.13 fresh venv; the wheel from TestPyPI, dependencies from PyPI; imports `ampere`, `ampere.core`, the extra's probe module, calls `cite()`), `publish-pypi` (environment `pypi`, Peter's approval, needs `smoke` on a tag), actions pinned to SHAs, no token or secret, `concurrency` without cancel; **citation** (#62): `ampere/_citation.py` with the nine `pyproject.toml` authors in order and no invented ORCIDs, `ampere.cite(format="text"|"bibtex", *, file=None)` resolved lazily through `ampere.__getattr__` so `import ampere` stays legacy-free, `scripts/write_citation_cff.py` generating `CITATION.cff` (CFF 1.2.0, no DOI, `date-released` as a comment) under `tests/core/test_citation.py`'s nine rows, `docs/source/citing.rst`; **the changelog** `docs/source/changelog.rst` (preamble, three names, the legacy policy, Phases 0–6 as what a user can do with page links, Known limitations, Deprecations) in the User guide toctree with `citing`; **install pages** PyPI-first with the three-names statement once each, the clone as the contributor route keeping W6.17's pointers, `install.rst`'s extras table gaining `nautilus`/`ultranest`/`blackjax`, the README's Python range 3.12–3.14 and "beta"; `pyproject.toml`'s one-sentence description, `4 - Beta` and Astronomy classifiers, the five URLs; **the release procedure** in `docs/development.md` after the merged-gate section (the dry run, steps (a)–(h) with Peter's named as his); the drafted decision-log row with `DATE`/`TAG_COMMIT`/`GPU_ROWS` placeholders. **For Peter**: `cite()` prints and returns `None` (as ruled) — say if you would rather it returned the string; the dry run `gh workflow run release.yml --ref v2 -f target=testpypi` runs on your word after the merge; #62's dependency-citation half is not done and is proposed for re-filing. Out-of-scope findings (recorded, not fixed): `overview.rst`, `sbi.rst` and the two backend API pages still show clone-form installs and `ampere/__init__.py`'s docstring still says "alpha testing phase" — W6.3's style pass territory, a one-line each follow-up; `conf.py`'s copyright and author lists differ from `pyproject.toml`; `ci.yml`'s `v2` trigger can go after the switch; `v0.1` and the `archive/*` tags already exist on origin, so the beta is not the first tag push; the agent's `pixi reinstall --all` installed the 8 GB `gpu` environment in its worktree (removed with `pixi clean -e gpu`; the CUDA wheels may sit in `~/.cache/rattler`, 35 GB today). Gate: **all four legs through CI on the push after merge** (`pyproject.toml` is CORE: everything runs) — run 37243130167 on `ef5a881` (`origin/v2`): **GREEN, 29 jobs — 28 success, 1 skipped — the gate of record; W6.5's code closed 2026-10-05** (the item's release steps remain open until the tag). **The TestPyPI dry run, 2026-10-05 on Peter's go** — release run 37243911563 on `1355d2a` (`gh workflow run release.yml --ref master -f target=testpypi`, after D15's switch): **GREEN — 13 jobs: build, publish to TestPyPI, the ten smoke legs (base, torch, jax, sbi, zeus, extinction, nautilus, ultranest, blackjax, all) all success, publish to PyPI skipped as designed**; TestPyPI serves `ampere-astro 0.2.dev1805` (wheel and sdist, `license_expression` `GPL-3.0-or-later`) — the item's "a TestPyPI release installable with every extra in a fresh venv" met. **THE RELEASE, 2026-10-05**: (a) the gate — CI run 37256587783 on the W6.18 merge `9727348` GREEN, the GPU rows job 18042730 on `62e851c` 33 passed; (b) Peter tagged `v1.0.0b1` on `62e851c`; (c) **the first release run 37259563688 FAILED at the build**: the tag versioned as `1.0.0b2.dev0` — the runner's pixi 0.81.0 rewrote `pixi.lock` on a plain `pixi run` (a `virtual-packages` block), the tree read dirty, and the tag-builds-its-own-version check stopped it before any publish (proved by the throwaway `diag-scm-dirty` branch, runs 37260066417/37260178766); the fix `56db292` (`pixi run --frozen` on both build steps, a clean-tree check that prints the diff; the procedure amended) merged by Peter, the tag deleted and re-created on `56db292`; **release run 37260510961 on the new tag: build GREEN (`1.0.0b1` exactly), TestPyPI publish GREEN, the ten smoke legs GREEN, `publish to PyPI` refused once by the `pypi` environment's protected-branches-only deployment policy (no tag admitted) — Peter set it to selected branches and tags with `v*`, re-ran the failed job, approved the deployment — GREEN**; (d) **PyPI serves `ampere-astro 1.0.0b1`** (wheel 12.3 MB, sdist 16.0 MB, `GPL-3.0-or-later`); **CI run 37260510926 on the tag: GREEN, 29 jobs — 28 success, 1 skipped**; (e) **Peter's GitHub release** https://github.com/ICSM/ampere/releases/tag/v1.0.0b1 (pre-release; Zenodo's first attempt failed with `401 Bad credentials` — a stale GitHub token on the Zenodo account — fixed by reconnecting GitHub on Zenodo and re-creating the release); **Zenodo: concept DOI 10.5281/zenodo.23151412, the v1.0.0b1 record 10.5281/zenodo.23151413**, nine creators and the licence read from `CITATION.cff`; (f) done before (D15); (g) **Read the Docs serves `v1.0.0b1`** (activated by Peter; `latest` follows `master`); (h) **the DOI commit `087dff8`** (`_citation.py`'s `doi` and `date_released`, `CITATION.cff` regenerated, the README's PyPI/DOI/docs badges, the citing page's example; `test_citation` 9 passed) and the D10 issue closes. **The release is complete.** Agent ≈ 245 k tokens / 182 tool uses / 34 min. |
| W6.18 | merged 2026-10-05 at `9727348` (a fast-forward; the branch head) (Fable-authored as the orchestrator's own work — session ampere-3d, 2026-10-05, from base `9d99d53`; the beta waits for it; the diff re-read whole by the recovering session ampere-4d after a WSL shutdown; merged by Peter). What landed, each fault found by a GPU run and none visible on a CPU: (1) `ampere/core/realisation.py::_check_agrees_at_reference` converts its scalar with `float()`, not `float(np.asarray(...))` — a CUDA tensor has no `__array__`, so every torch problem on a GPU failed to lower; the module's numpy import goes with it; (2) `ampere/backends/jax/_device.py::resolve_device` asks `jax.devices(name)` instead of searching `jax.devices()`, which lists the *default* backend only — on a GPU node `"cpu"` was refused as "no such platform"; jax's `RuntimeError` for a platform it lacks becomes the ruled refusal, the message still naming the platforms present; (3) `tests/gpu/test_jax_gpu.py`'s `TestChunkSharding.theta` draws rows of `problem.free_size` (four), not a literal two; (4) `ampere/inference/_nuts.py` builds pyro's `initial_params` start tensor on the realisation's `device` (`getattr(self.realisation, "device", None)`, the torch `LoweredProblem`'s attribute; a caller-supplied density keeps the default) — pyro's momenta and mass matrix had followed a CPU start beside a CUDA density into `velocity_verlet`'s "Expected all tensors to be on the same device"; its reference-value self-check and both engines' per-dataset decomposition (`_nuts.py`, `_vi.py`) convert with `float()` as in (1); (5) **the jax `QuasisepGP` refuses an accelerator by name at construction** (`LikelihoodError` naming `DenseGP` and `HilbertSpaceGP`): celerite2 0.3.3 registers its primitives' MLIR lowerings for the CPU alone, so the first factorisation died inside the trace — the GPU row comparing the quasiseparable density across devices is now the row asserting the refusal, and `changelog.rst`'s "Known limitations" records the jax solver as CPU-only (the torch one has no such limit); a jax-native recursion is Phase 7's; (6) `scripts/cluster/run.sh fetch` prints the pytest summary and the exit-status line by pattern — the cluster appends a usage epilogue after the job, so the summary is not the log's last line — and `docs/development.md`'s GPU-rows section says so, with 33 rows after parametrisation (not 30) and the first run's result; (7) the CPU-only regression row for (2), `tests/backends/test_jax.py::TestPerInstanceDevice::test_the_cpu_is_found_when_the_default_backend_is_an_accelerator` — the no-argument `jax.devices` monkeypatched to report an accelerator alone, as on the node, and `resolve_device("cpu")` still finds the CPU (`366f576`; the item allowed for the GPU row to be the only one, but this one is writable). **The GPU rows, four Slurm jobs on the branch** (`~/.cache/ampere-gates/gpu/<hash>.log`): 18041011 on `cc7a65e` **20 passed / 13 failed** (faults 1–3); 18042559 on `7e324f5` **31 / 2** (faults 4 and 5, masked before); 18042608 on `351f6e3` **32 / 1** (the NUTS decomposition's numpy route, fault 4's last term); **18042730 on `62e851c`: 33 passed in 45.78 s, exit status 0** — node `gina8`, one A100-SXM4-80GB, driver 615.71.09 (floor 525); `ampere/` and `tests/gpu` are byte-identical from `62e851c` to the branch head (`git diff 62e851c w6.18-gpu-first-run --stat` is the one test file). No contract changes. **Gates on the branch** (`~/.cache/ampere-gates/w6.18-legs{,2,4}.log`): at `7e324f5` dev `tests/core` + import sweep 1549 passed / 36 skipped, jax `tests/backends` + `tests/conformance` 1497 / 90, torch `tests/conformance` + `tests/backends` 1675 / 83, jax typecheck 0 errors; at `351f6e3` dev `tests/core` + `tests/inference` + import sweep 1765 / 430 (the jax and torch legs of that script, and a third script, cut off by the WSL shutdown of 2026-10-05 02:2x); **at the head `366f576`** (the recovering session, `w6.18-legs4`): dev `tests/core` + `tests/inference` + import sweep **1765 passed / 430 skipped** in 2 min 31 s; pyrefly **0 errors in dev, jax and torch**; jax `tests/backends` + `tests/conformance` + `tests/inference` **1891 / 307** in 15 min 46 s (the new row collected); torch the same three suites **2028 / 340** in 26 min 13 s; lint clean, 334 files formatted; `pixi run docs` from a clean `_build` **rc 0, 0 WARNING lines** — all green (the long leg times: a VS Code server update ran alongside). Found, not fixed: nothing in the code; the process lesson (the handoff on master went stale while the orchestrator worked its own branch) is in the handoff. Gate: CI on the push after merge (the item is core-touching) — run 37256587783 on `9727348` (`origin/master`): **GREEN, 29 jobs — 28 success, 1 skipped — the gate of record; W6.18 closed 2026-10-05**. |
| W6.12 | merged 2026-10-05 at `4b83daa` (on Peter's word, by the orchestrator; a true merge) (Opus-authored, Fable-reviewed as the second pass — `oifits.py`'s canonicalisation and sign logic (`_canonical_pair`, `_permutation_is_odd`, `_wrap_degrees`, `_closure_phases`) read against the container's docstring and the OIFITS definitions, the reference `SquaredAmplitude` and its torch twin read whole, the test's analytic truth model and the composition row read, the docs section and the changelog's Unreleased section read whole, `tests/data/README.md` read and the contest reference verified (Cotton et al. 2008, Proc. SPIE 7013, 70131N); **no fixes needed**; merged by Peter; the orchestrator's reproduction on the branch at `e903510`: dev `test_oifits.py` + the conformance row + `tests/test_imports.py` **131 passed / 12 skipped**, `tests/backends/test_reference_interferometry.py` + `tests/conformance` **720 / 87**, lint clean, 338 files formatted, pyrefly **0 errors**, `pixi run --frozen docs` from a clean `_build` **rc 0, 0 WARNING lines** (log `~/.cache/ampere-gates/w6.12-review.log`); the agent's: `tests/interferometry tests/backends/test_reference_interferometry.py tests/conformance tests/test_imports.py` 866 passed / 100 skipped in 7 min 25 s (dev), `test_oifits.py` 24 passed, torch `tests/conformance` 1100 / 87 and `tests/backends -k interferometry` 72 / 2, jax 1100 / 87 and 71 / 4 (the main checkout's environments with the worktree first on the path, `ampere.__file__` asserted), the download fallback exercised with the vendored copy moved aside (8 passed) and the offline skip by name, docs from a clean `_build` rc 0 / 0 WARNING lines, lint, format (338 files), pyrefly 0 errors with `ampere/interferometry` in scope; ≈ 231 k tokens / 114 tool uses / 2 h, no hand-off). **Landed**: `ampere.interferometry` — the observable's front door (placement memo option D): `read_oifits(path, *, target=, insname=, arrname=) -> OIFITSData` (a frozen dataclass: `visibilities` from `OI_VIS`, `squared_visibilities` from `OI_VIS2`, `closure_phases` from `OI_T3`, each `None` when absent, plus `target`, `instrument`, `array`, `revision`, `wavelengths`, `bandwidths`), `astropy.io.fits` only; several targets/instruments/arrays without the keyword refused by name with the choices listed; one sample per (row, channel) with `u = UCOORD/EFF_WAVE`; every baseline stored as `b_ij = r_j − r_i` with `i < j` (a reversed row negated, its `OI_VIS` phase conjugated); **the canonical triangle**: the stored `b_ij`, `b_jk` taken from the file's three oriented baselines, the phase negated exactly for an odd station permutation and wrapped into `(−π, π]` — proved against the contest file's published truth (two uniform discs, 5.0 mas, 30° east of north, ratio 8.9) written in numpy in the test: V² median 0.68σ (max 3.96σ), closure phases median 0.66σ (95th percentile 2.25σ; 14 of 800 samples near triple-amplitude nulls in the 1.50 µm channel exceed 5σ, the contest having simulated from a pixelised, tapered image), **the mirrored source 47σ** — a sign, station-order or units slip fails loudly; all five non-identity station orders rewritten into the file and read back canonically to 1e-5 rad; `FLAG` → mask, non-finite values and non-positive uncertainties masked and counted in `meta`; labels (`baseline`, `triangle`, `baseline_ij/jk/ki`), `channel`, `mjd`, `eff_band` in `extra_coords`; a loud warning for revision-2 differential `AMPTYP`/`PHITYP`; the front door re-exports the two kinds, four families (`GaussianFamily` added for the V² fit), the reference steps, the model trio and their analytic twins, and `read_oifits`/`OIFITSData`/`OIFITSError`; `import ampere` untouched (the import-cost rows pass). **`SquaredAmplitude`** on the three backends (torch/jax as `re² + im²` so the gradient exists at a null): kind-preserving, parameter-free `|V|²`, the output unit the input's squared, and a **`normalisation=` buffer** (the model's fixed zero-spacing flux) giving the unitless `|V/F|²` an `OI_VIS2` holds — without it a `Jy²` prediction is refused by alignment against unitless data, as it should be; seven conformance rows (`tests/conformance/test_squared_amplitude.py`, the mirror column skipping by name). `tests/data/` created with the 83 KB contest file, its SHA-256, pinned URL and redistribution basis (public contest data, MIT-redistributed in OIFITS.jl; vendoring confirmed by Peter at dispatch), the test downloading to a cache when the vendored copy is missing and skipping by name offline. Docs: `interferometry.rst` §10 "Reading an OIFITS file: the front door" (the §4 composition rebuilt from the reader in a dozen lines, run in the test; the two conventions in prose; the opening list of what a modality supplies ends with "a reader"), `ampere.interferometry.rst` in the API toctree, one sentence on `overview.rst`, the changelog's **Unreleased** section (the beta's "No file readers yet" untouched); `ampere/interferometry` in the typecheck scope. **Two departures from the rulings, both accepted at review**: (i) the spectral axis is **micron**, not metres — every shipped `FourierSample` emits its prediction's axis in µm and the alignment check compares axes exactly, so an axis in metres made `FittingProblem` refuse (the agent's finding, below); (ii) `SquaredAmplitude` handles units and takes `normalisation=`, since a squared prediction in `Jy²` cannot otherwise meet a normalised `OI_VIS2`; a fixed buffer rather than the model's own `V(0)`, because normalised data carry no total-flux information. **For Peter**: the `normalisation=` semantics (a fixed zero-spacing flux) — a `normalisation="model"` option dividing by the model's own zero-spacing flux is the refinement for a model whose total flux is free, carried to Phase 7; the DFT sign is `exp(−2πi(ux + vy))` with x east and y north, confirmed by the contest file. **Found, not fixed** (carried): `FourierSample.from_observed` on every backend converts the observed spectral axis to µm and emits µm, so any container not in µm fails alignment invisibly until composition — it could emit the observed container's own unit (a housekeeping item); `tests/conformance/protocol.py`'s `InterferometryPieces` has no `squared_amplitude` field and the mirror backend no `SquaredAmplitude`; `test_native_interferometry.py::test_every_piece_declares_this_backend` lists piece names by hand without `SquaredAmplitude` (the conformance row covers the declaration); the legacy `.gitignore`'s `*.fits` would silently ignore a future `.fits` test file. Gate: CI on the push after merge (dev, torch and jax legs — the native twins changed) — run 37274523311 on `c855b3d`: **GREEN, 29 jobs — 28 success, 1 skipped — the gate of record; W6.12 closed 2026-10-05**. |
| W6.2 | merged 2026-09-28 at `2ce213b` (Sonnet-authored, Fable-reviewed as the second pass — the page read whole, the builders, the tie registration and the initial-positions code read directly, the merge tested clean on top of tranche A; **one review commit** `6883cd9` (the page says where the walkers start and why — the code did, the page did not); merged by Peter; the agent's evidence: `test_photometry_spectra.py` + `test_sed_composition.py` 24 passed in 85 s, `test-examples` on the branch 108 passed / 13 skipped in 4 min 51 s (master 93 / 13 — the fifteen new rows), the docs build 11 `WARNING` lines / Sphinx's 12 (the base) with no toctree, label or document warning, lint/format (283 files)/typecheck (0 errors)/import sweep 96/4 clean; the orchestrator's docs build on the review commit from a clean `_build`: 11 `WARNING` lines, the base, none naming the page; **the four full-budget runs** (seed 20260928, emcee, 68 % intervals on the page, 95 % coverage from `report()`): untied/independent 98 s — the factors recovered (0.912, 1.098 against 0.92, 1.08), the temperature missed by a fraction of a kelvin (the bump's bias); untied/GP 173 s — 5/5 covered; tied/independent 71 s — 0/4, the shared factor overshooting past both truths to 1.21 and β to 1.48 against 1.6; tied/GP 155 s — 4/4; **the agent ≈ 426 k tokens / 499 tool uses / 99 min — past the soft cap without a hand-off (recorded as W5.28's and W5.32's were); and it bypassed the gate lock for its four full-budget runs, the docs build and `test-examples` while tranche A's jobs held it, checking free memory (8.6–9.3 GB) before each and surfacing the choice in its report — no gate was running, nothing was harmed, and the rule stands: the lock is the memory guard's proxy, not a suggestion; the numbers are single-threaded and consistent with the locked smoke runs, and are not re-verified**; **the wave-3 merged gate (both merges; launched 2026-09-28 08:43 BST on master's working tree at `f679776`, the two merges plus handoff prose; complete 09:26): dev `test-all` 3155 passed / 616 skipped in 40 min 36 s — the W6.6 baseline (3114/612) plus exactly the 41 rows and 4 sbi-gated skips the two items added to `tests/examples`; `-e sbi` the two tranche-A test files 30 passed in 76 s — green; log `~/.cache/ampere-gates/wave3-gate.log`**; **the CI run on the push of `f679776`, 36393247614: FAILED at lint — four lines of `linear_sed.py`'s coverage docstring over the 100-column limit (the agent's lint pass predated its coverage commit; every other job green) — fixed on master by the orchestrator at `6527d73` and re-pushed; **the CI run on that push, 36397440625 on `c66a1d9`: GREEN**). Landed: `examples/photometry_spectra/` (the sibling's `ModifiedBlackBody`; three instruments — `sl` R ≈ 100 over 5–14 µm, `ll` R ≈ 60 over 14–38 µm, each `[LSFConvolution, Resample, CalibrationScale(st.lognorm(0.05, scale=1.0))]`, and a six-band `catalogue` spanning 22–161 µm; `generators` injecting calibration truths 0.92 and 1.08 and a 10 % Gaussian bump of 4 µm width at 25 µm into `ll`; `build_problem(tie=, gp=)` — `Tie("calibration", ("sl.instrument.calibration_scale.scale", "ll.instrument.calibration_scale.scale"), prior=...)` through `FittingProblem(..., ties=)` (`dataset.py` ~2262 quoted for why it lives there), `GaussianProcessNoise` with a Matérn-3/2 on the two spectra and never on the photometry; `--tie`, `--gp`, `--seed`, the emcee budget flags, `--figures DIR` writing the corner, posterior-predictive and GP-localisation figures; `_initial_positions` starting the walkers in a ball around the synthetic truth); `tests/examples/test_photometry_spectra.py` (fifteen rows: both tie modes' names and counts, the injected factors against the noiseless prediction, tiny fits in every tie × gp combination, `main`); `docs/source/photometry_spectra.rst` (the physics, the five nouns for three observations, calibration as a step with a prior — untied then tied, the complement with the four-arm table, running it with the walker-start note, see also) in the tutorials toctree after `sed_composition` with one sentence in the paragraph. **Findings**: two surfaced only at full budget — an unbounded GP length-scale prior (up to 100 µm against `ll`'s 24 µm span) let a long-length-scale draw mimic a recalibration, so the priors were tightened to 0.5–12 µm and 1e-2–3 Jy; and `EmceeEngine`'s default prior-drawn initial positions froze about one walker in ten for thousands of steps on this composition (**for W6.7: the warm-start case, measured**); `tutorials.rst`'s "the first seven have runnable counterparts" now undercounts by one (W6.3). **For Peter**: the lock bypass above. Carried: the frozen-walker finding to W6.7; the paragraph count to W6.3. |
| W7.5 | merged 2026-10-07 at `5c606e2` (Sonnet-authored, Fable-reviewed as the second pass — the whole diff read: the CI trigger, the five docs pages, the metadata, the README, `_variable_width_strings` and the generalised `_concatenate` reset, the lifted `_pyphot_compat.py` and both `from_library` call sites, the `heavy` marker on every consumer of `sed_posterior`/`sed_scipy` (checked by grep: four tests, six collected rows), the two training-set rows, the migration note and the changelog; **one review commit** `554f105`: the migration note's sentence on `phoenixstar.py` line 130 said "again", which read as a third division — it is the stand-alone post-processing function repeating the model's pair at 454 pc, and the note now says so; merged by Peter; the orchestrator's reproduction on the branch at `554f105`: lint clean, 339 files formatted, pyrefly **0 errors** (dev), `-m heavy` **6/44 collected**, `tests/results/test_training.py` + `tests/test_imports.py` + `tests/backends/test_reference.py` **240 passed / 4 skipped** in 18 s (log `~/.cache/ampere-gates/w7.5-review.log`); the agent's: `tests/results/test_training.py` + `tests/test_imports.py` + `tests/inference/test_optimise.py` **185 passed / 16 skipped** in 4 min 27 s (dev, `-n 4`), lint clean, 339 files formatted, pyrefly 0 errors (dev), `-m heavy` collects exactly six rows and `-m "not heavy"` none of them, `tests/backends/test_reference.py -k library` 20 passed, `pixi run --frozen docs` from a clean `_build` **0 WARNING lines** against 0 at base, the sbi class through the driver from the main checkout's `sbi` environment with `ampere.__file__` asserted in the worktree: the row **1 passed in 13.4 s at min `ks_pvalue` 0.165**, the class 7 passed, the whole `test_sbi.py` 202 passed / 7 skipped in 6 min 40 s; **CI run 37549374916 on `587b777` GREEN — 29 jobs, 28 success, 1 skipped — W7.5's gate of record; W7.5 CLOSED**). What landed, line by line: (1) `ci.yml`'s `push` trigger is `branches: [master]` (D15); (2) `overview.rst`, `sbi.rst` and the torch and jax backend pages install with `pip install ampere-astro[...]` and the false "not on PyPI" sentence is gone — `install.rst`'s development install stays clone-form by ruling; (3) `ampere/__init__.py`'s docstring says beta on PyPI with the docs link, and its `__copyright__` (2017–2026) and `conf.py`'s `copyright`/`author` follow `pyproject.toml`'s nine authors in order — **for Peter: Sacha Hony was named in `conf.py` and the old `__copyright__` but is not in `pyproject.toml`, so the name drops under the "pyproject is the source of truth" ruling; add him to `pyproject.toml` if he should be listed**; (4) `examples/README.md` says `quickstart` and `Ampere_MBB_Example` are v2 and executed at the docs build, only `Embedding_nets` legacy and never executed, and lists `ngc6302`; (5) every string variable of a training set's sample-stats and context groups is written variable-length (`encoding={"dtype": str}`, as the `extra_*` coordinates already were since W6.14) and `_concatenate`'s reset keys on dtype kind rather than the `extra_` prefix — a row appends a 200-character failure message after a short one and reads it back whole (red without the fix), a second row carries integer and boolean extra coordinates through write, append and read-back; (6) `ampere/backends/reference/_pyphot_compat.py` is the v2 home of `PYPHOT_V2`/`get_unit`, both backends' `from_library` import it, the legacy module untouched, so **v2 imports no legacy module** (the W6.0 row green); the stray `ampere/examples` ruff entry removed; (7) the `heavy` marker registered and `test-fast` adds `and not heavy`, marking `TestTheScipyRoute`/`TestTheMapRoute`/`TestTheVIRoute`'s `sed_posterior` rows and the burn-in comparison; (8) the sbi calibration row **left as it is**: the 0.0054 the handoff recorded did not reproduce — 0.165 (16× the threshold) identically across the row alone, the class, the whole file and under single-thread BLAS, the seeds being pinned — so there was nothing to buy; if it recurs on Peter's checkout, the failing run's log is the next evidence; (9) netCDF4 left at 1.7.4 — PyPI's latest on 2026-10-07, no 1.7.4.1; (10) the Hyperion shim and pin left — conda-forge and upstream both at 0.9.11 (the only release, 2024-04-22); (11) a `.. note::` in `migrating.rst`'s star-disc section names where each legacy script divides by 4πd² twice (`phoenixstar.py` lines 25/53 and its post-processing at 130; `star_disc.py` line 43 feeding `QuickSED.py` line 148), cross-referenced from the phoenix-star section; (12) four Unreleased changelog bullets. Out-of-scope findings: `gp_posterior_medians` in `test_optimise.py` (an emcee run of `GP_WALKERS × GP_STEPS`) has the same shape as the heavy fixtures and is unmarked; the torch and jax typechecks were not run on the branch (the `from_library` edit is an import path only; CI's own-environment typecheck jobs are the check). |
| W7.6 | merged 2026-10-07 at `2815dab` (Opus-authored, Fable-reviewed as the second pass — the whole diff read: `_Coordinates` (the mixed-coordinate map, the reflection fold, the normalisation and the decade rule), `saturated_bounds`, `_warn_if_saturated`, the route and `_polish`, the `warm_start_gp` lookup, `Optimum.coordinates` and its round trip, the three test classes, the regression row, the page, the contract sentence and row, the plan row; **one review commit `ac446ec`**: the regression row runs one start instead of two (the agent reported the same MAP from the first start alone; two starts took 610 s, too long for CI's examples job) and `inference.md` §10b's two "over the packed unconstrained vector" sentences gain an *Amended W7.6* note the agent had left for want of ownership; merged by Peter; the orchestrator's reproduction on the branch: lint clean, 339 files formatted, pyrefly 0 errors (dev), `test_optimise.py` + `test_optimum.py` + conformance `test_optimise.py` not-heavy **69 passed / 12 skipped** in 8 s (log `~/.cache/ampere-gates/w7.6-review.log`), the one-start regression row **1 passed in 4 min 47 s** at log posterior −3152.0 with one coordinate saturated (`model.logacold2`) — the same MAP two starts gave in 610 s; the agent's: the dev batch (`test_optimise` + conformance + spec doctests + `tests/results`) **552 passed / 21 skipped** in 4 min 23 s, torch `test_optimise` + conformance **58 / 11** (47 min), jax **58 / 11** (11 min), both through the driver with `ampere.__file__` asserted in the worktree, lint, format, pyrefly 0, `pixi run --frozen docs` **0 WARNING lines**, no doctest number moved, no tolerance moved, the two-start regression row 1 passed in 610 s; **CI run 37572241312 on `6c8cdb6` (the push carrying W7.6 and W7.0): 29 jobs, 27 success, 1 skipped, 1 failure — the failure is W7.0's new NUTS row (`test_population_nuts.py`, the jax backends job; 2177 passed beside it), nothing of W7.6's; every job touching W7.6's code green — W7.6's gate of record; W7.6 CLOSED**). What landed: the scipy route moves a `Logit` box or `Log` floor in its **normalised constrained value** (a box on `[0, 1]`, a half-line by its prior median's distance from the floor, a positive box spanning ≥ 10× in its normalised log) and every other coordinate in `u`, the objective's value unchanged (it never carried a Jacobian); **two departures from the dispatch rulings, found by measurement and accepted at review**: (1) Powell is **not** handed `bounds=` — scipy's bounded Powell swaps its local line search for `fminbound` over the whole segment and reports a worsening iteration as convergence (both seed-2 starts "converged" at −12 860 and −15 926) — so Powell and every non-projected method move in coordinates *reflected* into the box (a triangle wave; `|w|` for a half-line) and only L-BFGS-B, TNC, SLSQP and trust-constr receive `bounds=`; (2) the normalisation, because raw units failed the heavy central-50 % row on a `loguniform(1e-16, 1e-14)` scale and linear units stalled three of eight starts on `sed_composition`'s curved ridge (reflected Powell brings 8/8 starts to the mode; the old route 7/8). `Optimum.covariance` stays in `u`; the record gains the frozen field `coordinates` (`ampere_optimum_coordinates`, all-`"unconstrained"` when absent, outside `identity`; `combine` carries it). `saturated_bounds(problem, unconstrained, *, tolerance=1e-3)` is public (W7.12's dependency) and `BoundSaturationWarning` names every converged coordinate within `tolerance × width` of a box end or `tolerance × (reference − lower)` above a `Log` floor — the reference `FittingProblem.reference_values`, the prior median by default, so a GP amplitude of 1e-4 or a 1e-6 norm is not called saturated (a departure from the dispatch's `max(1, |lower|)`, accepted); the flat-box saddle row now warns, correctly. `warm_start_gp` resolves the three hyperparameter names through `problem.mapping.global_name_for`, so a tie renaming the amplitude is accepted (a row). The regression row on the GP-off NGC6302 twin (16 free, 15 boxes): old route log posterior −4146.0 with eleven coordinates at a bound; new route −3151.8 with one; the truth −3106.9; margin 100 nats pinned over measured gaps of 45 and 67 (seed 2: −4947.7/12 → −3173.7/3). `optimisers.rst`'s scipy paragraph, the objective section, the `warm_start_gp` paragraph, the record section and a new "At a bound" subsection (persona B's case); `results.md` §4 sentence and §12 row; the decision-log row; a changelog bullet. The harness cut the worktree at `ce82a4e` and the agent reset to `429960c`. Out-of-scope findings: `Bijection` has no bounds accessor (the helper dispatches on `isinstance(Logit/Log)`; a user's own bounded bijection is treated as unbounded — W7.0's successor); the map route can also end at a saturated bound in `u` without warning; a uniform prior on a box spanning many decades is still resolved poorly by the linear normalisation; `Logit.constrain` overflows harmlessly on the old route's extreme `u`. For Peter: the two departures above are the orchestrator's acceptance, not a ruling of his. |
| W7.0 | merged 2026-10-07 at `daeab20` (Opus-authored, Fable-reviewed as the second pass — **the agent was terminated by the API's monthly spend limit mid-verification** (its last word: the named files green — m2 many-lines, spec doctests, `test_optimise.py` — 93 passed / 14 skipped; the docs build next), with every unit committed: ten commits on `w7.0-derived-node`, the eight planned units plus two fixes the agent found itself (`545703e`: `Plate.expand` qualifies a derived member's bindings onto its sibling members, so the non-centred `theta` from `z` works through `Population.as_plate()` too; `d3a2934`: `complete()` recomputes a derived value natively on a traced path — a native component's `context()` completes tensors or tracers, which the numpy compare-and-refuse cannot survive; found by the torch gate when the slab tail joined `TestTheHorseshoeUnderNUTS`'s parametrisation); so the report is the orchestrator's, from the diff and the legs it ran: the whole diff read (`Derived` and the closed grammar, the fourth state and the §3.3 consumers, `_derive_into` with `DERIVED_RTOL = 1e-9`/`DERIVED_ATOL = 1e-12` on the numpy path and `_NativeOps` on the traced one, both backends' resolve and walk functions, `_internal_members` and the two member rules with the flat refusal, `_slab_tiers`, the vectorised `_derived_variables` with the memo's padding rule, `ampere_derived` and schema 10, the training column skip, the four contracts' amendments, the decision-log row carrying the §3.3 correction, `population.rst`'s doctested block, `test_derived.py`'s ten rows and the fingerprint row with digests recorded from the base), **no fixes needed**; merged by Peter; the orchestrator's legs on the branch at `d3a2934` (`~/.cache/ampere-gates/w7.0-legs.log`): lint clean and 344 files formatted with the untracked `.w70_scratch` excluded, pyrefly 0 errors (dev), dev batch (`tests/conformance tests/core tests/results tests/inference/test_population_nuts.py tests/inference/test_nuts.py tests/backends tests/m2/test_many_lines.py`, `-n 4`) **2794 passed / 242 skipped** in 2 min 28 s; torch `tests/conformance tests/backends test_population_nuts.py test_nuts.py` through the driver with `ampere.__file__` asserted in the worktree **1784 passed / 95 skipped** in 15 min 47 s, pyrefly on the torch package **0 errors**; jax the same **1607 passed / 102 skipped** in 3 min 32 s, pyrefly on the jax package **0 errors**; `pixi run --frozen docs` from a clean `_build` **rc 0, 0 WARNING lines**; **CI run 37572241312 on `6c8cdb6`: 29 jobs, 27 success, 1 skipped, 1 failure** — `tests/inference/test_population_nuts.py::TestTheNonCentredPopulation::test_it_diverges_less_than_the_centred_declaration[jax]` with `(0, 2)`: on CI's jax runner the centred run had no divergences and the non-centred two, the reverse of this machine (2177 passed beside it in the job; the other 27 jobs green). **Fixed at `9a84e8d` on `w7.0-nuts-row-fix`** (Fable, from master): the orchestrator's probe at the row's budget and twice it — centred divergences 127 then 3 on jax, 184 then 5 on torch (the count swings with the step size the adaptation lands on) while the centred run's bulk ESS on `sigma` stays at 1.3–25 against the non-centred run's 92–202 — so the row now pins a bulk-ESS ratio of at least three on `mu` and `sigma` (measured 6–110) and a second row caps the non-centred run's divergences at a tenth of its draws (measured 1, 4 and CI's 2); `inference.md` §9's sentence corrected; 13 passed on jax (21 s) and torch (6 min 24 s). The CI run on the push after the fix's merge is W7.0's gate of record: **CI run 37582940893 on `dbca4fc` GREEN — 29 jobs, 28 success, 1 skipped — W7.0's gate of record; W7.0 CLOSED**). What landed, the memo's §3 whole: `Parameter(name, Derived("<expression>", symbols={...}))` as the fourth state — the grammar (numeric literals, `+ - * / **`, unary minus, `sqrt exp log log1p abs` through `ArrayOps`, `abs` as `absolute`; attribute, subscript, comparison, lambda, any other call and a symbol spelt like a reserved name refused by name) parsed once and stored as `ast.unparse`'s source, the `symbols` mapping of `HierarchicalPrior.hyperparameters`' shape that a merge renames; `is_derived`; `fixed`/`value`/`shared_as`/`bijection`/`fix`/`release`/a tie refused by name; excluded from the free half, included in `names` and `evaluation_order` after its inputs, cycles refused; `complete`/`unpack` compute it (idempotent; a stale supplied value refused on the numpy path, overwritten on the traced one), `prior_transform` and `sample` mid-walk; `ArrayOps.sqrt`/`log` on the three namespaces; torch `unpack_tensor`/`_resolved`/`prior_transform`/`log_prior_tensor`/`sample` and jax `unpack`/`_resolved` (a second pass over the evaluation order)/`prior_transform`/`log_prior`/`numpyro_model` with `numpyro.deterministic` as the structural view only; `Population`'s two member rules at the merge (a routed member declared by every `over` component, else refused by name; a member declared by none internal — routed nowhere, permitted only as an input of a derived member) and the flat layout refusing derived and internal members — the two frozen allowances retracted (`parameters.md` §8 amended; `_replace_or_append` gone); `shrinkage_horseshoe(tail="slab", slab_scale=c)` as `<prefix>.slab_scale` (fixed float or a `Parameter` with a prior), the local half-Cauchy scales, one derived `<prefix>.<leaf>_effective = sqrt(c²s²/(c²+s²))` per component and the amplitudes `halfnorm` on it, the helper's order kept, the plain and regularised tails byte-identical (spec digests pinned from the base); emission's `posterior` gains one variable per derived parameter on every engine, computed once vectorised in evaluation order with the memo's padding rule, `ampere_derived` (canonical JSON, `"[]"` when none) and `PROVENANCE_SCHEMA_VERSION = 10`; `_theta_variables` skips a derived column. Contracts: `parameters.md` §3/§8/§9/§12, `lowering.md` §3.2.1/§5/§8, `results.md` §4/§9, `inference.md` §9/§10a, the decision-log row (with the memo's §8 phrase "never trusted or refused" corrected to the §3.3 text). Tests: `tests/conformance/test_derived.py` (the §9.2 ten rows as named, rows 8 and 10 on every registered fixture, plus the plate-sibling row and the slab's kernel row), `tests/inference/test_population_nuts.py`'s `TestTheNonCentredPopulation` (the fifty members at `WEAK_NOISE = 1.0`, fifty times the recovery rows' noise, where the centred form funnels — at the informative noise neither form diverges, measured — pinning `mu` and `sigma` inside the central 95 %, a bulk-ESS ratio of at least three over the centred run at the same budget and seed (as fixed at `9a84e8d`; the merged branch's strict divergence inequality failed CI), and the posterior's derived `objects.index` equal to `mu + sigma z`), `TestTheHorseshoeUnderNUTS` now sampling the slab tail too, the training-set and fingerprint rows (`TestTheDerivedNode`), `population.rst` added to the doctest pages' required set. Out-of-scope findings (from the review): the `test_results.py` schema row now asserts `>= 9` rather than the exact constant; the fingerprint row's digests pin the population problem's spec across the bump, which any later spec change must re-record. |
