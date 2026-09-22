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


## Phase 5 — Scale-out & advanced inference (drafted 2026-09-15 by Fable; **approved by Peter 2026-09-15 with rulings D1–D8 as recommended — D1 as (a) and (b) together, carried by W5.20; D8 as two Opus items**; wave 1 dispatched 2026-09-15)

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

### W5.16 — RHMF exploratory trial [S; Sonnet]
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

### W5.18 — CI/CD Phase 5 [S; Sonnet]
The new extras (W5.14) as path-gated legs; the GPU job kept skip-clean; the
sbi leg's budget re-measured after W5.10; `test-fast` (the token-economy
proposal's last rule) as a pixi task excluding the `m2_full`,
`interferometry_full` and SBI-training rows, with each item's Accept line
naming its suites. **Depends:** W5.14. **Accept:** `actionlint` clean;
the path table in `path_filters.py` covering the new files; `test-fast`
under 4 min in dev; the docs' CI section updated; no five-suite gate.

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
updated to the new leg times. Minutes are free on the public repository,
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
or rank the derived groups too (W5.7's carried); (c) a `summarise`-side
helper reporting the evidence `B` (W5.0/W5.13); (d) `QuasisepGP.condition
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
first, then W5.26, W5.27 (before W5.19), W5.28 (after W5.27), W5.29,
with W5.21 and W5.24 as before. File ownership per wave in
the dispatch prompts; the sole shared file across waves is
`core/likelihood.py` (W5.4, W5.7, W5.9 each own a section — the solver,
the kernel, the noise-model — and merge in that order).


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
| W5.11 | merged 2026-09-18 at `c81f975` (Opus-authored, Fable-reviewed as the second pass — the thirteen-entry code table, the unit-to-code rule through `astropy.units.get_physical_type` with the one declared exception (spatial frequency, named by a kind putting `rad**-1` in its `AxisSpec.equivalent_units`, coded 8 whether the axis is given in wavelengths or in inverse radians, while a `u` in metres is honestly a length), `ENCODING_VERSION` in the hashed layout record, and the shared by-name refusal `axis_identity_complaint` at both doors read; the orchestrator's re-run of `tests/core/test_encoding.py` + `test_simulate.py` + `tests/results/test_training.py` in dev on the branch: 242 passed in 46 s (and `test_sbi.py` in sbi 179/2 in 5 min 20 s); **no fixes needed**; merged by Peter; the agent's runs: dev targeted 242 passed, sbi `tests/inference/test_sbi.py` 179/2 in 6 min 17 s, lint/format/typecheck (dev, sbi, torch) clean, import sweep 86/4; the merged-wave gate re-run from 2026-09-22 after the 18 September log was lost to a reboot — dev 2955 passed / 446 skipped in 28 min 17 s on `264fbd4`; sbi 4046 / 101 in 56 min 19 s; torch and jax to be recorded; the agent's usage is W5.10's row's — one agent did both items in sequence). Landed: `AXIS_TYPE_CODES` (contract: `absent` 0 for a padded column, `unknown` 1, then dimensionless, length, frequency, energy, time, angle, spatial frequency, wavenumber, speed, temperature, mass — a code may be appended, never renumbered), `AXIS_TYPE_NAMES`, `axis_type_code(unit, spec=)`; `DatasetLayout.axis_codes`/`axis_type_names` in the layout, its dict form, its hash and `compare`'s report; the `axis_identity` column group (`encoding.md` §3 group 4, width A, the code on every row of its dataset's block) exposed by `unpack` on `Unpacked` and every `DatasetView` and carried by `features`, so no embedding wrapper changed; `ENCODING_VERSION` 2 with `EncodingLayout.from_dict` and `append_training_set` refusing a version-1 record by name through one message. Measured: the module docstring's one-dataset `Spectrum` layout `e73db8a8032273f4c511142f1459c51c` (17 columns, version 1, computed from `dd79aa0`'s own code) → `822f4c900af7a854488b2044db598e26` (18 columns, version 2). Deliberately narrower than what it closes: the code names a physical type, not a role — two `length` axes in one collection (a wavelength and a baseline) still share a code, and a kind wanting them distinguished needs a new code, not a new mechanism. |
| W5.10 | merged 2026-09-18 at `d48c291` (Opus-authored, Fable-reviewed as the second pass — `ContextPrior` and the three shipped priors' draws, `Dataset.contextual_observed` as the single route of a context's σ into every family's `sample`, the `"<stream>.context"` sub-stream spawned by index (so every pre-existing budget's θ and noise are bit-for-bit unchanged), the three refusals by name, the FiLM branch's zero-initialised last layer and the pinned calibration inequalities read; **one fix at review**: the context prior did not enter the artefact cache key, so a network trained without a context would have been served for a context-amortised run with the same settings (§7's stale-artefact trap) — fixed as `artefact_key(context=<digest>)`, omitted when `None` so every pre-W5.10 digest is unchanged (verified byte-for-byte against master's own `artefacts.py`, `3af107dc…` either side; W3.12's pinned NPE digest untouched), with 7 store rows and 3 engine rows; the orchestrator's re-runs on the branch: dev `test_encoding.py` + `test_simulate.py` + `test_training.py` 242 passed in 46 s, sbi `test_sbi.py` 179/2 in 5 min 20 s (pre-fix), dev `test_artefacts.py` 63/3 in 4 s (post-fix); merged by Peter; the agent's runs: dev 242, sbi `test_sbi.py` + `test_artefacts.py` 248/2 in 6 min 14 s, torch `test_simulate.py` 128 in 56 s, lint/format/typecheck (dev, sbi, torch) clean, import sweep 86/4; the merged-wave gate re-run from 2026-09-22 (the 18 September log lost to a reboot) — dev 2955 / 446 in 28 min 17 s; sbi 4046 / 101 in 56 min 19 s; torch and jax to be recorded; the agent, W5.11 and W5.10 together over two segments: ≈ 566 k tokens / 462 tool uses / 2 h 38 min, most of it queued behind the wave-4 gate's four legs on the lock). Landed: `ampere.core.simulate.ObservationContext` (per-dataset σ arrays plus a JSON-safe record), the structural `ContextPrior` (`draw(rng, observed)`, `describe()`), `ScaledSigma`, `SigmaArchive`, `SignalToNoise`; `Dataset.contextual_observed(sigma)` and `draw_observation(..., sigma=)` — no `LikelihoodFamily.sample` signature moved; `Simulation.context` (failures included); `simulate_many(context=...)` drawing one context per simulation on its own sub-stream, materialised a chunk at a time; the training set's optional `context` group and `ampere_simulation_context`, append refusing a context mismatch in both directions; `SBIEngine(context=...)` with `ampere_sbi_context`/`ampere_sbi_context_hash` and the prior's digest in the cache key; `calibrate(context=...)` inheriting the run's prior by default; FiLM (`embedding={"type": "set", "film": True}`), off by default, the identity before training; `inference.md` §12/§13, `encoding.md` §3/§8, `results.md` §9/§11 amended, `sbi.rst` gains "Amortising over the observation context". `PROVENANCE_SCHEMA_VERSION` stays 8 and `TRAINING_SET_SCHEMA_VERSION` 1 (engine-specific extras; an optional group). Measured: a 1500-draw NPE under `ScaledSigma(0.5, 2.0)` — TARP area-to-curve −0.007 at a ×1.5 rescale the prior covers, −0.046 at ×12 it does not, pinned as inequalities at about half the separation; the whole of `test_sbi.py` in sbi in 6 min, so no marker. **For Peter**: a context on a problem with a joint noise group (W5.9) is refused by name rather than guessed (one shared context per group from its first dataset is the implementable alternative) — confirm. **Ruled by Peter 2026-09-22: confirmed** — the refusal stands and the shared-context-per-group route stays open until a problem needs it. Carried: the native batched path refuses a context by name (a per-draw σ does not reach the realisation; **scheduled as W5.29**, 2026-09-22); `SignalToNoise` shapes σ from the *observed* magnitudes because a context is drawn before the forward model runs (documented; a prior tracking each draw's own signal needs `draw` to receive the predicted containers — a protocol change); `_concatenate` in `training.py` indexes the addition by the existing tree's groups (a bare `KeyError` or a silent drop for any future optional group); `TrainingSet.contexts` defaults to `()`, so "written before W5.10" and "written without a context" are indistinguishable in Python (the file distinguishes them); sbi 0.27's outlier warning on the set layout under a wide prior (`z_score_x` guidance for the tutorial); a pre-existing wrong stream name in `_calibration_batch` fixed on the branch because W5.10 made it load-bearing. |
| W5.8 | merged 2026-09-18 at `2a1f596` (Opus-authored, Fable-reviewed as the second pass — the horseshoe's global-local chain as three `HierarchicalPrior` levels, the `with_shrinkage` hook (functional; verified directly: the original kernel's spec hash and parameter names unchanged, the shrunk kernel's hash moved, the scale levels registered ahead of the amplitudes, `terms` and `FAMILY` preserved), the pinned margins against the measured values, the SBC contrast re-pinned as the visibility precedent does, and the driver's own rows read; the orchestrator's re-runs on the branch through the lock: dev `tests/core/test_parameter.py` 158 passed in 3.6 s, the driver/arm rows of `tests/m2/test_many_lines.py` 17 passed in 2 min 14 s; **no fixes needed**; merged by Peter; the agent's runs: dev `tests/m2` 112/32 in 11 min 24 s, `test_many_lines.py` 29/2 in 4 min 41 s, `test_many_lines_calibration.py` 6/1 in 2 min 49 s, torch NUTS-over-the-knots 2/2 in 4 min 51 s, jax 2/2 in 1 min 1 s, lint/format/typecheck (dev, torch, jax) clean, import sweep 86/4, docs build clean; the merged-wave gate re-run from 2026-09-22 (the 18 September log lost to a reboot) — dev 2955 / 446 in 28 min 17 s; sbi 4046 / 101 in 56 min 19 s; torch and jax to be recorded; the agent: ≈ 503 k tokens / 82 tool uses over two segments, killed once by the session limit with two edits preserved as a WIP commit). Landed: `ampere.core.regularised_horseshoe(amplitudes, *, prefix, global_scale, tail, unit)` — `τ ~ C⁺(0, τ₀)`, `s_j | τ` a `gamma(½, scale=τ)` (`tail="regularised"`, the default, which lowers) or `C⁺(0, τ)` (`tail="cauchy"`, the plain horseshoe, which does not lower on torch/jax yet), `a_j | s_j ~ N⁺(0, s_j)`, all with `Log` bijections; `ampere.core.with_shrinkage(kernel, declaration)`, the one hook, appended to `kernels.py`; the `many_lines` scenario (kept out of `SCENARIOS`, so the milestone's four are untouched), the `warped` and `sum` arms in `build_kernel`, `examples/m2_misspecification/many_lines.py` and the driver's `--scenario many_lines --kernel {matern32,warped,sum}` (a trailing-block `KeyError` on the warped arm fixed in scope); `tests/m2/test_many_lines.py` (29 rows) and `test_many_lines_calibration.py` (SBC with a parameter-free `LineForestError` step in the simulating problem), the backend-agreement row for the warped arm; `likelihoods.md` §6 (line 889's promise discharged) and §13, `parameters.md` §9, `lowering.md` §3.2.1, `m2_misspecification.rst` (two sections) and `kernels.rst`. Measured at 200 points, bias in posterior widths: standard 35.70 (1/4 covered), stationary Matérn 2.84 (3/4; pinned ≥ 2.0 and > 1.5), warped 0.92 (4/4), `Sum` 0.78 (4/4; each pinned ≤ 1.5 and ≥ 2.0× better than stationary, measured 3.1× and 3.7×); all three localise inside the band (in-band contrast 2.76×, 1.90×, 3.02×, pinned ≥ 1.4×); the horseshoe row `min(a)/max(a)` 0.254 → 0.087 (2.93×, pinned ≥ 1.3) and mass below a tenth 0.257 → 0.537 (2.09×, pinned ≥ 1.2), thresholds a third below the worst of three seeds; SBC at nominal 0.90: standard [0.00, 0.80], warped [0.80, 0.90], gap pinned ≥ 0.4; the warp's fitted increments give an effective length scale of 0.028 µm under the smooth error and 0.0009 µm under the forest, the sign change at the band edge. **For Peter**: (1) `halfcauchy` is in neither backend's lowering table, so the plain-horseshoe tail and the global level are reference-only on NUTS — two lines per backend plus a §3.2 row and a conformance row, a Sonnet filler; (2) a `Derived` parameter node (a pure function of other parameters) would make P&V's slab declarable and let `WarpedKernel` declare its non-centring — a design question, not scheduled; (3) `tests/m2` in dev now costs 11 min 24 s per environment (was ≈ 100 s): `many_lines.SHRINKAGE_EMCEE` 3.5 min, the four arms 2.4 min, the calibration row 2.8 min, each set against a measured scatter — accept, or name the lever; (4) the helper is named `regularised_horseshoe` while its default tail is a gamma variant of the horseshoe rather than Piironen & Vehtari's slab — the docstring says so; keep the name, or rename. **Ruled by Peter 2026-09-22**: (1) scheduled as W5.25; (2) a `Derived` node goes to Phase 6 as a design-horizon note; (3) accepted for this gate, the trim is W5.26's; (4) renamed `shrinkage_horseshoe` with `regularised_horseshoe` a deprecated alias for one phase, W5.27. Carried: `with_shrinkage` rebuilds the parameter set through `copy.copy` and a `__dict__` write rather than a constructor path; `arviz_base`'s chain-longer-than-draw warning on every emcee run (pre-existing). |
