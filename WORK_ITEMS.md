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

## Phase 3 — The SBI layer (drafted 2026-09-08 by Fable; **approved by Peter 2026-09-08** with W3.1 revised into two slices from his three notes, and W3.8–W3.10 added from his rulings; W3.1's revised text approved 2026-09-09; dispatch open)

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
**Accept:** encode/decode round trip on every container fixture in
`tests/core` including a masked, an irregular and a complex one; two
datasets of different lengths encode into one padded tensor whose mask
column is exact; a layout mismatch is refused by name; `SBIEngine(...,
embedding="set")` trains and samples on a two-dataset toy problem in the
`sbi` environment; `dev` and `sbi` gates green; lint/format/pyrefly clean.

### W3.4 — TMNRE: reviving the swyft implementation [M; Opus; conditional]
Plan §5: "revive the swyft TMNRE implementation from
`docs/design/harvest/swyft/`". **Gate first**: the item begins with a
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
