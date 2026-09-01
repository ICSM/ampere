# Ampere v2 — Work Items, Phases 0–1

Companion to `DEVELOPMENT_PLAN.md` (the source of truth for architecture and
decisions — read it first). Each item below is sized for delegation to a
development agent. Phase 2+ items are written only once the Phase 1 spec
freezes (W1.13).

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

### W0.5 — Broken-module quarantine (issues #74–77) [M]
Fix the astropy `BlackBody` API usage (#77); fix `CCMExtinctionLaw` syntax
errors (#76); replace the abandoned `extinction` dependency with
`dust_extinction` or inline curves (#75); fix the `extinctionModels` import
error (#74). Make heavy/optional imports lazy so no broken or missing
optional corner can block `import ampere`.
**Accept:** `import ampere` succeeds in a minimal-deps env; a test imports
every non-quarantined module; quarantined corners raise informative errors
on use, not on import.

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

## Status

| Item | Status |
|---|---|
| W0.1 | done 2026-09-01 (committed directly to master — the pending fix existed only in the local working tree, so the branch-per-item rule was waived for it) |
| W0.2 | merged 2026-09-01 |
| W0.3 | merged 2026-09-01 |
| W0.4–W0.8 | not started |
| W1.1–W1.13 | not started |
