# Ampere v2 — development handoff notes

Audience: humans joining or reviewing the v2 redevelopment, and the person
dispatching agents. Agents themselves should start from `AGENTS.md`.

## The three documents

1. `DEVELOPMENT_PLAN.md` — decisions, architecture, phases, traps. The
   source of truth; changes to decisions go through its decision log.
2. `WORK_ITEMS.md` — agent-sized items for Phases 0–1, the working
   agreement, and the status table (updated at merge, not by agents).
3. `AGENTS.md` — ground rules and environment for any coding agent
   (Claude Code reads it via the `CLAUDE.md` symlink; Codex reads it
   natively).

## Dispatching work

- Each item is dispatched to one agent, on one branch (`w0.2-repo-hygiene`
  style), producing one PR. Respect the **Depends** lines in WORK_ITEMS.md;
  items with no unmet dependencies can run in parallel.
- Claude Code: point a session/agent at the item id — the `work-item`
  project skill carries the procedure. Codex: `codex exec` with a
  self-contained prompt; see `.claude/skills/delegate-codex/SKILL.md`.
- Agents do not push unless told to; review happens locally or on pushed
  branches at Peter's discretion.
- **Harness trap (hit and solved 2026-09-01)**: subagent worktrees may be
  created from the *session-start* commit rather than current HEAD. Every
  dispatch prompt must state the intended base commit and the recovery:
  verify the tree is clean and HEAD is a strict ancestor of the target
  (`git merge-base --is-ancestor HEAD <target>`), then branch directly
  from the target; STOP and report otherwise. Both W1.5/W1.6 agents
  recovered cleanly with this.
- **Harness trap (hit and solved 2026-09-03)**: `!`-prefixed commands
  typed by Peter run in the orchestrating session's *persistent shell
  cwd*. If that shell was last in an agent worktree, a `git merge` typed
  there silently no-ops ("Already up to date"). The orchestrator must
  return its shell to the main checkout before handing over any git
  command to run.

## Review checklist (per PR)

- Acceptance criteria evidenced in the report, not just claimed.
- Diff stays inside the item's scope; frozen legacy modules untouched
  (`ampere/data`, `ampere/models`, `ampere/infer`) unless the item says so.
- No run outputs, binaries, or checkpoints; British English in prose.
- From Phase 1 on: conformance suite green; contract changes carry a
  decision-log entry.
- On merge: update the WORK_ITEMS.md status table.

## Environment

- **pixi is the supported route** (W0.8, landed): from a clean clone with
  only `pixi` installed, `pixi install -e dev` sets up the daily-use
  environment (Python 3.13 + the `dev`, `zeus`, `extinction` extras).
  pyproject.toml's `[tool.pixi.*]` tables wrap `[project.dependencies]` /
  `[project.optional-dependencies]` — pixi never duplicates a dependency
  pyproject.toml already declares, it just adds Python (from conda-forge)
  and installs `ampere` itself as an editable PyPI dependency with the
  extras a given feature needs. `pixi.lock` is committed so environments
  are reproducible; `.pixi/` (the installed environments themselves) is
  not, and is gitignored.
  - Common tasks: `pixi run test` (fast: `tests/test_imports.py`),
    `pixi run test-phase1` (core + results + conformance in one pytest
    process — the main gate since W0.10; also `test-core`,
    `test-results` and `conformance` individually),
    `pixi run test-characterisation` (legacy still works),
    `pixi run lint`, `pixi run format-check`, `pixi run typecheck`,
    `pixi run docs`. Since W0.10, plain `pixi run <task>` equals
    `-e dev`.
  - Other environments (`pixi run -e <env> <task>`): `test-py311` /
    `test-py312` / `test-py313` (the CI matrix, one Python each); `sbi`
    (adds the `sbi` extra — torch is a large download, so it is not part
    of `dev` or the `test-py3*` environments); `torch` / `jax` (Phase 2
    backend placeholders, not exercised by anything yet).
- Plain-pip alternative: `pip install -e ".[dev]"`, Python ≥ 3.11.
- Conda env `ampere` (Python 3.13) has an editable install pointing at the
  main checkout — worktree-based work that needs importing its own changes
  should use pixi or `pip install -e .` into a fresh env instead.
- `pytest tests/characterisation` (or `pixi run test-characterisation`) is
  the "legacy still works" gate; run it before merging anything that
  touches shared files.

## ⚡ Pick up here (2026-09-07 session — the two W2.1 rulings landed; backend tracks next)

**In flight**: nothing. **W2.4 slice 1 (torch) and W2.5 slice 1 (jax) are MERGED** (`6964b80`, `3458657`; status rows in WORK_ITEMS.md carry the review records and every carried finding). Gates on merged master (`pixi run -e <env> test-all`, five suites in one process): **torch 1746 passed / 10 skipped, jax 1662 / 7, dev 1369 / 5** — every skip is a quasiseparable conformance row awaiting slice 2 or a backend-only test module declining without its extra; lint/format clean; pyrefly 0 errors in all three environments (each backend package checked with its real types in its own environment, excluded where its library is absent). Nothing pushed to origin.

**Rulings given by Peter 2026-09-07 (evening)**: the realisation registry approach is approved and **prototyped on branch `w2.13-realisation-prototype` at `9b26c3f`** (core `ampere/core/realisation.py`: `register_realisation`/`realise`/`registered_realisations`, checked on use against the contract path; jax registers `lower_problem` at import; `NUTSEngine(problem)` needs no density argument; 12 core rows + 4 NUTS rows; dev/jax suites and both typechecks green) — W2.13 turns it into the contract (§4.5 text, decision-log row, a conformance row per registered realisation, the torch realisation, and the fold-ins). `provenance_config` on `GPSolver` approved. CI stays CPU-only for torch/jax; **GPU smoke tests later** (API-works-level only, trusting CPU/GPU library parity). Dispatch order as proposed: W2.13 → both slice 2s in parallel → W2.7/W2.8/W2.9 as capacity allows.

**W2.13, W2.4 slice 2 and W2.5 slice 2 are ALL MERGED** (`0e7c97c`, `a529129`, `dd2c479`); status rows carry the review records. Gates on merged master (the reconciled tip, re-run by Fable): **torch 1924 passed / 19 skipped, jax 1806 / 19, dev 1419 / 69**, lint/format clean, pyrefly 0 errors in all three environments. **In flight: nothing.** Nothing pushed. Both modern backends now have: native QuasisepGP (celerite2 on both, by measurement), the realised path over every core family/censoring/both solvers, NUTS and VI through `realise`, batching, float32/device opt-ins, schema 5 provenance. Memory note for the orchestrator: two five-suite gates running concurrently exhaust this machine's 13 GB — run them one at a time, torch detached (its pyro NUTS rows are slow).

**Next session**: (1) **Peter's rulings** — W2.14 (the latent-GP fix, below; blocks both tracks' slice 3) and the noise-model papercut (below). (2) Dispatch **W2.14** (Opus, small, core + a conformance row). (3) **W2.7, W2.8, W2.9** — dispatchable now (prompts in this session's scratchpad; W2.7/W2.8 split `plots.py`, so sequential or with the split stated). (4) Then W2.11 (CI: torch/jax/no-extras jobs — the three environments exist; the benchmark harness choice) and W2.10 (M2). Slice 3 of the backend tracks (latent GPs once W2.14 lands, per-instance `device=` with the GPU smoke-test item, a torch/jax `conditional_loo` parity decision, `complex_gaussian` in Phase 4) waits behind those.

**RULING NEEDED — a science bug on the reference path (found by W2.4 slice 2, reproduced by Fable on master `eee491f`).** The latent-GP likelihood is **exactly flat in the kernel hyperparameters**: on a Poisson + `GaussianProcessNoise(Matern32, DenseGP)` problem, `problem.log_likelihood` is bit-identical for amplitude 0.5/5/50 and length scale 0.1/1/100 (`-23.785479914183670` every time), while both are free parameters. Cause: nothing on the scoring path calls `GPSolver.latent_transform` — `Dataset.log_likelihood_of` passes the whitened `z` straight through `GaussianProcessNoise.noise_params(latent=latent)` to the family, which uses it as `f`; the only caller of `latent_transform` is `GaussianFamily.sample`. `likelihoods.md` §~1414 reads "`latent_transform` [is] the declaration and the transform; sampling `f` is Phase 2's", which explains how it was missed: the transform was declared for the backends and never wired on the numpy path that is their oracle. A latent fit today converges to the prior on the GP hyperparameters and reports nothing wrong. **Proposed (W2.14, small, core):** `GaussianProcessNoise.noise_params` applies `f = solver.latent_transform(kernel, coordinates, z, values)` when `latent` is given, so `noise.latent` means `f` as every family already reads it — the solver and coordinates are already in hand there, and no family changes; a `tests/core` regression row asserting the likelihood *moves* with amplitude and length scale, and a conformance row per backend so the native latent paths (both refused by name today, torch's with a test that asserts the flatness so the refusal lifts when core is fixed) are held to the corrected oracle. A §4.4 clarification with a decision-log row. Dispatch to Opus before either track's slice 3.

**One user-facing consequence for Peter to confirm or soften**: fold-in 7 as ruled makes core's `IndependentNoise` declare `reference`, so `Dataset(observed)` inside a torch/jax problem is refused and the backend's own `IndependentNoise` must be passed. Consistent with requiring the backend's `PowerLaw`, and the message names the classes, but it is a papercut on the commonest problem. Softening options (both further §4.5 changes): a neutral `BACKEND` sentinel that `declared_capabilities` treats as agreeing with anything, or a backend-supplied default likelihood on `Dataset`. Fable's recommendation: leave as ruled until W2.10 shows whether it bites in practice.

**Next session**: (1) **Peter's rulings** below — the realisation surface (W2.13) is the one that unblocks NUTS on torch and VI on both, and both tracks' slice 2 should start only after it; the capability-parts widening folds into the same item. (2) Dispatch **W2.13** (Fable drafts the §4.5 contract text; Opus implements — it touches `ampere.core`, both backends and `ampere.inference`, so it is one agent, not two). (3) Then **W2.4 slice 2** and **W2.5 slice 2** in parallel (QuasisepGP library choice by measurement — torch: GPyTorch vs celerite2-torch; jax: tinygp vs celerite2.jax, noting celerite2.jax flips x64 on import — VI, batching/GPU, the native path widened to censoring/latent GPs/non-Gaussian families on jax, native kernels and noise models on torch). (4) **W2.7, W2.8, W2.9** are dispatchable now against merged master (Peter deferred them until the backends were ready — they are; prompts drafted in this session's scratchpad, re-derive from WORK_ITEMS.md if lost: W2.7 owns `plot_residuals`/`plot_gp_localisation`/`plot_anomaly_score`, W2.8 owns `plot_corner`/`plot_trace`/`plot_posterior_predictive` — both touch `ampere/results/plots.py`, so dispatch sequentially or with that split stated). (5) Spec corrections owed from the tracks, small enough to batch into W2.13's PR with decision-log rows: `lowering.md` §3.6's numpyro icdf table (Gamma/Beta need TensorFlow Probability); §4's `biject_to(lowered.support)` advice (torch's `TransformedDistribution.support` is the last transform's codomain — declare the true support); a note that numpyro's `AffineTransform` defaults `domain` to `real`.

**Open ruling for Peter (from W2.4, the larger one) — where the differentiable evaluation path lives.** The frozen spec never defined how a native backend evaluates the *composed* problem differentiably. Verified facts: (1) every `ampere.core` container coerces its values with `np.asarray` (`results_schema.py`), so `FittingProblem.log_prob` cuts a torch autograd graph at the first container — a tensor that requires grad cannot even be put into a `Spectrum`; (2) `Kernel.matrix` is numpy, so a native `DenseGP` receives a numpy covariance and no gradient reaches the GP hyperparameters (W2.4's torch `DenseGP` agrees with the oracle to 1e-9 but is differentiable only in the residual); (3) `ampere.inference` is AST-tested to import nothing but `ampere.core` and `ampere.results`, so a NUTS driver there cannot reach a backend's tensor path by import. `lowering.md` §0 says lowering happens "once, when a `FittingProblem` is *realised* on a backend" and that it "translates declarations, not user code" — the realisation step is named but never specified. W2.4 built every native piece with a tensor twin (`Model.evaluate_tensor`, `Transformation.apply_tensor`, `DenseGP.log_marginal_likelihood_tensor`, `TorchParameterSpace.log_prior_unconstrained_tensor`) and stopped there rather than ship a NUTS driver that finite-differences or differentiates only the prior. **Recommended ruling (Fable): a `Realisation` surface in §4.5, registered like lowering.** `ampere.core` gains `register_realisation(backend, factory)` / `realise(problem)` beside `register_lowering`; a backend registers at import a factory that turns a composed `FittingProblem` into a backend-native object exposing `log_prob_unconstrained(y_native) -> native scalar` (differentiable), `free_size`, `constrain`, and the failure flag for §4.5's non-strict path — built from the backend's own model/step/noise/solver twins and its lowered parameter space, **never** from core containers in the hot loop (which also removes the reference hot loop's 28% container-rebuild cost when the reference backend registers a realisation of its own). `ampere.inference`'s NUTS/VI drivers call `realise(problem)` through core and dispatch on `problem.backend`, keeping the import-graph test intact; the numpy `log_prob` remains the oracle the realisation must agree with (a conformance row). Kernels get the same treatment as models: a native kernel evaluation (`Kernel.matrix_native` or a per-backend kernel class registered by family) so GP hyperparameters are trainable. This is a §4.5 (+§4.4 for kernels) addition with a decision-log row, sized as one item (**W2.13**) that both tracks then consume in slice 2; it blocks NUTS/VI on both backends. The alternatives W2.4 named — a container that can hold a backend array, or narrowing what `differentiable` promises — are recorded and not recommended: the first drags torch/jax into the core container type, the second gives up the ladder's rung-2 promise.

**Smaller rulings from W2.4** — all four landed with W2.13 (shared `LoweringFallbackWarning`, public `evaluation_order()`, `provenance_config`, CPU-only CI ruled).

**Open ruling for Peter (from W2.12) — the capability flags stop at the instrument.** `Dataset.capability_parts` is `tuple(self.instrument.steps)` (`ampere/core/dataset.py:915`), and `FittingProblem` adds the compiled models; that is the whole set `declared_capabilities` reads. `NoiseModel`, `LikelihoodFamily`, `Likelihood` and `GPSolver` (`ampere/core/likelihood.py`) carry **none** of the four flags today — the reference noise models declare `BACKEND` (W2.12) but nothing reads it — so a problem whose models and steps are torch and whose GP solve is the numpy `DenseGP` reports `differentiable=True, backend="torch"` and both are false of the likelihood. The same hole exists for `DEVICE`. Two options; the first is recommended:
1. **Widen at the `GPSolver`/`NoiseModel` level (recommended).** Add the four ClassVars to `NoiseModel` and `GPSolver` with the reference defaults (`False`, `False`, `"cpu"`, `"reference"`); `Likelihood.capability_parts` returns its noise model and, when a GP is declared, its solver; `Dataset.capability_parts` appends `self.likelihood.capability_parts`. `LikelihoodFamily` stays flag-free — the families are closed-form algebra the backend evaluates, not objects that own arrays. Consequences: the reference `DenseGP`/`QuasisepGP` are declared non-differentiable, so a torch problem with a numpy solver correctly loses `differentiable` (which is what NUTS must see) and *fails composition* on `backend` disagreement unless the caller overrides — the loud outcome. A §4.4 + §4.5 change: `likelihoods.md` §7 and `inference.md`'s capability-flags section, one decision-log row, conformance rows for the identity check extended to the likelihood parts. Best landed as a small item **after** W2.4/W2.5 slice 1 merges (their native `DenseGP` classes then declare themselves) and **before** their slice 2.
2. **Leave the parts set as is** and document that the flags describe the model-and-instrument chain only; engines check the GP path separately. Cheaper, but it keeps `ampere_backend` a partial fact, which W2.12 was meant to end.

**State**: master has **W2.12 merged** (`53973cc`; the backend is §4.5's fourth capability flag, schema 4 — Peter ratifies at review), clean, gates green (`pixi run test-all` **1368 passed, zero skips**), lint/format/pyrefly clean. Nothing pushed. Earlier the same day (at 1344 tests) Peter ruled on the two
W2.1 findings and both rulings were **landed directly on master** (the handoff anticipated this
— small enough to skip a branch): (a) `requests` is declared in
`[project.dependencies]` on pyphot 2.1.1's behalf (`a02836d`; the
pyproject comment notes ampere makes no HTTP requests itself and that
httpx/aiohttp/niquests would be better if it ever does — drop the line
when pyphot fixes its metadata; plan pins-row addendum); (b)
`Instrument.__init__` now calls `configure_from` **last step first**
(`8a9534a`; `transformations.md` §5 clarified, decision-log row added
per ground rule 9, regression tests in `tests/core` and
`tests/backends`, both verified to fail under the old forward walk). The
minimal-install CI job should pass at the next push. Peter deferred
branch triage and the push to origin — both stay on the list.

**Next session**: dispatch the backend tracks **W2.4 (torch)** and
**W2.5 (jax)** — Opus, one agent per track, lockstep via the conformance
suite, per `docs/orchestration.md`. Settle the small **backend identity on
`FittingProblem`** addition first or as the tracks' opening step (drivers
currently *declare* `Engine(..., backend=)`; W2.2's finding). Fold into
the dispatch prompts: celerite2 returns quiet NaN where a Cholesky
raises — every backend's celerite2 path must guard preconditions as core
now does; the `eqx.partition`/x64-guard/icdf rulings in plan §2; a
backend that overrides chain construction must keep the last-step-first
`configure_from` order. Dispatchable in parallel to Sonnet: W2.7/W2.8
(follow W2.2), W2.9 (follows W2.1+W2.2). Smaller queued items: the zeus
re-ball helper, the `compile_for` template perf gap in the reference hot
loop, the `examples/` pyphot touch-up.

**Owed and pending**: the **retroactive terra passes** on the four merged
2026-09-05 ranges when the Codex quota returns (~2026-09-30) — W2.3's GP
algebra first — before the M2 milestone; branch-triage approval
(`docs/design/harvest/branch_triage.md`, deferred by Peter 2026-09-07); a
live-CI verification push (deferred likewise); pruning the stale
`.claude/worktrees/` checkouts and merged branches once Peter confirms
nothing on them is wanted.

## ⚡ SUPERSEDED session record (end of the 2026-09-05 session — W2.1/W2.2/W2.3/W2.6 MERGED)

**State**: the first Phase-2 block is **merged and integration-verified**.
Master's merge sequence: W2.1 (`0d8b9a0`), then W2.6 (`152b1fa`), W2.2
(`cea9413`), W2.3 (`c7eec38`) — each Fable-reviewed pre-merge (records in
the merge messages and the status table), with the pyproject dependency
union resolved and `pixi.lock` regenerated once at the W2.3 merge. Gates
on merged master: **`pixi run test-all` 1342 passed, zero skips** (five
suites — core, results, conformance, backends, inference — in one
process), lint/format/pyrefly clean. What now exists: the reference
backend (detector-aware photometry included), `Axis.locate` +
`describe()`/neutral model identity (provenance schema 3), the lowering
registry, emcee/dynesty/zeus drivers emitting DataTrees (arviz + h5netcdf
+ celerite2 are base deps now; the `arviz` extra is gone), and the exact
O(N) `QuasisepGP` (~57 ms at 10⁵ points; `conditional_loo` deferred with
the O(N) route recorded). Nothing pushed to origin.

**Next session**: dispatch the backend tracks **W2.4 (torch)** and
**W2.5 (jax)** — Opus, one agent per track, lockstep via the conformance
suite, per `docs/orchestration.md`. Fold into their dispatch prompts the
carried findings (all in the W2.2/W2.3 status rows): a backend identity
on `FittingProblem` — **done, W2.12**; celerite2 returns quiet NaN where a
Cholesky raises — every backend's celerite2 path must guard preconditions
as core now does; the `eqx.partition`/x64-guard/icdf rulings are already
in the plan §2. Smaller queued items: W2.7/W2.8 (follow W2.2), W2.9
(follows W2.1+W2.2), the zeus re-ball helper, the `compile_for` template
perf gap in the reference hot loop, the `examples/` pyphot touch-up.

**Owed and pending** (as recorded then; the two rulings below were given and landed 2026-09-07 — see the section above): the **retroactive terra passes** on all four merged
ranges when the Codex quota returns (~2026-09-30) — W2.3's GP algebra
first; **Peter's two rulings** — (a) pyphot 2.1.1's undeclared `requests`
(one line: add `requests` to deps or pin past it; the minimal-install CI
job fails at the next push until then) and (b) the chained-LSF
`configure_from` forward-order under-padding (`ampere/core/transform.py:944`;
proposed: reverse the call order + a regression test); branch-triage
approval (`docs/design/harvest/branch_triage.md`); a live-CI verification
push.

## ⚡ SUPERSEDED session record (2026-09-05, mid-session — kept as history)

**State**: **W2.1 is merged** — master is at `85976da` (merge commit
`0d8b9a0`; Opus-authored, Fable-reviewed twice per Peter's 2026-09-05
ruling, standing in for the cross-model pass while the **Codex quota is
exhausted until ~2026-09-30** — a retroactive `gpt-5.6-terra` review of
the merged W2.1 range is **owed** when it returns, before the M2
milestone). `ampere.backends.reference` is live and is the conformance
battery's `reference` fixture (adapter onto the shipped package;
Peter-confirmed, no test-local duplicate). Peter's same-day photometry
ruling is implemented: detector-aware weighting, `detector=` required,
photon `R dλ/λ` / energy `R dλ/λ²`, `from_library` reads pyphot's
per-filter `dtype`. Gates at merge: `pixi run test-all` 1133+85 in one
process, lint/format/pyrefly clean. Nothing pushed to origin.

**Nothing in flight — the whole W2.2/W2.3/W2.6 block is reviewed and
MERGE-READY pending Peter.** Proposed merge order: **W2.6 → W2.2 → W2.3**
(W2.6 touches no dependencies; W2.2 and W2.3 each add base deps in a
dedicated commit — resolve the trivial `[project.dependencies]` overlap,
regenerate `pixi.lock` once after the last merge, then re-run the full
gates on merged master before updating the status table).

**W2.3 (QuasisepGP via celerite2, Opus) is MERGE-READY pending Peter** —
`w2.3-quasisep-reference` at `7ea71bf` (3 commits). Fable-reviewed with
the GP algebra hand-checked line by line (in lieu of terra): the shipped
rank-2 semiseparable Matérn-3/2 representation is an algebraic identity
(verified on paper and by the test reconstructing the dense kernel matrix
to 1.7e-16); celerite2's own `Matern32Term` was rightly rejected (ε²-limit
approximation, misses `tolerances.cross_solver` at its default). Gates
independently re-verified: test-all 1236 with ZERO skips (the two
long-standing skips were the W1.10 solver-agreement skeleton, now live
and passing), conformance 294, lint/format/pyrefly clean; scaling
reproduced (57 ms at 10⁵ points, exponent ~0.94, >1500× over dense).
Review dispositions: placement in `ampere.core` accepted (decision-log
row + honest architecture.md clarifications); the conformance restructure
accepted (the new "not the dense solver in disguise" row is stronger than
the refusal it replaced, whose discipline moved to `tests/core`);
`conditional_loo` deferral properly recorded with the O(N) route named
for a future item; `tests/scaling/` + `pixi run scaling` outside the
gates, pending W2.11's benchmark-harness choice. Out-of-scope finding
worth carrying to W2.4/W2.5: celerite2 returns quiet NaN where a Cholesky
raises — the core solver now guards its preconditions, and the torch/jax
celerite2 paths must do the same.

**W2.2 (engine drivers, Opus) is MERGE-READY pending Peter** —
`w2.2-engine-drivers` at `1b71c10` (4 commits). Fable-reviewed, no fixes
needed; gates independently re-verified (test-all 1303/2 skipped — now
core+results+conformance+backends+inference in one process — lint/format/
pyrefly clean). Emcee/dynesty/zeus in `ampere.inference` against §4.5
only (backend neutrality is import-graph-tested); every run emits the
DataTree with per-dataset log-likelihood decomposition; the R1 promotion
landed with its decision-log row (arviz + h5netcdf base, `arviz` extra
REMOVED — loud failure preferred over an empty alias; pixi feature
renamed `netcdf`). Review dispositions for the PR description: EngineError
stays in `ampere.inference` (R5 moved ResultsError because the contract
text names it; no contract names EngineError — re-exported, movable in
two lines if ruled otherwise); the extras removal, matrix growth and
~4-minute suite all accepted. Out-of-scope findings worth carrying:
zeus draws from TWO process-global RNGs (numpy legacy + stdlib `random`)
— both seeded/restored by the driver, inherited by any future zeus/Pool
work; zeus needs a burn-in/re-ball helper (small future item); the
reference hot loop pays ~28% of 4 ms/eval rebuilding Spectrum containers
when no requirement was published (cheapest visible perf win, W2.1's
code); no FittingProblem surface names its backend, so `Engine(...,
backend=)` is declared — wants a small addition before W2.4/W2.5.
W2.2 and W2.3 both touch pyproject/pixi.lock — merge **sequentially**
(W2.6 → W2.2 → W2.3 natural), regenerating the lockfile at the collision.
All three ranges join the retroactive terra queue.

**W2.6 (lowering registry, Sonnet) is MERGE-READY pending Peter** —
`w2.6-lowering-registry` at `88b16c5` (4 commits; `ampere/core/lowering.py`,
20 tests). Fable-reviewed; two findings were fixed on the branch and
independently re-verified (test-all 1239/2 skipped, lint/format/pyrefly
clean): (1) `tests/core/conftest.py`'s autouse snapshot now restores
`ampere.core.lowering._REGISTRY` alongside `_FAMILIES` — the module-global
registry leaked across the single-process suite, W0.10 finding (c)'s class
again; (2) bijection rows are keyed on the module-qualified class name
(bare `cls.__name__` let two same-named classes silently resolve to each
other's lowering on lookup), bare names kept in messages. Accepted at
review, to record in the PR description: the shared `(kind, name, backend)`
store behind both slots, and provenance stamping via `provenance_attrs`'s
existing `extra=` hook — the first-class schema key is deferred until a
real backend (W2.4/W2.5) drives lowering end-to-end and can populate it
automatically, with the `PROVENANCE_SCHEMA_VERSION` bump taken then.
Merges are held for Peter's go-ahead (his "please merge" covered W2.1).

**Awaiting Peter's ruling** (both recorded in the W2.1 status row, small
enough to land directly on master between merges): (a) pyphot 2.1.1's
undeclared `requests` dependency breaks the minimal-install CI job — add
`requests` to `[project.dependencies]` or pin past the broken pyphot;
(b) `Instrument.__init__` runs `configure_from` in forward chain order
(`ampere/core/transform.py:944`), under-padding chained same-axis
convolutions — proposed fix is calling it in reverse order plus a
regression test. Also still outstanding from Phase 0–1: branch-triage
approval (`docs/design/harvest/branch_triage.md`); the `examples/`
scripts' standalone `pyphot.unit` usage (Phase 2 follow-up); a live-CI
verification push.

**After these three merge**: W2.4 (torch) and W2.5 (jax) backend tracks
become dispatchable (Opus, one per track, lockstep via the conformance
suite); W2.7/W2.8 follow W2.2; W2.9 follows W2.1+W2.2.

## ⚡ Pick up here (end of the 2026-09-03 session — PHASE 1 COMPLETE) — SUPERSEDED by the section above

**State**: **Phase 1 is done and the spec is frozen.** W1.13 merged
2026-09-03 (`a962b09`); the `spec-v1.0` tag sits on the freeze-content
commit `58daa86`. All gates green at the tip: **1109 tests** in a single
pytest process (`pixi run -e dev test-phase1`), pyrefly 0 errors, ruff
lint/format clean, and CI (W0.10) runs the Phase-1 suites on every PR.
Ground rule 9 is now in force: any §4 contract change needs a
decision-log entry in the same PR. Nothing is pushed to origin
(origin/master is far behind by design) — **a push is now worthwhile**,
to let the W0.10 CI prove itself on a live run.

**Next**: dispatch Phase 2 waves from WORK_ITEMS.md's W2.1–W2.11 per
`docs/orchestration.md` (Opus for the backend tracks W2.4/W2.5; the
reference slice W2.1 first — it carries the freeze-escalation landings:
`describe()` + the backend-neutral model identity, and `Axis.locate`,
each with its decision-log entry). **W0.9 merged 2026-09-03** — Phase 0
is entirely complete and W2.1's synthetic-photometry step is ungated:
new photometry code imports `get_unit` from
`ampere/utils/pyphot_compat.py` (the pyphot ≥2 idiom) from its first
line. Known follow-up for Phase 2: the `examples/` scripts still use the
removed `pyphot.unit` API standalone. Also outstanding: branch-triage
approval (`docs/design/harvest/branch_triage.md`). The paragraphs below
are the Phase 0–1 ledger, kept as history.

**R1–R4 ruled 2026-09-02 and implemented**: design B (nested merge)
ratified; `LikelihoodFamily.sample` landed (overridable, default refuses
specifically); dotted channel names landed (grouped data motivate it, not
just populations); the tie-based hierarchical pattern stands, with the
revisit design recorded in `inference.md` §17.2. **W1.8 and W1.10 are
unblocked.**

**W1.8 and W1.10 both merged 2026-09-02** (see WORK_ITEMS.md); W1.8's
seven ruling requests (R1–R7, `results.md` §15) join Peter's backlog.
**Only W1.13 (spec assembly & freeze) remains in Phase 1.** It owes the
batched sol review (still unavailable on this account — Peter
configures codex first) and resolves the decision backlog below.
**Findings worth W0.x items before the freeze** (all three **fixed by
W0.10, merged 2026-09-03**): (a) CI never runs
`tests/core`, `tests/results` or `tests/conformance` — the whole
Phase-1 suite is unexercised by `ci.yml` (the py3.11 import break that
survived until W1.10's review is the proof it bites; fixed `a8e6f54`);
(b) `pixi run <task>` in the *default* environment cannot import
ampere from a worktree — the tasks should default to `dev` or the
default env pin a supported Python; (c) found 2026-09-03: running
`tests/core` and `tests/conformance` in **one pytest process** fails
`test_the_table_covers_every_registered_family` —
a likelihood-family registration leaks into `list_families()`; each
suite alone passes, which is why nothing had caught it (see (a)).
*Attribution corrected at W0.10's review*: the leaking source is the
`likelihoods.md` worked example (a `"laplace"` family) run via
`test_spec_doctests.py` — `test_likelihood.py`'s own registration sites
use `test_*` names with their own cleanup, or fail before registering.
Fixed by W0.10: an autouse registry-snapshot fixture in
`tests/core/conftest.py` covers every source in the directory, and
ci.yml's `phase1-suites` job runs the three suites in one pytest
process to keep it honest.

**Ruled later on 2026-09-02, all landed**: lossless nesting +
`Binding.index` (implemented by the orchestrator after the Opus dispatch
hit its session limit — `parameters.md` §8; the problem's bindings now
compose to the leaves and plate elements route per component);
I-1/I-5/X-2/X-3 approved and landed; X-1 under iteration — Peter is
positive and asked for the full mechanics, now written as the "Detailed
design" subsection of `awkward_instrument.md` §6 X-1.

**Ruled 2026-09-03** (recorded in the spec preambles; implementation
routed to W1.13): W1.6 §17 Q3/Q4/Q6 accepted as recommended (Rice on the
complex value with an `Amplitude` step; `κ = 1/σ²` per sample; circular
complex GP declared `ANALYTIC`), Q5 ruled (the per-dataset scalar *is*
the per-instrument scale; cross-instrument relations via ties or
hierarchical priors, never a per-channel array), Q7 closed with W1.4 Q5,
Q8 superseded by a **consolidated cross-contract serialisation review**
at W1.13 (with `results.md` §15 R7 and `results_schema.md` §17 Q6);
W1.5 §15 Q2/Q3 accepted (`configure_from`, push-forward-and-raise;
`compile_for` refuses by default with an opt-out warning flag), Q5
accepted (`freeze()` — also closes `inference.md` §19 item 8), Q4 ruled
in substance (label and channel are distinct concepts, not one-to-one;
one residual: may the label still *default* to the channel name — see the
§15 preamble); W1.4 §17 Q3–Q6 ruled (`x`/`y` plus a future WCS carrier;
`(x, y, spectral)` stands; (u, v, λ) visibilities confirm the Phase 4
`Axis` extension for `extra_coords`; serialisation folded into the
consolidated review); lowering.md §12 items 4 and 6 ruled
(`LoweringError` into `ampere.core.exceptions`; the reference `icdf`
fallback sanctioned with a loud warning and a `strict` raise). **New
scope recorded** in `transformations.md` §15: image/1-D-spatial
standard-library transformations (PSF convolution, spatial resampling,
affine + WCS, Hankel, NUFFT) — W1.13 confirms the freeze precludes none.
**And later the same day, `results.md` §15 R1–R7 were all accepted as
recommended** (decision-log row in the plan's §2): R3's wording
corrections and R5's `ResultsError`-to-core move were implemented at
once; R2 (`pointwise_log_prob`, granted-not-stored) and R4
(`AnomalyScore` in core) land at W1.13; R1 (arviz + netCDF engine into
the base install) waits for Phase 2's engine drivers; R6 confirmed; R7
folds into the consolidated serialisation review.

**Third ruling batch, later on 2026-09-03 — everything but X-1.**
`inference.md` §19 items 5–7 ruled (substream ratified in core;
capability flags promoted into W1.5's ABCs at the freeze; narrow catch
set confirmed); `lowering.md` §12 item 1 ruled — option (a), **with the
door left ajar** for eventual discrete support ((variational) EM,
enumeration, SBI, nested sampling, Bayesian optimisation): the refusal
lives only in `default_bijection_for`, is a typed capability error, and
discreteness becomes queryable — and item 8 ruled (the hardened
registry, inert to jit/vmap/grad because it resolves before tracing);
`transformations.md` §15 Q4's residual settled (distinct labels required
only when >1 instrument reads a channel, checked at problem composition;
singleton default kept) and Q6 formally deferred, sketch-first; W1.10's
narrowed hash claim ratified (backend-neutral model identity routed to
the freeze beside `describe()`); the serialisation-review scope stands
as proposed; **W0.10 authorised** (CI runs the Phase-1 suites, with
different schedules per part if needed; the registry-leak and
pixi-default-env fixes ride along; py3.11 droppable if problematic).
Decision-log rows in the plan's §2; every ruling is in its spec's
§-preamble; implementation routes to W1.13 (and W0.10).

**X-1 accepted, 2026-09-03 — the review backlog is clear.** Peter
accepted `awkward_instrument.md` §6's detailed design as written; it
lands at the freeze (decision-log row updated; the sketch's §9 Q1 and Q2
close with it, Q4 dissolves). **W1.13 is now dispatchable**: its backlog
is the accepted-recommendation implementations recorded in the spec
preambles, the consolidated serialisation review, X-1, the standard
freeze duties (cross-review, tag, Phase 2 breakdown), and the Fable +
`gpt-5.6-terra` adversarial pass. `awkward_instrument.md` §9 Q3 was
ruled 2026-09-03 (decision-log row): ampere does not ship WStat; the
docs opinionatedly compare the profiled user-family workaround against
the two-dataset Bayesian formulation (pros, cons and results), landing
with Phase 2's engine drivers — W1.13 carries the obligation into the
Phase 2 breakdown. **The review backlog is now empty.** Housekeeping:
W0.9 scheduling; branch-triage approval
(`docs/design/harvest/branch_triage.md`).

**W0.10 dispatched, Fable-reviewed and merged 2026-09-03** (see the
status table): CI now runs `tests/core` and `tests/conformance` on every
matrix leg and all three Phase-1 suites in one pytest process in the
`dev` environment; the registry leak and the broken implicit-default
pixi environment are fixed. Plain `pixi run <task>` now behaves like
`-e dev`. **W1.13 is the only remaining Phase 1 item**, and the suites
it freezes are guarded from here on.

**Restart prompt for a new session** (paste as the opening message):

> Read docs/development.md's "Pick up here" section, then propose the
> first Phase 2 dispatch wave from WORK_ITEMS.md (W2.x) per
> docs/orchestration.md — W0.9 and W2.1 first — with the stale-worktree
> recovery and the obligations each item names. Review each at full
> Fable depth; merging is mine. Ground rule 9 is in force: any §4
> contract change carries a decision-log entry in the same PR.

### Consolidated obligations for the W1.7 dispatch prompt

W1.7 is the item most other contracts have been leaving notes for. Its
prompt must carry, at minimum:

- **Nested merging is an open design question Peter expressly kept open**
  (`parameters.md` §14 preamble): evaluate a nested `ParameterMapping` as
  a first-class option for `DatasetCollection`, on its merits — not
  dismissed because single-call merge is what exists. Ties crossing merge
  levels (`transformations.md` §14) and **nested result channels via
  qualified flat names** (`results_schema.md` §17 preamble — symmetrical
  question, routed here) are part of the same decision.
- `parameters.md` §13: build the joint space on
  `ParameterMapping`/`distribute()`; `merge` is not associative — one
  call, or the nested design above.
- `transformations.md` §14: a `Dataset` pairs an observed container with
  an `Instrument`; the instrument's `label` is its merge component; W1.7
  owns *when* `negotiate` and `compile_for` are called.
- `likelihoods.md` §16 (all load-bearing): call
  `Likelihood.check_alignment(predicted_template, observed)` once at
  construction — it is where unimplemented latent combinations are
  refused, so it is not optional; call `check_engine(..., observed=...)`
  so a censored sample the mask excludes is not counted against a
  gradient-free engine; `Likelihood.parameters` is one flat component for
  the single merge, and `latent_declaration(n)` (n = *retained* samples)
  joins that same merge; convert this contract's `LikelihoodError` on a
  non-positive-definite covariance into §4.5's −inf-with-recorded-reason
  failure signalling (Peter's pending W1.6 §17 Q1 ruling may add a
  `strict=False` alternative).
- Plan §4.5 verbatim: `log_prob`/`log_likelihood`/`log_prior` split,
  `prior_transform`, `simulate`, capability flags, failure signalling,
  RNG policy — `lowering.md` §9.2's `substream(seed, label)` design is
  the seed-derivation policy to adopt (its home in core is W1.9 §12.7's
  ratification item).
- Accept: a toy two-dataset joint problem with a tied parameter,
  end-to-end against a stub model.

### For the W1.11 dispatch prompt(s)

Sketches (a)–(f) per WORK_ITEMS.md, plus the specific claims the merged
specs ask it to check: the deliberately-awkward instrument stress test of
the §4.3/§4.4 split (`prior_art.md` 3M1/Tension 3); whether X-ray RMF/ARF
works as a matrix multiply on an energy-axis `Spectrum` *without a
per-sample exposure concept* (`likelihoods.md` §16); what Fourier sampling
can publish as requirements given no pull-back through chains
(`transformations.md` §13.1/§14); whether closure phases need
`VonMisesFamily` before the freeze, and input to Rice's parameterisation
(`likelihoods.md` §17 Q3/Q4).

**Sol review batched at the W1.13 freeze** (Peter's call): one adversarial
pass over the §4 contract code and the lowering spec. `codex exec -m
gpt-5.6-sol` is currently refused on this account ("not supported when
using Codex with a ChatGPT account"); Peter expects to configure the extra
codex steps himself before the freeze. Nothing blocks on it until then.
**Ruled 2026-09-03**: sol remains blocked, so the batched pass proceeds
with a combination of adversarial reviews from **Fable** and from
**`gpt-5.6-terra`** (via the codex CLI), chosen to capture the difference
in answer distributions between model families; sol can still be added if
it becomes available before the freeze.

## Phase 1 review outcomes (2026-09-01, all merged)

- **W1.1** (`w1.1-prior-art-memo`, Sonnet): merged as authored — no
  changes. Its two most load-bearing claims were re-verified against live
  sources during review: bilby's current core really does take an explicit
  `log_likelihood(parameters)` argument (checked against `bilby-dev/bilby`
  `main`), and gammapy issue #2859 / fix PR #2861 are as described. The
  memo's own §7 confidence table is accurate; its weakest section (3ML,
  tutorial-level sourcing) is honestly flagged there.
- **W1.2** (`w1.2-architecture-spec`, Sonnet + Fable amendments): one
  amendment commit on top of the draft. Core's dependency floor corrected
  to numpy/scipy/**astropy**/stdlib (the plan requires units on parameters;
  astropy is a required base dependency; this resolves W1.3's open
  question 1). §10 rewritten as the W1.1 reconciliation record — none of
  the four questions it posed required structural change. One genuine
  conflict with the plan was found (curated astropy→native translation:
  opt-in vs the plan's "silently restoring differentiability") and **ruled
  by Peter on 2026-09-01: opt-in, never silent** — recorded in
  `DEVELOPMENT_PLAN.md` §2's decision table with §4.7 amended, and
  architecture.md §1 updated to match (sketch: a backend-scoped
  `from_astropy()` constructor as the explicit opt-in, raising rather than
  silently falling back when a model is not fully translatable).
- **W1.3** (`w1.3-parameter-contract`, Opus + Fable fix commit): strong
  contract; three defects found and fixed on the branch, with tests
  (112 pass; pyrefly/ruff/format clean): (1) tying two hierarchical-prior
  parameters whose hyperparameters were not themselves tied silently wired
  the collapsed parameter to the first component's hyperparameters — now a
  `TyingError` with a hint; (2) `describe_prior` kept scipy's positional/
  keyword split, so `st.norm(1.5, 0.1)` and `st.norm(loc=1.5, scale=0.1)`
  counted as different priors when tied and serialised differently — the
  description is now canonical (keyword-only), which also gives W1.9 a
  single form to lower; (3) conflicting explicit bijections on tied sites
  silently took the first — now a `TyingError`. Checked against W1.1
  retrospectively: its declaration-based tying is exactly what memo lesson
  G2 prescribes.

- **W1.4** (`w1.4-results-schema`, Opus + Fable fix commit): merged
  2026-09-01. Strong contract; the agent's own self-review caught three
  unit-handling bugs, and the Fable pass found and fixed two more:
  (1) `from_unsorted` permuted the axis/values/uncertainty/mask but left
  `extra_coords` in declaration order — silent per-sample misalignment,
  now permuted with the data; (2) `to_unit` leaked astropy's raw
  `UnitConversionError` where every sibling path raises `SchemaError` —
  now wrapped. 247 core tests; W1.3 untouched.

- **W1.5** (`w1.5-transformation-contract`, Opus): merged 2026-09-01 with
  **no defects needing Fable fixes** — the agent's own self-review caught
  the deep one (a requirements union imposing one instrument's fine
  sampling on another's broad coverage; fixed with per-interval
  densities). Mask propagation is *enforced* (a step dropping its input's
  mask raises); the `Model` ABC landed here minimally with the
  engine-facing surface reserved for W1.7. Review probes confirmed loud
  partial-value failures, correct cross-unit unions, and negotiated grids
  passing `Spectrum`'s ordering validation.
- **W1.6** (`w1.6-likelihood-contract`, Opus): merged 2026-09-01, the
  strongest deliverable of the phase; no Fable fixes needed. The critical
  defect (Student-t/Cauchy/complex-Gaussian silently discarding the GP
  while every check reported a working latent problem) was caught by the
  adversarial review the agent commissioned itself and fixed with the
  `CONSUMES_LATENT_GP` opt-in (default `False`, so third-party families
  inherit the refusal). DenseGP anchor validated to ~1e-15; Fable
  verification covered the marginal-likelihood/conditioning/whitening
  mathematics, an independent Tobit-censoring probe (exact), and the
  masking-beats-censoring declaration behaviour. Key positions on record:
  masks are consumed via `weights()` + row/column excision (the
  infinite-variance limit diverges — spec §8 has the proof); the latent
  declaration is the whitened non-centred form, *not* `HierarchicalPrior`
  (wrong at any N — GP priors are not i.i.d.); the plan's "per-sample
  log_likelihood" ambiguity (per-draw vs per-observation) is flagged for
  W1.13, since a GP likelihood has no per-observation decomposition.

**Rulings recorded this session** (all implemented, not just noted):

- **W1.9**: *both rulings approved by Peter and merged 2026-09-01* —
  `eqx.partition` over `paramax.NonTrainable`, and x64 guard-and-raise —
  recorded in the plan's §2 decision table, §4.1 and `architecture.md` §5
  amended to match. Three follow-ups from Peter's review, all recorded:
  dependent priors (y ~ norm(x, 1)) are already first-class via
  `HierarchicalPrior` (verified live against the merged contract —
  `parameters.md` §9); `numpyro.enable_x64()` was source-verified as a
  thin wrapper over the same jax flag, so the construction guard covers
  numpyro with no separate switch (lowering.md §12.5); and
  backend-specific lowerings for user-defined priors/bijections are now
  lowering.md §12 item 8 (registration-hook sketch, W1.13 decides the
  principle, Phase 2 builds the plumbing). Remaining §12 items route to
  W1.13.
- **W1.12**: *approved by Peter and merged 2026-09-01* — the §10 defaults
  (peer `diagnostics` namespace, `AnomalyScore` in core, RHMF expert
  opt-in) stand as written. His review added one thing, recorded as the
  spec's §11: posterior calibration (SBC, expected-coverage/TARP) as a
  fourth, future diagnostic family — it validates the *inference* rather
  than the model, matters most for Phase 3's SBI layer, and needs no new
  contract surface (§4.5's `simulate()` is the hook). §11 also records
  the extension template for future families, since the area will evolve.
- W1.3 spec §14: *ruled by Peter 2026-09-01* — recommendations stand (flat
  tie labels; lone `shared_as` allowed; `npars` dropped;
  `OptionalDependencyError` to W1.13), **except that recursive merge is
  expressly kept open**: nested/hierarchical merging may be the natural
  approach for `DatasetCollection`, so W1.7 must evaluate a nested
  `ParameterMapping` as a first-class design option (recorded in
  `parameters.md` §14's preamble; W1.5 §14's ties-across-levels note is
  part of the same question). **W1.7's dispatch prompt must carry this.**
- W1.4 spec §17: *ruled by Peter 2026-09-01* — `"default"` stays
  unreserved (clash accepted, loudly documented); overlapping échelle
  orders as two channels confirmed; and θ-carrying results implemented:
  `ModelResult.parameters` exists and `Model.__call__` attaches the
  resolved values automatically (`65c997f`). His follow-up — nested
  result channels — assessed as feasible via qualified flat names
  (mirroring `ParameterSet.merge`) and **routed to W1.7 together with the
  symmetrical nested-merge question**; see `results_schema.md` §17's
  preamble. Questions 3–6 (`Cube` axis order, `extra_coords` units, axis
  naming, serialisation) remain open as written.
**Peter's remaining reading list** (nothing blocks dispatches; all are
freeze-relevant):

- W1.5 spec §15 — **three of eight ruled 2026-09-02** (recorded in the
  §15 preamble): the ANY mask rule stands as the default; `Model` stays
  in W1.5; the `max_step`/`min_resolving_power` scalars are demoted to
  constructor-only, with `segments()` the only public density form
  (implemented and merged the same day). Still open: Q2/Q3 now carry the
  interferometry sketch's recommendations (`configure_from` over
  `pull_back`; `compile_for` may raise — gaps I-3/I-4) and await ruling;
  Q4 (label collisions) has its concrete failure mode recorded in the
  spectrum+photometry sketch and is routed to W1.7; Q5/Q6 as written.
- W1.6 spec §17 — **the two sharpest ruled 2026-09-02** (recorded in the
  §17 preamble): Q1, a `strict` toggle — non-strict engine path
  (−inf + recorded reason, forced by jax tracing), strict raise for
  debugging, W1.7 owns conversion and recording; Q2, the `Dataset`
  resolves the effective mask once at construction, making the mask an
  evaluation-time invariant. Still awaiting ruling with W1.11
  recommendations on file: Rice's parameterisation and von Mises's
  concentration (Q3/Q4), the circular complex GP as `ANALYTIC` (Q6),
  `extra_coords` units eventually-not-freeze (Q7). Open as written:
  per-dataset vs per-channel `IndependentNoise.scale` (Q5);
  `Likelihood.to_spec()` for provenance (Q8, with W1.8).
- W1.4 spec §17 questions 3–6 (`Cube` axis order vs FITS — cheap to
  change only until the freeze; axis naming; `extra_coords` units;
  container serialisation, owed to W1.8).
- W1.9/`lowering.md` §12's remaining ratification items, all routed to
  W1.13 (discrete-family default bijection; `LoweringError` placement;
  reference-path `icdf` fallback; `substream` in core; backend-specific
  lowering registration — item 8, from Peter's review).

**Housekeeping for the next session:**

- **Nothing is pushed to origin.** All session work is on local `master`;
  pushing is Peter's call and Peter's action (agents never push).
- **Leftover worktrees and scratch branches**: several agent worktrees
  remain under `.claude/worktrees/` and the merged item branches
  (`w1.1-*` … `w1.6-*`, `w1.9-*`, `w1.12-*`) plus stale
  `worktree-agent-*` scratch branches still exist. All item branches are
  fully merged; cleanup (`git worktree remove` + branch deletion) awaits
  Peter's explicit approval per ground rule 3. Until then, note that
  `pixi run format-check` in the main checkout reports ~20 files to
  reformat — **all inside `.claude/worktrees/`** (nested copies of
  `tests/characterisation` that the task's top-level exclude does not
  match). Tracked files are format-clean; do not "fix" those.
- **Tooling trap (confirmed twice)**: `ruff format <explicit path>`
  bypasses `extend-exclude` and will rewrite the frozen W0.7 harvest
  snapshots under `docs/`. Format only `ampere/core tests/core` (or use
  the pixi task). Worth a guard in a small W0.x follow-up.
- A `SendFeedback` draft about the `/model` switch-back failure from the
  pre-handoff session remains queued; Peter can review/send it with
  `/feedback`.

## Current state (SUPERSEDED — see "⚡ Pick up here" above; kept as the end-of-2026-09-01 snapshot)

- **Phase 0: complete** (W0.1–W0.8 merged; CI green; issues #74–77
  closed; archival executed). **W0.9 not started** — the temporary
  `pyphot<2`/`sbi<0.28` pins stand, documented in `DEVELOPMENT_PLAN.md`
  §2; the pyphot half must land before Phase 2's synthetic-photometry
  Transformation is written.
- **Phase 1: eight of thirteen items merged** (W1.1–W1.6, W1.9, W1.12);
  W1.7 + W1.11 ready to dispatch in parallel; then W1.8, W1.10, W1.13.
  See "Pick up here" above for the restart procedure and consolidated
  dispatch obligations.
- **Decisions ruled this session** (all in `DEVELOPMENT_PLAN.md` §2's
  table): curated astropy→native translation is opt-in/never silent; jax
  non-trainables lower via `eqx.partition`, not paramax; jax x64 is
  guard-and-raise, never set-on-import.
- The AMPERE paper revision proceeds on the legacy code and takes priority
  in any conflict over `examples/examples_paper/`.
