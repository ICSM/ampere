# Handoff archive — session records superseded by `docs/development.md`'s "Pick up here"

Kept verbatim as history (rulings, review outcomes, dispatch obligations as they were written). The live restart point is always the "⚡ Pick up here" section of `docs/development.md`; nothing here is current.

## Archived 2026-09-08: the 2026-09-07/08 live section as it had accreted

## ⚡ Pick up here (2026-09-07 session — the two W2.1 rulings landed; backend tracks next)

**In flight**: nothing. **W2.4 slice 1 (torch) and W2.5 slice 1 (jax) are MERGED** (`6964b80`, `3458657`; status rows in WORK_ITEMS.md carry the review records and every carried finding). Gates on merged master (`pixi run -e <env> test-all`, five suites in one process): **torch 1746 passed / 10 skipped, jax 1662 / 7, dev 1369 / 5** — every skip is a quasiseparable conformance row awaiting slice 2 or a backend-only test module declining without its extra; lint/format clean; pyrefly 0 errors in all three environments (each backend package checked with its real types in its own environment, excluded where its library is absent). Nothing pushed to origin.

**Rulings given by Peter 2026-09-07 (evening)**: the realisation registry approach is approved and **prototyped on branch `w2.13-realisation-prototype` at `9b26c3f`** (core `ampere/core/realisation.py`: `register_realisation`/`realise`/`registered_realisations`, checked on use against the contract path; jax registers `lower_problem` at import; `NUTSEngine(problem)` needs no density argument; 12 core rows + 4 NUTS rows; dev/jax suites and both typechecks green) — W2.13 turns it into the contract (§4.5 text, decision-log row, a conformance row per registered realisation, the torch realisation, and the fold-ins). `provenance_config` on `GPSolver` approved. CI stays CPU-only for torch/jax; **GPU smoke tests later** (API-works-level only, trusting CPU/GPU library parity). Dispatch order as proposed: W2.13 → both slice 2s in parallel → W2.7/W2.8/W2.9 as capacity allows.

**W2.13, W2.4 slice 2 and W2.5 slice 2 are ALL MERGED** (`0e7c97c`, `a529129`, `dd2c479`); status rows carry the review records. Gates on merged master (the reconciled tip, re-run by Fable): **torch 1924 passed / 19 skipped, jax 1806 / 19, dev 1419 / 69**, lint/format clean, pyrefly 0 errors in all three environments. **In flight (2026-09-08)**: three agents dispatched from `6a0c128` — **W2.14** (`w2.14-latent-gp`, Opus: the latent-GP fix in `GaussianProcessNoise.noise_params`, both backends' latent paths, a conformance shape; ruled by Peter's 'continue' on the proposed plan), **W2.7** (`w2.7-diagnostics`, Opus: families B and C in `ampere.results` owning `plot_residuals`/`plot_gp_localisation`/`plot_anomaly_score`, RHMF re-verified), **W2.9** (`w2.9-wstat-example`, Sonnet: the docs example). **W2.7 MERGED** at `f57c72f` and **W2.9 MERGED** at `196b48f` (status rows; W2.9 found `pixi run docs` broken on master for pre-existing environment reasons — carried to W2.11); **W2.8 dispatched** from it (`w2.8-results-completion`, Opus: corner/trace/posterior-predictive, `add_posterior_predictive`, the pointwise group, the training-set writer). On return: Fable review each, merge in whatever order they land (W2.14 first if concurrent), gates one at a time. Nothing pushed. Both modern backends now have: native QuasisepGP (celerite2 on both, by measurement), the realised path over every core family/censoring/both solvers, NUTS and VI through `realise`, batching, float32/device opt-ins, schema 5 provenance. Memory note for the orchestrator: two five-suite gates running concurrently exhaust this machine's 13 GB — run them one at a time, torch detached (its pyro NUTS rows are slow).

**2026-09-08 midday: W2.7, W2.9, W2.14 and W2.8 are ALL MERGED** (`f57c72f`, `c5b4915`, `4bf1164`, `71b2e96`; status rows). Master gates at `4ec855a`: **dev 1542/71, jax 1933/21, torch 2052/21**, lint/format clean, pyrefly 0 errors ×3. **W2.11 MERGED** at `dbbc172` (status row; the lockfile changed — CPU torch, pandoc, ipykernel, pytest-benchmark — so `pixi install -e <env>` is needed in any stale checkout). Master gates at `205b84f` re-run by Fable: **dev 1555/71, jax 1946/21, torch 2065/21** (CPU wheel confirmed `2.13.0+cpu`), lint/types clean. **W2.10 MERGED — MILESTONE M2 REACHED** at `956e1ab` (status row carries the review dispositions and the numbers; the plan's Matérn row and §5 M2 paragraph record it). **Peter's ruling 2026-09-08 (afternoon): finish slice 3, then a documentation sweep, then Phase 3 in clean sessions.** Items written at `6fd8520` (W2.4 slice 3, W2.5 slice 3, W2.15). **In flight**: `w2.4-torch-slice3` and `w2.5-jax-slice3` (Opus, parallel, from `6fd8520`; ownership split — torch owns its package and fixture, jax owns its package, its fixture and `ampere/inference/engine.py` for the gradient-free fast path; shared conformance files only gain rows at the end). On return: Fable review each, merge sequentially, gates one at a time; then dispatch **W2.15** (docs sweep, Opus) and Fable restructures `docs/development.md` + `CLAUDE.md` for Phase 3. Phase 2's remaining items were the backend tracks' slice 3 (per-instance `device=` with the GPU smoke-test item, `complex_gaussian`, torch/jax `conditional_loo` parity, a jax fast path for the gradient-free engines — the contract path costs ~28 ms flat there) and Peter's open confirmations; then **Phase 3 (SBI)** opens — its items are not yet written in WORK_ITEMS.md (Fable drafts them from plan §5 Phase 3 before dispatch). The retroactive terra queue when the Codex quota returns (~2026-09-30): W2.1, W2.6, W2.2, W2.3 (W2.3's GP algebra first), then W2.13, W2.4/W2.5 slice 2 (the GP transcriptions), W2.14, W2.10. — its scoping note is on the item; both backends have QuasisepGP, the plots and diagnostics exist, and the harness will be chosen. Nothing pushed.

**2026-09-08 afternoon: W2.11 is IMPLEMENTED and awaiting review** on branch `w2.11-ci-phase2` (from `4ec855a`; nothing pushed, nothing merged). What is on it: `ci.yml` grew a `backend-suites` matrix (`torch`, `jax`, `fail-fast: false`, each running that environment's `typecheck` + `test-all` + `bench`), a `docs` job, and a `suites` job on `dev` — one environment per job, never two, on memory grounds as well as isolation; `minimal-install` stays a bare `pip install -e .` and now asserts that torch and jax are genuinely absent, so its import check proves something. `pixi run docs` **builds** (the `dev` pixi feature declares `pandoc` and `ipykernel`; `nbsphinx_execute = 'never'` with a per-notebook note; the dead `ampere.infer.ptemceesearch` autodoc entry is gone) — without `-W`, because the 21 residual warnings are legacy docstrings under frozen paths; that deviation from §5's Phase-1 line is in the decision-log row. `tests/examples` joined `test-all`. **The benchmark harness is pytest-benchmark, not asv** (decision-log row; §6's bullet struck): `tests/benchmarks/`, `pixi run bench`, `benchmark.json` uploaded per environment. Both `torch` and `sbi` now resolve torch from PyTorch's CPU index — all fifteen `nvidia-*` wheels are gone from `pixi.lock`. `tests/scaling` gained the torch rows it was missing. Branch gates: dev 1555/71, jax 1946/21, torch 2065/21; lint, format, `pixi run docs` and all three typechecks clean; the `minimal-install` job reproduced in a fresh venv. **For Peter at the merge review**: the `phase1-suites` job was renamed to `suites`, so any branch-protection required-check entry naming the old job needs updating.

**Next session**: (1) Peter's confirmations — the noise-model papercut (below), W2.7's pooling/`datasets=None` readings, the RHMF deferral row, and W2.11's benchmark-harness and docs-`-W` rows. (2) Review and merge **W2.11**. (3) Then **W2.10** (M2 — the misspecification study; its item notes name W2.11's harness choice as a prerequisite, now taken). Slice 3 of the backend tracks (per-instance `device=` with the GPU smoke-test item — CI is CPU-only by ruling and W2.11 now enforces that with CPU wheels — a torch/jax `conditional_loo` parity decision, `complex_gaussian` in Phase 4) waits behind those.

**The latent-GP science bug is FIXED — W2.14 merged 2026-09-08** (status row): `GaussianProcessNoise.noise_params` applies the whitening transform, both backends' native latent paths are live, and the conformance battery has a latent shape held to the corrected oracle. The earlier ruling text is superseded.

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


## Handoff as of the Phase 2 close (2026-09-08, superseded 2026-09-08 by the Phase 3 draft)

## ⚡ Pick up here — PHASE 2 CLOSED 2026-09-08; Phase 3 begins in a clean session

**Where the project is.** Phases 0–2 are complete. The core contracts are frozen (`spec-v1.0`); both modern backends implement them alongside the numpy reference; milestone M2 is reached (W2.10, the flagship misspecification validation on all three backends, with the science claim asserted as thresholds); slice 3 (devices, complex Gaussian, LOO parity, the gradient-free fast path) and the documentation sweep (W2.15) are merged. **Master is `fbd6970`+**, clean. Gates on the final Phase 2 code (`227d46a`, unchanged by the docs-only W2.15 apart from one spelling of a constant): **dev 1616 passed / 111 skipped, jax 2050 / 42, torch 2222 / 40**, lint/format clean, pyrefly 0 errors in all three environments; the docs build succeeds with 13 legacy-only warnings; the characterisation suite passes. **Nothing is pushed to origin** (Peter's call — the first live CI run is still owed).

**Peter's standing ruling (2026-09-08)**: Phase 3 begins in clean sessions. The procedure: (1) read this section, then WORK_ITEMS.md's status table (every merged item's row is the authoritative record of what landed and what it carried), then `DEVELOPMENT_PLAN.md` §5 Phase 3 and §2's decision log; (2) **Fable drafts the Phase 3 items** into WORK_ITEMS.md — the SBI layer: `simulate`'s batched form and the coordinate–value–mask tensor encoding for embedding networks (`inference.md` §13, limitation 17.5), the `sbi` package integration on torch and a jax-native route only if maturity warrants (plan §6), trained-artefact caching keyed on the spec hash (plan §7's trap), the training-set writer's unlimited-dimension append (`results.md` §13.9) when budgets need it, and SBC/coverage as `diagnostics.md` §11's fourth family — and gets **Peter's approval of the breakdown before dispatch**; (3) `docs/orchestration.md`'s policy is unchanged (Opus for judgement-within-spec, Sonnet for well-specified items, Fable for contracts and reviews; cross-model review via terra when the quota returns); (4) give every shared file exactly one owner when dispatching in parallel, and run five-suite gates one at a time.

**Small housekeeping for the first Phase 3 session** (one Sonnet item or done by hand): `pip install ampere` is an unrelated PyPI package — fix `OptionalDependencyError`'s remedy text (`ampere/core/exceptions.py`) and the six docstrings that quote it (W2.15's finding; the install page already warns); `SAMPLE_STATS_GROUP` is defined twice in `ampere.results`; rename `LSFConvolution.sigma_tensor` (torch) so it cannot be mistaken for the noise hook; add a native zero-uncertainty precondition to the jax realised GP path (the `realise` guard catches it today).

**Open confirmations for Peter** (none blocking; each recorded where it arose): the noise-model papercut from W2.13 fold-in 7 (recommendation: leave as ruled); W2.7's median pooling and `datasets=None` reading; W2.8's corner/trace caps and the exclusion of prior-rejected draws from corner plots; the "Peter to ratify" decision-log rows — RHMF deferral (W2.7), pytest-benchmark and the docs `-W` deferral (W2.11), the Matérn validation note (W2.10), the W2.15 annotation batch; `use_realisation=True` on torch gains nothing (W2.5 slice 3) and `ampere_realised = 1` on gradient-free realised runs is a visible provenance change; `:no-index:` on the legacy API pages; outside the repository, any branch-protection rule naming the old `phase1-suites` job now needs `suites`.

**Owed**: retroactive `gpt-5.6-terra` passes when the Codex quota returns (~2026-09-30), in this order — W2.3 (GP algebra), W2.1, W2.6, W2.2, then W2.13, W2.4/W2.5 slice 2 (the celerite2 transcriptions), W2.14, W2.10, W2.5 slice 3 (the LOO scan); branch-triage approval (`docs/design/harvest/branch_triage.md`); the push to origin and first live CI run; deletion of the auto-generated `worktree-agent-*` branches (harmless; Peter's action); GPU smoke tests actually run on hardware (`pixi run gpu` in the `torch`/`jax` environments on a machine with an accelerator).

**Operations notes for the orchestrator**: this machine has 13 GB — never run two five-suite gates concurrently (the harness kills background tasks on low memory); run torch's `test-all` detached (`nohup … &`) and watch the log with a Monitor, since it takes ~9–11 minutes; after a lockfile change, `pixi install -e <env>` in the main checkout; the merged-master gate is the one that counts, re-run after every merge; parallel agents appending decision-log rows conflict at the table's last line — keep both rows, keep the table contiguous; the two-track pattern (torch/jax) works when each shared file has one owner and the sibling's additions are appended.



## Handoff as of 2026-09-09 mid-Phase 3 (superseded 2026-09-09 by the clean rewrite; incremental record of the W3.0–W3.8 session)

## ⚡ Pick up here — W3.8 MERGED 2026-09-09; W3.6 in flight; merged-master gates queued

**Where the project is.** Phases 0–2 are complete (spec frozen at `spec-v1.0`; three backends; M2 reached; W2.15 docs sweep merged). Master is `09357b4`+ and clean apart from the untracked paper script under `examples/examples_paper/` (leave it alone). Gates on the final Phase 2 code (`227d46a`): **dev 1616 passed / 111 skipped, jax 2050 / 42, torch 2222 / 40**, lint/format clean, pyrefly 0 errors ×3; docs build succeeds with 13 legacy-only warnings; characterisation passes. **Nothing is pushed to origin** (origin/master is still at the W1.3 handoff, `b8e585b`; the first live CI run is owed on the whole of Phase 2). The `sbi` pixi environment was installed 2026-09-08 and resolves **sbi 0.27.0 + torch 2.13 (CPU)**.

**What happened 2026-09-08/09 (Fable).** Phase 3 W3.0–W3.10 drafted, approved and recorded (WORK_ITEMS.md "Phase 3"). **Merged**: W3.0 (`f0140dc`), W3.9 (`fcdd07d`), **W3.1 slice 1 (`11d01d7`)** — merged-master gates dev 1729/111, jax 2176/42, torch 2346/40, lint/format/pyrefly clean. Master is `11d01d7`+. **W3.2 merged** at `7ce3e81` with review fixes at `cbad8b4` (version attrs `str()`-wrapped; the `sbi` environment now includes the `torch` feature — `pixi install -e sbi` after pulling) and its decision-log row at `3d9cc2c`; merged-master gates sbi 2408/41, dev 1748/155, jax 2195/86, torch 2365/84 — all green. **W3.1 slice 2 merged** at `9e722ba` (review commit `3c034b9`); its four merged-master gates (dev, jax, torch, sbi) were queued behind the lock at the merge — results in the status row once recorded; if this session was cut off before then, re-run them one at a time. The encoding contract is drafted (`docs/design/contracts/encoding.md`, `1649a76`) and is W3.3's binding text. **Previously in flight (now merged)**: **W3.1 slice 2** (Opus, branch `w3.1-slice2-native`; owns `ampere/backends/{torch,jax}/*`, `ampere/core/simulate.py`, `ampere/core/likelihood.py`'s sampling dispatch, `inference.md`/`likelihoods.md`, the batched-simulation decision-log row). Gates are serialised through `flock /tmp/ampere-gate.lock` (the orchestrator's queued runs and the agent's alike). If this session was cut off: `git branch --list 'w3.*'`, `git worktree list`; review each branch with commits against its Accept line and merge `--no-ff` (W3.2 first if both are ready — slice 2 must then rebase), re-running merged-master gates one at a time; re-dispatch anything without commits from the item text.

**Approval and the W3.1 revision (Peter, 2026-09-08).** The breakdown is **approved**. Peter's three notes on W3.1 are recorded as a decision-log row and W3.1 is now two slices: slice 1 (numpy path) builds `simulate_many` on an order-preserving, partition-independent **executor protocol** with chunking, makes **external compiled simulators** the first-class case (process pool, timeouts, crash capture), and gives `BATCHABLE` a reference-backend meaning; slice 2 lands per-chunk `vmap` prediction and **native observation sampling on torch and jax** (numpy path the oracle, compared distributionally). Rows 85 and 90 ratified. W3.10 (automatic paging above the plot caps, loud warning; not urgent) added. **Dispatched**: W3.9 (Sonnet, docs-only, branch `w3.9-migration-page`). **Held**: W3.1 slice 1 — Fable asked Peter to glance at the revised text before an Opus session is spent on it; dispatch on his nod, from master, with the instruction to wait for any running gate (`pgrep -f pytest`) before starting its own.

**Horizon notes (2026-09-09)**: Peter's Phase 5+ look-ahead is in `docs/design/horizon_notes.md`; two things came out of it into the present — the **observation-context** reservation in the Phase 3 section (confirmed) and a **kernel-algebra item added to the plan's Phase 4** bullets (ruled).

**Phase 6 note (2026-09-09)**: a multi-observation composition tutorial (photometry + spectra, calibration uncertainty) is recorded in the plan's Phase 6 bullets — capability exists, example does not.

**Ruled by Peter 2026-09-09**: `ProcessExecutor` defaults to `forkserver` on POSIX (`fork` via `mp_context=`) — adopting it now also stops different Python versions seeing different behaviour and errors as the platform default changes; slice 2 implements and records it in the decision-log row.

**Rulings 2026-09-09 on slice 2's questions**: native sampling stays the default; further core `sample` implementations when a use case arrives, with the pathway kept open; the mesh sharder's padding accepted, to be tested and costed with the GPU item. The GPU item's scope now includes: run `tests/gpu` for real, and cost `MeshSharder`'s padding against refusal.

**For Peter — what is still open.** Only the repository actions he alone takes: push master to origin and watch the first live CI run (a branch-protection rule naming the old `phase1-suites` job needs `suites`); delete the stray local `worktree-agent-*` branches; run `pixi run gpu` in the `torch`/`jax` environments on a machine with an accelerator. Every ruling from the 2026-09-08 digest is recorded; W2.9's coverage study is scheduled inside W3.6.

**Owed (unchanged)**: retroactive `gpt-5.6-terra` passes when the Codex quota returns (~2026-09-30), in this order — W2.3 (GP algebra), W2.1, W2.6, W2.2, then W2.13, W2.4/W2.5 slice 2 (the celerite2 transcriptions), W2.14, W2.10, W2.5 slice 3 (the LOO scan); Phase 3 adds W3.1, W3.2, W3.6 to the queue as they merge.

**Resolved**: both agents resumed after the reset and reported; **W3.3 merged** at `1ca558a`, **W3.5 merged** at `6dba1f3` with the wiring commit after it; merged-master gates: dev 1867/191, jax 2338/122 green; **torch 2507 passed with 2 failures on the first run** — one a W3.5 test that assumed torch and sbi come together (fixed at `ed3eaaf`: checked per package), the other the backend-neutrality subprocess crashing while the orchestrator's `ruff format` rewrote `_sbi.py` under it (does not reproduce; the check passes standalone); both re-verified in the torch environment and the full torch re-run is **green: 2509/120**; **sbi 2574/55 green** — all four merged-master gates green at `ed3eaaf`+. Lesson for the operations notes: **never edit or format a tracked file while a gate is running on the main checkout.** **Dispatched 2026-09-09 from `2b3d621`, in parallel on disjoint files**: **W3.6** (Opus, branch `w3.6-calibration`; owns `ampere/results/*` except `provenance.py`, `ampere/inference/_sbi.py`, `tests/results/*`, `tests/inference/test_sbi.py`, the WStat example/page, `diagnostics.md`, `results.md`) and **W3.8** (Opus, branch `w3.8-foreign-parts`; owns `ampere/core/{dataset,likelihood,realisation,exceptions}.py`, `ampere/backends/*`, `tests/core|backends|conformance`, `inference.md`). Both gate through `flock /tmp/ampere-gate.lock`, behind the four merged-master gates queued at the W3.5 merge. **W3.8 merged** (see the status row; its decision-log row replaced with the landed text); four merged-master gates queued at the merge (`<scratchpad>/gates-w38.summary`, else re-run one at a time). W3.6 still in flight — merge it next (dev + sbi gates), then dispatch W3.4, W3.11, W3.12 (all touch `_sbi.py` or `results/`; W3.12 bumps the schema, so it goes alone or last), then W3.7. If cut off: `git branch --list 'w3.*'`, `git worktree list`; review branches with commits against their Accept lines. Then: **W3.4 is ruled (2026-09-09): TMNRE through sbi's own ratio estimators, no swyft** — dispatch after W3.6 merges (both touch `_sbi.py`); **W3.11** (embedding readouts) and **W3.12** (model identity hash, schema 6) are **ruled** (2026-09-09) and dispatch after W3.6 merges; an embedding *study* (width and readout vs data/model/structure, guidance for users) is Peter's rider, recorded in the deferred list; W3.7 (CI leg) last; W3.10 whenever convenient. The interruption record follows for history.

**Interruption (2026-09-09, early morning).** Both agents were killed by the account's session rate limit (resets **08:30 Europe/London**); the orchestrating session may be too. **All their work is committed** — nothing was lost — and the slice-2 merged-master gates finished green (dev 1780/169, jax 2251/100, torch 2422/98, sbi 2465/55), so master `4d6b315`+ is sound. State of the two branches, both from `3c034b9`, both worktrees still present under `.claude/worktrees/`:
- **`w3.3-encoding`** (Opus; worktree `agent-ac96abd2757708af3`): 4 commits — the encoding/layout/unpack, `layout=` on `SBIEngine` with the two masked embeddings, `encoding.md` bound and amended in place, a numpy-warning fix with a degenerate-masked-sample test. **Not yet done**: its `dev` and `sbi` gates and lint/format/typecheck, and the report (with the decision-log row amendment text). The orchestrator queued both gates in the background at the interruption (`<scratchpad>/gates-w33.summary` — the scratchpad is session-local, so if it is gone, re-run: from the worktree, `PYTHONPATH=$PWD ../../.pixi/envs/<env>/bin/python -m pytest tests/core tests/results tests/conformance tests/backends tests/inference tests/examples tests/m2 -q`, one env at a time).
- **`w3.5-artefact-cache`** (Sonnet; worktree `agent-a492de696162058b0`): 4 commits — `ampere/results/artefacts.py` (`ArtefactKey`/`artefact_key`/`ArtefactStore`, with a `model_hash` ingredient added beyond the item text so a kernel/solver swap invalidates), 36 tests including a real sbi round trip, the API entry, `examples/sbi/cached_fit.py` through the seam. Its **sbi gate was green (2501/55)** before the kill; its **dev gate is unknown**; **no `_sbi.py` wiring yet** (by design — it waits for W3.3). **Update after the resume**: W3.5 has **reported** — dev 1813/172 and sbi 2501/55 green, lint/format/pyrefly clean in both, docs unchanged; placement `ampere.results.artefacts` (it is provenance machinery, built on `provenance.py`'s hashes, beside `training.py`); a `model_hash` ingredient added beyond the item text so a kernel/solver/family swap with identical parameters is a miss (the item's `spec_hashes["spec"]` alone hashes only the parameter declaration); its `_sbi.py` wiring diff is in its report (transcript) and is ~16 new lines around the round loop — **the orchestrator applies it, adapted to W3.3's `layout=` and layout hash, when merging W3.5 after W3.3**. One real out-of-scope finding to schedule: `append_training_set` checks only `ampere_spec_hash` before appending, so a model change with identical parameter names could append onto a training set simulated under the old model — plan §7's trap in W2.8 code; a small follow-up item (compare `model_fingerprint` too).
**How to resume** (first session after 08:30): (1) `git worktree list`, `git branch --list 'w3.*'` — confirm the two branches above; (2) if the orchestrating session survives, `SendMessage` each agent to continue from its branch (they keep context); otherwise re-dispatch each with a prompt that says the branch exists and lists only what remains; (3) merge order is W3.3 first (review against its Accept line in WORK_ITEMS.md; all four gates on merged master since it touches `ampere/core`), then W3.5 (its wiring commit applied by the agent or by hand from its reported diff; gates dev + sbi); (4) the status table and this section at each merge.

**Dispatched 2026-09-09 from `3c034b9`, in parallel on disjoint files (record of the dispatch)**: **W3.3** (Opus, branch `w3.3-encoding`; owns `ampere/core/encoding.py`, `ampere/inference/_sbi.py` + `__init__.py`, `tests/core/test_encoding.py`, `tests/inference/test_sbi.py`, `encoding.md`) and **W3.5** (Sonnet, branch `w3.5-artefact-cache`; owns the new cache module, its tests, `examples/sbi/cached_fit.py`; wires `cache=` into `_sbi.py` only after W3.3 merges, else reports the diff). Both gate through `flock /tmp/ampere-gate.lock`. If cut off: `git branch --list 'w3.*'`, `git worktree list`; review branches with commits against their Accept lines; merge W3.3 first, then W3.5 (which may need its wiring commit applied by hand from its report). **Next for the orchestrator**: merge W3.3 then W3.5 as they report (gates one at a time: dev, sbi; jax/torch only if core files changed — W3.3 changes `ampere/core`, so all four) (disjoint: `ampere/core/encoding.py` + `_sbi.py`'s layout branch vs the cache module — W3.3 owns `_sbi.py`, W3.5 appends after W3.3 merges or works on a separate module and wires in at merge), W3.8 after slice 1; then W3.6; W3.4's maturity note any time after W3.2; W3.7 last; W3.10 whenever convenient. Policy unchanged (`docs/orchestration.md`): Opus for judgement-within-spec, Sonnet for well-specified items, Fable for contracts and reviews. One owner per shared file when dispatching in parallel; one five-suite gate at a time.

**Operations notes for the orchestrator**: this machine has 13 GB — never run two five-suite gates concurrently (the harness kills background tasks on low memory); run torch's `test-all` detached (`nohup … &`) and watch the log; after a lockfile change, `pixi install -e <env>` in the main checkout; the merged-master gate is the one that counts, re-run after every merge; **do not edit, format or check out files in the main checkout while a gate runs there** (a half-written module crashed a subprocess test on 2026-09-09); parallel agents appending decision-log rows conflict at the table's last line — keep both rows, keep the table contiguous; the two-track pattern works when each shared file has one owner and the sibling's additions are appended; the `sbi` gate (`pixi run -e sbi test-all`) joins the rotation from W3.2.

---

## Handoff as of the Phase 3 close (2026-09-10, superseded the same day by the clean rewrite; incremental record of the W3.12–W3.13 session)

## ⚡ Pick up here — PHASE 3 IN PROGRESS (2026-09-10): every item merged, W3.13 last; its merged-master gates running; then the close-out

**Where the project is.** Phases 0–2 complete (spec frozen at `spec-v1.0`; three backends; M2 reached). **Phase 3, the SBI layer, is more than half landed.** Merged, in order: W3.0 (housekeeping), W3.9 (migration page), W3.1 slice 1 (`simulate_many`, the executor protocol, external simulators), W3.2 (`SBIEngine` over sbi 0.27), W3.1 slice 2 (native batched prediction and sampling, sharding hooks, `forkserver`), W3.3 (the coordinate–value–mask encoding, `EncodingLayout`, the masked set/transformer embeddings), W3.5 (trained-artefact caching, wired into the engine), W3.8 (non-native parts refused by default, opt-in without gradients), **W3.6** (family D: SBC/TARP in `ampere.results.calibration`, `SBIEngine.calibrate`, the WStat coverage study — merged at `d25d8a5`; merged-master gates dev 1923/214, sbi 2652/70 green). Every merged item's row in WORK_ITEMS.md's status table is the authoritative record of what landed, what was accepted at review, and what is carried. **Master is `9e4bc7d`+**, clean apart from the untracked paper script under `examples/examples_paper/` (leave it). Last fully recorded gates: after the W3.4 merge (`ea7f1d7`, inference only) **dev 1943/248, sbi 2706/70**; after the W3.8 merge (`9e4bc7d`, the last to touch core) **jax 2375/136, torch 2546/134**, lint/format clean, pyrefly 0 errors in every environment. **Nothing is pushed to origin** (origin/master is still at the W1.3 handoff, `b8e585b`).

**Merged on 2026-09-09/10, in order**: **W3.12** at `449eda2` (Sonnet; one fix at review — `sample_with` joins the artefact key; merged-master gates dev 1982/258, jax 2477/191, torch 2665/187 green, sbi 2759/80 with one failure — a literal `== 5` schema pin in `tests/inference/test_sbi.py` that W3.12's worktree, lacking an sbi environment, could never run: the orchestrator's miss, fixed on W3.15's branch); **W3.7** at `a744009` (the blocking `new-namespace suites (sbi)` leg; no gates); **W3.10** at `1f4e8d3` (paging; merged-master dev 1994/258); **W3.15** at `1fb3097` (Sonnet; torch *and* numpy's legacy generator seeded per step from the problem and restored; branch gates dev 1982/264, sbi 2766/80; carries the schema-pin fix). **W3.15's merged-master gates: sbi 2778/80, dev 1994/264, both green** and recorded. **Last fully recorded gates on master (`1fb3097`)**: sbi 2778/80, dev 1994/264; jax 2477/191 and torch 2665/187 from `449eda2` (nothing touching core or a backend has merged since). Lint/format/pyrefly clean on every merged branch. **W3.13 merged** at `0ff223e` (Sonnet; one fix at review, `e0d0cf7`; branch gates dev 1994/268, sbi 2782/80; docs 13 warnings unchanged). **Its merged-master gates are running** — dev then sbi (its `tests/examples` addition), `<scratchpad>/gates-w313m.summary` — record them in W3.13's status row (the `MERGED_GATES_W313` placeholder). **Nothing is in flight.**

**Next, in order.** When W3.13 merges, write the **Phase 3 close-out**: the plan's §5 Phase 3 section marked complete with the date, the deferred-list state, the questions list below handed to Peter as one block, and this section rewritten for a clean start on the next phase (the plan's §5 says 3 and 4 were parallelisable; Phase 4 — one new modality end-to-end — is the next to draft items for, with `docs/design/modalities/` as its input). Then the Phase 4 item drafts, for Peter's approval before any dispatch. The owed small items below are Sonnet-sized fillers if a session has capacity while waiting on gates.

**Decisions recorded this phase** (all in the plan's §2 or the status rows): `forkserver` default; native sampling the default; further core `sample` families on a use case with the pathway kept open (structurally enforced by W3.8); the mesh sharder's padding to be costed with the GPU item; the observation-context reservation (W3.3's encoding carries σ by default; `context=` accepted as `None` and recorded); kernel algebra into Phase 4; the embedding *study* (width and readout vs data/model/structure, guidance for users) in the deferred list; W2.9's coverage study inside W3.6.
