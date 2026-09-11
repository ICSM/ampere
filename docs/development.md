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
    `pixi run test-all` (every new-namespace suite — core, results,
    conformance, backends, inference and, since W2.11, examples — in ONE
    pytest process; **this is the main gate**, and `test-phase1` remains
    the older core+results+conformance subset; the individual `test-core`,
    `test-results`, `test-backends`, `test-inference`, `test-examples` and
    `conformance` tasks exist too),
    `pixi run test-characterisation` (legacy still works),
    `pixi run lint`, `pixi run format-check`, `pixi run typecheck`,
    `pixi run docs` (repaired at W2.11 — it builds), `pixi run bench` (the
    tracked benchmarks; writes `benchmark.json`, which CI uploads as a run
    artefact) and `pixi run scaling` (the on-demand 10³–10⁵ GP scaling
    demonstration, deliberately in no gate). Since W0.10, plain
    `pixi run <task>` equals `-e dev`.
  - Other environments (`pixi run -e <env> <task>`): `test-py311` /
    `test-py312` / `test-py313` (the CI matrix, one Python each); `sbi`
    (since W3.2 the `dev` set plus the `torch` feature and the `sbi`
    extra — sbi 0.27, torch, pyro — and a real gate of its own:
    `pixi run -e sbi test-all`, plus `pixi run -e sbi typecheck`); `torch` /
    `jax` (the backend environments — real gates since W2.4/W2.5:
    `pixi run -e torch test-all`, `pixi run -e jax test-all`, and each
    backend package is typechecked with its real types only in its own
    environment). Since W2.11 both `torch` and `sbi` resolve torch from
    PyTorch's **CPU** index, so neither pulls ~2 GB of CUDA runtime onto a
    CPU-only runner; this changes nothing about what `pip install
    "ampere[torch]"` gives a user. **Run five-suite gates one at a time** —
    two concurrently exhaust a 13 GB machine. Torch's five-suite gate takes
    12–17 min here, sbi's 14–20 min, jax's 8–13 min, dev's 5–8 min.
  - **CI's `suites` matrix** (`.github/workflows/ci.yml`): one blocking
    job per new-namespace environment, each producing a check named
    "new-namespace suites (`<environment>`)" — `dev` (its own `suites` job:
    `test-all` + `bench`) and, in the `backend-suites` job's matrix,
    `torch`, `jax` and, since W3.7, `sbi` (each: `typecheck` then
    `test-all`; `torch`/`jax` also run `bench`, deliberately not repeated
    for `sbi` since it would only re-exercise torch's own benchmarks under
    another name — see the job's comment). `sbi`'s install is cached the
    same way as `torch`'s (`setup-pixi`'s `cache: true`). The weekly/
    on-demand `sbi-characterisation` job is unrelated: it exercises the
    *legacy* `ampere.infer.sbi` flow, not the new-namespace `sbi` gate, and
    stays non-blocking. Because the `sbi` leg's local reference time
    (14–20 min) sits above the `torch` leg's (12–17 min), a CI runner with
    torch's-leg headroom does not automatically have sbi's-leg headroom —
    see W3.7's report for the cheap levers considered (splitting
    `typecheck` out of the job, tightening the smoke budgets further)
    without widening that item's scope to touch test budgets.
- Plain-pip alternative: `pip install -e ".[dev]"`, Python ≥ 3.11.
- Conda env `ampere` (Python 3.13) has an editable install pointing at the
  main checkout — worktree-based work that needs importing its own changes
  should use pixi or `pip install -e .` into a fresh env instead.
- `pytest tests/characterisation` (or `pixi run test-characterisation`) is
  the "legacy still works" gate; run it before merging anything that
  touches shared files.

## ⚡ Pick up here — PHASE 4 IN PROGRESS (first wave dispatched 2026-09-11 from `63a319a`): W4.0 ∥ W4.1 ∥ W4.5 ∥ W4.6 in agent worktrees; nothing merged yet

**2026-09-11 (Fable, orchestrating session).** Peter ruled all four Phase 4 decisions (record: `docs/design/phase4_placement_memo.md`; applied at `63a319a`): **D1** the kind in core beside `VisibilitySet`, the steps in `backends/{reference,torch,jax}/interferometry.py`, no grouping namespace now (revisit after realistic usage); **D2** as revised after his chromatic-misspecification question — `VisibilitySet` gains a `spectral_axis`, `ClosurePhases` is `(u1, v1, u2, v2, spectral_axis)` with a canonical baseline ordering, kernels gain an `axes` selector so `Product` composes a (u, v) block with a spectral block; a GP on closure phases is latent and is **Phase 5's (W5.1)**; **D3** W4.9 runs; **D4** agreed, with the 2026-09-05 review-policy ruling extended to Phase 4 (W4.2/W4.5 merge on a second Fable pass while the Codex quota is blocked, terra owed). New items: **W4.11** (the photometry + spectrum composition as a documented example, from memo §7.1 — after W4.0), **W5.1**; the plan's Phase 5 notes matrix-free exact GPs (gpytorch/gpjax-style linear-operator solves) as a third solver strategy for the 3- and 5-axis kernels, not immediate.

**In flight (dispatched 2026-09-11, all from base `63a319a`, worktrees under `.claude/worktrees/`)**: `w4.0-phase4-housekeeping` (Sonnet), `w4.1-interferometry-reference` (Opus), `w4.5-kernel-algebra` (Opus), `w4.6-astropy-compat` (Opus). File ownership is in each prompt and in the items file's "Ordering" paragraph; W4.1 and W4.5 both edit `core/likelihood.py` (families vs kernels) — expect a hand merge there and on `core/__init__.py` exports. Gates serialise through `flock /tmp/ampere-gate.lock`; each agent runs dev, sbi, torch, jax. **A clean session does, in order**: (1) `git worktree list` and `git branch --list 'w4.*'` to see what survived; an agent cut off mid-flight is resumed from its branch (its committed work stands; re-dispatch with the same prompt and "continue from branch X"); (2) review each finished branch against its Accept line with the `docs/development.md` checklist, run the merged-master gates for every environment the item names, merge, update the status table, remove the worktree (`git worktree remove --force --force`); (3) after W4.0 merges dispatch **W4.11** (Sonnet); after W4.5 and W4.1 merge dispatch **W4.2** (Opus); after W4.6, **W4.7** (Sonnet); after W4.1 + W4.2, **W4.3** (Opus); W4.10 (Sonnet) any time; (4) refresh this section at every merge and dispatch. The scratchpad script behind memo §7.1 is not in the repo; W4.11 recreates it from the memo.

**Where the project is.** Phases 0–3 complete: the spec frozen at `spec-v1.0` (2026-09-03), three backends in lockstep, M2 reached, and **the SBI layer closed 2026-09-10** — every Phase 3 item W3.0–W3.15 merged, W3.13 (the documentation pass) last at `0ff223e`. Each item's row in WORK_ITEMS.md's status table is the authoritative record of what landed, what was accepted at review and what is carried; the plan's §5 Phase 3 section carries the landed summary and §2 the decision-log rows. **Gates on master, all green, all four after the last code merge (W3.16, `fea7ef9`)**: dev 1996/268, sbi 2784/80, jax 2491/201, torch 2667/197. Lint/format/pyrefly clean in every environment. Lint/format/pyrefly clean in every environment. **Nothing is in flight, no agent worktrees exist, and nothing is pushed to origin** (origin/master is still at the W1.3 handoff, `b8e585b`). Master is clean apart from Peter's untracked paper script under `examples/examples_paper/` (leave it).

**Environments.** `dev`, `torch`, `jax`, and **`sbi` (= dev + torch feature + sbi extra since W3.2; run `pixi install -e sbi` after pulling)**. The gate is `pixi run -e <env> test-all` in each; an item's gate set is dev + sbi, plus torch and jax when it touches `ampere/core` or a backend. **Gates are serialised through `flock /tmp/ampere-gate.lock`** — every agent and the orchestrator run `flock /tmp/ampere-gate.lock pixi run -e <env> test-all`, so two five-suite runs never overlap on this 13 GB machine (the harness kills background tasks on low memory). Torch's takes 12–17 min, sbi's 14–20, jax's 8–13, dev's 5–8.

**What a clean session does next, in order.** (1) Read the "Everything that needs Peter's attention" block below first; any ruling Peter has given since becomes either a direct edit (the one-liners in part A: the `sample_with` default, the zero-rate guard, `cauchy`) or a small item — do those before anything new. (2) **The Phase 4 items are drafted** (WORK_ITEMS.md "Phase 4", `5488e6d`: W4.0–W4.9 and four decisions D1–D4 — the `ampere/modalities/` namespace, the `ClosurePhases` signature, whether W4.9 runs, the model split) **and await Peter's approval before any dispatch**; once approved, dispatch W4.0 ∥ W4.1 ∥ W4.5 ∥ W4.6 from the standing paragraph. The inputs they were drafted from: the inputs are the plan's §5 Phase 4 bullets (interferometric visibilities end to end — decided; the astropy interop adapter, §4.7; kernel algebra and the public quasiseparable-term registry, ruled into Phase 4 on 2026-09-09) and `docs/design/modalities/interferometry.md`, whose §9 (interface gaps), §10 (requirements on the joint fit) and §11 (open questions) are the item list in all but name; `docs/design/modalities/astrometric_timeseries.md` and `ifu_cube.md` are the follow-on templates. Draft in W2/W3 style — one branch per item, an Accept line each, a conformance row per new container or transformation, the gate set named — and size them S/M for Sonnet/Opus per `docs/orchestration.md`. (3) If Peter approves the memo's §5 adaptations or §9 horizons (part B), draft those as Phase 5 items in the same pass. (4) The owed small items below are Sonnet-sized fillers while gates run.

**Phase 3 in one paragraph, for whoever reads only this.** `ampere.inference.SBIEngine` over sbi 0.27 (NPE/NLE/NRE, and TMNRE through sbi's own ratio estimators — no swyft, ruled) fits any `FittingProblem` through `simulate_many`'s executor protocol (`forkserver` default, native batched prediction and sampling on torch and jax); the coordinate–value–mask encoding and the masked set/transformer embeddings; trained-artefact caching keyed on spec, model and data hashes plus the run's settings; SBC/TARP calibration usable by every engine; non-native parts refused by default and opt-in without gradients; core `sample()` for Poisson, Student-t and complex Gaussian with native twins; `PROVENANCE_SCHEMA_VERSION` 6 with `ampere_model_hash` and the append check; runs reproducible bitwise from the problem's seed; the blocking sbi CI leg; plot paging; the SBI tutorial. Carried into the deferred list: the embedding study with the ε study folded in; jax-native SBI on maturity.

**Standing dispatch paragraph** (put in every prompt): base commit and the ancestor/descendant recovery recipe from `docs/development.md` "Dispatching work"; the `flock` gate rule and "never stop your turn to wait — poll with a foreground `until`/`sleep` loop"; explicit file ownership when two agents run in parallel; do not edit the WORK_ITEMS.md status table or DEVELOPMENT_PLAN.md (decision-log text goes in the report); the commit trailer; the report format from the `work-item` skill. Models: Opus for judgement-within-spec, Sonnet for well-specified items (`docs/orchestration.md`).

**Decisions recorded this phase** (all in the plan's §2 or the status rows): `forkserver` default; native sampling the default; further core `sample` families on a use case with the pathway kept open (structurally enforced by W3.8) — met at W3.14; the mesh sharder's padding to be costed with the GPU item; the observation-context reservation (W3.3's encoding carries σ by default; `context=` accepted as `None` and recorded); kernel algebra into Phase 4; the embedding *study* and the ε study in the deferred list; W2.9's coverage study inside W3.6; the `marginals` group stored by default (W3.4); the Poisson/Student-t/complex-Gaussian sampling forms (W3.14); schema 6 and the model hash (W3.12); plot paging (W3.10); the seeding rule (W3.15); the columnar store for population work waits but nothing may exclude it (2026-09-10, horizon (b)).

**Everything that needs Peter's attention, consolidated at the Phase 3 close (2026-09-10).** Rulings go in the status rows / decision log; nothing here blocks the close-out.

**Rulings received 2026-09-10 (Peter): parts A and B are all ruled and applied** — A1 and A6 landed as **W3.16** (merged `fea7ef9`; merged-master gates all green, recorded); A2–A4, A7, A8 confirmed in the decision-log rows; A5's principle in `likelihoods.md` §3; A9 as W4.10; A10 as W4.0 (1) and (6); B11–B14 and B17 in the memo's §10; B15's horizons (e)–(g) in the plan; B16 as W5.0. **Only part C remains open**, for Peter when he gets the chance. **Phase 3 is finished** (2026-09-10). The horizon notes are folded into the plan's §5 Phase 5 bullets and horizons (h)–(i) (`152d187`), so the whole horizon is visible in one place.

*A. Rulings on merged Phase 3 work (each is a one-line change or a "leave it")*
1. **W3.4 — `sample_with` default for TMNRE.** Current default `"rejection"` (i.i.d. draws). Measured 353.6 s against 43.8 s for `"mcmc"` on the example, and rejection's cost *rises* as truncation succeeds; `"mcmc"` is sbi's own default for ratio posteriors and what `method="nre"` uses. Recommendation: switch the default to `"mcmc"`, keep `"rejection"` selectable.
2. **W3.4 — the accuracy row.** The 0.5σ location tolerance is unreachable for a ratio estimator's joint on the toy problem; the suite asserts coverage plus the width band and holds the marginals to 1.0σ. Confirm, or state a tolerance.
3. **W3.4 — `ε = 1e-4`** gives a ±4σ box that stops shrinking after round 2. Recommendation: leave; the ε study is now in the deferred list (W3.13 entered it beside the embedding study).
4. **W3.4 — per-round pair estimators.** The 2-D marginal estimators train only on the final round; per-round pairs would cost budget for a box the 1-D marginals already set. Confirm the final-round-only choice.
5. **W3.14 — `cauchy`.** Refused only because the item listed three families; it is one branch away. Rule: a list (leave) or a principle (add it, and write the principle into `likelihoods.md` §3).
6. **W3.14 — the zero-rate asymmetry.** `PoissonFamily.sample` guards `rate < 0`, `log_prob` guards `rate <= 0`. Recommendation: align `sample` to `<= 0`.
7. **W3.12 — schema 6** (`ampere_model_hash` on every run and training set; a pre-schema-6 training set refused on append by name rather than guessed) and the `model_hash` name in place of the item's `model_identity_hash(problem)` (which already existed with other semantics). Confirm.
8. **W3.15** was drafted and dispatched on Fable's judgement: torch *and* numpy's legacy global generator seeded per step and restored. Confirm retroactively, or strike.
9. **W3.7 — the sbi CI leg** is likely the slowest (local 14–20 min against torch's 12–17). Accept, or apply a lever (typecheck split into its own job; tighter smoke budgets). Also: add `actionlint` to the dev toolchain?
10. **W3.13** — (a) `docs/source/overview.rst` is framed "as of M2" and still says five engines and "SBI is Phase 3": approve a small refresh item (W3.16). (b) No `examples/sbi` script does a bare in-process NPE fit; the tutorial covers the engine through the black-box case. Want a plain native example?

*B. The inference-extensions memo (`docs/design/inference_extensions_memo.md`)*
11. **§7.1 — the stored proposal density** as a rule for every approximate engine (VI today does not store it; SBI does). The hook with the widest consequences: importance correction, population reweighting at scale, retroactive density emulation. Yes/no.
12. **§7.2 — optimisation results**: a `DataTree` with an `optimum` group, or a typed object with `to_datatree()`.
13. **§7.3 — BO as an engine** or only as acquisition for emulation/multi-fidelity (recommendation: the latter; the *surrogate-posterior* family — VBMC, GPry — is filed separately at tier 2 as the route for simulators too expensive for SBI).
14. **§7.4 — dependency policy**: may `nautilus`/`ultranest` join the base install, or every new sampler behind an extra as `zeus` is?
15. **§9 — three new design horizons** proposed for the plan's list: (e) model comparison and averaging, (f) approximate-inference correction, (g) engines as a registry with uniform cost accounting. Approve and Fable adds them as one paragraph each; the per-horizon reservations in §9 (a fidelity column in training sets; the interim prior reconstructible from the archive; an emulator's identity including its training set; nested encodings) are sentences for the contracts they name, added by the first item that touches each.
16. **§5's six adaptations** — should they become drafted Phase 5 items now (the evidence attr, the weighted/approximate-draw rule, `ampere_approximation`, the `Optimum` result, batched `log_prob`, dependency tiers), or wait for the first engine that needs them? Recommendation: draft the first three as one small item, since they are results-contract text and the samplers of tier 1 want them.
17. **§8.6 / horizon (b) — the columnar store**: ruled "wait, but nothing may exclude it"; recorded as a constraint on horizon (b). Nothing further unless the ruling changes.

*C. Repository actions only Peter takes*
18. Push master (origin is still at the W1.3 handoff `b8e585b`) and watch the first live CI run: the branch-protection rule naming the old `phase1-suites` job needs `suites`, and **a required-check entry for `new-namespace suites (sbi)`** beside `(torch)`/`(jax)` (W3.7).
19. Delete the harness branches: `git branch --list 'worktree-agent-*' | xargs -r git branch -d`.
20. `pixi run gpu` on a machine with an accelerator, and cost the mesh sharder's padding there (W3.1 slice 2's ruling).
21. The retroactive `gpt-5.6-terra` reviews when the Codex quota returns (~2026-09-30): the Phase 2 queue (W2.3 first), then W3.1 slices 1–2, W3.2, W3.3, W3.8, W3.6's calibration statistics, W3.12's hash composition, W3.14's sampling forms.
22. `examples/examples_paper/flexible_likelihood_comparison.py` is untracked on master — his paper work; commit or ignore as he sees fit.

*D. Next phase*
23. **Phase 4 (one new modality end-to-end — interferometric visibilities, decided)** is the next to draft items for, with `docs/design/modalities/` as the input; Phase 5's inference items would follow from B. Fable drafts for approval before any dispatch.

**Owed / carried small items** (not scheduled; pick up when convenient, Sonnet-sized): a native sampler failure aborts a whole chunk where the numpy loop flags one draw (`_SimulateNatively.run`, W3.14's finding — wrap the sampler in the per-draw `try`); `Instrument._declarations()` fingerprints steps by `id()`, so a *bare* pickled `Instrument` still breaks (a value-based fingerprint in `transform.py`; `Dataset.__setstate__` works around it); torch `noise.py`'s `_check_one_device` partly redundant with the kernel-as-capability-part check; the "no uncertainty at all" branch of both backends' `_check_gp_uncertainty` is unreachable; upstream sbi 0.27 issues worth filing — `PermutationInvariantEmbedding` reads the first batch element's valid-row count for the whole batch, and `TransformerEmbedding` drops `attention_mask` unless causal; `z_score_x="none"` does not silence sbi's constant-column warning; the Phase 6 multi-observation composition tutorial (plan §5); the `tests/examples` smoke test of photometry + spectra composition.

**Operations notes for the orchestrator**: one five-suite gate at a time, always through the `flock`; run long gates detached (`nohup … &`) and wait with a background `until` loop rather than polling in the foreground; **never edit, format or check out source files in the main checkout while a gate runs there** (a half-written module crashed a subprocess test on 2026-09-09; markdown is safe); an agent that `cd`s to the main checkout can edit its files even though the isolation guard blocks git there — tell agents to stay in their worktree; after a lockfile change, `pixi install -e <env>` in the main checkout; the merged-master gate is the one that counts; parallel agents appending decision-log rows conflict at the table's last line — keep both rows, keep the table contiguous; agent worktrees under `.claude/worktrees/` may hold their own `.pixi` (2 GB) — remove a merged agent's worktree with `git worktree remove --force --force`; finished agents linger in the task list until dismissed with `TaskStop`; an agent that "stops to wait" for a gate is resumed with `SendMessage` and told to poll; a session rate limit kills agents mid-flight but their committed work survives — resume them from their branches.


---

Earlier session records (Phase 0–1 review outcomes, the 2026-09-01 to 2026-09-09 handoffs, including the incremental mid-Phase-3 record) are archived verbatim in `docs/handoff-archive.md`.
