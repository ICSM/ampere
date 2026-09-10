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

## ⚡ Pick up here — PHASE 3 IN PROGRESS (2026-09-10): every item merged, W3.13 last; its merged-master gates running; then the close-out

**Where the project is.** Phases 0–2 complete (spec frozen at `spec-v1.0`; three backends; M2 reached). **Phase 3, the SBI layer, is more than half landed.** Merged, in order: W3.0 (housekeeping), W3.9 (migration page), W3.1 slice 1 (`simulate_many`, the executor protocol, external simulators), W3.2 (`SBIEngine` over sbi 0.27), W3.1 slice 2 (native batched prediction and sampling, sharding hooks, `forkserver`), W3.3 (the coordinate–value–mask encoding, `EncodingLayout`, the masked set/transformer embeddings), W3.5 (trained-artefact caching, wired into the engine), W3.8 (non-native parts refused by default, opt-in without gradients), **W3.6** (family D: SBC/TARP in `ampere.results.calibration`, `SBIEngine.calibrate`, the WStat coverage study — merged at `d25d8a5`; merged-master gates dev 1923/214, sbi 2652/70 green). Every merged item's row in WORK_ITEMS.md's status table is the authoritative record of what landed, what was accepted at review, and what is carried. **Master is `9e4bc7d`+**, clean apart from the untracked paper script under `examples/examples_paper/` (leave it). Last fully recorded gates: after the W3.4 merge (`ea7f1d7`, inference only) **dev 1943/248, sbi 2706/70**; after the W3.8 merge (`9e4bc7d`, the last to touch core) **jax 2375/136, torch 2546/134**, lint/format clean, pyrefly 0 errors in every environment. **Nothing is pushed to origin** (origin/master is still at the W1.3 handoff, `b8e585b`).

**Environments.** `dev`, `torch`, `jax`, and **`sbi` (= dev + torch feature + sbi extra since W3.2; run `pixi install -e sbi` after pulling)**. The gate is `pixi run -e <env> test-all` in each; an item's gate set is dev + sbi, plus torch and jax when it touches `ampere/core` or a backend. **Gates are serialised through `flock /tmp/ampere-gate.lock`** — every agent and the orchestrator run `flock /tmp/ampere-gate.lock pixi run -e <env> test-all`, so two five-suite runs never overlap on this 13 GB machine (the harness kills background tasks on low memory). Torch's takes 12–17 min, sbi's 14–20, jax's 8–13, dev's 5–8.

**Merged on 2026-09-09/10, in order**: **W3.12** at `449eda2` (Sonnet; one fix at review — `sample_with` joins the artefact key; merged-master gates dev 1982/258, jax 2477/191, torch 2665/187 green, sbi 2759/80 with one failure — a literal `== 5` schema pin in `tests/inference/test_sbi.py` that W3.12's worktree, lacking an sbi environment, could never run: the orchestrator's miss, fixed on W3.15's branch); **W3.7** at `a744009` (the blocking `new-namespace suites (sbi)` leg; no gates); **W3.10** at `1f4e8d3` (paging; merged-master dev 1994/258); **W3.15** at `1fb3097` (Sonnet; torch *and* numpy's legacy generator seeded per step from the problem and restored; branch gates dev 1982/264, sbi 2766/80; carries the schema-pin fix). **W3.15's merged-master gates: sbi 2778/80, dev 1994/264, both green** and recorded. **Last fully recorded gates on master (`1fb3097`)**: sbi 2778/80, dev 1994/264; jax 2477/191 and torch 2665/187 from `449eda2` (nothing touching core or a backend has merged since). Lint/format/pyrefly clean on every merged branch. **W3.13 merged** at `0ff223e` (Sonnet; one fix at review, `e0d0cf7`; branch gates dev 1994/268, sbi 2782/80; docs 13 warnings unchanged). **Its merged-master gates are running** — dev then sbi (its `tests/examples` addition), `<scratchpad>/gates-w313m.summary` — record them in W3.13's status row (the `MERGED_GATES_W313` placeholder). **Nothing is in flight.**

**Next, in order.** When W3.13 merges, write the **Phase 3 close-out**: the plan's §5 Phase 3 section marked complete with the date, the deferred-list state, the questions list below handed to Peter as one block, and this section rewritten for a clean start on the next phase (the plan's §5 says 3 and 4 were parallelisable; Phase 4 — one new modality end-to-end — is the next to draft items for, with `docs/design/modalities/` as its input). Then the Phase 4 item drafts, for Peter's approval before any dispatch. The owed small items below are Sonnet-sized fillers if a session has capacity while waiting on gates.

**Standing dispatch paragraph** (put in every prompt): base commit and the ancestor/descendant recovery recipe from `docs/development.md` "Dispatching work"; the `flock` gate rule and "never stop your turn to wait — poll with a foreground `until`/`sleep` loop"; explicit file ownership when two agents run in parallel; do not edit the WORK_ITEMS.md status table or DEVELOPMENT_PLAN.md (decision-log text goes in the report); the commit trailer; the report format from the `work-item` skill. Models: Opus for judgement-within-spec, Sonnet for well-specified items (`docs/orchestration.md`).

**Decisions recorded this phase** (all in the plan's §2 or the status rows): `forkserver` default; native sampling the default; further core `sample` families on a use case with the pathway kept open (structurally enforced by W3.8); the mesh sharder's padding to be costed with the GPU item; the observation-context reservation (W3.3's encoding carries σ by default; `context=` accepted as `None` and recorded); kernel algebra into Phase 4; the embedding *study* (width and readout vs data/model/structure, guidance for users) in the deferred list; W2.9's coverage study inside W3.6.

**Everything that needs Peter's attention, consolidated at the Phase 3 close (2026-09-10).** Rulings go in the status rows / decision log; nothing here blocks the close-out.

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
