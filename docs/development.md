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

## ⚡ Pick up here — PHASE 3 IN PROGRESS (2026-09-09): W3.0–W3.6, W3.8, W3.11 merged; W3.4, W3.14 merged; W3.12 in flight; W3.15, W3.7, W3.10 queued

**Where the project is.** Phases 0–2 complete (spec frozen at `spec-v1.0`; three backends; M2 reached). **Phase 3, the SBI layer, is more than half landed.** Merged, in order: W3.0 (housekeeping), W3.9 (migration page), W3.1 slice 1 (`simulate_many`, the executor protocol, external simulators), W3.2 (`SBIEngine` over sbi 0.27), W3.1 slice 2 (native batched prediction and sampling, sharding hooks, `forkserver`), W3.3 (the coordinate–value–mask encoding, `EncodingLayout`, the masked set/transformer embeddings), W3.5 (trained-artefact caching, wired into the engine), W3.8 (non-native parts refused by default, opt-in without gradients), **W3.6** (family D: SBC/TARP in `ampere.results.calibration`, `SBIEngine.calibrate`, the WStat coverage study — merged at `d25d8a5`; merged-master gates dev 1923/214, sbi 2652/70 green). Every merged item's row in WORK_ITEMS.md's status table is the authoritative record of what landed, what was accepted at review, and what is carried. **Master is `9e4bc7d`+**, clean apart from the untracked paper script under `examples/examples_paper/` (leave it). Last fully recorded gates: after the W3.4 merge (`ea7f1d7`, inference only) **dev 1943/248, sbi 2706/70**; after the W3.8 merge (`9e4bc7d`, the last to touch core) **jax 2375/136, torch 2546/134**, lint/format clean, pyrefly 0 errors in every environment. **Nothing is pushed to origin** (origin/master is still at the W1.3 handoff, `b8e585b`).

**Environments.** `dev`, `torch`, `jax`, and **`sbi` (= dev + torch feature + sbi extra since W3.2; run `pixi install -e sbi` after pulling)**. The gate is `pixi run -e <env> test-all` in each; an item's gate set is dev + sbi, plus torch and jax when it touches `ampere/core` or a backend. **Gates are serialised through `flock /tmp/ampere-gate.lock`** — every agent and the orchestrator run `flock /tmp/ampere-gate.lock pixi run -e <env> test-all`, so two five-suite runs never overlap on this 13 GB machine (the harness kills background tasks on low memory). Torch's takes 12–17 min, sbi's 14–20, jax's 8–13, dev's 5–8.

**W3.4 merged** at `ea7f1d7` (merged-master gates dev 1943/248, sbi 2706/70 green). **W3.14 merged** at `4643f53` (four merged-master gates queued at the merge — `<scratchpad>/gates-w314.summary`; if unrecorded in its status row, re-run one at a time). **In flight.** **W3.12** (Sonnet, branch `w3.12-model-hash`, from `ea7f1d7`; `model_identity_hash`, schema 6, the append check, W3.4's `artefact_key` diff; owns `ampere/results/*`, the cache call in `_sbi.py`, `tests/results/*`; all four gates). Disjoint files. Merge whichever reports first; both need all four gates on merged master. *(Record of the W3.4 dispatch follows.)* **W3.4** (Opus, branch `w3.4-tmnre`, dispatched from the W3.11 merge): TMNRE through sbi's own `NRE` + `RestrictedPrior` rounds — `SBIEngine(method="tmnre", rounds=, marginals=, truncation_epsilon=)`, the truncation box per round in the attrs, the marginal posteriors emitted per the item's design paragraph; owns `ampere/inference/_sbi.py` (+ a new `_tmnre.py` if it prefers), `tests/inference/test_sbi.py`, `inference.md`/`results.md` amendments, `docs/source/ampere.inference.rst`. If cut off: `git branch --list 'w3.4*'`; review against its Accept line; merge with dev + sbi gates. (W3.11 merged before it at `531a518`; merged-master gates dev 1923/217, sbi 2655/70 green.)

**Next, in order** (all ruled unless marked; item texts in WORK_ITEMS.md "Phase 3"; dispatch prompts are regenerated from the item text plus the standing paragraph below): after W3.4 merges — ~~**W3.4**~~ (in flight); then (Opus; TMNRE through sbi's own `NRE` + `RestrictedPrior` rounds, *no swyft* — ruled 2026-09-09; the design paragraph on the item is the scope), **W3.11** (Sonnet; set-embedding default width `max(2·free_size, 32)` and a masked-mean transformer readout), **W3.12** (Sonnet; `model_identity_hash` promoted into `provenance.py`, recorded on runs and training sets, `PROVENANCE_SCHEMA_VERSION` → 6, `append_training_set` refuses a model change — goes **alone or last**, every run's attrs change); all three touch `_sbi.py` or `results/`, so W3.4 ∥ W3.11 is fine only with one owner of `_sbi.py` (give it to W3.4; W3.11 appends after) — the order is now W3.12 → W3.15 → W3.7 → W3.10, each a short session; W3.12 also lands W3.4's `artefact_key` diff. **W3.14** (core `sample()` for Poisson, Student-t, complex Gaussian with native twins) is **approved (2026-09-09) and dispatched** beside W3.4 (Opus, branch `w3.14-generative-forms`, from `89ab1b7`; owns `ampere/core/likelihood.py`, the backends' families/noise/problem sampling code, `tests/core|backends|conformance`, `inference.md` §13 and `likelihoods.md`, `examples/wstat_comparison.py`'s `CountingPoisson` deletion) — disjoint from W3.4's `ampere/inference/*`. Merge order: whichever reports first; W3.14 needs all four gates. Then **W3.7** (Sonnet; the blocking `sbi` CI leg) and **W3.10** (Sonnet; automatic paging above the plot caps, not urgent). Phase 3 closes with a documentation pass like W2.15 (the SBI tutorial, `examples/sbi/` on the docs site, the deferred-list review) — draft that as W3.13 when W3.7 lands.

**Standing dispatch paragraph** (put in every prompt): base commit and the ancestor/descendant recovery recipe from `docs/development.md` "Dispatching work"; the `flock` gate rule and "never stop your turn to wait — poll with a foreground `until`/`sleep` loop"; explicit file ownership when two agents run in parallel; do not edit the WORK_ITEMS.md status table or DEVELOPMENT_PLAN.md (decision-log text goes in the report); the commit trailer; the report format from the `work-item` skill. Models: Opus for judgement-within-spec, Sonnet for well-specified items (`docs/orchestration.md`).

**Decisions recorded this phase** (all in the plan's §2 or the status rows): `forkserver` default; native sampling the default; further core `sample` families on a use case with the pathway kept open (structurally enforced by W3.8); the mesh sharder's padding to be costed with the GPU item; the observation-context reservation (W3.3's encoding carries σ by default; `context=` accepted as `None` and recorded); kernel algebra into Phase 4; the embedding *study* (width and readout vs data/model/structure, guidance for users) in the deferred list; W2.9's coverage study inside W3.6.

**Questions collected for Peter's next session** (Peter, 2026-09-09: "continue sequentially through as many of the remaining tasks as you can … collect open questions for me to review later"; this list is appended to at every merge and is the first thing that session reads):
1. **W3.4 — `sample_with` default for TMNRE.** Ruled design and current default: `"rejection"` (i.i.d. draws, what the single "chain" promises). Measured: 353.6 s against 43.8 s for `"mcmc"` on the example, and rejection's cost *rises* as truncation succeeds. `"mcmc"` is sbi's own default for ratio posteriors and what `method="nre"` already uses. Fable's recommendation: switch the default to `"mcmc"`, keep `"rejection"` selectable. One line plus two test assertions.
2. **W3.4 — the accuracy row.** The 0.5σ location tolerance is unreachable for a ratio estimator's joint on the toy problem; the suite asserts coverage (emcee mean inside the central 95 %) plus the width band, and holds the marginals to 1.0σ. Fable accepted it as the stronger pair; confirm or ask for a stated tolerance.
3. **W3.4 — `ε = 1e-4`** gives a generous box (about ±4σ) and boxes stop shrinking after round 2; a larger ε concentrates harder at the risk of clipping tails. Recommendation: leave, and add an ε study to the deferred embedding study.
4. **W3.14 — `cauchy`.** It keeps the refusal only because the item listed three families; a Cauchy is as unambiguous as a Student-t (`predicted + σ·standard_cauchy`) and is one branch away. Rule: a list (leave) or a principle (add it, and write the principle into `likelihoods.md` §3).
5. **W3.14 — the zero-rate asymmetry.** `PoissonFamily.sample` guards `rate < 0` (the item's wording) while `log_prob` guards `rate <= 0`, so a zero rate is drawable but not scoreable. Recommendation: align `sample` to `<= 0`.
6. **W3.15 (drafted by Fable, dispatched after W3.12 on Fable's judgement)**: seeding torch from the problem so SBI runs are reproducible — confirm the scope or strike.

**For Peter — open, none blocking**: W3.12 will bump the schema — confirm at its merge review; the retroactive `gpt-5.6-terra` passes when the Codex quota returns (~2026-09-30): the Phase 2 queue (W2.3 first) then W3.1 slices 1–2 (the batched-equals-loop, partition-independence and native-sampling claims), W3.2 (the unconstrained prior and draw scoring), W3.3 (the encoding), W3.8 (the structural backend resolution), and W3.6's calibration statistics when it lands; the repository actions only he takes — push master and watch the first live CI run (a branch-protection rule naming the old `phase1-suites` job needs `suites`; the `sbi` leg becomes blocking at W3.7), delete the harness `worktree-agent-*` branches (`git branch --list 'worktree-agent-*' | xargs -r git branch -d`), run `pixi run gpu` on a machine with an accelerator (and cost the mesh sharder's padding there).

**Owed / carried small items** (not scheduled; pick up when convenient, Sonnet-sized): a native sampler failure aborts a whole chunk where the numpy loop flags one draw (`_SimulateNatively.run`, W3.14's finding — wrap the sampler in the per-draw `try`); `Instrument._declarations()` fingerprints steps by `id()`, so a *bare* pickled `Instrument` still breaks (a value-based fingerprint in `transform.py`; `Dataset.__setstate__` works around it); torch `noise.py`'s `_check_one_device` partly redundant with the kernel-as-capability-part check; the "no uncertainty at all" branch of both backends' `_check_gp_uncertainty` is unreachable; upstream sbi 0.27 issues worth filing — `PermutationInvariantEmbedding` reads the first batch element's valid-row count for the whole batch, and `TransformerEmbedding` drops `attention_mask` unless causal; `z_score_x="none"` does not silence sbi's constant-column warning; the Phase 6 multi-observation composition tutorial (plan §5); the `tests/examples` smoke test of photometry + spectra composition.

**Operations notes for the orchestrator**: one five-suite gate at a time, always through the `flock`; run long gates detached (`nohup … &`) and wait with a background `until` loop rather than polling in the foreground; **never edit, format or check out source files in the main checkout while a gate runs there** (a half-written module crashed a subprocess test on 2026-09-09; markdown is safe); an agent that `cd`s to the main checkout can edit its files even though the isolation guard blocks git there — tell agents to stay in their worktree; after a lockfile change, `pixi install -e <env>` in the main checkout; the merged-master gate is the one that counts; parallel agents appending decision-log rows conflict at the table's last line — keep both rows, keep the table contiguous; agent worktrees under `.claude/worktrees/` may hold their own `.pixi` (2 GB) — remove a merged agent's worktree with `git worktree remove --force --force`; finished agents linger in the task list until dismissed with `TaskStop`; an agent that "stops to wait" for a gate is resumed with `SendMessage` and told to poll; a session rate limit kills agents mid-flight but their committed work survives — resume them from their branches.


---

Earlier session records (Phase 0–1 review outcomes, the 2026-09-01 to 2026-09-09 handoffs, including the incremental mid-Phase-3 record) are archived verbatim in `docs/handoff-archive.md`.
