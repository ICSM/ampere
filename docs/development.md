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
    (adds the `sbi` extra); `torch` / `jax` (the backend environments —
    real gates since W2.4/W2.5: `pixi run -e torch test-all`,
    `pixi run -e jax test-all`, and each backend package is typechecked
    with its real types only in its own environment). Since W2.11 both
    `torch` and `sbi` resolve torch from PyTorch's **CPU** index, so
    neither pulls ~2 GB of CUDA runtime onto a CPU-only runner; this
    changes nothing about what `pip install "ampere[torch]"` gives a user.
    **Run five-suite gates one at a time** — two concurrently exhaust a
    13 GB machine.
- Plain-pip alternative: `pip install -e ".[dev]"`, Python ≥ 3.11.
- Conda env `ampere` (Python 3.13) has an editable install pointing at the
  main checkout — worktree-based work that needs importing its own changes
  should use pixi or `pip install -e .` into a fresh env instead.
- `pytest tests/characterisation` (or `pixi run test-characterisation`) is
  the "legacy still works" gate; run it before merging anything that
  touches shared files.

## ⚡ Pick up here — Phase 3 items DRAFTED 2026-09-08 (W3.0–W3.9), rulings recorded, awaiting Peter's approval; W3.0 in flight

**Where the project is.** Phases 0–2 are complete (spec frozen at `spec-v1.0`; three backends; M2 reached; W2.15 docs sweep merged). Master is `09357b4`+ and clean apart from the untracked paper script under `examples/examples_paper/` (leave it alone). Gates on the final Phase 2 code (`227d46a`): **dev 1616 passed / 111 skipped, jax 2050 / 42, torch 2222 / 40**, lint/format clean, pyrefly 0 errors ×3; docs build succeeds with 13 legacy-only warnings; characterisation passes. **Nothing is pushed to origin** (origin/master is still at the W1.3 handoff, `b8e585b`; the first live CI run is owed on the whole of Phase 2). The `sbi` pixi environment was installed 2026-09-08 and resolves **sbi 0.27.0 + torch 2.13 (CPU)**.

**What happened this session (2026-09-08, Fable).** (1) The Phase 3 breakdown is drafted into WORK_ITEMS.md as **W3.0–W3.7** plus a "deferred" list and a dispatch-order paragraph — read that section top to bottom; it carries the sbi 0.27 API facts and the environment facts every item needs. **Nothing beyond W3.0 may be dispatched until Peter approves the breakdown** (his standing ruling). (2) **W3.0** (the four housekeeping findings the previous handoff sanctioned as "one Sonnet item or by hand") is **dispatched to a Sonnet agent** on branch `w3.0-housekeeping` from `09357b4` in an agent worktree; it runs all three gates itself. If this session was cut off before the merge: check `git branch --list 'w3.0*'` and `git worktree list`; if the branch exists with commits, review it against W3.0's Accept line (the report is reproducible from the branch) and merge with `--no-ff`; if it does not, re-dispatch from the item text.

**Rulings recorded 2026-09-08 (Peter, on Fable's digest of 26 pending decisions).** Everything ratified or confirmed per the recommendations; the decision-log rows at §2 lines 91–98 now read "ratified by Peter 2026-09-08", the status rows' "Open for Peter" clauses are closed with the ruling, W0.7's stale branch-triage wording is struck. Four rulings carried consequences: (1) RHMF stays deferred but **early testing is wanted in a later phase** — a Phase 5 bullet in the plan's §5; (2) the corner/trace caps stand, with `var_names=`/`max_variables=` as the routes past them (no automatic paging); (3) a core numpy `Kernel` inside a native GP noise model is **refused by default, accepted under an explicit opt-in when no gradient is needed, always refused where a gradient is required** — decision-log row recorded, implementation drafted as **W3.8**; (4) the legacy API pages stay linked from their own section and gain a **"Migrating to the new API"** page — drafted as **W3.9** (Sonnet, docs-only), the full guide remaining Phase 6's.

**For Peter — what is still open.**
1. **Approve, amend or reject the Phase 3 breakdown**, now W3.0–W3.9 (WORK_ITEMS.md "Phase 3"). The load-bearing choices: W3.1 keeps the observation draw on the numpy path and vectorises only the noise-free prediction natively; W3.2 builds sbi's prior in the *unconstrained* coordinates and scores every stored draw on the numpy path so runs carry the true `log_prob`; W3.3 is a new §4 contract (`encoding.md`) with a layout hash; W3.4 (swyft) is gated on a maturity note; W3.5's cache-key ingredients; W3.6's placement in `ampere.results.calibration`; W3.8's opt-in spelling; the unlimited-dimension append deferred until a budget outgrows memory.
2. **Two decision-log rows the digest missed** still say "Peter to ratify": §2 line 85 (backend identity on `FittingProblem`, W2.12) and line 90 (celerite2.jax chosen over tinygp by measurement, W2.5 slice 2). Both were merged after Fable review; a one-word ratification closes them.
3. **Repository actions only Peter takes**: push master to origin and watch the first live CI run (a branch-protection rule naming the old `phase1-suites` job needs `suites`); delete the two stray local `worktree-agent-*` branches; run `pixi run gpu` in the `torch`/`jax` environments on a machine with an accelerator.

**Owed (unchanged)**: retroactive `gpt-5.6-terra` passes when the Codex quota returns (~2026-09-30), in this order — W2.3 (GP algebra), W2.1, W2.6, W2.2, then W2.13, W2.4/W2.5 slice 2 (the celerite2 transcriptions), W2.14, W2.10, W2.5 slice 3 (the LOO scan); Phase 3 adds W3.1, W3.2, W3.6 to the queue as they merge.

**Next for the orchestrator, once Peter approves**: merge W3.0; dispatch W3.1 (Opus) alone or in parallel with nothing else that touches `ampere/core/dataset.py`; after W3.1 merges, W3.2 (Opus) alone; then W3.3 ∥ W3.5; W3.8 after W3.1; W3.9 (docs-only) any time; then W3.6; W3.4's maturity note any time after W3.2; W3.7 last. Policy unchanged (`docs/orchestration.md`): Opus for judgement-within-spec, Sonnet for well-specified items, Fable for contracts and reviews. One owner per shared file when dispatching in parallel; one five-suite gate at a time.

**Operations notes for the orchestrator**: this machine has 13 GB — never run two five-suite gates concurrently (the harness kills background tasks on low memory); run torch's `test-all` detached (`nohup … &`) and watch the log; after a lockfile change, `pixi install -e <env>` in the main checkout; the merged-master gate is the one that counts, re-run after every merge; parallel agents appending decision-log rows conflict at the table's last line — keep both rows, keep the table contiguous; the two-track pattern works when each shared file has one owner and the sibling's additions are appended; the `sbi` gate (`pixi run -e sbi test-all`) joins the rotation from W3.2.


---

Earlier session records (Phase 0–1 review outcomes, the 2026-09-01 to 2026-09-08 handoffs) are archived verbatim in `docs/handoff-archive.md`.
