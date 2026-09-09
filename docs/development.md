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

## ⚡ Pick up here — Phase 3 APPROVED; W3.0 and W3.9 MERGED; W3.1 slice 1 dispatched 2026-09-09

**Where the project is.** Phases 0–2 are complete (spec frozen at `spec-v1.0`; three backends; M2 reached; W2.15 docs sweep merged). Master is `09357b4`+ and clean apart from the untracked paper script under `examples/examples_paper/` (leave it alone). Gates on the final Phase 2 code (`227d46a`): **dev 1616 passed / 111 skipped, jax 2050 / 42, torch 2222 / 40**, lint/format clean, pyrefly 0 errors ×3; docs build succeeds with 13 legacy-only warnings; characterisation passes. **Nothing is pushed to origin** (origin/master is still at the W1.3 handoff, `b8e585b`; the first live CI run is owed on the whole of Phase 2). The `sbi` pixi environment was installed 2026-09-08 and resolves **sbi 0.27.0 + torch 2.13 (CPU)**.

**What happened 2026-09-08/09 (Fable).** The Phase 3 breakdown W3.0–W3.10 is drafted, approved and recorded (WORK_ITEMS.md "Phase 3" — read it top to bottom; it carries the sbi 0.27 API and environment facts). **W3.0 merged** at `f0140dc` (dev gate on merged master 1616/111; branch gates jax 2052/42, torch 2222/40) and **W3.9 merged** at `fcdd07d` (docs build 13 legacy warnings, unchanged). Master is `fcdd07d`+. **W3.1 slice 1 is dispatched** to an Opus agent on branch `w3.1-simulate-many` from master, in an agent worktree, with the instruction to wait for any running gate before starting its own. If this session was cut off: `git branch --list 'w3.1*'`, `git worktree list`; a branch with commits is reviewed against the slice-1 Accept line and merged `--no-ff`; otherwise re-dispatch from the item text.

**Approval and the W3.1 revision (Peter, 2026-09-08).** The breakdown is **approved**. Peter's three notes on W3.1 are recorded as a decision-log row and W3.1 is now two slices: slice 1 (numpy path) builds `simulate_many` on an order-preserving, partition-independent **executor protocol** with chunking, makes **external compiled simulators** the first-class case (process pool, timeouts, crash capture), and gives `BATCHABLE` a reference-backend meaning; slice 2 lands per-chunk `vmap` prediction and **native observation sampling on torch and jax** (numpy path the oracle, compared distributionally). Rows 85 and 90 ratified. W3.10 (automatic paging above the plot caps, loud warning; not urgent) added. **Dispatched**: W3.9 (Sonnet, docs-only, branch `w3.9-migration-page`). **Held**: W3.1 slice 1 — Fable asked Peter to glance at the revised text before an Opus session is spent on it; dispatch on his nod, from master, with the instruction to wait for any running gate (`pgrep -f pytest`) before starting its own.

**Horizon notes (2026-09-09)**: Peter's Phase 5+ look-ahead is in `docs/design/horizon_notes.md`; two things came out of it into the present — the **observation-context** reservation in the Phase 3 section (confirmed) and a **kernel-algebra item added to the plan's Phase 4** bullets (ruled).

**For Peter — what is still open.** Only the repository actions he alone takes: push master to origin and watch the first live CI run (a branch-protection rule naming the old `phase1-suites` job needs `suites`); delete the stray local `worktree-agent-*` branches; run `pixi run gpu` in the `torch`/`jax` environments on a machine with an accelerator. Every ruling from the 2026-09-08 digest is recorded; W2.9's coverage study is scheduled inside W3.6.

**Owed (unchanged)**: retroactive `gpt-5.6-terra` passes when the Codex quota returns (~2026-09-30), in this order — W2.3 (GP algebra), W2.1, W2.6, W2.2, then W2.13, W2.4/W2.5 slice 2 (the celerite2 transcriptions), W2.14, W2.10, W2.5 slice 3 (the LOO scan); Phase 3 adds W3.1, W3.2, W3.6 to the queue as they merge.

**Next for the orchestrator**: review and merge W3.1 slice 1 when its report arrives (merged-master gates: dev, then jax and torch one at a time); after slice 1 merges, W3.2 (Opus) ∥ W3.1 slice 2 (Opus) — disjoint files; then W3.3 ∥ W3.5, W3.8 after slice 1; then W3.6; W3.4's maturity note any time after W3.2; W3.7 last; W3.10 whenever convenient. Policy unchanged (`docs/orchestration.md`): Opus for judgement-within-spec, Sonnet for well-specified items, Fable for contracts and reviews. One owner per shared file when dispatching in parallel; one five-suite gate at a time.

**Operations notes for the orchestrator**: this machine has 13 GB — never run two five-suite gates concurrently (the harness kills background tasks on low memory); run torch's `test-all` detached (`nohup … &`) and watch the log; after a lockfile change, `pixi install -e <env>` in the main checkout; the merged-master gate is the one that counts, re-run after every merge; parallel agents appending decision-log rows conflict at the table's last line — keep both rows, keep the table contiguous; the two-track pattern works when each shared file has one owner and the sibling's additions are appended; the `sbi` gate (`pixi run -e sbi test-all`) joins the rotation from W3.2.


---

Earlier session records (Phase 0–1 review outcomes, the 2026-09-01 to 2026-09-08 handoffs) are archived verbatim in `docs/handoff-archive.md`.
