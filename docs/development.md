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

## ⚡ Pick up here (2026-09-08 — Phase 2 closing: slice 3 in flight, then the docs sweep, then Phase 3)

**Where the project is.** Phase 2 is functionally complete and **milestone M2 is reached** (W2.10 merged at `956e1ab`): both modern backends implement the frozen contracts with native O(N) GP solvers, the realisation surface (`inference.md` §10a) gives NUTS and VI on torch and jax through `ampere.core.realise`, the latent-GP science bug is fixed (W2.14), diagnostics, the full plotting surface, the pointwise group, the training-set writer, the WStat example and the CI matrix are in. Master gates at `956e1ab`: **dev 1614 / 95 skipped, jax 2012 / 38, torch 2131 / 38**, lint/format clean, pyrefly 0 errors in all three environments. **Nothing is pushed to origin** (Peter's call; the first live CI run is still owed).

**Peter's standing ruling (2026-09-08 afternoon)**: finish slice 3 on both backend tracks, then a documentation sweep (**W2.15**), then **Phase 3 begins in clean sessions**. Items are written in WORK_ITEMS.md.

**Slice 3 is MERGED on both tracks** (`d6d9d01`, `ae1392a`; status rows). Both backends now have per-instance devices with API-level GPU smoke tests (`pixi run gpu`, outside every gate), the complex Gaussian family on the realised path (the circular complex GP refused until Phase 4 implements it in core), leave-one-out terms on both native quasiseparable solvers, and a gradient-free fast path through `realise` (44× on jax). **Phase 2's code is complete.** Merged-master gates: dev 1616/111, lint/format/pyrefly clean; jax and torch re-run at the merge (numbers in the next status commit). **In flight / next**: dispatch **W2.15** (the documentation sweep, Opus); the orchestrator's half (`CLAUDE.md`, this file, `docs/handoff-archive.md`) is done.

**Clean-session procedure for Phase 3**: (1) read this section, then WORK_ITEMS.md's status table (every merged item's row carries its review record and carried findings), then plan §5 Phase 3; (2) Fable drafts the Phase 3 items (SBI layer: the `simulate` batched form and the coordinate–value–mask tensor encoding for embedding nets, the `sbi` package integration on torch and a jax-native route if maturity warrants, trained-artefact caching keyed on the spec hash per plan §7, SBC/coverage as `diagnostics.md` §11's fourth family) with Peter's approval of the breakdown before dispatch; (3) the orchestration policy in `docs/orchestration.md` is unchanged.

**Open confirmations for Peter** (none blocking; each is recorded where it arose): the noise-model papercut from W2.13 fold-in 7 (a native problem must pass its backend's `IndependentNoise`; softenings are a neutral `BACKEND` sentinel or a backend default likelihood — recommendation: leave as ruled); W2.7's across-draw *median* pooling of the permutation test and its `gp_localisation(datasets=None)` = "every GP dataset" reading; W2.8's corner/trace caps (20/40) and unconditional exclusion of prior-rejected draws from corner plots; the decision-log rows written "Peter to ratify": RHMF deferral (W2.7), docs `-W` deferral and pytest-benchmark (W2.11), the Matérn validation note (W2.10); outside the repository, any branch-protection rule naming the old `phase1-suites` job now needs `suites`.

**Owed**: retroactive `gpt-5.6-terra` passes when the Codex quota returns (~2026-09-30), in this order — W2.3 (GP algebra), W2.1, W2.6, W2.2, then W2.13, W2.4/W2.5 slice 2 (the celerite2 transcriptions), W2.14, W2.10; branch-triage approval (`docs/design/harvest/branch_triage.md`); the push to origin and first live CI run; deletion of the twenty auto-generated `worktree-agent-*` branches (harmless; Peter's action).

**Operations notes for the orchestrator**: this machine has 13 GB — never run two five-suite gates concurrently (the harness kills background tasks on low memory); run torch's `test-all` detached (`nohup … &`) and watch the log, since its pyro NUTS rows push it to ~9 minutes; after a lockfile change, `pixi install -e <env>` in the main checkout; the merged-master gate is the one that counts, re-run after every merge; parallel agents appending decision-log rows conflict at the table's last line — keep both rows, keep the table contiguous; give every shared file exactly one owner when dispatching in parallel.


---

Earlier session records (Phase 0–1 review outcomes, the 2026-09-01 to 2026-09-08 handoffs) are archived verbatim in `docs/handoff-archive.md`.
