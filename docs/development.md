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

## ⚡ Pick up here — PHASE 2 CLOSED 2026-09-08; Phase 3 begins in a clean session

**Where the project is.** Phases 0–2 are complete. The core contracts are frozen (`spec-v1.0`); both modern backends implement them alongside the numpy reference; milestone M2 is reached (W2.10, the flagship misspecification validation on all three backends, with the science claim asserted as thresholds); slice 3 (devices, complex Gaussian, LOO parity, the gradient-free fast path) and the documentation sweep (W2.15) are merged. **Master is `fbd6970`+**, clean. Gates on the final Phase 2 code (`227d46a`, unchanged by the docs-only W2.15 apart from one spelling of a constant): **dev 1616 passed / 111 skipped, jax 2050 / 42, torch 2222 / 40**, lint/format clean, pyrefly 0 errors in all three environments; the docs build succeeds with 13 legacy-only warnings; the characterisation suite passes. **Nothing is pushed to origin** (Peter's call — the first live CI run is still owed).

**Peter's standing ruling (2026-09-08)**: Phase 3 begins in clean sessions. The procedure: (1) read this section, then WORK_ITEMS.md's status table (every merged item's row is the authoritative record of what landed and what it carried), then `DEVELOPMENT_PLAN.md` §5 Phase 3 and §2's decision log; (2) **Fable drafts the Phase 3 items** into WORK_ITEMS.md — the SBI layer: `simulate`'s batched form and the coordinate–value–mask tensor encoding for embedding networks (`inference.md` §13, limitation 17.5), the `sbi` package integration on torch and a jax-native route only if maturity warrants (plan §6), trained-artefact caching keyed on the spec hash (plan §7's trap), the training-set writer's unlimited-dimension append (`results.md` §13.9) when budgets need it, and SBC/coverage as `diagnostics.md` §11's fourth family — and gets **Peter's approval of the breakdown before dispatch**; (3) `docs/orchestration.md`'s policy is unchanged (Opus for judgement-within-spec, Sonnet for well-specified items, Fable for contracts and reviews; cross-model review via terra when the quota returns); (4) give every shared file exactly one owner when dispatching in parallel, and run five-suite gates one at a time.

**Small housekeeping for the first Phase 3 session** (one Sonnet item or done by hand): `pip install ampere` is an unrelated PyPI package — fix `OptionalDependencyError`'s remedy text (`ampere/core/exceptions.py`) and the six docstrings that quote it (W2.15's finding; the install page already warns); `SAMPLE_STATS_GROUP` is defined twice in `ampere.results`; rename `LSFConvolution.sigma_tensor` (torch) so it cannot be mistaken for the noise hook; add a native zero-uncertainty precondition to the jax realised GP path (the `realise` guard catches it today).

**Open confirmations for Peter** (none blocking; each recorded where it arose): the noise-model papercut from W2.13 fold-in 7 (recommendation: leave as ruled); W2.7's median pooling and `datasets=None` reading; W2.8's corner/trace caps and the exclusion of prior-rejected draws from corner plots; the "Peter to ratify" decision-log rows — RHMF deferral (W2.7), pytest-benchmark and the docs `-W` deferral (W2.11), the Matérn validation note (W2.10), the W2.15 annotation batch; `use_realisation=True` on torch gains nothing (W2.5 slice 3) and `ampere_realised = 1` on gradient-free realised runs is a visible provenance change; `:no-index:` on the legacy API pages; outside the repository, any branch-protection rule naming the old `phase1-suites` job now needs `suites`.

**Owed**: retroactive `gpt-5.6-terra` passes when the Codex quota returns (~2026-09-30), in this order — W2.3 (GP algebra), W2.1, W2.6, W2.2, then W2.13, W2.4/W2.5 slice 2 (the celerite2 transcriptions), W2.14, W2.10, W2.5 slice 3 (the LOO scan); branch-triage approval (`docs/design/harvest/branch_triage.md`); the push to origin and first live CI run; deletion of the auto-generated `worktree-agent-*` branches (harmless; Peter's action); GPU smoke tests actually run on hardware (`pixi run gpu` in the `torch`/`jax` environments on a machine with an accelerator).

**Operations notes for the orchestrator**: this machine has 13 GB — never run two five-suite gates concurrently (the harness kills background tasks on low memory); run torch's `test-all` detached (`nohup … &`) and watch the log with a Monitor, since it takes ~9–11 minutes; after a lockfile change, `pixi install -e <env>` in the main checkout; the merged-master gate is the one that counts, re-run after every merge; parallel agents appending decision-log rows conflict at the table's last line — keep both rows, keep the table contiguous; the two-track pattern (torch/jax) works when each shared file has one owner and the sibling's additions are appended.


---

Earlier session records (Phase 0–1 review outcomes, the 2026-09-01 to 2026-09-08 handoffs) are archived verbatim in `docs/handoff-archive.md`.
