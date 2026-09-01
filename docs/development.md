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

## Review checklist (per PR)

- Acceptance criteria evidenced in the report, not just claimed.
- Diff stays inside the item's scope; frozen legacy modules untouched
  (`ampere/data`, `ampere/models`, `ampere/infer`) unless the item says so.
- No run outputs, binaries, or checkpoints; British English in prose.
- From Phase 1 on: conformance suite green; contract changes carry a
  decision-log entry.
- On merge: update the WORK_ITEMS.md status table.

## Environment

- Conda env `ampere` (Python 3.13) has an editable install pointing at the
  main checkout — worktree-based work that needs importing its own changes
  should `pip install -e .` into a fresh env instead.
- Fresh setup: `pip install -e ".[dev]"`, Python ≥ 3.11. (Interim: once
  W0.8 lands, pixi becomes the supported route — `pixi install` /
  `pixi run test` — and environment.yml is retired.)
- `pytest tests/characterisation` is the "legacy still works" gate once
  W0.4 lands; run it before merging anything that touches shared files.

## Current state (2026-09-01)

- Phase 0 in progress: W0.1 done (committed directly to master — the fix
  existed only in the local working tree); W0.2, W0.3 dispatched.
- The AMPERE paper revision proceeds on the legacy code and takes priority
  in any conflict over `examples/examples_paper/`.
