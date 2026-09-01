# Ampere — agent guide

Ampere is a Bayesian fitting environment for heterogeneous astronomical data
(SEDs, spectra, and more), whose distinguishing feature is a flexible,
GP-based likelihood providing robustness to model misspecification.

The project is mid-way through a major redesign ("v2"). Before any
non-trivial work, read:

- **`DEVELOPMENT_PLAN.md`** — the source of truth: decisions taken, target
  architecture (backend-neutral core + reference/torch/jax backends), phased
  roadmap, and known traps (§7 — read it, they are easy to fall into).
- **`WORK_ITEMS.md`** — the agent-sized work items for the current phases,
  with the working agreement and per-item acceptance criteria.

## Ground rules

1. **Legacy is frozen.** Do not modify `ampere/data/`, `ampere/models/`, or
   `ampere/infer/` unless your work item explicitly says so.
2. **One work item = one branch = one PR.** Branch naming: `w0.2-repo-hygiene`
   style (item id + slug). Peter Scicluna reviews and merges everything.
3. **Never** push to `master`, push tags, delete branches, or publish
   anything externally. Do not push at all unless your dispatcher says to.
4. **Stay in scope.** Do only what the work item says. If you find other
   problems, record them in your report/PR description — do not fix them.
5. Acceptance criteria are the definition of done. Run them and report the
   evidence; never claim completion without it.
6. **British English** in all documentation and prose.
7. No binary artefacts, run outputs, logs, or checkpoints in git.
8. New namespaces (`ampere/core`, `ampere/backends`, `ampere/inference`,
   `ampere/results`) are typed from the first line (pyrefly); legacy is
   exempt and never gets typed.
9. After the Phase 1 spec freeze, any change to a §4 contract requires a
   decision-log entry in `DEVELOPMENT_PLAN.md` in the same PR.
10. End commit messages with:
    `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>` (or the
    equivalent for your agent/tool).

## Environment

- **pixi is the supported route** (W0.8): `pixi install -e dev` sets up the
  daily-use environment (Python 3.13 + `dev`, `zeus`, `extinction` extras)
  from a clean clone with only `pixi` installed; pyproject.toml's
  `[tool.pixi.*]` tables wrap `[project.dependencies]` /
  `[project.optional-dependencies]` rather than duplicating them —
  `pixi.lock` is committed, `.pixi/` is not. Common tasks: `pixi run test`
  (fast: `tests/test_imports.py`), `pixi run test-characterisation` (full
  suite), `pixi run lint`, `pixi run format-check`, `pixi run typecheck`,
  `pixi run docs`. Other environments: `test-py311`/`test-py312`/
  `test-py313` (the CI matrix), `sbi` (adds the `sbi` extra — torch is a
  large download, not part of `dev`), `torch`/`jax` (backend placeholders).
  Select one with `pixi run -e <env> <task>`.
- Plain-pip alternative: `pip install -e ".[dev]"` (Python ≥ 3.11; CI
  targets 3.11–3.13).
- Conda env **`ampere`** (`~/miniforge3/envs/ampere`, Python 3.13):
  `conda run -n ampere python ...`. Caution: that env has an editable
  install pointing at the main checkout — if you are working in a worktree
  and need to *import* your changed code, use pixi or `pip install -e .`
  into a fresh env/venv from your worktree instead.
- Tests: `pytest tests/` (characterisation suite: `pytest
  tests/characterisation`, slower; it defines "legacy still works") — or
  the equivalent pixi tasks above.
- Lint/format: `ruff check .` / `ruff format` (or `pixi run lint` /
  `pixi run format-check`). Types: `pyrefly check` (scoped to the new
  namespaces), or `pixi run typecheck`.

## Workflow for a work item

See the `work-item` skill (`.claude/skills/work-item/SKILL.md`) — it is the
step-by-step procedure. Codex agents: follow that file manually; it is
short.

## Delegation & model selection

`docs/orchestration.md` is the policy for which model tier gets which
work — read it before dispatching anything. In short: Sonnet for
well-specified items, Opus for judgement-within-spec, Haiku/luna for
mechanical bulk, Fable for contracts and integration, gpt-5.6-sol for
cross-model reviews.

Claude agents may delegate well-scoped, self-contained subtasks to the
OpenAI Codex CLI (installed: `codex`). See
`.claude/skills/delegate-codex/SKILL.md` (use `-m gpt-5.6-luna` for cheap
mechanical work), and `.claude/skills/sol-review/SKILL.md` for detailed or
adversarial reviews with gpt-5.6-sol. Codex reads this AGENTS.md
automatically. Always review Codex's diff before committing it; Codex never
pushes.

## Repository map

- `DEVELOPMENT_PLAN.md`, `WORK_ITEMS.md` — plan + items (see above).
- `docs/development.md` — human-facing onboarding/handoff notes.
- `docs/design/` — Phase 1 output: architecture spec, per-contract specs,
  lowering rules, modality sketches, prior-art memo, harvest of old
  branches.
- `ampere/` — the package. Legacy: `data/`, `models/`, `infer/`, `utils/`.
  New (as phases land): `core/`, `backends/{reference,torch,jax}/`,
  `inference/`, `results/`.
- `examples/` — legacy examples; `minimal_working_example*.py` are the
  characterisation-test anchors. `examples/examples_paper/` is
  paper-revision work in progress — leave it alone.
