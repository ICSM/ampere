# Ampere — agent guide

Ampere is a Bayesian fitting environment for heterogeneous astronomical data
(SEDs, spectra, and more), whose distinguishing feature is a flexible,
GP-based likelihood providing robustness to model misspecification.

The project is undergoing a major redesign ("v2"). **Phases 0–2 are
complete** (closed 2026-09-08): the core contracts are frozen
(`spec-v1.0`, 2026-09-03), both modern backends (`ampere.backends.torch`,
`.jax`) implement them alongside the numpy reference backend, milestone M2
(the flagship misspecification validation) is reached, and the
documentation reflects what landed. **Phase 3 (the SBI layer) is in
progress** (since 2026-09-08): `simulate_many` with its executor protocol,
the `SBIEngine` over sbi 0.27, the coordinate–value–mask encoding, artefact
caching and the foreign-parts opt-in are merged; the live state, what is in
flight and what comes next are in `docs/development.md`'s "⚡ Pick up
here". Before any non-trivial work, read:

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
9. **In force since the freeze (2026-09-03)**: any change to a §4
   contract requires a decision-log entry in `DEVELOPMENT_PLAN.md` in the
   same PR, and must keep the conformance suite green or update it in the
   same PR with justification.
10. End commit messages with:
    `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>` (or the
    equivalent for your agent/tool).
11. **Handoffs must survive an interruption.** Any session may be cut off
    without warning, so the resumable state lives in the repository, never
    only in a conversation. The orchestrating session refreshes
    `docs/development.md`'s "⚡ Pick up here" section at every state
    change — a merge, a dispatch, a ruling, a blocker — recording what is
    in flight (branches, agents, pending reviews, open rulings) and what a
    clean session should do next; the WORK_ITEMS.md status table is
    updated at merge as before. Working agents contribute by committing
    logical units as they complete them and keeping their branch's final
    report reproducible from the branch alone.

## Environment

- **pixi is the supported route** (W0.8): `pixi install -e dev` sets up the
  daily-use environment (Python 3.13 + `dev`, `zeus`, `extinction`,
  `netcdf` features; pandoc and ipykernel for the docs build) from a clean
  clone with only `pixi` installed; pyproject.toml's `[tool.pixi.*]`
  tables wrap `[project.dependencies]` / `[project.optional-dependencies]`
  rather than duplicating them — `pixi.lock` is committed, `.pixi/` is not.
  **Three real environments**: `dev` (no torch, no jax), `torch` and `jax`
  (each the `dev` feature set plus its backend; torch resolves from
  PyTorch's CPU index in pixi only). The gate is `pixi run -e <env>
  test-all` — core, results, conformance, backends, inference, examples
  and m2 in ONE pytest process — in **each** environment; a backend
  package is typechecked with its real types only in its own environment
  (`pixi run -e torch typecheck`), and excluded from the `dev` typecheck.
  Other tasks: `pixi run test` (fast import sweep), `test-core` /
  `test-results` / `conformance` / `test-inference` / `test-examples` /
  `test-m2` individually, `bench` (pytest-benchmark → `benchmark.json`),
  `scaling`, `test-characterisation` (legacy still works), `lint`,
  `format-check`, `typecheck`, `docs`. Plain `pixi run <task>` equals
  `pixi run -e dev <task>`. `test-py311`/`test-py312`/`test-py313` are the
  CI matrix; `sbi` adds the `sbi` extra. **Run five-suite gates one at a
  time** — two at once exhaust a 13 GB machine; torch's takes ~9 minutes.
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
mechanical bulk, Fable for contracts and integration, and cross-model
adversarial reviews via the Codex CLI — `gpt-5.6-sol` where available,
with **`gpt-5.6-terra` as the ruled stand-in** on this account (sol is
refused on ChatGPT-auth codex; the Fable + terra pairing was proven at
the W1.13 freeze).

Claude agents may delegate well-scoped, self-contained subtasks to the
OpenAI Codex CLI (installed: `codex`). See
`.claude/skills/delegate-codex/SKILL.md` (use `-m gpt-5.6-luna` for cheap
mechanical work), and `.claude/skills/sol-review/SKILL.md` for detailed or
adversarial reviews with gpt-5.6-sol. Codex reads this AGENTS.md
automatically. Always review Codex's diff before committing it; Codex never
pushes.

## Repository map

- `DEVELOPMENT_PLAN.md`, `WORK_ITEMS.md` — plan + items (see above).
- `docs/development.md` — onboarding + the live "⚡ Pick up here" restart
  point; `docs/handoff-archive.md` — superseded session records.
- `examples/m2_misspecification/` — the M2 study (generators, models per
  backend, driver, figures); `examples/wstat_comparison.py` — the WStat
  example; `docs/source/` — the Sphinx site (`pixi run docs`).
- `docs/design/` — Phase 1 output, **frozen at `spec-v1.0`**: architecture
  spec, per-contract specs (`contracts/`), lowering rules, the
  serialisation review, modality sketches, prior-art memo, harvest of old
  branches.
- `ampere/` — the package. Legacy (frozen): `data/`, `models/`, `infer/`,
  `utils/`. v2: `core/` (the frozen contracts, implemented — parameters,
  containers, transformations, likelihoods with `DenseGP`/`QuasisepGP`,
  datasets and `FittingProblem`, the lowering and realisation registries),
  `backends/{reference,torch,jax}/` (one name per backend everywhere; a
  native problem is composed entirely from one backend's pieces, noise
  models and solvers included), `inference/` (emcee/dynesty/zeus on the
  numpy path; NUTS and VI through `ampere.core.realise`; imports no
  backend), `results/` (ArviZ `DataTree` emission with provenance schema
  5, derived groups, diagnostics, the six plots, training sets).
  `ampere.diagnostics` (RHMF pre-fit screening) is deferred, not landed.
- `tests/` — `core`, `results`, `conformance` (the lockstep battery: one
  column per registered backend fixture, `backends/__init__.py` is the one
  place a backend is named), `backends`, `inference`, `examples`, `m2`
  (the M2 study's assertions), `benchmarks`, `scaling`, `gpu` (API-level
  smoke tests, skipped without an accelerator), `characterisation`
  (legacy).
- `examples/` — legacy examples; `minimal_working_example*.py` are the
  characterisation-test anchors. `examples/examples_paper/` is
  paper-revision work in progress — leave it alone.
