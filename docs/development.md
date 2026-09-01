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
    `pixi run test-characterisation` (the full suite), `pixi run lint`,
    `pixi run format-check`, `pixi run typecheck`, `pixi run docs`.
  - Other environments (`pixi run -e <env> <task>`): `test-py311` /
    `test-py312` / `test-py313` (the CI matrix, one Python each); `sbi`
    (adds the `sbi` extra — torch is a large download, so it is not part
    of `dev` or the `test-py3*` environments); `torch` / `jax` (Phase 2
    backend placeholders, not exercised by anything yet).
- Plain-pip alternative: `pip install -e ".[dev]"`, Python ≥ 3.11.
- Conda env `ampere` (Python 3.13) has an editable install pointing at the
  main checkout — worktree-based work that needs importing its own changes
  should use pixi or `pip install -e .` into a fresh env instead.
- `pytest tests/characterisation` (or `pixi run test-characterisation`) is
  the "legacy still works" gate; run it before merging anything that
  touches shared files.

## Phase 1: W1.1–W1.4 merged; W1.9 + W1.12 on Peter's review list (2026-09-01)

The Fable review the previous handoff asked for is done and, on Peter's
instruction, the reviewed branches were merged to local master
(W1.1 → W1.2 → W1.3, then W1.4), keeping master's `WORK_ITEMS.md` and
`docs/development.md` where branches carried stale snapshots. All gates
re-verified on the merged master after each merge (247 core tests as of
W1.4, the fast import suite, pyrefly, ruff lint and format). Nothing has
been pushed to origin. One correction to the earlier handoff text: W1.2 was
Sonnet-authored (not Fable, as previously recorded) and was therefore given
a full contract-tier review rather than a light reconciliation.

**Adversarial sol review is batched at the W1.13 freeze** (Peter's call,
2026-09-01): one pass over the §4 contract code (W1.3 parameters, W1.4
results schema, plus whatever lands by then) and the lowering spec (W1.9),
per `docs/orchestration.md` principle 3, rather than per-item passes now.

**Review outcomes (branches now in master's history):**

- **W1.1** (`w1.1-prior-art-memo`, Sonnet): merged as authored — no
  changes. Its two most load-bearing claims were re-verified against live
  sources during review: bilby's current core really does take an explicit
  `log_likelihood(parameters)` argument (checked against `bilby-dev/bilby`
  `main`), and gammapy issue #2859 / fix PR #2861 are as described. The
  memo's own §7 confidence table is accurate; its weakest section (3ML,
  tutorial-level sourcing) is honestly flagged there.
- **W1.2** (`w1.2-architecture-spec`, Sonnet + Fable amendments): one
  amendment commit on top of the draft. Core's dependency floor corrected
  to numpy/scipy/**astropy**/stdlib (the plan requires units on parameters;
  astropy is a required base dependency; this resolves W1.3's open
  question 1). §10 rewritten as the W1.1 reconciliation record — none of
  the four questions it posed required structural change. One genuine
  conflict with the plan was found (curated astropy→native translation:
  opt-in vs the plan's "silently restoring differentiability") and **ruled
  by Peter on 2026-09-01: opt-in, never silent** — recorded in
  `DEVELOPMENT_PLAN.md` §2's decision table with §4.7 amended, and
  architecture.md §1 updated to match (sketch: a backend-scoped
  `from_astropy()` constructor as the explicit opt-in, raising rather than
  silently falling back when a model is not fully translatable).
- **W1.3** (`w1.3-parameter-contract`, Opus + Fable fix commit): strong
  contract; three defects found and fixed on the branch, with tests
  (112 pass; pyrefly/ruff/format clean): (1) tying two hierarchical-prior
  parameters whose hyperparameters were not themselves tied silently wired
  the collapsed parameter to the first component's hyperparameters — now a
  `TyingError` with a hint; (2) `describe_prior` kept scipy's positional/
  keyword split, so `st.norm(1.5, 0.1)` and `st.norm(loc=1.5, scale=0.1)`
  counted as different priors when tied and serialised differently — the
  description is now canonical (keyword-only), which also gives W1.9 a
  single form to lower; (3) conflicting explicit bijections on tied sites
  silently took the first — now a `TyingError`. Checked against W1.1
  retrospectively: its declaration-based tying is exactly what memo lesson
  G2 prescribes.

- **W1.4** (`w1.4-results-schema`, Opus + Fable fix commit): merged
  2026-09-01. Strong contract; the agent's own self-review caught three
  unit-handling bugs, and the Fable pass found and fixed two more:
  (1) `from_unsorted` permuted the axis/values/uncertainty/mask but left
  `extra_coords` in declaration order — silent per-sample misalignment,
  now permuted with the data; (2) `to_unit` leaked astropy's raw
  `UnitConversionError` where every sibling path raises `SchemaError` —
  now wrapped. 247 core tests; W1.3 untouched.

**Peter's review list** (none block further Phase 1 work; the first two
gate their branches' merges):

- **W1.9** (`w1.9-lowering-spec`, Opus, Fable-reviewed, **not merged**):
  needs two decision-log rulings before merge — (a) jax non-trainable
  mechanism: the spec ranks `eqx.partition` above `paramax.NonTrainable`
  (whose freezing depends on `unwrap()` being called — forgetting it
  silently trains the buffers), narrowing plan §4.1's unranked wording;
  (b) x64 activation: guard-and-raise at construction instead of
  set-on-import, amending `architecture.md` §5's sketch (the policy —
  x64 always — is unchanged). The Fable review endorses both. The spec's
  §12 lists five further ratification items (most route to W1.13).
- **W1.12** (`w1.12-diagnostics-spec`, Sonnet, Fable-reviewed, **not
  merged**): the acceptance criterion is Peter's review pass. Key
  decisions to check: `ampere.diagnostics` as a new peer namespace behind
  a `diagnostics` extra; the shared `AnomalyScore` container with
  mandatory provenance metadata (Tension 5's resolution); RHMF
  hyperparameters as expert opt-in. Its §10 lists the open questions.
- W1.3 spec §14's remaining open questions — tie labels as a flat global
  namespace; whether a lone `shared_as` should raise; ratifying
  `OptionalDependencyError`'s shape at W1.13; confirming nothing needs
  the legacy `npars` alias. Review at
  `docs/design/contracts/parameters.md` §14, with §8 (tying semantics) and
  §12 (deliberate limitations) as the supporting context, and
  `ampere/core/parameter.py` (`_merge`/`_collapse`) /
  `ampere/core/exceptions.py` as the implementation. (Question 1,
  astropy-in-core, was resolved by the W1.2 amendment; question 2,
  recursive merge, is explicitly deferred to W1.7.)
- W1.4 spec §17's open questions — whether `"default"` should be a
  reserved channel name; overlapping échelle orders as two channels
  (duplicate coordinates refused); `Cube` axis order `(x, y, spectral)`
  vs FITS convention (cheap to change now); whether a `ModelResult`
  should carry the θ that produced it (emulator training sets). Review at
  `docs/design/contracts/results_schema.md` §17.

**Tooling note** (found by the W1.4 agent): `ruff format` invoked with an
explicit path bypasses `extend-exclude`, so it can silently rewrite the
frozen W0.7 harvest snapshots under `docs/design/harvest/`. Worth a guard
in a small W0.x follow-up.

**Obligations W1.3 places on later specs** are recorded in the spec's own
§13 (`docs/design/contracts/parameters.md`) — W1.4–W1.10 authors read that
section before starting; dispatch prompts should cite it.

A `SendFeedback` draft was queued in the downgraded session about the
`/model` switch-back failure (usage-credit downgrade not reversible via
`/model` even after `/login`) — the user can review and send it with
`/feedback` if they want to report it.

## Current state (2026-09-01)

- **Phase 0: complete.** W0.1–W0.8 merged to master; CI green on the first
  live run. Issues #74–77 closed. Branch archival executed (15 branches
  tagged `archive/*`; `jax` re-pushed filtered to drop ~50 MB of
  checkpoints; `small_silicates` kept live per Peter). W0.9 (pyphot ≥2 /
  current-sbi forward migration) not started — the temporary
  `pyphot<2`/`sbi<0.28` pins are in place and documented in
  `DEVELOPMENT_PLAN.md` §2.
- **Phase 1: in progress**, see the handoff section above for exact
  branch/review state of W1.1–W1.3. W1.4 onward not started.
- The AMPERE paper revision proceeds on the legacy code and takes priority
  in any conflict over `examples/examples_paper/`.
