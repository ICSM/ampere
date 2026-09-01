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

- **W1.9**: *both rulings approved by Peter and merged 2026-09-01* —
  `eqx.partition` over `paramax.NonTrainable`, and x64 guard-and-raise —
  recorded in the plan's §2 decision table, §4.1 and `architecture.md` §5
  amended to match. Three follow-ups from Peter's review, all recorded:
  dependent priors (y ~ norm(x, 1)) are already first-class via
  `HierarchicalPrior` (verified live against the merged contract —
  `parameters.md` §9); `numpyro.enable_x64()` was source-verified as a
  thin wrapper over the same jax flag, so the construction guard covers
  numpyro with no separate switch (lowering.md §12.5); and
  backend-specific lowerings for user-defined priors/bijections are now
  lowering.md §12 item 8 (registration-hook sketch, W1.13 decides the
  principle, Phase 2 builds the plumbing). Remaining §12 items route to
  W1.13.
- **W1.12**: *approved by Peter and merged 2026-09-01* — the §10 defaults
  (peer `diagnostics` namespace, `AnomalyScore` in core, RHMF expert
  opt-in) stand as written. His review added one thing, recorded as the
  spec's §11: posterior calibration (SBC, expected-coverage/TARP) as a
  fourth, future diagnostic family — it validates the *inference* rather
  than the model, matters most for Phase 3's SBI layer, and needs no new
  contract surface (§4.5's `simulate()` is the hook). §11 also records
  the extension template for future families, since the area will evolve.
- W1.3 spec §14: *ruled by Peter 2026-09-01* — recommendations stand (flat
  tie labels; lone `shared_as` allowed; `npars` dropped;
  `OptionalDependencyError` to W1.13), **except that recursive merge is
  expressly kept open**: nested/hierarchical merging may be the natural
  approach for `DatasetCollection`, so W1.7 must evaluate a nested
  `ParameterMapping` as a first-class design option (recorded in
  `parameters.md` §14's preamble; W1.5 §14's ties-across-levels note is
  part of the same question). **W1.7's dispatch prompt must carry this.**
- W1.4 spec §17: *ruled by Peter 2026-09-01* — `"default"` stays
  unreserved (clash accepted, loudly documented); overlapping échelle
  orders as two channels confirmed; and θ-carrying results implemented:
  `ModelResult.parameters` exists and `Model.__call__` attaches the
  resolved values automatically (`65c997f`). His follow-up — nested
  result channels — assessed as feasible via qualified flat names
  (mirroring `ParameterSet.merge`) and **routed to W1.7 together with the
  symmetrical nested-merge question**; see `results_schema.md` §17's
  preamble. Questions 3–6 (`Cube` axis order, `extra_coords` units, axis
  naming, serialisation) remain open as written.
- W1.5 spec §15's open questions (merged 2026-09-01) — most notably:
  whether the ANY mask-propagation rule is too conservative for real
  resampling (`min_valid_fraction` is the named extension); whether a
  model should be able to *refuse* a requirement rather than silently
  ignoring it; confirming the `Model` ABC's placement in W1.5 rather than
  W1.7 (§9 has the argument); and §15.8's `max_step`/`min_resolving_power`
  dual meaning (declaration vs post-union summary — the review leans
  towards demoting the scalars before the freeze, cheap now). Review at
  `docs/design/contracts/transformations.md` §15.

**Sol review (status 2026-09-01)**: `codex exec -m gpt-5.6-sol` is
currently refused — "The 'gpt-5.6-sol' model is not supported when using
Codex with a ChatGPT account" (reproduced independently of the W1.5
agent's report). Peter is aware and expects to configure the extra codex
steps to make sol available — one more reason the adversarial review is
batched at the W1.13 freeze rather than run per item. No Phase 1 work
blocks on it before then.

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
