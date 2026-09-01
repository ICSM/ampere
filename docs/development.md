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

## ⚠ Session handoff in effect (2026-09-01) — read this first

The orchestrating session hit a Claude usage-credit limit partway through
Phase 1 and was downgraded from Fable 5 to Sonnet 5; `/model` could not
switch it back even after re-login. Per Peter's instruction, from the
downgrade point onward that session stopped doing Fable-tier contract
authorship and contract-soundness review — it only finished mechanical
bookkeeping (staging/rebasing agent branches, updating this file and
`WORK_ITEMS.md`) so a fresh Fable session can pick up cleanly. **A new
Fable session should do the following before dispatching any further
Phase 1 work:**

1. Read `DEVELOPMENT_PLAN.md` and `WORK_ITEMS.md` (as always).
2. Read `docs/design/architecture.md` (**W1.2**) — drafted by Fable earlier
   in the downgraded session, *before* the downgrade, so its content is
   legitimate contract-tier work, not a Sonnet artefact. It is on local
   branch `w1.2-architecture-spec` (**not pushed to origin, not merged**).
   Its own §10 lists exactly what to re-check against W1.1's findings
   (bilby's capability-flag framing, gammapy's `Datasets` pattern, 3ML's
   instrument-owned-likelihood question, Starfish's kernel/cost lessons).
3. Read `docs/design/prior_art.md` (**W1.1**) — Sonnet-authored (on-policy;
   W1.1 was always meant for a cheaper tier), staged and rebased onto
   master on branch `w1.1-prior-art-memo` (**not merged**). It flags one
   correction to a dispatch-prompt error: gammapy's `Datasets`/tying
   material bears on `DEVELOPMENT_PLAN.md` **§4.5** (inference contracts,
   → W1.7), not §4.7 (astropy interop) — no committed doc currently
   contains that error, but keep it in mind when drafting W1.7.
4. Reconcile W1.2 against W1.1 (architecture.md §10), do a real Fable-tier
   review of both, merge what's sound, and only then continue.
5. **W1.3** (parameter & prior contract, Opus) **has completed** and is
   mechanically staged on branch `w1.3-parameter-contract` (rebased onto
   master, 4 commits: exceptions, `ampere/core/parameter.py`, the spec at
   `docs/design/contracts/parameters.md`, and `tests/core/`) — but it is
   **not content-reviewed or merged**. Treat it exactly like W1.1/W1.2. It
   is substantial (2,451-line contract, 854-line spec with 112 executable
   examples, 108 passing tests) and its own report names five open
   questions for Peter and per-downstream-spec obligations for W1.4–W1.10
   — captured below — before deciding what to merge as-is vs revise.
   Read `docs/design/contracts/parameters.md` (the spec itself, §12–14
   especially) for the full design reasoning; this section only carries
   the points that need a decision, not the whole rationale.

   **W1.3 ran without W1.1's findings** (the memo branch didn't exist yet
   when it started) — its tying design is the agent's own, not checked
   against 3ML's multi-instrument pattern; re-check when reconciling.

   **Five questions W1.3 raised for Peter specifically:** (1) it imports
   `astropy.units` at module level in `ampere/core`, which contradicts
   `architecture.md`'s literal "numpy/scipy/typing/stdlib" wording for
   core even though the *plan* requires units on parameters — the spec
   needs a wording fix, not a design change, but confirm; (2) tie labels
   are a flat global namespace (`"distance"`, not `"shared.distance"`);
   (3) a `shared_as` with nothing to merge with is currently allowed
   rather than raising; (4) `OptionalDependencyError`'s shape is pinned
   per `architecture.md` §9 and should be ratified or moved at W1.13;
   (5) the legacy `npars` alias was dropped in favour of `len(pset)` /
   `pset.free_size` — confirm nothing needs the old name.

   **Obligations it places on later specs** (its own §13): W1.5 must call
   the unconstrained-space slot `Bijection` (not `Transform`, which W1.5
   owns) and have parameterised transformations inherit its `Parameterised`
   mixin; W1.6 should treat GP hyperparameters as ordinary parameters with
   `Log` bijections and explicitly answer whether array-valued
   `HierarchicalPrior` declarations scale to ~10⁵ latent values (not
   assume they do); W1.7 should build `DatasetCollection` on
   `ParameterMapping`/`distribute()` and note that **`merge()` is not
   associative** — merge all components in one call, or extend to a
   nested mapping if `DatasetCollection`s must nest; W1.8 should hash
   `to_spec()` into provenance attrs and use `free_labels()` for ArviZ
   coordinates; W1.9 has an enumerated list (§13) of declaration forms
   needing a lowering row, plus `lnprior_unconstrained` as the reference
   semantics native backends must agree with; W1.10 already has seven
   conformance-suite rows implemented and testable straight from this
   contract; W1.4 should confirm channel names and parameter names are
   allowed to collide harmlessly as separate namespaces.

None of W1.1/W1.2/W1.3's branches were pushed to origin, so `git branch -v`
in the main checkout is the authoritative list of what exists locally.

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
