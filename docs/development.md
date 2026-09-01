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

## ⚡ Pick up here (end of the 2026-09-01 session)

**State**: Phase 1 is eight of thirteen items done — W1.1–W1.6, W1.9 and
W1.12 are all Fable-reviewed and **merged to local master** (nothing is
pushed to origin; origin/master is far behind by design). The working tree
is clean; all gates green at the tip: **498 core tests**, full suite
534 passed / 5 skipped, pyrefly 0 errors, ruff lint clean. Remaining items:
**W1.7** (Dataset/FittingProblem), **W1.11** (modality sketches) — both
fully unblocked and parallelisable — then W1.8 (needs W1.7), W1.10
(needs W1.3–W1.8), and the W1.13 freeze. W0.9 remains open in Phase 0.

**Next action**: dispatch W1.7 (Opus) and W1.11 (splittable, Sonnet/Opus
per `docs/orchestration.md`) in parallel, with self-contained prompts per
orchestration principle 4, carrying the consolidated obligations below and
the worktree-trap recovery above. Then Fable-review each result at the
depth established this session (independent gate verification, full spec
read, adversarial probes of the implementation, fix small defects on the
branch with tests) before presenting for Peter's merge decision.

**Restart prompt for a new session** (paste as the opening message):

> Read docs/development.md's "Pick up here" section and the W1.7/W1.11
> obligations checklist below it, then dispatch W1.7 and W1.11 in
> parallel per docs/orchestration.md, carrying every listed obligation in
> the prompts (including the stale-worktree recovery). Review each result
> at full Fable depth as in the 2026-09-01 session and present merge
> recommendations; merge only what needs no ruling from me, and put the
> rest on my review list.

### Consolidated obligations for the W1.7 dispatch prompt

W1.7 is the item most other contracts have been leaving notes for. Its
prompt must carry, at minimum:

- **Nested merging is an open design question Peter expressly kept open**
  (`parameters.md` §14 preamble): evaluate a nested `ParameterMapping` as
  a first-class option for `DatasetCollection`, on its merits — not
  dismissed because single-call merge is what exists. Ties crossing merge
  levels (`transformations.md` §14) and **nested result channels via
  qualified flat names** (`results_schema.md` §17 preamble — symmetrical
  question, routed here) are part of the same decision.
- `parameters.md` §13: build the joint space on
  `ParameterMapping`/`distribute()`; `merge` is not associative — one
  call, or the nested design above.
- `transformations.md` §14: a `Dataset` pairs an observed container with
  an `Instrument`; the instrument's `label` is its merge component; W1.7
  owns *when* `negotiate` and `compile_for` are called.
- `likelihoods.md` §16 (all load-bearing): call
  `Likelihood.check_alignment(predicted_template, observed)` once at
  construction — it is where unimplemented latent combinations are
  refused, so it is not optional; call `check_engine(..., observed=...)`
  so a censored sample the mask excludes is not counted against a
  gradient-free engine; `Likelihood.parameters` is one flat component for
  the single merge, and `latent_declaration(n)` (n = *retained* samples)
  joins that same merge; convert this contract's `LikelihoodError` on a
  non-positive-definite covariance into §4.5's −inf-with-recorded-reason
  failure signalling (Peter's pending W1.6 §17 Q1 ruling may add a
  `strict=False` alternative).
- Plan §4.5 verbatim: `log_prob`/`log_likelihood`/`log_prior` split,
  `prior_transform`, `simulate`, capability flags, failure signalling,
  RNG policy — `lowering.md` §9.2's `substream(seed, label)` design is
  the seed-derivation policy to adopt (its home in core is W1.9 §12.7's
  ratification item).
- Accept: a toy two-dataset joint problem with a tied parameter,
  end-to-end against a stub model.

### For the W1.11 dispatch prompt(s)

Sketches (a)–(f) per WORK_ITEMS.md, plus the specific claims the merged
specs ask it to check: the deliberately-awkward instrument stress test of
the §4.3/§4.4 split (`prior_art.md` 3M1/Tension 3); whether X-ray RMF/ARF
works as a matrix multiply on an energy-axis `Spectrum` *without a
per-sample exposure concept* (`likelihoods.md` §16); what Fourier sampling
can publish as requirements given no pull-back through chains
(`transformations.md` §13.1/§14); whether closure phases need
`VonMisesFamily` before the freeze, and input to Rice's parameterisation
(`likelihoods.md` §17 Q3/Q4).

**Sol review batched at the W1.13 freeze** (Peter's call): one adversarial
pass over the §4 contract code and the lowering spec. `codex exec -m
gpt-5.6-sol` is currently refused on this account ("not supported when
using Codex with a ChatGPT account"); Peter expects to configure the extra
codex steps himself before the freeze. Nothing blocks on it until then.

## Phase 1 review outcomes (2026-09-01, all merged)

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

- **W1.5** (`w1.5-transformation-contract`, Opus): merged 2026-09-01 with
  **no defects needing Fable fixes** — the agent's own self-review caught
  the deep one (a requirements union imposing one instrument's fine
  sampling on another's broad coverage; fixed with per-interval
  densities). Mask propagation is *enforced* (a step dropping its input's
  mask raises); the `Model` ABC landed here minimally with the
  engine-facing surface reserved for W1.7. Review probes confirmed loud
  partial-value failures, correct cross-unit unions, and negotiated grids
  passing `Spectrum`'s ordering validation.
- **W1.6** (`w1.6-likelihood-contract`, Opus): merged 2026-09-01, the
  strongest deliverable of the phase; no Fable fixes needed. The critical
  defect (Student-t/Cauchy/complex-Gaussian silently discarding the GP
  while every check reported a working latent problem) was caught by the
  adversarial review the agent commissioned itself and fixed with the
  `CONSUMES_LATENT_GP` opt-in (default `False`, so third-party families
  inherit the refusal). DenseGP anchor validated to ~1e-15; Fable
  verification covered the marginal-likelihood/conditioning/whitening
  mathematics, an independent Tobit-censoring probe (exact), and the
  masking-beats-censoring declaration behaviour. Key positions on record:
  masks are consumed via `weights()` + row/column excision (the
  infinite-variance limit diverges — spec §8 has the proof); the latent
  declaration is the whitened non-centred form, *not* `HierarchicalPrior`
  (wrong at any N — GP priors are not i.i.d.); the plan's "per-sample
  log_likelihood" ambiguity (per-draw vs per-observation) is flagged for
  W1.13, since a GP likelihood has no per-observation decomposition.

**Rulings recorded this session** (all implemented, not just noted):

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
**Peter's remaining reading list** (nothing blocks dispatches; all are
freeze-relevant):

- W1.5 spec §15's open questions — most notably: whether the ANY
  mask-propagation rule is too conservative for real resampling
  (`min_valid_fraction` is the named extension); whether a model should
  be able to *refuse* a requirement rather than silently ignoring it;
  confirming the `Model` ABC's placement in W1.5 rather than W1.7 (§9 has
  the argument); and §15.8's `max_step`/`min_resolving_power` dual
  meaning (declaration vs post-union summary — the review leans towards
  demoting the scalars before the freeze, cheap now). Review at
  `docs/design/contracts/transformations.md` §15.
- W1.6 spec §17's open questions — sharpest first: whether a
  non-positive-definite covariance should raise (current) or return
  `−inf` via a `strict=False` mode (Q1, affects every engine driver);
  whether the mask union lives in `Likelihood` or W1.7's `Dataset` (Q2);
  Rice's parameterisation and von Mises's concentration (Q3/Q4 — W1.11's
  interferometry sketch will inform both); per-dataset vs per-channel
  `IndependentNoise.scale` (Q5); the conservatively-`LATENT` correlated
  complex Gaussian (Q6, a Phase-4 decision that unblocks real
  functionality); `Likelihood.to_spec()` for provenance (Q8, with W1.8).
- W1.4 spec §17 questions 3–6 (`Cube` axis order vs FITS — cheap to
  change only until the freeze; axis naming; `extra_coords` units;
  container serialisation, owed to W1.8).
- W1.9/`lowering.md` §12's remaining ratification items, all routed to
  W1.13 (discrete-family default bijection; `LoweringError` placement;
  reference-path `icdf` fallback; `substream` in core; backend-specific
  lowering registration — item 8, from Peter's review).

**Housekeeping for the next session:**

- **Nothing is pushed to origin.** All session work is on local `master`;
  pushing is Peter's call and Peter's action (agents never push).
- **Leftover worktrees and scratch branches**: several agent worktrees
  remain under `.claude/worktrees/` and the merged item branches
  (`w1.1-*` … `w1.6-*`, `w1.9-*`, `w1.12-*`) plus stale
  `worktree-agent-*` scratch branches still exist. All item branches are
  fully merged; cleanup (`git worktree remove` + branch deletion) awaits
  Peter's explicit approval per ground rule 3. Until then, note that
  `pixi run format-check` in the main checkout reports ~20 files to
  reformat — **all inside `.claude/worktrees/`** (nested copies of
  `tests/characterisation` that the task's top-level exclude does not
  match). Tracked files are format-clean; do not "fix" those.
- **Tooling trap (confirmed twice)**: `ruff format <explicit path>`
  bypasses `extend-exclude` and will rewrite the frozen W0.7 harvest
  snapshots under `docs/`. Format only `ampere/core tests/core` (or use
  the pixi task). Worth a guard in a small W0.x follow-up.
- A `SendFeedback` draft about the `/model` switch-back failure from the
  pre-handoff session remains queued; Peter can review/send it with
  `/feedback`.

## Current state (end of 2026-09-01)

- **Phase 0: complete** (W0.1–W0.8 merged; CI green; issues #74–77
  closed; archival executed). **W0.9 not started** — the temporary
  `pyphot<2`/`sbi<0.28` pins stand, documented in `DEVELOPMENT_PLAN.md`
  §2; the pyphot half must land before Phase 2's synthetic-photometry
  Transformation is written.
- **Phase 1: eight of thirteen items merged** (W1.1–W1.6, W1.9, W1.12);
  W1.7 + W1.11 ready to dispatch in parallel; then W1.8, W1.10, W1.13.
  See "Pick up here" above for the restart procedure and consolidated
  dispatch obligations.
- **Decisions ruled this session** (all in `DEVELOPMENT_PLAN.md` §2's
  table): curated astropy→native translation is opt-in/never silent; jax
  non-trainables lower via `eqx.partition`, not paramax; jax x64 is
  guard-and-raise, never set-on-import.
- The AMPERE paper revision proceeds on the legacy code and takes priority
  in any conflict over `examples/examples_paper/`.
