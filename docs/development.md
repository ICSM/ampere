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
    (since W3.2 the `dev` set plus the `torch` feature and the `sbi`
    extra — sbi 0.27, torch, pyro — and a real gate of its own:
    `pixi run -e sbi test-all`, plus `pixi run -e sbi typecheck`); `torch` /
    `jax` (the backend environments — real gates since W2.4/W2.5:
    `pixi run -e torch test-all`, `pixi run -e jax test-all`, and each
    backend package is typechecked with its real types only in its own
    environment). Since W2.11 both `torch` and `sbi` resolve torch from
    PyTorch's **CPU** index, so neither pulls ~2 GB of CUDA runtime onto a
    CPU-only runner; this changes nothing about what `pip install
    "ampere[torch]"` gives a user. **Run five-suite gates one at a time** —
    two concurrently exhaust a 13 GB machine. Torch's five-suite gate takes
    12–17 min here, sbi's 14–20 min, jax's 8–13 min, dev's 5–8 min.
  - **CI's `suites` matrix** (`.github/workflows/ci.yml`): one blocking
    job per new-namespace environment, each producing a check named
    "new-namespace suites (`<environment>`)" — `dev` (its own `suites` job:
    `test-all` + `bench`) and, in the `backend-suites` job's matrix,
    `torch`, `jax` and, since W3.7, `sbi` (each: `test-all`; `torch`/`jax`
    also run `bench`, deliberately not repeated for `sbi` since it would
    only re-exercise torch's own benchmarks under another name — see the
    job's comment). `sbi`'s install is cached the same way as `torch`'s
    (`setup-pixi`'s `cache: true`). The weekly/on-demand
    `sbi-characterisation` job is unrelated: it exercises the *legacy*
    `ampere.infer.sbi` flow, not the new-namespace `sbi` gate, and stays
    non-blocking.
  - **W4.10 — `typecheck` split out of the backend legs.** Since W4.10,
    `pixi run -e <env> typecheck` for `torch`/`jax`/`sbi` runs in its own
    matrix job, `backend-typecheck` (check name
    "typecheck (pyrefly, `<environment>`)"), independent of
    `backend-suites` (no `needs:` between them) — the two run concurrently
    rather than one serialising in front of the other, so the slow `sbi`
    leg's critical path in `backend-suites` is `test-all` alone. This was
    W3.7's finding: the `sbi` leg's local reference time (14–20 min) sits
    above the `torch` leg's (12–17 min), and splitting `typecheck` out was
    one of the two cheap levers considered there without widening that
    item's scope to touch test budgets. `dev`'s own `typecheck` job is
    unchanged (it already ran on its own).
  - **W4.10 — path gating.** A `changes` job runs first and computes, from
    a pull request's changed files, which of five buckets are touched —
    `run_dev`, `run_torch`, `run_jax`, `run_sbi`, `run_docs` — using
    `.github/scripts/path_filters.py` as the single executable source of
    the rule table (reproduced as a comment at the top of `ci.yml` and
    below). Every job in the workflow always runs; a bucket controls
    whether that job's *steps* do real work or print a one-line skip note,
    so an unaffected leg's check still reports **success**, never
    "skipped" or "expected" — the latter is what a top-level
    `on.push.paths`/`on.pull_request.paths` filter would produce (the
    workflow never triggers at all on an unaffected PR, so a required
    check named in branch protection never gets a status and the merge
    button hangs waiting for it). Push to master, `workflow_dispatch` and
    `schedule` runs are never gated — path gating only narrows a
    *pull_request* diff's turnaround time, never what lands on master or
    what a human explicitly asked to run.

    The path-filter table (buckets are exclusive; a path is classified by
    the first rule that matches — see the script's module docstring for
    the authoritative version):

    | Path pattern | Runs |
    | --- | --- |
    | `pyproject.toml`, `pixi.lock`, `.github/**`, `ampere/core/**`, `ampere/results/**`, `ampere/inference/**`, and the shared test suites (`tests/{core,results,inference,conformance,m2,benchmarks,scaling,gpu,characterisation}/**` | everything: dev, torch, jax, sbi, docs |
    | `ampere/backends/torch/**`, `tests/backends/*torch*` | torch, sbi (sbi's environment installs torch too) |
    | `ampere/backends/jax/**`, `tests/backends/*jax*` | jax |
    | `ampere/backends/reference/**`, remaining `tests/backends/**` | dev only |
    | `examples/sbi/**` | dev, sbi |
    | `examples/interferometry/**` (W4.4, not landed at W4.10 — a marked slot) | dev, torch, jax, sbi (conservative; W4.4 should narrow this to the backend(s) each example file actually exercises) |
    | `examples/**` (remainder), `tests/examples/**` | dev only |
    | `docs/**`, any `*.md` | dev, docs |
    | anything unrecognised | everything (same as the first row) |

    Locally: `python .github/scripts/path_filters.py --paths <files...>`
    (or `--stdin` with a newline-separated list) prints the bucket for each
    path and the resulting job list — the same thing the `changes` job
    prints in its own log. This is also how the item's three synthetic
    diffs were evidenced (W4.10's report).
  - **W4.10 — `actionlint`.** `pixi run -e dev actionlint` (a conda-forge
    binary, alongside `pandoc` in the `dev` feature) validates
    `.github/workflows/*.yml`; CI's `actionlint` job runs the same task,
    ungated (cheap, and useful on every PR — not only ones touching a
    workflow file). Every `uses:` in `ci.yml` is pinned to a full commit
    SHA with a `# vX.Y` comment (pinact style), and
    `.github/dependabot.yml` enables the `github-actions` ecosystem
    (weekly) so a version bump arrives as a reviewable PR updating both the
    SHA and its comment together.
  - **Required-check names for branch protection** (unchanged names keep
    their existing branch-protection entries; new ones need adding):
    `lint + format-check`, `actionlint`,
    `typecheck (pyrefly, new namespaces)`, `test (py311)`, `test (py312)`,
    `test (py313)`, `new-namespace suites (dev)`,
    `typecheck (pyrefly, torch)`, `typecheck (pyrefly, jax)`,
    `typecheck (pyrefly, sbi)`, `new-namespace suites (torch)`,
    `new-namespace suites (jax)`, `new-namespace suites (sbi)`,
    `docs build`, `minimal install (no extras)`. (`path filters` — the
    `changes` job itself — and the weekly/on-demand
    `characterisation suite (with sbi extra)` are not required checks: the
    former is plumbing every other job depends on, not a signal in its own
    right, and the latter is deliberately non-blocking, as before W4.10.)
- Plain-pip alternative: `pip install -e ".[dev]"`, Python ≥ 3.11.
- Conda env `ampere` (Python 3.13) has an editable install pointing at the
  main checkout — worktree-based work that needs importing its own changes
  should use pixi or `pip install -e .` into a fresh env instead.
- `pytest tests/characterisation` (or `pixi run test-characterisation`) is
  the "legacy still works" gate; run it before merging anything that
  touches shared files.

## ⚡ Pick up here — PHASE 5 IN FLIGHT; **the wave-2 gate is green in all four environments and recorded (`d71be43`)**; **W5.5 merged by Peter at `614d7c9` and recorded (status row, three decision-log rows)**; **W5.13 fixed, re-reviewed and ready at `17aa2c4` — its merge is waiting for Peter** (rows below); **W5.6 (Opus, `w5.6-solver-bakeoff`, base `83ddb41`) dispatched 2026-09-16 into W5.5's worktree**; system crash of 2026-09-16 recovered

**State at the end of the recovery session (2026-09-16).** (1) The wave-2 gate (W5.20 + W5.4 + W5.3): dev 2594/407 at `12c6001`, sbi 3625/92, torch 3495/222, jax 3315/231 on `b7a14b2`/`db12df5` — all green, recorded in the three status rows at `d71be43`. (2) **W5.5 is merged** (`614d7c9`, merged by Peter after Fable's second-pass review found no fixes needed — the session's auto-mode permission classifier refuses `git merge` outright, so **Peter runs every merge from now on** with the orchestrator preparing the rows; status row and three decision-log rows written; two rulings for Peter sit in the rows — the two-edit `likelihood.py` change lifting the `Layout.GRID` refusals, and `--no-dense`). Its worktree `agent-a955f0a01fed1152b` (all four envs) is reused for W5.6. (3) **W5.13 is complete on `w5.13-population-reweighting` at `17aa2c4`** (five commits, base `12c6001`, merges clean — verified with `git merge-tree`). Fable's review of the first report found one defect — `_prepare_objects` divided the stored *joint* `log_prior` column out of every draw where the identity needs the named parameter's own *marginal* interim prior (the nuisance priors cancel), biasing any object with more than one parameter, invisible on the test's one-parameter model — and the agent fixed it in `17aa2c4`: the tree carries no `PriorSpec` (provenance records names and hashes only), so `fit_population` takes a **required `interim_prior: Prior | PriorSpec`** evaluated with `ampere.core.parameter.log_density` into `_PreparedObject.marginal_log_prior`, the proposal weights keep the joint column (correct there), a `TwoParameterModel` regression test (nuisance `phi ~ N(0, 0.05)` outside the likelihood) checks the two-parameter fit against a one-parameter control and the truth, and the docstring and `results.md` item 16 carry the derivation. Fable re-read the fix (`_log_object_evidence` divides by `marginal_log_prior`; `_self_normalised_log_weights` untouched): **accepted**. **Peter merges**: from the main checkout, `git merge --no-ff w5.13-population-reweighting`; then paste the status row and the decision-log row below, and remove the worktree (`git worktree remove --force --force .claude/worktrees/agent-a9660ead86a9c2bd3`). **Then launch the wave-3 gate** (all four environments, one detached `flock` chain on master) — note the new `tests/results/test_population.py` adds ≈10 min to the dev leg and ≈20 min to sbi's (measured under contention: 852 s dev, 1427 s sbi; ≈600 s dev standalone) — a `population_full` marker on the 200-object row, as `image_full`, is the obvious follow-up (**for Peter**: the Accept line says 200 fits, so the orchestrator did not cut it).

**Ready to paste at the W5.13 merge.** *Merge message*: "Merge W5.13: population inference by reweighting archived fits — ampere.results.population with the RunColumns protocol, the DataTree and netCDF-directory readers, GaussianPopulationModel, fit_population with the ESS and spec-hash refusals and the required interim_prior, results.md §13 item 16 (Sonnet-authored, Fable-reviewed; one defect found and fixed at review)". *Status row*:

> `| W5.13 | merged 2026-09-16 at <hash> (Sonnet-authored, Fable-reviewed; merged by Peter; **one defect found at review by reading the reweighting identity and fixed by the agent in 17aa2c4**: the ratio divided by the stored joint log_prior column where only the named parameter's marginal interim prior belongs — the nuisance priors cancel — so any object with more than one parameter was biased, invisible on the one-parameter test model; the fix takes a required interim_prior argument because a run's provenance stores no PriorSpec; the original agent was killed by the 2026-09-16 system crash after two commits with its test file uncommitted, preserved by the orchestrator as a WIP commit and finished by a fresh Sonnet agent; the agent's runs: tests/results/test_population.py dev 11 passed (852 s under contention; the VI row deselected), sbi 12 passed (1427 s) including the VI row; tests/results dev 409/4 before the fix; lint/format clean, pyrefly 0 errors in dev and sbi, import sweep 83/4; the wave-3 merged-master gate to be recorded with W5.5). Landed: ampere.results.population — RunColumns (a typing.Protocol: parameter_draws, log_prior, log_likelihood, proposal_log_density, attrs), DataTreeRunColumns and NetCDFRunColumns/runs_from_netcdf_directory, PopulationModel with GaussianPopulationModel (hyperparameters as ordinary Parameters so the hyperprior needs no bespoke method), fit_population (self-normalised importance reweighting after Hogg, Myers & Bovy 2010 — uniform weights for an exact sampler, exp(log_prior + log_likelihood − proposal_log_density) for an approximate engine, refused by name without W5.0's column; emcee over α; refusals when any object's ESS at the posterior mean falls below ess_floor, and when the runs do not share ampere_spec_hash); results.md §13 item 16 (the item text's "§13.15" is stale — W5.0 took 15). Measured (200 objects, truth μ = 1.0, τ = 0.3): μ 95 % (0.943, 1.049), τ (0.285, 0.372); the two-parameter regression at N = 50: μ (0.887, 1.118), τ (0.270, 0.466), agreeing with the one-parameter control to 0.15/0.1. The W5.12 half of the Accept line is not testable until W5.12 lands. Carried: the 200-object row costs ≈10 min in dev and ≈20 in sbi inside test-all (a population_full marker wants a ruling); the required interim_prior cannot be checked against the run because provenance stores no PriorSpec (storing each parameter's spec — a schema change — would let the reader verify it, proposed); the test's class-scoped fixture is an instance method (PytestRemovedIn10Warning); ruff format must never be pointed at a .md file (it rewrites doctest fences); results.md §14's Phase 5 bullet still reads prospectively; WORK_ITEMS.md's W5.13 text cites §13.15. |`

*Decision-log row*:

> `| Population inference by reweighting archived fits (W5.13) | ampere.results.population closes design horizon (b)'s "buildable entirely on stored runs" promise: RunColumns (a typing.Protocol) plus two implementations — DataTreeRunColumns (in memory) and NetCDFRunColumns/runs_from_netcdf_directory (a directory of archived .nc files) — read the per-draw log_prior/log_likelihood every run stores and, for approximate engines, W5.0's proposal_log_density; fit_population self-normalised-importance-reweights them (Hogg, Myers & Bovy 2010) under a declared PopulationModel (GaussianPopulationModel shipped) into a population-hyperparameter posterior sampled with emcee, refusing by name when any object's effective sample size at the posterior mean collapses below a stated floor or when the runs do not share ampere_spec_hash (merged 2026-09-16, <hash>; Sonnet-authored, Fable-reviewed). The identity divides each draw by the named parameter's marginal interim prior — the nuisance priors cancel — and that prior is a required interim_prior argument, because a run's provenance records parameter names and a hash of the declaration, never the PriorSpec (found and fixed at review; the joint log_prior column is used only where it is right, turning an approximate engine's proposal draws into posterior draws). Not a PROVENANCE_SCHEMA_VERSION change: the population posterior is a new derived DataTree, not a change to what a single-object run stores; results.md §13 item 16 carries the account, including what is deliberately left out (per-object nuisance re-sampling, multi-level populations, a selection function, and the joint-fit sibling route W5.12, which this item cross-checks against rather than merges with). Conformance suite untouched. Proposed for a later row: store each parameter's PriorSpec in provenance so the reader can verify the supplied interim prior. |`

(4) **W5.6 is in flight** (Opus, `w5.6-solver-bakeoff`, base `83ddb41`, dispatched 2026-09-16 into the reused worktree `agent-a955f0a01fed1152b`; owns new `ampere/core/efgp.py`/`vecchia.py`, `tests/benchmarks/` bake-off, `likelihoods.md` §7's table; told to use `examples/image/grid_gp.py`'s out-of-tree pattern and not to edit the layout refusals; gate dev only). **When W5.13 returns fixed and is merged**: launch the wave-3 gate (all four environments — W5.5 touched core and both backends) as one detached chain on master; then the second slot goes to **W5.7 ∥ W5.9** per the ordering, W5.6 being the first. Merges are Peter's (`git merge --no-ff <branch>` from the main checkout); the orchestrator prepares the rows in the handoff before asking. (5) The docs-warning note in the operations paragraph is stale: the `Binary`/`Product`/`Sum` warning pairs do not reproduce on `9bdba46` (W5.5's agent rebuilt with `-E` and saw 13 warnings, all legacy).

*(The W5.5 rows drafted here were pasted at the merge; see WORK_ITEMS.md's table and the plan's §2.)*

**Crash recovery (2026-09-16, session after the reboot).** The machine rebooted mid-session: `/tmp` (the gate log, the `flock` file) was lost, no gate or agent process survived, but every branch and worktree did. State found and actions taken: (a) master unchanged at `b7a14b2` (only handoff commits since the dev-green `12c6001`), so the dev leg stands and the **sbi, torch and jax legs were relaunched detached** (`nohup setsid`, script `wave2-gate-rest.sh` in this session's scratchpad, one line per leg appended to `wave2-gate.log` there, full output in `wave2-gate.log.<env>`); copy the three result lines into the W5.20, W5.4 and W5.3 status rows when green. (b) **W5.5** had five commits on `w5.5-image-psf` (HEAD `dac71b7`: the reference step and grid mask propagation, the torch/jax twins, the conformance rows and `image` capability, `examples/image/` + `docs/source/image.rst` + the contract amendments, the three `Layout.GRID` fixes in both `problem.py`s) and a clean tree; its final verification and report were lost, so a fresh Opus agent was dispatched into the same worktree (`.claude/worktrees/agent-a955f0a01fed1152b`, all four envs installed) to run the acceptance criteria and quick checks, fix what fails, and report — including the decision-log rows for the found blocker (core refuses correlated noise on `Layout.GRID`; `examples/image/grid_gp.py` is the out-of-tree stand-in for a two-edit library change to be ruled on). (c) **W5.13** had two commits plus an uncommitted `tests/results/test_population.py`, which the orchestrator committed as WIP `740983d` (the agent's to squash); a fresh Sonnet agent was dispatched into the same worktree (`agent-a9660ead86a9c2bd3`, dev + sbi envs) to finish the tests, amend `results.md` §13.15, run the quick checks and report. Both agents were told to prefix every command with a `cd` into their worktree and never to touch the main checkout while the gate runs. **If this session is cut off again**: check `git log master..<branch>` on both branches, re-dispatch the same way, and re-run any gate leg without a result line.

**Where the project is (2026-09-13).** Phases 0–4 complete. Phase 4 landed interferometry end to end (kinds, steps, native twins, the circular complex GP, the study), astrometry by the template with no core change, the kernel algebra and term registry, the astropy adapter and its native translations, the SED-composition and interferometry pages, path-gated CI. The plan's §5 Phase 4 paragraph is the landed summary; every item's status row records what landed, what was accepted and what is carried; §2 has the Phase 4 decision-log rows. **Gates on master**: the first wave (dev 2312/284, sbi 3145/112, torch 3023/234, jax 2847/238), the second (2380/308, 3250/114, 3128/236, 2951/241) and W4.3's (2392/371, 3349/90, 3222/217, 3048/219) all green; **the final wave's (W4.4 + W4.9, master at `318034e`): dev 2450/379, sbi 3427/90, torch 3300/217, jax 3125/220 — all green; W4.8's dev gate (master at `c8bf2c0`) green at 2450/379** — every Phase 4 wave gated green in all four environments. Nothing is in flight; no agent worktrees exist; nothing is pushed (origin is still at `b8e585b`). The `w4.*` branches and the `worktree-agent-*` harness branches are all merged and can be deleted (part C).

**What a clean session does next.** (1) **Wave 1 is merged and gated green** (counts in the header and in W5.2's status row). (2) **W5.20 is merged** (`acaad79`; decision-log row and status row written; its merged-master gate is taken with wave 2's, after W5.4 merges). **The weekly usage limit killed both Opus agents mid-item on 2026-09-16; W5.20 was resumed and finished; W5.4 is still in flight after its resume.** W5.4 on `w5.4-hsgp-solver`: two real commits (the Hilbert-space basis and the reference `HilbertSpaceGP`; the torch and jax solvers) plus a WIP commit of its in-progress conformance rows (`tests/conformance/protocol.py` and the four backend fixtures) — remaining: the conformance rows, the latent-path agreement row, NUTS on both backends, `likelihoods.md` §7 and `inference.md` 17.4, the report. W5.20 on `w5.20-parameter-namespace`: one real commit (the reserved set) plus a WIP commit of its in-progress renames across both backends' model modules, `problem.py`s, `torch/parameters.py`, tests, and a new `tests/backends/test_parameter_namespace.py` — remaining: finish the renames, the torch lowering guard's test, the ambiguity refusal, the docs, the report. **Resume each from its branch** (a fresh agent with the original prompt from this session's dispatch plus "start from the branch's HEAD; the WIP commit is yours to amend", or the same agent by message if the session still holds it); both agents' worktrees are under `.claude/worktrees/` with their own `.pixi` environments already installed. (2b) **W5.4 is merged** (`359386c`; decision-log and status rows written; two of its decisions — hashing `m`/`c`, `SpectralMixture` supported — are flagged for Peter in the row). (3) **W5.3 is in flight** (Sonnet, `w5.3-multi-axis-plots`, base `5af209f`, dispatched 2026-09-16; owns `ampere/results/{_plotting,plots,derived}.py`, the `PLOT_COORDINATE` attribute on the kinds, `tests/results`, the study's figures, `results.md` §4/§13). **W5.3 is merged** (`47988b4`; rows written). **The wave-2 merged-master gate** (dev, sbi, torch, jax, one at a time through the `flock`) was launched at the merge and covers W5.20, W5.4 and W5.3 together — **dev green at 2594/407 (master `12c6001`)**; the harness's memory guard then killed the gate's shell chain during the sbi leg, leaving an orphaned unlocked `pytest` that was killed by hand, and the sbi, torch and jax legs were relaunched **detached** (`nohup setsid`, script `wave2-gate-rest.sh` in the orchestrating session's scratchpad, results appended to `wave2-gate.log`) so the guard cannot take the chain down again; copy its four result lines into those three status rows when green (if this session was cut off: re-run `flock /tmp/ampere-gate.lock pixi run -e <env> test-all` per environment on master and record). **W5.5 is already in flight** (Opus, `w5.5-image-psf`, base `9bdba46`, dispatched into W5.4's freed slot; owns the reference/torch/jax PSF step, the grid-mask helper, `tests/conformance/test_image.py` and the `image` capability in `protocol.py`, `examples/image/`, `docs/source/image.rst`, `transformations.md` §10/§13.5). **W5.13 is in flight** (Sonnet, `w5.13-population-reweighting`, base `12c6001`; owns `ampere/results/population.py`, `tests/results/test_population.py`, `results.md`'s population item). After W5.5: W5.6 (the bake-off) and W5.7 ∥ W5.9 per the ordering. Note for targeted runs: `tests/backends/test_torch_device.py::TestTheComplexGaussianBody` fails outside `test-all` (`No module named 'tests'`, pre-existing) — not a regression. Then W5.5 ∥ W5.13 per WORK_ITEMS.md's ordering. **Lesson recorded in `docs/orchestration.md`'s token-economy section, as corrected by Peter**: the week was nearly spent before the session began and the remainder was used up deliberately before the reset (the session itself drew about 8 % of the weekly allowance); Peter budgets the week and the session start, the orchestrator the five-hour window; the incident's real lesson is that interruptions must be cheap — commit often, preserve WIP on a kill, keep the handoff resumable — and they were, since nothing was lost. The Phase 5 draft's inputs, for whoever re-reads it:

**Owed / carried from Phase 4** (beside the Phase 3 list further down): `NUTSEngine` tuning knobs (`target_accept`, `max_tree_depth`; W4.11); `examples/` excluded from ruff; `dataset.py`'s `part_name` docstring; a `Sum` containing a `Product` gets the generic O(N) refusal; W4.7's `BATCHABLE` unverified under `vmap`; four of the six plots refuse a point kind with several axes (Phase 5); `results/derived.py` casts complex to real with a warning in two places where a third refuses; `RiceFamily` declared-only with no schedule; `von_mises` not in `_TWINNED_FAMILIES` (native draws fall back to numpy); `encoding._sanitised`'s summary line over-promises; `QuasisepGP.conditional_loo` unimplemented on numpy; jax `bessel_jn` unusable across a realistic range and torch `bessel_j1` has no backward (documented on the uniform-disc models); two pre-existing docs-build warning pairs (`Binary`'s field list; duplicate `Product`/`Sum` descriptions); `interferometry.rst` attributes a non-verbatim quotation to the plan; the terra reviews owed on W4.1, W4.5 and W4.2 when the quota returns (~2026-09-30).

**The Phase 4 orchestration record** (for whoever reads only this): four agents dispatched in parallel on 2026-09-11 were killed twice by the session usage limit while polling gates with cold caches; the token-economy rules in `docs/orchestration.md` (agents never wait for gates; the orchestrator gates each merged wave once; bounded polling under the cache TTL) were adopted mid-phase and the second and later waves ran under them — every subsequent agent finished in 17–100 minutes. Two harness traps were re-hit and are in memory: a `cd` into a worktree persisting into a `git merge`, and a single Bash call exceeding ten minutes being pushed to the background.

**Where the project is.** Phases 0–3 complete: the spec frozen at `spec-v1.0` (2026-09-03), three backends in lockstep, M2 reached, and **the SBI layer closed 2026-09-10** — every Phase 3 item W3.0–W3.15 merged, W3.13 (the documentation pass) last at `0ff223e`. Each item's row in WORK_ITEMS.md's status table is the authoritative record of what landed, what was accepted at review and what is carried; the plan's §5 Phase 3 section carries the landed summary and §2 the decision-log rows. **Gates on master, all green, all four after the last code merge (W3.16, `fea7ef9`)**: dev 1996/268, sbi 2784/80, jax 2491/201, torch 2667/197. Lint/format/pyrefly clean in every environment. Lint/format/pyrefly clean in every environment. **Nothing is in flight, no agent worktrees exist, and nothing is pushed to origin** (origin/master is still at the W1.3 handoff, `b8e585b`). Master is clean apart from Peter's untracked paper script under `examples/examples_paper/` (leave it).

**Environments.** `dev`, `torch`, `jax`, and **`sbi` (= dev + torch feature + sbi extra since W3.2; run `pixi install -e sbi` after pulling)**. The gate is `pixi run -e <env> test-all` in each; an item's gate set is dev + sbi, plus torch and jax when it touches `ampere/core` or a backend. **Gates are serialised through `flock /tmp/ampere-gate.lock`** — every agent and the orchestrator run `flock /tmp/ampere-gate.lock pixi run -e <env> test-all`, so two five-suite runs never overlap on this 13 GB machine (the harness kills background tasks on low memory). Torch's takes 12–17 min, sbi's 14–20, jax's 8–13, dev's 5–8.

**What a clean session does next, in order.** (1) Read the "Everything that needs Peter's attention" block below first; any ruling Peter has given since becomes either a direct edit (the one-liners in part A: the `sample_with` default, the zero-rate guard, `cauchy`) or a small item — do those before anything new. (2) **The Phase 4 items are drafted** (WORK_ITEMS.md "Phase 4", `5488e6d`: W4.0–W4.9 and four decisions D1–D4 — the `ampere/modalities/` namespace, the `ClosurePhases` signature, whether W4.9 runs, the model split) **and await Peter's approval before any dispatch**; once approved, dispatch W4.0 ∥ W4.1 ∥ W4.5 ∥ W4.6 from the standing paragraph. The inputs they were drafted from: the inputs are the plan's §5 Phase 4 bullets (interferometric visibilities end to end — decided; the astropy interop adapter, §4.7; kernel algebra and the public quasiseparable-term registry, ruled into Phase 4 on 2026-09-09) and `docs/design/modalities/interferometry.md`, whose §9 (interface gaps), §10 (requirements on the joint fit) and §11 (open questions) are the item list in all but name; `docs/design/modalities/astrometric_timeseries.md` and `ifu_cube.md` are the follow-on templates. Draft in W2/W3 style — one branch per item, an Accept line each, a conformance row per new container or transformation, the gate set named — and size them S/M for Sonnet/Opus per `docs/orchestration.md`. (3) If Peter approves the memo's §5 adaptations or §9 horizons (part B), draft those as Phase 5 items in the same pass. (4) The owed small items below are Sonnet-sized fillers while gates run.

**Phase 3 in one paragraph, for whoever reads only this.** `ampere.inference.SBIEngine` over sbi 0.27 (NPE/NLE/NRE, and TMNRE through sbi's own ratio estimators — no swyft, ruled) fits any `FittingProblem` through `simulate_many`'s executor protocol (`forkserver` default, native batched prediction and sampling on torch and jax); the coordinate–value–mask encoding and the masked set/transformer embeddings; trained-artefact caching keyed on spec, model and data hashes plus the run's settings; SBC/TARP calibration usable by every engine; non-native parts refused by default and opt-in without gradients; core `sample()` for Poisson, Student-t and complex Gaussian with native twins; `PROVENANCE_SCHEMA_VERSION` 6 with `ampere_model_hash` and the append check; runs reproducible bitwise from the problem's seed; the blocking sbi CI leg; plot paging; the SBI tutorial. Carried into the deferred list: the embedding study with the ε study folded in; jax-native SBI on maturity.

**Standing dispatch paragraph** (put in every prompt): base commit and the ancestor/descendant recovery recipe from `docs/development.md` "Dispatching work"; the `flock` gate rule and "never stop your turn to wait — poll with a foreground `until`/`sleep` loop"; explicit file ownership when two agents run in parallel; do not edit the WORK_ITEMS.md status table or DEVELOPMENT_PLAN.md (decision-log text goes in the report); the commit trailer; the report format from the `work-item` skill. Models: Opus for judgement-within-spec, Sonnet for well-specified items (`docs/orchestration.md`).

**Decisions recorded this phase** (all in the plan's §2 or the status rows): `forkserver` default; native sampling the default; further core `sample` families on a use case with the pathway kept open (structurally enforced by W3.8) — met at W3.14; the mesh sharder's padding to be costed with the GPU item; the observation-context reservation (W3.3's encoding carries σ by default; `context=` accepted as `None` and recorded); kernel algebra into Phase 4; the embedding *study* and the ε study in the deferred list; W2.9's coverage study inside W3.6; the `marginals` group stored by default (W3.4); the Poisson/Student-t/complex-Gaussian sampling forms (W3.14); schema 6 and the model hash (W3.12); plot paging (W3.10); the seeding rule (W3.15); the columnar store for population work waits but nothing may exclude it (2026-09-10, horizon (b)).

**Rulings received 2026-09-15 (Peter), on the 2026-09-13 review block**: (1) W4.1's interferometric conventions **accepted**; (2) the four Fable-only decision-log rows **confirmed**; (3) W4.3's spelling: **(a) assumed provisionally**, and he asked for the consequences of the shadowing/duplication to be expanded — the analysis is in the plan's W4.3 row and W5.20's text (in one line: the reserved parameter names vary by backend because the check walks `type(self)`'s MRO, torch's `nn.Module` adds a third list, nothing reads a parameter as an attribute so the rule guards a convention, and the lookup order takes a class offering both pairs at the legacy name); (4) the token-economy rules: **continue as at present, re-evaluate later**; (5) W4.7's refusals and W4.10's path gating: he asked for more explanation — given in the session of 2026-09-15 (a tie is an opaque Python callable evaluated on floats, so on jax it cannot run on a tracer and on torch it would detach the gradient, and `|`/`&` have no elementwise reading as one channel's flux; results and inference are backend-neutral code that every engine on every backend writes through, so a change there skipping the torch or sbi leg is a change CI would not test) — awaiting his word; (6) W4.2's coverage pin **confirmed**. Phase 5: D2–D8 **as recommended** (D8: two Opus items, margins pinned); D1 provisional. **Closed later the same day**: D1 ruled as (a) and (b) together, W5.20 in; the W4.7/W4.10 explanations accepted ("the judgement calls make sense"). Nothing from the 2026-09-13 block remains open. The original block follows for the record.

**Awaiting Peter's review (2026-09-13; he will look "a little later")** — *ruled 2026-09-15, see above*: (1) W4.1's interferometric conventions — closure-phase sign and canonical ordering, baseline sign `r_j − r_i`, the DFT sign with x east and y north (`ClosurePhases`' docstring; `tests/conformance/oracles.py`'s header); (2) the four Phase 4 decision-log rows merged on Fable's second pass alone — W4.5 (Matérn-5/2 exactly quasiseparable, §15.3 withdrawn), W4.1's `VisibilitySet` amendment, the von Mises implementation, W4.2's `(n, k)` solver contract; (3) the token-economy proposal in `docs/orchestration.md` (three rules already in force); (4) W4.7's refusal of tied parameters and `|`/`&` on the native route, and W4.10's "results and inference count as core" for path gating; (5) W4.2's calibration row pinning coverage rather than the rank test where the GP over-covers. W4.11's far-infrared filter set was accepted the same day.

**Everything that needs Peter's attention, consolidated at the Phase 3 close (2026-09-10).** Rulings go in the status rows / decision log; nothing here blocks the close-out.

**Rulings received 2026-09-10 (Peter): parts A and B are all ruled and applied** — A1 and A6 landed as **W3.16** (merged `fea7ef9`; merged-master gates all green, recorded); A2–A4, A7, A8 confirmed in the decision-log rows; A5's principle in `likelihoods.md` §3; A9 as W4.10; A10 as W4.0 (1) and (6); B11–B14 and B17 in the memo's §10; B15's horizons (e)–(g) in the plan; B16 as W5.0. **Only part C remains open**, for Peter when he gets the chance. **Phase 3 is finished** (2026-09-10). The horizon notes are folded into the plan's §5 Phase 5 bullets and horizons (h)–(i) (`152d187`), so the whole horizon is visible in one place.

*A. Rulings on merged Phase 3 work (each is a one-line change or a "leave it")*
1. **W3.4 — `sample_with` default for TMNRE.** Current default `"rejection"` (i.i.d. draws). Measured 353.6 s against 43.8 s for `"mcmc"` on the example, and rejection's cost *rises* as truncation succeeds; `"mcmc"` is sbi's own default for ratio posteriors and what `method="nre"` uses. Recommendation: switch the default to `"mcmc"`, keep `"rejection"` selectable.
2. **W3.4 — the accuracy row.** The 0.5σ location tolerance is unreachable for a ratio estimator's joint on the toy problem; the suite asserts coverage plus the width band and holds the marginals to 1.0σ. Confirm, or state a tolerance.
3. **W3.4 — `ε = 1e-4`** gives a ±4σ box that stops shrinking after round 2. Recommendation: leave; the ε study is now in the deferred list (W3.13 entered it beside the embedding study).
4. **W3.4 — per-round pair estimators.** The 2-D marginal estimators train only on the final round; per-round pairs would cost budget for a box the 1-D marginals already set. Confirm the final-round-only choice.
5. **W3.14 — `cauchy`.** Refused only because the item listed three families; it is one branch away. Rule: a list (leave) or a principle (add it, and write the principle into `likelihoods.md` §3).
6. **W3.14 — the zero-rate asymmetry.** `PoissonFamily.sample` guards `rate < 0`, `log_prob` guards `rate <= 0`. Recommendation: align `sample` to `<= 0`.
7. **W3.12 — schema 6** (`ampere_model_hash` on every run and training set; a pre-schema-6 training set refused on append by name rather than guessed) and the `model_hash` name in place of the item's `model_identity_hash(problem)` (which already existed with other semantics). Confirm.
8. **W3.15** was drafted and dispatched on Fable's judgement: torch *and* numpy's legacy global generator seeded per step and restored. Confirm retroactively, or strike.
9. **W3.7 — the sbi CI leg** is likely the slowest (local 14–20 min against torch's 12–17). Accept, or apply a lever (typecheck split into its own job; tighter smoke budgets). Also: add `actionlint` to the dev toolchain?
10. **W3.13** — (a) `docs/source/overview.rst` is framed "as of M2" and still says five engines and "SBI is Phase 3": approve a small refresh item (W3.16). (b) No `examples/sbi` script does a bare in-process NPE fit; the tutorial covers the engine through the black-box case. Want a plain native example?

*B. The inference-extensions memo (`docs/design/inference_extensions_memo.md`)*
11. **§7.1 — the stored proposal density** as a rule for every approximate engine (VI today does not store it; SBI does). The hook with the widest consequences: importance correction, population reweighting at scale, retroactive density emulation. Yes/no.
12. **§7.2 — optimisation results**: a `DataTree` with an `optimum` group, or a typed object with `to_datatree()`.
13. **§7.3 — BO as an engine** or only as acquisition for emulation/multi-fidelity (recommendation: the latter; the *surrogate-posterior* family — VBMC, GPry — is filed separately at tier 2 as the route for simulators too expensive for SBI).
14. **§7.4 — dependency policy**: may `nautilus`/`ultranest` join the base install, or every new sampler behind an extra as `zeus` is?
15. **§9 — three new design horizons** proposed for the plan's list: (e) model comparison and averaging, (f) approximate-inference correction, (g) engines as a registry with uniform cost accounting. Approve and Fable adds them as one paragraph each; the per-horizon reservations in §9 (a fidelity column in training sets; the interim prior reconstructible from the archive; an emulator's identity including its training set; nested encodings) are sentences for the contracts they name, added by the first item that touches each.
16. **§5's six adaptations** — should they become drafted Phase 5 items now (the evidence attr, the weighted/approximate-draw rule, `ampere_approximation`, the `Optimum` result, batched `log_prob`, dependency tiers), or wait for the first engine that needs them? Recommendation: draft the first three as one small item, since they are results-contract text and the samplers of tier 1 want them.
17. **§8.6 / horizon (b) — the columnar store**: ruled "wait, but nothing may exclude it"; recorded as a constraint on horizon (b). Nothing further unless the ruling changes.

*C. Repository actions only Peter takes*
18. Push master (origin is still at the W1.3 handoff `b8e585b`) and watch the first live CI run: the branch-protection rule naming the old `phase1-suites` job needs `suites`, and **a required-check entry for `new-namespace suites (sbi)`** beside `(torch)`/`(jax)` (W3.7).
19. Delete the harness branches: `git branch --list 'worktree-agent-*' | xargs -r git branch -d`.
20. `pixi run gpu` on a machine with an accelerator, and cost the mesh sharder's padding there (W3.1 slice 2's ruling).
21. The retroactive `gpt-5.6-terra` reviews when the Codex quota returns (~2026-09-30): the Phase 2 queue (W2.3 first), then W3.1 slices 1–2, W3.2, W3.3, W3.8, W3.6's calibration statistics, W3.12's hash composition, W3.14's sampling forms.
22. `examples/examples_paper/flexible_likelihood_comparison.py` is untracked on master — his paper work; commit or ignore as he sees fit.

*D. Next phase*
23. **Phase 4 (one new modality end-to-end — interferometric visibilities, decided)** is the next to draft items for, with `docs/design/modalities/` as the input; Phase 5's inference items would follow from B. Fable drafts for approval before any dispatch.

**Owed / carried small items** (added 2026-09-13: `NUTSEngine` should expose `target_accept`/`max_tree_depth` — W4.11's modified-blackbody variant ran past 19 min at a 40/80 budget; `examples/` is excluded from ruff so example scripts escape the lint gate; `dataset.py:395`'s docstring cites `ampere.core.likelihood.Matern32`, now in `ampere.core.kernels`; a `Sum` containing a `Product` gets the generic O(N) refusal; W4.7's `BATCHABLE` unverified under `vmap`; **four of the six plots refuse a point kind with several axes** (`VisibilitySet`, `ClosurePhases`) — Phase 5's multi-axis diagnostics, W4.4's finding; `results/derived.py` casts complex to real with a warning in two places where a third refuses) (not scheduled; pick up when convenient, Sonnet-sized): a native sampler failure aborts a whole chunk where the numpy loop flags one draw (`_SimulateNatively.run`, W3.14's finding — wrap the sampler in the per-draw `try`); `Instrument._declarations()` fingerprints steps by `id()`, so a *bare* pickled `Instrument` still breaks (a value-based fingerprint in `transform.py`; `Dataset.__setstate__` works around it); torch `noise.py`'s `_check_one_device` partly redundant with the kernel-as-capability-part check; the "no uncertainty at all" branch of both backends' `_check_gp_uncertainty` is unreachable; upstream sbi 0.27 issues worth filing — `PermutationInvariantEmbedding` reads the first batch element's valid-row count for the whole batch, and `TransformerEmbedding` drops `attention_mask` unless causal; `z_score_x="none"` does not silence sbi's constant-column warning; the Phase 6 multi-observation composition tutorial (plan §5); the `tests/examples` smoke test of photometry + spectra composition.

**Operations notes for the orchestrator**: one five-suite gate at a time, always through the `flock`; run long gates detached (`nohup … &`) and wait with a background `until` loop rather than polling in the foreground; **never edit, format or check out source files in the main checkout while a gate runs there** (a half-written module crashed a subprocess test on 2026-09-09; markdown is safe); an agent that `cd`s to the main checkout can edit its files even though the isolation guard blocks git there — tell agents to stay in their worktree; after a lockfile change, `pixi install -e <env>` in the main checkout; the merged-master gate is the one that counts; parallel agents appending decision-log rows conflict at the table's last line — keep both rows, keep the table contiguous; agent worktrees under `.claude/worktrees/` may hold their own `.pixi` (2 GB) — remove a merged agent's worktree with `git worktree remove --force --force`; finished agents linger in the task list until dismissed with `TaskStop`; an agent that "stops to wait" for a gate is resumed with `SendMessage` and told to poll; a session rate limit kills agents mid-flight but their committed work survives — resume them from their branches.


---

Earlier session records (Phase 0–1 review outcomes, the 2026-09-01 to 2026-09-09 handoffs, including the incremental mid-Phase-3 record) are archived verbatim in `docs/handoff-archive.md`.
