# Ampere v2 — orchestration & model-selection policy

How to spend model capability where it matters. The orchestrating session
(Claude, holding the full project context) applies this when dispatching
work; humans dispatching agents directly should follow it too.

## The ladder

| Tier | What | When |
|---|---|---|
| Fable (orchestrator) | Frontier Claude, holds full project context | Architecture and contract decisions, plan changes, ambiguity resolution, integration, final review of agent output. Anything that changes a frozen agreement. |
| Opus subagent | Strong Claude | Implementation requiring design judgement *within* a frozen spec — Phase 2 backend tracks, contract implementations. |
| Sonnet subagent | Capable Claude, cheaper | The default for well-specified work items with tight acceptance criteria (most W0.x, harvest/memo work, test-writing from a spec). |
| Haiku subagent | Cheap Claude | Mechanical sweeps: formatting checks, bulk file edits from an exact recipe, status chores. |
| Codex + gpt-5.6-sol | Frontier OpenAI (separate budget) | Detailed and adversarial reviews — see the `sol-review` skill. Cross-provider review has uncorrelated failure modes with Claude-written code, which is the point. **Currently refused on this account** (ChatGPT-auth codex); the ruled stand-in (2026-09-03) is **`gpt-5.6-terra`**, run through the same skill's procedure with `-m gpt-5.6-terra` — proven at the W1.13 freeze, where it found two real pre-existing defects. |
| Codex + gpt-5.6-luna | Cheap OpenAI (separate budget) | Bulk mechanical work when preserving Claude quota, or independent second drafts of small pieces. |

Two mechanics worth knowing: Claude subagents **inherit the orchestrator's
model unless an override is given** — always set the model explicitly when
dispatching legwork; and a *fork* inherits both context and model — use it
only when full conversation context is genuinely required, never for cheap
work.

## Principles

1. **Match the model to judgement density, not task size.** A 14-line
   `.gitignore` item is Sonnet work; a 40-line contract definition is
   Fable work, because its blast radius is the whole project.
2. **Cheap generation, strong verification.** Anything produced by a
   cheaper model is reviewed by a stronger one before merge. Reviewer
   findings (from any model) are advisory: verify each against the code
   before acting; never auto-apply.
3. **Adversarial review before merging anything that touches the §4
   contracts, likelihood/GP mathematics, or lowering rules** — via the
   `sol-review` skill, in addition to the orchestrator's own pass.
4. **Self-contained dispatch prompts.** Embed the work-item text,
   acceptance criteria, and constraints; never rely on the agent's
   environment matching expectations (learned the hard way: agent
   worktrees have been created from `origin/master`, not the
   orchestrator's HEAD).
5. **One agent per work item; don't decompose below item granularity** —
   spawn overhead and context re-derivation dominate small tasks. The
   orchestrator does small things itself.
6. **Escalate ambiguity upward, not sideways.** An agent reporting a
   blocker gets an answer from the orchestrator (or Peter), not a
   re-prompt to a cheap model to guess.
7. **Budgets are separate pools.** Claude tiers draw on the Anthropic
   plan; Codex tiers on the OpenAI plan. Shifting mechanical bulk to luna
   preserves Claude quota for the work that needs it.

## Default phase mapping

- Phases 0–1: **complete** (spec frozen at `spec-v1.0`, 2026-09-03; the
  freeze's adversarial pass ran as Fable + `gpt-5.6-terra`, sol being
  blocked). Historical mapping: Phase 0 to Sonnet; contract specs
  Fable-drafted or Fable-reviewed; mechanical companions to Sonnet/Opus.
- Phase 2: **complete** (closed 2026-09-08). Backend tracks W2.4/W2.5 went
  to Opus, one per track, with the conformance suite as the cross-check;
  well-specified items to Sonnet or Opus by judgement density; the
  cross-model passes are owed retroactively (Codex quota) — the queue is
  in `docs/development.md`.
- Phase 3 (current, the SBI layer, W3.x): contracts and integration to
  Fable (the encoding contract was Fable-drafted before its item ran);
  items needing judgement within the spec (`simulate_many`, the engine,
  the encoding, native sampling, foreign parts, calibration, TMNRE) to
  Opus; well-specified items (housekeeping, docs pages, caching, CI, the
  embedding-readout and model-hash follow-ups) to Sonnet. Two agents run
  in parallel only with explicit file ownership in both prompts, and all
  gates serialise through `flock /tmp/ampere-gate.lock`. Cross-model
  review of the SBI mathematics (the batched-equals-loop and native
  sampling claims, the unconstrained prior, the calibration statistics)
  when the quota returns.
- Reviews of merged-candidate PRs: orchestrator (Fable) pass always;
  cross-model pass for contract/maths-touching changes.

## Token economy (proposed by Fable 2026-09-13 from Phase 4's first wave; **ruled by Peter 2026-09-15: the rules have worked well — continue as at present and re-evaluate a little later**; all seven stay in force as practised, none struck)

**What happened.** The first Phase 4 wave dispatched four agents at once
(one Sonnet, three Opus). Each coded for about an hour, then spent one to
three hours waiting for the five-suite gates, which are serialised through
one lock on a 13 GB machine (dev 8 min, jax 11, torch 15, sbi 20–25, so
sixteen gate runs ≈ 4 hours of wall clock behind a single lock). The
session hit the 5-hour usage limit twice in one day, both times more than
an hour before the window ended, each time killing all four agents
mid-flight (their committed work survived). Reported usage per agent:
W4.0 320 k tokens / 258 tool uses, W4.5 467 k / 314; the three read-only
surveys 123–152 k each on Opus.

**Where the tokens went.**

1. **Waiting is the dominant waste.** An agent polling a gate wakes with its
   whole context every poll. A subagent's prompt cache lives about five
   minutes; the polling loops slept eight, so every poll was an *uncached*
   re-read of 200–400 k tokens, ten to fifteen times per agent. That is
   several million uncached input tokens per agent spent watching a
   deterministic process that needs no model at all. Counter-intuitively,
   a shorter sleep (under five minutes) would have been ten times cheaper
   per poll, because cached reads are billed at a fraction of uncached.
2. **Parallelism that the lock defeats.** Four agents coding in parallel
   saved wall clock only during the coding hour; from then on they queued
   on one lock, each keeping a live context while idle.
3. **Context bloat.** 250–320 tool uses per agent, most of them reads of
   long modules (`likelihood.py` 3 400 lines, `dataset.py` 3 900) in
   full, several times, plus every gate log tail.
4. **Model tier for legwork.** The surveys, which only locate and
   summarise, ran on Opus 5.

**Rules proposed (the first three applied to the second wave on
2026-09-13 by message; the rest await the ruling).**

- **Agents never wait for gates.** An agent runs the quick checks (lint,
  format, typecheck, the import sweep) and the *targeted* test files its
  item touches, each command under **five** minutes (amended 2026-09-22
  from ten — the subagent cache lives five; a longer run goes detached
  with a keep-warm poll, see the re-evaluation below), then commits and
  reports.
  The five-suite gates run **once per merged wave** on master, launched by
  the orchestrator as a detached shell chain (`nohup … flock … &`) that
  needs no model; failures come back to a cheap fix-up agent with the log
  excerpt. Branch gates were a confidence measure; the merged-master gate
  is the one that counts and always was.
- **Nothing sleeps longer than the cache.** If anything must poll (the
  orchestrator on a gate it launched, say), it sleeps under the cache TTL
  — under five minutes for a subagent, under an hour for the orchestrating
  session — or it does not poll at all and lets the task notification
  wake it. The `Monitor` tool, which waits on a condition without a model
  turn, is the right primitive where available.
- **Two agents at a time, not four,** and staggered so their targeted
  test runs do not collide on the lock: dispatch the Opus item first, the
  Sonnet items an hour later. Use the last hour of a usage window for
  gates and reviews (cheap, cached, or model-free), not for fresh
  dispatches.
- **Read with a scalpel.** Prompts name file *and line range* for every
  read the item needs; agents are told to read sections, not modules, to
  use `git diff` rather than re-read a file they edited, and never to
  print more than the tail of a log. A survey that only locates code is
  Haiku work; one that summarises is Sonnet work; Opus reads only what it
  must judge.
- **Model tiers by judgement density, including the 4.6 generation.** If
  Opus 4.6 and Sonnet 4.6 draw less of the usage quota than the 5-series
  (to be checked on one item, not assumed), route: Sonnet 4.6 for
  well-specified items with tight acceptance lines (W4.7, W4.10, W4.11,
  documentation passes, CI); Opus 4.6 for judgement-within-spec where the
  contracts are settled and the tests define the answer (a native twin of
  an existing reference step, W4.3's kind); Opus 5 for GP mathematics,
  new contract surface and anything with a decision-log row (W4.2); Fable
  for integration, rulings and reviews. Measure the first such item's
  usage and record it here.
- **Make a fast gate.** `test-all` runs seven suites; most items touch
  two. Add a `test-fast` task that excludes the `m2_full` rows and the
  SBI training budgets, and a per-suite gate list in each item's Accept
  line, so the full four-environment gate is a wave-level event. Consider
  `pytest -x` on branch runs (stop at the first failure) and the
  `pytest-xdist` question only for `dev` (memory forbids it for torch and
  jax on this machine).
- **Orchestrator discipline.** Record state in the repository at each
  merge (done) but avoid turns whose only content is "still waiting":
  one bounded poll per gate set at most, otherwise wait for notifications.
  Batch the status-table, decision-log and handoff edits into one commit
  per wave.

**Second incident (2026-09-16), and Peter's correction of its reading.** The *weekly* usage limit ended a day of Phase 5 work: two Opus agents (W5.4 at about 40 % of its item, W5.20 at about 25 %) and two finished Sonnet agents were terminated at once, and the harness's low-memory guard killed the orchestrator's jax gate leg while the agents' own test runs were competing for the 13 GB. The orchestrator first read this as "Opus items burn the weekly allowance"; **that was wrong** (Peter, 2026-09-16): the session had used only about 8 % of the weekly allowance — the week was already nearly spent when the session began, and the remainder was *deliberately* used up in the five hours before the reset. The division of labour is therefore: **Peter budgets the week** and chooses when sessions start with the weekly window in view; **the orchestrator budgets the five-hour window** as the rules above already say. What the incident does teach is about *outages*, not budgets — the objective is that an interruption from either limit is short and loses nothing. Two rules follow. (8) **Make every interruption cheap**: agents commit logical units as they finish them (already the rule), the orchestrator preserves any uncommitted worktree edits as a WIP commit on the agent's branch the moment a kill is noticed, and the handoff names each branch's head and what remains so any session can resume it — this incident lost no work because all three held. (9) **Every pytest run goes through the lock while a gate runs** — dev runs included, not only torch and jax — because the memory guard does not distinguish who is at fault.

**Expected effect.** Per agent, the waiting cost (the bulk) goes to zero;
coding cost is unchanged; the gate wall clock is unchanged but no longer
holds a context open. The first wave's four agents would have finished
their reports within about ninety minutes each instead of three to four
hours, well inside one usage window.

### Re-evaluation (2026-09-22, Fable; **ruled by Peter the same day: all five adopted, the lifetime cap soft**)

**Method.** Every request in the project's 126 session transcripts since
31 August (22 419 requests; main sessions and subagents), weighted at the
API's relative rates — cache read 0.1, cache write 1.25, fresh input 1,
output 5. The script is not in the repository; the numbers are.

**Where the tokens went.** Cache reads 71 %, cache writes 22 %, output
6.5 %; subagents 89 % of the total, main sessions 11 %; Opus 61 %,
Sonnet 26 %, Fable 13 %. The bill is context size × request count, not
what is written.

**Finding 1 — the cache tier.** Every main session since 11 September ran
on the one-hour cache; **every subagent ran on the five-minute cache**
(the transcripts record the tier). The "under ten minutes" rule was
therefore wrong by construction for agents: since the 15 September
ruling, 16 agents suffered 68 cold restarts after gaps of five minutes
or more (median gap 9.8 min — a targeted test run), costing 16 % of all
agent usage; touching the cache every four minutes instead would have
cost a third of that.

**Finding 2 — context size is first-order.** Requests over 300 k of
context carried 63 % of all cost; capping every request at 200 k would
have left the total at 63 % of what it was. The heaviest agents made 500
to 850 requests at 280 to 390 k average context. Tool-result volume is
small (about 9 M tokens across every agent ever), so contexts are not
built from a few huge reads: they accumulate over hundreds of turns, and
the largest single items are still whole reads of 60–75 k-character
modules and design documents.

**Rules adopted (in force from wave 5).**

1. **Soft lifetime cap with a hand-off.** An agent past **250 k of
   context, or about 300 tool uses, or 2.5 hours**, whichever first, does
   not start a new unit of work: it finishes the unit in hand (never
   mid-edit; a WIP commit if the boundary cannot be reached cleanly),
   commits, and returns a report with a **Hand-off** section — what is
   done, what remains per Accept criterion, the exact next step, the
   files it owns, and anything learned that the item text lacks — the
   same text in its last commit message body, so the state survives an
   interruption. The orchestrator dispatches a successor on the same
   branch within the hour, same tier unless the remainder is mechanical.
   Two items are never sequenced through one agent (W5.11 → W5.10, 494
   requests, 12 cold restarts, 4.7 hours, is the counter-example).
2. **Five minutes, not ten.** Any command expected to exceed four minutes
   runs detached (`nohup setsid bash -c "flock /tmp/ampere-gate.lock pixi
   run -e <env> pytest <files> -q --tb=short > <log> 2>&1; echo DONE >>
   <log>" &`) and is polled with `sleep 240; tail -3 <log>` — one poll
   per turn, each a cached read; never a foreground sleep beyond the
   cache. Split a targeted run that cannot fit by file.
3. **Verify the tier on the next dispatch.** The orchestrator reads the
   first wave-5 agent's transcript (`grep -c 'ephemeral_1h_input_tokens":[1-9]'
   <transcript>` against the `5m` sibling) and records here whether
   subagents can be given the one-hour cache; if a harness setting does
   it, rule 2 becomes advice.
4. **Design documents by section.** A contract or design page (60 k
   characters each) is never read whole: the item's read-ranges are the
   default, and anything beyond them is located with `grep -n` and read
   as a range. The same for `likelihood.py` and `dataset.py`.
5. **Batch verification.** Lint, format-check, typecheck and the import
   sweep are one command; independent checks share a call; every
   command's output is trimmed (`-q --tb=short`, `tail`, `grep`) before
   it enters the context — never `cat` a log.

**Expected effect.** Rule 1 is the large one — a quarter to a third of
agent cost if the average context comes down toward 200 k; rule 2 removes
about two thirds of the 16 %; rules 4 and 5 slow the growth that rule 1
caps. Main-session cost is fine as it is. Re-measure after wave 6 with
the same method.
