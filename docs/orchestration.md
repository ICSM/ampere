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
