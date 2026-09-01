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
| Codex + gpt-5.6-sol | Frontier OpenAI (separate budget) | Detailed and adversarial reviews — see the `sol-review` skill. Cross-provider review has uncorrelated failure modes with Claude-written code, which is the point. |
| Codex + luna | Cheap OpenAI (separate budget) | Bulk mechanical work when preserving Claude quota, or independent second drafts of small pieces. |

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

- Remaining Phase 0 (W0.5, W0.6, W0.8): Sonnet.
- Phase 1 contract specs (W1.2–W1.9): Fable-drafted (or Fable-finalised),
  sol-reviewed adversarially at the W1.13 freeze; mechanical companions
  (prior-art memo W1.1, modality sketches W1.11) to Sonnet/Opus with
  Fable review.
- Phase 2 backend tracks: Opus per track (torch, jax), conformance suite
  as the cross-check, sol adversarial review at milestone M2.
- Reviews of merged-candidate PRs: orchestrator pass always; sol pass for
  contract/maths-touching changes.
