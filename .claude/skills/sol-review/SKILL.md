---
name: sol-review
description: Run a detailed or adversarial cross-model code/spec review using the Codex CLI with gpt-5.6-sol. Use when asked for a sol review, adversarial review, or second opinion — and before merging changes that touch the core contracts, likelihood/GP mathematics, or lowering rules.
---

# Cross-model review with gpt-5.6-sol

Reviews by a different model family have failure modes uncorrelated with
the (usually Claude-written) code under review — that is the value. Sol is
a frontier model on the user's OpenAI plan: use it for high-stakes reviews,
not routine lint.

## Invocation

Always read-only; never let a reviewer modify anything:

```bash
codex exec -m gpt-5.6-sol --sandbox read-only --cd <repo-or-worktree> \
  "<review prompt — see templates>"
```

- Scope the target precisely in the prompt: branch/diff range
  (`git diff master...<branch>`), or named files/sections for spec review.
- If the model name is rejected, check `codex exec --help` and available
  models, and report to the dispatcher — do not silently downgrade to a
  different model.
- `codex exec review` is a built-in alternative flow; prefer the explicit
  prompt templates below for adversarial work so the instructions are
  controlled.

## Template — detailed review

> Review the changes in `<range>` for the ampere project. Ground rules and
> architecture: AGENTS.md and DEVELOPMENT_PLAN.md (§4 contracts, §7 traps).
> Assess: correctness; contract compliance; edge cases (empty/masked data,
> irregular grids, float32 vs float64, units); test adequacy against the
> work item's acceptance criteria in WORK_ITEMS.md; scope creep; touches to
> frozen legacy modules. Report findings as: file:line — claim — concrete
> failure scenario (inputs → wrong behaviour) — severity. Label each
> CONFIRMED (you traced the failing path) or SPECULATIVE. British English.

## Template — adversarial review

> Adversarial review of `<range>` in the ampere project. Assume the code
> contains at least one serious defect and your job is to find it. Attack:
> mathematical correctness (likelihoods, GP algebra, determinants,
> normalisation, log-space handling); boundary behaviour (N=0/1, all-masked
> data, single-channel defaults, irregular/unsorted coordinates); dtype and
> precision (float32 Cholesky, jax x64); hidden assumptions (regular grids,
> sorted wavelengths, positive fluxes); contract violations against
> DEVELOPMENT_PLAN.md §4; failure signalling (NaN/crash paths); and
> whether the tests would actually catch the defects you hypothesise.
> Generate at least 10 candidate attack vectors, pursue each, then
> self-filter to the ones you can support with evidence. Report as:
> file:line — claim — concrete failure scenario — severity — CONFIRMED or
> SPECULATIVE. An empty verified-findings list is an acceptable outcome;
> do not pad.

## Handling the output

1. Findings are **advisory**. Verify every CONFIRMED finding against the
   code yourself before acting; spot-check SPECULATIVE ones worth the
   time. Never apply a reviewer's fix without understanding it.
2. Record the review and each finding's disposition (fixed / rejected,
   with reason) in the PR description.
3. Disagreement between sol and the orchestrator's own review is signal,
   not noise — escalate to Peter rather than quietly picking one.
