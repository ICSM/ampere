---
name: delegate-codex
description: Delegate a well-scoped, self-contained subtask to the OpenAI Codex CLI non-interactively and review its output. Use when parallel independent implementation work would help, or when explicitly asked to use Codex.
---

# Delegating to Codex CLI

The Codex CLI is installed (`codex`, v0.149+; auth is the user's — if it
demands login, stop and report rather than authenticating yourself).
Codex automatically reads `AGENTS.md`, so the ground rules travel with the
repo — but your prompt must still carry the task specifics.

## When to delegate

Good: mechanical, well-specified, verifiable subtasks (port this module per
this spec; write tests for X against these criteria). Bad: anything needing
project judgement, contract interpretation, or touching frozen legacy code.

## How

1. Prefer an isolated worktree so Codex cannot disturb your checkout:
   `git worktree add /tmp/codex-<task> <base-branch>`
2. Run non-interactively, sandboxed to the workspace:

   ```bash
   codex exec --cd /tmp/codex-<task> --sandbox workspace-write \
     "<self-contained task: goal, files, constraints, acceptance criteria, and 'do not commit'>"
   ```

   (Flags drift between versions — check `codex exec --help` if that
   errors. Never use a danger/full-access sandbox mode.)
3. **Review the entire diff yourself** (`git -C /tmp/codex-<task> diff`)
   before adopting any of it: correctness, scope creep, ground-rule
   violations, British English in prose. You own what you merge.
4. Adopt by committing the reviewed changes on your own branch (with the
   standard trailer), then `git worktree remove /tmp/codex-<task>`.
5. In your report, state what was delegated to Codex and that you reviewed
   the diff.

Codex must never commit to shared branches, push, or touch anything outside
its worktree.
