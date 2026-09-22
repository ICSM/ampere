---
name: work-item
description: Pick up and execute an ampere v2 work item (W0.x/W1.x) from WORK_ITEMS.md following the project working agreement. Use whenever asked to work on a numbered work item.
---

# Executing an ampere v2 work item

## Before writing anything

1. Read `AGENTS.md` (ground rules), then your item's entry in
   `WORK_ITEMS.md`, then every `DEVELOPMENT_PLAN.md` section your item
   references. Check the item's **Depends** line — if a dependency is not
   merged, stop and report that instead of working around it.
2. Create a branch from the base your dispatcher named (default: current
   HEAD): `w<item>-<slug>`, e.g. `w0.3-packaging-floor`.

## While working

- Scope: only what the item says. Discovered problems go in your final
  report, not in the diff.
- Verify claims against the actual code — the plan describes intent, the
  repository is the ground truth. If they conflict, say so in your report.
- Never touch frozen legacy modules unless the item explicitly says so;
  never commit run outputs or binaries; British English in prose.

## Budget (the token-economy rules, `docs/orchestration.md`; in force 2026-09-22)

- **Five-minute rule.** Your prompt cache lives five minutes. Any command
  expected to exceed four minutes (a targeted test file, a docs build)
  runs detached and is polled, never waited on in the foreground:
  `nohup setsid bash -c "flock /tmp/ampere-gate.lock pixi run -e <env> pytest <files> -q --tb=short > <log> 2>&1; echo DONE >> <log>" &`
  then, one per turn, `sleep 240; tail -3 <log>`. Never run a five-suite
  gate; that is the orchestrator's, on master, after the merge.
- **Soft lifetime cap, then hand off.** Past **250 k of context** (the
  harness's indicator when shown), **about 300 tool uses, or 2.5 hours**,
  whichever first: do not start a new unit of work. Finish the unit in
  hand — never mid-edit; a clearly labelled WIP commit if the boundary
  cannot be reached — commit, and return your report with a **Hand-off**
  section (below); put the same text in that last commit's message body.
  A successor agent continues on your branch.
- **Read by section.** Design documents and contract pages, and
  `likelihood.py` / `dataset.py`, are read as the ranges your prompt names
  or as ranges you locate with `grep -n`; never whole. Use `git diff`
  rather than re-reading a file you edited.
- **Batch and trim.** Lint, format-check, typecheck and the import sweep
  in one command; independent checks in one call; trim every output
  (`-q --tb=short`, `tail`, `grep`) — never `cat` a log.

## Definition of done

1. Every **Accept** criterion executed, with evidence (command + output)
   captured for the report.
2. Work committed on your branch in logical commits, messages explaining
   why, ending with the Co-Authored-By trailer (see AGENTS.md). Do not
   push unless the dispatcher said to.
3. Do not edit the WORK_ITEMS.md status table (merge-conflict magnet);
   status is updated at merge.

## Report format (returned to dispatcher / PR description)

- **Item**: id + one-line goal.
- **What changed**: files + why, in complete sentences.
- **Acceptance evidence**: per criterion, what you ran and what it showed.
- **Out-of-scope findings**: anything broken or suspicious you did NOT fix.
- **Open questions**: decisions you deferred to review.
- **Hand-off** (only when stopping under the lifetime cap): what is
  done; what remains, per Accept criterion; the exact next step; the
  files you own; anything learned that the item text lacks.
