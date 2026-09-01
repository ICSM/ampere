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
