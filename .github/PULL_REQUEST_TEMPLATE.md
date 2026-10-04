## Item or issue

<!-- A work item id (for example W6.17) or "Closes #123". -->

## What and why

<!-- What changed, and why. Keep to one topic. -->

## Acceptance evidence

<!-- Which tasks you ran and what they reported, with counts, for example:
`pixi run test-fast`: 1234 passed, 5 skipped. -->

## Checklist

- [ ] Tests added or updated, and `pixi run test-fast` passes
- [ ] `pixi run lint`, `pixi run format-check` and `pixi run typecheck` pass
- [ ] `pixi run docs` builds without warnings, if docs or docstrings changed
- [ ] British English in prose
- [ ] No run outputs, logs or binaries committed
- [ ] Nothing under `ampere/legacy/` changed; any change to a section 4
      contract carries a decision-log entry in `DEVELOPMENT_PLAN.md`
- [ ] `Co-Authored-By:` trailer on commits, if an agent assisted
