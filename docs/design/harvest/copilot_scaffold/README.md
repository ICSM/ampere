# Harvest: `origin/copilot/explore-implementation-plan` — core/backends scaffold

**Source branch:** `origin/copilot/explore-implementation-plan`, tip
`b7c62fe` ("feat: IMPLEMENTATION_PLAN.md + ampere/core and ampere/backends
scaffold"), 2026-04-01. Two commits, both by `copilot-swe-agent[bot]`, both
on 2026-04-01 (merge-base with `master` is `8d925b4`, i.e. this branch was
cut from a point very close to the current `master` tip — it is only 10
commits behind and 2 ahead).

**What it is:** a GitHub Copilot coding-agent exploration of the same
backend-neutral-core redesign this v2 effort is now doing independently. It
adds `IMPLEMENTATION_PLAN.md` (415 lines — Copilot's own phased plan and
issue breakdown, written without knowledge of this project's
`DEVELOPMENT_PLAN.md`/`WORK_ITEMS.md`) and a first-draft scaffold:
`ampere/core/parameter.py` (`Parameter` dataclass + `ParameterSet`:
named/prior-equipped parameters, pack/unpack, `lnprior`, `prior_transform`,
a `merge()` helper for composite models, and a stubbed `as_paramax()`
extension point), `ampere/core/model.py` (`Model` ABC), `ampere/core/data.py`
(`Data` ABC), `ampere/core/__init__.py`, plus `ampere/backends/__init__.py`,
`ampere/backends/legacy.py` (re-exports the existing `ampere.infer`
samplers under a stable path), and `ampere/backends/jax/__init__.py` (a
deliberate stub that raises `ImportError`/`NotImplementedError`).

**Files in this harvest:** exact copies, unmodified, of the eight files
above, under `ampere/` mirroring their source paths, plus
`IMPLEMENTATION_PLAN.md`.

**Status / caveats:**
- Unlike the `origin/jax` sketches, `ampere/core/parameter.py`,
  `model.py`, and `data.py` **do import successfully standalone** (checked
  with `importlib` in isolation against this repo's Python 3.13 env: only
  `numpy`/`dataclasses`/`abc`/stdlib `typing` are required at import time;
  `as_paramax()` lazily requires `paramax` and currently always raises
  `NotImplementedError`). `ampere/backends/legacy.py` imports the real
  `ampere.infer.{emceesearch,dynestysearch,zeussearch,nestedsearch}`
  classes, so it depends on legacy `ampere.infer` importing cleanly.
- This is a *scaffold*, not a finished contract: no unit tests, no
  docstring/behaviour review against `DEVELOPMENT_PLAN.md` §4, no
  `pyrefly` typing pass, and it was written by an agent working from a
  self-generated plan rather than this project's actual
  `DEVELOPMENT_PLAN.md`/`WORK_ITEMS.md`. Interfaces it assumes (e.g.
  `ampere.models.results.ModelResults`, `ParameterSet.as_paramax()`) are
  forward references (`TYPE_CHECKING`-guarded) to things that do not yet
  exist on `master`.
- Per `WORK_ITEMS.md` W1.3 ("Parameter & prior contract"), this scaffold —
  `ampere/core/parameter.py` in particular — is the designated **starting
  point** for that work item, not a drop-in final implementation: W1.3 must
  still produce `docs/design/contracts/parameters.md`, add unit tests, pass
  `pyrefly`, and reconcile the design against W1.2's architecture spec
  (transforms to unconstrained space, tying/sharing, plate-aware groups,
  units, buffer declaration — none of which this scaffold implements yet).
- Never merged; the branch remains unreviewed by Peter.
