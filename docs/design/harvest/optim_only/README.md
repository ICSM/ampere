# Harvest: `origin/optim_only` — optimiser-based inference module

**Source branch:** `origin/optim_only`, tip `6f0821e` ("Trying to implement
SAASBO - still figuring out some issues like how to pass the function to Ax
in this interface."), 2023-08-08. Twelve commits unique to this branch, all
by Peter Scicluna, spanning 2023-04-12 to 2023-08-08 (merge-base with
`master` is `4f0b1b5`).

**What it is:** a single new module, `ampere/infer/optim.py` (1,226 lines,
harvested in full — the branch touches no other file), adding a family of
point-estimate/optimisation-based "inference" backends alongside the
existing MCMC/nested-sampling searches: `OptimBase(BaseSearch, Logger)` as
the shared base, then thin wrappers around `scipy.optimize` routines
(`ScipyMinOpt`, `ScipyBasinOpt`, `ScipyDE`, `ScipyDualAnneal`, `ScipyShgo`,
`ScipyDirect`), an [Ax](https://ax.dev/) Bayesian-optimisation interface
(`AxOpt`, `AxBO`, `AxSAASBO` for sparse-axis-aligned-subspace BO,
`AxBOParallel`), and a Ray Tune-driven parallel BO variant
(`AxBOParallel` uses `from ray.tune import report` / `ray.tune.suggest.ax`
as lazy inline imports inside the method body). Extensive docstrings and
citations were added in the penultimate commit.

**Files in this harvest:**
- `optim.py` — full file as of `origin/optim_only` tip (`git show
  origin/optim_only:ampere/infer/optim.py`).

**Status / caveats:**
- The final commit message ("Trying to implement SAASBO - still figuring
  out some issues...") indicates the branch was left mid-debugging; `AxOpt`
  / `AxSAASBO` / `AxBOParallel` should be treated as work-in-progress, not
  finished, tested code.
- `scipy` is a core ampere dependency already; `ax-platform` and `ray[tune]`
  are not declared anywhere in `pyproject.toml` (including the `W0.3`
  extras skeleton) — reviving this module would need new optional extras.
- Not characterisation-tested; predates the W0.4/W0.5 hygiene work.
- Never merged; `ampere/infer` is frozen legacy, so reviving any of this is
  a Phase-2+ decision.
