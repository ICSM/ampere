# Harvest: `origin/swyft` — TMNRE (Swyft) simulation-based inference

**Source branch:** `origin/swyft`, tip `3f92195` ("Update gitignore for
swyft"), 2024-04-23. Eight commits unique to this branch, all by Peter
Scicluna, spanning 2024-02-23 to 2024-04-23 (merge-base with `master` is
`9a3d896`, the "Fixing plotting when SBI models are not cached" commit).

**What it is:** an implementation of Truncated Marginal Neural Ratio
Estimation (TMNRE) using the [swyft](https://github.com/undark-lab/swyft)
package, added on top of the existing `sbi`-package-based simulation-based
inference in `ampere/infer/sbi.py`. It adds `~300` lines to that file —
`SwyftNetwork` / `SwyftNetworkLinearEmbedding` (thin `swyft.SwyftModule`
wrappers) and a new `Swyft_TMNRE(LFIBase)` inference class implementing
`optimize()`/simulate/train/infer using `swyft.SwyftTrainer` — plus a new
example, `examples/minimal_working_example_tmnre.py` (257 lines), showing
end-to-end usage against a synthetic AGN power-law model.

**Files in this harvest:**
- `swyft-tmnre.patch` — `git diff master...origin/swyft` restricted to
  `ampere/infer/sbi.py` and `examples/minimal_working_example_tmnre.py`
  (the requested diff-only view).
- `sbi.py` — full file as of `origin/swyft` tip (`git show
  origin/swyft:ampere/infer/sbi.py`); includes the pre-existing
  `sbi`-package classes (`LFIBase`, `SBI_SNPE`, …) already on `master`, plus
  the new Swyft/TMNRE additions.
- `minimal_working_example_tmnre.py` — full file as of the branch tip.

**Status / caveats:**
- `swyft` is imported behind a `try`/`except (ModuleNotFoundError,
  ImportError)` guard at module top (`is_swift_installed`, sic), matching
  the existing lazy-import pattern for optional inference backends. However,
  `class SwyftNetwork(swyft.SwyftModule):` and
  `class SwyftNetworkLinearEmbedding(swyft.SwyftModule):` reference the bare
  `swyft` name unconditionally at class-definition time (module import
  time), so if `swyft` is not installed the guard does **not** prevent a
  `NameError` on import — the lazy-import guard is incomplete. This is a
  latent bug in the branch, not something introduced by this harvest.
- The base `sbi` package import (`from sbi.inference import SNPE,
  DirectPosterior`) is *already* unconditional on `master`'s
  `ampere/infer/sbi.py` — that is a pre-existing pattern, not something this
  branch changed.
- Not validated against a current `swyft` install; the branch predates the
  W0.5 quarantine/lazy-import work and has not been characterisation-tested.
- Never merged; `ampere/infer` is frozen legacy, so any revival of this work
  is a Phase-2+ decision (SBI is out of Phase 0/1 scope per
  `DEVELOPMENT_PLAN.md`).
