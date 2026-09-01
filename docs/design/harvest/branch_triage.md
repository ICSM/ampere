# Branch triage — W0.7

Every `origin/*` branch (excluding `origin/master`/`origin/HEAD`) as of this
harvest, enumerated with `git branch -r`. "Ahead" / "behind" are commit
counts relative to `master` (`git rev-list --count master..origin/<b>` /
`origin/<b>..master`), computed locally from already-fetched remote refs —
no network access used.

**This table is a recommendation only. No branch has been archived, tagged,
or deleted — that requires Peter's explicit approval (see AGENTS.md ground
rule 3).**

| Branch | Last commit | Ahead / behind master | Content summary | Recommendation |
|---|---|---|---|---|
| `copilot/explore-implementation-plan` | 2026-04-01 | 2 / 10 | Copilot coding-agent's own `IMPLEMENTATION_PLAN.md` + a first-draft `ampere/core`/`ampere/backends` scaffold (`Parameter`/`ParameterSet`, `Model`, `Data` ABCs; a legacy-backend re-export shim; a stubbed jax backend). | **Archive-tag.** Fully harvested to `docs/design/harvest/copilot_scaffold/`; it is the designated seed for W1.3. No reason to keep the branch active — tag for provenance, then it can be deleted. |
| `datamasking` | 2021-06-11 | 0 / 357 | Likelihood data-masking behaviour + a "relative abundances" variant of the `PowerLawAGN` model. | **Delete.** Zero commits ahead of `master` — every commit here is already an ancestor of `master`; content is fully subsumed. |
| `dynesty` | 2022-07-13 | 0 / 205 | Dynesty-sampler default-parameter fixes for a `Cstar` model, typo fixes. | **Delete.** Zero ahead; the dynesty backend itself lives on `master` as `ampere/infer/dynestysearch.py` already. |
| `fastresampling` | 2021-10-23 | 0 / 335 | Fast spectral resampling using `numpy.interp` instead of the slower default. | **Delete.** Zero ahead; fully subsumed by `master`. |
| `folder_selection` | 2023-07-31 | 0 / 21 | A Cstar SBI example + the switch to `pyproject.toml`/dynamic versioning. | **Delete.** Zero ahead; `master` already uses `pyproject.toml`. |
| `jax` | 2026-01-27 | 9 / 12 | Early design sketches for `ampere/jax/` and `ampere/torch/` backends (`Module`/`Data`/`InstrumentModel`/`Transformation`/`Likelihood` skeletons, several 0-byte placeholders). The final (2026-01-27) commit also added `examples/lightning_logs/` — ~50 MB of PyTorch-Lightning checkpoint binaries. | **Archive-tag, then consider a filtered re-push or delete.** Design sketches harvested to `docs/design/harvest/jax/` (annotated: does not import, must not be merged as-is). The branch itself is bulky (~50 MB of binaries never meant to be committed) — archive-tag preserves history for reference without keeping it live; deleting afterwards would reclaim remote storage, pending approval. |
| `mpire` | 2022-01-04 | 1 / 332 | Attempt to parallelise likelihood evaluation in `minimal_working_example.py` using the `mpire` package. Commit message: "Not currently working". | **Delete.** Single dead-end commit, explicitly marked broken by its own author; nothing to preserve. |
| `mwe` | 2021-10-18 | 0 / 336 | Dynesty functional-output fix + a minimal working example. | **Delete.** Zero ahead; fully subsumed by `master`. |
| `optim_only` | 2023-08-08 | 12 / 35 | `ampere/infer/optim.py` (1,226 lines): scipy/Ax/Ray-Tune point-estimate optimiser backends (`ScipyMinOpt`, `ScipyDE`, `AxOpt`, `AxSAASBO`, `AxBOParallel`, …). | **Archive-tag.** Fully harvested to `docs/design/harvest/optim_only/`. Left mid-debug (final commit: "still figuring out some issues" with SAASBO); reviving it is a Phase-2+ decision, not urgent to keep the branch live for. |
| `pyphot-integration` | 2022-05-24 | 0 / 261 | Integrating the `pyphot` package for photometric-filter handling, plus logging improvements. | **Delete.** Zero ahead; `pyphot` integration is live in `ampere/data/photometry.py` on `master` already. |
| `refactor_inference` | 2022-08-09 | 0 / 189 | SBI posterior-predictive caching work ("Now with caching working properly..."). | **Delete.** Zero ahead; superseded by `master`'s own SBI caching fixes (e.g. commit `9a3d896`, "Fixing plotting when SBI models are not cached"). |
| `sbi_embeddings` | 2023-07-23 | 0 / 35 | Embedding networks for SBI, with documentation and review fixes. | **Delete.** Zero ahead; fully subsumed by `master`. |
| `small_silicates` | 2022-08-17 | 10 / 306 | Personal science-analysis branch: crystalline-silicate line fitting for CygOB2 sources. Diverged from a pre-refactor point in history — includes now-superseded duplicate copies of `ampere/infer/*.py` (`basesearch.py`, `dynestysearch.py`, `emceesearch.py`, `mixins.py`, …) predating the current package layout, and predates `master`'s removal of the large `ampere/Opacities_Sascha/` opacity-data-file directory (tens of thousands of lines of `.q` files). | **Archive-tag; confirm with the author before deleting.** Not architecturally relevant to the v2 redesign (not part of this harvest's four required items), but may still hold value as a science record for the CygOB2 analysis — rebasing it onto current `master` would be non-trivial given how far it has diverged. |
| `swyft` | 2024-04-23 | 8 / 11 | TMNRE (Truncated Marginal Neural Ratio Estimation) implementation via the `swyft` package: `~300` new lines in `ampere/infer/sbi.py` + a new example, `minimal_working_example_tmnre.py`. | **Archive-tag.** Fully harvested to `docs/design/harvest/swyft/` (with a noted latent bug: the `swyft` lazy-import guard doesn't actually protect the new `SwyftNetwork*` class definitions). Reviving SBI/TMNRE work is out of Phase 0/1 scope. |
| `ysosilfit` | 2022-01-05 | 0 / 309 | NGC6302-specific silicate-fit model refinements (opacities list, `SpectrumNGC6302`, prior definitions). | **Delete.** Zero ahead; fully subsumed by `master`. |

## Summary

- **15** remote branches besides `master`.
- **9** are fully subsumed by `master` (0 commits ahead) → recommend **delete**: `datamasking`, `dynesty`, `fastresampling`, `folder_selection`, `mwe`, `pyphot-integration`, `refactor_inference`, `sbi_embeddings`, `ysosilfit`.
- **1** is a dead-end, explicitly-broken single commit → recommend **delete**: `mpire`.
- **5** carry content worth preserving beyond `master` → recommend **archive-tag** (not delete, pending review): `copilot/explore-implementation-plan`, `jax`, `optim_only`, `small_silicates`, `swyft`. Four of these five are the branches this work item was asked to harvest into `docs/design/harvest/`; `small_silicates` is flagged separately as it falls outside the harvest's required scope but still merits a recommendation.
- **0** are recommended to **keep active** — nothing here looks like ongoing, in-flight development as of 2026-09-01.

No tags were created, no branches were deleted, and nothing was pushed while producing this table.

## Execution record (2026-09-01, approved by Peter)

- All 15 branches archive-tagged as `archive/<branch>` and the tags pushed.
- `jax` re-pushed filtered: tip rewritten from `55e1266` to `803a780`,
  removing only `examples/lightning_logs/` (~50 MB of checkpoint binaries;
  the sole commit touching them was the tip). `archive/jax` points at the
  filtered tip, so the binaries are unreachable from any ref; GitHub
  reclaims unreachable objects on its own schedule (support can force it
  if storage matters sooner).
- 13 remote branches deleted (recoverable from their archive tags).
- Kept: `master`, `jax` (filtered), `small_silicates` (per Peter: must not
  be deleted; tagged and left live).
