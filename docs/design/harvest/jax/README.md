# Harvest: `origin/jax` — jax/torch backend design sketches

**Source branch:** `origin/jax`, tip `55e1266` ("Designing backends in torch
and jax (partial)"), 2026-01-27. Nine commits unique to this branch, all by
Peter Scicluna, in two separate sessions: eight commits from 2024-02-22 to
2024-02-23 sketching `ampere/jax/`, and a single much later commit
(2026-01-27) that revisited the branch to sketch `ampere/torch/data.py` and
tweak `ampere/jax/__init__.py`/`likelihoods.py`/`torch/module.py`
(merge-base with `master` is `060d2bf`).

**What it is:** early, deliberately minimal design sketches for a future
JAX backend (`ampere/jax/`) and a future torch backend (`ampere/torch/`),
each mirroring a `Module` / `Data` / `InstrumentModel` /
`Transformation` / `Likelihood` split. Several files are empty placeholders
(0 bytes: `ampere/jax/models.py`, `ampere/jax/noisemodels.py`,
`ampere/torch/likelihoods.py`, `ampere/torch/models.py`) that exist only to
reserve the module name. The docstrings describe an intended API (e.g.
`ampere.jax` as a "drop-in replacement" importable as `import ampere.jax as
ampere`) well beyond what any file actually implements.

**Files in this harvest** (mirrors `ampere/jax/` and `ampere/torch/` from
the source branch; every file carries a prepended banner comment — see
below):
- `ampere/jax/__init__.py`, `data.py`, `instrumentmodels.py`,
  `likelihoods.py`, `models.py`, `module.py`, `noisemodels.py`,
  `transformations.py`
- `ampere/torch/__init__.py`, `data.py` (source branch filename was
  `ampere/torch/data,py` — a comma typo — renamed here to `data.py` for
  clarity), `likelihoods.py`, `models.py`, `module.py`

**EXCLUDED from this harvest, deliberately:** the branch's final commit
(2026-01-27) also added `examples/lightning_logs/version_{0..7}/` —
8 PyTorch-Lightning checkpoint files (`*.ckpt`, ~6.5 MB each, ~50 MB total)
plus their `hparams.yaml` siblings. These are training-run binaries, not
design content, are already covered by the repository's `lightning_logs/`
`.gitignore` rule, and are **not present anywhere in this harvest** (see
acceptance evidence in the W0.7 report — `find … -size +200k` and a binary
scan both come up empty).

**CRITICAL — status, read before touching this code:**

Every file above carries a prepended, `grep`-able banner
(`ASPIRATIONAL SKETCH`, `DOES NOT IMPORT`, `MUST NOT BE MERGED AS-IS`).
Summary of why:

- `ampere/jax/module.py` does `from equinox import Dataset as EqxDataset`
  and `from equinox import Parameter as EqxParameter`. **Neither
  `equinox.Dataset` nor `equinox.Parameter` exists** in the equinox public
  API — this import fails immediately. Every other `ampere/jax/*.py` file
  imports `Module` from `.module` (directly, or transitively via
  `.transformations`), so the entire `ampere.jax` sketch package fails to
  import, full stop.
- `ampere/torch/data.py` does `from .likelihoods import Likelihood`, but
  `ampere/torch/likelihoods.py` in this same sketch is a 0-byte empty file
  with no `Likelihood` symbol — that import fails too.
  `ampere/torch/module.py` additionally depends on `gpytorch` and `pyro`,
  neither of which is declared anywhere in this project's dependencies or
  extras.
- These are reference sketches of intended structure/naming for a possible
  future JAX/torch backend (see `DEVELOPMENT_PLAN.md` §3's backend
  discussion) — **not working code, not a starting point to build from
  directly, and must never be merged as-is.** Any future JAX backend work
  should treat this only as prior design intent to consult, the same way
  `docs/design/prior_art.md` (W1.1) treats external libraries.
