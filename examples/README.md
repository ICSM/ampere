# `examples/`

This directory holds two generations of ampere example code side by side.

**v2 example packages**, each written against `ampere.core` and a modern
backend, each covered by its own test suite in `tests/examples/` so that none
can rot unnoticed: `sed_composition`, `m2_misspecification`, `interferometry`,
`astrometry`, `image`, `wstat_comparison.py`, `sbi/`, `linear_sed`,
`modified_blackbody`, `phoenix_star`, `star_disc` and `cstar` — see `docs/source/tutorials.rst` for the walk-throughs.
`population` has no example package of its own; its code is small enough to
live in the tutorial page inline.

`rhmf_trial` is W6.9's exploratory trial of pre-fit robust matrix factorisation
(`pixi run -e rhmf python -m examples.rhmf_trial --quick --out DIR`): a trial
behind the non-default `rhmf` extra, not a feature; its findings are in
`docs/design/contracts/diagnostics.md` §7.

**Legacy examples**, kept exactly as they are and run on legacy
(`ampere.legacy`) rather than touched — ampere v2's decision D1 (b). Six of
them have gained a v2 **twin**: the same model, the same data, the same
question, written once against `ampere.core` and the reference backend,
beside the original (W6.13).

| Legacy script(s) | Twin package | Status |
| --- | --- | --- |
| `minimal_working_example.py` and its `_dynesty`, `_zeus`, `_sbi`, `_sbi_embedding` variants | `examples/linear_sed/` | landed W6.13 (1) |
| `NGC6302.py`, `NGC6302_zeus.py`, `NGC6302-calculate-dust-mass.py` | `examples/ngc6302/` | landed W6.13 (2) |
| `examples_paper/modifiedblackbody.py` | `examples/modified_blackbody/` | landed W6.13 (3) |
| `examples_paper/phoenixstar.py` | `examples/phoenix_star/` | landed W6.13 (4), with ampere's own PHOENIX emulator (a PCA plus one GP per weight over 156 PHOENIX-ACES spectra, committed as `phoenix_emulator.npz`) in place of Starfish |
| `cstar_model_test_sbi_v2.py` and its `_embedding` variant | `examples/cstar/` | landed W6.13 (5), Hyperion kept as the physics with `miepython` in place of the unpackaged `bhmie` step; needs the pixi `hyperion` environment (`pixi install -e hyperion`), an example-only requirement absent from CI |
| `examples/star_disc.py` (the HD105 SED) | `examples/star_disc/` | landed W6.13 (6), `QuickSED` as a v2 `Model` on the same emulator |

`minimal_working_example.py` and its four variants are
**characterisation-test anchors**: `tests/characterisation` exercises them
directly, on legacy, to define "legacy still works" for the v2 redesign, and
they are not touched by `linear_sed`'s own work.
`examples_paper/modifiedblackbody.py` has no characterisation test of its
own; it too is kept byte-identical, untouched by `modified_blackbody`'s
work. `NGC6302-calculate-dust-mass.py`'s post-processing (equations 4 and 5
of Kemper et al. 2002) lives on as `examples/ngc6302/dust_mass.py`, a
function over a fit's posterior draws rather than a script with hard-coded
inputs; the original script itself is untouched, kept exactly as it is
beside its twin.

**Three scripts have no twin and never will**: `example.py`, `modbbtest.py`
and `modelClio.py` each import a path that no longer exists
(`ampere.emceesearch`, `ampere.PowerLawAGN` — neither is `ampere.infer.
emceesearch` or anything else legacy currently ships), so none of the three
has run in a long time. Ruled D12 (a): dead code, not migrated, not fixed.

**Notebooks** (`docs/source/notebooks/`: `quickstart`, `Ampere_MBB_Example`,
`Embedding_nets`) teach the legacy v1 API and are rendered from stored
output rather than executed by the documentation build; converting them is
W6.1's, not this item's.
