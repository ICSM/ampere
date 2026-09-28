# `examples/`

This directory holds two generations of ampere example code side by side.

**v2 example packages**, each written against `ampere.core` and a modern
backend, each covered by its own test suite in `tests/examples/` so that none
can rot unnoticed: `sed_composition`, `m2_misspecification`, `interferometry`,
`astrometry`, `image`, `wstat_comparison.py`, `sbi/`, `linear_sed` and
`modified_blackbody` — see `docs/source/tutorials.rst` for the walk-throughs.
`population` has no example package of its own; its code is small enough to
live in the tutorial page inline.

**Legacy examples**, kept exactly as they are and run on legacy
(`ampere.legacy`) rather than touched — ampere v2's decision D1 (b). Six of
them have gained a v2 **twin**: the same model, the same data, the same
question, written once against `ampere.core` and the reference backend,
beside the original (W6.13).

| Legacy script(s) | Twin package | Status |
| --- | --- | --- |
| `minimal_working_example.py` and its `_dynesty`, `_zeus`, `_sbi`, `_sbi_embedding` variants | `examples/linear_sed/` | landed W6.13 (1) |
| `NGC6302.py`, `NGC6302_zeus.py`, `NGC6302-calculate-dust-mass.py` | `examples/ngc6302/` | twin pending, W6.13 (2) |
| `examples_paper/modifiedblackbody.py` | `examples/modified_blackbody/` | landed W6.13 (3) |
| `examples_paper/phoenixstar.py` | `examples/phoenix_star/` | twin pending, W6.13 (4) |
| `cstar_model_test_sbi_v2.py` and its `_embedding` variant | `examples/cstar/` | twin pending, W6.13 (5) |
| `examples/star_disc.py` (the HD105 SED) | `examples/star_disc/` | twin pending, W6.13 (6) |

`minimal_working_example.py` and its four variants are
**characterisation-test anchors**: `tests/characterisation` exercises them
directly, on legacy, to define "legacy still works" for the v2 redesign, and
they are not touched by `linear_sed`'s own work.
`examples_paper/modifiedblackbody.py` has no characterisation test of its
own; it too is kept byte-identical, untouched by `modified_blackbody`'s
work.

**Three scripts have no twin and never will**: `example.py`, `modbbtest.py`
and `modelClio.py` each import a path that no longer exists
(`ampere.emceesearch`, `ampere.PowerLawAGN` — neither is `ampere.infer.
emceesearch` or anything else legacy currently ships), so none of the three
has run in a long time. Ruled D12 (a): dead code, not migrated, not fixed.

**Notebooks** (`docs/source/notebooks/`: `quickstart`, `Ampere_MBB_Example`,
`Embedding_nets`) teach the legacy v1 API and are rendered from stored
output rather than executed by the documentation build; converting them is
W6.1's, not this item's.
