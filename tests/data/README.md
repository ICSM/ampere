# Test data

Small, public, redistributable data files that tests read. Every file here
states where it came from, how to verify it and why ampere may redistribute
it. A file whose basis is ever ruled insufficient is removed with one
`git rm`: every test that reads one also knows how to download it to a cache
and skips by name when offline.

## `contest-2008-binary.oifits`

| | |
|---|---|
| What | The binary-star verification data set of the 2008 Optical/IR Interferometry Imaging Beauty Contest (Cotton et al. 2008, *Proc. SPIE* 7013, 70131N): simulated CHARA/MIRC H-band data of a binary, OIFITS revision 1 — `OI_ARRAY` (six stations), `OI_TARGET`, `OI_WAVELENGTH` (eight channels, 1.50–1.75 µm), `OI_VIS2` (75 rows) and `OI_T3` (100 rows); no `OI_VIS`. |
| Size | 83 520 bytes |
| SHA-256 | `2476bd412d25ddf9af3ee7002f8998a1e1e6f9fbbfbc60310b7412310c235005` |
| Pinned URL | <https://raw.githubusercontent.com/emmt/OIFITS.jl/0978576aeb42e25fa56223853997d9ddf79c83ac/test/contest-2008-binary.oifits> |
| Known truth | Separation 5.0 mas, position angle 30° east of north from the bright to the faint component, brightness ratio 8.9, uniform-disc diameters 1.2 mas (primary) and 0.75 mas (secondary); the field tapered with a 15 mas FWHM Gaussian (the contest's `CONTENTS.md`). |
| Read by | `tests/interferometry/test_oifits.py` (W6.12) |

**Redistribution basis.** The contest data were released publicly by the
organisers of the IAU/SPIE interferometry imaging beauty contest for anyone to
use in testing image-reconstruction and model-fitting software, and that is
what they are used for here. The copy vendored here is the one redistributed in
the MIT-licensed [OIFITS.jl](https://github.com/emmt/OIFITS.jl) repository
(the pinned commit above), whose licence permits redistribution. Vendoring
was confirmed by Peter Scicluna at the W6.12 dispatch (2026-10-05).

**Verifying.** `sha256sum tests/data/contest-2008-binary.oifits` must print
the checksum above; the test checks it too, on the vendored copy and on a
downloaded one.
