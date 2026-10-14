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

## `jw01536-o027_t008_miri_p750l_x1d.fits`

| | |
|---|---|
| What | A JWST MIRI low-resolution-spectrometer fixed-slit spectrum (`EXP_TYPE` `MIR_LRS-FIXEDSLIT`, detector `MIRIMAGE`, prism `P750L`) of the A-type flux standard BD+60 1753: programme 1536 (an absolute-flux-calibration programme, public with no exclusive-access period), observation 027, visit 001, observed 2022-07-08; the level-3 (`calwebb_spec3`) `x1d` product, calibrated with pipeline `CAL_VER` 3.0.0 and context `jwst_1584.pmap`. One `EXTRACT1D` table (`SRCTYPE` `POINT`, `SLTNAME` `MIR_LRS-FIXEDSLIT`) of 388 rows, 3.34–14.02 µm, written from long to short wavelength; `WAVELENGTH` in µm, `FLUX`/`FLUX_ERROR` in Jy, `SURF_BRIGHT`/`SB_ERROR` in MJy/sr; and the pipeline's `ASDF` extension. |
| Size | 195 840 bytes |
| SHA-256 | `16407b30bbf4bf86837e894cc39b702bf17207a0b7ca7bd56c26d7f23efba886` |
| Pinned URL | <https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:JWST/product/jw01536-o027_t008_miri_p750l_x1d.fits> (MAST's portal download, no token needed) |
| Known truth | The file's own structure: the first and last rows are padding (`DQ` = 1, `DO_NOT_USE`; a `NaN` flux; a zero uncertainty), the last repeating the wavelength of the row before it, and no other `DQ` bit is set; the flux falls from about 32 mJy at 3.3 µm to about 3.5 mJy at 14 µm, the Rayleigh–Jeans tail of an A star. |
| Read by | `tests/spectroscopy/test_jwst.py` (W7.3) |

**Redistribution basis.** JWST data are produced by NASA, ESA and CSA and
served by MAST at STScI; once a programme's exclusive-access period lapses
(calibration programmes have none) they are public. MAST's data-use policy
(<https://archive.stsci.edu/publishing/data-use>) states that "most data hosted
at MAST are in the public domain … and therefore do not have restrictions on
use", the copyrighted exceptions being the Digitized Sky Survey and the Guide
Star Catalogs, not JWST; it asks only that publications acknowledge the
mission and STScI. Vendoring a product under a megabyte was confirmed by Peter
Scicluna at the W7.3 dispatch (2026-10-09).

**Verifying.** `sha256sum tests/data/jw01536-o027_t008_miri_p750l_x1d.fits`
must print the checksum above; the test checks it too, on the vendored copy and
on a downloaded one.
