# Walkthrough scripts — the evidence behind `user_journeys_memo.md`

One script per persona, run on the beta with
`pixi run --frozen -e dev python docs/design/walkthroughs/persona_<x>.py`
from the repository root. Each prints a timestamped line per step, the
summary tables, and every refusal it meets; the memo's appendices quote
the numbers. Outputs (plots, netCDF runs) go under `$AMPERE_WALK_OUT`,
default `/tmp/ampere-walkthroughs` — never into the repository (ground
rule 7). These are measurements, not examples: they take no care to be
good science, and several deliberately do the wrong thing to see what the
package says.

- `persona_a.py`, `persona_a2.py` — the SED fitter: nine-band photometry,
  a modified blackbody, emcee from the prior and from the optimiser,
  dynesty, the plots, three unit and wavelength probes (Appendix A).
- `persona_b.py` — the spectroscopist: a real IRS spectrum read by hand,
  a power law, the flexible likelihood, the diagnostics (Appendix B).
- `persona_d.py` — the interferometrist: the contest OIFITS file through
  the reader and the binary fit (Appendix D).
- `persona_e.py` — the population person: twenty objects as a
  `Population`, then the same twenty as the catalogue loop a user writes
  today (Appendix E).
- Persona H runs the astrometry example as shipped
  (`python -m examples.astrometry --walkers 24 --steps 400 --burn-in 100`).
