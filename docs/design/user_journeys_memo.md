# User journeys — who uses ampere, how, and for what; what the beta gives each of them

Status: **assessment, drafted 2026-10-06 by Fable from Peter's answers of the
same day** (§1, recorded verbatim in substance) and the beta's own surfaces
(§2); the walkthroughs in the appendix are run on the code and are the
evidence. Its product is §6: the horizon's items re-ranked by the users
they serve rather than by the code's gaps, and the few items this view
adds. Companion to `horizon_beyond_phase7.md`; both fold into
`DEVELOPMENT_PLAN.md` §5 when Peter rules.

## 1. Peter's answers (2026-10-06)

The orchestrator asked seven groups of questions; the answers, condensed
but complete.

**Who.** About ten users are involved in development, and Peter is all of
the personas. Shares (a user may be several): **A, the SED fitter, 60 %**;
B, the spectroscopist, 40 %; C, the radiative-transfer modeller, 40 %;
E, the population or survey person, 40 %; G, the legacy user migrating,
30 %; D, the interferometrist, 20 %; F, the method developer, 20 %; H, the
time-domain or astrometry person, 20 %. Some personas combine naturally
into one example: **SED fitting and populations, with part of the SAGE AGB
sample or the sample of 2026MNRAS.545f2221M, which would also cover
radiative-transfer models.**

**How.** Effectively every user writes Python; a command line is probably
not useful. Most work in notebooks; Peter encourages scripts for
pipelining and reproducibility. Smaller problems run on laptops and
workstations; nearly everyone has a cluster for larger jobs, about 30 %
with GPUs; problems and results are shared with collaborators. Sizes span
the whole range: two parameters and one object at one end, hundreds of
parameters and millions of objects at the other.

**Where they start.** Users arrive from the paper. About half read the
concept page and half go straight to an example. Not enough real-user data
yet to know the trends.

**What they need out.** Must-have: residual and misspecification
diagnostics; a machine-readable summary table; evidence and model
comparison. Nice-to-have, because users adjust them per publication:
posterior predictive plots in publication form; a results file a referee
can rerun; derived physical quantities.

**What goes wrong.** Convergence is the single biggest problem, then units.
The rest is not tested enough to know. Users probably do not know what the
misspecification diagnostics mean, and would care if they did.

**Extend or use.** Most users, 60 % or more, will write a `Model` subclass
or use an astropy model; perhaps 30 % will write a step and 15–25 % a
kind. Radiative-transfer users will mostly wrap their codes themselves.
The analogy is numpyro, pymc and pyro users.

**Unattended and at scale.** Peter wants to loop over large catalogues as
soon as possible, expecting provenance and caching to make reruns cheap.
Not a general batch route: a dedicated example users can extrapolate to
their own clusters.

**Two further questions.** (1) Porting an existing reference-native model
to jax or torch, or using an astropy model on jax: possible, and is there a
guide? (2) Only NUTS and SVI are exposed from pyro and numpyro: should
other methods be, for instance numpyro's wrapper of jaxns?

## 2. What the beta gives each persona today

Weighted by Peter's shares. For each: the entry page and example, the
pieces the fit is made from, and where the persona has to do by hand what
the package should do. The appendix measures the top ones.

**A — the SED fitter (60 %).** Entry: the quickstart notebook (which is in
fact a spectrum-plus-photometry fit of a linear model) and
`photometry_spectra.rst`; `examples/modified_blackbody`. Pieces:
`PhotometricPoints` built from filter names, wavelengths and fluxes;
`SyntheticPhotometry.from_library` over the bundled 898-filter pyphot
library; `ModifiedBlackBody`, `BlackBody`, `PowerLaw` or a user model;
`IndependentNoise` or the GP; `EmceeEngine` or `DynestyEngine`; `summary`,
the corner, posterior-predictive and residual plots; `to_netcdf`. By
hand: the fluxes come from a catalogue the user queried themselves and the
effective wavelengths must be supplied (the appendix tests whether the
library's own are accepted); the filter names must be the library's; the
model beyond the three shipped is the user's or an example's. Blocked:
nothing for one object with a shipped model. Gaps: the catalogue reader
and filter library (horizon §5.1); the model zoo (§8).

**B — the spectroscopist (40 %).** Entry: `photometry_spectra.rst`,
`sed_composition.rst`, the M2 pages; the quickstart's IRS grids. Pieces:
`Spectrum`, `Resample`, `LSFConvolution`, `CalibrationScale`, the GP
likelihood with `QuasisepGP`, `plot_residuals` with the whiteness test and
the localisation plot. By hand: reading the file (the quickstart reads a
CASSIS IRS file in a dozen lines of astropy); the resolution curve. Blocked
until W7.3: a JWST product in one call; until W7.9: line fluxes. The
persona for whom the flexible likelihood was built, and the one Peter says
does not yet know what the diagnostics mean (§5).

**C — the radiative-transfer modeller (40 %).** Entry: `sbi.rst`,
`examples/cstar` (Hyperion), `examples/phoenix_star` (the emulator by
hand). Pieces: a `Model` subclass wrapping the code on the numpy path;
`simulate_many` under the process pool with timeouts; `SBIEngine` with
caching; or emcee and the nested samplers if the code is cheap enough. By
hand: the wrapper (Peter: they will write it themselves, which the foreign-
function item makes composable rather than replaces); an emulator if they
want gradients. Blocked: NUTS and VI on such a model (horizon §2), and
sharing one grid-based model across a population (the emulator again).

**E — the population or survey person (40 %).** Entry: `population.rst`,
`examples/population`? — no: the population example lives in
`tests/inference` and the docs page; `advanced.rst`'s section. Pieces:
`Plate`, `HierarchicalPrior`, `Population` over labelled components,
reweighting of archived fits, amortised SBI with context. By hand: the
loop over objects and its bookkeeping; the store. Blocked: a population
over each object's own GP amplitude until W7.1; nested populations;
anything at 10⁵ objects and above (horizon §4); the catalogue loop example
Peter asks for (§7 of his answers).

**G — the legacy user (30 %).** Entry: `migrating.rst` with the six twins.
Pieces: `ampere.legacy` unchanged; the twins as side-by-side material.
Blocked: nothing; the cost is reading. Not measured here.

**D — the interferometrist (20 %).** Entry: `interferometry.rst` with
`read_oifits`; `examples/interferometry`. Blocked: radial profiles and
uv-table data (horizon §3, §5.4).

**F — the method developer (20 %).** Entry: `advanced.rst` "Defining new
data types", the interferometry page as the worked template, the
conformance README. Blocked: nothing; but see §4 on porting a model, which
is the first thing a developer does.

**H — time domain and astrometry (20 %).** Entry: `astrometry.rst`,
`examples/astrometry`. Blocked: real-data entry (Gaia DR4, light curves).

## 3. The flagship example Peter named: AGB stars as a population

The combination he pointed at — SED fitting, a population, a
radiative-transfer model — is one example and it threads the three largest
personas plus the catalogue loop:

- **Data**: broadband photometry of a sample of AGB stars from SAGE (the
  Spitzer LMC survey: IRAC 3.6–8 µm, MIPS 24 µm, with 2MASS) or the sample
  of 2026MNRAS.545f2221M, read from a catalogue table into
  `PhotometricPoints` per object. This is horizon §5.1's reader and filter
  library, made concrete.
- **Model**: a dust-shell radiative-transfer grid (GRAMS-style: a
  precomputed grid over luminosity, optical depth, dust temperature and
  composition, interpolated or emulated) as a shipped model — the
  emulator route of horizon §2 (b) applied to a grid rather than a live
  code, which is the easier half of the emulator contract. A modified
  blackbody as the deliberately misspecified alternative.
- **Hierarchy**: the sample's mass-loss-rate or optical-depth distribution
  as a `Population` hyperprior; each object's GP amplitude drawn from a
  shared prior through W7.1's path populations; the non-centred
  declaration through W7.0.
- **Scale and loop**: the fits run per object on the cluster through the
  existing scripts, cached by spec hash, then reweighted as a population —
  the catalogue-loop example of Peter's answer 7, at a few hundred objects
  on the beta and at the survey's size once the store exists.
- **Outputs**: the summary table per object and for the population, the
  evidence contrast grid versus blackbody per object, the residual and
  localisation diagnostics on the objects the grid cannot explain.

The walkthrough plan (§7) builds this in stages so each stage is a real
example on its own. The open data question is in §8.

## 4. The two further questions, answered from the code

**(1) Porting a reference model to jax or torch, and astropy models on
jax.** Possible, by two routes, and documented only for modality authors.

- *Astropy on jax*: `ampere.backends.jax.from_astropy` (W4.7) computes a
  curated astropy model natively — six classes and their compound
  `+ - * /` forms — and refuses by name anything outside the table, as the
  opt-in-never-silent ruling requires; `ampere.core.from_astropy` wraps any
  astropy model as a black box on the numpy path. `astropy.rst` documents
  both. A user's own astropy model class is not in the table and cannot be,
  so for them the answer is route two.
- *A reference model to a native backend*: the inheriting pattern —
  subclass the reference model together with the backend's base, override
  the four capability flags and the arithmetic, keep the declaration
  (`compile_for`, negotiation, `from_observed`). `interferometry.rst` §8,
  `astrometry.rst` §8 and `image.rst` §9 state the rule and its clauses;
  the PHOENIX emulator example does it in thirty lines per backend by
  writing the model's arithmetic once against a small `Ops` namespace and
  supplying `jnp` or `torch` per subclass. **That is the pattern the
  package already uses for kernels**: `ampere.core.ArrayOps` (W4.5) is the
  protocol, `NumpyOps`, `JaxOps` and the torch twin the namespaces, and
  every kernel's closed form is written once. The jax astropy adapter has
  a private `_JaxOps` of its own, and the emulator example a third. So
  the mechanism exists three times and the user-facing guide zero times.
  The quickstart's `LinearModel` is written with `np.asarray` and
  `float(...)` casts, which is exactly what makes a model unportable.
- **The item this implies** (new; S code, M docs; recommended for Phase 7
  as a filler or Phase 8's first docs item): one public `ArrayOps` for
  models — the kernels' protocol widened by what a model needs (`where`,
  `interp`, `trapz`, `exp`, the dtype-and-device scalar) with the three
  namespaces exported from `ampere.core` and the two backends — and a
  guide page, "Writing a model once for three backends", whose example is
  the quickstart's model rewritten against the protocol and run under
  emcee on numpy and NUTS on jax in one notebook. The 60 % of users who
  write a `Model` subclass are its audience; the conformance row is the
  example model agreeing across the three fixtures.

**(2) Exposing more of pyro and numpyro.** `NUTSEngine` builds
`numpyro.infer.NUTS` and `pyro.infer.NUTS` directly; `VIEngine` the
autoguides; `BlackjaxEngine(method=)` is the one engine with a method
table. numpyro 0.22 offers, over the same `numpyro_model` ampere already
emits: `NestedSampler` (the jaxns wrapper — an evidence on the jax problem
without a new driver, which turns horizon §6.4's jaxns item from S to
XS), `SA` (sample-adaptive, gradient-free: the method a `ForeignModel`
without a Jacobian would run under on the native path, so it belongs with
horizon §2 (a)), `BarkerMH` (a robust gradient method, low priority),
`HMCECS` (energy-conserving subsampling, the only route to a million-object
likelihood under HMC; it needs a population plate to subsample over, so it
waits for Phase 9's scale work), and the discrete-variable kernels, which
do not apply. pyro offers little beyond HMC and a random-walk kernel.
**Recommendation**: a `NumpyroEngine(method=...)` on the blackjax pattern
with `nested`, `sa` and `barker`, S, Phase 8 beside the foreign function;
`HMCECS` recorded for Phase 9. The torch side gains nothing worth a driver.

## 5. What Peter's answers change

- **Convergence first.** It is the biggest problem and the beta's answer is
  `summary`'s R-hat and ESS columns, which a notebook user must know to
  read. Appendix A reproduces it on the simplest fit (R-hat 1.5, no
  warning) and shows the optimiser warm start curing it in 22 s. Three
  cheap items — the two below, and **the optimiser's mode as the ensemble
  engines' default start** (a behaviour change for Peter to rule): (i) a **convergence verdict** — `check_convergence(tree)`
  returning a typed verdict with the failing variables and the remedy in
  words (more steps, more walkers, a reparameterisation, the optimiser's
  start), and a loud `ResultsWarning` at emission when R-hat or ESS cross
  the thresholds, on the same guard `summary` uses for approximate
  engines; (ii) the FAQ's "Inference" section rewritten around it. S,
  Sonnet; recommended as a Phase 7 filler because every walkthrough will
  hit it.
- **Units second.** The appendix's probes show what a user sees for mJy
  fluxes and an Ångström axis; the refusals that exist are correct, and
  the gap, if any, is in the message naming the fix.
- **The must-have outputs are mostly there.** Residual and misspecification
  diagnostics: `plot_residuals` with the whiteness test, `plot_gp_localisation`,
  `gp_localisation_score`, `chi_square_pvalue`, the anomaly score. The
  summary table: `summary` (arviz's, guarded). Evidence: the nested
  engines' triple. **Model comparison across runs is the one must-have that
  is missing** — horizon §6.1's `compare` module moves from Phase 8 filler
  to a Phase 8 item.
- **The diagnostics need teaching, not more code.** Half the users read
  the concept page: it explains why the GP is there and not what the
  localisation plot shows. One page, "Reading the diagnostics", walking
  the M2 scenarios' plots with the sentence each one supports, is the
  docs item that makes persona B care. S.
- **Notebooks are the medium.** Every new example should ship as a
  notebook executed at docs build (the W6.1 pattern) with the script
  beside it, not the reverse.
- **No CLI, no problem file**: horizon §8's declarative problem file drops
  to the bottom. Sharing is by `DataTree` files and scripts, which exist.
- **The catalogue loop is an example, not a route**: horizon §4's store
  memo stays, and before it a `examples/catalogue_loop` running a few
  hundred SAGE objects through the cluster scripts with the artefact cache,
  as the first stage of §3.
- **Porting models is a guide plus one protocol** (§4 (1)); the foreign
  function stays as drafted, because the RT users will wrap their own code
  and need the composable primitive, not a wrapper.

## 6. The horizon re-ranked by users served

Weight is the share of the personas an item serves, from §1; cost from
the horizon memo. Items the beta already serves well are not listed.

| Item | Serves | Weight | Cost | Rank |
|---|---|---|---|---|
| Convergence verdict, warm start as the ensemble default, FAQ (§5, App. A 1) | all | 100 % | S | 1 |
| `ModifiedBlackBody`'s `scale` unit and photometric alignment by filter name (App. A 2–3) | A | 60 % | S + S | 1 |
| Portable-model protocol and guide (§4 (1)) | A, B, C, F — anyone writing a model | 60 % | S + M docs | 2 |
| Catalogue photometry reader and the SVO filter library (horizon §5.1) | A, E | 60 % | M | 3 |
| Model comparison module (horizon §6.1) | A, B, C | 60 % | M | 4 |
| Reading-the-diagnostics page (§5) | A, B | 60 % | S | 5 |
| The AGB-population flagship, in stages (§3) | A, C, E | 60 % | L across stages | 6 |
| Emulator memo and `EmulatedModel` (horizon §2 (b)) | C, E | 40 % | memo + L | 7 |
| `ForeignModel` and `NumpyroEngine(sa, nested)` (horizon §2 (a), §4 (2)) | C, F | 40 % | S + S | 8 |
| Catalogue-loop example with the cache (§5) | E | 40 % | S | 9 |
| JWST reader (W7.3, ruled) and `from_specutils` | B | 40 % | in Phase 7; S | — |
| Model zoo promotion (horizon §8) | A, C | 60 % | M | 10 |
| Nested populations memo and item (horizon §4) | E | 40 % | memo + L | 11 |
| Population at scale, the store (horizon §4) | E | 40 % | memo + L | 12 |
| Radial profiles (horizon §3) | D | 20 % | L | 13 |
| Correction module (horizon §6.2), pocoMC, jaxns | C, E | 40 % | M, S, XS | 14 |
| IFU cube (horizon §7) | B | 40 %, no user yet | L | 15 |
| Light-curve and Gaia readers | H | 20 % | S each | 16 |

What moves against the horizon memo's own order: the radial profile falls
from Phase 8's first wave to later, because it serves the smallest
persona, unless the paper needs it (horizon §11 Q7); the photometry reader
and filter library rise to Phase 8's first wave; three small docs and
results items that the horizon memo did not have rise to the top because
they serve everyone. The emulator keeps its place: it is the only route
to the flagship's model and to NUTS for 40 % of users.

## 7. The walkthrough plan

Each walkthrough is a script run on the beta, recording imports, lines
the user writes, wall time to the first posterior, the refusals met, and
the first point of hand work. Scripts live in the orchestrator's
scratchpad while the memo is drafted and move to `examples/walkthroughs/`
only if Peter wants them kept.

1. **A** — photometry of one object, a modified blackbody, emcee and
   dynesty, the table and three plots, three unit and wavelength probes.
   **Run; Appendix A.**
2. **B** — the quickstart's IRS spectrum with a two-component model, the
   GP likelihood, the diagnostics read as a user would. Next.
3. **C** — the carbon-star example under `SBIEngine` from a cold cache and
   from a warm one, the time to posterior both ways.
4. **E** — the population example on fifty synthetic objects, then the
   same fifty through the catalogue loop with the artefact cache on this
   machine, timed per object.
5. **D** — the contest OIFITS file through the reader and the binary fit,
   as the measure of what a reader saves.
6. **H** — the astrometry example as it stands.

G is reading, not running; F is §4 (1)'s guide, which the portable-model
item produces.

## 8. Questions for Peter, and his rulings of 2026-10-06

**Ruled**: (2) the ensemble engines start from the optimiser's mode by
default — W7.12; (3) the small items are Phase 7 fillers — W7.12, W7.13,
W7.14 drafted, W7.15 proposed; (4) the other walkthroughs are sketched and
run where the beta allows (Appendices B–H); (5) the radial profile is not
needed for the current paper but **is required for papers planned soon**,
so it keeps its place in Phase 8's first wave (horizon memo §3, §11 Q7).
On (1), the flagship's data, Peter's steer and the orchestrator's reading
of the paper are in §8.1; the choice is his.

### 8.1 The flagship's data: the paper read, and a recommendation

Peter's steer: analysing SAGE would be better without computing a new
grid, and could be a good example of amortised SBI. The paper he linked,
2026MNRAS.545f2221M (arXiv:2512.07573), is Marshall et al., "Systematic
determination of dust properties for a sample of 133 spatially resolved
debris discs" — Peter is a co-author — which fits simple analytical
radiative-transfer models of debris dust emission to multi-wavelength
photometry from the near-infrared to the millimetre, with the disc radius
from resolved imaging as an input, and derives per disc the minimum grain
size, the dust mass and the size-distribution exponent, finding a
population value q = 3.49 (+0.38, −0.33) and a trend of q with disc
radius. (Read from the arXiv abstract page; the data-availability
statement was not visible there.)

**Assessment.** It is a very good example, and a different one from
SAGE: it exercises A (one disc's SED with an analytic dust model — grain
size distribution, optical constants through Mie efficiencies, the
`miepython` route W6.13 (C2) already uses), E (the population of q, and
the q–radius trend as a hierarchical regression: a population whose
hyperprior location is a function of a per-object covariate, which is a
`Derived` on a buffer — W7.0's grammar over a per-object constant — under
a `HierarchicalPrior`, a construction the memo did not anticipate and
should be checked at W7.0), and the misspecification story (a smooth
analytic model against real SEDs with photospheric residuals and
silicate features, with the GP localising where). It does not need SBI:
the model is cheap. SAGE with a GRAMS-style grid is the amortised-SBI
example — the grid as the simulator, context amortisation over each
object's noise (W5.10), one trained posterior applied to thousands of
sources, reweighted as a population — and the catalogue-loop and store
work of Phase 9.

**Recommendation**: both, in sequence. The debris-disc sample as the
Phase 8 flagship (A + E + misspecification; the data are the
co-authors' and the model is public physics), staged as one disc, then the
133 as a population with the regression hyperprior, then the GP
localisation on the discs the analytic model fails; SAGE with GRAMS as
the Phase 9 flagship for amortised SBI at scale. Two questions remain for
Peter: whether the paper's photometry table and radii can be committed
under `examples/` or must be downloaded in the example, and whether the
analytic model is to be written fresh from the paper or ported from the
authors' code.

### 8.2 The questions as asked

1. **The flagship's data and model**: SAGE (public via VizieR/IRSA, with
   GRAMS grids public) or the 2026MNRAS.545f2221M sample (the orchestrator
   has not read it — which catalogue and model does it use, and is the
   sample public)? And is a GRAMS-style grid the model, or the carbon-star
   Hyperion setup emulated?
2. **The ensemble engines' default start** (Appendix A, finding 1): the
   optimiser's mode by default with `initial="prior"` as the escape, or
   the prior by default with the warm start documented on the quickstart?
   Recommendation: the mode by default — every walkthrough user will hit
   this first.
3. **The three small items at the top of §6** — the convergence verdict,
   the portable-model guide, the diagnostics page: as Phase 7 fillers now
   (recommended; each is S and serves everyone), or Phase 8?
4. **The walkthroughs**: run B–H now as this memo's appendix, or only A
   and the flagship's first stage?
5. **The radial profile's place** stands or falls on Q7 of the horizon
   memo: does the paper need it?

## Appendix A — persona A on the beta (2026-10-06)

Two scripts in the orchestrator's scratchpad (`walk/persona_a.py`,
`walk/persona_a2.py`), run with `pixi run --frozen -e dev` on this machine.
Synthetic photometry in nine bands (2MASS Ks, WISE W1–W4, MIPS 24 and 70,
PACS 100 and 160) from `ModifiedBlackBody(T = 180 K, β = 1.6, scale = 5)`
through `SyntheticPhotometry.from_library`, 8 % Gaussian noise; the fit
with `IndependentNoise`, three free parameters under log-uniform and
uniform priors.

**What the user writes**: 11 import lines, about 20 lines to the
`FittingProblem`, one line per engine, one per plot. Imports: `ampere`
0.5 s, the four subpackages 0.9 s in all.

| Step | Wall time | Result |
|---|---|---|
| `from_library` over nine filters | 0.7 s | — |
| `FittingProblem` | < 0.1 s | — |
| emcee, 24 walkers × 1 500 steps, burn-in 500, from the prior | 17 s | **R-hat 1.50, ESS 44** — not converged; T = 216 ± 245 K |
| `summary` | 0.1 s | the arviz table; no warning at R-hat 1.5 |
| `plot_corner` | 0.5 s | "too few points for contours" ×3 from the plotting library |
| `plot_posterior_predictive` after `add_posterior_predictive` | 8.4 s | slow for nine points |
| `plot_residuals` | — | refused: call `add_residuals` first (the message says so) |
| `to_netcdf` | 0.1 s | — |
| dynesty | 7 s | ln Z = −278.8 ± 0.4 |
| emcee with `GaussianProcessNoise` on nine points | 9 s | runs; R-hat 1.3; the GP amplitude and length scale unconstrained, as they should be |
| `optimise` (scipy route, eight starts) | 2.9 s | T = 180.4, β = 1.588, scale = 5.07 — the truth |
| emcee from the optimum, same budget | 19 s | **R-hat 1.06, ESS 300+**; T = 180.4 ± 0.9 K |
| emcee, 48 walkers × 6 000 steps from the prior | 197 s | R-hat 1.22 — still not converged |

**Findings, in order of what they cost a user.**

1. **The first fit a user would write does not converge, and nothing
   says so.** From prior draws over a log-uniform scale spanning six
   decades, the ensemble at a budget of 36 000 evaluations reaches R-hat
   1.5; ten times the budget reaches 1.22. The optimiser warm start
   (W6.7, `optimise` then `initial=`) reaches 1.06 at the first budget in
   22 s total — but it is opt-in and undocumented on the quickstart, so a
   user would not know. `summary` prints R-hat 1.5 without comment. This
   is Peter's "convergence is the biggest problem", reproduced on the
   simplest problem the package can be given. Two remedies, both cheap:
   the convergence verdict and emission-time warning (§5), and **the
   optimiser's mode as the ensemble engines' default start**, with the
   prior start kept behind `initial="prior"` — a behaviour change on
   `EmceeEngine` and `ZeusEngine` that needs Peter's ruling, since every
   recorded ensemble run so far started from the prior.
2. **The shipped `ModifiedBlackBody`'s `scale` is not a flux.** At
   `scale = 5` the model's fluxes are 10⁶–10¹⁵ Jy: `B_ν` carries its
   per-steradian magnitude, so `scale` is a solid angle in all but name
   and a user with catalogue fluxes in Jy needs a prior reaching 10⁻¹⁶.
   The docstring says "`scale` stays interpretable as the flux the source
   would have if it radiated as a pure blackbody at that wavelength",
   which reads as Jy. A persona-A user would spend an afternoon on this.
   Fix: a `scale` whose unit is the flux at the reference wavelength (the
   docstring's own promise), or a documented solid angle; S, with the
   examples and tests that pin the current value amended. The same check
   on `BlackBody` and `PowerLaw`.
3. **Photometric alignment is by wavelength, so a catalogue's own
   wavelengths are refused.** The step tabulates each filter on the
   user's grid and emits its own effective wavelength (2.16, 3.37, …
   160.93 µm); the library's pivot wavelengths differ in the second
   decimal (160.98) and `check_alignment` refuses the observed container
   by name. The only accepted route is to run the step once on a dummy
   prediction and copy its axis — hand work the quickstart hides because
   its data are simulated through the same step. A `PhotometricPoints`
   sample's identity is its filter name; alignment on that kind should
   key on the name and take the wavelength from the step. S in the core
   (a `results_schema.md` §8 amendment on the alignment rule for
   `PhotometricPoints`), or the reader item's first deliverable.
4. The refusals for a mJy flux and an Ångström axis are correct and name
   the fix (`.to_unit(...)`; "make the instrument chain produce the
   observed container's own coordinates") — the second message is the
   one finding 3 would make unnecessary. The hypothesis that a refused
   composition leaves a shared `Instrument` or model mutated is **refuted**:
   the log posterior at the truth is identical before and after.
5. The posterior-predictive plot's 8 s on nine points is the replicate
   draw over every stored draw; a default thinning would make it instant.
   The three "too few points to create valid contours" warnings on a
   corner plot of 24 000 draws are the plotting library's and should be
   silenced or explained.
6. The GP likelihood on nine photometric points runs and leaves its
   hyperparameters at the prior — correct behaviour, but a user who reads
   the concept page will try it, and a sentence on `photometry_spectra.rst`
   saying when the flexible likelihood has nothing to learn from would
   save the question.

**What the walkthrough did not test**: a real catalogue (the reader), a
user-written model, a notebook rather than a script.

## Appendix B — persona B on the beta (2026-10-06)

`walkthroughs/persona_b.py`: the real Spitzer IRS spectrum the quickstart
ships (`examples/test_data/cassis_yaaar_spcfw_14191360t.fits`, PG
1011-040, 360 unique samples over 5.2–37.4 µm after dropping the
overlapping orders' duplicates, median S/N 27), a `PowerLaw` continuum,
`IndependentNoise` then `GaussianProcessNoise(Matern32, QuasisepGP)`,
emcee from `optimise`'s mode, the diagnostics a user would reach for.

| Step | Wall time | Result |
|---|---|---|
| reading the CASSIS file by hand | 9 lines of astropy | the overlapping orders must be de-duplicated or `QuasisepGP` refuses the unsorted axis — the quickstart's `np.unique` trick, which a user must know |
| `optimise`, independent | 4 s | norm 0.00288, index 1.11 |
| emcee 24 × 1500 from it | 44 s | R-hat 1.03; a tight continuum |
| `optimise`, with the GP | 9 s | **norm pinned at its prior bound 0.001** (W6.7's Powell-at-the-bound finding, reproduced on real data) |
| emcee 24 × 1500 from that optimum | 73 s | **every walker sits on the bound**: norm sd 0, ESS 24 000, R-hat undefined; index 0.9 ± 1.0, unconstrained; the GP amplitude 0.033 Jy against a median flux 0.059 Jy — the GP is carrying the continuum |
| `add_residuals` + `plot_residuals` | 17 s | a correct warning that the whiteness family is scoped to standard-likelihood fits |
| `residual_whiteness` | 0.1 s | Q = 1578, p = 0.04 |
| `gp_localisation_score(run, problem)` | — | `TypeError`: takes the run alone — a signature the user guessed from `add_residuals(run, problem)` |
| `plot_posterior_predictive` | 39 s | slow: the replicate draw over every stored draw at 360 points |

**Findings.** (1) **The warm start can be worse than the prior start.**
When the optimiser lands on a bound, the ball around it is degenerate and
the ensemble never leaves: the fit is silently wrong with no diagnostic
except an ESS equal to the draw count. W7.12's default start must detect a
bound-saturated optimum (W7.6's check) and fall back to the prior with a
warning, and W7.6 should land before or with W7.12. (2) With the model at
a bound the GP absorbed the continuum entirely — the concept page's
"earthquake" paragraph in action, and the case the "Reading the
diagnostics" page (W7.14) must show. (3) The diagnostics' signatures are
not uniform (`add_residuals(run, problem)` against
`gp_localisation_score(run)`); a user guesses wrong once. (4) The
posterior-predictive plot's cost scales with draws × points; a default
thinning (W7.15 (3)) matters more here than on photometry.

## Appendix C — persona C, sketched from the rows (not run)

`examples/cstar` is the carbon star on Hyperion fitted with `SBIEngine`
(W6.13 (C2)); it needs the `hyperion` pixi environment. From its row: the
one-round quick budget left the envelope mass 0.14 dex short, the legacy
budget is about eleven hours, and the smoke row's forty simulations run
in the test. What the persona would meet, from the code rather than a
run: the wrapper is theirs to write (`Model` subclass, numpy path,
process pool, timeouts — the external-simulator example is the template);
`SBIEngine(cache=ArtefactStore)` makes a second run free; NUTS is refused
by name. The walkthrough to run when a slot allows: cold versus warm cache
at the quick budget, timed, and the calibration check. The emulator
(horizon §2 (b)) is what changes this persona's experience, not a
reader.

## Appendix D — persona D on the beta (2026-10-06)

`walkthroughs/persona_d.py`: the vendored contest OIFITS file through
`read_oifits`, V² and closure phases on one `sky` channel, the shipped
`Binary` on a 16 mas field, emcee from the optimiser.

| Step | Wall time | Result |
|---|---|---|
| `read_oifits` | 0.02 s | target, 600 V², 800 closure phases, 8 channels — one line |
| problem, 3 free parameters | 0.02 s | |
| `optimise` | 7 s | separation 3.886 mas, PA 126.4°, ratio 0.101 |
| emcee 24 × 800 | 39 s | R-hat 1.09; about 490 evaluations per second |

**Findings.** (1) The reader removes what was the whole cost of this
persona's entry: two containers in one call, the composition in a dozen
lines as the docs promise. (2) The log posterior at the *published*
geometry (5.0 mas at 30° east of north, ratio 1/8.9) is −456 346 while the
fit sits at 3.9 mas and 126° with sub-milliarcsecond precision. The
orchestrator did not resolve whether that is the `Binary` model's
position-angle convention, the four-pixel field it was built on, or the
walkthrough's own construction; W6.12's review verified the reader
against the published truth with an analytic model in numpy, not with
the shipped `Binary`. It is recorded as a question for the interferometry
page, not as a defect: a user porting a published geometry needs the
conventions stated where the model is documented, and a conformance row
fitting the shipped `Binary` to the contest file would settle it.

## Appendix E — persona E on the beta (2026-10-06)

`walkthroughs/persona_e.py`: the population page's example at twenty
objects on the numpy path (the page's own fit is NUTS on torch; the `dev`
environment has no torch), then the same twenty one by one as the loop a
user writes today.

| Step | Wall time | Result |
|---|---|---|
| the `Population` problem, 22 free | 0.02 s | |
| `optimise` on 22 dimensions | **79 s** | the scipy route's eight starts |
| emcee 48 × 3000 from it | 483 s | R-hat 1.23; μ = −1.49 ± 0.07 against −1.30, σ = 0.30 ± 0.05 against 0.35 |
| the loop, 20 single-object fits with `to_netcdf` | 43 s | 2.1 s per object; 0.7 MB per run file; provenance attrs present on reload |

**Findings.** (1) A population on the numpy path is the wrong tool and
nothing says so: eight minutes to R-hat 1.23 with a biased μ, where the
docs page fits fifty members by NUTS. The population page says it; a
`Population` problem handed to an ensemble engine could say it too (a
warning naming NUTS when `free_size` exceeds a threshold). (2) **The
optimiser's cost grows with dimension**: 79 s before the first draw on 22
parameters, which bounds what W7.12's default start may spend — one start
or a time cap for the default, the eight-start route on request. (3) The
loop a user writes today is fine at twenty objects and extrapolates to
six hours and seven gigabytes at ten thousand on one core: the
catalogue-loop example wants the cluster scripts and a per-object file
smaller than 0.7 MB for a one-parameter fit (the provenance and the
stored θ are most of it), and the store memo (horizon §4) is the ten
thousand case.

## Appendix H — persona H on the beta (2026-10-06)

`python -m examples.astrometry --walkers 24 --steps 400 --burn-in 100`:
6 s wall clock, the six orbit and proper-motion parameters printed with
their truths and the 95 % coverage flagged (one of six missed at this
budget, as a budget this small should). No friction; synthetic data by
construction, which is the persona's gap (the Gaia and light-curve
readers).

## Appendix F — the roll-up: what the walkthroughs change

Across A, B, D and E, the package's refusals were correct every time and
named the fix in four of five cases; nothing crashed; the imports cost a
second. What cost the user was never a missing feature: it was a default
(the prior start), a convention (the model's `scale`, the alignment
wavelength, a position angle), or a cost nobody capped (the optimiser's
starts, the replicate draw). The items that follow are all S, and all in
Phase 7's filler list: W7.12 with the bound-saturation fallback and a
budgeted default, W7.15 with the plotting defaults, W7.6 before W7.12. The
one thing no filler fixes is that a population fit needs the native path
and a `dev` user has no torch — which is a documentation and packaging
statement for `install.rst`: who needs which extra.

