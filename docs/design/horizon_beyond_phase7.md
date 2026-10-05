# Horizon beyond Phase 7 — what is blocked, what unblocks it, and a shape for Phases 8–10

Status: **assessment, not a plan.** Drafted 2026-10-05 by Fable at Peter's
request on the day Phase 7 was ruled, against the code at the beta
(`v1.0.0b1`) and the Phase 7 draft. Each section states what the code does
today in its own names, what a scientific use case needs that it does not
do, what it would take, and a first recommendation. Nothing here is
scheduled; §9 proposes a shape and §10 collects the questions. When Peter
rules, the entries fold into `DEVELOPMENT_PLAN.md` §5 the way
`horizon_notes.md` did on 2026-09-10, and this file stays as the reasoning.

## 0. The short answer

The beta can fit spectra, photometry, images, time series, visibilities and
closure phases, with a flexible likelihood on each, through eleven engines
on three backends, and read one file format. What it cannot yet do falls
into five groups, in the order a user would meet them:

1. **A model that calls foreign code cannot enter a native problem**, so
   any radiative-transfer code, Fortran routine or legacy numpy model is
   confined to the gradient-free engines and SBI. The emulator route the
   plan reserved as horizon (c) exists only as the hand-built PHOENIX
   example. (§2)
2. **Azimuthally symmetric sources have no home**: there is no radial
   profile kind, no Hankel transform, no deprojection, so a ring-and-gap
   disc, a dust shell or a galaxy's surface-brightness profile is fitted
   through a full image or not at all. The transformations contract named
   the slot at the freeze and nobody has filled it. (§3)
3. **Hierarchy stops at one level**, so a survey of objects each with its
   own repeated measurements, or a population whose members are themselves
   populations, cannot be declared; and population inference at the scale
   amortised SBI makes possible has no store. (§4)
4. **One reader.** Every other entry point is arrays the user prepared.
   (§5)
5. **Model comparison and approximate-inference correction are "roll your
   own"** although every hook was reserved and the per-draw densities are
   stored. (§6)

The recommendation in one line: Phase 8 takes the foreign function and the
radial profile, because they unblock the most science per item and neither
touches the parameter layer Phase 7 is changing; Phase 9 takes hierarchy
at depth and at scale, once Phase 7's path populations have been used in
anger; Phase 10 completes inference and closes with `1.0.0`. Readers are
fillers throughout.

## 1. What is true today

The inventory the rest of this memo reasons from, by the code's names.

- **Kinds**: `Spectrum`, `PhotometricPoints`, `TimeSeries`, `Image`,
  `Cube`, `VisibilitySet`, `ClosurePhases`; `register_kind` for a user's
  own; `Layout.GRID` refusals lifted in the core (W5.21).
- **Steps**: `CalibrationScale`, `Resample`, `LSFConvolution`,
  `SyntheticPhotometry`; `PSFConvolution`; `FourierSample`,
  `BandwidthSmearing`, `TimeSmearing`, `ClosurePhase`, `Amplitude`,
  `SquaredAmplitude`; `EpochSample`. Each on the three backends.
- **Models in the package**: `BlackBody`, `ModifiedBlackBody`, `PowerLaw`;
  `UniformDisc`, `GaussianSource`, `Binary` and their visibility twins;
  `ReflexOrbit`. In `examples/`, not the package: the PHOENIX emulator, the
  NGC6302 multi-temperature dust model with opacity tables, the star-plus-
  disc model, the carbon star on Hyperion, the astrometry generators.
  `from_astropy` wraps any `astropy.modeling` model as a black box and
  promotes six classes natively.
- **Likelihoods**: seven families (`Gaussian`, `StudentT`, `Cauchy`,
  `ComplexGaussian`, `Poisson`, `Rice`, `VonMises`); `IndependentNoise`,
  `FractionalModelNoise`, `GaussianProcessNoise` with the latent-GP opt-in;
  kernels `Matern12/32/52`, `SquaredExponential`, `SHO`, `RotationTerm`,
  `Sum`, `Product`, `SpectralMixture`, `WarpedKernel`; solvers `DenseGP`,
  `QuasisepGP`, `HilbertSpaceGP`, EFGP, Vecchia; joint noise over channels
  (`RotationCoupling`, the heteroscedastic routes); `Censoring`; the
  shrinkage helpers.
- **Hierarchy**: `Plate`, `HierarchicalPrior`, `Population` with the plate
  and flat layouts, `DatasetCollection.plate`, population reweighting of
  archived fits; **one plate per parameter**. Phase 7 adds `Derived` and
  populations over a dataset's own parameters.
- **Engines**: emcee, zeus, dynesty, nautilus, ultranest; NUTS (pyro,
  numpyro), VI (mean-field, full-rank, Laplace, flow), blackjax (MCLMC,
  Pathfinder); `SBIEngine` (NPE, NLE, NRE, TMNRE) with context
  amortisation and the native batched path; the optimisers (scipy, native
  MAP, the empirical-Bayes GP warm start). Every run a `DataTree` with
  per-draw `log_prior`/`log_likelihood`, the evidence triple where one
  exists, the proposal density for every approximate engine.
- **Results**: six plots, LOO through the pointwise group, SBC and
  coverage, the artefact store, training sets, provenance schema 9.
- **Readers**: `read_oifits`. Phase 7 adds `read_jwst`.
- **The capability ladder** (`architecture.md` §1): rung 0 black-box
  (gradient-free engines and SBI), rung 1 reference, rung 2 native. A
  native problem is composed entirely from one backend's pieces; a rung-0
  model cannot appear in a rung-2 problem.

## 2. The foreign function, and the emulator

**The use case.** A user has a model that is, or calls, code ampere cannot
trace: a radiative-transfer code (Hyperion, RADMC-3D, Dusty, MCFOST), a
Fortran opacity routine, a stellar-atmosphere interpolator in C, a legacy
numpy model. Today that model runs on the numpy path only. It gets emcee,
zeus, the nested samplers and SBI; it does not get NUTS, VI, the GPU, or
even the jax ensemble of a problem whose *other* pieces are native. The
docs' "Very slow models" section sends such a user to SBI, which is right
for the genuinely expensive case and wrong for the cheap foreign function
with few parameters, where the user merely wants their routine inside an
otherwise native problem.

**Route (a): a `ForeignModel` on each native backend.** jax's
`pure_callback` and torch's `autograd.Function` both admit a host callable
inside a traced computation, with the output's shape and dtype declared
up front. A `ForeignModel(fn, output=...)` on the jax backend wraps the
callable and declares `DIFFERENTIABLE = False`, so the capability flags do
what the ladder promises: NUTS refuses by name, while every gradient-free
engine on the native path, the batched path (`vmap` over the callback's
batch rule, which `pure_callback` supports with `vmap_method`), and the
GP likelihood on the accelerator all work. With an optional Jacobian —
user-supplied, or finite differences declared explicitly as
`jacobian="finite-difference"` with its step — a `custom_jvp` makes the
model differentiable at the cost of one call per parameter per gradient,
which for a five-parameter Fortran routine is a fine trade and for a
radiative-transfer code is not; the flag says which. Provenance records
the callable by qualified name and a source hash, as the external-simulator
example already does for pickling. Cost: S on each backend plus the
conformance rows (a foreign numpy blackbody agreeing with the native one
on every engine); a `lowering.md` §8 amendment stating what a foreign
piece is allowed to be. This is the cheap item that answers "can I use my
routine?" with yes.

**Route (b): the emulator as a first-class model** — horizon (c) of the
plan, reserved since 2026-09-01 and still unbuilt except by hand. The
PHOENIX emulator (W6.13 (4)) is the prototype: a model trained on
`(θ, ModelResult)` pairs that `simulate_many` produces, present on all
three backends, differentiable where the original was not, with its
training set's spec hash in provenance so a change to the simulator
invalidates it. What a package-level `EmulatedModel` adds over the
prototype: `train_emulator(model, prior, budget, architecture=...)` over
the existing training-set writer and artefact store; coordinate-conditioned
output (a neural-operator or SetConv decoder) so the emulator answers on
any sampling of its output function and can take part in requirements
negotiation, which the prototype's fixed grid cannot; and the one idea
that belongs to this project specifically — **the emulator's own error
enters the likelihood as a noise component.** An emulator is a
misspecified model by construction, and its predictive variance (an
ensemble's spread, a GP emulator's posterior variance, a held-out residual
model) is exactly a `FractionalModelNoise` term or a GP noise component
with a fitted amplitude, so the flexible likelihood absorbs emulation error
the way it absorbs physics error, and the M2 methodology measures whether
it did. Bayesian optimisation returns here as acquisition for the training
set, as the inference memo's §10 ruled. Cost: a design memo first (the
emulator contract: training set → model, the error model, the
coordinate-conditioned decoder, what provenance carries), then an L item.
This is the valuable item, and the one that turns every radiative-transfer
code into a NUTS-capable model.

**Recommendation.** Both, (a) before (b): (a) in Phase 8's first wave
(S, Opus, disjoint from Phase 7's parameter layer), the memo for (b) as
Fable's own work in Phase 8, the item itself Phase 8's long one. A worked
example for each: a Fortran-style opacity routine inside a native dust
model for (a); Hyperion's carbon star emulated and fitted with NUTS for
(b), beside the existing SBI fit of the same model so the two routes are
compared on one problem.

## 3. Radial profiles and azimuthal symmetry

**The use case**, in Peter's words, "radial profiles of observables (flux,
visibilities, etc.) assuming azimuthal symmetry". Two readings, both real:

- **The model is a radial profile.** A ring-and-gap protoplanetary disc
  observed by ALMA, a limb-darkened or shell-surrounded star observed by
  VLTI, a galaxy's Sérsic profile, the projected dust shell of an AGB star
  (which a 1D radiative-transfer code produces as `I(r)` directly). The
  observable is visibilities, an image, or photometry, and the symmetry
  makes the forward model one-dimensional: `V(ρ) = 2π ∫ I(r) J₀(2πρr) r dr`,
  a Hankel transform of order zero, with inclination, position angle and
  offset handled by deprojecting the `(u, v)` coordinates and a phase ramp
  before it.
- **The observable is a radial profile.** A surface-brightness profile
  extracted from an image by azimuthal averaging, or a visibility profile
  binned in deprojected baseline length — the form in which radio and
  optical interferometry data are routinely published and fitted
  (`frank`, Jennings et al. 2020, fits a non-parametric profile to
  deprojected visibilities with a GP prior on the profile, which is this
  project's own idiom).

**What exists.** The transformations contract's §10 table names "Hankel
transform: radial profile → visibility amplitudes, for radial profiles and
for visibilities as functions of u–v distance alone" as a standard-library
slot since the freeze, unfilled. `VisibilitySet` carries `(u, v)` per
sample; `FourierSample` samples an `Image` model on them; the interferometry
model trio has closed-form visibilities. There is no radial kind, no
Hankel step, no deprojection step, no azimuthal average. jax and torch
both ship `J₀` natively (`jax.scipy.special.bessel_jn`,
`torch.special.bessel_j0`), so the transform lowers without a foreign
function.

**What it takes.** A modality item by the template: a `RadialProfile`
kind (coordinate `r` in angular units, value a surface brightness, the
point layout); steps `Deproject` (`VisibilitySet` → `VisibilitySet`:
inclination, position angle, offset as parameters, the coordinates rotated
and scaled and the phase ramp applied), `HankelTransform` (`RadialProfile`
→ `VisibilitySet` at the deprojected baseline lengths, a dense `N_ρ × N_r`
quadrature at the sizes that matter, the discrete Hankel transform as the
large-N route), `Render` (`RadialProfile` → `Image` on a pixel grid, for
the image and photometry routes) and `AzimuthalAverage` (`Image` →
`RadialProfile`, the observed-side step for the second reading); models
`GaussianRing`, `PowerLawDisc` with gaps, `Sersic`, and a non-parametric
`LatentProfile` — the profile as a latent GP in `log r`, which the
`CONSUMES_LATENT_GP` opt-in already supports and which makes the
frank-style fit a declaration; the conformance rows (the uniform disc's
Hankel transform against `UniformDiscVisibilities`, the Gaussian's against
`GaussianSourceVisibilities`, a `Render`-then-`FourierSample` route against
the direct `HankelTransform`); and a worked example on real data (a
public ALMA continuum uv table — the DSHARP sets are public, and `frank`
ships a small one — or a VLTI set through the OIFITS reader). The
flexible likelihood then does what it is for: a parametric disc plus a GP
correction *in r* localises where the symmetric model fails, which is the
misspecification story on a new axis.

**Recommendation.** Phase 8, L, Opus, after the `Derived` node lands —
`Deproject`'s inclination and the ring parameters are natural `Derived`
customers (a ring's width from its radius and an aspect ratio) and the
example should use them. It unblocks four use cases at once (discs,
shells, galaxies, stellar diameters with limb darkening) and the
interferometry stack is the one most of the project's own data goes
through.

## 4. Hierarchy at depth and at scale

**Nested populations.** A parameter carries one plate. A survey of star
clusters each with member stars, a set of galaxies each with several
spectra, a population of discs each with multi-epoch visibilities — each
needs a population whose members are populations. The parameters contract
kept recursive merge open and then ratified the nested, lossless design
(one merge per level, every mapping retained, `Binding.index` for
per-element routing), so the *routing* substrate nests already; what does
not is the declaration (`Plate` is one-level), the plate layout (one
array-valued parameter per plate), the native resolve functions (one
`numpyro.plate`/torch loop level) and the emitted coordinates (one plate
dimension). Phase 7's W7.1 will show how far "a population over a path"
stretches, since a path into a dataset *is* a second level in disguise.
The memo-then-item pattern of W6.11 is right: a Fable memo after W7.1 has
been used, naming the tree layout, the per-level non-centring (the `Derived`
node generalises), the coordinates, and the two backends' nesting; then an
L Opus item with its conformance rows. **Phase 9**, not before — the
customer that forces it is most likely the IFU cube (§7) or the
population-at-scale work below, and designing it before either exists
would be designing in the abstract.

**Hierarchical SBI** (horizon (d)): plate-aware groups exist; what is
missing is a `simulate()` that draws hierarchically (hyperparameters, then
members) and an embedding that respects the plate (a set of sets). M, after
nested populations, since a one-level hierarchical SBI is already
expressible through `Population` and the set embedding.

**Population inference at scale** (horizon (b)'s constraint, Peter,
2026-09-10): reweighting archived fits works today one file per source;
at 10⁴–10⁹ sources from an amortised posterior it needs a columnar,
partitioned store (Parquet or Zarr — the per-draw columns `log_prior`,
`log_likelihood`, the proposal density are already separable from
provenance by rule) and a reweighting that streams. L, with a design memo
that chooses the store; the customer is a survey-scale SBI run, which the
cluster now makes possible. Phase 9.

## 5. Readers

The reader pattern (W6.12: a front door per observable, `astropy.io.fits`
only, a vendored public file with its provenance, a refusal by name for
every ambiguity) makes each reader an S item. In order of how many users
each unblocks:

1. **Photometry from the catalogues** — the SED-fitting entry point most
   users start from: `astroquery`'s VizieR photometry (the SED service) or
   a user's table of fluxes and filter names into `PhotometricPoints`,
   with the filter names resolved against a filter library. This is where
   issue #64 lands: a v2 `FilterLibrary` fetched from the SVO Filter
   Profile Service with a local cache and a pinned snapshot for the tests,
   replacing the bundled pyphot library that `SyntheticPhotometry.from_library`
   still reaches through legacy `pyphot_compat` (W6.1's carried finding).
   M, because the library is the larger half.
2. **Generic 1D spectra** through one optional adapter: `from_specutils`
   taking a `Spectrum1D` gives every format specutils reads (SDSS, ESO
   SDP, HST STIS/COS, generic FITS tables, ASCII) for one optional
   dependency, without ampere owning a format zoo. S.
3. **Spitzer IRS (CASSIS)** — Peter's own data path (the legacy star-disc
   example reads one); S, and it closes the missing-file question from
   W6.13 (C1).
4. **Radio and millimetre visibilities**: a uv-table reader (the
   `.npz`/text form `frank` and `galario` consume) and UVFITS through
   astropy into `VisibilitySet`; CASA measurement sets are out (the
   dependency is heavy and the conversion is a one-liner in CASA). S;
   pairs with §3's example.
5. **Light curves** (`lightkurve`'s `LightCurve` or a TESS/Kepler FITS
   table) into `TimeSeries`: S; it gives the time-series kind a real-data
   entry it lacks.
6. **Gaia epoch astrometry**: DR4's individual-epoch data (expected around
   the end of 2026) into the astrometry modality's `TimeSeries` pair, with
   the Hipparcos IAD format as the available stand-in today. S, timed to
   DR4.
7. **JWST cubes** (`s3d`) into `Cube` — with the IFU item (§7), not
   before.

X-ray PHA/RMF/ARF is a modality (the awkward-instrument sketch), not a
reader, and stays on the horizon.

## 6. Inference, model comparison and correction

The inference memo's ranking was executed to tier 1 at W5.14 and tier 2's
optimisers at W6.7. What remains, by value:

1. **Model comparison as a module** (horizon (e), the hooks all in place):
   `ampere.results.compare` — the learnt harmonic-mean estimator
   (`harmonic`) on archived draws, bridge sampling between two runs,
   Savage–Dickey for nested models, a Bayes-factor table across runs
   sharing a data hash, Bayesian model averaging of posterior predictives.
   The docs page says "roll your own" today. M; dependency-light.
2. **Approximate-inference correction** (horizon (f)): `ampere.results.correct`
   — Pareto-smoothed importance reweighting of VI, Laplace, Pathfinder
   and SBI runs toward the true posterior, using the stored proposal
   density and the problem's `log_prob`, with the `k̂` diagnostic stored;
   the same code is the population reweighting's core. M.
3. **pocoMC** (SMC with flow preconditioning, batched evaluations): the
   engine the native batched path was waiting for; multimodal and
   correlated posteriors with fewer calls, an evidence as a by-product.
   S–M behind an extra.
4. **jaxns**: a nested sampler that runs *on* the jax problem, so evidence
   on the GPU and under `vmap`; S behind an extra, once W7.10 makes the
   jax solver accelerator-capable.
5. **A noise-model selection study**: the question every user asks
   ("which kernel?") answered by evidence and LOO comparison across the
   kernel library on the M2 scenarios, as docs guidance with a measured
   table. M; depends on (1).
6. Lower: `snowline`, `eryn` parallel tempering, `nutpie`, ensemble Kalman
   inversion — on a use case only, as the memo ruled.

## 7. Modalities still on paper

- **IFU cubes.** The `Cube` kind, `PSFConvolution` on `Image` and the
  reduced-rank solvers exist; a correlated noise model on a `Cube` needs
  a separable spatial ⊗ spectral kernel whose Kronecker structure keeps
  the solve tractable (two small factorisations instead of one of size
  `N_x N_y N_λ`), the per-spaxel plate through `DatasetCollection.plate`
  as the alternative decomposition, the `s3d` reader, and an example on
  a JWST or MUSE cube fitting a line map with a misspecified continuum.
  L; the first real customer for nested populations (spaxels within
  cubes within a sample). Phase 9.
- **Line fluxes** are Phase 7's W7.9 (ruled in).
- **X-ray** (the awkward-instrument sketch): the response-matrix step,
  exposure, background as a second dataset, Poisson family — all
  expressible since X-1 was fixed at the freeze; an item when a user
  asks.
- **Polarimetry** as the vector-observable customer for the joint noise
  model (W5.9) — the Stokes kind and the rotation-structured coupling are
  a S–M item with no new contract.

## 8. The package, the docs and the release train

- **The model zoo.** The v2 package ships three spectral models; the
  examples hold the ones users want (multi-temperature dust with opacity
  tables, the PHOENIX emulator, star plus disc). Promoting them into an
  `ampere.models` namespace with their native twins, tests and a gallery
  page is M and mostly mechanical; the extinction feature joins them.
  Issues #65, #21–#23 (RT codes, Cloudy, LIME) stay closed: with §2's
  foreign function and emulator, an adapter is user-side.
- **A declarative problem file.** `to_spec` exists; a `from_spec` round
  trip plus `ampere fit problem.toml` would let a non-programmer run a
  fit from a file. S–M; useful for the paper's reproducibility claims,
  low priority against science.
- **Docs.** `advanced.rst` still ends sections with "more detailed
  tutorials will be available soon"; a cookbook of short recipes (one
  page per question a user asks) and a gallery of the examples' figures
  are the Phase 8–9 docs items; the accessibility pass is W7.7.
- **Python 3.15** joins the CI matrix when it ships (this month); 3.12
  stays the floor until arviz moves.
- **conda-forge**: a feedstock for `ampere-astro` is the install route
  most astronomers reach for first; every dependency is on conda-forge
  already. S, after `1.0.0b2`, Peter's maintainership.
- **The jax quasiseparable solver**: W7.10 lifts the accelerator refusal;
  `BATCHABLE` for the celerite2 provider stays false until a `vmap` rule
  exists, and the scan provider gives batching for free — note for W7.10's
  design.
- **celerite2 as a single point of failure** (plan §6): the scan provider
  of W7.10 is also the first in-house implementation of the recursion; if
  it holds on the CPU within a factor of a few, the dependency becomes
  optional and the risk closes.
- **`1.0.0` final**: the criteria to rule — the known-limitations list
  down to what is deliberate, conda-forge live, the paper's examples
  pinned to it, v2's deprecations removed, and the Zenodo record. The
  phase that closes with it is Phase 10 in §9's shape.

## 9. A shape for the next three phases

- **Phase 8 — any model, on any backend** (after Phase 7's beta):
  `ForeignModel` on torch and jax (S); the emulator memo (Fable) and
  `EmulatedModel` (L); the radial-profile modality with the disc example
  (L); the filter library and the photometry reader (M); the model zoo
  promotion (M); fillers: `from_specutils`, the IRS reader, the uv-table
  reader, `compare` (M). Theme: a user brings their own model and their
  own data and gets every engine.
- **Phase 9 — hierarchy at depth and at scale**: the nested-populations
  memo and item (L); hierarchical SBI (M); the columnar store and
  population inference at scale (L); the IFU cube with the `s3d` reader
  (L); fillers: the light-curve and Gaia readers, `correct` (M),
  polarimetry (S–M).
- **Phase 10 — inference completion and `1.0.0`**: pocoMC and jaxns (S
  each); the noise-model selection study (M); the declarative problem
  file (S–M); the cookbook and gallery (M); conda-forge; the final
  release with its criteria met.

Each phase keeps the Phase 5–7 working agreement: two agents at a time,
CI on the push as the gate, a memo before any §4 change, a decision list
for Peter at drafting.

## 10. Questions for Peter

1. **Radial profiles**: both readings (the model as a profile with the
   Hankel route, and the profile as an observable with `AzimuthalAverage`),
   or the model side first? Recommendation: both in one item; the second
   reading is two steps once the kind exists.
2. **The foreign function**: `ForeignModel` first and the emulator after
   its memo, as §2 recommends, or the emulator alone on the argument that
   it is the route that gives NUTS? Recommendation: both, (a) first.
3. **Nested populations' trigger**: Phase 9 by plan, or only when a
   customer appears? Recommendation: Phase 9 by plan, with the IFU cube as
   the customer written into the same phase.
4. **Readers' priority**: §5's order, or Peter's own data first (IRS,
   VLTI, ALMA)? Recommendation: §5's order with IRS as the first filler.
5. **The model zoo's home**: `ampere.models` in the package, or the
   examples as they are with a gallery page? Recommendation: the package,
   because a model in the package has native twins, conformance rows and
   provenance, and one in an example has none of those guarantees.
6. **`1.0.0`'s criteria**: §8's list, amended by Peter; and whether a
   software paper (JOSS or a journal) accompanies it, which fixes the
   docs items' priority.
7. **The paper**: which of the above the paper's revision needs — the
   radial profile and the emulator are the two with figures a paper would
   want — so the order can serve it.
