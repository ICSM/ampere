# Phase 4 decisions D1 and D2 — where a new observable lives, and what a closure phase is indexed by

Status: **discussion memo for Peter, 2026-09-11.** Written by Fable from three
read-only surveys of the tree (the v2 kinds and chain, the user's experience
of building a fit, and the legacy data classes), at Peter's request that the
placement question be settled from how the existing classes for spectra and
photometry actually behave rather than from taste. Nothing here is decided;
§6 lists the rulings asked for. Every claim about the code cites the file and
line it was checked against at `59454ea`.

The two questions from `WORK_ITEMS.md`'s Phase 4 header:

- **D1** — where a shipped observable's kind, steps and twins live, and what
  the namespace is called (Peter: not `modalities`; a shorter, simpler word,
  and decided only after looking at the existing interfaces).
- **D2** — what `ClosurePhases` is indexed by, and what the choice means for
  the dimensionality of the problem, bearing in mind that interferometric
  data are in general a function of wavelength as well.

D3 (W4.9 runs) and D4 (models and ordering; the cross-model review policy of
2026-09-05 extends to Phase 4) were ruled the same day and are recorded in
the items file.

---

## 1. How the existing interfaces behave

### 1.1 Legacy: the data object owns everything, and there is no instrument

The frozen API (`ampere/data`, `ampere/models`, `ampere/infer`) has no
instrument object of any kind. Instrument behaviour, the noise model, the
nuisance parameters *and their priors*, and the likelihood itself all live
on the data class; the model owns only its prediction and its own
wavelength grid.

| Concern | Legacy owner | Where |
|---|---|---|
| Likelihood | the data object | `data.lnlike(theta_slice, model_result)`, summed by `BaseSearch.lnlike` (`ampere/infer/basesearch.py:64-84`) |
| Filter integration (photometry) | the data object | pyphot `Filter`s on `self.filters`, the band integral inlined in `Photometry.lnlike` (`ampere/data/photometry.py:474-478`) and copied verbatim into `simulate` (`:586-603`) |
| Resampling (spectra) | the data object | a function-valued `self.resampler` (`spectres`, `np.interp` or a user callable; `ampere/data/spectrum.py:275-333`); no LSF or resolving-power model exists anywhere in legacy |
| Correlated noise | the data object | one hard-coded squared-exponential kernel, a dense N×N matrix inverted with `numpy.linalg.inv` (`spectrum.py:389-451`, `:572`) |
| Nuisance parameters and priors | the data object | `Spectrum.npars = 3` (`calVar`, cov weight, cov length), priors built in `__init__`; `Photometry.npars = 0` |
| The wavelength grid | the model | `self.wavelength` fixed at construction; spectra resample the model onto themselves, photometry is pushed onto the model's grid by the user calling `photometry.reloadFilters(model.wavelength)` out of band (`photometry.py:402-425`) |
| Negotiation | nobody | no mechanism tells the model what coverage or resolution the data need |

Three legacy facts bear directly on Phase 4:

- **The other observables were declared and never built.** `Image`,
  `Interferometry` and `Cube` exist in `ampere/data/data.py:200-350` as
  copy-paste stubs whose `__init__` lacks `self` (they cannot be
  instantiated) and whose `lnlike` raises before its dead placeholder
  arithmetic. `ModelResults.__init__` accepts `image`, `cube`,
  `visibilities` and `closure_phase` keywords and silently discards them
  (`ampere/models/results.py:14-18`). Legacy offers nothing to inherit for
  interferometry; the slot names are the only thing that survives.
- **The one attempt at a separate instrument object was on the abandoned
  `jax` branch**, harvested to `docs/design/harvest/jax/`: an
  `InstrumentModel` holding a list of transformations, a
  `DataGeneratingProcess` with `Scaling`, `OneDInterpolatingResampler` and
  `LineSpreadFunction`, and `Photometry`/`Spectroscopy`/`Interferometry`/
  `Imaging` as sibling `Data` subclasses. It never imported. Its structure
  is what v2 built properly.
- **The ordering trap the design exists to end.** A `Photometry` object on
  which `reloadFilters` has not been called raises inside `lnlike`, and
  `examples/star_disc.py` has that call commented out on both lines where
  it appears. v2's negotiation is the answer to this, which is why the
  placement of steps matters: a step's published requirements are what
  replaces the manual call.

### 1.2 v2: five nouns, grouped by layer and never by observable

The v2 tree groups everything by *layer* (neutral core; one package per
backend) and nothing by *observable*. For spectra and photometry the pieces
live in four places, and the pattern is the same for both:

| Piece | Where | Mechanism |
|---|---|---|
| The kind (`Spectrum`, `PhotometricPoints`, and also `TimeSeries`, `Image`, `Cube`, `VisibilitySet`) | `ampere/core/results_schema.py:1126-1553` | three `ClassVar`s (`AXES`, `LAYOUT`, `ALLOW_COMPLEX`) on a `FunctionSamples` subclass; recognised by `isinstance`, never by a registry |
| The reference-path steps (`CalibrationScale`, `Resample`, `LSFConvolution`, `SyntheticPhotometry`) | `ampere/backends/reference/instrument.py` | core ships **no** concrete step: "core is the vocabulary, and this is the library" (`instrument.py:5`; `transformations.md` §10) |
| The native twins | `ampere/backends/torch/instrument.py`, `ampere/backends/jax/instrument.py` | **no registry**: the same class names in parallel modules, a duck-typed `apply_flux` surface the realised problem checks (`backends/torch/problem.py:234`, `backends/jax/problem.py:196`), and the `BACKEND` flag that `foreign_parts` checks structurally (`core/dataset.py:409`) |
| The likelihood family (and the noise models, solvers) | `ampere/core/likelihood.py`, with family twins in `backends/torch/_families.py` and `backends/jax/families.py` | `register_family` keyed on `NAME` (`likelihood.py:1982-2019`) |
| Example models | `ampere/backends/{reference,torch,jax}/models.py` | one name per backend |

Two mechanisms that *look* like registries are not what a new observable
would register with: `core/lowering.py` has exactly two slots (prior
families and bijections; `lowering.py:20-42`), and `core/realisation.py`
registers whole problems per backend. `ampere.results.register_kind`
(`results/serialisation.py:151-196`) exists, but only so a kind can be
rebuilt by name from a stored run; an out-of-tree kind must call it or its
records cannot be read back. **W4.3's phrase "lowered through the registry
like every standard step" is therefore wrong**: there is no step registry,
and the twins associate by the three conventions above (§4 corrects it).

The out-of-tree extensibility claim is already *proven*, not promised:
`tests/core/thirdparty_polarimeter.py` defines a kind, a kind-preserving
step with a nuisance parameter, a kind-changing many-to-one step and an
instrument factory using only `ampere.core`'s public surface, and
`tests/core/test_transform.py::TestOutOfTreeExtension` runs it. Phase 4's
claim is a different one: that the *composition* holds for a kind-changing,
coordinate-changing step, complex data, wrapped angles and two datasets on
one channel. Where the shipped code lives does not change either claim.

One asymmetry between the backends is a live design fact for W4.3: the
torch twins **re-declare** everything (constructor validation, buffers,
published requirements) and are held to the reference by the conformance
battery (`backends/torch/instrument.py:10-31`), whereas the jax twins
**inherit** the reference class and override only the four flags and the
arithmetic (`backends/jax/instrument.py:8-31`), which is the pattern that
guarantees two backends publish identical requirements for one
declaration. The interferometry steps should say which pattern they follow.

### 1.3 How one builds the instruments for a fit today

The user-facing flow, taken from the only places it is written down
(`docs/source/overview.rst:129-163`, `examples/wstat_comparison.py:505-519`,
`tests/backends/test_reference.py:500-507`, `inference.md` §8):

1. **A model** is constructed with its own evaluation grid and a
   `channels=` keyword (default `"default"`), and returns one container per
   channel (`backends/reference/models.py:194-204`).
2. **A container** is built by hand from arrays with astropy units
   multiplied on at the call site. There is no reader, loader or bridge
   from `ampere.data` anywhere in v2; `docs/source/migrating.rst:30-34`
   says the two APIs do not interoperate.
3. **An instrument** is a list of steps: for a spectrum,
   `Instrument([LSFConvolution(resolving_power=R), Resample(observed.spectral_axis.values)], channel=..., label=...)`;
   for photometry, `Instrument([SyntheticPhotometry.from_library(filters, grid, detector=...)], channel="sed", label="phot")`.
   Kind-preserving steps must precede the kind-changing one, checked at
   construction (`core/transform.py:930-944`). A `Resample` target must be
   the observed container's own coordinates, never recomputed
   (`instrument.py:171-175`). A dataset with no instrument gets one that
   only binds.
4. **Two instruments on one channel need distinct labels**, because an
   instrument's label defaults to the channel name; the collision surfaces
   from `DatasetCollection`, whose comment names the photometry case as
   the motivating one (`core/dataset.py:1748-1761`).
5. **`Dataset(observed, instrument, likelihood=..., label=...)`**, a
   `DatasetCollection`, and `FittingProblem(model, datasets, seed=...)`,
   which binds, negotiates, compiles, merges and validates once in its
   constructor (`inference.md:556-575`). The user never calls
   `compile_for` or `configure_from`.

What the survey found missing, all of it relevant to the template page
W4.4 owes:

- **No shipped example composes photometry with a spectrum.** The step
  ships and is unit-tested; `docs/source/tutorials.rst:50-53` says the
  page "deserves a page of its own"; the only side-by-side appears in a
  `DatasetCollection` docstring with no instruments. The sketch
  `spectrum_photometry.md` uses a `SyntheticPhotometry` signature that
  differs from what shipped.
- **`configure_from` is already in use.** Its docstring
  (`core/transform.py:780-797`) says "the first real instances (the
  smearing steps) are Phase 4's", but `LSFConvolution` overrides it on all
  three backends (`backends/reference/instrument.py:342`,
  `backends/torch/instrument.py:480`, `backends/jax/instrument.py:243`) to
  learn its output range from the resampler after it. W4.1's "the first two
  `configure_from` instances" is wrong; the smearing steps are the first
  *cross-kind* uses, asking a kind-changing predecessor for extra samples.
- **No ingestion story at all.** Every example generates its data in
  process. Interferometry is the first observable with a standard archive
  format (OIFITS), which makes the question of where a reader would live
  part of D1.

---

## 2. D1 — where a shipped observable lives

### 2.1 What the question actually is

A shipped observable is four things: a **kind** (backend-neutral, three
class attributes), the **reference-path steps** (numpy), the **native
twins** (torch, jax), and optionally a **family** with its own twins and
some **example models**. The existing convention places each by layer.
The interferometry sketch (`interferometry.md:81-85`) was the first text to
propose grouping by observable instead: "this kind belongs in Phase 4's
interferometry package, not in `ampere.core`". The drafted D1 followed the
sketch and invented `ampere/modalities/` to hold the neutral kind and the
reference-path steps, with twins staying under `backends/`.

That draft has a defect the survey exposes: it moves the *reference*
backend's steps out of `backends/reference/`, so
`ampere.modalities.interferometry.FourierSample` and
`ampere.backends.torch.interferometry.FourierSample` are no longer
parallel names, breaking the rule the plan states as "one name per backend
everywhere; a native problem is composed entirely from one backend's
pieces". Whatever is decided, the reference steps belong in
`backends/reference/` beside the other three backends' twins.

### 2.2 The options

**A. Follow the layer convention exactly.** `ClosurePhases` joins
`core/results_schema.py` beside `VisibilitySet`; the steps go in a new
module per backend, `backends/{reference,torch,jax}/interferometry.py`; the
von Mises `sample()` and the circular GP closed form go where the families
already are. No new namespace, no naming question. The template for the
next observable is one sentence: "a kind in core, a module per backend".
The cost is honesty about the plan's "no change to `ampere.core`": one
kind is added (some forty lines, not a contract change — `results_schema.md`
§13 defines a kind as three attributes), and the module `results_schema.py`
grows by one class per observable.

**B. A grouping namespace by observable (the draft, corrected).**
`ampere/<word>/interferometry/` holds the kind and any neutral helpers;
steps and twins stay per backend under `backends/*/interferometry.py`.
What the neutral package then holds in Phase 4 is *one class*, plus
whatever conveniences are invented to justify it. The namespace is thin
until ingestion or observable-specific helpers arrive.

**C. Top-level packages named for the observable, no grouping word.**
`ampere/interferometry/` (later `ampere/astrometry/`, `ampere/ifu/`), the
way astropy has `astropy.timeseries` and `astropy.coordinates` rather than
a `modalities` layer. Same content as B's package; no word to choose. The
top level today is `core`, `backends`, `inference`, `results` plus the four
legacy packages, so one package per observable sits naturally beside them
and is the shortest possible import.

**D. A user-facing home layered over A.** Kinds in core and steps per
backend as in A, *and* a top-level `ampere.interferometry` as the
observable's front door: re-exports of the kind and the reference steps,
the analytic model trio, a reader (OIFITS), and the page. This is C's
package with A's placement of the contract-bound pieces; it is the only
option in which "where do I import interferometry from" has a one-word
answer while the backend rule still holds.

### 2.3 What separates them

- **The backend rule** (one name per backend; a native problem from one
  backend's pieces) is kept by A, C-as-stated-here and D, and was broken by
  the original draft. It should be treated as fixed.
- **Whether a kind may live in core.** The sketch said no on the analogy
  with steps, but the analogy does not hold: `transformations.md` §10 keeps
  *steps* out of core because they are implementations that differ per
  backend; a *kind* is declarative, backend-neutral, and the six existing
  ones, `VisibilitySet` and `Image` included, are in core already. The only
  argument for keeping `ClosurePhases` out is to make Phase 4 count as
  "out of tree", and §1.2 shows that claim is already carried by the
  polarimeter test.
- **The user's import experience.** Today a photometry user imports the
  kind from `ampere.core` and the step from `ampere.backends.reference`;
  A keeps that, which is consistent but not friendly. D is the option that
  improves it, and is the natural home for the reader when one exists.
- **The name.** Only B needs a grouping word. If B is nonetheless wanted,
  the candidates that are short, plain and not already taken by legacy:
  `obs` (observables, observations; `ampere.obs.interferometry`),
  `astro` (the sherpa convention, `sherpa.astro.data`), `kinds` (accurate
  for the kind, wrong for steps and readers), `probes`, `domains`. `data`
  is the right long-term word and is taken by legacy until Phase 6 retires
  it; `instruments` is wrong because a kind is not an instrument.

### 2.4 Recommendation

**A now, with D reserved.** Put `ClosurePhases` in core beside
`VisibilitySet`, the steps in `backends/*/interferometry.py` following the
jax inheritance pattern where torch's base class permits it, and record in
the plan that the "no change to core" sentence means "no contract change;
one kind added". Do not create a grouping namespace in Phase 4: it would
hold one class. When the first reader lands (OIFITS is the obvious
candidate; Phase 6's ingestion work is where readers belong), create
`ampere.interferometry` as the observable's front door per D, and let the
astrometry and IFU packages follow the same shape. That leaves no word to
choose, and it keeps the template for W4.9 and every later observable to
one sentence.

If Peter prefers a grouping namespace regardless, `obs` is the shortest
candidate that is not wrong, and the placement of steps under `backends/`
and the kind in core should stand either way; the namespace would then be
the front door of D under a grouping word, not a place the contract-bound
code moves to.

---

## 3. D2 — what a closure phase is indexed by, and the dimensionality

### 3.1 The object

A closure phase on the triangle of telescopes (i, j, k) at wavelength λ and
time t is the argument of the bispectrum,
Φ = arg(V_ij · V_jk · V_ki), with the three visibilities sampled at the
three baselines' spatial frequencies (u_ij, v_ij), (u_jk, v_jk),
(u_ki, v_ki), where u = B_x/λ, v = B_y/λ in wavelengths. The three
frequencies sum to zero, so two determine the third. Every one of those
frequencies is *already wavelength-divided*: in the frozen `VisibilitySet`
convention (`u`, `v` dimensionless, `Order.ANY`;
`core/results_schema.py:1491`) a spectrally dispersed observation of one
baseline is a set of points along a radial spoke through the origin, and
the wavelength is absorbed into the coordinate.

The two candidate signatures:

- **Four axes** `(u1, v1, u2, v2)`, dimensionless, `Layout.POINTS`, with a
  `triangle` label and (after W4.2 gives `extra_coords` units) a per-sample
  `frequency` in Hz in `extra_coords`. This is what W4.1 is drafted to.
- **No axes**, a `triangle` label in `extra_coords` in the manner of
  `PhotometricPoints.filters`, with the wavelength and time also as extra
  coordinates.

### 3.2 The fact the choice turns on

A GP kernel evaluates separation as a Euclidean distance over the
container's **declared axes**, requires them to share one unit, and never
sees `extra_coords` (`core/likelihood.py:720-745`; `_build_extra_coords`
at `results_schema.py:808-843` refuses even a `Quantity` there). Whatever
is put in the axes is what a noise model can correlate over; whatever is
put in `extra_coords` is visible to an instrument step and to the identity
check, and invisible to every solver. So the signature is not a labelling
choice. It decides whether `GaussianProcessNoise` can be composed on
closure phases at all.

### 3.3 Dimensionality, in the four senses the word can take here

**(a) Free parameters: none.** The (u, v) coverage is a buffer, never a
parameter, on both forms (`interferometry.md` §3).

**(b) Data count: none.** The container carries every measured closure
phase either way: N_cp = N_triangles × N_λ × N_t. What the count hides is
that only (N_tel − 1)(N_tel − 2)/2 of the C(N_tel, 3) triangles are
independent at each (λ, t) — 3 of 4 for a four-telescope array, 10 of 20
for six — and closure phases sharing a baseline are correlated. An
`IndependentNoise` likelihood ignores that under both forms. Only a joint
noise model can fix it (Phase 5, design horizon (h)); under the four-axis
form such a model can *discover* shared baselines from coordinate
equality, under the label form the baseline identity has to be carried
separately on both containers.

**(c) The coordinate dimension a noise model sees: 4 against 0.** This is
the decisive sense. Under the four-axis form a GP over closure phases
lives in a four-dimensional space whose coordinates all share one unit, so
an isotropic kernel is legal and `DenseGP` accepts it (the O(N)
`QuasisepGP` path refuses by name, as it does for `VisibilitySet`'s two
axes). Under the label form there are no axes, so no GP can be composed on
closure phases, and W4.4's flagship question — does the flexible
likelihood keep the binary calibrated when the sky model is wrong — could
be asked of the visibilities only.

The four-dimensional GP is also the *physically right* space for the
purpose the flagship serves. A missing or wrong sky component adds a
smooth function δI(x, y) to the image; its effect on the bispectrum is
smooth in the three (u, v) arguments, so a GP over (u1, v1, u2, v2)
describes sky misspecification on closure phases for exactly the reason
W4.2's (u, v) GP describes it on visibilities. Instrumental errors
(chromatic dispersion, piston, calibrator errors) are smooth in λ and t
rather than in (u, v); neither form exposes those as axes, and they are
the Phase 5 joint-noise question, not this one.

**(d) Wavelength.** *Superseded by §3.6 (Peter's question of 2026-09-11):
the claim below holds for an achromatic sky error only, and the chromatic
case is the common one.* Because the axes are already B/λ, the four-axis form
absorbs the spectral dimension without adding one:

- For an **achromatic sky** (the W4.4 study; an `Image` channel) nothing
  more is needed. The closure phases of one triangle across its spectral
  channels lie on a ray in the four-dimensional space, the GP sees them as
  neighbours along it, and the wavelength correlation of a sky-model error
  is captured automatically because along that ray the bispectrum *is* a
  function of one variable.
- For a **chromatic sky** (a `Cube` channel, sketch Q4) the Fourier step
  must know each sample's wavelength to pick the cube plane. That is the
  per-sample frequency in `extra_coords` with units (sketch Q2, in W4.2's
  scope): visible to the step, invisible to the kernel, which is the right
  division.
- Making λ a *fifth axis* would break `VisibilitySet`'s frozen convention
  for the visibilities beside it, and mix units, which the kernel refuses
  by name. Multi-window observations are one channel per spectral window
  (Q2's ruling), where each container is monochromatic in the instrument's
  sense and the (u, v) are still per-sample.

**(e) Cost, since a dense solver is the only one available in four
dimensions.** Illustrative sizes for a four-telescope array:

| Spectral channels | Exposures | Closure phases | Dense GP per evaluation |
|---|---|---|---|
| 6 (low resolution) | 10 | 240 | negligible |
| 200 (medium) | 10 | 8 000 | seconds; workable for emcee, slow for NUTS |
| 1 700 (high) | 10 | 68 000 | not feasible (a 37 GB matrix) |

The large case is handled by splitting into one dataset per exposure or
per spectral window (each negotiated onto the one `sky` channel, so the
model is still built once) or by the approximate-GP strategy slots of
Phase 5 (HSGP/EFGP); the four-axis form leaves both open, the label form
neither.

### 3.4 One requirement the four-axis form adds

A triangle has three representations by two of its baselines, and the
same closure phase listed as (u_ij, v_ij, u_jk, v_jk) in one file and as
(u_jk, v_jk, u_ki, v_ki) in another would look far apart to a GP. The
kind must fix a **canonical ordering** — baselines ij and jk for telescope
indices i < j < k, the third baseline implied as −(u1 + u2, v1 + v2), and
the sign convention of the phase stated — in its docstring, and the
`ClosurePhase` step and any reader must produce it. The Fourier step does
not care about the ordering; the kernel and the coordinate-equality
pairing with the visibility dataset (gap I-2's bit-identical rule, which
is also horizon (h)'s hook) both do.

### 3.5 Recommendation (as first written; revised in §3.6)

**Four axes, as drafted**, with the canonical ordering of §3.4 in the
kind's docstring and asserted by W4.1's conformance row; the `triangle`
label kept in `extra_coords` as the human-readable identity; the
per-sample frequency added there by W4.2. The dimensionality consequence
to record in the decision log: the kernel's coordinate space is
four-dimensional and single-unit, `DenseGP` only; wavelength is absorbed,
not added; per-axis (ARD) length scales are a W4.5 product-kernel option
on the dense path if isotropy proves too strong.

### 3.6 Correction: a chromatic sky error is correlated in wavelength, and (u, v) alone cannot express it

Peter's question (2026-09-11): if the model is missing spectral lines or
molecular bands that contribute at specific wavelengths in specific
patches of sky, is the residual correlated across *wavelength* rather than
across (u, v), or are the two so tightly coupled that it does not matter?

It matters, and the memo's §3.3(d) was wrong for this case. Write the
missing component as δI(x, y, λ) = P(x, y) S(λ) — a patch with a spectral
profile (a band of width Δλ at λ₀). Its effect on a baseline **B** at
wavelength λ is

δV(**B**, λ) = S(λ) · F(**B**/λ),

with F the Fourier transform of the patch, smooth in (u, v) on the scale
1/θ of the patch. The residual is therefore a product of two functions on
two *different* coordinates: sharp in λ, smooth in (u, v). Projecting onto
(u, v) alone destroys the first:

- Along one baseline's spoke, λ and |u| are monotonically related, so the
  band appears as a bump in |u| of width u₀ Δλ/λ₀ at radius B/λ₀.
- On a *different* baseline the same band sits at a different radius,
  B′/λ₀. Two points close in (u, v) — a 50 m and a 100 m baseline at
  similar position angle, at the same |u| — are at wavelengths a factor of
  two apart, one in the band and one out of it.

So an isotropic (u, v) kernel with one length scale is asked to do two
incompatible things: correlate along a spoke on the band's scale, and
correlate across the plane on the patch's scale, while *not* correlating
neighbouring points from different baselines that are at different
wavelengths. It cannot. The coupling u = B/λ absorbs wavelength only when
S(λ) is constant — a grey error, a missing continuum component — which is
the case the W4.4 study as drafted tests, and not the case Peter names,
which is the common one in spectro-interferometry (Brγ, CO band heads,
the silicate feature) and the case the M2 study was built around for
spectra.

**Consequence.** The kernel must see wavelength as a coordinate of its own,
separate from (u, v), and the natural covariance is a product,
k_uv(u, v) · k_λ(λ), with two length scales (patch size, band width). Given
§3.2 (a kernel sees axes only), that means:

1. **`VisibilitySet` gains a third axis, `spectral_axis`** (length,
   frequency or energy, `Order.ANY` since the layout is `POINTS` and each
   sample carries its own wavelength). Radio spectral windows and every
   modern optical beam combiner are dispersed, so the monochromatic case
   is the degenerate one (a constant column), not the other way round.
   This amends a frozen kind: a decision-log row and a conformance update
   in the same PR (ground rule 9). The blast radius measured at `59454ea`
   is small — the constructor is called positionally in about ten test
   sites (`tests/conformance/{composition,test_schema,test_results}.py`,
   `tests/results/{test_plots,test_results}.py`,
   `tests/backends/test_jax.py`) and one doctest in
   `ampere/results/serialisation.py`; no shipped step or model produces
   one yet. The alternative, a second dispersed kind beside a
   monochromatic one, doubles what every step and family must accept for
   no benefit.
2. **`ClosurePhases` is five axes**: `(u1, v1, u2, v2, spectral_axis)`,
   the same product structure, k(u1, v1, u2, v2) · k_λ(λ).
3. **A kernel must be able to act on a subset of a container's axes.**
   The single-unit rule (`likelihood.py:720-745`) refuses (u, v, λ) for an
   isotropic kernel, and its own message points at the fix: "declare a
   kernel that takes a length-scale per axis". `KernelSpec` (and every
   `Kernel`) gains an `axes` selector; `Matern32(axes=("u", "v"))` is the
   isotropic (u, v) kernel of W4.2 unchanged, `Matern32(axes=("spectral_axis",))`
   the band kernel, and W4.5's `Product` composes them. This is one field
   on the spec, checked in `check_compatible` against the subset's units,
   and it is what lets W4.2's achromatic flagship and the chromatic case
   coexist on one kind. It belongs to W4.5 (which owns the kernel section)
   and W4.2 consumes it, so the dependency W4.2 → W4.5 already drafted is
   the right order.
4. **Sketch Q2 lapses for this purpose.** The per-sample frequency was
   going into `extra_coords` with units so the Fourier step could pick a
   `Cube` plane; with wavelength an axis, the step reads it from the axis
   and the requirement it publishes on the cube's `spectral_axis` is
   `points=` the container's unique wavelengths. `extra_coords` units may
   still be wanted elsewhere but are no longer on Phase 4's critical
   path; W4.2's text drops that clause.
5. **Time stays out, deliberately.** The same argument applies to
   instrumental errors smooth in time (piston, coherence loss), and the
   same mechanism (an axis plus a product kernel) would serve; but those
   are not sky misspecification, they are Phase 5's joint-noise question
   (horizon (h)), and adding a `time` axis now would be a guess at a
   contract nobody has exercised. The `triangle`/baseline labels and an
   epoch label in `extra_coords` keep the structure recoverable.

**Dimensionality after the correction.** Kernel coordinates: 3 for
visibilities, 5 for closure phases, each split by the product into a
2-D (or 4-D) single-unit block and a 1-D spectral block; hyperparameters:
one extra length scale and amplitude per block. Data count unchanged.
Solver: `DenseGP` only, as before. **The scale route this opens** is worth
recording for Phase 5: in (B, λ) coordinates a dispersed observation is a
*product* of a set of baseline-epochs and a spectral axis, so a product
kernel's covariance is a Kronecker product K_B ⊗ K_λ, with K_B small and
K_λ quasiseparable over the ordered spectral axis — exact and close to
O(N) until masks break the structure. The kind must therefore keep the
baseline × channel structure recoverable (a baseline label per sample in
`extra_coords`, the wavelength as an axis), which the design above does.
This is the reduced-rank/structured-solver bullet of Phase 5 with a
concrete first customer; it is not Phase 4's.

**The W4.4 study should include the chromatic case.** Scenario (b), "the
disc omitted", is a grey error. Add a scenario in which the omitted
component is a compact patch with a band profile (the interferometric
analogue of M2's missing feature), fitted under (u, v)-only, λ-only and
product kernels, so the claim that the product is needed is measured
rather than argued. That is the study's real question for this
observable.

### 3.7 Revised recommendation for D2

Four baseline axes plus a spectral axis for `ClosurePhases` (five in all),
and `VisibilitySet` amended to three; the canonical baseline ordering of
§3.4; `triangle` and baseline labels in `extra_coords`; the `axes`
selector on kernels in W4.5 with `Product`; the amendment carried by W4.1
(the kind) with its decision-log row and conformance update; W4.2's
extra-coords-units clause dropped; the chromatic scenario added to W4.4.

---

## 4. Corrections to the drafted items from the surveys

To apply to `WORK_ITEMS.md` once D1 is ruled (they do not depend on it
unless noted):

1. **W4.1** — "the first two `configure_from` instances" → "the first
   cross-kind uses of `configure_from` (`LSFConvolution` already uses it
   within a kind on all three backends)". The four-axis `ClosurePhases`
   plus the `triangle` label is the assumption; state it once.
2. **W4.0** — add: the stale sentence in `configure_from`'s docstring
   (`core/transform.py:795-797`) corrected. Item (4), the architecture
   line, follows D1's ruling.
3. **W4.3** — "lowered through the registry like every standard step" →
   "associated with the reference step as every standard step is: the
   same class name in the backend's module, the `apply_flux` native
   surface, and `BACKEND` declared; there is no step registry". State
   which twin pattern (torch re-declares, jax inherits) the new steps
   follow, and why.
4. **W4.4** — the template page cites the photometry + spectrum
   composition as the simple case, and no such example exists. The owed
   "photometry + spectra composition smoke test" and Phase 6 tutorial
   should be pulled forward as a Sonnet filler before W4.4, or W4.4 should
   own it. Recommendation: a small item, W4.11, dispatched with W4.0.
5. **W4.8** — the annotation pass should include `spectrum_photometry.md`'s
   pre-shipping `SyntheticPhotometry` signature.
6. **Ingestion** — no reader in Phase 4 (W4.4's study is synthetic); note
   OIFITS as the first reader for Phase 6 and, under D, the first content
   of `ampere.interferometry`.

---

## 5. What does not change under any ruling

- The steps (`FourierSample`, `ClosurePhase`, the smearing steps) go under
  `backends/<backend>/interferometry.py` on all three backends.
- `VisibilitySet` is unchanged; amplitudes and squared visibilities are a
  real-valued `VisibilitySet`; the Rice, von Mises and circular-GP rulings
  of 2026-09-03 stand.
- The out-of-tree proof stays the polarimeter test; Phase 4 proves
  composition.

## 6. Rulings asked for

1. **D1 placement**: the kind in core (recommended), or in a separate
   namespace; and whether a per-observable top-level package
   (`ampere.interferometry`) is wanted now, later with the first reader
   (recommended), or not at all. If a grouping word is wanted after all,
   which of `obs`, `astro`, or another.
2. **D2** (revised in §3.6–3.7): `VisibilitySet` amended to `(u, v,
   spectral_axis)` and `ClosurePhases` as `(u1, v1, u2, v2,
   spectral_axis)` with the canonical ordering; kernels gain an `axes`
   selector so a product of a (u, v) block and a spectral block is
   expressible (recommended); or the label form with the consequence that
   no GP can be composed on closure phases.
3. Whether the photometry + spectrum composition example (§7.1, which
   runs today) becomes a small Phase 4 item ahead of W4.4 — the script,
   its smoke test and the owed tutorial page, plus the two reference-path
   findings it exposed (recommended).
4. (Added by §7.2) whether a latent GP on closure phases — `VonMisesFamily`
   consuming a latent GP on the native path, the Poisson pattern — is
   drafted as a Phase 4 item after W4.3 or left to Phase 5.

---

## 7. Two worked compositions: one model, two datasets (Peter, 2026-09-11)

Peter asked for a walk-through of a fit that targets a combined dataset with
a single model, both to make the third ruling concrete and to test whether
the D2 reasoning has gone down a blind alley. The first composition runs
today against `59454ea`; the second is the interferometric one written in
the same shape under §3.7, with every piece that does not yet exist marked.

### 7.1 A spectrum and photometry on one model channel — runs today

```python
import numpy as np, scipy.stats as st, astropy.units as u
from ampere.core import (Dataset, DatasetCollection, FittingProblem, Instrument, Likelihood,
                         GaussianFamily, IndependentNoise, Spectrum, PhotometricPoints, negotiate)
from ampere.backends.reference import ModifiedBlackBody, LSFConvolution, Resample, SyntheticPhotometry

# 1. One model, one channel. The grid it is given is a fallback; negotiation replaces it.
model = ModifiedBlackBody(np.geomspace(1.0, 100.0, 3000),
                          temperature=st.uniform(50.0, 400.0), scale=st.loguniform(1e-16, 1e-14),
                          beta=st.uniform(0.5, 2.5), channels="sed")

# 2. Two instruments reading the same channel. Coordinates that must match the
#    observed data (the spectrograph's grid) are the observed container's own.
obs_wave = np.linspace(5.0, 35.0, 120)
spectrograph = Instrument([LSFConvolution(resolving_power=100.0), Resample(obs_wave)],
                          channel="sed", label="irs")
camera = Instrument([SyntheticPhotometry.from_library(
                        ["2MASS_J", "2MASS_Ks", "WISE_RSR_W1", "WISE_RSR_W3", "WISE_RSR_W4", "SPITZER_MIPS_70"],
                        np.geomspace(1.0, 100.0, 500))],
                    channel="sed", label="catalogue")

# 3. Observed containers, built by hand from arrays (there is no reader in v2).
observed_spectrum = Spectrum(obs_wave * u.um, flux * u.Jy, uncertainty=sigma * u.Jy)
observed_photometry = PhotometricPoints(filters, pivots * u.um, phot * u.Jy, uncertainty=phot_sigma * u.Jy)

# 4. Datasets bind observed + instrument + likelihood. Two instruments on one
#    channel must carry distinct labels, or the collection refuses by name.
gauss = lambda: Likelihood(GaussianFamily(), IndependentNoise())
datasets = DatasetCollection({"irs": Dataset(observed_spectrum, spectrograph, likelihood=gauss()),
                              "catalogue": Dataset(observed_photometry, camera, likelihood=gauss())})

# 5. The problem negotiates once, compiles the model onto the union grid, validates.
problem = FittingProblem(model, datasets, seed=20260911)
problem.requirements["model"]["sed"].sources          # ('irs', 'catalogue')
problem.parameters.free_names                        # ('model.temperature', 'model.beta', 'model.scale')

# 6. Any engine.
from ampere.inference import EmceeEngine
run = EmceeEngine(problem, walkers=16).run(250, burn_in=125)
```

What it printed, and what each line teaches:

```
negotiated channels: ['sed']
sources asking of channel 'sed': ('irs', 'catalogue')
  axis spectral_axis: <AxisRequirement 'spectral_axis' micron [4.89383,35.7432]@R>=706.446 [5,35]@step<=0.12605 500 point(s)>
free parameters: ('model.temperature', 'model.beta', 'model.scale')
```

- The LSF asked for its padded range at three times its own resolving
  power (it learned the range from the `Resample` after it through
  `configure_from`); the resampler asked for its own interval at half its
  bin; the photometry step asked for the exact points it tabulated its
  responses on. One union grid served all three, and the model was built
  once per draw.
- The user wrote no grid arithmetic and no `reloadFilters`. The legacy
  ordering trap is gone by construction.
- Three free parameters: the instruments here have none. A calibration
  factor on the spectrum would be a `CalibrationScale` step and would appear
  as `irs.instrument.calibration_scale.scale`.

**Two findings from running it**, both out of this memo's scope and recorded
for the owed list:

1. Applying an instrument to a model that has *not* adopted the negotiated
   grid is refused by name — the photometry step asks for the points it
   published and finds the model's fallback grid instead. That is correct
   (the message even says "a negotiation defect rather than a usage
   error"), but it means synthetic data for an example must be generated
   through `negotiate` + `compile_for` (or through the problem), and no
   page says so. The template page needs that sentence.
2. **The reference `LSFConvolution` rebuilds its dense `(n, n)` matrix on
   every evaluation** (`backends/reference/instrument.py:407-425`; the
   torch twin builds it once as a buffer). On the union grid the
   photometry step's tabulation dominates *n*, so with a 4 000-point
   tabulation one likelihood evaluation took 2.5 s (5 647 grid points;
   2.2 s of it in `influence`), which is unusable even for an oracle. With
   500 tabulation points it is a quarter of a second. A cached operator
   keyed on the grid's identity is a small reference-backend fix
   (Sonnet-sized) and should precede any example that fits a spectrum
   through an LSF on the reference path.

### 7.2 Visibilities and closure phases on one sky channel — the same shape under §3.7

Pieces that exist today are unmarked; **[W4.x]** marks what the item builds;
**[sketch]** marks a signature invented here and not yet designed.

```python
from ampere.core import VisibilitySet, ClosurePhases                      # ClosurePhases [W4.1]; VisibilitySet + spectral axis [W4.1]
from ampere.core import ComplexGaussianFamily, VonMisesFamily, GaussianProcessNoise, DenseGP, Matern32, Product
from ampere.backends.reference import Binary, FourierSample, ClosurePhase, BandwidthSmearing   # [W4.1]

# 1. One model, one channel "sky", emitting an Image(x, y) in mas.
model = Binary(separation=st.uniform(2.0, 20.0), position_angle=st.uniform(0.0, 2 * np.pi),
               flux_ratio=st.loguniform(0.01, 1.0), channels="sky")          # [W4.1]

# 2. Observed containers from the OIFITS tables (by hand; a reader is Phase 6).
#    Every sample carries its own (u, v) in wavelengths AND its wavelength.
vis = VisibilitySet(u_pts, v_pts, wave * u.um, visibility,
                    uncertainty=sigma_vis, extra_coords={"baseline": baseline_labels})
t3 = ClosurePhases(u1, v1, u2, v2, wave3 * u.um, phase * u.rad,             # canonical ordering (§3.4)
                   uncertainty=sigma_phi * u.rad, extra_coords={"triangle": triangle_labels})

# 3. Two instruments on the one channel. The Fourier step takes its (u, v, λ)
#    buffers FROM the observed container (gap I-2), never recomputes them, and
#    publishes x/y intervals of ±fov/2 with max_step = 1/(2 u_max) (gap I-4).
vis_instrument = Instrument(
    [FourierSample.from_observed(vis, field_of_view=60 * u.mas),            # [sketch] constructor
     BandwidthSmearing(resolving_power=22.0)],                              # [W4.1]: asks FourierSample, via configure_from, for the extra samples it averages over
    channel="sky", label="gravity_vis")
t3_instrument = Instrument(
    [FourierSample.from_observed(t3, field_of_view=60 * u.mas),             # three baselines per triangle
     ClosurePhase()],                                                       # [W4.1]: three visibilities -> one angle in (-π, π], mask propagated
    channel="sky", label="gravity_t3")

# 4a. Independent likelihoods: the baseline fit.
independent = DatasetCollection({
    "vis": Dataset(vis, vis_instrument, likelihood=Likelihood(ComplexGaussianFamily(), IndependentNoise())),
    "t3":  Dataset(t3,  t3_instrument,  likelihood=Likelihood(VonMisesFamily(), IndependentNoise())),   # VonMises IMPLEMENTED flips [W4.1]
})
problem = FittingProblem(model, independent, seed=20260911)
problem.requirements["model"]["sky"].sources     # ('gravity_vis', 'gravity_t3'): the pairing horizon (h) wants, visible
problem.requirements["model"]["sky"].axes         # x, y: one image grid from the union of both Fourier steps' requirements

# 4b. The flexible likelihood on the visibilities — the chromatic case of §3.6.
sky_error = Product(                                                        # [W4.5] Product, axes= selector
    Matern32(axes=("u", "v"), amplitude=st.halfnorm(scale=0.05), length_scale=st.loguniform(1e5, 1e7)),   # smooth in spatial frequency: the patch
    Matern32(axes=("spectral_axis",), length_scale=st.loguniform(0.005, 0.1)),                            # sharp in wavelength: the band, µm
)
flexible = DatasetCollection({
    "vis": Dataset(vis, vis_instrument,
                   likelihood=Likelihood(ComplexGaussianFamily(), GaussianProcessNoise(sky_error, DenseGP()))),   # [W4.2] the circular closed form
    "t3":  Dataset(t3, t3_instrument, likelihood=Likelihood(VonMisesFamily(), IndependentNoise())),
})
```

What is the same as §7.1, which is the point of the phase: the model,
the two instruments on one channel with distinct labels, the observed
coordinates taken from the container, one negotiation, one compile, any
engine. What is new is exactly the list the plan wanted stressed: a
kind-changing and coordinate-changing step in the middle of the chain, a
complex container, a wrapped family, and `configure_from` used across a
kind change.

**Where the walk-through found the blind alley.** Step 4b puts the GP on
the *visibilities* only. The memo's D2 argument (§3.3(c)) was that the
axes make a GP on *closure phases* possible. Possible for the kernel, yes;
but a GP added to a **wrapped** observable is not analytically
marginalisable — the von Mises family with a GP is a `LATENT` composition,
and today only `PoissonFamily` declares `CONSUMES_LATENT_GP = True`
(`core/likelihood.py:2574`). So a GP on closure phases is reachable only by
giving `VonMisesFamily` the latent-conditional `log_prob` and running on
the native backends under NUTS or VI, the Poisson pattern — never on the
numpy path with emcee. That is a real and useful thing, and the five-axis
kind is what makes it *expressible*, but it is not what Phase 4's flagship
test runs. Corrected claim for D2: the axes decide whether a GP on closure
phases can ever be composed; the flagship GP of W4.2 and W4.4 is on the
visibilities; a latent GP on closure phases is a candidate item after
W4.3 (native path exists) and should be listed, not assumed.

**What the third ruling was asking.** Only this: §7.1 is the simplest
combined fit the package supports, nothing in `examples/` or `docs/source/`
shows it, and W4.4's template page will cite it as "the case you already
know" before showing §7.2. The ask is whether to write §7.1 up as a small
item (the example script, its `tests/examples` smoke test, and the tutorial
page `tutorials.rst` says is owed) *before* W4.4, so the template page has
something to point at. The two findings above are what running it turned
up, which is the argument for doing it early.

