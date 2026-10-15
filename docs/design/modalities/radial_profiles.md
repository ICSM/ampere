# Radial profiles and azimuthal symmetry — the design, explained for the D1 ruling

*Fable, 2026-10-15, for Peter's D1 of Phase 8. Written after his questions of
the same day: what exactly is proposed, why it is shaped like that, and what
the real covariance of an averaged profile involves. The item it informs is
W8.2 (`WORK_ITEMS.md`, Phase 8); the horizon memo's §3 is the one-paragraph
origin. The numbers in §4 come from a numpy probe reproduced in Appendix A.*

## 0. The short answer

A radially symmetric source is a one-dimensional model, `I(r)`, observed
through instruments that are not one-dimensional: a map with a beam, a set
of visibilities at `(u, v)`, a photometric aperture. The proposal is to
make `I(r)` a first-class model output — a kind, `RadialProfile` — and to
supply the *linear* steps that take it to each observable: `Render` to an
image (so the beam, the image likelihoods and the image GP solvers that
already exist apply unchanged), `HankelTransform` to visibilities (the
symmetric source's Fourier transform is one-dimensional), and, for the
case where someone hands you a profile rather than a map, `AzimuthalAverage`
from an image to a profile with **its covariance propagated exactly**
through the averaging. Nothing about the likelihoods, the engines or the
GP machinery changes; the GP goes where it always goes, and a second,
model-side use of a GP — a non-parametric profile — falls out of the
parameter contract as it stands.

The priority Peter gave (2026-10-15): continuum radial profiles of NESS
sub-mm maps and ALMA-ACA uv tables used through the uv distance alone,
first; elliptical averaging required but after circular; sector averaging
not a priority. §5 orders the item accordingly.

## 1. The use cases, concretely

| Case | The data | The model | The route (§2) |
|---|---|---|---|
| **NESS continuum maps** — AGB/RSG dust shells in JCMT SCUBA-2 450/850 µm maps (beam 8″/13″, 2–4″ pixels), with a variance map from the reduction | an `Image` with a known beam and per-pixel variance | a projected dust shell: a 1-D density law `ρ(r) ∝ r^-p` between inner and outer radii, projected along the line of sight to `I(r)` (closed form for a power law; numerical for a general law); or a 1-D RT code's `I(r)` | A (fit the map) or B (fit the averaged profile with its covariance) |
| **NESS / ACA uv tables, uv distance only** — ALMA-ACA Band 6/7 continuum, the source unresolved in azimuth, so only `|uv|` matters | a `VisibilitySet`; `(u, v)` per sample, but the model depends on `ρ = √(u²+v²)` alone | the same shell, Hankel-transformed | C (fit the visibilities) |
| **ALMA discs** — ring-and-gap protoplanetary discs at 1 mm, inclined | a `VisibilitySet` from a uv table | `PowerLawDisc` with gaps, `GaussianRing`s, or the non-parametric `LatentProfile` (frank's problem) | C with `Deproject` |
| **VLTI** — a star's angular diameter with limb darkening, or a shell around it | squared visibilities and closure phases through `read_oifits` | a limb-darkened disc profile; a star plus a thin shell | C (`HankelTransform` on `UniformDisc`'s territory, generalised) |
| **Galaxies** — a surface-brightness profile from an image | an `Image` or a published profile | `Sersic`, sums of them | A or B |

Two things these share. The *model* is one-dimensional in every case — that
is what "azimuthal symmetry" buys — and the *data* are whatever the
telescope produced. The design keeps those two facts separate: the kind
describes the model's output, the steps describe the telescope.

## 2. One NESS shell, three ways

Take a SCUBA-2 850 µm map of an AGB star: a 128 × 128 map at 2″ per pixel,
a 13″ beam, a variance map. The physics question is the density law of the
shell. There are three honest ways to fit it, and the proposal supports all
three because different users will have different inputs in hand.

**Route A — fit the map itself.** The model is `ProjectedShell(r)` →
`Render(x, y)` → `PSFConvolution(beam)` → the `Image` likelihood with the
variance map as `IndependentNoise`, and, if wanted, a 2-D GP on the image
for the background the model does not describe (W5.5's pattern; the
reduced-rank solvers for a large map). This is the statistically cleanest
route: every pixel is used with its own variance, the beam is applied
exactly, nothing is averaged, and a non-symmetric residual (a companion, a
bow shock) shows up in the image residual plot where it is. The radial
profile is then a *diagnostic* — the plotting side averages the data, the
model and the residual in annuli for the figure — not the fitted
observable. The cost is a likelihood over `N_pix` points, which for a
128 × 128 cutout is 16 384 — fine with `IndependentNoise`, and the HSGP /
EFGP solvers are what W5.5–W5.6 built for the GP case. **Why this route
still needs the kind.** Without `RadialProfile` and `Render`, the user must
write the shell as a 2-D `_ImageModel` evaluating `I(√(x²+y²))` on the pixel
grid — which works today, and is how `UniformDisc` is written — but every
symmetric model then re-implements the projection and the rendering, the
1-D structure is invisible to the chain (no Hankel route, no profile
figure), and the non-parametric profile of §3 cannot be declared at all.

**Route B — fit the azimuthally averaged profile.** The observed side
averages the map in annuli to a `RadialProfile` with `N_bins` points (15
bins of half a beam in the probe). The model side must do *the same
averaging to the same model*, so the chain is `ProjectedShell` → `Render`
→ `PSFConvolution` → `AzimuthalAverage(edges)`, or — cheaper and exactly
equivalent for a symmetric model and a symmetric beam — the 1-D
convolution `BeamConvolution(beam)` in `r` (a convolution of two radially
symmetric functions is a product in Hankel space, so it is two Hankel
transforms or one small dense matrix) followed by binning. Route B is what
a user does when the profile is what they have (a published table, a
collaborator's extraction), when the map is too large to fit directly,
or when they want a profile *figure* with error bars that mean something.
Its one trap is §4: the bins are correlated by the beam and their naive
errors are wrong by a large factor.

**Route C — fit the visibilities.** For the ACA table the model is
`ProjectedShell` → `HankelTransform` evaluated at each sample's
`ρ = √(u²+v²)` (`Deproject` with inclination zero is the identity, so an
unresolved-in-azimuth source needs no deprojection step at all; a disc
does). The data are the measured visibilities with their weights; no
averaging, no beam, the primary-beam attenuation a multiplicative factor
on the rendered model if the source is extended enough to need it. If the
table has already been binned in `ρ`, the bins are independent (thermal
noise is independent per visibility and binning averages independent
numbers), so the diagonal covariance from the weights is exact — the probe
shows this — and Route C stays clean either way.

**What the three share.** One model class, written once as `I(r)`. The
chain decides the observable. The GP likelihood is declared on the
observed container exactly as today: on the image (A), on the profile (B),
on `(u, v)` or `ρ` (C). The same shell can be fitted all three ways on one
object and the posteriors compared — which is the conformance row the item
needs (A against C on a synthetic shell agree; B against A agree once B's
covariance is right).

## 3. Why the pieces are shaped like this

**A kind, not just models.** `ampere.core`'s containers carry coordinates
and units and are what steps accept and produce (`ACCEPTS` / `PRODUCES`).
Making `RadialProfile` a kind — one axis `r` (angle), values a surface
brightness (Jy/sr or Jy/beam, the unit carried), `Layout.POINTS` — is what
lets `Render`, `HankelTransform`, `BeamConvolution` and `AzimuthalAverage`
be ordinary steps that the chain composes and negotiates
(`requirements` for the `r` grid the quadrature wants, as `Resample` and
`LineIntegration` negotiate theirs), and lets a profile be an *observed*
container with a mask, uncertainties and a covariance. Three class
attributes, as `LineFluxes` was at W7.9.

**Steps are linear operators with influence matrices.** `Resample`
exposes `influence(source)`, the `(n_out, n_in)` matrix taking the input
onto the output; `LineIntegration` has one fixed influence matrix;
`SyntheticPhotometry` the same. Every step here is linear in `I(r)`:
`Render` (interpolation of `I` at each pixel's radius — a sparse matrix),
`HankelTransform` (`V(ρ_k) = 2π Σ_j w_j I(r_j) J₀(2π ρ_k r_j) r_j` — a dense
`N_ρ × N_r` matrix, built once per grid and reused, W5.17's planning
pattern), `BeamConvolution` (a dense `N_r × N_r` matrix), `AzimuthalAverage`
(the `(N_bins, N_pix)` averaging matrix `A` of §4). Linearity is what makes
the covariance propagation of §4 exact and what makes the non-parametric
profile tractable: the whole chain from the knot values to the observable
is one matrix, so a Gaussian prior on the knots is a Gaussian predictive on
the data and the sampler sees a well-conditioned problem.

**`Deproject` is a step on the visibilities, not a parameter of the
model.** Inclination, position angle and offset describe how the symmetric
source sits on the sky, and they act on the *sampling* — rotate and scale
`(u, v)` to the disc's frame, apply the phase ramp of the offset — not on
`I(r)`. Written as a step they apply to any profile model and to
`LatentProfile`; written into each model they would be repeated in each.
Their parameters live on the step, as `CalibrationScale`'s does, and are
the first natural `Derived` customers: `cos i` from an axis ratio, a ring's
width from a radius and a fractional width.

**`HankelTransform` is a step, not a method of each model.** One quadrature
serves every profile, including the non-parametric one; a model that has a
closed-form transform (`UniformDisc` → the Airy pattern, which exists as
`UniformDiscVisibilities`) is the conformance oracle for the quadrature,
as `GaussianSourceVisibilities` is for the Gaussian. The quadrature itself
is unproblematic: on a 2 000-point `r` grid the trapezoid transform of a
Gaussian matches the closed form to 8 × 10⁻⁷ (Appendix A); the discrete
Hankel transform is the large-`N_r` route if a profile ever needs ten
thousand knots, and is not needed for the first cut.

**`Render` makes the image route free.** `PSFConvolution`, the `Image`
likelihoods, `propagate_mask_grid`, the 2-D GP solvers and the image plots
all exist. A symmetric model rendered onto the observed map's own `(x, y)`
axes inherits every one of them, which is why Route A costs the item one
small step rather than a modality.

**The GP has two homes, and only one is new.** On the observed container —
the image, the profile, the visibilities — `GaussianProcessNoise` is
unchanged: it absorbs what the parametric shell gets wrong and its
localisation plot says where (in `r` for Route B; in the image for A; in
`ρ` for C). The second home is the model: `LatentProfile` is a profile
whose knot values are `m(r) + L(θ) z`, with `z` an array-valued parameter of
shape `(N_knots,)` under a standard-normal prior and `L(θ)` the Cholesky
factor of a kernel on `log r` with its amplitude and length scale as
hyperparameters. This is the non-centred GP that `latent_transform` already
implements for latent likelihood families, used on the model side; the
parameter contract supports array-valued parameters with independent priors
(`parameters.md` §3, "array-valued parameters are supported"; a genuinely
multivariate prior is out of scope, which is exactly why the non-centred
form is the right one). It needs **no contract change**. Used alone it is
frank's fit as a declaration; added to a parametric shell it is a
misspecification correction *in the model's own coordinate* whose
posterior, pushed through the chain, is the "where does the symmetric
model fail, in radius" figure even when the data are visibilities.

## 4. The covariance an averaged profile really has (Peter's question 5)

**The mathematics.** Averaging is linear: `P = A m`, where `m` is the map
as a vector and `A` the `(N_bins, N_pix)` matrix with `1/N_i` in row `i` for
the pixels of annulus `i`. If the map's noise has covariance `C`, the
profile's covariance is **exactly** `Σ_P = A C Σ_pixel A^T`. For a map made
from white noise smoothed by a Gaussian beam of width `s` (the SCUBA-2 case
after the reduction's smoothing; interferometric images likewise), `C_pq =
σ_p σ_q exp(−d_pq² / 4s²)` with `d_pq` the pixel separation and `σ` from
the variance map — a closed form, and `Σ_P` is a 15 × 15 matrix that costs
a few hundred thousand exponentials to build. Per-pixel variance maps
(SCUBA-2 produces one), a mask and a non-Gaussian beam (`C` built from the
beam's autocorrelation) fit the same formula.

**What the naive error misses.** The probe (Appendix A: 128² map, 2″
pixels, 13″ beam, 6″ bins, unit smoothed rms):

| bin | r (″) | pixels | sd, naive (σ/√N) | sd, exact | ratio | corr with next bin | corr two bins on |
|---|---|---|---|---|---|---|---|
| 0 | 3 | 32 | 0.177 | 0.863 | 4.9 | 0.88 | 0.45 |
| 2 | 15 | 144 | 0.083 | 0.454 | 5.4 | 0.77 | 0.34 |
| 6 | 39 | 376 | 0.052 | 0.276 | 5.4 | 0.77 | 0.35 |
| 10 | 63 | 608 | 0.041 | 0.217 | 5.3 | 0.77 | 0.35 |
| 14 | 87 | 832 | 0.035 | 0.185 | 5.3 | — | — |

A 400-realisation Monte Carlo reproduces the exact column to a few per
cent. So with half-beam bins the per-bin error is **five times** the
independent-pixel estimate at every radius, and **adjacent bins are 77 %
correlated**. A fit that uses the naive diagonal would report a χ² of
order 25 per point on a perfect model, and a GP likelihood asked to make
sense of that would do the only thing it can: fit a kernel with the beam's
length scale and an amplitude that soaks up the instrument correlation —
after which the misspecification diagnostics are describing the beam, not
the physics. The exact covariance is cheap and removes the problem at its
source; the GP is then free to do its job.

**What ampere needs for it.** Two small things, both inside W8.2's scope.
(1) `AzimuthalAverage.covariance(variance_map, beam, mask=None)` — the
step already holds `A`; it returns `A C A^T` for the map's noise model
(white-plus-beam by default; a user-supplied `C` or beam autocorrelation
otherwise). (2) A way to declare that covariance on the observed profile.
Today's noise components are `IndependentNoise` (a diagonal),
`GaussianProcessNoise` (a kernel with free hyperparameters) and
`JointGaussianProcessNoise`; W5.8 made a *sum* of components the
likelihood's covariance (`S = Σ components + diag σ²`). The missing piece
is a **`FixedCovarianceNoise(matrix)`** component — a given positive
semi-definite matrix, no parameters, one Cholesky at construction — so
that a Route B likelihood reads `Likelihood(GaussianFamily(),
FixedCovarianceNoise(Σ_P) + GaussianProcessNoise(Matern32(...)))`: the
instrument's correlation stated, the GP for the physics. That is a §4.4
addition (one component, no change to any existing one), with its
decision-log row and conformance rows (a profile fitted with the exact
covariance recovering the injected shell's parameters and their coverage
where the naive diagonal does not — the probe's factor five as a test).
It is also useful well beyond this modality: any published table with a
covariance (a binned power spectrum, a Gaia solution's correlations) is the
same declaration.

**And for the uv case** (Peter's ACA tables): nothing. Unbinned visibilities
carry independent weights; binned ones average independent numbers, so the
per-bin `1/√Σw` is exact (Appendix A's second block) and `IndependentNoise`
is right. Only the *image* route creates correlated bins, because the beam
mixes neighbouring pixels before anything is averaged.

## 5. What a user writes — five declarations

Illustrative; names follow the reference backend's conventions (`as_parameter`,
`Dataset`, `Likelihood`); the units and grids are the observed container's.

**(1) A NESS shell, the map fitted directly (Route A).**
```python
shell = ProjectedShell(r_inner=..., r_outer=..., index=..., flux=...)   # I(r), closed form for a power law
chain = [Render(x=map.x, y=map.y), PSFConvolution(beam_fwhm=13 * u.arcsec)]
likelihood = Likelihood(GaussianFamily(), IndependentNoise())           # the variance map is map.uncertainty**2
ds = Dataset(observed=scuba2_map, model=shell, instrument=chain, likelihood=likelihood)
# diagnostics: plot_radial_profile(tree, edges=...) averages data, model and residual in annuli for the figure
```

**(2) The same shell from its averaged profile, with the covariance (Route B).**
```python
step = AzimuthalAverage(edges=np.arange(0, 90, 6) * u.arcsec, centre=(x0, y0))
profile = step.apply(scuba2_map)                                           # the observed RadialProfile
sigma_p = step.covariance(scuba2_map.uncertainty**2, beam_fwhm=13 * u.arcsec)
chain = [BeamConvolution(beam_fwhm=13 * u.arcsec), Bin(edges=step.edges)]  # the 1-D equivalent of Render→PSF→Average
likelihood = Likelihood(GaussianFamily(), FixedCovarianceNoise(sigma_p) + GaussianProcessNoise(Matern32(...)))
ds = Dataset(observed=profile, model=shell, instrument=chain, likelihood=likelihood)
```
With inclination and position angle on `AzimuthalAverage` (elliptical
annuli; Peter: required, lower priority) the same code fits a disc's
surface-brightness profile from an image.

**(3) The ACA uv table, uv distance only (Route C).**
```python
vis = read_uvtable("ness_target_aca.npz")                                  # W8.12: u, v, Re, Im, weights, wavelength
chain = [HankelTransform()]                                                # ρ = |uv|; no Deproject for an azimuthally unresolved source
likelihood = Likelihood(ComplexGaussianFamily(), IndependentNoise())       # weights → uncertainties; independent
ds = Dataset(observed=vis, model=shell, instrument=chain, likelihood=likelihood)
```
The shell's parameters shared with (1) by name, so the map and the uv
table fit **jointly** in one `FittingProblem` — the composition the SED
case already does, and the one that makes the two NESS data sets answer
one question.

**(4) A VLTI star with a shell.**
```python
star = LimbDarkenedDisc(diameter=..., coefficient=..., flux=...) + ThinShell(radius=..., flux=...)
chain = [HankelTransform(), SquaredAmplitude()]                           # V² as measured (W7.4)
ds = Dataset(observed=read_oifits("target.fits")["vis2"], model=star, instrument=chain, likelihood=...)
```

**(5) A disc, non-parametric, frank-style — and the hybrid.**
```python
disc = LatentProfile(knots=np.geomspace(0.01, 2.0, 100) * u.arcsec, kernel=Matern32(amplitude=..., scale=...))
chain = [Deproject(inclination=..., position_angle=..., offset=(..., ...)), HankelTransform()]
ds = Dataset(observed=read_uvtable("dsharp_disc.npz"), model=disc, instrument=chain,
             likelihood=Likelihood(ComplexGaussianFamily(), IndependentNoise()))
# the hybrid: a parametric disc plus a latent correction whose posterior localises the misfit in r
disc = PowerLawDisc(gaps=2, ...) + LatentProfile(knots=..., kernel=Matern32(...))
```
`Deproject`'s inclination as `Derived("arccos(b_over_a)")` from a fitted axis
ratio, and a gap's width as `Derived("r_gap * f_gap")`, are W7.0's grammar
on this item's parameters.

## 6. What this changes in W8.2, and the order Peter's priorities set

- **Order of work inside the item**: the kind; `Render` and the
  `ProjectedShell` model (Route A on a synthetic NESS-like map — the first
  end-to-end fit); `HankelTransform` and the uv route on the same shell
  (Route C; the A–C agreement row); `AzimuthalAverage` with `covariance`,
  `BeamConvolution`, `Bin` and `FixedCovarianceNoise` (Route B; the
  covariance row); `Deproject` and the disc models (`GaussianRing`,
  `PowerLawDisc`, `Sersic`) and `LatentProfile`; elliptical annuli on
  `AzimuthalAverage` last. Sector averaging is not in the item.
- **The example on real data** becomes a NESS object if Peter can name one
  whose SCUBA-2 map and ACA uv table are public (or releasable under the
  `tests/data` pattern for a small cutout); the disc example (DSHARP or
  frank's small set) is the second figure, not the first.
- **`FixedCovarianceNoise`** is a §4.4 addition inside the item (or an S
  filler before it) with its own decision-log row; the memo's §4 is its
  justification.
- **The conformance rows** named in §2 and §4: A against C on a synthetic
  shell; `HankelTransform` of `UniformDisc` against `UniformDiscVisibilities`
  and of a Gaussian against `GaussianSourceVisibilities`; `Render` →
  `FourierSample` against `HankelTransform`; `AzimuthalAverage` of a
  rendered profile recovering it; the covariance row (coverage with the
  exact `Σ_P` against the naive diagonal); `LatentProfile` recovering an
  injected ring.

## 7. Questions for Peter, narrowed

1. **The NESS example**: which object, and are its SCUBA-2 map (with the
   variance array) and ACA uv table public or releasable? If none can be,
   the example runs on a synthetic shell at NESS-like noise and the real
   fit is yours to run.
2. **The ACA tables' form**: CASA measurement sets exported how — `uvtable`
   text/npz as frank/galario take, or UVFITS? W8.12 reads what you have.
3. **`ProjectedShell`'s density laws for the first cut**: a single power
   law between two radii, or also a broken power law and a shell with a
   Gaussian-in-r density? Each is a few lines once the projection exists.
4. **Priority of `LatentProfile`** against the parametric models: in the
   first cut (recommended — it is the frank case and the hybrid's
   correction) or a follow-on?
5. **`FixedCovarianceNoise`** as described — in W8.2, or as its own S
   filler landing first so the SED and line-flux cases can use a covariance
   table too?

## Appendix A — the probe

```python
# numpy only. 128x128 map, 2"/pix, Gaussian beam FWHM 13" (s = 2.76 pix), unit smoothed rms;
# annuli of 6" to 90". A is the averaging matrix; the exact covariance is A C A^T with
# C_pq = exp(-d_pq^2 / (4 s^2)); a 400-draw Monte Carlo (white noise convolved with the beam
# by FFT, normalised to unit rms) checks it. Full script: the orchestrator's scratchpad,
# radial_cov_probe.py, 2026-10-15; the numbers are §4's table.
#
# Second block: 5 000 visibilities with weights U(0.5, 2) binned in 11 |uv| bins — the bin
# covariance is diagonal by construction, sd = 1/sqrt(sum of weights) = 0.040-0.044.
#
# Third block: the trapezoid Hankel transform of a 2-D Gaussian (sigma 0.5", 2 000 r points)
# against exp(-2 pi^2 sigma^2 rho^2) at 0-300 k-lambda: max |error| 7.7e-7.
```
