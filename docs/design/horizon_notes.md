# Horizon notes — questions for Phase 5 at the earliest

Status: **notes, not a plan.** Raised by Peter on 2026-09-09 while Phase 3
was in flight, with the instruction that none of this is to be considered
before Phase 5, and possibly not until Phase 6 is complete. Recorded so the
thinking is not re-derived. Where a note has a consequence for work that
is being built *now*, that consequence is stated explicitly at the end and
is the only part of this document with any claim on the present.
(The observation-context reservation was confirmed by Peter 2026-09-09.)

**Folded into `DEVELOPMENT_PLAN.md` on 2026-09-10 (Peter's instruction:
"so the whole horizon is visible in one place")**: §1 and §2 as Phase 5
bullets (warping with its degrees-of-freedom guard and the sparsity prior;
the approximate-GP bullet widened to HSGP/EFGP with its two contract
questions); §3 and its follow-up as a Phase 5 bullet (joint noise over a
tuple of channels) and design horizon (h); §4–5 as a Phase 5 bullet and
design horizon (i); the embedding follow-up's ConvCNP experiment into the
deferred embedding study (§6); the kernel follow-up was already Phase 4
(W4.5). The plan is the authority from here; this file stays as the
reasoning behind each entry.

Each section: the question, what the frozen design already provides, what
it would take, and a first recommendation. Citations are to the design
documents and to the literature by name; nothing here has been measured.

## 1. More flexible misspecification models: warping and deep kernels

**Question.** Could amplitude warping (splines or small networks) or Deep
Kernel Learning give a misspecification model that handles individual
spectral lines or molecular bands without stacking several non-stationary
kernels, and without inflating the degrees of freedom until anything fits?

**What exists.** `likelihoods.md` fixes `Kernel` as a core, backend-neutral
object with a `KernelSpec` that the reference, torch and jax backends
translate to celerite2 terms (W2.3, W2.4/W2.5 slice 2). The default is
Matérn-class, chosen for structured residuals (plan §2, "Flexible-likelihood
kernel"). The O(N) path (`QuasisepGP`) requires an *ordered* 1D coordinate
axis (`REQUIRES_ORDERED_1D`), not an evenly spaced one. M2's
`strong_sharp` scenario is already an unmodelled line, and the flexible
likelihood stays calibrated on it with a stationary kernel — the question
is about *efficiency and localisation*, not about whether the stationary
model fails.

**The two warps, and why they are cheap here.**

- **Input warping** (Snoek et al. 2014, "warped GP"): a monotone map
  `u = w(x)` of the coordinate, with a stationary kernel in `u`. A monotone
  warp preserves ordering, so `QuasisepGP`'s precondition survives and the
  O(N) solve is untouched — the kernel simply sees `w(x)` instead of `x`.
  Non-stationarity in `x` (a band where the correlation length shortens) is
  a region where `w` is steep. Parameterised as a monotone spline (knot
  increments as positive parameters) or a small monotone network.
- **Amplitude warping**: `K'(x, x') = a(x) K(x, x') a(x')` with `a > 0`.
  For a quasiseparable `K` this is `D K D` with `D` diagonal: the
  quasiseparable rank is unchanged and celerite2's factorisation applies
  unchanged (scale the `U`/`V` generators row-wise). So both warps keep the
  O(N) exact solve, which is the property that made the flexible likelihood
  usable at 20 000 points in M2.

Both compose as a **kernel wrapper**, `WarpedKernel(base, input_warp=,
amplitude_warp=)`, with the warp parameters ordinary `Parameter`s carrying
priors — the spec-hash, provenance and lowering machinery needs no new
concept. On torch/jax the warps are differentiable for free, so NUTS over
warp knots is available.

**Deep Kernel Learning** (Wilson et al. 2016) is the same idea with a
network as the feature map; in one dimension a learned 1D→1D map *is* an
input warp, so DKL adds little over a spline there except a worse
degrees-of-freedom problem. It becomes interesting for 2–3D (IFU, images),
where the feature map can learn anisotropy — Phase 5 territory with §2.

**Degrees of freedom — the real risk.** A warp with many knots can make
anything fit, which is exactly the failure mode the flexible likelihood is
meant to expose, not absorb. Three controls, in the order to try them:
(i) few knots with hierarchical shrinkage on the increments towards the
identity warp (so the stationary kernel is the prior mean); (ii) the
diagnostics families B/C (whiteness, localisation) applied to the *warped*
residuals as the check that the warp has not eaten the signal; (iii) the
posterior-predictive and the SBC/coverage machinery (W3.6) to show the
warped model's calibration on M2-style injected misspecification. The
M2 study is the natural harness: extend it with a "many lines / one band"
scenario and compare stationary Matérn, warped Matérn and a sum of two
kernels on bias, calibration and localisation.

**Alternatives worth listing so they are not forgotten.** Spectral-mixture
kernels (Wilson & Adams 2013) are sums of damped oscillators — celerite
terms already are, so a sum of SHO terms *is* a spectral mixture in this
codebase; change-point kernels for a band boundary; and the sparse
"line-list" prior (a fixed dictionary of Gaussian lines with sparse
amplitudes) which is not a GP at all and would live as a model component.

**Recommendation.** A Phase 5 item "non-stationary flexible likelihood:
input and amplitude warping as kernel wrappers preserving
quasiseparability", validated by an M2 extension. DKL proper waits for the
multi-dimensional GP work.

## 2. Beyond quasiseparable and sparse GPs: HSGP, EFGP and relatives

**Question.** Plan §5 Phase 5 names SVGP, SKI and Vecchia. Hilbert-space
GPs, equispaced Fourier GPs and similar reduced-rank spectral methods give
large speed-ups without sparsity or inducing points and work in 2–3D.
Should they be considered?

**What exists.** `GPSolver` (`ampere.core.likelihood`) is the extension
point and already carries `EXACT`, `REQUIRES_ORDERED_1D`,
`REQUIRES_QUASISEPARABLE`, `IMPLEMENTED`, `DIFFERENTIABLE` and `BATCHABLE`
flags — an approximate solver was anticipated. `latent_transform` is the
solver's whitening for the latent-GP and sampling paths. Multi-axis
containers exist (IFU cube sketch). The conformance suite compares solvers
against `DenseGP` bit-tightly, which an approximate solver cannot meet.

**The candidates.**

- **HSGP** (Solin & Särkkä 2020; Riutort-Mayol et al. 2023): the kernel's
  spectral density evaluated at the Laplacian eigenvalues of a bounded box
  gives `K ≈ Φ diag(S) Φᵀ` with `m` basis functions; cost O(N m + m³) via
  Woodbury, tensor-product bases in 2–3D (m grows as m₁·m₂·m₃, so it is a
  low-dimension method). Needs a stationary kernel with a closed-form
  spectral density (Matérn: yes) and a boundary factor chosen from the
  data extent. Trivially differentiable; the natural fit for NUTS on
  torch/jax and for the reduced-rank latent path (`m` latent variables
  rather than `N`).
- **EFGP** (Greengard, Rachh & Barnett 2023): equispaced Fourier features
  with the Toeplitz structure of the resulting normal equations exploited
  by FFT; O(N + m log m) with `m` set by the kernel's spectral decay;
  demonstrated to ~10⁸ points in 1–3D. Also stationary-kernel-only; the
  reference implementation is MATLAB/Python.
- **Random Fourier features / sparse spectrum**: the Monte Carlo cousin;
  cheaper to implement, noisier, rarely the right first choice.
- **SKI/KISS-GP**, **Vecchia**, **SVGP**: already named in the plan;
  Vecchia is the one that is *not* rank-limited and handles rough
  (Matérn-1/2, 3/2) processes in 2–3D well, where spectral methods need a
  large `m`.

**Where they fit.** Each is a `GPSolver` with `EXACT = False`, its own
`provenance_config()` (the approximation parameters — `m`, the box factor
— recorded, never hashed, per fold-in 10), a `latent_transform` of reduced
dimension, and `conditional_loo` either exact-in-the-approximation or
refused by name. Two contract questions to settle when the phase opens:
(a) the conformance tolerance class for `EXACT = False` solvers — compare
to `DenseGP` with an approximation-aware tolerance that tightens with
`m`, and assert the convergence, rather than a fixed number; (b) whether
the latent size (fixed at composition, `inference.md` limitation 17.4) may
be `m` instead of `N` — it is fixed at composition either way, so the
contract permits it, but `simulate(observe=True)` and the latent-GP
likelihood path must agree on which whitening they use.

**Recommendation.** Widen the Phase 5 bullet to "SVGP, SKI, Vecchia,
HSGP, EFGP — chosen by measurement, as celerite2 was", with the rule that
the first 2–3D solver landed is the one that wins a benchmark on the IFU
sketch's cube at realistic N. HSGP is the cheapest to land and the most
useful for NUTS; EFGP is the one to reach for at image scale.

## 3. Multi-task GPs for polarimetry and other vector observables

**Question.** Stokes components have physical cross-talk (instrumental
polarisation, Q/U rotating together). Is a multi-task GP the right
misspecification model for them?

**What exists.** `results_schema.md` §15.2 rules vector-valued observables
are *separate channels* (one value array per container); the astrometric
time-series sketch (§1, §8) identifies exactly this gap and its extension
point, and `likelihoods.md` §15 records it as a deliberate limitation:
"a `NoiseModel` bound to a tuple of channels rather than one, with a
coregionalisation-style kernel — Phase 5 territory". `inference.md`
limitation 17.3 says the joint likelihood is a plain sum with the
`DatasetCollection.contributions` override as the hook.

**Assessment.** Yes, worth it, and polarimetry is the sharpest motivation
but not the only one: vector astrometry, simultaneous multi-band light
curves with a shared calibration systematic, and IFU multi-line maps all
want the same thing. The model of choice is the **linear model of
coregionalisation** (Bonilla et al. 2008; Álvarez et al. 2012):
`K = Σ_q B_q ⊗ k_q(x, x')` with `B_q` small (T×T for T outputs) positive
semi-definite. Two properties matter here:

- When the outputs share one coordinate grid (Stokes I/Q/U at the same
  wavelengths — the common case), the single-`q` intrinsic model
  `B ⊗ K_x` is Kronecker-structured: diagonalise `B = V Λ Vᵀ` (T×T,
  trivial), rotate the outputs by `Vᵀ`, and the T decoupled problems are
  each a scalar GP in `x` with the *same* `K_x`. With `QuasisepGP` that is
  T O(N) solves — **exact, and it preserves the O(N) path**. The general
  LMC (several `q`) or mismatched grids loses the Kronecker trick and
  needs the dense or a reduced-rank solver from §2.
- For Q/U the physically motivated `B` is a rotation-structured one
  (leakage mixes Q and U by an angle), so `B` can be parameterised with
  two or three parameters rather than a free T(T+1)/2, which keeps the
  degrees of freedom honest.

**What it needs.** The `NoiseModel`-over-a-tuple-of-channels extension the
limitation names, which implies a *joint likelihood over several
datasets* — the first real use of `contributions` as something other than
a sum, so it touches `inference.md` §4 and the results emission (a joint
`log_likelihood` term spanning datasets; `results.md`'s per-dataset
decomposition needs a "joint" entry — the same vocabulary widening
`"mixed"` was). The diagnostics (whiteness, localisation) generalise per
rotated output.

**Recommendation.** A Phase 5 item alongside hierarchical inference,
scoped to the shared-grid intrinsic model first (exact, O(N)), with
polarimetry as the worked modality and the astrometric sketch as the
second test.

## 4. Amortising SBI over different observed noise

**Question.** Can the SBI layer being built now amortise over observations
whose noise differs — two SEDs from different instruments, two spectra
with different error bars — as astronomy needs?

**What exists, and the gap.** `simulate(observe=True)` draws noise
through `LikelihoodFamily.sample` using the *observed container's own
uncertainties* (`inference.md` §13), so every simulation in a training set
today shares one noise realisation's σ-pattern: the one attached to the
dataset. A posterior trained on it is amortised over θ, not over σ. Noise
*parameters* (a fitted `scale` or `jitter`) are drawn from the prior and
so are amortised, but the per-sample error bars are not.

**What it takes — and it is not hard.** Two things, both within Phase 3's
own vocabulary:

1. **Noise as a simulated input.** `simulate_many` needs a per-draw
   *observation context* — at minimum a per-dataset uncertainty array
   drawn from a *noise-realisation prior* (a distribution over σ-patterns:
   scaled copies of the observed one, a library of real error arrays from
   the archive, a parametric S/N model) — passed to `sample` in place of
   the container's σ. The `Simulation` records which context it used.
2. **The embedding must see σ.** The W3.3 encoding lists the uncertainty
   column as optional (`[coordinates…, value(s), uncertainty?, mask]`);
   for amortisation over noise it is *mandatory* — a network that cannot
   see the error bars cannot condition on them. With it, a set-based
   embedding conditions the posterior on the noise pattern of the
   observation at hand, which is exactly the amortisation asked for.

At inference time the observed dataset's actual σ enters the same column,
and the trained posterior is valid provided the training context prior
covered it — which SBC/coverage (W3.6) can check per observation. This is
the standard "noise-aware SBI" recipe in the astronomy SBI literature; the
codebase lacks only the context hook and the mandatory column.

## 5. Amortising over changes in sampling

**Question.** Can the current approach amortise over different sampling
(different wavelength grids, different resolutions), and does that belong
to the embedding networks?

**What exists, and the gap.** The coordinate–value–mask encoding (W3.3) is
*designed* for this: a set-based or attention embedding over
`(coordinate, value, σ, mask)` rows is invariant to how many rows there
are and where they sit, so it can consume any sampling — the plan's Phase
3 text says exactly that ("the encoding is what makes amortisation across
differently-sampled datasets possible"). Two things stand in the way:

1. **The simulator's output grid is fixed at composition.** The
   `Instrument` chain resamples the model to the observed container's
   coordinates once, at composition (`architecture.md` §4.3 requirements
   negotiation). Every simulation therefore lands on *one* grid. To
   amortise over sampling, the grid (and the resolution, the LSF width,
   the filter set) must be a per-draw input — the same **observation
   context** as §4, one level up: not just σ but the instrument's own
   settings. `FittingProblem` would compose per-context instruments, or
   the context would carry a grid the chain re-negotiates per chunk;
   `chunk_size` grouping by context keeps that cheap.
2. **Resolution is not just sampling.** Two spectra at different R differ
   in the *forward model* (the convolution), not only the grid; the
   context must carry the instrument parameters and the embedding must see
   them (a per-row or per-set feature: R, filter identity, exposure). This
   is what turns "amortised over sampling" into "amortised over
   instruments".

So: yes, the embedding is where the *invariance* lives, but the
*coverage* comes from the simulator — the training set has to contain the
variety the posterior is later asked about.

## Consequence for Phase 3 (the only present-tense claim here)

Points 4 and 5 share one hook — a per-draw **observation context**
(σ-pattern, grid, instrument settings) as an input to `simulate_many`,
recorded on the `Simulation`, drawn from a context prior, and surfaced to
the embedding through the encoding. Reserving it now is cheap; adding it
after W3.2/W3.3 have fixed the call signatures is not. Three small
consequences are recorded in `WORK_ITEMS.md`'s Phase 3 section rather than
as changes to the approved item bodies: the W3.3 encoding carries the
uncertainty column by default; W3.1 slice 2 and W3.2 leave a `context=`
slot in `simulate_many`/`SBIEngine` (accepting `None` today) and record it
in provenance; the full context machinery is a Phase 3 follow-on item
drafted when W3.3 lands, or a Phase 5 item if the budget is spent.

## Follow-ups from the second exchange (2026-09-09)

### 4/5. Embedding architectures suited to functional and process data

The W3.3 encoding gives any embedding rows of `(coordinates, value, σ,
mask, dataset id)`. What the network does with them decides how much
amortisation is real. In rough order of how well they fit astronomical
data, and all reachable through W3.2's `embedding=` slot (sbi 0.27 ships
the first two; the rest are user `nn.Module`s):

- **DeepSets / permutation-invariant pooling** (sbi's
  `PermutationInvariantEmbedding`): the baseline. Invariant to row order
  and count, so it handles irregular sampling and missing data; but the
  per-row network sees each point alone, so it learns *local* features and
  relies on pooling for the rest. Cheap, and the right first experiment.
- **Set transformers / attention over rows** (sbi's `TransformerEmbedding`
  with Fourier-feature positional encodings of the *coordinate*, not the
  index): attention lets rows interact, so line shapes, band edges and
  correlated residual structure are representable; masking is native. The
  natural default for spectra of varying length and grid.
- **Neural-process encoders**: the Conditional Neural Process family
  (Garnelo et al. 2018; Attentive NP, Kim et al. 2019; **Convolutional
  CNP**, Gordon et al. 2020; Transformer NP, Nguyen & Grover 2022) was
  built for exactly this — off-grid functional observations with varying
  sampling and missing data, treated as a *process*. The ConvCNP encoder
  ("SetConv" onto a fine internal grid, then a CNN) is translation
  equivariant along the coordinate, which is the right inductive bias for
  a spectrum whose features can appear anywhere, and it is explicitly
  discretisation-invariant. Any of these encoders is a valid embedding
  for NPE; the NP *decoder* is not needed because the density estimator
  plays that role.
- **Neural-operator branch networks** (DeepONet's branch net, Lu et al.
  2021; FNO encoders): coordinate-conditioned and discretisation-invariant
  by construction — the same ideas the plan's design horizon (c) names
  for coordinate-conditioned emulators, so an encoder written for
  emulation is reusable as an embedding and vice versa.
- **Noise awareness**: whichever architecture, feed `σ` (and `log σ`, and
  the whitened `y/σ`) as row features — a network cannot condition on
  error bars it cannot see — *and* make the σ-pattern a simulator input
  (the observation context of §4), because the encoder can only be
  invariant to what the training set varied.
- **Instrument context**: settings that are not per-row (resolution,
  filter identity, exposure) enter as a per-set conditioning vector, most
  simply through FiLM layers (Perez et al. 2018) or as extra tokens in an
  attention encoder. That is what turns "amortised over sampling" into
  "amortised over instruments".
- **Distributional data**: where the observation is itself a set of draws
  (a population of posterior samples, a Monte Carlo error budget), a
  DeepSets encoder over the draws is the standard treatment, and it
  connects to design horizon (b)'s population post-processing.

**What the W3.3 traps mean for every later embedding** (added
2026-09-09 after verifying sbi 0.27): sbi's interface is a module over a
single tensor `x`, so everything an encoder needs — coordinates, values,
error bars, masks, dataset identity, instrument context — has to travel
*inside* `x`. That makes the **packing format the contract** and the
network a free choice: `encoding.md` fixes named column groups and an
`unpack` helper, standardisation lives in the encoding with its
statistics in the layout (sbi's own z-scoring switched off), and each
embedding is a thin wrapper that unpacks and adapts — NaN rows for sbi's
set embedding, an attention mask for its transformer, a coordinate grid
for a ConvCNP encoder, a per-set context vector for FiLM conditioning.
The two features that matter most for amortisation are already in the
packing: whitened values with `log σ` (noise) and the coordinate with
Fourier features (sampling). What a future encoder adds is only its
inductive bias.

**Practical recommendation.** W3.3 lands set and transformer embeddings
because sbi ships them. A ConvCNP-style encoder is the first *experiment*
worth running once W3.6's calibration machinery exists to judge it — an
`examples/sbi/` script, not a contract change. Everything above consumes
the same encoding, so nothing in Phase 3 forecloses any of it.

### 3. Interferometry as a second motivation for multi-task noise

Peter's addition: interferometric data can also carry cross-talk between
flux, amplitude and phase, and — the stronger, more general point —
**misspecification is correlated across channels even when the data has
no intrinsic cross-talk**, because one wrong sky model produces coherent
errors in every derived quantity. That argument applies to any vector
observable derived from one model: Stokes components, RA/Dec, visibility
amplitude and phase (or real and imaginary parts), multi-band fluxes. It
makes a coregionalised misspecification model the *default* shape for
vector data rather than a special case for instrumental leakage. For
visibilities the coordinate is `(u, v)`, so the shared-grid Kronecker
trick of §3 still applies but the per-output solve is 2D — a dense solve
at thousands of baselines, or a §2 solver beyond that. Phase 4's
interferometry item (`docs/design/modalities/interferometry.md`) should
leave the channel pairing visible so Phase 5 can bind a joint noise model
to it.

### 1. Kernels: what exists, what is missing, and sparsity as the guard

**The kernel surface today is thinner than "already good".**
`ampere.core` has two kernels, `Matern32` and `SquaredExponential`; there
is **no kernel algebra** (no sum or product of kernels); and the celerite
translation table (`_QUASISEPARABLE_TERMS`) is private with one entry. A
user-defined `Kernel` subclass works on the dense path — the ABC is public
and `matrix`/`diagonal` are all a `DenseGP` needs — but it cannot reach
the O(N) path, because there is no public way to register its celerite
representation. So: users can define kernels; they cannot make them fast,
and they cannot combine them. Fringing in M2 is currently absorbed by a
stationary Matérn-3/2 (it stays calibrated), not modelled as a periodic
component.

**What a modest item would add** (**Phase 4**, ruled by Peter 2026-09-09:
it is extensibility work, letting users compose a more expressive noise
model; recorded in the plan's §5 Phase 4 bullets):
`Sum` and `Product` kernels with the obvious `KernelSpec` composition and
celerite translation (a sum of quasiseparable terms is quasiseparable; a
product is not in general and is refused on the O(N) path); a **damped
periodic / SHO term** — celerite's `SHOTerm` with `Q > 1/2`, or the
`RotationTerm` pair — which is exactly a quasi-periodic ripple and so the
right component for fringing; Matérn-1/2 and -5/2 as celerite-exact
siblings; a **public `register_quasiseparable_term`** so a user kernel
declares its own celerite builder; and, with sums in hand, a
spectral-mixture kernel is simply a sum of SHO terms with free
frequencies. Conformance rows: each new term against a dense evaluation
of its closed form, and the sum against the sum of matrices.

**Sparsity-inducing priors on the components, marginalised.** Peter's
point generalises to every noise model beyond `IndependentNoise`: once a
noise model is a *sum* of components (several kernels, plus warps), give
the component amplitudes a sparsity-inducing prior so the data switch off
what it does not need, and let the GP marginal likelihood integrate the
latent functions out — that marginalisation is already what the flexible
likelihood does; only the hyperparameters are sampled. The continuous
choice is the **regularised horseshoe** (Piironen & Vehtari 2017): a
half-Cauchy global scale and per-component local scales, expressible today
with `HierarchicalPrior` (`parameters.md` §9) without any new contract
surface; spike-and-slab is discrete and stays out (`lowering.md` §12 Q1
refuses discrete bijections). Two practical notes for when it is tried:
NUTS on horseshoe hierarchies wants the *non-centred* parameterisation, so
the lowering rules should offer it (a `lowering.md` §3 note, not a
contract change); and the same prior on warp-knot increments is the
cleanest answer to §1's degrees-of-freedom risk — the identity warp is the
sparse solution, and the data must pay to leave it.
