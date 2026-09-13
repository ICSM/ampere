# Modality sketch (c) — interferometric visibilities and closure phases

Status: **DRAFT for Peter's review.** W1.11 sketch (c). Validates
`DEVELOPMENT_PLAN.md` §4.2–§4.4 against the modality §2 names as the **Phase 4
proof modality**, and answers two claim-checks routed here by the merged
contracts:

- `transformations.md` §13.1 / §14 — what a Fourier-sampling `Transformation`
  can honestly publish as requirements, given that there is no pull-back
  through a chain;
- `likelihoods.md` §17 Q3/Q4 — whether closure phases need `VonMisesFamily`
  before the spec freeze, and how Rice should be parameterised.

**Landed at W4.1 (Phase 4).** This sketch is now a design record with an
implementation beside it, and where the two differ the implementation is the
truth. What shipped, and what moved:

- the kinds are `ampere.core.VisibilitySet` (amended to three axes) and
  `ampere.core.ClosurePhases` (five), both in core rather than in a modality
  namespace — D1's ruling of 2026-09-11;
- the steps are `ampere.backends.reference.interferometry`:
  `FourierSample`, `ClosurePhase`, `Amplitude`, `BandwidthSmearing`,
  `TimeSmearing`, with `UniformDisc`/`GaussianSource`/`Binary` emitting an
  `Image` and the same three emitting a `VisibilitySet` analytically;
- `VonMisesFamily` is implemented, with `sample()` (§5's verdict, plus the
  sampling-form principle added at W3.14); `RiceFamily` is still declared
  only, because its data type is polarimetry (`likelihoods.md` §3);
- §2's "`VisibilitySet` needs nothing added" is **superseded**: a chromatic
  sky error is sharp in wavelength and smooth in (u, v), and a kernel sees a
  container's axes only, so the kind gained a `spectral_axis`
  (`phase4_placement_memo.md` §3.6, ruled 2026-09-11). One consequence the
  sketch could not foresee: the axes now carry mixed units, so §7's GP over a
  real-valued `VisibilitySet` is refused until kernels gain an `axes`
  selector (W4.5);
- §4's `max_step = 1/(2 u_max)` shipped, with an `oversampling` multiplier
  beside it. The Nyquist step is the coarsest grid at which the *array's* own
  information is representable; it is not a statement about the accuracy of
  the quadrature, and a source with real power beyond the longest baseline
  has that power folded back by a sum at that step.

**Landed at W4.3 (Phase 4).** The native twins and the modality under every
engine, which is what the sketch was a proof *of*:

- the steps and the six source models exist on both modern backends,
  `ampere.backends.torch.interferometry` and
  `ampere.backends.jax.interferometry`, under the same class names. Both follow
  the **inheriting** twin pattern — each class derives from its
  `ampere.backends.reference` counterpart and overrides only the four capability
  flags and the arithmetic — because in this modality the declaration is the
  large and dangerous half: `FourierSample` publishes the field of view and the
  Nyquist `max_step` over the *expanded* coverage, and its `configure_from`
  fixes the flat sub-sample layout two other steps read back. Two backends that
  wrote that twice would eventually negotiate slightly different pixel scales
  for one declaration, and §4's own warning says why that would not look like an
  error;
- the native model surface is spelled `native_flux` / `native_grid` rather than
  `flux` / `grid`. Not a preference: every source model declares a *parameter*
  called `flux`, and `Parameterised._check_free_name` refuses a parameter whose
  name shadows a class attribute, so `flux` is unavailable as a method name on
  exactly these classes. Both backends' `problem.py` accept either spelling;
- §7's GP over a `VisibilitySet` **runs** now: W4.5's `axes=` selector landed,
  W4.2 implemented the circular complex closed form, and W4.3 samples it under
  NUTS on both backends. The GP on the *closure phases* is still refused by
  name — a GP added to a wrapped observable is a latent-variable model, which is
  W5.1's;
- `VonMisesFamily` gained native bodies on torch and jax. It had been refused on
  both on the stated grounds that `ampere.core` declared but did not implement
  it — true until W4.1 — and without them NUTS cannot see the closure-phase half
  of the joint fit at all;
- **two honest downgrades**, declared on the classes that own them:
  `UniformDisc` declares `DIFFERENTIABLE = False` on both backends, because a
  hard-edged disc rendered onto a grid is piecewise constant in its diameter and
  automatic differentiation returns the amplitude term while silently omitting
  the larger boundary term; and `UniformDiscVisibilities` declares it too, for a
  library reason each backend states — `torch.special.bessel_j1` has no
  backward at all, and `jax.scipy.special.bessel_jn` is accurate only inside a
  measured argument window. A uniform disc is therefore not a gradient-based
  model in ampere; the Gaussian and the binary are.

Nothing here was implemented when it was written. Every snippet below was
executed against the merged `ampere.core` at commit `8c4e99d` (see §8 for what
was checked and how);
the sketches are **not** wired into `tests/core/test_spec_doctests.py`, which
names its four contract specs explicitly. Whether modality sketches should be
executed by the suite is W1.13's call — **decided at the freeze: they stay
unexecuted.** The sketches are design records whose claims were verified by
execution at review time (§8); the *contracts* are the executed documents,
and freezing 600 lines of illustrative composition code as doctests would
pin exactly the code Phase 4's real implementation supersedes.

---

## 1. The composition, in one picture

```
DustyDisc (Model)
  └── channel "sky" : Image        (x, y in mas; surface brightness)
         │
         ├── Instrument "gravity"     [FourierSample]              → VisibilitySet (complex)
         │        └── Likelihood(ComplexGaussianFamily(), IndependentNoise())
         │
         └── Instrument "gravity_t3"  [FourierSample, ClosurePhase] → ClosurePhases (angles)
                  └── Likelihood(<wrapped family>, IndependentNoise())
```

Two instruments bind the **same** channel by name, publish requirements on the
same two axes, and `negotiate` unions them into one image grid the model builds
once. That is exactly the pattern §4.3 was designed for, and it holds here
without modification.

The choice of `Image` as the model channel — rather than the model emitting
visibilities directly — is deliberate and is what makes the modality a genuine
test of the contracts: it puts a kind-changing, coordinate-changing step in the
middle of the chain, which is the case `transformations.md` limitation 13.1
warns about.

A model that computes visibilities analytically (a uniform disc, a binary point
source) skips the Fourier step and emits a `VisibilitySet` channel directly,
with `Instrument([], channel="vis", input_kind=VisibilitySet)` — the "no steps"
instrument `transformations.md` §5 provides. Both routes are supported and
neither is privileged.

## 2. The containers

*Superseded in part at W4.1: `VisibilitySet` gained a `spectral_axis`, and
`ClosurePhases` gained one too (five axes, not four). See the status note at
the top.*

`VisibilitySet` needs nothing added. It is complex, its `(u, v)` axes are
dimensionless (baselines in wavelengths) with `Order.ANY`, its `uncertainty` is
real and means the per-component standard deviation of a circular complex
Gaussian, and its layout is `Layout.POINTS`, so a correlated noise model can
consume it. All four of those were checked against the merged code.

Closure phases need a **new kind**, and defining one out of tree is three class
attributes — the extensibility `results_schema.md` §13 promises:

```python
class ClosurePhases(FunctionSamples):
    """Closure phase per triangle, in radians, indexed by two of its baselines."""

    AXES = (
        AxisSpec("u1", physical_types=("dimensionless",), order=Order.ANY),
        AxisSpec("v1", physical_types=("dimensionless",), order=Order.ANY),
        AxisSpec("u2", physical_types=("dimensionless",), order=Order.ANY),
        AxisSpec("v2", physical_types=("dimensionless",), order=Order.ANY),
    )
    LAYOUT = Layout.POINTS
```

Two baselines fix the triangle (the third is their difference), so four
dimensionless axes identify a closure phase the way a filter name identifies a
photometric point. The values carry `unit=u.rad`. This kind belongs in Phase 4's
interferometry package, not in `ampere.core`, exactly as `transformations.md`
§10 keeps the standard library out of core.

**Squared visibilities and amplitudes** need no new kind at all: a
`VisibilitySet` whose values are *real* is legal (`ALLOW_COMPLEX` permits
complex, it does not require it), and its `(u, v)` axes are still the right
index. §6 explains why that turns out to be a trap the contracts do not
currently catch (gap I-1).

## 3. Fourier sampling as a `Transformation`

```python
RAD2MAS = (1.0 * u.rad).to_value(u.mas)


class FourierSample(Transformation):
    ACCEPTS = (Image,)
    PRODUCES = VisibilitySet

    def __init__(self, u_pts, v_pts, *, field_of_view, **kwargs):
        super().__init__(**kwargs)
        self.register_buffer("u_pts", u_pts)          # wavelengths, from the data
        self.register_buffer("v_pts", v_pts)
        self._fov = float(field_of_view)              # mas, the fibre/primary beam

    def requirements(self):
        half = self._fov / 2.0
        umax = float(np.max(np.abs(self.buffers["u_pts"].value)))
        vmax = float(np.max(np.abs(self.buffers["v_pts"].value)))
        return (
            AxisRequirement("x", unit=u.mas, intervals=(-half, half),
                            max_step=RAD2MAS / (2 * umax), source="gravity"),
            AxisRequirement("y", unit=u.mas, intervals=(-half, half),
                            max_step=RAD2MAS / (2 * vmax), source="gravity"),
        )

    def apply(self, samples, values):
        ...                                           # direct DFT; see §4
```

The `(u, v)` coverage is a **buffer**, not a parameter — `architecture.md` §6's
question ("would you ever put a prior on it?") answers itself — and it is read
from the same file as the observed container, which matters (gap I-2).

The step carries no parameters. Calibration factors, a coherence-loss term or a
fitted phase offset are separate one-step transformations later in the chain,
which is how `transformations.md` §5 says instrument-level nuisance parameters
are expressed.

## 4. Claim-check 1 — what a Fourier-sampling step can publish

> `transformations.md` §13.1: "a chain whose third step needs something only
> expressible in the second step's output coordinates cannot say so."
> §14: "whether interferometry's `Image` → `VisibilitySet` step can publish
> anything useful given limitation 13.1."

### It can, and what it publishes is exact

Everything the Fourier step needs of the image is fixed by facts the
*instrument* owns at construction, and every one of them is expressible as
coverage-plus-density on the image's own `x` and `y` axes:

| What the step needs | Where the number comes from | How it is published |
|---|---|---|
| pixel scale Δx ≤ 1/(2·max\|u\|) | the array's longest projected baseline | `max_step` on `x` |
| pixel scale Δy ≤ 1/(2·max\|v\|) | same, on the other axis | `max_step` on `y` |
| field of view ≥ the beam | the fibre / primary beam, an instrument property | `intervals` on `x` and `y` |

Two properties make this exact rather than approximate, and both are worth
stating because neither is obvious:

1. **The requirement is separable.** The DFT kernel is
   `exp(-2πi(ux + vy)) = exp(-2πiux)·exp(-2πivy)`, so the sampling constraint
   factorises into an independent constraint per image axis. A requirement
   language that only speaks about one axis at a time is therefore not a
   limitation here — it is the natural shape of the constraint. (A `max_step`
   from the *radial* baseline length would be conservative on both axes and
   also expressible; the factorised form is simply tighter.)
2. **The step is the chain's first, so no pull-back is needed.** Its
   requirements *are* statements about the channel's coordinates, which is what
   `AxisRequirement` means. §13.1's limitation does not bite for the case §14
   asked about.

Executed against the merged code, with a GRAVITY-like 4-telescope array at
2.2 µm (baselines 22–46 Mλ) and a 60 mas field of view:

```python
>>> negotiate([Instrument([FourierSample(uu, vv, field_of_view=60.0)],
...                       channel="sky", label="gravity")])["sky"]
<ChannelRequirements 'sky' Image axes=['x', 'y'] from ['gravity']>
# x: <AxisRequirement 'x' mas [-30,30]@step<=2.46621>  -> 26 coordinates
# y: <AxisRequirement 'y' mas [-30,30]@step<=2.98541>  -> 22 coordinates
```

Twenty-six by twenty-two pixels is the grid the negotiation asks for, and it is
the right one: finer buys nothing the array can measure, coarser aliases. The
closure-phase instrument publishes the identical requirement (it wraps the same
Fourier step), and `negotiate` unions the two to the same grid — the model
evaluates **once** and feeds both instruments, which is the whole argument for
§4.3.

### Three findings that do need recording

**(a) A step *after* the Fourier step can publish nothing.** Bandwidth smearing
and time smearing are the real cases: both want more `(u, v)` samples per output
visibility. Executed:

```
CompositionError: instrument 'alma2' published a requirement on axis 'u', but
the Image on channel 'sky' is indexed by ['x', 'y']. Requirements name the
container's own axes.
```

That message is correct and helpful, and the failure is loud rather than silent,
which is the important half. But it means the two effects cannot be written as
independent, reusable steps — see gap I-3. Note that the resolution is **not**
`pull_back`: what a smearing step wants is not a different *image* grid, it is
more `(u, v)` points, and those are chosen by the Fourier step, not by the
model. A pull-back through a kind-changing step would be the wrong mechanism for
the right problem.

**(b) Regularity can be taken but never required, so Phase 4 owes a DFT path.**
`AxisRequirement.coordinates()` happens to return a uniformly spaced grid for a
`max_step` requirement (verified), so an FFT-gridding implementation will
usually get what it wants. It cannot *require* it: `architecture.md` §7 and
`results_schema.md` §14 both make regularity advertised-only, and `compile_for`
may return any grid satisfying the requirement. A `FourierSample` that assumes
an FFT is therefore incorrect. The rule is: branch on `samples.x.regular and
samples.y.regular` for the gridded fast path, and always ship the direct DFT.
This is the contract behaving as designed, but it is an implementation
obligation that should be written down before somebody discovers it as a bug.

**(c) `compile_for` cannot refuse, and interferometry is where that is most
dangerous.** `transformations.md` limitation 13.9 and §15.3 leave open whether a
model may refuse a requirement. A radiative-transfer model with a fixed image
size that silently ignores `max_step` produces **aliased** visibilities: power
from outside the Nyquist limit folds back and looks exactly like real source
structure at the sampled baselines. That is not a bad fit that announces itself;
it is a plausible wrong answer. This modality is a concrete vote for §15.3, and
for the *loud* option rather than the report option — see gap I-4.

## 5. Claim-check 2, part one — closure phases and wrapped likelihoods

> `likelihoods.md` §16: "check whether closure phases need `VonMisesFamily`
> before the freeze, since it is currently declared-only."
> §17 Q4: "The natural mapping from a per-sample σ is `κ = 1/σ²`, exact only in
> the small-σ limit. … Closure-phase practice (W1.11) should decide."

### Closure phases cannot be fitted with any currently implemented family

Every implemented family computes an unwrapped `observed - predicted`.
`GaussianFamily` on phases is not merely imprecise near ±π, it is catastrophic,
and it is silent. Executed, with a true 2° error straddling the branch cut
(observed 179°, predicted −179°, σ = 5°):

```
GaussianFamily residual it actually uses (deg): [358.0, -2.0, -2.0]
GaussianFamily   log_prob = -2558.88
wrapped family   log_prob =     4.32
```

A 2° error is charged as a 358° error, a 5100-nat penalty, on a triangle that
fits perfectly. No sampler recovers from that; it simply avoids the region of
parameter space where the model phase is near ±π. So **the answer to "do closure
phases need a wrapped family" is yes, unconditionally** — this is not a
refinement, it is the difference between working and not.

### But that is an implementation, not an interface change

A wrapped family is six lines against the frozen `log_prob(predicted, observed,
noise)` signature, and needs nothing that is not already there:

```python
@register_family
class VonMisesFamily(LikelihoodFamily):
    NAME = "von_mises"

    def log_prob(self, predicted, observed, noise):
        kappa = 1.0 / _independent_sigma(noise, self.NAME) ** 2
        delta = np.angle(np.exp(1j * (observed - predicted)))     # wrap to (-pi, pi]
        return float(np.sum(kappa * np.cos(delta) - kappa
                            - np.log(2 * np.pi) - np.log(sp.i0e(kappa))))
```

Written that way it was composed, aligned and evaluated against the merged code
with no change to `ampere.core`. **Verdict on Q3: `VonMisesFamily` does not need
to be implemented before the spec freeze**, because implementing it needs no
contract amendment. It should be implemented with Phase 4, and `IMPLEMENTED`
flipped then. *(Amended W4.8: done — landed at W4.1, `IMPLEMENTED = True`, the
six-line sketch above essentially unchanged; see the class docstring.)* What
*does* need to land in the freeze is gap I-5 below — a
family-level composition-time hook — which is a genuine interface addition that
this family is the first to need.

### Ship von Mises, not a Gaussian on the wrapped residual

A Gaussian evaluated on the wrapped residual is the obvious cheap alternative,
and it is adequate at high signal-to-noise and wrong at low. It is not a
normalised density on the circle: its mass over (−π, π] is less than one, by an
amount that depends on σ, and σ varies per triangle in any real dataset, so the
error does not even cancel as a constant. Measured (log-density at zero
residual, per sample):

| σ | von Mises | unnormalised wrapped Gaussian | uniform on the circle |
|---|---|---|---|
| 5° | 1.5189 | 1.5199 | −1.8379 |
| 30° | −0.3135 | −0.2719 | −1.8379 |
| 90° | −1.4732 | −1.3705 | −1.8379 |
| 180° | −1.7391 | **−2.0637** | −1.8379 |

At σ = 5° the two agree to 0.001 nat and either would do. At σ = 180° the
wrapped Gaussian assigns *less* density to a perfect match than the uniform
distribution does, which is impossible for a density on the circle; von Mises
correctly tends to uniform. Since `likelihoods.md` §8 makes a point of the
distinction between "a statement about a chi-square term" and "a normalised
log-density", shipping the unnormalised form here would contradict the contract's
own reasoning for a few flops of saving. `scipy.special.i0e` makes the
normalisation cheap and numerically safe at large κ.

**Verdict on Q4: κ = 1/σ² per sample, derived from the container's own
uncertainties.** Three reasons, in order of weight:

1. Every closure-phase pipeline reports a per-triangle σ propagated from the
   constituent baselines' signal-to-noise, and those vary by an order of
   magnitude within one dataset. A single fitted κ throws that away and lets the
   best triangles be outvoted by the worst.
2. `κ = 1/σ²` is exact in the small-σ limit, which is the regime in which a
   closure phase carries information at all: at σ ≳ 1 rad the triangle is nearly
   uninformative under either convention, and both tend to the same uniform
   limit.
3. The case a fitted κ is really meant to serve — "the pipeline underestimates
   its closure-phase errors" — is already expressible without a new mechanism, as
   `IndependentNoise(scale=...)`, which is exactly the nuisance parameter it
   exists for. Composing `VonMisesFamily` with `IndependentNoise(scale=...)`
   gives `κ = 1/(s·σ)²` with `s` an ordinary fitted parameter, tied across
   datasets if wanted. No new declaration form is needed.

## 6. Claim-check 2, part two — Rice's parameterisation

> `likelihoods.md` §17 Q3: "Does the model predict the true amplitude, or the
> underlying complex value from which the amplitude is derived? The first is
> simpler; the second is what an interferometric model actually produces."

Both halves of that sentence are true, and the resolution is that they describe
different *stages*, not competing designs. **Recommendation: the model produces
the complex value; an `Amplitude` transformation in the instrument chain takes
the modulus; and `RiceFamily.log_prob` receives real arrays, with `predicted`
being the true amplitude ν and `noise.sigma` the per-component σ of the
underlying circular complex Gaussian.**

```
DiscModel → "sky": Image
    └── Instrument "amp"  [FourierSample, Amplitude]  → VisibilitySet (real |V|)
             └── Likelihood(RiceFamily(), IndependentNoise())
```

where `Amplitude` is `ACCEPTS = (VisibilitySet,)`, `PRODUCES = VisibilitySet`,
and `apply` returns `samples.with_values(np.abs(samples.values))` — a
kind-preserving, coordinate-preserving, parameter-free step whose mask
propagates for free.

Four reasons this is the right side of the split:

1. **The `log_prob(predicted, observed, noise)` signature survives contact with
   Rice.** The alternative — a family receiving a complex prediction and a real
   observation — breaks `Likelihood._check_kinds`, `check_alignment`'s complex
   check and `_retained`'s dtype handling all at once, for one family.
2. **`VisibilitySet`'s real `uncertainty` already means exactly what Rice's σ
   is.** `results_schema.md` §12 states it: the per-component standard deviation
   of a circular complex Gaussian. The amplitude of such a variable is Rice(ν, σ)
   with the same σ. No container change, no extra structure — which confirms the
   `RiceFamily` docstring's own claim.
3. **Discarding the phase becomes a visible step.** Fitting amplitudes rather
   than complex visibilities is a real scientific decision (you are throwing
   away the astrometric information because the phases are corrupted by the
   atmosphere). Making it a named transformation in the chain, rather than a
   convention buried in a family, is what §4.3 is for.
4. **The same step serves the plan's other Rice use case.** §4.4 lists Rice for
   *polarised intensity*, which arrives as a `Spectrum` or an `Image` of
   √(Q²+U²), not a `VisibilitySet`. If the modulus lived inside the family it
   would have to be re-solved per kind; as a transformation, `Amplitude` is
   written once and the family is kind-agnostic.

The cost, stated honestly: a fit to amplitudes alone still evaluates the full
complex visibility. That is free in practice (the DFT produces both) and is the
right structure anyway, because the common case is a **joint** fit to amplitudes
*and* closure phases, which wants the complex visibilities computed once and two
instruments reading the same channel — precisely the composition §1 shows.

One consequence worth flagging for Phase 4 rather than for the freeze: squared
visibilities V², the standard OIFITS observable, are the square of a Rice
variate. Whether the family should be `Rice` on |V| or a dedicated `V²`
formulation is a modelling choice with a Jacobian in it; it is expressible
either way (as a further `Square` transformation plus a family), and does not
change any interface.

## 7. The likelihood side, and one thing it cannot do

> **Status, W4.2 (2026-09-13): this section's recommendation is implemented**,
> and **W4.3** samples it: NUTS on a `complex_gaussian` + `GaussianProcessNoise`
> visibility problem runs on both modern backends, with the kernel bound through
> `GaussianProcessNoise.kernel_for(observed)` so that an `axes=("u", "v")`
> selector means the same thing natively as it does on the contract path.
> `complex_gaussian` + `GaussianProcessNoise` is `ANALYTIC`, with the circular
> GP as its meaning exactly as argued below; `GP_ANALYTIC_IMPLEMENTED` is
> `True`; and the implementation is the "one call on a stacked residual" this
> section predicted rather than the "two calls to
> `log_marginal_likelihood`" — `GPSolver`'s right-hand side may now carry `k`
> realisations sharing one covariance, so the factorisation and the
> log-determinant happen once (`likelihoods.md` §4 and §7). The prediction that
> `QuasisepGP` "will inherit it" is the one thing here that turned out **false**,
> and for a structural reason rather than a scheduling one: the O(N) path needs
> one ordered one-dimensional coordinate, and a visibility is a point of the
> `(u, v)` plane at a wavelength. It is refused by name. The section is kept as
> written below, as the record of the argument the ruling was made on.

The visibility likelihood is `ComplexGaussianFamily` plus `IndependentNoise`,
and it composes and evaluates cleanly on a `VisibilitySet` (verified). The
flexible likelihood, however, **was refused** for visibilities when this was
written:

```
LikelihoodError: the complex_gaussian family with a GaussianProcessNoise noise
model is a latent-variable model (see Marginalisation.LATENT), but
ComplexGaussianFamily does not implement the latent-conditional log_prob — its
CONSUMES_LATENT_GP is False.
```

That refusal is correct under the current declaration, and `likelihoods.md`
§17 Q6 explains why the declaration is conservative: "which complex GP" — equal
component covariances? a non-zero pseudo-covariance? — is a Phase 4 modelling
question the likelihood contract declined to answer by setting a flag, and it
routed the question here.

**Answer from the modality: declare `complex_gaussian` + `GaussianProcessNoise`
as `ANALYTIC`, with the model being one real kernel applied independently and
identically to the real and imaginary parts** (equal component covariances, zero
pseudo-covariance). Reasoning:

- It is the *only* complex GP consistent with what the container already
  encodes. `VisibilitySet.uncertainty` is a single real σ per sample precisely
  because the noise is circular; a GP with unequal component covariances or a
  non-zero pseudo-covariance is non-circular, and `results_schema.md` §16 already
  says such a model "must supply its own structure". Choosing the circular GP
  keeps the container and the noise model telling the same story.
- Under that model the likelihood is a product of two identical real Gaussian
  marginal likelihoods, one on `Re(residual)` and one on `Im(residual)`, sharing
  `K(θ) + diag(σ²)`. It closes in exactly the same algebra `DenseGP` already
  implements — two calls to `log_marginal_likelihood`, or one on a stacked
  residual — so the implementation cost is close to zero and `QuasisepGP` will
  inherit it.
- The alternative (staying `LATENT`) is not a conservative choice with a cost of
  nothing: because `CONSUMES_LATENT_GP` is `False`, the combination is currently
  *unusable*, so the flagship feature of the package is unavailable on the
  modality the plan chose to prove it with. §17 Q6 says as much.

The physical case for wanting it is strong: interferometric residuals are
correlated with baseline length (source structure the model gets wrong appears
at particular spatial frequencies, not independently per visibility), which is
the textbook situation the flexible likelihood exists for. Verified: a GP over
the `(u, v)` point set composes and evaluates today for a **real**-valued
`VisibilitySet` (both axes are dimensionless, so the mixed-unit check passes and
the isotropic length-scale is a baseline length in wavelengths, which is
physically the right thing). Only the complex case is blocked, and only by the
flag.

Two limitations inherited from `likelihoods.md` §15 apply and are acceptable
here: the kernel is isotropic in `(u, v)`, so an elongated source's anisotropic
residual correlation is not expressible (§15.2); and only one kernel per noise
model, so "short-baseline structure plus long-baseline structure" needs the
`SumKernel` extension point (§15.1).

*(W4.5 lifted the second of those — `Sum` and `Product` ship, and the chromatic
`Product(Matern32(axes=("u", "v")), Matern32(axes=("spectral_axis",)))` of
`phase4_placement_memo.md` §3.6 is one of W4.2's conformance rows. The first
stands: one isotropic length-scale per leaf, so an anisotropic `(u, v)`
correlation still has no expression. `likelihoods.md` §15.2 is where it lives.)*

## 8. What was verified, and how

Executed against `ampere.core` at `8c4e99d` in the pixi `dev` environment:

| Claim | Evidence |
|---|---|
| `VisibilitySet` axes are dimensionless with `Order.ANY`; `ALLOW_COMPLEX = True`; layout `POINTS` | read from `ampere/core/results_schema.py:1362-1434` |
| A user-defined `ClosurePhases` kind composes, and `Instrument([FourierSample, ClosurePhase])` builds and evaluates | ran the chain end to end; produced 4 closure phases from 6 visibilities |
| `propagate_mask` carries a masked baseline into every triangle that uses it | masking baseline 1 of 6 masked triangles 0 and 2 — the two containing it |
| Fourier requirements publish and negotiate | `26 × 22` mas grid from a 4-telescope array; two instruments unioned to one grid |
| A post-Fourier `(u, v)` requirement is refused | `CompositionError` quoted verbatim in §4(a) |
| `AxisRequirement.coordinates()` returns a uniform grid for `max_step` | checked `np.diff` constant |
| `GaussianFamily` on phases uses the unwrapped residual | 358° in place of 2°; log_prob −2558.88 vs 4.32 |
| A wrapped von Mises family needs no interface change | registered via `@register_family`, composed, aligned and evaluated |
| von Mises vs unnormalised wrapped Gaussian at four σ | table in §5 |
| `ComplexGaussianFamily` + `IndependentNoise` works on a complex `VisibilitySet` | `log_prob = 15.994799` on a 3-point set |
| `ComplexGaussianFamily` + `GaussianProcessNoise` is refused | error quoted in §7 |
| A GP over a real-valued `VisibilitySet`'s `(u, v)` composes and evaluates | `log_prob = 5.825558` |
| `check_alignment` does **not** compare value dtypes | see gap I-1 |
| `Axis.__eq__` uses `np.array_equal`, i.e. exact | `ampere/core/results_schema.py:392-399` |

## 9. Interface gaps

Five, none of them fatal. The composition closes: a model producing an `Image`,
two instruments negotiating one grid, a kind-changing Fourier step, a
user-defined closure-phase kind, mask propagation through a three-to-one step,
and two likelihoods — all of it runs today against merged `ampere.core` with no
change to the package. The gaps are one real defect, one flag, one missing hook,
and two documentation obligations.

### I-1 — `check_alignment` does not compare value dtypes, so a complex prediction can be fitted against real amplitudes

*Ruled 2026-09-02: **approved and landed** — `check_alignment` compares value
dtype kinds, with the proposed message.* *(Amended W4.8: exercised as
predicted — W4.2's own report confirms "`check_alignment`'s dtype check
(gap I-1) still catches an amplitude fit".)*

**Severity: defect.** This is the only kind in the schema for which both a real
and a complex value array are legal, so it is the only place this can happen —
and it happens silently. Verified:

```python
>>> lk = Likelihood(ComplexGaussianFamily(), IndependentNoise())
>>> lk.check_alignment(complex_prediction, amplitude_observation)   # passes
>>> lk.log_prob(complex_prediction, amplitude_observation)
-154.880201                     # fits |V| against V; no warning anywhere
```

`check_alignment` compares kind, shape, unit and axes, and checks that a complex
*observed* container is not handed to a real family — but never compares
`predicted.values.dtype.kind` with `observed.values.dtype.kind`.

**Proposed amendment** — `likelihoods.md` §13, and `Likelihood.check_alignment`:

> Add, after the unit check:
>
> ```python
> if predicted.values.dtype.kind != observed.values.dtype.kind:
>     raise LikelihoodError(
>         f"the predicted {type(predicted).__name__} holds "
>         f"{'complex' if predicted.values.dtype.kind == 'c' else 'real'} values but the "
>         f"observed one holds "
>         f"{'complex' if observed.values.dtype.kind == 'c' else 'real'} values. A likelihood "
>         f"compares like with like: if you mean to fit amplitudes, take the modulus in the "
>         f"instrument chain (W1.5) so both sides are real."
>     )
> ```
>
> and a row in §14: *`check_alignment` compares value dtypes as well as kinds* —
> *`VisibilitySet` is legal with real or complex values, so kind equality does
> not imply comparability; without this a fit to amplitudes silently accepts a
> complex prediction.*

### I-2 — the predicted and observed `(u, v)` axes must be bit-identical, and nothing says so

*Dispositioned at the freeze (W1.13): **landed as documentation** — the
rule is now stated in `transformations.md` §10 (every coordinate-reproducing
step takes the observed coordinates from the observed container, never
recomputes them), in `likelihoods.md` §16's W1.5 bullet, and in the
`Dataset` docstring on the caller side. No tolerance was added: the check
stays load-bearing and exact.* *(Amended W4.8: the concrete instances are
`FourierSample.from_observed` (W4.1) and `EpochSample.from_observed` (W4.9),
both reading the observed container's own coordinates rather than
recomputing them.)*

**Severity: documentation, with a real trap behind it.** `Axis.__eq__` uses
`np.array_equal`, so `check_alignment` requires the Fourier step's `(u, v)`
buffer to match the observed container's coordinates *exactly*. That is
satisfiable — and cheap — if the instrument takes its coverage from the dataset.
It is not satisfiable if `(u, v)` are recomputed from station coordinates,
hour angle and wavelength, which is what an observer would naturally write; the
last few bits will differ and `check_alignment` will raise a message about axes
differing that gives no hint why.

Note `transform.py` already defines `COORDINATE_RTOL` for requirement
arithmetic, so a tolerance would not be a new concept — but a tolerance on
`check_alignment` would weaken a check that is load-bearing for every other
modality. The cheaper and safer fix is the rule, not the tolerance.

**Proposed amendment** — `transformations.md` §10, in the Fourier-sampling row's
note, and `likelihoods.md` §16's W1.5 bullet:

> A transformation that reproduces the observed container's coordinates must
> take them *from* that container (or from the file it was read from), not
> recompute them: `check_alignment` compares axes with `np.array_equal`, and a
> recomputed coordinate differing in its last bits fails a check whose message
> is about axes rather than about arithmetic. This is the general form of the
> instruction `likelihoods.md` §16 already gives resampling steps
> ("negotiation should target the observed grid").

### I-3 — a step after a kind-changing step cannot publish requirements at all

*Ruled 2026-09-03: **approved and landed at the freeze** — the chain-internal
`configure_from(downstream)` is in the contract (`transformations.md` §5),
called once per step at `Instrument` construction; `pull_back` is not
adopted. The concrete smearing steps remain Phase 4's.* *(Amended W4.8: landed
at W4.1 — `BandwidthSmearing` and `TimeSmearing` in
`backends/reference/interferometry.py`, with `FourierSample` using
`configure_from` to learn how many extra `(u, v)` sub-samples they need,
exactly the mechanism this gap proposed.)*

**Severity: expressiveness; a real loss of reuse, with a workaround.** Bandwidth
smearing and time smearing are ordinary, reusable interferometric effects whose
needs are stated in `(u, v)`, and `negotiate` refuses them because the channel
is an `Image`. The workaround is to fold both into `FourierSample`, which works
and is what a Phase 4 implementation will do — at the cost of exactly the reuse
`transformations.md` §12 argues the factored design buys over 3ML's monolithic
plugin. A user who writes a smearing step for one array cannot compose it with
somebody else's Fourier step.

`transformations.md` limitation 13.1 names `pull_back` as the extension point.
**This modality's evidence says `pull_back` is the wrong mechanism**: the
smearing step does not want a different image grid, it wants more `(u, v)`
points, and those are the *Fourier step's* buffer, not the model's. A pull-back
would have to invert a Fourier transform to say something the model cannot act
on anyway.

**Proposed amendment** — `transformations.md` limitation 13.1 and §15.2, replace
the `pull_back` sketch with the honest statement:

> A step's `requirements()` are read as statements about the *channel's*
> coordinates. A step downstream of a kind-changing step therefore cannot
> publish at all, and `negotiate` says so. Two distinct needs hide behind that,
> and they want different mechanisms:
>
> - *"I need the model evaluated differently"* — the case `pull_back` would
>   serve, and no standard-library step has it, because every kind-changing step
>   in §10's table is the first in its chain.
> - *"I need the step before me configured differently"* (interferometric
>   bandwidth and time smearing, which want extra `(u, v)` samples the Fourier
>   step owns) — a **chain-internal** negotiation, between steps, which this
>   contract does not have and which `pull_back` would not provide. Until it
>   does, such a pair must be written as one step.
>
> The extension point for the second is an optional
> `Transformation.configure_from(downstream: Sequence[Transformation])` called
> once during `Instrument.__init__`, letting a step read its successors'
> declarations before the hot loop. Deferred to Phase 4, which is where the
> first real instance arrives.

### I-4 — `compile_for` cannot refuse, and here that produces a plausible wrong answer

*Ruled 2026-09-03: **approved and landed at the freeze** — `compile_for`
raises `CompositionError` when it cannot honour a requirement
(`transformations.md` §7), with `FittingProblem(lenient_compile=True)` the
explicit warning-and-proceed opt-out.* *(Amended W4.8: the concrete instance
is `FourierSample`'s Nyquist requirement, `max_step = 1/(2 s u_max)`, refused
by `compile_for` rather than aliased — landed at W4.1.)*

**Severity: this modality upgrades §15.3 from "would be nice" to "should be
decided".** `transformations.md` §15.3 asks whether a model should be able to
refuse a requirement. For spectra the failure mode is a bad fit. For
interferometry it is *aliasing*: a model that silently ignores
`max_step = 1/(2 u_max)` folds power from beyond the Nyquist limit back onto the
sampled baselines, where it is indistinguishable from real source structure.

**Proposed amendment** — `transformations.md` §15.3, ruling in favour of the
loud option:

> `compile_for` may raise `CompositionError` when it cannot honour a
> requirement, and a model that chooses to ignore requirements silently (the
> default `return self`) remains legal. The asymmetry is deliberate: ignoring
> negotiation entirely is a *declared* stance, visible in the model's code;
> accepting a requirement and quietly under-sampling it is not. W1.11's
> interferometry sketch is the case that makes the difference matter — an
> under-sampled image grid aliases, and aliased visibilities look like real
> source structure rather than like an error.

### I-5 — a `LikelihoodFamily` cannot declare a composition-time precondition on its data

*Ruled 2026-09-02: **approved and landed**, with one refinement found in
implementation: `check_observed` is called on the observed container only —
the unit-equality check already forces the two containers to agree on
everything a container carries, and value-range properties genuinely differ
between them (a Poisson rate is not an integer, and would fail the very check
its counts must pass).* *(Amended W4.8: the concrete instance predicted here
is `VonMisesFamily.check_observed`, landed at W4.1 — a family handed degrees
instead of radians is refused by name rather than scored as nonsense.)*

**Severity: interface addition; should land in the freeze.** `check_alignment`
calls `self._noise.check_compatible(family, observed)` — the *noise model* gets
a composition-time hook and the family does not. A circular family needs one:
its values must be angles in radians, and feeding it degrees is silently
accepted today. Verified:

```python
>>> w = Likelihood(VonMisesFamily(), IndependentNoise())
>>> w.check_alignment(pred_in_degrees, obs_in_degrees)     # passes
>>> w.log_prob(pred_in_degrees, obs_in_degrees)
-5.509                                                     # nonsense
```

The same hook would serve `PoissonFamily` (integer counts — currently checked in
the hot loop, on every evaluation, which §13's compile-once/evaluate-many split
says is the wrong place) and `RiceFamily` (non-negative amplitudes).

**Proposed amendment** — `likelihoods.md` §3 and §16, plus `LikelihoodFamily`:

> A family may override
>
> ```python
> def check_observed(self, observed: FunctionSamples) -> None:
>     """Composition-time precondition on the data. Default: no constraint."""
> ```
>
> called by `Likelihood.check_alignment` on both the predicted and the observed
> container, alongside `NoiseModel.check_compatible`. It is the family's half of
> the same obligation, and it exists because some preconditions are properties of
> the *sampling distribution* rather than of the noise: a circular family needs
> angles in radians, `PoissonFamily` needs integer counts, `RiceFamily` needs
> non-negative amplitudes. Checks that belong here are the ones currently forced
> either into the hot loop (`PoissonFamily`'s integrality test, which runs on
> every evaluation) or nowhere at all (a von Mises family handed degrees).

## 10. Requirements on W1.7

*Dispositioned at the freeze (W1.13): **all four discharged** — item 1 by
I-2's landed rule (the `Dataset` docstring now states the caller
obligation); items 2 and 4 by W1.7 as merged (`inference.md` §8: once,
jointly; `LikelihoodError` → −inf with a recorded reason); item 3 by
W1.7's dataset-label refusal plus the 2026-09-03 instrument-label ruling
landed at W1.13.* *(Amended W4.8: all four are now exercised by a real
two-dataset composition — W4.1's visibilities-plus-closure-phases fit on one
`sky` channel, with `label="vis"`/`label="t3"` naming item 3's collision
explicitly, as `interferometry.rst` §4 walks through.)*

None of these invents W1.7's API; they are properties its `Dataset` /
`DatasetCollection` must have for this modality to compose.

1. **A `Dataset` must be able to hand its own coordinates to its instrument.**
   Gap I-2's rule only works if the `(u, v)` buffer in `FourierSample` and the
   observed `VisibilitySet` come from one source. W1.7 should either construct
   the instrument from the dataset, or state that the caller must, and say so in
   the `Dataset` docstring.
2. **Two datasets sharing one model channel is the normal case here, not an edge
   case.** Visibilities and closure phases are two datasets, two instruments, two
   likelihoods, one `sky` channel, one negotiated grid, one model evaluation per
   sample. `DatasetCollection` must call `negotiate` across *all* datasets'
   instruments once and `compile_for` once — `transformations.md` §14 already
   assigns it the "when are they called" question; this modality is why the
   answer must be "once, jointly", not "per dataset".
3. **Instrument labels must be checked for collisions when datasets are
   assembled.** `Instrument.label` defaults to the channel name
   (`transformations.md` §15.4), so two instruments on `sky` both default to
   `"sky"` and silently become one component in the caller's merge dictionary.
   W1.7 should raise on duplicate instrument labels within a collection.
4. **The failure-signalling layer converts `LikelihoodError` to −inf.** A
   non-positive-definite `K + diag(σ²)` is reachable during sampling
   (`likelihoods.md` §17 Q1); interferometry adds a second reachable failure —
   a model image that is identically zero makes the normalised visibility
   undefined. Both should arrive at the engine as −inf with a recorded reason,
   per §4.5.

## 11. Open questions for review

*Dispositioned at the freeze (W1.13):* Q1 — **ruled 2026-09-03 as this
sketch recommends** (`likelihoods.md` §17 Q6) and landed: the declaration
is `ANALYTIC` with the circular GP as the fixed meaning, composition
refusing until Phase 4 implements the closed form. Q2 — ruled the same day
with `results_schema.md` §17 Q5: yes eventually, `extra_coords` gains
`Axis` support in **Phase 4**, not at the freeze. Q3 — deferred to
**Phase 4** as the sketch says (a user-kind signature, not a §4 contract).
Q4 — closed by `results_schema.md` §17 Q4's ruling: `(x, y, spectral)`
confirmed, this sketch's own analysis (the spatial axes stay adjacent
under a Fourier consumer) being part of the evidence.

*(Amended W4.8, checked against the code.* Q1 **landed at W4.2**:
`GP_ANALYTIC_IMPLEMENTED = True`, the closed form §7 describes. Q2 **landed
differently, at W4.1**: `extra_coords` did not gain `Axis` support at all;
instead `VisibilitySet` gained a third first-class axis, `spectral_axis`,
because the chromatic case (`phase4_placement_memo.md` §3.6) needed a
coordinate a kernel's `axes=` selector could see, not a label — see the
preamble above and `interferometry.rst` §1's "chromatic amendment" section.
Q3 **decided at W4.1**: `ClosurePhases` is the four-axis form plus
`spectral_axis` (five axes total), with the canonical ordering fixed in the
kind's own docstring, resolving the sketch's "four axes vs. a `triangle`
label" question in favour of the axes. Q4 is undisturbed by Phase 4.*)*

1. **Should `complex_gaussian` + `GaussianProcessNoise` be declared `ANALYTIC`
   (§7)?** This sketch recommends yes, with the circular (equal-component,
   zero-pseudo-covariance) GP as the fixed meaning. It unblocks the flagship
   feature on the plan's own proof modality. It is a flag plus a
   `CONSUMES_LATENT_GP`-free implementation path, and it answers
   `likelihoods.md` §17 Q6.
2. **Should `extra_coords` gain units, for per-visibility frequency?**
   `results_schema.md` §17.5 and `likelihoods.md` §17.7 both route this here.
   Answer: a multi-frequency observation is best expressed as **one channel per
   spectral window** (which needs nothing new), and within a window the
   per-sample frequency is genuinely a quantity in Hz rather than a label. So
   yes, eventually — `extra_coords` becoming a mapping of `Axis` is
   `results_schema.md` §15.4's own named extension point, and Phase 4 is when it
   is first needed. It is not needed for the freeze.
3. **Is `ClosurePhases` indexed by two baselines the right signature?** Four
   dimensionless axes is one option; a single `triangle` label in the manner of
   `PhotometricPoints.filters` is another, and would make the identity explicit
   rather than derived. The four-axis form has the advantage that a GP over
   closure phases (which nobody has asked for) would have coordinates to work
   with. This is a Phase 4 decision, not a freeze decision.
4. **A `Cube` channel for a frequency-dependent sky.** The multi-frequency
   Fourier step wants `ACCEPTS = (Cube,)` and publishes on `x`, `y` *and*
   `spectral_axis` — all three axes exist and the requirement language covers it.
   Worth confirming that `Cube`'s `(x, y, spectral)` ordering (`results_schema.md`
   §17.4) is still the right choice when the consumer is a Fourier transform over
   the two spatial axes; it is, since those stay adjacent.
