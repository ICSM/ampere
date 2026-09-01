# The deliberately-awkward instrument — a stress test of the §4.3/§4.4 split

Status: **DRAFT for Peter's review.** W1.11's third deliverable, mandated by
`prior_art.md` Tension 3 and lesson 3M1, and carrying the X-ray claim-check
`likelihoods.md` §16 routes here.

`prior_art.md` Tension 3 states the bet plainly: 3ML puts the instrument
response *and* the likelihood inside one opaque per-instrument plugin, which
"buys trivial extensibility for genuinely weird instruments at the cost of zero
code reuse across instruments"; ampere factors them, which "buys reuse at the
risk that some real instrument's physics won't decompose cleanly into the two
pieces". Lesson 3M1 asks for the specific test: **an instrument whose
calibration uncertainty and noise covariance are physically coupled**, to check
that the factoring does not force such users back into a monolithic plugin.

The instrument chosen is an **X-ray CCD spectrometer** (Chandra/ACIS,
XMM/EPIC, XRISM/Resolve — the family, not a particular one). It is the right
choice for three reasons: it carries the RMF/ARF claim-check the dispatch
requires; its noise is Poisson, so nothing about it fits the
Gaussian-residual assumptions the rest of the contracts are shaped around; and
its effective-area calibration systematic is precisely 3M1's coupled case, with
a large literature and no ambiguity about the physics.

**Verdict, stated first.** The split holds for the modality's *core* — response
folding, exposure, grouping and Poisson counts all decompose exactly as §4.3/§4.4
intend, and the flexible likelihood is available on top. It **breaks in one
specific and general place**: a `NoiseModel` never sees the model prediction, so
any noise whose magnitude depends on the prediction — the coupled calibration
case 3M1 asks about, and also the far more ordinary "add 10% model uncertainty"
— can only be written by collapsing the noise model into the family, which is
the monolithic outcome Tension 3 is about. That is gap X-1, it is a keyword
argument, and it should land before the freeze.

Every claim was executed against the merged `ampere.core` at commit `8c4e99d`;
§7 lists what was run.

---

## 1. The composition

```
PowerLaw (Model)
  └── channel "xray" : Spectrum        (spectral_axis in keV; photons/s/cm²/keV)
         │
         └── Instrument "acis"  [ResponseMatrix]        → Spectrum (counts per channel)
                  └── Likelihood(PoissonFamily(), IndependentNoise())
                          observed: Spectrum (channel energies in keV; integer counts)
```

That is the whole thing, and it works. `Spectrum`'s `spectral_axis` accepts the
`energy` physical type already (`AxisSpec(..., physical_types=("length",
"frequency", "energy"))`), so `results_schema.md` §16's question — *"whether
`Spectrum` with an energy axis plus a response matrix is expressible"* — is
answered yes, with nothing added.

## 2. Claim-check — RMF/ARF as a matrix multiply, and where exposure lives

> `likelihoods.md` §16: "The X-ray sketch should check that `PoissonFamily` plus
> a response matrix in the instrument chain is expressible without a per-sample
> exposure concept this contract lacks."

### The physics, and the one matrix

Expected counts in detector channel *i* are

```
mu_i  =  T · Σ_j  R(i, E_j) · A(E_j) · S(E_j) · ΔE_j
```

with `S` the model's photon flux density, `A` the ancillary response (effective
area, cm²), `R` the redistribution matrix (unit-normalised over channels), `ΔE`
the true-energy bin widths and `T` the exposure in seconds. Every factor except
`S` is a constant of the observation, so **all of them fold into one buffer**:

```python
class ResponseMatrix(Transformation):
    ACCEPTS = (Spectrum,)
    PRODUCES = Spectrum

    def __init__(self, channels, matrix, **kwargs):
        super().__init__(**kwargs)
        self.register_buffer("channels", channels, unit=u.keV)
        self.register_buffer("matrix", matrix)          # T · R · A · dE, precomputed once

    def requirements(self):
        return (AxisRequirement("spectral_axis", unit=u.keV,
                                points=self._true_energies, source="acis"),)

    def apply(self, samples, values):
        m = self.buffers["matrix"].value
        return Spectrum(self.buffers["channels"].value * u.keV,
                        m @ samples.values * u.ct,
                        mask=propagate_mask(samples, m))
```

Run against merged core with an 80-bin true-energy grid, a 40-channel response
and T = 12 000 s: the instrument composes, publishes its requirement as
`<AxisRequirement 'spectral_axis' keV 80 point(s)>`, negotiates, and produces a
40-sample `Spectrum` in `u.ct` totalling ~10 190 counts, which
`Likelihood(PoissonFamily(), IndependentNoise())` accepts and evaluates.

Three details that make it work and are worth naming, because each is a contract
decision paying off rather than an accident:

- **`points` is the right requirement, not `intervals` + `max_step`.** An RMF is
  tabulated on a fixed true-energy grid; asking the model for anything else
  would need re-interpolating the response every evaluation, which is issue
  #12's bug. `AxisRequirement.points` says "these coordinates exactly", and
  `transformations.md` §10's table already prescribes it for this row.
- **The value unit changes, so `with_values` cannot be used.** `with_values`
  refuses a unit change by design (`results_schema.md` §15.6), so the step
  constructs a fresh `Spectrum` — and therefore must propagate the mask itself.
  `Transformation.__call__` enforces that (`transformations.md` §6), and
  `propagate_mask(samples, m)` with the response matrix as the influence matrix
  is exactly right: a channel is masked if any true-energy bin that feeds it is.
- **`PoissonFamily.REQUIRES_UNCERTAINTY` is `False`**, so the counts container
  carries no `uncertainty` and none is invented. A count-based family defines
  its own dispersion; that is the declaration doing its job.

### Where the exposure lives: in the matrix, and the contracts already say so

**There is no missing concept.** The exposure is a constant of the observation,
so it belongs where every other constant of the observation belongs — in the
`Transformation`'s buffer. `likelihoods.md` is already unambiguous about the
direction, in `PoissonFamily.log_prob`'s own error message:

```
LikelihoodError: the poisson family needs non-negative integer counts, but the
observed values are not integral. Counts are counts; if the data are rates,
multiply by the exposure in the instrument chain (W1.5) rather than here.
```

That integrality check is what makes the design *force* the right answer rather
than merely permit it. A user who tries to work in count rates — the natural
instinct — is refused at the first evaluation with a message naming the fix.
Verified: a container of `counts × 0.93` raises exactly that.

The three variants all compose without a new concept:

| Situation | Where it goes |
|---|---|
| One observation, known exposure | folded into the response buffer |
| Several observations, different exposures and responses | several `Dataset`s, each with its own instrument and its own matrix |
| A *fitted* dead-time or livetime fraction | a one-step `Transformation` with an ordinary `Parameter` — verified: `Instrument([ResponseMatrix, ExposureScale])` gives `free_names == ('exposure_scale.livetime',)` |

**Answer to the claim-check: yes, RMF/ARF forward folding is a matrix multiply
on an energy-axis `Spectrum`, and no per-sample exposure concept is needed or
wanted.** The only thing missing is that neither spec says so, and someone will
reach for the likelihood; see gap X-4.

### One real constraint on the channel axis

`Spectrum` requires strictly increasing coordinates, so the channel axis must
use nominal channel energies that are strictly increasing (verified: a repeated
energy raises). Real RMFs satisfy this. A bare integer channel index would
*not* work, because the axis demands `length`, `frequency` or `energy` and
rejects `dimensionless`:

```
SchemaError: Spectrum's 'spectral_axis' coordinates are in  (physical type
'dimensionless'), which this axis does not accept.
```

That is arguably the wrong identity — a PHA channel is identified by its channel
number the way a photometric point is identified by its filter name — but it is
not a problem in practice and it buys a real benefit: a GP over the channel axis
(§4) has a physically meaningful length-scale in keV. Recorded as an
observation, not a gap.

## 3. Where the split holds under pressure

Three further awkward features of the modality, each of which the factoring
handles.

**Grouping.** X-ray spectra are almost always channel-grouped (adjacent channels
summed to a minimum count, or optimally binned). Grouping is another matrix
multiply and composes as a second step: `Instrument([ResponseMatrix, Grouping])`,
or by folding the grouping matrix into the response buffer. The grouped channel
coordinates must reproduce the observed container's exactly, because
`check_alignment` compares axes with `np.array_equal` — the same instruction
`likelihoods.md` §16 already gives resampling steps.

**Gain uncertainty.** If the detector's energy scale is uncertain, the response
matrix depends on a fitted parameter and cannot be a single precomputed buffer.
`apply` may do arbitrary array arithmetic over buffers and parameters, so this is
expressible — interpolate between tabulated gain settings inside `apply` — but it
does contradict `transformations.md` §10's flat claim that the response matrix
"is a buffer". One clause fixes the table (gap X-4).

**The flexible likelihood is available.** `PoissonFamily` +
`GaussianProcessNoise` composes, correctly declares `Marginalisation.LATENT`,
produces a whitened latent declaration (`Parameter('latent',
prior=norm(loc=0.0, scale=1.0), shape=(20,))`) and refuses a gradient-free
engine with a message naming the reason. That is the whole §4.4 machine working
on counts, and it is the case §4.4 singles out. Note that the kernel's
`length_scale` must then be declared in keV, enforced at composition:

```
LikelihoodError: the kernel's 'length_scale' is declared in um but the
Spectrum's coordinate axis is in keV.
```

## 4. Where it breaks — gap X-1, prediction-dependent noise

This is 3M1's coupled case, and it is the finding this document exists for.

### The physics

An X-ray effective-area calibration carries a systematic uncertainty of order
5–10 %, correlated smoothly across energy. It is **multiplicative on the
prediction**: the true expected counts are `mu_i · (1 + δ_i)` with `δ` a smooth,
correlated, zero-mean field. If `δ` is *sampled* — as a few spline knots or a GP
on the correction factor — it is an ordinary `Transformation` with ordinary
parameters, and the split holds perfectly; that is worth saying, because it is
the *right* way to do it and the contracts support it today.

If instead `δ` is **marginalised analytically**, which is what an observer who
wants a fast fit does, the marginal covariance of the counts is

```
Cov  =  diag(mu)  +  diag(mu) · C_delta · diag(mu)
```

— a covariance that depends on the **model prediction**, not on the data.

### The contract cannot express it

`NoiseModel` never sees the prediction. Verified by instrumenting a subclass:

```python
NoiseModel.sigma(self, observed, retain, values)
NoiseModel.noise_params(self, observed, retain, values, *, coordinates, latent, limits)
```

`Likelihood.log_prob` computes `predicted_values` and then calls
`self._noise.noise_params(observed, retain, resolved, ...)` — the prediction is
in scope, three lines away, and is not passed.

### The workaround is exactly the monolith Tension 3 warns about

It *is* expressible — as a **family**, because `log_prob(predicted, observed,
noise)` does receive the prediction. Verified end to end:

```python
@register_family
class FractionalModelError(LikelihoodFamily):
    NAME = "fractional"
    def __init__(self, f=0.1):
        self.register_parameter(Parameter("f", st.uniform(0.0, 0.5), value=f))
    def log_prob(self, predicted, observed, noise):
        s = np.sqrt(noise.sigma ** 2 + (noise.values["f"] * predicted) ** 2)
        r = observed - predicted
        return float(np.sum(-0.5 * ((r / s) ** 2 + np.log(2 * np.pi)) - np.log(s)))
```

That runs. And it is precisely the failure mode Tension 3 predicts: to express a
*noise model*, the user has had to reimplement the *sampling distribution*. The
consequences are not cosmetic:

- **It does not compose with the flexible likelihood.** Verified — composing it
  with `GaussianProcessNoise` is refused, because it is not `ANALYTIC_WITH_GP`
  and does not `CONSUME_LATENT_GP`. To get "10 % model error *and* a
  misspecification GP" the user would have to re-implement `GaussianFamily`'s
  delegation to `GPSolver.log_marginal_likelihood` inside their own family. That
  is the package's flagship feature, unreachable.
- **It is not reusable.** The same fractional-error idea is now welded to one
  family; a Student-t version is a second copy.
- **It is not the only case, or even the most common one.** "Add N % model
  uncertainty in quadrature" is a standard move in SED fitting — ampere's *own*
  primary use case — and today it requires writing a likelihood family.

### The amendment is a keyword argument

`Likelihood.log_prob` already holds `predicted_values` when it builds the noise
parameters. Passing it costs one argument, is backwards compatible for
out-of-tree noise models (keyword-only, defaulted `None`), and immediately makes
the case a ten-line `NoiseModel` that composes with every family and, being a
diagonal term, with `GaussianProcessNoise` too.

## 5. The other awkward case — gap X-2, a family with its own aligned array

X-ray fits with a background almost always use the Cash-with-background
("W") statistic: the source-region counts are Poisson in `mu_src + r·mu_bkg`,
the background-region counts are Poisson in `mu_bkg`, and the per-channel
background rates are profiled out analytically.

The Bayesian route — two datasets, a parametric background model, joint fit —
composes cleanly with the contracts as they stand and is the better science. But
the profiled statistic is what practitioners use, and it *should* be expressible
as a user-written family, since that is the extension surface §4.4 advertises.
It nearly is:

```python
@register_family
class WStat(LikelihoodFamily):
    NAME = "wstat"
    REQUIRES_UNCERTAINTY = False
    def __init__(self, background_counts, ratio):
        self.register_buffer("background", background_counts)   # aligned per channel
        self.register_buffer("ratio", ratio)
    def log_prob(self, predicted, observed, noise):
        ...
```

A family is `Parameterised`, so it can carry its own buffers — a genuinely nice
property, and it is what makes the "one method plus your own data" extension
work. But `Likelihood.log_prob` **excises** masked samples before calling the
family, and the family's buffer is full length. Verified:

```
unmasked : family sees 8 retained samples, buffer has 8   -> log_prob -24.6142
masked   : family sees 6 retained samples, buffer has 8   -> MISALIGNED
```

Masking two channels — an entirely ordinary thing to do — silently misaligns the
background array against the counts, or raises if the family is careful enough
to check its own shapes. Nothing in `NoiseParams` says which samples were kept:
its fields are `sigma, values, coordinates, kernel, solver, latent, limits`, it
is frozen, and only `limits` is excised on the family's behalf.

This is the same shape of problem as X-1 — the family/noise interface carries
what *ampere* knows about and nothing a user brings — and it has an equally small
fix (§6, gap X-2).

## 6. Interface gaps

The split **holds, with one break**. To be precise about the verdict Tension 3
asked for: of the four awkward features tested — response folding, exposure,
background, calibration systematics — three decompose exactly as §4.3/§4.4
intend, and the fourth (analytically-marginalised calibration error) does not,
for a reason that is one keyword argument deep rather than structural. Ampere is
not being pushed back towards 3ML's monolithic plugin by its architecture; it is
being pushed there by a signature.

### X-1 — a `NoiseModel` cannot see the prediction

**Severity: the stress test's finding. Should land in the freeze**, because it
changes an ABC signature that Phase 2's two backends will implement in lockstep.

**Proposed amendment** — `likelihoods.md` §5, §14 and §16, and `NoiseModel`:

> `NoiseModel.sigma` and `NoiseModel.noise_params` take the retained predicted
> values as a keyword-only argument:
>
> ```python
> def sigma(self, observed, retain, values, *, predicted=None) -> np.ndarray | None: ...
> def noise_params(self, observed, retain, values, *,
>                  predicted=None, coordinates=None, latent=None, limits=None) -> NoiseParams: ...
> ```
>
> `Likelihood.log_prob` passes the excised `predicted_values` it already holds.
> The argument is optional and defaults to `None` so that a noise model which
> does not need it — `IndependentNoise`, `GaussianProcessNoise` — is unchanged,
> and so out-of-tree subclasses written against the current signature keep
> working.
>
> Add to §14: *A noise model receives the prediction as well as the observation.*
> *A noise whose magnitude depends on the model — a fractional model uncertainty,*
> *an analytically marginalised multiplicative calibration systematic, a*
> *model-variance weighting of counts — is a `NoiseModel`, not a family. Without*
> *the argument the only way to express one is to re-implement the sampling*
> *distribution, which welds noise to family, cannot be reused, and cannot reach*
> *the GP path: exactly the monolithic collapse `prior_art.md` Tension 3 warns*
> *against.*
>
> Add to §15 as a limitation retained: the noise model sees the prediction but
> not the model's *parameters* beyond those the likelihood declares, so a noise
> term that depends on a physical parameter must be tied to a likelihood
> parameter or expressed as a transformation.

With the argument in place, the coupled case is:

```python
class FractionalModelNoise(NoiseModel):
    """sigma_eff² = (s·sigma_data)² + (f·predicted)²."""
    def __init__(self, f):
        self.register_parameter(_as_hyperparameter("f", f, None))
    def sigma(self, observed, retain, values, *, predicted=None):
        base = _observed_sigma(observed, retain, "FractionalModelNoise")
        f = self.context(values)["f"]
        return np.sqrt(base ** 2 + (f * predicted) ** 2)
```

which composes with `GaussianFamily`, `StudentTFamily` and, as a diagonal term
under the kernel, with the flexible GP.

### X-2 — a family cannot excise its own per-sample arrays

**Severity: correctness trap for any user-written family with aligned data.
Should land in the freeze** (it adds a `NoiseParams` field, which is frozen and
part of the family interface).

**Proposed amendment** — `likelihoods.md` §2, §8 and `NoiseParams`:

> `NoiseParams` gains
>
> ```python
> #: Boolean inclusion indicator over the *full* containers, so a family
> #: carrying its own aligned per-sample data can excise it the same way.
> retain: np.ndarray | None = None
> ```
>
> set by `Likelihood.log_prob` from the mask union it already computes. Add to
> §8: *Every array in `NoiseParams` covers the retained samples only. A family
> that carries per-sample data of its own — a background spectrum, an
> instrumental template, a per-sample weight from outside ampere — must excise it
> with `noise.retain`, because `Likelihood` cannot know about it. Failing to is
> a silent misalignment as soon as anything is masked.*

An alternative considered and rejected: a general `extra: Mapping[str,
np.ndarray]` on `NoiseParams`, excised alongside `limits`. It is more
convenient but forces the auxiliary data through the *noise model*, which is the
wrong owner when the data belong to the family. `retain` is smaller and puts the
excision where the knowledge is.

### X-3 — `PoissonFamily`'s integrality check runs in the hot loop

**Severity: minor; folds into interferometry gap I-5.** `PoissonFamily.log_prob`
runs `np.all(counts == np.round(counts))` on every evaluation. It is a
composition-time property of the data, and
`likelihoods.md` §13 is explicit that O(N) checks belong in `check_alignment`.
The interferometry sketch proposes a `LikelihoodFamily.check_observed` hook for
exactly this class of precondition (gap I-5); this is its second customer, and
`RiceFamily`'s non-negativity would be its third.

**Proposed amendment** — as I-5, plus: move `PoissonFamily`'s integrality test
into `check_observed`, keeping a cheap `rate > 0` guard in `log_prob` (which is
a property of the *prediction* and genuinely varies per evaluation).

### X-4 — the specs never say where the exposure goes, and overstate the response matrix

**Severity: documentation; two clauses.** The design works and nothing needs
changing in the code, but the reasoning lives only in an error message.

**Proposed amendment** — `transformations.md` §10, the response-matrix row and
its following note:

> | Response matrix (RMF/ARF) | `Spectrum` → `Spectrum` | usually none (the
> matrix is a buffer); a fitted gain or livetime is an ordinary parameter and
> then the matrix is rebuilt inside `apply` | `points` at the matrix's own
> tabulated energies |
>
> and, replacing the note's second sentence: *A response matrix is a matrix
> multiply, so X-ray forward folding is an ordinary `Transformation`.* **The
> exposure folds into that same matrix** — *`T · R · A · dE` is one buffer —
> because expected counts, not count rates, are what a Poisson likelihood
> compares against, and `PoissonFamily` enforces this by refusing non-integer
> observations. There is deliberately no per-sample exposure concept anywhere in
> the contracts: exposure is a constant of the observation and belongs with the
> observation's other constants. Several observations with different exposures
> are several `Dataset`s. A fitted livetime or dead-time fraction is a one-step
> `Transformation` with an ordinary `Parameter`. The value unit changes across
> this step (flux density in, counts out), so it constructs a fresh `Spectrum`
> rather than using `with_values`, and must therefore call `propagate_mask` with
> the response matrix as its influence matrix.*

and `likelihoods.md` §16's W1.11 bullet gains the answer: *confirmed — `Spectrum`
with a keV axis plus a response-matrix `Transformation` is expressible, and the
absence of a per-sample exposure concept is correct rather than a gap.*

## 7. What was verified, and how

Executed against `ampere.core` at `8c4e99d` in the pixi `dev` environment.

| Claim | Evidence |
|---|---|
| `Spectrum`'s `spectral_axis` accepts `energy` | `ampere/core/results_schema.py:1029-1034` |
| The RMF/ARF chain composes, negotiates and evaluates | 80 true bins → 40 channels; ~10 190 counts; `Likelihood(PoissonFamily(), IndependentNoise()).log_prob` = −159.522 |
| `points` requirement publishes and negotiates | `<AxisRequirement 'spectral_axis' keV 80 point(s)>` |
| The counts container needs no `uncertainty` | `PoissonFamily.REQUIRES_UNCERTAINTY` is `False`; `counts.uncertainty is None` |
| Non-integer counts are refused, with the exposure instruction | error quoted in §2 |
| A fitted livetime is an ordinary chain parameter | `free_names == ('exposure_scale.livetime',)` |
| A duplicated channel energy is refused | `SchemaError`, §2 |
| A dimensionless channel-index axis is refused | `SchemaError`, §2 |
| `PoissonFamily` + `GaussianProcessNoise` is `LATENT`, declares `(20,)` whitened latents, refuses emcee | all three run, §3 |
| `length_scale` unit is checked against the keV axis | error quoted in §3 |
| A `NoiseModel` receives `(observed, retain, values)` only | instrumented subclass, §4 |
| The prediction-dependent case works as a *family* and is refused with a GP | both run, §4 |
| A family's own aligned buffer misaligns under masking | 8 vs 6 retained, §5 |
| `NoiseParams` is frozen, with no user-supplied per-sample slot | `dataclasses.fields`, §5 |

## 8. Requirements on W1.7

1. **Two-region X-ray fits are a `DatasetCollection` case.** Source and
   background spectra are two `Dataset`s with a shared background model and
   different exposures. That is the Bayesian alternative to `WStat` and it needs
   nothing new — but W1.7 should confirm that two datasets may share a *model*
   channel while having different instruments, since here they share the
   background model but not the source model.
2. **`check_alignment` must be called with the observed container**, as
   `likelihoods.md` §16 already requires — here it is what catches a response
   matrix whose channel grid does not match the PHA file's.
3. **The failure-signalling layer sees a new reachable failure.** `PoissonFamily`
   raises when the model predicts a non-positive expected count, which happens
   for extreme parameter draws (a very hard spectral index, a zero
   normalisation). Like the non-positive-definite covariance of
   `likelihoods.md` §17 Q1, it should reach the engine as −inf with a recorded
   reason rather than as an exception.

## 9. Open questions for review

1. **Does X-1 belong in the freeze?** This sketch says yes: it is an ABC
   signature that both Phase 2 backends implement, and the conformance suite
   (W1.10) will encode it. Deferring it means either that "10 % model
   uncertainty" is unavailable through Phase 2, or that it arrives as a family
   and has to be un-picked later.
2. **Should ampere ship a `FractionalModelNoise` in the standard library?** Once
   X-1 lands it is ten lines, and it is the single most requested thing missing
   from legacy ampere's likelihood. It is not core's business (core ships no
   concrete transformations), so it would be a `backends/reference` /
   Phase 2 deliverable — but the *contract* should name it in
   `likelihoods.md` §5's list of noise models the way §10 of the transformations
   contract names the standard chain steps.
3. **Is `WStat`'s profiling acceptable at all?** It is a profile likelihood, not
   a marginal one, so it is not strictly Bayesian; ampere may reasonably decline
   to support it and point users at the two-dataset formulation. If so, X-2 is
   still worth fixing — the misalignment trap applies to any family carrying
   aligned data, not just this one — but the motivating example changes.
4. **Was a second awkward instrument needed?** Two others were considered and
   rejected as weaker tests: a heterodyne receiver with correlated gain drift
   (which is a `GaussianProcessNoise` on the residual, i.e. the contract's
   home ground), and an échelle spectrograph with per-order blaze (per-order
   calibration factors, which `parameters.md`'s array-valued parameters and
   `likelihoods.md` §17 Q5's per-channel scale question already cover). Neither
   strains the split as hard as X-ray, because neither has a noise term that
   depends on the prediction. If Peter wants a second, the strongest remaining
   candidate is a **microcalorimeter with pile-up**, where the response depends
   on the *count rate* and hence on the model — the same X-1 coupling, one level
   deeper, and worth checking only if X-1 is not fixed.
