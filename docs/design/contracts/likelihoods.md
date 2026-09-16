# Ampere v2 — Likelihood & NoiseModel Contract (W1.6)

Status: **frozen at `spec-v1.0`** (the tag created at the W1.13 merge, 2026-09; any later change to a §4 contract requires a decision-log entry in `DEVELOPMENT_PLAN.md` in the same PR — ground rule 9). Implements `DEVELOPMENT_PLAN.md` §4.4,
answers `results_schema.md` §16's and `parameters.md` §13's obligations on this
item, and closes issue #11's design question. Code:
`ampere/core/likelihood.py`, `ampere/core/exceptions.py`. Tests:
`tests/core/test_likelihood.py`.

Every worked example below is executed as a doctest by
`tests/core/test_spec_doctests.py`, so this document cannot drift from the
implementation without the suite going red. Examples share one namespace and
build on each other; read them in order.

Where this document and `DEVELOPMENT_PLAN.md` disagree, the plan wins and this
document is wrong — file it as a decision-log correction.

---

## 1. What this contract is for

This is where ampere says *how well a prediction matches an observation*, and
it is the home of the package's distinguishing feature: a flexible, GP-based
likelihood that is robust to model misspecification.

Legacy ampere fuses all of this into one place — a hardcoded RBF kernel, a
dense solve, one noise model, one family, and no way to say which parts you
want. This contract separates three things that vary independently:

| Piece | Question it answers |
|---|---|
| `LikelihoodFamily` | What is the sampling distribution of one datum? |
| `NoiseModel` | Is the noise independent, or correlated — and by what kernel? |
| `GPSolver` | How is the correlated algebra actually done? |

The separation is a deliberate bet. `prior_art.md` Tension 3 records the
alternative: 3ML puts instrument response *and* likelihood inside one opaque
per-instrument plugin, which buys trivial extensibility for genuinely weird
instruments at the price of no reuse at all between them. Ampere factors them,
which buys reuse at the risk that some real instrument does not decompose. That
risk is W1.11's to stress-test before the freeze; this document owns the split
once it is made.

It is **backend-neutral**: numpy, scipy, `astropy.units` and stdlib, with no
torch or jax import, lazily or otherwise (`architecture.md` §4 rule 1). Turning
these declarations into celerite2 terms, `torch.distributions` or numpyro sites
is W1.9's table and Phase 2's code.

*(**Amended W2.15**, 2026-09-08. The "no torch or jax import" half is
unchanged and is import-graph-tested. Two additions of fact. First, the
dependency list is now numpy, scipy, `astropy.units`, stdlib **and
celerite2** — which W2.3 made a base dependency of the whole distribution
rather than an extra, precisely so that the base install is not left with the
O(N³) problem, and which `ampere/core/likelihood.py` imports **lazily**,
inside `QuasisepGP`, so the module's import surface is exactly as described
above. Second, the `KernelSpec` → celerite2-term half of the last sentence
did **not** go to a backend: it landed here, in `ampere.core.QuasisepGP`,
at W2.3, over ampere's own exact rank-2 Matérn-3/2 representation. Only the
`torch.distributions` and numpyro halves are a backend's, and both landed at
W2.4/W2.5 through `ampere.core.lowering`'s registry.)*

### Setup for the examples

```pycon
>>> import numpy as np
>>> import scipy.stats as st
>>> import astropy.units as u
>>> from ampere.core import (
...     SHO, Censoring, ComplexGaussianFamily, DenseGP, GaussianFamily, GaussianProcessNoise,
...     HilbertSpaceGP, IndependentNoise, InducingPointGP, Likelihood, LikelihoodFamily,
...     LimitKind, Marginalisation, Matern12, Matern32, Matern52, NoiseParams, PhotometricPoints,
...     PoissonFamily, Product, QuasisepGP, RiceFamily, RotationTerm, SpectralMixture,
...     Spectrum, SquaredExponential, StationaryKernel, StudentTFamily, Sum, VisibilitySet,
...     WindowedSparseGP, family_named, list_families, quasiseparable_families,
...     register_family, register_quasiseparable_term,
... )
>>> from ampere.core.exceptions import LikelihoodError

```

## 2. The objects

| Object | Role |
|---|---|
| `LikelihoodFamily` | The sampling distribution of one datum. One method: `log_prob(predicted, observed, noise)` |
| `NoiseModel` | What the noise is. `Parameterised`, so its knobs are ordinary parameters |
| `IndependentNoise` | Uncorrelated per-sample uncertainties, optionally scaled and floored |
| `GaussianProcessNoise` | The flexible likelihood: a GP over the residuals |
| `Kernel`, `KernelSpec` | A covariance function, declared neutrally: family name plus `Parameter` hyperparameters |
| `Matern32` | The canonical flexible-likelihood kernel (plan §2) |
| `GPSolver` | How the GP algebra is done. A strategy: `DenseGP`, `QuasisepGP`, … |
| `Censoring`, `LimitKind` | Which observations are limits rather than detections (issue #11) |
| `Marginalisation` | `ANALYTIC` or `LATENT`: the declaration §4.4 requires |
| `LatentDeclaration` | The parameters a latent-GP formulation adds |
| `Likelihood` | The composition, and the object a `Dataset` (W1.7) holds |
| `LikelihoodError` | This contract's error; a `ContractError`, hence a `ValueError` |

## 3. Families: one method, and a registry

`DEVELOPMENT_PLAN.md` §4.4 fixes the family interface as
`log_prob(predicted, observed, noise_params)`, and `prior_art.md` lesson 3M3
says why it must stay that small: 3ML's `get_log_like()` is the right *shape*
of contract; its opacity is the part not to copy.

```pycon
>>> data = Spectrum(
...     [1.0, 2.0, 3.0, 4.0] * u.um,
...     [3.0, 2.5, 2.2, 1.4] * u.Jy,
...     uncertainty=[0.1, 0.1, 0.1, 0.1] * u.Jy,
... )
>>> model = data.with_values([3.05, 2.40, 2.25, 1.35])
>>> simple = Likelihood(GaussianFamily(), IndependentNoise())
>>> round(simple.log_prob(model, data), 6)
4.659586

```

`Likelihood(family)` alone means "independent noise", because the simple case
must stay simple:

```pycon
>>> Likelihood(GaussianFamily()).parameters.free_size
0

```

The registry is public and ordinary, so a user adds a family without touching
ampere (§4.4's requirement, verbatim):

```pycon
>>> @register_family
... class LaplaceFamily(LikelihoodFamily):
...     NAME = "laplace"
...     def log_prob(self, predicted, observed, noise):
...         z = (observed - predicted) / noise.sigma
...         return float(np.sum(st.laplace.logpdf(z) - np.log(noise.sigma)))
>>> family_named("laplace") is LaplaceFamily
True
>>> round(Likelihood(LaplaceFamily(), IndependentNoise()).log_prob(model, data), 6)
3.937752

```

### The optional generative half: `sample` (*Amended W3.14*)

**Ruled by Peter, 2026-09-02** (W1.7's ruling request R3, recorded in
`DEVELOPMENT_PLAN.md` §2): a family may also override

```
sample(predicted, noise, rng) -> np.ndarray
```

— one draw of the retained observed values from the same distribution
`log_prob` scores, given the same `NoiseParams`. It is what W1.7's
`simulate(observe=True)` delegates to. The **default refuses, specifically**:
it names the family and the override to provide, because a `log_prob` does
not determine an observation process and ampere will not guess one — a wrong
guess would silently train an SBI posterior on the wrong forward model.
`GaussianFamily` implements it for both noise models (including the GP draw
through the solver's own `latent_transform`, stabiliser included); a user
family supplies its own the same way `LaplaceFamily` above supplies
`log_prob`:

```pycon
>>> class SamplingLaplace(LaplaceFamily):
...     def sample(self, predicted, noise, rng):
...         return predicted + noise.sigma * rng.laplace(size=predicted.shape)

```

**Amended W3.14: which families implement it.** Peter ruled on 2026-09-09 that
further core `sample` families land "when we have a use case", and W3.6 supplied
one — a counting experiment could not `simulate(observe=True)`, so neither SBC
nor SBI on count data worked without the user subclassing `PoissonFamily`. The
refusal was written for families whose observation process is *genuinely
ambiguous*; three of the shipped families are not, each having one generative
form that its own `log_prob` already fixes:

| family | draw | `σ` means |
| --- | --- | --- |
| `poisson` | `counts ~ Poisson(rate)`, the rate being `predicted` (or `predicted · exp(f)` under a latent GP) | — the family defines its own dispersion |
| `student_t` | `x = μ + σ t_ν`, the family's own `ν` | the **scale**; the variance is `σ² ν/(ν−2)` |
| `complex_gaussian` | independent `Normal(0, σ²)` real and imaginary parts | the **per-component** standard deviation; the total variance is `2σ²` |

Each reads its parameters exactly as its `log_prob` reads them, which is the
whole content of "the generative counterpart": a fitted `scale` or `jitter` is
in the draw because it is in the `σ` the density uses, and the latent `f` a
Poisson draws at is the `f` the density scores at — `Dataset.draw_observation`
takes the whitened `z` out of θ and hands it to `noise_params`, exactly as
`Dataset.log_likelihood_of` does. `cauchy` and every declared-but-unimplemented
family keep the refusal, **word for word**; `student_t` refuses a correlated
noise model by name (a scale mixture of Gaussians does not commute with a GP
covariance) and `complex_gaussian` refuses one with the same
`GP_ANALYTIC_IMPLEMENTED` message its density refuses it with. `inference.md`
§13 carries the full statement, the reasoning and the doctests.

**Native twins (*Amended W3.1*, slice 2; *W3.14*).** Peter ruled on 2026-09-08
that every backend supports observation sampling natively, and W3.1 slice 2
landed it: a realisation may expose `sample_observations`, which draws a whole
chunk's retained values in the backend's own arithmetic. The ruling has a
ceiling, and it is this section's rule restated one level down: **a backend
samples exactly what `ampere.core` samples**. The default refusal is not a gap a
backend is free to fill — it is the contract declining to guess an observation
process — so `cauchy` refuses identically on torch and jax, and a family whose
`sample` a user has overridden is left to that override rather than replaced by
a native approximation of it. W3.14 moved the four families the core now samples
across that line: `gaussian`, `poisson`, `student_t` and `complex_gaussian` all
have twins on both backends, transcribing the same branches — GP draw and solver
stabiliser included, from the solver's own native `latent_transform`. The floor
is the other half of the same rule: a family the core samples but a backend has
no twin for falls back to the numpy path rather than refusing, because a backend
may be slower than the contract but never more restrictive than it.
`inference.md` §13's "sampling on a backend" carries the rest, including why the
numpy path stays the oracle, why the comparison is distributional, and why
torch's non-Gaussian twins are two-stage where jax's are not.

### The composition-time hook: `check_observed`

Also ruled 2026-09-02 (W1.11 gap I-5): a family may override
`check_observed(observed)` — the family's half of the obligation
`NoiseModel.check_compatible` already has. Some preconditions are properties
of the *sampling distribution* rather than of the noise: a circular family
needs angles in radians, `PoissonFamily` needs integer counts (its
integrality test now lives here rather than in the hot loop — gap X-3), a
Rice family needs non-negative amplitudes. `Likelihood.check_alignment`
calls it on the **observed** container only — the unit check has already
forced the two containers to agree on everything a container carries, and
value-range properties genuinely differ between them (a Poisson *rate* is
not an integer). Only retained samples are held to a precondition; a masked
sample carries zero information and cannot fail one.

```pycon
>>> counts = Spectrum([1.0, 2.0] * u.um, [4.5, 7.0])
>>> Likelihood(PoissonFamily(), IndependentNoise()).check_alignment(
...     counts.with_values([4.0, 7.0]), counts
... )
Traceback (most recent call last):
    ...
ampere.core.exceptions.LikelihoodError: the poisson family needs non-negative integer counts,
but the observed values are not integral. Counts are counts; if the data are rates, multiply
by the exposure in the instrument chain (W1.5) rather than here.

```

The families the plan's target scope needs are all *declared*, whether or not
they are implemented yet — a name in the registry is a commitment the spec
freeze can be reviewed against:

```pycon
>>> [name for name in list_families() if name != "laplace"]
['gaussian', 'student_t', 'cauchy', 'complex_gaussian', 'poisson', 'rice', 'von_mises']

```

A declared-but-unimplemented family refuses **composition**, not evaluation, so
the failure arrives when the problem is built rather than a thousand samples
into a run:

```pycon
>>> Likelihood(RiceFamily(), IndependentNoise())
Traceback (most recent call last):
    ...
ampere.core.exceptions.LikelihoodError: the rice family is declared but not implemented, so it cannot be composed into a Likelihood yet ...

```

Everything else a family declares is a class attribute, so declaring costs
nothing per evaluation and composition is checkable before a sampler starts:

```pycon
>>> GaussianFamily.ANALYTIC_WITH_GP, GaussianFamily.SUPPORTS_CENSORING
(True, True)
>>> PoissonFamily.REQUIRES_UNCERTAINTY, PoissonFamily.ANALYTIC_WITH_GP
(False, False)

```

A family may itself be `Parameterised` — Student-t's degrees of freedom are an
ordinary parameter, so they may be fitted, fixed or tied like anything else:

```pycon
>>> StudentTFamily(nu=st.loguniform(2.0, 50.0)).parameters.free_names
('nu',)
>>> StudentTFamily(nu=4.0).parameters["nu"].is_fixed
True

```

**The principle for a family's sampling form** (*ruled by Peter 2026-09-10*,
on W3.14's `cauchy` question; the rule §13 of `inference.md` applies).
`sample()` is added to a family when a data type needs that family's
likelihood — never speculatively, and never for a subset of the stack: once a
data type needs a likelihood, its sampling form is supported by **every
inference approach on every backend at once** (the numpy path as the oracle,
the native twins beside it, the conformance rows per backend), because a
likelihood that fits but cannot simulate breaks SBC and SBI for exactly the
users who brought the data type. So `cauchy` keeps its refusal until a data
type needs it, and polarimetry will bring Rice's sampling form together with
its likelihood, not after it.

## 4. The declaration that decides which engines can run: analytic or latent

This is the load-bearing paragraph of `DEVELOPMENT_PLAN.md` §4.4, and the
reason this contract exists as more than a function.

"The GP as a covariance matrix added to a Gaussian" only marginalises
analytically **because the family is Gaussian**. For Poisson, Student-t and the
rest, misspecification robustness needs a latent-GP formulation — `counts ~
Poisson(rate · exp(f))`, `f ~ GP` — marginalised numerically. So a
`Likelihood` declares which it is:

```pycon
>>> flexible = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(0.3, 2.0)))
>>> flexible.marginalisation
<Marginalisation.ANALYTIC: 'analytic'>
>>> counts = Spectrum([1.0, 2.0, 3.0, 4.0] * u.um, [4.0, 7.0, 2.0, 9.0])
>>> latent_poisson = Likelihood(PoissonFamily(), GaussianProcessNoise(Matern32(0.3, 2.0)))
>>> latent_poisson.marginalisation
<Marginalisation.LATENT: 'latent'>

```

Marginalisation is a property of the **combination**, not of either piece.
Poisson on its own is perfectly analytic; it is Poisson *plus a GP* that is not:

```pycon
>>> Likelihood(PoissonFamily(), IndependentNoise()).marginalisation
<Marginalisation.ANALYTIC: 'analytic'>
>>> Likelihood(StudentTFamily(), IndependentNoise()).marginalisation
<Marginalisation.ANALYTIC: 'analytic'>
>>> gp_noise = GaussianProcessNoise(Matern32(0.3, 2.0))
>>> StudentTFamily().marginalisation_with(gp_noise)
<Marginalisation.LATENT: 'latent'>

```

`marginalisation_with` is asked of the *family* there, and not of a composed
`Likelihood`, for a reason that is the sharpest edge in this contract.

### Declaring `LATENT` is not the same as implementing it

A family whose `log_prob` declares `LATENT` but never reads `noise.latent`
returns the **uncorrelated** likelihood. Every GP hyperparameter and every
latent value an engine sampled would then leave the log-probability untouched:
their posteriors would come back as their priors, and the physical parameters
would be biased exactly as they were under the rigid likelihood. The fit runs,
converges, and is wrong — which is the failure this whole contract exists to
make impossible, so it may not be reachable by composing two objects that each
look reasonable.

A family therefore opts in, and the default is refusal:

```pycon
>>> StudentTFamily.CONSUMES_LATENT_GP, PoissonFamily.CONSUMES_LATENT_GP
(False, True)
>>> Likelihood(StudentTFamily(), gp_noise)
Traceback (most recent call last):
    ...
ampere.core.exceptions.LikelihoodError: the student_t family with a GaussianProcessNoise noise model is a latent-variable model (see Marginalisation.LATENT), but StudentTFamily does not implement the latent-conditional log_prob — its CONSUMES_LATENT_GP is False. ...

```

`LikelihoodFamily.CONSUMES_LATENT_GP` is `False` on the base class, so a
third-party family inherits the refusal rather than the defect. This is the
same discipline the unimplemented `RiceFamily` gets, applied to a combination
rather than to a family.

### The circular complex GP: declared analytic, **implemented at W4.2**

**Ruled by Peter, 2026-09-03** (§17 Q6, the interferometry sketch's
recommendation accepted): `complex_gaussian` + `GaussianProcessNoise`
declares `ANALYTIC`, with the **circular complex GP** — one real kernel
applied independently to the real and imaginary parts: equal component
covariances, zero pseudo-covariance — as the fixed meaning. That is the
declaration under which the flexible likelihood reaches the plan's Phase-4
proof modality, and fixing it at the freeze is what let the freeze be reviewed
against it.

Between the freeze and W4.2 the combination was **refused at composition with
the schedule named** — the declared-but-staged discipline `RiceFamily` gets for
a whole family, applied to one combination. **W4.2 implements it**, so
`GP_ANALYTIC_IMPLEMENTED` is `True` and the composition is the flexible
likelihood on visibilities:

```pycon
>>> ComplexGaussianFamily().marginalisation_with(gp_noise)
<Marginalisation.ANALYTIC: 'analytic'>
>>> Likelihood(ComplexGaussianFamily(), gp_noise).marginalisation
<Marginalisation.ANALYTIC: 'analytic'>

```

The closed form, and **why it is one factorisation rather than two**. Write
`S = K(θ) + diag(σ²)` for the per-component covariance over the observed
points. Circularity says the real `2N` covariance of `(Re r, Im r)` is
`diag(S, S)` — the same block twice, and zero off the diagonal — so

```
log p = -½ [ rᵉᵀ S⁻¹ rᵉ  +  rⁱᵀ S⁻¹ rⁱ  +  2 log|S|  +  2N log 2π ]
```

with `rᵉ`, `rⁱ` the components of `observed − predicted`. That is **one**
Cholesky of an `N × N` matrix, **two** triangular solves against it, and
`log|S|` evaluated once and counted twice. Forming the `2N × 2N` matrix and
factorising it would cost eight times the arithmetic in order to carry a zero
block the model has already declared, and calling the real solver twice would
compute `log|S|` twice. §7's `(n, k)` right-hand side is the mechanism: the
family hands the solver two columns and one covariance.

`σ` is the **per-component** standard deviation, which is what
`results_schema.md` §16 says a `VisibilitySet`'s real-valued uncertainty
encodes, and the kernel's `amplitude` is a per-component marginal standard
deviation for the same reason — `K` appears once per component above, so a
correlated calibration error of RMS `a` in each component is the kernel with
`amplitude = a`. The total modulus variance `E|r|²` is `2 (K_ii + σ²)`.

The **draw** (`sample`, §3) follows the same model and is therefore *two* real
GP realisations rather than one complex one: `x = μ + (L z₁ + σ w₁) + i (L z₁′ +
σ w₁′)`, the two components sharing `L` and sharing nothing else. Using one
realisation for both would give a draw perfectly correlated between the real and
imaginary parts, whose modulus statistics this density does not score.

`LikelihoodFamily.GP_ANALYTIC_IMPLEMENTED` remains the staging flag (default
`True`), so a future family that declares `ANALYTIC_WITH_GP` before writing its
closed form inherits the discipline rather than reinventing it; its refusal no
longer names Phase 4, because the instance the wording was written for has
landed.

### Enforcing it against the engine

§4.4: "gradient-free samplers cannot realistically handle hundreds of latent
values, so non-Gaussian + flexible-GP robustness is effectively a
modern-backend capability."

```pycon
>>> latent_poisson.check_engine(differentiable=False, engine="emcee")
Traceback (most recent call last):
    ...
ampere.core.exceptions.LikelihoodError: emcee cannot run this likelihood: the poisson family with a GaussianProcessNoise noise model introduces N latent values ...

```

It must not be a blanket refusal, though: a gradient-based engine is fine, and
so is *any* engine on an analytic likelihood.

```pycon
>>> latent_poisson.check_engine(differentiable=True, engine="NUTS")
>>> flexible.check_engine(differentiable=False, engine="emcee")

```

W1.7 owns the capability record (`differentiable`, `batchable`, `device`); this
method is what it calls with it. It takes plain keyword flags rather than
defining a capability type here, so the two contracts do not have to agree on a
class. Pass `observed=` as well when the data are to hand — see §9, where a
censored sample the mask excludes must not count against a gradient-free
engine.

## 5. Noise models

`IndependentNoise()` declares nothing and does nothing but read the container's
own uncertainties. Two optional knobs cover what a real fit needs — a `scale`
multiplying every error bar (the "the catalogue underestimates its
uncertainties" nuisance parameter) and a `jitter` added in quadrature (an
unmodelled noise floor). Both are ordinary parameters:

```pycon
>>> noise = IndependentNoise(scale=st.loguniform(0.5, 5.0), jitter=0.05)
>>> noise.parameters.free_names, noise.parameters["jitter"].is_fixed
(('scale',), True)
>>> inflated = Likelihood(GaussianFamily(), noise)
>>> round(inflated.log_prob(model, data, {"scale": 2.0}), 6)
2.434866

```

`GaussianProcessNoise` is the flexible likelihood. It takes a kernel and a
solver strategy, and it adopts the kernel's hyperparameters as its own
parameters:

```pycon
>>> gp_noise = GaussianProcessNoise(
...     Matern32(st.loguniform(1e-3, 1e1), st.loguniform(0.1, 10.0))
... )
>>> gp_noise.parameters.free_names
('amplitude', 'length_scale')
>>> gp_noise.CORRELATED, gp_noise.solver
(True, DenseGP(jitter=0.0))

```

Those `Parameter` objects are the kernel's own, shared **by name**, not by
Python identity. `prior_art.md` Tension 1 is explicit that gammapy's
identity-based sharing broke downstream consumers that had to deduplicate it,
and that ampere has more tree-walking consumers than gammapy did; a
`Parameter` is a frozen dataclass, so sharing the object is safe and the
*declaration* is what travels.

A `Likelihood` holds one flat parameter namespace — the family's and the noise
model's together, with no nesting:

```pycon
>>> composed = Likelihood(
...     StudentTFamily(nu=st.loguniform(2.0, 50.0)),
...     IndependentNoise(scale=st.loguniform(0.5, 2.0), jitter=st.halfnorm(0.0, 1.0)),
... )
>>> composed.parameters.free_names
('scale', 'jitter', 'nu')

```

Flat, because `parameters.md` §12.4 records that `ParameterSet.merge` is **not
associative**: merging a `ParameterMapping.merged` again silently drops the
first merge's bindings. W1.7 merges every dataset in one call, so this contract
must not consume a merge level. A genuine collision therefore raises rather
than being silently qualified — see §14.

### Prediction-aware noise: `predicted` reaches the noise model

**Ruled by Peter, 2026-09-03** (W1.11 gap X-1 — `awkward_instrument.md` §6's
detailed design, accepted as written and landed at W1.13): `NoiseModel.sigma`
and `NoiseModel.noise_params` take the **retained predicted values** as a
keyword-only argument:

```
sigma(observed, retain, values, *, predicted=None)
noise_params(observed, retain, values, *,
             predicted=None, coordinates=None, latent=None, limits=None)
```

`predicted` is the identical, already-excised array the family's `log_prob`
receives as its first argument — float64, or complex128 for a complex family
(a noise model wanting an amplitude takes `np.abs(predicted)` itself; ampere
does not project on its behalf). Every call site passes it:
`Likelihood.log_prob`, `Likelihood.conditional` (so W1.12's diagnostics see
the same effective σ the fit used) and W1.7's `Dataset.draw_observation`,
where the `NoiseParams` are built from the noiseless prediction *before*
noise is added — the draw is σ(μ), not σ(x), the standard generative
reading. `IndependentNoise` and `GaussianProcessNoise` ignore it, so the
simple path is unchanged.

The case this exists for is a noise whose magnitude depends on the model.
The standard library names **`FractionalModelNoise`** —
`sigma_eff² = (s·σ_data)² + (f·predicted)²`, the single most requested thing
missing from legacy ampere's likelihood — the way `transformations.md` §10
names standard chain steps: the contract fixes the name and semantics here,
and the ten-line implementation lands with the reference backend in Phase 2
(with `f` an ordinary fitted parameter). **Landed W2.1** as
`ampere.backends.reference.FractionalModelNoise`, with `f` declared through
`register_parameter` — so a prior fits it and a number holds it fixed — and the
GP composition below as `FractionalModelGPNoise`. The doctest that follows
keeps `f` fixed as a plain attribute because it is illustrating the mechanism
in ten lines; the shipped class does not. A prediction-dependent σ is still
diagonal, so `GaussianFamily` plus a fractional noise stays `ANALYTIC` — the
marginalisation machinery never inspects *how* σ was computed — and the GP
composition ("10 % model error *and* a misspecification GP") is a
`GaussianProcessNoise` subclass overriding only `sigma`: the marginal
likelihood is `N(0, K + diag(σ_data² + (f·μ)²))` with no new mathematics.

```pycon
>>> class FractionalModelNoise(IndependentNoise):
...     """sigma_eff**2 = sigma_data**2 + (f * predicted)**2, with f fixed."""
...     def __init__(self, f):
...         super().__init__()
...         self.f = f
...     def sigma(self, observed, retain, values, *, predicted=None):
...         base = super().sigma(observed, retain, values)
...         return np.sqrt(base**2 + (self.f * predicted) ** 2)
>>> fractional = Likelihood(GaussianFamily(), FractionalModelNoise(0.1))
>>> fractional.marginalisation
<Marginalisation.ANALYTIC: 'analytic'>
>>> round(fractional.log_prob(model, data), 6)
1.842006

```

The boundary holds in both directions: the noise model sees the prediction,
but never the model's *parameters* beyond those the likelihood declares
(§15.10). A noise term that depends on a physical parameter is tied to a
likelihood parameter or expressed as a transformation.

## 6. Kernels: declared neutrally, hyperparameters are ordinary parameters

A kernel is a *declaration*, not an implementation detail: a neutral family
name plus its hyperparameters, which is exactly what W1.9's lowering table
needs to emit a celerite2 / tinygp / GPyTorch term and the most that translates
across all three.

```pycon
>>> kernel = Matern32(st.loguniform(1e-3, 1e1), st.loguniform(0.1, 10.0))
>>> kernel.spec()
KernelSpec(family='matern32', hyperparameters=('amplitude', 'length_scale'),
           quasiseparable=True, axes=None, terms=())
>>> kernel.spec().to_dict()
{'family': 'matern32', 'hyperparameters': ['amplitude', 'length_scale'], 'quasiseparable': True}

```

`parameters.md` §13 instructed this contract that "GP hyperparameters
(amplitude, length-scale) are ordinary parameters on a `Parameterised` noise
model, with `Log` bijections". Confirmed, and enforced in one place so no
kernel can quietly do it differently:

```pycon
>>> kernel.parameters.free_names
('amplitude', 'length_scale')
>>> kernel.parameters.bijections()
(Log(lower=0.0), Log(lower=0.0))

```

Nothing about them is special, which is the point: tying, fixing, plates,
serialisation and W1.9's lowering all work on them unchanged. A bare number
fixes one; a ready-made `Parameter` is accepted as-is (with its name checked,
because hyperparameter names are part of the `KernelSpec`):

```pycon
>>> Matern32(0.5, st.loguniform(0.1, 10.0)).parameters.free_names
('length_scale',)
>>> Matern32(st.halfnorm(0.0, 1.0), "a metre or so")
Traceback (most recent call last):
    ...
ampere.core.exceptions.LikelihoodError: kernel hyperparameter 'length_scale' must be a frozen scipy.stats distribution (a prior), a number (held fixed), or an ampere Parameter — got str.

```

### Matérn-3/2, and why

`k(r) = a² (1 + √3 r/ℓ) exp(−√3 r/ℓ)`, with `a` the marginal **standard
deviation**, so `k(0) == a²` — celerite2's `Matern32Term(sigma=…, rho=…)`
convention, chosen so a prior on the amplitude is a prior in the data's own
units.

```pycon
>>> float(kernel.matrix([[0.0]], [[0.0]], {"amplitude": 2.0, "length_scale": 1.0})[0, 0])
4.0

```

`DEVELOPMENT_PLAN.md` §2 makes it the default throughout, replacing legacy's
hardcoded RBF, for two independent reasons:

1. **It represents structured residuals better.** A Matérn-3/2 sample path is
   once differentiable, not analytic; real model deficiencies are not
   infinitely smooth. `prior_art.md` lesson S1 is independent, published,
   empirically-validated corroboration: Starfish (Czekala, Andrews, Mandel,
   Hogg & Green 2015, ApJ 812, 128) arrived at exactly `K_G(r) = w · a_G ·
   (1 + √3 r/ℓ) exp(−√3 r/ℓ)` in velocity separation, for exactly this purpose.
   The closest prior art to ampere's core scientific idea reached the same
   kernel for the same reason.
2. **It is exactly quasiseparable.** Matérn-3/2 has an exact representation as
   a sum of celerite/SHO terms, which is what makes `QuasisepGP` an *exact*
   O(N) solve rather than an approximation (§4.4, and §7 below). The squared
   exponential has no such form, and that asymmetry is the concrete reason the
   default changed:

```pycon
>>> Matern32(1.0, 1.0).spec().quasiseparable
True
>>> SquaredExponential(1.0, 1.0).spec().quasiseparable
False

```

The squared exponential is kept so milestone M2's misspecification study can
compare the new default against what legacy actually did, and so a user who
wants it can have it — on `DenseGP` only.

### The kernel family, and the registry that makes one fast (*Amended W4.5*)

*(W4.5, 2026-09-11. §15.1 and §15.3 recorded "one kernel per noise model" and
"only Matérn-3/2 and the squared exponential ship" as deliberate limitations
with named extension points. This is those extension points, built.)*

Seven families ship, and the table's last column is the one that matters: a
kernel reaches the O(N) path if and only if ampere holds an **exact**
semiseparable representation for it.

| Kernel | `k(τ)` | Rank | Sample path |
|---|---|---|---|
| `Matern12` | `a² e^{−τ/ℓ}` | 1 | continuous, nowhere differentiable |
| `Matern32` | `a² (1 + √3 τ/ℓ) e^{−√3 τ/ℓ}` | 2 | once differentiable |
| `Matern52` | `a² (1 + √5 τ/ℓ + 5τ²/3ℓ²) e^{−√5 τ/ℓ}` | 3 | twice differentiable |
| `SHO` | `a² e^{−τω₀/2Q}[cos ηω₀τ + (4Q²−1)^{−½} sin ηω₀τ]`, `η = √(1−1/4Q²)` | 2 | damped periodic, `Q > ½` |
| `RotationTerm` | two `SHO`s, at a damped period and its first harmonic | 4 | non-sinusoidal periodic |
| `SquaredExponential` | `a² e^{−τ²/2ℓ²}` | — | analytic; **no** exact form |
| `SpectralMixture` | `Σᵢ SHOᵢ` with free frequencies | Σ ranks | as its components |

`amplitude` is the marginal standard deviation in every one of them, so
`k(0) == a²` throughout and a prior on it is a prior in the data's own units.
`SHO`'s `period` is the *undamped natural* period `2π/ω₀` (celerite2's `rho`);
`RotationTerm`'s is the *damped* period of its fundamental, which is
celerite2's own parameterisation and exactly `period` by construction.

```pycon
>>> float(SHO(0.5, 2.0, 4.0).diagonal(np.zeros(1), {})[0])
0.25
>>> quasiseparable_families()
('matern12', 'matern32', 'matern52', 'rotation', 'sho', 'spectral_mixture', 'sum')

```

**Matérn-5/2 is exactly quasiseparable, and §15.3 said otherwise.** The
correction is the same one W2.3 made for Matérn-3/2 and is worth stating
plainly. There is no exact Matérn-5/2 in the **celerite basis**
`e^{−cτ}(a cos dτ + b sin dτ)`, which has no `τ² e^{−cτ}` member — that much
of §15.3 was right. But the solver underneath factorises any rank-J
**semiseparable** matrix `K[n,m] = Σⱼ U[n,j] V[m,j] e^{−cⱼ(tₙ−tₘ)}`, and a
degree-2 polynomial in `tₙ − tₘ` is a rank-3 bilinear form in `(1, t, t²)`:
with `f = √5/ℓ`,

```
1 + f(tₙ−tₘ) + ⅓f²(tₙ−tₘ)²
  = [1 + f tₙ + ⅓f² tₙ²]·1 + [−f − ⅔f² tₙ]·tₘ + [⅓f²]·tₘ²
```

so `Uₙ = a²(1 + f tₙ + ⅓f² tₙ², −f − ⅔f² tₙ, ⅓f²)`, `Vₘ = (1, tₘ, tₘ²)`,
`c = (f, f, f)`. Algebra, not a limit. The quadratic generators cancel one
order worse than Matérn-3/2's, which the midpoint re-referencing of §7 bounds;
it is the term to watch when data span very many length scales.

**`register_quasiseparable_term`: the table, opened.** Before W4.5 a user could
define a `Kernel` subclass and use it — the ABC is public and `matrix` /
`diagonal` are all `DenseGP` needs — but could not make it fast, because the
celerite translation was a private dict. It is now a registry of the
lowering/realisation registries' shape: one slot per kernel family, no silent
overwrite, `override=` to replace deliberately, `builtin=` rows distinguished,
and a provenance entry per user row (`term_provenance_entries()`, passed to
`provenance_attrs`'s `extra=` as `registered_quasiseparable_terms`). A builder
takes `(kernel, values, axis)` and returns a `CeleriteRepresentation` — `c`,
`U`, `V` and `k(0)`.

**One slot per family, not per (family, backend), and this is a decision.** A
semiseparable representation is generators built from the coordinate and the
hyperparameters with `exp`, `cos`, `sin` and `stack`: the same mathematics in
numpy, torch and jax. What differs between backends is the *array namespace*,
so that becomes an argument — `ArrayOps`, carried by the kernel instance
(`Kernel.ops`) and rebound by each solver to its own through
`Kernel.with_ops` — rather than a second registry key. Three consequences, all
of them wanted: a user registers **once** and reaches the O(N) path on all
three backends; the three backends cannot drift, where before W4.5 each carried
its own transcription with a comment saying the tables "must carry exactly the
same families" and nothing able to check it; and a backend that genuinely needs
to specialise can still register its own builder under `override=True`.
`StationaryKernel` is public for the same reason: a two-hyperparameter kernel
of the usual shape should not have to redeclare the units and bijections.

### Kernel algebra: `Sum` and `Product` (*Amended W4.5*)

```pycon
>>> composed = Sum(Matern32(st.loguniform(1e-3, 1e1), st.loguniform(0.1, 10.0)),
...                SHO(st.halfnorm(0.0, 0.1), 0.02, 30.0))
>>> composed.parameters.free_names
('term0.amplitude', 'term0.length_scale', 'term1.amplitude')
>>> composed.QUASISEPARABLE
True

```

A composite qualifies its children's hyperparameters with a **term label** —
`term0.amplitude` — which is the convention `ParameterSet.merge` already uses
for component parameters, so two Matérn terms in one sum cannot collide.
Positional indices (`terms.0.amplitude`) are not available: `parameters.md`'s
name rule requires every dot-separated segment to be a Python identifier,
because parameter values are handed to models as keyword arguments, and `0` is
not one. `labels=("broad", "narrow")` replaces the defaults where the names
carry meaning, and they do reach the sampler's coordinates and the stored
chains.

`KernelSpec` of a composite is a **tree**, `(label, spec)` pairs in declaration
order. Order-preserving rather than canonicalised, **even for `Sum`, where
addition commutes** — because the labels are what the children's parameters are
named with, so reordering a sum renames its parameters and is a different
declaration to everything downstream. Sorting the tree for the hash would let
two runs whose chains cannot be compared claim one identity, which is the thing
the hash exists to prevent.

**A sum of quasiseparable terms is quasiseparable**; the generators concatenate
and the ranks add, so an arbitrarily rich sum still costs O(N) with the
constant growing in the total rank. `Sum.QUASISEPARABLE` is therefore an
instance-level conjunction of its terms', not a class-level claim.

**A product is not**, and is refused on the O(N) path by name:

```pycon
>>> chromatic = Product(Matern32(1.0, 2.0, axes=("u", "v")),
...                     Matern32(0.3, 0.05, axes=("spectral_axis",)))
>>> chromatic.QUASISEPARABLE
False

```

## 7. Solver strategies

The solver is swappable behind the noise model precisely so the scaling story
can change without the science code changing. Every strategy computes the same
quantity — `log N(residual; 0, K(θ) + diag(σ²))` — and `DenseGP` is the
definition of the right answer.

| Strategy | Exact? | Applies to | Status |
|---|---|---|---|
| `DenseGP` | yes | anything, O(N³) | **implemented** — the correctness anchor |
| `QuasisepGP` | yes | one ordered axis, quasiseparable kernels and `Sum`s of them, O(N) | **implemented** (W2.3, celerite2; six families and `Sum` since W4.5) — `conditional_loo` deferred on the numpy solver; supplied by both differentiable backends' own |
| `HilbertSpaceGP` | **no** | stationary kernels with a closed-form spectral density (the three Matérns, `SquaredExponential`, `SHO`, `Sum`s of those), 1–3 axes, O(N m + m³) | **implemented** (W5.4) on all three backends — `conditional_loo` exact in the approximation; latent block of size `m` |
| `WindowedSparseGP` | no | any kernel, any dimension | slot; see below |
| `InducingPointGP` (SVGP) | no | 2D+ | slot; Phase 5 |
| `StructuredGridGP` (SKI) | no | gridded 2D+ | slot; Phase 5 |
| `VecchiaGP` | no | 2D+ | slot; Phase 5 |

```pycon
>>> (DenseGP.EXACT, QuasisepGP.EXACT, HilbertSpaceGP.EXACT, WindowedSparseGP.EXACT)
(True, True, False, False)
>>> (HilbertSpaceGP.IMPLEMENTED, WindowedSparseGP.IMPLEMENTED)
(True, False)

```

The two exact strategies compute the *same* number by genuinely different
recursions — a dense Cholesky against a rank-2 quasiseparable factorisation —
which is the equivalence `DenseGP` exists to define, and the conformance
battery's `DenseGP`↔`QuasisepGP` row:

```pycon
>>> dense = Likelihood(
...     GaussianFamily(), GaussianProcessNoise(Matern32(0.3, 2.0), DenseGP())
... )
>>> quasisep = Likelihood(
...     GaussianFamily(), GaussianProcessNoise(Matern32(0.3, 2.0), QuasisepGP())
... )
>>> abs(quasisep.log_prob(model, data) - dense.log_prob(model, data)) < 1e-9
True

```

A slot is a fixed interface with no implementation, and it says so rather than
pretending:

```pycon
>>> Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(0.3, 2.0), WindowedSparseGP()))\
...     .log_prob(model, data)
Traceback (most recent call last):
    ...
ampere.core.exceptions.LikelihoodError: WindowedSparseGP is a declared strategy slot with no implementation yet ...

```

Applicability is checked at **composition** time, and a permanent
incompatibility is reported ahead of a temporary one — a user who paired
`QuasisepGP` with a squared exponential needs to hear about the kernel, not
about the schedule:

```pycon
>>> GaussianProcessNoise(SquaredExponential(0.3, 2.0), QuasisepGP()).check_compatible(
...     GaussianFamily(), data
... )
Traceback (most recent call last):
    ...
ampere.core.exceptions.LikelihoodError: QuasisepGP needs a kernel with an exact quasiseparable representation, but SquaredExponential (squared_exponential) has none. ...

```

### The first approximate strategy (*Amended W5.4*)

`HilbertSpaceGP` is the first strategy in this table with `EXACT = False` and
`IMPLEMENTED = True`, and it makes concrete three things the contract had only
declared. It is the reduced-rank spectral method of Solin & Särkkä (2020) in
the practical form Riutort-Mayol et al. (2023) give: on a box containing the
data, the Dirichlet Laplacian's eigenfunctions diagonalise any stationary
kernel through its spectral density, so `K ≈ Φ diag(S) Φᵀ` with `m` basis
members and every quantity below follows from Woodbury at `O(N m + m³)`.

**What "approximate" means here, and what it does not.** It is the *exact*
marginal likelihood of a rank-`m` Gaussian process that converges to the
declared kernel as the basis fills the box. It is therefore not a different
quantity, and not a heuristic: `conditional_loo` is computed from the same
Woodbury factorisation and is exact in the approximation, rather than refused
the way `QuasisepGP`'s is; `condition` returns the approximation's own prior
variance, so the posterior variance it reports is non-negative at every `m`,
which mixing an exact prior with an approximate correction would not be.

**The approximation is a declaration, not configuration.** `basis_size` (the
members per axis; the total `m` is their product) and `boundary_factor` (the
box's half-width as a multiple of the data's own half-extent) are dataclass
fields, so `Likelihood.to_spec` records them and two runs at different `m` have
different spec hashes — they change the answer, and fold-in 10's test is
whether two backends may legitimately differ on a value, which here they may
not. They are repeated in `provenance_config()` because a stored run's attrs
are where a reader looks for "what approximation produced this".

```pycon
>>> solver = HilbertSpaceGP(basis_size=64, boundary_factor=3.0)
>>> approximate = Likelihood(
...     GaussianFamily(), GaussianProcessNoise(Matern32(0.3, 2.0), solver)
... )
>>> exact = Likelihood(
...     GaussianFamily(), GaussianProcessNoise(Matern32(0.3, 2.0), DenseGP())
... )
>>> abs(approximate.log_prob(model, data) - exact.log_prob(model, data)) < 1e-2
True
>>> HilbertSpaceGP(basis_size=(16, 16)).provenance_config()
{'basis_size': [16, 16], 'boundary_factor': 2.0}

```

The `boundary_factor=3.0` above is not decoration, and it is the trap worth
seeing once: `data` spans 1–5 µm, so its half-extent is 2 µm — the *same* as
the kernel's length scale — and at the default `boundary_factor` the error
stalls at about 0.14 nat however large `basis_size` grows, because what limits
it is the box rather than the basis. Widening the box to three half-extents
brings it to 8e-3 at `basis_size=32` and it barely improves after that. The
two parameters are one approximation, and refining only one of them converges
to the wrong answer rather than slowly.

The conformance battery holds it to a **convergence** rather than to a number,
which is the only honest assertion about an approximation: the disagreement
with `DenseGP` must fall as the basis is refined, inside an envelope anchored
on the coarsest setting of the same sweep and floored at the box's own
truncation error, and the finest setting must actually be close. See
`tests/conformance/README.md` §3's approximation class.

A kernel with no closed-form spectral density is refused **by name, at
composition**, and the refusal names the node rather than the leaf: a
`Product`'s factors each have a spectral density and the product does not,
because its transform is a convolution.

```pycon
>>> product = Product(Matern32(0.3, 2.0), Matern12(0.2, 0.5))
>>> GaussianProcessNoise(product, HilbertSpaceGP()).check_compatible(GaussianFamily(), data)
Traceback (most recent call last):
    ...
ampere.core.exceptions.LikelihoodError: HilbertSpaceGP needs a stationary kernel with a closed-form spectral density, but Product (product) has none. ...

```

The latent block this strategy declares is the **basis size**, not the sample
count; `inference.md` §17.4 carries that ruling and its consequences.

### The right-hand side may carry several realisations (*Amended W4.2*)

Every solver method that takes a `residual` accepts either an `(n,)` array —
one realisation, which is every call written before W4.2 — or an `(n, k)` block
of `k` realisations that are **independent of one another and share this one
covariance**. The declared quantities follow from that and from nothing else:

| method | `(n,)` | `(n, k)` |
|---|---|---|
| `log_marginal_likelihood` | the marginal | the **sum** of the `k` marginals: `log|K+Σ|` enters `k` times, the factorisation happens once |
| `conditional_loo` | one term per sample | still one term per **sample**, summing the `k` components, which share `A_ii` |
| `condition` | `(m,)` mean, `(m,)` variance | `(m, k)` mean, `(m,)` variance — a posterior variance does not depend on the data |
| `latent_transform` | `(n,)` → `(n,)` | `(n, k)` → `(n, k)` |

This exists for §4's circular complex GP, where `k = 2` and the two columns are
the real and imaginary parts of one complex residual. It is not a convenience
API: a two-column right-hand side is the only way to state "these two
realisations share this covariance" to a strategy whose whole job is to
factorise the covariance once.

`GPSolver.STACKED_RESIDUALS` is the declaration that a strategy can take the
extra columns. It is `False` by default, which is the honest answer for a
strategy written before the rule existed and for every declared slot, because
the two ways of getting it wrong are both silent: flattening `(n, 2)` scores
`2n` residuals against an `n × n` covariance, and taking the first column drops
the imaginary part of every visibility. `DenseGP` declares `True` on all three
backends.

```pycon
>>> (DenseGP.STACKED_RESIDUALS, QuasisepGP.STACKED_RESIDUALS)
(True, False)

```

### Why the O(N) path is refused on visibilities, word for word (*W4.2*)

```
the complex_gaussian family on a complex VisibilitySet is the circular complex
GP: the real and imaginary parts are two independent real processes sharing one
covariance, so the solve takes a two-column right-hand side, and QuasisepGP
declares STACKED_RESIDUALS = False. QuasisepGP could not take them even in
principle: it needs one ordered one-dimensional coordinate, and a visibility
lives at a point of the (u, v) plane at a wavelength, which no ordering reduces
to one coordinate. Use DenseGP, with the kernel selecting the axes it acts on --
Matern32(axes=("u", "v")) for an isotropic (u, v) kernel, or
Product(Matern32(axes=("u", "v")), Matern32(axes=("spectral_axis",))) for an
error that is smooth in (u, v) and sharp in wavelength.
```

Raised by `GaussianProcessNoise.check_compatible` **before** the solver's own
`check_compatible`, for the reason the `Product` refusal precedes the generic
quasiseparable one: the generic `REQUIRES_ORDERED_1D` message ends "select
exactly one axis to reach the O(N) path", and on this modality that is advice
nobody can take. `REQUIRES_ORDERED_1D` cannot hold for a visibility set,
whatever the kernel selects, so the refusal says so and names the solver that
does work.

A solver **may** branch on `Axis.regular` or `Axis.log_regular` for a fast path
and **must** have a path when both are false — `results_schema.md` §16's
instruction, and `architecture.md` §7's general rule that a fast path may be
taken but never required. `DenseGP` takes no such branch at all, which is the
trivial way to comply and one more reason it is the oracle. `QuasisepGP` takes
none either: it sorts the coordinates it is handed (a Gaussian density is
invariant under a simultaneous permutation of residuals, variances and
coordinates) and runs the same recursion whether or not the axis is regular.

### `QuasisepGP`: the exact representation, and what it does *not* do (W2.3)

The O(N) path is celerite2's numpy solver, and the term ampere hands it is
**its own**, not `celerite2.terms.Matern32Term`. celerite2's is an
approximation and says so: the celerite basis `e^{−cτ}(a cos dτ + b sin dτ)`
has no `τ e^{−cτ}` member, so its Matérn-3/2 is a limit in a parameter `eps`,
which at the default `eps=0.01` costs about `5e-3` in the log-likelihood — 
three orders of magnitude outside `tolerances.cross_solver`.

The *solver* underneath needs no such approximation: it factorises any rank-J
matrix of the form `K[n,m] = Σ_j U[n,j] V[m,j] exp(−c_j (t_n − t_m))` for
`n > m`, and Matérn-3/2 has an **exact** rank-2 representation in it. With
`f = √3/ℓ` and `Δ = t_n − t_m > 0`,

```
k(Δ) = a²(1 + fΔ)e^{−fΔ} = e^{−f(t_n − t_m)} [ a²(1 + f t_n)·1 + (−a² f)·t_m ]
```

so `U_n = (a²(1 + f t_n), −a² f)`, `V_m = (1, t_m)`, `c = (f, f)` and the
diagonal is `a² + σ_n²`. The bracket is `a²(1 + f t_n − f t_m) = a²(1 + fΔ)`:
an algebraic identity, not a limit. That is the sense in which §6's "exactly
quasiseparable" is cashed in, and the conformance row is what holds it.

Two honest caveats, both properties of the representation rather than of the
implementation:

* The generators grow linearly in the coordinate, so `U_n · V_m` is a
  difference of large numbers when `f t ≫ 1`. ampere re-references the
  coordinates to the midpoint of their own range (exact, since only
  differences enter), which bounds the cancellation by half the number of
  length scales the data span: measured against the dense Cholesky, ~5e-12
  over 10 length scales, ~3e-10 over 10², ~2e-8 over 10⁴.
* **`conditional_loo` is deferred on this (numpy) solver**
  (`DEVELOPMENT_PLAN.md` §2, 2026-09-05). Every leave-one-out term needs
  `A_ii` for `A = (K + diag(σ²))⁻¹`, and celerite2's public numpy interface
  exposes no O(N) route to that diagonal — its own conditional variance forms
  the cross-covariance densely. `ampere.core.QuasisepGP` therefore refuses,
  naming `DenseGP`, rather than quietly costing O(N²) under an O(N) name.
  `log_marginal_likelihood` and `latent_transform` are O(N); `condition` is
  O(N·M) for M evaluation points, because the cross-covariance block is dense
  by construction whatever the solver.

  **A backend that calls celerite2's compiled kernels directly is not bound by
  that deferral, and neither differentiable backend is** — the torch solver
  since W2.4 slice 2 (`celerite2.backprop`) and the jax one since W2.5 slice 3
  (`celerite2.jax.ops`), both by the same O(N) backward accumulation over the
  factorisation, both checked against `DenseGP`'s Cholesky. The deferral was
  always about a coupling to celerite2's internal factorisation convention
  rather than about the mathematics, and a backend that already holds `c, U,
  d, W` has paid it. A solver is free either way: the conformance row states
  the whole contract as "refuses by name **or** agrees with `DenseGP` at
  `tolerances.cross_solver`".

### Windowed-sparse truncation: a named slot, and why it is not the default

`prior_art.md` lesson S2 asked this contract either to give Starfish's actual
strategy a named slot or to write the "why we didn't use this" paragraph. It
gets both, because a Starfish-aware reviewer deserves the argument and not just
the class.

Starfish avoids the O(N³) solve by applying a Hann window that tapers the
kernel to **exactly zero** beyond a cutoff radius (`r₀ = 4ℓ` for the global
kernel), producing a genuinely sparse banded covariance a sparse Cholesky
factorises in roughly linear time. That is a real O(N) strategy and a
*qualitatively different* one from `QuasisepGP`: approximate rather than exact,
needing no state-space machinery, easy to reason about locally, and applicable
to kernels with no quasiseparable form and to more than one dimension. It is
therefore `WindowedSparseGP`, a named slot, not a footnote.

It is not the default for three reasons:

1. **`QuasisepGP` is exact for the kernel ampere recommends.** Matérn-3/2 has
   an exact celerite representation, so the windowed approximation buys nothing
   at the same asymptotic cost — you would be paying an approximation error for
   a speed you can have exactly.
2. **The cutoff is one more hyperparameter nobody chooses for the user.**
   `r₀ = 4ℓ` is a rule of thumb that interacts with a *fitted* `ℓ`: the
   sparsity pattern, and hence the answer, changes as the sampler moves. A
   strategy whose approximation error varies over the posterior is a bad
   default even when it is a good tool.
3. **The truncated tails are exactly what the kernel is for.** The Matérn tail
   is the statement "residuals stay correlated a long way out"; discarding it
   biases the length-scale posterior, which is itself a diagnostic output
   (§4.8, family C).

Where it earns its slot is the case `QuasisepGP` cannot reach: a
non-quasiseparable kernel, or 2D+ data where the state-space recursion has no
analogue. Phase 5 should implement it alongside the inducing-point strategies
rather than instead of them.

### The `axes` selector: which coordinates a kernel acts on (*Amended W4.5*)

*(Ruled by Peter 2026-09-11; `phase4_placement_memo.md` §3.6 item 3. This lifts
§15.2's second half — "mixed-unit axes are refused outright" — without lifting
its first: a kernel is still isotropic in the coordinates it sees, and there is
still one `length_scale`.)*

Every `Kernel` takes an `axes=` selector naming the container axes it acts on.
The default is `None`, meaning **every axis** — the pre-W4.5 reading, so no
declaration written before W4.5 changes, and `KernelSpec.to_dict` emits the
`axes` key only when a selection was made, so **no spec hash moves either**.

The selection changes two of `check_compatible`'s rules, and both changes are
the same change: the rule now applies per **leaf kernel**, over the axes that
leaf selects, rather than over the container's.

* **The single-unit rule.** A Euclidean separation across axes in different
  units is meaningless and is still refused. What the selector adds is that a
  container carrying `u`, `v` (dimensionless) and `spectral_axis` (micron) is
  refused for a bare isotropic kernel *exactly as before*, and accepted for
  `Matern32(axes=("u", "v"))`. The refusal's own message pointed at this fix
  before W4.5 existed: "use one axis, or declare a kernel that takes a
  length-scale per axis".
* **`REQUIRES_ORDERED_1D`.** The count is of the axes the kernel *tree* uses,
  not the container's, so a quasiseparable solve over the spectral axis of a
  three-axis container is expressible.

The hyperparameter-unit check follows the selection too: a `length_scale` is
compared against the unit of the axes **its own leaf** selects, so a product of
a dimensionless (u, v) kernel and a micron spectral one has two length scales
in two units and both are right.

Names become column indices once, at composition, in
`GaussianProcessNoise.kernel_for` — which is where the container still is; the
solver is handed a bare `(n, d)` block and no longer knows what the columns are
called. The binding is **functional**: `Kernel.for_axes` returns a bound copy,
so one kernel declaration may be shared between datasets whose axis orders
differ.

The case this exists for is `phase4_placement_memo.md` §3.6's chromatic sky
error. A missing patch of sky with a spectral profile contributes
`δV(B, λ) = S(λ)·F(B/λ)`: smooth in spatial frequency on the scale of the
patch, sharp in wavelength on the scale of the band. Its covariance is a
product of two kernels on **disjoint axis subsets** of one container, and no
isotropic kernel over all three axes can express it.

```python
sky_error = Product(
    Matern32(axes=("u", "v"), amplitude=..., length_scale=...),          # the patch
    Matern32(axes=("spectral_axis",), length_scale=...),                 # the band
)
```

A product's marginal variance is the product of its terms', so two terms each
declaring an `amplitude` over-parameterise it by one degree of freedom; fix all
but one (pass a number rather than a prior), as for any other multiplicative
redundancy.

### Why a product is refused on the O(N) path, word for word (*Amended W4.5*)

```
QuasisepGP cannot lower a Product: a product of quasiseparable kernels is not
quasiseparable. Where the factors act on different axes — which is what a
Product is for — the result is not a function of one ordered coordinate at all,
and where they act on the same one the semiseparable rank multiplies and is not
recoverable from the factors' own representations. Use DenseGP, or replace the
Product with a Sum, which is quasiseparable exactly when every term is.
```

The refusal is raised **before** the generic `QUASISEPARABLE` one, because
"products are not quasiseparable" is a sharper diagnosis than "this kernel is
not" — and because a user who reaches for a product is usually reaching for the
two-axis case, where no amount of registering would help.

## 8. Masks: `weights()`, and excision

`results_schema.md` §16 asked this contract to pick one of `weights()` and
`masked_uncertainty()` and state it. **It is `weights()`**, and masked samples
are removed by **row/column excision**.

The effective weight is the product of the predicted and observed containers'
weights, so a sample is used only if *both* call it valid; since those weights
are exactly 0 or 1, the product is an inclusion indicator, and retained samples
are the only rows and columns that enter the covariance at all.

```pycon
>>> masked = Spectrum(
...     [1.0, 2.0, 3.0, 4.0] * u.um,
...     [3.0, 2.5, 2.2, 1.4] * u.Jy,
...     uncertainty=[0.1, 0.1, 0.1, 0.1] * u.Jy,
...     mask=np.array([False, True, False, False]),
... )
>>> masked.weights().tolist()
[1.0, 0.0, 1.0, 1.0]
>>> excised = Spectrum(
...     [1.0, 3.0, 4.0] * u.um,
...     [3.0, 2.2, 1.4] * u.Jy,
...     uncertainty=[0.1, 0.1, 0.1] * u.Jy,
... )
>>> gp = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(0.4, 2.0)))
>>> a = gp.log_prob(masked.with_values([3.05, 2.40, 2.25, 1.35]), masked)
>>> b = gp.log_prob(excised.with_values([3.05, 2.25, 1.35]), excised)
>>> bool(np.isclose(a, b, rtol=0.0, atol=1e-12))
True

```

**Why excision and not the infinite-variance limit**, which
`results_schema.md` §7 offers as the equivalent statement.

The two are equivalent statements about a **chi-square term** — which is
exactly what §7 claims, and how it phrases `masked_uncertainty()`: "for
consumers whose algebra divides by σ² rather than multiplying by a weight".
`(r_i/σ_i)² → 0` as `σ_i → ∞`, and a zero weight deletes the same quantity.

They are **not** equivalent statements about a *normalised log-density*, and a
likelihood is a normalised log-density. Each Gaussian term carries a
`−log σ_i` alongside its chi-square, so

```
−½log(2π) − log σ_i − ½(r_i/σ_i)²  →  −∞   as σ_i → ∞,
```

which is not zero. Inflating an uncertainty does not remove a sample from a
log-likelihood; it makes that sample's contribution diverge. The GP case is the
same statement with a log-determinant instead of a `log σ`: letting `σ_i → ∞`
in `log N(r; 0, K + diag(σ²))` leaves the remaining samples' conditional
structure correct but adds `−½ log(2π σ_i²)`, so the limit approaches
*(excised value) − ½log(2π σ_i²)* and diverges rather than converging to the
excised value.

So the choice is not "the two rules agree on the diagonal and disagree under a
GP" — they disagree in both cases, by exactly the same divergent constant.
`results_schema.md` §7 defines a masked sample as carrying **exactly** zero
information; only excision delivers that, for any noise model. `weights()` is
the convention because a 0/1 weight multiplying a whole *term* (chi-square and
normalisation together) is excision written arithmetically, whereas
`masked_uncertainty()` is a statement about the chi-square alone. Both remain
correct for what §7 says they are; only one of them is a likelihood.

Consequences worth stating:

- A masked sample has no influence *whatever its value or its uncertainty*.
  Changing either changes nothing.

```pycon
>>> corrupted = Spectrum(
...     [1.0, 2.0, 3.0, 4.0] * u.um,
...     [3.0, 1e6, 2.2, 1.4] * u.Jy,
...     uncertainty=[0.1, 1e-9, 0.1, 0.1] * u.Jy,
...     mask=np.array([False, True, False, False]),
... )
>>> c = gp.log_prob(corrupted.with_values([3.05, 0.0, 2.25, 1.35]), corrupted)
>>> bool(np.isclose(a, c, rtol=0.0, atol=0.0))
True

```

- A fully masked pair returns `0.0`. That is the correct limit of the same
  rule, not a special case: no data, no information, no contribution.
  (W1.7's `FittingProblem` refuses to let a *parameter-dependent* mask reach
  this limit — `inference.md` §8's evaluation-invariance check — because a
  free `0.0` beats every finite log-likelihood.)

- **Every array in `NoiseParams` covers the retained samples only**, and a
  family that carries per-sample data of its own — a background spectrum, an
  instrumental template, a per-sample weight from outside ampere — must
  excise it with `noise.retain` (ruled 2026-09-02, W1.11 gap X-2): the
  boolean inclusion indicator over the *full* containers, set by this
  contract from the same weights product above, because `Likelihood` cannot
  know about data it was never handed. Failing to is a silent misalignment
  as soon as anything is masked.

```pycon
>>> nothing = Spectrum(
...     [1.0, 2.0] * u.um, [1.0, 1.0] * u.Jy,
...     uncertainty=[0.1, 0.1] * u.Jy, mask=np.array([True, True]),
... )
>>> gp.log_prob(nothing.with_values([0.0, 0.0]), nothing)
0.0

```

- **Masked coordinates remain available as evaluation locations.** The
  conditioned GP mean (§4.8's family C) is computed on the *full* axis,
  including the excluded samples, because "what would the GP have said here?"
  is exactly what a user asks about a region they masked.

- Ampere never threads sentinel NaNs through the arithmetic
  (`results_schema.md` §7, `prior_art.md` R1). A non-finite value in a
  *retained* sample is a loud error, not a quietly poisoned posterior:

```pycon
>>> broken = Spectrum(
...     [1.0, 2.0] * u.um, [1.0, np.nan] * u.Jy, uncertainty=[0.1, 0.1] * u.Jy
... )
>>> Likelihood(GaussianFamily()).log_prob(broken.with_values([1.0, 1.0]), broken)
Traceback (most recent call last):
    ...
ampere.core.exceptions.LikelihoodError: observed values contains non-finite entries. ...

```

### Per-observation terms: `pointwise_log_prob`

**Ruled by Peter, 2026-09-03** (`results.md` §15 R2 — a §4.4 addition,
recorded in the plan's decision log) and landed at the freeze:
`Likelihood.pointwise_log_prob(predicted, observed, values)` returns one
log-likelihood term per **retained** sample, on request — never stored by
default, and never called from the hot loop. Two decompositions, each under
the name `results.md` §6 reserves for it:

- **Independent noise** — `"factorised"`: the family's own `log_prob`
  evaluated pointwise, exactly; the terms sum to `log_prob`. Every family
  works, censoring included (a limit contributes its own Tobit term), and a
  user family carrying its own aligned data keeps X-2's excision contract
  per term: each single-sample call receives a full-length `retain`
  selecting exactly that sample.
- **A GP** — `"conditional_loo"`: the leave-one-out conditional terms
  `log N(y_i | μ_i^{-i}, σ_i^{2,-i})`, computed by
  `GPSolver.conditional_loo` from the same Cholesky the marginal likelihood
  forms (`DenseGP` implements the closed form; `QuasisepGP`'s O(N) recursion
  was deferred at W2.3 — `DEVELOPMENT_PLAN.md` §2, 2026-09-05 — and it
  refuses until one exists). These are what
  `arviz.loo`/`waic` consume, and they are a *different* decomposition: a
  GP joint has no per-observation factorisation, so the LOO terms
  deliberately do not sum to `log_prob`.

A latent combination is refused — its per-observation terms are conditional
on latent values that belong to inference.

```pycon
>>> terms = simple.pointwise_log_prob(model, data)
>>> terms.shape
(4,)
>>> bool(np.isclose(np.sum(terms), simple.log_prob(model, data)))
True
>>> gp_terms = flexible.pointwise_log_prob(model, data)
>>> bool(np.isclose(np.sum(gp_terms), flexible.log_prob(model, data)))
False

```

### The declarative spec: `to_spec`

Also landed at the freeze (`results.md` §15 R7, confirmed by the
consolidated serialisation review — `docs/design/serialisation_review.md`):
`Likelihood.to_spec()` returns the declarative description of the
composition as plain, JSON-able data — family name and class, noise-model
class, marginalisation, the parameters' spec, and for a GP the kernel spec
plus the solver's configuration. It is the one definition backends, the
conformance suite and provenance share, and it distinguishes exactly what
`ParameterSet.to_spec()` alone cannot:

```pycon
>>> matern = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(0.3, 2.0)))
>>> rbf = Likelihood(GaussianFamily(), GaussianProcessNoise(SquaredExponential(0.3, 2.0)))
>>> matern.parameters.to_spec() == rbf.parameters.to_spec()
True
>>> matern.to_spec()["kernel"]["family"], rbf.to_spec()["kernel"]["family"]
('matern32', 'squared_exponential')
>>> matern.to_spec()["solver"]
{'name': 'DenseGP', 'class': 'DenseGP', 'exact': True, 'config': {'jitter': 0.0}}

```

The spec describes the *declaration* only. Per-sample and bulk content — a
censoring declaration's code positions, a family's buffers — is
provenance's business: `ampere.results.describe_likelihood` composes this
mapping and adds the content fingerprints. That split (specs declare,
provenance fingerprints, storage carries values) is the review's one rule.

## 9. Censoring: what a limit is, and who consumes it (issue #11)

`results_schema.md` §8 draws the line and hands this side to this contract: a
masked sample carries zero information, whereas a 3σ non-detection says
something quite definite about the source. Discarding it throws that away;
recording it as a datum with a large error bar is simply the wrong likelihood.

A `Censoring` is a per-sample array of `LimitKind` codes, aligned index-by-index
with the observed container. It is declared **alongside** the container, not
inside it, because it is a statement about how an observation constrains a
model — a likelihood concern — and the same container may legitimately be
analysed with and without it.

```pycon
>>> list(LimitKind)
[<LimitKind.DETECTION: 0>, <LimitKind.UPPER_LIMIT: 1>, <LimitKind.LOWER_LIMIT: 2>]

```

`LimitKind` is an `IntEnum` so a declaration is an ordinary integer array — the
aligned-auxiliary-array hook §8 explicitly reserved:

```pycon
>>> photometry = PhotometricPoints(
...     ["WISE_W1", "WISE_W3", "WISE_W4"],
...     [3.4, 12.1, 22.2] * u.um,
...     [1.0, 0.40, 0.90] * u.Jy,
...     uncertainty=[0.05, 0.05, 0.30] * u.Jy,
...     extra_coords={"limit_kind": np.array([0, 0, 1])},
... )
>>> censoring = Censoring.from_extra_coord(photometry, "limit_kind")
>>> censoring
Censoring(DETECTION=2, UPPER_LIMIT=1)

```

A boolean flag array is not accepted directly, because "is a limit" does not
say *which* limit, and guessing is exactly the class of silent error this
project keeps designing out:

```pycon
>>> Censoring(np.array([False, False, True]))
Traceback (most recent call last):
    ...
ampere.core.exceptions.LikelihoodError: Censoring's kinds must be an integer array of LimitKind codes ...
>>> Censoring.upper_limits(np.array([False, False, True])).n_censored
1

```

**Which families consume it** is a class attribute, checked at composition:

```pycon
>>> sorted(name for name in list_families() if family_named(name).SUPPORTS_CENSORING)
['cauchy', 'gaussian', 'student_t']

```

The mathematics is the Tobit construction, written once for every family whose
standardised residual has a scipy distribution: a detection contributes
`log f(z) − log σ`; an upper limit contributes `log F(z)`, the probability that
the truth lies below the recorded value; a lower limit contributes
`log(1 − F(z))`.

```pycon
>>> censored = Likelihood(GaussianFamily(), IndependentNoise(), censoring=censoring)
>>> faint = photometry.with_values([1.0, 0.40, 0.10])
>>> bright = photometry.with_values([1.0, 0.40, 1.50])
>>> bool(censored.log_prob(faint, photometry) > censored.log_prob(bright, photometry))
True

```

Masking beats censoring: a masked sample contributes nothing regardless of its
limit kind. Declaring a limit on a masked sample is allowed and has no effect,
which is deliberate — masking a region for a test run should not require
editing the censoring array too.

That rule has to hold for the *declaration* as well as for the arithmetic, or a
problem that is analytic in fact would be refused a gradient-free engine.
Because whether a limit survives the mask is a property of the data,
`marginalisation` answers conservatively (it has seen no container) and
`marginalisation_for(observed)` answers for real data. W1.7, which has the
container, should use the second — and pass `observed=` to `check_engine` for
the same reason.

```pycon
>>> partly_masked = PhotometricPoints(
...     ["WISE_W1", "WISE_W3", "WISE_W4"],
...     [3.4, 12.1, 22.2] * u.um,
...     [1.0, 0.40, 0.90] * u.Jy,
...     uncertainty=[0.05, 0.05, 0.30] * u.Jy,
...     mask=np.array([False, False, True]),          # the only limit is masked
... )
>>> latent_by_declaration = Likelihood(
...     GaussianFamily(), GaussianProcessNoise(Matern32(0.3, 2.0)), censoring=censoring
... )
>>> latent_by_declaration.marginalisation
<Marginalisation.LATENT: 'latent'>
>>> latent_by_declaration.marginalisation_for(partly_masked)
<Marginalisation.ANALYTIC: 'analytic'>

```

**How the matrix algebra changes** — issue #11's own open question. It does not
close. A censored *multivariate* Gaussian likelihood is an orthant probability
of the multivariate normal, which has no closed form beyond a handful of
dimensions and needs Genz-style quadrature or a data-augmentation scheme over
truncated latent values. So censoring composed with a correlated noise model is
declared `LATENT`, and the analytic path refuses rather than approximating:

```pycon
>>> GaussianFamily().marginalisation_with(
...     GaussianProcessNoise(Matern32(0.3, 2.0)), censoring
... )
<Marginalisation.LATENT: 'latent'>

```

No family implements the truncated-latent form, so — by §4's rule that
declaring `LATENT` is not the same as implementing it — the combination is
refused as soon as data confirm a limit actually survives the mask:

```pycon
>>> observation = PhotometricPoints(
...     ["WISE_W1", "WISE_W3", "WISE_W4"],
...     [3.4, 12.1, 22.2] * u.um,
...     [1.0, 0.40, 0.90] * u.Jy,
...     uncertainty=[0.05, 0.05, 0.30] * u.Jy,
... )
>>> latent_by_declaration.check_alignment(observation.with_values([1.0, 0.4, 0.1]), observation)
Traceback (most recent call last):
    ...
ampere.core.exceptions.LikelihoodError: the gaussian family with a censoring declaration on correlated noise is a latent-variable model ...

```

Staged: the analytic (uncorrelated) case is implemented here, which is the case
issue #11 calls "fairly straightforward in other settings"; the correlated case
is declared and deferred to the modern backends along with every other latent
combination.

## 10. The latent-GP declaration, and whether it scales to 10⁵

`parameters.md` §13 asked this contract to answer, explicitly rather than by
assumption, whether an array-valued `HierarchicalPrior` scales as the latent-GP
declaration for `N ~ 10⁵`.

**The answer is that `HierarchicalPrior` is not the right declaration at any
N**, for a reason that has nothing to do with size. `parameters.md` §12.1
states that an array-valued parameter's prior is i.i.d. across its elements,
and a GP prior is precisely *not* i.i.d. — its entire content is the
correlation between elements. Declaring `f` as a hierarchical Normal
referencing the GP hyperparameters would describe white noise with a fitted
variance, and would do so silently.

The declaration is instead the standard non-centred (**whitened**)
parameterisation, which needs no extension to `parameters.md` at all:

- `z` is one array-valued parameter of shape `(N,)` with an i.i.d.
  standard-normal prior and an `Identity` bijection — exactly what §12.1
  already supports, and exactly the geometry a PPL wants for HMC;
- the covariance enters through `f = L(θ) z`, a deterministic transform owned by
  the `GPSolver`, where it can be a Cholesky factor (`DenseGP`) or a state-space
  recursion (`QuasisepGP`) without the declaration changing.

One honest caveat, since fixing the declaration fixes the geometry with it: the
non-centred form is the right default when the data constrain the latent
function weakly relative to its prior, which is the regime a misspecification
GP is usually in, but the *centred* form (`f` sampled directly, with the
hyperparameters entering its prior) has the better geometry when the data
constrain it strongly. Switching between them is a reparameterisation, not a
model change, so a backend may offer it; this contract fixes the default and
the transform, not the only possible parameterisation.

```pycon
>>> declaration = latent_poisson.latent_declaration(counts.n_samples)
>>> parameter = declaration.parameter
>>> parameter.shape, parameter.is_free, parameter.is_hierarchical, parameter.references
((4,), True, False, ())
>>> parameter.unconstraining_bijection()
Identity()
>>> declaration.as_parameter_set().free_size
4

```

Asking an analytic likelihood for one raises, because quietly returning an
empty declaration would let a caller build a sampler dimension that does
nothing:

```pycon
>>> flexible.latent_declaration(4)
Traceback (most recent call last):
    ...
ampere.core.exceptions.LikelihoodError: the gaussian family with a GaussianProcessNoise noise model marginalises analytically, so it declares no latent values. ...

```

The whitening transform is the mathematical content of the declaration, and it
is testable: `L L^T == K`.

```pycon
>>> x = np.array([[0.0], [1.0], [2.5]])
>>> values = {"amplitude": 0.6, "length_scale": 2.0}
>>> columns = [
...     DenseGP().latent_transform(Matern32(0.6, 2.0), x, e, values)
...     for e in np.eye(3)
... ]
>>> lower = np.column_stack(columns)
>>> reference = Matern32(0.6, 2.0).matrix(x, x, values)
>>> bool(np.allclose(lower @ lower.T, reference, atol=1e-8))
True

```

`latent=` is the **whitened** `z`, and `GaussianProcessNoise.noise_params`
applies the transform on the retained coordinates before the family ever sees
it (W2.14): the family reads `noise.latent` as `f`, so `f = L(θ) z` is imposed
once, in the one place that holds the solver, the kernel and the coordinates.
That is why the value below moves when the hyperparameters move — before
W2.14 it did not, and a latent fit sampled the amplitude and the length scale
against a flat likelihood.

```pycon
>>> rate = counts.with_values([3.5, 6.0, 2.5, 8.0])
>>> latent_poisson.log_prob(rate, counts)
Traceback (most recent call last):
    ...
ampere.core.exceptions.LikelihoodError: a Poisson likelihood with a GaussianProcessNoise model is a latent-variable model: ...
>>> z = np.array([0.1, -0.2, 0.05, 0.0])
>>> points = np.asarray(counts.axes[0].values).reshape(-1, 1)
>>> f = DenseGP().latent_transform(
...     Matern32(0.3, 2.0), points, z, {"amplitude": 0.3, "length_scale": 2.0}
... )
>>> bool(np.isclose(
...     latent_poisson.log_prob(rate, counts, latent=z),
...     np.sum(st.poisson.logpmf(counts.values, rate.values * np.exp(f))),
... ))
True

```

The same `z` under a wider kernel is a different number, which is the whole
content of the fix:

```pycon
>>> wider = Likelihood(PoissonFamily(), GaussianProcessNoise(Matern32(3.0, 2.0)))
>>> bool(
...     latent_poisson.log_prob(rate, counts, latent=z)
...     != wider.log_prob(rate, counts, latent=z)
... )
True

```

**On scale**, `N ~ 10⁵` is fine as a *declaration* — one `Parameter`, one shape
tuple, one `PriorSpec` — and pack/unpack is a single reshape. Three things
downstream are not fine, and are recorded as obligations in §16 rather than
left to be discovered:

1. `ParameterSet.free_labels()` materialises one string per element. At 10⁵
   that is a per-call allocation of ~10⁵ strings, and W1.8 must give ArviZ a
   *dimension and coordinate* for the latent block, not 10⁵ scalar names.
2. The sampler dimension becomes `N + k`, which no gradient-free engine can
   address. That is not a defect: it is `check_engine`'s reason to exist.
3. Forming `L` densely is O(N³) in time and O(N²) in memory. The latent path at
   that scale therefore requires `QuasisepGP` — the latent-GP promise is a
   modern-backend *and* an O(N)-solver promise, not one or the other.

`Plate` is not the right construct either: it expands to per-member parameters
for a *population*, and shares `free_labels`'s per-element cost without
supplying the correlation.

## 11. What ampere does instead of Starfish's local kernels

`prior_art.md` lesson S3 and Tension 4 both ask this contract to say the
following outright, so a future reader does not propose Starfish-style local
kernels as if they were a straightforward extension.

Beyond its global Matérn kernel, Starfish adds one local Gaussian "bump" kernel
per identified problem line, each with three hyperparameters, and the *number*
of them is itself chosen during fitting: `4·N_loc + 2` covariance
hyperparameters with `N_loc` not fixed. That is a trans-dimensional parameter
space, handled with a blocked Gibbs/Metropolis-Hastings sampler.

Ampere's §4.1 parameter contract is **statically declared and fixed-shape** by
design, because that is what lowers cleanly to `torch.distributions`/numpyro
plates (W1.9) and what NUTS and VI require. A Starfish-style local-kernel model
is not expressible in it without a reversible-jump/model-selection layer that
nothing in the current plan scopes.

**This is a deliberate, different answer to the same scientific problem, not an
oversight.** Ampere's single stationary Matérn kernel over the whole domain,
plus the GP-localisation diagnostics of §4.8 — the posteriors on amplitude and
length-scale, and the conditioned GP mean, which is signed and therefore shows
the *direction* of the local deficiency — is the substitute for explicit
per-line local kernels. The trade is: a global model that tells you where it is
locally wrong, rather than a model with explicit local components.

```pycon
>>> residual_data = Spectrum(
...     np.linspace(1.0, 10.0, 40) * u.um,
...     np.sin(np.linspace(1.0, 10.0, 40)) * u.Jy,
...     uncertainty=np.full(40, 0.01) * u.Jy,
... )
>>> localiser = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(2.0, 1.0)))
>>> conditioned = localiser.conditional(residual_data.with_values(np.zeros(40)), residual_data)
>>> bool(conditioned.mean.min() < 0.0 < conditioned.mean.max())   # signed, not |magnitude|
True
>>> bool(np.allclose(conditioned.mean, residual_data.values, atol=0.05))
True

```

A caveat is attached to that output, and it has two distinct sources that
should not be run together.

`prior_art.md` lesson S4 records Starfish's own finding: its global-kernel
amplitude and its explicit local kernels trade off against each other on real
data, so "a large fitted GP amplitude with a short length-scale can mean either
'genuinely misspecified locally' or 'the kernel's smooth global component is
under-amplitude and compensating locally'". That is a *global-versus-local*
degeneracy, and ampere's single-kernel design does not have it in that form —
there are no local components to trade against.

The caveat that does apply here, and which W1.12's family C carries as a
mandatory part of its output, is the more general one: a large fitted amplitude
localises a deficiency but does not say *why*. Model error, an underestimated
noise budget and a genuinely correlated astrophysical process are all
consistent with the same posterior. That is this contract's own statement, not
a quotation from S4, and it is why the localisation output is a pointer to
where to look rather than a diagnosis.

## 12. Position on the irregular-coordinate whiteness test

W1.12 §3.2 deliberately left open whether the post-fit whiteness test should
resample residuals onto a nominal regular grid or use a separation-binned
statistic native to irregular spacing, and asked W1.6 to decide with a
solver-level view of what is cheap. This contract decides: **the
separation-binned statistic**, and recommends against resampling.

Three reasons, two of them specific to this contract:

1. **Resampling manufactures the artefact the test looks for.** Interpolating
   residuals onto a grid applies a smoothing operator, which imprints
   correlation between neighbouring output samples that was not in the input.
   A whiteness test on resampled residuals is biased towards rejecting the null
   by construction — a false "you need the GP here" signal.
2. **The statistic and the model it motivates should live in the same
   coordinates.** The kernel this test exists to motivate is a function of
   separation alone: `k(r)`. A separation-binned empirical autocovariance of
   standardised residuals is the direct empirical estimate of the very object
   `Matern32` parameterises, so the diagnostic output and the fitted
   length-scale are comparable quantities rather than cousins.
3. **The cost objection is answerable from the solver side.** Binning does not
   need all `N²` separations: for lags out to `r_max` it needs only pairs
   within `r_max`, which is a windowed sweep over an ordered axis — `O(N·m)`
   with `m` the mean neighbour count inside `r_max`. That is the *same*
   neighbourhood structure `WindowedSparseGP` needs, so the two share
   machinery rather than each inventing their own.

One honest caveat to carry into Phase 2: bins share samples, so the classical
Ljung-Box asymptotic χ² degrees-of-freedom do not apply unchanged. The
reference implementation should calibrate the null by permutation or parametric
bootstrap of the standardised residuals rather than assuming the textbook
distribution. That is more work than calling `statsmodels`, and it is the price
of not assuming a grid that `results_schema.md` spent a whole contract refusing
to assume.

## 13. Worked example: the flexible likelihood earning its keep

A model that gets the continuum right and a feature wrong. Under an independent
Gaussian likelihood the residual structure is charged to the physical
parameters; under the flexible likelihood the GP absorbs it.

```pycon
>>> wavelength = np.linspace(5.0, 25.0, 60)
>>> continuum = 10.0 * (wavelength / 10.0) ** -1.5
>>> feature = 1.2 * np.exp(-0.5 * ((wavelength - 15.0) / 1.5) ** 2)
>>> observation = Spectrum(
...     wavelength * u.um,
...     (continuum + feature) * u.Jy,
...     uncertainty=np.full(60, 0.05) * u.Jy,
... )
>>> smooth_model = observation.with_values(continuum)      # misses the feature entirely

```

The independent-Gaussian likelihood is dominated by the unmodelled feature:

```pycon
>>> rigid = Likelihood(GaussianFamily(), IndependentNoise())
>>> bool(rigid.log_prob(smooth_model, observation) < -100.0)
True

```

The flexible likelihood, with a kernel whose length-scale is comparable to the
feature width, recognises the residual as correlated structure rather than 60
independent surprises:

```pycon
>>> flexible_fit = Likelihood(
...     GaussianFamily(), GaussianProcessNoise(Matern32(1.0, 2.0))
... )
>>> bool(flexible_fit.log_prob(smooth_model, observation)
...      > rigid.log_prob(smooth_model, observation))
True

```

And it says **where**, which is the whole point of §4.8's localisation family:

```pycon
>>> local = flexible_fit.conditional(smooth_model, observation)
>>> bool(abs(wavelength[np.argmax(local.mean)] - 15.0) < 1.5)
True

```

Composition is checked once, loudly, before any of that runs:

```pycon
>>> flexible_fit.check_alignment(smooth_model, observation)
>>> flexible_fit.check_alignment(smooth_model.to_unit(u.mJy), observation)
Traceback (most recent call last):
    ...
ampere.core.exceptions.LikelihoodError: the predicted Spectrum is in mJy but the observed one is in Jy. Convert once at composition time with .to_unit(...); units never reach the hot loop ...

```

Note what is *not* checked per evaluation: `check_alignment` compares
coordinates, units and kinds in O(N), and is a composition-time obligation
W1.7's `Dataset` discharges once. `log_prob` re-checks only shapes, mirroring
`results_schema.md` §10's compile-once/evaluate-many split.

## 14. Decisions and their reasoning

| Decision | Reasoning |
|---|---|
| Family interface is one method, `log_prob(predicted, observed, noise)` | `DEVELOPMENT_PLAN.md` §4.4's signature and `prior_art.md` 3M3's minimalism; everything declarative is a class attribute, so declaring costs nothing per evaluation |
| Family, noise model and solver are three objects, not one | `prior_art.md` Tension 3: 3ML's fused plugin buys extensibility at the cost of all reuse. The risk that a real instrument does not decompose is W1.11's stress test |
| `Marginalisation` is a property of the *combination* | Poisson alone is analytic; Poisson + GP is not. Declaring it on either piece separately would be wrong for half the pairs |
| `check_engine` refuses a gradient-free engine on a latent likelihood | §4.4's "effectively a modern-backend capability", enforced. A silent downgrade is a fit that runs and is wrong |
| `check_engine` takes keyword flags, not a capability object | W1.7 owns the capability record (§4.5); the two contracts should not have to agree on a class before either is written |
| Masks arrive as `weights()`; masked samples are **excised** | `results_schema.md` §16 asked for one convention. Zero weight and infinite uncertainty are equivalent statements about a chi-square term, which is what §7 claims for them — but not about a normalised log-density, whose `−log σ` (and, under a GP, log-determinant) diverges too. Excision is the only rule delivering §7's "exactly zero information", for *any* noise model, and a 0/1 weight on a whole term is excision written arithmetically |
| The effective mask is the union of predicted and observed | A sample needs both sides valid; this also composes with W1.5's obligation to propagate masks through transformations |
| A family opts in to the latent path via `CONSUMES_LATENT_GP`, default `False` | Declaring `LATENT` and implementing it are different things, and a family that declares without implementing returns the *uncorrelated* likelihood — so every GP hyperparameter and every latent value an engine sampled leaves the log-probability untouched. That is the "runs, converges, and is wrong" failure the contract exists to prevent, so it must not be reachable by composing two objects that each look reasonable |
| Family-and-noise latency is refused at construction; censoring-induced latency at `check_alignment` | Whether a family plus a noise model is latent is a fact about the pair. Whether *censoring* makes it latent depends on whether a limit survives the mask, which needs the data — so it is refused at the earliest point the question can honestly be asked |
| `marginalisation` is conservative; `marginalisation_for(observed)` is data-aware | §9's "masking beats censoring" has to hold for the declaration too, or an analytic problem gets refused a gradient-free engine. A property that has seen no container cannot know, so both answers exist and W1.7 uses the second |
| A bare 1-D array of coordinates means *n* 1-D points, never one *n*-D point | `np.atleast_2d` reads it the other way and then broadcasts silently through the kernel, turning a 101-point diagnostic grid into a one-point answer with no error. `_as_points` fixes the reading and raises on anything ambiguous |
| A fully masked pair returns `0.0` | The correct limit of the same rule, not a special case |
| The conditioned GP mean is evaluated on the *full* axis, masked points included | "What would the GP have said here?" is exactly the question a user asks about a region they excluded |
| Non-finite values in retained samples raise | `results_schema.md` §7 and `prior_art.md` R1: no sentinel NaNs threaded through arithmetic |
| GP hyperparameters are ordinary parameters with `Log` bijections | `parameters.md` §13's instruction, confirmed. Tying, fixing, plates, serialisation and lowering then work unchanged |
| Kernel amplitude is a standard deviation (`k(0) = a²`) | celerite2's `Matern32Term(sigma=…)` convention, so a prior on it is a prior in the data's own units |
| Kernels declared as family name + hyperparameter names | The minimum W1.9 needs and the maximum that translates across celerite2 / tinygp / GPyTorch |
| Hyperparameter units must match the data's exactly, not merely be convertible | Priors are numeric in the declared unit and rescaling a distribution correctly is family-specific — the same ruling `parameters.md` §8 makes for tying |
| Matérn-3/2 is the default kernel | `DEVELOPMENT_PLAN.md` §2; better residual representation *and* exact quasiseparability. `prior_art.md` S1 is independent published corroboration |
| Squared exponential retained but flagged non-quasiseparable | M2's comparison against legacy needs it; the asymmetry is the concrete argument for the default |
| Solver is a strategy with composition-time `check_compatible` | Scaling changes without science code changing; and a mismatch fails when the problem is built, not a thousand samples in |
| Declarative incompatibilities are reported before "not implemented yet" | A permanent fact about the choice is more useful than a temporary fact about the schedule |
| `WindowedSparseGP` gets a named slot *and* a "why not default" argument | `prior_art.md` S2 asked for one; a Starfish-aware reviewer deserves both. Exactness, a posterior-varying cutoff, and truncated tails that are themselves diagnostic output are the three reasons |
| `DenseGP` defaults to `jitter=0.0` | A covariance that will not factorise is a fact about the model. The error message names the cause and the remedy instead of hiding it behind a numerical fudge |
| A non-PD covariance raises rather than returning `−inf` | A silent `−inf` hides a mis-specified kernel. W1.7's failure-signalling layer (§4.5) is where a sampler-facing wrapper converts it, with a recorded reason |
| Latent values are declared **whitened**, not as a `HierarchicalPrior` | A GP prior is not i.i.d., and `parameters.md` §12.1 makes array-valued priors i.i.d. A hierarchical Normal would silently describe white noise |
| `Likelihood.parameters` is one flat `ParameterSet`, not a nested merge | `parameters.md` §12.4: `merge` is not associative, and W1.7 must merge every dataset in one call |
| A family/noise parameter-name collision raises | Silent qualification would make the flat namespace a lie; the fix (rename) is one call |
| `Censoring` lives beside the container, not inside it | It is a statement about how an observation constrains a model, and the same container may be analysed with and without it. `results_schema.md` §8 deliberately left `mask` strictly binary for this |
| `LimitKind` is an `IntEnum` | So a declaration is an ordinary aligned integer array — precisely §8's reserved hook |
| A boolean flag array is refused; `upper_limits`/`lower_limits` are explicit | "Is a limit" does not say which limit, and guessing is the error class this project keeps designing out |
| Censoring + correlated noise is `LATENT` | Issue #11's own open question, answered: it is a multivariate-normal orthant probability with no closed form beyond a few dimensions |
| Masking beats censoring on the same sample | Masking a region for a test run should not require editing the censoring array too |
| Complex data are the circular complex Gaussian only | `results_schema.md` §16: the container's real σ encodes exactly that. Non-circular noise supplies its own 2×2 structure and is not a container concern |
| `complex_gaussian` + GP is `ANALYTIC` — circular meaning fixed, implementation staged | Ruled 2026-09-03 (§17 Q6). The circular complex GP marginalises in closed form exactly as the real Gaussian does, and fixing the declaration now unblocks the flexible likelihood on the Phase-4 proof modality. `GP_ANALYTIC_IMPLEMENTED = False` keeps the pair a composition-time refusal until Phase 4 lands the closed form — declared-but-staged, never silently different. *(Amended W4.8: implemented at W4.2 — `GP_ANALYTIC_IMPLEMENTED` is now `True`; see §4's "declared analytic, implemented at W4.2" subsection for the closed form.)* |
| Rice takes amplitudes from an `Amplitude` chain step; von Mises takes `κ = 1/σ²` per sample | Ruled 2026-09-03 (§17 Q3/Q4). The model predicts what it physically produces — the complex value — and projection is the instrument chain's job; the concentration comes from the container's own uncertainties, exact in the small-σ limit where closure-phase practice lives. Families themselves are Phase 4's. *(Amended W4.8: `VonMisesFamily` implemented at W4.1, exactly to this interface. `RiceFamily` remains declared, not implemented — `Amplitude` landed at W4.1, but nothing consumes it as a Rician mean yet, and Phase 4 did not schedule it.)* |
| `check_alignment` is composition-time; `log_prob` re-checks only shapes | O(N) coordinate comparison is right once and wrong per evaluation — `results_schema.md` §10's split |
| A noise model receives the prediction as well as the observation | Ruled 2026-09-03 (X-1). A noise whose magnitude depends on the model — a fractional model uncertainty, an analytically marginalised multiplicative calibration systematic, a model-variance weighting of counts — is a `NoiseModel`, not a family. Without the `predicted` argument the only way to express one is to re-implement the sampling distribution, which welds noise to family, cannot be reused, and cannot reach the GP path: exactly the monolithic collapse `prior_art.md` Tension 3 warns against |

## 15. Deliberate limitations of v1.6

Each is a decision, not an oversight. Each has an extension point.

1. ~~**One kernel per noise model.**~~ **Lifted at W4.5.** `Sum` and
   `Product` compose kernels, their `spec()` is a tree of their children's, and
   a `Sum` of quasiseparable terms lowers to the O(N) path by concatenating
   their generators (§6, "Kernel algebra"). A `Product` stays on `DenseGP` and
   says so by name. The original note that this is *not* the Starfish
   local-kernel case (§11) stands: a fixed sum of global kernels is
   fixed-shape, which is why it fits the parameter contract at all.
2. **Isotropic kernels only** — *half lifted at W4.5*. Separation is still
   Euclidean in the coordinates a kernel sees and there is still one
   `length_scale`, so a genuinely anisotropic 2D kernel (two length scales in
   one plane) is still not expressible. What W4.5 adds is the `axes` selector
   (§7): a kernel acts on a named subset of its container's axes, the
   single-unit rule applies to that subset, and `Product` composes kernels on
   disjoint subsets — so a container with axes in different units is now
   *usable* rather than merely refused, and the two-length-scale case that
   motivated the limitation (a patch in (u, v) times a band in wavelength) has
   an expression. Mixed units **within one kernel's selection** are still
   refused outright rather than being silently wrong.
3. ~~**Only Matérn-3/2 and the squared exponential ship.**~~ **Lifted at
   W4.5**, which added `Matern12`, `Matern52`, `SHO`, `RotationTerm` and
   `SpectralMixture`, and opened the translation table as
   `register_quasiseparable_term` so a user kernel can reach the O(N) path
   too (§6). W2.3's rule is unchanged and now enforced through that registry:
   a declared `QUASISEPARABLE` with no registered representation is a loud
   refusal, never a silent approximation.

   **This clause's last sentence was wrong and is withdrawn.** Matérn-5/2 *is*
   exactly quasiseparable — rank 3, by the same argument W2.3 used for
   Matérn-3/2. It has no exact form in the celerite *basis*, which is what the
   sentence was really about; the solver factorises the wider *semiseparable*
   form, in which a degree-2 polynomial in `tₙ − tₘ` is a rank-3 bilinear form
   in `(1, t, t²)`. §6 carries the algebra. Nothing shipped on the strength of
   the wrong claim: no Matérn-5/2 existed to be flagged.
4. **`GaussianProcessNoise` only supports `Layout.POINTS`.** A gridded 2D GP is
   the SVGP/SKI/Vecchia slots' business (Phase 5); `IndependentNoise` works on
   any layout.
5. **Rice and von Mises are declared, not implemented.** Both raise on
   composition; the implementations are Phase 4's. Their parameterisations
   were fixed by the 2026-09-03 rulings (§17 Q3/Q4): the model predicts the
   complex value and an `Amplitude` chain step takes the modulus for Rice;
   `κ = 1/σ²` per sample for von Mises. *(Amended W4.8: `VonMisesFamily` was
   implemented at W4.1, to exactly this interface, with the interferometry
   modality that needed it — see the class docstring. `RiceFamily` remains
   declared, not implemented; W4.1 checked it as a third route for
   interferometric amplitudes and did not schedule it further.)*
5a. **`PoissonFamily` is the only family that consumes the latent path.**
   Student-t and Cauchy *declare* `LATENT` under a GP and neither implements
   it, so composing either with `GaussianProcessNoise` is refused (§4). That
   is a real capability gap — a heavy-tailed flexible likelihood is a
   reasonable thing to want — and it is a refusal rather than a wrong answer
   only because `CONSUMES_LATENT_GP` exists. Implementing one is a
   `log_prob` that reads `noise.latent` plus flipping the flag. (The complex
   Gaussian left this list at the freeze: under a GP it now declares a
   *staged* `ANALYTIC` — §4 — refused until Phase 4 implements the circular
   closed form.) *(Amended W4.8: Phase 4 did — W4.2 implements the circular
   complex GP in closed form, so `complex_gaussian` + `GaussianProcessNoise`
   is no longer staged; it marginalises analytically like every other
   `ANALYTIC` pair. Student-t and Cauchy are unaffected and remain refused
   under a GP.)*
6. **The latent path has no inference.** *(Amended W2.14, 2026-09-08.)* As
   written, this limitation also described a **defect**: `latent_transform`
   was "the transform" and nothing on the scoring path applied it, so
   `Likelihood.log_prob` read the whitened `z` as `f` and a latent
   log-likelihood was exactly flat in the kernel hyperparameters. That is
   fixed — `GaussianProcessNoise.noise_params` applies
   `f = solver.latent_transform(kernel, coordinates, z, values)` on the
   retained block, so `noise.latent` *is* `f` as every family already read it
   (§10; `DEVELOPMENT_PLAN.md` §2, 2026-09-08, a §4.4 clarification). What
   remains a limitation is only the inference: `ampere.core` supplies the
   declaration, the transform and the conditional density, and the *engines*
   that can sample an `N + k`-dimensional latent posterior are the torch and
   jax rungs' (both now lower the corrected path, and the conformance battery
   holds them to it). `DenseGP.latent_transform` still exists so the
   declaration is testable at small N and so a Laplace-type reference
   fallback (§4.4's "possible later fallback") is not blocked.
7. **Censoring is interval-free.** `LimitKind` has no `INTERVAL` member for a
   datum known only to lie between two values. Adding one is a code plus a
   second value array; nothing in the design forbids it, and nothing in the
   target scope needed it yet.
8. **No per-sample family mixing.** One family per `Likelihood`. A dataset that
   is Poisson in one region and Gaussian in another is two channels and two
   datasets (W1.7), not one likelihood with a per-sample switch.
9. **No log-likelihood *per sample*.** `log_prob` returns a scalar. W1.8's
   InferenceData requirement is per-*observation* log-likelihood, which for the
   GP case is not well defined anyway (the samples are not independent) — see
   §16's obligation on W1.8. *(Amended at the freeze: the named decomposition
   — pointwise terms for independent noise, leave-one-out conditionals for a
   GP — is now available on request through `pointwise_log_prob`, ruled
   2026-09-03, `results.md` §15 R2; it is not stored by default and
   `log_prob` still returns a scalar.)*
10. **The noise model sees the prediction, but not the model's parameters.**
   `predicted` (§5) carries the retained predicted *values* only; a noise
   term that depends on a physical parameter beyond those the likelihood
   declares must be tied to a likelihood parameter or expressed as a
   transformation. Retained deliberately (X-1's ruling keeps it): widening
   the argument to model internals would re-fuse the pieces this contract
   exists to separate.
11. **No cross-channel (vector-valued) correlated noise model.** A noise
   model belongs to one `Likelihood`, which belongs to one dataset, so a
   single correlated process over 2-vectors — astrometric (RA, Dec)
   residuals perturbed together, Stokes Q/U leakage, a calibration
   systematic shared across bands — is not expressible; two scalar GPs
   with tied hyperparameters is the nearest approximation and is a
   different model (it cannot express the cross-covariance). Recorded from
   the W1.11 astrometric sketch at the freeze; the extension point is a
   `JointGP` solver-strategy slot spanning datasets, deferred to
   **Phase 5** with the other advanced strategies.

## 16. What this contract hands to the specs downstream

- **W1.5 (Transformation / Instrument)** — the effective mask is the union of
  the predicted and observed containers' masks, so a transformation that drops
  a mask silently *adds* data to a fit. `results_schema.md` §16 already places
  the propagation obligation on you; this contract is the consumer that
  depends on it. Also: a resampling transformation must not produce coordinates
  that differ from the observed container's, because `check_alignment` compares
  axes for equality — negotiation should target the observed grid. (Its
  general form, landed at the freeze from W1.11 gap I-2: any step
  reproducing the observed coordinates takes them *from* the observed
  container, never recomputes them — `transformations.md` §10.)
- **W1.7 (Dataset / FittingProblem)** — `Dataset` should call
  `Likelihood.check_alignment(predicted_template, observed)` once at
  construction, and `Likelihood.check_engine(differentiable=…, engine=…,
  observed=…)` when an engine is selected. **Pass `observed`**: without it the
  censoring half of the marginalisation is answered conservatively, and a
  problem whose only limits are masked would be refused a gradient-free engine
  it can perfectly well use (§9). `check_alignment` is also where a latent
  combination whose family does not implement the latent path is refused, so
  calling it is not optional. `Likelihood.parameters` is a flat `ParameterSet` ready
  to be one component of a single `ParameterSet.merge` across datasets — do not
  merge it again on its own. For a latent combination, `latent_declaration(n)`
  returns the extra parameter to include in that same merge; `n` is the number
  of *retained* samples, which the dataset knows and this contract does not.
  Finally, this contract raises `LikelihoodError` on a non-positive-definite
  covariance rather than returning `−inf`: §4.5's failure-signalling layer is
  yours, and converting that exception to `−inf` *with a recorded reason* is the
  right place to do it.
- **W1.8 (Results)** — two things. (a) The latent block cannot be labelled with
  `free_labels()` at scale: 10⁵ scalar names is the wrong ArviZ representation,
  and a named dimension with a coordinate is the right one. (b) **"Per-sample
  `log_likelihood`" is ambiguous and the two readings differ for this
  contract.** `DEVELOPMENT_PLAN.md` §4.6 asks that "every run stores per-sample
  `log_likelihood` and `log_prior`", and justifies it by population-level
  importance reweighting (design horizon (b)) — which needs the *scalar joint*
  log-likelihood per posterior draw, and which this contract's `log_prob`
  returns directly. ArviZ's own `log_likelihood` group, however, is
  conventionally per-*observation*, because that is what LOO and WAIC consume;
  and W1.12 §7 reads §4.6 the second way when it asks for "signed per-point
  residuals, not just per-observation log-likelihood". Both are legitimate
  readings of the same sentence. It matters here because **a GP likelihood has
  no well-defined per-observation decomposition**: the samples are not
  independent, so the joint value does not factorise. If W1.8 wants the ArviZ
  convention it must choose a decomposition and name it — the leave-one-out
  conditional terms are the standard choice and are computable from the same
  Cholesky this contract already forms. W1.13 should reconcile the two readings
  in the plan's own wording. *(Done at the freeze: `results.md` §6 named both
  decompositions, `pointwise_log_prob` computes them (§8, ruled 2026-09-03
  R2), and `DEVELOPMENT_PLAN.md` §4.6's wording now states the per-draw
  scalar as what every run stores and the per-observation terms as available
  on request.)*
- **W1.9 (Lowering)** — the declaration forms needing a lowering row are:
  `KernelSpec` per family (`matern32` → `celerite2.terms.Matern32Term` /
  `tinygp.kernels.quasisep.Matern32` / a GPyTorch equivalent;
  `squared_exponential` → an RBF kernel with no quasiseparable path), the
  `GPSolver` strategies, the whitened latent parameter (an i.i.d.
  `Normal(0,1)` sample site of shape `(N,)` inside no plate, plus a
  deterministic `f = L z`), and the float64 requirement, which is
  non-negotiable here and is `architecture.md` §5's policy.
- **W1.10 (Conformance suite)** — candidate rows already implemented and
  testable here: DenseGP against `scipy.stats.multivariate_normal.logpdf`; the
  zero-amplitude reduction to the i.i.d. Gaussian; mask excision equalling
  deletion of the sample; DenseGP↔QuasisepGP agreement on Matérn-3/2 once the
  latter exists (§4.6 names it — live since W2.3); `L L^T == K` for the whitening transform; the
  Tobit censored likelihood against `scipy.stats.norm.logcdf`; and the
  marginalisation declaration for every family × noise-model pair. Added at
  the freeze (X-1's three rows, `awkward_instrument.md` §6 point 9): the
  σ_eff of a fractional noise against the manual quadrature at fixed θ;
  `simulate(observe=True)` draw variance growing with the prediction as
  `(f·μ)²`; and the GP composition against a `DenseGP` evaluation with a
  manually precomputed diagonal.
- **W1.11 (Modality sketches)** — `ComplexGaussianFamily` plus
  `IndependentNoise` is the visibility likelihood; check whether closure phases
  need `VonMisesFamily` before the freeze, since it is currently declared-only.
  The X-ray sketch should check that `PoissonFamily` plus a response matrix in
  the instrument chain is expressible without a per-sample exposure concept
  this contract lacks. *(Confirmed — `awkward_instrument.md` §2: a `Spectrum`
  with a keV axis plus a response-matrix `Transformation` is expressible, and
  the absence of a per-sample exposure concept is correct rather than a gap —
  the exposure folds into the response matrix, `transformations.md` §10.)* And `prior_art.md` Tension 3 asks for a deliberately
  awkward instrument as a stress test of the §4.3/§4.4 split — this contract
  is one half of what that test exercises.
- **W1.12 (Diagnostics)** — the three things you asked for: (a) GP
  hyperparameters are confirmed as ordinary parameters on a `Parameterised`
  noise model with `Log` bijections (§6); (b) the irregular-coordinate
  whiteness test is decided in favour of the separation-binned statistic, with
  the null calibrated by permutation rather than the classical χ² (§12); (c)
  the Starfish-substitute framing is on record here as well as in your §4.1
  (§11). Family C's input is `Likelihood.conditional`, which returns a signed
  mean and a variance on the full coordinate axis, masked samples included.

## 17. Open questions for review

**Ruled by Peter, 2026-09-02**: question 1 — a **`strict` mode** (or an
equivalent toggle) is the approach. The engine-facing path is non-strict:
`−inf` with a recorded reason — forced in the end by jax, whose traced hot
loops cannot use exception control flow, so the non-raising path must exist
by Phase 2 regardless. Strict raising remains the behaviour for direct use
and debugging, and the intended workflow is: a sampling run surfaces the
recorded-reason warnings, and the user re-runs with strict on to get the
raise at the offending draw with full details. W1.7 owns the conversion
layer and the reason recording (§16's obligation, extended the same day to
`PoissonFamily`'s non-positive-rate raise); the solver-level flag on
`DenseGP` follows once W1.7's recording mechanism exists, so try/except can
be replaced without contract change. Question 2 — **the `Dataset` resolves
the effective mask once at construction**, not `Likelihood` on every call.
Consequences: the effective mask becomes a declared evaluation-time
invariant (a transformation whose output mask depends on parameter values is
unsupported in a fitting problem, and W1.7 should check this loudly), which
also makes explicit the invariance that `latent_declaration(n)` — whose `n`
is fixed at composition — was already assuming. This contract's own
mechanics (§8's `weights()` product and excision) are unchanged: a
pre-resolved pair simply makes the internal union the identity. **Also ruled
2026-09-02** (W1.7's R3): `LikelihoodFamily` gained the optional generative
half, `sample(predicted, noise, rng)` — §3 has the contract; the default is
a specific refusal and `GaussianFamily` implements it. **And later the same
day, three W1.11 amendments were approved and landed**: `check_alignment`
compares value dtype kinds (gap I-1 — a complex prediction can no longer be
silently fitted against real amplitudes); families gained the
composition-time `check_observed(observed)` hook (gap I-5 — called by
`check_alignment` on the observed container, where `PoissonFamily`'s
integrality test now lives per gap X-3, masked samples exempt); and
`NoiseParams` carries `retain` (gap X-2 — §8), the caller's inclusion
indicator over the full containers, so a family with its own aligned
per-sample data excises it the same way. Questions
3, 4, 6 and 7 carry recommendations from the W1.11 interferometry sketch
(`docs/design/modalities/interferometry.md` §§5–7 and §11: the model
predicts the complex value with an `Amplitude` step taking the modulus;
`κ = 1/σ²` per sample; declare the circular complex GP `ANALYTIC`;
`extra_coords` units eventually, not for the freeze) and await Peter's
ruling; questions 5 and 8 remain open as written.

**Ruled by Peter, 2026-09-03**: questions 3, 4 and 6 — the sketch's
recommendations are accepted as written: the model predicts the complex
value and an `Amplitude` step takes the modulus; `κ = 1/σ²` per sample;
and the circular (equal-component, zero-pseudo-covariance) complex GP is
the fixed meaning under which `complex_gaussian` + `GaussianProcessNoise`
is declared `ANALYTIC` — unblocking the flexible likelihood on the plan's
proof modality. W1.13 lands the declaration change *(landed — §4's
"circular complex GP" subsection is the body record, with the
composition-time Phase-4 refusal)*; the implementations are Phase 4's. Question 5 — ruled as it stands: the use case is one
`scale` per instrument or survey, which the per-dataset scalar already
expresses (each dataset carries its own `Likelihood`); a per-channel scale
*within* one container is not wanted (one spectrum rarely carries the data
to constrain more), and where scales are physically related across
instruments (orders of an échelle spectrum, say) the relationship is
expressed by tying or a hierarchical prior over the per-dataset scale
parameters, not by an array-valued scale with a grouping declaration.
Question 7 — resolved with `results_schema.md` §17 Q5's same-day ruling:
per-visibility frequency is real (visibilities are functions of (u, v, λ)
in general), so `extra_coords` gains `Axis` support eventually — Phase 4,
not the freeze; this contract needs nothing. Question 8 — superseded:
serialisation is to be consolidated *once* across all the contracts (an
inventory of every `to_spec`/`to_dict`/emission mechanism, a gap
analysis, and one coherent approach) rather than settled piecemeal;
routed to W1.13 alongside `results.md` §15 R7. *(Closed at the freeze:
the review is `docs/design/serialisation_review.md`, and it confirmed
the promotion — `Likelihood.to_spec()` exists, §8 above documents it,
and `describe_likelihood` composes it.)* **And ruled later the same
day**: X-1 — the prediction-aware `NoiseModel`
(`awkward_instrument.md` §6's detailed design) — is **accepted as
written** and landed at the freeze (§5's "Prediction-aware noise"
subsection is the body record): `sigma` and `noise_params` gain a
keyword-only `predicted=None` (the retained predicted values, passed at
every call site including `conditional` and `draw_observation`; an
outright signature change, pre-freeze, no shims), §14 gains the
noise-model-not-family principle, §15 retains the
no-model-parameters limitation, and §5's list names
`FractionalModelNoise` (implemented with the reference backend).

1. **Should a non-positive-definite covariance be `−inf` instead of an
   exception?** This contract raises, on the argument that a silent `−inf`
   hides a mis-specified kernel and that W1.7 owns failure signalling. The
   counter-argument is real: during sampling this *will* happen for extreme
   hyperparameter draws, and every engine driver will then need the same
   try/except. If Peter prefers, the alternative is a `strict=False` mode on
   `DenseGP` that returns `−inf` and records the reason.
2. **Is `Likelihood` the right home for the mask union?** It currently reads
   both containers' `weights()`. The alternative is that W1.7's `Dataset`
   resolves the effective mask once and hands the likelihood a single
   pre-masked pair, which would make `log_prob` cheaper and move one concept
   out of this contract. That depends on whether an instrument's output mask
   can change between evaluations (it can, if a transformation masks on a
   parameter-dependent condition) — W1.5 should say.
3. **Rice's parameterisation.** Does the model predict the true amplitude, or
   the underlying complex value from which the amplitude is derived? The first
   is simpler; the second is what an interferometric model actually produces.
   Deciding fixes the interface, and the family can then be implemented in an
   afternoon.
4. **von Mises's concentration.** The natural mapping from a per-sample σ is
   `κ = 1/σ²`, exact only in the small-σ limit. The alternative is a fitted
   `κ` as an ordinary parameter, ignoring the container's uncertainties
   entirely. Closure-phase practice (W1.11) should decide.
5. **Should `IndependentNoise`'s `scale` be per-dataset or per-channel?** It is
   currently one scalar for the whole container. A per-instrument scale on a
   joint fit is the common case and is expressible today by giving each dataset
   its own `Likelihood`; a per-*channel* scale within one container would need
   an array-valued parameter and a grouping declaration this contract does not
   have.
6. **Should a correlated complex Gaussian be declared analytic?** One real
   kernel applied independently to the real and imaginary parts of a visibility
   residual *would* marginalise in closed form, so `complex_gaussian` +
   `GaussianProcessNoise` could be `ANALYTIC` rather than the `LATENT` it
   currently declares. It is declared conservatively because "which complex GP"
   is a genuine modelling question — equal component covariances? a non-zero
   pseudo-covariance? — that belongs with Phase 4's visibility modality (W1.11),
   not with this contract asserting an answer by setting a flag. The same
   applies, differently, to Student-t: a *multivariate* Student-t with scale
   matrix `K + D` is closed-form, but it is a single global scale mixture rather
   than the per-point outlier robustness the family is chosen for, so it is a
   different model and not a shortcut to `ANALYTIC`. Note the practical
   consequence: because `CONSUMES_LATENT_GP` is `False` for both, neither can
   currently be *composed* with a GP at all, so the conservative declaration is
   a refusal rather than a silently different model — which is the right way
   round, but does mean a Phase-4 decision unblocks real functionality.
7. **Does `extra_coords` need units after all?** `results_schema.md` §15.4 and
   §17.5 flag this and name W1.6 as the possible requester. This contract does
   not need them — `limit_kind` is a code, not a quantity — so the answer from
   here is no. A per-visibility frequency for a frequency-dependent noise model
   would want them, but that is a Phase 4 question.
8. **Where does a `Likelihood` get serialised?** `ParameterSet.to_spec()`
   covers its parameters, and `KernelSpec.to_dict()` its kernel, but there is
   no `Likelihood.to_spec()` for the family name plus solver plus censoring.
   W1.8 owns provenance hashing and should decide whether it wants one.
