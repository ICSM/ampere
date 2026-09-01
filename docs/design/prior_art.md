# Prior-art memo (W1.1)

Status: **draft, awaiting Peter's review**. This memo studies five prior-art
systems named in `DEVELOPMENT_PLAN.md` §4 and §4.8, and extracts concrete
interface lessons for the Phase 1 contract specs (W1.3–W1.9). It does not
propose contract text itself — it names what to copy or avoid, and why, so
the contract authors can cite it directly (e.g. "see `prior_art.md` §2,
lesson G2").

Each system section below is tagged with the `DEVELOPMENT_PLAN.md` §4
subsection(s) it bears on, per the assignment in `WORK_ITEMS.md` W1.1 and
`DEVELOPMENT_PLAN.md` §4's own prior-art pointer. One correction to the
dispatch brief: gammapy's `Datasets`/joint-fitting material maps to
`DEVELOPMENT_PLAN.md` **§4.5** (Inference contracts — the `Dataset` /
`DatasetCollection` containers that work item W1.7 turns into
`docs/design/contracts/inference.md`), not §4.7 (which is the Astropy
interop contract and is unrelated to dataset tying). Lessons below are
filed under §4.5, with a note where they said.

A confidence/sourcing summary is at the end (§7), ahead of the final report.

---

## 1. bilby — likelihood/prior/sampler decoupling (→ §4.5, §4.1, §4.6)

**What it is.** `bilby` (Ashton et al. 2019, ApJS 241, 27) is the standard
Bayesian inference library for gravitational-wave parameter estimation. Its
core (`bilby.core`) is deliberately domain-agnostic: `Likelihood`,
`Prior`/`PriorDict`, `Sampler`, and `Result` are separated so any of ~10
samplers (dynesty, emcee, ptemcee, nestle, pymultinest, …) can run any
likelihood against any prior through one contract.

**Sources consulted** (fetched 2026-09-01): `bilby.core.likelihood`,
`bilby.core.sampler.base_sampler`, `bilby.core.prior.dict` source at
`github.com/bilby-dev/bilby` (`main` branch — the project's current
canonical home, successor to the original `git.ligo.org/lscsoft/bilby`);
the bilby API docs (`lscsoft.docs.ligo.org/bilby`) via search; Ashton et al.
2019 (arXiv:1811.02042) via search summaries only, not fetched directly.

### Lesson B1 — explicit `parameters` argument, not mutable shared state (§4.5). **Copy.**

The current `bilby.core.likelihood.Likelihood` base class is:

```python
class Likelihood:
    def __init__(self):
        self._meta_data = None
        self._marginalized_parameters = []

    def log_likelihood(self, parameters):
        return np.nan

    def noise_log_likelihood(self):
        return np.nan

    def log_likelihood_ratio(self, parameters):
        return self.log_likelihood(parameters=parameters) - self.noise_log_likelihood()
```

Every concrete subclass checked (`GaussianLikelihood`,
`PoissonLikelihood`, `StudentTLikelihood`, and — importantly — the
flagship `GravitationalWaveTransient`) implements
`log_likelihood(self, parameters)` taking an explicit `dict`, not a
`log_likelihood(self)` reading a `self.parameters` attribute mutated by
the sampler beforehand. `base_sampler.Sampler` builds the parameter dict
itself each call and passes it in.

This is worth flagging precisely because it is *not* how bilby is usually
remembered from its 2019 paper and older tutorials, which document the
`self.parameters = {...}` pattern (sampler sets `likelihood.parameters`,
then calls `likelihood.log_likelihood()` with no arguments). Bilby's own
codebase has evidently moved away from that convention toward explicit
argument-passing — i.e. the project's own history is evidence in favour
of the pure-function shape. This is exactly the shape ampere's §4.5
`log_prob(params) -> float` contract already commits to, and it is a
precondition for jax's `jit`/`vmap`/`grad` (mutable instance state defeats
tracing and makes batched/vectorised evaluation fragile). **Verdict: keep
ampere's `log_prob(params)` as a pure function of an explicit parameter
argument; do not offer a `self.parameters`-mutation convenience path even
for interactive/notebook ergonomics** — 3ML and gammapy (§2, §3 below) show
where that convenience leads.

### Lesson B2 — `PriorDict` has no cross-dict tying mechanism (§4.1). **Gap to avoid inheriting.**

`PriorDict(dict)` computes `rescale(keys, theta)` (the nested-sampling
prior transform) and `ln_prob()` (joint log-prior) over its own keys, and
distinguishes `fixed_keys` / `non_fixed_keys` (via `DeltaFunction` priors
for constants). But GW parameter estimation rarely needs to tie a
parameter across two independent likelihood objects — each event is one
`PriorDict`. Consequently bilby has no first-class "this parameter in
model A is this parameter in model B" primitive; `ConditionalPriorDict`
only expresses *within-dict* conditional dependence (`condition_func`,
`required_variables`), not cross-object identity.

Ampere's target scope (joint SED + spectroscopy fits, population models,
DatasetCollection) needs tying as a first-class citizen (§4.1, §4.5) —
bilby simply never had to solve this problem, so there is nothing to copy
here, but it is worth naming explicitly so nobody assumes "do what bilby
does" answers the tying question. See §2 below for a system that *does*
attempt it, and the pitfall it hit.

### Lesson B3 — `Result` is a sampler-parametrised grab-bag, not a fixed schema (§4.6). **Avoid; validates ampere's decision.**

`bilby.core.result.Result` stores posterior samples as a `pandas.DataFrame`
plus `log_evidence`, `log_evidence_err`, `information_gain`, and
serialises via `save_to_file` to json/hdf5/pickle. `Sampler.__init__`
accepts a `result_class` override so different samplers (or users) can
subclass `Result` for extra fields. This is flexible but exactly the
"every engine grows its own output shape" failure mode that
`DEVELOPMENT_PLAN.md` §4.6 already rejects in favour of a single ArviZ
`InferenceData` format emitted by every engine. Nothing to port here; the
memo simply records that bilby is evidence *for* ampere's existing
decision, not a design to imitate. One incidental strength worth noting:
`GaussianLikelihood.log_likelihood` dispatches on
`xp = array_module(self.x)` before doing array ops — i.e. bilby already
does lightweight array-API-style backend dispatch inside a likelihood.
That is a small, working precedent for the `array_api_compat` aspiration
noted in `DEVELOPMENT_PLAN.md` §5 Phase 2.

---

## 2. gammapy — the `Datasets` container and joint fitting (→ §4.5 / W1.7, §4.1)

**What it is.** `gammapy` is the standard analysis package for
ground-based gamma-ray astronomy (Cherenkov telescopes) and Fermi-LAT.
Its `Datasets` container is the closest existing analogue to ampere's
planned `DatasetCollection`: joint likelihood fitting across
heterogeneous instruments (different IRFs, different data types — 3D
cube, 1D spectrum, flux points) with parameters shared between them.

**Sources consulted** (fetched 2026-09-01): `docs.gammapy.org/1.3`
user guide (`datasets/index.html`) and the multi-instrument joint-fit
tutorial (`analysis-3d/analysis_mwl.html`); GitHub issue
`gammapy/gammapy#2859` ("Linking parameters of two models").

### Lesson G1 — one `Datasets` container, additive joint log-likelihood (§4.5). **Copy the shape.**

A gammapy `Dataset` bundles data + IRFs/response + background model +
the `Models` (parameters) assigned to it; `Datasets` holds a list of
`Dataset`s. `Fit.run(datasets=datasets)` computes total log-likelihood as
the **sum of each dataset's own fit statistic** — explicitly *not* equal
to first stacking the data and fitting once ("the total fit statistic of
datasets is the sum of the fit statistic of each dataset... this is not
equal to the stacked fit statistic"). This is precisely the joint-vs-stack
distinction `DEVELOPMENT_PLAN.md`'s W1.7 `DatasetCollection` needs to get
right: joint fitting over per-dataset noise models is the general case;
stacking is a fast-path optimisation that changes the objective and must
never be silently substituted. **Copy**: `DatasetCollection.log_prob` as a
sum over member `Dataset.log_prob`, with any stacking/reduction path
opt-in and documented as changing the objective, not merely accelerating
it.

### Lesson G2 — parameter tying via Python object identity is a real footgun (§4.1). **Avoid this implementation strategy.**

Gammapy ties parameters by assigning one model's `Parameter` object onto
another (documented pattern: `pwl2.index = pwl.index`, or assigning the
same `Models` object to two `Dataset`s wholesale, per lesson G1). This is
ergonomic — no special API, just Python attribute assignment — but GitHub
issue #2859 shows the cost: a joint fit over datasets with linked
parameters produced a covariance-matrix shape mismatch (`ValueError:
... cannot assign 49 input values to the 144 output values`) because the
covariance-expansion code counted each *dataset's reference* to a shared
parameter rather than each *unique* parameter, double-counting tied
values. The bug was fixed downstream (PR #2861), but it is a clear
instance of a general risk: **when tying is implemented as "the same
Python object appears in two places," every piece of code that walks the
parameter tree (covariance matrices, serialisation, pytree flattening for
jax) must independently get deduplication-by-identity right, and it is
easy for one of them not to.**

Ampere's §4.1 tying/sharing design should therefore be **name- or
group-based** (an explicit tying declaration — "these two named
parameters, across these two models/datasets, are the same free
parameter" — resolved once into a canonical parameter list at
construction time) rather than relying on shared object identity
propagating correctly through every consumer, especially given jax's
pytree flattening is exactly the kind of code that would need to
special-case aliased leaves if identity-tying were used. This is doubly
important because ampere's tying must additionally survive lowering to
`torch.distributions`/numpyro plates (§4.1, §W1.9) — another consumer that
would need identity-awareness bilby and gammapy never had to give theirs.

### Lesson G3 — `Dataset` bundles IRF *and* background *and* data *and* model (§4.2/§4.3). **Contrast, not copy wholesale.**

A gammapy `Dataset` is a fairly monolithic unit: response function,
background model, counts, and the assigned `Models` all live on one
object, and `Datasets` is just a list of these. This works well for
gamma-ray analysis, where IRF and background genuinely are
per-observation and rarely reused verbatim elsewhere. Ampere's split is
more granular by design — Instrument-as-transformation-chain (§4.3) is
independent of Likelihood/NoiseModel (§4.4), and both are independent of
the `Dataset` that binds them to observed data (§4.5) — because ampere
explicitly wants Instruments (an LSF, a resampling operator, a
calibration factor) to be reusable building blocks composed differently
per channel, not baked into a single per-observation bundle. Gammapy's
choice is not wrong for its domain, but it is evidence that "bundle
everything into one Dataset object" is a viable simpler alternative if
ampere's finer-grained composition ever proves to add friction without
adding reuse — worth remembering as a fallback, not adopting now.

---

## 3. 3ML — per-instrument plugin likelihoods (→ §4.3, §4.4)

**What it is.** 3ML (the "Multi-Mission Maximum Likelihood framework") is
built for joint multi-messenger/multi-instrument fits (e.g. a GRB seen by
Fermi-GBM, Fermi-LAT, and a ground telescope simultaneously) where each
instrument's data format and response are wildly different. Its answer is
a `Plugin` per instrument, each fully opaque to the framework.

**Sources consulted** (fetched 2026-09-01): `threeml.readthedocs.io`
"Building Custom Plugins" tutorial (`v2.3.0/notebooks/custom_plugins.html`)
and search-summarised results for `PluginPrototype`, `DataList`,
`JointLikelihood`/`BayesianAnalysis` from the 3ML docs and API reference.
The main "Model construction" tutorial page returned HTTP 404 at the URL
tried; not independently re-verified — see confidence note in §7.

### Lesson 3M1 — the plugin owns its response *and* its likelihood, fully opaque to the framework (§4.3/§4.4). **The key tension this memo was asked to surface.**

3ML's `PluginPrototype` contract is minimal: a subclass implements
`set_model(model)` (receive the shared astrophysical model), `get_log_like()`
(return a float — "no restrictions are placed on how this number is
calculated allowing for it to be the product of complex instrument
software, mathematical formulas, etc."), and `inner_fit()` (profile out the
plugin's own nuisance parameters). Everything about how the instrument
response is applied to the model to produce a likelihood — resampling,
convolution, background subtraction, an entirely custom covariance
structure — happens *inside* the plugin, invisible to 3ML itself.

This is the opposite factoring from `DEVELOPMENT_PLAN.md` §4.3/§4.4, which
deliberately splits "Instrument as a chain of reusable Transformations"
from "Likelihood compares predicted vs observed given a swappable
NoiseModel." 3ML's monolithic-plugin design **trivially satisfies
extensibility** (bring an arbitrarily weird instrument, write one class,
done — nothing outside the plugin needs to understand it) at the cost of
**zero reuse**: no two instruments' plugins share a resampling routine or
a noise-model implementation unless their authors manually factor that out
themselves, and 3ML's own plugin library shows a lot of near-duplicated
response-handling code across instruments as a result.

**Verdict: keep ampere's factored §4.3/§4.4 split** — it is the right call
for a system whose stated goal is a *library* of reusable Transformations
(LSF, resampling, synthetic photometry, calibration) and NoiseModel
strategies (DenseGP, QuasisepGP) shared across every instrument, which
3ML's model does not give. **But treat 3ML's plugin as the acceptance
test for §4.3's extensibility requirement**: the out-of-tree extension
example that `DEVELOPMENT_PLAN.md` §4.3 requires ("a user-defined
Transformation imported as if third-party") should include at least one
case, in W1.11's modality sketches, of an instrument whose response
genuinely *doesn't* decompose cleanly into (generic Transformation chain)
+ (generic NoiseModel) — e.g. one where the calibration uncertainty and
the noise covariance are physically coupled — to prove ampere's factoring
doesn't quietly force such users back into 3ML-style monkey-patching a
single opaque class. If the factored contracts can't express that case
cleanly, that is exactly the kind of interface gap W1.11 exists to catch
before the freeze.

### Lesson 3M2 — `get_log_like()` reads a shared, externally-mutated `Model` (§4.5). **Avoid** — same failure mode as bilby's old convention.

The example plugin sets `self._model = model` in `set_model` and later
reads parameter values off it inside `get_log_like()`, which takes **no
arguments**. The optimiser/sampler mutates parameter values on the shared
`astromodels` `Model` object between calls. This is the stateful,
mutation-based pattern that bilby's own core has since moved away from
(Lesson B1) — and for the same jax-incompatibility reason. It is
mentioned here specifically because it is a second independent precedent
(3ML plus old-bilby) for the same anti-pattern, which strengthens the
case for ampere holding the line on explicit-argument `log_prob(params)`
even though it is a less immediately ergonomic API for quick interactive
use than "just set an attribute and call a no-arg method."

### Lesson 3M3 — `get_log_like()` as the entire likelihood contract (§4.4). **Copy the minimalism, not the opacity.**

The actual interface surface 3ML asks a plugin author to implement is
tiny: essentially one method returning a float. That minimalism is worth
copying in spirit for ampere's §4.4 family-registry interface
(`log_prob(predicted, observed, noise_params)`) — a user adding a new
likelihood family should not need to implement more than the equivalent
of "here is my log-probability." The difference is where the *inputs* to
that method come from: 3ML leaves "predicted" entirely to the plugin's own
internals (opaque), while ampere's contract gets `predicted` from the
generic Instrument chain (§4.3) and `observed`/`noise_params` from the
generic `ModelResult`/`Dataset` machinery (§4.2/§4.5) — so the user-written
part is smaller and more testable in isolation than a 3ML plugin, which
is the whole point of factoring instrument from likelihood.

### Lesson 3M4 — `DataList` + additive joint likelihood (§4.5). **Copy; consistent with gammapy G1.**

3ML's `DataList` combines multiple plugins for `JointLikelihood`/
`BayesianAnalysis` in the same "sum of per-plugin log-likelihoods" shape
gammapy uses. Two independent domains converging on the same joint-fit
shape (additive log-likelihood over a heterogeneous collection, one shared
parameter space) is reassuring confirmation that `DatasetCollection`'s
planned shape (§4.5, W1.7) is the right one — nothing novel to design
here, just corroboration.

---

## 4. Starfish — misspecification-robust GPs for spectra (→ §4.4, touches §4.8)

**What it is.** Starfish (Czekala, Andrews, Mandel, Hogg & Green 2015,
ApJ 812, 128, arXiv:1412.5177) is the direct scientific ancestor of
ampere's flexible-likelihood idea: it models correlated spectral residuals
from imperfect synthetic-spectrum models with a Gaussian-process
covariance term added to the noise, so the fit is robust to features the
underlying model gets wrong, rather than corrupting the physical
parameter estimates. This is the same idea ampere is generalising and
scaling.

**Sources consulted** (fetched 2026-09-01): the paper text via
`arxiv.org/html/1412.5177` (the PDF itself failed to parse as text — see
§7); the Starfish docs landing page (`starfish.readthedocs.io`, low
technical detail); the GitHub repository `github.com/iancze/Starfish`
(metadata only — no maintenance-recency or successor-project detail could
be confirmed from what was fetched).

### Lesson S1 — Matérn-3/2 as the global kernel, empirically validated (§4.4). **Copy; corroborates the plan's decision.**

Starfish's "global" covariance kernel is exactly Matérn-ν=3/2 in velocity
separation:

```
K_G(r_ij) = w_ij · a_G · (1 + √3 r_ij/ℓ) · exp(−√3 r_ij/ℓ)
```

This is independent, published (2015), empirically-validated evidence for
`DEVELOPMENT_PLAN.md`'s own decision (§2, §4.4) to make Matérn-3/2 the
default flexible-likelihood kernel: the closest prior art to ampere's core
scientific idea arrived at the same kernel family for the same reason
(better representation of structured residuals than squared-exponential).
Nothing to change here — this is reassurance, not a correction.

### Lesson S2 — sparsity via **windowed truncation**, not exact state-space recursion (§4.4/§4.6). **Note the difference; do not conflate the two O(N) strategies.**

Starfish avoids the dense O(N³) covariance solve, but not via a
quasiseparable/state-space recursion (celerite-class, what
`DEVELOPMENT_PLAN.md` §4.4 specifies for `QuasisepGP`). Instead it applies
a Hann window to *taper the kernel itself to exactly zero* beyond a cutoff
radius (`r₀ = 4ℓ` for the global kernel, `r₀ = 4σ_k` for each local
kernel), producing a genuinely **sparse banded** covariance matrix that a
sparse Cholesky factorisation can invert in roughly linear time for
well-separated bands. This is a real, different O(N) strategy from exact
quasiseparable solves: it is an *approximation* (the true covariance has
non-zero tails the window discards, traded for sparsity), whereas
celerite-class quasiseparable kernels give an *exact* O(N) solve for the
un-truncated kernel. `DEVELOPMENT_PLAN.md`'s `NoiseModel` strategy slot
(§4.4) already anticipates multiple solver strategies (`DenseGP`,
`QuasisepGP`, future SVGP/SKI/Vecchia for 2D+) — this memo's contribution
is simply to flag that **windowed-sparse is a third, qualitatively
different strategy family** (approximate-but-simple, no state-space
machinery required, easy to reason about locally) that is worth a named
slot or at least an explicit "why we didn't use this" line in the
likelihoods contract (W1.6), since it is what the field's own prior art
actually shipped and a reviewer familiar with Starfish will ask why
ampere didn't.

### Lesson S3 — variable-dimension local kernels do not fit a fixed-shape parameter contract (§4.1/§4.4). **Avoid porting directly; explains a design choice worth stating explicitly.**

Beyond the single global kernel, Starfish adds one **local** Gaussian
"bump" kernel per identified problem spectral line, each contributing 3
hyperparameters (amplitude, centre, width): "if there are N_loc local
covariance patches ... there are 4·N_loc + 2 elements in the set of
covariance hyperparameters." Critically, **N_loc is not fixed** — it is
itself chosen by the analyst/algorithm as part of fitting, i.e. this is a
trans-dimensional (variable-dimensionality) parameter space. Starfish
handles it with a blocked Gibbs/Metropolis-Hastings sampler, not a fixed
gradient-friendly parameter vector.

Ampere's §4.1 parameter contract is explicitly built around a
**statically-declared, plate-aware but fixed-shape** parameter space,
because that is what lowers cleanly to `torch.distributions`/numpyro
plates (§4.1, W1.9) and is required for NUTS/VI. A Starfish-style
trans-dimensional local-kernel model cannot be expressed in that contract
without a reversible-jump/model-selection layer nobody has scoped. This is
not a flaw to fix now — it is a **scope boundary worth stating
explicitly** in the likelihoods contract (W1.6): ampere's single
stationary Matérn kernel over the whole domain, plus the GP-localisation
diagnostic already planned in §4.8 ("the posteriors on GP amplitude and
length-scale, and the conditioned GP mean, already localise where the
model is deficient"), is the *substitute* for Starfish's explicit
per-line local kernels — a global model that tells you *where* it's
locally wrong, rather than a model with explicit local components. Naming
this trade-off in the spec pre-empts confusion later about why ampere
"only" has one kernel.

### Lesson S4 — global/local amplitude degeneracy is a real diagnostics concern (§4.8). **Feed forward to W1.12.**

The paper reports that its global kernel amplitude and its local kernels
trade off against each other on real data (Gl 51): "there is little
difference between the posteriors in the third and fourth tests," i.e.
either a larger global-kernel amplitude or extra local kernels explain
similar residual structure, and sensitivity is limited once amplitudes are
near-optimal. Since ampere's diagnostics spec (§4.8, W1.12) plans to
surface GP amplitude/length-scale posteriors as localisation output, it
should document this degeneracy risk explicitly (a large fitted GP
amplitude with a short length-scale can mean either "genuinely
misspecified locally" or "the kernel's smooth global component is
under-amplitude and compensating locally") so users reading the
localisation plots don't over-interpret them as pinpointing a single
physical cause.

### Lesson S5 (minor) — a concrete cost baseline (§4.4, benchmark expectations). **Informational.**

The paper reports "a fit of an R≈40,000 spectrum with >30 echelle orders
takes ~2 hours, parallelised on a cluster" using the windowed-sparse
approach on 2015-era hardware. This is a useful order-of-magnitude
baseline for what "acceptable" meant for the field's own prior art at
comparable N, worth keeping in mind (not as a target — ampere's whole
point is to beat this by orders of magnitude via exact O(N) quasiseparable
solves per M2's flagship benchmark) rather than any specific number to
design toward.

---

## 5. RHMF / Robusta-HMF — pre-fit misspecification screening (→ §4.8, touches §4.2)

**What it is.** Hilder, Hogg, Casey & Rix, "Robust Heteroskedastic Matrix
Factorization: A Generalization of PCA that Flags Outliers and Handles
Missing Data" (submitted to ApJ 2026-07-09; arXiv:2607.08081). RHMF is an
iteratively-reweighted low-rank factorisation with an implicit Student-t
likelihood, per-feature (heteroskedastic) uncertainties, native
missing-data support, and automatic per-feature/per-object anomaly scores.
Applied in the paper to Gaia DR3 RVS spectra to flag anomalous stars. This
is `DEVELOPMENT_PLAN.md` §4.8's named candidate for pre-fit, data-driven
misspecification screening across a *collection* of spectra (distinct from
Starfish/ampere's per-object, post-fit GP approach in §4).

**Sources consulted** (fetched 2026-09-01, in order): the arXiv abstract
page (`arxiv.org/abs/2607.08081`); the full paper text via
`arxiv.org/html/2607.08081v1`, which is where the software-availability
statement, algorithm detail, and licence information below came from; the
associated code repository, `github.com/TomHilder/robusta-hmf`, fetched
directly. A related blog post by D. W. Hogg
(`hoggresearch.blogspot.com/2025/07/robust-matrix-factorization.html`)
was also checked and links to an *earlier*, different repository
(`github.com/davidwhogg/SpectralAnomalies`, described there as where
"development is happening ... for now") — that appears to be superseded
by `TomHilder/robusta-hmf`, which is what the published paper's own
software-availability text points to and what carries the citation. Treat
`TomHilder/robusta-hmf` as authoritative.

### Adoptability assessment (front-loaded for W1.12)

| Question | Finding |
|---|---|
| **Package / repo** | `robusta-hmf` on PyPI (`pip install robusta-hmf` / `uv add robusta-hmf`); source at `github.com/TomHilder/robusta-hmf`. |
| **Code licence** | **MIT** on the code repository. (The arXiv paper's *text* is separately CC BY 4.0 — that licence governs the paper, not the code; do not conflate the two when writing the diagnostics spec's dependency note.) MIT is fully compatible with ampere's own licensing and with adding it as an optional extra. |
| **Language / dependency** | JAX. Adding `robusta-hmf` as an optional extra therefore pulls in a JAX dependency for users of the diagnostics module specifically, independent of which backend (`reference`/`torch`/`jax`) they use for fitting. Needs a lazy-import guard consistent with `DEVELOPMENT_PLAN.md`'s "never make `import ampere` require torch or jax" policy (§7) — the diagnostics module choosing to depend on JAX for one optional feature is a new, narrower instance of that same policy question, not a violation of it, provided the import stays lazy and scoped to the diagnostics extra. |
| **Maturity** | Small but real: 272 commits, 8 GitHub stars, 0 forks (as fetched); has a `tests/` directory and a GitHub Actions CI workflow; ships `examples_paper/` reproducing the paper's toy and Gaia-RVS pipelines. This is research-grade software from one active group (Hogg's), not a widely-adopted community package — expect API changes without a deprecation cycle, and budget review time rather than treating it as a stable dependency. |
| **API shape** | scikit-learn-esque: `Robusta(rank=K, robust=True, robust_scale=Q)`, `state, loss_history = model.fit(Y, W, max_iter=...)`, `model.synthesize()`. `Y` is the data matrix, `W` the per-element weight/uncertainty matrix (zero weight ≡ missing, per the paper's own equivalence of "setting weights to zero" and "setting uncertainties to infinity" for missing data). |
| **Verdict** | **Adoptable as an optional dependency**, not to be vendored or reimplemented from scratch. MIT licence and a real (if small) test suite clear the bar; the JAX dependency and small-team maturity are the costs to plan around. W1.12 should scope RHMF-style screening as "depend on `robusta-hmf` behind a lazy import in an optional diagnostics extra," matching how the plan already treats `sbi`/`torch`/`jax` as extras (§2 packaging decision), rather than porting the algorithm into ampere's own code. |

### Lesson R1 — masks/missing-data via zero-weight or infinite-uncertainty (§4.2). **Copy; already the direction ampere is going.**

RHMF's missing-data handling — "setting weights `w_ij` to zero or
(equivalently) setting the corresponding measurement uncertainties `σ_ij`
to infinity" — is precisely the semantics `DEVELOPMENT_PLAN.md` §4.2
already commits to for first-class masks ("masking... is distinct from
censoring... both are part of the likelihood contract"). This is
corroboration, not a new decision: when W1.4/W1.6 write the concrete mask
propagation rules, RHMF's weight-matrix convention is a working existing
implementation of the same idea and can be cited as precedent, and its
adoption as a diagnostics dependency will compose cleanly with ampere's
mask representation with no impedance mismatch (both express "missing" as
zero-weight/infinite-uncertainty, not as sentinel NaNs threaded through
arithmetic).

### Lesson R2 — per-element robust weights double as an anomaly/localisation signal (§4.8). **Copy the idea, not necessarily the exact statistic.**

RHMF's robust weights `w_ij^robust ∈ (0, 1]` from the IRLS reweighting
step are, by construction, a per-pixel "how much did this iteratively-
reweighted fit trust this data point" score; the paper aggregates them
per-object (e.g. a low quantile across features) to flag anomalous
objects. This is conceptually the same move `DEVELOPMENT_PLAN.md` §4.8
already plans for the *post-fit* GP path ("the posteriors on GP amplitude
and length-scale... already localise where the model is deficient") —
RHMF gives an analogous **pre-fit**, collection-level signal, cheaply,
before any physical model has even been evaluated. W1.12 should treat
RHMF's per-feature weight map and ampere's post-fit GP localisation as two
views of the same underlying question (where does a smooth/low-rank model
fail to explain this data) at two different stages of the pipeline, and
design the diagnostics module's output schema so both can be plotted with
the same conventions in `ampere.results` (shared colour scale / axis
conventions for "how anomalous is this point," even though the two
statistics are computed completely differently).

### Lesson R3 — two free hyperparameters need tuning, no automatic default (§4.8, UX note). **Anticipate this in the spec.**

Both the blog post and the paper are explicit that "the investigator has
to tune the rank of the low-rank part and also the soft outlier cutoff"
(rank `K` and `robust_scale`/`Q`), with cross-validation recommended for
reliable settings. This is a real cost for a *pre-fit screening* tool
meant to be cheap and automatic: if a user has to cross-validate two
hyperparameters before RHMF's flags are trustworthy, it stops being a
"quick look before you commit to the expensive fit" step and becomes its
own small research project. W1.12's diagnostics spec should either (a)
ship sane, documented defaults validated against a representative ampere
use case (SEDs/low-res spectra) so most users never need to tune
`rank`/`robust_scale` by hand, or (b) explicitly scope RHMF screening as
"opt-in, expert-tunable" rather than an always-on default step, so it
doesn't become friction in the common path.

---

## 6. Synthesis — cross-system tensions for W1.3 onward to resolve

The five systems above do not agree with each other, and the
disagreements are exactly the shape of decision the Phase 1 contract specs
need to make explicitly rather than by accident.

**Tension 1 — identity-based vs name-based parameter tying (§4.1).**
Gammapy ties parameters by Python object identity (`pwl2.index = pwl.index`,
or sharing a whole `Models` object across `Dataset`s) and hit a real bug
(§2, Lesson G2) because downstream code (covariance expansion) didn't
correctly deduplicate by identity. Bilby never needed cross-object tying at
all (§1, Lesson B2) — no lesson to borrow, just a gap. Ampere needs the
gammapy *use case* (assign one physical parameter to two models/datasets)
without the gammapy *implementation strategy* (bare object-identity
sharing) that broke downstream consumers. **Resolution for W1.3**: tying
must be an explicit, named declaration resolved once into a canonical
free-parameter list at `FittingProblem`/`DatasetCollection` construction
time (§4.1, §4.5) — never "the same Python object happens to appear twice
in the tree," specifically because ampere has *more* tree-walking
consumers than gammapy did (jax pytree flattening, torch parameter
registration, `torch.distributions`/numpyro plate lowering per W1.9, plus
InferenceData provenance) that would each need independent
identity-deduplication logic if tying stayed identity-based.

**Tension 2 — pure-function `log_prob(params)` vs mutable shared model state (§4.5).**
Bilby's current core (§1, Lesson B1) and ampere's own §4.5 contract both
commit to explicit-argument, side-effect-free `log_prob(params)`. Gammapy
(shared `Models` object mutated by the optimiser) and 3ML (`get_log_like()`
reading a `self._model` set once, mutated externally between calls) both
use the opposite, stateful convention (§2 Lesson G1's plumbing; §3 Lesson
3M2) — and it is the more common convention in the astronomy-fitting
ecosystem generally, not the minority one. This is worth stating plainly
because a contract author or reviewer coming from gammapy/3ML/older-bilby
experience will find ampere's insistence on pure functions unfamiliar and
may be tempted to add a mutable convenience layer "for interactive use."
**Resolution for W1.3/W1.7**: hold the line on explicit arguments in the
core contract (required for jax tracing and for the batched/vectorised
evaluation `DEVELOPMENT_PLAN.md` §3 promises on native backends) — if an
interactive convenience wrapper is wanted later, it must be a thin
adapter that constructs the parameter dict and calls the pure function,
never a code path the core contract itself depends on.

**Tension 3 — monolithic per-instrument likelihood vs factored Transformation/NoiseModel (§4.3/§4.4).**
3ML deliberately keeps instrument response and likelihood computation
fully inside one opaque `Plugin` (§3, Lesson 3M1); ampere deliberately
factors Instrument (chain of reusable Transformations) apart from
Likelihood/NoiseModel (swappable strategy), and gammapy sits in between
(response+background+data bundled per-`Dataset`, but the noise/fit
statistic is a shared, generic mechanism, not per-instrument code). 3ML's
choice buys trivial extensibility for genuinely weird instruments at the
cost of zero code reuse across instruments; ampere's choice buys reuse at
the risk that some real instrument's physics won't decompose cleanly into
the two pieces. **This is not resolvable from prior art alone** — it
depends on whether ampere's actual target instruments (SED photometry,
low-res + high-res spectral windows, eventually visibilities) do decompose
cleanly. **Resolution path**: W1.11's modality sketches must include at
least one deliberately-awkward instrument (per Lesson 3M1's suggestion) as
a stress test of the §4.3/§4.4 split before the spec freeze (W1.13), not
after.

**Tension 4 — fixed-shape parameter contract vs Starfish's variable-dimension local kernels (§4.1/§4.4).**
Starfish's local-kernel component (§4, Lesson S3) is trans-dimensional —
not expressible in a statically-shaped parameter vector at all. Ampere's
§4.1 contract is fixed-shape by design (required for the
torch.distributions/numpyro lowering in W1.9). There is no tension to
"resolve" here in the sense of making them compatible — the resolution is
to **document the boundary explicitly**: ampere's single global Matérn
kernel plus GP-localisation diagnostics (§4.8) is a deliberate, different
answer to the same scientific problem Starfish's local kernels solve, not
an oversight. W1.6 (likelihoods contract) should say this outright so a
future reader doesn't propose "let's add Starfish-style local kernels" as
if it were a straightforward extension — it would require a
model-selection/reversible-jump layer nothing in the current plan scopes.

**Tension 5 — where does pre-fit screening's output live relative to post-fit diagnostics (§4.8)?**
RHMF (§5) is a *collection-level, pre-fit* screening tool with its own
output shape (per-feature/per-object weight maps); Starfish/ampere's GP
localisation (§4, §4.8) is a *per-object, post-fit* diagnostic derived from
the fitted GP's own posterior. Both answer "where is a smooth model
inadequate," at different pipeline stages, with unrelated underlying
statistics. `DEVELOPMENT_PLAN.md` §4.8 already lists both but doesn't say
whether they share a visual/output convention. **Resolution for W1.12**:
decide explicitly whether pre-fit RHMF flags and post-fit GP localisation
share a common "anomaly/deficiency score" presentation in `ampere.results`
(Lesson R2) or are kept visually distinct because conflating two
differently-computed statistics under one colour scale would overstate
their comparability — either is defensible, but W1.12 should say which and
why.

---

## 7. Confidence and sourcing, by section

| Section | Confidence | Basis |
|---|---|---|
| §1 bilby | **High.** | Read actual current source (`likelihood.py`, `base_sampler.py`, `prior/dict.py`) from `bilby-dev/bilby` `main`, not just docs/search summaries. The one soft spot: I did not independently confirm the exact history/date of bilby's migration away from the `self.parameters`-mutation convention documented in the original 2019 paper — I infer the shift from what the current source shows, not from a changelog entry I read myself. |
| §2 gammapy | **High** for the `Datasets`/joint-fit shape and the G2 bug report (read the actual GitHub issue). **Medium** for whether PR #2861's fix changed the *public* tying API (e.g. added an explicit `link()`/`link_label` mechanism) — the fetched issue summary confirmed the bug and that it was closed via a PR, but I did not read the PR diff itself, so I cannot state precisely what gammapy's tying API looks like today. |
| §3 3ML | **Medium.** | The primary tutorial page I first tried (`Model_construction.html`) 404'd; I fell back to the "Building Custom Plugins" tutorial and search-engine summaries of the API reference for `DataList`/`JointLikelihood`/`BayesianAnalysis`, which is a shallower source than direct code reads. The core claims (plugin opacity, `get_log_like()` signature, `set_model` role) come from the tutorial's own example code, which I consider reliable, but I have not read `plugin_prototype.py` source directly the way I did for bilby's likelihood module. If W1.6/W1.11 need finer 3ML API detail, a direct source read is worth doing before relying further on this section. |
| §4 Starfish | **High** for the kernel mathematics, hyperparameter count, and the cost/degeneracy findings — these came from the paper's own HTML full text (`arxiv.org/html/1412.5177`), not an abstract or summary. **Low** for the repository's *current* maintenance status or whether a later Starfish version adopted celerite/celerite2 — the PDF fetch failed to parse as text and the GitHub repo fetch returned only metadata (no commit-recency or changelog detail); I could not confirm or rule out a more recent O(N) rewrite. Treat the "windowed-sparse, not quasiseparable" characterisation as accurate for the 2015 paper's method, which is what is cited in `DEVELOPMENT_PLAN.md` §4, but do not assume it is still true of the latest Starfish release without checking. |
| §5 RHMF/Robusta-HMF | **High.** | This is the best-sourced section: fetched the paper's own full HTML text (software-availability statement, algorithm, licence) and fetched the actual code repository directly (confirmed MIT licence, install instructions, API shape, test/CI presence, star/commit counts). The one gap: I did not clone/run the package or read its source beyond the README, so claims about exact internal behaviour (e.g. numerical edge cases) are paper-level, not code-verified. |
| §6 Synthesis | **High** for identifying that the tensions exist (each rests on the sourced material above). **N/A for "correctness" of the proposed resolutions** — these are recommendations for W1.3+ to weigh, not settled decisions; that weighing is explicitly out of this memo's scope. |

**Overall**: the two most load-bearing, least-covered-by-training-knowledge
sections — RHMF/Robusta-HMF and bilby's current source — are the two most
confidently sourced, because both were checked against live fetches rather
than recalled. The weakest section is 3ML (tutorial + search summaries
only, one dead link), and the Starfish section's *current* maintenance
status (as opposed to its 2015 method) is unconfirmed. Both are flagged
above rather than smoothed over.

---

## Sources consulted

- bilby: `github.com/bilby-dev/bilby` (`main` branch) — `bilby/core/likelihood.py`,
  `bilby/core/sampler/base_sampler.py`, `bilby/core/prior/dict.py`; Ashton
  et al. 2019, ApJS 241, 27 (arXiv:1811.02042, referenced via search, not
  fetched directly).
- gammapy: `docs.gammapy.org/1.3/user-guide/datasets/index.html`;
  `docs.gammapy.org/1.3/tutorials/analysis-3d/analysis_mwl.html`;
  `github.com/gammapy/gammapy` issue #2859.
- 3ML: `threeml.readthedocs.io/en/v2.3.0/notebooks/custom_plugins.html`;
  `threeml.readthedocs.io` API reference pages via search summary.
- Starfish: Czekala, Andrews, Mandel, Hogg & Green 2015, ApJ 812, 128
  (arXiv:1412.5177, full text via `arxiv.org/html/1412.5177`);
  `starfish.readthedocs.io`; `github.com/iancze/Starfish` (metadata only).
- RHMF/Robusta-HMF: Hilder, Hogg, Casey & Rix 2026, submitted ApJ
  (arXiv:2607.08081, full text via `arxiv.org/html/2607.08081v1`);
  `github.com/TomHilder/robusta-hmf`;
  `hoggresearch.blogspot.com/2025/07/robust-matrix-factorization.html`.

All fetches performed 2026-09-01.
