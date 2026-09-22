# Inference approaches beyond the budgeted six — what fits, what needs adaptation

*Drafted 2026-09-10 by Fable at Peter's request, while Phase 3 closes. A
brainstorm and a fit assessment, not a plan: nothing here is scheduled, and
the numbered adaptations at the end are the only things a future work item
would need to land first. Long-horizon — Phase 5 or later, depending on how
the pieces interact with emulation (design horizon (c)) and population
reweighting (horizon (b)).*

## 0. The question, and the short answer

Ampere budgets three black-box samplers (`emcee`, `zeus`, `dynesty`), two
gradient-based engines (NUTS, VI) and the SBI layer. The contract's claim is
that adding an engine in any of these categories is cheap. This memo checks
that claim against the code as it stands, category by category, and asks
what else is worth adding.

The short answer: **the claim holds for samplers, and holds best for
methods that produce draws.** The three surfaces an engine can consume are
already there and are enough for every sampler in this memo:

| Slot | Surface | What it gives an engine | Engines on it today |
|---|---|---|---|
| A | `inference.md` §10 — `log_prob`, `log_likelihood`, `log_prior`, `prior_transform`, `evaluate`, `initial_positions`, the seeded streams | scalar numpy density, any backend, no gradients | emcee, zeus, dynesty |
| B | §10a — `realise(problem).log_prob_unconstrained(θ)` in the backend's array type | a traceable density: gradients, Hessians and batching by autodiff | NUTS, VI |
| C | §13 — `simulate`, `simulate_many`, the encoding | draws of data, no density | SBIEngine |

What is **not** there, and where every non-sampling method runs into the
same wall, is the *output* side. `Engine.finish` assembles equal-weight
draws shaped `(chain, draw, free)` into a run; `ampere.results` stores draws.
Three kinds of method produce something else:

1. **Weighted particles with an evidence** (nested sampling, SMC, PMC,
   importance-corrected Laplace). Precedent exists: `dynesty` resamples to
   equal weight and records `ampere_dynesty_logz` in the attrs. It is a
   per-engine convention, not a contract rule.
2. **An approximating distribution** (Laplace, Pathfinder, richer VI
   families). Precedent exists: `VIEngine` draws from the fitted guide and
   emits one "chain" of i.i.d. draws with `ampere_vi_guide` naming the
   family. Again per-engine.
3. **A point estimate with curvature** (MAP, MINUIT-style profile errors,
   Bayesian optimisation). No precedent and no home in the results schema.

So the recommendation is: add samplers freely; before adding *approximate*
or *optimising* engines, take three small contract decisions (§5) so that
each of them is one engine rather than one engine plus one schema
extension.

## 1. Black-box samplers on slot A

Every candidate here consumes exactly what `dynesty` consumes. Each is
listed with what it wants, what it returns, and the dependency it drags in.

### 1.1 Nested-sampling family — the cheapest additions in this memo

- **`nautilus`** (Lange 2023; importance nested sampling with a neural
  network boundary). Wants `prior_transform` (or a `Prior` object) and a
  likelihood of θ; optionally a *vectorised* likelihood. Returns weighted
  points and `log Z` with an error. Dependencies: numpy, scipy,
  scikit-learn (the network). Reported to need 5–10× fewer likelihood calls
  than `dynesty` on typical astronomical posteriors. **Fits slot A
  verbatim**; the driver is `_dynesty.py`'s shape with the resample step
  and the evidence attr. Strong first candidate.
- **`ultranest`** (Buchner; MLFriends). Same inputs; returns weighted
  points, `log Z`, and a rich diagnostic bundle; pure Python + numpy; MPI
  parallelism built in; vectorised likelihood optional. Its
  `ReactiveNestedSampler` wants the parameter *names*, which
  `free_labels()` already gives. **Fits slot A verbatim.** More robust than
  `dynesty` on pathological posteriors, slower on easy ones.
- **`jaxns`** (jax-native nested sampling). Wants a jax log-likelihood and
  a prior model written in its own `tfp`-based DSL — the prior side does
  **not** map onto `prior_transform` without a wrapper, and the wrapper
  would have to re-express every §4.1 prior. Slot B on jax only. Defer:
  cost is in the prior bridge, not the density.
- `pymultinest`, `nessai`, `dynesty`'s own dynamic mode: nothing new to add
  (dynesty's dynamic mode is already reachable through `run_nested`
  kwargs); MultiNest needs a compiled library; `nessai` (flow-based nested
  sampling) needs torch and largely overlaps with pocoMC below.

**Adaptation needed**: one engine-neutral evidence attribute (§5.1) and
the weighted-draw rule written down once (§5.2). With those, the three
nested samplers share one `_nested.py` base beneath `_dynesty.py`.

### 1.2 Sequential / population Monte Carlo

- **`pocoMC`** (Karamanis et al.; preconditioned Monte Carlo — SMC with a
  normalising-flow preconditioner). Wants `log_likelihood(θ)`,
  `log_prior(θ)` and prior bounds or a prior object with `logpdf`/`rvs`;
  `prior_transform` of a uniform block gives `rvs` for free. Returns
  weighted particles and `log Z`. Handles multimodal and strongly
  correlated posteriors with far fewer evaluations than an ensemble
  sampler, and its evaluations are **batched** (the whole particle
  population at once), which is exactly where §5.5's batched `log_prob`
  would pay. Dependency: **torch** (via `zuko` for the flows) — so an
  extra, not base. Strong second candidate.
- **Generic SMC** (tempered, with MCMC mutations): `blackjax.smc` on jax
  (slot B, gradient kernels optional) or hand-rolled on slot A with
  emcee-style moves. pocoMC is the better-engineered version of the same
  idea; add SMC through pocoMC unless a jax-only user needs it.
- **Parallel tempering**: `ptemcee` is unmaintained; **`eryn`** (Katz et
  al.) is the maintained PT-ensemble sampler in the emcee idiom. It brings
  reversible-jump moves, which are out of scope (variable dimension does
  not fit the fixed-shape parameter contract — the same reason
  `prior_art.md` gave for Starfish's local kernels), but its fixed-dimension
  PT path fits slot A and yields an evidence by thermodynamic integration.
  Medium value: pocoMC and the nested samplers cover multimodality with
  less machinery.
- **Ensemble Kalman inversion** (EKI / EKS): gradient-free, embarrassingly
  batch-parallel, approximate (Gaussian in the limit). Interesting only
  for very expensive simulators where SBI is the alternative; mention, do
  not schedule.

### 1.3 Marginal additions

Adaptive Metropolis, differential-evolution MC (`pydream`), affine-invariant
variants: emcee and zeus already occupy this ground. Not worth a driver
each.

## 2. Gradient-based methods on slot B

The realisation exposes the density as a pure function of the unconstrained
vector in the backend's array type. Gradients, Hessians and batching come
from `jax.grad`/`jax.hessian`/`vmap` and `torch.func` — the engines take
them, the contract does not have to offer them. **The two backends are not
symmetric in what the ecosystem offers**, and that asymmetry is the main
finding of this section.

### 2.1 On jax: `blackjax` is one dependency for six methods

`blackjax` wants `logdensity_fn(θ)` and nothing else — the realisation is
that function. It supplies:

- **NUTS/HMC** (an alternative to numpyro's; not needed);
- **MCLMC** (microcanonical Langevin Monte Carlo; Robnik & Seljak) — the
  strongest new sampler in this memo for high-dimensional smooth
  posteriors (latent GPs, hierarchical plates), reported at several times
  NUTS's efficiency per gradient;
- **Pathfinder** (Zhang et al. 2022) — L-BFGS along the optimisation path,
  a mixture-of-Gaussians approximation with importance resampling. Two
  uses: a cheap approximate posterior in VI's emission shape, and the
  best available **initialiser** for NUTS/MCLMC chains (its draws feed
  `initial_positions`);
- **mean-field VI** and **SMC** with HMC mutation kernels;
- adaptive MALA and others.

Fit: one `BlackjaxEngine(problem, method=...)` over the jax realisation, or
routes inside the existing NUTS/VI engines. The jax `sbi` and `dev`
environments already carry jax; blackjax is small. **Cheapest slot-B
addition, jax only.**

### 2.2 On torch: less to borrow

Pyro offers HMC/NUTS and the autoguides; there is no maintained torch MCLMC
or Pathfinder (`hamiltorch` is dormant). Pyro's autoguides, however, are
where three items in Peter's list already live, almost for free:

- **Laplace approximation**: `pyro.infer.autoguide.AutoLaplaceApproximation`
  and numpyro's equivalent — MAP by optimisation, then the Hessian at the
  mode by autodiff, then a multivariate normal *in unconstrained space*
  mapped back through `constrain`. `VIEngine(guide="laplace")` is a new
  entry in `GUIDE_FAMILIES` and nothing else; the emission shape (one chain
  of i.i.d. draws, the family recorded) is already right.
- **MAP**: `AutoDelta` is the MAP point through the same optimiser. It
  *runs* through `VIEngine` today's machinery, but emitting a point as a
  run of identical draws is wrong — this is the §5.3 gap.
- **Normalising-flow guides**: `AutoNormalizingFlow` (pyro),
  `AutoBNAFNormal`/`AutoIAFNormal` (numpyro). `guide="flow"` closes VI's
  documented weakness ("neither captures a non-Gaussian posterior") at the
  cost of a fit that needs more steps and more judgement about
  architecture.

So on **both** backends, "other VI schemes" and "Laplace" are mostly guide
families, and the engine is already there.

### 2.3 Gradient-free Laplace: `snowline`

Buchner's `snowline` is MAP + Laplace + importance-sampling refinement
against a black-box likelihood and `prior_transform` (numerical Hessians via
`iminuit`). It fits **slot A**, so it is the Laplace route for reference-
backend problems and wrapped simulators. Output is importance-weighted
draws plus an evidence — §5.1/§5.2 again. Small dependency (`iminuit`).
Worth adding beside the autodiff Laplace, because the two answer the same
question on the two halves of the capability ladder.

### 2.4 Other gradient MCMC

- **`nutpie`** (Rust NUTS with a jax bridge and normalising-flow mass-matrix
  adaptation). Takes `logp` and gradient callables; fits slot B on jax.
  Fast, but overlaps numpyro's NUTS and blackjax's MCLMC; low priority.
- Stochastic-gradient MCMC (SGLD, SGHMC): needs minibatched data; ampere's
  likelihoods are whole-dataset GP marginals. No.
- Riemannian/manifold HMC: research territory; no.

## 3. Optimisation

This is the category where the *contract*, not the ecosystem, decides.

### 3.1 What fits today

- **Deterministic optimisers against slot A**: `scipy.optimize` (Nelder–Mead,
  Powell, L-BFGS-B with numerical gradients) on `log_prob_unconstrained`;
  **`iminuit`** (MINUIT2) — HESSE covariance and MINOS profile intervals,
  the astronomer's habitual tool, with `free_labels()` as the parameter
  names. Every one of these runs *now*, in a notebook, against the public
  surface. What is missing is a place to put the answer and a bridge from
  the answer to the samplers.
- **Against slot B**: `optax`/`optimistix` on jax (`jaxopt` is being
  retired in favour of optimistix), `torch.optim.LBFGS`, with exact
  gradients and autodiff Hessians. Same story: runs now, nowhere to store.

### 3.2 What an optimisation result needs (§5.3)

A run's `posterior` group is the wrong container for a point. An
`Optimum` result — the mode in constrained *and* unconstrained coordinates,
the value of the density there, the inverse Hessian as a covariance (or a
refusal when the Hessian is not positive definite, by name), profile
intervals when a profiler ran, the optimiser's convergence status and
evaluation count, and the full provenance attrs — either as its own small
`DataTree` (an `optimum` group, no `posterior`) or as a separate typed
object with `to_datatree()`. The consumers that make it worth doing:

1. **Warm starts**: `Engine.initial_positions(count, around=optimum)` — a
   ball at the mode scaled by the covariance — cuts burn-in for every
   slot-A engine and gives NUTS a mass matrix.
2. **Laplace as the two-line composition** `optimise → curvature → draw`,
   sharing the optimum with the samplers instead of each engine finding
   the mode again.
3. **Model comparison on the cheap**: a Laplace evidence from the same
   objects, recorded under §5.1's attribute.

### 3.3 Bayesian optimisation — a different animal

BO (`BoTorch`/`Ax` on torch; `optuna` backend-neutral) optimises an
expensive black-box objective by fitting a GP surrogate to its evaluations.
Two readings for ampere, and they belong in different places:

- **BO as a MAP finder** for an expensive simulator: legitimate, fits slot
  A, produces an `Optimum` (§5.3) and, as a by-product, a GP surrogate of
  `log_prob`. Low value on its own — a mode of an expensive simulator is
  rarely what an astronomer wants from it.
- **BO as an acquisition strategy** — choosing *where* to simulate next so
  a training set (§13's budget) or an emulator (horizon (c)) is informative
  rather than prior-random. This is the reading with a future: it is the
  active-learning half of emulation and of multi-fidelity modelling
  (horizon (a) names multi-fidelity BO explicitly). Its natural home is
  beside `simulate_many` as a *proposal* (the thing TMNRE's truncated prior
  already is, W3.4), not as an engine. **Recommendation**: keep BO out of
  `ampere.inference` and treat it as a proposal/acquisition item in the
  emulation line of work.

## 4. Approaches not on Peter's list, worth a line each

- **SBC for every engine, for free.** `ampere.results.calibration.sbc`
  takes a duck-typed `engine_factory` (W3.6). Nothing in it is SBI-specific:
  an MCMC or Laplace engine can be rank-calibrated on the same toy problems.
  Any engine added under this memo should land with an SBC row in
  `tests/inference` — a stronger acceptance test than "agrees with emcee
  within 0.5σ", and it is the calibration study that would expose an
  approximate engine's bias (Laplace on a skewed posterior, mean-field VI's
  under-dispersion) as a number.
- **Importance reweighting between engines** (horizon (b)'s machinery,
  already hooked by the per-draw `log_prior`/`log_likelihood`): an
  approximate run (Laplace, Pathfinder, VI) can be corrected toward the
  true posterior by importance weights computed from the stored per-draw
  terms plus the approximation's own density — which `SBIEngine` already
  stores per draw for exactly this reason. Making the approximation's
  log-density a *stored variable* for every approximate engine (§5.2)
  turns every one of them into the proposal of a self-correcting importance
  sampler. This is the single most useful cross-cutting hook in the memo.
- **Emulator-accelerated MCMC** (horizon (c)): once an emulator is a
  `Model`, every engine above runs on it unchanged, and a differentiable
  emulator of a black-box simulator upgrades it to slot B. No engine work;
  noted so the ordering is visible — emulation is worth more than any
  single sampler here.
- **Delayed-acceptance / multi-fidelity MCMC** (horizon (a)): a cheap model
  screens proposals, the expensive one accepts. Fits slot A once fidelity
  tags exist. Later.
- **Profile likelihoods and frequentist intervals** come with `iminuit`
  (§3.1) and with the `Optimum` type; some users will want them beside the
  Bayesian answer for the same fit.

## 5. The adaptations — what the contract should take once

Each is small, each is a §4 change needing a decision-log row, and each is
what turns "one engine" into "one engine, no schema work".

1. **An engine-neutral evidence attribute.** `ampere_log_evidence` and
   `ampere_log_evidence_err` in `results.md` §9, written by every engine
   that estimates a marginal likelihood (nested sampling, SMC/PMC, PT,
   Laplace, snowline), with `ampere_evidence_method` naming how. Today only
   `ampere_dynesty_logz` exists. A reader comparing two archived fits
   should not need to know which sampler made them.
2. **The weighted-draw and approximate-draw rule, stated once.** (a)
   Weighted output is resampled to equal weight for the `posterior` group,
   the original count recorded, and the raw weighted output kept on the
   engine object — dynesty's convention promoted to the contract. (b) Every
   engine whose draws are not from the target stores the *proposal's* own
   log-density per draw (`sample_stats.proposal_log_density`, or a named
   variable per engine) — VI, Laplace, Pathfinder, SBI alike — so §4's
   importance correction is a post-processing function over stored runs.
   `SBIEngine` does this already; VI does not yet.
3. **An approximation attribute**: `ampere_approximation` in the root
   attrs, `"none"` for exact samplers and the family name otherwise
   (`"mean_field"`, `"multivariate"`, `"laplace"`, `"pathfinder"`,
   `"flow"`, `"density_estimator"`). VI's `vi_guide` and SBI's attrs stay;
   this is the one key a plot or a summary can check before it reports an
   R-hat that means nothing.
4. **An `Optimum` result** (§3.2) and the warm-start bridge
   `initial_positions(count, around=...)`. This is the only adaptation with
   real design in it — where the object lives (`ampere.results`), whether
   it is a `DataTree` or a dataclass, how profile intervals are shaped.
5. **Batched `log_prob` on slot A for `BATCHABLE` problems** — nautilus,
   ultranest and pocoMC all accept a vectorised likelihood, and
   `simulate_batched`'s machinery (W3.1 slice 2) already lowers a batch of
   θ natively on torch and jax. An optional `log_prob_batched(thetas)` on
   the problem, present when the problem is batchable, would make the
   population samplers as fast on a native problem as NUTS is. On the
   reference backend it is a loop, and the samplers lose nothing.
6. **Dependency tiers**: base (numpy-only) — nautilus, ultranest, iminuit,
   snowline; `jax` extra — blackjax, jaxns; `torch` extra — pocoMC, BoTorch.
   None belongs in the base install by default; nautilus and ultranest are
   the two light enough to argue for.

## 6. A ranking, if the question is "what first"

| Tier | Item | Slot | Cost | Why |
|---|---|---|---|---|
| 1 | `nautilus` + `ultranest` behind a shared `_nested.py` | A | S each, after §5.1–5.2 | Evidence and multimodality with fewer calls than dynesty; no new deps of weight |
| 1 | `VIEngine` guides `laplace`, `flow` (pyro/numpyro autoguides) | B | S | Two of Peter's three "other inference" items, nearly free |
| 1 | `blackjax` route on jax: MCLMC, Pathfinder (as initialiser and as approximation) | B | M | The best new sampler and the best initialiser, one small dependency, jax only |
| 2 | `Optimum` + warm start + `iminuit`/scipy/optax optimisers | A/B | M (design) | Unblocks MAP, Laplace-as-composition, profile intervals |
| 2 | `pocoMC` | A | S–M | SMC with flow preconditioning; needs the torch extra and §5.5 to shine |
| 2 | `snowline` | A | S | Laplace on the gradient-free half of the ladder |
| 3 | `eryn` PT, `nutpie`, `jaxns`, EKI | A/B | S–M | Overlap with the above; add on a use case |
| — | Bayesian optimisation | — | — | Not an engine: a proposal strategy for training sets and emulators (horizon (a)/(c)) |

The engine battery that should come with any of them: every new engine on
the conformance toy problems, SBC-ranked through `calibration.sbc`, and its
evidence (where it has one) checked against the closed-form
linear-Gaussian case W3.6 already uses.

**Landed at W5.14 (2026-09-22): tier 1's nested-sampling pair.**
`nautilus` and `ultranest` ship behind one extra each, over a shared
`ampere/inference/_nested.py`, and the battery above is
`tests/inference/test_nested.py` — with `DynestyEngine` in it as the
*control*, because three nested samplers agreeing with each other proves
nothing if all three are wrong the same way. Two adjustments the table did
not anticipate, both recorded in that module: the closed-form
linear-Gaussian case is run with **two** free parameters rather than W3.6's
one, because nautilus refuses a one-dimensional problem; and "one cost
record per run" (horizon (g)) is discharged by `engine_evaluations`, which
`Engine.finish` already writes for every run, beside each sampler's own
count of likelihood calls under its own name — no new engine-neutral
attribute, as §10's ruling that (e) to (g) need no contract change
requires.

## 7. Questions for Peter

1. Is an approximate engine's *stored proposal density* (§5.2b) wanted as
   a rule, given it makes importance correction and population reweighting
   a post-processing step? It is the hook with the widest consequences.
2. Should optimisation results be a `DataTree` (so archiving, provenance
   and plotting reuse everything) or a lighter typed object with a
   `to_datatree()`? The former keeps one file format; the latter keeps the
   `posterior` group honest.
3. Is BO wanted as an *engine* at all, or only as acquisition for
   emulation and multi-fidelity work (this memo's recommendation)?
4. Dependency policy: may `nautilus`/`ultranest` join the base install, or
   does every new sampler go behind an extra as `zeus` did?

## 8. Follow-ups from the second exchange (2026-09-10)

Peter's response to §§0–7, and Fable's answers. Where an answer revises the
memo above, the revision is stated here rather than edited in, so the
reasoning stays visible.

### 8.1 Population inference by importance reweighting over amortised SBI

Peter's use case: an amortised estimator (over noise and sampling, the
`horizon_notes.md` §4–5 programme) run over 10⁴–10⁹ sources, then a
population model built by reweighting each object's draws under a
population prior and fitting the hyperparameters by VI. This is design
horizon (b) with the object count taken seriously, and it sharpens three
things above.

- **Which estimator to amortise.** For population work the quantity per
  object is the marginal likelihood under the population prior,
  `p(xᵢ|Λ) = ∫ p(xᵢ|θ) π(θ|Λ) dθ`. With an amortised *posterior*
  `q(θ|xᵢ)` trained under an interim prior `π₀`, that integral is the
  reweighting estimator `(1/N) Σⱼ π(θᵢⱼ|Λ)/π₀(θᵢⱼ)` over draws
  `θᵢⱼ ~ q(·|xᵢ)`, and its variance grows as the population prior narrows
  relative to `π₀` (the effective sample size collapse every hierarchical
  reweighting paper warns about). With an amortised *ratio* estimator
  `r(x, θ) = p(x|θ)/p(x)`, the same integral is a plain Monte Carlo average
  of `r(xᵢ, θ)` over draws from `π(θ|Λ)` — no interim prior to divide out,
  and draws that follow the population prior wherever VI moves it. **NRE is
  the natural amortised estimator for population inference**, and W3.4's
  marginal estimators are already ratios. Storing per draw the estimator's
  own log-density (§5.2b) is what keeps the NPE route usable when NRE is
  not the estimator in hand; the two are the same hook.
- **Where VI actually fits.** §2.4 dismissed stochastic-gradient methods
  because single-object likelihoods are whole-dataset GP marginals. The
  population objective `Σᵢ log p(xᵢ|Λ)` is the opposite shape — a sum over
  10⁹ independent terms — and is exactly what stochastic VI over
  minibatches of *sources* was made for. pyro/numpyro's `plate` with
  `subsample_size` is the mechanism; the per-object terms are the stored
  draws or the stored ratio network. This is the one place in ampere
  where SVI's minibatch machinery earns its place, and it is a `population`
  module over archived outputs, not an engine.
- **Storage is the real design question at 10⁹.** One `DataTree` per
  source does not scale; what scales is the *estimator* plus a columnar
  store of a few draws (or summary statistics) per source — the
  training-set writer's shape (§11 of `results.md`) but partitioned
  (zarr or parquet, one row group per chunk of sources), with the
  provenance attrs once per file. The per-draw `log_prior`, the
  estimator's log-density and, for NRE, the log-ratio are the columns
  that make reweighting possible without the network. A design item of
  its own, before any population code.

### 8.2 Surrogate-posterior methods — §3.3 revised

Peter's point that BO has use cases "where even SBI is too expensive"
identifies the family §3.3 mis-filed: Bayesian optimisation *of the
log-posterior itself*, where the product is a GP surrogate of the density
rather than a mode. The candidates, all slot A, all producing a surrogate
from tens to a few hundred evaluations:

- **VBMC** (Variational Bayesian Monte Carlo; Acerbi 2018, 2020;
  `pyvbmc`): a GP surrogate of the log-joint refined by active sampling,
  with a variational mixture fitted to the surrogate. Returns an
  approximate posterior *and* an ELBO estimate of the evidence; handles
  noisy likelihoods. Typically a few hundred evaluations in ≤ 10
  dimensions. The most complete package of the family.
- **GPry** (El Gammal, Schöneberg, Torrado, Fidler 2022–23): a GP
  surrogate of the log-posterior with a bespoke acquisition, then MCMC on
  the surrogate for contours; built for cosmological likelihoods,
  reported at roughly 10² evaluations for ≤ 10-parameter posteriors,
  callable with a plain log-posterior and bounds.
- **BAPE / Bayesian active learning of posteriors** (Kandasamy et al.
  2015 and successors) and **BOLFI** (Gutmann & Corander 2016, in ELFI):
  the same idea for likelihood-free discrepancies — the SBI-adjacent
  member of the family.
- **Bayesian quadrature** (WSABI, Gunter et al. 2014; `emukit`): the
  evidence as the primary product, the posterior as the by-product.

Fit: each consumes `log_prob_unconstrained` on slot A (no gradients), and
each returns a surrogate plus an evidence — so the output goes through
§5.1 (evidence), §5.2 (draws from the surrogate, resampled or i.i.d., with
the surrogate's own log-density stored) and §5.3 (`ampere_approximation =
"gp_surrogate"`). No new contract surface beyond the memo's six. **Revised
placement**: tier 2, as *the* route for a simulator too expensive for SBI
budgets, with VBMC first because it ships evidence, noise handling and a
mature package; "BO as acquisition for training sets" (§3.3) stands
separately and is unchanged.

On the specific memory — methods giving 1-D marginals from ~10 evaluations
and 2-D marginals from ~50 — **Fable could not identify the paper with
confidence.** The numbers fit two shapes: (i) slice-and-interpolate schemes
that evaluate the density along axis-aligned 1-D and 2-D slices through the
mode and fit a GP or spline per slice (a profile-likelihood construction
with Bayesian dressing, which is only exact when the posterior is close to
separable); (ii) the surrogate family above in low dimension, where
GPry-style budgets do come down to tens of points per pair. If Peter can
recall an author or a keyword, the question is worth settling, because a
method that is (i) is a diagnostic rather than an inference engine, and
should be filed with the profile-interval tools of §3.1.

#### 8.2.1 Identified (2026-09-10): Rizzato & Sellentin 2022, arXiv:2203.05009

Peter's recollection was "Extremely expensive likelihoods: a
variational-Bayes solution for precision cosmology" (Rizzato & Sellentin
2022; MNRAS 2023). Checked against the paper: the variational family is a
generalisation of the **DALI** expansion (Sellentin, Quartin & Amendola
2014) — a Taylor expansion of the log-posterior about the mode carried to
cubic and quartic terms in a form that stays positive-definite, so the
density is non-Gaussian but parametric; the fit minimises a **quadratic
loss** between the family and the posterior *values at a fixed set of
points*, not an ELBO; **no gradients**; the points are whatever is already
in hand (an old chain, a simulation grid) rather than adaptively chosen;
**no evidence** is produced. The numbers Peter remembered are the paper's
2-D marginals — **14–45 evaluations** — while the full 7-D KiDS-450
posterior took ~18 000 real evaluations plus 8 900 artificial zero-density
points (0.6 % of the original chain), so the handful-of-evaluations claim
is per-marginal and the full-dimensional cost is much larger.

Fit for ampere, and one observation that makes it more interesting than
the surrogate family for one use case:

- It is slot A in the weakest sense — it needs log-posterior *values*, at
  points it does not choose — and its output is a parametric density, so
  it goes through §5.3 (`ampere_approximation = "dali"`) and §5.2 (draws
  by sampling the cheap analytic form, with its log-density stored).
- **Its "only at pre-selected points" setting is ampere's training set.**
  A §11 training set stores `(θ, ModelResult)` for a whole budget, and the
  log-likelihood at every stored θ is then a *likelihood* evaluation on
  stored model outputs — no simulator call. So this method is a
  post-processing consumer of a budget plus an observation, in the same
  position as SBI but returning a parametric posterior rather than a
  trained network, and it can be run on a budget that was simulated for
  something else. That is the reading with a use: a quick,
  gradient-free, no-training posterior from an existing budget, and a
  cross-check on an SBI posterior trained from the same one.
- Its marginal-only mode is diagnostic-grade — a 2-D contour from a few
  dozen points is a look, not an inference — and belongs with the
  profile tools of §3.1 rather than the engines.
- Placement: tier 3, "budget-reuse posterior", beside the surrogate
  family; VBMC/GPry remain the choice when the points *can* be chosen.

### 8.3 `margarine` (Bevins et al. 2023, MNRAS 526, 4613)

What it does: trains a masked autoregressive flow (or a KDE) on posterior
samples — weighted nested-sampling chains included — for a chosen subset
of parameters, giving a marginal density with nuisance parameters
integrated out. From that density it computes marginal KL divergences and
Bayesian dimensionality (the `anesthetic`-style statistics, but on a
subset), marginal Bayes factors (the nuisance-marginalised evidence
comparison it was written for), and reusable densities: a posterior from
one experiment becomes a prior or an importance-reweighting proposal for
the next, and independent experiments' marginal posteriors can be combined
by multiplication.

Fit for ampere:

- It is a **results-layer consumer**, not an engine: it wants draws (and
  weights, if the raw weighted output is kept — §5.2a) for a subset of
  merged names, which is exactly the `posterior` group with `var_names=`.
- Three of its uses map onto hooks this memo already asks for: a fitted
  density on stored draws *is* the proposal log-density of §5.2b for a run
  whose engine did not store one (any MCMC run), which makes importance
  correction and §8.1's reweighting available retroactively on archives;
  marginal Bayes factors are a model-comparison route (§8.4); combining
  archived runs is a multi-dataset route that costs no refit.
- A flow trained on a posterior as the **next fit's prior** is the
  interesting contract question: ampere's priors are per-parameter
  (`parameters.md` §4.1), so a joint, correlated prior over several
  parameters is an extension. It is not a large one — a flow is an
  invertible map from a base distribution, so it supplies both `log_prior`
  and a multi-dimensional `prior_transform` (the nested samplers' need)
  directly — but it touches the frozen §4.1 and is a decision, not a
  driver. This is also the "arbitrary priors" page the tutorials list as
  still to be written.
- **The package itself is a poor fit**: it is built on TensorFlow
  Probability, a third array library beside numpy, torch and jax, with the
  install weight that implies. The *functions* are worth having and are
  small to write over what the project already carries: `anesthetic`
  (numpy/pandas; the reference implementation of the nested-sampling
  statistics, weighted samples native) for the KL/dimensionality half, and
  a torch or jax flow from the stack the `sbi` and `jax` extras already
  install for the density-emulation half. Recommendation: implement
  `ampere.results.marginals` (working name) in that shape, and cite
  margarine for the method.

### 8.4 Model comparison — a fuller picture than §5.1

What exists: `dynesty`'s evidence; `pointwise_log_likelihood` (W2.8) feeding
`arviz.loo` and `arviz.waic`, so LOO/WAIC comparison is already a solved
route wherever the solver offers `conditional_loo` (refused by name where
it does not — the quasiseparable numpy solver); posterior-predictive checks
(W2.8). What is worth adding, all in `ampere.results` and none an engine:

- **`harmonic`** (McEwen et al.; the learnt harmonic-mean estimator, with
  normalising-flow targets since v1.1): computes the evidence from
  *posterior samples plus their log-posterior values*, which every ampere
  run stores per draw (`log_likelihood + log_prior`). It wants several
  chains for its cross-validation split — emcee/zeus give them; a
  single-chain run (VI, SBI, nested) can be split by draw. Dependencies:
  jax and flax, so the `jax` extra. **Fits the results layer verbatim**,
  and it is the only route to an evidence for an MCMC run after the fact.
- **Laplace evidence** from §3.2's `Optimum`; **thermodynamic integration**
  if a PT sampler lands; **bridge sampling** on stored draws plus the
  `log_prob` callable (a short implementation, no dependency); the
  **Savage–Dickey ratio** for nested models, which needs the marginal
  density at the nested value — §8.3's fitted density supplies it.
- The engine-neutral evidence attribute (§5.1) is what makes
  `compare(runs)` a one-line table: a `results.comparison` module with
  `log_evidence(tree, method="attrs" | "harmonic" | "bridge" | "laplace")`,
  `bayes_factor`, `savage_dickey`, and a `compare` that lines up evidences
  and LOO/WAIC the way `arviz.compare` does for the latter. Every method
  records how it got its number and its error estimate, because a Bayes
  factor without its method is not a result.

### 8.5 `ChainConsumer`, `anesthetic`, and "new figures made easy"

`ChainConsumer` (v1, 2023 rewrite) takes a pandas `DataFrame` per chain
with optional `weight` and `log_posterior` columns and draws the
astronomy-style corner (1σ/2σ contours, several chains overlaid), summary
tables (LaTeX-ready), and walk plots. Its strength over the six ampere
plots is **comparison**: several runs — engines, models, data cuts — on one
figure, which `arviz.plot_pair` does awkwardly. `anesthetic` (Handley)
uses the same idiom — a weighted `DataFrame` subclass — and adds
nested-sampling-specific plots and statistics.

Fit: one bridge serves both, and pandas users generally —
`to_dataframe(tree, var_names=..., group="posterior", weights=...)`
flattening `(chain, draw)` into rows keyed by merged name, with
`log_posterior` from the stored per-draw split, the equal weights (or the
raw nested weights, §5.2a), and the run's provenance carried as
`DataFrame.attrs`. Since the groups are xarray, `.to_dataframe()` exists
already; the bridge is naming, weights, `log_posterior`, and a `name` per
run for overlay legends. Both packages are light (pandas, matplotlib,
scipy). Recommendation: keep the six contract plots in ArviZ/matplotlib —
they are the stable, tested surface — and add the bridge plus a
`chainconsumer` optional extra for publication figures and multi-run
comparison, with one worked example on the docs site. The bridge also
makes "new figures" a pandas exercise, which is what most users will
reach for anyway.

### 8.6 What this exchange adds to §5's adaptations

Nothing new in kind; two in emphasis. §5.2b (the stored proposal density)
is now the load-bearing hook for three things — importance correction,
population reweighting at scale, and retroactive density emulation on
archives — and should be the first of the six taken. And a seventh item:
**a columnar, partitioned output format for amortised runs over many
sources** (§8.1), designed before the population module rather than after
the first 10⁷-source run shows one file per source does not work. **Ruled by Peter 2026-09-10**: the columnar store waits for
drafting until population work begins, on the condition that nothing
landed in the meantime excludes it — concretely, no results or
training-set format may assume one file per source, and the per-draw
columns that reweighting needs (`log_prior`, the proposal's log-density,
a ratio where one exists) must stay separable from the per-run
provenance. That constraint is recorded in the plan's design horizon (b).

## 9. The design horizons, re-read against the extensions (2026-09-10)

Peter asked whether the plan's four design horizons (§5 of
`DEVELOPMENT_PLAN.md`) need anything further in the light of this memo.
Checked against the contracts first, so that the answer is about what is
*not* yet reserved rather than what is:

- `parameters.md` §12.1 already reserves a `MultivariatePrior` slot beside
  `Prior`, and `ParameterSet.prior_transform`/`unconstrain` are set-level
  methods evaluated in dependency order — so a joint prior (a flow trained
  on a previous posterior, §8.3; a correlated calibration prior) is an
  extension point that exists, not a reversal. A flow prior would be that
  slot's first concrete customer.
- `results_schema.md` §583 already carries a `fidelity` tag on every
  container (horizon (a)'s hook), and `likelihoods.md`'s ruling X-1 gives
  every noise model the prediction as an argument — so an emulator's own
  predictive uncertainty (horizon (c)) can enter the likelihood as a
  `NoiseModel` without touching the family contract.
- Nothing registers engines: they are classes with a `NAME`, and
  provenance records the string. Third-party engines and the "engine
  battery" of §4 therefore have nothing to enumerate.

What each horizon should say in addition, and three horizons worth adding:

**(a) Multi-fidelity.** The inference-side consumers are delayed-acceptance
MCMC (cheap model screens, expensive model accepts), multi-fidelity
surrogates (Kennedy–O'Hagan, and VBMC/GPry's family in §8.2 at two costs)
and multi-fidelity SBI budgets. All three need **two problems, or one
problem with two model variants, over one parameter space**, and the check
that they share it is already `ampere_spec_hash` equality. Two things to
reserve: the training-set format carries fidelity per sample (a column, not
an attr — a budget will mix fidelities), and cost accounting is per
fidelity (below).

**(b) Population.** Covered by §8.1 and Peter's ruling on the store. Two
reservations to add: the interim prior must be reconstructible *from the
archive alone* (the population module must evaluate `π₀(θ)` for stored draws
without the model — `ParameterSet.to_spec()` is what it has, so the spec
must round-trip the priors, which `serialisation_review.md` should be
checked against rather than assumed); and NRE's log-ratio joins the
per-draw columns as a first-class stored variable, since it is the
prior-free quantity population inference wants.

**(c) Emulation.** Three inference-side additions. An emulator's identity
for `model_hash`/the artefact cache must include the training set it was
fitted on (its `ampere_spec_hash` and `ampere_model_hash`, plus the budget
and the emulator's own architecture), so a retrained emulator moves every
key that depends on it. An emulated draw must be *distinguishable* from a
simulated one in a run's provenance — `ampere_model_identity_hashes` names
the model class, which is enough if the emulator is its own `Model`
subclass, and that should be stated. And the surrogate-posterior family
(§8.2) and emulation are the same machinery at two levels — a GP over the
*density* versus a network over the *model output* — so the acquisition
strategy (§3.3's BO reading) belongs to a shared "proposal" abstraction
that `simulate_many(values=)` and TMNRE's truncated prior already
foreshadow.

**(d) Hierarchical SBI.** The W3.3 encoding is per problem, with one frozen
layout hashed from the observed data. A population-level amortised
estimator sees a *set of objects*, so the encoding must be allowed to
nest — a set embedding over objects, each object a set embedding over
samples — with a layout whose object count varies. The masked set
embedding is the right primitive; what to reserve is that `EncodingLayout`
is not assumed to have a fixed leading dimension.

**Three horizons to add**, each with a hook that costs nothing now:

- **(e) Model comparison and averaging.** Hook: the engine-neutral
  evidence attributes (§5.1) and the rule that every run stores per-draw
  `log_prior` and `log_likelihood` (already decided). Consumers: nested-
  sampling and SMC evidences, `harmonic`, bridge sampling, Savage–Dickey,
  LOO/WAIC (already there), and Bayesian model averaging over archived runs
  with the same data hash.
- **(f) Approximate-inference correction.** Hook: every engine whose draws
  are not from the target stores the proposal's own log-density per draw
  (§5.2b), and `ampere_approximation` names the family (§5.3). Consumers:
  importance correction of VI/Laplace/Pathfinder/SBI runs, population
  reweighting (b), retroactive density emulation on archives (§8.3). This
  is the cross-cutting hook the memo keeps returning to, and it deserves
  to be named once at horizon level rather than inside three others.
- **(g) Engines as a registry, with uniform cost accounting.** Hook: an
  engine registry shaped like the realisation registry (`register_engine`,
  `registered_engines()`), so a third-party engine — a user's own, or one of
  §§1–3's behind an extra — is discoverable by name, stamped into
  provenance by the same code, and enumerated by the engine battery (§4);
  and one cost record every engine writes (`engine_evaluations` exists for
  slot A; SBI records simulated/usable; a surrogate method records surrogate
  evaluations) — per fidelity where (a) applies — so the tier claims in §6
  and multi-fidelity cost models are measurable rather than asserted.

None of these needs a contract change today. (e)–(g) are one paragraph
each in the plan's horizon list, and the per-horizon reservations above
are sentences in the contracts they name, all of which can wait for the
first item that touches them — provided the list exists so that nothing
lands against it in the meantime, which is the same condition Peter set
for the columnar store.

## 10. Rulings (Peter, 2026-09-10)

§7.1 **yes** — the stored proposal density is the rule for every approximate
engine (it will make VI results interpretable too); §7.2 **`DataTree`** with
an `optimum` group; §7.3 **acquisition only** — BO is not an engine; §7.4
**every new sampler behind its own extra** — the base install stays quick
to start with, a user upgrades for a specific problem; §9's three horizons
**approved** and added to the plan as (e), (f), (g); §5's first three
adaptations **drafted now** as W5.0 (the results contract for approximate
and evidence-producing engines); §8.6's columnar store confirmed as
recorded. Nothing in this memo is scheduled beyond W5.0.

