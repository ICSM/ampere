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
