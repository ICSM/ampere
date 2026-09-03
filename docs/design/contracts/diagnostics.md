# Ampere v2 — Diagnostics Design Spec (W1.12)

Status: **approved by Peter, 2026-09-01** (the §10 defaults stand as
written; §11 records the one addition from his review — posterior
calibration as a future family). Not frozen until W1.13. Implements
`DEVELOPMENT_PLAN.md` §4.8. This is a **design document only** — no code
lands with this item. It fixes where the diagnostic families live, what
they consume and produce in the §4.2 container vocabulary, their
dependency/extras story, and what is in scope for Phase 2 versus deferred.
Exact class/function APIs are sketched for orientation, not frozen; the
binding API is written when the code lands (Phase 2 for post-fit families;
Phase 2 also for the pre-fit family, per §4.8's "1D implementation lands
with Phase 2").

Where this document and `DEVELOPMENT_PLAN.md` disagree, the plan wins and
this document is wrong — file it as a decision-log correction.

This document builds on `docs/design/prior_art.md` §4 (Starfish, lesson S4),
§5 (RHMF/Robusta-HMF, the adoptability assessment) and §6 (Tension 5), and on
`architecture.md` §3–4 (namespace/extras/lazy-import policy) and §8 (this
item owns "Diagnostics API"). It does not repeat prior-art research already
done there; it cites and resolves the open questions those sections left for
this item.

---

## 1. What this document is for

`DEVELOPMENT_PLAN.md` §4.8 names three diagnostic families united by one
purpose: **tell the user where a flexible (GP) likelihood is needed, rather
than leaving it as an act of faith.** They sit at three different points in
the fitting pipeline:

| # | Family | Pipeline stage | Needs a fit? |
|---|---|---|---|
| A | RHMF-style screening | pre-fit, collection-level | No — runs on raw data |
| B | Residual whiteness / posterior-predictive checks | post-fit, standard-likelihood | Yes — a plain iid-Gaussian fit |
| C | GP-localisation | post-fit, flexible-likelihood | Yes — a GP-noise-model fit |
| D | Posterior calibration (SBC / coverage) — **future, Phase 3** (§11) | validation of the inference itself | Yes — many fits of simulated data |

B is explicitly the trigger for "should I turn the GP on"; C is what you get
once you have. A is a cheaper, earlier, collection-level version of the same
question, answerable before committing to any fit at all. This document
specs each in turn, then resolves the one question `prior_art.md` §6
Tension 5 deliberately left open: whether A and C should present as one
visual language or two.

---

## 2. Family A — pre-fit, data-driven screening (RHMF-style)

### 2.1 Purpose

Given a *collection* of comparable spectra (or SEDs), flag which
features/regions, and which objects, a smooth/low-rank model is going to
struggle with — before any physical model has been evaluated. This is a
cheap triage step: run it before committing to expensive per-object fits, to
decide up front which objects or regions need the flexible likelihood.

### 2.2 Adoptability assessment: Robusta-HMF

`prior_art.md` §5 already did this assessment (fetched 2026-09-01, the same
day as this item — treated as current, not re-fetched); it is reproduced
here because W1.12's acceptance criterion asks for it explicitly.

| Question | Finding | Source |
|---|---|---|
| Package | `robusta-hmf` on PyPI; source `github.com/TomHilder/robusta-hmf` | §5 |
| Licence | **MIT** on the code (the arXiv paper text is separately CC BY 4.0 — governs the paper, not the code; do not conflate the two) | §5 |
| Language | **JAX** | §5 |
| Maturity | Small but real: 272 commits, 8 stars, 0 forks, has `tests/` and CI, ships `examples_paper/` reproducing the paper. Research-grade, one active group (Hogg's) — not a widely-adopted community package | §5 |
| API shape | scikit-learn-esque: `Robusta(rank=K, robust=True, robust_scale=Q)`, `model.fit(Y, W, max_iter=...)`, `model.synthesize()` | §5 |
| Verdict | **Adoptable as an optional dependency.** MIT licence and a real (if small) test suite clear the bar. Do not vendor or reimplement | §5 |

**This item's addition to that assessment**: the "if Robusta-HMF proves
adoptable" condition in `DEVELOPMENT_PLAN.md` §5 Phase 2 scope is satisfied —
this family is in scope for Phase 2 (1D case only; see §2.6).

### 2.3 Inputs / outputs, in §4.2 vocabulary

**Input.** A collection of §4.2 containers sharing a comparable coordinate
axis (e.g. a `DatasetCollection` of `Spectrum`-kind containers) — coordinates
(the shared feature axis, e.g. wavelength), values, and per-element masks.

**A real integration point, not assumed away**: RHMF factorises a *matrix*
`Y` with a fixed feature axis. §4.2 commits to "no regular-grid assumption
anywhere in the contracts" and irregular per-object coordinates in general.
Applying RHMF therefore requires an explicit **alignment step** — resampling
or interpolating each member of the collection onto one common coordinate
grid before factorisation — which is not RHMF's problem to solve and is not
free of assumptions (interpolation error, choice of grid). This is scoped as
part of the `ampere.diagnostics` adapter, not of RHMF itself, and is called
out as an obligation in §7.

**Masks.** RHMF's own missing-data convention — zero weight, equivalently
infinite uncertainty (prior_art.md §5 Lesson R1) — is exactly ampere's §4.2
mask semantics. No impedance mismatch: an ampere mask converts to an RHMF
weight column directly (mask → `w=0`).

**Output.**
- The raw `robusta-hmf` state for expert users: per-feature/per-object
  robust weights, the low-rank reconstruction, the loss history.
- A converted **anomaly score** (§4, shared convention) for the common case:
  coordinate-indexed, one value per feature per object (or aggregated
  per-object), higher = less well explained by the low-rank model, derived
  from RHMF's robust weights (prior_art.md §5 Lesson R2: "the paper
  aggregates them per-object... a low quantile across features").

### 2.4 Where it lives: `ampere.diagnostics` (new namespace)

Not `ampere.results`. Three independent reasons, any one of which would be
sufficient:

1. **Nothing has been fit yet.** `ampere.results` is defined (`architecture.md`
   §3, `DEVELOPMENT_PLAN.md` §4.6) as "ArviZ `InferenceData` emission + all
   plotting" — there is no `InferenceData` at this pipeline stage for this
   family to consume.
2. **The dependency is orthogonal to backend choice.** `robusta-hmf` pulls
   JAX regardless of whether the user is fitting with `reference`, `torch`,
   or `jax`. Coupling it to `ampere.results` (which every backend's fits
   flow through) would mean every results-consuming user pays for a JAX
   import path unless it is scrupulously lazy; a dedicated namespace makes
   the optionality structural, not just a matter of import-order discipline.
3. **`architecture.md` §8 explicitly assigns "Diagnostics API" to this item**
   without presupposing it folds into an existing namespace — the namespace
   diagram in `architecture.md` §3 predates this item and does not list
   `diagnostics/`. This document adds it; see the obligation in §7 to
   reconcile that diagram.

`ampere.diagnostics.rhmf` is the proposed module. It depends on `ampere.core`
containers only (never on `ampere.results`, `ampere.backends.*`, or
`ampere.inference`) — it can run before, and independently of, any fit or
backend choice.

### 2.5 Defaults vs expert opt-in

Prior_art.md §5 Lesson R3 is explicit: RHMF has two free hyperparameters
(`rank`, `robust_scale`) that the paper itself says need cross-validation for
reliable settings, with no automatic default. It poses the choice this item
must make: ship validated defaults, or scope screening as expert opt-in.

**Decision: expert opt-in, not a silent default, for the Phase 2 landing.**

Reasoning:
- This is a *design* item; no empirical validation against representative
  ampere data (SEDs, low-resolution spectra) has been done, or can be done,
  in a document-only work item. Claiming "validated defaults" here would be
  false.
- RHMF's own authors say the settings need cross-validation for reliable
  results. Shipping unexamined defaults presented as "just works" risks the
  exact failure mode this feature exists to catch: misspecification going
  undetected because the screening step itself was silently miscalibrated.
- The tool's value proposition (`prior_art.md` §5 Lesson R3) is "quick look
  before you commit to the expensive fit" — that framing survives an
  opt-in API; it does not survive a validation study becoming a
  prerequisite for every use.

Concretely: `rank` and `robust_scale` are required arguments (or an explicit
`rank="heuristic"` / `robust_scale="heuristic"` sentinel using a documented,
clearly-labelled-as-unvalidated heuristic — e.g. an explained-variance elbow
for rank — that emits a one-time warning naming it as unvalidated). There is
no silent default that runs without the user seeing that a choice was made.

**Phase 2 obligation** (recorded here, not resolved here): once Phase 2 has
representative ampere SED/spectrum collections to test against, validate
`rank`/`robust_scale` choices empirically and promote a validated default at
that point — this is a natural W1.13-or-later decision-log entry, not a
Phase 1 one.

### 2.6 Phase 2 vs later

- **Phase 2**: 1D case only — spectra/SED collections, matching
  `DEVELOPMENT_PLAN.md` §4.8's own scoping ("1D implementation lands with
  Phase 2") and §5's Phase 2 diagnostics line. The alignment/resampling
  adapter in §2.3, the expert-opt-in API in §2.5, and the shared anomaly
  score conversion in §4 are all Phase 2 deliverables.
- **Later / out of scope for now**: 2D+ collections (images, cubes) — RHMF's
  matrix-factorisation formulation would need a genuinely different
  flattening/alignment strategy that has not been designed; no commitment is
  made here about whether or when that happens.

### 2.7 Dependency / extras story

New extra: `ampere[diagnostics]` → `robusta-hmf` → JAX. Per `architecture.md`
§4's lazy-import policy:
- `ampere.core` is untouched — the policy already forbids this.
- `ampere.diagnostics.rhmf` may `import robusta_hmf` (and transitively
  `jax`) at module top level, because nothing reaches that module without
  wanting it — same carve-out `architecture.md` §4 rule 2 grants
  `backends/torch` and `backends/jax`.
- `ampere.diagnostics/__init__.py` itself must not eagerly import
  `robusta_hmf` — only submodules that need it do, consistent with rule 2's
  "any module reachable without opting in... must import lazily."
- Failures raise `ampere.core.exceptions.OptionalDependencyError` (the one
  shared type, per `architecture.md` §4 rule 3) — not a bespoke message.
- **Risk to watch, not resolved here**: `robusta-hmf` pins its own JAX
  version range; a user with `ampere[jax]` (backends/jax) *and*
  `ampere[diagnostics]` installed together needs those two JAX pins to be
  mutually satisfiable. This is the same class of problem `pyphot`/`sbi`
  pins are already tracking (`DEVELOPMENT_PLAN.md` §2) — Phase 2
  implementation should check it, not assume it away.

---

## 3. Family B — post-fit residual whiteness and posterior-predictive checks

### 3.1 Purpose

After a **standard-likelihood** fit (plain iid Gaussian noise — no GP yet),
test whether the residuals still show structure the smooth model failed to
capture. Structure found here is the concrete, quantitative trigger for
"switch the GP on," replacing "turn on the flexible likelihood as an act of
faith" with a diagnosed decision — the framing `DEVELOPMENT_PLAN.md` §4.8
opens with. `DEVELOPMENT_PLAN.md` §4.8 scopes this family explicitly to
"standard-likelihood fits" — it is deliberately not run against a
GP-augmented fit, whose residuals are by construction whitened by the GP
term (that is family C's territory).

### 3.2 Concrete statistics

- **Whiteness / autocorrelation — Ljung-Box-style.** Tests the null that
  standardised residuals `(predicted − observed)/σ` are uncorrelated across
  the ordered coordinate (wavelength, time, …); a significant statistic
  means leftover structure a smooth model can't explain — exactly what a
  Matérn-kernel GP term is built to absorb.
- **Posterior-predictive checks (PPC).** Draw replicate data
  `y_rep ~ p(y | θ)` for posterior draws `θ`, compare a discrepancy
  statistic `T(y_rep, θ)` against the observed `T(y_obs, θ)` (a Bayesian
  p-value in the sense of Gelman et al.). `T` can be a chi-square-type
  standardised-residual sum, or the Ljung-Box statistic itself used as the
  discrepancy measure, tying the two checks together rather than treating
  them as unrelated.

**A note on irregular coordinates.** Classical Ljung-Box assumes a fixed lag
structure (evenly-spaced data). §4.2 forbids assuming a regular grid. Two
honest options, neither adopted here (this is a Phase-2/W1.6 design
question, not settled by this document): (a) resample residuals onto a
nominal regular grid before testing (simple, lossy — loses exactly the
irregular-spacing generality §4.2 insists on), or (b) use a
separation-binned autocorrelation / structure-function statistic native to
irregular coordinate spacing (more faithful to §4.2, more design work, no
off-the-shelf implementation to point to the way `statsmodels`' Ljung-Box
is). Flagged as an obligation in §7.

### 3.3 Inputs / outputs, in §4.2 vocabulary

**Input**: coordinates and masks from the original §4.2 container (unchanged
from the fit); **values** are the per-point standardised residual at each
retained posterior draw.

**What `DEVELOPMENT_PLAN.md` §4.6 already gives for free, and what it does
not.** §4.6 requires every run to store per-sample `log_likelihood` — in
ArviZ's own convention this is per-observation, not summed to one scalar per
draw (that per-observation shape is exactly what LOO/WAIC, the requirement's
stated purpose, need). That makes **posterior-predictive discrepancy
statistics built from log-likelihood contributions cheap already** — no new
`InferenceData` field, no re-running the forward model. It does **not**, by
itself, give the Ljung-Box test what it needs: Ljung-Box requires **signed**
per-point residuals (autocorrelation is a sign-sensitive statistic), and for
a Gaussian noise model `log_likelihood_i` recovers `|residual_i|` up to sign
but not the sign itself — the sign (equivalently, the raw predicted value or
signed residual) is a small addition, not something §4.6 already stores.
This is flagged precisely, not smoothed over, as an obligation in §7 — this
document should not claim §4.6 makes the whiteness test itself free when it
only makes the PPC half free.

**Output**: a (statistic, p-value) pair, optionally per-chain or pooled; a
diagnostic plot (residual-vs-coordinate panel plus an autocorrelation/PPC
panel) via `ampere.results`.

### 3.4 Where it lives: `ampere.results`

This family consumes `InferenceData` exclusively (residuals/log-likelihood
from a completed run) and produces a plot — precisely `ampere.results`'s
existing remit ("ArviZ `InferenceData` emission + all plotting," §4.6). It
introduces no new heavy/optional dependency beyond what `ampere.results`
already needs for ArviZ (Ljung-Box's statistic is simple autocorrelation
arithmetic; a `statsmodels` dependency is a possible but not required
implementation choice, deferred to Phase 2 code review). Unlike family A,
there is no reason to isolate this behind a separate extra.

### 3.5 Dependency / extras story

No new extra. Lives inside whatever `ampere.results` already requires
(ArviZ, folded into the base install per `architecture.md` §3's extras
table once `ampere.results` lands). If a Phase-2 implementation chooses
`statsmodels` for the Ljung-Box statistic specifically rather than a
hand-rolled version, that becomes a new small required (not lazy) dependency
of `ampere.results` — a Phase-2 implementation decision, not fixed here.

### 3.6 Phase 2 vs later

Phase 2, following directly behind `ampere.results`'s `InferenceData`
emission (§5 Phase 2: "Results: everything emits `InferenceData`... 
Diagnostics (§4.8): post-fit residual tests"). No later-deferred component —
once per-observation log-likelihood and (per the §7 obligation) signed
residuals are available, this family is cheap and complete for the Gaussian,
1D-ordered case; irregular-coordinate Ljung-Box generalisation (§3.2) is the
one open sub-question left to Phase 2/W1.6.

---

## 4. Family C — GP-localisation from the flexible fit itself

### 4.1 Purpose

Once the flexible (GP) likelihood is enabled, its own posterior already
answers "where is the model deficient": `DEVELOPMENT_PLAN.md` §4.8 states
this directly — "the posteriors on GP amplitude and length-scale, and the
conditioned GP mean, already localise where the model is deficient." This
family surfaces that as a standard diagnostic, not a new computation — it is
the closest thing ampere has to Starfish's explicit local kernels
(`prior_art.md` §4 Lesson S3/Tension 4), delivered as a diagnostic on a
single global kernel rather than as trans-dimensional model structure.

### 4.2 Inputs / outputs, in §4.2 vocabulary

**Input**: posterior draws of the GP hyperparameters (amplitude,
length-scale) — ordinary parameters on the noise model's `Parameterised`
interface (`parameters.md` §13's own note to W1.6: "GP hyperparameters...
are ordinary parameters on a `Parameterised` noise model") — plus the
conditioned GP mean function, an output of the fitted `NoiseModel`
(`DenseGP`/`QuasisepGP`, §4.4) evaluated at the data coordinates (or a finer
grid for visualisation).

**Output**: coordinates (the data's own coordinate axis, or the finer
evaluation grid), values = the conditioned GP mean (signed — shows the
direction of the model's local deficiency, not just its magnitude) together
with its posterior uncertainty band, mask inherited from the data. An
aggregated **anomaly score** (§4, shared convention) is also produced:
coordinate-indexed, amplitude-weighted magnitude of the local deviation,
for the same shared presentation family A uses.

### 4.3 The global/local amplitude degeneracy caveat (mandatory)

`prior_art.md` §4 Lesson S4, carried forward to this item by name: Starfish
found its global kernel amplitude and its explicit local kernels trade off
against each other on real data — "little difference between the posteriors
in the third and fourth tests." Ampere's single global kernel does not have
separate local components to trade off against, but the analogous ambiguity
survives in a different form: **a large fitted GP amplitude with a short
length-scale can mean either "genuinely misspecified locally" or "the
kernel's smooth global component is under-amplitude and compensating
locally."** A large fitted amplitude **localises** deficiency (it tells you
*where*); it does not, by itself, tell you *why* — model error, an
underestimated noise budget, or a genuinely correlated astrophysical
process are all consistent with the same posterior shape.

This caveat is not optional prose buried in this spec — it is a
**mandatory, permanently-attached piece of the output**, not just of the
documentation:
- Every GP-localisation plot `ampere.results` emits must carry this caveat
  in its caption/docstring by construction (a plotting-function requirement
  for W1.8, not a "please remember to mention this" note).
- The anomaly-score container itself (§4) carries a `caveat` or
  `interpretation_notes` field for this family, populated automatically,
  so a user extracting the score programmatically (not just viewing the
  plot) still receives the warning.

### 4.4 Where it lives: `ampere.results`

`DEVELOPMENT_PLAN.md` §4.8 already says this outright — "surfaced as
standard plots in `ampere.results`." Consistent with family B: this family
consumes only `InferenceData` (the fitted GP's posterior) and the fitted
`NoiseModel`'s conditioning output, and produces plots. No new namespace
question here; the plan already settled it.

### 4.5 Dependency / extras story

No new extra. Requires whichever backend (`torch`/`jax`/`reference`) the fit
was run on, already an established dependency of having a GP-augmented fit
at all — `ampere.results` consuming its posterior adds nothing new.

### 4.6 Phase 2 vs later

Phase 2 (§5: "Diagnostics (§4.8): post-fit residual tests and GP-localisation
plots"), gated on the `DenseGP`/`QuasisepGP` `NoiseModel` implementations
(§4.4) and `ampere.results`'s `InferenceData` emission both landing — the
last of the three families to become usable in practice, since it needs the
most machinery underneath it (a working flexible-GP backend, not just a
standard-likelihood fit). No later-deferred component once those land;
2D+ conditioned-mean localisation (images/cubes) inherits whatever staging
the underlying GP solver strategies use (`DEVELOPMENT_PLAN.md` §4.4:
"Future strategies for 2D+... implemented in Phase 5") — this document does
not commit family C to anything beyond what the underlying solver supports.

---

## 5. The shared anomaly-score question (Tension 5, resolved)

`prior_art.md` §6 Tension 5 posed the question this document is asked to
settle: do pre-fit RHMF flags (family A) and post-fit GP-localisation
(family C) share one visual/output convention in `ampere.results`, or stay
visually distinct?

### Decision: shared convention, via a common container type — not a shared computation

Both families answer the same underlying question — "where does a
smooth/low-rank model fail to explain this data" — at different pipeline
stages, computed by unrelated statistics (RHMF's IRLS robust weights vs a
fitted GP's conditioned mean/amplitude). `prior_art.md` §5 Lesson R2 already
recommends this: "design the diagnostics module's output schema so both can
be plotted with the same conventions... even though the two statistics are
computed completely differently."

**Mechanism**: a lightweight `AnomalyScore`-shaped container — coordinates,
values (a documented, comparable range; higher = more anomalous/deficient),
mask, plus two fields that are *not* optional: `provenance` (which family
produced it, e.g. `"rhmf_prefit"` / `"gp_localisation_postfit"`) and
`interpretation_notes` (free text; this is where family C's §4.3 caveat
lives programmatically). This type is proposed for wherever W1.4 lands the
§4.2 container definitions (`ampere.core.results_schema` or equivalent) —
see the W1.4 note in §7 — precisely so that `ampere.diagnostics` (family A)
and `ampere.results` (families B/C) can both produce it **without either
namespace depending on the other**: `ampere.diagnostics` never imports
`ampere.results`, and `ampere.results` is not required to import
`ampere.diagnostics` (the JAX-optional namespace) just to render a family-A
score someone handed it.

One plotting function in `ampere.results` renders any `AnomalyScore`
regardless of provenance, with a consistent colour/axis convention (a
sequential scale, documented range, same sense of "higher = worse" in both).

### Why shared, not distinct — and how the "overstating comparability" risk is actually addressed

Tension 5 flagged a real risk in the shared option: conflating two
differently-computed statistics under one colour scale could overstate their
comparability. This document's answer is not to ignore that risk but to
address it at the right layer: **the visual grammar is shared (so users get
one consistent "where's the problem" reading experience across the
pipeline, lowering the cognitive cost of moving from "should I screen this
collection" to "did my fit's residuals agree"), while comparability is
policed by the mandatory `provenance`/`interpretation_notes` metadata, not
by the plot style.** Two `AnomalyScore` panels are never allowed to be
captioned or overlaid as if interchangeable without their provenance shown —
that is a plotting-function requirement for W1.8 (§7), not a matter of
visual style. This is more robust than "keep them visually distinct," which
protects against conflation only until someone screenshots two panels
side by side without labels; provenance metadata travelling with the data
protects against conflation wherever the container goes, including outside
`ampere.results`'s own plots.

**What is deliberately *not* shared**: family A's raw RHMF outputs (weight
matrix, low-rank reconstruction) and family B's (statistic, p-value) pairs
remain in their own, unconverted form for expert users — the `AnomalyScore`
conversion is an additional, simplified view, not a replacement for the
underlying method-specific output. Family B is not folded into the
`AnomalyScore` convention at all: a Ljung-Box statistic/p-value is not a
coordinate-indexed deficiency map the way A and C's outputs are — forcing it
into the same container shape would be the conflation Tension 5 warned
about, applied to a case where it genuinely does not fit.

---

## 6. Module placement summary

| Family | Module | New extra? | New heavy dependency? |
|---|---|---|---|
| A — RHMF pre-fit screening | `ampere.diagnostics.rhmf` (new namespace) | `ampere[diagnostics]` | JAX (via `robusta-hmf`) |
| B — residual whiteness / PPC | `ampere.results` | none | none beyond `ampere.results`'s own (ArviZ; optionally `statsmodels`, Phase-2 decision) |
| C — GP-localisation | `ampere.results` | none | none beyond the fit's own backend |
| Shared `AnomalyScore` container | `ampere.core` (**landed at the freeze** — ruled 2026-09-03, `results.md` §15 R4; `ampere.core.results_schema.AnomalyScore`) | none | none — plain numpy/scipy, per `architecture.md` §4 rule 1 |

---

## 7. Obligations on downstream work

**On W1.6 (likelihoods contract)**:
- Confirm GP amplitude/length-scale are declared exactly as `parameters.md`
  §13 anticipates (ordinary parameters on a `Parameterised` `NoiseModel`),
  since family C's input is those posteriors.
- Decide the irregular-coordinate whiteness-test strategy (§3.2): resample
  onto a nominal grid, or a separation-binned native statistic. This
  document deliberately leaves it open rather than picking one without a
  §4.4 `NoiseModel`/solver-level view of what's cheap.
- State explicitly (as `prior_art.md` §4 Tension 4 already asks W1.6 to do)
  that GP-localisation is ampere's substitute for Starfish-style local
  kernels, not an oversight — family C's spec here depends on that framing
  being on record in the likelihoods contract too, not only here.

**On W1.8 (results contract / plot API)**:
- Reserve `InferenceData` fields (or a documented derivation path) for
  **signed** per-point residuals, not just per-observation log-likelihood —
  §3.3's precise gap: §4.6's log-likelihood requirement makes PPC cheap but
  does not by itself give Ljung-Box its sign-sensitive input.
- Decide whether posterior-predictive replicate draws (`y_rep`) are stored
  by default in every run's `InferenceData` (memory cost: N_draws × N_obs
  per run) or computed on demand when a user requests family B's PPC
  (cheaper by default, a small re-evaluation cost only when asked for) —
  this document's inputs/outputs assume the latter but does not bind W1.8
  to it.
- Own the `AnomalyScore` plotting function (§5) as a first-class stub: one
  renderer, both provenances, the mandatory provenance/caveat display
  requirement from §4.3 and §5 built in from the start, not bolted on.

**On Phase 2** (recorded here as a checklist, not resolved here):
- Family A: the alignment/resampling adapter (§2.3), the expert-opt-in
  hyperparameter API (§2.5), the JAX-pin compatibility check against
  `ampere[jax]` (§2.7), and — later, not at Phase-2 landing — the validation
  study that could eventually promote a default `rank`/`robust_scale`.
- Family B: the Ljung-Box irregular-coordinate strategy W1.6 picks (above),
  and the dependency choice (hand-rolled vs `statsmodels`) for the
  statistic itself.
- Family C: nothing beyond what W1.6/W1.8 already owe it — it is gated on
  their landing, not on new design work here.
- **Also flagged for W1.4 (results schema)**: the `AnomalyScore` container
  proposed in §5 needs a home in whatever W1.4 lands for §4.2 containers.
  This document proposes it belongs in `ampere.core` (plain numpy/scipy,
  no optional dependency, so both `ampere.diagnostics` and `ampere.results`
  can produce/consume it without a cross-namespace dependency) but does not
  bind W1.4's exact class hierarchy — reconcile at W1.13 spec assembly if
  W1.4 lands with a materially different container design.
- **Also flagged for `architecture.md`'s namespace diagram (§3)**: it does
  not currently list `ampere.diagnostics`. This document adds the namespace;
  the diagram should be updated to match at W1.13 spec assembly (a §4
  contract-affecting addition, per this repo's ground rule 9, belongs in a
  decision-log entry in `DEVELOPMENT_PLAN.md` in the same PR that actually
  adds the namespace in code — not asserted unilaterally here).

---

## 8. Decisions and their reasoning

| Decision | Reasoning |
|---|---|
| Family A lives in a new `ampere.diagnostics` namespace, not `ampere.results` | Pre-fit: no `InferenceData` exists yet to consume; its JAX dependency is orthogonal to backend choice, not to results consumption (§2.4) |
| Families B and C live in `ampere.results` | Both are pure `InferenceData`-in, plot-out — precisely `ampere.results`'s existing remit; introduces no dependency `ampere.results` doesn't already carry (§3.4, §4.4) |
| RHMF hyperparameters are expert opt-in, not a shipped default, for Phase 2 | No ampere-specific validation exists or can exist in a document-only item; RHMF's own authors say the settings need cross-validation; an unvalidated "just works" default risks silently miscalibrating the very screening step meant to catch misspecification (§2.5) |
| Pre-fit and post-fit anomaly scores share one container type and plot convention | Both answer the same question at different stages (prior_art.md §5 Lesson R2); a shared visual grammar lowers cognitive cost across the pipeline (§5) |
| Comparability risk (Tension 5) is policed by mandatory provenance metadata, not by keeping the two visually distinct | Metadata travels with the data wherever it's consumed; visual distinctness only protects the one plot function that respects it (§5) |
| Family B is explicitly scoped to standard-likelihood fits, never GP-augmented ones | `DEVELOPMENT_PLAN.md` §4.8's own wording; GP-augmented residuals are whitened by construction, so testing them for whiteness is close to circular (§3.1) |
| Family C's amplitude/length-scale posteriors carry a mandatory, machine-readable caveat field, not just prose | `prior_art.md` §4 Lesson S4's degeneracy is a real interpretation risk; a caption a user can silently crop out of a screenshot is not durable protection against over-interpretation (§4.3) |
| `AnomalyScore` proposed for `ampere.core`, not `ampere.results` or `ampere.diagnostics` | Neither namespace should depend on the other; a plain-numpy container in core is producible/consumable by both without a cross-dependency, and complies with `architecture.md` §4 rule 1 (core imports no optional dependency) (§5, §7) |

---

## 9. Deliberate limitations of v1 (this spec)

Each of these is a decision, not an oversight. Each has an extension point.

1. **Family A is 1D-only.** RHMF's matrix formulation and this spec's
   alignment adapter (§2.3) are not designed for images/cubes. Extension
   point: a 2D+ alignment strategy would be a new adapter module inside
   `ampere.diagnostics`, not a change to this document's namespace decision.
2. **No RHMF default hyperparameters ship in Phase 2.** §2.5's opt-in
   decision is explicit and intentional, not a gap to fill immediately.
   Extension point: a validated-default decision-log entry once Phase 2 has
   representative data to validate against.
3. **The irregular-coordinate Ljung-Box strategy is unresolved.** §3.2 names
   two honest options and picks neither. Extension point: W1.6 decides,
   informed by whatever the `NoiseModel`/solver contract already does for
   irregular 1D coordinates (quasiseparable solvers handle this natively —
   the whiteness test should reuse that thinking, not invent its own).
4. **`AnomalyScore` is a container shape, not a frozen API.** This document
   fixes what fields it must carry (coordinates, values, mask, provenance,
   interpretation notes) and where it lives, not its exact class definition
   — that is W1.4's to write when the §4.2 container hierarchy lands.
5. **Family B's PPC replicate-draw storage question is flagged, not
   answered** (§7). This document assumes on-demand computation for its own
   inputs/outputs description but does not bind W1.8 to that choice.
6. **No commitment on family C for 2D+ data.** It inherits whatever staging
   the underlying GP solver strategy uses (Phase 5 for 2D+ per §4.4); this
   document does not attempt to get ahead of that.

---

## 10. Open questions for review

1. Is `ampere.diagnostics` an acceptable new top-level namespace, or would
   Peter prefer it nested under an existing one (e.g.
   `ampere.results.diagnostics`, trading the "no InferenceData dependency"
   argument in §2.4 against fewer top-level namespaces)? This document's
   default is a peer namespace to `results/`, not nested under it, precisely
   because family A has no `InferenceData` dependency at all — but this is
   a naming/organisation call, not a load-bearing technical one, and is
   flagged for review rather than treated as settled.
2. Does the expert-opt-in framing for RHMF hyperparameters (§2.5) undersell
   the "quick triage" value proposition too far? An alternative not adopted
   here: ship a documented, clearly-labelled-experimental heuristic default
   that the user must explicitly request (rather than requiring explicit
   values with no default path at all) — closer to what §2.5 sketches as
   `rank="heuristic"`, but the exact API surface is Phase 2's to write.
3. Should the `AnomalyScore` shared-convention decision (§5) be revisited if
   W1.4 lands a container hierarchy that makes a plain-numpy `ampere.core`
   addition awkward (e.g. if all §4.2 containers turn out to carry more
   structure than this lightweight type wants)? Flagged for W1.13
   reconciliation rather than pre-empted here.

---

## 11. Family D — posterior calibration (recorded at review, 2026-09-01; future scope)

Raised by Peter at review: is simulation-based calibration, or other
posterior-calibration diagnostics, useful — especially for SBI? **Yes**, and
it is a genuinely distinct axis from families A–C, which is why it is
recorded as a fourth family rather than folded into one of them: A–C
diagnose the *model* (where is a smooth model, or the fit's noise budget,
inadequate for this data); calibration diagnoses the *inference machinery*
(does the posterior the engine produces actually have the coverage it
claims), independently of whether the model is right.

The candidates and their hooks:

- **Simulation-based calibration** (SBC; Talts et al. 2018,
  arXiv:1804.06788): draw parameters from the prior, simulate data, fit
  each simulated dataset, and check that the rank of each true parameter
  among its posterior draws is uniform. Its ingredients are exactly what
  the §4.5 contracts already provide — `simulate(params)` and
  `prior_transform`/`sample` — which is why it costs no new contract
  surface, only compute.
- **Expected-coverage / TARP-style tests** (e.g. Lemos et al. 2023,
  arXiv:2302.03026): the SBI-era refinements of the same question, testing
  whether credible regions contain the truth at their nominal rate.
- The `sbi` package ships SBC and coverage diagnostics in its own
  `diagnostics` module, so for the NPE/NLE/NRE path this family is largely
  an integration, not an implementation — the same
  depend-don't-reimplement posture §2.2 takes for RHMF.

**Why it matters most for SBI** (Phase 3): an amortised neural posterior
can be silently overconfident in a way no residual test detects — the fit
to the *observed* data can look excellent while the posterior's claimed
uncertainties are fiction. Calibration is the diagnostic that catches
this, and the SBI literature increasingly treats it as mandatory
reporting. For MCMC/nested paths it is a heavier, optional check (each
rank statistic costs a full fit), useful when validating a new likelihood
or noise-model configuration — for instance the flexible-GP likelihood
itself in milestone M2.

**Scope and placement**: lands with **Phase 3's SBI layer**, not Phase 2 —
it needs `simulate()` and the inference engines to exist first. Placement
follows this document's own rules: it consumes `InferenceData` (many of
them) plus `simulate()`, carries no new heavy dependency beyond what the
SBI extra already brings, and so belongs with `ampere.results` or the SBI
module's own diagnostics — decided when the code lands, per the §8 table's
logic. Its outputs (rank histograms, coverage curves) are not
coordinate-indexed deficiency maps, so — like family B, and for the same
reason — it does **not** adopt the `AnomalyScore` convention.

**On evolution generally** (Peter's second point): diagnostics are
expected to grow as the field produces new ones. This document's structure
is the extension template — a new family states its pipeline stage, its
inputs/outputs in contract vocabulary, its module placement per §8's
reasoning, its dependency story per `architecture.md` §4, and whether it
adopts or (with justification) declines the shared `AnomalyScore`
convention. Adding a family is an ordinary documentation-plus-code change,
not a contract change, so long as it follows that template; nothing in
§4.8's three named families was ever a closed list.
