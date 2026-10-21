# Amortised evidence and honest amortised posteriors: the ENFR/CF-NPE note against the code, and where it sits in the plan

*Drafted 2026-10-21 by Fable at Peter's request, reading "ENFR and CF-NPE
for ampere — design note" (Peter, 2026-10-20) against the merged code at
`432b256` (master = origin/master, PR #82 merged 2026-10-16) and against
`DEVELOPMENT_PLAN.md`, `WORK_ITEMS.md`'s ruled Phase 8, the horizon memo
(`horizon_beyond_phase7.md` §6, §9) and the user-journeys memo (§6). The
note itself is not in the repository; this memo quotes what it needs.
Nothing here changes a ruled item — §6 lists the decisions and §7 drafts
the items for Peter to lift once ruled.*

## 0. The short answer

Both proposals are compatible with ampere v2 as built, and more cheaply
than the note assumes. Every piece the note marks ASSUMED exists, three
of its "open questions" are already answered by the code (the simulator
*does* draw correlated GP residuals for every kernel, on the reference
path and natively; the locked `sbi 0.27` ships `FMPE` and the locked
`zuko 1.6` ships `CNF` and `NSF`, so no flow is hand-rolled; failures are
already counted by reason), and no frozen §4 contract changes for either
method. Two of the note's framings need correcting: "one more method in
an existing engine" understates CF-NPE, because `SBIEngine`'s `METHODS`
table dispatches to `sbi`'s own trainers and CF-NPE is a trainer of its
own; and `ampere.comparison` would be a new top-level namespace the
plan's §3 layout does not have, when the plan already owns model
comparison at W8.6 (`ampere.results.compare`, ruled D8) — the right shape
is an *evidence engine* beside `SBIEngine` in `ampere.inference`, writing
the engine-neutral evidence triple that W8.6's table then consumes.

Placement: **not Phase 8**, which is ruled, balanced and un-opened, and
whose customers for these methods (the 133-disc population of W8.9, the
`LatentProfile` of W8.2, `compare` of W8.6) land late in it; **Phase 9**,
whose theme — hierarchy at scale and the SAGE/GRAMS amortised-SBI
flagship — is exactly the regime both methods were built for: many
objects through one instrument, tens to hundreds of parameters, a fixed
bank of expensive simulations reused for many epochs. Two one-sentence
hooks go into Phase 8 now (W8.6's and W8.3's texts), one small filler is
offered (`reduced_chi2` in the calibration group), and the rest is a
memo-first workstream of four items drafted in §7: the feasibility memo
with a runnable reproduction of the paper's Gaussian-mixture benchmark
(Fable, between waves), `SBIEngine(method="cf-npe")`, the evidence
engine, and the "physics + GP" comparison study that merges with the
Phase 10 noise-model selection study the horizon memo already planned.
The cross-fitted ENFR variant is the note's own hypothesis with no
evidence anywhere; it stays a switch decided by the ablation.

## 1. Where the project stands (2026-10-21)

- Phase 7 closed 2026-10-15 with `ampere-astro 1.0.0b2`; Phase 8 is drafted
  and fully ruled (W8.0–W8.16, D1–D12) and **waits on Peter's word to
  open**. No agent, worktree, gate or detached run is in flight.
- Since the handoff: PR #82 (the CASSIS file for the star-disc example)
  merged on GitHub 2026-10-16, so **D10's trigger has fired** — W8.10 may
  be dispatched when Phase 8 opens; master was two commits behind and was
  fast-forwarded today. Two PRs are open: #83 (an updated `HD105_SED.vot`,
  2026-10-16 — the same example's data; W8.10 or W8.13 should take it the
  `tests/data` way) and #85 (dependabot, `setup-pixi` 0.11.2).
- `AGENTS.md` carries an uncommitted copy of `CLAUDE.md`'s "Phase 7
  complete / Phase 8 ruled" paragraph — the Codex mirror the closing
  session updated but did not commit; harmless, one commit away.
- `examples/examples_paper/flexible_likelihood_comparison.py` is Peter's
  paper work, untracked by design.

## 2. The note's assumptions against the code

| The note says (ASSUMED or asked) | What is true | Where |
|---|---|---|
| `simulate_many()` with `ProcessExecutor`, per-run timeouts, counted failures | Yes, on `FittingProblem`: `executor`, `chunk_size`, `as_chunks`, `native`, `sharder`, `context`, a `training_set=` sink on the engine | `ampere/core/dataset.py:3424`; `_sbi.py:1498` |
| Open question 4: does the simulator draw correlated GP residuals, for every kernel? | **Yes, already.** `draw_observation` delegates to the family's `sample` with the same `NoiseParams` the density scores; `GaussianFamily.sample` adds `L z` through the bound solver's `latent_transform` for any kernel the solver takes (`Sum`, `Product`, `SpectralMixture`, `WarpedKernel` included); the joint noise model draws its `(n, T)` block the same way; torch and jax draw natively and refuse by name what they cannot. Hyperparameters are drawn from their priors because the likelihood component is part of the joint prior. No `draw_residuals` hook is needed. | `dataset.py:1462`; `likelihood.py:4299` (lines 29–40 of `sample`), `:3732`; `backends/torch/problem.py:1780`, `backends/jax/problem.py:1557` |
| "Scales linearly like the likelihood — unconfirmed" | The draw goes through the same solver as the density (`latent_transform` is the whitening the latent declaration uses), so `QuasisepGP` draws are linear in N. Confirm with one timing row at W9.A, not by assertion. | `likelihood.py:4299` |
| Failed draws dropped and counted, by kind | Yes: `ampere_sbi_failures` and per-reason counts; the context path records which context a crash was in | `_sbi.py:1601`; contract `inference.md:2391` |
| `ampere/inference/sbi.py`, `ampere/results/calibration.py` | `_sbi.py` (3,343 lines); `calibration.py` exists | — |
| `method=` is "one more value" | `METHODS = {"npe": "NPE", "nle": "NLE", "nre": "NRE", "tmnre": "NRE"}` maps to **`sbi.inference` trainer classes**, with `_DEFAULT_ESTIMATORS`, `_NETWORK_KEYWORDS` per method and TMNRE as a wrapper around `NRE`. CF-NPE is a trainer of its own (K·M regressors plus a residual density) with its own posterior object that must satisfy `calibrate()`'s `sample`/`log_prob` surface. It is a **new family at the dispatch**, sized M–L, not S. | `_sbi.py:288–315`, `:2546` |
| `calibrate()` works unchanged for any amortised posterior | Yes if the posterior exposes `sample`/`log_prob`; returns `ranks`, `coverage`, `ks_pvalue`, `c2st_ranks`, `tarp_coverage`. `reduced_chi2` is absent and is a results-side addition to the `calibration` group (`results.md` §4 shapes the group; check whether an added variable bumps the schema). | `_sbi.py:2097–2352`; `results/calibration.py:273` |
| Seeds: torch seeded from the problem's sub-stream, recorded | Yes: `_seeded(torch, "sbi.torch")` records `ampere_sbi_torch_seed` | `_sbi.py:1898` |
| `ContextPrior` with `ScaledSigma`, `SigmaArchive`, `SignalToNoise` | Yes, a protocol with those three | `core/simulate.py:955–1224` |
| `layout="set"` carries log σ; `"flat"` does not; refuse the combination by name | Yes; `SET_EMBEDDINGS = ("set", "transformer")`; a context with a flat layout is the pattern already refused | `_sbi.py:116–121, 355–363` |
| `ArtefactStore` keyed on problem hashes; a miss names what changed | Yes; **but** the posterior key includes `data_hash` because a stored posterior is conditioned on its observation (W3.13's correction). An ENFR ensemble is *amortised* and must key on the training-set key (`spec_hash` + `model_hash` + context + budget), not the posterior key. | `results/artefacts.py:435, :603`; plan §2 row W3.13 |
| The evidence triple `ampere_log_evidence`, `_err`, `_method` | Yes, written by the three nested engines with `"nested_sampling"`; VI and Laplace deliberately do **not** write it because a bound is not an evidence. A learned estimate with an ensemble error is the case the second attribute was created to distinguish (`inference.md:1154–1159` asks exactly whether it should name the family). | `inference.md:1154, :1231` |
| Nested-sampling ground truth on held-out simulations needs re-observed problems | `replace_observations()` already builds a problem on a simulated observation (the generic SBC route); `validate()` is that plus the nested engine, per simulation | `results/calibration.py:404` |
| A hand-rolled flow-matching density with fixed-step Euler and per-dimension divergence | **Unnecessary.** The `sbi` extra's lock already carries `zuko 1.6.0` (`CNF` with exact or Hutchinson trace, `NSF`, `MAF`) and `sbi 0.27.0` ships `FMPE` and `NPSE`. The reference density is a `zuko` flow; the depend-don't-reimplement rule (`diagnostics.md` §2.2) applies. | `pixi.lock:8470`; checked in the `sbi` environment today |
| Torch only; a jax port "would mirror the other backends" | SBI is torch-only by ruling; jax-native SBI (sbijax/flowjax) is an open §6 deferral. Both methods live behind the `sbi` extra, as `SBIEngine` does. | plan §6 (`:1350`) |
| M2's "113 posterior standard deviations" and the 1-parameter lognormal's C2ST 0.575 → 0.875 | Both quoted correctly | `m2_misspecification.rst:219`; `sbi.rst:381–384` |
| "ampere compares models only through its three nested engines" | True today; **W8.6 is ruled** (D8: bridge sampling, Savage–Dickey, the Bayes-factor table, model averaging, `harmonic` behind an extra). ENFR is the simulator-only row of that table, not a replacement for it. | `WORK_ITEMS.md` W8.6; plan §2 row "Phase 8 opened" |

## 3. Compatibility, by proposal

### 3.1 CF-NPE as `SBIEngine(method="cf-npe")` — compatible; M–L, not S

The engine's shape fits: the encoding, embeddings, context, seeding,
failure accounting, the artefact store and `calibrate()` all carry over,
and the note's refusals (folds below 2, truncation options, a context
with the flat layout) follow the house pattern. What the note understates
is the trainer: `SBIEngine` delegates training to `sbi`'s classes and
CF-NPE has no `sbi` trainer, so the item writes `_cfnpe.py` (stage 1's
K·M regressors on a shared frozen feature map or per-stage encoders,
stage 2's diagonal Gaussian or `zuko.NSF` over standardised out-of-fold
residuals conditioned on the features and the stage-1 mean, the
`CFNPEPosterior` adapter), threads its provenance (`ampere_sbi_cf_folds`,
`_replicates`, `_variant`, `_stage_seconds`) into the cache key, and
`vmap`s small stage-1 networks with `torch.func` as the batched path
already does. The contract question is nil: `inference.md`'s SBI text
names the `sbi` families as what landed, not as a closed list, and
`ampere_sbi_method` gaining a value is provenance.

The honest caveat is the regime. The paper's gain is zero at two
parameters and only clear from about a hundred; a 1-parameter lognormal
or a five-parameter SED is not where CF-NPE pays. ampere's customers are
W8.9 stage 2 (133 members plus hyperparameters), `Population` fits in
general, W8.2's `LatentProfile`, images, and the Phase 9 flagship. The
item's acceptance must therefore measure on one of those — the paper's
Gaussian–Gaussian reproduction pins the mechanics; a many-parameter
ampere problem pins the benefit (D5).

### 3.2 ENFR — compatible; the namespace and the home are the decisions

Everything ENFR consumes exists (§2). The design questions are
ampere's, not the paper's:

- **Where it lives.** The note proposes `ampere.comparison`. The plan's
  §3 layout has no such namespace and treats new namespaces as
  deliberate; model comparison is already owned by W8.6 in
  `ampere.results.compare`, which is a *consumer* of evidences. ENFR is a
  *producer*: it simulates, trains networks, caches an ensemble and emits
  a `DataTree` carrying the evidence triple — the shape of `SBIEngine`,
  not of a results module. Recommended: an `EvidenceEngine` (name open)
  in `ampere.inference` taking a mapping of named problems sharing a data
  hash, with `evaluate`, `bayes_factor`, `gof_pvalue` and `validate`; its
  evidences land in W8.6's table with `ampere_evidence_method = "enfr"`
  rendered beside `"nested_sampling"` and `"bridge"`, and the rung of the
  validation ladder the ensemble reached recorded as an attribute the
  table shows (D2). `ampere.comparison` as a thin front door re-exporting
  both is possible later, on the `ampere.interferometry` pattern, if the
  user-facing story wants one name.
- **What it writes.** The triple, with `_err` the ensemble scatter;
  `ampere_enfr_*` for the rest, on the `ampere_sbi_*` convention
  (`results.md:788`: flat `ampere_`-prefixed root attributes). Writing
  the triple from a learned estimator is the first time
  `ampere_evidence_method` carries a second value, which is the question
  `inference.md:1156` left open; W8.6 should settle the vocabulary before
  ENFR arrives (§4.1).
- **The reference density.** `zuko` already in the lock; the memo of
  W9.A measures a flow-matching `CNF` against a maximum-likelihood `NSF`
  as the reference, because the paper's argument for flow matching is
  about sample quality and ampere needs the density (D4).
- **The data dimension.** The paper's accuracy degrades past about 500
  dimensions; a flat 20,000-point spectrum is far past it. The set
  embedding is the route, and then the evidence is of the summaries — and
  an embedding trained for posterior inference may discard what
  discriminates models (the note's own open question). An
  evidence-specific embedding trained with the classifier, as the paper's
  permutation-invariant networks were, is the default to try; the
  "summary-space caveat" is then a stated attribute, as the never-silent
  rules ask.
- **Cost.** The paper used about 10⁵ simulations per model and ensembles
  of five; rung 1 of the ladder is 50 nested runs per model. For the
  benchmarks this is cluster work on `scripts/cluster/run.sh`'s pattern
  and sized L with the cluster run on Peter's word, as W8.7 is (D9). For
  ampere's own problems the amortisation is the point: one trained
  ensemble per instrument serves a survey, which is the Phase 9 flagship's
  shape.

### 3.3 The GP inside the simulator — already true

The note's goal 3 ("treat physics + GP residuals as the generative
model") needs no new mechanism: `simulate(observe=True)` on a problem
whose likelihood carries `GaussianProcessNoise` draws correlated
residuals with hyperparameters from their priors, on every backend. The
three comparisons the note wants (A+GP against B+GP; physics+GP against
physics+white on one model; the bare model's goodness of fit) are three
`FittingProblem`s differing only in their noise component — and that is
also the *noise-model selection study* the horizon memo planned for
Phase 10 (§6 (5), "which kernel?", depending on W8.6). ENFR gives that
study an amortised route and nested sampling its ground truth; they are
one item (W9.D). The prior-sensitivity sweep on the amplitude prior is a
retraining per setting and must be budgeted as such.

### 3.4 `reduced_chi2`, and "several statistics, never one"

The cheapest and least contentious proposal: the paper's main calibration
statistic added to the `calibration` group beside `ks_pvalue`,
`c2st_ranks` and `tarp_coverage`, with the docs' calibration section
saying why one statistic misleads (TARP hard to read when the posterior
is far narrower than the prior; C2ST near 0.5 for every method on the
lensing problem). Results-side, S, no dependency; a Phase 8 filler if
Peter wants it now (D6), else the first row of W9.B.

### 3.5 Cross-fitted ENFR — a hypothesis, kept as a switch

The argument (the reference flow is inflated at its own training points,
so the classifier learns its correction where the flow is optimistic) is
sound by the CF-NPE paper's logic, but neither paper tests it, and the
price is each classifier seeing N/K real simulations. Agree with the
note: a `cross_fit=` switch whose default is decided by the in-fold
against cross-fitted ablation on the Gaussian mixture with known log Z,
run once at W9.A (D7). The note's default of `True` is not adopted until
then.

### 3.6 What the note gets wrong or leaves out

1. The "size of change" row (CF-NPE "one more method"; ENFR "new module
   and result attributes") — see §3.1, §3.2.
2. It ignores W8.6: ampere's model comparison is ruled and will exist
   before ENFR could; ENFR joins it.
3. It plans a hand-rolled flow; `zuko` and `sbi.FMPE` are already locked.
4. It leaves the ENFR ensemble's cache key as the posterior's; an
   amortised artefact keys on the bank (§2).
5. It does not say where the compute comes from (§3.2, last bullet).
6. It does not connect the goodness-of-fit p-value to the deferred
   `ampere.diagnostics` (family A, pre-fit screening; RHMF deferred at
   W5.16 with a revisit trigger checked at each phase's close). An
   amortised prior-predictive evidence test *is* a pre-fit screen, and a
   better-founded one than RHMF for a problem that already has a bank; it
   is a candidate for that namespace's first content (D8).
7. Its science hook is stronger than it says. The Phase 9 flagship (SAGE
   with GRAMS, NESS data) is a survey of evolved stars whose first
   question is the dust chemistry — oxygen-rich against carbon-rich
   grids — which is *amortised model comparison across a survey*, ENFR's
   exact use case; and its parameter count under a population model is
   CF-NPE's regime. That is the argument for Phase 9.

## 4. Where it sits in the plan

### 4.1 Phase 8: two hooks and one optional filler, nothing else

Phase 8 is ruled, its waves are balanced two-wide, and its three §4
additions already fill the contract budget of a phase. Adding two L
items would push `1.0.0b3` out and would be measured against customers
(W8.9's population, W8.2's latent profile, W8.6's table) that land in
waves 2–4. What *should* change now costs a sentence each and avoids
rework later:

- **W8.6 (`compare`)**: the table takes any archived run carrying the
  evidence triple, renders `ampere_evidence_method` beside the value, and
  the page's "when each estimator is trustworthy" section leaves a row
  for a learned, amortised estimator with a recorded validation rung.
  Settle the method vocabulary (`"nested_sampling"`, `"bridge"`,
  `"savage_dickey"`, `"harmonic"`, and the slot `"enfr"`) in that item.
- **W8.3 (the emulator memo, Fable's own work)**: the training-set-to-
  model contract it drafts should name the bank as a shareable artefact
  with three customers — the emulator, `SBIEngine`, and an evidence
  engine — so the simulation budget is paid once (the note's goal, and
  already the `--training-set` idiom).
- **Optional filler**: `reduced_chi2` in the calibration group (§3.4).

### 4.2 Phase 9: a workstream, memo first

Phase 9's drafted shape (horizon §9) is nested populations, hierarchical
SBI, the columnar store and population inference at scale, the IFU cube,
the SAGE/GRAMS amortised-SBI flagship, with readers and `correct` as
fillers. The proposals add a fourth strand, **amortised evidence and
honest amortised posteriors**, that serves the same flagship and the same
personas (C and E, the SBI users; A and B through `compare`'s table).
Order within the strand: W9.A (memo + reproduction) → W9.B (CF-NPE) ∥
W9.C (the evidence engine) → W9.D (the comparison study, merged with the
Phase 10 noise-model study) → W9.E (the goodness-of-fit screen as
diagnostics family A, if D8 says so). W9.A can run between Phase 8's
waves as Fable's own work, as W8.3 does, so that Phase 9's drafting has
measurements rather than a note.

### 4.3 What stays on the horizon

The note's further extensions — a set or transformer evidence network on
the encoding tokens, a flow-matching residual model for CF-NPE, masked
mean regressors, cross-fitting for NLE, Simformer as a multi-model
reference, a cross-fitted Simformer — are correctly marked future or
speculative. They join horizon §6 as a bullet, with one remark: the
cross-fitted NLE idea (a learned mean of x given θ plus a learned
residual density) is ampere's own flexible likelihood with the GP
replaced by a density learned across simulations, which makes it a
natural M2-style comparison rather than a product feature.

## 5. The validation ladder as acceptance criteria

The note's gates translate to the house style — pinned rows at the
per-PR budget, measured rows under `heavy` or on the cluster, and
"informational" where the number is a finding rather than a bar:

| Note's gate | Item | Pinned or measured |
|---|---|---|
| l-POP recovers a known ratio; the flow's density integrates to 1; the standardiser's Jacobian | W9.C | Pinned, small |
| Gaussian mixture with known log Z, within ensemble scatter; the cross-fit ablation | W9.A then W9.C | Measured at W9.A (the default decided); pinned at W9.C at a reduced budget |
| BAO five-model reproduction (RMSE ≤ 0.46–0.81 nats); DESI DR2 | W9.A | Measured once, cluster, informational — it is the reimplementation's check against a paper with no code release |
| Nested-sampling agreement on ampere problems: RMSE < 0.5, bias < 0.2 nats on ≥ 50 held-out simulations | W9.C | Measured on the M2 toy at `heavy`; the bar is Peter's to set (the note says so) |
| Gaussian–Gaussian 200-parameter reproduction (test LL near −128.4; χ² 0.99–1.01); in-fold ablation raises χ² to about 1.15 | W9.B | Pinned at a reduced size on every PR (the folds are wired); the full numbers under `heavy` |
| CF-NPE against NPE on a many-parameter ampere problem at matched budget | W9.B | Measured on W8.9 stage 2 or the flagship; pinned as a margin if cheap (W4.5's precedent) |
| M2 misspecified spectra: the bare model flagged; the evidence prefers physics + GP | W9.D | Measured; pinned where the sign is unambiguous |
| ≥ 300 calibration samples for any NPE-against-CF-NPE comparison | W9.B | A rule in the item text, from `sbi.rst`'s finding |
| Reproducibility: bitwise-identical runs from a seeded problem | W9.B, W9.C | Pinned (in-process identity, never draws' bytes across machines — the W7.12 lesson) |

## 6. Decisions for Peter

- **D1 — placement.** Phase 9 strand with W9.A between Phase 8's waves
  (recommended); or Phase 8 fillers; or a separate experimental track
  outside the phases. *Recommended: Phase 9, W9.A early.*
- **D2 — the home and the name.** An evidence engine in
  `ampere.inference` writing the triple, with `ampere.results.compare`
  (W8.6) the consumer (recommended); or the note's `ampere.comparison`
  namespace. *Recommended: `ampere.inference`; a front door later if
  wanted.*
- **D3 — memo first, with a runnable reproduction.** W9.A reproduces the
  ENFR Gaussian-mixture benchmark and the CF-NPE Gaussian–Gaussian
  numbers in a prototype outside the package before any item is dispatched
  (recommended, because ENFR has no code release and the paper's
  experimental code was machine-written); where the prototype lives
  (`examples/evidence/` kept, or scratch and only the numbers kept).
  *Recommended: yes; `examples/evidence/` so the measurement is
  reproducible.*
- **D4 — the reference density.** `zuko` (`CNF` by flow matching against
  `NSF` by maximum likelihood, measured) rather than the note's
  hand-rolled flow. *Recommended: `zuko`, measured at W9.A.*
- **D5 — CF-NPE's acceptance problem.** The 133-disc population (W8.9
  stage 2) at matched budget (recommended), or a synthetic many-parameter
  problem only. *Recommended: W8.9's population, with the synthetic
  reproduction pinned.*
- **D6 — the Phase 8 hooks.** Amend W8.6's and W8.3's texts now (one
  sentence each); `reduced_chi2` as a Phase 8 filler (S; Sonnet) or
  W9.B's first row. *Recommended: amend both; the filler if a slot
  opens, else W9.B.*
- **D7 — cross-fitted ENFR's default.** Decided by W9.A's ablation, not
  adopted from the note. *Recommended: so.*
- **D8 — the goodness-of-fit screen's home.** The deferred
  `ampere.diagnostics` (family A) with ENFR's p-value as its first
  content, re-opening the RHMF question with a better-founded candidate;
  or a method on the engine only. *Recommended: decide at W9.A with the
  numbers; the engine method either way.*
- **D9 — compute.** The BAO/DESI reproduction and rung 1 at 50 nested
  runs per model are cluster runs on Peter's word (W8.7's pattern);
  everything pinned per PR is CPU-small. *Recommended: so, budgeted in
  W9.A's and W9.C's text.*

## 7. Draft items (for Phase 9's drafting; lifted as ruled)

### W9.A — The amortised-evidence memo and reproduction [M; Fable's own work between Phase 8's waves]
`docs/design/amortised_evidence_memo.md` (this memo) extended with
measurements, on W8.3's pattern: a prototype in `examples/evidence/`
(torch, the `sbi` environment) that (a) trains a `zuko` reference —
`CNF` by flow matching and `NSF` by maximum likelihood — and the l-POP
evidence network on the ENFR paper's Gaussian-mixture benchmark at a
budget the cluster affords, reporting RMSE in log Z against the analytic
value for flow alone, ENFR in-fold and ENFR cross-fitted (the §3.5
ablation; decides D7 and D4); (b) runs the CF-NPE paper's
Gaussian–Gaussian problem at 200 parameters and 10,000 simulations
through a stand-alone CF-NPE, reproducing the test log-likelihood and
reduced χ² within tolerance, and the in-fold ablation's χ² rise; (c)
times a `QuasisepGP` residual draw against the likelihood at 2,000 and
20,000 points (the note's scaling question); (d) drafts W9.B–W9.D's final
text and the D2 namespace. **Depends:** nothing (the `sbi` environment).
**Accept:** the memo with the measured tables and the decision list; the
prototype runs from the repository; no package code.

### W9.B — `SBIEngine(method="cf-npe")` [M–L; Opus] (after W9.A)
`ampere/inference/_cfnpe.py`: stage 1's K·M mean regressors on the
encoded data (shared frozen feature map by default, a per-stage encoder
optional), `torch.func.vmap`-batched when small; out-of-fold residuals
standardised per coordinate; stage 2 a diagonal Gaussian or `zuko.NSF`
over residuals conditioned on the features and the stage-1 mean; the
mixture variant; folds from the problem's sub-stream; the
`CFNPEPosterior` satisfying `calibrate()`'s surface; `cf_folds`,
`cf_replicates`, `cf_residual`, `cf_mixture` with the note's defaults;
refusals by name (folds < 2, truncation options, a context with the flat
layout); `ampere_sbi_cf_*` provenance in the cache key; `reduced_chi2` in
the calibration group (unless taken as a Phase 8 filler) and `sbi.rst`'s
calibration section rewritten to report several statistics; the
hyperparameter rule "tune one member, reuse for all" in the docs.
**Depends:** W9.A; W8.9 for the acceptance problem. **Accept:** the
Gaussian–Gaussian reproduction pinned at reduced size with the in-fold
ablation's χ² rise; CF-NPE's test log-likelihood above NPE's at matched
budget on W8.9's population (measured; pinned as a margin if under the
per-PR budget); bitwise-identical seeded runs; gates `sbi` + dev.

### W9.C — The evidence engine: ENFR and MultiENFR [L; Opus] (after W9.A; D2)
`ampere.inference.EvidenceEngine` (name per D2) over named problems
sharing a data hash: the banks through `simulate_many` with the
`training_set=` sink shared with `SBIEngine`; the `zuko` reference per
D4; the (N+1)-class l-POP network; ensembles with the `cross_fit=`
switch at D7's default; `evaluate` → `DataTree` with the evidence triple
(`ampere_evidence_method = "enfr"` / `"multi-enfr"`), `ampere_enfr_*`
provenance, the validation rung reached; `bayes_factor` with the
"not distinguished" verdict; `gof_pvalue` with its floor;
`validate(engine=)` over `replace_observations` and a nested engine; the
never-silent rules (`ENFRExtrapolationWarning`, the summary-space
attribute, timeouts over 1 % of a bank warned by name, cache misses
naming the field); the artefact keyed on the bank. `compare` (W8.6)
renders the row. **Depends:** W9.A, W8.6. **Accept:** the small pinned
rows (l-POP ratio, density normalisation, Jacobian, disjoint folds,
seeded identity); the Gaussian-mixture evidence within the ensemble
scatter at a reduced budget; rung 1 on the M2 toy under `heavy` at the
bar Peter sets; the cluster reproduction of the BAO table recorded as a
measurement; gates `sbi` + dev.

### W9.D — Model comparison with the GP inside the simulator: the noise-model selection study [M; Opus] (after W9.C and W8.6; absorbs horizon §6 (5))
On the M2 scenarios and the debris-disc flagship: A+GP against B+GP,
physics+GP against physics+white, and the bare model's goodness of fit,
by ENFR amortised and by nested sampling as ground truth; the kernel
library compared by evidence and LOO ("which kernel?"); the amplitude-
prior sensitivity sweep; a docs page of guidance with the measured table
and the "when not to compare by evidence" caveat carried from
`advanced.rst`. **Depends:** W9.C, W8.6. **Accept:** the measured table;
the bare model flagged and physics+GP preferred on the misspecified
scenarios, pinned where the sign is unambiguous; gates `sbi` + nested.

### W9.E — The goodness-of-fit screen as diagnostics family A [S–M; Sonnet] (if D8 so rules)
`ampere.diagnostics` opened with the prior-predictive evidence p-value as
its first screen, over a trained evidence engine, beside the posterior-
side `chi_square_pvalue()` and the GP absorbed-power diagnostics; the
RHMF deferral re-ruled in the same decision-log row. **Depends:** W9.C.
**Accept:** the contract's family A row amended with its decision-log
entry; the screen on the M2 scenarios; gates dev + `sbi`.
