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

*Addendum, later on 2026-10-21.* §8 reads the proposals against the
reserved horizons and the Phase 9–10 outlines; §9 reads the state-of-
the-art survey Peter asked for (`sbi_state_of_the_art_2026.md`) against
both. The survey changes one thing in the short answer above: because
every ampere problem has an evaluable likelihood, the evidence route the
literature favours — importance sampling under the amortised posterior,
horizon (f)'s `correct` module — lands *first*, and ENFR is trialled
against it (§9.1, D10). CF-NPE stands (§9.3). The decisions are D1–D13
(§6 and §9.6); the draft items are W9.0–W9.D (§7 and §9.7).

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
  **Ruled 2026-10-22 (Peter): Phase 9 as recommended — the strand opens with the importance-sampling evidence (W9.0), W9.A early as Fable's own work between Phase 8's waves.**
- **D2 — the home and the name.** An evidence engine in
  `ampere.inference` writing the triple, with `ampere.results.compare`
  (W8.6) the consumer (recommended); or the note's `ampere.comparison`
  namespace. *Recommended: `ampere.inference`; a front door later if
  wanted.*
  **Ruled 2026-10-22 (Peter): `ampere.inference`; other options (a front door, a namespace) later if required or desired.**
- **D3 — memo first, with a runnable reproduction.** W9.A reproduces the
  ENFR Gaussian-mixture benchmark and the CF-NPE Gaussian–Gaussian
  numbers in a prototype outside the package before any item is dispatched
  (recommended, because ENFR has no code release and the paper's
  experimental code was machine-written); where the prototype lives
  (`examples/evidence/` kept, or scratch and only the numbers kept).
  *Recommended: yes; `examples/evidence/` so the measurement is
  reproducible.*
  **Ruled 2026-10-22 (Peter): yes — memo first, the reproduction kept in `examples/evidence/`, with §9.6's ENFR-against-importance-sampling comparison.**
- **D4 — the reference density.** `zuko` (`CNF` by flow matching against
  `NSF` by maximum likelihood, measured) rather than the note's
  hand-rolled flow. *Recommended: `zuko`, measured at W9.A.*
  **Ruled 2026-10-22 (Peter): `zuko`, CNF against NSF measured at W9.A.**
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

## 8. The fit with the design horizons and with the Phase 9–10 outlines (added 2026-10-21)

Peter asked for the elaboration before ruling: how the two proposals sit
against the horizons the plan reserves (`DEVELOPMENT_PLAN.md` §5, "Design
horizon", (a)–(i)), against the inference-extensions memo's rulings of
2026-09-10, the horizon notes of 2026-09-09, the user-journeys ranking,
and the Phase 9 and 10 shapes in `horizon_beyond_phase7.md` §9.

### 8.1 Against the nine reserved horizons

| Horizon | What it reserved | What ENFR / CF-NPE do with it |
|---|---|---|
| (a) Multi-fidelity | Fidelity tags on channels and registrations; "the training-set format carries fidelity per sample (a column, not an attr)" (inference memo §9) | A shared simulation bank mixing an emulator's draws with the simulator's is the first place that column is read: an ENFR member trained on the cheap fidelity and corrected on the expensive one is the method's own two-step shape (flow then classifier) applied across fidelities. Nothing to reserve beyond what §9 already asked; W8.3's bank contract should carry the column. |
| (b) Population by post-processing, at scale | Per-draw `log_prior`/`log_likelihood` stored; the columnar store; "NRE is the natural amortised estimator for population inference" because the per-object term is `p(xᵢ|Λ) = ∫ p(xᵢ|θ) π(θ|Λ) dθ` and a ratio needs no interim prior (§8.1) | Two things. CF-NPE's stage-2 density is a proper `log_prob`, so a CF-NPE run stores its own proposal density per draw and the reweighting route works as for NPE — with honest widths, which is what keeps the effective sample size from collapsing when the population prior narrows. More interesting: an evidence network **conditioned on the hyperparameters Λ** estimates `p(xᵢ|Λ)` directly — the per-object term the population objective sums — amortised over objects *and* over Λ, so the SVI over minibatches of sources that §8.1 wanted runs on a stored network with no draws at all. That is the ratio-estimator argument of §8.1 with ENFR's absolute normalisation; it is a horizon idea for Phase 9's hierarchical SBI, not an item. One question for W9.A: does an NRE classifier trained for many epochs on a fixed bank show CF-NPE's defect (its loss is also unbounded for a flexible network)? If so, cross-fitting belongs to the ratio route too. |
| (c) Emulation | Serialisable training sets; the emulator as a `Model`; its predictive variance as a noise component; the emulator's identity in every dependent cache key | The emulator is what makes ENFR's budget (10⁵ simulations per model in the paper) reachable for a radiative-transfer model: the pipeline is simulator → `train_emulator` (W8.4) → cheap emulated simulator with its error component → bank → evidence. Because the error component is part of the simulator, the evidence is honestly that of "emulated physics + emulator error + GP". The reservation to add: an ENFR ensemble's cache key includes the emulator's training-set identity, exactly as §9 (c) asks for any dependent artefact. |
| (d) Hierarchical SBI | Plate-aware groups; a `simulate()` drawing hierarchically; a set-of-sets embedding; "`EncodingLayout` is not assumed to have a fixed leading dimension" | CF-NPE's regime (tens to hundreds of parameters) is the population-level amortised posterior — the hyperparameters plus every member — not the per-object fit. Hierarchical SBI is therefore CF-NPE's Phase 9 customer, with W8.2's `LatentProfile` and images the customers it has before then. The encoding reservation stands unchanged. |
| (e) Model comparison and averaging | The evidence triple; consumers "nested sampling, SMC and Laplace, `harmonic`, bridge, Savage–Dickey, LOO/WAIC, BMA" | The list has no simulator-only producer. ENFR is it, and the first producer whose `ampere_evidence_method` is a learned estimator with a validation rung — the vocabulary question `inference.md:1156` left open. Add one line to the horizon: "an amortised evidence estimator for problems without a likelihood, consumed by the same table". |
| (f) Approximate-inference correction | The proposal's own log-density per draw; `ampere_approximation` names the family; PSIS reweighting with `k̂` stored | A CF-NPE run is a correctable proposal like any NPE run; its `k̂` under `correct` (Phase 9 filler) is a per-observation overconfidence flag that complements the ensemble scatter and the extrapolation flag. ENFR's reference flow `log q(x)` is itself a proposal density over *data*, which is what its extrapolation flag reads; nothing to reserve. |
| (g) Engines as a registry, cost accounting | A registry; one cost record per engine, per fidelity | An evidence engine's cost record is simulations per model, failures by kind, trainings per member — the fields the note's `ampere_enfr_failure_fraction` and `_budget` already name; they should be the shared cost record, not a private family. The engine battery (inference memo §4: "SBC for every engine, for free") gains a row that rank-calibrates model probabilities on held-out simulations — the ladder's rung 2. |
| (h) Vector observables with correlated misspecification | The joint noise model binds to the channel pairing | Untouched; a comparison that includes the joint GP draws from it like any other (§3.3). |
| (i) Amortisation over observation context | `context=` on `simulate_many`/`SBIEngine`; the σ column mandatory | ENFR trained under a `ContextPrior` is an evidence valid across a survey's noise levels, which is the whole case for amortising it; the engine takes the same `context=`. The horizon notes' further step — amortising over instruments (§5: the grid, the resolution, the filter set as per-draw inputs) — remains the open item it was, and applies to both methods equally. |

### 8.2 Against the inference-extensions memo's rulings

- **Every new sampler behind its own extra** (ruling 7.4). Both methods
  live behind `sbi`, which already resolves `zuko`; no new extra.
- **The stored proposal density is the rule** (ruling 7.1). CF-NPE
  satisfies it by construction; ENFR has no draws (§8.7).
- **BO is acquisition, not an engine** (ruling 7.3). The ENFR paper's
  finding that the classifier correction matters most at small budgets
  makes acquisition for the *bank* — where to simulate next so the
  evidence's error falls fastest — the same question W8.3's memo already
  owns for the emulator; one answer for both.
- **§8.4's `compare` as the one-line table.** ENFR's `evaluate` should
  return W8.6's typed comparison result, so that a nested-sampling
  comparison and an amortised one render identically and the method
  column tells them apart.

### 8.3 Against the horizon notes (embeddings, context)

The horizon notes' follow-up on embeddings (§4/5) ranked DeepSets, set
transformers, neural-process encoders (ConvCNP) and neural-operator
branch nets as embeddings for the W3.3 encoding, and recorded the trap
that makes the packing the contract. Both proposals consume that packing
unchanged. What they add to the embedding question is one requirement:
an embedding trained for *parameter* inference may discard what
*discriminates models*, so the evidence engine trains its own (or trains
a shared one jointly with the classifier), and the ConvCNP experiment the
notes recommended is as relevant to the evidence network as to the
posterior. The context programme (§4) is the mechanism that makes an
amortised evidence worth having at all.

### 8.4 Against the user-journeys ranking

By persona weight (user-journeys §6) the model-comparison module serves
A, B and C (60 %, rank 4) and the SBI-side items serve C and E (40 %,
rank 8–14). ENFR inherits the first weight *through* `compare` and the
second directly; CF-NPE serves C and E only. Neither is a top-five item
by that table, and the memo should not pretend otherwise. What lifts
them is the flagship and the papers: the Phase 9 flagship's first
question (oxygen-rich against carbon-rich dust across SAGE) is amortised
model comparison across a survey, the catalogue loop's listed output "the
evidence contrast grid versus blackbody per object, and the objects the
grid cannot explain" (user-journeys §3) is exactly `evaluate` plus
`gof_pvalue` at survey scale where nested sampling per object is too
expensive, and D12's component-led science papers on NESS data are
papers about population models and model comparison. The honest ranking:
below the convergence, units and reader items; level with the correction
module; above the IFU cube until a cube user appears.

### 8.5 The Phase 9 outline, re-read with the strand

Phase 9 as shaped (horizon §9; Q3 ruled "Phase 9 by plan") is: the
nested-populations memo and item (L, first, alone on the parameter
layer), hierarchical SBI (M, after it), the columnar store and population
inference at scale (L, with a store memo), the IFU cube (L), with the
light-curve and Gaia readers, `correct` (M) and polarimetry (S–M) as
fillers; SAGE with GRAMS as the amortised-SBI flagship; the catalogue
loop as an example before the store. Where each proposal attaches:

- **W9.A (the memo and reproduction)** has no dependency and is Fable's
  own work; it can run between Phase 8's waves, as W8.3 does, so that
  Phase 9's drafting has numbers. It should also answer the NRE question
  of §8.1 (b), because the store memo's choice of per-object columns
  (draws against a stored ratio or evidence network) depends on it.
- **W9.B (`cf-npe`)** has two customers in order: W8.2's `LatentProfile`
  under SBI (the misspecification-aware simulator the radial-profiles
  design named; a latent GP's coefficients are a many-parameter vector)
  the moment Phase 8 closes, and hierarchical SBI's population-level
  posterior once that lands. Draft it after the nested-populations item,
  in the slot beside hierarchical SBI, with the `LatentProfile` row as
  its first acceptance problem and the population row added when
  hierarchical SBI merges.
- **W9.C (the evidence engine)** serves the flagship and the catalogue
  loop. It needs W8.6's result type and the emulated GRAMS-style grid (the
  flagship's own model work) — so it sits in the flagship's wave, with the
  chemistry comparison and the per-object goodness-of-fit screen as its
  acceptance rows beside the Gaussian-mixture and M2 rows of §7.
- **`correct` (filler)** pairs with W9.B: the `k̂` flag on a CF-NPE run is
  its per-observation check, and the correction module is the route from
  an honest amortised posterior to the population reweighting the store
  serves. Keep them adjacent.
- **The store memo** should take the §8.1 (b) conditional-evidence idea
  as one of the per-object representations it weighs, beside draws and
  ratios.
- **W9.D (the comparison study)** closes the strand after W9.C and W8.6;
  it is the measured table the docs' "which kernel?" guidance needs and
  the first paper-shaped result of the strand.

Nothing in the strand touches the parameter layer, so none of it collides
with the nested-populations item; all of it is torch-side under the `sbi`
extra, so it runs beside jax-side Phase 9 work without sharing files.

### 8.6 The Phase 10 outline, re-read

Phase 10 as shaped is pocoMC and jaxns (S each, behind extras), the
noise-model selection study (M), the declarative problem file (S–M), the
cookbook and gallery (M), conda-forge, and `1.0.0` with its criteria.
Three consequences:

- **The noise-model selection study is absorbed** by W9.D and moves
  earlier, with the amortised route beside the nested one; Phase 10 loses
  an M item.
- **pocoMC and jaxns become evidence producers** for the same table (SMC
  gives an evidence as a by-product; jaxns runs on the jax problem on the
  GPU), which makes them the natural ground truth for the ladder's rung 1
  on problems where dynesty is slow — a reason to pull them into Phase 9's
  fillers rather than leave them for Phase 10, if slots allow.
- **`1.0.0`'s criteria** (horizon §8: the known-limitations list reduced to
  what is deliberate, conda-forge, the paper's examples pinned, v2's
  deprecations removed, the Zenodo record) are unaffected provided the
  evidence engine ships *experimental* behind the extra, with the
  validation rung it reached in every result and the docs saying what
  "experimental" means — a learned evidence that has not passed rung 1 on
  the user's own problem is a number with a caveat, and the known-
  limitations list says so. CF-NPE needs no such caveat: it is an
  ordinary amortised posterior that `calibrate()` checks.

### 8.7 One output-shape question the note missed

`Engine.finish` assembles equal-weight draws into a run, and every run in
`results.md` has a `posterior` group. An evidence engine produces no
draws: its output per observation is a table of evidences with errors, a
goodness-of-fit p-value and the provenance of the ensemble. The
optimisers met the same problem and were given the `optimum` group by
ruling (inference memo §7.2, W6.7). The two routes: (i) **reuse W8.6's
typed comparison result** — the Bayes-factor table that bridge sampling
and Savage–Dickey already return — so that an amortised comparison and a
nested-sampling one are the same object with a different method column
(recommended; this is the Phase 8 hook of §4.1 made precise: W8.6
defines the result type with a future producer in mind); or (ii) a new
group with its own `results.md` decision-log row at W9.C. Route (i)
keeps the "no §4 change" claim of §0 true; route (ii) makes it false by
one row. Either way the per-model evidences must be archivable with the
ensemble's provenance, so that `compare` can read them later without the
networks — the same "reconstructible from the archive alone" rule §9 (b)
set for the interim prior.

## 9. What the 2026 survey changes (added 2026-10-21, later the same day)

`sbi_state_of_the_art_2026.md` is the survey Peter asked for: robust SBI,
model comparison and evidence, SBI architectures and software, the
likelihood side, and the astrophysical software landscape, January 2025
to October 2026, with a candidate-tools table at its end. This section
reads it against §§3–8 and against the Phase 9–10 outlines. Four
findings change the framing above; the rest add items to the plan's
horizons or to items already ruled.

### 9.1 ampere always has a likelihood, and that reorders the evidence routes

The note, and §3.2 above, framed ENFR as the simulator-only route to an
evidence. The survey's one controlled study of amortised model
comparison (Kucharský et al. 2025–26) finds that **posterior-based
evidence routes beat direct evidence networks** in every case, that the
direct networks stay poorly calibrated out of distribution with or
without self-consistency training, and that self-consistency repairs the
posterior-based routes only when a likelihood is evaluable. ENFR was not
in that comparison and no paper in the window reports nats-level
accuracy for any direct evidence network; the ENFR paper's own 0.5–0.8
nats on BAO is the only number there is.

ampere is in the position most SBI users are not: every problem has an
evaluable likelihood, the GP marginal over the model's output. What SBI
amortises in ampere is an *expensive model*, not an intractable
likelihood. So the evidence route the survey favours is available here:
the amortised posterior as proposal, the GP marginal likelihood times
the prior as target, importance weights from the stored per-draw
proposal densities — which is **horizon (f)'s `correct` module made an
evidence producer**, with the importance-sampling efficiency per object
as the certificate the astronomy landscape now expects (DINGO's
`lvk-stable` line, BilbyFlow's 29 of 38 events above 1 % efficiency,
Gebhard et al.'s retrievals). The hooks were decided at W5.0; the stored
proposal density exists for every SBI run today; the module is a Phase 9
filler already. Its cost per object is the likelihood at the draws, which
for a radiative-transfer model without an emulator is the only thing
that makes it expensive — and W8.4's emulator removes that.

The consequence for ENFR is a narrower, more honest niche: **the survey
scale at which even a few thousand likelihood evaluations per object are
too many, and the genuinely likelihood-free case, which in ampere is
rare.** An ENFR trial must therefore be scored against the
importance-sampling route on the same bank, not only against nested
sampling; if it does not win in nats where the posterior route is
affordable, its place is the survey case alone. That is a cheaper and
better-founded gate than the note's ladder, and W9.A should run it.

### 9.2 Ground truth is now cheap on the jax backend

GPU-native nested sampling (blackjax nested slice sampling, TMLR 2026;
jaxns 3.0, both Apache-2.0) gives evidences "statistically identical" to
the CPU implementations at 20–40× speed. The ladder's rung 1 (50 nested
runs per model) stops being cluster work for any problem on the jax
backend, which is also where the flagship's emulated model will live.
Phase 10's jaxns item moves to Phase 9 as the strand's ground truth; the
unconfirmed point is whether nested slice sampling is upstream in
`blackjax` or in a fork, which W9.A checks. pocoMC, by contrast, has
had no release since September 2024 and should be struck from Phase 10.

### 9.3 CF-NPE stands, and has two cheaper neighbours to trial first

No cross-fitting work beyond Regier's and no calibration-aware training
loss since 2023 appeared, so CF-NPE has no published competitor on its
own ground. Two cheaper things sit beside it: CP4SBI's post-hoc local
conformal credible sets from the held-out simulations `calibrate()`
already draws (MIT; better conditional coverage on 8 of 10 benchmarks
for NPE flows), and the five-seed ensemble the practical guide found
"well calibrated" where each member was slightly overconfident, whose
pairwise-KL disagreement (Alvey et al.; MIT) is also a per-observation
out-of-distribution flag. Both are S items on the artefact cache with no
new mathematics, and both are honest comparisons for W9.B's acceptance
row rather than substitutes for it: conformal sets calibrate the
amortiser, not the model, and ensembles cost K trainings for a width that
is not learned. Two test-design findings bear on W9.B's rows: coverage on
average is not coverage per realisation (Alokda et al. 2026: variance and
kurtosis disagree with likelihood-based inference even when coverage
passes), and a population pipeline needed a post-hoc coverage
calibration before combining per-object ratios (Jarugula et al. 2026).
W9.B's validation should therefore include a fixed-parameter
per-realisation check, which the M2 harness can carry.

The self-consistency loss (Mishra et al.) is CF-NPE's natural companion
— it, too, exploits the evaluable likelihood — and it needs a
density-evaluable head, which pulls against the survey's other
recommendation that flow matching become the default head. The engine
should say which head is in use and why: spline flows where
self-consistency is on, FMPE (or VP-EDM NPSE above about fifty
parameters) elsewhere, all in the locked `sbi 0.27`.

### 9.4 The GP inside the simulator is unoccupied, and RVNP is the warning

No published work puts a GP residual model inside the simulator for
neural SBI (an absence from searches, not a review's statement). The
nearest neighbour is RVNP on 1,053 Gaia XP spectra, where a free
θ-dependent error covariance collapsed and lost the metallicity; the
survey reads that as independent evidence that the residual model needs
a kernel-structured prior with few hyperparameters, which is ampere's
construction. Rotman et al. (ApJ 2025) give the thesis its first
controlled number on real instruments: abundance offsets of 1.2–4.6 dex
fall to 0.2–0.4 dex with a global-plus-local GP, the misfit localised to
one bin, and the GP *not* preferred when no correlated noise was
injected. The comparison study W9.D proposes is, by that reading, a
paper: the kernel-structured residual against a free covariance on the
same XP-style task, with the bank shared with the evidence work. This is
the strand's strongest scientific claim and belongs in D12's list of
component-led papers.

### 9.5 What else the survey adds to the plan (by phase)

Items that are ruled or in flight, where one sentence now saves rework:

- **W8.5 (photometry and the filter library)**: `speclite 1.0` is a live
  filter library and `sedpy` is dead since 2020; the item's SVO fetch
  stands, and speclite is the fallback bundle candidate D11 wants struck at
  `1.0.0`.
- **W8.16 (`FixedCovarianceNoise`)**: when the covariance is *estimated*
  (an azimuthal average's bins, an emulator's residual covariance), the
  Sellentin–Heavens correction and Heavens, Whiteway and Sellentin's
  (2026) Student-t with extra Gaussian error sources are the right family;
  a `family="student_t"` route on the same component is the hook, and it
  is the concrete recipe for W8.3's "emulator error as a noise component".
- **W8.3 (the emulator memo)**: Różański and Ting's scaling law (10× compute
  for ~7× lower MSE) and the asteroseismology rule "residuals below the
  survey's uncertainty or model them" are the budgeting rules; multi-
  fidelity NPE (ICLR 2026, up to 100× fewer high-fidelity simulations) is
  horizon (a)'s consumer and wants the fidelity column in the bank.
- **W8.6 (`compare`)**: `harmonic 1.3.1` now carries flow-based
  Savage–Dickey as well as the learnt harmonic mean, so one dependency
  covers two of its five estimators — subject to reading its second
  licence file; the importance-sampling evidence of §9.1 is a sixth
  estimator the table should name from the start.
- **W5.16 (the RHMF trial, deferred)**: Hilder, Hogg, Casey and Rix's
  robust heteroskedastic matrix factorisation is accepted in ApJ (2026)
  with a JAX implementation applied to Gaia RVS; the deferral's revisit
  trigger has fired, and D8's "family A's first content" question now has
  two candidates.
- **W7.10's jax scan provider**: `smolgp` (AJ 2026; associative-scan
  state-space GP, exposure integration, quasiperiodic kernels,
  tinygp-compatible) is its direct successor and the benchmark to run at
  W8.7's solver table.
- **The docs**: the `sbi` practical guide's validation chain (prior-
  predictive check, pre-training misspecification test, posterior
  predictive checks, global then local coverage, "necessary but not
  sufficient") is the checklist the SBI page should state; `sbi 0.27`'s
  two misspecification tests are wired in a filler, with the documented
  caveat that one high-dimensional observation has weak power.
- **Kernels and M2**: Rotman et al.'s localised Gaussian-bump kernels
  with free centres are a kernel-sum recipe the horseshoe already selects
  over, and their MSE (bias² + variance) metric is a reporting form the M2
  rows can add.

Phase 9 and 10 additions and removals:

- Phase 9 gains, as fillers: the importance-sampling evidence and
  per-object efficiency (S–M, on `correct`), CP4SBI and the ensemble flag
  (S each), the `sbi` misspecification tests (S), FMPE/NPSE as selectable
  heads with the head recorded in provenance (S), jaxns or blackjax nested
  slice sampling as the GPU evidence engine (S, from Phase 10), and an
  online fresh-batch training loop for emulated models (S–M).
- Phase 9's population work takes two findings: per-object draw budgets
  far above 2 × 10³ before reweighting (two 2026 papers quantify the
  Monte Carlo failure), and a learned population prior as a prior node
  (pop-cosmos' diffusion prior, with its admitted overconfidence as the
  calibration target). The store memo should weigh a stored evidence
  network (§8.1 (b)) beside draws and ratios.
- Phase 10 loses pocoMC and the noise-model study (absorbed by W9.D), and
  keeps jaxns only if Phase 9 did not take it.
- Not adopted: RoPE and minimum-distance summaries (need labelled pairs
  or many i.i.d. observations), direct evidence regressors, LtU-ILI
  (pins `sbi` ≤ 0.22), swyft, lampe, `nbi`, `sedpy`; BayesFlow 2 on jax
  stays the watch item for jax-native SBI, which answers the plan's §6
  deferral with "not yet, and not sbijax/flowjax".

### 9.6 The decisions, revised

D1–D9 of §6 stand with these changes:

- **D1** (placement): unchanged — Phase 9 strand, W9.A early; but the
  strand opens with the importance-sampling evidence on `correct`, not
  with ENFR.
- **D3** (memo first): W9.A's reproduction gains the comparison that
  decides ENFR's place — ENFR against the importance-sampling route on
  the same bank, scored in nats against GPU nested sampling.
- **D5** (CF-NPE's acceptance): add the per-realisation check and the
  CP4SBI and ensemble rows as comparisons.
- **D8** (the goodness-of-fit screen's home): two candidates for family
  A now — the amortised p-value and the landed RHMF.

Added:

- **D10 — the evidence route's order.** The importance-sampling evidence
  under the amortised posterior (horizon (f), `correct`) lands first and
  is the default amortised evidence; ENFR is trialled against it and
  takes the survey-scale and likelihood-free cases only if it wins.
  *Recommended: so.*
- **D11 — ground truth.** jaxns or blackjax nested slice sampling pulled
  from Phase 10 into Phase 9 as the jax-backend evidence engine and the
  ladder's rung-1 reference; pocoMC struck. *Recommended: so, the
  upstream question checked at W9.A.*
- **D12 — the paper.** The kernel-structured residual against a free
  covariance on an XP-style task (W9.D's design) as one of the
  component-led papers of the Phase 8 ruling's D12. *Recommended: yes,
  with W9.D drafted to produce it.*
- **D13 — the small hooks now.** The sentences in §9.5 for W8.5, W8.16,
  W8.3 and W8.6 added to their texts before dispatch. *Recommended:
  yes.*

### 9.7 The draft items, revised

- **W9.0 — Importance-sampling evidence and per-object efficiency** [S–M;
  Opus] (new; before W9.A): `ampere.results.correct` (horizon (f), the
  Phase 9 filler) returning, for any run with a stored proposal density,
  the Pareto-smoothed importance weights, `k̂`, the efficiency, and the
  evidence with its error, written as the triple with
  `ampere_evidence_method = "importance_sampling"`; the efficiency and
  `k̂` as the per-object certificate on every SBI run's summary; the
  M2 toy's evidence against dynesty as the pinned row. Depends on
  nothing. This is the item §9.1 makes first.
- **W9.A** gains: the ENFR-against-W9.0 comparison in nats on the
  Gaussian mixture and the M2 toy; the NRE-overfitting probe of §8.1 (b);
  the `blackjax`/jaxns check of D11; the head question of §9.3.
- **W9.B** gains: CP4SBI and the five-seed ensemble as comparison rows;
  the per-realisation variance and kurtosis check; the head recorded in
  provenance.
- **W9.C** is conditional on W9.A's comparison and narrows to the survey
  and likelihood-free cases if ENFR does not win where W9.0 is
  affordable.
- **W9.D** is drafted as the paper of D12.
