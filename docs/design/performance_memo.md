# The benchmark-driven optimisation pass — profile and levers

*W5.17, drafted 2026-09-27 by Opus. The plan's last Phase 5 bullet
(`DEVELOPMENT_PLAN.md` §5): "profile against the CI benchmark baselines
established in Phase 2, then attack the levers in evidence order". This memo
is the profile (§§1–4), the ranked lever list it produces (§5), the
proposals that would need a contract change (§6), and — once the levers have
been tried — what landed and what did not (§7). No lever lands without its
before/after row; nothing in `tests/benchmarks` asserts a time, and nothing
here changes that.*

## 0. The short answer

The lifecycle docstring's claim holds: **nothing in the hot loop
re-negotiates, re-compiles or re-merges at the problem level** (measured, §3.1).
Construction — negotiation, compilation, the merge and the one validating
evaluation — costs 6–16 ms per problem, once. The time goes elsewhere, and
where it goes depends on the load:

| Load (the study's own small budget) | Where the time goes, in order |
|---|---|
| M2, 200 points, flexible, reference | the six scalar prior densities (scipy frozen `logpdf`), ~45 %; the quasiseparable solve, ~17 %; building the model's `Spectrum`, ~15 % |
| M2 on jax, the numpy-facing contract path | **XLA recompiling celerite2's `lax.cond` on every call**, ~95 % (31 ms per call against 0.28 ms realised) |
| Interferometry, flexible arm, reference | `FourierSample.apply` re-deriving its DFT matrices and cell solid angles every call, ~33 %; the priors, ~25 %; the dense GP, ~13 % |
| Image, 24², flexible, reference | the dense kernel matrix, ~60 %; Cholesky factor and solve, ~15 % |
| Image, 64², flexible, reference | `cho_factor` ~40 %, **`cho_solve` ~31 % — a layout copy, not arithmetic**, the kernel matrix ~25 % |
| Any spectrum through `Resample` (issue #12) | the weight matrix, **re-planned every call**: 47 % of the step at 2 000 → 200 points, 92 % at 20 000 → 2 000 |

Five levers are on the plan's candidate list, are exactly-zero changes,
and landed (§5, §7): the jax contract path compiles once (up to 6.8×), the
dense solve drops a layout copy (2.6× at 64²), `Resample` and
`FourierSample` plan once per grid (30–65× on the step; up to 2.1× on the
interferometry `log_prob`), and the M2 models reuse their container
(1.2–1.35×). The biggest remaining item on the M2 reference load — the
scalar prior densities, 37 % — is measured with an 11× floor under it but
needs either a scipy private-API coupling or a non-bit-identical closed
form, so it is a proposal (§6, P1), with four others.

## 1. The machine and the loads

**Machine.** AMD Ryzen 9 6900HS (8 cores, 16 threads), 13 GB RAM, WSL2 on
Linux 6.6. Python 3.13.15, numpy 2.5.2, scipy 1.18.1 (scipy-openblas
0.3.34, 64-bit ints, dynamic arch), celerite2 0.3.3, torch 2.13.0+cpu,
jax/jaxlib 0.11.1. No accelerator (GPU placement is out of scope). A second
agent (W5.15) shared the machine throughout; every timing run went through
`flock /tmp/ampere-gate.lock`, so no two measurements overlapped, but the
machine was not otherwise idle, and the benchmark rows' inter-quartile
ranges (§2) show it. A lever therefore has to beat that scatter, not merely
move a median.

**Base commit.** `3320eef` (master at dispatch).

**The loads**, each at the study's own small budget:

* **M2** — `examples.m2_misspecification.study.build_problem("strong_smooth",
  size=200, likelihood="flexible")`: the Matérn-3/2 `QuasisepGP` at 200
  points, the problem `python -m examples.m2_misspecification --quick` runs
  first. On the reference backend, torch and jax.
* **Interferometry** — `examples.interferometry.study.build_problem(backend,
  "flexible")`: visibilities under the complex GP (`DenseGP`, 2-column
  right-hand side) plus closure phases under von Mises noise.
* **Image** — `examples.image.study.build_problem(backend, "flexible",
  pixels=...)` at `SMALL_PIXELS` (24², N = 576) and `--pixels 64` (N = 4 096),
  `DenseGP` over a two-axis Matérn-3/2.

**The loop.** 300 proposals drawn from the problem's own prior
(`problem.sample_prior`, seed 1), scored by `problem.log_prob` (the
contract path every gradient-free engine uses) or, for the modern
backends, by `realise(problem).log_prob_unconstrained` with its gradient
(`backward()` on torch; `jax.jit(jax.value_and_grad(...))` on jax — one
NUTS leapfrog step's work). One warm-up call first. Each loop is run twice:
once under `cProfile` (the tables, sorted by cumulative time, top 40), once
bare (the ms/call below). The harness is reproduced in Appendix A.

## 2. The baselines (pytest-benchmark, base commit)

`pixi run -e <env> bench` for dev, torch and jax, one at a time through the
lock; the JSON files are kept outside git at
`~/.cache/ampere-gates/w5.17/baseline-<env>.json`. Medians and IQR in ms;
only the rows a lever below touches, plus the headline rows. The full
tables are in the JSON.

| Row | dev median (IQR) | torch median (IQR) | jax median (IQR) |
|---|---|---|---|
| `test_reference_dense[200]` | 1.32 (0.35) | 1.40 (0.37) | 1.46 (0.52) |
| `test_reference_dense[2000]` | 308 (269) | 497 (66) | 248 (202) |
| `test_reference_dense` (GP marginal, n1000) | 47.1 (22.9) | 67.7 (37.2) | 48.7 (10.5) |
| `test_dense_anchor` (image `DenseGP`, n1024) | 208 (30) | 199 (19) | 63.7 (86.3) |
| `test_route_cholesky` (EFGP normal equations, m1089) | 326 (574) | 564 (75) | 79.9 (67.2) |
| `test_reference_quasisep[200]` | 0.722 (0.101) | 0.735 (0.108) | 0.717 (0.107) |
| `test_reference_quasisep[20000]` | 2.99 (0.74) | 2.49 (0.34) | 2.68 (0.53) |
| `test_reference_emcee[flexible]` | 1716 (14) | 1773 (26) | 1977 (421) |
| `test_backend_quasisep_contract[200-jax]` | — | — | 30.2 (3.3) |
| `test_backend_quasisep_contract[20000-jax]` | — | — | 31.7 (1.4) |
| `test_one_proposal[contract-QuasisepGP-jax]` | — | — | 28.5 (3.3) |
| `test_one_proposal[realised-QuasisepGP-jax]` | — | — | 0.671 (0.174) |
| `test_backend_quasisep_realised[200-jax]` | — | — | 0.286 (0.043) |
| `test_backend_nuts[flexible-torch/jax]` | — | 7962 (110) | 3167 (2231) |

Two readings. First, the scatter: several dense rows have an IQR of the
same order as their median, which is the shared machine, and a lever that
moves such a row by less than a factor of about two cannot be told from
noise on it. Second, and already a finding: **the jax contract row is flat
at ~30 ms from 200 to 20 000 points** while the realised row is 0.29 ms at
200 — a per-call constant that has nothing to do with N (§3.2).

## 3. Where the time goes, per load

Per-call times, bare loop (profiled in brackets), base commit:

| Load | reference contract | torch | jax |
|---|---|---|---|
| M2, n = 200 | 1.05 ms (1.23) | contract 1.19; realised v+g 3.50 | **contract 31.3**; realised v+g 0.283 |
| Interferometry, flexible | 1.19 ms (2.05) | realised v+g 8.10 | realised v+g 0.473 |
| Image 24², flexible | 18.2 ms | realised v+g 34.6 | contract 24.0; realised v+g 45.9 |
| Image 64², flexible | 1 789 ms | — | — |
| Construction (`FittingProblem`) | M2 5.8 ms; interferometry 16.2 ms; image 64² 10.2 ms | | |

### 3.1 The lifecycle claim

`ampere/core/dataset.py`'s docstring says negotiation, compilation, the merge
and alignment checking happen once, in `FittingProblem.__init__`. Counting
calls over 50 `log_prob` evaluations of each load: **zero** calls to
`negotiate`, `compile_for`, `ParameterSet.merge`, `check_alignment` or the
instrument's `_declarations` fingerprint. The claim holds. What *is* redone
per call is below the problem level, inside steps and solvers, and is what
§5's levers are about. Construction itself is dominated by building scipy
frozen distributions (M2) and by the one validating evaluation
(interferometry): 6–16 ms, never worth a lever.

### 3.2 M2

* **Reference, 1.05 ms.** `ParameterSet.lnprior` is 47 % of the profiled
  loop: six scalar parameters, each through a scipy *frozen* distribution's
  generic `logpdf` (`_distn_infrastructure.py:2093` — argument parsing,
  `argsreduce`, broadcasting, support masking), ~40 µs apiece for one
  number. The quasiseparable solve (`QuasisepGP.log_marginal_likelihood`,
  celerite2) is 17 %. The model's own `evaluate` is 20 %, three quarters of
  it constructing and validating a fresh `Spectrum` (`results_schema.py`
  `_build_axes`/`_spacing`/`allclose`).
* **torch realised, 3.50 ms.** `backward()` is 34 %; the prior
  (`log_prior_unconstrained_tensor`, torch.distributions plus the
  bijector's Jacobian) 31 %; the GP 19 %. Small-tensor dispatch overhead,
  not arithmetic; nothing in ampere's code is redone.
* **jax contract, 31.3 ms.** 300 calls, **300 XLA compilations**
  (`compiler.py:backend_compile_and_load`, 20 ms each, plus lowering).
  The cause is inside celerite2.jax: `GaussianProcess._do_compute` builds
  `lax.cond(bad, _bad, _good)` whose `_good` branch closes over the freshly
  computed log-determinant, a concrete array. Called eagerly — which is what
  `ampere.backends.jax.QuasisepGP.log_marginal_likelihood`, the numpy-facing
  contract path, does — the closure's constant is baked into the branch
  jaxpr, so every call presents jax with a new jaxpr and a compile-cache
  miss. The realised path is unaffected because it runs under `jit`, where
  the closure captures a tracer. The fit is **98 % compile**.
* **jax realised, 0.283 ms.** Nothing to take out at this size.

### 3.3 Interferometry

* **Reference, 1.19 ms.** `FourierSample.apply` (`ampere/backends/reference/
  interferometry.py`) is 33 % of the profiled loop. Per call it recomputes
  `cell_solid_angle(x, y)` (the grid's quadrature weights: widths, an outer
  product, a unit-free scale), and the two DFT factors
  `exp(-2πi · outer(v, y))` and `exp(-2πi · outer(u, x))` — every one a
  function of the model grid and the (u, v) sampling alone, both fixed after
  compilation — and it forms a composite astropy unit (`unit * u.sr`) and
  compares it with the template's. The priors are 25 %, the complex GP's
  dense solve 13 %.
* **torch realised, 8.10 ms**: `backward()` 50 %, `native_flux`
  (the torch model's DFT) 17 %. **jax realised, 0.473 ms**: compiled.

### 3.4 Image

* **24², 18.2 ms.** `DenseGP.log_marginal_likelihood` is 78 %: the kernel
  matrix 60 % (`kernels.py` `separation` 30 %, `exp` 22 %, `_covariance`
  7 %), the Cholesky factor 9 %, the solve 5 %. `PSFConvolution.apply` 8 %,
  the priors 5 %.
* **64², 1 789 ms.** `cho_factor` 41 %, **`cho_solve` 31 %**, the kernel
  matrix 25 %. A solve with one right-hand side is O(N²) against the
  factor's O(N³), so a solve costing three quarters of the factorisation is
  not arithmetic. It is a layout copy: scipy 1.18's batched `cho_factor`
  returns a **C-ordered** factor, and `cho_solve` hands it to LAPACK
  `potrs`, which f2py must first copy into Fortran order — a transposing
  copy of an N × N matrix, 128 MB at N = 4 096, which the measurement puts at
  ~0.9 s. Handing `potrs` the transpose instead (an F-ordered view of the
  same memory, with `lower=False`) does the same solve with no copy
  (§5, L2).
* **torch realised 34.6 ms, jax realised 45.9 ms** at 24² (value and
  gradient of a dense 576-point GP): backward/XLA arithmetic, no ampere
  overhead worth a lever.

### 3.5 Resampling — issue #12

None of the three loads resamples, so #12 is measured directly: the
reference `Resample` on an evenly sampled `Spectrum`,

| Input → output | per call | of which the weight matrix | the product |
|---|---|---|---|
| 2 000 → 200 | 8.52 ms | 3.97 ms (47 %) | 0.042 ms |
| 20 000 → 2 000 | 471 ms | 433 ms (92 %) | 13.5 ms |

`Resample.apply` calls `influence(samples.spectral_axis.values)` on every
evaluation, rebuilding the dense `(n_out, n_in)` overlap matrix from bin
edges although, after negotiation, the input grid is the compiled model's
and does not change between draws. The v2 core has **not** inherited
legacy's SpectRes loop, but it has inherited the issue's shape: the
expensive half of resampling is planning, and it is re-planned per call.
The jax `Resample` already avoids this — its `_influence_jax` caches the
matrix keyed on the input grid's exact bytes, citing W2.3's measurement of
exactly this waste in the reference backend — so the reference class is
the one that inherited the problem. The torch `Resample` builds its matrix
in torch (`influence_tensor`) on every call too (§6, P4).

### 3.6 The three issues, against the v2 core

* **#12 (resampling).** Inherited in shape, not in implementation: see §3.5.
  Lever L3.
* **#29 (covariance solves).** Largely answered by construction —
  `QuasisepGP` is O(N) and exact for the quasiseparable kernels (the M2
  rows: 0.57 ms at 10⁴ points, 3.0 ms at 2 × 10⁴), and W5.4/W5.6 ship
  `HilbertSpaceGP`, EFGP and Vecchia for the rest. What the profile adds is
  that the *dense* anchor was paying a layout copy as large as its
  arithmetic (§3.4). Lever L2.
* **#67 (loops pretending to be vectorised; distribution evaluation).**
  Two halves. The engines do not pretend: emcee, zeus and dynesty call the
  scalar contract path once per proposal and none claims otherwise (the
  nested driver passes `vectorized=False` explicitly), and the realised path
  is where batching lives (`vmap` on jax). What *is*
  inherited is the issue's other observation — "a lot of time is spent …
  evaluating probability distributions": the six scalar scipy `logpdf`
  calls are the largest single item on the M2 reference load. §6, P1.

## 4. Candidate levers not taken, and why

* **Solver selection.** The image study's flexible arm uses `DenseGP`, the
  exact anchor, by choice; the bake-off (W5.6) already measures the
  approximate solvers against it. Switching a default to an approximate
  solver changes answers, which no optimisation pass may do. M2 is already
  on `QuasisepGP`.
* **Precision policy.** Float64 stays the default for likelihood linear
  algebra (§7's trap). The opt-out is already reachable (the float32 route
  recorded at the x64 decision); no load spends measurable time on
  precision bookkeeping, so there is nothing to make cheaper.
* **GPU placement.** Out of scope: no accelerator on this machine.

## 5. The levers, ranked by measured share

Each lever is an **exactly-zero** change: the same floating-point
operations on the same operands, so every conformance row keeps its value
bit for bit. That is checked per lever, not assumed.

| # | Lever (plan category) | Where | Measured share | Rows it should move |
|---|---|---|---|---|
| L1 | jax contract path: compile once, not per call (compilation caching) | `ampere/backends/jax/gp.py`, `QuasisepGP.log_marginal_likelihood` / `_factorise` | ~98 % of the jax contract loop | `test_backend_quasisep_contract[*-jax]`, `test_one_proposal[contract-QuasisepGP-jax]` |
| L2 | Dense solve without the layout copy (#29) | `ampere/core/likelihood.py`, `DenseGP` and every other `cho_solve` on a real factor | 31 % of image 64², 5 % of image 24² | `test_reference_dense[2000]`, `test_reference_dense`, `test_dense_anchor`, not `test_route_cholesky` (EFGP's normal matrix is complex Hermitian, which falls through) |
| L3 | `Resample`: plan the weights once per input grid (#12) | `ampere/backends/reference/instrument.py`, `Resample.apply` | 47–92 % of the step | a timed loop (§3.5); no benchmark row resamples |
| L4 | `FourierSample`: derive the DFT factors and solid angles once per grid (compilation caching; the units trap) | `ampere/backends/reference/interferometry.py`, `FourierSample.apply` | ~33 % of the interferometry loop | a timed loop (§3.3) |
| L5 | The M2 models: build the container template once on the un-negotiated path (re-allocated per call) | `examples/m2_misspecification/model*.py`, `evaluate` | ~15 % of the M2 reference loop | `test_reference_standard[*]`, `test_reference_quasisep[*]`, `test_reference_emcee[*]` |

## 6. Proposals — measured, not landed

Measured, on the candidate list, and not landed — each for the reason
given. The decision-log row each would need is drafted in the W5.17 report.

**P1 — scalar prior evaluation (issue #67's "evaluating probability
distributions").** The largest single item on the M2 reference load:
`ParameterSet.lnprior` is **262 µs of a 701 µs `log_prob`** (37 %; six
scalar priors, four uniform and two half-normal, through scipy frozen
`logpdf`). The arithmetic underneath — the distribution's private
`_logpdf` on the standardised value, minus `log(scale)` — costs **20.5 µs**
for all six, and agrees with `logpdf` bit for bit on all 300 × 6 draws
measured. So an 11× faster prior, and ~1.5× on this `log_prob`, is there.
Not landed because the only exactly-zero route calls scipy's private
`_parse_args`, `_argcheck`, `_support_mask` and `_logpdf`, reproducing
`rv_continuous.logpdf`'s wrapper; a private-API coupling across the CI
matrix's three Pythons (and whatever scipy each resolves) is a maintenance
decision, not an optimisation. The public alternative — ampere's own
closed forms for its common priors — is not bit-identical, so it would move
conformance values within tolerance, which this pass may not do. Either
needs a ruling; `parameters.md` §4 (a prior is any object satisfying the
`Prior` protocol, evaluated through `logpdf`/`logpmf`) and §13's
conformance row (`lnprior` against summed scipy `logpdf`) are the texts a
closed-form route would touch.

**P2 — batched evaluation (issue #67's shape).** The contract path is
scalar by design (`log_prob(values) -> float`), every problem on it
declares `batchable=False`, and the engines call it once per proposal
without pretending otherwise. Batching exists where the plan put it: the
realised path, where jax `vmap`s. Two facts bound it: the reference
backend has no batched evaluation to call, and `ampere.backends.jax.
QuasisepGP` declares `BATCHABLE = False` because celerite2's primitives
register no `vmap` rule (`DEVELOPMENT_PLAN.md` §2, the jax quasiseparable
row). A batched contract path would be a §4.5 addition (a stacked
`log_prob`); nothing here measured what it would buy on a numpy model, so
it is recorded, not ranked.

**P3 — solver selection on the image.** The image study's flexible arm
defaults to `DenseGP`, and at 64² that is the whole cost (§3.4). The
baseline's own rows say what the alternatives cost at N = 2 304:
`HilbertSpaceGP` (m = 256) 80 ms, Vecchia (k = 30) 207 ms, EFGP
(m = 289) 249 ms. Choosing an approximate solver by default changes
answers, so it is a science decision (W5.6's bake-off is where it is
made), not an optimisation.

**P4 — the torch `Resample`.** `influence_tensor` rebuilds the weights in
torch on every call, as the reference class did before L3. The same
exact-bytes cache would apply, but the matrix there can sit inside an
autograd graph, and whether a cached tensor may be reused across graphs
without `detach` is a torch-backend question this pass did not measure.

**P5 — the jax contract path's remaining overhead.** After L1 the M2 jax
contract `log_prob` is ~4–6 ms against 0.28 ms realised: eager dispatch of
~30 small jax operations per call, plus `float()` synchronisations.
Jitting the whole contract-path solve takes it to ~2 ms but is **not**
bit-identical (XLA fuses, measured max difference 8e-11 on the same
draws), so it cannot land under this pass's rule. The realised path is
already the fast route for anyone using jax.

## 7. Results

Five levers landed, one commit each, every one **bit-for-bit**: each
before/after pair below was measured in one process on the same draws,
and the values compared with `np.array_equal`, not a tolerance. After
each, `tests/conformance` in dev (and the lever's own test files) stayed
green; torch and jax conformance ran once at the end (§7.3).

### 7.1 The timed loops (one process, same machine, values bitwise equal)

| Lever | Load | before | after | × |
|---|---|---|---|---|
| L1 jax cond compiled once | M2 jax contract, Matérn-3/2 (200 draws) | 30.89 ms | 4.51 ms | 6.8 |
| | M2 jax contract, `Sum` | 31.80 ms | 6.32 ms | 5.0 |
| | M2 jax contract, warped | 36.30 ms | 11.61 ms | 3.1 |
| L2 dense solve, no layout copy (#29) | image 24², flexible (100 draws) | 18.90 ms | 17.46 ms | 1.08 |
| | image 48², flexible (10 draws) | 261.8 ms | 232.7 ms | 1.12 |
| | image 64², flexible (6 draws) | 1 970 ms | 760 ms | 2.59 |
| | interferometry, flexible (200 draws) | 1.59 ms | 1.12 ms | 1.42 |
| L3 `Resample` planned once (#12) | 2 000 → 200 | 8.856 ms | 0.137 ms | 65 |
| | 20 000 → 2 000 | 450.0 ms | 15.0 ms | 30 |
| L4 `FourierSample` planned once | interferometry correct (300 draws) | 1.328 ms | 0.641 ms | 2.07 |
| | interferometry incomplete | 0.808 ms | 0.536 ms | 1.51 |
| | interferometry flexible | 1.201 ms | 0.878 ms | 1.37 |
| | interferometry chromatic, product kernel (100 draws) | 7.855 ms | 5.924 ms | 1.33 |
| L5 M2 template on the simple path (min of 7 interleaved rounds) | M2 reference standard | 0.436 ms | 0.322 ms | 1.35 |
| | M2 reference flexible | 0.691 ms | 0.557 ms | 1.24 |
| | M2 torch standard / flexible | 0.504 / 1.054 ms | 0.390 / 0.924 ms | 1.29 / 1.14 |
| | M2 jax standard / flexible | 0.887 / 3.774 ms | 0.767 / 3.858 ms | 1.16 / — |

L5's jax flexible row is the one "tried, no gain" in the table: after L1
the jax contract path's remaining eager overhead (P5) dwarfs the
container, and the difference is inside the scatter. It is kept because
the same one-line change is a measured gain on the other five rows and the
three M2 models share it by construction (`model._template_for`).

The M2 reference `log_prob` at 200 points, the study's own first rung,
went from ~0.70 ms to ~0.56 ms with L5; the rest of it is now the priors
(P1) and the quasiseparable solve.

### 7.2 The benchmark rows

`pixi run -e <env> bench` at the branch head, same machine, same lock
discipline as the baselines (§2); JSON at
`~/.cache/ampere-gates/w5.17/after-<env>.json`. Medians in ms, baseline →
head. The rows each lever targets:

| Row | Lever | dev | torch | jax |
|---|---|---|---|---|
| `test_backend_quasisep_contract[200-jax]` | L1 | — | — | 30.2 → **3.65** |
| `test_backend_quasisep_contract[2000-jax]` | L1 | — | — | 33.0 → **6.38** |
| `test_backend_quasisep_contract[20000-jax]` | L1 | — | — | 31.7 → **7.45** |
| `test_one_proposal[contract-QuasisepGP-jax]` | L1 | — | — | 28.5 → **4.61** |
| `test_backend_quasisep_realised[*-jax]` | (untouched) | — | — | 0.29 / 0.68 / 5.08 → 0.29 / 0.63 / 4.82 |
| `test_dense_anchor` (image `DenseGP`, n1024) | L2 | 208 → **83.8** | 199 → **39.1** | 63.7 → 256 |
| `test_reference_dense[2000]` | L2 | 308 → 209 | 497 → 315 | 248 → 171 |
| `test_reference_standard[200]` | L5 | 0.476 → 0.339 | 0.437 → 0.317 | 0.450 → 0.323 |
| `test_reference_quasisep[200]` | L5 | 0.722 → 0.590 | 0.735 → 0.605 | 0.717 → 0.578 |
| `test_reference_emcee[flexible]` | L5 | 1716 → **1435** | 1773 → **1391** | 1977 → 1517 |

**How far to believe them.** The jax contract rows move by 4–8× against an
IQR of 1–3 ms: L1 is unambiguous. The M2 reference rows move by 1.2–1.4×
in all three environments, consistently, against IQRs of 0.07–0.2 ms, and
`test_reference_emcee[flexible]` — the one sampling row with a tight IQR
(14–26 ms) — moves by 1.2–1.3×: L5 is real on the benchmark too. The
dense rows are where the shared machine shows: `test_dense_anchor` runs
the *same numpy code* in all three environments, and its baseline alone
spans 64–208 ms between them; its head moved 2.5× and 5.1× down in dev and
torch and 4× *up* in jax. Rows no lever touches moved just as far in the
same runs (the legacy dense RBF row by anything from 0.06× to 9×; EFGP,
whose normal matrix is complex and so falls through L2, by up to 3×).
L2's evidence is therefore the in-process timed loop (§7.1), not these
rows.

**A regression that was not one.** The torch `test_one_proposal` rows
came back 1.4–4× slower at the head. None of the five levers is on the
torch realised path, so this was re-measured as an A/B in one session:
`tests/benchmarks/test_engine_fast_path.py` run alternately against an
export of the base commit's tree (shadowing the editable install through
`PYTHONPATH`) and the head, twice each, in torch and jax. torch
`contract-DenseGP`: base 8.00 / 6.47 ms, head 4.88 / 7.89 ms;
`realised-DenseGP`: base 8.18 / 5.69, head 6.19 / 8.56 — the base itself
had slowed by the same factor, so it was the machine, not the code. The
same A/B confirms L1 in isolation: jax `contract-QuasisepGP` base
35.6 / 33.4 ms, head 5.39 / 5.81 ms.
`test_the_fast_path_beats_the_contract_path` (a ratio, realised against
contract, floor 3×) still passes on jax: the contract path is now ~7× the
realised one rather than ~40×, comfortably above the floor.

### 7.3 Conformance and the gates

At the head, every conformance row passes in all three environments —
dev 669 passed / 75 skipped, torch 1 073 / 75, jax 1 073 / 75 — and no
row's value moved: every lever was checked bit-for-bit against its
unchanged twin on the loads above before it was committed. Collected
tests: 3 869 before and after (plus the one pre-existing legacy collection
error, `ampere/models/test_Dusty.py`, which needs the absent `Dusty`
module and is not this item's). Lint, format-check, the import sweep and
pyrefly (dev, torch and jax: 0 errors) are clean. No §4 contract and no
frozen `docs/design/contracts/` page changed; float64 remains the default
everywhere; no GPU placement was attempted.

**Effect on W5.15's measurements.** L4 speeds up the reference
interferometry forward model that any fit sharing it would use; values
are bit-identical, so only timings move.

## Appendix A. The harness

The profiling loop, run as `pixi run -e <env> python prof.py <load> <backend>
[n] [pixels] [--realised] [--build]` from the repository root:

```python
import cProfile, io, pstats, sys, time
sys.path[:0] = [".", "tests/backends"]          # examples/, interferometry_fixtures
import numpy as np

def build(load, backend, pixels=None):
    if load == "m2":
        from examples.m2_misspecification import study
        return study.build_problem("strong_smooth", size=200, backend=backend,
                                   likelihood="flexible")
    if load == "itf":
        from examples.interferometry import study
        return study.build_problem(backend, "flexible")
    from examples.image import study
    solver = None if backend == "reference" else study._backend_module(backend).DenseGP()
    return study.build_problem(backend, "flexible",
                               pixels=pixels or study.SMALL_PIXELS, solver=solver)

problem = build(load, backend, pixels)
rng = np.random.default_rng(1)
draws = [problem.sample_prior(rng) for _ in range(n)]
problem.log_prob(draws[0])                       # warm-up
profiler = cProfile.Profile()
profiler.enable(); values = [problem.log_prob(d) for d in draws]; profiler.disable()
pstats.Stats(profiler).strip_dirs().sort_stats("cumulative").print_stats(40)
# --realised: theta = problem.unconstrain(d); torch: backward() on
# realise(problem).log_prob_unconstrained(tensor); jax: jit(value_and_grad(...)).
```

(The image study's flexible arm defaults its solver to `ampere.core.DenseGP`,
which a native torch or jax problem refuses as a foreign part; the harness
passes the backend's own `DenseGP`. Recorded as an out-of-scope finding in
the W5.17 report.)
