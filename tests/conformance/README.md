# The conformance suite

`DEVELOPMENT_PLAN.md` §4.6 calls this "the contract that keeps two lockstep
backends aligned across agent tracks", and §7 records why it has to exist:
*lockstep backends drift without a mechanical forcing function*. A feature is
done when **both** backends pass this battery, not when one of them works.

Run it with:

```console
$ pixi run conformance          # or: pixi run -e dev python -m pytest tests/conformance -q
```

Every test in this directory takes the `backend` fixture and runs once per
registered backend. No test body names a concrete backend, and none may.

---

## 1. Adding a backend

Two steps.

1. Write a class satisfying `ConformanceBackend` (§2 below). Put it in
   `tests/conformance/backends/`.
2. Add one entry to `_REGISTRY` in `tests/conformance/backends/__init__.py`,
   wrapped in `_optional(...)` if it needs a dependency the base install does
   not have.

That is the whole of it. The battery then runs against your backend with no
other change — which is W1.10's acceptance criterion, and the reason there are
two in-repo fixtures rather than one (see §5).

---

## 2. The fixture protocol

`tests/conformance/protocol.py` is the authority; this section is the prose.
Seven members.

### `name: str`

The pytest fixture id. Short, lower-case, stable across runs — it appears in
every test id and in cross-backend row names (`reference-vs-mirror`).

**It is also the backend's name in the sense W2.12 fixed**: one name per
backend, everywhere. Every model and instrument step this fixture builds must
declare `BACKEND` equal to this string, because `TestBackendIdentity` in
`test_inference.py` asserts that a composed problem reports it — and because
that is the same string `lowering.md` §12.8's registry is keyed on and the same
string a run records as `ampere_backend`. A fixture that reuses another
backend's classes must subclass to change the declaration (`MirrorResample` and
its siblings are the worked example); reusing them under a different fixture
name would make the row a tautology and would put a name in provenance that
nothing else in the project answers to.

### `capabilities: BackendCapabilities`

A frozen record:

| field | meaning |
|---|---|
| `differentiable`, `batchable`, `device` | three of the four flags `inference.md` §18 says a backend declares. Mirror `ampere.core.Capabilities`. The fourth, `backend`, is not repeated here: it *is* `name` above. |
| `float64` | whether the likelihood linear algebra runs in double precision. `architecture.md` §5 makes float64 the policy for GP solves; a backend that opts out for GPU throughput says so here and widens `tolerances.cross_backend`. |
| `solvers` | the `SolverKind`s `gp_solver` can return an *implemented* solver for. Rows for absent kinds skip with a reason naming what is owed. |
| `complex_models` | whether `model()` can realise `ModelKind.COMPLEX` — a channel of complex values, which the `complex_gaussian` rows need. `False` by default; those rows then skip with a reason. Added at W2.4 slice 3, and made a capability rather than a required protocol member so that one track's slice is not a change to every other track's fixture. |
| `tolerances` | the per-comparison table (§3). |

### `model(spec: ModelSpec) -> Model`

Realise a model declaration. The returned object must be an
`ampere.core.Model` **and** a `CountingModel` — that is, it must expose an
`evaluations: int` counter incremented once per `evaluate` call, and a
`reset_evaluations()` that zeroes it. (`inference.md` §18 asks for a row
proving an out-of-support θ is refused *without* running the model; there is no
way to assert that from outside without the model saying so.)

`ModelSpec` carries:

* **`kind`** — one of three closed forms, so every downstream oracle is
  analytic:
  * `LINEAR`: `f(x) = offset + slope * x`, with `offset ~ Normal(0, 1)` and
    `slope ~ Normal(1, 0.5)` — both unbounded, so both lower to `Identity`.
  * `POWER_LAW`: `f(x) = norm * (x / x_ref) ** index`, with
    `norm ~ LogUniform(0.1, 10)` (bounded both ways → `Logit`) and
    `index ~ Normal(-1, 0.5)` (`Identity`).
  * `COMPLEX`: `V(x) = norm * exp(i * index * x)`, the same two parameters
    under the same priors, emitted as a **`VisibilitySet`** whose second (`v`)
    axis comes from `protocol.complex_axes` — the same rule the observed
    container uses, because `check_alignment` compares the two for equality.
    Optional: only asked for if you declare `capabilities.complex_models`. Note
    that the modulus does not depend on `index`, deliberately — a backend that
    dropped the imaginary part would score a density independent of one of its
    own parameters, which the row then catches by nats rather than by digits.
  Declare exactly those parameters, under exactly those names. The `lnprior`,
  `prior_transform` and spec-hash rows compare by name across backends, so a
  renamed parameter is a failure, not a detail.
* **`channels`** — one `Spectrum` per named channel in the returned
  `ModelResult` (a `VisibilitySet` for `COMPLEX`).
* **`coordinates`** — the model's own grid, in **micron**; emitted flux is in
  **Jy**. Units are fixed here so containers stay comparable across backends
  without making unit negotiation a backend author's problem.
* **`reference_coordinate`** — `x_ref` for the power law.
* **`plated`** — when `True`, additionally declare a `Plate` named `objects`
  of size `len(channels)`, with hyperparameters `mu ~ Normal(0, 1)` and
  `sigma ~ HalfNormal(1)` and one array-valued member `offsets` drawn from
  `Normal(mu, sigma)`; channel *i* adds `offsets[i]` to its flux. Register the
  *expanded* parameters (`Plate.expand()`), so `parameter.plate` survives to
  the emitted `InferenceData` and names the dimension.
* **`compiled`** — when `False` (the default), leave `compile_for` inherited:
  it returns `self` and the model evaluates on its own grid. When `True`,
  honour the request — build one container per channel from
  `requirements[channel].coordinates()`, keep it, and refill it with
  `with_values` on every evaluation, so successive results **share their axes
  by identity**. Both behaviours are conformance rows.

### `transformation(spec: TransformationSpec) -> Transformation`

Three kinds:

* **`SCALE`** — a multiplicative calibration. `ACCEPTS = (Spectrum,)`,
  `PRODUCES = None`, one free parameter `scale ~ LogNormal(0.2)` (support
  `(0, ∞)`, so a `Log` bijection). Publishes no requirements.
* **`REBIN`** — resamples a `Spectrum` onto `spec.target` through an explicit
  influence matrix, exposed as `.influence(source_coordinates) -> (n_out,
  n_in)`. Publishes an `AxisRequirement` on `spectral_axis` covering the target
  grid, and **must** propagate the input mask with `propagate_mask`.
* **`PHOTOMETRY`** — `Spectrum → PhotometricPoints`, integrating onto
  `spec.target` under the names in `spec.filters`. `PRODUCES =
  PhotometricPoints`. It exists so a row can build a chain whose kinds do not
  compose.

### `kernel(spec: CovarianceSpec) -> Kernel`

`MATERN32` or `SQUARED_EXPONENTIAL`, with `amplitude` (the marginal standard
deviation: `k(0) == amplitude²`) and `length_scale` in the coordinate's own
units. Both hyperparameters arrive as plain numbers and are held fixed.

The kernel's `QUASISEPARABLE` class flag must be truthful — the
`DenseGP`↔`QuasisepGP` row selects on it, and `GPSolver.check_compatible`
refuses a quasiseparable solver a kernel that has no such representation.

### `gp_solver(kind: SolverKind) -> GPSolver`

`DENSE` or `QUASISEP`. Only called for a kind present in
`capabilities.solvers`; you may raise for anything else.

### `parameter_space(declaration: ParameterSet) -> ParameterSpace`

Your engine-facing view of a parameter declaration. `ampere.core.ParameterSet`
already satisfies every method, so the reference adapter is a thin wrapper; a
backend that lowers a declaration into native sample sites (numpyro, pyro)
returns its own object, and `test_parameters.py` then holds it to the same
behaviour — which is the point.

Required surface: `declaration`, `free_size`, `free_labels()`, `pack`,
`unpack`, `prior_transform`, `lnprior`, `constrain`, `unconstrain`,
`lnprior_unconstrained`. Every method takes and returns plain numpy on the
boundary; convert at the edge if you work in another array type.

### `to_numpy(values) -> np.ndarray`

Bring a backend array back to numpy for comparison against an oracle. Only the
edges call it, so no row has to know your array type.

### What a backend does **not** supply

`Instrument`, `Likelihood`, `Dataset`, `DatasetCollection` and
`FittingProblem` are reused unchanged, assembled by
`tests/conformance/composition.py` from your pieces. That is the claim
`inference.md` §18 and `results.md` §14 both make — "a backend supplies models
and transformations […] and nothing else", "a backend supplies draws and
`Evaluation`s and nothing else" — and this suite is where it gets tested.

---

## 3. Tolerances

`architecture.md` §2 point 1 commits the project to tolerances "documented and
specific per comparison, never a single blanket tolerance", and hands the table
to W1.10. It is `Tolerances`, and it lives on the backend so that a float32 or
different-accumulation-order backend can widen its own entry without loosening
anyone else's.

| field | default | what it covers |
|---|---|---|
| `exact` | `1e-12` | identities that hold in exact arithmetic: `pack`∘`unpack`, `constrain`∘`unconstrain`, a spec round trip, mask excision equalling deletion. Anything looser here is a bug, not a tolerance. |
| `analytic` | `1e-9` | a closed form on the *same* arithmetic path: the i.i.d. Gaussian written out, `scipy.stats.norm.logcdf` for a Tobit limit, a prior quantile from `ppf`. |
| `linear_algebra` | `1e-8` | a closed form through a *different* factorisation: `DenseGP` against `scipy.stats.multivariate_normal.logpdf`, `L Lᵀ` against `K`. |
| `cross_solver` | `1e-6` | two `GPSolver` strategies by genuinely different recursions — dense Cholesky against the quasiseparable state-space solve. |
| `cross_backend` | `1e-9` | two backends' independent arithmetic, and in Phase 2 a different autodiff accumulation order. Cross-backend rows use the **looser** of the pair. |
| `monte_carlo_sigmas` | `5.0` | not a tolerance but a multiple: how many standard errors an empirical moment may sit from its analytic value. The standard error is computed from the estimator, so the assertion stays honest as the draw count changes. |

---

## 4. Layout

```
tests/conformance/
  protocol.py            the fixture protocol and the declaration types
  composition.py         ProblemSpec -> FittingProblem, plus the deterministic data
  oracles.py             the closed forms every row compares against
  conftest.py            the `backend`, `tolerances` and `backends` fixtures
  backends/
    __init__.py          the registry — the one place a backend is named
    reference.py         an adapter onto the shipped ampere.backends.reference
    mirror.py            a second fixture (§5)
  test_parameters.py     parameters.md §13's rows
  test_transformations.py transformations.md §14's rows
  test_likelihoods.py    likelihoods.md §16's rows
  test_inference.py      inference.md §18's rows
  test_results.py        results.md §14's rows (arviz-gated)
  test_schema.py         results_schema.md §16's rows
  test_cross_backend.py  the rows that compare two backends
```

Two rules for writing a row:

* **Oracles are analytic or `scipy`, never a second copy of ampere's own
  arithmetic.** If you find yourself importing an ampere function to compute
  the expected value, the row has stopped being a test.
* **A row whose second implementation does not exist yet is a skeleton, not an
  absence.** Present, parametrised, and `skip`ped or `xfail`ed with a reason
  string naming exactly what is owed — so `pytest tests/conformance -q` shows
  the shape of the debt.

---

## 5. Why there are two in-repo backends

`mirror` is not a backend in `architecture.md` §1's sense and must not be
mistaken for one. It exists because a battery with a single fixture cannot
demonstrate that it is parametrised at all: the machinery would be untested,
and the first Phase-2 author would discover the hard way which rows had quietly
baked in an assumption about the reference path.

It shares the reference backend's transformations, kernels and solvers, because
it has nothing different to offer there. What it supplies of its own is exactly
the two things a real backend owns: **models** — the same closed forms under
the same parameter names, computed by a deliberately different route (Horner's
rule; `exp` of a logarithm) in classes of their own — and **a parameter space**
that round-trips the declaration through `to_spec`/`from_spec` before
delegating, so every parameter row is also a serialisation row on that fixture.

---

## 6. What is currently owed

Run the suite to see it; as of W1.10:

* ~~**`DenseGP`↔`QuasisepGP` agreement on Matérn-3/2** — `DEVELOPMENT_PLAN.md`
  §4.6 names it and `ampere.core.QuasisepGP` is a declared strategy slot with
  no implementation.~~ **Discharged by W2.3**: `ampere.core.QuasisepGP` is
  celerite2's numpy solver over an exact rank-2 representation of Matérn-3/2,
  the in-repo fixtures declare `SolverKind.QUASISEP`, and the marginal- and
  conditional-agreement rows run at `tolerances.cross_solver`. The
  no-fallback row now asserts the positive claim — a declared quasiseparable
  strategy must be implemented, exact, and a *different* strategy from the
  dense one, since a backend that declared `QUASISEP` and handed back
  `DenseGP` would satisfy every agreement row and prove nothing. Its old half
  (asking for the empty slot and asserting it refused) moved to `tests/core`:
  that is a property of `ampere.core`'s declared-slot discipline, not of a
  backend, and `gp_solver` is contractually called only for a kind the
  backend declares. A backend that has not yet supplied one (torch, jax)
  skips these rows with a reason naming what to supply (`tinygp`'s
  `QuasisepSolver`, GPyTorch or celerite2-torch on the torch side).
* **`ampere_problem_hash` agreement across backends** — a recorded
  contradiction, not a skip. `results.md` §14 promises all three hashes are
  functions of the declaration and the data, "not of the arithmetic", so a
  cross-backend disagreement can only be a lowering bug. The spec and data
  hashes hold that promise and are live rows. The problem hash does not:
  `provenance.model_fingerprint` deliberately records the model's class and
  module — with a docstring arguing, correctly, that two models of different
  classes can declare the same parameters and compute completely different
  things — so two backends implementing one declaration disagree by
  construction.

  Both behaviours are defensible and W1.10 has no standing to change either,
  so `test_the_problem_hash_tracks_the_model_class_not_only_the_declaration`
  asserts the *actual* rule as an equivalence: the problem hash agrees exactly
  when the recorded model identity does. It is green today, green for a
  backend that reuses another's model classes, and fails loudly the moment
  either the contract or the fingerprint changes. A ruling is owed on whether
  §14's backend-spanning claim should be narrowed to the spec and data hashes,
  or the fingerprint should carry a backend-neutral model identity.

  **Ruled 2026-09-03: the narrowed claim stands.** The spec and data hashes
  are the backend-spanning promise; `ampere_problem_hash` remains
  deliberately backend-variant (a lowered model is a different
  implementation, and different implementations must not share a cache
  key). Whether Phase 2 additionally wants a backend-neutral model identity
  — so an emulator trained on one backend can be *offered*, never silently
  served, to a fit on another — stays routed to the freeze beside
  `results.md` §13.13's `describe()` hook.

  **Ruled at the freeze's escalations (2026-09-03): adopted for early
  Phase 2**, as one mechanism with the `describe()` hook — both land with
  W2.1 (decision-log entry there). The derived neutral identity (the
  fingerprint minus class/module) offers, never serves; the problem hash
  itself stays backend-variant, so this suite's equivalence row is
  unchanged. (W1.13's report carries a
  recommendation; the ruling is Peter's.)

  **Landed W2.1, and the debt is discharged rather than removed.** The
  equivalence row is untouched, as the ruling required.
  `test_the_neutral_model_identity_agrees_across_backends` is its
  counterpart: `ampere_model_identity_hashes` must agree across every
  registered backend, unconditionally. Together they state the whole rule —
  the neutral identity is what licenses *offering* a cross-backend emulator,
  and the backend-variant problem hash is what stops one being served. That
  the pair is a real assertion rather than a tautology is why `mirror`
  exists: same declaration, different classes, different modules.

Added at the freeze (W1.13):

* ~~**The X-1 rows run live** against an in-repo fractional-noise double
  (`test_likelihoods.py::TestPredictionAwareNoise`); Phase 2 points them at
  the shipped reference-backend `FractionalModelNoise` when it lands.~~
  **Discharged by W2.1**: the three rows now exercise
  `ampere.backends.reference.FractionalModelNoise` and
  `FractionalModelGPNoise` themselves, and the doubles are deleted — a
  double that shadows a shipped class only tests itself.
* **`complex_gaussian` + GP is declared `ANALYTIC` and staged**
  (`TestStagedAnalyticCombination`): the declaration row runs, the
  composition-refusal row runs, and Phase 4 replaces the refusal with
  agreement rows against the circular closed form.

  **The uncorrelated half landed at W2.4 slice 3** and the staged half is
  unchanged. `ModelKind.COMPLEX` and the `COMPLEX` shape in
  `test_inference.py` give the family a *realised* row —
  `test_the_complex_realisation_agrees_with_the_numpy_path`, at
  `tolerances.cross_backend`, on any backend declaring
  `capabilities.complex_models` — and
  `test_a_complex_gp_is_refused_by_name_on_every_path` pins the invariant that
  makes the staging safe: the contract path refuses the composition, so no
  backend can be quietly computing something for it. When Phase 4 implements
  the circular complex GP in `ampere.core`, that refusal row becomes the
  agreement row and each backend's `refuse_family` drops its `correlated`
  branch.
* ~~**`GPSolver.conditional_loo`** outlived the QuasisepGP debt: `DenseGP`
  implements the leave-one-out terms and `QuasisepGP` refuses, W2.3 having
  **deferred** the O(N) recursion with a decision-log entry
  (`DEVELOPMENT_PLAN.md` §2, 2026-09-05) — celerite2's public numpy interface
  exposes no O(N) route to the diagonal of `(K + diag(σ²))⁻¹`. The refusal is
  a live row (`TestSolverAgreement`); the eventual agreement row mirrors the
  marginal one.~~ **Stated in full by W2.4 slice 2**: the row is now
  `test_the_leave_one_out_terms_either_agree_exactly_or_refuse_by_name`, and
  it says the whole contract — a quasiseparable strategy either refuses by
  name or computes exactly the decomposition `DenseGP` does, at
  `tolerances.cross_solver`. `ampere.core.QuasisepGP` still refuses (the
  deferral stands where it was taken: celerite2's *public numpy* interface has
  no route to that diagonal), and `ampere.backends.torch.QuasisepGP` supplies
  the terms, because it calls celerite2's compiled kernels directly and
  therefore already holds the factorisation the O(N) backward recursion needs.
  The decision-log row of 2026-09-07 records the change of circumstances.
  Backends are free to differ here, and the row is what keeps a *wrong*
  implementation from passing for either choice.

  **All three implemented solvers have now made that choice** (W2.5 slice 2
  completes the entry). `ampere.core.QuasisepGP` refuses;
  `ampere.backends.torch.QuasisepGP` supplies the terms;
  `ampere.backends.jax.QuasisepGP` **refuses**, and the reason is worth
  recording rather than leaving as a gap, because it is not the same reason
  the numpy path gives. The jax solver goes through `celerite2.jax`'s public
  `GaussianProcess`, whose surface is `log_likelihood`, `apply_inverse`,
  `dot_tril`, `predict`, `condition` and `sample`: the two members that could
  give the diagonal form the cross-covariance densely and cost O(N·M), and the
  O(N) route would need the private `_d`/`_W` — the coupling the torch backend
  already pays for and tests, and which would be a *new* coupling here. The
  alternative was measured rather than assumed: tinygp's quasiseparable factor
  exposes the diagonal exactly (`L.inv().transpose() @ L.inv()`, matching a
  dense inverse to 7e-15), and tinygp was rejected on a factor of 200 in the
  marginal likelihood itself. So the debt is closed as a *statement* — every
  shipped solver either supplies the terms or refuses by name, and the row
  asserts exactly that — with one open path recorded: a jax
  `conditional_loo` becomes cheap the day this backend calls celerite2's
  kernels directly, as torch does.
