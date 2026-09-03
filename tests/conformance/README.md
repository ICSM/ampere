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

### `capabilities: BackendCapabilities`

A frozen record:

| field | meaning |
|---|---|
| `differentiable`, `batchable`, `device` | the three flags `inference.md` §18 says a backend declares. Mirror `ampere.core.Capabilities`. |
| `float64` | whether the likelihood linear algebra runs in double precision. `architecture.md` §5 makes float64 the policy for GP solves; a backend that opts out for GPU throughput says so here and widens `tolerances.cross_backend`. |
| `solvers` | the `SolverKind`s `gp_solver` can return an *implemented* solver for. Rows for absent kinds skip with a reason naming what is owed. |
| `tolerances` | the per-comparison table (§3). |

### `model(spec: ModelSpec) -> Model`

Realise a model declaration. The returned object must be an
`ampere.core.Model` **and** a `CountingModel` — that is, it must expose an
`evaluations: int` counter incremented once per `evaluate` call, and a
`reset_evaluations()` that zeroes it. (`inference.md` §18 asks for a row
proving an out-of-support θ is refused *without* running the model; there is no
way to assert that from outside without the model saying so.)

`ModelSpec` carries:

* **`kind`** — one of two closed forms, so every downstream oracle is
  analytic:
  * `LINEAR`: `f(x) = offset + slope * x`, with `offset ~ Normal(0, 1)` and
    `slope ~ Normal(1, 0.5)` — both unbounded, so both lower to `Identity`.
  * `POWER_LAW`: `f(x) = norm * (x / x_ref) ** index`, with
    `norm ~ LogUniform(0.1, 10)` (bounded both ways → `Logit`) and
    `index ~ Normal(-1, 0.5)` (`Identity`).
  Declare exactly those parameters, under exactly those names. The `lnprior`,
  `prior_transform` and spec-hash rows compare by name across backends, so a
  renamed parameter is a failure, not a detail.
* **`channels`** — one `Spectrum` per named channel in the returned
  `ModelResult`.
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
    reference.py         ampere.core's numpy path
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

* **`DenseGP`↔`QuasisepGP` agreement on Matérn-3/2** — `DEVELOPMENT_PLAN.md`
  §4.6 names it and `ampere.core.QuasisepGP` is a declared strategy slot with
  no implementation. Two rows skip with a reason naming what Phase 2 must
  supply (celerite2 on the numpy and jax sides, `tinygp`'s `QuasisepSolver`,
  GPyTorch or celerite2-torch on the torch side). A third row runs today and
  asserts the empty slot *refuses* rather than silently falling back to the
  dense path — without it the agreement row would be vacuous the day someone
  forgot to implement it.
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
