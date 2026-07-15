# Ampere v1 Implementation Plan

## Background and Motivation

The current ampere codebase provides a working SED/spectrum fitting framework
built around:

- `ampere.models` – physical models (blackbodies, power laws, RT wrappers, …)
- `ampere.data` – data containers with built-in likelihood evaluation
- `ampere.infer` – inference backends (emcee, dynesty, zeus, SBI/PyTorch)

Two structural problems block the next stage of development:

1. **Parameter management is entirely manual.** Each model counts its own
   `npars`, slices raw `theta` numpy arrays positionally, and delegates
   `lnprior` / `prior_transform` through the same positional slicing. This
   makes models hard to compose, hard to introspect, and essentially
   incompatible with gradient-based / probabilistic-programming backends
   (numpyro, blackjax, …).

2. **`Data` classes conflate too many concerns.** The current `Data.lnlike`
   signature is `lnlike(self, synWave, synFlux)`, which bakes in both the
   data-generating process and the likelihood into a single monolithic method.
   This was identified in [issue #73][i73] as the primary barrier to making
   data classes flexible enough to support, e.g., IFU cubes, censored data, or
   user-supplied noise models.

Ampere is already structured around a **dependency-injection** (DI) pattern
(models and data objects are injected into inference objects), but this can be
taken further. By extending DI into the `Data` class itself — injecting the
data-generating process, noise model, and likelihood as separate objects — we
gain far greater flexibility without sacrificing backward compatibility.

[i73]: https://github.com/ICSM/ampere/issues/73

### Goals for v1

1. Introduce a **thin parameter layer** (`ampere.core.parameter`) so that
   parameters are first-class objects with names, priors, and optional
   constraints.
2. Establish `ampere.core` as a stable package of base abstractions that all
   other sub-packages depend on.
3. **Extend the DI pattern into `Data` classes** so that data-generating
   process, noise model, and likelihood are independently replaceable
   components (see [issue #73][i73]).
4. Re-organise the inference layer so the existing numpy samplers keep working
   while gradient-based / numpyro-based backends can be added later *without*
   breaking the public API.
5. Leave a clean extension point for hierarchical inference.

---

## Design Decisions

| Topic | Decision |
|---|---|
| Parameter management | `ampere.core.parameter.ParameterSet`; parameters are named with scipy priors |
| DI pattern for Data | Extend DI to `Data`: inject `data_generating_process`, `noise_model`, `likelihood` (issue #73) |
| `Data.lnlike` signature | Change from `lnlike(synWave, synFlux)` to `lnlike(theta, result)`, also accepting a param dict |
| Getters / setters | Use `@property` throughout for all important instance attributes (issue #68) |
| True vectorisation | Replace loop-based `lnprob_vector` with `np.vectorize` / `map()` where appropriate (issue #67) |
| JAX | Optional; only required for JAX backends; numpy path must remain functional |
| Gradient-based inference | Not in scope for v1 but must not be architecturally blocked |
| Hierarchical inference | Not in scope; extension point (`hyperpriors` slot) must be kept clean |
| numpyro | Future backend; model API must be numpyro-friendly (named params, no raw slicing) |
| Existing backends | emcee, dynesty, zeus, SBI remain; refactored to use new parameter API |
| New backends (v1) | Snowline (variational Bayes, issue #41); L-M pre-optimiser via Scipy (issue #40) |
| `CompositeModel` | Complete arithmetic composition (blocked by parameter layer) |
| Censored data | Upper/lower limits supported in all likelihood classes (issue #11) |
| Covariance optimisation | Expose hooks for efficient covariance structures (issue #29); full implementation deferred |
| Public API compatibility | Best-effort; breaking changes allowed in `ampere.infer` internals |

---

## Target Package Structure

```
ampere/
├── core/                        # NEW – stable base abstractions
│   ├── __init__.py
│   ├── parameter.py             # ParameterSet / Parameter (Issue 2)
│   ├── model.py                 # Abstract Model base class (Issue 3)
│   └── data.py                  # Abstract Data base class (Issue 4)
├── models/                      # Physical models (refactored in Issue 5)
│   ├── __init__.py
│   ├── models.py
│   ├── blackbodies.py
│   ├── powerlaws.py
│   └── …
├── data/                        # Data containers (Issue 4)
│   ├── __init__.py
│   ├── data.py
│   ├── photometry.py
│   └── spectrum.py
├── infer/                       # Inference backends (Issue 6)
│   ├── __init__.py
│   ├── basesearch.py
│   ├── emceesearch.py
│   ├── dynestysearch.py
│   ├── zeussearch.py
│   ├── nestedsearch.py
│   ├── snowlinesearch.py        # NEW – Issue 7a
│   ├── scipyminsearch.py        # NEW – Issue 7b
│   └── sbi.py
├── backends/                    # NEW – optional heavy backends (Issue 7)
│   ├── __init__.py
│   ├── legacy.py                # Re-exports of numpy-based samplers
│   └── jax/
│       └── __init__.py          # JAX/numpyro stub (future)
└── utils/
    └── …
```

---

## Issues / Work Items

### Issue 0 – Immediate bug fixes (pre-work)
**Branch:** `fix/immediate-bugs`

Fix the blockers that prevent the existing codebase from running cleanly, as a
prerequisite for all other work:

- **Issue #77** – Update `ampere/models/blackbodies.py` to use the current
  astropy `BlackBody` model API:
  ```python
  from astropy.modeling.models import BlackBody
  bb = BlackBody(temperature=t * u.K)
  flux = bb(freq)
  ```
- **Issue #76** – Fix syntax errors in `CCMExtinctionLaw` module.
- **Issue #75** – Replace dependency on the abandoned `extinction` package;
  migrate to `dust_extinction` or inline the required curves.
- **Issue #74** – Fix import error in `ampere/models/extinctionModels.py`.

No architectural changes in this issue.


### Issue 1 – Project setup, CI, and dependencies
**Branch:** `feature/project-setup`
**Depends on:** Issue 0

- Add optional extras to `pyproject.toml`: `jax = ["jax", "numpyro"]`,
  `all = [<everything>]`.
- Set up `tests/` directory with a minimal `pytest` suite.
- Add GitHub Actions CI (`.github/workflows/ci.yml`): lint + tests on push/PR.
- Update `[tool.setuptools] packages` to include `ampere.core` and
  `ampere.backends`.
- Pin minimum Python to 3.9 (3.7/3.8 are EOL).


### Issue 2 – `ampere.core.parameter`: ParameterSet
**Branch:** `feature/core-parameter`
**Depends on:** Issue 1

Implement `ampere/core/parameter.py` containing:

```python
class Parameter:
    name: str
    prior: Any          # scipy.stats frozen rv, or numpyro distribution
    constraint: Any     # optional bijector slot (reserved for JAX path)
    value: float | None
    hyperpriors: ParameterSet | None   # extension point for hierarchical inference

class ParameterSet:
    def pack(self, values: dict) -> np.ndarray: ...
    def unpack(self, theta: np.ndarray) -> dict: ...
    def lnprior(self, theta: np.ndarray) -> float: ...
    def prior_transform(self, u: np.ndarray) -> np.ndarray: ...
    def as_paramax(self): ...   # requires jax; raises ImportError otherwise
    @staticmethod
    def merge(sets, prefixes) -> "ParameterSet": ...  # for CompositeModel
```

Key requirements:
- scipy.stats priors are supported out-of-the-box; numpyro distributions are an
  optional extra.
- `as_paramax()` raises `ImportError` when JAX is absent; otherwise returns a
  Paramax pytree for use with optax/blackjax/numpyro.
- `hyperpriors` slot exists but defaults to `None`; needed for Issue 8.


### Issue 3 – `ampere.core.model`: abstract Model base class
**Branch:** `feature/core-model`
**Depends on:** Issue 2

Split `ampere/models/models.py::Model` into two layers:

- `ampere.core.model.Model` – pure abstract base; interface only.
- `ampere.models.models.Model` – inherits from core; adds arithmetic composition
  and delegates `lnprior`/`prior_transform` to a `ParameterSet`.

New `Model.__call__` signature accepts **keyword arguments** (unpacked from
theta) rather than positional `*args`, making models introspectable and
numpyro-compatible.

Add a `LegacyModel` mixin that accepts positional `*args` for backward
compatibility.

Use `@property` for all public attributes (e.g. `parameters`, `npars`)
so that derived classes can validate or process values on assignment
([issue #68][i68]).

[i68]: https://github.com/ICSM/ampere/issues/68

Complete `CompositeModel` (currently raises `NotImplementedError`):
- Merge parameter sets with name-scoping (`component1.param`, `component2.param`).
- Route keyword arguments correctly in `__call__`.


### Issue 4 – `ampere.core.data`: DI-pattern Data base class
**Branch:** `feature/core-data`
**Depends on:** Issue 2
**Ref:** [issue #73](https://github.com/ICSM/ampere/issues/73)

Ampere's `Data` class should follow a full **dependency-injection** pattern,
with three independently replaceable components injected at construction time:

1. **`data_generating_process`** – converts a `ModelResults` object into a
   noise-free observable in the same space as the measurements (e.g. integrates
   a spectrum through filter curves for photometry).
2. **`noise_model`** – characterises the covariance structure of the data
   (diagonal Gaussian, correlated noise, student-t, …).  Its nuisance
   parameters are declared as a `ParameterSet`.
3. **`likelihood`** – given the noise-free prediction and the noise model,
   returns the scalar log-likelihood.

The minimum interface for the refactored `Data` base class:

```python
class Data:
    parameters: ParameterSet          # nuisance parameters of the noise model
    data_generating_process: callable # (result: ModelResults) -> prediction
    noise_model: object               # declares covariance; owns nuisance params
    likelihood: callable              # (observed, predicted, noise_model) -> float

    def lnlike(self, theta, result: ModelResults) -> float: ...
    # theta may be np.ndarray (legacy) or dict (new API)
```

Additional changes:
- Change `lnlike` signature from the current `lnlike(self, synWave, synFlux)`
  to `lnlike(self, theta, result)`.  The `theta` slice passed to a `Data`
  object covers only *that object's* nuisance parameters.
- Declare each data object's nuisance parameters as a `ParameterSet`.
- Add `@property` getters/setters for all public attributes ([issue #68][i68]).
- Add first-class support for **upper and lower limits** (censored data) in the
  likelihood component ([issue #11](https://github.com/ICSM/ampere/issues/11)).
- Expose a covariance-structure hook so that future work can plug in efficient
  sparse/Toeplitz solvers for large datasets ([issue #29][i29]).
- No changes to `Photometry` or `Spectrum` computation logic in v1; only the
  architecture changes.

[i29]: https://github.com/ICSM/ampere/issues/29


### Issue 5 – Refactor physical models to use `ParameterSet`
**Branch:** `feature/models-parameter-refactor`
**Depends on:** Issues 2, 3

For each concrete model, replace manual `npars` counting and
`lnprior`/`prior_transform` with a `ParameterSet` declaration, and change
`__call__` to accept keyword arguments.

Priority order:
1. `BlackBody` / `ModifiedBlackBody`
2. `PowerLaw`
3. `StarPlusScreen`
4. Dusty / RT wrappers (more complex; can be a sub-issue)

Add `LegacyModel` shim so existing user scripts keep working.


### Issue 6 – Refactor inference backends to use `ParameterSet`
**Branch:** `feature/infer-refactor`
**Depends on:** Issues 3, 4

Update `BaseSearch` and all subclasses so that:
- Parameter counting and slicing is driven by `model.parameters` and
  `data.parameters`, not manual `npars` attributes.
- `lnprior`, `lnlike`, `prior_transform` use the new dict-based API internally.
- The external calling convention of `EmceeSearch`, `DynestySearch`,
  `ZeusSearch` is **unchanged** for users.

Replace:
```python
theta[:self.nparsMod]                                    # positional slice (old)
model.parameters.unpack(theta[:model.parameters.size])  # named dict (new)
```

Also in this issue ([issue #67](https://github.com/ICSM/ampere/issues/67)):
- Replace the loop-based `lnprob_vector` / `lnprior_vector` / `lnlike_vector`
  in `BaseSearch` with genuinely vectorised implementations using
  `map()` / `np.vectorize` / `jax.vmap` (optional).  The current
  implementation simply iterates over samples in Python, defeating the purpose
  of the `vectorize=True` flag passed to emcee.


### Issue 7 – New inference backends
**Branch:** `feature/new-backends`
**Depends on:** Issue 6

Add two new backends that complement the existing emcee / dynesty / zeus set:

#### 7a – Snowline (variational Bayes)
[issue #41](https://github.com/ICSM/ampere/issues/41)
- Implement `SnowlineSearch` inheriting from `BaseSearch`.
- Snowline is a near drop-in replacement for dynesty; use `prior_transform`
  from `ParameterSet`.
- Note: Snowline is optimised for N < 10 parameters; document this limitation.

#### 7b – Scipy L-M pre-optimiser
[issue #40](https://github.com/ICSM/ampere/issues/40)
- Implement `ScipyMinSearch` (or extend existing `ScipyMinMixin`).
- Wrap `scipy.optimize.minimize` (L-M / Nelder-Mead) with the ampere `Model`
  and `Data` interface.
- Primarily useful as a **pre-optimiser** to seed MCMC walkers near the MAP
  solution.


### Issue 8 – `ampere.backends`: optional heavy backends
**Branch:** `feature/backends-stub`
**Depends on:** Issue 6

The `ampere/backends/` package (already scaffolded) provides:
- `legacy.py` – thin re-exports of `EmceeSearch`, `DynestySearch`, `ZeusSearch`
  under the new stable import path.
- `jax/__init__.py` – stub with `NotImplementedError` and a description of the
  intended numpyro integration pattern.

No new sampler logic in this issue; that lives in Issue 7.


### Issue 9 – Hierarchical inference: extension point
**Branch:** `feature/hierarchical-stub`
**Depends on:** Issues 2, 3

**Not a full implementation.** Ensure the architecture does not block future
hierarchical inference:

- `ParameterSet.hyperpriors` slot must exist (defaults to `None`).
- `Model` must expose a `population_model()` hook (returns `None` by default)
  that numpyro can later use to build a plate model.
- Add notes in `ampere/backends/jax/__init__.py` about the intended numpyro
  integration pattern.


### Issue 10 – Documentation and examples
**Branch:** `feature/docs-v1`
**Depends on:** Issues 1–7

- Update `README.md` with new package layout and installation instructions.
- Add `examples/01_basic_fitting.ipynb` using the new API.
- Add `docs/migration.md` explaining how to convert existing models and data
  classes to the DI / `ParameterSet` API.
- Update docstrings throughout `ampere.core`.


### Issue 11 – Test suite
**Branch:** `feature/tests`
**Depends on:** Issues 2–7 (developed alongside each)

Minimum coverage:
- `tests/test_parameter.py` – pack/unpack/lnprior/prior_transform.
- `tests/test_model.py` – `BlackBody` call, `CompositeModel` arithmetic.
- `tests/test_data.py` – `Photometry.lnlike`, `Spectrum.lnlike`; censored data
  upper/lower limits.
- `tests/test_infer.py` – smoke test: construct `EmceeSearch`, run 10 steps.
- `tests/test_infer_snowline.py` – smoke test: `SnowlineSearch`, run 1 iteration.

---

## Recommended Implementation Order

```
Issue 0  (bug fixes)
    └─> Issue 1  (setup)
            └─> Issue 2  (ParameterSet)          ← most critical; unblocks everything
                    ├─> Issue 3  (Model base)
                    │       ├─> Issue 5  (physical models)
                    │       └─> Issue 6  (inference refactor)
                    │               ├─> Issue 7  (new backends: Snowline, L-M)
                    │               └─> Issue 8  (backends stub: legacy + JAX)
                    ├─> Issue 4  (Data DI base)
                    │       └─> Issue 6  (inference refactor)
                    └─> Issue 9  (hierarchical stub)
Issue 10 (docs)    – after Issues 1–7
Issue 11 (tests)   – alongside Issues 2–7
```

---

## Backward Compatibility Strategy

- All existing user-facing classes (`EmceeSearch`, `DynestySearch`, `BlackBody`,
  etc.) must continue to work after each issue is merged.
- `LegacyModel` mixin is the main shim for models still using positional `*args`.
- The `Data.lnlike(synWave, synFlux)` → `Data.lnlike(theta, result)` signature
  change is the most disruptive; a shim wrapper will be provided in the
  `Data` base class to accept both call signatures during the transition.
- Users are encouraged (not forced) to migrate to the `ParameterSet` API.
- The SBI backend uses a different calling convention; a compatibility wrapper
  will be assessed during Issue 6.

---

## Open Questions (deferred)

- Should `ParameterSet` support correlated priors (e.g. multivariate Gaussian)?
  Deferred to the hierarchical inference work (Issue 9).
- Efficient covariance matrix inversion for large spectra / megapixel images
  ([issue #29](https://github.com/ICSM/ampere/issues/29)) – hooks exposed in
  Issue 4; full implementation is a separate milestone.
- Image / interferometry / IFU data classes are out of scope for this plan.
