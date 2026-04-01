# Ampere v1 Implementation Plan

## Background and Motivation

The current ampere codebase provides a working SED/spectrum fitting framework
built around:

- `ampere.models` – physical models (blackbodies, power laws, RT wrappers, …)
- `ampere.data` – data containers with built-in likelihood evaluation
- `ampere.infer` – inference backends (emcee, dynesty, zeus, SBI/PyTorch)

The parameter-management story is entirely manual: each model counts its own
`npars`, slices raw `theta` numpy arrays positionally, and delegates `lnprior` /
`prior_transform` through the same positional slicing. This makes models hard to
compose, hard to introspect, and essentially incompatible with gradient-based /
probabilistic-programming backends (numpyro, blackjax, …).

### Goals for v1

1. Introduce a **thin parameter layer** (`ampere.core.parameter`) so that
   parameters are first-class objects with names, priors, and optional
   constraints.
2. Establish `ampere.core` as a stable package of base abstractions that all
   other sub-packages depend on.
3. Re-organise the inference layer so the existing numpy samplers keep working
   while gradient-based / numpyro-based backends can be added later *without*
   breaking the public API.
4. Leave a clean extension point for hierarchical inference.

---

## Design Decisions

| Topic | Decision |
|---|---|
| Parameter management | `ampere.core.parameter.ParameterSet`; parameters are named with scipy priors |
| JAX | Optional; only required for JAX backends; numpy path must remain functional |
| Gradient-based inference | Not in scope for v1 but must not be architecturally blocked |
| Hierarchical inference | Not in scope; extension point (`hyperpriors` slot) must be kept clean |
| numpyro | Future backend; model API must be numpyro-friendly (named params, no raw slicing) |
| Existing backends | emcee, dynesty, zeus, SBI remain; refactored to use new parameter API |
| `CompositeModel` | Complete arithmetic composition (blocked by parameter layer) |
| Data classes | Minimal changes; `lnlike` accepts both raw arrays and parameter dicts |
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

### Issue 1 – Project setup, CI, and dependencies
**Branch:** `feature/project-setup`

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

Complete `CompositeModel` (currently raises `NotImplementedError`):
- Merge parameter sets with name-scoping (`component1.param`, `component2.param`).
- Route keyword arguments correctly in `__call__`.


### Issue 4 – `ampere.core.data`: abstract Data base class
**Branch:** `feature/core-data`
**Depends on:** Issue 2

- Move abstract `Data` base class to `ampere/core/data.py`.
- Update `lnlike(self, theta, result)` to also accept a parameter **dict**.
- Declare each data object's nuisance parameters as a `ParameterSet`.
- No changes to `Photometry` or `Spectrum` computation logic.


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
theta[:self.nparsMod]                         # positional slice (old)
model.parameters.unpack(theta[:model.parameters.size])  # named dict (new)
```


### Issue 7 – `ampere.backends`: optional backend stubs
**Branch:** `feature/backends-stub`
**Depends on:** Issue 6

Create `ampere/backends/` with:
- `legacy.py` – thin re-exports of `EmceeSearch`, `DynestySearch`, `ZeusSearch`
  under the new stable import path.
- `jax/__init__.py` – stub with `NotImplementedError` and a description of the
  intended numpyro integration pattern.


### Issue 8 – Hierarchical inference: extension point
**Branch:** `feature/hierarchical-stub`
**Depends on:** Issues 2, 3

**Not a full implementation.** Ensure the architecture does not block future
hierarchical inference:

- `ParameterSet.hyperpriors` slot must exist (defaults to `None`).
- `Model` must expose a `population_model()` hook (returns `None` by default)
  that numpyro can later use to build a plate model.
- Add notes in `ampere/backends/jax/__init__.py` about the intended numpyro
  integration pattern.


### Issue 9 – Documentation and examples
**Branch:** `feature/docs-v1`
**Depends on:** Issues 1–6

- Update `README.md` with new package layout and installation instructions.
- Add `examples/01_basic_fitting.ipynb` using the new API.
- Add `docs/migration.md` explaining how to convert existing models.
- Update docstrings throughout `ampere.core`.


### Issue 10 – Test suite
**Branch:** `feature/tests`
**Depends on:** Issues 2–6 (developed alongside each)

Minimum coverage:
- `tests/test_parameter.py` – pack/unpack/lnprior/prior_transform.
- `tests/test_model.py` – `BlackBody` call, `CompositeModel` arithmetic.
- `tests/test_data.py` – `Photometry.lnlike`, `Spectrum.lnlike`.
- `tests/test_infer.py` – smoke test: construct `EmceeSearch`, run 10 steps.

---

## Recommended Implementation Order

```
Issue 1  (setup)
    └─> Issue 2  (ParameterSet)          ← most critical; unblocks everything
            ├─> Issue 3  (Model base)
            │       ├─> Issue 5  (physical models)
            │       └─> Issue 6  (inference refactor)
            │               └─> Issue 7  (backends stub)
            ├─> Issue 4  (Data base)
            │       └─> Issue 6  (inference refactor)
            └─> Issue 8  (hierarchical stub)
Issue 9  (docs)    – after Issues 1–6
Issue 10 (tests)   – alongside Issues 2–6
```

---

## Backward Compatibility Strategy

- All existing user-facing classes (`EmceeSearch`, `DynestySearch`, `BlackBody`,
  etc.) must continue to work after each issue is merged.
- `LegacyModel` mixin is the main shim for models still using positional `*args`.
- Users are encouraged (not forced) to migrate to the `ParameterSet` API.
- The SBI backend uses a different calling convention; a compatibility wrapper
  will be assessed during Issue 6.

---

## Open Questions (deferred)

- Should `ParameterSet` support correlated priors (e.g. multivariate Gaussian)?
  Deferred to the hierarchical inference work.
- Image / interferometry / IFU data classes are out of scope for this plan.
