# ampere

Ampere is a Bayesian fitting environment for astronomers: a tool for modelling
several kinds of astronomical data at once — SEDs, spectra, and more — **even
when the model cannot explain everything in the data**.

That last clause is the point of the package. Deficiencies in a physical model
show up as structure in the residuals, and structure in the residuals is
indistinguishable from correlated noise. So ampere gives each dataset a
flexible likelihood: a Gaussian process over the residuals, marginalised over
while the physical parameters are fitted. Whatever the model cannot explain is
absorbed by the GP instead of biasing the parameters you care about, and no
one has to identify the offending regions by hand.

It works, and the cost is measured. On a deliberately misspecified 20 000-point
spectrum, an ordinary chi-square fit reports a parameter **113 of its own
posterior standard deviations away from the truth** — a precise answer that is
wrong, with nothing in the fit to say so — while the flexible likelihood stays
within 0.6 and keeps the truth inside its 68 % interval. The error of the
standard fit *grows* as √N along the size ladder; the flexible fit's does not,
which is why misspecification matters more, not less, for real spectra with
thousands of pixels. That is only reachable because the Matérn-3/2 kernel is
exactly quasiseparable: ampere solves the GP in **O(N)** rather than O(N³), and
on the jax backend takes a value *and* a gradient of a 20 000-point GP
likelihood in **4.7 ms**. The whole study is
`docs/source/m2_misspecification.rst`, and it is re-run by `pixi run test-m2`.

Ampere is in **alpha**, undergoing a v2 redesign. What is there today is a
backend-neutral core of frozen contracts, three backends that implement it
(pure numpy/scipy, torch and jax), and five inference engines written once
against the contracts and run on any of them. If you are interested, please
get in touch.

## Installation

Ampere needs Python 3.11–3.13 and is not on PyPI yet, so start from a clone.

> **Note:** there is an unrelated project called `ampere` on PyPI (a battery
> modelling package). `pip install ampere` installs *that*, not this.

```bash
git clone https://github.com/ICSM/ampere.git
cd ampere
```

**With [pixi](https://pixi.sh) — the supported route.** It builds the whole
environment, interpreter included, from `pyproject.toml` and the committed
lock file:

```bash
pixi install -e dev
pixi run test-all        # everything the v2 namespaces own
pixi run docs            # build the documentation
```

**With pip**, in a clean virtual environment:

```bash
pip install -e .
```

That is the base install, and it is a complete fitting environment rather than
a stub: the reference backend, the flexible GP likelihood with its exact O(N)
solver, emcee and dynesty, and the ArviZ results format. celerite2, arviz and
h5netcdf are base dependencies, not extras, for exactly that reason.

| Extra | Adds | Unlocks |
|---|---|---|
| *(none)* | numpy, scipy, astropy, matplotlib, spectres, tqdm, pyphot, emcee, dynesty, corner, celerite2, arviz, h5netcdf | the base install; must always work |
| `torch` | torch, pyro-ppl | `ampere.backends.torch`, and NUTS/VI on a torch problem |
| `jax` | jax, numpyro, equinox | `ampere.backends.jax`, and NUTS/VI on a jax problem |
| `zeus` | zeus-mcmc | `ampere.inference.ZeusEngine` |
| `sbi` | torch, sbi | neural posterior estimation (legacy `ampere.infer.sbi`) |
| `extinction` | dust_extinction | legacy `ampere.models.extinctionModels.F99Extinction` |
| `dev` | pytest, ruff, pyrefly, sphinx, … | contributor tooling |
| `all` | every feature extra above | — |

```bash
pip install -e ".[torch]"     # or [jax], [zeus], [all], …
```

The corresponding pixi environments are `pixi run -e torch …`,
`pixi run -e jax …`, `pixi run -e sbi …`. No environment carries both torch
and jax; nothing in ampere needs them together.

## Where to look next

- **Documentation**: `pixi run docs`, then `docs/_build/html/index.html`.
  Start at the architecture overview (`docs/source/overview.rst`), then the
  misspecification study (`docs/source/m2_misspecification.rst`).
- **Runnable examples**: `examples/m2_misspecification` (the study above) and
  `examples/wstat_comparison.py` (registering your own likelihood family).
  Both are covered by tests, so neither can rot unnoticed.
- **Design**: `docs/design/` holds the frozen contract specifications
  (`spec-v1.0`) the v2 API implements — start with `architecture.md`.
- **Contributing**: `DEVELOPMENT_PLAN.md` is the source of truth for where the
  redesign is going; `WORK_ITEMS.md` is the current work, item by item; and
  `docs/development.md` is the onboarding note.

Legacy ampere — `ampere.data`, `ampere.models`, `ampere.infer` — is the v1
code that produced the published science. It still runs and is frozen in
place; it implements none of the v2 contracts and the two halves do not
interoperate.
