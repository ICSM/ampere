Changelog
=========

This page is maintained by hand, one section per release. It says what you
can do with each release, not which pull requests made it.

1.0.0b1 — the first beta of ampere v2
-------------------------------------

*Released from the* ``v1.0.0b1`` *tag.* This is the first release of the
redesigned ampere ("v2") and the first ampere published on PyPI. The major
version marks the redesign: v2 is a new package grown beside the old one, not
an upgrade of it. The final release will be ``1.0.0``.

**Three names.** The distribution on PyPI is ``ampere-astro``
(``pip install ampere-astro``; the name ``ampere`` on PyPI belongs to an
unrelated package); the import name is ``ampere``; the documentation is at
https://ampere.readthedocs.io/. Extras are spelt the same way:
``pip install "ampere-astro[jax]"``. See :doc:`install`.

**The legacy code is kept.** The v1 code that produced ampere's published
science — ``ampere.data``, ``ampere.models``, ``ampere.infer`` — now lives in
``ampere.legacy``, still answers to its old import names, and is kept
indefinitely: it is frozen, receives only critical fixes, and will not be
removed. It implements none of the v2 contracts, and the two halves do not
interoperate. :doc:`legacy` is the policy; :doc:`migrating` maps every legacy
class and call onto its v2 counterpart.

What follows is what the beta lets you do, phase by phase of the redesign.

Phase 0 — a safety net under the old code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The legacy examples are pinned by characterisation tests with fixed seeds, so
the v1 code keeps producing the numbers it always did. ``import ampere`` no
longer fails because of a broken corner of the legacy code: the modules that
did not import (issues #74–#77) are fixed or import lazily, the abandoned
``extinction`` package is replaced by ``dust_extinction`` (the ``extinction``
extra), and pyphot 2 is supported. Python 3.12, 3.13 and 3.14 are supported
and tested.

Phase 1 — the contracts
~~~~~~~~~~~~~~~~~~~~~~~

ampere v2 is built on a small set of frozen, backend-neutral contracts in
:mod:`ampere.core`: parameters and priors (with ties, fixed values and
buffers), result containers indexed by their coordinates with first-class
masks and no assumption of a regular grid, instruments as chains of
transformations, likelihood families and noise models, datasets, and the
:class:`~ampere.core.FittingProblem` that composes them. Every model,
instrument and likelihood you write against them runs on every backend.
:doc:`overview` is the map and :doc:`concept` the idea behind the flexible
likelihood.

Phase 2 — three backends, and the flexible likelihood at scale
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Three backends implement the contracts and are held to one another by a
shared conformance suite: a pure numpy/scipy **reference** backend, which is
the base install and needs no extra; **torch** (``pip install
"ampere-astro[torch]"``); and **jax** (``"ampere-astro[jax]"``). The flexible
likelihood — a Gaussian process over the residuals, marginalised while the
physical parameters are fitted — defaults to a Matérn-3/2 kernel solved
exactly in **O(N)** by a quasiseparable solver (``QuasisepGP``), with the
dense solver beside it. Fits run under emcee and dynesty on any backend, and
under NUTS and variational inference on torch and jax. Every run is an ArviZ
``DataTree`` with provenance — the problem's hash, the library versions, the
seed — that can be saved, reloaded and plotted with :mod:`ampere.results`'s
plots (corner, trace, posterior predictive, residuals, GP localisation).
:doc:`m2_misspecification` is the flagship measurement: on a deliberately
misspecified 20 000-point spectrum, a chi-square fit lands 113 posterior
standard deviations from the truth and the flexible likelihood within 0.6.
:doc:`sed_composition` and :doc:`photometry_spectra` are the first fits to
read.

Phase 3 — simulation-based inference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:class:`~ampere.inference.SBIEngine` fits a model with no likelihood to write
down — a compiled radiative-transfer code behind a Python call, say — by
neural posterior, likelihood or ratio estimation, and by truncated marginal
ratio estimation for a tighter fit (the ``sbi`` extra). Simulations run in
batches under a process pool with timeouts and crash capture; trained
posteriors are cached and reused; runs are reproducible from the problem's
seed; and calibration checks (SBC, TARP, coverage) tell you whether to trust
the result. Embedding networks read any container through one
coordinate–value–mask encoding, so irregular sampling and missing data need
no special handling. :doc:`sbi` is the tutorial and :doc:`wstat_comparison`
a worked calibration study.

Phase 4 — new observables, kernels and astropy models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Interferometric **visibilities and closure phases** fit end to end, with
Fourier sampling of an image model, bandwidth and time smearing, a von Mises
likelihood for phases and a circular complex GP for correlated visibility
residuals (:doc:`interferometry`, which is also the template for adding an
observable of your own). **Astrometric time series** were then added by
following that template (:doc:`astrometry`). The flexible likelihood's
**kernel algebra** grew to seven families, sums, products, spectral mixtures,
a damped oscillator for periodic residuals, an ``axes=`` selector for
multi-axis data and a public registry for your own quasiseparable term
(:doc:`kernels`). Any **astropy.modeling** model, compound models included,
wraps as an ampere model with :func:`~ampere.core.from_astropy`, with an
opt-in differentiable translation for six common models (:doc:`astropy`).

Phase 5 — scale and advanced inference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Approximate GP solvers** for data too large or too high-dimensional for the
exact ones — reduced-rank (``HilbertSpaceGP``, ``EquispacedFourierGP``),
Vecchia, inducing-point and structured-grid solvers — each chosen by
measurement and documented with its error envelope (:doc:`solvers`).
**Images**, with PSF convolution and correlated noise over the grid
(:doc:`image`). **Non-stationary noise** through a warped kernel, and a
shrinkage prior that switches off noise components a fit does not need
(:doc:`kernels`). **Joint noise** over several channels, for astrometric and
similar data. **More engines** behind their own extras: nautilus and
UltraNest nested sampling (``nautilus``, ``ultranest``) and blackjax's MCLMC
and Pathfinder (``blackjax``, with ``jax``), plus more variational guide
families — nine engines in all, written once against the contracts.
**Amortised SBI** over the observation context, so one trained network
serves differently-sampled observations. **Populations**: many objects fitted
hierarchically in one joint fit, or by reweighting archived single-object
fits (:doc:`population`). And a measured optimisation pass that made the
common paths several times faster without changing a number.

Phase 6 — documentation, migration and this release
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Optimisers** for a fast point estimate and a warm start:
:func:`~ampere.inference.optimise` (multi-start scipy minimisation on any
backend, gradient MAP on torch and jax, or a Laplace approximation) and
:func:`~ampere.inference.warm_start_gp` for a GP's hyperparameters in
milliseconds; any sampling engine can start from the result
(:doc:`optimisers`). The **documentation** was rebuilt around v2 and is
published per version on Read the Docs, with every public name documented;
the **migration guide** (:doc:`migrating`) and v2 twins of the legacy
examples' models sit beside the originals, and :doc:`photometry_spectra`
shows calibration uncertainty in a combined fit. SBI training sets accept
non-numeric coordinates, and scoring the stored draws in
:class:`~ampere.inference.SBIEngine` is optional. The package is on PyPI as
``ampere-astro``, built and published by a release workflow; it can be cited
with :func:`ampere.cite` (:doc:`citing`); and the repository has a
contributing guide, a code of conduct and a security policy.

Known limitations at the beta
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **No pre-fit diagnostics module.** ``ampere.diagnostics`` (RHMF screening
  of a dataset before fitting) has not landed: the exploratory trial did not
  meet the maturity gate set for it. Post-fit residual diagnostics are in
  :mod:`ampere.results`.
* **One level of hierarchy.** A parameter carries at most one plate, so
  nested populations (objects within surveys) cannot be declared.
* **Two design items are planned for the next phase, not shipped**: a
  population over a dataset's own parameters (each spectrum's GP amplitude
  drawn from one shared prior), and a ``Derived`` parameter node (a parameter
  that is a function of others). Until then, a per-dataset model component
  that owns the nuisance, or a tie, is the workaround.
* **GPU support is exercised by smoke tests only**, run by hand on a cluster
  at release time, never in continuous integration; the hosted CI is
  CPU-only. The first such run, before this release, found and fixed three
  device-placement faults; one limitation it found stands: **the jax
  ``QuasisepGP`` solver is CPU-only**, because celerite2's jax primitives have
  CPU lowerings alone — on an accelerator it refuses by name, and ``DenseGP``
  or ``HilbertSpaceGP`` is the solver to use there. The torch quasiseparable
  solver has no such limit.
* **No file readers yet.** Observations are built from arrays you read
  yourself; an OIFITS reader for interferometric data is planned after the
  beta.

Deprecations
~~~~~~~~~~~~

* ``ampere.core.regularised_horseshoe`` is a deprecated alias of
  :func:`~ampere.core.shrinkage_horseshoe` (same signature, same return
  value). It emits a ``DeprecationWarning`` and will be removed in ``1.0.0``.
