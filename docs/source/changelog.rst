Changelog
=========

This page is maintained by hand, one section per release. It says what you
can do with each release, not which pull requests made it.

Unreleased
----------

* **Read OIFITS files.** :func:`ampere.interferometry.read_oifits` turns an
  OIFITS file's squared visibilities, closure phases and complex visibilities
  into the interferometric containers, ready for a fit, and the new
  :mod:`ampere.interferometry` is the one place to import interferometry
  from; :class:`~ampere.backends.reference.SquaredAmplitude` fits squared
  visibilities as they were measured. See :doc:`interferometry` §10.
* **Derived parameters.** ``Parameter(name, Derived("mu + sigma * z"))``
  declares a parameter computed from others — no sampler dimension, no prior
  term, handed to the model like any other value and stored in the
  ``posterior`` beside the sampled ones on every engine (named in the run's
  ``ampere_derived`` attribute; provenance schema 10). It makes the
  non-centred population declarable (a member no model declares is internal
  to the population, see :doc:`population`) and gives
  :func:`~ampere.core.shrinkage_horseshoe` Piironen & Vehtari's slab as
  ``tail="slab"``. A population member must now be declared by every
  component it is routed to; one that is not is refused when the problem is
  composed, rather than failing at the first model evaluation.

* **Install from PyPI throughout the documentation.** The overview, the SBI
  guide and the torch and jax backend pages say ``pip install
  ampere-astro[...]`` rather than an install from a clone; the development
  install stays on :doc:`install`.
* **Training sets keep long failure messages and extra-coordinate dtypes.**
  Appending a batch whose failure message (or context record) is longer than
  any in the file no longer truncates it on write: every string variable of a
  training set is stored variable-length. Integer and boolean extra
  coordinates are covered by a test through write, append and read-back.
* **Migrating note: the legacy star scripts are** :math:`4\pi` **too faint.**
  :doc:`migrating` now says that the legacy ``phoenixstar.py`` and
  ``star_disc.py``/``QuickSED`` divide by :math:`4\pi d^2` twice; the twins
  use the correct law.
* **For developers.** A ``heavy`` pytest marker keeps the optimiser rows that
  need a 40-second emcee run out of ``pixi run test-fast``; ``test-all`` and CI
  still run them. ``SyntheticPhotometry.from_library`` no longer depends on a
  legacy module.
* **The scipy optimiser no longer stops at a corner of the prior box.**
  ``optimise(method="scipy")`` now moves a parameter with a box or half-line
  prior in its own normalised value, kept inside the box, rather than
  through a saturating sigmoid. On the NGC6302 twin's sixteen
  parameters this takes the MAP from eleven coordinates at a bound to one.
  :func:`ampere.inference.saturated_bounds` names a converged parameter
  that still sits at a bound, and the route warns with
  :class:`~ampere.inference.BoundSaturationWarning`. The
  :class:`~ampere.results.Optimum` records the coordinates its minimiser
  moved in (``coordinates``). :func:`~ampere.inference.warm_start_gp`
  accepts a GP hyperparameter that a tie has renamed. See
  :doc:`optimisers`.
* **Reading the diagnostics.** A new page, :doc:`reading_the_diagnostics`,
  says what the residual, GP-localisation, posterior-predictive and
  anomaly-score plots show and what they do not, opening with a power law
  fitted to an IRS spectrum (the GP carrying the silicate emission), and
  closing with how to choose a kernel's priors against the model's own
  scales. Its figures are the first to be built with the documentation, by
  matplotlib's plot directive from ``docs/source/plots/``, at a reduced
  budget.

1.0.0b1 — the first beta of ampere v2
-------------------------------------

*Released from the* ``v1.0.0b1`` *tag on 2026-10-05.*

What this release is
~~~~~~~~~~~~~~~~~~~~

ampere v2 is a ground-up redesign of the package, and this beta is its first
release and the first ampere on PyPI. The idea is the same as before — fit a
physical model to heterogeneous astronomical data under a likelihood that is
robust to the model being imperfect — but the package around that idea is
new: you *declare* a fit from named parts (parameters with priors, a model,
an instrument, a dataset with its likelihood, an engine) instead of
subclassing a framework, and the same declaration runs unchanged on numpy,
torch or jax. The old code is not gone: it lives in ``ampere.legacy``, still
answers to its old import names, and is kept indefinitely (see
:doc:`legacy`). The two halves do not interoperate; :doc:`migrating` is the
bridge, with every legacy name mapped and the same fit shown side by side.
The final release of this line will be ``1.0.0``.

Installing it: ``pip install ampere-astro`` (the name ``ampere`` on PyPI
belongs to an unrelated package); the import name is still ``ampere``. The
base install is a complete fitting environment — the reference backend, the
flexible likelihood with its exact O(N) solver, the gradient-free engines and
the results tools. Extras add a backend or an engine:
``pip install "ampere-astro[torch]"``, ``"[jax]"``, ``"[sbi]"``,
``"[nautilus]"``, ``"[ultranest]"``, ``"[blackjax]"``, ``"[zeus]"``,
``"[extinction]"`` or ``"[all]"``. Python 3.12 to 3.14. See :doc:`install`.

What you can fit
~~~~~~~~~~~~~~~~

* **Spectra and SEDs**, singly or combined, including photometry and spectra
  in one fit with the calibration uncertainty of each instrument as a
  parameter (:doc:`sed_composition`, :doc:`photometry_spectra`). Synthetic
  photometry comes from a bundled filter library or one of your own.
* **Images**, with PSF convolution and noise correlated across the grid
  (:doc:`image`).
* **Interferometric visibilities and closure phases**, with Fourier sampling
  of an image model, bandwidth and time smearing, a von Mises likelihood for
  the phases and a complex Gaussian process for correlated visibility
  residuals (:doc:`interferometry`).
* **Astrometric time series** (:doc:`astrometry`), and other time series
  through the same container.
* **Several objects at once**: a population fitted hierarchically in one
  joint run, or built afterwards by reweighting archived single-object fits
  (:doc:`population`).
* **Anything with a simulator and no likelihood** — a radiative-transfer code
  behind a Python call, say — by simulation-based inference (below).

Data arrive as plain arrays in containers (:class:`~ampere.core.Spectrum`,
:class:`~ampere.core.Image`, :class:`~ampere.core.VisibilitySet`,
:class:`~ampere.core.TimeSeries`) indexed by their coordinates, with masks
for missing values and no assumption of a regular grid. There are no file
readers yet (see the limitations below).

The flexible likelihood
~~~~~~~~~~~~~~~~~~~~~~~

The reason to use ampere. A Gaussian process over the residuals absorbs the
structure a wrong or incomplete model leaves behind and is marginalised
while the physical parameters are fitted, so the posterior on those
parameters stays honest. On a deliberately misspecified 20 000-point
spectrum, a chi-square fit lands 113 posterior standard deviations from the
truth; the flexible likelihood lands within 0.6 (:doc:`m2_misspecification`;
:doc:`concept` explains why). Attaching it to a dataset is one line —
:class:`~ampere.core.GaussianProcessNoise` in place of
:class:`~ampere.core.IndependentNoise` — and from there:

* **Kernels** (:doc:`kernels`): seven families, sums and products, spectral
  mixtures, a damped oscillator for periodic residuals, a selector for which
  axis of multi-dimensional data a kernel acts on, and a registry for a
  kernel of your own. Non-stationary noise through a warped kernel, and a
  shrinkage prior that switches off noise components a fit does not need.
* **Solvers** (:doc:`solvers`): the default Matérn-3/2 kernel is solved
  exactly in O(N) by the quasiseparable solver, with the dense solver beside
  it; for data too large or too high-dimensional for either there are
  reduced-rank, Vecchia, inducing-point and structured-grid solvers, each
  documented with its error envelope.
* **Joint noise** across several channels, for astrometric and similar data.
* **Other likelihood families** where a Gaussian is wrong: Student-t,
  Cauchy, Poisson, Rice, von Mises and complex Gaussian, all with the same
  noise models.

Models and instruments
~~~~~~~~~~~~~~~~~~~~~~

A model is a class with named parameters and an ``evaluate`` method that
returns a container; written once, it runs on every backend. Priors are any
scipy distribution, with ties, fixed values and conditional priors
(:doc:`arbitrary_priors`, :doc:`conditional_priors`). Any
**astropy.modeling** model, compound models included, becomes an ampere model
with :func:`~ampere.core.from_astropy` (:doc:`astropy`). An **instrument** is
a chain of transformations — resampling, convolution, synthetic photometry,
a calibration scale — applied to the model's prediction before it meets the
data, so one physical model serves several instruments.

Inference and results
~~~~~~~~~~~~~~~~~~~~~

Nine engines share one interface, so changing sampler is changing one name:
:class:`~ampere.inference.EmceeEngine`, :class:`~ampere.inference.ZeusEngine`
and :class:`~ampere.inference.DynestyEngine` on any backend;
:class:`~ampere.inference.NautilusEngine` and
:class:`~ampere.inference.UltranestEngine` for nested sampling with an
evidence; :class:`~ampere.inference.NUTSEngine` and
:class:`~ampere.inference.VIEngine` (mean-field, full-rank, Laplace and
normalising-flow guides) with gradients on torch and jax;
:class:`~ampere.inference.BlackjaxEngine` (MCLMC, Pathfinder) on jax.
:func:`~ampere.inference.optimise` gives a fast point estimate or Laplace
approximation and :func:`~ampere.inference.warm_start_gp` the GP's
hyperparameters in milliseconds; any engine can start from either
(:doc:`optimisers`).

:class:`~ampere.inference.SBIEngine` fits a model with no likelihood by
neural posterior, likelihood or ratio estimation, or truncated marginal
ratio estimation; simulations run in a process pool with timeouts and crash
capture, trained networks are cached and reused, one network can be
amortised over differently-sampled observations, and calibration checks
(simulation-based calibration, TARP, coverage) say whether to trust the
result (:doc:`sbi`, :doc:`wstat_comparison`).

Every run returns an ArviZ ``DataTree``: the posterior, the sample
statistics, the posterior predictive and residuals, and provenance — the
problem's hash, the seed, the library versions — so a result can be saved
(:func:`~ampere.results.to_netcdf`), reloaded, compared and reproduced.
:mod:`ampere.results` plots it: corner, trace, posterior predictive,
residuals, and a view of where the GP absorbed structure the model missed.

What changed in the interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you have code against the old ampere, this is what to expect.

* **The name on PyPI** is ``ampere-astro``; ``import ampere`` is unchanged.
* **The old modules moved**: ``ampere.data``, ``ampere.models``,
  ``ampere.infer``, ``ampere.utils`` and ``ampere.logger`` are now
  ``ampere.legacy.data`` and so on. The old names still import, as aliases
  of the same modules, with no warning and no removal date, so old scripts
  keep running (:doc:`legacy`).
* **A fit is declared, not subclassed.** Where the old code asked for a
  ``Model`` subclass with ``__call__``, ``lnprior`` and ``prior_transform``,
  a ``Spectrum`` data object carrying its own noise parameters, and an
  ``EmceeSearch`` holding both, v2 separates them:

  .. code-block:: python

      model = ASimpleModel(wavelength)                  # parameters and priors declared inside
      observed = Spectrum(wave, flux, uncertainty=uncertainty)
      dataset = Dataset(observed, likelihood=Likelihood(GaussianFamily(), GaussianProcessNoise(...)))
      problem = FittingProblem(model, [dataset], seed=1)
      run = EmceeEngine(problem, walkers=32).run(steps=2000, burn_in=1000)

  Priors attach to named :class:`~ampere.core.Parameter` declarations; the
  noise model is a named choice on the dataset; the engine is bound to the
  problem rather than mixed into the model. The same fit is shown in full,
  old and new, in :doc:`migrating`.
* **Results are a** ``DataTree``, not the sampler's raw chain and a
  ``postProcess`` call; plots are functions in :mod:`ampere.results` that
  take the tree.
* **Filters and photometry** are the reference backend's
  :class:`~ampere.backends.reference.SyntheticPhotometry` over a bundled
  pyphot library; the legacy filter-building helpers are carried unchanged
  under ``ampere.legacy.utils``.
* **Backends are explicit.** The base install is the numpy reference
  backend; torch and jax are extras, chosen when a model or problem is
  built, and a device (a GPU) is asked for by name and never detected.

The complete name-by-name table is :ref:`migrating-every-name`.

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
