Ampere v2: the architecture
===========================

Ampere v2 is one backend-neutral core plus a small ladder of backends that
implement it. This page is the map: what the pieces are, which environment
each of them needs, how a problem is put together, and what comes out of a
run. :doc:`api` is the reference for every name mentioned here.

The state of play, as of milestone **M2** and the close of Phase 3: the
contracts are frozen and implemented, three backends ship, six engines
ship — simulation-based inference (:class:`~ampere.inference.SBIEngine`)
landed alongside the five likelihood-based engines — and the flexible
likelihood has been measured against a deliberately misspecified problem at
three data sizes and on all three backends. :doc:`m2_misspecification` is
that measurement; :doc:`sbi` is the SBI tutorial.

.. note::

   Everything on this page is the **v2** API — ``ampere.core``,
   ``ampere.backends``, ``ampere.inference``, ``ampere.results``. Legacy
   ampere (``ampere.data``, ``ampere.models``, ``ampere.infer``) still works
   and is frozen; the two do not interoperate. See :doc:`ampere`.

The capability ladder
---------------------

The axis that actually matters is **differentiability**, not which array
library you like. So the backends are not peers; they are rungs, and each
rung can do everything the ones below it can:

.. list-table::
   :header-rows: 1
   :widths: 12 34 54

   * - Rung
     - What runs here
     - Inference available
   * - 0 — black box
     - Any callable ampere cannot see inside: an external radiative-transfer
       code, a legacy numpy model, a wrapped binary
     - Gradient-free sampling (emcee, zeus, dynesty), plus simulation-based
       inference (:class:`~ampere.inference.SBIEngine`) where ``log_prob``
       cannot be evaluated at all.
   * - 1 — reference
     - :mod:`ampere.backends.reference`, pure numpy/scipy
     - The same, plus the exact O(N) GP likelihood as a correctness anchor
   * - 2 — native
     - :mod:`ampere.backends.torch`, :mod:`ampere.backends.jax`
     - Everything above, plus NUTS, variational inference, batching and GPU
       placement

A model does not choose a rung. Its **capability flags** — ``DIFFERENTIABLE``,
``BATCHABLE``, ``DEVICE`` and ``BACKEND`` — are facts about how it was
written, and :func:`ampere.core.declared_capabilities` aggregates them over
every part of a composed problem so an engine can self-select. An engine
asked for something the problem cannot do refuses by name, at construction,
before any sampling starts; it never silently downgrades and it never fails
deep inside an autodiff error.

Three environments
------------------

.. list-table::
   :header-rows: 1
   :widths: 22 30 48

   * - Install
     - Environment
     - What you get

   * - ``pip install -e .``
     - ``pixi run -e dev …``
     - The reference backend, the whole of :mod:`ampere.core` — including the
       flexible GP likelihood and its **exact O(N)** solver — the
       gradient-free engines, and :mod:`ampere.results`. This is a complete
       fitting environment on its own, not a stub: celerite2, arviz and
       h5netcdf are base dependencies, not extras.

   * - ``pip install -e ".[torch]"``
     - ``pixi run -e torch …``
     - Adds :mod:`ampere.backends.torch` and pyro, and therefore
       :class:`~ampere.inference.NUTSEngine`,
       :class:`~ampere.inference.VIEngine`, ``vmap`` batching and per-instance
       device placement.

   * - ``pip install -e ".[jax]"``
     - ``pixi run -e jax …``
     - Adds :mod:`ampere.backends.jax`, numpyro and equinox, with the same
       consequences. The fastest path in the repository: a value **and** a
       gradient of a 20 000-point GP likelihood in 4.7 ms.

No environment has both torch and jax; the two are separately installable by
design, and nothing in ampere requires them together. ``ampere.core``,
``ampere.inference`` and ``ampere.results`` import neither, ever — importing
a backend is your explicit opt-in to its dependency, and
``ampere/backends/__init__.py`` imports nothing at all.

See :doc:`install` for the pixi and pip routes in full.

Composing a problem
-------------------

A fit is built out of five kinds of thing, and the same five on every
backend:

**A model** produces a prediction from parameters. It declares its own
parameters (each with a prior) and its own buffers (fixed arrays such as a
wavelength grid), and returns a container — a :class:`~ampere.core.Spectrum`,
:class:`~ampere.core.PhotometricPoints`, an :class:`~ampere.core.Image`, a
:class:`~ampere.core.TimeSeries`, a :class:`~ampere.core.VisibilitySet` or
:class:`~ampere.core.ClosurePhases`, and so on — per named channel. A
container **kind** is three class attributes (axes, layout, whether values
may be complex), extensible out of tree with no change to ``ampere.core``;
:doc:`interferometry` is the worked template for adding one, and
:doc:`astrometry` the second modality built by following it.

**An instrument** is a chain of transformations from what the model produces
to what a particular dataset observed: a calibration scale, a resampling onto
the observed grid, a line-spread-function convolution, synthetic photometry
through filter curves. A dataset with no instrument gets an implicit one that
only does the binding.

**A likelihood** is a family (the sampling distribution) plus a noise model
(what the covariance is). :class:`~ampere.core.IndependentNoise` is the
ordinary chi-square; :class:`~ampere.core.GaussianProcessNoise` is ampere's
flexible likelihood, a GP over the residuals that absorbs structure the model
cannot explain instead of letting it bias the physical parameters.

**A dataset** binds observed data to an instrument and a likelihood.

**A fitting problem** binds one or more models to a collection of datasets,
carries the seed, and exposes the surface every engine consumes: ``log_prob``,
the ``log_likelihood``/``log_prior`` split, ``prior_transform``, ``simulate``
and the capability flags.

The snippet below composes one model with one instrument and one dataset;
:doc:`sed_composition` is the same five nouns with **two** instruments
bound to one model channel — a spectrum and a photometric catalogue — which
is the ordinary shape of a real SED fit and the case where a dataset's label
stops being optional.

On the reference backend
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    import astropy.units as u
    import numpy as np
    import scipy.stats as st

    from ampere.backends.reference import BlackBody
    from ampere.core import (
        Dataset, FittingProblem, GaussianFamily, GaussianProcessNoise,
        Likelihood, Matern32, QuasisepGP, Spectrum,
    )
    from ampere.inference import EmceeEngine

    wavelength = np.linspace(5.0, 30.0, 200)

    model = BlackBody(
        wavelength,
        temperature=st.loguniform(100.0, 3000.0),
        scale=st.loguniform(0.1, 10.0),
    )

    observed = Spectrum(
        wavelength * u.um, flux * u.Jy, uncertainty=uncertainty * u.Jy
    )

    noise = GaussianProcessNoise(
        Matern32(amplitude=st.halfnorm(scale=0.1),
                 length_scale=st.loguniform(0.5, 10.0)),
        QuasisepGP(),          # exact, and O(N)
    )

    dataset = Dataset(observed, likelihood=Likelihood(GaussianFamily(), noise))
    problem = FittingProblem(model, [dataset], seed=20260908)

    run = EmceeEngine(problem, walkers=32).run(steps=2000, burn_in=1000)

Two details worth pointing at. The priors are ordinary
``scipy.stats`` distributions, and the half-normal on the GP amplitude is
load-bearing: its mass piles up at zero, so "this model is not misspecified"
is a conclusion the fit can reach rather than a corner it is pushed away from.
And :class:`~ampere.core.QuasisepGP` computes the *same* likelihood as
:class:`~ampere.core.DenseGP` — the Matérn-3/2 kernel is exactly
quasiseparable — in O(N) rather than O(N³), which is the difference between
200 points and 20 000.

On torch
~~~~~~~~

The same shape, with the parts taken from the backend instead:

.. code-block:: python

    from ampere.backends.torch import (
        BlackBody, GaussianProcessNoise, Matern32, QuasisepGP,
    )
    from ampere.core import Dataset, FittingProblem, GaussianFamily, Likelihood
    from ampere.inference import NUTSEngine

    model = BlackBody(wavelength, temperature=..., scale=...)
    noise = GaussianProcessNoise(Matern32(...), QuasisepGP())
    dataset = Dataset(observed, likelihood=Likelihood(GaussianFamily(), noise))
    problem = FittingProblem(model, [dataset], seed=20260908)

    problem.capabilities        # differentiable=True, backend='torch'
    run = NUTSEngine(problem).run(draws=1000, warmup=1000)

On jax
~~~~~~

Identical, with one thing to do first — jax's ``jax_enable_x64`` is
process-global state that its own documentation addresses to the application
author, so ampere never flips it for you, and every piece of the jax backend
raises at construction if it is off:

.. code-block:: python

    from ampere.backends.jax import configure_x64
    configure_x64()             # idempotent; before any jax work

    from ampere.backends.jax import (
        BlackBody, GaussianProcessNoise, Matern32, QuasisepGP,
    )

Raising rather than warning is deliberate: a warning in a notebook scrolls
away, and what follows a float32 GP solve is a plausible-looking wrong answer
rather than a crash. float64 for all likelihood and GP linear algebra is a
policy on every backend, not a default.

The one-backend rule
~~~~~~~~~~~~~~~~~~~~

**Every part of one problem must come from one backend, and one device.**
:func:`~ampere.core.declared_capabilities` aggregates ``DIFFERENTIABLE`` and
``BATCHABLE`` conjunctively — one black-box step is enough to stop a gradient
— but ``BACKEND`` and ``DEVICE`` are identities rather than promises, so
disagreement is refused at composition. There is no conservative answer to
"half of this problem is torch and half is jax", and there is no silent
answer either: converting arrays between libraries behind your back is how a
run loses its gradients, and moving them between devices behind your back is
how a run becomes mysteriously slow.

The parts that count are the model, every instrument step, **the noise model,
and the GP solver when a GP is declared**. The last two joined the list in
W2.13, and that is the trap this section exists to name: composing
``ampere.core.IndependentNoise`` or ``ampere.core.DenseGP`` into an otherwise
torch problem used to leave a numpy island in the middle of a differentiable
chain. It is now a backend disagreement, reported by name at composition.

Two things are deliberately *not* in that list. A likelihood **family** is a
declaration of a sampling distribution, evaluated by whichever path the
problem is on, so it has no backend of its own — ``GaussianFamily`` above is
imported from :mod:`ampere.core` on all three backends, and that is correct.
A **kernel** is consumed by its solver, which is the piece that decides
whether the covariance is built in numpy or natively, and the solver is
already on the list.

Realisation, and the engines
----------------------------

:mod:`ampere.inference` is written once against the fitting-problem surface
and **imports no backend**, in any module, at any depth. Six engines ship:

.. list-table::
   :header-rows: 1
   :widths: 20 22 58

   * - Engine
     - Needs
     - What it is for

   * - :class:`~ampere.inference.EmceeEngine`
     - base install
     - Affine-invariant ensemble MCMC. The general-purpose default.
   * - :class:`~ampere.inference.DynestyEngine`
     - base install
     - Nested sampling. Multimodal posteriors, and the only one of the six
       that yields a marginal likelihood — so model comparison.
   * - :class:`~ampere.inference.ZeusEngine`
     - ``zeus`` extra
     - Ensemble slice sampling.
   * - :class:`~ampere.inference.NUTSEngine`
     - ``torch`` or ``jax``
     - The No-U-Turn sampler: numpyro's on a jax problem, pyro's on a torch
       one. Gradient-based, so it scales to many more parameters.
   * - :class:`~ampere.inference.VIEngine`
     - ``torch`` or ``jax``
     - Stochastic variational inference. **Approximate** — the guide family
       is recorded in the run — and the honest use is a first look, or the
       only tractable route when the space is too large for MCMC.
   * - :class:`~ampere.inference.SBIEngine`
     - ``sbi`` extra
     - Simulation-based inference (NPE, NLE, NRE, truncated variants) —
       trains a neural density estimator on simulated ``(theta, x)`` pairs
       instead of consuming ``log_prob``, so it is the route for a rung-0
       black box with no tractable likelihood at all. **Approximate**, with
       calibration and caching built in. :doc:`sbi` is the tutorial.

The first three consume only the neutral surface, which is the architectural
bet of the whole redesign, cashed: one driver, any backend, and a black-box
model behind a thin adapter is not a special case. :class:`~ampere.inference.SBIEngine`
consumes the same surface through :meth:`~ampere.core.dataset.FittingProblem.simulate_many`
rather than ``log_prob``, so it too runs on any rung, black box included.

The last two cannot be, and the reason is precise: ``log_prob`` is not
*traceable*. It runs the chain through containers that coerce with
``numpy.asarray``, and a tensor carrying an autograd graph cannot be put in
one, so the graph is cut at the first container on every backend. The answer
is a **realisation**: the backend's own registered, differentiable, native
form of the same problem. Importing a backend registers its realisation;
:func:`ampere.core.realise` hands it back, and checks it against the contract
path before anyone samples with it. So ``NUTSEngine(problem)`` takes no
density argument, imports no backend, and dispatches on the problem's own
``backend`` flag.

Failures, and your priors
~~~~~~~~~~~~~~~~~~~~~~~~~

A proposal the model or likelihood cannot score becomes ``-inf`` with a
recorded reason rather than an exception, and the aggregate is surfaced at
the end of the run three ways: a warning, an ``ampere_failure_summary``
attribute on the run, and
:attr:`~ampere.inference.engine.Engine.last_failure_summary`. Run non-strict,
read the summary, then re-run with ``FittingProblem(..., strict=True)`` to get
the raise at the offending draw with a traceback.

**The catch set is narrow, on purpose, and the sharp edge is your priors.**
Only :class:`~ampere.core.LikelihoodError` and the types you name in
``simulator_failures=`` are caught; everything else propagates, because a
composition bug turned into ``-inf`` is a fit that runs, converges and is
wrong. A model that refuses a physically meaningless value raises
``ValueError`` — the shipped ``BlackBody`` does, for a non-positive
temperature — and ``ValueError`` is not caught. So a ``scipy.stats.norm``
prior on a positivity-constrained parameter will kill a run mid-flight the
first time the sampler proposes a negative value. The remedy is nearly always
a prior whose support *is* the parameter's support: ``loguniform``,
``lognorm``, ``truncnorm``, ``halfnorm``.

Results, and diagnostics
------------------------

There is one results format, and every engine emits it: an ArviZ
``DataTree``, returned by ``run()`` and built by :func:`ampere.results.emit`.
It carries the posterior, per-draw ``log_prior`` and ``log_likelihood``, the
per-dataset log-likelihood decomposition, the observed data, and provenance
attributes — the parameter-spec hash, the data hashes, the backend-neutral
model identity, the seed, package versions, failure counts. Write it with
:func:`~ampere.results.to_netcdf` and read it back with
:func:`~ampere.results.from_netcdf`; the hashes survive the round trip by
value, which is what makes an archived run identifiable years later.

All plotting is written once against that format —
:func:`~ampere.results.plot_corner`, :func:`~ampere.results.plot_trace`,
:func:`~ampere.results.plot_posterior_predictive`,
:func:`~ampere.results.plot_residuals` — so "sampler plotting parity keeps
regressing" stops being a class of bug.

Two post-fit diagnostics answer the two questions a flexible likelihood
raises, and each is scoped to the fit it makes sense on:

:func:`~ampere.results.residual_whiteness`
    *Did a standard-likelihood fit leave structure behind?* A
    separation-binned, permutation-calibrated whiteness test, with a
    chi-square posterior-predictive *p*-value beside it. Asking it of a GP fit
    would be close to circular — those residuals are whitened by construction.

:func:`~ampere.results.gp_localisation`
    *Where did a flexible fit need its GP?* The conditioned GP mean, converted
    into an :class:`~ampere.core.AnomalyScore` over the observed coordinate —
    which is a map of where your model is failing.

In the M2 study, the second of these was told nothing about a 12 % emission
line injected at 0.86300 µm and put its peak 0.05 nm away.

Where to go next
----------------

* :doc:`m2_misspecification` — the evidence that the flexible likelihood
  works, what it costs, and the size ladder that shows why misspecification
  matters *more* as spectra get larger.
* :doc:`interferometry` — the template for adding a new observable kind,
  proved end to end on interferometric visibilities and closure phases;
  :doc:`astrometry` is the second modality built by following it.
* :doc:`kernels` — the kernel algebra in full: the seven quasiseparable
  families, ``Sum``/``Product``/``SpectralMixture``, the ``axes=`` selector,
  and registering your own term.
* :doc:`astropy` — wrapping an ``astropy.modeling`` model as an ampere
  model, the capability consequence, and the opt-in native route for the
  common cases.
* :doc:`wstat_comparison` — a worked example of registering your own
  likelihood family.
* :doc:`api` — the reference for every name above.
* ``docs/design/`` in the repository — the frozen contract specifications
  (``spec-v1.0``) these APIs implement, if you are extending ampere rather
  than using it.
