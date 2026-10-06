Migrating to the new API
=========================

Ampere has two APIs, and this page is about the gap between them: what maps
to what, what "frozen" means for the old one, and the same small fit shown
both ways, then routes **every public legacy name** to its v2 equivalent (or
says that it is carried), sets the six legacy examples beside their v2 twins,
and maps the post-processing onto :mod:`ampere.results`. The legacy API stays
**linked**, not removed: there is no deprecation timetable. See
:doc:`overview` for how the v2 pieces fit together and :doc:`ampere.legacy`
for the legacy reference.

What "frozen" means
--------------------

**Legacy ampere** — ``ampere.legacy`` — ``ampere.legacy.data``, ``.models``,
``.infer``, ``.utils``, still importable under the old top-level names — is
frozen: it still runs, the characterisation suite exists to keep it running,
and it is not being extended or refactored. It stays because the published
science was produced with it, and because rewriting every existing script
the day v2 lands would be a worse outcome than living with two APIs for a
while.

Three things follow from "frozen":

* **It keeps working.** Nothing in this migration removes it, and nothing in
  v2's roadmap plans to.
* **It implements none of the v2 contracts.** A legacy ``Model`` is not an
  :class:`ampere.core.Model`; a legacy search is not an
  :class:`ampere.inference.Engine`. Neither side's classes substitute for
  the other's.
* **The two APIs do not interoperate.** You cannot hand a legacy
  ``Photometry`` object to a v2 :class:`~ampere.core.Dataset`, or point an
  :class:`~ampere.inference.EmceeEngine` at a legacy model. A fit is built
  entirely on one side or entirely on the other.

The policy — kept indefinitely, never removed, and where it lives — is
:doc:`legacy`.

Concept map
-----------

The two APIs solve the same problems with different shapes. This table lines
up each legacy piece with the v2 piece that plays the equivalent role; it is
a map for orientation, not a claim that either side is a drop-in replacement
for the other.

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Legacy (``ampere.infer``/``ampere.data``/``ampere.models``)
     - Ampere v2
     - Notes
   * - A ``Model`` subclass (``ampere.models``), returning a flux array
     - A :class:`~ampere.core.Model` subclass publishing named
       :class:`~ampere.core.ModelResult` channels
     - v2 models declare their own parameters and buffers, and return typed
       containers per channel rather than a bare array.
   * - ``Photometry``/``Spectrum`` (``ampere.data``)
     - :class:`~ampere.core.Dataset` binding an observed container to an
       :class:`~ampere.core.Instrument` chain
     - The instrument (calibration, resampling, synthetic photometry, LSF
       convolution) is a composable chain in v2 rather than built into the
       data class.
   * - The likelihood's built-in GP switches (correlated-noise strength and
       scale-length parameters carried on ``Spectrum``)
     - :class:`~ampere.core.Likelihood` with
       :class:`~ampere.core.GaussianProcessNoise` and an explicit
       :class:`~ampere.core.Kernel`
     - The flexible likelihood is opt-in and composable in v2: pick a noise
       model and a kernel, rather than accepting the hardcoded
       squared-exponential term.
   * - ``EmceeSearch``
     - :class:`~ampere.inference.EmceeEngine`
     -
   * - ``ZeusSearch``
     - :class:`~ampere.inference.ZeusEngine`
     -
   * - ``DynestyNestedSampler``/``DynestyDynamicNestedSampler``
     - :class:`~ampere.inference.DynestyEngine`
     - One engine covers both nested-sampling modes: ``dynamic=True`` is the
       dynamic one.
   * - *(gradient-free only; no legacy equivalent)*
     - :class:`~ampere.inference.NUTSEngine`,
       :class:`~ampere.inference.VIEngine`
     - Native (torch/jax) samplers with no legacy counterpart — see the
       capability ladder in :doc:`overview`.
   * - ``SBI_SNPE``
     - :class:`~ampere.inference.SBIEngine`
     - Landed (W3.2–W3.15, Phase 3). Simulation-based inference over ``sbi``
       0.27, fitting the same :class:`~ampere.core.FittingProblem` every
       other v2 engine fits by training on simulated pairs from
       :attr:`~ampere.core.FittingProblem.simulate_many` rather than
       consuming ``log_prob``. Beyond a name change: NLE and NRE join NPE,
       plus truncated marginal ratio estimation (``method="tmnre"``); a
       trained posterior can be cached (``cache=``) and checked for
       calibration (:meth:`~ampere.inference.SBIEngine.calibrate`); and a
       seeded problem's run repeats bitwise, network included. See
       :doc:`sbi` for the worked tutorial.

The same fit, side by side
---------------------------

Both snippets below fit the shape of
``examples/minimal_working_example.py``'s toy model — a straight line in
wavelength, fit to a synthetic spectrum with emcee. They are abridged (no
synthetic-data generation, no plotting, priors elided) to show how the
composition differs, not to run standalone.

Legacy:

.. code-block:: python

    from ampere.data import Spectrum
    from ampere.infer.emceesearch import EmceeSearch
    from ampere.models import Model

    class ASimpleModel(Model):
        def __call__(self, slope, intercept, **kwargs):
            self.modelFlux = slope * self.wavelength + intercept
            return {"spectrum": {"wavelength": self.wavelength, "flux": self.modelFlux}}
        # lnprior / prior_transform also required — see the full example

    model = ASimpleModel(wavelengths)
    spectrum = Spectrum(wave, flux, uncertainty, "um", "Jy")

    optimizer = EmceeSearch(model=model, data=[spectrum], nwalkers=100)
    optimizer.optimise(nsamples=150, burnin=100, guess=guess)
    optimizer.postProcess()

v2, on the reference backend:

.. code-block:: python

    import scipy.stats as st

    from ampere.core import (
        Dataset, FittingProblem, GaussianFamily, IndependentNoise, Likelihood,
        Model, Parameter, Spectrum,
    )
    from ampere.inference import EmceeEngine

    class ASimpleModel(Model):
        def __init__(self, wavelength):
            self.register_buffer("wavelength", wavelength)
            self.register_parameter(Parameter("slope", st.uniform(-10, 20)))
            self.register_parameter(Parameter("intercept", st.uniform(-10, 20)))

        def evaluate(self, **values):
            ctx = self.context(values)
            flux = ctx["slope"] * ctx["wavelength"] + ctx["intercept"]
            return Spectrum(ctx["wavelength"], flux)

    model = ASimpleModel(wavelength)
    observed = Spectrum(wave, flux, uncertainty=uncertainty)
    dataset = Dataset(observed, likelihood=Likelihood(GaussianFamily(), IndependentNoise()))
    problem = FittingProblem(model, [dataset], seed=20260908)

    run = EmceeEngine(problem, walkers=32).run(steps=2000, burn_in=1000)

The shapes that carry over: a model class producing a prediction, a data
container, an emcee-driven sampler you configure and run. What moves: priors
attach to named :class:`~ampere.core.Parameter` declarations instead of a
hand-written ``lnprior``/``prior_transform`` pair; the noise model is a
named, composable choice (:class:`~ampere.core.IndependentNoise` here —
:doc:`overview`'s own worked example swaps in
:class:`~ampere.core.GaussianProcessNoise` for the flexible likelihood)
rather than parameters built into the data class; and the run is driven by
binding a :class:`~ampere.inference.EmceeEngine` to a
:class:`~ampere.core.FittingProblem` rather than mixing the sampler and the
model together in one ``EmceeSearch`` object.

.. _migrating-every-name:

Every legacy name, and where it goes
-------------------------------------

One table per legacy package, one row per public name, enumerated from the
code rather than from memory. The **v2 route** is a real v2 name, a worked
twin (the v2 versions of six legacy examples, described in
:ref:`migrating-side-by-side`), or **no equivalent; carried**. "Carried"
means v2 does not do the thing and has no plan to; the code keeps working,
indefinitely, on ``ampere.legacy`` under the policy in :doc:`legacy`
(D1 (b)) — it is not deprecated and nothing is scheduled for removal. The
legacy names are written as literals because the legacy reference is
rendered without cross-reference targets; every v2 name is a link.

Several legacy names are not what their package's ``__init__`` exports. The
``data`` package exports only ``Photometry`` and ``Spectrum``; the ``models``
package star-imports ``models``, ``blackbodies`` and ``powerlaws``; the
``infer`` and ``utils`` packages export nothing (``__all__`` is empty), so
most names below are reached by their module path, as the minimal working
examples do (``from ampere.infer.emceesearch import EmceeSearch``).

Data (``ampere.legacy.data``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 24 40 36

   * - Legacy name
     - v2 route
     - Note
   * - ``Photometry``
     - :class:`~ampere.core.Dataset` binding a
       :class:`~ampere.core.PhotometricPoints` observation to an
       :class:`~ampere.core.Instrument` of
       :meth:`SyntheticPhotometry.from_library <ampere.backends.reference.SyntheticPhotometry.from_library>`
     - The filter names move from the data class to the instrument step;
       ``reloadFilters`` has no counterpart because the model's grid is
       negotiated, not set by hand. ``Photometry.fromFile`` has no shipped
       reader: the twins read their VOTables with astropy in each package's
       ``generators.py`` (see :ref:`migrating-twin-modified-blackbody`).
   * - ``Spectrum``
     - :class:`~ampere.core.Spectrum` observed in a
       :class:`~ampere.core.Dataset`, with an
       :class:`~ampere.core.Instrument` chain of
       :class:`~ampere.backends.reference.Resample`,
       :class:`~ampere.backends.reference.LSFConvolution` and
       :class:`~ampere.backends.reference.CalibrationScale`
     - ``calUnc`` becomes the calibration scale's prior; ``scaleLengthPrior``
       and ``covWeightPrior`` become the kernel of a
       :class:`~ampere.core.GaussianProcessNoise` (see
       :ref:`migrating-twin-linear-sed`). ``Spectrum.fromFile`` has no
       shipped reader either.
   * - ``LineStrengths``
     - no equivalent; carried
     - A subclass of ``Photometry`` whose own docstring calls it "a
       placeholder": it adds nothing to ``Photometry``. v2 has no
       line-strength container.
   * - ``Data``
     - no equivalent; carried
     - The base class, not exported by the package. v2's counterpart of the
       idea is the container family (:class:`~ampere.core.Spectrum`,
       :class:`~ampere.core.PhotometricPoints`, :class:`~ampere.core.Image`
       and the rest) plus a :class:`~ampere.core.Dataset`.
   * - ``Image``
     - :class:`~ampere.core.Image`, fitted as in :doc:`image`
     - The legacy class is a stub: its methods raise ``NotImplementedError``.
       The v2 container is real, with a PSF-convolution step and a
       Hilbert-space GP over the grid.
   * - ``Interferometry``
     - :class:`~ampere.core.VisibilitySet` and
       :class:`~ampere.core.ClosurePhases`, fitted as in :doc:`interferometry`
     - Also a stub on the legacy side. v2 fits visibilities and closure
       phases end to end.
   * - ``Cube``
     - :class:`~ampere.core.Cube`
     - Also a stub on the legacy side. The v2 container exists; there is
       no tutorial for a cube fit.

Models (``ampere.legacy.models``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 24 40 36

   * - Legacy name
     - v2 route
     - Note
   * - ``Model``
     - :class:`~ampere.core.Model`
     - Priors attach to declared :class:`~ampere.core.Parameter` objects
       instead of a hand-written ``lnprior``/``prior_transform`` pair.
   * - ``CompositeModel``
     - one :class:`~ampere.core.Model` whose ``evaluate`` sums the
       components
     - Legacy combined two models with an operator object. In v2 you write
       the sum: ``StarDisc.evaluate`` in :ref:`migrating-twin-star-disc` is
       a photosphere plus a dust term.
   * - ``AnalyticalModel``, ``RTModel``
     - no equivalent; carried
     - Thin base classes. v2 draws no analytic/radiative-transfer
       distinction; a black-box simulator is an ordinary
       :class:`~ampere.core.Model` (``CarbonStarShell`` in
       :ref:`migrating-twin-cstar`).
   * - ``ModelResults``
     - :class:`~ampere.core.ModelResult`
     - A dictionary of named channels, each a typed container, in place of
       a ``{"spectrum": {"wavelength": ..., "flux": ...}}`` dictionary.
   * - ``SingleModifiedBlackBody``
     - :class:`~ampere.backends.reference.ModifiedBlackBody`; the paper's
       four-parameter form is ``examples/modified_blackbody``
     - Defined twice, in ``blackbodies.py`` and ``PowerLawAGN.py``; the
       ``ampere.legacy.models`` namespace gets the first.
   * - ``DualBlackBodyDust``
     - no equivalent as a class; carried. The pieces are
       :class:`~ampere.backends.reference.BlackBody` and
       :class:`~ampere.backends.reference.ModifiedBlackBody`
     - Three different classes share this name, in ``blackbodies.py``,
       ``starScreen.py`` and ``PowerLawAGN.py``. A model summing two
       blackbodies is the ``StarDisc`` pattern.
   * - ``PowerLawContinuumAbsoluteAbundances``,
       ``PowerLawContinuumRelativeAbundances``
     - no equivalent; carried. The nearest v2 pieces are
       :class:`~ampere.backends.reference.PowerLaw` and the opacity-table
       model of ``examples/ngc6302``
     - Dust-composition models over tabulated opacities; neither is
       shipped as a v2 class.
   * - ``PowerLawAGN``, ``PowerLawAGNRelativeAbundances``,
       ``OpacitySpectrum``
     - no equivalent; carried
     - The same family, in ``PowerLawAGN.py``.
   * - ``PolynomialSource``, ``PolynomialSource2``
     - no equivalent; carried
     - Polynomial stellar-screen models in ``starScreen.py``.
   * - ``QuickSEDModel``
     - ``StarDisc`` in ``examples/star_disc``
     - See :ref:`migrating-twin-star-disc`. The photosphere comes from
       ampere's own PHOENIX emulator rather than Starfish.
   * - ``HyperionCStarRTModel``
     - ``CarbonStarShell`` in ``examples/cstar``
     - See :ref:`migrating-twin-cstar`. Needs the pixi ``hyperion``
       environment.
   * - ``DustySpectrum``
     - no equivalent; carried
     - Wraps the external DUSTY code through a template file.
   * - ``BaseExtinctionLaws``, ``F99Extinction``, ``CCMExtinctionLaw``
     - no shipped equivalent; carried. v2 models apply extinction
       themselves, as ``examples/phoenix_star`` does with its own CCM89
     - ``F99Extinction`` is what the ``extinction`` extra
       (:doc:`install`) is for. v2 ships no extinction-law library.
       ``extinction_models.py`` is an empty file beside the live
       ``extinctionModels.py``, which defines all three classes.

Inference (``ampere.legacy.infer``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 24 40 36

   * - Legacy name
     - v2 route
     - Note
   * - ``BaseSearch``
     - :class:`~ampere.inference.Engine`
     - A search owned its model and data; an engine is bound to a
       :class:`~ampere.core.FittingProblem` and returns its run.
   * - ``MCMCSampler``, ``EnsembleSampler``
     - :class:`~ampere.inference.Engine`
     - Intermediate bases with no v2 counterpart: each engine is concrete.
   * - ``EmceeSearch``
     - :class:`~ampere.inference.EmceeEngine`
     - ``optimise(nsamples=, burnin=, guess=)`` becomes ``run(steps,
       burn_in=, initial=)``: ``initial`` is a ``(walkers, n_dim)`` array of
       start positions or an :class:`~ampere.results.Optimum` from
       :func:`ampere.inference.optimise` (the walkers then start in a ball
       at its mode); left out, the walkers start from the prior.
   * - ``ZeusSearch``
     - :class:`~ampere.inference.ZeusEngine`
     -
   * - ``BaseNestedSampler``
     - :class:`~ampere.inference.Engine`
     - The v2 nested samplers are :class:`~ampere.inference.DynestyEngine`,
       :class:`~ampere.inference.NautilusEngine` and
       :class:`~ampere.inference.UltranestEngine`.
   * - ``DynestyNestedSampler``
     - :class:`~ampere.inference.DynestyEngine`
     -
   * - ``DynestyDynamicNestedSampler``
     - :class:`~ampere.inference.DynestyEngine` with ``dynamic=True``
     - One class, a flag.
   * - ``LFIBase``
     - :class:`~ampere.inference.SBIEngine`
     - The base class of the legacy SBI search.
   * - ``SBI_SNPE``
     - :class:`~ampere.inference.SBIEngine`
     - See the concept map above and :doc:`sbi`. The embedding-network
       dictionary documented in the legacy ``Embedding_nets`` notebook is
       read verbatim by ``embedding=``.
   * - ``CustomPriorWrapper``
     - no equivalent; carried
     - Adapted a prior for the ``sbi`` package. v2 priors are declared on
       the parameters and :class:`~ampere.inference.SBIEngine` builds what
       ``sbi`` needs.
   * - ``ScipyMinMixin``
     - :func:`ampere.inference.optimise`
     - Point estimates and starting values; see :doc:`optimisers`.
   * - ``SimpleMCMCPostProcessor``, ``SBIPostProcessor``,
       ``ArvizPostProcessor``
     - :mod:`ampere.results`
     - See :ref:`migrating-post-processing`. ``ArvizPostProcessor`` is an
       empty subclass of the shared base: ArviZ is the one format of v2.
   * - ``SimulatorMixin``
     - :attr:`FittingProblem.simulate_many <ampere.core.FittingProblem.simulate_many>`
     - The legacy mixin's ``simulate`` is a stub; the v2 method is the
       real thing, and the training-set writer builds on it.
   * - ``AnotherMixin``
     - no equivalent; carried
     - An empty class (``pass``) in ``mixins.py``.
   * - ``create_worker_init``, ``lnlike_pool``, ``lnprob_pool``
     - :class:`~ampere.core.ProcessExecutor`, for pooled simulation
     - Module-level helpers in ``basesearch.py`` that give a
       multiprocessing pool its worker state. v2 pools simulations through
       an executor rather than through the search.

Utilities (``ampere.legacy.utils``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 24 40 36

   * - Legacy name
     - v2 route
     - Note
   * - ``pyphot_compat`` (``get_unit``)
     - no equivalent; carried
     - Not only legacy: v2's
       :meth:`SyntheticPhotometry.from_library <ampere.backends.reference.SyntheticPhotometry.from_library>`
       names this module's ``get_unit`` as the way to attach units for
       pyphot 2, so it stays because the v2 library depends on it.
   * - ``makeFilterSet`` (two copies: ``makeFilterSet.py`` and
       ``makeFilterSet_mod.py``), ``getFilterList``, ``get_filter_svo``,
       ``get_filter_file``
     - no equivalent; carried
     - Build a pyphot filter library from a CSV or from the SVO service.
       v2 uses the bundled library
       (:func:`~ampere.backends.reference.bundled_filter_library`) and
       accepts a ``library=`` of your own in ``from_library``.
   * - ``make_alma_filter``
     - no equivalent; carried
     - A filter built from frequency windows. For a top-hat in v2, the
       ``examples/star_disc`` twin builds ALMA and ATCA bands with
       ``detector="energy"`` (:ref:`migrating-twin-star-disc`).
   * - ``read_eso_filt``, ``eso_filt_to_csv``
     - no equivalent; carried
     - ESO filter-file converters, for building a library by hand.

Logger (``ampere.legacy.logger``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 24 40 36

   * - Legacy name
     - v2 route
     - Note
   * - ``Logger``, ``handle_exception``
     - no equivalent; carried
     - A mixin that sets up the ``"ampere logger"`` for a legacy search
       (file and terminal handlers), and an uncaught-exception hook for it.
       v2 engines do not write log files.

.. _migrating-side-by-side:

Side by side: the six examples
-------------------------------

Under D11 the legacy example scripts stay as they are and keep working on
``ampere.legacy``; each of six has a v2 **twin** beside it, the same model,
the same data and the same question written once against
:mod:`ampere.core` and the reference backend (``examples/README.md`` has the
pairs table). Each subsection below puts the legacy lines next to the twin's
own code, abridged, and says what changed in the translation. The full
version of every snippet is in the named file and function; the twins are
covered by ``tests/examples/`` and are not modified by this page. Run them as
``python -m examples.<twin>`` from the repository root; :doc:`tutorials`
describes them alongside the other worked examples.

.. _migrating-twin-linear-sed:

Linear SED: ``minimal_working_example.py`` and ``examples/linear_sed``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Legacy ``minimal_working_example.py`` and its ``_dynesty``, ``_zeus``,
``_sbi`` and ``_sbi_embedding`` variants; twin ``examples/linear_sed/``, whose
``--engine`` flag chooses between them.

.. code-block:: python

    # legacy
    spec0 = Spectrum(irs.wavelength, flux0, unc0, "um", "Jy",
                     calUnc=0.0025, scaleLengthPrior=0.01)
    photometry = Photometry(filterName=names, value=values,
                            uncertainty=unc, photunits="Jy", libName=libname)
    photometry.reloadFilters(wavelengths)
    optimizer = EmceeSearch(model=model, data=[spec0, spec1, photometry],
                            nwalkers=100, moves=m, vectorize=False)
    optimizer.optimise(nsamples=150, burnin=100, guess=guess)
    optimizer.postProcess()

.. code-block:: python

    # v2: linear_sed.build_instruments / _likelihood / build_problem / fit
    sl = Instrument([Resample(sl_grid), CalibrationScale(st.lognorm(0.0025, scale=1.0))],
                    channel="sed", label="sl")
    catalogue = Instrument([SyntheticPhotometry.from_library(FILTERS, GRID)],
                           channel="sed", label="catalogue")
    kernel = Matern32(GP_AMPLITUDE_PRIOR, st.halfnorm(scale=0.01),
                      amplitude_unit=u.Jy, length_scale_unit=u.um,
                      axes=("spectral_axis",))
    gp = Likelihood(GaussianFamily(), GaussianProcessNoise(kernel, DenseGP()))
    datasets = DatasetCollection({
        "catalogue": Dataset(obs_phot, catalogue,
                             likelihood=Likelihood(GaussianFamily(), IndependentNoise())),
        "sl": Dataset(obs_sl, sl, likelihood=gp),
    })
    problem = FittingProblem(model, datasets, seed=20260928)
    run = EmceeEngine(problem, walkers=100).run(150, burn_in=100)

The legacy noise pair is the thing to translate: ``calUnc`` ("the sigma for a
LogNormal distribution" on the calibration scale) becomes each spectrum's
:class:`~ampere.backends.reference.CalibrationScale` prior, and
``scaleLengthPrior`` (a half-normal on the correlated noise's length scale)
becomes the length-scale prior of the kernel of a
:class:`~ampere.core.GaussianProcessNoise`. ``covWeightPrior`` has no number
to carry over (its units differ from a Matern amplitude's), so the kernel
amplitude gets a weakly informative prior of its own. The instruments carry
distinct labels because several of them share one channel.
:doc:`photometry_spectra` is the tutorial for this composition.

.. _migrating-twin-ngc6302:

NGC 6302: ``NGC6302.py`` and ``examples/ngc6302``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Legacy ``NGC6302.py``, ``NGC6302_zeus.py`` and
``NGC6302-calculate-dust-mass.py``; twin ``examples/ngc6302/``, with the
dust-mass calculation as ``dust_mass.py``, a function over a fit's posterior
draws.

.. code-block:: python

    # legacy
    class SpectrumNGC6302(Model):
        def lnprior(self, theta, **kwargs):
            # flat box on every parameter, plus the ordering
            # theta[12] > theta[11] and theta[14] > theta[13]
            ...
    spec = Spectrum(wave, flux, unc, "um", "Jy",
                    calUnc=1e-10, scalelengthPrior=0.1)
    optimizer = EmceeSearch(model=model, data=dataset, nwalkers=50, moves=m)
    optimizer.optimise(nsamples=50000, burnin=40000, guess=guess)
    optimizer.postProcess()

.. code-block:: python

    # v2: ngc6302.build_model / build_problem / fit
    model = KemperTwoShell(
        DEFAULT_GRID,
        logacold0=st.uniform(-6.0, 6.0),             # eleven of these
        Tcold0=st.triang(c=0, loc=10.0, scale=70.0),
        Tcold_fraction=st.uniform(0.0, 1.0),
        Twarm0=st.triang(c=0, loc=80.0, scale=100.0),
        Twarm_fraction=st.uniform(0.0, 1.0),
    )
    instrument = Instrument([Resample(observed_wavelength),
                             CalibrationScale(CALIBRATION_PRIOR)],
                            channel="sed", label="iso")
    likelihood = Likelihood(GaussianFamily(),
                            GaussianProcessNoise(kernel, QuasisepGP()))
    problem = FittingProblem(model, DatasetCollection(
        {"iso": Dataset(observed, instrument, likelihood=likelihood)}), seed=seed)
    run = EmceeEngine(problem, walkers=50).run(steps, burn_in=burn_in)

Two things changed. The legacy hard rejection that keeps the inner
temperature above the outer is reproduced **exactly** and without a
rejection step: ``T0`` gets a triangular prior, a fraction gets a uniform
one, and ``T1 = T0 + fraction * (T_max - T0)`` is computed in ``evaluate``;
the Jacobian cancels the triangular density, so the joint prior is flat on
the ordered triangle, as in legacy. And the Matern-3/2 GP runs on the
O(N) :class:`~ampere.core.QuasisepGP` solver rather than a dense one, which
is exact for that kernel on sorted one-dimensional coordinates. The module
docstring of ``examples/ngc6302/ngc6302.py`` has the derivation.

.. _migrating-twin-modified-blackbody:

Modified blackbody: ``modifiedblackbody.py`` and ``examples/modified_blackbody``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Legacy ``examples/examples_paper/modifiedblackbody.py``; twin
``examples/modified_blackbody/``, where one ``--engine`` flag (or ``--all``)
replaces the four blocks the script ran by hand. The ``Ampere_MBB_Example``
notebook is this example in notebook form (see :doc:`tutorials`).

.. code-block:: python

    # legacy
    class ModifiedBlackBody(Model):
        def __init__(self, wavelength, kappa=10, kappawave=250, lims=None, **kw):
            self.priors = [uniform(lim[0], lim[1] - lim[0]) for lim in lims]
        def __call__(self, t, logm, beta, d, **kwargs):
            ...
            return {"spectrum": {"wavelength": self.wavelength, "flux": flux}}
        # lnprior and prior_transform loop over self.priors
    optimizer_e = EmceeSearch(model=model, data=[photometry], nwalkers=100, moves=m)
    optimizer_e.optimise(nsamples=1000, burnin=900, guess="None")
    optimizer_e.postProcess()

.. code-block:: python

    # v2: modified_blackbody.ModifiedBlackBody / build_model / fit
    class ModifiedBlackBody(Model):
        def __init__(self, wavelength, *, temperature, logmass, beta, distance, ...):
            self.register_buffer("wavelength", grid, unit=u.um)
            self.register_parameter(_as_parameter("temperature", temperature, unit=u.K))
            ...
        def evaluate(self, **values):
            ctx = self.context(values)
            ...
            return ModelResult({self.channel: spectrum})

    instrument = Instrument([SyntheticPhotometry.from_library(FILTERS, GRID)],
                            channel="sed", label="catalogue")
    run = EmceeEngine(problem, walkers=100).run(1000, burn_in=900)

The formula is transcribed exactly, including two quirks that the twin keeps
rather than corrects (the script's "kpc" distance is converted as parsecs,
and astropy's ``BlackBody`` is called without units, so the flux is not
physically in Jy; the recovery under noise is what the example shows). The
priors and the ten-band catalogue carry over unchanged. The shipped
:class:`~ampere.backends.reference.ModifiedBlackBody` is the three-parameter
form that folds distance and mass into one scale; the twin keeps the paper's
four on purpose.

.. _migrating-twin-phoenix-star:

PHOENIX star: ``phoenixstar.py`` and ``examples/phoenix_star``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Legacy ``examples/examples_paper/phoenixstar.py``; twin
``examples/phoenix_star/``, which also runs NUTS on torch or jax through
native twins of the emulator.

.. code-block:: python

    # legacy
    class StarfishStellarModel(Model):
        def __call__(self, logl, teff, logg, av, rv):
            fl = self.emulator.load_flux([teff, logg, 0.0])
            ...
            fl = fl / (4 * np.pi * self.distance**2)
            ext = ccm89(self.wavelength * 1e4, av, rv)
            return {"spectrum": {"wavelength": self.wavelength,
                                 "flux": apply(ext, fint)}}
    spec0 = Spectrum(specwaves, spec0, unc0, "um", "Jy",
                     calUnc=0.0025, scaleLengthPrior=0.01)
    optimizer = SBI_SNPE(model=model, data=dataset, name="star_test_sbi_spec")
    optimizer.optimise(nsamples=10000, nsamples_post=10000)

.. code-block:: python

    # v2: phoenix_star.build_instruments / fit
    catalogue = Instrument([SyntheticPhotometry.from_library(FILTERS, GRID)],
                           channel="sed", label="catalogue")
    rvs = Instrument([Resample(data.segment_b),
                      CalibrationScale(CALIBRATION_PRIOR)],
                     channel="rvs", label="rvs")
    ...
    store = ArtefactStore(CACHE_DIR)
    run = SBIEngine(problem, budget=10_000, rounds=1, cache=store).run(10_000)

Three things changed, ruled and recorded in the twin's docstring. There is
**no Starfish**: the emulator is ampere's own, a PCA plus one Gaussian
process per weight, committed as ``phoenix_emulator.npz``, because Starfish
on PyPI still needs Python below 3.10. The legacy divides by :math:`4\pi d^2`
a second time after a flux that already carries the :math:`4\pi`, so its
fluxes are :math:`4\pi` too faint; the twin does not reproduce that, and a
luminosity fitted with the legacy script is :math:`4\pi` times the twin's for
the same photometry (the note under :ref:`migrating-twin-star-disc` has the
detail). And CCM89 is implemented inside the model, so it is
differentiable on every backend, instead of calling the ``extinction``
package. The trained posterior is cached by an
:class:`~ampere.results.ArtefactStore`, keyed on the problem's hashes.

.. _migrating-twin-cstar:

Carbon star: ``cstar_model_test_sbi_v2.py`` and ``examples/cstar``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Legacy ``examples/cstar_model_test_sbi_v2.py`` and its ``_embedding``
variant; twin ``examples/cstar/``. It needs Hyperion, which is an
example-only requirement in the pixi ``hyperion`` environment
(``pixi install -e hyperion``), not an ampere extra.

.. code-block:: python

    # legacy
    p = data.Photometry.fromFile("Observed_SED.vot", format="votable", libName=lib)
    s = data.Spectrum.fromFile("SPEC_OGLE_CAGB_IRS.csv", "User-Defined", ...)
    optimizer_sbi = SBI_SNPE(model=model, data=[p, *s], n_rounds=2,
                             name="Cstar_test_sbi_v2", nproc=70)
    optimizer_sbi.optimise(nsamples=10000, nsamples_post=10000)

.. code-block:: python

    # v2: cstar.build_problem / fit
    problem = FittingProblem(
        build_model(photons=photons, carbon=carbon),
        DatasetCollection(datasets),
        seed=seed,
        validate=False,
        simulator_failures=(SimulatorFailed,),
    )
    run = SBIEngine(problem, budget=10_000, rounds=2,
                    executor=ProcessExecutor(4), cache=store).run(10_000)

The legacy model exposed two abundances with a Dirichlet(1, 1) prior and a
prior transform that normalised two gamma quantiles. A Dirichlet(1, 1) over
two fractions that sum to one *is* a uniform on either, so the twin has one
parameter, ``sic_fraction``, uniform on [0, 1]: the legacy's eight
parameters were seven degrees of freedom. The Hyperion run is a black-box
:class:`~ampere.core.Model` that declares its failures through
``simulator_failures=``, and the pooled simulations go through a
:class:`~ampere.core.ProcessExecutor`. The twin uses ``miepython`` where the
legacy used an unpackaged ``bhmie`` step, and does not reproduce the
embedding variant's emcee block: a likelihood sampler needs far more Hyperion
runs than anyone can wait for.

.. _migrating-twin-star-disc:

Star and debris disc: ``star_disc.py`` and ``examples/star_disc``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Legacy ``examples/star_disc.py`` (the HD 105 SED); twin
``examples/star_disc/``. ``QuickSEDModel`` was the one legacy model class
with no v2 counterpart; ``StarDisc`` is that counterpart.

.. code-block:: python

    # legacy
    model = QuickSED.QuickSEDModel(wavelengths, lims=np.array([...]),
                                   dstar=1000.0 / parallax,
                                   starfish_dir=path + "../",
                                   fbol_1l1p=fbol_1l1p)
    optimizer = EmceeSearch(model=model, data=dataSet, nwalkers=40, vectorize=False)
    optimizer.optimise(nsamples=4000, nburnin=1000, guess=guess)
    optimizer.postProcess(logx=True, logy=True)

.. code-block:: python

    # v2: star_disc.StarDisc.evaluate / fit
    def evaluate(self, **values):
        ctx = self.context(values)
        emitted = {}
        for channel in self.CHANNELS:           # "sed" and "rvs"
            total = self.star(channel, ctx) + self.dust(channel, ctx)
            emitted[channel] = self.templates[channel].with_values(total)
        return ModelResult(emitted)

    run = EmceeEngine(problem, walkers=40).run(4000, burn_in=1000)

The star is the PHOENIX emulator of :ref:`migrating-twin-phoenix-star`, with
the corrected inverse-square law and no extinction (``QuickSED`` has none);
the dust is transcribed exactly from the legacy ``__call__``. The model
emits two channels, ``"sed"`` and ``"rvs"``, so the Gaia-RVS spectrum is
modelled on its own grid. Two photometric bands that the bundled filter
library lacks (ALMA band 6 and ATCA 9 mm) are top-hats with
``detector="energy"``.

.. note::

   **The legacy scripts divide by** :math:`4\pi d^2` **twice, so their fluxes
   are** :math:`4\pi` **too faint.** ``examples/examples_paper/phoenixstar.py``
   divides the emulator's flux by :math:`4\pi d^2` in the model (line 53) after
   scaling the bolometric flux by :math:`1/4\pi (1\,\mathrm{pc})^2` (line 25),
   and again in its post-processing (line 130). ``examples/star_disc.py``
   forms ``fbol_1l1p / (4π (1 pc)²)`` (line 43) and passes it to
   ``ampere/legacy/models/QuickSED.py``, which divides by :math:`4\pi d_\star^2`
   a second time (line 148), so for ``QuickSED`` the second division is the
   example feeding the model as much as the model itself. A luminosity fitted
   with either legacy script is :math:`4\pi` times the twin's for the same
   photometry. The twins use the correct law. The legacy code is frozen (D11),
   so this note is the only fix.

.. _migrating-post-processing:

The post-processing
--------------------

Every legacy search ended in ``postProcess()``, which called the plotting
and summary methods its post-processor mixin supplied
(``print_summary``, ``plot_corner``, ``plot_trace`` and
``plot_posteriorpredictive``; the nested sampler's version also drew a
summary plot). Those mixins — ``SimpleMCMCPostProcessor`` for the ensemble
samplers, ``SBIPostProcessor`` for the SBI search, ``ArvizPostProcessor`` as
an empty subclass — hung figures and summaries on a sampler class, and
``ModelResults`` carried a model's output. In v2 the run **is** the result:
an engine returns an ArviZ ``DataTree`` with provenance attached, and
everything that reads it is a function in :mod:`ampere.results`, written once
against that format and independent of the sampler.

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Legacy
     - v2
   * - ``print_summary``, ``summary``, ``get_credible_interval``
     - :func:`ampere.results.summary` over the run
   * - ``get_ess``, ``get_rhat``
     - :func:`ampere.results.effective_sample_sizes` and the convergence
       columns of :func:`~ampere.results.summary`
   * - ``plot_corner``
     - :func:`ampere.results.plot_corner`
   * - ``plot_trace``
     - :func:`ampere.results.plot_trace`
   * - ``plot_posteriorpredictive``
     - :func:`ampere.results.add_posterior_predictive` then
       :func:`ampere.results.plot_posterior_predictive`
   * - ``plot_covmats``
     - no equivalent plot; the GP is a property of the likelihood, and
       :func:`ampere.results.plot_residuals` and
       :func:`ampere.results.residual_whiteness` test what it left
   * - ``get_map``
     - :func:`ampere.inference.optimise`, which returns an
       :class:`~ampere.results.Optimum`
   * - pickling the whole search object
     - :func:`ampere.results.to_netcdf` and
       :func:`ampere.results.from_netcdf`
   * - the trained network of ``SBI_SNPE``
     - cached by an :class:`~ampere.results.ArtefactStore`

.. code-block:: python

    # legacy
    optimizer.optimise(nsamples=150, burnin=100, guess=guess)
    optimizer.postProcess()

    # v2
    from ampere.results import (
        add_posterior_predictive, plot_corner, plot_posterior_predictive,
        plot_trace, summary,
    )

    run = EmceeEngine(problem, walkers=100).run(150, burn_in=100)
    print(summary(run))
    plot_corner(run)
    plot_trace(run)
    plot_posterior_predictive(add_posterior_predictive(run, problem))

The plotting functions each return a matplotlib figure and write nothing to
disk. The ``quickstart`` notebook (:doc:`tutorials`) runs this sequence on a
real fit.

What has no legacy counterpart
-------------------------------

Nothing in the legacy API corresponds to these; each has its own page.

* **Gradient samplers.** :class:`~ampere.inference.NUTSEngine` and
  :class:`~ampere.inference.VIEngine` on the torch and jax backends: see the
  capability ladder in :doc:`overview`.
* **The optimisers.** :func:`ampere.inference.optimise` with its multi-start
  and warm-start options: :doc:`optimisers`.
* **Approximate GP solvers.** The sparse and structured solvers behind the
  flexible likelihood: :doc:`solvers`, and the kernels in :doc:`kernels`.
* **Populations.** :class:`~ampere.core.Population`, fitted jointly or by
  reweighting archived fits: :doc:`population`.
* **Other data.** Images, visibilities, closure phases, astrometry and
  time series as fitted modalities: :doc:`image`, :doc:`interferometry`,
  :doc:`astrometry`.
* **Non-Gaussian likelihoods.** A registered likelihood family such as a
  profiled Cash statistic: :doc:`wstat_comparison`.
* **Training sets.** :func:`ampere.results.write_training_set` and
  :func:`ampere.results.read_training_set`, for amortised inference: see
  :doc:`sbi`.
* **The misspecification study.** The flagship validation of the flexible
  likelihood: :doc:`m2_misspecification`.

The policy for the legacy code, and where it lives, is :doc:`legacy`; the
legacy reference is :doc:`ampere.legacy`.
