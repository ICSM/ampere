Migrating to the new API
=========================

Ampere has two APIs, and this page is about the gap between them: what maps
to what, what "frozen" means for the old one, and the same small fit shown
both ways. It is a seed, not the guide: the legacy API stays **linked**, not
removed, and a full migration guide and a deprecation timetable are planned
for a later phase. See
:doc:`overview` for how the v2 pieces fit together and :doc:`ampere` for the
legacy reference.

What "frozen" means
--------------------

**Legacy ampere** — ``ampere.data``, ``ampere.models``, ``ampere.infer`` — is
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

A deprecation policy — if and when legacy pieces are retired, and on what
notice — is Phase 6 work and is not written yet.

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
     - One engine covers both the static and dynamic nested-sampling modes.
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
       :meth:`~ampere.core.dataset.FittingProblem.simulate_many` rather than
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

What comes next
----------------

This page grows into the full migration guide, and gains the Phase 6
deprecation policy for the legacy API, once that phase is scoped. Until
then, :doc:`overview` is the complete reference for composing a v2 problem,
and :doc:`ampere` is where the legacy API stays documented and supported as
frozen code.
