Advanced usage
==============

The previous guides described how to do fairly simple things with AMPERE like defining a model. However, you will often want to do something more complex than that. These tutorials are intended to help you achieve that.

Model comparison
----------------

Often our objective with inference is not only to determine which distribution of parameter values is supported by our observations, but instead to determine which model (out of some set of 2 or more models) is most probable, given our observations.
This is the *model selection* problem, rather than parameter estimation.

AMPERE natively supports model selection with Nested Sampling - use :class:`~ampere.inference.DynestyEngine` to do your inference instead of any of the other engines, and you will get an estimate of the *model evidence* at the end of the run.
Do this for multiple different models, and compare the logarithm of the evidence at the end - whichever model has the highest evidence is the preferred model (assuming all models are equally probable *a priori*).
This is very convenient, since it automatically penalises models with different numbers of parameters, and doesn't require that models be nested.
However, this is only well-justified if the model parameters have physical meaning.

At present, if you want to do model comparison for other approaches than Nested Sampling, you will have to roll your own.
For computing the evidence from *emcee* or *zeus* results, the `harmonic <https://astro-informatics.github.io/harmonic/index.html>`_ package may be effective.
Alternatively you can do model comparison with `arviz <https://www.arviz.org/en/latest/>`_'s WAIC and/or LOO-PIT methods.
These methods have the advantage of giving meaningful results even if the parameters are ad hoc without real physical meaning, and nothing needs exporting: every v2 run **is** an ArviZ ``DataTree`` already, carrying the per-draw and per-dataset log-likelihood decomposition those methods need (see :doc:`overview`).

More detailed tutorials will be available soon!


Noise models: beyond the default kernel
------------------------------------------

:doc:`concept` introduces the flexible likelihood with its default
Matérn-3/2 kernel; :doc:`kernels` is the full reference for
:class:`~ampere.core.GaussianProcessNoise`'s algebra, for the cases the
default does not cover.

**A residual with a period** — interference fringing, an instrumental
ripple — wants :class:`~ampere.core.SHO`, celerite's damped
simple-harmonic-oscillator term, rather than a stationary Matérn kernel that
has no notion of periodicity at all. Composed with the default through
:class:`~ampere.core.Sum`, a broad Matérn plus a narrow ``SHO`` still costs
O(N) on :class:`~ampere.core.QuasisepGP` — a sum of quasiseparable terms is
quasiseparable, because the semiseparable generators concatenate and the
ranks simply add. ``examples/m2_misspecification/fringing.py`` demonstrates
exactly this: refitting M2's ``fringing`` scenario with ``Matern32 + SHO``
against the stationary default, with the bias and calibration improvement
reported.

**A residual correlated in two different ways at once** — a missing patch
of sky with a spectral profile, smooth across spatial frequency and sharp
across wavelength — wants a :class:`~ampere.core.Product` of two kernels,
each acting on a named subset of the container's axes through the
``axes=`` selector, rather than one isotropic kernel that cannot express
two structures on two coordinates simultaneously. :doc:`interferometry` §6
measures this case end to end.

**A kernel ampere does not ship** works on the dense solver the moment you
write its covariance function — :class:`~ampere.core.Kernel` is a public
ABC — and reaches the O(N) path once you register its celerite
representation with :func:`~ampere.core.register_quasiseparable_term`, out
of tree, with no change to ``ampere.core``. :doc:`kernels` §5 walks through
the registration with the same trivial worked example
``tests/core/test_kernels.py`` uses to prove the route.

Very slow models
----------------

AMPERE allows you to use a wide variety of models to interpret your data, which may include models which take a very long time to compute, or whose likelihood cannot be written down at all — a compiled radiative-transfer code behind a Python call, say.
In such cases, simulation-based inference is the right tool: :class:`~ampere.inference.SBIEngine` trains a neural network on simulated ``(theta, x)`` pairs, drawn through :meth:`~ampere.core.dataset.FittingProblem.simulate_many`, instead of consuming ``log_prob``.
It runs under a process pool with per-simulation timeouts and crash capture, so a slow or occasionally-crashing external simulator is the case it is built for rather than an edge case it tolerates.
:doc:`sbi` is the full tutorial — a black-box simulator fitted end to end, caching the trained posterior, truncated marginal ratio estimation for a tighter fit, and checking the result is calibrated; the runnable scripts are ``examples/sbi/``.

.. note::

   This is the **current** route. ``ampere.infer.sbi`` is legacy and frozen
   (``pip install "ampere[sbi]"`` unlocks either): it still runs and is the
   fullest worked example of an SED fit in this repository's legacy
   notebooks, but it gains nothing new and the two APIs do not interoperate
   — see :doc:`migrating`'s ``SBI_SNPE`` row. :class:`~ampere.inference.SBIEngine`
   is where new SBI work lands, over ``sbi`` 0.27 (NPE, NLE, NRE and TMNRE),
   consuming the same :class:`~ampere.core.FittingProblem` every other v2
   engine does.


Embedding networks for automatic summary statistics
-----------------------------------------------------

Neural posterior/likelihood/ratio estimation needs a fixed-size summary of the observed data to condition the network on. When your data is high dimensional, comparing raw simulated and observed data directly makes training slow and can make a poor summary; a network that learns the summary statistics — an embedding network — usually does better than one hand-picked, and generalises across problems a hand-picked one would not.

:class:`~ampere.inference.SBIEngine` supports this two ways, chosen with ``layout=``/``embedding=``. The default, ``layout="flat"``, is a fixed-size vector — each dataset's observed values, masked samples dropped, concatenated in ``datasets`` order — which is the simple, sufficient choice for a single fitting problem, since the data layout never changes between simulation and inference. Where the layout itself can vary between simulated draws — irregular sampling, missing data, or amortising across differently-configured instruments — ``layout="set"`` packs every dataset through the coordinate–value–mask **encoding** (:class:`~ampere.core.encoding.EncodingLayout`) and reads it with a masked permutation-invariant (``embedding="set"``) or attention (``embedding="transformer"``) network, so the summary is learned over data whose shape is not fixed in advance. Both are ``sbi`` 0.27 nets behind an ampere wrapper that handles the mask column and pools every retained token, at a default output width of ``max(2 * free_size, 32)`` — a default, not a finding: which width and readout suit which combination of data, model and structure is a deferred study (``DEVELOPMENT_PLAN.md`` §6). :doc:`sbi`'s "Embedding choices" section works through both, and the legacy dict-based embedding vocabulary this section used to describe is documented, under the legacy warning, in :doc:`notebooks/Embedding_nets`.

Parallelised evaluation
-----------------------

There is deliberately **no multiprocessing pool** in :mod:`ampere.inference`.
The failure history a run records is per-process, so a pooled run would leave every worker with its own counts and the run's aggregate silently incomplete; aggregating them properly is real work rather than a keyword argument, and it is not done yet.

Two things do already parallelise the evaluation, and neither needs a pool.
On the torch backend, ``log_prob_unconstrained_batched`` evaluates a whole stack of parameter vectors in one ``torch.func.vmap`` call, for any problem whose parts all report ``BATCHABLE``.
On either differentiable backend, the gradient itself is what buys the speed: a NUTS draw guided by a gradient is worth many blind ensemble evaluations, and :doc:`m2_misspecification` measures both sides of that trade.
Failing those, you can parallelise your model itself if that makes sense.


Defining new data types
-----------------------

AMPERE packages a selection of container types suitable for the most common astronomical datasets - spectra, photometry, images, cubes, time series, interferometric visibilities and closure phases - but these might not always cover what you need.

In v2 a container **kind** is a registered, extensible thing rather than a fixed list: :func:`ampere.results.register_kind` adds yours, :func:`~ampere.results.registered_kinds` lists what is known, and the registration is what lets a run store, hash and reload data of your kind alongside everything else.
:doc:`interferometry` is now the worked guide, not a sketch — a kind-changing step, a complex container, and the two composition shapes ("two instruments on one channel" and, in :doc:`astrometry`, "one model, two channels") that ship — written as a template a new modality follows section by section. The design sketches in ``docs/design/modalities/`` cover what has not shipped yet - IFU cubes, awkward instruments, and the joint 2-vector GP over correlated channels deferred to Phase 5 - and are the place to start for those.
