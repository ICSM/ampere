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


Very slow models
----------------

AMPERE allows you to use a wide variety of models to interpret your data, which may include models which take a very long time to compute.
In such cases, Neural Posterior Estimation (NPE) is probably a good bet!

AMPERE uses `sbi <https://www.mackelab.org/sbi/>`_ under the hood to do NPE. You can see a few examples in the tutorials, but a more complete guide will appear here in the future.

.. note::

   NPE currently lives in the **legacy** half of the package
   (``ampere.infer.sbi``, ``pip install "ampere[sbi]"``). The v2 SBI layer is
   Phase 3 of the redesign and will live in :mod:`ampere.inference` beside the
   other engines; this section describes what exists today.


Embedding Networks for automatic summary statistics with NPE
------------------------------------------------------------

NPE is a powerful tool for speeding up inference with models where the likelihood is difficult to evaluate. 
However, *because* the likelihood is difficult to evaluate, we can find ourselves simply comparing the simulated data to the real ones.
When your data is high dimensional, this can make the comparison difficult, and even worse, it makes training the neural network for the posterior _very_ slow.

In these cases, it is better to define some summary statistics that can reduce the dimensionality of the data with minimal loss of useful information. 
However, in the case of astronomical data it can be difficult to define a good statistics that are also easy to transfer to other problems.
Hence, we can use a neural network to learn the best summary statistics for our problem, and then use these to train the NPE network.
The interface provided by SBI is exposed and instructions for how to use it can be found in :doc:`notebooks/Embedding_nets`.

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

AMPERE packages a selection of container types suitable for the most common astronomical datasets - spectra, photometry, images, cubes and visibilities - but these might not always cover what you need.

In v2 a container **kind** is a registered, extensible thing rather than a fixed list: :func:`ampere.results.register_kind` adds yours, :func:`~ampere.results.registered_kinds` lists what is known, and the registration is what lets a run store, hash and reload data of your kind alongside everything else.
The design sketches in ``docs/design/modalities/`` work several new modalities through the contracts end to end - interferometric visibilities, IFU cubes, astrometric time series, awkward instruments - and are the place to start.
A guide will appear here in future!
