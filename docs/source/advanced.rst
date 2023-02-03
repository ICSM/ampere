Advanced usage
==============

The previous guides described how to do fairly simple things with AMPERE like defining a model. However, you will often want to do something more complex than that. These tutorials are intended to help you achieve that.

Model comparison
----------------

Often our objective with inference is not only to determine which distribution of parameter values is supported by our observations, but instead to determine which model (out of some set of 2 or more models) is most probable, given our observations.
This is the *model selection* problem, rather than parameter estimation.

AMPERE natively supports model selection with Nested Sampling - simply use the interface to Dynesty to do your inference instead of any of the other approaches, and you will get an estimate of the *model evidence* at the end of run.
Do this for multiple different models, and compare the logarithm of the evidence at the end - whichever model has the highest evidence is the preferred model (assuming all models are equally probable *a priori*).
This is very convenient, since it automatically penalises models with different numbers of parameters, and doesn't require that models be nested.
However, this is only well-justified if the model parameters have physical meaning.

At present, if you want to do model comparison for other approaches than Nested Sampling, you will have to roll your own.
For computing the evidence from *emcee* or *zeus* results, the `harmonic <https://astro-informatics.github.io/harmonic/index.html>`_ package may be effective.
Alternatively you can do model comparison by exporting AMPERE's chains to `arviz <https://www.arviz.org/en/latest/>`_ and using their WAIC and/or LOO-PIT methods.
These methods have the advantage of giving meaningful results even if the parameters are ad hoc without real physical meaning.

More detailed tutorials will be available soon!


Very slow models
----------------

AMPERE allows you to use a wide variety of models to interpret your data, which may include models which take a very long time to compute.
In such cases, Neural Posterior Estimation (NPE) is probably a good bet!

AMPERE uses `sbi <https://>`_ under the hood to do NPE. You can see a few examples in the tutorials, but a more complete guide will appear here in the future. 


Parallelised evaluation
-----------------------

Automatic parallel evaluation of paramater sets will be implemented in AMPERE soon! In the meantime, you can parallelise your model itself if that makes sense.


Defining new data types
-----------------------

AMPERE packages a selection of data objects suitable for the most common astronomical datasets - photometry and spectra (with more to come).
However, these might not always cover what you need, as you might need to specialise these for the features of your specific dataset.
In that case, you will want to create your own data classes that encapsulate these features. A guide will appear here in future!
