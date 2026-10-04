Advanced usage
==============

The previous guides described how to do fairly simple things with AMPERE like defining a model. However, you will often want to do something more complex than that. These tutorials are intended to help you achieve that.

Model comparison
----------------

Often our objective with inference is not only to determine which distribution of parameter values is supported by our observations, but instead to determine which model (out of some set of 2 or more models) is most probable, given our observations.
This is the *model selection* problem, rather than parameter estimation.

AMPERE natively supports model selection with Nested Sampling - use one of the three nested samplers, :class:`~ampere.inference.DynestyEngine`, :class:`~ampere.inference.NautilusEngine` or :class:`~ampere.inference.UltranestEngine`, instead of any of the other engines, and you will get an estimate of the *model evidence* at the end of the run.
Whichever you use, the estimate is recorded under the same three names — ``ampere_log_evidence``, ``ampere_log_evidence_err`` and ``ampere_evidence_method`` — so comparing two archived runs never requires knowing which sampler produced them.
See :doc:`ampere.inference` for which of the three to reach for; dynesty is in the base install and the other two are behind an extra each.
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
of tree, with no change to ``ampere.core``. :doc:`kernels` §6 walks through
the registration with the same trivial worked example
``tests/core/test_kernels.py`` uses to prove the route.

Shrinkage priors on kernel amplitudes
-------------------------------------

A ``Sum`` of noise kernels is as easy to write as one kernel, and that freedom needs a guard: a component the data do not need will otherwise be fitted to whatever the others leave over. :func:`~ampere.core.shrinkage_horseshoe` declares the guard — one global scale shared by every component, one local scale per component under it, and each component's amplitude under its own local scale — and :func:`~ampere.core.with_shrinkage` puts it on a kernel without changing anything else about it. :doc:`kernels` §1 is the reference for both, with the declaration written out; here is the short form, which names the amplitudes of the two components by their labels:

.. code-block:: python

    kernel = with_shrinkage(
        Sum(Matern32(...), Matern32(...), labels=("broad", "narrow")),
        shrinkage_horseshoe(("broad.amplitude", "narrow.amplitude")),
    )

The shared global scale is what makes this a *sparsity* prior rather than several independent shrinkage priors: a component the data insist on pulls the global scale up, and every other component is then shrunk against that same scale. Reach for it whenever you add a second or third term to a ``Sum`` — a broad term for the smooth residual plus a narrow one for the structure the model has not got, say — and you do not know in advance that the data support all of them. It is the recommended prior for any ``Sum`` of noise components, and the same hierarchical shape is what the warped kernels' default priors give their knots (:doc:`kernels` §3).

It is not the right tool when the question is *which structure*, rather than *how many components*. The prior shrinks amplitudes towards zero; it does not tell two components which of them should survive, so if the data cannot distinguish the terms (two nearly degenerate Matérn terms are the extreme case), the prior decides the split, and an answer that depends on the prior is a statement about the prior. Nor does it help when a single component has the wrong *shape*: shrinking its amplitude only hides the misfit, and the remedy is a different kernel (:doc:`kernels` §3's warps, for a deviation whose length scale changes across the band). For a lone component with a sensible amplitude prior of its own there is nothing to share a scale across, and the three levels only add a funnel to the posterior, which an ensemble sampler explores unevenly (:doc:`m2_misspecification` makes the same point about its own assertions).

What the M2 study measured is on :doc:`m2_misspecification`, under "A sum of noise components needs a sparsity guard". Two nearly degenerate Matérn terms were fitted to a spectrum whose deviation has exactly one smooth component, so the likelihood pins their total and says almost nothing about how it is divided. Under a flat prior the larger amplitude had a posterior median of 0.0118 Jy and the median ratio of the smaller amplitude to the larger was 0.254, with 0.257 of the posterior mass below a ratio of 0.1; under the horseshoe the figures were 0.0077 Jy, 0.087 and 0.537. The redundant component's share fell by a factor of 2.9 and the posterior mass at "the fit chose one component" rose by 2.1, while the component the truth does have survived — which is the check that matters, since a prior that shrank everything would move both numbers together and be useless. ``python -m examples.m2_misspecification.many_lines --shrinkage`` prints that table.

**What is and is not the regularised horseshoe.** The default ``tail="regularised"`` is Piironen and Vehtari's idea in the form this version of ampere can declare, not their exact construction. Their slab bounds the local scale through :math:`\tilde\lambda_j^2 = c^2\lambda_j^2/(c^2 + \tau^2\lambda_j^2)`, which is a deterministic function of two sampled parameters, and the parameter contract declares parameters and priors rather than deterministic nodes. What it declares instead is a gamma-tailed local scale, ``gamma(a=1/2, scale=tau)``: the same spike at zero, where all the shrinkage comes from, with an exponential tail in place of the Cauchy one, which bounds the tails and calms the geometry as the slab is meant to. The exact slab waits on a ``Derived`` node (see :doc:`kernels` §1 and ``likelihoods.md`` §6 for the account); ``tail="cauchy"`` gives the plain horseshoe if you want the heavier tail.

Populations: fitting many objects together
------------------------------------------

Sometimes the objects are not independent: a hundred discs, each with its own
spectrum and its own fit, whose spectral indices you believe are drawn from
*one* distribution. Fitting them separately throws that belief away; fitting
them with a shared index asserts something much stronger than you believe.
The middle is a hierarchical, or population, model — each object keeps its own
value, and the values are draws from a parent distribution whose parameters
are fitted too.

:class:`~ampere.core.Population` is how you say that, and the point of it is
that you say it **once, where you compose the fit**, not inside each object's
model:

.. code-block:: python

   from ampere.core import FittingProblem, HierarchicalPrior, Parameter, Population
   import scipy.stats as st

   discs = Population(
       "discs",
       members=[Parameter("index", HierarchicalPrior("norm", {"loc": "mu", "scale": "sigma"}))],
       hyperpriors=[Parameter("mu", st.norm(-1.0, 1.0)),
                    Parameter("sigma", st.halfnorm(0.0, 1.0))],
       over=[f"disc{i}" for i in range(50)],
   )
   problem = FittingProblem(models, datasets, populations=[discs])

``models`` here is fifty ordinary models — library models, written by someone
who had never heard of your survey. None of them declares ``mu`` or ``sigma``,
and none of them is modified: the declaration replaces each one's own prior on
``index`` with its draw from the population, and hands it back a plain scalar
under the name it already uses. ``problem.parameters`` then holds ``discs.mu``,
``discs.sigma`` and one array-valued ``discs.index`` of fifty elements — 52
dimensions, three parameter objects.

That last sentence is the reason to prefer this over writing the hierarchy out
by hand. The fifty draws are **one** parameter, so they are one sample site,
and on the torch and jax backends they lower to a real ``pyro``/``numpyro``
plate; NUTS then fits the whole population jointly. The hand-written form —
fifty scalar parameters sharing a prior — is still available as
``Population(..., layout="flat")``, declares exactly the same density, and is
refused above 128 members, because the prior evaluation is linear in the
number of parameter *objects*: measured at 46 ms per evaluation at a hundred
members and 373 ms at a thousand, against about 1.3 ms for the same structure
as one plate. A user who has read that cost and wants the flat layout anyway
can turn the refusal into a loud warning with
``ampere.core.settings.override(flat_population_cap="warn")`` — a one-field
setting, not a ``Population`` argument, and it does not change the limit.

A hierarchical prior's default unconstraining bijection is inferred, not just
for location-scale families: a shape family whose support does not move with
its shape — ``gamma``, ``beta``, ``lognorm`` — infers correctly provided at
least one shape argument is given as a constant, while a family whose support
does move with its shape (``truncnorm``) still asks for an explicit
``Parameter(bijection=...)``.

Two companions are worth knowing about:

- :meth:`DatasetCollection.plate <ampere.core.DatasetCollection.plate>` builds
  the population from a list of datasets in order, which is the natural form
  when each member *is* a dataset — a spaxel of an IFU cube, an échelle order,
  one catalogue entry. The draws are routed to the models those datasets name.
- :func:`ampere.results.fit_population` gets population hyperparameters out of
  fits you have **already run**, by importance reweighting, with no joint fit
  at all. It cannot shrink the members towards each other the way a joint fit
  does, but it costs nothing beyond the archive and scales to any number of
  objects. Fit jointly when the members are being fitted anyway or when the
  shrinkage is part of the answer; reweight when the archive already exists.

Very slow models
----------------

AMPERE allows you to use a wide variety of models to interpret your data, which may include models which take a very long time to compute, or whose likelihood cannot be written down at all — a compiled radiative-transfer code behind a Python call, say.
In such cases, simulation-based inference is the right tool: :class:`~ampere.inference.SBIEngine` trains a neural network on simulated ``(theta, x)`` pairs, drawn through :meth:`~ampere.core.dataset.FittingProblem.simulate_many`, instead of consuming ``log_prob``.
It runs under a process pool with per-simulation timeouts and crash capture, so a slow or occasionally-crashing external simulator is the case it is built for rather than an edge case it tolerates.
:doc:`sbi` is the full tutorial — a black-box simulator fitted end to end, caching the trained posterior, truncated marginal ratio estimation for a tighter fit, and checking the result is calibrated; the runnable scripts are ``examples/sbi/``.

.. note::

   This is the **current** route. ``ampere.infer.sbi`` is legacy and frozen
   (``pip install "ampere-astro[sbi]"`` unlocks either): it still runs and is the
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
