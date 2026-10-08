Frequently Asked Questions
==========================



Models
------



Data
----



Inference
---------

**summary says my R-hat is 1.5 — what now?** Ask for the verdict:
:func:`ampere.results.check_convergence` holds every posterior variable to an
R-hat below 1.05 and a bulk effective sample size of at least 100 (both
keywords you can change) and returns the failing variables with their numbers
and a remedy in words: more steps when only the ESS is short, more walkers and
more steps when most variables have not mixed, the optimiser's start if the run
began from the prior, a reparameterisation when one variable fails alone.
``emit`` and :func:`~ampere.results.summary` print the same verdict as a
:class:`~ampere.results.plots.ResultsWarning` whenever a chain-based run fails
it; a VI or SBI run, or a nested sampler's resample, is not a chain and is
never warned about. The usual first cause is the start, the next answer's
subject.

**Which start does run use?** For :class:`~ampere.inference.EmceeEngine` and
:class:`~ampere.inference.ZeusEngine` the default is the optimiser's mode:
``run()`` calls ``optimise(problem, method="scipy", starts=1)`` once and puts
the walkers in a tight ball around it (:doc:`optimisers`). ``initial="prior"``
starts them from prior draws instead — the default in 1.0.0b1, reproduced
bit for bit — and ``initial=optimise(problem)`` hands over the full eight-start
route, or any :class:`~ampere.results.Optimum` you already have. When the
optimiser finds no start it can score, or its mode sits on a prior bound (a
ball there is degenerate and every walker stays on the bound), the run falls
back to the prior with a :class:`~ampere.inference.DefaultStartWarning` that
names the bound and the remedies. The run records which in
``ampere_start_kind`` (``"optimum"``, ``"prior"`` or ``"supplied"``) and any
fallback in ``ampere_start_fallback``. NUTS, VI and the nested samplers keep
their own starts.

**How many walkers and steps?** The default walker count is four per free
parameter, at least eight; more walkers buy a better-sampled ensemble, more
steps a longer chain. From the optimiser's start a few-parameter SED fit
converges in the low thousands of steps (the
:doc:`notebooks/quickstart` runs a deliberately small 24 walkers by 400);
let :func:`~ampere.results.check_convergence` decide rather than a rule of
thumb, and double the steps until it passes.

**When do I move to NUTS?** When the problem has many free parameters. An
ensemble sampler's mixing time grows with the dimension, and emcee and zeus
warn once per run (:class:`~ampere.inference.EnsembleSizeWarning`) from
sixteen free coordinates: a twenty-member :doc:`population <population>` took
eight minutes on the numpy path to reach R-hat 1.23, where
:class:`~ampere.inference.NUTSEngine` on a torch or jax problem fits fifty.
Composing the problem from a native backend is the move; :doc:`overview`'s
"Realisation, and the engines" says which engines run where.



Other
-----
