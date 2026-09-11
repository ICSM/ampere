How does AMPERE work?
=====================



Diverse datasets
----------------

Astronomy increasingly faces the problem of combining multiple datasets to answer one question.
This results in large, complex datasets that are also rich in information.


Model mis-specification
-----------------------

As we build more powerful and precise instruments, our datasets also get richer.
This often results in features that our models are unable to explain (progress!) but this 'model mis-specification' leads to problems in inference;
unexplainable features create structure in the residuals, which in turn tend to distract the optimiser and distort parameter estimates.
For example, if your data is a spectrum that contains emission lines, but your model only predicts continuum, the posterior parameter estimates
will tend to slightly overestimate the continuum flux to compromise between line and continuum levels.

Flexible noise models
---------------------

This problem can be at least partially mitigated by realising that the structure residuals created by a mis-specified model are indistinguishable from correlated noise in the data.
Therefore, by explicitly modelling an additional component of correlated noise in the dataset and marginalising over it in our inference, we can add a certain amount of resistance
to mis-specification, and get less-biased posterior estimates of the parameters we are really interested in.

.. note::

   There is no *mis-specification proof* model here, just like there is no such thing as an earthquake-proof building. A big enough mis-specification
   will still break your inference, in this case by increasing the correlated noise component to arbitrary large levels so that the likelihood remains
   reasonable, just like a sufficiently energetic earthquake will bring down any building.

AMPERE achieves this by modelling the covariance matrix of the data as the sum of the identity matrix and an additional correlation matrix [#1]_, all multiplied by the variances of the data.
The correlation matrix is composed from a set of *stationary* kernel functions, as in a Gaussian Process.
This provides a high degree of flexibility in the types of noise that get generated, without adding large numbers of additional free parameters.

The default kernel is Matérn-3/2 rather than the squared exponential more familiar from the Gaussian-process literature, for two independent reasons.
A Matérn-3/2 sample path is once differentiable rather than analytic, which is a better description of a real model deficiency than infinite smoothness.
And it is *exactly* quasiseparable, which is what makes :class:`ampere.core.QuasisepGP` an **exact** O(N) solve rather than an approximation - the difference between fitting two hundred points and fitting twenty thousand.

A default is not a restriction.
:class:`~ampere.core.Matern12` and :class:`~ampere.core.Matern52` sit either side of it in smoothness; :class:`~ampere.core.SHO` is a damped oscillator, which is the right shape for a residual with a *period* - interference fringing, an instrumental ripple - and :class:`~ampere.core.RotationTerm` is the pair of them celerite2 uses for a non-sinusoidal one.
:class:`~ampere.core.Sum` and :class:`~ampere.core.Product` compose them, and a sum of quasiseparable terms is still quasiseparable, so a broad Matérn plus a narrow oscillator still costs O(N).
:class:`~ampere.core.SpectralMixture` is that sum with free frequencies.
Your own kernel works on the dense solver the moment you write its covariance function, and reaches the O(N) path once you register its celerite representation with :func:`~ampere.core.register_quasiseparable_term`.

A kernel may also act on a *subset* of a container's coordinates, named with ``axes=``.
That matters wherever a residual is correlated in two different ways at once: a missing patch of sky with a spectral profile is smooth across spatial frequency and sharp across wavelength, and its covariance is a :class:`~ampere.core.Product` of one kernel on ``("u", "v")`` and another on ``("spectral_axis",)``.
A single isotropic length scale over axes in different units is meaningless, and ampere refuses it rather than being silently wrong.

:doc:`m2_misspecification` measures all of this on a controlled problem, and :doc:`overview` shows how to compose it.

.. rubric:: Footnotes

.. [#1] This effectively separates the noise into a correlated and an uncorrelated component.
