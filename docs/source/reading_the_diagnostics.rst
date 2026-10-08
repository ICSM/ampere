Reading the diagnostics
=======================

You have fitted something, perhaps with a flexible likelihood, and ampere has
drawn you a residual plot, a ribbon labelled "GP localisation" or a
posterior-predictive histogram. This page says what each of them shows, the one
sentence each can support, the one it cannot, and what you would do next.
:doc:`concept` explains why the Gaussian process (GP) is there; this page is
about reading what it did.

**Where the figures and numbers come from.** The figures are built when the
documentation is, by the scripts in ``docs/source/plots/``, each of which
names in its docstring the driver it reproduces. They run at a **docs
budget**, short enough that the whole page builds in a couple of minutes, so
a figure's chain is shorter than the one behind the numbers in the prose. The
prose quotes the full-budget numbers, from the tables of
:doc:`m2_misspecification` (the "M2 page") and, for the first case, from the
user-journeys memo's Appendix B.2 ("B.2"); the captions say which is which.

The gross case: what "absorbing a lot of power" looks like
-----------------------------------------------------------

:doc:`concept` says that there is no misspecification-proof model, "just like
there is no such thing as an earthquake-proof building". This is the
earthquake. The data are a *Spitzer*/IRS spectrum of the quasar PG 1011-040
(``examples/test_data/cassis_yaaar_spcfw_14191360t.fits``, 360 points from 5.2
to 37.4 µm), and the model is a single power law,
:math:`F_\nu = \mathrm{norm}\cdot\lambda^{\mathrm{index}}`. The spectrum has
silicate emission at 10 and 18 µm, which a power law cannot have, so the model
is wrong in a way that is hard to miss.

Fitted with independent noise, the power law reports an index of
1.111 ± 0.003 (B.2). That is a precise answer to a question the data do not
support, and the residuals say so:

.. plot:: plots/gross_residuals.py
   :alt: Two panels. Top, the standardised residuals of a power-law fit to an
      infrared spectrum run smoothly from minus ten to plus ten standard
      deviations: negative below 9 microns, positive from 10 to 24, negative
      again beyond 25. Bottom, the binned autocorrelation of those residuals
      is between 0.6 and 0.85 at every separation up to 0.5 microns, and the
      panel's title reports a whiteness statistic that rejects white noise at
      the test's floor, p = 0.005.
   :caption: The power law's standardised residuals, with the whiteness test
      (``plot_residuals``). B.2's full-budget run: Q = 1490, p = 0.005, the
      residuals "thread the silicate hump at ±10σ". This figure is the same fit
      at the docs budget (24 walkers, 300 steps).

*It supports:* these residuals are not white, they are ten standard deviations
from zero, and the model leaves the silicate bumps in the data. *It does not
support:* the index being wrong by any particular amount, or anything about
what the structure is.

Now give the same power law a GP, with a **permissive** prior: a half-normal on
the amplitude of scale 0.05 Jy and one on the length scale of scale 5 µm.
(The median flux of this spectrum is 0.059 Jy, so that amplitude prior says "let
the GP carry whatever it must".) The GP does exactly that:

.. plot:: plots/gross_localisation.py
   :alt: The conditioned Gaussian-process mean of a power-law fit with a
      permissive kernel prior, plotted against wavelength from 5 to 37 microns
      with its 68 per cent band. It is a smooth curve of several microns'
      length scale that peaks near 20 microns and is negative below 9 microns,
      the outline of the silicate emission the power law lacks. The caveat text
      is printed below the axis.
   :caption: The conditioned GP mean from the same fit with the permissive
      kernel prior (``plot_gp_localisation``), at the docs budget. B.2's
      full-budget runs: a GP of 0.033 Jy at a length scale of 9.4 µm from the
      optimiser, 0.022 Jy at 7.1 µm from the prior (the chains were not
      converged, R-hat 1.4); the conditioned mean is the 10 and 18 µm silicate
      emission, +0.02 Jy at 18-20 µm, with a trough below 9 µm. The caveat
      quoted under the axis is the plot's own.

Read B.2's numbers against the independent fit's. The GP carries about a third
of the flux, at a length scale of several microns (it is not noise in any
sense; it is the spectrum's shape). The power law's index loses two orders of
magnitude of precision, from ±0.003 to about ±1 (0.9 ± 1.0, from the
optimiser). And under the permissive prior the normalisation runs to its bound
while the GP carries the lot: the optimiser's answer sits on ``model.norm``'s
edge, which is why ampere warns that the optimum is at a bound. A GP this
large is not a robust correction to your model. It is the model.

*The decision here is to fix the model*, in this case by adding the silicate
features (a dust component, a template), and not to report the GP's amplitude as
a result. The GP told you where to look, and what it found at 9.4 and 18 µm is
what the model is missing. Ampere cannot make a wrong model right; it can show
you how wrong it is.

The mild case, where the GP is doing its job, is the rest of this page; the
four plots below are read on the M2 study's scenarios, a four-parameter toy
(a linear continuum with two Gaussian absorption lines, 200 points, 1 % noise)
into which a known deficiency is injected. :doc:`m2_misspecification` describes
it. (A second mild case, a real spectrum, joins this page with the NGC 6302
notebook.)

The four plots
--------------

The residual plot and its whiteness test
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``plot_residuals`` draws the signed standardised residuals ((data - model) /
uncertainty) of a fit, and below them the autocorrelation of the residuals
binned by separation, with the statistic *Q* and the permutation *p* of
:func:`~ampere.results.residual_whiteness` in its title. It is for fits with
**independent** noise. A GP fit's residuals are whitened by construction, so the
plot warns if you give it one, and points at the localisation plot instead.

First, what passing looks like. The control (``none``: no deficiency injected):

.. plot:: plots/m2_residuals_none.py
   :alt: Standardised residuals of a correct model scatter about zero, within
      about plus or minus two, with no pattern. The binned autocorrelation
      below is small, at most 0.15 in magnitude and changing sign, and the
      panel title reports a p-value of about 0.2, consistent with white noise.
   :caption: ``none``, standard likelihood, at the docs budget. Full budget (M2
      page, "The diagnostics agree"): p = 0.185.

*It supports:* these residuals are consistent with independent noise of the
stated size. *It does not support:* the model being right, because a test that
fails to reject is not an acceptance, and a model error smaller than the noise
or one that varies faster than the sampling is invisible to it.

Now the mild scenario, in which 2.5 % sinusoidal fringing (a period of
0.0028 µm) is injected that the model lacks:

.. plot:: plots/m2_residuals_mild.py
   :alt: Standardised residuals of the mild scenario oscillate regularly
      between about minus three and plus three, about eleven cycles across the band.
      The binned autocorrelation below falls from about 0.8 at the shortest
      separation to about minus 0.65 at the longest, a sinusoid's signature,
      and the panel title reports a p-value at its floor of 0.005.
   :caption: ``mild``, standard likelihood, at the docs budget. Full budget (M2
      page): p = 0.005, the floor of a 199-permutation test, Q = 376.

*It supports:* the fit left structure behind, and the model, noise budget or
both are incomplete. *It does not support:* what the structure is, nor that the
parameters are wrong (they are biased here, but the plot does not measure
that). **The decision:** this is the signal to do something. Either fix the
model, or let a GP take it (the next plot); and if the structure is real but
you cannot model it, carry on with the flexible likelihood instead of the
independent one.

The GP localisation plot, and its caveat
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With a flexible likelihood the question changes. The GP has absorbed whatever
it could, so the interesting thing is *where* it did.
``plot_gp_localisation`` draws the conditioned GP mean (what the GP says the
model is missing, given the data) with its posterior band on the data's own
axis, and the number behind it is :func:`~ampere.results.gp_localisation_score`.

.. plot:: plots/m2_localisation_mild.py
   :alt: The conditioned Gaussian-process mean for the mild scenario is a small
      oscillation across the whole band, at the same period as the injected
      fringing, with a band that excludes zero at many wavelengths. The caveat
      text is printed below the axis.
   :caption: ``mild``, flexible likelihood, at the docs budget. Full budget (M2
      page): peak score 3.70, GP amplitude 0.0248 ± 0.0050 Jy, about 2.5 % of a
      continuum of about 1 Jy, at a length scale of 0.00089 µm.

The GP's amplitude is the size of the fringing the model lacks, and nothing
else: this is the case that the flexible likelihood is for, and you can read it
off the plot. *It supports:* the GP absorbed structure of about this size,
at about this wavelength scale, and where it did. *It does not support:* "the
model is wrong". The same shape is also what an underestimated noise budget
gives, and the GP cannot tell them apart.

When one place stands out:

.. plot:: plots/m2_localisation_sharp.py
   :alt: The conditioned Gaussian-process mean for the strong_sharp scenario
      is close to zero across the band except for one narrow peak of about 0.12
      Jy at 0.863 microns, where the injected emission line sits, whose band
      excludes zero. The caveat text is printed below the axis.
   :caption: ``strong_sharp``, flexible likelihood, at the docs budget. Full
      budget (M2 page): peak score 12.8 at 0.86295 µm; GP amplitude 0.0215 ±
      0.0034 Jy at a length scale of 0.00070 µm.

*It supports:* "the GP absorbed structure at 0.863 µm; the model is missing a
feature there", which is what an unmodelled 12 % emission line is, at 0.8630 µm.
*It does not support:* "the model is wrong" in general, nor that the feature is
an emission line (it could be a blend, a calibration artefact or a correlated
noise source at that wavelength). **The decision:** look at the spectrum there,
and if it is a real feature, put it in the model. If you cannot, accept the
GP's correction and *report its amplitude* with the result, because that is the
size of what you did not model.

The plot prints a caveat under the axis, in the figure and in its metadata
(``show_caveat=False`` removes the drawing but keeps the text). Its clauses,
in turn:

* "A large fitted GP amplitude localises where the model is deficient; it does
  not say why." The plot says *where*, never *why*.
* "Model error, an underestimated noise budget and a genuinely correlated
  astrophysical process are all consistent with the same posterior shape." Three
  different problems, with three different fixes (the model, the uncertainties,
  and nothing at all, because the correlation is part of the sky), look the same
  here.
* "A large amplitude at a short length-scale may equally mean that the kernel's
  smooth global component is under-amplitude and compensating locally." A
  spike can be a real narrow feature, or a smooth model deficiency that the
  kernel was too stiff to describe in one piece. Compare the length scale with
  the width of the feature you see.

The posterior-predictive check
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``plot_posterior_predictive`` asks a different question: if the fitted model
were true, would data like mine be typical? Panel one shows the observations
against the replicate median and band; panel two the distribution, over
replicates, of a statistic *T* (by default the sum of squared standardised
residuals), with the observed *T* marked and the p-value.

.. plot:: plots/m2_posterior_predictive_sharp.py
   :alt: Two panels for the strong_sharp scenario's standard fit. On the left,
      the observed spectrum follows the replicate band everywhere except a
      group of points near 0.863 microns that rise to 1.17 Jy against a band at
      about 1.0, where the emission line was injected. On the right, the histogram of the
      replicate statistic is centred near 200 and the observed value, near
      800, lies far to its right, giving a p-value of zero.
   :caption: ``strong_sharp``, standard likelihood, at the docs budget. The same
      fit's whiteness statistic at full budget is Q = 272, p = 0.005 (M2 page,
      "The diagnostics agree").

*It supports:* the model, with the noise you gave it, does not reproduce this
data set as a whole: a very small p-value here means the observed misfit is
well outside what the model would produce. *It does not support:* where the
problem is (use the localisation plot or the residuals) or which part of the
model is at fault, and a p-value near 0.5 does not say the model is right, only
that this statistic did not notice. **The decision:** a clear failure sends you
to the localisation plot or the residuals; do not tune the model to the
statistic.

The anomaly score
~~~~~~~~~~~~~~~~~

``plot_anomaly_score`` draws one number per coordinate on a single axis, so you
can compare it between data sets. For a fit with a GP it is the size of the
conditioned GP mean in units of its own uncertainty; before a fit, the same
renderer draws the pre-fit RHMF map.

.. plot:: plots/m2_anomaly_sharp.py
   :alt: A line plot of the anomaly score against wavelength for the
      strong_sharp scenario: below about 2 everywhere except for one tall peak
      of about 14 at 0.863 microns, with the provenance legend and the
      interpretation notes drawn on the figure.
   :caption: ``strong_sharp``, flexible likelihood, at the docs budget. Full
      budget (M2 page): peak score 12.8 at 0.86295 µm. The control, ``none``,
      peaks at 0.28; ``mild`` at 3.70; ``strong_smooth`` at 2.08.

*It supports:* where, relative to its uncertainty, the GP is working hardest.
*It does not support:* a significance in the usual sense (it is a ratio, not a
test, and it inherits the caveat above), or a comparison between fits with different kernels or priors. Note what the M2 table says
about the scale: the 7 % smooth fringing scores only 2.08 and the 2.5 % mild
fringing 3.70, while the 12 % sharp line scores 12.8, so the score ranks
*places* within one fit and is not a measure of how large the deficiency is. **The decision:** use the score to decide where to
look, not whether to worry; the amplitude and the length scale decide that.

When the GP has several things to absorb
----------------------------------------

One stationary length scale suits one scale of deficiency. The M2 study's
``many_lines`` scenario has two: five narrow lines (a correlation length of
0.00035 µm) over a smooth continuum error (a hundred times longer):

.. plot:: plots/m2_many_lines.py
   :alt: Two panels sharing a wavelength axis. Top, the injected deviation,
      five narrow lines of up to 0.11 Jy between 0.860 and 0.870 microns on a
      smooth trend of about 0.02 Jy, with the observed residuals. Bottom, the
      conditioned Gaussian-process means of three kernels, which at this short
      budget all recover the five lines and differ mainly in their bands, the
      warped kernel's being much the widest.
   :caption: The ``many_lines`` scenario with three kernels, at the docs budget
      (26 walkers, 260 steps), where the three conditioned means are hard to
      tell apart by eye: what separates the kernels is the parameters they
      recover, and that is in the M2 page's table, not this figure. Full budget
      (M2 page, "When one length scale is not enough"): worst offsets 2.84
      widths for the single Matern-3/2, 0.92 for the warped kernel and 0.78 for
      the sum of two terms.

A **sum** of terms gives each scale a component, and also gives the model the
freedom to add a component the data do not need. The M2 page's shrinkage
section ("A sum of noise components needs a sparsity guard") measures what the
horseshoe prior does about that on two nearly degenerate terms fitted to a
one-component truth:

============ ==================== ==================== ====================
prior        larger amplitude     median min/max       P(min/max < 0.1)
============ ==================== ==================== ====================
flat         0.0118 Jy            0.254                0.257
horseshoe    0.0077 Jy            0.087                0.537
============ ==================== ==================== ====================

*It supports:* the redundant component's share falls by a factor of 2.9 under
the horseshoe while the real one survives. *It does not support:* that every
component the horseshoe keeps is real; it supports that the ones the data do
not need are discouraged. **The decision:** for any ``Sum`` of noise terms, put
the shrinkage prior on it (``ampere.core.with_shrinkage``) and read the
components' amplitudes as upper limits on what each term absorbed, not as
detections. (The figure is the study's ``many_lines`` figure, not a plot of the
shrinkage table, whose numbers the page quotes.)

Choosing the kernel's priors
----------------------------

The GP can only absorb what its priors allow. They are the one place where
*you* tell it how much of the data to give up, and they are best chosen against
what the model itself knows.

**The amplitude, as a fraction of the flux.** A half-normal whose scale is a few
per cent of the median flux says "the model should be nearly right; a GP bigger
than that is a surprise", and a large fitted amplitude then stands out. A scale
comparable to the flux says "let the GP carry whatever it must", and it will:
it is the permissive prior of the gross case above (0.05 Jy against a median
flux of 0.059 Jy). A good default is the size of the deficiency you would be
willing to tolerate without comment: the M2 study's fringing was 2.5 % of the
flux and the GP found 0.0248 Jy on a ~1 Jy continuum, which a prior at a few
per cent of the flux contains with room to spare. If your prior is much wider
than the amplitude you find, the amplitude is not constrained by the prior, and
the fit will tell you so; if the fitted amplitude is near the *prior's* scale
or hugs a bound, the GP is doing the model's job, which is the gross case.

**The length scale, against the features' widths and the data's span.** Three
regimes.

* *Below the grid spacing*, the GP is white noise: it can only be correlated
  between points closer than a length scale, and there are none. It absorbs
  nothing structured and the diagnostics then say nothing. The quickstart's
  half-normal of scale 0.01 µm is on an IRS grid whose spacing here is 0.03 to
  0.17 µm (median 0.06 µm), so that prior is below the sampling over the whole
  band. That is the case in the quickstart, and it is harmless there only
  because the quickstart's model fits, so there is nothing to absorb; to catch
  structure at the 10 µm silicate feature a scale of a micron or so would do.
* *Above the span*, a length scale longer than the wavelength range is a
  continuum offset, degenerate with the model's normalisation, and the fit will
  trade one against the other. The gross case's permissive prior (5 µm on a
  32 µm span) is already near this, which is part of why the GP and the
  normalisation traded freely.
* *In between*, the length scale should be about the width of the structure you
  expect the model to miss. The M2 study's fringing has a 0.0028 µm period and
  the fitted length scale settled at a fraction of it (0.00089 µm); the 12 %
  line had a σ of 0.35 nm (0.00035 µm) and the fitted length scale was 0.00070 µm.

A worked example on this page's IRS spectrum: a grid spacing of about 0.06 µm
and a span of 32 µm give a range from 0.06 to 32 µm. The silicate emission
that the power law misses spans several microns (the gross case's GP found
length scales of 7 to 9 µm), so a half-normal of scale 1 to 3 µm covers the
features of that width and everything narrower, with an amplitude scale of about 0.003 Jy (5 % of the
median flux of 0.059 Jy). If the fit then reports an amplitude near 0.03 Jy,
the model is missing something large, and you should go back and fix it. A
spectrum whose structure spans several scales is the case for a sum of two
kernels, or for a warped kernel, whose input map is built from the data's own
quantiles with :func:`~ampere.core.quantile_knots`.

What none of the plots can tell you
-----------------------------------

All four plots say *where* and *how much*, never *why*. An underestimated noise
budget, a model error and a real, correlated astrophysical process all produce
the same GP, and no plot here separates them: that takes physics, a second
instrument or a better calibration. What you can always do is report the GP's
amplitude and length scale next to your result, so that a reader knows how much
the answer leaned on the correction. For where to go next, :doc:`kernels` lists
the kernels, :doc:`m2_misspecification` measures what the correction costs and
buys, and :doc:`concept` explains the idea.
