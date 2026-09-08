WStat, and why ampere does not ship it
=======================================

X-ray astronomers fitting a source spectrum together with its background
almost always reach for the profiled Cash-with-background statistic — the
"W statistic" in XSPEC's terminology. It needs no background *model* at
all: the background rate in each channel is replaced by its
maximum-likelihood value, profiled out analytically, so the fit only ever
varies the source parameters.

**Ampere does not ship a WStat likelihood family**, and this is a
deliberate, documented decision (ruled 2026-09-03; see
``DEVELOPMENT_PLAN.md`` §2's "WStat / profile likelihoods" row and
``docs/design/modalities/awkward_instrument.md`` §5 and §9 question 3). A
profiled statistic is not a marginal likelihood — the background's own
uncertainty is replaced by a point estimate rather than integrated out — so
it does not produce a Bayesian posterior in the sense the rest of ampere
does. Ampere's documentation takes the opinionated line that the correct
approach is the **two-dataset Bayesian formulation**: source and background
as two datasets, sharing one background model, fit jointly.

That said, WStat is expressible in ampere, because the likelihood contract's
extension surface is exactly "a user adds a family without touching ampere"
(``docs/design/contracts/likelihoods.md`` §3). The worked example below
builds it as a user family, builds the two-dataset joint fit alongside it,
runs both, and compares them — pros, cons, and results, not false balance.

The full, runnable code is :download:`examples/wstat_comparison.py
<../../examples/wstat_comparison.py>`. It is included here as a literal
script rather than an executed notebook deliberately: this repository's
policy is that no run outputs, figures or other binary artefacts are
committed to git, and Sphinx's ``nbsphinx_allow_errors = True`` setting
means a notebook that failed to execute would not fail this documentation
build either — so a script the reader (and
``tests/examples/test_wstat_comparison.py``) can actually run, rather than a
page whose apparent success proves nothing, is the honest choice here. Run
it yourself with::

    python examples/wstat_comparison.py

which prints the comparison below, deterministically (every random draw —
the synthetic data and both engines' own streams — is seeded), in around
ten to twenty seconds on the reference backend.

The two routes
---------------

**Route 1 — the profiled statistic, as a user family**
(``ProfiledCashWithBackground``). One dataset (the source-region counts);
the background-region counts are data the family carries as its own
buffer, aligned per channel, together with the exposure/area ratio between
the two regions. For channel :math:`i`, the profiled background estimate is

.. math::

    \hat b_i = \frac{C_i + \sqrt{C_i^2 + 4\rho(\rho+1) B_i m_i}}{2\rho(\rho+1)},
    \qquad C_i = \rho (S_i + B_i) - (\rho + 1) m_i,

with :math:`\rho` the ratio, :math:`m_i` the model's predicted source
counts, and :math:`S_i`, :math:`B_i` the observed source- and
background-region counts. ``log_prob`` sums
``log Poisson(S_i; m_i + ratio * b_hat_i) + log Poisson(B_i; b_hat_i)``
over channels.

Three things this family gets right, on purpose, because they are exactly
what the likelihood contract asks of a user family carrying its own aligned
data:

* **Safe under masking.** The background buffer is excised with
  ``NoiseParams.retain`` — the boolean indicator over the *full*,
  unmasked containers every ampere noise model sets — so a masked source
  channel cannot silently pair against the wrong background channel.
* **Its ``sample()`` refuses**, and for a specific, structural reason
  rather than only the contract's general one: the profiled background
  estimate is a function of the very counts a draw would need to produce,
  so there is no forward generative model to sample from without already
  having the data it profiles over.
* **Its per-sample log-likelihood terms are flagged as degraded.** Each
  channel's term already used that channel's own counts to fix its own
  nuisance parameter, so — even though the terms sum exactly to
  ``log_prob``, unlike a GP's — they are not valid predictive densities.
  Handing them to ``arviz.loo``/``waic`` would be circular.

**Route 2 — the two-dataset Bayesian joint fit** (``XraySourceAndBackground``,
``build_joint_problem``). Source and background regions as two
``Dataset``\ s in one ``DatasetCollection``, observing two channels of one
model that shares the background parameter between them (`inference.md` §8's
pattern, specialised to the case where the shared physics needs no channel
of its own). Both channels are ordinary ``PoissonFamily`` likelihoods — no
profiling, no user family, nothing bespoke.

The comparison
---------------

Both routes are run with ``emcee`` on the same deterministic, seeded
synthetic data — 32 channels, a faint power-law source over a flat
background, chosen so several channels have very few or zero counts, which
is precisely the regime where a profiled nuisance and a marginalised one can
part company. A representative run's output::

    WStat (profiled, user family) vs the two-dataset Bayesian joint fit
    ======================================================================

    parameter          truth          wstat median [68%]          joint median [68%]
    src_norm           6.000        6.534 [5.727, 7.442]        6.592 [5.666, 7.488]
    src_index         -1.500     -1.701 [-1.911, -1.514]     -1.803 [-2.024, -1.614]

    On the source index, this run's two posteriors agree: each median falls inside
    the other approach's 68% interval; the two intervals are close to the same width.

**Recommendation: prefer the two-dataset Bayesian joint fit.** This holds
regardless of how closely the two routes' numbers happen to agree in any one
run. The joint fit is a proper marginal likelihood: its posterior is a
posterior, its per-sample terms are valid predictive densities for
``arviz.loo``/``waic``, and both regions can be simulated forward for
posterior-predictive checks. WStat can do none of that structurally, not as
a matter of a particular run's luck.

**The trade-off, stated plainly rather than as false balance.** The joint
fit pays for those guarantees with a background *model* — here, one
flat-shape parameter. If that shape assumption is wrong, the joint fit is
biased in a way WStat cannot be, because WStat's per-channel nuisance needs
no shape assumption at all; that is its one genuine advantage, and the
reason it remains a reasonable quick-look tool. The literature on profile
likelihoods in counting experiments (Cash 1979 and its X-ray descendants)
documents that profiling a background out, rather than marginalising it,
tends to bias and overstate the precision of the parameters of interest
specifically in the low-count regime this example uses. A single seeded run
— including the one printed above — cannot itself demonstrate that bias
reliably; only a repeated-trial coverage study can, which this fast worked
example deliberately does not attempt. The recommendation above does not
depend on it.
