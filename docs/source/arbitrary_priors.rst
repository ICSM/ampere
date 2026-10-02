Arbitrary priors
================

The canonical way to declare a prior is a frozen ``scipy.stats`` distribution:
``st.norm(0.0, 1.0)``, ``st.loguniform(1e-3, 1e3)``. It is the only form ampere
can describe neutrally, and therefore the only form that can be lowered to the
torch and jax backends or recorded in a run's provenance. But a prior is
really a *protocol*, and anything that satisfies it will evaluate on the
reference path: the numpy engines (:class:`~ampere.inference.EmceeEngine`,
:class:`~ampere.inference.ZeusEngine`, :class:`~ampere.inference.DynestyEngine`
and the other nested samplers) and a simulation budget for
:class:`~ampere.inference.SBIEngine`. This page says what the protocol is, writes
three priors that are not frozen scipy distributions — your own class, one of
scipy's newer distribution objects, and an empirical prior built from the draws
of a previous run — and says what each costs. The contract is
``docs/design/contracts/parameters.md`` §4 and §6; the measurements behind the
section on scipy's new objects are ``docs/design/performance_memo.md`` §8.
:doc:`conditional_priors` is the companion page for priors that depend on other
parameters.

Every block below is run by the documentation's doctest runner
(``tests/core/test_spec_doctests.py``).

The protocol
------------

A prior is any object with three things:

* ``logpdf(x)`` — the log density (a discrete family exposes ``logpmf`` instead,
  and :func:`~ampere.core.parameter.log_density` accepts either);
* ``ppf(q)`` — the inverse CDF, which nested samplers use to map the unit cube
  onto the prior (``prior_transform``);
* ``support()`` — a *method* returning ``(lower, upper)``, which ampere reads to
  choose the map to unconstrained space for gradient-based engines (a finite
  interval gets a logit, a half-line a log, the whole line the identity).

and, if you want ``ParameterSet.sample`` to work (initialisation of an ensemble,
the simulation budget of an SBI fit), an ``rvs(size=None, random_state=None)``
method as scipy's own distributions have. The structural type is
:class:`~ampere.core.parameter.Prior`; a frozen scipy distribution satisfies all
four, which is why none of this is visible when you use one.

A frozen scipy distribution is the case to prefer whenever it exists, because it
alone can be *described*: :func:`~ampere.core.describe_prior` reads its family
name and numbers, and that description is what lowering and provenance use.
Anything else is a **duck-typed** prior, which evaluates but cannot be lowered or
serialised, and says so by name the moment you ask:

.. code-block:: pycon

    >>> import numpy as np
    >>> import scipy.stats as st
    >>> from ampere.core import describe_prior
    >>> describe_prior(object())
    Traceback (most recent call last):
        ...
    ampere.core.exceptions.ParameterError: prior <object object at ...> is not a frozen scipy.stats distribution, so it cannot be described neutrally (and therefore cannot be lowered to torch/numpyro or recorded in a run's provenance). Declare priors as e.g. scipy.stats.norm(0, 1); a custom prior object may still be evaluated on the reference path, but it must not be serialised or lowered.

That is the deliberate middle path: a user with a genuinely custom prior is never
blocked, and a prior that cannot travel never travels silently as far as a
backend before failing there.

Your own class
--------------

A truncated power law — a Salpeter-like initial mass function, say — has no
frozen scipy form, but its CDF inverts in closed form, which makes it a ten-line
class:

.. code-block:: pycon

    >>> class TruncatedPowerLaw:
    ...     """p(x) proportional to x**-index on [lower, upper]."""
    ...     def __init__(self, index, lower, upper):
    ...         self.index, self.lower, self.upper = index, lower, upper
    ...         self._k = 1.0 - index
    ...         self._norm = (upper**self._k - lower**self._k) / self._k
    ...     def logpdf(self, x):
    ...         x = np.asarray(x, dtype=float)
    ...         inside = (x >= self.lower) & (x <= self.upper)
    ...         with np.errstate(divide="ignore"):
    ...             return np.where(inside, -self.index * np.log(x) - np.log(self._norm), -np.inf)
    ...     def ppf(self, q):
    ...         q = np.asarray(q, dtype=float)
    ...         return (self.lower**self._k + q * self._k * self._norm) ** (1.0 / self._k)
    ...     def support(self):
    ...         return (self.lower, self.upper)
    ...     def rvs(self, size=None, random_state=None):
    ...         return self.ppf(np.random.default_rng(random_state).uniform(size=size))
    >>> from ampere.core import Parameter, ParameterSet
    >>> from ampere.core.parameter import Prior
    >>> mass = Parameter("mass", TruncatedPowerLaw(2.35, 0.1, 100.0))
    >>> isinstance(mass.prior, Prior)
    True
    >>> mass_set = ParameterSet([mass])
    >>> bool(np.isclose(mass_set.lnprior(np.array([1.0])), -np.log(mass.prior._norm)))
    True
    >>> bool(np.isclose(mass_set.prior_transform([0.0])[0], 0.1))
    True
    >>> bool(0.1 <= mass_set.sample(np.random.default_rng(0))["mass"] <= 100.0)
    True

Outside the support ``logpdf`` returns ``-inf``, which is what every engine
reads as "rejected"; it must never raise. ``support()`` is what lets the
gradient-based engines treat the parameter, and here it chooses a logit
between the bounds:

.. code-block:: pycon

    >>> from ampere.core import default_bijection_for
    >>> default_bijection_for(mass.prior)
    Logit(lower=0.1, upper=100.0)

What the class cannot do is travel. ``to_spec`` — and with it the torch and jax
lowering and a run's provenance record — needs a description, and a custom class
has none:

.. code-block:: pycon

    >>> mass_set.to_spec()
    Traceback (most recent call last):
        ...
    ampere.core.exceptions.ParameterError: prior <...TruncatedPowerLaw object at ...> is not a frozen scipy.stats distribution, so it cannot be described neutrally (and therefore cannot be lowered to torch/numpyro or recorded in a run's provenance). Declare priors as e.g. scipy.stats.norm(0, 1); a custom prior object may still be evaluated on the reference path, but it must not be serialised or lowered.

In practice: a duck-typed prior is for the numpy engines and for nested
sampling. If you need NUTS or variational inference on a torch or jax backend,
the prior must be a frozen scipy distribution, or a reparameterisation of one
(see the ordering example in :doc:`conditional_priors`).

scipy's new distribution objects
--------------------------------

scipy 1.15 added a new distribution infrastructure — ``scipy.stats.Normal``,
``Uniform``, ``Logistic`` and ``make_distribution`` for the rest — which is much
faster to evaluate (the memo measures 7.2 times on the six M2 priors, and
bit-identical results). It is natural to ask whether these objects are priors.
They are not, quite: they satisfy **half** the protocol.

.. code-block:: pycon

    >>> prior = st.Normal(mu=0.0, sigma=1.0)
    >>> isinstance(prior, Prior)
    False
    >>> pset = ParameterSet([Parameter("x", prior)])
    >>> bool(np.isclose(pset.lnprior(np.array([0.5])), st.norm(0.0, 1.0).logpdf(0.5)))
    True
    >>> pset.prior_transform([0.5])
    Traceback (most recent call last):
        ...
    AttributeError: 'Normal' object has no attribute 'ppf'...
    >>> pset.sample(np.random.default_rng(0))
    Traceback (most recent call last):
        ...
    AttributeError: 'Normal' object has no attribute 'rvs'...

``logpdf`` and ``support()`` are there, so constrained-space evaluation and the
map to unconstrained space work; but the new objects call the inverse CDF
``icdf`` and sampling ``sample``, so everything that draws from the prior or
needs ``ppf`` — nested sampling, ``ParameterSet.sample``, an SBI budget —
fails by name. Objects made by ``make_distribution`` have a second problem:
they are all instances of one anonymous ``CustomDistribution`` class, with no
family name for :func:`~ampere.core.describe_prior` to read, so they cannot be
lowered or recorded even with an adapter. The recommendation in the memo is not
to adopt them as priors; if you want one anyway, a four-method adapter makes it
a duck-typed prior, with every consequence of the previous section:

.. code-block:: pycon

    >>> class FromNewInfrastructure:
    ...     def __init__(self, dist):
    ...         self._dist = dist
    ...     def logpdf(self, x):
    ...         return self._dist.logpdf(x)
    ...     def ppf(self, q):
    ...         return self._dist.icdf(q)
    ...     def support(self):
    ...         return self._dist.support()
    ...     def rvs(self, size=None, random_state=None):
    ...         rng = np.random.default_rng(random_state)
    ...         return self._dist.sample(() if size is None else size, rng=rng)
    >>> adapted = ParameterSet([Parameter("x", FromNewInfrastructure(prior))])
    >>> adapted.prior_transform([0.5])
    array([0.])
    >>> bool(np.isfinite(adapted.sample(np.random.default_rng(0))["x"]))
    True

An empirical prior from a previous run's draws
----------------------------------------------

The most useful arbitrary prior is the posterior of an earlier analysis, used as
the prior of the next: a temperature measured from one band, carried into a fit
of another. Summarising the draws as a normal throws away what made the
posterior worth reusing, so use a kernel density estimate. It needs the three
protocol methods and a decision about the support, because a KDE has
unbounded support but a prior must say where its mass is.

The class below fits a ``scipy.stats.gaussian_kde`` to the draws, tabulates its
CDF on a grid, and inverts that by interpolation to get ``ppf``. Its
``support()`` is the range of the draws padded by three kernel bandwidths — the
region outside which the KDE has negligible mass — and ``logpdf`` is ``-inf``
beyond it, so the density is renormalised over exactly the interval the
protocol advertises:

.. code-block:: pycon

    >>> class EmpiricalPrior:
    ...     """A Gaussian-KDE prior over one parameter, from posterior draws."""
    ...     def __init__(self, draws, n_grid=2048, pad=3.0):
    ...         draws = np.asarray(draws, dtype=float)
    ...         self._kde = st.gaussian_kde(draws)
    ...         bandwidth = float(np.sqrt(self._kde.covariance[0, 0]))
    ...         self._lower = float(draws.min() - pad * bandwidth)
    ...         self._upper = float(draws.max() + pad * bandwidth)
    ...         self._grid = np.linspace(self._lower, self._upper, n_grid)
    ...         density = self._kde(self._grid)
    ...         self._norm = float(np.trapezoid(density, self._grid))
    ...         cdf = np.cumsum(density)
    ...         self._cdf = (cdf - cdf[0]) / (cdf[-1] - cdf[0])
    ...     def logpdf(self, x):
    ...         x = np.asarray(x, dtype=float)
    ...         inside = (x >= self._lower) & (x <= self._upper)
    ...         density = self._kde(np.atleast_1d(x)).reshape(x.shape) / self._norm
    ...         with np.errstate(divide="ignore"):
    ...             return np.where(inside, np.log(density), -np.inf)
    ...     def ppf(self, q):
    ...         return np.interp(np.asarray(q, dtype=float), self._cdf, self._grid)
    ...     def support(self):
    ...         return (self._lower, self._upper)
    ...     def rvs(self, size=None, random_state=None):
    ...         return self.ppf(np.random.default_rng(random_state).uniform(size=size))

Used on 4 000 draws standing in for an earlier posterior on a temperature, it
reproduces the density where the draws were and is zero where they were not:

.. code-block:: pycon

    >>> earlier = np.random.default_rng(7).normal(180.0, 12.0, 4000)    # a previous run's draws
    >>> empirical = EmpiricalPrior(earlier)
    >>> lower, upper = empirical.support()
    >>> bool(lower < earlier.min() and earlier.max() < upper)
    True
    >>> bool(abs(empirical.logpdf(180.0) - st.norm(180.0, 12.0).logpdf(180.0)) < 0.1)
    True
    >>> empirical.logpdf(np.array([0.0, 400.0]))
    array([-inf, -inf])
    >>> bool(abs(empirical.ppf(0.5) - 180.0) < 1.5)
    True
    >>> temperature = ParameterSet([Parameter("temperature", empirical)])
    >>> bool(abs(temperature.prior_transform([0.5])[0] - 180.0) < 1.5)
    True
    >>> default_bijection_for(empirical).lower == lower
    True

A few cautions apply to any prior built this way. The draws are only as good as
the sampler that produced them, so check their effective sample size first; a
KDE of fifty correlated draws is a bumpy prior. A posterior is a legitimate
prior only for data that are *independent* of the new ones — reusing the same
data twice counts them twice. And an empirical prior is a duck-typed one: it
runs on the numpy engines and the nested samplers, and it cannot be lowered, so
a gradient-based fit on torch or jax needs a frozen scipy distribution fitted to
the draws instead (a log-normal or a skew-normal often does well enough to
carry the information that matters).

A checklist for your own prior
------------------------------

* ``logpdf`` returns ``-inf`` outside the support and never raises.
* ``ppf`` is monotonic and maps ``[0, 1]`` onto the support; nested sampling
  depends on it.
* ``support()`` is a method, and the bounds it returns are the true bounds of
  the mass — a logit is built from them.
* ``rvs(size=, random_state=)`` is there if you use ``ParameterSet.sample``
  or an SBI budget.
* Declare it in the parameter's own unit; ampere does not convert priors.
* If it must be lowered or recorded, it is a frozen scipy distribution, or it
  is a reparameterisation of one.
