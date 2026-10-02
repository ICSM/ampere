Conditional priors
==================

Most priors in ampere are independent: each parameter carries its own
distribution and the joint prior is their product. Sometimes the prior on one
parameter *depends on another* — the members of a population are drawn around a
shared mean, the upper of two temperatures must lie above the lower — and this
page is what v2 offers for that, and what it does not. There are two routes,
and they suit different cases:

* **a hierarchical prior**, for a parameter whose distribution has
  hyperparameters that are themselves fitted (:class:`~ampere.core.HierarchicalPrior`,
  and :class:`~ampere.core.Plate` for a population of members);
* **a reparameterisation**, for a constraint such as an ordering, where the
  trick is to fit different parameters and compute the physical ones from them.

Every block below is run by the documentation's doctest runner
(``tests/core/test_spec_doctests.py``), so what you read is what the code does.
The contract behind both routes is ``docs/design/contracts/parameters.md`` §9;
:doc:`arbitrary_priors` is the companion page about what a prior *is*.

A parameter whose prior has a fitted hyperparameter
---------------------------------------------------

A :class:`~ampere.core.HierarchicalPrior` names a distribution family, exactly
as a frozen ``scipy.stats`` object does, but refers to its parameters *by the
name of another parameter* in the same set. That is what a numpyro model does
anyway (``dist.Normal(mu, sigma)`` with ``mu`` and ``sigma`` earlier sample
sites), and it is why the construct lowers to the torch and jax backends like
any other prior. Every argument is a reference: a quantity you want held
fixed is declared as a fixed :class:`~ampere.core.Parameter` and referenced.

.. code-block:: pycon

    >>> import numpy as np
    >>> import scipy.stats as st
    >>> from ampere.core import HierarchicalPrior, Parameter, ParameterSet
    >>> pset = ParameterSet([
    ...     Parameter("mu", st.norm(0.0, 5.0)),
    ...     Parameter("sigma", value=0.5, fixed=True),
    ...     Parameter("theta", HierarchicalPrior("norm", {"loc": "mu", "scale": "sigma"})),
    ... ])
    >>> pset.free_names
    ('mu', 'theta')
    >>> bool(np.isclose(
    ...     pset.lnprior(np.array([1.0, 1.2])),
    ...     st.norm(0.0, 5.0).logpdf(1.0) + st.norm(1.0, 0.5).logpdf(1.2),
    ... ))
    True

``lnprior`` is the joint density written as a chain: ``mu`` is evaluated under
its own prior and ``theta`` under a normal centred on the value of ``mu`` it was
given. A set orders its evaluations topologically, so a hyperparameter is always
resolved before anything that refers to it, and a reference that does not
resolve, or a cycle, is refused when the set is built rather than part-way
through a fit:

.. code-block:: pycon

    >>> ParameterSet([Parameter("theta", HierarchicalPrior("norm", {"loc": "mu"}))])
    Traceback (most recent call last):
        ...
    ampere.core.exceptions.ParameterError: parameter 'theta' has a hierarchical prior referencing 'mu', which is not in this set (it declares ['theta']). Hierarchical references must resolve within the set that will evaluate them.

The shape of the dependence is the only thing a hierarchical prior can say: the
*family* is fixed and its arguments are other parameters, one for one. It cannot
say "location ``a``, scale ``80 - a``", because an argument is a name, not an
arithmetic expression — which is the reason the ordering case below is done the
other way.

A population: N members sharing hyperparameters
-----------------------------------------------

The commonest conditional structure in astronomy is a population — *N* objects
whose parameter is each drawn from one distribution that is itself unknown.
:class:`~ampere.core.Plate` declares it in one place: the hyperparameters, and
the member parameters that refer to them.

.. code-block:: pycon

    >>> from ampere.core import Plate
    >>> objects = Plate(
    ...     "objects",
    ...     size=4,
    ...     hyperparameters=[
    ...         Parameter("mu", st.norm(0.0, 5.0)),
    ...         Parameter("sigma", st.halfnorm(0.0, 2.0)),
    ...     ],
    ...     members=[
    ...         Parameter("theta", HierarchicalPrior("norm", {"loc": "mu", "scale": "sigma"})),
    ...     ],
    ... )
    >>> population = ParameterSet([Parameter("background", st.norm(0.0, 1.0))], plates=[objects])
    >>> population.names
    ('background', 'objects.mu', 'objects.sigma', 'objects.theta')
    >>> population.free_size            # 1 + 1 + 1 + 4
    7

A plate is a constructor, not a container: it expands into ordinary parameters,
the hyperparameters qualified by the plate's name and each member a single
array-valued parameter of shape ``(size,)``, so nothing downstream needs to know
it was there. The members are drawn from the population distribution implied by
the hyperparameters currently in force, which is what a sampler needs:

.. code-block:: pycon

    >>> drawn = population.unpack(population.prior_transform(np.full(7, 0.6)))
    >>> mu, sigma = drawn["objects.mu"], drawn["objects.sigma"]
    >>> bool(np.allclose(drawn["objects.theta"], st.norm(mu, sigma).ppf(0.6)))
    True

For a whole fitting problem — one model per member, a dataset per member —
:class:`~ampere.core.Population` is the composition-time form of the same idea,
and :doc:`population` is the tutorial for it: a joint fit on the native path,
and reweighting of fits you have already run.

An ordering: reparameterise, do not constrain
---------------------------------------------

The other conditional prior people reach for is an *ordering*: two temperatures
with ``T0 < T1``, or two line centres that must not swap. ampere has no prior
that is a density on the ordered triangle, and it does not need one. The
construction in the NGC 6302 twin (``examples/ngc6302/ngc6302.py``) gives that
joint density exactly, with no rejection step. Fit instead a *triangular*
prior on the lower temperature and a *fraction* of the remaining interval for
the upper one,

.. math::

    T_0 \sim \mathrm{Triangular}(10, 80) \;\propto\; 80 - T_0, \qquad
    f \sim \mathrm{Uniform}(0, 1), \qquad
    T_1 = T_0 + f\,(80 - T_0),

and compute :math:`T_1` inside the model. The map :math:`(T_0, f) \to (T_0,
T_1)` has Jacobian :math:`\partial T_1/\partial f = 80 - T_0`, so the joint
density of the pair is

.. math::

    p(T_0, T_1) = \frac{p(T_0)\,p(f)}{|\partial T_1 / \partial f|}
                = \frac{k\,(80 - T_0)\cdot 1}{80 - T_0} = k,

constant on the ordered triangle :math:`\{10 \le T_0 < T_1 \le 80\}` — the flat
ordered prior, reproduced from two independent one-dimensional priors. The
cancellation can be checked numerically, at any points of the triangle:

.. code-block:: pycon

    >>> T0_prior = st.triang(c=0, loc=10.0, scale=70.0)     # density proportional to 80 - T0
    >>> f_prior = st.uniform(0.0, 1.0)
    >>> T0 = np.linspace(10.5, 79.5, 7)
    >>> f = np.linspace(0.05, 0.95, 7)
    >>> log_p = T0_prior.logpdf(T0) + f_prior.logpdf(f) - np.log(80.0 - T0)
    >>> bool(np.ptp(log_p) < 1e-12)
    True
    >>> bool(np.isclose(np.exp(log_p[0]), 2.0 / 70.0**2))
    True

The declaration is two ordinary parameters per pair, and the model reads them:

.. code-block:: pycon

    >>> pair = ParameterSet([
    ...     Parameter("T0", T0_prior),
    ...     Parameter("T_fraction", f_prior),
    ... ])
    >>> pair.free_names
    ('T0', 'T_fraction')
    >>> T0_draw, fraction = pair.prior_transform(np.array([0.5, 0.25]))
    >>> T1_draw = T0_draw + fraction * (80.0 - T0_draw)
    >>> bool(10.0 <= T0_draw < T1_draw <= 80.0)
    True

The same trick covers any constraint that can be written as "``y`` between
``x`` and a bound": a fraction of the remaining interval, with a prior on ``x``
chosen so that the Jacobian cancels. It buys the ordering for every engine at
once, because every engine sees only independent priors; the price is that
:math:`T_1` is a *derived* quantity of the model rather than a parameter, so it
does not appear in the posterior unless you compute it from the draws (the twin
does, in its ``report``).

What you cannot declare today
-----------------------------

An arbitrary joint density on two parameters — a prior that correlates them
because a previous analysis did, a banana-shaped constraint, a density known
only up to a normalising constant — is not a :class:`~ampere.core.Prior`. The
protocol is one-dimensional by design (a parameter's prior applies
independently to every element of an array-valued parameter, and a
multivariate prior is out of scope for the parameter contract, §5), and a
hierarchical prior's arguments are names, not expressions. There is no
``Derived`` node in the contract either, which is the missing piece for
declaring a prior over a function of other parameters; its design is the
subject of a memo (W6.11), not something this page can promise.

There are two workarounds, in the order to try them:

1. **Reparameterise** so that the joint is a product of one-dimensional
   priors, as above. This is the only route that keeps every engine working.
2. **Fit under a simpler prior and reweight.** Fit with independent priors on
   the two parameters, then multiply each posterior draw by the ratio of the
   density you wanted to the density you used. This is exact in the limit and
   inefficient when the two disagree strongly, because few draws then carry
   the weight; it needs no change to the fit. Check the effective sample size
   of the weights before trusting the result.
