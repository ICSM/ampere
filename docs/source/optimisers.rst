Point estimates and warm starts
===============================

Every sampling engine starts, by default, from draws of the joint prior that
merely score finitely (:meth:`~ampere.inference.Engine.initial_positions`).
That is honest and needs no tuning, and it is also where a hard problem goes
wrong first. ``examples/photometry_spectra`` measured it (W6.2): on that
composition about one emcee walker in ten, started from a technically
scoreable but astronomically improbable prior draw — ``model.scale``'s prior
alone spans two decades — never accepts a single proposal in thousands of
steps, because the stretch move cannot climb back from so far off, and every
flattened statistic of the run is corrupted by it. The tutorial worked around
it by starting near the known synthetic truth, which a real fit does not
have.

:func:`ampere.inference.optimise` is the general answer: a cheap point
estimate of the posterior mode, stored as an :class:`ampere.results.Optimum`
that every sampler can start from. :func:`ampere.inference.warm_start_gp` does
the same for a GP likelihood's hyperparameters, in milliseconds.
``inference.md`` §10b is the contract.

What is maximised
-----------------

Every route maximises the **constrained-space posterior density**,
``log p(θ) + log p(D | θ)`` at ``θ = constrain(u)``, over the packed
unconstrained vector ``u`` — so every iterate stays inside the prior's
support, and the answer does not depend on which bijection maps a bounded
parameter to the real line. (``log_prob_unconstrained``, the density NUTS
samples, adds the change-of-variables term, and its maximum moves with the
bijection.) The :class:`~ampere.results.Optimum` records both numbers, each
named for what it is.

The three routes
----------------

**scipy** runs on every backend, without gradients: multi-start
:func:`scipy.optimize.minimize` (Powell by default; ``minimiser="L-BFGS-B"``
selectable) from ``starts`` prior draws, the best kept, and the covariance
from a central-difference Hessian.

.. code-block:: pycon

    >>> import numpy as np, scipy.stats as st, astropy.units as u
    >>> from ampere.backends.reference import PowerLaw
    >>> from ampere.core import Dataset, FittingProblem, Spectrum
    >>> from ampere.inference import optimise
    >>> grid = np.geomspace(1.0, 10.0, 20)
    >>> rng = np.random.default_rng(7)
    >>> flux = 2.0 / grid + rng.normal(0.0, 0.2, grid.size)
    >>> observed = Spectrum(grid * u.um, flux * u.Jy, uncertainty=np.full(20, 0.2) * u.Jy)
    >>> model = PowerLaw(grid, norm=st.norm(2.0, 0.5), index=st.uniform(-3.0, 3.0),
    ...                  reference_wavelength=1.0)
    >>> problem = FittingProblem(model, [Dataset(observed)], seed=20260929)
    >>> optimum = optimise(problem, method="scipy", starts=4)
    >>> optimum.route, optimum.converged, len(optimum.starts)
    ('scipy', True, 4)
    >>> {name: round(value, 3) for name, value in optimum.constrained.items()}
    {'model.norm': 2.014, 'model.index': -1.103}

**map** is the gradient route on a torch or jax problem, through
:func:`ampere.core.realise`: ``torch.optim.LBFGS`` or
``jax.scipy.optimize.minimize(method="BFGS")``, with Adam (``optax.adam`` on
jax) as the fallback — a quasi-Newton step from a prior draw can leap into a
saturated bijection, so the fallback walks into the bulk with Adam and
finishes with the quasi-Newton method again. The covariance is an autodiff
Hessian. It finds the same mode as the scipy route to ``1e-3`` in every
unconstrained coordinate (a conformance row says so on every backend with a
realisation), and on the reference backend it refuses by name:

.. code-block:: pycon

    >>> optimise(problem, method="map")
    Traceback (most recent call last):
        ...
    ampere.inference.exceptions.EngineError: optimise(method='map') needs a gradient through ampere.core.realise, ...

On a native problem ``method="auto"`` (the default) takes this route, and
the scipy route everywhere else.

**vi** fits :class:`~ampere.inference.VIEngine`'s ``laplace`` guide and
returns its mean and covariance. The Laplace guide centres on the mode of
the density NUTS samples, so on a problem with bounded parameters it sits a
little apart from the other two routes: it is the alternative *start*, not a
third estimate of the same number.

**warm_start_gp** handles a :class:`~ampere.core.GaussianProcessNoise`
likelihood's ``amplitude``, ``length_scale`` and ``scale``. On a
reduced-rank solver (:class:`~ampere.core.HilbertSpaceGP`,
:class:`~ampere.core.EquispacedFourierGP`) the covariance is linear in its
features at a fixed length scale, so one eigendecomposition of the whitened
feature Gram matrix per length scale gives the marginal likelihood in closed
form as a function of the noise-to-amplitude ratio; the amplitude profiles
out, a bracketed root find gives the ratio, and a twelve-point log grid over
the length-scale prior's central 99 % gives the rest. A
:class:`~ampere.core.DenseGP` dataset is searched on a temporary
``HilbertSpaceGP(basis_size=32)``, which the result's ``message`` says.

.. code-block:: pycon

    >>> from ampere.core import GaussianFamily, GaussianProcessNoise, HilbertSpaceGP
    >>> from ampere.core import Likelihood, Matern32
    >>> from ampere.inference import warm_start_gp
    >>> gp = Likelihood(GaussianFamily(), GaussianProcessNoise(
    ...     Matern32(st.loguniform(0.05, 5.0), st.loguniform(0.3, 20.0)),
    ...     HilbertSpaceGP(basis_size=24), scale=st.loguniform(0.3, 3.0)))
    >>> gp_problem = FittingProblem(model, [Dataset(observed, likelihood=gp)], seed=20260929)
    >>> hyper = warm_start_gp(gp_problem)["default"]
    >>> hyper.route, hyper.free_names
    ('empirical_bayes', ('default.likelihood.amplitude', 'default.likelihood.length_scale', 'default.likelihood.scale'))

The bridge to the samplers
--------------------------

An :class:`~ampere.results.Optimum` is a start. ``initial_positions(count,
around=optimum)`` draws a ball at the mode — ``u* + 0.5 L z`` in the
unconstrained coordinates, ``L`` the covariance's Cholesky factor, so the
ball is tighter than the posterior and no walker starts in a tail — and
every sampling engine's ``run(initial=)`` takes the optimum directly: emcee
and zeus through that ball, NUTS and blackjax with every chain at the mode
plus a small jitter, VI at the mode exactly. The nested samplers refuse it by
name: they draw their live points from the prior transform and have no start.
:meth:`Optimum.combine <ampere.results.Optimum.combine>` joins a model's
optimum and its GP's hyperparameters into one start.

.. code-block:: pycon

    >>> from ampere.inference import EmceeEngine
    >>> start = optimise(gp_problem, method="scipy", starts=2)
    >>> run = EmceeEngine(gp_problem, walkers=12).run(50, initial=start)
    >>> run.attrs["ampere_start_route"], run.attrs["ampere_schema_version"]
    ('scipy', 10)

What an Optimum records, and saving it
--------------------------------------

The mode as the packed unconstrained vector and as constrained values by
qualified name; the density in both conventions; the covariance (inverse
Hessian, unconstrained coordinates) or a ``covariance_refusal`` naming why
the Hessian was not positive definite — never a silently regularised matrix;
whether it converged, the optimiser's message and the evaluation count; one
summary per start; and the problem's provenance at the time. A run seeded
from it records ``ampere_start`` (the route, the optimum's identity hash, the
density, the evaluations and convergence) and every run records
``ampere_start_route`` (``"prior"`` by default) — provenance schema 9.
``to_datatree()`` writes an ``optimum`` group and **no** ``posterior``, so it
goes through the ordinary netCDF route and back:

.. code-block:: pycon

    >>> from ampere.results import Optimum, from_netcdf, to_netcdf
    >>> path = to_netcdf(optimum.to_datatree(), "optimum.nc")
    >>> Optimum.from_datatree(from_netcdf(path)).identity == optimum.identity
    True
    >>> print(optimum.summary())
    Optimum by route 'scipy' on backend 'reference': converged after ... evaluations over 4 start(s)
    log p (constrained) = 6.42836; log p (unconstrained) = 6.06802
    parameter       constrained   unconstrained     sd (unc.)
    model.norm          2.01405         2.01405      0.125458
    model.index         -1.1029         0.54238      0.157001
    message: scipy.optimize.minimize('Powell'): Optimization terminated successfully.

The honest limits
-----------------

* **A mode is not a posterior.** The covariance is the local curvature at one
  point; it says nothing about skew, heavy tails or a second mode. Use it to
  start a sampler, not to replace one.
* **Multimodality.** Multi-start keeps the best of ``starts`` local optima
  and summarises the others; it does not find every mode, and a best-by-
  objective start inside a minor mode is possible. The per-start summaries
  are there to be read: starts that ended at very different objectives are a
  warning.
* **The harvest's warning stands** (``docs/design/harvest/optim_only/``): a
  local optimiser on a badly conditioned or ridge-shaped posterior converges
  somewhere, reports success and is wrong about where the mass is. The
  covariance refusal catches a flat direction; it cannot catch a banana.

Bayesian optimisation is deliberately not offered. Its surrogate assumes a
high-dimensional standard-normal latent and a locally linear objective, which
a posterior of ten to fifty curved, often multimodal parameters does not
supply; an expensive simulator is the SBI layer's case first.
