Fitting a population
======================

:doc:`advanced`'s "Populations: fitting many objects together" section
introduces :class:`~ampere.core.Population` — this page is not that
introduction again, it links to it. What follows is a tutorial: one small
population, declared once, fitted two different ways — jointly, on the
native path, and by reweighting archived single-object fits — because
Phase 5 shipped both routes (``hierarchical_population.md`` §11 Q5 left the
choice open between them) and a real project's data usually decides which
one is affordable.

The population: fifty objects, one shared index
--------------------------------------------------

The running example is the one :doc:`advanced` sketches and
``tests/inference/test_population_nuts.py`` fits for real: fifty power-law
spectra, each observed at three wavelengths, whose indices are believed to
be draws from one Gaussian population rather than fifty independent
numbers.

.. code-block:: python

    import numpy as np
    import scipy.stats as st
    import astropy.units as u

    from ampere.core import (
        Dataset, FittingProblem, GaussianFamily, HierarchicalPrior,
        Instrument, Likelihood, Parameter, Population, Spectrum,
    )
    from ampere.backends.torch import IndependentNoise, PowerLaw

    GRID = np.array([1.0, 3.0, 9.0])
    NOISE = 0.02
    rng = np.random.default_rng(20260919)
    indices = rng.normal(-1.30, 0.35, 50)

    models, datasets = {}, []
    for i, index in enumerate(indices):
        label = f"obj{i}"
        truth = 2.0 * GRID**index
        observed = Spectrum(
            GRID * u.micron,
            (truth + rng.normal(0.0, NOISE, GRID.size)) * u.Jy,
            uncertainty=np.full(GRID.size, NOISE) * u.Jy,
        )
        models[label] = PowerLaw(
            GRID, norm=Parameter("norm", value=2.0, fixed=True), index=st.norm(-1.3, 1.0)
        )
        datasets.append(Dataset(
            observed,
            Instrument([], channel="default", input_kind=Spectrum, label=f"scope{i}"),
            Likelihood(GaussianFamily(), IndependentNoise()),
            model=label, label=f"d{i}",
        ))

    population = Population(
        "objects",
        members=[Parameter("index", HierarchicalPrior("norm", {"loc": "mu", "scale": "sigma"}))],
        hyperpriors=[Parameter("mu", st.norm(-1.0, 1.0)), Parameter("sigma", st.halfnorm(0.0, 1.0))],
        over=[f"obj{i}" for i in range(50)],
    )
    problem = FittingProblem(models, datasets, populations=[population], seed=20260919)

Nothing about ``models`` or ``datasets`` knows a population exists — each
``PowerLaw`` still declares its own ``index`` prior, and ``Population``
replaces it at merge time. ``problem.free_size`` is 52: fifty member draws
plus ``objects.mu`` and ``objects.sigma``, not fifty-two separately-costed
parameter objects — that distinction is what makes the fit tractable at
this scale, and :doc:`advanced` has the measured cost of getting it wrong
(46 ms per ``lnprior`` evaluation at a hundred hand-written members, 373 ms
at a thousand, against about 1.3 ms for the same structure as one plate).

**Two layouts, one declaration.** ``layout="plate"`` is the default used
above: the fifty draws are **one** array-valued parameter, which is what
lowers to a real ``pyro``/``numpyro`` plate on the torch and jax backends.
``layout="flat"`` declares exactly the same density as fifty individually
tied scalar parameters — the pattern ``parameters.md`` §9 describes for a
handful of objects — and is refused above ``MAX_FLAT_MEMBERS`` (128) by
default, because the flat layout's cost is linear in the number of
parameter *objects* while
the plate's is not. A caller who has read that cost and wants the flat
layout anyway past the limit can turn the refusal into a loud warning with
``ampere.core.settings.override(flat_population_cap="warn")`` — a
process-wide setting, not a ``Population`` argument, so it does not change
the limit itself.

When each member *is* a dataset — a spaxel of an IFU cube, one catalogue
entry — :meth:`DatasetCollection.plate <ampere.core.DatasetCollection.plate>`
is the more direct route: it builds the same ``Population`` declaration
from a list of datasets in plate order, routing draws to the model label
each dataset already names, so the loop above collapses into one call.

The same population, non-centred
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The declaration above is *centred*: each ``index`` is sampled directly from
``Normal(mu, sigma)``. The same density can be sampled in independent
coordinates instead — ``z ~ Normal(0, 1)`` per member and
``index = mu + sigma * z`` — which removes the funnel a gradient sampler
meets when the data say little about each member. ``z`` is declared by no
member model, so it is an *internal* member: sampled on the population's own
component and routed nowhere. ``index`` is a :class:`~ampere.core.Derived`
member, computed from the others wherever values are formed and routed to each
model exactly as the centred draw was. A small version, on plain parameter
sets:

.. code-block:: pycon

    >>> import numpy as np
    >>> import scipy.stats as st
    >>> from ampere.core import Derived, Parameter, ParameterSet, Population
    >>> non_centred = Population(
    ...     "objects",
    ...     members=[Parameter("z", st.norm(0.0, 1.0)),
    ...              Parameter("index", Derived("mu + sigma * z"))],
    ...     hyperpriors=[Parameter("mu", st.norm(-1.0, 1.0)),
    ...                  Parameter("sigma", st.halfnorm(0.0, 1.0))],
    ...     over=["obj0", "obj1", "obj2"],
    ... )
    >>> objects = {
    ...     f"obj{i}": ParameterSet([Parameter("index", st.norm(-1.3, 1.0))]) for i in range(3)
    ... }
    >>> merged = ParameterSet.merge(objects, populations=[non_centred]).merged
    >>> merged.free_names
    ('objects.mu', 'objects.sigma', 'objects.z')
    >>> values = merged.complete({"objects.mu": -1.3, "objects.sigma": 0.35,
    ...                           "objects.z": np.array([-1.0, 0.0, 2.0])})
    >>> values["objects.index"]
    array([-1.65, -1.3 , -0.6 ])

The sampler sees ``mu``, ``sigma`` and the three ``z``; ``index`` is no
sampler dimension, yet every model receives it, and the ``posterior`` of a
run carries it as a variable of its own beside the sampled ones (named in the
run's ``ampere_derived`` attribute). For the fifty-object problem above the
change is the ``members=`` list alone. Which form samples better depends on
the data: the centred form suits members the data pin down individually, the
non-centred one members the data barely constrain (``inference.md`` §9).

The joint fit, on the native path
-------------------------------------

The population's own hyperparameters and every member's index are one
posterior, fitted jointly by NUTS through :func:`~ampere.core.realise` —
Phase 5's first evidence that the joint route in
``hierarchical_population.md``'s design horizon (a) actually works at a
scale past what the flat layout tolerates:

.. code-block:: python

    from ampere.inference import NUTSEngine

    run = NUTSEngine(problem).run(draws=300, warmup=300, chains=1)
    mu = run.posterior["objects.mu"].values.reshape(-1)
    sigma = run.posterior["objects.sigma"].values.reshape(-1)

This needs a differentiable backend (torch or jax) to lower the plate and
run NUTS — it is not reachable on the reference backend, which is why the
code above imports ``ampere.backends.torch`` for the model rather than
``ampere.core``'s own ``PowerLaw``. ``tests/inference/test_population_nuts.py``
is this exact fit, run on both installed backends, and it is the acceptance
evidence for the number worth stating precisely: at fifty well-measured
members, the truth used to generate the sample sits inside the posterior's
central 95 % interval on both ``objects.mu`` and ``objects.sigma``, and the
posterior concentrates on the **sample's** mean and scatter rather than
drifting toward the hyperprior's own location.

Reweighting archived fits, when refitting jointly is not affordable
------------------------------------------------------------------------

Sometimes the fifty (or two hundred) single-object fits already exist —
run independently, on whatever engine and schedule suited each one, long
before anyone asked a population question of them. Refitting all of them
jointly to get a population posterior throws that work away.
:func:`~ampere.results.fit_population` (**W5.13**) does not: it recovers
the population hyperparameters from the archive alone, by self-normalised
importance reweighting (Hogg, Myers & Bovy 2010), with no joint fit at all.

The example below is deliberately small — twelve objects, each a single
free parameter read straight back from one datum, fitted with
:class:`~ampere.inference.EmceeEngine` on the reference backend — so that
it runs in seconds rather than minutes, but the call is the one a two
-hundred-object archive makes too:

.. code-block:: python

    import numpy as np
    import scipy.stats as st
    import astropy.units as u

    from ampere.core import Dataset, FittingProblem, Model, ModelResult, Parameter, Spectrum
    from ampere.inference import EmceeEngine
    from ampere.results import DataTreeRunColumns, GaussianPopulationModel, fit_population

    class ConstantModel(Model):
        """The smallest model this contract can fit: theta, read straight back."""

        def __init__(self):
            self.register_buffer("grid", np.array([1.0]), unit=u.micron)
            self.register_parameter(Parameter("theta", st.norm(0.0, 10.0)))

        def evaluate(self, **values):
            ctx = self.context(values)
            return ModelResult(
                Spectrum(ctx["grid"] * u.micron, np.full_like(ctx["grid"], ctx["theta"]) * u.Jy)
            )

    sigma = 0.2
    rng = np.random.default_rng(28)
    truths = rng.normal(1.0, 0.3, 12)
    data = truths + rng.normal(0.0, sigma, 12)

    def fit_object(datum, seed):
        observed = Spectrum(
            np.array([1.0]) * u.micron, np.array([datum]) * u.Jy,
            uncertainty=np.array([sigma]) * u.Jy,
        )
        problem = FittingProblem(ConstantModel(), [Dataset(observed)], seed=seed)
        return EmceeEngine(problem, walkers=8).run(steps=300, burn_in=100)

    runs = [fit_object(float(datum), seed=28 + i) for i, datum in enumerate(data)]
    columns = [DataTreeRunColumns(run) for run in runs]

    population_model = GaussianPopulationModel(
        Parameter("mu", st.norm(0.0, 5.0)), Parameter("tau", st.halfnorm(scale=2.0)),
    )
    result = fit_population(
        columns, "model.theta", population_model, st.norm(0.0, 10.0), steps=1500, burn_in=500, seed=1,
    )

Run on the branch, the twelve-object archive above (truth ``mu = 1.0``,
``tau = 0.3``) gives a 95 % interval of ``(0.77, 1.26)`` on ``mu`` and
``(0.15, 0.61)`` on ``tau``, both covering the truth, with a per-object
effective sample size no worse than 989 — well clear of
:func:`~ampere.results.fit_population`'s default refusal floor, which
exists so a population fit never silently trusts an archived run whose own
posterior is not well sampled. The whole cell — twelve ``EmceeEngine`` fits
plus the reweighting — runs in under thirty seconds.

The interim prior matters here, and it is not a free choice: it must be the
**marginal** prior each archived object was actually fitted under
(``st.norm(0.0, 10.0)`` above), because the reweighting divides each
object's stored posterior density by it. Since **W5.22**, supplying it is
optional — every run's provenance stores its own free parameters' priors
(``ampere_free_priors``, schema 8), so ``fit_population`` reads the stored
one back when none is given, and refuses by name rather than guessing
whenever a caller-supplied prior disagrees with the archive's own record.

What a larger archive looks like: ``population_full``
----------------------------------------------------------

The twelve-object archive above is sized for a tutorial, not for a claim.
``tests/results/test_population.py``'s own two-hundred-object archive —
registered under the ``population_full`` marker, and not run by default
because building it costs several minutes — is the row that actually
validates the method at a realistic size: it recovers ``mu`` and ``sigma``
to intervals consistent with a fifty-object reduction of the same draw, and
is where to look for the calibrated behaviour this page's toy only
gestures at.

Simulation-based calibration on a population
-------------------------------------------------

A population declared on a :class:`~ampere.core.FittingProblem`
survives an SBC replica exactly as an ordinary problem's parameters do:
since **W5.30 (c)**, :func:`~ampere.results.replace_observations`
passes ``populations=problem.populations`` into the rebuilt replica, so a
population-level SBC run is not silently independent draws for each
member — a gap the item closed after :doc:`astrometry`'s own joint-noise
calibration study found the analogous one for ``joint=``.

See also
------------

* :doc:`advanced` — the populations section, for the concept, the two
  layouts' cost trade-off and the flat-population-cap setting in full.
* ``docs/design/modalities/hierarchical_population.md`` §9, §11 — the
  design horizons (a) (joint) and (b) (reweighting) this page's two routes
  answer, and gap H-1's disposition.
* ``docs/design/contracts/parameters.md`` §8–§9 — the frozen ``Population``
  and ``Plate`` contract.
* ``docs/design/contracts/inference.md`` §9 — ``DatasetCollection.plate``
  and the native lowering.
