Photometry and several spectra, with calibration uncertainty
================================================================

:doc:`sed_composition` is the one-model, two-instrument case: a spectrum and
a photometric catalogue. Real SED fits are usually richer than that — more
than one spectrograph, each with its own uncertain absolute flux
calibration — and that is the case ``spectrum_photometry.md`` §4 leaves to
an example and no example shows, until now. This page is the worked case:
the runnable example is ``examples/photometry_spectra``, a small package
(:download:`__init__.py <../../examples/photometry_spectra/__init__.py>`,
:download:`__main__.py <../../examples/photometry_spectra/__main__.py>`,
:download:`generators.py <../../examples/photometry_spectra/generators.py>`,
:download:`photometry_spectra.py <../../examples/photometry_spectra/photometry_spectra.py>`)
that fits one :class:`~ampere.backends.reference.ModifiedBlackBody` to
**three observations of one source** — a photometric catalogue and two
spectra of different resolving power — with each spectrum's calibration
factor either free or tied to the other, and the flexible likelihood as
the complement for what tying a calibration factor cannot explain (W6.2).

Run it yourself::

    python -m examples.photometry_spectra                  # untied, independent noise
    python -m examples.photometry_spectra --tie
    python -m examples.photometry_spectra --gp
    python -m examples.photometry_spectra --tie --gp

``tests/examples/test_photometry_spectra.py`` is this page's own coverage —
see that module's docstring for exactly what it does and does not check.

The physics, briefly
---------------------

A spectrograph measures a shape well and an absolute level poorly: its
wavelength-dependent throughput is calibrated against a standard, and that
calibration carries its own uncertainty, so the flux it reports is the
truth times an unknown multiplicative factor close to one — the "calibration
factor" this page is about. A photometric catalogue is not immune to this
either, but by convention it is treated as the calibration reference a
spectrum's factor is measured relative to, which is why the catalogue in
this composition carries none. Two spectrographs covering adjacent
wavelength ranges (here, the Spitzer IRS "short-low" and "long-low" modules'
shape, synthetic) each carry their *own* factor, and there is no reason the
two should agree — :mod:`examples.photometry_spectra.generators` gives them
different truths, 0.92 and 1.08, on purpose.

The five nouns, for three observations
-----------------------------------------

The negotiate-then-``compile_for`` dance, the label-collision lesson, and
why synthetic data needs the negotiated grid rather than the model's own are
all :doc:`sed_composition`'s ground to cover, not this page's — read that
page first if you have not. What is new here is a third observation on the
same channel, ``"sed"``: **one model**
(:class:`~ampere.backends.reference.ModifiedBlackBody`, the sibling
example's exact declaration — same priors, same ``channels="sed"``);
**three instruments**, all bound to channel ``"sed"`` with distinct labels
— ``"sl"``, a slit spectrograph at resolving power 100 over 5-14 micron;
``"ll"``, one at resolving power 60 over 14-38 micron; ``"catalogue"``, six
bundled mid/far-infrared filters (:data:`~examples.photometry_spectra.photometry_spectra.FILTERS`,
spanning about 22-161 micron) via
:meth:`~ampere.backends.reference.SyntheticPhotometry.from_library`; **three
datasets**, one per instrument, each with its own
:class:`~ampere.core.Likelihood`; **one fitting problem**, negotiating all
three instruments' requirements onto one union grid exactly as the sibling
page's two-instrument case does, just with a third source added to the
union.

Each spectrograph's chain is
:class:`~ampere.backends.reference.LSFConvolution`, then
:class:`~ampere.backends.reference.Resample`, then
:class:`~ampere.backends.reference.CalibrationScale` — the same three-step
shape the sibling page's ``"irs"`` instrument uses, once per spectrograph.

Calibration as a step with a prior
-------------------------------------

:class:`~ampere.backends.reference.CalibrationScale`'s own docstring gives
the natural prior: log-normal about one, since a calibration factor is
positive and multiplicative. Declared once and reused for both
spectrographs' own steps:

.. code-block:: python

    CALIBRATION_PRIOR = st.lognorm(0.05, scale=1.0)

    sl = Instrument(
        [LSFConvolution(resolving_power=100.0), Resample(SL_OBSERVED_WAVELENGTH),
         CalibrationScale(CALIBRATION_PRIOR)],
        channel="sed", label="sl",
    )
    ll = Instrument(
        [LSFConvolution(resolving_power=60.0), Resample(LL_OBSERVED_WAVELENGTH),
         CalibrationScale(CALIBRATION_PRIOR)],
        channel="sed", label="ll",
    )

Left alone, this is the **untied** case: ``sl.instrument.calibration_scale.scale``
and ``ll.instrument.calibration_scale.scale`` are two independent free
parameters, each with its own posterior.

**Tying** the two together — one shared factor rather than two — is a
composition-time decision, not a change to either step, and it is declared
on the :class:`~ampere.core.dataset.FittingProblem`, not on the collection:
its own docstring is explicit about why (``ampere/core/dataset.py`` line
2262) — *"Declared here because this is the only level at which every site
is visible; a tie names sites by their full merged path"*:

.. code-block:: python

    tie = Tie(
        "calibration",
        ("sl.instrument.calibration_scale.scale", "ll.instrument.calibration_scale.scale"),
        prior=CALIBRATION_PRIOR,
    )
    problem = FittingProblem(model, datasets, ties=(tie,), seed=seed)

``problem.parameters.free_names`` then carries one ``"calibration"`` name in
place of the two qualified ones — the tied problem has exactly one fewer
free parameter than the untied one (4 against 5, or 8 against 9 with the
GP arm below), which
``tests/examples/test_photometry_spectra.py::TestTheCompositionBuilds::test_tying_removes_exactly_one_free_parameter``
asserts directly.

But the truth gives the two spectrographs *different* factors, 0.92 and
1.08, and one shared parameter cannot recover two different numbers. The
measured cost, at 68 % (:func:`~examples.photometry_spectra.photometry_spectra.fit`'s
own defaults, seed :data:`~examples.photometry_spectra.generators.SEED`) is
the first two rows of the table below: **without** the flexible likelihood,
the tied factor does not even land at the two truths' mean (1.00) — it
settles at 1.21, pulled past both truths by the "ll" spectrum's unexplained
residual (below), and drags ``model.beta`` and ``model.scale`` off with it.
**With** the flexible likelihood, the tied factor lands at 0.96, close to
the mean the untied case's own two numbers bracket, and the physical
parameters recover.

The flexible likelihood as the complement
---------------------------------------------

The generator also injects a smooth Gaussian bump — 10 % amplitude, 4
micron width, centred at 25 micron — into the ``"ll"`` spectrum alone
(:data:`~examples.photometry_spectra.generators.BUMP_FRACTION`,
:data:`~examples.photometry_spectra.generators.BUMP_WIDTH`,
:data:`~examples.photometry_spectra.generators.BUMP_CENTRE`): the classic
"the dust is not one blackbody" feature, and something no three-parameter
greybody can express regardless of temperature, emissivity index or scale.
``build_problem(gp=True)`` keeps :class:`~ampere.core.IndependentNoise` on
the photometry (three or four broadband points carry no exploitable
correlation structure at this resolution) but gives each spectrum its own
:class:`~ampere.core.GaussianProcessNoise` — a Matérn-3/2 in wavelength with
log-uniform ("shrinkage-free") priors on the amplitude and length scale,
following :class:`~ampere.core.GaussianProcessNoise`'s own docstring
example; a prior with mass piled up near zero would pull a genuinely
present feature back towards "no GP needed", which is the opposite of what
this page needs to show.

When the data are a handful of photometric points with no shared structure
between them, the GP has nothing to learn from: its hyperparameters stay at
their prior and the fit is the independent one at extra cost. The flexible
likelihood earns its place on a spectrum, or on a catalogue dense enough to
carry correlated residuals — which is why the photometry here keeps
:class:`~ampere.core.IndependentNoise`.

Four combinations, four fits, the physical parameters and the calibration
factor(s) at 68 % central credible intervals (truth in brackets; the full
95 % coverage check and every free parameter are in the branch report and
``tests/examples/test_photometry_spectra.py``):

.. list-table::
   :header-rows: 1

   * - Arm
     - Temperature (K) [180]
     - :math:`\beta` [1.6]
     - Calibration, sl [0.92]
     - Calibration, ll [1.08]
   * - Untied, independent
     - 179.07 [178.70, 179.42]
     - 1.616 [1.605, 1.628]
     - 0.912 [0.893, 0.931]
     - 1.098 [1.077, 1.119]
   * - Untied, GP
     - 177.74 [176.50, 178.93]
     - 1.624 [1.594, 1.653]
     - 0.934 [0.902, 0.966]
     - 1.057 [1.012, 1.099]
   * - Tied, independent
     - 178.99 [178.64, 179.33]
     - 1.484 [1.473, 1.495]
     - 1.210 [1.184, 1.236] (shared)
     - 1.210 [1.184, 1.236] (shared)
   * - Tied, GP
     - 178.24 [177.00, 179.56]
     - 1.606 [1.576, 1.637]
     - 0.958 [0.922, 0.994] (shared)
     - 0.958 [0.922, 0.994] (shared)

Reading the table is the point. **Untied, independent** looks almost fine —
both calibration factors land close to their own truths — but the fit is
precise enough (226 points across three datasets) that the unmodelled bump
still shows up as a small, real bias: ``model.temperature``'s 68 % interval
does not reach 180 K, and neither does its 95 % one (a MISS in the branch
report, by a fraction of a kelvin). **Untied, GP** absorbs the bump instead
of biasing the physical parameters; both calibration factors and the
temperature move back within reach of their truths (95 % coverage: all
five parameters). **Tied, independent** is the sharpest illustration in the
table: forced to explain both spectrographs' calibration *and* the "ll"
residual with one number, the shared factor overshoots past both truths to
1.21, and ``model.beta`` (1.48 against a truth of 1.6) and ``model.scale``
are pulled along with it — every one of the four parameters misses at 95 %
in the branch report, not only at 68 %. **Tied, GP** recovers: with the
bump absorbed separately, the one shared factor is free to settle near the
two truths' honest compromise (0.96, against a mean of 1.00), and the
physical parameters recover alongside it.

.. _photometry-spectra-limits:

An upper limit is a likelihood statement, not a datum
-----------------------------------------------------

A non-detection is not a measurement of zero with a small error bar, and
fitting it as one pulls the model towards a flux nobody observed. HD 105's
catalogue (``HD105_SED.csv``) flags two far-infrared points as upper limits,
SPIRE 500 micron and LABOCA 870 micron, and the star-disc twin
(:doc:`migrating`) declares them as censored observations on its photometry:

.. code-block:: python

    names, _, limit_flux, limit_error = generators.read_limits()  # the two flagged rows
    names = [*votable_names, *names]                             # nineteen, then the limits
    flags = np.arange(len(names)) >= len(names) - 2             # True for the last two
    photometry = Likelihood(
        GaussianFamily(),
        IndependentNoise(),
        censoring=Censoring.upper_limits(flags),
    )

Each flagged point then contributes the log of the probability that the true
flux lies *below* its recorded value, :math:`\log \Phi(z)` with :math:`z =
(\text{recorded} - \text{predicted}) / \sigma`, in place of the Gaussian
density; the other nineteen are untouched. This is the Tobit term of the
likelihood contract (``docs/design/contracts/likelihoods.md``, section 9). Two
rules go with it. A masked sample contributes nothing whatever its limit
kind, so masking beats censoring. And a limit under a correlated-noise model
is refused when the problem is composed, not approximated, so declare limits
on an independent-noise dataset: which is why the twin puts its Gaussian
process on the RVS spectrum and not on the photometry. The recorded value of
a limit is whatever the catalogue holds; here that is a flux below three
times its error, flagged by the catalogue, and a reader who wants a stated
:math:`3\sigma` limit declares that number as the recorded value instead.

Running it
-------------

.. code-block:: text

    python -m examples.photometry_spectra                  # untied, independent noise
    python -m examples.photometry_spectra --tie             # one shared calibration factor
    python -m examples.photometry_spectra --gp              # the flexible likelihood on both spectra
    python -m examples.photometry_spectra --tie --gp

Each finishes in one to three minutes on the reference backend with
:class:`~ampere.inference.EmceeEngine` (71-173 s across the four arms in one
specific run — untied/independent 98 s, untied/GP 173 s, tied/independent
71 s, tied/GP 155 s; the GP arms are slower, both for the extra dimensions
and for the Cholesky solve :class:`~ampere.core.GaussianProcessNoise`'s
default :class:`~ampere.core.DenseGP` does every evaluation). Add
``--figures DIR`` to any of them to write ``corner.png``,
``posterior_predictive.png`` and (for
the ``--gp`` arms) ``gp_localisation.png`` into *DIR* — the last one is the
GP's own conditioned mean on the ``"ll"`` dataset, showing the bump it
localised. No figure is committed (AGENTS.md ground rule 7).

**Where the walkers start.** The example does not draw its initial
positions from the prior, which is :class:`~ampere.inference.EmceeEngine`'s
default; it starts every walker in a tight ball around the known synthetic
truth (``_initial_positions`` in ``photometry_spectra.py``). Writing this
page found out why that matters: with ``model.scale``'s prior spanning two
decades, about one walker in ten drawn from the joint prior lands on a
finitely-scored but astronomically improbable point and then never accepts
a proposal in thousands of steps, silently corrupting every flattened
statistic. Starting near a good point is ordinary MCMC practice; on a
synthetic problem the truth is that point, and on real data it is a
preliminary optimum — the warm start Phase 6's optimisers module is for.

See also
---------

* :doc:`sed_composition` — the two-instrument groundwork this page builds
  on: the five nouns, the negotiated grid, and the label-collision lesson.
* :doc:`m2_misspecification` — the flagship misspecification study the same
  flexible likelihood underlies, at three data sizes and on all three
  backends.
* :doc:`astrometry` — another worked example of
  :class:`~ampere.core.GaussianProcessNoise` against an injected systematic,
  on a very different kind of axis (time, not wavelength).
