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

.. _photometry-spectra-jwst:

Reading a JWST spectrum: the front door
-----------------------------------------

Everything above starts from synthetic arrays. A real spectrum starts from a
file, and for JWST that file is the pipeline's one-dimensional product: an
``x1d`` (one ``EXTRACT1D`` table per extracted source or slit) or a ``c1d``
(one ``COMBINE1D`` table per combined source).
:func:`~ampere.spectroscopy.read_jwst` reads one table into a
:class:`~ampere.core.Spectrum` in one line, and :mod:`ampere.spectroscopy` is
the front door it lives behind — the kinds, the families, the reference
instrument steps and the reader in one import, in the shape
:mod:`ampere.interferometry` set for interferometry (:doc:`interferometry`
§10). ``import ampere`` does not import it.

.. code-block:: python

    from ampere.spectroscopy import read_jwst

    spectrum = read_jwst("tests/data/jw01536-o027_t008_miri_p750l_x1d.fits")

The file is the test suite's own: a MIRI low-resolution-spectrometer
fixed-slit spectrum of the flux standard BD+60 1753 from calibration
programme 1536, 196 kB, public domain (``tests/data/README.md`` has its
provenance). What the reader did with it is in ``spectrum.meta``:

* **where it came from** — ``instrument`` ``"MIRI"``, ``detector``,
  ``filter`` ``"P750L"``, ``grating`` (``None``: MIRI has none),
  ``exp_type`` ``"MIR_LRS-FIXEDSLIT"``, ``program``, ``observation``,
  ``target`` ``"BD+60 1753"``, ``date_obs``, ``pipeline_version`` (the
  file's ``CAL_VER``, ``"3.0.0"``), and the table read: ``product`` ``"x1d"``,
  ``extension`` (its ``EXTVER``), ``slit`` and ``srctype`` ``"POINT"``;
* **what it changed** — the table's 388 rows run from long to short
  wavelength, so they were reversed (``reordered = "reversed"``); the last row
  repeats the wavelength of the one before it, and an axis cannot repeat a
  coordinate, so the padding row of the two was dropped (``n_duplicate = 1``),
  leaving 387 samples;
* **what it masked** — ``mask`` is the pipeline's ``DO_NOT_USE`` data-quality
  bit (``DQ & 1``), together with any non-finite value or uncertainty and any
  uncertainty that is not positive; ``n_do_not_use``, ``n_nonfinite``,
  ``n_nonpositive_error`` and ``n_masked`` count them (one padding row here,
  at the long-wavelength end). The other ``DQ`` bits mask nothing and are kept
  raw in ``spectrum.extra_coords["dq"]``, beside ``npixels``.

The axis is in micron, as the reference steps take it; the flux is in the
unit the file declares (Jy here, read from the column's ``TUNIT``, never
assumed). ``read_jwst(path, surface_brightness=True)`` reads ``SURF_BRIGHT``
and ``SB_ERROR`` in MJy/sr instead — the quantity to fit for an extended
source. A file holding several spectra (a NIRSpec fixed-slit exposure writes
one table per slit) is refused, with every table listed, until ``source=``
(a slit name such as ``"S200A1"``, or a source id) or ``extension=`` (the
table's ``EXTVER``) says which; the reader never picks the first for you.

**The line-spread function is yours to give.** The pipeline does not write a
resolving power into its products, so the reader carries none. Supply the
mode's resolution to :class:`~ampere.backends.reference.LSFConvolution` from
the instrument's documentation: today the step takes one constant
``resolving_power`` (or one ``fwhm``), so take your resolution curve's value
at the band's centre — about 100 near 7.5 µm for the LRS prism. A
wavelength-dependent :math:`R(\lambda)` is the step's next extension; the LRS
curve rises from about 40 at 5 µm to about 160 at 10 µm (the MIRI
documentation's figures), and
a fit that cares about line shapes across the band should know that it is
using one value.

**The fit.** The ``Spectrum`` composes into this page's pattern unchanged — a
model on a grid covering the band, the instrument chain built on the
observed axis, and the flexible likelihood:

.. code-block:: python

    import astropy.units as u
    import numpy as np
    import scipy.stats as st

    from ampere.backends.reference import PowerLaw
    from ampere.core import Dataset, FittingProblem, GaussianProcessNoise, Instrument
    from ampere.core import Likelihood, Matern32
    from ampere.inference import EmceeEngine
    from ampere.spectroscopy import GaussianFamily, LSFConvolution, Resample, read_jwst

    spectrum = read_jwst("tests/data/jw01536-o027_t008_miri_p750l_x1d.fits")
    model = PowerLaw(np.geomspace(3.0, 15.0, 400), norm=st.loguniform(1e-3, 1e-1),
                     index=st.uniform(-4.0, 4.0), reference_wavelength=7.5)
    instrument = Instrument([LSFConvolution(resolving_power=100.0),
                             Resample(spectrum.spectral_axis.values)])
    kernel = Matern32(st.loguniform(1e-5, 1e-2), st.loguniform(0.1, 10.0),
                      amplitude_unit=u.Jy, length_scale_unit=u.um)
    likelihood = Likelihood(GaussianFamily(), GaussianProcessNoise(kernel))
    problem = FittingProblem(model, {"miri": Dataset(spectrum, instrument, likelihood)})
    run = EmceeEngine(problem, walkers=16).run(1000, burn_in=500)

The block is not executed when the documentation builds; it runs as
``tests/spectroscopy/test_jwst.py::TestTheFit``, at a token budget of eight
steps and with ``initial="prior"`` in place of the engine's default optimiser
warm start, since the row asserts that the reader's ``Spectrum`` composes and
fits, not what the fit finds. The ``model`` and ``instrument`` lines are
this page's example with one spectrograph and no photometry. A power law is a deliberately crude model of a star's
spectrum — the Rayleigh–Jeans tail of an A star falls about as
:math:`\lambda^{-2}`, but the calibrated spectrum is not a pure power law
across a factor of four in wavelength — and the Matérn-3/2 process is what
absorbs the difference while the index is fitted: the page's theme, on a
real file.

At the block's budget the fit took about twelve minutes on one laptop core
(the optimiser's warm start, then 16 000 evaluations of a 387-sample dense
Gaussian process) and returned a median index of −1.60 with the GP's
amplitude at about 1.2 mJy and its length scale at about 0.25 µm; emcee's
diagnostics flagged the index (R-hat 1.11), so a fit meant to be quoted
runs longer. That the index is shallower than −2 with a short-length-scale
process beside it is the flexible likelihood doing its job, not a stellar
measurement: the residual structure is the instrument's and the model's, and
this page does not try to tell them apart.

See also
---------

* :doc:`sed_composition` — the two-instrument groundwork this page builds
  on: the five nouns, the negotiated grid, and the label-collision lesson.
* :doc:`m2_misspecification` — the flagship misspecification study the same
  flexible likelihood underlies, at three data sizes and on all three
  backends.
* :mod:`ampere.spectroscopy` — the front door's API page, and
  :func:`~ampere.spectroscopy.read_jwst`'s full statement of what it reads,
  masks and refuses.
* :doc:`astrometry` — another worked example of
  :class:`~ampere.core.GaussianProcessNoise` against an injected systematic,
  on a very different kind of axis (time, not wavelength).
