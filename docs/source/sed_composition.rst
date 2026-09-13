Combining a spectrum and photometry
=====================================

Real sources are rarely observed once. A dust-dominated infrared source
might have a spectrograph scan and a handful of catalogue magnitudes, both
telling you about the same physical emission — and ampere's answer is not a
special "SED fitting" mode, but the same five nouns :doc:`overview` already
introduced, composed twice over one model channel. This page is the worked
case: the runnable example is ``examples/sed_composition``, a small package
(:download:`generators.py <../../examples/sed_composition/generators.py>`,
:download:`sed_composition.py <../../examples/sed_composition/sed_composition.py>`,
:download:`__main__.py <../../examples/sed_composition/__main__.py>`) that
fits one :class:`~ampere.backends.reference.ModifiedBlackBody` to a spectrum *and* a
photometric catalogue at once, on the reference backend with
:class:`~ampere.inference.EmceeEngine` or on ``torch``/``jax`` with
:class:`~ampere.inference.NUTSEngine` — the same script, one flag. It is
memo §7.1's worked example (``docs/design/phase4_placement_memo.md``),
expanded into a runnable example and this page (W4.11).

Run it yourself::

    python -m examples.sed_composition                 # reference, emcee, ~1-2 minutes
    python -m examples.sed_composition --backend torch  # NUTS

``tests/examples/test_sed_composition.py`` is this page's own coverage — see
that module's docstring for exactly what it does and does not check.

The physics, briefly
---------------------

The injected source (``examples.sed_composition.generators.TRUTH``) is a
180 K, optically thin dust greybody with emissivity index 1.6 — cool enough
that it has essentially nothing to say in the near infrared and most of its
flux beyond 20 micron, which is realistic for the instrument combination
below (a debris disk or an embedded protostar, not a stellar photosphere).
``scale`` is dimensionless — it is what turns :math:`B_\nu(T)`, in Jy/sr,
into an observed flux density, i.e. the source's solid angle in effect — so
its natural size is of order :math:`10^{-15}`, not a flux; the flux the
instruments actually see comes out at a few tenths of a Jy to a few Jy, a
sensible size for a catalogued mid/far-infrared source. The catalogue's five
filters, below, are mid/far-infrared only for exactly this reason — a first
draft of this example included 2MASS J/Ks and found the source at
:math:`10^{-19}` Jy there: harmless for emcee, but the extreme dynamic range
across the dataset made the NUTS variant's warmup pathologically slow, for a
point that constrains nothing. A filter with no real signal is not free
information; it is a numerical liability.

The five nouns, for this case
-------------------------------

**One model.** :class:`~ampere.backends.reference.ModifiedBlackBody` (``ampere.backends.reference``,
or its ``torch``/``jax`` twin), on channel ``"sed"``, given its own fallback
grid at construction. That grid is a starting point only — the next
paragraph is why it never actually gets used once an instrument is bound.

**Two instruments, one channel.** ``"irs"`` — a slit spectrograph:
:class:`~ampere.backends.reference.LSFConvolution` at constant resolving power 100, then
:class:`~ampere.backends.reference.Resample` onto a 120-point observed grid from 5 to 35
micron, then a :class:`~ampere.backends.reference.CalibrationScale` nuisance parameter for
the instrument's own flux calibration. ``"catalogue"`` —
:meth:`~ampere.backends.reference.SyntheticPhotometry.from_library`, five
bundled mid/far-infrared filters (WISE W3/W4, IRAS 60, Spitzer MIPS 70,
Herschel PACS 100) tabulated on a 500-point grid. Both bind channel ``"sed"``.

**Two datasets**, one per instrument, each with its own
:class:`~ampere.core.Likelihood` — :class:`~ampere.core.GaussianFamily` over
:class:`~ampere.core.IndependentNoise` (the ordinary chi-square; there is no
GP here, so no flexible-likelihood story on this page).

**One fitting problem.** :class:`~ampere.core.dataset.FittingProblem` binds
the model and the :class:`~ampere.core.DatasetCollection`, and does
everything below in its constructor: negotiates the channel's requirements,
compiles the model onto the union grid, merges the parameters, validates.

Why synthetic data needs the negotiated grid
----------------------------------------------

There is no reader in v2 (:doc:`migrating`), so an example builds its
"observed" data by hand: evaluate the model at the truth, push the result
through each instrument, add noise. The tempting way to write that —
``model(**truth)`` on the model's own declared grid, then straight into the
photometry step — fails, and it is worth seeing the failure once, because it
is the reason ``examples.sed_composition.generators.synthetic_data``
does something less obvious:

.. code-block:: pycon

    >>> raw_result = model(**generators.TRUTH)           # the model's own grid
    >>> camera(raw_result)
    Traceback (most recent call last):
        ...
    ampere.core.exceptions.SchemaError: axis 'spectral_axis' has no coordinate
    within COORDINATE_RTOL (1e-12) of 1.0092715146305713, 1.0186289902446874,
    1.0280732238308652, 1.0376050197669118, 1.0472251898884348 (and 494 more).
    The axis spans [1, 200] over 2000 point(s). A step should only ask for
    coordinates it published as a points= requirement, so this is a
    negotiation defect rather than a usage error.

The photometry step tabulated its filter responses on 500 particular
wavelengths and asks for exactly those back
(:meth:`~ampere.core.Axis.locate`); the model's own 2000-point
fallback grid does not contain them, so the lookup fails, correctly — "a
negotiation defect rather than a usage error" is the message being explicit
that this is not a bounds check with a sensible fallback, because a silent
nearest-neighbour substitution here is exactly the bug class
(``spectrum_photometry.md`` Gap 1) that a *published* requirement exists to
prevent.

The fix is memo §7.1 finding 1: do what
:class:`~ampere.core.dataset.FittingProblem` itself does before evaluating
anything —

.. code-block:: python

    requirements = negotiate([spectrograph, camera])
    compiled = model.compile_for(requirements)
    truth = compiled(**generators.TRUTH)
    spectrum_truth = spectrograph(truth, {"calibration_scale.scale": generators.CALIBRATION_TRUTH})
    photometry_truth = camera(truth)

— :func:`~ampere.core.negotiate` collects what both instruments ask of
channel ``"sed"``, :meth:`~ampere.core.transform.Model.compile_for` adopts
the union grid as the model's template (once; every later evaluation refills
values on the same axis object), and only then does an instrument have
something it is safe to read from. Passing the same, now-compiled ``model``
on to :class:`~ampere.core.dataset.FittingProblem` costs nothing extra: its
constructor calls ``compile_for`` again on the same requirements, which is
idempotent.

Requirements, printed and explained
--------------------------------------

Once the problem is built, ``problem.requirements`` is the negotiation
record — what each source asked of the channel, and what the union of those
demands turned into. Running the example prints exactly this::

    negotiated channels: ['sed']
    sources asking of channel 'sed': ('irs', 'catalogue')
      axis spectral_axis: <AxisRequirement 'spectral_axis' micron [4.89383,35.7432]@R>=706.446 [5,35]@step<=0.12605 500 point(s)>
    free parameters: ('model.temperature', 'model.beta', 'model.scale', 'irs.instrument.calibration_scale.scale')

One requirement, three clauses, each contributed by a different step:

- ``R>=706.446`` over ``[4.89, 35.74]`` micron is :class:`~ampere.backends.reference.LSFConvolution`'s:
  it learned its output range from :class:`~ampere.backends.reference.Resample` *after* it
  in the chain (:meth:`~ampere.core.transform.Transformation.configure_from`,
  gap I-3) — padded by several kernel widths beyond ``[5, 35]`` so the
  convolution is never evaluated against an edge — and asked to be sampled
  several times finer than its own resolving power (100) so the kernel is not
  itself under-resolved.
- ``[5,35]@step<=0.12605`` is :class:`~ampere.backends.reference.Resample`'s own density
  requirement: cover its target range, sampled finer than half its coarsest
  output bin.
- ``500 point(s)`` is :class:`~ampere.backends.reference.SyntheticPhotometry`'s exact
  tabulation — the same 500 wavelengths the previous section's failure was
  about.

Negotiation takes the union: the coverage from every clause, the sampling
density at its strictest wherever the clauses overlap, and the photometry
step's exact points folded in regardless. The model never sees three
requests — it is compiled once, onto the single grid that satisfies all
three, and the four free parameters are the physics (temperature, beta,
scale) plus the one nuisance parameter the "irs" instrument's calibration
scale contributes.

The label collision, and the fix
-----------------------------------

An :class:`~ampere.core.Instrument`'s label defaults to its channel name, so
two instruments bound to one channel collide **unless you name them**. Build
the two instruments the way a first attempt might, without ``label=``:

.. code-block:: pycon

    >>> spectrograph = Instrument([LSFConvolution(resolving_power=100.0), Resample(grid)], channel="sed")
    >>> camera = Instrument([SyntheticPhotometry.from_library(["WISE_RSR_W3"], tabulation)], channel="sed")
    >>> spectrograph.label, camera.label
    ('sed', 'sed')
    >>> DatasetCollection([Dataset(observed_spectrum, spectrograph), Dataset(observed_photometry, camera)])
    Traceback (most recent call last):
        ...
    ampere.core.exceptions.DatasetError: two datasets are labelled 'sed'. Labels
    are the component names of the joint parameter space, so they must be
    unique. Both labels came from an Instrument's own label, which itself
    defaults to the channel name — so two instruments bound to one channel
    collide unless you name them. Pass label='...' to the Dataset (or to the
    Instrument), or build the collection from a mapping whose keys are the
    labels.

The message names the fix, and ``examples.sed_composition.sed_composition.build_instruments``
takes it: ``label="irs"`` and ``label="catalogue"`` at construction. This is
not a corner case — ``core/dataset.py``'s own comment on the check names
"the photometry case" as its motivating example — it is what happens the
first time a second instrument reads a channel a first one already does, and
:doc:`the composition paragraph in the overview page </overview>` links here
for exactly this reason.

Fitting it: reference backend, emcee
----------------------------------------

.. code-block:: python

    from examples.sed_composition.sed_composition import build_problem, fit, report

    problem = build_problem("reference")
    run = fit(problem, backend="reference")   # EmceeEngine, the module's own defaults
    print(report(run))

Sizing the budget is measurement, not guesswork: twenty calls to
``problem.log_prob`` at the truth (this is W4.11's own acceptance-sizing
step) took in the high teens of milliseconds each under ampere's default,
multi-threaded BLAS — the LSF convolution's ``(n, n)`` operator is now built
once and cached across draws (W4.0), but the matrix-vector multiply against
it still runs every evaluation, on a roughly 2 100-point union grid. Small
matrices like that do not repay multi-threaded BLAS's own thread-pool
overhead, so ``python -m examples.sed_composition`` restricts every BLAS
thread pool to one thread before numpy is imported (see
``examples/sed_composition/__main__.py``'s docstring) — measured, that
roughly halves the per-evaluation time. With that restriction, 16 walkers x
450 steps (180 discarded as burn-in — 7 200 evaluations) finished in one to
under two minutes on the reference backend (65-100 s across a few runs on a
machine shared with other work — see the branch report for one specific
run's numbers), and recovered all four parameters inside their central 95 %
interval:

.. code-block:: text

    emcee on reference: 16 chain(s) x 270 draw(s)
      posterior (truth in brackets, 95 % coverage flagged):
        irs.instrument.calibration_scale.scale        +0.997 +- 0.034   95%[+0.937, +1.080]  (truth +1.02)  ok
        model.beta                                    +1.526 +- 0.431   95%[+0.508, +2.202]  (truth +1.6)  ok
        model.scale                                   +2.24e-15 +- 2.6e-15   95%[+4.30e-16, +8.14e-15]  (truth +1e-15)  ok
        model.temperature                             +177.9 +- 26.6   95%[+87.2, +216.4]  (truth +180)  ok

(Figures are rounded here; the branch report for W4.11 has the exact values
from one specific run — every random draw is seeded, so re-running the
script reproduces them exactly.) ``model.temperature``, ``model.beta`` and
``model.scale`` are the classic modified-blackbody degeneracy — wide, and
correlated with each other — which is exactly why the interval, not the
point estimate, is what this page's acceptance criterion is about.

The NUTS variant: one flag
------------------------------

.. code-block:: python

    problem = build_problem("torch")            # or "jax"
    run = fit(problem, backend="torch")          # NUTSEngine

Nothing about the *composition* changes: the same ``ModifiedBlackBody``,
the same two instruments on the same channel, the same two labelled
datasets. What changes is which module
supplied every backend-declared piece —
``examples.sed_composition.sed_composition.backend_module`` reaches
into ``ampere.backends.torch`` instead of ``ampere.backends.reference`` for
the model and the instrument steps, and
``examples.sed_composition.sed_composition.noise_module`` reaches into
the same module for :class:`~ampere.backends.torch.IndependentNoise` rather
than ``ampere.core.IndependentNoise`` — composing the reference
backend's ``IndependentNoise`` into a torch problem is refused as a backend
disagreement (:doc:`overview`'s one-backend rule), because it would leave a
numpy island in an otherwise differentiable chain. :class:`~ampere.core.Instrument`,
:class:`~ampere.core.Dataset` and :class:`~ampere.core.dataset.FittingProblem`
are backend-neutral, so none of that code changes at all — only
``examples.sed_composition.sed_composition.build_model`` and
``examples.sed_composition.sed_composition.build_instruments``, which
this example already routes through one backend-name flag.

See also
---------

* :doc:`overview` — the five nouns, and the one-backend rule this page's NUTS
  variant relies on.
* :doc:`m2_misspecification` — the flexible likelihood this page's
  composition does not use; the GP over the residuals is a property of a
  dataset's own :class:`~ampere.core.Likelihood`, so it composes with
  everything on this page unchanged.
* :doc:`sbi` — for a model with no likelihood to write down at all, rather
  than the two ordinary ones here.
