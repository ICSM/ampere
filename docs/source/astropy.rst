Using an astropy model
========================

`astropy.modeling` is where a great deal of astronomy's analytic models
already live, and a great many people arrive at ampere with one in hand.
This page is the reference for :mod:`ampere.core.astropy_compat` (landed at
W4.6) and its curated native translations (W4.7): the casual route that
wraps any astropy model as an ampere ``Model``, the capability consequence
of doing so, and the opt-in translation to a differentiable native
equivalent for the six common cases. ``tests/core/test_astropy_compat.py``,
``test_astropy_engines.py``, ``test_astropy_backend_hook.py``,
``test_astropy_translations.py`` and ``test_astropy_translations_backends.py``
are this page's own coverage; ``docs/design/contracts/astropy_compat.md``
is the binding contract this page summarises.

1. The casual route: ``from_astropy``
-----------------------------------------

One function, one class, both exported from :mod:`ampere.core`:

.. code-block:: python

    from_astropy(model, *, kind=None, priors=None, channel="default",
                 grid=None, output_unit=None, equivalencies=()) -> AdaptedAstropyModel

It wraps **any** ``astropy.modeling`` model — compound models included — as
an ordinary ampere :class:`~ampere.core.Model`, by translating a declaration
that already exists on the astropy side rather than asking for it again:

.. code-block:: pycon

    >>> import astropy.units as u, numpy as np, scipy.stats as st
    >>> from astropy.modeling.models import PowerLaw1D
    >>> from ampere.core import from_astropy
    >>> grid = np.array([1.0, 2.0, 4.0]) * u.micron
    >>> model = from_astropy(
    ...     PowerLaw1D(amplitude=2.0 * u.Jy, x_0=1.0 * u.micron, alpha=1.0),
    ...     grid=grid,
    ...     priors={"amplitude": st.loguniform(0.1, 10.0), "x_0": 1.0, "alpha": st.norm(1.0, 0.5)},
    ...     output_unit=u.Jy,
    ... )
    >>> model.kind.__name__, model.parameters.free_names
    ('Spectrum', ('amplitude', 'alpha'))

Each astropy ``Parameter`` becomes exactly one of three things, tried in
this order: named in ``priors=`` wins outright (a frozen ``scipy.stats``
distribution to fit, a number to hold fixed, or a ready-made
:class:`~ampere.core.Parameter`); a ``tied`` callable becomes an
``AstropyTie`` — a derived quantity, not an ampere :class:`~ampere.core.Tie`,
because astropy's ``tied`` is an arbitrary Python function of the whole
model and belongs in none of ``Parameter``'s three states; ``fixed=True``
becomes a frozen parameter at its value; finite ``bounds`` become a uniform
prior over them. **A free, unbounded parameter is refused rather than given
a default** — ``BlackBody.temperature`` ships with ``bounds=(0, None)``,
which is not a uniform prior over anything, and an improper default here
would be a fit silently becoming a different fit.

astropy's own parameter names are kept exactly — a compound model's
``temperature_0``/``temperature_1`` are the posterior's names too — which is
what keeps a tie callable working, since it reads the astropy model by
attribute:

.. code-block:: pycon

    >>> from astropy.modeling.models import Gaussian1D, Const1D
    >>> compound = Gaussian1D(3.0, 6.0, 1.5) + Const1D(0.3)
    >>> compound.mean_0.tied = lambda m: float(m.amplitude_0.value) * 2.0
    >>> for name in ("amplitude_0", "stddev_0", "amplitude_1"):
    ...     getattr(compound, name).bounds = (0.01, 12.0)
    >>> from ampere.core import Spectrum
    >>> wrapped = from_astropy(compound, grid=grid, kind=Spectrum, output_unit=u.Jy)
    >>> wrapped.parameters.free_names          # mean_0 is derived, not sampled
    ('amplitude_0', 'stddev_0', 'amplitude_1')
    >>> sorted(wrapped.ties)
    ['mean_0']

The kind is the **caller's declaration**, inferred only where the grid's own
units settle it unambiguously (a length/frequency/energy axis is a
``Spectrum``, a time axis a ``TimeSeries``, two axes an ``Image``), and
never from the model's class name — a wrong kind guessed from a class name
would be a fit that runs and is silently wrong.

2. The solid-angle rule: no steradian is ever invented
------------------------------------------------------------

astropy's own ``BlackBody`` emits a **surface brightness** — everything per
steradian — and a surface brightness reaches a flux density only through a
solid angle. ``from_astropy`` will not supply one:

.. code-block:: pycon

    >>> from astropy.modeling.models import BlackBody
    >>> from ampere.core.exceptions import CompositionError
    >>> try:
    ...     from_astropy(
    ...         BlackBody(temperature=3000.0 * u.K, scale=1.0 * u.Jy / u.sr),
    ...         grid=grid, output_unit=u.Jy,
    ...         priors={"temperature": st.uniform(100.0, 9900.0), "scale": st.loguniform(0.1, 10.0)},
    ...     )
    ... except CompositionError as exc:
    ...     print(exc)
    from_astropy() cannot express this model's output, which astropy returns in Jy / sr,
    in the requested output_unit=Jy. No equivalency in force makes the conversion, and
    ampere will not invent one...

The remedy states the convention rather than hiding it:
``equivalencies=[astropy.units.dimensionless_angles()]`` declares a solid
angle of **exactly one steradian**, which is the same convention
:class:`ampere.backends.reference.BlackBody`'s dimensionless ``scale``
carries — the factor that absorbs the solid angle and the distance dilution
together. The two are the same physics written twice, and
``TestItAgreesWithTheReferenceModels`` compares their ``log_prob`` at the
conformance suite's exact tolerance rather than approximately. Both
conversions (input units, output units) happen **once, at configuration
time**, never per evaluation — the units trap ``DEVELOPMENT_PLAN.md`` §7
names.

3. The capability consequence
-----------------------------------

``AdaptedAstropyModel`` declares its four capability flags honestly:

.. code-block:: text

    DIFFERENTIABLE = False    BATCHABLE = False    DEVICE = "cpu"    BACKEND = "reference"

All four declared rather than inherited, as every piece of the reference
path declares them. The consequences follow directly from
:doc:`overview`'s capability ladder:

- **Reachable**: :class:`~ampere.inference.EmceeEngine`,
  :class:`~ampere.inference.DynestyEngine`, :class:`~ampere.inference.ZeusEngine`
  and :class:`~ampere.inference.SBIEngine` — SBI is the engine a wrapped
  external model exists for, and it is the one this adapter most obviously
  unlocks.
- **Not reachable, ever, on the black-box route**:
  :class:`~ampere.inference.NUTSEngine` and :class:`~ampere.inference.VIEngine`.
  There is no gradient through an arbitrary Python callable, so
  :func:`~ampere.core.realise` refuses the problem by name. That is the
  contract, not a gap to be closed later — closing it is what §4 below is
  for, and it closes it by using a **different model**.

``BATCHABLE = False`` deserves its own precision: astropy models *are*
vectorised over their input grid (one evaluation covers a whole spectrum),
but nothing in ``astropy.modeling`` takes a stack of *parameter* vectors in
one call, and ``BATCHABLE`` is a claim about theta. ``simulate_many``
therefore parallelises an adapted model by **process**, not by
vectorisation, which is what a black box wants.

4. The opt-in native route
--------------------------------

**Never silent, in either direction** — the 2026-09-01 ruling this whole
module implements. ``ampere.core.from_astropy`` never substitutes anything:
it always evaluates the caller's actual astropy model. A native equivalent
is reachable only through a **backend-scoped** hook —
``ampere.backends.torch.from_astropy``, ``ampere.backends.jax.from_astropy``
— which decomposes a compound model into its leaves against a curated table
and **raises for any leaf it has no row for**, rather than falling back to
the black box. Both directions of silence are forbidden and for symmetric
reasons: falling back silently hands somebody who asked for a
differentiable model one that is not, discovered several composition steps
later when :class:`~ampere.inference.NUTSEngine` refuses; substituting
silently hands somebody a model they did not write, and a curated
``BlackBody`` is not guaranteed numerically identical to astropy's.

W4.7 filled the table with six curated classes, and their compound sums,
products, differences and ratios (astropy's ``+ - * /``):

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - astropy class
     - Native translation
   * - ``BlackBody``
     - Each backend's own tested ``planck_jy``, with astropy's raw-unit
       correction folded in once.
   * - ``PowerLaw1D``
     - Direct: ``amplitude * (x / x_0) ** -alpha``.
   * - ``BrokenPowerLaw1D``
     - A ``where``-selected pair of power laws about the break wavelength.
   * - ``Polynomial1D``
     - Horner's method over the coefficients astropy already carries.
   * - ``Gaussian1D``
     - Direct: an amplitude, a mean and a standard deviation.
   * - ``Const1D``
     - Direct: one amplitude, broadcast over the grid.

The physics is written **once**, backend-neutral, in
:mod:`ampere.core.astropy_translations` — a small ``TranslationOps``
namespace (``exp``, ``where``, ``blackbody``) each backend implements in its
own array type — so "one table serves both backends" rather than two
tables that could silently drift apart. Parameters, priors, frozen-ness,
the channel and the kind all come from the **same**
:func:`~ampere.core.astropy_compat.translate_astropy_parameters` the
black-box route uses, so a native and a black-box fit of the same astropy
model cannot disagree about what a bound or a fixed value means:

.. code-block:: pycon

    >>> from ampere.backends.torch import from_astropy as native_from_astropy
    >>> compound = Gaussian1D(3.0, 6.0, 1.5) + Const1D(0.3)
    >>> native = native_from_astropy(compound, grid=grid, output_unit=u.Jy,
    ...     priors={"amplitude_0": st.loguniform(0.1, 10.0), "mean_0": 6.0,
    ...             "stddev_0": st.loguniform(0.5, 5.0), "amplitude_1": st.loguniform(0.01, 1.0)})
    >>> native.DIFFERENTIABLE
    True

NUTS runs directly on a ``Gaussian1D + Const1D`` translated this way, on
both torch and jax.

5. The two refusals
--------------------------

The native route refuses two things by name rather than attempting a
compromise:

**A tied parameter.** A compound model carrying an astropy ``tied=``
parameter is refused on both native backends: a tie is an arbitrary Python
callable evaluated on plain floats, so there is no gradient through it and
no way to evaluate it at all while jax is tracing. The black-box route is
unaffected — it still applies a tie exactly as astropy defines it, at the
cost of the capability consequences of §3.

**``|`` and ``&``.** astropy's other two composition operators have no
elementwise meaning as one channel's flux — ``|`` is functional composition
and ``&`` stacks independent axes — so a compound model built with either
is refused by name on the native route, the same refusal an untranslated
leaf gets:

.. code-block:: pycon

    >>> from ampere.core import astropy_components, translation_refusal
    >>> [type(part).__name__ for part in astropy_components(compound)]
    ['Gaussian1D', 'Const1D']
    >>> print(translation_refusal("torch", compound, {}))
    ampere.backends.torch.from_astropy() has no native translation for Const1D, Gaussian1D.
    Translation is opt-in and never silent (DEVELOPMENT_PLAN.md §2, ruled 2026-09-01), so
    this refuses rather than falling back to the black-box adapter...

The black-box route consumes both operators, and every model outside the
curated table, without complaint — the refusal is specifically about the
promise of a gradient, not about whether ampere can evaluate the model at
all.

See also
------------

* :doc:`overview` — the capability ladder this page's §3 is an instance of.
* :doc:`sed_composition` — a worked composition in the same shape as one
  built with ``from_astropy``.
* ``docs/design/contracts/astropy_compat.md`` — the frozen contract this
  page summarises, with the full parameter-translation table and every
  limitation.
