"""M2 over many spectra: the per-dataset GP amplitude as a population (**W7.2**).

The milestone study (:mod:`examples.m2_misspecification`) fits one spectrum at
a time, so the flexible likelihood's GP amplitude is one number with one prior.
A survey is N spectra of one object class, observed through one instrument,
and the question it raises is how those N amplitudes should be declared. This
package is the measured answer for the three declarations a user can write:

**independent**
    Each spectrum's ``likelihood.amplitude`` keeps the study's own half-normal
    prior — N independent amplitudes, the status quo.
**tied**
    One :class:`~ampere.core.Tie` over every ``d*.likelihood.amplitude`` — one
    amplitude for all N spectra, the wrong model when the misspecification
    differs from object to object.
**population**
    The design memo's §4 (``docs/design/nuisance_populations_and_derived_memo.md``):
    a :meth:`~ampere.core.DatasetCollection.plate` over ``d*.likelihood``,
    non-centred — an internal ``z`` per spectrum, the amplitude derived as
    ``exp(mu + sigma * z)`` and routed by element to each spectrum's noise model
    — so the amplitudes share one prior whose location and spread are fitted.

The scenario is :data:`~examples.m2_populations.populations.N_SPECTRA` spectra
of the M2 toy, the ``strong_sharp`` deviation injected in three of them and one
clean member observed at five times the noise, so the population of amplitudes
has a genuine spread (most members near zero, a few large) and one member whose
own data barely constrain its amplitude.

The modules
-----------
:mod:`~examples.m2_populations.populations`
    The scenario, the three arms, the driver (:func:`~examples.m2_populations.populations.compare`)
    and the reductions ``tests/m2/test_populations.py`` asserts on.
:mod:`~examples.m2_populations.figures`
    The page figure ``docs/source/population.rst`` shows.

``python -m examples.m2_populations`` runs the comparison and prints its table.
Nothing it writes is committed. The study's own pieces — the generator, the
model, the likelihood, the engine dispatch, the budgets' types and the
:class:`~examples.m2_misspecification.study.Summary` reduction — are imported
from :mod:`examples.m2_misspecification` rather than copied, and no new
scenario joins the milestone's fixed set.
"""
