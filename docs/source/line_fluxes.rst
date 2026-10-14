Adding a modality: emission-line fluxes
========================================

The modality template applied to the measurement a spectroscopist most often
has (W7.9, issue #70): not a spectrum but a **table of integrated line
fluxes** — one flux and its uncertainty per transition, each transition
identified by a label and a rest wavelength — fitted with a model that
produces a spectrum, through an instrument that integrates that spectrum over
each line's window. A planetary nebula observed with JWST's mid-infrared
instrument is usually reduced to exactly such a table. Section order follows
:doc:`interferometry` ("the template") and :doc:`astrometry`, section for
section, so a gap shows as a gap.

The runnable example is ``examples/line_fluxes``
(:download:`generators.py <../../examples/line_fluxes/generators.py>`,
:download:`model.py <../../examples/line_fluxes/model.py>`,
:download:`study.py <../../examples/line_fluxes/study.py>`): a synthetic,
JWST-like list of ten mid-infrared nebular lines drawn from a two-temperature
excitation ladder, fitted with a one-temperature ladder with and without the
flexible likelihood. Run it yourself::

    python -m examples.line_fluxes                 # reference, emcee
    python -m examples.line_fluxes --quick         # the per-PR budget
    python -m examples.line_fluxes --backend jax   # NUTS on the native path

``tests/conformance/test_line_fluxes.py``,
``tests/backends/test_native_line_fluxes.py``,
``tests/core/test_results_schema.py::TestLineFluxes``,
``tests/core/test_dataset.py::TestLineFluxAlignmentByTransition`` and
``tests/examples/test_line_fluxes.py`` are this page's own coverage.

1. The kind: an identity, a rest wavelength, an energy flux
------------------------------------------------------------

:class:`~ampere.core.LineFluxes` is :class:`~ampere.core.PhotometricPoints`'
template with three things changed:

- **an identity.** A line flux *is* its transition, so the transition labels
  must be unique and are carried as ``extra_coords["transitions"]``. The kind
  declares ``IDENTITY = "transitions"`` (photometry declares ``"filters"``),
  which is what tells composition to align an observed table with the
  chain's prediction **by label**, in any order (§4);
- **a rest wavelength** as the one axis — the transition's laboratory
  wavelength, not where the line was observed. Order is arbitrary, as for
  photometry; the axis is what the lines are plotted against and what the
  flexible likelihood's kernel measures distances along (§5);
- **an energy-flux unit**: W m⁻², not a flux density. A published table in
  erg s⁻¹ cm⁻² is converted once with ``.to_unit`` — the likelihood's unit
  check is exact, as it is for every kind.

Building one by hand from a real table is one call per table — a label, a
rest wavelength, a flux and its uncertainty per row:

.. code-block:: python

    import astropy.units as u
    from ampere.spectroscopy import LineFluxes

    erg = u.erg / u.s / u.cm**2
    table = LineFluxes(
        ["[Ne II] 12.81", "[Ne III] 15.56", "[S IV] 10.51"],
        [12.8136, 15.5551, 10.5105] * u.um,
        [3.1e-13, 1.2e-13, 0.8e-13] * erg,
        uncertainty=[0.1e-13, 0.05e-13, 0.04e-13] * erg,
    ).to_unit(u.W / u.m**2)

2. The step: what it takes from the observed table, and what it publishes
---------------------------------------------------------------------------

``LineIntegration`` (``Spectrum -> LineFluxes``, in each backend's
``instrument`` module beside synthetic photometry) integrates the model's
flux density over a window around each line:

.. math::

    F_k = \int_{\mathrm{window}_k} \left(F_\nu - C_\nu\right)\, d\nu
    \quad [\mathrm{W\,m^{-2}}],

from a model in **Jy** on a **micron** grid, with ``dν = c dλ / λ²``
(``∫F_λ dλ`` is the same number and is not offered twice). Two declarations
shape it:

- **the window**: one half-width per line (or one for all), in km/s — the
  default, the spectroscopist's unit, so the window scales with the line's
  wavelength — or in µm;
- **the continuum**: ``"none"`` integrates the model as it is (a model of the
  lines alone, as here); ``"linear"`` subtracts, per draw, the straight line
  through the model's own values at the window's two edges. A line table was
  continuum-subtracted by whoever measured it, so a model that carries a
  continuum must be continuum-subtracted the same way before the two are
  compared.

Both are linear in the model's samples, so **the whole step is one fixed
``(n_lines, n_samples)`` influence matrix** built once — synthetic
photometry's shape exactly. Each window is tabulated on
``samples_per_window`` (default 64) evenly spaced wavelengths, edge to edge;
the union of those grids is the step's published ``points=`` requirement, and
the columns are placed on whatever grid negotiation hands back through
:meth:`~ampere.core.Axis.locate`, never positionally, so a second instrument
on the same channel cannot move the fluxes. The quadrature is the trapezoid
rule on each window's own grid — :func:`~ampere.backends.reference.bin_edges`'
midpoint bins clipped at the window's edges, so no window borrows a bin from
the gap beside it — accurate to order ``(h/λ)²``. Overlapping windows are
refused by name: a sample inside two windows would have two owners.

The transitions and their rest wavelengths are **taken from the observed
table**, never retyped (gap I-2's rule):

.. code-block:: python

    from ampere.backends.reference import LineIntegration
    from ampere.core import Instrument

    step = LineIntegration.from_observed(table, window=300.0)   # km/s
    instrument = Instrument([step], channel="lines", label="lines")

A table that carries its own integration windows as
``extra_coords["window"]`` (half-widths in ``window_unit``) needs no
``window=`` at all.

3. ``configure_from``: not needed
----------------------------------

The step needs nothing from any other step: its tabulation is fixed by the
line list and the windows. A chain that resampled or convolved the spectrum
first would put those steps before it, exactly as for synthetic photometry.

4. The composition: one model channel, one instrument
------------------------------------------------------

One model channel (a spectrum), one instrument, one dataset:

.. code-block:: python

    from ampere.core import Dataset, FittingProblem, GaussianFamily, IndependentNoise, Likelihood

    dataset = Dataset(table, instrument, likelihood=Likelihood(GaussianFamily(), IndependentNoise()))
    problem = FittingProblem(model, {"lines": dataset})

The table may list its lines in any order and on its own rounding of the rest
wavelengths: composition re-keys it onto the step's transitions and adopts
the step's rest wavelengths, once. An observed line the step does not
integrate, and a step line without an observation, are refused by name. This
is W7.15's photometric rule, generalised rather than copied — see §9.

5. The likelihood: the family default, and a kernel over rest wavelength
-------------------------------------------------------------------------

A line flux with an uncertainty is a Gaussian measurement, so the family is
:class:`~ampere.core.GaussianFamily`, the default for every kind (nothing in
code declares a family per kind). :class:`~ampere.core.PoissonFamily` is wrong
here: an integrated flux is a calibrated energy flux, not a count.

The flexible likelihood's kernel runs over the **rest wavelength**:

.. code-block:: python

    import scipy.stats as st
    from ampere.core import DenseGP, GaussianProcessNoise, Matern32

    kernel = Matern32(st.halfnorm(scale=0.5 * median_flux), st.uniform(1.0, 29.0),
                      axes=("spectral_axis",))
    noise = GaussianProcessNoise(kernel, DenseGP())

What a correlated residual across neighbouring transitions means is the
point of the modality: a model whose **excitation ladder** is wrong — one
temperature where the gas has two, a missing ionisation zone — mispredicts
lines in groups, not one at a time, because neighbouring transitions share
their excitation. The kernel absorbs that grouped residual instead of letting
it drag the physical parameters, and
:func:`~ampere.results.plot_gp_localisation` shows where along the ladder the
model failed. Ten lines is a small dense problem; ``QuasisepGP`` agrees with
``DenseGP`` on it (it sorts the axis itself), and is the one to reach for
with a few hundred lines.

6. What the conformance suite owes a new modality
--------------------------------------------------

``tests/conformance/test_line_fluxes.py`` runs on every registered backend
(behind a ``line_fluxes`` capability flag):

- **the step against a closed form** — a power law integrated over each
  window, ``∫ norm (λ/λ₀)^index c/λ² dλ``, which is analytic; the
  quadrature's error bound is the trapezoid rule's own,
  ``(b − a) h² max|g''| / 12``, derived in the row from
  ``samples_per_window`` (about 10⁻⁶ of the flux at 64 samples on a 0.6 µm
  window), and a second row shows four times the samples tightening the
  agreement more than tenfold;
- **the linear continuum is exact on a linear spectrum** — a straight-line
  model with ``continuum="linear"`` integrates to zero, to ``exact``
  relative to the uncorrected integral;
- **the published requirement** is ``points=`` at the tabulation, and the
  Gap-1 placement row: a resampling instrument beside the step on the same
  channel widens the negotiated grid and the fluxes do not change;
- **the mask**: a masked model sample inside one window masks that line and
  no other;
- **the composition**: a table in another order and in erg s⁻¹ cm⁻² aligns
  by transition label, and the two refusals name the transitions;
- **a draw row** — ``simulate`` produces a ``LineFluxes``;
- **the flexible likelihood** — ``DenseGP`` and ``QuasisepGP`` agree at
  ``cross_solver`` on a ``LineFluxes`` container whose axis is in no
  particular order;
- **the four diagnostics' plots** run on a stored run of the composition (a
  smoke row on the reference column);
- **the native realisation** — on torch and jax, the realised density equals
  the contract path's at ``cross_backend``, independent noise and the dense
  GP.

7. The native twins: inherited
-------------------------------

Both twins **inherit** the reference step, and override only the capability
flags, their placement, and the one ``matmul``:

.. code-block:: python

    class LineIntegration(_JaxStep, _ReferenceLineIntegration):
        ...  # apply_flux rebuilds the Axis and multiplies; the windows are inherited

The template's §8 rule decides it in its first clause: the declaration — the
per-line grids, the overlap refusal, the trapezoid weights, the continuum's
edge columns, whose union is the published ``points=`` — is the substantial,
dangerous half, and two backends tabulating it twice could negotiate
different grids; the arithmetic is one matrix product. On torch this departs
from ``ampere.backends.torch.instrument``'s own habit of re-declaring every
step, deliberately: the weights are built once in numpy by the inherited
constructor and registered a second time as a torch buffer, which is what a
device move and ``apply_flux`` see. The model in the example is written once
against :class:`~ampere.core.PortableModel`, so its twins are one line each.

8. The study: does the flexible likelihood survive a misspecified ladder?
--------------------------------------------------------------------------

The truth is a two-temperature ladder: a cool component (600 K) feeding every
line and a hot one (6000 K) feeding the high-excitation, short-wavelength
lines — two-thirds of [Ar II] 6.99 µm, a third of [Ar III] 8.99 µm, a few per
cent of the long-wavelength lines. The uncertainties are 3 per cent. The
line list is real (vacuum rest wavelengths from the fine-structure
compilation of Feuchtgruber et al. 1997 and the NIST database); **the ladder
is a toy** — the upper-level energies are ``30000 K µm / λ`` and the
strengths are chosen to make the lines comparably bright — because the
inference pattern, not the nebular physics, is what this page is about.
Three arms, as :doc:`m2_misspecification` has them:

- ``correct``: the two-temperature ladder with the hot component at its true,
  fixed value, independent noise — the control;
- ``misspecified``: one temperature, independent noise;
- ``flexible``: one temperature under the flexible likelihood.

Measured at the per-PR budget (emcee, 16 walkers, 600 steps with 300 burned;
about eleven seconds per arm on the reference backend, thirty-five for the
three) over three seeds — the posterior mean ± standard deviation of the cool
component's ``temperature`` (truth 600 K) and ``column`` (truth 3000, in
units of 10⁻¹⁶ W m⁻²), and the bias in posterior widths:

============  ========  ==================  ==================
arm           seed      temperature (K)     column
============  ========  ==================  ==================
correct       20261014  593 ± 5 (1.44)      3147 ± 80 (1.85)
correct       1         601 ± 6 (0.18)      3008 ± 84 (0.09)
correct       2         598 ± 6 (0.36)      3031 ± 78 (0.40)
misspecified  20261014  703 ± 6 (16.3)      2356 ± 55 (11.7)
misspecified  1         714 ± 7 (16.8)      2261 ± 53 (13.9)
misspecified  2         710 ± 7 (16.0)      2272 ± 57 (12.7)
flexible      20261014  720 ± 77 (1.57)     2017 ± 511 (1.93)
flexible      1         722 ± 66 (1.84)     2160 ± 485 (1.73)
flexible      2         733 ± 71 (1.87)     2179 ± 553 (1.49)
============  ========  ==================  ==================

The one-temperature ladder under independent noise is confidently wrong:
sixteen widths off in temperature, its interval six kelvin wide. Under the
flexible likelihood the same model's interval widens tenfold and the truth
sits within two widths at every seed, and the GP localisation peaks at
[Ar II] 6.99 µm — the highest-excitation line, where the hot component
contributes most — at all three. ``tests/examples/test_line_fluxes.py`` pins
those claims with the study's margin convention (each threshold a third
beyond the worst seed): the flexible arm and the control within 2.6 widths,
the misspecified arm beyond 12, the localisation's peak among the three
highest-excitation lines. NUTS on torch and jax runs the same flexible arm
(``heavy`` rows).

.. plot:: plots/line_fluxes.py
   :alt: The conditioned Gaussian-process mean of the one-temperature ladder's
      fit under the flexible likelihood, plotted against rest wavelength from 7
      to 26 microns at the ten transitions, with its 68 per cent band. It is
      largest at the shortest-wavelength line, [Ar II] 6.99, at about
      0.75e-16 W per square metre, drops to near zero between 10 and 14
      microns, and sits at a smaller positive offset of 0.2 to 0.4e-16 across
      the long-wavelength lines.

   The flexible arm's GP localisation along rest wavelength (seed 20261014):
   the residual the missing hot component leaves is largest at the
   high-excitation [Ar II] line; the smaller offset at long wavelengths is
   the other half of the one-temperature compromise, a slightly warm
   temperature and a low column under-predicting the cool lines.

9. What the template did not say
---------------------------------

- **The continuum is a declaration, and the template has no word for it.**
  A line flux is defined against a continuum someone chose. The step offers
  the two choices a linear operator can express (none; a straight line
  through the model's own edge values); anything cleverer — a fitted
  polynomial over side-bands — would make the step nonlinear in the samples
  and lose the one-matrix shape. The observer's choice should match the
  step's; the step cannot check that.
- **Identity, not position, was already a rule — for photometry only.**
  W7.15 had written "align by filter name" for one kind. A second named kind
  made it a property of kinds: ``FunctionSamples.IDENTITY``, with
  ``align_by_identity`` applying the rule to any kind that declares one and
  ``align_by_filter`` kept as the photometric spelling. A user-defined kind
  whose samples are named (a set of spectral indices, a list of molecular
  bands) gets the same alignment by declaring the attribute.
- **The window is part of the measurement.** A line table usually records
  the integration window; ``from_observed`` reads it from
  ``extra_coords["window"]`` when the table carries one, so the model is
  integrated over the window the data were.
- **The input unit is checked on one path only.** The native surface sees
  bare numbers. The weights carry the jansky's conversion, so the contract
  path refuses a model not in Jy at composition — where every problem runs
  it once — rather than letting the two paths integrate different numbers.
- **The real-file example the item hoped for does not exist yet.** W7.3's
  vendored JWST file is a MIRI LRS spectrum of a flux standard — a star with
  no emission lines — so there is no line list to read from it. A real fit
  needs a line-list reader (a table of labels, rest wavelengths, fluxes,
  uncertainties and windows) — a reader of the kind
  ``docs/design/horizon_beyond_phase7.md`` §5 surveys, none of which reads
  line tables yet — and until one lands, the hand-built table of §1 is the
  route.
- **A misspecification that projects onto a parameter is not fully
  rescued.** The hot component's excess falls off smoothly along the ladder,
  and part of it looks like a slightly warmer single temperature; the
  flexible arm absorbs the grouped residual and widens the interval, but its
  temperature still leans towards the warmer value (§8's table). The
  flexible likelihood helps against residual a parameter cannot imitate; it
  cannot undo the part of a misspecification that a parameter can.
