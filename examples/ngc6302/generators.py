"""The real spectrum, the opacity-table locations, and the synthetic truth.

``examples.ngc6302``'s data-loading half (W6.13 (2)) -- see :mod:`.ngc6302`
for the model and the fitting problem.

The real data
-------------
:func:`load_observed_spectrum` reads the tracked
``examples/NGC6302/NGC6302_100.tab`` (two columns, wavelength in micron and
flux in Jy, a two-line header) with :func:`numpy.loadtxt` -- no legacy
reader involved -- and applies legacy's own 25-120 micron selection
(``examples/NGC6302.py`` line ~380, ``spec.selectWaves(low=25, up=120)``).

The uncertainty rule
--------------------
The ``.tab`` file carries no uncertainty column. Legacy's own ``__main__``
sets it at construction (``examples/NGC6302.py`` lines ~360-368)::

    unc = specdata[1][:] * input_noise_spec   # input_noise_spec = 0.01
    spec = Spectrum(specdata[0][:], specdata[1][:] + randn(...) * unc,
                    specdata[1][:] * 0.05, "um", "Jy", ...)

The **third positional argument is the stored uncertainty**: ``0.05`` of the
flux, five per cent -- a different number from the ``0.01`` (one per cent)
used only to jitter the demonstration flux that precedes it. This twin
reproduces the *stored* rule (:data:`FRACTIONAL_UNCERTAINTY`, five per cent
of the flux) for both the real spectrum's uncertainty and, in
:func:`synthetic_data`, the synthetic noise draw -- the legacy script's own
separate one-per-cent jitter has no continuing role once the twin generates
its own synthetic truth from the model rather than perturbing the file's
real flux.

The opacity tables and the synthetic truth
-------------------------------------------
:mod:`.ngc6302` resolves the eight opacity tables and the observed data file
by walking up from ``__file__`` (:data:`OPACITY_DIRECTORY`,
:data:`OPACITY_FILE_LIST`), the same convention
:mod:`examples.linear_sed.generators` uses for its own tracked CASSIS file,
rather than ``importlib.resources`` -- ``examples/NGC6302/`` is tracked data
belonging to the legacy example and is not moved, copied or turned into a
package.

:data:`TRUTH` is the "2002 solution" ruling 3 asks for: not literally Kemper
et al. (2002)'s published table, but ``examples/NGC6302-calculate-dust-mass.py``'s
own hard-coded ``n0``/``Tin``/``Tout`` -- a previous fit's result, reported
in ``examples/NGC6302/table.tex`` as "present work" -- qualified onto the
model's own fifteen parameter names via the species correspondence below
(the dust-mass script's own ``names``/``rhod`` order does not match the
model's species index order, so the mapping is stated explicitly rather than
assumed):

======================  =================================  ==================
species                 dust-mass script's ``n0`` index     model parameter
======================  =================================  ==================
calcite                 cold[2]                              ``logacold0``
enstatite               cold[5]                               ``logacold1``
forsterite              cold[1]                               ``logacold2``
diopside                cold[4]                               ``logacold3``
dolomite                cold[6]                               ``logacold4``
ice                     cold[3]                               ``logacold6``
olivine                 cold[0]                               ``logacold7``
enstatite               warm[5]                               ``logawarm1``
forsterite              warm[1]                               ``logawarm2``
iron                    warm[7]                               ``logawarm5``
olivine                 warm[0]                               ``logawarm7``
======================  =================================  ==================

and ``Tcold0/Tcold1 = Tout/Tin[0]``, ``Twarm0/Twarm1 = Tout/Tin[1]`` directly.
:mod:`.dust_mass` uses the identical correspondence in the other direction.
"""

from __future__ import annotations

from pathlib import Path

__all__ = [
    "DATA_FILE",
    "FRACTIONAL_UNCERTAINTY",
    "OPACITY_DIRECTORY",
    "OPACITY_FILE_LIST",
    "SEED",
    "SPECIES",
    "TRUTH",
    "WAVELENGTH_SELECTION",
]

#: Reproducible everywhere this example is run, in the convention
#: :mod:`examples.linear_sed` and :mod:`examples.modified_blackbody` also use.
#: Legacy seeds nothing (``np.random.randn`` unseeded in ``NGC6302.py``).
SEED = 20260928

#: The eight opacity tables' own order (``examples/NGC6302-opacities.txt``),
#: and the model's own species index (0-7) throughout this package and
#: :mod:`.ngc6302`.
SPECIES = ("calcite", "enstatite", "forsterite", "diopside", "dolomite", "iron", "ice", "olivine")

#: Resolved by walking up from this file, as
#: :data:`examples.linear_sed.generators.IRS_FILE` resolves its own tracked
#: data file -- not ``importlib.resources``, since ``examples/NGC6302/`` is
#: the legacy example's own tracked directory and stays exactly where it is
#: (the item's own instruction: "not moved or copied").
OPACITY_DIRECTORY = Path(__file__).resolve().parents[1] / "NGC6302"
OPACITY_FILE_LIST = Path(__file__).resolve().parents[1] / "NGC6302-opacities.txt"
DATA_FILE = OPACITY_DIRECTORY / "NGC6302_100.tab"

#: Legacy's own selection (``examples/NGC6302.py`` line ~380).
WAVELENGTH_SELECTION = (25.0, 120.0)

#: The stored-uncertainty rule -- five per cent of the flux (see the module
#: docstring for the one-per-cent jitter this is *not*).
FRACTIONAL_UNCERTAINTY = 0.05

#: The "2002 solution" -- see the module docstring for the species
#: correspondence and its source
#: (``examples/NGC6302-calculate-dust-mass.py`` lines 20-34, 87-88).
TRUTH: dict[str, float] = {
    "logacold0": -4.44527,  # calcite
    "logacold1": -5.03878,  # enstatite
    "logacold2": -2.95189,  # forsterite
    "logacold3": -4.23599,  # diopside
    "logacold4": -4.89447,  # dolomite
    "logacold6": -4.01175,  # ice
    "logacold7": -1.07680,  # olivine
    "logawarm1": -5.49659,  # enstatite
    "logawarm2": -5.57701,  # forsterite
    "logawarm5": -3.70948,  # iron
    "logawarm7": -3.72738,  # olivine
    "Tcold0": 36.13010,
    "Tcold1": 57.03042,
    "Twarm0": 105.19534,
    "Twarm1": 122.69678,
}
