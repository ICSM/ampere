"""Line fluxes as an observable: a table of integrated fluxes, fitted (W7.9).

The modality template applied to the measurement a spectroscopist most often
has — a set of integrated emission-line fluxes with uncertainties, one per
transition — expanded into a runnable example and :doc:`the line-fluxes
tutorial page </line_fluxes>`. :mod:`.generators` holds the synthetic,
JWST-like line list and its two-temperature truth, :mod:`.model` the
excitation-ladder model (with its one-line torch and jax twins), and
:mod:`.study` the three arms, the fit and the thresholds the tests assert.

Deliberately minimal: ``__main__.py`` restricts every BLAS thread pool to one
thread before numpy is imported, so this file imports neither submodule.
"""
