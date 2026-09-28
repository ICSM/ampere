"""One ``ModifiedBlackBody``, three datasets: two spectrographs and a photometric
catalogue on the same model channel, with a per-spectrograph calibration factor.

``docs/design/modalities/spectrum_photometry.md`` §4 leaves the "several
observations" case to an example; this package is that example (W6.2), the
sibling :mod:`examples.sed_composition` extends to three observations of one
source instead of two, and :doc:`the tutorial page </photometry_spectra>`
walks through it. See :mod:`.generators` for how the synthetic data is made
(the same negotiate-then-``compile_for`` dance
:mod:`examples.sed_composition.generators` uses, extended to three
instruments and an injected smooth residual one of them cannot express), and
:mod:`.photometry_spectra` for the problem, the tie, the flexible likelihood
and the CLI.

Deliberately minimal, for the same reason as the sibling package:
``__main__.py`` restricts every BLAS thread pool to one thread before numpy
is imported anywhere in the process, and this file is imported *before*
``__main__.py`` runs — so it does not import :mod:`.generators` or
:mod:`.photometry_spectra` (both of which import numpy) eagerly. Import the
submodules directly: ``from examples.photometry_spectra import generators``
or ``from examples.photometry_spectra.photometry_spectra import build_problem``.
"""

from __future__ import annotations
