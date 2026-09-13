"""``ReflexOrbit`` on two ``TimeSeries`` channels: the second worked modality (W4.9).

Memo-free this time — the design sketch is
``docs/design/modalities/astrometric_timeseries.md`` and the template it
tests is :doc:`the interferometry page </interferometry>` — expanded into a
runnable example and :doc:`the astrometry tutorial page </astrometry>`. See
:mod:`.generators` for the synthetic data and :mod:`.astrometry` for the
problem, the fit, and the backend flag.

Deliberately minimal, in :mod:`examples.sed_composition`'s own style:
``__main__.py`` restricts every BLAS thread pool to one thread before numpy
is imported anywhere in the process, so this file does not eagerly import
either submodule — ``from examples.astrometry import generators`` or
``from examples.astrometry.astrometry import build_problem``.
"""

from __future__ import annotations
