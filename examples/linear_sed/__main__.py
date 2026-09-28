"""Run the linear SED example: ``python -m examples.linear_sed``.

    python -m examples.linear_sed                            # emcee, reference
    python -m examples.linear_sed --engine zeus
    python -m examples.linear_sed --engine dynesty
    python -m examples.linear_sed --engine sbi --embedding    # pixi run -e sbi

Restricts every BLAS thread pool to one thread **before numpy is imported
anywhere in this process** -- the same reasoning as
:mod:`examples.sed_composition`'s ``__main__.py``: this composition's grids
are a couple of thousand points, too small for multi-threaded BLAS to pay for
its own thread-pool overhead. Scoped to a subprocess invocation of this
module, never imposed on a caller who imports :mod:`.linear_sed` as a
library.
"""

from __future__ import annotations

import os

for _threads in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_threads, "1")

from .linear_sed import main  # noqa: E402 -- must follow the thread-pool env vars above

if __name__ == "__main__":
    raise SystemExit(main())
