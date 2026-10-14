"""Run the line-fluxes study: ``python -m examples.line_fluxes``.

    python -m examples.line_fluxes                    # reference, emcee, the docs budget
    python -m examples.line_fluxes --quick            # the per-PR budget
    python -m examples.line_fluxes --backend jax      # NUTS on the native path
    python -m examples.line_fluxes --seed 7

Restricts every BLAS thread pool to one thread **before numpy is imported
anywhere in this process**, as :mod:`examples.astrometry`'s own
``__main__.py`` does: a ten-line problem is far too small for multi-threaded
BLAS to pay for its thread-pool overhead.
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

from .study import main  # noqa: E402 -- must follow the thread-pool env vars above

if __name__ == "__main__":
    raise SystemExit(main())
