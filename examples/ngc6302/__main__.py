"""Run the NGC 6302 example: ``python -m examples.ngc6302``.

    python -m examples.ngc6302                          # emcee, real data, full legacy budget
    python -m examples.ngc6302 --quick                   # emcee, real data, 2000/1000
    python -m examples.ngc6302 --synthetic --engine zeus  # a look at the other engine

Restricts every BLAS thread pool to one thread **before numpy is imported
anywhere in this process**, the same reasoning as
:mod:`examples.linear_sed` and :mod:`examples.modified_blackbody`'s
``__main__.py`` modules. Scoped to a subprocess invocation of this module,
never imposed on a caller who imports :mod:`.ngc6302` as a library.
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

from .ngc6302 import main  # noqa: E402 -- must follow the thread-pool env vars above

if __name__ == "__main__":
    raise SystemExit(main())
