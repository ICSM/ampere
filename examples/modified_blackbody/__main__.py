"""Run the modified blackbody example: ``python -m examples.modified_blackbody``.

    python -m examples.modified_blackbody                              # emcee
    python -m examples.modified_blackbody --engine dynesty
    python -m examples.modified_blackbody --all --figures /tmp/mbb-figures

Restricts every BLAS thread pool to one thread **before numpy is imported
anywhere in this process**, the same reasoning as
:mod:`examples.sed_composition` and :mod:`examples.linear_sed`'s
``__main__.py`` modules. Scoped to a subprocess invocation of this module,
never imposed on a caller who imports :mod:`.modified_blackbody` as a
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

from .modified_blackbody import main  # noqa: E402 -- must follow the thread-pool env vars above

if __name__ == "__main__":
    raise SystemExit(main())
