"""Run the SED composition example: ``python -m examples.sed_composition``.

    python -m examples.sed_composition                      # reference, emcee
    python -m examples.sed_composition --backend torch       # NUTS
    python -m examples.sed_composition --walkers 32 --steps 200

Restricts every BLAS thread pool to one thread **before numpy is imported
anywhere in this process**, which is why this is a separate, minimal module
rather than a ``if __name__ == "__main__":`` guard at the foot of
:mod:`.sed_composition`. The composition's grids are a couple of thousand
points — too small for multi-threaded BLAS to pay for its own thread-pool
overhead — and on this example, measured with twenty calls to
``problem.log_prob`` at the truth (this docstring's own design guidance),
restricting to one thread roughly halves the per-evaluation time (measured in
the high teens of milliseconds under the default, multi-threaded BLAS;
roughly half that restricted to one thread — see :doc:`the tutorial page
</sed_composition>` for the numbers from one specific run), which is the
difference between comfortably and uncomfortably inside W4.11's two-minute
budget. It has no effect on the *answer*, only on how many CPU threads
compute it, and it is scoped to a subprocess invocation of this module, never
imposed on a caller who imports :mod:`.sed_composition` as a library.
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

from .sed_composition import main  # noqa: E402 -- must follow the thread-pool env vars above

if __name__ == "__main__":
    raise SystemExit(main())
