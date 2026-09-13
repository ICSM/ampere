"""Run the interferometry study: ``python -m examples.interferometry``.

    python -m examples.interferometry                  # 3 arms, reference, emcee
    python -m examples.interferometry --chromatic       # arm (d) instead
    python -m examples.interferometry --backend torch   # NUTS, one arm
    python -m examples.interferometry --calibration     # the SBC coverage row
    python -m examples.interferometry --figures /tmp/itf

Restricts every BLAS thread pool to one thread before numpy is imported,
``examples/sed_composition/__main__.py``'s own reason: this study's image
grids are a few pixels across (a handful of ``(u, v)`` points DFT'd against
them), too small for multi-threaded BLAS to repay its own overhead.
"""

from __future__ import annotations

import os

for _threads in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_threads, "1")

import argparse  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from collections.abc import Sequence  # noqa: E402

from . import study  # noqa: E402


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="python -m examples.interferometry")
    parser.add_argument("--backend", default="reference", choices=list(study.BACKENDS))
    parser.add_argument("--arms", nargs="+", default=list(study.ARMS), choices=list(study.ARMS))
    parser.add_argument("--chromatic", action="store_true", help="run arm (d) instead of (a)-(c)")
    parser.add_argument("--calibration", action="store_true", help="run the SBC coverage row instead")
    parser.add_argument("--doc-budget", action="store_true", help="use the longer documentation budget")
    parser.add_argument("--figures", default=None, help="write the figures into this directory")
    parser.add_argument("--seed", type=int, default=study.gen.SEED, help="the data/run seed")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the study, print the recovery table (or the calibration row), return an exit code."""
    arguments = _parser().parse_args(argv)
    started = time.perf_counter()

    if arguments.calibration:
        print(f"# interferometry calibration: {arguments.backend} backend, arms {arguments.arms}")
        calibrations = {arm: study.run_calibration(arm, backend=arguments.backend, seed=arguments.seed) for arm in arguments.arms}
        print(f"# ran in {time.perf_counter() - started:.1f} s")
        for arm, calibration in calibrations.items():
            coverage = study.coverage_at(calibration, 0.9)
            print(f"{arm:<12s} coverage@0.9 = {coverage} failures={int(calibration.attrs['ampere_calibration_failures'])}")
        if arguments.figures is not None:
            from . import figures

            written = figures.save_calibration_figures(calibrations, arguments.figures)
            print(f"\nWrote {len(written)} figures to {arguments.figures}")
        return 0

    if arguments.chromatic:
        budget = study.DOC_BUDGET if arguments.doc_budget else study.CHROMATIC_BUDGET
        print(f"# interferometry chromatic arm: {arguments.backend} backend, {budget}")
        results = study.run_chromatic_arm(backend=arguments.backend, budget=budget, seed=arguments.seed)
        print(f"# sampled in {time.perf_counter() - started:.1f} s")
        _print_table(results)
        if arguments.figures is not None:
            from . import figures

            written = figures.save_arm_figures(results, arguments.figures)
            print(f"\nWrote {len(written)} figures to {arguments.figures}")
        return 0

    budget = study.DOC_BUDGET if arguments.doc_budget else study.CI_BUDGET
    print(f"# interferometry study: {arguments.backend} backend, arms {arguments.arms}, {budget}")
    results = study.run_study(backend=arguments.backend, arms=arguments.arms, budget=budget, seed=arguments.seed)
    print(f"# sampled in {time.perf_counter() - started:.1f} s")
    _print_table(results)
    if arguments.figures is not None:
        from . import figures

        written = figures.save_arm_figures(results, arguments.figures)
        print(f"\nWrote {len(written)} figures to {arguments.figures}")
    return 0


def _print_table(results: dict[str, dict]) -> None:
    names = list(study.gen.TRUTH)
    header = " ".join(f"{name:>16s}" for name in names)
    print(f"\n{'arm/kernel':<12s}{header}")
    print("-" * (12 + 17 * len(names)))
    for key, entry in results.items():
        summaries = study.summarise(entry["run"], names=names)
        row = " ".join(f"{summaries[name].bias_in_widths:16.2f}" for name in names)
        print(f"{key:<12s}{row}")
    print("\nEach column is |median - truth| in units of the posterior's own 68 % half-width.")


if __name__ == "__main__":
    sys.exit(main())
