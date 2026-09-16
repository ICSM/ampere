"""Run the image study: ``python -m examples.image``.

    python -m examples.image                        # 3 arms, reference, emcee
    python -m examples.image --calibration          # the SBC coverage row
    python -m examples.image --benchmark            # DenseGP vs HSGP at three N
    python -m examples.image --benchmark --no-dense # ... without the ~10 GiB cell
    python -m examples.image --backend torch        # NUTS
    python -m examples.image --figures /tmp/img

Restricts every BLAS thread pool to one thread before numpy is imported, as
every other example's runner here does — except under ``--benchmark``, where
the whole point is the cost of a large dense solve and pinning BLAS to one
thread would measure something nobody runs. ``--benchmark`` therefore leaves
the environment alone, and says so in its own header line.
"""

from __future__ import annotations

import os
import sys

_BENCHMARKING = "--benchmark" in sys.argv[1:]
if not _BENCHMARKING:
    for _threads in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ.setdefault(_threads, "1")

import argparse  # noqa: E402
import time  # noqa: E402
from collections.abc import Sequence  # noqa: E402

from . import study  # noqa: E402


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="python -m examples.image")
    parser.add_argument("--backend", default="reference", choices=list(study.BACKENDS))
    parser.add_argument("--arms", nargs="+", default=list(study.ARMS), choices=list(study.ARMS))
    parser.add_argument(
        "--calibration", action="store_true", help="run the SBC coverage row instead"
    )
    parser.add_argument(
        "--benchmark", action="store_true", help="run the solver cost table instead"
    )
    parser.add_argument(
        "--sizes",
        nargs="+",
        type=int,
        default=list(study.BENCHMARK_SIZES),
        help="image sizes (pixels on a side) for --benchmark",
    )
    parser.add_argument(
        "--no-dense",
        action="store_true",
        help=(
            "skip the exact DenseGP cells of --benchmark: its peak is five copies of an "
            "N x N covariance, which at 128x128 is about 10 GiB"
        ),
    )
    parser.add_argument(
        "--doc-budget", action="store_true", help="use the longer documentation budget"
    )
    parser.add_argument("--pixels", type=int, default=study.SMALL_PIXELS, help="fitted image size")
    parser.add_argument("--figures", default=None, help="write the figures into this directory")
    parser.add_argument("--seed", type=int, default=study.gen.SEED, help="the data/run seed")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the study, print the table it produces, and return an exit code."""
    arguments = _parser().parse_args(argv)
    started = time.perf_counter()

    if arguments.benchmark:
        print(f"# image solver benchmark: {arguments.backend} backend, BLAS threads unrestricted")
        costs = study.benchmark_solvers(
            sizes=arguments.sizes,
            backend=arguments.backend,
            seed=arguments.seed,
            include_dense=not arguments.no_dense,
        )
        print(f"# measured in {time.perf_counter() - started:.1f} s\n")
        print(study.benchmark_table(costs))
        print(
            "\nInformational (DEVELOPMENT_PLAN.md §5's 'chosen by measurement'): one log_prob "
            "each,\nthe same likelihood and the same image, so the difference is the solve."
        )
        return 0

    if arguments.calibration:
        print(f"# image calibration: {arguments.backend} backend, arms {arguments.arms}")
        calibrations = {
            arm: study.run_calibration(
                arm, backend=arguments.backend, pixels=arguments.pixels, seed=arguments.seed
            )
            for arm in arguments.arms
        }
        print(f"# ran in {time.perf_counter() - started:.1f} s")
        for arm, calibration in calibrations.items():
            coverage = study.coverage_at(calibration, 0.9)
            failures = int(calibration.attrs["ampere_calibration_failures"])
            print(f"{arm:<12s} coverage@0.9 = {coverage} failures={failures}")
        if arguments.figures is not None:
            from . import figures

            written = figures.save_calibration_figures(calibrations, arguments.figures)
            print(f"\nWrote {len(written)} figures to {arguments.figures}")
        return 0

    budget = study.DOC_BUDGET if arguments.doc_budget else study.CI_BUDGET
    print(
        f"# image study: {arguments.backend} backend, arms {arguments.arms}, "
        f"{arguments.pixels}x{arguments.pixels} pixels, {budget}"
    )
    results = study.run_study(
        backend=arguments.backend,
        arms=arguments.arms,
        budget=budget,
        pixels=arguments.pixels,
        seed=arguments.seed,
    )
    print(f"# sampled in {time.perf_counter() - started:.1f} s")
    _print_table(results)
    if arguments.figures is not None:
        from . import figures

        written = figures.save_arm_figures(results, arguments.figures)
        written += figures.save_image_panels(results, arguments.figures)
        print(f"\nWrote {len(written)} figures to {arguments.figures}")
    return 0


def _print_table(results: dict[str, dict]) -> None:
    names = list(study.gen.TRUTH)
    header = " ".join(f"{name:>16s}" for name in names)
    print(f"\n{'arm':<12s}{header}")
    print("-" * (12 + 17 * len(names)))
    for key, entry in results.items():
        summaries = study.summarise(entry["run"], names=names)
        row = " ".join(f"{summaries[name].bias_in_widths:16.2f}" for name in names)
        print(f"{key:<12s}{row}")
    print("\nEach column is |median - truth| in units of the posterior's own 68 % half-width.")


if __name__ == "__main__":
    sys.exit(main())
