"""Run the study and print its tables: ``python -m examples.m2_misspecification``.

One command, the same code the tests and the benchmarks call, and nothing
written anywhere unless ``--figures`` names a directory.

    python -m examples.m2_misspecification                 # 200 points, reference
    python -m examples.m2_misspecification --size 20000    # the top of the ladder
    python -m examples.m2_misspecification --backend jax   # NUTS, if the extra is there
    python -m examples.m2_misspecification --figures /tmp/m2

The default budget is the milestone one
(:data:`~examples.m2_misspecification.study.MILESTONE_EMCEE`), which is what
the numbers in ``docs/source/m2_misspecification.rst`` were measured at and
takes a few minutes at 200 points; ``--quick`` swaps in ``tests/m2``'s much
shorter budget for a first look.
"""

from __future__ import annotations

import argparse
import sys
import time
from collections.abc import Sequence

from . import study
from .generators import SIZES
from .model import PARAMETER_NAMES


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="python -m examples.m2_misspecification")
    parser.add_argument("--size", type=int, default=SIZES[0], help="points in the spectrum")
    parser.add_argument("--backend", default="reference", choices=list(study.BACKENDS))
    parser.add_argument(
        "--solver",
        default="quasisep",
        choices=list(study.SOLVERS),
        help="GP solver for the flexible likelihood; 'dense' is O(N^3)",
    )
    parser.add_argument(
        "--kernel", default="matern32", choices=list(study.KERNELS), help="GP kernel"
    )
    parser.add_argument("--quick", action="store_true", help="use tests/m2's short budget")
    parser.add_argument("--figures", default=None, help="write the figures into this directory")
    parser.add_argument("--seed", type=int, default=study.SEED, help="the run seed")
    return parser


def _budget(backend: str, quick: bool) -> study.EmceeBudget | study.NutsBudget:
    if backend == "reference":
        return study.TEST_EMCEE if quick else study.MILESTONE_EMCEE
    return study.TEST_NUTS if quick else study.MILESTONE_NUTS


def main(argv: Sequence[str] | None = None) -> int:
    """Run the study, print the recovery and diagnostics tables, return an exit code."""
    arguments = _parser().parse_args(argv)
    budget = _budget(arguments.backend, arguments.quick)
    print(
        f"# M2: {arguments.backend} backend, {arguments.size} points, "
        f"{arguments.kernel}/{arguments.solver}, {budget}"
    )
    started = time.perf_counter()
    results = study.run_study(
        size=arguments.size,
        backend=arguments.backend,
        budget=budget,
        kernel=arguments.kernel,
        solver=arguments.solver,
        seed=arguments.seed,
    )
    print(f"# sampled in {time.perf_counter() - started:.1f} s")
    diagnoses = study.prepare(results)

    header = " ".join(f"{name:>9s}" for name in PARAMETER_NAMES)
    print(f"\n{'scenario':<16s}{'likelihood':<11s}{header}   covers  diagnostic")
    print("-" * 96)
    for (key, kind), entry in results.items():
        summaries = study.summarise(entry["run"], names=study.PHYSICAL_NAMES)
        offsets = " ".join(
            f"{summaries[name].bias_in_widths:9.2f}" for name in study.PHYSICAL_NAMES
        )
        covers = all(summaries[name].covers_truth for name in study.PHYSICAL_NAMES)
        diagnosis = diagnoses[(key, kind)]
        if diagnosis.whiteness_p_value is not None:
            note = (
                f"whiteness p = {diagnosis.whiteness_p_value:.4f} "
                f"(Q = {diagnosis.whiteness_statistic:.1f})"
            )
        else:
            note = (
                f"GP peak at {diagnosis.localisation_peak:.5f} um, "
                f"score {diagnosis.localisation_max:.2f}"
            )
        print(f"{key:<16s}{kind:<11s}{offsets}   {covers!s:<7s} {note}")
    print(
        "\nThe four columns are |median - truth| in units of the posterior's own 68 % half-width."
    )

    print("\nGP hyperparameters, flexible likelihood:")
    for (key, kind), entry in results.items():
        if kind != "flexible":
            continue
        names = [
            f"{study.DATASET_LABEL}.likelihood.amplitude",
            f"{study.DATASET_LABEL}.likelihood.length_scale",
        ]
        summaries = study.summarise(entry["run"], names=names)
        amplitude, length = (summaries[name] for name in names)
        print(
            f"  {key:<16s} amplitude = {amplitude.median:.5g} +- {amplitude.width:.3g} Jy, "
            f"length scale = {length.median:.5g} +- {length.width:.3g} um"
        )

    if arguments.figures is not None:
        from . import figures

        written = figures.save_paper_figures(results, arguments.figures)
        written += figures.save_result_figures(results, arguments.figures)
        print(f"\nWrote {len(written)} figures to {arguments.figures}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
