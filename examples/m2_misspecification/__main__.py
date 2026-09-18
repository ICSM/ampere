"""Run the study and print its tables: ``python -m examples.m2_misspecification``.

One command, the same code the tests and the benchmarks call, and nothing
written anywhere unless ``--figures`` names a directory.

    python -m examples.m2_misspecification                 # 200 points, reference
    python -m examples.m2_misspecification --size 20000    # the top of the ladder
    python -m examples.m2_misspecification --backend jax   # NUTS, if the extra is there
    python -m examples.m2_misspecification --figures /tmp/m2

W5.8's extension is reached from here too — one scenario, one kernel at a time::

    python -m examples.m2_misspecification --scenario many_lines --kernel matern32
    python -m examples.m2_misspecification --scenario many_lines --kernel warped
    python -m examples.m2_misspecification --scenario many_lines --kernel sum

and the comparison of the three, side by side in one table, is
``python -m examples.m2_misspecification.many_lines``.

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
from .generators import EXTENDED_SCENARIOS, SIZES
from .model import PARAMETER_NAMES


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="python -m examples.m2_misspecification")
    parser.add_argument("--size", type=int, default=SIZES[0], help="points in the spectrum")
    parser.add_argument("--backend", default="reference", choices=list(study.BACKENDS))
    parser.add_argument(
        "--scenario",
        default=None,
        choices=[scenario.key for scenario in EXTENDED_SCENARIOS],
        help="run one scenario rather than the four; 'many_lines' is W5.8's",
    )
    parser.add_argument(
        "--solver",
        default="quasisep",
        choices=list(study.SOLVERS),
        help="GP solver for the flexible likelihood; 'dense' is O(N^3)",
    )
    parser.add_argument(
        "--kernel",
        default="matern32",
        choices=list(study.KERNELS),
        help="GP kernel; 'warped' is W5.7's and 'sum' is W4.5's",
    )
    parser.add_argument("--quick", action="store_true", help="use tests/m2's short budget")
    parser.add_argument("--figures", default=None, help="write the figures into this directory")
    parser.add_argument("--seed", type=int, default=study.SEED, help="the run seed")
    return parser


def _budget(backend: str, quick: bool, kernel: str) -> study.EmceeBudget | study.NutsBudget:
    """The sampler budget for this backend, this kernel and this ``--quick``.

    W5.8's two arms are why *kernel* is an argument. ``warped`` carries six
    more free parameters than the stationary kernel and ``sum`` two, and an
    emcee ensemble needs more than twice the dimension in walkers — so
    ``TEST_EMCEE``'s twenty cannot run either of them, and
    :data:`~examples.m2_misspecification.study.MANY_LINES_EMCEE` is the budget
    sized for the widest of them. It is also short enough to serve as the quick
    look, so both arms get it either way.
    """
    if backend != "reference":
        return study.TEST_NUTS if quick else study.MILESTONE_NUTS
    if kernel in ("warped", "sum"):
        return study.MANY_LINES_EMCEE
    return study.TEST_EMCEE if quick else study.MILESTONE_EMCEE


def main(argv: Sequence[str] | None = None) -> int:
    """Run the study, print the recovery and diagnostics tables, return an exit code."""
    arguments = _parser().parse_args(argv)
    budget = _budget(arguments.backend, arguments.quick, arguments.kernel)
    chosen = None if arguments.scenario is None else [arguments.scenario]
    print(
        f"# M2: {arguments.backend} backend, {arguments.size} points, "
        f"{arguments.kernel}/{arguments.solver}, "
        f"{'the four scenarios' if chosen is None else arguments.scenario}, {budget}"
    )
    started = time.perf_counter()
    results = study.run_study(
        size=arguments.size,
        backend=arguments.backend,
        scenarios=chosen,
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
    prefix = f"{study.DATASET_LABEL}.likelihood."
    for (key, kind), entry in results.items():
        if kind != "flexible":
            continue
        # Read the names off the problem rather than spelling them. W5.8's
        # ``warped`` and ``sum`` arms declare a different set from the
        # stationary kernel's amplitude-and-length-scale pair — six warp
        # variables, or two labelled terms — and a hard-coded pair is a
        # KeyError the moment --kernel names one of them.
        declared = entry["problem"].parameters
        names = [name for name in declared.free_names if name.startswith(prefix)]
        summaries = study.summarise(entry["run"], names=names)
        print(f"  {key}:")
        for name in names:
            summary = summaries[name]
            unit = declared[name].unit
            suffix = "" if unit is None else f" {unit}"
            print(
                f"    {name[len(prefix) :]:<26s} "
                f"{summary.median:>12.5g} +- {summary.width:.3g}{suffix}"
            )

    if arguments.figures is not None:
        from . import figures

        written = figures.save_paper_figures(results, arguments.figures)
        written += figures.save_result_figures(results, arguments.figures)
        print(f"\nWrote {len(written)} figures to {arguments.figures}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
