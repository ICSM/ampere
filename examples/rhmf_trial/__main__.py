"""Run the RHMF trial: ``python -m examples.rhmf_trial --out DIR [--quick]``.

    pixi run -e rhmf python -m examples.rhmf_trial --quick --out /tmp/rhmf
    pixi run -e rhmf python -m examples.rhmf_trial --out /tmp/rhmf
    pixi run -e rhmf python -m examples.rhmf_trial --only image --out /tmp/rhmf

Prints the per-scenario tables and writes them, with the figures, into ``DIR``
(created if need be, and never inside git: ground rule 7). ``--quick`` is a
two-minute smoke preset; the default is the preset the findings in
``diagnostics.md`` §7 were taken from.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import os
import pathlib
import sys
import time
from collections.abc import Sequence
from typing import Any

for _threads in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_threads, "4")


@dataclasses.dataclass(frozen=True)
class Preset:
    """What one run covers."""

    ranks: tuple[int, ...]
    scales: tuple[float, ...]
    controls: int
    copies: tuple[int, ...]  # deviated rows per scenario, one collection each
    headline_copies: int
    max_iter: int
    image_pixels: int
    image_count: int
    backgrounds: tuple[float, ...]  # multiples of W5.5's background flux


QUICK = Preset((1, 2, 3), (1.0, 3.0), 6, (2,), 2, 150, 24, 8, (1.0, 5.0))
FULL = Preset((1, 2, 3, 4), (1.0, 2.0, 3.0, 5.0), 12, (1, 3, 6), 3, 500, 24, 24, (1.0, 5.0))


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="python -m examples.rhmf_trial", description=__doc__)
    parser.add_argument("--out", required=True, help="run directory for tables and figures")
    parser.add_argument("--quick", action="store_true", help="the under-two-minutes preset")
    parser.add_argument("--only", choices=("spectra", "image"), default=None)
    parser.add_argument("--seed", type=int, default=0)
    return parser


def _write_csv(path: pathlib.Path, records: Sequence[dict[str, Any]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def _print_table(rows: Sequence[dict[str, Any]], columns: Sequence[str]) -> None:
    widths = [max(len(c), 9) for c in columns]
    print(" ".join(f"{c:>{w}s}" for c, w in zip(columns, widths, strict=True)))
    for row in rows:
        cells = []
        for column, width in zip(columns, widths, strict=True):
            value = row[column]
            cells.append(
                f"{value:>{width}.3f}" if isinstance(value, float) else f"{value!s:>{width}s}"
            )
        print(" ".join(cells))


def run_spectra(preset: Preset, out: pathlib.Path, seed: int) -> list[dict[str, Any]]:
    """Trial (a); returns every record."""
    from . import figures, trial

    all_records: list[dict[str, Any]] = []
    for copies in preset.copies:
        collection = trial.spectra_collection(controls=preset.controls, copies=copies, seed=seed)
        records, fits = trial.scan_spectra(
            collection, ranks=preset.ranks, scales=preset.scales, max_iter=preset.max_iter
        )
        for record in records:
            record["copies"] = copies
        all_records += records
        best = trial.best_by_scenario(records)
        print(
            f"\n== spectra: {preset.controls} controls + {copies} row(s) of each of "
            f"{len(collection.bands)} scenarios ({len(collection.labels)} rows x "
            f"{collection.Y.shape[1]} features); best grid point per scenario "
            f"(oracle choice: the deviation's location is known)"
        )
        _print_table(
            list(best.values()),
            (
                "scenario",
                "rank",
                "robust_scale",
                "band_points",
                "ratio",
                "control_ratio",
                "excess",
                "lift",
                "auc_rows",
            ),
        )
        if copies == preset.headline_copies:
            for path in figures.save_spectra_figures(collection, fits, best, out):
                print(f"   wrote {path}")
    _write_csv(out / "spectra_grid.csv", all_records)
    return all_records


def run_image(preset: Preset, out: pathlib.Path, seed: int) -> list[dict[str, Any]]:
    """Trial (b), at each background multiple; returns every record."""
    from . import figures, trial

    all_records: list[dict[str, Any]] = []
    for background in preset.backgrounds:
        collection = trial.image_collection(
            pixels=preset.image_pixels, count=preset.image_count, seed=seed, background=background
        )
        print(
            f"\n== image: {preset.image_count} images of {preset.image_pixels}x"
            f"{preset.image_pixels} (half with the omitted background at {background:g}x W5.5's "
            f"flux; peak {collection.peak_signal_to_noise:.2f} sigma per pixel); band = "
            f"{int(collection.band.sum())} of {collection.band.size} pixels"
        )
        started = time.perf_counter()
        records_a, fits_a = trial.scan_image_collection(
            collection, ranks=preset.ranks, scales=preset.scales, max_iter=preset.max_iter
        )
        records_b, fits_b = trial.scan_single_image(
            collection, ranks=preset.ranks, scales=preset.scales, max_iter=preset.max_iter
        )
        records = records_a + records_b
        columns = (
            "flattening", "rank", "robust_scale", "ratio", "control_ratio", "excess",
            "auc_rows", "iterations", "seconds",
        )  # fmt: skip
        _print_table(records, columns)
        print(f"   ({time.perf_counter() - started:.1f} s for both flattenings, JIT included)")
        for record in records:
            record["peak_snr"] = collection.peak_signal_to_noise
        all_records += records
        best_a = max(records_a, key=lambda r: r["excess"])
        best_b = max(records_b, key=lambda r: r["excess"])
        path = figures.save_image_figure(
            collection,
            fits_a[(best_a["rank"], best_a["robust_scale"])],
            fits_b[(best_b["rank"], best_b["robust_scale"])],
            out,
            name=f"image_scores_x{background:g}.png",
        )
        print(f"   wrote {path}")
    _write_csv(out / "image_grid.csv", all_records)
    return all_records


def main(argv: Sequence[str] | None = None) -> int:
    """Run the trial and return an exit code."""
    arguments = _parser().parse_args(argv)
    from . import trial

    version = trial.require_rhmf() and trial.installed_version()
    preset = QUICK if arguments.quick else FULL
    out = pathlib.Path(arguments.out)
    out.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    print(
        f"# RHMF trial ({'quick' if arguments.quick else 'full'}): robusta-hmf {version}, "
        f"commit {trial.COMMIT[:12]}; ranks {preset.ranks}, scales {preset.scales}; "
        f"output in {out}"
    )
    print("# no default is promoted: the grid is scanned and reported, and 'best' knows the answer")
    if arguments.only in (None, "spectra"):
        run_spectra(preset, out, arguments.seed)
    if arguments.only in (None, "image"):
        run_image(preset, out, arguments.seed)
    print(f"\n# done in {time.perf_counter() - started:.1f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
