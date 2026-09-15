#!/usr/bin/env python
"""An external simulator as a first-class ampere model, run under a process pool.

``DEVELOPMENT_PLAN.md`` §2's *Batched simulation* row, note (2): a compiled
Fortran/C/C++/Rust routine behind a Python call is **the** canonical SBI
simulator in this field, not a fallback. This example is that case end to end —
a black-box :class:`~ampere.core.transform.Model` on the reference backend whose
``evaluate`` marshals arguments, runs a program in a working directory of its
own, reads the answer back off disk, and turns a non-zero exit status into a
flagged failure rather than an exception — driven by
:meth:`~ampere.core.dataset.FittingProblem.simulate_many` over a
:class:`~ampere.core.simulate.ProcessExecutor`.

Run it::

    python examples/sbi/external_simulator.py            # 24 draws, 2 workers
    python examples/sbi/external_simulator.py 200 4      # 200 draws, 4 workers

The four things it is here to demonstrate
-----------------------------------------
**Argument marshalling.** ``ExternalPowerlaw.evaluate`` writes the evaluation
grid to a file, builds a command line, and parses the output file. Nothing about
the simulator knows what ampere is, and nothing about ampere knows what the
simulator is: the contract between them is a command line and two files, which
is exactly the contract a wrapped compiled code offers.

**A working directory per worker.** The scratch directory is created lazily,
keyed on the process id, and is therefore *not* part of the model's pickled
state: the model travels to a worker, the directory is made there, and two
workers never write to the same files. Simulators that insist on fixed
filenames in the current directory — and there are many — need exactly this.

**stdout and stderr into the failure detail.** When the program exits non-zero,
the wrapper raises :class:`SimulatorFailed` carrying the tail of stderr. That
class is declared to the problem through ``simulator_failures=``, which is what
turns it into ``inference.md`` §11's flagged
:class:`~ampere.core.dataset.Failure` — reason ``model_failed``, exception type
``SimulatorFailed``, and the simulator's own diagnostic as the message. The
budget carries on; the count is on the problem afterwards.

**A timeout that actually stops the work.** ``ProcessExecutor(timeout=...)``
kills the worker whose draw overruns and flags that draw ``execution_failed``,
which is a different count from a simulator that reported its own failure. Pass
``--slow`` to see it: a fraction of the draws are told to sleep past the
deadline.

What this example is *not*
--------------------------
It is not an SBI fit. It stops at the ``(θ, x)`` pairs, because that is where
W3.1 slice 1 stops; W3.2's ``SBIEngine`` consumes exactly what this produces,
with the same ``executor=`` and ``chunk_size=`` arguments passed straight
through.
"""

from __future__ import annotations

import atexit
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

import astropy.units as u
import numpy as np
import scipy.stats as st

from ampere.core import (
    Dataset,
    FailureReason,
    FittingProblem,
    Model,
    Parameter,
    ProcessExecutor,
    SimulationBatch,
    Spectrum,
)

#: The "compiled" program. A separate file, run as a subprocess, exactly as a
#: real binary would be — this example never imports it.
SIMULATOR = Path(__file__).resolve().parent / "toy_powerlaw.py"

#: How much of the simulator's stderr is kept on a failure. A crashing code can
#: produce megabytes of diagnostics and the failure record is kept per draw, so
#: the tail is what goes in: it is where the reason is.
STDERR_TAIL = 400

GRID = np.geomspace(1.0, 10.0, 12)

#: One scratch directory per worker process, made on first use *in that process*.
#: Deliberately module state rather than an attribute: the model is pickled to
#: the worker, and a path made in the parent would have every worker writing to
#: one directory.
_WORKING_DIRECTORIES: dict[int, Path] = {}


class SimulatorFailed(RuntimeError):
    """The external program declined to produce a result.

    Declared to :class:`~ampere.core.dataset.FittingProblem` through
    ``simulator_failures=``, which is the whole of what a user has to do to turn
    their simulator's own error class into a flagged failure
    (``inference.md`` §11: the catch set is deliberately narrow, and it is the
    user who widens it).
    """


@atexit.register
def _clean_working_directories() -> None:
    """Remove this process's scratch directory when it exits.

    Registered per process, so each worker tidies up after itself. A worker
    killed by :class:`~ampere.core.simulate.ProcessExecutor`'s timeout never
    runs this — a killed process runs no ``atexit`` handlers — which is one more
    honest cost of a simulator that has to be stopped rather than asked to stop.
    """
    for path in _WORKING_DIRECTORIES.values():
        shutil.rmtree(path, ignore_errors=True)


def working_directory() -> Path:
    """This process's scratch directory, made once, on first use."""
    pid = os.getpid()
    existing = _WORKING_DIRECTORIES.get(pid)
    if existing is None:
        existing = Path(tempfile.mkdtemp(prefix=f"ampere-external-{pid}-"))
        _WORKING_DIRECTORIES[pid] = existing
    return existing


class ExternalPowerlaw(Model):
    """A black-box model whose ``evaluate`` is a subprocess call.

    Conservative capability flags, and they are honest ones: a subprocess is not
    differentiable, and it does not take a stack of θ — this simulator is a
    one-parameter-set-per-run program, which is the common case. A code that
    *does* take a table would implement
    :meth:`~ampere.core.transform.Model.evaluate_batch` and set
    ``BATCHABLE = True``.

    Parameters
    ----------
    grid
        Wavelengths to evaluate on, in microns.
    program
        The executable to run. Defaults to :data:`SIMULATOR`.
    sleep
        Seconds to make the simulator sleep for draws whose ``index`` is below
        ``-1.9`` — how the example produces a draw that overruns a timeout.
    """

    def __init__(
        self, grid: np.ndarray = GRID, *, program: Path = SIMULATOR, sleep: float = 0.0
    ) -> None:
        self.register_buffer("grid", np.asarray(grid, dtype=float), unit=u.micron)
        self.register_parameter(Parameter("index", st.uniform(-2.5, 2.5)))
        self.register_parameter(Parameter("norm", st.uniform(0.5, 9.5)))
        self.program = Path(program)
        self.sleep = float(sleep)

    def evaluate(self, **values: Any) -> Spectrum:
        context = self.context(values)
        grid = np.asarray(context["grid"], dtype=float)
        index = float(context["index"])
        norm = float(context["norm"])

        # (1) Marshal. A file for the grid, a command line for the scalars, and
        #     a named output file: the interface a compiled code actually has.
        scratch = working_directory()
        grid_path = scratch / "grid.txt"
        out_path = scratch / "flux.txt"
        grid_path.write_text("\n".join(repr(float(x)) for x in grid) + "\n")
        command = [
            sys.executable,
            str(self.program),
            "--index",
            repr(index),
            "--norm",
            repr(norm),
            "--grid",
            str(grid_path),
            "--out",
            str(out_path),
        ]
        if self.sleep and index < -1.9:
            command += ["--sleep", repr(self.sleep)]

        # (2) Run it, in its own directory, capturing both streams.
        completed = subprocess.run(
            command,
            cwd=scratch,
            capture_output=True,
            text=True,
            check=False,
        )

        # (3) A non-zero exit is a *failure*, not a bug: raise the declared
        #     class, carrying the simulator's own diagnostic, and let
        #     `simulate` flag it. Nothing here catches anything: flagging is the
        #     problem's job, and doing it here would hide the draw's theta.
        if completed.returncode != 0 or not out_path.exists():
            raise SimulatorFailed(
                f"{self.program.name} exited {completed.returncode} at "
                f"index={index:.6g}, norm={norm:.6g}; "
                f"stderr: {completed.stderr.strip()[-STDERR_TAIL:] or '(empty)'}; "
                f"stdout: {completed.stdout.strip()[-STDERR_TAIL:] or '(empty)'}"
            )

        # (4) Read the answer back and put it in a container.
        flux = np.loadtxt(out_path, ndmin=1)
        out_path.unlink()
        return Spectrum(grid * u.micron, flux * u.Jy)


def build_problem(*, seed: int = 20260909, sleep: float = 0.0) -> FittingProblem:
    """One dataset observed by the external simulator, with its failure class declared."""
    truth = 2.0 * GRID**-1.0
    observed = Spectrum(
        GRID * u.micron,
        truth * u.Jy,
        uncertainty=np.full(GRID.size, 0.05) * u.Jy,
        mask=None,
    )
    return FittingProblem(
        ExternalPowerlaw(GRID, sleep=sleep),
        [Dataset(observed, label="sed")],
        seed=seed,
        # The one line that makes this simulator's crashes flagged rather than
        # fatal. Without it, the first bad draw ends the budget.
        simulator_failures=(SimulatorFailed,),
    )


def run(
    count: int = 24,
    workers: int = 2,
    *,
    timeout: float | None = None,
    sleep: float = 0.0,
    chunk_size: int | None = None,
) -> tuple[FittingProblem, SimulationBatch]:
    """Simulate *count* draws through a pool of *workers* subprocesses."""
    problem = build_problem(sleep=sleep)
    executor = ProcessExecutor(workers, timeout=timeout)
    batch = problem.simulate_many(count, observe=True, executor=executor, chunk_size=chunk_size)
    return problem, batch


def report(problem: FittingProblem, batch: SimulationBatch) -> str:
    """What a driver should print after a budget: usable pairs, and why the rest are not."""
    lines = [
        (
            f"{len(batch)} draw(s) simulated, {len(batch.usable)} usable "
            f"({int(batch.failed.sum())} failed)"
        ),
    ]
    for reason, number in sorted(problem.failure_counts.items()):
        lines.append(f"  {number:>4} x {reason}")
    first = next((failure for failure in batch.failures if failure is not None), None)
    if first is not None:
        lines.append(f"  first failure: {first}")
    if batch.observations is not None:
        lines.append(f"  observations stacked as {batch.observations['sed'].shape}")
    lines.append(f"  theta stacked as {batch.theta.shape}")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    arguments = list(sys.argv[1:] if argv is None else argv)
    slow = "--slow" in arguments
    positional = [item for item in arguments if not item.startswith("-")]
    count = int(positional[0]) if positional else 24
    workers = int(positional[1]) if len(positional) > 1 else 2

    problem, batch = run(
        count,
        workers,
        timeout=1.0 if slow else None,
        sleep=5.0 if slow else 0.0,
    )
    print(report(problem, batch))
    if slow:
        expired = problem.failure_counts.get(FailureReason.EXECUTION_FAILED, 0)
        print(f"  {expired} draw(s) were killed by the per-simulation timeout")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
