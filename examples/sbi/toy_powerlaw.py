#!/usr/bin/env python
"""A stand-in for a compiled simulator, with a compiled simulator's manners.

This script is **not** ampere code and imports nothing from ampere. It exists to
be run as a subprocess by ``external_simulator.py``, and it behaves the way the
Fortran, C++ and Rust radiative-transfer codes ampere is asked to wrap actually
behave:

* parameters arrive as **command-line arguments and an input file**, not as
  Python objects;
* output goes to a **file in the working directory**, named on the command line;
* progress chatter goes to **stdout** and diagnostics to **stderr**;
* a parameter combination it cannot handle is a **non-zero exit status** with a
  message on stderr, not an exception.

The physics is a power law, ``f(x) = norm * x ** index``, so that the wrapper
around it can be checked against a closed form.

Usage
-----
    python toy_powerlaw.py --index -1.0 --norm 2.0 --grid grid.txt --out flux.txt

Exit codes
----------
0
    Success; ``--out`` holds one flux per grid point.
2
    "Failed to converge": raised for ``norm`` above 9.5, which is how this toy
    injects the crash rate a real budget has.
"""

from __future__ import annotations

import argparse
import sys
import time


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="toy power-law simulator")
    parser.add_argument("--index", type=float, required=True)
    parser.add_argument("--norm", type=float, required=True)
    parser.add_argument("--grid", required=True, help="input file: one coordinate per line")
    parser.add_argument("--out", required=True, help="output file: one flux per line")
    parser.add_argument("--sleep", type=float, default=0.0, help="pretend to be slow")
    args = parser.parse_args(argv)

    print(f"toy_powerlaw: index={args.index:g} norm={args.norm:g}", flush=True)

    if args.norm > 9.5:
        print(
            f"toy_powerlaw: FATAL: failed to converge at norm={args.norm:.6g} "
            f"(the opacity table does not extend this far)",
            file=sys.stderr,
            flush=True,
        )
        return 2

    if args.sleep > 0.0:
        time.sleep(args.sleep)

    with open(args.grid) as handle:
        grid = [float(line) for line in handle if line.strip()]
    with open(args.out, "w") as handle:
        for coordinate in grid:
            handle.write(f"{args.norm * coordinate**args.index!r}\n")

    print(f"toy_powerlaw: wrote {len(grid)} point(s) to {args.out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
