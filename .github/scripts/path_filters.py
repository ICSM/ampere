#!/usr/bin/env python3
"""W4.10: the single source of truth for CI's path-gating rules.

``.github/workflows/ci.yml``'s ``changes`` job calls this script to decide
which of the path-gated jobs actually do their work on a given pull
request; this file — not a duplicated set of glob patterns pasted into the
workflow's ``if:`` conditions — is the one place the rule table lives, so a
change to the rules is a one-file diff and a local run reproduces exactly
what CI decided. The path-filter table below is also reproduced in
``ci.yml``'s header comment and in ``docs/development.md``'s CI section
(three copies of the same table would drift; this script is the one that
is actually executed, the other two are read by humans).

Usage
-----
In CI (``ci.yml``'s ``changes`` job)::

    python .github/scripts/path_filters.py \\
        --event-name "$EVENT_NAME" --base "$BASE_SHA" --head "$HEAD_SHA" \\
        --github-output "$GITHUB_OUTPUT"

Locally, against a synthetic diff (this is the acceptance evidence for
W4.10 — "the report shows the job list for three synthetic diffs")::

    printf 'docs/source/index.rst\\nREADME.md\\n' | \\
        python .github/scripts/path_filters.py --stdin

    python .github/scripts/path_filters.py --paths ampere/backends/torch/problem.py

Both forms print the same two things: which bucket each changed path fell
into, and which named CI jobs/checks would run as a result. Nothing here
calls the GitHub API or needs network access — it is pure path
classification, which is what makes it usable outside CI.

The path-filter table
----------------------
Buckets, in priority order (a path is classified by the first rule that
matches it):

  BLACKJAX -- ``ampere/inference/_blackjax.py``,
              ``tests/inference/test_blackjax.py`` (W5.14). Checked BEFORE
              CORE, for the reason the NESTED entry below gives, and gated
              to the **jax** leg rather than a leg of its own: blackjax is
              a jax library (memo §2.2: there is no torch counterpart to
              borrow) and it installs into the existing ``jax``
              environment, so there is no new environment to run and
              nothing for a fourth leg to do. Runs dev -- whose
              ``test-all`` still imports the namespace and parses its
              import graph -- and the jax leg.
  NESTED   -- ``ampere/inference/_nested.py``,
              ``tests/inference/test_nested.py`` (W5.14). Checked BEFORE
              CORE, which would otherwise claim them: these two files are
              the whole of the nautilus/ultranest drivers, they are reached
              only through ``ampere/inference/__init__.py`` (itself CORE, so
              a change there still runs everything), and nothing under
              ``ampere/backends`` can be affected by them. Runs dev -- whose
              ``test-all`` still imports the namespace and parses its import
              graph -- and the nested leg.
  CORE     -- ``pyproject.toml``, ``pixi.lock``, ``.github/**``,
              ``ampere/core/**``, ``ampere/results/**``,
              ``ampere/inference/**``, and the shared/backend-neutral test
              suites (``tests/core``, ``tests/results``, ``tests/inference``,
              ``tests/conformance``, ``tests/m2``, ``tests/interferometry``,
              ``tests/astrometry`` (both W5.26), ``tests/benchmarks``,
              ``tests/scaling``, ``tests/gpu``, ``tests/characterisation``).
              Runs EVERYTHING: dev, torch, jax, sbi, nested and the docs
              build.
              (A registry-leak regression in shared code is exactly what
              ``test-all``'s single-process suites exist to catch across
              every backend -- see ``pyproject.toml``'s ``test-all`` task
              comment -- so a shared-code change cannot skip any leg.)
  TORCH    -- ``ampere/backends/torch/**``, ``tests/backends/*torch*``.
              Runs the torch leg, and the sbi leg (sbi's environment
              installs torch and typechecks ``ampere/backends/torch`` too
              -- see ``backend-typecheck``'s job comment in ci.yml).
  JAX      -- ``ampere/backends/jax/**``, ``tests/backends/*jax*``.
              Runs the jax leg only.
  DEV_ONLY -- ``ampere/backends/reference/**`` and anything else under
              ``tests/backends/`` not caught above (the reference backend's
              own tests, ``tests/backends/__init__.py``). The reference
              backend is exercised inside the *dev* environment's
              ``test-all`` (tests/conformance's reference fixture,
              tests/backends/test_reference*.py); it shares no code path
              with the torch/jax legs' own fixtures, so it need not run
              them. Runs dev only.
  EX_SBI   -- ``examples/sbi/**``. These need the ``sbi`` extra (several
              are skip-marked without it) and are exercised by
              ``tests/examples`` inside ``test-all`` on every environment,
              but only actually run on the sbi leg. Runs dev and sbi.
  EX_IFM   -- ``examples/interferometry/**`` -- W4.4's study, not landed at
              W4.10. MARKED SLOT: until W4.4 defines which model file
              belongs to which backend, this bucket conservatively runs
              dev plus every backend leg (torch, jax, sbi); W4.4 should
              narrow it to the backend(s) each example file actually
              exercises rather than widen it further.
  EXAMPLES -- ``examples/**`` (anything not caught above),
              ``tests/examples/**``. Runs dev only.
  DOCS     -- ``docs/**``, any ``*.md`` file anywhere. Runs dev and the
              docs build.
  OTHER    -- anything not matched above (legacy ``ampere/{data,models,
              infer,utils}``, root files with no extension, etc.) -- an
              unrecognised path is exactly the case this table cannot
              afford to get wrong, so it is treated like CORE and runs
              everything.

``run_dev``    = CORE | DEV_ONLY | EXAMPLES | EX_SBI | EX_IFM | DOCS | NESTED
                 | BLACKJAX | OTHER
``run_torch``  = CORE | TORCH | EX_IFM | OTHER
``run_jax``    = CORE | JAX | EX_IFM | BLACKJAX | OTHER
``run_sbi``    = CORE | TORCH | EX_SBI | EX_IFM | OTHER
``run_nested`` = CORE | NESTED | OTHER
``run_docs``   = CORE | DOCS | OTHER

Non-pull_request events (push to master, workflow_dispatch, schedule) and
any failure to compute a diff (missing base/head, git error, or a genuinely
empty diff) fall back to "everything", by design -- see the module
docstring's "Safety" note below repeated at the call site. Path gating only
narrows *pull_request* runs; a push to master (post-merge confidence) and
the weekly/manual runs always exercise every leg.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import PurePosixPath

# --- job/check names, exactly as they appear in ci.yml's `name:` fields ---
ALWAYS_JOBS = [
    "lint + format-check",
    "actionlint",
]
# W5.26: `suites`/`backend-suites` each gained a `group` matrix dimension
# (`core`, `backends`, `studies` -- pyproject.toml's `test-group-*` tasks),
# so the one "new-namespace suites (<environment>)" check per environment
# below is now three, one per group.
_SUITE_GROUPS = ("core", "backends", "studies")


def _suite_jobs(environment: str) -> list[str]:
    return [f"new-namespace suites ({environment}, {group})" for group in _SUITE_GROUPS]


DEV_JOBS = [
    "typecheck (pyrefly, new namespaces)",
    "test (py311)",
    "test (py312)",
    "test (py313)",
    *_suite_jobs("dev"),
    "minimal install (no extras)",
]
TORCH_JOBS = ["typecheck (pyrefly, torch)", *_suite_jobs("torch")]
JAX_JOBS = ["typecheck (pyrefly, jax)", *_suite_jobs("jax")]
SBI_JOBS = ["typecheck (pyrefly, sbi)", *_suite_jobs("sbi")]
# W5.14: the nautilus/ultranest environment. A small, numpy-only
# environment (no torch, no jax), which is what makes a leg of its own
# affordable rather than folding two nested samplers into `dev`.
NESTED_JOBS = ["typecheck (pyrefly, nested)", *_suite_jobs("nested")]
DOCS_JOBS = ["docs build"]
# Note: `characterisation suite (with sbi extra)` (sbi-characterisation) is
# not listed here -- it is never path-gated (see ci.yml's job comment); it
# keeps its own schedule/workflow_dispatch-only `if:`, independent of this
# script.

CORE = "CORE"
NESTED = "NESTED"
BLACKJAX = "BLACKJAX"
TORCH = "TORCH"
JAX = "JAX"
DEV_ONLY = "DEV_ONLY"
EX_SBI = "EX_SBI"
EX_IFM = "EX_IFM"
EXAMPLES = "EXAMPLES"
DOCS = "DOCS"
OTHER = "OTHER"

_CORE_EXACT = {"pyproject.toml", "pixi.lock"}
_NESTED_EXACT = {"ampere/inference/_nested.py", "tests/inference/test_nested.py"}
_BLACKJAX_EXACT = {"ampere/inference/_blackjax.py", "tests/inference/test_blackjax.py"}
_CORE_PREFIXES = (
    ".github/",
    "ampere/core/",
    "ampere/results/",
    "ampere/inference/",
    "tests/core/",
    "tests/results/",
    "tests/inference/",
    "tests/conformance/",
    "tests/m2/",
    # W5.26: joined `test-all` (and the `test-group-studies` cell of the new
    # CI matrix); both are backend-neutral shared suites in the same sense
    # as `tests/m2` above (their per-backend rows skip themselves via
    # `pytest.importorskip`/`needs_torch`/`needs_jax`, exactly like the rest
    # of this bucket), so they run in every environment's gate.
    "tests/interferometry/",
    "tests/astrometry/",
    "tests/benchmarks/",
    "tests/scaling/",
    "tests/gpu/",
    "tests/characterisation/",
)


def classify(path: str) -> str:
    """Classify one changed path into a bucket. See the module docstring's table."""
    p = PurePosixPath(path)
    posix = p.as_posix()

    # Before CORE, deliberately: see the module docstring's NESTED entry.
    if posix in _NESTED_EXACT:
        return NESTED
    # Same rule, gated to the jax leg: see the BLACKJAX entry.
    if posix in _BLACKJAX_EXACT:
        return BLACKJAX

    if posix in _CORE_EXACT or posix.startswith(_CORE_PREFIXES):
        return CORE

    if posix.startswith("ampere/backends/torch/"):
        return TORCH
    if posix.startswith("tests/backends/") and "torch" in p.name:
        return TORCH

    if posix.startswith("ampere/backends/jax/"):
        return JAX
    if posix.startswith("tests/backends/") and "jax" in p.name:
        return JAX

    if posix.startswith("ampere/backends/reference/") or posix.startswith("tests/backends/"):
        return DEV_ONLY

    # W4.4 marked slot -- see EX_IFM in the module docstring.
    if posix.startswith("examples/interferometry/"):
        return EX_IFM

    if posix.startswith("examples/sbi/"):
        return EX_SBI

    if posix.startswith("examples/") or posix.startswith("tests/examples/"):
        return EXAMPLES

    if posix.startswith("docs/") or posix.endswith(".md"):
        return DOCS

    return OTHER


def compute_flags(buckets: set[str]) -> dict[str, bool]:
    return {
        "run_dev": bool(
            buckets & {CORE, DEV_ONLY, EXAMPLES, EX_SBI, EX_IFM, DOCS, NESTED, BLACKJAX, OTHER}
        ),
        "run_torch": bool(buckets & {CORE, TORCH, EX_IFM, OTHER}),
        "run_jax": bool(buckets & {CORE, JAX, EX_IFM, BLACKJAX, OTHER}),
        "run_sbi": bool(buckets & {CORE, TORCH, EX_SBI, EX_IFM, OTHER}),
        "run_nested": bool(buckets & {CORE, NESTED, OTHER}),
        "run_docs": bool(buckets & {CORE, DOCS, OTHER}),
    }


def job_list(flags: dict[str, bool]) -> tuple[list[str], list[str]]:
    """Return (jobs that run, jobs skipped) given the computed flags."""
    running = list(ALWAYS_JOBS)
    skipped: list[str] = []
    buckets_and_jobs = (
        ("run_dev", DEV_JOBS),
        ("run_torch", TORCH_JOBS),
        ("run_jax", JAX_JOBS),
        ("run_sbi", SBI_JOBS),
        ("run_nested", NESTED_JOBS),
        ("run_docs", DOCS_JOBS),
    )
    for key, jobs in buckets_and_jobs:
        if flags[key]:
            running.extend(jobs)
        else:
            skipped.extend(jobs)
    return running, skipped


def git_diff_paths(base: str, head: str) -> list[str]:
    result = subprocess.run(
        ["git", "diff", "--name-only", base, head],
        capture_output=True,
        text=True,
        check=True,
    )
    return [line for line in result.stdout.splitlines() if line.strip()]


def run_all_flags() -> dict[str, bool]:
    return {
        "run_dev": True,
        "run_torch": True,
        "run_jax": True,
        "run_sbi": True,
        "run_nested": True,
        "run_docs": True,
    }


def report(paths: list[str] | None, flags: dict[str, bool], reason: str) -> None:
    print(f"# path-filter decision ({reason})", file=sys.stderr)
    if paths is None:
        print("  (no diff computed -- running everything)", file=sys.stderr)
    else:
        if not paths:
            print("  (empty diff)", file=sys.stderr)
        for path in paths:
            print(f"  {path}  -> {classify(path)}", file=sys.stderr)
    print(
        "  flags: " + " ".join(f"{k}={'true' if v else 'false'}" for k, v in flags.items()),
        file=sys.stderr,
    )
    running, skipped = job_list(flags)
    print("  jobs that run:", file=sys.stderr)
    for name in running:
        print(f"    - {name}", file=sys.stderr)
    print("  jobs skipped (step-level; the job itself still reports success):", file=sys.stderr)
    for name in skipped:
        print(f"    - {name}", file=sys.stderr)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--event-name",
        default="",
        help="github.event_name, e.g. pull_request, push, schedule, workflow_dispatch",
    )
    parser.add_argument("--base", default="", help="base SHA to diff from (pull_request only)")
    parser.add_argument("--head", default="", help="head SHA to diff to (pull_request only)")
    parser.add_argument(
        "--stdin",
        action="store_true",
        help="read a newline-separated list of changed paths from stdin (synthetic-diff testing)",
    )
    parser.add_argument(
        "--paths",
        nargs="*",
        default=None,
        help="explicit changed paths (synthetic-diff testing)",
    )
    parser.add_argument(
        "--github-output",
        default="",
        help="path to $GITHUB_OUTPUT; if given, the flags are appended there",
    )
    args = parser.parse_args()

    paths: list[str] | None
    reason: str

    if args.stdin:
        paths = [line.strip() for line in sys.stdin if line.strip()]
        reason = "synthetic diff via --stdin"
    elif args.paths is not None:
        paths = args.paths
        reason = "synthetic diff via --paths"
    elif args.event_name == "pull_request" and args.base and args.head:
        try:
            paths = git_diff_paths(args.base, args.head)
            reason = f"pull_request diff {args.base[:12]}..{args.head[:12]}"
        except subprocess.CalledProcessError as exc:
            msg = f"warning: git diff failed ({exc}); falling back to 'everything'"
            print(msg, file=sys.stderr)
            paths = None
            reason = "git diff failed -- fallback"
    else:
        event = args.event_name or "(unset)"
        paths = None
        reason = f"event_name={event} -- not a pull_request diff, running everything"

    if paths is None:
        flags = run_all_flags()
    elif not paths:
        # An empty diff on a pull_request is unusual (an empty commit) and not
        # worth a special case that skips everything -- run everything rather
        # than nothing.
        flags = run_all_flags()
    else:
        buckets = {classify(p) for p in paths}
        flags = compute_flags(buckets)

    report(paths, flags, reason)

    if args.github_output:
        with open(args.github_output, "a", encoding="utf-8") as fh:
            for key, value in flags.items():
                fh.write(f"{key}={'true' if value else 'false'}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
