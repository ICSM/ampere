# Contributing to ampere

Thank you for considering a contribution. Ampere is a Bayesian fitting
environment for heterogeneous astronomical data; this page is the short
version of how to work on it. The longer onboarding note is
[`docs/development.md`](docs/development.md).

Please read the [Code of Conduct](CODE_OF_CONDUCT.md) first. To report a
security problem, follow [`SECURITY.md`](SECURITY.md) rather than opening a
public issue.

## Setting up

The supported route is [pixi](https://pixi.sh). From a clean clone:

```bash
pixi install -e dev
```

There are three real environments: `dev` (the daily-use one, with neither
torch nor jax), `torch` and `jax` (each `dev` plus its backend). Run a task
in a given environment with `pixi run -e <env> <task>`; plain
`pixi run <task>` means `-e dev`. If you prefer pip,
`pip install -e ".[dev]"` works on Python 3.12 or later.

## Before you open a pull request

```bash
pixi run test-fast      # core, results, conformance, backends, inference; about 3 minutes
pixi run lint           # ruff check
pixi run format-check   # ruff format --check
pixi run typecheck      # pyrefly, scoped to the v2 namespaces
pixi run docs           # only if you touched docs or docstrings
```

If your change touches a backend, also run `pixi run -e torch test-all` or
`pixi run -e jax test-all` for that backend. The full gate (`test-all` in
every environment) is CI's job, not yours. Please say in the pull request
which tasks you ran and what they reported.

## Branches and pull requests

- One topic per branch and per pull request. Name the branch for what it
  does, for example `fix-closure-phase-wrap` or `w6.17-community-health`
  for a numbered work item.
- Open the pull request against the repository's default branch. Do not push to it directly.
- The maintainer, Peter Scicluna, reviews and merges every change. Expect
  questions; a review is a conversation, not a verdict.
- Stay in scope. If you notice another problem while working, record it in
  the pull request description or open an issue rather than fixing it in the
  same change.

## What is frozen, and why

- **Legacy is frozen.** Nothing under `ampere/legacy/` changes (the old
  top-level names `ampere.data`, `ampere.models` and `ampere.infer` are
  aliases for it). It produced published science and stays as it was. The
  policy is in [`docs/source/legacy.rst`](docs/source/legacy.rst).
- **The core contracts are frozen** at `spec-v1.0`. Any change to a contract
  in [`DEVELOPMENT_PLAN.md`](DEVELOPMENT_PLAN.md) section 4 needs a
  decision-log entry in that file in the same pull request, and must keep the
  conformance suite green (or update it, with justification). This is ground
  rule 9 in [`AGENTS.md`](AGENTS.md).

[`DEVELOPMENT_PLAN.md`](DEVELOPMENT_PLAN.md) is the source of truth for the
architecture and its known traps (section 7); read it before any non-trivial
change.

## Conventions

- British English in all prose: documentation, docstrings, comments, commit
  messages.
- Code in the v2 namespaces (`ampere/core`, `backends`, `inference`,
  `results`) is typed from the first line.
- No run outputs, logs, checkpoints or binary artefacts in git.
- Commit messages explain why. When an agent assisted, end the message with a
  `Co-Authored-By:` trailer naming it, for example
  `Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>`.

## Where to talk

- **Issues** for bugs and feature requests; the forms ask for what is needed
  to reproduce a problem (the version, the backend, the smallest
  `FittingProblem` that shows it, the traceback).
- **Discussions** for questions, ideas and "how do I fit...?".

## How v2 was built

Much of the v2 code was written by orchestrated coding agents working under
the working agreement in [`AGENTS.md`](AGENTS.md) and
[`WORK_ITEMS.md`](WORK_ITEMS.md): one work item per branch, explicit
acceptance criteria, and a report with evidence. Every change was reviewed
and merged by the maintainer. You are welcome to work the same way, with or
without an agent of your own, or by hand alone; the checks above and the
review are the same either way.
