# Security policy

## Supported versions

Security fixes are made for the `1.0.0b1` line onward. Legacy ampere
(`ampere.legacy`, and the pre-v2 top-level names `ampere.data`,
`ampere.models` and `ampere.infer`) is frozen and receives no security
fixes; see [`docs/source/legacy.rst`](docs/source/legacy.rst).

## Reporting a vulnerability

Please report privately, not in a public issue or discussion.

1. Use GitHub's private vulnerability reporting:
   <https://github.com/ICSM/ampere/security/advisories/new>.
2. If that is unavailable to you, email peter.scicluna@eso.org.

Please include the ampere version
(`python -c "import ampere; print(ampere.__version__)"`), what you did, and
what happened.

## What to expect

- An acknowledgement within fourteen days of your report.
- A fix on a best-effort basis, with no stated deadline. Ampere is a
  scientific library maintained by one person and has no network surface of
  its own.
- No embargo, CVE assignment or coordinated-disclosure process. If a fix is
  made, it will be described in the release notes.

## What ampere does with files you give it

Ampere runs your Python and reads your files with your privileges. These
behaviours are documented and by design; take the same care as you would with
any pickle.

- **The SBI artefact cache unpickles.** `ArtefactStore` (in
  `ampere.results.artefacts`) loads cached posteriors with `pickle.load` from
  digest-named `.pkl` files. Do not point it at a cache directory from a
  source you do not trust.
- **`ProcessExecutor` pickles whole problems.** It sends a `FittingProblem`,
  user models and simulators included, into worker processes. Do not run a
  problem, or load a pickled one, that you did not build or do not trust.
- **Runs and training sets are netCDF.** They are read through xarray and the
  HDF5/netCDF-C libraries. Do not open run or training-set files from an
  untrusted source; a malformed file is a risk in those libraries as much as
  in ampere.
- **Models and simulators are user code.** Anything you register or pass as a
  model, simulator or likelihood is ordinary Python and runs with your
  privileges.

A report that only describes one of these documented behaviours is not a
vulnerability report. A way to make ampere do something beyond them is.
