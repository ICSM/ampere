"""W4.5's demonstration: fringing modelled, not merely absorbed.

M2's two fringing scenarios deviate from the model by a **sinusoid** of known
period — 0.0028 µm, an optical path difference, the same shape a real etalon
leaves in a real spectrum. The study fits them with a stationary Matérn-3/2 and
the flexible likelihood stays calibrated: that is M2's result and it is not in
question here. But a stationary kernel is the wrong *shape* for a ripple. It
has one length scale to spend, and it must spend it either short enough to
follow each oscillation — in which case it learns nothing about the ripple's
coherence — or long enough to see the coherence, in which case it cannot follow
the oscillation. The horizon note's §1 follow-up named the consequence and the
remedy in one sentence: "fringing in M2 is currently absorbed by a stationary
Matérn-3/2 (it stays calibrated), not modelled as a periodic component", and
the component for it is a damped harmonic oscillator.

This module runs that comparison. Same data, same physical model, same sampler
budget; two flexible likelihoods, one with `Matern32` and one with
`Matern32 + SHO`. It reports, per physical parameter, the bias in posterior
widths and whether the 68 % interval covers the truth — beside the stationary
fit, and beside the standard likelihood, so the reader can see all three.

**This is informational.** M2's pinned assertions are untouched and this module
adds none: the sum kernel is a *better-shaped* model for this deviation, not a
better-calibrated one — the flexible likelihood was already calibrated, which
was the whole point of M2 — and the number that should move is the GP's own
account of the residual, not the coverage. Reporting a pass/fail here would be
inventing a threshold for a comparison nobody has a prediction for.

Run it as::

    python -m examples.m2_misspecification.fringing
    python -m examples.m2_misspecification.fringing --scenario mild --size 400

Nothing here writes a file (ground rule 7): the table goes to stdout.
"""

from __future__ import annotations

import argparse
import dataclasses
from typing import Any

import astropy.units as u
import scipy.stats as st

from ampere.core import Dataset, FittingProblem, GaussianFamily, Kernel, Likelihood, Sum

from .generators import SyntheticSpectrum, generate
from .study import (
    DATASET_LABEL,
    GP_AMPLITUDE_SCALE,
    GP_LENGTH_SCALE,
    PHYSICAL_NAMES,
    SEED,
    TEST_EMCEE,
    EmceeBudget,
    Summary,
    _backend_module,
    build_likelihood,
    model_for,
    run,
    summarise,
)

__all__ = [
    "FRINGE_PERIOD",
    "FRINGE_QUALITY",
    "Comparison",
    "build_fringing_kernel",
    "build_fringing_problem",
    "compare",
    "main",
]

#: The fringe period the two smooth scenarios were generated with, micron.
#: Given to the SHO term as a *prior scale*, not as a fixed value: a real
#: instrument's fringe period is known to a factor rather than exactly, and a
#: demonstration that handed the model the generating number would be
#: demonstrating nothing.
FRINGE_PERIOD = 0.0028

#: The SHO's quality factor, held fixed. ``Q`` is how many periods the ripple
#: stays coherent for, to within a factor of pi; 15 is "coherent across the
#: band but not forever", which is what an etalon fringe is. Fixed rather than
#: fitted because the demonstration is about the *shape* of the component, and
#: a third free hyperparameter would make the comparison against the two-parameter
#: stationary kernel a comparison of flexibility as well as of shape.
FRINGE_QUALITY = 15.0


def build_fringing_kernel(backend: str = "reference") -> Kernel:
    """``Matern32 + SHO``: a broad stationary term plus a damped ripple.

    The sum is what W4.5 added and the reason it is the right answer here: a
    sum of quasiseparable terms is quasiseparable, so this kernel reaches the
    same O(N) solver the stationary one does, at rank 4 instead of rank 2. The
    demonstration would not be affordable at M2's 20 000 points otherwise.

    Both terms are this backend's own classes, reached by name rather than by a
    branch — every backend's namespace spells them identically (W2.12) — so the
    comparison runs on the reference, torch and jax paths without a
    backend-specific line.
    """
    module = _backend_module(backend)
    return Sum(
        module.Matern32(
            st.halfnorm(scale=GP_AMPLITUDE_SCALE),
            st.halfnorm(scale=GP_LENGTH_SCALE),
            amplitude_unit=u.Jy,
            length_scale_unit=u.micron,
        ),
        module.SHO(
            st.halfnorm(scale=GP_AMPLITUDE_SCALE),
            st.halfnorm(scale=FRINGE_PERIOD),
            FRINGE_QUALITY,
            amplitude_unit=u.Jy,
            period_unit=u.micron,
        ),
        labels=("smooth", "fringe"),
    )


def build_fringing_problem(
    data: SyntheticSpectrum,
    *,
    backend: str = "reference",
    seed: int = SEED,
) -> FittingProblem:
    """The flexible likelihood over :func:`build_fringing_kernel`, on *backend*."""
    module = _backend_module(backend)
    likelihood = Likelihood(
        GaussianFamily(),
        module.GaussianProcessNoise(build_fringing_kernel(backend), module.QuasisepGP()),
    )
    dataset = Dataset(data.container(), likelihood=likelihood, label=DATASET_LABEL)
    return FittingProblem(model_for(backend, data.wavelength), [dataset], seed=seed)


def _stationary_problem(
    data: SyntheticSpectrum, kind: str, *, backend: str, seed: int
) -> FittingProblem:
    likelihood = build_likelihood(kind, backend=backend, kernel="matern32", solver="quasisep")
    dataset = Dataset(data.container(), likelihood=likelihood, label=DATASET_LABEL)
    return FittingProblem(model_for(backend, data.wavelength), [dataset], seed=seed)


@dataclasses.dataclass(frozen=True)
class Comparison:
    """Three fits of one spectrum, reduced to what the demonstration reports.

    Attributes
    ----------
    scenario, size
        Which spectrum was fitted, and at how many points.
    standard, stationary, composed
        ``{parameter name: Summary}`` for the standard likelihood, the flexible
        likelihood with a stationary Matérn-3/2, and the flexible likelihood
        with ``Matern32 + SHO``.
    """

    scenario: str
    size: int
    standard: dict[str, Summary]
    stationary: dict[str, Summary]
    composed: dict[str, Summary]

    @property
    def worst_bias(self) -> dict[str, float]:
        """The largest |bias| in posterior widths, per likelihood."""
        return {
            name: max(
                float(summary.bias_in_widths or 0.0)
                for summary in getattr(self, name).values()
                if summary.truth is not None
            )
            for name in ("standard", "stationary", "composed")
        }

    @property
    def coverage(self) -> dict[str, int]:
        """How many of the physical parameters' 68 % intervals cover the truth."""
        return {
            name: sum(
                1
                for summary in getattr(self, name).values()
                if summary.truth is not None and summary.covers_truth
            )
            for name in ("standard", "stationary", "composed")
        }

    def table(self) -> str:
        """The comparison as plain text, one row per physical parameter."""
        header = (
            f"M2 {self.scenario!r}, {self.size} points — bias in posterior widths "
            f"(68 % coverage in brackets)\n"
            f"{'parameter':<22}{'standard':>14}{'Matern32':>14}{'Matern32 + SHO':>18}"
        )
        rows = [header, "-" * 68]
        for name in PHYSICAL_NAMES:
            cells = []
            for which in ("standard", "stationary", "composed"):
                summary = getattr(self, which)[name]
                offset = summary.bias_in_widths
                covered = "y" if summary.covers_truth else "n"
                cells.append(f"{offset:.2f} [{covered}]" if offset is not None else "-")
            rows.append(f"{name:<22}{cells[0]:>14}{cells[1]:>14}{cells[2]:>18}")
        worst = self.worst_bias
        covered = self.coverage
        rows.append("-" * 68)
        rows.append(
            f"{'worst bias':<22}{worst['standard']:>14.2f}"
            f"{worst['stationary']:>14.2f}{worst['composed']:>18.2f}"
        )
        rows.append(
            f"{'covered / 4':<22}{covered['standard']:>14d}"
            f"{covered['stationary']:>14d}{covered['composed']:>18d}"
        )
        return "\n".join(rows)


def compare(
    scenario: str = "strong_smooth",
    *,
    size: int = 200,
    backend: str = "reference",
    budget: EmceeBudget | None = None,
    seed: int = SEED,
    progress: bool = False,
) -> Comparison:
    """Fit one spectrum three ways and reduce the posteriors to :class:`Comparison`.

    One spectrum, generated once and shared, so the three fits differ in their
    likelihood and in nothing else — which is the same discipline
    ``study.run_study`` keeps, for the same reason.
    """
    spectrum = generate(scenario, size=size)
    chosen = budget or TEST_EMCEE
    fits = {
        "standard": _stationary_problem(spectrum, "standard", backend=backend, seed=seed),
        "stationary": _stationary_problem(spectrum, "flexible", backend=backend, seed=seed),
        "composed": build_fringing_problem(spectrum, backend=backend, seed=seed),
    }
    reduced: dict[str, dict[str, Summary]] = {}
    for key, problem in fits.items():
        stored = run(problem, chosen, progress=progress)
        reduced[key] = summarise(stored, names=PHYSICAL_NAMES)
    return Comparison(
        scenario=scenario,
        size=size,
        standard=reduced["standard"],
        stationary=reduced["stationary"],
        composed=reduced["composed"],
    )


def main(argv: Any = None) -> int:
    """Print the comparison table. Informational: nothing here asserts."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--scenario", default="strong_smooth")
    parser.add_argument("--size", type=int, default=200)
    parser.add_argument("--backend", default="reference")
    parser.add_argument("--walkers", type=int, default=TEST_EMCEE.walkers)
    parser.add_argument("--steps", type=int, default=TEST_EMCEE.steps)
    parser.add_argument("--burn-in", type=int, default=TEST_EMCEE.burn_in)
    parser.add_argument("--progress", action="store_true")
    args = parser.parse_args(argv)
    comparison = compare(
        args.scenario,
        size=args.size,
        backend=args.backend,
        budget=EmceeBudget(args.walkers, args.steps, args.burn_in),
        progress=args.progress,
    )
    print(comparison.table())
    return 0


if __name__ == "__main__":  # pragma: no cover - CLI
    raise SystemExit(main())
