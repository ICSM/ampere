"""The adapter and the two trials of W6.9 (see the package docstring).

Everything that touches ``robusta_hmf`` is here, behind :func:`require_rhmf`,
so importing this module (and so the package) needs neither JAX nor the
``rhmf`` extra; the adapter functions raise
:class:`~ampere.core.exceptions.OptionalDependencyError` on *use*.

The adapter — :func:`to_matrix`, :func:`fit_rhmf`, :func:`anomaly_score` — is
what ``ampere.diagnostics.rhmf`` would have been (``diagnostics.md`` §2.3,
§2.4); it is kept in the example because the maturity gate of
``DEVELOPMENT_PLAN.md`` §2 is unmet and nothing under ``ampere/`` may import
the dependency. It adds no default for ``rank`` or ``robust_scale``: both are
required keyword arguments.
"""

from __future__ import annotations

import contextlib
import dataclasses
import importlib
import importlib.metadata
import io
import time
from collections.abc import Sequence
from typing import Any

import numpy as np

from ampere.core import AnomalyScore, Spectrum
from ampere.core.exceptions import OptionalDependencyError

#: The ``robusta-hmf`` commit the trial was run against: the upstream ``paper``
#: tag (2026-09-15). The pin lives in ``pyproject.toml``'s ``rhmf`` pixi feature.
COMMIT = "cf2fcffa10bc6b150b6cdee76cb18980bbe908ef"

#: ``AnomalyScore.provenance`` for everything this adapter produces.
PROVENANCE = "rhmf_prefit"


def require_rhmf() -> Any:
    """Import and return the ``robusta_hmf`` package, or name the extra that is missing.

    Also switches JAX to float64: the weights here are inverse variances of
    order 1e4 against flux of order one, which is uncomfortable in float32.
    A script may do that on its own behalf; ``ampere`` itself never flips
    ``jax_enable_x64`` for a user (``lowering.md`` §10.2(a)).
    """
    try:
        module = importlib.import_module("robusta_hmf")
        import jax
    except ImportError as error:
        raise OptionalDependencyError(
            "robusta_hmf", extra="rhmf", context="running the RHMF exploratory trial"
        ) from error
    jax.config.update("jax_enable_x64", True)
    return module


def installed_version() -> str:
    """The installed ``robusta-hmf`` version string (a dev version for a git pin)."""
    try:
        return importlib.metadata.version("robusta-hmf")
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


# ---------------------------------------------------------------------------
# The adapter
# ---------------------------------------------------------------------------


def to_matrix(spectra: Sequence[Spectrum]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Aligned spectra to ``(Y, W, coordinates)``: values, inverse variances, the shared axis.

    ``diagnostics.md`` §2.3's alignment step is **checked, not performed**: the
    spectra must already share one coordinate axis (the M2 generators put every
    scenario on one grid at a given size), and a mismatch raises rather than
    resampling — interpolation is a choice the caller should see, and this trial
    has no use for it. An ampere mask becomes zero weight (Lesson R1); a missing
    or non-positive uncertainty is refused, because RHMF is *heteroskedastic* and
    a unit weight would be a silent fabrication.
    """
    if len(spectra) == 0:
        raise ValueError("to_matrix needs at least one spectrum.")
    coordinates = np.asarray(spectra[0].spectral_axis.values, dtype=float)
    values, weights = [], []
    for index, spectrum in enumerate(spectra):
        axis = np.asarray(spectrum.spectral_axis.values, dtype=float)
        if axis.shape != coordinates.shape or not np.allclose(axis, coordinates):
            raise ValueError(
                f"spectrum {index} is not on the first spectrum's coordinate axis; RHMF "
                f"factorises a matrix with one feature axis, so resample the collection onto a "
                f"common grid first (diagnostics.md §2.3)."
            )
        if spectrum.uncertainty is None:
            raise ValueError(f"spectrum {index} carries no uncertainty; RHMF needs weights.")
        sigma = np.asarray(spectrum.uncertainty, dtype=float)
        if not np.all(sigma > 0.0):
            raise ValueError(f"spectrum {index} has non-positive uncertainties.")
        weight = 1.0 / sigma**2
        if spectrum.mask is not None:
            weight = np.where(np.asarray(spectrum.mask, dtype=bool), 0.0, weight)
        values.append(np.asarray(spectrum.values, dtype=float))
        weights.append(weight)
    return np.vstack(values), np.vstack(weights), coordinates


@dataclasses.dataclass(frozen=True)
class RHMFFit:
    """One fitted factorisation, with what :func:`anomaly_score` needs to read it."""

    model: Any
    state: Any
    loss_history: np.ndarray
    Y: np.ndarray
    W: np.ndarray
    rank: int
    robust_scale: float
    robust_nu: float
    seconds: float

    @property
    def iterations(self) -> int:
        """Iterations run (the loss history's length is one more than ``state.it``)."""
        return int(self.state.it)

    def weights(self) -> np.ndarray:
        """The IRLS robust weights, in ``(0, 1]``, shape of ``Y``."""
        return np.asarray(self.model.robust_weights(self.Y, self.W, self.state), dtype=float)

    def score(self) -> np.ndarray:
        """``1 - weights``: higher means less well explained. Zero where ``W == 0``."""
        return 1.0 - self.weights()


def fit_rhmf(
    Y: np.ndarray,
    W: np.ndarray,
    *,
    rank: int,
    robust_scale: float,
    robust_nu: float = 1.0,
    max_iter: int = 500,
    seed: int = 0,
) -> RHMFFit:
    """Fit ``Robusta(rank=..., robust=True, robust_scale=...)`` to ``(Y, W)``.

    ``rank`` and ``robust_scale`` are required: ``diagnostics.md`` §2.5 rules out
    a silent default, and this trial has none to offer. ``robust_scale`` is in
    units of the declared sigma (the likelihood sees ``W * r**2``). ``robust_nu``
    is left at the package's own default of 1 (a Cauchy-like tail) and is not
    scanned. RHMF prints its progress to standard output; it is captured here.
    """
    module = require_rhmf()
    model = module.Robusta(
        rank=int(rank), robust=True, robust_scale=float(robust_scale), robust_nu=robust_nu
    )
    started = time.perf_counter()
    with contextlib.redirect_stdout(io.StringIO()):
        state, history = model.fit(np.asarray(Y), np.asarray(W), max_iter=int(max_iter), seed=seed)
    seconds = time.perf_counter() - started
    return RHMFFit(
        model=model,
        state=state,
        loss_history=np.asarray(history),
        Y=np.asarray(Y),
        W=np.asarray(W),
        rank=int(rank),
        robust_scale=float(robust_scale),
        robust_nu=float(robust_nu),
        seconds=seconds,
    )


def anomaly_score(
    fit: RHMFFit,
    coordinates: np.ndarray,
    *,
    row: int | None = None,
    aggregate: float | None = None,
) -> AnomalyScore:
    """An :class:`~ampere.core.AnomalyScore` from the robust weights.

    * ``row=i`` (the default view): one object's per-feature score ``1 - w``,
      indexed by *coordinates*, masked where the weight was zero.
    * ``aggregate=q``: one score per object, ``1 - quantile_q(w)`` across its
      features — Lesson R2's "a low quantile across features" — indexed by the
      object's position in the collection. ``coordinates`` is ignored.

    Higher is less well explained by the rank-``K`` model *of the collection*.
    ``provenance`` is ``"rhmf_prefit"`` and the notes carry the rank, the scale,
    the commit and the caveat that nothing has been fitted.
    """
    if (row is None) == (aggregate is None):
        raise ValueError("pass exactly one of row= (per feature) or aggregate= (per object).")
    weights = fit.weights()
    notes = (
        f"RHMF pre-fit screen (robusta-hmf {installed_version()}, commit {COMMIT[:12]}): "
        f"rank {fit.rank}, robust_scale {fit.robust_scale:g} sigma, Student-t nu "
        f"{fit.robust_nu:g}. Score = 1 - IRLS robust weight; higher = less well explained by a "
        f"rank-{fit.rank} model of this collection. Nothing has been fitted: it says where a "
        f"low-rank model of the collection struggles, not where any physical model will, and "
        f"rank and robust_scale have no validated default (diagnostics.md 2.5)."
    )
    if aggregate is not None:
        if not 0.0 <= aggregate <= 1.0:
            raise ValueError(f"aggregate is a quantile in [0, 1], got {aggregate!r}.")
        values = 1.0 - np.quantile(weights, aggregate, axis=1)
        return AnomalyScore(
            coordinates=np.arange(weights.shape[0], dtype=float),
            values=values,
            provenance=PROVENANCE,
            interpretation_notes=notes
            + f" Per object: the {aggregate:g} quantile across features.",
        )
    assert row is not None
    coordinates = np.asarray(coordinates, dtype=float)
    if coordinates.shape[0] != weights.shape[1]:
        raise ValueError(f"{coordinates.shape[0]} coordinates for {weights.shape[1]} features.")
    return AnomalyScore(
        coordinates=coordinates,
        values=1.0 - weights[row],
        mask=fit.W[row] == 0.0,
        provenance=PROVENANCE,
        interpretation_notes=notes,
    )


# ---------------------------------------------------------------------------
# Contrast: the trial's one statistic
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class Contrast:
    """Score inside versus outside the injected deviation, deviated rows versus controls."""

    inside: float  # mean score inside the band, deviated rows
    outside: float  # mean score outside it, deviated rows
    ratio: float  # inside / outside, deviated rows
    control_ratio: float  # the same ratio for the control rows, same band
    lift: float  # inside, deviated rows / inside, control rows
    excess: float  # ratio / control_ratio: how much the deviation adds to the contrast


def contrast(
    score: np.ndarray,
    deviated: Sequence[int],
    controls: Sequence[int],
    band: np.ndarray,
    exclude: np.ndarray | None = None,
) -> Contrast:
    """The three ratios of the trial's tables, for a ``(rows, features)`` score matrix.

    *band* is the boolean feature mask of the injected deviation. The control
    ratio is **not** 1: noise alone gives every feature a score of order a half
    at ``robust_scale=1``, and a band can sit where the factorisation happens to
    fit less well, so the deviated rows' ratio is only meaningful beside it.
    *exclude* removes features from the "outside" as well (the image trial drops
    the source's core, which every flattening fits badly whatever the deviation).
    """
    rest = ~band if exclude is None else ~band & ~exclude
    dev = score[np.asarray(deviated)]
    ctl = score[np.asarray(controls)]
    tiny = 1e-12
    inside, outside = float(dev[:, band].mean()), float(dev[:, rest].mean())
    ctl_in, ctl_out = float(ctl[:, band].mean()), float(ctl[:, rest].mean())
    ratio = inside / max(outside, tiny)
    control_ratio = ctl_in / max(ctl_out, tiny)
    return Contrast(
        inside=inside,
        outside=outside,
        ratio=ratio,
        control_ratio=control_ratio,
        lift=inside / max(ctl_in, tiny),
        excess=ratio / max(control_ratio, tiny),
    )


def auc(scores_deviated: np.ndarray, scores_control: np.ndarray) -> float:
    """P(a deviated object's score exceeds a control's): the rank statistic, ties as half."""
    a = np.asarray(scores_deviated, dtype=float)[:, None]
    b = np.asarray(scores_control, dtype=float)[None, :]
    return float((a > b).mean() + 0.5 * (a == b).mean())


# ---------------------------------------------------------------------------
# Trial (a): the M2 spectra
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class SpectraCollection:
    """A collection of M2 spectra: the controls, then each deviated scenario's copies."""

    spectra: list[Spectrum]
    labels: list[str]  # the scenario key of each row; "none" for a control
    wavelength: np.ndarray
    bands: dict[str, np.ndarray]  # scenario key -> boolean mask of the injected deviation
    Y: np.ndarray
    W: np.ndarray

    def rows(self, key: str) -> list[int]:
        """Row indices of the spectra of scenario *key*."""
        return [i for i, label in enumerate(self.labels) if label == key]


def deviation_band(deviation: np.ndarray, *, fraction: float = 0.5) -> np.ndarray:
    """Where ``|delta|`` reaches *fraction* of its peak: the "inside" of the contrast.

    For a line that is its core; for a fringe it is the half of each cycle near
    the peaks, which is a *region of the axis* but not a localised one — the
    fringe has no location, as ``Scenario.localised_at`` says.
    """
    magnitude = np.abs(np.asarray(deviation, dtype=float))
    return magnitude >= fraction * magnitude.max()


def spectra_collection(
    *, controls: int, copies: int, size: int = 200, seed: int = 0
) -> SpectraCollection:
    """*controls* undeviated spectra, then *copies* of each deviated M2 scenario.

    Every row has its own noise seed (``seed`` plus its position), so the copies
    of a scenario share the deviation but not the noise. They share one
    wavelength grid at a given size, so ``diagnostics.md`` §2.3's alignment step
    is trivial here and :func:`to_matrix` only checks it.
    """
    from examples.m2_misspecification.generators import EXTENDED_SCENARIOS, generate

    spectra: list[Spectrum] = []
    labels: list[str] = []
    bands: dict[str, np.ndarray] = {}
    wavelength = None
    order = [("none", controls)] + [(s.key, copies) for s in EXTENDED_SCENARIOS if s.key != "none"]
    for key, count in order:
        for _ in range(count):
            data = generate(key, size=size, seed=seed + len(spectra))
            spectra.append(data.container())
            labels.append(key)
            wavelength = data.wavelength
            if key != "none" and key not in bands:
                bands[key] = deviation_band(data.deviation)
    assert wavelength is not None
    Y, W, coordinates = to_matrix(spectra)
    return SpectraCollection(spectra, labels, coordinates, bands, Y, W)


def scan_spectra(
    collection: SpectraCollection,
    *,
    ranks: Sequence[int],
    scales: Sequence[float],
    max_iter: int,
) -> tuple[list[dict[str, Any]], dict[tuple[int, float], RHMFFit]]:
    """Fit the collection at every ``(rank, robust_scale)`` and score each scenario's band.

    Returns the records (one per grid point and scenario) and the fits, keyed by
    grid point, for the figures.
    """
    controls = collection.rows("none")
    records: list[dict[str, Any]] = []
    fits: dict[tuple[int, float], RHMFFit] = {}
    for rank in ranks:
        for scale in scales:
            fit = fit_rhmf(
                collection.Y, collection.W, rank=rank, robust_scale=scale, max_iter=max_iter
            )
            fits[(rank, scale)] = fit
            score = fit.score()
            for key, band in collection.bands.items():
                rows = collection.rows(key)
                c = contrast(score, rows, controls, band)
                records.append(
                    {
                        "scenario": key,
                        "rank": rank,
                        "robust_scale": scale,
                        "rows": len(rows),
                        "band_points": int(band.sum()),
                        "inside": c.inside,
                        "outside": c.outside,
                        "ratio": c.ratio,
                        "control_ratio": c.control_ratio,
                        "lift": c.lift,
                        "excess": c.excess,
                        "auc_rows": auc(
                            np.quantile(fit.weights()[rows][:, band], 0.1, axis=1) * -1.0,
                            np.quantile(fit.weights()[controls][:, band], 0.1, axis=1) * -1.0,
                        ),
                        "iterations": fit.iterations,
                        "seconds": fit.seconds,
                    }
                )
    return records, fits


def best_by_scenario(records: Sequence[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    """The grid point with the largest ``excess`` for each scenario (an oracle choice).

    It is made *knowing where the deviation is*, so it is an upper bound on what
    a user could pick without that knowledge, not a recommended setting.
    """
    best: dict[str, dict[str, Any]] = {}
    for record in records:
        key = record["scenario"]
        if key not in best or record["excess"] > best[key]["excess"]:
            best[key] = record
    return best


# ---------------------------------------------------------------------------
# Trial (b): W5.5's image
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class ImageCollection:
    """Images of W5.5's truth at varied source parameters, with and without the background.

    ``deviated[i]`` is ``True`` where image *i* carries the smooth background
    the study's misspecified arms omit. ``injected`` is the noiseless background
    (after the PSF) of a reference image, and ``band`` the pixels where it
    exceeds half its peak. Images, sigma and ``injected`` are in units of the
    median per-pixel sigma.
    """

    images: np.ndarray  # (n, pixels, pixels)
    sigma: np.ndarray  # (n,)
    deviated: np.ndarray  # (n,) bool
    injected: np.ndarray  # (pixels, pixels)
    band: np.ndarray  # (pixels, pixels) bool: background-dominated, off the source core
    core: np.ndarray  # (pixels, pixels) bool: the source's core, excluded from "outside"
    pixels: int
    background: float  # the multiple of W5.5's background flux

    @property
    def peak_signal_to_noise(self) -> float:
        """The injected deviation's peak in units of the per-pixel sigma."""
        return float(self.injected.max() / np.median(self.sigma))

    def rows(self, deviated: bool) -> list[int]:
        """Indices of the deviated (or the control) images."""
        return [i for i, flag in enumerate(self.deviated) if bool(flag) == deviated]


def _noiseless_image(pixels: int, *, flux: float, fwhm: float, background: float) -> np.ndarray:
    """One noiseless PSF-convolved image, through W5.5's own model and camera.

    *background* is a multiple of W5.5's background flux; ``0`` omits it.
    """
    from ampere.backends import reference
    from ampere.backends.reference import interferometry
    from ampere.core import negotiate
    from examples.image import generators as gen
    from examples.image.model import SourceWithBackground

    grid = gen.seed_grid()
    kwargs: dict[str, Any] = {"flux": flux, "fwhm": fwhm}
    if background:
        kwargs["background_flux"] = background * gen.BACKGROUND["flux"]
        kwargs["background_fwhm"] = gen.BACKGROUND["fwhm"]
    instrument = gen.camera(reference, gen.observed_template(pixels))
    compiled = SourceWithBackground(interferometry, grid, grid, **kwargs).compile_for(
        negotiate([instrument])
    )
    return np.asarray(instrument(compiled.evaluate()).values, dtype=float)


def image_collection(
    *, pixels: int, count: int, seed: int = 0, spread: float = 0.25, background: float = 1.0
) -> ImageCollection:
    """*count* images, alternately controls and deviated, source flux and FWHM varied by *spread*.

    The source parameters are drawn uniformly within ``+-spread`` of W5.5's
    truth, so the images are not identical copies (a rank-1 factorisation of
    identical copies would be a trivially good model of everything but noise).
    The noise is W5.5's: Gaussian, uniform, ``SIGMA_FRACTION`` of the brightest
    pixel of each noiseless image. *background* multiplies the omitted
    component's flux (``1`` is W5.5's study; a larger value is the sensitivity
    row that separates "the flattening cannot see it" from "it is below the
    noise"). Values are returned in units of the median per-pixel sigma, because
    W5.5's surface brightness is ~1e15 Jy/sr and a weight of ~1e-29 is not a
    number to hand a factorisation.
    """
    from examples.image import generators as gen

    rng = np.random.default_rng(seed)
    images, sigmas, flags = [], [], []
    for index in range(count):
        flux = gen.SOURCE["flux"] * (1.0 + spread * rng.uniform(-1.0, 1.0))
        fwhm = gen.SOURCE["fwhm"] * (1.0 + spread * rng.uniform(-1.0, 1.0))
        deviated = index % 2 == 1
        clean = _noiseless_image(
            pixels, flux=flux, fwhm=fwhm, background=background if deviated else 0.0
        )
        sigma = gen.SIGMA_FRACTION * float(np.max(np.abs(clean)))
        images.append(clean + rng.normal(0.0, sigma, clean.shape))
        sigmas.append(sigma)
        flags.append(deviated)
    reference = gen.SOURCE
    injected = _noiseless_image(
        pixels, flux=reference["flux"], fwhm=reference["fwhm"], background=background
    ) - _noiseless_image(pixels, flux=reference["flux"], fwhm=reference["fwhm"], background=0.0)
    source_only = _noiseless_image(
        pixels, flux=reference["flux"], fwhm=reference["fwhm"], background=0.0
    )
    core = source_only >= 0.05 * source_only.max()
    unit = float(np.median(sigmas))
    return ImageCollection(
        images=np.array(images) / unit,
        sigma=np.array(sigmas) / unit,
        deviated=np.array(flags),
        injected=injected / unit,
        band=(injected >= 0.5 * injected.max()) & ~core,
        core=core,
        pixels=pixels,
        background=background,
    )


def scan_image_collection(
    collection: ImageCollection, *, ranks: Sequence[int], scales: Sequence[float], max_iter: int
) -> tuple[list[dict[str, Any]], dict[tuple[int, float], RHMFFit]]:
    """Flattening 1: each image is one row, each pixel one feature."""
    n, pixels = collection.images.shape[0], collection.pixels
    Y = collection.images.reshape(n, pixels * pixels)
    W = np.repeat((1.0 / collection.sigma**2)[:, None], pixels * pixels, axis=1)
    band, core = collection.band.reshape(-1), collection.core.reshape(-1)
    deviated, controls = collection.rows(True), collection.rows(False)
    records, fits = [], {}
    for rank in ranks:
        for scale in scales:
            fit = fit_rhmf(Y, W, rank=rank, robust_scale=scale, max_iter=max_iter)
            fits[(rank, scale)] = fit
            c = contrast(fit.score(), deviated, controls, band, core)
            low = np.quantile(fit.weights(), 0.1, axis=1)
            records.append(
                {
                    "flattening": "images-as-rows",
                    "background": collection.background,
                    "rank": rank,
                    "robust_scale": scale,
                    "inside": c.inside,
                    "outside": c.outside,
                    "ratio": c.ratio,
                    "control_ratio": c.control_ratio,
                    "excess": c.excess,
                    "auc_rows": auc(-low[deviated], -low[controls]),
                    "iterations": fit.iterations,
                    "seconds": fit.seconds,
                }
            )
    return records, fits


def scan_single_image(
    collection: ImageCollection, *, ranks: Sequence[int], scales: Sequence[float], max_iter: int
) -> tuple[list[dict[str, Any]], dict[tuple[int, float], dict[str, RHMFFit]]]:
    """Flattening 2: one image's rows are the objects and its columns the features.

    One deviated image and one control image (the first of each) are factorised
    separately; the contrast compares the deviated image's score inside the
    background's band with the control image's score in the same pixels. There
    are no "rows of deviated versus control objects" here: the objects are image
    rows, and the band cuts across all of them.
    """
    deviated, controls = collection.rows(True)[0], collection.rows(False)[0]
    band, rest = collection.band, ~collection.band & ~collection.core
    records, fits = [], {}
    for rank in ranks:
        for scale in scales:
            scores, pair = {}, {}
            for name, index in (("deviated", deviated), ("control", controls)):
                Y = collection.images[index]
                W = np.full_like(Y, 1.0 / collection.sigma[index] ** 2)
                fit = fit_rhmf(Y, W, rank=rank, robust_scale=scale, max_iter=max_iter)
                pair[name] = fit
                scores[name] = fit.score()
            fits[(rank, scale)] = pair
            inside_dev, outside_dev = (
                scores["deviated"][band].mean(),
                scores["deviated"][rest].mean(),
            )
            inside_ctl, outside_ctl = scores["control"][band].mean(), scores["control"][rest].mean()
            ratio, control_ratio = inside_dev / outside_dev, inside_ctl / outside_ctl
            records.append(
                {
                    "flattening": "single-image",
                    "background": collection.background,
                    "rank": rank,
                    "robust_scale": scale,
                    "inside": float(inside_dev),
                    "outside": float(outside_dev),
                    "ratio": float(ratio),
                    "control_ratio": float(control_ratio),
                    "excess": float(ratio / control_ratio),
                    "auc_rows": float("nan"),
                    "iterations": pair["deviated"].iterations,
                    "seconds": pair["deviated"].seconds + pair["control"].seconds,
                }
            )
    return records, fits
