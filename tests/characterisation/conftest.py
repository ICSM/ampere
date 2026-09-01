"""Shared fixtures/helpers for the ampere v2 characterisation suite (W0.4).

This suite wraps the frozen legacy flows demonstrated in
``examples/minimal_working_example*.py`` (emcee, dynesty, zeus, SBI) with
fixed seeds and small iteration budgets, so that later PRs have a
"legacy still works" regression check. Per ``AGENTS.md`` / ``WORK_ITEMS.md``
W0.4, this suite must NOT modify anything under ``ampere/`` or the example
scripts themselves -- it only imports/exercises the real, frozen legacy
classes (``ampere.data.Spectrum``/``Photometry``,
``ampere.infer.emceesearch.EmceeSearch``,
``ampere.infer.dynestysearch.DynestyNestedSampler``,
``ampere.infer.zeussearch.ZeusSearch``, ``ampere.infer.sbi.SBI_SNPE``).

The model class and synthetic-data recipe below are a lightly parametrised
*copy* of the ``ASimpleModel`` / data-synthesis code shared by all four
``examples/minimal_working_example*.py`` scripts (allowed per the W0.4 work
item: "you may lightly parametrise copies of the example flows inside
tests/"). Differences from the original examples, and why:

* The spectrum is subsampled (every ``SPECTRUM_STRIDE``-th point of the
  CASSIS example spectrum) purely to keep the dense-GP noise-model
  likelihood (``ampere/data/spectrum.py``, O(N^3) per evaluation) cheap
  enough for a fast CI-friendly test budget; this does not touch the
  covariance-matrix code path itself, only the number of points fed into
  it.
* ``ASimpleModel.__init__`` here takes its wavelength grid from its own
  argument rather than (as in the example scripts, where it is a bug of no
  practical consequence when run as ``__main__``) implicitly reading the
  module-global ``wavelengths`` variable.
"""

from __future__ import annotations

import importlib.util
import itertools
from pathlib import Path

import numpy as np
import pyphot
import pytest
from spectres import spectres

import ampere
from ampere.data import Photometry, Spectrum
from ampere.models import Model

# ---------------------------------------------------------------------------
# Fixed problem definition
# ---------------------------------------------------------------------------

#: Base seed for all characterisation tests. Fixed so golden values are
#: reproducible; each test derives its own sub-seeds from this so that
#: tests remain independent of execution order.
SEED = 20260901

#: Every N-th point of the CASSIS example spectrum is kept (see module
#: docstring for why).
SPECTRUM_STRIDE = 3

TRUE_SLOPE = 1.0
TRUE_INTERCEPT = 1.0

REPO_ROOT = Path(__file__).resolve().parents[2]
SPECTRUM_FILE = (
    REPO_ROOT / "examples" / "test_data" / "cassis_yaaar_spcfw_14191360t.fits"
)
FILTER_NAMES = np.array(["WISE_RSR_W1", "SPITZER_MIPS_70"])


class LinearModel(Model):
    """A lightly-cleaned copy of ``examples/minimal_working_example*.py``'s
    ``ASimpleModel``: flux linear in wavelength, two free parameters.
    """

    def __init__(
        self,
        wavelengths,
        flatprior=True,
        lims=np.array([[-10, 10], [-10, 10]]),
    ):
        self.wavelength = wavelengths
        self.npars = 2
        self.npars_ptform = 2
        self.lims = lims
        self.flatprior = flatprior
        self.parLabels = ["slope", "intercept"]

    def __call__(self, slope, intercept, **kwargs):
        self.modelFlux = slope * self.wavelength + intercept
        return {"spectrum": {"wavelength": self.wavelength, "flux": self.modelFlux}}

    def lnprior(self, theta, **kwargs):
        if not self.flatprior:
            raise NotImplementedError()
        if (self.lims[0, 0] < theta[0] < self.lims[0, 1]) and (
            self.lims[1, 0] < theta[1] < self.lims[1, 1]
        ):
            return 0
        return -np.inf

    def prior_transform(self, u, **kwargs):
        if not self.flatprior:
            raise NotImplementedError()
        return (self.lims[:, 1] - self.lims[:, 0]) * u + self.lims[:, 0]


def build_linear_sed_problem(seed=SEED, stride=SPECTRUM_STRIDE):
    """Build the (model, dataset) pair shared by all four characterisation
    flows, following the recipe in ``examples/minimal_working_example*.py``.

    Seeds ``numpy``'s legacy global RNG (``np.random.seed``) so the
    synthetic photometry/spectrum noise realisation is reproducible.
    Returns a dict with the model, the [photometry, spectrum] dataset, and
    the true parameter values used to generate the synthetic data.

    Raises ``AttributeError`` (uncaught) if the installed ``pyphot`` is
    incompatible with the legacy filter-convolution API used here and by
    ``ampere.data.photometry`` -- see ``PYPHOT_LEGACY_COMPATIBLE`` below.
    """
    np.random.seed(seed)

    wavelengths = 10 ** np.linspace(0.0, 1.9, 2000)
    model = LinearModel(wavelengths)
    model(TRUE_SLOPE, TRUE_INTERCEPT)
    model_flux = model.modelFlux

    libdir = ampere.__file__.strip("__init__.py")
    libname = f"{libdir}ampere_allfilters.hd5"
    filter_library = pyphot.get_library(fname=libname)
    filters = filter_library.load_filters(
        FILTER_NAMES, interp=True, lamb=wavelengths * pyphot.unit["micron"]
    )
    flam = model_flux / wavelengths**2
    mod_sed = []
    for f in filters:
        lp = f.lpivot.to("micron").value
        fphot = f.get_flux(
            wavelengths * pyphot.unit["micron"], flam * pyphot.unit["flam"], axis=-1
        ).value
        mod_sed.append(fphot * lp**2)
    mod_sed = np.array(mod_sed)
    phot_unc = 0.1 * mod_sed
    mod_sed = mod_sed + np.random.randn(len(FILTER_NAMES)) * phot_unc

    irs_example = Spectrum.fromFile(str(SPECTRUM_FILE), format="SPITZER-YAAAR")
    wave1 = irs_example[1].wavelength[::stride]
    spec1_flux = spectres(wave1, wavelengths, model_flux)
    unc1 = 0.1 * spec1_flux
    spec1_flux = spec1_flux + np.random.randn(len(spec1_flux)) * unc1
    spectrum = Spectrum(
        wave1,
        spec1_flux,
        unc1,
        "um",
        "Jy",
        calUnc=0.0025,
        scaleLengthPrior=0.01,
    )
    spectrum.setResampler(resampleMethod="fast")

    photometry = Photometry(
        filterName=FILTER_NAMES,
        value=mod_sed,
        uncertainty=phot_unc,
        photunits="Jy",
        libName=libname,
    )
    photometry.reloadFilters(wavelengths)

    return {
        "model": model,
        "dataset": [photometry, spectrum],
        "wavelengths": wavelengths,
        "true_slope": TRUE_SLOPE,
        "true_intercept": TRUE_INTERCEPT,
    }


def seed_default_rng(monkeypatch, seed=SEED):
    """Make every call to ``numpy.random.default_rng()`` inside legacy
    ampere code deterministic for the duration of a test.

    Several ``ampere.infer`` / ``ampere.data`` code paths call
    ``np.random.default_rng()`` with *no* seed -- e.g. the pre-optimisation
    starting-point selection in
    ``ampere/infer/emceesearch.py::EmceeSearch.optimise`` and
    ``ampere/infer/zeussearch.py::ZeusSearch.optimise``
    (both around their ``preopt`` branch), and
    ``ampere/infer/sbi.py::SBI_SNPE.sample`` /
    ``_check_prior_normalisation`` / ``mixins.py::SBIPostProcessor``. That
    makes repeated runs of those code paths non-reproducible even with
    ``np.random.seed()`` set (see W0.4 report: out-of-scope finding). We
    intercept the *global* ``numpy.random.default_rng`` function (not any
    ampere code -- this is standard ``pytest`` ``monkeypatch`` usage) so
    that two independent test runs are reproducible.
    """
    counter = itertools.count()
    original_default_rng = np.random.default_rng

    def _fixed_default_rng(*_args, **_kwargs):
        return original_default_rng(seed + next(counter))

    monkeypatch.setattr(np.random, "default_rng", _fixed_default_rng)


# ---------------------------------------------------------------------------
# Optional / broken dependency handling
# ---------------------------------------------------------------------------

#: True if the installed pyphot still exposes the pint-based ``pyphot.unit``
#: API that ``ampere.data.photometry`` (legacy, frozen) calls directly in
#: ``Photometry.reloadFilters`` and ``Photometry.lnlike``
#: (``ampere/data/photometry.py`` lines ~421-422, 476-477, 590-602).
#: pyphot >= ~2.0 removed this attribute as part of a unit-adapter rework;
#: pip resolves the unpinned ``pyphot`` dependency in ``pyproject.toml`` to
#: the latest release, so a plain ``pip install -e ".[dev]"`` currently
#: installs an incompatible pyphot. See the W0.4 report for full details.
PYPHOT_LEGACY_COMPATIBLE = hasattr(pyphot, "unit")

xfail_if_pyphot_incompatible = pytest.mark.xfail(
    condition=not PYPHOT_LEGACY_COMPATIBLE,
    reason=(
        "ampere.data.photometry.Photometry (legacy, frozen) calls "
        "pyphot.unit['micron'/'flam'] directly in reloadFilters()/lnlike() "
        "(ampere/data/photometry.py:421-422,476-477,590-602). The installed "
        f"pyphot ({getattr(pyphot, '__VERSION__', 'unknown version')}) "
        "removed the module-level 'unit' attribute in its unit-adapter "
        "rework (pyphot>=~2.0), which is what 'pip install -e \".[dev]\"' "
        "resolves for the unpinned pyproject.toml dependency. This blocks "
        "every SED-fitting flow that uses Photometry, not just this one. "
        "Not fixed here (frozen legacy code; pin/adapter fix is out of "
        "scope for W0.4) -- see the W0.4 report."
    ),
    strict=False,
    raises=AttributeError,
)


def assert_means_within_tolerance(means, golden, tolerances, labels):
    """Assert each entry of ``means`` is within the matching ``tolerances``
    entry of the matching ``golden`` entry, raising a single readable
    ``AssertionError`` listing every parameter that is out of tolerance
    (rather than stopping at the first one).
    """
    means = np.asarray(means, dtype=float)
    golden = np.asarray(golden, dtype=float)
    tolerances = np.asarray(tolerances, dtype=float)
    diffs = np.abs(means - golden)
    bad = [
        f"{label}: got {m:.5f}, golden {g:.5f} +/- {t:.5f} (|diff|={d:.5f})"
        for label, m, g, d, t in zip(labels, means, golden, diffs, tolerances)
        if d > t
    ]
    assert not bad, "Posterior mean(s) outside golden tolerance:\n" + "\n".join(bad)


HAVE_ZEUS = importlib.util.find_spec("zeus") is not None
HAVE_TORCH = importlib.util.find_spec("torch") is not None
HAVE_SBI = HAVE_TORCH and importlib.util.find_spec("sbi") is not None

requires_zeus = pytest.mark.skipif(
    not HAVE_ZEUS,
    reason="zeus-mcmc is not installed (optional extra: ampere[zeus])",
)
requires_sbi = pytest.mark.skipif(
    not HAVE_SBI,
    reason="torch and/or sbi are not installed (optional extra: ampere[sbi])",
)
