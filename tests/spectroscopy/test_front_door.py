"""``ampere.spectroscopy``: the front door's exports, and ``import ampere`` untouched (W7.3)."""

from __future__ import annotations

import subprocess
import sys

import ampere.backends.reference as reference
import ampere.core as core
import ampere.spectroscopy as spectroscopy
from ampere.spectroscopy import jwst

#: Where each re-export is defined; the door re-exports, it does not redefine.
HOMES = {
    "Spectrum": core,
    "PhotometricPoints": core,
    "LineFluxes": core,
    "GaussianFamily": core,
    "StudentTFamily": core,
    "CauchyFamily": core,
    "PoissonFamily": core,
    "Resample": reference,
    "LSFConvolution": reference,
    "CalibrationScale": reference,
    "SyntheticPhotometry": reference,
    "LineIntegration": reference,
    "read_jwst": jwst,
    "JWSTError": jwst,
}


def test_the_exports_are_the_list_and_the_same_objects() -> None:
    assert sorted(spectroscopy.__all__) == sorted(HOMES)
    for name, home in HOMES.items():
        assert getattr(spectroscopy, name) is getattr(home, name), name


def test_the_interferometric_families_and_the_noise_models_stay_out() -> None:
    for name in ("ComplexGaussianFamily", "RiceFamily", "VonMisesFamily", "GaussianProcessNoise"):
        assert not hasattr(spectroscopy, name), name


def test_import_ampere_does_not_import_the_door() -> None:
    probe = subprocess.run(
        [
            sys.executable,
            "-c",
            (
                "import sys\n"
                "import ampere\n"
                "print(sorted(m for m in sys.modules if m.startswith('ampere.spectroscopy')))\n"
            ),
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    assert probe.stdout.strip() == "[]"
