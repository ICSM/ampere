"""``ampere.cite()`` and ``CITATION.cff`` (W6.5, issue #62).

The citation's facts live once, in ``ampere/_citation.py``; ``CITATION.cff``
is generated from them by ``scripts/write_citation_cff.py``. These rows hold
the committed file to the generator's output and the printed citation to the
same facts, so neither can drift from the other or from ``pyproject.toml``.
"""

from __future__ import annotations

import importlib.util
import io
import tomllib
from pathlib import Path
from types import ModuleType

import pytest

import ampere
from ampere._citation import CITATION, DOI_PENDING

ROOT = Path(__file__).resolve().parents[2]


def _generator() -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        "write_citation_cff", ROOT / "scripts" / "write_citation_cff.py"
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _printed(**kwargs: object) -> str:
    buffer = io.StringIO()
    ampere.cite(file=buffer, **kwargs)  # type: ignore[arg-type]
    return buffer.getvalue()


class TestCitationFile:
    def test_the_committed_file_is_the_generated_one(self) -> None:
        committed = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
        assert committed == _generator().render(), (
            "CITATION.cff is stale: run `pixi run python scripts/write_citation_cff.py`"
        )

    def test_the_authors_are_pyproject_s_in_order(self) -> None:
        pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
        declared = [a["name"] for a in pyproject["project"]["authors"]]
        cited = [f"{a['given_names']} {a['family_names']}" for a in CITATION["authors"]]
        assert cited == declared

    def test_the_licence_is_pyproject_s(self) -> None:
        pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
        assert CITATION["license"] == pyproject["project"]["license"]
        assert CITATION["distribution"] == pyproject["project"]["name"]

    def test_it_is_cff_1_2_software_with_the_concept_doi(self) -> None:
        text = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
        assert "cff-version: 1.2.0\n" in text
        assert "type: software\n" in text
        # Zenodo minted it at the first GitHub release (2026-10-05); the concept
        # DOI, which always resolves to the latest version, is the one cited.
        assert 'doi: "10.5281/zenodo.23151412"\n' in text
        assert 'date-released: "2026-10-05"\n' in text


class TestCite:
    def test_text_names_the_version_the_repository_and_the_doi(self) -> None:
        text = _printed()
        assert text.startswith("Scicluna, P., Kemper, F., Srinivasan, S.")
        assert "(2026). ampere (version 1.0.0b1) [software]" in text
        assert "https://github.com/ICSM/ampere" in text
        assert "pip install ampere-astro" in text
        assert "DOI: https://doi.org/10.5281/zenodo.23151412" in text
        assert DOI_PENDING not in text

    def test_bibtex_is_one_software_entry(self) -> None:
        entry = _printed(format="bibtex")
        assert entry.startswith("@software{ampere,\n")
        assert entry.rstrip().endswith("}")
        assert "author = {Scicluna, Peter and Kemper, Francisca and" in entry
        assert "version = {1.0.0b1}" in entry
        assert "date = {2026-10-05}" in entry
        assert "doi = {10.5281/zenodo.23151412}" in entry

    def test_it_prints_to_standard_output_by_default(
        self, capsys: pytest.CaptureFixture[str]
    ) -> None:
        assert ampere.cite() is None
        assert capsys.readouterr().out == _printed()

    def test_an_unknown_format_is_refused(self) -> None:
        with pytest.raises(ValueError, match="'text' or 'bibtex'"):
            ampere.cite(format="ris")  # type: ignore[arg-type]

    def test_the_paper_is_named_once_there_is_one(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setitem(CITATION, "paper", "Scicluna et al. (in prep.)")
        assert "Please also cite the paper: Scicluna et al. (in prep.)" in _printed()
        assert "Scicluna et al. (in prep.)" in _printed(format="bibtex")
