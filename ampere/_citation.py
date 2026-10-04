"""How to cite ampere: the one place the citation's facts live.

:data:`CITATION` holds them; :func:`cite` prints them as text or BibTeX, and
``scripts/write_citation_cff.py`` writes the repository's ``CITATION.cff``
(which GitHub's "Cite this repository" button and Zenodo read) from the same
dictionary. ``tests/core/test_citation.py`` asserts that the committed file
equals the generated one, so the two cannot drift: change a fact here, then
run the script.

Two facts are not known yet and are left empty on purpose rather than
guessed. The **DOI** is minted by Zenodo when the first GitHub release is made
from the ``v1.0.0b1`` tag; the concept DOI then goes into ``"doi"`` here, into
``CITATION.cff`` (by the script) and into a badge on the README, in one commit.
The **paper** reference goes into ``"paper"`` once there is a paper.
"""

from __future__ import annotations

import sys
from typing import Literal, TextIO, TypedDict


class Author(TypedDict):
    """One author, in the Citation File Format's name fields."""

    given_names: str
    family_names: str
    email: str | None


class Citation(TypedDict):
    """The facts a citation of ampere is built from."""

    title: str
    message: str
    version: str
    date_released: str | None
    doi: str | None
    license: str
    repository_code: str
    url: str
    distribution: str
    authors: tuple[Author, ...]
    paper: str | None


def _author(given: str, family: str, email: str | None = None) -> Author:
    return {"given_names": given, "family_names": family, "email": email}


#: The citation's facts. The authors are ``pyproject.toml``'s
#: ``[[project.authors]]``, in that order; no ORCID is given because none is
#: stated anywhere in the repository. ``date_released`` and ``doi`` are filled
#: on the tag day and after the first GitHub release respectively (the release
#: procedure in ``docs/development.md``); ``paper`` when there is one.
CITATION: Citation = {
    "title": "ampere",
    "message": (
        "If you use ampere in your research, please cite the software as below "
        "(and the paper, once there is one)."
    ),
    "version": "1.0.0b1",
    "date_released": None,
    "doi": None,
    "license": "GPL-3.0-or-later",
    "repository_code": "https://github.com/ICSM/ampere",
    "url": "https://ampere.readthedocs.io/",
    "distribution": "ampere-astro",
    "authors": (
        _author("Peter", "Scicluna", "peter.scicluna@eso.org"),
        _author("Francisca", "Kemper"),
        _author("Sundar", "Srinivasan"),
        _author("Jonathan", "Marshall"),
        _author("Oscar", "Morata"),
        _author("Alfonso", "Trejo"),
        _author("Sascha", "Zeegers"),
        _author("Lapo", "Fanciullo"),
        _author("Thavisha", "Dharmawardena"),
    ),
    "paper": None,
}

#: What the DOI line says until Zenodo has minted one.
DOI_PENDING = "minted at the first release — see the badge on the README"


def _text(citation: Citation) -> str:
    authors = ", ".join(f"{a['family_names']}, {a['given_names'][0]}." for a in citation["authors"])
    doi = f"https://doi.org/{citation['doi']}" if citation["doi"] else DOI_PENDING
    year = f" ({citation['date_released'][:4]})" if citation["date_released"] else ""
    # The last initial already ends in a full stop; a year in brackets does not.
    stop = "." if year else ""
    lines = [
        (
            f"{authors}{year}{stop} {citation['title']} (version {citation['version']}) "
            f"[software]. {citation['repository_code']}"
        ),
        f"Documentation: {citation['url']}",
        f"Install: pip install {citation['distribution']}",
        f"DOI: {doi}",
    ]
    if citation["paper"]:
        lines.append(f"Please also cite the paper: {citation['paper']}")
    return "\n".join(lines)


def _bibtex(citation: Citation) -> str:
    authors = " and ".join(f"{a['family_names']}, {a['given_names']}" for a in citation["authors"])
    fields = [
        ("author", authors),
        ("title", citation["title"]),
        ("version", citation["version"]),
        ("url", citation["repository_code"]),
        ("license", citation["license"]),
    ]
    if citation["date_released"]:
        fields.append(("date", citation["date_released"]))
    if citation["doi"]:
        fields.append(("doi", citation["doi"]))
    else:
        fields.append(("note", f"DOI {DOI_PENDING}"))
    body = ",\n".join(f"  {key} = {{{value}}}" for key, value in fields)
    entry = f"@software{{ampere,\n{body}\n}}"
    if citation["paper"]:
        entry += f"\n\n% Please also cite the paper:\n% {citation['paper']}"
    return entry


def cite(format: Literal["text", "bibtex"] = "text", *, file: TextIO | None = None) -> None:  # noqa: A002
    """Print how to cite ampere.

    Parameters
    ----------
    format
        ``"text"`` (the default) for a reference to paste into a paper's
        bibliography, or ``"bibtex"`` for a biblatex ``@software`` entry.
    file
        Where to print it; standard output by default.

    Raises
    ------
    ValueError
        If ``format`` is neither ``"text"`` nor ``"bibtex"``.

    Examples
    --------
    >>> import ampere
    >>> ampere.cite()  # doctest: +ELLIPSIS
    Scicluna, P., Kemper, F., ... ampere (version 1.0.0b1) [software]. https://github.com/ICSM/ampere
    Documentation: https://ampere.readthedocs.io/
    Install: pip install ampere-astro
    DOI: ...
    """
    if format == "text":
        rendered = _text(CITATION)
    elif format == "bibtex":
        rendered = _bibtex(CITATION)
    else:
        raise ValueError(f"format must be 'text' or 'bibtex', not {format!r}")
    print(rendered, file=sys.stdout if file is None else file)
