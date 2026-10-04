Citing ampere
=============

If ampere contributes to work you publish, please cite the software. The
citation is built into the package, so it always names the version you have:

.. code-block:: python

    >>> import ampere
    >>> ampere.cite()                  # a reference for your bibliography
    >>> ampere.cite(format="bibtex")   # a biblatex @software entry

For the ``1.0.0b1`` beta it prints:

.. code-block:: text

    Scicluna, P., Kemper, F., Srinivasan, S., Marshall, J., Morata, O., Trejo, A.,
    Zeegers, S., Fanciullo, L., Dharmawardena, T. ampere (version 1.0.0b1)
    [software]. https://github.com/ICSM/ampere
    Documentation: https://ampere.readthedocs.io/
    Install: pip install ampere-astro
    DOI: minted at the first release — see the badge on the README

**The DOI.** The software is archived on `Zenodo <https://zenodo.org>`_, which
mints a DOI for each GitHub release of the repository and a *concept* DOI
that always resolves to the latest one. Until the first release has been made
the DOI line says so, as above; from then on :func:`ampere.cite` prints the
concept DOI and the README carries it as a badge. Cite the concept DOI unless
you need to pin the exact version you used, in which case the version's own
DOI is on the Zenodo record.

**The paper.** There is no ampere paper yet. When there is one, :func:`ampere.cite`
will ask you to cite it too and print its reference.

**On GitHub**, the repository's *Cite this repository* button reads the same
facts from ``CITATION.cff`` at the repository's root. That file is generated
from the dictionary :func:`ampere.cite` reads (``scripts/write_citation_cff.py``),
and a test fails if the two disagree, so the button, Zenodo and the function
always say the same thing.

The software is licensed under the GNU General Public License, version 3 or
later (``GPL-3.0-or-later``).

.. autofunction:: ampere.cite
