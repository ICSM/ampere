"""
Welcome to AMPERE

Ampere is an attempt to produce a fitting environment that can natively
handle multiple different kinds of astronomical data with differing
information content, even when the model being applied to them might be
missing some of the key processes and be unable to actually reproduce all
aspects of the observations properly.
This is motivated by the need within out team for a tool to simultaneously
fit the SEDs and spectra of dusty objects to constrain, among other things,
the dust properties and mineralogy.
However, the final product will be more general, such that it can be applied
to a wide range of astronomical questions.

To achieve this, we include a simple parametric model of correlated noise for
each dataset which is marginalised over when fitting the parameters of the
models.
This approach was chosen because, typically, deficiencies in the model result
in structured residuals, and structured residuals are equivalent to
correlated noise.
This effectively downweights parts of the data which the model doesn't
represent well without having to manually identify these regions.

ampere is in beta, released on PyPI as ``ampere-astro``; the documentation is
at https://ampere.readthedocs.io.
"""

from importlib.metadata import PackageNotFoundError, version as _version

try:
    # The distribution is `ampere-astro` (D13: PyPI's `ampere` is an unrelated
    # package); the import name stays `ampere`.
    __version__ = _version("ampere-astro")
except PackageNotFoundError:  # a checkout on sys.path without an install
    __version__ = "0+unknown"

__copyright__ = """ Copyright (C) 2017-2026  P. Scicluna, F. Kemper, S. Srinivasan,
J.P. Marshall, O. Morata, A. Trejo, S. Zeegers, L. Fanciullo,
T. Dharmawardena

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>
"""
# Add a statement like this to each file/module/subpackage we include.
# Exactly how it should be included is a matter of debate... but including
# it in variables makes it both part of the source and part of the program
# itself, meaning it is accessible in python.

__license__ = "GNU Public License v3"

from ._legacy_aliases import LEGACY_NAMES as _LEGACY_NAMES, install as _install_legacy_aliases

_install_legacy_aliases()


def __getattr__(name: str):
    # `ampere.cite()`, how to cite ampere (W6.5), resolved on first use so that
    # `import ampere` alone still imports nothing but the alias finder.
    if name == "cite":
        from ._citation import cite

        return cite
    if name in _LEGACY_NAMES or name == "legacy":
        import importlib

        return importlib.import_module(f"ampere.{name}")
    raise AttributeError(f"module 'ampere' has no attribute {name!r}")
