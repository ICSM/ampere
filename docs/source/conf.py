# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))

from datetime import date

# -- Project information -----------------------------------------------------

project = "AMPERE"
copyright = (
    f"{date.today().year}, Peter Scicluna, Francisca Kemper, Sundar"
    " Srinivasan, Jonathan Marshall, Sacha Hony, Sascha Zeegers, "
    "Lapo Fanciullo"
)
author = (
    "Peter Scicluna, Francisca Kemper, Sundar Srinivasan, Jonathan"
    " Marshall, Sacha Hony, Sascha Zeegers, Lapo Fanciullo"
)

from importlib.metadata import version

# 'pgmuvi' was another project's package name, left over from when this
# configuration was templated from it; it has never been installable here, so
# `pixi run docs` failed at configuration time regardless of content
# (W2.9's docs-build gate found this). The distribution actually installed
# for this documentation build is 'ampere' (pyproject.toml's [project].name).
release = version("ampere")
# for example take major/minor
version = ".".join(release.split(".")[:2])

master_doc = "index"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "nbsphinx",  # 'sphinx.ext.imgmath'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["pyphot*", "test*", "old*"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "alabaster"

# No custom static files. This entry was `['_static']`, a directory that has
# never existed in this repository, so every build emitted
# "html_static_path entry '_static' does not exist" (W2.15). Restore the entry
# together with the directory, when there is a stylesheet to put in it.
html_static_path = []


# -- autodoc: what is faked, and why (W2.15) ---------------------------------
#
# The documentation is built in the `dev` environment, which deliberately
# carries neither the `torch` nor the `jax` extra (pyproject.toml's
# [tool.pixi.environments]; `architecture.md` §4 rule 2 is why -- `dev` must
# never resolve either array library). A backend subpackage is the one place
# in ampere that imports its array library at module top level, so without
# mocks `ampere.backends.torch` and `ampere.backends.jax` cannot be imported
# here at all, and their API pages would be two autodoc "failed to import"
# warnings rather than documentation.
#
# The alternative -- rendering those two pages only under
# `pixi run -e torch docs` / `pixi run -e jax docs`, with a "needs the extra"
# note otherwise -- was rejected. No single environment has both libraries, so
# the *published* site, which CI builds in `dev`, would then carry no API
# reference at all for the two flagship Phase 2 backends. Mocking costs
# nothing here, because these pages document signatures and docstrings rather
# than runtime values.
#
# `bs4`, `requests`, `astropy` and `emcee` used to be on this list and are all
# gone. astropy and emcee are base dependencies and genuinely installed, and
# mocking an installed package replaces it wholesale -- which would have
# rendered `ampere.core`'s astropy.units-typed signatures as mock objects the
# moment those pages existed; bs4 and requests are imported nowhere in the
# package.
autodoc_mock_imports = [
    "torch",
    "pyro",
    "jax",
    "jaxlib",
    "numpyro",
    "equinox",
]


# Methods and attributes in the order their module declares them. The new
# namespaces' `__all__` lists are already alphabetical, so this only affects
# the inside of a class, where the order it was written in is the order it is
# meant to be read in.
autodoc_member_order = "bysource"

# Render a numpydoc `Attributes` section as `:ivar:` fields rather than as
# standalone `.. attribute::` directives (W2.15). Without this, every
# documented dataclass field and every property that is also listed in its
# class's Attributes section is described twice -- once by napoleon and once
# by autodoc -- which Sphinx reports as "duplicate object description of ...".
# The v2 namespaces document their dataclasses that way throughout, so this
# accounts for several dozen warnings on its own.
napoleon_use_ivar = True

# W2.11: notebooks are rendered from what they contain; the docs build never
# executes them. nbsphinx's default ('auto') executes any notebook with no
# stored outputs, and this repository stores none by policy (AGENTS.md ground
# rule 7 — no run outputs or figures in git), so 'auto' means "execute every
# notebook here except Ampere_MBB_Example.ipynb". Neither of the two it would
# reach can run, and both were part of why `pixi run docs` was red before this
# item:
#
#   * notebooks/quickstart.ipynb reads
#     'PGQuasars/PG1011-040/cassis_yaaar_spcfw_14191360t.fits' relative to the
#     *working directory*. No such file is in the repository (it is a Spitzer
#     CASSIS spectrum, i.e. exactly the kind of binary artefact that is not
#     committed), so the notebook cannot run from a clean clone anywhere. It
#     then fits it with a 100-walker emcee run and calls postProcess(), which
#     writes figures.
#   * notebooks/Embedding_nets.ipynb imports torch and sbi — neither is in the
#     `dev` environment the docs are built in, deliberately (pyproject.toml's
#     [tool.pixi.environments]) — and trains an SNPE posterior on 10 000
#     simulations.
#
# So this is the "explicitly excluded, with a note" half of that repair rather
# than a claim that they work. The gate that actually proves ampere's shipped
# example code runs is tests/examples (W2.9), which is in `pixi run test-all`
# as of W2.11 — a real pytest run, not a docs build that
# `nbsphinx_allow_errors` would let pass regardless.
#
# `ipykernel` is nonetheless declared in the `dev` pixi feature alongside
# `pandoc`: without a registered `python3` kernelspec the failure mode for any
# notebook that *is* executed is a build-aborting NoSuchKernel rather than a
# cell error, and this line is one word away from being turned back on for a
# notebook that earns it.
nbsphinx_execute = "never"

# Belt and braces with the line above: a cell that raises is reported, not
# fatal. tests/examples/test_wstat_comparison.py's docstring cites this
# setting as the reason the docs build is not a trustworthy gate for example
# code — that reasoning stands.
nbsphinx_allow_errors = True

# imgmath_latex = "latex"
