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

import os
from datetime import date

# -- Project information -----------------------------------------------------

project = "AMPERE"
copyright = (
    f"{date.today().year}, Peter Scicluna, Francisca Kemper, Sundar"
    " Srinivasan, Jonathan Marshall, Oscar Morata, Alfonso Trejo, "
    "Sascha Zeegers, Lapo Fanciullo, Thavisha Dharmawardena"
)
author = (
    "Peter Scicluna, Francisca Kemper, Sundar Srinivasan, Jonathan"
    " Marshall, Oscar Morata, Alfonso Trejo, Sascha Zeegers, Lapo"
    " Fanciullo, Thavisha Dharmawardena"
)

from importlib.metadata import version

# 'pgmuvi' was another project's package name, left over from when this
# configuration was templated from it; it has never been installable here, so
# `pixi run docs` failed at configuration time regardless of content
# (W2.9's docs-build gate found this). The distribution actually installed
# for this documentation build is 'ampere-astro' (pyproject.toml's
# [project].name since W6.5, D13); the import name is still 'ampere'.
release = version("ampere-astro")
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
    "matplotlib.sphinxext.plot_directive",
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
#
# `sbi` (W5.32 (d)): the legacy `ampere.legacy.infer.sbi` module imports `sbi` at
# module level (on top of `torch`, already mocked above), and the `dev` docs
# environment does not carry the `sbi` extra -- adding it here would pull in
# torch too, for one page, which the two backend pages above already reject
# for the same reason. `sbi` is import-only in the docs build, exactly like
# `torch`/`jax` here: it is never installed in the environment this page is
# built in, so mocking it costs nothing and unblocks the autodoc import.
autodoc_mock_imports = [
    "torch",
    "pyro",
    "jax",
    "jaxlib",
    "numpyro",
    "equinox",
    "sbi",
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

# W6.1: the two v2 notebooks are executed at docs-build time; the legacy one
# never is. nbsphinx's default ('auto') executes any notebook that stores no
# outputs, and this repository stores none by policy (AGENTS.md ground rule
# 7 -- no run outputs or figures in git), so 'auto' means "execute every
# notebook here except one that opts out in its own metadata":
#
#   * notebooks/quickstart.ipynb and notebooks/Ampere_MBB_Example.ipynb are
#     written on v2, import nothing but ampere and the `dev` environment's own
#     packages, read the one data file they need (a tracked Spitzer spectrum
#     under examples/test_data) by a path found from the notebook's own
#     directory, and run in about a minute each. Their last cells state their
#     sampling budgets and wall times.
#   * notebooks/Embedding_nets.ipynb teaches the legacy `SBI_SNPE` embedding
#     vocabulary, imports torch and sbi -- neither is in the `dev`
#     environment the docs are built in, deliberately (pyproject.toml's
#     [tool.pixi.environments]) -- and trains an SNPE posterior on 10 000
#     simulations. It carries "nbsphinx": {"execute": "never"} in its own
#     metadata, so it is rendered from what it contains and never run.
#
# `ipykernel` is declared in the `dev` pixi feature alongside `pandoc`:
# without a registered `python3` kernelspec a notebook that is executed
# aborts the build with NoSuchKernel rather than a cell error.
nbsphinx_execute = "auto"

# A cell that raises fails the build. That is the point of executing the
# notebooks: a tutorial that no longer runs is a bug the docs build must see.
# (tests/examples/test_wstat_comparison.py's docstring used to cite the
# opposite setting as the reason the docs build was no gate for example code;
# for the two executed notebooks it now is one.)
nbsphinx_allow_errors = False

# nbsphinx's default is 30 s per cell; the slowest cell here is a zeus run of
# about half a minute, and a loaded machine is slower still.
nbsphinx_timeout = 300

# The kernel logs "Kernel is running over TCP without encryption" at WARNING
# on every start; it is about the local loopback connection the build opens,
# not about the notebooks, and it would otherwise add one WARNING line to the
# build log per executed notebook.
nbsphinx_execute_arguments = ["--IPKernelApp.log_level=ERROR"]

# imgmath_latex = "latex"

# W7.14: the figures of reading_the_diagnostics.rst are built at docs-build
# time by matplotlib's own plot directive, from the scripts under
# docs/source/plots/ -- this repository commits no binary artefacts
# (AGENTS.md ground rule 7), and a figure that is a script is a figure that
# cannot drift from the code it teaches. Each script runs the M2 study's (or
# persona_b2.py's) public calls at a reduced, stated budget and names the
# driver it reproduces in its docstring. The directive re-runs only a script
# whose own text changed, so a local rebuild after editing the page is fast;
# edit a shared helper (plots/_m2.py, plots/_irs.py) and `touch` the scripts
# that import it, or delete docs/_build/plot_directive.
# Build cost, `time pixi run --frozen docs` from a clean docs/_build:
# before (the same build without the page) 196 s; after, 281 s, 318 s and
# 341 s on three clean builds at a load average of about 3 on a 16-core
# machine (another agent's tests share it), so +85 to +145 s, under the
# three-minute cap; the nine figure scripts' CPU
# time is about 150 s (a shared fit is run once and reused across figures).
# The figure scripts' samplers are small-matrix, many-call workloads, where a
# 16-thread BLAS pool is slower than one thread (measured on this 16-core
# machine: the power-law fit took 100 s with the pool and 15 s without), so
# the build pins the pools before numpy loads. The notebooks' kernels inherit
# the setting.
for _pool in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_pool, "1")
plot_include_source = False
plot_html_show_source_link = False
plot_formats = [("png", 120)]
plot_html_show_formats = False
# The repository root, so `examples.m2_misspecification` and the IRS file's
# relative path resolve; conf.py lives in docs/source.
plot_working_directory = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
