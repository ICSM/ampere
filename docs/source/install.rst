Installing AMPERE
=================

Ampere runs on **Python 3.11, 3.12 and 3.13**; those are the versions CI
tests. It is not yet published on PyPI, so every route below starts from a
clone.

.. warning::

   There is an unrelated project called ``ampere`` on PyPI (a battery
   modelling package). ``pip install ampere`` will install **that**, not this.
   Install from the repository.

.. code-block:: console

    $ git clone https://github.com/ICSM/ampere.git
    $ cd ampere

With pixi (recommended)
-----------------------

`pixi <https://pixi.sh>`_ is the supported route: it builds the whole
environment — the interpreter included — from ``pyproject.toml`` and the
committed lock file, so everyone gets the same package set. With only pixi
installed:

.. code-block:: console

    $ pixi install -e dev

That is the daily-use environment: Python 3.13, ampere as an editable
install, the developer tooling, and the ``zeus`` and ``extinction`` extras.
Then run anything through ``pixi run``:

.. code-block:: console

    $ pixi run test-phase1          # the core, results and conformance suites
    $ pixi run test-all             # everything the v2 namespaces own
    $ pixi run lint                 # ruff
    $ pixi run typecheck            # pyrefly, scoped to the v2 namespaces
    $ pixi run docs                 # build this documentation

A plain ``pixi run <task>`` is the same as ``pixi run -e dev <task>``. Other
environments are selected with ``-e``:

.. code-block:: console

    $ pixi run -e torch test-all    # the torch backend and its suites
    $ pixi run -e jax test-all      # the jax backend and its suites
    $ pixi run -e sbi  test-characterisation

``pixi.lock`` is committed; ``.pixi/`` is not.

With pip
--------

If you would rather manage the environment yourself, use a clean virtual
environment (conda, venv, uv — whichever you prefer) and install the clone:

.. code-block:: console

    $ python -m venv .venv && source .venv/bin/activate
    $ pip install -e .

That gives you the **base install**, which is a complete fitting environment
rather than a stub: the reference backend, the whole of :mod:`ampere.core`
including the flexible GP likelihood and its exact O(N) solver, the
gradient-free engines, and :mod:`ampere.results`.

Add extras in the usual way:

.. code-block:: console

    $ pip install -e ".[torch]"
    $ pip install -e ".[jax]"
    $ pip install -e ".[all]"

Extras
------

.. list-table::
   :header-rows: 1
   :widths: 16 34 50

   * - Extra
     - Adds
     - Unlocks

   * - *(none)*
     - numpy, scipy, astropy, matplotlib, spectres, tqdm, pyphot, emcee,
       dynesty, corner, **celerite2**, **arviz**, **h5netcdf**
     - The base install. The reference backend, the flexible likelihood, the
       O(N) solver, emcee and dynesty, and the ArviZ results format. This must
       always work; a CI job installs exactly this and imports the package.

   * - ``torch``
     - torch, pyro-ppl
     - :mod:`ampere.backends.torch`, and therefore
       :class:`~ampere.inference.NUTSEngine` and
       :class:`~ampere.inference.VIEngine` on a torch problem

   * - ``jax``
     - jax, numpyro, equinox
     - :mod:`ampere.backends.jax`, and the same two engines on a jax problem

   * - ``zeus``
     - zeus-mcmc
     - :class:`~ampere.inference.ZeusEngine`

   * - ``sbi``
     - torch, sbi
     - Neural posterior estimation, in the **legacy** ``ampere.infer.sbi``

   * - ``extinction``
     - dust_extinction
     - Legacy ``ampere.models.extinctionModels.F99Extinction``

   * - ``dev``
     - pytest, pytest-cov, pytest-benchmark, coverage, ruff, pyrefly, sphinx,
       nbsphinx, nbconvert
     - Contributor tooling. Building this documentation also needs ``pandoc``
       and ``ipykernel``, which are not Python packages — the pixi ``dev``
       environment installs them for you.

   * - ``all``
     - every feature extra above (not ``dev``)
     - —

Two notes on the table. **celerite2 is a base dependency, not an extra**: the
reference backend's whole point is being a complete, scalable environment
with no heavy dependencies, and a base install that still had the O(N³)
problem would defeat it. **arviz and h5netcdf are base dependencies too**, as
of the engine drivers — a base install that could build a run but not save it
would be a base install that cannot sample. There is deliberately no
``arviz`` extra any more, not even as an empty alias.

Nothing outside ``ampere/backends/torch`` and ``ampere/backends/jax`` imports
torch or jax, at any depth, even lazily. Importing a backend is your explicit
opt-in to its dependency, and a missing optional dependency raises
:class:`~ampere.core.OptionalDependencyError` on *use*, naming the package and
the extra, rather than a bare ``ModuleNotFoundError`` three frames down.

Checking it worked
------------------

.. code-block:: console

    $ python -c "import ampere; from ampere.core import FittingProblem; print(ampere.__version__)"

and, with the ``dev`` tooling present:

.. code-block:: console

    $ pixi run test-all
