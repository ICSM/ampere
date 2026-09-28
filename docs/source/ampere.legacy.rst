ampere.legacy package (legacy, frozen)
=======================================

.. warning::

   **This is legacy ampere — the v1 code, frozen in place.** It still runs,
   the characterisation suite exists to keep it running, and the published
   science came out of it. Nothing in ``ampere.legacy.data``,
   ``ampere.legacy.models`` or ``ampere.legacy.infer`` is being changed, and
   none of it implements the v2 contracts: a legacy model is not an
   :class:`ampere.core.Model`, a legacy search is not an
   :class:`ampere.inference.Engine`, and the two halves do not interoperate.

   The old top-level names (``ampere.data`` and the rest) still import, as
   aliases of these modules — see :doc:`legacy`.

   New work belongs on the v2 API — :mod:`ampere.core`,
   :mod:`ampere.backends`, :mod:`ampere.inference`, :mod:`ampere.results`.
   See :doc:`overview`.

   Docstrings here pre-date the project's numpydoc conventions, so this half
   of the documentation renders with warnings that are deliberately left
   alone rather than fixed in place.

.. automodule:: ampere.legacy
    :members:
    :undoc-members:
    :show-inheritance:
    :no-index:

Subpackages
-----------

.. toctree::

    ampere.legacy.data
    ampere.legacy.infer
    ampere.legacy.models
    ampere.legacy.utils
