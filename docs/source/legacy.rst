The legacy API: kept, frozen, and where it lives
=================================================

Ampere's first API — the ``Spectrum``/``Photometry`` data classes, the
``Model`` base class and its subclasses, the ``ampere.infer`` searches and
their post-processors — produced the published science, and it stays.
This page is the policy: what "kept" and "frozen" mean, where the code now
lives, and the one rule that applies to v2's own deprecations instead.

Where it lives
--------------

Since W6.0 the legacy code is one subpackage, :mod:`ampere.legacy`:
``ampere.legacy.data``, ``ampere.legacy.models``, ``ampere.legacy.infer``,
``ampere.legacy.utils`` and ``ampere.legacy.logger``. A legacy script gets
its whole old surface under one name with a single substitution::

    from ampere import legacy as ampere

    spec = ampere.data.Spectrum(...)        # was ampere.data.Spectrum
    search = ampere.infer.emceesearch.EmceeSearch(...)

**The old top-level names still import.** ``import ampere.data``,
``from ampere.models import Model``, ``from ampere.infer.sbi import SBI_SNPE``
and ``ampere.utils.pyphot_compat`` all keep working, with no warning and no
removal date: they are aliases, resolved on first import to the module under
:mod:`ampere.legacy`. The alias and the target are the *same* module object,
so ``ampere.data is ampere.legacy.data`` and a class's ``__module__``,
its repr and any pickle of it say ``ampere.legacy.data``. Prefer the
``ampere.legacy`` spelling in new writing; the old names are kept for the
scripts that already exist.

``import ampere`` itself no longer loads any of this. It binds the version
and licence attributes and installs the aliases; the legacy modules and
their dependencies are imported the first time one of them is named. v2's
namespaces — :mod:`ampere.core`, :mod:`ampere.backends`,
:mod:`ampere.inference`, :mod:`ampere.results` — are likewise imported
only when asked for.

What "kept" means
------------------

* **Indefinitely.** There is no removal scheduled, no deprecation warning,
  and no date after which a legacy import stops working. Ruled at Phase 6
  (D1 (b), 2026-09-28) over the alternatives of a deprecation interval or
  a removal at the beta: the published science was produced with this
  code, and every script that reproduces it should keep running.
* **Documented as legacy.** The reference stays on the site
  (:doc:`ampere.legacy`), under its own heading, with the banner that says
  what it is. Its docstrings pre-date the project's conventions and render
  with warnings that are left alone rather than fixed in place.
* **Tested.** The characterisation suite (``tests/characterisation``, the
  ``test-characterisation`` pixi task) runs the minimal working examples
  end to end on the frozen code, and is what "legacy still works" means.
  It runs on every release.

What "frozen" means
--------------------

* **Unchanged.** No new features, no refactors, no renames, no style or
  typing passes. The only edits are the mechanical relocation above and
  bug fixes that are **a line** — never a behaviour change, and each with
  a characterisation row that would have caught the bug.
* **It implements none of the v2 contracts.** A legacy ``Model`` is not an
  :class:`ampere.core.Model`; a legacy search is not an
  :class:`ampere.inference.Engine`; a legacy ``Photometry`` is not a
  :class:`~ampere.core.Dataset`. The two halves do not interoperate, and a
  fit is built entirely on one side or the other. :doc:`migrating` maps
  each legacy piece to the v2 piece that plays its role.
* **Its dependencies are its own.** A package only the legacy code needs
  is installed by the extra that names it and imported only when the
  legacy module that uses it is; a missing optional dependency never blocks
  ``import ampere`` or any v2 namespace (``tests/test_imports.py`` holds
  the package to this).

v2's own deprecations are a separate rule
-------------------------------------------

The policy above is about the legacy code. Names that **v2** deprecates —
a helper renamed, a keyword replaced — follow the ordinary rule instead:
the old name keeps working through a :class:`DeprecationWarning` that
names the replacement and the version the old name is removed in, and it
is removed at **1.0.0 final**. Deprecations made during the beta are
therefore gone at 1.0.0; nothing deprecated after 1.0.0 is removed before
the next major version. The current list:

* :func:`ampere.core.regularised_horseshoe` → :func:`ampere.core.shrinkage_horseshoe`
  (renamed at W5.27; removed at 1.0.0 final).

The examples
------------

The legacy examples under ``examples/`` — the minimal working examples and
their engine variants, ``NGC6302``, the modified blackbody, the PHOENIX
star, the carbon star, the star-plus-disc model — stay exactly as they
are and keep running on the legacy code; the minimal working examples are
the characterisation suite's anchors. Each gains a **v2 twin** beside it
(the same model, data and question written against :mod:`ampere.core`),
and the migration guide links the pairs as its side-by-side material.
