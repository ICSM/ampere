"""
ampere.backends.legacy
======================

Thin re-exports of the existing numpy-based inference backends under a stable
``ampere.backends`` import path.

These classes are unchanged from ``ampere.infer``; this module simply provides
a cleaner public API for code that should not depend on the internal layout of
``ampere.infer``.

Classes
-------
EmceeSearch
    Affine-invariant ensemble sampler (emcee).
DynestySearch
    Dynamic nested sampling (dynesty).
ZeusSearch
    Ensemble slice sampler (zeus-mcmc).
NestedSearch
    Generic nested sampling wrapper.
"""

from ampere.infer.emceesearch import EmceeSearch
from ampere.infer.dynestysearch import DynestySearch
from ampere.infer.zeussearch import ZeusSearch
from ampere.infer.nestedsearch import NestedSearch

__all__ = ["EmceeSearch", "DynestySearch", "ZeusSearch", "NestedSearch"]
