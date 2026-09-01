"""
ampere.core
===========

Stable base abstractions shared by all ampere sub-packages.

This package must not import from ``ampere.models``, ``ampere.infer``, or
``ampere.data``.  Everything else may depend on it.

Contents
--------
parameter
    :class:`Parameter` and :class:`ParameterSet` – named, prior-equipped
    parameter containers that replace manual ``npars`` / positional ``theta``
    array slicing.
model
    :class:`Model` – abstract base class for all physical models.
data
    :class:`Data` – abstract base class for all data objects.
"""

from .parameter import Parameter, ParameterSet
from .model import Model
from .data import Data

__all__ = ["Parameter", "ParameterSet", "Model", "Data"]
