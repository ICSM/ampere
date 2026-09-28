"""Legacy ampere (v1), kept and frozen — see the legacy page of the docs.

``from ampere import legacy as ampere`` gives a legacy script its whole old
surface under one name. The old top-level names (``ampere.data`` and the
rest) remain importable as aliases of these modules.
"""

from . import data, infer, logger, models, utils

__all__ = ["data", "infer", "logger", "models", "utils"]
