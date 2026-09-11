"""One recursive kernel builder, shared by every fixture (W4.5).

A :class:`~tests.conformance.protocol.CovarianceSpec` became a *tree* at W4.5,
and four fixtures translating that tree independently would have been four
places for a row to be subtly different. They are not different: each fixture
hands this module its own class table and nothing else, which is the same
arrangement ``backends/__init__.py`` makes for the fixtures themselves.

The user term is here for a sharper reason. The acceptance criterion is that
**one** registration of a user kernel's celerite representation carries it onto
the O(N) path on all three backends, and the registry is keyed on the kernel
*family* — so every backend's user kernel must declare the **same** family name
while subclassing its **own** ``Matern12``. :func:`user_kernel_type` builds
exactly that, one class per backend class, all sharing
:data:`USER_FAMILY`.
"""

from __future__ import annotations

import functools
from typing import Any

from ampere.core import Kernel, Product, Sum

from ..protocol import CovarianceSpec, KernelFamily

__all__ = ["USER_FAMILY", "build_kernel", "user_kernel_type"]

#: The family name every backend's user kernel declares. One name, so one
#: ``register_quasiseparable_term`` call reaches all three O(N) paths — which
#: is the property the registry row exists to demonstrate.
USER_FAMILY = "conformance_user_term"


@functools.cache
def user_kernel_type(matern12: type[Kernel]) -> type[Kernel]:
    """A user kernel: *matern12* under a family name ampere does not know.

    Deliberately a re-labelling rather than new mathematics. What the row
    proves is the *route* — that a kernel declared outside ampere reaches
    ``QuasisepGP`` once its representation is registered — and a term whose
    algebra also had to be checked would make a failure ambiguous.
    """

    return type(
        f"UserTerm{matern12.__module__.rsplit('.', 2)[-2].title()}",
        (matern12,),
        {
            "FAMILY": USER_FAMILY,
            "__doc__": (
                f"A conformance user kernel over {matern12.__module__}.{matern12.__name__}."
            ),
        },
    )


def build_kernel(spec: CovarianceSpec, classes: dict[KernelFamily, type[Kernel]]) -> Kernel:
    """The kernel *spec* declares, from one backend's classes.

    ``classes`` maps each leaf family this backend supplies to its class;
    :class:`~ampere.core.Sum` and :class:`~ampere.core.Product` are taken from
    ``ampere.core`` unconditionally, because a composite has no arithmetic of
    its own — it combines its children's matrices and adopts their namespace,
    device and capability flags.
    """
    if spec.family in (KernelFamily.SUM, KernelFamily.SPECTRAL_MIXTURE):
        return Sum(*(build_kernel(term, classes) for term in spec.terms))
    if spec.family is KernelFamily.PRODUCT:
        return Product(*(build_kernel(term, classes) for term in spec.terms))
    keywords: dict[str, Any] = {}
    if spec.axes is not None:
        keywords["axes"] = spec.axes
    if spec.family is KernelFamily.SHO:
        return classes[KernelFamily.SHO](spec.amplitude, spec.period, spec.quality, **keywords)
    if spec.family is KernelFamily.ROTATION:
        return classes[KernelFamily.ROTATION](
            spec.amplitude,
            spec.period,
            spec.quality,
            spec.delta_quality,
            spec.fraction,
            **keywords,
        )
    if spec.length_scale_unit is not None:
        keywords["length_scale_unit"] = spec.length_scale_unit
    if spec.family is KernelFamily.USER:
        return user_kernel_type(classes[KernelFamily.MATERN12])(
            spec.amplitude, spec.length_scale, **keywords
        )
    return classes[spec.family](spec.amplitude, spec.length_scale, **keywords)
