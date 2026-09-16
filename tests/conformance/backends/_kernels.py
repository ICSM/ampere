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

from ampere.core import Kernel, Product, Sum, WarpedKernel

from ..protocol import CovarianceSpec, KernelFamily

__all__ = [
    "USER_FAMILY",
    "WARP_AMPLITUDE_KNOTS",
    "WARP_AMPLITUDE_LEVELS",
    "WARP_INCREMENTS",
    "WARP_INPUT_KNOTS",
    "build_kernel",
    "user_kernel_type",
    "warped_kernel",
]

#: The warp of the W5.7 rows, fixed here so every fixture declares the *same*
#: model and the rows can be compared with each other. The input knots span
#: ``composition.GP_GRID`` (1.0 to 11.5) and the increments are deliberately
#: not all equal — a warp with equal slopes is a rescaling of the coordinate,
#: which a Matérn's own ``length_scale`` already expresses, so it would be a
#: weak row.
WARP_INPUT_KNOTS: tuple[float, ...] = (1.0, 4.5, 8.0, 11.5)
#: Fixed values for the three input-warp knot variables, held rather than
#: sampled: the battery's GP rows compare numbers against ``scipy`` and a
#: fitted warp would only add sampling dimensions they do not use.
WARP_INCREMENTS: tuple[float, ...] = (0.7, -0.5, 0.3)
#: The amplitude warp's knot locations and its fixed log-amplitude levels.
WARP_AMPLITUDE_KNOTS: tuple[float, ...] = (1.0, 6.25, 11.5)
WARP_AMPLITUDE_LEVELS: tuple[float, ...] = (0.25, -0.35, 0.4)

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
    if spec.family is KernelFamily.WARPED:
        return warped_kernel(build_kernel(spec.terms[0], classes), axes=spec.axes)
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


def warped_kernel(
    base: Kernel,
    *,
    input_warp: bool = True,
    amplitude_warp: bool = True,
    identity: bool = False,
    axes: tuple[str, ...] | None = None,
) -> WarpedKernel:
    """The battery's warp over *base*, with every knot variable held fixed.

    One place, so that the conformance rows, their oracles and the
    identity-warp row all mean the same warp. ``identity=True`` sets every knot
    variable to zero, which is the point at which a :class:`WarpedKernel`'s
    covariance must be **bit-identical** to its base's, not merely close to it.
    """
    keywords: dict[str, Any] = {"non_centred": True}
    if input_warp:
        keywords["input_warp"] = WARP_INPUT_KNOTS
        keywords["increments"] = 0.0 if identity else WARP_INCREMENTS
        keywords["input_scale"] = 1.0
    if amplitude_warp:
        keywords["amplitude_warp"] = WARP_AMPLITUDE_KNOTS
        keywords["levels"] = 0.0 if identity else WARP_AMPLITUDE_LEVELS
        keywords["amplitude_scale"] = 1.0
    if axes is not None:
        keywords["axes"] = axes
    return WarpedKernel(base, **keywords)
