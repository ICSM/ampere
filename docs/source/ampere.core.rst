ampere.core
===========

The backend-neutral contracts, and everything shared between backends: the
parameter and prior model, the container schema, transformations and
instruments, likelihood families, noise models, kernels and GP solvers, the
dataset and fitting-problem surface, batched simulation with its executor
protocol, the coordinate-value-mask encoding embedding networks consume, the
lowering registry and the realisation registry.

Nothing here imports torch, jax or any other optional dependency, lazily or
otherwise — that is what makes the three backends interchangeable consumers
of one vocabulary. The one exception in spirit but not in policy is
:class:`~ampere.core.QuasisepGP`, whose celerite2 dependency is a *base*
dependency of the whole distribution and is imported lazily on first use.

.. automodule:: ampere.core
   :members:
   :imported-members:
   :undoc-members:
   :show-inheritance:

Constants
---------

Every name below is re-exported from ``ampere.core`` — ``from ampere.core
import SEPARATOR`` — and is shown here under the module that defines it,
which is where its documentation lives.

.. autodata:: ampere.core.parameter.SEPARATOR
.. autodata:: ampere.core.results_schema.DEFAULT_CHANNEL
.. autodata:: ampere.core.results_schema.COORDINATE_RTOL
.. autodata:: ampere.core.results_schema.REGULARITY_RTOL
.. autodata:: ampere.core.dataset.MODEL_COMPONENT
.. autodata:: ampere.core.dataset.INSTRUMENT_COMPONENT
.. autodata:: ampere.core.dataset.LIKELIHOOD_COMPONENT
.. autodata:: ampere.core.dataset.LATENT_COMPONENT
.. autodata:: ampere.core.dataset.SHARED_COMPONENT
.. autodata:: ampere.core.dataset.LATENT_NAME
.. autodata:: ampere.core.dataset.DEFAULT_FAILURE_HISTORY
.. autodata:: ampere.core.kernels.NUMPY_OPS
.. autodata:: ampere.core.kernels.TermBuilder
.. autodata:: ampere.core.realisation.RealisationFactory
.. autodata:: ampere.core.simulate.ChunkHook
.. autodata:: ampere.core.encoding.ENCODING_VERSION
.. autodata:: ampere.core.encoding.DEFAULT_FOURIER_BANDS
.. autodata:: ampere.core.encoding.COLUMN_GROUPS
.. autodata:: ampere.core.encoding.SET_KIND
.. autodata:: ampere.core.encoding.FLAT_KIND

