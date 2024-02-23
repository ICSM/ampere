""" This module defines the base class for all `ampere.jax` likelihoods. 

The `Likelihood` class is a subclass of `Module`, and is intended to be used as
a base class for all likelihoods in the `ampere.jax` namespace. It provides a
`forward` method that should be implemented in a subclass, and a `log_prob`
method that should be implemented in a subclass.

This module also provides a few of the most important astrophyiscal likelihoods
as subclasses of `Likelihood`, including `GaussianLikelihood` and
`PoissonLikelihood`. These are intended to be used in data objects, and to be
extended to include more complex likelihoods as needed.
"""