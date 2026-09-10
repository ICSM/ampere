.. AMPERE documentation master file, created by
   sphinx-quickstart on Mon Mar  1 09:58:54 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to AMPERE's documentation!
==================================

**AMPERE** is a Bayesian fitting environment for astronomers. It exists to
make it possible to model complex, heterogeneous datasets — spectra and
photometry together, and more besides — *even when your model cannot explain
everything in the data*.

That last clause is the point of the package. A model deficiency shows up as
structure in the residuals, and structure in the residuals is
indistinguishable from correlated noise. So ampere gives each dataset a
flexible likelihood — a Gaussian process over the residuals, marginalised
over while the physical parameters are fitted — which absorbs what the model
cannot explain instead of letting it bias the answer. The
:doc:`m2_misspecification` page is the measurement of what that buys and what
it costs: on a deliberately misspecified 20 000-point spectrum, an ordinary
chi-square fit reports a parameter **113 posterior standard deviations** away
from the truth, while the flexible likelihood stays within 0.6 and keeps the
truth inside its 68 % interval.

Ampere is in its **v2** redesign. The current release is a backend-neutral
core of frozen contracts with three backends implementing it — pure
numpy/scipy, torch, and jax — and six inference engines written once against
the contracts and run on any of them, the sixth being
:class:`~ampere.inference.SBIEngine` for models with no likelihood to write
down at all (see :doc:`sbi`). :doc:`overview` is the map.

User guide
----------

.. toctree::
   :maxdepth: 2

   install
   overview
   concept
   tutorials
   advanced
   faqs


API reference
-------------

.. toctree::
   :maxdepth: 2

   api


Contributing
------------

We very much welcome contributions to AMPERE! Please take a look at our `github repository <https://github.com/ICSM/ampere/>`_ for more information on how to contribute!

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
