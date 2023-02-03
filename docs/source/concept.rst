How does AMPERE work?
=====================



Diverse datasets
----------------

Astronomy increasingly faces the problem of combining multiple datasets to answer one question.
This results in large, complex datasets that are also rich in information.


Model mis-specification
-----------------------

As we build more powerful and precise instruments, our datasets also get richer.
This often results in features that our models are unable to explain (progress!) but this 'model mis-specification' leads to problems in inference;
unexplainable features create structure in the residuals, which in turn tend to distract the optimiser and distort parameter estimates.
For example, if your data is a spectrum that contains emission lines, but your model only predicts continuum, the posterior parameter estimates
will tend to slightly overestimate the continuum flux to compromise between line and continuum levels.

Flexible noise models
---------------------

This problem can be at least partially mitigated by realising that the structure residuals created by a mis-specified model are indistinguishable from correlated noise in the data.
Therefore, by explicitly modelling an additional component of correlated noise in the dataset and marginalising over it in our inference, we can add a certain amount of resistance
to mis-specification, and get less-biased posterior estimates of the parameters we are really interested in.

.. note::

   There is no *mis-specification proof* model here, just like there is no such thing as an earthquake-proof building. A big enough mis-specification
   will still break your inference, in this case by increasing the correlated noise component to arbitrary large levels so that the likelihood remains
   reasonable, just like a sufficiently energetic earthquake will bring down any building.


 
 AMPERE achieves this modelling the covariance matrix of the data as the sum of the identity matrix and an additional correlation matrix [#1]_, all multiplied by the variances of the data.
 The correlation matrix is composed from a set of *stationary* kernel functions, as in a Gaussian Process.
 This provides a high degree of flexibility in the types of noise that get generated, without adding large numbers of additional free parameters.

 A more detailed mathematical description will be added soon!


 .. rubric:: Footnotes
 ..  [#1] This effectively separates the noise into a correlated and an uncorrelated component.
