# ampere
Ampere is an attempt to produce a fitting environment that can natively handle multiple different kinds of astronomical data with differing information content, even when the model being applied to them might be missing some of the key processes and be unable to actually reproduce all aspects of the observations properly.
This is motivated by the need within out team for a tool to simultaneously fit the SEDs and spectra of dusty objects to constrain, among other things, the dust properties and mineralogy. 
However, the final product will be more general, such that it can be applied to a wide range of astronomical questions.

To achieve this, we include a simple parametric model of correlated noise for each dataset which is marginalised over when fitting the parameters of the models. 
This approach was chosen because, typically, deficiencies in the model result in structured residuals, and structured residuals are equivalent to correlated noise.
This effectively downweights parts of the data which the model doesn't represent well without having to manually identify these regions.

At present, ampere is in the alpha testing phase, but we anticipate a beta release in the near future. If you are interested, please get in touch with us!


## Installation

git clone repository and install with pip (ideally in a new conda environment):

> git clone git@github.com:ICSM/ampere.git
> 
> cd ampere
> 
> pip install -e .

This will install the minimal package, but leave out some extra features. For all features, instead do

> pip install -e .[all]
