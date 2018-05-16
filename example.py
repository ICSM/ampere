import numpy as np
from ampere.data import Spectrum, Photometry
from ampere.emceesearch import EmceeSearch
from ampere.PowerLawAGN import PowerLawAGN
from ampere.extinction import CCMExtinctionLaw
import corner
import matplotlib.pyplot as plt


if __name__=="__main__":
    """
    This script is a demonstration of how ampere (should) work(s). 
    """


    """ Read in some data """
    irs = Spectrum.fromArray(filename)
    #If we have more than one dataset we need multiple data objects

    """ Define the model """
    model = PowerLawAGN(irs.wavelengths) #We should probably define a wavelength grid separate from the data wavelengths for a number of reasons, primarily because the photometry needs an oversampled spectrum to compute synthetic photometry
    
    """ and the extinction law (if necessary) (could be moved inside model as only some models will need this) """
    ext = CCMExtinctionLaw()

    model.extinction = ext #combination of multiple different models is something we need to consider how to implement. 

    """ Connect the dots: """
    """ Hook it up to an optimiser """
    opt = EmceeSearch(model = model, data = [irs]) #will need to either define some info about the parameters here or write some code to determine how the model works by introspection

    """ if you want, overload the default priors with a new function """
    def lnprior(self, inputs, **kwargs):
        return 0 #flat prior
    
    model.lnprior = lnprior
    
    """ Run it """
    pos = [
           [
               0 + np.random.randn() for i in range(model.npars)
           ]
           for j in range(nwalkers)
          ]
    opt.optimise(nsamples = 10000, nwalkers = 100,burnin=1000,pos=pos)

    """save optimiser state for later """
    opt.dump(filename="output_file",pickle=True) #Save as a python object
    
    """ Sort out the results """
    fig = corner.corner(opt.samples,labels=opt.labels)
    plt.show()

    """ Figure out what it all means """