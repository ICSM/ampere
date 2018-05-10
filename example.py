import numpy as np
from ampere.data import Spectrum, Photometry
from ampere.emceesearch import EmceeSearch
from ampere.ClassRT import RadTransWrap
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
    model = RadTransWrap(irs.wavelengths)
    """ and the extinction law (if necessary) (could be moved inside model as only some models will need this) """
    ext = CCMExtinctionLaw()

    model.extinction = ext

    """ Connect the dots: """
    """ Hook it up to an optimiser """
    opt = EmceeSearch(model = model, data = irs)

    """ if you want, overload the default priors with a new function """
    def lnprior(self, inputs, **kwargs):
        return 0 #flat prior
    
    opt.lnprior = lnprior
    
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
