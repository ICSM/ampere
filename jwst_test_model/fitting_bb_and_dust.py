import sys
sys.path.insert(1, '/home/zeegers/git_ampere/ampere/')
import numpy as np
import os
import ampere
from ampere.data import Spectrum, Photometry
from ampere.emceesearch import EmceeSearch
from ampere.models import Model
from spectres import spectres
import pyphot
from emcee import moves
from astropy.modeling import models
from astropy import units as u
from astropy.modeling.models import BlackBody
import matplotlib.pyplot as plt
from scipy import interpolate


#First we will define a rather simple model

class Blackbody_dust(Model):
    '''This is a blackbody model based on the very simple model from Peter 

    We will first try to fit a blackbody model to a dataset with dust 

    '''
    def __init__(self, wavelengths, flatprior=True,
                 lims=np.array([[10., 300000.],
                                [1.,6.]])):
        '''The model constructor, which will set everything up

        This method does essential setup actions, primarily things that 
        may change from one fit to another, but will stay constant throughout 
        the fit. This may be things like the grid of wavelengths to calculate
        the model output on, or establishing the dust opacities if involved.
        There are also several important variables it *MUST* define here
        '''
        self.wavelength = wavelengths
        
        # Getting the opacities from the folder 
        #opacityDirectory = os.path.dirname(__file__)+'/Opacities/'
        opacityDirectory = os.path.dirname(os.path.realpath('__file__'))+'/Opacities/'
        print("Directory:", opacityDirectory)
        opacityFileList = os.listdir(opacityDirectory)
        opacityFileList = np.array(opacityFileList)[['sub.q' in zio for zio in opacityFileList]] # Only files ending in sub.q are valid (for now). At the moment there are 6 files that meet this criteria
        nSpecies = opacityFileList.__len__()
        opacity_array = np.zeros((wavelengths.__len__(), nSpecies))
        
        
        for j in range(nSpecies):
            tempData = np.loadtxt(opacityDirectory + opacityFileList[j], comments = '#')
            print(opacityFileList[j])
            tempWl = tempData[:, 0]
            tempOpac = tempData[:, 2]            
            

            f = interpolate.interp1d(tempWl, tempOpac, assume_sorted = False)
            opacity_array[:,j] = f(wavelengths)#wavelengths)
            
        self.opacity_array = opacity_array
        self.nSpecies = nSpecies
        
        self.npars = nSpecies + 2 #Number of free parameters for the model (__call__()). For some models this can be determined through introspection, but it is still strongly recommended to define this explicitly here. Introspection will only be attempted if self.npars is not defined.
        self.npars_ptform = nSpecies + 2 #Sometimes the number of free parameters is different when using the prior transform instead of the prior. In that case, self.npars_ptform should also be defined.
        #You can do any other set up you need in this method.
        #For example, we could define some cases to set up different priors
        #But that's for a slightly more complex example.
        #Here we'll just use a simple flat prior
        self.lims = lims
        self.flatprior = flatprior

    def __call__(self, temp, radius_sol, *args, **kwargs):
        '''The model itself, using the callable class functionality of python.

        This is an essential method. It should do any steps required to 
        calculate the output fluxes. Once done, it should stop the output fluxes
        in self.modelFlux.
        '''
        dustAbundances = np.array(args) # instead of 10**np.array(args)
        fModel = (np.matmul(self.opacity_array, dustAbundances))
        fModel2 = np.exp(-fModel)
        
        #plt.plot(wavelengths,self.opacity_array[:,0], label = "0")
        #plt.plot(wavelengths,self.opacity_array[:,1], label = "1")
        #plt.plot(wavelengths,self.opacity_array[:,2], label = "2")
        #plt.plot(wavelengths,self.opacity_array[:,3], label = "3")
        #plt.plot(wavelengths,self.opacity_array[:,4], label = "4")
        #plt.plot(wavelengths,self.opacity_array[:,5], label = "5")
        #plt.legend()
        #plt.show()

        wavelengths_aa = (self.wavelength*u.micron).to(u.AA)
        bb =  BlackBody(temperature=temp*u.K)
        pc  = 3.086e16 # m
        rsol = 696340e3 # m
        
        distance_pc=1100.
        
        Rstar = radius_sol*rsol
        r = distance_pc*pc
        
        flux = bb(wavelengths_aa)
        flux_mjy = flux.to(u.mJy / u.sr)*(Rstar/r)**2.*fModel2
        
        # here we need to put dust models ...
        dustAbundances = np.array(args) # instead of 10**np.array(args)
        
        self.modelFlux = flux_mjy #slope*self.wavelength + intercept        

    def lnprior(self, theta, **kwargs):
        temp = theta[0]
        #distance_pc  = theta[1]
        radius_sol= theta[1]
        if self.flatprior:
            if np.all(theta[2:] > 0.) and (self.lims[0,0] < theta[0] < self.lims[0,1]) and (self.lims[1,0] < theta[1] < self.lims[1,1]):
                return 0
            else:
                return -np.inf
        else:
            raise NotImplementedError()
        

    def prior_transform(self, u, **kwargs):
        '''The prior transform, which takes samples from the Uniform(0,1) distribution to the desired distribution.

        This is only included for completeness and to demonstrate how a prior 
        transform function should look. This example only uses emcee for 
        fitting, which uses the lnprior function instead. Prior transforms are 
        required by nested-sampling codes and similar approaches.
        '''
        if self.flatprior:
            theta = np.zeros_like(u)
            return (lims[:,1] - lims[:,0]) * u - lims[:,0]
        else:
            raise NotImplementedError()


if __name__ == "__main__":
    """ Set up the inputs for the model """
    """ wavelength grid """
    wavelengths = np.linspace(2.0,28, 900)
    
    # constants and definition of R factor
    pc  = 3.086e16 # m
    rsol = 696340e3 # m

    """ Choose some model parameters """
    temp = 18000.#Keep it super simple for now
    radius_sol = 5.0
    distance_pc = 1100. # perhaps you should fix the distance!!

    #Now init the model:
    model = Blackbody_dust(wavelengths)
    #And call it to produce the fluxes for our chosen parameters
    model(temp, radius_sol, 1.0, 0.1, 0.01, 0.01, 0.01, 0.01)
    model_flux = model.modelFlux
    
    # Here starts the data part 
    
    keywords = {'waveCol': 'wavelength', #replace with the name of the column with the wavelength information
                'fluxCol': 'flux',
                'waveUnit': 'um',
                'fluxUnit': 'Jy',
                'uncsCol': 'uncertainty'}    
    
    dataDir = os.getcwd()
    specFile = 'blackbody_jwst.txt'
    #photFile = 'vizier_votable_pruned_no2MASS.vot'
    irs = Spectrum.fromFile(dataDir+'/'+specFile,format='User-Defined', filetype='text', keywords=keywords)
    #irs = [irs.selectWaves(low = 1., up=18.) for wavelength in irs] # alternative to next two lines, doesn't work
    irs[0].selectWaves(low = 5.0, up = 18.) #following Srinivasan et al. 2017
    #irs[1].selectWaves(low = 0.6, up = 18.) #following Srinivasan et al. 2017
    
    print(irs[0])   
    
    dataSet = irs
    
    fig = plt.figure()    
    ax = fig.add_subplot(111)
    
    for i in irs:
        ax.plot(i.wavelength, i.value, '-',color='red')
    
    
    ax.plot(wavelengths, model_flux)
    plt.show()
    
    #Ampere exposes acces to emcee's moves interface. This can be useful if the posterior turns out to not be well behaved - the default move only deals well with posteriors that are monomodal and approximately Gaussian. Here's an example that usually deals a bit better with posteriors that don't meet these criteria:
    m = [(moves.DEMove(), 0.8),
        (moves.DESnookerMove(), 0.2),
         ]

    #Now we set up the optimizer object:
    optimizer = EmceeSearch(model=model, data=dataSet, nwalkers=100, moves=m)
    
    optimizer.optimise(nsamples = 2000, burnin=1500, guess=[
        [20000., 5.0, 1.0, 0.1, 0.01, 0.01, 0.01, 0.01, #The parameters of the model
         #1.0, 0.1, 0.1, #Each Spectrum object contains a noise model with three free parameters
         #The first one is a calibration factor which the observed spectrum will be multiplied by
         #The second is the fraction of correlated noise assumed
         #And the third is the scale length (in microns) of the correlated component of the noise
         1.0 ,0.1, 0.1
        ] #
        + np.random.rand(optimizer.npars)*[1,1,1,1,1,1,1,1,
                                           #1,1,1,
                                           1,1,1
                                           ]
        for i in range(optimizer.nwalkers)])

    optimizer.postProcess() #now we call the postprocessing to produce some figures

