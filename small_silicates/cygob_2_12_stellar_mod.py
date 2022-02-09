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
from astropy.constants import c

#First we will define a rather simple model

class Blackbody_dust(Model):
    '''This is a blackbody model based on the very simple model from Peter 

    We will first try to fit a blackbody model to a dataset with dust 

    '''
    def __init__(self, wavelengths, flatprior=True,
                 lims=np.array([[10., 300000.],
                                [1.,6.],
                                [0.001,10.]])):
        '''The model constructor, which will set everything up

        This method does essential setup actions, primarily things that 
        may change from one fit to another, but will stay constant throughout 
        the fit. This may be things like the grid of wavelengths to calculate
        the model output on, or establishing the dust opacities if involved.
        There are also several important variables it *MUST* define here
        '''
        self.wavelength = wavelengths
        
        stellar_modDir = os.path.dirname(os.path.realpath('__file__'))
        stellar_mod = np.loadtxt(stellar_modDir + 'ssap.txt', comments = '#')
        temp_wl_mod = stellar_mod[:,0]
        temp_flux_mod = stellar_mod[:,1] #units ERG/CM2/S/A
        wavelengths_aa=wavelengths
        f1 = interpolate.interp1d(temp_wl_mod, temp_flux_mod, assume_sorted = False)
        stellar_flux = f1(wavelengths_aa) #units ERG/CM2/S/A
        
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
            

            f2 = interpolate.interp1d(tempWl, tempOpac, assume_sorted = False)
            opacity_array[:,j] = f2(wavelengths)#wavelengths)
            
        self.stellar_flux = stellar_flux
        self.opacity_array = opacity_array
        self.nSpecies = nSpecies
        
        self.npars = nSpecies + 3 #Number of free parameters for the model (__call__()). For some models this can be determined through introspection, but it is still strongly recommended to define this explicitly here. Introspection will only be attempted if self.npars is not defined.
        self.npars_ptform = nSpecies + 3 #Sometimes the number of free parameters is different when using the prior transform instead of the prior. In that case, self.npars_ptform should also be defined.
        #You can do any other set up you need in this method.
        #For example, we could define some cases to set up different priors
        #But that's for a slightly more complex example.
        #Here we'll just use a simple flat prior
        self.lims = lims
        self.flatprior = flatprior

    def __call__(self, temp, radius_sol, scaling, *args, **kwargs):
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
        pc  = 3.086e16 # m
        rsol = 696340e3 # m
        
        distance_pc=1100.
        
        Rstar = radius_sol*rsol
        r = distance_pc*pc
        # from angstrom to micron 
        # [Y Jy]                       = 3.33564095E+04 * [X1 erg/cm^2/s/A] * [X2 A]^2
        flux = self.stellar_flux*wavelengths_aa**2 * 3.33564095E+04 * (Rstar/r)**2   # conversion to jansky's 
        flux_freefree = (((c.to(u.cm/u.s))/wavelengths_aa.to(u.cm))**0.6).value
        flux_mjy = (flux.to(u.Jy / u.sr).value*(Rstar/r)**2.+scaling*flux_freefree)*fModel2
        # here we need to put dust models ...
        dustAbundances = np.array(args) # instead of 10**np.array(args)
        
        self.modelFlux = flux_mjy #slope*self.wavelength + intercept        

    def lnprior(self, theta, **kwargs):
        temp = theta[0]
        #distance_pc  = theta[1]
        radius_sol= theta[1]
        scaling = theta[2]
        if self.flatprior:
            if np.all(theta[3:] > 0.) and (self.lims[0,0] < theta[0] < self.lims[0,1]) and (self.lims[1,0] < theta[1] < self.lims[1,1])and (self.lims[2,0] < theta[2] < self.lims[2,1]):
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
    
    # here should be the model input
    wavelengths = np.linspace(5.0,38, 1000)
    
    # constants and definition of R factor
    pc  = 3.086e16 # m
    rsol = 696340e3 # m

    """ Choose some model parameters """
    temp = 12000.#Keep it super simple for now
    radius_sol = 2.9
    distance_pc = 1100. # perhaps you should fix the distance!!
    scaling=1.0

    #Now init the model:
    model = Blackbody_dust(wavelengths)
    #And call it to produce the fluxes for our chosen parameters
    model(temp, radius_sol, scaling, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
    model_flux = model.modelFlux
    
    # Here starts the data part 
    

    #now we'll create a synthetic spectrum from the model fluxes, using a Spitzer IRS observation to get the wavelength sampling
    dataDir = ampere.__file__.strip('__init__.py') + 'Testdata/'
    specFile1 = 'cassis_yaaar_spcfw_27570176t.fits'
    irsEx_1 = Spectrum.fromFile(dataDir+specFile1,format='SPITZER-YAAAR')
    
    specFile2 = 'cassis_yaaar_spcfw_9834496t.fits'
    irsEx_2 = Spectrum.fromFile(dataDir+specFile2,format='SPITZER-YAAAR')
    
    specFile3 = 'cassis_yaaar_optdiff_9834496.fits'
    irsEx_3 = Spectrum.fromFile(dataDir+specFile3,format='SPITZER-YAAAR_OPTDIFFHR')
    
    dataSet = irsEx_1
    
    for s in irsEx_2:             #include the next two lines when appending spectroscopy to photometry
        dataSet.append(s)

#    for s in irsEx_3:             #include the next two lines when appending spectroscopy to photometry
#        dataSet.append(s)
                
   # print(dataSet[0])   
        
    fig = plt.figure()    
    ax = fig.add_subplot(111)
    ax.set_yscale('log')
    
    for i in dataSet:
        ax.plot(irsEx_1[0].wavelength, irsEx_1[0].value, '-',color='red')
        ax.plot(irsEx_1[1].wavelength, irsEx_1[1].value, '-',color='red')        
        ax.plot(irsEx_2[0].wavelength, irsEx_2[0].value, '-',color='red')
        ax.plot(irsEx_3[0].wavelength, irsEx_3[0].value, '-',color='blue')
        #ax.plot(irsEx_3[1].wavelength, irsEx_3[1].value, '-',color='blue')
    
    ax.plot(wavelengths, model_flux)

    plt.show()
    
    #Ampere exposes acces to emcee's moves interface. This can be useful if the posterior turns out to not be well behaved - the default move only deals well with posteriors that are monomodal and approximately Gaussian. Here's an example that usually deals a bit better with posteriors that don't meet these criteria:
    m = [(moves.DEMove(), 0.8),
        (moves.DESnookerMove(), 0.2),
         ]

    #Now we set up the optimizer object:
    optimizer = EmceeSearch(model=model, data=dataSet, nwalkers=100, moves=m)
            
    optimizer.optimise(nsamples = 2, burnin=1, guess=[
        [12000., 2.9, scaling, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, #The parameters of the model
         #1.0, 0.1, 0.1, #Each Spectrum object contains a noise model with three free parameters
         #The first one is a calibration factor which the observed spectrum will be multiplied by
         #The second is the fraction of correlated noise assumed
         #And the third is the scale length (in microns) of the correlated component of the noise
         1.0 ,0.1, 0.1,1.0 ,0.1, 0.1,1.0 ,0.1, 0.1
        ] #
        + np.random.rand(optimizer.npars)*[1,1,1,1,1,1,1,1,1,
                                           1,1,1,
                                           1,1,1,
                                           1,1,1
                                           ]
        for i in range(optimizer.nwalkers)])

    optimizer.postProcess() #now we call the postprocessing to produce some figures
    
    
    
        
