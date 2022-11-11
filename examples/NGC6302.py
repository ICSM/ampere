import numpy as np
import os
import math
import ampere
from ampere.data import Spectrum, Photometry
from ampere.infer.emceesearch import EmceeSearch
from ampere.models import Model
from spectres import spectres
import pyphot
from emcee import moves
import matplotlib.pyplot as plt
from astropy.io import ascii


class SpectrumNGC6302(Model):
    '''This model fits a modfied blackbody multiplied by sum of opacties,
        consisting of a warm and cold component, not necessarily of the same
        composition, over two temperature ranges, following Kemper et al. 2002.
        We use the following parameters:
       wavelengths : the array of wavelengths considered
       opacityFileList
       acold : (relative) abundances for the cold component 
       awarm : (relative) abundances for the warm component
       Tcold : temperature range of the cold component (low, high)
       Twarm : temperature range of the warm component (low, high)
       indexp : index p of the density distribution 
       indexq : index q of the temperature distribution
       multfact : multiplication factor that the sum of the modified black
                  bodies has to be multiplied with to fit the spectrum
       Output is an array of model fluxes (fluxes), to match wavelengths
    '''
    def __init__(self, wavelengths, flatprior=True,
                 opacityFileName='NGC6302-opacities.txt', lims=None):
        '''The model constructor, which will set everything up
        This method does essential setup actions, primarily things that
        may change from one fit to another, but will stay constant throughout
        the fit. This may be things like the grid of wavelengths to calculate
        the model output on, or establishing the dust opacities if involved.
        There are also several important variables it *MUST* define here
        '''
        self.wavelength = wavelengths
        import os
        from scipy import interpolate
        opacityDirectory = os.getcwd()+'/NGC6302/'
        opacityFileList = np.loadtxt(opacityFileName, dtype='str')
        nSpecies = opacityFileList.__len__() 
        opacity_array = np.zeros((wavelengths.__len__(), nSpecies))
        for j in range(nSpecies):
            tempData = np.loadtxt(opacityDirectory + opacityFileList[j],
                                  comments='#')
            tempWl = tempData[:, 0]
            tempOpac = tempData[:, 1]
            f = interpolate.interp1d(tempWl, np.log10(tempOpac), assume_sorted=False, fill_value = "extrapolate") 
            opacity_array[:, j] = np.power(np.full(self.wavelength.shape,10),f(self.wavelength))
            #extrapolate is needed because some opacities don't cover the
            #entire ISO SWS-LWS range, from 2.4 to 198 micron
            #the extrapolation is done in log space for the opacities for
            #better results. 
            #plt.xscale('log')
            #plt.plot(self.wavelength, opacity_array[:,j])
            #plt.plot(tempWl,tempOpac)
            #plt.title(opacityFileList[j])
            #plt.show()
            print("Reading in species: ", j, " : ", opacityFileList[j])

 
        self.opacity_array = opacity_array
        self.nSpecies = nSpecies
        self.npars = 2*nSpecies + 6
        # the number of free parameters for the model (__call__()). For
        # some models this can be determined through introspection, but it is
        # still strongly recommended to define this explicitly here.
        # Introspection will only be attempted if self.npars is not defined.
        self.npars_ptform = 2
        # Sometimes the number of free parameters is different when using the
        # prior transform instead of the prior. In that case, self.npars_ptform
        # should also be defined.
        # You can do any other set up you need in this method.
        # For example, we could define some cases to set up different priors
        # But that's for a slightly more complex example.
        # Here we'll just use a simple flat prior
        if lims is None:
            self.lims = np.zeros((self.npars,2)) 
            # acold = theta[0:8] 0, 1
            # awarm = theta[8:16] 0, 1
            # Tcold = theta[16:18] 10, 80; 10, 80 
            # Twarm = theta[18:20] 80, 200; 80, 200
            # indexp = theta[20] 0, 3;
            # multfact = theta[21] 0, np.inf
            self.lims[:, 1] = 1
            self.lims[16:18, 0] = 10.
            self.lims[16:18, 1] = 80.
            self.lims[18:20, 0] = 80.
            self.lims[18:20, 1] = 200.
            self.lims[20, 1] = 3.
            self.lims[21, 1] = np.inf
        else:            
            self.lims = lims
        self.flatprior = flatprior


    def __call__(self, acold0, acold1, acold2, acold3, acold4, acold5,
                 acold6, acold7, 
                 awarm0, awarm1, awarm2, awarm3, awarm4, awarm5,
                 awarm6, awarm7, 
                 Tcold0, Tcold1, Twarm0, Twarm1, indexp, multfact,
                 *args, **kwargs):
        '''The model itself, using the callable class functionality of python.
        This is an essential method. It should do any steps required to
        calculate the output fluxes. Once done, it should stop the output
        fluxes in self.modelFlux.
        '''
        # number of dust species in opacityDirectory are all going to be fitted
        # for each we will consider a hot and a cold component. We will allow
        # only 2 temperature ranges, the same for all dust species.
        # the module translated from ck_modbb calculates the flux for each
        # component at a single temperature
        # ckmodbb and shbb are now written. Now I just need to figure out how
        # to do the fit over 8 dust species and 2 temperature ranges.
        # let's first do this over fixed temperature ranges, instead of
        # allowing them to be free the temperature ranges are 118-100 K and
        # 60-30 K, following Kemper et al. 2002

#        print(acold0, acold1, acold2, acold3, acold4, acold5,
#                 acold6, acold7, 
#                 awarm0, awarm1, awarm2, awarm3, awarm4, awarm5,
#                 awarm6, awarm7, 
#                 Tcold0, Tcold1, Twarm0, Twarm1, indexp, multfact)

        
        coldcomponent = np.zeros((2, wavelengths.__len__()))
        coldcomponent[0, :] = self.wavelength
        warmcomponent = np.zeros((2, wavelengths.__len__()))
        warmcomponent[0, :] = self.wavelength
        acold = [acold0, acold1, acold2, acold3, acold4, acold5, acold6,
                 acold7]
        awarm = [awarm0, awarm1, awarm2, awarm3, awarm4, awarm5, awarm6,
                 awarm7] 
 
#        print("acold: ",acold)
        for i, a in enumerate(acold):
            onespeciescold = self.ckmodbb(self.opacity_array[:, i],
                                          tin=Tcold[1], tout=Tcold[0],
                                          n0=a, index=indexp)
            coldcomponent[1, :] = coldcomponent[1, :] + onespeciescold[1, :]
            #print(max(coldcomponent[1,:]))
            #plt.plot(coldcomponent[1,:])
            
        for i, a in enumerate(awarm):
            onespecieswarm = self.ckmodbb(self.opacity_array[:, i],
                                          tin=Twarm[1], tout=Twarm[0],
                                          n0=a, index=indexp)
            warmcomponent[1, :] = warmcomponent[1, :] + onespecieswarm[1, :]
            #print(max(warmcomponent[1,:]))
            #plt.plot(warmcomponent[1,:])
        fModel = np.full_like(coldcomponent, 1.)
        fModel[1, :] = coldcomponent[1,:] + warmcomponent[1, :]
        self.modelFlux = fModel[1, :]
        
        #plt.show()
        

        return {"spectrum":{"wavelength":self.wavelength, "flux": self.modelFlux}}

    def ckmodbb(self, q, tin, tout, n0, index=0.5, r0=1e15, distance=910.,
                grainsize=0.1, steps=10):
        d = distance * 3.0857e18  # convert distance from pc to cm
        a = grainsize * 1e-4  # convert grainsize from micron to cm
        fnu = np.zeros((2, self.wavelength.__len__()))
        fnu[0, :] = self.wavelength  #TODO: make fnu (and therefore onespecies)
                                     #an array with just 1 column, remove wl

        for j in range(steps - 1):
            t = tin - j * (tin-tout)/steps
            power = (t/tin)**(2*index - 6)
            bb = self.shbb(fnu, t, 0.)
          #  print(min(bb[1,:]), max(bb[1,:]))
          #  print(power*((tin-tout)/steps))
            fnu[1, :] = np.add(fnu[1, :], 
                               np.multiply(q,
                                           bb[1, :])*(power*((tin-tout)/steps)))
        extra = r0/d
        factor = 4 * math.pi * a * a * r0 * n0 * extra * extra / (3-index)
      #  print("pi: ", math.pi)
      #  print("a: ", a)
      #  print("r0: ", r0)
      #  print("n0: ", n0)
      #  print("extra: ", extra)
      #  print("index: ", index)
#        print("factor: ", factor)

        fnu[1, :] = fnu[1, :] * factor
        
        return fnu

    def shbb(self, aar, temp, pinda):
        a1 = 3.97296e19
        a2 = 1.43875e4
        mbb = np.copy(aar)
        bbflux = a1/(mbb[0, :]**3)/(np.exp(a2/(mbb[0, :]*temp))-1)
        mbb[1, :] = bbflux * mbb[0, :]**pinda
        return mbb

    def lnprior(self, theta, **kwargs): 
        # acold = theta[0:8]
        # awarm = theta[8:16]
        # Tcold = theta[16:18]; theta[17] > theta[16]
        # Twarm = theta[18:20]; theta[19] > theta[18]
        # indexp = theta[20]
        # multfact = theta[21]
        
        if self.flatprior:
            if np.all([self.lims[i,0] <= theta[i] <= self.lims[i,1] for i in
                      range(len(self.lims[:,0]))] and (theta[17] > theta[16])
                      and (theta[19] > theta[18])):
                return 0
            else:
                return -np.inf
        else:
            raise NotImplementedError()

    def prior_transform(self, u, **kwargs):  
        '''The prior transform, which takes samples from the Uniform(0,1)
        distribution to the desired distribution.
        This is only included for completeness and to demonstrate how a prior 
        transform function should look. This example only uses emcee for 
        fitting, which uses the lnprior function instead. Prior transforms are 
        required by nested-sampling codes and similar approaches.
        '''
        if self.flatprior:
            theta = np.zeros_like(u)
            return (self.lims[:,1] - self.lims[:,0]) * u + self.lims[:,0]
        else:
            raise NotImplementedError()

if __name__ == "__main__": 
    """ Set up the inputs for the model """
    """ wavelength grid """
    wavelengths = np.linspace(2.4,198.,1956)

    """ Choose some model parameters """
    acold = [5e-9, 1e-8,  2e-3, 3e-9, 3e-9, 0,    1e-3, 4e-2]
    awarm = [0,    1e-10, 1e-5, 0,    0,    1e-3, 0,    1e-4]
    Tcold = [30, 60]
    Twarm = [100, 118]
    indexp = 0.5
    multfact = 1.

    #Now init the model:
    model = SpectrumNGC6302(wavelengths)
    #And call it to produce the fluxes for our chosen parameters
    model(acold[0], acold[1], acold[2], acold[3], acold[4], acold[5],
                 acold[6], acold[7], 
                 awarm[0], awarm[1], awarm[2], awarm[3], awarm[4], awarm[5],
                 awarm[6], awarm[7], 
          Tcold[0], Tcold[1], Twarm[0], Twarm[1],
          indexp=indexp, multfact=multfact)
    model_flux = model.modelFlux
    plt.plot(wavelengths, model.modelFlux)
    plt.show()
    
    #Now we create synthetic data:

    #now we'll create a synthetic spectrum from the model fluxes, using the to be fitted spectrum to get the wavelength sampling
    dataDir = os.getcwd() + '/NGC6302/'
    specFileExample = 'NGC6302_100.tab'
    specdata = ascii.read(dataDir+specFileExample,data_start=2)
    input_noise_spec = 0.05 #assume a 5% error on the flux measurements. check with what I did in 2002
    unc = specdata[1][:]*input_noise_spec
    spec = Spectrum(specdata[0][:],specdata[1][:] +
                    np.random.randn(len(specdata[1][:]))*unc,specdata[1][:]*0.05,"um","Jy")
#    spec0 = spectres(spec[0].wavelength,wavelengths,model_flux)
#    spec1 = spectres(spec[1].wavelength,wavelengths,model_flux)
    
    #And again, add some noise to it
#    input_noise_spec = 0.05
#    unc0 = input_noise_spec*spec0
#    unc1 = input_noise_spec*spec1
#    spec0 = spec0 + np.random.randn(len(spec0))*unc0
#    spec1 = spec1 + np.random.randn(len(spec1))*unc1
    
#    spec0 = Spectrum(irsEx[0].wavelength, spec0, unc0,"um", "Jy",calUnc=0.0025, scaleLengthPrior = 0.01) #, resampleMethod=resmethod)
#    spec1 = Spectrum(irsEx[1].wavelength, spec1, unc1,"um", "Jy",calUnc=0.0025, scaleLengthPrior = 0.01) #, resampleMethod=resmethod)

    #Now let's try changing the resampling method so it's faster
    #This model is very simple so exact flux conservation is not important
    resmethod = "fast" #"exact"#"fast"#
    spec.setResampler(resampleMethod=resmethod)
#    spec1.setResampler(resampleMethod=resmethod)

    """ now set up ampere to try and fit the same stuff """

    dataset = [               
               spec   
               ]


    #Ampere exposes acces to emcee's moves interface. This can be useful if the posterior turns out to not be well behaved - the default move only deals well with posteriors that are monomodal and approximately Gaussian. Here's an example that usually deals a bit better with posteriors that don't meet these criteria:
    m = [(moves.DEMove(), 0.8),
        (moves.DESnookerMove(), 0.2),
         ]

    #Now we set up the optimizer object:
    #optimizer = EmceeSearch(model=model, data=dataset, nwalkers=100, moves=m)

    optimizer = EmceeSearch(model=model, data=dataset, nwalkers=100, moves=m)
    guess = [
        [5e-9,  1e-8,  2e-3, 3e-9,  3e-9,  1e-10, 1e-3,  4e-2,
         1e-10, 1e-10, 1e-5, 1e-10, 1e-10, 1e-3,  1e-10, 1e-4,
         30, 60,
         100, 118,
         0.5,
         1., #The parameters of the model
             #Each Spectrum object contains a noise model with three free parameters
             #The first one is a calibration factor which the observed spectrum will be multiplied by
             #The second is the fraction of correlated noise assumed
             #And the third is the scale length (in microns) of the correlated component of the noise
         1.0 ,0.1, 0.1
        ]
        + np.random.rand(optimizer.npars)*[1e-9, 1e-8, 1e-3, 1e-9, 1e-9, 1e-10, 1e-3, 1e-2,
                                           1e-10, 1e-10, 1e-5, 1e-10, 1e-10, 1e-3, 1e-10, 1e-4, 
                                           10, 10, 10, 10, 0.5, 1,
                                           1,0.1,0.1
        ]
        for i in range(optimizer.nwalkers)]
    #guess = "None"

    #Then we tell it to explore the parameter space

    #optimizer.optimise(nsamples = 1500, burnin=1000, guess=guess)
    optimizer.optimise(nsamples = 50, burnin=10, guess=guess)


    optimizer.postProcess() #now we call the postprocessing to produce some figures
