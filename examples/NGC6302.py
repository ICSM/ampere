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
       logacold : log of (relative) abundances for the cold component 
       logawarm : log of (relative) abundances for the warm component
       Tcold : temperature range of the cold component (low, high)
       Twarm : temperature range of the warm component (low, high)
       indexp : index p of the density distribution 
       indexq : index q of the temperature distribution
       Kemper et al. 2002 use p and q of 0.5.
       Output is an array of model fluxes (fluxes), to match wavelengths
    '''
    # Kemper et al. 2002 talk about p and q being -1/2, but what they actually
    # mean is that -p and -q are -1/2, or p and q are 1/2. 
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
            #entire spectral range, from 2.36 to 196.6 micron
            #the extrapolation is done in log space for the opacities for
            #better results.

            # calcite and dolomite are M.A.C. (kappa) values
            # enstatite and diopside are Q/a values
            # we need to convert those to Q values
                        
            #plt.xscale('log')
            #plt.plot(self.wavelength, opacity_array[:,j])
            #plt.plot(tempWl,tempOpac)
            #plt.title(opacityFileList[j])
            #plt.xlim(2,200)
            #plt.show()
            print("Reading in species: ", j, " : ", opacityFileList[j])

        opacity_array[:,1] = opacity_array[:,1]*1e-4 # convert enstatite to Q
        opacity_array[:,3] = opacity_array[:,3]*1e-4 # convert diopside to Q

        opacity_array[:,0] = opacity_array[:,0]*2.71*(4/3)*1e-4 # convert calcite from kappa to Q. Calcite has a density of 2.71 g cm^-3.
        opacity_array[:,4] = opacity_array[:,4]*2.87*(4/3)*1e-4 # convert calcite from kappa to Q. Calcite has a density of 2.71 g cm^-3. 


 
        self.opacity_array = opacity_array
        self.nSpecies = nSpecies
        self.npars = 15
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
            # logacold = theta[0:8] -6, 0
            # logawarm = theta[8:16] -6, 0
            # Tcold = theta[16:18] 10, 80; 10, 80 
            # Twarm = theta[18:20] 80, 200; 80, 180
            # indexp = theta[20] 0, 1;
            self.lims[:, 0] = -6.
            self.lims[11:13, 0] = 10.
            self.lims[11:13, 1] = 80.
            self.lims[13:15, 0] = 80.
            self.lims[13:15, 1] = 180.
#            self.lims[15, 0] = 0. 
#            self.lims[15, 1] = 1.
        else:            
            self.lims = lims
        print("Limits: ")
        print(self.lims)
        self.flatprior = flatprior


    def __call__(self, logacold0, logacold1, logacold2, logacold3, logacold4,
                 logacold6, logacold7, 
                 logawarm1, logawarm2, logawarm5,
                 logawarm7, 
                 Tcold0, Tcold1, Twarm0, Twarm1, #indexp, 
                 *args, **kwargs):
        '''The model itself, using the callable class functionality of python.
        This is an essential method. It should do any steps required to
        calculate the output fluxes. Once done, it should stop the output
        fluxes in self.modelFlux.
        '''
        # number of dust species in opacityDirectory.
        # for each we will consider a hot and a cold component. We will allow
        # only 2 temperature ranges, the same for all dust species.
        # we only fit the species that were non-zero in the original model fit,
        # e.g. 7 species in the cold component, and 3 species in the hot
        # component
        # the module translated from ck_modbb calculates the flux for each
        # component at a single temperature
        # ckmodbb and shbb are now written.

        coldcomponent = np.zeros((2, wavelengths.__len__()))
        coldcomponent[0, :] = self.wavelength
        warmcomponent = np.zeros((2, wavelengths.__len__()))
        warmcomponent[0, :] = self.wavelength
        acold = [10**logacold0, 10**logacold1, 10**logacold2, 10**logacold3, 10**logacold4, 0, 10**logacold6,
                 10**logacold7] #0 values correspond to dust species that were not included in the fit by Kemper et al. 2002
        awarm = [0, 10**logawarm1, 10**logawarm2, 0, 0, 10**logawarm5, 0,
                 10**logawarm7] #same here       

        for ii, aa in enumerate(acold):
            onespeciescold = self.ckmodbb(self.opacity_array[:, ii],
                                          tin=Tcold1, tout=Tcold0,
                                          n0=aa)#, index=indexp)
            coldcomponent[1, :] = coldcomponent[1, :] + onespeciescold[1, :]
            #plt.plot(coldcomponent[0,:],coldcomponent[1,:])
            
        for jj, bb in enumerate(awarm):
            onespecieswarm = self.ckmodbb(self.opacity_array[:, jj],
                                          tin=Twarm1, tout=Twarm0,
                                          n0=bb)#, index=indexp)
            warmcomponent[1, :] = warmcomponent[1, :] + onespecieswarm[1, :]
            #plt.plot(warmcomponent[0,:],warmcomponent[1,:])
        fModel = np.zeros((2, wavelengths.__len__()))
        fModel[0, :] = self.wavelength
        fModel[1, :] = coldcomponent[1,:] + warmcomponent[1, :]
        self.modelFlux = fModel[1, :]
        #plt.plot(fModel[0,:],fModel[1,:])
        #plt.show()

        return {"spectrum":{"wavelength":self.wavelength, "flux": self.modelFlux}}

    def ckmodbb(self, q, tin, tout, n0, index=0.5, r0=1e15, distance=910.,
                grainsize=0.1, steps=15):
        ## ckmodbb calculates the flux of either the warm or the cold component
        ## according to equation 6 and 7 from Kemper et al. (2002, A&A 394, 679)
        ## translated from ck_modbb.pro (IDL script)
        d = distance * 3.0857e18  # convert distance from pc to cm
        a = grainsize * 1e-4  # convert grainsize from micron to cm
        fnu = np.zeros((2, self.wavelength.__len__()))
        fnu[0, :] = self.wavelength  #maybe it is better to  make fnu
                                     #(and therefore onespecies)
                                     #an array with just 1 column, remove wl
        pindex = index
        qindex = index

        for j in range(steps - 1):
            t = tin - j * (tin-tout)/steps # integrating. This is how it was
            # done by Kemper et al. 2002. This can be done better, but I'll
            # leave it as is for now.
            power = (t/tin)**(-1*(3-pindex)/qindex)
            bb = self.shbb(fnu, t, 0.)
            fnu[1, :] = np.add(fnu[1, :], 
                               np.multiply(q,
                                           bb[1, :])*(power*((tin-tout)/steps)))
            # adding the fnu together calculated in the different steps. 
        
        factor = (4*math.pi*a*a * r0**3 * n0)/((3-pindex)*d*d)
        # this factor from eq. 6 (Kemper et al. 2002) can be placed outside the
        # integral

        fnu[1, :] = fnu[1, :] * factor
        return fnu

    def shbb(self, aar, temp, pinda):
        ## This calculates the blackbody flux, if necessary multiplied with a
        ## powerlaw opacity.
        ## We call it with the index 0, so no power-law opacity
        ## based on sh_bb.pro (IDL script) by Sacha Hony. His original comments:
        ## ;(SH Jan 27 1999)
        ## ;wavelength dependent emissivity powerindex:power
        ## ; fit = bb(T)*H*l**power        
        a1 = 3.97296e19
        a2 = 1.43875e4
        mbb = np.copy(aar)
        bbflux = a1/(mbb[0, :]**3)/(np.exp(a2/(mbb[0, :]*temp))-1)
        mbb[1, :] = bbflux * mbb[0, :]**pinda
        return mbb

    def lnprior(self, theta, **kwargs): 
        # acold = theta[0:7]
        # awarm = theta[7:10]
        # Tcold = theta[10:12]; theta[11] > theta[10]
        # Twarm = theta[12:14]; theta[13] > theta[12]
        # indexp = theta[14]
 

        if self.flatprior: 
            if (np.all([self.lims[i,0] <= theta[i] <= self.lims[i,1] for i in
                       range(len(self.lims[:,0]))]) and (theta[12] > theta[11])
                      and (theta[14] > theta[13])):
                # the temperature of the inner radius needs to be higher than
                # the temperature of the outer radius for both temperature
                # components. Hence theta[12] > theta[11] and theta[14] >
                # theta[13].
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
    wave1 = np.linspace(2.3603,35.0603,327)
    wave2 = np.linspace(1/196.6261,1/35.1,117)
    wavelengths = np.concatenate((wave1,1/wave2[::-1]))
    #print(wavelengths[:])

    """ Choose some model parameters """
    acold0 = -4.9  # Cold calcite (converted to Q)
    acold1 = -4 # Cold enstatite (converted to Q)
    acold2 = -2.7 # Cold forsterite (Q)
    acold3 = -4.5 # Cold diopside (converted to Q)
    acold4 = -5.1 # Cold dolomite (converted to Q)
#    acold5 = -10 # Cold metallic iron. Not used by Kemper et al. 2002
    acold6 = -3 # Cold crystalline water ice (Q)
    acold7 = -1.4 # Cold amorphous olivine (Q)
#    awarm0 = -10 # Warm calcite. Not used by Kemper et al. 2002
    awarm1 = -5 # Warm enstatite. (converted to Q)
    awarm2 = -5 # Warm forsterite (Q)
#    awarm3 = -10 # Warm diopside. Not used by Kemper et al. 2002
#    awarm4 = -10 # Warm dolomite. Not used by Kemper et al. 2002
    awarm5 = -3 # Warm metallic iron (Q).
#    awarm6 = -10 # Warm crystalline water ice. Not used in Kemper et al. 2002. Would have been evaporated at these temperatures anyway
    awarm7 = -4 # Warm amorphous olivine (Q)
    Tcold0 = 30 # Temperature outer radius cold dust shell
    Tcold1 = 60 # Temperature inner radius cold dust shell
    Twarm0 = 100 # Temperature outer radius warm dust shell
    Twarm1 = 118 # Temperature inner radius warm dust shell
#    indexp = 0.5 # Index p from Kemper et al. 2002. Kept fixed for now to reduce the number of parameters 


    #Now init the model:
    model = SpectrumNGC6302(wavelengths)
    #And call it to produce the fluxes for our chosen parameters
    model(acold0, acold1, acold2, acold3, acold4, 
                 acold6, acold7, 
                 awarm1, awarm2, awarm5,
                 awarm7, #zero-value dust species are removed from the call
          Tcold0, Tcold1, Twarm0, Twarm1)
#          indexp=indexp)
    model_flux = model.modelFlux
    #plt.plot(wavelengths, model.modelFlux) #sanity check plot. Can be commented
    #plt.show()                             # out.
    
    #Now we create synthetic data:

    #now we'll create a synthetic spectrum from the model fluxes, using the to be fitted spectrum to get the wavelength sampling
    dataDir = os.getcwd() + '/NGC6302/'
    specFileExample = 'NGC6302_100.tab'
    specdata = ascii.read(dataDir+specFileExample,data_start=2)
        #And again, add some noise to it
    input_noise_spec = 0.01 #assume a 1% error on the spectrum
    unc = specdata[1][:]*input_noise_spec
    spec = Spectrum(specdata[0][:],specdata[1][:] +
                    np.random.randn(len(specdata[1][:]))*unc,specdata[1][:]*0.05,"um","Jy", calUnc=1e-10, scalelengthPrior=0.1) #added extra keyword calUnc with a small value, and scalelengthPrior to restrict cov. scale length
    #plt.plot(spec.wavelength, spec.value) #another sanity check
    #plt.show()                            #comment out if it is bothersome


    #Now let's try changing the resampling method so it's faster
    #This model is very simple so exact flux conservation is not important
    resmethod = "fast" #"exact"#"fast"#
    spec.setResampler(resampleMethod=resmethod)
#    spec1.setResampler(resampleMethod=resmethod)

    spec.selectWaves(low=25, up=120) #limit wavelength range over which fit is performed -- similar to 2002 result
    
    """ now set up ampere to try and fit the same stuff """

    dataset = [               
               spec   
               ]


    #Ampere exposes access to emcee's moves interface. This can be useful if the posterior turns out to not be well behaved - the default move only deals well with posteriors that are monomodal and approximately Gaussian. Here's an example that usually deals a bit better with posteriors that don't meet these criteria:
    m = [(moves.DEMove(), 0.8),
        (moves.DESnookerMove(), 0.2),
         ]

    #Now we set up the optimizer object:
    #optimizer = EmceeSearch(model=model, data=dataset, nwalkers=100, moves=m)

    optimizer = EmceeSearch(model=model, data=dataset, nwalkers=50, moves=m)
    guess = [
        [-4.9, -4, -2.7, -4.5, -5.1, -3, -1.4, #values are essentially the 
         -5, -5,  -3,   -4,                        #same as fitted values from
         30, 60,                               # 2002.
         100, 118,
         #The parameters of the model
             #Each Spectrum object contains a noise model with three free parameters
             #The first one is a calibration factor which the observed spectrum will be multiplied by
             #The second is the fraction of correlated noise assumed
             #And the third is the scale length (in microns) of the correlated component of the noise
         1.0 ,0.1, 0.1
        ]
        + np.random.rand(optimizer.npars)*[0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                                           0.1, 0.1, 0.1, 0.1, #these values are
                                           10, 10,          #adjusted so that 
                                           10, 10,          #the first guesses
                                           #0.1,       #do not fall outside
                                           1, 1, 1]       #the prior range
        for i in range(optimizer.nwalkers)]
    #guess = "None"

    #Then we tell it to explore the parameter space

    optimizer.optimise(nsamples = 50000, burnin=40000, guess=guess)
#    optimizer.optimise(nsamples = 50, burnin=10, guess=guess) #short run for tests


    optimizer.postProcess() #now we call the postprocessing to produce some figures
