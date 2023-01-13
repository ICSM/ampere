#
# Ampere Stuff
#
import copy
import numpy as np
from astropy import constants as const
from astropy import units
try:
    from astropy.modeling import blackbody
except ImportError:
    from astropy.modeling.physical_models import BlackBody
from .models import AnalyticalModel

from scipy.stats import dirichlet

# 
# Hyperion Stuff
# 

import hyperion
from hyperion.model import Model, AnalyticalYSOModel
from hyperion.model import ModelOutput
from hyperion.dust import BHDust
#from hyperion.util.constants import lsun, rsun, au, msun, sigma, pc, pi #get these from astropy units 

import subprocess

# Silly string to link Ampere and Hyperion together
class HyperionRTModel(Model):
    '''
    Class to link Ampere and Hyperion.

    '''

    def __init__(self,flatprior=True,                 
                 #RT parameters - use in __init__ and store as self.XXX
                        niter=5,
                        nph_initial=1e3,
                        nph_imging=1e4,
                        nph_rtsrcs=1e4,
                        nph_rtdust=1e4,
                        track_mode='detailed', #one of no/basic/detailed/scatterings
                        raytracing=True,
                        modrndwlk=True,
                        mrw_gamma=2,
                        lmin=0.1,lmax=1000.0,nl=101,
                 #Peel photons to get images - use in __init__ and store as self.XXX
                        api_sed=True,
                        api_img=False,
                        view_angles=np.linspace(0., 90., 10),
                        view_repeats=np.repeat(45., 10),
                 #Parallel processing - use in __init__ and store as self.XXX
                        useMPI=True,
                        nproc=60,
                        npars=0,
                        **kwargs):
        
        #Assign keyword variables to the object
        #Assign all the inputs to __init__ to instance variables with the same name as above
        #this is equivalent to many lines of self.niter = niter etc
        l = locals()
        for key, value in l.items():
            if key == "self": #skip self to avoid possible recursion
                continue
            setattr(self, key, value) #builtin that assigns to an attribute with name = key of self with a value of value
        
        
        #self.flatprior = flatprior
        
        self.model = Model()
        
        #Use raytracing to improve s/n of thermal/source emission
        self.model.set_raytracing(self.raytracing)
        
        #Use the modified random walk
        self.model.set_mrw(self.modrndwlk, gamma=self.mrw_gamma)
        
        #Set up SED for 10 viewing angles
        sed = self.model.add_peeled_images(sed=self.api_sed, image=self.api_img)
        sed.set_viewing_angles(self.view_angles,self.view_repeats)
        sed.set_wavelength_range(nl, lmin, lmax)
        sed.set_track_origin(self.track_mode)

        self.sed = sed

        #Set number of photons
        self.model.set_n_photons(initial=self.nph_initial, imaging=self.nph_imging,
                        raytracing_sources=self.nph_rtsrcs, raytracing_dust=self.nph_rtdust)
        
        #Set number of temperature iterations
        self.model.set_n_initial_iterations(self.niter)

        print("Hyperion RT Model setup complete.")
        #self.dictionary of fixed parameters that are needed for modelling in __call__

    #Only give __call__ the numbers that emcee is going to change.
    def __call__(self,dust="astrosilicate",fileformat=2,amin=[0.5],amax=[1000.],na=[101],
                        nang=91,nanx=11,
                        nchem=1,gtd=0,
                        tmin=2.7,tmax=1000.0,nt=101,
                        massfrac=[1.0],
                        density=[3.3],
                        disttype=['power'],
                        optconst=["silicate_d03.lnk"],
                        q=[3.5],
                 #Source parameters
                        #sources=[['spherical',1.0,1.0,1.0,5784,[0.0,0.0,0.0]]],(type,lstar,rstar,mstar,tstar,position[x,y,z],spectrum file)
                        sources =[['spherical',1.0,1.0,1.0,5784,[0,0,0]]],
                 #Disc parameters
                 #type, rin, rout, alpha
                        distribution=[['power_law_shell',10.,1000.,-1]],
                        gridtype='cartesian',
                        rmax= 100.,
                        rho0= 1.5e-19,
                        ngrid= 251,
                #output to be added in SED
                        components = ['total'],
                        **kwargs):
        #Same as in __init__!
        l = locals()
        for key, value in l.items():
            if key == "self": #skip self to avoid possible recursion
                continue
            setattr(self, key, value)
        
        if nchem != len(massfrac):
            print("Number of chemical components does not equal mass fractions")
        
        if np.sum(massfrac) != 1.0:
            print("Mass fraction of all chemical components must equal 1")
        

        #Define opacities for use in model calculation - probably only going to be a single species in most cases, but
        #should make it general for multiple species per component, and multiple components (see above)    
        #Read in opacities - do this in __init__, check they exist during __call__

        #Generating the dust optical constants can be part of __init__, too
        #Convenience function to write dust parameter file '<dust>.params' for Hyperion BHDust calculator (separate program)
        
        f=open(str(dust)+'.params','w')
        self.param_file = str(dust)+'.params'
        f.write(str(dust)+'_'+str(amin[0])+'\n')
        f.write(str(int(fileformat))+'\n')
        f.write(str(np.min(amin))+'\n')
        f.write(str(np.max(amax))+'\n')
        f.write(str(np.max(na))+'\n')
        f.write(str(nang)+'\n')
        f.write(str(nanx)+'\n')
        f.write(str(nchem)+'\n')
        f.write(str(gtd)+'\n')
        f.write(str(self.lmin)+' '+str(self.lmax)+' '+str(self.nl)+'\n')
        f.write(str(self.nproc)+'\n') #parallel processing version of BHMie
        for i in range(0,len(self.optconst)):
            f.write(''+'\n')
            f.write(str(massfrac[i])+'\n')
            f.write(str(density[i])+'\n')
            f.write(str(optconst[i])+'\n')
            f.write(str(disttype[i])+'\n')
            f.write(str(amin[i])+' '+str(amax[i])+' '+str(q[i])+'\n')
            f.close()
            
            self.npars += 5 #3 for size distribution, 1 for material, 1 for mass fraction
            
        print("BHMie dust input file created.")
    
        subprocess.run(['bhmie',self.param_file])

        print("BHMie dust output file created")
        #need a way to iteratively add dust models to the Model object so that they can be called later by name
        self.d = BHDust(str(dust)+'_'+str(amin))
        # self.d.optical_properties.extrapolate_nu(5e7, 5e16)
        self.d.optical_properties.extrapolate_wav(0.95*self.lmin, 1.05*self.lmax)
        self.d.set_lte_emissivities(n_temp=self.nt,temp_min=self.tmin,temp_max=self.tmax)

        print("Read in optical constants.")            
               
        #Add the grid(s) to the model - can specify multiple density grids for different species, physical components
        #Set up spatial grid
        #The grid limits here should be dynamic - set up to fit the largest extent of observations or a 
        #specific field of view or some such
        #maybe have standard models in here - power law shell = 3 parameters (S0, Rin, alpha_out), 
        #                                     power law disc = 5 parameters (S0, Rin, Rout, h, alpha_out),
        #                                     2 power law disc = 5 parameters (S0, R0, h, alpha_in, alpha_out), 
        #                                     Gaussian annulus = 4 parameters (S0, R0, dR, h)
        self.npars += len(self.distribution) - 1

        for distro in self.distribution:
        
            if gridtype == 'cartesian':
                self.x = np.linspace(-1.*self.rmax*units.au, self.rmax*units.au, ngrid)
                self.y = np.linspace(-1.*self.rmax*units.au, self.rmax*units.au, ngrid)
                self.z = np.linspace(-1.*self.rmax*units.au, self.rmax*units.au, ngrid)
                self.model.set_cartesian_grid(self.x,self.y,self.z)
            #Set up density grid
                self.rr = np.sqrt(self.model.grid.gx**2 + self.model.grid.gy**2 + self.model.grid.gz**2)
            if gridtype == 'polar':
                print("Not implemented.")
            if gridtype == 'spherical':
                print("Not implemented")

            print("Spatial grid set up.")

        #define density grid
            self.density = eval(distro[0] + "()")
            self.model.add_density_grid(self.density, self.d)

            print("Density grid set up.")

        #Set central source position, add source(s) to the grid
        for source in sources: #source with temperature -> blackbody
            if source[0] == 'spherical':
                if type(source[-1] == 'list'):
                    self.model.add_spherical_source(luminosity  = source[1] * units.lsun,
                                                    radius      = source[2] * units.rsun,
                                                    mass        = source[3] * units.msun,
                                                    temperature = source[4] * units.K,
                                                    position    = source[5])
                elif type(source[-1] == 'str'):
                    data = np.loadtxt(source[6], dtype=[('wav', float), ('fnu', float)])
        
                    # Convert to nu, fnu from erg/cm2/A/s
                    freq = const.c / (np.array(data['wav'].data[1:])*1e-8)
                    angs = np.array(data['wav'].data[1:])
                    flux = np.array(data['fnu'].data[1:])
                    flux = flux*3.33564095E+04*angs**2
                    nu = []
                    fnu = []
                    for i in range(0,len(freq)-1):
                        if freq[i] != freq[i-1]:
                            nu.append(freq[i])
                            fnu.append(flux[i])
                    
                    nu  = np.array(nu)
                    fnu = np.array(fnu)
        
                    # Set up the source
                    self.model.add_point_source(luminosity  = source[1] * units.lsun,
                                                mass        = source[3] * units.msun,
                                                temperature = source[4] * units.K,
                                                position    = source[5])
                    self.source.luminosity = source[1] * units.lsun # [ergs/s]
                    self.source.spectrum = (nu, fnu)
            
            elif source[0] == 'point':
                if type(source[-1]) == 'list':
                    self.model.add_spherical_source(luminosity  = source[1] * units.lsun,
                                                    radius      = source[2] * units.rsun,
                                                    mass        = source[3] * units.msun,
                                                    temperature = source[4] * units.K,
                                                    position    = source[5])
                elif type(source[-1]) == 'str':
                    #The source spectrum files need to be passed from the data class
                    data = np.loadtxt(source[6], dtype=[('wav', float), ('fnu', float)])
        
                    # Convert to nu, fnu from erg/cm2/A/s
                    freq = const.c / (np.array(data['wav'].data[1:])*1e-8)
                    angs = np.array(data['wav'].data[1:])
                    flux = np.array(data['fnu'].data[1:])
                    flux = flux*3.33564095E+04*angs**2
                    nu = []
                    fnu = []
                    for i in range(0,len(freq)-1):
                        if freq[i] != freq[i-1]:
                            nu.append(freq[i])
                            fnu.append(flux[i])
                    
                    nu  = np.array(nu)
                    fnu = np.array(fnu)
        
                    # Set up the source
                    self.source = self.model.add_point_source(position = source[4])
                    self.source.luminosity = source[1] * units.lsun # [ergs/s]
                    self.source.spectrum = (nu, fnu)
            else:
                print("Source must be point or spherical.")
                
        print("Source(s) set up.")

        self.model.write('HyperionRT_sed.rtin')
        print("Hyperion RT model created.")

        self.model.run('HyperionRT_sed.rtout',mpi=self.useMPI,n_processes=self.nproc,overwrite=True)
        print("Hyperion RT model executed.")

        self.result = ModelOutput('HyperionRT_sed.rtout')
        #forcing code to return total (scattered + emitted) component here - might want to avoid this if self-scattering is not desired, for example
        self.HyperionRTSED = self.result.get_sed(inclination=0, aperture=-1, distance=self.dstar * units.pc,component='total',units='mJy')
        
        self.wave = self.HyperionRTSED.wav
        self.flux = self.HyperionRTSED.val

    def lnprior(self, theta, **kwargs):
        if self.flatprior:
            if (self.lims[0,0] < theta[0] < self.lims[0,1]) and \
               (self.lims[1,0] < theta[1] < self.lims[1,1]) and \
               (self.lims[2,0] < theta[2] < self.lims[2,1]) and \
               (self.lims[3,0] < theta[3] < self.lims[3,1]) and \
                np.sum(10**theta[4:]) <= 1. and np.all(theta[4:] < 0.): 
                return 0
            else:
                return -np.inf
        else:
            raise NotImplementedError()

    def __str__(self, **kwargs):
        raise NotImplementedError()
    #This will provide a summary of the properties of the model
    def __repr__(self, **kwargs):
        raise NotImplementedError()

    def power_law_shell(self):
        density = np.zeros(self.rr)
        density = self.rho0 * (self.rr/self.distribution[1])**self.distribution[3] #* np.exp(-((abs(self.model.grid.gz)/self.rr)**2/scaleheight**2)/2.0)
        density[np.where(self.rr <= self.distribution[1])] = 0.0
        density[np.where(self.rr >= self.distribution[2])] = 0.0

        return density

 
# a test class for carbon-star fitting
class HyperionCStarRTModel(Model):
    '''
    Class to link Ampere and Hyperion and to set up a
    spherically-symmetric RT model for a carbon star.
    '''

    def __init__(self,wavelength, #This is a grid of wavelengths onto which we will interpolate the RT output. 
                                   #This is kinda dangerous, but necessary at least for now.
                 flatprior=True,
                 lims = None,
                 parLabels = None,
                 #RT parameters - use in __init__ and store as self.XXX
                        niter=5,
                        nph_initial=1e4,
                        nph_imging=1e5,
                        nph_rtsrcs=1e5,
                        nph_rtdust=1e5,
                        track_mode='detailed', #one of no/basic/detailed/scatterings
                        raytracing=True,
                        modrndwlk=True,
                        mrw_gamma=2,
                        lmin=0.1,lmax=200.0,nl=101,
                 #Peel photons to get images - use in __init__ and store as self.XXX
                        api_sed=True,
                        api_img=False,
                        view_angles=[45],
                        view_repeats=[45],
                 #Parallel processing - use in __init__ and store as self.XXX
                        useMPI=True,
                        nproc=1,
                        npars=0,
                 #Dust parameters that remain fixed
                        dust="AmC_plus_SiC",fileformat=2,amin=[0.5,0.5],amax=[1000.,1000.],na=[101,101],
                        nang=91,nanx=11,
                        nchem=2,gtd=0,
                        tmin=2.7,tmax=2000.0,nt=101, density=[1.80, 3.22],
                        disttype=['power', 'power'],
                        optconst=["zubko96_ac_acar.optc", "SiC_Pegourie1988.dat"],
                        q=[3.5, 3.5],
                 #Source parameters #(type,lstar,rstar,mstar,tstar,position[x,y,z],spectrum file)
                        stellar_spectrum = 'photosphere_model.csv',
                        stellar_luminosity = 6.165950E+03 * units.solLum,
                        stellar_temperature = 3000 * units.K,
                        stellar_radius = 2.056479E+13 * units.cm,
                        stellar_mass = 2 * units.solMass,
                        stellar_distance = 50.0 * units.kpc,
                        # sources=[['spherical',1.0,1.0,1.0,5784,[0.0,0.0,0.0],'photosphere_model.csv']],
                 #Disc parameters
                        envelope_mass = 6.985718e-6, #in Solar masses
                        envelope_rin = 4.4859, # in stellar radii
                        envelope_rout = 1000.0, # in inner radii
                        envelope_r0 = 4.4859, # in stellar radii
                        envelope_power = -2,
                        envelope_rho0 = 1.41e-18 * units.g / units.cm**3, # only required if mass not given. Write code to compute one given the other
                        envelope_velocity = 10 * units.km / units.s,
                        envelope_ngrid = 251,
                        envelope_gridtype = 'spherical',
                        # distribution=[['power_law_shell',3.0,1000.0,-2]], #type,rin,rout,alpha # define radii in terms of the stellar radius
                        # gridtype='spherical',
                        # rmax= 2000.0, # we need this grid to be logarithmic
                        # rho0= 1.5e-19,
                        # ngrid= 251,
                #output to be added in SED
                        components = ['total'],
                        **kwargs):
        
        #Assign keyword variables to the object
        #Assign all the inputs to __init__ to instance variables with the same name as above
        #this is equivalent to many lines of self.niter = niter etc
        l = locals()
        for key, value in l.items():
            if key == "self": #skip self to avoid possible recursion
                continue
            setattr(self, key, value) #builtin that assigns to an attribute with name = key of self with a value of value
        
        self.model = AnalyticalYSOModel()
        
        # Set up stellar parameters
        if stellar_spectrum is not None:
            nu, fnu = np.loadtxt(stellar_spectrum, delimiter = ',', skiprows = 1, unpack = True)
            self.model.star.spectrum = (nu, fnu)
        elif temperature is not None:
            self.model.star.temperature = stellar_temperature.value
        else: 
            raise ValueError("Either the stellar temperature must be specified or a template spectrum must be specified")
        self.model.star.luminosity = stellar_luminosity.to(units.erg / units.s).value
        self.model.star.mass = stellar_mass.to('g').value
        self.model.star.radius = stellar_radius.to('cm').value
        self.dstar = stellar_distance.to('cm').value
        
        # Set up envelope parameters
        self.envelope_shell = self.model.add_power_law_envelope()
        # Convert to appropriate units and then feed in the magnitudes
        self.envelope_shell.mass = (envelope_mass * units.solMass).to('g').value
        self.envelope_shell.rmin = (envelope_rin * self.model.star.radius)
        self.envelope_shell.rmax = envelope_rout * self.envelope_shell.rmin
        self.envelope_shell.r_0 = envelope_r0 * self.envelope_shell.rmin
        self.envelope_shell.power = envelope_power
        # self.envelope_shell.dust This should be done by Hyperion (see Jonty's code)
        # self.envelope_shell.ngrid = envelope_ngrid # where is this used?
        # self.envelope_shell.gridtype = envelope_gridtype # where is this used?
        
        # A spherical grid with sampling along r but not along theta or phi
        self.model.set_spherical_polar_grid_auto(envelope_ngrid, 1, 1)
        
        self.nSpecies = nchem
        
        #Use raytracing to improve s/n of thermal/source emission
        self.model.set_raytracing(self.raytracing)
        
        #Use the modified random walk
        self.model.set_mrw(self.modrndwlk, gamma=self.mrw_gamma)
        
        #Set up SED for 10 viewing angles
        sed = self.model.add_peeled_images(sed=self.api_sed, image=False)
        sed.set_viewing_angles(self.view_angles,self.view_repeats)
        sed.set_uncertainties(uncertainties = True) #This is currently ignored by ampere
        sed.set_wavelength_range(nl, lmin, lmax)
        sed.set_track_origin(self.track_mode)

        self.sed = sed

        #Set number of photons
        self.model.set_n_photons(initial=self.nph_initial, imaging=self.nph_imging,
                        raytracing_sources=self.nph_rtsrcs, raytracing_dust=self.nph_rtdust)
        
        #Set number of temperature iterations
        self.model.set_n_initial_iterations(self.niter)
        self.model.set_convergence(True, percentile=99.0, absolute=2.0, relative=1.1)

        print("Hyperion RT Model setup for carbon star complete.")
        #self.dictionary of fixed parameters that are needed for modelling in __call__
        
        self.npars = 8
        self.npars_ptform = 8

    #Only give __call__ the numbers that emcee is going to change.
    #labels = ['envelope_mass','envelope_rin','envelope_rout','envelope_r0','stellar_mass','stellar_luminosity']
    def __call__(self,envelope_mass,envelope_rin,envelope_rout,envelope_r0,stellar_mass,stellar_luminosity, *args, **kwargs):
        """
        The parameters passed to *args at the moment are the two
        mass fractions.
        """
        l = locals()
        for key, value in l.items():
            if key == "self": #skip self to avoid possible recursion
                continue
            setattr(self, key, value)
            
        if len(args) == self.nSpecies-1: #emcee-like codes need to end up in this pathway
            args = np.append(np.array(args),1.-np.sum(np.array(args)))
        elif len(args) == self.nSpecies: #codes which use a prior transform need to end up in this branch
            args = args
        
        if self.nchem != len(args):
            print("Number of chemical components does not equal mass fractions")
        
        if np.sum(args) != 1.0:
            print("Total mass fraction of all chemical components must equal 1")
        
        #Define opacities for use in model calculation - probably only going to be a single species in most cases, but
        #should make it general for multiple species per component, and multiple components (see above)    
        #Read in opacities - do this in __init__, check they exist during __call__

        #Generating the dust optical constants can be part of __init__, too
        #Convenience function to write dust parameter file '<dust>.params' for Hyperion BHDust calculator (separate program)
        
        with open(str(self.dust)+'.params','w') as f:
            self.param_file = str(self.dust)+'.params'
            f.write(str(self.dust)+'_'+str(self.amin[0])+'\n')
            f.write(str(int(self.fileformat))+'\n')
            f.write(str(np.min(self.amin))+'\n')
            f.write(str(np.max(self.amax))+'\n')
            f.write(str(np.max(self.na))+'\n')
            f.write(str(self.nang)+'\n')
            f.write(str(self.nanx)+'\n')
            f.write(str(self.nchem)+'\n')
            f.write(str(self.gtd)+'\n')
            f.write(str(self.lmin)+' '+str(self.lmax)+' '+str(self.nl)+'\n')
            f.write(str(self.nproc)+'\n') #parallel processing version of BHMie
            
            print(args,len(args))
            for i in range(0,len(self.optconst)):
                f.write(''+'\n')
                f.write(str(args[i])+'\n')
                f.write(str(self.density[i])+'\n')
                f.write(str(self.optconst[i])+'\n')
                f.write(str(self.disttype[i])+'\n')
                f.write(str(self.amin[i])+' '+str(self.amax[i])+' '+str(self.q[i])+'\n')
            f.close()
            
            self.npars += 5 #3 for size distribution, 1 for material, 1 for mass fraction
            
        print("BHMie dust input file created.")
    
        subprocess.run(['bhmie',self.param_file])

        print("BHMie dust output file created")
        #need a way to iteratively add dust models to the Model object so that they can be called later by name
        self.d = BHDust(str(self.dust)+'_'+str(self.amin[0]))
        self.d.optical_properties.extrapolate_wav(0.05, 1200)#0.95*self.lmin, 1.05*self.lmax)
        # self.d.optical_properties.extrapolate_nu(5e7, 5e16)
        self.d.set_lte_emissivities(n_temp=self.nt,temp_min=self.tmin,temp_max=self.tmax)

        print("Read in optical constants.")   
        self.envelope_shell.dust = self.d
               
        #Add the grid(s) to the model - can specify multiple density grids for different species, physical components
        #Set up spatial grid
        self.distribution = ['power_law_shell',self.envelope_rin,self.envelope_rout,self.envelope_power]
        
        self.npars += len(self.distribution) - 1
        
        self.model = AnalyticalYSOModel()
        
        # Set up stellar parameters
        if self.stellar_spectrum is not None:
            nu, fnu = np.loadtxt(self.stellar_spectrum, delimiter = ',', skiprows = 1, unpack = True)
            self.model.star.spectrum = (nu, fnu)
        elif temperature is not None:
            self.model.star.temperature = self.stellar_temperature.value
        else: 
            raise ValueError("Either the stellar temperature must be specified or a template spectrum must be specified")
        self.model.star.luminosity = (stellar_luminosity*units.solLum).to(units.erg / units.s).value
        self.model.star.mass = (stellar_mass*units.solMass).to('g').value
        self.model.star.radius = self.stellar_radius.value
        #self.dstar = (stellar_distance*units.kpc).to('cm').value
        
        # Set up envelope parameters
        self.envelope_shell = self.model.add_power_law_envelope()
        # Convert to appropriate units and then feed in the magnitudes
        self.envelope_shell.mass = ((10**envelope_mass)*units.solMass).to('g').value
        self.envelope_shell.rmin = (10**envelope_rin) * self.model.star.radius
        self.envelope_shell.rmax = (10**envelope_rout) * self.model.star.radius
        self.envelope_shell.r_0 = (10**envelope_r0) * self.model.star.radius
        self.envelope_shell.power = self.envelope_power
        self.envelope_shell.dust = self.d #This should be done by Hyperion (see Jonty's code)
        # self.envelope_shell.ngrid = envelope_ngrid # where is this used?
        # self.envelope_shell.gridtype = envelope_gridtype # where is this used?
        
        # A spherical grid with sampling along r but not along theta or phi
        self.model.set_spherical_polar_grid_auto(self.envelope_ngrid, 1, 1)
        
        self.nSpecies = self.nchem
        
        #Use raytracing to improve s/n of thermal/source emission
        self.model.set_raytracing(self.raytracing)
        
        #Use the modified random walk
        self.model.set_mrw(self.modrndwlk, gamma=self.mrw_gamma)
        
        #Set up SED for 10 viewing angles
        sed = self.model.add_peeled_images(sed=self.api_sed, image=False)
        sed.set_viewing_angles(self.view_angles,self.view_repeats)
        sed.set_uncertainties(uncertainties = True) #This is currently ignored by ampere
        sed.set_wavelength_range(self.nl, self.lmin, self.lmax)
        sed.set_track_origin(self.track_mode)

        self.sed = sed

        #Set number of photons
        self.model.set_n_photons(initial=self.nph_initial, imaging=self.nph_imging,
                        raytracing_sources=self.nph_rtsrcs, raytracing_dust=self.nph_rtdust)
        
        #Set number of temperature iterations
        self.model.set_n_initial_iterations(self.niter)
        self.model.set_convergence(True, percentile=99.0, absolute=2.0, relative=1.1)
                
        print("Source(s) set up.")       

        self.model.write('HyperionRT_sed.rtin')
        print("Hyperion RT model created.")

        self.model.run('HyperionRT_sed.rtout',mpi=self.useMPI,n_processes=self.nproc,overwrite=True)
        print("Hyperion RT model executed.")

        self.result = ModelOutput('HyperionRT_sed.rtout')
        #forcing code to return total (scattered + emitted) component here - might want to avoid this if self-scattering is not desired, for example
        self.HyperionRTSED = self.result.get_sed(inclination=0, aperture=-1, distance=self.dstar,component='total',units='Jy')
        
        self.wave = copy.copy(self.HyperionRTSED.wav)
        self.flux = copy.copy(self.HyperionRTSED.val)
        
        #Interpolate model fluxes onto fixed wavelength grid
        self.modelFlux = np.interp(self.wavelength, np.flip(self.wave), np.flip(self.flux))
        
        return {"spectrum":{"wavelength":self.wavelength, "flux": self.modelFlux}}

    def lnprior(self, theta, **kwargs): #theta will only have nSpecies-1 entries for this case
        
        #Some of our random variables are relative abundances - this means that they lie on a simplex and we have to sample from a dirichlet distribution to get it right:
        if theta[-self.nSpecies::].all() > 0. and theta[-self.nSpecies::].all() < 1.:
            try:
                lp = dirichlet((1.,)*self.nSpecies).logpdf(theta) #When you have a set of random variates whose sum = 1 by definition, the correct prior is a dirichlet distribution (which is a generalisation of the beta distribution).
            except ValueError:
                return -np.inf
        else:
            return -np.inf
        #print(theta)
        #print(self.lims)
        #The rest of our parameters are uniform random variables, so we need to draw those appropriately
        lp+=np.sum([uniform(self.lims[i,0], self.lims[i,1] - self.lims[i,0]).logpdf(theta[i]) for i in range(self.npars_ptform - self.nSpecies)])
        return lp
    
    def prior_transform(self, u, **kwargs): #In this case, u will contain nSpecies entries, rather than nSpecies-1 as above
        theta = np.zeros_like(u)
        #First we deal with all the non-abundance parameters - for simplicity, let's say they're linearly distributed in the range of the limits
        theta[:-self.nSpecies] = self.lims[:-self.nSpecies,0] + np.asarray(u[:-self.nSpecies]) * (self.lims[:-self.nSpecies,1] - self.lims[:-self.nSpecies,0])
        #Now we come to the abundances, which should be dirichlet-distributed:
        gamma_quantiles = -np.log(u[-self.nSpecies:])
        theta[-self.nSpecies:] = gamma_quantiles/gamma_quantiles.sum()
        
        return theta
    
    #def __fuck__(self):
    #    return  str(self.__class__) + '\n'+ '\n'.join(('{} = {}'.format(item, type(self.__dict__[item]).__mro__) for item in self.__dict__))
    def __str__(self, **kwargs):
        return  str(self.__class__) + '\n' + '\n'.join((str(item) + ' = ' + str(self.__dict__[item]) for item in sorted(self.__dict__)))
    #This will provide a summary of the properties of the model
    def __repr__(self, **kwargs):
        from pprint import pformat
        return "<" + type(self).__name__ + "> " + pformat(vars(self), indent=4, width=1)

    def power_law_shell(self):
        density = np.zeros(self.rr)
        density = self.rho0 * (self.rr/self.distribution[1])**self.distribution[3] #* np.exp(-((abs(self.model.grid.gz)/self.rr)**2/scaleheight**2)/2.0)
        density[np.where(self.rr <= self.distribution[1])] = 0.0
        density[np.where(self.rr >= self.distribution[2])] = 0.0
        return density
