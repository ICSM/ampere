#
# Ampere Stuff
#

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
                        nph_initial=1e4,
                        nph_imging=1e5,
                        nph_rtsrcs=1e5,
                        nph_rtdust=1e5,
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
                        nproc=1,
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
            
        
        ##Assign keyword variables to the object - there has to be a better way
        ##to do this!
        #self.niter = niter
        #self.nph_initial = nph_initial
        #self.nph_imging = nph_imging
        #self.nph_rtsrcs = nph_rtsrcs
        #self.nph_rtdust = nph_rtdust
        #self.track_mode = track_mode
        #self.raytracing = raytracing
        #self.modrndwlk = modrndwlk
        #self.mrw_gamma = mrw_gamma
        #self.api_sed = api_sed
        #self.api_img = api_img
        #self.view_angles = view_angles
        #self.view_repeats = view_repeats
        #self.useMPI = useMPI
        #self.nproc = nproc
        #self.lmin=lmin
        #self.lmax=lmax
        #self.nl=nl
        
        
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
                        sources=[['spherical',1.0,1.0,1.0,5784,[0.0,0.0,0.0]]], #(type,lstar,rstar,mstar,tstar,position[x,y,z],spectrum file)
                 #Disc parameters
                        distribution=[['power_law_shell',30,70,-1]], #type,rin,rout,alpha
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
        
        ##Assign keyword variables to the object - there has to be a better way
        ##to do this!
        #self.dust=dust
        #self.fileformat=fileformat
        #self.amin=amin
        #self.amax=amax
        #self.na=na
        #self.nang=nang
        #self.nanx=nanx
        #self.nchem=nchem
        #self.gtd=gtd
        #self.tmin=tmin
        #self.tmax=tmax
        #self.nt=nt
        #self.massfrac=massfrac
        #self.density=density
        #self.optconst=optconst
        #self.disttype=disttype
        #self.q=q
#Source parameters
        #self.sources=sources #(type,lstar,rstar,mstar,tstar,position[x,y,z],spectrum file)
#Disc parameters
        #self.distribution=distribution #type,rin,rout,alpha
        #self.gridtype=gridtype
        #self.rmax=rmax
        #self.rho0=rho0
        #self.ngrid=ngrid
        #self.components = components
        
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

    def __init__(self,flatprior=True,                 
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
                        lmin=0.1,lmax=1000.0,nl=101,
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
                        stellar_luminosity = 3000 * units.solLum,
                        stellar_temperature = 3000 * units.K,
                        stellar_radius = 1 * units.au,
                        stellar_mass = 1 * units.solMass,
                        # sources=[['spherical',1.0,1.0,1.0,5784,[0.0,0.0,0.0],'photosphere_model.csv']],
                 #Disc parameters
                        envelope_mass = 6.985718e-6 * units.solMass,
                        envelope_rin = 4.4859,
                        envelope_rout = 1000.0,
                        envelope_r0 = 4.4859,
                        envelope_power = -2,
                        envelope_rho0 = 1.5e-19 * units.g / units.cm**3, # only required if mass not given. Write code to compute one given the other
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
        
        # Set up envelope parameters
        self.envelope_shell = self.model.add_power_law_envelope()
        # Convert to appropriate units and then feed in the magnitudes
        self.envelope_shell.mass = envelope_mass.to('g').value
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

    #Only give __call__ the numbers that emcee is going to change.
    def __call__(self, *args, **kwargs):
        """
        The only parameters passed to *args at the moment are the two
        mass fractions.
        """
        #Same as in __init__!
        l = locals()
        for key, value in l.items():
            if key == "self": #skip self to avoid possible recursion
                continue
            setattr(self, key, value)
            
        if len(args) == self.nSpecies-1: #emcee-like codes need to end up in this pathway
            args = np.append(np.array(args),1.-np.sum(np.array(args)))
        elif len(args) == self.nSpecies: #codes which use a prior transform need to end up in this branch
            args = args
        
        if nchem != len(args):
            print("Number of chemical components does not equal mass fractions")
        
        if np.sum(args) != 1.0:
            print("Total mass fraction of all chemical components must equal 1")
        

        #Define opacities for use in model calculation - probably only going to be a single species in most cases, but
        #should make it general for multiple species per component, and multiple components (see above)    
        #Read in opacities - do this in __init__, check they exist during __call__

        #Generating the dust optical constants can be part of __init__, too
        #Convenience function to write dust parameter file '<dust>.params' for Hyperion BHDust calculator (separate program)
        
        f=open(str(self.dust)+'.params','w')
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
        self.d.optical_properties.extrapolate_wav(0.95*self.lmin, 1.05*self.lmax)
        # self.d.optical_properties.extrapolate_nu(5e7, 5e16)
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
        
        # for distro in self.distribution:
        # 
        #     if gridtype == 'cartesian':
        #         self.x = np.linspace(-1.*self.rmax*units.au, self.rmax*units.au, self.ngrid)
        #         self.y = np.linspace(-1.*self.rmax*units.au, self.rmax*units.au, self.ngrid)
        #         self.z = np.linspace(-1.*self.rmax*units.au, self.rmax*units.au, self.ngrid)
        #         self.model.set_cartesian_grid(self.x,self.y,self.z)
        #     #Set up density grid
        #         self.rr = np.sqrt(self.model.grid.gx**2 + self.model.grid.gy**2 + self.model.grid.gz**2)
        #     if gridtype == 'polar':
        #         print("Not implemented.")
        #     if gridtype == 'spherical':
        #         print("Not implemented")
        #     print("Spatial grid set up.")

        # #define density grid
        #     self.density = eval(distro[0] + "()")
        #     self.model.add_density_grid(self.density, self.d)
        #     print("Density grid set up.")

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
        self.HyperionRTSED = self.result.get_sed(inclination=0, aperture=-1, distance=self.dstar * units.pc,component='total',units='Jy')
        
        self.wave = self.HyperionRTSED.wav
        self.flux = self.HyperionRTSED.val

    def lnprior(self, theta, **kwargs): #theta will only have nSpecies-1 entries for this case
        
        if theta > 0. and theta < 1.:
            try:
                return dirichlet((1.,)*self.nSpecies).logpdf(theta) #When you have a set of random variates whose sum = 1 by definition, the correct prior is a dirichlet distribution (which is a generalisation of the beta distribution).
            except ValueError:
                return -np.inf
        else:
            return -np.inf
        #if self.flatprior:
        #    if (self.lims[0,0] < theta[0] < self.lims[0,1]) and \
        #       (self.lims[1,0] < theta[1] < self.lims[1,1]) and \
        #       (self.lims[2,0] < theta[2] < self.lims[2,1]) and \
        #       (self.lims[3,0] < theta[3] < self.lims[3,1]) and \
        #        np.sum(10**theta[4:]) <= 1. and np.all(theta[4:] < 0.): 
        #        return 0
        #    else:
        #        return -np.inf
        #else:
        #    raise NotImplementedError()
        
    def prior_transform(self, u, **kwargs): #In this case, u will contain nSpecies entries, rather than nSpecies-1 as above
        gamma_quantiles = -np.log(u)
        theta = gamma_quantiles/gamma_quantiles.sum()
        return theta

    
    def __str__(self, **kwargs):
        raise NotImplementedError()
    #This will provide a summary of the properties of the model
    def __repr__(self, **kwargs):
        raise NotImplementedError()

    # def power_law_shell(self):
    #     density = np.zeros(self.rr)
    #     density = self.rho0 * (self.rr/self.distribution[1])**self.distribution[3] #* np.exp(-((abs(self.model.grid.gz)/self.rr)**2/scaleheight**2)/2.0)
    #     density[np.where(self.rr <= self.distribution[1])] = 0.0
    #     density[np.where(self.rr >= self.distribution[2])] = 0.0
    #     return density
