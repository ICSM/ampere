#
# Ampere Stuff
#

import numpy as np
from astropy import constants as const
from astropy import units as u
from astropy.modeling import blackbody
from .models import AnalyticalModel

# 
# Hyperion Stuff
# 

import hyperion
from hyperion.model import Model
from hyperion.model import ModelOutput
from hyperion.dust import BHDust
#from hyperion.util.constants import lsun, rsun, au, msun, sigma, pc, pi #get these from astropy units 

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
                        track_mode='origin',
                        raytracing=True,
                        modrndwlk=True,
                        mrw_gamma=2,
                 #Peel photons to get images - use in __init__ and store as self.XXX
                        api_sed=True,
                        api_img=False,
                        view_angles=np.linspace(0., 90., 10),
                        view_repeats=np.repeat(45., 10),
                 #Parallel processing - use in __init__ and store as self.XXX
                        useMPI=True,
                        nproc=1,
                        **kwargs):
        
        #Assign keyword variables to the object - there has to be a better way
        #to do this!
        self.niter = niter
        self.nph_initial = nph_initial
        self.nph_imging = nph_imging
        self.nph_rtsrcs = nph_rtsrcs
        self.nph_rtdust = nph_rtdust
        self.track_mode = track_mode
        self.raytracing = raytracing
        self.modrndwlk = modrndwlk
        self.mrw_gamma = mrw_gamma
        self.api_sed = api_sed
        self.api_img = api_img
        self.view_angles = view_angles
        self.view_repeats = view_repeats
        self.useMPI = useMPI
        self.nproc = nproc

        self.flatprior = flatprior
        self.model = Model()
        
        #Use raytracing to improve s/n of thermal/source emission
        self.model.set_raytracing(self.raytracing)
        
        #Use the modified random walk
        self.model.set_mrw(self.modrndwlk, gamma=self.mrw_gamma)
        
        #Set up SED for 10 viewing angles
        self.model.add_peeled_images(sed=self.api_sed, image=self.api_img)
        self.model.set_viewing_angles(self.view_angles,self.view_repeats)
        self.model.set_wavelength_range(nl, lmin, lmax)
        self.model.set_track_origin(self.track_mode)
        
        #Set number of photons
        self.model.set_n_photons(initial=self.nph_initial, imaging=self.nph_imging,
                        raytracing_sources=self.nph_rtsrcs, raytracing_dust=self.nph_rtdust)
        
        #Set number of temperature iterations
        self.model.set_n_initial_iterations(self.niter)

        print("Hyperion RT Model setup complete.")
        #self.dictionary of fixed parameters that are needed for modelling in __call__

    #Only give __call__ the numbers that emcee is going to change.
    def __call__(self,dust="astrosilicate",fileformat=2,amin=0.5,amax=1000.,na=101,
                        nang=91,nanx=11,
                        nchem=1,gtd=0,
                        lmin=0.1,lmax=1000.0,nl=101,
                        tmin=2.7,tmax=1000.0,nt=101,
                        massfrac=1.0,
                        rho0=1.5e-19,
                        optconst="silicate_d03.lnk",
                        q=3.5,
                 #Source parameters
                        sources=[['spherical',1.0,1.0,1.0,5784,[0.0,0.0,0.0]]], #(type,lstar,rstar,mstar,tstar,position[x,y,z],spectrum file)
                 #Disc parameters
                        distribution=['example',30,70,-1], #type,rin,rout,alpha
                        gridtype='cartesian',
                        rmax= 100.,
                        ngrid= 251,
                        **kwargs):
        
        #Assign keyword variables to the object - there has to be a better way
        #to do this!
        self.dust=dust
        self.fileformat=fileformat
        self.amin=amin
        self.amax=amax
        self.na=na
        self.nang=nang
        self.nanx=nanx
        self.nchem=nchem
        self.gtd=gtd
        self.lmin=lmin
        self.lmax=lmax
        self.nl=nl
        self.tmin=tmin
        self.tmax=tmax
        self.nt=nt
        self.massfrac=massfrac
        self.rho0=rho0
        self.optconst=optconst
        self.q=q
#Source parameters
        self.sources=sources #(type,lstar,rstar,mstar,tstar,position[x,y,z],spectrum file)
#Disc parameters
        self.distribution=distribution #type,rin,rout,alpha
        self.gridtype=gridtype
        self.rmax=rmax
        self.ngrid=ngrid
        try:
            self.components = len(massfrac) #number of sub-models to run to generate complete model
        except:
            self.components = 1
        #maybe have standard models in here - power law shell = 3 parameters (S0, Rin, alpha_out), 
        #                                     power law disc = 5 parameters (S0, Rin, Rout, h, alpha_out),
        #                                     2 power law disc = 5 parameters (S0, R0, h, alpha_in, alpha_out), 
        #                                     Gaussian annulus = 4 parameters (S0, R0, dR, h)
        self.npars += len(distribution) - 1

        #Define opacities for use in model calculation - probably only going to be a single species in most cases, but
        #should make it general for multiple species per component, and multiple components (see above)    
        #Read in opacities - do this in __init__, check they exist during __call__
        import os
        from scipy import interpolate
        opacityDirectory = os.path.dirname(__file__)+'/Opacities/'
        opacityFileList = np.array(opacityFileList)
        #opacityFileList = np.array(opacityFileList)[['.q' in zio for zio in opacityFileList]] # Only files ending in .q are valid (for now)
        nSpecies = opacityFileList.__len__()
        #print(opacityFileList,nSpecies)
        opacity_array = np.zeros((wavelengths.__len__(), nSpecies))
        for j in range(0,nchem):
            tempData = np.loadtxt(opacityDirectory + opacityFileList[j], comments = '#')
            tempWl = tempData[:, 0]
            tempOpac = tempData[:, 1]
            f = interpolate.interp1d(tempWl, tempOpac, assume_sorted = False)
            opacity_array[:,j] = f(self.wavelength)
        self.opacity_array = opacity_array
        self.nSpecies = nchem
        self.npars += nSpecies - 1 + 4       

        #Generating the dust optical constants can be part of __init__, too
        #Convenience function to write dust parameter file '<dust>.params' for Hyperion BHDust calculator (separate program)
        f=open(str(dust)+'.params','w')
        f.write(str(dust)+'_'+str(amin)+'\n')
        f.write(str(int(fileformat))+'\n')
        f.write(str(amin)+'\n')
        f.write(str(amax)+'\n')
        f.write(str(na)+'\n')
        f.write(str(nang)+'\n')
        f.write(str(nanx)+'\n')
        f.write(str(nchem)+'\n')
        f.write(str(gtd)+'\n')
        f.write(str(lmin)+' '+str(lmax)+' '+str(nl)+'\n')
        f.write(str(self.nproc)+'\n') #parallel processing version of BHMie
        f.write(''+'\n')
        f.write(str(massfrac)+'\n')
        f.write(str(rho0)+'\n')
        f.write(str(optconst)+'\n')
        f.write(str(self.disttype)+'\n')
        f.write(str(amin)+' '+str(amax)+' '+str(q)+'\n')
        f.close()

        print("BHMie dust input file created.")

        import subprocess
        subprocess.run(['bhmie',self.param_file])

        print("BHMie dust output file created")

        self.d = BHDust(dust+'_'+str(amin))
        self.d.optical_properties.extrapolate_nu(5e7, 5e16)
        self.d.set_lte_emissivities(n_temp=self.nt,temp_min=self.tmin,temp_max=self.tmax)

        print("Read in optical constants.")
        
        #Add the grid to the model
        #Set up spatial grid
        #The grid limits here should be dynamic - set up to fit the largest extent of observations or a 
        #specific field of view or some such
        if gridtype == 'cartesian':
            self.x = np.linspace(-1.*rmax*u.au, rmax*u.au, ngrid)
            self.y = np.linspace(-1.*rmax*u.au, rmax*u.au, ngrid)
            self.z = np.linspace(-1.*rmax*u.au, rmax*u.au, ngrid)
            self.model.set_cartesian_grid(self.x,self.y,self.z)
        #Set up density grid
            self.rr = np.sqrt(self.model.grid.gx**2 + self.model.grid.gy**2 + self.model.grid.gz**2)
        if gridtype == 'polar':
            print("Not implemented.")
        if gridtype == 'spherical':
            print("Not implemented")

        print("Spatial grid set up.")

        #define density grid
        self.density = eval(self.distribution[0] + "()")
        self.model.add_density_grid(self.density, self.d)

        print("Density grid set up.")

        #Set central source position, add source(s) to the grid
        for source in sources: #source with temperature -> blackbody
            if source[0] == 'spherical':
                if type(source[-1] == 'list'):
                    self.model.add_spherical_source(luminosity  = source[1] * u.lsun,
                                                    radius      = source[2] * u.rsun,
                                                    mass        = source[3] * u.msun,
                                                    temperature = source[4] * u.K,
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
                    self.model.add_point_source(luminosity  = source[1] * u.lsun,
                                                mass        = source[3] * u.msun,
                                                temperature = source[4] * u.K,
                                                position    = source[5])
                    self.source.luminosity = source[1] * u.lsun # [ergs/s]
                    self.source.spectrum = (nu, fnu)
            
            elif source[0] == 'point'
                if type(source[-1]) == 'list':
                    self.model.add_spherical_source(luminosity  = source[1] * u.lsun,
                                                    radius      = source[2] * u.rsun,
                                                    mass        = source[3] * u.msun,
                                                    temperature = source[4] * u.K,
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
                    self.source.luminosity = source[1] * u.lsun # [ergs/s]
                    self.source.spectrum = (nu, fnu)
            else:
                print("Source must be point or spherical.")
                
        print("Source(s) set up.")

        self.model.write('HyperionRT_sed.rtin')
        print("Hyperion RT model created.")

        self.model.run('HyperionRT_sed.rtout',mpi=self.useMPI,n_processes=self.nproc,overwrite=True)
        print("Hyperion RT model executed.")

        self.result = ModelOutput('HyperionRT_sed.rtout')

        #Extract fluxes and wavelengths into objects that can be passed to fitting or plotting
        self.HyperionRTWave = somethingelse
        self.HyperionRTFlux = something


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

    def example(self):
        density = np.zeros(self.rr)
        density = self.rho0 * (self.rr/self.distribution[1])**self.distribution[3] #* np.exp(-((abs(self.model.grid.gz)/self.rr)**2/scaleheight**2)/2.0)
        density[np.where(self.rr <= self.distribution[1])] = 0.0
        density[np.where(self.rr >= self.distribution[2])] = 0.0

        return density