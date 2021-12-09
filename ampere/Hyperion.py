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

    def __init__(self,flatprior=True,opacityFileList=None,parameters=par,**kwargs):

        self.components = components #number of sub-models to run to generate complete model
        self.distribution = distribution #dust density distribution + its properties, somehow?
        #maybe have standard models in here - power law shell = 3 parameters (S0, Rin, alpha_out), 
        #                                     power law disc = 5 parameters (S0, Rin, Rout, h, alpha_out),
        #                                     2 power law disc = 5 parameters (S0, R0, h, alpha_in, alpha_out), 
        #                                     Gaussian annulus = 4 parameters (S0, R0, dR, h)
        self.npars += distribution.nParameters

        self.flatprior = flatprior
        #Define opacities for use in model calculation - probably only going to be a single species in most cases, but
        #should make it general for multiple species per component, and multiple components (see above)
        import os
        from scipy import interpolate
        opacityDirectory = os.path.dirname(__file__)+'/Opacities/'
        opacityFileList = np.array(opacityFileList)
        #opacityFileList = np.array(opacityFileList)[['.q' in zio for zio in opacityFileList]] # Only files ending in .q are valid (for now)
        nSpecies = opacityFileList.__len__()
        #print(opacityFileList,nSpecies)
        opacity_array = np.zeros((wavelengths.__len__(), nSpecies))
        for j in range(nSpecies):
            tempData = np.loadtxt(opacityDirectory + opacityFileList[j], comments = '#')
            tempWl = tempData[:, 0]
            tempOpac = tempData[:, 1]
            f = interpolate.interp1d(tempWl, tempOpac, assume_sorted = False)
            opacity_array[:,j] = f(self.wavelength)
        self.opacity_array = opacity_array
        self.nSpecies = nSpecies
        self.npars += nSpecies - 1 + 4       
    
        model = Model()

        return model

    def HyperionRTAddDust(self,par,**kwargs):
        #Set up dust model for Hyperion
        # Initalize the model
        
        from astropy.units import au
        #Add the grid to the model
        #Set up spatial grid
        #The grid limits here should be dynamic - set up to fit the largest extent of observations or a 
        #specific field of view or some such
        x = np.linspace(-1.5*par.get('rmax',100.*u.au), 1.5*par.get('rmax',100.*u.au), par.get('ngrid',257))
        y = np.linspace(-1.5*par.get('rmax',100.*u.au), 1.5*par.get('rmax',100.*u.au), par.get('ngrid',257))
        z = np.linspace(-1.5*par.get('rmax',100.*u.au), 1.5*par.get('rmax',100.*u.au), par.get('ngrid',257))
        model.set_cartesian_grid(x,y,z)

        print("Spatial grid set up.")

        #Set up density grid
        rho0 = par.get('rho0',1.5e-19)
        alpha_in = par.get('alpha_in',-5.)
        alpha_out = par.get('alpha_out',5.)
        scaleheight= par.get('scaleheight',0.1)
        rmin = par.get('rmin',70.*u.au)
        rmax = par.get('rmax',100.*u.au)
        rmid = (rmax + rmin) / 2
        
        rr = np.sqrt(model.grid.gx ** 2 + model.grid.gy ** 2 + model.grid.gz ** 2)
        
        #define density grid
        density = eval(distribution + "(model,par)") #np.zeros(model.grid.shape)
        
        model.add_density_grid(density, d)

        print("Density grid set up.")

    def HyperionRTAddSource(self,par,**kwargs):
        #Add the [central] source to the model - should be either spherical or point source
        #Can also specify spectrum as part of the inputs...    
        if par.any('source'):
            #Set central source position
            model.add_spherical_source(luminosity  = par.get('lsource',1.0*lsun),
                                        radius      = par.get('rsource',1.0*rsun),
                                        mass        = par.get('msource',1.0*msun),
                                        temperature = par.get('tsource',5784.0),
                                        position    = par.get('position',(0.0,0.0,0.0))

        print("Source set up.")


    def HyperionRTSetup(self,par,**kwargs):
        #This is a catch all for the bits and pieces that actually run the model
        #The number of photons, ray tracing, wavelengths, all that jazz

        
        #Setup wavelengths
        

        #Use raytracing to improve s/n of thermal/source emission
        model.set_raytracing(par.get('raytracing',True))
        
        #Use the modified random walk
        model.set_mrw(par.get('modrndwlk',True), gamma=par.get('mrw_gamma',2))
        
        #Set up SED for 10 viewing angles
        model.add_peeled_images(sed=par.get('api_sed',True), image=par.get('api_img',False))
        model.set_viewing_angles(np.linspace(0., 90., 10), np.repeat(45., 10))
        model.set_wavelength_range(par.get('nl',101), par.get('lmin',0.1), par.get('lmax',1000.))
        model.set_track_origin('basic')
        
        #Set number of photons
        model.set_n_photons(initial=par.get('nph_initial',1e4), imaging=par.get('nph_imging',1e5),
                        raytracing_sources=par.get('nph_rtsrcs',1e5), raytracing_dust=par.get('nph_rtdust',1e5))
        
        #Set number of temperature iterations
        model.set_n_initial_iterations(par.get('niter',5))

        print("Model setup complete.")

    def __call__(self, *args, **kwargs):

        model.write('HyperionRT_sed.rtin')
        print("Hyperion RT model created.")
        model.run('HyperionRT_sed.rtout', mpi=True,n_processes=6,overwrite=True)
        print("Hyperion RT model executed.")

        result = ModelOutput('HyperionRT_sed.rtout')

        #Save both the model and the results
        self.model = model
        self.result = result

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

    def __repr__(self, **kwargs):
        raise NotImplementedError()

    def inputFile(self, **kwargs):
        #Create a dictionary with which to setup and run Hyperion RT models
        #Dust parameters
        par =   {'dust':'"astrosilicate"',
                 'format':2,
                 'size':0.5,
                 'amin':0.5,
                 'amax':1000.,
                 'na':101,
                 'nang':91,
                 'nanx':11,
                 'nchem':1,
                 'gtd':0,
                 'lmin':0.1,
                 'lmax':1000.0,
                 'nl':101,
                 'massfrac':1.0,
                 'rho0':1.5e-19,
                 'optconst':'"silicate_d03.lnk"',
                 'disttype':'power',
                 'q':3.5,
                 #Source parameters
                 'lsource': 1.0*u.lsun,
                 'rsource': 1.0*u.rsun,
                 'msource': 1.0*u.msun,
                 'tsource': 5784.,
                 'position':[0.0,0.0,0.0],
                 #Disc parameters
                 'rmin': 70.*u.au,
                 'rmax': 100.*u.au,
                 'alpha_in': -5.,
                 'alpha_out': 5.,
                 'scaleheight': 0.1,
                 #RT parameters
                 'niter':5,
                 'nph_initial':1e4,
                 'nph_imging':1e5,
                 'nph_rtsrcs':1e5,
                 'nph_rtdust':1e5,
                 #Peel photons to get images
                 'api_sed':True,
                 'api_img':False,
                 }

        return par

    def HyperionRTPlot(self):
        #
        # Insert code here to call the plotting routines for the respective data classes
        # and Hyperion model.
        #


#convenience function to write dust parameter file '<dust>.params' for Hyperion BHDust calculator (separate program)
    def HyperionRTBHMie(self,par):
        f=open(par.get('dust','astrosilicate')+'.params','w')
        f.write('"'+par['dust']+'_'+str(par['size'])+'"'+'\n')
        f.write(str(par['format'])+'\n')
        f.write(str(par['amin'])+'\n')
        f.write(str(par['amax'])+'\n')
        f.write(str(par['na'])+'\n')
        f.write(str(par['nang'])+'\n')
        f.write(str(par['nanx'])+'\n')
        f.write(str(par['nchem'])+'\n')
        f.write(str(par['gtd'])+'\n')
        f.write(str(par['lmin'])+' '+str(par['lmax'])+' '+str(par['nl'])+'\n')
        f.write(''+'\n')
        f.write(str(par['massfrac'])+'\n')
        f.write(str(par['rho0'])+'\n')
        f.write(str(par['optconst'])+'\n')
        f.write(str(par['disttype'])+'\n')
        f.write(str(par['amin'])+' '+str(par['amax'])+' '+str(par['q'])+'\n')
        f.close()

        print("BHMie dust input file created.")

        import subprocess
        subprocess.run(['bhmie',param_file])

        print("BHMie dust output file created")