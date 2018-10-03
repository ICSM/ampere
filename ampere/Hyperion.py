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

    def __init__(self, par, **kwargs):

        #Set up dust model for Hyperion
        # Initalize the model
        m = Model()
        
        #Set central source position
        m.add_spherical_source(luminosity  = par.get('lsource',1.0*lsun),
                               radius      = par.get('rsource',1.0*rsun),
                               mass        = par.get('msource',1.0*msun),
                               temperature = par.get('tsource',5784.0),
                               position    = par.get('position',(0.0,0.0,0.0))
        
        # Use raytracing to improve s/n of thermal/source emission
        m.set_raytracing(par.get('raytracing',True))
        
        # Use the modified random walk
        m.set_mrw(par.get('modrndwlk',True), gamma=par.get('mrw_gamma',2))

        # Set up spatial grid
        x = np.linspace(-1.5*par.get('rmax',100.*u.au), 1.5*par.get('rmax',100.*u.au), par.get('ngrid',257))
        y = np.linspace(-1.5*par.get('rmax',100.*u.au), 1.5*par.get('rmax',100.*u.au), par.get('ngrid',257))
        z = np.linspace(-1.5*par.get('rmax',100.*u.au), 1.5*par.get('rmax',100.*u.au), par.get('ngrid',257))
        m.set_cartesian_grid(x,y,z)

        print("Spatial grid set up.")

        #Set up density grid
        rho0 = par.get('rho0',1.5e-19)
        alpha_in = par.get('alpha_in',-5.)
        alpha_out = par.get('alpha_out',5.)
        scaleheight= par.get('scaleheight',0.1)
        rmin = par.get('rmin',70.*u.au)
        rmax = par.get('rmax',100.*u.au)
        rmid = (rmax + rmin) / 2
        rr = np.sqrt(m.grid.gx ** 2 + m.grid.gy ** 2 + m.grid.gz ** 2)
        
        #define density grid
        density = np.zeros(m.grid.shape)
        density = rho0 * ( (rr/rmid)**(2.*alpha_in) + (rr/rmid)**(2.*alpha_out) )**(-0.5) * np.exp(-((abs(m.grid.gz)/rr)**2/scaleheight**2)/2.0)
        m.add_density_grid(density, d)

        print("Density grid set up.")
        
        # Set up SED for 10 viewing angles
        sed = m.add_peeled_images(sed=par.get('api_sed',False), image=par.get('api_img',False))
        sed.set_viewing_angles(np.linspace(0., 90., 10), np.repeat(45., 10))
        sed.set_wavelength_range(par.get('nl',101), par.get('lmin',0.1), par.get('lmax',1000.))
        sed.set_track_origin('basic')
        
        # Set number of photons
        m.set_n_photons(initial=par.get('nph_initial',1e4), imaging=par.get('nph_imging',1e5),
                        raytracing_sources=par.get('nph_rtsrcs',1e5), raytracing_dust=par.get('nph_rtdust',1e5))
        
        # Set number of temperature iterations
        m.set_n_initial_iterations(par.get('niter',5))
        
        # Write out file
        m.write('HyperionRT_sed.rtin')
        print("Hyperion RT model created.")

    def __call__(self, *args, **kwargs):
        m.run('HyperionRT_sed.rtout', mpi=True,n_processes=6,overwrite=True)
        print("Hyperion RT model executed.")

        m = ModelOutput('HyperionRT_sed.rtout')

        self.modelFlux = sed.value

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
                 'optconst':'"silicated03.lnk"',
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
                 'api_sed':False,
                 'api_img':False,
                 }
#convenience function to plot the SED of the Hyperion RT output
    def plot_sed(par):
        #Set up figure
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        
        #Read in Hyperion model output
        try:
            m = ModelOutput('HyperionRT_sed.rtout')

            #Total SED
            sed = m.get_sed(inclination=0, aperture=-1, distance=par['distance'],
                                    component='total',units='Jy')        
            ax.loglog(sed.wav, sed.val, color='black', lw=3, alpha=0.5)
            
            # Direct stellar photons
            sed = m.get_sed(inclination=0, aperture=-1, distance=par['distance'],
                                   component='source_emit',units='Jy')
            ax.loglog(sed.wav, sed.val, color='blue')
            # Scattered stellar photons
            sed = m.get_sed(inclination=0, aperture=-1, distance=par['distance'],
                                   component='source_scat',units='Jy')
            ax.loglog(sed.wav, sed.val, color='teal')
            # Direct dust photons
            sed = m.get_sed(inclination=0, aperture=-1, distance=par['distance'],
                                   component='dust_emit',units='Jy')
            ax.loglog(sed.wav, sed.val, color='red')        
            # Scattered dust photons
            sed = m.get_sed(inclination=0, aperture=-1, distance=par['distance'],
                                   component='dust_scat',units='Jy')
            ax.loglog(sed.wav, sed.val, color='orange')

            ax.set_xlabel(r'$\lambda$ [$\mu$m]')
            ax.set_ylabel(r'Flux [Jy]')
            ax.set_xlim(par.get('lmin',0.1), par.get('lmax',1000.))
            ax.set_ylim(1e-6*np.max(sed.value),1.1*np.max(sed.value))
            fig.savefig('HyperionRT_sed_plot_components.png')
            plt.close(fig)

        except IOError:
            print("No Hyperion RT output found, SED not plotted.")

#convenience function to write dust parameter file 'BHDust.input' for Hyperion BHDust calculator (separate program)
    def write_bhmie_file(par):
        f=open('BHDust.input','w')
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