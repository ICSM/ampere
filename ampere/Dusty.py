# I have started with version 4.2 of dusty (https://github.com/ivezic/dusty)

# the model is called in the following way:
# dusty model.inp (the name is actually flexible, the extension is not)

# this generates a number of outputs in the current directory:
# model.NNN = runtime messages for each model iteration
# model.out = summery of the model inputs and its outputs (like A_V and dust temps)
# model.rNNN = radial profiles for each model iteration
# model.sNNN = output spectrum for each model iteration
# model.Xsec = effective extinction and absorption curves

# for the moment the most interesting output is model.sNNN where NNN is the largest available.
# but model.out might also contain interesting info
# perhaps in the future it would be interestingto be able to fit the radial profiles


# the input file looks like this:

#  ----------------------------------------------------------------------
#     Input data for DUSTY                                               
#  ---------------------------------------------------------------------- 
#  This is an input file for radiative transfer code DUSTY, version 4.0. 
#  NOTE: this input file is not compatible with old versions of Dusty 
#  due to the added new input options. Examples of each input option are 
#  given at the end of this file. For a more detailed description please 
#  refer to the Manual. 
#  
# 
#  The input file has a free format, text and empty lines can be entered
#  arbitrarily. All lines that start with the '*' sign are copied to the
#  output, and can be used to print out notes and comments. This option 
#  can also be useful when the program fails for some mysterious reason 
#  and you want to compare its output with an exact copy of the input line
#  as it was read in before processing by DUSTY. The occurrence of relevant 
#  numerical input, which is entered in standard FORTRAN conventions, is 
#  flagged by the equal sign `='. The only restrictions are that all required 
#  input entries must be specified, and in the correct order; the most likely 
#  source of an input error is failure to comply with these requirements. 
#  Recall, also, that FORTRAN requires a carriage return termination of the 
#  file's last line if it contains relevant input. Single entries are always 
#  preceded by the equal sign, `=', and must be padded by blanks on both sides; 
#  the terminating blank can be optionally preceded with a comma. For example: 
#  T = 10,000 K as well as Temperature = 1.E4 degrees and simply T = 10000.00 
#  are all equivalent, legal input entries (note that comma separations of long 
#  numbers are permitted).  Some input is entered as a list, in which case the 
#  first member is preceded by `=' and each subsequent member must be preceded 
#  by a blank (an optional comma can be entered before the blank for additional 
#  separation); for example, Temperatures  = 1E4, 2E4 30,000. Because of the 
#  special role of '=' as a flag for input entry, care must be taken not to 
#  introduce any '=' except when required.  All text following the  '%' sign 
#  is ignored (as in TeX) and this can be used to comment out material that 
#  includes '=' signs. For example, different options for the same physical 
#  property may require a different number of input entries. By commenting out 
#  with '%', all options may be retained in the input file with only the 
#  relevant one switched on.
# >
# 
# * ----------------------------------------------------------------------
# * NOTES:                                                                
# * Sample input file (external.inp)
# * Spherical dust distribution with
# * constant density profile
# * heated by external radiation with a temperature of 5000K
# * for composite dust grain 70% silicates 30% carbon
# * for 2 optical depth (10,100)
# * ----------------------------------------------------------------------
# 
# I. GEOMETRY %(available options: sphere\slab)
# 
#      geometry = sphere
# 
# II. PHYSICAL PARAMETERS                                                 
#      1) Central Source %(available options: on\off)
# 
#         	central = off
# 
#      2) External source  %(available options: on\off)
#                
# 		external = on
# 
#      2.1) Shape: %(available options: black_body\engelkd_marengo\power_law\file_lambda_f_lambda\file_f_lambda\file_f_nu)
# 
#       	        Spectral shape = black_body  
#                 Number of BB = 1
#                 Temperature = 5000 K 
# 
#      2.2) Scale: %(available options: flux\Lum_r1\energy_den\dilutn_fac\T1)
#         
# 		Scale = energy_den   % Td at the inner boundary
# 		flux = 1e4
# 
#      3) Dust Properties 
#      
#      3.1 Chemical composition %(available options: common_grain_composite\common_and_addl_grain\tabulated)
# 
#                 optical properties index = common_grain_composite
#      		Abundances for supported grain types:
#                	Sil-Ow  Sil-Oc  Sil-DL  grf-DL  amC-Hn  SiC-Pg 
#            x =  0.00    0.70    0.00    0.00    0.30    0.00
# 	        SIZE DISTRIBUTION = MRN
# 	        Tsub = 1500.
# 	   
#      4) Density Distribution %(available options: powd\expd\rdw\rdwa\usr_suppld)
# 
#          	 density type = POWD
#         	 number of powers = 1                
#         	 shell's relative thickness = 100.
#         	 power = 0.
# 
# 
#      5) Optical Depth: %(available options: linear\logarithmic\user_supplied)
#    
# 		 grid type = logarithmic % log grid
#         	 lambda0 = 0.55 micron   % fiducial wavelength
# 		 % minimum optical depth @ fiducial wavelength
#         	 tau(min) = 10.0; 
# 		 % maximum optical depth @ fiducial wavelength
# 		 tau(max) = 100.0  
#         	 number of models = 2
# 
#   ----------------------------------------------------------------------
#                                                                         
#   III. NUMERICS                                                           
#       
#      - accuracy for flux conservation = 0.10
#  
#   ----------------------------------------------------------------------
#                                                                         
#   IV. OUTPUT PARAMETERS                                                 
# 
#   	The flags governing file production are as follows: 
#   	If flag.eq.0 the particular file(s) is not produced. If flag.eq.1
# 	all model results are in corresponding files with extensions 'spp' 
# 	(sp.properties), 'stb' (spectra), 'itb' (images and visibilities, 
# 	if chosen), 'rtb' (radial profiles) and 'mtb' (messages).  If 
# 	flag.eq.2 each model result is in a separate corresponding file, 
# 	with visibilities contained in fname.i##. If the images flag.eq.3 
# 	the visibilities will be in separate files fname.v## (the flag for 
# 	visibilities has to be the same as for images).
# 	Note that choosing imaging output requires additional input data 
# 	(please refer to the exmaples below or to the Manual).
# 
# 
#         FILE DESCRIPTION                               FLAG        
#        ------------------------------------------------------------     
#        - detailed spectra for each model;           fname.s### = 2
#        - images at specified wavelengths;           fname.i### = 0
#        - en.density at specified radii;             fname.j### = 0
#        - radial profiles for each model;            fname.r### = 2
#        - detailed run-time messages;                fname.m### = 2
#        ------------------------------------------------------------- 
#  
# 
#   The end of the input parameters listing
#   **********************************************************************
#  
 

# in order to create a valid input file we will read a template input-file (dusty_inp.tmpl)
# with $ values for those parameters we want to control.
# use the following code to substitute them:

# from string import Template
# # dictionary of variables with values to be replaced
# var_dict=dict( 'geometry'=this_geom, 'shape'=this_shape, ... )
# dusty_inp = Template( open( 'dusty_inp.tmpl','r' ).read() )
# #do the substitution
# dusty_inp.substitute(var_dict)
# dusty_inp_file = open('model.inp','w')
# dusty_inp_file.write(dusty_inp)

# We have to think a bit how to do this because the sub-options sometimes change depending on the top-level option.
# for example a powerlaw density profile has different keywords than an exponential dependency
# perhaps it is better to construct bit by bit the input file rather than with a global template







# Input: wavelengths from obs., dust model (q_ij), parameters A, B, C_j
# j = 0, nmaterial-1
# Wavelengths need include both pivotal wavelengths for photometry and wavelengths for (multiple) spectra
# Regrid the q values onto the wavelengths
# Execution: calculate the model flux; return the result

import numpy as np
from astropy import constants as const
from astropy import units as u
from astropy.modeling import blackbody
from .models import AnalyticalModel

class DustySpectrum(AnalyticalModel):
    def __init__(self, wavelengths, flatprior=True,
                 normWave = 1., sigmaNormWave = 1.,opacityFileList=opacities,
                 redshift = False, lims=np.array([[0,1e6],[-100,100],[-10,10],[0,np.inf]]),
                 **kwargs):
        self.wavelength = wavelengths #grid of observed wavelengths to calculate BB for
        self.flatprior = flatprior #whether to assume flat priors
        self.normWave = normWave #wavelength at which opacity is normalised
        self.sigmaNormWave = sigmaNormWave #value to which the opacity is normalised at wavelength normWave
        self.lims = lims #in same order as inputs to __call__
        print(self.lims)
        #Define opacities for use in model calculation
        import os
        from scipy import interpolate
        opacityDirectory = os.path.dirname(__file__)+'/Opacities/'
        opacityFileList = opacities
        #opacityFileList = np.array(opacityFileList)[['.q' in zio for zio in opacityFileList]] # Only files ending in .q are valid (for now)
        nSpecies = opacities.__len__()
        opacity_array = np.zeros((wavelengths.__len__(), nSpecies))
        for j in range(nSpecies):
            tempData = np.loadtxt(opacityDirectory + opacityFileList[j], comments = '#')
            print(opacityFileList[j])
            tempWl = tempData[:, 0]
            tempOpac = tempData[:, 1]
            f = interpolate.interp1d(tempWl, tempOpac, assume_sorted = False)
            opacity_array[:,j] = f(self.restwaves)#wavelengths)
        self.opacity_array = opacity_array
        self.nSpecies = nSpecies
        self.redshift = redshift
        if redshift:
            from astropy.cosmology import FlatLambdaCDM
            self.cosmo=FlatLambdaCDM(H0=70, Om0=0.3)

    def __call__(self, t = 1., scale = 1., index = 1., dist=1.,
                 weights=np.array([0.]),
                 **kwargs):
        if self.redshift:
            z = dist
            dist = cosmo.luminosity_distance(z).to(u.m)
            freq = const.c.value / ((self.wavelength/(1.+z))*1e-6)
        else:
            dist = dist*u.pc.to(u.m)
            freq = const.c.value / (self.wavelength*1e-6)
        #Simple bug catches for inputs to opacity spectrum
        if len(weights) != self.nSpecies :
            print('Number of weights must be same as number of species')
        #if scale <= 0.0:
        #    print('Scale factor must be positive and non-zero.') #not relevant - scale is a log
        if t <= 0.0:
            print('Temperature must be positive and in Kelvin.')


    def dusty_write_input_file(self):
        # we constuct a dusty input file piece by piece
        dusty_inp_file = open("dusty_model.inp","w")

        # geometry
        if self.geometry.lower() == "slab":
            dusty_inp_file.write("geometry = SLAB")

            # slab needs also angular distribution = isotropic/directional
            if self.source_angular_distribution.lower() == "isotropic":
                dusty_inp_file.write("anugular distribution = ISOTROPIC")
            elif self.source_angular_distribution.lower() == "directional":
                dusty_inp_file.write("anugular distribution = DIRECTIONAL")
                dusty_inp_file.write("illumination angle = "+str(self.source_illumination_angle)+" degrees")
            else:
                raise ValueError('DustySpectrum: no valid source angular distribution is specfied.')
                
        elif self.geometry.lower() == "sphere":
            dusty_inp_file.write("geometry = sphere")
            if self.source_central.lower() == "on" and self.source_external.lower() == "on":
                raise ValueError('DustySpectrum: the spherical geometry does not allow external and central irradiation.')
            if self.source_central.lower() == "off":
                dusty_inp_file.write("central = off")
                dusty_inp_file.write("external = on")
            else:
                dusty_inp_file.write("central = on")
                dusty_inp_file.write("external = off")

        elif  self.geometry.lower() == "sphere_matrix":
            dusty_inp_file.write("geometry = sphere_matrix")
            if self.source_central.lower() == "on" and self.source_external.lower() == "on":
                raise ValueError('DustySpectrum: the spherical geometry does not allow external and central irradiation.')
            if self.source_central.lower() == "off":
                dusty_inp_file.write("central = off")
                dusty_inp_file.write("external = on")
            else:
                dusty_inp_file.write("central = on")
                dusty_inp_file.write("external = off")
        else:
            raise ValueError('DustySpectrum: no valid geometry specfied.')
            
        # spectral shape
        if self.spectralshape.lower() == "blackbody":
            dusty_inp_file.write("spectral shape = BLACK_BODY")
            dusty_inp_file.write("Number of BB = "+str(len(self.blackbody_temperatures)))
            dusty_inp_file.write("Temperatures = "+', '.join(map(str, self.blackbody_temperatures))+' K')
            dusty_inp_file.write("Luminosities = "+', '.join(map(str, self.blackbody_luminosities)))
        elif self.spectralshape.lower() == "engelke":
            dusty_inp_file.write("spectral shape = ENGELKE_MARENGO")
            dusty_inp_file.write("Temperature = "+str(self.engelke_temperature)+' K')
            dusty_inp_file.write("SiO absorption depth = "+str(self.engelke_sio_depth)+' percents')
        elif self.spectralshape.lower() == "powerlaw":
            dusty_inp_file.write("spectral shape = POWERLAW")
            dusty_inp_file.write("N = "+str(len(self.powerlaw_lambda)-1))
            dusty_inp_file.write("lambda = "+', '.join(map(str, self.powerlaw_lambda))+' micron')
            dusty_inp_file.write("k = "+', '.join(map(str, self.powerlaw_k)))
        elif self.spectralshape.lower() == "file_lambda_f_lambda":
            dusty_inp_file.write("spectral shape = FILE_LAMBDA_F_LAMBDA")
            dusty_inp_file.write("filename = "+self.spectralshape_filename)
        elif self.spectralshape.lower() == "file_f_lambda":
            dusty_inp_file.write("spectral shape = FILE_F_LAMBDA")
            dusty_inp_file.write("filename = "+self.spectralshape_filename)
        elif self.spectralshape.lower() == "file_f_nu":
            dusty_inp_file.write("spectral shape = FILE_F_NU")
            dusty_inp_file.write("filename = "+self.spectralshape_filename)
        else:
            raise ValueError('DustySpectrum: no valid spectral shape specfied.')

        # scale of the input spectrum (FLUX/LUM_R1/ENERGY_DEN/DILUTN_FAC/T1)
        # note that the doc lists several combinations of shape,geometry and scaling that are either required or
        # mutually exclusive. Best to catch those in __init__
        if self.spectralscale.lower() == "flux":
            dusty_inp_file.write("Scale: type of entry = FLUX")
            dusty_inp_file.write("Fe = "+str(self.spectralscale_flux_entering)+" W/m^2")
        elif self.spectralscale.lower() == "luminosity":
            dusty_inp_file.write("Scale: type of entry = LUM_R1")
            dusty_inp_file.write("L = "+str(self.spectralscale_luminosity)+" % in L_sun")
            dusty_inp_file.write("d = "+str(self.spectralscale_distance)+" cm")
        elif self.spectralscale.lower() == "energy-density":
            dusty_inp_file.write("Scale: type of entry = ENERGY_DEN")
            dusty_inp_file.write("J = "+str(self.spectralscale_energy_density)+" W/m^2")
        elif self.spectralscale.lower() == "dilution-factor":
            dusty_inp_file.write("Scale: type of entry = DILUTN_FAC")
            dusty_inp_file.write("W = "+str(self.spectralscale_dilution_factor))
        elif self.spectralscale.lower() == "inner-temperature":
            dusty_inp_file.write("Scale: type of entry = T1")
            dusty_inp_file.write("Td = "+str(self.spectralscale_inner_temperature)+" K")
        else:
            raise ValueError('DustySpectrum: no valid spectral scaling specfied.')

# repeat a section if slab geometry and right is on
        # geometry
        if (self.geometry.lower() == "slab") and (self.geometry_slab_right.lower() == "on"):
            dusty_inp_file.write("right = ON")

            # slab needs also angular distribution = isotropic/directional
            if self.source_right_angular_distribution.lower() == "isotropic":
                dusty_inp_file.write("anugular distribution = ISOTROPIC")
            elif self.source_right_angular_distribution.lower() == "directional":
                dusty_inp_file.write("anugular distribution = DIRECTIONAL")
                dusty_inp_file.write("illumination angle = "+str(self.source_right_illumination_angle)+" degrees")
            else:
                raise ValueError('DustySpectrum: no valid right source angular distribution is specfied.')
            
            # spectral shape
            if self.right_spectralshape.lower() == "blackbody":
                dusty_inp_file.write("spectral shape = BLACK_BODY")
                dusty_inp_file.write("Number of BB = "+str(len(self.right_blackbody_temperatures)))
                dusty_inp_file.write("Temperatures = "+', '.join(map(str, self.right_blackbody_temperatures))+' K')
                dusty_inp_file.write("Luminosities = "+', '.join(map(str, self.right_blackbody_luminosities)))
            elif self.right_spectralshape.lower() == "engelke":
                dusty_inp_file.write("spectral shape = ENGELKE_MARENGO")
                dusty_inp_file.write("Temperature = "+str(self.right_engelke_temperature)+' K')
                dusty_inp_file.write("SiO absorption depth = "+str(self.right_engelke_sio_depth)+' percents')
            elif self.right_spectralshape.lower() == "powerlaw":
                dusty_inp_file.write("spectral shape = POWERLAW")
                dusty_inp_file.write("N = "+str(len(self.right_powerlaw_lambda)-1))
                dusty_inp_file.write("lambda = "+', '.join(map(str, self.right_powerlaw_lambda))+' micron')
                dusty_inp_file.write("k = "+', '.join(map(str, self.right_powerlaw_k)))
            elif self.right_spectralshape.lower() == "file_lambda_f_lambda":
                dusty_inp_file.write("spectral shape = FILE_LAMBDA_F_LAMBDA")
                dusty_inp_file.write("filename = "+self.right_spectralshape_filename)
            elif self.right_spectralshape.lower() == "file_f_lambda":
                dusty_inp_file.write("spectral shape = FILE_F_LAMBDA")
                dusty_inp_file.write("filename = "+self.right_spectralshape_filename)
            elif self.right_spectralshape.lower() == "file_f_nu":
                dusty_inp_file.write("spectral shape = FILE_F_NU")
                dusty_inp_file.write("filename = "+self.right_spectralshape_filename)
            else:
                raise ValueError('DustySpectrum: no valid right spectral shape specfied.')

            # scale of the input spectrum (FLUX/LUM_R1/ENERGY_DEN/DILUTN_FAC/T1)
            # note that the doc lists several combinations of shape,geometry and scaling that are either required or
            # mutually exclusive. Best to catch those in __init__
            if self.right_spectralscale.lower() == "flux":
                dusty_inp_file.write("Scale: type of entry = FLUX")
                dusty_inp_file.write("Fe = "+str(self.right_spectralscale_flux_entering)+" W/m^2")
            elif self.right_spectralscale.lower() == "luminosity":
                dusty_inp_file.write("Scale: type of entry = LUM_R1")
                dusty_inp_file.write("L = "+str(self.right_spectralscale_luminosity)+" % in L_sun")
                dusty_inp_file.write("d = "+str(self.right_spectralscale_distance)+" cm")
            elif self.right_spectralscale.lower() == "energy-density":
                dusty_inp_file.write("Scale: type of entry = ENERGY_DEN")
                dusty_inp_file.write("J = "+str(self.right_spectralscale_energy_density)+" W/m^2")
            elif self.right_spectralscale.lower() == "dilution-factor":
                dusty_inp_file.write("Scale: type of entry = DILUTN_FAC")
                dusty_inp_file.write("W = "+str(self.right_spectralscale_dilution_factor))
            elif self.right_spectralscale.lower() == "inner-temperature":
                dusty_inp_file.write("Scale: type of entry = T1")
                dusty_inp_file.write("Td = "+str(self.right_spectralscale_inner_temperature)+" K")
            else:
                raise ValueError('DustySpectrum: no valid right spectral scaling specfied.')



        
#      3.1 Chemical composition %(available options: common_grain_composite\common_and_addl_grain\tabulated)

        if self.grain_composition.lower() == "common_grain_composite":
            dusty_inp_file.write("optical properties index = COMMON_GRAIN_COMPOSITE")
            dusty_inp_file.write("Abundances for supported grain types:")
            dusty_inp_file.write("Sil-Ow  Sil-Oc  Sil-DL  grf-DL  amC-Hn  SiC-Pg")
            dusty_inp_file.write("x = "+map(str, self.grain_fractional_abundances))
        elif self.grain_composition.lower() == "common_and_addl_grain":
            dusty_inp_file.write("optical properties index = COMMON_AND_ADDL_GRAIN")
            dusty_inp_file.write("Abundances for supported grain types:")
            dusty_inp_file.write("Sil-Ow  Sil-Oc  Sil-DL  grf-DL  amC-Hn  SiC-Pg")
            dusty_inp_file.write("x = "+map(str, self.grain_fractional_abundances[0:6]))
            dusty_inp_file.write("Number of additional components = 3, propperties listed in files")
            dusty_inp_file.write("amC-zb1.nk")
            dusty_inp_file.write("amC-zb2.nk")
            dusty_inp_file.write("amC-zb3.nk")
            dusty_inp_file.write("Abundances for these components = "+map(str, self.grain_fractional_abundances[6:9]))
        elif self.grain_composition.lower() == "tabulated":
            dusty_inp_file.write("optical properties index = TABULATED")
            dusty_inp_file.write("X-sections input file = "+self.grain_optical_properties_filename)
        else:
            raise ValueError('DustySpectrum: no valid grain composition specfied.')

        if self.grainsize_distribution.lower() == "mrn":
            dusty_inp_file.write("Size distribution = MRN")
        elif self.grainsize_distribution.lower() == "modified mrn":
            dusty_inp_file.write("Size distribution = MODIFIED_MRN")
            dusty_inp_file.write("q = "+str(self.grainsize_q))
            dusty_inp_file.write("a(min) = "+str(self.grainsize_amin)+" micron")
            dusty_inp_file.write("a(max) = "+str(self.grainsize_amax)+" micron")
        elif self.grain_size_distribution.lower() == "kmh":
            dusty_inp_file.write("Size distribution = KMH")
            dusty_inp_file.write("q = "+str(self.grainsize_q))
            dusty_inp_file.write("a(min) = "+str(self.grainsize_amin)+" micron")
            dusty_inp_file.write("a0 = "+str(self.grainsize_a0)+" micron")
        else:
            raise ValueError('DustySpectrum: no valid grain size distribution specfied.')

# 	   
#      4) Density Distribution %(available options: powd\expd\rdw\rdwa\usr_suppld)
# 
#          	 density type = POWD
#         	 number of powers = 1                
#         	 shell's relative thickness = 100.
#         	 power = 0.
# 
# 
#      5) Optical Depth: %(available options: linear\logarithmic\user_supplied)
#    
# 		 grid type = logarithmic % log grid
#         	 lambda0 = 0.55 micron   % fiducial wavelength
# 		 % minimum optical depth @ fiducial wavelength
#         	 tau(min) = 10.0; 
# 		 % maximum optical depth @ fiducial wavelength
# 		 tau(max) = 100.0  
#         	 number of models = 2




        bb = blackbody.blackbody_nu(freq,t).to(u.Jy / u.sr).value
        bb = bb / dist**2
        bb = bb * 10**(scale) * self.sigmaNormWave * ((self.wavelength / self.normWave)**index)
        #Subtract the sum of opacities from the blackbody continuum to calculate model spectrum
        fModel = bb * (1.0 - (np.matmul(self.opacity_array, weights)))
        self.modelFlux = fModel
        #return (blackbody.blackbody_nu(const.c.value*1e6/self.wavelengths,t).to(u.Jy / u.sr).value / (dist_lum.value)**2 * kappa230.value * ((wave/230.)**betaf) * massf) #*M_sun.cgs.value

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



