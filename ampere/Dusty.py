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

# is it better to use all these explicit keyword arguments or use the kwargs mechanism in python?
class DustySpectrum(AnalyticalModel):
    def __init__(self, wavelengths, flatprior=True,
                 normWave = 1., sigmaNormWave = 1.,opacityFileList=opacities,
                 redshift = False, lims=np.array([[0,1e6],[-100,100],[-10,10],[0,np.inf]])
                 # these are all the possible input values. Some have defaults most are not set
                 # the defaults are set to what I think might be a reasonable basic setup
                 #
                 # The checks below will complain about needed options that are set to None
                 # If these are going to be variables set them to "variable"
                 #
                 ,geometry="sphere"
                 ,geometry_illumination_angle=None
                 ,geometry_angular_distribution=None
                 ,geometry_central="on"
                 ,geometry_external ="off"
                 ,spectralshape="blackbody"
                 ,spectralshape_blackbody_temperatures=None
                 ,spectralshape_blackbody_luminosities=None
                 ,spectralshape_powerlaw_lambda=None
                 ,spectralshape_powerlaw_k=None
                 ,spectralshape_engelke_temperature=None
                 ,spectralshape_engelke_sio_depth=None
                 ,spectralshape_filename=None
                 ,spectralscale="inner temperature"
                 ,spectralscale_flux_entering=None
                 ,spectralscale_luminosity=None
                 ,spectralscale_distance=None
                 ,spectralscale_dilution_factor=None
                 ,spectralscale_energy_density=None
                 ,spectralscale_inner_temperature=None
                 ,geometry_toggle_right=None
                 ,geometry_right_angular_distribution=None
                 ,geometry_right_illumination_angle=None
                 ,right_spectralshape=None
                 ,right_spectralshape_blackbody_luminosities=None
                 ,right_spectralshape_blackbody_temperatures=None
                 ,right_spectralshape_engelke_sio_depth=None
                 ,right_spectralshape_engelke_temperature=None
                 ,right_spectralshape_powerlaw_k=None
                 ,right_spectralshape_powerlaw_lambda=None
                 ,right_spectralshape_filename=None
                 ,right_spectralscale=None
                 ,right_spectralscale_dilution_factor=None
                 ,right_spectralscale_distance=None
                 ,right_spectralscale_energy_density=None
                 ,right_spectralscale_flux_entering=None
                 ,right_spectralscale_inner_temperature=None
                 ,right_spectralscale_luminosity=None
                 ,density_distribution=None
                 ,density_transition_radii=None
                 ,density_powers=None
                 ,density_outer_radius=None
                 ,density_falloff_rate=None
                 ,density_filename=None
                 ,grain_composition="common_grain_composite"
                 ,grain_fractional_abundances=None
                 ,grain_optical_properties_filename=None
                 ,grain_sublimation_temperature=1500.
                 ,grainsize_distribution="mrn"
                 ,grainsize_amax=None
                 ,grainsize_amin=None
                 ,grainsize_q=None
                 ,grainsize_a0=None
                 ,tau_grid="linear"
                 ,tau_filename=None
                 ,tau_max=None
                 ,tau_min=None
                 ,tau_nmodels=1
                 ,tau_wavelength=0.55
                 ,flux_conservation_accuracy=0.10
                 ,**kwargs):

        # make sure to lowercase all model input parameters such that we do not have to use .lower() when checking later
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

        # check the validity of the inputs.
        # make all sure that all self.xxx string values are lowercase

        # GEOMETRY:
        # slab needs angular distribution
        ## Note that for the slab geometry the density distribution is irrelevant
        ## Should we warn if densities are specief for slab?
        
        # directional distribution needs illumination angle < 85

        # toggle_right is only valid in the slab case

        # sphere can have only internal or external illumination

        # sphere_matrix can have only central illumination

        # SPECTRAL SHAPES:
        # black body need temperature(s)
        # if more than 1 temperature then relative luminosities are needed

        # engelke_marengo needs temperature and SiO depth

        # GEOMETRY BLOCK
        if geometry.lower() == "sphere":
            self.geometry="sphere"
            spherical=True
        elif geometry.lower() == "sphere_matrix":
            self.geometry="sphere_matrix"
            spherical=True
        elif geometry.lower() == "slab":
            self.geometry="slab"
            spherical=False
        else:
            raise ValueError('Dusty.__init__: geometry must be set (slab/sphere/sphere_matrix)')

        if spherical:
            if geometry_central.lower()=="on":
                self.geometry_central=True
            else:
                self.geometry_central=False
            if geometry_external.lower()=="on":
                self.geometry_external=True
            else:
                self.geometry_external=False
            if not self.geometry_external and not self.geometry_central:
                raise ValueError('Dusty.__init__: spherical needs either geometry_central or geometry_external set to on')
            if self.geometry_external and self.geometry_central:
                raise ValueError('Dusty.__init__: spherical cannot handle central AND external radiation')
        else:
            if geometry_angular_distribution.lower() == "isotropic":
                self.geometry_angular_distribution="isotropic"
            elif geometry_angular_distribution.lower() == "directional":
                self.geometry_angular_distribution="directional"
                if geometry_illumination_angle == None:
                    raise ValueError('Dusty.__init__: slab with directional needs an illumination angle')
                elif str(geometry_illumination_angle).lower() == "variable":
                    #
                elif geometry_illumination_angle > 85.0:
                    raise ValueError('Dusty.__init__: slab directional: illumination angle needs to be smaller than 85')
                else:
                    self.geometry_illumination_angle = geometry_illumination_angle
            else:
                raise ValueError('Dusty.__init__: slab geometry:specify angular distribution of the radiation (isotropic/directional)')
                
            if geometry_toggle_right.lower() == "on":
                self.geometry_toggle_right=True
                if geometry_right_angular_distribution.lower() == "isotropic":
                    self.geometry_right_angular_distribution="isotropic"
                elif geometry_right_angular_distribution.lower() == "directional":
                    self.geometry_right_angular_distribution="directional"
                    if geometry_right_illumination_angle == None:
                        raise ValueError('Dusty.__init__: slab-right with directional needs an illumination angle')
                    elif str(geometry_right_illumination_angle).lower() == "variable":
                        #
                    elif geometry_right_illumination_angle > 85.0:
                        raise ValueError('Dusty.__init__: slab-right directional: illumination angle needs to be smaller than 85')
                    else:
                        self.geometry_right_illumination_angle = geometry_right_illumination_angle
            else:
                raise ValueError('Dusty.__init__: slab-right geometry:specify angular distribution of the radiation (isotropic/directional)')
            elif geometry_toggle_right.lower() == None:
                #
            else:
                raise ValueError('Dusty.__init__: value passed to geometry_toggle_right not understood')
                    
            
        # SPECTRALSHAPE BLOCK
        needfile=False
        if spectralshape.lower() == "blackbody":
            self.spectralshape="blackbody"
            if spectralshape_blackbody_temperatures=None:
                raise ValueError('Dusty.__init__: spectralshape blackbody needs temperature(s) in K')
            elif str(spectralshape_blackbody_temperatures).lower()="variable":
                #
            else:
                self.spectralshape_blackbody_temperatures = spectralshape_blackbody_temperatures
                if len(spectralshape_blackbody_temperatures) == 1:
                    spectralshape_blackbody_luminosities = 1.
                else:
                    if spectralshape_blackbody_luminosities=None:
                        raise ValueError('Dusty.__init__: spectralshape multiple blackbody needs luminosities(s) in K')
                    elif str(spectralshape_blackbody_luminosities).lower()="variable":
                        #
                    else:
                        # in principle we should check that temperatures and luminosities are equal length
                        self.spectralshape_blackbody_luminosities = spectralshape_blackbody_luminosities
                
        elif spectralshape.lower() == "engelke":
            self.spectralshape="engelke"
            if spectralshape_engleke_temperatures=None:
                raise ValueError('Dusty.__init__: spectralshape engleke needs a temperature in K')
            elif str(spectralshape_engelke_temperature).lower()="variable":
                #
            else:
                self.spectralshape_engelke_temperature = spectralshape_engelke_temperature
            if spectralshape_engelke_sio-depth=None:
                raise ValueError('Dusty.__init__: spectralshape engelke needs an SiO depth')
            elif str(spectralshape_engelke_sio_depth).lower()="variable":
                #
            else:
                self.spectralshape_engelke_sio_depth = spectralshape_engelke_sio_depth

        elif spectralshape.lower() == "powerlaw":
            self.spectralshape="powerlaw"
            if spectralshape_powerlaw_lambda=None:
                raise ValueError('Dusty.__init__: spectralshape powerlaw needs wavelengths in microns')
            elif str(spectralshape_powerlaw_lambda).lower()="variable":
                #
            else:
                self.spectralshape_powerlaw_lambda = spectralshape_powerlaw_lambda
            if spectralshape_powerlaw_k=None:
                raise ValueError('Dusty.__init__: spectralshape multiple powerlaw needs k (spectral indices)')
            elif str(spectralshape_powerlaw_k).lower()="variable":
                #
            else:
                # in principle we should check that len(lambda) = len(k)+1
                self.spectralshape_powerlaw_k = spectralshape_powerlaw_k
                
        elif spectralshape.lower() == "file_lambda_f_lambda":
            self.spectralshape="file_lambda_f_lambda"
            needfile=True
        elif spectralshape.lower() == "file_f_lambda":
            self.spectralshape="file_f_lambda"
            needfile=True
        elif spectralshape.lower() == "file_f_nu":
            self.spectralshape="file_f_nu"
            needfile=True
        else:
            raise ValueError('Dusty.__init__: spectralshape must be set (blackbody/engelke/powerlaw/file_lambda_f_lambda/file_f_lambda/file_f_nu)')

        if needfile:
            if spectralshape_filename=None:
                raise ValueError('Dusty.__init__: spectralshape file_* needs a filename')
            else:
                # in principle we should check that the file exists
                self.spectralshape_filename = spectralshape_filename


        # RIGHT_SPECTRALSHAPE BLOCK
        if self.geometry_toggle_right:
            needfile=False
            if right_spectralshape.lower() == "blackbody":
                self.right_spectralshape="blackbody"
                if right_spectralshape_blackbody_temperatures=None:
                    raise ValueError('Dusty.__init__: right_spectralshape blackbody needs temperature(s) in K')
                elif str(right_spectralshape_blackbody_temperatures).lower()="variable":
                    #
                else:
                    self.right_spectralshape_blackbody_temperatures = right_spectralshape_blackbody_temperatures
                    if len(right_spectralshape_blackbody_temperatures) == 1:
                        right_spectralshape_blackbody_luminosities = 1.
                    else:
                        if right_spectralshape_blackbody_luminosities=None:
                            raise ValueError('Dusty.__init__: right_spectralshape multiple blackbody needs luminosities(s) in K')
                        elif str(right_spectralshape_blackbody_luminosities).lower()="variable":
                            #
                        else:
                            # in principle we should check that temperatures and luminosities are equal length
                            self.right_spectralshape_blackbody_luminosities = right_spectralshape_blackbody_luminosities
                                    
            elif right_spectralshape.lower() == "engelke":
                self.right_spectralshape="engelke"
                if right_spectralshape_engleke_temperatures=None:
                    raise ValueError('Dusty.__init__: right_spectralshape engleke needs a temperature in K')
                elif str(right_spectralshape_engelke_temperature).lower()="variable":
                    #
                else:
                    self.right_spectralshape_engelke_temperature = right_spectralshape_engelke_temperature
                if right_spectralshape_engelke_sio-depth=None:
                    raise ValueError('Dusty.__init__: right_spectralshape engelke needs an SiO depth')
                elif str(right_spectralshape_engelke_sio_depth).lower()="variable":
                    #
                else:
                    self.right_spectralshape_engelke_sio_depth = right_spectralshape_engelke_sio_depth
                    
            elif right_spectralshape.lower() == "powerlaw":
                self.right_spectralshape="powerlaw"
                if right_spectralshape_powerlaw_lambda=None:
                    raise ValueError('Dusty.__init__: right_spectralshape powerlaw needs wavelengths in microns')
                elif str(right_spectralshape_powerlaw_lambda).lower()="variable":
                    #
                else:
                    self.right_spectralshape_powerlaw_lambda = right_spectralshape_powerlaw_lambda
                if right_spectralshape_powerlaw_k=None:
                    raise ValueError('Dusty.__init__: right_spectralshape multiple powerlaw needs k (spectral indices)')
                elif str(right_spectralshape_powerlaw_k).lower()="variable":
                    #
                else:
                    # in principle we should check that len(lambda) = len(k)+1
                    self.right_spectralshape_powerlaw_k = right_spectralshape_powerlaw_k
                
            elif right_spectralshape.lower() == "file_lambda_f_lambda":
                self.right_spectralshape="file_lambda_f_lambda"
                needfile=True
            elif right_spectralshape.lower() == "file_f_lambda":
                self.right_spectralshape="file_f_lambda"
                needfile=True
            elif right_spectralshape.lower() == "file_f_nu":
                self.right_spectralshape="file_f_nu"
                needfile=True
            else:
                raise ValueError('Dusty.__init__: right_spectralshape must be set (blackbody/engelke/powerlaw/file_lambda_f_lambda/file_f_lambda/file_f_nu)')

            if needfile:
                if right_spectralshape_filename=None:
                    raise ValueError('Dusty.__init__: right_spectralshape file_* needs a filename')
                else:
                    # in principle we should check that the file exists
                    self.right_spectralshape_filename = right_spectralshape_filename


        ,spectralscale = "inner temperature"
        ,spectralscale_flux_entering
        ,spectralscale_luminosity = 1e4
        ,spectralscale_distance = 1.4959787e+13 # 1 AU in cm
        ,spectralscale_dilution_factor
        ,spectralscale_energy_density
        ,spectralscale_inner_temperature

        ,right_spectralscale
        ,right_spectralscale_dilution_factor
        ,right_spectralscale_distance
        ,right_spectralscale_energy_density
        ,right_spectralscale_flux_entering
        ,right_spectralscale_inner_temperature
        ,right_spectralscale_luminosity

        ,density_distribution
        ,density_transition_radii
        ,density_powers
        ,density_outer_radius = 100
        ,density_falloff_rate
        ,density_filename

        ,grain_composition = "common_grain_composite"
        ,grain_fractional_abundances = [0,0,1,0,0,0]
        ,grain_optical_properties_filename
        ,grain_sublimation_temperature = 1500.

        ,grainsize_distribution = "mrn"
        ,grainsize_amax
        ,grainsize_amin
        ,grainsize_q
        ,grainsize_a0

        ,tau_grid
        ,tau_filename
        ,tau_max
        ,tau_min
        ,tau_nmodels = 1
        ,tau_wavelength = 0.55

        ,flux_conservation_accuracy = 0.10
        
    def __call__(self, *arg,
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

    def dusty_write_input_file_radiation_field(self,
                                               spectralshape=spectralshape,
                                               blackbody_temperatures=blackbody_temperatures,
                                               blackbody_luminosities=blackbody_luminosities,
                                               engelke_sio_depth=engelke_sio_depth,
                                               engelke_temperature=engelke_temperature,
                                               powerlaw_k=powerlaw_k,
                                               powerlaw_lambda=powerlaw_lambda,
                                               filename=filename,
                                               spectralscale=spectralscale,
                                               dilution_factor=dilution_factor,
                                               distance=distance,
                                               energy_density=energy_density,
                                               flux_entering=flux_entering,
                                               inner_temperature=inner_temperature,
                                               luminosity=luminosity)

        # spectral shape
        if spectralshape == "blackbody":
            dusty_inp_file.write("spectral shape = BLACK_BODY")
            dusty_inp_file.write("Number of BB = "+str(len(blackbody_temperatures)))
            dusty_inp_file.write("Temperatures = "+', '.join(map(str, blackbody_temperatures))+' K')
            dusty_inp_file.write("Luminosities = "+', '.join(map(str, blackbody_luminosities)))
        elif spectralshape == "engelke":
            dusty_inp_file.write("spectral shape = ENGELKE_MARENGO")
            dusty_inp_file.write("Temperature = "+str(engelke_temperature)+' K')
            dusty_inp_file.write("SiO absorption depth = "+str(engelke_sio_depth)+' percents')
        elif spectralshape == "powerlaw":
            dusty_inp_file.write("spectral shape = POWERLAW")
            dusty_inp_file.write("N = "+str(len(powerlaw_lambda)-1))
            dusty_inp_file.write("lambda = "+', '.join(map(str, powerlaw_lambda))+' micron')
            dusty_inp_file.write("k = "+', '.join(map(str, powerlaw_k)))
        elif spectralshape == "file_lambda_f_lambda":
            dusty_inp_file.write("spectral shape = FILE_LAMBDA_F_LAMBDA")
            dusty_inp_file.write("filename = "+filename)
        elif spectralshape == "file_f_lambda":
            dusty_inp_file.write("spectral shape = FILE_F_LAMBDA")
            dusty_inp_file.write("filename = "+filename)
        elif spectralshape == "file_f_nu":
            dusty_inp_file.write("spectral shape = FILE_F_NU")
            dusty_inp_file.write("filename = "+filename)
        else:
            raise ValueError('DustySpectrum: no valid spectral shape specfied.')

        # scale of the input spectrum (FLUX/LUM_R1/ENERGY_DEN/DILUTN_FAC/T1)
        if spectralscale == "flux":
            dusty_inp_file.write("Scale: type of entry = FLUX")
            dusty_inp_file.write("Fe = "+str(flux_entering)+" W/m^2")
        elif spectralscale == "luminosity":
            dusty_inp_file.write("Scale: type of entry = LUM_R1")
            dusty_inp_file.write("L = "+str(luminosity)+" % in L_sun")
            dusty_inp_file.write("d = "+str(distance)+" cm")
        elif spectralscale == "energy-density":
            dusty_inp_file.write("Scale: type of entry = ENERGY_DEN")
            dusty_inp_file.write("J = "+str(energy_density)+" W/m^2")
        elif spectralscale == "dilution-factor":
            dusty_inp_file.write("Scale: type of entry = DILUTN_FAC")
            dusty_inp_file.write("W = "+str(dilution_factor))
        elif spectralscale == "inner-temperature":
            dusty_inp_file.write("Scale: type of entry = T1")
            dusty_inp_file.write("Td = "+str(inner_temperature)+" K")
        else:
            raise ValueError('DustySpectrum: no valid spectral scaling specfied.')


    def dusty_write_input_file_radiation_field_left_right(self,side='left')
    if side == left:
        dusty_write_input_file_radiation_field(self,
                                               spectralshape=self.spectralshape,
                                               blackbody_temperatures=self.spectralshape_blackbody_temperatures,
                                               blackbody_luminosities=self.spectralshape_blackbody_luminosities,
                                               engelke_sio_depth=self.spectralshape_engelke_sio_depth,
                                               engelke_temperature=self.spectralshape_engelke_temperature,
                                               powerlaw_k=self.spectralshape_powerlaw_k,
                                               powerlaw_lambda=self.spectralshape_powerlaw_lambda,
                                               filename=self.spectralshape_filename,
                                               spectralscale=self.spectralscale,
                                               dilution_factor=self.spectralscale_dilution_factor,
                                               distance=self.spectralscale_distance,
                                               energy_density=self.spectralscale_energy_density,
                                               flux_entering=self.spectralscale_flux_entering,
                                               inner_temperature=self.spectralscale_inner_temperature,
                                               luminosity=self.spectralscale_luminosity)
    else:
        dusty_write_input_file_radiation_field(self,
                                               spectralshape=self.right_spectralshape,
                                               blackbody_temperatures=self.right_spectralshape_blackbody_temperatures,
                                               blackbody_luminosities=self.right_spectralshape_blackbody_luminosities,
                                               engelke_sio_depth=self.right_spectralshape_engelke_sio_depth,
                                               engelke_temperature=self.right_spectralshape_engelke_temperature,
                                               powerlaw_k=self.right_spectralshape_powerlaw_k,
                                               powerlaw_lambda=self.right_spectralshape_powerlaw_lambda,
                                               filename=self.right_spectralshape_filename,
                                               spectralscale=self.right_spectralscale,
                                               dilution_factor=self.right_spectralscale_dilution_factor,
                                               distance=self.right_spectralscale_distance,
                                               energy_density=self.right_spectralscale_energy_density,
                                               flux_entering=self.right_spectralscale_flux_entering,
                                               inner_temperature=self.right_spectralscale_inner_temperature,
                                               luminosity=self.right_spectralscale_luminosity)

    def dusty_write_input_file_geometry(self)
        # geometry
        if self.geometry == "slab":
            dusty_inp_file.write("geometry = SLAB")

            # slab needs also angular distribution = isotropic/directional
            if self.geometry_angular_distribution == "isotropic":
                dusty_inp_file.write("anugular distribution = ISOTROPIC")
            elif self.geometry_angular_distribution == "directional":
                dusty_inp_file.write("anugular distribution = DIRECTIONAL")
                dusty_inp_file.write("illumination angle = "+str(self.geometry_illumination_angle)+" degrees")
            else:
                raise ValueError('DustySpectrum: no valid source angular distribution is specfied.')

            dusty_write_input_file_radiation_field_left_right(self,'left')

            # repeat a section if slab geometry and right is on
            # geometry
            if (self.geometry == "slab") and (self.geometry_toggle_right == "on"):
                dusty_inp_file.write("right = ON")
                
                # slab needs also angular distribution = isotropic/directional
                if self.geometry_right_angular_distribution == "isotropic":
                    dusty_inp_file.write("anugular distribution = ISOTROPIC")
                elif self.geometry_right_angular_distribution == "directional":
                    dusty_inp_file.write("anugular distribution = DIRECTIONAL")
                    dusty_inp_file.write("illumination angle = "+str(self.geometry_right_illumination_angle)+" degrees")
                else:
                    raise ValueError('DustySpectrum: no valid right source angular distribution is specfied.')
                
                dusty_write_input_file_radiation_field_left_right(self,'right')
                
        elif self.geometry == "sphere" or self.geometry == "sphere_matrix":
            if self.geometry == "sphere":
                dusty_inp_file.write("geometry = SPHERE")
            elif self.geometry == "sphere_matrix":
                dusty_inp_file.write("geometry = SPHERE_MATRIX")

            if self.geometry_central == "on" and self.geometry_external == "on":
                raise ValueError('DustySpectrum: the spherical geometry does not allow external and central irradiation.')

            if self.geometry_central == "off":
                dusty_inp_file.write("central = OFF")
                dusty_inp_file.write("external = ON")
                dusty_write_input_file_radiation_field_left_right(self,'left')
            else:
                dusty_inp_file.write("central = ON")
                dusty_write_input_file_radiation_field_left_right(self,'left')
                dusty_inp_file.write("external = OFF")
        else:
            raise ValueError('DustySpectrum: no valid geometry specfied.')

        
    def dusty_write_input_file_grain_composition(self)
        # Chemical composition %(available options: common_grain_composite\common_and_addl_grain\tabulated)
        if self.grain_composition == "common_grain_composite":
            dusty_inp_file.write("optical properties index = COMMON_GRAIN_COMPOSITE")
            dusty_inp_file.write("Abundances for supported grain types:")
            dusty_inp_file.write("Sil-Ow  Sil-Oc  Sil-DL  grf-DL  amC-Hn  SiC-Pg")
            dusty_inp_file.write("x = "+map(str, self.grain_fractional_abundances))
        elif self.grain_composition == "common_and_addl_grain":
            dusty_inp_file.write("optical properties index = COMMON_AND_ADDL_GRAIN")
            dusty_inp_file.write("Abundances for supported grain types:")
            dusty_inp_file.write("Sil-Ow  Sil-Oc  Sil-DL  grf-DL  amC-Hn  SiC-Pg")
            dusty_inp_file.write("x = "+map(str, self.grain_fractional_abundances[0:6]))
            dusty_inp_file.write("Number of additional components = 3, propperties listed in files")
            dusty_inp_file.write("amC-zb1.nk")
            dusty_inp_file.write("amC-zb2.nk")
            dusty_inp_file.write("amC-zb3.nk")
            dusty_inp_file.write("Abundances for these components = "+map(str, self.grain_fractional_abundances[6:9]))
        elif self.grain_composition == "tabulated":
            dusty_inp_file.write("optical properties index = TABULATED")
            dusty_inp_file.write("X-sections input file = "+self.grain_optical_properties_filename)
        else:
            raise ValueError('DustySpectrum: no valid grain composition specfied.')

        dusty_inp_file.write("Sublimation Temperature = "+self.grain_sublimation_temperature+" K")
        
    def dusty_write_input_file_grainsize_distribution(self)
        if self.grainsize_distribution == "mrn":
            dusty_inp_file.write("Size distribution = MRN")
        elif self.grainsize_distribution == "modified mrn":
            dusty_inp_file.write("Size distribution = MODIFIED_MRN")
            dusty_inp_file.write("q = "+str(self.grainsize_q))
            dusty_inp_file.write("a(min) = "+str(self.grainsize_amin)+" micron")
            dusty_inp_file.write("a(max) = "+str(self.grainsize_amax)+" micron")
        elif self.grainsize_distribution == "kmh":
            dusty_inp_file.write("Size distribution = KMH")
            dusty_inp_file.write("q = "+str(self.grainsize_q))
            dusty_inp_file.write("a(min) = "+str(self.grainsize_amin)+" micron")
            dusty_inp_file.write("a0 = "+str(self.grainsize_a0)+" micron")
        else:
            raise ValueError('DustySpectrum: no valid grain size distribution specfied.')

    def dusty_write_input_file_density_distribution(self)
        # Density Distribution %(available options: powd\expd\rdw\rdwa\usr_suppld)
        if self.density_distribution == "powerlaw":
            dusty_inp_file.write("density type = POWD")
            dusty_inp_file.write("N = "+str(len(self.density_transition_radii)))
            dusty_inp_file.write("transition radii = "+', '.join(map(str, self.density_transition_radii)))
            dusty_inp_file.write("power indices = "+', '.join(map(str, self.density_powers)))
        elif self.density_distribution == "exponential":
            dusty_inp_file.write("density type = EXPD")
            dusty_inp_file.write("Y = "+str(self.density_outer_radius))
            dusty_inp_file.write("sigma = "+str(self.density_falloff_rate))
        elif self.density_distribution == "radiation driven wind":
            dusty_inp_file.write("density type = RDW")
            dusty_inp_file.write("Y = "+str(self.density_outer_radius))
        elif self.density_distribution == "radiation driven wind analytic":
            dusty_inp_file.write("density type = RDWA")
            dusty_inp_file.write("Y = "+str(self.density_outer_radius))
        elif self.density_distribution == "user supplied":
            dusty_inp_file.write("density type = USR_SUPPLD")
            dusty_inp_file.write("profile filename = "+self.density_filename)
        else:
            raise ValueError('DustySpectrum: no valid density distribution specfied.')

    def dusty_write_input_file_tau_grid(self)
        #      5) Optical Depth: %(available options: linear\logarithmic\user_supplied)
        # (SH Oct  4 2018)
        #      Some thought needs to go into this. Dusty can calculate many models that stop at different depths in the dust layer.
        #      For the sampler it is simplest to calculate only one model. But in terms of overheads it might be interestin to calculate several
        if self.tau_grid == "linear":
            dusty_inp_file.write("grid type = LINEAR")
            dusty_inp_file.write("lambda0 = "+str(self.tau_wavelength)+" micron")
            dusty_inp_file.write("tau(min) = "+str(self.tau_min))
            dusty_inp_file.write("tau(max) = "+str(self.tau_max))
            dusty_inp_file.write("number of models = "+str(self.tau_nmodels))
        elif self.tau_grid == "logarithmic":
            dusty_inp_file.write("grid type = LOGARITHMIC")
            dusty_inp_file.write("lambda0 = "+str(self.tau_wavelength)+" micron")
            dusty_inp_file.write("tau(min) = "+str(self.tau_min))
            dusty_inp_file.write("tau(max) = "+str(self.tau_max))
            dusty_inp_file.write("number of models = "+str(self.tau_nmodels))
        elif self.tau_grid == "user supplied":
            dusty_inp_file.write("grid type = USER_SUPPLIED")
            # not sure about the following
            dusty_inp_file.write("tau values filename = "+self.tau_filename)
        else:
            raise ValueError('DustySpectrum: no valid tau grid type specfied.')

    def dusty_write_input_file(self):
        # we constuct a dusty input file piece by piece
        dusty_inp_file = open("dusty_model.inp","w")
        dusty_write_input_file_geometry(self)
        dusty_write_input_file_grain_composition(self)
        dusty_write_input_file_grainsize_distribution(self)
        dusty_write_input_file_density_distribution(self)
        dusty_write_input_file_tau_grid(self)
        # NUMERICS                                                           
        dusty_inp_file.write("accuracy for flux conservation = "+self.flux_conservation_accuracy)

        #   IV. OUTPUT PARAMETERS                                                 
        #         FILE DESCRIPTION                               FLAG        
        #        ------------------------------------------------------------     
        #        - detailed spectra for each model;           fname.s### = 2
        #        - images at specified wavelengths;           fname.i### = 0
        #        - en.density at specified radii;             fname.j### = 0
        #        - radial profiles for each model;            fname.r### = 2
        #        - detailed run-time messages;                fname.m### = 2
        #        ------------------------------------------------------------- 
        dusty_inp_file.write("- detailed spectra for each model;           fname.s### = 2")
        dusty_inp_file.write("- images at specified wavelengths;           fname.i### = 0")
        dusty_inp_file.write("- en.density at specified radii;             fname.j### = 0")
        dusty_inp_file.write("- radial profiles for each model;            fname.r### = 0")
        dusty_inp_file.write("- detailed run-time messages;                fname.m### = 0")
        
        close(dusty_inp_file)
        

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




