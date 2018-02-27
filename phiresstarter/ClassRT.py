
# Input: wavelengths from obs., dust model (q_ij), parameters A, B, C_j
# j = 0, nmaterial-1
# Wavelengths need include both pivotal wavelengths for photometry and wavelengths for (multiple) spectra
# Regrid the q values onto the wavelengths
# Execution: calculate the model flux; return the result

import numpy as np

class RadTransWrap(object):
    '''Input: fit parameters (multiplicationFactor, powerLawIndex, relativeAbundances), 
              opacities (opacity_matrix): q_ij where i = wavelength, j = species
              wavelength grid from data (wavelengths)
    Output: model fluxes (modelFlux)'''

	def __init__(self, **kwargs):
            import os#, os.path
            opacityDirectory = './Opacities/'
            opacityFileList = os.listdir(opacityDirectory)
            opacityFileList.sort()
            nSpecies = opacityFileList.__len__() # Warning: need to check all files found are ordinary files (not directory or such). Look up os.path.isfile()
            opacityData = np.loadtxt(opacityDirectory+'ss_Dorschneretal1995_Olivine_0.10.q', comments='#')
            # To finish: Construct opacity_matrix
           pass
        
        def __call__(self, multiplicationFactor = 1, powerLawIndex = 2, relativeAbundances = np.ones(nSpecies)/nSpecies, **kwargs):
            fModel = (np.matmul(opacity_matrix, relativeAbundances)+1)*wavelengths**powerLawIndex*multiplicationFactor
            self.modelFlux = fModel
