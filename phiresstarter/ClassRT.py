
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
            pass
        
        def __call__(self, multiplicationFactor = 1, powerLawIndex = 2, relativeAbundances = np.ones(nSpecies)/nSpecies, **kwargs):
            fModel = (np.matmul(opacity_matrix, relativeAbundances)+1)*wavelengths**powerLawIndex*multiplicationFactor
            self.modelFlux = fModel
