from __future__ import print_function

import numpy as np



class BaseExtinctionLaws(object):
	"""
	"""

	def __init__(self, **kwargs):
		raise NotImplimentedError("")

	def __str__(self, **kwargs):
		raise NotImplimentedError("")

	def __repr__(self, **kwargs):
		raise NotImplimentedError("")

	def extinctionModel(self, **kwargs): 
		raise NotImplimentedError("")



class CCMExtinctionLaw(BaseExtinctionLaws):
	"""
	"""

        def __init__(self, **kwargs):
		raise NotImplimentedError("")

	def extinctionModel(self, extinctionModel, **kwargs):
                return CCMModel(self, wavelength, R_A, A_V, **kwargs)

	def a(self, wavelength, type_Wavelength, **kwargs):

		inverse_wavelength = 1 / wavelength

		if type_Wavelength == "IR":

			a = 0.574 * (inverse_wavelength**1.61)

		elif type_Wavelength == "NIR" or "Optical":

			y = inverse_wavelength - 1.82

			a = 1 + (0.17699*y) - (0.5044*(y**2)) - (0.02427*(y**3)) + (0.72085*(y**4)) + (0.01979*(y**5)) - (0.77530*(y**6)) + (0.32999*(y**7))

		elif type_Wavelength == "UV":

			if 8 >= inverse_wavelength >= 5.9: 

				P_a = -(0.04473((inverse_wavelength - 5.9)**2)) - (0.009779((inverse_wavelength - 5.9)**3))

			elif inverse_wavelength < 5.9:

				P_a = 0

			a = 1.752 - (0.316*inverse_wavelength) - (0.104 / ( ((inverse_wavelength - 4.67)**2) + 0.341) ) + P_a

		elif type_Wavelength == "FUV":

			a = -1.073 - (0.628*(inverse_wavelength - 8)) + (0.137 * ((inverse_wavelength - 8)**2)) - (0.070 * ((inverse_wavelength - 8)**3))

		else: 
			raise NotImplementedError("Incorrect type_Wavelength string") 


		return a



	def b(self, wavelength, type_Wavelength, **kwargs):

		inverse_wavelength = 1 / wavelength

		if type_Wavelength == "IR":

			b = 0.527 * (inverse_wavelength**1.61)

		elif type_Wavelength == "NIR" or "Optical":

			y = inverse_wavelength - 1.82

			b = (1.41338*y) - (2.28305*(y**2)) - (1.07233*(y**3)) - (5.38434*(y**4)) - (0.62251*(y**5)) + (5.30260*(y**6)) - (2.09002*(y**7))

		elif type_Wavelength == "UV":

			if 8 >= inverse_wavelength >= 5.9: 

				P_b = (0.2130((inverse_wavelength - 5.9)**2)) + (0.1207((inverse_wavelength - 5.9)**3))

			elif inverse_wavelength < 5.9:

				P_b = 0

			b = (-3.090) + (1.825*inverse_wavelength) + (1.206 / ( ((inverse_wavelength - 4.62)**2) + 0.263) ) + P_b

		elif type_Wavelength == "FUV":

			b = 13.670 + (4.257*(inverse_wavelength - 8)) - (0.420 * ((inverse_wavelength - 8)**2)) - (0.374 * ((inverse_wavelength - 8)**3))

		else: 
			raise NotImplementedError("Incorrect type_Wavelength string")

		return b

	def CCMModel(self, wavelength, R_A, A_V, **kwargs):

		a_lambda = self.a(wavelength, type_Wavelength, **kwargs)
		b_lambda = self.b(wavelength, type_Wavelength, **kwargs)

		A_lambda = ( a_lambda + (b_lambda / R_A) ) * A_V
                return A_Lambda

