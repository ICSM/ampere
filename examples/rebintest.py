import numpy as np
import os
import math
import ampere
from ampere.data import Spectrum, Photometry
from ampere.infer.emceesearch import EmceeSearch
from ampere.models import Model
from spectres import spectres
import pyphot
from emcee import moves
import matplotlib.pyplot as plt
from astropy.io import ascii
from scipy import interpolate

dataDir = os.getcwd() + '/NGC6302/'
specFileExample = 'NGC6302_100.tab'
ngc100 = ascii.read(dataDir+specFileExample,data_start=2)

specFileExample = 'NGC6302.tab'
ngc = ascii.read(dataDir+specFileExample,data_start=2)

input_noise_spec = 0.01 #assume a 1% error on the spectrum
unc = ngc[1][:]*input_noise_spec
spec = Spectrum(ngc[0][:],ngc[1][:] +
                    np.random.randn(len(ngc[1][:]))*unc,ngc[1][:]*0.05,"um","Jy", calUnc=1e-10, scalelengthPrior=0.1)

input_noise_spec = 0.01 #assume a 1% error on the spectrum
unc100 = ngc100[1][:]*input_noise_spec
spec100 = Spectrum(ngc100[0][:],ngc100[1][:] +
                    np.random.randn(len(ngc100[1][:]))*unc100,ngc100[1][:]*0.05,"um","Jy", calUnc=1e-10, scalelengthPrior=0.1) #added extra keyword calUnc with

wavelengths = np.linspace(2.3603,196.6261,1956)

opacityFileName = 'NGC6302-opacities.txt'
opacityDirectory = os.getcwd()+'/NGC6302/'
opacityFileList = np.loadtxt(opacityFileName, dtype='str')
j = 3
tempData = np.loadtxt(opacityDirectory + opacityFileList[j],
                                  comments='#')
tempWl = tempData[:, 0]
tempOpac = tempData[:, 1]
regrid = np.arange(1/max(tempWl),1/min(tempWl),2e-4)
print(regrid)
print(1/tempWl[::-1])
print(tempOpac[::-1])
diop = spectres(regrid,1/tempWl[::-1],tempOpac[::-1])

plt.xscale('log')
#plt.plot(1/wavelengths, f(wavelengths))
plt.plot(tempWl,tempOpac)
plt.plot(1/regrid,diop)
plt.title(opacityFileList[j])
#plt.xlim(2,130)
#plt.ylim(0,5000)
plt.show()




