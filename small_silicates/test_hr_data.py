import sys
sys.path.insert(1, '/home/zeegers/git_ampere/ampere/')
import numpy as np
import os
import ampere
from ampere.data import Spectrum, Photometry
from ampere.emceesearch import EmceeSearch
from ampere.models import Model
from spectres import spectres
import pyphot
from emcee import moves
from astropy.modeling import models
from astropy import units as u
from astropy.modeling.models import BlackBody
import matplotlib.pyplot as plt
from scipy import interpolate

from numpy import ma
from astropy.table import Table
from astropy.table import vstack
from astropy import constants as const
import sys
sys.path.insert(1, '/home/zeegers/git_ampere/ampere/')
import numpy as np
import os
import ampere
from ampere.data import Spectrum, Photometry
from ampere.emceesearch import EmceeSearch
from ampere.models import Model
from spectres import spectres

from emcee import moves
from astropy.modeling import models
from astropy import units as u
from astropy.modeling.models import BlackBody
import matplotlib.pyplot as plt
from scipy import interpolate
from astropy.constants import c

import pyphot
from astropy.io import fits
# for reading VO table formats into a single table
from astropy.io.votable import parse_single_table
from astropy.io import ascii
from scipy.stats import rv_continuous, norm, halfnorm
#from scipy.linalg import inv
from numpy.linalg import inv

from astropy.table import Table
from astropy.table import vstack

dataDir = ampere.__file__.strip('__init__.py') + 'Testdata/'
specFile3 = 'cassis_yaaar_optdiff_9834496.fits'

filename = dataDir+specFile3


hdul = fits.open(filename)
hdu = hdul[0]
header=hdu.header
data = hdu.data
table = Table(data,names=[header['COL01DEF'],header['COL02DEF'],header['COL03DEF'],header['COL04DEF'],header['COL05DEF'],header['COL06DEF'],header['COL07DEF'],header['COL08DEF']])
table['wavelength'].unit='um'
table['flux'].unit='Jy'
table['flux_error'].unit='Jy'
table.rename_column('flux_error','uncertainty')
#table['error (RMS)'].unit='Jy'
#table['error (SYS)'].unit='Jy'
 #table['offset uncertainty (CAL)'].unit='Jy'
#table['sky'].unit='Jy'
#table['sky error'].unit='Jy'
mask = np.logical_and.reduce([np.isfinite(c) for c in table.columns.values()]) #require the elements to be non-NaNs

table = table[mask]
chunks = np.zeros_like(table['module'].data)

wavelengths = np.array(table['wavelength'].data)
orders = np.array(table['IRS_order'].data)
diff_wav=np.diff(wavelengths)
split=np.where(diff_wav > 5)
splits=np.asarray(split[0])
order_array=np.arange(11,21,1)


#for i in range(len(wavelengths)-1):
    #for j in range(0,10):
        ##print(i,j)
        
        #if((orders[i]==order_array[j]) and (i<splits[j])):
            #print(i)
            ##chunks[i] = 1.
            ##print(i)
        #elif ((orders[i]==order_array[j]) and (i>=splits[j])):
            #chunks[i] = 2.
            ##print(i)
            

chunks2=np.zeros(len(wavelengths))
for i_order,order in enumerate(order_array):
    wl_sel=wavelengths[order==orders]
    nums=np.where(order==orders)[0]
    print(nums)
    print(wl_sel)
    print((np.diff(wl_sel)))
    print((np.diff(wl_sel)>5))
    breaknum=np.where(np.diff(wl_sel)>5)[0][0]+1
    print(breaknum)
    chunks2[nums[0:breaknum]]=1
    chunks2[nums[breaknum::]]=2
plt.plot(wavelengths,chunks2,"o")
plt.show()
    
#teller gaat alleen omhoog als orders[i]-orders[i-1]>0.
    

teller=0
for i in range(1, len(wavelengths)):
    
    if((orders[i-1]==order_array[teller]) & (i-1<=splits[teller])):
        chunks[i-1]=1.
    elif((orders[i-1]==order_array[teller]) & (i-1>splits[teller])):
        chunks[i-1]=2.
    if(i==len(wavelengths)-1.):
        chunks[i]=2.
        
    if((orders[i]-orders[i-1]>=0.5)):
        teller=teller+1
        #print(i)
        #print(teller)
        
plt.plot(wavelengths,chunks,"o")
plt.show()


#for i in range(0,10):
    #chunks[sh] = 
    #chunks = 1
    #chunks = 2



table['chunk'] = chunks    


    #sh = 
    #lh = 
    ##should we divide this in four chunks: module 0 = SL2, module 1 = SL1, module 2 = LL2, module 3 is LL1?
    #chunks[sh] = 1.
    #chunks[lh] = 2.
    #table['chunk'] = chunks

    ##normalise the column names (wavelength,flux,uncertainty)
    #table.rename_column('flux_error','uncertainty')
                
    #''' If I'm interpreting the CASSIS data correctly, Module(SL) = 0, 1; Module(LL) = 2, 3 '''
    ##a = table['module'] > 1.
    ##tablell = table[a]#.sort(keys='wavelength')
    ##tablell.sort(keys='wavelength')
    ##tablell.pprint()
    ##table = tablell
    #table.sort(keys='wavelength')            
