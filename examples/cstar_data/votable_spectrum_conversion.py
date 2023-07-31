from astropy.io import votable
from astropy.table import Table
from astropy import constants as c
import numpy as np

direc   = '/Users/jonty/mydata/sundar/hyperion/'

spectrum = votable.parse_single_table("OGLE_LMC_LPV_28579_photosphere.vot").to_table()

lam = np.asarray(spectrum['LAM']).flatten()
nu = c.c.value/(lam*1e-6)
flx = np.asarray(spectrum['SPEC']).flatten()

output_file = 'photosphere_original.csv'

data = np.column_stack([nu,flx]) 
np.savetxt(direc+output_file , data, delimiter=',') 



spectrum = votable.parse_single_table("OGLE_LMC_LPV_28579_photosphere_interpolated.vot").to_table()

lam = np.asarray(spectrum['LAMINTERP']).flatten()
nu = c.c.value/(lam*1e-6)
flx = np.asarray(spectrum['SPECINTERP']).flatten()

output_file = 'photosphere_interpolated.csv'

data = np.column_stack([nu,flx]) 
np.savetxt(direc+output_file , data, delimiter=',') 