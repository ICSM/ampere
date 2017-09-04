""" 
Some functions to interface with pyphot
"""

import pyphot #Needed for the synthphot package.
import numpy as np

from astropy.table import Table
from astropy.io import fits #Only used by getGRAMSCgrid at the moment.

#from general import ss_physconst as pc #what is this, and what functions do we need from your "general" module? Why can't it be done with astropy.constants? Oh, I see, it doesn't get called

def synthphot(inlam, influx, filtlist='SPITZER_IRAC_36', fnu=True, **kwargs):
    """
    Function to compute synthetic photometry using pyphot routines.  
    
    For now, assumes that INLAM is in micron
    FILTLIST is the list of filters. The script should check if these exist in the library.

    List of inputs go here with description

    List of outputs go here with description
    """
    lib = pyphot.get_library()
    f = lib['SPITZER_IRAC_36']
    if fnu:
        out = f.get_flux(1e4*inlam[0], influx*(f.lpivot._magnitude/(inlam*1e4))**2)
    else:
        out = f.get_flux(1e4*inlam, influx)
    return out

def getGRAMSCgrid(**kwargs):
    """
    Read in the carbon-star GRAMS grid. 

    This is a test module that need not appear in the fitter.
    """
    hdulist = fits.open('/Users/sundar/work/PRO/sundar_misc/GRAMS/grams_c.fits')
    cgrid = hdulist[1].data; hdulist.close()
    return cgrid

def dostuff(**kwargs):
    """
    Does some stuff (I think?)
    """
    cg = getGRAMSCgrid()
    fphot = synthphot(cg['Lspec'], cg['Fspec'])
    return fphot

def makemeanphot(invot,outfits, **kwargs):
    """

    """
    tab = Table.read(invot+".vot",format="votable")
    u, ui = np.unique(tab["sed_filter"], return_index=True)
    favg = []; e_favg = []; fmed = []; e_fmed = []
    for t in u:
        tab2 = tab[tab["sed_filter"]==t]
        flux = tab2["sed_flux"]
        dflux = tab2["sed_eflux"]
        #All the while ignoring nans
        favg.append(np.nansum(flux)/len(flux[flux!=np.nan]))
        e_favg.append(np.sqrt(np.nansum(dflux)**2/len(dflux[dflux!=np.nan]) + (0.5*(np.nanmax(flux)-np.nanmin(flux)))**2))
        fmed.append(np.nanmedian(flux))
        e_fmed.append(np.nanmedian(np.abs(flux - np.nanmedian(flux))))

    filterNames = tab[ui]["sed_filter"]
    filterNames = np.array(filterNames); favg = np.array(favg); e_favg = np.array(e_favg); fmed = np.array(fmed); e_fmed = np.array(e_fmed)
    tab2 = Table([filterNames, favg, e_favg, fmed, e_fmed], names=['sed_filter','favg','e_favg','fmed','e_fmed'], meta=
                 {'name': 'mean photometry'}, dtype=('S30', np.float32, np.float32, np.float32, np.float32))
    tab2.write(outfits+'_meanphot.fits', format='fits')
    return 1

if __name__ == "__main__":
    print "Main"
