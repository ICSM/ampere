"""Download filter transmission curves from the Spanish Virtual Observatory.
Save these as an hd5 file to be ingested into pyphot.
INPUT: filters.csv, a two-column file containing the name of the magnitude column and the name of the SVO/VOSA filter file to download.
At the moment, the header in the output files only contains the DetectorType specification if the filter
is a photon counter (DetectorType value = "1"). We read in the entire file and check for the occurrence
of this line and set the detector type accordingly.
"""

import numpy as np
import subprocess
from astropy.table import Table
import pyphot as pyp
import h5py
#from astropy import units as u

def makeFilterSet(infile = 'filters.csv', outfile = 'filters.hd5'):
    #infile = 'filters.csv'
    tin = Table.read(infile, format = 'csv', \
                     names = ('column', 'filtername'))
    url = 'http://svo2.cab.inta-csic.es/theory/fps/fps.php?ID='
    filters = []
    fnames = []
    for t in tin:
        filt = t['filtername']
        try:
            g, fname = get_filter_svo(filt, url = url)
        except ValueError:
            try:
                g, fname = get_filter_file(filt)
            except OSError:
                logging.info("Filter %s could not be found on the SVO or in the local directory. \n Please check where to find this filter. \n Continuing to the next filter "%(filt))
                continue
        filters.append(g)
        fnames.append(fname)
    _ = subprocess.call(['rm', 'temp.vot'])
    h = h5py.File(outfile, 'w')
    h.create_group('filters')
    h.close()
    h = pyp.HDF_Library(source = outfile)
    #Don't try to add a filter if it's already in there
    _, u = np.unique([f.name for f in filters], return_index = True)
    for f in list(np.array(filters)[u]):
        h.add_filter(f)




def getFilterList(filterList, outfile='filters.hd5'):
    url = 'http://svo2.cab.inta-csic.es/theory/fps/fps.php?ID='
    filters = []
    fnames=[]
    for filt in filterList:
        try:
            g, fname = get_filter_svo(filt, url = url)
        except ValueError:
            try:
                g, fname = get_filter_file(filt)
            except OSError:
                print("Filter %s could not be found on the SVO or in the local directory. \nPlease check where to find this filter. \nContinuing to the next filter "%(filt))
                continue
        filters.append(g)
        fnames.append(fname)
    _ = subprocess.call(['rm', 'temp.vot'])
    h = h5py.File(outfile, 'w')
    h.create_group('filters')
    h.close()
    h = pyp.HDF_Library(source = outfile)
    #Don't try to add a filter if it's already in there
    _, u = np.unique([f.name for f in filters], return_index = True)
    for f in list(np.array(filters)[u]):
        h.add_filter(f)
    return fnames


def get_filter_svo(filt, url = 'http://svo2.cab.inta-csic.es/theory/fps/fps.php?ID='):
    """
    This function retrieves the filter file specified by the input filt from the URL url.
    
    Parameters
    ----------
    filt (str): The file name or path of the filter file.
    url (str): The URL to fetch the filter file from. Default: SVO url (http://svo2.cab.inta-csic.es/theory/fps/fps.php?ID=)

    Returns
    -------
    g (pyphot.Filter): An object representing the filter.
    fname (str): The file name of the filter file, with '/' replaced by '_'.
    """
    _ = subprocess.call(['curl', '-o', 'temp.vot', url + filt])
    with open('temp.vot') as f:
        content = f.readlines()
    if any("DetectorType" in c for c in content):
        det_type = 'photon'
    else:
        det_type = 'energy'
    #try:
    temp = Table.read('temp.vot', format = 'votable')
    #except ValueError:
    #    print("Filter ",filt," could not be found on the SVO. Please check the filter name. \n Continuing with the next filter")
    #    continue
    fname = filt.replace('/','_')
    g = pyp.Filter(np.array(temp['Wavelength']), np.array(temp['Transmission']), \
                   name = filt.replace('/','_'), unit = temp['Wavelength'].unit.name, \
                   dtype = det_type)
    return g, fname






def get_filter_file(filt, comment='#'):
    """
    This function gets the filter file specified by the input filt.

    Args:
    filt (str): The file name or path of the filter file.
    comment (str): The comment character used in the filter file.

    Returns:
    g (pyphot.Filter): An object representing the filter.
    fname (str): The file name of the filter file, with '/' replaced by '_'.
    """
    header = np.loadtxt(filt + '.csv', max_rows = 2, delimiter = ',', dtype = 'str')
    data = np.loadtxt(filt + '.csv', skiprows = 2, delimiter = ',')
    fname = filt.replace('/','_')
    #Note: dtype in the following line is the detector type of the filter
    g = pyp.Filter(data[:, 0], data[:, 1], name = filt.replace('/', '_'), unit = header[1, 0], \
                   dtype = header[1, 1])
    #filters.append(g)
    return g, fname
    pass


if __name__=="__main__":
    makeFilterSet()
