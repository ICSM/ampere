#!/usr/bin/env python
"""Download filter transmission curves from the Spanish Virtual Observatory and write into an HDF file."""
__author__ = "Sundar Srinivasan"
__copyright__ = "Copyright 2021"

import numpy as np
import subprocess
from astropy.table import Table
import pyphot as pyp
import h5py

def makeFilterSet(infile = 'filters.csv', outfile = 'filters.hd5'):
    """Download filter transmission curves from the Spanish Virtual Observatory and write into an HDF file.

    Keyword arguments:
    infile -- (default 'filters.csv')
        name of two-column CSV file, with the first column containing a unique name for the filter,
        and the second column containing the name of the filter as specified on the SVO/VOSA webpage. For
        example, a row "IRAC45, Spitzer/IRAC.I2" in infile will download the IRAC 4.5 µm filter curve.

        If the name in the second column does not match any filters on the SVO database, the code will
        search the local machine. This allows the user to incorporate filter files not available on the SVO.
        In this case, the filter profile must be a CSV, with two header lines providing the column content
        as well as the units for the wavelength, and the detector type. For example:
                Wavelength, Transmission
                um, energy
                0.5, 0.00
                0.6, 0.10

    outfile -- (default 'filters.hd5') name of output HDF file.

    Notes: The current SVO data only specify the DetectorType if the filter is a photon counter, by setting
    DetectorType = "1". This script reads in the entire file and checks for the occurrence of this line
    to set the detector type in the output file accordingly. The SVO has said that they will improve
    the implementation of this keyword in the future.
    """

    tin = Table.read(infile, format = 'csv', \
                     names = ('column', 'filtername'))
    url = 'http://svo2.cab.inta-csic.es//theory/fps3/fps.php?ID='
    filters = []
    for t in tin:
        try:
            #Look for the filter in the SVO database
            _ = subprocess.call(['curl', '-o', 'temp.vot', url + t['filtername']])
            with open('temp.vot') as f:
                content = f.readlines()
            if any("DetectorType" in c for c in content):
                det_type = 'photon'
            else:
                det_type = 'energy'
            temp = Table.read('temp.vot', format = 'votable')
            g = pyp.Filter(np.array(temp['Wavelength']), np.array(temp['Transmission']), \
                           name = t['filtername'].replace('/','_'), unit = temp['Wavelength'].unit.name, \
                           dtype = det_type)
            logging.info("Generated filter %s from Spanish Virtual Observatory"%(t))
        except:
            logging.info("Generating filter %s from Spanish Virtual Observatory failed, searching for local file"%(t))
            #Look for the filter on the local machine
            header = np.loadtxt(t['filtername'] + '.csv', max_rows = 2, delimiter = ',', dtype = 'str')
            data = np.loadtxt(t['filtername'] + '.csv', skiprows = 2, delimiter = ',')
            #Note: dtype in the following line is the detector type of the filter
            g = pyp.Filter(data[:, 0], data[:, 1], name = t['filtername'].replace('/', '_'), unit = header[1, 0], \
                           dtype = header[1, 1])
        filters.append(g)
    _ = subprocess.call(['rm', 'temp.vot'])
    h = h5py.File(outfile, 'w')
    h.create_group('filters')
    h.close()
    h = pyp.HDF_Library(source = outfile)
    #Don't try to add a filter if it's already in there
    _, u = np.unique([f.name for f in filters], return_index = True)
    for f in list(np.array(filters)[u]):
        h.add_filter(f)





