import numpy as np
import matplotlib.pyplot as plt
from astropy.constants import c
import astropy.units as u

def make_alma_filter(windows, unit = "GHz",**kwargs):
    """ 
    Make an effective filter transmission curve for a set of ALMA spectral windows, ignoring the effects of the atmosphere

    windows is an array of start and end points of the spectral windows, with shape 2x n_windows, in GHz or micron.
    
    """

    """ Find the outer bounds of the whole set of spectral windows and expand the coverage a little to make sure we have a few zeroes outside the windows """

    if unit != "GHz":
        raise NotImplementedError("Use of units other than GHz hasn't been implemented yet!")
    print(windows)
    start_freq = 0.95 * np.min(windows) 
    end_freq = 1.05 * np.max(windows)
    print(start_freq, end_freq)

    freqs = np.linspace(start_freq, end_freq, 200)
    print(freqs)
    trans = np.zeros_like(freqs)
    for window in windows:
        print(window)
        winmask = np.logical_and(freqs > np.min(window),
                                 freqs < np.max(window)
                                 )
        trans[winmask] = 1.
    return freqs, trans
    


if __name__=="__main__":

    #  GHz -- 8 spws!
    windows_B10 = [[849.3, 851.2],[851.3, 853.3],[853.3, 855.2],[855.3, 857.2],[865.3, 867.2],[867.3, 869.2],[869.3, 871.2],[871.3, 873.2]]
    #  GHz -- 8 spws!
    windows_B9 = [[659,661],[661, 663],[663, 665],[665,667],[675,677],[677, 679],[679, 681],[681,683]]
    # [397,399],[399,401],[409,411],[411,413] GHz
    windows_B8 = [[397,399],[399,401],[409,411],[411,413]]
    # 335.5 -- 337.5, 337.5 -- 339.5, 347.5 -- 349.5, 349.5 -- 351.5 GHz
    windows_B7 = [[335.5, 337.5],[337.5, 339.5],[347.5, 349.5],[349.5, 351.5]]
    #213--214, 214--216, 228--230, 230--232 GHz.
    windows_B6 = [[213,214],[214,216],[228,230],[230,232]]
    #195.05--196.95, 197.05--198.95, 207.05--208.95, 209.05--210.95
    windows_B5 = [[195,197],[197,199],[207,209],[209,211]]
    #137.05--138.95, 139.05--140.95, 149.05--150.95, 151.05--152.95
    windows_B4 = [[137,139],[139,141],[149,151],[151,153]]
    #85--87, 87--89, 97--99, 99--101 GHz.
    windows_B3 = [[85,87],[87,89],[97,99],[99,101]]
    #35--37, 37--39, 39--41, 41--43 GHz
    windows_B1 = [[35,37],[37,39],[39,41],[41,43]]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    filternames = ['ALMA/ALMA.B10','ALMA/ALMA.B9','ALMA/ALMA.B8','ALMA/ALMA.B7','ALMA/ALMA.B6','ALMA/ALMA.B5','ALMA/ALMA.B4','ALMA/ALMA.B3','ALMA/ALMA.B1']
    bands = [windows_B10, windows_B9, windows_B8, windows_B7, windows_B6, windows_B5, windows_B4, windows_B3, windows_B1]

    for i in range(len(bands)):
        #print(band)
        band = bands[i]
        filtername = filternames[i]
        filt = make_alma_filter(band)
        ax.plot(filt[0],filt[1], '-')
        freq = filt[0]
        trans = filt[1]
        wave = c.to(u.um/u.s).value / ((freq*u.GHz).to(u.Hz)).value
        print(wave)
        wave = wave[::-1]
        trans = trans[::-1]
        with open(filtername+".csv", "w") as f:
            f.write("wavelength, transmission\n")
            f.write("um, energy\n")
            for i in range(len(wave)):
                f.write(str(wave[i])+","+str(trans[i])+"\n")
        
    plt.show()
