import numpy as np


def read_eso_filt(filename):
    with open(filename, 'r') as f:
        a = np.loadtxt(f)
        waves = a[:,1]
        trans = a[:,2]
        try:
            sig_trans = a[:,3]
            #print(waves, trans, sig_trans)
            return waves, trans, sig_trans
        except:
            return waves, trans

def eso_filt_to_csv(wave, trans, outfile, sigtrans=None, dettype='energy', unit="um",**kwargs):
    with open(outfile, 'w') as f:
        """ Write the header that Sundar wants """
        outstr = "wavelength, transmission"
        if sigtrans is not None:
            outstr += ", uncertainty"
        outstr += "\n"
        f.write(outstr)
        outstr = unit + ', ' + dettype
        if sigtrans is not None:
            outstr += ", "

        outstr += "\n"
        f.write(outstr)


        """ Now loop over the filter curve and output it """
        for i in range(len(wave)):
            outstr = str(wave[i]) + ', ' + str(trans[i])
            if sigtrans is not None:
                outstr += ', '
                outstr += str(sigtrans[i])
            outstr += "\n"
            f.write(outstr)

if __name__=="__main__":
    filts = ['Paranal/VISIR.SIV_2']#['Paranal_VISIR.J89','Paranal_VISIR.M']

    for f in filts:
        wave, trans, sig_trans = read_eso_filt(f)
        outf = f + ".csv"
        eso_filt_to_csv(wave, trans, outf, sigtrans = sig_trans)
