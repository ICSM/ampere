import numpy as np
import sys
sys.path.insert(1, '/home/zeegers/git_ampere/ampere/')
import ampere
from ampere.data import Spectrum, Photometry
from ampere.emceesearch import EmceeSearch
from ampere.PowerLawAGN import PowerLawAGN, SingleModifiedBlackBody
from ampere.PowerLawAGN import PowerLawAGN, PowerLawAGN
#from extinction import CCMExtinctionLaw
#from extinction import apply, fitzpatrick99
import corner
import matplotlib.pyplot as plt
import os

from astropy import constants as const
from astropy import units as u
from spectres import spectres
import pyphot
from astropy.io.votable import parse_single_table

import pdb

if __name__=="__main__":
    """
    This script should contain
    """
    #ardila_2008=readfits('/home/zeegers/small_silicates_crystal/small_silicates/infrared_spectra/cygob2_12/cassis_yaaar_spcfw_27570176t.fits',ardila_header) 
    #evans_2004_HR=readfits('/home/zeegers/small_silicates_crystal/small_silicates/infrared_spectra/cygob2_12/cassis_yaaar_optdiff_9834496.fits',evans_header) 
    #evans_2004_LR=readfits('/home/zeegers/small_silicates_crystal/small_silicates/infrared_spectra/cygob2_12/cassis_yaaar_spcfw_9834496t.fits',evans_header2)

    """ Read in some data """
    filename = ampere.__file__.strip('__init__.py')+'Testdata/cassis_yaaar_spcfw_27570176t.fits' # path to data of cyg ob 2 12
    print(filename)
    irs = Spectrum.fromFile(filename,format='SPITZER-YAAAR')
    #irs[0].selectWaves(low = 5.4, up = 22.) #following Srinivasan et al. 2017, check if limits are ok
    #irs[1].selectWaves(low = 5.4, up = 22.) #following Srinivasan et al. 2017, check if limiets are ok
    
    #for s in irs:
    #    print(s)
    #exit()
    #If we have more than one dataset we need multiple data objects
    #print(irs.wavelength)

    filename2 = ampere.__file__.strip('__init__.py')+'Testdata/cassis_yaaar_spcfw_9834496t.fits' # path to data of cyg ob 2 12
    
    irs2 = Spectrum.fromFile(filename2,format='SPITZER-YAAAR')
    
    photFile = ampere.__file__.strip('__init__.py')+'Testdata/vizier_votable_cygob212_time.vot'
    
    """ Define the model """
    modwaves = 10**np.linspace(0.,2., 2000)
    #print(modwaves)
    #print(10**modwaves)
    #exit()

    
    ''' Now we get the photometry too *** look for the SED information of Schulte 12 ''' 
    #filters=[
             ##'SPITZER_IRAC_36',
             #'WISE_RSR_W1',
             #'WISE_RSR_W2',
             #'WISE_RSR_W3',
             #'WISE_RSR_W4',
             #'HERSCHEL_PACS_BLUE'
             #]
    #libDir = os.getcwd() + '/ampere/'
#    libname = libDir + 'ampere_allfilters.hd5'
    libname = '/home/zeegers/ampere/ampere/ampere_allfilters.hd5'

    table = parse_single_table(photFile)
    
    table = table.to_table()    
    #data = table.array
    
    # https://numpy.org/doc/stable/reference/generated/numpy.allclose.html
    # Sascha: read in real photometry data ....
    desired_filters=[b'2MASS:J', b'2MASS:H', b'2MASS:Ks', b'Spitzer/MIPS:24'] #these are the filters we're after
    mask = np.isin(table['sed_filter'], desired_filters) #np.isin() is true for each element of table['filter'] that matches one of the elements of desired_filters
    phot_new = table[mask] #now we make a new table which is just the rows that have the filters we want
    desired_filters_again=[b'II/328/allwise']
    mask_again = np.isin(phot_new['_tabname'], desired_filters_again)
    #phot_table = table([phot])
    phot_again = phot_new[mask_again]
    
    #table.write('new_table.vot', format='votable')
    phot_again.write('new_table.vot', format='votable', overwrite=True)
    
    
    
    
    phot = Photometry.fromFile('new_table.vot', libName = libname) # gaat hier nog steeds iets niet goed
    phot.selectWaves(low = 35., interval = "right-open") #using only MIPS-70 and PACS, following Srinivasan et al. 2017

    print('did we get here?')

    phot.reloadFilters(modwaves)
    dataSet = [phot]          #use this line when using photometry
    print(phot)

    for s in irs:             #include the next two lines when appending spectroscopy to photometry
        dataSet.append(s)

    for s in irs2:             #include the next two lines when appending spectroscopy to photometry
        dataSet.append(s)

    ''' now turn the numbers into a Photometry object '''
#    phot = Photometry(np.array(filters), np.array(sed), np.array(0.1*sed), np.array(['Jy', 'Jy','Jy','Jy','Jy']), bandUnits='um',libName = libname)



    #model = PowerLawAGN(irs.wavelengths) #We should probably define a wavelength grid separate from the data wavelengths for a number of reasons, primarily because the photometry needs an oversampled spectrum to compute synthetic photometry

    model = PowerLawAGN(modwaves, flatprior=True, redshift=None)
                                    
                                                                        
    print(model.npars)                                
    
    """ and the extinction law (if necessary) (could be moved inside model as only some models will need this) """
    #ext = CCMExtinctionLaw()

    #model.extinction = ext #combination of multiple different models is something we need to consider how to implement. 

    """ Connect the dots: """
    """ Hook it up to an optimiser """
    opt = EmceeSearch(model = model, data = dataSet, nwalkers = 100) #introspection is used inside the optimiser to determine the number of parameters that each part of the model has
    
    print(opt.npars)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in irs:
        ax.loglog(i.wavelength, i.value, '-',color='blue')
    for j in irs2:
        ax.loglog(j.wavelength, j.value, '-',color='blue')        
    ax.loglog(phot.wavelength, phot.value, 'o',color='blue')
    #ax.set_ylim(0., 1.5*np.max([np.max(i.value) for i in dataSet]))
    #fig.savefig("sed_test.png")

    model(-2.,-0.5,
                    -0.5,-0.5,-0.5,
                    -0.5,-0.5)
    
    modspec = model.modelFlux
    print(modspec)
    ax.plot(modwaves,modspec, color="red")
    plt.show() #this plots the spectrum and photometry plus the shape of the model SED using the input parameters


    """ if you want, overload the default priors with a new function """
#    def lnprior(self, inputs, **kwargs):
#        return 0 #flat prior
    
    #model.lnprior = lnprior
    
    #20 + np.random.randn() for i in range(np.int(opt.npars))
    
    pdb.set_trace()
    print(opt.npars,np.int(opt.npars))
    """ Run it """
    pos = [
           [
               -2., -0.78, -0.78, -0.78, -0.78, -0.78, -0.78, 1., 0.5, 1., 1., 0.5, 1., 1., 0.5, 1.
               #1000., 17., -1.0, 2., 1., 0.5, 1., 1., 0.5, 1., 1., 0.5, 1.
               #1000., 17., -1.0, 2., 0., 0.,  0., 0., 0., 0.
               #100. + np.random.randn() for i in range(np.int(opt.npars))
           ]
           + np.random.randn(np.int(opt.npars)) for j in range(opt.nwalkers)
          ]
    print(pos[0])
    print(np.max(pos, axis=0))
    print(np.min(pos, axis=0))
    #pdb.set_trace()
    
    opt.optimise(nsamples = 1000, burnin=100,guess=pos)
    lnprob = opt.sampler.lnprobability

    """save optimiser state for later """
    #opt.save(filename="output_file",pickle=True) #Save as a python object
    
    """ Sort out the results """
    """ First produce some terminal output for the numbers and confidence intervals """

    """ Then produce a couple of plots """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    print(opt.samples.shape)
    print(np.max(opt.samples, axis=0))
    print(np.min(opt.samples, axis=0))
    nneg=0
    for i in range(0,opt.samples.shape[0],100):
        #print(opt.samples[i,:])
        if opt.samples[i,0] > 0.:
            opt.model(opt.samples[i,0],opt.samples[i,1],opt.samples[i,2],opt.samples[i,3])
            ax.plot(modwaves,opt.model.modelFlux, '.', 'k', alpha=0.02) #irs.wavelength
        else:
            nneg += 1
            
    #    ax.plot(irs.wavelength, model(opt.samples[i,:]), '-', alpha = 0.1)
    print(nneg)
    for i in irs:
        ax.plot(i.wavelength, i.value, '.','b')
        ax.plot(phot.wavelength, phot.value, 'o',color='b') 
    for i in irs2:
        ax.plot(i.wavelength, i.value, '.','b')
    fig2 = corner.corner(opt.samples)#,labels=opt.labels)
    plt.show()

    """ Figure out what it all means """
