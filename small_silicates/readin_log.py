# this program reads in the .log file from ampere and recreates the plotting from the program 

import sys
sys.path.insert(1, '/home/zeegers/git_ampere/ampere/')
import numpy as np
import os
import ampere
from ampere.data import Spectrum, Photometry
from ampere.infer.emceesearch import EmceeSearch
from ampere.models import Model
from spectres import spectres
import pyphot
from emcee import moves
from astropy.modeling import models
from astropy import units as u
from astropy.modeling.models import BlackBody
import matplotlib.pyplot as plt
from scipy import interpolate
from astropy.constants import c
from astropy.io.votable import parse_single_table
import pdb
import re


class Blackbody_dust(Model): 
    '''This is a blackbody model based on the very simple model from Peter 

    We will first try to fit a blackbody model to a dataset with dust 

    '''
    def __init__(self, wavelengths, flatprior=True,
                 lims=np.array([[10., 300000.],
                                [1.,6.],
                                [0.001,10.]])):
        '''The model constructor, which will set everything up

        This method does essential setup actions, primarily things that 
        may change from one fit to another, but will stay constant throughout 
        the fit. This may be things like the grid of wavelengths to calculate
        the model output on, or establishing the dust opacities if involved.
        There are also several important variables it *MUST* define here
        '''
        self.wavelength = wavelengths
        
        # Getting the opacities from the folder 
        #opacityDirectory = os.path.dirname(__file__)+'/Opacities/'
        opacityDirectory = os.path.dirname(os.path.realpath('__file__'))+'/optical_const_bulk/'
        print("Directory:", opacityDirectory)
        opacityFileList = os.listdir(opacityDirectory)
        opacityFileList = np.array(opacityFileList)[['.q' in zio for zio in opacityFileList]] # Only files ending in sub.q are valid (for now). At the moment there are 6 files that meet this criteria
        print(opacityFileList)
        nSpecies = opacityFileList.__len__()
        opacity_array = np.zeros((wavelengths.__len__(), nSpecies))
        
        
        for j in range(nSpecies):
            tempData = np.loadtxt(opacityDirectory + opacityFileList[j], comments = '#')
            print(opacityFileList[j])
            tempWl = tempData[:, 0]
            tempOpac = tempData[:, 1]            
            

            f = interpolate.interp1d(tempWl, tempOpac, assume_sorted = False)
            opacity_array[:,j] = f(wavelengths)#wavelengths)
            
        self.opacity_array = opacity_array
        self.nSpecies = nSpecies
        
        self.npars = nSpecies + 3 #Number of free parameters for the model (__call__()). For some models this can be determined through introspection, but it is still strongly recommended to define this explicitly here. Introspection will only be attempted if self.npars is not defined.
        self.npars_ptform = nSpecies + 3 #Sometimes the number of free parameters is different when using the prior transform instead of the prior. In that case, self.npars_ptform should also be defined.
        #You can do any other set up you need in this method.
        #For example, we could define some cases to set up different priors
        #But that's for a slightly more complex example.
        #Here we'll just use a simple flat prior
        self.lims = lims
        self.flatprior = flatprior
        labels = ["Temperature", "radius star", "Scaling parameter"]
        labels2 = labels+opacityFileList.tolist()
        print("Labels:", labels2)
        self.parLabels = labels2

    def __call__(self, temp, radius_sol, scaling, *args, **kwargs):
        '''The model itself, using the callable class functionality of python.

        This is an essential method. It should do any steps required to 
        calculate the output fluxes. Once done, it should stop the output fluxes
        in self.modelFlux.
        '''
        dustAbundances = np.array(args) # instead of 10**np.array(args)
        fModel = (np.matmul(self.opacity_array, dustAbundances))
        fModel2 = np.exp(-fModel)
        
        #plt.plot(wavelengths,self.opacity_array[:,0], label = "0")
        #plt.plot(wavelengths,self.opacity_array[:,1], label = "1")
        #plt.plot(wavelengths,self.opacity_array[:,2], label = "2")
        #plt.plot(wavelengths,self.opacity_array[:,3], label = "3")
        #plt.plot(wavelengths,self.opacity_array[:,4], label = "4")
        #plt.plot(wavelengths,self.opacity_array[:,5], label = "5")
        #plt.legend()
        #plt.show()

        wavelengths_aa = (self.wavelength*u.micron).to(u.AA)
        bb =  BlackBody(temperature=temp*u.K)
        pc  = 3.086e16 # m
        rsol = 696340e3 # m
        
        distance_pc=1100.
        
        Rstar = radius_sol*rsol
        r = distance_pc*pc
        # from angstrom to micron 
        flux = bb(wavelengths_aa)
        flux_freefree = (((c.to(u.cm/u.s))/wavelengths_aa.to(u.cm))**0.6).value
        flux_mjy = (flux.to(u.mJy / u.sr).value*(Rstar/r)**2.+scaling*1.e-7*flux_freefree)*fModel2
        # here we need to put dust models ...
        dustAbundances = np.array(args) # instead of 10**np.array(args)
        
        self.modelFlux = flux_mjy #slope*self.wavelength + intercept 
        
def tracking_logfile(string1):
    # opening the log file
    file1 = open("ampere.log", "r")

    # setting flag and index to 0
    flag = 0
    teller = 0
    index = 0
    index_array = np.array([])

    # Loop through the file line by line
    for line in file1:  
        index = index + 1 
        
        # checking string is present in line or not
        if string1 in line:
            index_array=np.append(index_array, index)
            #teller = teller + 1
    
    # checking condition for string found or not
    #if flag == 0: 
        #print('String', string1 , 'Not Found') 
    #else: 
        #print('String', string1, 'Found In Line', index)
        #index_array[teller]
        #teller = teller + 1
    
    # closing text file    
    file1.close() 
    
    return(index_array)

def finding_number(string_no):
    x = string_no
    f = open('ampere.log')
    for line in f:
        if x in line:
            print(line[line.find(x)+len(x):])
            parameter_string=line[line.find(x)+len(x):]

    emp_lis = []
    for z in parameter_string.split():
        if z.isdigit():
            emp_lis.append(int(z))
      
    f.close()
    return(emp_lis)


# first receate the plotting of the data 

wavelengths = np.linspace(1.0,38, 1000)

pc  = 3.086e16 # m
rsol = 696340e3 # m

""" Choose some model parameters """

#Now init the model:

# Spitzer Spectra 

dataDir = ampere.__file__.strip('__init__.py') + 'Testdata/'
specFile1 = 'cassis_yaaar_spcfw_27570176t.fits'
irsEx_1 = Spectrum.fromFile(dataDir+specFile1,format='SPITZER-YAAAR')
    
specFile2 = 'cassis_yaaar_spcfw_9834496t.fits'
irsEx_2 = Spectrum.fromFile(dataDir+specFile2,format='SPITZER-YAAAR')
    
specFile3 = 'cassis_yaaar_optdiff_9834496.fits'
irsEx_3 = Spectrum.fromFile(dataDir+specFile3,format='SPITZER-YAAAR_OPTDIFFHR')

# Photometry 

libname = '/home/zeegers/ampere/ampere/ampere_allfilters.hd5'

photFile = ampere.__file__.strip('__init__.py')+'Testdata/vizier_votable_cygob212_time_again.vot'
table1 = parse_single_table(photFile)
   
table1 = table1.to_table()    

    # Sascha: read in real photometry data ....
desired_filters=['2MASS:J', '2MASS:H', '2MASS:Ks', 'Spitzer/MIPS:24', 'WISE:W4','WISE:W3'] #these are the filters we're after
mask = np.isin(table1[u'sed_filter'], desired_filters) #np.isin() is true for each element of table['filter'] that matches one of the elements of desired_filters
phot_new = table1[mask] #now we make a new table which is just the rows that have the filters we want
desired_filters_again=['II/328/allwise','I/ApJS/191/301/table1']
mask_again = np.isin(phot_new['_tabname'], desired_filters_again)
phot_again = phot_new[mask_again]
    
# I don't think you really need to write the votable again though... so for now this is commented     
# phot_again.write('new_table.vot', format='votable', overwrite=True)

phot = Photometry.fromFile('new_table.vot', libName = libname)
    
phot.reloadFilters(wavelengths)
    
dataSet = [phot]

fig = plt.figure()    
ax = fig.add_subplot(111)
ax.set_yscale('log')
ax.set_xscale('log')

#ax.plot(irsEx_1[0].wavelength, irsEx_1[0].value, '-',color='red')
#ax.plot(irsEx_1[1].wavelength, irsEx_1[1].value, '-',color='red')        
#ax.plot(irsEx_2[0].wavelength, irsEx_2[0].value, '-',color='red')
#ax.plot(irsEx_3[0].wavelength, irsEx_3[0].value, '-',color='blue')
    #ax.plot(irsEx_3[1].wavelength, irsEx_3[1].value, '-',color='blue')
    
#ax.plot(wavelengths, model_flux)
ax.errorbar(phot.wavelength, phot.value, yerr = phot.uncertainty,marker='o',color='green', ls='none')
#ax.plot(phot.wavelength, phot.value, 'o',color='green')

#Need to know how many parameters there are
#Need to put the names of these parameters in an array
#when the names are recognized after the line with the MAP solution they should be read in


# string to search for 
string1 = 'MAP solution:'
string2 = 'This model has'
string3 = 'There are also'

index1 = tracking_logfile(string1) # line number of the MAP solution 

emp_lis = finding_number(string2)

print("Find number in string:",emp_lis)

parameter_array=np.empty(emp_lis[0])

lines_to_print = np.arange(emp_lis[0])+index1[len(index1)-1] 

file_again = open('ampere.log')

teller = 0 

for index, line_again in enumerate(file_again):
    if (index in lines_to_print):
        print(line_again)
        print('index', index)
        numbers = (re.findall(r"[-+]?(?:\d*\.\d+|\d+)", line_again))
        parameter_array[teller] = numbers[-1]
        teller=teller+1

file_again.close()        

# now we need to fill the parameter 

""" get the model parameters from the log file """

model = Blackbody_dust(wavelengths)

model(parameter_array[0], parameter_array[1], parameter_array[2], parameter_array[3], parameter_array[4], parameter_array[5], parameter_array[6], parameter_array[7], parameter_array[8], parameter_array[9])

model_flux = model.modelFlux

ax.plot(wavelengths, model_flux)

parameter_names=model.parLabels # we also need an array with all the parameter names 

# now we need to overplot also the shifted data 

# x10, x13 and x16 give the shifts in the datasets. Find them in the dataset and multiply by this number 

# take the line number of the MAP, the number of parameters and search how many parameters there are for the noise model
# Every first out of three is the scaling parameter 

noise_parameters = finding_number(string3)

scaling_parameters = int(noise_parameters[0]/3)

index_noise_param = np.empty(scaling_parameters)


for i in range(len(index_noise_param)):
    print(i)
    index_noise_param[i]= emp_lis[0]+index1[len(index1)-1]+i*3

file_again2 = open('ampere.log')

teller = 0 

noise_param = np.empty(scaling_parameters)

for index2, line_again2 in enumerate(file_again2):
    if (index2 in index_noise_param):
        print(line_again2)
        print('index', index2)
        numbers2 = (re.findall(r"[-+]?(?:\d*\.\d+|\d+)", line_again2))
        noise_param[teller] = numbers2[-1]
        teller=teller+1

file_again2.close() 

for i in dataSet:
    ax.errorbar(irsEx_1[0].wavelength, irsEx_1[0].value*noise_param[0], yerr=irsEx_1[0].uncertainty, color='red', label='Spitzer spectrum')
    ax.errorbar(irsEx_1[1].wavelength, irsEx_1[1].value*noise_param[1], yerr=irsEx_1[1].uncertainty,color='red')    
    ax.errorbar(irsEx_2[0].wavelength, irsEx_2[0].value*noise_param[2], yerr=irsEx_2[0].uncertainty,color='red')
    #ax.plot(irsEx_3[0].wavelength, irsEx_3[0].value, '-',color='blue')
    
ax.set_xlabel(r"Wavelength ($\mu$m)")
ax.set_ylabel(r"Flux density (mJy)")

plt.show()
