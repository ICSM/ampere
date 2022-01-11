from astropy.modeling import models
from astropy import units as u
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from extinction import apply, fitzpatrick99
from scipy import interpolate
import math
import glob
from astropy.io import fits
from spectres import spectres
from astropy.modeling.models import BlackBody

# constants
pc  = 3.086e16 # m
rsol = 696340e3 # m


#Define stellar spectrum: blackbody, star with 10 solar radii and at 1 kpc 
R = 3.9*rsol
r = 2000.*pc

bb = BlackBody(temperature=15000*u.K)
wav = np.arange(100,1000000,10) * u.AA     # wavelength from 0.1 to 100 micron
flux = bb(wav)

# transform flux to mJy
flux_mjy = flux.to(u.mJy / u.sr)*(R/r)**2.


#
# Here we want a generalised function to read in a stellar spectrum from the subfolder.
#
#B0V T = 30000 K, log g = 4, Rstar = 10.0 Rsol
#B5V T = 15200 K, log g = 4, Rstar = 3.9 Rsol
#B9V T = 10600 K, log g = 4, Rstar = 2.7 Rsol
#B2V T = 20600 K, log g = 4, Rstar = 3.85 Rsol
#O8V T = 35000 K, log g = 4, Rstar = 8.75 Rsol



# get list of stellar spectra filenames from subfolder

# print list of filenames

dstar = 2000. # pc

# read in stellar spectrum

wav_aa = wav
wav_um = wav/1e4


spectrum_total_int = 0

# Set the depletion of silicon
depletion_si = 0.98


#
# Here we want a function to select the dust component(s) we want to add to the extincted stellar model
# 
def read_dust_qext(species='mg05fe05sio3_am_cde01.q'): 
    # read in the Q values 
    data_dust=np.loadtxt(species)

    wavelengths_q=data_dust[:,0]
    Qext=data_dust[:,1]

    return wavelengths_q, Qext

# Here we read the properties, e.g. density and molmass of the dust species 

def read_dust_prop(species='mg05fe05sio3_am_cde01.q'): 
    # read in the Q values 
    data_dust=np.loadtxt(species)
    sample = open('./dust_qext/silicate_sample.txt','r')
    linessample=sample.readlines()
    sample_column_number = 0
    resultsample=[]
    for x in linessample:
        resultsample.append(x.split()[sample_column_number])
    sample.close()
    #print(resultsample)
    
    index = resultsample.index(species)
    name,density,molmass=linessample[index].split()

    return density, molmass
#
# Calculate number of atoms along line of sight based on extinction
#
def calculate_Nd(Av,rho,gmol,fraction, depletion_si):
    NH=1.9e21*Av
    si_abundance=3.4e-5
    # total amount of silicon atoms available
    si_los=NH*si_abundance*depletion_si
    si_los_specie=si_los*fraction
    # calculate the silicon atoms per particle
    #rho=3.2 # check per composition
    rho=rho
    V=4./3.*math.pi*(0.1*1e-4)**3
    grams=V*rho
    gmol=gmol
    #gmol=100. # check per composition
    mol=grams/gmol
    Avogadro=6.0221409e+23
    si_atoms_per_particle=mol*Avogadro
    Nd=si_los/si_atoms_per_particle
    
    # this needs to be adapted, because it will be different for different species
    # also has the assumption that all the silicon is locked up in dust 
    
    print('Nd',Nd)
    return Nd


#generate noise for a spectrum
def spec_noise(spectrum,noise_para=[0.0,1.0],noise_type='Poisson'):
    
    if noise_type == 'Poisson':
        noise = np.zeros((len(spectrum)))
                
        for k in range(len(noise)):
            noise[k] = 1e-3*np.sqrt(1e3*spectrum[k])*np.random.uniform(-1,1)

    if noise_type == 'Gaussian':
        mu=noise_para[0]
        sigma = noise_para[1]
        noise = np.zeros((len(spectrum)))

        for k in range(len(noise)):
           noise[k] = spectrum[k]*np.random.normal(mu,sigma)

    return noise

# get list of stellar spectra filenames from subfolder

dust_list = glob.glob("./dust_qext/*.q")
print(dust_list)
# print list of filenames
for i in range(len(dust_list)):
    print(i,dust_list[i])

# ask user for desired spectrum, distance

dust_index  = "2,10" #input("Please enter index(s) of desired dust component(s), separated by commas (no spaces e.g. '1,2,3'):")

dust_split = dust_index.split(',')

species=[]
species_string = ""
for i in range(len(dust_split)):
    species.append(dust_list[int(dust_split[i])])
    species_string += species[i][12:]+"_"
Nspecies = len(species)

print("Using ",Nspecies," dust species to generate features on reddened stellar spectrum.")

##################################
dust_fraction  = "0.9,0.1" #input("Please enter the fractions(s) of the desired dust component(s), separated by commas (no spaces e.g. '0.2,0.2,0.6', must add up to 1):")

dust_split_fraction = dust_fraction.split(',')
fraction=np.empty(len(dust_split_fraction))

for i in range(len(dust_split_fraction)):
    fraction[i]=float(dust_split_fraction[i])

# Here we define the extinction properties of the ISM
#Avs = [1.0,3.0,5.0,6.5,8.0,10.0,15.0] # in mags
#Avs = [4.0,5.5,7.0] # in mags
Avs = [7.0] # in mags
Av = 7.0

Rv = 3.1 # for MW

# create list of filenames based on material composition and extinction
filename_array = []
for i in range(len(Avs)): 
    filename_array.append("extinction_"+species_string+"av"+str(Avs[i])+"_"+str(int(dstar/1000))+"kpc.txt")

#filename_array=["extinction_enstatite_av1_1kpc.txt","extinction_enstatite_av3_1kpc.txt","extinction_enstatite_av5_1kpc.txt","extinction_enstatite_av8_1kpc.txt","extinction_enstatite_av10_1kpc.txt","extinction_enstatite_av15_1kpc.txt"]

spectrum = flux_mjy
    
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel(r'Wavelength [$\mu$m]',fontsize=16)
ax.set_ylabel(r'Flux density [Jy]',fontsize=16)
ax.set_ylim(1e-5,1e3)
ax.set_xlim(1e-2,1e3)
ax.loglog(wav_um,spectrum,color='blue')
ax.loglog(wav_um,flux_mjy,color='red')
plt.show()
plt.close()

print(np.min(wav_um),np.max(wav_um))
# wavelength for interpolation
wavelength_int=np.arange(0.1,40.0,0.001)
interpfunc_spectrum=interpolate.interp1d(wav_um,spectrum)
spectrum3_int=interpfunc_spectrum(wavelength_int)

interpfunc_model=interpolate.interp1d(wav_um,flux_mjy)
model_int=interpfunc_model(wavelength_int)
# optical depth array

tau_cont=np.log(spectrum3_int/model_int)*-1.
wavelength_int_q=np.arange(2,40.0,0.001)
tau_int_func=interpolate.interp1d(wavelength_int,tau_cont)
tau_int=tau_int_func(wavelength_int_q)    
model_int2=interpfunc_model(wavelength_int_q)
    
tau_lab=np.empty([len(species),len(wavelength_int_q)])

teller=0

for specie in species:
        
        #def read_dust_qext(species='mg05fe05sio3_am_cde01.q'): 
     ## read in the Q values 
    #data_dust=np.loadtxt(species)

    #wavelengths_q=data_dust[:,0]
    #Qext=data_dust[:,1]

    #return wavelengths_q, Qext
    wavelengths_q, Qext = read_dust_qext(specie)
       
    rho, molmass = read_dust_prop(specie)
    rho=float(rho)
    molmass=float(molmass)
    # interpolating lab data
    interpfunc_Qext=interpolate.interp1d(wavelengths_q,Qext)
    Qext_int=interpfunc_Qext(wavelength_int_q)

    # How many dust particles do we have along the line of sight for a certain Av value?
    #Nd = calculate_Nd(Av)
    Nd=calculate_Nd(Av, rho, molmass, fraction, depletion_si)
       
    # Now we need the lab spectrum part of the intensity
    tau_lab[teller,:]=Qext_int*math.pi*(0.1*1e-4)**2*Nd*fraction[teller]

    teller=teller+1
    
tau_lab=np.sum(tau_lab, axis=0)

    # interpolating lab data
    
# now we need to bring back the original stellar spectrum      
tau=tau_int+tau_lab
intensity=model_int2*np.exp(-1*(tau))
wavelength_int_where=np.where(wavelength_int<2) 
wavelength_int_again=wavelength_int[wavelength_int_where]
spectrum_again=spectrum3_int[wavelength_int_where]

spectrum_total=np.concatenate((spectrum_again,intensity))
wavelength_total=np.concatenate((wavelength_int_again,wavelength_int_q))

spectrum_total_func=interpolate.interp1d(wavelength_total,spectrum_total,fill_value="extrapolate")
spectrum_total_int = spectrum_total_func(wavelength_int)
# Adding noise to the spectrum
noise = spec_noise(spectrum_total_int,noise_para=None,noise_type='Poisson') #Noise = sqrt(Signal)
#noise = spec_noise(spectrum_total_int,noise_para=[0.0,0.10],noise_type='Gaussian') #Noise = 10% Gaussian uncertainty
    
spectrum_total_int_noisy = spectrum_total_int + noise
    

# we now also need to consider binning and the 
# rebinning the spectrum
#JWST MIRSPEC resolutions by filter: 
#https://jwst-docs.stsci.edu/mid-infrared-instrument/miri-observing-modes/miri-medium-resolution-spectroscopy#MIRIMediumResolutionSpectroscopy-wavelengthMRSwavelengthcoverage
 
##lam_min, lam_max, avg. spec. resln.
#spectral_resolution = {'nirspec_prism' :[0.6,5.0,100],
                       #'nirspec_lores' :[0.6,5.0,1000],
                       #'nirspec_hires' :[0.6,5.0,2700],
                       #'mirispec_ch1_s':[4.88,5.75,3515],
                       #'mirispec_ch1_m':[5.63,6.63,3470],
                       #'mirispec_ch1_l':[6.41,7.52,3355],
                       #'mirispec_ch2_s':[7.48,8.76,3050],
                       #'mirispec_ch2_m':[8.71 ,10.23,2960],
                       #'mirispec_ch2_l':[10.02,11.75,3080],
                       #'mirispec_ch3_s':[11.52,13.49,2705],
                       #'mirispec_ch3_m':[13.36,15.65,2215],
                       #'mirispec_ch3_l':[15.43,18.08,2385],
                       #'mirispec_ch4_s':[17.65,20.94,1695],
                       #'mirispec_ch4_m':[20.41,24.22,1725],
                       #'mirispec_ch4_l':[23.88,28.34,1495]}
# mock JWST spectral data    
    
spectral_resolution = {'nirspec_hires' :[0.6,4.87,270],
                       'mirispec_ch1_s':[4.88,5.63,352],
                       'mirispec_ch1_m':[5.63,6.40,347],
                       'mirispec_ch1_l':[6.41,7.47,336],
                       'mirispec_ch2_s':[7.48,8.76,305],
                       'mirispec_ch2_m':[8.77 ,10.23,296],
                       'mirispec_ch2_l':[10.25,11.75,308],
                       'mirispec_ch3_s':[11.76,13.49,271],
                       'mirispec_ch3_m':[13.50,15.65,222],
                       'mirispec_ch3_l':[15.66,18.08,239],
                       'mirispec_ch4_s':[18.09,20.94,170],
                       'mirispec_ch4_m':[20.95,24.22,173],
                       'mirispec_ch4_l':[24.23,28.34,150]}
     
    
#empty arrays for the final spectrum
jwst_spec_wave = np.asarray([])
jwst_spec_flux = np.asarray([])
 
for key in spectral_resolution:
    channel = spectral_resolution[key]
     
    new_waves = np.arange(channel[0],channel[1],
                           (channel[1]-0.5*(channel[1]-channel[0]))/channel[2])
    new_waves = np.append(new_waves,channel[1])
    new_fluxs = spectres(new_waves,wavelength_int,spectrum_total_int_noisy)
      
    jwst_spec_wave = np.append(jwst_spec_wave,new_waves)
    jwst_spec_flux = np.append(jwst_spec_flux,new_fluxs)
        
jwst_spex_error = spec_noise(jwst_spec_flux,noise_para=None,noise_type='Poisson') #Noise = sqrt(Signal)
jwst_spex_error2 = np.sqrt(jwst_spex_error**2)
    
    #plot figures
    
    # writing data to ascii text file, which 
    
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel(r'Wavelength [$\mu$m]',fontsize=16)
ax.set_ylabel(r'Flux density [mJy]',fontsize=16)
ax.set_ylim(1e-1,1e3)
ax.set_xlim(1,40)
ax.loglog(wavelength_int,spectrum_total_int,color='blue')
ax.loglog(wav_um,flux_mjy,color='red')
ax.text(20,5e2,r'$A_{V}$ = '+str(Av),fontsize=16)

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel(r'Wavelength [$\mu$m]',fontsize=16)
ax.set_ylabel(r'Relative intenstity (compared to photosphere)',fontsize=16)
ax.set_ylim(1e-1,1.1e0)
ax.set_xlim(1,40) 
ax.loglog(wavelength_int,spectrum_total_int_noisy/spectrum3_int,color='green')
ax.loglog(wavelength_int,spectrum_total_int/spectrum3_int,color='blue')
ax.loglog(wavelength_int,model_int/model_int,color='red')
ax.text(20,0.7,r'$A_{V}$ = '+str(Av),fontsize=16)
plt.show()
plt.close()
    
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel(r'Wavelength [$\mu$m]',fontsize=16)
ax.set_ylabel(r'Flux density [mJy]',fontsize=16)
ax.set_ylim(1e-1,1e3)
ax.set_xlim(1,40)
ax.set_xscale('log')
ax.set_yscale('log')
ax.loglog(wavelength_int,spectrum3_int,color='black')
cval = np.arange(0,1,1/len(spectral_resolution))
nk = 0
for keys in spectral_resolution:
    #cval = (nk+1)/len(spectral_resolution)
    lam_lo = spectral_resolution[keys][0]
    lam_hi = spectral_resolution[keys][1]
    filt = np.where((jwst_spec_wave >= lam_lo) & (jwst_spec_wave <= lam_hi))
    ax.scatter(jwst_spec_wave[filt],jwst_spec_flux[filt], cmap='rainbow',linestyle='-',marker='.',linewidth=0.5)
    nk += 1
ax.text(20,5e2,r'$A_{V}$ = '+str(Av),fontsize=16)

plt.show()
   
fig = plt.figure(figsize=(8,6))    
ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel(r'Wavelength [$\mu$m]',fontsize=16)
ax.set_ylabel(r'Flux density [mJy]',fontsize=16)
ax.set_ylim(1e-1,1e3)
ax.set_xlim(1,40)
ax.set_xscale('log')
ax.set_yscale('log')
ax.loglog(wavelength_int,spectrum3_int,color='black')    
ax.errorbar(jwst_spec_wave, jwst_spec_flux, yerr=jwst_spex_error2,fmt='-o')
plt.show()
plt.close()
 
zeroes=np.zeros(len(jwst_spec_flux))

np.savetxt(filename_array[i], np.c_[jwst_spec_wave, jwst_spec_flux, jwst_spex_error2, zeroes], delimiter=' ')

