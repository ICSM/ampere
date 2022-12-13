import math
import numpy as np

# rho_d and n_0 for specific dust species, rho_d from literature, n_0 from solution

#rho_d:
#am. oli: 3.71
#cr. oli: 3.33
#calcite: 2.71
#cr. ice: 1.00
#diopside: 3.4
#cr. pyr: 2.80
#dolomite: 2.84
#iron: 7.874

rhod = np.array([3.71, 3.33, 2.71, 1.00, 3.4, 2.80, 2.84, 7.874])
names = np.array(['am. oliv.', 'forst.', 'calcite', 'ice', 'diopside', 'c-enst.', 'dolomite', 'iron'])
n0 = np.array([(10**(-1.47491), 10**(-3.12641), 10**(-4.95898), 10**(-4.78796), 10**(-4.38619), 10**(-4.69366), 10**(-5.03824), 0), #cold
               (10**(-4.86979),10**(-5.22728), 0, 0, 0, 0, 0, 10**(-2.54621))])
n0pos = np.array([(10**(-1.47491+0.04872), 10**(-3.12641+0.10695), 10**(-4.95898+0.12598), 10**(-4.78796+0.74259), 10**(-4.38619+0.06672), 10**(-4.69366+0.23300), 10**(-5.03824+0.26593), 0), #cold
               (10**(-4.86979+0.50563),10**(-5.22728+0.26170), 0, 0, 0, 0, 0, 10**(-2.54621+0.19764))])
n0neg = np.array([(10**(-1.47491-0.05432), 10**(-3.12641-0.10895), 10**(-4.95898-0.11465), 10**(-4.78796-0.70013), 10**(-4.38619-0.08870), 10**(-4.69366-0.38323), 10**(-5.03824-0.33590), 0), #cold
               (10**(-4.86979-0.54834),10**(-5.22728-0.27771), 0, 0, 0, 0, 0, 10**(-2.54621-0.16419))])

#warm

# Tout and Tin from solution
Tout = np.array([27.83584, 108.80126]) 
Tin = np.array([61.06154, 122.91945])

#other constants/values
p = 0.5 #slope density distribution
q = 0.5 #slope temperature distribution
a = 0.1e-4 #grain size (radius) in cm 
msun = 1.9885e33 #in g

# calculate r0

#r0 = 1e15 #placeholder value, to be calculated
r0 = np.zeros((2))
r0cubed = 4.7e-2*msun / (math.pi/(3-p) * ((30/60)**(-4)-1) * (4/3) * math.pi * a**3 * 3.71 * 3.9e-2) #Kemper et al, 2002, using cold amorphous olivine
# m = 4.7e-2 Msun, Tin = 60, Tout = 30 K, n0 = 3.9e-2, rhod = 3.71
r02002 = r0cubed**(1/3)

# calculate our r0 from equation (4) from Kemper et al. 2002, A&A
r0[0] = r02002 * (Tin[0]/60)**(-1/q)

# same for warm component
r0cubed = 6.1e-6*msun / (math.pi/(3-p) * ((100/118)**(-4)-1) * (4/3) * math.pi * a**3 * 3.71 * 1.2e-4) #Kemper et al, 2002, using cold amorphous olivine
# m = 6.1e-6 Msun, Tin = 118, Tout = 100 K, n0 = 1.2e-4, rhod = 3.71
r02002 = r0cubed**(1/3)

# calculate our r0 from equation (4) from Kemper et al. 2002, A&A
r0[1] = r02002 * (Tin[1]/118)**(-1/q)


print('r0 (cold): ',r0[0],' cm')
print('r0 (warm): ',r0[1],' cm')

# Use Eq(5) from Kemper et al. 2002, A&A, integrated over the range Tin - Tout
# to calculate dust mass in c.g.s.

nspecies = names.__len__()
mdust = np.zeros((2,nspecies))
mdustup = np.zeros((2,nspecies))
mdustlow = np.zeros((2,nspecies))

              
for i in range(nspecies):
    mdust[0,i] = ((math.pi * r0[0]**3)/(3-p) * ((Tout[0]/Tin[0])**(-4)-1) * (4/3) * math.pi * a**3 * rhod[i] * n0[0,i])/msun # in solar masses
    mdust[1,i] = ((math.pi * r0[1]**3)/(3-p) * ((Tout[1]/Tin[1])**(-4)-1) * (4/3) * math.pi * a**3 * rhod[i] * n0[1,i])/msun # in solar masses
    mdustup[0,i] = ((math.pi * r0[0]**3)/(3-p) * ((Tout[0]/Tin[0])**(-4)-1) * (4/3) * math.pi * a**3 * rhod[i] * n0pos[0,i])/msun # in solar masses
    mdustup[1,i] = ((math.pi * r0[1]**3)/(3-p) * ((Tout[1]/Tin[1])**(-4)-1) * (4/3) * math.pi * a**3 * rhod[i] * n0pos[1,i])/msun # in solar masses
    mdustlow[0,i] = ((math.pi * r0[0]**3)/(3-p) * ((Tout[0]/Tin[0])**(-4)-1) * (4/3) * math.pi * a**3 * rhod[i] * n0neg[0,i])/msun # in solar masses
    mdustlow[1,i] = ((math.pi * r0[1]**3)/(3-p) * ((Tout[1]/Tin[1])**(-4)-1) * (4/3) * math.pi * a**3 * rhod[i] * n0neg[1,i])/msun # in solar masses
    print(names[i],": ",mdust[0,i]," Msun in cold dust; ",mdust[1,i]," Msun in warm dust")
    print("cold dust mass solutions range from ", mdustlow[0,i], " to ", mdustup[0,i])
    print("warm dust mass solutions range from ", mdustlow[1,i], " to ", mdustup[1,i])


#(3.14*(1e15)**3)/2.5 * ((27.83/61.06)**(-4) -1) * (4/3) * 3.14 * (0.1e-4)**3 * rho_d * n0 in c.g.s units




#another option is to calibrate the cold dust mass to what is derived in our 2002 paper for the total cold dust mass: 0.050 Msun

