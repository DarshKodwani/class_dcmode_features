import matplotlib.pyplot as plt
import numpy as np
from classy import Class

params = { 
    'output': 'tCl pCl tpCl',
    'lensing': 'no',
    'A_s' : 2.3e-9,
    'h' : 0.6711,
    'omega_b' : 0.022068,
    'omega_cdm' : 0.12029}

lmax=2500
l = np.array(range(2,lmax+1))
factor = l*(l+1)/(2*np.pi)

# Defining detector noise                                                                                                                                                                               
#Using the defintiion of detector noise used in "Forecasting isocurvature models with CMB lensing information" by Santos et al (2012). See Eq A9 in this paper.                                            # Defining the noise paramaters - all quantities taken from the paper given above Table VII - we only take the 143 GHz channel                                                                            
arcmin_rad = 3437.75 
thetaarcmin = 7.1
thetaarcmin100 = 9.5
thetaarcmin217 = 5.0
thetaarcmin353 = 5.0
thetarad = thetaarcmin/arcmin_rad
thetarad100 = thetaarcmin100/arcmin_rad
thetarad217 = thetaarcmin217/arcmin_rad
thetarad353 = thetaarcmin353/arcmin_rad

sigmaT = 6.0016*(10**(-6))
sigmaT100 = 6.82*(10**(-6))
sigmaT217 = 13.0944*(10**(-6))
sigmaT353 = 40.1016*(10**(-6))

#sigmaT = 6.0016
#sigmaT100 = 6.82
#sigmaT217 = 13.0944
#sigmaT353 = 40.1016


sigmaP = 11.4576*(10**(-6))
sigmaP100 = 10.9120*(10**(-6))
sigmaP217 = 26.7644*(10**(-6))
sigmaP353 = 81.2944*(10**(-6))

#sigmaP = 11.4576
#sigmaP100 = 10.9120
#sigmaP217 = 26.7644
#sigmaP353 = 81.2944

def dnoise(l):                                                                                                                                                                                             
    return ((thetaarcmin*sigmaT)**2)*np.exp(l*(l+1)*thetarad**2/(8*np.log(2)))                                                                                                                             
#                                                                                                                                                                                                           
def dnoiseP(l):                                                                                                                                                                                            
    return ((thetaarcmin*sigmaP)**2)*np.exp(l*(l+1)*thetarad**2/(8*np.log(2)))                                                                                                                             

def dnoise_full(l):
    return ( ((thetaarcmin*sigmaT)**(-2))*np.exp(-l*(l+1)*(thetarad**2)/(8*np.log(2)))
    + ((thetaarcmin100*sigmaT100)**(-2))*np.exp(-l*(l+1)*(thetarad100**2)/(8*np.log(2)))
    + ((thetaarcmin217*sigmaT217)**(-2))*np.exp(-l*(l+1)*(thetarad217**2)/(8*np.log(2)))
    + ((thetaarcmin353*sigmaT353)**(-2))*np.exp(-l*(l+1)*(thetarad353**2)/(8*np.log(2))))**(-1)

def dnoiseP_full(l):
    return (  ((thetaarcmin*sigmaP)**(-2))*np.exp(-l*(l+1)*(thetarad**2)/(8*np.log(2)))
    + ((thetaarcmin100*sigmaP)**(-2))*np.exp(-l*(l+1)*(thetarad100**2)/(8*np.log(2)))
    + ((thetaarcmin217*sigmaP217)**(-2))*np.exp(-l*(l+1)*(thetarad217**2)/(8*np.log(2)))
    + ((thetaarcmin353*sigmaP353)**(-2))*np.exp(-l*(l+1)*(thetarad353**2)/(8*np.log(2))))**(-1)

# Theory Cls

cosmo = Class()
cosmo.set(params)
cosmo.compute()
fullcl = cosmo.raw_cl(lmax)

#plt.loglog(fullcl['ell'][2:], dnoise(fullcl['ell'][2:]), label="Noise")
#plt.loglog(fullcl['ell'][2:], dnoise_full(fullcl['ell'][2:]), label="Noise full")
plt.loglog(fullcl['ell'][2:], factor * fullcl['tt'][2:], label="Temp")
plt.legend()
plt.show()
#quit()
TT=factor*fullcl['tt'][2:] + dnoise(fullcl['ell'][2:])
EE=factor*fullcl['ee'][2:] + dnoiseP(fullcl['ell'][2:])
TE=factor*fullcl['te'][2:]
TT2 = TT**2
EE2 = EE**2
TE2 = TE**2

#### Fisher matrix ####

numns = 3
ns_arr = np.linspace(0.8, 1.2, numns)
dns = 0.05
numAs = 4
As_arr = np.linspace(2.1e-9, 2.4e-9, numAs)
dAs = 10e-10

numks = 3
ks_arr = np.linspace(10**(-4), 10**(0), numks)
dks = 10**(-5)


fish_sum = []
full_fish = np.zeros((numAs,numAs))

for i in np.arange(numAs) : 
    cosmop = Class(); cosmop.set(params) ; cosmop.set({'A_s' : As_arr[i]*(1 + dAs)}) ; cosmop.compute() ; cosmop.empty()
    cosmom = Class(); cosmom.set(params) ; cosmom.set({'A_s' : As_arr[i]*(1 - dAs)}) ; cosmom.compute() ; cosmom.empty()
    dTTi = (np.array(factor*cosmop.raw_cl(lmax)['tt'][2:]) - np.array(factor*cosmom.raw_cl(lmax)['tt'][2:]))/(2*dAs)
    dEEi = (np.array(factor*cosmop.raw_cl(lmax)['ee'][2:]) - np.array(factor*cosmom.raw_cl(lmax)['ee'][2:]))/(2*dAs)
    dTEi = (np.array(factor*cosmop.raw_cl(lmax)['te'][2:]) - np.array(factor*cosmom.raw_cl(lmax)['te'][2:]))/(2*dAs)
    for j in range(i, numAs) : 
        cosmop = Class(); cosmop.set(params) ; cosmop.set({'A_s' : As_arr[j]*(1 + dAs)}) ; cosmop.compute() ; cosmop.empty()
        cosmom = Class(); cosmom.set(params) ; cosmom.set({'A_s' : As_arr[j]*(1 - dAs)}) ; cosmom.compute() ; cosmom.empty()
        dTTj = (np.array(factor*cosmop.raw_cl(lmax)['tt'][2:]) - np.array(factor*cosmom.raw_cl(lmax)['tt'][2:]))/(2*dAs)
        dEEj = (np.array(factor*cosmop.raw_cl(lmax)['ee'][2:]) - np.array(factor*cosmom.raw_cl(lmax)['ee'][2:]))/(2*dAs)
        dTEj = (np.array(factor*cosmop.raw_cl(lmax)['te'][2:]) - np.array(factor*cosmom.raw_cl(lmax)['te'][2:]))/(2*dAs)
        for l in range(2,lmax-1):
            fsky=0.5
            fish_factor = fsky*((2*l+1)/2.)*( dTTi[l]*EE2[l]*dTTj[l]/( (TE2[l] - EE[l]*TT[l])**2) 
                                              + 2*dTTi[l]*TE2[l]*dTEj[l]/( (TE2[l] - EE[l]*TT[l])**2) 
                                              - 4*dTTi[l]*EE[l]*TE[l]*dEEj[l]/( (TE2[l] - EE[l]*TT[l])**2)
                                              + dTEi[l]*TT2[l]*dTEj[l]/( (TE2[l] - EE[l]*TT[l])**2) 
                                              - 4*dTEi[l]*TE[l]*TT[l]*dEEj[l]/( (TE2[l] - EE[l]*TT[l])**2)
                                              + 2*dEEi[l]*(TE2[l]+EE[l]*TT[l])*EE[l]/( (TE2[l] - EE[l]*TT[l])) ) 
            fish_sum.append(fish_factor)
        full_fish[i,j] = np.sum(fish_sum)
        if i!=j : 
            full_fish[j,i] = full_fish[i,j]

#print full_fish
#plt.imshow(full_fish)
#plt.show()

plt.figure()
CS = plt.contour(ks_arr,ks_arr,full_fish,100)
cbar = plat.colorbar(CS)
plt.title('Fisher information $I(k_1, k_2)$')
plt.yscale('log')
plt.xscale('log')
plt.show()
