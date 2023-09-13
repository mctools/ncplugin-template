import NCrystal as NC
import numpy as np
import matplotlib.pyplot as plt
import pandas
import random 

A1=132.869
b1=-1.33605
A2=0.0519763
b2=-3.97314
Q0=0.0510821
sigma0=5.551

def calcCrossSection(neutron_ekin):
    
    lamb = NC.ekin2wl(neutron_ekin) #wavelength
    k =  2*np.arccos(-1)/lamb
    total_sigma = (sigma0/(2*k*k))*(A1/(b1+2)*pow(Q0,b1+2) + A2/(b2+2)*pow(2*k,b2+2) - A2/(b2+2)*pow(Q0,b2+2))
    return total_sigma


def sampleScatteringVector(rand, neutron_ekin ): 
    lamb = NC.ekin2wl(neutron_ekin) #wavelength
    k =  2*np.arccos(-1)/lamb
    #sample a random scattering vector Q from the inverse CDF (see plugin readme)
    ratio_sigma = (sigma0/(2*k*k))/calcCrossSection(neutron_ekin) #cross section over total cross section ratio
    CDF_Q0 = (A1*pow(Q0, b1+2)/(b1+2))*ratio_sigma
    if rand < CDF_Q0 :
        Q = pow(((b1+2)*rand/A1)/ratio_sigma, 1/(b1+2))
    else:
        Q = pow((rand/ratio_sigma - (A1/(b1+2))*pow(Q0,b1+2) + (A2/(b2+2))*pow(Q0,b2+2))*(b2+2)/A2, 1/(b2+2))

    return Q

def calcCDF(Q, neutron_ekin):
    
    lamb = NC.ekin2wl(neutron_ekin) #wavelength
    k =  2*np.arccos(-1)/lamb
    ratio_sigma = (sigma0/(2*k*k))/calcCrossSection(neutron_ekin) 
    if Q < Q0:
        CDF = ratio_sigma*A1/(b1+2)*pow(Q, b1+2)
    else:
        CDF = ratio_sigma*( A1/(b1+2)*pow(Q0, b1+2) + A2/(b2+2)*pow(Q, b2+2) - A2/(b2+2)*pow(Q0, b2+2))
    return CDF

def calcP(Q, neutron_ekin):
    
    lamb = NC.ekin2wl(neutron_ekin) #wavelength
    k =  2*np.arccos(-1)/lamb
    ratio_sigma = (sigma0/(2*k*k))/calcCrossSection(neutron_ekin) 
    if Q < Q0:
        Prob = ratio_sigma*Q*A1*pow(Q, b1)
    else:
        Prob = ratio_sigma*Q*A2*pow(Q, b2)
    return Prob

#rand = random.random()
#plt.figure()
#x = np.linspace(0,1,100) # random
#y = np.linspace(0.0002,0.002,100) # energy
#xx, yy = np.meshgrid(x,y)
#z = sampleScatteringVector(xx,yy)
#plt.title('Radicand ')
#plt.ylabel('Energy (meV)')
#plt.xlabel('random number')
#plt.pcolor(xx, yy*1000, z, cmap='coolwarm', vmin=-0.001, vmax=0.001)
#plt.text(0.2,1.2, "> 0", size=30)
#plt.text(0.8,0.4, "< 0", size=30)

plt.figure()
x = np.linspace(0.001,1.4,10000) # Qs
problist = []
for qs in x:
    problist.append(calcP(qs,0.001))
plt.plot(x,problist)
histPlist = []
for i in range(0,100000):
    rand = random.random()
    histPlist.append(sampleScatteringVector(rand,0.001))
plt.hist(histPlist, bins=280, density=True)
plt.title('P(Q)')
plt.ylabel('P(Q)')
plt.xlabel('Q ($Å^{-1}$)')
plt.show()

plt.figure()
x = np.linspace(0.001,1.4,10000) # Qs
CDFlist = []
for qs in x:
    CDFlist.append(calcCDF(qs,0.001))
plt.plot(x,CDFlist)
plt.title('CDF(Q)')
plt.ylabel('p(Q)')
plt.xlabel('Q ($Å^{-1}$)')
plt.show()