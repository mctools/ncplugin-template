import NCrystal as NC
import numpy as np
import matplotlib.pyplot as plt
import pickle

with open("data/muvals", "rb") as f: 
    zip_ = pickle.load(f)
f.close()

with open("data/muvals_files", "rb") as f: 
    zip_files = pickle.load(f)
f.close()

nbins = 200
nsample = 1e7

muvals=[]
muvals2=[]
muvals_files=[]
muvals2_files=[]
for mu, mu2 in zip_:
    muvals.append(mu)
    muvals2.append(mu2)

for mu_files, mu2_files in zip_files:
    muvals_files.append(mu_files)
    muvals2_files.append(mu2_files)

plt.hist(muvals, bins=nbins,alpha=0.5, label=f'GP',range=(0.8,1),color="k", fill=False,histtype="step")
plt.hist(muvals2, bins=nbins,alpha=0.5, label=f'Simple',range=(0.8,1),color="r", fill=False,histtype="step")

plt.hist(muvals_files, bins=nbins,alpha=0.5, label=f'teshi_xs',range=(0.8,1),color="g", fill=False,histtype="step")
plt.hist(muvals2_files, bins=nbins,alpha=0.5, label=f'Ersez_data_corr',range=(0.8,1),color="b", fill=False,histtype="step")

plt.title(r'$\mu=\cos(\theta)$' + f' distribution in {float(nsample):g} scattering events'+ r' $\lambda=1.98\AA$')
plt.legend()
plt.xlabel(r'$\mu=\cos(\theta)$')
plt.yscale('log')
plt.show()