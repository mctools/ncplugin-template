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

with open("data/muvals_dist", "rb") as f: 
    zip_dist = pickle.load(f)
f.close()

with open("data/muvals_mono_R", "rb") as f: 
    zip_mono_R = pickle.load(f)
f.close()

nbins = 100
nsample = 1e6

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

muvals_dist=[]
muvals2_dist=[]
muvals_mono_R=[]
muvals2_mono_R=[]
for mu, mu2 in zip_dist:
    muvals_dist.append(mu)
    muvals2_dist.append(mu2)

for mu_mono_R, mu2_mono_R in zip_mono_R:
    muvals_mono_R.append(mu_mono_R)
    muvals2_mono_R.append(mu2_mono_R)   

fig, ax = plt.subplots(1,1, figsize=(7, 5), dpi=200) 

plt.hist(muvals, bins=nbins,alpha=0.8, label=f'PPF',range=(0.99999,1),color="k", fill=False,histtype="step",density=True)
print(f"Average mu in PPF: {sum(muvals)/len(muvals)}")
plt.hist(muvals2, bins=nbins,alpha=0.8, label=f'GPF',range=(0.99999,1),color="tab:red", fill=False,histtype="step",density=True)
print(f"Average mu in GPF: {sum(muvals2)/len(muvals2)}")

plt.hist(muvals_files, bins=nbins,alpha=0.8, label=f'Teshi',range=(0.99999,1),color="tab:green", fill=False,histtype="step",density=True)
print(f"Average mu in Teshi: {sum(muvals_files)/len(muvals_files)}")
plt.hist(muvals2_files, bins=nbins,alpha=0.8, label=f'Grammer Sample 1',range=(0.99999,1),color="tab:blue", fill=False,histtype="step",density=True)
print(f"Average mu in Grammer Sample 1: {sum(muvals2_files)/len(muvals2_files)}")

plt.hist(muvals_mono_R, bins=nbins,alpha=0.8, label=f'R=50AA',range=(0.99999,1),color="tab:cyan", fill=False,histtype="step",density=True)
print(f"Average mu in HS 50 AA: {sum(muvals_mono_R)/len(muvals_mono_R)}")
plt.hist(muvals2_mono_R, bins=nbins,alpha=0.8, label=f'R=25AA',range=(0.99999,1),color="tab:purple", fill=False,histtype="step",density=True)
print(f"Average mu in HS 25 AA: {sum(muvals2_mono_R)/len(muvals2_mono_R)}")

plt.hist(muvals_dist, bins=nbins,alpha=0.8, label=f'Teshi R dist',range=(0.99999,1),color="tab:orange", fill=False,histtype="step",density=True)
print(f"Average mu in Teshi R dist: {sum(muvals_dist)/len(muvals_dist)}")
plt.hist(muvals2_dist, bins=nbins,alpha=0.8, label=f'Alksenskii R dist',range=(0.99999,1),color="tab:brown", fill=False,histtype="step",density=True)
print(f"Average mu in Alksenskii R dist: {sum(muvals2_dist)/len(muvals2_dist)}")

plt.title(r'$\mu=\cos(\theta)$' + f' PDF for'+ r' $\lambda=20\AA$')
plt.legend(loc="upper left")
plt.xlabel(r'$\mu=\cos(\theta)$')
#plt.yscale('log')
plt.tight_layout()
plt.show()