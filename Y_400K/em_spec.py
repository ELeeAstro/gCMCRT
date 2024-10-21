import numpy as np
import matplotlib.pylab as plt
from scipy.stats import binned_statistic
import seaborn as sns

fname = 'Em_001.txt'
head = np.loadtxt(fname,max_rows=1)
nwl = int(head[0])
Rp_b = float(head[1])
Rp_t = float(head[2])
Fp0 = np.zeros(nwl)
data = np.loadtxt(fname,skiprows=1)
wl0 = data[:,0]
frac = data[:,1]
Ltot = data[:,2]
Fp0[:] = (frac[:] * Ltot[:]) / (Rp_t**2)

fig = plt.figure()

c = sns.color_palette('colorblind')

plt.plot(wl0,Fp0,c=c[0],label=r'Y 400 K')

plt.xscale('log')
plt.yscale('log')

plt.xlabel(r'$\lambda$ [$\mu$m]',fontsize=14)
plt.ylabel(r'F$_{\rm bd}$ [erg s$^{-1}$ cm$^{-2}$ cm$^{-1}$]',fontsize=14)


plt.legend()

#plt.xlim(1,30)
#plt.ylim(1e6,1e13)

plt.tick_params(axis='both',which='major',labelsize=12)

plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)

#plt.savefig('eq_pole_comp.png',dpi=144,bbox_inches='tight')


plt.show()
