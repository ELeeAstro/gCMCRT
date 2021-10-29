import numpy as np
import matplotlib.pylab as plt


pc = 648000.0/np.pi
au = 1.495978707e13
c_s = 2.99792458e10


data = np.loadtxt('Em.txt',skiprows=1)
wl = data[:,0]
frac = data[:,1]
Ltot = data[:,2] 

L = frac * Ltot


R_jup = 6.9911e9
R_p = 1.25 * R_jup
#R_p = 8684071391.53085
R = 8905787324.09163

Fp = L / (np.pi *  R**2) * (np.pi * R_p**2) 


fname = 'seA_15___.tsv'
data = np.loadtxt(fname)

wl2 = data[:,0]
Flx2 = data[:,1]

fname = 'seE_15___.tsv'
data = np.loadtxt(fname)

wl3 = data[:,0]
Flx3 = data[:,1]

fname = 'sep_15___.tsv'
data = np.loadtxt(fname)

wl4 = data[:,0]
Flx4 = data[:,1]


dist = 10.0
De = dist * pc * au
freq = c_s / (wl[:] * 1.0e-4)

Fp = (Fp *  (wl[:] * 1.0e-4) / freq) / ( De**2)

fig = plt.figure()

plt.plot(wl,Fp,c='darkorange',lw=1,label='gCMCRT',zorder=4)
plt.plot(wl2,Flx2,c='darkred',label=r'ATMO',zorder=2,lw=2)
plt.plot(wl3,Flx3,c='darkblue',label='Exo-REM',zorder=3,lw=2)
plt.plot(wl4,Flx4,c='darkgreen',label='petitCODE',zorder=1,lw=2)

plt.xscale('log')
#plt.yscale('log')
plt.xlim(0.3,30)
plt.ylim(0,2.0e-25)

plt.legend()

plt.xlabel('Wavelength [um]')
plt.ylabel('Fp [erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$]')


plt.show()
