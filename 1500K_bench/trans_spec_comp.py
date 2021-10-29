import numpy as np
import matplotlib.pylab as plt

R_sun = 6.957e10
R_jup = 6.9911e9
#R_jup = 7.1492e9

R_star = 0.78 * R_sun
R_p = 1.25 * R_jup
R_p = 8684071391.53085 

fname = 'Transmission.txt'
data = np.loadtxt(fname,skiprows=1)

wl = data[:,0]
sumf = data[:,1]

Rp = np.sqrt(R_p**2 + 2.0 * sumf)/R_jup

fname = 'stA_15___.tsv'
data = np.loadtxt(fname)

wl2 = data[:,0]
Rp2 = data[:,1]

fname = 'stp_15___.tsv'
data = np.loadtxt(fname)

wl4 = data[:,0]
Rp4 = data[:,1]


fig = plt.figure()

plt.plot(wl,Rp,label='gCMCRT',c='darkorange',lw=1,zorder=4)
plt.plot(wl2,Rp2,c='darkred',label=r'ATMO',zorder=2,lw=2)
plt.plot(wl4,Rp4,c='darkgreen',label='petitCODE',zorder=3,lw=2)

plt.legend()

plt.xscale('log')
plt.xlim(0.3,30.0)
plt.ylim(1.255,1.266)


plt.show()
