import numpy as np
import matplotlib.pylab as plt
from scipy.stats import binned_statistic

plt.rc('font', family='sans-serif')
plt.rc('font', serif='Helvetica Neue')
plt.rc('text', usetex='false')

n_ph = 2
nwl = 503

# Calculate planetary flux from Em_* files
Fp = np.zeros((n_ph,nwl))

for n in range(n_ph):
  fname = 'Em_'+str(n+1).zfill(3)+'.txt'
  print(fname)
  head = np.loadtxt(fname,max_rows=1)
  ph = int(head[3])
  Rp = float(head[1])
  Rp2 = float(head[2])
  data = np.loadtxt(fname,skiprows=1)
  wl = data[:,0]
  frac = data[:,1]
  Ltot = data[:,2]
  Fp[n,:] = (frac[:] * Ltot[:]) / (Rp2**2) # erg s-1 cm-2 cm-1


# Calculate the brightness temperature from inverting the Planck function
c0 = 2.99792458e10
h = 6.626176e-27
kb = 1.3080662e-16

wl_cm = wl * 1e-4

BT = np.zeros((n_ph,nwl))
for n in range(n_ph):
  BT[n,:] = (h*c0)/(kb*wl_cm) / np.log(1.0 +  ((2.0*h*c0**2)/(Fp[n,:]/np.pi * wl_cm**5)))

# Calculate the Fp/Fs

data = np.loadtxt('wavelengths_R100_UV_edge.txt',skiprows=1)
wl_ed = data[:,1]

Rsun = 6.957e10
Rs = 1.444 * Rsun

fname = 'WASP-33_stellar_spectrum.txt'
data = np.loadtxt(fname,skiprows=1)
wl_s = data[:,0]
Fs = data[:,1] * 1e4  #in erg/s/cm2/um to erg/s/cm2/cm

Fs_binned, bin_edges, binnumber = binned_statistic(wl_s,Fs,bins=wl_ed)


# Now you have the option of using the observed Rp/Rs or model Rp/Rs scaling factor
# Here we use the observed value since it fits better
FpFs = np.zeros((n_ph,nwl))
for n in range(n_ph):
    #FpFs[n,:] = (Rp/Rs)**2 * (Fp[n,:]/Fs_binned)
    FpFs[n,:] = 0.103**2 * (Fp[n,:]/Fs_binned)

# Planetary flux plot
fig = plt.figure()

plt.plot(wl,Fp[0,:],label='Dayside')
plt.plot(wl,Fp[1,:],label='Nightside')

plt.legend()

plt.xscale('log')
plt.yscale('log')

plt.xlabel(r'$\lambda$ [$\mu$m]',fontsize=14)
plt.ylabel(r'F$_{\rm p}$ [erg s$^{-1}$ cm$^{-2}$ cm$^{-1}$]',fontsize=14)

plt.xlim(0.2,31.0)

xticks = [0.2,0.5,1,2,3,4,5,6,8,10,20,30]
xticks_lab = ['0.2','0.5','1','2','3','4','5','6','8','10','20','30']

plt.xticks(xticks,xticks_lab)

plt.tick_params(axis='both',which='major',labelsize=12)

plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)

#plt.savefig('Fp.png',dpi=300,bbox_inches='tight')

# Brightness temperature plots
fig = plt.figure()

plt.plot(wl,BT[0,:],label='Dayside')
plt.plot(wl,BT[1,:],label='Nightside')

plt.legend()

plt.xscale('log')

plt.xlabel(r'$\lambda$ [$\mu$m]',fontsize=14)
plt.ylabel(r'T$_{\rm b}$ [K]',fontsize=14)

plt.xlim(0.2,31.0)

xticks = [0.2,0.5,1,2,3,4,5,6,8,10,20,30]
xticks_lab = ['0.2','0.5','1','2','3','4','5','6','8','10','20','30']

plt.xticks(xticks,xticks_lab)

plt.tick_params(axis='both',which='major',labelsize=12)

plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)

#plt.savefig('BT.png',dpi=300,bbox_inches='tight')

# Fp/Fs plots

# Observational data
S36 = 3506.0/1e6 
S36_err = 173.0/1e6

S45 = 4250.0/1e6
S45_err = 160.0/1e6

fig = plt.figure()

plt.plot(wl,FpFs[0,:],label='Dayside')
plt.plot(wl,FpFs[1,:],label='Nightside')

plt.scatter(3.6,S36,c='darkcyan',marker='o',s=12.0,label='Zhang et al. (2013)',zorder=1)
plt.errorbar(3.6,S36,yerr=S36_err,fmt='none',c='darkcyan',zorder=1)

plt.scatter(4.5,S45,c='darkcyan',marker='o',s=12.0,zorder=1)
plt.errorbar(4.5,S45,yerr=S45_err,fmt='none',c='darkcyan',zorder=1)

plt.legend()

plt.xscale('log')

plt.xlabel(r'$\lambda$ [$\mu$m]',fontsize=14)
plt.ylabel(r'F$_{\rm p}$/F$_{\star}$',fontsize=14)

plt.xlim(0.2,31.0)

xticks = [0.2,0.5,1,2,3,4,5,6,8,10,20,30]
xticks_lab = ['0.2','0.5','1','2','3','4','5','6','8','10','20','30']

plt.xticks(xticks,xticks_lab)

plt.tick_params(axis='both',which='major',labelsize=12)

plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)

#plt.savefig('FpFs.png',dpi=300,bbox_inches='tight')

plt.show()
