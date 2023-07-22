# Import needed packages
import numpy as np
import matplotlib.pylab as plt
from scipy.stats import binned_statistic

# Nice fints (not required)
plt.rc('font', family='sans-serif')
plt.rc('font', serif='Helvetica Neue')
plt.rc('text', usetex='false')

# Number of phases and wavelength points
n_ph = 2
nwl = 503

# Calculate planetary flux from Em_* files
Fp = np.zeros((n_ph,nwl))

# Cycle through each phase, reading the emission data
for n in range(n_ph):
  fname = 'Em_'+str(n+1).zfill(3)+'.txt'
  print(fname)
  head = np.loadtxt(fname,max_rows=1)
  # Phase (longitude degree)
  ph = int(head[3])
  # Radius of planet at surface
  Rp = float(head[1])
  # Radius of planet + atmosphere top
  Rp2 = float(head[2])
  data = np.loadtxt(fname,skiprows=1)
  # Read wavelength
  wl = data[:,0]
  # Read the escaped energy fraction
  frac = data[:,1]
  # Read total energy
  Ltot = data[:,2]
  # Calculate the flux of the planet toward phase n
  # Remember, flux is defined as the energy passing through a surface, so here the surface can be 
  # scaled between Rp and Rp2 at will, but we typically use Rp2 or Rp
  Fp[n,:] = (frac[:] * Ltot[:]) / (Rp2**2) # erg s-1 cm-2 cm-1


# Calculate the brightness temperature from inverting the Planck function
c0 = 2.99792458e10
h = 6.626176e-27
kb = 1.3080662e-16

wl_cm = wl * 1e-4

# Brightness temperature calculation
BT = np.zeros((n_ph,nwl))
for n in range(n_ph):
  BT[n,:] = (h*c0)/(kb*wl_cm) / np.log(1.0 +  ((2.0*h*c0**2)/(Fp[n,:]/np.pi * wl_cm**5)))

# Calculate the Fp/Fs
data = np.loadtxt('wavelengths_R100_UV_edge.txt',skiprows=1)
wl_ed = data[:,1]

# Stellar properties
Rsun = 6.957e10
Rs = 1.444 * Rsun

# Read in stellar fluxes and convert to same units
fname = 'WASP-33_stellar_spectrum.txt'
data = np.loadtxt(fname,skiprows=1)
wl_s = data[:,0]
Fs = data[:,1] * 1e4  #in erg/s/cm2/um to erg/s/cm2/cm

# Find the averaged stellar flux in each wavelength bin
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

data = np.loadtxt('obs_data.txt')
wl_obs = data[:,0]
FpFs_obs = data[:,1]/100.0
FpFs_obs_err = data[:,2]/100.0

fig = plt.figure()

plt.plot(wl,FpFs[0,:],label='Dayside')
plt.plot(wl,FpFs[1,:],label='Nightside')

plt.scatter(wl_obs,FpFs_obs,c='darkcyan',marker='o',s=12.0,label='Obs. data',zorder=1)
plt.errorbar(wl_obs,FpFs_obs,yerr=FpFs_obs_err,fmt='none',c='darkcyan',zorder=1)

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
