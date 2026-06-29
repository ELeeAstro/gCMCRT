# Import needed packages
import numpy as np
import matplotlib.pylab as plt
from scipy.stats import binned_statistic
from pathlib import Path

# Nice fints (not required)
plt.rc('font', family='sans-serif')
plt.rc('font', serif='Helvetica Neue')
plt.rc('text', usetex='false')

def read_emission_file(path):
  head = np.atleast_1d(np.loadtxt(path,max_rows=1))
  data = np.atleast_2d(np.loadtxt(path,skiprows=1))

  wl = data[:,0]
  frac = data[:,1]
  if data.shape[1] >= 4:
    frac_occ = data[:,2]
    Ltot = data[:,3]
  else:
    frac_occ = frac
    Ltot = data[:,2]

  phase = float(head[5]) if len(head) > 5 else float(head[3])

  return {
    'wl': wl,
    'frac': frac,
    'frac_occ': frac_occ,
    'Ltot': Ltot,
    'Rp': float(head[1]),
    'Rp2': float(head[2]),
    'viewphi': float(head[3]) if len(head) > 3 else np.nan,
    'viewthet': float(head[4]) if len(head) > 4 else np.nan,
    'phase': phase,
    'occ_state': int(head[6]) if len(head) > 6 else 0,
    'xstar': float(head[7]) if len(head) > 7 else np.nan,
    'ystar': float(head[8]) if len(head) > 8 else np.nan,
    'zstar': float(head[9]) if len(head) > 9 else np.nan,
    'separation': float(head[10]) if len(head) > 10 else np.nan,
  }


def wavelength_edges_for(wl):
  edges = np.loadtxt('wavelengths_R100_UV_edge.txt',skiprows=1)[:,1]
  centers = np.loadtxt('wavelengths.wl',skiprows=1)[:,1]
  start = int(np.argmin(np.abs(centers - wl[0])))
  return edges[start:start+len(wl)+1]


em_files = sorted(Path('.').glob('Em_*.txt'))
if not em_files:
  raise FileNotFoundError('No Em_*.txt files found')

first = read_emission_file(em_files[0])
n_ph = len(em_files)
nwl = len(first['wl'])

# Calculate planetary flux from Em_* files
Fp = np.zeros((n_ph,nwl))
Fp_occ = np.zeros((n_ph,nwl))
viewphi = np.zeros(n_ph)
phase = np.zeros(n_ph)
occ_state = np.zeros(n_ph, dtype=int)

# Cycle through each phase, reading the emission data
for n, fname in enumerate(em_files):
  print(fname)
  em = read_emission_file(fname)
  wl = em['wl']
  Rp = em['Rp']
  Rp2 = em['Rp2']
  viewphi[n] = em['viewphi']
  phase[n] = em['phase']
  occ_state[n] = em['occ_state']
  # Calculate the flux of the planet toward phase n
  # Remember, flux is defined as the energy passing through a surface, so here the surface can be 
  # scaled between Rp and Rp2 at will, but we typically use Rp2 or Rp
  Fp[n,:] = (em['frac'][:] * em['Ltot'][:]) / (Rp2**2) # erg s-1 cm-2 cm-1
  Fp_occ[n,:] = (em['frac_occ'][:] * em['Ltot'][:]) / (Rp2**2) # erg s-1 cm-2 cm-1


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
wl_ed = wavelength_edges_for(wl)

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
FpFs_occ = np.zeros((n_ph,nwl))
for n in range(n_ph):
    #FpFs[n,:] = (Rp/Rs)**2 * (Fp[n,:]/Fs_binned)
    FpFs[n,:] = 0.103**2 * (Fp[n,:]/Fs_binned)
    FpFs_occ[n,:] = 0.103**2 * (Fp_occ[n,:]/Fs_binned)

# Planetary flux plot
fig = plt.figure()

for n in range(n_ph):
  plt.plot(wl,Fp[n,:],label=f'Phase {phase[n]:.4f} (occ {occ_state[n]})')

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

for n in range(n_ph):
  plt.plot(wl,BT[n,:],label=f'Phase {phase[n]:.4f} (occ {occ_state[n]})')

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

for n in range(n_ph):
  plt.plot(wl,FpFs[n,:],label=f'Phase {phase[n]:.4f} (occ {occ_state[n]})')

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
