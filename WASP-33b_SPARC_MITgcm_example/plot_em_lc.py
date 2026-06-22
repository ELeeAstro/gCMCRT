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
  head = np.loadtxt(path,max_rows=1)
  data = np.atleast_2d(np.loadtxt(path,skiprows=1))

  wl = data[:,0]
  frac = data[:,1]
  if data.shape[1] >= 4:
    frac_occ = data[:,2]
    Ltot = data[:,3]
  else:
    frac_occ = frac
    Ltot = data[:,2]

  return {
    'wl': wl,
    'frac': frac,
    'frac_occ': frac_occ,
    'Ltot': Ltot,
    'Rp': float(head[1]),
    'Rp2': float(head[2]),
    'viewphi': float(head[3]) if len(head) > 3 else np.nan,
    'viewthet': float(head[4]) if len(head) > 4 else np.nan,
    'phase': float(head[5]) if len(head) > 5 else float(head[3]),
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
phase = np.zeros(n_ph)

# Cycle through each phase, reading the emission data
for n, fname in enumerate(em_files):
  print(fname)
  em = read_emission_file(fname)
  wl = em['wl']
  Rp = em['Rp']
  Rp2 = em['Rp2']
  phase[n] = em['phase']
  # Calculate the flux of the planet toward phase n
  # Remember, flux is defined as the energy passing through a surface, so here the surface can be 
  # scaled between Rp and Rp2 at will, but we typically use Rp2 or Rp
  Fp[n,:] = (em['frac'][:] * em['Ltot'][:]) / (Rp2**2) # erg s-1 cm-2 cm-1
  Fp_occ[n,:] = (em['frac_occ'][:] * em['Ltot'][:]) / (Rp2**2) # erg s-1 cm-2 cm-1


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


mean_FpFs = np.zeros(n_ph)
mean_FpFs_occ = np.zeros(n_ph)

for n in range(n_ph):
    mean_FpFs[n] = np.mean(FpFs[n,:])
    mean_FpFs_occ[n] = np.mean(FpFs_occ[n,:])

order = np.argsort(phase)
phase = phase[order]
mean_FpFs = mean_FpFs[order]
mean_FpFs_occ = mean_FpFs_occ[order]

# Planetary flux plot
fig = plt.figure()

plt.plot(phase,mean_FpFs[:],label='Unocculted')
plt.scatter(phase,mean_FpFs_occ[:],label='Occulted')

plt.legend()

#plt.yscale('log')

plt.xlabel(r'phase',fontsize=14)
plt.ylabel(r'Band mean F$_{\rm p}$/F$_{\star}$',fontsize=14)


plt.tick_params(axis='both',which='major',labelsize=12)

plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)

#plt.savefig('Fp.png',dpi=300,bbox_inches='tight')
plt.show()
