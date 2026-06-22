# Import required packages
import numpy as np
import matplotlib.pylab as plt
from pathlib import Path

try:
    import xarray as xr
except ImportError:
    xr = None

# Nice fonts (not required)
plt.rc('font', family='sans-serif')
plt.rc('font', serif='Helvetica Neue')
plt.rc('text', usetex='false')

# Stellar radius
Rsun = 6.957e10
Rs =  0.932 * Rsun

# Scaling factor 
# NOTE: using this factor is usually required to fit observations, it is theoretically sound due to the radius-pressure degeneracy
# i.e. you cannot know the exact radius -> pressure relation for the atmosphere and the `base' radius
scalfac = -0.0068

# Import the transmission spectra data
fname = 'Transmission.txt'
head = np.loadtxt(fname,max_rows=1)
# Radius of planet
Rp = head[1]
data = np.atleast_2d(np.loadtxt(fname,skiprows=1))
if data.shape[1] < 4:
    raise ValueError(
        'Transmission.txt must contain wavelength, total, east, and west columns. '
        'Re-run the current 3D_sph_trans code to produce east/west limb output.'
    )
# Wavelength
wl = data[:,0]
# The transit equation integral approximation
sumf = data[:,1]
sumf_east = data[:,2]
sumf_west = data[:,3]
# RpRs^2 calculation
RpRs = (Rp**2 + 2.0 * sumf)/Rs**2 + scalfac
RpRs_east = (Rp**2 + 2.0 * (2.0 * sumf_east))/Rs**2 + scalfac
RpRs_west = (Rp**2 + 2.0 * (2.0 * sumf_west))/Rs**2 + scalfac

has_g395h = xr is not None and Path('ERS_NIRSpec_G395H_weighted-mean-transmission-spectrum_using-DGNELFMANWSBPAR.nc').exists()
if has_g395h:
    # Import the G395H observational data
    fname = 'ERS_NIRSpec_G395H_weighted-mean-transmission-spectrum_using-DGNELFMANWSBPAR.nc'
    data = xr.open_dataset(fname)
    wl_N = data.coords['central_wavelength']
    dwl_N = data.coords['bin_half_width']
    wl_N = wl_N.values
    dwl_N = dwl_N.values
    wl_e = np.zeros(len(wl_N)+1)
    wl_e[:-1] = wl_N[:] - dwl_N[:]
    wl_e[-1] = wl_N[-1] + dwl_N[-1]
    RpRs_N = data.variables['transit_depth']
    RpRs_N_err = data.variables['transit_depth_error']

has_soss = Path('NIRISS_SOSS_FINAL.txt').exists()
if has_soss:
    # Import the NIRISS SOSS data
    fname = 'NIRISS_SOSS_FINAL.txt'
    data = np.loadtxt(fname,delimiter=',')
    wl_S = data[:,0]
    dwl_S = data[:,1]
    RpRs_S = data[:,2]/1e6
    RpRs_S_err = data[:,3]/1e6


# Plot the RpRs^2 figure, comparing the GCM result to the obs data
fig = plt.figure()

plt.plot(wl,RpRs,c='darkmagenta',label='10x Exo-FMS GCM model')

if has_g395h:
    d1 = plt.scatter(wl_N,RpRs_N,c='grey',marker='o',s=3.0,zorder=1,alpha=0.7)
    plt.errorbar(wl_N,RpRs_N,xerr=dwl_N,yerr=RpRs_N_err,fmt='none',c='grey',zorder=1,alpha=0.7,label='G395H')

if has_soss:
    d2 = plt.scatter(wl_S,RpRs_S,c='darkcyan',marker='s',s=3.0,zorder=1,alpha=0.7)
    plt.errorbar(wl_S,RpRs_S,xerr=dwl_S,yerr=RpRs_S_err,fmt='none',c='darkcyan',zorder=1,alpha=0.7,label='SOSS')

plt.xlabel(r'$\lambda$ [$\mu$m]',fontsize=14)
plt.ylabel(r'(R$_{\rm p}$/R$_{\star}$)$^{2}$',fontsize=14)

plt.xscale('log')
plt.xlim(0.3,6.0)

plt.legend()

xticks = [0.3,0.5,1,2,3,4,5,6]
xticks_lab = ['0.3','0.5','1','2','3','4','5','6']

plt.xticks(xticks,xticks_lab)

plt.tick_params(axis='both',which='major',labelsize=12)

plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)

#plt.savefig('Exo-FMS_WASP-39b_trans_spec.png',dpi=300,bbox_inches='tight')

fig = plt.figure()

plt.plot(wl,RpRs,c='black',label='Total')
plt.plot(wl,RpRs_east,c='firebrick',label='East limb')
plt.plot(wl,RpRs_west,c='royalblue',label='West limb')

plt.xlabel(r'$\lambda$ [$\mu$m]',fontsize=14)
plt.ylabel(r'(R$_{\rm p}$/R$_{\star}$)$^{2}$',fontsize=14)

plt.xscale('log')
plt.xlim(0.3,6.0)

plt.legend()

plt.xticks(xticks,xticks_lab)

plt.tick_params(axis='both',which='major',labelsize=12)

plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)

#plt.savefig('Exo-FMS_WASP-39b_trans_spec_limbs.png',dpi=300,bbox_inches='tight')

plt.show()
