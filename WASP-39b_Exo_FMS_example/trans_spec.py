# Import required packages
import numpy as np
import matplotlib.pylab as plt
import xarray as xr

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
data = np.loadtxt(fname,skiprows=1)
# Wavelength
wl = data[:,0]
# The transit equation integral approximation
sumf = data[:,1]
# RpRs^2 calculation
RpRs = (Rp**2 + 2.0 * sumf)/Rs**2 + scalfac

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

d1 = plt.scatter(wl_N,RpRs_N,c='grey',marker='o',s=3.0,zorder=1,alpha=0.7)
plt.errorbar(wl_N,RpRs_N,xerr=dwl_N,yerr=RpRs_N_err,fmt='none',c='grey',zorder=1,alpha=0.7,label='G395H')

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

plt.show()

