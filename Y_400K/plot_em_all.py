import numpy as np
import matplotlib.pylab as plt
import glob, os
import seaborn as sns
from astropy.timeseries import LombScargle

# Catch all file with same name

files = glob.glob('./Em_*')
files = sorted(files,key=os.path.getmtime)

files.pop(-1)

print(files)

ns = len(files)

print(ns)

nwl = 503 
wl = np.zeros(nwl)
Fp = np.zeros((ns,nwl))

for i in range(ns):
  fname = files[i]
  head = np.loadtxt(fname,max_rows=1)
  Rp_b = float(head[1])
  Rp_t = float(head[2])
  data = np.loadtxt(fname,skiprows=1)
  wl[:] = data[:,0]
  frac = data[:,1]
  Ltot = data[:,2]
  Fp[i,:] = (frac[:] * Ltot[:]) / (Rp_b**2)

hrs = np.arange(0,ns,1)


# Read spectral data
fname = 'Beiler_2023.txt'
data = np.loadtxt(fname)
wl_B = data[:,0]
Jy_B = data[:,1]
Jy_B_err = data[:,2]


c = 29979245800.0
Hz = c/(wl*1e-4)

#flx_B = Jy_B * 1.0e23 * Hz / (wl_B*1e-4) #erg/cm^2/s/cm
#flx_B_err = Jy_B_err * 1.0e23 *  Hz/ (wl_B*1e-4)  #erg/cm^2/s/cm


D = 13.57 * 3.0857e18 # cm

m_wl = np.array([15.0, 18.0, 21.0])
m_wl_err = np.array([1.80,2.95,4.58])
m_Jy = np.array([0.1231e-3, 0.1019e-3, 0.0731e-3])
m_Jy_err = np.array([0.0062e-3, 0.0051e-3, 0.0037e-3])

# Convert model to Jy
for i in range(ns):
    Fp[i,:] = Fp[i,:] * Rp_b**2 / D**2 * 1e23 * (wl[:]*1e-4) / Hz[:] * 6.5e-1


fig = plt.subplots(figsize=(12.8, 4.8))

col = sns.color_palette('colorblind')

i = 0
plt.plot(wl[:],Fp[0,:]*1e3,label=r'GCM',c=col[1])
for i in range(1,ns):
    plt.plot(wl[:],Fp[i,:]*1e3,c=col[1],lw=1)

plt.plot(wl_B[:],Jy_B[:]*1e3,ls='solid',c=col[0],lw=1,label=r'WISE 0359-54')

plt.errorbar(m_wl,m_Jy*1e3,yerr=m_Jy_err*1e3,xerr=m_wl_err,fmt='o',c=col[0])

plt.xscale('log')
#plt.yscale('log')

plt.xlabel(r'$\lambda$ [$\mu$m]',fontsize=16)
plt.ylabel(r'$f_{\lambda}$ [mJy]',fontsize=16)


plt.legend()

plt.xlim(1.0,25)
plt.ylim(0.0,0.30)


xticks = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]
xticks_lab = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','','18','','20','','22','','24','']
plt.xticks(xticks,xticks_lab)

plt.tick_params(axis='both',which='major',labelsize=14)

plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)

plt.savefig('obs_comp.pdf',dpi=144,bbox_inches='tight')


# Integrate flux in wavelength
Tbol = np.zeros(ns)
sb = 5.670374419e-5
for i in range(ns):
  Tbol[i] = (np.trapz(Fp[i,:],wl*1e-4))
  print(i,(Tbol[i]/sb)**0.25)


Fp_av = np.zeros(nwl)
for i in range(nwl):
  Fp_av[i] = np.average(Fp[:,i])
  #Fp_av[i] = np.median(Fp[:,i])


# Planetary average spectral variation J, H, Ks and Spitzer
fig = plt.figure()
b_file = ['Spitzer_1.txt','Spitzer_2.txt']
b_lab = [r'Spitzer [3.6]',r'Spitzer [4.5]']
nb = len(b_file)
col = sns.color_palette("husl", nb)
# Read band information and convolve
for n in range(nb):
  print(n, nb)
  fname = 'filters/'+b_file[n]
  data = np.loadtxt(fname)
  wl_b = data[:,0]/1e4
  t_b = data[:,1]
  Fp_b = np.zeros(ns)
  for i in range(ns):
    # Now perform filter convolution
    i_Fp = np.interp(wl_b,wl,Fp[i,:])
    Fp_b[i] = np.trapz(i_Fp*t_b,wl_b)/np.trapz(t_b,wl_b)
  Fp_b[:] = Fp_b[:]/np.average(Fp_b[:])
  plt.plot(hrs,Fp_b[:],label=b_lab[n],c=col[n])

plt.hlines(1.0,hrs[0],hrs[-1],colors='black',ls='dashed')
plt.legend(title='Photometry Band',loc='lower left',ncol=3)
plt.xlabel(r'Time [hrs]',fontsize=14)
plt.ylabel(r'Relative Flux Density',fontsize=14)
plt.xlim(hrs[0],hrs[-1])
#plt.ylim(0.9,1.1)
#xticks = [0.2,0.5,1,2,3,4,5,6,8,10,20,30]
#xticks_lab = ['0.2','0.5','1','2','3','4','5','6','8','10','20','30']
#plt.xticks(xticks,xticks_lab)
plt.tick_params(axis='both',which='major',labelsize=12)
plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)
plt.savefig('Y_400_var.pdf',dpi=144,bbox_inches='tight')

plt.show()

quit()

plt.figure()

w = np.logspace(-6, -4, 10000)

# Try get periodigram out
for n in range(0):
  frequency, power = LombScargle(hrs, Fp_b).autopower()
  frequency.unit
  power.unit
  plt.plot(frequency, power)

plt.xlabel('Angular frequency [rad/s]')
plt.ylabel('Normalized amplitude')

plt.show()


