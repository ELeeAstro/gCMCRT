import numpy as np
import matplotlib.pylab as plt
from scipy.io import FortranFile
from matplotlib import ticker, cm
import seaborn as sns


# Read wavelengths from the Em file!!!
wl = np.loadtxt('Em_001.txt',skiprows=1)
wl = wl[:,0]

nwl = len(wl)
nlat = 30
nlon = 64
nlay = 53

ni = nlat * nlon * nlay

f = FortranFile('cf_001.txt','r')

cf = np.zeros((nwl,nlay,nlon,nlat))
for l in range(nwl):
  #cf_1D = np.fromfile('cf.out',count=ni,dtype=np.float32)
  cf_1D = f.read_reals(np.float32)
  cf[l,:,:,:] = np.reshape(cf_1D,[nlay,nlon,nlat],order='F')/np.sum(cf_1D[:])

prf = np.loadtxt('SPARC-MITgcm.prf',skiprows=27)
P_1D = prf[:,1]

P = np.zeros((nlay,nlon,nlat))
P[:,:,:] = np.reshape(P_1D,[nlay,nlon,nlat],order='F')


fig = plt.figure()

lwl = [0.5,1.0,5.0,10.0,20.0]
iwl = np.searchsorted(wl,lwl)
print(iwl)
ilat = int(nlat/2)
ilon = 0

for i in range(len(iwl)):
  plt.plot(cf[iwl[i],:,ilon,ilat],P[:,ilon,ilat],label=str(wl[iwl[i]])+' um')

plt.legend()

plt.ylim(1e-6,100)
plt.xscale('log')
plt.yscale('log')

plt.gca().invert_yaxis()

plt.ylabel(r'P [bar]',fontsize=14)
plt.xlabel(r'Fractional Contribution',fontsize=14)

yticks = [100,10,1,0.1,0.01,1e-3,1e-4,1e-5,1e-6]
yticks_lab = ['100','10','1','0.1','0.01','10$^{-3}$','10$^{-4}$','10$^{-5}$','10$^{-6}$']

plt.yticks(yticks,yticks_lab)

plt.tick_params(axis='both',which='major',labelsize=12)

plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)

#plt.savefig('3D_cont_sub_stellar_0G_wl.png',dpi=300,bbox_inches='tight')

fig = plt.figure()

#cm = sns.color_palette("rocket", as_cmap=True)
cm = sns.color_palette("mako", as_cmap=True)

cf_swap = np.swapaxes(cf[:,:,ilon,ilat],0,1)
lev = np.linspace(-7,np.log10(np.max(cf_swap)),20)

con = plt.contourf(wl,P[:,ilon,ilat],np.log10(cf_swap[:,:]),levels=lev,extend='both',cmap=cm)
for c in con.collections:
    c.set_edgecolor("face")


cb = fig.colorbar(con)
cb.solids.set_rasterized(True)


plt.ylim(1e-6,100)
plt.xscale('log')
plt.yscale('log')

plt.gca().invert_yaxis()

plt.ylabel(r'P [bar]',fontsize=14)
plt.xlabel(r'$\lambda$ [$\mu$m]',fontsize=14)

yticks = [100,10,1,0.1,0.01,1e-3,1e-4,1e-5,1e-6]
yticks_lab = ['100','10','1','0.1','0.01','10$^{-3}$','10$^{-4}$','10$^{-5}$','10$^{-6}$']

plt.yticks(yticks,yticks_lab)


xticks = [0.2,0.5,1,2,3,4,5,6,7,10,20,30]
xticks_lab = ['0.2','0.5','1','2','3','4','5','6','7','10','20','30']

plt.xticks(xticks,xticks_lab)


plt.tick_params(axis='both',which='major',labelsize=12)

plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)

#plt.savefig('3D_cont_sub_stellar_0G.png',dpi=300,bbox_inches='tight')

plt.show()
