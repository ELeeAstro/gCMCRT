import numpy as np
import matplotlib.pylab as plt
from scipy.io import FortranFile
from matplotlib import cm

plt.rc('font', family='sans-serif')
plt.rc('font', serif='Helvetica Neue')
plt.rc('text', usetex='false')

n_ph = 2
nwl = 503

Fp = np.zeros((n_ph,nwl))
Ltot = np.zeros((n_ph,nwl))
ph = np.zeros(n_ph)

for n in range(n_ph):
  fname = 'Em_'+str(n+1).zfill(3)+'.txt'
  print(fname)
  head = np.loadtxt(fname,max_rows=1)
  ph[n] = int(head[3])
  Rp = float(head[1])
  Rp2 = float(head[2])
  data = np.loadtxt(fname,skiprows=1)
  wl = data[:,0]
  frac = data[:,1]
  Ltot[n,:] = data[:,2]
  Fp[n,:] = (frac[:] * Ltot[n,:]) / (Rp2**2)

xpix = 300
ypix = 300

# Find area of a pixel
area = ((Rp2*1.01)/xpix) * ((Rp2*1.01)/ypix)

print(area,4.0*np.pi*Rp2**2,area/(4.0*np.pi*Rp2**2))

c0 = 2.99792458e10
h = 6.626176e-27
kb = 1.3080662e-16

#read in Spitzer 3.6 um throughput
fname = 'Spitzer_36_response.txt'
data = np.loadtxt(fname)
wl_S = data[:,0]
T1_S = data[:,1]

T1_S_int = np.interp(wl,wl_S,T1_S)
T1_S_bot = np.trapz(T1_S_int,wl)

# Read in the pixel data

f_pix = np.zeros((xpix, ypix))
wf_pix = np.zeros((nwl,xpix,ypix))
BT = np.zeros((n_ph, xpix, ypix))

for n in range(n_ph):
  fname = 'f_im_'+str(n+1).zfill(3)+'.txt'
  print(fname)
  f = FortranFile(fname, 'r')
  for l in range(nwl):
    # Flux in each pixel
    f_pix[:,:] = (f.read_record(np.float32).reshape((xpix,ypix),order='F') * Ltot[n,l])/area 
    # Calculate brightness temperature in each pixel
    wl_cm = wl[l] * 1e-4 
    wf_pix[l,:,:] = (h*c0)/(kb*wl_cm) / np.log(1.0 +  ((2.0*h*c0**2)/(f_pix[:,:]/np.pi * wl_cm**5))) 
  # Integrate brightness temperature through Spizter 4.5um bandpass
  for x in range(xpix):
      for y in range(ypix):
          BT[n,x,y] = np.trapz(wf_pix[:,x,y]*T1_S_int[:],wl)/T1_S_bot


minval = np.min(BT[:,:,:][np.nonzero(BT[:,:,:])])
maxval = np.max(BT[:,:,:][np.nonzero(BT[:,:,:])])
print(minval,maxval)
for n in range(n_ph):
  fig = plt.figure()
  cmap = cm.get_cmap('RdYlBu_r').copy()
  cmap.set_under(color='black')
  pmap = plt.imshow(BT[n,:,:],vmin=500, vmax=4000,cmap=cmap)
  cbar = plt.colorbar(pmap,extend='min')
  cbar.ax.tick_params(labelsize=14)
  cbar.set_label(r'T$_{\rm b}$ [K]',fontsize=14)
  plt.ylabel(r'y pixel',fontsize=14)
  plt.xlabel(r'x pixel',fontsize=14)
  plt.title(r'Longitude: ' + str(ph[n]),fontsize=14)

  plt.tick_params(axis='both', which='major', labelsize=12)
  plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)

  #plt.savefig('W33b_Tb_'+str(n)+'.png',dpi=300,bbox_inches='tight')

plt.show()
