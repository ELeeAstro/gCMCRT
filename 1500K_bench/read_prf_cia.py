import numpy as np
import matplotlib.pylab as plt
from matplotlib.pyplot import cm
from scipy.io import FortranFile

wl = np.loadtxt('wavelengths.wl',skiprows=1)

wl = wl[:,1]

k_tab = FortranFile('CIA.cmcrt', 'r')

fig = plt.figure()
ax = fig.add_subplot(111)


k_ints =  k_tab.read_ints(dtype=np.int32)

NX = k_ints[0]
n_bins = k_ints[1]

k_data = np.zeros((n_bins,NX))
for n in range(n_bins):
    k_data[n,:] = k_tab.read_reals(dtype=np.float32)

color=iter(cm.rainbow(np.linspace(0,1,NX)))
for x in range(0,NX):
   c=next(color)
   ax.plot(wl[:n_bins],k_data[:,x],c=c)
   ax.plot(wl[:n_bins],k_data[:,x],c=c,ls='dashed')



plt.yscale('log')
plt.xscale('log')
plt.show()

