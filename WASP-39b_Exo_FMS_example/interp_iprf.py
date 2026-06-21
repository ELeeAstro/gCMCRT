import numpy as np
from scipy.interpolate import RegularGridInterpolator


def make_regular_interpolator(y_grid, x_grid, values):
    return RegularGridInterpolator(
        (y_grid, x_grid),
        values,
        method='linear',
        bounds_error=False,
        fill_value=None,
    )

# Read in iprf
iprf = 'FMS.iprf'
data = np.loadtxt(iprf)
idx = data[:,0]
P = np.log10(data[:,1])
T = np.log10(data[:,2])
ni = len(idx)

# Read in processed interpolation file

ifile = '../data/CE_interp/interp_FastChem_10x_cond.txt'
f = open(ifile,'r')
dim = f.readline().split()

nT = int(dim[0])
nP = int(dim[1])
nL = int(dim[2])
nsp = int(dim[3])

print(nT, nP, nL, nsp)

sp = f.readline().split()

print(sp)

iT = f.readline().split()
iT = np.array(iT,dtype=float)
iP = f.readline().split()
iP = np.array(iP,dtype=float)

print(iT)
print(iP)

lT = np.log10(iT)
lP = np.log10(iP)

# Read in mu and VMR data using loadtxt
data = np.loadtxt(ifile,skiprows=4)
mu_1D = np.log10(data[:,0])
VMR_1D = np.zeros((nL,nsp))
VMR_1D = np.log10(data[:,1:])

# Convert 1D to regular P-T grids. The interpolation file has T as the
# inner loop, so a plain reshape gives arrays indexed as [P, T].
mu = mu_1D.reshape(nP, nT)
VMR = VMR_1D.reshape(nP, nT, nsp)

# Interpolate all profile points in one vectorized call. VMR has species as a
# trailing dimension, so RegularGridInterpolator returns every species at once.
points = np.column_stack((P, T))
mu_1D = make_regular_interpolator(lP, lT, mu)(points)
VMR_1D = make_regular_interpolator(lP, lT, VMR)(points)
#for s in range(nsp):
#  VMR_1D[:,s] = np.log10(10.0**VMR_1D[:,s]/np.sum(10.0**VMR_1D,axis=1))

sp = ['H2','H','H-','e-','He','H2O','CO','CO2','CH4','TiO','VO','Fe','Fe+','FeH','SiO','Na','K','NH3','H2S','HCN','HCl','PH3','OH','C2H2','HF','SH','V','V+','Ti','Ti+']

# Create .prf file
head = open('../data/header.txt','r')
lines = head.readlines()

fname = 'FMS.prf'
prf = open(fname,'w')
prf.write(lines[0])
prf.write(lines[1])
prf.write(str(ni) + '\n')
prf.write(lines[2])
prf.write(str(nsp) + '\n')
for n in range(nsp):
    prf.write(sp[n] + '\n')
prf.write(lines[3])
prf.write(lines[4])
for n in range(ni):
    prf.write(str(n+1) + ' ' + str(10.0**(P[n])) + ' ' + str(10.0**(T[n])) + ' ' + str(10.0**(mu_1D[n])) + ' ' + " ".join(str(l) for l in 10.0**(VMR_1D[n,:])) + '\n')
