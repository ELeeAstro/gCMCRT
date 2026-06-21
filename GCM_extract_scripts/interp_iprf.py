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

simname = 'SPARC-MITgcm' 

# Read in iprf
iprf = simname+'.iprf'
data = np.loadtxt(iprf)
idx = data[:,0]
P = data[:,1]
T = data[:,2]
ni = len(idx)

# Read in processed interpolation file

ifile = '../data/CE_interp/UHJ_interp_co_1x.txt'
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

# Read in mu and VMR data using loadtxt
data = np.loadtxt(ifile,skiprows=4)
mu_1D = data[:,0]
VMR_1D = np.zeros((nL,nsp))
VMR_1D = data[:,1:]

# Convert 1D to regular P-T grids. The interpolation file has T as the
# inner loop, so a plain reshape gives arrays indexed as [P, T].
mu = mu_1D.reshape(nP, nT)
VMR = VMR_1D.reshape(nP, nT, nsp)

# Interpolate all profile points in one vectorized call. VMR has species as a
# trailing dimension, so RegularGridInterpolator returns every species at once.
points = np.column_stack((P, T))
mu_1D = make_regular_interpolator(iP[::-1], iT[::-1], mu)(points)
VMR_1D = 10.0**make_regular_interpolator(iP[::-1], iT[::-1], VMR)(points)

# Create .prf file
head = open('../data/header.txt','r')
lines = head.readlines()

fname = simname+'.prf'
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
    prf.write(str(n+1) + ' ' + str(P[n]) + ' ' + str(T[n]) + ' ' + str(mu_1D[n]) + ' ' + " ".join(str(l) for l in VMR_1D[n,:]) + '\n')
