import numpy as np
from scipy import interpolate

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

# Convert 1D to 3D arrays
mu = np.zeros((nP, nT))
VMR = np.zeros((nP, nT, nsp))

# Loop in the GGchem way to get 3D arrays in correct T and P values (T inner loop)
n = 0
for i in range(nP):
    for j in range(nT):
      mu[i,j] = mu_1D[n]
      VMR[i,j,:] = VMR_1D[n,:]
      n = n + 1

# Create a scipy interpolation function for mu
f_mu = interpolate.interp2d(lT[:], lP[:], mu, kind='linear')

# Create a scipy interpolation function for the VMR for each species
f_VMR = []
for i in range(nsp):
    f_VMR.append(interpolate.interp2d(lT[:], lP[:], VMR[:,:,i], kind='linear'))

# Now we can loop across the whole 1D profile and interpolate to each P-T point in the GCM
mu_1D = np.zeros(ni)
VMR_1D = np.zeros((ni,nsp))
for n in range(ni):
    # interpolate mu
    mu_1D[n] = f_mu(T[n],P[n])
    #print(n,mu_1D[n],P[n],T[n])
    # interpolate VMRs
    for s in range(nsp):
        VMR_1D[n,s] = f_VMR[s](T[n],P[n])
        #print(n,s,10.0**VMR_1D[n,s],P[n],T[n])
    #for s in range(nsp):
    #  VMR_1D[n,s] = np.log10(10.0**VMR_1D[n,s]/np.sum(10.0**VMR_1D[n,:]))

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
