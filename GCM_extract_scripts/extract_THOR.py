import numpy as np
import h5py
import matplotlib.pylab as plt


# Get surface radius in m from log file (manually)
radius = 79698540.000000

# Filename of regridded data
filename = "regrid_height_hd189b_744_826.h5"

# Read data from HDF5 file
with h5py.File(filename, "r") as f:
    # List all groups
    print("Keys: %s" % f.keys())
    alt_key = list(f.keys())[0]
    cp_key = list(f.keys())[1]
    lat_key = list(f.keys())[2]
    lon_key = list(f.keys())[3]
    Rd_key = list(f.keys())[9]
    Rho_key = list(f.keys())[10]
    T_key = list(f.keys())[12]
    U_key = list(f.keys())[14]
    V_key = list(f.keys())[16]
    W_key = list(f.keys())[18]

    # Get the data
    alt_mid = np.array(list(f[alt_key]))
    cp = np.array(list(f[cp_key]))
    lats = np.array(list(f[lat_key]))
    lons = np.array(list(f[lon_key]))
    Rd = np.array(list(f[Rd_key]))
    Rho = np.array(list(f[Rho_key]))
    T = np.array(list(f[T_key]))
    U = np.array(list(f[U_key]))
    V = np.array(list(f[V_key]))
    W = np.array(list(f[W_key]))


print(np.shape(T))

# THOR is on a height based grid - so outputting is fairly straightforward
# But CMCRT needs the height at the cell levels, so we create a height grid averaging the cell layer points
# and setting the upermost level = upermost layer
# THOR output goes (lat,lon,vert) and vert index starts at alt = 0

nlay = len(alt_mid)
nlev = nlay + 1
nlon = len(lons)
nlat = len(lats)

print(nlay, nlon, nlat)

# Find the altitude mid points
alt = np.zeros(nlev)
for k in range(1,nlay):
  alt[k] = (alt_mid[k] + alt_mid[k-1])/2.0  

alt[-1] = alt_mid[-1]


# Now output gCMCRT height grid in cm
fname = 'THOR.hprf'
print('Outputting height profile: ', fname)
f = open(fname,'w')
for k in range(nlev):
  f.write(str(k+1) + ' ' + str((radius + alt[k]) * 100.0) + '\n')
f.close()


# Find the pressure at each layer from ideal gas law, convert to bar
P = np.zeros((nlat,nlon,nlay))
for i in range(nlat):
    for j in range(nlon):
        for k in range(nlay):
            P[i,j,k] = (Rho[i,j,k]*Rd[i,j,k]*T[i,j,k])/1e5


# Now output variables in CMCRT 1D iprf format (Interpolatable file)
fname = 'THOR.iprf'
print('Outputting interpolatable T-p profile: ', fname)
f = open(fname,'w')
n = 0
for i in range(nlat):
    for j in range(nlon):
        for k in range(nlay):
            f.write(str(n+1) + ' ' +  str(P[i,j,k]) + ' ' + str(T[i,j,k]) + '\n')
            n = n + 1
f.close()

# Output the 3D wind profiles in CMCRT layout -> (lat,lon,lay) in units cm s-1
fname = 'THOR.wprf'
print('Outputting wind profile: ', fname)
f = open(fname,'w')
n = 0
for i in range(nlat):
    for j in range(nlon):
        for k in range(nlay):
            f.write(str(n+1) + ' ' +  str(U[i,j,k] * 100.0) + ' ' + str(V[i,j,k] * 100.0) + ' ' + str(W[i,j,k] * 100.0)  + '\n')
            n = n + 1
f.close()



#fig = plt.figure()

#print(lats[int(nlat/2)],lons[0])
#plt.plot(T[0,0,:],P[0,0,:])
#plt.plot(T[int(nlat/2),int(nlon/2),:],P[int(nlat/2),0,:])
#plt.yscale('log')
#plt.show()
