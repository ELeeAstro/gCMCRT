import numpy as np
from netCDF4 import Dataset
import matplotlib.pylab as plt

# Filename to import
fname = 'MITgcm_HD20-like.nc'

# Indexing goes from 0 (lowest height) - nlev (highest hight) in the ExoRad netCDF

# Hypsometric equation o find height grid (constant gravity version)
def calc_alt_col(nlev,Plev,T,grav,Rd):
    alt = np.zeros(nlev)
    alt[0] = 0.0
    for k in range(nlay):
        alt[k+1] = alt[k] + Rd/grav * T[k] * np.log(Plev[k]/Plev[k+1])
    return alt


it = -1

nc = Dataset(fname)

# Print the keys
print (nc.file_format)
print (nc.dimensions.keys())
print (nc.variables.keys())

# Print the time
time = nc.variables['time'][:]
print ('time: ',time)

# Get the longitudes, latitudes and pressure levels
lons = nc.variables['XC'][:]
lats = nc.variables['YC'][:]
lay = nc.variables['Z'][:]

# Get dimensions of system
nt = len(time)
nlat = len(lats)
nlon = len(lons)
nlay = len(lay)
nlev = nlay + 1

print('lats: ', nlat, lats[:])
print('lons: ', nlon, lons[:])
print('lay: ', nlay, lay[:])

# Get the global attributes
p_ref = nc.getncattr('p_ref')
cp = nc.getncattr('cp')
R = nc.getncattr('R')
g = nc.getncattr('g')
radius = nc.getncattr('radius')

print('p_ref: ',p_ref)
print('cp: ',cp)
print('R: ',R)
print('g: ',g)
print('radius: ',radius)

# Extract temperature
T = np.zeros((nlay,nlat,nlon))
T = nc.variables['T'][it,:,:,:]

# Extract winds
U = np.zeros((nlay,nlat,nlon))
U = nc.variables['U'][it,:,:,:]
V = np.zeros((nlay,nlat,nlon))
V = nc.variables['V'][it,:,:,:]
W = np.zeros((nlay,nlat,nlon))
W = nc.variables['W'][it,:,:,:]

# Shift the latitude coordinate to start from 0-180, the gCMCRT coordinate system
lats[:] = lats[:] + 90.0
# Now we've got to shift longitudes to make the 0th index start at 0 degrees
# gCMCRT works from 0-360 degree longitude coordinates

nr = int(nlon/2)
lons[:] = np.roll(lons[:],nr)
lons[:] = np.where(lons[:] > 0, lons[:], 360.0-abs(lons[:]))

print(lons[:])

# Now we've got to shift all arrays by nr in the longitude axis for each pressure level and latitude
for j in range(nlat):
    for k in range(nlay):
        T[k,j,:] = np.roll(T[k,j,:],nr)
        U[k,j,:] = np.roll(U[k,j,:],nr)
        V[k,j,:] = np.roll(V[k,j,:],nr)
        W[k,j,:] = np.roll(W[k,j,:],nr)


# All variables are now in the CMCRT coordinate system - now perform height calculation and interpolation
# In the netCDF file, we are given the pressure at the layers, for accurate height calculation
# we need the pressure at the levels - it looks like it's a stright average at the lower boundary
# so for now use a stright average and assume lowest pressure = lowest layer pressure
Play = lay
Plev = np.zeros(nlev)
Plev[0] = p_ref

for k in range(1,nlay):
    Plev[k] = (Play[k-1] + Play[k])/2.0
Plev[-1] = Play[-1]

# All variables are ready for the height calculation
# perform altitude calculation for all lat and lon
alt_3D = np.zeros((nlev,nlat,nlon))
alt_3D_mid = np.zeros((nlay,nlat,nlon))
for j in range(nlat):
    for i in range(nlon):
      alt_3D[:,j,i] = calc_alt_col(nlev,Plev[:],T[:,j,i],g,R)
      alt_3D_mid[:,j,i] = (alt_3D[0:-1,j,i] + alt_3D[1:,j,i])/2.0
#find the maximum altitude profile and use that as the vertical grid
imax = np.unravel_index(np.argmax(alt_3D, axis=None), alt_3D.shape)
print(imax)
# alt is the altitude grid used for gCMCRT
alt = alt_3D[:,imax[1],imax[2]]
# Find midpoint altitudes to interpolate to 
alt_mid = np.zeros(nlay)
alt_mid = (alt[0:-1] + alt[1:])/2.0

# Now find the interpolated temperature and pressure at each altitude midpoint
# We assume that the temperature is constant beyond the uppermost alitude boundary
iT = np.zeros((nlay,nlat,nlon))
iP = np.zeros((nlay,nlat,nlon))
for j in range(nlat):
    for i in range(nlon):
        iT[:,j,i] = np.interp(alt_mid[:],alt_3D_mid[:,j,i],T[:,j,i])
        iP[:,j,i] = 10.0**np.interp(alt_mid[:],alt_3D_mid[:,j,i],np.log10(Play[:]),right=-7,left=-7)
        # Now perform pressure extrapolation in height for when max alt < alt_grid in a column
        # approximate this as a hydrostatic layer
        for k in range(nlay):
          if (iP[k,j,i] == 1e-7):
            p0 = iP[k-1,j,i]
            z0 = alt_mid[k-1]
            H0 = R * iT[k,j,i] / g
            iP[k,j,i] = p0 * np.exp(-(alt_mid[k] - z0)/H0)
            if (iP[k,j,i] < 1e-7):
              iP[k,j,i] = 1e-7

# Now output T-p profile in bar and K to an interpolation file (iprf) - after which we can interpolate to GGChem values to get VMRs
fname = 'MITgcm.iprf'
print('Outputting interpolatable T-p profile: ', fname)
f = open(fname,'w')
n = 0
for j in range(nlat):
    for i in range(nlon):
        for k in range(nlay):
            f.write(str(n+1) + ' ' +  str(iP[k,j,i]/1e5) + ' ' + str(iT[k,j,i]) + '\n')
            n = n + 1
f.close()

# Now output Height grid in cm
fname = 'MITgcm.hprf'
print('Outputting height profile: ', fname)
f = open(fname,'w')
for k in range(nlev):
  f.write(str(k+1) + ' ' + str((radius + alt[k]) * 100.0) + '\n')
f.close()


# Now interpolate the wind profiles to the new height grid
iU = np.zeros((nlay,nlat,nlon))
iV = np.zeros((nlay,nlat,nlon))
iW = np.zeros((nlay,nlat,nlon))
for j in range(nlat):
    for i in range(nlon):
        iU[:,j,i] = np.interp(alt_mid[:],alt_3D_mid[:,j,i],U[:,j,i])
        iV[:,j,i] = np.interp(alt_mid[:],alt_3D_mid[:,j,i],V[:,j,i])
        iW[:,j,i] = np.interp(alt_mid[:],alt_3D_mid[:,j,i],W[:,j,i])


#Output the 3D wind profiles in CMCRT layout -> (lat,lon,lay) in units cm s-1
fname = 'MITgcm.wprf'
print('Outputting wind profile: ', fname)
f = open(fname,'w')
n = 0
for j in range(nlat):
    for i in range(nlon):
        for k in range(nlay):
            f.write(str(n+1) + ' ' +  str(iU[k,j,i] * 100.0) + ' ' + str(iV[k,j,i] * 100.0) + ' ' + str(iW[k,j,i] * 100.0)  + '\n')
            n = n + 1
f.close()

