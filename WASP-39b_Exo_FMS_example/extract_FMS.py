import numpy as np
from netCDF4 import Dataset

g = 4.26
Rd = 3221.0
R = 8.31446261815324
Rp = 91438268.0

def calc_alt_col_var_g(nlev,Plev,T,gs,Rd,Rp):
    alt = np.zeros(nlev)
    alt_new = np.zeros(nlev)

    # Initial guess with constant gravity
    alt[0] = 0.0
    for k in range(nlev-1):
        alt[k+1] = alt[k] + (Rd*T[k])/gs * np.log(Plev[k]/Plev[k+1])

    alt_new[:] = alt[:]
    # Converge using alt[k] for gravity
    itera = 0
    converge = False
    while (converge == False) :

      # Current delta alt
      atdepth = alt[-1] - alt[0]

      # Perform hydrostatic calcuation
      alt_new[0] = 0.0
      for k in range(nlev-1):
        grav = gs * (Rp / (Rp + alt[k]))**2
        alt_new[k+1] = alt_new[k] + (Rd*T[k])/grav * np.log(Plev[k]/Plev[k+1])

      # New delta alt
      atdepth1 = alt_new[-1] - alt_new[0]

      # Convergence test
      itera = itera + 1
      xdepth  = 100.0 * abs((atdepth1-atdepth)/atdepth)
      #print(itera, atdepth, atdepth1, xdepth)
      if (xdepth < 1e-8):
        converge = True

      alt[:] = alt_new[:]

    return alt

def calc_alt_col(nlev,Plev,T,grav,Rd):
    alt = np.zeros(nlev)
    alt[0] = 0.0
    for k in range(nlay):
        alt[k+1] = alt[k] + (Rd*T[k])/grav * np.log(Plev[k]/Plev[k+1])
    return alt


nc = Dataset('atmos_average_W39_10x_1010.nc')

time = nc.variables['time'][:]
lats = nc.variables['grid_yt'][:]
lons = nc.variables['grid_xt'][:]
lev_f = nc.variables['pfull'][:]

pk = nc.variables['pk'][:]
bk = nc.variables['bk'][:]

nt = len(time)
nlat = len(lats)
nlon = len(lons)
nlay = len(lev_f)
nlev = nlay + 1

ni = nlat * nlon * nlay
print(nlay,nlat,nlon,ni)


it = -1

# Extract temperature
T = np.zeros((nlay,nlat,nlon))
T = nc.variables['temp'][it,:,:,:]

# Get level pressure grid from bk and pk
Ps = np.zeros((nlat,nlon))
Ps = nc.variables['ps'][it,:,:]

Plev = np.zeros((nlev,nlat,nlon))
for j in range(nlat):
  for i in range(nlon):
    for k in range(nlev):
      Plev[k,j,i] = Ps[j,i] * bk[k] + pk[k]

#Invert T and P
for j in range(nlat):
    for i in range(nlon):
      T[:,j,i] = T[::-1,j,i] 
      Plev[:,j,i] = Plev[::-1,j,i]

# Calculate layer P
Play = np.zeros((nlay,nlat,nlon))
for j in range(nlat):
    for i in range(nlon):
        for k in range(nlay):
          Play[k,j,i] = (Plev[k,j,i] - Plev[k+1,j,i])/(np.log(Plev[k,j,i]/Plev[k+1,j,i]))

# All variables are ready for the height calculation
# perform altitude calculation for all lat and lon
alt_3D = np.zeros((nlev,nlat,nlon))
alt_3D_mid = np.zeros((nlay,nlat,nlon))
for j in range(nlat):
    for i in range(nlon):
      alt_3D[:,j,i] = calc_alt_col(nlev,Plev[:,j,i],T[:,j,i],g,Rd)
      #alt_3D[:,j,i] = calc_alt_col_var_g(nlev,Plev[:,j,i],T[:,j,i],g,Rd,Rp)
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
        iP[:,j,i] = 10.0**np.interp(alt_mid[:],alt_3D[:,j,i],np.log10(Plev[:,j,i]),right=-7,left=-7)
        # Now perform pressure extrapolation in height for when max alt < alt_grid in a column
        # approximate this as a hydrostatic layer
        for k in range(nlay):
          if (iP[k,j,i] == 1e-7):
            p0 = iP[k-1,j,i]
            z0 = alt_mid[k-1]
            H0 = Rd * iT[k,j,i] / g
            iP[k,j,i] = p0 * np.exp(-(alt_mid[k] - z0)/H0)
            if (iP[k,j,i] < 1e-7):
              iP[k,j,i] = 1e-7


# Now output T-p profile in bar and K to an interpolation file (iprf) - after which we can interpolate to GGChem values to get VMRs
fname = 'FMS.iprf'
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
fname = 'FMS.hprf'
print('Outputting height profile: ', fname)
f = open(fname,'w')
for k in range(nlev):
  f.write(str(k+1) + ' ' + str((Rp + alt[k]) * 100.0) + '\n')
f.close()



# Extract winds
U = np.zeros((nlay,nlat,nlon))
U = nc.variables['ucomp'][it,:,:,:]
V = np.zeros((nlay,nlat,nlon))
V = nc.variables['vcomp'][it,:,:,:]
W = np.zeros((nlay,nlat,nlon))
W = nc.variables['omega'][it,:,:,:]

# Invert and change omega to vertical velocity
rho = np.zeros(nlay)
for j in range(nlat):
    for i in range(nlon):
        U[:,j,i] = U[::-1,j,i]
        V[:,j,i] = V[::-1,j,i]
        W[:,j,i] = W[::-1,j,i]
        rho[:] = Play[:,j,i]/(Rd * T[:,j,i])
        W[:,j,i] = -W[:,j,i]/(rho[:]*g)

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
fname = 'FMS.wprf'
print('Outputting wind profile: ', fname)
f = open(fname,'w')
n = 0
for j in range(nlat):
    for i in range(nlon):
        for k in range(nlay):
            f.write(str(n+1) + ' ' +  str(iU[k,j,i] * 100.0) + ' ' + str(iV[k,j,i] * 100.0) + ' ' + str(iW[k,j,i] * 100.0)  + '\n')
            n = n + 1
f.close()
