import sys
import numpy as np
from netCDF4 import Dataset

g = 177.8
Rd = 3568.0
R = 8.31446261815324
Rp = 71492000.0

p0 = 1.0 * 1e5

#print('Number of arguments:', len(sys.argv), 'arguments.')
#print('Argument List:', str(sys.argv))

print(sys.argv[1])

t = sys.argv[1]


def calc_alt_col_var_g_p0(nlev,Plev,Play,T,gs,Rd,Rp,p0):
    alt = np.zeros(nlev)
    alt_new = np.zeros(nlev)

    # Initial guess with constant gravity
    ia = np.searchsorted(Plev[::-1],p0)
    ia = len(Plev) - ia - 1
    alt0 = 0.0

    T0 = np.interp(p0, Play[::-1], T[::-1])

    # Do downward sweep
    alt[ia] = alt0 + (Rd*T0)/gs * np.log(p0/Plev[ia])
    for k in range(ia,-1,-1):
        alt[k] = alt[k+1] + (Rd*T[k])/gs * np.log(Plev[k+1]/Plev[k])
        #print(k,alt[k],Plev[k])

    # Do upward sweep
    alt[ia+1] = alt0 + (Rd*T0)/gs * np.log(p0/Plev[ia+1])
    for k in range(ia+1,nlev-1):
        alt[k+1] = alt[k] + (Rd*T[k])/gs * np.log(Plev[k]/Plev[k+1])
        #print(k,alt[k],Plev[k])

    alt_new[:] = alt[:]
    # Converge using alt[k] for gravity
    itera = 0
    converge = False
    while (converge == False) :

      # Current delta alt
      atdepth = alt[-1] - alt[0]

      # Do downward sweep
      alt_new[ia] = alt0 + (Rd*T0)/gs * np.log(p0/Plev[ia])
      for k in range(ia,-1,-1):
        grav = gs * (Rp / (Rp + alt[k]))**2
        alt[k] = alt[k+1] + (Rd*T[k])/grav * np.log(Plev[k+1]/Plev[k])
        #print(k,alt[k],Plev[k])

      # Do upward sweep
      alt[ia+1] = alt0 + (Rd*T0)/gs * np.log(p0/Plev[ia+1])
      for k in range(ia+1,nlev-1):
        grav = gs * (Rp / (Rp + alt[k]))**2
        alt[k+1] = alt[k] + (Rd*T[k])/grav * np.log(Plev[k]/Plev[k+1])
        #print(k,alt[k],Plev[k])

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

def calc_alt_col_p0(nlev,Plev,Play,T,gs,Rd,Rp,p0):
    alt = np.zeros(nlev)
    alt_new = np.zeros(nlev)

    # Initial guess with constant gravity
    ia = np.searchsorted(Plev[::-1],p0)
    ia = len(Plev) - ia - 1 
    alt0 = 0.0

    T0 = np.interp(p0, Play[::-1], T[::-1])

    # Do downward sweep
    alt[ia] = alt0 + (Rd*T0)/gs * np.log(p0/Plev[ia])
    for k in range(ia,-1,-1):
        alt[k] = alt[k+1] + (Rd*T[k])/gs * np.log(Plev[k+1]/Plev[k])
        #print(k,alt[k],Plev[k])    

    # Do upward sweep
    alt[ia+1] = alt0 + (Rd*T0)/gs * np.log(p0/Plev[ia+1])
    for k in range(ia+1,nlev-1):
        alt[k+1] = alt[k] + (Rd*T[k])/gs * np.log(Plev[k]/Plev[k+1])
        #print(k,alt[k],Plev[k])

    return alt

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

nc = Dataset('atmos_average_T_dwarf.nc')

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
print(time,nlay,nlat,nlon,ni)

it = int(t)

# Extract temperature
T = np.zeros((nlay,nlat,nlon))
T = abs(nc.variables['temp'][it,:,:,:])

# Extract VMR of species
sp = ['OH', 'H2', 'H2O', 'H', 'CO', 'CO2', 'O', 'CH4', 'C2H2', 'NH3', 'N2', 'HCN','He']
nsp = len(sp)
# Molecular weights of species
mw = [17.00734, 2.01588, 18.01528, 1.007940, 28.0101, 44.0095, 15.99940, 16.0425, 26.0373, 17.03052, 28.01340, 27.0253, 4.0026020]
mw = np.array(mw)
VMR = np.zeros((nlay,nlat,nlon,nsp))
for s in range(nsp):
    VMR[:,:,:,s] = np.maximum(abs(nc.variables[sp[s]][it,:,:,:]),1e-36)

# Rescale VMR to total 1
for i in range(nlon):
    for j in range(nlat):
        for k in range(nlay):
          for s in range(nsp):
            VMR[k,j,i,s] = VMR[k,j,i,s]/np.sum(VMR[k,j,i,:])

q_v = np.zeros((nlay,nlat,nlon))
q_0 = np.zeros((nlay,nlat,nlon))
q_1 = np.zeros((nlay,nlat,nlon))
q_2 = np.zeros((nlay,nlat,nlon))
q_v[:,:,:] = abs(nc.variables['q_v'][t,:,:,:])
q_0[:,:,:] = abs(nc.variables['q_0'][t,:,:,:])
q_1[:,:,:] = abs(nc.variables['q_1'][t,:,:,:])
q_2[:,:,:] = abs(nc.variables['q_2'][t,:,:,:])

mu = np.zeros((nlay,nlat,nlon))
mu[:,:,:] = abs(nc.variables['mu'][t,:,:,:])


# Get level pressure grid from bk and pk
Ps = np.zeros((nlat,nlon))
Ps = nc.variables['ps'][it,:,:]

Plev = np.zeros((nlev,nlat,nlon))
for j in range(nlat):
  for i in range(nlon):
    for k in range(nlev):
      Plev[k,j,i] = Ps[j,i] * bk[k] + pk[k]


#Invert T, P, VMR and q_c
for j in range(nlat):
    for i in range(nlon):
      T[:,j,i] = T[::-1,j,i] 
      Plev[:,j,i] = Plev[::-1,j,i]
      q_v[:,j,i] = q_v[::-1,j,i]
      q_0[:,j,i] = q_0[::-1,j,i]
      q_1[:,j,i] = q_1[::-1,j,i]
      q_2[:,j,i] = q_2[::-1,j,i]
      mu[:,j,i] = mu[::-1,j,i]
      for s in range(nsp):
          VMR[:,j,i,s] = VMR[::-1,j,i,s]

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
      #alt_3D[:,j,i] = calc_alt_col(nlev,Plev[:,j,i],T[:,j,i],g,Rd)
      #alt_3D[:,j,i] = calc_alt_col_var_g(nlev,Plev[:,j,i],T[:,j,i],g,Rd,Rp)
      #alt_3D[:,j,i] = calc_alt_col_p0(nlev,Plev[:,j,i],Play[:,j,i],T[:,j,i],g,Rd,Rp,p0) + Rp
      alt_3D[:,j,i] = calc_alt_col_var_g_p0(nlev,Plev[:,j,i],Play[:,j,i],T[:,j,i],g,Rd,Rp,p0) + Rp
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
iVMR = np.zeros((nlay,nlat,nlon,nsp))
imu = np.zeros((nlay,nlat,nlon))
iq_v =  np.zeros((nlay,nlat,nlon))
iq_0 = np.zeros((nlay,nlat,nlon))
iq_1 = np.zeros((nlay,nlat,nlon))
iq_2 = np.zeros((nlay,nlat,nlon))

for j in range(nlat):
    for i in range(nlon):
        iT[:,j,i] = np.interp(alt_mid[:],alt_3D_mid[:,j,i],T[:,j,i])
        iP[:,j,i] = 10.0**np.interp(alt_mid[:],alt_3D[:,j,i],np.log10(Plev[:,j,i]),right=-7)
        imu[:,j,i] = np.interp(alt_mid[:],alt_3D_mid[:,j,i],mu[:,j,i])
        iq_v[:,j,i] = 10.0**np.interp(alt_mid[:],alt_3D_mid[:,j,i],np.log10(q_v[:,j,i]),right=-30,left=-30)
        iq_0[:,j,i] = 10.0**np.interp(alt_mid[:],alt_3D_mid[:,j,i],np.log10(q_0[:,j,i]),right=-30,left=-30)
        iq_1[:,j,i] = 10.0**np.interp(alt_mid[:],alt_3D_mid[:,j,i],np.log10(q_1[:,j,i]),right=-30,left=-30)
        iq_2[:,j,i] = 10.0**np.interp(alt_mid[:],alt_3D_mid[:,j,i],np.log10(q_2[:,j,i]),right=-30,left=-30)
        for s in range(nsp):
          iVMR[:,j,i,s] = 10.0**np.interp(alt_mid[:],alt_3D_mid[:,j,i],np.log10(VMR[:,j,i,s]))
        # Now perform pressure extrapolation in height for when max alt < alt_grid in a column
        # approximate this as a hydrostatic layer
        for k in range(nlay):
          if (iP[k,j,i] == 1e-7):
            p0 = iP[k-1,j,i]
            z0 = alt_mid[k-1]
            H0 = Rd * iT[k,j,i] / g
            iP[k,j,i] = p0 * np.exp(-(alt_mid[k] - z0)/H0)
            if (iP[k,j,i] < 1e-12):
              iP[k,j,i] = 1e-12

# Calculate the cloud size distribution for gCMCRT log-normal distribution - do everything in cgs
kb = 1.380649e-16
R_gas = 8.31446261815324e7
amu = 1.66053906892e-24
rho_d = 1.99

nd_atm = np.zeros((nlay,nlat,nlon))
nd_atm[:,:,:] = (iP[:,:,:]*1e1)/(kb*iT[:,:,:])
rho = np.zeros((nlay,nlat,nlon))
rho[:,:,:] = (iP[:,:,:]*1e1*imu[:,:,:]*amu)/(kb * iT[:,:,:])

N0 = np.zeros((nlay,nlat,nlon))
N0[:,:,:] = iq_0[:,:,:] * nd_atm[:,:,:]


m_c = np.zeros((nlay,nlat,nlon))
m_c[:,:,:] = (iq_1[:,:,:]*rho[:,:,:])/(N0[:,:,:])
r_c = np.zeros((nlay,nlat,nlon))
r_c[:,:,:] = ((3.0*m_c[:,:,:])/(4.0*np.pi*rho_d))**(1.0/3.0)

nclsp = 1
nmode = 1
VMR_cl = 1.0
sig = 1.0
cl_sp = 'KCl'
fname = 'FMS.clprf'
print('Outputting cloud profile: ', fname)
f = open(fname,'w')
f.write('\n')
f.write('\n')
f.write(str(ni) + ' ' + str(nmode) + '\n')
f.write('\n')
f.write(str(nclsp) + '\n')
f.write(cl_sp + '\n')
f.write('\n')
f.write('\n')
n = 0
for j in range(nlat):
    for i in range(nlon):
        for k in range(nlay):
            f.write(str(n+1) + ' ' + str(r_c[k,j,i]*1e4) + ' ' + str(sig) + ' ' + str(N0[k,j,i]) + ' ' +  str(VMR_cl) +  '\n')
            n = n + 1
f.close()

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
  f.write(str(k+1) + ' ' + str((alt[k]) * 100.0) + '\n')
f.close()

# Output the 3D profiles in CMCRT format
head = open('../data/header.txt','r')
lines = head.readlines()

fname = 'FMS.prf'
print('Outputting main profile: ', fname)
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
n = 0
for j in range(nlat):
    for i in range(nlon):
        for k in range(nlay):
            prf.write(str(n+1) + ' ' + str(iP[k,j,i]/1e5) + ' ' + str(iT[k,j,i]) + ' ' + str(imu[k,j,i]) + ' ' + " ".join(str(l) for l in iVMR[k,j,i,:]) + '\n')
            n = n + 1
prf.close()
