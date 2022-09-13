import numpy as np
import matplotlib.pylab as plt

kb = 1.380649e-23
amu = 1.66053906660e-27
R = 8.31446261815324

# Script to convert V. Parmentier format SPARC/MITgcm output
# You will need the 'delr.txt' file which denotes the pressure levels

# Filename to output
fname = 'PTprofiles-with-winds-and-abundances-WASP-33b-TiO-fix-3-NEW-OPA_nit_518400av_more_species_2_C16_res.dat'

# Manually add gravity, Rd and Radius (now in cgs)
grav = 20.25 * 100.0
Rd = 3714.3 * 1e4
Rp = 1.16e8 * 100.0
Rspec = Rd / 1e4 # mks Rd

# Have to manually add the number of layers, latitude points and longitude points
nlay = 53
nlev = nlay + 1
nlat = 30
nlon = 64

# Have to manually give the bottom pressure boundary (level) value (reference pressure)
p0 = 200.0

# Species in the GCM output file
sp = ['H2','H-','H2O','CO','TiO','VO','e-','Fe','FeH','CH4','Na','K','NH3','CO2','HCN','H','e-','H+','H3+','He']
nsp = len(sp)

# Hypsometric equation to calculate height grid
def calc_alt_col(nlev,p,T,grav,Rd):
    alt = np.zeros(nlev)
    alt[0] = 0.0
    for l in range(nlay):
        alt[l+1] = alt[l] + Rd/grav * T[l] * np.log(p[l]/p[l+1])
        #print(l,alt[l],T[l],p[l],p[l+1])
    return alt

data = np.loadtxt(fname,skiprows=6,comments=' ',delimiter=',')

ni = nlay * nlat * nlon

lons = np.unique(data[:,1])
lats = np.unique(data[:,2])
P = np.unique(data[:,3])

print(nlon,len(lons),lons)
print(nlat,len(lats),lats)
print(nlay,len(P),P)

# Find the pressure at the level interfaces
P[:] = P[::-1]
fname = 'delr.txt'
dP = np.loadtxt(fname)

print(dP)
print(len(dP))

Plev = np.zeros(nlev)
Plev[0] = p0
for i in range(nlay):
    Plev[i+1] = Plev[i] - dP[i]/1e5  

print(Plev)
print(len(Plev))

print(len(data[:,0]),nlay*nlat*nlon)


# Roll the longitude array
rl = int(len(lons)/2)
print(rl)
lons = np.roll(lons,rl)
lons = np.where(lons > 0.0, lons, 360 + lons)
print(lons)

# Save T and VMR
T = np.zeros((nlon,nlat,nlay))
u = np.zeros((nlon,nlat,nlay))
v = np.zeros((nlon,nlat,nlay))
w =  np.zeros((nlon,nlat,nlay))
rho = np.zeros((nlon,nlat,nlay))
mu = np.zeros((nlon,nlat,nlay))
VMR = np.zeros((nlon,nlat,nlay,nsp))
P_3D = np.zeros((nlon,nlat,nlay))

P_3D[:,:,:] = data[:,3].reshape((nlon,nlat,nlay))
T[:,:,:] = data[:,4].reshape((nlon,nlat,nlay))
u[:,:,:] = data[:,5].reshape((nlon,nlat,nlay))
v[:,:,:] = data[:,6].reshape((nlon,nlat,nlay))
w[:,:,:] = data[:,7].reshape((nlon,nlat,nlay))
rho[:,:,:] = data[:,8].reshape((nlon,nlat,nlay))
for i in range(nsp):
    VMR[:,:,:,i] = data[:,9+i].reshape((nlon,nlat,nlay))

#Roll all the arrays
for k in range(nlat):
  for i in range(nlay):
      T[:,k,i] = np.roll(T[:,k,i],rl)
      u[:,k,i] = np.roll(u[:,k,i],rl)
      v[:,k,i] = np.roll(v[:,k,i],rl)
      w[:,k,i] = np.roll(w[:,k,i],rl)
      rho[:,k,i] = np.roll(rho[:,k,i],rl)
      for l in range(nsp):
          VMR[:,k,i,l] = np.roll(VMR[:,k,i,l],rl)

mu[:,:,:] = R/Rspec * 1000.0


# Reverse the i index - convert to cgs units
for j in range(nlon):
  for k in range(nlat):
      T[j,k,:] = T[j,k,::-1]
      u[j,k,:] = u[j,k,::-1] * 100.0
      v[j,k,:] = v[j,k,::-1] * 100.0
      w[j,k,:] = w[j,k,::-1] * 100.0
      #rho[j,k,:] = rho[j,k,::-1]
      for l in range(nsp):
          VMR[j,k,:,l] = 10.0**VMR[j,k,::-1,l]


#find altitude and midpoint values
alt = np.zeros((nlon,nlat,nlev))
alt_mid = np.zeros((nlon,nlat,nlay))
for j in range(nlon):
  for k in range(nlat):
      alt[j,k,:] = calc_alt_col(nlev,Plev[:],T[j,k,:],grav,Rd)
      alt_mid[j,k,:] = (alt[j,k,0:nlay] + alt[j,k,1:nlev])/2.0

imax = np.unravel_index(np.argmax(alt, axis=None), alt.shape)
print(imax)
alt_grid = alt[imax[0],imax[1],:]
alt_grid_mid = (alt_grid[0:nlay] + alt_grid[1:nlev])/2.0

print(alt_grid)
print(alt_grid + Rp) 


# Interpolate p, T VMR and windws to alt_grid midpoints
iT = np.zeros((nlon,nlat,nlay))
iP = np.zeros((nlon,nlat,nlay))
iu = np.zeros((nlon,nlat,nlay))
iv = np.zeros((nlon,nlat,nlay))
iw = np.zeros((nlon,nlat,nlay))
iVMR = np.zeros((nlon,nlat,nlay,nsp))
for j in range(nlon):
  for k in range(nlat):
      iT[j,k,:] = np.interp(alt_grid_mid[:],alt_mid[j,k,:],T[j,k,:])
      iP[j,k,:] = 10.0**np.interp(alt_grid_mid[:],alt_mid[j,k,:],np.log10(P[:]),right=-12)
      iu[j,k,:] = np.interp(alt_grid_mid[:],alt_mid[j,k,:],u[j,k,:])
      iv[j,k,:] = np.interp(alt_grid_mid[:],alt_mid[j,k,:],v[j,k,:])
      iw[j,k,:] = np.interp(alt_grid_mid[:],alt_mid[j,k,:],w[j,k,:])
      for l in range(nsp):
          iVMR[j,k,:,l] = 10.0**np.interp(alt_grid_mid[:],alt_mid[j,k,:],np.log10(VMR[j,k,:,l]))
      # Now perform extrapolation in height for max alt < alt_grid
      for i in range(nlay):
          if (iP[j,k,i] == 1e-12):
            p0 = iP[j,k,i-1]
            z0 = alt_grid_mid[i-1]
            H0 = Rd * iT[j,k,i] / grav
            iP[j,k,i] = p0 * np.exp(-(alt_grid_mid[i] - z0)/H0)
            if (iP[j,k,i] < 1e-12):
              iP[j,k,i] = 1e-12
          #print(l,iP[l,y,x]/1e5,iT[l,y,x],alt_grid_mid[l])


# Output the 3D profiles in CMCRT format
head = open('../data/header.txt','r')
lines = head.readlines()

fname = 'SPARC-MITgcm.prf'
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
for k in range(nlat):
    for j in range(nlon):
        for i in range(nlay):
            prf.write(str(n+1) + ' ' + str(iP[j,k,i]) + ' ' + str(iT[j,k,i]) + ' ' + str(mu[j,k,i]) + ' ' + " ".join(str(l) for l in iVMR[j,k,i,:]) + '\n')
            n = n + 1
prf.close()

# Now output T-p profile in bar and K to an interpolation file (iprf) - after which we can interpolate to GGChem values to get VMRs
fname = 'SPARC-MITgcm.iprf'
print('Outputting interpolatable T-p profile: ', fname)
f = open(fname,'w')
n = 0
for k in range(nlat):
    for j in range(nlon):
        for i in range(nlay):
            f.write(str(n+1) + ' ' +  str(iP[j,k,i]) + ' ' + str(iT[j,k,i]) + '\n')
            n = n + 1
f.close()


fname = 'SPARC-MITgcm.hprf'
print('Outputting height profile: ', fname)
f = open(fname,'w')
for i in range(nlev):
    f.write(str(i+1) + ' ' + str(alt_grid[i] + Rp) + '\n')

f.close()


# Output the 3D profiles in CMCRT format
fname = 'SPARC-MITgcm.wprf'
print('Outputting wind profile profile: ', fname)
f = open(fname,'w')
n = 0
for k in range(nlat):
    for j in range(nlon):
        for i in range(nlay):
            f.write(str(n+1) + ' ' +  str(iu[j,k,i]) + ' ' + str(iv[j,k,i]) + ' ' + str(iw[j,k,i])  + '\n')
            n = n + 1
f.close()



