import numpy as np
import matplotlib.pylab as plt


Rd = 3523.0 
Rp = 9.65e7 
grav = 10.0


nlat = 48
nlon = 96
nlay = 50
nlev = nlay + 1
nlines = nlat * nlon * nlay

#fname = 'cloudreport_T2250g10_10pct_nuclim.txt'
fname = 'cloudreport_T2250g10_nuc_cor_mass.txt'
f = open(fname,'r')
for l in range(5):
  line = f.readline()

# KCl, ZnS, Na2S, MnS, Cr2O3, SiO2, Mg2SiO4, VO, Ni, Fe, Ca2SiO4, CaTiO3, Al2O3

cld_sp = ['KCl', 'ZnS', 'Na2S', 'MnS', 'Cr2O3', 'SiO2', 'Mg2SiO4', 'VO', 'Ni', 'Fe', 'Ca2SiO4', 'CaTiO3', 'Al2O3']

rhod = [1.98e3,409e3,1.86e3,4.0e3,5.22e3,2.65e3,3.27e3,5.76e3,8.9e3,7.9e3,3.34e3,3.98e3,3.95e3]

fname = 'radius.txt'
rd = np.loadtxt(fname, delimiter=',')
rd[:] = rd[::-1]

# Make containers for data
# lat , lon , level , altitude(m), pressure(bars) , temp(k) , temp_stdev(k), EW vel(m/s) , NS vel , vert vel ,
# Cloud 1-13 vis tau (650nm), cloud total VIS tau, cloud 1-13 IR tau (5um), cloud total IR tau

lats = np.zeros(nlines)
lons =  np.zeros(nlines)
lays = np.zeros(nlines)
alt_1D =  np.zeros(nlines)
P_1D =  np.zeros(nlines)
T_1D =  np.zeros(nlines)
T_st_1D =  np.zeros(nlines)
u_1D =  np.zeros(nlines)
v_1D =  np.zeros(nlines)
w_1D =  np.zeros(nlines)
cld_mass_1D = np.zeros((13,nlines))
cld_mass_tot_1D = np.zeros(nlines)
cld_IR_1D = np.zeros((13,nlines))
cld_IR_tot_1D = np.zeros(nlines)


for l in range(nlines):
  line1 = f.readline().split()
  line2 = f.readline().split()
  lats[l] = float(line1[0])
  lons[l] = float(line1[1])
  lays[l] = int(line1[2])
  alt_1D[l] = float(line1[3])
  P_1D[l] = float(line1[4])
  T_1D[l] = float(line1[5])
  T_st_1D[l] = float(line1[6])
  u_1D[l] = float(line1[7])
  v_1D[l] = float(line1[8])
  w_1D[l] = float(line1[9])
  for c in range(13):
      cld_mass_1D[c,l] = float(line1[10+c])
  cld_mass_tot_1D[l] = float(line1[23])
  for c in range(11):
      cld_IR_1D[c,l] = float(line1[24+c])
  cld_IR_1D[-2,l] = float(line2[0])
  cld_IR_1D[-1,l] = float(line2[1])
  cld_IR_tot_1D[l] = float(line2[1])
  #print(P_1D[l],cld_mass_tot_1D[l],cld_IR_tot_1D[l])

## First deal with generating the iprf
lons = np.unique(lons)
lats = np.unique(lats)

print(nlon,len(lons),lons)
print(nlat,len(lats),lats)

## Convert 1D arrays to 3D

P = np.zeros((nlat,nlon,nlay))
T = np.zeros((nlat,nlon,nlay))
u = np.zeros((nlat,nlon,nlay))
v = np.zeros((nlat,nlon,nlay))
w =  np.zeros((nlat,nlon,nlay))
alt_mid = np.zeros((nlat,nlon,nlay))
cld_mass = np.zeros((13,nlat,nlon,nlay))
cld_IR = np.zeros((13,nlat,nlon,nlay))

P[:,:,:] = P_1D.reshape((nlat,nlon,nlay))
T[:,:,:] = T_1D.reshape((nlat,nlon,nlay))
u[:,:,:] = u_1D.reshape((nlat,nlon,nlay))
v[:,:,:] = v_1D.reshape((nlat,nlon,nlay))
w[:,:,:] = w_1D.reshape((nlat,nlon,nlay))
alt_mid[:,:,:] = alt_1D.reshape((nlat,nlon,nlay))
for i in range(13):
  cld_mass[i,:,:,:] = cld_mass_1D[i,:].reshape((nlat,nlon,nlay))
  cld_IR[i,:,:,:] = cld_IR_1D[i,:].reshape((nlat,nlon,nlay))

# Get the delta cloud IR
cld_IR_lev = np.zeros((13,nlat,nlon,nlev))
dcld_IR = np.zeros((13,nlat,nlon,nlay))
for j in range(nlat):
  for k in range(nlon):
    for z in range(nlay-1):
      for i in range(13):
        dcld_IR[i,j,k,z] = abs(cld_IR[i,j,k,z+1] - cld_IR[i,j,k,z])

# Reverse the z axis
for j in range(nlat):
  for k in range(nlon):
      P[j,k,:] = P[j,k,::-1]
      T[j,k,:] = T[j,k,::-1]
      u[j,k,:] = u[j,k,::-1] 
      v[j,k,:] = v[j,k,::-1] 
      w[j,k,:] = w[j,k,::-1] 
      alt_mid[j,k,:] = alt_mid[j,k,::-1]
      for i in range(13):
          cld_mass[i,j,k,:] = cld_mass[i,j,k,::-1]
          dcld_IR[i,j,k,:] = dcld_IR[i,j,k,::-1]

# Plot intermediate results 

fig = plt.figure()
ilat = int(nlat/2)
ilon = int(nlon/2)

#plt.plot(T[ilat,0,:],P[ilat,0,:])
#plt.plot(T[ilat,ilon,:],P[ilat,ilon,:])

for i in range(13):
  plt.plot(dcld_IR[i,ilat,ilon,:],P[ilat,ilon,:])
  #plt.plot(cld_IR[i,ilat,ilon,:],P[ilat,ilon,:])
  #plt.plot(cld_mass[i,ilat,ilon,:],P[ilat,ilon,:])

plt.yscale('log')
plt.xscale('log')

plt.show()

#quit()

# Find maximum altitude profile
imax = np.unravel_index(np.argmax(alt_mid, axis=None), alt_mid.shape)
alt_grid_mid = alt_mid[imax[0],imax[1],:]

# Find and approximate the altitude levels
alt_grid = np.zeros(nlev)
alt_grid[0] = 0.0
for i in range(nlay-1):
    alt_grid[i+1] = (alt_grid_mid[i+1] + alt_grid_mid[i])/2.0
# Approximate uppermost level as midpoint layer altitude
alt_grid[-1] = alt_grid_mid[-1]

# Now output Height grid in cm
fname = 'Roman.hprf'
print('Outputting height profile: ', fname)
f = open(fname,'w')
for k in range(nlev):
  f.write(str(k+1) + ' ' + str((Rp + alt_grid[k]) * 100.0) + '\n')

# Find number density and VMR of particles
cld_nd = np.zeros((nlat,nlon,nlay))
cld_VMR = np.zeros((13,nlat,nlon,nlay))
for j in range(nlat):
  for k in range(nlon):
    for i in range(nlay):
      for s in range(13):
        cld_nd[j,k,i] = cld_nd[j,k,i] + dcld_IR[s,j,k,i]/(np.pi*rd[i]**2*(alt_grid[i+1] - alt_grid[i]))
        cld_VMR[s,j,k,i] = np.maximum(cld_mass[s,j,k,i]/np.sum(cld_mass[:,j,k,i]),1e-30)
        if (np.isnan(cld_VMR[s,j,k,i]) == True or np.isinf(cld_VMR[s,j,k,i]) == True):
          cld_VMR[s,j,k,i] = 1e-30
        cld_nd[j,k,i] = np.maximum(cld_nd[j,k,i],1e-30)

# Interpolate variables to max profile altitude
# Interpolate p, T VMR and windws to alt_grid midpoints
iT = np.zeros((nlat,nlon,nlay))
iP = np.zeros((nlat,nlon,nlay))
iu = np.zeros((nlat,nlon,nlay))
iv = np.zeros((nlat,nlon,nlay))
iw = np.zeros((nlat,nlon,nlay))
icld_nd = np.zeros((nlat,nlon,nlay))
icld_VMR = np.zeros((13,nlat,nlon,nlay))
for j in range(nlat):
  for k in range(nlon):
      iT[j,k,:] = np.interp(alt_grid_mid[:],alt_mid[j,k,:],T[j,k,:])
      iP[j,k,:] = 10.0**np.interp(alt_grid_mid[:],alt_mid[j,k,:],np.log10(P[j,k,:]),right=-12)
      iu[j,k,:] = np.interp(alt_grid_mid[:],alt_mid[j,k,:],u[j,k,:])
      iv[j,k,:] = np.interp(alt_grid_mid[:],alt_mid[j,k,:],v[j,k,:])
      iw[j,k,:] = np.interp(alt_grid_mid[:],alt_mid[j,k,:],w[j,k,:])
      icld_nd[j,k,:] = np.interp(alt_grid_mid[:],alt_mid[j,k,:],cld_nd[j,k,:],right=1e-30,left=1e-30)
      for s in range(13):
          icld_VMR[s,j,k,:] = np.interp(alt_grid_mid[:],alt_mid[j,k,:],cld_VMR[s,j,k,:],right=1e-30,left=1e-30)
      # Now perform extrapolation in height for max alt < alt_grid
      for i in range(nlay):
          if (iP[j,k,i] == 1e-12):
            p0 = iP[j,k,i-1]
            z0 = alt_grid_mid[i-1]
            H0 = Rd * iT[j,k,i] / grav
            iP[j,k,i] = p0 * np.exp(-(alt_grid_mid[i] - z0)/H0)
            if (iP[j,k,i] < 1e-12):
              iP[j,k,i] = 1e-12

# Now output T-p profile in bar and K to an interpolation file (iprf) - after which we can interpolate to GGChem values to get VMRs
fname = 'Roman.iprf'
print('Outputting interpolatable T-p profile: ', fname)
f = open(fname,'w')
n = 0
for j in range(nlat):
    for i in range(nlon):
        for k in range(nlay):
            f.write(str(n+1) + ' ' +  str(iP[j,i,k]) + ' ' + str(iT[j,i,k]) + '\n')
            n = n + 1
f.close()

#Output the 3D wind profiles in CMCRT layout -> (lat,lon,lay) in units cm s-1
fname = 'Roman.wprf'
print('Outputting wind profile: ', fname)
f = open(fname,'w')
n = 0
for j in range(nlat):
    for i in range(nlon):
        for k in range(nlay):
            f.write(str(n+1) + ' ' +  str(iu[j,i,k] * 100.0) + ' ' + str(iv[j,i,k] * 100.0) + ' ' + str(iw[j,i,k] * 100.0)  + '\n')
            n = n + 1
f.close()


#Output the 3D cloud profiles
fname = 'Roman.clprf'
print('Outputting cloud profile: ', fname)
f = open(fname,'w')
f.write('\n')
f.write('\n')
f.write(str(nlines) + ' ' + '1' + '\n')
f.write('\n')
f.write('13' + '\n')
for s in range(13):
  f.write(cld_sp[s] + '\n')
f.write('\n')
f.write('\n')
n = 0
for j in range(nlat):
    for i in range(nlon):
        for k in range(nlay):
            f.write(str(n+1) + ' ' + str(rd[k]*1e6) + ' ' + str(icld_nd[j,i,k]/1e6) + ' ' + " ".join(str(l) for l in icld_VMR[:,j,i,k]) + '\n')
            n = n + 1
f.close()



