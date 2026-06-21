# gCMCRT - 3D Monte Carlo Radiative Transfer for exoplanet atmospheres with GPUs

WARNING - this code is not a black box and requires some getting used to, that said it typically takes a student a few hours of tinkering to get a first spectra. Once the model works the first time, further tinkering with options/physics etc becomes more simple.
The best way to learn the code is to use it! And feel free to contact Elsie or other experienced users if stuck or confused.

! If you need the code to do more fancy things than explained here, contact an experienced user or Elsie, usually most fancy things are possible to perform or calculate.

More extensive documentation is in the making, for now here are some bare bones instructions.

Two example models are provided:

1. A 3D WASP-33b SPARC/MITgcm GCM model emission example - which extracts chemical abundances at CE from the file
2. A 3D 10x Solar WASP-39b Exo-FMS GCM model transmission example - uses chemical equilibirum abundance table interpolation
3. A 1D WASP-39b 1D VULCAN + CARMA output example - example extracting 1D chemical and cloud data
4. Several MALBEC 1D benchmark models - examples for using the code in 1D
5. Y-dwarf example with KCl cloud opacity calculation

## k-tables and CE data

k-tables, chemical equilibirum interpolation tables and Y dwarf example can be found here:
https://drive.google.com/drive/folders/1HVa9xWK_GqOqknIErcXVszzhhztACExw?usp=sharing

CIA tables can be downloaded from HITRAN:
https://hitran.org/cia/

! There will be periodic updates to the k-tables and CE interpolation tables. A more flexible option to extract data from the tables is under development.

## To compile

You will need a Nvidia GPU card on your system.

You will need to install the latest drivers from Nvidia: https://www.nvidia.com/Download/index.aspx

(optional) install the CUDA toolkit: https://developer.nvidia.com/cuda-toolkit

You will need to install the CUDA hpc sdk: https://developer.nvidia.com/hpc-sdk

# Required gCMCRT formatted files

.prf file

.hprf file

.clprf (for clouds)

.iprf (if CE interpolation required)

wavelengths.wl (central band wavelengths of calculation)

# How to operate optools

To compile cd to src_optools_V2 and enter 'make'.
To de-compile enter 'make clean'.

Compile options can be altered in the Makefile

optools uses a fortran namelist (.nml) and parameter (.par) file to communicate with the code.

## optools.par file

Is quite self-explanatory - fill in the number of species followed by the number of species (See examples)

## optools.nml file

Is more difficult to fill out correctly:

### &CK_nml - corr-k namelist

pre_mixed - Does a pre-mixed table interpolation (.True.), otherwise random overlap (.False.)

interp_wl - Interpolate to the wavelengths.wl file (.True.) (.False. if exact wavelengths of the ck table are used)

iopts - Integer option number (dev-only)

form - 1 (NEMESIS format), 2 (gCMCRT format) Note, for multiple k-tables this is a comma separated list, e.g. 2,2,2

nG - number of g-ordinances in k-table 

gmin1, gmax1, gmin2, gmax2 - g-ordinate split limits used when reading/rebinning two-part k-tables.

rebin - Rebin the k-table g-ordinates (.True.) or keep the table as supplied (.False.)

nrebin - Number of g-ordinates after rebinning.

paths - list of path to the k-table data
NOTE: THESE PATHS MUST BE IN THE SAME SPECIES ORDER AS THE SPECIES IN THE optools.par FILE !!!!

### &lbl_nml - line-by-line namelist

interp_wl - interpolate to wavelengths.wl file (.True.) or use wavelengths.wl file directly (.False.)

iopts - Integer option number (dev-only)

form - 0 (Joost's format), 1 (gCMCRT format)

paths - list of paths to the lbl data
NOTE: THESE PATHS MUST BE IN THE SAME SPECIES ORDER AS THE SPECIES IN THE optools FILE !!!

### &CIA_nml - CIA namelist

iopts - Integer option number (dev-only)

form - 0 (Special CIA species), 1 (NEMESIS CIA table), 4 (HITRAN CIA table)

paths - list of paths to the CIA data
NOTE: put a dummy path (e.g. './' ) for special species (H-, He- etc.)
NOTE: THESE PATHS MUST BE IN THE SAME SPECIES ORDER AS THE SPECIES IN THE optools FILE !!!


### &Rayleigh_nml - Rayleigh opacity namelist

iopts - Integer option number (dev-only)

### &xsec_nml - cross-section opacity namelist

iopts - Integer option number (dev-only)

form - input cross-section table format option

paths - list of paths to the cross-section data
NOTE: THESE PATHS MUST BE IN THE SAME SPECIES ORDER AS THE SPECIES IN THE optools FILE !!!

### &cl_nml - clouds namelist

iopts - Integer option number (dev-only)

imix = 1 (Bruggeman optical constant mixing), 2 (LLL method)

idist - 0 (Read bin model results), 1 (single particle size), 2 (3 size peaked near mean size), 3 (log-normal), 4 (Gamma), 5 (Inv. Gamma), 6 (Rayleigh), 7 (Hansen), 8 (Exponential)

ndist - number of size distibution points (log-spaced between amin and amax)

idist_int - 1 (Trapezium rule integration)

imie - 0 (Size limiting method), 1 (MieX), 2 (MieExt), 3 (BHMIE), 4 (DHS), 5 (BHCOAT), 6 (LX-MIE)
[Note, some of these are experimental, 0 or 6 typically recommended]

form - 5 (DIHRT format nk-tables)

paths - list of paths to the nk data
\NOTE: THESE PATHS MUST BE IN THE SAME SPECIES ORDER AS THE SPECIES IN THE optools FILE !!!

sig -  (log-normal sigma value (note, not ln(sigma) the actual sigma))

eff_fac - effective mean size variance

veff - effective variance (Hansen)

amin - Minimum distribution particle size (um)

amax - Maximum distribution particle size (um)

fmax - parameter for DHS theory

cld_tab_read - .True. = read already tabulated cloud optical properties from cld_ext.txt, cld_a.txt and cld_g.txt instead of calculating them from optical constants.


# How to operate gCMCRT

To compile cd to src_gCMCRT and enter 'make'.
To de-compile enter 'make clean'.

Compile options can be altered in the Makefile

gCMCRT uses a fortran namelist (.nml) file to communicate with the code.

## Options for gCMCRT in CMCRT.nml

###  &main

exp_name - Name of experiment; used as the prefix for profile files such as `<exp_name>.prf`, `<exp_name>.hprf` and related inputs.

xper - Required mode of gCMCRT. Active values in `gpuCMCRT.f90` are `3D_sph_tests`, `3D_sph_pol`, `3D_sph_alb`, `3D_sph_trans`, `3D_sph_em`, `3D_sph_em_hi` and `3D_sph_trans_hi`. `1D_pp`, `1D_sph` and `2D_car_gal` are listed in the dispatcher but currently commented out.

do_trans - .True. = Transmission limb sampling (for transit spectra), .False. = Normal sampling

oneD - Use a 1D profile as input

threeD - Use a 3D profile as input

do_infslab - Legacy infinite slab switch; not used by the active example experiments.

do_diffuse - Legacy diffuse/cartesian experiment switch; not used by the active example experiments.

do_cart_3D - Legacy cartesian 3D switch; not used by the active example experiments.

do_moments - experimental not in use

do_images - produce pixel images of the output (warning, can be a large file)

do_cf - calculate the contribution functions for emission or transmission (transmission is experimental)

inc_ck - .True. = Read in CK.cmcrt file (corr-k mode)

inc_lbl - .True. = Read in lbl.cmcrt file (lbl mode)

inc_CIA - .True. = Read in CIA.cmcrt file (inc. CIA opacity)

inc_Ray - .True. = Read in Rayleigh.cmcrt file (inc. Rayleigh opacity)

inc_cld - .True. = Read in cl_k.cmcrt, cl_a.cmcrt and cl_g.cmcrt files (inc. Cloud opacity and scattering properties)

inc_xsec - .True. = Read in xsec.cmcrt file (inc. xsec opacity)

do_scat_loop - .True. = Allow multiple scattering (usually .False. unless needed)

do_g_bias - .True. = Use biasing when sampling g-ordinance (typically .True.)

wght_deg - experimental do not use

do_BB_band - experimental do not use

lbl - .True. (lbl mode), .False. (corr-k mode)

ck - .True. (corr-k mode), .False. (lbl mode)

LHS - .True. (use Latin Hypercube Sampling - currently for transmission limbs only)

! Wind parameters for hi-res los velocity (lbl mode only)

doppler_on - .True. = Apply doppler shifting to local opacity

winds_on - .True. = Apply doppler shifting due to winds

rotation_on - .True. = Apply doppler shifting due to rotation

orbit_on - .True. = Apply doppler shifting due to orbital motion

orbital_period - X  orbital period (days)

systemic_velocity -  X  Systematic velocity (km s-1)

sm_ax -  X semi-major axis (AU)

ecc - Orbital eccentricity, used by the phase/eclipse geometry.

xpix -  number of pixels in x direction for pixel maps

ypix - number of pixels in y direction for pixel maps

do_Draine - .True. = Perform the Draine G value calculation

Draine_alp - Alpha value for the Draine G calculation (typically around 0.5)

do_LD - .True. = Applying the limb darkening scheme (lbl mode only)

ilimb - Integer limb darkening law selection. Options are 1 linear, 2 quadratic, 3 square-root, 4 logarithmic, 5 exponential, 6 three-parameter, 7 four-parameter non-linear and 8 power-2.

LD_c - Limb darkening coefficients

Rs - Radius of star in Solar units. Used for limb darkening and phase/eclipse geometry.

inc - Inclination in degrees. Used for limb darkening and phase/eclipse geometry.

phase - Phase value used by the limb-darkening setup.

do_phase - .True. = calculate an emission phase curve using `n_phase`; .False. = use the single `viewphi`/`viewthet` value from the experiment namelist.

do_eclipse - .True. = include eclipse/occultation geometry during phase-curve emission calculations.

n_phase - Number of phases to calculate (= 1 for single phase)

n_ecl - Number of eclipse points used when `do_phase` and `do_eclipse` are both .True.

do_surf - Surface switch for surface-enabled experiments.

T_surf - Surface temperature.

emis_surf - Surface emissivity.

alb_surf - Surface albedo.

### &sph_3D_em

Nph_tot = Number of photon packets per phase

s_wl = Start wavelength integer

n_wl = End wavelength integer

pl = 0.51 (polarisation parameter; experimental)

pc = 0.39 (polarisation parameter; experimental)

sc = 1.0 (polarisation parameter; experimental)

iscat = Scattering phase function choice (See mc_k_scatt.f90)

n_theta = Number of latitudes + 1 in GCM (i.e. number of latitude edges)

n_phi = Number of longitudes + 1 in GCM (i.e. number of latitude edges)

n_lay = Number of layers in GCM (NOTE: not levels)

viewthet = Viewing angle in latitude (typically 90 = edge on)

viewphi = Viewing angle(s) in longitude (typically 0 = dayside, 180 = nightside)

Emission only option:

xi_emb = Emission biasing (typically ~0.99)

### &sph_3D_trans

Nph = Number of photon packets per wavelength

s_wl = Start wavelength integer

n_wl = End wavelength integer

pl = 0.51 (polarisation parameter; experimental)

pc = 0.39 (polarisation parameter; experimental)

sc = 1.0 (polarisation parameter; experimental)

iscat = Scattering phase function choice (See mc_k_scatt.f90)

n_theta = Number of latitudes + 1 in GCM

n_phi = Number of longitudes + 1 in GCM

n_lay = Number of layers in GCM (NOTE: not levels)

viewthet = Viewing angle in latitude (typically 90 = edge on)

viewphi = Viewing angle in longitude (typically 180 for transmission)

nb_cf = number of layers for interpolation to calculate contribution function

### &sph_3D_alb

Used by `xper = '3D_sph_alb'`.

Nph = Number of photon packets per wavelength and phase

s_wl = Start wavelength integer

n_wl = End wavelength integer

pl = 0.51 (polarisation parameter; experimental)

pc = 0.39 (polarisation parameter; experimental)

sc = 1.0 (polarisation parameter; experimental)

iscat = Scattering phase function choice (See mc_k_scatt.f90)

n_theta = Number of latitudes + 1 in GCM

n_phi = Number of longitudes + 1 in GCM

n_lay = Number of layers in GCM (NOTE: not levels)

viewthet = Viewing angle in latitude

viewphi = Viewing angle(s) in longitude. For phase curves this should have `n_phase` values unless using the `-vphi` command-line override.

### &sph_3D_alb for 3D_sph_pol

The polarization test experiment `xper = '3D_sph_pol'` also reads a namelist called `&sph_3D_alb`, but its accepted options are different:

Nph = Number of photon packets per wavelength

n_wl = End wavelength integer

pl = 0.51 (polarisation parameter; experimental)

pc = 0.39 (polarisation parameter; experimental)

sc = 1.0 (polarisation parameter; experimental)

n_theta = Number of latitudes + 1

n_phi = Number of longitudes + 1

n_lay = Number of layers

viewthet = Viewing angle in latitude

viewphi = Viewing angle in longitude

Do not include `s_wl` or `iscat` in this block when using `xper = '3D_sph_pol'`.

### &sph_3D_em_hires

Used by `xper = '3D_sph_em_hi'`.

Nph_tot = Number of photon packets per phase

n_wl = End wavelength integer

pl = 0.51 (polarisation parameter; experimental)

pc = 0.39 (polarisation parameter; experimental)

sc = 1.0 (polarisation parameter; experimental)

iscat = Scattering phase function choice (See mc_k_scatt.f90)

n_theta = Number of latitudes + 1 in GCM

n_phi = Number of longitudes + 1 in GCM

n_lay = Number of layers in GCM

viewthet = Viewing angle in latitude

viewphi = Viewing angle(s) in longitude; array length should match `n_phase`.

xi_emb = Emission biasing (typically ~0.99)

### &sph_3D_trans_hires

Used by `xper = '3D_sph_trans_hi'`.

Nph = Number of photon packets per wavelength and phase

n_wl = End wavelength integer

pl = 0.51 (polarisation parameter; experimental)

pc = 0.39 (polarisation parameter; experimental)

sc = 1.0 (polarisation parameter; experimental)

iscat = Scattering phase function choice (See mc_k_scatt.f90)

n_theta = Number of latitudes + 1 in GCM

n_phi = Number of longitudes + 1 in GCM

n_lay = Number of layers in GCM

viewthet = Viewing angle in latitude

viewphi = Viewing angle(s) in longitude; array length should match `n_phase`.

### &diffuse_nml

Legacy cartesian diffuse experiment namelist. The dispatcher currently lists `2D_car_gal`, but that call is commented out.

Nph_tot = Total number of photon packets

kappa = Opacity

albedo = Single-scattering albedo

hgg = Henyey-Greenstein asymmetry parameter

pl = 0.51 (polarisation parameter; experimental)

pc = 0.39 (polarisation parameter; experimental)

sc = 1.0 (polarisation parameter; experimental)

xmax, ymax, zmax = Cartesian domain extents

rimage = Image radius

viewthet = Viewing angle in latitude

viewphi = Viewing angle in longitude

nxim, nyim = Image pixel dimensions

nxg, nyg, nzg = Cartesian grid dimensions

# What is output and how do I make synthetic observations?

# Example tutorials

## WASP-33b GCM dayside and nightside emission spectrum (with details on phase curve operation)

Example using a SPARC/MITgcm WASP-33b model. In this example we produce a dayside and nightside spectra. 
1. Use the extract script to extract the GCM data into the gCMCRT .hprf and .prf format, chemical abundances are extracted alongside.
2. Run goptools to produce the corr-k, CIA and Rayleigh opacity files.
3. Run gCMCRT to produce Em_001.txt (dayside) and Em_002.txt (nightside) output files
4. Run em_spec.py to convert the output to synthetic observations, Fp/Fs, Fp and Tb. This file contains useful information on how to produce emission spectra.

This example can be extended easily to produce phase curves.

## WASP-39b GCM transmision spectrum (with details on using CE interpolation tables)

Example using an Exo-FMS WASP-39b model. In this example we produce a transmission spectrum.
1. Use the extract script to extract the GCM data into the gCMCRT .hprf and .iprf format
2. Use interp_iprf.py to interpolate the CE abundances to the T,p of the GCM and produce the .prf file.
3. Run goptools to produce the corr-k, CIA and Rayleigh opacity files.
4. Run gCMCRT to produce the Transmission.txt file 
5. Run trans_spec.py to convert the output to synthetic observations, Rp/Rs and compare to G395H and SOSS observations. (This file contains useful information on transmission spectrum fitting etc)

## WASP-39b 1D VULCAN + CARMA model (with details on using the code in 1D)

## MALBEC benchmarks (further examples on 1D code use and comparing to other codes)
