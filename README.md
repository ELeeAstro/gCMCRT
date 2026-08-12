# gCMCRT - 3D Monte Carlo Radiative Transfer for exoplanet atmospheres with GPUs

WARNING - this code is not a black box and requires some getting used to, that said it typically takes a student a few hours of tinkering to get a first spectra. Once the model works the first time, further tinkering with options/physics etc becomes more simple.
The best way to learn the code is to use it! And feel free to contact Elsie or other experienced users if stuck or confused.

! If you need the code to do more fancy things than explained here, contact an experienced user or Elsie, usually most fancy things are possible to perform or calculate.

More extensive documentation is in the making, for now here are some bare bones instructions.

The repository provides the following example models and benchmarks:

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

xper - Required mode of gCMCRT. Active values in `gpuCMCRT.f90` are `3D_sph_alb`, `3D_sph_trans`, `3D_sph_trans_lc`, `3D_sph_em`, `3D_sph_em_hi` and `3D_sph_trans_hi`.

The selected experiment determines its source sampling automatically.
`3D_sph_trans` and `3D_sph_trans_hi` use atmospheric-annulus sampling.
`3D_sph_trans_lc` uses its dedicated full-disk/overlap-lens sampler.

oneD - Use a 1D profile as input

threeD - Use a 3D profile as input

Exactly one of `oneD` and `threeD` must be `.True.`.

do_images - produce pixel images of the output (warning, can be a large file)

do_cf - calculate the contribution functions for emission or transmission (transmission is experimental)

inc_ck - .True. = Read in CK.cmcrt file (corr-k mode)

inc_lbl - .True. = Read in lbl.cmcrt file (lbl mode)

inc_CIA - .True. = Read in CIA.cmcrt file (inc. CIA opacity)

inc_Ray - .True. = Read in Rayleigh.cmcrt file (inc. Rayleigh opacity)

inc_cld - .True. = Read in cl_k.cmcrt, cl_a.cmcrt and cl_g.cmcrt files (inc. Cloud opacity and scattering properties)

inc_xsec - .True. = Read in xsec.cmcrt file (inc. xsec opacity)

inc_refrac - `.True.` reads `refrac.cmcrt`, containing gas refractivity
`nu = n - 1` for each wavelength and atmospheric cell. This input is independent
of `inc_Ray`. `mc_k_refraction` supplies the local spherical-interface Snell and
total-internal-reflection physics; the forward ray tracer does not yet invoke it.

do_scat_loop - .True. = Allow multiple scattering (usually .False. unless needed)

do_g_bias - .True. = Use biasing when sampling g-ordinance (typically .True.)

lbl - .True. (lbl mode), .False. (corr-k mode)

ck - .True. (corr-k mode), .False. (lbl mode)

Exactly one of `ck` and `lbl` must be `.True.`. `inc_ck` is only valid in `ck` mode, and `inc_lbl` is only valid in `lbl` mode. Both gas-inclusion flags may be `.False.` for a scattering/cloud-only calculation.

LHS - `.True.` enables Latin-hypercube projected-position sampling in the albedo and transmission source modes. Storage is allocated once and samples are regenerated for each requested phase/wavelength.

random_seed - Non-negative base seed used by both the CPU/LHS and GPU random-number generators. GPU packet `id - 1` selects a distinct CURAND sequence under this common seed. Reusing a value reproduces the same Monte Carlo realization; change it for an independent run. The default is `12345`.

! Wind parameters for hi-res los velocity (lbl mode only)

doppler_on - `.True.` applies Doppler shifting to local opacity. It is accepted only by `3D_sph_em_hi` and `3D_sph_trans_hi`; all other experiment modes reject it.

winds_on - .True. = Apply doppler shifting due to winds

rotation_on - .True. = Apply doppler shifting due to rotation

orbit_on - .True. = Apply doppler shifting due to orbital motion

orbital_period - X  orbital period (days)

systemic_velocity - Systemic velocity (km s-1). This observer-frame shift is applied independently of `orbit_on`.

sm_ax -  X semi-major axis (AU)

ecc - Orbital eccentricity. Only circular orbits are currently supported, so this must be 0.

xpix -  number of pixels in x direction for pixel maps

ypix - number of pixels in y direction for pixel maps

Draine_alp - Alpha parameter for the Draine phase function (typically around 0.5). When an experiment uses `iscat = 6`, the corresponding Draine `Dg` field is calculated automatically from `Draine_alp` and the cloud asymmetry parameter.

do_LD - `.True.` applies the selected stellar limb-darkening law to transit source sampling and transit-light-curve normalisation. It is independent of the opacity mode and is not valid for albedo illumination.

ilimb - Integer limb darkening law selection. Options are 1 linear, 2 quadratic, 3 square-root, 4 logarithmic, 5 exponential, 6 three-parameter, 7 four-parameter non-linear and 8 power-2.

LD_c - Limb darkening coefficients

Rs - Radius of star in Solar units. Used for limb darkening and phase/eclipse geometry.

inc - Inclination in degrees. Used for limb darkening and phase/eclipse geometry.

do_phase - .True. = calculate an emission phase curve using `n_phase`; .False. = use the single `viewphi`/`viewthet` value from the experiment namelist.

do_eclipse - .True. = include secondary eclipse/occultation geometry during phase-curve emission calculations. Partial occultations are masked in the per-phase image plane.

n_phase - Number of phases to calculate (= 1 for single phase)

n_ecl - Number of eclipse points used when `do_phase` and `do_eclipse` are both .True.

For all spherical 3D experiments, `n_theta` is the number of theta
boundaries and must be odd and at least 3, which places a boundary exactly at
the equator. `n_phi` is the number of periodic longitude boundaries and must
be at least 4.

The standard and high-resolution transmission, transit-light-curve, albedo,
standard emission and high-resolution emission experiment namelists accept
`use_block_accum`.
It defaults to `.True.` when omitted. In this mode, integrated transport
quantities are accumulated into separate CUDA-block slots and reduced on the
host after the kernel completes, reducing contention on single global scalar
atomics. Set it to `.False.` to restore the original scalar-atomic path for
A/B testing. The two paths should be statistically equivalent, but floating-
point summation order means results are not expected to be bit-for-bit
identical. Image pixels, contribution functions and failure diagnostics retain
their existing accumulation paths.

### &sph_3D_em

Nph_tot = Number of photon packets per phase

s_wl = Start wavelength integer

n_wl = End wavelength integer

pl = 0.51 (polarisation parameter; experimental)

pc = 0.39 (polarisation parameter; experimental)

sc = 1.0 (polarisation parameter; experimental)

iscat = Volumetric cloud scattering phase function: 1 isotropic, 3 Rayleigh, 4 Henyey-Greenstein, 5 two-term Henyey-Greenstein, or 6 Draine. Value 2 is a Lambertian surface law and is rejected for atmospheric volume scattering.

n_theta = Number of latitudes + 1 in GCM (i.e. number of latitude edges)

n_phi = Number of longitudes + 1 in GCM (i.e. number of longitude edges)

n_lay = Number of layers in GCM (NOTE: not levels)

viewthet = Viewing angle in latitude (typically 90 = edge on)

viewphi = Viewing angle(s) in longitude (typically 0 = dayside, 180 = nightside)

Emission only option:

xi_emb = Emission biasing (typically ~0.99)

use_block_accum = `.True.` uses per-CUDA-block accumulation for integrated flux, Stokes values, occulted flux and scattering counts; `.False.` selects the legacy scalar atomics.

### &sph_3D_trans

Nph = Number of photon packets per wavelength

s_wl = Start wavelength integer

n_wl = End wavelength integer

pl = 0.51 (polarisation parameter; experimental)

pc = 0.39 (polarisation parameter; experimental)

sc = 1.0 (polarisation parameter; experimental)

iscat = Volumetric cloud scattering phase function: 1, 3, 4, 5, or 6 as defined above.

n_theta = Number of latitudes + 1 in GCM

n_phi = Number of longitudes + 1 in GCM

n_lay = Number of layers in GCM (NOTE: not levels)

viewthet = Viewing angle in latitude (typically 90 = edge on)

viewphi = Viewing angle in longitude (typically 180 for transmission)

nb_cf = number of radial bins used for the transmission contribution function

use_block_accum = `.True.` uses per-CUDA-block accumulation for the transmission estimators, integrated flux/Stokes values and scattering counts; `.False.` selects the legacy scalar atomics.

### &sph_3D_alb

Used by `xper = '3D_sph_alb'`.

Albedo requires `do_scat_loop = .True.` because the signal is produced by scattered packets. The `viewphi` array must contain `n_phase` values.

Nph = Number of photon packets per wavelength and phase

s_wl = Start wavelength integer

n_wl = End wavelength integer

pl = 0.51 (polarisation parameter; experimental)

pc = 0.39 (polarisation parameter; experimental)

sc = 1.0 (polarisation parameter; experimental)

iscat = Volumetric cloud scattering phase function: 1, 3, 4, 5, or 6 as defined above.

n_theta = Number of latitudes + 1 in GCM

n_phi = Number of longitudes + 1 in GCM

n_lay = Number of layers in GCM (NOTE: not levels)

viewthet = Viewing angle in latitude

viewphi = Viewing angle(s) in longitude. For phase curves this should have `n_phase` values unless using the `-vphi` command-line override.

use_block_accum = `.True.` uses per-CUDA-block accumulation for integrated flux/Stokes values and scattering counts; `.False.` selects the legacy scalar atomics.

### &sph_3D_em_hires

Used by `xper = '3D_sph_em_hi'`.

This mode requires `lbl = .True.` and `ck = .False.`. If `doppler_on = .True.`, `inc_lbl` must also be `.True.` and `inc_xsec` must be `.False.`.

Nph_tot = Number of photon packets per phase

n_wl = End wavelength integer

pl = 0.51 (polarisation parameter; experimental)

pc = 0.39 (polarisation parameter; experimental)

sc = 1.0 (polarisation parameter; experimental)

iscat = Volumetric cloud scattering phase function: 1, 3, 4, 5, or 6 as defined above.

n_theta = Number of latitudes + 1 in GCM

n_phi = Number of longitudes + 1 in GCM

n_lay = Number of layers in GCM

viewthet = Viewing angle in latitude

viewphi = Viewing angle(s) in longitude; array length should match `n_phase`.

xi_emb = Emission biasing (typically ~0.99)

use_block_accum = `.True.` uses per-CUDA-block accumulation for integrated flux, Stokes values, occulted flux and scattering counts; `.False.` selects the legacy scalar atomics.

### &sph_3D_trans_hires

Used by `xper = '3D_sph_trans_hi'`.

This mode requires `lbl = .True.` and `ck = .False.`. If `doppler_on = .True.`, `inc_lbl` must also be `.True.` and `inc_xsec` must be `.False.`.

Nph = Number of photon packets per wavelength and phase

n_wl = End wavelength integer

pl = 0.51 (polarisation parameter; experimental)

pc = 0.39 (polarisation parameter; experimental)

sc = 1.0 (polarisation parameter; experimental)

iscat = Volumetric cloud scattering phase function: 1, 3, 4, 5, or 6 as defined above.

n_theta = Number of latitudes + 1 in GCM

n_phi = Number of longitudes + 1 in GCM

n_lay = Number of layers in GCM

viewthet = Viewing angle in latitude

viewphi = Viewing angle(s) in longitude; array length should match `n_phase`.

use_block_accum = `.True.` uses per-CUDA-block accumulation for the transmission estimator, integrated flux/Stokes values and scattering counts; `.False.` selects the legacy scalar atomics.

### &sph_3D_trans_lc

Used by `xper = '3D_sph_trans_lc'`. This is the phase-resolved primary-transit light-curve mode. It launches direct stellar packets through the phase-dependent projected stellar disk and ignores multiple-scattering loops.

The transit geometry is derived self-consistently from the system parameters in `&main`. For a synchronously rotating, zero-obliquity planet the GCM grid is not rotated; instead, the observer/image basis changes with phase. The projected stellar centre is stored in that same per-phase image basis before ray rejection, so ingress/egress samples the correct backlit limb sector.

Nph = Number of photon packets per wavelength and transit phase

s_wl = Start wavelength integer

n_wl = End wavelength integer

n_trans = Number of sampled transit phases across the primary-transit window

n_ingress = Number of samples assigned to each ingress/egress band when the transit has a full-transit section

pl = 0.51 (polarisation parameter; experimental)

pc = 0.39 (polarisation parameter; experimental)

sc = 1.0 (polarisation parameter; experimental)

iscat = Scattering phase-function choice retained for opacity setup consistency; this direct-only mode does not sample scattering angles.

n_theta = Number of latitudes + 1 in GCM

n_phi = Number of longitudes + 1 in GCM

n_lay = Number of layers in GCM

use_block_accum = `.True.` uses per-CUDA-block accumulation for the total, east-limb and west-limb atmospheric transmission estimators; `.False.` selects the legacy scalar atomics.

This mode derives `viewphi`, `viewthet`, projected star position and transit state from `Rs`, `inc`, `sm_ax`, `orbital_period`, `do_LD`, `ilimb` and `LD_c` in `&main`.

The code writes one file per sampled transit phase:

`Transit_001.txt`, `Transit_002.txt`, ...

The first row is:

`n_wl_out H_base H_top viewphi viewthet phase trans_state xstar ystar zstar`

where `n_wl_out = n_wl - s_wl + 1` is the number of wavelength rows in the file, `H_base = H(1)` is the opaque/reference radius, `H_top = H(n_lev)` is the top-of-atmosphere radius, and `xstar/ystar/zstar` are projected into the per-phase image basis. `trans_state` is 0 out of transit, 1 partial transit and 2 full atmospheric transit.

The data rows are:

`wavelength depth depth_atm depth_atm_east depth_atm_west depth_opaque`

`depth` is the total blocked fraction from the opaque body plus atmosphere. `depth_atm`, `depth_atm_east` and `depth_atm_west` are the atmosphere-only contribution and its east/west limb split. `depth_opaque` is the solid-body contribution alone, including the active limb-darkening law.

# What is output and how do I make synthetic observations?

The WASP-33b and WASP-39b examples use the following plotting/post-processing scripts; the interactive tool is shared from `tools/`:

`plot_trans.py` - transmission spectrum from `Transmission.txt`

`plot_trans_lc.py` - primary-transit light curve from `Transit_*.txt`

`tools/plot_trans_lc_interactive.py` - interactive primary-transit light-curve explorer from `Transit_*.txt`

`plot_em.py` - emission spectra from `Em_*.txt`

`plot_em_lc.py` - emission phase/eclipsing light curve from `Em_*.txt`

`plot_pixel.py` - pixel brightness-temperature maps from `f_im_*.txt`

`plot_3D_cf.py` - 3D contribution functions from `cf_*.txt`

# Example tutorials

## WASP-33b GCM emission, eclipse and transit-light-curve example

Example using a SPARC/MITgcm WASP-33b model. The default `CMCRT.nml` is prepared for emission phase-curve calculations with secondary-eclipse/occultation sampling enabled. The same namelist also contains dormant transmission-spectrum, albedo and primary-transit light-curve blocks for experimentation.
1. Use the extract script to extract the GCM data into the gCMCRT .hprf and .prf format, chemical abundances are extracted alongside.
2. Run goptools to produce the corr-k, CIA and Rayleigh opacity files.
3. Run gCMCRT to produce `Em_*.txt` output files.
4. Run `plot_em.py` to convert the output to synthetic observations, Fp/Fs, Fp and Tb. This file contains useful information on how to produce emission spectra.

For phase curves and eclipses, `plot_em_lc.py` reads the unocculted and occulted columns in the current `Em_*.txt` output and plots the band-mean Fp/Fs versus phase. For image outputs, `plot_pixel.py` reads the current `Em_*.txt` header and `f_im_*.txt` files.

To test the primary-transit light-curve mode in the WASP-33b directory, change `xper` to `3D_sph_trans_lc`. The light-curve mode writes `Transit_*.txt`; run `plot_trans_lc.py` from the WASP-33b example directory to plot transit spectra, normalised flux versus time from mid-transit, and the wavelength-time normalised-flux map.

For interactive exploration, run `python ../tools/plot_trans_lc_interactive.py --pattern "Transit_*.txt" --period-days 1.21987` from the WASP-33b example directory. The sliders select a wavelength bin and averaging width, and clicking a light-curve point shows the spectrum for that transit phase.

## WASP-39b GCM transmission spectrum and transit light curve (with details on using CE interpolation tables)

Example using an Exo-FMS WASP-39b model. The default `CMCRT.nml` produces a transmission spectrum, and the same namelist also contains dormant blocks for emission, secondary-eclipse experiments and a `3D_sph_trans_lc` primary-transit light-curve test setup.
1. Use the extract script to extract the GCM data into the gCMCRT .hprf and .iprf format
2. Use interp_iprf.py to interpolate the CE abundances to the T,p of the GCM and produce the .prf file.
3. Run goptools to produce the corr-k, CIA and Rayleigh opacity files.
4. Run gCMCRT to produce the Transmission.txt file 
5. Run `plot_trans.py` to convert the output to synthetic observations, Rp/Rs and compare to G395H and SOSS observations. (This file contains useful information on transmission spectrum fitting etc)

To test the transit light-curve mode, change `xper` in `WASP-39b_Exo_FMS_example/CMCRT.nml` from `3D_sph_trans` to `3D_sph_trans_lc`. The example light-curve block uses all 503 wavelength bins, 26 transit samples and a modest photon count; increase `Nph`, `n_trans` or `n_ingress` for production runs.

The light-curve mode writes `Transit_*.txt` files. Run `plot_trans_lc.py` from the WASP-39b example directory to plot:

1. phase-resolved transmission spectra
2. band-mean transmission diagnostics
3. normalised flux versus time from mid-transit
4. wavelength-time normalised-flux map

For interactive exploration, run `python ../tools/plot_trans_lc_interactive.py --pattern "Transit_*.txt" --period-days 4.055` from the WASP-39b example directory. The sliders select a wavelength bin and averaging width, and clicking a light-curve point shows the spectrum for that transit phase.

To experiment with WASP-39b emission, change `xper` to `3D_sph_em`. The provided `&sph_3D_em` block is a single-view setup by default. To run an emission phase curve, set `do_phase = .True.`. To include secondary-eclipse/occultation sampling, set both `do_phase = .True.` and `do_eclipse = .True.`; `n_phase` and `n_ecl` are already populated with example values.

## WASP-39b 1D VULCAN + CARMA model (with details on using the code in 1D)

## MALBEC benchmarks (further examples on 1D code use and comparing to other codes)
