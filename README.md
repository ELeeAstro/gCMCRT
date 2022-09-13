# gCMCRT - 3D Monte Carlo Radiative Transfer for exoplanet atmospheres with GPUs

WARNING - this code is not a black box and requires some getting used to, that said it typically takes a student a few hours of tinkering to get a first spectra. Once the model works the first time, further tinkering with options/physics etc becomes more simple.
The best way to learn the code is to use it! And feel free to contact Elsie or other experianced users if stuck or confused.

More extensive documentation is in the making, for now here are some bare bones instructions.

Two example models are provided:

1. A WASP-33b SPARC/MITgcm GCM model emission example - which includes chemical abundances
2. A WASP-39b Exo-FMS GCM model transmission example - uses chemical equilibirum abundance table interpolation
3. A WASP-39b 1D VULCAN + CARMA output example - uses VULCAN

## k-tables and CE data

k-tables and chemical equilibirum interpolation tables can be found here:
https://drive.google.com/drive/folders/1HVa9xWK_GqOqknIErcXVszzhhztACExw?usp=sharing


## To compile

You will need a Nvidia GPU card on your system.

You will need to install the latest drivers from Nvidia: https://www.nvidia.com/Download/index.aspx

(optional) install the CUDA toolkit: https://developer.nvidia.com/cuda-toolkit

You will need to install the CUDA hpc sdk: https://developer.nvidia.com/hpc-sdk

# Required gCMCRT formated files

.prf file

.hprf file

.clprf (for clouds)

.iprf (if CE interpolation required)

# How to operate optools

To compile cd to src_optoools_V2 and enter 'make'.
To de-compile enter 'make clean'.

optools uses a fortran namelist and parameter file to communicate with the code.

## optools.par file

Is quite self-explainatory - fill in the number of species followed by the number of species (See examples)

## optools.nml file

Is more difficult to fill out correctly:

### &CK_nml - corr-k namelist

pre_mixed - Does a pre-mixed table interpolation (.True.), otherwise random overlap (.False.)

interp_wl - Interpolate to the wavelengths.wl file (.True.) (.False. if exact wavelengths of the ck table are used)

iopts - Integer option number (dev-only)

form - 1 (NEMESIS format), 2 (gCMCRT format) Note, for multiple k-tables this is a comma separated list, e.g. 2,2,2

nG - number of g-ordinances in k-table 

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

### &cl_nml - clouds namelist

iopts - Integer option number (dev-only)

imix = 1

idist - 0 (Read bin model results), 1 (single particle size), 2 (3 size peaked near mean size), 3 (log-normal), 4 (Gamma), 5 (Inv. Gamma), 6 (Rayleigh), 7 (Hansen)

ndist - number of size distibution points (log-spaced between amin and amax)

idist_int - 1 (Trapezium rule integration)

imie - 1 (MieX solver), 2 (DHS solver)

form - 5 (DIHRT format nk-tables)

paths - list of paths to the nk data
\NOTE: THESE PATHS MUST BE IN THE SAME SPECIES ORDER AS THE SPECIES IN THE optools FILE !!!

sig - ln sigma value (log-normal)

eff_fac - effetive mean size varience

veff - effective varience (Hansen)

amin - Minimum distribution particle size (um)

amax - Maximum distribution particle size (um)

fmax - parameter for DHS theory


# How to operate gCMCRT
