# gCMCRT - 3D Monte Carlo Radiative Transfer for exoplanet atmospheres with GPUs

More extensive documentation is in the making, for now here are some bare bones instructions.

Two example models are provided:

1. 1D 1500K benchmark case from Baundino et al. (2017).
2. A WASP-33b SPARC/MITgcm GCM model output example.

## k-tables and CE data

k-tables and chemical equilibirum interpolation tables can be found here:
https://drive.google.com/drive/folders/1HVa9xWK_GqOqknIErcXVszzhhztACExw?usp=sharing

## To compile

You will need a Nvidia GPU card on your system.

You will need to install the latest drivers from Nvidia: https://www.nvidia.com/Download/index.aspx

You will need to install the CUDA toolkit: https://developer.nvidia.com/cuda-toolkit

You will need to install the CUDA hpc sdk: https://developer.nvidia.com/hpc-sdk


# How to operate optools

To compile cd to src_optoools_V2 and enter 'make'.
To de-compile enter 'make clean'

# How to operate gCMCRT
