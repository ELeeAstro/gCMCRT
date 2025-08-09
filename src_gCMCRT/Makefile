# Compiler settings
#===================
F90 = nvfortran
F90LINKER = $(F90)

# Default flags for release build
FFLAGS = -fast -O3 -cuda -cudalib=curand

# Debug flags for debugging purposes
DEBUG_FLAGS = -O0 -g -C -traceback -cuda -cudalib=curand

# Include directories and libraries
INCLUDES =
LFLAGS = $(FFLAGS)
LIBS =

# Source files
OBJECTS = \
    mc_precision.o \
    mc_data_mod.o \
    mc_class_pac.o \
    mc_class_grid.o \
    mc_class_imag.o \
    mc_k_aux.o \
    mc_k_moments.o \
    mc_k_vol_samp.o \
    mc_k_limb_dark.o \
    mc_k_scatt_mat.o \
    mc_k_scatt.o \
    mc_k_findcell.o \
    mc_k_findwall_cart.o \
    mc_k_findwall_sph.o \
    mc_k_raytrace.o \
    mc_k_peeloff_emit.o \
    mc_k_peeloff_scatt.o \
    mc_k_source_pac_inc.o \
    mc_k_emit_iso.o \
    mc_k_tauint.o \
    mc_k_tau_samp.o \
    mc_k_gord_samp.o \
    mc_k_RR.o \
    mc_read_namelist.o \
    mc_read_prf.o \
    mc_Draine_G.o \
    mc_set_opac.o \
    mc_set_em.o \
    mc_los_velocity.o \
    mc_LHS_sampling.o \
    exp_3D_sph_atm_tests.o \
    exp_3D_sph_atm_pol.o \
    exp_3D_sph_atm_alb.o \
    exp_3D_sph_atm_trans.o \
    exp_3D_sph_atm_trans_hires.o \
    exp_3D_sph_atm_em.o \
    exp_3D_sph_atm_em_hires.o \
    gpuCMCRT.o

# Executable
EXECS = ../gCMCRT

#====================
# Phony targets
.PHONY: all clean debug release

# Default target (release build by default)
all: $(EXECS)

# Debug build (overrides FFLAGS with debug flags)
debug: FFLAGS = $(DEBUG_FLAGS)
debug: clean all

# Release build (uses default release flags)
release: clean all

# Linking the executable
$(EXECS): $(OBJECTS)
	$(F90LINKER) $(FFLAGS) $(OBJECTS) -o $(EXECS) $(LIBS)

# Pattern rule for compiling Fortran 90 files
%.o: %.f90
	$(F90) $(FFLAGS) $(DEFS) -c $<

# Clean target (removes object files and the executable)
clean:
	rm -f *.o *.mod $(EXECS)
