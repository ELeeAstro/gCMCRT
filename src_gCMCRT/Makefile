F90          = nvfortran
F90LINKER    = nvfortran

# for nvfortran Compiler
#===================
DEFS      =
#FFLAGS   = -O0 -g -C -traceback -cuda -cudalib=curand
#FFLAGS   = -fast -g -traceback -cuda
FFLAGS    = -fast -O3 -cuda -cudalib=curand 
INCLUDES  =
LFLAGS    = $(FFLAGS)

OBJECTS = mc_precision.o \
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
exp_3D_sph_atm_tests.o \
exp_3D_sph_atm_pol.o \
exp_3D_sph_atm_alb.o \
exp_3D_sph_atm_trans.o \
exp_3D_sph_atm_trans_hires.o \
exp_3D_sph_atm_em.o \
exp_3D_sph_atm_em_hires.o \
gpuCMCRT.o

# exp_1D_pp_atm.o \
# exp_1D_sph_atm.o \
# exp_3D_cart_galaxy.o \
# exp_3D_sph_atm.o \

# executable statement
EXECS  = ../gCMCRT

.SUFFIXES : .o .f90

default: gCMCRT

gCMCRT:  $(OBJECTS)
	 $(F90LINKER) $(LFLAGS) $(OBJECTS) -o $(EXECS)

clean:
	rm -f *.o *.mod  $(EXECS)

.f90.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<
