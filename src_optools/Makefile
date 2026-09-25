# optools CPU/GPU build configuration
# ===================================
#
# Both executables are built from the same source list. GPU_OFFLOAD is only
# defined for the GPU build so source files can select OpenMP target regions
# as they are introduced.

NVFORTRAN ?= nvfortran
CPU_FC ?= $(NVFORTRAN)
GPU_FC ?= $(NVFORTRAN)

# Raven's NVIDIA A100 GPUs use CUDA compute capability 8.0.
GPU_ARCH ?= cc80
GPU_LBL_BLOCK ?= 64
GPU_CK_BLOCK ?= 8

COMMON_FLAGS ?= -fast -O3 -Mpreprocess
CPU_FLAGS ?= -mp=multicore $(COMMON_FLAGS)
GPU_FLAGS ?= -mp=gpu -gpu=$(GPU_ARCH) -Minfo=mp $(COMMON_FLAGS)

CPU_CPPFLAGS ?=
GPU_CPPFLAGS ?= -DGPU_OFFLOAD -DLBL_GPU_BLOCK=$(GPU_LBL_BLOCK) \
	-DCK_GPU_BLOCK=$(GPU_CK_BLOCK)

CPU_MODULE_FLAGS ?= -module $(CPU_BUILD) -I$(CPU_BUILD)
GPU_MODULE_FLAGS ?= -module $(GPU_BUILD) -I$(GPU_BUILD)

CPU_DEBUG_FLAGS ?= -mp=multicore -O0 -g -traceback -Mbounds -Mpreprocess
GPU_DEBUG_FLAGS ?= -mp=gpu -gpu=$(GPU_ARCH),debug,lineinfo -O0 -g -traceback -Mbounds -Mpreprocess

INCLUDES ?=
LDFLAGS ?=
LIBS ?=

BUILD_ROOT ?= .build
CPU_BUILD := $(BUILD_ROOT)/cpu
GPU_BUILD := $(BUILD_ROOT)/gpu

CPU_EXE ?= ../goptools_cpu
GPU_EXE ?= ../goptools_gpu

# Source order also supplies the Fortran module compilation order. Parallel
# make is disabled until explicit generated module dependencies are added.
OBJECT_NAMES = \
    MIEX.o \
    mieext_mod.o \
    bhmie_mod.o \
    bhcoat_mod.o \
    DMiLay.o \
    dhs_mod.o \
    lxmie_mod.o \
    mie_approx_mod.o \
    optools_table_class.o \
    optools_data_mod.o \
    optools_aux.o \
    optools_read_io.o \
    CK_tables_RO.o \
    CK_tables_interp.o \
    CK_tables_read.o \
    CK_tables_mod.o \
    lbl_tables_read.o \
    lbl_tables_interp.o \
    lbl_tables_combine.o \
    lbl_tables_mod.o \
    xsec_tables_read.o \
    xsec_tables_interp.o \
    xsec_tables_mod.o \
    cloud_tables_read.o \
    cloud_tables_interp.o \
    cloud_tables_emt.o \
    cloud_tables_mie.o \
    cloud_tables_dist.o \
    cloud_tables_carma.o \
    cloud_tables_mod.o \
    CIA_tables_Hminus.o \
    CIA_tables_Heminus.o \
    CIA_tables_H2minus.o \
    CIA_tables_fake_H2O_special.o \
    CIA_tables_interp.o \
    CIA_tables_read.o \
    CIA_tables_mod.o \
    Ray_tables_mod.o \
    gas_refractivity_mod.o \
    refrac_tables_mod.o \
    optools.o

CPU_OBJECTS := $(addprefix $(CPU_BUILD)/,$(OBJECT_NAMES))
GPU_OBJECTS := $(addprefix $(GPU_BUILD)/,$(OBJECT_NAMES))

.NOTPARALLEL:

.PHONY: all cpu gpu release debug cpu-debug gpu-debug clean clean-cpu clean-gpu

# Keep the default build usable on systems without a GPU.
all: cpu

cpu: $(CPU_EXE)

gpu: $(GPU_EXE)

release: cpu

debug: cpu-debug

cpu-debug: clean-cpu
	$(MAKE) CPU_FLAGS="$(CPU_DEBUG_FLAGS)" cpu

gpu-debug: clean-gpu
	$(MAKE) GPU_FLAGS="$(GPU_DEBUG_FLAGS)" gpu

$(CPU_EXE): $(CPU_OBJECTS)
	$(CPU_FC) $(CPU_FLAGS) $(LDFLAGS) $^ -o $@ $(LIBS)

$(GPU_EXE): $(GPU_OBJECTS)
	$(GPU_FC) $(GPU_FLAGS) $(LDFLAGS) $^ -o $@ $(LIBS)

$(CPU_BUILD)/%.o: %.f90 | $(CPU_BUILD)
	$(CPU_FC) $(CPU_CPPFLAGS) $(CPU_FLAGS) $(INCLUDES) \
		$(CPU_MODULE_FLAGS) -c $< -o $@

$(CPU_BUILD)/%.o: %.f | $(CPU_BUILD)
	$(CPU_FC) $(CPU_CPPFLAGS) $(CPU_FLAGS) $(INCLUDES) \
		$(CPU_MODULE_FLAGS) -c $< -o $@

$(GPU_BUILD)/%.o: %.f90 | $(GPU_BUILD)
	$(GPU_FC) $(GPU_CPPFLAGS) $(GPU_FLAGS) $(INCLUDES) \
		$(GPU_MODULE_FLAGS) -c $< -o $@

$(GPU_BUILD)/%.o: %.f | $(GPU_BUILD)
	$(GPU_FC) $(GPU_CPPFLAGS) $(GPU_FLAGS) $(INCLUDES) \
		$(GPU_MODULE_FLAGS) -c $< -o $@

$(CPU_BUILD) $(GPU_BUILD):
	mkdir -p $@

clean-cpu:
	$(RM) -r $(CPU_BUILD) $(CPU_EXE)

clean-gpu:
	$(RM) -r $(GPU_BUILD) $(GPU_EXE)

clean: clean-cpu clean-gpu
	-rmdir $(BUILD_ROOT) 2>/dev/null
