# Compiler settings
#===================
F90 = nvfortran
F90LINKER = $(F90)

# Default flags for release build
FFLAGS = -mp -fast -O3

# Debug flags for debugging purposes
DEBUG_FLAGS = -O0 -g -traceback

# Libraries and includes (define as needed)
INCLUDES =
LFLAGS =

# Define additional libraries or definitions as needed
LIBS =

# Source files (.f and .f90 mix)
OBJECTS = \
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
    optools.o

# Executable
EXECS = ../goptools

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
	$(F90) $(FFLAGS) $(INCLUDES) -c $<

# Pattern rule for compiling Fortran 77 (.f) files
%.o: %.f
	$(F90) $(FFLAGS) $(INCLUDES) -c $<

# Clean target (removes object files and the executable)
clean:
	rm -f *.o *~ *.mod $(EXECS)
