F90          = nvfortran
F90LINKER    = nvfortran

# for nvfortran Compiler
#===================
DEFS      =
#FFLAGS    = -O0 -g -traceback
FFLAGS    = -mp -fast -O3 -Mipa=fast
#FFLAGS    = -fast -O3 -Mipa=fast
#FFLAGS = -O2 -C -g -traceback
INCLUDES  =
LFLAGS    =

#F90          = gfortran
#F90LINKER    = gfortran

# for gfortran Compiler
#===================
DEFS      =
#FFLAGS    = -Og -g -pipe -Wall -Wextra -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow
#FFLAGS    = 
#FFLAGS   = 
INCLUDES  =
LFLAGS    =

OBJECTS = \
MIEX.o \
mieext_mod.o \
bhmie_mod.o \
bhcoat_mod.o \
DMiLay.o \
dhs_mod.o \
lxmie_mod.o \
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

#dhs_mod.o: FFLAGS = -mp -O1

# executable statement
EXECS  = ../goptools

.SUFFIXES : .o .f .f90

default: goptools

goptools: $(OBJECTS)
	$(F90LINKER) $(INCLUDES) $(LFLAGS) $(OBJECTS) -o $(EXECS)

clean:
	rm -f *.o *~ *.mod $(EXECS)

.f.o:
	$(F90) $(FFLAGS) $(INCLUDES)  $(DEFS) -c $< $(LFLAGS)
.f90.o:
	$(F90) $(FFLAGS) $(INCLUDES)  $(DEFS) -c $< $(LFLAGS)
