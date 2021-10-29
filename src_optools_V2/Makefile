F90          = nvfortran
F90LINKER    = nvfortran

# for nvfortran Compiler
#===================
DEFS      =
#FFLAGS    = -O0 -g -traceback
FFLAGS    = -mp -fast -Mipa=fast
INCLUDES  =
LFLAGS    =

OBJECTS = \
MIEX.o \
q_dhs.o \
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
cloud_tables_read.o \
cloud_tables_interp.o \
cloud_tables_emt.o \
cloud_tables_mie.o \
cloud_tables_dist.o \
cloud_tables_mod.o \
CIA_tables_Hminus.o \
CIA_tables_Heminus.o \
CIA_tables_interp.o \
CIA_tables_read.o \
CIA_tables_mod.o \
Ray_tables_mod.o \
optools.o


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
