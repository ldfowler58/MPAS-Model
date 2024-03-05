.SUFFIXES: .F90 .o

all: dummy chemistry_gocart2G

dummy:
	echo "****** compiling chemistry_gocart2G ******"

OBJS = \
	GOCART2G_MieMod.o \
	GOCART2G_Process.o

chemistry_gocart2G: $(OBJS)
	ar -ru ./../libchem.a $(OBJS)


# DEPENDENCIES:
GOCART2G_Process.o: \
	GOCART2G_MieMod.o

clean:
	$(RM) *.f90 *.o *.mod
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.F90.o:
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS) -c $*.F90 $(CPPINCLUDES) $(FCINCLUDES) -I.. -I../../../framework -I../../../external/esmf_time_f90
