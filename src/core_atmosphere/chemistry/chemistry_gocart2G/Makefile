.SUFFIXES: .F .o

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

.F.o:
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I.. -I../../../framework -I../../../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I.. -I../../../framework -I../../../external/esmf_time_f90
endif
