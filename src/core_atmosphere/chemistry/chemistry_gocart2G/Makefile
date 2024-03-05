.SUFFIXES: .F90 .o

all: dummy chemistry_gocart2G

dummy:
	echo "****** compiling chemistry_gocart2G ******"

OBJS = \
	CA2G_instance_CA_bc.o \
	CA2G_instance_CA_br.o \
	CA2G_instance_CA_oc.o \
	DU2G_instance_DU.o \
	NI2G_instance_NI.o \
	SS2G_instance_SS.o \
	SU2G_instance_SU.o \
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
