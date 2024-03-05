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
	CA2G_bc_GridCompMod.o \
	CA2G_br_GridCompMod.o \
	CA2G_oc_GridCompMod.o \
	DU2G_GridCompMod.o \
	NI2G_GridCompMod.o \
	SS2G_GridCompMod.o \
	SU2G_GridCompMod.o \
	CA2G_StateSpecs.o \
	DU2G_StateSpecs.o \
	NI2G_StateSpecs.o \
	SS2G_StateSpecs.o \
	SU2G_StateSpecs.o \
	GA_EnvironmentMod.o \
	GOCART2G_MieMod.o \
	GOCART2G_Process.o

chemistry_gocart2G: $(OBJS)
	ar -ru ./../libchem.a $(OBJS)


# DEPENDENCIES:
CA2G_bc_GridCompMod.o: \
	GA_EnvironmentMod.o

CA2G_br_GridCompMod.o: \
	GA_EnvironmentMod.o

CA2G_oc_GridCompMod.o: \
	GA_EnvironmentMod.o

DU2G_GridCompMod.o: \
	GA_EnvironmentMod.o

NI2G_GridCompMod.o: \
	GA_EnvironmentMod.o

SS2G_GridCompMod.o: \
	GA_EnvironmentMod.o

SU2G_GridCompMod.o: \
	GA_EnvironmentMod.o

GA_EnvironmentMod.o: \
	GOCART2G_MieMod.o

GOCART2G_Process.o: \
	GOCART2G_MieMod.o

clean:
	$(RM) *.f90 *.o *.mod
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.F90.o:
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS) -c $*.F90 $(CPPINCLUDES) $(FCINCLUDES) -I.. -I../../../framework -I../../../external/esmf_time_f90
