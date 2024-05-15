.SUFFIXES: .F90 .o

all: dummy chemistry_gocart2G

dummy:
	echo "****** compiling chemistry_gocart2G ******"

OBJS = \
	CA2G_bc_instance_CA.o \
	CA2G_br_instance_CA.o \
	CA2G_oc_instance_CA.o \
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
	CA2G_bc_StateSpecs.o \
	CA2G_br_StateSpecs.o \
	CA2G_oc_StateSpecs.o \
	DU2G_StateSpecs.o \
	NI2G_StateSpecs.o \
	SS2G_StateSpecs.o \
	SU2G_StateSpecs.o \
	GA_EnvironmentMod.o \
	GOCART2G_AeroGeneric.o \
	GOCART2G_GridComp.o \
	GOCART2G_MieMod_smiol.o \
	GOCART2G_Process.o \
	GOCART2G_StateSpecs.o


chemistry_gocart2G: $(OBJS)
	ar -ru ./../libchem.a $(OBJS)


# DEPENDENCIES:
CA2G_bc_GridCompMod.o: \
	CA2G_bc_instance_CA.o \
	CA2G_bc_StateSpecs.o \
	GA_EnvironmentMod.o \
	GOCART2G_AeroGeneric.o \
	GOCART2G_MieMod_smiol.o \
	GOCART2G_Process.o

CA2G_br_GridCompMod.o: \
	CA2G_br_instance_CA.o \
	CA2G_br_StateSpecs.o \
	GA_EnvironmentMod.o \
	GOCART2G_AeroGeneric.o \
	GOCART2G_MieMod_smiol.o \
	GOCART2G_Process.o

CA2G_oc_GridCompMod.o: \
	CA2G_oc_instance_CA.o \
	CA2G_oc_StateSpecs.o \
	GA_EnvironmentMod.o \
	GOCART2G_AeroGeneric.o \
	GOCART2G_MieMod_smiol.o \
	GOCART2G_Process.o

DU2G_GridCompMod.o: \
	DU2G_instance_DU.o \
	DU2G_StateSpecs.o \
	GA_EnvironmentMod.o \
	GOCART2G_AeroGeneric.o \
	GOCART2G_MieMod_smiol.o \
	GOCART2G_Process.o

NI2G_GridCompMod.o: \
	NI2G_instance_NI.o \
	NI2G_StateSpecs.o \
	GA_EnvironmentMod.o \
	GOCART2G_AeroGeneric.o \
	GOCART2G_MieMod_smiol.o \
	GOCART2G_Process.o

SS2G_GridCompMod.o: \
	SS2G_instance_SS.o \
	SS2G_StateSpecs.o \
	GA_EnvironmentMod.o \
	GOCART2G_AeroGeneric.o \
	GOCART2G_MieMod_smiol.o \
	GOCART2G_Process.o

SU2G_GridCompMod.o: \
	SU2G_instance_SU.o \
	SU2G_StateSpecs.o \
	GA_EnvironmentMod.o \
	GOCART2G_AeroGeneric.o \
	GOCART2G_MieMod_smiol.o \
	GOCART2G_Process.o

CA2G_bc_StateSpecs.o: \
	CA2G_bc_instance_CA.o

CA2G_br_StateSpecs.o: \
	CA2G_br_instance_CA.o

CA2G_oc_StateSpecs.o: \
	CA2G_oc_instance_CA.o

DU2G_StateSpecs.o: \
	DU2G_instance_DU.o

SU2G_StateSpecs.o: \
	GOCART2G_GridComp.o \
	SU2G_instance_SU.o

GA_EnvironmentMod.o: \
	GOCART2G_MieMod_smiol.o

GOCART2G_MieMod_smiol.o: \
	GOCART2G_MieMod_smiol.o


clean:
	$(RM) *.f90 *.o *.mod
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.F90.o:
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS) -c $*.F90 $(CPPINCLUDES) $(FCINCLUDES) -I.. -I../../../framework -I../../../external/esmf_time_f90
