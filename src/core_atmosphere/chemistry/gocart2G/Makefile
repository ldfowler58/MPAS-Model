.SUFFIXES: .F90 .o


all: dummy gocart2G


dummy:
	echo "****** compiling gocart2G ******"


OBJS = \
	CA2G_bc_instance.o \
	CA2G_br_instance.o \
	CA2G_oc_instance.o \
	DU2G_instance.o \
	NI2G_instance.o \
	SS2G_instance.o \
	SU2G_instance.o \
	CA2G_bc_GridCompMod.o \
	CA2G_br_GridCompMod.o \
	CA2G_oc_GridCompMod.o \
	DU2G_GridCompMod.o \
	NI2G_GridCompMod.o \
	SOA2G_GridCompMod.o \
	SS2G_GridCompMod.o \
	SU2G_GridCompMod.o \
	CA2G_bc_StateSpecs.o \
	CA2G_br_StateSpecs.o \
	CA2G_oc_StateSpecs.o \
	DU2G_StateSpecs.o \
	NI2G_StateSpecs.o \
	SOA2G_StateSpecs.o \
	SS2G_StateSpecs.o \
	SU2G_StateSpecs.o \
	GA_EnvironmentMod.o \
	GOCART2G_instance.o \
	GOCART2G_GridCompMod.o \
	GOCART2G_MieMod_smiol.o \
	GOCART2G_Process.o \
	GOCART2G_StateSpecs.o \
	MAPL.o


gocart2G: $(OBJS)
	ar -ru ./../libchem.a $(OBJS)


# DEPENDENCIES:
CA2G_bc_GridCompMod.o: \
	CA2G_bc_instance.o \
	CA2G_bc_StateSpecs.o \
	GA_EnvironmentMod.o \
	GOCART2G_AeroGeneric.o \
	GOCART2G_MieMod_smiol.o \
	GOCART2G_Process.o \
	MAPL.o

CA2G_br_GridCompMod.o: \
	CA2G_br_instance.o \
	CA2G_br_StateSpecs.o \
	GA_EnvironmentMod.o \
	GOCART2G_AeroGeneric.o \
	GOCART2G_MieMod_smiol.o \
	GOCART2G_Process.o

CA2G_oc_GridCompMod.o: \
	CA2G_oc_instance.o \
	CA2G_oc_StateSpecs.o \
	GA_EnvironmentMod.o \
	GOCART2G_AeroGeneric.o \
	GOCART2G_MieMod_smiol.o \
	GOCART2G_Process.o

DU2G_GridCompMod.o: \
	DU2G_instance.o \
	DU2G_StateSpecs.o \
	GA_EnvironmentMod.o \
	GOCART2G_AeroGeneric.o \
	GOCART2G_MieMod_smiol.o \
	GOCART2G_Process.o

GOCART2G_GridCompMod.o: \
	GA_EnvironmentMod.o \
	GOCART2G_instance.o \
	GOCART2G_StateSpecs.o \
	CA2G_bc_StateSpecs.o \
	CA2G_br_StateSpecs.o \
	CA2G_oc_StateSpecs.o \
	DU2G_StateSpecs.o \
	NI2G_StateSpecs.o \
	SS2G_StateSpecs.o \
	SU2G_StateSpecs.o

NI2G_GridCompMod.o: \
	NI2G_instance.o \
	NI2G_StateSpecs.o \
	GA_EnvironmentMod.o \
	DU2G_GridCompMod.o \
	SS2G_GridCompMod.o \
	GOCART2G_AeroGeneric.o \
	GOCART2G_MieMod_smiol.o \
	GOCART2G_Process.o

SOA2G_GridCompMod.o: \
	SOA2G_StateSpecs.o

SS2G_GridCompMod.o: \
	SS2G_instance.o \
	SS2G_StateSpecs.o \
	GA_EnvironmentMod.o \
	GOCART2G_AeroGeneric.o \
	GOCART2G_MieMod_smiol.o \
	GOCART2G_Process.o

SU2G_GridCompMod.o: \
	SU2G_instance.o \
	SU2G_StateSpecs.o \
	GA_EnvironmentMod.o \
	GOCART2G_AeroGeneric.o \
	GOCART2G_MieMod_smiol.o \
	GOCART2G_Process.o \
	MAPL.o

CA2G_bc_StateSpecs.o: \
	CA2G_bc_instance.o \
	GOCART2G_instance.o

CA2G_br_StateSpecs.o: \
	CA2G_br_instance.o \
	GOCART2G_instance.o

CA2G_oc_StateSpecs.o: \
	CA2G_oc_instance.o \
	GOCART2G_instance.o

DU2G_StateSpecs.o: \
	DU2G_instance.o \
	GOCART2G_instance.o

GOCART2G_Process.o: \
	GOCART2G_MieMod_smiol.o

GOCART2G_StateSpecs.o : \
	GOCART2G_instance.o

NI2G_StateSpecs.o: \
	NI2G_instance.o \
	GOCART2G_instance.o

SS2G_StateSpecs.o: \
	SS2G_instance.o \
	GOCART2G_instance.o

SU2G_StateSpecs.o: \
	SU2G_instance.o \
	GOCART2G_instance.o

GA_EnvironmentMod.o: \
	GOCART2G_instance.o \
	GOCART2G_MieMod_smiol.o


clean:
	$(RM) *.f90 *.o *.mod
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i


.F90.o:
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) $< > $*.f90 -I../../../external/SMIOL
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I.. -I../../../framework -I../../../external/esmf_time_f90 -I../../../external/SMIOL
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS) -c $*.F90 $(CPPINCLUDES) $(FCINCLUDES) -I.. -I../../../framework -I../../../external/esmf_time_f90 -I../../../external/SMIOL
endif
