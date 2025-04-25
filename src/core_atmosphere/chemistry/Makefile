.SUFFIXES: .F90 .o


ifeq ($(CORE),atmosphere)
COREDEF = -Dmpas
endif


all: dummy core_gocart2G core_chemistry


dummy:
	echo "****** compiling chemistry ******"


OBJS = \
	mpas_chemistry_gocart2G_date_time.o \
	mpas_chemistry_gocart2G_diagnostics.o \
	mpas_chemistry_gocart2G_driver.o \
	mpas_chemistry_gocart2G_emissions.o \
	mpas_chemistry_gocart2G_emissions_init.o \
	mpas_chemistry_gocart2G_emissions_update.o \
	mpas_chemistry_gocart2G_finalize.o \
	mpas_chemistry_gocart2G_forMPASphysics.o \
	mpas_chemistry_gocart2G_forMPASrestart.o \
	mpas_chemistry_gocart2G_init.o \
	mpas_chemistry_gocart2G_interface.o \
	mpas_chemistry_gocart2G_manager.o \
	mpas_chemistry_gocart2G_update.o \
	mpas_chemistry_gocart2G_vars.o


core_gocart2G:
	( cd gocart2G; $(MAKE) all COREDEF="$(COREDEF)" )

core_chemistry: core_gocart2G
	( cd ../physics; $(MAKE) core_physics_init )
	( $(MAKE) chem_interface COREDEF="$(COREDEF)" )
	ar -ru libchem.a $(OBJS)

chem_interface: $(OBJS)


# DEPENDENCIES:
mpas_chemistry_gocart2G_driver.o: \
	mpas_chemistry_gocart2G_diagnostics.o \
	mpas_chemistry_gocart2G_interface.o \
	mpas_chemistry_gocart2G_manager.o \
	mpas_chemistry_gocart2G_vars.o

mpas_chemistry_gocart2G_emissions_init.o: \
	mpas_chemistry_gocart2G_vars.o

mpas_chemistry_gocart2G_finalize.o: \
	mpas_chemistry_gocart2G_vars.o

mpas_chemistry_gocart2G_forMPASrestart.o: \
	mpas_chemistry_gocart2G_vars.o

mpas_chemistry_gocart2G_init.o: \
	mpas_chemistry_gocart2G_vars.o

mpas_chemistry_gocart2G_interface.o: \
	mpas_chemistry_gocart2G_update.o

mpas_chemistry_gocart2G_manager.o: \
	mpas_chemistry_gocart2G_emissions_init.o \
	mpas_chemistry_gocart2G_emissions_update.o \
	mpas_chemistry_gocart2G_update.o

mpas_chemistry_gocart2G_update.o: \
	mpas_chemistry_gocart2G_date_time.o

mpas_chemistry_gocart2G_vars.o: \
	mpas_chemistry_gocart2G_emissions.o \
	mpas_chemistry_gocart2G_interface.o \
	mpas_chemistry_gocart2G_forMPASphysics.o


clean:
	$(RM) *.f90 *.o *.mod libchem.a
	( cd gocart2G; $(MAKE) clean )
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i


.F90.o:
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I. -I./gocart2G -I.. -I../dynamics -I../physics -I../../framework -I../../external/esmf_time_f90 -I../../external/SMIOL
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(HYDROSTATIC) $(PHYSICS) $(FFLAGS) -c $*.F90 $(CPPINCLUDES) $(FCINCLUDES) -I. -I./gocart2G -I.. -I../dynamics -I../physics -I../../framework -I../../external/esmf_time_f90  -I../../external/SMIOL
endif
