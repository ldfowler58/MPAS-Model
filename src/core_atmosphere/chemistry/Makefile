.SUFFIXES: .F .o

ifeq ($(CORE),atmosphere)
COREDEF = -Dmpas
endif

all: dummy core_chemistry_wrf core_chemistry_gocart2G core_chemistry

dummy:
	echo "****** compiling chemistry ******"

OBJS = \
	module_mynn_get_pblh.o \
	mpas_atmphys_chems.o \
	mpas_atmphys_chems_physics.o \
	mpas_atmphys_tracers.o \
	mpas_atmphys_tracers_physics.o \
	mpas_chemistry_gocart_driver_aerosols.o \
	mpas_chemistry_gocart_driver_chem.o \
	mpas_chemistry_gocart_driver_co.o \
	mpas_chemistry_gocart_driver_dms.o \
	mpas_chemistry_gocart_driver_dust.o \
	mpas_chemistry_gocart_driver_opt.o \
	mpas_chemistry_gocart_driver_seasalt.o \
	mpas_chemistry_gocart_driver.o \
	mpas_chemistry_manager.o \
	mpas_chemistry_gocart_interface.o \
	mpas_chemistry_gocart2G_driver.o \
	mpas_chemistry_gocart2G_init.o \
	mpas_chemistry_gocart2G_interface.o \
	mpas_chemistry_gocart2G_vars.o


core_chemistry_wrf:
	(cd chemistry_wrf; $(MAKE) all COREDEF="$(COREDEF)")

core_chemistry_gocart2G:
	(cd chemistry_gocart2G; $(MAKE) all COREDEF="$(COREDEF)")

core_chemistry: core_chemistry_wrf
	($(MAKE) chem_interface COREDEF="$(COREDEF)")
	ar -ru libchem.a $(OBJS)

chem_interface: $(OBJS)


# DEPENDENCIES:
mpas_atmphys_tracers.o: \
	module_mynn_get_pblh.o

mpas_chemistry_gocart_driver_aerosols.o: \
	mpas_chemistry_gocart_interface.o

mpas_chemistry_gocart_driver_chem.o: \
	mpas_chemistry_gocart_interface.o

mpas_chemistry_gocart_driver_co.o: \
	mpas_chemistry_gocart_interface.o

mpas_chemistry_gocart_driver_dms.o: \
	mpas_chemistry_gocart_interface.o

mpas_chemistry_gocart_driver_dust.o: \
	mpas_chemistry_gocart_interface.o

mpas_chemistry_gocart_driver_seasalt.o: \
	mpas_chemistry_gocart_interface.o

mpas_chemistry_gocart_driver.o: \
	mpas_chemistry_gocart_driver_aerosols.o \
	mpas_chemistry_gocart_driver_chem.o \
	mpas_chemistry_gocart_driver_co.o \
	mpas_chemistry_gocart_driver_dms.o \
	mpas_chemistry_gocart_driver_dust.o \
	mpas_chemistry_gocart_driver_opt.o \
	mpas_chemistry_gocart_driver_seasalt.o \
	mpas_chemistry_gocart_interface.o

mpas_chemistry_gocart2G_driver.o:\
	mpas_chemistry_gocart2G_interface.o \
	mpas_chemistry_gocart2G_vars.o

mpas_chemistry_gocart2G_init.o:\
	mpas_chemistry_gocart2G_interface.o \
	mpas_chemistry_gocart2G_vars.o


clean:
	$(RM) *.o *.mod *.f90 libchem.a
	( cd chemistry_wrf; $(MAKE) clean )
	( cd chemistry_gocart2G; $(MAKE) clean )
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(COREDEF) $(HYDROSTATIC) $(PHYSICS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I. -I./chemistry_wrf -I./chemistry_gocart2G -I.. -I../dynamics -I../physics -I../../framework -I../../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(HYDROSTATIC) $(PHYSICS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I. -I./chemistry_wrf -I./chemistry_gocart2G -I.. -I../dynamics -I../physics -I../../framework -I../../external/esmf_time_f90
endif
