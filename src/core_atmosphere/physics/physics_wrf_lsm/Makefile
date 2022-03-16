.SUFFIXES: .F .o

all: dummy physics_wrf_lsm

dummy:
	echo "****** compiling physics_wrf_lsm ******"

OBJS = \
	module_bep_bem_helper.o \
	module_sf_gecros.o \
	module_sf_bem.o \
	module_sf_bep_bem.o \
	module_sf_bep.o \
	module_sf_noahdrv.o \
	module_sf_noahlsm.o \
	module_sf_noahlsm_glacial_only.o \
	module_sf_noah_seaice_drv.o \
	module_sf_noah_seaice.o \
	module_sf_noahmpdrv.o \
	module_sf_noahmp_glacier.o \
	module_sf_noahmp_groundwater.o \
	module_sf_noahmplsm.o \
	module_sf_urban.o

physics_wrf_lsm: $(OBJS)
	ar -ru ./../libphys.a $(OBJS)


#DEPENDENCIES;
module_sf_bep_bem.o: \
	module_bep_bem_helper.o \
	module_sf_bem.o \
	module_sf_urban.o

module_sf_bep.o: \
	module_bep_bem_helper.o \
	module_sf_urban.o

module_sf_noahdrv.o: \
	module_sf_bep_bem.o \
	module_sf_bep.o \
	module_sf_noahlsm.o \
	module_sf_noahlsm_glacial_only.o \
	module_sf_urban.o

module_sf_noahmpdrv.o: \
	module_sf_bep_bem.o \
	module_sf_bep.o \
	module_sf_gecros.o \
	module_sf_noahmp_glacier.o \
	module_sf_noahmp_groundwater.o \
	module_sf_noahmplsm.o \
	module_sf_urban.o

module_sf_noahmp_groundwater.o: \
	module_sf_noahmplsm.o

module_sf_noahmplsm.o: \
	module_sf_gecros.o

module_sf_noahlsm_glacial_only.o: \
	module_sf_noahlsm.o

module_sf_noah_seaice_drv.o: \
	module_sf_noah_seaice.o

module_sf_noah_seaice.o: \
	module_sf_noahlsm.o


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
