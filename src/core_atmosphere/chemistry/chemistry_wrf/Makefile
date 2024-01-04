.SUFFIXES: .F .o

all: dummy chemistry_wrf

dummy:
	echo "****** compiling chemistry_wrf ******"

OBJS = \
	module_data_gocart_chem.o \
	module_data_gocart_dust.o \
	module_data_gocart_seas.o \
	module_data_mosaic_asect.o \
	module_data_radm2.o \
	module_data_rrtmg_aeropt.o \
	module_data_sorgam.o \
	module_gocart_aerosols.o \
	module_gocart_chem.o \
	module_gocart_dms_emission.o \
	module_gocart_dust_afwa.o \
	module_gocart_dust.o \
	module_gocart_mie.o \
	module_gocart_opt.o \
	module_gocart_seasalt.o

chemistry_wrf: $(OBJS)
	ar -ru ./../libchem.a $(OBJS)


# DEPENDENCIES:
module_data_sorgam.o: \
	module_data_radm2.o

module_gocart_dms_emission.o: \
	module_data_gocart_chem.o

module_gocart_dust.o: \
	module_data_gocart_dust.o

module_gocart_mie.o: \
	module_data_rrtmg_aeropt.o

module_gocart_opt.o: \
	module_data_gocart_chem.o \
	module_data_gocart_dust.o \
	module_data_gocart_seas.o \
	module_data_mosaic_asect.o \
	module_data_rrtmg_aeropt.o \
	module_data_sorgam.o \
	module_gocart_mie.o

module_gocart_seas.o: \
	module_data_gocart_seas.o


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
