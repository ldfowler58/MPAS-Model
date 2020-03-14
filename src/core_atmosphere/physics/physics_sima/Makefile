.SUFFIXES: .F .o

all: dummy physics_sima

dummy:
	echo "****** compiling physics_sima ******"

OBJS = \
	bl_ysu.o             \
	module_bl_ysu.o      \
	ra_rrtmg_mcica.o     \
	ra_rrtmg_vinterp.o   \
	ra_rrtmg_shared.o    \
	ra_rrtmg_sw_mpas.o   \
	ra_rrtmg_sw.o        \
        ra_sw_radconst.o

physics_sima: $(OBJS)
	ar -ru ./../libphys.a $(OBJS)

# DEPENDENCIES:
module_bl_ysu.o: \
	bl_ysu.o

ra_rrtmg_mcica.o: \
	ra_rrtmg_sw.o

ra_rrtmg_sw_mpas.o: \
	ra_rrtmg_mcica.o \
	ra_rrtmg_shared.o \
	ra_rrtmg_vinterp.o \
	ra_sw_radconst.o

ra_rrtmg_sw.o: \
	ra_rrtmg_shared.o

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
