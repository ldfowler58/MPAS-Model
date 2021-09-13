.SUFFIXES: .F90 .o

all: dummy physics_cam6_utils physics_cam6

dummy:
	echo "****** compiling physics_cam6 ******"

OBJS = \
	mp_pumas_mpas.o \
	mp_pumas_vars.o \
	micro_pumas_utils.o \
	micro_mg3_0.o \
	wv_sat_methods.o

physics_cam6_utils:
	(cd utils; $(MAKE))

physics_cam6: $(OBJS)
	ar -ru ./../libphys.a $(OBJS)


# DEPENDENCIES:
micro_mg3_0.o: \
	micro_pumas_utils.o \
	wv_sat_methods.o

mp_pumas_mpas.o: \
	mp_pumas_vars.o

clean:
	$(RM) *.f90 *.o *.mod
	( cd utils; $(MAKE) clean )
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.F90.o:
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I.. -I../../../framework -I./utils
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS) -c $*.F90 $(CPPINCLUDES) $(FCINCLUDES) -I.. -I../../../framework -I./utils
endif
