.SUFFIXES: .F90 .o

all: dummy physics_cam6_utils

dummy:
	echo "****** compiling physics_cam6_utils ******"

OBJS = \
	shr_const_mod.o \
	shr_kind_mod.o \
	shr_spfn_mod.o

physics_cam6_utils: $(OBJS)
	ar -ru ./../../libphys.a $(OBJS)


# DEPENDENCIES:
shr_const_mod.o: \
	shr_kind_mod.o

shr_spfn_mod.o: \
	shr_const_mod.o \
	shr_kind_mod.o


clean:
	$(RM) *.f90 *.o *.mod
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.F90.o:
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS) -c $*.F90
endif
