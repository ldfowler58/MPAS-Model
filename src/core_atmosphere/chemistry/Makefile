.SUFFIXES: .F .o

ifeq ($(CORE),atmosphere)
COREDEF = -Dmpas
endif

all: dummy core_chemistry

dummy:
	echo "****** compiling chemistry ******"

OBJS = \
	module_mynn_get_pblh.o \
	mpas_atmphys_advection.o \
	mpas_atmphys_chems.o \
	mpas_atmphys_chems_physics.o \
	mpas_atmphys_tracers.o \
	mpas_atmphys_tracers_physics.o

core_chemistry: $(OBJS)
	ar -ru libchem.a $(OBJS)

# DEPENDENCIES:
mpas_atmphys_tracers.o: \
	module_mynn_get_pblh.o

clean:
	$(RM) *.o *.mod *.f90 libchem.a
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(COREDEF) $(HYDROSTATIC) $(PHYSICS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I. -I.. -I../dynamics -I../physics -I../../framework -I../../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(HYDROSTATIC) $(PHYSICS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I. -I.. -I../dynamics -I../physics -I../../framework -I../../external/esmf_time_f90
endif
