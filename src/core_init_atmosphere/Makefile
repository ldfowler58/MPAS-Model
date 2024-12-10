.SUFFIXES: .F .o

ifeq "$(GOCART2G)" "true"
    CHEMISTRY = -DDO_GOCART2G
    CHEMCORE = core_gocart2G
    CHEM_OBJS = libchem/*.o
endif


OBJS_INIT = \
	mpas_init_atm_hinterp.o \
	mpas_init_atm_llxy.o \
	mpas_init_atm_read_met.o \
	read_geogrid.o

OBJS =  \
	mpas_init_atm_core.o \
	mpas_init_atm_core_interface.o \
	mpas_init_atm_cases.o \
	mpas_atm_advection.o \
	mpas_init_atm_bitarray.o \
 	mpas_init_atm_queue.o \
	mpas_init_atm_static.o \
	mpas_init_atm_gwd.o \
	mpas_gsl_oro_data_sm_scale.o \
	mpas_gsl_oro_data_lg_scale.o \
	mpas_init_atm_gwd_gsl.o \
	mpas_init_atm_surface.o \
	mpas_init_atm_vinterp.o \
	mpas_init_atm_thompson_aerosols.o \
	read_geogrid.o \
	mpas_atmphys_constants.o \
	mpas_atmphys_date_time.o \
	mpas_atmphys_functions.o \
	mpas_atmphys_initialize_real.o \
	mpas_atmphys_utilities.o \
	mpas_stack.o \
	mpas_kd_tree.o \
	mpas_parse_geoindex.o \
	mpas_geotile_manager.o

all: core_hyd_init $(CHEMCORE) core_hyd

core_hyd_init: $(OBJS_INIT)
	ar -ru libdycore.a $(OBJS_INIT)

core_hyd: $(OBJS)
	ar -ru libdycore.a $(OBJS)

core_reg:
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $(CHEMISTRY) Registry.xml > Registry_processed.xml

core_input_gen:
	if [ ! -e default_inputs ]; then  mkdir default_inputs; fi
	( cd default_inputs; $(NL_GEN) ../Registry_processed.xml namelist.init_atmosphere in_defaults=true )
	( cd default_inputs; $(ST_GEN) ../Registry_processed.xml streams.init_atmosphere stream_list.init_atmosphere. listed )

core_gocart2G:
	( cd ../core_atmosphere/chemistry; $(MAKE) -f Makefile.init_gocart2G all )

gen_includes:
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) Registry.xml > Registry_processed.xml
	(if [ ! -d inc ]; then mkdir -p inc; fi) # To generate *.inc files
	(cd inc; $(REG_PARSE) ../Registry_processed.xml $(CPPFLAGS) )

post_build:
	if [ ! -e $(ROOT_DIR)/default_inputs ]; then mkdir $(ROOT_DIR)/default_inputs; fi
	cp default_inputs/* $(ROOT_DIR)/default_inputs/.
	( cd $(ROOT_DIR)/default_inputs; for FILE in `ls -1`; do if [ ! -e ../$$FILE ]; then cp $$FILE ../.; fi; done )

mpas_init_atm_cases.o: \
	read_geogrid.o \
	mpas_atm_advection.o \
	mpas_init_atm_read_met.o \
	mpas_init_atm_llxy.o \
	mpas_init_atm_hinterp.o \
	mpas_init_atm_static.o \
	mpas_init_atm_gwd.o \
	mpas_init_atm_gwd_gsl.o \
	mpas_init_atm_surface.o \
	mpas_init_atm_thompson_aerosols.o \
	mpas_init_atm_vinterp.o \
	mpas_atmphys_constants.o \
	mpas_atmphys_functions.o \
	mpas_atmphys_initialize_real.o

mpas_init_atm_hinterp.o: mpas_init_atm_queue.o mpas_init_atm_bitarray.o

mpas_init_atm_thompson_aerosols.o: \
	mpas_init_atm_read_met.o \
	mpas_init_atm_hinterp.o \
	mpas_init_atm_llxy.o \
	mpas_init_atm_vinterp.o \
	mpas_atmphys_date_time.o \
	mpas_atmphys_utilities.o

mpas_advection.o: 

mpas_init_atm_read_met.o:

read_geogrid.o:

mpas_kd_tree.o:

mpas_init_atm_llxy.o:

mpas_init_atm_core_interface.o: mpas_init_atm_core.o

mpas_init_atm_core.o: mpas_advection.o mpas_init_atm_cases.o

mpas_stack.o:

mpas_parse_geoindex.o:

mpas_geotile_manager.o: mpas_parse_geoindex.o

mpas_init_atm_static.o: \
	mpas_atm_advection.o \
	mpas_init_atm_hinterp.o \
	mpas_init_atm_llxy.o \
	mpas_atmphys_utilities.o \
	mpas_stack.o \
	mpas_kd_tree.o \
	mpas_parse_geoindex.o \
	mpas_geotile_manager.o

mpas_init_atm_surface.o: \
	mpas_init_atm_hinterp.o  \
	mpas_init_atm_llxy.o     \
	mpas_init_atm_read_met.o

mpas_atmphys_initialize_real.o:  \
	mpas_init_atm_surface.o  \
	mpas_atmphys_date_time.o \
	mpas_atmphys_utilities.o

mpas_init_atm_gwd_gsl.o: \
	mpas_gsl_oro_data_sm_scale.o \
	mpas_gsl_oro_data_lg_scale.o

clean:
	(cd ../core_atmosphere/chemistry; $(MAKE) -f Makefile.init_gocart2G clean )
	$(RM) *.o *.mod *.f90 libdycore.a
	$(RM) Registry_processed.xml
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i
	$(RM) -r default_inputs

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) -I./inc $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../core_atmosphere/chemistry -I../framework -I../operators  -I../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I./inc -I../core_atmosphere/chemistry -I../framework -I../operators  -I../external/esmf_time_f90
endif

.c.o:
	$(CC) $(CFLAGS) $(CPPFLAGS) $(CPPINCLUDES) -c $<
