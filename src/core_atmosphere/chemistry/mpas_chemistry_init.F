!=================================================================================================================
 module mpas_chemistry_init
 use mpas_kind_types
 use mpas_log,only          : mpas_log_write
 use mpas_derived_types,only: mpas_pool_type
 use mpas_pool_routines,only: mpas_pool_get_config

 use mpas_init_chem2DOPP,only: init_chem2DOPP

 implicit none
 private
 public:: chemistry_init2DOPP


!initialization and interpolation to MPAS mesh of input data needed to run the CHEM2D Ozone Photochemistry
!Parameterization (CHEM2D-OPP) of McCormack et al. (2006).


 contains


!=================================================================================================================
 subroutine chemistry_init2DOPP(mesh,chemopp_input)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh


!inout arguments:
 type(mpas_pool_type),intent(inout):: chemopp_input

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine chemistry_init2DOPP:')
 
 call init_chem2DOPP(mesh,chemopp_input)

 call mpas_log_write('--- end subroutine chemistry_init2DOPP:')
 call mpas_log_write(' ')

 end subroutine chemistry_init2DOPP

!=================================================================================================================
 end module mpas_chemistry_init
!=================================================================================================================
