!=================================================================================================================
 module mpas_chemistry_init
 use mpas_kind_types
 use mpas_log,only          : mpas_log_write
 use mpas_derived_types,only: mpas_pool_type
 use mpas_pool_routines,only: mpas_pool_get_config

 use mpas_chemistry_opp,only: init_opp


 implicit none
 private
 public:: chemistry_init


!initialization and interpolation to MPAS mesh of input data needed to run the CHEM2D Ozone Photochemistry
!Parameterization (CHEM2D-OPP) of McCormack et al. (2006).


 contains


!=================================================================================================================
 subroutine chemistry_init(configs,mesh,chem_input)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh


!inout arguments:
 type(mpas_pool_type),intent(inout):: chem_input


!local variables and pointers:
 logical,pointer:: o3prognostic

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine chemistry_init:')

 call mpas_pool_get_config(configs,'config_o3prognostic',o3prognostic)


!--- initializes the input data need to run the CHEM2D-OPP of McCormack et al. (2006):
 if(o3prognostic) then
    call init_opp(configs,mesh,chem_input)
 endif


 call mpas_log_write('--- end subroutine chemistry_init:')
 call mpas_log_write(' ')

 end subroutine chemistry_init

!=================================================================================================================
 end module mpas_chemistry_init
!=================================================================================================================
