!=================================================================================================================
 module mpas_init_atm_chems
 use mpas_log
 use mpas_kind_types
 use mpas_pool_routines

 implicit none
 private
 public:: init_atm_chems


 contains


!=================================================================================================================
 subroutine init_atm_chems(block,mesh,configs,state)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh

!inout arguments:
 type(block_type),intent(inout),target:: block
 type(mpas_pool_type),intent(inout)   :: state

!local variables and pointers:
 character (len=StrKIND),pointer:: config_start_time

 logical:: lexist

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine init_atm_chems:')

!--- lexist is initialized to true (for now):
 lexist = .true.

 if(lexist) then
    call mpas_pool_get_config(configs,'config_start_time',config_start_time)
 
    !--- horizontal interpolation:
    call init_hinterp_chems(block,mesh,state)

!   !--- vertical interpolation:
!   initial_date = trim(config_start_time)
!   call init_vinterp_gocart(initial_date,mesh,diag,state)
!else
!   call mpas_log_write('QNWFA_QNIFA_Monthly_Sigma.dat was not found in local directory:')
!   call mpas_log_write('nwfa and nifa are set to zero and not interpolated from climatological data.')
 endif

 call mpas_log_write('--- end subroutine init_atm_chems.')
!call mpas_log_write(' ')
 
 end subroutine init_atm_chems

!=================================================================================================================
 subroutine init_hinterp_chems(block,mesh,state)
!=================================================================================================================

!input arguments:
 type (mpas_pool_type),intent(in)  :: mesh

!inout arguments:
 type(block_type),intent(inout),target:: block
 type (mpas_pool_type),intent(inout)  :: state

!local variables and pointers:
 integer,pointer:: nCells,nVertLevels
 integer,pointer:: num_CO,num_chems
 integer,pointer:: CO_emissions_start,CO_emissions_end
 integer,pointer:: CO_chems_start,CO_chems_end
 integer:: iCell,k,n

 real(kind=RKIND),dimension(:,:),pointer:: CO_emissions

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine init_hinterp_chems:')

 call mpas_pool_get_dimension(mesh,'nCells',nCells)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 call mpas_pool_get_dimension(state,'num_chems',num_chems)
 call mpas_pool_get_dimension(mesh,'CO_emissions_start',CO_emissions_start)
 call mpas_pool_get_dimension(mesh,'CO_emissions_end'  ,CO_emissions_end  )
 call mpas_pool_get_dimension(state,'CO_chems_start'   ,CO_chems_start    )
 call mpas_pool_get_dimension(state,'CO_chems_end'     ,CO_chems_end      )
 call mpas_log_write('--- num_chems          = $i',intArgs=(/num_chems/))
 call mpas_log_write('--- CO_emissions_start = $i',intArgs=(/CO_emissions_start/))
 call mpas_log_write('--- CO_emissions_end   = $i',intArgs=(/CO_emissions_end/))
 call mpas_log_write('--- CO_chems_start     = $i',intArgs=(/CO_chems_start/))
 call mpas_log_write('--- CO_chems_end       = $i',intArgs=(/CO_chems_end/))

 call mpas_pool_get_array(mesh,'CO',CO_emissions)

!--- initialization:
 do n = CO_emissions_start, CO_emissions_end
    do iCell = 1, nCells
       CO_emissions(n,iCell) = 0_RKIND
    enddo
 enddo

 call mpas_log_write('--- end subroutine init_hinterp_chems:')

 end subroutine init_hinterp_chems

!=================================================================================================================
 end module mpas_init_atm_chems
!=================================================================================================================
