!=================================================================================================================
 module mpas_atmphys_chems
 use mpas_log
 use mpas_kind_types
 use mpas_pool_routines


 implicit none
 private
 public:: chems_get_tend


 contains


!=================================================================================================================
 subroutine chems_get_tend(block,configs,mesh,state,tend,tend_chems)
!=================================================================================================================

!input variables:
 type(block_type),intent(in),target:: block
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: state

 type(mpas_pool_type),intent(in):: tend_chems

!inout variables:
 type(mpas_pool_type),intent(inout):: tend

!local variables and pointers:
 character(len=StrKIND),pointer:: pbl_scheme,       &
                                  convection_scheme

 integer,pointer:: nCells,nCellsSolve,nVertLevels
 integer,pointer:: num_bl_chems,num_cu_chems,num_chems

 real(kind=RKIND),dimension(:,:),pointer  :: mass !time level 2 rho_zz
 real(kind=RKIND),dimension(:,:,:),pointer:: bl_chems,cu_chems
 real(kind=RKIND),dimension(:,:,:),pointer:: chems

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine chemistry_get_tend:')

 call mpas_pool_get_config(configs,'config_pbl_scheme'       ,pbl_scheme       )
 call mpas_pool_get_config(configs,'config_convection_scheme',convection_scheme)

 call mpas_pool_get_dimension(mesh,'nCells'     ,nCells     )
 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 call mpas_pool_get_dimension(tend,'num_chems_tend',num_chems)
 call mpas_pool_get_dimension(tend_chems,'num_bl_chems',num_bl_chems)
 call mpas_pool_get_dimension(tend_chems,'num_cu_chems',num_cu_chems)
 call mpas_log_write('--- num_chems_tend = $i',intArgs=(/num_chems/))
 call mpas_log_write('--- num_bl_chems   = $i',intArgs=(/num_bl_chems/))
 call mpas_log_write('--- num_cu_chems   = $i',intArgs=(/num_cu_chems/))

 call mpas_pool_get_array(state,'rho_zz',mass,2)

 call mpas_pool_get_array(tend,'chems_tend',chems)
 call mpas_pool_get_array(tend_chems,'bl_chems',bl_chems)
 call mpas_pool_get_array(tend_chems,'cu_chems',cu_chems)

!in case some variables are not allocated due to their associated packages. We need to make their pointers
!associated here to avoid triggering run-time. checks when calling physics_get_tend_work:
 if(.not. associated(chems))    allocate(chems(0,0,0))
 if(.not. associated(bl_chems)) allocate(bl_chems(0,0,0))
 if(.not. associated(cu_chems)) allocate(cu_chems(0,0,0))


 call chems_get_tend_work(block,nCells,nCellsSolve,nVertLevels,pbl_scheme,convection_scheme, &
                           num_bl_chems,num_cu_chems,num_chems,bl_chems,cu_chems,chems,mass)


!clean up any pointers that were allocated with zero size before the call to physics_get_tend_work:
 if(size(chems) == 0)    deallocate(chems)
 if(size(bl_chems) == 0) deallocate(bl_chems)
 if(size(cu_chems) == 0) deallocate(cu_chems)


 call mpas_log_write('--- end subroutine chemistry_get_tend:')
!call mpas_log_write(' ')

 end subroutine chems_get_tend

!=================================================================================================================
 subroutine chems_get_tend_work(block,nCells,nCellsSolve,nVertLevels,pbl_scheme,convection_scheme, &
                                 num_bl_chems,num_cu_chems,num_chems,bl_chems,cu_chems,chems,mass)
!=================================================================================================================

!input arguments:
 type(block_type),intent(in):: block

 character(len=StrKIND),intent(in):: pbl_scheme
 character(len=StrKIND),intent(in):: convection_scheme

 integer,intent(in):: nCells,nCellsSolve,nVertLevels
 integer,intent(in):: num_chems,num_bl_chems,num_cu_chems

 real(kind=RKIND),intent(in),dimension(nVertLevels,nCells+1):: mass
 real(kind=RKIND),intent(in),dimension(num_bl_chems,nVertLevels,nCells+1):: bl_chems
 real(kind=RKIND),intent(in),dimension(num_cu_chems,nVertLevels,nCells+1):: cu_chems

!inout arguments:
 real(kind=RKIND),intent(inout),dimension(num_chems,nVertLevels,nCells+1):: chems

!local variables:
 integer:: i,ic,k

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write('--- enter subroutine chems_get_tend_work:')

!add coupled tendencies due to PBL processes:
 if(pbl_scheme .ne. 'off') then
    do ic = 1,num_chems
       do k = 1,nVertLevels
          do i = 1,nCellsSolve
             chems(ic,k,i) = chems(ic,k,i) + mass(k,i)*bl_chems(ic,k,i)
          enddo
       enddo
    enddo
 endif

!add coupled tendencies due to subgrid-scale convection:
 if(convection_scheme .ne. 'off') then
    do ic = 1,num_chems
       do k = 1,nVertLevels
          do i = 1,nCellsSolve
             chems(ic,k,i) = chems(ic,k,i) + mass(k,i)*cu_chems(ic,k,i)
          enddo
       enddo
    enddo
 endif

!call mpas_log_write('--- end subroutine chems_get_tend_work:')

 end subroutine chems_get_tend_work

!=================================================================================================================
 end module mpas_atmphys_chems
!=================================================================================================================
