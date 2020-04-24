!=================================================================================================================
 module mpas_atmphys_tracers
 use mpas_derived_types
 use mpas_dmpar
 use mpas_log
 use mpas_kind_types
 use mpas_pool_routines

 use mpas_atmphys_constants,only: R_d,R_v
 use module_mynn_get_pblh

 implicit none
 private
 public:: update_tracers,  &
          tracers_get_tend


 contains


!=================================================================================================================
 subroutine update_tracers(itimestep,configs,mesh,state,time_lev,diag,diag_physics,sfc_input)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: diag
 type(mpas_pool_type),intent(in):: diag_physics
 type(mpas_pool_type),intent(in):: sfc_input

 integer,intent(in):: itimestep,time_lev

!inout arguments:
 type(mpas_pool_type),intent(inout):: state

!local variables and pointers:
 logical,pointer:: config_do_restart

 integer,pointer:: nCells,nCellsSolve,nVertLevels,nVertLevelsP1
 integer,pointer:: num_tracers,passive_start,passive_end
 integer,pointer:: index_tr17_7,index_tr17_8
 integer,pointer:: index_qv
 integer,dimension(:),pointer:: kpbl

 integer:: iCell,k
 integer:: count_pbl,kpblh

 real(kind=RKIND),pointer:: dt
 real(kind=RKIND),dimension(:),pointer:: hpbl
 real(kind=RKIND),dimension(:),pointer:: xland
 real(kind=RKIND),dimension(:,:),pointer:: exner,zgrid
 real(kind=RKIND),dimension(:,:),pointer:: qv,theta_m
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars,tracers,tracers2

 real(kind=RKIND):: rvd,rvdi
 real(kind=RKIND):: factor_decay,ht,landsea,pblh
 real(kind=RKIND),dimension(:),allocatable:: qke,sqv,th,thv,dzw,zw

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine update_tracers:')

 call mpas_pool_get_config(configs,'config_do_restart',config_do_restart)
 call mpas_pool_get_config(configs,'config_dt',dt)

 call mpas_pool_get_dimension(mesh,'nCells',nCells)
 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
 call mpas_pool_get_dimension(mesh,'nVertLevelsP1',nVertLevelsP1)

 call mpas_pool_get_array(mesh,'zgrid',zgrid)

 call mpas_pool_get_array(diag,'exner',exner)

 call mpas_pool_get_array(diag_physics,'hpbl',hpbl)
 call mpas_pool_get_array(diag_physics,'kpbl',kpbl)
 call mpas_pool_get_array(sfc_input,'xland',xland)

 call mpas_pool_get_dimension(state,'num_tracers'  ,num_tracers  )
 call mpas_pool_get_dimension(state,'passive_start',passive_start)
 call mpas_pool_get_dimension(state,'passive_end'  ,passive_end  )
 call mpas_log_write('--- num_tracers        = $i',intArgs=(/num_tracers/)  )
 call mpas_log_write('--- passive_start      = $i',intArgs=(/passive_start/))
 call mpas_log_write('--- passive_end        = $i',intArgs=(/passive_end/)  )

 call mpas_pool_get_dimension(state,'index_tr17_7',index_tr17_7)
 call mpas_pool_get_dimension(state,'index_tr17_8',index_tr17_8)
 call mpas_pool_get_array(state,'tracers',tracers ,time_lev)
 call mpas_pool_get_array(state,'tracers',tracers2,2)

 call mpas_pool_get_dimension(state,'index_qv',index_qv)
 call mpas_pool_get_array(state,'scalars',scalars,time_lev)
 call mpas_pool_get_array(state,'theta_m',theta_m,time_lev)

 rvd  = R_v/R_d
 rvdi = 1._RKIND/rvd

 if(.not. config_do_restart .and. itimestep == 1) then
    !--- for a cold start and on the first timestep, calculate the height of the PBL to set passive tracers
    !    tr17_7 and tr17_8 to zero above the PBL top.
!   call mpas_log_write('--- compute initial PBL-top height.')

    !calculate the PBL top:
    qv => scalars(index_qv,:,:)

    if(.not.allocated(dzw)) allocate(dzw(nVertLevels) )
    if(.not.allocated(qke)) allocate(qke(nVertLevels) )
    if(.not.allocated(sqv)) allocate(sqv(nVertLevels) )
    if(.not.allocated(th) ) allocate(th(nVertLevels)  )
    if(.not.allocated(thv)) allocate(thv(nVertLevels) )
    if(.not.allocated(zw) ) allocate(zw(nVertLevels+1))

    do iCell = 1, nCellsSolve
       landsea = xland(iCell)
       k = 1
       zw(k) = 0._RKIND
       do k = 1, nVertLevels
          qv(k,iCell) = max(qv(k,iCell),0._RKIND)
          sqv(k)  = qv(k,iCell)/(1._RKIND+qv(k,iCell))
          th(k)   = theta_m(k,iCell)/(1._RKIND+rvd*qv(k,iCell))
          thv(k)  = th(k)*(1._RKIND+0.61*sqv(k)) ! as in subroutine mynn_bl_driver but should use rvdi instead.
          dzw(k)  = zgrid(k+1,iCell)-zgrid(k,iCell)
          zw(k+1) = zw(k) + dzw(k)
          qke(k)  = 0.1-MIN(zw(k)*0.001,0.) !for initial PBLH calc only
       enddo
       call mynn_get_pblh(1,nVertLevels,pblh,thv,qke,zw,dzw,landsea,kpblh)
       hpbl(iCell) = pblh
       kpbl(iCell) = kpblh
    enddo

    if(allocated(dzw)) deallocate(dzw)
    if(allocated(qke)) deallocate(qke)
    if(allocated(sqv)) deallocate(sqv)
    if(allocated(th) ) deallocate(th )
    if(allocated(thv)) deallocate(thv)
    if(allocated(zw) ) deallocate(zw )

    !remove initial passive tracers tr17_7 and tr17_8 above the layer containing the PBL-top height:
    do iCell = 1, nCells
       count_pbl = 0
       k = 1
       ht = zgrid(k,iCell)
       do k = 1, nVertLevelsP1
          if((zgrid(k,iCell) - ht) .le. hpbl(iCell)) then
              count_pbl = count_pbl + 1
          endif
       enddo
       if(count_pbl .ge. 1) then
          do k = count_pbl+1, nVertLevels
             !--- time-level 1:
             tracers(index_tr17_7,k,iCell) = 0._RKIND
             tracers(index_tr17_8,k,iCell) = 0._RKIND

             !--- time-level 2:
             tracers2(index_tr17_7,k,iCell) = 0._RKIND
             tracers2(index_tr17_8,k,iCell) = 0._RKIND
          enddo
       endif
    enddo

 else

    !--- calculate factor_decay. replenish the passive tracers tr17_7 and tr17_8 below the PBL top.
    factor_decay = 1._RKIND/(86400._RKIND/dt)
    factor_decay = 1._RKIND - factor_decay

    call mpas_log_write('--- itimestep    = $i',intArgs=(/itimestep/))
    call mpas_log_write('--- factor_decay = $r',realArgs=(/factor_decay/))

    do iCell = 1, nCells
       !apply decay to scalars_chem tr17_8:
       do k = 1, nVertLevels
          tracers(index_tr17_8,k,iCell) = factor_decay*tracers(index_tr17_8,k,iCell)
       enddo

       !replenish scalars_chem tr17_7 and tr17_8 in the boundary layer:
       count_pbl = 0
       k = 1
       ht = zgrid(k,iCell)
       do k = 1, nVertLevelsP1
          if((zgrid(k,iCell) - ht) .le. hpbl(iCell)) then
              count_pbl = count_pbl + 1
          endif
       enddo
       if(count_pbl .ge. 1) then
          do k = 1, count_pbl
             tracers(index_tr17_7,k,iCell) = 1._RKIND
             tracers(index_tr17_8,k,iCell) = 1._RKIND
          enddo
       endif
    enddo
 endif

 call mpas_log_write('--- end subroutine update_tracers.')
!stop

 end subroutine update_tracers

!=================================================================================================================
 subroutine tracers_get_tend(block,configs,mesh,state,tend,tend_tracers)
!=================================================================================================================

!input variables:
 type(block_type),intent(in),target:: block
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: state

 type(mpas_pool_type),intent(in):: tend_tracers

!inout variables:
 type(mpas_pool_type),intent(inout):: tend

!local variables and pointers:
 character(len=StrKIND),pointer:: pbl_scheme,       &
                                  convection_scheme

 integer,pointer:: nCells,nCellsSolve,nVertLevels
 integer,pointer:: num_bl_tracers,num_cu_tracers,num_tracers

 real(kind=RKIND),dimension(:,:),pointer  :: mass !time level 2 rho_zz
 real(kind=RKIND),dimension(:,:,:),pointer:: bl_tracers,cu_tracers
 real(kind=RKIND),dimension(:,:,:),pointer:: tracers

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine tracers_get_tend:')

 call mpas_pool_get_config(configs,'config_pbl_scheme'       ,pbl_scheme       )
 call mpas_pool_get_config(configs,'config_convection_scheme',convection_scheme)

 call mpas_pool_get_dimension(mesh,'nCells'     ,nCells     )
 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 call mpas_pool_get_dimension(tend,'num_tracers_tend',num_tracers)
 call mpas_pool_get_dimension(tend_tracers,'num_bl_tracers',num_bl_tracers)
 call mpas_pool_get_dimension(tend_tracers,'num_cu_tracers',num_cu_tracers)
 call mpas_log_write('--- num_tracers_tend = $i',intArgs=(/num_tracers/))
 call mpas_log_write('--- num_bl_tracers   = $i',intArgs=(/num_bl_tracers/))
 call mpas_log_write('--- num_cu_tracers   = $i',intArgs=(/num_cu_tracers/))

 call mpas_pool_get_array(state,'rho_zz',mass,2)

 call mpas_pool_get_array(tend,'tracers_tend',tracers)
 call mpas_pool_get_array(tend_tracers,'bl_tracers',bl_tracers)
 call mpas_pool_get_array(tend_tracers,'cu_tracers',cu_tracers)

!in case some variables are not allocated due to their associated packages. We need to make their pointers
!associated here to avoid triggering run-time. checks when calling physics_get_tend_work:
 if(.not. associated(tracers))    allocate(tracers(0,0,0))
 if(.not. associated(bl_tracers)) allocate(bl_tracers(0,0,0))
 if(.not. associated(cu_tracers)) allocate(cu_tracers(0,0,0))


 call tracers_get_tend_work(block,nCells,nCellsSolve,nVertLevels,pbl_scheme,convection_scheme, &
                 num_bl_tracers,num_cu_tracers,num_tracers,bl_tracers,cu_tracers,tracers,mass)


!clean up any pointers that were allocated with zero size before the call to physics_get_tend_work:
 if(size(tracers) == 0)    deallocate(tracers)
 if(size(bl_tracers) == 0) deallocate(bl_tracers)
 if(size(cu_tracers) == 0) deallocate(cu_tracers)

 call mpas_log_write('--- end subroutine tracers_get_tend:')
 call mpas_log_write(' ')

 end subroutine tracers_get_tend

!=================================================================================================================
 subroutine tracers_get_tend_work(block,nCells,nCellsSolve,nVertLevels,pbl_scheme,convection_scheme, &
                       num_bl_tracers,num_cu_tracers,num_tracers,bl_tracers,cu_tracers,tracers,mass)
!=================================================================================================================

!input arguments:
 type(block_type),intent(in):: block

 character(len=StrKIND),intent(in):: pbl_scheme
 character(len=StrKIND),intent(in):: convection_scheme

 integer,intent(in):: nCells,nCellsSolve,nVertLevels
 integer,intent(in):: num_tracers,num_bl_tracers,num_cu_tracers

 real(kind=RKIND),intent(in),dimension(nVertLevels,nCells+1):: mass
 real(kind=RKIND),intent(in),dimension(num_bl_tracers,nVertLevels,nCells+1):: bl_tracers
 real(kind=RKIND),intent(in),dimension(num_cu_tracers,nVertLevels,nCells+1):: cu_tracers

!inout arguments:
 real(kind=RKIND),intent(inout),dimension(num_tracers,nVertLevels,nCells+1):: tracers

!local variables:
 integer:: i,it,k

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine tracers_get_tend_work:')

!add coupled tendencies due to PBL processes:
 if(pbl_scheme .ne. 'off') then
    do it = 1,num_tracers
       do k = 1,nVertLevels
          do i = 1,nCellsSolve
             tracers(it,k,i) = tracers(it,k,i) + mass(k,i)*bl_tracers(it,k,i)
          enddo
       enddo
    enddo
 endif

!add coupled tendencies due to subgrid-scale convection:
 if(convection_scheme .ne. 'off') then
    do it = 1,num_tracers
       do k = 1,nVertLevels
          do i = 1,nCellsSolve
             tracers(it,k,i) = tracers(it,k,i) + mass(k,i)*cu_tracers(it,k,i)
          enddo
       enddo
    enddo
 endif

!it = 1
!do i = 1,nCellsSolve
!   do k = 1,nVertLevels
!      call mpas_log_write('$i $i $r $r $r $r',intArgs=(/i,k/),realArgs=(/mass(k,i),bl_tracers(it,k,i), &
!                          cu_tracers(it,k,i),tracers(it,k,i)/))
!   enddo
!   call mpas_log_write(' ')
!enddo

 call mpas_log_write('--- end subroutine tracers_get_tend_work:')

 end subroutine tracers_get_tend_work

!=================================================================================================================

 end module mpas_atmphys_tracers
!=================================================================================================================
