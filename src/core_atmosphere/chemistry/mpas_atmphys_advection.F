!=================================================================================================================
 module mpas_atmphys_advection
 use mpas_log
 use mpas_kind_types
 use mpas_pool_routines
 use mpas_timer

 use mpas_atm_boundaries,only: nRelaxZone

 implicit none
 private

 public:: atmphys_advection


 contains


!=================================================================================================================
 subroutine atmphys_advection(domain,rk_timestep)
!=================================================================================================================

!input arguments:
 real(kind=RKIND),intent(in),dimension(3):: rk_timestep

!inout arguments:
 type(domain_type),intent(inout):: domain

!local variables and arrays:
 type(block_type),pointer:: block
 type(mpas_pool_type),pointer:: mesh
 type(mpas_pool_type),pointer:: diag
 type(mpas_pool_type),pointer:: state
 type(mpas_pool_type),pointer:: tend

 logical,pointer:: scalar_advection,         &
                   split_dynamics_transport, &
                   apply_lbcs,               &
                   monotonic,                &
                   positive_definite

 integer,pointer:: nThreads
 integer,pointer:: nCells,nCellsSolve,nEdges,nVertLevels
 integer,pointer:: time_integration_order
 integer,pointer:: num_scalars
 integer,dimension(:),pointer:: cellThreadEnd,cellThreadStart,              &
                                cellSolveThreadEnd,cellSolveThreadStart,    &
                                edgeThreadEnd,edgeThreadStart,              &
                                edgeSolveThreadEnd,edgeSolveThreadStart,    &
                                vertexThreadEnd,vertexThreadStart,          &
                                vertexSolveThreadEnd,vertexSolveThreadStart
 integer,dimension(:),pointer:: nAdvCellsForEdge,nEdgesOnCell
 integer,dimension(:),pointer:: bdyMaskCell,bdyMaskEdge
 integer,dimension(:,:),pointer:: advCellsForEdge,cellsOnEdge,edgesOnCell

 integer:: i,k,ithread,rk_step

 real(kind=RKIND),dimension(:,:),allocatable:: scalar_old_arr,     &
                                               scalar_new_arr,     &
                                               s_max_arr,          &
                                               s_min_arr,          &
                                               flux_array,         &
                                               wdtn_arr,           &
                                               rho_zz_int,         &
                                               flux_tmp_arr,       &
                                               flux_upwind_tmp_arr

 real(kind=RKIND),dimension(:,:,:),allocatable:: scale_array,       &
                                                 scalar_tend_array, &
                                                 horiz_flux_arr

 real(kind=RKIND),pointer:: dt
 real(kind=RKIND),pointer:: coef_3rd_order
 real(kind=RKIND),dimension(:),pointer:: invAreaCell
 real(kind=RKIND),dimension(:),pointer:: fnm,fnp,rdnw,meshScalingDel2,meshScalingDel4
 real(kind=RKIND),dimension(:),pointer:: dvEdge, qv_init
 real(kind=RKIND),dimension(:,:),pointer:: adv_coefs,adv_coefs_3rd,edgesOnCell_sign
 real(kind=RKIND),dimension(:,:),pointer:: rho_edge,zgrid
 real(kind=RKIND),dimension(:,:),pointer:: rho_zz_old,rho_zz_new
 real(kind=RKIND),dimension(:,:),pointer:: kdiff,uhAvg,wwAvg

 real(kind=RKIND),dimension(:,:,:),pointer:: deriv_two
 real(kind=RKIND),dimension(:,:,:),pointer:: scalar_old,scalar_new,scalar_tend_save

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
!call mpas_log_write('--- enter subroutine atmphys_advection:')

 RK3_SPLIT_TRANSPORT : do rk_step = 1,3
    call mpas_log_write('--- rk_step = $i',intArgs=(/rk_step/))
    block => domain%blocklist
    do while(associated(block))
       call mpas_pool_get_config(block%configs,'config_scalar_advection'        ,scalar_advection        )
       call mpas_pool_get_config(block%configs,'config_split_dynamics_transport',split_dynamics_transport)
       call mpas_pool_get_config(block%configs,'config_time_integration_order'  ,time_integration_order  )
       call mpas_pool_get_config(block%configs,'config_apply_lbcs'              ,apply_lbcs              )
       call mpas_pool_get_config(block%configs,'config_monotonic'               ,monotonic               )
       call mpas_pool_get_config(block%configs,'config_positive_definite'       ,positive_definite       )
!      call mpas_log_write('--- config_monotonic = $l',logicArgs=(/monotonic/))
!      call mpas_log_write('--- config_positive_definite = $l',logicArgs=(/positive_definite/))

       if(rk_step < 3 .or. (.not. monotonic .and. .not. positive_definite)) then
          call mpas_timer_start('atm_advance_chemistry')
       else
          call mpas_timer_start('atm_advance_chemistry')
       endif

       call mpas_pool_get_config(block%configs,'config_dt',dt)
       call mpas_pool_get_config(block%configs,'config_coef_3rd_order',coef_3rd_order)

       call mpas_pool_get_subpool(block%structs,'tend' ,tend)
       call mpas_pool_get_subpool(block%structs,'state',state)
       call mpas_pool_get_subpool(block%structs,'diag' ,diag)
       call mpas_pool_get_subpool(block%structs,'mesh' ,mesh)

       call mpas_pool_get_dimension(mesh ,'nCells'     ,nCells     )
       call mpas_pool_get_dimension(mesh ,'nCellsSolve',nCellsSolve)
       call mpas_pool_get_dimension(mesh ,'nEdges'     ,nEdges     )
       call mpas_pool_get_dimension(mesh ,'nVertLevels',nVertLevels)

       call mpas_pool_get_dimension(state,'num_scalars',num_scalars)

       call mpas_pool_get_dimension(block%dimensions,'nThreads',nThreads)
       call mpas_pool_get_dimension(block%dimensions,'cellThreadStart'     ,cellThreadStart     )
       call mpas_pool_get_dimension(block%dimensions,'cellThreadEnd'       ,cellThreadEnd       )
       call mpas_pool_get_dimension(block%dimensions,'cellSolveThreadStart',cellSolveThreadStart)
       call mpas_pool_get_dimension(block%dimensions,'cellSolveThreadEnd'  ,cellSolveThreadEnd  )

       call mpas_pool_get_dimension(block%dimensions,'vertexThreadStart'     ,vertexThreadStart     )
       call mpas_pool_get_dimension(block%dimensions,'vertexThreadEnd'       ,vertexThreadEnd       )
       call mpas_pool_get_dimension(block%dimensions,'vertexSolveThreadStart',vertexSolveThreadStart)
       call mpas_pool_get_dimension(block%dimensions,'vertexSolveThreadEnd'  ,vertexSolveThreadEnd  )

       call mpas_pool_get_dimension(block%dimensions,'edgeThreadStart'     ,edgeThreadStart     )
       call mpas_pool_get_dimension(block%dimensions,'edgeThreadEnd'       ,edgeThreadEnd       )
       call mpas_pool_get_dimension(block%dimensions,'edgeSolveThreadStart',edgeSolveThreadStart)
       call mpas_pool_get_dimension(block%dimensions,'edgeSolveThreadEnd'  ,edgeSolveThreadEnd  )

       call mpas_pool_get_array(mesh,'nEdgesOnCell'    ,nEdgesOnCell    )
       call mpas_pool_get_array(mesh,'nAdvCellsForEdge',nAdvCellsForEdge)
       call mpas_pool_get_array(mesh,'cellsOnEdge'     ,cellsOnEdge     )
       call mpas_pool_get_array(mesh,'edgesOnCell'     ,edgesOnCell     )
       call mpas_pool_get_array(mesh,'advCellsForEdge' ,advCellsForEdge )

       call mpas_pool_get_array(mesh,'bdyMaskCell',bdyMaskCell)
       call mpas_pool_get_array(mesh,'bdyMaskEdge',bdyMaskEdge)

       call mpas_pool_get_array(mesh,'fzm' ,fnm )
       call mpas_pool_get_array(mesh,'fzp' ,fnp )
       call mpas_pool_get_array(mesh,'rdzw',rdnw)
       call mpas_pool_get_array(mesh,'adv_coefs'       ,adv_coefs       )
       call mpas_pool_get_array(mesh,'adv_coefs_3rd'   ,adv_coefs_3rd   )
       call mpas_pool_get_array(mesh,'edgesOnCell_sign',edgesOnCell_sign)
       call mpas_pool_get_array(mesh,'deriv_two'       ,deriv_two       )
       call mpas_pool_get_array(mesh,'dvEdge'          ,dvEdge          )
       call mpas_pool_get_array(mesh,'qv_init'         ,qv_init         )
       call mpas_pool_get_array(mesh,'zgrid'           ,zgrid           )
       call mpas_pool_get_array(mesh,'invAreaCell'     ,invAreaCell     )
       call mpas_pool_get_array(mesh,'meshScalingDel2' ,meshScalingDel2 )
       call mpas_pool_get_array(mesh,'meshScalingDel4' ,meshScalingDel4 )

       call mpas_pool_get_array(state,'rho_zz' ,rho_zz_old,1)
       call mpas_pool_get_array(state,'rho_zz' ,rho_zz_new,2)
       call mpas_pool_get_array(state,'scalars',scalar_old,1)
       call mpas_pool_get_array(state,'scalars',scalar_new,2)

       call mpas_pool_get_array(diag,'rho_edge',rho_edge)
       call mpas_pool_get_array(diag,'kdiff'   ,kdiff   )
       call mpas_pool_get_array(diag,'ruAvg'   ,uhAvg   )
       call mpas_pool_get_array(diag,'wwAvg'   ,wwAvg   )

       call mpas_pool_get_array(tend,'scalars_tend',scalar_tend_save)

!--- allocate local arrays:
       allocate(scalar_old_arr(nVertLevels,nCells+1))
       allocate(scalar_new_arr(nVertLevels,nCells+1))
       allocate(s_max_arr(nVertLevels,nCells+1))
       allocate(s_min_arr(nVertLevels,nCells+1))
       allocate(scale_array(nVertLevels,2,nCells+1))
       allocate(flux_array(nVertLevels,nEdges+1))
       allocate(wdtn_arr(nVertLevels+1,nCells+1))
       allocate(rho_zz_int(nVertLevels,nCells+1))
       allocate(scalar_tend_array(num_scalars,nVertLevels,nCells+1))

       scalar_old_arr(:,nCells+1)      = 0.0_RKIND
       scalar_new_arr(:,nCells+1)      = 0.0_RKIND
       s_max_arr(:,nCells+1)           = 0.0_RKIND
       s_min_arr(:,nCells+1)           = 0.0_RKIND
       scale_array(:,:,nCells+1)       = 0.0_RKIND
       flux_array(:,nEdges+1)          = 0.0_RKIND
       wdtn_arr(:,nCells+1)            = 0.0_RKIND
       rho_zz_int(:,nCells+1)          = 0.0_RKIND
       scalar_tend_array(:,:,nCells+1) = 0.0_RKIND

       if(rk_step < 3 .or. (.not. monotonic .and. .not. positive_definite)) then
          allocate(horiz_flux_arr(num_scalars,nVertLevels,nEdges+1))
          horiz_flux_arr(:,:,nEdges+1) = 0.0_RKIND
       else
          allocate(flux_upwind_tmp_arr(nVertLevels,nEdges+1))
          flux_upwind_tmp_arr(:,nEdges+1) = 0.0_RKIND
          allocate(flux_tmp_arr(nVertLevels,nEdges+1))
          flux_tmp_arr(:,nEdges+1) = 0.0_RKIND
       endif

       call mpas_log_write('--- before subroutine atm_advance_chemistry:')
       do i = 1,1
          do k = 1,nVertLevels
             call mpas_log_write('$i $i $r $r',intArgs=(/i,k/),realArgs= &
                           (/scalar_old(1,k,i),scalar_new(1,k,i)/))
          enddo
          call mpas_log_write(' ')
       enddo

!---
!$OMP PARALLEL DO
      do ithread = 1,nThreads
         if(rk_step < 3 .or. (.not. monotonic .and. .not. positive_definite)) then
            call atm_advance_chemistry(num_scalars,nCells,nCellsSolve,nEdges,nVertLevels,  &
                        cellThreadStart(ithread),cellThreadEnd(ithread),                   &
                        cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread),         &
                        edgeThreadStart(ithread),edgeThreadEnd(ithread),                   &
                        edgeSolveThreadStart(ithread),edgeSolveThreadEnd(ithread),         &
                        vertexThreadStart(ithread),vertexThreadEnd(ithread),               &
                        vertexSolveThreadStart(ithread),vertexSolveThreadEnd(ithread),     &
                        nEdgesOnCell,nAdvCellsForEdge,bdyMaskCell,bdyMaskEdge,cellsOnEdge, &
                        edgesOnCell,advCellsForEdge,rk_step,time_integration_order,        &
                        apply_lbcs,fnm,fnp,rdnw,meshScalingDel2,meshScalingDel4,           &
                        adv_coefs, adv_coefs_3rd,qv_init,rho_edge,zgrid,invAreaCell,       &
                        edgesOnCell_sign,dvEdge,coef_3rd_order,rk_timestep(rk_step),kdiff, &
                        uhAvg,wwAvg,deriv_two,rho_zz_old,rho_zz_new,scalar_old,scalar_new, &
                        scalar_tend_save,horiz_flux_arr,advance_density=.true.)
         else
            call mpas_log_write('--- rk_step = $i - do not call atm_advance_chemistry:',intArgs=(/rk_step/))
         endif
      enddo
!$OMP END PARALLEL DO

       call mpas_log_write('--- after subroutine atm_advance_chemistry:')
       do i = 1,1
          do k = 1,nVertLevels
             call mpas_log_write('$i $i $r $r',intArgs=(/i,k/),realArgs= &
                           (/scalar_old(1,k,i),scalar_new(1,k,i)/))
          enddo
          call mpas_log_write(' ')
       enddo

!--- deallocate local arrays:
       deallocate(scalar_old_arr)
       deallocate(scalar_new_arr)
       deallocate(s_max_arr)
       deallocate(s_min_arr)
       deallocate(scale_array)
       deallocate(flux_array)
       deallocate(wdtn_arr)
       deallocate(rho_zz_int)
       deallocate(scalar_tend_array)

       if(rk_step < 3 .or. (.not. monotonic .and. .not. positive_definite)) then
          deallocate(horiz_flux_arr)
       else
          deallocate(flux_upwind_tmp_arr)
          deallocate(flux_tmp_arr)
       endif

       if(rk_step < 3 .or. (.not. monotonic .and. .not. positive_definite)) then
          call mpas_timer_stop('atm_advance_chemistry')
       else
          call mpas_timer_stop('atm_advance_chemistry')
       endif
       call mpas_log_write('--- end subroutine atmphys_advection:')
       call mpas_log_write(' ')
    block => block % next
    end do

 enddo RK3_SPLIT_TRANSPORT

!call mpas_log_write('--- end subroutine atmphys_advection:')
!call mpas_log_write(' ')

 end subroutine atmphys_advection

!=================================================================================================================
 subroutine atm_advance_chemistry(num_scalars,nCells,nCellsSolve,nEdges,nVertLevels,     &
                   cellStart,cellEnd,cellSolveStart,cellSolveEnd,                        &
                   edgeStart,edgeEnd,edgeSolveStart,edgeSolveEnd,                        &
                   vertexStart,vertexEnd,vertexSolveStart,vertexSolveEnd,                &
                   nEdgesOnCell,nAdvCellsForEdge,bdyMaskCell,bdyMaskEdge,cellsOnEdge,    &
                   edgesOnCell, advCellsForEdge,rk_step,time_integration_order,          &
                   apply_lbcs,fnm,fnp,rdnw,meshScalingDel2,meshScalingDel4,adv_coefs,    &
                   adv_coefs_3rd,qv_init,rho_edge,zgrid,invAreaCell,edgesOnCell_sign,    &
                   dvEdge,coef_3rd_order,dt,kdiff,uhAvg,wwAvg,deriv_two,rho_zz_old,      &
                   rho_zz_new,scalar_old,scalar_new,scalar_tend_save,horiz_flux_arr,     &
                   rho_zz_int,scalar_tend,advance_density)
!=================================================================================================================

!input arguments:
 logical,intent(in):: apply_lbcs
 logical,intent(in),optional:: advance_density

 integer,intent(in):: num_scalars
 integer,intent(in):: nCells,nCellsSolve,nEdges,nVertLevels
 integer,intent(in):: cellEnd,cellStart,cellSolveEnd,cellSolveStart,        &
                      edgeEnd,edgeStart,edgeSolveEnd,edgeSolveStart,        &
                      vertexEnd,vertexStart,vertexSolveEnd,vertexSolveStart
 integer,intent(in):: rk_step,time_integration_order

 integer,intent(in),dimension(:):: nEdgesOnCell,nAdvCellsForEdge
 integer,intent(in),dimension(:):: bdyMaskCell,bdyMaskEdge
 integer,intent(in),dimension(:,:):: cellsOnEdge,edgesOnCell
 integer,intent(in),dimension(:,:):: advCellsForEdge

 real(kind=RKIND),intent(in):: coef_3rd_order,dt

 real(kind=RKIND),intent(in),dimension(:):: invAreaCell
 real(kind=RKIND),intent(in),dimension(:):: dvEdge
 real(kind=RKIND),intent(in),dimension(:):: fnm,fnp,rdnw 
 real(kind=RKIND),intent(in),dimension(:):: qv_init
 real(kind=RKIND),intent(in),dimension(:):: meshScalingDel2,meshScalingDel4

 real(kind=RKIND),intent(in),dimension(:,:):: adv_coefs,adv_coefs_3rd
 real(kind=RKIND),intent(in),dimension(:,:):: rho_edge,zgrid
 real(kind=RKIND),intent(in),dimension(:,:):: edgesOnCell_sign
 real(kind=RKIND),intent(in),dimension(:,:):: kdiff,uhAvg,wwAvg
 real(kind=RKIND),intent(in),dimension(:,:):: rho_zz_old,rho_zz_new

 real(kind=RKIND),intent(in),dimension(:,:,:):: deriv_two
 real(kind=RKIND),intent(in),dimension(:,:,:):: scalar_old

!inout arguments:
 real(kind=RKIND),intent(inout),dimension(:,:,:):: scalar_new
 real(kind=RKIND),intent(inout),dimension(num_scalars,nVertLevels,nCells+1):: scalar_tend_save
 real(kind=RKIND),intent(inout),dimension(num_scalars,nVertLevels,nEdges+1):: horiz_flux_arr

 real(kind=RKIND),intent(inout),dimension(nVertLevels,nCells+1),optional:: rho_zz_int
 real(kind=RKIND),intent(inout),dimension(num_scalars,nVertLevels,nCells+1),optional :: scalar_tend

!local variables and arrays:
 logical:: local_advance_density

 integer:: iAdvCell,iCell,iEdge,iScalar,cell1,cell2,i,j,k
 integer,dimension(10):: ica

 real(kind=RKIND):: u_direction,u_negative,u_positive
 real(kind=RKIND):: rho_zz_new_inv,scalar_weight,weight_time_old,weight_time_new
 real(kind=RKIND),dimension(nVertLevels,10):: scalar_weight2
 real(kind=RKIND),dimension(num_scalars,nVertLevels):: scalar_tend_column! local storage to accumulate tendencies
 real(kind=RKIND),dimension(num_scalars,nVertLevels+1 ):: wdtn

 real(kind=RKIND):: flux3,flux4
 real(kind=RKIND):: q_im2,q_im1,q_i,q_ip1,ua,coef3
 flux4(q_im2,q_im1,q_i,q_ip1,ua) = ua*(7.*(q_i+q_im1) - (q_ip1+q_im2))/12.0
 flux3(q_im2,q_im1,q_i,q_ip1,ua,coef3) = flux4(q_im2,q_im1,q_i,q_ip1, ua) &
                                       + coef3*abs(ua)*((q_ip1-q_im2)-3.*(q_i-q_im1))/12.0


!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine atm_advance_chemistry:')
!call mpas_log_write('--- num_scalars      = $i',intArgs=(/num_scalars/))
!call mpas_log_write('--- nCells           = $i',intArgs=(/nCells/))
!call mpas_log_write('--- nCellsSolve      = $i',intArgs=(/nCellsSolve/))
!call mpas_log_write('--- nEdges           = $i',intArgs=(/nEdges/))
!call mpas_log_write('--- nVertLevels      = $i',intArgs=(/nVertLevels/))

 if(present(advance_density)) then
    local_advance_density = advance_density
 else
    local_advance_density = .true.
 endif

!weights for the time interpolation of the input density:
!
 if(.not. advance_density) then
    weight_time_new = 1.
 else
    if((rk_step == 1) .and. time_integration_order == 3) weight_time_new = 1./3
    if((rk_step == 1) .and. time_integration_order == 2) weight_time_new = 1./2
    if(rk_step == 2) weight_time_new = 1./2
    if(rk_step == 3) weight_time_new = 1.
 endif
 weight_time_old = 1. - weight_time_new
!call mpas_log_write('--- local_advance_density = $l',logicArgs=(/local_advance_density/))
!call mpas_log_write('--- rk_step    = $i',intArgs=(/rk_step/))
!call mpas_log_write('--- nRelaxZone = $i',intArgs=(/nRelaxZone/))
!call mpas_log_write('--- weight_time_new = $r',realArgs=(/weight_time_new/))
!call mpas_log_write('--- weight_time_old = $r',realArgs=(/weight_time_old/))

 do iEdge = edgeStart,edgeEnd
    if((.not.apply_lbcs) .or. (bdyMaskEdge(iEdge).lt.nRelaxZone-1)) then  ! full flux calculation

       select case(nAdvCellsForEdge(iEdge))
          case(10)
             do j = 1,10
!DIR$ IVDEP
                do k = 1,nVertLevels
                   scalar_weight2(k,j) = adv_coefs(j,iEdge) + sign(1.0_RKIND,uhAvg(k,iEdge))*adv_coefs_3rd(j,iEdge)
               enddo
             enddo
             do j = 1,10
                ica(j) = advCellsForEdge(j,iEdge)
             enddo
!DIR$ IVDEP
             do k = 1,nVertLevels
!DIR$ IVDEP
                do iScalar = 1,num_scalars
                   horiz_flux_arr(iScalar,k,iEdge) = &
                        scalar_weight2(k,1)  * scalar_new(iScalar,k,ica(1)) + &
                        scalar_weight2(k,2)  * scalar_new(iScalar,k,ica(2)) + &
                        scalar_weight2(k,3)  * scalar_new(iScalar,k,ica(3)) + &
                        scalar_weight2(k,4)  * scalar_new(iScalar,k,ica(4)) + &
                        scalar_weight2(k,5)  * scalar_new(iScalar,k,ica(5)) + &
                        scalar_weight2(k,6)  * scalar_new(iScalar,k,ica(6)) + &
                        scalar_weight2(k,7)  * scalar_new(iScalar,k,ica(7)) + &
                        scalar_weight2(k,8)  * scalar_new(iScalar,k,ica(8)) + &
                        scalar_weight2(k,9)  * scalar_new(iScalar,k,ica(9)) + &
                        scalar_weight2(k,10) * scalar_new(iScalar,k,ica(10))
                enddo
             enddo

          case default

             horiz_flux_arr(:,:,iEdge) = 0.0
             do j = 1,nAdvCellsForEdge(iEdge)
                iAdvCell = advCellsForEdge(j,iEdge)
!DIR$ IVDEP
                do k = 1, nVertLevels
                   scalar_weight = adv_coefs(j,iEdge) + sign(1.0_RKIND,uhAvg(k,iEdge))*adv_coefs_3rd(j,iEdge)
!DIR$ IVDEP
                   do iScalar = 1, num_scalars
                      horiz_flux_arr(iScalar,k,iEdge) = horiz_flux_arr(iScalar,k,iEdge) &
                                                      + scalar_weight*scalar_new(iScalar,k,iAdvCell)
                   enddo
                   
                enddo
             enddo

       end select

!      do k = 1,nVertLevels
!         call mpas_log_write('$i $i $r',intArgs=(/iEdge,k/),realArgs=(/horiz_flux_arr(1,k,iEdge)/))
!      enddo

    elseif(apply_lbcs .and. (bdyMaskEdge(iEdge).ge.nRelaxZone-1) .and. (bdyMaskEdge(iEdge).le.nRelaxZone)) then

       !upwind flux evaluation for outermost 2 edges in specified zone
       cell1 = cellsOnEdge(1,iEdge)
       cell2 = cellsOnEdge(2,iEdge)
!DIR$ IVDEP
       do k = 1,nVertLevels
          u_direction = sign(0.5_RKIND,uhAvg(k,iEdge))
          u_positive = dvEdge(iEdge)*abs(u_direction + 0.5_RKIND)
          u_negative = dvEdge(iEdge)*abs(u_direction - 0.5_RKIND)
!DIR$ IVDEP
          do iScalar = 1,num_scalars
             horiz_flux_arr(iScalar,k,iEdge) = u_positive*scalar_new(iScalar,k,cell1) &
                                             + u_negative*scalar_new(iScalar,k,cell2)
          enddo
       enddo

    endif ! end of regional MPAS test
 enddo

!$OMP BARRIER
!scalar update,for each column sum fluxes over horizontal edges, add physics tendency, and add vertical flux
!divergence in update.

 do iCell = cellSolveStart,cellSolveEnd
    if(bdyMaskCell(iCell) <= nRelaxZone) then !specified zone for regional_MPAS is not updated in this routine

#ifndef DO_PHYSICS
       scalar_tend_save(:,:,iCell) = 0.0  !  testing purposes - we have no sources or sinks
#endif
       scalar_tend_column(1:num_scalars,1:nVertlevels) = 0.

       do i = 1,nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i,iCell)

          !here we add the horizontal flux divergence into the scalar tendency.
          !note that the scalar tendency is modified.
!DIR$ IVDEP
          do k = 1,nVertLevels
!DIR$ IVDEP
             do iScalar = 1,num_scalars
                scalar_tend_column(iScalar,k) = scalar_tend_column(iScalar,k) &
                      - edgesOnCell_sign(i,iCell) * uhAvg(k,iEdge)*horiz_flux_arr(iScalar,k,iEdge)
             enddo
          enddo
       enddo
!      do k = 1, nVertLevels
!         call mpas_log_write('$i $i $r $r',intArgs=(/iCell,k/),realArgs=(/scalar_tend_column(1,k),scalar_tend_save(1,k,iCell)/))
!      enddo
!      call mpas_log_write(' ')

!DIR$ IVDEP
       do k = 1,nVertLevels
!DIR$ IVDEP
          do iScalar = 1,num_scalars
             scalar_tend_column(iScalar,k) = scalar_tend_column(iScalar,k)*invAreaCell(iCell) &
                                           + scalar_tend_save(iScalar,k,iCell)
          enddo
       enddo
!      do k = 1, nVertLevels
!         call mpas_log_write('$i $i $r $r',intArgs=(/iCell,k/),realArgs=(/scalar_tend_column(1,k),scalar_tend_save(1,k,iCell)/))
!      enddo
!      call mpas_log_write(' ')

!
!vertical flux divergence and update of the scalars
!
       wdtn(:,1) = 0.0
       wdtn(:,nVertLevels+1) = 0.0

       k = 2
       do iScalar = 1,num_scalars
          wdtn(iScalar,k) = wwAvg(k,iCell)*(fnm(k)*scalar_new(iScalar,k,iCell) &
                          + fnp(k)*scalar_new(iScalar,k-1,iCell))
       enddo
!DIR$ IVDEP
       do k = 3,nVertLevels-1
!DIR$ IVDEP
          do iScalar =1,num_scalars
             wdtn(iScalar,k) = flux3(scalar_new(iScalar,k-2,iCell),scalar_new(iScalar,k-1,iCell),  &
                                     scalar_new(iScalar,k  ,iCell),scalar_new(iScalar,k+1,iCell),  &
                                     wwAvg(k,iCell),coef_3rd_order )
          enddo
       enddo
       k = nVertLevels
       do iScalar = 1, num_scalars
          wdtn(iScalar,k) = wwAvg(k,iCell)*(fnm(k)*scalar_new(iScalar,k,iCell) &
                          + fnp(k)*scalar_new(iScalar,k-1,iCell))
       enddo

!      do k = 1,nVertLevels
!         call mpas_log_write('$i $i $r $r',intArgs=(/iCell,k/),realArgs=(/wdtn(1,k),scalar_new(1,k,iCell)/))
!      enddo

!DIR$ IVDEP
       do k = 1,nVertLevels
          rho_zz_new_inv = 1.0_RKIND / (weight_time_old*rho_zz_old(k,iCell) &
                         + weight_time_new*rho_zz_new(k,iCell))
!DIR$ IVDEP
          do iScalar =1,num_scalars
             scalar_new(iScalar,k,iCell) = (scalar_old(iScalar,k,iCell)*rho_zz_old(k,iCell) &
                                         + dt*( scalar_tend_column(iScalar,k)               &
                                         - rdnw(k)*(wdtn(iScalar,k+1)-wdtn(iScalar,k))))    &
                                         * rho_zz_new_inv
          enddo
!         call mpas_log_write('$i $i $r $r $r $r $r $r',intArgs=(/iCell,k/),realArgs=(/dt,rho_zz_new_inv, &
!                             scalar_tend_column(1,k),wdtn(1,k),scalar_old(1,k,iCell),scalar_new(1,k,iCell)/))
       enddo

    endif ! specified zone regional_MPAS test

 enddo

!call mpas_log_write('--- end subroutine atm_advance_chemistry:')
!call mpas_log_write(' ')

 end subroutine atm_advance_chemistry

!=================================================================================================================
 end module mpas_atmphys_advection
!=================================================================================================================

