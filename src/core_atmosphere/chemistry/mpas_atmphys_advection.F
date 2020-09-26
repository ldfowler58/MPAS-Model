!=================================================================================================================
 module mpas_atmphys_advection
 use mpas_log
 use mpas_kind_types
 use mpas_dmpar
 use mpas_pool_routines
 use mpas_timer

 use mpas_atm_boundaries,only: nRelaxZone,nSpecZone


 implicit none
 private

 public:: atmphys_advection, &
          atmphys_advection_intro, &
          atmphys_rk_integration_setup


 contains


!=================================================================================================================
 subroutine atmphys_advection_intro(domain,rk_timestep,string_in)
!=================================================================================================================

!input arguments:
 character(len=*),intent(in):: string_in
 real(kind=RKIND),intent(in),dimension(3):: rk_timestep

!inout arguments:
 type(domain_type),intent(inout):: domain

!local variables and arrays:
 type(block_type),pointer:: block
 type(mpas_pool_type),pointer:: state
 type(mpas_pool_type),pointer:: tend

 character(len=StrKind):: num_string,tend_string

 integer,pointer:: num_species
 real(kind=RKIND),dimension(:,:,:),pointer:: species_old,species_new
 real(kind=RKIND),dimension(:,:,:),pointer:: species_tend

 integer:: iCell

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine atmphys_advection_intro:')

 block => domain%blocklist
 do while(associated(block))
    num_string  = 'num_'//trim(string_in)
    tend_string = trim(string_in)//'_tend' 

    call mpas_pool_get_subpool(block%structs,'state',state)
    call mpas_pool_get_subpool(block%structs,'tend' ,tend )

    call mpas_pool_get_dimension(state,trim(num_string),num_species)
    call mpas_log_write('--- num_species = $i',intArgs=(/num_species/))

    call mpas_pool_get_array(state,trim(string_in),species_old,1)
    call mpas_pool_get_array(state,trim(string_in),species_new,2)

    call mpas_pool_get_array(tend,trim(tend_string),species_tend)

    block => block % next
 end do

 call mpas_log_write('--- end subroutine atmphys_advection_intro:')

 end subroutine atmphys_advection_intro

!=================================================================================================================
 subroutine atmphys_rk_integration_setup(domain)
!=================================================================================================================

!inout arguments:
 type(domain_type),intent(inout):: domain

!local variables and arrays:
 type(block_type),pointer:: block
 type(mpas_pool_type),pointer:: mesh
 type(mpas_pool_type),pointer:: state

 integer,pointer:: nThreads
 integer,dimension(:),pointer:: cellThreadEnd,cellThreadStart
 integer:: ithread
 integer:: cellStart,cellEnd

 real(kind=RKIND),dimension(:,:,:),pointer:: chems_1,chems_2
 real(kind=RKIND),dimension(:,:,:),pointer:: tracers_1,tracers_2

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter atmphys_rk_integration_setup:')

 block => domain%blocklist
 do while(associated(block))

    call mpas_pool_get_subpool(block%structs,'state',state)

    call mpas_pool_get_dimension(block%dimensions,'nThreads',nThreads)
    call mpas_pool_get_dimension(block%dimensions,'cellThreadStart'     ,cellThreadStart     )
    call mpas_pool_get_dimension(block%dimensions,'cellThreadEnd'       ,cellThreadEnd       )

!$OMP PARALLEL DO
    do ithread=1,nThreads
       cellStart = cellThreadStart(ithread)
       cellEnd   = cellThreadEnd(ithread)

       call mpas_pool_get_array(state,'chems',chems_1,1)
       call mpas_pool_get_array(state,'chems',chems_2,2)
       call mpas_pool_get_array(state,'tracers',tracers_1,1)
       call mpas_pool_get_array(state,'tracers',tracers_2,2)

       chems_2(:,:,cellStart:cellEnd)   = chems_1(:,:,cellStart:cellEnd)
       tracers_2(:,:,cellStart:cellEnd) = tracers_1(:,:,cellStart:cellEnd)
    enddo
!$OMP END PARALLEL DO

 block => block%next
 enddo

 call mpas_log_write('--- end atmphys_rk_integration_setup')

 end subroutine atmphys_rk_integration_setup

!=================================================================================================================
 subroutine atmphys_advection(domain,rk_timestep,string_in)
!=================================================================================================================

!input arguments:
 character(len=*),intent(in):: string_in
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

 character(len=StrKind):: num_string,tend_string

 integer,pointer:: nThreads
 integer,pointer:: nCells,nCellsSolve,nEdges,nVertLevels
 integer,pointer:: time_integration_order
!integer,pointer:: num_scalars
 integer,pointer:: num_species
 integer,dimension(:),pointer:: cellThreadEnd,cellThreadStart,           &
                                cellSolveThreadEnd,cellSolveThreadStart, &
                                edgeThreadEnd,edgeThreadStart,           &
                                edgeSolveThreadEnd,edgeSolveThreadStart
 integer,dimension(:),pointer:: nAdvCellsForEdge,nEdgesOnCell
 integer,dimension(:),pointer:: bdyMaskCell,bdyMaskEdge
 integer,dimension(:,:),pointer:: advCellsForEdge,cellsOnCell,cellsOnEdge,edgesOnCell

 integer:: ithread,rk_step
 integer:: i,k,n

 real(kind=RKIND),pointer:: dt
 real(kind=RKIND),pointer:: coef_3rd_order
 real(kind=RKIND),dimension(:),pointer:: invAreaCell
 real(kind=RKIND),dimension(:),pointer:: fnm,fnp,rdnw
 real(kind=RKIND),dimension(:),pointer:: dvEdge
 real(kind=RKIND),dimension(:,:),pointer:: adv_coefs,adv_coefs_3rd,edgesOnCell_sign
 real(kind=RKIND),dimension(:,:),pointer:: rho_zz_old,rho_zz_new
 real(kind=RKIND),dimension(:,:),pointer:: uhAvg,wwAvg

!real(kind=RKIND),dimension(:,:,:),pointer:: scalar_old,scalar_new,scalar_tend

!type(field3DReal),pointer:: scalar_old_field
!type(field3DReal),pointer:: scalar_field

 real(kind=RKIND),dimension(:,:,:),pointer:: species_old,species_new
 real(kind=RKIND),dimension(:,:,:),pointer:: species_tend
 type(field3DReal),pointer:: species_old_field
 type(field3DReal),pointer:: species_field

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine atmphys_advection:')

 num_string  = 'num_'//trim(string_in)
 tend_string = trim(string_in)//'_tend'

 RK3_SPLIT_TRANSPORT : do rk_step = 1,3
    block => domain%blocklist
    do while(associated(block))
       call mpas_pool_get_config(block%configs,'config_scalar_advection'        ,scalar_advection        )
       call mpas_pool_get_config(block%configs,'config_split_dynamics_transport',split_dynamics_transport)
       call mpas_pool_get_config(block%configs,'config_time_integration_order'  ,time_integration_order  )
       call mpas_pool_get_config(block%configs,'config_apply_lbcs'              ,apply_lbcs              )
       call mpas_pool_get_config(block%configs,'config_monotonic'               ,monotonic               )
       call mpas_pool_get_config(block%configs,'config_positive_definite'       ,positive_definite       )

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

!      call mpas_pool_get_dimension(state,'num_scalars',num_scalars)

       call mpas_pool_get_dimension(block%dimensions,'nThreads',nThreads)
       call mpas_pool_get_dimension(block%dimensions,'cellThreadStart'     ,cellThreadStart     )
       call mpas_pool_get_dimension(block%dimensions,'cellThreadEnd'       ,cellThreadEnd       )
       call mpas_pool_get_dimension(block%dimensions,'cellSolveThreadStart',cellSolveThreadStart)
       call mpas_pool_get_dimension(block%dimensions,'cellSolveThreadEnd'  ,cellSolveThreadEnd  )

       call mpas_pool_get_dimension(block%dimensions,'edgeThreadStart'     ,edgeThreadStart     )
       call mpas_pool_get_dimension(block%dimensions,'edgeThreadEnd'       ,edgeThreadEnd       )
       call mpas_pool_get_dimension(block%dimensions,'edgeSolveThreadStart',edgeSolveThreadStart)
       call mpas_pool_get_dimension(block%dimensions,'edgeSolveThreadEnd'  ,edgeSolveThreadEnd  )

       call mpas_pool_get_array(mesh,'nEdgesOnCell'    ,nEdgesOnCell    )
       call mpas_pool_get_array(mesh,'nAdvCellsForEdge',nAdvCellsForEdge)
       call mpas_pool_get_array(mesh,'cellsOnCell'     ,cellsOnCell     )
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
       call mpas_pool_get_array(mesh,'dvEdge'          ,dvEdge          )
       call mpas_pool_get_array(mesh,'invAreaCell'     ,invAreaCell     )

       call mpas_pool_get_array(state,'rho_zz' ,rho_zz_old,1)
       call mpas_pool_get_array(state,'rho_zz' ,rho_zz_new,2)
!      call mpas_pool_get_array(state,'scalars',scalar_old,1)
!      call mpas_pool_get_array(state,'scalars',scalar_new,2)

       call mpas_pool_get_array(diag,'ruAvg',uhAvg)
       call mpas_pool_get_array(diag,'wwAvg',wwAvg)

!      call mpas_pool_get_array(tend,'scalars_tend',scalar_tend)

       call mpas_pool_get_dimension(state,trim(num_string),num_species)
       call mpas_pool_get_array(state,trim(string_in),species_old,1 )
       call mpas_pool_get_array(state,trim(string_in),species_new,2 )
       call mpas_pool_get_array(tend ,trim(tend_string),species_tend)

!---
!$OMP PARALLEL DO
       do ithread = 1,nThreads
          call mpas_log_write('--- rk_step = $i',intArgs=(/rk_step/))

          if(rk_step < 3 .or. (.not. monotonic .and. .not. positive_definite)) then

             call atm_advance_chemistry(num_species,nCells,nEdges,nVertLevels,            &
                      cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread),          &
                      edgeThreadStart(ithread),edgeThreadEnd(ithread),                    &
                      nEdgesOnCell,nAdvCellsForEdge,bdyMaskCell,bdyMaskEdge,cellsOnEdge,  &
                      edgesOnCell,advCellsForEdge,rk_step,time_integration_order,         &
                      apply_lbcs,fnm,fnp,rdnw,adv_coefs,adv_coefs_3rd,invAreaCell,        &
                      edgesOnCell_sign,dvEdge,coef_3rd_order,rk_timestep(rk_step),        &
                      uhAvg,wwAvg,rho_zz_old,rho_zz_new,species_old,species_new,          &
                      species_tend,advance_density=.true.)

          else

             call mpas_pool_get_field(state,trim(string_in),species_old_field,1)
             call atm_advance_chemistry_mono(block,num_species,nCells,nCellsSolve,        &
                      nEdges,nVertLevels,cellThreadStart(ithread),cellThreadEnd(ithread), &
                      cellSolveThreadStart(ithread),cellSolveThreadEnd(ithread),          &
                      edgeThreadStart(ithread),edgeThreadEnd(ithread),                    &
                      nEdgesOnCell,nAdvCellsForEdge,bdyMaskCell,bdyMaskEdge,cellsOnCell,  &
                      cellsOnEdge,edgesOnCell,advCellsForEdge,adv_coefs,adv_coefs_3rd,    &
                      apply_lbcs,rk_timestep(rk_step),fnm,fnp,edgesOnCell_sign,dvEdge,    &
                      coef_3rd_order,invAreaCell,rdnw,uhAvg,wwAvg,rho_zz_old,rho_zz_new,  &
                      species_old,species_new,species_tend,species_old_field,             &
                      advance_density=.true.)

          endif
       enddo
!$OMP END PARALLEL DO

       if(rk_step < 3 .or. (.not. monotonic .and. .not. positive_definite)) then
          call mpas_timer_stop('atm_advance_chemistry')
       else
          call mpas_timer_stop('atm_advance_chemistry')
       endif

    block => block%next
    enddo

    if(rk_step < 3) then
       call mpas_pool_get_field(state,trim(string_in),species_field,2)
       call mpas_dmpar_exch_halo_field(species_field)
    endif

 enddo RK3_SPLIT_TRANSPORT

 call mpas_log_write('--- end subroutine atmphys_advection:')

 end subroutine atmphys_advection

!=================================================================================================================
 subroutine atm_advance_chemistry(num_scalars,nCells,nEdges,nVertLevels,           &
                cellSolveStart,cellSolveEnd,edgeStart,edgeEnd,                     &
                nEdgesOnCell,nAdvCellsForEdge,bdyMaskCell,bdyMaskEdge,cellsOnEdge, &
                edgesOnCell,advCellsForEdge,rk_step,time_integration_order,        &
                apply_lbcs,fnm,fnp,rdnw,adv_coefs,adv_coefs_3rd,invAreaCell,       &
                edgesOnCell_sign,dvEdge,coef_3rd_order,dt,uhAvg,wwAvg,rho_zz_old,  &
                rho_zz_new,scalar_old,scalar_new,scalar_tend,advance_density)
!=================================================================================================================

!input arguments:
 logical,intent(in):: apply_lbcs
 logical,intent(in),optional:: advance_density

 integer,intent(in):: num_scalars
 integer,intent(in):: nCells,nEdges,nVertLevels
 integer,intent(in):: cellSolveEnd,cellSolveStart,edgeEnd,edgeStart
 integer,intent(in):: rk_step,time_integration_order

 integer,intent(in),dimension(:):: nEdgesOnCell,nAdvCellsForEdge
 integer,intent(in),dimension(:):: bdyMaskCell,bdyMaskEdge
 integer,intent(in),dimension(:,:):: cellsOnEdge,edgesOnCell
 integer,intent(in),dimension(:,:):: advCellsForEdge

 real(kind=RKIND),intent(in):: coef_3rd_order,dt

 real(kind=RKIND),intent(in),dimension(:):: invAreaCell
 real(kind=RKIND),intent(in),dimension(:):: dvEdge
 real(kind=RKIND),intent(in),dimension(:):: fnm,fnp,rdnw 

 real(kind=RKIND),intent(in),dimension(:,:):: adv_coefs,adv_coefs_3rd
 real(kind=RKIND),intent(in),dimension(:,:):: edgesOnCell_sign
 real(kind=RKIND),intent(in),dimension(:,:):: uhAvg,wwAvg
 real(kind=RKIND),intent(in),dimension(:,:):: rho_zz_old,rho_zz_new

 real(kind=RKIND),intent(in),dimension(:,:,:):: scalar_old
 real(kind=RKIND),intent(in),dimension(num_scalars,nVertLevels,nCells+1):: scalar_tend

!inout arguments:
 real(kind=RKIND),intent(inout),dimension(:,:,:):: scalar_new

!local variables and arrays:
 logical:: local_advance_density

 integer:: iAdvCell,iCell,iEdge,iScalar,cell1,cell2,i,j,k
 integer,dimension(10):: ica

 real(kind=RKIND):: u_direction,u_negative,u_positive
 real(kind=RKIND):: rho_zz_new_inv,scalar_weight,weight_time_old,weight_time_new
 real(kind=RKIND),dimension(nVertLevels,10):: scalar_weight2
 real(kind=RKIND),dimension(num_scalars,nVertLevels  ):: scalar_tend_column! local storage to accumulate tendencies
 real(kind=RKIND),dimension(num_scalars,nVertLevels+1):: wdtn
 real(kind=RKIND),dimension(num_scalars,nVertLevels,nEdges+1):: horiz_flux_arr

 real(kind=RKIND):: flux3,flux4
 real(kind=RKIND):: q_im2,q_im1,q_i,q_ip1,ua,coef3
 flux4(q_im2,q_im1,q_i,q_ip1,ua) = ua*(7.*(q_i+q_im1) - (q_ip1+q_im2))/12.0
 flux3(q_im2,q_im1,q_i,q_ip1,ua,coef3) = flux4(q_im2,q_im1,q_i,q_ip1, ua) &
                                       + coef3*abs(ua)*((q_ip1-q_im2)-3.*(q_i-q_im1))/12.0

!-----------------------------------------------------------------------------------------------------------------
 if(present(advance_density)) then
    local_advance_density = advance_density
 else
    local_advance_density = .true.
 endif

 do iScalar = 1,num_scalars
    do k = 1,nVertLevels
       horiz_flux_arr(iScalar,k,nEdges+1) = 0._RKIND
    enddo
 enddo

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

!DIR$ IVDEP
       do k = 1,nVertLevels
!DIR$ IVDEP
          do iScalar = 1,num_scalars
             scalar_tend_column(iScalar,k) = scalar_tend_column(iScalar,k)*invAreaCell(iCell) &
                                           + scalar_tend(iScalar,k,iCell)
          enddo
       enddo

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
       enddo

    endif ! specified zone regional_MPAS test
 enddo

 end subroutine atm_advance_chemistry

!=================================================================================================================
 subroutine atm_advance_chemistry_mono(block,num_scalars,nCells,nCellsSolve,nEdges,         &
                nVertLevels,cellStart,cellEnd,cellSolveStart,cellSolveEnd,edgeStart,        &
                edgeEnd,nEdgesOnCell,nAdvCellsForEdge,bdyMaskCell,bdyMaskEdge,cellsOnCell,  &
                cellsOnEdge,edgesOnCell,advCellsForEdge,adv_coefs,adv_coefs_3rd,apply_lbcs, &
                dt,fnm,fnp,edgesOnCell_sign,dvEdge,coef_3rd_order,invAreaCell,rdnw,uhAvg,   &
                wwAvg,rho_zz_old,rho_zz_new,scalar_old,scalar_new,scalar_tend,              &
                scalar_old_field,advance_density)
!=================================================================================================================

!input arguments:
 logical,intent(in):: apply_lbcs
 logical,intent(in),optional:: advance_density

 integer,intent(in):: num_scalars
 integer,intent(in):: nCells,nCellsSolve,nEdges,nVertLevels
 integer,intent(in):: cellStart,cellEnd,cellSolveEnd,cellSolveStart
 integer,intent(in):: edgeStart,edgeEnd

 integer,intent(in),dimension(:):: nEdgesOnCell,nAdvCellsForEdge
 integer,intent(in),dimension(:):: bdyMaskCell,bdyMaskEdge
 integer,intent(in),dimension(:,:):: cellsOnCell,cellsOnEdge,edgesOnCell
 integer,intent(in),dimension(:,:):: advCellsForEdge

 real(kind=RKIND),intent(in):: coef_3rd_order,dt

 real(kind=RKIND),intent(in),dimension(:):: invAreaCell
 real(kind=RKIND),intent(in),dimension(:):: dvEdge
 real(kind=RKIND),intent(in),dimension(:):: fnm,fnp,rdnw

 real(kind=RKIND),intent(in),dimension(:,:):: adv_coefs,adv_coefs_3rd
 real(kind=RKIND),intent(in),dimension(:,:):: edgesOnCell_sign
 real(kind=RKIND),intent(in),dimension(:,:):: uhAvg,wwAvg
 real(kind=RKIND),intent(in),dimension(:,:):: rho_zz_old,rho_zz_new
!real(kind=RKIND),intent(in),dimension(nVertLevels,nCells+1):: rho_zz_old,rho_zz_new

!inout arguments:
 type(block_type),intent(inout),target:: block
 type(field3DReal),intent(inout),pointer:: scalar_old_field

 real(kind=RKIND),intent(inout),dimension(num_scalars,nVertLevels,nCells+1):: scalar_old,scalar_new
 real(kind=RKIND),intent(inout),dimension(num_scalars,nVertLevels,nCells+1):: scalar_tend

!local variables and arrays:
 logical:: local_advance_density

 integer,parameter:: &
    scale_in  = 1, &
    scale_out = 2

 integer:: iCell,iEdge,iScalar,i,ii,k,kk,jj
 integer:: cell1,cell2
 integer,dimension(10):: ica

 real(kind=RKIND),parameter:: eps=1.e-20
!real(kind=RKIND),parameter:: eps=1.e-15

 real(kind=RKIND):: scale_factor,scalar_weight
 real(kind=RKIND):: scmax,scmin
 real(kind=RKIND):: flux,flux_upwind
 real(kind=RKIND),dimension(10,2):: swa
 real(kind=RKIND),dimension(nVertLevels):: flux_upwind_arr
 real(kind=RKIND),dimension(nVertLevels,nCells+1)  :: s_min,s_max
 real(kind=RKIND),dimension(nVertLevels,nCells+1)  :: rho_zz_int

 real(kind=RKIND),dimension(nVertLevels,nCells+1)  :: temp_old,temp_new
 real(kind=RKIND),dimension(nVertLevels,2,nCells+1),target:: scale_arr
 real(kind=RKIND),dimension(nVertLevels+1,nCells+1):: wdtn

 real(kind=RKIND),dimension(nVertLevels,nEdges+1):: flux_arr,flux_upwind_tmp, flux_tmp

 real(kind=RKIND):: flux3,flux4
 real(kind=RKIND):: q_im2,q_im1,q_i,q_ip1,ua,coef3
 flux4(q_im2,q_im1,q_i,q_ip1,ua) = ua*(7.*(q_i+q_im1) - (q_ip1+q_im2))/12.0
 flux3(q_im2,q_im1,q_i,q_ip1,ua,coef3) = flux4(q_im2,q_im1,q_i,q_ip1, ua) &
                                       + coef3*abs(ua)*((q_ip1-q_im2)-3.*(q_i-q_im1))/12.0

 real(kind=RKIND),dimension(:,:,:),pointer:: scalar_old_arr
 type(field3DReal),pointer:: tempField
 type(field3DReal),target :: tempFieldTarget

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write('--- enter subroutine atm_advance_chemistry_mono: dt = $r',realArgs=(/dt/))

 if(present(advance_density)) then
    local_advance_density = advance_density
 else
    local_advance_density = .true.
 endif

!initialization:
 do k = 1,nVertLevels
    do iCell = 1,nCells+1
       s_min(k,iCell)      = 0._RKIND
       s_max(k,iCell)      = 0._RKIND
       temp_old(k,iCell)   = 0._RKIND
       temp_new(k,iCell)   = 0._RKIND
       rho_zz_int(k,iCell) = 0._RKIND
       do kk = 1,2
          scale_arr(k,kk,iCell) = 0._RKIND
       enddo
    enddo
 enddo
 do k = 1,nVertLevels+1
    do iCell = 1,nCells+1
       wdtn(k,iCell) = 0._RKIND
    enddo
 enddo

 do k = 1,nVertLevels
    do iEdge = 1,nEdges+1
       flux_arr(k,iEdge) = 0._RKIND
       flux_tmp(k,iEdge) = 0._RKIND
       flux_upwind_tmp(k,iEdge) = 0._RKIND
    enddo
 enddo

!
!for positive-definite or monotonic option, we first update scalars using the tendency from sources other than
!the resolved transport (these should constitute a positive definite update).  
!Note, however, that we enforce positive-definiteness in this update.
!The transport will maintain this positive definite solution and optionally, shape preservation (monotonicity).
!
 do iCell = cellSolveStart,cellSolveEnd
!DIR$ IVDEP
    do k = 1,nVertLevels
!DIR$ IVDEP
       do iScalar = 1,num_scalars
#ifndef DO_PHYSICS
          scalar_tend(iScalar,k,iCell) = 0.0  !  testing purposes - we have no sources or sinks
#endif
          scalar_old(iScalar,k,iCell) = scalar_old(iScalar,k,iCell) &
                                      + dt*scalar_tend(iScalar,k,iCell)/rho_zz_old(k,iCell)
          scalar_tend(iScalar,k,iCell) = 0.0
       enddo
    enddo
 enddo

!$OMP BARRIER
!$OMP MASTER
 call mpas_dmpar_exch_halo_field(scalar_old_field)
 scalar_old_arr => scalar_old_field%array
!enddo
!$OMP END MASTER
!$OMP BARRIER

!
!Runge Kutta integration, so we compute fluxes from scalar_new values, update starts from scalar_old
!
 if(local_advance_density) then
    !begin with update of density
    do iCell = cellStart,cellEnd
       rho_zz_int(:,iCell) = 0.0
    end do
!$OMP BARRIER
    do iCell = cellSolveStart,cellSolveEnd
       do i = 1,nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i,iCell)
!DIR$ IVDEP
          do k = 1,nVertLevels
             rho_zz_int(k,iCell) = rho_zz_int(k,iCell) - edgesOnCell_sign(i,iCell) * uhAvg(k,iEdge) &
                                 * dvEdge(iEdge) * invAreaCell(iCell)
          enddo
       enddo
    enddo
    do iCell = cellSolveStart,cellSolveEnd
!DIR$ IVDEP
       do k = 1,nVertLevels
          rho_zz_int(k,iCell) = rho_zz_old(k,iCell) + dt*( rho_zz_int(k,iCell) &
                              - rdnw(k)*(wwAvg(k+1,iCell)-wwAvg(k,iCell)) )
       enddo
    enddo
!$OMP BARRIER
 endif

!
!next, do one scalar at a time:
!
 do iScalar = 1,num_scalars
    do iCell = cellStart,cellEnd
!DIR$ IVDEP
       do k = 1,nVertLevels
          temp_old(k,iCell) = scalar_old(iScalar,k,iCell)
          temp_new(k,iCell) = scalar_new(iScalar,k,iCell)
       enddo
    enddo
    do k = 1,nVertLevels
       temp_old(k,nCells+1) = 0.
       temp_new(k,nCells+1) = 0.
    enddo

!$OMP BARRIER
#ifdef DEBUG_TRANSPORT
    scmin = temp_old(1,1)
    scmax = temp_old(1,1)
    do iCell = cellStart,cellEnd
       do k = 1,nVertLevels
          scmin = min(scmin,temp_old(k,iCell))
          scmax = max(scmax,temp_old(k,iCell))
       enddo
    enddo
    call mpas_log_write(' scmin, scmin old in $r $r', realArgs=(/scmin,scmax/))

    scmin = temp_new(1,1)
    scmax = temp_new(1,1)
    do iCell = cellStart,cellEnd
       do k = 1,nVertLevels
          scmin = min(scmin,temp_new(k,iCell))
          scmax = max(scmax,temp_new(k,iCell))
       enddo
    enddo
    call mpas_log_write(' scmin, scmin new in $r $r', realArgs=(/scmin,scmax/))
#endif

    !
    !vertical flux divergence, and min and max bounds for flux limiter:
    !
    do iCell = cellSolveStart,cellSolveEnd

       !zero flux at top and bottom:
       wdtn(1,iCell) = 0.0
       wdtn(nVertLevels+1,iCell) = 0.0

       k = 1
       s_max(k,iCell) = max(temp_old(1,iCell),temp_old(2,iCell))
       s_min(k,iCell) = min(temp_old(1,iCell),temp_old(2,iCell))

       k = 2
       wdtn(k,iCell) = wwAvg(k,iCell)*(fnm(k)*temp_new(k,iCell)+fnp(k)*temp_new(k-1,iCell))
       s_max(k,iCell) = max(temp_old(k-1,iCell),temp_old(k,iCell),temp_old(k+1,iCell))
       s_min(k,iCell) = min(temp_old(k-1,iCell),temp_old(k,iCell),temp_old(k+1,iCell))

!DIR$ IVDEP
       do k = 3,nVertLevels-1
          wdtn(k,iCell) = flux3(temp_new(k-2,iCell),temp_new(k-1,iCell), &
                                temp_new(k  ,iCell),temp_new(k+1,iCell), &
                                wwAvg(k,iCell),coef_3rd_order)
          s_max(k,iCell) = max(temp_old(k-1,iCell),temp_old(k,iCell),temp_old(k+1,iCell))
          s_min(k,iCell) = min(temp_old(k-1,iCell),temp_old(k,iCell),temp_old(k+1,iCell))
       enddo

       k = nVertLevels
       wdtn(k,iCell) = wwAvg(k,iCell)*(fnm(k)*temp_new(k,iCell)+fnp(k)*temp_new(k-1,iCell))
       s_max(k,iCell) = max(temp_old(k,iCell),temp_old(k-1,iCell))
       s_min(k,iCell) = min(temp_old(k,iCell),temp_old(k-1,iCell))

       !
       !pull s_min and s_max from the (horizontal) surrounding cells:
       !

       !speclal treatment of calculations involving hexagonal cells
       !original code retained in select "default" case
       select case(nEdgesOnCell(iCell))
          case(6)
!DIR$ IVDEP
             do k = 1,nVertLevels
                s_max(k,iCell) = max(s_max(k,iCell), &
                      temp_old(k,cellsOnCell(1,iCell)), &
                      temp_old(k,cellsOnCell(2,iCell)), &
                      temp_old(k,cellsOnCell(3,iCell)), &
                      temp_old(k,cellsOnCell(4,iCell)), &
                      temp_old(k,cellsOnCell(5,iCell)), &
                      temp_old(k,cellsOnCell(6,iCell)))
                s_min(k,iCell) = min(s_min(k,iCell), &
                      temp_old(k,cellsOnCell(1,iCell)), &
                      temp_old(k,cellsOnCell(2,iCell)), &
                      temp_old(k,cellsOnCell(3,iCell)), &
                      temp_old(k,cellsOnCell(4,iCell)), &
                      temp_old(k,cellsOnCell(5,iCell)), &
                      temp_old(k,cellsOnCell(6,iCell)))
             enddo

          case default
             do i = 1,nEdgesOnCell(iCell)
!DIR$ IVDEP
                do k = 1,nVertLevels
                   s_max(k,iCell) = max(s_max(k,iCell),temp_old(k,cellsOnCell(i,iCell)))
                   s_min(k,iCell) = min(s_min(k,iCell),temp_old(k,cellsOnCell(i,iCell)))
                enddo
             enddo
       end select
    enddo

    !
    !horizontal flux divergence:
    !
    do iEdge = edgeStart,edgeEnd

       cell1 = cellsOnEdge(1,iEdge)
       cell2 = cellsOnEdge(2,iEdge)

       if(cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then  ! only for owned cells

          !special treatment of calculations involving edges between hexagonal cells
          !original code retained in select "default" case
          !be sure to see additional declarations near top of subroutine
          select case(nAdvCellsForEdge(iEdge))
             case(10)
                do jj = 1,10
                   ica(jj)   = advCellsForEdge(jj,iEdge)
                   swa(jj,1) = adv_coefs(jj,iEdge) + adv_coefs_3rd(jj,iEdge)
                   swa(jj,2) = adv_coefs(jj,iEdge) - adv_coefs_3rd(jj,iEdge)
                enddo
                do k = 1,nVertLevels
                   ii = merge(1,2,uhAvg(k,iEdge) > 0)
                   flux_arr(k,iEdge) = uhAvg(k,iEdge)*( &
                        swa(1,ii)*temp_new(k,ica(1)) + swa(2,ii)*temp_new(k,ica(2)) + &
                        swa(3,ii)*temp_new(k,ica(3)) + swa(4,ii)*temp_new(k,ica(4)) + &
                        swa(5,ii)*temp_new(k,ica(5)) + swa(6,ii)*temp_new(k,ica(6)) + &
                        swa(7,ii)*temp_new(k,ica(7)) + swa(8,ii)*temp_new(k,ica(8)) + &
                        swa(9,ii)*temp_new(k,ica(9)) + swa(10,ii)*temp_new(k,ica(10)))
                enddo

             case default
                do k = 1,nVertLevels
                   flux_arr(k,iEdge) = 0.0_RKIND
                enddo
                do i =1,nAdvCellsForEdge(iEdge)
                   iCell = advCellsForEdge(i,iEdge)
!DIR$ IVDEP
                   do k = 1,nVertLevels
                      scalar_weight = uhAvg(k,iEdge)*(adv_coefs(i,iEdge) &
                                    + sign(1.0_RKIND,uhAvg(k,iEdge))*adv_coefs_3rd(i,iEdge))
                      flux_arr(k,iEdge) = flux_arr(k,iEdge) + scalar_weight*temp_new(k,iCell)
                   enddo
                enddo
          end select
       else
          flux_arr(:,iEdge) = 0.0_RKIND
       endif
    enddo

!$OMP BARRIER
    !
    !vertical flux divergence for upwind update, we will put upwind update into scalar_new,
    !and put factor of dt in fluxes:
    !
    do iCell = cellSolveStart,cellSolveEnd
       k = 1
       temp_new(k,iCell) = temp_old(k,iCell)*rho_zz_old(k,iCell)
!DIR$ IVDEP
       do k = 2,nVertLevels
          temp_new(k,iCell)  = temp_old(k,iCell)*rho_zz_old(k,iCell)
          flux_upwind_arr(k) = dt*(max(0.0_RKIND,wwAvg(k,iCell))*temp_old(k-1,iCell) &
                             + min(0.0_RKIND,wwAvg(k,iCell))*temp_old(k,iCell))
       enddo
       do k = 1,nVertLevels-1
          temp_new(k,iCell) = temp_new(k,iCell) - flux_upwind_arr(k+1)*rdnw(k)
       enddo
!DIR$ IVDEP
       do k = 2,nVertLevels
          temp_new(k,iCell) = temp_new(k  ,iCell) + flux_upwind_arr(k)*rdnw(k)
          wdtn(k,iCell) = dt*wdtn(k,iCell) - flux_upwind_arr(k)
       enddo
       !
       !scale_arr(SCALE_IN,:,:) and scale_arr(SCALE_OUT:,:) are used here to store the incoming and
       !outgoing perturbation flux contributions to the update:  first the vertical flux component,
       !then the horizontal.
       !
!DIR$ IVDEP
       do k = 1,nVertLevels
          scale_arr(k,SCALE_IN,iCell)  = - rdnw(k)*(min(0.0_RKIND,wdtn(k+1,iCell))-max(0.0_RKIND,wdtn(k,iCell)))
          scale_arr(k,SCALE_OUT,iCell) = - rdnw(k)*(max(0.0_RKIND,wdtn(k+1,iCell))-min(0.0_RKIND,wdtn(k,iCell)))
       enddo
    enddo

    !
    !horizontal flux divergence for upwind update:
    !
    !upwind flux computation
    do iEdge = edgeStart,edgeEnd
       cell1 = cellsOnEdge(1,iEdge)
       cell2 = cellsOnEdge(2,iEdge)
!DIR$ IVDEP
       do k =1,nVertLevels
          flux_upwind_tmp(k,iEdge) = dvEdge(iEdge)*dt & 
                                   * (max(0.0_RKIND,uhAvg(k,iEdge))*temp_old(k,cell1) &
                                   +  min(0.0_RKIND,uhAvg(k,iEdge))*temp_old(k,cell2))
          flux_tmp(k,iEdge) = dt*flux_arr(k,iEdge) - flux_upwind_tmp(k,iEdge)
       enddo

       if(apply_lbcs .and. (bdyMaskEdge(iEdge) == nRelaxZone) .or. &
          (bdyMaskEdge(iEdge) == nRelaxZone-1) ) then
           flux_tmp(:,iEdge) = 0.
           flux_arr(:,iEdge) = flux_upwind_tmp(:,iEdge)
       endif
    enddo

!$OMP BARRIER
    do iCell = cellSolveStart,cellSolveEnd
       do i = 1,nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i,iCell)
!DIR$ IVDEP
          do k = 1,nVertLevels
             temp_new(k,iCell) = temp_new(k,iCell) - edgesOnCell_sign(i,iCell) &
                               * flux_upwind_tmp(k,iEdge)*invAreaCell(iCell)

             scale_arr(k,SCALE_OUT,iCell) = scale_arr(k,scale_out,iCell) &
                               - max(0.0_RKIND,edgesOnCell_sign(i,iCell)*flux_tmp(k,iEdge))*invAreaCell(iCell)
             scale_arr(k,SCALE_IN, iCell) = scale_arr(k,scale_in, iCell) &
                               - min(0.0_RKIND,edgesOnCell_sign(i,iCell)*flux_tmp(k,iEdge))*invAreaCell(iCell)
          enddo
       enddo
    enddo

    !
    ! next, the limiter:
    !

    ! simplification of limiter calculations
    ! worked through algebra and found equivalent form
    ! added benefit that it should address ifort single prec overflow issue
    if(local_advance_density) then
       do iCell = cellSolveStart,cellSolveEnd
!DIR$ IVDEP
          do k = 1,nVertLevels
             scale_factor = (s_max(k,iCell)*rho_zz_int(k,iCell) - temp_new(k,iCell)) &
                          / (scale_arr(k,SCALE_IN,iCell)  + eps)
             scale_arr(k,SCALE_IN,iCell) = min( 1.0_RKIND, max( 0.0_RKIND, scale_factor))

             scale_factor = (s_min(k,iCell)*rho_zz_int(k,iCell) - temp_new(k,iCell)) &
                          / (scale_arr(k,SCALE_OUT,iCell) - eps)
             scale_arr(k,SCALE_OUT,iCell) = min( 1.0_RKIND, max( 0.0_RKIND, scale_factor) )
          enddo
       enddo
    else
       do iCell = cellSolveStart,cellSolveEnd
!DIR$ IVDEP
          do k = 1,nVertLevels
             scale_factor = (s_max(k,iCell)*rho_zz_new(k,iCell) - temp_new(k,iCell)) &
                          / (scale_arr(k,SCALE_IN,iCell)  + eps)
             scale_arr(k,SCALE_IN,iCell) = min( 1.0_RKIND, max( 0.0_RKIND, scale_factor) )

             scale_factor = (s_min(k,iCell)*rho_zz_new(k,iCell) - temp_new(k,iCell)) &
                          / (scale_arr(k,SCALE_OUT,iCell) - eps)
             scale_arr(k,SCALE_OUT,iCell) = min( 1.0_RKIND, max( 0.0_RKIND, scale_factor) )
          enddo
       enddo
    endif

    !
    ! communicate scale factors here:
    ! communicate only first halo row in these next two exchanges:
    !
!$OMP BARRIER
!$OMP MASTER
    tempField => tempFieldTarget

    tempField % block => block
    tempField % dimSizes(1) = nVertLevels
    tempField % dimSizes(2) = 2
    tempField % dimSizes(3) = nCells
    tempField % sendList => block % parinfo % cellsToSend
    tempField % recvList => block % parinfo % cellsToRecv
    tempField % copyList => block % parinfo % cellsToCopy
    tempField % prev => null()
    tempField % next => null()
    tempField % isActive = .true.

    tempField % array => scale_arr
    call mpas_dmpar_exch_halo_field(tempField, (/ 1 /))
!$OMP END MASTER
!$OMP BARRIER

    do iEdge = edgeStart,edgeEnd
       cell1 = cellsOnEdge(1,iEdge)
       cell2 = cellsOnEdge(2,iEdge)
       if(cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then  ! only for owned cells
!DIR$ IVDEP
          do k = 1,nVertLevels
             flux_upwind = dvEdge(iEdge)*dt &
                         * (max(0.0_RKIND,uhAvg(k,iEdge))*temp_old(k,cell1) &
                         +  min(0.0_RKIND,uhAvg(k,iEdge))*temp_old(k,cell2))
             flux_arr(k,iEdge) = dt*flux_arr(k,iEdge) - flux_upwind
          enddo

          if(apply_lbcs .and. (bdyMaskEdge(iEdge) == nRelaxZone) .or. &
             (bdyMaskEdge(iEdge) == nRelaxZone-1) ) then
             flux_arr(:,iEdge) = 0.
          endif
       endif
    enddo

    !
    ! rescale the fluxes:
    !

    ! moved assignment to temp_new from separate loop (see commented code below)
    ! into the following loops. Avoids having to save elements of flux array
    do iEdge = edgeStart,edgeEnd
       cell1 = cellsOnEdge(1,iEdge)
       cell2 = cellsOnEdge(2,iEdge)
       if(cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then
!DIR$ IVDEP
          do k = 1,nVertLevels
             flux = flux_arr(k,iEdge)
             flux = max(0.0_RKIND,flux) * min(scale_arr(k,SCALE_OUT,cell1), scale_arr(k,SCALE_IN, cell2)) &
                  + min(0.0_RKIND,flux) * min(scale_arr(k,SCALE_IN, cell1), scale_arr(k,SCALE_OUT,cell2))
             flux_arr(k,iEdge) = flux
          enddo
       endif
    enddo

    !
    ! rescale the vertical flux
    !
!$OMP BARRIER
    do iCell = cellSolveStart,cellSolveEnd
!DIR$ IVDEP
       do k = 2,nVertLevels
          flux = wdtn(k,iCell)
          flux = max(0.0_RKIND,flux)*min(scale_arr(k-1,SCALE_OUT,iCell),scale_arr(k  ,SCALE_IN,iCell)) &
               + min(0.0_RKIND,flux)*min(scale_arr(k  ,SCALE_OUT,iCell),scale_arr(k-1,SCALE_IN,iCell))
          wdtn(k,iCell) = flux
       enddo
    enddo

    !
    ! do the scalar update now that we have the fluxes:
    !
    do iCell = cellSolveStart,cellSolveEnd
       do i = 1,nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i,iCell)
!DIR$ IVDEP
          do k = 1,nVertLevels
             temp_new(k,iCell) = temp_new(k,iCell) - edgesOnCell_sign(i,iCell)*flux_arr(k,iEdge) &
                               * invAreaCell(iCell)
          enddo
       enddo

       if(local_advance_density) then
!DIR$ IVDEP
          do k = 1,nVertLevels
             temp_new(k,iCell) = (temp_new(k,iCell) + (-rdnw(k)*(wdtn(k+1,iCell)-wdtn(k,iCell)))) &
                               / rho_zz_int(k,iCell)
          enddo
       else
!DIR$ IVDEP
          do k = 1,nVertLevels
             temp_new(k,iCell) = (temp_new(k,iCell) + (-rdnw(k)*(wdtn(k+1,iCell)-wdtn(k,iCell)))) &
                               /  rho_zz_new(k,iCell)
          enddo
       endif
    enddo

#ifdef DEBUG_TRANSPORT
    scmin = scalar_new(1,1)
    scmax = scalar_new(1,1)
    do iCell = 1,nCellsSolve
       do k =1,nVertLevels
          scmax = max(scmax,scalar_new(k,iCell))
          scmin = min(scmin,scalar_new(k,iCell))
          if(s_max(k,iCell) < scalar_new(k,iCell)) then
             call mpas_log_write(' over - k,iCell,s_min,s_max,scalar_new ',intArgs=(/k,iCell/), &
                       realArgs=(/s_min(k,iCell),s_max(k,iCell),scalar_new(k,iCell)/))
          endif
          if(s_min(k,iCell) > scalar_new(k,iCell)) then
             call mpas_log_write(' under - k,iCell,s_min,s_max,scalar_new ',intArgs=(/k,iCell/), &
                       realArgs=(/s_min(k,iCell),s_max(k,iCell),scalar_new(k,iCell)/))
          endif
       enddo
    enddo
    call mpas_log_write(' scmin, scmax new out $r $r',realArgs=(/scmin,scmax/))
    call mpas_log_write(' icell_min, k_min ',intArgs=(/icellmax, kmax/))
#endif

    ! the update should be positive definite. but roundoff can sometimes leave small negative values
    ! hence the enforcement of PD in the copy back to the model state.
!$OMP BARRIER
    do iCell = cellStart,cellEnd
       if(bdyMaskCell(iCell) <= nSpecZone) then ! regional_MPAS does spec zone update after transport.
          do k =1,nVertLevels
             scalar_new(iScalar,k,iCell) = max(0.0_RKIND,temp_new(k,iCell))
          enddo
       endif
    enddo

 enddo

!call mpas_log_write('--- end subroutine atm_advance_chemistry_mono:")
 end subroutine atm_advance_chemistry_mono

!=================================================================================================================
 end module mpas_atmphys_advection
!=================================================================================================================

