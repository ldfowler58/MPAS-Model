! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================================
 module mpas_chemistry_gocart2G_update
 use mpas_log
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines,only: mpas_pool_get_dimension,mpas_pool_get_array

 use mpas_chemistry_gocart2G_date_time


 implicit none
 private
 public:: update_gocart2G_climatology, &
          vinterp_backgrounds


 contains


!=================================================================================================================
 subroutine update_gocart2G_climatology(current_date,mesh,gocart2G_met,gocart2G_backgrounds)
!=================================================================================================================

!--- input arguments:
 type(mpas_pool_type),intent(in):: mesh
 character(len=StrKIND),intent(in):: current_date

!--- inout arguments:
 type(mpas_pool_type),intent(inout):: gocart2G_met
 type(mpas_pool_type),intent(inout):: gocart2G_backgrounds

!--- local variables and pointers:


!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine update_gocart2G_climatology: '//trim(current_date))


!--- updates the climatological monthly-averaged GOCART-2G tropopause pressure to the current forecast time:
 call update_tinterp_background_ptrop(current_date,mesh,gocart2G_met) 


!--- updates the climatological monthly-averaged GOCART-2G DMS to the current forecast time:
 call update_tinterp_background_dms(current_date,mesh,gocart2G_backgrounds)


!--- updates the climatological monthly-averaged GOCART-2G H2O2, OH, and NO3 to the current forecast time:
 call update_tinterp_backgrounds(current_date,mesh,gocart2G_backgrounds)


 call mpas_log_write('--- end subroutine update_gocart2G_climatology.')

 end subroutine update_gocart2G_climatology

!=================================================================================================================
 subroutine update_tinterp_background_ptrop(current_date,mesh,gocart2G_met)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 character(len=StrKIND),intent(in):: current_date

!inout arguments:
 type(mpas_pool_type),intent(inout):: gocart2G_met

!local variables and pointers:
 integer,pointer:: nCells,nMonths
 real(kind=RKIND),dimension(:),pointer:: ptrop
 real(kind=RKIND),dimension(:,:),pointer:: ptrop_clim

 real(kind=RKIND),dimension(:),allocatable::   dummy2
 real(kind=RKIND),dimension(:,:),allocatable:: dummy1

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write('--- enter subroutine update_tinterp_background_ptrop: '//trim(current_date))

 call mpas_pool_get_dimension(mesh,'nCells' ,nCells )
 call mpas_pool_get_dimension(mesh,'nMonths',nMonths)

 call mpas_pool_get_array(mesh,'ptrop_gocart2G_clim',ptrop_clim)
 call mpas_pool_get_array(gocart2G_met,'background_ptrop',ptrop)


!--- interpolation of the monthly-averaged pressure at tropopause to the current date:
 if(.not.allocated(dummy2)) allocate(dummy2(nCells))
 if(.not.allocated(dummy1)) allocate(dummy1(nMonths,nCells))

 dummy2(1:nCells) = 0._RKIND
 dummy1(1:nMonths,1:nCells) = ptrop_clim(1:nMonths,1:nCells)
 call monthly_interp_to_date(nCells,current_date,dummy1,dummy2)
 ptrop(1:nCells) = dummy2(1:nCells)

 if(allocated(dummy2)) deallocate(dummy2)
 if(allocated(dummy1)) deallocate(dummy1)


!call mpas_log_write('--- end subroutine update_tinterp_background_ptrop.')

 end subroutine update_tinterp_background_ptrop

!=================================================================================================================
 subroutine update_tinterp_background_dms(current_date,mesh,gocart2G_backgrounds)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 character(len=StrKIND),intent(in):: current_date

!inout arguments:
 type(mpas_pool_type),intent(inout):: gocart2G_backgrounds

!local variables and pointers:
 integer,pointer:: nCells,nMonths
 real(kind=RKIND),dimension(:),pointer:: dms
 real(kind=RKIND),dimension(:,:),pointer:: dms_clim

 real(kind=RKIND),dimension(:),allocatable::   dummy2
 real(kind=RKIND),dimension(:,:),allocatable:: dummy1

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write('--- enter subroutine update_tinterp_background_dms: '//trim(current_date))

 call mpas_pool_get_dimension(mesh,'nCells' ,nCells )
 call mpas_pool_get_dimension(mesh,'nMonths',nMonths)

 call mpas_pool_get_array(mesh,'dms_gocart2G_clim',dms_clim)
 call mpas_pool_get_array(gocart2G_backgrounds,'background_dms',dms)


!--- interpolation of the monthly-averaged DMS to the current date:
 if(.not.allocated(dummy2)) allocate(dummy2(nCells))
 if(.not.allocated(dummy1)) allocate(dummy1(nMonths,nCells))

 dummy2(1:nCells) = 0._RKIND
 dummy1(1:nMonths,1:nCells) = dms_clim(1:nMonths,1:nCells)
 call monthly_interp_to_date(nCells,current_date,dummy1,dummy2)
 dms(1:nCells) = dummy2(1:nCells)

 if(allocated(dummy2)) deallocate(dummy2)
 if(allocated(dummy1)) deallocate(dummy1)


!call mpas_log_write('--- end subroutine update_tinterp_background_dms.')

 end subroutine update_tinterp_background_dms

!=================================================================================================================
 subroutine update_tinterp_backgrounds(current_date,mesh,gocart2G_backgrounds)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 character(len=StrKIND),intent(in):: current_date

!inout arguments:
 type(mpas_pool_type),intent(inout):: gocart2G_backgrounds

!local variables and pointers:
 integer,pointer:: nCells,nBCKLevels,nMonths
 integer:: k

 real(kind=RKIND),dimension(:,:),pointer:: oh,h2o2,no3
 real(kind=RKIND),dimension(:,:),pointer:: pres,dpres

 real(kind=RKIND),dimension(:,:,:),pointer:: oh_clim,h2o2_clim,no3_clim
 real(kind=RKIND),dimension(:,:,:),pointer:: pres_clim,dpres_clim

 real(kind=RKIND),dimension(:),allocatable::   dummy2
 real(kind=RKIND),dimension(:,:),allocatable:: dummy1

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write('--- enter subroutine update_tinterp_backgrounds: '//trim(current_date))

 call mpas_pool_get_dimension(mesh,'nCells'     ,nCells     )
 call mpas_pool_get_dimension(mesh,'nBCKLevels' ,nBCKLevels )
 call mpas_pool_get_dimension(mesh,'nMonths'    ,nMonths    )

 call mpas_pool_get_array(mesh,'oh_gocart2G_clim'   ,oh_clim   )
 call mpas_pool_get_array(mesh,'h2o2_gocart2G_clim' ,h2o2_clim )
 call mpas_pool_get_array(mesh,'no3_gocart2G_clim'  ,no3_clim  )
 call mpas_pool_get_array(mesh,'pres_gocart2G_clim' ,pres_clim )
 call mpas_pool_get_array(mesh,'dpres_gocart2G_clim',dpres_clim)

 call mpas_pool_get_array(gocart2G_backgrounds,'oh_gocart2G'   ,oh   )
 call mpas_pool_get_array(gocart2G_backgrounds,'h2o2_gocart2G' ,h2o2 )
 call mpas_pool_get_array(gocart2G_backgrounds,'no3_gocart2G'  ,no3  )
 call mpas_pool_get_array(gocart2G_backgrounds,'pres_gocart2G' ,pres )
 call mpas_pool_get_array(gocart2G_backgrounds,'dpres_gocart2G',dpres)


!--- interpolation of the monthly-averaged H2O2, OH, and NO3 volume mixing ratios to the current date:
 if(.not.allocated(dummy2)) allocate(dummy2(nCells))
 if(.not.allocated(dummy1)) allocate(dummy1(nMonths,nCells))

 do k = 1, nBCKLevels
    dummy2(1:nCells) = 0._RKIND
    dummy1(1:nMonths,1:nCells) = pres_clim(1:nMonths,k,1:nCells)
    call monthly_interp_to_date(nCells,current_date,dummy1,dummy2)
    pres(k,1:nCells) = dummy2(1:nCells)

    dummy2(1:nCells) = 0._RKIND
    dummy1(1:nMonths,1:nCells) = dpres_clim(1:nMonths,k,1:nCells)
    call monthly_interp_to_date(nCells,current_date,dummy1,dummy2)
    dpres(k,1:nCells) = dummy2(1:nCells)

    dummy2(1:nCells) = 0._RKIND
    dummy1(1:nMonths,1:nCells) = h2o2_clim(1:nMonths,k,1:nCells)
    call monthly_interp_to_date(nCells,current_date,dummy1,dummy2)
    h2o2(k,1:nCells) = dummy2(1:nCells)

    dummy2(1:nCells) = 0._RKIND
    dummy1(1:nMonths,1:nCells) = oh_clim(1:nMonths,k,1:nCells)
    call monthly_interp_to_date(nCells,current_date,dummy1,dummy2)
    oh(k,1:nCells) = dummy2(1:nCells)

    dummy2(1:nCells) = 0._RKIND
    dummy1(1:nMonths,1:nCells) = no3_clim(1:nMonths,k,1:nCells)
    call monthly_interp_to_date(nCells,current_date,dummy1,dummy2)
    no3(k,1:nCells) = dummy2(1:nCells)
 enddo

 if(allocated(dummy2)) deallocate(dummy2)
 if(allocated(dummy1)) deallocate(dummy1)


!call mpas_log_write('--- end subroutine update_tinterp_backgrounds.')

 end subroutine update_tinterp_backgrounds

!=================================================================================================================
 subroutine vinterp_backgrounds(mesh,diag,gocart2G_backgrounds)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: diag

!inout arguments:
 type(mpas_pool_type),intent(inout):: gocart2G_backgrounds

!local variables and pointers:
 integer,pointer:: nCells,nBCKLevels,nVertLevels,nMonths
 integer:: iCell,k,kk,n

 real(kind=RKIND),dimension(:,:),pointer:: pressure
 real(kind=RKIND),dimension(:,:),pointer:: oh,h2o2,no3
 real(kind=RKIND),dimension(:,:),pointer:: oh_clim,h2o2_clim,no3_clim
 real(kind=RKIND),dimension(:,:),pointer:: pres_clim,dpres_clim

 real(kind=RKIND):: target_p
 real(kind=RKIND),dimension(:,:),allocatable:: sorted_arr

!------------------------------------------------------------------------------------------------------------------
!call mpas_log_write('--- enter subroutine init_vinterp_backgrounds:')

 call mpas_pool_get_dimension(mesh,'nCells'     ,nCells     )
 call mpas_pool_get_dimension(mesh,'nBCKLevels' ,nBCKLevels )
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 call mpas_pool_get_array(diag,'pressure_base',pressure)

 call mpas_pool_get_array(gocart2G_backgrounds,'oh_gocart2G'   ,oh_clim   )
 call mpas_pool_get_array(gocart2G_backgrounds,'h2o2_gocart2G' ,h2o2_clim )
 call mpas_pool_get_array(gocart2G_backgrounds,'no3_gocart2G'  ,no3_clim  )
 call mpas_pool_get_array(gocart2G_backgrounds,'pres_gocart2G' ,pres_clim )
 call mpas_pool_get_array(gocart2G_backgrounds,'dpres_gocart2G',dpres_clim)

 call mpas_pool_get_array(gocart2G_backgrounds,'background_oh'  ,oh  )
 call mpas_pool_get_array(gocart2G_backgrounds,'background_h2o2',h2o2)
 call mpas_pool_get_array(gocart2G_backgrounds,'background_no3' ,no3 )


!--- allocation of local arrays:
 if(.not.allocated(sorted_arr)) allocate(sorted_arr(2,nBCKLevels))


!--- H2O2:
 do iCell = 1, nCells
    sorted_arr(1,1:nBCKLevels) = 0._RKIND
    sorted_arr(2,1:nBCKLevels) = 0._RKIND
    do k = 1, nBCKLevels
       kk = nBCKLevels + 1 -k
       sorted_arr(1,kk) = pres_clim(k,iCell)
       sorted_arr(2,kk) = h2o2_clim(k,iCell)
    enddo
    do k = nVertLevels, 1, -1
       target_p = pressure(k,iCell)
       h2o2(k,iCell) = vertical_interp(target_p,nBCKLevels-1, &
                          sorted_arr(:,1:nBCKLevels-1),order=1,extrap=0)
       if(target_p >= pres_clim(1,iCell) .and. k < nVertLevels) h2o2(k,iCell) = h2o2(k+1,iCell)
    enddo
 enddo


!--- OH:
 do iCell = 1, nCells
    sorted_arr(1,1:nBCKLevels) = 0._RKIND
    sorted_arr(2,1:nBCKLevels) = 0._RKIND
    do k = 1, nBCKLevels
       kk = nBCKLevels + 1 -k
       sorted_arr(1,kk) = pres_clim(k,iCell)
       sorted_arr(2,kk) = oh_clim(k,iCell)
    enddo
    do k = nVertLevels, 1, -1
       target_p = pressure(k,iCell)
       oh(k,iCell) = vertical_interp(target_p,nBCKLevels-1, &
                        sorted_arr(:,1:nBCKLevels-1),order=1,extrap=0)
       if(target_p >= pres_clim(1,iCell) .and. k < nVertLevels) oh(k,iCell) = oh(k+1,iCell)
    enddo
 enddo


!--- NO3:
 do iCell = 1, nCells
    sorted_arr(1,1:nBCKLevels) = 0._RKIND
    sorted_arr(2,1:nBCKLevels) = 0._RKIND
    do k = 1, nBCKLevels
       kk = nBCKLevels + 1 -k
       sorted_arr(1,kk) = pres_clim(k,iCell)
       sorted_arr(2,kk) = no3_clim(k,iCell)
    enddo
    do k = nVertLevels, 1, -1
       target_p = pressure(k,iCell)
       no3(k,iCell) = vertical_interp(target_p,nBCKLevels-1, &
                         sorted_arr(:,1:nBCKLevels-1),order=1,extrap=0)
       if(target_p >= pres_clim(1,iCell) .and. k < nVertLevels) no3(k,iCell) = no3(k+1,iCell)
    enddo
 enddo


!--- deallocation of local arrays:
 if(allocated(sorted_arr)) deallocate(sorted_arr)


!call mpas_log_write('--- end subroutine init_vinterp_backgrounds.')

 end subroutine vinterp_backgrounds

!=================================================================================================================
 real(kind=RKIND) function vertical_interp(target_z,nz,zf,order,extrap,surface_val,sealev_val)
!=================================================================================================================

 implicit none

 real(kind=RKIND),intent(in):: target_z
 integer,intent(in):: nz
 real(kind=RKIND),dimension(2,nz),intent(in):: zf
 integer,intent(in),optional:: order
 integer,intent(in),optional:: extrap
 real(kind=RKIND),intent(in),optional:: surface_val
 real(kind=RKIND),intent(in),optional:: sealev_val

 integer:: k, lm, lp
 real(kind=RKIND):: wm, wp
 real(kind=RKIND):: slope

 integer:: interp_order, extrap_type
 real(kind=RKIND):: surface, sealevel

!------------------------------------------------------------------------------------------------------------------

 if(present(order)) then
    interp_order = order
 else
    interp_order = 2
 end if

 if(present(extrap)) then
    extrap_type = extrap
 else
    extrap_type = 1
 end if

 if(present(surface_val)) then
    surface = surface_val
 else
    surface = 200100.0
 end if

 if(present(sealev_val)) then
    sealevel = sealev_val
 else
    sealevel = 201300.0
 end if

!
! Extrapolation required
!
 if(target_z < zf(1,1)) then
    if(extrap_type == 0) then
       vertical_interp = zf(2,1)
    else if(extrap_type == 1) then
       slope = (zf(2,2) - zf(2,1)) / (zf(1,2) - zf(1,1))
       vertical_interp = zf(2,1) + slope * (target_z - zf(1,1))
    end if
    return
 end if
 if(target_z >= zf(1,nz)) then
    if(extrap_type == 0) then
       vertical_interp = zf(2,nz)
    else if(extrap_type == 1) then
       slope = (zf(2,nz) - zf(2,nz-1)) / (zf(1,nz) - zf(1,nz-1))
       vertical_interp = zf(2,nz) + slope * (target_z - zf(1,nz))
    end if
    return
 end if


!
! No extrapolation required
!
 do k=1,nz-1
    if(target_z >= zf(1,k) .and. target_z < zf(1,k+1)) then
       lm = k
       lp = k+1
       wm = (zf(1,k+1) - target_z) / (zf(1,k+1) - zf(1,k))
       wp = (target_z - zf(1,k)) / (zf(1,k+1) - zf(1,k))
       exit
    end if
 end do

 vertical_interp = wm*zf(2,lm) + wp*zf(2,lp)

 return

 end function vertical_interp

!==================================================================================================================
 end module mpas_chemistry_gocart2G_update
!==================================================================================================================
