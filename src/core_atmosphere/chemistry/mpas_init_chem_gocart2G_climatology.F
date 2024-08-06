! Copyright (c) 2024 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_init_chem_gocart2G_climatology
 use mpas_derived_types
 use mpas_dmpar,only: mpas_dmpar_abort
 use mpas_pool_routines
 use mpas_kind_types

 use init_atm_read_met
 use init_atm_hinterp
 use init_atm_llxy
 use init_atm_vinterp
 use mpas_atmphys_date_time

 implicit none
 private
 public:: init_chem_gocart2G_climatology 


 contains


!=================================================================================================================
 subroutine init_chem_gocart2G_climatology(block,configs,mesh,diag,gocart2G_met,gocart2G_backgrounds)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: diag

!inout arguments:
 type(block_type),intent(inout),target:: block
 type(mpas_pool_type),intent(inout):: gocart2G_met
 type(mpas_pool_type),intent(inout):: gocart2G_backgrounds

!local variables and pointers:
 character (len=StrKIND),pointer:: config_start_time
 character(len=StrKIND):: initial_date
 character(len=StrKIND):: filename_dms,filename_sulf,filename_ptrop

 logical:: lexist

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine init_chem_gocart2G_climatology:')

!inquire if the input file of the monthly-mean climatological tropopause pressure exists:
 lexist = .false.
 filename_ptrop = "BACKGROUND_ptrop.dat"
 
 inquire(file=filename_ptrop,exist=lexist)
 if(lexist) then

    call mpas_pool_get_config(configs,'config_start_time',config_start_time)

    !--- horizontal interpolation of the climatological monthly-averaged GOCART-2G tropopause pressure
    !    data to the MPAS mesh:
    call init_hinterp_background_ptrop(filename_ptrop,block,mesh,gocart2G_met)

    !--- interpolation of the climatological monthly-averaged GOCART-2G tropopause pressure to the
    !    initial forecast time:
    initial_date = trim(config_start_time)
    call init_tinterp_background_ptrop(initial_date,mesh,gocart2G_met)
 else
    call mpas_log_write(trim(filename_ptrop)//' was not found in local directory: ',messageType=MPAS_LOG_CRIT)
 endif


!inquire if the input file of the monthly-mean climatological background DMS exists:
 lexist = .false.
 filename_dms = "BACKGROUND_dms.dat"

 inquire(file=filename_dms,exist=lexist)
 if(lexist) then

    call mpas_pool_get_config(configs,'config_start_time',config_start_time)

    !--- horizontal interpolation of the climatological monthly-averaged GOCART-2G DMS
    !    data to the MPAS mesh:
    call init_hinterp_background_dms(filename_dms,block,mesh,gocart2G_backgrounds)

    !--- interpolation of the climatological monthly-averaged GOCART-2G DMS to the initial forecast time:
    initial_date = trim(config_start_time)
    call init_tinterp_background_dms(initial_date,mesh,gocart2G_backgrounds)
 else
    call mpas_log_write(trim(filename_dms)//' was not found in local directory:',messageType=MPAS_LOG_CRIT)
 endif


!inquire if the input file of the monthly-mean climatology background H2O2, NO3, AND OH exist:
 lexist = .false.
 filename_sulf = "BACKGROUND_sulf.dat"

 inquire(file=filename_sulf,exist=lexist)
 if(lexist) then

    call mpas_pool_get_config(configs,'config_start_time',config_start_time)

    !--- horizontal interpolation of the climatological monthly-averaged GOCART-2G H2O2, OH, and NO3
    !    data to the MPAS mesh:
    call init_hinterp_background_sulf(filename_sulf,block,mesh,gocart2G_backgrounds)

    !--- horizontal interpolation of the climatological monthly-averaged GOCART-2G H2O2, OH, and NO3
    !    to the initial forecast time:
    initial_date = trim(config_start_time)
    call init_vinterp_backgrounds(initial_date,mesh,diag,gocart2G_backgrounds)
 else
    call mpas_log_write(trim(filename_sulf)//' was not found in local directory:',messageType=MPAS_LOG_CRIT)
 endif

 call mpas_log_write('--- end subroutine init_chem_gocart2G_climatology.')

 end subroutine init_chem_gocart2G_climatology

!=================================================================================================================
 subroutine init_tinterp_background_ptrop(initial_date,mesh,gocart2G_met)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 character(len=StrKIND),intent(in):: initial_date

!inout arguments:
 type(mpas_pool_type),intent(inout):: gocart2G_met

!local variables and pointers:
 integer,pointer:: nCells,nMonths
 real(kind=RKIND),dimension(:),pointer:: ptrop
 real(kind=RKIND),dimension(:,:),pointer:: ptrop_clim

 real(kind=RKIND),dimension(:),allocatable::   dummy2
 real(kind=RKIND),dimension(:,:),allocatable:: dummy1

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine init_tinterp_background_ptrop:')

 call mpas_pool_get_dimension(mesh,'nCells' ,nCells )
 call mpas_pool_get_dimension(mesh,'nMonths',nMonths)

 call mpas_pool_get_array(gocart2G_met,'ptrop_gocart2G_clim',ptrop_clim)
 call mpas_pool_get_array(gocart2G_met,'background_ptrop',ptrop)


 if(.not.allocated(dummy2)) allocate(dummy2(nCells))
 if(.not.allocated(dummy1)) allocate(dummy1(nMonths,nCells))

 dummy2(1:nCells) = 0._RKIND
 dummy1(1:nMonths,1:nCells) = ptrop_clim(1:nMonths,1:nCells)
 call monthly_interp_to_date(nCells,initial_date,dummy1,dummy2)
 ptrop(1:nCells) = dummy2(1:nCells)

 if(allocated(dummy2)) deallocate(dummy2)
 if(allocated(dummy1)) deallocate(dummy1)


 call mpas_log_write('--- end subroutine init_tinterp_background_ptrop.')

 end subroutine init_tinterp_background_ptrop

!=================================================================================================================
 subroutine init_tinterp_background_dms(initial_date,mesh,gocart2G_backgrounds)
!=================================================================================================================
  
!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 character(len=StrKIND),intent(in):: initial_date

!inout arguments:
 type(mpas_pool_type),intent(inout):: gocart2G_backgrounds

!local variables and pointers:
 integer,pointer:: nCells,nMonths
 real(kind=RKIND),dimension(:),pointer:: dms
 real(kind=RKIND),dimension(:,:),pointer:: dms_clim

 real(kind=RKIND),dimension(:),allocatable::   dummy2
 real(kind=RKIND),dimension(:,:),allocatable:: dummy1

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine init_tinterp_background_dms:')

 call mpas_pool_get_dimension(mesh,'nCells' ,nCells )
 call mpas_pool_get_dimension(mesh,'nMonths',nMonths)

 call mpas_pool_get_array(gocart2G_backgrounds,'dms_gocart2G_clim',dms_clim)
 call mpas_pool_get_array(gocart2G_backgrounds,'background_dms',dms)


 if(.not.allocated(dummy2)) allocate(dummy2(nCells))
 if(.not.allocated(dummy1)) allocate(dummy1(nMonths,nCells))

 dummy2(1:nCells) = 0._RKIND
 dummy1(1:nMonths,1:nCells) = dms_clim(1:nMonths,1:nCells)
 call monthly_interp_to_date(nCells,initial_date,dummy1,dummy2)
 dms(1:nCells) = dummy2(1:nCells)

 if(allocated(dummy2)) deallocate(dummy2)
 if(allocated(dummy1)) deallocate(dummy1)


 call mpas_log_write('--- end subroutine init_tinterp_background_dms.')

 end subroutine init_tinterp_background_dms

!=================================================================================================================
 subroutine init_vinterp_backgrounds(initial_date,mesh,diag,gocart2G_backgrounds)
!=================================================================================================================

!input arguments:
 character(len=StrKIND),intent(in):: initial_date
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: diag

!inout arguments:
 type(mpas_pool_type),intent(inout):: gocart2G_backgrounds

!local variables and pointers:
 integer,pointer:: nCells,nBCKLevels,nVertLevels,nMonths
 integer:: iCell,k,kk,n

 real(kind=RKIND),dimension(:,:),pointer:: pressure
 real(kind=RKIND),dimension(:,:),pointer:: oh,h2o2,no3
 real(kind=RKIND),dimension(:,:,:),pointer:: oh_clim,h2o2_clim,no3_clim
 real(kind=RKIND),dimension(:,:,:),pointer:: pres_clim,dpres_clim


 real(kind=RKIND):: target_p
 real(kind=RKIND),dimension(:,:),allocatable:: oh_int,h2o2_int,no3_int,pres_int,sorted_arr

 real(kind=RKIND),dimension(:),allocatable::   dummy2
 real(kind=RKIND),dimension(:,:),allocatable:: dummy1

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine init_vinterp_backgrounds:')

 call mpas_pool_get_dimension(mesh,'nCells'     ,nCells     )
 call mpas_pool_get_dimension(mesh,'nBCKLevels' ,nBCKLevels )
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
 call mpas_pool_get_dimension(mesh,'nMonths'    ,nMonths    )

 call mpas_pool_get_array(diag,'pressure_base',pressure)

 call mpas_pool_get_array(gocart2G_backgrounds,'oh_gocart2G_clim'   ,oh_clim   )
 call mpas_pool_get_array(gocart2G_backgrounds,'h2o2_gocart2G_clim' ,h2o2_clim )
 call mpas_pool_get_array(gocart2G_backgrounds,'no3_gocart2G_clim'  ,no3_clim  )
 call mpas_pool_get_array(gocart2G_backgrounds,'pres_gocart2G_clim' ,pres_clim )
 call mpas_pool_get_array(gocart2G_backgrounds,'dpres_gocart2G_clim',dpres_clim)

 call mpas_pool_get_array(gocart2G_backgrounds,'background_oh'  ,oh  )
 call mpas_pool_get_array(gocart2G_backgrounds,'background_h2o2',h2o2)
 call mpas_pool_get_array(gocart2G_backgrounds,'background_no3' ,no3 )


!--- allocation of local arrays:
 if(.not.allocated(oh_int)  ) allocate(oh_int(nBCKLevels,nCells)  )
 if(.not.allocated(h2o2_int)) allocate(h2o2_int(nBCKLevels,nCells))
 if(.not.allocated(no3_int) ) allocate(no3_int(nBCKLevels,nCells) )
 if(.not.allocated(pres_int)) allocate(pres_int(nBCKLevels,nCells))
 if(.not.allocated(sorted_arr)) allocate(sorted_arr(2,nBCKLevels))


!--- interpolation of the monthly-averaged GOCART data to the initial date, and vertical interpolation to the
!    MPAS levels:
 if(.not.allocated(dummy2)) allocate(dummy2(nCells))
 if(.not.allocated(dummy1)) allocate(dummy1(nMonths,nCells))

 do k = 1, nBCKLevels
    dummy2(1:nCells) = 0._RKIND
    dummy1(1:nMonths,1:nCells) = pres_clim(1:nMonths,k,1:nCells)
    call monthly_interp_to_date(nCells,initial_date,dummy1,dummy2)
    pres_int(k,1:nCells) = dummy2(1:nCells)
 enddo

!--- H2O2:
 do k = 1, nBCKLevels
    dummy2(1:nCells) = 0._RKIND
    dummy1(1:nMonths,1:nCells) = h2o2_clim(1:nMonths,k,1:nCells)
    call monthly_interp_to_date(nCells,initial_date,dummy1,dummy2)
    h2o2_int(k,1:nCells) = dummy2(1:nCells)
 enddo
 do iCell = 1, nCells
    sorted_arr(1,1:nBCKLevels) = 0._RKIND
    sorted_arr(2,1:nBCKLevels) = 0._RKIND
    do k = 1, nBCKLevels
       kk = nBCKLevels + 1 -k
       sorted_arr(1,kk) = pres_int(k,iCell)
       sorted_arr(2,kk) = h2o2_int(k,iCell)
    enddo
    do k = nVertLevels, 1, -1
       target_p = pressure(k,iCell)
       h2o2(k,iCell) = vertical_interp(target_p,nBCKLevels-1, &
                          sorted_arr(:,1:nBCKLevels-1),order=1,extrap=0)
       if(target_p >= pres_int(1,iCell) .and. k < nVertLevels) h2o2(k,iCell) = h2o2(k+1,iCell)
    enddo
 enddo


!--- OH:
 do k = 1, nBCKLevels
    dummy2(1:nCells) = 0._RKIND
    dummy1(1:nMonths,1:nCells) = oh_clim(1:nMonths,k,1:nCells)
    call monthly_interp_to_date(nCells,initial_date,dummy1,dummy2)
    oh_int(k,1:nCells) = dummy2(1:nCells)
 enddo
 do iCell = 1, nCells
    sorted_arr(1,1:nBCKLevels) = 0._RKIND
    sorted_arr(2,1:nBCKLevels) = 0._RKIND
    do k = 1, nBCKLevels
       kk = nBCKLevels + 1 -k
       sorted_arr(1,kk) = pres_int(k,iCell)
       sorted_arr(2,kk) = oh_int(k,iCell)
    enddo
    do k = nVertLevels, 1, -1
       target_p = pressure(k,iCell)
       oh(k,iCell) = vertical_interp(target_p,nBCKLevels-1, &
                        sorted_arr(:,1:nBCKLevels-1),order=1,extrap=0)
       if(target_p >= pres_int(1,iCell) .and. k < nVertLevels) oh(k,iCell) = oh(k+1,iCell)
    enddo
 enddo


!--- NO3:
 do k = 1, nBCKLevels
    dummy2(1:nCells) = 0._RKIND
    dummy1(1:nMonths,1:nCells) = no3_clim(1:nMonths,k,1:nCells)
    call monthly_interp_to_date(nCells,initial_date,dummy1,dummy2)
    no3_int(k,1:nCells) = dummy2(1:nCells)
 enddo
 do iCell = 1, nCells
    sorted_arr(1,1:nBCKLevels) = 0._RKIND
    sorted_arr(2,1:nBCKLevels) = 0._RKIND
    do k = 1, nBCKLevels
       kk = nBCKLevels + 1 -k
       sorted_arr(1,kk) = pres_int(k,iCell)
       sorted_arr(2,kk) = no3_int(k,iCell)
    enddo
    do k = nVertLevels, 1, -1
       target_p = pressure(k,iCell)
       no3(k,iCell) = vertical_interp(target_p,nBCKLevels-1, &
                         sorted_arr(:,1:nBCKLevels-1),order=1,extrap=0)
       if(target_p >= pres_int(1,iCell) .and. k < nVertLevels) no3(k,iCell) = no3(k+1,iCell)
    enddo
 enddo


!--- deallocation of local arrays:
 if(allocated(dummy1)  ) deallocate(dummy1  )
 if(allocated(dummy2)  ) deallocate(dummy2  )
 if(allocated(oh_int)  ) deallocate(oh_int  )
 if(allocated(h2o2_int)) deallocate(h2o2_int)
 if(allocated(no3_int) ) deallocate(no3_int )
 if(allocated(pres_int)) deallocate(pres_int)
 if(allocated(sorted_arr)) deallocate(sorted_arr)


 call mpas_log_write('--- end subroutine init_vinterp_backgrounds.')

 end subroutine init_vinterp_backgrounds

!=================================================================================================================
 subroutine init_hinterp_background_ptrop(filename_ptrop,block,mesh,gocart2G_met)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 character(len=StrKIND),intent(in):: filename_ptrop

!inout arguments:
 type(block_type),intent(inout),target:: block
 type(mpas_pool_type),intent(inout):: gocart2G_met

!local variables:
 type(dm_info),pointer:: dminfo
 type(met_data) :: field !real*4 meteorological data.
 type(proj_info):: proj

 logical:: have_landmask

 integer:: iCell,istatus,k,masked,nmonths,nInterpPoints
 integer,pointer:: nCells
 integer,dimension(5):: interp_list
 integer,dimension(:),pointer:: landmask
 integer,dimension(:),pointer:: mask_array

 real(kind=RKIND):: fillval,maskval,msgval
 real(kind=RKIND):: lat,lon,x,y
 real(kind=RKIND),dimension(:),pointer:: latCell,lonCell
 real(kind=RKIND),dimension(:),pointer:: latPoints,lonPoints
 real(kind=RKIND),dimension(:,:),pointer:: destField2d

 real(kind=RKIND),dimension(:,:),allocatable:: maskslab,rslab

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine init_hinterp_ptrop:')

 dminfo => block%domain%dminfo

 call mpas_pool_get_dimension(mesh,'nCells',nCells)

 call mpas_pool_get_array(mesh,'landmask',landmask)
 call mpas_pool_get_array(mesh,'latCell' ,latCell )
 call mpas_pool_get_array(mesh,'lonCell' ,lonCell )


!open intermediate file:
 istatus = 0
 call read_met_init(trim(filename_ptrop),.true.,'not needed',istatus)
 if(istatus /= 0) then
    call mpas_log_write('********************************************************************************')
    call mpas_log_write('Error opening gocart2G file '//trim(filename_ptrop))
    call mpas_log_write('********************************************************************************')
    call mpas_dmpar_abort(dminfo)
 else
    call mpas_log_write('Processing file '//trim(filename_ptrop))
 end if


!scan through all fields in the file, looking for the LANDSEA field:
 have_landmask = .false.
 call read_next_met_field(field,istatus)
 do while (istatus == 0)
    if(index(field % field, 'LANDSEA') /= 0) then
       have_landmask = .true.
       if(.not.allocated(maskslab)) allocate(maskslab(-2:field%nx+3,field%ny))

       maskslab(1:field%nx,1:field%ny) = field%slab(1:field%nx,1:field%ny)
       maskslab(0 ,1:field%ny) = field%slab(field%nx  ,1:field%ny)
       maskslab(-1,1:field%ny) = field%slab(field%nx-1,1:field%ny)
       maskslab(-2,1:field%ny) = field%slab(field%nx-2,1:field%ny)
       maskslab(field%nx+1,1:field%ny) = field%slab(1,1:field%ny)
       maskslab(field%nx+2,1:field%ny) = field%slab(2,1:field%ny)
       maskslab(field%nx+3,1:field%ny) = field%slab(3,1:field%ny)
!      call mpas_log_write('minval, maxval of LANDSEA = $r $r',realArgs=(/minval(maskslab),maxval(maskslab)/))
    end if
    deallocate(field%slab)
    call read_next_met_field(field,istatus)
 end do
 call read_met_close()

 if(.not. have_landmask) then
    call mpas_log_write('********************************************************************************')
    call mpas_log_write('Landsea mask not available from the surface file')
    call mpas_log_write('********************************************************************************')
 end if


!read climatological tropopause pressure data:
 istatus = 0
 call read_met_init(trim(filename_ptrop),.true.,'not needed',istatus)
 if(istatus /= 0) then
    call mpas_log_write('********************************************************************************')
    call mpas_log_write('Error opening file '// trim(filename_ptrop))
    call mpas_log_write('********************************************************************************')
    call mpas_dmpar_abort(dminfo)
 endif
 call read_next_met_field(field,istatus)


!horizontally interpolate climatological tropopause pressure  data:
 do while(istatus == 0)

    interp_list(1) = FOUR_POINT
    interp_list(2) = W_AVERAGE4
    interp_list(3) = W_AVERAGE16
    interp_list(4) = SEARCH
    interp_list(5) = 0

    maskval = -1.0
    masked  = -1
    fillval = 0.0
    msgval  = 0.0

    mask_array => landmask

    if(index(field % field, 'PTROP_JAN') /= 0 .or. &
       index(field % field, 'PTROP_FEB') /= 0 .or. &
       index(field % field, 'PTROP_MAR') /= 0 .or. &
       index(field % field, 'PTROP_APR') /= 0 .or. &
       index(field % field, 'PTROP_MAY') /= 0 .or. &
       index(field % field, 'PTROP_JUN') /= 0 .or. &
       index(field % field, 'PTROP_JUL') /= 0 .or. &
       index(field % field, 'PTROP_AUG') /= 0 .or. &
       index(field % field, 'PTROP_SEP') /= 0 .or. &
       index(field % field, 'PTROP_OCT') /= 0 .or. &
       index(field % field, 'PTROP_NOV') /= 0 .or. &
       index(field % field, 'PTROP_DEC') /= 0) then

       !
       !set up projection:
       !
       call map_init(proj)

       if(field%iproj == PROJ_LATLON) then
          call map_set(PROJ_LATLON,proj, &
                       latinc = real(field%deltalat,RKIND), &
                       loninc = real(field%deltalon,RKIND), &
                       knowni = 1.0_RKIND, &
                       knownj = 1.0_RKIND, &
                       lat1   = real(field%startlat,RKIND), &
                       lon1   = real(field%startlon,RKIND))
       elseif(field%iproj == PROJ_GAUSS) then
          call map_set(PROJ_GAUSS,proj, &
                       nlat = nint(field%deltalat), &
                       loninc = 360.0_RKIND / real(field%nx,RKIND), &
                       lat1 = real(field%startlat,RKIND), &
                       lon1 = real(field%startlon,RKIND))
       endif

       !
       !horizontally interpolate field at level k:
       !
       if(index(field%field,'PTROP_JAN') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating PTROP_JAN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_met,'ptrop_gocart2G_clim',destField2d)
          nmonths = 1
       elseif(index(field%field,'PTROP_FEB') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating PTROP_FEB at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_met,'ptrop_gocart2G_clim',destField2d)
          nmonths = 2
       elseif(index(field%field,'PTROP_MAR') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating PTROP_MAR at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_met,'ptrop_gocart2G_clim',destField2d)
          nmonths = 3
       elseif(index(field%field,'PTROP_APR') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating PTROP_APR at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_met,'ptrop_gocart2G_clim',destField2d)
          nmonths = 4
       elseif(index(field%field,'PTROP_MAY') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating PTROP_MAY at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_met,'ptrop_gocart2G_clim',destField2d)
          nmonths = 5
       elseif(index(field%field,'PTROP_JUN') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating PTROP_JUN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_met,'ptrop_gocart2G_clim',destField2d)
          nmonths = 6
       elseif(index(field%field,'PTROP_JUL') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating PTROP_JUL at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_met,'ptrop_gocart2G_clim',destField2d)
          nmonths = 7
       elseif(index(field%field,'PTROP_AUG') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating PTROP_AUG at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_met,'ptrop_gocart2G_clim',destField2d)
          nmonths = 8
       elseif(index(field%field,'PTROP_SEP') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating PTROP_SEP at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_met,'ptrop_gocart2G_clim',destField2d)
          nmonths = 9
       elseif(index(field%field,'PTROP_OCT') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating PTROP_OCT at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_met,'ptrop_gocart2G_clim',destField2d)
          nmonths = 10
       elseif(index(field%field,'PTROP_NOV') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating PTROP_NOV at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_met,'ptrop_gocart2G_clim',destField2d)
          nmonths = 11
       elseif(index(field%field,'PTROP_DEC') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating PTROP_DEC at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_met,'ptrop_gocart2G_clim',destField2d)
          nmonths = 12
       endif

       allocate(rslab(-2:field%nx+3,field%ny))
       rslab(1:field%nx,1:field%ny) = field%slab(1:field%nx,1:field%ny)
       rslab(0,1:field%ny)  = field%slab(field%nx  ,1:field%ny)
       rslab(-1,1:field%ny) = field%slab(field%nx-1,1:field%ny)
       rslab(-2,1:field%ny) = field%slab(field%nx-2,1:field%ny)
       rslab(field%nx+1,1:field%ny) = field%slab(1,1:field%ny)
       rslab(field%nx+2,1:field%ny) = field%slab(2,1:field%ny)
       rslab(field%nx+3,1:field%ny) = field%slab(3,1:field%ny)

       do iCell = 1, nInterpPoints
          if(mask_array(iCell) /= masked) then
             lat = latPoints(iCell)*DEG_PER_RAD
             lon = lonPoints(iCell)*DEG_PER_RAD
             call latlon_to_ij(proj,lat,lon,x,y)
             if(x < 0.5) then
                lon = lon + 360.0
                call latlon_to_ij(proj,lat,lon,x,y)
             elseif(x > real(field%nx,kind=RKIND)+ 0.5) then
                lon = lon - 360.0
                call latlon_to_ij(proj,lat,lon,x,y)
             endif

             if(maskval /= -1.0) then
                destField2d(nmonths,iCell) = interp_sequence(x,y,1,rslab,-2,field%nx+3,1,field%ny,1,1,msgval, &
                                                  interp_list,1,maskval=maskval,mask_array=maskslab)
             else
                destField2d(nmonths,iCell) = interp_sequence(x,y,1,rslab,-2,field%nx+3,1,field%ny,1,1,msgval, &
                                                  interp_list,1)
             endif
          else
             destField2d(nmonths,iCell) = fillval
          endif
       enddo
       deallocate(rslab)

    endif
    deallocate(field%slab)
    call read_next_met_field(field,istatus)

 enddo

 call read_met_close()

 call mpas_log_write('--- end subroutine init_hinterp_ptrop.')

 end subroutine init_hinterp_background_ptrop

!=================================================================================================================
 subroutine init_hinterp_background_dms(filename_dms,block,mesh,gocart2G_backgrounds)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 character(len=StrKIND),intent(in):: filename_dms

!inout arguments:
 type(block_type),intent(inout),target:: block
 type(mpas_pool_type),intent(inout):: gocart2G_backgrounds

!local variables:
 type(dm_info),pointer:: dminfo
 type(met_data) :: field !real*4 meteorological data.
 type(proj_info):: proj

 logical:: have_landmask

 integer:: iCell,istatus,k,masked,nmonths,nInterpPoints
 integer,pointer:: nCells
 integer,dimension(5):: interp_list
 integer,dimension(:),pointer:: landmask
 integer,dimension(:),pointer:: mask_array

 real(kind=RKIND):: fillval,maskval,msgval
 real(kind=RKIND):: lat,lon,x,y
 real(kind=RKIND),dimension(:),pointer:: latCell,lonCell
 real(kind=RKIND),dimension(:),pointer:: latPoints,lonPoints
 real(kind=RKIND),dimension(:,:),pointer:: destField2d

 real(kind=RKIND),dimension(:,:),allocatable:: maskslab,rslab

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine init_hinterp_dms:')

 dminfo => block%domain%dminfo

 call mpas_pool_get_dimension(mesh,'nCells',nCells)

 call mpas_pool_get_array(mesh,'landmask',landmask)
 call mpas_pool_get_array(mesh,'latCell' ,latCell )
 call mpas_pool_get_array(mesh,'lonCell' ,lonCell )


!open intermediate file:
 istatus = 0
 call read_met_init(trim(filename_dms),.true.,'not needed',istatus)
 if(istatus /= 0) then
    call mpas_log_write('********************************************************************************')
    call mpas_log_write('Error opening gocart2G2G file '//trim(filename_dms))
    call mpas_log_write('********************************************************************************')
    call mpas_dmpar_abort(dminfo)
 else
    call mpas_log_write('Processing file '//trim(filename_dms))
 end if


!scan through all fields in the file, looking for the LANDSEA field:
 have_landmask = .false.
 call read_next_met_field(field,istatus)
 do while (istatus == 0)
    if(index(field % field, 'LANDSEA') /= 0) then
       have_landmask = .true.
       if(.not.allocated(maskslab)) allocate(maskslab(-2:field%nx+3,field%ny))


       maskslab(1:field%nx,1:field%ny) = field%slab(1:field%nx,1:field%ny)
       maskslab(0 ,1:field%ny) = field%slab(field%nx  ,1:field%ny)
       maskslab(-1,1:field%ny) = field%slab(field%nx-1,1:field%ny)
       maskslab(-2,1:field%ny) = field%slab(field%nx-2,1:field%ny)
       maskslab(field%nx+1,1:field%ny) = field%slab(1,1:field%ny)
       maskslab(field%nx+2,1:field%ny) = field%slab(2,1:field%ny)
       maskslab(field%nx+3,1:field%ny) = field%slab(3,1:field%ny)
!      call mpas_log_write('minval, maxval of LANDSEA = $r $r',realArgs=(/minval(maskslab),maxval(maskslab)/))
    end if
    deallocate(field%slab)
    call read_next_met_field(field,istatus)
 end do
 call read_met_close()

 if(.not. have_landmask) then
    call mpas_log_write('********************************************************************************')
    call mpas_log_write('Landsea mask not available from the surface file')
    call mpas_log_write('********************************************************************************')
 end if


!read climatological DMS data:
 istatus = 0
 call read_met_init(trim(filename_dms),.true.,'not needed',istatus)
 if(istatus /= 0) then
    call mpas_log_write('********************************************************************************')
    call mpas_log_write('Error opening file '// trim(filename_dms))
    call mpas_log_write('********************************************************************************')
    call mpas_dmpar_abort(dminfo)
 endif
 call read_next_met_field(field,istatus)


!horizontally interpolate DMS data:
 do while(istatus == 0)

!   interp_list(1) = FOUR_POINT
!   interp_list(2) = W_AVERAGE4
!   interp_list(3) = W_AVERAGE16
!   interp_list(4) = SEARCH
!   interp_list(5) = 0

!   maskval = -1.0
!   masked  = -1
!   fillval = 0.0
!   msgval  = 0.0

    interp_list(1) = FOUR_POINT
    interp_list(2) = W_AVERAGE4
    interp_list(3) = SEARCH
    interp_list(4) = 0
    interp_list(5) = 0

    maskval = 1.0
    masked  = 1
    fillval = 0.0
    msgval  = 0.0

    mask_array => landmask

    if(index(field % field, 'DMSO_JAN') /= 0 .or. &
       index(field % field, 'DMSO_FEB') /= 0 .or. &
       index(field % field, 'DMSO_MAR') /= 0 .or. &
       index(field % field, 'DMSO_APR') /= 0 .or. &
       index(field % field, 'DMSO_MAY') /= 0 .or. &
       index(field % field, 'DMSO_JUN') /= 0 .or. &
       index(field % field, 'DMSO_JUL') /= 0 .or. &
       index(field % field, 'DMSO_AUG') /= 0 .or. &
       index(field % field, 'DMSO_SEP') /= 0 .or. &
       index(field % field, 'DMSO_OCT') /= 0 .or. &
       index(field % field, 'DMSO_NOV') /= 0 .or. &
       index(field % field, 'DMSO_DEC') /= 0) then

       !
       !set up projection:
       !
       call map_init(proj)


       if(field%iproj == PROJ_LATLON) then
          call map_set(PROJ_LATLON,proj, &
                       latinc = real(field%deltalat,RKIND), &
                       loninc = real(field%deltalon,RKIND), &
                       knowni = 1.0_RKIND, &
                       knownj = 1.0_RKIND, &
                       lat1   = real(field%startlat,RKIND), &
                       lon1   = real(field%startlon,RKIND))
       elseif(field%iproj == PROJ_GAUSS) then
          call map_set(PROJ_GAUSS,proj, &
                       nlat = nint(field%deltalat), &
                       loninc = 360.0_RKIND / real(field%nx,RKIND), &
                       lat1 = real(field%startlat,RKIND), &
                       lon1 = real(field%startlon,RKIND))
       endif

       !
       !horizontally interpolate field at level k:
       !
       if(index(field%field,'DMSO_JAN') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DMSO_JAN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_backgrounds,'dms_gocart2G_clim',destField2d)
          nmonths = 1
       elseif(index(field%field,'DMSO_FEB') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DMSO_FEB at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_backgrounds,'dms_gocart2G_clim',destField2d)
          nmonths = 2
       elseif(index(field%field,'DMSO_MAR') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DMSO_MAR at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_backgrounds,'dms_gocart2G_clim',destField2d)
          nmonths = 3
       elseif(index(field%field,'DMSO_APR') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DMSO_APR at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_backgrounds,'dms_gocart2G_clim',destField2d)
          nmonths = 4
       elseif(index(field%field,'DMSO_MAY') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DMSO_MAY at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_backgrounds,'dms_gocart2G_clim',destField2d)
          nmonths = 5
       elseif(index(field%field,'DMSO_JUN') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DMSO_JUN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_backgrounds,'dms_gocart2G_clim',destField2d)
          nmonths = 6
       elseif(index(field%field,'DMSO_JUL') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DMSO_JUL at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_backgrounds,'dms_gocart2G_clim',destField2d)
          nmonths = 7
       elseif(index(field%field,'DMSO_AUG') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DMSO_AUG at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_backgrounds,'dms_gocart2G_clim',destField2d)
          nmonths = 8
       elseif(index(field%field,'DMSO_SEP') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DMSO_SEP at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_backgrounds,'dms_gocart2G_clim',destField2d)
          nmonths = 9
       elseif(index(field%field,'DMSO_OCT') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DMSO_OCT at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_backgrounds,'dms_gocart2G_clim',destField2d)
          nmonths = 10
       elseif(index(field%field,'DMSO_NOV') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DMSO_NOV at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_backgrounds,'dms_gocart2G_clim',destField2d)
          nmonths = 11
       elseif(index(field%field,'DMSO_DEC') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DMSO_DEC at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_backgrounds,'dms_gocart2G_clim',destField2d)
          nmonths = 12
       endif

       allocate(rslab(-2:field%nx+3,field%ny))
       rslab(1:field%nx,1:field%ny) = field%slab(1:field%nx,1:field%ny)
       rslab(0,1:field%ny)  = field%slab(field%nx  ,1:field%ny)
       rslab(-1,1:field%ny) = field%slab(field%nx-1,1:field%ny)
       rslab(-2,1:field%ny) = field%slab(field%nx-2,1:field%ny)
       rslab(field%nx+1,1:field%ny) = field%slab(1,1:field%ny)
       rslab(field%nx+2,1:field%ny) = field%slab(2,1:field%ny)
       rslab(field%nx+3,1:field%ny) = field%slab(3,1:field%ny)

       do iCell = 1, nInterpPoints
          if(mask_array(iCell) /= masked) then
             lat = latPoints(iCell)*DEG_PER_RAD
             lon = lonPoints(iCell)*DEG_PER_RAD
             call latlon_to_ij(proj,lat,lon,x,y)
             if(x < 0.5) then
                lon = lon + 360.0
                call latlon_to_ij(proj,lat,lon,x,y)
             elseif(x > real(field%nx,kind=RKIND)+ 0.5) then
                lon = lon - 360.0
                call latlon_to_ij(proj,lat,lon,x,y)
             endif

             if(maskval /= -1.0) then
                destField2d(nmonths,iCell) = interp_sequence(x,y,1,rslab,-2,field%nx+3,1,field%ny,1,1,msgval, &
                                                  interp_list,1,maskval=maskval,mask_array=maskslab)
             else
                destField2d(nmonths,iCell) = interp_sequence(x,y,1,rslab,-2,field%nx+3,1,field%ny,1,1,msgval, &
                                                  interp_list,1)
             endif
          else
             destField2d(nmonths,iCell) = fillval
          endif
       enddo
       deallocate(rslab)

    endif
    deallocate(field%slab)
    call read_next_met_field(field,istatus)

 enddo

 call read_met_close()

 call mpas_log_write('--- end subroutine init_hinterp_dms.')

 end subroutine init_hinterp_background_dms

!=================================================================================================================
 subroutine init_hinterp_background_sulf(filename_sulf,block,mesh,gocart2G_backgrounds)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 character(len=StrKIND),intent(in):: filename_sulf

!inout arguments:
 type(block_type),intent(inout),target:: block
 type(mpas_pool_type),intent(inout):: gocart2G_backgrounds

!local variables:
 type(dm_info),pointer:: dminfo
 type(met_data) :: field !real*4 meteorological data.
 type(proj_info):: proj

 logical:: have_landmask

 integer:: iCell,istatus,k,masked,nmonths,nInterpPoints
 integer,pointer:: nCells
 integer,dimension(5):: interp_list
 integer,dimension(:),pointer:: landmask
 integer,dimension(:),pointer:: mask_array

 real(kind=RKIND):: fillval,maskval,msgval
 real(kind=RKIND):: lat,lon,x,y
 real(kind=RKIND),dimension(:),pointer:: latCell,lonCell
 real(kind=RKIND),dimension(:),pointer:: latPoints,lonPoints
 real(kind=RKIND),dimension(:,:,:),pointer:: destField3d

 real(kind=RKIND),dimension(:,:),allocatable:: maskslab,rslab

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine init_hinterp_sulf:')

 dminfo => block%domain%dminfo

 call mpas_pool_get_dimension(mesh,'nCells',nCells)

 call mpas_pool_get_array(mesh,'landmask',landmask)
 call mpas_pool_get_array(mesh,'latCell' ,latCell )
 call mpas_pool_get_array(mesh,'lonCell' ,lonCell )


!open intermediate file:
 istatus = 0
 call read_met_init(trim(filename_sulf),.true.,'not needed',istatus)
 if(istatus /= 0) then
    call mpas_log_write('********************************************************************************')
    call mpas_log_write('Error opening gocart2G file '//trim(filename_sulf))
    call mpas_log_write('********************************************************************************')
    call mpas_dmpar_abort(dminfo)
 else
    call mpas_log_write('Processing file '//trim(filename_sulf))
 end if


!scan through all fields in the file, looking for the LANDSEA field:
 have_landmask = .false.
 call read_next_met_field(field,istatus)
 do while (istatus == 0)
    if(index(field % field, 'LANDSEA') /= 0) then
       have_landmask = .true.
       if(.not.allocated(maskslab)) allocate(maskslab(-2:field%nx+3,field%ny))


       maskslab(1:field%nx,1:field%ny) = field%slab(1:field%nx,1:field%ny)
       maskslab(0 ,1:field%ny) = field%slab(field%nx  ,1:field%ny)
       maskslab(-1,1:field%ny) = field%slab(field%nx-1,1:field%ny)
       maskslab(-2,1:field%ny) = field%slab(field%nx-2,1:field%ny)
       maskslab(field%nx+1,1:field%ny) = field%slab(1,1:field%ny)
       maskslab(field%nx+2,1:field%ny) = field%slab(2,1:field%ny)
       maskslab(field%nx+3,1:field%ny) = field%slab(3,1:field%ny)
!      call mpas_log_write('minval, maxval of LANDSEA = $r $r',realArgs=(/minval(maskslab),maxval(maskslab)/))
    end if
    deallocate(field%slab)
    call read_next_met_field(field,istatus)
 end do
 call read_met_close()

 if(.not. have_landmask) then
    call mpas_log_write('********************************************************************************')
    call mpas_log_write('Landsea mask not available from the surface file')
    call mpas_log_write('********************************************************************************')
 end if


!read climatological SULF data:
 istatus = 0
 call read_met_init(trim(filename_sulf),.true.,'not needed',istatus)
 if(istatus /= 0) then
    call mpas_log_write('********************************************************************************')
    call mpas_log_write('Error opening file '// trim(filename_sulf))
    call mpas_log_write('********************************************************************************')
    call mpas_dmpar_abort(dminfo)
 endif
 call read_next_met_field(field,istatus)


!horizontally interpolate DMS data:
 do while(istatus == 0)

    interp_list(1) = FOUR_POINT
    interp_list(2) = W_AVERAGE4
    interp_list(3) = W_AVERAGE16
    interp_list(4) = SEARCH
    interp_list(5) = 0

    maskval = -1.0
    masked  = -1
    fillval = 0.0
    msgval  = 0.0

    mask_array => landmask

    if(index(field % field, 'H2O2_JAN')  /= 0 .or. &
       index(field % field, 'H2O2_FEB')  /= 0 .or. &
       index(field % field, 'H2O2_MAR')  /= 0 .or. &
       index(field % field, 'H2O2_APR')  /= 0 .or. &
       index(field % field, 'H2O2_MAY')  /= 0 .or. &
       index(field % field, 'H2O2_JUN')  /= 0 .or. &
       index(field % field, 'H2O2_JUL')  /= 0 .or. &
       index(field % field, 'H2O2_AUG')  /= 0 .or. &
       index(field % field, 'H2O2_SEP')  /= 0 .or. &
       index(field % field, 'H2O2_OCT')  /= 0 .or. &
       index(field % field, 'H2O2_NOV')  /= 0 .or. &
       index(field % field, 'H2O2_DEC')  /= 0 .or. &
       index(field % field, 'NO3_JAN')   /= 0 .or. &
       index(field % field, 'NO3_FEB')   /= 0 .or. &
       index(field % field, 'NO3_MAR')   /= 0 .or. &
       index(field % field, 'NO3_APR')   /= 0 .or. &
       index(field % field, 'NO3_MAY')   /= 0 .or. &
       index(field % field, 'NO3_JUN')   /= 0 .or. &
       index(field % field, 'NO3_JUL')   /= 0 .or. &
       index(field % field, 'NO3_AUG')   /= 0 .or. &
       index(field % field, 'NO3_SEP')   /= 0 .or. &
       index(field % field, 'NO3_OCT')   /= 0 .or. &
       index(field % field, 'NO3_NOV')   /= 0 .or. &
       index(field % field, 'NO3_DEC')   /= 0 .or. &
       index(field % field, 'OH_JAN')    /= 0 .or. &
       index(field % field, 'OH_FEB')    /= 0 .or. &
       index(field % field, 'OH_MAR')    /= 0 .or. &
       index(field % field, 'OH_APR')    /= 0 .or. &
       index(field % field, 'OH_MAY')    /= 0 .or. &
       index(field % field, 'OH_JUN')    /= 0 .or. &
       index(field % field, 'OH_JUL')    /= 0 .or. &
       index(field % field, 'OH_AUG')    /= 0 .or. &
       index(field % field, 'OH_SEP')    /= 0 .or. &
       index(field % field, 'OH_OCT')    /= 0 .or. &
       index(field % field, 'OH_NOV')    /= 0 .or. &
       index(field % field, 'OH_DEC')    /= 0 .or. &
       index(field % field, 'PRES_JAN')  /= 0 .or. &
       index(field % field, 'PRES_FEB')  /= 0 .or. &
       index(field % field, 'PRES_MAR')  /= 0 .or. &
       index(field % field, 'PRES_APR')  /= 0 .or. &
       index(field % field, 'PRES_MAY')  /= 0 .or. &
       index(field % field, 'PRES_JUN')  /= 0 .or. &
       index(field % field, 'PRES_JUL')  /= 0 .or. &
       index(field % field, 'PRES_AUG')  /= 0 .or. &
       index(field % field, 'PRES_SEP')  /= 0 .or. &
       index(field % field, 'PRES_OCT')  /= 0 .or. &
       index(field % field, 'PRES_NOV')  /= 0 .or. &
       index(field % field, 'PRES_DEC')  /= 0 .or. &
       index(field % field, 'DPRES_JAN') /= 0 .or. &
       index(field % field, 'DPRES_FEB') /= 0 .or. &
       index(field % field, 'DPRES_MAR') /= 0 .or. &
       index(field % field, 'DPRES_APR') /= 0 .or. &
       index(field % field, 'DPRES_MAY') /= 0 .or. &
       index(field % field, 'DPRES_JUN') /= 0 .or. &
       index(field % field, 'DPRES_JUL') /= 0 .or. &
       index(field % field, 'DPRES_AUG') /= 0 .or. &
       index(field % field, 'DPRES_SEP') /= 0 .or. &
       index(field % field, 'DPRES_OCT') /= 0 .or. &
       index(field % field, 'DPRES_NOV') /= 0 .or. &
       index(field % field, 'DPRES_DEC') /= 0 ) then

       !
       !set up projection:
       !
       call map_init(proj)


       if(field%iproj == PROJ_LATLON) then
          call map_set(PROJ_LATLON,proj, &
                       latinc = real(field%deltalat,RKIND), &
                       loninc = real(field%deltalon,RKIND), &
                       knowni = 1.0_RKIND, &
                       knownj = 1.0_RKIND, &
                       lat1   = real(field%startlat,RKIND), &
                       lon1   = real(field%startlon,RKIND))
       elseif(field%iproj == PROJ_GAUSS) then
          call map_set(PROJ_GAUSS,proj, &
                       nlat = nint(field%deltalat), &
                       loninc = 360.0_RKIND / real(field%nx,RKIND), &
                       lat1 = real(field%startlat,RKIND), &
                       lon1 = real(field%startlon,RKIND))
       endif

       !
       !horizontally interpolate field at level k:
       !
       if(index(field%field,'H2O2_JAN') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating H2O2_JAN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_backgrounds,'h2o2_gocart2G_clim',destField3d)
          nmonths = 1
       elseif(index(field%field,'H2O2_FEB') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating H2O2_FEB at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_backgrounds,'h2o2_gocart2G_clim',destField3d)
          nmonths = 2
       elseif(index(field%field,'H2O2_MAR') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating H2O2_MAR at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_backgrounds,'h2o2_gocart2G_clim',destField3d)
          nmonths = 3
       elseif(index(field%field,'H2O2_APR') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating H2O2_APR at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_backgrounds,'h2o2_gocart2G_clim',destField3d)
          nmonths = 4
       elseif(index(field%field,'H2O2_MAY') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating H2O2_MAY at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_backgrounds,'h2o2_gocart2G_clim',destField3d)
          nmonths = 5
       elseif(index(field%field,'H2O2_JUN') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating H2O2_JUN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_backgrounds,'h2o2_gocart2G_clim',destField3d)
          nmonths = 6
       elseif(index(field%field,'H2O2_JUL') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating H2O2_JUL at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_backgrounds,'h2o2_gocart2G_clim',destField3d)
          nmonths = 7
       elseif(index(field%field,'H2O2_AUG') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating H2O2_AUG at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_backgrounds,'h2o2_gocart2G_clim',destField3d)
          nmonths = 8
       elseif(index(field%field,'H2O2_SEP') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating H2O2_SEP at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_backgrounds,'h2o2_gocart2G_clim',destField3d)
          nmonths = 9
       elseif(index(field%field,'H2O2_OCT') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating H2O2_OCT at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_backgrounds,'h2o2_gocart2G_clim',destField3d)
          nmonths = 10
       elseif(index(field%field,'H2O2_NOV') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating H2O2_NOV at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_backgrounds,'h2o2_gocart2G_clim',destField3d)
          nmonths = 11
       elseif(index(field%field,'H2O2_DEC') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating H2O2_DEC at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_backgrounds,'h2o2_gocart2G_clim',destField3d)
          nmonths = 12
       elseif(index(field%field,'NO3_JAN') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating NO3_JAN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_backgrounds,'no3_gocart2G_clim',destField3d)
          nmonths = 1
       elseif(index(field%field,'NO3_FEB') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating NO3_FEB at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'no3_gocart2G_clim',destField3d)
           nmonths = 2
       elseif(index(field%field,'NO3_MAR') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating NO3_MAR at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'no3_gocart2G_clim',destField3d)
           nmonths = 3
       elseif(index(field%field,'NO3_APR') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating NO3_APR at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'no3_gocart2G_clim',destField3d)
           nmonths = 4
       elseif(index(field%field,'NO3_MAY') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating NO3_MAY at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'no3_gocart2G_clim',destField3d)
           nmonths = 5
       elseif(index(field%field,'NO3_JUN') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating NO3_JUN at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'no3_gocart2G_clim',destField3d)
           nmonths = 6
       elseif(index(field%field,'NO3_JUL') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating NO3_JUL at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'no3_gocart2G_clim',destField3d)
           nmonths = 7
       elseif(index(field%field,'NO3_AUG') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating NO3_AUG at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'no3_gocart2G_clim',destField3d)
           nmonths = 8
       elseif(index(field%field,'NO3_SEP') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating NO3_SEP at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'no3_gocart2G_clim',destField3d)
           nmonths = 9
       elseif(index(field%field,'NO3_OCT') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating NO3_OCT at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'no3_gocart2G_clim',destField3d)
           nmonths = 10
       elseif(index(field%field,'NO3_NOV') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating NO3_NOV at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'no3_gocart2G_clim',destField3d)
           nmonths = 11
       elseif(index(field%field,'NO3_DEC') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating NO3_DEC at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'no3_gocart2G_clim',destField3d)
           nmonths = 12
       elseif(index(field%field,'OH_JAN') /= 0) then 
          k = field%xlvl
          call mpas_log_write('Interpolating OH_JAN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_backgrounds,'oh_gocart2G_clim',destField3d)
          nmonths = 1
       elseif(index(field%field,'OH_FEB') /= 0) then 
           k = field%xlvl
           call mpas_log_write('Interpolating OH_FEB at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'oh_gocart2G_clim',destField3d)
           nmonths = 2
       elseif(index(field%field,'OH_MAR') /= 0) then 
           k = field%xlvl
           call mpas_log_write('Interpolating OH_MAR at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'oh_gocart2G_clim',destField3d)
           nmonths = 3
       elseif(index(field%field,'OH_APR') /= 0) then 
           k = field%xlvl
           call mpas_log_write('Interpolating OH_APR at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'oh_gocart2G_clim',destField3d)
           nmonths = 4
       elseif(index(field%field,'OH_MAY') /= 0) then 
           k = field%xlvl
           call mpas_log_write('Interpolating OH_MAY at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'oh_gocart2G_clim',destField3d)
           nmonths = 5
       elseif(index(field%field,'OH_JUN') /= 0) then 
           k = field%xlvl
           call mpas_log_write('Interpolating OH_JUN at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'oh_gocart2G_clim',destField3d)
           nmonths = 6
       elseif(index(field%field,'OH_JUL') /= 0) then 
           k = field%xlvl
           call mpas_log_write('Interpolating OH_JUL at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'oh_gocart2G_clim',destField3d)
           nmonths = 7
       elseif(index(field%field,'OH_AUG') /= 0) then 
           k = field%xlvl
           call mpas_log_write('Interpolating OH_AUG at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'oh_gocart2G_clim',destField3d)
           nmonths = 8
       elseif(index(field%field,'OH_SEP') /= 0) then 
           k = field%xlvl
           call mpas_log_write('Interpolating OH_SEP at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'oh_gocart2G_clim',destField3d)
           nmonths = 9
       elseif(index(field%field,'OH_OCT') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating OH_OCT at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'oh_gocart2G_clim',destField3d)
           nmonths = 10
       elseif(index(field%field,'OH_NOV') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating OH_NOV at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'oh_gocart2G_clim',destField3d)
           nmonths = 11
       elseif(index(field%field,'OH_DEC') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating OH_DEC at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'oh_gocart2G_clim',destField3d)
           nmonths = 12
       elseif(index(field%field,'PRES_JAN') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating PRES_JAN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_backgrounds,'pres_gocart2G_clim',destField3d)
          nmonths = 1
       elseif(index(field%field,'PRES_FEB') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating PRES_FEB at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'pres_gocart2G_clim',destField3d)
           nmonths = 2
       elseif(index(field%field,'PRES_MAR') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating PRES_MAR at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'pres_gocart2G_clim',destField3d)
           nmonths = 3
       elseif(index(field%field,'PRES_APR') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating PRES_APR at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'pres_gocart2G_clim',destField3d)
           nmonths = 4
       elseif(index(field%field,'PRES_MAY') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating PRES_MAY at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'pres_gocart2G_clim',destField3d)
           nmonths = 5
       elseif(index(field%field,'PRES_JUN') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating PRES_JUN at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'pres_gocart2G_clim',destField3d)
           nmonths = 6
       elseif(index(field%field,'PRES_JUL') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating PRES_JUL at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'pres_gocart2G_clim',destField3d)
           nmonths = 7
       elseif(index(field%field,'PRES_AUG') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating PRES_AUG at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'pres_gocart2G_clim',destField3d)
           nmonths = 8
       elseif(index(field%field,'PRES_SEP') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating PRES_SEP at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'pres_gocart2G_clim',destField3d)
           nmonths = 9
       elseif(index(field%field,'PRES_OCT') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating PRES_OCT at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'pres_gocart2G_clim',destField3d)
           nmonths = 10
       elseif(index(field%field,'PRES_NOV') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating PRES_NOV at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'pres_gocart2G_clim',destField3d)
           nmonths = 11
       elseif(index(field%field,'PRES_DEC') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating PRES_DEC at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'pres_gocart2G_clim',destField3d)
           nmonths = 12
       elseif(index(field%field,'DPRES_JAN') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DPRES_JAN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(gocart2G_backgrounds,'dpres_gocart2G_clim',destField3d)
          nmonths = 1
       elseif(index(field%field,'DPRES_FEB') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating DPRES_FEB at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'dpres_gocart2G_clim',destField3d)
           nmonths = 2
       elseif(index(field%field,'DPRES_MAR') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating DPRES_MAR at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'dpres_gocart2G_clim',destField3d)
           nmonths = 3
       elseif(index(field%field,'DPRES_APR') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating DPRES_APR at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'dpres_gocart2G_clim',destField3d)
           nmonths = 4
       elseif(index(field%field,'DPRES_MAY') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating DPRES_MAY at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'dpres_gocart2G_clim',destField3d)
           nmonths = 5
       elseif(index(field%field,'DPRES_JUN') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating DPRES_JUN at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'dpres_gocart2G_clim',destField3d)
           nmonths = 6
       elseif(index(field%field,'DPRES_JUL') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating DPRES_JUL at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'dpres_gocart2G_clim',destField3d)
           nmonths = 7
       elseif(index(field%field,'DPRES_AUG') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating DPRES_AUG at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'dpres_gocart2G_clim',destField3d)
           nmonths = 8
       elseif(index(field%field,'DPRES_SEP') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating DPRES_SEP at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'dpres_gocart2G_clim',destField3d)
           nmonths = 9
       elseif(index(field%field,'DPRES_OCT') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating DPRES_OCT at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'dpres_gocart2G_clim',destField3d)
           nmonths = 10
       elseif(index(field%field,'DPRES_NOV') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating DPRES_NOV at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'dpres_gocart2G_clim',destField3d)
           nmonths = 11
       elseif(index(field%field,'DPRES_DEC') /= 0) then
           k = field%xlvl
           call mpas_log_write('Interpolating DPRES_DEC at $i',intArgs=(/k/))
           nInterpPoints = nCells
           latPoints => latCell
           lonPoints => lonCell
           call mpas_pool_get_array(gocart2G_backgrounds,'dpres_gocart2G_clim',destField3d)
           nmonths = 12
       endif

       allocate(rslab(-2:field%nx+3,field%ny))
       rslab(1:field%nx,1:field%ny) = field%slab(1:field%nx,1:field%ny)
       rslab(0,1:field%ny)  = field%slab(field%nx  ,1:field%ny)
       rslab(-1,1:field%ny) = field%slab(field%nx-1,1:field%ny)
       rslab(-2,1:field%ny) = field%slab(field%nx-2,1:field%ny)
       rslab(field%nx+1,1:field%ny) = field%slab(1,1:field%ny)
       rslab(field%nx+2,1:field%ny) = field%slab(2,1:field%ny)
       rslab(field%nx+3,1:field%ny) = field%slab(3,1:field%ny)

       do iCell = 1, nInterpPoints
          if(mask_array(iCell) /= masked) then
             lat = latPoints(iCell)*DEG_PER_RAD
             lon = lonPoints(iCell)*DEG_PER_RAD
             call latlon_to_ij(proj,lat,lon,x,y)
             if(x < 0.5) then
                lon = lon + 360.0
                call latlon_to_ij(proj,lat,lon,x,y)
             elseif(x > real(field%nx,kind=RKIND)+ 0.5) then
                lon = lon - 360.0
                call latlon_to_ij(proj,lat,lon,x,y)
             endif

             if(maskval /= -1.0) then
                destField3d(nmonths,k,iCell) = interp_sequence(x,y,1,rslab,-2,field%nx+3,1,field%ny,1,1,msgval, &
                                                               interp_list,1,maskval=maskval,mask_array=maskslab)
             else
                destField3d(nmonths,k,iCell) = interp_sequence(x,y,1,rslab,-2,field%nx+3,1,field%ny,1,1,msgval, &
                                                               interp_list,1)
             endif
          else
             destField3d(nmonths,k,iCell) = fillval
          endif
       enddo
       deallocate(rslab)

    endif
    deallocate(field%slab)
    call read_next_met_field(field,istatus)

 enddo

 call read_met_close()

 call mpas_log_write('--- end subroutine init_hinterp_sulf.')

 end subroutine init_hinterp_background_sulf

!=================================================================================================================
 end module mpas_init_chem_gocart2G_climatology
!=================================================================================================================

