! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
 module mpas_init_atm_gocart
 use mpas_derived_types
 use mpas_kind_types
 use mpas_log
 use mpas_dmpar
 use mpas_pool_routines

 use init_atm_read_met
 use init_atm_hinterp
 use init_atm_llxy
 use init_atm_vinterp
 use mpas_atmphys_date_time
 use mpas_atmphys_utilities

 implicit none
 private
 public:: init_atm_gocart


!mpas_init_atm_gocart contains the subroutines needed for the interpolation of climatological
!monthly-averaged hygroscopic ("water friendly") and nonhygroscopic ("ice friendly") aerosols
!used in the Thompson parameterization of cloud microphysics with CCN and IN nucleation (use
!physics option mp_thompson_aerorosols=true).
!Laura D. Fowler (laura@ucar.edu) / 2016-05-27.

! add-ons and modifications to sourcecode:
! ----------------------------------------
! * in subroutines init_hinterp_gocart and init_vinterp_gocart, added the horizontal and vertical interpolations
!   of the individual aerosols that constitute the "water-friendly" (OC,SO4,SEAS1,and SEAS2) and "ice_friendly"
!   (DUST1,DUST2,DUST3,DUST4,and DUST5) aerosols so that we can compute activation and wet scavenging processes
!   during the convective transport of aerosols in the Grell-Freitas convection scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2021-02-17.  


 contains


!=================================================================================================================
 subroutine init_atm_gocart(block,mesh,configs,diag,state)
!=================================================================================================================

!input arguments:
 type (mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: diag

!inout arguments:
 type(block_type),intent(inout),target:: block 
 type(mpas_pool_type),intent(inout)   :: mesh
 type(mpas_pool_type),intent(inout)   :: state
!local variables and pointers:
 character (len=StrKIND),pointer:: config_start_time
 character(len=StrKIND):: filename_gocart
 character(len=StrKIND):: initial_date,mess

 logical:: lexist

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine init_atm_gocart:')

!inquire if the GOCART input file exists:
 filename_gocart = "QNWFA_QNIFA_Monthly_Sigma.dat"
!filename_gocart = "QNWFA_QNIFA_Monthly_Sigma.withAerosols.dat"

 inquire(file=filename_gocart,exist=lexist)
 if(lexist) then

    call mpas_pool_get_config(configs,'config_start_time',config_start_time)

    !--- horizontal interpolation of the climatological monthly-averaged GOCART data to the MPAS mesh:  
    call init_hinterp_gocart(block,mesh)

    !--- interpolation of the monthly-averaged GOCART data to the initial date, and vertical interpolation to
    !    the MPAS levels:
    initial_date = trim(config_start_time)
    call init_vinterp_gocart(initial_date,mesh,diag,state)
 else
    call mpas_log_write('QNWFA_QNIFA_Monthly_Sigma.dat was not found in local directory:')
!   call mpas_log_write('QNWFA_QNIFA_Monthly_Sigma.withAerosols.dat was not found in local directory:')
    call mpas_log_write('nwfa and nifa are set to zero and not interpolated from climatological data.')
 endif

 call mpas_log_write('--- end subroutine init_atm_gocart.')
 call mpas_log_write(' ')

 end subroutine init_atm_gocart

!=================================================================================================================
 subroutine init_vinterp_gocart(initial_date,mesh,diag,state)
!=================================================================================================================

!input arguments:
 character(len=StrKIND),intent(in):: initial_date
 type(mpas_pool_type),intent(in):: diag

!inout arguments:
 type(mpas_pool_type),intent(inout):: mesh
 type(mpas_pool_type),intent(inout):: state

!local variables and pointers:
 integer,pointer:: nCells,nGocartLevels,nVertLevels,nMonths
 integer,pointer:: index_nifa,index_nwfa
 integer,pointer:: num_aerosols
 integer:: iCell,k,kk,n

 real(kind=RKIND),dimension(:,:),pointer  :: nifa,nwfa,pressure
 real(kind=RKIND),dimension(:,:,:),pointer:: nifa_clim,nwfa_clim,pwif_clim
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars
 real(kind=RKIND),dimension(:,:,:),pointer:: aerosols,aerosols_weight
 real(kind=RKIND),dimension(:,:,:,:),pointer:: aerosols_clim

 real(kind=RKIND):: target_p
 real(kind=RKIND),dimension(:,:),allocatable:: nifa_int,nwfa_int,pwif_int,aers_int,sorted_arr

 real(kind=RKIND),dimension(:),allocatable::   dummy2
 real(kind=RKIND),dimension(:,:),allocatable:: dummy1

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine init_vinterp_gocart:')

 call mpas_pool_get_dimension(mesh,'nCells'       ,nCells       )
 call mpas_pool_get_dimension(mesh,'nGocartLevels',nGocartLevels)
 call mpas_pool_get_dimension(mesh,'nVertLevels'  ,nVertLevels  )
 call mpas_pool_get_dimension(mesh,'nMonths'      ,nMonths      )

 call mpas_pool_get_dimension(state,'index_nifa',index_nifa)
 call mpas_pool_get_dimension(state,'index_nwfa',index_nwfa)
    
 call mpas_pool_get_array(diag,'pressure_base',pressure)

 call mpas_pool_get_array(mesh,'nifa_clim',nifa_clim)
 call mpas_pool_get_array(mesh,'nwfa_clim',nwfa_clim)
 call mpas_pool_get_array(mesh,'pwif_clim',pwif_clim)

 call mpas_pool_get_array(state,'scalars',scalars)
 nifa => scalars(index_nifa,:,:)
 nwfa => scalars(index_nwfa,:,:)

 if(.not.allocated(nifa_int)  ) allocate(nifa_int(nGocartLevels,nCells))
 if(.not.allocated(nwfa_int)  ) allocate(nwfa_int(nGocartLevels,nCells))
 if(.not.allocated(pwif_int)  ) allocate(pwif_int(nGocartLevels,nCells))
 if(.not.allocated(sorted_arr)) allocate(sorted_arr(2,nGocartLevels))

!--- interpolation of the monthly-averaged GOCART data to the initial date, and vertical interpolation to the
!    MPAS levels:
 if(.not.allocated(dummy2)) allocate(dummy2(nCells))
 if(.not.allocated(dummy1)) allocate(dummy1(nMonths,nCells))

 do k = 1, nGocartLevels
    dummy2(1:nCells) = 0._RKIND
    dummy1(1:nMonths,1:nCells) = pwif_clim(1:nMonths,k,1:nCells)
    call monthly_interp_to_date(nCells,initial_date,dummy1,dummy2)
    pwif_int(k,1:nCells) = dummy2(1:nCells)
 enddo

!--- nifa:
 do k = 1, nGocartLevels
    dummy2(1:nCells) = 0._RKIND
    dummy1(1:nMonths,1:nCells) = nifa_clim(1:nMonths,k,1:nCells)
    call monthly_interp_to_date(nCells,initial_date,dummy1,dummy2)
    nifa_int(k,1:nCells) = dummy2(1:nCells)
 enddo
 do iCell = 1, nCells
    sorted_arr(1,1:nGocartLevels) = 0._RKIND
    sorted_arr(1,1:nGocartLevels) = 0._RKIND
    do k = 1, nGocartLevels
       kk = nGocartLevels + 1 -k
       sorted_arr(1,kk) = pwif_int(k,iCell)
       sorted_arr(2,kk) = nifa_int(k,iCell)
    enddo
    do k = nVertLevels, 1, -1
       target_p = pressure(k,iCell)
       nifa(k,iCell) = vertical_interp(target_p,nGocartLevels-1, &
                                sorted_arr(:,1:nGocartLevels-1),order=1,extrap=0)
       if(target_p >= pwif_int(1,iCell) .and. k < nVertLevels) nifa(k,iCell) = nifa(k+1,iCell)
    enddo
 enddo

!--- nwfa:
 do k = 1, nGocartLevels
    dummy2(1:nCells) = 0._RKIND
    dummy1(1:nMonths,1:nCells) = nwfa_clim(1:nMonths,k,1:nCells)
    call monthly_interp_to_date(nCells,initial_date,dummy1,dummy2)
    nwfa_int(k,1:nCells) = dummy2(1:nCells)
 enddo
 do iCell = 1, nCells
    sorted_arr(1,1:nGocartLevels) = 0._RKIND
    sorted_arr(1,1:nGocartLevels) = 0._RKIND
    do k = 1, nGocartLevels
       kk = nGocartLevels + 1 -k
       sorted_arr(1,kk) = pwif_int(k,iCell)
       sorted_arr(2,kk) = nwfa_int(k,iCell)
    enddo
    do k = nVertLevels, 1, -1
       target_p = pressure(k,iCell)
       nwfa(k,iCell) = vertical_interp(target_p,nGocartLevels-1, &
                                sorted_arr(:,1:nGocartLevels-1),order=1,extrap=0)
       if(target_p >= pwif_int(1,iCell) .and. k < nVertLevels) nwfa(k,iCell) = nwfa(k+1,iCell)
    enddo
 enddo

!--- contribution of individual species of aerosols to nwfa and nifa, and calculation of their
!    respective weight to nwfa and nifa:
 call mpas_pool_get_dimension(mesh,'num_aerosols',num_aerosols)
 call mpas_pool_get_array(mesh,'aerosols_clim'  ,aerosols_clim  )
 call mpas_pool_get_array(mesh,'aerosols'       ,aerosols       )
 call mpas_pool_get_array(mesh,'aerosols_weight',aerosols_weight)
 call mpas_log_write(' ')
 call mpas_log_write('--- num_aerosols = $i',intArgs=(/num_aerosols/))

 if(.not.allocated(aers_int)) allocate(aers_int(nGocartLevels,nCells))

 do n = 1, num_aerosols
    do k = 1, nGocartLevels
       dummy2(1:nCells) = 0._RKIND
       dummy1(1:nMonths,1:nCells) = aerosols_clim(n,1:nMonths,k,1:nCells)
       call monthly_interp_to_date(nCells,initial_date,dummy1,dummy2)
       aers_int(k,1:nCells) = dummy2(1:nCells)
    enddo

    do iCell = 1, nCells
       sorted_arr(1,1:nGocartLevels) = 0._RKIND
       sorted_arr(1,1:nGocartLevels) = 0._RKIND
       do k = 1, nGocartLevels
          kk = nGocartLevels + 1 -k
          sorted_arr(1,kk) = pwif_int(k,iCell)
          sorted_arr(2,kk) = aers_int(k,iCell)
       enddo
       do k = nVertLevels, 1, -1
          target_p = pressure(k,iCell)
          aerosols(n,k,iCell) = vertical_interp(target_p,nGocartLevels-1, &
                                    sorted_arr(:,1:nGocartLevels-1),order=1,extrap=0)
          if(target_p >= pwif_int(1,iCell) .and. k < nVertLevels) &
             aerosols(n,k,iCell) = aerosols(n,k+1,iCell)
       enddo
    enddo
 enddo

 do n = 1, num_aerosols
    if(n .le. 5) then
       do k = 1, nVertLevels
          do iCell = 1, nCells
             if(nifa(k,iCell) .ge. 1.e-10) then
                aerosols_weight(n,k,iCell) = aerosols(n,k,iCell)/nifa(k,iCell)
             else
                aerosols_weight(n,k,iCell) = 0._RKIND
             endif
          enddo
       enddo
    elseif(n .gt. 5) then
       do k = 1, nVertLevels
          do iCell = 1, nCells
             if(nwfa(k,iCell) .ge. 1.e-10) then
                aerosols_weight(n,k,iCell) = aerosols(n,k,iCell)/nwfa(k,iCell)
             else
                aerosols_weight(n,k,iCell) = 0._RKIND
             endif
          enddo
       enddo
    endif

    do k = 1, nVertLevels
       call mpas_log_write('$i $i $r $r',intArgs=(/n,k/),realArgs=(/maxval(aerosols_weight(n,k,1:nCells)), &
                           minval(aerosols_weight(n,k,1:nCells))/))
    enddo
    call mpas_log_write(' ')
 enddo

!--- deallocation of local arrays:
 if(allocated(dummy1)    ) deallocate(dummy1    )
 if(allocated(dummy2)    ) deallocate(dummy2    )
 if(allocated(aers_int)  ) deallocate(aers_int  )
 if(allocated(nifa_int)  ) deallocate(nifa_int  )
 if(allocated(nwfa_int)  ) deallocate(nwfa_int  )
 if(allocated(pwif_int)  ) deallocate(pwif_int  )
 if(allocated(sorted_arr)) deallocate(sorted_arr)
 
 call mpas_log_write('--- enter subroutine init_vinterp_gocart:')

 end subroutine init_vinterp_gocart

!=================================================================================================================
 subroutine init_hinterp_gocart(block,mesh)
!=================================================================================================================

!inout arguments:
 type(block_type),intent(inout),target:: block 
 type (mpas_pool_type),intent(inout)  :: mesh

!local variables:
 type(dm_info),pointer:: dminfo
 type(met_data) :: field !real*4 meteorological data.
 type(proj_info):: proj

 character(len=StrKIND):: filename_gocart
 logical:: have_landmask

 integer,pointer:: nCells
 integer:: i,j
 integer:: iCell,istatus,k,masked,nmonths,nInterpPoints
 integer,dimension(5):: interp_list
 integer,pointer:: index_dust1,index_dust2,index_dust3,index_dust4,index_dust5
 integer,pointer:: index_oc,index_so4,index_seas1,index_seas2
 integer,dimension(:),pointer:: landmask
 integer,dimension(:),pointer:: mask_array

 real(kind=RKIND):: fillval,maskval,msgval
 real(kind=RKIND):: lat,lon,x,y
 real(kind=RKIND),dimension(:),pointer    :: latCell,lonCell
 real(kind=RKIND),dimension(:),pointer    :: latPoints,lonPoints
 real(kind=RKIND),dimension(:,:,:),pointer:: nifa_clim,nwfa_clim,pwif_clim
 real(kind=RKIND),dimension(:,:,:),pointer:: dust1_clim,dust2_clim,dust3_clim,dust4_clim,dust5_clim
 real(kind=RKIND),dimension(:,:,:),pointer:: oc_clim,so4_clim,seas1_clim,seas2_clim
 real(kind=RKIND),dimension(:,:,:),pointer:: destField3d
 real(kind=RKIND),dimension(:,:,:,:),pointer:: aerosols_clim

 real(kind=RKIND),dimension(:,:),allocatable:: maskslab,rslab

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine init_hinterp_gocart:')

 dminfo => block%domain%dminfo

 filename_gocart = "QNWFA_QNIFA_Monthly_Sigma.dat"
!filename_gocart = "QNWFA_QNIFA_Monthly_Sigma.withAerosols.dat"

 call mpas_pool_get_dimension(mesh,'nCells',nCells)

 call mpas_pool_get_array(mesh,'landmask',landmask)
 call mpas_pool_get_array(mesh,'latCell' ,latCell )
 call mpas_pool_get_array(mesh,'lonCell' ,lonCell )

 call mpas_pool_get_array(mesh,'nifa_clim',nifa_clim)
 call mpas_pool_get_array(mesh,'nwfa_clim',nwfa_clim)
 call mpas_pool_get_array(mesh,'pwif_clim',pwif_clim)

 call mpas_pool_get_dimension(mesh,'index_dust1_clim',index_dust1)
 call mpas_pool_get_dimension(mesh,'index_dust2_clim',index_dust2)
 call mpas_pool_get_dimension(mesh,'index_dust3_clim',index_dust3)
 call mpas_pool_get_dimension(mesh,'index_dust4_clim',index_dust4)
 call mpas_pool_get_dimension(mesh,'index_dust5_clim',index_dust5)
 call mpas_pool_get_dimension(mesh,'index_oc_clim'   ,index_oc   )
 call mpas_pool_get_dimension(mesh,'index_so4_clim'  ,index_so4  )
 call mpas_pool_get_dimension(mesh,'index_seas1_clim',index_seas1)
 call mpas_pool_get_dimension(mesh,'index_seas2_clim',index_seas2)
 call mpas_pool_get_array(mesh,'aerosols_clim',aerosols_clim)

 dust1_clim => aerosols_clim(index_dust1,:,:,:)
 dust2_clim => aerosols_clim(index_dust2,:,:,:)
 dust3_clim => aerosols_clim(index_dust3,:,:,:)
 dust4_clim => aerosols_clim(index_dust4,:,:,:)
 dust5_clim => aerosols_clim(index_dust5,:,:,:)
 oc_clim    => aerosols_clim(index_oc,:,:,:   )
 so4_clim   => aerosols_clim(index_so4,:,:,:  )
 seas1_clim => aerosols_clim(index_seas1,:,:,:)
 seas2_clim => aerosols_clim(index_seas2,:,:,:)

!open intermediate file:
 istatus = 0
 call read_gocart_init(trim(filename_gocart),istatus)
 if(istatus /= 0) then
    call mpas_log_write('********************************************************************************')
    call mpas_log_write('Error opening gocart file '//trim(filename_gocart))
    call mpas_log_write('********************************************************************************')
    call mpas_dmpar_abort(dminfo)
 else
    call mpas_log_write('Processing file '//trim(filename_gocart))
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


!read gocart data:
 istatus = 0
 call read_gocart_init(trim(filename_gocart),istatus)
 if(istatus /= 0) then
    call mpas_log_write('********************************************************************************')
    call mpas_log_write('Error opening gocart file '// trim(filename_gocart))
    call mpas_log_write('********************************************************************************')
    call mpas_dmpar_abort(dminfo)
 endif
 call read_next_met_field(field, istatus)

!horizontally interpolate GOCART data:
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

    if(index(field % field, 'QNIFA_JAN') /= 0 .or. &
       index(field % field, 'QNIFA_FEB') /= 0 .or. &
       index(field % field, 'QNIFA_MAR') /= 0 .or. &
       index(field % field, 'QNIFA_APR') /= 0 .or. &
       index(field % field, 'QNIFA_MAY') /= 0 .or. &
       index(field % field, 'QNIFA_JUN') /= 0 .or. &
       index(field % field, 'QNIFA_JUL') /= 0 .or. &
       index(field % field, 'QNIFA_AUG') /= 0 .or. &
       index(field % field, 'QNIFA_SEP') /= 0 .or. &
       index(field % field, 'QNIFA_OCT') /= 0 .or. &
       index(field % field, 'QNIFA_NOV') /= 0 .or. &
       index(field % field, 'QNIFA_DEC') /= 0 .or. &
       index(field % field, 'QNWFA_JAN') /= 0 .or. &
       index(field % field, 'QNWFA_FEB') /= 0 .or. &
       index(field % field, 'QNWFA_MAR') /= 0 .or. &
       index(field % field, 'QNWFA_APR') /= 0 .or. &
       index(field % field, 'QNWFA_MAY') /= 0 .or. &
       index(field % field, 'QNWFA_JUN') /= 0 .or. &
       index(field % field, 'QNWFA_JUL') /= 0 .or. &
       index(field % field, 'QNWFA_AUG') /= 0 .or. &
       index(field % field, 'QNWFA_SEP') /= 0 .or. &
       index(field % field, 'QNWFA_OCT') /= 0 .or. &
       index(field % field, 'QNWFA_NOV') /= 0 .or. &
       index(field % field, 'QNWFA_DEC') /= 0 .or. &
       index(field % field, 'DU1_JAN'  ) /= 0 .or. &
       index(field % field, 'DU1_FEB'  ) /= 0 .or. &
       index(field % field, 'DU1_MAR'  ) /= 0 .or. &
       index(field % field, 'DU1_APR'  ) /= 0 .or. &
       index(field % field, 'DU1_MAY'  ) /= 0 .or. &
       index(field % field, 'DU1_JUN'  ) /= 0 .or. &
       index(field % field, 'DU1_JUL'  ) /= 0 .or. &
       index(field % field, 'DU1_AUG'  ) /= 0 .or. &
       index(field % field, 'DU1_SEP'  ) /= 0 .or. &
       index(field % field, 'DU1_OCT'  ) /= 0 .or. &
       index(field % field, 'DU1_NOV'  ) /= 0 .or. &
       index(field % field, 'DU1_DEC'  ) /= 0 .or. &
       index(field % field, 'DU2_JAN'  ) /= 0 .or. &
       index(field % field, 'DU2_FEB'  ) /= 0 .or. &
       index(field % field, 'DU2_MAR'  ) /= 0 .or. &
       index(field % field, 'DU2_APR'  ) /= 0 .or. &
       index(field % field, 'DU2_MAY'  ) /= 0 .or. &
       index(field % field, 'DU2_JUN'  ) /= 0 .or. &
       index(field % field, 'DU2_JUL'  ) /= 0 .or. &
       index(field % field, 'DU2_AUG'  ) /= 0 .or. &
       index(field % field, 'DU2_SEP'  ) /= 0 .or. &
       index(field % field, 'DU2_OCT'  ) /= 0 .or. &
       index(field % field, 'DU2_NOV'  ) /= 0 .or. &
       index(field % field, 'DU2_DEC'  ) /= 0 .or. &
       index(field % field, 'DU3_JAN'  ) /= 0 .or. &
       index(field % field, 'DU3_FEB'  ) /= 0 .or. &
       index(field % field, 'DU3_MAR'  ) /= 0 .or. &
       index(field % field, 'DU3_APR'  ) /= 0 .or. &
       index(field % field, 'DU3_MAY'  ) /= 0 .or. &
       index(field % field, 'DU3_JUN'  ) /= 0 .or. &
       index(field % field, 'DU3_JUL'  ) /= 0 .or. &
       index(field % field, 'DU3_AUG'  ) /= 0 .or. &
       index(field % field, 'DU3_SEP'  ) /= 0 .or. &
       index(field % field, 'DU3_OCT'  ) /= 0 .or. &
       index(field % field, 'DU3_NOV'  ) /= 0 .or. &
       index(field % field, 'DU3_DEC'  ) /= 0 .or. &
       index(field % field, 'DU4_JAN'  ) /= 0 .or. &
       index(field % field, 'DU4_FEB'  ) /= 0 .or. &
       index(field % field, 'DU4_MAR'  ) /= 0 .or. &
       index(field % field, 'DU4_APR'  ) /= 0 .or. &
       index(field % field, 'DU4_MAY'  ) /= 0 .or. &
       index(field % field, 'DU4_JUN'  ) /= 0 .or. &
       index(field % field, 'DU4_JUL'  ) /= 0 .or. &
       index(field % field, 'DU4_AUG'  ) /= 0 .or. &
       index(field % field, 'DU4_SEP'  ) /= 0 .or. &
       index(field % field, 'DU4_OCT'  ) /= 0 .or. &
       index(field % field, 'DU4_NOV'  ) /= 0 .or. &
       index(field % field, 'DU4_DEC'  ) /= 0 .or. &
       index(field % field, 'DU5_JAN'  ) /= 0 .or. &
       index(field % field, 'DU5_FEB'  ) /= 0 .or. &
       index(field % field, 'DU5_MAR'  ) /= 0 .or. &
       index(field % field, 'DU5_APR'  ) /= 0 .or. &
       index(field % field, 'DU5_MAY'  ) /= 0 .or. &
       index(field % field, 'DU5_JUN'  ) /= 0 .or. &
       index(field % field, 'DU5_JUL'  ) /= 0 .or. &
       index(field % field, 'DU5_AUG'  ) /= 0 .or. &
       index(field % field, 'DU5_SEP'  ) /= 0 .or. &
       index(field % field, 'DU5_OCT'  ) /= 0 .or. &
       index(field % field, 'DU5_NOV'  ) /= 0 .or. &
       index(field % field, 'DU5_DEC'  ) /= 0 .or. &
       index(field % field, 'OC_JAN'   ) /= 0 .or. &
       index(field % field, 'OC_FEB'   ) /= 0 .or. &
       index(field % field, 'OC_MAR'   ) /= 0 .or. &
       index(field % field, 'OC_APR'   ) /= 0 .or. &
       index(field % field, 'OC_MAY'   ) /= 0 .or. &
       index(field % field, 'OC_JUN'   ) /= 0 .or. &
       index(field % field, 'OC_JUL'   ) /= 0 .or. &
       index(field % field, 'OC_AUG'   ) /= 0 .or. &
       index(field % field, 'OC_SEP'   ) /= 0 .or. &
       index(field % field, 'OC_OCT'   ) /= 0 .or. &
       index(field % field, 'OC_NOV'   ) /= 0 .or. &
       index(field % field, 'OC_DEC'   ) /= 0 .or. &
       index(field % field, 'SO4_JAN'  ) /= 0 .or. &
       index(field % field, 'SO4_FEB'  ) /= 0 .or. &
       index(field % field, 'SO4_MAR'  ) /= 0 .or. &
       index(field % field, 'SO4_APR'  ) /= 0 .or. &
       index(field % field, 'SO4_MAY'  ) /= 0 .or. &
       index(field % field, 'SO4_JUN'  ) /= 0 .or. &
       index(field % field, 'SO4_JUL'  ) /= 0 .or. &
       index(field % field, 'SO4_AUG'  ) /= 0 .or. &
       index(field % field, 'SO4_SEP'  ) /= 0 .or. &
       index(field % field, 'SO4_OCT'  ) /= 0 .or. &
       index(field % field, 'SO4_NOV'  ) /= 0 .or. &
       index(field % field, 'SO4_DEC'  ) /= 0 .or. &
       index(field % field, 'SS1_JAN'  ) /= 0 .or. &
       index(field % field, 'SS1_FEB'  ) /= 0 .or. &
       index(field % field, 'SS1_MAR'  ) /= 0 .or. &
       index(field % field, 'SS1_APR'  ) /= 0 .or. &
       index(field % field, 'SS1_MAY'  ) /= 0 .or. &
       index(field % field, 'SS1_JUN'  ) /= 0 .or. &
       index(field % field, 'SS1_JUL'  ) /= 0 .or. &
       index(field % field, 'SS1_AUG'  ) /= 0 .or. &
       index(field % field, 'SS1_SEP'  ) /= 0 .or. &
       index(field % field, 'SS1_OCT'  ) /= 0 .or. &
       index(field % field, 'SS1_NOV'  ) /= 0 .or. &
       index(field % field, 'SS1_DEC'  ) /= 0 .or. &
       index(field % field, 'SS2_JAN'  ) /= 0 .or. &
       index(field % field, 'SS2_FEB'  ) /= 0 .or. &
       index(field % field, 'SS2_MAR'  ) /= 0 .or. &
       index(field % field, 'SS2_APR'  ) /= 0 .or. &
       index(field % field, 'SS2_MAY'  ) /= 0 .or. &
       index(field % field, 'SS2_JUN'  ) /= 0 .or. &
       index(field % field, 'SS2_JUL'  ) /= 0 .or. &
       index(field % field, 'SS2_AUG'  ) /= 0 .or. &
       index(field % field, 'SS2_SEP'  ) /= 0 .or. &
       index(field % field, 'SS2_OCT'  ) /= 0 .or. &
       index(field % field, 'SS2_NOV'  ) /= 0 .or. &
       index(field % field, 'SS2_DEC'  ) /= 0 .or. &
       index(field % field, 'P_WIF_JAN') /= 0 .or. &
       index(field % field, 'P_WIF_FEB') /= 0 .or. &
       index(field % field, 'P_WIF_MAR') /= 0 .or. &
       index(field % field, 'P_WIF_APR') /= 0 .or. &
       index(field % field, 'P_WIF_MAY') /= 0 .or. &
       index(field % field, 'P_WIF_JUN') /= 0 .or. &
       index(field % field, 'P_WIF_JUL') /= 0 .or. &
       index(field % field, 'P_WIF_AUG') /= 0 .or. &
       index(field % field, 'P_WIF_SEP') /= 0 .or. &
       index(field % field, 'P_WIF_OCT') /= 0 .or. &
       index(field % field, 'P_WIF_NOV') /= 0 .or. &
       index(field % field, 'P_WIF_DEC') /= 0) then

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
       if(index(field%field,'QNIFA_JAN') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating QNIFA_JAN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nifa_clim',destField3d)
          nmonths = 1
       elseif(index(field%field,'QNIFA_FEB') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating QNIFA_FEB at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nifa_clim',destField3d)
          nmonths = 2
       elseif(index(field%field,'QNIFA_MAR') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating QNIFA_MAR at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nifa_clim',destField3d)
          nmonths = 3
       elseif(index(field%field,'QNIFA_APR') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating QNIFA_APR at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nifa_clim',destField3d)
          nmonths = 4
       elseif(index(field%field,'QNIFA_MAY') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating QNIFA_MAY at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nifa_clim',destField3d)
          nmonths = 5
       elseif(index(field%field,'QNIFA_JUN') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating QNIFA_JUN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nifa_clim',destField3d)
          nmonths = 6
       elseif(index(field%field,'QNIFA_JUL') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating QNIFA_JUL at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nifa_clim',destField3d)
          nmonths = 7
       elseif(index(field%field,'QNIFA_AUG') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating QNIFA_AUG at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nifa_clim',destField3d)
          nmonths = 8
       elseif(index(field%field,'QNIFA_SEP') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating QNIFA_SEP at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nifa_clim',destField3d)
          nmonths = 9
       elseif(index(field%field,'QNIFA_OCT') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating QNIFA_OCT at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nifa_clim',destField3d)
          nmonths = 10
       elseif(index(field%field,'QNIFA_NOV') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating QNIFA_NOV at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nifa_clim',destField3d)
          nmonths = 11
       elseif(index(field%field,'QNIFA_DEC') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating QNIFA_DEC at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nifa_clim',destField3d)
          nmonths = 12
       elseif(index(field%field,'QNWFA_JAN') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating QNWFA_JAN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nwfa_clim',destField3d)
          nmonths = 1
       elseif(index(field%field,'QNWFA_FEB') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating QNWFA_FEB at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nwfa_clim',destField3d)
          nmonths = 2
       elseif(index(field%field,'QNWFA_MAR') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating QNWFA_MAR at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nwfa_clim',destField3d)
          nmonths = 3
       elseif(index(field%field,'QNWFA_APR') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating QNWFA_APR at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nwfa_clim',destField3d)
          nmonths = 4
       elseif(index(field%field,'QNWFA_MAY') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating QNWFA_MAY at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nwfa_clim',destField3d)
          nmonths = 5
       elseif(index(field%field,'QNWFA_JUN') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating QNWFA_JUN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nwfa_clim',destField3d)
          nmonths = 6
       elseif(index(field%field,'QNWFA_JUL') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating QNWFA_JUL at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nwfa_clim',destField3d)
          nmonths = 7
       elseif(index(field%field,'QNWFA_AUG') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating QNWFA_AUG at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nwfa_clim',destField3d)
          nmonths = 8
       elseif(index(field%field,'QNWFA_SEP') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating QNWFA_SEP at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nwfa_clim',destField3d)
          nmonths = 9
       elseif(index(field%field,'QNWFA_OCT') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating QNWFA_OCT at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nwfa_clim',destField3d)
          nmonths = 10
       elseif(index(field%field,'QNWFA_NOV') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating QNWFA_NOV at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nwfa_clim',destField3d)
          nmonths = 11
       elseif(index(field%field,'QNWFA_DEC') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating QNWFA_DEC at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'nwfa_clim',destField3d)
          nmonths = 12
       elseif(index(field%field,'DU1_JAN') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU1_JAN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust1_clim
          nmonths = 1
       elseif(index(field%field,'DU1_FEB') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU1_FEB at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust1_clim
          nmonths = 2
       elseif(index(field%field,'DU1_MAR') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU1_MAR at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust1_clim
          nmonths = 3
       elseif(index(field%field,'DU1_APR') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU1_APR at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust1_clim
          nmonths = 4
       elseif(index(field%field,'DU1_MAY') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU1_MAY at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust1_clim
          nmonths = 5
       elseif(index(field%field,'DU1_JUN') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU1_JUN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust1_clim
          nmonths = 6
       elseif(index(field%field,'DU1_JUL') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU1_JUL at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust1_clim
          nmonths = 7
       elseif(index(field%field,'DU1_AUG') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU1_AUG at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust1_clim
          nmonths = 8
       elseif(index(field%field,'DU1_SEP') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU1_SEP at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust1_clim
          nmonths = 9
       elseif(index(field%field,'DU1_OCT') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU1_OCT at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust1_clim
          nmonths = 10
       elseif(index(field%field,'DU1_NOV') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU1_NOV at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust1_clim
          nmonths = 11
       elseif(index(field%field,'DU1_DEC') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU1_DEC at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust1_clim
          nmonths = 12
       elseif(index(field%field,'DU2_JAN') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU2_JAN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust2_clim
          nmonths = 1
       elseif(index(field%field,'DU2_FEB') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU2_FEB at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust2_clim
          nmonths = 2
       elseif(index(field%field,'DU2_MAR') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU2_MAR at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust2_clim
          nmonths = 3
       elseif(index(field%field,'DU2_APR') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU2_APR at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust2_clim
          nmonths = 4
       elseif(index(field%field,'DU2_MAY') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU2_MAY at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust2_clim
          nmonths = 5
       elseif(index(field%field,'DU2_JUN') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU2_JUN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust2_clim
          nmonths = 6
       elseif(index(field%field,'DU2_JUL') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU2_JUL at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust2_clim
          nmonths = 7
       elseif(index(field%field,'DU2_AUG') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU2_AUG at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust2_clim
          nmonths = 8
       elseif(index(field%field,'DU2_SEP') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU2_SEP at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust2_clim
          nmonths = 9
       elseif(index(field%field,'DU2_OCT') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU2_OCT at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust2_clim
          nmonths = 10
       elseif(index(field%field,'DU2_NOV') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU2_NOV at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust2_clim
          nmonths = 11
       elseif(index(field%field,'DU2_DEC') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU2_DEC at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust2_clim
          nmonths = 12
       elseif(index(field%field,'DU3_JAN') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU3_JAN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust3_clim
          nmonths = 1
       elseif(index(field%field,'DU3_FEB') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU3_FEB at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust3_clim
          nmonths = 2
       elseif(index(field%field,'DU3_MAR') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU3_MAR at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust3_clim
          nmonths = 3
       elseif(index(field%field,'DU3_APR') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU3_APR at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust3_clim
          nmonths = 4
       elseif(index(field%field,'DU3_MAY') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU3_MAY at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust3_clim
          nmonths = 5
       elseif(index(field%field,'DU3_JUN') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU3_JUN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust3_clim
          nmonths = 6
       elseif(index(field%field,'DU3_JUL') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU3_JUL at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust3_clim
          nmonths = 7
       elseif(index(field%field,'DU3_AUG') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU3_AUG at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust3_clim
          nmonths = 8
       elseif(index(field%field,'DU3_SEP') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU3_SEP at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust3_clim
          nmonths = 9
       elseif(index(field%field,'DU3_OCT') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU3_OCT at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust3_clim
          nmonths = 10
       elseif(index(field%field,'DU3_NOV') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU3_NOV at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust3_clim
          nmonths = 11
       elseif(index(field%field,'DU3_DEC') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU3_DEC at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust3_clim
          nmonths = 12
        elseif(index(field%field,'DU4_JAN') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU4_JAN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust4_clim
          nmonths = 1
       elseif(index(field%field,'DU4_FEB') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU4_FEB at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust4_clim
          nmonths = 2
       elseif(index(field%field,'DU4_MAR') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU4_MAR at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust4_clim
          nmonths = 3
       elseif(index(field%field,'DU4_APR') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU4_APR at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust4_clim
          nmonths = 4
       elseif(index(field%field,'DU4_MAY') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU4_MAY at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust4_clim
          nmonths = 5
       elseif(index(field%field,'DU4_JUN') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU4_JUN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust4_clim
          nmonths = 6
       elseif(index(field%field,'DU4_JUL') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU4_JUL at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust4_clim
          nmonths = 7
       elseif(index(field%field,'DU4_AUG') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU4_AUG at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust4_clim
          nmonths = 8
       elseif(index(field%field,'DU4_SEP') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU4_SEP at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust4_clim
          nmonths = 9
       elseif(index(field%field,'DU4_OCT') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU4_OCT at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust4_clim
          nmonths = 10
       elseif(index(field%field,'DU4_NOV') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU4_NOV at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust4_clim
          nmonths = 11
       elseif(index(field%field,'DU4_DEC') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU4_DEC at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust4_clim
          nmonths = 12
        elseif(index(field%field,'DU5_JAN') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU5_JAN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust5_clim
          nmonths = 1
       elseif(index(field%field,'DU5_FEB') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU5_FEB at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust5_clim
          nmonths = 2
       elseif(index(field%field,'DU5_MAR') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU5_MAR at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust5_clim
          nmonths = 3
       elseif(index(field%field,'DU5_APR') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU5_APR at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust5_clim
          nmonths = 4
       elseif(index(field%field,'DU5_MAY') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU5_MAY at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust5_clim
          nmonths = 5
       elseif(index(field%field,'DU5_JUN') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU5_JUN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust5_clim
          nmonths = 6
       elseif(index(field%field,'DU5_JUL') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU5_JUL at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust5_clim
          nmonths = 7
       elseif(index(field%field,'DU5_AUG') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU5_AUG at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust5_clim
          nmonths = 8
       elseif(index(field%field,'DU5_SEP') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU5_SEP at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust5_clim
          nmonths = 9
       elseif(index(field%field,'DU5_OCT') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU5_OCT at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust5_clim
          nmonths = 10
       elseif(index(field%field,'DU5_NOV') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU5_NOV at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust5_clim
          nmonths = 11
       elseif(index(field%field,'DU5_DEC') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating DU5_DEC at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => dust5_clim
          nmonths = 12
       elseif(index(field%field,'OC_JAN') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating OC_JAN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => oc_clim
          nmonths = 1
       elseif(index(field%field,'OC_FEB') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating OC_FEB at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => oc_clim
          nmonths = 2
       elseif(index(field%field,'OC_MAR') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating OC_MAR at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => oc_clim
          nmonths = 3
       elseif(index(field%field,'OC_APR') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating OC_APR at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => oc_clim
          nmonths = 4
       elseif(index(field%field,'OC_MAY') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating OC_MAY at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => oc_clim
          nmonths = 5
       elseif(index(field%field,'OC_JUN') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating OC_JUN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => oc_clim
          nmonths = 6
       elseif(index(field%field,'OC_JUL') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating OC_JUL at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => oc_clim
          nmonths = 7
       elseif(index(field%field,'OC_AUG') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating OC_AUG at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => oc_clim
          nmonths = 8
       elseif(index(field%field,'OC_SEP') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating OC_SEP at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => oc_clim
          nmonths = 9
       elseif(index(field%field,'OC_OCT') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating OC_OCT at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => oc_clim
          nmonths = 10
       elseif(index(field%field,'OC_NOV') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating OC_NOV at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => oc_clim
          nmonths = 11
       elseif(index(field%field,'OC_DEC') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating OC_DEC at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => oc_clim
          nmonths = 12
       elseif(index(field%field,'SO4_JAN') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SO4_JAN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => so4_clim
          nmonths = 1
       elseif(index(field%field,'SO4_FEB') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SO4_FEB at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => so4_clim
          nmonths = 2
       elseif(index(field%field,'SO4_MAR') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SO4_MAR at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => so4_clim
          nmonths = 3
       elseif(index(field%field,'SO4_APR') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SO4_APR at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => so4_clim
          nmonths = 4
       elseif(index(field%field,'SO4_MAY') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SO4_MAY at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => so4_clim
          nmonths = 5
       elseif(index(field%field,'SO4_JUN') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SO4_JUN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => so4_clim
          nmonths = 6
       elseif(index(field%field,'SO4_JUL') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SO4_JUL at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => so4_clim
          nmonths = 7
       elseif(index(field%field,'SO4_AUG') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SO4_AUG at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => so4_clim
          nmonths = 8
       elseif(index(field%field,'SO4_SEP') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SO4_SEP at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => so4_clim
          nmonths = 9
       elseif(index(field%field,'SO4_OCT') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SO4_OCT at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => so4_clim
          nmonths = 10
       elseif(index(field%field,'SO4_NOV') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SO4_NOV at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => so4_clim
          nmonths = 11
       elseif(index(field%field,'SO4_DEC') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SO4_DEC at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => so4_clim
          nmonths = 12
        elseif(index(field%field,'SS1_JAN') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SS1_JAN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => seas1_clim
          nmonths = 1
       elseif(index(field%field,'SS1_FEB') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SS1_FEB at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => seas1_clim
          nmonths = 2
       elseif(index(field%field,'SS1_MAR') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SS1_MAR at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => seas1_clim
          nmonths = 3
       elseif(index(field%field,'SS1_APR') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SS1_APR at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => seas1_clim
          nmonths = 4
       elseif(index(field%field,'SS1_MAY') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SS1_MAY at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => seas1_clim
          nmonths = 5
       elseif(index(field%field,'SS1_JUN') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SS1_JUN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => seas1_clim
          nmonths = 6
       elseif(index(field%field,'SS1_JUL') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SS1_JUL at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => seas1_clim
          nmonths = 7
       elseif(index(field%field,'SS1_AUG') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SS1_AUG at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => seas1_clim
          nmonths = 8
       elseif(index(field%field,'SS1_SEP') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SS1_SEP at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => seas1_clim
          nmonths = 9
       elseif(index(field%field,'SS1_OCT') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SS1_OCT at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => seas1_clim
          nmonths = 10
       elseif(index(field%field,'SS1_NOV') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SS1_NOV at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => seas1_clim
          nmonths = 11
       elseif(index(field%field,'SS1_DEC') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SS1_DEC at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => seas1_clim
          nmonths = 12
       elseif(index(field%field,'SS2_JAN') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SS2_JAN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => seas2_clim
          nmonths = 1
       elseif(index(field%field,'SS2_FEB') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SS2_FEB at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => seas2_clim
          nmonths = 2
       elseif(index(field%field,'SS2_MAR') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SS2_MAR at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => seas2_clim
          nmonths = 3
       elseif(index(field%field,'SS2_APR') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SS2_APR at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => seas2_clim
          nmonths = 4
       elseif(index(field%field,'SS2_MAY') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SS2_MAY at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => seas2_clim
          nmonths = 5
       elseif(index(field%field,'SS2_JUN') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SS2_JUN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => seas2_clim
          nmonths = 6
       elseif(index(field%field,'SS2_JUL') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SS2_JUL at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => seas2_clim
          nmonths = 7
       elseif(index(field%field,'SS2_AUG') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SS2_AUG at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => seas2_clim
          nmonths = 8
       elseif(index(field%field,'SS2_SEP') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SS2_SEP at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => seas2_clim
          nmonths = 9
       elseif(index(field%field,'SS2_OCT') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SS2_OCT at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => seas2_clim
          nmonths = 10
       elseif(index(field%field,'SS2_NOV') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SS2_NOV at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => seas2_clim
          nmonths = 11
       elseif(index(field%field,'SS2_DEC') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating SS2_DEC at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField3d => seas2_clim
          nmonths = 12
      elseif(index(field%field,'P_WIF_JAN') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating P_WIF_JAN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'pwif_clim',destField3d)
          nmonths = 1
       elseif(index(field%field,'P_WIF_FEB') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating P_WIF_FEB at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'pwif_clim',destField3d)
          nmonths = 2
       elseif(index(field%field,'P_WIF_MAR') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating P_WIF_MAR at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'pwif_clim',destField3d)
          nmonths = 3
       elseif(index(field%field,'P_WIF_APR') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating P_WIF_APR at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'pwif_clim',destField3d)
          nmonths = 4
       elseif(index(field%field,'P_WIF_MAY') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating P_WIF_MAY at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'pwif_clim',destField3d)
          nmonths = 5
       elseif(index(field%field,'P_WIF_JUN') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating P_WIF_JUN at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'pwif_clim',destField3d)
          nmonths = 6
       elseif(index(field%field,'P_WIF_JUL') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating P_WIF_JUL at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'pwif_clim',destField3d)
          nmonths = 7
       elseif(index(field%field,'P_WIF_AUG') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating P_WIF_AUG at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'pwif_clim',destField3d)
          nmonths = 8
       elseif(index(field%field,'P_WIF_SEP') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating P_WIF_SEP at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'pwif_clim',destField3d)
          nmonths = 9
       elseif(index(field%field,'P_WIF_OCT') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating P_WIF_OCT at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'pwif_clim',destField3d)
          nmonths = 10
       elseif(index(field%field,'P_WIF_NOV') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating P_WIF_NOV at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'pwif_clim',destField3d)
          nmonths = 11
       elseif(index(field%field,'P_WIF_DEC') /= 0) then
          k = field%xlvl
          call mpas_log_write('Interpolating P_WIF_DEC at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          call mpas_pool_get_array(mesh,'pwif_clim',destField3d)
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
                destField3d(nmonths,k,iCell) = interp_sequence(x,y,1,rslab,-2,field%nx+3,1,field%ny,1,1,msgval, \
                                                  interp_list,1,maskval=maskval,mask_array=maskslab)
             else
                destField3d(nmonths,k,iCell) = interp_sequence(x,y,1,rslab,-2,field%nx+3,1,field%ny,1,1,msgval, \
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

 call mpas_log_write('--- end subroutine init_hinterp_gocart:')

 end subroutine init_hinterp_gocart
 
!=================================================================================================================
 subroutine read_gocart_init(fg_source,istatus)
!=================================================================================================================

!input arguments:
 character(len=*),intent(in):: fg_source

!output arguments:
 integer,intent(out):: istatus

!local variables:
 logical:: is_used
 integer:: io_status

!-----------------------------------------------------------------------------------------------------------------

 istatus = 0

 do input_unit = 10, 100
    inquire(unit=input_unit,opened=is_used)
    if (.not. is_used) exit
 end do
 if(input_unit > 100) &
    call mpas_log_write('Error: In read_gocart_init(), couldn''t find an available Fortran unit.')
 open(unit=input_unit, file=trim(fg_source), status='old', form='unformatted', iostat=io_status)

 if (io_status > 0) istatus = 1

 end subroutine read_gocart_init

!=================================================================================================================
 end module mpas_init_atm_gocart
!=================================================================================================================
