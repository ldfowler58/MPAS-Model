! Copyright (c) 2024 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================================
 module mpas_chemistry_init_gocart2G_interp
 use mpas_dmpar
 use mpas_log
 use mpas_kind_types
 use mpas_pool_routines
 use mpas_stream_manager
 use mpas_timekeeping,only: sub_t_t
 use mpas_derived_types,only     : mpas_time_type,mpas_timeinterval_type
 use mpas_timekeeping,only       : mpas_get_clock_time,mpas_set_time,mpas_get_timeInterval
 use init_atm_read_met,only      : met_data,read_met_init,read_met_close,read_next_met_field
 use init_atm_hinterp,only       : interp_sequence,FOUR_POINT,SEARCH,W_AVERAGE4,W_AVERAGE16
 use init_atm_llxy,only          : latlon_to_ij,map_init,map_set,proj_info,PROJ_LATLON,PROJ_GAUSS,DEG_PER_RAD


 implicit none
 private
 public:: init_gocart2G_aerosols


 contains


!=================================================================================================================
 subroutine init_gocart2G_aerosols(configs,mesh,fg,diag,state,gocart2G_backgrounds)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: diag

!inout arguments:
 type(mpas_pool_type),intent(inout):: fg
 type(mpas_pool_type),intent(inout):: state
 type(mpas_pool_type),intent(inout):: gocart2G_backgrounds

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine init_gocart2G_aerosols:')

 call init_hinterp_gocart2G(configs,mesh,fg)
 call init_vinterp_gocart2G(configs,mesh,fg,diag,state)
 call init_vinterp_gocart2G_hno3(configs,mesh,fg,diag,gocart2G_backgrounds)

 call mpas_log_write('--- end subroutine init_gocart2G_aerosols.')
 call mpas_log_write(' ')

 end subroutine init_gocart2G_aerosols

!==================================================================================================================
 subroutine init_hinterp_gocart2G(configs,mesh,fg)
!==================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh

!inout arguments:
 type(mpas_pool_type),intent(inout):: fg

!local variables and arrays:
 type(met_data) :: field !real*4 meteorological data.
 type(proj_info):: proj

 character(len=StrKIND),pointer:: prefix,start_time

 logical:: have_landmask

!First-guess gocart2G aerosol mixing ratios:
 integer:: iCell,istatus,k,masked,nInterpPoints
 integer,dimension(5):: interp_list

 integer,pointer:: index_qbcphobic,index_qbcphilic
 integer,pointer:: index_qbrphobic,index_qbrphilic
 integer,pointer:: index_qocphobic,index_qocphilic
 integer,pointer:: index_qdust1,index_qdust2,index_qdust3,index_qdust4,index_qdust5
 integer,pointer:: index_qni1,index_qni2,index_qni3
 integer,pointer:: index_qseas1,index_qseas2,index_qseas3,index_qseas4,index_qseas5
 integer,pointer:: index_qso2,index_qso2v,index_qso4,index_qso4v
 integer,pointer:: index_qdms,index_qmsa
 integer,pointer:: index_qnh3,index_qnh4a,index_qhno3

 integer,pointer:: nCells,nAerLevels
 integer,pointer:: num_scalars_fg
 integer,pointer:: gocart2G_start,gocart2G_end
 integer,dimension(:),pointer:: landmask,mask_array

 real(kind=RKIND):: fillval,maskval,msgval
 real(kind=RKIND):: lat,lon,x,y
 real(kind=RKIND),dimension(:),pointer:: latCell,lonCell
 real(kind=RKIND),dimension(:),pointer:: latPoints,lonPoints

 real(kind=RKIND),dimension(:,:),pointer:: dpgoc,pgoc

 real(kind=RKIND),dimension(:,:),pointer:: qbcphobic,qbcphilic
 real(kind=RKIND),dimension(:,:),pointer:: qbrphobic,qbrphilic
 real(kind=RKIND),dimension(:,:),pointer:: qocphobic,qocphilic
 real(kind=RKIND),dimension(:,:),pointer:: qdust1,qdust2,qdust3,qdust4,qdust5
 real(kind=RKIND),dimension(:,:),pointer:: qni1,qni2,qni3
 real(kind=RKIND),dimension(:,:),pointer:: qseas1,qseas2,qseas3,qseas4,qseas5
 real(kind=RKIND),dimension(:,:),pointer:: qso2,qso2v,qso4,qso4v
 real(kind=RKIND),dimension(:,:),pointer:: qdms,qmsa
 real(kind=RKIND),dimension(:,:),pointer:: qnh3,qnh4a,qhno3
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars_fg

 real(kind=RKIND),dimension(:,:),pointer:: destField2d
 real(kind=RKIND),dimension(:,:),allocatable:: maskslab,rslab

 real(kind=RKIND):: rmax

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine init_hinterp_gocart2G:')

 call mpas_pool_get_config(configs,'config_aerosolFG_prefix',prefix)
 call mpas_pool_get_config(configs,'config_start_time',start_time)

 call mpas_pool_get_dimension(mesh,'nCells',nCells)
 call mpas_pool_get_dimension(mesh,'nAerLevels',nAerLevels)

 call mpas_pool_get_array(mesh,'landmask',landmask)
 call mpas_pool_get_array(mesh,'latCell' ,latCell )
 call mpas_pool_get_array(mesh,'lonCell' ,lonCell )


!--- list of input gocart2G aerosol species from intermediate binary file:
 call mpas_pool_get_dimension(fg,'index_qbcphobic',index_qbcphobic)
 call mpas_pool_get_dimension(fg,'index_qbcphilic',index_qbcphilic)
 call mpas_pool_get_dimension(fg,'index_qbrphobic',index_qbrphobic)
 call mpas_pool_get_dimension(fg,'index_qbrphilic',index_qbrphilic)
 call mpas_pool_get_dimension(fg,'index_qocphobic',index_qocphobic)
 call mpas_pool_get_dimension(fg,'index_qocphilic',index_qocphilic)
 call mpas_pool_get_dimension(fg,'index_qdust1'   ,index_qdust1   )
 call mpas_pool_get_dimension(fg,'index_qdust2'   ,index_qdust2   )
 call mpas_pool_get_dimension(fg,'index_qdust3'   ,index_qdust3   )
 call mpas_pool_get_dimension(fg,'index_qdust4'   ,index_qdust4   )
 call mpas_pool_get_dimension(fg,'index_qdust5'   ,index_qdust5   )
 call mpas_pool_get_dimension(fg,'index_qni1'     ,index_qni1     )
 call mpas_pool_get_dimension(fg,'index_qni2'     ,index_qni2     )
 call mpas_pool_get_dimension(fg,'index_qni3'     ,index_qni3     )
 call mpas_pool_get_dimension(fg,'index_qso2'     ,index_qso2     )
 call mpas_pool_get_dimension(fg,'index_qso2v'    ,index_qso2v    )
 call mpas_pool_get_dimension(fg,'index_qso4'     ,index_qso4     )
 call mpas_pool_get_dimension(fg,'index_qso4v'    ,index_qso4v    )
 call mpas_pool_get_dimension(fg,'index_qseas1'   ,index_qseas1   )
 call mpas_pool_get_dimension(fg,'index_qseas2'   ,index_qseas2   )
 call mpas_pool_get_dimension(fg,'index_qseas3'   ,index_qseas3   )
 call mpas_pool_get_dimension(fg,'index_qseas4'   ,index_qseas4   )
 call mpas_pool_get_dimension(fg,'index_qseas5'   ,index_qseas5   )
 call mpas_pool_get_dimension(fg,'index_qdms'     ,index_qdms     )
 call mpas_pool_get_dimension(fg,'index_qmsa'     ,index_qmsa     )
 call mpas_pool_get_dimension(fg,'index_qnh3'     ,index_qnh3     )
 call mpas_pool_get_dimension(fg,'index_qnh4a'    ,index_qnh4a    )
 call mpas_pool_get_dimension(fg,'index_qhno3'    ,index_qhno3    )
 call mpas_pool_get_array(fg,'scalars_fg',scalars_fg)
 scalars_fg = 0._RKIND

 call mpas_log_write('--- index_qbcphobic = $i',intArgs=(/index_qbcphobic/))
 call mpas_log_write('--- index_qbcphilic = $i',intArgs=(/index_qbcphilic/))
 call mpas_log_write('--- index_qbrphobic = $i',intArgs=(/index_qbrphobic/))
 call mpas_log_write('--- index_qbrphilic = $i',intArgs=(/index_qbrphilic/))
 call mpas_log_write('--- index_qocphobic = $i',intArgs=(/index_qocphobic/))
 call mpas_log_write('--- index_qocphilic = $i',intArgs=(/index_qocphilic/))
 call mpas_log_write('--- index_qdust1    = $i',intArgs=(/index_qdust1/))
 call mpas_log_write('--- index_qdust2    = $i',intArgs=(/index_qdust2/))
 call mpas_log_write('--- index_qdust3    = $i',intArgs=(/index_qdust3/))
 call mpas_log_write('--- index_qdust4    = $i',intArgs=(/index_qdust4/))
 call mpas_log_write('--- index_qdust5    = $i',intArgs=(/index_qdust5/))
 call mpas_log_write('--- index_qni1      = $i',intArgs=(/index_qni1/))
 call mpas_log_write('--- index_qni2      = $i',intArgs=(/index_qni2/))
 call mpas_log_write('--- index_qni3      = $i',intArgs=(/index_qni3/))
 call mpas_log_write('--- index_qso2      = $i',intArgs=(/index_qso2/))
 call mpas_log_write('--- index_qso2v     = $i',intArgs=(/index_qso2v/))
 call mpas_log_write('--- index_qso4      = $i',intArgs=(/index_qso4/))
 call mpas_log_write('--- index_qso4v     = $i',intArgs=(/index_qso4v/))
 call mpas_log_write('--- index_qseas1    = $i',intArgs=(/index_qseas1/))
 call mpas_log_write('--- index_qseas2    = $i',intArgs=(/index_qseas2/))
 call mpas_log_write('--- index_qseas3    = $i',intArgs=(/index_qseas3/))
 call mpas_log_write('--- index_qseas4    = $i',intArgs=(/index_qseas4/))
 call mpas_log_write('--- index_qseas5    = $i',intArgs=(/index_qseas5/))
 call mpas_log_write('--- index_qdms      = $i',intArgs=(/index_qdms/))
 call mpas_log_write('--- index_qmsa      = $i',intArgs=(/index_qmsa/))
 call mpas_log_write('--- index_qnh3      = $i',intArgs=(/index_qnh3/))
 call mpas_log_write('--- index_qnh4a     = $i',intArgs=(/index_qnh4a/))
 call mpas_log_write('--- index_qhno3     = $i',intArgs=(/index_qhno3/))

 qbcphobic => scalars_fg(index_qbcphobic,:,:)
 qbcphilic => scalars_fg(index_qbcphilic,:,:)
 qbrphobic => scalars_fg(index_qbrphobic,:,:)
 qbrphilic => scalars_fg(index_qbrphilic,:,:)
 qocphobic => scalars_fg(index_qocphobic,:,:)
 qocphilic => scalars_fg(index_qocphilic,:,:)
 qdust1    => scalars_fg(index_qdust1,:,:)
 qdust2    => scalars_fg(index_qdust2,:,:)
 qdust3    => scalars_fg(index_qdust3,:,:)
 qdust4    => scalars_fg(index_qdust4,:,:)
 qdust5    => scalars_fg(index_qdust5,:,:)
 qni1      => scalars_fg(index_qni1,:,:)
 qni2      => scalars_fg(index_qni2,:,:)
 qni3      => scalars_fg(index_qni3,:,:)
 qso2      => scalars_fg(index_qso2,:,:)
 qso2v     => scalars_fg(index_qso2v,:,:)
 qso4      => scalars_fg(index_qso4,:,:)
 qso4v     => scalars_fg(index_qso4v,:,:)
 qseas1    => scalars_fg(index_qseas1,:,:)
 qseas2    => scalars_fg(index_qseas2,:,:)
 qseas3    => scalars_fg(index_qseas3,:,:)
 qseas4    => scalars_fg(index_qseas4,:,:)
 qseas5    => scalars_fg(index_qseas5,:,:)
 qdms      => scalars_fg(index_qdms,:,:)
 qmsa      => scalars_fg(index_qmsa,:,:)
 qnh3      => scalars_fg(index_qnh3,:,:)
 qnh4a     => scalars_fg(index_qnh4a,:,:)
 qhno3     => scalars_fg(index_qhno3,:,:)

 call mpas_pool_get_dimension(fg,'num_scalars_fg' ,num_scalars_fg )
 call mpas_pool_get_dimension(fg,'gocart2G_start',gocart2G_start)
 call mpas_pool_get_dimension(fg,'gocart2G_end'  ,gocart2G_end  )
 call mpas_log_write('--- num_scalars_fg     = $i',intArgs=(/num_scalars_fg/) )
 call mpas_log_write('--- gocart2G_start     = $i',intArgs=(/gocart2G_start/))
 call mpas_log_write('--- gocart2G_end       = $i',intArgs=(/gocart2G_end/)  )

 call mpas_pool_get_array(fg,'pgoc',pgoc)
 call mpas_pool_get_array(fg,'dpgoc',dpgoc)


!--- open intermediate binary file:
 istatus = 0
 call read_met_init(trim(prefix),.false.,start_time(1:13),istatus)
 if(istatus /= 0) then
    call mpas_log_write('**************************************************',messageType=MPAS_LOG_ERR)
    call mpas_log_write('Error opening intermediate input data file ' &
                                       //trim(prefix)//':'//start_time(1:13),messageType=MPAS_LOG_ERR)
    call mpas_log_write('**************************************************',messageType=MPAS_LOG_CRIT)
 endif


!scan through all fields in the file, looking for the LANDSEA field:
 have_landmask = .false.
 call read_next_met_field(field,istatus)
 do while (istatus == 0)
    if(index(field%field, 'LANDSEA') /= 0) then
       have_landmask = .true.
       if(.not.allocated(maskslab)) allocate(maskslab(-2:field%nx+3,field%ny))

       maskslab(1:field%nx,1:field%ny) = field%slab(1:field%nx,1:field%ny)
       maskslab(0 ,1:field%ny) = field%slab(field%nx  ,1:field%ny)
       maskslab(-1,1:field%ny) = field%slab(field%nx-1,1:field%ny)
       maskslab(-2,1:field%ny) = field%slab(field%nx-2,1:field%ny)
       maskslab(field%nx+1,1:field%ny) = field%slab(1,1:field%ny)
       maskslab(field%nx+2,1:field%ny) = field%slab(2,1:field%ny)
       maskslab(field%nx+3,1:field%ny) = field%slab(3,1:field%ny)
       call mpas_log_write('minval,maxval LANDSEA = $r $r',realArgs=(/minval(maskslab),maxval(maskslab)/))
    endif
    !note that field%slab is initialized in subroutine read_next_met_field but deallocated here:
    deallocate(field%slab)
    call read_next_met_field(field,istatus)
 enddo
 call read_met_close()

 if(.not. have_landmask) then
    call mpas_log_write('**************************************************')
    call mpas_log_write('Landsea mask not available from the intermediate CAM-Chem data file ' &
                                       //trim(prefix)//':'//start_time(1:13))
    call mpas_log_write('**************************************************')
    call mpas_log_write(' ')
 endif


!horizontally interpolate first-guess data:
 istatus = 0
 call read_met_init(trim(prefix),.false.,start_time(1:13),istatus)
 if(istatus /= 0) then
    call mpas_log_write('**************************************************',messageType=MPAS_LOG_ERR)
    call mpas_log_write('Error opening intermediate input data file ' &
                                       //trim(prefix)//':'//start_time(1:13),messageType=MPAS_LOG_ERR)
    call mpas_log_write('**************************************************',messageType=MPAS_LOG_CRIT)
 endif
 call read_next_met_field(field,istatus)


 do while(istatus == 0)

!--- use the same values as the default values for meteorological fields in mpas_init_case_gfs:
!   interp_list(1) = FOUR_POINT
!   interp_list(2) = W_AVERAGE4
!   interp_list(3) = W_AVERAGE16
!   interp_list(4) = SEARCH
!   interp_list(5) = 0
    interp_list(1) = FOUR_POINT
    interp_list(2) = SEARCH
    interp_list(3) = 0
!--- end use.

    maskval = -1.0
    masked  = -1
    fillval = 0.0
    msgval  = 1.e15

    mask_array => landmask

    if(trim(field%field) == 'BCPHOBIC' .or. &
       trim(field%field) == 'BCPHILIC' .or. &
       trim(field%field) == 'BRPHOBIC' .or. &
       trim(field%field) == 'BRPHILIC' .or. &
       trim(field%field) == 'OCPHOBIC' .or. &
       trim(field%field) == 'OCPHILIC' .or. &
       trim(field%field) == 'DU001'    .or. &
       trim(field%field) == 'DU002'    .or. &
       trim(field%field) == 'DU003'    .or. &
       trim(field%field) == 'DU004'    .or. &
       trim(field%field) == 'DU005'    .or. &
       trim(field%field) == 'NI001'    .or. &
       trim(field%field) == 'NI002'    .or. &
       trim(field%field) == 'NI003'    .or. &
       trim(field%field) == 'SO2'      .or. &
       trim(field%field) == 'SO2V'     .or. &
       trim(field%field) == 'SO4'      .or. &
       trim(field%field) == 'SO4V'     .or. &
       trim(field%field) == 'SS001'    .or. &
       trim(field%field) == 'SS002'    .or. &
       trim(field%field) == 'SS003'    .or. &
       trim(field%field) == 'SS004'    .or. &
       trim(field%field) == 'SS005'    .or. &
       trim(field%field) == 'DMS'      .or. &
       trim(field%field) == 'MSA'      .or. &
       trim(field%field) == 'HNO3'     .or. &
       trim(field%field) == 'NH3'      .or. &
       trim(field%field) == 'AIRDENS'  .or. &
       trim(field%field) == 'RH'       .or. &
       trim(field%field) == 'DPRES'    .or. &
       trim(field%field) == 'PRES'    ) then

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
       if(trim(field%field) == 'BCPHOBIC') then
          k = field%xlvl
          call mpas_log_write('Interpolating BCPHOBIC at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qbcphobic
       elseif(trim(field%field) == 'BCPHILIC') then
          k = field%xlvl
          call mpas_log_write('Interpolating BCPHILIC at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qbcphilic
       elseif(trim(field%field) == 'BRPHOBIC') then
          k = field%xlvl
          call mpas_log_write('Interpolating BRPHOBIC at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qbrphobic
       elseif(trim(field%field) == 'BRPHILIC') then
          k = field%xlvl
          call mpas_log_write('Interpolating BRPHILIC at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qbrphilic
       elseif(trim(field%field) == 'OCPHOBIC') then
          k = field%xlvl
          call mpas_log_write('Interpolating OCPHOBIC at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qocphobic
       elseif(trim(field%field) == 'OCPHILIC') then
          k = field%xlvl
          call mpas_log_write('Interpolating OCPHILIC at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qocphilic
       elseif(trim(field%field) == 'DU001') then
          k = field%xlvl
          call mpas_log_write('Interpolating DU001 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qdust1
       elseif(trim(field%field) == 'DU002') then
          k = field%xlvl
          call mpas_log_write('Interpolating DU002 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qdust2
       elseif(trim(field%field) == 'DU003') then
          k = field%xlvl
          call mpas_log_write('Interpolating DU003 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qdust3
       elseif(trim(field%field) == 'DU004') then
          k = field%xlvl
          call mpas_log_write('Interpolating DU004 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qdust4
       elseif(trim(field%field) == 'DU005') then
          k = field%xlvl
          call mpas_log_write('Interpolating DU005 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qdust5
       elseif(trim(field%field) == 'NI001') then
          k = field%xlvl
          call mpas_log_write('Interpolating NI001 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qni1
       elseif(trim(field%field) == 'NI002') then
          k = field%xlvl
          call mpas_log_write('Interpolating NI002 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qni2
       elseif(trim(field%field) == 'NI003') then
          k = field%xlvl
          call mpas_log_write('Interpolating NI003 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qni3
       elseif(trim(field%field) == 'SO2') then
          k = field%xlvl
          call mpas_log_write('Interpolating SO2 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qso2
       elseif(trim(field%field) == 'SO2V') then
          k = field%xlvl
          call mpas_log_write('Interpolating SO2V at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qso2v
       elseif(trim(field%field) == 'SO4') then
          k = field%xlvl
          call mpas_log_write('Interpolating SO4 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qso4
       elseif(trim(field%field) == 'SO4V') then
          k = field%xlvl
          call mpas_log_write('Interpolating SO4V at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qso4v
       elseif(trim(field%field) == 'SS001') then
          k = field%xlvl
          call mpas_log_write('Interpolating SS001 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qseas1
       elseif(trim(field%field) == 'SS002') then
          k = field%xlvl
          call mpas_log_write('Interpolating SS002 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qseas2
       elseif(trim(field%field) == 'SS003') then
          k = field%xlvl
          call mpas_log_write('Interpolating SS003 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qseas3
       elseif(trim(field%field) == 'SS004') then
          k = field%xlvl
          call mpas_log_write('Interpolating SS004 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qseas4
       elseif(trim(field%field) == 'SS005') then
          k = field%xlvl
          call mpas_log_write('Interpolating SS005 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qseas5
       elseif(trim(field%field) == 'DMS') then
          k = field%xlvl
          call mpas_log_write('Interpolating DMS at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qdms
       elseif(trim(field%field) == 'MSA') then
          k = field%xlvl
          call mpas_log_write('Interpolating MSA at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qmsa
       elseif(trim(field%field) == 'HNO3') then
          k = field%xlvl
          call mpas_log_write('Interpolating HNO3 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qhno3
       elseif(trim(field%field) == 'NH3') then
          k = field%xlvl
          call mpas_log_write('Interpolating HN3 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qnh3
       elseif(trim(field%field) == 'DPRES') then
          k = field%xlvl
          call mpas_log_write('Interpolating DPRES at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => dpgoc
       elseif(trim(field%field) == 'PRES') then
          k = field%xlvl
          call mpas_log_write('Interpolating PRES at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => pgoc
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
                destField2d(k,iCell) = interp_sequence(x,y,1,rslab,-2,field%nx+3,1,field%ny,1,1,msgval, \
                                              interp_list,1,maskval=maskval,mask_array=maskslab)
             else
                destField2d(k,iCell) = interp_sequence(x,y,1,rslab,-2,field%nx+3,1,field%ny,1,1,msgval, \
                                              interp_list,1)
             endif
          else
             destField2d(k,iCell) = fillval
          endif
       enddo
       deallocate(rslab)

    endif
    deallocate(field%slab)
    call read_next_met_field(field,istatus)

 enddo
 call read_met_close()


 call mpas_log_write('--- end subroutine init_hinterp_gocart2G.')

 end subroutine init_hinterp_gocart2G

!==================================================================================================================
 subroutine init_vinterp_gocart2G(configs,mesh,fg,diag,state)
!==================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: fg
 type(mpas_pool_type),intent(in):: diag

!inout arguments:
 type(mpas_pool_type),intent(inout):: state

!local variables and arrays:
 integer:: k,iCell,n,nn
 integer,pointer:: nCells,nAerLevels,nVertLevels
 integer,pointer:: num_scalars,gocart2G_start,gocart2G_end
 integer,pointer:: num_scalars_fg,gocart2G_fg_start,gocart2G_fg_end
 integer,pointer:: index_qnh3,index_qnh4a,index_qso4

 real(kind=RKIND),dimension(:,:),pointer:: pgoc,pressure
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars_fg
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars

 real(kind=RKIND):: fMassAir,fMassNH3,fMassNH4a,fMassSO4,fmult
 real(kind=RKIND):: target_p
 real(kind=RKIND),dimension(:,:),allocatable:: sorted_arr

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine init_vinterp_gocart2G:')

 call mpas_pool_get_dimension(mesh,'nCells',nCells)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
 call mpas_pool_get_dimension(mesh,'nAerLevels' ,nAerLevels )

 call mpas_pool_get_array(diag,'pressure',pressure)

 call mpas_pool_get_dimension(fg,'num_scalars_fg',num_scalars_fg   )
 call mpas_pool_get_dimension(fg,'gocart2G_start',gocart2G_fg_start)
 call mpas_pool_get_dimension(fg,'gocart2G_end'  ,gocart2G_fg_end  )
 call mpas_log_write('--- num_scalars_fg     = $i',intArgs=(/num_scalars_fg/)   )
 call mpas_log_write('--- gocart2G_fg_start  = $i',intArgs=(/gocart2G_fg_start/))
 call mpas_log_write('--- gocart2G_fg_end    = $i',intArgs=(/gocart2G_fg_end/)  )
 call mpas_log_write(' ')

 call mpas_pool_get_dimension(state,'num_scalars',num_scalars)
 call mpas_pool_get_dimension(state,'gocart2G_start',gocart2G_start)
 call mpas_pool_get_dimension(state,'gocart2G_end'  ,gocart2G_end  )
 call mpas_log_write('--- num_scalars    = $i',intArgs=(/num_scalars/)   )
 call mpas_log_write('--- gocart2G_start = $i',intArgs=(/gocart2G_start/))
 call mpas_log_write('--- gocart2G_end   = $i',intArgs=(/gocart2G_end/)  )
 call mpas_log_write(' ')

 call mpas_pool_get_array(fg,'pgoc',pgoc)
 call mpas_pool_get_array(fg,'scalars_fg',scalars_fg)
 call mpas_pool_get_array(state,'scalars',scalars)
 do nn = gocart2G_start,gocart2G_end
    scalars(nn,:,:) = 0._RKIND
 enddo


 if(.not.allocated(sorted_arr)) allocate(sorted_arr(2,nAerLevels))
 n = 0
 do nn = gocart2G_start,gocart2G_end
    n = n+1
    do iCell = 1,nCells
       sorted_arr(1,1:nAerLevels) = 0._RKIND
       sorted_arr(2,1:nAerLevels) = 0._RKIND
       do k = 1,nAerLevels
          sorted_arr(1,k) = pgoc(k,iCell)
          sorted_arr(2,k) = scalars_fg(n,k,iCell)
       enddo
       do k = nVertLevels,1,-1
          target_p = pressure(k,iCell)
          scalars(nn,k,iCell) = pressure_interp(iCell,k,target_p,nAerLevels,sorted_arr(:,1:nAerLevels))
          if(target_p.gt.sorted_arr(1,1)) scalars(nn,k,iCell) = scalars(nn,k+1,iCell)
       enddo
    enddo
 enddo
 if(allocated(sorted_arr)) deallocate(sorted_arr)


!--- convert qnh3 from mole per mole to kg per kg. initialize qnh4a as a function of qso4 (we assume
!    that the number of moles of NH4a is equal to the number of moles of SO4, and then convert mole
!    per mole to kg per kg:
 fMassNH3  = 62._RKIND   ! as defined in NI2G_GridCompMod.F90
 fMassNH4a = 18._RKIND   ! as defined in NI2G_GridCompMod.F90
 fMassAir  = 28.97_RKIND ! as defined in NI2G_GridCompMod.F90
 fMassSO4  = 96._RKIND   ! as defined in SU2G_GridCompMod.F90

 call mpas_pool_get_dimension(state,'index_qnh3' ,index_qnh3 )
 call mpas_pool_get_dimension(state,'index_qnh4a',index_qnh4a)
 call mpas_pool_get_dimension(state,'index_qso4' ,index_qso4 )

 fmult = fMassNH3/fMassAir
 scalars(index_qnh3,:,:) = fmult*scalars(index_qnh3,:,:)

 fmult = fMassNH4a/fMassSO4
 scalars(index_qnh4a,:,:) = fmult*scalars(index_qso4,:,:)


 call mpas_log_write('--- end subroutine init_vinterp_gocart2G.')

 end subroutine init_vinterp_gocart2G

!==================================================================================================================
 subroutine init_vinterp_gocart2G_hno3(configs,mesh,fg,diag,gocart2G_backgrounds)
!==================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: fg
 type(mpas_pool_type),intent(in):: diag

!inout arguments:
 type(mpas_pool_type),intent(inout):: gocart2G_backgrounds

!local variables and arrays:
 integer:: k,iCell,n,nn
 integer,pointer:: nCells,nAerLevels,nVertLevels
 integer,pointer:: index_qhno3

 real(kind=RKIND),dimension(:,:),pointer:: background_hno3
 real(kind=RKIND),dimension(:,:),pointer:: pgoc,pressure
 real(kind=RKIND),dimension(:,:),pointer:: qhno3_fg
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars_fg

 real(kind=RKIND):: target_p
 real(kind=RKIND),dimension(:,:),allocatable:: sorted_arr

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine init_vinterp_gocart2G_hno3:')

 call mpas_pool_get_dimension(mesh,'nCells',nCells)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
 call mpas_pool_get_dimension(mesh,'nAerLevels' ,nAerLevels )

 call mpas_pool_get_array(diag,'pressure',pressure)

 call mpas_pool_get_dimension(fg,'index_qhno3',index_qhno3)
 call mpas_pool_get_array(fg,'pgoc',pgoc)
 call mpas_pool_get_array(fg,'scalars_fg',scalars_fg)
 qhno3_fg => scalars_fg(index_qhno3,:,:)

 call mpas_pool_get_array(gocart2G_backgrounds,'background_hno3',background_hno3)

 if(.not.allocated(sorted_arr)) allocate(sorted_arr(2,nAerLevels))
 do iCell = 1,nCells
    sorted_arr(1,1:nAerLevels) = 0._RKIND
    sorted_arr(2,1:nAerLevels) = 0._RKIND
    do k = 1,nAerLevels
       sorted_arr(1,k) = pgoc(k,iCell)
       sorted_arr(2,k) = qhno3_fg(k,iCell)
    enddo
    do k = nVertLevels,1,-1
       target_p = pressure(k,iCell)
       background_hno3(k,iCell) = pressure_interp(iCell,k,target_p,nAerLevels,sorted_arr(:,1:nAerLevels))
       if(target_p.gt.sorted_arr(1,1)) background_hno3(k,iCell) = background_hno3(k+1,iCell)
    enddo
 enddo
 if(allocated(sorted_arr)) deallocate(sorted_arr)

 call mpas_log_write('--- end subroutine init_vinterp_gocart2G_hno3:')

 end subroutine init_vinterp_gocart2G_hno3

!=================================================================================================================
 real(kind=RKIND) function pressure_interp(ii,kk,target_z,nz,zf)
!=================================================================================================================

!input arguments:
 integer,intent(in):: ii,kk
 integer,intent(in):: nz

 real(kind=RKIND),intent(in):: target_z
 real(kind=RKIND),intent(in),dimension(2,nz):: zf

!local variables:
 integer:: k,lm,lp
 real(kind=RKIND):: wm,wp

!-----------------------------------------------------------------------------------------------------------------

 do k = 1,nz-1
    if(target_z <= zf(1,k) .and. target_z > zf(1,k+1)) then
       lm = k
       lp = k+1
       wm = (zf(1,k+1) - target_z)/(zf(1,k+1) - zf(1,k))
       wp = (target_z - zf(1,k))/(zf(1,k+1) - zf(1,k))
       exit
    else
       lm = nz-1
       lp = nz
       wm = 0.
       wp = 0.
    endif
 enddo
 pressure_interp = wm*zf(2,lm) + wp*zf(2,lp)

 return

 end function pressure_interp

!==================================================================================================================
 end module mpas_chemistry_init_gocart2G_interp
!==================================================================================================================
