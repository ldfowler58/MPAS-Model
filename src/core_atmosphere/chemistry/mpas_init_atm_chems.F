!=================================================================================================================
 module mpas_init_atm_chems
 use mpas_dmpar
 use mpas_log
 use mpas_kind_types
 use mpas_pool_routines
 use mpas_stream_manager
 use mpas_derived_types,only     : mpas_time_type,mpas_timeinterval_type
 use mpas_timekeeping,only       : mpas_get_clock_time,mpas_set_time,mpas_get_timeInterval
 use init_atm_read_met,only      : met_data,read_met_init,read_met_close,read_next_met_field
 use init_atm_hinterp,only       : interp_sequence,FOUR_POINT,SEARCH,W_AVERAGE4,W_AVERAGE16
 use init_atm_llxy,only          : latlon_to_ij,map_init,map_set,proj_info,PROJ_LATLON,PROJ_GAUSS,DEG_PER_RAD


 implicit none
 private
 public:: init_atm_CAMS_emissions,init_atm_camchem


!initialization of CAMS emissions and interpolation of CAM-Chem input data to the MPAS mesh.
!Laura D. Fowler (laura@ucar.edu) / 2022-02-08.


 contains


!=================================================================================================================
 subroutine init_atm_CAMS_emissions(clock,stream_manager,mesh,CAMS_emissions)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in),pointer:: mesh
 type(mpas_Clock_type),intent(in),pointer:: clock

!inout arguments:
 type(MPAS_streamManager_type),intent(inout):: stream_manager
 type(mpas_pool_type),intent(inout),pointer:: CAMS_emissions

!local variables and arrays:
 type(mpas_time_type):: beforeTime,afterTime,currTime
 type(mpas_timeinterval_type):: beforeDelta,afterDelta,totalDelta

 character(len=StrKIND):: actualTimeStamp
 integer,pointer:: nCells
 real(kind=RKIND):: total_dt,before_dt,after_dt

 real(kind=RKIND),dimension(:),pointer:: bc1_before,co_before,oc1_before,nh3_before,so2_before
 real(kind=RKIND),dimension(:),pointer:: bc1,co,oc1,nh3,so2

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine init_atm_CAMS_emissions:')

 call mpas_pool_get_dimension(mesh,'nCells',nCells)

 call mpas_pool_get_array(CAMS_emissions,'bc1_em_anthro',bc1)
 call mpas_pool_get_array(CAMS_emissions,'co_em_anthro' ,co )
 call mpas_pool_get_array(CAMS_emissions,'oc1_em_anthro',oc1)
 call mpas_pool_get_array(CAMS_emissions,'nh3_em_anthro',nh3)
 call mpas_pool_get_array(CAMS_emissions,'so2_em_anthro',so2)

!read the latest time slice from the file that is before (or equal to) the current time:
 call mpas_stream_mgr_read(stream_manager,'emissions',rightNow=.true.,whence=MPAS_STREAM_LATEST_BEFORE, &
                           actualWhen=actualTimestamp)
 call mpas_log_write('latest time before is '//trim(actualTimestamp))
 call mpas_log_write('maxval(bc1) = $r', realArgs=[maxval(bc1)])
 call mpas_log_write('maxval(co)  = $r', realArgs=[maxval(co )])
 call mpas_log_write('maxval(oc1) = $r', realArgs=[maxval(oc1)])
 call mpas_log_write('maxval(nh3) = $r', realArgs=[maxval(nh3)])
 call mpas_log_write('maxval(so2) = $r', realArgs=[maxval(so2)])
 call mpas_set_time(beforeTime,dateTimeString=trim(actualTimestamp))

!before reading the next time from the file,copy all arrays to temporary arrays for later use in interpolation:
 call mpas_pool_get_dimension(mesh,'nCells',nCells)
 if(.not.associated(bc1_before)) allocate(bc1_before(nCells+1))  !allocate the nCells+1 garbage cell,too.
 if(.not.associated(co_before) ) allocate(co_before(nCells+1) )
 if(.not.associated(oc1_before)) allocate(oc1_before(nCells+1))
 if(.not.associated(nh3_before)) allocate(nh3_before(nCells+1))
 if(.not.associated(so2_before)) allocate(so2_before(nCells+1))
 bc1_before(:) = bc1(:)
 co_before(:)  = co(:)
 oc1_before(:) = oc1(:)
 nh3_before(:) = nh3(:)
 so2_before(:) = so2(:)

!read the latest time slice from the file that is after (or equal to) the current time:
 call mpas_stream_mgr_read(stream_manager,'emissions',rightNow=.true.,whence=MPAS_STREAM_EARLIEST_AFTER, &
                           actualWhen=actualTimestamp)
 call mpas_log_write('earliest time after is '//trim(actualTimestamp))
 call mpas_log_write('maxval(bc1) = $r', realArgs=[maxval(bc1)])
 call mpas_log_write('maxval(co)  = $r', realArgs=[maxval(co )])
 call mpas_log_write('maxval(oc1) = $r', realArgs=[maxval(oc1)])
 call mpas_log_write('maxval(nh3) = $r', realArgs=[maxval(nh3)])
 call mpas_log_write('maxval(so2) = $r', realArgs=[maxval(so2)])
 call mpas_set_time(afterTime,dateTimeString=trim(actualTimestamp))

!get current time:
 currTime = mpas_get_clock_time(clock,MPAS_NOW)

!calculate time deltas between the times that were actually read and the current time:
 totalDelta  = afterTime - beforeTime
 beforeDelta = currTime - beforeTime
 afterDelta  = afterTime - currTime

!retrieve time deltas as real values:
 call mpas_get_timeInterval(totalDelta ,dt=total_dt)
 call mpas_get_timeInterval(beforeDelta,dt=before_dt)
 call mpas_get_timeInterval(afterDelta , dt=after_dt)
 call mpas_log_write(' ')
 call mpas_log_write('--- totalDelta  = $r',realArgs=(/total_dt/))
 call mpas_log_write('--- beforeDelta = $r',realArgs=(/before_dt/))
 call mpas_log_write('--- afterDelta  = $r',realArgs=(/after_dt/))

!interpolation of surface emissions to the current time:
 if(total_dt > 0.0_RKIND) then
    bc1(:) = (after_dt/total_dt)*bc1_before(:) + (before_dt/total_dt)*bc1(:)
    co(:)  = (after_dt/total_dt)*co_before(:)  + (before_dt/total_dt)*co(:)
    oc1(:) = (after_dt/total_dt)*oc1_before(:) + (before_dt/total_dt)*oc1(:)
    nh3(:) = (after_dt/total_dt)*nh3_before(:) + (before_dt/total_dt)*nh3(:)
    so2(:) = (after_dt/total_dt)*so2_before(:) + (before_dt/total_dt)*so2(:)
    call mpas_log_write(' ')
    call mpas_log_write('maxval(bc1) = $r', realArgs=[maxval(bc1)])
    call mpas_log_write('maxval(co)  = $r', realArgs=[maxval(co )])
    call mpas_log_write('maxval(oc1) = $r', realArgs=[maxval(oc1)])
    call mpas_log_write('maxval(nh3) = $r', realArgs=[maxval(nh3)])
    call mpas_log_write('maxval(so2) = $r', realArgs=[maxval(so2)])
 endif

 if(associated(bc1_before)) deallocate(bc1_before)  !allocate the nCells+1 garbage cell,too.
 if(associated(co_before) ) deallocate(co_before )
 if(associated(oc1_before)) deallocate(oc1_before)
 if(associated(nh3_before)) deallocate(nh3_before)
 if(associated(so2_before)) deallocate(so2_before)

 call mpas_log_write('--- end subroutine init_atm_CAMS_emissions.')

 end subroutine init_atm_CAMS_emissions

!=================================================================================================================
 subroutine init_atm_camchem(configs,mesh,fg,diag,chems_input,state)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: diag

!inout arguments:
 type(mpas_pool_type),intent(inout):: fg
 type(mpas_pool_type),intent(inout):: chems_input
 type(mpas_pool_type),intent(inout):: state

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine init_atm_camchem:')

 call init_hinterp_camchem(configs,mesh,fg)
 call init_vinterp_camchem(configs,mesh,fg,diag,chems_input)
 call init_camchem_tothompson(mesh,chems_input,state)

 call mpas_log_write('--- end subroutine init_atm_camchem.')
 call mpas_log_write(' ')

 end subroutine init_atm_camchem

!=================================================================================================================
 subroutine init_hinterp_camchem(configs,mesh,fg)
!=================================================================================================================

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

 integer,pointer:: nCells,nCAMLevels
 integer:: iCell,istatus,k,kk,masked,nInterpPoints
 integer,dimension(5):: interp_list

 integer,pointer:: num_chems_fg
 integer,pointer:: gocart_fg_start,gocart_fg_end
 integer,pointer:: index_qco,index_qnh3,index_qso2,index_qso4a1,index_qso4a2,index_qso4a3,index_qbc1,    &
                   index_qbc2,index_qoc1,index_qoc2,index_qseas1,index_qseas2,index_qseas3,index_qseas4, &
                   index_qdust1,index_qdust2,index_qdust3,index_qdust4,index_qdust5
 integer,dimension(:),pointer:: landmask,mask_array

 real(kind=RKIND):: fillval,maskval,msgval
 real(kind=RKIND):: lat,lon,x,y
 real(kind=RKIND),dimension(:),pointer:: latCell,lonCell
 real(kind=RKIND),dimension(:),pointer:: latPoints,lonPoints
 real(kind=RKIND),dimension(:,:),pointer:: qco_fg,qnh3_fg,qso2_fg,qso4a1_fg,qso4a2_fg,qso4a3_fg,qbc1_fg,    &
                                           qbc2_fg,qoc1_fg,qoc2_fg,qseas1_fg,qseas2_fg,qseas3_fg,qseas4_fg, &
                                           qdust1_fg,qdust2_fg,qdust3_fg,qdust4_fg,qdust5_fg,pcam_fg
 real(kind=RKIND),dimension(:,:),pointer:: destField2d
 real(kind=RKIND),dimension(:,:,:),pointer:: chems_fg

 real(kind=RKIND),dimension(:,:),pointer:: pcam_r
 real(kind=RKIND),dimension(:,:,:),pointer:: chems_r

 real(kind=RKIND),dimension(:,:),allocatable:: maskslab,rslab

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine init_hinterp_camchem:')

 call mpas_pool_get_config(configs,'config_chems_prefix',prefix    )
 call mpas_pool_get_config(configs,'config_start_time'  ,start_time)

 call mpas_pool_get_dimension(mesh,'nCells'    ,nCells    )
 call mpas_pool_get_dimension(mesh,'nCAMLevels',nCAMLevels)

 call mpas_pool_get_array(mesh,'landmask',landmask)
 call mpas_pool_get_array(mesh,'latCell' ,latCell )
 call mpas_pool_get_array(mesh,'lonCell' ,lonCell )

 call mpas_pool_get_dimension(fg,'num_chems_fg'   ,num_chems_fg  )
 call mpas_pool_get_dimension(fg,'gocart_fg_start',gocart_fg_start)
 call mpas_pool_get_dimension(fg,'gocart_fg_end'  ,gocart_fg_end  )
 call mpas_log_write('--- num_chems_fg    = $i',intArgs=(/num_chems_fg/)   )
 call mpas_log_write('--- gocart_fg_start = $i',intArgs=(/gocart_fg_start/))
 call mpas_log_write('--- gocart_fg_end   = $i',intArgs=(/gocart_fg_end/)  )
 call mpas_log_write(' ')

 call mpas_pool_get_dimension(fg,'index_qco_fg'   ,index_qco   )
 call mpas_pool_get_dimension(fg,'index_qnh3_fg'  ,index_qnh3  )
 call mpas_pool_get_dimension(fg,'index_qso2_fg'  ,index_qso2  )
 call mpas_pool_get_dimension(fg,'index_qso4a1_fg',index_qso4a1)
 call mpas_pool_get_dimension(fg,'index_qso4a2_fg',index_qso4a2)
 call mpas_pool_get_dimension(fg,'index_qso4a3_fg',index_qso4a3)
 call mpas_pool_get_dimension(fg,'index_qbc1_fg'  ,index_qbc1  )
 call mpas_pool_get_dimension(fg,'index_qbc2_fg'  ,index_qbc2  )
 call mpas_pool_get_dimension(fg,'index_qoc1_fg'  ,index_qoc1  )
 call mpas_pool_get_dimension(fg,'index_qoc2_fg'  ,index_qoc2  )
 call mpas_pool_get_dimension(fg,'index_qseas1_fg',index_qseas1)
 call mpas_pool_get_dimension(fg,'index_qseas2_fg',index_qseas2)
 call mpas_pool_get_dimension(fg,'index_qseas3_fg',index_qseas3)
 call mpas_pool_get_dimension(fg,'index_qseas4_fg',index_qseas4)
 call mpas_pool_get_dimension(fg,'index_qdust1_fg',index_qdust1)
 call mpas_pool_get_dimension(fg,'index_qdust2_fg',index_qdust2)
 call mpas_pool_get_dimension(fg,'index_qdust3_fg',index_qdust3)
 call mpas_pool_get_dimension(fg,'index_qdust4_fg',index_qdust4)
 call mpas_pool_get_dimension(fg,'index_qdust5_fg',index_qdust5)
 call mpas_log_write('--- index_qbc1_fg   = $i',intArgs=(/index_qbc1/)  )
 call mpas_log_write('--- index_qbc2_fg   = $i',intArgs=(/index_qbc2/)  )
 call mpas_log_write('--- index_qoc1_fg   = $i',intArgs=(/index_qoc1/)  )
 call mpas_log_write('--- index_qoc2_fg   = $i',intArgs=(/index_qoc2/)  )
 call mpas_log_write('--- index_qseas1_fg = $i',intArgs=(/index_qseas1/))
 call mpas_log_write('--- index_qseas2_fg = $i',intArgs=(/index_qseas2/))
 call mpas_log_write('--- index_qseas3_fg = $i',intArgs=(/index_qseas3/))
 call mpas_log_write('--- index_qseas4_fg = $i',intArgs=(/index_qseas4/))
 call mpas_log_write('--- index_qdust1_fg = $i',intArgs=(/index_qdust1/))
 call mpas_log_write('--- index_qdust2_fg = $i',intArgs=(/index_qdust2/))
 call mpas_log_write('--- index_qdust3_fg = $i',intArgs=(/index_qdust3/))
 call mpas_log_write('--- index_qdust4_fg = $i',intArgs=(/index_qdust4/))
 call mpas_log_write('--- index_qdust5_fg = $i',intArgs=(/index_qdust5/))
 call mpas_log_write('--- index_qco_fg    = $i',intArgs=(/index_qco/)   )
 call mpas_log_write('--- index_qnh3_fg   = $i',intArgs=(/index_qnh3/)  )
 call mpas_log_write('--- index_qso2_fg   = $i',intArgs=(/index_qso2/)  )
 call mpas_log_write('--- index_qso4a1_fg = $i',intArgs=(/index_qso4a1/))
 call mpas_log_write('--- index_qso4a2_fg = $i',intArgs=(/index_qso4a2/))
 call mpas_log_write('--- index_qso4a3_fg = $i',intArgs=(/index_qso4a3/))
 call mpas_log_write(' ')

 call mpas_pool_get_array(fg,'chems_fg',chems_fg)
 qco_fg    => chems_fg(index_qco,:,: )
 qnh3_fg   => chems_fg(index_qnh3,:,:)
 qso2_fg   => chems_fg(index_qso2,:,:)
 qso4a1_fg => chems_fg(index_qso4a1,:,:)
 qso4a2_fg => chems_fg(index_qso4a2,:,:)
 qso4a3_fg => chems_fg(index_qso4a3,:,:)
 qbc1_fg   => chems_fg(index_qbc1,:,:)
 qbc2_fg   => chems_fg(index_qbc2,:,:)
 qoc1_fg   => chems_fg(index_qoc1,:,:)
 qoc2_fg   => chems_fg(index_qoc2,:,:)
 qseas1_fg => chems_fg(index_qseas1,:,:)
 qseas2_fg => chems_fg(index_qseas2,:,:)
 qseas3_fg => chems_fg(index_qseas3,:,:)
 qseas4_fg => chems_fg(index_qseas4,:,:)
 qdust1_fg => chems_fg(index_qdust1,:,:)
 qdust2_fg => chems_fg(index_qdust2,:,:)
 qdust3_fg => chems_fg(index_qdust3,:,:)
 qdust4_fg => chems_fg(index_qdust4,:,:)
 qdust5_fg => chems_fg(index_qdust5,:,:)
 chems_fg(:,:,:) = 0._RKIND

 call mpas_pool_get_array(fg,'pcam_fg',pcam_fg)


!open intermediate file:
 istatus = 0
 call read_met_init(trim(prefix),.false.,start_time(1:13),istatus)
 if(istatus /= 0) then
    call mpas_log_write('**************************************************',messageType=MPAS_LOG_ERR)
    call mpas_log_write('Error opening intermediate CAM-Chem data file ' &
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


!horizontally interpolate CAM-CHEM data:
 istatus = 0
 call read_met_init(trim(prefix),.false.,start_time(1:13),istatus)
 if(istatus /= 0) then
    call mpas_log_write('**************************************************',messageType=MPAS_LOG_ERR)
    call mpas_log_write('Error opening intermediate CAM-Chem data file ' &
                                       //trim(prefix)//':'//start_time(1:13),messageType=MPAS_LOG_ERR)
    call mpas_log_write('**************************************************',messageType=MPAS_LOG_CRIT)
 endif
 call read_next_met_field(field,istatus)


 do while(istatus == 0)

    interp_list(1) = FOUR_POINT
    interp_list(2) = W_AVERAGE4
    interp_list(3) = W_AVERAGE16
    interp_list(4) = SEARCH
    interp_list(5) = 0

    maskval = -1.0
    masked  = -1
    fillval = 0.0
    msgval  = -1.e30

    mask_array => landmask

    if(index(field % field, 'CO'    ) /= 0 .or. &
       index(field % field, 'NH3'   ) /= 0 .or. &
       index(field % field, 'SO2'   ) /= 0 .or. &
       index(field % field, 'SO4A1' ) /= 0 .or. &
       index(field % field, 'SO4A2' ) /= 0 .or. &
       index(field % field, 'SO4A3' ) /= 0 .or. &
       index(field % field, 'BC1'   ) /= 0 .or. &
       index(field % field, 'BC2'   ) /= 0 .or. &
       index(field % field, 'OC1'   ) /= 0 .or. &
       index(field % field, 'OC2'   ) /= 0 .or. &
       index(field % field, 'SEAS_1') /= 0 .or. &
       index(field % field, 'SEAS_2') /= 0 .or. &
       index(field % field, 'SEAS_3') /= 0 .or. &
       index(field % field, 'SEAS_4') /= 0 .or. &
       index(field % field, 'DUST_1') /= 0 .or. &
       index(field % field, 'DUST_2') /= 0 .or. &
       index(field % field, 'DUST_3') /= 0 .or. &
       index(field % field, 'DUST_4') /= 0 .or. &
       index(field % field, 'DUST_5') /= 0 .or. &
       index(field % field, 'PCAM'  ) /= 0 ) then

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
       if(index(field%field,'CO') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating CO at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qco_fg
       elseif(index(field%field,'NH3') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating NH3 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qnh3_fg
       elseif(index(field%field,'SO2') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating SO2 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qso2_fg
       elseif(index(field%field,'SO4A1') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating SO4A1 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qso4a1_fg
       elseif(index(field%field,'SO4A2') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating SO4A1 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qso4a2_fg
       elseif(index(field%field,'SO4A3') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating SO4A1 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qso4a3_fg
       elseif(index(field%field,'BC1') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating BC1 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qbc1_fg
       elseif(index(field%field,'BC2') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating BC2 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qbc2_fg
       elseif(index(field%field,'OC1') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating OC1 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qoc1_fg
       elseif(index(field%field,'OC2') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating OC2 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qoc2_fg
       elseif(index(field%field,'SEAS_1') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating SEAS_1 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qseas1_fg
       elseif(index(field%field,'SEAS_2') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating SEAS_2 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qseas2_fg
       elseif(index(field%field,'SEAS_3') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating SEAS_3 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qseas3_fg
       elseif(index(field%field,'SEAS_4') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating SEAS_4 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qseas4_fg
       elseif(index(field%field,'DUST_1') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating DUST_1 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qdust1_fg
       elseif(index(field%field,'DUST_2') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating DUST_2 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qdust2_fg
       elseif(index(field%field,'DUST_3') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating DUST_3 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qdust3_fg
       elseif(index(field%field,'DUST_4') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating DUST_4 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qdust4_fg
       elseif(index(field%field,'DUST_5') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating DUST_5 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qdust5_fg
       elseif(index(field%field,'PCAM') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating PCAM at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => pcam_fg
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

 call mpas_log_write(' ')
 call mpas_log_write('--- check for negative qbc1_fg:')
 do iCell = 1,nCells
    do k = 1,nCAMLevels
       if(qbc1_fg(k,iCell) .lt. 0._RKIND) then
          call mpas_log_write('--- negative qbc1_fg: $i $i $r',intArgs=(/iCell,k/), &
                              realArgs=(/qbc1_fg(k,iCell)/))
       endif
    enddo
 enddo
 call mpas_log_write('--- check for negative qbc2_fg:')
 do iCell = 1,nCells
    do k = 1,nCAMLevels
       if(qbc2_fg(k,iCell) .lt. 0._RKIND) then
          call mpas_log_write('--- negative qbc2_fg: $i $i $r',intArgs=(/iCell,k/), &
                              realArgs=(/qbc2_fg(k,iCell)/))
       endif
    enddo
 enddo
 call mpas_log_write('--- check for negative qoc1_fg:')
 do iCell = 1,nCells
    do k = 1,nCAMLevels
       if(qoc1_fg(k,iCell) .lt. 0._RKIND) then
          call mpas_log_write('--- negative qoc1_fg: $i $i $r',intArgs=(/iCell,k/), &
                              realArgs=(/qoc1_fg(k,iCell)/))
       endif
    enddo
 enddo
 call mpas_log_write('--- check for negative qoc2_fg:')
 do iCell = 1,nCells
    do k = 1,nCAMLevels
       if(qoc2_fg(k,iCell) .lt. 0._RKIND) then
          call mpas_log_write('--- negative qoc2_fg: $i $i $r',intArgs=(/iCell,k/), &
                              realArgs=(/qoc2_fg(k,iCell)/))
       endif
    enddo
 enddo
 call mpas_log_write('--- check for negative qseas1_fg:')
 do iCell = 1,nCells
    do k = 1,nCAMLevels
       if(qseas1_fg(k,iCell) .lt. 0._RKIND) then
          call mpas_log_write('--- negative qseas1_fg: $i $i $r',intArgs=(/iCell,k/), &
                              realArgs=(/qseas1_fg(k,iCell)/))
       endif
    enddo
 enddo
 call mpas_log_write('--- check for negative qseas2_fg:')
 do iCell = 1,nCells
    do k = 1,nCAMLevels
       if(qseas2_fg(k,iCell) .lt. 0._RKIND) then
          call mpas_log_write('--- negative qseas2_fg: $i $i $r',intArgs=(/iCell,k/), &
                              realArgs=(/qseas2_fg(k,iCell)/))
       endif
    enddo
 enddo
 call mpas_log_write('--- check for negative qseas3_fg:')
 do iCell = 1,nCells
    do k = 1,nCAMLevels
       if(qseas3_fg(k,iCell) .lt. 0._RKIND) then
          call mpas_log_write('--- negative qseas3_fg: $i $i $r',intArgs=(/iCell,k/), &
                              realArgs=(/qseas3_fg(k,iCell)/))
       endif
    enddo
 enddo
 call mpas_log_write('--- check for negative qseas4_fg:')
 do iCell = 1,nCells
    do k = 1,nCAMLevels
       if(qseas4_fg(k,iCell) .lt. 0._RKIND) then
          call mpas_log_write('--- negative qseas4_fg: $i $i $r',intArgs=(/iCell,k/), &
                              realArgs=(/qseas4_fg(k,iCell)/))
       endif
    enddo
 enddo
 call mpas_log_write('--- check for negative qdust1_fg:')
 do iCell = 1,nCells
    do k = 1,nCAMLevels
       if(qdust1_fg(k,iCell) .lt. 0._RKIND) then
          call mpas_log_write('--- negative qdust1_fg: $i $i $r',intArgs=(/iCell,k/), &
                              realArgs=(/qdust1_fg(k,iCell)/))
       endif
    enddo
 enddo
 call mpas_log_write('--- check for negative qdust2_fg:')
 do iCell = 1,nCells
    do k = 1,nCAMLevels
       if(qdust2_fg(k,iCell) .lt. 0._RKIND) then
          call mpas_log_write('--- negative qdust2_fg: $i $i $r',intArgs=(/iCell,k/), &
                              realArgs=(/qdust2_fg(k,iCell)/))
       endif
    enddo
 enddo
 call mpas_log_write('--- check for negative qdust3_fg:')
 do iCell = 1,nCells
    do k = 1,nCAMLevels
       if(qdust3_fg(k,iCell) .lt. 0._RKIND) then
          call mpas_log_write('--- negative qdust3_fg: $i $i $r',intArgs=(/iCell,k/), &
                              realArgs=(/qdust3_fg(k,iCell)/))
       endif
    enddo
 enddo
 call mpas_log_write('--- check for negative qdust4_fg:')
 do iCell = 1,nCells
    do k = 1,nCAMLevels
       if(qdust4_fg(k,iCell) .lt. 0._RKIND) then
          call mpas_log_write('--- negative qdust4_fg: $i $i $r',intArgs=(/iCell,k/), &
                              realArgs=(/qdust4_fg(k,iCell)/))
       endif
    enddo
 enddo
 call mpas_log_write('--- check for negative qdust5_fg:')
 do iCell = 1,nCells
    do k = 1,nCAMLevels
       if(qdust5_fg(k,iCell) .lt. 0._RKIND) then
          call mpas_log_write('--- negative qdust5_fg: $i $i $r',intArgs=(/iCell,k/), &
                              realArgs=(/qdust5_fg(k,iCell)/))
       endif
    enddo
 enddo
 call mpas_log_write('--- check for negative qco_fg:')
 do iCell = 1,nCells
    do k = 1,nCAMLevels
       if(qco_fg(k,iCell) .lt. 0._RKIND) then
          call mpas_log_write('--- negative qco_fg: $i $i $r',intArgs=(/iCell,k/), &
                              realArgs=(/qco_fg(k,iCell)/))
       endif
    enddo
 enddo
 call mpas_log_write('--- check for negative qnh3_fg:')
 do iCell = 1,nCells
    do k = 1,nCAMLevels
       if(qnh3_fg(k,iCell) .lt. 0._RKIND) then
          call mpas_log_write('--- negative qnh3_fg: $i $i $r',intArgs=(/iCell,k/), &
                              realArgs=(/qnh3_fg(k,iCell)/))
       endif
    enddo
 enddo
 call mpas_log_write('--- check for negative qso2_fg:')
 do iCell = 1,nCells
    do k = 1,nCAMLevels
       if(qso2_fg(k,iCell) .lt. 0._RKIND) then
          call mpas_log_write('--- negative qso2_fg: $i $i $r',intArgs=(/iCell,k/), &
                              realArgs=(/qso2_fg(k,iCell)/))
       endif
    enddo
 enddo
 call mpas_log_write('--- check for negative qso4a1_fg:')
 do iCell = 1,nCells
    do k = 1,nCAMLevels
       if(qso4a1_fg(k,iCell) .lt. 0._RKIND) then
          call mpas_log_write('--- negative qso4a1_fg: $i $i $r',intArgs=(/iCell,k/), &
                              realArgs=(/qso4a1_fg(k,iCell)/))
       endif
    enddo
 enddo
 call mpas_log_write('--- check for negative qso4a2_fg:')
 do iCell = 1,nCells
    do k = 1,nCAMLevels
       if(qso4a2_fg(k,iCell) .lt. 0._RKIND) then
          call mpas_log_write('--- negative qso4a2_fg: $i $i $r',intArgs=(/iCell,k/), &
                              realArgs=(/qso4a2_fg(k,iCell)/))
       endif
    enddo
 enddo
 call mpas_log_write('--- check for negative qso4a3_fg:')
 do iCell = 1,nCells
    do k = 1,nCAMLevels
       if(qso4a3_fg(k,iCell) .lt. 0._RKIND) then
          call mpas_log_write('--- negative qso4a3_fg: $i $i $r',intArgs=(/iCell,k/), &
                              realArgs=(/qso4a3_fg(k,iCell)/))
       endif
    enddo
 enddo
 call mpas_log_write('--- check for negative pcam_fg:')
 do iCell = 1,nCells
    do k = 1,nCAMLevels
       if(pcam_fg(k,iCell) .lt. 0._RKIND) then
          call mpas_log_write('--- negative pcam_fg: $i $i $r',intArgs=(/iCell,k/), &
                              realArgs=(/pcam_fg(k,iCell)/))
       endif
    enddo
 enddo

 call mpas_log_write('--- exit subroutine init_hinterp_camchem.')
 call mpas_log_write(' ')

 end subroutine init_hinterp_camchem

!=================================================================================================================
 subroutine init_vinterp_camchem(configs,mesh,fg,diag,chems_input)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: fg
 type(mpas_pool_type),intent(in):: diag

!inout arguments:
 type(mpas_pool_type),intent(inout):: chems_input

!local variables and arrays:
 integer:: iCell,k,kk,n

 integer,pointer:: nCells,nCAMLevels,nVertLevels
 integer,pointer:: num_chems,num_chems_fg
 integer,pointer:: gocart_start,gocart_end,gocart_fg_start,gocart_fg_end
 integer,pointer:: index_qbc1,index_qbc2,index_qoc1,index_qoc2,index_qseas1,index_qseas2,index_qseas3, &
                   index_qseas4,index_qdust1,index_qdust2,index_qdust3,index_qdust4,index_qdust5

 real(kind=RKIND),dimension(:,:),pointer:: pcam_fg,pressure
 real(kind=RKIND),dimension(:,:,:),pointer:: chems,chems_fg

 real(kind=RKIND):: target_p
 real(kind=RKIND),dimension(:,:),allocatable:: sorted_arr

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine init_vinterp_camchem:')

 call mpas_pool_get_dimension(mesh,'nCells',nCells)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
 call mpas_pool_get_dimension(mesh,'nCAMLevels' ,nCAMLevels )

 call mpas_pool_get_dimension(chems_input,'num_chems_gocart',num_chems   )
 call mpas_pool_get_dimension(chems_input,'gocart_start'    ,gocart_start)
 call mpas_pool_get_dimension(chems_input,'gocart_end'      ,gocart_end  )
 call mpas_log_write('--- num_chems_gocart = $i',intArgs=(/num_chems/)   )
 call mpas_log_write('--- gocart_start     = $i',intArgs=(/gocart_start/))
 call mpas_log_write('--- gocart_end       = $i',intArgs=(/gocart_end/)  )
 call mpas_log_write(' ')

 call mpas_pool_get_dimension(fg,'num_chems_fg'   ,num_chems_fg   )
 call mpas_pool_get_dimension(fg,'gocart_fg_start',gocart_fg_start)
 call mpas_pool_get_dimension(fg,'gocart_fg_end'  ,gocart_fg_end  )
 call mpas_log_write('--- num_chems_fg    = $i',intArgs=(/num_chems_fg/)   )
 call mpas_log_write('--- gocart_fg_start = $i',intArgs=(/gocart_fg_start/))
 call mpas_log_write('--- gocart_fg_end   = $i',intArgs=(/gocart_fg_end/)  )
 call mpas_log_write(' ')

 call mpas_pool_get_array(diag,'pressure_base',pressure)

 call mpas_pool_get_array(fg,'pcam_fg' ,pcam_fg )
 call mpas_pool_get_array(fg,'chems_fg',chems_fg) 
 call mpas_pool_get_array(chems_input,'chems_gocart',chems)
 chems(:,:,:) = 0._RKIND

 if(.not.allocated(sorted_arr)) allocate(sorted_arr(2,nCAMLevels))
 do n = 1,num_chems_fg
    do iCell = 1,nCells
       sorted_arr(1,1:nCAMLevels) = 0._RKIND
       sorted_arr(1,1:nCAMLevels) = 0._RKIND
       do k = 1,nCAMLevels
          sorted_arr(1,k) = pcam_fg(k,iCell)
          sorted_arr(2,k) = chems_fg(n,k,iCell)
       enddo
       do k = 1,nVertLevels
          target_p = pressure(k,iCell)
          chems(n,k,iCell) = pressure_interp(target_p,nCAMLevels,sorted_arr(:,1:nCAMLevels), &
                                             order=1,extrap=0)
       enddo
    enddo
 enddo
 if(allocated(sorted_arr)) deallocate(sorted_arr)

 call mpas_log_write('--- negative chem variables:')
 do n = 1,num_chems
    do iCell = 1,nCells
       do k = 1,nVertLevels
          if(chems(n,k,iCell) .lt. 0._RKIND) then
             call mpas_log_write('$i $i $i $r',intArgs=(/n,iCell,k/),realArgs=(/chems(n,k,iCell)/))
          endif
       enddo
    enddo
 enddo

 call mpas_log_write('--- end subroutine init_vinterp_camchem.')

 end subroutine init_vinterp_camchem

!=================================================================================================================
 subroutine init_camchem_tothompson(mesh,chems_input,state)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh

!inout arguments:
 type(mpas_pool_type),intent(inout):: chems_input
 type(mpas_pool_type),intent(inout):: state

!local variables and arrays:
 integer:: iCell,k
 integer,pointer:: nCells,nVertLevels
 integer,pointer:: index_nifa,index_nwfa,index_qifa,index_qwfa
 integer,pointer:: index_qoc2,index_qseas1,index_qseas2,index_qseas3,index_qseas4,index_qso4a1, &
                   index_qso4a2,index_qso4a3
 integer,pointer:: index_qdust1,index_qdust2,index_qdust3,index_qdust4,index_qdust5

 real(kind=RKIND),parameter:: &
    conv_dust1 = 1.58e-15, & ! conversion factor from mixing ratio to number concentration for dust1.
    conv_dust2 = 1.68e-14, & ! conversion factor from mixing ratio to number concentration for dust2.
    conv_dust3 = 8.46e-14, & ! conversion factor from mixing ratio to number concentration for dust3.
    conv_dust4 = 5.58e-13, & ! conversion factor from mixing ratio to number concentration for dust4.
    conv_dust5 = 3.13e-13    ! conversion factor from mixing ratio to number concentration for dust5.
 real(kind=RKIND),parameter:: &
    conv_oc2   = 6.61e-18    ! conversion factor from mixing ratio to number concentration for oc2.
 real(kind=RKIND),parameter:: &
    conv_seas1 = 1.04e-15, & ! conversion factor from mixing ratio to number concentration for seas1.
    conv_seas2 = 3.88e-13    ! conversion factor from mixing ratio to number concentration for sea2.
 real(kind=RKIND),parameter:: &
    conv_so4   = 1.07e-17    ! conversion factor from mixing ratio to number concentration for so4.
 real(kind=RKIND):: iconv_dust1,iconv_dust2,iconv_dust3,iconv_dust4,iconv_dust5
 real(kind=RKIND):: iconv_oc2
 real(kind=RKIND):: iconv_seas1,iconv_seas2
 real(kind=RKIND):: iconv_so4

 real(kind=RKIND),dimension(:,:),pointer:: nifa,nwfa,qifa,qwfa
 real(kind=RKIND),dimension(:,:),pointer:: qoc2,qseas1,qseas2,qseas3,qseas4,qso4a1,qso4a2,qso4a3
 real(kind=RKIND),dimension(:,:),pointer:: qdust1,qdust2,qdust3,qdust4,qdust5
 real(kind=RKIND),dimension(:,:,:),pointer:: chems,scalars

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine init_camchem_tothompson:')

 call mpas_pool_get_dimension(mesh,'nCells'     ,nCells     )
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)

 call mpas_pool_get_dimension(chems_input,'index_qoc2'  ,index_qoc2  )
 call mpas_pool_get_dimension(chems_input,'index_qseas1',index_qseas1)
 call mpas_pool_get_dimension(chems_input,'index_qseas2',index_qseas2)
 call mpas_pool_get_dimension(chems_input,'index_qseas3',index_qseas3)
 call mpas_pool_get_dimension(chems_input,'index_qseas4',index_qseas4)
 call mpas_pool_get_dimension(chems_input,'index_qso4a1',index_qso4a1)
 call mpas_pool_get_dimension(chems_input,'index_qso4a2',index_qso4a2)
 call mpas_pool_get_dimension(chems_input,'index_qso4a3',index_qso4a3)
 call mpas_pool_get_dimension(chems_input,'index_qdust1',index_qdust1)
 call mpas_pool_get_dimension(chems_input,'index_qdust2',index_qdust2)
 call mpas_pool_get_dimension(chems_input,'index_qdust3',index_qdust3)
 call mpas_pool_get_dimension(chems_input,'index_qdust4',index_qdust4)
 call mpas_pool_get_dimension(chems_input,'index_qdust5',index_qdust5)
 call mpas_pool_get_dimension(chems_input,'index_qifa'  ,index_qifa  )
 call mpas_pool_get_dimension(chems_input,'index_qwfa'  ,index_qwfa  )
 call mpas_log_write('--- index_qoc2   = $i',intArgs=(/index_qoc2/)  )
 call mpas_log_write('--- index_qseas1 = $i',intArgs=(/index_qseas1/))
 call mpas_log_write('--- index_qseas2 = $i',intArgs=(/index_qseas2/))
 call mpas_log_write('--- index_qseas3 = $i',intArgs=(/index_qseas3/))
 call mpas_log_write('--- index_qseas4 = $i',intArgs=(/index_qseas4/))
 call mpas_log_write('--- index_qso4a1 = $i',intArgs=(/index_qso4a1/))
 call mpas_log_write('--- index_qso4a2 = $i',intArgs=(/index_qso4a2/))
 call mpas_log_write('--- index_qso4a3 = $i',intArgs=(/index_qso4a3/))
 call mpas_log_write('--- index_qdust1 = $i',intArgs=(/index_qdust1/))
 call mpas_log_write('--- index_qdust2 = $i',intArgs=(/index_qdust2/))
 call mpas_log_write('--- index_qdust3 = $i',intArgs=(/index_qdust3/))
 call mpas_log_write('--- index_qdust4 = $i',intArgs=(/index_qdust4/))
 call mpas_log_write('--- index_qdust5 = $i',intArgs=(/index_qdust5/))
 call mpas_log_write('--- index_qifa   = $i',intArgs=(/index_qifa/)  )
 call mpas_log_write('--- index_qwfa   = $i',intArgs=(/index_qwfa/)  )

 call mpas_pool_get_array(chems_input,'chems_gocart',chems)
 qoc2   => chems(index_qoc2,:,:)
 qseas1 => chems(index_qseas1,:,:)
 qseas2 => chems(index_qseas2,:,:)
 qseas3 => chems(index_qseas3,:,:)
 qseas4 => chems(index_qseas4,:,:)
 qso4a1 => chems(index_qso4a1,:,:)
 qso4a2 => chems(index_qso4a2,:,:)
 qso4a3 => chems(index_qso4a3,:,:)
 qdust1 => chems(index_qdust1,:,:)
 qdust2 => chems(index_qdust2,:,:)
 qdust3 => chems(index_qdust3,:,:)
 qdust4 => chems(index_qdust4,:,:)
 qdust5 => chems(index_qdust5,:,:)
 qifa   => chems(index_qifa,:,:)
 qwfa   => chems(index_qwfa,:,:)

 call mpas_pool_get_dimension(state,'index_nifa',index_nifa)
 call mpas_pool_get_dimension(state,'index_nwfa',index_nwfa)
 call mpas_log_write(' ')
 call mpas_log_write('--- index_nifa   = $i',intArgs=(/index_nifa/))
 call mpas_log_write('--- index_nwfa   = $i',intArgs=(/index_nwfa/))

 call mpas_pool_get_array(state,'scalars',scalars)
 nifa => scalars(index_nifa,:,:)
 nwfa => scalars(index_nwfa,:,:)


!--- initialization:
 iconv_dust1 = 1._RKIND/conv_dust1
 iconv_dust2 = 1._RKIND/conv_dust2
 iconv_dust3 = 1._RKIND/conv_dust3
 iconv_dust4 = 1._RKIND/conv_dust4
 iconv_dust5 = 1._RKIND/conv_dust5

 iconv_oc2   = 1._RKIND/conv_oc2
 iconv_seas1 = 1._RKIND/conv_seas1
 iconv_seas2 = 1._RKIND/conv_seas2

 iconv_so4   = 1._RKIND/conv_so4

 do iCell = 1,nCells
    do k = 1,nVertLevels
       !"ice-friendly" aerosols:
       qifa(k,iCell) = qdust1(k,iCell) + qdust2(k,iCell) + qdust3(k,iCell) &
                     + qdust4(k,iCell) + qdust5(k,iCell)

       nifa(k,iCell) = iconv_dust1*qdust1(k,iCell) + iconv_dust2*qdust2(k,iCell) &
                     + iconv_dust3*qdust3(k,iCell) + iconv_dust4*qdust4(k,iCell) &
                     + iconv_dust5*qdust5(k,iCell)

       !"water-friendly" aerosols:
       qwfa(k,iCell) = 8*qoc2(k,iCell) + 5*qso4a1(k,iCell) + 5*qso4a2(k,iCell) &
                     + 9*qseas1(k,iCell) + 9*qseas2(k,iCell)

       nwfa(k,iCell) = 8*iconv_oc2*qoc2(k,iCell) + 5*iconv_so4*qso4a1(k,iCell) + 5*iconv_so4*qso4a2(k,iCell) &
                     + 9*iconv_seas1*qseas1(k,iCell) + 9*iconv_seas2*qseas2(k,iCell)
    enddo
 enddo

 call mpas_log_write('--- end subroutine init_camchem_tothompson:')

 end subroutine init_camchem_tothompson

!=================================================================================================================
 real(kind=RKIND) function pressure_interp(target_z,nz,zf,order,extrap,surface_val,sealev_val)
!=================================================================================================================

!input arguments:
 integer,intent(in):: nz
 integer,intent(in),optional:: order
 integer,intent(in),optional:: extrap

 real(kind=RKIND),intent(in):: target_z
 real(kind=RKIND),intent(in),optional:: surface_val
 real(kind=RKIND),intent(in),optional:: sealev_val
 real(kind=RKIND),intent(in),dimension(2,nz):: zf

!local variables:
 integer:: k,lm,lp
 integer:: interp_order,extrap_type

 real(kind=RKIND):: wm,wp
 real(kind=RKIND):: slope
 real(kind=RKIND):: surface,sealevel

!-----------------------------------------------------------------------------------------------------------------

 if(present(order)) then
    interp_order = order
 else
    interp_order = 2
 endif
 if(present(extrap)) then
    extrap_type = extrap
 else
    extrap_type = 1
 endif
 if(present(surface_val)) then
    surface = surface_val
 else
    surface = 200100.0
 endif
 if(present(sealev_val)) then
    sealevel = sealev_val
 else
    sealevel = 201300.0
 endif

!--- Extrapolation required:
 if(target_z > zf(1,1)) then
    if(extrap_type == 0) then
       pressure_interp = zf(2,1)
    elseif(extrap_type == 1) then
       slope = (zf(2,2) - zf(2,1))/(zf(1,2) - zf(1,1))
       pressure_interp = zf(2,1) + slope * (target_z - zf(1,1))
    endif
    return
 endif
 if(target_z <= zf(1,nz)) then
    if(extrap_type == 0) then
       pressure_interp = zf(2,nz)
    elseif(extrap_type == 1) then
       slope = (zf(2,nz) - zf(2,nz-1)) / (zf(1,nz) - zf(1,nz-1))
       pressure_interp = zf(2,nz) + slope * (target_z - zf(1,nz))
    endif
    return
 endif

!--- No extrapolation required:
 do k =1,nz-1
    if(target_z <= zf(1,k) .and. target_z > zf(1,k+1)) then
       lm = k
       lp = k+1
       wm = (zf(1,k+1) - target_z)/(zf(1,k+1) - zf(1,k))
       wp = (target_z - zf(1,k))/(zf(1,k+1) - zf(1,k))
       exit
    endif
 enddo
 pressure_interp = wm*zf(2,lm) + wp*zf(2,lp)

 return

 end function pressure_interp

!=================================================================================================================
 subroutine read_camchem_init(fg_source,istatus)
!=================================================================================================================

!input arguments:
 character(len=*),intent(in):: fg_source

!output arguments:
 integer,intent(out):: istatus

!local variables:
 logical:: is_used
 integer:: io_status,input_unit

!-----------------------------------------------------------------------------------------------------------------

 istatus = 0

 do input_unit = 10, 100
    inquire(unit=input_unit,opened=is_used)
    if (.not. is_used) exit
 end do
 if(input_unit > 100) &
    call mpas_log_write('Error: In read_camchem_init(), couldn''t find an available Fortran unit.')
 open(unit=input_unit, file=trim(fg_source), status='old', form='unformatted', iostat=io_status)

 if (io_status > 0) istatus = 1

 end subroutine read_camchem_init

!=================================================================================================================
 end module mpas_init_atm_chems
!=================================================================================================================
