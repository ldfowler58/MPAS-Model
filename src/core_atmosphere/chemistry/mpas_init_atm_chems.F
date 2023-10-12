!=================================================================================================================
 module mpas_init_atm_chems
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
 public:: init_atm_CAMS_emissions,init_atm_camchem_toGOCART


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
 totalDelta  = sub_t_t(afterTime,beforeTime)
 beforeDelta = sub_t_t(currTime,beforeTime)
 afterDelta  = sub_t_t(afterTime,currTime)

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
 subroutine init_atm_camchem_toGOCART(configs,mesh,fg,diag,state)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: diag

!inout arguments:
 type(mpas_pool_type),intent(inout):: fg
 type(mpas_pool_type),intent(inout):: state

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine init_atm_camchem:')

 call init_hinterp_camchem(configs,mesh,fg)
 call init_vinterp_camchem(configs,mesh,fg,diag,state)

 call mpas_log_write('--- end subroutine init_atm_camchem.')
 call mpas_log_write(' ')

 end subroutine init_atm_camchem_toGOCART

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

 integer,dimension(:),pointer:: landmask,mask_array

 real(kind=RKIND):: fillval,maskval,msgval
 real(kind=RKIND):: lat,lon,x,y
 real(kind=RKIND),dimension(:),pointer:: latCell,lonCell
 real(kind=RKIND),dimension(:),pointer:: latPoints,lonPoints


!First-guess CAM-Chem aerosol mixing ratios and number concentrations:
 integer,pointer:: num_chems_fg
 integer,pointer:: camchem_start,camchem_end

 integer,pointer:: index_qbca1,index_qbca4
 integer,pointer:: index_qpoma1,index_qpoma4
 integer,pointer:: index_qdsta1,index_qdsta2,index_qdsta3
 integer,pointer:: index_qncla1,index_qncla2,index_qncla3
 integer,pointer:: index_qso4a1,index_qso4a2,index_qso4a3
 integer,pointer:: index_numa1,index_numa2,index_numa3
 integer,pointer:: index_qco,index_qh2o2,index_qoh,index_qso2,index_qno3
 integer,pointer:: index_qdms,index_qmsa

 real(kind=RKIND),dimension(:,:),pointer:: qbca1_fg,qbca4_fg
 real(kind=RKIND),dimension(:,:),pointer:: qpoma1_fg,qpoma4_fg
 real(kind=RKIND),dimension(:,:),pointer:: qdsta1_fg,qdsta2_fg,qdsta3_fg
 real(kind=RKIND),dimension(:,:),pointer:: qncla1_fg,qncla2_fg,qncla3_fg
 real(kind=RKIND),dimension(:,:),pointer:: qso4a1_fg,qso4a2_fg,qso4a3_fg
 real(kind=RKIND),dimension(:,:),pointer:: numa1_fg,numa2_fg,numa3_fg
 real(kind=RKIND),dimension(:,:),pointer:: qco_fg,qh2o2_fg,qoh_fg,qso2_fg,qno3_fg
 real(kind=RKIND),dimension(:,:),pointer:: qdms_fg,qmsa_fg

 real(kind=RKIND),dimension(:,:),pointer:: pcam_fg,dpcam_fg
 real(kind=RKIND),dimension(:,:),pointer:: destField2d
 real(kind=RKIND),dimension(:,:,:),pointer:: chems_fg

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


!--- CAMchem black carbon (a1: hydrophilic - a4: hydrophobic)
 call mpas_pool_get_dimension(fg,'index_qbca1_fg',index_qbca1)
 call mpas_pool_get_dimension(fg,'index_qbca4_fg',index_qbca4)
 call mpas_log_write('--- index_qbca1_fg  = $i',intArgs=(/index_qbca1/))
 call mpas_log_write('--- index_qbca4_fg  = $i',intArgs=(/index_qbca4/))


!--- CAMchem organic carbon (a1: hydrophilic - a4: hydrophobic)
 call mpas_pool_get_dimension(fg,'index_qpoma1_fg',index_qpoma1)
 call mpas_pool_get_dimension(fg,'index_qpoma4_fg',index_qpoma4)
 call mpas_log_write('--- index_qpoma1_fg = $i',intArgs=(/index_qpoma1/))
 call mpas_log_write('--- index_qpoma4_fg = $i',intArgs=(/index_qpoma4/))


!--- CAMchem dust:
 call mpas_pool_get_dimension(fg,'index_qdsta1_fg',index_qdsta1)
 call mpas_pool_get_dimension(fg,'index_qdsta2_fg',index_qdsta2)
 call mpas_pool_get_dimension(fg,'index_qdsta3_fg',index_qdsta3)
 call mpas_log_write('--- index_qdsta1_fg = $i',intArgs=(/index_qdsta1/))
 call mpas_log_write('--- index_qdsta2_fg = $i',intArgs=(/index_qdsta2/))
 call mpas_log_write('--- index_qdsta3_fg = $i',intArgs=(/index_qdsta3/))


!--- CAMchem seasalt:
 call mpas_pool_get_dimension(fg,'index_qncla1_fg',index_qncla1)
 call mpas_pool_get_dimension(fg,'index_qncla2_fg',index_qncla2)
 call mpas_pool_get_dimension(fg,'index_qncla3_fg',index_qncla3)
 call mpas_log_write('--- index_qncla1_fg = $i',intArgs=(/index_qncla1/))
 call mpas_log_write('--- index_qncla2_fg = $i',intArgs=(/index_qncla2/))
 call mpas_log_write('--- index_qncla3_fg = $i',intArgs=(/index_qncla3/))


!--- CAMchem sulfate:
 call mpas_pool_get_dimension(fg,'index_qso4a1_fg',index_qso4a1)
 call mpas_pool_get_dimension(fg,'index_qso4a2_fg',index_qso4a2)
 call mpas_pool_get_dimension(fg,'index_qso4a3_fg',index_qso4a3)
 call mpas_log_write('--- index_qso4a1_fg = $i',intArgs=(/index_qso4a1/))
 call mpas_log_write('--- index_qso4a2_fg = $i',intArgs=(/index_qso4a2/))
 call mpas_log_write('--- index_qso4a3_fg = $i',intArgs=(/index_qso4a3/))


!--- CAMchem number concentration:
 call mpas_pool_get_dimension(fg,'index_numa1_fg',index_numa1)
 call mpas_pool_get_dimension(fg,'index_numa2_fg',index_numa2)
 call mpas_pool_get_dimension(fg,'index_numa3_fg',index_numa3)
 call mpas_log_write('--- index_numa1_fg  = $i',intArgs=(/index_numa1/))
 call mpas_log_write('--- index_numa2_fg  = $i',intArgs=(/index_numa2/))
 call mpas_log_write('--- index_numa3_fg  = $i',intArgs=(/index_numa3/))


!--- CAMchem DMS,MSA:
 call mpas_pool_get_dimension(fg,'index_qdms_fg',index_qdms)
 call mpas_pool_get_dimension(fg,'index_qmsa_fg',index_qmsa)
 call mpas_log_write('--- index_qdms _fg  = $i',intArgs=(/index_qdms/))
 call mpas_log_write('--- index_qmsa_fg   = $i',intArgs=(/index_qmsa/))


!--- CAMchem OH,CO,H2O2,SO2,NO3:
 call mpas_pool_get_dimension(fg,'index_qco_fg'  ,index_qco  )
 call mpas_pool_get_dimension(fg,'index_qh2o2_fg',index_qh2o2)
 call mpas_pool_get_dimension(fg,'index_qoh_fg'  ,index_qoh  )
 call mpas_pool_get_dimension(fg,'index_qso2_fg' ,index_qso2 )
 call mpas_pool_get_dimension(fg,'index_qno3_fg' ,index_qno3 )
 call mpas_log_write('--- index_qco_fg    = $i',intArgs=(/index_qco/)  )
 call mpas_log_write('--- index_qh2o2_fg  = $i',intArgs=(/index_qh2o2/))
 call mpas_log_write('--- index_qoh_fg    = $i',intArgs=(/index_qoh/)  )
 call mpas_log_write('--- index_qso2_fg   = $i',intArgs=(/index_qso2/) )
 call mpas_log_write('--- index_qno3_fg   = $i',intArgs=(/index_qno3/) )
 call mpas_log_write(' ')


!--- read CAMchem data:
 call mpas_pool_get_dimension(fg,'num_chems_fg' ,num_chems_fg )
 call mpas_pool_get_dimension(fg,'camchem_start',camchem_start)
 call mpas_pool_get_dimension(fg,'camchem_end'  ,camchem_end  )
 call mpas_log_write('--- num_chems_fg    = $i',intArgs=(/num_chems_fg/) )
 call mpas_log_write('--- camchem_start   = $i',intArgs=(/camchem_start/))
 call mpas_log_write('--- camchem_end     = $i',intArgs=(/camchem_end/)  )

 call mpas_pool_get_array(fg,'chems_fg',chems_fg)
 chems_fg(:,:,:) = 0._RKIND

 qbca1_fg  => chems_fg(index_qbca1,:,:)
 qbca4_fg  => chems_fg(index_qbca4,:,:)

 qpoma1_fg => chems_fg(index_qpoma1,:,:)
 qpoma4_fg => chems_fg(index_qpoma4,:,:)

 qdsta1_fg => chems_fg(index_qdsta1,:,:)
 qdsta2_fg => chems_fg(index_qdsta2,:,:)
 qdsta3_fg => chems_fg(index_qdsta3,:,:)

 qncla1_fg => chems_fg(index_qncla1,:,:)
 qncla2_fg => chems_fg(index_qncla2,:,:)
 qncla3_fg => chems_fg(index_qncla3,:,:)

 qso4a1_fg => chems_fg(index_qso4a1,:,:)
 qso4a2_fg => chems_fg(index_qso4a2,:,:)
 qso4a3_fg => chems_fg(index_qso4a3,:,:)

 numa1_fg  => chems_fg(index_numa1,:,:)
 numa2_fg  => chems_fg(index_numa2,:,:)
 numa3_fg  => chems_fg(index_numa3,:,:)

 qdms_fg   => chems_fg(index_qdms,:,:)
 qmsa_fg   => chems_fg(index_qmsa,:,:)

 qco_fg    => chems_fg(index_qco,:,:)
 qh2o2_fg  => chems_fg(index_qh2o2,:,:)
 qoh_fg    => chems_fg(index_qoh,:,:)
 qso2_fg   => chems_fg(index_qso2,:,:)
 qno3_fg   => chems_fg(index_qno3,:,:)

 call mpas_pool_get_array(fg,'pcam_fg' ,pcam_fg )
 call mpas_pool_get_array(fg,'dpcam_fg',dpcam_fg)


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
    msgval  = -1.e30

    mask_array => landmask

    if(index(field % field, 'BCA1'  ) /= 0 .or. &
       index(field % field, 'BCA4'  ) /= 0 .or. &
       index(field % field, 'POMA1' ) /= 0 .or. &
       index(field % field, 'POMA4' ) /= 0 .or. &
       index(field % field, 'DSTA1' ) /= 0 .or. &
       index(field % field, 'DSTA2' ) /= 0 .or. &
       index(field % field, 'DSTA3' ) /= 0 .or. &
       index(field % field, 'NCLA1' ) /= 0 .or. &
       index(field % field, 'NCLA2' ) /= 0 .or. &
       index(field % field, 'NCLA3' ) /= 0 .or. &
       index(field % field, 'SO4A1' ) /= 0 .or. &
       index(field % field, 'SO4A2' ) /= 0 .or. &
       index(field % field, 'SO4A3' ) /= 0 .or. &
       index(field % field, 'NUM_A1') /= 0 .or. &
       index(field % field, 'NUM_A2') /= 0 .or. &
       index(field % field, 'NUM_A3') /= 0 .or. &
       index(field % field, 'DMS'   ) /= 0 .or. &
       index(field % field, 'MSA'   ) /= 0 .or. &
       index(field % field, 'CO'    ) /= 0 .or. &
       index(field % field, 'H2O2'  ) /= 0 .or. &
       index(field % field, 'OH'    ) /= 0 .or. &
       index(field % field, 'SO2'   ) /= 0 .or. &
       index(field % field, 'NO3'   ) /= 0 .or. &
       index(field % field, 'PCAM'  ) /= 0 .or. &
       index(field % field, 'DPRES' ) /= 0) then

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
       if(index(field%field,'BCA1') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating BCA1 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qbca1_fg
       elseif(index(field%field,'BCA4') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating BCA4 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qbca4_fg
       elseif(index(field%field,'POMA1') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating POMA1 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qpoma1_fg
       elseif(index(field%field,'POMA4') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating POMA4 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qpoma4_fg
       elseif(index(field%field,'DSTA1') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating DSTA1 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qdsta1_fg
       elseif(index(field%field,'DSTA2') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating DSTA2 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qdsta2_fg
       elseif(index(field%field,'DSTA3') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating DSTA3 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qdsta3_fg
       elseif(index(field%field,'NCLA1') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating NCLA1 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qncla1_fg
       elseif(index(field%field,'NCLA2') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating NCLA2 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qncla2_fg
       elseif(index(field%field,'NCLA3') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating NCLA3 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qncla3_fg
       elseif(index(field%field,'SO4A1') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating SO4A1 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qso4a1_fg
       elseif(index(field%field,'SO4A2') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating SO4A2 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qso4a2_fg
       elseif(index(field%field,'SO4A3') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating SO4A3 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qso4a3_fg
       elseif(index(field%field,'NUM_A1') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating NUM_A1 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => numa1_fg
       elseif(index(field%field,'NUM_A2') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating NUM_A2 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => numa2_fg
       elseif(index(field%field,'NUM_A3') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating NUM_A3 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => numa3_fg
       elseif(index(field%field,'DMS') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating DMS at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qdms_fg
       elseif(index(field%field,'MSA') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating MSA at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qmsa_fg
       elseif(index(field%field,'CO') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating CO at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qco_fg
       elseif(index(field%field,'H2O2') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating H2O2 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qh2o2_fg
       elseif(index(field%field,'OH') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating OH at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qoh_fg
       elseif(index(field%field,'SO2') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating SO2 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qso2_fg
       elseif(index(field%field,'NO3') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating NO3 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qno3_fg
       elseif(index(field%field,'PCAM') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating PCAM at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => pcam_fg
       elseif(index(field%field,'DPRES') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating DPRES at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => dpcam_fg
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


 end subroutine init_hinterp_camchem


!=================================================================================================================
 subroutine init_vinterp_camchem(configs,mesh,fg,diag,state)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: fg
 type(mpas_pool_type),intent(in):: diag

!inout arguments:
 type(mpas_pool_type),intent(inout):: state

!local variables and arrays:
 integer:: iCell,k,kk,n,nn
 integer,pointer:: nCells,nCAMLevels,nVertLevels

 integer,pointer:: num_chems_fg,camchem_start,camchem_end
 integer,pointer:: index_qbca1,index_qbca4
 integer,pointer:: index_qpoma1,index_qpoma4
 integer,pointer:: index_qdsta1,index_qdsta2,index_qdsta3
 integer,pointer:: index_qncla1,index_qncla2,index_qncla3
 integer,pointer:: index_qso4a1,index_qso4a2,index_qso4a3
 integer,pointer:: index_qdms_fg,index_qmsa_fg,index_qco_fg,index_qh2o2_fg,index_qoh_fg, &
                   index_qso2_fg,index_qno3_fg

 integer,pointer:: num_chems,gocart_start,gocart_end
 integer,pointer:: index_qbc1,index_qbc2
 integer,pointer:: index_qoc1,index_qoc2
 integer,pointer:: index_qdust1,index_qdust2,index_qdust3,index_qdust4,index_qdust5
 integer,pointer:: index_qseas1,index_qseas2,index_qseas3,index_qseas4
 integer,pointer:: index_qsulf
 integer,pointer:: index_qdms,index_qmsa,index_qco,index_qh2o2,index_qoh,index_qso2,index_qno3

 real(kind=RKIND),dimension(:,:),pointer:: pcam_fg,pressure
 real(kind=RKIND),dimension(:,:,:),pointer:: chems_fg
 real(kind=RKIND),dimension(:,:,:),pointer:: chems

 real(kind=RKIND):: target_p
 real(kind=RKIND),dimension(:,:),allocatable:: sorted_arr

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine init_vinterp_camchem:')

 call mpas_pool_get_dimension(mesh,'nCells',nCells)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
 call mpas_pool_get_dimension(mesh,'nCAMLevels' ,nCAMLevels )

 call mpas_pool_get_array(diag,'pressure',pressure)

 call mpas_pool_get_array(fg,'pcam_fg' ,pcam_fg )
 call mpas_pool_get_array(fg,'chems_fg',chems_fg) 
 call mpas_pool_get_array(state,'chems',chems)
 chems(:,:,:) = 0._RKIND

 call mpas_pool_get_dimension(fg,'num_chems_fg' ,num_chems_fg )
 call mpas_pool_get_dimension(fg,'camchem_start',camchem_start)
 call mpas_pool_get_dimension(fg,'camchem_end'  ,camchem_end  )
 call mpas_log_write('--- num_chems_fg  = $i',intArgs=(/num_chems_fg/) )
 call mpas_log_write('--- camchem_start = $i',intArgs=(/camchem_start/))
 call mpas_log_write('--- camchem_end   = $i',intArgs=(/camchem_end/)  )
 call mpas_log_write(' ')


 call mpas_pool_get_dimension(state,'num_chems' ,num_chems )
 call mpas_pool_get_dimension(state,'gocart_start',gocart_start)
 call mpas_pool_get_dimension(state,'gocart_end'  ,gocart_end  )
 call mpas_log_write('--- num_chems     = $i',intArgs=(/num_chems/)   )
 call mpas_log_write('--- gocart_start  = $i',intArgs=(/gocart_start/))
 call mpas_log_write('--- gocart_end    = $i',intArgs=(/gocart_end/)  )
 call mpas_log_write(' ')

 if(.not.allocated(sorted_arr)) allocate(sorted_arr(2,nCAMLevels))

!--- match CAMchem primary black (organic) carbon to Gocart hydrophobic black (organic) carbon, and match
!    CAMchem black (organic) carbon from accumulation mode to hydrophilic black (organic) carbon:
!    QBCA4_FG  -> QBC1; QBCA1_FG  -> QBC2.
!    QPOMA4_FG -> QOC2; QPOMA1_FG -> QOC2.

 call mpas_pool_get_dimension(fg,'index_qbca1_fg',index_qbca1)
 call mpas_pool_get_dimension(fg,'index_qbca4_fg',index_qbca4)
 call mpas_pool_get_dimension(state,'index_qbc1',index_qbc1)
 call mpas_pool_get_dimension(state,'index_qbc2',index_qbc2)
 call mpas_log_write('--- index_qbca1   = $i',intArgs=(/index_qbca1/))
 call mpas_log_write('--- index_qbca4   = $i',intArgs=(/index_qbca4/))
 call mpas_log_write('--- index_qbc1    = $i',intArgs=(/index_qbc1/) )
 call mpas_log_write('--- index_qbc2    = $i',intArgs=(/index_qbc2/) )
 call mpas_log_write(' ')

 do n = camchem_start,camchem_end
    do nn = gocart_start,gocart_end
       if(n == index_qbca4 .and. nn == index_qbc1) then
          call mpas_log_write('--- BLACK CARBON - HYDROPHOBIC: n = $i  nn = $i',intArgs=(/n,nn/))
          do iCell = 1,nCells
             sorted_arr(1,1:nCAMLevels) = 0._RKIND
             sorted_arr(2,1:nCAMLevels) = 0._RKIND
             do k = 1,nCAMLevels
                sorted_arr(1,k) = pcam_fg(k,iCell)
                sorted_arr(2,k) = chems_fg(n,k,iCell)
             enddo
             do k = nVertLevels,1,-1
                target_p = pressure(k,iCell)
                chems(nn,k,iCell) = pressure_interp(iCell,k,target_p,nCAMLevels,sorted_arr(:,1:nCAMLevels))
                if(target_p.gt.sorted_arr(1,1)) chems(nn,k,iCell) = chems(nn,k+1,iCell)
             enddo
          enddo
       elseif(n == index_qbca1 .and. nn == index_qbc2) then
          call mpas_log_write('--- BLACK CARBON - HYDROPHILIC: n = $i  nn = $i',intArgs=(/n,nn/))
          do iCell = 1,nCells
             sorted_arr(1,1:nCAMLevels) = 0._RKIND
             sorted_arr(2,1:nCAMLevels) = 0._RKIND
             do k = 1,nCAMLevels
                sorted_arr(1,k) = pcam_fg(k,iCell)
                sorted_arr(2,k) = chems_fg(n,k,iCell)
             enddo
             do k = nVertLevels,1,-1
                target_p = pressure(k,iCell)
                chems(nn,k,iCell) = pressure_interp(iCell,k,target_p,nCAMLevels,sorted_arr(:,1:nCAMLevels))
                if(target_p.gt.sorted_arr(1,1)) chems(nn,k,iCell) = chems(nn,k+1,iCell)
             enddo
          enddo
       endif
    enddo
 enddo


 call mpas_pool_get_dimension(fg,'index_qpoma1_fg',index_qpoma1)
 call mpas_pool_get_dimension(fg,'index_qpoma4_fg',index_qpoma4)
 call mpas_pool_get_dimension(state,'index_qoc1',index_qoc1)
 call mpas_pool_get_dimension(state,'index_qoc2',index_qoc2)
 call mpas_log_write('--- index_qpoma1   = $i',intArgs=(/index_qpoma1/))
 call mpas_log_write('--- index_qpoma4   = $i',intArgs=(/index_qpoma4/))
 call mpas_log_write('--- index_qoc1     = $i',intArgs=(/index_qoc1/)  )
 call mpas_log_write('--- index_qoc2     = $i',intArgs=(/index_qoc2/)  )
 call mpas_log_write(' ')

 do n = camchem_start,camchem_end
    do nn = gocart_start,gocart_end
       if(n == index_qpoma4 .and. nn == index_qoc1) then
          call mpas_log_write('--- ORGANIC CARBON - HYDROPHOBIC: n = $i  nn = $i',intArgs=(/n,nn/))
          do iCell = 1,nCells
             sorted_arr(1,1:nCAMLevels) = 0._RKIND
             sorted_arr(2,1:nCAMLevels) = 0._RKIND
             do k = 1,nCAMLevels
                sorted_arr(1,k) = pcam_fg(k,iCell)
                sorted_arr(2,k) = chems_fg(n,k,iCell)
             enddo
             do k = nVertLevels,1,-1
                target_p = pressure(k,iCell)
                chems(nn,k,iCell) = pressure_interp(iCell,k,target_p,nCAMLevels,sorted_arr(:,1:nCAMLevels))
                if(target_p.gt.sorted_arr(1,1)) chems(nn,k,iCell) = chems(nn,k+1,iCell)
             enddo
          enddo
       elseif(n == index_qpoma1 .and. nn == index_qoc2) then
          call mpas_log_write('--- ORGANIC CARBON - HYDROPHILIC: n = $i  nn = $i',intArgs=(/n,nn/))
          do iCell = 1,nCells
             sorted_arr(1,1:nCAMLevels) = 0._RKIND
             sorted_arr(2,1:nCAMLevels) = 0._RKIND
             do k = 1,nCAMLevels
                sorted_arr(1,k) = pcam_fg(k,iCell)
                sorted_arr(2,k) = chems_fg(n,k,iCell)
             enddo
             do k = nVertLevels,1,-1
                target_p = pressure(k,iCell)
                chems(nn,k,iCell) = pressure_interp(iCell,k,target_p,nCAMLevels,sorted_arr(:,1:nCAMLevels))
                if(target_p.gt.sorted_arr(1,1)) chems(nn,k,iCell) = chems(nn,k+1,iCell)
             enddo
          enddo
       endif
    enddo
 enddo


!--- match CAMchem dust from the accumulation mode, aitken mode, and coarse mode to Gocart dust in the
!    five dust bins:
 call mpas_pool_get_dimension(fg,'index_qdsta1_fg',index_qdsta1)
 call mpas_pool_get_dimension(fg,'index_qdsta2_fg',index_qdsta2)
 call mpas_pool_get_dimension(fg,'index_qdsta3_fg',index_qdsta3)
 call mpas_pool_get_dimension(state,'index_qdust1',index_qdust1)
 call mpas_pool_get_dimension(state,'index_qdust2',index_qdust2)
 call mpas_pool_get_dimension(state,'index_qdust3',index_qdust3)
 call mpas_pool_get_dimension(state,'index_qdust4',index_qdust4)
 call mpas_pool_get_dimension(state,'index_qdust5',index_qdust5)
 call mpas_log_write(' ')
 call mpas_log_write('--- index_qdsta1   = $i',intArgs=(/index_qdsta1/))
 call mpas_log_write('--- index_qdsta2   = $i',intArgs=(/index_qdsta2/))
 call mpas_log_write('--- index_qdsta3   = $i',intArgs=(/index_qdsta3/))
 call mpas_log_write('--- index_qdust1   = $i',intArgs=(/index_qdust1/))
 call mpas_log_write('--- index_qdust2   = $i',intArgs=(/index_qdust2/))
 call mpas_log_write('--- index_qdust3   = $i',intArgs=(/index_qdust3/))
 call mpas_log_write('--- index_qdust4   = $i',intArgs=(/index_qdust4/))
 call mpas_log_write('--- index_qdust5   = $i',intArgs=(/index_qdust5/))
 call mpas_log_write(' ')

 do nn = gocart_start,gocart_end
    if(nn == index_qdust1) then
       do iCell = 1,nCells
          sorted_arr(1,1:nCAMLevels) = 0._RKIND
          sorted_arr(2,1:nCAMLevels) = 0._RKIND
          do k = 1,nCAMLevels
             sorted_arr(1,k) = pcam_fg(k,iCell)
             sorted_arr(2,k) = 0.02*chems_fg(index_qdsta3,k,iCell)
          enddo
          do k = nVertLevels,1,-1
             target_p = pressure(k,iCell)
             chems(nn,k,iCell) = pressure_interp(iCell,k,target_p,nCAMLevels,sorted_arr(:,1:nCAMLevels))
             if(target_p.gt.sorted_arr(1,1)) chems(nn,k,iCell) = chems(nn,k+1,iCell)
          enddo
       enddo
    elseif(nn == index_qdust2) then
       do iCell = 1,nCells
          sorted_arr(1,1:nCAMLevels) = 0._RKIND
          sorted_arr(2,1:nCAMLevels) = 0._RKIND
          do k = 1,nCAMLevels
             sorted_arr(1,k) = pcam_fg(k,iCell)
             sorted_arr(2,k) = 0.93*chems_fg(index_qdsta3,k,iCell)
          enddo
          do k = nVertLevels,1,-1
             target_p = pressure(k,iCell)
             chems(nn,k,iCell) = pressure_interp(iCell,k,target_p,nCAMLevels,sorted_arr(:,1:nCAMLevels))
             if(target_p.gt.sorted_arr(1,1)) chems(nn,k,iCell) = chems(nn,k+1,iCell)
          enddo
       enddo
    elseif(nn == index_qdust3) then
       do iCell = 1,nCells
          sorted_arr(1,1:nCAMLevels) = 0._RKIND
          sorted_arr(2,1:nCAMLevels) = 0._RKIND
          do k = 1,nCAMLevels
             sorted_arr(1,k) = pcam_fg(k,iCell)
             sorted_arr(2,k) = 0.05*chems_fg(index_qdsta3,k,iCell)
          enddo
          do k = nVertLevels,1,-1
             target_p = pressure(k,iCell)
             chems(nn,k,iCell) = pressure_interp(iCell,k,target_p,nCAMLevels,sorted_arr(:,1:nCAMLevels))
             if(target_p.gt.sorted_arr(1,1)) chems(nn,k,iCell) = chems(nn,k+1,iCell)
          enddo
       enddo
    elseif(nn == index_qdust4 .or. nn == index_qdust5) then
       do iCell = 1,nCells
          do k = nVertLevels,1,-1
             chems(nn,k,iCell) = 0._RKIND
          enddo
       enddo
    endif
 enddo


!--- match CAMchem seasalt from the accumulation mode, aitken mode, and coarse mode to Gocart dust in the
!    four seasalt bins:
 call mpas_pool_get_dimension(fg,'index_qncla1_fg',index_qncla1)
 call mpas_pool_get_dimension(fg,'index_qncla2_fg',index_qncla2)
 call mpas_pool_get_dimension(fg,'index_qncla3_fg',index_qncla3)
 call mpas_pool_get_dimension(state,'index_qseas1',index_qseas1)
 call mpas_pool_get_dimension(state,'index_qseas2',index_qseas2)
 call mpas_pool_get_dimension(state,'index_qseas3',index_qseas3)
 call mpas_pool_get_dimension(state,'index_qseas4',index_qseas4)
 call mpas_log_write('--- index_qncla1   = $i',intArgs=(/index_qncla1/))
 call mpas_log_write('--- index_qncla2   = $i',intArgs=(/index_qncla2/))
 call mpas_log_write('--- index_qncla3   = $i',intArgs=(/index_qncla3/))
 call mpas_log_write('--- index_qseas1   = $i',intArgs=(/index_qseas1/))
 call mpas_log_write('--- index_qseas2   = $i',intArgs=(/index_qseas2/))
 call mpas_log_write('--- index_qseas3   = $i',intArgs=(/index_qseas3/))
 call mpas_log_write('--- index_qseas4   = $i',intArgs=(/index_qseas4/))
 call mpas_log_write(' ')

 do nn = gocart_start,gocart_end
    if(nn == index_qseas1) then
       do iCell = 1,nCells
          sorted_arr(1,1:nCAMLevels) = 0._RKIND
          sorted_arr(2,1:nCAMLevels) = 0._RKIND
          do k = 1,nCAMLevels
             sorted_arr(1,k) = pcam_fg(k,iCell)
             sorted_arr(2,k) = chems_fg(index_qncla1,k,iCell) + chems_fg(index_qncla2,k,iCell)
          enddo
          do k = nVertLevels,1,-1
             target_p = pressure(k,iCell)
             chems(nn,k,iCell) = pressure_interp(iCell,k,target_p,nCAMLevels,sorted_arr(:,1:nCAMLevels))
             if(target_p.gt.sorted_arr(1,1)) chems(nn,k,iCell) = chems(nn,k+1,iCell)
          enddo
       enddo
    elseif(nn == index_qseas2) then
       do iCell = 1,nCells
          sorted_arr(1,1:nCAMLevels) = 0._RKIND
          sorted_arr(2,1:nCAMLevels) = 0._RKIND
          do k = 1,nCAMLevels
             sorted_arr(1,k) = pcam_fg(k,iCell)
             sorted_arr(2,k) = 0.5*chems_fg(index_qncla3,k,iCell)
          enddo
          do k = nVertLevels,1,-1
             target_p = pressure(k,iCell)
             chems(nn,k,iCell) = pressure_interp(iCell,k,target_p,nCAMLevels,sorted_arr(:,1:nCAMLevels))
             if(target_p.gt.sorted_arr(1,1)) chems(nn,k,iCell) = chems(nn,k+1,iCell)
          enddo
       enddo
    elseif(nn == index_qseas3) then
       do iCell = 1,nCells
          sorted_arr(1,1:nCAMLevels) = 0._RKIND
          sorted_arr(2,1:nCAMLevels) = 0._RKIND
          do k = 1,nCAMLevels
             sorted_arr(1,k) = pcam_fg(k,iCell)
             sorted_arr(2,k) = 0.5*chems_fg(index_qncla3,k,iCell)
          enddo
          do k = nVertLevels,1,-1
             target_p = pressure(k,iCell)
             chems(nn,k,iCell) = pressure_interp(iCell,k,target_p,nCAMLevels,sorted_arr(:,1:nCAMLevels))
             if(target_p.gt.sorted_arr(1,1)) chems(nn,k,iCell) = chems(nn,k+1,iCell)
          enddo
       enddo
    elseif(nn == index_qseas4) then
       do iCell = 1,nCells
          do k = nVertLevels,1,-1
             chems(nn,k,iCell) = 0._RKIND
          enddo
       enddo
    endif
 enddo


!SULFATES:
 call mpas_pool_get_dimension(fg,'index_qso4a1_fg',index_qso4a1)
 call mpas_pool_get_dimension(fg,'index_qso4a2_fg',index_qso4a2)
 call mpas_pool_get_dimension(fg,'index_qso4a3_fg',index_qso4a3)
 call mpas_pool_get_dimension(state,'index_qsulf',index_qsulf)
 call mpas_log_write('--- index_qso4a1   = $i',intArgs=(/index_qso4a1/))
 call mpas_log_write('--- index_qso4a2   = $i',intArgs=(/index_qso4a2/))
 call mpas_log_write('--- index_qso4a3   = $i',intArgs=(/index_qso4a3/))
 call mpas_log_write('--- index_qsulf    = $i',intArgs=(/index_qsulf/))
 call mpas_log_write(' ')

 do nn = gocart_start,gocart_end
    if(nn == index_qsulf) then
       do iCell = 1,nCells
          sorted_arr(1,1:nCAMLevels) = 0._RKIND
          sorted_arr(2,1:nCAMLevels) = 0._RKIND
          do k = 1,nCAMLevels
             sorted_arr(1,k) = pcam_fg(k,iCell)
             sorted_arr(2,k) = chems_fg(index_qso4a1,k,iCell) + chems_fg(index_qso4a2,k,iCell)
          enddo
          do k = nVertLevels,1,-1
             target_p = pressure(k,iCell)
             chems(nn,k,iCell) = pressure_interp(iCell,k,target_p,nCAMLevels,sorted_arr(:,1:nCAMLevels))
             if(target_p.gt.sorted_arr(1,1)) chems(nn,k,iCell) = chems(nn,k+1,iCell)
          enddo
       enddo
    endif
 enddo


!DMS,MSA,QCO,QH2O,QOH,QSO2,NO3:
 call mpas_pool_get_dimension(fg,'index_qdms_fg' ,index_qdms_fg )
 call mpas_pool_get_dimension(fg,'index_qmsa_fg' ,index_qmsa_fg )
 call mpas_pool_get_dimension(fg,'index_qco_fg'  ,index_qco_fg  )
 call mpas_pool_get_dimension(fg,'index_qh2o2_fg',index_qh2o2_fg)
 call mpas_pool_get_dimension(fg,'index_qoh_fg'  ,index_qoh_fg  )
 call mpas_pool_get_dimension(fg,'index_qso2_fg' ,index_qso2_fg )
 call mpas_pool_get_dimension(fg,'index_qno3_fg' ,index_qno3_fg )
 call mpas_pool_get_dimension(state,'index_qdms' ,index_qdms )
 call mpas_pool_get_dimension(state,'index_qmsa' ,index_qmsa )
 call mpas_pool_get_dimension(state,'index_qco'  ,index_qco  )
 call mpas_pool_get_dimension(state,'index_qh2o2',index_qh2o2)
 call mpas_pool_get_dimension(state,'index_qoh'  ,index_qoh  )
 call mpas_pool_get_dimension(state,'index_qso2' ,index_qso2 )
 call mpas_pool_get_dimension(state,'index_qno3' ,index_qno3 )
 call mpas_log_write('--- index_qdms_fg  = $i',intArgs=(/index_qdms_fg/) )
 call mpas_log_write('--- index_qmsa_fg  = $i',intArgs=(/index_qmsa_fg/) )
 call mpas_log_write('--- index_qco_fg   = $i',intArgs=(/index_qco_fg/)  )
 call mpas_log_write('--- index_qh2o2_fg = $i',intArgs=(/index_qh2o2_fg/))
 call mpas_log_write('--- index_qoh_fg   = $i',intArgs=(/index_qoh_fg/)  )
 call mpas_log_write('--- index_qso2_fg  = $i',intArgs=(/index_qso2_fg/) )
 call mpas_log_write('--- index_qno3_fg  = $i',intArgs=(/index_qno3_fg/) )
 call mpas_log_write('--- index_qdms     = $i',intArgs=(/index_qdms/)    )
 call mpas_log_write('--- index_qmsa     = $i',intArgs=(/index_qmsa/)    )
 call mpas_log_write('--- index_qco      = $i',intArgs=(/index_qco/)     )
 call mpas_log_write('--- index_qh2o2    = $i',intArgs=(/index_qh2o2/)   )
 call mpas_log_write('--- index_qoh      = $i',intArgs=(/index_qoh/)     )
 call mpas_log_write('--- index_qso2     = $i',intArgs=(/index_qso2/)    )
 call mpas_log_write('--- index_qno3     = $i',intArgs=(/index_qno3/)    )
 call mpas_log_write(' ')

 do nn = gocart_start,gocart_end
    if(nn == index_qdms) then
       do iCell = 1,nCells
          sorted_arr(1,1:nCAMLevels) = 0._RKIND
          sorted_arr(2,1:nCAMLevels) = 0._RKIND
          do k = 1,nCAMLevels
             sorted_arr(1,k) = pcam_fg(k,iCell)
             sorted_arr(2,k) = chems_fg(index_qdms_fg,k,iCell)
          enddo
          do k = nVertLevels,1,-1
             target_p = pressure(k,iCell)
             chems(nn,k,iCell) = pressure_interp(iCell,k,target_p,nCAMLevels,sorted_arr(:,1:nCAMLevels))
             if(target_p.gt.sorted_arr(1,1)) chems(nn,k,iCell) = chems(nn,k+1,iCell)
          enddo
       enddo
    elseif(nn == index_qmsa) then
       do iCell = 1,nCells
          sorted_arr(1,1:nCAMLevels) = 0._RKIND
          sorted_arr(2,1:nCAMLevels) = 0._RKIND
          do k = 1,nCAMLevels
             sorted_arr(1,k) = pcam_fg(k,iCell)
             sorted_arr(2,k) = chems_fg(index_qmsa_fg,k,iCell)
          enddo
          do k = nVertLevels,1,-1
             target_p = pressure(k,iCell)
             chems(nn,k,iCell) = pressure_interp(iCell,k,target_p,nCAMLevels,sorted_arr(:,1:nCAMLevels))
             if(target_p.gt.sorted_arr(1,1)) chems(nn,k,iCell) = chems(nn,k+1,iCell)
          enddo
       enddo
    elseif(nn == index_qco) then
       do iCell = 1,nCells
          sorted_arr(1,1:nCAMLevels) = 0._RKIND
          sorted_arr(2,1:nCAMLevels) = 0._RKIND
          do k = 1,nCAMLevels
             sorted_arr(1,k) = pcam_fg(k,iCell)
             sorted_arr(2,k) = chems_fg(index_qco_fg,k,iCell)
          enddo
          do k = nVertLevels,1,-1
             target_p = pressure(k,iCell)
             chems(nn,k,iCell) = pressure_interp(iCell,k,target_p,nCAMLevels,sorted_arr(:,1:nCAMLevels))
             if(target_p.gt.sorted_arr(1,1)) chems(nn,k,iCell) = chems(nn,k+1,iCell)
          enddo
       enddo
    elseif(nn == index_qh2o2) then
       do iCell = 1,nCells
          sorted_arr(1,1:nCAMLevels) = 0._RKIND
          sorted_arr(2,1:nCAMLevels) = 0._RKIND
          do k = 1,nCAMLevels
             sorted_arr(1,k) = pcam_fg(k,iCell)
             sorted_arr(2,k) = chems_fg(index_qh2o2_fg,k,iCell)
          enddo
          do k = nVertLevels,1,-1
             target_p = pressure(k,iCell)
             chems(nn,k,iCell) = pressure_interp(iCell,k,target_p,nCAMLevels,sorted_arr(:,1:nCAMLevels))
             if(target_p.gt.sorted_arr(1,1)) chems(nn,k,iCell) = chems(nn,k+1,iCell)
          enddo
       enddo
    elseif(nn == index_qoh) then
       do iCell = 1,nCells
          sorted_arr(1,1:nCAMLevels) = 0._RKIND
          sorted_arr(2,1:nCAMLevels) = 0._RKIND
          do k = 1,nCAMLevels
             sorted_arr(1,k) = pcam_fg(k,iCell)
             sorted_arr(2,k) = chems_fg(index_qoh_fg,k,iCell)
          enddo
          do k = nVertLevels,1,-1
             target_p = pressure(k,iCell)
             chems(nn,k,iCell) = pressure_interp(iCell,k,target_p,nCAMLevels,sorted_arr(:,1:nCAMLevels))
             if(target_p.gt.sorted_arr(1,1)) chems(nn,k,iCell) = chems(nn,k+1,iCell)
          enddo
       enddo
    elseif(nn == index_qso2) then
       do iCell = 1,nCells
          sorted_arr(1,1:nCAMLevels) = 0._RKIND
          sorted_arr(2,1:nCAMLevels) = 0._RKIND
          do k = 1,nCAMLevels
             sorted_arr(1,k) = pcam_fg(k,iCell)
             sorted_arr(2,k) = chems_fg(index_qso2_fg,k,iCell)
          enddo
          do k = nVertLevels,1,-1
             target_p = pressure(k,iCell)
             chems(nn,k,iCell) = pressure_interp(iCell,k,target_p,nCAMLevels,sorted_arr(:,1:nCAMLevels))
             if(target_p.gt.sorted_arr(1,1)) chems(nn,k,iCell) = chems(nn,k+1,iCell)
          enddo
       enddo
    elseif(nn == index_qno3) then
       do iCell = 1,nCells
          sorted_arr(1,1:nCAMLevels) = 0._RKIND
          sorted_arr(2,1:nCAMLevels) = 0._RKIND
          do k = 1,nCAMLevels
             sorted_arr(1,k) = pcam_fg(k,iCell)
             sorted_arr(2,k) = chems_fg(index_qno3_fg,k,iCell)
          enddo
          do k = nVertLevels,1,-1
             target_p = pressure(k,iCell)
             chems(nn,k,iCell) = pressure_interp(iCell,k,target_p,nCAMLevels,sorted_arr(:,1:nCAMLevels))
             if(target_p.gt.sorted_arr(1,1)) chems(nn,k,iCell) = chems(nn,k+1,iCell)
          enddo
       enddo
    endif
 enddo

 if(allocated(sorted_arr)) deallocate(sorted_arr)

 call mpas_log_write('--- end subroutine init_vinterp_camchem.')

 end subroutine init_vinterp_camchem

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
