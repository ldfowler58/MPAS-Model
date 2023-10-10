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
 public:: init_atm_CAMS_emissions,init_atm_camchem


!initialization of CAMS emissions and interpolation of CAM-Chem input data to the MPAS mesh.
!Laura D. Fowler (laura@ucar.edu) / 2022-02-08.

!add-ons and modifications to sourcecode:
!----------------------------------------
! * in subroutine init_vinterp_camchem, changed the initialization of qco from a diagnostic variable to a
!   prognostic state variable so that carbon monoxide can be used as a passive tracer.
!   Laura D. Fowler (laura@ucar.edu) / 2022-06-23.


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
 subroutine init_atm_camchem(configs,mesh,fg,diag,diag_physics,chems_input,state)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: diag

!inout arguments:
 type(mpas_pool_type),intent(inout):: fg
 type(mpas_pool_type),intent(inout):: chems_input
 type(mpas_pool_type),intent(inout):: state
 type(mpas_pool_type),intent(inout):: diag_physics

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine init_atm_camchem:')

 call init_hinterp_camchem(configs,mesh,fg)
 call init_vinterp_camchem(configs,mesh,fg,diag,diag_physics,chems_input,state)
!call init_camchem_tothompson(mesh,chems_input,state)

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
 integer,dimension(:),pointer:: landmask,mask_array

 real(kind=RKIND):: fillval,maskval,msgval
 real(kind=RKIND):: lat,lon,x,y
 real(kind=RKIND),dimension(:),pointer:: latCell,lonCell
 real(kind=RKIND),dimension(:),pointer:: latPoints,lonPoints

!CAM-Chem aerosol mixing ratios and number concentrations:
 integer,pointer:: index_qco
 integer,pointer:: index_numa1,index_qbca1,index_qdsta1,index_qncla1,index_qpoma1,index_qso4a1
 integer,pointer:: index_numa2,index_qdsta2,index_qncla2,index_qso4a2
 integer,pointer:: index_numa3,index_qdsta3,index_qncla3,index_qso4a3
 integer,pointer:: index_qbca4,index_qpoma4
 integer,pointer:: index_qsoa1_a1,index_qsoa2_a1,index_qsoa3_a1,index_qsoa4_a1,index_qsoa5_a1
 integer,pointer:: index_qsoa1_a2,index_qsoa2_a2,index_qsoa3_a2,index_qsoa4_a2,index_qsoa5_a2

 real(kind=RKIND),dimension(:,:),pointer:: qco_fg
 real(kind=RKIND),dimension(:,:),pointer:: qbca1_fg,qdsta1_fg,qncla1_fg,qpoma1_fg,qso4a1_fg
 real(kind=RKIND),dimension(:,:),pointer:: qdsta2_fg,qncla2_fg,qso4a2_fg
 real(kind=RKIND),dimension(:,:),pointer:: qdsta3_fg,qncla3_fg,qso4a3_fg
 real(kind=RKIND),dimension(:,:),pointer:: qbca4_fg,qpoma4_fg
 real(kind=RKIND),dimension(:,:),pointer:: numa1_fg,numa2_fg,numa3_fg

 real(kind=RKIND),dimension(:,:),pointer:: qsoa1_a1_fg,qsoa2_a1_fg,qsoa3_a1_fg,qsoa4_a1_fg,qsoa5_a1_fg
 real(kind=RKIND),dimension(:,:),pointer:: qsoa1_a2_fg,qsoa2_a2_fg,qsoa3_a2_fg,qsoa4_a2_fg,qsoa5_a2_fg

!GOCART aerosol mixing ratios:
 integer,pointer:: index_qbc1,index_qbc2,index_qoc1,index_qoc2
 integer,pointer:: index_qdust1,index_qdust2,index_qdust3,index_qdust4,index_qdust5
 integer,pointer:: index_qseas1,index_qseas2,index_qseas3,index_qseas4
 real(kind=RKIND),dimension(:,:),pointer:: qbc1_fg,qbc2_fg,qoc1_fg,qoc2_fg
 real(kind=RKIND),dimension(:,:),pointer:: qdust1_fg,qdust2_fg,qdust3_fg,qdust4_fg,qdust5_fg
 real(kind=RKIND),dimension(:,:),pointer:: qseas1_fg,qseas2_fg,qseas3_fg,qseas4_fg

!Ice-friendly and water-friendly Thompson aerosols:
 integer,pointer:: index_nifa,index_nwfa
 real(kind=RKIND),dimension(:,:),pointer:: nifa_fg,nwfa_fg

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

 call mpas_pool_get_dimension(fg,'num_chems_fg'   ,num_chems_fg  )
 call mpas_pool_get_dimension(fg,'gocart_fg_start',gocart_fg_start)
 call mpas_pool_get_dimension(fg,'gocart_fg_end'  ,gocart_fg_end  )
 call mpas_log_write('--- num_chems_fg    = $i',intArgs=(/num_chems_fg/)   )
 call mpas_log_write('--- gocart_fg_start = $i',intArgs=(/gocart_fg_start/))
 call mpas_log_write('--- gocart_fg_end   = $i',intArgs=(/gocart_fg_end/)  )
 call mpas_log_write(' ')

 call mpas_pool_get_dimension(fg,'index_qco_fg'   ,index_qco   )

 call mpas_pool_get_dimension(fg,'index_qbca1_fg' ,index_qbca1 )
 call mpas_pool_get_dimension(fg,'index_qdsta1_fg',index_qdsta1)
 call mpas_pool_get_dimension(fg,'index_qncla1_fg',index_qncla1)
 call mpas_pool_get_dimension(fg,'index_qpoma1_fg',index_qpoma1)
 call mpas_pool_get_dimension(fg,'index_qso4a1_fg',index_qso4a1)
 call mpas_pool_get_dimension(fg,'index_qsoa1_a1_fg',index_qsoa1_a1)
 call mpas_pool_get_dimension(fg,'index_qsoa2_a1_fg',index_qsoa2_a1)
 call mpas_pool_get_dimension(fg,'index_qsoa3_a1_fg',index_qsoa3_a1)
 call mpas_pool_get_dimension(fg,'index_qsoa4_a1_fg',index_qsoa4_a1)
 call mpas_pool_get_dimension(fg,'index_qsoa5_a1_fg',index_qsoa5_a1)
 call mpas_pool_get_dimension(fg,'index_numa1_fg' ,index_numa1 )

 call mpas_pool_get_dimension(fg,'index_qdsta2_fg',index_qdsta2)
 call mpas_pool_get_dimension(fg,'index_qncla2_fg',index_qncla2)
 call mpas_pool_get_dimension(fg,'index_qso4a2_fg',index_qso4a2)
 call mpas_pool_get_dimension(fg,'index_qsoa1_a2_fg',index_qsoa1_a2)
 call mpas_pool_get_dimension(fg,'index_qsoa2_a2_fg',index_qsoa2_a2)
 call mpas_pool_get_dimension(fg,'index_qsoa3_a2_fg',index_qsoa3_a2)
 call mpas_pool_get_dimension(fg,'index_qsoa4_a2_fg',index_qsoa4_a2)
 call mpas_pool_get_dimension(fg,'index_qsoa5_a2_fg',index_qsoa5_a2)
 call mpas_pool_get_dimension(fg,'index_numa2_fg' ,index_numa2 )

 call mpas_pool_get_dimension(fg,'index_qdsta3_fg',index_qdsta3)
 call mpas_pool_get_dimension(fg,'index_qncla3_fg',index_qncla3)
 call mpas_pool_get_dimension(fg,'index_qso4a3_fg',index_qso4a3)
 call mpas_pool_get_dimension(fg,'index_numa3_fg' ,index_numa3 )

 call mpas_pool_get_dimension(fg,'index_qbca4_fg' ,index_qbca4 )
 call mpas_pool_get_dimension(fg,'index_qpoma4_fg',index_qpoma4)

 call mpas_pool_get_dimension(fg,'index_qbc1_fg'  ,index_qbc1  )
 call mpas_pool_get_dimension(fg,'index_qbc2_fg'  ,index_qbc2  )
 call mpas_pool_get_dimension(fg,'index_qoc1_fg'  ,index_qoc1  )
 call mpas_pool_get_dimension(fg,'index_qoc2_fg'  ,index_qoc2  )
 call mpas_pool_get_dimension(fg,'index_qdust1_fg',index_qdust1)
 call mpas_pool_get_dimension(fg,'index_qdust2_fg',index_qdust2)
 call mpas_pool_get_dimension(fg,'index_qdust3_fg',index_qdust3)
 call mpas_pool_get_dimension(fg,'index_qdust4_fg',index_qdust4)
 call mpas_pool_get_dimension(fg,'index_qdust5_fg',index_qdust5)
 call mpas_pool_get_dimension(fg,'index_qseas1_fg',index_qseas1)
 call mpas_pool_get_dimension(fg,'index_qseas2_fg',index_qseas2)
 call mpas_pool_get_dimension(fg,'index_qseas3_fg',index_qseas3)
 call mpas_pool_get_dimension(fg,'index_qseas4_fg',index_qseas4)

 call mpas_pool_get_dimension(fg,'index_nifa_fg',index_nifa)
 call mpas_pool_get_dimension(fg,'index_nwfa_fg',index_nwfa)

 call mpas_log_write('--- index_qco_fg    = $i',intArgs=(/index_qco/)   )

 call mpas_log_write('--- accumulation mode:')
 call mpas_log_write('--- index_qbca1_fg  = $i',intArgs=(/index_qbca1/) )
 call mpas_log_write('--- index_qdsta1_fg = $i',intArgs=(/index_qdsta1/))
 call mpas_log_write('--- index_qncla1_fg = $i',intArgs=(/index_qncla1/))
 call mpas_log_write('--- index_qpoma1_fg = $i',intArgs=(/index_qpoma1/))
 call mpas_log_write('--- index_qso4a1_fg = $i',intArgs=(/index_qso4a1/))
 call mpas_log_write('--- index_qsoa1_a1_fg = $i',intArgs=(/index_qsoa1_a1/))
 call mpas_log_write('--- index_qsoa2_a1_fg = $i',intArgs=(/index_qsoa2_a1/))
 call mpas_log_write('--- index_qsoa3_a1_fg = $i',intArgs=(/index_qsoa3_a1/))
 call mpas_log_write('--- index_qsoa4_a1_fg = $i',intArgs=(/index_qsoa4_a1/))
 call mpas_log_write('--- index_qsoa5_a1_fg = $i',intArgs=(/index_qsoa5_a1/))
 call mpas_log_write('--- index_numa1_fg  = $i',intArgs=(/index_numa1/) )
 call mpas_log_write(' ')
 call mpas_log_write('--- aitken mode:')
 call mpas_log_write('--- index_qdsta2_fg = $i',intArgs=(/index_qdsta2/))
 call mpas_log_write('--- index_qncla2_fg = $i',intArgs=(/index_qncla2/))
 call mpas_log_write('--- index_qso4a2_fg = $i',intArgs=(/index_qso4a2/))
 call mpas_log_write('--- index_qsoa1_a2_fg = $i',intArgs=(/index_qsoa1_a2/))
 call mpas_log_write('--- index_qsoa2_a2_fg = $i',intArgs=(/index_qsoa2_a2/))
 call mpas_log_write('--- index_qsoa3_a2_fg = $i',intArgs=(/index_qsoa3_a2/))
 call mpas_log_write('--- index_qsoa4_a2_fg = $i',intArgs=(/index_qsoa4_a2/))
 call mpas_log_write('--- index_qsoa5_a2_fg = $i',intArgs=(/index_qsoa5_a2/))
 call mpas_log_write('--- index_numa2_fg  = $i',intArgs=(/index_numa2/) )
 call mpas_log_write(' ')
 call mpas_log_write('--- coarse mode:')
 call mpas_log_write('--- index_qdsta3_fg = $i',intArgs=(/index_qdsta3/))
 call mpas_log_write('--- index_qncla3_fg = $i',intArgs=(/index_qncla3/))
 call mpas_log_write('--- index_qso4a3_fg = $i',intArgs=(/index_qso4a3/))
 call mpas_log_write('--- index_numa3_fg  = $i',intArgs=(/index_numa3/) )
 call mpas_log_write(' ')
 call mpas_log_write('--- primary carbon mode:')
 call mpas_log_write('--- index_qbca4_fg  = $i',intArgs=(/index_qbca4/) )
 call mpas_log_write('--- index_qpoma4_fg = $i',intArgs=(/index_qpoma4/))
 call mpas_log_write(' ')

 call mpas_log_write('--- gocart aerosols:')
 call mpas_log_write('--- index_qbc1_fg   = $i',intArgs=(/index_qbc1/)  )
 call mpas_log_write('--- index_qbc2_fg   = $i',intArgs=(/index_qbc2/)  )
 call mpas_log_write('--- index_qoc1_fg   = $i',intArgs=(/index_qoc1/)  )
 call mpas_log_write('--- index_qoc2_fg   = $i',intArgs=(/index_qoc2/)  )
 call mpas_log_write('--- index_qdust1_fg = $i',intArgs=(/index_qdust1/))
 call mpas_log_write('--- index_qdust2_fg = $i',intArgs=(/index_qdust2/))
 call mpas_log_write('--- index_qdust3_fg = $i',intArgs=(/index_qdust3/))
 call mpas_log_write('--- index_qdust4_fg = $i',intArgs=(/index_qdust4/))
 call mpas_log_write('--- index_qdust5_fg = $i',intArgs=(/index_qdust5/))
 call mpas_log_write('--- index_qseas1_fg = $i',intArgs=(/index_qseas1/))
 call mpas_log_write('--- index_qseas2_fg = $i',intArgs=(/index_qseas2/))
 call mpas_log_write('--- index_qseas3_fg = $i',intArgs=(/index_qseas3/))
 call mpas_log_write('--- index_qseas4_fg = $i',intArgs=(/index_qseas4/))
 call mpas_log_write(' ')
 call mpas_log_write('--- Thompson ice-friendly and water-friendly aerosols:')
 call mpas_log_write('--- index_nifa_fg   = $i',intArgs=(/index_nifa/)  )
 call mpas_log_write('--- index_nwfa_fg   = $i',intArgs=(/index_nwfa/)  )


 call mpas_pool_get_array(fg,'chems_fg',chems_fg)
 chems_fg(:,:,:) = 0._RKIND

 qco_fg    => chems_fg(index_qco,:,: )

 qbca1_fg  => chems_fg(index_qbca1,:,:)
 qdsta1_fg => chems_fg(index_qdsta1,:,:)
 qncla1_fg => chems_fg(index_qncla1,:,:)
 qpoma1_fg => chems_fg(index_qpoma1,:,:)
 qso4a1_fg => chems_fg(index_qso4a1,:,:)
 qsoa1_a1_fg => chems_fg(index_qsoa1_a1,:,:)
 qsoa2_a1_fg => chems_fg(index_qsoa2_a1,:,:)
 qsoa3_a1_fg => chems_fg(index_qsoa3_a1,:,:)
 qsoa4_a1_fg => chems_fg(index_qsoa4_a1,:,:)
 qsoa5_a1_fg => chems_fg(index_qsoa5_a1,:,:)
 numa1_fg  => chems_fg(index_numa1,:,:)

 qdsta3_fg => chems_fg(index_qdsta3,:,:)
 qncla3_fg => chems_fg(index_qncla3,:,:)
 qso4a3_fg => chems_fg(index_qso4a3,:,:)
 numa3_fg  => chems_fg(index_numa3,:,:)

 qdsta2_fg => chems_fg(index_qdsta2,:,:)
 qncla2_fg => chems_fg(index_qncla2,:,:)
 qso4a2_fg => chems_fg(index_qso4a2,:,:)
 qsoa1_a2_fg => chems_fg(index_qsoa1_a2,:,:)
 qsoa2_a2_fg => chems_fg(index_qsoa2_a2,:,:)
 qsoa3_a2_fg => chems_fg(index_qsoa3_a2,:,:)
 qsoa4_a2_fg => chems_fg(index_qsoa4_a2,:,:)
 qsoa5_a2_fg => chems_fg(index_qsoa5_a2,:,:)
 numa2_fg  => chems_fg(index_numa2,:,:)

 qbca4_fg  => chems_fg(index_qbca4,:,:)
 qpoma4_fg => chems_fg(index_qpoma4,:,:)

 qbc1_fg   => chems_fg(index_qbc1,:,:)
 qbc2_fg   => chems_fg(index_qbc2,:,:)
 qoc1_fg   => chems_fg(index_qoc1,:,:)
 qoc2_fg   => chems_fg(index_qoc2,:,:)
 qdust1_fg => chems_fg(index_qdust1,:,:)
 qdust2_fg => chems_fg(index_qdust2,:,:)
 qdust3_fg => chems_fg(index_qdust3,:,:)
 qdust4_fg => chems_fg(index_qdust4,:,:)
 qdust5_fg => chems_fg(index_qdust5,:,:)
 qseas1_fg => chems_fg(index_qseas1,:,:)
 qseas2_fg => chems_fg(index_qseas2,:,:)
 qseas3_fg => chems_fg(index_qseas3,:,:)
 qseas4_fg => chems_fg(index_qseas4,:,:)

 nifa_fg   => chems_fg(index_nifa,:,:)
 nwfa_fg   => chems_fg(index_nwfa,:,:)


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

    if(index(field % field, 'CO'    ) /= 0 .or. &
       index(field % field, 'BCA1'  ) /= 0 .or. &
       index(field % field, 'DSTA1' ) /= 0 .or. &
       index(field % field, 'NCLA1' ) /= 0 .or. &
       index(field % field, 'POMA1' ) /= 0 .or. &
       index(field % field, 'SO4A1' ) /= 0 .or. &
       index(field % field, 'SOA1A1') /= 0 .or. &
       index(field % field, 'SOA2A1') /= 0 .or. &
       index(field % field, 'SOA3A1') /= 0 .or. &
       index(field % field, 'SOA4A1') /= 0 .or. &
       index(field % field, 'SOA5A1') /= 0 .or. &
       index(field % field, 'DSTA2' ) /= 0 .or. &
       index(field % field, 'NCLA2' ) /= 0 .or. &
       index(field % field, 'SO4A2' ) /= 0 .or. &
       index(field % field, 'SOA1A2') /= 0 .or. &
       index(field % field, 'SOA2A2') /= 0 .or. &
       index(field % field, 'SOA3A2') /= 0 .or. &
       index(field % field, 'SOA4A2') /= 0 .or. &
       index(field % field, 'SOA5A2') /= 0 .or. &
       index(field % field, 'DSTA3' ) /= 0 .or. &
       index(field % field, 'NCLA3' ) /= 0 .or. &
       index(field % field, 'SO4A3' ) /= 0 .or. &
       index(field % field, 'BCA4'  ) /= 0 .or. &
       index(field % field, 'POMA4' ) /= 0 .or. &
       index(field % field, 'NUM_A1') /= 0 .or. &
       index(field % field, 'NUM_A2') /= 0 .or. &
       index(field % field, 'NUM_A3') /= 0 .or. &
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
       if(index(field%field,'CO') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating CO at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qco_fg
       elseif(index(field%field,'BCA1') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating BCA1 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qbca1_fg

       elseif(index(field%field,'DSTA1') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating DSTA1 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qdsta1_fg
       elseif(index(field%field,'NCLA1') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating NCLA1 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qncla1_fg
       elseif(index(field%field,'POMA1') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating POMA1 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qpoma1_fg
       elseif(index(field%field,'SO4A1') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating SO4A1 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qso4a1_fg
       elseif(index(field%field,'SOA1A1') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating SOA1A1 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qsoa1_a1_fg
       elseif(index(field%field,'SOA2A1') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating SOA2A1 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qsoa2_a1_fg
       elseif(index(field%field,'SOA3A1') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating SOA3A1 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qsoa3_a1_fg
       elseif(index(field%field,'SOA4A1') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating SOA4A1 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qsoa4_a1_fg
       elseif(index(field%field,'SOA5A1') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating SOA5A1 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qsoa5_a1_fg
       elseif(index(field%field,'DSTA2') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating DSTA2 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qdsta2_fg
       elseif(index(field%field,'NCLA2') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating NCLA2 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qncla2_fg
       elseif(index(field%field,'SO4A2') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating SO4A2 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qso4a2_fg
       elseif(index(field%field,'SOA1A2') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating SOA1A2 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qsoa1_a2_fg
       elseif(index(field%field,'SOA2A2') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating SOA2A2 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qsoa2_a2_fg
       elseif(index(field%field,'SOA3A2') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating SOA3A2 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qsoa3_a2_fg
       elseif(index(field%field,'SOA4A2') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating SOA4A2 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qsoa4_a2_fg
       elseif(index(field%field,'SOA5A2') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating SOA5A2 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qsoa5_a2_fg
       elseif(index(field%field,'DSTA3') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating DSTA3 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qdsta3_fg
       elseif(index(field%field,'NCLA3') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating NCLA3 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qncla3_fg
       elseif(index(field%field,'SO4A3') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating SO4A3 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qso4a3_fg
       elseif(index(field%field,'BCA4') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating BCA4 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qbca4_fg
       elseif(index(field%field,'POMA4') /= 0) then
          k = nCAMLevels+1-field%xlvl
          call mpas_log_write('Interpolating POMA4 at $i',intArgs=(/k/))
          nInterpPoints = nCells
          latPoints => latCell
          lonPoints => lonCell
          destField2d => qpoma4_fg
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


!--- check that all mixing ratios and other variables are all positive:
 call mpas_log_write(' ')
 call mpas_log_write('--- check for negative qco_fg:')
 do iCell = 1,nCells
    do k = 1,nCAMLevels
       if(qco_fg(k,iCell) .lt. 0._RKIND) then
          call mpas_log_write('--- negative qco_fg: $i $i $r',intArgs=(/iCell,k/), &
                              realArgs=(/qco_fg(k,iCell)/))
       endif
    enddo
 enddo

!coarse mode:
 call mpas_log_write('--- check for negative qdst3_fg:')
 do iCell = 1,nCells
    do k = 1,nCAMLevels
       if(qdsta3_fg(k,iCell) .lt. 0._RKIND) then
          call mpas_log_write('--- negative qdst3_fg: $i $i $r',intArgs=(/iCell,k/), &
                              realArgs=(/qdsta3_fg(k,iCell)/))
       endif
    enddo
 enddo
 call mpas_log_write('--- check for negative qncla3_fg:')
 do iCell = 1,nCells
    do k = 1,nCAMLevels
       if(qncla3_fg(k,iCell) .lt. 0._RKIND) then
          call mpas_log_write('--- negative qncla3_fg: $i $i $r',intArgs=(/iCell,k/), &
                              realArgs=(/qncla3_fg(k,iCell)/))
       endif
    enddo
 enddo
 call mpas_log_write('--- check for negative qso4a3_fg:')
 do iCell = 1,nCells
    do k = 1,nCAMLevels
       if(qso4a3_fg(k,iCell) .lt. 0._RKIND) then
          call mpas_log_write('--- negative qdst3_fg: $i $i $r',intArgs=(/iCell,k/), &
                              realArgs=(/qso4a3_fg(k,iCell)/))
       endif
    enddo
 enddo
 call mpas_log_write('--- check for negative numa3_fg:')
 do iCell = 1,nCells
    do k = 1,nCAMLevels
       if(numa3_fg(k,iCell) .lt. 0._RKIND) then
          call mpas_log_write('--- negative numa3_fg: $i $i $r',intArgs=(/iCell,k/), &
                              realArgs=(/numa3_fg(k,iCell)/))
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
 call mpas_log_write('--- check for negative dpcam_fg:')
 do iCell = 1,nCells
    do k = 1,nCAMLevels
       if(dpcam_fg(k,iCell) .lt. 0._RKIND) then
          call mpas_log_write('--- negative dpcam_fg: $i $i $r',intArgs=(/iCell,k/), &
                              realArgs=(/dpcam_fg(k,iCell)/))
       endif
    enddo
 enddo
!call mpas_log_write(' ')
!do iCell = 1, 5
!   do k = 1,nCAMLevels
!      call mpas_log_write('$i $r $r',intArgs=(/k/),realArgs=(/pcam_fg(k,iCell),dpcam_fg(k,iCell)/))
!   enddo
!   call mpas_log_write(' ')
!enddo


!--- convert mass mixing ratios in the accumulation, aitken, and coarse modes to number concentrations
!    and rescale the number concentrations:
 call init_scale_camchem(nCells,nCAMLevels,numa1_fg,numa2_fg,numa3_fg,qbca1_fg,qdsta1_fg,qncla1_fg,   &
                   qpoma1_fg,qso4a1_fg,qsoa1_a1_fg,qsoa2_a1_fg,qsoa3_a1_fg,qsoa4_a1_fg,qsoa5_a1_fg,   &
                   qncla2_fg,qso4a2_fg,qsoa1_a2_fg,qsoa2_a2_fg,qsoa3_a2_fg,qsoa4_a2_fg,qsoa5_a2_fg,   &
                   qdsta3_fg,qncla3_fg,qso4a3_fg,qbca4_fg,qpoma4_fg,qdust1_fg,qdust2_fg,qdust3_fg,    &
                   qdust4_fg,qdust5_fg,qbc1_fg,qbc2_fg,qoc1_fg,qoc2_fg,qseas1_fg,qseas2_fg,qseas3_fg, &
                   qseas4_fg,nifa_fg,nwfa_fg)

 call mpas_log_write('--- exit subroutine init_hinterp_camchem.')
 call mpas_log_write(' ')

 end subroutine init_hinterp_camchem


!=================================================================================================================
 subroutine init_scale_camchem(nCells,nCAMLevels,numa1_fg,numa2_fg,numa3_fg,qbca1_fg,qdsta1_fg,qncla1_fg, &
                       qpoma1_fg,qso4a1_fg,qsoa1_a1_fg,qsoa2_a1_fg,qsoa3_a1_fg,qsoa4_a1_fg,qsoa5_a1_fg,   &
                       qncla2_fg,qso4a2_fg,qsoa1_a2_fg,qsoa2_a2_fg,qsoa3_a2_fg,qsoa4_a2_fg,qsoa5_a2_fg,   &
                       qdsta3_fg,qncla3_fg,qso4a3_fg,qbca4_fg,qpoma4_fg,qdust1_fg,qdust2_fg,qdust3_fg,    &
                       qdust4_fg,qdust5_fg,qbc1_fg,qbc2_fg,qoc1_fg,qoc2_fg,qseas1_fg,qseas2_fg,qseas3_fg, &
                       qseas4_fg,nifa_fg,nwfa_fg)
!=================================================================================================================

!input arguments:
 integer,intent(in):: nCells,nCAMLevels
 real(kind=RKIND),dimension(nCAMLevels,nCells+1),intent(in):: numa1_fg,numa2_fg,numa3_fg

!inout arguments:
 real(kind=RKIND),intent(inout),dimension(nCAMLevels,nCells+1):: qbca1_fg,qdsta1_fg,qncla1_fg,qpoma1_fg,qso4a1_fg
 real(kind=RKIND),intent(inout),dimension(nCAMLevels,nCells+1):: qncla2_fg,qso4a2_fg
 real(kind=RKIND),intent(inout),dimension(nCAMLevels,nCells+1):: qso4a3_fg,qncla3_fg,qdsta3_fg
 real(kind=RKIND),intent(inout),dimension(nCAMLevels,nCells+1):: qbca4_fg,qpoma4_fg

 real(kind=RKIND),intent(inout),dimension(NCAMLevels,nCells+1):: qbc1_fg,qbc2_fg,qoc1_fg,qoc2_fg
 real(kind=RKIND),intent(inout),dimension(nCAMLevels,nCells+1):: qdust1_fg,qdust2_fg,qdust3_fg,qdust4_fg,qdust5_fg
 real(kind=RKIND),intent(inout),dimension(nCAMLevels,nCells+1):: qseas1_fg,qseas2_fg,qseas3_fg,qseas4_fg

 real(kind=RKIND),intent(inout),dimension(nCAMLevels,nCells+1):: qsoa1_a1_fg,qsoa2_a1_fg,qsoa3_a1_fg, &
                                                           qsoa4_a1_fg,qsoa5_a1_fg
 real(kind=RKIND),intent(inout),dimension(nCAMLevels,nCells+1):: qsoa1_a2_fg,qsoa2_a2_fg,qsoa3_a2_fg, &
                                                           qsoa4_a2_fg,qsoa5_a2_fg

 real(kind=RKIND),intent(inout),dimension(nCAMLevels,nCells+1):: nifa_fg,nwfa_fg


!local variables and arrays:
 integer:: iCell,k

!these parameters should be moved somewhere else:
 real(kind=RKIND),parameter:: pii = 3.141592653589793

 real(kind=RKIND),parameter:: dens_bc       = 1.7e03  ! density of black carbon (kg/m3).
 real(kind=RKIND),parameter:: dens_oc       = 1.0e03  ! density of organic carbon (kg/m3).
 real(kind=RKIND),parameter:: dens_soa      = 1.0e03  ! density of secondary organic carbon (kg/m3)
 real(kind=RKIND),parameter:: dens_so4      = 1.77e03 ! density of sulfates (kg/m3).
 real(kind=RKIND),parameter:: dens_ncl      = 1.9e03  ! density of seasalts (kg/m3).
 real(kind=RKIND),parameter:: dens_dst      = 2.6e03  ! density of dust (kg/m3).

 real(kind=RKIND),parameter:: dgnum_accum   = 1.1e-07 ! (m).
 real(kind=RKIND),parameter:: sigmag_accum  = 1.6     ! (-).
 real(kind=RKIND):: alnsg_accum,voltonumb_accum
 real(kind=RKIND):: sum_a1,ratio_a1

 real(kind=RKIND),parameter:: dgnum_aitken  = 2.6e-08 ! (m).
 real(kind=RKIND),parameter:: sigmag_aitken = 1.6     ! (-).
 real(kind=RKIND):: alnsg_aitken,voltonumb_aitken
 real(kind=RKIND):: sum_a2,ratio_a2

 real(kind=RKIND),parameter:: dgnum_coarse  = 9.0e-07 ! (m).
 real(kind=RKIND),parameter:: sigmag_coarse = 1.2     ! (-).
 real(kind=RKIND):: alnsg_coarse,voltonumb_coarse
 real(kind=RKIND):: sum_a3,ratio_a3
 real(kind=RKIND):: mr_to_nb_bc,mr_to_nb_oc,mr_to_nb_so4,mr_to_nb_ncl,mr_to_nb_dst,mr_to_nb_soa

 real(kind=RKIND),dimension(nCAMLevels,nCells+1):: nbca1,ndsta1,nncla1,npoma1,nso4a1
 real(kind=RKIND),dimension(nCAMLevels,nCells+1):: nsoa1_a1,nsoa2_a1,nsoa3_a1,nsoa4_a1,nsoa5_a1
 real(kind=RKIND),dimension(nCAMLevels,nCells+1):: nncla2,nso4a2
 real(kind=RKIND),dimension(nCAMLevels,nCells+1):: nsoa1_a2,nsoa2_a2,nsoa3_a2,nsoa4_a2,nsoa5_a2
 real(kind=RKIND),dimension(nCAMLevels,nCells+1):: ndsta3,nncla3,nso4a3

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine init_scale_camchem:')

!--- accumulation mode: conversion from mass mixing ratio to number concentrations and rescaling of numbers:
 alnsg_accum = log(sigmag_accum)
 voltonumb_accum = 1./((pii/6.)*(dgnum_accum**3.)*exp(4.5*alnsg_accum**2.))

 mr_to_nb_bc  = voltonumb_accum/dens_bc
 mr_to_nb_oc  = voltonumb_accum/dens_oc
 mr_to_nb_so4 = voltonumb_accum/dens_so4
 mr_to_nb_ncl = voltonumb_accum/dens_ncl
 mr_to_nb_dst = voltonumb_accum/dens_dst
 mr_to_nb_soa = voltonumb_accum/dens_soa

 call mpas_log_write('--- alnsg_accum      = $r',realArgs=(/alnsg_accum/))
 call mpas_log_write('--- voltonumb_accum  = $r',realArgs=(/voltonumb_accum/))
 call mpas_log_write('--- mr_to_nb_bc      = $r',realArgs=(/mr_to_nb_bc/))
 call mpas_log_write('--- mr_to_nb_dst     = $r',realArgs=(/mr_to_nb_dst/))
 call mpas_log_write('--- mr_to_nb_ncl     = $r',realArgs=(/mr_to_nb_ncl/))
 call mpas_log_write('--- mr_to_nb_oc      = $r',realArgs=(/mr_to_nb_oc/))
 call mpas_log_write('--- mr_to_nb_so4     = $r',realArgs=(/mr_to_nb_so4/))
 call mpas_log_write('--- mr_to_nb_soa     = $r',realArgs=(/mr_to_nb_soa/))
 call mpas_log_write(' ')

 do iCell = 1,nCells
    do k = 1,nCAMLevels
       nbca1(k,iCell)  = mr_to_nb_bc*qbca1_fg(k,iCell)
       ndsta1(k,iCell) = mr_to_nb_dst*qdsta1_fg(k,iCell)
       nncla1(k,iCell) = mr_to_nb_ncl*qncla1_fg(k,iCell)
       npoma1(k,iCell) = mr_to_nb_oc*qpoma1_fg(k,iCell)
       nso4a1(k,iCell) = mr_to_nb_so4*qso4a1_fg(k,iCell)

       nsoa1_a1(k,iCell) = mr_to_nb_soa*qsoa1_a1_fg(k,iCell)
       nsoa2_a1(k,iCell) = mr_to_nb_soa*qsoa2_a1_fg(k,iCell)
       nsoa3_a1(k,iCell) = mr_to_nb_soa*qsoa3_a1_fg(k,iCell)
       nsoa4_a1(k,iCell) = mr_to_nb_soa*qsoa4_a1_fg(k,iCell)
       nsoa5_a1(k,iCell) = mr_to_nb_soa*qsoa5_a1_fg(k,iCell)
    enddo

    do k = 1,nCAMLevels
       sum_a1 = nbca1(k,iCell) + ndsta1(k,iCell) + nncla1(k,iCell) + npoma1(k,iCell) + nso4a1(k,iCell) &
              + nsoa1_a1(k,iCell) + nsoa2_a1(k,iCell) + nsoa3_a1(k,iCell) + nsoa4_a1(k,iCell)          &
              + nsoa5_a1(k,iCell)
       ratio_a1 = numa1_fg(k,iCell)/sum_a1

       nbca1(k,iCell)  = ratio_a1*nbca1(k,iCell)
       ndsta1(k,iCell) = ratio_a1*ndsta1(k,iCell)
       nncla1(k,iCell) = ratio_a1*nncla1(k,iCell)
       npoma1(k,iCell) = ratio_a1*npoma1(k,iCell)
       nso4a1(k,iCell) = ratio_a1*nso4a1(k,iCell)

       nsoa1_a1(k,iCell) = ratio_a1*nsoa1_a1(k,iCell)
       nsoa2_a1(k,iCell) = ratio_a1*nsoa2_a1(k,iCell)
       nsoa3_a1(k,iCell) = ratio_a1*nsoa3_a1(k,iCell)
       nsoa4_a1(k,iCell) = ratio_a1*nsoa4_a1(k,iCell)
       nsoa5_a1(k,iCell) = ratio_a1*nsoa5_a1(k,iCell)
    enddo
 enddo


!--- aitken mode: conversion from mass mixing ratio to number concentrations and rescaling of numbers:
 alnsg_aitken = log(sigmag_aitken)
 voltonumb_aitken = 1./((pii/6.)*(dgnum_aitken**3.)*exp(4.5*alnsg_aitken**2.))

 mr_to_nb_ncl = voltonumb_aitken/dens_ncl
 mr_to_nb_so4 = voltonumb_aitken/dens_so4
 mr_to_nb_soa = voltonumb_aitken/dens_soa

 call mpas_log_write('--- alnsg_aitken     = $r',realArgs=(/alnsg_aitken/))
 call mpas_log_write('--- voltonumb_aitken = $r',realArgs=(/voltonumb_aitken/))
 call mpas_log_write('--- mr_to_nb_ncl     = $r',realArgs=(/mr_to_nb_ncl/))
 call mpas_log_write('--- mr_to_nb_so4     = $r',realArgs=(/mr_to_nb_so4/))
 call mpas_log_write('--- mr_to_nb_soa     = $r',realArgs=(/mr_to_nb_soa/))
 call mpas_log_write(' ')

 do iCell = 1,nCells
    do k = 1,nCAMLevels
       nncla2(k,iCell) = mr_to_nb_ncl*qncla2_fg(k,iCell)
       nso4a2(k,iCell) = mr_to_nb_so4*qso4a2_fg(k,iCell)

       nsoa1_a2(k,iCell) = mr_to_nb_soa*qsoa1_a2_fg(k,iCell)
       nsoa2_a2(k,iCell) = mr_to_nb_soa*qsoa2_a2_fg(k,iCell)
       nsoa3_a2(k,iCell) = mr_to_nb_soa*qsoa3_a2_fg(k,iCell)
       nsoa4_a2(k,iCell) = mr_to_nb_soa*qsoa4_a2_fg(k,iCell)
       nsoa5_a2(k,iCell) = mr_to_nb_soa*qsoa5_a2_fg(k,iCell)
    enddo

    do k = 1,nCAMLevels
       sum_a2 = nncla2(k,iCell) + nso4a2(k,iCell) + nsoa1_a2(k,iCell) + nsoa2_a2(k,iCell) &
              + nsoa3_a2(k,iCell) + nsoa4_a2(k,iCell) + nsoa5_a2(k,iCell)
       ratio_a2 = numa2_fg(k,iCell)/sum_a2

       nncla2(k,iCell) = ratio_a2*nncla2(k,iCell)
       nso4a2(k,iCell) = ratio_a2*nso4a2(k,iCell)

       nsoa1_a2(k,iCell) = ratio_a2*nsoa1_a2(k,iCell)
       nsoa2_a2(k,iCell) = ratio_a2*nsoa2_a2(k,iCell)
       nsoa3_a2(k,iCell) = ratio_a2*nsoa3_a2(k,iCell)
       nsoa4_a2(k,iCell) = ratio_a2*nsoa4_a2(k,iCell)
       nsoa5_a2(k,iCell) = ratio_a2*nsoa5_a2(k,iCell)
     enddo
 enddo


!--- coarse mode: conversion from mass mixing ratio to number concentrations and rescaling of numbers
!    in the layer adjacent to the surface:
 alnsg_coarse = log(sigmag_coarse)
 voltonumb_coarse = 1./((pii/6.)*(dgnum_coarse**3.)*exp(4.5*alnsg_coarse**2.))

 mr_to_nb_ncl = voltonumb_coarse/dens_ncl
 mr_to_nb_dst = voltonumb_coarse/dens_dst
 mr_to_nb_so4 = voltonumb_coarse/dens_so4

 call mpas_log_write('--- alnsg_coarse     = $r',realArgs=(/alnsg_coarse/))
 call mpas_log_write('--- voltonumb_coarse = $r',realArgs=(/voltonumb_coarse/))
 call mpas_log_write('--- mr_to_nb_so4     = $r',realArgs=(/mr_to_nb_so4/))
 call mpas_log_write('--- mr_to_nb_ncl     = $r',realArgs=(/mr_to_nb_ncl/))
 call mpas_log_write('--- mr_to_nb_dst     = $r',realArgs=(/mr_to_nb_dst/))
 call mpas_log_write(' ')

 do iCell = 1,nCells
    do k = 1,nCAMLevels
       ndsta3(k,iCell) = mr_to_nb_dst*qdsta3_fg(k,iCell)
       nncla3(k,iCell) = mr_to_nb_ncl*qncla3_fg(k,iCell)
       nso4a3(k,iCell) = mr_to_nb_so4*qso4a3_fg(k,iCell)
    enddo

    do k = 1,nCAMLevels
       sum_a3 = ndsta3(k,iCell) + nncla3(k,iCell) + nso4a3(k,iCell)
       ratio_a3 = numa3_fg(k,iCell)/sum_a3

       ndsta3(k,iCell) = ratio_a3*ndsta3(k,iCell)
       nncla3(k,iCell) = ratio_a3*nncla3(k,iCell)
       nso4a3(k,iCell) = ratio_a3*nso4a3(k,iCell)
     enddo
 enddo


!--- remapping of the CAM-Chem aerosols to the GOCART aerosols:
 do iCell = 1,nCells
    do k = 1,nCAMLevels
       qbc1_fg(k,iCell)   = qbca4_fg(k,iCell)
       qbc2_fg(k,iCell)   = qbca1_fg(k,iCell)

       qoc1_fg(k,iCell)   = qpoma4_fg(k,iCell)
       qoc2_fg(k,iCell)   = qpoma1_fg(k,iCell)

       qdust1_fg(k,iCell) = 0.02*qdsta3_fg(k,iCell)
       qdust2_fg(k,iCell) = 0.93*qdsta3_fg(k,iCell)
       qdust3_fg(k,iCell) = 0.05*qdsta3_fg(k,iCell)
       qdust4_fg(k,iCell) = 0.0_RKIND
       qdust5_fg(k,iCell) = 0.0_RKIND

       qseas1_fg(k,iCell) = 1.0*qncla1_fg(k,iCell) + 1.0*qncla2_fg(k,iCell)
       qseas2_fg(k,iCell) = 0.5*qncla3_fg(k,iCell)
       qseas3_fg(k,iCell) = 0.5*qncla3_fg(k,iCell)
       qseas4_fg(k,iCell) = 0.0_RKIND

       nifa_fg(k,iCell)   = ndsta3(k,iCell)
       nwfa_fg(k,iCell)   = numa1_fg(k,iCell)
    enddo
 enddo


 call mpas_log_write('--- end subroutine init_scale_camchem:')

 end subroutine init_scale_camchem

!=================================================================================================================
 subroutine init_vinterp_camchem(configs,mesh,fg,diag,diag_physics,chems_input,state)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: fg
 type(mpas_pool_type),intent(in):: diag

!inout arguments:
 type(mpas_pool_type),intent(inout):: chems_input
 type(mpas_pool_type),intent(inout):: state
 type(mpas_pool_type),intent(inout):: diag_physics

!local variables and arrays:
 integer:: iCell,k,kk,n,nn

 integer,pointer:: nCells,nCAMLevels,nVertLevels
 integer,pointer:: num_chems,num_chems_gocart,num_chems_fg
 integer,pointer:: gas_start,gas_end,gocart_start,gocart_end,gocart_fg_start,gocart_fg_end
 integer,pointer:: index_qco,index_qco_fg
 integer,pointer:: index_nwfa,index_nwfa_fg,index_nifa,index_nifa_fg
 integer,pointer:: index_chems,index_chems_fg
 integer,dimension(:),pointer:: indexToCellID

 real(kind=RKIND),dimension(:),pointer:: nifa2d,nwfa2d
 real(kind=RKIND),dimension(:,:),pointer:: pcam_fg,pressure
 real(kind=RKIND),dimension(:,:,:),pointer:: chems,chems_gocart,chems_fg
 real(kind=RKIND),dimension(:,:,:),pointer:: scalars

 real(kind=RKIND),parameter:: m_CO  = 28.01 ! molecular weight of carbon monoxide
 real(kind=RKIND),parameter:: m_air = 28.97 ! molecular weight of dry air.
 real(kind=RKIND):: conv
 real(kind=RKIND):: target_p
 real(kind=RKIND),dimension(:,:),allocatable:: sorted_arr

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine init_vinterp_camchem:')

 call mpas_pool_get_dimension(mesh,'nCells',nCells)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
 call mpas_pool_get_dimension(mesh,'nCAMLevels' ,nCAMLevels )

 call mpas_pool_get_array(mesh,'indexToCellID',indexToCellID)

 call mpas_pool_get_dimension(state,'num_chems' ,num_chems )
 call mpas_pool_get_dimension(state,'gas_start' ,gas_start )
 call mpas_pool_get_dimension(state,'gas_end'   ,gas_end   )
 call mpas_pool_get_dimension(state,'index_qco' ,index_qco )
 call mpas_pool_get_dimension(state,'index_nifa',index_nifa)
 call mpas_pool_get_dimension(state,'index_nwfa',index_nwfa)
 call mpas_log_write('--- num_chems        = $i',intArgs=(/num_chems/) )
 call mpas_log_write('--- gas_start        = $i',intArgs=(/gas_start/) )
 call mpas_log_write('--- gas_end          = $i',intArgs=(/gas_end/)   )
 call mpas_log_write('--- index_qco        = $i',intArgs=(/index_qco/) )
 call mpas_log_write('--- index_nifa       = $i',intArgs=(/index_nifa/))
 call mpas_log_write('--- index_nwfa       = $i',intArgs=(/index_nwfa/))
 call mpas_log_write(' ')

 call mpas_pool_get_dimension(chems_input,'num_chems_gocart',num_chems_gocart)
 call mpas_pool_get_dimension(chems_input,'gocart_start'    ,gocart_start    )
 call mpas_pool_get_dimension(chems_input,'gocart_end'      ,gocart_end      )
 call mpas_log_write('--- num_chems_gocart = $i',intArgs=(/num_chems_gocart/))
 call mpas_log_write('--- gocart_start     = $i',intArgs=(/gocart_start/)    )
 call mpas_log_write('--- gocart_end       = $i',intArgs=(/gocart_end/)      )
 call mpas_log_write(' ')

 call mpas_pool_get_dimension(fg,'num_chems_fg'   ,num_chems_fg   )
 call mpas_pool_get_dimension(fg,'gocart_fg_start',gocart_fg_start)
 call mpas_pool_get_dimension(fg,'gocart_fg_end'  ,gocart_fg_end  )
 call mpas_pool_get_dimension(fg,'index_qco_fg'   ,index_qco_fg   )
 call mpas_pool_get_dimension(fg,'index_nifa_fg'  ,index_nifa_fg  )
 call mpas_pool_get_dimension(fg,'index_nwfa_fg'  ,index_nwfa_fg  )
 call mpas_log_write('--- num_chems_fg    = $i',intArgs=(/num_chems_fg/)   )
 call mpas_log_write('--- gocart_fg_start = $i',intArgs=(/gocart_fg_start/))
 call mpas_log_write('--- gocart_fg_end   = $i',intArgs=(/gocart_fg_end/)  )
 call mpas_log_write('--- index_qco_fg    = $i',intArgs=(/index_qco_fg/)   )
 call mpas_log_write(' ')

 call mpas_pool_get_array(diag,'pressure',pressure)

 call mpas_pool_get_array(diag_physics,'nifa2d',nifa2d)
 call mpas_pool_get_array(diag_physics,'nwfa2d',nwfa2d)

 call mpas_pool_get_array(fg,'pcam_fg' ,pcam_fg )
 call mpas_pool_get_array(fg,'chems_fg',chems_fg) 
 call mpas_pool_get_array(chems_input,'chems_gocart',chems_gocart)
 call mpas_pool_get_array(state,'chems',chems)
 call mpas_pool_get_array(state,'scalars',scalars)
 chems_gocart(:,:,:) = 0._RKIND
 chems(:,:,:)        = 0._RKIND

 if(.not.allocated(sorted_arr)) allocate(sorted_arr(2,nCAMLevels))
!--- interpolation of the first-guess CAM-Chem aerosols (in Registry_init_chems.xml, interpolation of first-
!    guess fields in chems_fg (qbc1_fg to qsoa4a3_fg) to pressure-interpolated fields in chems_gocart (qbc1
!    to qsoa4a3, i.e. first 18 fields):
 call mpas_log_write('--- vertical interpolation of aerosols: num_chems_gocart = $i',intArgs=(/num_chems_gocart/))
 do n = 1,num_chems_gocart
    do iCell = 1,nCells
       sorted_arr(1,1:nCAMLevels) = 0._RKIND
       sorted_arr(2,1:nCAMLevels) = 0._RKIND
       do k = 1,nCAMLevels
          sorted_arr(1,k) = pcam_fg(k,iCell)
          sorted_arr(2,k) = chems_fg(n,k,iCell)
       enddo
!      if(pressure(1,iCell) .ge. sorted_arr(1,1)) then
!         call mpas_log_write('--- before vertical interpolation: PCAM levels')
!         do k = 1,nCAMLevels
!            call mpas_log_write('$i $i $i $i $r $r',intArgs=(/n,iCell,indexToCellID(iCell),k/), &
!                                                    realArgs=(/sorted_arr(1,k),sorted_arr(2,k)/))
!         enddo
!         call mpas_log_write('------:')
!      endif
       do k = nVertLevels,1,-1
          target_p = pressure(k,iCell)
          chems_gocart(n,k,iCell) = pressure_interp(iCell,k,target_p,nCAMLevels,sorted_arr(:,1:nCAMLevels))
          if(target_p.gt.sorted_arr(1,1)) chems_gocart(n,k,iCell) = chems_gocart(n,k+1,iCell)
       enddo
!      if(pressure(1,iCell) .ge. sorted_arr(1,1)) then
!         call mpas_log_write('--- after vertical interpolation: PRESSURE levels')
!         do k = 1,nVertLevels
!            call mpas_log_write('$i $i $i $r $r',intArgs=(/iCell,indexToCellID(iCell),k/), &
!                                realArgs=(/pressure(k,iCell),chems_gocart(n,k,iCell)/))
!         enddo
!         call mpas_log_write(' ')
!      endif
    enddo
 enddo
 call mpas_log_write('--- end vertical interpolation of aerosols:')

!--- interpolation of the CO first-guess CAM-Chem mixing ratio (in Registry_init_chems.mxl, index_qco_fg
!    equals 29, meaning that we have to loop beyond num_chems_gocart to find it). This is the reason why
!    the length of the loop below equals num_chems_fg instead of num_chems_gocart:
 do n = 1,num_chems_fg
    if(n .eq. index_qco_fg) then
       call mpas_log_write(' ')
       call mpas_log_write('--- vertical interpolation of CO mixing ratio: index_qco_fg = $i',intArgs=(/n/))
       call mpas_log_write('--- index_qco = $i',intArgs=(/index_qco/))
       conv = m_CO / m_air
       do iCell = 1,nCells
          sorted_arr(1,1:nCAMLevels) = 0._RKIND
          sorted_arr(2,1:nCAMLevels) = 0._RKIND
          do k = 1,nCAMLevels
             sorted_arr(1,k) = pcam_fg(k,iCell)
             sorted_arr(2,k) = chems_fg(n,k,iCell)
          enddo
!         if(pressure(1,iCell) .ge. sorted_arr(1,1)) then
!            call mpas_log_write('--- before vertical interpolation: PCAM levels')
!            do k = 1,nCAMLevels
!               call mpas_log_write('$i $i $i $i $r $r',intArgs=(/n,iCell,indexToCellID(iCell),k/), &
!                                                       realArgs=(/sorted_arr(1,k),sorted_arr(2,k)/))
!            enddo
!            call mpas_log_write('------:')
!         endif
          do k = nVertLevels,1,-1
             target_p = pressure(k,iCell)
             chems(index_qco,k,iCell) = pressure_interp(iCell,k,target_p,nCAMLevels,sorted_arr(:,1:nCAMLevels))
             chems(index_qco,k,iCell) = conv*chems(index_qco,k,iCell)
             if(target_p.gt.sorted_arr(1,1)) chems(index_qco,k,iCell) = chems(index_qco,k+1,iCell)
          enddo
!         if(pressure(1,iCell) .ge. sorted_arr(1,1)) then
!            call mpas_log_write('--- after vertical interpolation: PRESSURE levels')
!            do k = 1,nVertLevels
!               call mpas_log_write('$i $i $i $r $r',intArgs=(/iCell,indexToCellID(iCell),k/), &
!                                   realArgs=(/pressure(k,iCell),chems(index_qco,k,iCell)/))
!            enddo
!            call mpas_log_write(' ')
!         endif
       enddo
       call mpas_log_write('--- end vertical interpolation of CO mixing ratio:')
       call mpas_log_write(' ')
    endif
 enddo

!--- interpolation of the first-guess Thompson ice-friendly and water-friendly aerosols:
 call mpas_log_write('--- vertical interpolation of water- and ice-friendly Thompson aerosols:')
 nifa2d(:) = 0._RKIND
 nwfa2d(:) = 0._RKIND
 scalars(index_nifa,:,:) = 0._RKIND
 scalars(index_nwfa,:,:) = 0._RKIND
 do iCell = 1,nCells
    !nifa:
    sorted_arr(1,1:nCAMLevels) = 0._RKIND
    sorted_arr(2,1:nCAMLevels) = 0._RKIND
    do k = 1,nCAMLevels
       sorted_arr(1,k) = pcam_fg(k,iCell)
       sorted_arr(2,k) = chems_fg(index_nifa_fg,k,iCell)
    enddo
!   if(pressure(1,iCell) .ge. sorted_arr(1,1)) then
!      call mpas_log_write('--- before vertical interpolation: PCAM levels')
!      do k = 1,nCAMLevels
!         call mpas_log_write('$i $i $i $i $r $r',intArgs=(/n,iCell,indexToCellID(iCell),k/), &
!                                                 realArgs=(/sorted_arr(1,k),sorted_arr(2,k)/))
!      enddo
!      call mpas_log_write('------:')
!   endif
    do k = nVertLevels,1,-1
       target_p = pressure(k,iCell)
       scalars(index_nifa,k,iCell) = pressure_interp(iCell,k,target_p,nCAMLevels,sorted_arr(:,1:nCAMLevels))
       if(target_p.gt.sorted_arr(1,1)) scalars(index_nifa,k,iCell) = scalars(index_nifa,k+1,iCell)
    enddo
!   if(pressure(1,iCell) .ge. sorted_arr(1,1)) then
!      call mpas_log_write('--- after vertical interpolation: PRESSURE levels')
!      do k = 1,nVertLevels
!         call mpas_log_write('$i $i $i $r $r',intArgs=(/iCell,indexToCellID(iCell),k/), &
!                             realArgs=(/pressure(k,iCell),scalars(index_nifa,k,iCell)/))
!      enddo
!      call mpas_log_write(' ')
!   endif
    k = 1
    nifa2d(iCell) = scalars(index_nifa,k,iCell)

    !nwfa:
    sorted_arr(1,1:nCAMLevels) = 0._RKIND
    sorted_arr(2,1:nCAMLevels) = 0._RKIND
    do k = 1,nCAMLevels
       sorted_arr(1,k) = pcam_fg(k,iCell)
       sorted_arr(2,k) = chems_fg(index_nwfa_fg,k,iCell)
    enddo
!   if(pressure(1,iCell) .ge. sorted_arr(1,1)) then
!      call mpas_log_write('--- before vertical interpolation: PCAM levels')
!      do k = 1,nCAMLevels
!         call mpas_log_write('$i $i $i $i $r $r',intArgs=(/n,iCell,indexToCellID(iCell),k/), &
!                                                 realArgs=(/sorted_arr(1,k),sorted_arr(2,k)/))
!      enddo
!      call mpas_log_write('------:')
!   endif
    do k = nVertLevels,1,-1
       target_p = pressure(k,iCell)
       scalars(index_nwfa,k,iCell) = pressure_interp(iCell,k,target_p,nCAMLevels,sorted_arr(:,1:nCAMLevels))
       if(target_p.gt.sorted_arr(1,1)) scalars(index_nwfa,k,iCell) = scalars(index_nwfa,k+1,iCell)
    enddo
!   if(pressure(1,iCell) .ge. sorted_arr(1,1)) then
!      call mpas_log_write('--- after vertical interpolation: PRESSURE levels')
!      do k = 1,nVertLevels
!         call mpas_log_write('$i $i $i $r $r',intArgs=(/iCell,indexToCellID(iCell),k/), &
!                             realArgs=(/pressure(k,iCell),scalars(index_nwfa,k,iCell)/))
!      enddo
!      call mpas_log_write(' ')
!   endif
    k = 1
    nwfa2d(iCell) = scalars(index_nwfa,k,iCell)
 enddo

 if(allocated(sorted_arr)) deallocate(sorted_arr)


 call mpas_log_write('--- negative chem variables:')
 do n = 1,num_chems_gocart
    do iCell = 1,nCells
       do k = 1,nVertLevels
          if(chems_gocart(n,k,iCell) .lt. 0._RKIND) then
             call mpas_log_write('$i $i $i $r',intArgs=(/n,iCell,k/),realArgs=(/chems_gocart(n,k,iCell)/))
          endif
       enddo
    enddo
 enddo
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

       !"water-friendly" aerosols so that nwfa resembles the global distribution of the number
       !concentration of aerosols in the accumulation mode:
       qwfa(k,iCell) = 3*qoc2(k,iCell) + 0.4*qso4a1(k,iCell) &
                     + 50*qseas1(k,iCell) + 50*qseas2(k,iCell)

       nwfa(k,iCell) = 3*iconv_oc2*qoc2(k,iCell) + 0.4*iconv_so4*qso4a1(k,iCell) &
                     + 50*iconv_seas1*qseas1(k,iCell) + 50*iconv_seas2*qseas2(k,iCell)
    enddo
 enddo

 call mpas_log_write('--- end subroutine init_camchem_tothompson:')

 end subroutine init_camchem_tothompson

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
