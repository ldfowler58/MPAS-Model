! Copyright (c) 2024 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================================
 module mpas_chemistry_init_gocart2G_emissions
 use mpas_dmpar
 use mpas_log
 use mpas_kind_types
 use mpas_pool_routines
 use mpas_stream_manager
 use mpas_timekeeping,only: sub_t_t
 use mpas_derived_types,only     : mpas_time_type,mpas_timeinterval_type
 use mpas_timekeeping,only       : mpas_get_clock_time,mpas_set_time,mpas_get_timeInterval


 implicit none
 private
 public:: init_CAMS_emissions


!initialization of CAMS emissions.
!Laura D. Fowler (laura@ucar.edu) / 2022-02-08.


 contains


!==================================================================================================================
 subroutine init_CAMS_emissions(clock,stream_manager,mesh,CAMS_emissions)
!==================================================================================================================

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

!------------------------------------------------------------------------------------------------------------------
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

 end subroutine init_CAMS_emissions

!==================================================================================================================
 end module mpas_chemistry_init_gocart2G_emissions
!==================================================================================================================
