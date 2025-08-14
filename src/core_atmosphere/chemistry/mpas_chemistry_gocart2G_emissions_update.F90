! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================================
 module mpas_chemistry_gocart2G_emissions_update
 use mpas_kind_types
 use mpas_log
 use mpas_pool_routines
 use mpas_stream_manager
 use mpas_timekeeping,only  : sub_t_t
 use mpas_derived_types,only: mpas_time_type,mpas_timeinterval_type
 use mpas_timekeeping,only  : mpas_get_clock_time,mpas_set_time,mpas_get_timeInterval


 implicit none
 private
 public:: update_anth_emissions_bc,  &
          update_anth_emissions_co,  &
          update_anth_emissions_oc,  &
          update_anth_emissions_nh3, &
          update_anth_emissions_su,  &
          update_anth_emissions_iso, &
          update_anth_emissions_mnt


!update of anthropogenic emissions.
!laura D. Fowler (laura@ucar.edu) / 2025-08-12.


 contains


!=================================================================================================================
 subroutine update_anth_emissions_bc(clock,stream_manager,mesh,CAMS_anth_emissions,anth_emissions)
!=================================================================================================================

!--- input arguments:
 type(mpas_Clock_type),intent(in):: clock
 type(mpas_pool_type),intent(in),pointer:: mesh
 type(mpas_pool_type),intent(in),pointer:: CAMS_anth_emissions

!--- inout arguments:
 type(MPAS_streamManager_type),intent(inout):: stream_manager
 type(mpas_pool_type),intent(inout),pointer:: anth_emissions
 
!local variables and arrays:
 type(mpas_time_type):: beforeTime,afterTime,currTime
 type(mpas_timeinterval_type):: beforeDelta,afterDelta,totalDelta

 character(len=StrKIND):: actualTimeStamp
 integer,pointer:: nCells
 real(kind=RKIND):: total_dt,before_dt,after_dt

 real(kind=RKIND),dimension(:),pointer:: anth_less100m,anth_biofuel

!CAMS anthropogenic emissions:
 real(kind=RKIND),dimension(:),pointer:: bc_anth_res,bc_anth_sum

!gocart2G anthropogenic emissions:
 real(kind=RKIND),dimension(:),pointer:: bc_anth_less100m,bc_anth_biofuel

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine update_anth_emissions_bc:')


 call mpas_pool_get_dimension(mesh,'nCells',nCells)
 call mpas_log_write('    nCells      = $i',intArgs=(/nCells/))

 allocate(anth_less100m(nCells+1))
 allocate(anth_biofuel(nCells+1) )


!--- BC anthropogenic emissions:
 call mpas_pool_get_array(anth_emissions,'bc_anth_less100m',bc_anth_less100m)
 call mpas_pool_get_array(anth_emissions,'bc_anth_biofuel' ,bc_anth_biofuel )
 bc_anth_less100m(:) = 0._RKIND
 bc_anth_biofuel(:)  = 0._RKIND

 call mpas_pool_get_array(CAMS_anth_emissions,'bc_anth_res',bc_anth_res)
 call mpas_pool_get_array(CAMS_anth_emissions,'bc_anth_sum',bc_anth_sum)

 call mpas_stream_mgr_read(stream_manager,'anth_bc_emissions',rightNow=.true., &
                  whence=MPAS_STREAM_LATEST_BEFORE,actualWhen=actualTimestamp)
 call mpas_log_write('    latest time before is  '//trim(actualTimestamp))
 call mpas_set_time(beforeTime,dateTimeString=trim(actualTimestamp))

 anth_biofuel(1:nCells)  = bc_anth_res(1:nCells)
 anth_less100m(1:nCells) = bc_anth_sum(1:nCells) - bc_anth_res(1:nCells)

 call mpas_stream_mgr_read(stream_manager,'anth_bc_emissions',rightNow=.true., &
                  whence=MPAS_STREAM_EARLIEST_AFTER,actualWhen=actualTimestamp)
 call mpas_log_write('    earliest time after is '//trim(actualTimestamp))
 call mpas_set_time(afterTime,dateTimeString=trim(actualTimestamp))


!get current time and calculate time deltas between the times that were actually read and the current time.
!retrieve time deltas as real values:
 currTime = mpas_get_clock_time(clock,MPAS_NOW)
 totalDelta  = sub_t_t(afterTime,beforeTime)
 beforeDelta = sub_t_t(currTime,beforeTime)
 afterDelta  = sub_t_t(afterTime,currTime)

 call mpas_get_timeInterval(totalDelta ,dt=total_dt)
 call mpas_get_timeInterval(beforeDelta,dt=before_dt)
 call mpas_get_timeInterval(afterDelta , dt=after_dt)
 call mpas_log_write('    totalDelta  = $r',realArgs=(/total_dt/))
 call mpas_log_write('    beforeDelta = $r',realArgs=(/before_dt/))
 call mpas_log_write('    afterDelta  = $r',realArgs=(/after_dt/))


!interpolation of surface emissions to the current time:
 if(total_dt > 0.0_RKIND) then
    bc_anth_biofuel(:)  = (after_dt/total_dt)*anth_biofuel(:) + (before_dt/total_dt)*bc_anth_res(:)
    bc_anth_less100m(:) = (after_dt/total_dt)*anth_less100m(:) &
                        + (before_dt/total_dt)*(bc_anth_sum(:)-bc_anth_res(:))
 endif


 deallocate(anth_less100m)
 deallocate(anth_biofuel )

 call mpas_log_write('--- end subroutine update_anth_emissions_bc.')

 end subroutine update_anth_emissions_bc

!=================================================================================================================
 subroutine update_anth_emissions_oc(clock,stream_manager,mesh,CAMS_anth_emissions,anth_emissions)
!=================================================================================================================

!--- input arguments:
 type(mpas_Clock_type),intent(in):: clock
 type(mpas_pool_type),intent(in),pointer:: mesh
 type(mpas_pool_type),intent(in),pointer:: CAMS_anth_emissions

!--- inout arguments:
 type(MPAS_streamManager_type),intent(inout):: stream_manager
 type(mpas_pool_type),intent(inout),pointer:: anth_emissions

!local variables and arrays:
 type(mpas_time_type):: beforeTime,afterTime,currTime
 type(mpas_timeinterval_type):: beforeDelta,afterDelta,totalDelta

 character(len=StrKIND):: actualTimeStamp
 integer,pointer:: nCells
 real(kind=RKIND):: total_dt,before_dt,after_dt

 real(kind=RKIND),dimension(:),pointer:: anth_less100m,anth_biofuel

!CAMS anthropogenic emissions:
 real(kind=RKIND),dimension(:),pointer:: oc_anth_res,oc_anth_sum

!gocart2G anthropogenic emissions:
 real(kind=RKIND),dimension(:),pointer:: oc_anth_less100m,oc_anth_biofuel

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine update_anth_emissions_oc:')


 call mpas_pool_get_dimension(mesh,'nCells',nCells)
 call mpas_log_write('    nCells      = $i',intArgs=(/nCells/))

 allocate(anth_less100m(nCells+1))
 allocate(anth_biofuel(nCells+1) )


!--- OC anthropogenic emissions:
 call mpas_pool_get_array(anth_emissions,'oc_anth_less100m',oc_anth_less100m)
 call mpas_pool_get_array(anth_emissions,'oc_anth_biofuel' ,oc_anth_biofuel )
 oc_anth_less100m(:) = 0._RKIND
 oc_anth_biofuel(:)  = 0._RKIND

 call mpas_pool_get_array(CAMS_anth_emissions,'oc_anth_res',oc_anth_res)
 call mpas_pool_get_array(CAMS_anth_emissions,'oc_anth_sum',oc_anth_sum)

 call mpas_stream_mgr_read(stream_manager,'anth_oc_emissions',rightNow=.true., &
                  whence=MPAS_STREAM_LATEST_BEFORE,actualWhen=actualTimestamp)
 call mpas_log_write('    latest time before is  '//trim(actualTimestamp))
 call mpas_set_time(beforeTime,dateTimeString=trim(actualTimestamp))

 anth_biofuel(1:nCells)  = oc_anth_res(1:nCells)
 anth_less100m(1:nCells) = oc_anth_sum(1:nCells) - oc_anth_res(1:nCells)

 call mpas_stream_mgr_read(stream_manager,'anth_oc_emissions',rightNow=.true., &
                  whence=MPAS_STREAM_EARLIEST_AFTER,actualWhen=actualTimestamp)
 call mpas_log_write('    earliest time after is '//trim(actualTimestamp))
 call mpas_set_time(afterTime,dateTimeString=trim(actualTimestamp))


!get current time and calculate time deltas between the times that were actually read and the current time.
!retrieve time deltas as real values:
 currTime = mpas_get_clock_time(clock,MPAS_NOW)
 totalDelta  = sub_t_t(afterTime,beforeTime)
 beforeDelta = sub_t_t(currTime,beforeTime)
 afterDelta  = sub_t_t(afterTime,currTime)

 call mpas_get_timeInterval(totalDelta ,dt=total_dt)
 call mpas_get_timeInterval(beforeDelta,dt=before_dt)
 call mpas_get_timeInterval(afterDelta , dt=after_dt)
 call mpas_log_write('    totalDelta  = $r',realArgs=(/total_dt/))
 call mpas_log_write('    beforeDelta = $r',realArgs=(/before_dt/))
 call mpas_log_write('    afterDelta  = $r',realArgs=(/after_dt/))


!interpolation of surface emissions to the current time:
 if(total_dt > 0.0_RKIND) then
    oc_anth_biofuel(:)  = (after_dt/total_dt)*anth_biofuel(:) + (before_dt/total_dt)*oc_anth_res(:)
    oc_anth_less100m(:) = (after_dt/total_dt)*anth_less100m(:) &
                        + (before_dt/total_dt)*(oc_anth_sum(:)-oc_anth_res(:))
 endif


 deallocate(anth_less100m)
 deallocate(anth_biofuel )

 call mpas_log_write('--- end subroutine update_anth_emissions_oc.')

 end subroutine update_anth_emissions_oc

!=================================================================================================================
 subroutine update_anth_emissions_su(clock,stream_manager,mesh,CAMS_anth_emissions,anth_emissions)
!=================================================================================================================

!--- input arguments:
 type(mpas_Clock_type),intent(in):: clock
 type(mpas_pool_type),intent(in),pointer:: mesh
 type(mpas_pool_type),intent(in),pointer:: CAMS_anth_emissions

!--- inout arguments:
 type(MPAS_streamManager_type),intent(inout):: stream_manager
 type(mpas_pool_type),intent(inout),pointer:: anth_emissions
 
!local variables and arrays:
 type(mpas_time_type):: beforeTime,afterTime,currTime
 type(mpas_timeinterval_type):: beforeDelta,afterDelta,totalDelta

 character(len=StrKIND):: actualTimeStamp
 integer,pointer:: nCells
 real(kind=RKIND):: total_dt,before_dt,after_dt

 real(kind=RKIND),dimension(:),pointer:: anth_less100m,anth_biofuel

!CAMS anthropogenic emissions:
 real(kind=RKIND),dimension(:),pointer:: so2_anth_sum

!gocart2G anthropogenic emissions:
 real(kind=RKIND),dimension(:),pointer:: su_anth_less100m

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine update_anth_emissions_su:')


 call mpas_pool_get_dimension(mesh,'nCells',nCells)
 call mpas_log_write('    nCells      = $i',intArgs=(/nCells/))

 allocate(anth_less100m(nCells+1))


!--- SU anthropogenic emissions:
 call mpas_pool_get_array(anth_emissions,'su_anth_less100m',su_anth_less100m)
 su_anth_less100m(:) = 0._RKIND

 call mpas_pool_get_array(CAMS_anth_emissions,'so2_anth_sum',so2_anth_sum)

 call mpas_stream_mgr_read(stream_manager,'anth_so2_emissions',rightNow=.true., &
                  whence=MPAS_STREAM_LATEST_BEFORE,actualWhen=actualTimestamp)
 call mpas_log_write('    latest time before is  '//trim(actualTimestamp))
 call mpas_set_time(beforeTime,dateTimeString=trim(actualTimestamp))

 anth_less100m(1:nCells) = so2_anth_sum(1:nCells)

 call mpas_stream_mgr_read(stream_manager,'anth_so2_emissions',rightNow=.true., &
                  whence=MPAS_STREAM_EARLIEST_AFTER,actualWhen=actualTimestamp)
 call mpas_log_write('    earliest time after is '//trim(actualTimestamp))
 call mpas_set_time(afterTime,dateTimeString=trim(actualTimestamp))


!get current time and calculate time deltas between the times that were actually read and the current time.
!retrieve time deltas as real values:
 currTime = mpas_get_clock_time(clock,MPAS_NOW)
 totalDelta  = sub_t_t(afterTime,beforeTime)
 beforeDelta = sub_t_t(currTime,beforeTime)
 afterDelta  = sub_t_t(afterTime,currTime)

 call mpas_get_timeInterval(totalDelta ,dt=total_dt)
 call mpas_get_timeInterval(beforeDelta,dt=before_dt)
 call mpas_get_timeInterval(afterDelta , dt=after_dt)
 call mpas_log_write('    totalDelta  = $r',realArgs=(/total_dt/))
 call mpas_log_write('    beforeDelta = $r',realArgs=(/before_dt/))
 call mpas_log_write('    afterDelta  = $r',realArgs=(/after_dt/))


!interpolation of surface emissions to the current time:
 if(total_dt > 0.0_RKIND) then
    su_anth_less100m(:) = (after_dt/total_dt)*anth_less100m(:) + (before_dt/total_dt)*so2_anth_sum(:)
 endif


 deallocate(anth_less100m)

 call mpas_log_write('--- end subroutine update_anth_emissions_su.')

 end subroutine update_anth_emissions_su

!==================================================================================================================
 subroutine update_anth_emissions_nh3(clock,stream_manager,mesh,CAMS_anth_emissions,anth_emissions)
!==================================================================================================================

!input arguments:
 type(mpas_Clock_type),intent(in):: clock
 type(mpas_pool_type),intent(in),pointer:: mesh
 type(mpas_pool_type),intent(in),pointer:: CAMS_anth_emissions

!inout arguments:
 type(MPAS_streamManager_type),intent(inout):: stream_manager
 type(mpas_pool_type),intent(inout),pointer:: anth_emissions

!local variables and arrays:
 type(mpas_time_type):: beforeTime,afterTime,currTime
 type(mpas_timeinterval_type):: beforeDelta,afterDelta,totalDelta

 character(len=StrKIND):: actualTimeStamp
 integer,pointer:: nCells
 real(kind=RKIND):: total_dt,before_dt,after_dt

 real(kind=RKIND),dimension(:),pointer:: anth_ag,anth_en,anth_in,anth_oc,anth_re,anth_tr

!CAMS anthropogenic emissions:
 real(kind=RKIND),dimension(:),pointer:: nh3_anth_agl,nh3_anth_ags,nh3_anth_awb,nh3_anth_ene,nh3_anth_fef, &
                                         nh3_anth_ind,nh3_anth_res,nh3_anth_shp,nh3_anth_slv,nh3_anth_sum, &
                                         nh3_anth_swd,nh3_anth_tnr,nh3_anth_tro

!gocart2G anthropogenic emissions:
 real(kind=RKIND),dimension(:),pointer:: nh3_anth_ag,nh3_anth_en,nh3_anth_in,nh3_anth_oc,nh3_anth_re,nh3_anth_tr

!------------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine update_anth_emissions_nh3:')


 call mpas_pool_get_dimension(mesh,'nCells',nCells)
 call mpas_log_write('    nCells = $i',intArgs=(/nCells/))

 allocate(anth_ag(nCells+1))
 allocate(anth_en(nCells+1))
 allocate(anth_in(nCells+1))
 allocate(anth_oc(nCells+1))
 allocate(anth_re(nCells+1))
 allocate(anth_tr(nCells+1))


!--- NH3 anthropogenic emissions:
 call mpas_pool_get_array(anth_emissions,'nh3_anth_ag',nh3_anth_ag)
 call mpas_pool_get_array(anth_emissions,'nh3_anth_en',nh3_anth_en)
 call mpas_pool_get_array(anth_emissions,'nh3_anth_in',nh3_anth_in)
 call mpas_pool_get_array(anth_emissions,'nh3_anth_oc',nh3_anth_oc)
 call mpas_pool_get_array(anth_emissions,'nh3_anth_re',nh3_anth_re)
 call mpas_pool_get_array(anth_emissions,'nh3_anth_tr',nh3_anth_tr)
 nh3_anth_ag(:) = 0._RKIND
 nh3_anth_en(:) = 0._RKIND
 nh3_anth_in(:) = 0._RKIND
 nh3_anth_oc(:) = 0._RKIND
 nh3_anth_re(:) = 0._RKIND
 nh3_anth_tr(:) = 0._RKIND

 call mpas_pool_get_array(CAMS_anth_emissions,'nh3_anth_agl',nh3_anth_agl)
 call mpas_pool_get_array(CAMS_anth_emissions,'nh3_anth_ags',nh3_anth_ags)
 call mpas_pool_get_array(CAMS_anth_emissions,'nh3_anth_awb',nh3_anth_awb)
 call mpas_pool_get_array(CAMS_anth_emissions,'nh3_anth_ene',nh3_anth_ene)
 call mpas_pool_get_array(CAMS_anth_emissions,'nh3_anth_fef',nh3_anth_fef)
 call mpas_pool_get_array(CAMS_anth_emissions,'nh3_anth_ind',nh3_anth_ind)
 call mpas_pool_get_array(CAMS_anth_emissions,'nh3_anth_res',nh3_anth_res)
 call mpas_pool_get_array(CAMS_anth_emissions,'nh3_anth_shp',nh3_anth_shp)
 call mpas_pool_get_array(CAMS_anth_emissions,'nh3_anth_slv',nh3_anth_slv)
 call mpas_pool_get_array(CAMS_anth_emissions,'nh3_anth_sum',nh3_anth_sum)
 call mpas_pool_get_array(CAMS_anth_emissions,'nh3_anth_swd',nh3_anth_swd)
 call mpas_pool_get_array(CAMS_anth_emissions,'nh3_anth_tnr',nh3_anth_tnr)
 call mpas_pool_get_array(CAMS_anth_emissions,'nh3_anth_tro',nh3_anth_tro)

 call mpas_stream_mgr_read(stream_manager,'anth_nh3_emissions',rightNow=.true., &
                  whence=MPAS_STREAM_LATEST_BEFORE,actualWhen=actualTimestamp)
 call mpas_log_write('    latest time before is  '//trim(actualTimestamp))
 call mpas_set_time(beforeTime,dateTimeString=trim(actualTimestamp))


!agriculture emissions: agriculture livestock (agl) plus agriculture soils (ags) plus
!agriculture waste burning (awv) sectors:
 anth_ag(1:nCells) = nh3_anth_agl(1:nCells) + nh3_anth_ags(1:nCells) + nh3_anth_awb(1:nCells)

!energy emissions:
 anth_en(1:nCells) = nh3_anth_ene(1:nCells)

!industry emissions:
 anth_in(1:nCells) = nh3_anth_ind(1:nCells)

!ocean emissions:
 anth_oc(1:nCells) = 0._RKIND

!residential emissions:
 anth_re(1:nCells) = nh3_anth_res(1:nCells)

!transport emissions:
 anth_tr(1:nCells) = nh3_anth_tnr(1:nCells) + nh3_anth_tro(1:nCells)


 call mpas_stream_mgr_read(stream_manager,'anth_nh3_emissions',rightNow=.true., &
                  whence=MPAS_STREAM_EARLIEST_AFTER,actualWhen=actualTimestamp)
 call mpas_log_write('    earliest time after is '//trim(actualTimestamp))
 call mpas_set_time(afterTime,dateTimeString=trim(actualTimestamp))

 nh3_anth_ag(1:nCells) = nh3_anth_agl(1:nCells) + nh3_anth_ags(1:nCells) + nh3_anth_awb(1:nCells)
 nh3_anth_en(1:nCells) = nh3_anth_ene(1:nCells)
 nh3_anth_in(1:nCells) = nh3_anth_ind(1:nCells)
 nh3_anth_oc(1:nCells) = 0._RKIND
 nh3_anth_re(1:nCells) = nh3_anth_res(1:nCells)
 nh3_anth_tr(1:nCells) = nh3_anth_tnr(1:nCells) + nh3_anth_tro(1:nCells)


!get current time and calculate time deltas between the times that were actually read and the current time.
!retrieve time deltas as real values:
 currTime = mpas_get_clock_time(clock,MPAS_NOW)
 totalDelta  = sub_t_t(afterTime,beforeTime)
 beforeDelta = sub_t_t(currTime,beforeTime)
 afterDelta  = sub_t_t(afterTime,currTime)

 call mpas_get_timeInterval(totalDelta ,dt=total_dt)
 call mpas_get_timeInterval(beforeDelta,dt=before_dt)
 call mpas_get_timeInterval(afterDelta , dt=after_dt)
 call mpas_log_write('    totalDelta  = $r',realArgs=(/total_dt/))
 call mpas_log_write('    beforeDelta = $r',realArgs=(/before_dt/))
 call mpas_log_write('    afterDelta  = $r',realArgs=(/after_dt/))


!interpolation of surface emissions to the current time:
 if(total_dt > 0.0_RKIND) then
    nh3_anth_ag(:) = (after_dt/total_dt)*anth_ag(:) + (before_dt/total_dt)*nh3_anth_ag(:)
    nh3_anth_en(:) = (after_dt/total_dt)*anth_en(:) + (before_dt/total_dt)*nh3_anth_en(:)
    nh3_anth_in(:) = (after_dt/total_dt)*anth_in(:) + (before_dt/total_dt)*nh3_anth_in(:)
    nh3_anth_oc(:) = (after_dt/total_dt)*anth_oc(:) + (before_dt/total_dt)*nh3_anth_oc(:)
    nh3_anth_re(:) = (after_dt/total_dt)*anth_re(:) + (before_dt/total_dt)*nh3_anth_re(:)
    nh3_anth_tr(:) = (after_dt/total_dt)*anth_tr(:) + (before_dt/total_dt)*nh3_anth_tr(:)
 endif


 deallocate(anth_ag)
 deallocate(anth_en)
 deallocate(anth_in)
 deallocate(anth_oc)
 deallocate(anth_re)
 deallocate(anth_tr)

 call mpas_log_write('--- end subroutine update_anth_emissions_nh3.')

 end subroutine update_anth_emissions_nh3

!==================================================================================================================
 subroutine update_anth_emissions_co(clock,stream_manager,mesh,CAMS_anth_emissions,anth_emissions)
!==================================================================================================================

!input arguments:
 type(mpas_Clock_type),intent(in):: clock
 type(mpas_pool_type),intent(in),pointer:: mesh
 type(mpas_pool_type),intent(in),pointer:: CAMS_anth_emissions

!inout arguments:
 type(MPAS_streamManager_type),intent(inout):: stream_manager
 type(mpas_pool_type),intent(inout),pointer:: anth_emissions

!local variables and arrays:
 type(mpas_time_type):: beforeTime,afterTime,currTime
 type(mpas_timeinterval_type):: beforeDelta,afterDelta,totalDelta

 character(len=StrKIND):: actualTimeStamp
 integer,pointer:: nCells
 real(kind=RKIND):: total_dt,before_dt,after_dt

 real(kind=RKIND),dimension(:),pointer:: anth_em

!CAMS anthropogenic emissions:
 real(kind=RKIND),dimension(:),pointer:: co_anth_sum

!gocart2G anthropogenic emissions:
 real(kind=RKIND),dimension(:),pointer:: co_anth_em

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine update_anth_emissions_co:')


 call mpas_pool_get_dimension(mesh,'nCells',nCells)
 call mpas_log_write('    nCells      = $i',intArgs=(/nCells/))

 allocate(anth_em(nCells+1))


!--- CO anthropogenic emissions:
 call mpas_pool_get_array(anth_emissions,'co_anth_em',co_anth_em)
 co_anth_em(:) = 0._RKIND

 call mpas_pool_get_array(CAMS_anth_emissions,'co_anth_sum',co_anth_sum)

 call mpas_stream_mgr_read(stream_manager,'anth_co_emissions',rightNow=.true., &
                  whence=MPAS_STREAM_LATEST_BEFORE,actualWhen=actualTimestamp)
 call mpas_log_write('    latest time before is  '//trim(actualTimestamp))
 call mpas_set_time(beforeTime,dateTimeString=trim(actualTimestamp))

 anth_em(1:nCells) = co_anth_sum(1:nCells)

 call mpas_stream_mgr_read(stream_manager,'anth_co_emissions',rightNow=.true., &
                  whence=MPAS_STREAM_EARLIEST_AFTER,actualWhen=actualTimestamp)
 call mpas_log_write('    earliest time after is '//trim(actualTimestamp))
 call mpas_set_time(afterTime,dateTimeString=trim(actualTimestamp))


!get current time and calculate time deltas between the times that were actually read and the current time.
!retrieve time deltas as real values:
 currTime = mpas_get_clock_time(clock,MPAS_NOW)
 totalDelta  = sub_t_t(afterTime,beforeTime)
 beforeDelta = sub_t_t(currTime,beforeTime)
 afterDelta  = sub_t_t(afterTime,currTime)

 call mpas_get_timeInterval(totalDelta ,dt=total_dt)
 call mpas_get_timeInterval(beforeDelta,dt=before_dt)
 call mpas_get_timeInterval(afterDelta , dt=after_dt)
 call mpas_log_write('    totalDelta  = $r',realArgs=(/total_dt/))
 call mpas_log_write('    beforeDelta = $r',realArgs=(/before_dt/))
 call mpas_log_write('    afterDelta  = $r',realArgs=(/after_dt/))


!interpolation of surface emissions to the current time:
 if(total_dt > 0.0_RKIND) then
    co_anth_em(:) = (after_dt/total_dt)*anth_em(:) + (before_dt/total_dt)*co_anth_sum(:)
 endif


 deallocate(anth_em)

 call mpas_log_write('--- end subroutine update_anth_emissions_co.')

 end subroutine update_anth_emissions_co

!==================================================================================================================
 subroutine update_anth_emissions_iso(clock,stream_manager,mesh,CAMS_anth_emissions,anth_emissions)
!==================================================================================================================

!input arguments:
 type(mpas_Clock_type),intent(in):: clock
 type(mpas_pool_type),intent(in),pointer:: mesh
 type(mpas_pool_type),intent(in),pointer:: CAMS_anth_emissions

!inout arguments:
 type(MPAS_streamManager_type),intent(inout):: stream_manager
 type(mpas_pool_type),intent(inout),pointer:: anth_emissions

!local variables and arrays:
 type(mpas_time_type):: beforeTime,afterTime,currTime
 type(mpas_timeinterval_type):: beforeDelta,afterDelta,totalDelta

 character(len=StrKIND):: actualTimeStamp
 integer,pointer:: nCells
 real(kind=RKIND):: total_dt,before_dt,after_dt

 real(kind=RKIND),dimension(:),pointer:: anth_em

!CAMS anthropogenic emissions:
 real(kind=RKIND),dimension(:),pointer:: iso_anth_sum

!gocart2G anthropogenic emissions:
 real(kind=RKIND),dimension(:),pointer:: iso_anth_em

!------------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine update_anth_emissions_iso:')


 call mpas_pool_get_dimension(mesh,'nCells',nCells)
 call mpas_log_write('    nCells      = $i',intArgs=(/nCells/))

 allocate(anth_em(nCells+1))


!--- ISOPRENE anthropogenic emissions:
 call mpas_pool_get_array(anth_emissions,'iso_anth_em',iso_anth_em)
 iso_anth_em(:) = 0._RKIND

 call mpas_pool_get_array(CAMS_anth_emissions,'iso_anth_sum',iso_anth_sum)

 call mpas_stream_mgr_read(stream_manager,'anth_iso_emissions',rightNow=.true., &
                  whence=MPAS_STREAM_LATEST_BEFORE,actualWhen=actualTimestamp)
 call mpas_log_write('    latest time before is  '//trim(actualTimestamp))
 call mpas_set_time(beforeTime,dateTimeString=trim(actualTimestamp))

 anth_em(1:nCells) = iso_anth_sum(1:nCells)

 call mpas_stream_mgr_read(stream_manager,'anth_iso_emissions',rightNow=.true., &
                  whence=MPAS_STREAM_EARLIEST_AFTER,actualWhen=actualTimestamp)
 call mpas_log_write('    earliest time after is '//trim(actualTimestamp))
 call mpas_set_time(afterTime,dateTimeString=trim(actualTimestamp))


!get current time and calculate time deltas between the times that were actually read and the current time.
!retrieve time deltas as real values:
 currTime = mpas_get_clock_time(clock,MPAS_NOW)
 totalDelta  = sub_t_t(afterTime,beforeTime)
 beforeDelta = sub_t_t(currTime,beforeTime)
 afterDelta  = sub_t_t(afterTime,currTime)

 call mpas_get_timeInterval(totalDelta ,dt=total_dt)
 call mpas_get_timeInterval(beforeDelta,dt=before_dt)
 call mpas_get_timeInterval(afterDelta , dt=after_dt)
 call mpas_log_write('    totalDelta  = $r',realArgs=(/total_dt/))
 call mpas_log_write('    beforeDelta = $r',realArgs=(/before_dt/))
 call mpas_log_write('    afterDelta  = $r',realArgs=(/after_dt/))


!interpolation of surface emissions to the current time:
 if(total_dt > 0.0_RKIND) then
    iso_anth_em(:) = (after_dt/total_dt)*anth_em(:) + (before_dt/total_dt)*iso_anth_sum(:)
 endif


 deallocate(anth_em)

 call mpas_log_write('--- end subroutine update_anth_emissions_iso.')

 end subroutine update_anth_emissions_iso

!==================================================================================================================
 subroutine update_anth_emissions_mnt(clock,stream_manager,mesh,CAMS_anth_emissions,anth_emissions)
!==================================================================================================================

!input arguments:
 type(mpas_Clock_type),intent(in):: clock
 type(mpas_pool_type),intent(in),pointer:: mesh
 type(mpas_pool_type),intent(in),pointer:: CAMS_anth_emissions

!inout arguments:
 type(MPAS_streamManager_type),intent(inout):: stream_manager
 type(mpas_pool_type),intent(inout),pointer:: anth_emissions

!local variables and arrays:
 type(mpas_time_type):: beforeTime,afterTime,currTime
 type(mpas_timeinterval_type):: beforeDelta,afterDelta,totalDelta

 character(len=StrKIND):: actualTimeStamp
 integer,pointer:: nCells
 integer:: iCell
 real(kind=RKIND):: total_dt,before_dt,after_dt

 real(kind=RKIND),dimension(:),pointer:: anth_em

!CAMS anthropogenic emissions:
 real(kind=RKIND),dimension(:),pointer:: mnt_anth_sum

!gocart2G anthropogenic emissions:
 real(kind=RKIND),dimension(:),pointer:: mnt_anth_em

!------------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine update_anth_emissions_mnt:')


 call mpas_pool_get_dimension(mesh,'nCells',nCells)
 call mpas_log_write('    nCells      = $i',intArgs=(/nCells/))

 allocate(anth_em(nCells+1))


!--- MONOTERPENES anthropogenic emissions:
 call mpas_pool_get_array(anth_emissions,'mnt_anth_em',mnt_anth_em)
 mnt_anth_em(:) = 0._RKIND

 call mpas_pool_get_array(CAMS_anth_emissions,'mnt_anth_sum',mnt_anth_sum)

 call mpas_stream_mgr_read(stream_manager,'anth_mnt_emissions',rightNow=.true., &
                  whence=MPAS_STREAM_LATEST_BEFORE,actualWhen=actualTimestamp)
 call mpas_log_write('    latest time before is  '//trim(actualTimestamp))
 call mpas_set_time(beforeTime,dateTimeString=trim(actualTimestamp))

 anth_em(1:nCells) = mnt_anth_sum(1:nCells)

 call mpas_stream_mgr_read(stream_manager,'anth_mnt_emissions',rightNow=.true., &
                  whence=MPAS_STREAM_EARLIEST_AFTER,actualWhen=actualTimestamp)
 call mpas_log_write('    earliest time after is '//trim(actualTimestamp))
 call mpas_set_time(afterTime,dateTimeString=trim(actualTimestamp))


!get current time and calculate time deltas between the times that were actually read and the current time.
!retrieve time deltas as real values:
 currTime = mpas_get_clock_time(clock,MPAS_NOW)
 totalDelta  = sub_t_t(afterTime,beforeTime)
 beforeDelta = sub_t_t(currTime,beforeTime)
 afterDelta  = sub_t_t(afterTime,currTime)

 call mpas_get_timeInterval(totalDelta ,dt=total_dt)
 call mpas_get_timeInterval(beforeDelta,dt=before_dt)
 call mpas_get_timeInterval(afterDelta , dt=after_dt)
 call mpas_log_write('    totalDelta  = $r',realArgs=(/total_dt/))
 call mpas_log_write('    beforeDelta = $r',realArgs=(/before_dt/))
 call mpas_log_write('    afterDelta  = $r',realArgs=(/after_dt/))


!interpolation of surface emissions to the current time:
 if(total_dt > 0.0_RKIND) then
    mnt_anth_em(:) = (after_dt/total_dt)*anth_em(:) + (before_dt/total_dt)*mnt_anth_sum(:)
 endif

 deallocate(anth_em)


 call mpas_log_write('--- end subroutine update_anth_emissions_mnt.')

 end subroutine update_anth_emissions_mnt

!=================================================================================================================
 end module mpas_chemistry_gocart2G_emissions_update
!=================================================================================================================
