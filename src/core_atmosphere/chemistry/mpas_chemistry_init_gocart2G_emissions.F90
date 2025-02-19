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
 public:: init_anth_emissions_bc,         &
          init_anth_emissions_oc,         &
          init_anth_emissions_nh3,        &
          init_anth_emissions_su,         &
          init_biomass_burning_emissions, &
          init_BIOG_emissions


!initialization of anthropogenic emissions.
!Laura D. Fowler (laura@ucar.edu) / 2022-02-08.


 contains


!==================================================================================================================
 subroutine init_anth_emissions_bc(clock,stream_manager,mesh,CAMS_anth_emissions,anth_emissions)
!==================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in),pointer:: mesh
 type(mpas_pool_type),intent(in),pointer:: CAMS_anth_emissions
 type(mpas_Clock_type),intent(in),pointer:: clock

!inout arguments:
 type(MPAS_streamManager_type),intent(inout):: stream_manager
 type(mpas_pool_type),intent(inout),pointer:: anth_emissions

!local variables and arrays:
 type(mpas_time_type):: beforeTime,afterTime,currTime
 type(mpas_timeinterval_type):: beforeDelta,afterDelta,totalDelta

 character(len=StrKIND):: actualTimeStamp
 integer,pointer:: nCells,nVertLevels
 integer:: iCell
 real(kind=RKIND):: total_dt,before_dt,after_dt

 real(kind=RKIND),dimension(:),pointer:: anth_less100m,anth_less500m,anth_ship,anth_aviation_lto, &
                                         anth_aviation_cds,anth_aviation_crs
 real(kind=RKIND),dimension(:,:),pointer:: anth_aircraft

!CAMS anthropogenic emissions:
 real(kind=RKIND),dimension(:),pointer:: bc_anth_agl,bc_anth_ags,bc_anth_awb,bc_anth_ene,bc_anth_fef, &
                                         bc_anth_ind,bc_anth_res,bc_anth_shp,bc_anth_slv,bc_anth_sum, &
                                         bc_anth_swd,bc_anth_tnr,bc_anth_tro

!gocart2G anthropogenic emissions:
 real(kind=RKIND),dimension(:),pointer:: bc_anth_less100m,bc_anth_less500m,bc_anth_ship,bc_anth_aviation_lto, &
                                         bc_anth_aviation_cds,bc_anth_aviation_crs
 real(kind=RKIND),dimension(:,:),pointer:: bc_anth_aircraft

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine init_anth_emissions_bc:')


 call mpas_pool_get_dimension(mesh,'nCells',nCells)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
 call mpas_log_write('    nCells      = $i',intArgs=(/nCells/))
 call mpas_log_write('    nVertLevels = $i',intArgs=(/nVertLevels/))


 allocate(anth_less100m(nCells+1)    )
 allocate(anth_less500m(nCells+1)    )
 allocate(anth_ship(nCells+1)        )
 allocate(anth_aviation_lto(nCells+1))
 allocate(anth_aviation_cds(nCells+1))
 allocate(anth_aviation_crs(nCells+1))
 allocate(anth_aircraft(nVertLevels,1:nCells+1))


!--- BC anthropogenic emissions:
!    here, we assume that anthropogenic emissions from ships and aircrafts are not included. they do not need to
!    be interpolated to the initial time of the forecast and are set to zero. furthermore, we assume that surface
!    emissions are all non-energy related emissions and simply interpolated to layers below 100 meters.

 call mpas_pool_get_array(anth_emissions,'bc_anth_less100m'    ,bc_anth_less100m    )
 call mpas_pool_get_array(anth_emissions,'bc_anth_less500m'    ,bc_anth_less500m    )
 call mpas_pool_get_array(anth_emissions,'bc_anth_ship'        ,bc_anth_ship        )
 call mpas_pool_get_array(anth_emissions,'bc_anth_aviation_lto',bc_anth_aviation_lto)
 call mpas_pool_get_array(anth_emissions,'bc_anth_aviation_cds',bc_anth_aviation_cds)
 call mpas_pool_get_array(anth_emissions,'bc_anth_aviation_crs',bc_anth_aviation_crs)
 call mpas_pool_get_array(anth_emissions,'bc_anth_aircraft'    ,bc_anth_aircraft    )
 bc_anth_less100m(:)     = 0._RKIND
 bc_anth_less500m(:)     = 0._RKIND
 bc_anth_ship(:)         = 0._RKIND
 bc_anth_aviation_lto(:) = 0._RKIND
 bc_anth_aviation_cds(:) = 0._RKIND
 bc_anth_aviation_crs(:) = 0._RKIND
 bc_anth_aircraft(:,:)   = 0._RKIND


 call mpas_pool_get_array(CAMS_anth_emissions,'bc_anth_agl',bc_anth_agl)
 call mpas_pool_get_array(CAMS_anth_emissions,'bc_anth_ags',bc_anth_ags)
 call mpas_pool_get_array(CAMS_anth_emissions,'bc_anth_awb',bc_anth_awb)
 call mpas_pool_get_array(CAMS_anth_emissions,'bc_anth_ene',bc_anth_ene)
 call mpas_pool_get_array(CAMS_anth_emissions,'bc_anth_fef',bc_anth_fef)
 call mpas_pool_get_array(CAMS_anth_emissions,'bc_anth_ind',bc_anth_ind)
 call mpas_pool_get_array(CAMS_anth_emissions,'bc_anth_res',bc_anth_res)
 call mpas_pool_get_array(CAMS_anth_emissions,'bc_anth_shp',bc_anth_shp)
 call mpas_pool_get_array(CAMS_anth_emissions,'bc_anth_slv',bc_anth_slv)
 call mpas_pool_get_array(CAMS_anth_emissions,'bc_anth_sum',bc_anth_sum)
 call mpas_pool_get_array(CAMS_anth_emissions,'bc_anth_swd',bc_anth_swd)
 call mpas_pool_get_array(CAMS_anth_emissions,'bc_anth_tnr',bc_anth_tnr)
 call mpas_pool_get_array(CAMS_anth_emissions,'bc_anth_tro',bc_anth_tro)

 call mpas_stream_mgr_read(stream_manager,'anth_bc_emissions',rightNow=.true., &
                  whence=MPAS_STREAM_LATEST_BEFORE,actualWhen=actualTimestamp)
 call mpas_log_write('    latest time before is  '//trim(actualTimestamp))
 call mpas_set_time(beforeTime,dateTimeString=trim(actualTimestamp))


 anth_less100m(1:nCells) = bc_anth_sum(1:nCells)


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
    bc_anth_less100m(:) = (after_dt/total_dt)*anth_less100m(:) + (before_dt/total_dt)*bc_anth_sum(:)
 endif


 deallocate(anth_less100m    )
 deallocate(anth_less500m    )
 deallocate(anth_ship        )
 deallocate(anth_aviation_lto)
 deallocate(anth_aviation_cds)
 deallocate(anth_aviation_crs)
 deallocate(anth_aircraft    )


 call mpas_log_write('--- end subroutine init_anth_emissions_bc.')

 end subroutine init_anth_emissions_bc

!==================================================================================================================
 subroutine init_anth_emissions_oc(clock,stream_manager,mesh,CAMS_anth_emissions,anth_emissions)
!==================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in),pointer:: mesh
 type(mpas_pool_type),intent(in),pointer:: CAMS_anth_emissions
 type(mpas_Clock_type),intent(in),pointer:: clock

!inout arguments:
 type(MPAS_streamManager_type),intent(inout):: stream_manager
 type(mpas_pool_type),intent(inout),pointer:: anth_emissions

!local variables and arrays:
 type(mpas_time_type):: beforeTime,afterTime,currTime
 type(mpas_timeinterval_type):: beforeDelta,afterDelta,totalDelta

 character(len=StrKIND):: actualTimeStamp
 integer,pointer:: nCells,nVertLevels
 integer:: iCell
 real(kind=RKIND):: total_dt,before_dt,after_dt

 real(kind=RKIND),dimension(:),pointer:: anth_less100m,anth_less500m,anth_ship,anth_aviation_lto, &
                                         anth_aviation_cds,anth_aviation_crs
 real(kind=RKIND),dimension(:,:),pointer:: anth_aircraft

!CAMS anthropogenic emissions:
 real(kind=RKIND),dimension(:),pointer:: oc_anth_agl,oc_anth_ags,oc_anth_awb,oc_anth_ene,oc_anth_fef, &
                                         oc_anth_ind,oc_anth_res,oc_anth_shp,oc_anth_slv,oc_anth_sum, &
                                         oc_anth_swd,oc_anth_tnr,oc_anth_tro

!gocart2G anthropogenic emissions:
 real(kind=RKIND),dimension(:),pointer:: oc_anth_less100m,oc_anth_less500m,oc_anth_ship,oc_anth_aviation_lto, &
                                         oc_anth_aviation_cds,oc_anth_aviation_crs
 real(kind=RKIND),dimension(:,:),pointer:: oc_anth_aircraft

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine init_anth_emissions_oc:')


 call mpas_pool_get_dimension(mesh,'nCells',nCells)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
 call mpas_log_write('    nCells      = $i',intArgs=(/nCells/))
 call mpas_log_write('    nVertLevels = $i',intArgs=(/nVertLevels/))


 allocate(anth_less100m(nCells+1)    )
 allocate(anth_less500m(nCells+1)    )
 allocate(anth_ship(nCells+1)        )
 allocate(anth_aviation_lto(nCells+1))
 allocate(anth_aviation_cds(nCells+1))
 allocate(anth_aviation_crs(nCells+1))
 allocate(anth_aircraft(nVertLevels,1:nCells+1))


!--- OC anthropogenic emissions:
!    here, we assume that anthropogenic emissions from ships and aircrafts are not included. they do not need to
!    be interpolated to the initial time of the forecast and are set to zero. furthermore, we assume that surface
!    emissions are all non-energy related emissions and simply interpolated to layers below 100 meters.

 call mpas_pool_get_array(anth_emissions,'oc_anth_less100m'    ,oc_anth_less100m    )
 call mpas_pool_get_array(anth_emissions,'oc_anth_less500m'    ,oc_anth_less500m    )
 call mpas_pool_get_array(anth_emissions,'oc_anth_ship'        ,oc_anth_ship        )
 call mpas_pool_get_array(anth_emissions,'oc_anth_aviation_lto',oc_anth_aviation_lto)
 call mpas_pool_get_array(anth_emissions,'oc_anth_aviation_cds',oc_anth_aviation_cds)
 call mpas_pool_get_array(anth_emissions,'oc_anth_aviation_crs',oc_anth_aviation_crs)
 call mpas_pool_get_array(anth_emissions,'oc_anth_aircraft'    ,oc_anth_aircraft    )
 oc_anth_less100m(:)     = 0._RKIND
 oc_anth_less500m(:)     = 0._RKIND
 oc_anth_ship(:)         = 0._RKIND
 oc_anth_aviation_lto(:) = 0._RKIND
 oc_anth_aviation_cds(:) = 0._RKIND
 oc_anth_aviation_crs(:) = 0._RKIND
 oc_anth_aircraft(:,:)   = 0._RKIND


 call mpas_pool_get_array(CAMS_anth_emissions,'oc_anth_agl',oc_anth_agl)
 call mpas_pool_get_array(CAMS_anth_emissions,'oc_anth_ags',oc_anth_ags)
 call mpas_pool_get_array(CAMS_anth_emissions,'oc_anth_awb',oc_anth_awb)
 call mpas_pool_get_array(CAMS_anth_emissions,'oc_anth_ene',oc_anth_ene)
 call mpas_pool_get_array(CAMS_anth_emissions,'oc_anth_fef',oc_anth_fef)
 call mpas_pool_get_array(CAMS_anth_emissions,'oc_anth_ind',oc_anth_ind)
 call mpas_pool_get_array(CAMS_anth_emissions,'oc_anth_res',oc_anth_res)
 call mpas_pool_get_array(CAMS_anth_emissions,'oc_anth_shp',oc_anth_shp)
 call mpas_pool_get_array(CAMS_anth_emissions,'oc_anth_slv',oc_anth_slv)
 call mpas_pool_get_array(CAMS_anth_emissions,'oc_anth_sum',oc_anth_sum)
 call mpas_pool_get_array(CAMS_anth_emissions,'oc_anth_swd',oc_anth_swd)
 call mpas_pool_get_array(CAMS_anth_emissions,'oc_anth_tnr',oc_anth_tnr)
 call mpas_pool_get_array(CAMS_anth_emissions,'oc_anth_tro',oc_anth_tro)

 call mpas_stream_mgr_read(stream_manager,'anth_oc_emissions',rightNow=.true., &
                  whence=MPAS_STREAM_LATEST_BEFORE,actualWhen=actualTimestamp)
 call mpas_log_write('    latest time before is  '//trim(actualTimestamp))
 call mpas_set_time(beforeTime,dateTimeString=trim(actualTimestamp))


 anth_less100m(1:nCells) = oc_anth_sum(1:nCells)


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
    oc_anth_less100m(:) = (after_dt/total_dt)*anth_less100m(:) + (before_dt/total_dt)*oc_anth_sum(:)
 endif


 deallocate(anth_less100m    )
 deallocate(anth_less500m    )
 deallocate(anth_ship        )
 deallocate(anth_aviation_lto)
 deallocate(anth_aviation_cds)
 deallocate(anth_aviation_crs)
 deallocate(anth_aircraft    )


 call mpas_log_write('--- end subroutine init_anth_emissions_oc.')

 end subroutine init_anth_emissions_oc

!==================================================================================================================
 subroutine init_anth_emissions_su(clock,stream_manager,mesh,CAMS_anth_emissions,anth_emissions)
!==================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in),pointer:: mesh
 type(mpas_pool_type),intent(in),pointer:: CAMS_anth_emissions
 type(mpas_Clock_type),intent(in),pointer:: clock

!inout arguments:
 type(MPAS_streamManager_type),intent(inout):: stream_manager
 type(mpas_pool_type),intent(inout),pointer:: anth_emissions

!local variables and arrays:
 type(mpas_time_type):: beforeTime,afterTime,currTime
 type(mpas_timeinterval_type):: beforeDelta,afterDelta,totalDelta

 character(len=StrKIND):: actualTimeStamp
 integer,pointer:: nCells,nVertLevels
 integer:: iCell
 real(kind=RKIND):: total_dt,before_dt,after_dt

 real(kind=RKIND),dimension(:),pointer:: anth_less100m,anth_less500m,anth_shipso2,anth_shipso4, &
                                         anth_aviation_lto,anth_aviation_cds,anth_aviation_crs
 real(kind=RKIND),dimension(:,:),pointer:: anth_aircraft

!CAMS anthropogenic emissions:
 real(kind=RKIND),dimension(:),pointer:: so2_anth_agl,so2_anth_ags,so2_anth_awb,so2_anth_ene,so2_anth_fef, &
                                         so2_anth_ind,so2_anth_res,so2_anth_shp,so2_anth_slv,so2_anth_sum, &
                                         so2_anth_swd,so2_anth_tnr,so2_anth_tro

!gocart2G anthropogenic emissions:
 real(kind=RKIND),dimension(:),pointer:: su_anth_less100m,su_anth_less500m,su_anth_shipso2,su_anth_shipso4, &
                                         su_anth_aviation_lto,su_anth_aviation_cds,su_anth_aviation_crs
 real(kind=RKIND),dimension(:,:),pointer:: su_anth_aircraft

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine init_anth_emissions_su:')


 call mpas_pool_get_dimension(mesh,'nCells',nCells)
 call mpas_pool_get_dimension(mesh,'nVertLevels',nVertLevels)
 call mpas_log_write('    nCells      = $i',intArgs=(/nCells/))
 call mpas_log_write('    nVertLevels = $i',intArgs=(/nVertLevels/))


 allocate(anth_less100m(nCells+1)    )
 allocate(anth_less500m(nCells+1)    )
 allocate(anth_shipso2(nCells+1)     )
 allocate(anth_shipso4(nCells+1)     )
 allocate(anth_aviation_lto(nCells+1))
 allocate(anth_aviation_cds(nCells+1))
 allocate(anth_aviation_crs(nCells+1))
 allocate(anth_aircraft(nVertLevels,1:nCells+1))


!--- SU anthropogenic emissions:
!    here, we assume that anthropogenic emissions from ships and aircrafts are not included. they do not need to
!    be interpolated to the initial time of the forecast and are set to zero. furthermore, we assume that surface
!    emissions are all non-energy related emissions and simply interpolated to layers below 100 meters.

 call mpas_pool_get_array(anth_emissions,'su_anth_less100m'    ,su_anth_less100m    )
 call mpas_pool_get_array(anth_emissions,'su_anth_less500m'    ,su_anth_less500m    )
 call mpas_pool_get_array(anth_emissions,'su_anth_shipso2'     ,su_anth_shipso2     )
 call mpas_pool_get_array(anth_emissions,'su_anth_shipso4'     ,su_anth_shipso4     )
 call mpas_pool_get_array(anth_emissions,'su_anth_aviation_lto',su_anth_aviation_lto)
 call mpas_pool_get_array(anth_emissions,'su_anth_aviation_cds',su_anth_aviation_cds)
 call mpas_pool_get_array(anth_emissions,'su_anth_aviation_crs',su_anth_aviation_crs)
 call mpas_pool_get_array(anth_emissions,'su_anth_aircraft'    ,su_anth_aircraft    )
 su_anth_less100m(:)     = 0._RKIND
 su_anth_less500m(:)     = 0._RKIND
 su_anth_shipso2(:)      = 0._RKIND
 su_anth_shipso4(:)      = 0._RKIND
 su_anth_aviation_lto(:) = 0._RKIND
 su_anth_aviation_cds(:) = 0._RKIND
 su_anth_aviation_crs(:) = 0._RKIND
 su_anth_aircraft(:,:)   = 0._RKIND


 call mpas_pool_get_array(CAMS_anth_emissions,'so2_anth_agl',so2_anth_agl)
 call mpas_pool_get_array(CAMS_anth_emissions,'so2_anth_ags',so2_anth_ags)
 call mpas_pool_get_array(CAMS_anth_emissions,'so2_anth_awb',so2_anth_awb)
 call mpas_pool_get_array(CAMS_anth_emissions,'so2_anth_ene',so2_anth_ene)
 call mpas_pool_get_array(CAMS_anth_emissions,'so2_anth_fef',so2_anth_fef)
 call mpas_pool_get_array(CAMS_anth_emissions,'so2_anth_ind',so2_anth_ind)
 call mpas_pool_get_array(CAMS_anth_emissions,'so2_anth_res',so2_anth_res)
 call mpas_pool_get_array(CAMS_anth_emissions,'so2_anth_shp',so2_anth_shp)
 call mpas_pool_get_array(CAMS_anth_emissions,'so2_anth_slv',so2_anth_slv)
 call mpas_pool_get_array(CAMS_anth_emissions,'so2_anth_sum',so2_anth_sum)
 call mpas_pool_get_array(CAMS_anth_emissions,'so2_anth_swd',so2_anth_swd)
 call mpas_pool_get_array(CAMS_anth_emissions,'so2_anth_tnr',so2_anth_tnr)
 call mpas_pool_get_array(CAMS_anth_emissions,'so2_anth_tro',so2_anth_tro)

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


 deallocate(anth_less100m    )
 deallocate(anth_less500m    )
 deallocate(anth_shipso2     )
 deallocate(anth_shipso4     )
 deallocate(anth_aviation_lto)
 deallocate(anth_aviation_cds)
 deallocate(anth_aviation_crs)
 deallocate(anth_aircraft    )


 call mpas_log_write('--- end subroutine init_anth_emissions_su.')

 end subroutine init_anth_emissions_su

!==================================================================================================================
 subroutine init_anth_emissions_nh3(clock,stream_manager,mesh,CAMS_anth_emissions,anth_emissions)
!==================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in),pointer:: mesh
 type(mpas_pool_type),intent(in),pointer:: CAMS_anth_emissions
 type(mpas_Clock_type),intent(in),pointer:: clock

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

 real(kind=RKIND),dimension(:),pointer:: anth_ag,anth_en,anth_in,anth_oc,anth_re,anth_tr

!CAMS anthropogenic emissions:
 real(kind=RKIND),dimension(:),pointer:: nh3_anth_agl,nh3_anth_ags,nh3_anth_awb,nh3_anth_ene,nh3_anth_fef, &
                                         nh3_anth_ind,nh3_anth_res,nh3_anth_shp,nh3_anth_slv,nh3_anth_sum, &
                                         nh3_anth_swd,nh3_anth_tnr,nh3_anth_tro

!gocart2G anthropogenic emissions:
 real(kind=RKIND),dimension(:),pointer:: nh3_anth_ag,nh3_anth_en,nh3_anth_in,nh3_anth_oc,nh3_anth_re,nh3_anth_tr

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine init_anth_emissions_nh3:')


 call mpas_pool_get_dimension(mesh,'nCells',nCells)
 call mpas_log_write('    nCells = $i',intArgs=(/nCells/))


 allocate(anth_ag(nCells+1))
 allocate(anth_en(nCells+1))
 allocate(anth_in(nCells+1))
 allocate(anth_oc(nCells+1))
 allocate(anth_re(nCells+1))
 allocate(anth_tr(nCells+1))


!--- NH3 anthropogenic emissions:
 call mpas_pool_get_array(anth_emissions,'nh3_anth_ag'    ,nh3_anth_ag)
 call mpas_pool_get_array(anth_emissions,'nh3_anth_en'    ,nh3_anth_en)
 call mpas_pool_get_array(anth_emissions,'nh3_anth_in'    ,nh3_anth_in)
 call mpas_pool_get_array(anth_emissions,'nh3_anth_oc'    ,nh3_anth_oc)
 call mpas_pool_get_array(anth_emissions,'nh3_anth_re'    ,nh3_anth_re)
 call mpas_pool_get_array(anth_emissions,'nh3_anth_tr'    ,nh3_anth_tr)
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


 call mpas_log_write('--- end subroutine init_anth_emissions_nh3.')

 end subroutine init_anth_emissions_nh3

!==================================================================================================================
 subroutine init_BIOG_emissions(clock,stream_manager,mesh,BIOG_emissions)
!==================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in),pointer:: mesh
 type(mpas_Clock_type),intent(in),pointer:: clock

!inout arguments:
 type(MPAS_streamManager_type),intent(inout):: stream_manager
 type(mpas_pool_type),intent(inout),pointer:: BIOG_emissions

!local variables and arrays:
 type(mpas_time_type):: beforeTime,afterTime,currTime
 type(mpas_timeinterval_type):: beforeDelta,afterDelta,totalDelta

 character(len=StrKIND):: actualTimeStamp
 integer,pointer:: nCells
 real(kind=RKIND):: total_dt,before_dt,after_dt

 real(kind=RKIND),dimension(:),pointer:: bc_before,br_before,oc_before
 real(kind=RKIND),dimension(:),pointer:: bc,br,oc

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine init_BIOG_emissions:')

 call mpas_log_write('--- end subroutine init_BIOG_emissions.')

 end subroutine init_BIOG_emissions

!==================================================================================================================
 subroutine init_biomass_burning_emissions(clock,stream_manager,mesh,FINN_biob_emissions,biob_emissions)
!==================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in),pointer:: mesh
 type(mpas_pool_type),intent(in),pointer:: FINN_biob_emissions
 type(mpas_Clock_type),intent(in),pointer:: clock

!inout arguments:
 type(MPAS_streamManager_type),intent(inout):: stream_manager
 type(mpas_pool_type),intent(inout),pointer:: biob_emissions

!local variables and arrays:
 type(mpas_time_type):: beforeTime,afterTime,currTime
 type(mpas_timeinterval_type):: beforeDelta,afterDelta,totalDelta

 character(len=StrKIND):: actualTimeStamp
 integer:: iCell
 integer,pointer:: nCells
 real(kind=RKIND):: total_dt,before_dt,after_dt

 real(kind=RKIND),dimension(:),pointer:: biob_modis

!FINN biomass burning emissions:
 real(kind=RKIND),dimension(:),pointer:: bc_biob_modis
 real(kind=RKIND),dimension(:),pointer:: oc_biob_modis
 real(kind=RKIND),dimension(:),pointer:: nh3_biob_modis
 real(kind=RKIND),dimension(:),pointer:: so2_biob_modis

!gocart2G biomass burning emissions:
 real(kind=RKIND),dimension(:),pointer:: bc_biob_em
 real(kind=RKIND),dimension(:),pointer:: br_biob_em
 real(kind=RKIND),dimension(:),pointer:: oc_biob_em
 real(kind=RKIND),dimension(:),pointer:: ni_biob_em
 real(kind=RKIND),dimension(:),pointer:: su_biob_em

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine init_biomass_burning_emissions:')


 call mpas_pool_get_dimension(mesh,'nCells',nCells)
 call mpas_log_write('    nCells = $i',intArgs=(/nCells/))


 allocate(biob_modis(nCells+1))


!--- FINN biomass burning emissions:
 call mpas_pool_get_array(biob_emissions,'bc_biob_em',bc_biob_em)
 call mpas_pool_get_array(biob_emissions,'br_biob_em',br_biob_em)
 call mpas_pool_get_array(biob_emissions,'oc_biob_em',oc_biob_em)
 call mpas_pool_get_array(biob_emissions,'ni_biob_em',ni_biob_em)
 call mpas_pool_get_array(biob_emissions,'su_biob_em',su_biob_em)
 bc_biob_em(:) = 0._RKIND
 br_biob_em(:) = 0._RKIND
 oc_biob_em(:) = 0._RKIND
 ni_biob_em(:) = 0._RKIND
 su_biob_em(:) = 0._RKIND


 call mpas_pool_get_array(FINN_biob_emissions,'bc_biob_modis' ,bc_biob_modis )
 call mpas_pool_get_array(FINN_biob_emissions,'oc_biob_modis' ,oc_biob_modis )
 call mpas_pool_get_array(FINN_biob_emissions,'nh3_biob_modis',nh3_biob_modis)
 call mpas_pool_get_array(FINN_biob_emissions,'so2_biob_modis',so2_biob_modis)


!biomass burning of black carbon:
 call mpas_stream_mgr_read(stream_manager,'biob_bc_emissions',rightNow=.true., &
                  whence=MPAS_STREAM_LATEST_BEFORE,actualWhen=actualTimestamp)
 call mpas_log_write('    latest time before is  '//trim(actualTimestamp))
 call mpas_set_time(beforeTime,dateTimeString=trim(actualTimestamp))

 biob_modis(1:nCells) = bc_biob_modis(1:nCells)

 call mpas_stream_mgr_read(stream_manager,'biob_bc_emissions',rightNow=.true., &
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

!interpolation of black carbon biomass burning emissions to the current time:
 if(total_dt > 0.0_RKIND) then
    bc_biob_em(:) = (after_dt/total_dt)*biob_modis(:) + (before_dt/total_dt)*bc_biob_modis(:)
 endif


!biomass burning of organic carbon:
 call mpas_stream_mgr_read(stream_manager,'biob_oc_emissions',rightNow=.true., &
                  whence=MPAS_STREAM_LATEST_BEFORE,actualWhen=actualTimestamp)
 call mpas_set_time(beforeTime,dateTimeString=trim(actualTimestamp))

 biob_modis(1:nCells) = oc_biob_modis(1:nCells)

 call mpas_stream_mgr_read(stream_manager,'biob_oc_emissions',rightNow=.true., &
                  whence=MPAS_STREAM_EARLIEST_AFTER,actualWhen=actualTimestamp)
 call mpas_set_time(afterTime,dateTimeString=trim(actualTimestamp))

!interpolation of organic carbon biomass burning emissions to the current time:
 if(total_dt > 0.0_RKIND) then
    oc_biob_em(:) = (after_dt/total_dt)*biob_modis(:) + (before_dt/total_dt)*oc_biob_modis(:)
 endif


!biomass burning of nitrate:
 call mpas_stream_mgr_read(stream_manager,'biob_ni_emissions',rightNow=.true., &
                  whence=MPAS_STREAM_LATEST_BEFORE,actualWhen=actualTimestamp)
 call mpas_set_time(beforeTime,dateTimeString=trim(actualTimestamp))

 biob_modis(1:nCells) = nh3_biob_modis(1:nCells)

 call mpas_stream_mgr_read(stream_manager,'biob_ni_emissions',rightNow=.true., &
                  whence=MPAS_STREAM_EARLIEST_AFTER,actualWhen=actualTimestamp)
 call mpas_set_time(afterTime,dateTimeString=trim(actualTimestamp))

!interpolation of organic carbon biomass burning emissions to the current time:
 if(total_dt > 0.0_RKIND) then
    ni_biob_em(:) = (after_dt/total_dt)*biob_modis(:) + (before_dt/total_dt)*nh3_biob_modis(:)
 endif


!biomass burning of sulphur dioxide:
 call mpas_stream_mgr_read(stream_manager,'biob_su_emissions',rightNow=.true., &
                  whence=MPAS_STREAM_LATEST_BEFORE,actualWhen=actualTimestamp)
 call mpas_set_time(beforeTime,dateTimeString=trim(actualTimestamp))

 biob_modis(1:nCells) = so2_biob_modis(1:nCells)

 call mpas_stream_mgr_read(stream_manager,'biob_su_emissions',rightNow=.true., &
                  whence=MPAS_STREAM_EARLIEST_AFTER,actualWhen=actualTimestamp)
 call mpas_set_time(afterTime,dateTimeString=trim(actualTimestamp))

!interpolation of organic carbon biomass burning emissions to the current time:
 if(total_dt > 0.0_RKIND) then
    su_biob_em(:) = (after_dt/total_dt)*biob_modis(:) + (before_dt/total_dt)*so2_biob_modis(:)
 endif


 deallocate(biob_modis)


 call mpas_log_write('--- end subroutine init_biomass_burning_emissions.')

 end subroutine init_biomass_burning_emissions

!==================================================================================================================
 end module mpas_chemistry_init_gocart2G_emissions
!==================================================================================================================
