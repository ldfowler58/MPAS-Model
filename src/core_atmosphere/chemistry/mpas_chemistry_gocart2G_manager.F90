! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================================
 module mpas_chemistry_gocart2G_manager
 use mpas_kind_types
 use mpas_pool_routines,only: mpas_pool_get_config,mpas_pool_get_subpool
 use mpas_timekeeping
 use mpas_stream_manager

 use mpas_chemistry_gocart2G_emissions_init
 use mpas_chemistry_gocart2G_emissions_update
 use mpas_chemistry_gocart2G_update


 implicit none
 private
 public:: gocart2G_timetracker,init_gocart2G_timetracker

 character(len=*),parameter:: gocart2GAlarmID     = 'gocart2G'
 character(len=*),parameter:: gocart2GanthAlarmID = 'gocart2G_anth'
 character(len=*),parameter:: gocart2GbiobAlarmID = 'gocart2G_biob'
 character(len=*),parameter:: gocart2GbiogAlarmID = 'gocart2G_biog'

 integer,public:: iyear
 integer,public:: imonth
 integer,public:: iday
 integer,public:: ihour
 integer,public:: iminute
 integer,public:: isecond


 contains


!==================================================================================================================
 subroutine gocart2G_timetracker(domain,clock,stream_manager)
!==================================================================================================================

!--- inout arguments:
 type(MPAS_Clock_type),intent(inout):: clock
 type(domain_type),intent(inout):: domain
 type(MPAS_streamManager_type),intent(inout):: stream_manager

!--- local variables:
 type(block_type),pointer:: block
 type(MPAS_Time_Type):: currTime
 type(mpas_pool_type),pointer:: mesh
 type(mpas_pool_type),pointer:: gocart2G_met
 type(mpas_pool_type),pointer:: gocart2G_backgrounds
 type(mpas_pool_type),pointer:: anth_emissions
 type(mpas_pool_type),pointer:: biob_emissions
 type(mpas_pool_type),pointer:: biog_emissions
 type(mpas_pool_type),pointer:: CAMS_anth_emissions
 type(mpas_pool_type),pointer:: CAMS_biog_emissions
 type(mpas_pool_type),pointer:: FINN_biob_emissions

 character(len=StrKIND):: timeStamp
 character(len=StrKIND),pointer:: backg_interval
 character(len=StrKIND),pointer:: anth_interval
 character(len=StrKIND),pointer:: biob_interval
 character(len=StrKIND),pointer:: biog_interval

 integer:: ierr

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine gocart2G_timetracker:')

 currTime = mpas_get_clock_time(clock,MPAS_NOW,ierr)
 call mpas_get_time(curr_time=currTime,dateTimeString=timeStamp,YYYY=iyear,MM=imonth, &
                    DD=iday,H=ihour,M=iminute,S=isecond,ierr=ierr)
 call mpas_log_write('--- YEAR   = $i',intArgs=(/iyear/))
 call mpas_log_write('--- MONTH  = $i',intArgs=(/imonth/))
 call mpas_log_write('--- DAY    = $i',intArgs=(/iday/))
 call mpas_log_write('--- HOUR   = $i',intArgs=(/ihour/))
 call mpas_log_write('--- MINUTE = $i',intArgs=(/iminute/))
 call mpas_log_write('--- SECOND = $i',intArgs=(/isecond/))


!--- check to see if it is time to update the background climatological fields:
 block => domain%blocklist
 do while(associated(block))

    call mpas_pool_get_subpool(block%structs,'mesh'                ,mesh                )
    call mpas_pool_get_subpool(block%structs,'gocart2G_met     '   ,gocart2G_met        )
    call mpas_pool_get_subpool(block%structs,'gocart2G_backgrounds',gocart2G_backgrounds)

    call mpas_pool_get_config(domain%blocklist%configs,'config_gocart2G_backg_interval',backg_interval)
    if(mpas_is_alarm_ringing(clock,gocart2GAlarmID,ierr=ierr)) then
       call mpas_reset_clock_alarm(clock,gocart2GAlarmID,ierr=ierr)
       call mpas_log_write('--- time to update gocart2G climatology:')
       call update_gocart2G_climatology(timeStamp,mesh,gocart2G_met,gocart2G_backgrounds)
    endif

    block => block%next

 end do


!--- check to see if it is time to update anthropogenic, biomass burning, and biogenic emissions:
 block => domain%blocklist
 do while(associated(block))

    call mpas_pool_get_config(domain%blocklist%configs,'config_gocart2G_anthemis_interval',anth_interval)
    call mpas_pool_get_config(domain%blocklist%configs,'config_gocart2G_biobemis_interval',biob_interval)
    call mpas_pool_get_config(domain%blocklist%configs,'config_gocart2G_biogemis_interval',biog_interval)

    call mpas_log_write(' ')
    call mpas_log_write('--- config_gocart2G_anthemis_interval = '//anth_interval)
    call mpas_log_write('--- config_gocart2G_biobemis_interval = '//biob_interval)
    call mpas_log_write('--- config_gocart2G_biogemis_interval = '//biog_interval)

    call mpas_pool_get_subpool(block%structs,'mesh'                ,mesh              )
    call mpas_pool_get_subpool(block%structs,'CAMS_anth_emissions',CAMS_anth_emissions)
    call mpas_pool_get_subpool(block%structs,'CAMS_biog_emissions',CAMS_biog_emissions)
    call mpas_pool_get_subpool(block%structs,'FINN_biob_emissions',FINN_biob_emissions)
    call mpas_pool_get_subpool(block%structs,'anth_emissions'     ,anth_emissions     )
    call mpas_pool_get_subpool(block%structs,'biob_emissions'     ,biob_emissions     )
    call mpas_pool_get_subpool(block%structs,'biog_emissions'     ,biog_emissions     )

    if(mpas_is_alarm_ringing(clock,gocart2GanthAlarmID,ierr=ierr)) then
       call mpas_reset_clock_alarm(clock,gocart2GanthAlarmID,ierr=ierr)
       call mpas_log_write(' ')
       call mpas_log_write('--- time to update gocart2G anthropogenic surface emissions:')
       call update_anth_emissions_bc(clock,stream_manager,mesh,CAMS_anth_emissions,anth_emissions)
       call update_anth_emissions_oc(clock,stream_manager,mesh,CAMS_anth_emissions,anth_emissions)
       call update_anth_emissions_nh3(clock,stream_manager,mesh,CAMS_anth_emissions,anth_emissions)
       call update_anth_emissions_su(clock,stream_manager,mesh,CAMS_anth_emissions,anth_emissions)
       call update_anth_emissions_iso(clock,stream_manager,mesh,CAMS_anth_emissions,anth_emissions)
       call update_anth_emissions_mnt(clock,stream_manager,mesh,CAMS_anth_emissions,anth_emissions)
    endif

    if(mpas_is_alarm_ringing(clock,gocart2GbiobAlarmID,ierr=ierr)) then
       call mpas_reset_clock_alarm(clock,gocart2GbiobAlarmID,ierr=ierr)
       call mpas_log_write(' ')
       call mpas_log_write('--- time to update gocart2G biomass burning surface emissions:')
       call update_biomass_burning_emissions(clock,stream_manager,mesh,FINN_biob_emissions,biob_emissions)
    endif

    if(mpas_is_alarm_ringing(clock,gocart2GbiogAlarmID,ierr=ierr)) then
       call mpas_reset_clock_alarm(clock,gocart2GbiogAlarmID,ierr=ierr)
       call mpas_log_write(' ')
       call mpas_log_write('--- time to update gocart2G biogenic surface emissions:')
       call update_biog_emissions(clock,stream_manager,mesh,CAMS_biog_emissions,biog_emissions)
    endif

    !--- updates local arrays:
    call init_gocart2G_emissions(mesh,anth_emissions,biob_emissions,biog_emissions)

    block => block%next

 end do


 call mpas_log_write('--- end subroutine gocart2G_timetracker.')
 call mpas_log_write(' ')

 end subroutine gocart2G_timetracker

!==================================================================================================================
 subroutine init_gocart2G_timetracker(configs,clock)
!==================================================================================================================

!--- input arguments:
 type(mpas_pool_type),intent(in):: configs

!--- local variables:
 type(MPAS_Clock_type),intent(inout):: clock
 type(MPAS_Time_Type):: startTime,alarmStartTime
 type(MPAS_TimeInterval_Type):: alarmTimeStep

 character(len=StrKIND),pointer:: backg_interval
 character(len=StrKIND),pointer:: anth_interval
 character(len=StrKIND),pointer:: biob_interval
 character(len=StrKIND),pointer:: biog_interval

 integer:: ierr

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine gocart2G_timetracker_init:')

!--- computes startTime:
 ierr = 0
 startTime = mpas_get_clock_time(clock,MPAS_START_TIME,ierr)
 if(ierr /=0) &
    call mpas_log_write('--- gocart2G_timetracker_init: error getting startTime', &
                        messageType=MPAS_LOG_CRIT)

 ierr = 0
 call mpas_get_time(curr_time=startTime,YYYY=iyear,MM=imonth,DD=iday,H=ihour,M=iminute,S=isecond,ierr=ierr)
 if(ierr /=0) then
    call mpas_log_write('--- gocart2G_timetracker_init: error deriving iyear,imonth,etc...', &
                        messageType=MPAS_LOG_CRIT)
 else
    call mpas_log_write('--- YEAR   = $i',intArgs=(/iyear/))
    call mpas_log_write('--- MONTH  = $i',intArgs=(/imonth/))
    call mpas_log_write('--- DAY    = $i',intArgs=(/iday/))
    call mpas_log_write('--- HOUR   = $i',intArgs=(/ihour/))
    call mpas_log_write('--- MINUTE = $i',intArgs=(/iminute/))
    call mpas_log_write('--- SECOND = $i',intArgs=(/isecond/))
 endif


!--- initializes alarm to update the background climatological fields:
 call mpas_pool_get_config(configs,'config_gocart2G_backg_interval',backg_interval)

 ierr = 0
 call mpas_set_timeInterval(alarmTimeStep,timeString=backg_interval,ierr=ierr)
 alarmStartTime = startTime
 call mpas_add_clock_alarm(clock,gocart2GAlarmID,alarmStartTime,alarmTimeStep,ierr=ierr)
 if(ierr /= 0) &
 call mpas_log_write('--- gocart2G_timetracker_init: error creating gocart2GAlarmID', &
                          messageType=MPAS_LOG_CRIT)


!--- initializes alarm to update anthropogenic, biomass burning, and biogenic surface emissions:
 call mpas_pool_get_config(configs,'config_gocart2G_anthemis_interval',anth_interval)
 call mpas_pool_get_config(configs,'config_gocart2G_biobemis_interval',biob_interval)
 call mpas_pool_get_config(configs,'config_gocart2G_biogemis_interval',biog_interval)

 ierr = 0
 call mpas_set_timeInterval(alarmTimeStep,timeString=anth_interval,ierr=ierr)
 alarmStartTime = startTime
 call mpas_add_clock_alarm(clock,gocart2GanthAlarmID,alarmStartTime,alarmTimeStep,ierr=ierr)
 if(ierr /= 0) &
 call mpas_log_write('--- gocart2G_timetracker_init: error creating gocart2GanthAlarmID', &
                          messageType=MPAS_LOG_CRIT)


 ierr = 0
 call mpas_set_timeInterval(alarmTimeStep,timeString=biob_interval,ierr=ierr)
 alarmStartTime = startTime
 call mpas_add_clock_alarm(clock,gocart2GbiobAlarmID,alarmStartTime,alarmTimeStep,ierr=ierr)
 if(ierr /= 0) &
 call mpas_log_write('--- gocart2G_timetracker_init: error creating gocart2GbiobAlarmID', &
                          messageType=MPAS_LOG_CRIT)


 ierr = 0
 call mpas_set_timeInterval(alarmTimeStep,timeString=biog_interval,ierr=ierr)
 alarmStartTime = startTime
 call mpas_add_clock_alarm(clock,gocart2GbiogAlarmID,alarmStartTime,alarmTimeStep,ierr=ierr)
 if(ierr /= 0) &
 call mpas_log_write('--- gocart2G_timetracker_init: error creating gocart2GbiogAlarmID', &
                          messageType=MPAS_LOG_CRIT)


 call mpas_log_write('--- end subroutine gocart2G_timetracker_init.')

 end subroutine init_gocart2G_timetracker

!==================================================================================================================
 end module mpas_chemistry_gocart2G_manager
!==================================================================================================================
