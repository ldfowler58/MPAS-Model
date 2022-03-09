!=================================================================================================================
 module mpas_chemistry_manager
 use mpas_log
 use mpas_kind_types
 use mpas_derived_types
 use mpas_pool_routines
 use mpas_timekeeping

 use mpas_atmphys_utilities


 implicit none
 private
 public:: chemistry_run_init,chemistry_timetracker


!defines alarm to update the coefficients used in the CHEM2D Ozone Photochemistry Parameterization (CHEM2D-OPP)
!of McCormack et al.,2006).
 character(len=*), parameter:: oppAlarmID = 'opp'


 contains


!=================================================================================================================
 subroutine chemistry_timetracker(domain,clock)
!=================================================================================================================

!inpout arguments:
 type(MPAS_Clock_type),intent(inout):: clock
 type(domain_type),intent(inout)    :: domain


!local pointers:
 type(block_type),pointer    :: block
 type(mpas_pool_type),pointer:: mesh
 type(mpas_pool_type),pointer:: chem_input


!local variables:
 type(MPAS_Time_Type):: currTime

 logical:: Leapyear
 character(len=StrKIND):: timeStamp
 integer:: h,m,s,s_n,s_d,DoY,yr
 integer:: ierr,julday,year

 real(kind=RKIND):: curr_julday,gmt
 real(kind=RKIND):: utc_h,utc_s

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine chemistry_timetracker:')


!--- update the current julian day and current year:
 currTime = mpas_get_clock_time(clock,MPAS_NOW,ierr)
 call mpas_get_time(curr_time=currTime,dateTimeString=timeStamp,YYYY=yr,H=h,M=m, &
                    S=s,S_n=s_n,S_d=s_d,DoY=DoY,ierr=ierr)

 utc_h = real(h) + real(m) / 60.0 + real(s + s_n / s_d) / 3600.0
 utc_s = real(s + s_n + s_d)
 gmt   = utc_h
 year   = yr
 julday = DoY
 curr_julday = real(julday-1) + utc_h / 24.0
 LeapYear = isLeapYear(year)
 call mpas_log_write('     YEAR        = $i', intArgs=(/year/))
 call mpas_log_write('     JULDAY      = $i', intArgs=(/julday/))
 call mpas_log_write('     GMT         = $r', realArgs=(/gmt/))
 call mpas_log_write('     UTC_H       = $r', realArgs=(/utc_h/))
 call mpas_log_write('     CURR_JULDAY = $r', realArgs=(/curr_julday/))
 call mpas_log_write('     LEAP_YEAR   = $l', logicArgs=(/LeapYear/))
 call mpas_log_write('     TIME STAMP  = '//trim(timeStamp))


 block => domain%blocklist
 do while(associated(block))
    call mpas_pool_get_subpool(block%structs,'mesh'      ,mesh      )
    call mpas_pool_get_subpool(block%structs,'chem_input',chem_input)

    !update the coefficients for the GFS linearized CHEM2D-OPP (Mc Cormack et al. 2006):
!   if(mpas_is_alarm_ringing(clock,oppAlarmID,ierr=ierr)) then
!      call mpas_reset_clock_alarm(clock,greenAlarmID,ierr=ierr)
       call mpas_log_write('--- time to update CHEM2D-OPP coefficients.')
       call chemistry_update_opp(julday,curr_julday,mesh,chem_input)
!   endif

    block => block%next
 end do


 call mpas_log_write('--- end subroutine chemistry_timetracker:')

 end subroutine chemistry_timetracker

!=================================================================================================================
 subroutine chemistry_update_opp(julian,curr_julday,mesh,chem_input)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 integer,intent(in):: julian
 real(kind=RKIND):: curr_julday

!inout arguments:
 type(mpas_pool_type),intent(inout):: chem_input

!local pointers:
 logical:: finddate

 integer,pointer:: nCells,nO3Coeffs,nO3Levels,nMonths
 integer:: iCell,imonth,k,n

 real(kind=RKIND),dimension(:),pointer      :: O3pres
 real(kind=RKIND),dimension(:,:,:),pointer  :: O3coeffs
 real(kind=RKIND),dimension(:,:,:,:),pointer:: O3data

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine update_o3chem2dopp')

 call mpas_pool_get_dimension(mesh,'nCells'    ,nCells   )
 call mpas_pool_get_dimension(mesh,'nMonths'   ,nMonths  )
 call mpas_pool_get_dimension(mesh,'nOPPLevels',nO3Levels)
 call mpas_pool_get_dimension(mesh,'nOPPCoeffs',nO3Coeffs)

 call mpas_pool_get_array(chem_input,'OPPpres'  ,O3pres)
 call mpas_pool_get_array(chem_input,'OPPdata'  ,O3data)
 call mpas_pool_get_array(chem_input,'OPPcoeffs',O3coeffs)


!--- for now, fix the array GFSO3coeffs to GFSO3data for July:
 imonth = 7
 O3coeffs(:,:,:) = O3data(imonth,:,:,:)

!do iCell = 1,nCells
!   do k = 1,nO3Levels
!      call mpas_log_write('$i $i $r $r $r $r $r $r $r',intArgs=(/iCell,k/), &
!                          realArgs=(/O3pres(k),(O3coeffs(n,k,iCell),n=1,nO3Coeffs)/))
!   enddo
!   call mpas_log_write(' ')
!enddo

 call mpas_log_write('--- end subroutine update_o3chem2dopp')

 end subroutine chemistry_update_opp

!=================================================================================================================
 subroutine chemistry_run_init(configs,clock)
!=================================================================================================================

!input arguments:
 type(MPAS_Clock_type):: clock
 type(mpas_pool_type),intent(in):: configs


!local pointers:
 character(len=StrKIND),pointer:: config_OPP_update


!local variables:
 type(MPAS_Time_Type):: startTime,alarmStartTime
 type(MPAS_TimeInterval_Type):: alarmTimeStep

 integer:: ierr

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine chemistry_run_init:')

 call mpas_pool_get_config(configs,'config_OPP_update',config_OPP_update)


!--- get startTime prior to initializating alarms:
 startTime = mpas_get_clock_time(clock,MPAS_START_TIME,ierr)


!--- set alarm for updating the coefficients in the GFS linearized CHEM2D-OPP:
 call mpas_set_timeInterval(alarmTimeStep,timeString=config_OPP_update,ierr=ierr)
 alarmStartTime = startTime
 call mpas_add_clock_alarm(clock,oppAlarmID,alarmStartTime,alarmTimeStep,ierr=ierr)
    if(ierr /= 0) &
       call physics_error_fatal('subroutine chemistry_init: error creating alarm GFS CHEM2D-OPP')

 call mpas_log_write('--- end subroutine chemistry_run_init:')

 end subroutine chemistry_run_init

!=================================================================================================================
 end module mpas_chemistry_manager
!=================================================================================================================


