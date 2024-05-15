!=================================================================================================================
 module mpas_chemistry_manager
 use mpas_kind_types
 use mpas_timekeeping


 implicit none
 private
 public:: chemistry_timetracker

 integer,public:: iyear
 integer,public:: imonth
 integer,public:: iday
 integer,public:: ihour
 integer,public:: iminute
 integer,public:: isecond


 contains


!=================================================================================================================
 subroutine chemistry_timetracker(clock)
!=================================================================================================================

!input arguments:
 type(MPAS_Clock_type),intent(in):: clock

!local variables:
 type(MPAS_Time_Type):: currTime

 character(len=StrKIND):: timeStamp
 integer:: ierr

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine chemistry_timetracker:')

 currTime = mpas_get_clock_time(clock,MPAS_NOW,ierr)
 call mpas_get_time(curr_time=currTime,dateTimeString=timeStamp,YYYY=iyear,MM=imonth, &
                    DD=iday,H=ihour,M=iminute,S=isecond,ierr=ierr)
 call mpas_log_write('--- YEAR   = $i',intArgs=(/iyear/))
 call mpas_log_write('--- MONTH  = $i',intArgs=(/imonth/))
 call mpas_log_write('--- DAY    = $i',intArgs=(/iday/))
 call mpas_log_write('--- HOUR   = $i',intArgs=(/ihour/))
 call mpas_log_write('--- MINUTE = $i',intArgs=(/iminute/))
 call mpas_log_write('--- SECOND = $i',intArgs=(/isecond/))

 call mpas_log_write('--- end subroutine chemistry_timetracker:')
 call mpas_log_write(' ')

 end subroutine chemistry_timetracker

!=================================================================================================================
 end module mpas_chemistry_manager
!=================================================================================================================
