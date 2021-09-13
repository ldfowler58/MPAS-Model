!======================================================================================================================
 module pumas_mpas
!======================================================================================================================
 use mpas_kind_types
 use mpas_log

 use pumas_vars,only: mgncol,nlev,cpair,gravit,latice,latvap,rair,rh2o,tmelt


 implicit none
 private
 public:: pumas_mpas_init,          &
          pumas_mpas_timestep_init, &
          pumas_mpas_timestep_final


 contains


!======================================================================================================================
 subroutine pumas_mpas_init(ncells,nlevels,con_cp,con_grav,con_rd,con_rv,con_tmelt,con_xlf,con_xlv,errmsg,errflg)
!======================================================================================================================

!input arguments:
 integer,intent(in):: &
      ncells,    &
      nlevels

 real(kind=RKIND),intent(in):: &
      con_cp,    &
      con_grav,  &
      con_rd,    &
      con_rv,    &
      con_tmelt, &
      con_xlf,   &
      con_xlv


!output arguments:
 character(len=*),intent(out):: errmsg
 integer,intent(out):: errflg

!----------------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine mp_pumas_init:')

 mgncol = ncells
 nlev   = nlevels

 cpair  = con_cp
 gravit = con_grav
 rair   = con_rd
 rh2o   = con_rv
 tmelt  = con_tmelt
 latvap = con_xlv
 latice = con_xlf
 
 errmsg = ' '
 errflg = 0

 call mpas_log_write('--- end subroutine mp_pumas_init:')

 end subroutine pumas_mpas_init

!======================================================================================================================
 subroutine pumas_mpas_run()
!======================================================================================================================

!----------------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
!call mpas_log_write('--- enter subroutine mp_pumas_run:')

 end subroutine pumas_mpas_run

!======================================================================================================================
 subroutine pumas_mpas_timestep_init()
!======================================================================================================================
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine pumas_mpas_timestep_init:')

 

 call mpas_log_write('--- end subroutine pumas_mpas_timestep_init:')

 end subroutine pumas_mpas_timestep_init

!======================================================================================================================
 subroutine pumas_mpas_timestep_final()
!======================================================================================================================

 end subroutine pumas_mpas_timestep_final

!======================================================================================================================
 end module pumas_mpas
!======================================================================================================================
