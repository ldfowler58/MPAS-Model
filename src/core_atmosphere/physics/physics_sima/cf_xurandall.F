!=================================================================================================================
! cf_xurandall     : calculates the cloud fraction as a function of the relative humidity (Xu-Randall scheme)
!
! cf_xurandall.F created  Jihyeon Jang / 2020-03-18
!=================================================================================================================
!
module cf_xurandall
  use ccpp_kinds,only: kind_phys

  IMPLICIT NONE

  PUBLIC :: cf_xurandall_run
  PUBLIC :: cf_xurandall_init
  PUBLIC :: cf_xurandall_final
  PUBLIC :: cf_xurandall_timestep_init
  PUBLIC :: cf_xurandall_timestep_final

contains
!
!-------------------------------------------------------------------------------
!
   subroutine cf_xurandall_run(cldfrac,t_p,pres_p,qv,qc,qi,qs,                 &
                               svp1, svp2, svpi2, svp3,                        &
                               svpi3, svpt0, ep_2,                             &
                               its,ite,kte,                                    &
                               errmsg,errflg                                   &
                               )
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
!  parameters:
!
   real(kind=kind_phys),parameter    ::  alpha0  = 100.
   real(kind=kind_phys),parameter    ::  gamma   = 0.49
   real(kind=kind_phys),parameter    ::  qcldmin = 1.e-12
   real(kind=kind_phys),parameter    ::  pexp    = 0.25
   real(kind=kind_phys),parameter    ::  rhgrid  = 1.0
!
   integer,parameter :: kts=1
!
!  input arguments:
!
   integer,  intent(in)      :: its,ite,kte
!
   real(kind=kind_phys),                                                       &
             intent(in)      ::                                          svp1, &
                                                                         svp2, &
                                                                        svpi2, &
                                                                         svp3, &
                                                                        svpi3, &
                                                                        svpt0, &
                                                                         ep_2
!
   real(kind=kind_phys),     dimension( its:ite, kts:kte ),                    &
             intent(in)      ::                                            qv, &
                                                                           qc, &
                                                                           qi, &
                                                                           qs
!
   real(kind=kind_phys),     dimension( its:ite, kts:kte ),                    &
             intent(in)      ::                                           t_p, &
                                                                       pres_p
!
!  output arguments:
!
   real(kind=kind_phys),     dimension( its:ite, kts:kte ),                    &
             intent(out)     ::                                       cldfrac
!
   character(len=*), intent(out)   ::                                  errmsg
   integer,          intent(out)   ::                                  errflg
!
!  local variables:
!
   integer           :: i,k
   real(kind=kind_phys)              :: esi,esw,qvsi,qvsw
   real(kind=kind_phys)              :: arg,denom,qcld,qvs,rhum,subsat,weight
!
!-------------------------------------------------------------------------------
!
   do k = kts,kte
     do i = its,ite

!... calculation of the saturation mixing ratios over water and over ice (Murray, 1966):
       esw = 1000. * svp1 * exp(svp2 * (t_p(i,k) - svpt0) / (t_p(i,k) - svp3))
       esi = 1000. * svp1 * exp(svpi2 * (t_p(i,k) - svpt0) / (t_p(i,k) - svpi3))

       qvsw = ep_2 * esw / (pres_p(i,k) - esw)
       qvsi = ep_2 * esi / (pres_p(i,k) - esi)

       qcld = qc(i,k) + qi(i,k) + qs(i,k)
       if(qcld .lt. qcldmin) then
         weight = 0.
       else
         weight = (qi(i,k) + qs(i,k)) / qcld
       endif

       qvs = (1-weight) * qvsw + weight * qvsi
       rhum = qv(i,k) / qvs

       if(qcld .lt. qcldmin) then
         !assume that the cloud fraction is equal to 0. when the cloudy mixing ratio equals 0.
         cldfrac(i,k) = 0.

       elseif(rhum .ge. rhgrid) then
         !assume that the cloud fraction is equal to 1. when the relative humidity equal 100%.
         cldfrac(i,k) = 1.

       else
         !computation of the cloud fraction:
         subsat = max(1.e-10,rhgrid*qvs-qv(i,k))
         denom  = subsat**gamma
         arg    = max(-6.9,-alpha0*qcld/denom) ! exp(-6.9) = 0.001

         rhum = max(1.e-10,rhum)
         cldfrac(i,k) = (rhum/rhgrid)**pexp*(1.-exp(arg))
         if(cldfrac(i,k) .lt. 0.01) cldfrac(i,k) = 0.

       endif

     enddo
   enddo
!
   errmsg = 'cf_xurandall_run OK'
   errflg = 0
!
   end subroutine cf_xurandall_run
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   subroutine cf_xurandall_init (errmsg, errflg)

    character(len=*),        intent(out)   :: errmsg
    integer,                 intent(out)   :: errflg

    ! This routine currently does nothing

    errmsg = ''
    errflg = 0

   end subroutine cf_xurandall_init
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   subroutine cf_xurandall_final (errmsg, errflg)

    character(len=*),        intent(out)   :: errmsg
    integer,                 intent(out)   :: errflg

    ! This routine currently does nothing

    errmsg = ''
    errflg = 0

   end subroutine cf_xurandall_final
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine cf_xurandall_timestep_init (errmsg, errflg)

    character(len=*),        intent(out)   :: errmsg
    integer,                 intent(out)   :: errflg

    ! This routine currently does nothing

    errmsg = ''
    errflg = 0

   end subroutine cf_xurandall_timestep_init
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   subroutine cf_xurandall_timestep_final (errmsg, errflg)

    character(len=*),        intent(out)   :: errmsg
    integer,                 intent(out)   :: errflg

    ! This routine currently does nothing

    errmsg = ''
    errflg = 0

   end subroutine cf_xurandall_timestep_final
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
 end module cf_xurandall
!=================================================================================================================
