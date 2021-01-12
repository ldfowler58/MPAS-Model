#define NEED_B4B_DURING_CCPP_TESTING 1
!=================================================================================================================
! module_cld_xurandall.F created  Jihyeon Jang / 2020-03-18
!=================================================================================================================
!

module module_cf_xurandall

#if (NEED_B4B_DURING_CCPP_TESTING == 1)
  use mpas_atmphys_constants, only: ep_2
#elif (NEED_B4B_DURING_CCPP_TESTING != 1)
  use mpas_atmphys_constants, only: ep_2, svp1, svp2, svpt0
#endif
  use mpas_log
  use cf_xurandall, only: cf_xurandall_run

  IMPLICIT NONE

  PUBLIC :: xurandall

contains
!
!-------------------------------------------------------------------------------
   subroutine xurandall(cldfrac_3d,t_p_3d,pres_p_3d,                           &
                        qv_3d,qc_3d,qi_3d,qs_3d,                               &
                        ids,ide, jds,jde, kds,kde,                             &
                        ims,ime, jms,jme, kms,kme,                             &
                        its,ite, jts,jte, kts,kte,                             &
                        errmsg,errflg                                          &
                       )
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
!  input arguments:
!
   integer,  intent(in)      ::      ids,ide, jds,jde, kds,kde,                &
                                     ims,ime, jms,jme, kms,kme,                &
                                     its,ite, jts,jte, kts,kte

   real,     dimension( ims:ime, kms:kme, jms:jme )                          , &
             intent(in)      ::                                         qv_3d, &
                                                                        qc_3d, &
                                                                        qi_3d, &
                                                                        qs_3d
!
   real,     dimension( ims:ime, kms:kme, jms:jme )                          , &
             intent(in)      ::                                        t_p_3d, &
                                                                    pres_p_3d
!
!  output arguments:
!
   real,     dimension( ims:ime, kms:kme, jms:jme )                          , &
             intent(out)     ::                                    cldfrac_3d
!
   character(len=*), intent(out)   ::                                  errmsg
   integer,          intent(out)   ::                                  errflg
!
!  local variables:
!
   real,     dimension( its:ite, kts:kte )                 ::                  &
                                                                           qv, &
                                                                           qc, &
                                                                           qi, &
                                                                           qs, &
                                                                          t_p, &
                                                                       pres_p, &
                                                                      cldfrac
!
#if (NEED_B4B_DURING_CCPP_TESTING == 1)
   real,parameter    ::  svp1  = 0.61078
   real,parameter    ::  svp2  = 17.2693882
   real,parameter    ::  svpt0 = 273.15
#endif
   real,parameter    ::  svpi2 = 21.8745584
   real,parameter    ::  svp3  = 35.86
   real,parameter    ::  svpi3 = 7.66
!
   integer :: i,j,k
!-------------------------------------------------------------------------------
!
   do j = jts,jte
!
      do k = kts, kte
         do i = its, ite
            t_p(i,k)     = t_p_3d(i,k,j)
            pres_p(i,k)  = pres_p_3d(i,k,j)
            qv(i,k)      = qv_3d(i,k,j)
            qc(i,k)      = qc_3d(i,k,j)
            qi(i,k)      = qi_3d(i,k,j)
            qs(i,k)      = qs_3d(i,k,j)
         end do
      end do
!
      call cf_xurandall_run(cldfrac,t_p,pres_p,qv,qc,qi,qs,   &
                            svp1, svp2, svpi2, svp3,          &
                            svpi3, svpt0, ep_2,               &
                            its,ite,kte,                      &
                            errmsg,errflg)
!
      do k = kts, kte
         do i = its, ite
            cldfrac_3d(i,k,j) = cldfrac(i,k)
         end do
      end do
!
   end do

   end subroutine xurandall
!=================================================================================================================
 end module module_cf_xurandall
!=================================================================================================================
