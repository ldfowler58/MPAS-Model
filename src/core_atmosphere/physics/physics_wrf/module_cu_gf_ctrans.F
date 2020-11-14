!=================================================================================================================
 module module_cu_gf_ctrans
 use mpas_log
 use mpas_timer,only: mpas_timer_start,mpas_timer_stop
 use mpas_kind_types,only: kind_phys=>RKIND

 use mpas_atmphys_constants,only: grav=>gravity,xlv
!use module_cu_gf_fct,only      : fct1d3
 use module_cu_gf_shared,only   : get_cloud_bc

 implicit none
 private
 public:: cup_transport_clev,    &
          cup_transport_dn,      &
          cup_transport_tends,   &
          cup_up_moisture_ctrans

!subroutine cup_transport_clev computes the number concentration of cloud liquid water, cloud ice, "ice-friendly"
!and "water-friendly" aerosols at cloud levels, i.e. at the interface between model layers. cup_transport_clev
!also computes the mixing ratio (number concentration) of passive tracers and chemical species at cloud levels.

!subroutine cup_up_moisture_ctrans is a cleaned-up version of subroutine cup_up_moisture.

!subroutine cup_transport_dn computes the number concentration of cloud liquid water, cloud ice, "ice-friendly"
!and "water-friendly" aerosols at cloud levels in downdrafts. also computes mixing ratio and number concentration
!of passive tracers and chemical species in downdrafts.

!subroutine cup_transport_tends computes the tendencies for the number concentration of cloud liquid water, cloud
!ice, "ice-friendly" and "water-friendly" aerosols. also computes the tendencies of passive tracers and chemical
!species.

!subroutine cup_transport_up_moisture is the cleaned-up subroutine cup_up_moisture. Results are exactly the same
!as in module_cu_gf_deep.F when using the "-g" or "-O0 -fp-model precise" options. It does not work with autoconv
!is equal to 2. It includes the convective transport for the number concentration of cloud liquid water, cloud
!ice, "ice-friendly" and "water-friendly" aerosols from the Thompson cloud microphysics scheme. It also includes
!the convective transport of passive tracers. Note that the convective transport of passive tracers could be in
!its own subroutine since they do not interact with moisture, it was more convenient to put them here.

!add-ons and modifications to sourcecode:
!----------------------------------------
! * corrected how to limit the cloud base mass flux.
! * added the convective transport of "ice-friendly" and "water-friendly" Thompson aerosols.
!   Laura D. Fowler (laura@ucar.edu) / 2020-11-14.


 contains


!=================================================================================================================
 subroutine cup_transport_clev(its,ite,itf,kts,kte,ktf,nc,ni,nifa,nwfa,ncl2,nil2,nifal2,nwfal2, &
                               nchem,chem,cheml2,ntracer,tracer,tracerl2)
!=================================================================================================================

!input arguments:
 integer,intent(in):: its,ite,itf,kts,kte,ktf
 integer,intent(in):: nchem,ntracer

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte):: nc,ni,nifa,nwfa
 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte,nchem)  :: chem
 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte,ntracer):: tracer

!inout arguments:
 real(kind=kind_phys),intent(inout),dimension(its:ite,kts:kte):: ncl2,nil2,nifal2,nwfal2
 real(kind=kind_phys),intent(inout),dimension(its:ite,kts:kte,nchem)  :: cheml2
 real(kind=kind_phys),intent(inout),dimension(its:ite,kts:kte,ntracer):: tracerl2

!local variables and arrays:
 integer:: i,ic,it,k

 real(kind=kind_phys),parameter:: half = 0.5_kind_phys

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
!call mpas_log_write('--- enter subroutine cup_transport_clev:')

 k = kts
 do i = its, itf
    ncl2(i,k) = nc(i,k)
    nil2(i,k) = ni(i,k)
    nifal2(i,k) = nifa(i,k)
    nwfal2(i,k) = nwfa(i,k)
 enddo
 do k = kts+1, ktf
    do i = its, itf
       ncl2(i,k) = half*(nc(i,k-1) + nc(i,k))
       nil2(i,k) = half*(ni(i,k-1) + ni(i,k))
       nifal2(i,k) = half*(nifa(i,k-1) + nifa(i,k))
       nwfal2(i,k) = half*(nwfa(i,k-1) + nwfa(i,k))
    enddo
 enddo

 if(nchem > 0) then
    do ic = 1, nchem
       k = kts
       do i = its, itf
          cheml2(i,k,ic) = chem(i,k,ic)
       enddo
       do k = kts+1, ktf
          do i = its, itf
             cheml2(i,k,ic) = half*(chem(i,k-1,ic) + chem(i,k,ic))
          enddo
       enddo
    enddo
 endif

 if(ntracer > 0) then
    do it = 1, ntracer
       k = kts
       do i = its, itf
          tracerl2(i,k,it) = tracer(i,k,it)
       enddo
       do k = kts+1, ktf
          do i = its, itf
!            tracerl2(i,k,it) = half*(tracer(i,k-1,it) + tracer(i,k,it))
             tracerl2(i,k,it) = tracer(i,k-1,it)
          enddo
       enddo
    enddo
 endif

!call mpas_log_write('--- enter subroutine cup_transport_clev:')

 end subroutine cup_transport_clev

!=================================================================================================================
 subroutine cup_up_moisture_ctrans(name,ierr,k22,kbcon,ktop,autoconv,z_cup,t_cup,qe_cup,qes_cup, &
                     gamma_cup,q,dby,c1d,up_massentr,up_massdetr,zu,qc,qrc,pw,clw_all,pwav,psum, &
#ifdef ATM_CHEMISTRY
                     is_mp,nifa,nwfa,nifal2,nwfal2,nifa_up,nwfa_up,                              &
                     nchem,chem,cheml2,chem_up,ntracer,tracer,tracerl2,tracer_up,                &
#endif
                     its,ite,itf,kts,kte,ktf)
!=================================================================================================================

!input arguments:
 character *(*),intent(in):: &
    name            ! type of convection (shallow,mid-convection,deep).

 integer,intent(in):: &
    its,ite,itf,   &! indices for horizontal dimensions.
    kts,kte,ktf     ! indices for vertical dimensions.

 integer,intent(in):: &
    autoconv        !

 integer,intent(in),dimension(its:ite):: &
    k22,           &!
    kbcon,         &! index of convective cloud base.
    ktop            ! index of convective cloud top.

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte):: &
    up_massentr,   &! normalized entrainment mass flux.
    up_massdetr,   &! normalized detrainment mass flux.
    zu              ! normalized updraft mass flux.

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte):: &
    t_cup,         &!
    qe_cup,        &!
    qes_cup,       &!
    z_cup           !

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte):: &
    q               ! layer water vapor mixing ratio.

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte):: &
    gamma_cup,     &!
    dby,           &!
    c1d             !

!inout arguments:
 integer,intent(inout),dimension(its:ite) :: &
    ierr            !

!output arguments:
 real(kind=kind_phys),intent(out),dimension(its:ite):: &
    pwav,&          !
    psum            !

 real(kind=kind_phys),intent(out),dimension(its:ite,kts:kte):: &
    clw_all,&       !
    pw,&            !
    qc,&            !
    qrc

!local variables and arrays:
 integer,parameter:: iall = 0
 integer:: i,k,k2,kb,kt

!real(kind=kind_phys),parameter:: c0 = 0.004
 real(kind=kind_phys):: top,bottom
 real(kind=kind_phys):: c0,denom,dz,qaver,qrch

 real(kind=kind_phys),dimension(kts:kte):: dummy1
 real(kind=kind_phys),dimension(its:ite,kts:kte):: c0t,clw_allh,pwh,qch,qrcb

#ifdef ATM_CHEMISTRY
!input arguments:
 logical,intent(in):: &
    is_mp           ! if is_mp is true, then compute convective transport of nifa and nwfa.

 integer,intent(in):: &
    nchem,         &! number of chemical species.
    ntracer         ! number of passive tracers.

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte):: &
    nifa,          &! layer environmental "ice-friendly" aerosol number concentration from Thompson microphysics.
    nwfa            ! layer environmental "water-friendly" aerosol number concentration from Thompson microphysics.
 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte):: &
    nifal2,        &! level environmental "ice-friendly" aerosol number concentration from Thompson microphysics.
    nwfal2          ! level environmental "water-friendly" aerosol number concentration from Thompson microphysics.
 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte,nchem):: &
    chem,          &!
    cheml2          !
 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte,ntracer):: &
    tracer,        &!
    tracerl2        !

!output arguments:
 real(kind=kind_phys),intent(inout),dimension(its:ite,kts:kte):: &
    nifa_up,&       !
    nwfa_up
 real(kind=kind_phys),intent(inout),dimension(its:ite,kts:kte,nchem):: &
    chem_up         !
 real(kind=kind_phys),intent(inout),dimension(its:ite,kts:kte,ntracer):: &
    tracer_up       !

!local variables:
 integer:: ic,it
 real(kind=kind_phys):: nifa_aver,nwfa_aver,chem_aver,tracer_aver
#endif

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine cup_up_moisture_ctrans:')

!--- initialization:
 do i = its, itf
    pwav(i)  = 0._kind_phys
    psum(i)  = 0._kind_phys
 enddo
 do k = kts, ktf
    do i = its, itf
       pw(i,k)  = 0._kind_phys
       pwh(i,k) = 0._kind_phys

       clw_all(i,k)  = 0._kind_phys
       clw_allh(i,k) = 0._kind_phys

       qrc(i,k)  = 0._kind_phys
       qrcb(i,k) = 0._kind_phys
       qc(i,k)   = 0._kind_phys
       qch(i,k)  = 0._kind_phys
       if(ierr(i) .eq. 0) qc(i,k)  = qe_cup(i,k)
       if(ierr(i) .eq. 0) qch(i,k) = qe_cup(i,k)

       !--- conversion to rain:
       !if(t_cup(i,k) > 273.16) then
       !   c0t(i,k) = c0
       !else
       !   c0t(i,k) = c0*exp(0.07*(t_cup(i,k)-273.16))
       !endif
    enddo
 enddo
 do i = its, itf
    if(ierr(i) .eq. 0) then
       k2 = k22(i)
       do k = kts, kte
          dummy1(k) = qe_cup(i,k)
       enddo
       call get_cloud_bc(kte,dummy1,qaver,k2)
       qc(i,k2)  = qaver
       qch(i,k2) = qaver
       do k = kts, k2-1
          qc(i,k)  = qe_cup(i,k)
          qch(i,k) = qe_cup(i,k)
       enddo
    endif
 enddo

#ifdef ATM_CHEMISTRY
 call mpas_timer_start('atm_chem: physics - GF')
!--- number concentrations:
 if(is_mp) then
    do i = its, itf
       if(ierr(i) .eq. 0) then
          k = kts
          nifa_up(i,k) = nifa(i,k)
          nwfa_up(i,k) = nwfa(i,k)
          do k = kts+1, ktf
             nifa_up(i,k) = nifa(i,k-1)
             nwfa_up(i,k) = nwfa(i,k-1)
          enddo
       endif
    enddo
 endif
!--- chemical species:
 do ic = 1, nchem
    do i = its, itf
       if(ierr(i) .eq. 0) then
          k = kts
          chem_up(i,k,ic) = chem(i,k,ic)
          do k = kts+1, ktf
             chem_up(i,k,ic) = chem(i,k-1,ic)
          enddo
       endif
    enddo
 enddo
!--- passive tracers:
 do it = 1, ntracer
    do i = its, itf
       if(ierr(i) .eq. 0) then
          k = kts
          tracer_up(i,k,it) = tracer(i,k,it)
          do k = kts+1, ktf
             tracer_up(i,k,it) = tracer(i,k-1,it)
          enddo
       endif
    enddo
 enddo
 call mpas_timer_stop('atm_chem: physics - GF')
#endif

!---
 do i = its, itf
    if(ierr(i) .eq. 0) then
       k2 = k22(i)
       kb = kbcon(i)
       kt = ktop(i)

       !--- below cloud base:
       do k = k2+1, kb
          dz = z_cup(i,k)-z_cup(i,k-1)
          c0 = 0.004
          if(t_cup(i,k) > 273.16) then
             c0t(i,k) = c0
          else
             c0t(i,k) = c0*exp(0.07*(t_cup(i,k)-273.16))
          endif

          qc(i,k)  = (qc(i,k-1)*zu(i,k-1) - .5*up_massdetr(i,k-1)*qc(i,k-1)   &
                     + up_massentr(i,k-1)*q(i,k-1))                           &
                   / (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
          qrch = qes_cup(i,k) + 1./xlv*(gamma_cup(i,k)/(1.+gamma_cup(i,k)))*dby(i,k)

          if(k .lt. kb) qrch = qc(i,k)
          if(qc(i,k) .gt. qrch) then
             qrc(i,k)     = (qc(i,k)-qrch)/(1.+c0t(i,k)*dz)
             pw(i,k)      = c0t(i,k)*dz*qrc(i,k)*zu(i,k)
             qc(i,k)      = qrch+qrc(i,k)
             clw_all(i,k) = qrc(i,k)
          endif

#ifdef ATM_CHEMISTRY
          call mpas_timer_start('atm_chem: physics - GF')
          !--- number concentrations:
          if(is_mp) then
             nifa_up(i,k) = (nifa_up(i,k-1)*zu(i,k-1) - .5*up_massdetr(i,k-1)*nifa_up(i,k-1) &
                          + up_massentr(i,k-1)*nifa(i,k-1))                                 &
                          / (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
             nwfa_up(i,k) = (nwfa_up(i,k-1)*zu(i,k-1) - .5*up_massdetr(i,k-1)*nwfa_up(i,k-1) &
                          + up_massentr(i,k-1)*nwfa(i,k-1))                                 &
                          / (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
          endif
!         !--- chemical species:
!         if(nchem > 0) then
!            do ic = 1,nchem
!               chem_up(i,k,ic) = (chem_up(i,k-1,ic)*zu(i,k-1) - .5*up_massdetr(i,k-1)*chem_up(i,k-1,ic) &
!                               + up_massentr(i,k-1)*chem(i,k-1,ic))                                     &
!                               / (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
!            enddo
!         endif
          !--- passive tracers:
          if(ntracer > 0) then
             do it = 1,ntracer
                tracer_up(i,k,it) = (tracer_up(i,k-1,it)*zu(i,k-1) - .5*up_massdetr(i,k-1)*tracer_up(i,k-1,it) &
                                  + up_massentr(i,k-1)*tracer(i,k-1,it))                                       &
                                  / (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
             enddo
          endif
          call mpas_timer_stop('atm_chem: physics - GF')
#endif
       enddo

       !--- above cloud base:
       do k = kb+1, kt

          !--- I think that this check on the value of denom can be moved in the initialization:
          denom = zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1)
          if(denom .lt. 1.e-16) then
             ierr(i) = 51
             exit
          endif

          c0 = 0.004
          if(t_cup(i,k) .lt. 270.) c0 = 0.002
          if(t_cup(i,k) > 273.16) then
             c0t(i,k) = c0
          else
             c0t(i,k) = c0*exp(0.07*(t_cup(i,k)-273.16))
          endif

          qc(i,k)  = (qc(i,k-1)*zu(i,k-1) - .5*up_massdetr(i,k-1)*qc(i,k-1)   &
                   + up_massentr(i,k-1)*q(i,k-1))                             &
                   / (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))

          qch(i,k) = (qch(i,k-1)*zu(i,k-1) - .5*up_massdetr(i,k-1)*qch(i,k-1) &
                   + up_massentr(i,k-1)*q(i,k-1))                             &
                   / (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))

          qrch = qes_cup(i,k)+(1./xlv)*(gamma_cup(i,k)/(1.+gamma_cup(i,k)))*dby(i,k)
          if(qc(i,k) .le. qrch)  qc(i,k)  = qrch
          if(qch(i,k) .le. qrch) qch(i,k) = qrch

          !--- total condensed water before rainout:
          clw_all(i,k)  = max(0.,qc(i,k)-qrch)
          clw_allh(i,k) = max(0.,qch(i,k)-qrch)
          qrc(i,k)      = max(0.,(qc(i,k)-qrch))  ! /(1.+c0*dz*zu(i,k))
          qrcb(i,k)     = max(0.,(qch(i,k)-qrch)) ! /(1.+c0*dz*zu(i,k))

          dz = z_cup(i,k)-z_cup(i,k-1)
          if(autoconv .eq. 2) then

          else
             if(iall .eq. 0) then
                qrc(i,k) = (qc(i,k)-qrch)/(1.+(c1d(i,k)+c0t(i,k))*dz)
                pw(i,k)  = c0t(i,k)*dz*qrc(i,k)*zu(i,k)
                if(qrc(i,k) .lt. 0.) then
                   qrc(i,k) = 0.
                   pw(i,k)  = 0.
                endif
             else
                qrc(i,k) = 0.
                pw(i,k)  = (qc(i,k)-qrch)*zu(i,k)
                if(pw(i,k).lt.0.)pw(i,k)=0.
             endif
             qc(i,k) = qrc(i,k) + qrch
          endif
          pwav(i) = pwav(i) + pw(i,k)
          psum(i) = psum(i) + clw_all(i,k)*zu(i,k) *dz

#ifdef ATM_CHEMISTRY
          !--- number concentrations:
          if(is_mp) then
             nifa_up(i,k) = (nifa_up(i,k-1)*zu(i,k-1) - .5*up_massdetr(i,k-1)*nifa_up(i,k-1) &
                          + up_massentr(i,k-1)*nifa(i,k-1))                                  &
                          / (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
             nwfa_up(i,k) = (nwfa_up(i,k-1)*zu(i,k-1) - .5*up_massdetr(i,k-1)*nwfa_up(i,k-1) &
                          + up_massentr(i,k-1)*nwfa(i,k-1))                                  &
                          / (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
          endif
!         !--- chemical species:
!         if(nchem > 0) then
!            do ic = 1, nchem
!               chem_up(i,k,ic) = (chem_up(i,k-1,ic)*zu(i,k-1) - .5*up_massdetr(i,k-1)*chem_up(i,k-1,ic) &
!                               + up_massentr(i,k-1)*chem(i,k-1,ic))                                     &
!                               / (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
!            enddo
!         endif

          !--- passive tracers:
          if(ntracer > 0) then
             do it = 1, ntracer
                tracer_up(i,k,it) = (tracer_up(i,k-1,it)*zu(i,k-1) - .5*up_massdetr(i,k-1)*tracer_up(i,k-1,it) &
                                  + up_massentr(i,k-1)*tracer(i,k-1,it))                                       &
                                  / (zu(i,k-1)-.5*up_massdetr(i,k-1)+up_massentr(i,k-1))
             enddo
          endif
#endif
       enddo

       !--- do not include liquid/ice in qc
       do  k = k22(i)+1,ktop(i)
           qc(i,k) = qc(i,k)-qrc(i,k)
       enddo

    endif
 enddo

 call mpas_log_write('--- end subroutine cup_up_moisture_ctrans:')

 end subroutine cup_up_moisture_ctrans

!=================================================================================================================
 subroutine cup_transport_dn(its,ite,itf,kts,kte,ktf,ierr,kdn,zd,dn_massentr,dn_massdetr,is_mp,nifa,nwfa,     &
                             nifal2,nwfal2,nifa_dn,nwfa_dn,nchem,chem,cheml2,chem_dn,ntracer,tracer,tracerl2, &
                             tracer_dn)
!=================================================================================================================

!input arguments:
 logical,intent(in):: &
    is_mp           !

 integer,intent(in):: &
    its,ite,itf,   &! indices for horizontal dimensions.
    kts,kte,ktf     ! indices for vertical dimensions.

 integer,intent(in):: &
    nchem,         &! number of chemical species.
    ntracer         ! number of passive tracers.

 integer,intent(in),dimension(its:ite):: &
    kdn             ! index of downdraft originating level.

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte):: &
    dn_massentr,   &! normalized entrainment mass flux.
    dn_massdetr,   &! normalized detrainment mass flux.
    zd              ! normalized downdraft mass flux.

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte):: &
    nifa,          &! layer environmental "ice-friendly" aerosol number concentration from Thompson microphysics.
    nwfa            ! layer environmental "water-friendly" aerosol number concentration from Thompson microphysics.

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte):: &
    nifal2,&        ! level environmental "ice-friendly" aerosol number concentration from Thompson microphysics.
    nwfal2          ! level environmental "water-friendly" aerosol number concentration from Thompson microphysics.

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte,nchem):: &
    chem,          &!
    cheml2          !

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte,ntracer):: &
    tracer,        &!
    tracerl2        !


!inout arguments:
 integer,intent(inout),dimension(its:ite) :: &
    ierr            !

 real(kind=kind_phys),intent(inout),dimension(its:ite,kts:kte):: &
    nifa_dn,       &!
    nwfa_dn         !

 real(kind=kind_phys),intent(inout),dimension(its:ite,kts:kte,nchem):: &
    chem_dn         !

 real(kind=kind_phys),intent(inout),dimension(its:ite,kts:kte,ntracer):: &
    tracer_dn       !


!local variables and arrays:
 integer:: i,ic,it,k,kd

 real(kind=kind_phys):: denom

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine cup_transport_dn:')

!--- initialization:
!--- number concentrations:
 if(is_mp) then
    do i = its,itf
       if(ierr(i) .eq. 0) then
          do k = kts,ktf
             nifa_dn(i,k) = nifa(i,k)
             nwfa_dn(i,k) = nwfa(i,k)
          enddo
       endif
    enddo
 endif
!--- passive tracers:
 do it = 1,ntracer
    do i = its,itf
       if(ierr(i) .eq. 0) then
          do k = kts,ktf
             tracer_dn(i,k,it) = tracer(i,k,it)
          enddo
       endif
    enddo
 enddo

!--- computation of the in-cloud cloud liquid water and ice number concentrations,"ice-friendly" and
!    "water-friendly" number concentrations from the Thompson cloud microphysics scheme,and in-cloud
!    mass mixing ratios of passive tracers and chemical species due to cloud downdrafts:
 do i = its,itf
    if(ierr(i) .eq. 0) then
       kd = kdn(i)
       do k = kd-1,kts,-1
          denom = zd(i,k+1) -.5*dn_massdetr(i,k) + dn_massentr(i,k)
          if(denom .lt. 1.e-16) ierr(i) = 51
       enddo
    endif

    if(ierr(i) .eq. 0) then
       kd = kdn(i)
       if(is_mp) then
          do k = kd-1,kts,-1
             nifa_dn(i,k) = (nifa_dn(i,k+1)*zd(i,k+1) -.5*dn_massdetr(i,k)*nifa_dn(i,k+1) &
                          +  dn_massentr(i,k)*nifa(i,k))                                  &
                          / (zd(i,k+1)-.5*dn_massdetr(i,k)+dn_massentr(i,k))
             nwfa_dn(i,k) = (nwfa_dn(i,k+1)*zd(i,k+1) -.5*dn_massdetr(i,k)*nwfa_dn(i,k+1) &
                          +  dn_massentr(i,k)*nwfa(i,k))                                  &
                          / (zd(i,k+1)-.5*dn_massdetr(i,k)+dn_massentr(i,k))
          enddo
       endif
!      if(nchem > 0) then
!         do ic = 1,nchem
!            chem_dn(i,kd,ic) = cheml2(i,kd,ic)
!            do k = kd-1,kts,-1
!               chem_dn(i,k,ic) = (chem_dn(i,k+1,ic)*zd(i,k+1) -.5*dn_massdetr(i,k)*chem_dn(i,k+1,ic) &
!                                +  dn_massentr(i,k)*chem(i,k,ic))                                    &
!                                / (zd(i,k+1)-.5*dn_massdetr(i,k)+dn_massentr(i,k))
!            enddo
!         enddo
!      endif
       if(ntracer > 0) then
          do it = 1,ntracer
             do k = kd-1,kts,-1
                tracer_dn(i,k,it) = (tracer_dn(i,k+1,it)*zd(i,k+1) -.5*dn_massdetr(i,k)*tracer_dn(i,k+1,it) &
                                  +  dn_massentr(i,k)*tracer(i,k,it))                                       &
                                  / (zd(i,k+1)-.5*dn_massdetr(i,k)+dn_massentr(i,k))
             enddo
          enddo
       endif
    endif
 enddo

 call mpas_log_write('--- exit subroutine cup_transport_dn:')

 end subroutine cup_transport_dn

!=================================================================================================================
 subroutine cup_transport_tends(its,ite,itf,kts,kte,ktf,ierr,k22,kbcon,ktop,kdn,kpbl,dtime,xmb,edt,zu,zd,pres, &
                               is_mp,nifa,nifal2,nwfa,nwfal2,nifa_up,nwfa_up,nifa_dn,nwfa_dn,outnifa,outnwfa,  &
                               nchem,chem,cheml2,chem_up,chem_dn,outchem,ntracer,tracer,tracerl2,tracer_up,    &
                               tracer_dn,outtracer)
!=================================================================================================================

!input arguments:
 logical,intent(in):: &
    is_mp           !

 integer,intent(in):: &
    its,ite,itf,   &! indices for horizontal dimensions.
    kts,kte,ktf     ! indices for vertical dimensions.

 integer,intent(in):: &
    nchem,         &!
    ntracer         !

 integer,intent(in),dimension(its:ite) :: &
    ierr,          &!
    k22,           &!
    kbcon,         &!
    ktop,          &!
    kdn,           &!
    kpbl

 real(kind=kind_phys),intent(in):: &
    dtime           !

 real(kind=kind_phys),intent(in),dimension(its:ite):: &
    edt,           &!
    xmb             !

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte):: &
    pres,          &!
    zd,            &!
    zu

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte):: &
    nifa,          &!
    nwfa,          &!
    nifal2,        &! level environmental "ice-friendly" aerosol number concentration from Thompson microphysics.
    nwfal2          ! level environmental "water-friendly" aerosol number concentration from Thompson microphysics.

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte,nchem):: &
    chem,          &!
    cheml2          ! level environmental chemical species.

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte,ntracer):: &
    tracer,        &!
    tracerl2        ! level environmental passive tracers.

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte):: &
    nifa_dn,       &!
    nwfa_dn,       &!
    nifa_up,       &!
    nwfa_up         !

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte,nchem):: &
    chem_up,       &!
    chem_dn         !

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte,ntracer):: &
    tracer_up,     &!
    tracer_dn       !

!inout arguments:
 real(kind=kind_phys),intent(inout),dimension(its:ite,kts:kte):: &
    outnifa,       &!
    outnwfa         !

 real(kind=kind_phys),intent(inout),dimension(its:ite,kts:kte,nchem):: &
    outchem         !

 real(kind=kind_phys),intent(inout),dimension(its:ite,kts:kte,ntracer):: &
    outtracer       !

!local variables and arrays:
 integer:: i,ic,it,k,k2,kb,kd,kt

 real(kind=kind_phys):: sumi,sumf

 real(kind=kind_phys),parameter:: epsilon = 1.e-10
 real(kind=kind_phys),dimension(its:ite,kts:kte):: ratio

 real(kind=kind_phys),dimension(its:ite):: cflmax,xmbl
 real(kind=kind_phys),dimension(kts:kte):: cfl
 real(kind=kind_phys),dimension(its:ite,kts:kte):: aterm,bterm,esubs
 real(kind=kind_phys),dimension(its:ite,kts:kte,ntracer):: tracerp,tracerpu_expl,tracerpd_expl

 real(kind=kind_phys):: numerator,denominator
 real(kind=kind_phys),dimension(kts:kte):: aterm2,bterm2,cterm2,dterm2,esubs2,x2
 real(kind=kind_phys),dimension(kts:kte):: aterm3,bterm3,cterm3,dterm3,esubs3,x3
 real(kind=kind_phys),dimension(kts:kte):: num,den,fm,fp,znet
 real(kind=kind_phys),dimension(its:ite,kts:kte,ntracer):: outt

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- subroutine cup_transport_tends: itf = $i  ktf = $i',intArgs=(/itf,ktf/))

!--- CFL (maybe I only need to limit the CFL conditions between the cloud-base and cloud-top levels):
 do i = its, ite
    xmbl(i) = xmb(i)
    do k = kts, kte
       ratio(i,k) = 0._kind_phys
    enddo
 enddo

 do i = its, itf
    if(ierr(i) .eq. 0) then
       k2 = k22(i)
       kb = kbcon(i)
       kt = ktop(i)
       kd = kdn(i)
       cflmax(i) = 0._kind_phys

!      do k = kts, ktf-1
!         ratio(i,k) = 0.01_kind_phys*grav/(pres(i,k) - pres(i,k+1))
!         cfl(k) = zu(i,k)/(ratio(i,k)*dtime)
!         cfl(k) = xmbl(i)*cfl(k)
!         cflmax(i) = max(cflmax(i),cfl(k))
!      enddo
!      if(xmbl(i).gt.cflmax(i)) xmbl(i) = cflmax(i)
       do k = kb,kt
          ratio(i,k) = 0.01_kind_phys*grav/(pres(i,k) - pres(i,k+1))
          cflmax(k)  = 1._kind_phys/(ratio(i,k)*dtime)
          znet(k)    = zu(i,k)-edt(i)*zd(i,k)
          cfl(k)     = xmb(i)*znet(k)
          if(cfl(k).gt.cflmax(k)) xmbl(i) = min(xmbl(i),cflmax(k)/znet(k))
       enddo
       do k = kts, ktf
          ratio(i,k) = ratio(i,k)*xmbl(i)
       enddo

!      do k = kts, ktf
!         call mpas_log_write('$i $i $i $i $i $i $r $r $r $r $r $r',intArgs=(/i,kd,k2,kb,kt,k/),realArgs= &
!              (/nifa_up(i,k),nifa_dn(i,k),outnifa(i,k),nwfa_up(i,k),nwfa_dn(i,k),outnwfa(i,k)/))
!      enddo
!      do k = kts, ktf
!         call mpas_log_write('$i $i $i $i $i $i $r $r $r $r $r $r $r $r $r',intArgs=(/i,kd,k2,kb,kt,k/),realArgs= &
!              (/(chem_up(i,k,ic),chem_dn(i,k,ic),outchem(i,k,ic),ic=1,nchem)/))
!      enddo
!      do k = kts, ktf
!         call mpas_log_write('$i $i $i $i $i $i $r $r $r $r $r $r $r $r',intArgs=(/i,kd,k2,kb,kt,k/),realArgs= &
!              (/(tracer(i,k,it),tracer_up(i,k,it),tracer_dn(i,k,it),outtracer(i,k,it),it=1,ntracer)/))
!      enddo
    endif
 enddo

!--- convective transport of "ice-friendly" and "water-friendly" aerosols from the Thompson cloud microphysics:
 if(is_mp) then
    do i = its, itf
       if(ierr(i) .eq. 0) then
          k2 = k22(i)
          kb = kbcon(i)
          kt = ktop(i)
          kd = kdn(i)

          do k = kts, ktf
             znet(k) = zu(i,k)-edt(i)*zd(i,k)
!            znet(k) = zu(i,k)
             fp(k) = 0.5*(znet(k)+abs(znet(k)))
             fm(k) = 0.5*(znet(k)-abs(znet(k)))
          enddo
          do k = kts, ktf-1
             aterm3(k) = ratio(i,k)*dtime*fm(k)
             bterm3(k) = 1.+ratio(i,k)*dtime*(fp(k)-fm(k+1))
             cterm3(k) = -ratio(i,k)*dtime*fp(k+1)
         enddo
         !nifa:
         do k = kts, ktf-1
             dterm3(k) = (zu(i,k)*nifa_up(i,k)-zu(i,k+1)*nifa_up(i,k+1)) &
                       - edt(i)*(zd(i,k)*nifa_dn(i,k)-zd(i,k+1)*nifa_dn(i,k+1))
             dterm3(k) = nifa(i,k) + ratio(i,k)*dtime*dterm3(k)
          enddo
          k = ktf
          dterm3(k) = nifa(i,k)
          call tridiag(ktf-1,aterm3(kts:ktf-1),bterm3(kts:ktf-1),cterm3(kts:ktf-1),dterm3(kts:ktf-1))
          do k = kts, ktf
             outnifa(i,k) = (dterm3(k)-nifa(i,k))/dtime
          enddo
          !nwfa:
          do k = kts, ktf-1
             dterm3(k) = (zu(i,k)*nwfa_up(i,k)-zu(i,k+1)*nwfa_up(i,k+1)) &
                       - edt(i)*(zd(i,k)*nwfa_dn(i,k)-zd(i,k+1)*nwfa_dn(i,k+1))
             dterm3(k) = nwfa(i,k) + ratio(i,k)*dtime*dterm3(k)
          enddo
          k = ktf
          dterm3(k) = nwfa(i,k)
          call tridiag(ktf-1,aterm3(kts:ktf-1),bterm3(kts:ktf-1),cterm3(kts:ktf-1),dterm3(kts:ktf-1))
          do k = kts, ktf
             outnwfa(i,k) = (dterm3(k)-nwfa(i,k))/dtime
          enddo
!         do k = kts, ktf
!            call mpas_log_write('$i $i $i $i $i $i $r $r $r $r',intArgs=(/i,kd,k2,kb,kt,k/), &
!                      realArgs=(/nifa(i,k),outnifa(i,k),nwfa(i,k),outnwfa(i,k)/))
!         enddo
!         call mpas_log_write(' ')
       endif
    enddo
 endif

 do i = its, itf
    if(ierr(i) .eq. 0) then
       k2 = k22(i)
       kb = kbcon(i)
       kt = ktop(i)
       kd = kdn(i)

       do k = kts, ktf
          znet(k) = zu(i,k)-edt(i)*zd(i,k)
!         znet(k) = zu(i,k)
          fp(k) = 0.5*(znet(k)+abs(znet(k)))
          fm(k) = 0.5*(znet(k)-abs(znet(k)))

          den(k)    = 0._kind_phys
          num(k)    = 0._kind_phys

          aterm2(k) = 0._kind_phys
          bterm2(k) = 0._kind_phys
          cterm2(k) = 0._kind_phys
          dterm2(k) = 0._kind_phys

          x3(k) = 0._kind_phys
          aterm3(k) = 0._kind_phys
          bterm3(k) = 0._kind_phys
          cterm3(k) = 0._kind_phys
          dterm3(k) = 0._kind_phys
       enddo

       do it = 1, ntracer
          do k = kts, ktf
             x2(k) = tracer(i,k,it)
          enddo

!--- implicit formulation for updrafts only (downward sweep):
          do k = kts, ktf-1
             aterm(i,k) = ratio(i,k)*dtime*zu(i,k)
             esubs(i,k) = zu(i,k)*tracer_up(i,k,it)-zu(i,k+1)*tracer_up(i,k+1,it)
             esubs(i,k) = ratio(i,k)*dtime*esubs(i,k)
          enddo

          !updated tracers in layers above layer index k = kt:
          do k = ktf, kt+1, -1
             tracerp(i,k,it) = tracer(i,k,it)
          enddo
          !updated tracer in layer above the cloud-top:
          k = kt
          numerator   = zu(i,k)*tracer_up(i,k,it) - zu(i,k+1)*tracer_up(i,k+1,it)
          numerator   = ratio(i,k)*dtime*numerator + tracer(i,k,it)
          denominator = 1. + ratio(i,k)*dtime*zu(i,k)
          tracerp(i,k,it) = numerator/denominator
          num(k) = numerator
          den(k) = denominator
          !updated tracers in layers below the cloud-top:
          do k = kt-1, kb-1, -1
             numerator   = zu(i,k)*tracer_up(i,k,it) - zu(i,k+1)*tracer_up(i,k+1,it)
             numerator   = numerator + zu(i,k+1)*tracerp(i,k+1,it)
             numerator   = ratio(i,k)*dtime*numerator + tracer(i,k,it)
             denominator = 1. + ratio(i,k)*dtime*zu(i,k)
             tracerp(i,k,it) = numerator/denominator
             num(k) = numerator
             den(k) = denominator
          enddo
          do k = kb-2, kts, -1
             tracerp(i,k,it) = tracer(i,k,it)
          enddo
          !tendencies:
          do k = kts, ktf
             outt(i,k,it) = (tracerp(i,k,it)-tracer(i,k,it))/dtime
          enddo

!         do k = kts, ktf-1
!            aterm2(k) = 0.
!            bterm2(k) = 1.+ratio(i,k)*dtime*zu(i,k)
!            cterm2(k) = -ratio(i,k)*dtime*zu(i,k+1)
!            dterm2(k) = tracer(i,k,it) + esubs(i,k)
!         enddo
!         call tridiag2(ktf-1,aterm2,bterm2,cterm2,dterm2,x2)

!--- implicit formulation of updrafts plus downdrafts:
          do k = kts, ktf-1
             aterm3(k) = ratio(i,k)*dtime*fm(k)
             bterm3(k) = 1.+ratio(i,k)*dtime*(fp(k)-fm(k+1))
             cterm3(k) = -ratio(i,k)*dtime*fp(k+1)
             dterm3(k) = (zu(i,k)*tracer_up(i,k,it)-zu(i,k+1)*tracer_up(i,k+1,it)) &
                       - edt(i)*(zd(i,k)*tracer_dn(i,k,it)-zd(i,k+1)*tracer_dn(i,k+1,it))
             dterm3(k) = tracer(i,k,it) + ratio(i,k)*dtime*dterm3(k)
          enddo
          k = ktf
          dterm3(k) = tracer(i,k,it)
          call tridiag(ktf-1,aterm3(kts:ktf-1),bterm3(kts:ktf-1),cterm3(kts:ktf-1),dterm3(kts:ktf-1))
          do k = kts, ktf
             outtracer(i,k,it) = (dterm3(k)-tracer(i,k,it))/dtime
          enddo
!         if(it .eq. 1) then
!            do k = kts, ktf
!               call mpas_log_write('$i $i $r $r $r $r $r',intArgs=(/i,k/),realArgs=(/tracer_up(i,k,it), &
!                    tracer_dn(i,k,it),tracer(i,k,it),dterm3(k),(dterm3(k)-tracer(i,k,it))/dtime/))
!            enddo
!            call mpas_log_write(' ')
!         endif
       enddo

    endif
 enddo

 call mpas_log_write('--- end subroutine cup_transport_tends:')
!call mpas_log_write(' ')

 end subroutine cup_transport_tends

!=================================================================================================================
 subroutine tridiag(m,a,b,c,f)
!=================================================================================================================
!copied as is from module ConvPar_GF_GEOS5.F90, as written by Saulo Freita.
!this routine solves the problem: aa*f(k-1,t+1) + bb*f(k,t+1) + cc*f(k+1,t+1) = dd
!an updated "f" at time t+1 is the output
 implicit none
 integer, intent(in):: m
 real(kind=kind_phys),dimension(m),intent(inout):: a,b,c
 real(kind=kind_phys),dimension(m),intent(inout):: f
!locals:
 integer:: k
 real(kind=kind_phys):: p
 real(kind=kind_phys),dimension(m):: q

 c(m)=0.
 q(1)=-c(1)/b(1)
 f(1)= f(1)/b(1)
 do k = 2,m
    p    = 1./(b(k)+a(k)*q(k-1))
    q(k) = -c(k)*p
    f(k) = p*(f(k)-a(k)*f(k-1))
 enddo
 do k = m-1,1,-1
    f(k) = f(k)+q(k)*f(k+1)
 enddo

 end subroutine tridiag

!=================================================================================================================
 subroutine tridiag2(n,a,b,c,d,x)
!=================================================================================================================
!copied as is from module_bl_mynn.F, as written by Joe Olson.
!a - sub-diagonal (means it is the diagonal below the main diagonal)
!b - the main diagonal
!c - sup-diagonal (means it is the diagonal above the main diagonal)
!d - right part
!x - the answer
!n - number of unknowns (levels)
 implicit none
 integer,intent(in):: n
 real(kind=kind_phys),dimension(n),intent(in):: a,b,c,d
 real(kind=kind_phys),dimension(n),intent(out):: x
 real(kind=kind_phys),dimension(n):: cp,dp
 real(kind=kind_phys):: m
 integer:: i

!initialize c-prime and d-prime
 cp(1) = c(1)/b(1)
 dp(1) = d(1)/b(1)
!solve for vectors c-prime and d-prime
 do i = 2,n
    m = b(i)-cp(i-1)*a(i)
    cp(i) = c(i)/m
    dp(i) = (d(i)-dp(i-1)*a(i))/m
 enddo
!initialize x
 x(n) = dp(n)
!solve for x from the vectors c-prime and d-prime
 do i = n-1,1,-1
    x(i) = dp(i)-cp(i)*x(i+1)
 end do

 end subroutine tridiag2

!=================================================================================================================
 end module module_cu_gf_ctrans
!=================================================================================================================
