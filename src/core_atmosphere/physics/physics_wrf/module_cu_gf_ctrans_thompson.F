!================================================================================================================
 module module_cu_gf_ctrans_thompson
!=================================================================================================================
 use mpas_kind_types
 use mpas_log
 use module_cu_gf_ctrans,only:tridiag

 implicit none
 private
 public:: ctrans_thompson


 contains


!=================================================================================================================
 subroutine ctrans_thompson(its,ite,kts,kte,ierr,k22,kbcon,ktop,kdn,nthom,use_thom_scaven,use_thom_evap,       &
                        use_thom_implicit,dtime,grav,edt,xmb,p,p_cup,z_cup,nifa,nwfa,pwavo,pwevo,pwo,pwdo,     &
                        zuo,up_massentro,up_massdetro,zdo,dn_massentro,dn_massdetro,nifa_dn_evap,nifa_up_scav, &
                        nwfa_dn_evap,nwfa_up_scav,nifa_dn_evapt,nifa_up_scavt,nwfa_dn_evapt,nwfa_up_scavt,     &
                        outnifa,outnwfa)
!=================================================================================================================

!--- input arguments:
 integer,intent(in):: its,ite,kts,kte

 integer,intent(in):: &
    nthom      ! number of aerosols that describe "water-friendly" and "ice-friendly" aerosols.

 integer,intent(in):: &
    use_thom_scaven,  &!
    use_thom_evap,    &!
    use_thom_implicit  !

 integer,intent(in),dimension(its:ite):: &
    ierr,     &!
    k22,      &!
    kbcon,    &!
    ktop,     &!
    kdn

 real(kind=RKIND),intent(in):: &
    dtime,    &!
    grav

 real(kind=RKIND),intent(in),dimension(its:ite):: &
    edt,      &!
    xmb,      &!
    pwavo,    &!
    pwevo      !

 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte):: &
    p,           &!
    p_cup,       &!
    z_cup         !

 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte):: &
    pwdo,        &!
    pwo           !

 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte):: &
    zuo,          &! normalized mass flux in updrafts               (-).
    up_massentro, &! normalized entrainment mass flux in updrafts   (-).
    up_massdetro, &! normalized detraiment mass flux in updrafts    (-).
    zdo,          &! normalized mass flux in downdrafts (-).
    dn_massentro, &! normalized entrainment mass flux in downdrafts (-).
    dn_massdetro   ! normalized detrainment mass flux in downdrafts (-).

 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte):: &
    nifa,         &! number concentration of "ice-friendly" aerosols   (# kg^{-1}).
    nwfa           ! number concentration of "water-friendly" aerosols (# kg^{-1}).


!--- inout arguments:
 real(kind=RKIND),intent(inout),dimension(its:ite):: &
    nifa_dn_evapt, &!
    nifa_up_scavt, &!
    nwfa_dn_evapt, &!
    nwfa_up_scavt   !

 real(kind=RKIND),intent(inout),dimension(its:ite,kts:kte):: &
    nifa_dn_evap, &!
    nifa_up_scav, &!
    nwfa_dn_evap, &!
    nwfa_up_scav, &!
    outnifa,      &! tendency of number concentration of "ice-friendly" aerosols   (# kg^{-1} s^{-1}).
    outnwfa        ! tendency of number concentration of "water-friendly" aerosols (# kg^{-1} s^{-1}).


!--- local variables and arrays:
 integer:: i,k,n
 integer:: k2,kb,kc,kd,kt

 real(kind=RKIND),parameter:: c0 = 0.004
 real(kind=RKIND):: denom,dp
 real(kind=RKIND),dimension(kts:kte):: fm,fp,znet
 real(kind=RKIND),dimension(its:ite):: pwav,pwev
 real(kind=RKIND),dimension(its:ite):: tot_nifa,tot_nifa_evap,tot_nifa_scav,tot_nwfa,tot_nwfa_evap,tot_nwfa_scav
 real(kind=RKIND),dimension(its:ite,kts:kte):: pw,pwd
 real(kind=RKIND),dimension(its:ite,kts:kte):: up_entr,up_detr,zu,dn_entr,dn_detr,zd
 real(kind=RKIND),dimension(its:ite,kts:kte):: nifal2,nifa_dn,nifa_up
 real(kind=RKIND),dimension(its:ite,kts:kte):: nifa_out_evap,nifa_out_scav,nifa_out
 real(kind=RKIND),dimension(its:ite,kts:kte):: nwfal2,nwfa_dn,nwfa_up
 real(kind=RKIND),dimension(its:ite,kts:kte):: nwfa_out_evap,nwfa_out_scav,nwfa_out

 real(kind=RKIND),dimension(nthom):: thom_wetscav,thom_conv

 real(kind=RKIND),dimension(:),allocatable:: mult,aterm,bterm,cterm,dterm

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine ctrans_thompson: nthom = $i',intArgs=(/nthom/))
 call mpas_log_write('--- use_thom_implicit = $i',intArgs=(/use_thom_implicit/))
 call mpas_log_write('--- use_thom_scaven   = $i',intArgs=(/use_thom_scaven/))
 call mpas_log_write('--- use_thom_evap     = $i',intArgs=(/use_thom_evap/))


!--- initialization:
 do k = kts,kte
    do i = its,ite
       nifa_dn(i,k)      = 0._RKIND
       nifa_up(i,k)      = 0._RKIND
       nifa_out(i,k)     = 0._RKIND
       nwfa_dn(i,k)      = 0._RKIND
       nwfa_up(i,k)      = 0._RKIND
       nwfa_out(i,k)     = 0._RKIND

       nifa_dn_evap(i,k) = 0._RKIND
       nifa_up_scav(i,k) = 0._RKIND
       nwfa_dn_evap(i,k) = 0._RKIND
       nwfa_up_scav(i,k) = 0._RKIND

       nifa_out_evap(i,k) = 0._RKIND
       nifa_out_scav(i,k) = 0._RKIND
       nwfa_out_evap(i,k) = 0._RKIND
       nwfa_out_scav(i,k) = 0._RKIND
    enddo
 enddo
 do i = its,ite
    nifa_dn_evapt(i) = 0._RKIND
    nifa_up_scavt(i) = 0._RKIND
    nwfa_dn_evapt(i) = 0._RKIND
    nwfa_up_scavt(i) = 0._RKIND

    tot_nifa(i)      = 0._RKIND
    tot_nifa_evap(i) = 0._RKIND
    tot_nifa_scav(i) = 0._RKIND
    tot_nwfa(i)      = 0._RKIND
    tot_nwfa_evap(i) = 0._RKIND
    tot_nwfa_scav(i) = 0._RKIND
 enddo


!--- conversion of all normalized mass fluxes, etc:
 do i = its,ite
    if(ierr(i) .eq. 0) then
       pwav(i) = xmb(i)*pwavo(i)
       pwev(i) = xmb(i)*pwevo(i)

       do k = kts,kte
          zu(i,k) = xmb(i)*zuo(i,k)
          zd(i,k) = xmb(i)*zdo(i,k)
          up_entr(i,k) = xmb(i)*up_massentro(i,k)
          up_detr(i,k) = xmb(i)*up_massdetro(i,k)
          dn_entr(i,k) = xmb(i)*dn_massentro(i,k)
          dn_detr(i,k) = xmb(i)*dn_massdetro(i,k)

          pw(i,k)  = xmb(i)*pwo(i,k)
          pwd(i,k) = xmb(i)*pwdo(i,k)
       enddo
    endif
 enddo


!--- fractional removal of cloud-borne aerosols for the individual aerosol species that are included in
!    the "water-friendly" and "ice-friendly" aerosols:
 thom_wetscav(1) = 0.2   ! DUST1 (ice-friendly).
 thom_wetscav(2) = 0.2   ! DUST2 (ice-friendly).
 thom_wetscav(3) = 0.2   ! DUST3 (ice-friendly).
 thom_wetscav(4) = 0.2   ! DUST4 (ice-friendly).
 thom_wetscav(5) = 0.2   ! DUST5 (ice-friendly).
 thom_wetscav(6) = 0.4   ! OC    (water-friendly).
 thom_wetscav(7) = 0.4   ! SO4   (water-friendly).
 thom_wetscav(8) = 0.4   ! SEAS1 (water-friendly).
 thom_wetscav(9) = 0.4   ! SEAS2 (water-friendly).


!--- computes the "ice-friendly", "water-friendly" aerosolos, and separate aerosols species at cloud
!    levels using an upstream approximation:
 do i = its,ite
    if(ierr(i) .eq. 0) then
       !environmental "ice-frienly" and "water-friendly" aerosols at cloud levels:
       k = kts
       nifal2(i,k) = nifa(i,k)
       nwfal2(i,k) = nwfa(i,k)
       do k = kts+1,kte
          nifal2(i,k) = 0.5*(nifa(i,k-1)+nifa(i,k))
          nwfal2(i,k) = 0.5*(nwfa(i,k-1)+nwfa(i,k))
       enddo
    endif
 enddo


!--- convective transport in convective_updrafts:
 call ctrans_thompson_up(its,ite,kts,kte,ierr,k22,kbcon,ktop,nthom,use_thom_scaven,c0,grav,thom_wetscav,p, &
                      p_cup,z_cup,nifa,nifal2,nwfa,nwfal2,zu,up_entr,up_detr,nifa_up,nwfa_up,nifa_up_scav, &  
                      nwfa_up_scav)


!--- convective transport in convective downdrafts:
 call ctrans_thompson_dn(its,ite,kts,kte,ierr,k22,kbcon,ktop,kdn,use_thom_evap,'deep',z_cup,nifa, &
                         nifa_up_scav,nwfa,nwfa_up_scav,zd,dn_entr,dn_detr,pwav,pwev,pwd,nifa_dn, &
                         nifa_dn_evap,nwfa_dn,nwfa_dn_evap)


!--- computation budget terms due to wet scavenging in updrafts and evaporation in downdrafts:
 do i = its,ite
    if(ierr(i) .eq. 0) then
       k2 = k22(i)
       kb = kbcon(i)
       kt = ktop(i)
       kd = kdn(i)

       if(use_thom_scaven .gt. 0) then
          do k = kts,kt
             nifa_up_scavt(i) = nifa_up_scavt(i) - zu(i,k)*nifa_up_scav(i,k)
             nwfa_up_scavt(i) = nwfa_up_scavt(i) - zu(i,k)*nwfa_up_scav(i,k)

             dp = 0.01_RKIND*grav/(p_cup(i,k)-p_cup(i,k+1))
             nifa_out_scav(i,k) = - 0.5*dp*(zu(i,k)*nifa_up_scav(i,k)+zu(i,k+1)*nifa_up_scav(i,k+1))
             nifa_out(i,k) = nifa_out_scav(i,k)

             nwfa_out_scav(i,k) = - 0.5*dp*(zu(i,k)*nwfa_up_scav(i,k)+zu(i,k+1)*nwfa_up_scav(i,k+1))
             nwfa_out(i,k) = nwfa_out_scav(i,k)
          enddo
       endif
       if(use_thom_evap .gt. 0) then
          do k = kts,kt
             nifa_dn_evapt(i) = nifa_dn_evapt(i) + edt(i)*zd(i,k)*nifa_dn_evap(i,k)
             nwfa_dn_evapt(i) = nwfa_dn_evapt(i) + edt(i)*zd(i,k)*nwfa_dn_evap(i,k)

             dp = 0.01_RKIND*grav/(p_cup(i,k)-p_cup(i,k+1))
             nifa_out_evap(i,k) = - 0.5*dp*edt(i)*(zd(i,k)*nifa_dn_evap(i,k)+zd(i,k+1)*nifa_dn_evap(i,k+1))
             nifa_out(i,k) = nifa_out(i,k) + nifa_out_evap(i,k)

             nwfa_out_evap(i,k) = - 0.5*dp*edt(i)*(zd(i,k)*nwfa_dn_evap(i,k)+zd(i,k+1)*nwfa_dn_evap(i,k+1))
             nwfa_out(i,k) = nwfa_out(i,k) + nwfa_out_evap(i,k)
          enddo
       endif
!      do k = kts,kte
!         call mpas_log_write('$i $i $i $i $i $i $r $r $r $r $r $r $r $r $r $r $r $r',intArgs=(/i,k2,kb,kd,kt,k/), &
!                            realArgs=(/nifal2(i,k),nifa_up(i,k),nifa_dn(i,k),nifa_up_scav(i,k),nifa_dn_evap(i,k), &
!                            nifa_out(i,k),nwfal2(i,k),nwfa_up(i,k),nwfa_dn(i,k),nwfa_up_scav(i,k),nwfa_dn_evap(i,k), &
!                            nwfa_out(i,k)/))
!      enddo
!      call mpas_log_write(' ')
!      exit

    endif
 enddo


!--- computation of the tendencies of "ice-friendly" and "water-friendly" aerosols:
 if(use_thom_implicit .eq. 0) then
    do i = its,ite
       if(ierr(i) .eq. 0) then
          kt = ktop(i)

          do k = kts,kt
             dp = 0.01_RKIND*grav/(p_cup(i,k)-p_cup(i,k+1))
             outnifa(i,k) = - dp*(zu(i,k+1)*(nifa_up(i,k+1)-nifal2(i,k+1))-zu(i,k)*(nifa_up(i,k)-nifal2(i,k))) &
                     + edt(i)*dp*(zd(i,k+1)*(nifa_dn(i,k+1)-nifal2(i,k+1))-zd(i,k)*(nifa_dn(i,k)-nifal2(i,k)))
             outnifa(i,k) = outnifa(i,k) + nifa_out(i,k)

             outnwfa(i,k) = - dp*(zu(i,k+1)*(nwfa_up(i,k+1)-nwfal2(i,k+1))-zu(i,k)*(nwfa_up(i,k)-nwfal2(i,k))) &
                     + edt(i)*dp*(zd(i,k+1)*(nwfa_dn(i,k+1)-nwfal2(i,k+1))-zd(i,k)*(nwfa_dn(i,k)-nwfal2(i,k)))
             outnwfa(i,k) = outnwfa(i,k) + nwfa_out(i,k)
          enddo
!         do k = kts,kte
!            call mpas_log_write('$i $i $i $r $r $r $r',intArgs=(/i,kt,k/),realArgs=(/nifa_out(i,k),outnifa(i,k), &
!                                nwfa_out(i,k),outnwfa(i,k)/))
!         enddo
!         exit

       endif
    enddo

 elseif(use_thom_implicit .gt. 0) then
    do i = its,ite
       if(ierr(i) .eq. 0) then
          k2 = k22(i)
          kb = kbcon(i)
          kt = ktop(i)
          kd = kdn(i)
          kc = kt-kts+1 
          do k = kts,kte
             znet(k) = zu(i,k)-edt(i)*zd(i,k)
             fp(k)   = 0.5*(znet(k)+abs(znet(k)))
             fm(k)   = 0.5*(znet(k)-abs(znet(k)))
          enddo

          if(.not.allocated(mult) ) allocate(mult(kts:kt) )
          if(.not.allocated(aterm)) allocate(aterm(kts:kt))
          if(.not.allocated(bterm)) allocate(bterm(kts:kt))
          if(.not.allocated(cterm)) allocate(cterm(kts:kt))
          if(.not.allocated(dterm)) allocate(dterm(kts:kt))
          do k = kts,kt
             mult(k)  = 0.01_RKIND*grav/(p_cup(i,k)-p_cup(i,k+1))
             aterm(k) = mult(k)*dtime*fm(k)
             bterm(k) = 1.+mult(k)*dtime*(fp(k)-fm(k+1))
             cterm(k) = -mult(k)*dtime*fp(k+1)
          enddo

          !--- "ice-friendly" aerosols:
          do k = kts,kt
             dterm(k) = (zu(i,k)*nifa_up(i,k)-zu(i,k+1)*nifa_up(i,k+1)) &
                      - edt(i)*(zd(i,k)*nifa_dn(i,k)-zd(i,k+1)*nifa_dn(i,k+1))
             dterm(k) = nifa(i,k) + mult(k)*dtime*dterm(k)
          enddo
          if(use_thom_scaven .gt. 0) then
             do k = kts,kt
                dterm(k) = dterm(k) + dtime*nifa_out_scav(i,k)
             enddo
          endif
          if(use_thom_evap) then
             do k = kts,kt
                dterm(k) = dterm(k) + dtime*nifa_out_evap(i,k)
             enddo
          endif
          call tridiag(kc,aterm,bterm,cterm,dterm)
          do k = kts,kt
             outnifa(i,k) = (dterm(k)-nifa(i,k))/dtime
          enddo

          !--- "water-friendly" aerosols:
          do k = kts,kt
             dterm(k) = (zu(i,k)*nwfa_up(i,k)-zu(i,k+1)*nwfa_up(i,k+1)) &
                      - edt(i)*(zd(i,k)*nwfa_dn(i,k)-zd(i,k+1)*nwfa_dn(i,k+1))
             dterm(k) = nwfa(i,k) + mult(k)*dtime*dterm(k)
          enddo
          if(use_thom_scaven > 0) then
             do k = kts,kt
                dterm(k) = dterm(k) + dtime*nwfa_out_scav(i,k)
             enddo
          endif
          if(use_thom_evap) then
             do k = kts,kt
                dterm(k) = dterm(k) + dtime*nwfa_out_evap(i,k)
             enddo
          endif
          call tridiag(kc,aterm,bterm,cterm,dterm)
          do k = kts,kt
             outnwfa(i,k) = (dterm(k)-nwfa(i,k))/dtime
          enddo

          if(allocated(mult) ) deallocate(mult )
          if(allocated(aterm)) deallocate(aterm)
          if(allocated(bterm)) deallocate(bterm)
          if(allocated(cterm)) deallocate(cterm)
          if(allocated(dterm)) deallocate(dterm)

       endif
    enddo
 endif


!--- conservation of "ice-friendly" and "water-friendly" Thompson aerosols:
!call mpas_log_write(' ')
!call mpas_log_write('--- conservation of number concentration:')
!do i = its,ite
!   if(ierr(i) .eq. 0) then
!      kt = ktop(i)
!      do k = kts,kt
!         dp = 100_RKIND*(p_cup(i,k)-p_cup(i,k+1))/grav
!         tot_nifa(i) = tot_nifa(i) + dp*outnifa(i,k)
!         tot_nwfa(i) = tot_nwfa(i) + dp*outnwfa(i,k)

!         tot_nifa_scav(i) = tot_nifa_scav(i) &
!                          - 0.5*(zu(i,k)*nifa_up_scav(i,k)+zu(i,k+1)*nifa_up_scav(i,k+1))
!         tot_nwfa_scav(i) = tot_nwfa_scav(i)-0.5*(zu(i,k)*nwfa_up_scav(i,k)+zu(i,k+1)*nwfa_up_scav(i,k+1))

!         tot_nifa_evap(i) = tot_nifa_evap(i) &
!                          - 0.5*edt(i)*(zd(i,k)*nifa_dn_evap(i,k)+zd(i,k+1)*nifa_dn_evap(i,k+1))
!         tot_nwfa_evap(i) = tot_nwfa_evap(i) &
!                          - 0.5*edt(i)*(zd(i,k)*nwfa_dn_evap(i,k)+zd(i,k+1)*nwfa_dn_evap(i,k+1))
!      enddo
!      call mpas_log_write('$i  --- tot_nifa = $r  tot_nifa_scav = $r  tot_nifa_evap = $r',intArgs=(/i/), &
!                          realArgs=(/tot_nifa(i),tot_nifa_scav(i),tot_nifa_evap(i)/))
!      call mpas_log_write('$i  --- tot_nwfa = $r  tot_nwfa_scav = $r  tot_nwfa_evap = $r',intArgs=(/i/), &
!                          realArgs=(/tot_nwfa(i),tot_nwfa_scav(i),tot_nwfa_evap(i)/))
!   endif
!enddo


 call mpas_log_write('--- end subroutine ctrans_thompson:')
 call mpas_log_write(' ')

 end subroutine ctrans_thompson

!=================================================================================================================
 subroutine ctrans_thompson_up(its,ite,kts,kte,ierr,k22,kbcon,ktop,nthom,use_thom_scaven,c0,grav,fscav,p,p_cup, &
                               z_cup,nifa,nifal2,nwfa,nwfal2,zu,up_entr,up_detr,nifa_up,nwfa_up,nifa_up_scav,   &
                               nwfa_up_scav)
!=================================================================================================================

!--- input arguments:
 integer,intent(in):: its,ite,kts,kte

 integer,intent(in):: &
    nthom,         &! number of aerosols that describe "water-friendly" and "ice-friendly" aerosols.
    use_thom_scaven !
                              
 integer,intent(in),dimension(its:ite):: &
    ierr,          &!
    k22,           &!
    kbcon,         &!
    ktop            !

 real(kind=RKIND),intent(in):: &
    c0,            &!
    grav            !

 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte):: &
    p,             &!
    p_cup,         &!
    z_cup

 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte):: &
    up_entr,       &!
    up_detr,       &!
    zu              !

 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte):: &
    nifa,          &!
    nifal2,        &!
    nwfa,          &!
    nwfal2          !

 real(kind=RKIND),intent(in),dimension(nthom):: &
    fscav           !


!--- inout arguments:
 real(kind=RKIND),intent(inout),dimension(its:ite,kts:kte):: &
    nifa_up,     &!
    nwfa_up,     &!
    nifa_up_scav,&!
    nwfa_up_scav  !


!--- local variables and arrays:
 integer:: i,k,n
 integer:: k2,kb,kt

 real(kind=RKIND):: fscav_i,fscav_w
 real(kind=RKIND):: denom,dp,dz,thom_b,thom_c,thom_g,thom_l,thom_r
 real(kind=RKIND):: coef1,coef2

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine ctrans_thompson_up: nthom = $i',intArgs=(/nthom/))

 fscav_i = maxval(fscav(1:5))
 fscav_w = maxval(fscav(6:nthom))
 call mpas_log_write('--- fscav_i = $r',realArgs=(/fscav_i/))
 call mpas_log_write('--- fscav_w = $r',realArgs=(/fscav_w/))

 do i = its,ite
    if(ierr(i) .eq. 0) then
       k2 = k22(i)
       kb = kbcon(i)
       kt = ktop(i)

       !--- convective transport of "ice-friendly" and "water-friendly" aerosols in updrafts (use upstream
       !    method to define the in-cloud variables at the k2 levels):
       k = k2
       nifa_up(i,k) = nifa(i,k-1)
       nwfa_up(i,k) = nwfa(i,k-1)
       do k = k2+1,kt+1
          dz = z_cup(i,k)-z_cup(i,k-1)
          denom = zu(i,k-1)-.5*up_detr(i,k-1)+up_entr(i,k-1)
          if(denom .lt. 1.e-16) exit

          !transport in updrafts:
          nifa_up(i,k) = (nifa_up(i,k-1)*zu(i,k-1) - .5*up_detr(i,k-1)*nifa_up(i,k-1) &
                       + up_entr(i,k-1)*nifa(i,k-1))/denom
          nwfa_up(i,k) = (nwfa_up(i,k-1)*zu(i,k-1) - .5*up_detr(i,k-1)*nwfa_up(i,k-1) &
                       + up_entr(i,k-1)*nwfa(i,k-1))/denom

          !wet scavenging:
          !formulation as in GOCART (using sourcecode in ConvPar_GF_GEOS5 provided by Saulo Freitas):
          if(use_thom_scaven == 1) then
             nifa_up_scav(i,k) = max(0.,nifa_up(i,k)*(1.-exp(-fscav_i*(dz/1000.))))
             nwfa_up_scav(i,k) = max(0.,nwfa_up(i,k)*(1.-exp(-fscav_w*(dz/1000.))))

             nifa_up(i,k) = nifa_up(i,k) - nifa_up_scav(i,k)
             nwfa_up(i,k) = nwfa_up(i,k) - nwfa_up_scav(i,k)
          endif
       enddo
    endif
 enddo

 call mpas_log_write('--- end subroutine ctrans_thompson_up:')

 end subroutine ctrans_thompson_up

!=================================================================================================================
 subroutine ctrans_thompson_dn(its,ite,kts,kte,ierr,k22,kbcon,ktop,kdn,use_thom_evap,cumulus,z_cup,nifa, &
                               nifa_up_scav,nwfa,nwfa_up_scav,zd,dn_entr,dn_detr,pwav,pwev,pwd,nifa_dn,  &
                               nifa_dn_evap,nwfa_dn,nwfa_dn_evap)
!=================================================================================================================

!--- input arguments:
 character(len=*),intent(in):: &
    cumulus        !

 integer,intent(in):: its,ite,kts,kte

 integer,intent(in):: &
    use_thom_evap  !                       

 integer,intent(in),dimension(its:ite):: &
    ierr,         &!
    k22,          &!
    kbcon,        &!
    ktop,         &!
    kdn            !

 real(kind=RKIND),intent(in),dimension(its:ite):: &
    pwav,         &!total rain in each updraft.
    pwev           !total evaporation in each downdraft.

 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte):: &
    z_cup ,       &!
    pwd            !evaporation in each layer of updraft.

 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte):: &
    dn_entr,      &!
    dn_detr,      &!
    zd             !

 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte):: &
    nifa,         &!
    nwfa,         &!
    nifa_up_scav, &!
    nwfa_up_scav   !


!--- inout arguments:
 real(kind=RKIND),intent(inout),dimension(its:ite,kts:kte):: &
    nifa_dn,      &!
    nwfa_dn,      &!
    nifa_dn_evap, &!
    nwfa_dn_evap   !


!--- local variables and arrays:
 integer:: i,k,k2,kb,kd,kt

 real(kind=RKIND):: frac_evap,pwdper
 real(kind=RKIND):: denom,dz

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine ctrans_thompson_dn:')

 if(cumulus == 'shallow') return

 do i = its,ite
    if(ierr(i) /= 0) cycle
    k2 = k22(i)
    kb = kbcon(i)
    kd = kdn(i)
    kt = ktop(i)


    !--- fraction of the total rain that was evaporated:
    frac_evap = - pwev(i)/(1.e-16+pwav(i))

    !--- initializes downdrafts at kd level with a mixture of 90% environmental and 10% in-cloud mixing ratios,
    !    or initializes downdrafts using upstream approximation:
    k = kd
!   nifa_dn(i,k) = 0.9*nifa(i,k) + 0.1*nifa_up(i,k)
!   nwfa_dn(i,k) = 0.9*nwfa(i,k) + 0.1*nwfa_up(i,k)
    nifa_dn(i,k) = nifa(i,k)
    nwfa_dn(i,k) = nwfa(i,k)

    pwdper = pwd(i,k)/(1.e-16+pwev(i))*frac_evap ! > 0
!   pwdper = pwd(i,k)/(1.e-16+pwav(i))
    if(use_thom_evap .eq. 0) pwdper = 0._RKIND
    nifa_dn_evap(i,k) = - pwdper*nifa_up_scav(i,k)
    nwfa_dn_evap(i,k) = - pwdper*nwfa_up_scav(i,k)
    nifa_dn(i,k) = nifa_dn(i,k) - nifa_dn_evap(i,k)
    nwfa_dn(i,k) = nwfa_dn(i,k) - nwfa_dn_evap(i,k)
!   call mpas_log_write('$i $i $i $i $i $i $r $r $r $r $r $r $r',intArgs=(/i,k2,kb,kd,kt,k/),realArgs= &
!            (/pwdper,nifa_dn(i,k),nwfa_dn(i,k),nifa_up_scav(i,k),nwfa_up_scav(i,k),nifa_dn_evap(i,k), &
!            nwfa_dn_evap(i,k)/))

    !--- calculates the downdraft mass terms:
    do k = kd-1,kts,-1
       dz = z_cup(i,k+1)-z_cup(i,k)
       denom = zd(i,k+1)-.5*dn_detr(i,k)+dn_entr(i,k)
       if(denom > 0._RKIND) then
          nifa_dn(i,k) = (nifa_dn(i,k+1)*zd(i,k+1) - .5*dn_detr(i,k)*nifa_dn(i,k+1) &
                       + dn_entr(i,k)*nifa(i,k))/denom

          nwfa_dn(i,k) = (nwfa_dn(i,k+1)*zd(i,k+1) - .5*dn_detr(i,k)*nwfa_dn(i,k+1) &
                       +  dn_entr(i,k)*nwfa(i,k))/denom
       else
          nifa_dn(i,k) = nifa(i,k+1)
          nwfa_dn(i,k) = nwfa(i,k+1)
       endif

       !--- evaporation:
       pwdper = pwd(i,k)/(1.e-16+pwev(i))*frac_evap ! > 0
!      pwdper = pwd(i,k)/(1.e-16+pwav(i))
       if(use_thom_evap .eq. 0) pwdper = 0._RKIND
       nifa_dn_evap(i,k) = - pwdper*nifa_up_scav(i,k)
       nwfa_dn_evap(i,k) = - pwdper*nwfa_up_scav(i,k)
       nifa_dn(i,k) = nifa_dn(i,k) - nifa_dn_evap(i,k)
       nwfa_dn(i,k) = nwfa_dn(i,k) - nwfa_dn_evap(i,k)
!      call mpas_log_write('$i $i $i $i $i $i $r $r $r $r $r $r $r',intArgs=(/i,k2,kb,kd,kt,k/),realArgs= &
!               (/pwdper,nifa_dn(i,k),nwfa_dn(i,k),nifa_up_scav(i,k),nwfa_up_scav(i,k),nifa_dn_evap(i,k), &
!               nwfa_dn_evap(i,k)/))
    enddo
!   call mpas_log_write(' ')
 enddo

 call mpas_log_write('--- end subroutine ctrans_thompson_dn:')

 end subroutine ctrans_thompson_dn

!=================================================================================================================
 end module module_cu_gf_ctrans_thompson
!=================================================================================================================
