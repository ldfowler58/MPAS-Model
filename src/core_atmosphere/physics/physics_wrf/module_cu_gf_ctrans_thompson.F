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
 subroutine ctrans_thompson(its,ite,kts,kte,ierr,k22,kbcon,ktop,kdn,nthom,use_thom_scaven,use_thom_evap,dtime, &
                        grav,edt,xmb,p,p_cup,z_cup,thom_w,nifa,nwfa,thom_scav,nifa_scav,nwfa_scav,pwavo,pwevo, &
                        pwo,pwdo,zu,up_entr,up_detr,zd,dn_entr,dn_detr,outnifa,outnwfa)
!=================================================================================================================

!--- input arguments:
 integer,intent(in):: its,ite,kts,kte

 integer,intent(in):: &
    nthom      ! number of aerosols that describe "water-friendly" and "ice-friendly" aerosols.

 integer,intent(in):: &
    use_thom_scaven, &!
    use_thom_evap     !

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
    p,        &!
    p_cup,    &!
    z_cup      !

 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte):: &
    pwdo,     &!
    pwo        !

 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte):: &
    up_entr,  &!
    up_detr,  &!
    zu,       &!
    dn_entr,  &!
    dn_detr,  &!
    zd         !

 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte):: &
    nifa,     &! "ice-friendly" aerosols   (# per kg).
    nwfa       ! "water-friendly" aerosols (# per kg).

 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte,nthom):: &
    thom_w     ! weight of individual aerosols to "water-friendly" and "ice-friendly" aerosols (-).


!--- inout arguments:
 real(kind=RKIND),intent(inout),dimension(its:ite,kts:kte):: &
    outnifa,  &! "ice-friendly" aerosols   (# per kg).
    outnwfa    ! "water-friendly" aerosols (# per kg).

 real(kind=RKIND),intent(inout),dimension(its:ite,kts:kte):: &
    nifa_scav,&!
    nwfa_scav  !

 real(kind=RKIND),intent(inout),dimension(its:ite,kts:kte,nthom):: &
    thom_scav  !


!--- local variables and arrays:
 integer:: i,k,n
 integer:: k2,kb,kc,kd,kt

 real(kind=RKIND),parameter:: c0 = 0.004
 real(kind=RKIND):: denom
 real(kind=RKIND),dimension(kts:kte):: fm,fp,znet
 real(kind=RKIND),dimension(kts:kte):: aeri,aeri_dn,aeri_up,aerw,aerw_dn,aerw_up
 real(kind=RKIND),dimension(its:ite):: nifa_up_scav_tot,nwfa_up_scav_tot
 real(kind=RKIND),dimension(its:ite,kts:kte):: nifal2,nifa_dn,nifa_up,nifa_dn_scav,nifa_up_scav,nifa_out
 real(kind=RKIND),dimension(its:ite,kts:kte):: nwfal2,nwfa_dn,nwfa_up,nwfa_dn_scav,nwfa_up_scav,nwfa_out

 real(kind=RKIND),dimension(nthom):: thom_wetscav,thom_conv
 real(kind=RKIND),dimension(its:ite,nthom):: thom_up_scav_tot
 real(kind=RKIND),dimension(its:ite,kts:kte,nthom):: thom,thoml2
 real(kind=RKIND),dimension(its:ite,kts:kte,nthom):: thom_dn,thom_dn_scav,thom_up,thom_up_scav,thom_out
 real(kind=RKIND),dimension(its:ite,kts:kte,nthom):: thom_mult,thom_multo

 real(kind=RKIND),dimension(:),allocatable:: mult,aterm,bterm,cterm,dterm

 real(kind=RKIND),dimension(kts:kte):: trashi,trashw

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine ctrans_thompson: nthom = $i',intArgs=(/nthom/))


!--- allocate the different aerosol components of "water-friendly" and "ice-friendly" aerosols:
 do n = 1,nthom
    do k = kts,kte
       do i = its,ite
          thom_up(i,k,n)   = 0._RKIND
          thom_dn(i,k,n)   = 0._RKIND
          thom_out(i,k,n)  = 0._RKIND
          thom_scav(i,k,n) = 0._RKIND
          thom_dn_scav(i,k,n) = 0._RKIND
          thom_up_Scav(i,k,n) = 0._RKIND
       enddo
    enddo
    do i = its,ite
       thom_up_scav_tot(i,n) = 0._RKIND
    enddo
 enddo

 do k = kts,kte
    do i = its,ite
       nifa_dn(i,k)   = 0._RKIND
       nifa_up(i,k)   = 0._RKIND
       nifa_out(i,k)  = 0._RKIND
       nwfa_dn(i,k)   = 0._RKIND
       nwfa_up(i,k)   = 0._RKIND
       nwfa_out(i,k)  = 0._RKIND

       nifa_dn_scav(i,k) = 0._RKIND
       nifa_up_scav(i,k) = 0._RKIND
       nwfa_dn_scav(i,k) = 0._RKIND
       nwfa_up_scav(i,k) = 0._RKIND
    enddo
 enddo
 do i = its,ite
    nifa_up_scav_tot(i) = 0._RKIND
    nwfa_up_scav_tot(i) = 0._RKIND
 enddo


!--- conversion coefficients from kg/kg to #/kg:
 thom_conv(1) = 1.58e-15 ! DUST1.
 thom_conv(2) = 1.68e-14 ! DUST2.
 thom_conv(3) = 8.46e-14 ! DUST3.
 thom_conv(4) = 5.58e-13 ! DUST4.
 thom_conv(5) = 3.13e-13 ! DUST5.
 thom_conv(6) = 6.61e-18 ! OC.
 thom_conv(7) = 1.07e-17 ! SO4.
 thom_conv(8) = 1.04e-15 ! SEAS1.
 thom_conv(9) = 3.88e-13 ! SEAS2.


!--- fractional removal of cloud-borne aerosols for the individual aerosol sepcies that are included in
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


!--- computes the different "environmental" aerosols that describe "water-friendly" and "ice-friendly"
!    aerosols:
 do i = its,ite
    do k = kts,kte
       do n = 1,nthom
          thom_mult(i,k,n)  = thom_w(i,k,n)
!         thom_multo(i,k,n) = 1./thom_mult(i,k,n)
       enddo
    enddo
    do k = kts,kte
       !aerosols that describe "ice-friendly" aerosols:
       do n = 1,5
          thom(i,k,n) = thom_mult(i,k,n)*nifa(i,k)
       enddo
       !aerosols that describe "water-friendly" aerosols:
       do n = 6,nthom
          thom(i,k,n) = thom_mult(i,k,n)*nwfa(i,k)
       enddo
    enddo
 enddo


!--- computes the "ice-friendly", "water-friendly" aerosolos, and separate aerosols species at cloud
!    levels using an upstream approximation:
 do i = its,ite
    if(ierr(i) .eq. 0) then
       !environmental "ice-frienly" and "water-friendly" aerosols at cloud levels:
       k = kts
       nifal2(i,k) = nifa(i,k)
       nwfal2(i,k) = nwfa(i,k)
       do k = kts+1,kte
          nifal2(i,k) = nifa(i,k-1)
          nwfal2(i,k) = nwfa(i,k-1)
       enddo
       !environmental aerosols at cloud levels:
       k = kts
       do n = 1,nthom
          thoml2(i,k,n) = thom(i,k,n)
       enddo
       do k = kts+1,kte
          do n = 1,nthom
             thoml2(i,k,n) = thom(i,k-1,n)
          enddo
       enddo
    endif
 enddo


!--- convective transport in convective_updrafts:
 call ctrans_thompson_up(its,ite,kts,kte,ierr,k22,kbcon,ktop,nthom,use_thom_scaven,c0,grav,thom_wetscav,p,   &
                         p_cup,z_cup,thom,thoml2,nifa,nifal2,nwfa,nwfal2,zu,up_entr,up_detr,thom_up,nifa_up, &
                         nwfa_up,thom_up_scav,nifa_up_scav,nwfa_up_scav,thom_up_scav_tot,nifa_up_scav_tot,   &
                         nwfa_up_scav_tot)


!--- convective transport in convective downdrafts:
 call ctrans_thompson_dn(its,ite,kts,kte,ierr,k22,kbcon,ktop,kdn,nthom,use_thom_evap,'deep',grav,thom_wetscav,    &
                         p,p_cup,z_cup,thom,thoml2,thom_up,nifa,nifal2,nifa_up,nwfa,nwfal2,nwfa_up,zd,dn_entr,    &
                         dn_detr,pwavo,pwevo,pwo,pwdo,thom_up_scav_tot,nifa_up_scav_tot,nwfa_up_scav_tot,thom_dn, &
                         thom_dn_scav,nifa_dn,nifa_dn_scav,nwfa_dn,nwfa_dn_scav)


!--- updates the "water-friendly" and "ice-friendly" aerosols, and computes tendencies of "water-friendly"
!    and "ice-friendly" aerosols due transport,activation,and wet scavenging in convective updrafts and
!    downdrafts.
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

          aeri(k)    = 0._RKIND
          aerw(k)    = 0._RKIND
          aeri_up(k) = 0._RKIND
          aerw_up(k) = 0._RKIND
          aeri_dn(k) = 0._RKIND
          aerw_dn(k) = 0._RKIND
          do n = 1,5
             aeri_dn(k) = aeri_dn(k) + thom_dn(i,k,n)
             aeri_up(k) = aeri_up(k) + thom_up(i,k,n)
          enddo
          do n = 6,nthom
             aerw_dn(k) = aerw_dn(k) + thom_dn(i,k,n)
             aerw_up(k) = aerw_up(k) + thom_up(i,k,n)
          enddo
       enddo

       if(.not.allocated(mult) ) allocate(mult(kts:kt) )
       if(.not.allocated(aterm)) allocate(aterm(kts:kt))
       if(.not.allocated(bterm)) allocate(bterm(kts:kt))
       if(.not.allocated(cterm)) allocate(cterm(kts:kt))
       if(.not.allocated(dterm)) allocate(dterm(kts:kt))
       do k = kts,kt
          mult(k)  = 0.01_RKIND*grav/(p_cup(i,k)-p_cup(i,k+1))
          mult(k)  = xmb(i)*mult(k)

          aterm(k) = mult(k)*dtime*fm(k)
          bterm(k) = 1.+mult(k)*dtime*(fp(k)-fm(k+1))
          cterm(k) = -mult(k)*dtime*fp(k+1)
       enddo

       !--- "ice-friendly" aerosols:
       do k = kts,kt
          dterm(k) = (zu(i,k)*aeri_up(k)-zu(i,k+1)*aeri_up(k+1)) &
                   - edt(i)*(zd(i,k)*aeri_dn(k)-zd(i,k+1)*aeri_dn(k+1))
          dterm(k) = nifa(i,k) + mult(k)*dtime*dterm(k)
       enddo
       
       call tridiag(kc,aterm,bterm,cterm,dterm)
       do k = kts,kt
          nifa_scav(i,k) = dterm(k)
       enddo
       do k = kt+1,kte
          nifa_scav(i,k) = nifa(i,k)
       enddo

       do k = kts,kt
          dterm(k) = (zu(i,k)*nifa_up(i,k)-zu(i,k+1)*nifa_up(i,k+1)) &
                   - edt(i)*(zd(i,k)*nifa_dn(i,k)-zd(i,k+1)*nifa_dn(i,k+1))
          dterm(k) = nifa(i,k) + mult(k)*dtime*dterm(k)
       enddo
       if(use_thom_scaven == 1) then
          do k = kts,kt
             nifa_out(i,k) = nifa_out(i,k) &
                 - 0.5*(zu(i,k)*nifa_up_scav(i,k)+zu(i,k+1)*nifa_up_scav(i,k+1))*mult(k)*dtime
             dterm(k) = dterm(k) + nifa_out(i,k)
          enddo
       endif
       if(use_thom_evap == 1) then
          do k = kts,kt
             nifa_out(i,k) = nifa_out(i,k) &
                 - 0.5*(zd(i,k)*nifa_dn_scav(i,k)+zd(i,k+1)*nifa_dn_scav(i,k+1))*edt(i)**mult(k)*dtime
             dterm(k) = dterm(k) + nifa_out(i,k)
          enddo
       endif
       call tridiag(kc,aterm,bterm,cterm,dterm)
       do k = kts,kt
          outnifa(i,k) = (dterm(k)-nifa(i,k))/dtime
       enddo

       !--- "water-friendly" aerosols:
       do k = kts,kt
          dterm(k) = (zu(i,k)*aerw_up(k)-zu(i,k+1)*aerw_up(k+1)) &
                   - edt(i)*(zd(i,k)*aerw_dn(k)-zd(i,k+1)*aerw_dn(k+1))
          dterm(k) = nwfa(i,k) + mult(k)*dtime*dterm(k)
       enddo
       call tridiag(kc,aterm,bterm,cterm,dterm)
       do k = kts,kt
          nwfa_scav(i,k) = dterm(k)
       enddo
       do k = kt+1,kte
          nwfa_scav(i,k) = nwfa(i,k)
       enddo

       do k = kts,kt
          dterm(k) = (zu(i,k)*nwfa_up(i,k)-zu(i,k+1)*nwfa_up(i,k+1)) &
                   - edt(i)*(zd(i,k)*nwfa_dn(i,k)-zd(i,k+1)*nwfa_dn(i,k+1))
          dterm(k) = nwfa(i,k) + mult(k)*dtime*dterm(k)
       enddo
       if(use_thom_scaven == 1) then
          do k = kts,kt
             nwfa_out(i,k) = nwfa_out(i,k) &
                 - 0.5*(zu(i,k)*nwfa_up_scav(i,k)+zu(i,k+1)*nwfa_up_scav(i,k+1))*mult(k)*dtime
             dterm(k) = dterm(k) + nwfa_out(i,k)
          enddo
       endif
       if(use_thom_evap == 1) then
          do k = kts,kt
             nwfa_out(i,k) = nwfa_out(i,k) &
                 - 0.5*(zd(i,k)*nwfa_dn_scav(i,k)+zd(i,k+1)*nwfa_dn_scav(i,k+1))*edt(i)**mult(k)*dtime
             dterm(k) = dterm(k) + nwfa_out(i,k)
          enddo
       endif
       call tridiag(kc,aterm,bterm,cterm,dterm)
       do k = kts,kt
          outnwfa(i,k) = (dterm(k)-nwfa(i,k))/dtime
       enddo

       !--- separate "Thompson" aerosols:
       do n = 1,nthom
          do k = kts,kt
             dterm(k) = (zu(i,k)*thom_up(i,k,n)-zu(i,k+1)*thom_up(i,k+1,n)) &
                      - edt(i)*(zd(i,k)*thom_dn(i,k,n)-zd(i,k+1)*thom_dn(i,k+1,n))
             dterm(k) = thom(i,k,n) + mult(k)*dtime*dterm(k)
          enddo
          !include scavenging:
          if(use_thom_scaven == 1) then
             do k = kts,kt
                thom_out(i,k,n) = thom_out(i,k,n) &
                    - 0.5*(zu(i,k)*thom_up_scav(i,k,n)+zu(i,k+1)*thom_up_scav(i,k+1,n))*mult(k)*dtime
                dterm(k) = dterm(k) + thom_out(i,k,n)
             enddo
          endif
          if(use_thom_scaven == 1) then
             do k = kts,kt
                thom_out(i,k,n) = thom_out(i,k,n) &
                    - 0.5*(zd(i,k)*thom_dn_scav(i,k,n)+zd(i,k+1)*thom_dn_scav(i,k+1,n))*edt(i)**mult(k)*dtime
                dterm(k) = dterm(k) + thom_out(i,k,n)
             enddo
          endif
          call tridiag(kc,aterm,bterm,cterm,dterm)
          do k = kts,kt
             thom_scav(i,k,n) = dterm(k)
          enddo
          do k = kt+1,kte
             thom_scav(i,k,n) = thom(i,k,n)
          enddo
       enddo

       if(allocated(mult) ) deallocate(mult )
       if(allocated(aterm)) deallocate(aterm)
       if(allocated(bterm)) deallocate(bterm)
       if(allocated(cterm)) deallocate(cterm)
       if(allocated(dterm)) deallocate(dterm)

    endif
 enddo


!--- conservation of mass:
!if(.not.allocated(mult)) allocate(mult(kts:kte))
!do i = its,ite
!   kt = ktop(i)
!   if(ierr(i) .eq. 0) then
!      trash = 0._RKIND
!      do k = kts,kt
!         mult(k)  = 0.01_RKIND*grav/(p_cup(i,k)-p_cup(i,k+1))
!         trash = trash + outnifa(i,k)/mult(k) 
!      enddo
!      call mpas_log_write('$i $r',intArgs=(/i/),realArgs=(/trash/))
!   endif
!enddo
!if(allocated(mult)) deallocate(mult)


 call mpas_log_write('--- end subroutine ctrans_thompson:')
 call mpas_log_write(' ')

 end subroutine ctrans_thompson

!=================================================================================================================
 subroutine ctrans_thompson_up(its,ite,kts,kte,ierr,k22,kbcon,ktop,nthom,use_thom_scaven,c0,grav,fscav,p,p_cup, &
                       z_cup,thom,thoml2,nifa,nifal2,nwfa,nwfal2,zu,up_entr,up_detr,thom_up,nifa_up,nwfa_up,    &
                       thom_up_scav,nifa_up_scav,nwfa_up_scav,thom_up_scav_tot,nifa_up_scav_tot,nwfa_up_scav_tot)
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

 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte,nthom):: &
    thom,          &!
    thoml2          !


!--- inout arguments:
 real(kind=RKIND),intent(inout),dimension(its:ite):: &
    nifa_up_scav_tot, &!
    nwfa_up_scav_tot   !

 real(kind=RKIND),intent(inout),dimension(its:ite,nthom):: &
    thom_up_scav_tot   !

 real(kind=RKIND),intent(inout),dimension(its:ite,kts:kte):: &
    nifa_up,     &!
    nwfa_up,     &!
    nifa_up_scav,&!
    nwfa_up_scav  !

 real(kind=RKIND),intent(inout),dimension(its:ite,kts:kte,nthom):: &
    thom_up,     &!
    thom_up_scav  !

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

             !do we assume that nifa_up_scav_tot is at in-cloud levels or enviromnental levels?
             dp = 100._RKIND*(p_cup(i,k) - p_cup(i,k+1))/grav
             nifa_up_scav_tot(i) = nifa_up_scav_tot(i) + nifa_up_scav(i,k)*dp
             nwfa_up_scav_tot(i) = nwfa_up_scav_tot(i) + nwfa_up_scav(i,k)*dp

             nwfa_up(i,k) = nwfa_up(i,k) - nwfa_up_scav(i,k)
             nifa_up(i,k) = nifa_up(i,k) - nifa_up_scav(i,k)
          endif

!         thom_b = nwfa_up(i,k)
!         thom_c = nwfa_wetr*thom_b        ! nb of activated (liquid water+rain) water-friendly aerosols.
!         thom_g = thom_b - thom_c         ! nb of non-activated (gas) water-friendly aerosols.
!         thom_l = thom_c/(1.+c0*dz)       ! nb of activated (liquid) water-friendly aerosols.
!         thom_r = c0*dz*thom_c/(1.+c0*dz) ! nb of activated (rain) water-friendly aerosols.
!         nwfa_up(i,k)   = thom_g + thom_l
!         nwfa_scav(i,k) = thom_r
!         coef1 = 1./(1.+c0*dz)
!         coef2 = c0*dz/(1.+c0*dz)
       enddo

       !--- convective transport of individual aerosols in updrafts (use upstream method to define the
       !    in-cloud variables at the k2 levels):
       k = k2
       do n = 1, nthom
          thom_up(i,k,n) = thoml2(i,k-1,n)
       enddo
       do k = k2+1,kte
          dz = z_cup(i,k)-z_cup(i,k-1)
          denom = zu(i,k-1)-.5*up_detr(i,k-1)+up_entr(i,k-1)
          if(denom .lt. 1.e-16) exit

          !transport in updrafts:
          do n = 1,nthom
             thom_up(i,k,n) = (thom_up(i,k-1,n)*zu(i,k-1) - .5*up_detr(i,k-1)*thom_up(i,k-1,n) &
                            + up_entr(i,k-1)*thom(i,k-1,n))/denom
          enddo

          !wet scavenging:
          !formulation as in GOCART (using sourcecode in ConvPar_GF_GEOS5 provided by Saulo Freitas):
          if(use_thom_scaven == 1) then
             do n = 1,nthom
                thom_up_scav(i,k,n) = max(0.,thom_up(i,k,n)*(1.-exp(-fscav(n)*(dz/1000.))))

                !do we assume that nifa_up_scav_tot is at in-cloud levels or enviromnental levels?
                dp = 100._RKIND*(p_cup(i,k) - p_cup(i,k+1))/grav
                thom_up_scav_tot(i,n) = thom_up_scav_tot(i,n) + thom_up_scav(i,k,n)*dp

                thom_up(i,k,n) = thom_up(i,k,n) - thom_up_scav(i,k,n)
             enddo
          endif
       enddo
    endif
 enddo

 call mpas_log_write('--- end subroutine ctrans_thompson_up:')

 end subroutine ctrans_thompson_up

!=================================================================================================================
 subroutine ctrans_thompson_dn(its,ite,kts,kte,ierr,k22,kbcon,ktop,kdn,nthom,use_thom_evap,cumulus,grav,fscav, &
                      p,p_cup,z_cup,thom,thoml2,thom_up,nifa,nifal2,nifa_up,nwfa,nwfal2,nwfa_up,zd,dn_entr,    &
                      dn_detr,pwavo,pwevo,pwo,pwdo,thom_up_scav_tot,nifa_up_scav_tot,nwfa_up_scav_tot,thom_dn, &
                      thom_dn_scav,nifa_dn,nifa_dn_scav,nwfa_dn,nwfa_dn_scav)
!=================================================================================================================

!--- input arguments:
 character(len=*),intent(in):: &
    cumulus            !

 integer,intent(in):: its,ite,kts,kte

 integer,intent(in):: &
    nthom,            &! number of aerosols that describe "water-friendly" and "ice-friendly" aerosols.
    use_thom_evap      !                       

 integer,intent(in),dimension(its:ite):: &
    ierr,             &!
    k22,              &!
    kbcon,            &!
    ktop,             &!
    kdn                !

 real(kind=RKIND):: &
    grav               !

 real(kind=RKIND),intent(in),dimension(its:ite):: &
    pwavo,            &!total rain in each updraft.
    pwevo,            &!total evaporation in each downdraft.
    nifa_up_scav_tot, &!
    nwfa_up_scav_tot   !

 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte):: &
    p,                &!
    p_cup,            &!
    z_cup ,           &!
    pwo,              &!rain in each layer of updraft.
    pwdo               !evaporation in each layer of updraft.

 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte):: &
    dn_entr,          &!
    dn_detr,          &!
    zd                 !

 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte):: &
    nifa,             &!
    nifal2,           &!
    nifa_up,          &!
    nwfa,             &!
    nwfal2,           &!
    nwfa_up            !

 real(kind=RKIND),intent(in),dimension(nthom):: &
    fscav              !

 real(kind=RKIND),intent(in),dimension(its:ite,nthom):: &
    thom_up_scav_tot   !

 real(kind=RKIND),intent(in),dimension(its:ite,kts:kte,nthom):: &
    thom,             &!
    thoml2,           &!
    thom_up            !


!--- inout arguments:
 real(kind=RKIND),intent(inout),dimension(its:ite,kts:kte):: &
    nifa_dn,          &!
    nwfa_dn,          &!
    nifa_dn_scav,     &!
    nwfa_dn_scav       !

 real(kind=RKIND),intent(inout),dimension(its:ite,kts:kte,nthom):: &
    thom_dn,          &!
    thom_dn_scav       !


!--- local variables and arrays:
 integer:: i,k,k2,kb,kd,kt,n

 real(kind=RKIND):: frac_evap,pwdper
 real(kind=RKIND):: denom,dp,dz

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine ctrans_thompson_dn: nthom = $i',intArgs=(/nthom/))

 if(cumulus == 'shallow') return

 do i = its,ite
    if(ierr(i) /= 0) cycle
    k2 = k22(i)
    kb = kbcon(i)
    kd = kdn(i)
    kt = ktop(i)

    !--- fraction of the total rain that was evaporated:
    frac_evap = - pwevo(i)/(1.e-16+pwavo(i))

    !--- initializes downdrafts at kd level with a mixture of 90% environmental and 10% in-cloud mixing ratios,
    !    or initializes downdrafts using upstream approximation:
    k = kd
!   nifa_dn(i,k) = 0.9*nifa(i,k) + 0.1*nifa_up(i,k)
!   nwfa_dn(i,k) = 0.9*nwfa(i,k) + 0.1*nwfa_up(i,k)
!   do n = 1,nthom
!      thom_dn(i,k,n) = 0.9*thom(i,k,n) + 0.1*thom_up(i,k,n)
!   enddo
    nifa_dn(i,k) = nifa(i,k)
    nwfa_dn(i,k) = nwfa(i,k)
    do n = 1, nthom
       thom_dn(i,k,n) = thom(i,k,n)
    enddo

    pwdper = pwdo(i,k)/(1.e-16+pwevo(i))*frac_evap ! > 0
    if(use_thom_evap == 0) pwdper = 0._RKIND

    dp= 100.*(p_cup(i,k)-p_cup(i,k+1))
    nifa_dn_scav(i,k) = -pwdper*nifa_up_scav_tot(i)*grav/dp
    nifa_dn(i,k) = nifa_dn(i,k) - nifa_dn_scav(i,k)
    nwfa_dn_scav(i,k) = -pwdper*nwfa_up_scav_tot(i)*grav/dp
    nwfa_dn(i,k) = nwfa_dn(i,k) - nwfa_dn_scav(i,k)

    do n = 1,nthom
       thom_dn_scav(i,k,n) = -pwdper*thom_up_scav_tot(i,n)*grav/dp
       thom_dn(i,k,n) = thom_dn(i,k,n) - thom_dn_scav(i,k,n)
    enddo

    !--- calculates the downdraft mass terms:
    do k = kd-1,kts,-1
       dz = z_cup(i,k+1)-z_cup(i,k)
       denom = zd(i,k+1)-.5*dn_detr(i,k)+dn_entr(i,k)
       if(denom > 0._RKIND) then
          nifa_dn(i,k) = (nifa_dn(i,k+1)*zd(i,k+1) - .5*dn_detr(i,k)*nifa_dn(i,k+1) &
                       + dn_entr(i,k)*nifa(i,k))/denom

          nwfa_dn(i,k) = (nwfa_dn(i,k+1)*zd(i,k+1) - .5*dn_detr(i,k)*nwfa_dn(i,k+1) &
                       +  dn_entr(i,k)*nwfa(i,k))/denom

          do n = 1,nthom
             thom_dn(i,k,n) = (thom_dn(i,k+1,n)*zd(i,k+1) - .5*dn_detr(i,k)*thom_dn(i,k+1,n) &
                            +  dn_entr(i,k)*thom(i,k,n))/denom
          enddo
       else
          nifa_dn(i,k) = nifa(i,k+1)
          nwfa_dn(i,k) = nwfa(i,k+1)

          do n = 1,nthom
             thom_dn(i,k,n) = thom_dn(i,k+1,n)
          enddo
       endif

       !--- evaporation:
       pwdper = pwdo(i,k)/(1.e-16+pwevo(i))*frac_evap ! > 0
       if(use_thom_evap == 0) cycle

       dp = 100.*(p_cup(i,k)-p_cup(i,k+1))
       nifa_dn_scav(i,k) = -pwdper*nifa_up_scav_tot(i)*grav/dp
       nifa_dn(i,k) = nifa_dn(i,k) -nifa_dn_scav(i,k)
       nwfa_dn_scav(i,k) = -pwdper*nwfa_up_scav_tot(i)*grav/dp
       nwfa_dn(i,k) = nwfa_dn(i,k) -nwfa_dn_scav(i,k)

       do n = 1,nthom
          thom_dn_scav(i,k,n) = -pwdper*thom_up_scav_tot(i,n)*grav/dp
          thom_dn(i,k,n) = thom_dn(i,k,n) - thom_dn_scav(i,k,n)
       enddo
    enddo
 enddo

 call mpas_log_write('--- end subroutine ctrans_thompson_dn:')

 end subroutine ctrans_thompson_dn

!=================================================================================================================
 end module module_cu_gf_ctrans_thompson
!=================================================================================================================
