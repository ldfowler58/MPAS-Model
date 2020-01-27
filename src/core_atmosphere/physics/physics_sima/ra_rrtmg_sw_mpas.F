!=================================================================================================================
 module rrtmg_sw_mpas
 use mpas_kind_types,only: kind_phys => RKIND
 use mpas_log

 use module_ra_rrtmg_vinterp,only: vinterp_ozn
 use module_ra_rrtmg_lw,only: inirad,reicalc,relcalc,retab
!use parrrsw,only: nbndsw,ngptsw,naerec

 implicit none
 public:: rrtmg_sw_mpas_run,           &
          rrtmg_sw_mpas_timestep_init, &
          rrtmg_sw_mpas_timestep_final


 contains


!   call rrtmg_sw_run &
!              (dorrsw    , ncol     , nlay      , icld      , play      , plev      , &
!               tlay      , tlev     , tsfc      , h2ovmr    , o3vmr     , co2vmr    , &
!               ch4vmr    , n2ovmr   , o2vmr     , asdir     , asdif     , aldir     , &
!               aldif     , coszen   , adjes     , dyofyr    , scon      , inflgsw   , &
!               iceflgsw  , liqflgsw , cldfmcl   , taucmcl   , ssacmcl   , asmcmcl   , &
!               fsfcmcl   , ciwpmcl  , clwpmcl   , cswpmcl   , reicmcl   , relqmcl   , &
!               resnmcl   , tauaer   , ssaaer    , asmaer    , ecaer     , swuflx    , &
!               swdflx    , swhr     , swuflxc   , swdflxc   , swhrc     , swuflxcln , &
!               swdflxcln , aer_opt  , sibvisdir , sibvisdif , sibnirdir , sibnirdif , &
!               swdkdir   , swdkdif  , swdkdirc  , calc_clean_atm_diag   , errmsg    , &
!               errflg)


!   call rrtmg_sw_run &
!               cldfmcl   , taucmcl  , ssacmcl   , asmcmcl   ,                         &
!               fsfcmcl   , ciwpmcl  , clwpmcl   , cswpmcl   , reicmcl   , relqmcl   , &
!               resnmcl   , tauaer   , ssaaer    , asmaer    , ecaer     , swuflx    , &
!               swdflx    , swhr     , swuflxc   , swdflxc   , swhrc     , swuflxcln , &
!               swdflxcln , aer_opt  , sibvisdir , sibvisdif , sibnirdir , sibnirdif , &
!               swdkdir   , swdkdif  , swdkdirc  , calc_clean_atm_diag   , errmsg    , &
!               errflg)

!=================================================================================================================
 subroutine rrtmg_sw_mpas_run(errmsg,errflg)
!=================================================================================================================

!--- output arguments:
 character(len=*),intent(out):: errmsg
 integer,intent(out)::errflg

 errmsg = 'rrtmg_sw_mpas_run OK'
 errflg = 0

 end subroutine rrtmg_sw_mpas_run

!=================================================================================================================
 subroutine rrtmg_sw_mpas_timestep_init &
                  (ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte, &
                   f_qc     , f_qr       , f_qi       , f_qs       , icloud     , &
                   has_reqc , has_reqi   , has_reqs   , cldovrlp   , gravity    , &
                   gmt      , declin     , degrad     , radt       , xtime      , &
                   xlat     , xlon       , tsk        , albedo     , xland      , &
                   xice     , snow       , p3d        , t3d        , qv3d       , &
                   qc3d     , qr3d       , qi3d       , qs3d       , cldfra3d   , &
                   re_cloud , re_ice     , re_snow    , dz8w       , p8w        , &
                   t8w      , alswvisdir , alswvisdif , alswnirdir , alswnirdif , &
                   o3input  , noznlevels , pin        , o3clim     , ncol       , &
                   nlay     , icld       , inflgsw    , iceflgsw   , liqflgsw   , &
                   dyofyr   , adjes      , dorrsw     , coszen     , tsfc       , &
                   asdir    , asdif      , aldir      , aldif      , play       , &
                   tlay     , pdel       , h2ovmr     , o3vmr      , co2vmr     , &
                   o2vmr    , ch4vmr     , n2ovmr     , plev       , tlev       , &
                   reliq    , reice      , resnow     , cliqwp     , cicewp     , &
                   csnowp   , hgt        , cldfrac)
!=================================================================================================================

!input arguments:
 logical,intent(in),optional:: &
    f_qc,       &!
    f_qr,       &!
    f_qi,       &!
    f_qs         !

 integer,intent(in):: ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte

 integer,intent(in):: &
    icloud,     &!
    has_reqc,   &!
    has_reqi,   &!
    has_reqs,   &!
    cldovrlp,   &!
    noznlevels, &!
    o3input      !

 real(kind=kind_phys),intent(in):: &
    declin,     &!
    degrad,     &!
    gmt,        &!
    gravity,    &!
    radt,       &!
    xtime        !

 real(kind=kind_phys),intent(in),dimension(ims:ime,jms:jme):: &
    xlat,       &!
    xlon,       &!
    xland,      &!
    xice,       &!
    snow,       &!
    albedo,     &! surface albedo (-).
    tsk          ! surface temperature (K).

 real(kind=kind_phys),intent(in),dimension(ims:ime,kms:kme,jms:jme):: &
    t3d,        &!
    p3d,        &!
    qv3d,       &!
    re_cloud,   &!
    re_ice,     &!
    re_snow      !

 real(kind=kind_phys),intent(in),dimension(ims:ime,jms:jme),optional:: &
   alswvisdir,  &!
   alswvisdif,  &!
   alswnirdir,  &!
   alswnirdif    !

 real(kind=kind_phys),intent(in),dimension(ims:ime,kms:kme,jms:jme):: &
    dz8w,       &!
    p8w,        &!
    t8w          !

 real(kind=kind_phys),intent(in),dimension(ims:ime,kms:kme,jms:jme),optional:: &
    cldfra3d,   &!
    qc3d,       &!
    qr3d,       &!
    qi3d,       &!
    qs3d         !

 real(kind=kind_phys),intent(in),dimension(1:noznlevels),optional:: &
    pin          !

 real(kind=kind_phys),intent(in),dimension(ims:ime,1:noznlevels,jms:jme),optional:: &
    o3clim       !

!output arguments:
 logical,intent(out),dimension(its:ite):: &
    dorrsw     !

 integer,intent(out):: &
    ncol,     &!
    nlay,     &!
    dyofyr     !

 integer,intent(out):: &
    icld,     &!
    inflgsw,  &!
    iceflgsw, &!
    liqflgsw   !

 real(kind=kind_phys),intent(out):: &
    adjes

 real(kind=kind_phys),intent(out),dimension(its:ite):: &
    coszen,   &!
    tsfc,     &!
    asdir,    &!
    asdif,    &!
    aldir,    &!
    aldif      !

 real(kind=kind_phys),intent(out),dimension(its:ite,kts:kte):: &
    pdel

 real(kind=kind_phys),intent(out),dimension(its:ite,kts:kte):: &
    cicewp,   &!
    cliqwp,   &!
    csnowp,   &!
    reliq,    &!
    reice,    &!
    resnow     !

 real(kind=kind_phys),intent(out),dimension(its:ite,kts:kte+1):: &
    play,     &!
    tlay,     &!
    h2ovmr,   &!

    o3vmr,    &!
    co2vmr,   &!
    o2vmr,    &!
    ch4vmr,   &!
    n2ovmr     !

 real(kind=kind_phys),intent(out),dimension(its:ite,kts:kte+1):: &
    hgt,      &!
    cldfrac    !

 real(kind=kind_phys),intent(out),dimension(its:ite,kts:kte+2):: &
    plev,     &!
    tlev

!local variables and arrays:
 integer:: i,j,k,n
 integer:: pcols,pver
 integer:: idx_rei

 real(kind=kind_phys):: coszrs,hrang,tloctm,xt24,xxlat
 real(kind=kind_phys):: dz,dzsum
 real(kind=kind_phys):: corr
 real(kind=kind_phys):: gliqwp,gicewp,gsnowp,snow_mass_factor
 real(kind=kind_phys),dimension(1):: landfrac,landm,snowh,icefrac
 real(kind=kind_phys),dimension(kts:kte):: qc1d,qr1d,qi1d,qs1d
 real(kind=kind_phys),dimension(kts:kte+1):: o3mmr,p1dlay
 real(kind=kind_phys),dimension(kts:kte+2):: p1dlev
 real(kind=kind_phys),dimension(1,kts:kte):: t1d,rei1d,rel1d
 real(kind=kind_phys),dimension(1:noznlevels):: o3clim1d

!real(kind=kind_phys),dimension(its:ite,kts:kte+1):: hgt,cldfrac

!--- set trace gas volume mixing ratios, 2005 values, IPCC (2007):
!carbon dioxide (379 ppmv)
 real(kind=kind_phys):: co2
 data co2 / 379.e-6 /
!methane (1774 ppbv)
 real(kind=kind_phys):: ch4
 data ch4 / 1774.e-9 /
!nitrous oxide (319 ppbv)
 real(kind=kind_phys):: n2o
 data n2o / 319.e-9 /
!--- set oxygen volume mixing ratio (for o2mmr=0.23143):
 real(kind=kind_phys):: o2
 data o2 / 0.209488 /

!--- molecular weights and ratios for converting mmr to vmr units
!real(kind=kind_phys):: amd       ! Effective molecular weight of dry air (g/mol)
!real(kind=kind_phys):: amw       ! Molecular weight of water vapor (g/mol)
!real(kind=kind_phys):: amo       ! Molecular weight of ozone (g/mol)
!real(kind=kind_phys):: amo2      ! Molecular weight of oxygen (g/mol)
!Atomic weights for conversion from mass to volume mixing ratios
!data amd   /  28.9660   /
!data amw   /  18.0160   /
!data amo   /  47.9998   /
!data amo2  /  31.9999   /

 real(kind=kind_phys):: amdw     ! Molecular weight of dry air / water vapor
 real(kind=kind_phys):: amdo     ! Molecular weight of dry air / ozone
 real(kind=kind_phys):: amdo2    ! Molecular weight of dry air / oxygen
 data amdw /  1.607793 /
 data amdo /  0.603461 /
 data amdo2 / 0.905190 /

!-----------------------------------------------------------------------------------------------------------------

!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine rrtmg_sw_mpas_timestep_init:')

!--- number of columns and number of layers between the surface and the top-of-the-atmosphere (i.e. number of
!    MPAS layers plus one):
 ncol = ite-its+1
 nlay = (kte-kts+1) + 1

!--- initialize option for the calculation of the cloud optical properties:
 icld     = cldovrlp
 inflgsw  = 2
 iceflgsw = 3
 liqflgsw = 1

!--- set solar constant:
 dyofyr = 0
 adjes  = 1.0

!--- calculation of the local input soundings for daytime columns:
 j_loop: do j = jts,jte
    i_loop: do i = its,ite

       !calculation of the solar zenith angle:
       xt24 = mod(xtime+radt*0.5,1440.)
       tloctm = gmt + xt24/60. + xlon(i,j)/15.
       hrang = 15. * (tloctm-12.) * degrad
       xxlat = xlat(i,j) * degrad
       coszrs = sin(xxlat) * sin(declin) + cos(xxlat) * cos(declin) * cos(hrang)
       coszen(i) = coszrs

       dorrsw(i) = .true.
       if(coszen(i) .le. 0._kind_phys) dorrsw(i) = .false.

       if(dorrsw(i)) then !begin if statement checking for daytime columns.

          !initialization of the surface albedo:
          if(present(alswvisdir) .and. present(alswvisdir) .and. &
             present(alswnirdir) .and. present(alswnirdif)) then
             asdir(i) = alswvisdir(i,j)
             asdif(i) = alswvisdif(i,j)
             aldir(i) = alswnirdir(i,j)
             aldif(i) = alswnirdif(i,j)
          else
             asdir(i) = albedo(i,j)
             asdif(i) = albedo(i,j)
             aldir(i) = albedo(i,j)
             aldif(i) = albedo(i,j)
          endif

          !initialization of the local input soundings:
          k = kts
          plev(i,kts) = p8w(i,k,j) / 100.
          tlev(i,kts) = t8w(i,k,j)
          tsfc(i)     = tsk(i,j)
          do k = kts,kte
             plev(i,k+1) = p8w(i,k+1,j) / 100.
             tlev(i,k+1) = t8w(i,k+1,j)
             play(i,k)   = p3d(i,k,j) / 100.
             pdel(i,k)   = plev(i,k) - plev(i,k+1)
             tlay(i,k)   = t3d(i,k,j)
             h2ovmr(i,k) = amax1(qv3d(i,k,j),1.e-12)
             h2ovmr(i,k) = h2ovmr(i,k)*amdw
             co2vmr(i,k) = co2
             o2vmr(i,k)  = o2
             ch4vmr(i,k) = ch4
             n2ovmr(i,k) = n2o
          enddo
          plev(i,kte+2)   = 1.0e-5
          tlev(i,kte+2)   = tlev(i,kte+1)
          play(i,kte+1)   = 0.5 * plev(i,kte+1)
          tlay(i,kte+1)   = tlev(i,kte+1)
          h2ovmr(i,kte+1) = h2ovmr(i,kte)
          co2vmr(i,kte+1) = co2vmr(i,kte)
          o2vmr(i,kte+1)  = o2vmr(i,kte)
          ch4vmr(i,kte+1) = ch4vmr(i,kte)
          n2ovmr(i,kte+1) = n2ovmr(i,kte)

          !initialization of the ozone volume mixing ratio:
          pcols = 1
          do k = kts,kte+1
             p1dlay(k) = play(i,k)
             p1dlev(k) = plev(i,k)
          enddo
          p1dlev(kte+2) = plev(i,kte+2)
          call inirad(o3mmr,p1dlev,kts,kte)

          if(o3input .eq. 2) then
             do k = 1, noznlevels
                o3clim1d(k) = o3clim(i,k,j)
             enddo
             call vinterp_ozn(1,pcols,pcols,kte+1,p1dlay,pin,noznlevels,o3clim1d,o3mmr)
             do k = kts,kte+1
                o3vmr(i,k) = o3mmr(k)
             enddo
          else
             do k = kts,kte+1
                o3vmr(i,k) = o3mmr(k)*amdo
             enddo
          endif
       endif !end if statement checking for daytime columns.

    enddo i_loop
 enddo j_loop

!--- calculation of the cloud optical properties for daytime columns:
 do i = its,ite
 do k = kts,kte
    reliq(i,k)  = 0._kind_phys
    reice(i,k)  = 0._kind_phys
    resnow(i,k) = 0._kind_phys
    cliqwp(i,k) = 0._kind_phys
    cicewp(i,k) = 0._kind_phys
    csnowp(i,k) = 0._kind_phys
 enddo
 enddo

 if(icloud .gt. 0) then !begin cloud effect to shortwave radiation.
    j_loop2: do j = jts,jte
       i_loop2: do i = its,ite

          if(dorrsw(i)) then !begin if statement checking for daytime columns.

             dzsum = 0.0
             do k = kts, kte
                dz = dz8w(i,k,j)
                hgt(i,k) = dzsum + 0.5*dz
                dzsum = dzsum + dz
             enddo
             hgt(i,kte+1) = 0.

             do k = kts,kte
                qc1d(k) = 0.
                qr1d(k) = 0.
                qi1d(k) = 0.
                qs1d(k) = 0.
             enddo
             if(present(f_qc) .and. present(qc3d)) then
                do k = kts,kte
                   qc1d(k) = amax1(qc3d(i,k,j),0.)
                enddo
             endif
             if(present(f_qr) .and. present(qr3d)) then
                do k = kts,kte
                   qr1d(k) = amax1(qr3d(i,k,j),0.)
                enddo
             endif
             if(present(f_qi) .and. present(qi3d)) then
                do k = kts,kte
                   qi1d(k) = amax1(qi3d(i,k,j),0.)
                enddo
             endif
             if(present(f_qs) .and. present(qs3d)) then
                do k = kts,kte
                   qs1d(k) = amax1(qs3d(i,k,j),0.)
                enddo
             endif
             if(present(f_qc).and.present(qc3d) .and. present(f_qi).and.present(qi3d)  .and. &
                present(f_qs).and.present(qs3d) .and. present(cldfra3d)) then
                do k = kts,kte
                   cldfrac(i,k) = cldfra3d(i,k,j)
                enddo
             endif

             if(inflgsw .gt. 0) then
                !--- calculation of the effective radii for cloud liquid water, cloud ice, and snow:
                if(has_reqc .ne. 0) then
                   !fill the effective radius for cloud water with that calculated in the Thompson cloud
                   !microphysics parameterization:
                   inflgsw = 3
                   do k = kts,kte
                      reliq(i,k) = max(2.5,re_cloud(i,k,j)*1.e6)
                      if(reliq(i,k).le.2.5 .and. cldfrac(i,k).gt.0. .and. &
                         (xland(i,j)-1.5).gt.0.) then        !--- ocean.
                         reliq(i,k) = 10.5
                      elseif(reliq(i,k).le.2.5 .and. cldfrac(i,k).gt.0. .and. &
                         (xland(i,j)-1.5).lt.0.) then        !--- land.
                         reliq(i,k) = 7.5
                      endif
                   enddo
                else
                   !calculate the effective radius for cloud water in layers below the model top. we added the
                   !dimensions pcols, pver so that the calls to subroutines relcalc and reicalc remain the same
                   !as in WRF although the two dimensions are duplicate and not needed inside the subroutines:
                   pcols = 1
                   pver  = nlay-1
                   landfrac(pcols) = 2.-xland(i,j)
                   landm(pcols)    = landfrac(pcols)
                   snowh(pcols)    = 0.001*snow(i,j)
                   icefrac(pcols)  = xice(i,j)
                   do k = kts,kte
                      t1d(pcols,k) = tlay(i,k)
                   enddo
                   call relcalc(pcols,pcols,pver,t1d,landfrac,landm,icefrac,rel1d,snowh)
                   do k = kts,kte
                      reliq(i,k) = rel1d(pcols,k)
                   enddo
                endif
                if(has_reqi .ne. 0) then
                   !fill the effective radius for cloud ice with that calculated in the Thompson cloud
                   !microphysics parameterization:
                   inflgsw  = 4
                   iceflgsw = 4
                   do k = kts,kte
                      reice(i,k) = max(5.,re_ice(i,k,j)*1.e6)
                      if(reice(i,k).le.5. .and. cldfrac(i,k).gt.0.) then
                         idx_rei = int(tlay(i,k)-179.)
                         idx_rei = min(max(idx_rei,1),75)
                         corr = tlay(i,k) - int(tlay(i,k))
                         reice(i,k) = retab(idx_rei)*(1.-corr) + retab(idx_rei+1)*corr
                         reice(i,k) = max(reice(i,k),5.0)
                      endif
                   enddo
                else
                   !calculate the effective radius for cloud ice in layers below the model top, following
                   !Kristjansson and Mitchell. limit upper bound of reice for Fu ice parameterization and
                   !convert effective radius to generalized effective size (*1.0315; Fu, 1996):
                   pcols = 1
                   pver  = nlay-1
                   do k = kts,kte
                      t1d(pcols,k) = tlay(i,k)
                   enddo
                   call reicalc(pcols,pcols,pver,t1d,rei1d)
                   do k = kts,kte
                      reice(i,k) = rei1d(pcols,k)
                   enddo
                   if(iceflgsw .eq. 3) then
                      do k = kts, kte
                         reice(i,k) = reice(i,k) * 1.0315
                         reice(i,k) = min(140.0,reice(i,k))
                      enddo
                   endif
                endif
                if(has_reqs .ne. 0) then
                   !fill the effective radius for snow with that calculated in the Thompson cloud
                   !microphysics parameterization:
                   inflgsw  = 5
                   iceflgsw = 5
                   do k = kts, kte
                      resnow(i,k) = max(10.,re_snow(i,k,j)*1.e6)
                   enddo
                else
                   do k = kts, kte
                      resnow(i,k) = 10.
                   enddo
                endif
             endif

             !--- calculation of the cloud liquid water, cloud ice, and snow paths (units are g/m2):
             !cloud liquid water path:
             do k = kts, kte
                gliqwp = qc1d(k) * pdel(i,k)*100.0 / gravity * 1000.0              !grid box liquid water path.
                cliqwp(i,k) = gliqwp / max(0.01,cldfrac(i,k))                      !in-cloud liquid water path.
             enddo
             !cloud ice path:
             if(iceflgsw .eq. 3) then
                do k = kts, kte
                   gicewp = (qi1d(k)+qs1d(k)) * pdel(i,k)*100.0 / gravity * 1000.0 !grid box ice water path.
                   cicewp(i,k) = gicewp / max(0.01,cldfrac(i,k))                   !in-cloud ice water path.
                enddo
             elseif(iceflgsw .ge. 4) then
                do k = kts, kte
                   gicewp = qi1d(k) * pdel(i,k)*100.0 / gravity * 1000.0           !grid box ice water path.
                   cicewp(i,k) = gicewp / max(0.01,cldfrac(i,k))                   !in-cloud ice water path.
                enddo
             endif
             !snow path:
             if(iceflgsw.eq.5) then
                do k = kts, kte
                   snow_mass_factor = 1.0
                   if(resnow(i,k) .gt. 130.)then
                      snow_mass_factor = (130.0/resnow(i,k))*(130.0/resnow(i,k))
                      resnow(i,k)   = 130.0
                   endif
                   gsnowp = qs1d(k) * pdel(i,k)*100.0 / gravity * 1000.0          !grid box snow path.
                   csnowp(i,k) = snow_mass_factor * gsnowp / max(0.01,cldfrac(i,k))
                enddo
             endif

          endif !end if statement checking for daytime columns.

       enddo i_loop2
    enddo j_loop2
 endif

 call mpas_log_write('--- end subroutine rrtmg_sw_mpas_timestep_init:')
!call mpas_log_write(' ')

 end subroutine rrtmg_sw_mpas_timestep_init

!=================================================================================================================
 subroutine rrtmg_sw_mpas_timestep_final &
                  (dorrsw,pi3d,alswvisdir,alswvisdif,alswnirdir,alswnirdif,coszen,                 &
                   coszr,gsw,swhr,swhrc,swuflx,swuflxc,swdflx,swdflxc,swuflxcln,swdflxcln,         &
                   sibvisdir,sibvisdif,sibnirdir,sibnirdif,swdkdir,swdkdirc,swdkdif,               &
                   rthratensw,swcf,swupt,swuptc,swdnt,swdntc,swupb,swupbc,swdnb,swdnbc,            &
                   swuptcln,swdntcln,swupbcln,swdnbcln,swvisdir,swvisdif,swnirdir,swnirdif,        &
                   swddir,swddni,swddif,swddirc,swddnic,swdownc,swupflx,swupflxc,swdnflx,swdnflxc, &
                   ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)
!=================================================================================================================

!input arguments:
 logical,intent(in),dimension(its:ite):: dorrsw

 integer,intent(in):: ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte

 real(kind=kind_phys),intent(in),dimension(its:ite):: &
    coszen       !

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte+1):: &
    swhr,       &! all-sky shortwave heating rate (K/day).
    swhrc        ! clear-sky shortwave heating rate (K/day).

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte+2):: &
    swuflx,     &! all-sky upward shortwave radiation flux (W/m2).
    swuflxc,    &! clear-sky upward shortwave radiation flux (W/m2)
    swdflx,     &! all-sky downward shortwave radiation flux (W/m2).
    swdflxc      ! clear-sky downward shortwave radiation flux (W/m2).

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte+2):: &
    swuflxcln,  &!
    swdflxcln    !

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte+2):: &
    sibvisdir,  &!
    sibvisdif,  &!
    sibnirdir,  &!
    sibnirdif    !

 real(kind=kind_phys),dimension(its:ite,kts:kte+2):: &
    swdkdir,    &!
    swdkdirc,   &!
    swdkdif      !

 real(kind=kind_phys),intent(in),dimension(ims:ime,kms:kme,jms:jme):: &
    pi3d         ! exner function (-).

 real(kind=kind_phys),intent(in),dimension(ims:ime,jms:jme),optional:: &
    alswvisdir, &!
    alswvisdif, &!
    alswnirdir, &!
    alswnirdif   !

!output arguments:
 real(kind=kind_phys),intent(out),dimension(ims:ime,jms:jme):: &
    coszr,      &! cosine of the solar zenith angle (-).
    gsw          ! all-sky net surface radiation (W/m2).

 real(kind=kind_phys),intent(out),dimension(ims:ime,jms:jme),optional:: &
    swcf,       &! all-sky shortwave cloud radiative forcing (W/m2).
    swupt,      &! all-sky shortwave top-of-the-atmosphere upward radiation (W/m2).
    swuptc,     &! clear-sky shortwave top-of-the-atmosphere upward radiation (W/m2).
    swdnt,      &! all-sky shortwave top-of-the-atmosphere downward radiation (W/m2).
    swdntc,     &! clear-sky shortwave top-of-the-atmosphere downward radiation (W/m2).
    swupb,      &! all-sky shortwave surface upward radiation (W/m2).
    swupbc,     &! clear-sky shortwave surface upward radiation (W/m2).
    swdnb,      &! all-sky shortwave surface downward radiation (W/m2).
    swdnbc       ! clear-sky shortwave surface downward radiation (W/m2).

 real(kind=kind_phys),intent(out),dimension(ims:ime,jms:jme),optional:: &
    swuptcln,   &!
    swdntcln,   &!
    swupbcln,   &!
    swdnbcln     !

 real(kind=kind_phys),intent(out),dimension(ims:ime,jms:jme),optional:: &
    swvisdir,   &!
    swvisdif,   &!
    swnirdir,   &!
    swnirdif     !

 real(kind=kind_phys),intent(out),dimension(ims:ime,jms:jme),optional:: &
      swddir,   &!
      swddni,   &!
      swddif,   &!
      swddirc,  &!
      swddnic,  &!
      swdownc    !

 real(kind=kind_phys),intent(out),dimension(ims:ime,kms:kme,jms:jme):: &
    rthratensw   !

 real(kind=kind_phys),intent(out),dimension(ims:ime,kms:kme+1,jms:jme ),optional:: &
    swupflx,    &!
    swupflxc,   &!
    swdnflx,    &!
    swdnflxc     !

!local variables and arrays:
 integer:: i,k,j
 integer:: calc_clean_atm_diag
 integer:: sf_surface_physics

!-----------------------------------------------------------------------------------------------------------------

!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine rrtmg_sw_mpas_timestep_final:')

!--- initialization of calc_clean_atm_diag:
 calc_clean_atm_diag = 0
 if(present(swuptcln) .and. present(swdntcln) .and. present(swupbcln) .and. present(swdnbcln)) &
    calc_clean_atm_diag = 1

!--- initialization sf_surface_physics:
 sf_surface_physics = 0
 if(present(alswvisdir) .and. present(alswvisdif) .and. present(alswnirdir) .and. present(alswnirdif)) &
    sf_surface_physics = 8

!--- latitude loop:
 j_loop: do j = jts,jte

    !--- longitude loop:
    i_loop: do i = its,ite
       coszr(i,j) = coszen(i)
       if(dorrsw(i)) then 
          gsw(i,j) = swdflx(i,1) - swuflx(i,1)
          do k = kts,kte
             rthratensw(i,k,j) = swhr(i,k)/86400.
             rthratensw(i,k,j) = rthratensw(i,k,j)/pi3d(i,k,j)
          enddo

          !optional arguments:
          if(present(swupt)) then
             swcf(i,j) = (swdflx(i,kte+2)-swuflx(i,kte+2)) - (swdflxc(i,kte+2)-swuflxc(i,kte+2))
             swupt(i,j)  = swuflx(i,kte+2)
             swuptc(i,j) = swuflxc(i,kte+2)
             swdnt(i,j)  = swdflx(i,kte+2)
             swdntc(i,j) = swdflxc(i,kte+2)
             swupb(i,j)  = swuflx(i,1)
             swupbc(i,j) = swuflxc(i,1)
             swdnb(i,j)  = swdflx(i,1)
             swdnbc(i,j) = swdflxc(i,1)
          endif
          if(present(swupflx)) then
             do k = kts, kte+1
                swupflx(i,k,j)  = swuflx(i,k)
                swupflxc(i,k,j) = swuflxc(i,k)
                swdnflx(i,k,j)  = swdflx(i,k)
                swdnflxc(i,k,j) = swdflxc(i,k)
             enddo
          endif
          if(calc_clean_atm_diag .gt. 0)then
             swuptcln(i,j) = swuflxcln(i,kte+2)
             swdntcln(i,j) = swdflxcln(i,kte+2)
             swupbcln(i,j) = swuflxcln(i,1)
             swdnbcln(i,j) = swdflxcln(i,1)
          endif
          if(sf_surface_physics .eq. 8) then
             swvisdir(i,j) = sibvisdir(i,1)
             swvisdif(i,j) = sibvisdif(i,1)
             swnirdir(i,j) = sibnirdir(i,1)
             swnirdif(i,j) = sibnirdif(i,1)
          endif
          if(present(swddir) .and. present(swddif) .and. present(swddni) .and. &
             present(swddirc) .and. present(swddnic) .and. present(swdownc)) then
             swddir(i,j)   = swdkdir(i,1)
             swddni(i,j)   = swddir(i,j) / coszen(i)
             swddif(i,j)   = swdkdif(i,1)
             swdownc(i, j) = swdflxc(i,1)
             swddirc(i,j)  = swdkdirc(i,1)
             swddnic(i,j)  = swddirc(i,j) / coszen(i)
          endif
 
       else

          gsw(i,j) = 0._kind_phys
          do k = kts,kte
             rthratensw(i,k,j) = 0._kind_phys
          enddo

          !optional arguments:
          if(present(swupt)) then
             swcf(i,j)   = 0._kind_phys
             swupt(i,j)  = 0._kind_phys
             swuptc(i,j) = 0._kind_phys
             swdnt(i,j)  = 0._kind_phys
             swdntc(i,j) = 0._kind_phys
             swupb(i,j)  = 0._kind_phys
             swupbc(i,j) = 0._kind_phys
             swdnb(i,j)  = 0._kind_phys
             swdnbc(i,j) = 0._kind_phys
          endif
          if(present(swupflx)) then
             do k = kts, kte+1
                swupflx(i,k,j)  = 0._kind_phys
                swupflxc(i,k,j) = 0._kind_phys
                swdnflx(i,k,j)  = 0._kind_phys
                swdnflxc(i,k,j) = 0._kind_phys
             enddo
          endif
          if(calc_clean_atm_diag .gt. 0) then
             swuptcln(i,j) = 0._kind_phys
             swdntcln(i,j) = 0._kind_phys
             swupbcln(i,j) = 0._kind_phys
             swdnbcln(i,j) = 0._kind_phys
          endif
          if(sf_surface_physics .eq. 8) then
             swvisdir(i,j) = 0._kind_phys
             swvisdif(i,j) = 0._kind_phys
             swnirdir(i,j) = 0._kind_phys
             swnirdif(i,j) = 0._kind_phys
          endif
          if(present(swddir) .and. present(swddif) .and. present(swddni) .and. &
             present(swddirc) .and. present(swddnic) .and. present(swdownc)) then
             swddir(i,j)   = 0._kind_phys
             swddif(i,j)   = 0._kind_phys
             swdownc(i, j) = 0._kind_phys
             swddirc(i,j)  = 0._kind_phys
             swddnic(i,j)  = 0._kind_phys
          endif

       endif !endif daytime grid columns.

    enddo i_loop !end longitude loop.

 enddo j_loop !end latitude loop.

 call mpas_log_write('--- end subroutine rrtmg_sw_mpas_timestep_final:')
!call mpas_log_write(' ')

 end subroutine rrtmg_sw_mpas_timestep_final

!=================================================================================================================
 end module rrtmg_sw_mpas
!=================================================================================================================



