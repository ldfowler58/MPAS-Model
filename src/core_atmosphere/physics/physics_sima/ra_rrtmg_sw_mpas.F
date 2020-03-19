!=================================================================================================================
 module rrtmg_sw_mpas
 use mpas_kind_types,only: kind_phys => RKIND
 use mpas_log

 use ra_rrtmg_shared,only    : inirad,reicalc,relcalc,retab
 use ra_rrtmg_vinterp,only   : vinterp_ozn
 use ra_sw_radconst,only     : radconst
 use mcica_subcol_gen_sw,only: mcica_subcol_sw

 implicit none
 private
 public:: rrtmg_sw_mpas_run,           &
          rrtmg_sw_mpas_timestep_init, &
          rrtmg_sw_mpas_timestep_final


 contains


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
                  (its        , ite        , kts        , kte        , kme        , &
                   nbndsw     , ngptsw     , naerec     , f_qc       , f_qr       , &
                   f_qi       , f_qs       , icloud     , has_reqc   , has_reqi   , &
                   has_reqs   , icld       , gravity    , gmt        , julian     , &
                   declin     , degrad     , dpd        , solcon_0   , radt       , &
                   xtime      , xlat       , xlon       ,                           &
                   albedo     , xland      , xice       , snow       , p2d        , &
                   t2d        , qv2d       , qc2d       , qr2d       , qi2d       , &
                   qs2d       , cldfra2d   , re_cloud   , re_ice     , re_snow    , &
                   dz8w2d     , p2di       , t2di       , alswvisdir , alswvisdif , &
                   alswnirdir , alswnirdif , tsfc       , o3input    , noznlevels , &
                   pin        , o3clim     , ncol       , nlay       , inflgsw    , &
                   iceflgsw   , liqflgsw   , aer_opt    , dyofyr     , adjes      , &
                   dorrsw     , coszen     , asdir      , asdif      , aldir      , &
                   aldif      , play       , tlay       , h2ovmr     , o3vmr      , &
                   co2vmr     , o2vmr      , ch4vmr     , n2ovmr     , plev       , &
                   tlev       , cldfmcl    , clwpmcl    , ciwpmcl    , cswpmcl    , &
                   relqmcl    , reicmcl    , resnmcl    , taucmcl    , ssacmcl    , &
                   asmcmcl    , fsfcmcl    , tauaer     , ssaaer     , asmaer     , &
                   ecaer      , errmsg     , errflg)
!=================================================================================================================

!input arguments:
 logical,intent(in),optional:: &
    f_qc,       &! if true, the physics package includes the cloud liquid water mixing ratio.
    f_qr,       &! if true, the physics package includes the rain mixing ratio.
    f_qi,       &! if true, the physics package includes the cloud ice mixing ratio.
    f_qs         ! if true, the physics package includes the snow mixing ratio.

 integer,intent(in):: &
    its,ite,    &! first and last index for horizontal dimension.
    kts,kte,    &! first and last index for vertical dimension.
    kme          ! defined as kte+1.

 integer,intent(in):: &
    nbndsw,     &! number of shortwave spectral bands.
    ngptsw,     &! number of reduced g-intervals for rrtmg_sw.
    naerec       ! number of aerosols for the option iaer=6 (ecmwf aerosols option)

 integer,intent(in):: &
    icloud,     &! option that defines the cloud effect on the SW optical depth.
    has_reqc,   &! = 1 if effective radius of cloud liquid water computed from cloud microphysics, = 0 otherwise.
    has_reqi,   &! = 1 if effective radius of cloud ice computed from cloud microphysics, = 0 otherwise.
    has_reqs,   &! = 1 if effective radius of snow computed from cloud microphysics, = 0 otherwise.
    icld,       &! option that defines the type of cloud overlap.
    noznlevels, &! number of levels of climatological ozone data.
    o3input      ! ozone input option for shortwave radiation.

 real(kind=kind_phys),intent(in):: &
    declin,     &! solar declination (radian).
    dpd,        &! ratio between two fixed julian days (-).
    degrad,     &! conversion factor from degree to radian (-).
    gmt,        &! GMT time.
    gravity,    &! acceleration due to gravity (m/s2).
    julian,     &! current julian day (-).
    radt,       &! shortwave radiation timestep (s).
    solcon_0,   &! initial solar constant (W/m2).
    xtime        ! MPAS run time.

 real(kind=kind_phys),intent(in),dimension(its:ite):: &
    xlat,       &! latitude (radian).
    xlon,       &! longitude (radian).
    xland,      &! land-ocean mask (1=land including sea-ice ; 2=ocean)
    xice,       &! fractional seaice amount (-).
    snow,       &! snow water equivalent (kg/m2).
    albedo,     &! surface albedo (-).
    tsfc         ! surface temperature (K).

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte):: &
    t2d,        &! temperature at middle of layers (K).
    p2d,        &! pressure at middle of layers (Pa).
    qv2d,       &! water vapor mixing ratio at middle of layers (kg/kg).
    dz8w2d,     &! thickness of layers (m).
    re_cloud,   &! effective radius of cloud liquid droplets at middle of layers (m).
    re_ice,     &! effective radius of cloud ice crystals at middle of layers (m). 
    re_snow      ! effective radius of snow crystals at middle of layers (m).

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kme):: &
    p2di,       &! pressure at interface of layers (Pa).
    t2di         ! temperature at interface of layers (K).

 real(kind=kind_phys),intent(in),dimension(its:ite),optional:: &
    alswvisdir, &! visible direct component of surface albedo when using SIB LSM (-).
    alswvisdif, &! visible diffuse component of surface albedo when using SIB LSM (-).
    alswnirdir, &! near-IR direct component of surface albedo when using SIB LSM (-).
    alswnirdif   ! near-IR diffuse component of surface albedo when using SIB LSM (-).

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte),optional:: &
    cldfra2d,   &! cloud fraction at middle of layers (-).
    qc2d,       &! cloud liquid water mixing ratio at middle of layers (kg/kg).
    qr2d,       &! rain mixing ratio at middle of layers (kg/kg).
    qi2d,       &! cloud ice mixing ratio at middle of layers (kg/kg).
    qs2d         ! snow mixing ratio (kg/kg).

 real(kind=kind_phys),intent(in),dimension(1:noznlevels),optional:: &
    pin          ! fixed pressures at which climatological ozone data are defined (hPa).

 real(kind=kind_phys),intent(in),dimension(its:ite,1:noznlevels):: &
    o3clim       ! monthly-mean climatological ozone data at fixed pressure levels.


!output arguments:
 character(len=*),intent(out):: &
    errmsg       ! output error message (-).

 logical,intent(out),dimension(its:ite):: &
    dorrsw       ! logical set to true for daytime grid-points (-).

 integer,intent(out):: &
    errflg       ! output error flag (-).

 integer,intent(out):: &
    ncol,       &! number of columns in each MPI domain.
    nlay,       &! number of layers plus one (between the top of the model and top-of-the-atmosphere).
    dyofyr       ! day of the year.

 integer,intent(out):: &
    aer_opt      ! local option for shortwave optical properties of aerosols.

 integer,intent(out):: &
    inflgsw,    &! local flag for calculation of shortwave cloud optical properties.
    iceflgsw,   &! local flag for calculation of shortwave cloud ice particle properties.
    liqflgsw     ! local flag for calculation of shortwave cloud liquid droplet properties.

 real(kind=kind_phys),intent(out):: &
    adjes        ! flux adjustment for Earth/Sun distance (-).

 real(kind=kind_phys),intent(out),dimension(its:ite):: &
    coszen,     &! local cosine of solar zenith angle (-).
    asdir,      &! local visible direct surface albedo (-).
    asdif,      &! local visible diffuse surface albedo (-).
    aldir,      &! local near-IR direct surface albedo (-).
    aldif        ! local near-IR diffuse surface albedo (-).

 real(kind=kind_phys),intent(out),dimension(its:ite,kts:kme):: &
    play,       &! local pressure at middle of layers (hPa).
    tlay,       &! local temperature at middle of layers (K).
    h2ovmr,     &! local volume mixing ratio of water vapor ().
    o3vmr,      &! local volume mixing ratio of ozone ().
    co2vmr,     &! local volume mixing ratio of carbon dioxide ().
    o2vmr,      &! local volume mixing ratio of oxygen ().
    ch4vmr,     &! local volume mixing ratio of methane ().
    n2ovmr       ! local volume mixing ratio of nitrous oxide ().

 real(kind=kind_phys),intent(out),dimension(its:ite,kts:kme):: &
    relqmcl,    &! local cloud liquid water particle size (microns).
    reicmcl,    &! local cloud ice particle size (microns).
    resnmcl      ! local snow particle size (microns).

 real(kind=kind_phys),intent(out),dimension(its:ite,kts:kte+2):: &
    plev,       &! local pressure at interface between layers (hPa).
    tlev         ! local temperature at interface between layers (K).

 real(kind=kind_phys),intent(out),dimension(ngptsw,its:ite,kts:kme):: &
    cldfmcl,    &! in-cloud cloud fraction (-).
    clwpmcl,    &! in-cloud liquid water path (g/m2).
    ciwpmcl,    &! in-cloud ice water path (g/m2).
    cswpmcl,    &! in-cloud snow path (g/m2).
    taucmcl,    &! in-cloud optical depth (-).
    ssacmcl,    &! in-cloud single scattering albedo (-).
    asmcmcl,    &! in-cloud asymmetry factor (-).
    fsfcmcl      ! in-cloud forward scattering fraction (-).

 real(kind=kind_phys),intent(out),dimension(its:ite,kts:kme,nbndsw):: &
    tauaer,     &! aerosol optical depth (-).
    ssaaer,     &! aerosol single scattering albedo (-).
    asmaer       ! aerosol asymmetry factor (-).

 real(kind=kind_phys),intent(out),dimension(its:ite,kts:kme,naerec):: &
    ecaer        ! aerosol optical depth for ECMWF aerosols ().


!local variables and arrays:
 integer:: i,k,n
 integer:: pcols,pver
 integer:: idx_rei
 integer:: iplon,irng,permuteseed

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

 real(kind=kind_phys),dimension(its:ite,kts:kte):: pdel
 real(kind=kind_phys),dimension(its:ite,kts:kte+1):: hgt,cldfrac
 real(kind=kind_phys),dimension(its:ite,kts:kte+1):: ciwpth,clwpth,cswpth,rei,rel,res
 real(kind=kind_phys),dimension(nbndsw,its:ite,kts:kte+1):: taucld,ssacld,asmcld,fsfcld

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

!---
 real(kind=kind_phys):: declin0,scon0

!-----------------------------------------------------------------------------------------------------------------

!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine rrtmg_sw_mpas_timestep_init:')

!--- number of columns and number of layers between the surface and the top-of-the-atmosphere (i.e. number of
!    MPAS layers plus one):
 ncol = ite-its+1
 nlay = (kte-kts+1) + 1

!--- initialize option for the calculation of the cloud optical properties:
 inflgsw  = 2
 iceflgsw = 3
 liqflgsw = 1

!--- initialize option for the calculation of the aerosol optical properties:
 aer_opt = 0

!--- calculate solar declination and set solar constant:
 dyofyr = 0
 adjes  = 1.0

 call radconst(declin0,scon0,julian,degrad,dpd,solcon_0)
!call mpas_log_write('--- declin = $r  declin0 = $r  solcon_0 = $r',realArgs=(/declin,declin0,solcon_0/)) 

!--- calculation of the local input soundings for daytime columns:
 i_loop: do i = its,ite

    !calculation of the solar zenith angle:
    xt24 = mod(xtime+radt*0.5,1440.)
    tloctm = gmt + xt24/60. + xlon(i)/15.
    hrang = 15. * (tloctm-12.) * degrad
    xxlat = xlat(i) * degrad
!   coszrs = sin(xxlat) * sin(declin) + cos(xxlat) * cos(declin) * cos(hrang)
    coszrs = sin(xxlat) * sin(declin0) + cos(xxlat) * cos(declin0) * cos(hrang)
    coszen(i) = coszrs

    dorrsw(i) = .true.
    if(coszen(i) .le. 0._kind_phys) dorrsw(i) = .false.

    if(dorrsw(i)) then !begin if statement checking for daytime columns.

       !initialization of the surface albedo:
       if(present(alswvisdir) .and. present(alswvisdir) .and. &
          present(alswnirdir) .and. present(alswnirdif)) then
          asdir(i) = alswvisdir(i)
          asdif(i) = alswvisdif(i)
          aldir(i) = alswnirdir(i)
          aldif(i) = alswnirdif(i)
       else
          asdir(i) = albedo(i)
          asdif(i) = albedo(i)
          aldir(i) = albedo(i)
          aldif(i) = albedo(i)
       endif

       !initialization of the local input soundings:
       k = kts
       plev(i,kts) = p2di(i,k) / 100.
       tlev(i,kts) = t2di(i,k)
       do k = kts,kte
          plev(i,k+1) = p2di(i,k+1) / 100.
          tlev(i,k+1) = t2di(i,k+1)
          play(i,k)   = p2d(i,k) / 100.
          pdel(i,k)   = plev(i,k) - plev(i,k+1)
          tlay(i,k)   = t2d(i,k)
          h2ovmr(i,k) = amax1(qv2d(i,k),1.e-12)
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
             o3clim1d(k) = o3clim(i,k)
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

       !initialization of the aerosol cloud optical properties:
       do k = kts,kte+1
          do n = 1,nbndsw
             tauaer(i,k,n) = 0._kind_phys
             ssaaer(i,k,n) = 1._kind_phys
             asmaer(i,k,n) = 0._kind_phys
          enddo
          do n = 1,naerec
             ecaer(i,k,n) = 0._kind_phys
          enddo
       enddo
    endif !end if statement checking for daytime columns.

 enddo i_loop


!--- calculation of the cloud optical properties for daytime columns:
 do i = its,ite
    !initialization of input variables:
    do k = kts,kte
       cldfrac(i,k)  = 0._kind_phys
       ciwpth(i,k)   = 0._kind_phys
       clwpth(i,k)   = 0._kind_phys
       cswpth(i,k)   = 0._kind_phys
       rei(i,k)      = 0._kind_phys
       rel(i,k)      = 0._kind_phys
       res(i,k)      = 0._kind_phys
       do n = 1, nbndsw
          taucld(n,i,k) = 0._kind_phys
          ssacld(n,i,k) = 1._kind_phys
          asmcld(n,i,k) = 0._kind_phys
          fsfcld(n,i,k) = 0._kind_phys
       enddo
    enddo
    k = kte+1
    cldfrac(i,k) = 0._kind_phys
    ciwpth(i,k)  = 0._kind_phys
    clwpth(i,k)  = 0._kind_phys
    cswpth(i,k)  = 0._kind_phys
    rei(i,k)     = 10._kind_phys
    rel(i,k)     = 10._kind_phys
    res(i,k)     = 10._kind_phys
    do n = 1,nbndsw
       taucld(n,i,k) = 0._kind_phys
       ssacld(n,i,k) = 1._kind_phys
       asmcld(n,i,k) = 0._kind_phys
       fsfcld(n,i,k) = 0._kind_phys
    enddo
    !initialization of output variables:
    do k = kts,kte+1
       relqmcl(i,k) = 0._kind_phys
       reicmcl(i,k) = 0._kind_phys
       resnmcl(i,k) = 0._kind_phys
       do n = 1,ngptsw
          cldfmcl(n,i,k) = 0._kind_phys
          clwpmcl(n,i,k) = 0._kind_phys
          ciwpmcl(n,i,k) = 0._kind_phys
          cswpmcl(n,i,k) = 0._kind_phys
          taucmcl(n,i,k) = 0._kind_phys
          ssacmcl(n,i,k) = 0._kind_phys
          asmcmcl(n,i,k) = 0._kind_phys
          fsfcmcl(n,i,k) = 0._kind_phys
       enddo
    enddo
 enddo

 if(icloud .gt. 0) then !begin cloud effect to shortwave radiation.
    i_loop2: do i = its,ite

       if(dorrsw(i)) then !begin if statement checking for daytime columns.
          dzsum = 0.0
          do k = kts, kte
             dz = dz8w2d(i,k)
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
          if(present(f_qc) .and. present(qc2d)) then
             do k = kts,kte
                qc1d(k) = amax1(qc2d(i,k),0.)
             enddo
          endif
          if(present(f_qr) .and. present(qr2d)) then
             do k = kts,kte
                qr1d(k) = amax1(qr2d(i,k),0.)
             enddo
          endif
          if(present(f_qi) .and. present(qi2d)) then
             do k = kts,kte
                qi1d(k) = amax1(qi2d(i,k),0.)
             enddo
          endif
          if(present(f_qs) .and. present(qs2d)) then
             do k = kts,kte
                qs1d(k) = amax1(qs2d(i,k),0.)
             enddo
          endif
          if(present(f_qc).and.present(qc2d) .and. present(f_qi).and.present(qi2d)  .and. &
             present(f_qs).and.present(qs2d) .and. present(cldfra2d)) then
             do k = kts,kte
                cldfrac(i,k) = cldfra2d(i,k)
             enddo
          endif

          if(inflgsw .gt. 0) then
             !--- calculation of the effective radii for cloud liquid water, cloud ice, and snow:
             if(has_reqc .ne. 0) then
                !fill the effective radius for cloud water with that calculated in the Thompson cloud
                !microphysics parameterization:
                inflgsw = 3
                do k = kts,kte
                   rel(i,k) = max(2.5,re_cloud(i,k)*1.e6)
                   if(rel(i,k).le.2.5 .and. cldfrac(i,k).gt.0. .and. &
                      (xland(i)-1.5).gt.0.) then        !--- ocean.
                      rel(i,k) = 10.5
                   elseif(rel(i,k).le.2.5 .and. cldfrac(i,k).gt.0. .and. &
                      (xland(i)-1.5).lt.0.) then        !--- land.
                      rel(i,k) = 7.5
                   endif
                enddo
             else
                !calculate the effective radius for cloud water in layers below the model top. we added the
                !dimensions pcols, pver so that the calls to subroutines relcalc and reicalc remain the same
                !as in WRF although the two dimensions are duplicate and not needed inside the subroutines:
                pcols = 1
                pver  = nlay-1
                landfrac(pcols) = 2.-xland(i)
                landm(pcols)    = landfrac(pcols)
                snowh(pcols)    = 0.001*snow(i)
                icefrac(pcols)  = xice(i)
                do k = kts,kte
                   t1d(pcols,k) = tlay(i,k)
                enddo
                call relcalc(pcols,pcols,pver,t1d,landfrac,landm,icefrac,rel1d,snowh)
                do k = kts,kte
                   rel(i,k) = rel1d(pcols,k)
                enddo
             endif
             if(has_reqi .ne. 0) then
                !fill the effective radius for cloud ice with that calculated in the Thompson cloud
                !microphysics parameterization:
                inflgsw  = 4
                iceflgsw = 4
                do k = kts,kte
                   rei(i,k) = max(5.,re_ice(i,k)*1.e6)
                   if(rei(i,k).le.5. .and. cldfrac(i,k).gt.0.) then
                      idx_rei = int(tlay(i,k)-179.)
                      idx_rei = min(max(idx_rei,1),75)
                      corr = tlay(i,k) - int(tlay(i,k))
                      rei(i,k) = retab(idx_rei)*(1.-corr) + retab(idx_rei+1)*corr
                      rei(i,k) = max(rei(i,k),5.0)
                   endif
                enddo
             else
                !calculate the effective radius for cloud ice in layers below the model top, following
                !Kristjansson and Mitchell. limit upper bound of rei for Fu ice parameterization and
                !convert effective radius to generalized effective size (*1.0315; Fu, 1996):
                pcols = 1
                pver  = nlay-1
                do k = kts,kte
                   t1d(pcols,k) = tlay(i,k)
                enddo
                call reicalc(pcols,pcols,pver,t1d,rei1d)
                do k = kts,kte
                   rei(i,k) = rei1d(pcols,k)
                enddo
                if(iceflgsw .eq. 3) then
                   do k = kts, kte
                      rei(i,k) = rei(i,k) * 1.0315
                      rei(i,k) = min(140.0,rei(i,k))
                   enddo
                endif
             endif
             if(has_reqs .ne. 0) then
                !fill the effective radius for snow with that calculated in the Thompson cloud
                !microphysics parameterization:
                inflgsw  = 5
                iceflgsw = 5
                do k = kts, kte
                   res(i,k) = max(10.,re_snow(i,k)*1.e6)
                enddo
             else
                do k = kts, kte
                   res(i,k) = 10.
                enddo
             endif
          endif

          !--- calculation of the cloud liquid water, cloud ice, and snow paths (units are g/m2):
          !cloud liquid water path:
          do k = kts,kte
             gliqwp = qc1d(k) * pdel(i,k)*100.0 / gravity * 1000.0              !grid box liquid water path.
             clwpth(i,k) = gliqwp / max(0.01,cldfrac(i,k))                      !in-cloud liquid water path.
          enddo
          !cloud ice path:
          if(iceflgsw .eq. 3) then
             do k = kts,kte
                gicewp = (qi1d(k)+qs1d(k)) * pdel(i,k)*100.0 / gravity * 1000.0 !grid box ice water path.
                ciwpth(i,k) = gicewp / max(0.01,cldfrac(i,k))                   !in-cloud ice water path.
             enddo
          elseif(iceflgsw .ge. 4) then
             do k = kts,kte
                gicewp = qi1d(k) * pdel(i,k)*100.0 / gravity * 1000.0           !grid box ice water path.
                ciwpth(i,k) = gicewp / max(0.01,cldfrac(i,k))                   !in-cloud ice water path.
             enddo
          endif
          !snow path:
          if(iceflgsw.eq.5) then
             do k = kts,kte
                snow_mass_factor = 1.0
                if(res(i,k) .gt. 130.)then
                   snow_mass_factor = (130.0/res(i,k))*(130.0/res(i,k))
                   res(i,k)   = 130.0
                endif
                gsnowp = qs1d(k) * pdel(i,k)*100.0 / gravity * 1000.0          !grid box snow path.
                cswpth(i,k) = snow_mass_factor * gsnowp / max(0.01,cldfrac(i,k))
             enddo
          endif
       endif !end if statement checking for daytime columns.

    enddo i_loop2

    !--- call to subroutine mcica_subcol_sw which outputs the stochastic cloud optical properties:
    iplon = 1
    irng  = 0
    permuteseed = 1
    call mcica_subcol_sw &
               (iplon   , ncol    , nlay    , icld    , permuteseed , irng    , play    , &
                hgt     , cldfrac , ciwpth  , clwpth  , cswpth      , rei     , rel     , &
                res     , taucld  , ssacld  , asmcld  , fsfcld      , cldfmcl , ciwpmcl , &
                clwpmcl , cswpmcl , reicmcl , relqmcl , resnmcl     , taucmcl , ssacmcl , &
                asmcmcl , fsfcmcl)
 endif !end cloud effect to shortwave radiation.

 errmsg = 'rrtmg_sw_mpas_timestep_init OK'
 errflg = 0

 call mpas_log_write('--- end subroutine rrtmg_sw_mpas_timestep_init:')

 end subroutine rrtmg_sw_mpas_timestep_init

!=================================================================================================================
 subroutine rrtmg_sw_mpas_timestep_final &
                  (its        , ite       , kts        , kte        , kme        , &
                   dorrsw     , pi2d      , alswvisdir , alswvisdif , alswnirdir , &
                   alswnirdif , coszen    , swhr       , swhrc      , swuflx     , &
                   swuflxc    , swdflx    , swdflxc    , swuflxcln  , swdflxcln  , &
                   sibvisdir  , sibvisdif , sibnirdir  , sibnirdif  , swdkdir    , &
                   swdkdirc   , swdkdif   , gsw        , rthratensw , swcf       , &
                   swupt      , swuptc    , swdnt      , swdntc     , swupb      , &
                   swupbc     , swdnb     , swdnbc     , swuptcln   , swdntcln   , &
                   swupbcln   , swdnbcln  , swvisdir   , swvisdif   , swnirdir   , &
                   swnirdif   , swddir    , swddni     , swddif     , swddirc    , &
                   swddnic    , swdownc   , swupflx    , swupflxc   , swdnflx    , &
                   swdnflxc   , errmsg    , errflg)
!=================================================================================================================

!input arguments:
 logical,intent(in),dimension(its:ite):: &
    dorrsw       ! logical set to true for daytime grid-points (-).

 integer,intent(in):: &
    its,ite,    &! first and last index for horizontal dimension.
    kts,kte,    &! first and last index for vertical dimension.
    kme          ! defined as kte+1.

 real(kind=kind_phys),intent(in),dimension(its:ite):: &
    coszen       ! local cosine of the solar zenith angle (-). 

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte+1):: &
    swhr,       &! local instantaneous all-sky shortwave heating rate (K/day).
    swhrc        ! local instantaneous clear-sky shortwave heating rate (K/day).

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte+2):: &
    swuflx,     &! local instantaneous all-sky upward shortwave radiation flux (W/m2).
    swuflxc,    &! local instantaneous clear-sky upward shortwave radiation flux (W/m2)
    swdflx,     &! local instantaneous all-sky downward shortwave radiation flux (W/m2).
    swdflxc      ! local instantaneous clear-sky downward shortwave radiation flux (W/m2).

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte+2):: &
    swuflxcln,  &! local instantaneous clean-sky upward shortwave radiation flux (W/m2).
    swdflxcln    ! local instantaneous clean-sky downward shortwave radiation flux (W/m2).

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte+2):: &
    sibvisdir,  &! local instantaneous direct visible shortwave radiation flux when using SIB LSM (W/m2).
    sibvisdif,  &! local instantaneous diffuse visible shortwave radiation flux when using SIB LSM (W/m2).
    sibnirdir,  &! local instantaneous direct near-IR shortwave radiation flux when using SIB LSM (W/m2).
    sibnirdif    ! local instantaneous diffuse near-IR shortwave radiation flux when using SIB LSM (W/m2).

 real(kind=kind_phys),dimension(its:ite,kts:kte+2):: &
    swdkdir,    &! local instantaneous all-sky downward direct shortwave irradiance (W/m2).
    swdkdirc,   &! local instantaneous clear-sky downward direct shortwave irradiance (W/m2).
    swdkdif      ! local instantaneous all-sky downward diffuse shortwave irradiance (W/m2).

 real(kind=kind_phys),intent(in),dimension(its:ite,kts:kte):: &
    pi2d         ! exner function (-).

 real(kind=kind_phys),intent(in),dimension(its:ite),optional:: &
    alswvisdir, &! visible direct component of surface albedo when using SIB LSM (-).
    alswvisdif, &! visible diffuse component of surface albedo when using SIB LSM (-).
    alswnirdir, &! near-IR direct component of surface albedo when using SIB LSM (-).
    alswnirdif   ! near-IR diffuse component of surface albedo when using SIB LSM (-).


!output arguments:
 character(len=*),intent(out):: &
    errmsg       ! output error message (-).

 integer,intent(out):: &
    errflg       ! output error flag (-).

 real(kind=kind_phys),intent(out),dimension(its:ite):: &
    gsw          ! all-sky net surface radiation (W/m2).

 real(kind=kind_phys),intent(out),dimension(its:ite),optional:: &
    swcf,       &! instantaneous all-sky shortwave top-of-the-atmosphere cloud radiative forcing (W/m2).
    swupt,      &! instantaneous all-sky shortwave top-of-the-atmosphere upward radiation (W/m2).
    swuptc,     &! instantaneous clear-sky shortwave top-of-the-atmosphere upward radiation (W/m2).
    swdnt,      &! instantaneous all-sky shortwave top-of-the-atmosphere downward radiation (W/m2).
    swdntc,     &! instantaneous clear-sky shortwave top-of-the-atmosphere downward radiation (W/m2).
    swupb,      &! instantaneous all-sky shortwave surface upward radiation (W/m2).
    swupbc,     &! instantaneous clear-sky shortwave surface upward radiation (W/m2).
    swdnb,      &! instantaneous all-sky shortwave surface downward radiation (W/m2).
    swdnbc       ! instantaneous clear-sky shortwave surface downward radiation (W/m2).

 real(kind=kind_phys),intent(out),dimension(its:ite),optional:: &
    swuptcln,   &! instantaneous clean-sky shortwave top-of-the-atmosphere upward radiation (W/m2).
    swdntcln,   &! instantaneous clean-sky shortwave top-of-the-atmosphere downward radiation (W/m2).
    swupbcln,   &! instantaneous clean-sky shortwave surface upward radiation (W/m2).
    swdnbcln     ! instantaneous clean-sky shortwave surface downward radiation (W/m2).

 real(kind=kind_phys),intent(out),dimension(its:ite),optional:: &
    swvisdir,   &! instantaneous visible direct all-sky shortwave surface radiation when using SIB LSM (W/m2).
    swvisdif,   &! instantaneous visible diffuse all-sky shortwave surface radiation when using SIB LSM (W/m2).
    swnirdir,   &! instantaneous near-IR direct all-sky shortwave surface radiation when using SIB LSM (W/m2).
    swnirdif     ! instantaneous near-IR diffuse all-sky shortwave surface radiation when using SIB LSM (W/m2).

 real(kind=kind_phys),intent(out),dimension(its:ite),optional:: &
    swddir,     &! instantaneous direct all-sky shortwave surface downward irradiance (W/m2).
    swddni,     &! instantaneous direct all-sky shortwave surface downward normal irradiance (W/m2).
    swddif,     &! instantaneous diffuse all-sky shortwave surface downward irradiance (W/m2).
    swddirc,    &! instantaneous direct clear-sky shortwave surface downward irradiance (W/m2).
    swddnic,    &! instantaneous direct clear-sky shortwave surface downward normal irradiance (W/m2).
    swdownc      ! instantaneous clear-sky shortwave surface downward radiation (W/m2).

 real(kind=kind_phys),intent(out),dimension(its:ite,kts:kte):: &
    rthratensw   ! all-sky shortwave radiation heating rate (K/s).

 real(kind=kind_phys),intent(out),dimension(its:ite,kts:kme+1),optional:: &
    swupflx,    &! instantaneous all-sky shortwave upward radiation flux profile (W/m2).
    swupflxc,   &! instantaneous clear-sky shortwave upward radiation flux profile (W/m2).
    swdnflx,    &! instantaneous all-sky shortwave downward radiation flux profile (W/m2).
    swdnflxc     ! instantaneous clear-sky shortwave downward radiation flux profile (W/m2).


!local variables and arrays:
 integer:: i,k

!-----------------------------------------------------------------------------------------------------------------

 call mpas_log_write('--- enter subroutine rrtmg_sw_mpas_timestep_final:')

 i_loop: do i = its,ite

    if(dorrsw(i)) then
       gsw(i) = swdflx(i,1) - swuflx(i,1)
       do k = kts,kte
          rthratensw(i,k) = swhr(i,k)/86400.
          rthratensw(i,k) = rthratensw(i,k)/pi2d(i,k)
       enddo

       !optional arguments:
       if(present(swupt)) then
         swcf(i)   = (swdflx(i,kte+2)-swuflx(i,kte+2)) - (swdflxc(i,kte+2)-swuflxc(i,kte+2))
         swupt(i)  = swuflx(i,kte+2)
         swuptc(i) = swuflxc(i,kte+2)
         swdnt(i)  = swdflx(i,kte+2)
         swdntc(i) = swdflxc(i,kte+2)
         swupb(i)  = swuflx(i,1)
         swupbc(i) = swuflxc(i,1)
         swdnb(i)  = swdflx(i,1)
         swdnbc(i) = swdflxc(i,1)
       endif
       if(present(swupflx)) then
          do k = kts, kme+1
             swupflx(i,k)  = swuflx(i,k)
             swupflxc(i,k) = swuflxc(i,k)
             swdnflx(i,k)  = swdflx(i,k)
             swdnflxc(i,k) = swdflxc(i,k)
          enddo
       endif
       if(present(swuptcln) .and. present(swdntcln) .and. &
          present(swupbcln) .and. present(swdnbcln)) then
          swuptcln(i) = swuflxcln(i,kte+2)
          swdntcln(i) = swdflxcln(i,kte+2)
          swupbcln(i) = swuflxcln(i,1)
          swdnbcln(i) = swdflxcln(i,1)
       endif
       if(present(alswvisdir) .and. present(alswvisdif) .and. &
          present(alswnirdir) .and. present(alswnirdif)) then
          swvisdir(i) = sibvisdir(i,1)
          swvisdif(i) = sibvisdif(i,1)
          swnirdir(i) = sibnirdir(i,1)
          swnirdif(i) = sibnirdif(i,1)
       endif
       if(present(swddir) .and. present(swddif) .and. present(swddni) .and. &
          present(swddirc) .and. present(swddnic) .and. present(swdownc)) then
          swddir(i)  = swdkdir(i,1)
          swddni(i)  = swddir(i) / coszen(i)
          swddif(i)  = swdkdif(i,1)
          swdownc(i) = swdflxc(i,1)
          swddirc(i) = swdkdirc(i,1)
          swddnic(i) = swddirc(i) / coszen(i)
       endif
 
    else

       gsw(i) = 0._kind_phys
       do k = kts,kte
          rthratensw(i,k) = 0._kind_phys
       enddo

       !optional arguments:
       if(present(swupt)) then
          swcf(i)   = 0._kind_phys
          swupt(i)  = 0._kind_phys
          swuptc(i) = 0._kind_phys
          swdnt(i)  = 0._kind_phys
          swdntc(i) = 0._kind_phys
          swupb(i)  = 0._kind_phys
          swupbc(i) = 0._kind_phys
          swdnb(i)  = 0._kind_phys
          swdnbc(i) = 0._kind_phys
       endif
       if(present(swupflx)) then
          do k = kts, kte+1
             swupflx(i,k)  = 0._kind_phys
             swupflxc(i,k) = 0._kind_phys
             swdnflx(i,k)  = 0._kind_phys
             swdnflxc(i,k) = 0._kind_phys
          enddo
       endif
       if(present(swuptcln) .and. present(swdntcln) .and. &
          present(swupbcln) .and. present(swdnbcln)) then
          swuptcln(i) = 0._kind_phys
          swdntcln(i) = 0._kind_phys
          swupbcln(i) = 0._kind_phys
          swdnbcln(i) = 0._kind_phys
       endif
       if(present(alswvisdir) .and. present(alswvisdif) .and. &
          present(alswnirdir) .and. present(alswnirdif)) then
          swvisdir(i) = 0._kind_phys
          swvisdif(i) = 0._kind_phys
          swnirdir(i) = 0._kind_phys
          swnirdif(i) = 0._kind_phys
       endif
       if(present(swddir) .and. present(swddif) .and. present(swddni) .and. &
          present(swddirc) .and. present(swddnic) .and. present(swdownc)) then
          swddir(i)  = 0._kind_phys
          swddif(i)  = 0._kind_phys
          swdownc(i) = 0._kind_phys
          swddirc(i) = 0._kind_phys
          swddnic(i) = 0._kind_phys
       endif
    endif !endif daytime grid columns.

 enddo i_loop

 errmsg = 'rrtmg_sw_mpas_timestep_final OK'
 errflg = 0

 call mpas_log_write('--- end subroutine rrtmg_sw_mpas_timestep_final:')

 end subroutine rrtmg_sw_mpas_timestep_final

!=================================================================================================================
 end module rrtmg_sw_mpas
!=================================================================================================================



