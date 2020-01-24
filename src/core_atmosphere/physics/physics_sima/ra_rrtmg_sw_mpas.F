!=================================================================================================================
 module rrtmg_sw_mpas
 use mpas_kind_types,only: kind_phys => RKIND
 use mpas_log

use module_ra_rrtmg_vinterp,only: vinterp_ozn
use module_ra_rrtmg_lw,only: inirad
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
                   cldovrlp , gmt        , declin     , degrad     , radt       , &
                   xtime    , xlat       , xlon       , tsk        , albedo     , &
                   p3d      , t3d        , qv3d       , dz8w       , p8w        , &
                   t8w      , alswvisdir , alswvisdif , alswnirdir , alswnirdif , &
                   o3input  , noznlevels , pin        , o3clim     , ncol       , &
                   nlay     , icld       , inflgsw    , iceflgsw   , liqflgsw   , &
                   dyofyr   , adjes      , dorrsw     , coszen     , tsfc       , &
                   asdir    , asdif      , aldir      , aldif      , play       , &
                   tlay     , pdel       , h2ovmr     , o3vmr      , co2vmr     , &
                   o2vmr    , ch4vmr     , n2ovmr     , plev       , tlev)
!=================================================================================================================

!input arguments:
 integer,intent(in):: ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte

 integer,intent(in):: &
    cldovrlp,   &!
    noznlevels, &!
    o3input      !

 real(kind=kind_phys),intent(in):: &
    declin,     &!
    degrad,     &!
    gmt,        &!
    radt,       &!
    xtime        !

 real(kind=kind_phys),intent(in),dimension(ims:ime,jms:jme):: &
    xlat,       &!
    xlon,       &!
    albedo,     &! surface albedo (-).
    tsk          ! surface temperature (K).

 real(kind=kind_phys),intent(in),dimension(ims:ime,kms:kme,jms:jme):: &
    t3d,        &!
    p3d,        &!
    qv3d         !

 real(kind=kind_phys),intent(in),dimension(ims:ime,jms:jme),optional:: &
   alswvisdir,  &!
   alswvisdif,  &!
   alswnirdir,  &!
   alswnirdif    !

 real(kind=kind_phys),intent(in),dimension(ims:ime,kms:kme,jms:jme):: &
    dz8w,       &!
    p8w,        &!
    t8w          !

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

 real(kind=kind_phys),intent(out),dimension(its:ite,kts:kte+1):: &
    play,     &!
    tlay,     &!
    h2ovmr,   &!

    o3vmr,    &!
    co2vmr,   &!
    o2vmr,    &!
    ch4vmr,   &!
    n2ovmr     !

 real(kind=kind_phys),intent(out),dimension(its:ite,kts:kte+2):: &
    plev,     &!
    tlev

!local variables and arrays:
 integer:: i,j,k
 integer:: pcols

 real(kind=kind_phys):: coszrs,hrang,tloctm,xt24,xxlat
 real(kind=kind_phys),dimension(kts:kte+1):: o3mmr,p1dlay
 real(kind=kind_phys),dimension(kts:kte+2):: p1dlev
 real(kind=kind_phys),dimension(1:noznlevels):: o3clim1d

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



