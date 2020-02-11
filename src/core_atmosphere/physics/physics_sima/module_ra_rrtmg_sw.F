!=================================================================================================================
 module module_ra_rrtmg_sw

 use mpas_kind_types,only: kind_phys => RKIND
 use mpas_log
 use mpas_atmphys_constants,only: gravity

 use parrrsw,only      : nbndsw,ngptsw,naerec
 use rrtmg_sw_mpas,only: rrtmg_sw_mpas_timestep_init, &
                         rrtmg_sw_mpas_timestep_final
 use rrtmg_sw_rad,only : rrtmg_sw_run

 implicit none
 public:: rrtmg_swrad


 contains


!=================================================================================================================
 subroutine rrtmg_swrad &
                  (ims        , ime        , jms        , jme        , kms        , kme      , &
                   its        , ite        , jts        , jte        , kts        , kte      , &
                   f_qc       , f_qr       , f_qi       , f_qs       , icloud     , cldovrlp , &
                   has_reqc   , has_reqi   , has_reqs   , noznlevels , o3input    , xtime    , &
                   declin     , gmt        , julday     , radt       , scon       , degrad   , &
                   xland      , xice       , snow       , xlat       , xlon       , tsk      , &
                   albedo     , alswvisdir , alswvisdif , alswnirdir , alswnirdif , pin      , &
                   o3clim     , p3d        , pi3d       , t3d        , qv3d       , qc3d     , &
                   qr3d       , qi3d       , qs3d       , cldfra3d   , o33d       , re_cloud , &
                   re_ice     , re_snow    , dz8w       , p8w        , t8w        , coszr    , &
                   rthratensw , gsw        , swcf       , swupt      , swuptc     , swdnt    , &
                   swdntc     , swupb      , swupbc     , swdnb      , swdnbc     , swupflx  , &
                   swupflxc   , swdnflx    , swdnflxc   , swuptcln   , swdntcln   , swupbcln , &
                   swdnbcln   , swddir     , swddni     , swddif     , swddirc    , swddnic  , &
                   swdownc    , swvisdir   , swvisdif   , swnirdir   , swnirdif   ,            &
                   tauaer300  , tauaer400  , tauaer600  , tauaer999  , gaer300    , gaer400  , &
                   gaer600    , gaer999    , waer300    , waer400    , waer600    , waer999  , &
                   errmsg     , errflg)
!=================================================================================================================

!--- input arguments:
 logical,intent(in),optional:: &
    f_qc,       &! if true, the physics package includes the cloud liquid water mixing ratio.
    f_qr,       &! if true, the physics package includes the rain mixing ratio.
    f_qi,       &! if true, the physics package includes the cloud ice mixing ratio.
    f_qs         ! if true, the physics package includes the snow mixing ratio.

 integer,intent(in):: ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte

 integer,intent(in):: &
    icloud,     &! option that defines the cloud effect on the SW optical depth.
    has_reqc,   &! = 1 if effective radius of cloud liquid water computed from cloud microphysics, = 0 otherwise.
    has_reqi,   &! = 1 if effective radius of cloud ice computed from cloud microphysics, = 0 otherwise.
    has_reqs,   &! = 1 if effective radius of snow computed from cloud microphysics, = 0 otherwise.
    cldovrlp,   &! option that defines the type of cloud overlap.
    noznlevels, &! number of levels of climatological ozone data.
    o3input,    &! ozone input option for shortwave radiation.
    julday       !julian day.

 real(kind=kind_phys),intent(in):: &
    declin,     &! solar declination (radian).
    degrad,     &! conversion factor from degree to radian (-).
    gmt,        &! GMT time.
    radt,       &! shortwave radiation timestep (s).
    xtime        ! MPAS run time.

 real(kind=kind_phys),intent(in)::scon

 real(kind=kind_phys),intent(in),dimension(ims:ime,jms:jme):: &
    xlat,       &! latitude (radian).
    xlon,       &! longitude (radian).
    xland,      &! land-ocean mask (1=land including sea-ice ; 2=ocean)
    xice,       &! fractional seaice amount (-).
    snow,       &! snow water equivalent (kg/m2).
    albedo,     &! surface albedo (-).
    tsk          ! surface temperature (K).

 real(kind=kind_phys),intent(in),dimension(ims:ime,kms:kme,jms:jme):: &
    t3d,        &! temperature at middle of layers (K).
    p3d,        &! pressure at middle of layers (Pa).
    pi3d,       &! exner function (-).
    qv3d,       &! water vapor mixing ratio at middle of layers (kg/kg).
    re_cloud,   &! effective radius of cloud liquid droplets at middle of layers (m).
    re_ice,     &! effective radius of cloud ice crystals at middle of layers (m). 
    re_snow      ! effective radius of snow crystals at middle of layers (m).

 real(kind=kind_phys),intent(in),dimension(ims:ime,kms:kme,jms:jme):: &
    dz8w,       &! thickness of layers (m).
    p8w,        &! pressure at interface of layers (Pa).
    t8w          ! temperature at interface of layers (K).

 real(kind=kind_phys),intent(in),dimension(ims:ime,jms:jme),optional:: &
    alswvisdir, &! visible direct component of surface albedo when using SIB LSM (-).
    alswvisdif, &! visible diffuse component of surface albedo when using SIB LSM (-).
    alswnirdir, &! near-IR direct component of surface albedo when using SIB LSM (-).
    alswnirdif   ! near-IR diffuse component of surface albedo when using SIB LSM (-).

 real(kind=kind_phys),intent(in),dimension(ims:ime,kms:kme,jms:jme),optional:: &
    cldfra3d,   &! cloud fraction at middle of layers (-).
    qc3d,       &! cloud liquid water mixing ratio at middle of layers (kg/kg).
    qr3d,       &! rain mixing ratio at middle of layers (kg/kg).
    qi3d,       &! cloud ice mixing ratio at middle of layers (kg/kg).
    qs3d,       &! snow mixing ratio (kg/kg).
    o33d         ! ozone mixing ratio (kg/kg).

 real(kind=kind_phys),intent(in),dimension(ims:ime,kms:kme,jms:jme),optional:: &
    tauaer300,  &! aerosol optical depth at ??? nm from WRF-chem.
    tauaer400,  &! aerosol optical depth at ??? nm from WRF-chem.
    tauaer600,  &! aerosol optical depth at ??? nm from WRF-chem.
    tauaer999,  &! aerosol optical depth at ??? nm from WRF-chem.
    gaer300,    &! asymmetry factor at ??? nm from WRF-chem.
    gaer400,    &! asymmetry factor at ??? nm from WRF-chem.
    gaer600,    &! asymmetry factor at ??? nm from WRF-chem.
    gaer999,    &! asymmetry factor at ??? nm from WRF-chem.
    waer300,    &! single scattering albedo at ??? nm from WRF-chem.
    waer400,    &! single scattering albedo at ??? nm from WRF chem.
    waer600,    &! single scattering albedo at ??? nm from WRF chem.
    waer999      !  single scattering albedo at ??? nm from WRF chem.

 real(kind=kind_phys),intent(in),dimension(1:noznlevels),optional:: &
    pin          ! fixed pressures at which climatological ozone data are defined (hPa).

 real(kind=kind_phys),intent(in),dimension(ims:ime,1:noznlevels,jms:jme),optional:: &
    o3clim       ! monthly-mean climatological ozone data at fixed pressure levels.


!--- output arguments:
 character(len=*),intent(out):: &
    errmsg       ! output error message (-).

 integer,intent(out):: &
    errflg       ! output error flag (-).

 real(kind=kind_phys),intent(out),dimension(ims:ime,jms:jme):: &
    coszr,      &! cosine of the solar zenith angle (-).
    gsw          ! all-sky net surface radiation (W/m2).

 real(kind=kind_phys),intent(out),dimension(ims:ime,kms:kme,jms:jme):: &
    rthratensw   ! all-sky shortwave radiation heating rate (K/s).

 real(kind=kind_phys),intent(out),dimension(ims:ime,jms:jme),optional:: &
    swcf,       &! instantaneous all-sky shortwave top-of-the-atmosphere cloud radiative forcing (W/m2).
    swupt,      &! instantaneous all-sky shortwave top-of-the-atmosphere upward radiation (W/m2).
    swuptc,     &! instantaneous clear-sky shortwave top-of-the-atmosphere upward radiation (W/m2).
    swdnt,      &! instantaneous all-sky shortwave top-of-the-atmosphere downward radiation (W/m2).
    swdntc,     &! instantaneous clear-sky shortwave top-of-the-atmosphere downward radiation (W/m2).
    swupb,      &! instantaneous all-sky shortwave surface upward radiation (W/m2).
    swupbc,     &! instantaneous clear-sky shortwave surface upward radiation (W/m2).
    swdnb,      &! instantaneous all-sky shortwave surface downward radiation (W/m2).
    swdnbc       ! instantaneous clear-sky shortwave surface downward radiation (W/m2).

 real(kind=kind_phys),intent(out),dimension(ims:ime,jms:jme),optional:: &
    swuptcln,   &! instantaneous clean-sky shortwave top-of-the-atmosphere upward radiation (W/m2).
    swdntcln,   &! instantaneous clean-sky shortwave top-of-the-atmosphere downward radiation (W/m2).
    swupbcln,   &! instantaneous clean-sky shortwave surface upward radiation (W/m2).
    swdnbcln     ! instantaneous clean-sky shortwave surface downward radiation (W/m2).

 real(kind=kind_phys),intent(out),dimension(ims:ime,jms:jme),optional:: &
    swvisdir,   &! instantaneous visible direct all-sky shortwave surface radiation when using SIB LSM (W/m2).
    swvisdif,   &! instantaneous visible diffuse all-sky shortwave surface radiation when using SIB LSM (W/m2).
    swnirdir,   &! instantaneous near-IR direct all-sky shortwave surface radiation when using SIB LSM (W/m2).
    swnirdif     ! instantaneous near-IR diffuse all-sky shortwave surface radiation when using SIB LSM (W/m2).

 real(kind=kind_phys),intent(out),dimension(ims:ime,kms:kme+1,jms:jme),optional:: &
    swupflx,    &! instantaneous all-sky shortwave upward radiation flux profile (W/m2).
    swupflxc,   &! instantaneous clear-sky shortwave upward radiation flux profile (W/m2).
    swdnflx,    &! instantaneous all-sky shortwave downward radiation flux profile (W/m2).
    swdnflxc     ! instantaneous clear-sky shortwave downward radiation flux profile (W/m2).

 real(kind=kind_phys),intent(out),dimension(ims:ime,jms:jme),optional:: &
    swddir,     &! instantaneous direct all-sky shortwave surface downward irradiance (W/m2).
    swddni,     &! instantaneous direct all-sky shortwave surface downward normal irradiance (W/m2).
    swddif,     &! instantaneous diffuse all-sky shortwave surface downward irradiance (W/m2).
    swddirc,    &! instantaneous direct clear-sky shortwave surface downward irradiance (W/m2).
    swddnic,    &! instantaneous direct clear-sky shortwave surface downward normal irradiance (W/m2).
    swdownc      ! instantaneous clear-sky shortwave surface downward radiation (W/m2).


!local variables only needed for conversion from ims->its,ime->ite,kms->kts,kme->kte,etc::
 real(kind=kind_phys),dimension(its:ite):: &
    xlat_hv,xlon_hv,xland_hv,xice_hv,snow_hv,albedo_hv,tsk_hv
 real(kind=kind_phys),dimension(its:ite,kts:kte):: &
    t3d_hv,p3d_hv,qv3d_hv,re_cloud_hv,re_ice_hv,re_snow_hv,dz8w_hv
 real(kind=kind_phys),dimension(its:ite,kts:kte+1):: &
    p8w_hv,t8w_hv
 real(kind=kind_phys),dimension(its:ite,kts:kte):: cldfra3d_hv,qc3d_hv,qr3d_hv,qi3d_hv,qs3d_hv


!local variables and arrays:
 logical,dimension(its:ite):: dorrsw

 integer:: i,j,k
 integer:: ncol,nlay,icld,inflgsw,iceflgsw,liqflgsw
 integer:: dyofyr
 integer:: calc_clean_atm_diag
 integer:: sf_surface_physics
 integer:: aer_ra_feedback
 integer:: aer_opt

 real(kind=kind_phys):: adjes

 real(kind=kind_phys),dimension(its:ite):: coszen,tsfc
 real(kind=kind_phys),dimension(its:ite):: asdir,asdif,aldir,aldif

 real(kind=kind_phys),dimension(its:ite,kts:kte+1):: play,tlay,h2ovmr,o3vmr,co2vmr,o2vmr,ch4vmr,n2ovmr
 real(kind=kind_phys),dimension(its:ite,kts:kte+1):: relqmcl,reicmcl,resnmcl
 real(kind=kind_phys),dimension(its:ite,kts:kte+1):: swhr,swhrc

 real(kind=kind_phys),dimension(its:ite,kts:kte+2):: plev,tlev
 real(kind=kind_phys),dimension(its:ite,kts:kte+2):: swuflx,swdflx,swuflxc,swdflxc
 real(kind=kind_phys),dimension(its:ite,kts:kte+2):: swuflxcln,swdflxcln
 real(kind=kind_phys),dimension(its:ite,kts:kte+2):: sibvisdir,sibvisdif,sibnirdir,sibnirdif
 real(kind=kind_phys),dimension(its:ite,kts:kte+2):: swdkdir,swdkdif,swdkdirc

 real(kind=kind_phys),dimension(its:ite,kts:kte+1,nbndsw):: tauaer,ssaaer,asmaer
 real(kind=kind_phys),dimension(its:ite,kts:kte+1,naerec):: ecaer

 real(kind=kind_phys),dimension(nbndsw,its:ite,kts:kte+1):: taucld,ssacld,asmcld,fsfcld
 real(kind=kind_phys),dimension(ngptsw,its:ite,kts:kte+1):: cldfmcl,clwpmcl,ciwpmcl,cswpmcl,taucmcl, &
                                                            ssacmcl,asmcmcl,fsfcmcl

!-----------------------------------------------------------------------------------------------------------------

 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine rrtmg_swrad:')


 do j = jts,jte
    do i = its,ite
       xlat_hv(i)   = xlat(i,j)
       xlon_hv(i)   =  xlon(i,j)
       xland_hv(i)  = xland(i,j)
       xice_hv(i)   = xice(i,j)
       snow_hv(i)   = snow(i,j)
       albedo_hv(i) = albedo(i,j)
       tsk_hv(i)    = tsk(i,j)
    enddo

    do k = kts,kte
       do i = its,ite
          t3d_hv(i,k)  = t3d(i,k,j)
          p3d_hv(i,k)  = p3d(i,k,j)
          qv3d_hv(i,k) = qv3d(i,k,j)
          dz8w_hv(i,k) = dz8w(i,k,j)
          p8w_hv(i,k)  = p8w(i,k,j)
          t8w_hv(i,k)  = t8w(i,k,j)
          re_cloud_hv(i,k) = re_cloud(i,k,j)
          re_ice_hv(i,k)   = re_ice(i,k,j)
          re_snow_hv(i,k)  = re_snow(i,k,j)
       enddo
    enddo
    k = kte+1
    do i = its,ite
       p8w_hv(i,k)  = p8w(i,k,j)
       t8w_hv(i,k)  = t8w(i,k,j)
    enddo

    !--- optional arrays:
    if(present(f_qc) .and. present(qc3d)) then
       do k = kts,kte
          do i = its,ite
             qc3d_hv(i,k) = qc3d(i,k,j)
          enddo
       enddo
    endif
    if(present(f_qr) .and. present(qr3d)) then
       do k = kts,kte
          do i = its,ite
             qr3d_hv(i,k) = qr3d(i,k,j)
          enddo
       enddo
    endif
    if(present(f_qi) .and. present(qi3d)) then
       do k = kts,kte
          do i = its,ite
             qi3d_hv(i,k) = qi3d(i,k,j)
          enddo
       enddo
    endif
    if(present(f_qs) .and. present(qs3d)) then
       do k = kts,kte
          do i = its,ite
             qs3d_hv(i,k) = qs3d(i,k,j)
          enddo
       enddo
    endif
    if(present(cldfra3d)) then
       do k = kts,kte
          do i = its,ite
             cldfra3d_hv(i,k) = cldfra3d(i,k,j)
          enddo
       enddo
    endif
 enddo


!--- initialization of all local variables needed to run rrtmg_sw_run:
 call rrtmg_sw_mpas_timestep_init &
                  (ims         , ime        , jms        , jme        , kms         , &
                   kme         , its        , ite        , jts        , jte         , &
                   kts         , kte        , nbndsw     , ngptsw     , naerec      , &
                   f_qc        , f_qr       , f_qi       , f_qs       , icloud      , &
                   has_reqc    , has_reqi   , has_reqs   , cldovrlp   , gravity     , &
                   gmt         , declin     , degrad     , radt       , xtime       , &
                   xlat_hv     , xlon_hv    , tsk_hv     , albedo_hv  , xland_hv    , &
                   xice_hv     , snow_hv    , p3d_hv     , t3d_hv     , qv3d_hv     , &
                   qc3d_hv     , qr3d_hv    , qi3d_hv    , qs3d_hv    , cldfra3d_hv , &
                   re_cloud_hv , re_ice_hv  , re_snow_hv , dz8w_hv    , p8w_hv      , &
                   t8w_hv      , alswvisdir , alswvisdif , alswnirdir , alswnirdif  , &
                   o3input     , noznlevels , pin        , o3clim     , ncol        , &
                   nlay        , icld       , inflgsw    , iceflgsw   , liqflgsw    , &
                   aer_opt     , dyofyr     , adjes      , dorrsw     , coszen      , &
                   tsfc        , asdir      , asdif      , aldir      , aldif       , &
                   play        , tlay       , h2ovmr     , o3vmr      , co2vmr      , &
                   o2vmr       , ch4vmr     , n2ovmr     , plev       , tlev        , &
                   cldfmcl     , clwpmcl    , ciwpmcl    , cswpmcl    , relqmcl     , &
                   reicmcl     , resnmcl    , taucmcl    , ssacmcl    , asmcmcl     , &
                   fsfcmcl     , tauaer     , ssaaer     , asmaer     , ecaer       , &
                   errmsg      , errflg)


!--- initialize calc_clean_atm_diag:
 calc_clean_atm_diag = 0
 if(present(swuptcln) .and. present(swdntcln) .and. present(swupbcln) .and. present(swdnbcln)) &
    calc_clean_atm_diag = 1

!--- initialize sf_surface_physics:
 sf_surface_physics = 0
 if(present(alswvisdir) .and. present(alswvisdif) .and. present(alswnirdir) .and. present(alswnirdif)) &
    sf_surface_physics = 8


!--- call the short-wave radiation code:
 call rrtmg_sw_run &
            (dorrsw    , ncol     , nlay      , icld      , play      , plev      , &
             tlay      , tlev     , tsfc      , h2ovmr    , o3vmr     , co2vmr    , &
             ch4vmr    , n2ovmr   , o2vmr     , asdir     , asdif     , aldir     , &
             aldif     , coszen   , adjes     , dyofyr    , scon      , inflgsw   , &
             iceflgsw  , liqflgsw , cldfmcl   , taucmcl   , ssacmcl   , asmcmcl   , &
             fsfcmcl   , ciwpmcl  , clwpmcl   , cswpmcl   , reicmcl   , relqmcl   , &
             resnmcl   , tauaer   , ssaaer    , asmaer    , ecaer     , swuflx    , &
             swdflx    , swhr     , swuflxc   , swdflxc   , swhrc     , swuflxcln , &
             swdflxcln , aer_opt  , sibvisdir , sibvisdif , sibnirdir , sibnirdif , &
             swdkdir   , swdkdif  , swdkdirc  , calc_clean_atm_diag   , errmsg    , &
             errflg)


!--- output the all-sky shortwave heating rate, shortwave net surface radiation for land surface scheme, and
!    optional shortwave surface and top-of-the-atmosphere all-sky and clear-sky fluxes (mainly diagnostics):
 call rrtmg_sw_mpas_timestep_final &
            (ims        , ime       , jms        , jme        , kms        , kme        , &
             its        , ite       , jts        , jte        , kts        , kte        , &
             dorrsw     , pi3d      , alswvisdir , alswvisdif , alswnirdir , alswnirdif , &
             coszen     , swhr      , swhrc      , swuflx     , swuflxc    , swdflx     , &       
             swdflxc    , swuflxcln , swdflxcln  , sibvisdir  , sibvisdif  , sibnirdir  , &
             sibnirdif  , swdkdir   , swdkdirc   , swdkdif    , coszr      , gsw        , &
             rthratensw , swcf      , swupt      , swuptc     , swdnt      , swdntc     , &
             swupb      , swupbc    , swdnb      , swdnbc     , swuptcln   , swdntcln   , &
             swupbcln   , swdnbcln  , swvisdir   , swvisdif   , swnirdir   , swnirdif   , &
             swddir     , swddni    , swddif     , swddirc    , swddnic    , swdownc    , &
             swupflx    , swupflxc  , swdnflx    , swdnflxc   , errmsg     , errflg)

 call mpas_log_write('--- end subroutine rrtmg_swrad:')
 call mpas_log_write(' ')

 end subroutine rrtmg_swrad

!=================================================================================================================
 end module module_ra_rrtmg_sw
!=================================================================================================================
