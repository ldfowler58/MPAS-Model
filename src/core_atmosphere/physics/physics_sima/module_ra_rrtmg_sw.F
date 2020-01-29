!=================================================================================================================
 module module_ra_rrtmg_sw

 use mpas_kind_types,only: kind_phys => RKIND
 use mpas_log
 use mpas_atmphys_constants,only : cp,gravity
 use module_ra_rrtmg_vinterp,only: vinterp_ozn
 use mpas_atmphys_utilities, only: physics_message,physics_error_fatal
#define FATAL_ERROR(M) call physics_error_fatal(M)
#define WRITE_MESSAGE(M) call physics_message(M)

 use rrtmg_sw_mpas,only: rrtmg_sw_mpas_timestep_init, &
                         rrtmg_sw_mpas_timestep_final

 use parrrsw,only: nbndsw, ngptsw, naerec
 use rrtmg_sw_rad,only: rrtmg_sw_run
 use mcica_subcol_gen_sw, only: mcica_subcol_sw

 use module_ra_rrtmg_lw, only : inirad, o3data, relcalc, reicalc, retab


 contains


!=================================================================================================================
 subroutine rrtmg_swrad( &
                       p3d,p8w,pi3d,t3d,t8w,dz8w,qv3d,qc3d,qr3d,       &
                       qi3d,qs3d,cldfra3d,o33d,tsk,albedo,             &
                       xland,xice,snow,coszr,xtime,gmt,julday,radt,    &
                       degrad,declin,scon,xlat,xlong,icloud,o3input,   &
                       noznlevels,pin,o3clim,gsw,swcf,rthratensw,      &
                       has_reqc,has_reqi,has_reqs,re_cloud,            &
                       re_ice,re_snow,                                 &
                       swupt,swuptc,swdnt,swdntc,                      &
                       swupb,swupbc,swdnb,swdnbc,                      &
                       swupflx,swupflxc,swdnflx,swdnflxc,              &
                       swuptcln,swdntcln,swupbcln,swdnbcln,            &
                       swddir,swddni,swddif,swddirc,swddnic,swdownc,   &
                       f_qc,f_qr,f_qi,f_qs,                            &
                       alswvisdir,alswvisdif,alswnirdir,alswnirdif,    &
                       swvisdir,swvisdif,swnirdir,swnirdif,            &
                       tauaer300,tauaer400,tauaer600,tauaer999,        &
                       gaer300,gaer400,gaer600,gaer999,                &
                       waer300,waer400,waer600,waer999,                &
                       cldovrlp,                                       &
                       errmsg,errflg,                                  &
                       ids,ide, jds,jde, kds,kde,                      & 
                       ims,ime, jms,jme, kms,kme,                      &
                       its,ite, jts,jte, kts,kte                       &
                       )

 implicit none

!=================================================================================================================

!--- input arguments:
 logical,intent(in),optional:: f_qc,f_qr,f_qi,f_qs

 integer,intent(in):: ids,ide, jds,jde, kds,kde, &
                      ims,ime, jms,jme, kms,kme, &
                      its,ite, jts,jte, kts,kte
 integer,intent(in):: icloud,has_reqc,has_reqi,has_reqs
 integer,intent(in):: julday
 integer,intent(in):: cldovrlp
 integer,intent(in):: o3input

 real(kind=kind_phys),intent(in):: radt,degrad,xtime,declin,scon,gmt
 real(kind=kind_phys),intent(in),dimension(ims:ime,jms:jme):: xlat,xlong
 real(kind=kind_phys),intent(in),dimension(ims:ime,jms:jme):: albedo,tsk,snow,xice,xland
 real(kind=kind_phys),intent(in),dimension(ims:ime,kms:kme,jms:jme):: t3d,p3d,pi3d,qv3d
 real(kind=kind_phys),intent(in),dimension(ims:ime,kms:kme,jms:jme):: dz8w,p8w,t8w
 real(kind=kind_phys),intent(in),dimension(ims:ime,kms:kme,jms:jme):: re_cloud,re_ice,re_snow

 real(kind=kind_phys),intent(in),dimension(ims:ime,jms:jme),optional:: &
    alswvisdir,alswvisdif,alswnirdir,alswnirdif
 real(kind=kind_phys),intent(in),dimension(ims:ime,kms:kme,jms:jme),optional:: cldfra3d,qc3d,qr3d,qi3d,qs3d,o33d

 real(kind=kind_phys),intent(in),dimension(ims:ime,kms:kme,jms:jme),optional:: &
    tauaer300,tauaer400,tauaer600,tauaer999,gaer300,gaer400,gaer600,gaer999,waer300,waer400,waer600,waer999

!--- additional input arguments to use the CAM ozone climatology:
 integer,intent(in):: noznlevels
 real(kind=kind_phys),intent(in),dimension(1:noznlevels),optional:: pin
 real(kind=kind_phys),intent(in),dimension(ims:ime,1:noznlevels,jms:jme),optional:: o3clim

!--- output arguments:
 character(len=*),intent(out):: errmsg
 integer,intent(out)::errflg

 real(kind=kind_phys),intent(out),dimension(ims:ime,jms:jme):: coszr,gsw
 real(kind=kind_phys),intent(out),dimension(ims:ime,kms:kme,jms:jme):: rthratensw

 real(kind=kind_phys),intent(out),dimension(ims:ime,jms:jme),optional:: &
      swcf,swupt,swuptc,swdnt,swdntc,swupb,swupbc,swdnb,swdnbc
 real(kind=kind_phys),intent(out),dimension(ims:ime,jms:jme),optional:: swuptcln,swdntcln,swupbcln,swdnbcln
 real(kind=kind_phys),intent(out),dimension(ims:ime,jms:jme),optional:: &
      swddir,swddni,swddif,swddirc,swddnic,swdownc
 real(kind=kind_phys),intent(out),dimension(ims:ime,jms:jme),optional:: &
      swvisdir,swvisdif,swnirdir,swnirdif

 real(kind=kind_phys),intent(out),dimension(ims:ime,kms:kme+2,jms:jme ),optional:: &
      swupflx,swupflxc,swdnflx,swdnflxc

!local variables and arrays:
 logical,dimension(its:ite):: dorrsw

 integer:: na,nb,ncol,nlay,icld,inflgsw,iceflgsw,liqflgsw
 integer:: dyofyr
 integer:: iplon,irng,permuteseed
 integer:: pcols,pver
 integer:: idx_rei
 integer:: i,j,k,n

 integer:: calc_clean_atm_diag
 integer:: sf_surface_physics
 integer:: aer_ra_feedback

 real(kind=kind_phys):: coszrs,xt24,tloctm,hrang,xxlat,adjes
 real(kind=kind_phys):: ro,dz
 real(kind=kind_phys):: corr
 real(kind=kind_phys):: gliqwp,gicewp,gsnowp
 real(kind=kind_phys):: snow_mass_factor
 real(kind=kind_phys):: dzsum

 real(kind=kind_phys),dimension(its:ite):: coszen,tsfc
 real(kind=kind_phys),dimension(its:ite,kts:kte):: pdel,cicewp,cliqwp,csnowp,reliq,reice,resnow
 real(kind=kind_phys),dimension(its:ite,kts:kte+1):: hgt

 real(kind=kind_phys),dimension(1):: landfrac,landm,snowh,icefrac
 real(kind=kind_phys),dimension(its:ite):: asdir,asdif,aldir,aldif

 real(kind=kind_phys),dimension(kts:kte):: tten1d,cldfra1d,qc1d,qr1d,qi1d,qs1d,o31d
 real(kind=kind_phys),dimension(kts:kte+1):: o3mmr,p1dlay
 real(kind=kind_phys),dimension(kts:kte+2):: p1dlev

 real(kind=kind_phys),dimension(1,kts:kte):: t1d,rei1d,rel1d

!--- additional local variables and arrays needed to include additional layers between the model top
!    and the top of the atmosphere:
 real(kind=kind_phys),dimension(its:ite,kts:kte+1):: play,tlay,h2ovmr,o3vmr,co2vmr,o2vmr,ch4vmr,n2ovmr
 real(kind=kind_phys),dimension(its:ite,kts:kte+1):: &
      clwpth,ciwpth,cswpth,rel,rei,res,cldfrac,relqmcl,reicmcl,resnmcl
 real(kind=kind_phys),dimension(its:ite,kts:kte+1):: swhr,swhrc

 real(kind=kind_phys),dimension(its:ite,kts:kte+2):: plev,tlev
 real(kind=kind_phys),dimension(its:ite,kts:kte+2):: swuflx,swdflx,swuflxc,swdflxc
 real(kind=kind_phys),dimension(its:ite,kts:kte+2):: swuflxcln,swdflxcln
 real(kind=kind_phys),dimension(its:ite,kts:kte+2):: sibvisdir,sibvisdif,sibnirdir,sibnirdif
 real(kind=kind_phys),dimension(its:ite,kts:kte+2):: swdkdir,swdkdif,swdkdirc

 real(kind=kind_phys),dimension(its:ite,kts:kte+1,nbndsw):: tauaer,ssaaer,asmaer

 real(kind=kind_phys),dimension(nbndsw,its:ite,kts:kte+1):: taucld,ssacld,asmcld,fsfcld
 real(kind=kind_phys),dimension(ngptsw,its:ite,kts:kte+1):: &
      cldfmcl,clwpmcl,ciwpmcl,cswpmcl,taucmcl,ssacmcl,asmcmcl,fsfcmcl

!--- additional local variables and arrays needed for the CAM ozone climatologyL
 real(kind=kind_phys),dimension(1:noznlevels):: o3clim1d

!--- additional local variables related to the implementation of aerosols in rrtmg_swrad in WRF 3.8.
 integer:: aer_opt 
 real(kind=kind_phys),dimension(its:ite,kts:kte+1,naerec):: ecaer

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

!--- additional data:
 character(len=200):: msg
 real(kind=kind_phys),parameter:: thresh=1.e-9
 real(kind=kind_phys):: ang,slope
 real(kind=kind_phys),dimension(nbndsw),save :: wavemin ! minimum wavelength (um) of 14 intervals.
 real(kind=kind_phys),dimension(nbndsw),save :: wavemax ! maximum wavelength (um) of 14 intervals.
 real(kind=kind_phys),dimension(nbndsw):: wavemid
 data wavemin /3.077,2.500,2.150,1.942,1.626,1.299,1.242,0.778,0.625,0.442,0.345,0.263,0.200,3.846/
 data wavemax /3.846,3.077,2.500,2.150,1.942,1.626,1.299,1.242,0.778,0.625,0.442,0.345,0.263,12.195/

!-----------------------------------------------------------------------------------------------------------------

 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine rrtmg_swrad:')

 call rrtmg_sw_mpas_timestep_init &
                  (ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte, &
                   nbndsw   , ngptsw     ,                                        &
                   f_qc     , f_qr       , f_qi       , f_qs       , icloud     , &
                   has_reqc , has_reqi   , has_reqs   , cldovrlp   , gravity    , &
                   gmt      , declin     , degrad     , radt       , xtime      , &
                   xlat     , xlong      , tsk        , albedo     , xland      , &
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
                   rel      , rei        , res        , clwpth     , ciwpth     , &
                   cswpth   , hgt        , cldfrac    , cldfmcl    , clwpmcl    , &
                   ciwpmcl  , cswpmcl    , relqmcl    , reicmcl    , resnmcl    , &
                   taucmcl  , ssacmcl    , asmcmcl    , fsfcmcl)

!--- initialize calc_clean_atm_diag:
 calc_clean_atm_diag = 0
 if(present(swuptcln) .and. present(swdntcln) .and. present(swupbcln) .and. present(swdnbcln)) &
    calc_clean_atm_diag = 1

!--- initialize sf_surface_physics:
 sf_surface_physics = 0
 if(present(alswvisdir) .and. present(alswvisdif) .and. present(alswnirdir) .and. present(alswnirdif)) &
    sf_surface_physics = 8

!--- latitude loop:
 j_loop: do j = jts,jte

!--- longitude loop:

    i_loop: do i = its,ite

       !--- initialization of outputs:
       coszr(i,j) = 0.
       gsw(i,j)   = 0.
       do k = kts,kte
          rthratensw(i,k,j) = 0.
       enddo

       if(present(swupt)) then
          swcf(i,j) = 0.
          swupt(i,j)  = 0.
          swuptc(i,j) = 0.
          swdnt(i,j)  = 0.
          swdntc(i,j) = 0.
          swupb(i,j)  = 0.
          swupbc(i,j) = 0.
          swdnb(i,j)  = 0.
          swdnbc(i,j) = 0.
       endif
       if(calc_clean_atm_diag .gt. 0)then
          swuptcln(i,j) = 0.
          swdntcln(i,j) = 0.
          swupbcln(i,j) = 0.
          swdnbcln(i,j) = 0.
       endif
       if(sf_surface_physics .eq. 8) then
          swvisdir(i,j) = 0.
          swvisdif(i,j) = 0.
          swnirdir(i,j) = 0.
          swnirdif(i,j) = 0.
       endif
       if(present(swddir) .and. present(swddif) .and. present(swddni) .and. &
          present(swddirc) .and. present(swddnic) .and. present(swdownc)) then
          swddir(i,j)   = 0.
          swddni(i,j)   = 0.
          swddif(i,j)   = 0.
          swdownc(i, j) = 0.
          swddirc(i,j)  = 0.
          swddnic(i,j)  = 0.
       endif
       if(present(swupflx)) then
          do k = kts, kte+2
             swupflx(i,k,j)  = 0.
             swupflxc(i,k,j) = 0.
             swdnflx(i,k,j)  = 0.
             swdnflxc(i,k,j) = 0.
          enddo
       endif

       !--- initialization of input and output arrays in call to subroutine rrtmg_sw_run. since shortwave
       !    radiation computations are done for daytime columns only, these arrays are initialized when
       !    dorrsw is set to true:
       aer_opt = 0

       if(dorrsw(i)) then    !begin if statement checking for daytime columns.
          do k = kts,kte+1
             !inputs:
             do nb = 1,nbndsw
                tauaer(i,k,nb) = 0.
                ssaaer(i,k,nb) = 1.
                asmaer(i,k,nb) = 0.
             enddo
             do na = 1,naerec
                ecaer(i,k,na) = 0.
             enddo

             !outputs:
             swuflx(i,k)    = 0.
             swuflxc(i,k)   = 0.
             swdflx(i,k)    = 0.
             swdflxc(i,k)   = 0.
             swhr(i,k)      = 0.
             swhrc(i,k)     = 0.
             swuflxcln(i,k) = 0.
             swdflxcln(i,k) = 0.
             sibvisdir(i,k) = 0.
             sibvisdif(i,k) = 0.
             sibnirdir(i,k) = 0.
             sibnirdif(i,k) = 0.
             swdkdir(i,k)   = 0.
             swdkdirc(i,k)  = 0.
             swdkdif(i,k)   = 0.
          enddo
          k = kte+2
          swuflx(i,k)    = 0.
          swuflxc(i,k)   = 0.
          swdflx(i,k)    = 0.
          swdflxc(i,k)   = 0.
          swuflxcln(i,k) = 0.
          swdflxcln(i,k) = 0.
          sibvisdir(i,k) = 0.
          sibvisdif(i,k) = 0.
          sibnirdir(i,k) = 0.
          sibnirdif(i,k) = 0.
          swdkdir(i,k)   = 0.
          swdkdirc(i,k)  = 0.
          swdkdif(i,k)   = 0.
       endif

    end do i_loop !end longitude loop.

 end do j_loop !end latitude loop.                                           

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
!    shortwave surface and top-of-the-atmosphere all-sky and clear-sky fluxes (mainly diagnostics):
 call rrtmg_sw_mpas_timestep_final &
            (dorrsw,pi3d,alswvisdir,alswvisdif,alswnirdir,alswnirdif,coszr,coszen,           &
             gsw,swhr,swhrc,swuflx,swuflxc,swdflx,swdflxc,swuflxcln,swdflxcln,               &
             sibvisdir,sibvisdif,sibnirdir,sibnirdif,swdkdir,swdkdirc,swdkdif,               &
             rthratensw,swcf,swupt,swuptc,swdnt,swdntc,swupb,swupbc,swdnb,swdnbc,            &
             swuptcln,swdntcln,swupbcln,swdnbcln,swvisdir,swvisdif,swnirdir,swnirdif,        &
             swddir,swddni,swddif,swddirc,swddnic,swdownc,swupflx,swupflxc,swdnflx,swdnflxc, &
             ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)

 call mpas_log_write('--- end subroutine rrtmg_swrad:')
 call mpas_log_write(' ')

 end subroutine rrtmg_swrad

!=================================================================================================================
 end module module_ra_rrtmg_sw
!=================================================================================================================
