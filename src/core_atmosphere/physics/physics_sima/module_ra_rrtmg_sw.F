!=================================================================================================================
 module module_ra_rrtmg_sw

 use mpas_log
 use mpas_atmphys_constants,only : cp,gravity
 use module_ra_rrtmg_vinterp,only: vinterp_ozn
 use mpas_atmphys_utilities, only: physics_message,physics_error_fatal
#define FATAL_ERROR(M) call physics_error_fatal(M)
#define WRITE_MESSAGE(M) call physics_message(M)

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
                       degrad,declin,solcon,xlat,xlong,icloud,o3input, &
                       noznlevels,pin,o3clim,gsw,swcf,rthratensw,      &
                       has_reqc,has_reqi,has_reqs,re_cloud,            &
                       re_ice,re_snow,                                 &
                       swupt,swuptc,swdnt,swdntc,                      &
                       swupb,swupbc,swdnb,swdnbc,                      &
                       swupflx,swupflxc,swdnflx,swdnflxc,              &
                       swuptcln,swdntcln,swupbcln,swdnbcln,            & !--- CPF development.
                       swddir,swddni,swddif,swddirc,swddnic,swdownc,   &
                       f_qc,f_qr,f_qi,f_qs,                            &
                       alswvisdir,alswvisdif,alswnirdir,alswnirdif,    &
                       swvisdir,swvisdif,swnirdir,swnirdif,            &
                       tauaer300,tauaer400,tauaer600,tauaer999,        &
                       gaer300,gaer400,gaer600,gaer999,                &
                       waer300,waer400,waer600,waer999,                &
                       cldovrlp,                                       & !--- CPF development end.
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
 integer,intent(in),optional:: o3input

 real,intent(in):: radt,degrad,xtime,declin,solcon,gmt
 real,intent(in),dimension(ims:ime,jms:jme):: xlat,xlong
 real,intent(in),dimension(ims:ime,jms:jme):: albedo,tsk,snow,xice,xland
 real,intent(in),dimension(ims:ime,kms:kme,jms:jme):: t3d,p3d,pi3d,qv3d
 real,intent(in),dimension(ims:ime,kms:kme,jms:jme):: dz8w,p8w,t8w
 real,intent(in),dimension(ims:ime,kms:kme,jms:jme):: re_cloud,re_ice,re_snow

 real,intent(in),dimension(ims:ime,jms:jme),optional:: alswvisdir,alswvisdif,alswnirdir,alswnirdif
 real,intent(in),dimension(ims:ime,kms:kme,jms:jme),optional:: cldfra3d,qc3d,qr3d,qi3d,qs3d,o33d

 real,intent(in),dimension(ims:ime,kms:kme,jms:jme),optional:: &
    tauaer300,tauaer400,tauaer600,tauaer999,gaer300,gaer400,gaer600,gaer999,waer300,waer400,waer600,waer999

!--- additional input arguments to use the CAM ozone climatology:
 integer,intent(in):: noznlevels
 real,intent(in),dimension(1:noznlevels),optional:: pin
 real,intent(in),dimension(ims:ime,1:noznlevels,jms:jme),optional:: o3clim

!--- output arguments:
 character(len=*),intent(out):: errmsg
 integer,intent(out)::errflg

 real,intent(out),dimension(ims:ime,jms:jme):: coszr,gsw
 real,intent(out),dimension(ims:ime,kms:kme,jms:jme):: rthratensw

 real,intent(out),dimension(ims:ime,jms:jme),optional:: swcf,swupt,swuptc,swdnt,swdntc,swupb,swupbc,swdnb,swdnbc
 real,intent(out),dimension(ims:ime,jms:jme),optional:: swuptcln,swdntcln,swupbcln,swdnbcln
 real,intent(out),dimension(ims:ime,jms:jme),optional:: swddir,swddni,swddif,swddirc,swddnic,swdownc
 real,intent(out),dimension(ims:ime,jms:jme),optional:: swvisdir,swvisdif,swnirdir,swnirdif

 real,intent(out),dimension(ims:ime,kms:kme+2,jms:jme ),optional:: swupflx,swupflxc,swdnflx,swdnflxc

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
 integer:: cldovrlp

 real:: coszrs,xt24,tloctm,hrang,xxlat,adjes,scon
 real:: ro,dz
 real:: corr
 real:: gliqwp,gicewp,gsnowp
 real:: snow_mass_factor
 real:: dzsum

 real,dimension(its:ite):: coszen,tsfc
 real,dimension(its:ite,kts:kte):: pdel,cicewp,cliqwp,csnowp,reliq,reice,resnow
 real,dimension(its:ite,kts:kte+1):: hgt

 real,dimension(1):: landfrac,landm,snowh,icefrac
 real,dimension(its:ite):: asdir,asdif,aldir,aldif

 real,dimension(kts:kte):: tten1d,cldfra1d,qc1d,qr1d,qi1d,qs1d,o31d
 real,dimension(kts:kte+1):: o3mmr,p1dlay
 real,dimension(kts:kte+2):: p1dlev

 real,dimension(1,kts:kte):: t1d,rei1d,rel1d

!--- additional local variables and arrays needed to include additional layers between the model top
!    and the top of the atmosphere:
 real,dimension(its:ite,kts:kte+1):: play,tlay,h2ovmr,o3vmr,co2vmr,o2vmr,ch4vmr,n2ovmr
 real,dimension(its:ite,kts:kte+1):: clwpth,ciwpth,cswpth,rel,rei,res,cldfrac,relqmcl,reicmcl,resnmcl
 real,dimension(its:ite,kts:kte+1):: swhr,swhrc

 real,dimension(its:ite,kts:kte+2):: plev,tlev
 real,dimension(its:ite,kts:kte+2):: swuflx,swdflx,swuflxc,swdflxc
 real,dimension(its:ite,kts:kte+2):: swuflxcln,swdflxcln
 real,dimension(its:ite,kts:kte+2):: sibvisdir,sibvisdif,sibnirdir,sibnirdif
 real,dimension(its:ite,kts:kte+2):: swdkdir,swdkdif,swdkdirc

 real,dimension(its:ite,kts:kte+1,nbndsw):: tauaer,ssaaer,asmaer

 real,dimension(nbndsw,its:ite,kts:kte+1):: taucld,ssacld,asmcld,fsfcld
 real,dimension(ngptsw,its:ite,kts:kte+1):: cldfmcl,clwpmcl,ciwpmcl,cswpmcl,taucmcl,ssacmcl,asmcmcl,fsfcmcl

!--- additional local variables and arrays needed for the CAM ozone climatologyL
 real,dimension(1:noznlevels):: o3clim1d

!--- additional local variables related to the implementation of aerosols in rrtmg_swrad in WRF 3.8.
!    In WRF 3.8, these variables are in the argument list of subroutine rrtmg_swrad, but are made
!    local here:
 integer:: aer_opt 
 real,dimension(its:ite,kts:kte+1,naerec):: ecaer

!--- set trace gas volume mixing ratios, 2005 values, IPCC (2007):
!carbon dioxide (379 ppmv)
 real :: co2
 data co2 / 379.e-6 / 
!methane (1774 ppbv)
 real :: ch4
 data ch4 / 1774.e-9 / 
!nitrous oxide (319 ppbv)
 real :: n2o
 data n2o / 319.e-9 / 

!--- set oxygen volume mixing ratio (for o2mmr=0.23143):
 real :: o2
 data o2 / 0.209488 /

!--- molecular weights and ratios for converting mmr to vmr units
! real :: amd       ! Effective molecular weight of dry air (g/mol)  
! real :: amw       ! Molecular weight of water vapor (g/mol)        
! real :: amo       ! Molecular weight of ozone (g/mol)              
! real :: amo2      ! Molecular weight of oxygen (g/mol)              
! Atomic weights for conversion from mass to volume mixing ratios                
! data amd   /  28.9660   /                                                  
! data amw   /  18.0160   /                                                  
! data amo   /  47.9998   /                                                  
! data amo2  /  31.9999   /
                                                                                 
  real :: amdw     ! Molecular weight of dry air / water vapor  
  real :: amdo     ! Molecular weight of dry air / ozone
  real :: amdo2    ! Molecular weight of dry air / oxygen
  data amdw /  1.607793 /                                                    
  data amdo /  0.603461 /
  data amdo2 / 0.905190 /

!--- additional data:
 character(len=200):: msg
 real,parameter:: thresh=1.e-9
 real:: ang,slope
 real,dimension(nbndsw),save :: wavemin ! minimum wavelength (um) of 14 intervals.
 real,dimension(nbndsw),save :: wavemax ! maximum wavelength (um) of 14 intervals.
 real,dimension(nbndsw):: wavemid
 data wavemin /3.077,2.500,2.150,1.942,1.626,1.299,1.242,0.778,0.625,0.442,0.345,0.263,0.200,3.846/
 data wavemax /3.846,3.077,2.500,2.150,1.942,1.626,1.299,1.242,0.778,0.625,0.442,0.345,0.263,12.195/

!-----------------------------------------------------------------------------------------------------------------

!--- all fields are ordered vertically from bottom to top (pressures are in mb):
 ncol = ite-its+1
 nlay = (kte-kts+1) + 1

!--- initialize option for the calculation of the cloud optical properties:
 icld     = cldovrlp
 inflgsw  = 2
 iceflgsw = 3
 liqflgsw = 1

!--- initialize calc_clean_atm_diag:
 calc_clean_atm_diag = 0
 if(present(swuptcln) .and. present(swdntcln) .and. present(swupbcln) .and. present(swdnbcln)) &
    calc_clean_atm_diag = 1

!--- initialize sf_surface_physics:
 sf_surface_physics = 0
 if(present(alswvisdir) .and. present(alswvisdif) .and. present(alswnirdir) .and. present(alswnirdif)) &
    sf_surface_physics = 8

!--- set solar constant:
 scon   = solcon
 dyofyr = 0
 adjes  = 1.0 

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

       !--- initialization of variables in call to subroutine mcica_subcol_sw:
       do k = kts,kte+1
          !input arrays:
          cldfrac(i,k) = 0.
          clwpth(i,k)  = 0.
          ciwpth(i,k)  = 0.
          cswpth(i,k)  = 0.
          rel(i,k)     = 10.
          rei(i,k)     = 10.
          res(i,k)     = 10.
          !--- note: these variables are only calculated if inflgsw is set to zero (inactive option at present):
          do n = 1,nbndsw
             taucld(n,i,k) = 0.
             ssacld(n,i,k) = 1.
             asmcld(n,i,k) = 0.
             fsfcld(n,i,k) = 0.
          enddo
 
          !output arrays:
          relqmcl(i,k) = 0.
          reicmcl(i,k) = 0.
          resnmcl(i,k) = 0.
          do n = 1,ngptsw
             cldfmcl(n,i,k) = 0.
             clwpmcl(n,i,k) = 0.
             ciwpmcl(n,i,k) = 0.
             cswpmcl(n,i,k) = 0.
             taucmcl(n,i,k) = 0.
             ssacmcl(n,i,k) = 0.
             asmcmcl(n,i,k) = 0.
             fsfcmcl(n,i,k) = 0.
          enddo
       enddo

       !--- calculate the cosine of the solar zenith angle at the current time step to determine if the sun is
       !    above or below the horizon (xt24 is the fractional part of simulation days plus half of radt in 
       !    units of minutes, julian is in days, and radt is in minutes). do not call rrtmg_sw is night-time:
       xt24 = mod(xtime+radt*0.5,1440.)
       tloctm = gmt + xt24/60. + xlong(i,j)/15.
       hrang = 15. * (tloctm-12.) * degrad
       xxlat = xlat(i,j) * degrad
       coszrs = sin(xxlat) * sin(declin) + cos(xxlat) * cos(declin) * cos(hrang)
       coszr(i,j) = coszrs
       coszen(i)  = coszrs

       dorrsw(i) = .true.
       if (coszrs.le.0.0) dorrsw(i) = .false.


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

          !initialization of surface albedo:
          if(sf_surface_physics.eq.8 .and. xland(i,j).lt.1.5) then
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

          !initialization of input sounding:
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

          !--- initialization of ozone volume mixing ratio:
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

       endif    !begin if statement checking for daytime columns.

    end do i_loop !end longitude loop.


!--- calculate the cloud optical properties:
    i_loop2: do i = its,ite

       if(dorrsw(i)) then    !begin if statement checking for daytime columns.

          if(icloud .gt. 0) then !begin cloud effect to shortwave radiation.

             !---  mji - derive height of each layer mid-point from layer thickness. needed for exponential
             !     (icld=4) and exponential-random overlap option (icld=5) only.
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
                      gsnowp = qs1d(k) * pdel(n,k)*100.0 / gravity * 1000.0          !grid box snow path.
                      csnowp(i,k) = snow_mass_factor * gsnowp / max(0.01,cldfrac(i,k))
                   enddo
                endif

                !--- initialize input variables to subroutine mcica_subcol_sw. do not allow clouds between
                !    the model-top and the top of the atmosphere:
                do k = kts, kte
                   clwpth(i,k) = cliqwp(i,k)
                   ciwpth(i,k) = cicewp(i,k)
                   cswpth(i,k) = csnowp(i,k)
                   rel(i,k) = reliq(i,k)
                   rei(i,k) = reice(i,k)
                   res(i,k) = resnow(i,k)
                enddo
                clwpth(i,kte+1)  = 0.
                ciwpth(i,kte+1)  = 0.
                cswpth(i,kte+1)  = 0.
                rel(i,kte+1)     = 10.
                rei(i,kte+1)     = 10.
                res(i,kte+1)     = 10.

             endif !end inflgsw.

          endif !end cloud effect to shortwave radiation.

       endif    !end if statement checking for daytime columns.

    end do i_loop2 !end longitude loop.


!--- call the short-wave radiation code:
    iplon = 1
    irng  = 0
    permuteseed = 1
    call mcica_subcol_sw &
               (iplon   , ncol    , nlay    , icld    , permuteseed , irng    , play    , &
                hgt     , cldfrac , ciwpth  , clwpth  , cswpth      , rei     , rel     , &
                res     , taucld  , ssacld  , asmcld  , fsfcld      , cldfmcl , ciwpmcl , &
                clwpmcl , cswpmcl , reicmcl , relqmcl , resnmcl     , taucmcl , ssacmcl , &
                asmcmcl , fsfcmcl)

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


!--- outputs:
    i_loop3: do i = its,ite


       if(dorrsw(i)) then    !begin if statement checking for daytime columns.
          gsw(i,j)  = swdflx(i,1) - swuflx(i,1)

          if(present(swupt)) then
             swcf(i,j) = (swdflx(i,kte+2)-swuflx(i,kte+2)) - (swdflxc(i,kte+2)-swuflxc(i,kte+2))
             !output up and down toa fluxes for total and clear sky:
             swupt(i,j)  = swuflx(i,kte+2)
             swuptc(i,j) = swuflxc(i,kte+2)
             swdnt(i,j)  = swdflx(i,kte+2)
             swdntc(i,j) = swdflxc(i,kte+2)
             !output up and down surface fluxes for total and clear sky:
             swupb(i,j)  = swuflx(i,1)
             swupbc(i,j) = swuflxc(i,1)
             swdnb(i,j)  = swdflx(i,1)
             swdnbc(i,j) = swdflxc(i,1)
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
             swddir(i,j)   = swdkdir(i,1)          ! jararias 2013/08/10
             swddni(i,j)   = swddir(i,j) / coszrs  ! jararias 2013/08/10
             swddif(i,j)   = swdkdif(i,1)          ! jararias 2013/08/10
             swdownc(i, j) = swdflxc(i,1)          ! PAJ: clear-sky GHI
             swddirc(i,j)  = swdkdirc(i,1)         ! PAJ: clear-sky direct normal irradiance
             swddnic(i,j)  = swddirc(i,j) / coszrs ! PAJ: clear-sky direct normal irradiance
          endif
          if(present(swupflx)) then
             do k = kts, kte+2
                swupflx(i,k,j)  = swuflx(i,k)
                swupflxc(i,k,j) = swuflxc(i,k)
                swdnflx(i,k,j)  = swdflx(i,k)
                swdnflxc(i,k,j) = swdflxc(i,k)
             enddo
          endif

          !--- output heating rate tendency; convert heating rate from K/d to K/s:
          do k = kts, kte
             tten1d(k) = swhr(i,k)/86400.
             rthratensw(i,k,j) = tten1d(k)/pi3d(i,k,j)
          enddo

       endif    !end if statement checking for daytime columns.

    end do i_loop3

 end do j_loop !end latitude loop.                                           


 end subroutine rrtmg_swrad

!=================================================================================================================
 end module module_ra_rrtmg_sw
!=================================================================================================================
