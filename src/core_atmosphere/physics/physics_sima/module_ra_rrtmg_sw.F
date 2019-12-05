!=================================================================================================================
 module module_ra_rrtmg_sw

 use mpas_atmphys_constants,only : cp,g=>gravity
 use module_ra_rrtmg_vinterp,only: vinterp_ozn
 use mpas_atmphys_utilities, only: physics_message,physics_error_fatal
#define FATAL_ERROR(M) call physics_error_fatal(M)
#define WRITE_MESSAGE(M) call physics_message(M)

 use parrrsw,only: nbndsw, ngptsw, naerec
 use rrtmg_sw_rad,only: rrtmg_sw
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
                       ids,ide, jds,jde, kds,kde,                      & 
                       ims,ime, jms,jme, kms,kme,                      &
                       its,ite, jts,jte, kts,kte                       &
                       )

 implicit none

!=================================================================================================================

!--- input arguments:
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
 real,intent(in),dimension(ims:ime,kms:kme,jms:jme),optional:: &
    cldfra3d,qc3d,qr3d,qi3d,qs3d,o33d

!--- additional input arguments to use the CAM ozone climatology:
 integer,intent(in):: noznlevels
 real,intent(in),dimension(1:noznlevels),optional:: pin
 real,intent(in),dimension(ims:ime,1:noznlevels,jms:jme),optional:: o3clim

!--- inout arguments:
 real,intent(inout),dimension(ims:ime,jms:jme):: coszr,gsw,swcf
 real,intent(inout),dimension(ims:ime,jms:jme),optional:: &
    swupt,swuptc,swdnt,swdntc,swupb,swupbc,swdnb,swdnbc

 real,intent(inout),dimension(ims:ime,kms:kme,jms:jme):: rthratensw

!--- output arguments:
  real,intent(out),dimension(ims:ime,kms:kme+2,jms:jme ),optional:: &
     swupflx,swupflxc,swdnflx,swdnflxc

!local variables and arrays:
 logical:: dorrsw

 integer:: na,nb,ncol,nlay,icld,inflgsw,iceflgsw,liqflgsw
 integer:: dyofyr
 integer:: iplon,irng,permuteseed
 integer:: pcols,pver
 integer:: idx_rei
 integer:: i,j,k,n

 real:: coszrs,xt24,tloctm,hrang,xxlat,adjes,scon
 real:: ro,dz
 real:: corr
 real:: gliqwp,gicewp,gsnowp,gravmks
 real:: snow_mass_factor
 real,dimension(1):: tsfc,landfrac,landm,snowh,icefrac
 real,dimension(1):: asdir,asdif,aldir,aldif,coszen
 real,dimension(1,1:kte-kts+1):: pdel,cicewp,cliqwp,csnowp,reliq,reice,resnow
 real,dimension(kts:kte):: &
    tten1d,cldfra1d,dz1d,p1d,t1d,qv1d,qc1d,qr1d,qi1d,qs1d,o31d
 real,dimension(kts:kte+1):: pw1d,tw1d
 real,dimension(kts:kte+1):: o3mmr

!--- additional local variables and arrays needed to include additional layers between the model top
!    and the top of the atmosphere:
 real,dimension(1,kts:kte+1):: play,tlay,h2ovmr,o3vmr,co2vmr,o2vmr,ch4vmr,n2ovmr
 real,dimension(1,kts:kte+1):: clwpth,ciwpth,cswpth,rel,rei,res,cldfrac,relqmcl,reicmcl,resnmcl
 real,dimension(1,kts:kte+1):: swhr,swhrc

 real,dimension(1,kts:kte+2):: plev,tlev
 real,dimension(1,kts:kte+2):: swuflx,swdflx,swuflxc,swdflxc
 real,dimension(1,kts:kte+2):: swuflxcln,swdflxcln
 real,dimension(1,kts:kte+2):: sibvisdir,sibvisdif,sibnirdir,sibnirdif
 real,dimension(1,kts:kte+2):: swdkdir,swdkdif,swdkdirc

 real,dimension(1,kts:kte+1,nbndsw):: tauaer,ssaaer,asmaer

 real,dimension(nbndsw,1,kts:kte+1):: taucld,ssacld,asmcld,fsfcld
 real,dimension(ngptsw,1,kts:kte+1):: cldfmcl,clwpmcl,ciwpmcl,cswpmcl,taucmcl,ssacmcl,asmcmcl,fsfcmcl

!--- additional local variables and arrays needed for the CAM ozone climatologyL
 real,dimension(1:noznlevels):: o3clim1d

!--- additional local variables related to the implementation of aerosols in rrtmg_swrad in WRF 3.8.
!    In WRF 3.8, these variables are in the argument list of subroutine rrtmg_swrad, but are made
!    local here:
 integer:: aer_opt 
 real,dimension(1,kts:kte+1,naerec):: ecaer

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

!--- CPF development:
 logical,intent(in),optional:: f_qc,f_qr,f_qi,f_qs
 real,intent(in),dimension(ims:ime,jms:jme),optional:: alswvisdir,alswvisdif,alswnirdir,alswnirdif
 real,intent(in),dimension(ims:ime,kms:kme,jms:jme),optional:: &
    tauaer300,tauaer400,tauaer600,tauaer999,gaer300,gaer400,gaer600,gaer999, &
    waer300,waer400,waer600,waer999

 real,intent(out),dimension(ims:ime,jms:jme),optional:: swuptcln,swdntcln,swupbcln,swdnbcln
 real,intent(out),dimension(ims:ime,jms:jme),optional:: swddir,swddni,swddif,swddirc,swddnic,swdownc
 real,intent(out),dimension(ims:ime,jms:jme),optional:: swvisdir,swvisdif,swnirdir,swnirdif

 integer:: calc_clean_atm_diag
 integer:: sf_surface_physics
 integer:: aer_ra_feedback
 real:: dzsum
 real,dimension(1,kts:kte+1):: hgt

 character(len=200):: msg
 real,parameter:: thresh=1.e-9
 real:: ang,slope
 real,dimension(nbndsw),save :: wavemin ! minimum wavelength (um) of 14 intervals.
 real,dimension(nbndsw),save :: wavemax ! maximum wavelength (um) of 14 intervals.
 real,dimension(nbndsw):: wavemid
 data wavemin /3.077,2.500,2.150,1.942,1.626,1.299,1.242,0.778,0.625,0.442,0.345,0.263,0.200,3.846/
 data wavemax /3.846,3.077,2.500,2.150,1.942,1.626,1.299,1.242,0.778,0.625,0.442,0.345,0.263,12.195/

 integer:: cldovrlp

!-----------------------------------------------------------------------------------------------------------------

!--- all fields are ordered vertically from bottom to top (pressures are in mb):
 ncol = 1

!--- initialize option for the calculation of the cloud optical properties:
 icld     = cldovrlp
 inflgsw  = 2
 iceflgsw = 3
 liqflgsw = 1

!--- CPF development:
!--- initialize calc_clean_atm_diag:
 calc_clean_atm_diag = 0
 if(present(swuptcln) .and. present(swdntcln) .and. present(swupbcln) .and. present(swdnbcln)) &
    calc_clean_atm_diag = 1

!--- initialize sf_surface_physics:
 sf_surface_physics = 0
 if(present(alswvisdir) .and. present(alswvisdif) .and. present(alswnirdir) .and. present(alswnirdif)) &
    sf_surface_physics = 8

!--- initialize aer_ra_feedback:
 aer_ra_feedback = 0
 if(present(tauaer300) .and. present(gaer300) .and. present(waer300) .and. &
    present(tauaer400) .and. present(gaer400) .and. present(waer400) .and. &
    present(tauaer600) .and. present(gaer600) .and. present(waer600) .and. &
    present(tauaer999) .and. present(gaer999) .and. present(waer999)) &
    aer_ra_feedback = 1
!--- CPF development end.

!--- latitude loop:
 j_loop: do j = jts,jte

!--- longitude loop:
    i_loop: do i = its,ite

       !--- calculate the cosine of the solar zenith angle at the current time step to determine if the sun is
       !    above or below the horizon (xt24 is the fractional part of simulation days plus half of radt in 
       !    units of minutes, julian is in days, and radt is in minutes). do not call rrtmg_sw is night-time:
       dorrsw = .true.
       xt24 = mod(xtime+radt*0.5,1440.)
       tloctm = gmt + xt24/60. + xlong(i,j)/15.
       hrang = 15. * (tloctm-12.) * degrad
       xxlat = xlat(i,j) * degrad
       coszrs = sin(xxlat) * sin(declin) + cos(xxlat) * cos(declin) * cos(hrang)
       coszr(i,j) = coszrs

       if (coszrs.le.0.0) dorrsw = .false.
       
       if(dorrsw) then

          !--- INITIALIZE COLUMN SOUNDING (the call to the short wave radiation code is done one column
          !     at a time):
          do k = kts, kte+1
             pw1d(k) = p8w(i,k,j) / 100.
             tw1d(k) = t8w(i,k,j)
          enddo

          do k = kts, kte
             p1d(k)      = p3d(i,k,j) / 100.
             dz1d(k)     = dz8w(i,k,j)
             t1d(k)      = t3d(i,k,j)
             qv1d(k)     = amax1(qv3d(i,k,j),1.e-12)
             o31d(k)     = 0.
             qc1d(k)     = 0.
             qr1d(k)     = 0.
             qi1d(k)     = 0.
             qs1d(k)     = 0.
             cldfra1d(k) = 0.
             if(present(o33d)) o31d(k) = o33d(i,k,j)
          enddo

          !--- initialize the local arrays containing the different cloud water and ice condenstates:
          if(icloud .gt. 0) then
             if(present(f_qc) .and. present(qc3d)) then
                do k = kts, kte
                   qc1d(k) = amax1(qc3d(i,k,j),0.)
                enddo
             endif
             if(present(f_qr) .and. present(qr3d)) then
                do k = kts, kte
                   qr1d(k) = amax1(qr3d(i,k,j),0.)
                enddo
             endif
             if(present(f_qi) .and. present(qi3d)) then
                do k = kts, kte
                   qi1d(k) = amax1(qi3d(i,k,j),0.)
                enddo
             endif
             if(present(f_qs) .and. present(qs3d)) then
                do k = kts, kte
                   qs1d(k) = amax1(qs3d(i,k,j),0.)
                enddo
             endif
             if(present(f_qc).and.present(qc3d) .and. present(f_qi).and.present(qi3d)  .and. &
                present(f_qs).and.present(qs3d) .and. present(cldfra3d)) then
                do k = kts, kte
                   cldfra1d(k) = cldfra3d(i,k,j)
                enddo
             endif
          endif

          !--- initialize the local radiative heating rate:
          do k = kts, kte
             tten1d(k) = 0.
          enddo

          !--- add extra layer to include absorption between the top of the model and the top of the atmosphere:
          nlay = (kte-kts+1) + 1

          !--- initialize local arrays called in the calculation of the cloud optical properties and radiative
          !    fluxes:
          do n = 1, ncol
             do k = kts, kte
                reliq(ncol,k)  = 10.
                reice(ncol,k)  = 10.
                resnow(ncol,k) = 10.
                cliqwp(ncol,k) = 0.
                cicewp(ncol,k) = 0.
                csnowp(ncol,k) = 0.
             enddo

             do k = 1, nlay
                clwpth(n,k)    = 0.
                ciwpth(n,k)    = 0.
                cswpth(n,k)    = 0.
                rel(n,k)       = 0.
                rei(n,k)       = 0.
                res(n,k)       = 0.
                cldfrac(n,k)   = 0.
                relqmcl(n,k)   = 0.
                reicmcl(n,k)   = 0.
                resnmcl(n,k)   = 0.
                swuflx(n,k)    = 0.
                swuflxc(n,k)   = 0.
                swdflx(n,k)    = 0.
                swdflxc(n,k)   = 0.
                swhr(n,k)      = 0.
                swhrc(n,k)     = 0.
                swuflxcln(n,k) = 0.
                swdflxcln(n,k) = 0.
                taucld(1:nbndsw,n,k)  = 0.
                tauaer(n,k,1:nbndsw)  = 0.
                ssaaer(n,k,1:nbndsw)  = 0.
                asmaer(n,k,1:nbndsw)  = 0.
                cldfmcl(1:ngptsw,n,k) = 0.
                clwpmcl(1:ngptsw,n,k) = 0.
                ciwpmcl(1:ngptsw,n,k) = 0.
                cswpmcl(1:ngptsw,n,k) = 0.
                taucmcl(1:ngptsw,n,k) = 0.
             enddo
             swuflx(n,nlay+1)    = 0.
             swuflxc(n,nlay+1)   = 0.
             swdflx(n,nlay+1)    = 0.
             swdflxc(n,nlay+1)   = 0.
             swuflxcln(n,nlay+1) = 0.
             swdflxcln(n,nlay+1) = 0.

             do k = 1, nlay + 1
                sibvisdir(ncol,k) = 0.
                sibvisdif(ncol,k) = 0.                
                sibnirdir(ncol,k) = 0.
                sibnirdif(ncol,k) = 0.                
                swdkdir(n,k)      = 0.
                swdkdirc(n,k)     = 0.
                swdkdif(n,k)      = 0.
             enddo
          enddo

          !--- initialization of aerosol optical properties:
          aer_opt = 0
          do n = 1, ncol
             do k = 1, nlay
                do na = 1, naerec
                   ecaer(n,k,na)   = 0.
                enddo
             enddo
          enddo

          !--- fill local arrays with input sounding. convert water vapor mass mixing ratio to
          !    volume mixing ratio:
          plev(ncol,1) = pw1d(1)
          tlev(ncol,1) = tw1d(1)
          tsfc(ncol)   = tsk(i,j)
          do k = kts, kte
             plev(ncol,k+1)   = pw1d(k+1)
             tlev(ncol,k+1)   = tw1d(k+1)
             play(ncol,k)     = p1d(k)
             pdel(ncol,k)     = plev(ncol,k) - plev(ncol,k+1)
             tlay(ncol,k)     = t1d(k)
             h2ovmr(ncol,k)   = qv1d(k) * amdw
             co2vmr(ncol,k)   = co2
             o2vmr(ncol,k)    = o2
             ch4vmr(ncol,k)   = ch4
             n2ovmr(ncol,k)   = n2o
          enddo
          plev(ncol,kte+2)   = 1.0e-5
          tlev(ncol,kte+2)   = tlev(ncol,kte+1) + 0.0
          play(ncol,kte+1)   = 0.5 * plev(ncol,kte+1)
          tlay(ncol,kte+1)   = tlev(ncol,kte+1) + 0.0
          h2ovmr(ncol,kte+1) = h2ovmr(ncol,kte) 
          co2vmr(ncol,kte+1) = co2vmr(ncol,kte) 
          o2vmr(ncol,kte+1)  = o2vmr(ncol,kte) 
          ch4vmr(ncol,kte+1) = ch4vmr(ncol,kte) 
          n2ovmr(ncol,kte+1) = n2ovmr(ncol,kte) 

          !---  mji - derive height of each layer mid-point from layer thickness. needed for exponential
          !     (icld=4) and exponential-random overlap option (icld=5) only.
          dzsum = 0.0
          do k = kts, kte
             dz = dz1d(k)
             hgt(ncol,k) = dzsum + 0.5*dz
             dzsum = dzsum + dz
          enddo

          !--- initialize the ozone voume mixing ratio:
          call inirad(o3mmr,plev,kts,kte)
          if(o3input .eq. 2) then
             do k = 1, noznlevels
                o3clim1d(k) = o3clim(i,k,j)
             enddo
             call vinterp_ozn(1,ncol,ncol,kte+1,play,pin,noznlevels,o3clim1d,o3mmr)
             do k = kts,kte+1
                o3vmr(ncol,k) = o3mmr(k)
             enddo
          else
             do k = kts,kte+1
                o3vmr(ncol,k) = o3mmr(k)*amdo
             enddo
          endif

          !--- initialize the surface albedo:
          !CPF development:
          if(sf_surface_physics.eq.8 .and. xland(i,j).lt.1.5) then
             asdir(ncol) = alswvisdir(i,j)
             asdif(ncol) = alswvisdif(i,j)
             aldir(ncol) = alswnirdir(i,j)
             aldif(ncol) = alswnirdif(i,j)
          else
             asdir(ncol) = albedo(i,j)
             asdif(ncol) = albedo(i,j)
             aldir(ncol) = albedo(i,j)
             aldif(ncol) = albedo(i,j)
          endif
          !CPF development end.

          !--- set solar constant:
          scon = solcon
          !--- set cosine of solar zenith angle:
          coszen(ncol) = coszrs
          !--- in MPAS, the solar constant is already provided with eccentricity adjustment, so do not do
          !    this here:
          dyofyr = 0
          adjes  = 1.0 



          !--- CALCULATE CLOUD OPTICAL PROPERTIES:
          if(inflgsw .gt. 0) then 

             do k = kts, kte
                cldfrac(ncol,k) = cldfra1d(k)
             enddo

             !--- zero out cloud optical properties here (not used when passing physical properties to rrtmg_sw):
             do k = kts, kte
                do nb = 1, nbndsw
                   taucld(nb,ncol,k) = 0.0
                   ssacld(nb,ncol,k) = 1.0
                   asmcld(nb,ncol,k) = 0.0
                   fsfcld(nb,ncol,k) = 0.0
                enddo
             enddo

             pcols = ncol
             pver  = kte - kts + 1
             gravmks = g

             if(has_reqc .ne. 0) then
                !--- fill the effective radius for cloud water with that calculated in the Thompson cloud
                !    microphysics parameterization:
                inflgsw = 3
                do k = kts, kte
                   reliq(ncol,k) = max(2.5,re_cloud(i,k,j)*1.e6)
                   if(reliq(ncol,k).le.2.5 .and. cldfrac(ncol,k).gt.0. .and. &
                      (xland(i,j)-1.5).gt.0.) then        !--- ocean.
                      reliq(ncol,k) = 10.5
                   elseif(reliq(ncol,k).le.2.5 .and. cldfrac(ncol,k).gt.0. .and. &
                      (xland(i,j)-1.5).lt.0.) then        !--- land.
                      reliq(ncol,k) = 7.5
                   endif
                enddo

             else

                !---  calculate the effective radius for cloud water in layers below the model top. we added the
                !     dimensions pcols, pver so that the calls to subroutines relcalc and reicalc remain the same
                !     as in WRF although the two dimensions are duplicate and not needed inside the subroutines:
                landfrac(ncol) = 2.-xland(i,j)
                landm(ncol)    = landfrac(ncol)
                snowh(ncol)    = 0.001*snow(i,j)
                icefrac(ncol)  = xice(i,j)

                call relcalc(ncol,pcols,pver,tlay,landfrac,landm,icefrac,reliq,snowh)
             endif

             if(has_reqi .ne. 0) then
                !--- fill the effective radius for cloud ice with that calculated in the Thompson cloud
                !    microphysics parameterization:
                inflgsw  = 4
                iceflgsw = 4
                do k = kts, kte
                   reice(ncol,k) = max(5.,re_ice(i,k,j)*1.e6)
                   if(reice(ncol,k).le.5. .and. cldfrac(ncol,k).gt.0.) then
                      idx_rei = int(t3d(i,k,j)-179.)
                      idx_rei = min(max(idx_rei,1),75)
                      corr = t3d(i,k,j) - int(t3d(i,k,j))
                      reice(ncol,k) = retab(idx_rei)*(1.-corr) + retab(idx_rei+1)*corr
                      reice(ncol,k) = max(reice(ncol,k),5.0)
                   endif
                enddo

             else

                !--- calculate the effective radius for cloud ice in layers below the model top, following
                !    Kristjansson and Mitchell. limit upper bound of reice for Fu ice parameterization and
                !    convert effective radius to generalized effective size (*1.0315; Fu, 1996):

                call reicalc(ncol,pcols,pver,tlay,reice)

                if(iceflgsw .eq. 3) then
                   do k = kts, kte
                      reice(ncol,k) = reice(ncol,k) * 1.0315
                      reice(ncol,k) = min(140.0,reice(ncol,k))
                   enddo
                endif
             endif

             if(has_reqs .ne. 0) then
                !--- fill the effective radius for snow with that calculated in the Thompson cloud
                !    microphysics parameterization:
                inflgsw  = 5
                iceflgsw = 5
                do k = kts, kte
                   resnow(ncol,k) = max(10.,re_snow(i,k,j)*1.e6)
                enddo
             else
                do k = kts, kte
                   resnow(ncol,k) = 10.
                enddo
             endif

             !--- calculate the cloud liquid water path in units of g/m2 below the model top:
             do k = kts, kte
                gliqwp = qc1d(k) * pdel(ncol,k)*100.0 / gravmks * 1000.0              !grid box liquid water path.
                cliqwp(ncol,k) = gliqwp / max(0.01,cldfrac(ncol,k))                   !in-cloud liquid water path.
             enddo

             !--- calculate the cloud ice path in units of g/m2 below the model top:
             if(iceflgsw .eq. 3)then
                do k = kts, kte
                   gicewp = (qi1d(k)+qs1d(k)) * pdel(ncol,k)*100.0 / gravmks * 1000.0 !grid box ice water path.
                   cicewp(ncol,k) = gicewp / max(0.01,cldfrac(ncol,k))                !in-cloud ice water path.
                enddo
             elseif(iceflgsw .ge. 4) then
                do k = kts, kte
                   gicewp = qi1d(k) * pdel(ncol,k)*100.0 / gravmks * 1000.0           !grid box ice water path.
                   cicewp(ncol,k) = gicewp / max(0.01,cldfrac(ncol,k))                !in-cloud ice water path.
                enddo
             endif

             !--- calculate the snow path in units of g/m2 below the model top:
             if(iceflgsw.eq.5)then
                do k = kts, kte
                   snow_mass_factor = 1.0
                   if(resnow(ncol,k) .gt. 130.)then
                      snow_mass_factor = (130.0/resnow(ncol,k))*(130.0/resnow(ncol,k))
                      resnow(ncol,k)   = 130.0
                   endif
                   gsnowp = qs1d(k) * pdel(ncol,k)*100.0 / gravmks * 1000.0           !grid box snow path.
                   csnowp(ncol,k) = snow_mass_factor * gsnowp / max(0.01,cldfrac(ncol,k))
                enddo
             endif

             !--- set the cloud physical property between the model top and the top of the atmosphere.
             !    do not allow clouds between the model top and the top of the atmosphere: 
             do k = kts, kte
                clwpth(ncol,k) = cliqwp(ncol,k)
                ciwpth(ncol,k) = cicewp(ncol,k)
                cswpth(ncol,k) = csnowp(ncol,k)
                rel(ncol,k) = reliq(ncol,k)
                rei(ncol,k) = reice(ncol,k)
                res(ncol,k) = resnow(ncol,k)
             enddo
             !---- no clouds are allowed in the extra layer from model top to TOA:
             clwpth(ncol,kte+1)  = 0.
             ciwpth(ncol,kte+1)  = 0.
             cswpth(ncol,kte+1)  = 0.
             rel(ncol,kte+1)     = 10.
             rei(ncol,kte+1)     = 10.
             res(ncol,kte+1)     = 10.
             cldfrac(ncol,kte+1) = 0.
             do nb = 1, nbndsw
                taucld(nb,ncol,kte+1) = 0.
                ssacld(nb,ncol,kte+1) = 1.
                asmcld(nb,ncol,kte+1) = 0.
                fsfcld(nb,ncol,kte+1) = 0.
             enddo

          endif

          !--- sub-column generator for McICA:
          iplon = 1
          irng  = 0
          permuteseed = 1

          call mcica_subcol_sw &
                     (iplon   , ncol    , nlay    , icld    , permuteseed , irng    , play    , &
                      hgt     , cldfrac , ciwpth  , clwpth  , cswpth      , rei     , rel     , &
                      res     , taucld  , ssacld  , asmcld  , fsfcld      , cldfmcl , ciwpmcl , &
                      clwpmcl , cswpmcl , reicmcl , relqmcl , resnmcl     , taucmcl , ssacmcl , &
                      asmcmcl , fsfcmcl)

          !--- initialization of aerosol optical properties:
          do nb = 1, nbndsw
             do k = kts, kte+1
                tauaer(ncol,k,nb) = 0.
                ssaaer(ncol,k,nb) = 1.
                asmaer(ncol,k,nb) = 0.
             enddo
          enddo

          do na = 1, naerec
             do k = kts, kte+1
                ecaer(ncol,k,na) = 0.
             enddo
          enddo

          !CPF development:
          if(aer_ra_feedback .eq. 1) then
             do nb = 1, nbndsw
                wavemid(nb)=0.5*(wavemin(nb)+wavemax(nb)) ! um

                !convert optical properties at 300,400,600,and 999 to conform to band-wavelengths:
                !tauaer - use angstrom exponent.
                !ssaer - linear interpolation/extrapolation.
                !g - linear interpolation/extrapolation.
                do k = kts, kte
                   if(tauaer300(i,k,j).gt.thresh .and. tauaer999(i,k,j).gt.thresh) then
                      ang = alog(tauaer300(i,k,j)/tauaer999(i,k,j))/alog(999./300.)
                      tauaer(ncol,k,nb)=tauaer400(i,k,j)*(0.4/wavemid(nb))**ang

                      slope = (waer600(i,k,j)-waer400(i,k,j))/.2
                      ssaaer(ncol,k,nb) = slope*(wavemid(nb)-.6)+waer600(i,k,j)
                      if(ssaaer(ncol,k,nb).lt.0.4) ssaaer(ncol,k,nb) = 0.4
                      if(ssaaer(ncol,k,nb).ge.1.0) ssaaer(ncol,k,nb) = 1.0

                      slope=(gaer600(i,k,j)-gaer400(i,k,j))/.2
                      asmaer(ncol,k,nb) = slope*(wavemid(nb)-.6)+gaer600(i,k,j)
                      if(asmaer(ncol,k,nb).lt.0.5) asmaer(ncol,k,nb) = 0.5
                      if(asmaer(ncol,k,nb).ge.1.0) asmaer(ncol,k,nb) = 1.0
                   endif
                enddo
             enddo

             do nb = 1, nbndsw
                slope = 0.  !use slope as a sum holder
                do k = kts, kte
                   slope = slope + tauaer(ncol,k,nb)
                enddo
                if(slope < 0.) then
                   write(msg,'("ERROR: negative total optical depth of ",f8.2," at point i,j,nb=",3i5)') slope,i,j,nb
                   FATAL_ERROR(msg)
                elseif(slope > 6.) then
                   WRITE_MESSAGE("-------------------------")
                   write(msg,'("WARNING: Large total sw optical depth of ",f8.2," at point i,j,nb=",3i5)') slope,i,j,nb
                   WRITE_MESSAGE(msg)
                   WRITE_MESSAGE("Diagnostics 1: k,tauaer300,tauaer400,tauaer600,tauaer999,tauaer")
                   do k = kts, kte
                      write(msg,'(i4,5f8.2)') k,tauaer300(i,k,j),tauaer400(i,k,j), &
                            tauaer600(i,k,j),tauaer999(i,k,j),tauaer(ncol,k,nb)
                      WRITE_MESSAGE(msg)
                      tauaer(ncol,k,nb)=tauaer(ncol,k,nb)*6.0/slope ! up-limit to avoid segmentation fault from AOD (czhao).
                   enddo
                   WRITE_MESSAGE("Diagnostics 2: k,gaer300,gaer400,gaer600,gaer999")
                   do k = kts, kte
                      write(msg,'(i4,4f8.2)') k,gaer300(i,k,j),gaer400(i,k,j), &
                            gaer600(i,k,j),gaer999(i,k,j)
                      WRITE_MESSAGE(msg)
                   enddo
                   WRITE_MESSAGE("Diagnostics 3: k,waer300,waer400,waer600,waer999")
                   do k = kts, kte
                      write(msg,'(i4,4f8.2)') k,waer300(i,k,j),waer400(i,k,j), &
                            waer600(i,k,j),waer999(i,k,j)
                      WRITE_MESSAGE(msg)
                   enddo
                   WRITE_MESSAGE("Diagnostics 4: k,ssaal,asyal,taual")
                   do k = kts-1, kte
                      write(msg,'(i4,3f8.2)') k, ssaaer(i,k,nb),asmaer(i,k,nb),tauaer(i,k,nb)
                      WRITE_MESSAGE(msg)
                   enddo
                   WRITE_MESSAGE("-------------------------")
                endif
             enddo
          endif
          !CPF development end.



          !--- CALL TO THE RRTMG SHORT WAVE RADIATION MODEL:
          call rrtmg_sw &
                     (ncol      , nlay      , icld      , play      , plev      , tlay      , &
                      tlev      , tsfc      , h2ovmr    , o3vmr     , co2vmr    , ch4vmr    , &
                      n2ovmr    , o2vmr     , asdir     , asdif     , aldir     , aldif     , &
                      coszen    , adjes     , dyofyr    , scon      , inflgsw   , iceflgsw  , &
                      liqflgsw  , cldfmcl   , taucmcl   , ssacmcl   , asmcmcl   , fsfcmcl   , &
                      ciwpmcl   , clwpmcl   , cswpmcl   , reicmcl   , relqmcl   , resnmcl   , &
                      tauaer    , ssaaer    , asmaer    , ecaer     , swuflx    , swdflx    , &
                      swhr      , swuflxc   , swdflxc   , swhrc     , swuflxcln , swdflxcln , &
                      aer_opt   , sibvisdir , sibvisdif , sibnirdir , sibnirdif , swdkdir   , &
                      swdkdif   , swdkdirc  , calc_clean_atm_diag)



          !--- OUTPUTS:
          gsw(i,j)  = swdflx(1,1) - swuflx(1,1)
          swcf(i,j) = (swdflx(1,kte+2)-swuflx(1,kte+2)) - (swdflxc(1,kte+2)-swuflxc(1,kte+2))

          if(present(swupt)) then 
             !output up and down toa fluxes for total and clear sky:
             swupt(i,j)  = swuflx(1,kte+2)
             swuptc(i,j) = swuflxc(1,kte+2)
             swdnt(i,j)  = swdflx(1,kte+2)
             swdntc(i,j) = swdflxc(1,kte+2)
             !output up and down surface fluxes for total and clear sky:
             swupb(i,j)  = swuflx(1,1)
             swupbc(i,j) = swuflxc(1,1)
             swdnb(i,j)  = swdflx(1,1)
             swdnbc(i,j) = swdflxc(1,1)

             !CPF development:
             if(calc_clean_atm_diag .gt. 0)then
                swuptcln(i,j) = swuflxcln(1,kte+2)
                swdntcln(i,j) = swdflxcln(1,kte+2)
                swupbcln(i,j) = swuflxcln(1,1)
                swdnbcln(i,j) = swdflxcln(1,1)
             endif
             !added by Zhenxin for 4 components of swdown radiation:
             if(sf_surface_physics .eq. 8) then
                swvisdir(i,j)  = sibvisdir(1,1)
                swvisdif(i,j)  = sibvisdif(1,1)
                swnirdir(i,j)  = sibnirdir(1,1)
                swnirdif(i,j)  = sibnirdif(1,1)
             endif
             !CPF development end.
          endif
          !CPF development:
          if(present(swddir) .and. present(swddif) .and. present(swddni) .and. &
             present(swddirc) .and. present(swddnic) .and. present(swdownc)) then
             swddir(i,j)   = swdkdir(1,1)          ! jararias 2013/08/10
             swddni(i,j)   = swddir(i,j) / coszrs  ! jararias 2013/08/10
             swddif(i,j)   = swdkdif(1,1)          ! jararias 2013/08/10
             swdownc(i, j) = swdflxc(1,1)          ! PAJ: clear-sky GHI
             swddirc(i,j)  = swdkdirc(1,1)         ! PAJ: clear-sky direct normal irradiance
             swddnic(i,j)  = swddirc(i,j) / coszrs ! PAJ: clear-sky direct normal irradiance
          endif
          !CPF development end.

          if(present(swupflx)) then
             do k = kts, kte+2
                swupflx(i,k,j)  = swuflx(1,k)
                swupflxc(i,k,j) = swuflxc(1,k)
                swdnflx(i,k,j)  = swdflx(1,k)
                swdnflxc(i,k,j) = swdflxc(1,k)
             enddo
          endif

          !--- output heating rate tendency; convert heating rate from K/d to K/s:
          do k = kts, kte 
             tten1d(k) = swhr(ncol,k)/86400.
             rthratensw(i,k,j) = tten1d(k)/pi3d(i,k,j)
          enddo

       else

          if(present(swupt)) then 
             !output up and down toa fluxes for total and clear sky:
             swupt(i,j)     = 0.
             swuptc(i,j)    = 0.
             swdnt(i,j)     = 0.
             swdntc(i,j)    = 0.
             !output up and down surface fluxes for total and clear sky:
             swupb(i,j)     = 0.
             swupbc(i,j)    = 0.
             swdnb(i,j)     = 0.
             swdnbc(i,j)    = 0.

             !CPF development:
             if(calc_clean_atm_diag .gt. 0)then
                swuptcln(i,j) = 0.
                swdntcln(i,j) = 0.
                swupbcln(i,j) = 0.
                swdnbcln(i,j) = 0.
             endif
             !added by Zhenxin for 4 components of swdown radiation:
             if(sf_surface_physics .eq. 8) then
                swvisdir(i,j)  = 0.
                swvisdif(i,j)  = 0.
                swnirdir(i,j)  = 0.
                swnirdif(i,j)  = 0.
             endif
             !CPF development end.
          endif
          !CPF development:
          if(present(swddir) .and. present(swddif) .and. present(swddni) .and. &
             present(swddirc) .and. present(swddnic) .and. present(swdownc)) then
             swddir(i,j)   = 0. ! jararias 2013/08/10
             swddni(i,j)   = 0. ! jararias 2013/08/10
             swddif(i,j)   = 0. ! jararias 2013/08/10
             swdownc(i, j) = 0. ! PAJ: clear-sky GHI
             swddirc(i,j)  = 0. ! PAJ: clear-sky direct normal irradiance
             swddnic(i,j)  = 0. ! PAJ: clear-sky direct normal irradiance
          endif
          !CPF development end.

       endif

    end do i_loop !end longitude loop.

 end do j_loop !end latitude loop.                                           

 end subroutine rrtmg_swrad

!=================================================================================================================
 end module module_ra_rrtmg_sw
!=================================================================================================================
