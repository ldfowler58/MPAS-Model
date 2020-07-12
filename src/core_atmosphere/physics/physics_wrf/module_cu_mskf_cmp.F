!
!module  mskf_microphysics
! Adapted to WRF3.8 by Kiran Alapaty
!  !ckay = !dkay = Kiran Alapaty, EPA
! PSH : Sep 2015: copuled with CESM global climatological aerosol data
! TWG : Jun 2016: porting to WRFV3.8 
! TWG & cKAY: Feb 2017: replaced sheet cloud microphysics with that of cumulus clouds

! Purpose:
!!!!#define WRF_PORT
!   CAM Interface for cumulus microphysics
!
! Initial Authors: Xiaoliang Song and Guang Jun Zhang, June 2010
! MSKF adaptation authors: Kiran Alapaty, Tim Glotfelty, and Patrick Hawbecker
!
! Adapted to MSKF scheme: Kiran Alapaty at EPA March 2013 (WRF version)
!

!=================================================================================================================
 module module_cu_mskf_cmp

#if defined(mpas)
use mpas_kind_types,only: r8=>R8KIND
use mpas_atmphys_utilities, only: physics_message,physics_error_fatal
#define FATAL_ERROR(M) call physics_error_fatal( M )
#define WRITE_MESSAGE(M) call physics_message( M )
#else
use shr_kind_mod,only: r8=>shr_kind_r8
use module_wrf_error
#define FATAL_ERROR(M) call wrf_error_fatal( M )
#define WRITE_MESSAGE(M) call wrf_message( M )
#endif

 use error_function, only: erf,erfc

 implicit none
 private
 public:: mskf_mphyi, mskf_mphy, mskf_GAMMA, mskf_polysvp

! Private module data

  integer, parameter :: naer_cu = 10        
  integer, parameter :: pcols = 1

!constants remaped
  real(r8), private::  g              !gravity
  real(r8), private::  mw             !molecular weight of water
  real(r8), private::   r        !Dry air Gas constant
  real(r8), private::   rv       !water vapor gas contstant
  real(r8), private::   rr   !universal gas constant
  real(r8), private::   cpp                  !specific heat of dry air
  real(r8), private::   rhow               !density of liquid water
  real(r8), private::  xlf !latent heat of freezing

!from physconst 
  real(r8), private, parameter ::  gravit = 9.80616_r8      ! acceleration of gravity ~ m/s^2
  real(r8), private, parameter ::  rair   = 287.04239_r8    ! Dry air gas constant     ~ J/K/kg 
  real(r8), private, parameter ::  tmelt  = 273.15_r8       ! freezing T of fresh water  ~ K
  real(r8), private, parameter ::  cpair  = 1.00464e3_r8    ! specific heat of dry air   ~ J/kg/K
  real(r8), private, parameter ::  rh2o   = 461.915_r8      ! Water vapor gas constant ~ J/K/kg
  real(r8), private, parameter ::  r_universal = 8.31447e3_r8  ! Universal gas constant ~ J/K/kmole
  real(r8), private, parameter ::  mwh2o  = 18._r8          ! molecular weight h2o
  real(r8), private, parameter ::  rhoh2o = 1.000e3_R8      ! density of fresh water     ~ kg/m^3
  real(r8), private, parameter ::  latvap = 2.501e6_r8      ! latent heat of evaporation ~ J/kg
  real(r8), private, parameter ::  latice = 3.337e5_r8      ! latent heat of fusion      ~ J/kg
  real(r8), private, parameter ::  epsilo = 0.622_r8        ! ratio of h2o to dry air molecular weights         

!from 'microconstants'
  real(r8), private:: rhosn  ! bulk density snow
  real(r8), private:: rhoi   ! bulk density ice

  real(r8), private:: ac,bc,as,bs,ai,bi,ar,br  !fall speed parameters 
  real(r8), private:: ci,di    !ice mass-diameter relation parameters
  real(r8), private:: cs,ds    !snow mass-diameter relation parameters
  real(r8), private:: cr,dr    !drop mass-diameter relation parameters
  real(r8), private:: Eii      !collection efficiency aggregation of ice
  real(r8), private:: Ecc      !collection efficiency
  real(r8), private:: Ecr      !collection efficiency cloud droplets/rain
  real(r8), private:: DCS      !autoconversion size threshold
  real(r8), private:: F14      !Ferrier (1994) Time scale parameter
  real(r8), private:: qsmall   !min mixing ratio 
  real(r8), private:: bimm,aimm !immersion freezing
  real(r8), private:: rhosu     !typical 850mn air density
  real(r8), private:: mi0       ! new crystal mass
  real(r8), private:: rin       ! radius of contact nuclei
  real(r8), private:: pi       ! pi

  real(r8), private:: rn_dst1, rn_dst2, rn_dst3, rn_dst4  !dust number mean radius for contact freezing
!..........................................................................

!needed for findsp
real(r8), private:: t0       ! Freezing temperature

! activate parameters

      integer, private:: psat
      parameter (psat=6) ! number of supersaturations to calc ccn concentration
      real(r8), private:: aten
!      
      real(r8), private:: alogsig(naer_cu) ! natl log of geometric standard dev of aerosol
      real(r8), private:: exp45logsig(naer_cu)
      real(r8), private:: argfactor(naer_cu)
      real(r8), private:: amcube(naer_cu) ! cube of dry mode radius (m)
      real(r8), private:: smcrit(naer_cu) ! critical supersatuation for activation
      real(r8), private:: lnsm(naer_cu) ! ln(smcrit)
      real(r8), private:: amcubesulfate(pcols) ! cube of dry mode radius (m) of sulfate
      real(r8), private:: smcritsulfate(pcols) ! critical supersatuation for activation of sulfate
      real(r8), private:: amcubefactor(naer_cu) ! factors for calculating mode radius
      real(r8), private:: smcritfactor(naer_cu) ! factors for calculating critical supersaturation
      real(r8), private:: super(psat)
      real(r8), private:: alogten,alog2,alog3,alogaten
      real(r8), private, parameter :: supersat(psat)= &! supersaturation (%) to determine ccn concentration
               (/0.02,0.05,0.1,0.2,0.5,1.0/)
      real(r8), private:: ccnfact(psat,naer_cu)

      real(r8), private:: f1(naer_cu),f2(naer_cu) ! abdul-razzak functions of width
      real(r8), private:: third, sixth,zero
      real(r8), private:: sq2, sqpi


!wrf      integer :: naer_all    ! number of aerosols affecting climate
!xsong 2013-08-22---------------
      integer :: idxsul = 1 ! index in aerosol list for sulfate  
      integer :: idxdst1 = 3 ! index in aerosol list for dust1
      integer :: idxdst2 = 4 ! index in aerosol list for dust2
      integer :: idxdst3 = 5 ! index in aerosol list for dust3
      integer :: idxdst4 = 6 ! index in aerosol list for dust4
      integer :: idxbcphi = 10 ! index in aerosol list for Soot (BCPHI)
!xsong 2013-08-22---------------
      ! aerosol properties
      character(len=20)  aername(naer_cu)
      real(r8) dryrad_aer(naer_cu)
      real(r8) density_aer(naer_cu)
      real(r8) hygro_aer(naer_cu)
      real(r8) dispersion_aer(naer_cu)
      real(r8) num_to_mass_aer(naer_cu)

!xsong 2013-08-22--------------------
   data aername /"SULFATE","SEASALT2","DUST1","DUST2","DUST3","DUST4","OCPHO","BCPHO",   &
                 "OCPHI","BCPHI"/
   data dryrad_aer /0.695E-7_r8,0.200E-5_r8,0.151E-5_r8,0.151E-5_r8,0.151E-5_r8,0.151E-5_r8,     &
                    0.212E-7_r8,0.118E-7_r8,0.212E-7_r8, 0.118E-7_r8/
   data density_aer /1770._r8,2200._r8,2600._r8,2600._r8,2600._r8,2600._r8,1800._r8,  &
                     1000._r8,2600._r8,1000._r8/
   data hygro_aer /0.507_r8,1.160_r8,0.140_r8,0.140_r8,0.140_r8,0.140_r8,0.100_r8,0.100_r8,  &
                   0.140_r8,0.100_r8/
   data dispersion_aer /2.030_r8,1.3732_r8,1.900_r8,1.900_r8,1.900_r8,1.900_r8,2.240_r8,  &
                        2.000_r8,2.240_r8,2.000_r8/
   data num_to_mass_aer /42097098109277080._r8,8626504211623._r8,3484000000000000._r8,213800000000000._r8,&
                         22050000000000._r8,3165000000000._r8,0.745645E+18_r8,0.167226E+20_r8,&
                         0.516216E+18_r8,0.167226E+20_r8/
!xsong 2013-08-22-----------------------


 contains


!=================================================================================================================
 subroutine mskf_mphyi

!----------------------------------------------------------------------- 
! 
! Purpose:
! initialize constants for the cumulus microphysics
! called from zm_conv_init() in zm_conv_intr.F90
!
! Author: Xiaoliang Song, June 2010
! 
!-----------------------------------------------------------------------

!     save    ! sep6
!wrf   use pmgrid, only: plev, plevp
      integer k

      integer l,m, iaer
      real(r8) surften       ! surface tension of water w/respect to air (N/m)
!      real(r8) arg

! hm modify to use my error function


!declarations for morrison codes (transforms variable names)

!   g= gravit                  !gravity
!   mw = mwh2o / 1000._r8      !molecular weight of water
!   r= rair                   !Dry air Gas constant: note units(phys_constants
!   are in J/K/kmol)
!   rv= rh2o                   !water vapor gas contstant
!   rr = r_universal           !universal gas constant
!   cpp = cpair                 !specific heat of dry air
!   rhow = rhoh2o              !density of liquid water

!NOTE:
! latent heats should probably be fixed with temperature 
! for energy conservation with the rest of the model
! (this looks like a +/- 3 or 4% effect, but will mess up energy balance)

   xlf = latice          ! latent heat freezing


! from microconstants

! parameters below from Reisner et al. (1998)
! density parameters (kg/m3)

      rhosn = 100._r8    ! bulk density snow
      rhoi = 500._r8     ! bulk density ice
      rhow = 1000._r8    ! bulk density liquid

! fall speed parameters, V = aD^b
! V is in m/s

! droplets
        ac = 3.e7_r8
        bc = 2._r8

! snow
        as = 11.72_r8
        bs = 0.41_r8

! cloud ice
        ai = 700._r8
        bi = 1._r8

! rain
        ar = 841.99667_r8
        br = 0.8_r8

! particle mass-diameter relationship
! currently we assume spherical particles for cloud ice/snow
! m = cD^d

        pi= 3.1415927_r8

! cloud ice mass-diameter relationship

        ci = rhoi*pi/6._r8
        di = 3._r8

! snow mass-diameter relationship

        cs = rhosn*pi/6._r8
        ds = 3._r8

! drop mass-diameter relationship

        cr = rhow*pi/6._r8
        dr = 3._r8

! collection efficiency, aggregation of cloud ice and snow

        Eii = 0.1_r8

! collection efficiency, accretion of cloud water by rain

        Ecr = 1.0_r8

! autoconversion size threshold for cloud ice to snow (m)

!          Dcs = 100.e-6_r8
        Dcs = 200.e-6_r8

! Ferrier [1994] time period parameter ! TWG Feb17
         F14 = 100.0 !180.0 Original
         

! smallest mixing ratio considered in microphysics

        qsmall = 1.e-28_r8 !Shaocai   !1.e-18_r8  

! immersion freezing parameters, bigg 1953

        bimm = 100._r8
        aimm = 0.66_r8

! contact freezing due to dust
! dust number mean radius (m), Zender et al JGR 2003 assuming number mode radius
! of 0.6 micron, sigma=2

        rn_dst1=0.258e-6_r8
        rn_dst2=0.717e-6_r8
        rn_dst3=1.576e-6_r8
        rn_dst4=3.026e-6_r8

! typical air density at 850 mb

        rhosu = 85000._r8/(rair * tmelt)

! mass of new crystal due to aerosol freezing and growth (kg)

        mi0 = 4._r8/3._r8*pi*rhoi*(10.e-6_r8)*(10.e-6_r8)*(10.e-6_r8)

! radius of contact nuclei aerosol (m)

        rin = 0.1e-6_r8

! freezing temperature
        t0=273.15_r8

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! set parameters for droplet activation, following abdul-razzak and ghan 2000,
! JGR

!      mathematical constants

      zero=0._r8
      third=1./3._r8
      sixth=1./6._r8
      sq2=sqrt(2._r8)
      pi=4._r8*atan(1.0_r8)
      sqpi=sqrt(pi)

      surften=0.076_r8
      aten=2.*mwh2o*surften/(r_universal*t0*rhoh2o)
      alogaten=log(aten)
      alog2=log(2._r8)
      alog3=log(3._r8)
      super(:)=0.01*supersat(:)

      do m=1,naer_cu
!         use only if width of size distribution is prescribed
          alogsig(m)=log(dispersion_aer(m))
          exp45logsig(m)=exp(4.5*alogsig(m)*alogsig(m))
          argfactor(m)=2./(3.*sqrt(2.)*alogsig(m))
          f1(m)=0.5*exp(2.5*alogsig(m)*alogsig(m))
          f2(m)=1.+0.25*alogsig(m)
          amcubefactor(m)=3._r8/(4._r8*pi*exp45logsig(m)*density_aer(m))
          smcritfactor(m)=2._r8*aten*sqrt(aten/(27._r8*max(1.e-10_r8,hygro_aer(m))))
!         use only if mode radius of size distribution is prescribed
          amcube(m)=amcubefactor(m)/(num_to_mass_aer(m))
!         use only if only one component per mode
          if(hygro_aer(m).gt.1.e-10) then
             smcrit(m)=smcritfactor(m)/sqrt(amcube(m))
          else
             smcrit(m)=100.
          endif
          lnsm(m)=log(smcrit(m))

      end do

   return
 end subroutine mskf_mphyi

!===============================================================================

subroutine mskf_mphy(su,    qu,   mu,   du,   cmel, cmei, zf,  pm,  te,   qe, eps0,    &
                   jb,    jt,   jlcl, msg,  il2g, grav, cp,  rd,  qc,   qi, qr, qni, & ! TWG
                    rprd,  wu,    eu,   nc,   ni, nr, ns, dum2l, sprd, frz, aer_mmr, deltat, & !TWG
                   Pver,PverP,gamhat,qsatzm,wu_mskf_act,qc_mskf_act,qi_mskf_act,effc,effi,effs)

! Purpose:
! microphysic parameterization for Zhang-McFarlane convection scheme
! called from cldprp() in zm_conv.F90
!
! Author: Xiaoliang Song, June 2010
!
! Adaptation: Kiran Alapaty at EPA 2013 for MSKF convection scheme in WRF

!wrf  use time_manager,    only: get_nstep, get_step_size

! variable declarations

  implicit none

! input variables

  integer, parameter :: naer_cu = 10   
  integer, parameter :: pcols = 1

  integer :: pver                  ! number of vertical levels(mid-layer)
  integer :: pverp                 ! number of vertical levels(interface)      
  real(r8) :: su(pcols,pver)        ! normalized dry stat energy of updraft
  real(r8) :: qu(pcols,pver)        ! spec hum of updraft
  real(r8) :: mu(pcols,pver)        ! updraft mass flux
  real(r8) :: du(pcols,pver)        ! detrainement rate of updraft
  real(r8) :: cmel(pcols,pver)      ! condensation rate of updraft
  real(r8) :: cmei(pcols,pver)      ! condensation rate of updraft
  real(r8) :: zf(pcols,pverp)       ! height of interfaces
  real(r8) :: pm(pcols,pver)        ! pressure of env
  real(r8) :: te(pcols,pver)        ! temp of env
  real(r8) :: qe(pcols,pver)        ! spec. humidity of env
  real(r8) :: eps0(pcols)
  real(r8) :: eu(pcols,pver)        ! entrainment rate of updraft
!ckay  real(r8) :: aer_mmr(:,:,:)        ! aerosol mass mixing ratio
  real(r8) :: aer_mmr(Pcols,Pver,naer_cu)        ! aerosol mass mixing ratio
!  real(r8) :: gamhat(pcols,pver)    ! kf_GAMMA=L/cp(dq*/dT) at interface
!ckay
  real(r8) :: qsatzm(pcols,pver)        ! spec hum of updraft
  real(r8) :: wu_mskf_act(pver)        ! KF incloud updraft velocity
  real(r8) :: qc_mskf_act(pver)        ! KF incloud updraft velocity
  real(r8) :: qi_mskf_act(pver)        ! KF incloud updraft velocity

  integer :: jb(pcols)              ! updraft base level
  integer :: jt(pcols)              ! updraft plume top
  integer :: jlcl(pcols)            ! updraft lifting cond level
  integer :: msg                    ! missing moisture vals
  integer :: il2g                   ! CORE GROUP REMOVE

  real(r8) grav                                 ! gravity
  real(r8) cp                                   ! heat capacity of dry air
  real(r8) rd                                   ! gas constant for dry air

! output variables
!ckay
  real(r8) qc(pcols,pver)       ! cloud water mixing ratio (kg/kg)
  real(r8) qi(pcols,pver)       ! cloud ice mixing ratio (kg/kg)
  real(r8) nc(pcols,pver)       ! cloud water number conc (1/kg)
  real(r8) ni(pcols,pver)       ! cloud ice number conc (1/kg)
  real(r8)  qni(pcols,pver)      ! snow mixing ratio
  real(r8)  qr(pcols,pver)       ! rain mixing ratio
  real(r8)  ns(pcols,pver)       ! snow number conc
  real(r8)  nr(pcols,pver)       ! rain number conc
  real(r8) rprd(pcols,pver)     ! rate of production of precip at that layer
!ckay  real(r8), intent(out) :: rprd(pcols,pver)     ! rate of production of
!precip at that layer
  real(r8) sprd(pcols,pver)     ! rate of production of snow at that layer
  real(r8) frz(pcols,pver)      ! rate of freezing

! tendency for output

  real(r8) :: autolm(pcols,pver)    !mass tendency due to autoconversion of droplets to rain
  real(r8) :: accrlm(pcols,pver)    !mass tendency due to accretion of droplets by rain
  real(r8) :: bergnm(pcols,pver)    !mass tendency due to Bergeron process
  real(r8) :: fhtimm(pcols,pver)    !mass tendency due to immersion freezing
  real(r8) :: fhtctm(pcols,pver)    !mass tendency due to contact freezing
  real(r8) :: fhmlm (pcols,pver)    !mass tendency due to homogeneous freezing
  real(r8) :: hmpim (pcols,pver)    !mass tendency due to HM process
  real(r8) :: accslm(pcols,pver)    !mass tendency due to accretion of droplets by snow
  real(r8) :: dlfm  (pcols,pver)    !mass tendency due to detrainment of droplet
  real(r8) :: autoln(pcols,pver)    !num tendency due to autoconversion of droplets to rain
  real(r8) :: accrln(pcols,pver)    !num tendency due to accretion of droplets by rain
  real(r8) :: bergnn(pcols,pver)    !num tendency due to Bergeron process
  real(r8) :: fhtimn(pcols,pver)    !num tendency due to immersion freezing
  real(r8) :: fhtctn(pcols,pver)    !num tendency due to contact freezing
  real(r8) :: fhmln (pcols,pver)    !num tendency due to homogeneous freezing
  real(r8) :: accsln(pcols,pver)    !num tendency due to accretion of droplets by snow
  real(r8) :: activn(pcols,pver)    !num tendency due to droplets activation
  real(r8) :: dlfn  (pcols,pver)    !num tendency due to detrainment of droplet
  real(r8) :: autoim(pcols,pver)    !mass tendency due to autoconversion of cloud ice to snow
  real(r8) :: accsim(pcols,pver)    !mass tendency due to accretion of cloud ice by snow
  real(r8) :: difm  (pcols,pver)    !mass tendency due to detrainment of cloud ice 
  real(r8) :: nuclin(pcols,pver)    !num tendency due to ice nucleation
  real(r8) :: nuclim(pcols,pver)    !mass tendency due to ice nucleation
  real(r8) :: collrm(pcols,pver)    !mass tendency due to rain-ice collection
  real(r8) :: collrn(pcols,pver)    !number tendency due to rain-ce collection
  real(r8) :: fhtcrm(pcols,pver)    !mass tendency due to rain freezing to snow
  real(r8) :: fhtcrn(pcols,pver)    !num tendency due to rain freezing to snow
  real(r8) :: autorn(pcols,pver)    !num tendency for autoconversion of clouds to rain (rain term)
  real(r8) :: aggrn(pcols,pver)     !num tendency for self collection of rain
  real(r8) :: aggsn(pcols,pver)     !num tendency for self collection of snow
  real(r8) :: autoin(pcols,pver)    !num tendency due to autoconversion of cloud ice to snow
  real(r8) :: accsin(pcols,pver)    !num tendency due to accretion of cloud ice by snow
  real(r8) :: hmpin (pcols,pver)    !num tendency due to HM process
  real(r8) :: difn  (pcols,pver)    !num tendency due to detrainment of cloud ice
  real(r8) :: trspcm(pcols,pver)    !LWC tendency due to convective transport
  real(r8) :: trspcn(pcols,pver)    !droplet num tendency due to convective transport
  real(r8) :: trspim(pcols,pver)    !IWC tendency due to convective transport
  real(r8) :: trspin(pcols,pver)    !ice crystal num tendency due to convective transport

  real(r8) :: ncadj(pcols,pver)     !droplet num tendency due to adjustment
  real(r8) :: niadj(pcols,pver)     !ice crystal num tendency due to adjustment
  real(r8) :: qcadj(pcols,pver)     !droplet mass tendency due to adjustment
  real(r8) :: qiadj(pcols,pver)     !ice crystal mass tendency due to adjustment

! output for ice nucleation
  real(r8) :: nimey(pcols,pver)     !output number conc of ice nuclei due to meyers deposition (1/m3)
  real(r8) :: nihf(pcols,pver)      !output number conc of ice nuclei due to heterogenous freezing (1/m3)
  real(r8) :: nidep(pcols,pver)     !output number conc of ice nuclei due to deoposion nucleation (hetero nuc) (1/m3)
  real(r8) :: niimm(pcols,pver)     !output number conc of ice nuclei due to immersion freezing (hetero nuc) (1/m3)
  real(r8) :: effc(pcols,pver)    ! droplet effective radius (micron)
  real(r8) :: effi(pcols,pver)    ! cloud ice effective radius (micron)
  real(r8) :: effs(pcols,pver)    ! snow effective radius (micron)

!................................................................................
! local workspace
! all units mks unless otherwise stated
  real(r8) :: deltat                ! time step (s)
  real(r8) :: omsm                  ! number near unity for round-off issues
  real(r8) :: dum                   ! temporary dummy variable
  real(r8) :: arg                   ! argument of erfc
  real(r8) :: dum1                  ! temporary dummy variable 
  real(r8) :: dum2                  ! temporary dummy variable

  real(r8) :: q(pcols,pver)         ! water vapor mixing ratio (kg/kg)
  real(r8) :: t(pcols,pver)         ! temperature (K)
  real(r8) :: rho(pcols,pver)       ! air density (kg m-3)
  real(r8) :: dz(pcols,pver)        ! height difference across model verticallevel

  real(r8) :: qcic(pcols,pver)      ! in-cloud cloud liquid mixing ratio
  real(r8) :: qiic(pcols,pver)      ! in-cloud cloud ice mixing ratio
!dkay
  real (r8) :: tot_qc_qi
  real(r8) :: qniic(pcols,pver)     ! in-precip snow mixing ratio
  real(r8) :: qric(pcols,pver)      ! in-precip rain mixing ratio
  real(r8) :: ncic(pcols,pver)      ! in-cloud droplet number conc
  real(r8) :: niic(pcols,pver)      ! in-cloud cloud ice number conc
  real(r8) :: nsic(pcols,pver)      ! in-precip snow number conc
  real(r8) :: nric(pcols,pver)      ! in-precip rain number conc

  real(r8) :: lami(pver)            ! slope of cloud ice size distr
  real(r8) :: n0i(pver)             ! intercept of cloud ice size distr
  real(r8) :: lamc(pver)            ! slope of cloud liquid size distr
  real(r8) :: n0c(pver)             ! intercept of cloud liquid size distr
  real(r8) :: lams(pver)            ! slope of snow size distr
  real(r8) :: n0s(pver)             ! intercept of snow size distr
  real(r8) :: lamr(pver)            ! slope of rain size distr
  real(r8) :: n0r(pver)             ! intercept of rain size distr
  real(r8) :: cdist1(pver)          ! size distr parameter to calculate droplet freezing
  real(r8) :: pgam(pver)            ! spectral width parameter of droplet size distr
  real(r8) :: lammax                ! maximum allowed slope of size distr
  real(r8) :: lammin                ! minimum allowed slope of size distr

  real(r8) :: mnuccc(pver)          ! mixing ratio tendency due to freezing of cloud water
  real(r8) :: nnuccc(pver)          ! number conc tendency due to freezing of cloud water
  real(r8) :: mnucct(pver)          ! mixing ratio tendency due to contact freezing of cloud water
  real(r8) :: nnucct(pver)          ! number conc tendency due to contact freezing of cloud water
  real(r8) :: msacwi(pver)          ! mixing ratio tendency due to HM ice multiplication
  real(r8) :: nsacwi(pver)          ! number conc tendency due to HM ice multiplication
  real(r8) :: prf(pver)             ! mixing ratio tendency due to fallout of rain
  real(r8) :: psf(pver)             ! mixing ratio tendency due to fallout of snow
  real(r8) :: pnrf(pver)            ! number conc tendency due to fallout of rain
  real(r8) :: pnsf(pver)            ! number conc tendency due to fallout of snow
  real(r8) :: prc(pver)             ! mixing ratio tendency due to autoconversion of cloud droplets
  real(r8) :: nprc(pver)            ! number conc tendency due to autoconversion of cloud droplets
  real(r8) :: nprc1(pver)           ! qr tendency due to autoconversion of cloud droplets
  real(r8) :: nsagg(pver)           ! ns tendency due to self-aggregation of snow
  real(r8) :: dc0                   ! mean size droplet size distr
  real(r8) :: ds0                   ! mean size snow size distr (area weighted)
  real(r8) :: eci                   ! collection efficiency for riming of snow by droplets
  real(r8) :: dv(pcols,pver)        ! diffusivity of water vapor in air
  real(r8) :: mua(pcols,pver)       ! viscocity of air
  real(r8) :: psacws(pver)          ! mixing rat tendency due to collection of droplets by snow
  real(r8) :: npsacws(pver)         ! number conc tendency due to collection of droplets by snow
  real(r8) :: pracs(pver)           ! mixing rat tendency due to collection of rain by snow
  real(r8) :: npracs(pver)          ! number conc tendency due to collection of rain by snow
  real(r8) :: mnuccr(pver)          ! mixing rat tendency due to freezing of rain
  real(r8) :: nnuccr(pver)          ! number conc tendency due to freezing of rain
  real(r8) :: pra(pver)             ! mixing rat tendnency due to accretion of droplets by rain
  real(r8) :: npra(pver)            ! nc tendnency due to accretion of droplets by rain
  real(r8) :: nragg(pver)           ! nr tendency due to self-collection of rain
  real(r8) :: prci(pver)            ! mixing rat tendency due to autoconversion of cloud ice to snow
  real(r8) :: nprci(pver)           ! number conc tendency due to autoconversion of cloud ice to snow
  real(r8) :: prai(pver)            ! mixing rat tendency due to accretion of cloud ice by snow
  real(r8) :: nprai(pver)           ! number conc tendency due to accretion of cloud ice by snow
  real(r8) :: prb(pver)             ! rain mixing rat tendency due to Bergeron process
  real(r8) :: nprb(pver)            ! number conc tendency due to Bergeron process


! fall speed
  real(r8) :: arn(pcols,pver)       ! air density corrected rain fallspeed
  real(r8) :: asn(pcols,pver)       ! air density corrected snow fallspeed
  real(r8) :: acn(pcols,pver)       ! air density corrected cloud droplet fallspeed parameter
  real(r8) :: ain(pcols,pver)       ! air density corrected cloud ice fallspeed
  real(r8) :: uns(pver)             ! number-weighted snow fallspeed
  real(r8) :: ums(pver)             ! mass-weighted snow fallspeed
  real(r8) :: unr(pver)             ! number-weighted rain fallspeed
  real(r8) :: umr(pver)             ! mass-weighted rain fallspeed

! conservation check
  real(r8) :: qce                   ! dummy qc for conservation check
  real(r8) :: qie                   ! dummy qi for conservation check
  real(r8) :: nce                   ! dummy nc for conservation check
  real(r8) :: nie                   ! dummy ni for conservation check
  real(r8) :: qre                   ! dummy qr for conservation check
  real(r8) :: nre                   ! dummy nr for conservation check
  real(r8) :: qnie                  ! dummy qni for conservation check
  real(r8) :: nse                   ! dummy ns for conservation check      
  real(r8) :: ratio                 ! parameter for conservation check

! sum of source/sink terms for cloud hydrometeor
  real(r8) :: qctend(pcols,pver)    ! microphysical tendency qc (1/s)
  real(r8) :: qitend(pcols,pver)    ! microphysical tendency qi (1/s)
  real(r8) :: nctend(pcols,pver)    ! microphysical tendency nc (1/(kg*s))
  real(r8) :: nitend(pcols,pver)    ! microphysical tendency ni (1/(kg*s))
  real(r8) :: qnitend(pcols,pver)   ! snow mixing ratio source/sink term
  real(r8) :: nstend(pcols,pver)    ! snow number concentration source/sink term
  real(r8) :: qrtend(pcols,pver)    ! rain mixing ratio source/sink term
  real(r8) :: nrtend(pcols,pver)    ! rain number concentration source/sink term

! terms for Bergeron process
  real(r8) :: bergtsf               !bergeron timescale to remove all liquid
  real(r8) :: plevap                ! cloud liquid water evaporation rate

! aerosol variables
  real(r8) :: naermod(naer_cu)      ! aerosol number concentration (/m3)
  real(r8) :: naer2(pcols,pver,naer_cu)    ! new aerosol number concentration (/m3)
  real(r8) :: naer2h(pcols,pver,naer_cu)   ! new aerosol number concentration (/m3) 
  real(r8) :: maerosol(1,naer_cu)   ! aerosol mass conc (kg/m3)
  real(r8) naer(pcols)

! droplet activation
  real(r8) :: dum2l(pcols,pver)     ! number conc of CCN (1/kg)
  real(r8) :: npccn(pver)           ! droplet activation rate
  real(r8) :: ncmax
  real(r8) :: mtimec                ! factor to account for droplet activation timescale

! ice nucleation
  real(r8) :: dum2i(pcols,pver)     ! number conc of ice nuclei available (1/kg)
  real(r8) :: qs(pcols,pver)        ! liquid-ice weighted sat mixing rat (kg/kg)
  real(r8) :: es(pcols,pver)        ! sat vapor press (pa) over water
  real(r8) :: relhum(pcols,pver)    ! relative humidity
  real(r8) :: esi(pcols,pver)       ! sat vapor press (pa) over ice
  real(r8) :: nnuccd(pver)          ! ice nucleation rate from deposition/cond.-freezing
  real(r8) :: mnuccd(pver)          ! mass tendency from ice nucleation
  real(r8) :: nimax
  real(r8) :: mtime                 ! factor to account for ice nucleation timescale
  real(r8) :: gamhat(pcols,pver)    ! kf_GAMMA=L/cp(dq*/dT) at interface


! loop array variables
  integer i,k,nstep,n, l
  integer ii,kk, m

! loop variables for iteration solution
  integer iter,it,ltrue(pcols)

! used in contact freezing via dust particles
  real(r8)  tcnt, viscosity, mfp
  real(r8)  slip1, slip2, slip3, slip4
  real(r8)  dfaer1, dfaer2, dfaer3, dfaer4
  real(r8)  nacon1,nacon2,nacon3,nacon4

! used in immersion freezing via soot
  real(r8) ttend(pver)
  real(r8) naimm
  real(r8) :: ntaer(pcols,pver)
  real(r8) :: ntaerh(pcols,pver)

! used in secondary ice production
  real(r8) ni_secp

! used in vertical velocity calculation
  real(r8) th(pcols,pver)
  real(r8) qh(pcols,pver)
  real(r8) wu(pcols,pver)
  real(r8) zkine(pcols,pver)
  real(r8) zbuo(pcols,pver)
  real(r8) zfacbuo, cwdrag, cwifrac, retv,  zbuoc
  real(r8) zbc, zbe,  zdkbuo, zdken
  real(r8) arcf(pcols,pver)
  real(r8) p(pcols,pver)
  real(r8) ph(pcols,pver)

  real(r8) :: rhoh(pcols,pver)    ! air density (kg m-3) at interface 
  real(r8) :: rhom(pcols,pver)    ! air density (kg m-3) at mid-level
  real(r8) :: tu(pcols,pver)      ! temperature in updraft (K)

  real(r8) :: fhmrm (pcols,pver)  !mass tendency due to homogeneous freezing of rain

  real(r8) ncorg,niorg,qcorg,qiorg

  integer  kqi(pcols),kqc(pcols)
  logical  lcbase(pcols), libase(pcols)

!ckay introduced save sep6
! save

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! initialization
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


! parameters for scheme
        omsm=0.99999_r8
        zfacbuo = 0.5_r8/(1._r8+0.5_r8)
        cwdrag  = 1.875_r8*0.506_r8
        cwifrac = 0.5_r8
        retv    = 0.608_r8
        bergtsf = 1800._r8

! initialize multi-level fields
        do i=1,il2g
          do k=1,pver
            q(i,k)= qu(i,k)
            tu(i,k)= su(i,k) - grav/cp*zf(i,k)
            t(i,k)= su(i,k) - grav/cp*zf(i,k)
            p(i,k) = 100._r8*pm(i,k)
            wu(i,k)  = 0._r8
            zkine(i,k)= 0._r8
            arcf(i,k) = 0._r8
            zbuo(i,k) = 0._r8
            nc(i,k) = 0._r8
            ni(i,k) = 0._r8
            qc(i,k) = 0._r8
            qi(i,k) = 0._r8
            qcic(i,k) = 0._r8
            qiic(i,k) = 0._r8
            ncic(i,k) = nc(i,k)
            niic(i,k) = ni(i,k)
            qr(i,k) = 0._r8
            qni(i,k)= 0._r8
            nr(i,k) = 0._r8
            ns(i,k) = 0._r8
            qric(i,k) = qr(i,k)
            qniic(i,k) = qni(i,k)
            nric(i,k) = nr(i,k)
            nsic(i,k) = ns(i,k)
            nimey(i,k) = 0._r8
            nihf(i,k)  = 0._r8
            nidep(i,k) = 0._r8
            niimm(i,k) = 0._r8

            autolm(i,k) = 0._r8
            accrlm(i,k) = 0._r8
            bergnm(i,k) = 0._r8
            fhtimm(i,k) = 0._r8
            fhtctm(i,k) = 0._r8
            fhmlm (i,k) = 0._r8
            hmpim (i,k) = 0._r8
            accslm(i,k) = 0._r8
            dlfm  (i,k) = 0._r8
            collrm(i,k) = 0._r8
            collrn(i,k) = 0._r8
            fhtcrm(i,k) = 0._r8
            fhtcrn(i,k) = 0._r8
            autorn(i,k) = 0._r8
             aggrn(i,k) = 0._r8
             aggsn(i,k) = 0._r8

            autoln(i,k) = 0._r8

            accrln(i,k) = 0._r8
            bergnn(i,k) = 0._r8
            fhtimn(i,k) = 0._r8
            fhtctn(i,k) = 0._r8
            fhmln (i,k) = 0._r8
            accsln(i,k) = 0._r8
            activn(i,k) = 0._r8
            dlfn  (i,k) = 0._r8
            ncadj (i,k) = 0._r8
            qcadj (i,k) = 0._r8
!cloud ice------------------------
           autoim(i,k) = 0._r8
            accsim(i,k) = 0._r8
            difm  (i,k) = 0._r8
            nuclin(i,k) = 0._r8
            nuclim(i,k) = 0._r8
            autoin(i,k) = 0._r8
            accsin(i,k) = 0._r8
            hmpin (i,k) = 0._r8
            difn  (i,k) = 0._r8
            niadj (i,k) = 0._r8
            qiadj (i,k) = 0._r8

            trspcm(i,k) = 0._r8
            trspcn(i,k) = 0._r8
            trspim(i,k) = 0._r8
            trspin(i,k) = 0._r8

            effc(i,k) = 0._r8
            effi(i,k) = 0._r8
            effs(i,k) = 0._r8

            fhmrm (i,k) = 0._r8
          end do
        end do


! initialize time-varying parameters
        do k=1,pver
          do i=1,il2g
!-------------Shaocai Yu
             if (k .eq.1) then
                rhoh(i,k) = p(i,k)/(t(i,k)*rd)
                rhom(i,k) = p(i,k)/(t(i,k)*rd)
                th (i,k) = te(i,k)
                qh (i,k) = qe(i,k)
                dz (i,k)  = zf(i,k) - zf(i,k+1)
                ph(i,k)   = p(i,k)
             else

               rhoh(i,k) = 0.5_r8*(p(i,k)+p(i,k-1))/(t(i,k)*rd)
                if (k .eq. pver) then
                  rhom(i,k) = p(i,k)/(rd*t(i,k))
                else
                  rhom(i,k) = 2.0_r8*p(i,k)/(rd*(t(i,k)+t(i,k+1)))
                end if
                th (i,k) = 0.5_r8*(te(i,k)+te(i,k-1))
                qh (i,k) = 0.5_r8*(qe(i,k)+qe(i,k-1))
                dz(i,k)  = zf(i,k-1) - zf(i,k)
                ph(i,k)  = 0.5_r8*(p(i,k) + p(i,k-1))
             end if

            dv(i,k) = 8.794E-5_r8*t(i,k)**1.81_r8/ph(i,k)
            mua(i,k) = 1.496E-6_r8*t(i,k)**1.5_r8/ &
                     (t(i,k)+120._r8)

            rho(i,k) = rhoh(i,k)

! air density adjustment for fallspeed parameters
! add air density correction factor to the power of 
! 0.54 following Heymsfield and Bansemer 2006

            arn(i,k)=ar*(rhosu/rho(i,k))**0.54
            asn(i,k)=as*(rhosu/rho(i,k))**0.54
            acn(i,k)=ac*(rhosu/rho(i,k))**0.54
            ain(i,k)=ai*(rhosu/rho(i,k))**0.54

          end do
        end do

! initialize aerosol number
        do k=1,pver
          do i=1,il2g
            naer2(i,k,:)=0._r8
            naer2h(i,k,:)=0._r8
            dum2l(i,k)=0._r8
            dum2i(i,k)=0._r8
          end do
        end do

        do k=1,pver
          do i=1,il2g
            ntaer(i,k) = 0.0_r8
            ntaerh(i,k) = 0.0_r8
            do m=1,naer_cu

              maerosol(1,m)=aer_mmr(i,k,m)*rhom(i,k)

!------------------------------------------------------------------           
        
! set number nucleated for sulfate based on Lohmann et al. 2000 (JGR) Eq.2
!    Na=340.*(massSO4)^0.58  where Na=cm-3 and massSO4=ug/m3
! convert units to Na [m-3] and SO4 [kgm-3]
!    Na(m-3)= 1.e6 cm3 m-3 Na(cm-3)=340. *(massSO4[kg/m3]*1.e9ug/kg)^0.58
!    or Na(m-3)= 1.e6* 340.*(1.e9ug/kg)^0.58 * (massSO4[kg/m3])^0.58
              if(m .eq. idxsul) then
                naer2(i,k,m)= 5.64259e13_r8 * maerosol(1,m)**0.58
              else
                naer2(i,k,m)=maerosol(1,m)*num_to_mass_aer(m)
              endif
                ntaer(i,k) = ntaer(i,k) + naer2(i,k,m)
            enddo
          end do ! i loop
        end do ! k loop

        do i=1,il2g
          ltrue(i)=0
          do k=1,pver
            if (qc(i,k).ge.qsmall.or.qi(i,k).ge.qsmall.or.cmel(i,k).ge.qsmall.or.cmei(i,k).ge.qsmall) ltrue(i)=1
!            print *,'qc flag =',ltrue(i)
          end do
        end do

! skip microphysical calculations if no cloud water
      do i=1,il2g
        if (ltrue(i).eq.0) then
          do k=1,pver
            qctend(i,k)=0._r8
            qitend(i,k)=0._r8
            qnitend(i,k)=0._r8
            qrtend(i,k)=0._r8
            nctend(i,k)=0._r8
            nitend(i,k)=0._r8
            nrtend(i,k)=0._r8
            nstend(i,k)=0._r8
            qniic(i,k)=0._r8
            qric(i,k)=0._r8
            nsic(i,k)=0._r8
            nric(i,k)=0._r8
            qni(i,k)=0._r8
            qr(i,k)=0._r8
            ns(i,k)=0._r8
            nr(i,k)=0._r8
            qc(i,k) = 0._r8
            qi(i,k) = 0._r8
            nc(i,k) = 0._r8
            ni(i,k) = 0._r8
            rprd(i,k) = 0._r8
            sprd(i,k) = 0._r8
            frz(i,k) = 0._r8
          end do
          goto 300
        end if

        kqc(i) = 1
        kqi(i) = 1
        lcbase(i) = .true.
        libase(i) = .true.

! assign number of steps for iteration
! use 2 steps following Song and Zhang, 2011, J. Clim.
        iter = 2  !5 !Shaocai Yu !2

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  iteration
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        do it=1,iter
! initialize sub-step microphysical tendencies
         do k=1,pver
           qctend(i,k)=0._r8
           qitend(i,k)=0._r8
           qnitend(i,k)=0._r8
           qrtend(i,k)=0._r8
           nctend(i,k)=0._r8
           nitend(i,k)=0._r8
           nrtend(i,k)=0._r8
           nstend(i,k)=0._r8
           rprd(i,k) = 0._r8
           sprd(i,k) = 0._r8
           frz(i,k)  = 0._r8
           qniic(i,k)=0._r8
           qric(i,k)=0._r8
           nsic(i,k)=0._r8
           nric(i,k)=0._r8
           qiic(i,k)=0._r8
           qcic(i,k)=0._r8
           niic(i,k)=0._r8
           ncic(i,k)=0._r8
!<songxl 2012-01-06---------------
            accrlm(i,k) = 0._r8
            bergnm(i,k) = 0._r8
            fhtimm(i,k) = 0._r8
            fhtctm(i,k) = 0._r8
            fhmlm (i,k) = 0._r8
            hmpim (i,k) = 0._r8
            accslm(i,k) = 0._r8
            dlfm  (i,k) = 0._r8

            autoln(i,k) = 0._r8
            accrln(i,k) = 0._r8
            bergnn(i,k) = 0._r8
            fhtimn(i,k) = 0._r8
            fhtctn(i,k) = 0._r8
            fhmln (i,k) = 0._r8
            accsln(i,k) = 0._r8
            activn(i,k) = 0._r8
            dlfn  (i,k) = 0._r8
            ncadj (i,k) = 0._r8
            qcadj (i,k) = 0._r8

            autoim(i,k) = 0._r8
            accsim(i,k) = 0._r8
            difm  (i,k) = 0._r8

            nuclin(i,k) = 0._r8
            nuclim(i,k) = 0._r8
            autoin(i,k) = 0._r8
            accsin(i,k) = 0._r8
            hmpin (i,k) = 0._r8
            difn  (i,k) = 0._r8
            niadj (i,k) = 0._r8
            qiadj (i,k) = 0._r8

            trspcm(i,k) = 0._r8
            trspcn(i,k) = 0._r8
            trspim(i,k) = 0._r8
            trspin(i,k) = 0._r8

            effc(i,k) = 0._r8
            effi(i,k) = 0._r8
            effs(i,k) = 0._r8

            fhmrm (i,k) = 0._r8
!songxl 2012-01-06>---------------
!songxl 2012-01-06>---------------
         end do

!---------------------Shaocai
!        goto 9910  !Stop2
        
         do k = pver,msg+2,-1

!within the cloud processing...
            if (k > jt(i) .and. k <= jb(i) .and. eps0(i) > 0._r8      &
              .and.mu(i,k).gt.0._r8 .and. mu(i,k-1).gt.0._r8) then

! initialize precip fallspeeds to zero
            ums(k)=0._r8
            uns(k)=0._r8
            umr(k)=0._r8
            unr(k)=0._r8
            prf(k)=0._r8
            pnrf(k)=0._r8
            psf(k) =0._r8
            pnsf(k) = 0._r8
            ttend(k)=0._r8
            nnuccd(k)=0._r8
            npccn(k)=0._r8

!************************************************************************************
! obtain values of cloud water/ice mixing ratios and number concentrations in
! updraft
! for microphysical process calculations
! units are kg/kg for mixing ratio, 1/kg for number conc
!************************************************************************************

! limit values to 0.005 kg/kg
!dkay      qc(i,k)=min(qc(i,k),5.e-3_r8)
!dkay      qi(i,k)=min(qi(i,k),5.e-3_r8)
           nc(i,k)=max(nc(i,k),0._r8)
           ni(i,k)=max(ni(i,k),0._r8)
           if (it.eq.1) then
             qcic(i,k) = qc(i,k)
             qiic(i,k) = qi(i,k)
!             print *,'at it=1',qcic(i,k),k,it
!dkay
!            qcic(i,k) = qc_kf_act(k)
!            qiic(i,k) = qi_kf_act(k)
!dkay
             ncic(i,k) = nc(i,k)
             niic(i,k) = ni(i,k)
             qniic(i,k)= qni(i,k)
             qric(i,k) = qr(i,k)
             nsic(i,k) = ns(i,k)
             nric(i,k) = nr(i,k)
           else  ! for it 
             if (k.le.kqc(i)) then
                qcic(i,k) = qc(i,k)
                ncic(i,k) = nc(i,k)
                if (k.eq.kqc(i)) then
                  qcic(i,k) = qc(i,k-1)
                  ncic(i,k) = nc(i,k-1)
                end if
! consider rain falling from above
                do kk= k,jt(i)+2,-1
                   qric(i,k) = qr(i,k) + max(0._r8, qr(i,kk-1)-qr(i,kk-2) )
                   if (qr(i,kk-1) .gt. 0._r8)  &
                   nric(i,k) = nr(i,k) + max(0._r8,qr(i,kk-1)-qr(i,kk-2))/qr(i,kk-1)*nr(i,kk-1)
                end do
             end if
             if(k.le.kqi(i)) then
                qiic(i,k) = qi(i,k)
                niic(i,k) = ni(i,k)
                if(k.eq.kqi(i)) then
                  qiic(i,k) = qi(i,k-1)
                  niic(i,k) = ni(i,k-1)
                end if
! consider snow falling from above
                do kk= k,jt(i)+2,-1
                  qniic(i,k) = qni(i,k) + max(0._r8, qni(i,kk-1)-qni(i,kk-2) )
                  if (qni(i,kk-1) .gt. 0._r8)  &
                  nsic(i,k) = ns(i,k) + max(0._r8,qni(i,kk-1)-qni(i,kk-2))/qni(i,kk-1)*ns(i,kk-1)
                end do
             end if
           end if

            if(it.eq.1) then
!              print
!              *,'qcic,qiic=',qcic(i,k),qiic(i,k),i,k,cmel(i,k),cmei(i,k),tu(i,k),it
            end if
!**********************************************************************
! boundary condition for cloud liquid water and cloud ice
!***********************************************************************

! boundary condition for provisional cloud water
        if (cmel(i,k-1).gt.qsmall .and. lcbase(i) .and. it.eq.1 ) then
             kqc(i) = k
             lcbase(i) = .false.
             qcic(i,k) = dz(i,k)*cmel(i,k-1)/(mu(i,k-1)+dz(i,k)*du(i,k-1))
             if(qcic(i,k).eq.0.0) then
              if(it.eq.1) then
!              print *,'dz,cmel...',
!              dz(i,k),cmel(i,k+1),mu(i,k+1),dz(i,k),du(i,k+1)
              end if
             end if
             ncic(i,k) = qcic(i,k)/(4._r8/3._r8*pi*11.e-6_r8**3*rhow)
         end if

! boundary condition for provisional cloud ice
         if (qiic(i,k).gt.qsmall .and. libase(i) .and. it.eq.1 ) then
             kqi(i) = k
             libase(i) = .false.
         else if ( cmei(i,k-1).gt.qsmall .and.   &
             cmei(i,k).lt.qsmall .and. k.lt.jb(i) .and. libase(i) .and. it.eq.1) then
             kqi(i)=k
             libase(i) = .false.
             qiic(i,k) = dz(i,k)*cmei(i,k-1)/(mu(i,k-1)+dz(i,k)*du(i,k-1))
             niic(i,k) = qiic(i,k)/(4._r8/3._r8*pi*25.e-6_r8**3*rhoi)
         end if

!***************************************************************************
! get size distribution parameters based on in-cloud cloud water/ice 
! these calculations also ensure consistency between number and mixing ratio
!***************************************************************************
! cloud ice
           if (qiic(i,k).ge.qsmall) then
! add upper limit to in-cloud number concentration to prevent numerical error
              niic(i,k)=min(niic(i,k),qiic(i,k)*1.e20_r8)
              lami(k) = (mskf_GAMMA(1._r8+di)*ci* &
                  niic(i,k)/qiic(i,k))**(1._r8/di)
              n0i(k) = niic(i,k)*lami(k)
! check for slope
              lammax = 1._r8/10.e-6_r8
              lammin = 1._r8/(2._r8*dcs)
! adjust vars
              if (lami(k).lt.lammin) then
                lami(k) = lammin
                n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*mskf_GAMMA(1._r8+di))
                niic(i,k) = n0i(k)/lami(k)
              else if (lami(k).gt.lammax) then
                lami(k) = lammax
                n0i(k) = lami(k)**(di+1._r8)*qiic(i,k)/(ci*mskf_GAMMA(1._r8+di))
                niic(i,k) = n0i(k)/lami(k)
              end if
           else
              lami(k) = 0._r8
              n0i(k) = 0._r8
           end if

!cloud water
           if (qcic(i,k).ge.qsmall) then

! add upper limit to in-cloud number concentration to prevent numerical error
              ncic(i,k)=min(ncic(i,k),qcic(i,k)*1.e20_r8)

! get pgam from fit to observations of martin et al. 1994

              pgam(k)=0.0005714_r8*(ncic(i,k)/1.e6_r8*rho(i,k))+0.2714_r8  ! TWG 2017 change / to * for consistenct with Morrison
              pgam(k)=1._r8/(pgam(k)**2)-1._r8
              pgam(k)=max(pgam(k),2._r8)
              pgam(k)=min(pgam(k),15._r8)

! calculate lamc
              lamc(k) = (pi/6._r8*rhow*ncic(i,k)*mskf_GAMMA(pgam(k)+4._r8)/ &
                 (qcic(i,k)*mskf_GAMMA(pgam(k)+1._r8)))**(1._r8/3._r8)

! lammin, 50 micron diameter max mean size
              lammin = (pgam(k)+1._r8)/50.e-6_r8
              lammax = (pgam(k)+1._r8)/2.e-6_r8

              if (lamc(k).lt.lammin) then
                 lamc(k) = lammin
                 ncic(i,k) = 6._r8*lamc(k)**3*qcic(i,k)* &
                      mskf_GAMMA(pgam(k)+1._r8)/ &
                      (pi*rhow*mskf_GAMMA(pgam(k)+4._r8))
              else if (lamc(k).gt.lammax) then
                 lamc(k) = lammax
                 ncic(i,k) = 6._r8*lamc(k)**3*qcic(i,k)* &
                       mskf_GAMMA(pgam(k)+1._r8)/ &
                      (pi*rhow*mskf_GAMMA(pgam(k)+4._r8))
              end if

! parameter to calculate droplet freezing

              cdist1(k) = ncic(i,k)/mskf_GAMMA(pgam(k)+1._r8)
           else
              lamc(k) = 0._r8
              cdist1(k) = 0._r8
           end if
! boundary condition for cloud liquid water
         if ( kqc(i) .eq. k  ) then
              qc(i,k) =  0._r8
              nc(i,k) = 0._r8
          end if
! boundary condition for cloud ice
          if (kqi(i).eq.k  ) then
             qi(i,k) = 0._r8
             ni(i,k) = 0._r8
          end if

!**************************************************************************
! begin micropysical process calculations 
!**************************************************************************

!.................................................................
! autoconversion of cloud liquid water to rain
! formula from Khrouditnov and Kogan (2000)
! minimum qc of 1 x 10^-8 prevents floating point error

           if (qcic(i,k).ge.1.e-8_r8) then

! nprc is increase in rain number conc due to autoconversion
! nprc1 is decrease in cloud droplet conc due to autoconversion

!TWG Feb 2017 Update from Khrouditnov and Kogan (2000) to Kogan (2013) for
!convection

!              prc(k) = 1350._r8*qcic(i,k)**2.47_r8*    &
!                    (ncic(i,k)/1.e6_r8*rho(i,k))**(-1.79_r8)
              prc(k) = 7.98E10_r8*qcic(i,k)**4.22_r8*    &
                    (ncic(i,k)/1.e6_r8*rho(i,k))**(-3.01_r8)

              nprc(k) = prc(k)/(4._r8/3._r8*pi*rhow*(25.e-6_r8)**3)
              nprc1(k) = prc(k)/(qcic(i,k)/ncic(i,k))
           else
              prc(k)=0._r8
              nprc(k)=0._r8
              nprc1(k)=0._r8
           end if

! provisional rain mixing ratio and number concentration (qric and nric)
! at boundary are estimated via autoconversion

         if (k.eq.kqc(i) .and. it.eq.1) then
             qric(i,k) = prc(k)*dz(i,k)/0.55_r8
             nric(i,k) = nprc(k)*dz(i,k)/0.55_r8
             qr(i,k) = 0.0_r8
             nr(i,k) = 0.0_r8
         end if
!          print *,'qric,nric,qr,nr afer autoconversion cld water to rain'
!          print *, 'qric=',qric
!          print *,
!          'nric=',nric(i,15),i,nprc(15),prc(15),ncic(i,15),rhow,qcic(i,15)
!          print *, 'qr=',qr
!          print *, 'qr=',qr

!.......................................................................

! similar to Ferrier (1994)

           if (t(i,k).le.273.15_r8.and.qiic(i,k).ge.qsmall) then

! note: assumes autoconversion timescale of 180 sec  !TWG Feb17 adjust
! autoconversion time scale
              ! nprci(k) = n0i(k)/(lami(k)*180._r8)*exp(-lami(k)*dcs)
              nprci(k) = n0i(k)/(lami(k)*F14)*exp(-lami(k)*dcs)
             ! prci(k) = pi*rhoi*n0i(k)/(6._r8*180._r8)* &
               prci(k) = pi*rhoi*n0i(k)/(6._r8*F14)* &  
                  (dcs**3/lami(k)+3._r8*dcs**2/lami(k)**2+ &
                  6._r8*dcs/lami(k)**3+6._r8/lami(k)**4)*exp(-lami(k)*dcs)
           else
              prci(k)=0._r8
              nprci(k)=0._r8
           end if

! provisional snow mixing ratio and number concentration (qniic and nsic) 
! at boundary are estimated via autoconversion

           if (k.eq.kqi(i) .and. it.eq.1) then
              qniic(i,k)= prci(k)*dz(i,k)*0.25_r8
              nsic(i,k)= nprci(k)*dz(i,k)*0.25_r8
              qni(i,k)= 0.0_r8
              ns(i,k)= 0.0_r8
           end if
! if precip mix ratio is zero so should number concentration
           if (qniic(i,k).lt.qsmall) then
              qniic(i,k)=0._r8
              nsic(i,k)=0._r8
           end if
           if (qric(i,k).lt.qsmall) then
              qric(i,k)=0._r8
              nric(i,k)=0._r8
           end if

! make sure number concentration is a positive number to avoid 
! taking root of negative later
           nric(i,k)=max(nric(i,k),0._r8)
           nsic(i,k)=max(nsic(i,k),0._r8)

!**********************************************************************
! get size distribution parameters for precip
!**********************************************************************
! rain

           if (qric(i,k).ge.qsmall) then
             lamr(k) = (pi*rhow*nric(i,k)/qric(i,k))**(1._r8/3._r8)
             n0r(k) = nric(i,k)*lamr(k)
! check for slope
             lammax = 1._r8/20.e-6_r8
             lammin = 1._r8/500.e-6_r8
! adjust vars
             if (lamr(k).lt.lammin) then
               lamr(k) = lammin
               n0r(k) = lamr(k)**4*qric(i,k)/(pi*rhow)
               nric(i,k) = n0r(k)/lamr(k)
             else if (lamr(k).gt.lammax) then
               lamr(k) = lammax
               n0r(k) = lamr(k)**4*qric(i,k)/(pi*rhow)
               nric(i,k) = n0r(k)/lamr(k)
             end if

! provisional rain number and mass weighted mean fallspeed (m/s)
! Eq.18 of Morrison and Gettelman, 2008, J. Climate
             unr(k) = min(arn(i,k)*mskf_GAMMA(1._r8+br)/lamr(k)**br,10._r8)
             umr(k) = min(arn(i,k)*mskf_GAMMA(4._r8+br)/(6._r8*lamr(k)**br),10._r8)
           else
             lamr(k) = 0._r8
             n0r(k) = 0._r8
             umr(k) = 0._r8
             unr(k) = 0._r8
           end if
!......................................................................
! snow
           if (qniic(i,k).ge.qsmall) then
             lams(k) = (mskf_GAMMA(1._r8+ds)*cs*nsic(i,k)/ &
                       qniic(i,k))**(1._r8/ds)
             n0s(k) = nsic(i,k)*lams(k)

! check for slope
             lammax = 1._r8/10.e-6_r8
             lammin = 1._r8/2000.e-6_r8
! adjust vars
             if (lams(k).lt.lammin) then
               lams(k) = lammin
               n0s(k) = lams(k)**(ds+1._r8)*qniic(i,k)/(cs*mskf_GAMMA(1._r8+ds))
               nsic(i,k) = n0s(k)/lams(k)
             else if (lams(k).gt.lammax) then
               lams(k) = lammax
               n0s(k) = lams(k)**(ds+1._r8)*qniic(i,k)/(cs*mskf_GAMMA(1._r8+ds))
               nsic(i,k) = n0s(k)/lams(k)
             end if

! provisional snow number and mass weighted mean fallspeed (m/s)
             ums(k) = min(asn(i,k)*mskf_GAMMA(4._r8+bs)/(6._r8*lams(k)**bs),3.6_r8)
             uns(k) = min(asn(i,k)*mskf_GAMMA(1._r8+bs)/lams(k)**bs,3.6_r8)
           else
             lams(k) = 0._r8
             n0s(k) = 0._r8
             ums(k) = 0._r8
             uns(k) = 0._r8
           end if

!.......................................................................
! snow self-aggregation from passarelli, 1978, used by Reisner(1998,Eq.A.35)
! this is hard-wired for bs = 0.4 for now
! ignore self-collection of cloud ice

          if (qniic(i,k).ge.qsmall .and. t(i,k).le.273.15_r8) then
              nsagg(k) = -1108._r8*asn(i,k)*Eii* &
                   pi**((1._r8-bs)/3._r8)*rhosn**((-2._r8-bs)/3._r8)*rho(i,k)** &
                   ((2._r8+bs)/3._r8)*qniic(i,k)**((2._r8+bs)/3._r8)* &
                   (nsic(i,k)*rho(i,k))**((4._r8-bs)/3._r8)/ &
                   (4._r8*720._r8*rho(i,k))
           else
              nsagg(k)=0._r8
           end if

!.......................................................................
! accretion of cloud droplets onto snow/graupel
! here use continuous collection equation with
! simple gravitational collection kernel
! ignore collisions between droplets/cloud ice

! ignore collision of snow with droplets above freezing

           if (qniic(i,k).ge.qsmall .and. t(i,k).le.273.15_r8 .and. &
              qcic(i,k).ge.qsmall) then

! put in size dependent collection efficiency
! mean diameter of snow is area-weighted, since
! accretion is function of crystal geometric area
! collection efficiency is from stoke's law (Thompson et al. 2004)

              dc0 = (pgam(k)+1._r8)/lamc(k)
              ds0 = 1._r8/lams(k)
              dum = dc0*dc0*uns(k)*rhow/(9._r8*mua(i,k)*ds0)
              eci = dum*dum/((dum+0.4_r8)*(dum+0.4_r8))
              eci = max(eci,0._r8)
              eci = min(eci,1._r8)

              psacws(k) = pi/4._r8*asn(i,k)*qcic(i,k)*rho(i,k)* &
                  n0s(k)*Eci*mskf_GAMMA(bs+3._r8)/ &
                  lams(k)**(bs+3._r8)   
              npsacws(k) = pi/4._r8*asn(i,k)*ncic(i,k)*rho(i,k)* &
                  n0s(k)*Eci*mskf_GAMMA(bs+3._r8)/ &
                  lams(k)**(bs+3._r8)
           else
              psacws(k)=0._r8
              npsacws(k)=0._r8
           end if

! secondary ice production due to accretion of droplets by snow 
! (Hallet-Mossop process) (from Cotton et al., 1986)
           if((t(i,k).lt.270.16_r8) .and. (t(i,k).ge.268.16_r8)) then
              ni_secp   = 3.5e8_r8*(270.16_r8-t(i,k))/2.0_r8*psacws(k)
              nsacwi(k) = ni_secp
              msacwi(k) = min(ni_secp*mi0,psacws(k))
           else if((t(i,k).lt.268.16_r8) .and. (t(i,k).ge.265.16_r8)) then
              ni_secp   = 3.5e8_r8*(t(i,k)-265.16_r8)/3.0_r8*psacws(k)
              nsacwi(k) = ni_secp
              msacwi(k) = min(ni_secp*mi0,psacws(k))
           else
              ni_secp   = 0.0_r8
              nsacwi(k) = 0.0_r8
              msacwi(k) = 0.0_r8
           endif
           psacws(k) = max(0.0_r8,psacws(k)-ni_secp*mi0)

!.......................................................................
! accretion of rain water by snow
! formula from ikawa and saito, 1991, used by reisner et al., 1998

           if (qric(i,k).ge.1.e-8_r8 .and. qniic(i,k).ge.1.e-8_r8 .and. &
              t(i,k).le.273.15_r8) then

              pracs(k) = pi*pi*ecr*(((1.2_r8*umr(k)-0.95_r8*ums(k))**2+ &
                  0.08_r8*ums(k)*umr(k))**0.5_r8*rhow*rho(i,k)* &
                  n0r(k)*n0s(k)* &
                  (5._r8/(lamr(k)**6*lams(k))+ &
                  2._r8/(lamr(k)**5*lams(k)**2)+ &
                  0.5_r8/(lamr(k)**4*lams(k)**3)))

              npracs(k) = pi/2._r8*rho(i,k)*ecr*(1.7_r8*(unr(k)-uns(k))**2+ &
                  0.3_r8*unr(k)*uns(k))**0.5_r8*n0r(k)*n0s(k)* &
                  (1._r8/(lamr(k)**3*lams(k))+ &
                  1._r8/(lamr(k)**2*lams(k)**2)+ &
                  1._r8/(lamr(k)*lams(k)**3))
           else
              pracs(k)=0._r8
              npracs(k)=0._r8
           end if

!.......................................................................
! heterogeneous freezing of rain drops
! follows from Bigg (1953)

           if (t(i,k).lt.269.15_r8 .and. qric(i,k).ge.qsmall) then

              mnuccr(k) = 20._r8*pi*pi*rhow*nric(i,k)*bimm* &
                  exp(aimm*(273.15_r8-t(i,k)))/lamr(k)**3 &
                  /lamr(k)**3

              nnuccr(k) = pi*nric(i,k)*bimm* &
                   exp(aimm*(273.15_r8-t(i,k)))/lamr(k)**3
           else
              mnuccr(k)=0._r8
              nnuccr(k)=0._r8
           end if

!.......................................................................
! accretion of cloud liquid water by rain
! formula from Khrouditnov and Kogan (2000)
! gravitational collection kernel, droplet fall speed neglected

           if (qric(i,k).ge.qsmall .and. qcic(i,k).ge.qsmall) then
              pra(k) = 67._r8*(qcic(i,k)*qric(i,k))**1.15_r8
              npra(k) = pra(k)/(qcic(i,k)/ncic(i,k))
           else
              pra(k)=0._r8
              npra(k)=0._r8
           end if

!.......................................................................
! Self-collection of rain drops
! from Beheng(1994)

           if (qric(i,k).ge.qsmall) then
              nragg(k) = -8._r8*nric(i,k)*qric(i,k)*rho(i,k)
           else
              nragg(k)=0._r8
           end if

!.......................................................................
! Accretion of cloud ice by snow
! For this calculation, it is assumed that the Vs >> Vi
! and Ds >> Di for continuous collection

           if (qniic(i,k).ge.qsmall.and.qiic(i,k).ge.qsmall &
              .and.t(i,k).le.273.15_r8) then
              prai(k) = pi/4._r8*asn(i,k)*qiic(i,k)*rho(i,k)* &
                   n0s(k)*Eii*mskf_GAMMA(bs+3._r8)/ &
                   lams(k)**(bs+3._r8)  
              nprai(k) = pi/4._r8*asn(i,k)*niic(i,k)* &
                   rho(i,k)*n0s(k)*Eii*mskf_GAMMA(bs+3._r8)/ &
                   lams(k)**(bs+3._r8)
           else
              prai(k)=0._r8
              nprai(k)=0._r8
           end if

!.......................................................................
! fallout term
        prf(k)  = -umr(k)*qric(i,k)/dz(i,k)
        pnrf(k) = -unr(k)*nric(i,k)/dz(i,k)
        psf(k)  = -ums(k)*qniic(i,k)/dz(i,k)
        pnsf(k) = -uns(k)*nsic(i,k)/dz(i,k)

!........................................................................
! calculate vertical velocity in cumulus updraft

     if (k.eq.jb(i)) then
       zkine(i,jb(i)) = 0.5_r8
       wu   (i,jb(i)) = 1._r8
       zbuo (i,jb(i)) = (tu(i,jb(i))*(1._r8+retv*qu(i,jb(i)))-    &
                     th(i,jb(i))*(1._r8+retv*qh(i,jb(i))))/   &
                     (th(i,jb(i))*(1._r8+retv*qh(i,jb(i))))
     else
       if (.true.) then
!          print *,'before ecmwf qcs=',qc(i,k),qi(i,k),qr(i,k),k
! ECMWF formula
  !          print *,'using ecmwrf CKE, retv=',retv
           zbc = tu(i,k)*(1._r8+retv*qu(i,k)-qr(i,k)-qni(i,k)-qi(i,k)-qc(i,k))
           zbe = th(i,k)*(1._r8+retv*qh(i,k))
           zbuo(i,k) = (zbc-zbe)/zbe
           zbuoc= (zbuo(i,k)+zbuo(i,k+1))*0.5_r8
           zdkbuo = dz(i,k+1)*grav*zfacbuo*zbuoc
           zdken = min(.99_r8,(1._r8+cwdrag)*max(du(i,k),eu(i,k))*dz(i,k+1)/ &
                      max(1.e-10_r8,mu(i,k+1)))
           zkine(i,k) = (zkine(i,k+1)*(1._r8-zdken)+zdkbuo)/      &
                      (1._r8+zdken)
!        print *,'zkine=',(zkine(i,k)),dz(i,k),k
        else
! Gregory formula
           write(*,*) "Gregory vertical velocity"
           zbc = tu(i,k)*(1._r8+retv*qu(i,k))
           zbe = th(i,k)*(1._r8+retv*qh(i,k))
           zbuo(i,k) = (zbc-zbe)/zbe-qr(i,k)-qni(i,k)-qi(i,k)-qc(i,k)
           zbuoc= (zbuo(i,k)+zbuo(i,k+1))*0.5_r8
           zdkbuo = dz(i,k+1)*grav*zbuoc*(1.0-0.25)/6.
           zdken = du(i,k)*dz(i,k+1)/max(1.e-10_r8,mu(i,k+1))
           zkine(i,k) = (zkine(i,k+1)*(1._r8-zdken)+zdkbuo)/      &
                      (1._r8+zdken)
         end if
              wu(i,k) = min(15._r8,sqrt(2._r8*max(0.1_r8,zkine(i,k) )))
!dkay         wu(i,k) = wu_kf_act(k)
       end if
!             print *,'wu from cke= & kf',wu(i,k),wu_kf_act(k),k
!ckay
       arcf(i,k)= mu(i,k)/wu(i,k)

!............................................................................
! droplet activation
! calculate potential for droplet activation if cloud water is present
! formulation from Abdul-Razzak and Ghan (2000) and Abdul-Razzak et al. (1998),
! AR98

       naer2h(i,k,:) = 0.5_r8*(naer2(i,k,:) + naer2(i,k+1,:))
       ntaerh(i,k)   = 0.5_r8*(ntaer(i,k) + ntaer(i,k+1))

!      write(*,*)'naer2h(i,k,:)',naer2h(i,k,:)

       if (qcic(i,k).ge.qsmall.or.cmel(i,k+1).ge.qsmall ) then

!dkay
! added qsatzm
!         print *, 'before activate'
         call mskf_activate(wu(i,k),t(i,k),rho(i,k), &
                 naer2h(i,k,:), naer_cu,naer_cu, maerosol,  &
                 dispersion_aer,hygro_aer, density_aer, dum2,qsatzm(i,k))
!             print *,'ccn, massmixing ratio of aerosols='
!            print *, dum2, maerosol
         dum2l(i,k) = dum2
       else
         dum2l(i,k) = 0._r8
       end if

! get droplet activation rate
       if (qcic(i,k).ge.qsmall .and. t(i,k).gt.238.15_r8 .and. k.gt.jt(i)+2) then

! assume aerosols already activated are equal number of existing droplets for
! simplicity
         if (k.eq.kqc(i))  then
              npccn(k) = dum2l(i,k)/deltat
         else
              npccn(k) = (dum2l(i,k)-ncic(i,k))/deltat
         end if
! make sure number activated > 0
         npccn(k) = max(0._r8,npccn(k))
         ncmax = dum2l(i,k)
       else
         npccn(k)=0._r8
         ncmax = 0._r8
       end if

!..............................................................................
!ice nucleation

       esi(i,k)= mskf_polysvp(t(i,k),1)      ! over ice          
       es(i,k) = mskf_polysvp(t(i,k),0)
       qs(i,k) = 0.622_r8*es(i,k)/(ph(i,k) - (1.0_r8-0.622_r8)*es(i,k))
       qs(i,k) = min(1.0_r8,qs(i,k))
       if (qs(i,k) < 0.0_r8)  qs(i,k) = 1.0_r8

       relhum(i,k)= 1.0_r8

       if (t(i,k).lt.tmelt ) then
         if (.true.) then

! Liu et al.,J. climate, 2007
!         print *, 'before ice nuke'

            call mskf_nucleati(wu(i,k),t(i,k),p(i,k),q(i,k),qcic(i,k),rho(i,k),  & ! TWG add p and replace relhum with q
                         naer2h(i,k,:),naer_cu,dum2i(i,k) &
                        , nihf(i,k),     &
                        niimm(i,k),nidep(i,k),nimey(i,k))

             nihf(i,k)=nihf(i,k)*rho(i,k)           !  convert from #/kg -> #/m3)
             niimm(i,k)=niimm(i,k)*rho(i,k)
             nidep(i,k)=nidep(i,k)*rho(i,k)
             nimey(i,k)=nimey(i,k)*rho(i,k)
          else

! cooper curve (factor of 1000 is to convert from L-1 to m-3)
            dum2i(i,k)=0.005_r8*exp(0.304_r8*(273.15_r8-t(i,k)))*1000._r8
! put limit on number of nucleated crystals, set to number at T=-30 C
! cooper (limit to value at -35 C)
            dum2i(i,k)=min(dum2i(i,k),208.9e3_r8)/rho(i,k) ! convert from m-3 to kg-1
          endif
        else
          dum2i(i,k)=0._r8
        end if
!ckay
!       print *,'nucleated ccn=',dum2i(i,k),k

! ice nucleation if activated nuclei exist at t<0C 

        if (dum2i(i,k).gt.0._r8.and.t(i,k).lt.tmelt.and. &
           relhum(i,k)*es(i,k)/esi(i,k).gt. 1.05_r8  .and. k.gt.jt(i)+1) then

           if (k.eq.kqi(i)) then
                nnuccd(k)=dum2i(i,k)/deltat
           else
                nnuccd(k)=(dum2i(i,k)-niic(i,k))/deltat
           end if
           nnuccd(k)=max(nnuccd(k),0._r8)
           nimax = dum2i(i,k)

!Calc mass of new particles using new crystal mass...
!also this will be multiplied by mtime as nnuccd is...
           mnuccd(k) = nnuccd(k) * mi0
         else
           nnuccd(k)=0._r8
           nimax = 0._r8
           mnuccd(k) = 0._r8
         end if
!................................................................................
! Bergeron process
! If 0C< T <-40C and both ice and liquid exist

         if (t(i,k).le.273.15_r8 .and. t(i,k).gt.233.15_r8 .and.  &
              qiic(i,k).gt.0.5e-6_r8 .and. qcic(i,k).gt. qsmall)  then
              plevap = qcic(i,k)/bergtsf
              prb(k) = max(0._r8,plevap)
              nprb(k) = prb(k)/(qcic(i,k)/ncic(i,k))
         else
              prb(k)=0._r8
              nprb(k)=0._r8
         end if

!................................................................................
! heterogeneous freezing of cloud water (-5C < T < -35C)

        if (qcic(i,k).ge.qsmall .and.ncic(i,k).gt.0._r8 .and. ntaerh(i,k).gt.0._r8 .and.  &
              t(i,k).le.268.15_r8 .and. t(i,k).gt.238.15_r8 ) then

          if (.false.)  then
! immersion freezing (Diehl and Wurzler, 2004)
              ttend(k) = -grav*wu(i,k)/cp/(1.0_r8+gamhat(i,k))
              naimm = (0.00291_r8*naer2h(i,k,idxbcphi)+32.3_r8*(naer2h(i,k,idxdst1)  &
                      +naer2h(i,k,idxdst2)+naer2h(i,k,idxdst3)+              &
                       naer2h(i,k,idxdst4)))/ntaerh(i,k)             !m-3
              if (ttend(k) .lt. 0._r8) then
                 nnuccc(k) = -naimm*exp(273.15_r8-t(i,k))*ttend(k)*qcic(i,k)/rhow   ! kg-1s-1                        
                 mnuccc(k) = nnuccc(k)*qcic(i,k)/ncic(i,k)
              end if
          else


! immersion freezing (Bigg, 1953)
              mnuccc(k) = pi*pi/36._r8*rhow* &
                    cdist1(k)*mskf_GAMMA(7._r8+pgam(k))* &
                    bimm*exp(aimm*(273.15_r8-t(i,k)))/ &
                    lamc(k)**3/lamc(k)**3

              nnuccc(k) = pi/6._r8*cdist1(k)*mskf_GAMMA(pgam(k)+4._r8) &
                    *bimm*exp(aimm*(273.15_r8-t(i,k)))/lamc(k)**3
           end if

! contact freezing (Young, 1974) with hooks into simulated dust

           tcnt=(270.16_r8-t(i,k))**1.3_r8
           viscosity=1.8e-5_r8*(t(i,k)/298.0_r8)**0.85_r8    ! Viscosity (kg/m/s)          
           mfp=2.0_r8*viscosity/(ph(i,k)  &                  ! Mean free path (m)
               *sqrt(8.0_r8*28.96e-3_r8/(pi*8.314409_r8*t(i,k))))

           slip1=1.0_r8+(mfp/rn_dst1)*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rn_dst1/mfp))))! Slip correction factor
           slip2=1.0_r8+(mfp/rn_dst2)*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rn_dst2/mfp))))
           slip3=1.0_r8+(mfp/rn_dst3)*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rn_dst3/mfp))))
           slip4=1.0_r8+(mfp/rn_dst4)*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rn_dst4/mfp))))

           dfaer1=1.381e-23_r8*t(i,k)*slip1/(6._r8*pi*viscosity*rn_dst1)  !aerosol diffusivity (m2/s)
           dfaer2=1.381e-23_r8*t(i,k)*slip2/(6._r8*pi*viscosity*rn_dst2)
           dfaer3=1.381e-23_r8*t(i,k)*slip3/(6._r8*pi*viscosity*rn_dst3)
           dfaer4=1.381e-23_r8*t(i,k)*slip4/(6._r8*pi*viscosity*rn_dst4)

           nacon1=0.0_r8
           nacon2=0.0_r8
           nacon3=0.0_r8
           nacon4=0.0_r8


           if (idxdst1.gt.0) then
              nacon1=naer2(i,k,idxdst1)*tcnt *0.0_r8
           endif
           if (idxdst2.gt.0) then
              nacon2=naer2(i,k,idxdst2)*tcnt ! 1/m3
           endif
           if (idxdst3.gt.0) then
              nacon3=naer2(i,k,idxdst3)*tcnt
           endif
           if (idxdst4.gt.0) then
              nacon4=naer2(i,k,idxdst4)*tcnt
           endif

           mnucct(k) = (dfaer1*nacon1+dfaer2*nacon2+dfaer3*nacon3+dfaer4*nacon4)*pi*pi/3._r8*rhow* &
                       cdist1(k)*mskf_GAMMA(pgam(k)+5._r8)/lamc(k)**4

           nnucct(k) = (dfaer1*nacon1+dfaer2*nacon2+dfaer3*nacon3+dfaer4*nacon4)*2._r8*pi*  &
                       cdist1(k)*mskf_GAMMA(pgam(k)+2._r8)/lamc(k)

!              if (nnuccc(k).gt.nnuccd(k)) then
!                 dum=nnuccd(k)/nnuccc(k)
! scale mixing ratio of droplet freezing with limit
!                 mnuccc(k)=mnuccc(k)*dum
!                 nnuccc(k)=nnuccd(k)
!              end if

           else
             mnuccc(k) = 0._r8
             nnuccc(k) = 0._r8
             mnucct(k) = 0._r8
             nnucct(k) = 0._r8
           end if

!****************************************************************************************
! conservation to ensure no negative values of cloud water/precipitation
! in case microphysical process rates are large
! note: for check on conservation, processes are multiplied by omsm
! to prevent problems due to round off error

! since activation/nucleation processes are fast, need to take into account
! factor mtime = mixing timescale in cloud / model time step
! for now mixing timescale is assumed to be 15 min
!*****************************************************************************************

       mtime=deltat/900._r8
       mtimec=deltat/900._r8

       mtime = AMAX1(1.0,mtime)   !TWG remove time scale limitation from CAM5
       mtimec = AMAX1(1.0,mtimec)    

! conservation of qc

        qce = mu(i,k)*qc(i,k)+dz(i,k)*(cmel(i,k-1)-du(i,k-1)*qc(i,k))
        dum = arcf(i,k)*(pra(k)+prc(k)+prb(k)+mnuccc(k)+mnucct(k)+msacwi(k)+   &
                         psacws(k)  )*dz(i,k)
        if( qce.lt.0._r8)  then
          prc(k) = 0._r8
          pra(k) = 0._r8
          prb(k) = 0._r8
          mnuccc(k) = 0._r8
          mnucct(k) = 0._r8
          msacwi(k) = 0._r8
          psacws(k) = 0._r8
        else  if (dum.gt.qce) then
          ratio = qce/dum*omsm
          prc(k) = prc(k)*ratio
          pra(k) = pra(k)*ratio
          prb(k) = prb(k)*ratio
          mnuccc(k) = mnuccc(k)*ratio
          mnucct(k) = mnucct(k)*ratio
          msacwi(k) = msacwi(k)*ratio
          psacws(k) = psacws(k)*ratio
        end if

! conservation of nc
        nce = mu(i,k)*nc(i,k)+(arcf(i,k)*npccn(k)*mtimec-du(i,k-1)*nc(i,k))*dz(i,k)
        dum = arcf(i,k)*dz(i,k)*(nprc1(k)+npra(k)+nnuccc(k)+nnucct(k)+ &
              npsacws(k)+ nprb(k) )
        if (nce.lt.0._r8) then
          nprc1(k) = 0._r8
!          nprc(k) = 0._r8
          npra(k) = 0._r8
          nnuccc(k) = 0._r8
          nnucct(k) = 0._r8
          npsacws(k) = 0._r8
          nprb(k) = 0._r8
        else if (dum.gt.nce) then
          ratio = nce/dum*omsm
          nprc1(k) = nprc1(k)*ratio
          npra(k) = npra(k)*ratio
          nnuccc(k) = nnuccc(k)*ratio
          nnucct(k) = nnucct(k)*ratio
          npsacws(k) = npsacws(k)*ratio
          nprb(k) = nprb(k)*ratio
        end if

! conservation of qi
        qie = mu(i,k)*qi(i,k)+dz(i,k)*(cmei(i,k-1)-du(i,k-1)*qi(i,k)+  &
                   ( mnuccc(k)+mnucct(k)+msacwi(k)+prb(k))*arcf(i,k) )
        dum = arcf(i,k)*(prci(k)+ prai(k))*dz(i,k)
        if (qie.lt.0._r8) then
          prci(k) = 0._r8
          prai(k) = 0._r8
        else if (dum.gt.qie) then
          ratio = qie/dum*omsm
          prci(k) = prci(k)*ratio
          prai(k) = prai(k)*ratio
        end if

! conservation of ni
         nie = mu(i,k)*ni(i,k)+dz(i,k)*(nnuccd(k)*mtime*arcf(i,k)-du(i,k-1)*ni(i,k)  &
                       + nnucct(k)*arcf(i,k) )
         dum = arcf(i,k)*dz(i,k)*(-nsacwi(k)+nprci(k)+ &
               nprai(k))
         if( nie.lt.0._r8) then
           nsacwi(k)= 0._r8
           nprci(k) = 0._r8
           nprai(k) = 0._r8
         else  if (dum.gt.nie) then
           ratio = nie/dum*omsm
           nsacwi(k)= nsacwi(k)*ratio
           nprci(k) = nprci(k)*ratio
           nprai(k) = nprai(k)*ratio
         end if

! conservation of qr

        qre = mu(i,k)*qr(i,k)+dz(i,k)*(pra(k)+prc(k))*arcf(i,k)
        dum = arcf(i,k)*dz(i,k)*(pracs(k)+ mnuccr(k)-prf(k))
        if (qre.lt.0._r8) then
           prf(k) = 0._r8
           pracs(k) = 0._r8
           mnuccr(k) = 0._r8
        else if (dum.gt.qre) then
           ratio = qre/dum*omsm
           prf(k) = prf(k)*ratio
           pracs(k) = pracs(k)*ratio
           mnuccr(k) = mnuccr(k)*ratio
        end if

! conservation of nr
         nre = mu(i,k)*nr(i,k)
         dum = arcf(i,k)*dz(i,k)*(-nprc(k)+npracs(k)+nnuccr(k) &
                   -nragg(k)-pnrf(k))
         if(nre.lt.0._r8) then
           nprc(k) = 0._r8
           npracs(k)= 0._r8
           nnuccr(k)= 0._r8
           nragg(k) = 0._r8
           pnrf(k) = 0._r8
         else if (dum.gt.nre) then
           ratio = nre/dum*omsm
           nprc(k) = nprc(k)*ratio
           npracs(k)= npracs(k)*ratio
           nnuccr(k)= nnuccr(k)*ratio
           nragg(k) = nragg(k)*ratio
           pnrf(k) = pnrf(k)*ratio
         end if

! conservation of qni

        qnie = mu(i,k)*qni(i,k)+dz(i,k)*( (prai(k)+psacws(k)+prci(k)+     &
                   pracs(k)+mnuccr(k))*arcf(i,k) )
        dum = arcf(i,k)*dz(i,k)*(-psf(k))

        if(qnie.lt.0._r8) then
           psf(k) = 0._r8
        else if (dum.gt.qnie) then
           ratio = qnie/dum*omsm
           psf(k) = psf(k)*ratio
        end if

! conservation of ns
        nse = mu(i,k)*ns(i,k)+dz(i,k)*(nprci(k)+nnuccr(k))*arcf(i,k)
        dum = arcf(i,k)*dz(i,k)*(-nsagg(k)-pnsf(k))
        if (nse.lt.0._r8) then
           nsagg(k) = 0._r8
           pnsf(k) = 0._r8
        else if (dum.gt.nse) then
           ratio = nse/dum*omsm
           nsagg(k) = nsagg(k)*ratio
           pnsf(k) = pnsf(k)*ratio
        end if

!*****************************************************************************
! get tendencies due to microphysical conversion processes
!*****************************************************************************

      if (k.le.kqc(i))   then
        qctend(i,k) = qctend(i,k)+  &
                 (-pra(k)-prc(k)-prb(k)-mnuccc(k)-mnucct(k)-msacwi(k)- &
                  psacws(k))

!             print *,'qctend components=',qctend(i,k),pra(k),prc(k), &
!             mnuccc(k)-mnucct(k)-msacwi(k),psacws(k)

        qitend(i,k) = qitend(i,k)+  &
                  (prb(k)+mnuccc(k)+mnucct(k)+msacwi(k)-prci(k)- &
                  prai(k)+mnuccd(k)*mtimec) !TWG ice nucleation change

        qrtend(i,k) = qrtend(i,k)+ &
                 (pra(k)+prc(k))+(-pracs(k)- &
                  mnuccr(k))


        qnitend(i,k) = qnitend(i,k)+ &
                (prai(k)+psacws(k)+prci(k))+( &
                   pracs(k)+mnuccr(k))

! multiply activation/nucleation by mtime to account for fast timescale

        nctend(i,k) = nctend(i,k)+ npccn(k)*mtimec+&
                  (-nnuccc(k)-nnucct(k)-npsacws(k) &
                  -npra(k)-nprc1(k)-nprb(k))

        nitend(i,k) = nitend(i,k)+ nnuccd(k)*mtime+&
                  (nnuccc(k)+ nnucct(k)+nsacwi(k)-nprci(k)- &
                  nprai(k))

        nstend(i,k) = nstend(i,k)+( &
                  nsagg(k)+nnuccr(k))+nprci(k)

        nrtend(i,k) = nrtend(i,k)+ &
                  nprc(k)+(-npracs(k)-nnuccr(k) +nragg(k))

! for output
! cloud liquid water-------------
        autolm(i,k) = -prc(k)*arcf(i,k)
        accrlm(i,k) = -pra(k)*arcf(i,k)
        bergnm(i,k) = -prb(k)*arcf(i,k)
        fhtimm(i,k) = -mnuccc(k)*arcf(i,k)
        fhtctm(i,k) = -mnucct(k)*arcf(i,k)
        hmpim (i,k) = -msacwi(k)*arcf(i,k)
        accslm(i,k) = -psacws(k)*arcf(i,k)
        collrm(i,k) = -pracs(k)*arcf(i,k)
        collrn(i,k) = -npracs(k)*arcf(i,k)
        fhtcrm(i,k) = -mnuccr(k)*arcf(i,k)
        fhtcrn(i,k) = -nnuccr(k)*arcf(i,k)
        dlfm  (i,k) = -du(i,k)*qc(i,k)

        autoln(i,k) = -nprc1(k)*arcf(i,k)*rho(i,k)
        autorn(i,k) = -nprc(k)*arcf(i,k)*rho(i,k)
        aggrn(i,k) =  nragg(k)*arcf(i,k)*rho(i,k)
        aggsn(i,k) =  nsagg(k)*arcf(i,k)*rho(i,k)
        accrln(i,k) = -npra(k)*arcf(i,k)*rho(i,k)
        bergnn(i,k) = -nprb(k)*arcf(i,k)*rho(i,k)
        fhtimn(i,k) = -nnuccc(k)*arcf(i,k)*rho(i,k)
        fhtctn(i,k) = -nnucct(k)*arcf(i,k)*rho(i,k)
        accsln(i,k) = -npsacws(k)*arcf(i,k)*rho(i,k)
        activn(i,k) = npccn(k)*mtimec*arcf(i,k)*rho(i,k)
        dlfn  (i,k) = -du(i,k)*nc(i,k)*rho(i,k)
!cloud ice------------------------        
        autoim(i,k) = -prci(k)*arcf(i,k)
        accsim(i,k) = -prai(k)*arcf(i,k)
        difm  (i,k) = -du(i,k)*qi(i,k)             !TWG 2017 change -du(i,k+1)*qi(i,k)

        nuclin(i,k) = nnuccd(k)*mtime*arcf(i,k)*rho(i,k)
        nuclim(i,k) = mnuccd(k)*mtime*arcf(i,k)*rho(i,k)
        autoin(i,k) = -nprci(k)*arcf(i,k)*rho(i,k)
        accsin(i,k) = -nprai(k)*arcf(i,k)*rho(i,k)
        hmpin (i,k)  = nsacwi(k)*arcf(i,k)*rho(i,k)
        difn  (i,k) = -du(i,k)*ni(i,k)*rho(i,k)
      else
        qctend(i,k) = 0._r8
        qitend(i,k) = 0._r8
        qrtend(i,k) = 0._r8
        qnitend(i,k) = 0._r8
        nctend(i,k) = 0._r8
        nitend(i,k) = 0._r8
        nstend(i,k) = 0._r8
        nrtend(i,k) = 0._r8
      end if

!********************************************************************************
!  vertical integration
!********************************************************************************
! snow
        if ( k.le.kqi(i) ) then
          qni(i,k-1) = 1._r8/mu(i,k-1)*                                    &
                   (mu(i,k)*qni(i,k)+dz(i,k)*(qnitend(i,k)+psf(k))*arcf(i,k) )

          ns(i,k-1) = 1._r8/mu(i,k-1)*                                    &
                   (mu(i,k)*ns(i,k)+dz(i,k)*(nstend(i,k)+pnsf(k))*arcf(i,k) )

         else
           qni(i,k-1)=0._r8
           ns(i,k-1)=0._r8
         end if

         if (qni(i,k-1).le.0._r8) then
          qni(i,k-1)=0._r8
          ns(i,k-1)=0._r8
         end if

! rain
         if (k.le.kqc(i) ) then
          qr(i,k-1) = 1._r8/mu(i,k-1)*                                    &
                   (mu(i,k)*qr(i,k)+dz(i,k)*(qrtend(i,k)+prf(k))*arcf(i,k) )

          nr(i,k-1) = 1._r8/mu(i,k-1)*                                    &
                   (mu(i,k)*nr(i,k)+dz(i,k)*(nrtend(i,k)+pnrf(k))*arcf(i,k) )

        else
          qr(i,k-1)=0._r8
          nr(i,k-1)=0._r8
        end if

        if( qr(i,k-1) .le. 0._r8) then
          qr(i,k-1)=0._r8
          nr(i,k-1)=0._r8
        end if

! freeze rain homogeneously at -40 C

         if (t(i,k-1) < 233.15_r8 .and. qr(i,k-1) > 0._r8) then

! make sure freezing rain doesn't increase temperature above threshold
          dum = xlf/cp*qr(i,k-1)
          if (t(i,k-1)+dum.gt.233.15_r8) then
              dum = -(t(i,k-1)-233.15_r8)*cp/xlf
!bugfix 2012-01-06              dum = dum/(xlf/cp*qr(i,k-1))
              dum = dum/qr(i,k-1)
              dum = max(0._r8,dum)
              dum = min(1._r8,dum)
          else
              dum = 1._r8
          end if
          qni(i,k-1)=qni(i,k-1)+dum*qr(i,k-1)
          ns(i,k-1)=ns(i,k-1)+dum*nr(i,k-1)
          qr(i,k-1)=(1._r8-dum)*qr(i,k-1)
          nr(i,k-1)=(1._r8-dum)*nr(i,k-1)
          fhmrm(i,k-1) = -mu(i,k-1)*dum*qr(i,k-1)/dz(i,k)
        end if

!        if( qr(i,k-1) .le. 0._r8) then
!          qr(i,k-1)=0._r8
!          nr(i,k-1)=0._r8
!        end if  

! cloud water
         if ( k.le.kqc(i) ) then
          qc(i,k-1) = 1._r8/mu(i,k-1)*                                    &
                   (mu(i,k)*qc(i,k)-dz(i,k)*du(i,k-1)*qc(i,k)             &
                    +dz(i,k)*qctend(i,k)*arcf(i,k)+dz(i,k)*cmel(i,k-1) )

          nc(i,k-1) = 1._r8/mu(i,k-1)*                                    &
                   (mu(i,k)*nc(i,k)-dz(i,k)*du(i,k-1)*nc(i,k)             &
                    +dz(i,k)*nctend(i,k)*arcf(i,k) )

        else
          qc(i,k-1)=0._r8
          nc(i,k-1)=0._r8
        end if

        qcorg = qc(i,k-1)
        ncorg = nc(i,k-1)
        if (qc(i,k-1).le. 0._r8) then
          qc(i,k-1)=0._r8
          nc(i,k-1)=0._r8
        end if
        qcadj(i,k-1)= (qc(i,k-1)- qcorg)*mu(i,k-1)/dz(i,k)*rho(i,k)
        ncadj(i,k-1)= (nc(i,k-1)- ncorg)*mu(i,k-1)/dz(i,k)*rho(i,k)

! cloud ice
         if( k.le.kqi(i)) then
           qi(i,k-1) = 1._r8/mu(i,k-1)*                                    &
                   (mu(i,k)*qi(i,k)-dz(i,k)*du(i,k-1)*qi(i,k)             &
                    +dz(i,k)*qitend(i,k)*arcf(i,k)+dz(i,k)*cmei(i,k-1) )

           ni(i,k-1) = 1._r8/mu(i,k-1)*                                    &
                   (mu(i,k)*ni(i,k)-dz(i,k)*du(i,k-1)*ni(i,k)             &
                    +dz(i,k)*nitend(i,k)*arcf(i,k) )

         else
          qi(i,k-1)=0._r8
          ni(i,k-1)=0._r8
         end if

        qiorg = qi(i,k-1)
        niorg = ni(i,k-1)
        if (qi(i,k-1).le. 0._r8) then
          qi(i,k-1)=0._r8
          ni(i,k-1)=0._r8
        end if
        qiadj(i,k-1)= (qi(i,k-1)- qiorg)*mu(i,k-1)/dz(i,k)*rho(i,k)
        niadj(i,k-1)= (ni(i,k-1)- niorg)*mu(i,k-1)/dz(i,k)*rho(i,k)

!        trspcm(i,k-1) = (mu(i,k)*qc(i,k) - mu(i,k-1)*qc(i,k-1))/dz(i,k)
!        trspcn(i,k-1) = (mu(i,k)*nc(i,k) -
!        mu(i,k-1)*nc(i,k-1))/dz(i,k)*rho(i,k)
!        trspim(i,k-1) = (mu(i,k)*qi(i,k) - mu(i,k-1)*qi(i,k-1))/dz(i,k)
!        trspin(i,k-1) = (mu(i,k)*ni(i,k) -
!        mu(i,k-1)*ni(i,k-1))/dz(i,k)*rho(i,k)


! freeze rain homogeneously at -38 C
         if (t(i,k-1) < 233.15_r8 .and. qc(i,k-1) > 0._r8) then
! make sure freezing rain doesn't increase temperature above threshold
          dum = xlf/cp*qc(i,k-1)
          if (t(i,k-1)+dum.gt.233.15_r8) then
              dum = -(t(i,k-1)-233.15_r8)*cp/xlf
!bugfix 2012-01-06      dum = dum/(xlf/cp*qc(i,k-1))
              dum = dum/qc(i,k-1)
              dum = max(0._r8,dum)
              dum = min(1._r8,dum)
          else
              dum = 1._r8
          end if
          qi(i,k-1)=qi(i,k-1)+dum*qc(i,k-1)
          ni(i,k-1)=ni(i,k-1)+dum*nc(i,k-1)
          fhmlm(i,k-1) = -mu(i,k-1)*dum*qc(i,k-1)/dz(i,k)
          fhmln(i,k-1) = -mu(i,k-1)*dum*nc(i,k-1)/dz(i,k)*rho(i,k)
          qc(i,k-1)=(1._r8-dum)*qc(i,k-1)
          nc(i,k-1)=(1._r8-dum)*nc(i,k-1)
        end if

        frz(i,k-1) = cmei(i,k-1) + arcf(i,k)*(prb(k)+mnuccc(k)+mnucct(k)+msacwi(k)+   &
                     pracs(k)+mnuccr(k)+psacws(k) )-fhmlm(i,k-1)-fhmrm(i,k-1)

!******************************************************************************
! get size distribution parameters based on in-cloud cloud water/ice
! these calculations also ensure consistency between number and mixing ratio

! following equation(2,3,4) of Morrison and Gettelman, 2008, J. Climate.
! Gamma(n)= (n-1)! 
! lamc <-> lambda for cloud liquid water
! pgam <-> meu    for cloud liquid water
! meu=0 for ice,rain and snow         
!*******************************************************************************
!songxl 2011-12-31

           niorg = ni(i,k-1)

! cloud ice
           if (qi(i,k-1).ge.qsmall) then
! add upper limit to in-cloud number concentration to prevent numerical error
              ni(i,k-1)=min(ni(i,k-1),qi(i,k-1)*1.e20_r8)
              lami(k-1) = (mskf_gamma(1._r8+di)*ci* &
                  ni(i,k-1)/qi(i,k-1))**(1._r8/di)
              n0i(k-1) = ni(i,k-1)*lami(k-1)
! check for slope
              lammax = 1._r8/10.e-6_r8
              lammin = 1._r8/(2._r8*dcs)
! adjust vars
              if (lami(k-1).lt.lammin) then
                lami(k-1) = lammin
                n0i(k-1) = lami(k-1)**(di+1._r8)*qi(i,k-1)/(ci*mskf_gamma(1._r8+di))
                ni(i,k-1) = n0i(k-1)/lami(k-1)
              else if (lami(k-1).gt.lammax) then
                lami(k-1) = lammax
                n0i(k-1) = lami(k-1)**(di+1._r8)*qi(i,k-1)/(ci*mskf_gamma(1._r8+di))
                ni(i,k-1) = n0i(k-1)/lami(k-1)
              end if
              effi(i,k-1) = 1.5_r8/lami(k-1)*1.e6_r8
           else
              lami(k-1) = 0._r8
              n0i(k-1) = 0._r8
              effi(i,k-1) = 0._r8
           end if

!songxl 2011-12-31-----
           niadj(i,k-1)= niadj(i,k-1)+(ni(i,k-1)-niorg)*mu(i,k-1)/dz(i,k)*rho(i,k)
!................................................................................
!songxl 2011-12-31
              ncorg = nc(i,k-1)

!cloud water
           if (qc(i,k-1).ge.qsmall) then

! add upper limit to in-cloud number concentration to prevent numerical error
              nc(i,k-1)=min(nc(i,k-1),qc(i,k-1)*1.e20_r8)

! get pgam from fit to observations of martin et al. 1994

              pgam(k-1)=0.0005714_r8*(nc(i,k-1)/1.e6_r8*rho(i,k-1))+0.2714_r8 !TWG 2017 change / to * in front of rho 
              pgam(k-1)=1._r8/(pgam(k-1)**2)-1._r8
              pgam(k-1)=max(pgam(k-1),2._r8)
              pgam(k-1)=min(pgam(k-1),15._r8)
! calculate lamc

              lamc(k-1) = (pi/6._r8*rhow*nc(i,k-1)*mskf_gamma(pgam(k-1)+4._r8)/ &
                 (qc(i,k-1)*mskf_gamma(pgam(k-1)+1._r8)))**(1._r8/3._r8)

! lammin, 50 micron diameter max mean size
              lammin = (pgam(k)+1._r8)/50.e-6_r8
              lammax = (pgam(k-1)+1._r8)/2.e-6_r8

              if (lamc(k-1).lt.lammin) then
                 lamc(k-1) = lammin
                 nc(i,k-1) = 6._r8*lamc(k-1)**3*qc(i,k-1)* &
                      mskf_gamma(pgam(k-1)+1._r8)/ &
                      (pi*rhow*mskf_gamma(pgam(k-1)+4._r8))
              else if (lamc(k-1).gt.lammax) then
                 lamc(k-1) = lammax
                 nc(i,k-1) = 6._r8*lamc(k-1)**3*qc(i,k-1)* &
                       mskf_gamma(pgam(k-1)+1._r8)/ &
                      (pi*rhow*mskf_gamma(pgam(k-1)+4._r8))
              end if
              effc(i,k-1) = mskf_gamma(pgam(k-1)+4._r8)/ &
                            mskf_gamma(pgam(k-1)+3._r8)/lamc(k-1)/2._r8*1.e6_r8
! parameter to calculate droplet freezing

              cdist1(k-1) = nc(i,k-1)/mskf_gamma(pgam(k-1)+1._r8)
           else
              lamc(k-1) = 0._r8
              cdist1(k-1) = 0._r8
              effc(i,k-1) = 0._r8
           end if

!songxl 2011-12-31-----
           ncadj(i,k-1) = ncadj(i,k-1)+ (nc(i,k-1)-ncorg)*mu(i,k-1)/dz(i,k)*rho(i,k)

           trspcm(i,k-1) = (mu(i,k)*qc(i,k) - mu(i,k-1)*qc(i,k-1))/dz(i,k)
           trspcn(i,k-1) = (mu(i,k)*nc(i,k) - mu(i,k-1)*nc(i,k-1))/dz(i,k)*rho(i,k)
           trspim(i,k-1) = (mu(i,k)*qi(i,k) - mu(i,k-1)*qi(i,k-1))/dz(i,k)
           trspin(i,k-1) = (mu(i,k)*ni(i,k) - mu(i,k-1)*ni(i,k-1))/dz(i,k)*rho(i,k)

           if (k-1 .eq. jt(i)+1)  then
             trspcm(i,k-2) =  mu(i,k-1)*qc(i,k-1)/dz(i,k)
             trspcn(i,k-2) =  mu(i,k-1)*nc(i,k-1)/dz(i,k)*rho(i,k)
             trspim(i,k-2) =  mu(i,k-1)*qi(i,k-1)/dz(i,k)
             trspin(i,k-2) =  mu(i,k-1)*ni(i,k-1)/dz(i,k)*rho(i,k)
             dlfm  (i,k-2) = -du(i,k-2)*qc(i,k-1)
             dlfn  (i,k-2) = -du(i,k-2)*nc(i,k-1)*rho(i,k)
             difm  (i,k-2) = -du(i,k-2)*qi(i,k-1)
             difn  (i,k-2) = -du(i,k-2)*ni(i,k-1)*rho(i,k)
           end if
!.......................................................................
! get size distribution parameters for precip
!......................................................................
! rain
           if (qr(i,k-1).ge.qsmall) then

             lamr(k-1) = (pi*rhow*nr(i,k-1)/qr(i,k-1))**(1._r8/3._r8)
             n0r(k-1) = nr(i,k-1)*lamr(k-1)
! check for slope
             lammax = 1._r8/20.e-6_r8
             lammin = 1._r8/500.e-6_r8
! adjust vars
             if (lamr(k-1).lt.lammin) then
               lamr(k-1) = lammin
               n0r(k-1) = lamr(k-1)**4*qr(i,k-1)/(pi*rhow)
               nr(i,k-1) = n0r(k-1)/lamr(k-1)
             else if (lamr(k-1).gt.lammax) then
               lamr(k-1) = lammax
               n0r(k-1) = lamr(k-1)**4*qr(i,k-1)/(pi*rhow)
               nr(i,k-1) = n0r(k-1)/lamr(k-1)
             end if
           else
             lamr(k-1) = 0._r8
             n0r(k-1) = 0._r8
           end if
!......................................................................
! snow
           if (qni(i,k-1).ge.qsmall) then
             lams(k-1) = (mskf_gamma(1._r8+ds)*cs*ns(i,k-1)/ &
                       qni(i,k-1))**(1._r8/ds)
             n0s(k-1) = ns(i,k-1)*lams(k-1)

! check for slope
             lammax = 1._r8/10.e-6_r8
             lammin = 1._r8/2000.e-6_r8
! adjust vars
             if (lams(k-1).lt.lammin) then
               lams(k-1) = lammin
               n0s(k-1) = lams(k-1)**(ds+1._r8)*qni(i,k-1)/(cs*mskf_gamma(1._r8+ds))
               ns(i,k-1) = n0s(k-1)/lams(k-1)
             else if (lams(k-1).gt.lammax) then
               lams(k-1) = lammax
               n0s(k-1) = lams(k-1)**(ds+1._r8)*qni(i,k-1)/(cs*mskf_gamma(1._r8+ds))
               ns(i,k-1) = n0s(k-1)/lams(k-1)
             end if
             effs(i,k-1) = 1.5_r8/lams(k-1)*1.e6_r8
           else
             lams(k-1) = 0._r8
             n0s(k-1) = 0._r8
             effs(i,k-1) = 0._r8
           end if

!dkay : since KF treats rain and snow separately, no need to add snow to the
!rprd (kg/kg/m)
!dkay        rprd(i,k-1)= (qnitend(i,k) + qrtend(i,k))*arcf(i,k)    ! original
        rprd(i,k-1)=  qrtend(i,k)  *arcf(i,k)
        sprd(i,k-1)=  qnitend(i,k) *arcf(i,k)

!dkay
!dkay       print *,'k,rprd,qrtend,qcic
!=',k,rprd(i,k-1),qrtend(i,k-1),qcic(i,k-1)
!dkay 
     end if  ! k<jlcl

! if rain/snow mix ratio is zero so should number concentration

         if (qni(i,k-1).lt.qsmall) then
           qni(i,k-1)=0._r8
           ns(i,k-1)=0._r8
         end if

         if (qr(i,k-1).lt.qsmall) then
           qr(i,k-1)=0._r8
           nr(i,k-1)=0._r8
         end if
         if (qi(i,k-1).lt.qsmall) then
           qi(i,k-1)=0._r8
           ni(i,k-1)=0._r8
         end if

         if (qc(i,k-1).lt.qsmall) then
           qc(i,k-1)=0._r8
           nc(i,k-1)=0._r8
         end if

! make sure number concentration is a positive number to avoid 
! taking root of negative

         nr(i,k-1)=max(nr(i,k-1),0._r8)
         ns(i,k-1)=max(ns(i,k-1),0._r8)
         ni(i,k-1)=max(ni(i,k-1),0._r8)
         nc(i,k-1)=max(nc(i,k-1),0._r8)
!!......................................................................
!dkay
  !      frz(i,k) = amin1(1.E-07, frz(i,k))  ! constrain frz 
       end do ! k loop
!dkay          print *,'jb, jt=',jb(i), jt(i)
       end do ! it loop, iteration
300    continue  ! continue if no cloud water
       end do ! i loop
!........................................................................

!        deallocate( &
!         naermod,  &
!         naer2,    &
!         naer2h,    &
!         maerosol)

!        deallocate( &
!         naermod,  &
!         naer2,    &
!         naer2h,    &
!         maerosol)

return
end subroutine mskf_mphy

!##############################################################################


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      FUNCTION mskf_GAMMA(X)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!D    DOUBLE PRECISION FUNCTION DGAMMA(X)
!----------------------------------------------------------------------
!
! THIS ROUTINE CALCULATES THE GAMMA FUNCTION FOR A REAL ARGUMENT X.
!   COMPUTATION IS BASED ON AN ALGORITHM OUTLINED IN REFERENCE 1.
!   THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE GAMMA
!   FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS.  COEFFICIENTS
!   FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPUBLISHED.
!   THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM REFERENCE 2.
!   THE ACCURACY ACHIEVED DEPENDS ON THE ARITHMETIC SYSTEM, THE
!   COMPILER, THE INTRINSIC FUNCTIONS, AND PROPER SELECTION OF THE
!   MACHINE-DEPENDENT CONSTANTS.
!
!
!*******************************************************************
!*******************************************************************
!
! EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
!
! BETA   - RADIX FOR THE FLOATING-POINT REPRESENTATION
! MAXEXP - THE SMALLEST POSITIVE POWER OF BETA THAT OVERFLOWS
! XBIG   - THE LARGEST ARGUMENT FOR WHICH GAMMA(X) IS REPRESENTABLE
!          IN THE MACHINE, I.E., THE SOLUTION TO THE EQUATION
!                  GAMMA(XBIG) = BETA**MAXEXP
! XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER;
!          APPROXIMATELY BETA**MAXEXP
! EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1.0+EPS .GT. 1.0
! XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1/XMININ IS MACHINE REPRESENTABLE
!
!     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
!
!                            BETA       MAXEXP        XBIG
!
! CRAY-1         (S.P.)        2         8191        966.961
! CYBER 180/855
!   UNDER NOS    (S.P.)        2         1070        177.803
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)        2          128        35.040
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)        2         1024        171.624
! IBM 3033       (D.P.)       16           63        57.574
! VAX D-FORMAT   (D.P.)        2          127        34.844
! VAX G-FORMAT   (D.P.)        2         1023        171.489
!
!                            XINF         EPS        XMININ
!
! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
! CYBER 180/855
!   UNDER NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
! VAX D-FORMAT   (D.P.)   1.70D+38     1.39D-17    5.88D-39
! VAX G-FORMAT   (D.P.)   8.98D+307    1.11D-16    1.12D-308
!
!*******************************************************************
!*******************************************************************
!
! ERROR RETURNS
!
!  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
!     WHEN OVERFLOW WOULD OCCUR.  THE COMPUTATION IS BELIEVED
!     TO BE FREE OF UNDERFLOW AND OVERFLOW.
!
!
!  INTRINSIC FUNCTIONS REQUIRED ARE:
!
!     INT, DBLE, EXP, LOG, REAL, SIN
!
!
! REFERENCES:  AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL
!              FUNCTIONS   W. J. CODY, LECTURE NOTES IN MATHEMATICS,
!              506, NUMERICAL ANALYSIS DUNDEE, 1975, G. A. WATSON
!              (ED.), SPRINGER VERLAG, BERLIN, 1976.
!
!              COMPUTER APPROXIMATIONS, HART, ET. AL., WILEY AND
!              SONS, NEW YORK, 1968.
!
!  LATEST MODIFICATION: OCTOBER 12, 1989
!
!  AUTHORS: W. J. CODY AND L. STOLTZ
!           APPLIED MATHEMATICS DIVISION
!           ARGONNE NATIONAL LABORATORY
!           ARGONNE, IL 60439
!
!----------------------------------------------------------------------
      INTEGER I,N
      LOGICAL PARITY

      real(r8) mskf_GAMMA
      REAL(r8) &
!D    DOUBLE PRECISION
         C,CONV,EPS,FACT,HALF,ONE,P,PI,Q,RES,SQRTPI,SUM,TWELVE, &
         TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
      DIMENSION C(7),P(8),Q(8)
!----------------------------------------------------------------------
!  MATHEMATICAL CONSTANTS
!----------------------------------------------------------------------
      DATA ONE,HALF,TWELVE,TWO,ZERO/1.0E0_r8,0.5E0_r8,12.0E0_r8,2.0E0_r8,0.0E0_r8/, &
          SQRTPI/0.9189385332046727417803297E0_r8/, &
          PI/3.1415926535897932384626434E0_r8/
!D    DATA ONE,HALF,TWELVE,TWO,ZERO/1.0D0,0.5D0,12.0D0,2.0D0,0.0D0/,
!D   1     SQRTPI/0.9189385332046727417803297D0/,
!D   2     PI/3.1415926535897932384626434D0/
!----------------------------------------------------------------------
!  MACHINE DEPENDENT PARAMETERS
!----------------------------------------------------------------------
      DATA XBIG,XMININ,EPS/35.040E0_r8,1.18E-38_r8,1.19E-7_r8/, &
          XINF/3.4E38_r8/
!D    DATA XBIG,XMININ,EPS/171.624D0,2.23D-308,2.22D-16/,
!D   1     XINF/1.79D308/
!----------------------------------------------------------------------
!  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
!     APPROXIMATION OVER (1,2).
!----------------------------------------------------------------------
      DATA P/-1.71618513886549492533811E+0_r8,2.47656508055759199108314E+1_r8,&
            -3.79804256470945635097577E+2_r8,6.29331155312818442661052E+2_r8,&
            8.66966202790413211295064E+2_r8,-3.14512729688483675254357E+4_r8,&
            -3.61444134186911729807069E+4_r8,6.64561438202405440627855E+4_r8/
      DATA Q/-3.08402300119738975254353E+1_r8,3.15350626979604161529144E+2_r8,&
           -1.01515636749021914166146E+3_r8,-3.10777167157231109440444E+3_r8,&
             2.25381184209801510330112E+4_r8,4.75584627752788110767815E+3_r8,&
           -1.34659959864969306392456E+5_r8,-1.15132259675553483497211E+5_r8/
!D    DATA P/-1.71618513886549492533811D+0,2.47656508055759199108314D+1,
!D   1       -3.79804256470945635097577D+2,6.29331155312818442661052D+2,
!D   2       8.66966202790413211295064D+2,-3.14512729688483675254357D+4,
!D   3       -3.61444134186911729807069D+4,6.64561438202405440627855D+4/
!D    DATA Q/-3.08402300119738975254353D+1,3.15350626979604161529144D+2,
!D   1      -1.01515636749021914166146D+3,-3.10777167157231109440444D+3,
!D   2        2.25381184209801510330112D+4,4.75584627752788110767815D+3,
!D   3      -1.34659959864969306392456D+5,-1.15132259675553483497211D+5/
!----------------------------------------------------------------------
!  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
!----------------------------------------------------------------------
      DATA C/-1.910444077728E-03_r8,8.4171387781295E-04_r8, &
          -5.952379913043012E-04_r8,7.93650793500350248E-04_r8,&
          -2.777777777777681622553E-03_r8,8.333333333333333331554247E-02_r8,&
           5.7083835261E-03_r8/
!D    DATA C/-1.910444077728D-03,8.4171387781295D-04,
!D   1     -5.952379913043012D-04,7.93650793500350248D-04,
!D   2     -2.777777777777681622553D-03,8.333333333333333331554247D-02,
!D   3      5.7083835261D-03/
!----------------------------------------------------------------------
!  STATEMENT FUNCTIONS FOR CONVERSION BETWEEN INTEGER AND FLOAT
!----------------------------------------------------------------------
      CONV(I) = REAL(I,r8)
!D    CONV(I) = DBLE(I)
      PARITY=.FALSE.
      FACT=ONE
      N=0
      Y=X
      IF(Y.LE.ZERO)THEN
!----------------------------------------------------------------------
!  ARGUMENT IS NEGATIVE
!----------------------------------------------------------------------
        Y=-X
        Y1=AINT(Y)
        RES=Y-Y1
        IF(RES.NE.ZERO)THEN
          IF(Y1.NE.AINT(Y1*HALF)*TWO)PARITY=.TRUE.
          FACT=-PI/SIN(PI*RES)
          Y=Y+ONE
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ENDIF
!----------------------------------------------------------------------
!  ARGUMENT IS POSITIVE
!----------------------------------------------------------------------
      IF(Y.LT.EPS)THEN
!----------------------------------------------------------------------
!  ARGUMENT .LT. EPS
!----------------------------------------------------------------------
        IF(Y.GE.XMININ)THEN
          RES=ONE/Y
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ELSEIF(Y.LT.TWELVE)THEN
        Y1=Y
        IF(Y.LT.ONE)THEN
!----------------------------------------------------------------------
!  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
          Z=Y
          Y=Y+ONE
        ELSE
!----------------------------------------------------------------------
!  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
!----------------------------------------------------------------------
          N=INT(Y)-1
          Y=Y-CONV(N)
          Z=Y-ONE
        ENDIF
!----------------------------------------------------------------------
!  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
!----------------------------------------------------------------------
        XNUM=ZERO
        XDEN=ONE
        DO 260 I=1,8
          XNUM=(XNUM+P(I))*Z
          XDEN=XDEN*Z+Q(I)
  260   CONTINUE
        RES=XNUM/XDEN+ONE
        IF(Y1.LT.Y)THEN
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
          RES=RES/Y1
        ELSEIF(Y1.GT.Y)THEN
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  2.0 .LT. ARGUMENT .LT. 12.0
!----------------------------------------------------------------------
          DO 290 I=1,N
            RES=RES*Y
            Y=Y+ONE
  290     CONTINUE
        ENDIF
      ELSE
!----------------------------------------------------------------------
!  EVALUATE FOR ARGUMENT .GE. 12.0,
!----------------------------------------------------------------------
       IF(Y.LE.XBIG)THEN
          YSQ=Y*Y
          SUM=C(7)
          DO 350 I=1,6
            SUM=SUM/YSQ+C(I)
  350     CONTINUE
          SUM=SUM/Y-Y+SQRTPI
          SUM=SUM+(Y-HALF)*LOG(Y)
          RES=EXP(SUM)
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ENDIF
!----------------------------------------------------------------------
!  FINAL ADJUSTMENTS AND RETURN
!----------------------------------------------------------------------
      IF(PARITY)RES=-RES
      IF(FACT.NE.ONE)RES=FACT/RES
  900 mskf_GAMMA=RES
!D900 DGAMMA = RES
      RETURN
! ---------- LAST LINE OF mskf_GAMMA ----------
      END function mskf_GAMMA

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! error function in single precision
!
!    Copyright(C) 1996 Takuya OOURA (email: ooura@mmm.t.u-tokyo.ac.jp).
!    You may use, copy, modify this code for any purpose and
!    without fee. You may distribute this ORIGINAL package.

      function derf(x)
      implicit real (a - h, o - z)
      real(r8) a,b,x
      dimension a(0 : 64), b(0 : 64)
      integer i,k
      data (a(i), i = 0, 12) / &
         0.00000000005958930743d0, -0.00000000113739022964d0, &
         0.00000001466005199839d0, -0.00000016350354461960d0, &
         0.00000164610044809620d0, -0.00001492559551950604d0, &
         0.00012055331122299265d0, -0.00085483269811296660d0, &
         0.00522397762482322257d0, -0.02686617064507733420d0, &
         0.11283791670954881569d0, -0.37612638903183748117d0, &
         1.12837916709551257377d0 /
      data (a(i), i = 13, 25) / &
         0.00000000002372510631d0, -0.00000000045493253732d0, &
         0.00000000590362766598d0, -0.00000006642090827576d0, &
         0.00000067595634268133d0, -0.00000621188515924000d0, &
         0.00005103883009709690d0, -0.00037015410692956173d0, &
         0.00233307631218880978d0, -0.01254988477182192210d0, &
         0.05657061146827041994d0, -0.21379664776456006580d0, &
         0.84270079294971486929d0 /
      data (a(i), i = 26, 38) / &
         0.00000000000949905026d0, -0.00000000018310229805d0, &
         0.00000000239463074000d0, -0.00000002721444369609d0, &
         0.00000028045522331686d0, -0.00000261830022482897d0, &
         0.00002195455056768781d0, -0.00016358986921372656d0, &
         0.00107052153564110318d0, -0.00608284718113590151d0, &
         0.02986978465246258244d0, -0.13055593046562267625d0, &
         0.67493323603965504676d0 /
      data (a(i), i = 39, 51) / &
         0.00000000000382722073d0, -0.00000000007421598602d0, &
         0.00000000097930574080d0, -0.00000001126008898854d0, &
         0.00000011775134830784d0, -0.00000111992758382650d0, &
         0.00000962023443095201d0, -0.00007404402135070773d0, &
         0.00050689993654144881d0, -0.00307553051439272889d0, &
         0.01668977892553165586d0, -0.08548534594781312114d0, &
         0.56909076642393639985d0 /
      data (a(i), i = 52, 64) / &
         0.00000000000155296588d0, -0.00000000003032205868d0, &
         0.00000000040424830707d0, -0.00000000471135111493d0, &
         0.00000005011915876293d0, -0.00000048722516178974d0, &
         0.00000430683284629395d0, -0.00003445026145385764d0, &
         0.00024879276133931664d0, -0.00162940941748079288d0, &
         0.00988786373932350462d0, -0.05962426839442303805d0, &
         0.49766113250947636708d0 /
     data (b(i), i = 0, 12) / &
         -0.00000000029734388465d0, 0.00000000269776334046d0, &
         -0.00000000640788827665d0, -0.00000001667820132100d0, &
         -0.00000021854388148686d0, 0.00000266246030457984d0, &
         0.00001612722157047886d0, -0.00025616361025506629d0, &
         0.00015380842432375365d0, 0.00815533022524927908d0, &
         -0.01402283663896319337d0, -0.19746892495383021487d0,&
         0.71511720328842845913d0 /
      data (b(i), i = 13, 25) / &
         -0.00000000001951073787d0, -0.00000000032302692214d0, &
         0.00000000522461866919d0, 0.00000000342940918551d0, &
         -0.00000035772874310272d0, 0.00000019999935792654d0, &
         0.00002687044575042908d0, -0.00011843240273775776d0, &
         -0.00080991728956032271d0, 0.00661062970502241174d0, &
         0.00909530922354827295d0, -0.20160072778491013140d0, &
         0.51169696718727644908d0 /
      data (b(i), i = 26, 38) / &
         0.00000000003147682272d0, -0.00000000048465972408d0, &
         0.00000000063675740242d0, 0.00000003377623323271d0, &
         -0.00000015451139637086d0, -0.00000203340624738438d0,&
         0.00001947204525295057d0, 0.00002854147231653228d0, &
         -0.00101565063152200272d0, 0.00271187003520095655d0, &
         0.02328095035422810727d0, -0.16725021123116877197d0, &
         0.32490054966649436974d0 /
      data (b(i), i = 39, 51) / &
         0.00000000002319363370d0, -0.00000000006303206648d0, &
         -0.00000000264888267434d0, 0.00000002050708040581d0, &
         0.00000011371857327578d0, -0.00000211211337219663d0, &
         0.00000368797328322935d0, 0.00009823686253424796d0, &
         -0.00065860243990455368d0, -0.00075285814895230877d0,&
         0.02585434424202960464d0, -0.11637092784486193258d0, &
         0.18267336775296612024d0 /
      data (b(i), i = 52, 64) / &
         -0.00000000000367789363d0, 0.00000000020876046746d0, &
         -0.00000000193319027226d0, -0.00000000435953392472d0, &
         0.00000018006992266137d0, -0.00000078441223763969d0, &
         -0.00000675407647949153d0, 0.00008428418334440096d0, &
         -0.00017604388937031815d0, -0.00239729611435071610d0, &
         0.02064129023876022970d0, -0.06905562880005864105d0, &
         0.09084526782065478489d0 /
      w = abs(x)
      if (w .lt. 2.2d0) then
          t = w * w
          k = int(t)
          t = t - k
          k = k * 13
          y = ((((((((((((a(k) * t + a(k + 1)) * t + &
             a(k + 2)) * t + a(k + 3)) * t + a(k + 4)) * t + &
             a(k + 5)) * t + a(k + 6)) * t + a(k + 7)) * t + &
             a(k + 8)) * t + a(k + 9)) * t + a(k + 10)) * t + &
             a(k + 11)) * t + a(k + 12)) * w
      else if (w .lt. 6.9d0) then
          k = int(w)
          t = w - k
          k = 13 * (k - 2)
          y = (((((((((((b(k) * t + b(k + 1)) * t + &
             b(k + 2)) * t + b(k + 3)) * t + b(k + 4)) * t + &
             b(k + 5)) * t + b(k + 6)) * t + b(k + 7)) * t + &
             b(k + 8)) * t + b(k + 9)) * t + b(k + 10)) * t + &
             b(k + 11)) * t + b(k + 12)
          y = y * y
          y = y * y
          y = y * y
          y = 1 - y * y
      else
          y = 1
      end if
      if (x .lt. 0) y = -y
      derf = y
      end function derf

!-----------------------------------------------------------------------
        real function erfc_num_recipes( x )
!
!   from press et al, numerical recipes, 1990, page 164
!
        implicit none
        real x
        double precision erfc_dbl, dum, t, zz

        zz = abs(x)
        t = 1.0/(1.0 + 0.5*zz)

!       erfc_num_recipes =
!     &   t*exp( -zz*zz - 1.26551223 + t*(1.00002368 + t*(0.37409196 +
!     &   t*(0.09678418 + t*(-0.18628806 + t*(0.27886807 +
!     &                                    t*(-1.13520398 +
!     &   t*(1.48851587 + t*(-0.82215223 + t*0.17087277 )))))))))

        dum =  ( -zz*zz - 1.26551223 + t*(1.00002368 + t*(0.37409196 +   &
          t*(0.09678418 + t*(-0.18628806 + t*(0.27886807 +   &
                                           t*(-1.13520398 +   &
          t*(1.48851587 + t*(-0.82215223 + t*0.17087277 )))))))))

        erfc_dbl = t * exp(dum)
        if (x .lt. 0.0) erfc_dbl = 2.0d0 - erfc_dbl

        erfc_num_recipes = erfc_dbl

        return
        end function erfc_num_recipes

!-----------------------------------------------------------------------
    real function erf_alt( x )

    implicit none

    real,intent(in) :: x

    erf_alt = 1. - erfc_num_recipes(x)

    end function erf_alt
      subroutine mskf_activate(wbar, tair, rhoair,  &
                 na, pmode, nmode, ma, sigman, hygro, rhodry, nact,qs)
!      calculates number, surface, and mass fraction of aerosols activated as
!      CCN
!      calculates flux of cloud droplets, surface area, and aerosol mass into
!      cloud
!      assumes an internal mixture within each of up to pmode multiple aerosol
!      modes
!      a gaussiam spectrum of updrafts can be treated.

!      mks units

!      Abdul-Razzak and Ghan, A parameterization of aerosol activation.
!      2. Multiple aerosol types. J. Geophys. Res., 105, 6837-6844.

!      use physconst, only: rair, epsilo, cpair, rh2o, latvap, gravit,   &
!                                 rhoh2o, mwh2o, r_universal
!ckay      use wv_saturation, only: estblf, epsqs

      implicit none

!      save  ! sep6

!      input

      integer pmode,ptype ! dimension of modes, types in modes
      real(r8) wbar          ! grid cell mean vertical velocity (m/s)
      real(r8) tair          ! air temperature (K)
      real(r8) rhoair        ! air density (kg/m3)
      real(r8) na(pmode)           ! aerosol number concentration (/m3)
      integer nmode      ! number of aerosol modes
      real(r8) ma(pmode)     ! aerosol mass concentration (kg/m3)
      real(r8) rhodry(pmode) ! density of aerosol material
      real(r8) sigman(pmode)  ! geometric standard deviation of aerosol size distribution
      real(r8) hygro(pmode)  ! hygroscopicity of aerosol mode


!      output

      real(r8) nact      ! number fraction of aerosols activated

!      local
#if (defined AIX)
#define ERF erf
#define ERFC erfc
#else
#define ERF derf
#define ERFC derfc
#define ERF_ALT erf_alt
      real(r8) derf,derfc, erf_alt
#endif

      integer, parameter:: nx=200
      integer :: maxmodes

      real(r8) surften       ! surface tension of water w/respect to air (N/m)
      data surften/0.076/
      save surften
      real(r8) p0     ! reference pressure (Pa)
      data p0/1013.25e2/
      save p0

      real(r8) :: volc(naer_cu) ! total aerosol volume  concentration (m3/m3)
      real(r8) tmass ! total aerosol mass concentration (g/cm3)
      real(r8) rm ! number mode radius of aerosol at max supersat (cm)
      real(r8) pres ! pressure (Pa)
      real(r8) path ! mean free path (m)
      real(r8) diff ! diffusivity (m2/s)
      real(r8) conduct ! thermal conductivity (Joule/m/sec/deg)
      real(r8) diff0,conduct0
      real(r8) qs ! water vapor saturation mixing ratio
      real(r8) dqsdt ! change in qs with temperature
      real(r8) dqsdp ! change in qs with pressure
      real(r8) gloc ! thermodynamic function (m2/s)
      real(r8) zeta
      real(r8) :: eta(naer_cu)
      real(r8) :: smc(naer_cu)
      real(r8) lnsmax ! ln(smax)
      real(r8) alpha
      real(r8) gammaloc
      real(r8) beta
      real(r8) sqrtg
      real(r8) alogam
      real(r8) rlo,rhi,xint1,xint2,xint3,xint4
      real(r8) w,wnuc,wb
      real(r8) alw,sqrtalw
      real(r8) smax
      real(r8) x,arg
      real(r8) xmincoeff,xcut,volcut,surfcut
      real(r8) z,z1,z2,wf1,wf2,zf1,zf2,gf1,gf2,gf
      real(r8) :: etafactor1,etafactor2max
      real(r8) :: etafactor2(naer_cu)
      real(r8) es
      integer m,n

      real(r8) :: amcubeloc(naer_cu)
      real(r8) :: lnsmloc(naer_cu)
      maxmodes = naer_cu

      if(maxmodes<pmode)then
!         write(*,*)'maxmodes,pmode in activate =',maxmodes,pmode
!         call endrun('kf_activate')
      endif

      nact=0._r8

      if(nmode.eq.1.and.na(1).lt.1.e-20)return

      if(wbar.le.0.)return

      pres=rair*rhoair*tair
      diff0=0.211e-4*(p0/pres)*(tair/t0)**1.94
      conduct0=(5.69+0.017*(tair-t0))*4.186e2*1.e-5 ! convert to J/m/s/deg
!ckay      es = estblf(tair)
!ckay      qs = epsilo*es/(pres-(1.0_r8 - epsqs)*es)
!        print *,'rh2o=',rh2o
      dqsdt=latvap/(rh2o*tair*tair)*qs
      alpha=gravit*(latvap/(cpair*rh2o*tair*tair)-1./(rair*tair))
      gammaloc=(1+latvap/cpair*dqsdt)/(rhoair*qs)
!     growth coefficent Abdul-Razzak & Ghan 1998 eqn 16
!     should depend on mean radius of mode to account for gas kinetic effects
      gloc=1./(rhoh2o/(diff0*rhoair*qs)                                    &
          +latvap*rhoh2o/(conduct0*tair)*(latvap/(rh2o*tair)-1.))
      sqrtg=sqrt(gloc)
      beta=4.*pi*rhoh2o*gloc*gammaloc
      etafactor2max=1.e10/(alpha*wbar)**1.5 ! this should make eta big if na is very small.

      do m=1,nmode
!         internal mixture of aerosols
          volc(m)=ma(m)/(rhodry(m)) ! only if variable size dist
         if(volc(m).gt.1.e-39_r8.and.na(m).gt.1.e-39_r8)then
            etafactor2(m)=1./(na(m)*beta*sqrtg)  !fixed or variable size dist
!            number mode radius (m)
            amcubeloc(m)=(3.*volc(m)/(4.*pi*exp45logsig(m)*na(m)))  ! only if variable size dist
            smc(m)=smcrit(m) ! only for prescribed size dist

!May30,2014
                 if(hygro(m).gt.1.e-10)then   ! loop only if variable size dist
                    smc(m)=2.*aten*sqrt(aten/(27.*hygro(m)*amcubeloc(m)))
                 else
                   smc(m)=100.
                 endif
         else
            smc(m)=1.
            etafactor2(m)=etafactor2max ! this should make eta big if na is very small.
         endif
         lnsmloc(m)=log(smc(m)) ! only if variable size dist
      enddo

!         single  updraft
         wnuc=wbar
!        write(iulog,*)'uniform updraft =',wnuc

            w=wbar
            alw=alpha*wnuc
            sqrtalw=sqrt(alw)
            zeta=2.*sqrtalw*aten/(3.*sqrtg)
            etafactor1=2.*alw*sqrtalw

            do m=1,nmode
               eta(m)=etafactor1*etafactor2(m)
            enddo

!             print *,' kf_maxsat '
            call mskf_maxsat(zeta,eta,nmode,smc,smax)

            lnsmax=log(smax)
!           print *,'smc,smax=',smc,smax
            xmincoeff=alogaten-2.*third*(lnsmax-alog2)-alog3

            nact=0._r8
            do m=1,nmode
               x=2*(lnsmloc(m)-lnsmax)/(3*sq2*alogsig(m))
!original ghan code
!               nact=nact+0.5*(1.-ERF(x))*na(m)
!++ag replace sg erf with hm derf pre 1.68
!               nact=nact+0.5*(1.-derf(x))*na(m)
!++ag 1.68 new error function
                nact=nact+0.5*(1.-erf(x))*na(m)  
!               nact=nact+0.5*(1.-derf(x))*na(m)  
!       write(*,*)'nact',nact,derf(x),na(m),m
!       write(*,*) 'lnsmloc(m)',lnsmloc(m),lnsmax,alogsig(m)
!                write(*,*) 'wbar=',wbar
            enddo
            nact=nact/rhoair ! convert from #/m3 to #/kg

!      write(*,*)'na(m),qs',na(m),m,qs
!      write(*,*)'nact',nact
!      deallocate( &
!         volc,       &
!         eta,        &
!         smc,        &
!         etafactor2, &
!         amcubeloc,  &
!         lnsmloc     )

      return
      end subroutine mskf_activate

      subroutine mskf_maxsat(zeta,eta,nmode,smc,smax)

!      calculates maximum supersaturation for multiple
!      competing aerosol modes.

!      Abdul-Razzak and Ghan, A parameterization of aerosol activation.
!      2. Multiple aerosol types. J. Geophys. Res., 105, 6837-6844.

      implicit none
!     save ! sep6
      integer nmode ! number of modes
      real(r8) :: smc(:) ! critical supersaturation for number mode radius
      real(r8) zeta
      real(r8) :: eta(:)
      real(r8) smax ! maximum supersaturation
      integer m  ! mode index
      real(r8) sum, g1, g2

      do m=1,nmode
         if(zeta.gt.1.e5*eta(m).or.smc(m)*smc(m).gt.1.e5*eta(m))then
!            weak forcing. essentially none activated
            smax=1.e-20
         else
!            significant activation of this mode. calc activation all modes.
            go to 1
         endif
      enddo

      return

  1   continue

      sum=0
      do m=1,nmode
         if(eta(m).gt.1.e-20)then
            g1=sqrt(zeta/eta(m))
            g1=g1*g1*g1
            g2=smc(m)/sqrt(eta(m)+3*zeta)
            g2=sqrt(g2)
            g2=g2*g2*g2
            sum=sum+(f1(m)*g1+f2(m)*g2)/(smc(m)*smc(m))
!    write(*,*)'f1(m)',f1(m),m
         else
            sum=1.e20
         endif
      enddo

      smax=1./sqrt(sum)

      return

      end subroutine mskf_maxsat

subroutine mskf_nucleati(wbar, tair, pair, qv,  qc,  rhoair, & ! TWG add pair and replace relhum with qv
       na,  naer_all, nuci  &
       , onihf, oniimm, onidep, onimey)

!---------------------------------------------------------------
! Purpose:
!  The parameterization of ice nucleation.
!
! Method: The current method is based on Liu & Penner (2005)
!  It related the ice nucleation with the aerosol number, temperature and the
!  updraft velocity. It includes homogeneous freezing of sulfate, immersion
!  freezing of soot, and Meyers et al. (1992) deposition nucleation
!
! Authors: Xiaohong Liu, 01/2005, modifications by A. Gettelman 2009-2010
!----------------------------------------------------------------
! Input Arguments
!
!    save   ! sep6
  integer  naer_all
  real(r8) :: wbar                ! grid cell mean vertical velocity (m/s)
  real(r8) :: tair                ! temperature (K)
  real(r8) :: qv                  ! water vapor mixing ratio (kg/kg)
  real(r8) :: pair                ! pressure (Pa)

  real(r8) :: qc                  ! liquid water mixing ratio (kg/kg)
  real(r8) :: rhoair              ! air density (kg/m3)
  real(r8) :: na(naer_all)        ! aerosol number concentration (/m3)

!
! Output Arguments
!
  real(r8) :: nuci               ! ice number nucleated (#/kg)
  real(r8) :: onihf              ! nucleated number from homogeneous freezing of so4
  real(r8) :: oniimm             ! nucleated number from immersion freezing
  real(r8) :: onidep             ! nucleated number from deposition nucleation
  real(r8) :: onimey             ! nucleated number from deposition nucleation (meyers: mixed phase)
!
! Local workspace
!
  real(r8)  so4_num                                      ! so4 aerosol number (#/cm^3)
  real(r8)  soot_num                                     ! soot (hydrophilic) aerosol number (#/cm^3)
  real(r8)  dst1_num,dst2_num,dst3_num,dst4_num          ! dust aerosol number (#/cm^3)
  real(r8)  dst_num                                      ! total dust aerosol number (#/cm^3)
  real(r8)  nihf                                         ! nucleated number from homogeneous freezing of so4
  real(r8)  niimm                                        ! nucleated number from immersion freezing
  real(r8)  nidep                                        ! nucleated number from deposition nucleation
  real(r8)  nimey                                        ! nucleated number from deposition nucleation (meyers)
  real(r8)  n1,ni                                        ! nucleated number
  real(r8)  tc,A,B,C,regm                                ! work variable
  real(r8)  esl,esi,deles,qsi,qsl,relhum                 ! work variable
  real(r8)  dst_scale
  real(r8)  subgrid
  real(r8)  dmc,ssmc         ! variables for modal scheme.

    so4_num=0.0_r8
    soot_num=0.0_r8
    dst_num=0.0_r8
    dst1_num = 0.0_r8
    dst2_num = 0.0_r8
    dst3_num = 0.0_r8
    dst4_num = 0.0_r8

!For modal aerosols, assume for the upper troposphere:
! soot = accumulation mode
! sulfate = aiken mode
! dust = coarse mode
! since modal has internal mixtures.

    if(idxsul .gt. 0) then
      so4_num=na(idxsul)*1.0e-6_r8 ! #/cm^3
    end if

    if(idxbcphi .gt. 0) then
      soot_num=na(idxbcphi)*1.0e-6_r8 !#/cm^3
    end if

    if(idxdst1 .gt. 0) then
       dst1_num=na(idxdst1)*1.0e-6_r8 !#/cm^3
    end if

    if(idxdst2 .gt. 0) then
       dst2_num=na(idxdst2)*1.0e-6_r8 !#/cm^3
    end if

    if(idxdst3 .gt. 0) then
       dst3_num=na(idxdst3)*1.0e-6_r8 !#/cm^3
    end if

    if(idxdst4 .gt. 0) then
       dst4_num=na(idxdst4)*1.0e-6_r8 !#/cm^3
    end if

    dst_num =dst1_num+dst2_num+dst3_num+dst4_num
! no soot nucleation for now.
   ! soot_num=0.0_r8

    ni=0._r8
    tc=tair-273.15_r8

!TWG calculate Ice saturation ratio
    esi = 611.2*exp(21.87*(tair-273.16)/(tair-7.66))
    esl = 611.2*exp(17.67*(tair-273.16)/(243.5+tair-273.16))
    qsi = 0.622*esi/(pair-esi)
    qsl = 0.622*esl/(pair-esl)
    deles = qv/qsi
    relhum = qv/qsl

    ! initialize
    niimm=0._r8
    nidep=0._r8
    nihf=0._r8

    if(so4_num.ge.1.0e-10_r8 .and. (soot_num+dst_num).ge.1.0e-10_r8 ) then

      subgrid = 1.0_r8

!          print *,'nucleiate wbar=', wbar
     if((wbar.lt.4.0_r8) .and. (tc.le.-35.0_r8) .and.((deles*subgrid).ge.1.0_r8)) then !TWG Make RHi consistent
!< regm => T in Eq.10 of Liu et al., J. Climate, 2007>
       print*,'Aerosol Ice Nucleation is Doing Something'
       A = -1.4938_r8 * log(soot_num+dst_num) + 12.884_r8
       B = -10.41_r8  * log(soot_num+dst_num) - 67.69_r8
       regm = A * log(wbar) + B

!         print *,'before bunch of hetero'
       if(tc.gt.regm) then    ! heterogeneous nucleation only
         if(tc.lt.-40._r8 .and. wbar.gt.1._r8) then ! exclude T<-40 & W>1m/s from hetero. nucleation
           call mskf_hf(tc,wbar,relhum,subgrid,so4_num,nihf)
           niimm=0._r8
           nidep=0._r8
           n1=nihf
         else
           call mskf_hetero(tc,wbar,soot_num+dst_num,niimm,nidep)
           nihf=0._r8
           n1=niimm+nidep
         endif
       elseif (tc.lt.regm-5._r8) then ! homogeneous nucleation only
         call mskf_hf(tc,wbar,relhum,subgrid,so4_num,nihf)
         niimm=0._r8
         nidep=0._r8
         n1=nihf
       else        ! transition between homogeneous and heterogeneous: interpolate in-between
         if(tc.lt.-40._r8 .and. wbar.gt.1._r8) then ! exclude T<-40 & W>1m/s from hetero. nucleation
           call mskf_hf(tc,wbar,relhum,subgrid,so4_num,nihf)
           niimm=0._r8
           nidep=0._r8
           n1=nihf
         else
           call mskf_hf(regm-5._r8,wbar,relhum,subgrid,so4_num,nihf)
           call mskf_hetero(regm,wbar,soot_num+dst_num,niimm,nidep)
           if(nihf.le.(niimm+nidep)) then
             n1=nihf
           else
              n1=(niimm+nidep)*((niimm+nidep)/nihf)**((tc-regm)/5._r8)
           endif
         endif
       endif

       ni=n1

    endif
    endif
1100  continue

! deposition/condensation nucleation in mixed clouds (-37<T<0C) (Meyers, 1992)
!<Eq.12 of Liu et al., J. Climate, 2007
! Nid(L-1)*1.e-3 => Nid(m-3)
! Question:  RHi=RHw*esl/esi

    if(tc.lt.0._r8 .and. tc.gt.-37._r8 .and. qc.gt.1.e-12_r8) then
!     if(tc.lt.0._r8 .and. tc.gt.-37._r8) then  ! TWG  remove cloud water constraint
!      esl = kf_polysvp(tair,0)     ! over water in mixed clouds
!      esi = kf_polysvp(tair,1)     ! over ice
!songxl      deles = (esl - esi)
!      deles = (relhum*esl - esi)
       if (deles.gt.1.5) THEN
             deles = 1.5
        end if
      nimey=1.e-3_r8*exp(12.96_r8*(deles-1.0_r8) - 0.639_r8) ! TWG fix Meyers formulation
    else
      nimey=0._r8
    endif

    nuci=ni+nimey
    if(nuci.gt.9999._r8.or.nuci.lt.0._r8) then
       write(*, *) 'incorrect ice nucleation number'
       write(*, *) ni, tair, relhum, wbar, nihf, niimm,nidep,deles,esi,dst2_num,dst3_num,dst4_num
       nuci=0._r8
       FATAL_ERROR('Incorrect Ice Nucleation Number, diags')
    endif

    nuci=nuci*1.e+6_r8/rhoair    ! change unit from #/cm3 to #/kg
    onimey=nimey*1.e+6_r8/rhoair
    onidep=nidep*1.e+6_r8/rhoair
    oniimm=niimm*1.e+6_r8/rhoair
    onihf=nihf*1.e+6_r8/rhoair

!     print *,'inputs=',wbar, tair, relhum,  qc,  rhoair, &
!      na,  naer_all, nuci,onimey,onidep,oniimm,onihf 
!     print *,'na,tari,nuci.. =', na,tair,nuci,onimey,onidep,oniimm,onihf
  return
  end subroutine mskf_nucleati

  subroutine mskf_hetero(T,ww,Ns,Nis,Nid)

    real(r8) :: T, ww, Ns
    real(r8) :: Nis, Nid

    real(r8) A11,A12,A21,A22,B11,B12,B21,B22
    real(r8) A,B,C

!    save    ! spe6
!---------------------------------------------------------------------
! parameters

      A11 = 0.0263_r8
      A12 = -0.0185_r8
      A21 = 2.758_r8
      A22 = 1.3221_r8
      B11 = -0.008_r8
      B12 = -0.0468_r8
      B21 = -0.2667_r8
      B22 = -1.4588_r8
!<Eq.11 of Liu et al., J. Climate, 2007>
!     ice from immersion nucleation (cm-3)

      B = (A11+B11*log(Ns)) * log(ww) + (A12+B12*log(Ns))
      C =  A21+B21*log(Ns)

      Nis = exp(A22) * Ns**B22 * exp(B*T) * ww**C
      Nis = min(Nis,Ns)

      Nid = 0.0_r8    ! don't include deposition nucleation for cirrus clouds when T<-37C

      return
  end subroutine mskf_hetero

 subroutine mskf_hf(T,ww,RH,subgrid,Na,Ni)

      real(r8) :: T, ww, RH, subgrid, Na
      real(r8), intent(out) :: Ni

      real(r8)    A1_fast,A21_fast,A22_fast,B1_fast,B21_fast,B22_fast
      real(r8)    A2_fast,B2_fast
      real(r8)    C1_fast,C2_fast,k1_fast,k2_fast
      real(r8)    A1_slow,A2_slow,B1_slow,B2_slow,B3_slow
      real(r8)    C1_slow,C2_slow,k1_slow,k2_slow
      real(r8)    regm
      real(r8)    A,B,C
      real(r8)    RHw

!     save   ! sep6
!---------------------------------------------------------------------
!<Table 1 of  Liu et al., J. Climate, 2007>
! parameters

      A1_fast  =0.0231_r8
      A21_fast =-1.6387_r8  !(T>-64 deg)
      A22_fast =-6.045_r8   !(T<=-64 deg)
      B1_fast  =-0.008_r8
      B21_fast =-0.042_r8   !(T>-64 deg)
      B22_fast =-0.112_r8   !(T<=-64 deg)
      C1_fast  =0.0739_r8
      C2_fast  =1.2372_r8

      A1_slow  =-0.3949_r8
      A2_slow  =1.282_r8
      B1_slow  =-0.0156_r8
      B2_slow  =0.0111_r8
      B3_slow  =0.0217_r8
      C1_slow  =0.120_r8
      C2_slow  =2.312_r8

      Ni = 0.0_r8

!----------------------------
!<Eq.6 of Liu et al., J. Climate, 2007 
! w~m/s, T~degree C, RHw~% => RHw*0.01~fraction  >
!RHw xiaohong's parameter
      A = 6.0e-4_r8*log(ww)+6.6e-3_r8
      B = 6.0e-2_r8*log(ww)+1.052_r8
      C = 1.68_r8  *log(ww)+129.35_r8
      RHw=(A*T*T+B*T+C)*0.01_r8

      if((T.le.-37.0_r8) .and. ((RH*subgrid).ge.RHw)) then

!<Eq.9 of Liu et al., J. Climate, 2007>
        regm = 6.07_r8*log(ww)-55.0_r8

        if(T.ge.regm) then    ! fast-growth regime

          if(T.gt.-64.0_r8) then
            A2_fast=A21_fast
            B2_fast=B21_fast
          else
            A2_fast=A22_fast
            B2_fast=B22_fast
          endif
!<Eq.7 of Liu et al., J. Climate, 2007> 
          k1_fast = exp(A2_fast + B2_fast*T + C2_fast*log(ww))
          k2_fast = A1_fast+B1_fast*T+C1_fast*log(ww)

          Ni = k1_fast*Na**(k2_fast)
          Ni = min(Ni,Na)
        else       ! slow-growth regime
!<Eq.7 of Liu et al., J. Climate, 2007>
          k1_slow = exp(A2_slow + (B2_slow+B3_slow*log(ww))*T + C2_slow*log(ww))
          k2_slow = A1_slow+B1_slow*T+C1_slow*log(ww)

          Ni = k1_slow*Na**(k2_slow)
          Ni = min(Ni,Na)
        endif
      end if

      return
  end subroutine mskf_hf

      function mskf_polysvp (T,type)
!  Compute saturation vapor pressure by using
! function from Goff and Gatch (1946)

!  Polysvp returned in units of pa.
!  T is input in units of K.
!  type refers to saturation with respect to liquid (0) or ice (1)

      real(r8) dum

      real(r8) T,mskf_polysvp

      integer type

! ice

      if (type.eq.1) then

! Goff Gatch equation (good down to -100 C)

         mskf_polysvp = 10._r8**(-9.09718_r8*(273.16_r8/t-1._r8)-3.56654_r8* &
          log10(273.16_r8/t)+0.876793_r8*(1._r8-t/273.16_r8)+ &
          log10(6.1071_r8))*100._r8

      end if


! Goff Gatch equation, uncertain below -70 C

      if (type.eq.0) then
         mskf_polysvp = 10._r8**(-7.90298_r8*(373.16_r8/t-1._r8)+ &
             5.02808_r8*log10(373.16_r8/t)- &
             1.3816e-7_r8*(10._r8**(11.344_r8*(1._r8-t/373.16_r8))-1._r8)+ &
             8.1328e-3_r8*(10._r8**(-3.49149_r8*(373.16_r8/t-1._r8))-1._r8)+ &
             log10(1013.246_r8))*100._r8
         end if


      end function mskf_polysvp

!=================================================================================================================
 end module module_cu_mskf_cmp
!=================================================================================================================
