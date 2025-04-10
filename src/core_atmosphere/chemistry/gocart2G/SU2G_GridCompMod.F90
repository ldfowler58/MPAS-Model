! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================================
 module SU2G_GridCompMod
 use mpas_kind_types,only: RKIND,StrKIND
 use mpas_derived_types,only: MPAS_LOG_CRIT
 use mpas_log

 use GA_EnvironmentMod
 use GOCART2G_AeroGeneric,only: findKlid
 use GOCART2G_MieMod_smiol
 use GOCART2G_Process,only: Chem_BiomassDiurnal,        &
                            Chem_Settling,              &
                            DMSemission,                &
                            SulfateChemDriver,          &
                            SulfateDistributeEMissions, &
                            SulfateUpdateOxidants,      &
                            SU_Compute_Diags,           &
                            SU_Wet_Removal

 use MAPL,only: MAPL_PackTime

 use SU2G_instance,only: nbins,particle_radius_microns,particle_density,fscav,molecular_weight, &
                         fnum,rhFlag,pressure_lid_in_hPa,sigma,SO4_anthropogenic_fraction,      &
                         aircraft_fuel_emission_factor,aviation_vertical_layers
 use SU2G_StateSpecs,only: SU2G_State


 implicit none
 private


!--- constants (these parameters need to be accessed from MPAS physics instead of redefined here):
 real(kind=RKIND),parameter:: cpd      = 1003.0_RKIND
 real(kind=RKIND),parameter:: grav     = 9.80616_RKIND
 real(kind=RKIND),parameter:: pi       = 3.141592653589793_RKIND
 real(kind=RKIND),parameter:: karman   = 0.4_RKIND
 real(kind=RKIND),parameter:: airmw    = 28.97_RKIND
 real(kind=RKIND),parameter:: radTodeg = 180._RKIND/pi
 real(kind=RKIND),parameter:: Avogadro = 6.02214076e23
 real(kind=RKIND),parameter:: undefval = 1.0e15 


!--- relative position of sulfate tracers
 integer,parameter:: nDMS = 1, &
                     nSO2 = 2, &
                     nSO4 = 3, &
                     nMSA = 4


!--- molecular weights of sulfate species (grams):
 real(kind=RKIND),parameter:: fMassSulfur = 32._RKIND, &
                              fMassSO2    = 64._RKIND, &
                              fMassSO4    = 96._RKIND, &
                              fMassDMS    = 62._RKIND, &
                              fMassMSA    = 96._RKIND


!--- land/ocean/sea-ice mask (these parameters needs to be accessed from MPAS phys instead of redefined here):
 real(kind=RKIND),parameter:: OCEAN   = 2._RKIND, &
                              LAND    = 1._RKIND, &
                              SEA_ICE = 1._RKIND


!--- types needed to define SU2G:
 type:: ThreadWorkspace
    logical:: firstRun = .true.
    logical:: recycle_H2O2 = .false.

    integer:: nVolc = 0
    integer:: nPts = -1
    integer:: nymd_oxidants = -1 ! update the oxidant files?
    integer:: nymd_last = -1     ! previous nymd. updated daily.
    integer,dimension(:),allocatable:: pstart,pend
    integer,dimension(:),allocatable:: vStart,vEnd

    real(kind=RKIND),dimension(:),allocatable:: vLat,vLon,vSO2,vElev,vCloud
    real(kind=RKIND),dimension(:),allocatable:: pLat,pLon,pBase,pTop,pEmis
 end type ThreadWorkspace

 type,extends(GA_Environment),public:: SU2G_GridComp
    logical:: diurnal_bb                              ! diurnal biomass burning
    integer:: myDOW = -1                              ! day of the week: Sun=1, Mon=2,...,Sat=7

    real(kind=RKIND):: eAircraftFuel                  ! aircraft emission factor: go from kg fuel to kg SO2
    real(kind=RKIND):: aviation_layers(4)             ! heights of the LTO, CDS and CRS layers
    real(kind=RKIND):: fSO4anth                       ! fraction of anthropogenic emissions that are SO4
    real(kind=RKIND),dimension(:),allocatable:: sigma ! sigma of lognormal number distribution

    !special handling for volcanic emissions:
    character(len=strKIND):: volcano_srcfilen

    !workspace for point emissions:
    character(len=StrKIND):: point_emissions_srcfilen ! filename for pointwise emissions
    logical:: doing_point_emissions = .false.
    type(ThreadWorkspace),dimension(:),allocatable:: workspaces

    contains
       procedure:: emissions_GridComp => emissions_SU2G_GridComp
       procedure:: load_GridComp      => load_SU2G_GridComp
       procedure:: processes_GridComp => processes_SU2G_GridComp
 end type SU2G_GridComp

 type wrap_
    type(SU2G_GridComp),pointer:: PTR !=> null()
 end type wrap_


 contains


!==================================================================================================================
 subroutine load_SU2G_GridComp(self,kts,kte)
!==================================================================================================================

!--- input arguments:
 integer,intent(in):: kts,kte

!--- inout arguments:
 class(SU2G_GridComp),intent(inout) :: self

!--- local variables:
 integer:: n

!------------------------------------------------------------------------------------------------------------------
!call mpas_log_write('--- enter subroutine load_SU2G_GridComp:')


!--- initialization of vertical index:
 self%klid = kts
 self%km = kte-kts+1


!--- initialization using parameters defined in SU2G_instance_SU:
 self%diurnal_bb = .false.

 call self%load_from_config(nbins,particle_radius_microns,particle_density,fscav,molecular_weight,fnum, &
                            rhFlag,pressure_lid_in_hPa)

 if(.not.allocated(self%sigma)) allocate(self%sigma(self%nbins))
 self%eAirCraftFuel = aircraft_fuel_emission_factor
 self%fSO4anth      = SO4_anthropogenic_fraction
 do n = 1,self%nbins
    self%sigma(n) = sigma(n)
 enddo
 do n = 1,4
    self%aviation_layers(n) = aviation_vertical_layers(n)
 enddo

!call mpas_log_write('--- nbins = $i',intArgs=(/self%nbins/))
!call mpas_log_write('--- radius,rhop,fscav,molwght,fnum:')
!do n = 1,self%nbins
!   call mpas_log_write('$i $r $r $r $r $r $r',intArgs=(/n/),realArgs=(/self%radius(n),self%rhop(n), &
!                       self%fscav(n),self%molwght(n),self%fnum(n),self%sigma(n)/))
!enddo


!call mpas_log_write('--- end subroutine load_SU2G_GridCOMP.')

 end subroutine load_SU2G_GridComp

!==================================================================================================================
 subroutine emissions_SU2G_GridComp(self_params,self,its,ite,jts,jte,kts,kte,iyr,imm,idd,ihr,imn,isc)
!==================================================================================================================

!--- input arguments:
 integer,intent(in):: its,ite,jts,jte,kts,kte
 integer,intent(in):: iyr,imm,idd,ihr,imn,isc

 class(SU2G_GridComp),intent(in):: self_params

!--- inout arguments:
 class(SU2G_State),intent(inout):: self

!--- local variables and arrays:
 integer:: i,j,k
 integer:: nymd,nhms
 integer:: istat
 integer:: nVolc

 real(kind=RKIND),dimension(:,:),allocatable:: so2biomass_src_

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine emissions_SU2G_GridComp:')


!--- extract nymd and nhms from clock:
 call MAPL_PackTime(nymd,iyr,imm,idd)
 call MAPL_PackTime(nhms,ihr,imn,isc)
!call mpas_log_write('--- nymd = $i',intArgs=(/nymd/))
!call mpas_log_write('--- nhms = $i',intArgs=(/nhms/))


!--- for now, we do not support volcanic emissions and point emissions:
 nVolc = 0


!--- apply diurnal cycle to biomass burning if needed:
 if(self_params%diurnal_bb ) then
!   call mpas_log_write('--- enter subroutine Chem_BiomassDiurnal:')
    so2biomass_src_ = self%su_biomass
    call Chem_BiomassDiurnal( &
       cdt  = self_params%cdt,    &
       nhms = nhms,               &
       eout = self%su_biomass,    &
       ein  = so2biomass_src_,    &
       lons = self%lons*radTodeg, &
       lats = self%lats*radTodeg  &
                            )
!   call mpas_log_write('--- end subroutine Chem_BiomassDiurnal.')
 endif


!--- apply emissions to SO2 and SO4:
!call mpas_log_write('--- enter subroutine SulfateDistributeEmissions:')
 istat = 0
 call SulfateDistributeEmissions( &
    km                = self_params%km,              &
    cdt               = self_params%cdt,             &
    fSO4ant           = self_params%fSO4anth,        &
    eAircraftFuel     = self_params%eAircraftFuel,   &
    aviation_layers   = self_params%aviation_layers, &
    nbins             = self_params%nbins,           &
    grav              = grav,                        &
    nymd              = nymd,                        &
    nhms              = nhms,                        &
    fMassSO2          = fMassSO2,                    &
    fMassSO4          = fMassSO4,                    &
    nSO2              = nSO2,                        &
    nSO4              = nSO4,                        &
    so2biomass_src    = self%su_biomass,             &
    so2anthro_l1_src  = self%su_anthrol1,            &
    so2anthro_l2_src  = self%su_anthrol2,            &
    so2ship_src       = self%su_shipso2,             &
    so4ship_src       = self%su_shipso4,             &
    aviation_lto_src  = self%su_aviation_lto,        &
    aviation_cds_src  = self%su_aviation_cds,        &
    aviation_crs_src  = self%su_aviation_crs,        &
    aircraft_fuel_src = self%su_aircraft,            &
    SO2               = self%so2,                    &
    SO4               = self%so4,                    &
    oro               = self%lwi,                    &
    u10m              = self%u10m,                   &
    v10m              = self%v10m,                   &
    hghte             = self%zle,                    &
    pblh              = self%zpbl,                   &
    tmpu              = self%t,                      &
    rhoa              = self%airdens,                &
    delp              = self%delp,                   &
    nVolc             = nVolc,                       &
    su_emis           = self%suem,                   &
    su_so4eman        = self%so4eman,                &
    su_so2eman        = self%so2eman,                &
    su_so2embb        = self%so2embb,                &
    rc                = istat                        &
                                )
 if(istat /=0) then
    call mpas_log_write('--- SU2G_bc_GridComp: error in subroutine SulfateDistributeEmissions.', &
                        messageType=MPAS_LOG_CRIT)
 else
!   call mpas_log_write('--- end subroutine SulfateDistributionEmissions.')
 endif


 if(associated(self%DMS)) then
!   call mpas_log_write('--- enter subroutine DMSemission:')
    istat = 0
    call DMSemission( &
       km        = self_params%km,  &
       cdt       = self_params%cdt, &
       nDMS      = nDMS,            &
       fMassDMS  = fMassDMS,        &
       grav      = grav,            &
       tmpu      = self%t,          &
       u10m      = self%u10m,       &
       v10m      = self%v10m,       &
       oro       = self%lwi,        &
       delp      = self%delp,       &
       dmso_conc = self%su_dmso,    &
       dms       = self%dms,        &
       su_emis   = self%suem,       &
       rc        = istat            &
                    )
    if(istat /=0) then
       call mpas_log_write('--- SU2G_bc_GridComp: error in subroutine DMSEmission.', &
                           messageType=MPAS_LOG_CRIT)
    else
!      call mpas_log_write('--- end subroutine DMSEmission.')
    endif
 endif


 call mpas_log_write('--- end subroutine emissions_SU2G_GridComp.')

 end subroutine emissions_SU2G_GridComp

!==================================================================================================================
 subroutine processes_SU2G_GridComp(self_params,self,its,ite,jts,jte,kts,kte,iyr,imm,idd,ihr,imn,isc)
!==================================================================================================================

!--- input arguments:
 integer,intent(in):: its,ite,jts,jte,kts,kte
 integer,intent(in):: iyr,imm,idd,ihr,imn,isc

!--- inout arguments:
 class(SU2G_GridComp),intent(inout):: self_params
 class(SU2G_State),intent(inout):: self

!--- local variables:
 logical:: KIN
 logical:: correctionMaring
 logical:: recycle_h2o2
 integer:: istat
 integer:: i,i1,i2,j,j1,j2,k,ibin,km,n
 integer:: nymd,nhms
 integer:: nymd_last
 integer:: nw_profile,nw_vertint

!real(kind=RKIND),dimension(:,:,:),pointer:: rh20,rh80
 real(kind=RKIND),dimension(:,:),allocatable:: drydepf
 real(kind=RKIND),dimension(:,:,:),allocatable:: h2o2_init
 real(kind=RKIND),dimension(:,:,:),allocatable:: xh2o2,xoh,xno3
 real(kind=RKIND),dimension(:,:,:),allocatable:: dms_init,so2_init,so4_init,msa_init
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: qsu2G

 real(kind=RKIND),allocatable,dimension(:,:,:),target:: rh20,rh80

!------------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine processes_SU2G_GridComp:')

 nw_profile = size(self_params%wavelengths_profile)
 nw_vertint = size(self_params%wavelengths_vertint)


!--- extract nymd and nhms from clock:
 call MAPL_PackTime(nymd,iyr,imm,idd)
 call MAPL_PackTime(nhms,ihr,imn,isc)
!call mpas_log_write('--- nymd = $i',intArgs=(/nymd/))
!call mpas_log_write('--- nhms = $i',intArgs=(/nhms/))


!--- SU2G oxidants:
 if(.not.allocated(xoh)      ) allocate(xoh(its:ite,jts:jte,kts:kte)      )
 if(.not.allocated(xh2o2)    ) allocate(xh2o2(its:ite,jts:jte,kts:kte)    )
 if(.not.allocated(xno3)     ) allocate(xno3(its:ite,jts:jte,kts:kte)     )
 if(.not.allocated(h2o2_init)) allocate(h2o2_init(its:ite,jts:jte,kts:kte))
 xoh(:,:,:)       = 0._RKIND
 xh2o2(:,:,:)     = 0._RKIND
 xno3(:,:,:)      = 0._RKIND
 h2o2_init(:,:,:) = 0._RKIND

 if(.not.allocated(dms_init)) allocate(dms_init(its:ite,jts:jte,kts:kte))
 if(.not.allocated(so2_init)) allocate(so2_init(its:ite,jts:jte,kts:kte))
 if(.not.allocated(so4_init)) allocate(so4_init(its:ite,jts:jte,kts:kte))
 if(.not.allocated(msa_init)) allocate(msa_init(its:ite,jts:jte,kts:kte))
 do j = jts,jte
    do i = its,ite
       do k = kts,kte
          dms_init(i,j,k) = self%dms(i,j,k)
          so2_init(i,j,k) = self%so2(i,j,k)
          so4_init(i,j,k) = self%so4(i,j,k)
          msa_init(i,j,k) = self%msa(i,j,k)
       enddo
    enddo
 enddo

!call mpas_log_write('--- enter subroutine SulfateUpdateOxidants:')
 recycle_h2o2 = .true.
 nymd_last    = -1
 istat = 0
 call SulfateUpdateOxidants( &
              nymd_current   = nymd            , &
              nymd_last      = nymd_last       , &
              nhms_current   = nhms            , &
              km             = self_params%km  , &
              cdt            = self_params%cdt , &
              undefval       = undefval        , &
              radToDeg       = radToDeg        , &
              nAvogadro      = Avogadro/1000.  , &
              pi             = pi              , &
              airMolWght     = airmw           , &
              lonRad         = self%lons       , &
              latRad         = self%lats       , &
              rhoa           = self%airdens    , &
              oh_clim        = self%su_oh      , & !climatological OH source.
              no3_clim       = self%su_no3     , & !climatological NO3 source.
              h2o2_clim      = self%su_h2o2    , & !climatological H2O2 source.
              xoh            = xoh             , &
              xno3           = xno3            , &
              xh2o2          = xh2o2           , &
              recycle_h2o2   = recycle_h2o2    , &
              rc             = istat             &
                           )
 if(istat /=0) then
    call mpas_log_write('--- SU2G_GridComp: error in subroutine SulfateUpdateOxidants.', &
                        messageType=MPAS_LOG_CRIT)
 else
!   call mpas_log_write('--- end subroutine SulfateChemOxidants.')
 endif


!--- SU2G settling:
!call mpas_log_write('--- enter subroutine Chem_Settling:')
 if(.not.allocated(qsu2G)) allocate(qsu2G(its:ite,jts:jte,kts:kte,self_params%nbins))
 do j = jts,jte
    do i = its,ite
       do k = kts,kte
          qsu2G(i,j,k,nDMS) = self%dms(i,j,k)
          qsu2G(i,j,k,nSO2) = self%so2(i,j,k)
          qsu2G(i,j,k,nSO4) = self%so4(i,j,k)
          qsu2G(i,j,k,nMSA) = self%msa(i,j,k)
       enddo
    enddo
 enddo

 istat = 0
 do ibin = 1, self_params%nbins
    !if radius == 0, then we're dealing with a gas which has no settling losses:
    if(self_params%radius(ibin) == 0.0) then
       if(associated(self%susd)) self%susd(:,:,ibin) = 0.0
       cycle
    endif
    call Chem_Settling( &
              km        = self_params%km                 , &
              klid      = self_params%klid               , &
              flag      = self_params%rhFlag             , &
              cdt       = self_params%cdt                , &
              radiusInp = self_params%radius(ibin)*1.e-6 , &
              rhopInp   = self_params%rhop(ibin)         , &
              bin       = ibin                           , &
              grav      = grav                           , &
              int_qa    = qsu2G(:,:,:,ibin)              , &
              tmpu      = self%t                         , &
              rhoa      = self%airdens                   , &
              rh        = self%rh2                       , &
              hghte     = self%zle                       , &
              delp      = self%delp                      , &
              fluxout   = self%susd                      , &
              rc        = istat                            &
                      )
 enddo
 if(istat /=0) then
    call mpas_log_write('--- SU2G_GridComp: error in subroutine Chem_Settling:', &
                        messageType=MPAS_LOG_CRIT)
 else
    do j = jts,jte
       do i = its,ite
          do k = kts,kte
             self%dms(i,j,k) = qsu2G(i,j,k,nDMS)
             self%so2(i,j,k) = qsu2G(i,j,k,nSO2)
             self%so4(i,j,k) = qsu2G(i,j,k,nSO4)
             self%msa(i,j,k) = qsu2G(i,j,k,nMSA)
          enddo
       enddo
    enddo
    if(allocated(qsu2G)) deallocate(qsu2G)
!   call mpas_log_write('--- end subroutine Chem_Settling:')
 endif


!--- SU2G chem driver:
!call mpas_log_write('--- enter subroutine SulfateChemDriver:')
 istat = 0
 if(associated(self%suvdep)) self%suvdep(:,:) = 0._RKIND
 if(.not.allocated(drydepf)) allocate(drydepf(its:ite,jts:jte))
 call SulfateChemDriver( &
              km             = self_params%km   , &
              klid           = self_params%klid , &
              cdt            = self_params%cdt  , &
              pi             = pi               , &
              radToDeg       = radToDeg         , &
              von_karman     = karman           , &
              airMolWght     = airmw            , &
              nAvogadro      = Avogadro/1000.   , &
              cpd            = cpd              , &
              grav           = grav             , &
              fMassMSA       = fMassMSA         , &
              fMassDMS       = fMassDMS         , &
              fMassSO2       = fMassSO2         , &
              fMassSO4       = fMassSO4         , &
              nymd           = nymd             , &
              nhms           = nhms             , &
              lonRad         = self%lons        , &
              latRad         = self%lats        , &
              dms            = self%dms         , &
              so2            = self%so2         , &
              so4            = self%so4         , &
              msa            = self%msa         , &
              nDMS           = nDMS             , &
              nSO2           = nSO2             , &
              nSO4           = nSO4             , &
              nMSA           = nMSA             , &
              xoh            = xoh              , &
              xno3           = xno3             , &
              xh2o2          = xh2o2            , &
              h2o2_init      = h2o2_init        , &
              delp           = self%delp        , &
              tmpu           = self%t           , &
              cloud          = self%fcld        , &
              rhoa           = self%airdens     , &
              hghte          = self%zle         , &
              ustar          = self%ustar       , &
              shflux         = self%sh          , &
              oro            = self%lwi         , &
              pblh           = self%zpbl        , &
              z0h            = self%z0h         , &
              su_dep         = self%sudp        , &
              su_pso2        = self%supso2      , &
              su_pmsa        = self%supmsa      , &
              su_pso4        = self%supso4      , &
              su_pso4g       = self%supso4g     , &
              su_pso4aq      = self%supso4aq    , &
              pso2           = self%pso2        , &
              pmsa           = self%pmsa        , &
              pso4           = self%pso4        , &
              pso4g          = self%pso4g       , &
              pso4aq         = self%pso4aq      , &
              drydepositionfrequency = drydepf  , &
              rc             = istat              &
                       )
 if(associated(self%suvdep)) then
    self%suvdep(:,:) = drydepf(:,:)*self%delz(:,:,self_params%km)
 endif
 if(istat /=0) then
    call mpas_log_write('--- SU2G_GridComp: error in subroutine SulfateChemDriver.', &
                        messageType=MPAS_LOG_CRIT)
 else
    if(allocated(drydepf)) deallocate(drydepf)
!   call mpas_log_write('--- end subroutine SulfateChemDriver.')
 endif


!--- SU2G wet removal:
!call mpas_log_write('--- enter subroutine SU_Wet_Removal:')
 istat = 0
 KIN = .true.
 call SU_Wet_Removal( &
              km              = self_params%km    , &
              nbins           = self_params%nbins , &
              klid            = self_params%klid  , &
              cdt             = self_params%cdt   , &
              kin             = KIN               , & 
              grav            = grav              , &
              airMolWght      = airmw             , &
              delp            = self%delp         , &
              fMassSO4        = fMassSO4          , &
              fMassSO2        = fMassSO2          , &
              h2o2_int        = h2o2_init         , &
              ple             = self%ple          , &
              rhoa            = self%airdens      , &
              precc           = self%cn_prcp      , &
              precl           = self%ncn_prcp     , &
              pfllsan         = self%pfl_lsan     , &
              pfilsan         = self%pfi_lsan     , &
              tmpu            = self%t            , &
              nDMS            = nDMS              , &
              nSO2            = nSO2              , &
              nSO4            = nSO4              , &
              nMSA            = nMSA              , &
              dms             = self%dms          , &
              so2             = self%so2          , &
              so4             = self%so4          , &
              msa             = self%msa          , &
              fluxout         = self%suwt         , &
              pso4_colflux    = self%supso4       , &
              pso4wet_colflux = self%supso4wt     , &
              pso4            = self%pso4         , &
              pso4wet         = self%pso4wet      , &
              rc              = istat               &
                    )
 if(istat /=0) then
    call mpas_log_write('--- SU2G_GridComp: error in subroutine SU_Wet_Removal.', &
                        messageType=MPAS_LOG_CRIT)
 else
!   call mpas_log_write('--- end subroutine SU_Wet_Removal.')
 endif


!--- SU2G diagnostics:
!Certain variables are multiplied by 1.0e-9 to convert from nanometers to meters
!call mpas_log_write('--- enter subroutine SU_Compute_Diags:')
!call mpas_log_write('--- nw_profile = $i',intArgs=(/nw_profile/))
!call mpas_log_write('--- nw_vertint = $i',intArgs=(/nw_vertint/))
 if(associated(self%dmssmass)  ) self%dmssmass   = 0._RKIND
 if(associated(self%dmscmass)  ) self%dmscmass   = 0._RKIND
 if(associated(self%msasmass)  ) self%msasmass   = 0._RKIND
 if(associated(self%msacmass)  ) self%msacmass   = 0._RKIND
 if(associated(self%so2smass)  ) self%so2smass   = 0._RKIND
 if(associated(self%so2cmass)  ) self%so2cmass   = 0._RKIND
 if(associated(self%so4smass)  ) self%so4smass   = 0._RKIND
 if(associated(self%so4cmass)  ) self%so4cmass   = 0._RKIND
 if(associated(self%so4mass)   ) self%so4mass    = 0._RKIND
 if(associated(self%suconc )   ) self%suconc     = 0._RKIND
 if(associated(self%so4snum)   ) self%so4snum    = 0._RKIND
 if(associated(self%so4sarea)  ) self%so4sarea   = 0._RKIND
 if(associated(self%suexttau)  ) self%suexttau   = 0._RKIND
 if(associated(self%sustexttau)) self%sustexttau = 0._RKIND
 if(associated(self%suscatau)  ) self%suscatau   = 0._RKIND
 if(associated(self%sustscatau)) self%sustscatau = 0._RKIND
 if(associated(self%suextcoef) ) self%suextcoef  = 0._RKIND
 if(associated(self%suscacoef) ) self%suscacoef  = 0._RKIND
 if(associated(self%subckcoef) ) self%subckcoef  = 0._RKIND
 if(associated(self%suangstr)  ) self%suangstr   = 0._RKIND
 if(associated(self%sufluxu)   ) self%sufluxu    = 0._RKIND
 if(associated(self%sufluxv)   ) self%sufluxv    = 0._RKIND
 istat = 0
 call SU_Compute_Diags( &
              km         = self_params%km                                  , &
              klid       = self_params%klid                                , &
              rmed       = self_params%radius(nSO4)                        , &
              sigma      = self_params%sigma(nSO4)                         , &
              rhop       = self_params%rhop(nSO4)                          , &
              wavelengths_profile = self_params%wavelengths_profile*1.0e-9 , &
              wavelengths_vertint = self_params%wavelengths_vertint*1.0e-9 , &
              mie        = self_params%diag_Mie                            , &
              grav       = grav                                            , &
              pi         = pi                                              , &
              nSO4       = nSO4                                            , &
              tmpu       = self%t                                          , &
              rhoa       = self%airdens                                    , &
              delp       = self%delp                                       , &
              ple        = self%ple                                        , &
              tropp      = self%tropp                                      , &
              rh         = self%rh2                                        , &
              u          = self%u                                          , &
              v          = self%v                                          , &
              dms        = self%dms                                        , &
              so2        = self%so2                                        , &
              so4        = self%so4                                        , &
              msa        = self%msa                                        , &
              dmssfcmass = self%dmssmass                                   , &
              dmscolmass = self%dmscmass                                   , &
              msasfcmass = self%msasmass                                   , &
              msacolmass = self%msacmass                                   , &
              so2sfcmass = self%so2smass                                   , &
              so2colmass = self%so2cmass                                   , &
              so4sfcmass = self%so4smass                                   , &
              so4colmass = self%so4cmass                                   , &
              so4mass    = self%so4mass                                    , &
              so4conc    = self%suconc                                     , &
              snum       = self%so4snum                                    , &
              sarea      = self%so4sarea                                   , &
              exttau     = self%suexttau                                   , &
              scatau     = self%suscatau                                   , &
!             stexttau   = self%sustexttau                                 , &
!             stscatau   = self%sustscatau                                 , &
              extcoef    = self%suextcoef                                  , &
              scacoef    = self%suscacoef                                  , &
              bckcoef    = self%subckcoef                                  , &
              angstrom   = self%suangstr                                   , &
              fluxu      = self%sufluxu                                    , &
              fluxv      = self%sufluxv                                    , &
              rc         = istat                                             &
                      )
 if(istat /=0) then
    call mpas_log_write('--- SU2G_GridComp: error in subroutine SU_Compute_Diags.', &
                        messageType=MPAS_LOG_CRIT)
 else
!   call mpas_log_write('--- end subroutine SU_Compute_Diags.')
 endif


 i1 = lbound(self%rh2,1); i2 = ubound(self%rh2,1) 
 j1 = lbound(self%rh2,2); j2 = ubound(self%rh2,2) 
 km = ubound(self%rh2,3)                      

!call mpas_log_write('--- enter subroutine SU_Compute_Diags RH20:')
 if(associated(self%suextcoefrh20)) self%suextcoefrh20(:,:,:,:) = 0._RKIND
 if(associated(self%suscacoefrh20)) self%suscacoefrh20(:,:,:,:) = 0._RKIND
 if(.not.allocated(rh20)) allocate(rh20(i1:i2,j1:j2,km))
 rh20(:,:,:) = 0.20
 call SU_Compute_Diags( &
              km                  = self_params%km,                         &
              klid                = self_params%klid,                       &
              rmed                = self_params%radius(nSO4),               &
              sigma               = self_params%sigma(nSO4),                &
              rhop                = self_params%rhop(nSO4),                 &
              wavelengths_profile = self_params%wavelengths_profile*1.0e-9, &
              wavelengths_vertint = self_params%wavelengths_vertint*1.0e-9, &
              mie                 = self_params%diag_Mie,                   &
              grav                = grav,                                   &
              pi                  = pi,                                     &
              rh                  = rh20,                                   &
              nSO4                = nSO4,                                   &
              tmpu                = self%t,                                 &
              rhoa                = self%airdens,                           &
              delp                = self%delp,                              &
              ple                 = self%ple,                               &
              tropp               = self%tropp,                             &
              u                   = self%u,                                 &
              v                   = self%v,                                 &
              dms                 = self%dms,                               &
              so2                 = self%so2,                               &
              so4                 = self%so4,                               &
              msa                 = self%msa,                               &
              extcoef             = self%suextcoefrh20,                     &
              scacoef             = self%suscacoefrh20,                     &
              rc = istat                                                    &
                      )
 if(istat /=0) then
    call mpas_log_write('--- SU2G_GridComp: error in subroutine SU_Compute_Diags RH20.', &
                        messageType=MPAS_LOG_CRIT)
 else
!   call mpas_log_write('--- end subroutine SU_Compute_Diags RH20.')
 endif


!call mpas_log_write('--- enter subroutine SU_Compute_Diags RH80:')
 if(associated(self%suextcoefrh80)) self%suextcoefrh80(:,:,:,:) = 0._RKIND
 if(associated(self%suscacoefrh80)) self%suscacoefrh80(:,:,:,:) = 0._RKIND
 if(.not.allocated(rh80)) allocate(rh80(i1:i2,j1:j2,km))
 rh80(:,:,:) = 0.80
 call SU_Compute_Diags( &
              km                  = self_params%km,                         &
              klid                = self_params%klid,                       &
              rmed                = self_params%radius(nSO4),               &
              sigma               = self_params%sigma(nSO4),                &
              rhop                = self_params%rhop(nSO4),                 &
              wavelengths_profile = self_params%wavelengths_profile*1.0e-9, &
              wavelengths_vertint = self_params%wavelengths_vertint*1.0e-9, &
              mie                 = self_params%diag_Mie,                   &
              grav                = grav,                                   &
              pi                  = pi,                                     &
              rh                  = rh80,                                   &
              nSO4                = nSO4,                                   &
              tmpu                = self%t,                                 &
              rhoa                = self%airdens,                           &
              delp                = self%delp,                              &
              ple                 = self%ple,                               &
              tropp               = self%tropp,                             &
              u                   = self%u,                                 &
              v                   = self%v,                                 &
              dms                 = self%dms,                               &
              so2                 = self%so2,                               &
              so4                 = self%so4,                               &
              msa                 = self%msa,                               &
              extcoef             = self%suextcoefrh80,                     &
              scacoef             = self%suscacoefrh80,                     &
              rc = istat                                                    &
                      )
 if(istat /=0) then
    call mpas_log_write('--- SU2G_GridComp: error in subroutine SU_Compute_Diags RH80.', &
                        messageType=MPAS_LOG_CRIT)
 else
!   call mpas_log_write('--- end subroutine SU_Compute_Diags RH80.')
 endif
 if(allocated(rh20)) deallocate(rh20)
 if(allocated(rh80)) deallocate(rh80)


!--- deallocate SU2G oxidants:
 if(allocated(xoh)      ) deallocate(xoh      )
 if(allocated(xh2o2)    ) deallocate(xh2o2    )
 if(allocated(xno3)     ) deallocate(xno3     )
 if(allocated(h2o2_init)) deallocate(h2o2_init)

 if(allocated(dms_init)) deallocate(dms_init)
 if(allocated(so2_init)) deallocate(so2_init)
 if(allocated(so4_init)) deallocate(so4_init)
 if(allocated(msa_init)) deallocate(msa_init)

 call mpas_log_write('--- end subroutine processes_SU2G_GridComp.')

 end subroutine processes_SU2G_GridComp

!=================================================================================================================
 end module SU2G_GridCompMod
!=================================================================================================================

