! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================================
 module CA2G_bc_GridCompMod
 use mpas_kind_types,only: RKIND
 use mpas_derived_types,only: MPAS_LOG_CRIT
 use mpas_log

 use GA_EnvironmentMod
 use GOCART2G_AeroGeneric,only: findKlid
 use GOCART2G_MieMod_smiol
 use GOCART2G_Process,only: Aero_Compute_Diags,  &
                            Chem_BiomassDiurnal, &
                            Chem_Settling,       &
                            CAEmission,          &
                            DryDeposition,       &
                            phobicTophilic,      &
                            WetRemovalGOCART2G
 use MAPL,only: MAPL_PackTime

 use CA2G_bc_instance,only: nbins,particle_radius_microns,particle_density,fscav,molecular_weight, &
                            fnum,rhFlag,pressure_lid_in_hPa,sigma,hydrophobic_fraction,            &
                            aircraft_fuel_emission_factor,aviation_vertical_layers,                &
                            point_emissions_srcfilen
 use CA2G_bc_StateSpecs,only: CA2G_bc_State


 implicit none
 private


!--- constants (these parameters need to be accessed from MPAS physics instead of redefined here):
 real(kind=RKIND),parameter:: cpd      = 1003.0_RKIND
 real(kind=RKIND),parameter:: grav     = 9.80616_RKIND
 real(kind=RKIND),parameter:: pi       = 3.141592653589793_RKIND
 real(kind=RKIND),parameter:: karman   = 0.4_RKIND
 real(kind=RKIND),parameter:: radTodeg = 180._RKIND/pi
 real(kind=RKIND),parameter:: undefval = 1.0e15


!--- types needed to define CA2G_bc:
 type :: ThreadWorkspace
    integer:: nPts = -1
    integer,dimension(:),allocatable:: pstart,pend

    real(kind=RKIND),dimension(:),allocatable:: pLat,pLon,pBase,pTop,pEmis
 end type ThreadWorkspace

 type,extends(GA_Environment),public:: CA2G_bc_GridComp
    logical:: diurnal_bb                              ! diurnal biomass burning
    integer:: myDOW = -1                              ! day of the week: Sun=1, Mon=2,...,Sat=7

    real(kind=RKIND):: eAircraftFuel                  ! aircraft emission factor: go from kg fuel to kg SO2
    real(kind=RKIND):: aviation_layers(4)             ! heights of the LTO, CDS and CRS layers
    real(kind=RKIND):: fMonoterpenes = 0.0            ! fraction of monoterpene emissions -> aerosol
    real(kind=RKIND):: fIsoprene = 0.0                ! fraction of isoprene emissions -> aerosol
    real(kind=RKIND):: fHydrophobic                   ! initially hydrophobic portion
    real(kind=RKIND):: ratPOM = 1.0                   ! ratio of POM to OC mass
    real(kind=RKIND),dimension(:),allocatable:: sigma ! sigma of lognormal number distribution

    !workspace for point emissions:
    logical:: doing_point_emissions = .false.
    character(len=255):: point_emissions_srcfilen  ! filename for pointwise emissions

    type(ThreadWorkspace),dimension(:),allocatable :: workspaces

    contains
       procedure:: emissions_GridComp => emissions_CA2G_bc_GridComp
       procedure:: load_GridComp      => load_CA2G_bc_GridComp
       procedure:: processes_GridComp => processes_CA2G_bc_GridComp
 end type CA2G_bc_GridComp

 type wrap_
    type(CA2G_bc_GridComp),pointer:: PTR !=> null()
 end type wrap_


 contains


!==================================================================================================================
 subroutine load_CA2G_bc_GridComp(self,kts,kte)
!==================================================================================================================

!--- input arguments:
 integer,intent(in):: kts,kte

!--- inout arguments:
 class(CA2G_bc_GridComp),intent(inout) :: self

!local variables:
 integer:: n

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine load_CA2G_bc_GridComp:')


!--- initialization of vertical index:
 self%klid = kts
 self%km = kte-kts+1


!--- initialization using parameters defined in GA_Environment:
 self%diurnal_bb = .false.

 call self%load_from_config(nbins,particle_radius_microns,particle_density,fscav,molecular_weight,fnum, &
                            rhFlag,pressure_lid_in_hPa)

 call mpas_log_write('--- nbins = $i',intArgs=(/self%nbins/))
 call mpas_log_write('--- radius,rhop,fscav,molwght,fnum,sigma:')
 do n = 1,self%nbins
    call mpas_log_write('$i $r $r $r $r $r',intArgs=(/n/),realArgs=(/self%radius(n),self%rhop(n), &
                        self%fscav(n),self%molwght(n),self%fnum(n)/))
 enddo


!--- initialization of all other variables in CA2G_bc_GridComp:
 self%fHydrophobic = hydrophobic_fraction
 self%point_emissions_srcfilen = trim(point_emissions_srcfilen)

 if(.not.allocated(self%sigma)) allocate(self%sigma(self%nbins))
 do n = 1,self%nbins
    self%sigma(n) = sigma(n)
 enddo

 self%eAirCraftFuel = aircraft_fuel_emission_factor
 do n = 1,4
    self%aviation_layers(n) = aviation_vertical_layers(n)
 enddo


 call mpas_log_write('--- end subroutine load_CA2G_bc_GridCOMP:')

 end subroutine load_CA2G_bc_GridComp

!==================================================================================================================
 subroutine emissions_CA2G_bc_GridComp(self_params,self,its,ite,jts,jte,kts,kte,iyr,imm,idd,ihr,imn,isc)
!==================================================================================================================

!--- input arguments:
 integer,intent(in):: its,ite,jts,jte,kts,kte
 integer,intent(in):: iyr,imm,idd,ihr,imn,isc

 class(CA2G_bc_GridComp),intent(in):: self_params

!--- inout arguments:
 class(CA2G_bc_State),intent(inout):: self

!--- local variables and arrays:
 integer:: nymd,nhms
 integer:: i,j,istat

 real(kind=RKIND),dimension(:,:),allocatable:: biomass_src         !
 real(kind=RKIND),dimension(:,:),allocatable:: biofuel_src         !
 real(kind=RKIND),dimension(:,:),allocatable:: eocant1_src         !
 real(kind=RKIND),dimension(:,:),allocatable:: eocant2_src         !
 real(kind=RKIND),dimension(:,:),allocatable:: oc_ship_src         !
 real(kind=RKIND),dimension(:,:),allocatable:: aviation_lto_src    !
 real(kind=RKIND),dimension(:,:),allocatable:: aviation_cds_src    !
 real(kind=RKIND),dimension(:,:),allocatable:: aviation_crs_src    !
 real(kind=RKIND),dimension(:,:,:),allocatable:: aircraft_fuel_src !

 real(kind=RKIND),dimension(:,:),allocatable:: biomass_src_
 real(kind=RKIND),dimension(:,:),allocatable:: biogvoc_src

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine emissions_CA2G_bc_GridComp:')


!--- extract nymd and nhms from clock:
 call MAPL_PackTime(nymd,iyr,imm,idd)
 call MAPL_PackTime(nhms,ihr,imn,isc)
 call mpas_log_write('--- nymd = $i',intArgs=(/nymd/))
 call mpas_log_write('--- nhms = $i',intArgs=(/nhms/))


!--- initialize black carbon emissions: at this time, only anthropogenic surface emission is available as an
!    input stream and put into eocant1_src. all other emission types are set to zero.
 biomass_src       = self%bc_biomass
 biofuel_src       = self%bc_biofuel
 eocant1_src       = self%bc_antebc1
 eocant2_src       = self%bc_antebc2
 oc_ship_src       = self%bc_ship
 aviation_lto_src  = self%bc_aviation_lto
 aviation_cds_src  = self%bc_aviation_cds
 aviation_crs_src  = self%bc_aviation_crs
 aircraft_fuel_src = self%bc_aircraft

 biomass_src(:,:)         = 0._RKIND
 biofuel_src(:,:)         = 0._RKIND
!eocant1_src(:,:)         = 0._RKIND
 eocant2_src(:,:)         = 0._RKIND
 oc_ship_src(:,:)         = 0._RKIND
 aviation_lto_src(:,:)    = 0._RKIND
 aviation_cds_src(:,:)    = 0._RKIND
 aviation_crs_src(:,:)    = 0._RKIND
 aircraft_fuel_src(:,:,:) = 0._RKIND

 if(.not.allocated(biogvoc_src)) allocate(biogvoc_src(its:ite,jts:jte))
 biogvoc_src(:,:) = 0._RKIND


!--- as a safety check, all undefined values are set to zero. this may be needed when all emission types become
!    available.
!where(1.01*biomass_src > undefval) biomass_src = 0._RKIND
!where(1.01*biogvoc_src > undefval) biogvoc_src = 0._RKIND
!where(1.01*biofuel_src > undefval) biofuel_src = 0._RKIND
 where(1.01*eocant1_src > undefval) eocant1_src = 0._RKIND
!where(1.01*eocant2_src > undefval) eocant2_src = 0._RKIND
!where(1.01*oc_ship_src > undefval) oc_ship_src = 0._RKIND
!where(1.01*aviation_lto_src  > undefval) aviation_lto_src  = 0._RKIND
!where(1.01*aviation_cds_src  > undefval) aviation_cds_src  = 0._RKIND
!where(1.01*aviation_crs_src  > undefval) aviation_crs_src  = 0._RKIND
!where(1.01*aircraft_fuel_src > undefval) aircraft_fuel_src = 0._RKIND


!--- apply diurnal cycle to biomass burning if needed:
 if(self_params%diurnal_bb ) then
    call mpas_log_write('--- enter subroutine Chem_BiomassDiurnal:')
    biomass_src_      = self%bc_biomass
    call Chem_BiomassDiurnal( &
       cdt  = self_params%cdt,    &
       nhms = nhms,               &
       eout = biomass_src,        &
       ein  = biomass_src_,       &
       lons = self%lons*radTodeg, &
       lats = self%lats*radTodeg  &
                            )
    call mpas_log_write('--- end subroutine Chem_BiomassDiurnal:')
 endif


!--- apply emissions to CA2G_bc:
 istat = 0
 call mpas_log_write('--- enter subroutine CAEmission:')
 call CAEmission( &
    mie               = self_params%diag_Mie,        &
    km                = self_params%km,              &
    cdt               = self_params%cdt,             &
    ratPOM            = self_params%ratPOM,          &
    fHydrophobic      = self_params%fHydrophobic,    &
    eAircraftfuel     = self_params%eAircraftfuel,   &
    aviation_layers   = self_params%aviation_layers, &
    nbins             = self_params%nbins,           &
    grav              = grav,                        &
    prefix            = 'BC',                        &
    biomass_src       = biomass_src,                 &
    biofuel_src       = biofuel_src,                 &
    terpene_src       = biogvoc_src,                 &
    eocant1_src       = eocant1_src,                 &
    eocant2_src       = eocant2_src,                 &
    oc_ship_src       = oc_ship_src,                 &
    aircraft_fuel_src = aircraft_fuel_src,           &
    aviation_lto_src  = aviation_lto_src,            &
    aviation_cds_src  = aviation_cds_src,            &
    aviation_crs_src  = aviation_crs_src,            &
    pblh              = self%zpbl,                   &
    tmpu              = self%t,                      &
    rhoa              = self%airdens,                &
    rh                = self%rh2,                    &
    delp              = self%delp,                   &
    aerosolPhilic     = self%bcPHOBIC,               &
    aerosolPhobic     = self%bcPHILIC,               &
    oc_emis           = self%bcEM,                   &
    oc_emisan         = self%bcEMAN,                 &
    oc_emisbb         = self%bcEMBB,                 &
    oc_emisbf         = self%bcEMBF,                 &
    oc_emisbg         = self%bcEMBG,                 &
    rc                = istat                        &
                )
 if(istat /=0) then
    call mpas_log_write('--- CA2G_bc_GridComp: error in subroutine CAEmission:', &
                        messageType=MPAS_LOG_CRIT)
 else
    call mpas_log_write('--- end subroutine CAEmission:')
 endif


!--- for now, we do not support point emissions:


 call mpas_log_write('--- end subroutine emissions_CA2G_bc_GridComp:')

 end subroutine emissions_CA2G_bc_GridComp

!==================================================================================================================
 subroutine processes_CA2G_bc_GridComp(self_params,self,its,ite,jts,jte,kts,kte)
!==================================================================================================================
!--- input arguments:
 integer,intent(in):: its,ite,jts,jte,kts,kte

!--- inout arguments:
 class(CA2G_bc_GridComp),intent(inout):: self_params
 class(CA2G_bc_State),intent(inout):: self

!--- local variables and arrays:
 logical:: KIN

 integer:: i,i1,i2,j,j1,j2,k,km,ibin,n
 integer:: istat
 integer:: n_profile,n_vertint

 real(kind=RKIND):: fwet
 real(kind=RKIND),dimension(:,:,:),pointer:: rh20,rh80
 real(kind=RKIND),dimension(:,:),allocatable:: drydepf,dqa
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: qca2G

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine processes_CA2G_bc_GridComp:')


!--- add SOA from anthropogenic VOC oxidation: no oxidation for black carbon:


!--- add hoc transfer of hydrophobic to hydrophilic aerosols following Chin's parameterization:
!    the rate constant is k = 4.63e-6 s-1 (.4 day-1; e-folding time = 2.5 days)
 call mpas_log_write('--- enter subroutine phobicTophilic:')
 if(associated(self%bchyphil)) self%bchyphil(:,:) = 0._RKIND
 istat = 0
 call phobicTophilic( &
           aerosol_phobic        = self%bcPHOBIC   , &
           aerosol_philic        = self%bcPHILIC   , &
           aerosol_toHydrophilic = self%bcHYPHIL   , &
           km                    = self_params%km  , &
           cdt                   = self_params%cdt , &
           grav                  = grav            , &
           delp                  = self%delp       , &
           rc = istat                                &
                    )
 if(istat /=0) then
    call mpas_log_write('--- CA2G_bc_GridComp: error in subroutine phobicTophilic', &
                        messageType=MPAS_LOG_CRIT)
 else
    call mpas_log_write('--- end subroutine phobicTophilic:')
 endif


!--- CA2G_bc settling:
 call mpas_log_write('--- enter subroutine Chem_Settling:')
 if(.not.allocated(qca2G)) allocate(qca2G(its:ite,jts:jte,kts:kte,self_params%nbins))
 do j = jts,jte
    do i = its,ite
       do k = kts,kte
          qca2G(i,j,k,1) = self%bcphobic(i,j,k)
          qca2G(i,j,k,2) = self%bcphilic(i,j,k)
       enddo
    enddo
 enddo

 do ibin = 1, self_params%nbins
    !if radius == 0, then we're dealing with a gas which has no settling losses:
    if(self_params%radius(ibin) == 0.0) then
       if(associated(self%bcSD)) self%bcSD(:,:,ibin) = 0._RKIND
       cycle
    endif
    istat = 0
    call Chem_Settling( &
              km        = self_params%km                 , &
              klid      = self_params%klid               , &
              bin       = ibin                           , &
              flag      = self_params%rhFlag             , &
              cdt       = self_params%cdt                , &
              grav      = grav                           , &
              radiusInp = self_params%radius(ibin)*1.e-6 , &
              rhopInp   = self_params%rhop(ibin)         , &
              int_qa    = qca2G(:,:,:,ibin)              , &
              tmpu      = self%t                         , &
              rhoa      = self%airdens                   , &
              rh        = self%rh2                       , &
              hghte     = self%zle                       , &
              delp      = self%delp                      , &
              fluxout   = self%bcSD                      , &
              rc        = istat                            &
                      )
 enddo
 if(istat /=0) then
    call mpas_log_write('--- CA2G_bc_GridComp: error in subroutine Chem_Settling', &
                        messageType=MPAS_LOG_CRIT)
 else
    call mpas_log_write('--- end subroutine Chem_Settling:')
 endif


!--- CA2G_bc dry deposition:
 call mpas_log_write('--- enter subroutine DryDeposition:')
 if(associated(self%bcdp)  ) self%bcdp(:,:,:) = 0._RKIND
 if(.not.allocated(dqa)    ) allocate(dqa(its:ite,jts:jte)    )
 if(.not.allocated(drydepf)) allocate(drydepf(its:ite,jts:jte))
 drydepf = 0.
 istat   = 0
 call DryDeposition( &
              km         = self_params%km , &
              tmpu       = self%t         , &
              rhoa       = self%airdens   , &
              hghte      = self%zle       , &
              oro        = self%lwi       , &
              ustar      = self%ustar     , &
              pblh       = self%zpbl      , &
              shflux     = self%sh        , &
              von_karman = karman         , &
              cpd        = cpd            , &
              grav       = grav           , &
              z0h        = self%z0h       , &
              drydepf    = drydepf        , &
              rc         = istat            &
                   )
 do ibin = 1, self_params%nbins
    dqa = 0.
    dqa = max(0.0,qca2G(:,:,self_params%km,ibin)*(1.-exp(-drydepf*self_params%cdt)))
    qca2G(:,:,self_params%km,ibin) = qca2G(:,:,self_params%km,ibin) - dqa
    if(associated(self%bcDP)) then
       self%bcdp(:,:,ibin) = dqa*self%delp(:,:,self_params%km)/grav/self_params%cdt
    end if
 enddo
 if(istat /=0) then
    call mpas_log_write('--- CA2G_bc_GridComp: error in subroutine DryDeposition', &
                        messageType=MPAS_LOG_CRIT)
 else
    if(allocated(dqa)    ) deallocate(dqa    )
    if(allocated(drydepf)) deallocate(drydepf)
    call mpas_log_write('--- end subroutine DryDeposition:')
 endif


 do j = jts,jte
    do i = its,ite
       do k = kts,kte
          self%bcPHOBIC(i,j,k) = qca2G(i,j,k,1)
          self%bcPHILIC(i,j,k) = qca2G(i,j,k,2)
       enddo
    enddo
 enddo


!--- CA2G_bc large-scale wet removal (hydrophilic mode is removed):
 call mpas_log_write('--- enter subroutine WetRemovalGOCART2G:')
 if(associated(self%bcWT)) self%bcWT(:,:,:) = 0._RKIND
 KIN   = .true.
 fwet  = 1._RKIND
 istat = 0
 call WetRemovalGOCART2G( &
              km        = self_params%km    , &
              klid      = self_params%klid  , &
              n1        = self_params%nbins , &
              n2        = self_params%nbins , &
              bin_ind   = 2                 , &
              cdt       = self_params%cdt   , &
              aero_type = 'BC'              , &
              kin       = KIN               , &
              grav      = grav              , &
              fwet      = fwet              , &
              aerosol   = self%bcPHILIC     , &
              ple       = self%ple          , &
              tmpu      = self%t            , &
              rhoa      = self%airdens      , &
              pfllsan   = self%pfl_lsan     , &
              pfilsan   = self%pfi_lsan     , &
              precc     = self%cn_prcp      , &
              precl     = self%ncn_prcp     , &
              fluxout   = self%bcWT         , &
              rc        = istat               &
                        )
 if(istat /=0) then
    call mpas_log_write('--- CA2G_bc_GridComp: error in subroutine WetRemovalGOCART2G', &
                        messageType=MPAS_LOG_CRIT)
 else
    call mpas_log_write('--- end subroutine WetRemovalGOCART2G:')
 endif


!--- CA2G_bc diagnostics:
 do j = jts,jte
    do i = its,ite
       do k = kts,kte
          qca2G(i,j,k,1) = self%bcphobic(i,j,k)
          qca2G(i,j,k,2) = self%bcphilic(i,j,k)
       enddo
    enddo
 enddo
 n_profile = size(self_params%wavelengths_profile)
 n_vertint = size(self_params%wavelengths_vertint)
 call mpas_log_write('--- enter subroutine Aero_Compute_Diags:')
 call mpas_log_write('--- nbins     = $i',intArgs=(/nbins/))
 call mpas_log_write('--- n_profile = $i',intArgs=(/n_profile/))
 call mpas_log_write('--- n_vertint = $i',intArgs=(/n_vertint/))
 if(associated(self%bcsmass)   ) self%bcSMASS(:,:)       = 0._RKIND
 if(associated(self%bccmass)   ) self%bcCMASS(:,:)       = 0._RKIND
 if(associated(self%bcmass )   ) self%bcMASS(:,:,:)      = 0._RKIND
 if(associated(self%bcexttau)  ) self%bcEXTTAU(:,:,:)    = 0._RKIND
 if(associated(self%bcstexttau)) self%bcSTEXTTAU(:,:,:)  = 0._RKIND
 if(associated(self%bcscatau)  ) self%bcSCATAU(:,:,:)    = 0._RKIND
 if(associated(self%bcstscatau)) self%bcSTSCATAU(:,:,:)  = 0._RKIND
 if(associated(self%bcfluxu)   ) self%bcFLUXU(:,:)       = 0._RKIND
 if(associated(self%bcfluxv)   ) self%bcFLUXV(:,:)       = 0._RKIND
 if(associated(self%bcconc)    ) self%bcCONC(:,:,:)      = 0._RKIND
 if(associated(self%bcextcoef) ) self%bcEXTCOEF(:,:,:,:) = 0._RKIND
 if(associated(self%bcscacoef) ) self%bcSCACOEF(:,:,:,:) = 0._RKIND
 if(associated(self%bcbckcoef) ) self%bcBCKCOEF(:,:,:,:) = 0._RKIND
 if(associated(self%bcangstr)  ) self%bcANGSTR(:,:)      = 0._RKIND
 if(associated(self%bcaeridx)  ) self%bcAERIDX(:,:)      = 0._RKIND
 istat = 0
 call Aero_Compute_Diags( &
              mie                 = self_params%diag_Mie                   , &
              km                  = self_params%km                         , &
              klid                = self_params%klid                       , &
              nbegin              = 1                                      , &
              nbins               = 2                                      , &
              wavelengths_profile = self_params%wavelengths_profile*1.0e-9 , &
              wavelengths_vertint = self_params%wavelengths_vertint*1.0e-9 , &
              aerosol             = qca2G                                  , &
              grav                = grav                                   , &
              tmpu                = self%t                                 , &
              rhoa                = self%airdens                           , &
              rh                  = self%rh2                               , &
              u                   = self%u                                 , &
              v                   = self%v                                 , &
              delp                = self%delp                              , &
              ple                 = self%ple                               , &
              tropp               = self%tropp                             , &
              sfcmass             = self%bcSMASS                           , &
              colmass             = self%bcCMASS                           , &
              mass                = self%bcMASS                            , &
              exttau              = self%bcEXTTAU                          , &
              scatau              = self%bcSCATAU                          , &
!             stexttau            = self%bcSTEXTTAU                        , &
!             stscatau            = self%bcSTSCATAU                        , &
              fluxu               = self%bcFLUXU                           , &
              fluxv               = self%bcFLUXV                           , &
              conc                = self%bcCONC                            , &
              extcoef             = self%bcEXTCOEF                         , &
              scacoef             = self%bcSCACOEF                         , &
              bckcoef             = self%bcBCKCOEF                         , &
              angstrom            = self%bcANGSTR                          , &
              aerindx             = self%bcAERIDX                          , &
              NO3nFlag            = .false.                                , &
              rc                  = istat                                    &
                        )
 if(istat /=0) then
    call mpas_log_write('--- CA2G_bc_GridComp: error in subroutine Aero_Compute_Diags', &
                        messageType=MPAS_LOG_CRIT)
 else
    call mpas_log_write('--- end subroutine Aero_Compute_Diags:')
 endif


 i1 = lbound(self%rh2,1); i2 = ubound(self%rh2,1)
 j1 = lbound(self%rh2,2); j2 = ubound(self%rh2,2)
 km = ubound(self%rh2,3)

 call mpas_log_write('--- enter subroutine Aero_Compute_Diags RH20:')
 if(.not.associated(rh20)) allocate(rh20(i1:i2,j1:j2,km))
 if(associated(self%bcEXTCOEFRH20)) self%bcEXTCOEFRH20(:,:,:,:) = 0._RKIND
 if(associated(self%bcSCACOEFRH20)) self%bcEXTCOEFRH20(:,:,:,:) = 0._RKIND
 rh20(:,:,:) = 0.20
 istat = 0
 call Aero_Compute_Diags( &
           mie                 = self_params%diag_Mie                   , &
           km                  = self_params%km                         , &
           klid                = self_params%klid                       , &
           nbegin              = 1                                      , &
           nbins               = self_params%nbins                      , &
           wavelengths_profile = self_params%wavelengths_profile*1.0e-9 , &
           wavelengths_vertint = self_params%wavelengths_vertint*1.0e-9 , &
           aerosol             = qca2G                                  , &
           grav                = grav                                   , &
           tmpu                = self%t                                 , &
           rhoa                = self%airdens                           , &
           rh                  = rh20                                   , &
           u                   = self%u                                 , &
           v                   = self%v                                 , &
           delp                = self%delp                              , &
           ple                 = self%ple                               , &
           tropp               = self%tropp                             , &
           extcoef             = self%bcEXTCOEFRH20                     , &
           scacoef             = self%bcSCACOEFRH20                     , &
           NO3nFlag            = .false.                                , &
           rc                  = istat                                    &
                        )
 if(istat /=0) then
    call mpas_log_write('--- CA2G_bc_GridComp: error in subroutine Aero_Compute_Diags RH20', &
                        messageType=MPAS_LOG_CRIT)
 else
    call mpas_log_write('--- end subroutine Aero_Compute_Diags RH20:')
 endif


 call mpas_log_write('--- enter subroutine Aero_Compute_Diags RH80:')
 if(.not.associated(rh80)) allocate(rh80(i1:i2,j1:j2,km))
 if(associated(self%bcEXTCOEFRH80)) self%bcEXTCOEFRH80(:,:,:,:) = 0._RKIND
 if(associated(self%bcSCACOEFRH80)) self%bcSCACOEFRH80(:,:,:,:) = 0._RKIND
 rh80(:,:,:) = 0.80
 istat = 0
 call Aero_Compute_Diags( &
           mie                 = self_params%diag_Mie                   , &
           km                  = self_params%km                         , &
           klid                = self_params%klid                       , &
           nbegin              = 1                                      , &
           nbins               = 2                                      , &
           wavelengths_profile = self_params%wavelengths_profile*1.0e-9 , &
           wavelengths_vertint = self_params%wavelengths_vertint*1.0e-9 , &
           aerosol             = qca2G                                  , &
           grav                = grav                                   , &
           tmpu                = self%t                                 , &
           rhoa                = self%airdens                           , &
           rh                  = rh80                                   , &
           u                   = self%u                                 , &
           v                   = self%v                                 , &
           delp                = self%delp                              , &
           ple                 = self%ple                               , &
           tropp               = self%tropp                             , &
           extcoef             = self%bcEXTCOEFRH80                     , &
           scacoef             = self%bcSCACOEFRH80                     , &
           NO3nFlag            = .false.                                , &
           rc                  = istat                                    &
                        )
 if(istat /=0) then
    call mpas_log_write('--- CA2G_bc_GridComp: error in subroutine Aero_Compute_Diags RH80', &
                        messageType=MPAS_LOG_CRIT)
 else
    call mpas_log_write('--- end subroutine Aero_Compute_Diags RH80:')
 endif
 if(associated(rh20)) deallocate(rh20)
 if(associated(rh80)) deallocate(rh80)

 if(allocated(qca2G)) deallocate(qca2G)


 call mpas_log_write('--- end subroutine processes_CA2G_bc_GridComp:')

 end subroutine processes_CA2G_bc_GridComp

!==================================================================================================================
 end module CA2G_bc_GridCompMod
!==================================================================================================================
