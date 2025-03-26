! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module CA2G_oc_GridCompMod
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
                            WetRemovalGocart2G
 use MAPL,only: MAPL_PackTime

 use CA2G_oc_instance,only: nbins,particle_radius_microns,particle_density,fscav,molecular_weight,   &
                            fnum,rhFlag,pressure_lid_in_hPa,sigma,hydrophobic_fraction,pom_ca_ratio, &
                            monoterpenes_emission_fraction,isoprene_emission_fraction,               &
                            aircraft_fuel_emission_factor,aviation_vertical_layers,                  &
                            point_emissions_srcfilen
 use CA2G_oc_StateSpecs,only: CA2G_oc_State


 implicit none
 private


!--- constants (these parameters need to be accessed from MPAS physics instead of redefined here):
 real(kind=RKIND),parameter:: cpd      = 1003.0_RKIND
 real(kind=RKIND),parameter:: grav     = 9.80616_RKIND
 real(kind=RKIND),parameter:: pi       = 3.141592653589793_RKIND
 real(kind=RKIND),parameter:: karman   = 0.4_RKIND
 real(kind=RKIND),parameter:: radTodeg = 180._RKIND/pi
 real(kind=RKIND),parameter:: undefval = 1.0e15


!--- types needed to define CA2G_oc:
 type :: ThreadWorkspace
    integer:: nPts = -1
    integer,dimension(:),allocatable:: pstart,pend

    real(kind=RKIND),dimension(:),allocatable:: pLat,pLon,pBase,pTop,pEmis
 end type ThreadWorkspace

 type,extends(GA_Environment),public:: CA2G_oc_GridComp
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
       procedure:: emissions_GridComp => emissions_CA2G_oc_GridComp
       procedure:: load_GridComp      => load_CA2G_oc_GridComp
       procedure:: processes_GridComp => processes_CA2G_oc_GridComp
 end type CA2G_oc_GridComp

 type wrap_
    type(CA2G_oc_GridComp),pointer:: PTR !=> null()
 end type wrap_


 contains


!=================================================================================================================
 subroutine load_CA2G_oc_GridComp(self,kts,kte)
!=================================================================================================================

!--- input arguments:
 integer,intent(in):: kts,kte

!--- inout arguments:
 class(CA2G_oc_GridComp),intent(inout) :: self

!local variables:
 integer:: n

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine load_CA2G_oc_GridComp:')

!--- initialization of vertical index:
 self%klid = kts
 self%km = kte-kts+1


!--- initialization using parameters defined in GA_Environment:
 self%diurnal_bb = .false.

 call self%load_from_config(nbins,particle_radius_microns,particle_density,fscav,molecular_weight,fnum, &
                            rhFlag,pressure_lid_in_hPa)

 call mpas_log_write('--- nbins = $i',intArgs=(/self%nbins/))
 call mpas_log_write('--- radius,rhop,fscav,molwght,fnum:')
 do n = 1,self%nbins
    call mpas_log_write('$i $r $r $r $r $r',intArgs=(/n/),realArgs=(/self%radius(n),self%rhop(n), &
                        self%fscav(n),self%molwght(n),self%fnum(n)/))
 enddo


!--- initialization of all other variables in CA2G_oc_GridComp:
 self%fMonoterpenes = monoterpenes_emission_fraction
 self%fIsoprene     = isoprene_emission_fraction
 self%fHydrophobic  = hydrophobic_fraction
 self%ratPOM        = pom_ca_ratio
 self%point_emissions_srcfilen = trim(point_emissions_srcfilen)

 if(.not.allocated(self%sigma)) allocate(self%sigma(self%nbins))
 do n = 1,self%nbins
    self%sigma(n) = sigma(n)
 enddo

 self%eAirCraftFuel = aircraft_fuel_emission_factor
 do n = 1,4
    self%aviation_layers(n) = aviation_vertical_layers(n)
 enddo


 call mpas_log_write('--- end subroutine load_CA2G_oc_GridCOMP:')

 end subroutine load_CA2G_oc_GridComp

!==================================================================================================================
 subroutine emissions_CA2G_oc_GridComp(self_params,self,its,ite,jts,jte,kts,kte,iyr,imm,idd,ihr,imn,isc)
!==================================================================================================================

!--- input arguments:
 integer,intent(in):: its,ite,jts,jte,kts,kte
 integer,intent(in):: iyr,imm,idd,ihr,imn,isc

 class(CA2G_oc_GridComp),intent(in):: self_params

!--- inout arguments:
 class(CA2G_oc_State),intent(inout):: self

!--- local variables and arrays:
 integer:: nymd,nhms
 integer:: i,j,istat

 real(kind=RKIND),dimension(:,:),allocatable:: biomass_src_
 real(kind=RKIND),dimension(:,:),allocatable:: biogvoc_src

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine emissions_CA2G_oc_GridComp:')


!--- extract nymd and nhms from clock:
 call MAPL_PackTime(nymd,iyr,imm,idd)
 call MAPL_PackTime(nhms,ihr,imn,isc)
 call mpas_log_write('--- nymd = $i',intArgs=(/nymd/))
 call mpas_log_write('--- nhms = $i',intArgs=(/nhms/))


 if(.not.allocated(biogvoc_src)) allocate(biogvoc_src(its:ite,jts:jte))
 biogvoc_src(:,:) = 0._RKIND


!--- apply diurnal cycle to biomass burning if needed:
 if(self_params%diurnal_bb) then
    call mpas_log_write('--- enter subroutine Chem_Biomass Diurnal:')
    biomass_src_ = self%oc_biomass
    call Chem_BiomassDiurnal( &
       cdt  = self_params%cdt,    &
       nhms = nhms,               &
       eout = self%oc_biomass,    &
       ein  = biomass_src_,       &
       lons = self%lons*radTodeg, &
       lats = self%lats*radTodeg  &
                            )
    call mpas_log_write('--- end subroutine Chem_Biomass Diurnal.')
 endif


!--- apply emissions to CA2G_oc: 
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
    prefix            = 'OC',                        &
    terpene_src       = biogvoc_src,                 &
    biomass_src       = self%oc_biomass,             &
    biofuel_src       = self%oc_biofuel,             &
    eocant1_src       = self%oc_anteoc1,             &
    eocant2_src       = self%oc_anteoc2,             &
    oc_ship_src       = self%oc_ship,                &
    aircraft_fuel_src = self%oc_aircraft,            &
    aviation_lto_src  = self%oc_aviation_lto,        &
    aviation_cds_src  = self%oc_aviation_cds,        &
    aviation_crs_src  = self%oc_aviation_crs,        &
    pblh              = self%zpbl,                   &
    tmpu              = self%t,                      &
    rhoa              = self%airdens,                &
    rh                = self%rh2,                    &
    delp              = self%delp,                   &
    aerosolPhilic     = self%ocphilic,               &
    aerosolPhobic     = self%ocphobic,               &
    oc_emis           = self%ocem,                   &
    oc_emisan         = self%oceman,                 &
    oc_emisbb         = self%ocembb,                 &
    oc_emisbf         = self%ocembf,                 &
    oc_emisbg         = self%ocembg,                 &
    rc                = istat                        &
                )
 if(istat /=0) then
    call mpas_log_write('--- CA2G_oc_GridComp: error in subroutine CAEmission.', &
                        messageType=MPAS_LOG_CRIT)
 else
    call mpas_log_write('--- end subroutine CAEmission.')
 endif


!--- for now, we do not support point emissions:


 call mpas_log_write('--- end subroutine emissions_CA2G_oc_GridComp.')

 end subroutine emissions_CA2G_oc_GridComp

!=================================================================================================================
 subroutine processes_CA2G_oc_GridComp(self_params,self,its,ite,jts,jte,kts,kte)
!=================================================================================================================
!--- input arguments:
 integer,intent(in):: its,ite,jts,jte,kts,kte

!--- inout arguments:
 class(CA2G_oc_GridComp),intent(inout):: self_params
 class(CA2G_oc_State),intent(inout):: self

!--- local variables:
 logical:: KIN

 integer:: i,i1,i2,j,j1,j2,k,km,ibin,n
 integer:: istat
 integer:: n_profile,n_vertint

 real(kind=RKIND):: fwet
!real(kind=RKIND),dimension(:,:,:),pointer:: rh20,rh80
 real(kind=RKIND),dimension(:,:),allocatable:: drydepf,dqa
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: qca2G

 real(kind=RKIND),allocatable,dimension(:,:,:),target:: rh20,rh80

!-----------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine processes_CA2G_oc_GridComp:')


 if(associated(self%ocpsoa)) self%ocpsoa(:,:) = 0._RKIND

 if(associated(self%psoa_anthro_voc)) then
    self%ocphilic(:,:,:) = self%ocphilic(:,:,:) + self%psoa_anthro_voc(:,:,:)*self_params%cdt/self%airdens(:,:,:)
    if(associated(self%ocpsoa)) then
       do k = kts,kte
          self%ocpsoa(:,:) = self%ocpsoa(:,:) + self%psoa_anthro_voc(:,:,k)*self%delp(:,:,k)/grav &
                           / self%airdens(:,:,k)
       enddo
    endif
 endif

 if(associated(self%psoap_biog_voc)) then
    self%ocphilic(:,:,:) = self%ocphilic(:,:,:) + self%psoap_biog_voc(:,:,:)*self_params%cdt/self%airdens(:,:,:)
    if(associated(self%ocpsoa)) then
       do k = kts,kte
          self%ocpsoa(:,:) = self%ocpsoa(:,:) + self%psoap_biog_voc(:,:,k)*self%delp(:,:,k)/grav &
                           / self%airdens(:,:,k)
       enddo
    endif
 endif

 if(associated(self%psoas_biog_voc)) then
    self%ocphilic(:,:,:) = self%ocphilic(:,:,:) + self%psoas_biog_voc(:,:,:)*self_params%cdt/self%airdens(:,:,:)
    if(associated(self%ocpsoa)) then
       do k = kts,kte
          self%ocpsoa(:,:) = self%ocpsoa(:,:) + self%psoa_anthro_voc(:,:,k)*self%delp(:,:,k)/grav &
                           / self%airdens(:,:,k)
       enddo
    endif
 endif


!--- add hoc transfer of hydrophobic to hydrophilic aerosols following Chin's parameterization:
!    the rate constant is k = 4.63e-6 s-1 (.4 day-1; e-folding time = 2.5 days)
 call mpas_log_write('--- enter subroutine phobicTophilic:')
 if(associated(self%ochyphil)) self%ochyphil(:,:) = 0._RKIND
 istat = 0
 call phobicTophilic( &
           aerosol_phobic        = self%ocphobic   , &
           aerosol_philic        = self%ocphilic   , &
           aerosol_toHydrophilic = self%ochyphil   , &
           km                    = self_params%km  , &
           cdt                   = self_params%cdt , &
           grav                  = grav            , &
           delp                  = self%delp       , &
           rc = istat                                &
                    )
 if(istat /=0) then
    call mpas_log_write('--- CA2G_oc_GridComp: error in subroutine phobicTophilic.', &
                        messageType=MPAS_LOG_CRIT)
 else
    call mpas_log_write('--- end subroutine phobicTophilic.')
 endif


!--- CA2G_oc settling:
 call mpas_log_write('--- enter subroutine Chem_Settling:')
 if(.not.allocated(qca2G)) allocate(qca2G(its:ite,jts:jte,kts:kte,self_params%nbins))
 do j = jts,jte
    do i = its,ite
       do k = kts,kte
          qca2G(i,j,k,1) = self%ocphobic(i,j,k)
          qca2G(i,j,k,2) = self%ocphilic(i,j,k)
       enddo
    enddo
 enddo

 do ibin = 1, self_params%nbins
    !if radius == 0, then we're dealing with a gas which has no settling losses:
    if(self_params%radius(ibin) == 0.0) cycle
    if(associated(self%ocSD)) self%ocSD(:,:,ibin) = 0._RKIND
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
              fluxout   = self%ocsd                      , &
              rc        = istat                            &
                      )
 enddo
 if(istat /=0) then
    call mpas_log_write('--- CA2G_oc_GridComp: error in subroutine Chem_Settling.', &
                        messageType=MPAS_LOG_CRIT)
 else
    call mpas_log_write('--- end subroutine Chem_Settling.')
 endif


!--- CA2G_oc dry deposition:
 call mpas_log_write('--- enter subroutine DryDeposition:')
 if(associated(self%ocDP)  ) self%ocdp(:,:,:) = 0._RKIND
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
    if(associated(self%ocDP)) then
       self%ocDP(:,:,ibin) = dqa*self%delp(:,:,self_params%km)/grav/self_params%cdt
    end if
 enddo
 if(istat /=0) then
    call mpas_log_write('--- CA2G_oc_GridComp: error in subroutine DryDeposition.', &
                        messageType=MPAS_LOG_CRIT)
 else
    if(allocated(drydepf)) deallocate(drydepf)
    call mpas_log_write('--- end subroutine DryDeposition.')
 endif


 do j = jts,jte
    do i = its,ite
       do k = kts,kte
          self%ocphobic(i,j,k) = qca2G(i,j,k,1)
          self%ocphilic(i,j,k) = qca2G(i,j,k,2)
       enddo
    enddo
 enddo


!--- CA2G_oc large-scale wet removal (hydrophobic mode is removed):
 call mpas_log_write('--- enter subroutine WetRemovalGOCART2G:')
 if(associated(self%ocWT)) self%ocwt(:,:,:) = 0._RKIND
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
           aero_type = 'OC'              , &
           kin       = KIN               , &
           grav      = grav              , &
           fwet      = fwet              , &
           aerosol   = self%ocphilic     , &
           ple       = self%ple          , &
           tmpu      = self%t            , &
           rhoa      = self%airdens      , &
           pfllsan   = self%pfl_lsan     , &
           pfilsan   = self%pfi_lsan     , &
           precc     = self%cn_prcp      , &
           precl     = self%ncn_prcp     , &
           fluxout   = self%ocwt         , &
           rc        = istat               &
                        )
 if(istat /=0) then
    call mpas_log_write('--- CA2G_oc_GridComp: error in subroutine WetRemovalGOCART2G.', &
                        messageType=MPAS_LOG_CRIT)
 else
    call mpas_log_write('--- end subroutine WetRemovalGOCART2G.')
 endif
!return


!--- CA2G_oc diagnostics:
 do j = jts,jte
    do i = its,ite
       do k = kts,kte
          qca2G(i,j,k,1) = self%ocphobic(i,j,k)
          qca2G(i,j,k,2) = self%ocphilic(i,j,k)
       enddo
    enddo
 enddo
 n_profile = size(self_params%wavelengths_profile)
 n_vertint = size(self_params%wavelengths_vertint)
 call mpas_log_write('--- enter subroutine Aero_Compute_Diags:')
 call mpas_log_write('--- nbins     = $i',intArgs=(/nbins/))
 call mpas_log_write('--- n_profile = $i',intArgs=(/n_profile/))
 call mpas_log_write('--- n_vertint = $i',intArgs=(/n_vertint/))
 if(associated(self%ocsmass)   ) self%ocsmass(:,:)       = 0._RKIND
 if(associated(self%occmass)   ) self%occmass(:,:)       = 0._RKIND
 if(associated(self%ocmass)    ) self%ocmass(:,:,:)      = 0._RKIND
 if(associated(self%ocexttau)  ) self%ocexttau(:,:,:)    = 0._RKIND
 if(associated(self%ocstexttau)) self%ocstexttau(:,:,:)  = 0._RKIND
 if(associated(self%ocscatau)  ) self%ocscatau(:,:,:)    = 0._RKIND
 if(associated(self%ocstscatau)) self%ocstscatau(:,:,:)  = 0._RKIND
 if(associated(self%ocfluxu)   ) self%ocfluxu(:,:)       = 0._RKIND
 if(associated(self%ocfluxv)   ) self%ocfluxv(:,:)       = 0._RKIND
 if(associated(self%occonc)    ) self%occonc(:,:,:)      = 0._RKIND
 if(associated(self%ocextcoef) ) self%ocextcoef(:,:,:,:) = 0._RKIND
 if(associated(self%ocscacoef) ) self%ocscacoef(:,:,:,:) = 0._RKIND
 if(associated(self%ocbckcoef) ) self%ocbckcoef(:,:,:,:) = 0._RKIND
 if(associated(self%ocangstr)  ) self%ocangstr(:,:)      = 0._RKIND
 if(associated(self%ocaeridx)  ) self%ocaeridx(:,:)      = 0._RKIND
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
              sfcmass             = self%ocsmass                           , &
              colmass             = self%occmass                           , &
              mass                = self%ocmass                            , &
              exttau              = self%ocexttau                          , &
              scatau              = self%ocscatau                          , &
!             stexttau            = self%ocstexttau                        , &
!             stscatau            = self%ocstscatau                        , &
              fluxu               = self%ocfluxu                           , &
              fluxv               = self%ocfluxv                           , &
              conc                = self%occonc                            , &
              extcoef             = self%ocextcoef                         , &
              scacoef             = self%ocscacoef                         , &
              bckcoef             = self%ocbckcoef                         , &
              angstrom            = self%ocangstr                          , &
              aerindx             = self%ocaeridx                          , &
              NO3nFlag            = .false.                                , &
              rc                  = istat                                    &
                        )
 if(istat /=0) then
    call mpas_log_write('--- CA2G_bc_GridComp: error in subroutine Aero_Compute_Diags.', &
                        messageType=MPAS_LOG_CRIT)
 else
    call mpas_log_write('--- end subroutine Aero_Compute_Diags.')
 endif
 if(allocated(qca2G)) deallocate(qca2G)
 return


 i1 = lbound(self%rh2,1); i2 = ubound(self%rh2,1)
 j1 = lbound(self%rh2,2); j2 = ubound(self%rh2,2)
 km = ubound(self%rh2,3)

 call mpas_log_write('--- enter subroutine Aero_Compute_Diags RH20:')
 if(associated(self%ocextcoefrh20)) self%ocextcoefrh20(:,:,:,:) = 0._RKIND
 if(associated(self%ocscacoefrh20)) self%ocscacoefrh20(:,:,:,:) = 0._RKIND
 if(.not.allocated(rh20)) allocate(rh20(i1:i2,j1:j2,km))
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
           extcoef             = self%ocextcoefrh20                     , &
           scacoef             = self%ocscacoefrh20                     , &
           NO3nFlag            = .false.                                , &
           rc                  = istat                                    &
                        )
 if(istat /=0) then
    call mpas_log_write('--- CA2G_bc_GridComp: error in subroutine Aero_Compute_Diags RH20.', &
                        messageType=MPAS_LOG_CRIT)
 else
    call mpas_log_write('--- end subroutine Aero_Compute_Diags RH20.')
 endif


 call mpas_log_write('--- enter subroutine Aero_Compute_Diags RH80:')
 if(associated(self%ocextcoefrh80)) self%ocextcoefrh80(:,:,:,:) = 0._RKIND
 if(associated(self%ocscacoefrh80)) self%ocscacoefrh80(:,:,:,:) = 0._RKIND
 if(.not.allocated(rh80)) allocate(rh80(i1:i2,j1:j2,km))
 rh80(:,:,:) = 0.80
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
           rh                  = rh80                                   , &
           u                   = self%u                                 , &
           v                   = self%v                                 , &
           delp                = self%delp                              , &
           ple                 = self%ple                               , &
            tropp               = self%tropp                            , &
           extcoef             = self%ocextcoefrh80                     , &
           scacoef             = self%ocscacoefRH80                     , &
           NO3nFlag            = .false.                                , &
           rc                  = istat                                    &
                        )
 if(istat /=0) then
    call mpas_log_write('--- CA2G_bc_GridComp: error in subroutine Aero_Compute_Diags RH80', &
                        messageType=MPAS_LOG_CRIT)
 else
    call mpas_log_write('--- end subroutine Aero_Compute_Diags RH80.')
 endif
 if(allocated(rh20)) deallocate(rh20)
 if(allocated(rh80)) deallocate(rh80)

 if(allocated(qca2G)) deallocate(qca2G)


 call mpas_log_write('--- end subroutine processes_CA2G_oc_GridComp.')

 end subroutine processes_CA2G_oc_GridComp

!=================================================================================================================
 end module CA2G_oc_GridCompMod
!=================================================================================================================
