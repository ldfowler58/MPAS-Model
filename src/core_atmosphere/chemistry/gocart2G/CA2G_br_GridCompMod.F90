! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!==================================================================================================================
 module CA2G_br_GridCompMod
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

 use CA2G_br_instance,only: nbins,particle_radius_microns,particle_density,fscav,molecular_weight,   &
                            fnum,rhFlag,pressure_lid_in_hPa,sigma,hydrophobic_fraction,pom_ca_ratio, &
                            aircraft_fuel_emission_factor,aviation_vertical_layers,                  &
                            point_emissions_srcfilen
 use CA2G_br_StateSpecs,only: CA2G_br_State


 implicit none
 private


!--- constants (these parameters need to be accessed from MPAS physics instead of redefined here):
 real(kind=RKIND),parameter:: cpd      = 1003.0_RKIND
 real(kind=RKIND),parameter:: grav     = 9.80616_RKIND
 real(kind=RKIND),parameter:: pi       = 3.141592653589793_RKIND
 real(kind=RKIND),parameter:: karman   = 0.4_RKIND
 real(kind=RKIND),parameter:: radTodeg = 180._RKIND/pi
 real(kind=RKIND),parameter:: undefval = 1.0e15


!--- types needed to define CA2G_br:
 type :: ThreadWorkspace
    integer:: nPts = -1
    integer,dimension(:),allocatable:: pstart,pend

    real(kind=RKIND),dimension(:),allocatable:: pLat,pLon,pBase,pTop,pEmis
 end type ThreadWorkspace

 type,extends(GA_Environment),public:: CA2G_br_GridComp
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
       procedure:: emissions_GridComp => emissions_CA2G_br_GridComp
       procedure:: load_GridComp      => load_CA2G_br_GridComp
       procedure:: processes_GridComp => processes_CA2G_br_GridComp
 end type CA2G_br_GridComp

 type wrap_
    type(CA2G_br_GridComp),pointer:: PTR !=> null()
 end type wrap_


 contains


!==================================================================================================================
 subroutine load_CA2G_br_GridComp(self,kts,kte)
!==================================================================================================================

!--- input arguments:
 integer,intent(in):: kts,kte

!--- inout arguments:
 class(CA2G_br_GridComp),intent(inout) :: self

!local variables:
 integer:: n

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write('--- enter subroutine load_CA2G_br_GridComp:')


!--- initialization of vertical index:
 self%klid = kts
 self%km = kte-kts+1


!--- initialization of variables in GA_Environment:
 self%diurnal_bb = .false.

 call self%load_from_config(nbins,particle_radius_microns,particle_density,fscav,molecular_weight,fnum, &
                            rhFlag,pressure_lid_in_hPa)

 call mpas_log_write('--- nbins = $i',intArgs=(/self%nbins/))
 call mpas_log_write('--- radius,rhop,fscav,molwght,fnum:')
 do n = 1,self%nbins
    call mpas_log_write('$i $r $r $r $r $r',intArgs=(/n/),realArgs=(/self%radius(n),self%rhop(n), &
                        self%fscav(n),self%molwght(n),self%fnum(n)/))
 enddo


!--- initialization of all other variables in CA2G_br_GridComp:
 self%fHydrophobic = hydrophobic_fraction
 self%ratPOM       = pom_ca_ratio
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

 end subroutine load_CA2G_br_GridComp

!==================================================================================================================
 subroutine emissions_CA2G_br_GridComp(self_params,self,its,ite,jts,jte,kts,kte,iyr,imm,idd,ihr,imn,isc)
!==================================================================================================================

!--- input arguments:
 integer,intent(in):: its,ite,jts,jte,kts,kte
 integer,intent(in):: iyr,imm,idd,ihr,imn,isc

 class(CA2G_br_GridComp),intent(in):: self_params

!--- inout arguments:
 class(CA2G_br_State),intent(inout):: self

!--- local variables and arrays:
 integer:: nymd,nhms
 integer:: istat

 real(kind=RKIND),dimension(:,:),pointer:: biomass_src         !
 real(kind=RKIND),dimension(:,:),pointer:: biofuel_src         !
 real(kind=RKIND),dimension(:,:),pointer:: eocant1_src         ! 
 real(kind=RKIND),dimension(:,:),pointer:: eocant2_src         !
 real(kind=RKIND),dimension(:,:),pointer:: oc_ship_src         !
 real(kind=RKIND),dimension(:,:),pointer:: aviation_lto_src    !
 real(kind=RKIND),dimension(:,:),pointer:: aviation_cds_src    !
 real(kind=RKIND),dimension(:,:),pointer:: aviation_crs_src    !
 real(kind=RKIND),dimension(:,:,:),pointer:: aircraft_fuel_src !

 real(kind=RKIND),dimension(:,:),pointer:: biogvoc_src  => null() !
 real(kind=RKIND),dimension(:,:),pointer:: biomass_src_ => null() !

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine emissions_CA2G_br_GridComp:')


!--- extract nymd and nhms from clock:
 call MAPL_PackTime(nymd,iyr,imm,idd)
 call MAPL_PackTime(nhms,ihr,imn,isc)
 call mpas_log_write('--- nymd = $i',intArgs=(/nymd/))
 call mpas_log_write('--- nhms = $i',intArgs=(/nhms/))


!--- initialize brown carbon emissions: all other emission types are set to zero. 
 biomass_src       => self%br_biomass
 biofuel_src       => self%br_biofuel
 eocant1_src       => self%br_antebr1
 eocant2_src       => self%br_antebr2
 oc_ship_src       => self%br_ship
 aviation_lto_src  => self%br_aviation_lto
 aviation_cds_src  => self%br_aviation_cds
 aviation_crs_src  => self%br_aviation_crs
 aircraft_fuel_src => self%br_aircraft

 biomass_src(:,:)         = 0._RKIND
 biofuel_src(:,:)         = 0._RKIND
 eocant1_src(:,:)         = 0._RKIND
 eocant2_src(:,:)         = 0._RKIND
 oc_ship_src(:,:)         = 0._RKIND
 aviation_lto_src(:,:)    = 0._RKIND
 aviation_cds_src(:,:)    = 0._RKIND
 aviation_crs_src(:,:)    = 0._RKIND
 aircraft_fuel_src(:,:,:) = 0._RKIND

 if(.not.associated(biogvoc_src)) allocate(biogvoc_src(its:ite,jts:jte))
 biogvoc_src(:,:) = self%br_terpene


!--- as a safety check, all undefined values are set to zero. this may be needed when all emission types become
!    available.
 where(1.01*biomass_src > undefval) biomass_src = 0._RKIND
 where(1.01*biogvoc_src > undefval) biogvoc_src = 0._RKIND
 where(1.01*biofuel_src > undefval) biofuel_src = 0._RKIND
 where(1.01*eocant1_src > undefval) eocant1_src = 0._RKIND
 where(1.01*eocant2_src > undefval) eocant2_src = 0._RKIND
 where(1.01*oc_ship_src > undefval) oc_ship_src = 0._RKIND
 where(1.01*aviation_lto_src  > undefval) aviation_lto_src  = 0._RKIND
 where(1.01*aviation_cds_src  > undefval) aviation_cds_src  = 0._RKIND
 where(1.01*aviation_crs_src  > undefval) aviation_crs_src  = 0._RKIND
 where(1.01*aircraft_fuel_src > undefval) aircraft_fuel_src = 0._RKIND


!--- apply diurnal cycle to biomass burning if needed:
 if(self_params%diurnal_bb ) then
    call mpas_log_write('--- enter subroutine Chem_Biomass Diurnal:')
    biomass_src_ => biomass_src
    call Chem_BiomassDiurnal( &
       cdt  = self_params%cdt,    &
       nhms = nhms,               &
       eout = biomass_src,        &
       ein  = biomass_src_,       &
       lons = self%lons*radTodeg, &
       lats = self%lats*radTodeg  &
                            )
    call mpas_log_write('--- end subroutine Chem_Biomass Diurnal:')
 endif


!--- apply emissions to CA2G_br: 
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
    aerosolPhilic     = self%brphobic,               &
    aerosolPhobic     = self%brphilic,               &
    oc_emis           = self%brem,                   &
    oc_emisan         = self%breman,                 &
    oc_emisbb         = self%brembb,                 &
    oc_emisbf         = self%brembf,                 &
    oc_emisbg         = self%brembg,                 &
    rc                = istat                        &
                )
 if(istat /=0) then
    call mpas_log_write('--- CA2G_br_GridComp: error in subroutine CAEmission', &
                        messageType=MPAS_LOG_CRIT)
 else
    call mpas_log_write('--- end subroutine CAEmission:')
 endif


!--- for now, we do not support point emissions:


 call mpas_log_write('--- end subroutine emissions_CA2G_br_GridComp:')

 end subroutine emissions_CA2G_br_GridComp

!==================================================================================================================
 subroutine processes_CA2G_br_GridComp(self_params,self,its,ite,jts,jte,kts,kte)
!==================================================================================================================
!--- input arguments:

 integer,intent(in):: its,ite,jts,jte,kts,kte

!--- inout arguments:
 class(CA2G_br_GridComp),intent(inout):: self_params
 class(CA2G_br_State),intent(inout):: self

!------------------------------------------------------------------------------------------------------------------
 call mpas_log_write(' ')
 call mpas_log_write('--- enter subroutine processes_CA2G_br_GridComp:')


 call mpas_log_write('--- end subroutine processes_CA2G_br_GridComp:')

 end subroutine processes_CA2G_br_GridComp

!==================================================================================================================
 end module CA2G_br_GridCompMod
!==================================================================================================================
