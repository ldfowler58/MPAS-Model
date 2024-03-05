!=================================================================================================================
 module NI2G_state_SS
 use mpas_kind_types,only: RKIND
 public
 save

!this module is the state variable specification file for nitrate parameters. it is the same as NI2G_StateSpecs.rc
!in the GOCART-2G directory ./GOCART-2G/ESMF/GOCART2G_GridComp/NI2G_GridComp.

!schema_version: 2.0.0
!component: NI


 type NI2G_State

!category: IMPORT
 real(kind=RKIND),dimension(:,:),allocatable    :: LWI           ! land-ocean-ice_mask (-)
 real(kind=RKIND),dimension(:,:),allocatable    :: TROPP         ! tropopause_pressure_based_on_blended_estimate (Pa)
 real(kind=RKIND),dimension(:,:),allocatable    :: USTAR         ! surface_velocity_scale (m s-1)
 real(kind=RKIND),dimension(:,:),allocatable    :: ZPBL          ! planetary_boundary_layer_height (m)
 real(kind=RKIND),dimension(:,:),allocatable    :: SH            ! sensible_heat_flux_from_turbulence (W m-2)
 real(kind=RKIND),dimension(:,:),allocatable    :: Z0H           !  surface_roughness_for_heat (m)
 real(kind=RKIND),dimension(:,:),allocatable    :: CN_PRCP       ! surface_conv._rain_flux_needed_by_land (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),allocatable    :: NCN_PRCP      ! Non-convective precipitation (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),allocatable    :: AREA          ! agrid_cell_area (^2)
!.................................................................................................................
 real(kind=RKIND),dimension(:,:,:),allocatable  :: AIRDENS       ! moist_air_density (kg/m^3)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: DELP          ! pressure_thickness (Pa)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: T             ! air_temperature (K)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: RH2           ! Rel_Hum_after_moist (-)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: ZLE           ! geopotential_height (m)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: PLE           ! air_pressure (Pa)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: PFL_LSAN      ! 3D_flux_of_liquid_nonconvective_precipitation (kg/m2/s)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: PFI_LSAN      ! 3D_flux_of_ice_nonconvective_precipitation (kg/m2/s)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: U             ! eastward_wind (m s-1)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: V             ! northward_wind (m s-1)
!.................................................................................................................
 real(kind=RKIND),dimension(:,:),allocatable    :: EMI_NH3_AG    ! agriculture emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),allocatable    :: EMI_NH3_BB    ! biomass burning emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),allocatable    :: EMI_NH3_EN    ! energy emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),allocatable    :: EMI_NH3_IN    ! industry emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),allocatable    :: EMI_NH3_OC    ! ocean emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),allocatable    :: EMI_NH3_RE    ! resedential emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),allocatable    :: EMI_NH3_TR    ! transport emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),allocatable    :: NITRATE_HNO3  ! nitrate hno3 emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: DU            ! Dust Mixing Ratio all bins (kg kg-1)
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: SS            ! Sea Salt Mixing Ratio all bins (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: SO4           ! Sulfate Mixing Ratio (kg kg-1)

!category: EXPORT
 real(kind=RKIND),dimension(:,:,:),allocatable  :: NH3MASS       ! Ammonia Mass Mixing Ratio (kg/kg)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: NH4MASS       ! Ammonium Aerosol Mass Mixing Ratio (kg/kg)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: NIMASS        ! Nitrate Mass Mixing Ratio (kg/kg)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: NIMASS25      ! Nitrate Mass Mixing Ratio [PM2.5] (kg/kg)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: HNO3CONC      ! Nitric Acid Mass Concentration (kg m-3)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: NH3CONC       ! Ammonia Mass Concentration (kg m-3)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: NH4CONC       ! Ammonium Mass Concentration (kg m-3)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: NICONC        ! Nitrate Mass Concentration (kg m-3)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: NICONC25      ! Nitrate Mass Concentration [PM2.5] (kg m-3)
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: NIEXTCOEF     ! Nitrate Extinction Coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: NIEXTCOEFRH20 ! Nitrate Extinction Coefficient - fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: NIEXTCOEFRH80 ! Nitrate Extinction Coefficient - fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: NISCACOEF     ! Nitrate Scattering Coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: NISCACOEFRH20 ! Nitrate Scattering Coefficient - fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: NISCACOEFRH80 ! Nitrate Scattering Coefficient - fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: NIBCKCOEF     ! Nitrate Backscatter Coefficient (m-1 sr-1)
!.................................................................................................................
 real(kind=RKIND),dimension(:,:),allocatable    :: NIPNO3AQ      ! Nitrate Production from Aqueous Chemistry (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),allocatable    :: NIPNH4AQ      ! Ammonium Production from Aqueous Chemistry (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),allocatable    :: NIPNH3AQ      ! Ammonia Change from Aqueous Chemistry (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: NIHT          ! Nitrate Production from Het Chem (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: NISD          ! Nitrate Sedimentation (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: NIDP          ! Nitrate Dry Deposition (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: NIWT          ! Nitrate Wet Deposition (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: NISV          ! Nitrate Convective Scavenging (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),allocatable    :: NH3EM         ! Ammonia Emission (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),allocatable    :: NH3DP         ! Ammonia Dry Deposition (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),allocatable    :: NH3WT         ! Ammonia Wet Deposition (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),allocatable    :: NH3SV         ! Ammonia Convective Scavenging (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),allocatable    :: NH4SD         ! Ammonium Settling (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),allocatable    :: NH4DP         ! Ammonium Dry Deposition (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),allocatable    :: NH4WT         ! Ammonium Wet Deposition (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),allocatable    :: NH4SV         ! Ammonium Convective Scavenging (kg m-2 s-1)

 real(kind=RKIND),dimension(:,:),allocatable    :: HNO3SMASS     ! Nitric Acid Surface Mass Concentration (kg m-3)
 real(kind=RKIND),dimension(:,:),allocatable    :: NH3SMASS      ! Ammonia Surface Mass Concentration (kg m-3)
 real(kind=RKIND),dimension(:,:),allocatable    :: NH4SMASS      ! Ammonium Surface Mass Concentration (kg m-3)
 real(kind=RKIND),dimension(:,:),allocatable    :: NISMASS       ! Nitrate Surface Mass Concentration (kg m-3)
 real(kind=RKIND),dimension(:,:),allocatable    :: NISMASS25     ! Nitrate Surface Mass Concentration [PM2.5] (kg m-3)
 real(kind=RKIND),dimension(:,:),allocatable    :: HNO3CMASS     ! Nitric Acid Column Mass Density (kg m-3)
 real(kind=RKIND),dimension(:,:),allocatable    :: NH3CMASS      ! Ammonia Column Mass Density (kg m-3)
 real(kind=RKIND),dimension(:,:),allocatable    :: NH4CMASS      ! Ammonium Column Mass Density (kg m-3)
 real(kind=RKIND),dimension(:,:),allocatable    :: NICMASS       ! Nitrate Column Mass Density (kg m-2)
 real(kind=RKIND),dimension(:,:),allocatable    :: NICMASS25     ! Nitrate Column Mass Density [PM2.5] (kg m-2)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: NIEXTTFM      ! Nitrate Extinction AOT - PM 1.0 um (-)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: NISCATFM      ! Nitrate Scattering AOT - PM 1.0 um (-)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: NIEXTT25      ! Nitrate Extinction AOT - PM 2.5 um (-)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: NISCAT25      ! Nitrate Scattering AOT - PM 2.5 um (-)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: NIEXTTAU      ! Nitrate Extinction AOT (-)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: NISTEXTTAU    ! Nitrate Extinction AOT Stratosphere (-)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: NISCATAU      ! Nitrate Scattering AOT (-)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: NISTSCATAU    ! Nitrate Scattering AOT Stratosphere (-)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: NIANGSTR      ! Nitrate Angstrom parameter [470-870 nm] (-)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: NIFLUXU       ! Nitrate column u-wind mass flux (kg m-1 s-1)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: NIFLUXV       ! Nitrate column v-wind mass flux (kg m-1 s-1)

!category: INTERNAL
 real(kind=RKIND),dimension(:,:,:),allocatable  :: NH3           ! Ammonia (NH3, gas phase) (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: NH4a          ! Ammonium ion (NH4+, aerosol phase) (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: NO3an1        ! Nitrate size bin 001 (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: NO3an2        ! Nitrate size bin 002 (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: NO3an3        ! Nitrate size bin 003 (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: XHNO3         ! buffer for NITRATE_HNO3 (kg m-2 s-1)

 end type NI2G_State

!=================================================================================================================
 end module NI2G_state_SS
!=================================================================================================================
