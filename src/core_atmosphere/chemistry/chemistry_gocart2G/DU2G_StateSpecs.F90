!=================================================================================================================
 module DU2G_StateSpecs
 use mpas_kind_types,only: RKIND
 public
 save

!this module is the state variable specification file for dust parameters. it is the same as DU2G_StateSpecs.rc
!in the GOCART-2G directory ./GOCART-2G/ESMF/GOCART2G_GridComp/DU2G_GridComp.

!schema_version: 2.0.0
!component: DU


 type DU2G_State

 real(kind=RKIND),dimension(:,:),allocatable  :: DU_SRC     ! erod - dust emissions (-)
 real(kind=RKIND),dimension(:,:),allocatable  :: DU_Z0      ! aerodynamic_surface_roughness_for_aeolian_processes (-)
 real(kind=RKIND),dimension(:,:),allocatable  :: DU_GVF     ! GVF (-)
 real(kind=RKIND),dimension(:,:),allocatable  :: DU_SAND    ! volume_fraction_of_sand_in_soil (-)
 real(kind=RKIND),dimension(:,:),allocatable  :: DU_SILT    ! volume_fraction_of_silt_in_soil (-)
 real(kind=RKIND),dimension(:,:),allocatable  :: DU_CLAY    ! volume_fraction_of_clay_in_soil (-)
 real(kind=RKIND),dimension(:,:),allocatable  :: DU_RDRAG   ! drag_partition (m-1)
 real(kind=RKIND),dimension(:,:),allocatable  :: DU_SSM     ! sediment_supply_map (-)
 real(kind=RKIND),dimension(:,:),allocatable  :: DU_UTHRES  ! surface_dry_threshold_velocity (m s-1)
 real(kind=RKIND),dimension(:,:),allocatable  :: FRSNOW     ! surface_snow_area_fraction (-)
 real(kind=RKIND),dimension(:,:),allocatable  :: SLC        ! liquid_water_content_of_soil_layer (-)
 real(kind=RKIND),dimension(:,:),allocatable  :: DU_TEXTURE ! soil_texture (-)
 real(kind=RKIND),dimension(:,:),allocatable  :: DU_VEG     ! vegetation_type (-)
 real(kind=RKIND),dimension(:,:),allocatable  :: FRLAKE     ! fraction_of_lake (-)
 real(kind=RKIND),dimension(:,:),allocatable  :: FRLAND     ! fraction_of_land (-)
 real(kind=RKIND),dimension(:,:),allocatable  :: ASNOW      ! snow_covered_fraction_of_land (-)
 real(kind=RKIND),dimension(:,:),allocatable  :: WET1       ! surface_soil_wetness (-)
 real(kind=RKIND),dimension(:,:),allocatable  :: LWI        ! land-ocean-ice_mask (-)
 real(kind=RKIND),dimension(:,:),allocatable  :: TROPP      ! tropopause_pressure_based_on_blended_estimate (Pa)
 real(kind=RKIND),dimension(:,:),allocatable  :: U10M       ! 10-meter_eastward_wind (m s-1)
 real(kind=RKIND),dimension(:,:),allocatable  :: V10M       ! 10-meter_northward_wind (m s-1)
 real(kind=RKIND),dimension(:,:),allocatable  :: U10N       ! equivalent_neutral_10-meter_eastward_wind (m s-1)
 real(kind=RKIND),dimension(:,:),allocatable  :: V10N       ! equivalent_neutral_10-meter_northward_wind (m s-1)
 real(kind=RKIND),dimension(:,:),allocatable  :: AREA       ! agrid_cell_area (m^2)
 real(kind=RKIND),dimension(:,:),allocatable  :: USTAR      ! equivalent_neutral_10-meter_northward_wind (m s-1)
 real(kind=RKIND),dimension(:,:),allocatable  :: CN_PRCP    ! surface_conv._rain_flux_needed_by_land (kg/m^2/s)
 real(kind=RKIND),dimension(:,:),allocatable  :: NCN_PRCP   ! Non-convective precipitation (kg/m^2/s)
 real(kind=RKIND),dimension(:,:),allocatable  :: ZPBL       ! planetary_boundary_layer_height (m)
 real(kind=RKIND),dimension(:,:),allocatable  :: SH         ! sensible_heat_flux_from_turbulence (W m-2)
 real(kind=RKIND),dimension(:,:),allocatable  :: Z0H        ! surface_roughness_for_heat (m)
 real(kind=RKIND),dimension(:,:),allocatable  :: WCSF       ! water_surface_layer (m3 m-3)
 real(kind=RKIND),dimension(:,:),allocatable  :: TSOIL1     ! soil_temperatures_layer_1 (k)
 real(kind=RKIND),dimension(:,:),allocatable  :: RHOS       ! air_density_at_surface (kg m-3)
!.................................................................................................................
 real(kind=RKIND),dimension(:,:,:),allocatable:: AIRDENS    ! moist_air_density (kg/m^3)
 real(kind=RKIND),dimension(:,:,:),allocatable:: DELP       ! pressure_thickness (Pa)
 real(kind=RKIND),dimension(:,:,:),allocatable:: RH2        ! Rel_Hum_after_moist (-)
 real(kind=RKIND),dimension(:,:,:),allocatable:: T          ! air_temperature (K)
 real(kind=RKIND),dimension(:,:,:),allocatable:: ZLE        ! geopotential_height (m)
 real(kind=RKIND),dimension(:,:,:),allocatable:: PLE        ! air_pressure (Pa)
 real(kind=RKIND),dimension(:,:,:),allocatable:: PFL_LSAN   ! 3D_flux_of_liquid_nonconvective_precipitation (kg/m2/s)
 real(kind=RKIND),dimension(:,:,:),allocatable:: PFI_LSAN   ! 3D_flux_of_ice_nonconvective_precipitation (kg/m2/s)
 real(kind=RKIND),dimension(:,:,:),allocatable:: U          ! eastward_wind (m s-1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: V          ! northward_wind (m s-1)

!category: EXPORT
 real(kind=RKIND),dimension(:,:,:),allocatable  ::DUMASS        ! Dust Mass Mixing Ratio (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),allocatable  ::DUMASS25      ! Dust Mass Mixing Ratio (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),allocatable  ::DUCONC        ! Dust Mass Concentration (kg m-3)
 real(kind=RKIND),dimension(:,:,:,:),allocatable::DUEXTCOEF     ! Dust Extinction Coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),allocatable::DUEXTCOEFRH20 ! Dust Extinction Coefficient - Fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),allocatable::DUEXTCOEFRH80 ! Dust Extinction Coefficient - Fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),allocatable::DUSCACOEF     ! Dust Scattering Coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),allocatable::DUSCACOEFRH20 ! Dust Scattering Coefficient - Fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),allocatable::DUSCACOEFRH80 ! Dust Scattering Coefficient - Fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),allocatable::DUBCKCOEF     ! Dust Backscatter Coefficient (m-1 sr-1)
!.................................................................................................................
 real(kind=RKIND),dimension(:,:),allocatable  :: DUSMASS        ! Dust Surface Mass Concentration (kg m-3)
 real(kind=RKIND),dimension(:,:),allocatable  :: DUCMASS        ! Dust Column Mass Density (kg m-2)
 real(kind=RKIND),dimension(:,:,:),allocatable:: DUEXTTAU       ! Dust Extinction AOT (-)
 real(kind=RKIND),dimension(:,:,:),allocatable:: DUSTEXTTAU     ! Dust Extinction AOT Stratosphere (-)
 real(kind=RKIND),dimension(:,:,:),allocatable:: DUSCATAU       ! Dust Scattering AOT (-)
 real(kind=RKIND),dimension(:,:,:),allocatable:: DUSTSCATAU     ! Dust Scattering AOT Stratosphere (-)
 real(kind=RKIND),dimension(:,:),allocatable  :: DUSMASS25      ! Dust Surface Mass Concentration - PM 2.5 (kg m-3)
 real(kind=RKIND),dimension(:,:),allocatable  :: DUCMASS25      ! Dust Column Mass Density - PM 2.5 (kg m-2)
 real(kind=RKIND),dimension(:,:,:),allocatable:: DUEXTT25       ! Dust Extinction AOT - PM 2.5 (-)
 real(kind=RKIND),dimension(:,:,:),allocatable:: DUSCAT25       ! Dust Scattering AOT - PM 2.5 (-)
 real(kind=RKIND),dimension(:,:),allocatable  :: DUAERIDX       ! Dust TOMS UV Aerosol Index (-)
 real(kind=RKIND),dimension(:,:),allocatable  :: DUFLUXU        ! Dust column u-wind mass flux (kg m-1 s-1)
 real(kind=RKIND),dimension(:,:),allocatable  :: DUFLUXV        ! Dust column v-wind mass flux (kg m-1 s-1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: DUEXTTFM       ! Dust Extinction AOT - PM 1.0 um (-)
 real(kind=RKIND),dimension(:,:,:),allocatable:: DUSCATFM       ! Dust Scattering AOT - PM 1.0 um (-)
 real(kind=RKIND),dimension(:,:),allocatable  :: DUANGSTR       ! Dust Angstrom parameter [470-870 nm] (-)
 real(kind=RKIND),dimension(:,:,:),allocatable:: DUEM           ! Dust Emission (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: DUSD           ! Dust Sedimentation (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: DUDP           ! Dust Dry Deposition (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: DUWT           ! Dust Wet Deposition (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: DUSV           ! Dust Convective Scavenging (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),allocatable  :: DU_UST         ! aeolian_friction_velocity (-)
 real(kind=RKIND),dimension(:,:),allocatable  :: DU_UST_T       ! aeolian_threshold_friction_velocity (-)
 real(kind=RKIND),dimension(:,:),allocatable  :: DU_UST_TS      ! aeolian_threshold_friction_velocity_over_smooth_surface (-)
 real(kind=RKIND),dimension(:,:),allocatable  :: DU_DPC         ! aeolian_drag_partition_correction (-)
 real(kind=RKIND),dimension(:,:),allocatable  :: DU_SMC         ! aeolian_soil_moisture_correction (-)
 real(kind=RKIND),dimension(:,:),allocatable  :: DU_EROD        ! aeolian_erodibilitiy (-)

!category: INTERNAL
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: DU           ! Dust Mixing Ratio (Bin %d) (kg kg-1)

 end type DU2G_State

!=================================================================================================================
 end module DU2G_StateSpecs
!=================================================================================================================

