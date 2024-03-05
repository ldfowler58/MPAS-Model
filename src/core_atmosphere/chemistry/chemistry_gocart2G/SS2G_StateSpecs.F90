!=================================================================================================================
 module SS2G_state_SS
 use mpas_kind_types,only: RKIND
 public
 save

!this module is the state variable specification file for sea-salt parameters. it is the same as SS2G_StateSpecs.rc
!in the GOCART-2G directory ./GOCART-2G/ESMF/GOCART2G_GridComp/SG2G_GridComp.

!schema_version: 2.0.0
!component: SS


 type SS2G_State

!category: IMPORT
 real(kind=RKIND),dimension(:,:),allocatable  :: frocean  !fraction_of_ocean (1)
 real(kind=RKIND),dimension(:,:),allocatable  :: fraci    !ice_covered_fraction_of_tile (1)
 real(kind=RKIND),dimension(:,:),allocatable  :: lwi      !land-ocean-ice_mask (1)
 real(kind=RKIND),dimension(:,:),allocatable  :: tropp    !tropopause_pressure_based_on_blended_estimate (Pa)
 real(kind=RKIND),dimension(:,:),allocatable  :: u10m     !10-meter_eastward_wind (m s-1)
 real(kind=RKIND),dimension(:,:),allocatable  :: v10m     !10-meter_northward_wind (m s-1)
 real(kind=RKIND),dimension(:,:),allocatable  :: ustar    !surface_velocity_scale (m s-1)
 real(kind=RKIND),dimension(:,:),allocatable  :: ts       !surface skin temperature (K)
 real(kind=RKIND),dimension(:,:),allocatable  :: dz       !surface_layer_height (m)
 real(kind=RKIND),dimension(:,:),allocatable  :: frlake   !fraction_of_lake (1)
 real(kind=RKIND),dimension(:,:),allocatable  :: area     !grid_cell_area (m^2)
 real(kind=RKIND),dimension(:,:),allocatable  :: zpbl     !planetary_boundary_layer_height (m)
 real(kind=RKIND),dimension(:,:),allocatable  :: sh       !sensible_heat_flux_from_turbulence (W m-2)
 real(kind=RKIND),dimension(:,:),allocatable  :: z0h      !surface_roughness_for_heat(m)
 real(kind=RKIND),dimension(:,:),allocatable  :: cn_prcp  !surface_conv._rain_flux_needed_by_land (kg/m^2/s)
 real(kind=RKIND),dimension(:,:),allocatable  :: ncn_prcp !non-convective precipitation (kg/m^2/s)
!.................................................................................................................
 real(kind=RKIND),dimension(:,:,:),allocatable:: airdens  !moist_air_density (kg/m^3)
 real(kind=RKIND),dimension(:,:,:),allocatable:: delp     !pressure_thickness (Pa)
 real(kind=RKIND),dimension(:,:,:),allocatable:: t        !air_temperature (K)
 real(kind=RKIND),dimension(:,:,:),allocatable:: rh2      !rel_hum_after_moist (1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: zle      !geopotential_height (m)
 real(kind=RKIND),dimension(:,:,:),allocatable:: ple      !air_pressure (Pa)
 real(kind=RKIND),dimension(:,:,:),allocatable:: pfl_lsan !3d_flux_of_liquid_nonconvective_precipitation (kg/m2s)
 real(kind=RKIND),dimension(:,:,:),allocatable:: pfi_lsan !3d_flux_of_ice_nonconvective_precipitation (kg/m2/s)
 real(kind=RKIND),dimension(:,:,:),allocatable:: u        !eastward_wind (m s-1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: v        !northward_wind (m s-1)

!EXPORT: category
 real(kind=RKIND),dimension(:,:,:),allocatable  :: SSMASS        !Sea Salt Mass Mixing Ratio (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: SSMASS25      !Sea Salt Mass Mixing Ratio - PM 2.5 (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: SSCONC        !Sea Salt Mass Concentration (kg m-3)
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: SSEXTCOEF     !Sea Salt Extinction Coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: SSEXTCOEFRH20 !Sea Salt Extinction Coefficient - Fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: SSEXTCOEFRH80 !Sea Salt Extinction Coefficient - Fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: SSSCACOEF     !Sea Salt Scattering Coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: SSSCACOEFRH20 !Sea Salt Scattering Coefficient - Fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: SSSCACOEFRH80 !Sea Salt Scattering Coefficient - Fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: SSBCKCOEF     !Sea Salt Backscatter Coefficient (m-1 sr-1)
!.................................................................................................................
 real(kind=RKIND),dimension(:,:,:),allocatable  :: SSEM          !Sea Salt Emission (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: SSSD          !Sea Salt Sedimentation (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: SSDP          !Sea Salt Dry Deposition (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: SSWT          !Sea Salt Wet Deposition (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: SSSV          !Sea Salt Convective Scavenging (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),allocatable    :: SSSMASS       !Sea Salt Surface Mass Concentration (kg m-3)
 real(kind=RKIND),dimension(:,:),allocatable    :: SSCMASS       !Sea Salt Column Mass Density (kg m-2)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: SSEXTTAU      !Sea Salt Extinction AOT (-)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: SSSTEXTTAU    !Sea Salt Extinction AOT Stratosphere (-)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: SSSCATAU      !Sea Salt Scattering AOT (-)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: SSSTSCATAU    !Sea Salt Scattering AOT Stratosphere (-)
 real(kind=RKIND),dimension(:,:),allocatable    :: SSSMASS25     !Sea Salt Surface Mass Concentration - PM 2.5 (kg m-3)
 real(kind=RKIND),dimension(:,:),allocatable    :: SSCMASS25     !Sea Salt Column Mass Density - PM 2.5 (kg m-2)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: SSEXTT25      !Extinction AOT - PM 2.5 (-)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: SSSCAT25      !Sea Salt Scattering AOT - PM 2.5 (-)
 real(kind=RKIND),dimension(:,:),allocatable    :: SSAERIDX      !Sea Salt TOMS UV Aerosol Index (-)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: SSEXTTFM      !Sea Salt Extinction AOT [550 nm] - PM 1.0 um (-)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: SSSCATFM      !Sea Salt Scattering AOT [550 nm] - PM 1.0 um (-)
 real(kind=RKIND),dimension(:,:),allocatable    :: SSANGSTR      !Sea Salt Angstrom parameter [470-870 nm] (-)
 real(kind=RKIND),dimension(:,:),allocatable    :: SSFLUXU       !Sea Salt column u-wind mass flux (kg m-1 s-1)
 real(kind=RKIND),dimension(:,:),allocatable    :: SSFLUXV       !Sea Salt column v-wind mass flux (kg m-1 s-1)

!category: INTERNAL
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: SS              !Sea Salt Mixing Ratio (bin %d) (kg kg-1)
 real(kind=RKIND),dimension(:,:),allocatable    :: DEEP_LAKES_MASK !Deep Lakes Mask

 end type SS2G_State

!=================================================================================================================
 end module SS2G_state_SS
!=================================================================================================================
