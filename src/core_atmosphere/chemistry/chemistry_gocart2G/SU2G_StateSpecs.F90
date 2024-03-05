!=================================================================================================================
 module SU2G_StateSpecs
 use mpas_kind_types,only: RKIND
 public
 save

!this module is the state variable specification file for sulfur parameters. it is the same as SU2G_StateSpecs.rc
!in the GOCART-2G directory ./GOCART-2G/ESMF/GOCART2G_GridComp/SU2G_GridComp.

!schema_version: 2.0.0
!component: SU


 type SU2G_State
 
!category: IMPORT
 real(kind=RKIND),dimension(:,:),allocatable:: frocean    !fraction_of_ocean (1)
 real(kind=RKIND),dimension(:,:),allocatable:: lwi        !land-ocean-ice_mask (1)
 real(kind=RKIND),dimension(:,:),allocatable:: tropp      !tropopause_pressure_based_on_blended_estimate (Pa)
 real(kind=RKIND),dimension(:,:),allocatable:: u10m       !10-meter_eastward_wind (m s-1)
 real(kind=RKIND),dimension(:,:),allocatable:: v10m       !10-meter_northward_wind (m s-1)
 real(kind=RKIND),dimension(:,:),allocatable:: area       !grid_cell_area (m^2)
 real(kind=RKIND),dimension(:,:),allocatable:: zpbl       !planetary_boundary_layer_height (m)
 real(kind=RKIND),dimension(:,:),allocatable:: ustar      !surface_velocity_scale (m s-1)
 real(kind=RKIND),dimension(:,:),allocatable:: sh         !sensible_heat_flux_from_turbulence (W m-2)
 real(kind=RKIND),dimension(:,:),allocatable:: z0h        !surface_roughness_for_heat (m)
 real(kind=RKIND),dimension(:,:),allocatable:: cn_prcp    !surface_conv._rain_flux_needed_by_land (kg/m^2/s)
 real(kind=RKIND),dimension(:,:),allocatable:: ncn_prcp   !non-convective precipitation (kg/m^2/s)
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
 real(kind=RKIND),dimension(:,:,:),allocatable:: fcld     !cloud fraction for radiation (1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: pSO2_OCS      !source species (1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: SU_AIRCRAFT   !fuel source species (1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: SU_NO3        !climatological NO3 source (1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: SU_OH         !climatological OH source (1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: SU_H2O2       !climatological H2O2 source (1)
!.................................................................................................................
 real(kind=RKIND),dimension(:,:),allocatable:: SU_BIOMASS      !biomass burning emissions (1)
 real(kind=RKIND),dimension(:,:),allocatable:: SU_ANTHROL1     !anthropogenic BF emissions (1)
 real(kind=RKIND),dimension(:,:),allocatable:: SU_ANTHROL2     !anthropogenic FF emissions (1)
 real(kind=RKIND),dimension(:,:),allocatable:: SU_SHIPSO2      !SO2 ship emissions (1)
 real(kind=RKIND),dimension(:,:),allocatable:: SU_SHIPSO4      !SO4 ship emissions (1)
 real(kind=RKIND),dimension(:,:),allocatable:: SU_DMSO         !DMS emissions (1)
 real(kind=RKIND),dimension(:,:),allocatable:: SU_AVIATION_LTO !Landing/Take-off aircraft source species (1)
 real(kind=RKIND),dimension(:,:),allocatable:: SU_AVIATION_CDS !Climb/Descent aircraft source species (1)
 real(kind=RKIND),dimension(:,:),allocatable:: SU_AVIATION_CRS !Cruise aircraft source species (1)

!category: EXPORT
 real(kind=RKIND),dimension(:,:,:),allocatable:: SUEM   !Sulfur Emission (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: SUDP   !Sulfate Dry Deposition (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: SUSD   !Sulfate Settling (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: SUWT   !Sulfate Wet Deposition (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: SUSV   !Sulfate Convective Scavenging (Bin %d) (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),allocatable:: SO4EMAN  !SO4 Anthropogenic Emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),allocatable:: SO2EMAN  !SO2 Anthropogenic Emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),allocatable:: SO2EMBB  !SO2 Biomass Burning Emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),allocatable:: SO2EMVN  !SO2 Volcanic (non-explosive) Emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),allocatable:: SO2EMVE  !SO2 Volcanic (explosive) Emissions (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),allocatable:: SUPSO2   !SO2 Prod from DMS Oxidation [column] (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),allocatable:: SUPSO4   !SO4 Prod from All SO2 Oxidation [column] (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:),allocatable:: SUPSO4G  !SO4 Prod from Gaseous SO2 Oxidation [column] (kg m-2 s-1) 
 real(kind=RKIND),dimension(:,:),allocatable:: SUPSO4AQ !SO4 Prod from Aqueous SO2 Oxidation [column] (kg m-2 s-1) 
 real(kind=RKIND),dimension(:,:),allocatable:: SUPSO4WT !SO4 Prod from Aqueous SO2 Oxidation (kg m-2 s-1) 
 real(kind=RKIND),dimension(:,:),allocatable:: SUPMSA   !MSA Prod from DMS Oxidation [column] (kg m-2 s-1) 
 real(kind=RKIND),dimension(:,:),allocatable:: SO2SMASS !SO2 Surface Mass Concentration (kg m-3)     
 real(kind=RKIND),dimension(:,:),allocatable:: SO2CMASS !SO2 Column Mass Density (kg m-2)     
 real(kind=RKIND),dimension(:,:),allocatable:: SO4SMASS !SO4 Surface Mass Concentration (kg m-3)     
 real(kind=RKIND),dimension(:,:),allocatable:: SO4CMASS !SO4 Column Mass Density (kg m-2)     
 real(kind=RKIND),dimension(:,:),allocatable:: DMSSMASS !DMS Surface Mass Concentration (kg m-3)     
 real(kind=RKIND),dimension(:,:),allocatable:: DMSCMASS !DMS Column Mass Density (kg m-2)     
 real(kind=RKIND),dimension(:,:),allocatable:: MSASMASS !MSA Surface Mass Concentration (kg m-3)     
 real(kind=RKIND),dimension(:,:),allocatable:: MSACMASS !MSA Column Mass Density (kg m-3)
 real(kind=RKIND),dimension(:,:,:),allocatable:: PSO2   !SO2 Prod from DMS oxidation (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: PMSA   !MSA Prod from DMS oxidation (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: PSO4   !SO4 Prod from all SO2 oxidation (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: PSO4G  !SO4 Prod from gaseous SO2 oxidation (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: PSO4WET!SO4 Prod from wet SO2 oxidation (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: PSO4AQ !SO4 Prod from aqueous SO2 oxidation (kg m-2 s-1)
 real(kind=RKIND),dimension(:,:,:),allocatable  :: SUCONC        !SO4 Aerosol Mass Concentration (kg m-3)
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: SUEXTCOEF     !SO4 Extinction Coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: SUEXTCOEFRH20 !SO4 Extinction Coefficient - Fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: SUEXTCOEFRH80 !SO4 Extinction Coefficient - Fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: SUSCACOEF     !SO4 Scattering Coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: SUSCACOEFRH20 !SO4 Scattering Coefficient - Fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: SUSCACOEFRH80 !SO4 Scattering Coefficient - Fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: SUBCKCOEF     !SO4 Backscatter Coefficient (m-1 sr-1) 

 real(kind=RKIND),dimension(:,:),allocatable:: SUANGSTR          !SO4 Angstrom parameter [470-870 nm] (1)
 real(kind=RKIND),dimension(:,:),allocatable:: SUFLUXU           !SO4 column u-wind mass flux (kg m-1 s-1)
 real(kind=RKIND),dimension(:,:),allocatable:: SUFLUXV           !SO4 column v-wind mass flux (kg m-1 s-1)
 real(kind=RKIND),dimension(:,:),allocatable:: SO4MASS           !SO4 Aerosol Mass Mixing Ratio (kg kg-1)
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: SUEXTTAU      !SO4 Extinction AOT (1)
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: USTEXTTAU     !SO4 Extinction AOT Stratosphere (1)
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: SUSCATAU      !SO4 Scattering AOT (1)
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: SUSTSCATAU    !SO4 Scattering AOT Stratosphere (1)
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: SO4SAREA      !SO4 Surface Area Density (m2 m-3 )
 real(kind=RKIND),dimension(:,:,:,:),allocatable:: SO4SNUM       !SO4 Number Density (m-3)

 real(kind=RKIND),dimension(:,:,:),allocatable:: DMS             !Dimethylsulphide (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: SO2             !Sulphur dioxide (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: SO4             !Sulphate aerosol (kg kg-1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: MSA             !Methanesulphonic acid (kg kg-1) 
 real(kind=RKIND),dimension(:,:,:),allocatable:: H2O2_INIT       !private H2O2 (kg kg-1)

 end type SU2G_State

!=================================================================================================================
 end module SU2G_StateSpecs
!=================================================================================================================
