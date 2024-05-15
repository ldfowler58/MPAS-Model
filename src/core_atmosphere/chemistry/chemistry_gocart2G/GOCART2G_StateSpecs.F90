!=================================================================================================================
 module GOCART2G_StateSpecs
 use mpas_kind_types,only: RKIND
 public
 save

!this module is the state variable specification file for AODs. it is the same as GOCART2G_StateSpecs.rc
!in the GOCART-2G directory ./GOCART-2G/ESMF/GOCART2G_GridComp.

!schema_version: 2.0.0
!component: GOCART2G


 type GOCART2G_State

!category: IMPORT
 real(kind=RKIND),dimension(:,:,:),allocatable:: DELP         !pressure_thickness (Pa)
!real(kind=RKIND),dimension(:,:,:),allocatable:: RH2          !rel_Hum_after_moist (1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: AIRDENS      !moist air density  (kg/m^3)
 real(kind=RKIND),dimension(:,:,:),allocatable:: T            !air_temperature (K)
 real(kind=RKIND),dimension(:,:,:),allocatable:: PLE          !air pressure (Pa)

!category: EXPORT
 real(kind=RKIND),dimension(:,:,:),allocatable:: PSO4TOT        !Total Sulfate Produced in GOCART (kg m-2 s-1)
!........................................................................................
 real(kind=RKIND),dimension(:,:),allocatable:: TOTEXTTAU        !Total Aerosol Extinction AOT [550 nm]
 real(kind=RKIND),dimension(:,:),allocatable:: TOTSTEXTTAU      !Total Aerosol Extinction AOT [550 nm] Stratosphere
 real(kind=RKIND),dimension(:,:),allocatable:: TOTSCATAU        !Total Aerosol Scattering AOT [550 nm]
 real(kind=RKIND),dimension(:,:),allocatable:: TOTSTSCATAU      !Total Aerosol Scattering AOT [550 nm] Stratosphere
 real(kind=RKIND),dimension(:,:),allocatable:: TOTEXTT25        !Total Aerosol Extinction AOT [550 nm] - PM2.5
 real(kind=RKIND),dimension(:,:),allocatable:: TOTSCAT25        !Total Aerosol Extinction AOT [550 nm] - PM2.5
 real(kind=RKIND),dimension(:,:),allocatable:: TOTEXTTFM        !Total Aerosol Extinction AOT [550 nm] - PM1.0
 real(kind=RKIND),dimension(:,:),allocatable:: TOTSCATFM        !Total Aerosol Extinction AOT [550 nm] - PM1.0
 real(kind=RKIND),dimension(:,:),allocatable:: TOTANGSTR        !Total Aerosol Angstrom parameter [470-870 nm]
 real(kind=RKIND),dimension(:,:,:),allocatable:: TOTEXTCOEF     !Total Aerosol Extinction coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: TOTEXTCOEFRH20 !Total Aerosol Extinction coefficient - Fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: TOTEXTCOEFRH80 !Total Aerosol Extinction coefficient - Fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: TOTSCACOEF     !Total Aerosol Scattering coefficient (m-1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: TOTSCACOEFRH20 !Total Aerosol Scattering coefficient - Fixed RH=20% (m-1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: TOTSCACOEFRH80 !Total Aerosol Scattering coefficient - Fixed RH=80% (m-1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: TOTBCKCOEF     !Total Aerosol Single Scattering Backscatter coefficient (m-1 sr-1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: TOTABCKTOA     !Total Attenuated Backscatter Coefficient from TOA [532nm] (m-1 sr-1)
 real(kind=RKIND),dimension(:,:,:),allocatable:: TOTABCKSFC     !Total Attenuated Backscatter Coefficient from surface [532nm] (m-1 sr-1)
 real(kind=RKIND),dimension(:,:),allocatable:: PM               !Total reconstructed PM (kg m-3)
 real(kind=RKIND),dimension(:,:),allocatable:: PM_RH35          !Total reconstructed PM(RH=35%) (kg m-3)
 real(kind=RKIND),dimension(:,:),allocatable:: PM_RH50          !Total reconstructed PM(RH=50%) (kg m-3)
 real(kind=RKIND),dimension(:,:),allocatable:: PM25             !Total reconstructed PM2.5 (kg m-3)
 real(kind=RKIND),dimension(:,:),allocatable:: PM25_RH35        !Total reconstructed PM2.5(RH=35%) (kg m-3)
 real(kind=RKIND),dimension(:,:),allocatable:: PM25_RH50        !Total reconstructed PM2.5(RH=50%) (kg m-3)

!category: INTERNAL

 end type GOCART2G_State

!=================================================================================================================
 end module GOCART2G_StateSpecs
!=================================================================================================================
